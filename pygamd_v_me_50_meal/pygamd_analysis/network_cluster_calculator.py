import os
import numpy as np
import networkx as nx
from sklearn.cluster import SpectralClustering, AgglomerativeClustering
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
from tqdm import tqdm

class NetworkClusterCalculator:
    def __init__(self, path, data, node_mol_type, edge_mol_type, distance_threshold: float=10.0):
        self.path = path
        self.data = data
        self.node_mol_type = node_mol_type
        self.edge_mol_type = edge_mol_type
        self.distance_threshold = distance_threshold

        self.chain_path = os.path.join(self.path, "chain_xyz/")

        self.node_mol_count = self.data.mol_class_dict[self.node_mol_type][0]
        self.edge_mol_count = self.data.mol_class_dict[self.edge_mol_type][0]

        self.output_path = os.path.join(self.path, "network_cluster/")
        os.makedirs(self.output_path, exist_ok=True)

    def build_network(self, name):
        x_mat = eval(open(os.path.join(self.chain_path, name), 'r').read())

        G = nx.Graph()

        for i in range(self.node_mol_count):
            G.add_node(f"{self.node_mol_type}_{i}")

        node_coords = x_mat[self.node_mol_type]
        edge_coords = x_mat[self.edge_mol_type]

        for i, node_coord in enumerate(node_coords):
            for j, edge_coord in enumerate(edge_coords):
                dist = np.linalg.norm(np.mean(node_coord, axis=0) - np.mean(edge_coord, axis=0))
                if dist < self.distance_threshold:
                    G.add_edge(f"{self.node_mol_type}_{i}", f"{self.edge_mol_type}_{j}", weight=dist)

        return G

    def cluster_network(self, G, n_clusters=None, method='spectral'):
        adj_matrix = nx.to_numpy_array(G)

        if n_clusters is None:
            n_clusters = self._estimate_clusters(adj_matrix)

        if method == 'spectral':
            clustering = SpectralClustering(
                n_clusters=n_clusters,
                affinity='precomputed',
                random_state=42
            )
        elif method == 'hierarchical':
            clustering = AgglomerativeClustering(
                n_clusters=n_clusters,
                metric='euclidean',
                linkage='ward'
            )
        else:
            raise ValueError(f"Unknown clustering method: {method}")

        labels = clustering.fit_predict(adj_matrix)

        return labels

    def _estimate_clusters(self, adj_matrix):
        best_k = 2
        best_score = -1

        for k in range(2, min(10, self.node_mol_count)):
            try:
                clustering = SpectralClustering(
                    n_clusters=k,
                    affinity='precomputed',
                    random_state=42
                )
                labels = clustering.fit_predict(adj_matrix)
                if len(set(labels)) > 1:
                    score = silhouette_score(adj_matrix, labels, metric='precomputed')
                    if score > best_score:
                        best_score = score
                        best_k = k
            except Exception:
                continue

        return best_k

    def calculate_parallel(self):
        file_list = [f for f in os.listdir(self.chain_path) if f.endswith('.npy') or f.endswith('.txt') or f.endswith('.xml')]

        results = {}
        for name in tqdm(file_list, desc="Building networks"):
            G = self.build_network(name)
            labels = self.cluster_network(G)

            node_labels = {f"{self.node_mol_type}_{i}": labels[i] for i in range(self.node_mol_count)}
            results[name] = {
                'graph': G,
                'labels': labels,
                'node_labels': node_labels
            }

        self._save_results(results)
        self._plot_clusters(results)

        return results

    def _save_results(self, results):
        summary_path = os.path.join(self.output_path, "clustering_summary.txt")

        with open(summary_path, 'w') as f:
            f.write("Network Clustering Summary\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"Node molecule type: {self.node_mol_type}\n")
            f.write(f"Edge molecule type: {self.edge_mol_type}\n")
            f.write(f"Distance threshold: {self.distance_threshold}\n\n")

            for name, data in results.items():
                f.write(f"\nFile: {name}\n")
                f.write(f"Number of clusters: {len(set(data['labels']))}\n")
                f.write(f"Cluster distribution:\n")
                unique, counts = np.unique(data['labels'], return_counts=True)
                for cluster_id, count in zip(unique, counts):
                    f.write(f"  Cluster {cluster_id}: {count} nodes\n")

    def _plot_clusters(self, results):
        n_files = len(results)
        fig, axes = plt.subplots(1, min(n_files, 4), figsize=(5 * min(n_files, 4), 5))

        if n_files == 1:
            axes = [axes]
        else:
            axes = axes.flatten()

        for idx, (name, data) in enumerate(results.items()):
            if idx >= 4:
                break

            G = data['graph']
            labels = data['labels']

            if G.number_of_nodes() == 0:
                continue

            pos = nx.spring_layout(G, k=2, iterations=50, seed=42)

            nx.draw_networkx_nodes(G, pos, ax=axes[idx],
                                   node_color=labels,
                                   cmap=plt.cm.tab10,
                                   node_size=100)
            nx.draw_networkx_edges(G, pos, ax=axes[idx], alpha=0.3)

            axes[idx].set_title(f"{name[:15]}...\nClusters: {len(set(labels))}")
            axes[idx].axis('off')

        plt.tight_layout()
        plt.savefig(os.path.join(self.output_path, "network_clusters.png"), dpi=150)
        plt.close()

    def get_cluster_statistics(self, results):
        all_labels = []
        for data in results.values():
            all_labels.extend(data['labels'])

        unique, counts = np.unique(all_labels, return_counts=True)

        stats = {}
        for cluster_id, count in zip(unique, counts):
            stats[f"Cluster_{cluster_id}"] = {
                'count': int(count),
                'percentage': float(count / len(all_labels) * 100)
            }

        return stats
