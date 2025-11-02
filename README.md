# pygamd_v_me_50_meal: 一款用于 pygamd 分子动力学模拟软件的脚本。

由于 PYGAMD 未提供适用于生物大分子凝聚体模拟的数据处理软件，因此我整合了我所使用过的脚本，集合成一个脚本，并提供了丰富的选项以供选择。

本软件分为两部分：分子动力学模拟前的文件生成、分子动力学模拟后的文件处理。在拿到 PDB 文件后，可以使用本软件进行模拟体系的构建、运行脚本的生成、分子动力学模拟的运行。再得到结果后，可以继续使用本软件进行模拟所得到的 XML 文件的处理。

---
## 安装
可直接使用 `pip` 安装：

    pip install git+https://gitee.com/lyxlyxlyxxx/pygamd_v_me_50_meal.git
    pip install git+https://github.com/xcslyx/pygamd_v_me_50_meal

或者

    git clone https://gitee.com/lyxlyxlyxxx/pygamd_v_me_50_meal.git
    cd pygamd_v_me_50_meal
    pip install -e .

---
## 分子动力学模拟前的文件生成

### pygamd_v_me_50_meal/simulate_creation/xml_generator.py
XMLGenerator 类用于从 PDB 文件生成模拟所需要的粗粒化 XML 文件，目前支持蛋白质与 DNA 粗粒化模型的生成。
#### 使用
    v50 -pdb2xml /path/to/pdb_file.pdb
    v50 -p /path/to/system -pdb2xml filename.pdb

例如，要将 `/home/protein` 文件夹下的 **1kx5.pdb** 转化为 **xml** 文件，输入 

    v50 -pdb2xml /home/protein/1kx5.pdb
    或 v50 -p /home/protein -pdb2xml 1kx5.pdb

即可生成 **1kx5.xml** 文件。

##### 下列功能直接在运行过程中以交互式选项呈现
设置弹性网络、选择 DNA 粗粒化模型、生成模拟脚本、设置盒子大小（盒子为正方体盒子，大小默认为 100nm，如需要设置其他长度，请使用-box_size选项）

##### 蛋白质粗粒化模型：一个氨基酸粗粒化为一个粒子：Dignon G L, Zheng W, Kim Y C, Best R B, Mittal J. PLOS Computational Biology, 2018, 14(1): e1005941
##### DNA 粗粒化模型：目前支持两种模型
1. 3SPN 模型
> Knotts T A, Rathore N, Schwartz D C, De Pablo J J. A coarse grain model for DNA[J]. _The Journal of Chemical Physics_, **2007**, 126(8): 084901.
2. Mittal 2 Bead 模型：
>Kapoor U, Kim Y C, Mittal J. Coarse-Grained Models to Study Protein–DNA Interactions and Liquid–Liquid Phase Separation[J]. _Journal of Chemical Theory and Computation_, **2023**: acs.jctc.3c00525.

注意：目前“Mittal 2 Bead 模型”还未能成功在 PYGAMD 中实现。
    
    
## 分子动力学模拟后的文件处理

### 1. 坐标提取
注意：在进行任何数据处理前，都需要先需要提取坐标。

1. 由于 PYGAMD 未提供拓扑文件，因此需要标注您的体系信息。
请确保提供路径的最后一级为 "数字+分子名称-分子长度" 的格式，并且按照 xml 文件中的顺序排好。  
例如，体系中含有 40 个长度为 256 的 A 分子、20 个长度为 512 的 B 分子、30 个长度为 729 的 C 分子，
须将体系命名为"40A-256+20B-512+30C-729" (本手册将以`40A-256+20B-512+30C-729`文件夹为例). 
2. 请将分子动力学模拟得到的 XML 文件放置在提供的目录下或者在该目录下的 xml 文件夹内。
3. 进行坐标提取：`v50 -p 40A-256+20B-512+30C-729 -xyz t`
4. #### 可选参数 
   `-remove_enm`：如果加入了弹性网络，可用此选项移除弹性键。 如：`v50 -p 40A-256+20B-512+30C-729 -xyz t -remove_enm t`
   `-remove_condensate_pbc`：若设置为 t，则会去除 PBC，并将最大的凝聚体移动到盒子中央。如：`v50 -p 40A-256+20B-512+30C-729 -xyz t -remove_condensate_pbc t`


[//]: # (### GetSequence)

[//]: # (GetSequence 类用于提取 xml 或 pdb 文件的序列，根据文件后缀自动检测文件类型。)

[//]: # (#### 使用)

[//]: # (你可以直接使用 `-get_seq` 选项：)

[//]: # ()
[//]: # (    v50 -get_seq /path/to/file.xml)

[//]: # (    v50 -get_seq /path/to/file.pdb)

[//]: # ()
[//]: # (也可以通过 `-p` 提供路径，然后使用 `-get_seq` 提供文件名：)

[//]: # ()
[//]: # (    v50 -p /path/to/system -get_seq filename.xml)

[//]: # (    v50 -p /path/to/system -get_seq filename.pdb)

### 2. 计算接触图 (Contact Map)
注意，在计算分子间的 Contact Map 时，已经排除了单个分子内的相互作用。 

    v50 -p 40A-256+20B-512+30C-729 -cm t

我们提供了以下选项：
1. `-r_cut`：设置计算 Contact Map 的 $r_{cut}$。

    `v50 -p 40A-256+20B-512+30C-729 -cm t -r_cut 4.0`

### 3. 计算 Rg、RMSD、RMSF

    v50 -p 40A-256+20B-512+30C-729 -rg t
    v50 -p 40A-256+20B-512+30C-729 -rmsd t -ref reference_file.xml
    v50 -p 40A-256+20B-512+30C-729 -rmsf t -ref reference_file.xml

在运行过程中会询问需要计算的分子，ref 提供的哪个就计算哪个

### 4. 计算体系的质量数密度分布  
注意：在使用前请先进行对凝聚体PBC的去除。

    v50 -mass_density t

