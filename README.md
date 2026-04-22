# pygamd_v_me_50_meal: 一款用于 pygamd 分子动力学模拟软件的脚本。

由于 PYGAMD 未提供适用于生物大分子凝聚体模拟的数据处理软件，因此我整合了我所使用过的脚本，集合成一个脚本，并提供了丰富的选项以供选择。

本软件分为两部分：分子动力学模拟前的文件生成、分子动力学模拟后的文件处理。在拿到 PDB 文件后，可以使用本软件进行模拟体系的构建、运行脚本的生成、分子动力学模拟的运行。再得到结果后，可以继续使用本软件进行模拟所得到的 XML 文件的处理。

---
## 安装
可直接使用 `pip` 安装：

    pip install git+https://gitee.com/lyxlyxlyxxx/pygamd_v_me_50_meal.git

或者

    git clone https://gitee.com/lyxlyxlyxxx/pygamd_v_me_50_meal.git
    cd pygamd_v_me_50_meal
    pip install -e .

---
## 分子动力学模拟前的文件生成

### 生成 XML 文件
从 PDB 文件生成模拟所需要的粗粒化 XML 文件，目前支持蛋白质、DNA 和 RNA 粗粒化模型的生成。
#### 使用
    v50 -pdb2xml /path/to/pdb_file.pdb
    v50 -p /path/to/system -pdb2xml filename.pdb

例如，要将 `/home/protein` 文件夹下的 **1kx5.pdb** 转化为 **xml** 文件，输入 

    v50 -pdb2xml /home/protein/1kx5.pdb
    或 v50 -p /home/protein -pdb2xml 1kx5.pdb

即可生成 **1kx5.xml** 文件。

##### 下列功能直接在运行过程中以交互式选项呈现
- 设置弹性网络（ENM）
- 选择 DNA/RNA 粗粒化模型
- 生成模拟脚本
- 设置盒子大小（盒子为正方体盒子，大小默认为 100nm，如需要设置其他长度，请使用-box_size选项）
- 设置刚体结构域
- 对结构域单独设置粒子类型

##### 蛋白质粗粒化模型：一个氨基酸粗粒化为一个粒子：Dignon G L, Zheng W, Kim Y C, Best R B, Mittal J. PLOS Computational Biology, 2018, 14(1): e1005941
- **HPS 模型**：适用于 HPS-Urry、CALVADOS 系列等力场
- **Mpipi 模型**：另一种蛋白质粗粒化模型

##### DNA 粗粒化模型：目前支持两种模型
1. 3SPN 模型 
   > Knotts T A, Rathore N, Schwartz D C, De Pablo J J. A coarse grain model for DNA[J]. _The Journal of Chemical Physics_, **2007**, 126(8): 084901.

2. Mittal 2 Bead 模型：
   > Kapoor U, Kim Y C, Mittal J. Coarse-Grained Models to Study Protein–DNA Interactions and Liquid–Liquid Phase Separation[J]. _Journal of Chemical Theory and Computation_, **2023**: acs.jctc.3c00525.

注意：目前“Mittal 2 Bead 模型”还未能成功在 PYGAMD 中实现。

##### RNA 粗粒化模型
- **3SPN 模型**：基于 DNA 的 3SPN 模型扩展
- **CAVADOS-RNA 模型**：（未提供完整实现）

#### 高级选项

**参数说明：**

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `-box_size` | 设置盒子大小（正方体盒子，单位：nm） | `100` |
| `-add_enm_bond` | 预定义弹性网络结构域，格式如 "1-50,51-100" | `None` |
| `-add_rigid_body` | 预定义刚体结构域，格式如 "1-50,51-100" | `None` |
| `-gen_run_file` | 是否自动生成运行文件，设为 `t` 启用 | `unset` |
| `-dna_model` | 指定 DNA 模型，1=3SPN, 2=2BeadMittal | `None` |

使用示例：
```bash
# 设置盒子大小
v50 -pdb2xml /path/to/pdb_file.pdb -box_size 150

# 预定义弹性网络
v50 -pdb2xml /path/to/pdb_file.pdb -add_enm_bond "1-50,51-100"

# 预定义刚体结构域
v50 -pdb2xml /path/to/pdb_file.pdb -add_rigid_body "1-50,51-100"

# 自动生成运行文件
v50 -pdb2xml /path/to/pdb_file.pdb -gen_run_file t

# 指定 DNA 模型
v50 -pdb2xml /path/to/pdb_file.pdb -dna_model 1  # 1=3SPN, 2=2BeadMittal
```

---
## 分子动力学模拟后的文件处理

### 1. 坐标提取
注意：在进行任何数据处理前，都需要先需要提取坐标。

1. 由于 PYGAMD 未提供拓扑文件，因此需要标注您的体系信息。
请确保提供路径的最后一级为 "数字+分子名称-分子长度" 的格式，并且按照 xml 文件中的顺序排好。  
例如，体系中含有 40 个长度为 256 的 A 分子、20 个长度为 512 的 B 分子、30 个长度为 729 的 C 分子，
须将体系命名为"40A-256+20B-512+30C-729" (本手册将以`40A-256+20B-512+30C-729`文件夹为例). 
2. 请将分子动力学模拟得到的 XML 文件放置在提供的目录下或者在该目录下的 xml 文件夹内。
3. 进行坐标提取：`v50 -p 40A-256+20B-512+30C-729 -xyz`
4. #### 可选参数 

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `-remove_enm` | 是否移除弹性键 | 无 |
| `-remove_condensate_pbc` | 是否去除 PBC，并将最大的凝聚体移动到盒子中央 | 无 |

使用示例：
```
v50 -p 40A-256+20B-512+30C-729 -xyz -remove_enm
v50 -p 40A-256+20B-512+30C-729 -xyz -remove_condensate_pbc
```


### 2. 计算接触图 (Contact Map)
注意，在计算分子间的 Contact Map 时，已经排除了单个分子内的相互作用。 

    v50 -p 40A-256+20B-512+30C-729 -cm

**参数说明：**

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `-cm` | 是否计算接触图（contact map），计算后会自动绘图 | 无 |
| `-r_cut` | 设置计算 Contact Map 的 $r_{cut}$ | 4.0 |

使用示例：
```
v50 -p 40A-256+20B-512+30C-729 -cm -r_cut 4.0
```

### 3. 计算 Rg、RMSD、RMSF

**参数说明：**

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `-rg` | 是否计算回旋半径 (Rg) | 无 |
| `-rmsd` | 是否计算均方根偏差 (RMSD) | 无 |
| `-rmsf` | 是否计算均方根涨落 (RMSF) | 无 |
| `-ref` | 参考结构文件路径，用于 RMSD 和 RMSF 计算 | `None` |

使用示例：
```
v50 -p 40A-256+20B-512+30C-729 -rg
v50 -p 40A-256+20B-512+30C-729 -rmsd -ref reference_file.xml
v50 -p 40A-256+20B-512+30C-729 -rmsf -ref reference_file.xml
```

在运行过程中会询问需要计算的分子，ref 提供的哪个就计算哪个

### 4. 计算体系的质量数密度分布  
注意：在使用前请先进行对凝聚体PBC的去除。

**参数说明：**

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `-mass_density` | 是否计算质量数密度分布 | 无 |

使用示例：
```
v50 -mass_density
```

## 其他分析

### 序列分析

序列分析模块提供了两种分析功能：NCPR（净电荷分析）和芳香性分析。

**功能特点：**
- 统一的命令行接口，通过交互式界面选择分析类型
- 支持直接传入序列字符串或序列文件路径
- 边缘自适应补全算法，保持输出长度等于输入长度
- 分段着色（NCPR：红色正电荷，蓝色负电荷；芳香性：紫色）
- 面积填充，增强可视化效果
- 自动调整 x 轴刻度间隔，适应不同长度的序列

**使用方法：**

1. 直接传入序列字符串：
```
v50 -seq_analysis -seq "MKVDELVQGLLKQISAEELKKARNEIARQHLEKTHQDLKKDILTYLTDRQIKQLEDAFQKLLAEKTEENKLAQAVENSLGQLEEKLKEA"
```

2. 通过文件路径传入序列：
```
v50 -seq_analysis -seq pro_sequence.txt
```

3. 指定滑动窗口大小（默认 15）：
```
v50 -seq_analysis -seq pro_sequence.txt -seq_window 20
```

4. 指定输出文件路径：
```
v50 -seq_analysis -seq pro_sequence.txt -seq_output pro_analysis.png
```

**参数说明：**

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `-seq_analysis` | 是否进行序列分析 | 无 |
| `-seq` | 序列字符串或序列文件路径 | `None` |
| `-seq_window` | 滑动窗口大小 | `15` |
| `-seq_output` | 输出图片路径 | `自动生成` |

**使用流程：**
1. 运行命令后，会提示选择分析类型
2. 输入 `1` 选择 NCPR 分析，输入 `2` 选择芳香性分析
3. 分析完成后会生成相应的 PNG 格式图片

## 故障排除

### 常见问题及解决方法

1. **PDB 文件转换失败**
   - 检查 PDB 文件格式是否正确
   - 确保 PDB 文件中包含完整的原子信息
   - 对于多链结构，确保使用 "TER" 标签分隔

2. **弹性网络设置错误**
   - 确保结构域范围设置正确（从 1 开始计数）
   - 检查 PDB 文件中的残基编号是否连续

3. **序列分析失败**
   - 确保输入的序列只包含标准氨基酸代码
   - 对于文件输入，确保文件编码为 UTF-8
   - 检查滑动窗口大小是否合理（建议 10-30）

4. **坐标提取失败**
   - 确保目录命名格式正确（如 "40A-256+20B-512"）
   - 检查 XML 文件是否存在于正确位置
   - 对于大型体系，可能需要增加内存限制

## 版本信息

### v1.0.0
- 初始版本
- 支持蛋白质、DNA、RNA 粗粒化模型生成
- 实现序列分析功能（NCPR 和芳香性分析）
- 提供接触图、Rg、RMSD、RMSF 计算
- 支持质量密度分布分析

## 示例

### 示例 1：生成蛋白质 XML 文件并添加弹性网络

```bash
# 生成 XML 文件并添加弹性网络
v50 -pdb2xml /path/to/protein.pdb

# 运行过程中选择：
# 1. 蛋白质模型：HPS
# 2. 是否设置弹性网络：y
# 3. 结构域范围：1-100,101-200
# 4. 是否生成运行文件：y
```

### 示例 2：分析蛋白质序列的 NCPR 和芳香性

```bash
# 分析序列文件
v50 -seq_analysis t -seq protein_sequence.txt -seq_window 20

# 运行过程中选择：
# 1. NCPR 分析
# 2. 芳香性分析

# 生成两个分析结果图片
```

### 示例 3：分析模拟结果

```bash
# 提取坐标
v50 -p 40A-256+20B-512 -xyz t -remove_condensate_pbc t

# 计算接触图
v50 -p 40A-256+20B-512 -cm t -r_cut 4.5

# 计算 Rg
v50 -p 40A-256+20B-512 -rg t
```

## 致谢

本软件参考了以下研究工作：

1. Dignon G L, Zheng W, Kim Y C, Best R B, Mittal J. PLOS Computational Biology, 2018, 14(1): e1005941
2. Knotts T A, Rathore N, Schwartz D C, De Pablo J J. A coarse grain model for DNA[J]. _The Journal of Chemical Physics_, **2007**, 126(8): 084901.
3. Kapoor U, Kim Y C, Mittal J. Coarse-Grained Models to Study Protein–DNA Interactions and Liquid–Liquid Phase Separation[J]. _Journal of Chemical Theory and Computation_, **2023**: acs.jctc.3c00525.

---

**注意**：本软件仅供研究使用，如有任何问题或建议，请联系开发者。