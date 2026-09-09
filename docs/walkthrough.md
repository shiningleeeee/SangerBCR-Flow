# 使用教程（WSL / Linux 环境）

> 本教程由仓库最初的单文件 `code` 整理而来。脚本已拆分为本仓库 `scripts/` 下的三个独立文件（代码内容**未做任何修改**），安装与运行步骤保持原样。
>
> 运行前提：三个脚本放在 `~/bcr/scripts`（脚本用相对路径 `../raw_data`、`../results` 访问同级文件夹），因此目录结构为：
>
> ```
> bcr/
> ├── igblast_tool/   # 02 自动下载的 IgBLAST（含参考库）
> ├── raw_data/       # 放入 .ab1 原始测序文件
> ├── references/     # 02 自动准备的 IMGT 参考基因
> ├── results/        # 01/02/03 的输出
> └── scripts/        # 本仓库 scripts/ 下的三个脚本
> ```

## 1. WSL 的安装及关闭

```bash
#安装wsl至D盘wsl文件夹
wsl --install --web-download --location D:\wsl
#查看所有发行版详细信息,简写wsl -l -v
wsl --list --verbose
#关闭特定发行版
wsl --terminate Ubuntu
#关闭所有
wsl --shutdown
```

## 2. 文件夹的建立及包的安装

```bash
# 进入 Ubuntu/WSL 后，先更新系统软件源并安装 IgBLAST 可能需要的系统运行库
sudo apt update
sudo apt install -y libgomp1
#创建文件夹
mkdir bcr
cd bcr
mkdir igblast_tool raw_data references results scripts
#将原始的测序文件ab1放入raw_data文件夹下
#安装miniconda3
cd ~
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
#装完后执行
source ~/.bashrc
#成功后运行下面代码，看到版本说明安装成功
conda --version
#创建分析环境
conda create -n bcr python=3.9
#创建完成后激活
conda activate bcr
#安装必要的Python包
cd ~/bcr
pip install biopython pandas seaborn matplotlib numpy
#（等价方式：仓库附带的 requirements.txt，与上面命令行效果一致）
# pip install -r requirements.txt
```

> 说明：02 脚本会用 Perl 运行 IgBLAST 自带的 `edit_imgt_file.pl`，Ubuntu 默认自带 perl；若提示缺少 perl，先 `sudo apt install -y perl`。

## 3. 脚本的准备

把仓库 `scripts/` 下的三个文件放进 `~/bcr/scripts`：

```bash
# 方式一（推荐）：把仓库克隆/下载到 ~/bcr 下，scripts/ 即仓库中的 scripts/
cd ~/bcr
git clone https://github.com/shiningleeeee/SangerBCR-Flow.git .

# 方式二：手动复制
# 将仓库中的 scripts/01_ab1_to_fasta.py、scripts/02_run_igblast.sh、
# scripts/03_analyze_clones.py 复制到 ~/bcr/scripts/
```

三个脚本的作用与数据流：

| 脚本 | 作用 | 主要输出（均在 `results/`） |
| --- | --- | --- |
| `01_ab1_to_fasta.py` | 读取 `raw_data/*.ab1`，按质量修剪（abi-trim），过滤过短序列 | `01_cleaned_sequences.fasta`、`01_raw_sequences.fasta`、`01_sequence_summary.tsv` |
| `02_run_igblast.sh` | 首次运行自动下载 IgBLAST 1.22.0 与 IMGT 人源 V/D/J（+C 区）参考并建库，随后对清洗后序列做 V(D)J 注释 | `02_igblast_results.tsv`（IgBLAST AIRR TSV，outfmt 19） |
| `03_analyze_clones.py` | 过滤 productive 序列；基因使用 / SHM / CDR3 长度分布图；heavy-only 克隆分型（同 V、同 J、同 junction 长度 + junction nt 相似度）与克隆饼图 | `analysis_output/` 下的 PNG 图与 `Heavy_Clonotypes.tsv`、`Clone_Summary.tsv` |

## 4. 执行脚本

```bash
cd ~/bcr/scripts
python3 01_ab1_to_fasta.py
bash 02_run_igblast.sh
python3 03_analyze_clones.py
```

## 5. 注意事项

```bash
#results文件夹下的结果请及时转移，避免新运行产生的结果混淆
#后续只需把原始测序文件放入raw_data文件夹下，先source ~/.bashrc , conda activate bcr并且cd ~/bcr/scripts依次执行三个脚本即可
```

- 可调参数（如 `MIN_LEN`、`MAX_DISTANCE`、`TOP_N_CLONES_IN_PIE` 等）在 `03_analyze_clones.py` 顶部「路径和参数」区域修改（此版重组未改动任何代码）。
