# 历史 v3 离线转换

新增 `tools/fpd2vtk_v3.py`，配套只读格式模块 `tools/fpd_format_v3.py`。
历史布局固定自提交 `42c279a`：40 字节头部，3 个流体数组，12 个粒子数组，8 字节校验和。
文件长度为 `48+24*size+96*N`。当前 v4 工具和生产程序没有改动。

新入口在独立进程中将共用转换模块的格式读取替换为 v3；顶层替换也在 spawn worker 中执行。
内存估算、并行、原子写入、增量索引及写出实现共用当前工具。
转换不输出历史受力，只保留粒子位置、velocity、unwrapped。
现有 `verify_vtk.py` 仍仅支持 v4；本次 v3 数值往返验收使用此目录的 `check.py`。

## 验收

- pvpython 下的 `tools/fpd2vtk_v3.py --self-test` 通过：轴序、插值、无 force 的 VTK 字段。
  共用自检中的“v4 粒子字段”表示当前写出字段集合，输入读取采用 v3 布局。
- `check.py` 独立构造 N=0/2、8×4×2 的已知值 v3 文件，全部字段、偏移、长度及版本隔离通过。
- CLI `-j 2` 完整转换 3 个小文件，通过 VTK 读回的流体、粒子、空粒子及无 force 检查；
  VTM 引用、PVD 时间轴、重跑跳过和索引合并通过，见 `check.log`。
- 新模块读取实际目录的全部 401 个文件头部；布局末尾与各文件长度一致，见 `input_scan.log`。
  未执行大文件转换或完整校验和验证。

复现小文件检查：

```bash
PV=/home/doll/Software/ParaView-6.2.0-RC1-MPI-Linux-Python3.12-x86_64
"$PV/bin/pvpython" tools/fpd2vtk_v3.py --self-test
"$PV/bin/pvpython" baseline/vtk_v3_20261027/check.py
```

## 用户转换命令

仓库根目录执行；使用全新输出目录，以免跳过已有旧产物。默认自动选择并行度，也可用 `-j` 指定。

```bash
PV=/home/doll/Software/ParaView-6.2.0-RC1-MPI-Linux-Python3.12-x86_64
DATA="$HOME/mnt/Storage/FPD_data/out/test-20261021-g=0.06"
"$PV/bin/pvpython" tools/fpd2vtk_v3.py "$DATA"/512x512x128_[0-9]*.fpd \
  --config "$DATA/512x512x128.config.used" \
  -o vis/test-20261021-g=0.06-v3/
```

ParaView 打开 `vis/test-20261021-g=0.06-v3/512x512x128.pvd`。
