# v4 检查点删除粒子受力

任务开始于 2026-10-25，验证完成于 2026-10-26。旧实现基于 `42c279a`，改动前重新构建并复制
为 `/tmp/fpd-v3-before`；新实现只读写 v4。GPU 为 RTX 3060，运行显式设置 `ACC_DEVICE_TYPE=nvidia`。
GPU 在沙箱内不可见，数值检查在获准的沙箱外完成，无 CPU 回退。

## 实现与覆盖

文件只保存流体 vx/vy/vz 和粒子 R/Ru/V；长度为 `48+24*size+72*N`，每份减少 `24*N` 字节。
C++/Python 同步删除受力字段；删除终态及 PCG 失败落盘前的受力补算、VTK force、
`--verify-forces`。运行时力数组及 ST_PARTICLE 的上传/下载语义保持原样。

覆盖减少：无法再用存档受力验证势参数与文件的一致性；不再输出可视化受力。
覆盖保留：独立势导数、GPU/CPU 力装配、近壁受力、力投影、重叠相场、壁面噪声及 PCG 噪声判据。
新增常驻 CPU 入口 `fpd_check --check-io <目录>`，由 `interop.py` 准备已知输入。

## 已执行

- 默认生产构建及 fpd_check、fpd_tool、fpd_calib 构建通过（build.log）。
- check-config、check、check-potential、check-wall、pcg-noise 均返回 0，各自同名日志。
- Python→C++ 独立已知值核对和 C++→Python 整文件比较：N=0/2、全部保留字段、文件长度、
  校验和、v3 拒绝、损坏/截断/尾部多余数据拒绝均通过；索引 (3,1,1) 的 vx=3001001（interop.log）。
- wall：32³ 单粒子，z=4.7，WCA 壁面受力 max|F|≈12；第 0/5/10 步新旧共同状态逐位一致。
- pair_disjoint：64×32×32 双粒子，中心距 20，模板半宽 8，支撑域不重叠；WCA sigma=20、eps=0.01，
  max|F|≈0.012。第 0/5/10 步新旧共同状态逐位一致（compare.log）。
- wall 连续 10 步与第 5 步重启到 10 步：整个 v4 文件逐位相同。
- 零步运行两版均不落盘；强制 pcg_tol=1e-30，两版均以 3 退出，step=1 的共同落盘状态逐位一致。
- Python 语法检查、git diff --check 通过。

初始 pair 测试中心距 6，支撑域重叠，第 5 步出现跨进程末位差异，保留
`compare_overlap_failed.log` 及 before_pair/after_pair 日志。严格对照改用 pair_disjoint，
未降低逐位判据。重叠情况下的独立算法覆盖仍由 --check 提供。

二进制 .fpd 和构建产物按仓库规则不提交；当前本地留存。`prepare.py` 可重新生成测试初态。
重跑命令（仓库根目录；旧二进制需从上述提交单独构建并放回 /tmp/fpd-v3-before）：

```bash
python3 baseline/ckpt_v4_20261025/prepare.py
ACC_DEVICE_TYPE=nvidia /tmp/fpd-v3-before baseline/ckpt_v4_20261025/wall.cfg
bash baseline/ckpt_v4_20261025/validate_extra.sh
bash baseline/ckpt_v4_20261025/validate_gpu.sh
python3 baseline/ckpt_v4_20261025/interop.py
python3 baseline/ckpt_v4_20261025/compare.py
```

## 待用户执行：VTK

已更新自检（粒子 velocity/unwrapped 存在且 force 不存在）和往返验证逻辑；未运行 pvpython。
使用全新输出目录避免复用已有转换产物：

```bash
PV=/home/doll/Software/ParaView-6.2.0-RC1-MPI-Linux-Python3.12-x86_64
"$PV/bin/pvpython" tools/fpd2vtk.py --self-test
"$PV/bin/pvpython" tools/fpd2vtk.py baseline/ckpt_v4_20261025/after/wall_0000010.fpd --config baseline/ckpt_v4_20261025/wall.cfg -o vis/ckpt_v4_20261025/
"$PV/bin/pvpython" tools/verify_vtk.py vis/ckpt_v4_20261025/wall.pvd baseline/ckpt_v4_20261025/after/wall_0000010.fpd --config baseline/ckpt_v4_20261025/wall.cfg
```

原有用户初态/输出未转换或删除，tools/make_init.py 顶部数据块未改动。新程序不能读取 v3；
用户后续运行前按目标配置重新生成初态。未修改 doc/ 或物理公式。
