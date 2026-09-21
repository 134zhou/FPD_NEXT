# AGENTS.md

本仓库是个人使用的 FPD 流体粒子模拟程序。优先让代码直接、易读、便于修改。
物理背景与历史验证见 [README.md](README.md)，进度与决策见 [PROGRESS.md](PROGRESS.md)。

## 工作方式

- 实现前明确计划、基线和数值验收标准；用户已批准的计划直接执行。
- 发现方案有误，先更新计划，并在 PROGRESS.md 写清原因，再继续实现。
- 每阶段完成后用中文提交；实验配置与日志归档至 baseline/。
- 不以目测或“推导应该如此”代替数值判据；不确定的机制先用 spike 验证。
- 精简以合法输入下结果一致为准，不追求删行比例，不顺手改物理算法。

## 目录与构建

- src/ 与 src/include/：生产实现；fpd_core 只能包含生产代码。
- tests/：fpd_check、fpd_tool；tests/spike/：不参与常规构建的实验。
- tools/：Python 初态生成与离线可视化；config/：运行配置。
- 新验证代码加入 CMakeLists.txt 的 fpd_check 目标。

```bash
cmake -B build && cmake --build build -j
./build/fpd config/smoke.cfg --set n_steps=200
./build/fpd_check --check
./build/fpd_check --check-potential
./build/fpd_check --check-config
./build/fpd_check --check-poisson
./build/fpd_check --check-wall
./build/fpd_check --lambda 32
./build/fpd_tool --dump-ckpt out/run_0000200.fpd
./build/fpd_tool --diff-ckpt before.fpd after.fpd
```

使用 nvc++、OpenACC、CUDA/cuFFT/cuRAND；GPU 为 RTX 3060 12 GB。
GPU 验证设置 `ACC_DEVICE_TYPE=nvidia`，不能把 CPU 回退结果当作 GPU 验证。
不要在运行中的二进制上重建；长跑前先测吞吐，日志用 std::endl 或 fflush 刷新。

## 可信输入与精简原则

- 用户负责数值范围、物理可行性、时间步稳定性、模板盒与截断半径限制。
- 必须提供 init_file；Nx/Ny/Nz 只来自其头部，配置和 --set 显式给出就报错。
- 保留基础必填项和所选势模型的必填参数；基础项和势参数在处理完 --set 后检查。
- 保留语法、数值转换、未知键、势名称与移位名称的基本报错；不做拼写建议。
- 重复键最后赋值生效；已知但当前模型不用的参数允许存在；配置归档输出全部可配置字段。
- 不添加范围/组合检查、参数预测警告、严格/宽松模式或面向通用用户的防御层。
- 保留实际穿墙、NaN/Inf、初态边界一致性、文件完整性、读写失败及 GPU 库错误检查。
- 数值判据验证算法正确性，继续保留；CPU/GPU 独立参考不能互调。
- 不引入新框架或泛化抽象；优先删除重复逻辑，保持合法输入的接口和行为。

## 代码约定

- 中文注释；函数 snake_case，文件 PascalCase.cpp，Allman 括号。
- 小写 vx/fx/eta 表示流体场；大写 Vx/Fx/Rx 表示粒子量。
- 沿用裸指针、POD、手动 OpenACC 映射；不引入 class/智能指针/异常体系。
- GPU 函数的配置 POD 按值传入；acc routine seq 的引用参数会要求设备对象。
- 新增数组参数必须核对每个 present()；遗漏可能静默给出垃圾数据。
- 不在 lambda 内使用 acc update；模板的 acc routine seq 放在 template 之前。
- 多粒子写网格用 atomic；每粒子求和用 gang + vector reduction。
- cuFFT/cuRAND 调用使用 Check.h 包装。

## 数值不变量

- 边界固定为 x/y 周期、z 上下无滑移；不要恢复已删除的全周期路径。
- 索引 `i + j*Nx + k*Nx*Ny`，x 最快；dx、密度、溶剂粘度均为 1。
- Stencil.h 是半格偏移、支撑域截断和模板回绕的唯一实现。
  模板盒 n_range = 2*range+1；投影与归一化使用同一位置和同一权重，各方向独立。
- Wall.h 是 z 棱边存储、法向面钉死和 RNG 长度的唯一实现；棱边为 Nx*Ny*(Nz+1)，
  调用方保证索引范围。Stokes.cpp 按体心、内部 z 棱、两片壁面三个 kernel 计算通量，
  壁面 kernel 使用已经代入 ghost 后的显式 ±2vx/±2vy 公式。
- Poisson.cpp 与 doc/PressurePoisson.md §12 共同定义 div∘grad。
  变换归一化为 Nx*Ny；(0,0) 奇异列用 d_0 -= 1 定规，不能清空整列。
- 速度/噪声边界属于 Wall.h；粒子壁面力属于 Potential.h，二者分开维护。
  壁面力使用表面间隙 h，要求 h>0；不对负间隙 clamp，保留运行中越界中止。
- 壁面剪切噪声幅度 ×√2 只作用于两壁的 EDGE_YZ/EDGE_ZX；公式见 §14。
- Philox 每步按 `(2*step+c)*size*3` 设置 offset；不从存档恢复 RNG 状态。
- 固定粒子数和 GPU 配置下保持力求和顺序；不引入粒子力 atomic 求和。
- 相场支撑域重叠时粘度 atomic 导致跨进程末位差异；逐位重启只承诺不重叠情况。
- 不把旧代码作为定量真值；独立参考不得复制待测实现的模板盒等假设。
- 离散常量使用模拟相同的求和；不能用解析球体积代替扩散界面体积。
- 平衡态误差用分块平均：b≥10τ_int、块数≥20，冲突时报 FAIL，不能降低门槛。

## 按修改内容选择验证

| 修改内容 | 必须执行 |
|---|---|
| Stencil、Viscosity、Force、Velocity、Potential.h | --check；相场加权须覆盖 N=2 重叠均匀流 |
| Poisson.cpp | --check-poisson；同步 PressurePoisson.md §12 |
| Wall.h、Stokes.cpp、Viscosity.cpp | --check-wall；同步 PressurePoisson.md §14 |
| Potential.h/.cpp、Wall.h、main.cpp 越界逻辑 | --check-wall，含 W8；壁面力改动同步 §15 |
| 势公式 | --check-potential，保持能量与导数独立实现 |
| 配置/主流程精简 | --check-config、合法输入检查点前后对照、单粒子重启 |
| 二进制 IO 任一侧 | 下述 C++/Python 互操作检查 |

--check-wall 覆盖 W1–W6、W3b、W4s、W5b 及 W8；W3b 直接核对两片壁面剪切通量，
W8e 是近壁 GPU/CPU 力装配对照。
关闭壁面势时的回归锚是原壁面判据输出不变，不能用 W8e 单项替代。
粒子侧壁面 FDT/冻结粒子能量均分尚未完成；历史 W7 重启记录不代表它已验证。
已删除的 --noise/--equipart 不属于现有覆盖；详见 PROGRESS.md 已知隐患。

## IO 与可视化

C++ 只读写 .fpd；init_file 必填，同一文件格式用于初态和重启，文件 step 决定起点。
n_steps 是绝对目标步数。网格尺寸、N、step 来自文件，其余参数来自配置；逐位重启须保持同一配置和 seed。
格式由 src/include/IOBin.h 与 tools/fpd_format.py 共同定义，当前版本为 3；v2 明确拒绝。
头部只有尺寸、N、step 等布局元数据；.config.used 不含尺寸键，离线工具用 --config 读取物理量。
numpy 形状必须为 (Nz, Ny, Nx)；版本、字段顺序与 FNV-1a 两侧同步。

```bash
python3 tools/make_init.py --grid 8 4 2 --pattern index --empty -o /tmp/t.fpd
./build/fpd_tool --dump-ckpt /tmp/t.fpd --at 3 1 1  # 必须为 3001001
PV=/home/doll/Software/ParaView-6.2.0-RC1-MPI-Linux-Python3.12-x86_64
$PV/bin/pvpython tools/fpd2vtk.py --self-test
$PV/bin/pvpython tools/fpd2vtk.py out/run_*.fpd --config config/production.cfg -o vis/
```

make_init.py 仅依赖 stdlib；转换和 verify_vtk.py 必须使用 pvpython。
.fpd 不保存压力，压力仅在求解器内部及数值判据中使用；沉降统计需 --config。
用户自行转换可视化，提供命令即可；ParaView 打开生成的 .pvd。
沉降数据选择 out/sed_[0-9]*.fpd，避免包含或删除 sed_init.fpd。
