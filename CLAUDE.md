# CLAUDE.md

给在本仓库工作的 Claude Code 实例的指引。

本文只写**怎么在这里干活**。物理背景、缺陷清单、验证结果见 `README.md`；
进度与决策记录见 `PROGRESS.md`。

## 构建与测试

**代码布局**：`src/` = 生产代码（只有 `fpd`），`tests/` = 全部验证与测试代码
（`fpd_check`、`fpd_tool`、以及不参与构建的 `tests/spike/`）。
`fpd_core` 静态库里**只有 src/** —— 加新文件时先想清楚它属于哪一侧。

```bash
cmake -B build && cmake --build build -j     # nvc++ 26.3 + CUDA 13.1 + CMake 4.2.3
                                              # 产出 build/fpd / fpd_check / fpd_tool

./build/fpd config/smoke.cfg [--set k=v ...]  # 生产运行；--set 覆盖配置项
./build/fpd_check --check                     # 自检：力守恒 / 亚格点 / 混叠 / N=2 重叠 / 力链路
./build/fpd_check --check-potential           # 势函数自检（纯 CPU，无卡可跑）
./build/fpd_check --check-wall                # z 向壁面 W1–W6 + W5'a/W5b（需 GPU）
./build/fpd_check --lambda [L]                # 常量表 λ^T, M_i（纯 CPU，秒级）
./build/fpd_check --noise [L] [dt] [kT] [steps]   # 测试 3A：纯流体噪声谱
./build/fpd_check --equipart <ghost|frozen|moving> [L] [dt] [kT] [steps] [seed]   # 测试 3B

./build/fpd_tool --dump-ckpt <f> [--at i j k] # 看 .fpd 头部 / 指定格点的值
./build/fpd_tool --diff-ckpt <a> <b>          # 逐位比较两个检查点

# Python 侧（见下方「IO 与可视化」）
python3 tools/make_init.py --grid 128 64 32 --single -o out/init.fpd
$PV/bin/pvpython tools/fpd2vtk.py out/*.fpd -o vis/

cmake --build build -j 2>&1 | grep -i error   # 编译错误
nvc++ -acc -O3 -Isrc/include -Minfo=accel -c src/X.cpp -o /tmp/x.o   # 看 kernel 生成
```

## IO 与可视化

**C++ 只写二进制 `.fpd`，绝不产出可视化格式。** 初始构型、检查点、重启文件都是同一个格式，
文件里的 `step` 决定从哪继续 —— 所以 C++ 除配置文件外没有任何文本解析。
格式规范在 `src/include/IOBin.h`，Python 侧镜像实现在 `tools/fpd_format.py`。

| 脚本 | 依赖 | 跑在哪 |
|---|---|---|
| `tools/make_init.py` | **纯 stdlib** | 任何 python3（含裸 conda base） |
| `tools/fpd2vtk.py` | numpy + vtk | **必须 pvpython**（系统 python3 无 numpy） |
| `tools/verify_vtk.py` | numpy + vtk + paraview | 必须 pvpython |

`PV=/home/doll/Software/ParaView-6.2.0-RC1-MPI-Linux-Python3.12-x86_64`

转换后在 ParaView 里**打开一个 `vis/<name>.pvd`** 就得到完整时间动画 + 流体/粒子双 block。
`fpd2vtk.py` 不手写 XML，全部交给 VTK 自己的 writer。

> ⚠️ **C++ 与 Python 共同定义同一个二进制格式**，这是新的漂移点。改动任一侧后必须重跑：
> ```bash
> python3 tools/make_init.py --grid 8 4 2 --pattern index --empty -o /tmp/t.fpd
> ./build/fpd_tool --dump-ckpt /tmp/t.fpd --at 3 1 1     # 必须报 3001001
> $PV/bin/pvpython tools/fpd2vtk.py --self-test          # 轴序与插值方向
> ```
> `IDX = i + j*Nx + k*Nx*Ny` 意味着 numpy 必须 reshape 成 `(Nz, Ny, Nx)`，**x 是最后一维**。
> 搞反不报错，只得到一个转置的场。FNV-1a 校验和同样是两处独立实现。

`init_file` 与「重启文件」是同一个 key，文件里的 `step` 决定从哪继续。
`n_steps` 是**绝对目标步数**，不是增量。

## 环境与机时

RTX 3060 12 GB。实测吞吐：**153 步/秒**（128×64×32）、**1128 步/秒**（32³）。
注意 32³ 只有 1/8 的格点却只快 7.4 倍 —— 那个规模已接近**启动延迟主导**，
「小盒子会快很多」的直觉不成立。排期前先按实测数算机时；
Phase 4 的 L=192 立方盒是 27 倍网格量。检查点写入 3.9 ms/次，可忽略。

> ⚠️ **别在后台运行时重建** —— 正在执行的二进制被写会触发 `ETXTBSY`。
> 长跑输出用 `std::endl` 或 `fflush`，否则重定向到文件时会被缓冲，只能盲等。
> 长跑前先小规模测吞吐，别盲目启动几小时的任务。

## 不可违反的约束

### 1. `Stencil.h` 是交错网格的唯一真值源

`src/include/Stencil.h` 的 `stencil_point<Loc L>()` 是**唯一**实现半格偏移约定、
球形截断判据、周期回绕的地方。`Viscosity` / `Force` / `Velocity` 三个模块必须共用它。

**不要在任何地方重新手写模板盒遍历循环。** 缺陷 C1/C2 的成因正是这段循环被复制粘贴了
三份然后各自漂移 —— 力投影用 x-face 的权重、归一化却用 cell-center 的求和。

规则一句话：**做投影用的权重函数、和做归一化的求和，必须是同一个函数在同一组点上的求和。
三个方向各一套，不能共用。**

**`n_range` 必须是 `2*range + 1`，不能是 `2*range`。** 盒取 `[kn-range, kn+range]`。
半径 `range` 的球心落在格胞任意位置时能触及的整数层最多 `2*range+1` 个，窄一号的盒
会在负方向**静默漏掉一层**（Phase 7-B 修掉的缺陷，复现脚本 `tests/spike/spike_stencil_box.py`）。
漏掉的权重随亚格点位置与 `Loc` 变化 ⇒ 方向相关的伪偏差，且**力守恒判据抓不到**
（分子分母一起漏，比值自洽）。判它的手段只有两个：暴力枚举全空间对照，或与连续积分对照。

### 2. 改完必须重跑 `--check`

动过 `Stencil.h` / `Viscosity.cpp` / `Force.cpp` / `Velocity.cpp` / `Potential.h`
之后**必须**跑 `./build/fpd_check --check`，全绿才算数。`stencil_point` 是单点故障，
它错了三个模块一起错。

### 2b. 随机数流按 step 定位，**不要改成累计记账**

`Stokes.cpp` 每步生成前显式 `curandSetGeneratorOffset` 到
`offset(step,c) = (2*step+c) * size*3`。这看起来绕，但是实测逼出来的：

- 只有 **Philox** 支持有意义的 offset；XORWOW（原来的 `CURAND_RNG_PSEUDO_DEFAULT`）
  和 MRG32K3A 都对不上，**无法逐位重启**
- **「生成 n 个数后序列前进 n」是错的**：`offset=1*n` 能对上第 2 批，
  但 `2n`、`3n` 都对不上连续生成的对应批次。所以任何累计记账都会失配
- 但同一 offset 生成同样数量必定逐位可复现，不同 slot 之间零共同值、
  互相关 ~1/√n（统计独立）

证据在 `tests/spike/spike_rng_offset.cpp`、`spike_rng_slice.cpp`、`spike_rng_advance.cpp`。
`.fpd` 里的 `rng_draws` 字段**仅供人读**，读回时不采信，offset 一律由 step 重算。

### 3. 两个测试盲区（都真实骗过人）

- **力守恒判据对 `stencil_point` 的内部公式错误是盲的** —— 投影和归一化会一起错，
  `Σf == F` 照样成立。需要独立数值对照，见 `tests/spike/spike_template_routine.cpp`
  （用独立写的 Python 复核过 `∫φ`）。⚠️ 那份 Python 参考**当时复制了同一个模板盒
  尺寸**，所以两处一起错、没能发现模板盒缺陷（Phase 7-B 才发现，见约束 1）。
  **共享假设的参考实现不是独立参考** —— 真正的判决者是连续球坐标积分与暴力枚举。
- **单粒子测试对一整类缺陷是盲的** —— C2 在 N=1 时数学上恰好抵消，骗过了 20 万步模拟。
  **任何涉及相场加权的改动，都必须补 N=2 重叠粒子的用例。**
  判决性判据：均匀流场 `v ≡ 1` 下 `V_i` 必须精确等于 1，与重叠无关。

### 4. 旧代码不是真值

`doc/old_code/` 仅供参考。已确认它在四处是错的或有瑕疵（cellList 的 `nnIndex[12]`、
cell 数不足 3、z 向边界 hack、重力硬编码），详见 `README.md`。
「和旧代码对齐」只能作定性 sanity check，**定量真值来自文献公式和自洽性测试**。

### 5. 粒子间力的求和顺序必须确定（逐位重启的生命线）

`Potential.cpp` 的 `compute_particle_forces` 是 GPU 全矩阵（gang-over-i +
vector-reduction-over-j，零 atomic），求和顺序在固定 N 与固定 gang/vector 配置下是
确定的 —— 与 `Viscosity.cpp` 算 `sum_phi`、`Velocity.cpp` 算 `V` 是同一机制
（`spike_pair_bench` 已实测逐位可复现）。

**但 FPD 的逐位重启只在粒子支撑域不重叠时成立**（间距 > 2·range = 16.8）。
重叠时 `eta[ijk] += d_eta*w` 的 atomic 累加顺序跨进程不确定（~1e-16 差异），
这是 Viscosity.cpp 的既有行为，不是粒子间力引入的。而「粒子间力」需要粒子靠近
（间距 < potential rcut ≤ 15），靠近必然重叠 —— 所以「有力 + 逐位重启」在当前盒子下
**原理上不可兼得**。逐位重启的常驻回归用 `smoke.cfg`（N=1 或 N=2 不重叠）。

### 6. 粒子间力的 CPU/GPU 对照是判据，不是生产代码

`compute_particle_forces_cpu`（串行 i<j 半矩阵）只用于判据对照（`--check-potential`
的 N=3、`--check` 的 J5、`fpd_tool --verify-forces`），生产热路径走 GPU 版。
两份实现必须独立维护（不互调），否则「GPU vs CPU 对照」判据退化成自洽的、盲的。

### 7. 改 `Poisson.cpp` 必须同步改 `doc/PressurePoisson.md` 并重跑 `--check-poisson`

`src/Poisson.cpp` 与 `doc/PressurePoisson.md` **共同定义** `div∘grad` 算子（Phase 7-A
引入的新漂移点，与 C++/Python 共同定义 `.fpd` 是同一类风险）。改求解器后必须：
同步更新文档 §12 的公式 ↔ 代码对照表，并重跑 `./build/fpd_check --check-poisson`。
归一化因子（周期 `size`、壁面 `Nx·Ny`）在 `solve_pressure` 里**只出现一次**，别散落。
壁面 `(0,0)` 奇异列用 `d_0 -= 1` 定规，**不许整列清零**（那会删掉支撑粒子重量的静压）。

### 8. `Wall.h` 是 z 向边界访问的唯一真值源

`src/include/Wall.h` 是**唯一**实现 z 向邻居访问、壁面 ghost、棱边存储映射、
壁面噪声因子的地方（Phase 7-B 引入的新漂移点，与 `Stencil.h`、`Poisson.cpp`
是同一类风险）。

**不要在 `Stokes.cpp` 里重新手写 `k` 或 `k-1` 的 z 向访问。** 逻辑 `k` 的物理含义
在周期/壁面两种模式下完全相同，只有存储映射不同 —— 手写一次就破坏这个不变量。

改完 `Wall.h` / `Stokes.cpp` / `Viscosity.cpp` 后**必须**跑
`./build/fpd_check --check-wall`（含 W1/W2/W3/W4/W6/W5），并同步
`doc/PressurePoisson.md` §14（公式 ↔ 代码对照表在那里）。

几个不能忘的数：
- 壁面剪切噪声因子 **γ=2（幅度 ×√2）**，只有 `EDGE_YZ`/`EDGE_ZX` 的两片壁面棱边需要。
  `Π_zz` 与 `Π_xy` **不需要**。推导见 §14.3，判据 W5b 用 γ≡1 对照证明它被数据选中
- `dim = n_v − size + 1`（能量均分的目标 `⟨Σ|v|²⟩ = kT·dim`），周期/壁面通用
- 棱边数组在壁面模式下是 `Nx·Ny·(Nz+1)`

## 代码约定

- **命名**：小写 = 流体场（`vx`、`fx`、`eta`），大写 = 粒子量（`Vx`、`Fx`、`Rx`）。贯穿全代码，别破坏
- **注释用中文**，物理术语可混排英文。函数 `snake_case`，文件名 `PascalCase.cpp`，Allman 括号
- **索引** `IDX(i,j,k) = i + j*Nx + k*Nx*Ny`，x 最快。与 `cufftPlan3d(Nz,Ny,Nx)` 匹配，**别动**
- **单位**：格距 `dx ≡ 1`，密度 `ρ ≡ 1`，溶剂粘度 `η_ℓ ≡ 1`
- 风格接近「带 vector 的 C」：裸指针 + OpenACC 手动映射，不用 class/智能指针/异常。
  新代码沿用，别引入现代 C++ 抽象（`Stencil.h` 的模板是唯一例外，且已验证不影响 kernel 生成）
- **放新文件**：生产代码进 `src/`（头文件进 `src/include/`），验证代码进 `tests/`。
  `tests/` 里加文件要在 `CMakeLists.txt` 的 `fpd_check` 目标里登记
- POD 结构体（`NS_Config`、`PhiParams`、未来的 `PotentialParams`）**按值传进 GPU 函数**。
  `acc routine seq` 里用引用参数会要求对象位于设备内存

## nvc++ / OpenACC 的坑

- **`present()` 遗漏不报错，只给垃圾数据**。nvc++ 对未列出的裸指针走隐式查找，通常「碰巧能跑」。
  每次给函数加数组参数，必须同步检查所有 `present()`
- **`#pragma acc update` 不能放在 lambda 里** —— 闭包捕获让运行时查不到设备映射，
  运行时报 `data in update host clause was not found on device`。用普通函数或直接内联
- `#pragma acc routine seq` 作用于**模板函数**可行（已验证），但 pragma 要放在 `template<...>` **之前**
- 写到**网格**上必须用 `atomic`（多粒子支撑域会重叠）；写到**每粒子标量**上用
  `parallel loop gang` + 内层 `loop vector reduction`，不要用 atomic 打同一地址
- cuFFT/cuRAND 调用一律用 `Check.h` 的 `CUFFT_CHECK` / `CURAND_CHECK` 包住
- 重定向输出时 `std::cout << "\n"` 不刷新，长跑要用 `std::endl`，否则盲等

## 数值定则

- **稳定性**：显式粘性项要求 `dt ≲ ρdx²/(2d·η_c)`。η_c=50、d=3 时是 `dt < 0.0033`
- **噪声的 O(dt) 偏差**：纯流体 η=1 实测 `bias(%) ≈ 335 × dt`，要 kT 精度 < 1% 需 `dt < 0.0033`。
  **粒子侧（η_c=50）已由测试 3B 证明不会按扩散数放大 50 倍** ——
  实测 +0.5~0.7%，与流体式预期 0.67% 吻合。原因是粒子速度是 φ 支撑域上的低通平均，
  高 k 被滤掉，而偏差 `~ν k² dt/2` 对高 k 最大。**生产可用 dt = 0.002**
- **平衡态统计的误差棒必须来自分块平均**。把相关样本当独立样本会低估约 7.5 倍。
  三条硬约束：`b ≥ 10τ_int`、`n_b ≥ 20`、冲突则**报 FAIL 要求延长，不许降标准接受**。
  采样间隔要在**物理时间**上固定（`samp_every = 0.5/dt`），否则 dt 越小样本越相关
- **`Σφ_α²` 三方向一致性是交错混叠的判据**（`Σφ_α` 只是 k=0 分量）。已并入 `--check`
- **平衡时间**：最慢模式 `τ = (L/2π)²/ν`。做平衡态统计时盒子越大越慢，
  纯流体测试务必用小立方盒（L=128 要 41 万步，L=32 只要 2.6 万步）
- **λ^T、M_i 等常量必须用与模拟相同的离散求和**（`./build/fpd_check --lambda`）。
  理由是**自洽性**，不是格点不准 —— 格点其实精确到 5e-6。
  但**绝不能用解析球体 `(4/3)πa³`**：真值比它大 24.1%，那是扩散界面的曲率项
  （`+8πaξ²π²/24`），误用会让 `M_i` 差 24%，直接毁掉验证
- **势截断半径必须 `rcut < min(Nx,Ny,Nz)/2`**（严格小于，等号会双重计数周期镜像）。
  旧代码用 `cutoff = 3·Re = 22.2` 在它自己的 128×128×64 里合法，但**照搬到
  `production.cfg` 的 128×64×32（min/2=16）就违反**。`validate_config` 已硬性拦截。
  势残余力比 `|U'(rcut)|/max|U'| > 1e-6`（峰值取物理可达区 `[2a, rcut]`）会报警

## 粒子间势与外场

势函数在 `src/include/Potential.h`（enum + POD + `routine seq` 的 switch，无函数指针）。
支持 `none / wca / morse / lj126`；配置项 `potential` + `pot_eps/pot_sigma/pot_De/
pot_alpha/pot_r_eq/pot_rcut/pot_shift`。

- **WCA 不是独立分支**：WCA ≡ LJ12-6 + `rcut = 2^{1/6}σ`（派生）+ 移位，
  且此时能量移位与力移位恒等（U'(rc)=0）。switch 只有 3 个真实分支。
- **默认 `shift=energy`**（U 移位保证 U(rc)=0，力不变）。`force` 会改阱深
  （Morse α=1, rcut=15 时 50→49.57），只在残余力比不可接受时用。
- **参数查表**：`validate_config` 检查「用不到的势参数报错、需要的必须给」，不给默认值
  兜底；启动时打印势的完整解释（含残余力比）。Morse 宽度参数叫 `pot_alpha` 不叫 `a`
  （`a` 是粒子半径，README.md:21 的相场公式里就是）。
- **外场** `gravity_x/y/z`（每粒子恒力）+ `gravity_compensate`（默认 1）：外场均匀时
  Σ_n F = N·g ≠ 0，流体 k=0 被线性加速、整盒漂移，`bg = −ΣF/size` 抵消它。
  **z 向无滑移壁面（Phase 7）上线后此开关默认转 0**（壁面提供真实动量汇）。

`pair_energy_bare` 与 `pair_force_over_r_bare` 是**两份独立实现**（势 vs 解析导数），
`--check-potential` 用四阶数值微分 + Python 黄金表交叉检验；`min_image` 在
`Potential.h` 里是唯一真值源。势参数**不进 `.fpd` header**（Python 侧用不到），
自洽性由 `fpd_tool --verify-forces` 判据兜底。

## 当前进度（详见 `PROGRESS.md`）

Phase 0/1/2/3/5 完成，**Phase 7 完成**（7-A 求解器 + 7-B 接线）。流体求解器、
交错网格封装、噪声标定（3A + 3B）、配置系统、逐位断点重启、
**粒子间相互作用力（WCA/Morse/LJ126 + 外场）**、
**z 向无滑移壁面（`boundary_z = noslip`：切向 ghost + 壁面剪切噪声 ×√2 +
相场按 `Loc` 截断 + 棱边数组 Nz+1 层）**均已通过数值判据。
`FpdState` 已实现，自检拆成独立可执行 `fpd_check`。

**下一步是 Phase 6（生产算例 + 旧代码对照）或 Phase 4（L 外推）**，无阻塞项。
Phase 7-B 明确**未做**粒子-壁面排斥势（粒子靠初始构型远离壁面，越界报错中止）。

## 工作方式

- 每个阶段完成后 git 提交，commit message 用中文
- 实验日志归档到 `baseline/`
- **判据必须是数值判据**，不是「看起来对」。本项目已经有三次「实测推翻了计划里的假设」
  （模板 `routine seq` 可行性、k=0 修正系数、RNG offset 语义），
  凡是「文档说」「推导应该」的地方都先写个 spike 测掉
- 发现计划有误时**先改计划文件再动手**，并把「为什么原方案不成立」写进 `PROGRESS.md`
