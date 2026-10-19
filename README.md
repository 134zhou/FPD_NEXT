# FPD_NEXT

胶体悬浮液的**流体粒子动力学**（Fluid Particle Dynamics, FPD）GPU 模拟器。

刚性胶体粒子表示为高粘度的流体区域（相场 φ），流体力学相互作用由此自然涌现，
无需显式追踪移动边界。实现为 OpenACC + cuFFT + cuRAND，单 GPU，双精度。

边界条件固定：**x/y 周期，z 上下两片无滑移硬壁**（不是配置项）。
进度、决策与历史记录见 [PROGRESS.md](PROGRESS.md)，工程约定见 [AGENTS.md](AGENTS.md)。

## 参考文献与推导笔记

- **[PRL]** Tanaka & Araki, *Simulation Method of Colloidal Suspensions with Hydrodynamic
  Interactions*, PRL **85**, 1338 (2000)
- **[SM]** Furukawa, Tateno & Tanaka, *Physical foundation of the fluid particle dynamics
  method for colloid dynamics simulation*, Soft Matter **14**, 3738 (2018)

两篇的 PDF 与 md 都在 `doc/`。本实现以 **[SM] 2018** 为准（[PRL] 是简化版：无噪声、
无转动、力密度未归一化、丢对流项），并在其基础上扩展了**变密度**（见下）。
转动自由度未实现。

`doc/` 下其余文件是推导笔记，**原理与代码一致，但编号、实测数值、接口名等细节可能
滞后于代码，以代码和 `fpd_check` 输出为准**：

| 文件 | 内容 |
|---|---|
| `doc/ns-equation.md` | 从 F=ma 到 FPD 形式的动量方程（含变密度） |
| `doc/fourier-projection.md` | 投影法、2D FFT 求解、变密度 PCG |
| `doc/discretization.md` | MAC 网格、模板盒、z 向三对角的推导 |
| `doc/PressurePoisson.md` | 压力泊松与壁面的完整推导（§12 算子对照、§14 壁面、§15 壁面势、§16 变系数） |
| `doc/PoissonCodeWalkthrough.md` | `solve_pressure` 代码阅读记录 |
| `doc/Stokes_line_by_line.md` | `Stokes.cpp` 逐行解释（常密度分支） |
| `doc/old_code/` | 上一代 OpenMP 实现，**仅供参考，不是真值** |

## 控制方程

```
相场      φ_i(r) = ½[tanh((a - |r-r_i|)/ξ) + 1]          球形截断半径 range = ⌊2(a+ξ)⌋
粘度      η(r)   = 1 + (ratio_eta - 1) Σ_i φ_i(r)         η_l ≡ 1，ratio_eta = η_c/η_l
密度      ρ(r)   = 1 + (ratio_rho - 1) Σ_i φ_i(r)         ρ_l ≡ 1，ratio_rho = ρ_p/ρ_l
动量      ρ(∂_t v + v·∇v) = -∇p + ∇·[η(∇v+∇vᵀ)] + ∇·σ̃ + f + (ρ-1)g
不可压    ∇·v = 0
力密度    f_α(r) = Σ_i F_iα φ_i^α(r) / S_i^α,   S_i^α = Σ_r φ_i^α(r)
粒子速度  V_iα   = Σ_r v_α φ_i^α / S_i^α,       dR_i/dt = V_i
```

- `φ_i^α` 是同一个相场在 α 方向面心上的采样（三套错开的格点），`S_i^α` 是该方向的离散积分；
  投影与平均共用同一组位置和权重，因此 `Σ_r f_α = Σ_i F_iα` 恒成立。
- `σ̃` 是涨落应力：对角分量幅度 `√(2η)·W`，非对角 `√η·W`，`W = √(2kT/dt)`（`noise_on=0` 时 W=0）。
- 重力的均匀部分 `ρ_l g` 被静压吸收，只剩净体力 `(ρ-1)g`，粒子净重力自动是
  `(ρ_p-ρ_l)·∫φ·g`。`gravity_*` 是**加速度**：早期版本里它是直接施加在每个粒子上的力，
  旧配置按 `g = F_old/((ratio_rho-1)·∫φ)` 换算（∫φ = 170.306，例如 `-10 → -0.058718`）。
- 代码里没有「粒子质量」：惯性全部由流体密度场 ρ(r) 承担，粒子速度只是流场的相场加权平均。
- 论文 [SM] 式 (3) 假设液体与粒子等密度；`ratio_rho = 1` 时本实现精确退回该形式。

### 一步时间推进（`src/main.cpp`）

```
compute_particle_forces      F = 粒子对势 + 壁面势（重力不在这里）
[写检查点 / 打日志]           按 interval_ckpt / interval_log，含越界、NaN 检查
update_viscosity_fields      η（体心）、η_XY/YZ/ZX（棱心）、S_i^x/y/z
update_mass_fields           面心 a = 1/ρ（仅 ratio_rho ≠ 1）
update_force_field           f ← F（面心）
apply_gravity                f += (1/a - 1)·g（仅 ratio_rho ≠ 1 且 g ≠ 0）
step_navier_stokes           噪声 → Π → v* → b = ∇·v*/dt → 解压力 → 投影
update_particle_velocity     V ← v
update_particle_position     R += dt·V（x/y 回绕；z 不折叠；Ru 为不折叠坐标）
[PCG 收敛检查]                每步读主机侧 diag[1]，未收敛则写出 step+1 并中止（rc=3）
```

`step_navier_stokes`（`src/Stokes.cpp`）的显式步与投影：

```
ratio_rho = 1:   v* = v + dt·(-∇·Π + f)                            Π = vv - η(∇v+∇vᵀ) - σ̃
ratio_rho ≠ 1:   v* = v + dt·[ -∇·(vv) + a·(-∇·(Π - vv) + f) ]      a = 1/ρ（面心）
投影:            ∇·(a∇p) = ∇·v*/dt,     v = v* - dt·a·∇p
```

**1/ρ 只乘力密度（粘性、噪声、外力、压力梯度），不乘对流加速度**。代码里总通量散度
`d_pix` 被拆成对流部分 `adv`（由 `conv_*` 唯一实现）与 `d_pix - adv`。
`ratio_rho = 1` 时走合并写法与常系数直解。

## 空间离散

| 量 | 位置 | 存储 |
|---|---|---|
| `p`、`η`、`Π_xx/yy/zz` | 体心 | `IDX(i,j,k) = i + j*Nx + k*Nx*Ny`，x 最快 |
| `vx/fx`、`vy/fy`、`vz/fz`、`a_x/y/z` | 面心 | `IDX`；`vz` 的 `k=Nz-1` 是上壁，恒 0 |
| `η_XY`、`Π_xy` | xy 棱 | `IDX` |
| `η_YZ/ZX`、`Π_yz/zx` | yz/zx 棱 | `Nx*Ny*(Nz+1)`，逻辑 `k ∈ [-1, Nz-1]`，两片壁面各占一层 |

格距、溶剂粘度、溶剂密度均为 1。几个单一真值源：

- **`src/include/Stencil.h`**：`stencil_point<Loc>()` 是半格偏移、球形截断、x/y 回绕和
  z 向截断的唯一实现。模板盒边长 `n_range = 2*range+1`（盒取 `[kn-range, kn+range]`，
  窄一格会在负方向静默漏层）。Viscosity / Mass / Force / Velocity 共用它。
- **`src/include/Wall.h`**：z 棱边存储映射、法向面钉死（`v*z`、`vz` 在上壁为 0；下壁不在数组里）
  和 RNG slot 长度。
- **`src/Stokes.cpp`** 的应力装配分三个 kernel：体心与 XY 棱（`k∈[0,Nz-1]`）、内部 YZ/ZX 棱
  （`k∈[0,Nz-2]`）、两片壁面棱（`k=-1, Nz-1`）。切向无滑移由 ghost `vα[-1]=-vα[0]`、
  `vα[Nz]=-vα[Nz-1]` 代入得到：壁面对流为 0、粘性项为 `∓2η vα`；
  壁面棱的剪切噪声幅度 ×√2（控制体只有一半，推导见 `doc/PressurePoisson.md` §14.3）。

### 压力求解（`src/Poisson.cpp`）

- **常系数**（`ratio_rho = 1`）：xy 向批量 2D FFT（`cufftPlanMany`，batch=Nz）+ 每个 `(kx,ky)`
  在 z 向解三对角（Thomas，前推系数 `tri_w` 预算一次）。算子是修正步所用梯度与散度的
  **精确复合** `div∘grad`；z 向 Neumann 由「壁面法向面不修正」导出。`(0,0)` 奇异列先做相容投影
  （减均值，`|Σb̂|` 记入 `diag[0]`），再以 `d_0 -= 1` 定规，不能整列清零。逆变换归一化为 `Nx*Ny`。
- **变系数**（`ratio_rho ≠ 1`，默认 `pressure_solver = pcg`）：`apply_laplacian_var` 给出同构离散的
  `div(a∇p)`；`solve_pressure_var` 用 PCG 求解，预条件器为常系数求解器 M（乘 1/ā，标量缩放不改变
  迭代序列）。初值固定 `p0 = 0`（保证逐位重启，不做 warm start），相对残差容差 `pcg_tol`
  （默认 `1e-12`），上限 200 次。`diag[1]` 是递推残差，`diag[2]` 是投影后从速度场直接量的真残差
  `‖div vⁿ⁺¹‖/‖div v*‖`；此路径的 `∇·v=0` 只到容差，常系数路径为机器精度。
- **PCG 为什么已接近最优**：对任意 x，`xᵀAx = Σ_f a_f(∇x)_f²`，所以预条件后的条件数
  `κ ≤ ρ_max/ρ_min`，**与网格尺寸和粒子数无关**。CG 误差界 `2q^k`、`q=(√κ-1)/(√κ+1)`：
  `ratio_rho=2` 时每次迭代缩小约 0.17 倍，到 1e-12 约 16 次——32³ 到 512×512×128 实测都是 16 次。
  多重网格同样是 O(1) 次迭代，每次不比一次 FFT 便宜；只有 `ratio_rho` 很大（κ≳10）时才值得考虑。
  所以真正的杠杆是**容差**，见下面「压力求解精度」。
- **压力分裂法**（`pressure_solver = split`，可选）：Dodd & Ferrante (2014) 的外推法，
  `∇²pⁿ⁺¹ = ρ₀b + ∇·[(1-ρ₀a)∇p̂]`、`vⁿ⁺¹ = v* - dt[a₀∇pⁿ⁺¹ + (a-a₀)∇p̂]`，`p̂ = 2pⁿ-pⁿ⁻¹`，
  `ρ₀ = min ρ`。每步一次常系数直解，`∇·v=0` 到机器精度，前 2 步用 PCG 建立压力历史。
  **有热噪声时不推荐**：压力在相邻步之间几乎不相关，外推误差集中在粒子内部（见「压力求解精度」）。

## 粒子、势与壁面

- **粒子间势**（`src/include/Potential.h`）：`none / wca / morse / lj126`，enum + POD + switch。
  WCA ≡ LJ12-6 + 派生截断 `rcut = 2^{1/6}σ`；截断移位 `pot_shift = none | energy | force`，默认
  `energy`。`pair_energy_bare` 与 `pair_force_over_r_bare` 两份独立实现（供 `F=-dU/dr` 判据）。
  x/y 最小镜像，z 不做镜像；使用者须保证 `rcut < min(Nx,Ny)/2`。
- **力装配**：GPU 版 gang-over-i + vector reduction over j 全矩阵（零 atomic，逐位确定），
  CPU 版 i<j 半矩阵，两份独立维护。O(N²)，无 cell list。
- **壁面势**：`wallpotential` 与粒子间势同一族函数，但自变量是**表面间隙**
  `h_bottom = (z+½) - a`、`h_top = (Nz-½-z) - a`，参数 `wall_*` 独立。这是让粒子停在可控间隙上的
  模型，不是第一性原理推导。要求 `h > 0`，不 clamp。
- **越界**：z 从不折叠。主程序在检查点/日志步检查 `Rz`：壁面势关闭时要求 `[-½, Nz-½]`，
  开启时要求 `[-½+a, Nz-½-a]`；越界即落盘并以退出码 3 中止。
- **质量场**（`src/Mass.cpp`）：用与 `S_i^α` 相同的 `stencil_point<FACE_*>` 先建 ρ 再取倒数得
  `a = 1/ρ`；`apply_gravity` 是重力的唯一入口。密度是派生量，不入存档。

## 构建与运行

```bash
cmake -B build && cmake --build build -j            # 默认只产出 build/fpd
cmake --build build --target fpd_check -j           # 验证程序
cmake --build build --target fpd_tool -j            # 纯 CPU 检查工具
./build/fpd config/smoke.cfg [--set key=value ...]
```

环境：nvc++ 26.3（NVIDIA HPC SDK）、CUDA 13.1、CMake 4.2.3、RTX 3060 12 GB。
GPU 验证须设 `ACC_DEVICE_TYPE=nvidia`。

| 可执行 | 来源 | 说明 |
|---|---|---|
| `build/fpd` | `src/main.cpp` | 生产程序（需 GPU） |
| `build/fpd_check` | `tests/Check*.cpp`、`tests/Analysis.cpp` | 数值判据（`--check-config/-potential/-tridiag`、`--lambda` 纯 CPU，其余需 GPU） |
| `build/fpd_tool` | `tests/tool_main.cpp`、`tests/Analysis.cpp` | 纯 CPU：dump/diff 检查点、重算力、沉降统计 |

`fpd_core` 静态库只含 `src/` 生产代码，三个可执行共用它。`tests/spike/` 是一次性风险验证程序，
不参与构建。退出码：0 完成；2 配置/IO 错误；3 物理异常中止（NaN/Inf、穿墙、PCG 未收敛，均已落盘）。

## 配置

`key = value`，`#` 起注释；`--set key=value` 在文件之后应用，重复键最后一次生效。
`./build/fpd` 不带参数会打印全部字段。

| 类别 | 键 |
|---|---|
| 必填 | `init_file`、`dt`、`n_steps`（**绝对**目标步数）、`kT`、`radius`、`xi`、`ratio_eta`、`ratio_rho` |
| 可选物理 | `gravity_x/y/z`（加速度，默认 0）、`noise_on`（默认 1） |
| 压力求解 | `pressure_solver`（`pcg`/`split`，默认 pcg）、`pcg_tol`（默认 1e-12；有噪声时推荐 1e-6） |
| 粒子间势 | `potential`、`pot_eps/sigma/De/alpha/r_eq/rcut/shift` |
| 壁面势 | `wallpotential`、`wall_eps/sigma/De/alpha/r_eq/rcut/shift` |
| IO | `out_dir`（out）、`run_name`（run）、`interval_ckpt`、`interval_log`（各 10000）、`seed`（1234） |

所选势的参数必须显式给出：wca 要 `eps, sigma`；morse 要 `De, alpha, r_eq, rcut`；
lj126 要 `eps, sigma, rcut`（壁面势同理，前缀 `wall_`）。`Nx/Ny/Nz` 只来自 `init_file` 头部，
在配置或 `--set` 中出现即报错。

配置层只检查语法、数值转换、未知键、必填项和势/移位/压力求解器名称；参数范围、时间步稳定性、模板盒与截断
限制由使用者保证（显式粘性项大致要求 `dt ≲ min(ρ)/(6η_c)`，`ratio_rho < 1` 须相应减小 dt）。
每次运行写出 `<out_dir>/<run_name>.config.used`：全部字段（不含尺寸键）+ 以 `#` 注释的派生量，
可直接作为输入重跑。

| 配置 | 用途 |
|---|---|
| `config/smoke.cfg` | 5 步小算例，回归与重启验证（`ratio_rho=1`） |
| `config/production.cfg` | 128×64×32 单粒子（`ratio_rho=1`） |
| `config/sed.cfg` | 128×64×32、N=95 集体沉降（`ratio_rho=2`，WCA 粒子间势 + WCA 壁面势） |
| `config/lattice512.cfg` | 512×512×64 粒子阵列（`ratio_rho=2`） |

## 输入输出

C++ 只读写二进制 `.fpd`；初态、检查点、续跑文件是同一格式，文件里的 `step` 决定起点。

```
magic "FPDCKPT\0" | endian 0x01020304 | version=4 | flags (uint32) | Nx Ny Nz N (int32) | step (int64)
vx vy vz [size]  |  Rx Ry Rz Rux Ruy Ruz Vx Vy Vz Fx Fy Fz [N]   （f64，IDX 顺序连续）
[p p_prev [size]]                                                 （仅 flags bit0）
checksum (FNV-1a，覆盖此前全部字节)
```

- 写端写 v4，读端接受 v3（无 flags）与 v4；v2 及更早版本明确拒绝。
- 压力只在 `pressure_solver = split` 且历史完整时入档（flags bit0），用于逐位重启；PCG 路径冷启动
  `p0 = 0`，不存压力。粘度、密度、`S_i`、力密度、RNG 状态和物理参数都不入档。
- 文件大小恒为 v3 `48 + 24*size + 96*N`、v4 `52 + 24*size + 96*N [+16*size]` 字节，
  `fpd_tool --dump-ckpt` 校验这一点。
- 格式由 `src/include/IOBin.h` 与 `tools/fpd_format.py` 共同定义；numpy 形状为 `(Nz, Ny, Nx)`。
- 载入后检查上壁 `vz[:,:,Nz-1]` 逐位为 0（否则泊松相容性失效），不满足即拒绝。
- 检查点路径 `<out_dir>/<run_name>_<step 补零 7 位>.fpd`；终态在 `n_steps` 写出，文件里的 F 与 R 对应。

```bash
# 1. 生成初态：编辑 tools/make_init.py 顶部数据块（网格、半径、POINTS、OUT）后运行，仅依赖 stdlib
#    POINTS = None 时按 PHI 随机生成离壁 ≥ WALL_MARGIN 的不重叠构型；也可直接写坐标列表
python3 tools/make_init.py

# 2. 运行与续跑（同一配置、同一 seed；把检查点当 init_file）
./build/fpd config/sed.cfg
./build/fpd config/sed.cfg --set init_file=out/sed_0010000.fpd --set n_steps=400000

# 3. 离线可视化（必须 pvpython；--config 只用来取 dt）
PV=/home/doll/Software/ParaView-6.2.0-RC1-MPI-Linux-Python3.12-x86_64
$PV/bin/pvpython tools/fpd2vtk.py out/sed_[0-9]*.fpd --config config/sed.cfg -o vis/
# ParaView 打开 vis/sed.pvd：流体 .vti（速度插值到格心）+ 粒子 .vtp（V/F/Ru），时间轴完整
```

`fpd2vtk.py` 可重复运行：产物齐全的文件跳过，`.pvd` 增量合并，产物先写 `.tmp` 再原子改名；
`-j` 默认按可用内存自动定并行度，运行前后各打印一行内存估算与实测。
`tools/verify_vtk.py <pvd> <fpd> --config <cfg>` 把 VTK 读回与源文件逐点比对。
沉降数据用 `out/sed_[0-9]*.fpd` 选择，避免带上 `sed_init.fpd`。

`fpd_tool`：

```bash
./build/fpd_tool --dump-ckpt <f> [--at i j k]              # 头部、尺寸校验；--at 打印格点速度
./build/fpd_tool --diff-ckpt <a> <b>                       # 逐位比较
./build/fpd_tool --verify-forces <ckpt> <config>           # 按配置重算 F（势 + 壁面）与文件对照
./build/fpd_tool --sed-stats out/sed_[0-9]*.fpd --config config/sed.cfg [--window A B] [--bins K]
```

`--sed-stats` 给出 `<Vz>(t)`、质心、最小壁面间隙、φ(z) 剖面（解析球与扩散界面两种体积分数），
误差棒来自 `tests/Analysis.cpp` 的分块平均；找不到平台时拒绝背书，要求延长模拟。

### 逐位重启

同一配置、同一 seed 下，从检查点续跑与一口气跑完**逐位相同**（N=1/N=2 不重叠、常系数与
`ratio_rho=2` 的 PCG 与分裂法三条路径均实测）。依据：

- cuRAND Philox 每步显式定位：`base = step·(slotD+slotN)`，`randD` 从 `base`、`randN` 从
  `base+slotD` 生成，`slotD = 3·size`、`slotN = 3·Nx·Ny·(Nz+1)`（均向上取偶）。不需要恢复 RNG 状态。
  （实测只有 Philox 的 offset 有意义，且「生成 n 个数前进 n」的累计模型不成立。）
- 密度场、粘度场由 R 确定性重算；PCG 固定 `p0=0`；分裂法的压力历史在 v4 检查点里。

粒子支撑域重叠（间距 < 2·range）时，网格 atomic 累加顺序不确定，跨进程会有 ~1e-16 差异，
逐位重启不再成立；有粒子间力的致密体系通常处于这种情况。

## 数值判据（`fpd_check`）

```bash
export ACC_DEVICE_TYPE=nvidia
./build/fpd_check --check              # 力守恒、亚格点、N=2 重叠、J5 力链路（128×64×32）
./build/fpd_check --check-config       # 配置解析、必填、覆盖顺序、归档往返（CPU）
./build/fpd_check --check-potential    # J1–J6 势函数（CPU）
./build/fpd_check --check-tridiag      # Thomas vs 稠密 Gauss、奇异列定规（CPU）
./build/fpd_check --check-poisson      # P5/P6 常系数往返、V0–V5 变系数投影、V6 分裂法
./build/fpd_check --check-wall         # W1–W4、W3b、W4s、W6、W8a–e、W5b（约数分钟）
./build/fpd_check --check-mass         # MA1–MA7 质量场与重力
./build/fpd_check --debug-pcg [Nx Ny Nz]   # PCG 到双精度地板的迭代数（调试，不改生产容差）
./build/fpd_check --pcg-noise [Nx Ny Nz]   # PCG 误差 vs 同一步热噪声增量（PN1），含分裂法一行
./build/fpd_check --check-density-fdt [tol ...] [split]   # W5c 变密度能量均分 + 容差成对比较（约 4 分钟/档）
./build/fpd_check --lambda [L] [ratio_rho] # 相场离散常量表（CPU）
```

当前实测（日志 `baseline/mass_20261017/*.log`、`baseline/pcg_precision_20261019.log`、
`baseline/pcg_tol_20261019/`）：

| 判据 | 内容 | 实测 |
|---|---|---|
| `--check` | 力守恒 `Σf = ΣF`（三个亚格点） | ≤ 6.2e-15 |
| | Σφ² 三方向散布（交错混叠）/ (Σφ)²/Σφ² 的亚格点散布 | 9.185e-5 / 0.0348% |
| | N=2 重叠：力守恒 / 均匀流 `V_i = 1` | ≤ 5.3e-15 / 0 |
| | J5：GPU vs CPU 力 / GPU 两次逐位 / 投影守恒 | 0 / 逐位 / 3.3e-14 |
| P5/P6 | 常系数壁面算子往返（偶/奇 Nz、16×8×8）/ 不相容右端诊断 | ≤ 5.6e-16 / 报 3.811 |
| V0 | 均匀 a=0.5 vs 1.0、η=f=0：`v*` 与 a 无关 | 逐位 0 |
| V1–V3 | `L·1=0` / 对称性 / 变系数往返（16³，三种 a 场） | 逐位 / ≤ 1.3e-15 / ≤ 1.2e-11 |
| V4/V5 | 一步投影后 div v 相对最大速度 / PCG 迭代与残差 | ≤ 8.1e-12 / 14–18 次、≤ 3.2e-13 |
| V6 | 分裂法一步后 div v / a≡1 时与 PCG 之差 | 1.3e-15 / 3.5e-16 |
| W1/W2 | 壁面 `div v`（含两壁层）/ 上壁 `vz` | 3.4e-16 / 逐位 0 |
| W3 | 平面 Poiseuille 与**离散**闭式 `vx[k]=(g/2η)[(Nz²+1)/4-(k-(Nz-1)/2)²]` | 0 |
| W3b | 两片壁面 `Π_yz/Π_zx` 与独立 CPU 展开式 | 4.4e-17 |
| W4/W4s | 贴壁截断下力守恒、均匀流 V；上下壁镜像 `S^z/S^x` | ≤ 3.1e-15、0；0 |
| W6 | z 向动量收支恒等式逐步成立 | 1.3e-15 |
| W8a–e | 壁面力 `=-dU/dz`、镜像反对称、体相零力、平衡点、GPU/CPU 装配 | 2.7e-10、0、0、7.3e-14、0 |
| W5b | 4³ 壁面盒能量均分：γ=2 vs γ≡1 | 1.0079±0.0030 vs 0.9427±0.0028（15.8σ） |
| MA1–MA7 | 面心 a vs 暴力枚举、望远镜求和、重力总和、N=0/等密度恒等、壁面力平衡 | ≤ 3.3e-16、≤ 3.9e-15、≤ 7.4e-15、逐位、6.3e-14 |
| J1–J6 | 力=-dU/dr、Python 黄金表、最小镜像、N=3 三角形、特征点 | ≤ 1.1e-11、≤ 8.2e-16、通过 |
| P0 | `--debug-pcg`：a≡1 时 PCG 1 步到 <1e-15；粒子场 15 步到 1e-12、19 步到 1e-15 | 地板 ~4e-16 |
| PN1 | `--pcg-noise`：默认容差下求解误差 / 同一步热噪声增量（ρ 加权能量 ε_E、粒子速度 ε_V） | 7e-13 / 4e-12 |
| W5c | 4³ 壁面盒、ρ∈[1,2]：`<Σρv²>/(kT·dim)`；pcg_tol 1e-8…1e-3 与 1e-12 同 seed 成对差 | 1.0047±0.0035；≤ 4e-9 |

判据设计要点：

- 力守恒对 `stencil_point` 内部公式错误是**盲的**（投影与归一化一起错）。独立对照是连续球坐标
  积分与暴力枚举（`--lambda`、MA1）；N=2 均匀流 `V=1` 是速度平均误用全场 φ 的判决性判据。
- W3 必须用离散闭式：连续抛物线带 `O(1/Nz²)` 假偏差。
- W5b 是差分判据（γ=2 必须显著优于 γ≡1），绝对比值受 4³ 小盒与显式 Euler 偏差影响。
- 靠近壁面时 `--lambda` 的各向同性判据会 FAIL，这是截断的预期结果；测试粒子须离两壁 ≥ range。

相场离散常量（`--lambda 32`，a=3.2、ξ=1、range=8、模板盒 17³，FACE_X）：

```
∫φ = 170.306086   ∫φ² = 100.841136   三方向散布 2.3e-7 / 9.2e-5
离散 vs 连续球坐标积分：∫φ 差 4.7e-7
```

`∫φ` 比解析球体 `(4/3)πa³ = 137.26` 大 24.1%，这是扩散界面的曲率项
（`+8πaξ²π²/24`），不是格点误差。凡涉及粒子体积的量（重力换算、体积分数）必须用与模拟相同的
离散求和。`--lambda` 另外打印的 λ^T、M_i、M_eq、M_eff 只是 ∫φ、∫φ² 的组合，代码不使用。

## 压力求解精度：误差要比热噪声小多少

投影后的速度误差有精确表达式 `δv = -dt·a·∇e`（e 为压力误差），其 ρ 加权能量
`Σρ|δv|² = dt²·eᵀAe`；下一步右端会把残留的 `div v` 重新投影掉，所以误差**不累积**。
合适的对照尺度是**同一步注入的噪声速度增量**（而非平衡态热速度，后者大 5 倍以上，会把容差放得过松）。
判据门槛取 1e-3：`ε_E = sqrt(Σρ|δv|²/Σρ|PΔv*|²)` 与 `ε_V = max|δV_i|/rms|ΔV_i|` 都 ≤ 1e-3。

32³、N=4、`ratio_rho=2`、kT=0.25、dt=0.002（`--pcg-noise`；右端 `b` 几乎 100% 是噪声）：

| pcg_tol | 迭代 | ε_E | ε_V | 结论 |
|---|---|---|---|---|
| 1e-3 | 4 | 9.0e-4 | 7.9e-3 | 不满足（粒子速度） |
| 1e-4 | 6 | 2.7e-5 | 2.5e-4 | 满足 |
| 1e-6 | 8 | 8.0e-7 | 7.6e-6 | 满足 |
| 1e-12 | 16 | 7.2e-13 | 4.0e-12 | 满足 |
| 分裂法 | 0 | 9.3e-2 | 1.8 | 不满足 |

经验关系：`ε_E ≈ 1.2×` 相对残差，`ε_V ≈ 10×` 相对残差（粒子速度的噪声增量被相场平均压低，
分母小）。分裂法的误差项 `(a-a₀)∇(p̂-p)` 在 ρ=ρ₀ 的流体里为零，只在粒子内部出现。
端到端（`sed.cfg`，5000 步，同 seed，以 1e-14 为参考；不同 seed 的差为 max|ΔR| = 0.41）：

| 方案 | 粒子位置最大偏差 ΔR | 步/秒 |
|---|---|---|
| pcg_tol 1e-10 / 1e-8 / 1e-6 | 3.6e-13 / 3.6e-11 / 3.1e-8 | 38.6 / 44.7 / 53.3 |
| pcg_tol 1e-4 / 1e-3 | 1.8e-7 / 7.5e-6 | 61.1 / 71.5 |
| 分裂法 | 8.3e-3（`<Vz>` 偏 0.55%） | 99.2 |
| 默认 1e-12 | — | 35.9 |

W5c 显示各容差下温度与参考的差 ≤ 4e-9；分裂法的温度偏 +1.9e-4（统计显著但小），但其单步误差在粒子
速度上与噪声增量同量级，不满足 1e-3 门槛。**推荐：有噪声的生产算例用 `pcg_tol = 1e-6`**
（满足门槛的最松档 1e-4 再收紧 100 倍，吞吐 +48%）；确定性算例（`noise_on=0`）的对照尺度不再是噪声，
需要按所测量的精度要求保留紧容差（E2 单粒子在 1e-4 时 ΔR = 2.7e-11）。默认值仍为 1e-12。

## 性能

RTX 3060，`ratio_rho ≠ 1` 时每步一次 PCG（约 16 次迭代，每次一对 2D FFT + Thomas）：

| 算例 | ratio_rho | 步/秒 | 来源 |
|---|---|---|---|
| 128×64×32，N=95（`sed.cfg`） | 2.0 | **35.9**（PCG 16 次/步；20 万步 ≈ 1.5 h） | 2026-10-19 |
| 同上，`pcg_tol=1e-6` | 2.0 | **53.3**（8 次/步；20 万步 ≈ 1 h） | 2026-10-19 |
| 同上 | 1.0 | **134.8**（20 万步 ≈ 25 min） | 2026-10-19 |
| 128×64×32，N=1 | 1.0 | 153 | Phase 2 |
| 32³，N=1 | 1.0 / 2.0 | 1128 / ~100–200 | Phase 2 / 2026-10-17 |
| 512×512×128，N=28611 | 2.0 | ~0.1（前 3 步） | 2026-10-17 |

小盒子已接近 kernel 启动延迟主导（32³ 格点少 8 倍只快 7.4 倍）。PCG 迭代数已与网格无关，
进一步提速的空间在容差（见上节）与非压力部分（粒子模板盒散写、O(N²) 力）。

## 能说什么、不能说什么

| 能说 | 不能说 |
|---|---|
| 同一盒子内的**相对**比较：`U(φ)/U(φ₀)`、终端速度比、浓度剖面演化 | 绝对 Stokes 阻力 / 迁移率 |
| 净重力由密度差给出：`U(2.0)/U(1.5) = 1.9985`（期望 2，E2） | 复现了某个真实胶体的壁面势 |
| 壁面势托住了沉积层（最小间隙 > 0） | 粒子速度涨落满足 FDT（W7 未做） |

`a/L_z = 0.1` 时壁面受限与周期镜像是 10% 量级偏差；单粒子绝对终端速度比解析 Stokes 小约 2.8 倍
（扩散界面的有效半径问题）。绝对值需要 L 外推（Phase 4）。

## 已知局限

- **粒子侧 FDT（W7）未做**：周期专属的 3A/3B 能量均分判据已随 z 周期路径删除，壁面盒下的目标值需
  重新推导；现在只有流体侧 W5b。
- **变密度近似**：`∇·v=0` 加上由 R 预设的 ρ 场，`∂_tρ + v·∇ρ` 的残差为 O(ξ/a)，与 FPD 界面近似同阶。
- **PCG 成本**：变密度路径默认容差下比常系数慢约 3.8 倍（`pcg_tol=1e-6` 时约 2.5 倍）。
- **分裂法不保证 FDT**：它不是 ρ 加权的正交投影；W5c 的温度偏差 2e-4 很小，但粒子速度的单步误差
  与噪声同量级，只适合无噪声或对粒子动力学精度要求低的场合。
- **重叠支撑域**：`η`、`a` 的 atomic 累加使逐位重启只对不重叠构型成立；重叠区 `Σφ > 1`，
  η 可超过 η_c（N=2 间距 4 时 η_max = 90.5），文献按字面求和、未 clamp。
- `stencil_point` 的引用出参在 `-Minfo` 下报 local memory，目前不是瓶颈。

## 历史缺陷与教训

详细过程见 PROGRESS.md。最值得记住的几条：

- **C1/C2（交错归一化塌缩）**：旧代码把三套错开的相场塌缩成一个体心场，力投影与归一化位置不一致；
  速度平均用全场 Σφ 除以单粒子归一化，**N=1 时恰好抵消，骗过 20 万步模拟**。根因是同一段遍历被
  复制三份后漂移，解法是 `Stencil.h` 去重 + N=2 重叠判据。
- **模板盒少一格**：盒宽 `2·range` 在负方向静默漏一层，力守恒抓不到（分子分母一起漏）。由上下壁
  镜像判据 W4s 抓出；当时与之「全位一致」的 Python 参考复制了同一盒尺寸——共享假设的参考不独立。
- **势函数 `s6` 写成三次方**：力差 200 倍，被四阶数值微分 + 独立黄金表抓出。
- **噪声标定**：旧代码 `W=10, dt=0.001` 等于 kT=0.05（噪声几乎关着）；周期盒 3A 测出偏差为
  显式 Euler 的 O(dt) 项（`bias% ≈ 335·dt`，外推为零），3B 测出该偏差不被 η_c=50 放大，
  这是生产用 `dt=0.002` 的依据。两项判据已删除，结论无常驻判据守护。
- **对流项误乘 1/ρ**：变密度初版把整个 `∇·Π` 乘 a，均匀密度下对流加速度错一个因子；由 V0 守住。
- 旧代码中不可照抄的 bug：cell list 邻居表漏赋值导致每个 cell 都把 0 号 cell 当邻居；z 向反射 hack
  在投影之后改 `vz`；重力硬编码 `-10`；读写 header 字段顺序不一致。
