# FPD_NEXT

胶体悬浮液的**流体粒子动力学**（Fluid Particle Dynamics, FPD）GPU 模拟器。

把刚性胶体粒子表示为高粘度的流体区域（相场 φ），流体力学相互作用因此自然涌现，
无需显式求解移动边界。实现基于 OpenACC + cuFFT + cuRAND，单 GPU。

## 参考文献

- **[PRL]** Tanaka & Araki, *Simulation Method of Colloidal Suspensions with Hydrodynamic
  Interactions*, PRL **85**, 1338 (2000) —— `doc/` 下有 PDF 和 md
- **[SM]** Furukawa, Tateno & Tanaka, *Physical foundation of the fluid particle dynamics
  method for colloid dynamics simulation*, Soft Matter **14**, 3738 (2018)

本实现以 **[SM] 2018** 为准（[PRL] 是简化版：无噪声、无转动、力密度未归一化、丢对流项）。
`doc/old_code/` 是上一代 OpenMP 实现，**仅供参考，不是真值**（已知错误见下）。

## 控制方程

```
相场      φ_i(r) = ½[tanh((a - |r-r_i|)/ξ) + 1]
粘度      η(r)   = η_ℓ + (η_c - η_ℓ) Σ_i φ_i(r)
动量      ρ(∂_t + v·∇)v = f - ∇·[Π - σ],   Π = pI - η(∇v + ∇v†)
不可压    ∇·v = 0
力密度    f(r)   = Σ_i φ_i(r) F_i / ∫φ_i dr
粒子速度  V_i    = ∫ v φ_i dr / ∫ φ_i dr
```

## 构建

```bash
cmake -B build && cmake --build build -j
./build/main
```

环境：nvc++ 26.3（NVIDIA HPC SDK）、CUDA 13.1、CMake 4.2.3、RTX 3060 12 GB。

> 注：`build/` 里原有的构建产物来自另一台机器（缓存路径为 `/home/rebecca/VSCode_File/...`），
> 已清理并重新配置。`build/result/` 保留了上一代码的输出，仅供参考。

## 当前状态

| 模块 | 状态 |
|---|---|
| 流体求解器（MAC 交错网格 + FFT 压力泊松） | ✅ 完成且已验证自洽 |
| 相场 / 粘度场构造 | ⚠️ 有缺陷，见 C1-C3 |
| 力投影 / 粒子速度 | ⚠️ 有缺陷，见 C1-C2、C4-C5 |
| 粒子间相互作用力 | ❌ 未实现（`F[n]` 恒为 0） |
| 多粒子 / 初始构型读取 | ❌ 未实现（硬编码 N=1） |
| 热噪声按 kT 标定 | ❌ 未标定（见 C6） |
| 配置文件 / checkpoint | ❌ 未实现 |
| 定量验证 | ❌ 未做 |

**流体求解器为什么可信**：`src/Stokes.cpp` 的泊松解用的是 7 点差分格式的**精确离散本征值**
`0.5/(cos kx + cos ky + cos kz - 3)`，而不是连续谱的 `-k²`，与第 3 阶段的有限差分离散严格自洽。
`cufftPlan3d(Nz,Ny,Nx)` 与 `IDX(i,j,k)=i+j*Nx+k*Nx*Ny` 的搭配也已核对无误。

## 已知缺陷清单

以下每条都已逐行对照源码确认。**修复前不要相信任何物理结果。**

### 致命：交错网格归一化被塌缩

旧代码维护**三套错开的相场** `phiX/phiY/phiZ`（`doc/old_code/main.cpp:326-328`）各配
`sum_phix/y/z`，因此 `Σ_grid fx == Fx[n]` 恒等成立。新代码把三套塌缩成了单个 cell-center 的
`phi_grid`，且存的是 `RATIO_ETA*order()`（50φ 而非 φ）：

- **C1** `Force.cpp:61-63` vs `Viscosity.cpp:79-89` —— 分子用 x-face 错开的 `order(dx+0.5,dy,dz)`，
  分母却是 cell-center 的 `50*order(dx,dy,dz)`。双重错误：50 倍因子 + 交错位置不匹配。
  即使去掉 50，`Σf ≠ F_n` 仍成立，且误差随粒子亚格点位置起伏 → **位置相关的伪力**
- **C2** `Velocity.cpp:39-48` —— `phi_grid` 是**所有粒子 φ 的求和场**，分母却是粒子 n 自己的
  `sum_phi[n]`。**N=1 时 50 和 φ 恰好抵消，所以跑 20 万步也看不出来**；多粒子重叠时直接错

### 其余

- **C3** `Viscosity.cpp:79` —— 用 `RATIO_ETA*order()` = 50，应为 `(RATIO_ETA-1)*order()` = 49
  （文献 `η = η_ℓ + (η_c-η_ℓ)φ`）。当前粒子内 η = 51 而非 50。旧代码此处是对的
- **C4** `Velocity.cpp:42` —— `Vx[n] +=` 前未清零。*实际数值影响极小*：
  `V_new = (V_old + Σvφ)/sum_phi` 而 `sum_phi ≈ 6850`，残留被压掉 4 个数量级。仍须修，
  但它**不是**任何可见漂移的成因，别误判
- **C5** `Velocity.cpp:74` —— PBC 用 `fmod()`，负数得负值，周期边界破了
- **C6** `Stokes.cpp:53` —— `W=10, dt=0.001` 隐含 `kT = W²dt/2 = 0.05`，文献用 **kT=14.3**。
  **当前等于噪声关着**，"20 万步无 NaN"没有验证过噪声通路
- **C7** `Viscosity.cpp:88-89`、`Velocity.cpp:41-48` —— `collapse(4)` 下 atomic 累加到
  `sum_phi[n]`，N=1 时 4096 线程排队打同一地址，GPU 退化成串行
- **C8** `Viscosity.cpp:47`、`Force.cpp:33`、`Velocity.cpp:65` —— `present()` 遗漏 `sum_phi`
- **C9** `Stokes.cpp:176-187` —— 用旧索引约定的死代码
- **C10** `Common.h` —— `dt` 与 `inv_dt` 靠人工维护互为倒数
- **C11** 全部 —— cuFFT/cuRAND 返回码无一检查

### 旧代码里不要照抄的 bug

- `doc/old_code/src/cellList.cpp:141-142` —— `nnIndex[12]` 未赋值、`nnIndex[13]` 赋两次。
  `std::vector<int>` 零初始化 → 每个 cell 都把 **0 号 cell** 当邻居，产生跨盒伪相互作用
- 同上：14 半邻居方案要求每维 ≥ 3 个 cell。`CELL_SIZE=32` 时 `init128x128x64` 得 4×4×2，
  **z 向只有 2** → 自相邻重复计数
- `doc/old_code/main.cpp:798-806` —— z 向反射 hack 在投影**之后**改 `vz`，破坏 ∇·v=0，
  且与全周期 FFT 求解器不自洽
- `doc/old_code/main.cpp:475` —— 重力硬编码 `-10`，命令行 `argv[5]` 读了不用
- `savePositions3D_FPD` 与 `readPositions3D_FPD` 的 header 字段顺序不一致

## 设计约定

- **命名**：小写 = 流体场（`vx`、`fx`），大写 = 粒子量（`Vx`、`Fx`、`Rx`）。贯穿全代码
- **索引**：`IDX(i,j,k) = i + j*Nx + k*Nx*Ny`，**x 变化最快**
- **交错网格**：`v` 在面心，`p`/`eta`/对角应力在体心，`etaXY/YZ/ZX` 与非对角应力在棱边
- **单位**：格距 `dx ≡ 1`，密度 `ρ ≡ 1`，溶剂粘度 `η_ℓ ≡ 1`
- **稳定性**：显式粘性项要求 `dt ≲ ρdx²/(2d·η_c) = 1/(6×50) ≈ 0.0033`

### 交错网格的单一真值源

`src/include/Stencil.h` 的 `stencil_point<Loc L>()` 是**唯一**实现半格偏移约定、球形截断判据、
周期回绕的地方。`Viscosity`/`Force`/`Velocity` 三个模块共用它 —— C1/C2 的成因正是这段循环
被复制粘贴了三份然后各自漂移。

> ⚠️ **代价**：`stencil_point` 成了单点故障，它错了三个模块一起错。
> 而"力守恒 `Σfx == Fx[n]`"这个判据对它的内部公式错误是**盲的**（投影和归一化会一起错）。
> 所以必须同时有独立的数值对照，见 `spike/`。

## 已完成的验证

**spike（`spike/spike_template_routine.cpp`）** —— nvc++ 26.3 下：

- `#pragma acc routine seq` 作用于**模板函数**可行，`-Minfo=accel` 确认
  `stencil_point<(Loc)1>` 生成了 device routine。不需要退化成宏
- 力守恒 `|Σfx - 1.0| = 3.3e-16`（机器精度），粒子在非整数亚格点 (64.3, 32.7, 16.5) 上仍精确
- `∫φ`（FACE_X，a=3.2, ξ=1, 截断半径 8）= **170.3051225223**，与独立写的 Python 实现全位一致，
  与连续球坐标积分（170.306）相对差 5e-4

> 📌 **`∫φ` 比解析球体 `(4/3)πa³ = 137.26` 大 24.1%。**
> 所以 λ^T、M_i 等常量**必须用与模拟相同的离散求和**，用解析公式会带来 24% 误差，
> 直接毁掉能量均分和 Stokes-Einstein 的验证。

已知的次要问题：`-Minfo` 报 `Local memory used` —— `stencil_point` 的引用出参
（`int& ijk, double& w`）溢出到 local memory。当前不是瓶颈，若后续 profiling 显示有影响，
改成返回小结构体。

截断半径的差异：本实现用 `range2 = range² = 64`（半径 8）；旧代码开 `_8BLOCKLOOP_` 时用
`RANGE_CORE² = 49`（半径 7），二者 `∫φ` 差约 0.3%。不追求与旧结果逐位一致。

## 基线

`baseline/phase0_particle_traj.txt` —— 修改前的原始代码跑 19 万步的单粒子轨迹
（`W=10` 即 kT=0.05，dt=0.001）。Phase 1 修完缺陷后数值结果**必然改变**，此文件仅供定性对照。

## 遗留文件

`storage/test_Stokes_only.cpp` —— **已失效，不参与构建**。三处与当前代码不符：
调用了已重命名的 `save_vtk()`；`cufftPlan3d(Nx,Ny,Nz)` 维度顺序与 `main.cpp` 的 `(Nz,Ny,Nx)` 相反；
第 31 行用的是旧索引约定 `(i*Ny+j)*Nz+k`。

保留是因为它有一个有价值的骨架：**无粒子的纯流体驱动**（`eta≡1` + 初始速度块）——
正是 Phase 3 测试 3A（纯流体噪声谱 `⟨|v|²⟩ = 2kT`）需要的结构。改造时以它为起点，
但上述三处必须先修。
