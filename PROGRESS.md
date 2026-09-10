# 进度记录

最后更新：2026-09-07（Phase 7-A 完成）

## 一句话状态

流体求解器、交错网格封装、噪声标定（流体侧+粒子侧）、配置系统与断点重启均已验证；
**粒子间相互作用力（WCA/Morse/LJ126 + 外场）已实现并通过数值判据**；
自检拆成独立可执行 `fpd_check`，`FpdState` 消除了 4 份分配/映射样板；
**z 向壁面泊松求解器（xy 2D FFT + z 向 Thomas）已实现并通过算子往返判据**，
数学推导全文见 `doc/PressurePoisson.md`（尚未接线到 `Stokes.cpp`）。

## 总体规划

8 个阶段。顺序原则：**先修缺陷建立可信基线 → 钉死噪声 → 定量验证 → 才上生产算例**。

| 阶段 | 内容 | 状态 |
|---|---|---|
| 0 | 基线固化、git、README | ✅ 完成 |
| 1 | 交错网格封装 + 缺陷修复 | ✅ 完成 |
| 2 | 配置文件 + 二进制 IO/重启 + CMake | ✅ 完成（FpdState 已在 Phase 5 一并做掉） |
| 3 | λ^T + 噪声标定 + dt 扫描 | ✅ 3A + 3B + 常量表完成 |
| 4 | Stokes 阻力 + VAF/MSD + L 外推 | ⬜ 未开始（下一步） |
| 5 | 粒子间力 + 外场 + fpd_check 拆分 | ✅ 完成（本次） |
| 6 | 生产算例 + 旧代码对照 | ⬜ 未开始 |
| 7 | z 向无滑移壁面（重力沉降） | 🔶 7-A 求解器完成；7-B 接线未开始 |

完整计划（含各阶段的详细设计决策 D0-D8）：`~/.claude/plans/doc-doc-old-code-fuzzy-blum.md`

---

## 已完成

### Phase 0 — 基线固化

- `git init`，仓库 6.9 MB（`doc/old_code/init/` 的 1.6 GB 已排除）
- `README.md` 记录 11 条已核实缺陷（C1-C11）+ 旧代码 5 处不可照抄的 bug
- `baseline/phase0_particle_traj.txt` 归档修改前的 19 万步单粒子轨迹
- 确认工具链：nvc++ 26.3 + CUDA 13.1 + CMake 4.2.3 + RTX 3060 12 GB
  （原 `build/` 是另一台机器的产物，已重新配置）

### spike — 前置风险验证

`spike/spike_template_routine.cpp` 在动核心代码前验证了三件事：

1. `#pragma acc routine seq` 作用于**模板函数**在 nvc++ 26.3 下可行，
   `-Minfo=accel` 确认 `stencil_point<(Loc)1..3>` 生成了 device routine。**不需要退化成宏**
2. 力守恒 `|Σfx − 1.0| = 3.3e-16`（机器精度），非整数亚格点位置上仍精确
3. `∫φ = 170.3051225223` 与独立写的 Python 实现**全位一致**，
   与连续球坐标积分（170.306）相对差 5e-4

> ⚠️ **Phase 7-B 更正**：上面那个 `∫φ` 与它的「5e-4 相对差」都是**模板盒缺陷**
> 的产物（见「已知隐患」第 8 条）。修掉后 `∫φ = 170.3060863408`，与连续球坐标
> 积分差 **4.7e-7**。当时「与 Python 全位一致」之所以成立，是因为那份 Python 参考
> **复制了同一个盒尺寸** —— 共享假设的参考实现不是独立参考。

> 📌 副产品：`∫φ` 比解析球体 `(4/3)πa³ = 137.26` 大 **24.1%**。
> **归因后来被更正**：这是扩散界面的曲率项（`+8πaξ²π²/24 = 33.07`），
> 不是格点离散化误差 —— 格点其实精确到 5e-6。结论（必须用离散求和）不变，
> 但理由是**自洽性**而非格点不准。详见「测试 3B」一节。

### Phase 1 — 交错网格封装 + 缺陷修复

**新增** `src/include/Stencil.h`：`Loc` / `LocOffset<L>` / `stencil_point<L>()`，
把半格偏移约定、球形截断判据、周期回绕收进唯一一份实现，三个模块共用。

**新增** `src/include/Check.h`：`CUFFT_CHECK` / `CURAND_CHECK`。

**修复**：C1（归一化塌缩）、C2（速度平均误用全场 φ）、C3（粘度系数差 1）、
C4（V 未清零）、C5（fmod 破坏 PBC）、C7（atomic 归约热点）、C8（present 遗漏）、
C9（死代码）、C11（返回码未检查）。

**接口变化**：
- `NS_Config` 不再手工维护 `inv_dt`；新增 `make_ns_config()` 作为**唯一**计算
  `inv_dt = 1/dt` 和 `W = sqrt(2kT/dt)` 的地方
- 新增 POD `PhiParams`，`RADIUS`/`XI`/`range` 从编译期常量改为运行时成员
- `order()` 参数化为 `order(dx,dy,dz,radius,inv_xi)`
- `sum_phi` 拆成 `sum_phix/y/z` 三套；`phi_grid` 这个场**已删除**
- `update_particale_VandR` 拆成 `update_particle_velocity` + `update_particle_position`
- 新增不折叠位置 `Rux/Ruy/Ruz`（MSD 用，与 `R` 同 kernel 积分）

**验收结果**（`./build/fpd_check --check`）：

```
力守恒 |err| ~ 1e-15 ... 2e-14        F=(1,2,3)，三个亚格点位置
sum_phi 各向同性离散度 < 0.001%
N=2 重叠：力守恒 8.4e-15，均匀流场下 V 偏差 = 0.000e+00
eta_max 分离 49.85 ≈ 50（C3 修复前是 ~50.85）；重叠 90.5（FPD 已知行为）
20 万步无 NaN，轨迹与 Phase 0 基线高度吻合
```

> N=2 均匀流场检验是 C2 的**判决性判据**：`v ≡ 1` 时 `V_i = ∫vφ_i/∫φ_i = 1`
> 必须精确成立。旧实现用全场 `Σφ_n` 除以单粒子归一化，重叠时 `V > 1`。

### 测试 3A — 纯流体噪声谱（提前做，只依赖 Phase 1）

32³ 立方盒，无粒子（η ≡ 1），检验 `⟨|v|²⟩ = 2kT`。
**完全绕开**相场、力投影、速度平均的全部机制，只检验噪声通路。

dt 扫描（各档总物理时间相同 T=600）：

| dt | 偏差 | 相邻比值 |
|---|---|---|
| 0.02 | +6.6134% | — |
| 0.01 | +3.1844% | 2.08 |
| 0.005 | +1.5602% | 2.04 |
| 0.0025 | +0.7726% | 2.02 |

线性拟合 `bias(%) = 334.6 × dt − 0.104`，残差 < 0.06%，**dt→0 截距 −0.10%（统计误差内为零）**。

> ✅ **噪声实现正确，不需要任何修正因子。**
> `W = sqrt(2kT/dt)` 标定正确；缺失的 `-2/3 δδ` 迹项确实被不可压投影消掉；FFT 投影正确。
> 观测偏差**纯粹是显式 Euler 的 O(dt) 离散误差**，随 dt 线性消失。

日志：`baseline/phase3a_dt_scan.log`

### 常量表 + 测试 3B — 粒子能量均分

**新增** `src/Analysis.{h,cpp}`：`compute_fpd_constants`（含 4³ 亚格点扫描）、
`equipartition_target`、`blocking_analysis` / `pick_plateau`（分块平均）。
**新增** `main.cpp` 的 `--lambda` 和 `--equipart <ghost|frozen|moving>` 两个模式。

#### 常量表（`--lambda`，纯 CPU 秒级）

```
∫φ = 170.3061  ∫φ² = 100.841   λ^T = 1.688855   M_i = 287.6223   M_eff = 431.4335
Σφ_α² 三方向散布 9.185e-5     4³ 亚格点扫描散布 0.0348%
```

三种**互相独立**的方法交叉验证：离散求和 287.6407 / 连续球坐标积分 287.6409 ——
**差 7e-7** / 扩散界面解析展开。（Phase 7-B 修掉模板盒缺陷前，离散与连续差 0.008%，
见「已知隐患」第 8 条。）
**0.0348% 是整个方法的系统误差下限** —— 比早先用 4 个手挑位置估的 0.009% 大，
因为 4³ 扫描覆盖了整个格胞，这才是诚实的数。

#### 目标值与三处推导更正

```
⟨|V|²⟩ = 2kT/M_i − 2kT/(ρL³)
```

**更正 1（实质错误）**：k=0 修正系数是 **−2** 不是 −3。基线式 `2kT/M_i` 来自 Parseval
对**全部** k 求和、每个 k 权重 `tr P = 2` —— k=0 已按权重 2 计入。代码把它冻结为 0
（`Σ_r v_α` 逐步严格守恒于 0），故减 2 个自由度。L=32 时是 −0.878%，不是 −1.317%。

**更正 2**：「φ_x 是 φ 的平移，故 `|φ̂_x|²=|φ̂|²`」的论证**错了** ——
那是**重新采样**不是平移，混叠分支带不同相位。但实测 `Σφ_α²` 散布 9.19e-5，
结论保住，理由换掉。已把这条判据并入 `--check`。

**更正 3**：`∫φ` 超出解析球体 24.1% 是**扩散界面的曲率项**，不是格点误差（见上）。

#### 结果（32³，kT=0.25，T=2000）

| 模式 | η_c | 位置 | dt | seed | 比值 | σ | 分块平台 |
|---|---|---|---|---|---|---|---|
| ghost | 1 | 冻结 | 0.01 | 1 | 0.9952 ± 0.0349 | 0.1 | b=64（勉强） |
| frozen | 50 | 冻结 | 0.002 | 1 | **1.0050 ± 0.0175** | 0.3 | b=128 ✅ |
| frozen | 50 | 冻结 | 0.002 | 2 | 1.0072 ± 0.0150 | 0.5 | ❌ 无平台 |
| moving | 50 | 更新 | 0.002 | 1 | 0.9712 ± 0.0149 | 1.9 | b=16（过早） |

> ✅ **核心结论：O(dt) 偏差【不会】被 η_c=50 放大。**
> frozen 两个种子给 +0.50% / +0.72%，与「流体式」预期 **+0.67%** 高度吻合；
> 「按扩散数放大 50 倍」情景预期 **+33.5%**，被排除约 19σ。
>
> 物理原因：粒子速度是 φ 支撑域上的**低通平均**，高 k 被滤掉，
> 而偏差 `~ν k² dt/2` 对高 k 最大。低通抑制（约 16×）与界面粘度放大（约 25×）大致抵消。
>
> **实际影响：生产可继续用 dt = 0.002，无需降到 0.0005，Phase 4 的 L=192 仍然可行。**

日志：`baseline/phase3b_equipart.log`、`baseline/phase3b_constants.txt`

#### ⚠️ 3B 遗留的未解决问题

1. **`moving` 比 `frozen` 低约 3.4%（约 2σ）**。它的分块平台在 b=16 就"找到"了
   （frozen 是 b=128），很可能是噪声造成的**假平台、误差棒被低估**。
   物理上 `R += dt·V` 是否严格保测度（`∇_R·V = 0`）只在离散模板恰好匹配时成立，
   所以**可能是真实的小系统效应，也可能只是涨落**。需 2~4 倍长运行 + 第二个种子定论
2. **frozen seed 2 未找到分块平台**，误差棒不可信（中心值与 seed 1 一致是好迹象）。
   分块守卫按设计拒绝背书 —— **这是功能不是缺陷**
3. **本次无法验证亚 1% 的偏差**：`T = (4/3)τ₂/ε²` 给出 3σ 分辨 0.5% 需约 7500 万步

### Phase 2 — 配置 + 二进制 IO + 断点重启

**新增**：`src/Config.{h,cpp}`（key=value 解析，解析/覆盖/dump 共用同一张字段表）、
`src/IOBin.{h,cpp}`（`.fpd` 格式）、`src/tool_main.cpp`（`fpd_tool`）、
`tools/fpd_format.py`、`tools/make_init.py`、`tools/fpd2vtk.py`、`tools/verify_vtk.py`。
**删除**：`src/IO.cpp` + `IO.h`（ASCII VTK 整个移除）。
**CMake**：拆 `fpd_core` 静态库 + `fpd` + `fpd_tool`。

**核心简化**：初始构型直接用 `.fpd` 检查点格式，于是 C++ 只有一条输入路径、
除配置文件外没有任何文本解析；「init 文件」与「重启文件」是同一个 key 同一段代码。

**验收结果**：

| 判据 | 结果 |
|---|---|
| 配置：未知 key / 必填缺失 / 类型错 / init 与 config 维度冲突 | 全部非零退出并指出错误项（未知 key 还给拼写建议） |
| C++/Python 轴序互操作 | `--at 3 1 1` 报 `3001001` ✓ |
| **逐位重启 N=1** | 5 步 vs 2+3 步，全部数组**逐位相同** ✓ |
| **逐位重启 N=2** | 同样**逐位相同**（粒子不重叠，无 atomic 竞争） |
| 检查点自校验 | 改一字节 → 校验和不匹配；截断 → 数据区截断；均非零退出 |
| VTK 往返（pvpython 读回比对） | 时间轴 6 步、双 block、所有数组相对差 **0.000e+00** |
| 插值方向与轴序自检 | PASS |
| Phase 1 / Phase 3 回归 | 全绿，`--lambda` 与 baseline 完全一致 |
| 热路径开销 | **3.9 ms/检查点**（8.39 MB）；旧 ASCII 约 100 ms，快 26 倍 |

#### RNG：实测推翻了原方案

计划原打算「累计 `rng_draws` + `curandSetGeneratorOffset`」。实测发现两件事：

1. **只有 Philox 支持有意义的 offset**。XORWOW（原来的 `CURAND_RNG_PSEUDO_DEFAULT`）
   和 MRG32K3A 试了系数 1/2/4 都对不上，**根本无法逐位重启**
2. **「生成 n 个数后序列前进 n」这个模型是错的**：`offset=1*n` 能对上第 2 批，
   但 `2n`、`3n` 都对不上连续生成的对应批次。所以**累计记账从前提上就不成立**

改用 **slot 方案**：每次生成前显式定位到 `offset(step,c) = (2*step+c)*size*3`。
已验证不同 slot 之间零共同值、互相关 ~1/√n（统计独立），且同一 offset 逐位可复现。
于是重启无需恢复任何 RNG 内部状态，只需 seed 相同 —— 逐位相等是构造上保证的。

证据：`spike/spike_rng_offset.cpp`、`spike_rng_slice.cpp`、`spike_rng_advance.cpp`。

#### 顺手拆掉的三个地雷

1. `Fx/Fy/Fz` 输出前漏 `update host`（现在无害只因力恒为 0，Phase 5 后会写陈旧零值）
2. NaN 早退的 `return 2` 跳过所有 `exit data` 和 cuFFT/cuRAND 清理 → 改为先落盘再正常清理
3. 输出目录不存在时静默丢数据 → `ensure_dir` + 全路径错误检查

#### 新踩的坑

**`#pragma acc update` 不能放在 lambda 里** —— 闭包捕获让运行时查不到设备映射，
报 `data in update host clause was not found on device`。已改成普通结构体 + 内联 pragma。

---

### Phase 5 — 粒子间力 + 外场 + fpd_check 拆分

本次按用户要求**提前做了 Phase 5 的力**（原排期在 Phase 4 之后），理由：拆分自检是
Phase 4/5 的共同前置（Phase 4 要加 VAF/MSD 分析路径，不能再往 943 行的 main.cpp 塞），
`FpdState` 正是长跑路径要用的样板去重，两个阶段互不阻塞。

**三个前置 spike**（见 `baseline/spike_state_and_pair.log`）：

- `spike_state_map`：运行时 API（acc_copyin/acc_create/...）建的映射能被另一个 TU 的
  `present()` 认账；**成员函数里的 pragma 也实测可行**（推翻「成员函数= lambda 会挂」
  的担忧，CLAUDE.md 只记录过 lambda）。FpdState 最终用运行时 API，`acc_is_present`
  直接支撑 `require()` 运行时断言。
- `spike_pair_bench`：**推翻了「粒子间力在 CPU 侧算」的决策记录**。实测 CPU 串行在
  N=363 时 2.7 ms（占 42%！），GPU 全矩阵 129 µs（占 2%）—— 决策记录漏算了 O(N²)
  计算本身，且「CPU 才能逐位重启」不成立（GPU reduction 逐位可复现）。这是本项目
  **第四次「实测推翻了计划里的假设」**。结论：力在 GPU 算（gang-over-i +
  vector-reduction-over-j 全矩阵，零 atomic）。
- `spike_potential_ref`：独立 Python 实现三种势，产出黄金表。

**交付**：

| 模块 | 文件 | 说明 |
|---|---|---|
| `FpdState` | `src/State.{h,cpp}` | 4 份分配/映射/释放样板收成 1 份，`require()` 把 present 遗漏从「跑出垃圾数」变成「启动即 abort」 |
| 势函数 | `src/include/Potential.h` + `src/Potential.cpp` | WCA/Morse/LJ126 + none；WCA 不是独立分支（≡ LJ + 派生 rcut + 移位）；pair_energy 与 pair_force 两份独立实现 |
| 外场 | `Config.h` 的 `gravity_x/y/z` + `gravity_compensate` | 背景力密度 bg=−ΣF/size 抵消 k=0 漂移 |
| 自检拆分 | `src/CheckMain/CheckStencil/CheckNoise/CheckPotential.cpp` | 独立可执行 `fpd_check`，`fpd` 旧命令打印迁移提示 |
| RNG 统一 | — | 3A/3B 从 XORWOW 改 Philox（与生产一致，CMakeLists 注释早已要求） |
| 判据 J1-J7 | — | 见下 |

**验证结果**：

```
--check-potential（纯 CPU，无卡可跑）：全绿
  J1 力=-dU/dr 四阶差分    max 相对差 ~1e-12   （抓到真实 bug：s6 算成三次方，力差 200 倍）
  J2 Python 黄金表对照     max 相对差 ~1e-16
  J3 最小镜像 vs 暴力枚举    J4 N=3 等边三角形   J6 特征点 + 零判据
--check（含 J5）：全绿
  GPU vs CPU 对照 0 相对差   GPU 两次逐位  力守恒 Σf==ΣF 及 bg 补偿 Σf==0（~1e-14）
fpd_tool --verify-forces：正向 PASS（F 自洽），负向 FAIL（配错参数检测到）
smoke.cfg 逐位回归：potential=none 默认逐位相同
N=4 不重叠逐位：跨进程逐位相同
吞吐：N=100 加力 111.2 vs 无势 110.9 步/秒（力代价 <0.3%）
3B frozen seed=1（Philox）：0.9747 ± 0.0337（0.7σ PASS），与 baseline 1.005 统计一致
```

**关键发现：FPD 逐位重启只在粒子支撑域不重叠时成立**。实测 N=8 间距 10.5（< 2·range=16.8，
重叠）跨进程有 ~1e-16 差异（v 场），N=4 间距 17.2（不重叠）逐位相同。根源是
`Viscosity.cpp` 的 `eta[ijk] += d_eta*w` atomic 累加顺序跨进程不确定 —— 这是 FPD 的
固有特性，**不是粒子间力引入的 bug**（potential=none 的 N=8 重叠构型同样不逐位）。
而「粒子间力」需要粒子靠近（间距 < rcut ≤ 15），靠近必然重叠，所以「有力 + 逐位重启」
在当前盒子下**原理上不可兼得**。J8 因此改为「GPU 力确定性（J5 的 GPU vs CPU 对照 +
GPU 两次逐位）+ 不重叠构型的完整逐位重启」，而非「N=8 重叠粒子的逐位重启」。

---

### Phase 7-A — 压力泊松求解器（xy 2D FFT + z 向三对角）

Phase 7 要做的 z 向无滑移壁面，挡在最前面的是 `Stokes.cpp` 的压力泊松求解器整体是
3D 全周期 FFT，z 不再周期后它整个失效。本轮按用户要求**只做求解器 + 判据 + 文档**，
用合成右端的算子往返判据验证，`Stokes.cpp` 一行不动（生产路径零风险）。数学推导全文
见 **`doc/PressurePoisson.md`**（与 `Poisson.cpp` 共同定义同一个 `div∘grad` 算子）。

**交付**：

| 模块 | 文件 | 说明 |
|---|---|---|
| 求解器 | `src/Poisson.{h,cpp}` | `solve_pressure` 按 `cfg.wall_z` 分派：周期 3D FFT（算术与 Stokes.cpp 一致）/ 壁面 xy 2D 批量 FFT + z 向 Thomas；`build_tridiag_coeffs` 预算前推系数 |
| 配置 | `Common.h` 的 `NS_Config` 加 `wall_z` | 默认 0，`make_ns_config` 加默认参数，现有 4 处调用点零改动 |
| 接线 | `State.{h,cpp}` | `tri_w` 数组 + `plan_xy` 句柄（`cufftPlanMany` batch=Nz） |
| 判据 | `src/CheckPoisson.cpp` | `--check-tridiag`（纯 CPU）/ `--check-poisson`（需 GPU） |
| 文档 | `doc/PressurePoisson.md` | 从离散 MAC 网格起的完整推导 + 实测数值 |

**关键设计**（详见文档）：Neumann 边界**不是假设**，而是从「壁面法向面不修正」导出；
`(0,0)` 奇异列用 `d_0 -= 1` 定规（有证明，非罚函数、非近似），**不能整列清零**
（那会删掉支撑粒子重量的静压）；归一化因子是 `Nx·Ny` 不是 `size`（2D 变换点数）。

**验证结果**（实测，非推导值）：

```
spike_fft2d_batch   P1 轴序 1.0e-14  P2 归一化 == Nx*Ny  P3 inembed 逐位相同
spike_tridiag       T1 Thomas vs Gauss 2.8e-16  T2 定规 |x0|=2.2e-16  T3 不相容漂移 496
--check-tridiag     三对角全列对照 2.8e-16；奇异列定规 |x0|=8.2e-15
--check-poisson     壁面往返 4.3e-16（8x4x6）、2.7e-15（128x64x32）；周期往返 2.0e-15
                    相容性诊断对不相容右端报 |Σb|=3.811（非 0）
```

> 算子往返抓到过一个真实 bug：`thomas_sweep` 里 z 向步长写成 1 而非 `Nx·Ny`，
> 往返误差 O(1)。这是「手写独立算子 + 求解器」对照的判决性体现。

日志：`baseline/phase7a_poisson.log`

---

## 待办：下一步

**建议直接进 Phase 4（Stokes 阻力 + VAF/MSD）。** Phase 0/1/2/3/5/7-A 全部完成，
没有阻塞项：3B 解除了 dt 的不确定性（生产可用 0.002），Phase 2 给了配置系统和断点重启，
Phase 5 给了粒子间力和 fpd_check 拆分，Phase 7-A 把壁面泊松求解器提前落地（未接线，
不占生产路径）。Phase 7-B（壁面接线）是独立轨道，与 Phase 4 互不阻塞，可并行。

Phase 4 的关键设计已经想清楚：**必须做 L∈{64,96,128,192} 外推** ——
a/L=0.05 时 Hasimoto 有限尺寸修正是 **14% 偏差**，远大于要验证的几个百分点，
不外推就得不出可信的 λ^T。而从阻力反解出的 λ^T 与 `--lambda` 的 1.68885 对照，
是**静态平衡 vs 稳态耗散两条完全不同的物理路径**交叉验证同一个常量，证据力很强。

### 其余待办

- **Phase 4**（Stokes 阻力 + VAF/MSD + L 外推）：下一步，无阻塞项。
  `FpdState` 已就位，加 VAF/MSD 分析路径时直接复用。
- **Phase 7-B**（z 向无滑移壁面接线）：7-A 已把求解器（`Poisson.cpp`）做出来并通过
  算子往返判据。7-B 是把壁面 BC 接进生产：`Stokes.cpp` 重构为显式 `v*`（现在的源项是
  `div(v*)/dt` 的展开式，边界上要扣的 `v*z[−1]` 分散在三处、打不了补丁）、`Wall.h`
  作 z 向邻居访问的唯一真值源、相场在壁面截断（`FACE_Z` 有效范围 `k ∈ [0, Nz−2]`）、
  粒子侧壁面排斥势与 `min_image` 分方向、配置项 `boundary_z`、端到端 W1–W8 判据。
  **3A/3B 的能量均分目标值要重新推导**（模式基不再是纯 Fourier）。旧代码在这里是
  错的（z 向反射 hack 破坏 ∇·v=0），无照抄对象。壁面上线后 `gravity_compensate` 默认转 0。
- **cell list**：363 粒子 O(N²) 足够。GPU 方案下迁移阈值比原来的 `N > 2×10⁴` 还高
  （S0-C 实测 N=1000 仍只占 11%），具体以 N 标度实测为准。
- **Phase 6**：生产算例，边界条件三选一（先复现旧行为建立对照，别一次改两个变量）
- **性能优化未立项**：128×64×32 只有 153 步/秒。Phase 4 的 L=192 是 27 倍网格量，
  若排期吃紧，瓶颈在 `Stokes.cpp` 散度计算的约 18 次邻居访问（`k±1` 跨步 64 KB）

### 备选池（有价值但未立项）

- **离线精确预测器** `C = kT(I − Ã·dt/2)⁻¹`（逐本征值精确）。
  已验证该框架能重现 3A 四个数据点到 2.1~2.5%（精确式 6.46/3.11/1.53/0.757%
  vs 实测 6.61/3.18/1.56/0.773%，残差近乎恒定比例，很可能是对流项）。
  它能把「跑几小时测偏差再外推」变成「秒算精确目标 + 一次确认」，
  且是验证亚 1% 偏差的**唯一**可行路径
- **`adv_on` 开关**（`Stokes.cpp` 6 处对流表达式各乘一个系数，6 行）。
  关掉对流后系统严格线性高斯，比值**必须**是 1.000 —— 能造出来的最强单元测试。
  **若将来出现无法解释的偏差，这是第一个该加的东西**
- **粘性项改 Crank–Nicolson**：平稳协方差在**任意 dt** 下精确等于 `kT·I`，
  且所需噪声幅度 `BBᵀ = 2kT·A·dt` **恰好就是代码现在用的**，一个字不用改。
  代价是每步内嵌变系数 CG。若将来 dt 成为瓶颈，这是原理性的解

---

## 关键决策记录

| 决策 | 选择 | 理由 |
|---|---|---|
| 交错网格管理 | **自建轻量封装**，不用 AMReX/PETSc DMStag | 核对 C1-C11 后发现库的类型系统只能挡住 C5（`fmod`）。C1/C2 根因是 `sum_phi[n]` 是 per-particle **标量**，`Field<FACE_X> / double` 在任何库里都合法。真正成因是重复代码漂移，解药是去重不是换框架。AMReX 的 GPU 模型是 nvcc+lambda，与 OpenACC 不是一条路，采用等于重写已验证正确的求解器 |
| 封装形式 | `stencil_point<L>()` 返回 `(ijk, w)`，**不传 lambda** | OpenACC 的 `#pragma acc loop` 藏不进模板，nvc++ 对 compute region 里的 lambda 支持有边界。pragma 结构留在调用点，不赌编译器行为 |
| `Field<L>` 类型标签 | **降级为可选**，暂未实现 | 挡不住 C1（标量除法合法）。若与 `present()` 配合不顺就放弃，不为类型漂亮牺牲 GPU 映射可控性 |
| 转动自由度 | 暂不实现 | 用户决定。保持 [PRL 2000] / 旧代码的层次 |
| 势函数 | enum + POD 参数 + `routine seq` 的 switch | 用户要可切换。同一份源码 CPU/GPU 都能编译，无函数指针（OpenACC 设备端不可靠） |
| 粒子间力 | **GPU 侧算**，gang-over-i + vector-reduction-over-j 全矩阵（零 atomic） | `spike_pair_bench` 实测**推翻**了原「CPU 侧算」决策：CPU 串行 N=363 时 2.7 ms（占 42%！），GPU 全矩阵 129 µs（占 2%）。原论据漏算了 O(N²) 计算本身，且「CPU 才能逐位重启」不成立（GPU reduction 逐位可复现，与 sum_phi/V 同机制） |
| 参数输入 | 配置文件 + `.fpd` 初始构型（由 Python 生成） | 用户决定弃用旧 init 文本格式。结果：C++ 只剩一条输入路径、零文本解析，旧格式的「三重循环顺序写反不报错」陷阱直接消失 |
| 可视化 | C++ 只写二进制，Python 离线转 VTK | 用户决定。热路径开销从 ~100 ms/次降到 3.9 ms；且 Python 侧用 VTK 自带 writer，不手写 XML，schema 出错的风险类别整个消失 |
| RNG offset | slot 方案（按 step 定位），**不做累计记账** | 实测推翻原方案：XORWOW 无法逐位重启；且「生成 n 前进 n」的模型是错的 |
| 构建系统 | CMake（用户中途装上了） | 曾因本机无 cmake 临时写过 Makefile，已删除，避免两套构建系统漂移 |
| 长程稳定性判据 | 跑 20 万步而非计划的 200 万步 | 噪声在 kT=0.05 下等于关着，信息量低。真正考验是 Phase 3 在 kT=14.3 下的 dt 扫描。GPU 机时花在测试 3A 上更值 |
| 3B 范围 | **最小版本**：只回答「偏差是否被 η 放大」 | 用户决定。两种情景差 50 倍，10% 精度足够判别。代价是无法验证亚 1% 偏差（已记入遗留问题） |
| 3B 的 kT | 0.25 而非 1.0 | 信噪比与 kT 无关（信号和涨落同比例）、τ 与 kT 无关，所以降 kT **不花任何代价**，却把对流非线性风险降 4 倍。3A 在 kT=1 下相对精确式有 ~2.4% 残差，很可能就是对流项 |
| 3B 的 ghost 模式 | 靠 `ratio_eta=1.0` 实现，**不写特例分支** | C3 修复留下的 `d_eta = ratio_eta − 1` 让它零成本，且走**完全相同**的代码路径（`sum_phi` 照算、`stencil_point` 照跑），只是 η 场是平的 |
| 采样间隔 | 在**物理时间**上固定（`samp_every = 0.5/dt`） | 固定步数会让 dt 越小样本越相关，同样物理时间得到更多但更相关的样本，分块平台反而更难出现 —— pilot 实测踩到了这个坑 |
| 自检代码去处 | **独立可执行 `fpd_check`** | 用户决定。main.cpp 涨到 943 行成瓶颈，自检拆出后 `fpd` 只剩生产路径；`fpd` 旧命令打印迁移提示而非静默 |
| 势参数进 `.fpd` header | **不进** | Python 侧用不到（现有 5 个物理量都是 make_init/可视化要用的）；加进去制造「文件说 Morse 配置说 WCA」的静默冲突。自洽性由 `fpd_tool --verify-forces` 兜底 |
| WCA 实现 | **不是独立分支**（≡ LJ + 派生 rcut=2^{1/6}σ + 移位） | U'(rc)=0 使能量移位与力移位恒等，switch 只剩 3 分支，消除旧代码 cul_WCA/cul_LJ6 两份复制粘贴漂移 |
| LJ/Morse 截断移位 | 默认 `shift=energy`（U 移位，力不变） | `force` 会改阱深（Morse α=1 rcut=15 时 50→49.57），是物理改变。残余力比 `|U'(rc)|/max|U'| > 1e-6` 报警，把「照抄 22.2」变成有数字支撑的决策 |
| RNG 发生器 | 3A/3B 统一 Philox | CMakeLists 注释白纸黑字「验证与生产走同一代码路径」，且 XORWOW 无法逐位重启。代价是 3A/3B 数字变（统计等价），重跑 3B frozen 确认比值落在误差棒内 |
| 外场 k=0 漂移 | `gravity_compensate`（默认 1）+ 背景力密度 bg=−ΣF/size | 均匀外场 ΣF=N·g≠0 使流体 k=0 线性加速、整盒漂移，废掉 3B 目标公式。壁面 Phase 7 上线后此开关默认转 0（壁面提供真实动量汇） |
| `FpdState` 映射机制 | **运行时 API**（acc_copyin/acc_create/...） | spike_state_map 实测：运行时 API 建的映射能被另一 TU 的 present() 认账，成员函数 pragma 也可行；运行时 API 可查询（acc_is_present 支撑 require()） |
| 壁面泊松 z 向求解 | **xy 2D 批量 FFT + 逐 (kx,ky) Thomas 三对角**，不用 DCT-II | cuFFT 无 DCT；2N 延拓/twiddle 技巧是「看起来对、归一化悄悄错」的典型。Thomas 约 20 行、无归一化歧义，且可扩展到非均匀 z 网格或 Robin 壁面。前推系数与右端无关 → 预算一次 |
| 壁面泊松奇异列 (0,0) | **`d_0 -= 1` 定规**，而非整列清零 | 有证明（`𝟙ᵀA'x = −x₀` → `x₀=0` 且精确满足原方程），非罚函数、非近似。整列清零会删掉支撑粒子重量的静压（该列的 z 结构是水平均匀的竖直压力分布） |

---

## 已知隐患

0. **FPD 逐位重启只在粒子支撑域不重叠时成立**（间距 > 2·range = 16.8）。重叠时
   `eta[ijk] += d_eta*w` 的 atomic 累加顺序跨进程不确定（~1e-16 差异）。这是 FPD 的
   固有特性，不是粒子间力引入的 bug（potential=none 的 N=8 重叠同样不逐位）。而
   「粒子间力」需要粒子靠近（间距 < rcut ≤ 15），靠近必然重叠 —— 「有力 + 逐位重启」
   当前盒子下原理上不可兼得。逐位重启回归用 N=1/N=2 不重叠构型。
1. **性能**：128×64×32 是 **153 步/秒**、32³ 是 **1128 步/秒**（GPU 均 99% 饱和）。
   注意 32³ 格点少 8 倍却只快 7.4 倍 —— 已接近**启动延迟主导**，
   所以「小盒子会快很多」的直觉不成立。瓶颈在 `Stokes.cpp` 散度计算的约 18 次邻居访问，
   `k±1` 跨步 `Nx*Ny` = 64 KB 对缓存不友好。
   Phase 4 的 L=192 是 27 倍网格量，按实测数排期。**优化不在当前范围**
2. **`moving` 模式比 `frozen` 低约 2σ，未解决**（见 3B 遗留问题）
3. **`stencil_point` 是单点故障** —— 去重的代价是它错了三个模块一起错，
   且力守恒判据对其内部公式错误是盲的。力守恒必须作为常驻回归测试
4. `-Minfo` 报 `Local memory used` —— `stencil_point` 的引用出参溢出到 local memory。
   当前不是瓶颈，若 profiling 显示有影响，改成返回小结构体
5. **截断半径与旧代码不同**：本实现 `range2 = range² = 64`（半径 8），
   旧代码开 `_8BLOCKLOOP_` 时用 `49`（半径 7），`∫φ` 差约 0.3%。不追求逐位一致
6. **C++ 与 Python 共同定义 `.fpd` 格式**（Phase 2 引入的新漂移点）。
   两处独立实现同一个布局和同一个 FNV-1a 哈希，改动任一侧必须重跑互操作判据 ——
   `--dump-ckpt --at 3 1 1` 应报 `3001001`，`fpd2vtk.py --self-test` 应全 PASS。
   轴序（numpy 要 reshape 成 `(Nz,Ny,Nx)`）搞反不报错，只得到一个转置的场
7. **`doc/PressurePoisson.md` 与 `Poisson.cpp` 共同定义 `div∘grad` 算子**（Phase 7-A
   引入的新漂移点）。改任一侧必须同步更新另一侧的 §12 对照表，并重跑
   `--check-poisson`（`CLAUDE.md` 已列为第 7 条不可违反约束）。往返判据对
   「算子与修正步不匹配」是盲的 —— 真正的 `∇·v=0` 端到端检验要等 7-B 接线后的 W1。
8. ✅ **模板盒少一格（Phase 7-B 发现并【已修复】）**。`stencil_point` 的模板盒原本取
   `[kn-range+1, kn+range]`（宽 `2·range`），而半径 `range` 的球心落在格胞任意位置时
   能触及的整数层最多 `2·range+1` 个 ⇒ 负方向静默漏一层。**先于壁面存在**，周期模式下
   远离壁面就能复现（FACE_X @ Rz=16.5：盒遍历 2140 点 vs 暴力枚举 2152 点）。
   后果：`∫φ`/`M_i`/`λ^T` 带方向相关、位置相关的伪偏差（C1/C2 同一失败类别）；
   **力守恒判据抓不到**（分子分母一起漏，比值自洽）。
   修法：`n_range = 2*range+1`、盒 `[kn-range, kn+range]`，已用暴力枚举对全部
   `Loc × 位置` 验证零漏点。代价：模板盒 16³→17³，吞吐 −4%（N=8@32³ 实测）。
   判决者：连续球坐标积分（离散 vs 连续从 0.008% → **7e-7**）与暴力枚举。
   抓它的是 Phase 7-B 新增的 W4「上下壁对称」判据 —— 只测一侧的判据对它完全盲。

---

## 提交历史

```
ce2ef34  Phase 7-A: 壁面压力泊松求解器(xy 2D FFT + z 向 Thomas) + 算子往返判据 + 数学推导文档
3aa0f97  Phase 5: 主循环接线粒子间力 + 背景力密度补偿 + 端到端判据
d7643fb  Phase 5: 配置接线势参数 + 最小镜像硬约束 + 参数查表
e1a747e  Phase 5: 势函数模块 + 纯 CPU 判据 (--check-potential)
9fa2853  重构: 提取 FpdState + 拆分 fpd_check + RNG 统一 Philox + 修 config.used 覆盖 bug
9ed91b5  Phase 2: 配置文件 + 二进制 .fpd IO + 逐位重启, C++ 不再产出可视化格式
1574432  测试 3B: 粒子能量均分验证 + 常量表, 证实 O(dt) 偏差不被 eta 放大
48b8fb3  补充 CLAUDE.md 与 PROGRESS.md
9695701  测试 3A: 纯流体噪声谱验证, dt 扫描证实偏差为 O(dt) 且外推为零
c0c0e53  Phase 1: 交错网格封装 + 修复 C1-C5/C7-C9/C11
ae13f10  Phase 0: 基线固化 - README 记录缺陷清单, spike 验证模板 routine seq 可行
6b2a9c7  初始提交：FPD_NEXT 现状基线
```

## 文件地图

```
src/            Stencil.h(交错网格唯一真值源) Common.h(POD+工厂) Check.h
                Stokes/Viscosity/Force/Velocity(物理)  Analysis(常量表+分块平均)
                Config(key=value) IOBin(.fpd) State(分配/映射/释放唯一持有者)
                Potential(势函数+最小镜像+力装配)  Poisson(压力泊松: 周期 3D FFT / 壁面 2D FFT+Thomas)
                Tests.h(自检声明)
                main.cpp(生产)  CheckMain/CheckStencil/CheckNoise/CheckPotential/CheckPoisson.cpp(fpd_check)
                tool_main.cpp(fpd_tool)
tools/          fpd_format.py(格式镜像) make_init.py(纯stdlib) fpd2vtk.py verify_vtk.py(pvpython)
spike/          template_routine(模板+routine seq) rng_offset/rng_slice/rng_advance(RNG 语义)
                state_map(映射机制) pair_bench(力算 GPU 判决) potential_ref(黄金表)
                fft2d_batch(2D 批量 FFT 轴序/归一化) tridiag(Thomas vs 稠密 Gauss + 定规)
config/         smoke.cfg production.cfg
baseline/       各阶段实验日志与参考轨迹
doc/            PressurePoisson.md(泊松数学推导，与 Poisson.cpp 共同定义算子)  old_code/(仅供参考)
```
