# 进度记录

最后更新：2026-10-20（blog 修改清单；待用户处理事项汇总）

## 当前状态

生产路径完整可用：MAC 交错网格流体求解器、x/y 周期 + z 无滑移壁面、相场粘度与力投影、
粒子间势与粒子-壁面势、变密度（`ratio_rho`）惯性与重力体力、`.fpd` v3 初态/检查点/逐位重启、
离线 VTK 转换。变密度算例的 `pcg_tol` 由用户在跑 `fpd` 前用独立程序 `fpd_calib` 标定
（判据：残差导致的粒子速度误差 ≤ 同一步热噪声增量 × 1e-3）；默认仍为 1e-12。

尚未完成的物理验证：粒子侧 FDT（W7）、绝对 Stokes 阻力（Phase 4 的 L 外推）、
沉降运行到稳态给出 `U(φ)`、旧代码对照。

## 阶段

| 阶段 | 内容 | 状态 |
|---|---|---|
| 0 | 基线固化、git、缺陷清单 | ✅ |
| 1 | 交错网格封装 `Stencil.h` + 缺陷 C1–C11 修复 | ✅ |
| 2 | 配置文件 + `.fpd` 二进制 IO + 逐位重启 + CMake | ✅ |
| 3 | 相场离散常量表 + 噪声标定（3A）+ 粒子能量均分（3B） | ✅ 结论保留；3A/3B 判据于 8-A 删除 |
| 4 | Stokes 阻力 + VAF/MSD + L 外推 | ⬜ 未开始 |
| 5 | 粒子间势 + `FpdState` + `fpd_check` 拆分 | ✅ |
| 6 | 粒子-壁面势 + 集体沉降算例 + 旧代码对照 | 🟡 壁面势 ✅、沉降能跑通；稳态 `U(φ)` 与旧代码对照未做 |
| 7 | z 向无滑移壁面（7-A 求解器，7-B 接线） | ✅ |
| 8-A | 删除三维 FFT 与 z 周期路径 | ✅ |
| 密度 | `ratio_rho` 变密度惯性 + 重力体力 + 变系数投影（PCG） | ✅ 2026-10-17/18 |
| 压力精度 | PCG 容差判据、真残差与每步收敛检查、`fpd_calib` 标定程序；分裂法评估后移出 | ✅ 2026-10-19 |

## 判据覆盖

当前判据与实测值见 README「数值判据」。按模块：

| 模块 | 判据 |
|---|---|
| 相场投影/平均（Stencil、Force、Velocity、Viscosity） | `--check`：力守恒、亚格点、Σφ² 混叠、N=2 重叠均匀流；W4/W4s 壁面截断 |
| 势与力装配 | `--check-potential` J1–J6；`--check` 的 J5（GPU/CPU、逐位、投影守恒） |
| 常系数压力 | `--check-tridiag`；`--check-poisson` P5/P6 |
| 变系数压力与变密度显式步 | `--check-poisson` V0–V5；`--debug-pcg`（P0）；`--pcg-noise`（PN1，fpd_calib 同一核心）；`--check-density-fdt`（W5c：变密度能量均分与容差成对差） |
| 壁面流体力学 | `--check-wall` W1、W2、W3、W3b、W6；噪声 √2 由 W5b 统计差分判据负责 |
| 壁面势 | W8a–W8e |
| 质量场与重力 | `--check-mass` MA1–MA7 |
| 配置 | `--check-config`（34 项） |
| IO | `--dump-ckpt` 尺寸守卫；`--at 3 1 1` = 3001001；`fpd2vtk.py --self-test`；`verify_vtk.py` 往返 |

不在常驻覆盖内（历史结论或缺口）：3A 的 `bias% ≈ 335·dt`、3B 的「O(dt) 偏差不被 η_c 放大」、
W5'a（FDT 矩阵恒等式，随 `adv_on/noise_skip` 删除）、粒子速度涨落（W7）。

## 已知缺口与隐患

1. **粒子侧 FDT 是空洞**。3B 删除后没有任何粒子速度涨落统计；W4 的 `v≡1 ⇒ V=1` 是运动学，
   不是热力学。补 W7 需要：用自由度计数（`doc/PressurePoisson.md` §14.4 的模板）重新推导壁面盒
   目标值，不能沿用周期盒的 `-2kT/(ρL³)`；逐分量 `⟨Vx²⟩/⟨Vy²⟩/⟨Vz²⟩`；只跑冻结粒子（贴壁时
   `S^z < S^x`，且移动粒子会随机游走出界）；粒子离两壁 ≥ range ⇒ `Nz ≳ 24`。变密度下惯性来自
   ρ(r) 场，目标值要按变密度流体推导。补 W7 时顺带恢复「O(dt) 偏差不被 η_c 放大」的守护。
2. **dt 偏差无常驻判据**。`dt=0.002` 的依据是 3A/3B 的历史测量（`baseline/phase3a_dt_scan.log`、
   `phase3b_equipart.log`）。可行替代：W5b 跑 `dt` 与 `dt/2`，断言 `|ratio-1|` 约减半。
3. **W5b 只是差分判据**。4³ 壁面盒 γ=2 为 1.0079±0.0030（2.6σ 偏离 1），γ≡1 为 0.9427±0.0028。
   绝对一致需更大盒子 + dt 外推 + 更长采样。（早期记录的 γ=2 = 0.9750、8.0σ 是 Stokes 三 kernel
   重构前的旧基线。）
4. **逐位重启只对不重叠支撑域成立**（间距 > 2·range）。重叠时 `eta`、`inv_rho_*` 的 atomic 累加
   顺序跨进程不确定（~1e-16）。有粒子间力的致密体系通常重叠，长跑不能靠重启获得逐位可复现。
5. **`stencil_point` 是单点故障**：Viscosity/Mass/Force/Velocity 共用，力守恒判据对其内部公式
   错误是盲的。独立对照：连续球坐标积分（`--lambda`）、暴力枚举（MA1、spike）。
6. **两侧共同定义的格式/算子**：`.fpd` 由 `IOBin.h` 与 `fpd_format.py` 共同定义（轴序搞反不报错，
   只得到转置场）；常系数/变系数算子与 `doc/PressurePoisson.md` §12/§16 对应。
7. **`wz_edge_idx` 没有越界兜底**：范围是调用方前置条件，`acc routine seq` 里做不了运行时检查。
8. **壁面附近各向异性判据会 FAIL 是设计行为**：`Analysis.cpp` 的 `spread_int_phi2`/`spread_subgrid`
   假设三方向等价，测试粒子须离两壁 ≥ range；检验截断对称性用 W4s。
9. **PCG 成本**：迭代数与网格无关（κ ≤ ρ_max/ρ_min），`sed.cfg` 默认容差 35.9 步/秒、
   `pcg_tol=1e-6` 53.3、`1e-4` 61.1 步/秒、等密度 134.8 步/秒。`ratio_rho` 很大时迭代数按 √κ 增长。
10. **确定性算例的容差**：`fpd_calib` 的门槛以噪声为尺度，`noise_on=0` 时不适用（程序直接说明并退出）；
    这类算例按所测量的精度保留紧容差。
11. **绝对阻力**：单粒子终端速度比解析 Stokes 小约 2.8 倍（扩散界面有效半径）；
    `a/L_z = 0.1` 时壁面与镜像偏差 10% 量级。
12. `-Minfo` 报 `stencil_point` 引用出参用 local memory；目前不是瓶颈。
13. 截断半径 `range = ⌊2(a+ξ)⌋ = 8`，旧代码用 7，`∫φ` 差约 0.3%；不追求与旧结果一致。

## 待用户处理（2026-10-20 汇总）

代码侧的工作都已提交（HEAD `0b91b9b` 及之前），以下事项等用户决定或手动处理：

1. **blog 核对**：`doc/` 下推导笔记与代码的出入已逐条列在根目录 `tmp.md`
   （不进 git，经 `.git/info/exclude` 排除）。重点：ns-equation §3.3 力投影/速度平均的分母写成了
   Σ_i φ_i（即旧 C2 缺陷）、§1 把 `pi_*` 当成 Cauchy 应力、§5 守恒写法 warning 与第二篇矛盾；
   λ^T/M_i/M_eff 相关内容建议删除；「ā 缩放降迭代数」说法错误；discretization §2.2 只令 z 分量散度
   为 0、`p[K]` 笔误；PressurePoisson §16.1 噪声符号前后不一致；W5b、W3b、V0–V5、PCG 停止条件、
   `conv_*`、行号等过时数值与接口。
2. **`pcg_tol`**：默认仍为 1e-12。有噪声的变密度配置跑前先 `fpd_calib <配置> --state <检查点> --warm 0`，
   把推荐值（已测算例均为 1e-5）手动写入配置；是否改默认值待定。
3. **未入库的配置**（`config/512x512x128*.cfg`、`lattice512_2/3.cfg`、`mass_destribution.cfg`）：
   - `512x512x128_2.cfg` 的 `init_file` 是第 600000 步检查点而 `n_steps = 600000`，启动即「无事可做」；
   - `lattice512.cfg`、`lattice512_2.cfg`、`512x512x128.cfg` 的 `gravity_z = +0.58716` 指向 +z（沿用旧配置的正号），需确认是否有意；
   - `gravity_z = -100` 已按新语义换算为 `-0.58716`；首行构型注释已改正。
4. **`doc/` 已移出版本控制**（用户提交 `8a0913b`，`.gitignore` 忽略整个 `doc/`）：笔记、论文、`old_code/`
   与 `doc/pressure_split_failed.md` 只在本机。README/AGENTS 里对 `doc/` 的引用指向本地文件；
   分裂法的代码仍可从 `git show a424dfd`、`0f22da7` 找回。

## 下一步

1. **W7：冻结粒子在壁面盒中的能量均分**（缺口 1），同一统计框架补 dt 守护（缺口 2）。
2. **沉降延长到稳态**，给出同盒相对 `U(φ)/U(φ₀)`；400 时间单位时分块平均找不到平台。
   注意 `sed.cfg` 现为 `ratio_rho=2`，20 万步约 1.5 h。
3. **Phase 4**：L ∈ {64, 96, 128, 192} 外推 Stokes 阻力（a/L=0.05 时 Hasimoto 修正约 14%）。
   壁面上线后可直接用沉降终端速度测（背景力密度恒 0）。
4. **旧代码对照**：`doc/old_code/` 只在少数处可对照；先复现旧行为建立对照，别一次改两个变量。

备选池（有价值但未立项）：

- **大 `ratio_rho` 的预条件**：κ ≳ 10 时考虑两级预条件或多重网格。
- **离线精确预测器** `C = kT(I - Ã·dt/2)⁻¹`：曾重现 3A 四个数据点到 2.1–2.5%，是验证亚 1%
  偏差的唯一可行路径。
- **`adv_on` 开关**：关掉对流后系统线性高斯，能量均分比值必须是 1.000；出现无法解释的偏差时先加它。
- **粘性项 Crank–Nicolson**：任意 dt 下平稳协方差精确为 kT·I，所需噪声幅度与现有一致；代价是每步
  变系数 CG。
- **cell list**：GPU 全矩阵在 N=1000 时力仅占 11%；N=28611 的 512 网格算例需按实测再评估。

## 关键决策记录

| 决策 | 选择 | 理由 |
|---|---|---|
| 交错网格管理 | 自建 `stencil_point<L>()`，不用 AMReX/PETSc | C1/C2 根因是重复代码漂移，库的类型系统挡不住（`Field / double` 在任何库里都合法）；AMReX 的 GPU 模型与 OpenACC 不是一条路 |
| 封装形式 | 返回 `(ijk, w)`，不传 lambda | pragma 结构留在调用点，不赌 nvc++ 对 compute region 内 lambda 的支持 |
| 模板盒 | `2*range+1` | 半径 range 的球触及 `2*range+1` 层；窄盒负方向漏层，代价约 4% 吞吐 |
| 转动自由度 | 不实现 | 用户决定 |
| 势函数 | enum + POD + `routine seq` switch；WCA ≡ LJ + 派生 rcut | 设备端不用函数指针；WCA 下能量与力移位恒等，消除两份复制 |
| 截断移位 | 默认 `energy` | `force` 改变阱深（Morse α=1、rcut=15 时 50→49.57），是物理改变 |
| 粒子间力 | GPU 全矩阵 gang/vector reduction，零 atomic | spike 实测 CPU 串行 N=363 占 42%，GPU 占 2%；reduction 逐位可复现 |
| 输入 | 配置文件 + Python 生成的 `.fpd` | C++ 只有一条输入路径、零文本解析 |
| 可视化 | C++ 只写二进制，Python 离线转 VTK | 热路径开销从 ~100 ms 降到 ms 级；用 VTK writer，不手写 XML |
| `.fpd` 内容 | v3 只存续跑状态（速度 + 粒子数组） | 压力、派生量、物理参数、RNG 状态都可重算或来自配置；势参数入档会制造「文件与配置冲突」 |
| RNG | Philox + 每步显式 offset（slot 方案） | 实测 XORWOW/MRG 无法逐位定位，「生成 n 前进 n」的累计模型错误 |
| `FpdState` 映射 | OpenACC 运行时 API | 映射能被另一 TU 的 `present()` 认账；`acc_is_present` 支撑 `require()` |
| 验证程序 | 独立 `fpd_check`/`fpd_tool`，共用只含生产代码的 `fpd_core` | 验证与生产走同一代码路径；默认构建只编生产目标 |
| z 边界 | 只保留 x/y 周期 + z 无滑移 | 两套语义 = 两倍漂移面；代价是 3A/3B 判据 |
| 壁面泊松 | xy 2D 批量 FFT + 逐列 Thomas，不用 DCT | cuFFT 无 DCT；2N 延拓归一化易错；Thomas 无歧义且可扩展 |
| 奇异列 (0,0) | 相容投影 + `d_0 -= 1` 定规 | 有证明、非罚函数；整列清零会删掉支撑粒子重量的静压 |
| 壁面切向 BC | ghost 代入后的显式壁面 kernel | 壁面上没有 vx/vy 自由度；三 kernel 直接对应离散定义域 |
| 棱边布局 | `Nz+1` 层，逻辑 k 与存储解耦 | 下壁 k=-1 与上壁 k=Nz-1 在 Nz 层布局下撞同一槽位 |
| 壁面剪切噪声 | ×√2（γ=2） | 离散 FDT 要求 `γ·w² = w`，壁面控制体 w=1/2；W5b 用 γ≡1 对照证明被数据选中 |
| 粒子越界 | 不 clamp，落盘中止（rc=3） | clamp 把穿墙变成看似正常的轨迹；检查只在检查点/日志步做，避免每步同步 |
| 壁面势自变量 | 表面间隙 h | 改 radius 时平衡间隙不变 |
| 重力 | 加速度；体力 `(ρ-1)g` 经 `apply_gravity` 进入面心力场 | 浮力由密度差给出；「力 = -dU/dr」与「重力 = 质量 × 加速度」两条路不互相污染 |
| 背景力密度 | `bg` 形参保留，生产恒传 0 | 壁面是真实动量汇；全周期才需要的 `gravity_compensate` 已删除 |
| 变密度显式步 | 1/ρ 只乘力密度，对流项不乘 | ρDv/Dt 除以 ρ 后对流系数为 1；V0 守护 |
| 变系数投影 | PCG，预条件器 `-ā·M` 复用常系数求解器，p0=0 | 不复制 FFT/Thomas；条件数与网格无关（≤ ρ_max/ρ_min），16 次/步已接近最优；warm start 会破坏逐位重启 |
| 压力求解精度判据 | 误差与**同一步噪声增量**比，门槛 1e-3（ε_E 与 ε_V） | 用户提出「误差比热噪声小即可」；对照尺度取单步增量而非平衡态热速度（后者大 5 倍，会放松过头）；误差不累积（下一步重新投影） |
| `pcg_tol` 怎么定 | 独立程序 `fpd_calib`，用户跑 `fpd` 前手动运行，推荐值手动写入配置 | 判据依赖真实配置与状态（噪声占比、ε_V/残差）；不做运行时自适应（逐位重启与判停逻辑更简单）；实测 τ* 与 kT 无关 |
| 压力分裂法 | 移出主路径，记录在 `doc/pressure_split_failed.md` | 吞吐 ×2.8，但噪声主导的压力步间几乎不相关，外推误差在粒子速度上比噪声增量大（ε_V = 1.8） |
| `.fpd` 版本 | 回退 v3 | v4 只为分裂法存压力历史；PCG 路径冷启动 p0=0，不需要压力 |
| 主循环 | 四段抽到 `src/Step.cpp` | `fpd` 与 `fpd_calib` 共用，标定测到的就是生产路径 |
| 配置防御 | 只保留必填、语法、未知键、名称检查 | 个人使用、输入可信；不做范围/组合检查与预测警告 |
| 统计 | 分块平均，找不到平台拒绝背书；采样间隔按物理时间 | 相关样本当独立会低估误差约 7.5 倍；固定步数采样会随 dt 变相关 |
| 判据形态 | W3 用离散闭式；W5b 用差分 | 连续抛物线带 O(1/Nz²) 假偏差；小盒量测机器自带系统偏差 |

## 文件地图

```
src/              生产代码（fpd_core + fpd）
  main.cpp        配置解析、run_production 主循环、越界/NaN 检查、每步 PCG 收敛检查、检查点
  Config.cpp      字段表（解析/--set/dump 共用）、必填与势参数检查、POD 工厂
  Stokes.cpp      噪声、Π 三 kernel、v*、泊松右端、投影、真残差
  Step.cpp        主循环四段（forces / fields / fluid / particles），fpd 与 fpd_calib 共用
  Poisson.cpp     常系数 2D FFT + Thomas；变系数算子 + PCG（可选 PcgTrace）
  Viscosity.cpp   η 场与 S_i^x/y/z；Mass.cpp 面心 1/ρ 与重力体力
  Force.cpp       F → f 面心投影；Velocity.cpp v → V 与位置更新
  Potential.cpp   势参数派生量、GPU/CPU 两份力装配
  State.cpp       FpdState：主机存储、设备映射、FFT plan、Philox
  IOBin.cpp       .fpd v3 读写、FNV-1a
  include/        Common.h(POD/IDX/工厂) Stencil.h Wall.h Potential.h Check.h 及各模块头文件
tests/            fpd_check（Check*.cpp + Analysis.cpp + PcgCalib.cpp）、fpd_tool（tool_main.cpp + Analysis.cpp）、
                  fpd_calib（calib_main.cpp + PcgCalib.cpp）
  Analysis.*      常量表、亚格点扫描、分块平均（fpd_check 与 --sed-stats 共用）
  spike/          一次性风险验证（RNG offset、模板 routine、FFT 批量、三对角、状态映射、势黄金表等）
tools/            fpd_format.py(格式镜像) make_init.py(初态) fpd_config.py(读配置数值)
                  fpd2vtk.py verify_vtk.py(pvpython)
config/           smoke / production / sed / lattice512
baseline/         各阶段日志；当前判据基线 mass_20261017/、pcg_precision_20261019.log、pcg_tol_20261019/、pcg_calib_20261019/
doc/              推导笔记与参考文献（用户维护），pressure_split_failed.md（失败尝试记录），old_code/ 仅供参考
```

---

## 历史记录（倒序）

### 2026-10-20 blog 修改清单

对照 HEAD `0b91b9b` 逐篇核对 `doc/` 下六篇推导笔记（ns-equation、fourier-projection、discretization、
PressurePoisson、PoissonCodeWalkthrough、Stokes_line_by_line），把与代码矛盾、过时或建议补充的点
汇总到根目录 `tmp.md`（全局 9 条 + 每篇一张表，标 [错]/[旧]/[补]，附位置、应改内容与依据）。
`tmp.md` 不进 git；按 AGENTS.md 约定 `doc/` 由用户自行修改。本次无代码改动。

### 2026-10-19 `fpd_calib` 容差标定程序；压力分裂法与 `.fpd` v4 移出主路径

用户决定：容差由独立程序在跑主代码前标定并手动写入配置；分裂法作为失败的尝试移出主路径。

- **移除**：`pressure_solver` 键、`NS_Config.split`、`FpdState.p_prev`、`step_navier_stokes` 的分裂分支、
  `apply_laplacian_affine`、V6、W5c/`--pcg-noise` 的分裂行；`.fpd` 回退 v3（v4 只为分裂法存压力历史）。
  记录写入 `doc/pressure_split_failed.md`（方法、代码、数据、失败原因、恢复方法）。
  PCG 路径 50 步检查点与移除前逐位相同，判据日志与 `mass_20261017` 逐字节相同。
- **`src/Step.cpp`**：主循环四段抽成共用函数，`main.cpp` 与 `fpd_calib` 共用；改动前后检查点逐位相同。
- **`fpd_calib`**（`tests/calib_main.cpp` + `tests/PcgCalib.cpp`）：读真实配置与状态，推进生产轨迹，
  在测量点上对 1e-2…1e-12 每半量级一档按生产判停求解，量 ε_V（残差导致的粒子速度误差 / 同一步
  热噪声导致的粒子速度增量）与 ε_E，推荐 τ*/margin（默认 10）。`--pcg-noise` 改为调用同一核心，
  PN1 与旧实现数值一致。
- **实测**：`sed.cfg` 在 kT = 0.01/0.25/1 下 τ* 都是 1e-4（rms|ΔV_noise| 严格 ∝ √kT，噪声占比不变），
  证实 τ* 与 kT 无关；kT=0.01 且重力 10 倍仍为 1e-4（预热 200 步时沉降流尚未发展，生产标定应用
  运行中的检查点）。`lattice512.cfg`（512×512×64，N=7803）τ* = 1e-4，耗时 3.5 min、内存 6.8 GB。
  推荐值（余量 10）均为 1e-5。

### 2026-10-19 变密度压力求解的精度验收

用户问：变密度后 FFT 只能当预条件器，PCG 是不是最好的；怎么检验精度够；「误差比热噪声小」对不对。

**PCG 是否最优**：预条件后条件数 `κ ≤ ρ_max/ρ_min`（`xᵀAx = Σa_f(∇x)²`），与网格、粒子数无关；
CG 界 `2q^k`、`q = (√κ-1)/(√κ+1)`，`ratio_rho=2` 时到 1e-12 约 16 次，与 32³ 到 512×512×128 实测一致。
多重网格同为 O(1) 迭代且每次不比 FFT 便宜，不值得。`ā` 缩放在精确算术下不改变迭代序列（原注释
「降一个量级」是错的，已改）。真正的杠杆是容差。

**判据**（门槛 1e-3，用户确认）：投影后速度误差 `δv = -dt·a·∇e`，`Σρ|δv|² = dt²eᵀAe`，与同一步
噪声增量（同一 RNG slot 下 W 正常与 W=0 的 v* 之差）投影后的 ρ 加权能量比 → ε_E；粒子平均后 → ε_V。

- S1：`pcg_tol` 配置键；`diag[2]` 真残差 `‖div vⁿ⁺¹‖/‖div v*‖`；主循环每步按主机侧 `diag[1]` 判收敛
  （原先只在日志步看 `> 1e-6`，放松容差后会误报，且中间步失败看不见），未收敛写出 step+1 并 rc=3。
  默认路径与旧二进制 50 步检查点逐位相同。
- S2 `--pcg-noise`（32³、N=4、`ratio_rho=2`、kT=0.25）：`b` 中噪声占 99.8–100%。ε_E ≈ 1.2×残差、
  ε_V ≈ 10×残差；满足门槛的最松容差 1e-4（6 次迭代，ε_V 2.5e-4）；1e-3 不满足（ε_V 7.9e-3）。
- S3 W5c `--check-density-fdt`（4³、ρ∈[1,2]、300k 步）：`<Σρv²>/(kT·dim)` = 1.0047±0.0035——
  变密度路径的噪声/FDT 此前没有任何判据。同 seed 成对差：1e-8/1e-6/1e-4/1e-3 分别 3e-13/6e-12/2e-10/4e-9。
- S4 `sed.cfg` 5000 步同 seed（参考 1e-14，seed 差 max|ΔR| = 0.41）：1e-6 时 ΔR 3.1e-8、1e-3 时 7.5e-6；
  吞吐 1e-12 35.9、1e-6 53.3、1e-4 61.1、1e-3 71.5 步/秒。E2（noise_on=0）1e-4 时 ΔR 2.7e-11。
- S6 压力分裂法（Dodd & Ferrante 2014，`pressure_solver=split`）：V6 一步 div v 1.3e-15，a≡1 时与 PCG
  差 3.5e-16；吞吐 99.2 步/秒。但历史取自参考轨迹时 ε_E = 9.3e-2、ε_V = 1.8；端到端 ΔR 8.3e-3
  （seed 差的 2%）、`<Vz>` 偏 0.55%；W5c 温度偏 +1.9e-4。原因：噪声主导的压力在步间几乎不相关，
  外推 `2pⁿ-pⁿ⁻¹` 误差与压力本身同量级，且误差项 `(a-a₀)∇(p̂-p)` 只在粒子内部非零。无噪声时（E2）
  ΔR 2.4e-8，可用。为逐位重启新增 `.fpd` v4（flags bit0 存 p、p_prev），PCG 路径载荷与 v3 逐字节相同。

**结论**（当时；之后分裂法与 v4 已移出主路径、容差改由 `fpd_calib` 标定，见上一条）：
PCG 保留为默认；有噪声的生产算例推荐 `pcg_tol = 1e-6`（门槛最松档 1e-4 的 100 倍余量），
默认值是否改由用户决定。日志 `baseline/pcg_tol_20261019/`。需用户核对的 doc 章节：
`doc/fourier-projection.md` §6（PCG、ā 缩放的作用、判据）、`doc/PressurePoisson.md` §16。

### 2026-10-19 消除 p 导出；PCG 精度调试；文档同步

- `FpdState::download` 不再回传 p；`fpd_tool --dump-ckpt` 校验文件大小 `48 + 24*size + 96*N`，
  防止整场数组回流。W6 改为显式 `update self(p)`。
- `Poisson.h` 新增可选 `PcgTrace`（生产传 null，tol/上限/算术逐位不变）；`fpd_check --debug-pcg`
  在 a≡1、高斯、球、真实粒子场上测到 1e-8…1e-16 的迭代数。32³ 粒子场 15 步到 1e-12、19 步到 1e-15，
  误差地板约 4e-16；a≡1 时 PCG 1 步，FFT 直解地板约 1e-15。生产容差维持 1e-12。
- 文档同步：README/AGENTS/PROGRESS 按当前代码重写。纠正的过时内容：RNG offset 公式
  （旧文档写 `(2*step+c)*size*3`，代码是 `step*(slotD+slotN)`）、W5b 数值、W3b 已验证、
  `make_init.py` 没有固定默认算例、`fpd --noise` 提示早已不存在、沉降算例吞吐（变密度后约 1.5 h/20 万步）。
  重新构建后九个 `fpd_check` 模式输出与归档逐字节相同；`sed.cfg` 吞吐实测见 README「性能」，
  记录在 `baseline/docsync_20261019/`。
- 删除文档中关于 λ^T / M_i / M_eq / M_eff 的讨论（用户确认）：代码里没有粒子质量，惯性只来自
  ρ(r)，这些量用不到。`--lambda` 仍打印它们，但只当 ∫φ、∫φ² 的组合；`Mass.h`、`Analysis.*` 的
  相应注释同步删改。
- 代码注释同步：`Stokes.h` 的 RNG offset、`State.h` 的数组分组、`Config.cpp` 的 `init_file` 帮助、
  `make_init.py` 的「默认算例」、`CheckMain.cpp` 用法、`CheckWall.cpp` 的 W5b 数值，以及引用已不存在的
  「CLAUDE.md 约束 N」「tmp.txt」的注释。只改注释与帮助文字，数值路径不变。
- 未入库的 512 网格配置：`gravity_z = -100`（旧语义下直接施加在每个粒子上的力）按
  `-100/170.306` 换算为加速度 `-0.58716`；删除过时的「最低一层位于 z=0」注释并改正首行构型说明。

### 2026-10-17/18 密度对比：`ratio_rho`、变密度惯性、重力体力

**动机**：论文 [SM] 式 (3) 写明液体与粒子等密度；原实现把它实现成 ρ ≡ 1，重力却作为每粒子恒力
绕过质量直接加到粒子上。改为：

```
rho(r) = 1 + (ratio_rho - 1) Σφ          a = 1/rho（面心）
v* = v + dt·[ -∇·(vv) + a·(∇·τ + ∇·σ̃ + f + (rho-1)g) ]
∇·(a∇p) = ∇·v*/dt,   v = v* - dt·a·∇p
```

**实现**：新增 `Mass.cpp`（`update_mass_fields` 先建 ρ 再取倒数；`apply_gravity` 加 `(1/a-1)g`）；
`Poisson.cpp` 新增 `apply_laplacian_var` + `solve_pressure_var`；`compute_particle_forces` 删除外场形参，
`ExternalField` 移到 `Common.h` 且语义改为加速度；`ratio_rho` 必填，`density_on = (ratio_rho != 1)`；
`FpdState` 常驻三个面心 a，`density_on` 时另分配三个 PCG 暂存，`diag` 扩到 2。

**修正（10-18）**：初版把整个 `∇·Π` 乘 a，均匀密度下对流加速度会小一个因子（用户审 blog 时给出反例）。
改为把总通量散度拆成对流（`conv_*` 唯一实现，不乘 a）与其余部分（乘 a）；新增 V0：均匀 a=0.5 vs 1.0、
η=f=0 时 `v*` 逐位相同。`density_on=0` 保留合并写法。

**判据**：`--check-mass` MA1–MA7 与 `--check-poisson` V0–V5 全绿（数值见 README）。密度关闭路径与
`ratio_rho=2`、非零重力、N=1 的 5 步 vs 2+3 步重启都逐位相同。

**E2 终端速度比**：单粒子 32³、无势、`noise_on=0`、`g_z=-0.02`，30000 步：
`U(2.0)/U(1.5) = 1.9985`（期望 2，阻力相同 ⇒ 比值等于净重力比）。

**配置迁移**：`smoke/production` 用 `ratio_rho=1`；`sed.cfg` 用 `ratio_rho=2`、
`gravity_z = -10/170.306 = -0.058718`（旧每粒子恒力 -10 的等净力换算；贴壁 z 向 ∫φ 小 0.47%，
平衡点 4.7189 → 4.7191）；512 网格配置按 `100/170.306 = 0.58716` 换算。旧 `.config.used` 缺 `ratio_rho`
会响亮失败；旧检查点里的 F 含重力，`--verify-forces` 对旧文件会 FAIL。

**旧基线不再作为门槛**（用户确认）：保留独立解析判据与新判据，不要求与历史日志逐位一致。

### 2026-10-10 转换器增量化与并行化

`fpd2vtk.py` 为边下载 808 MB/个存档边转换而改：产物齐全即跳过、`.pvd` 增量合并、先写 `.tmp` 再原子
改名（半截文件不会被误判完成）；`imap_unordered` + `spawn`（numpy 的 OpenBLAS 线程池使 fork 不安全）。
首版并行度只按「可用内存 ÷ 单 worker 峰值」算，512×512×128 单 worker 峰值 2.73 GB、开 7 个进程压垮
31 GB 机器。修正：逐分量插值并及时释放源场、`numpy_to_vtk(deep=0)`，峰值降到 1.55 GB；默认只用一半可用
内存并留 4 GiB；运行前后打印内存估算与实测。`MEM_PER_CELL = 50` 与 `write_vti` 写法绑定。
808 MB 真实档案与 32³ 数据改前改后产物逐字节一致，`--self-test` 与 `verify_vtk.py` 通过。

### 2026-09-27 CMake 按目标拆分

`src/CMakeLists.txt` 定义 `fpd_core`/`fpd`，`tests/CMakeLists.txt`（`EXCLUDE_FROM_ALL`）定义
`fpd_check`/`fpd_tool`；默认构建只编生产程序，可执行仍输出到 `build/` 根目录。
同期开始 512×512×64 粒子阵列大算例（`config/lattice512.cfg`）。

### 2026-09-23 初态生成脚本精简

`make_init.py` 删除命令行接口，改为直接编辑顶部数据块（网格、半径、`POINTS`、速度图样、输出路径）。
`POINTS = None` 按 `PHI` 随机生成离壁 ≥ `WALL_MARGIN` 的不重叠构型（128×64×32 时 N=95，即 `sed.cfg`
的初态）；也可直接写坐标列表。脚本同时打印标称与扩散界面两种体积分数。

### 2026-09-22 Stokes 应力三 kernel；泊松代码阅读

应力装配拆成体心/XY 棱、内部 YZ/ZX 棱、两片壁面三个 kernel；壁面段直接写代入 ghost 后的 `∓2ηv`，
删除 `wz_vx/wz_vy/wz_vz` 与带上壁 off-by-one 的 `wz_edge_at_wall/wz_noise_gamma`。用户决定删除只服务
W5'a 的 `adv_on/noise_skip`（W5'a 随之删除，覆盖损失）。新增 W3b 直接核对两片壁面剪切通量。
当时 GPU 驱动不可用，W3b/W5b 待验证；之后已在 `baseline/mass_20261017/check_wall.log` 通过
（W3b 4.4e-17，W5b 1.0079 vs 0.9427）。同日新增 `doc/PoissonCodeWalkthrough.md`，`solve_wall`
包装并入 `solve_pressure`。记录在 `baseline/stokes_3kernel_20260922/`。

### 2026-09-21 `.fpd` v3

文件头只留 magic/endian/version、Nx/Ny/Nz/N、step；数据区只留速度和粒子数组；删除压力、flags、物理量
与 RNG 记录。配置或 `--set` 出现尺寸键即报错；离线工具用 `--config` 取物理量。两种算例 v2/v3 共用载荷
逐位相同，v3 的 50 步与 25+25 步重启逐位相同。记录在 `baseline/io_v3_20260921/`。

### 2026-09-20 配置精简

按个人使用、输入可信处理：保留必填、语法、未知键与势/移位名称检查及实际穿墙/文件/GPU 错误检查；
删除参数范围、组合限制、重复键拒绝、拼写建议、预测警告与 `h_rest/z_rest` 归档。`Config.cpp` 674 → 405 行；
新增 `--check-config`。修正截断上界说明为 `min(Nx,Ny)/2`（z 非周期）。前后判据日志逐字节相同，检查点与
重启逐位相同。记录在 `baseline/simplify_20260920/`。

### 2026-09-14 Phase 6：粒子-壁面排斥势与集体沉降

无滑移速度边界不托粒子，而 z 不 clamp、越界即中止，所以没有壁面势的重力沉降必然以 rc=3 收场。

- 自变量取表面间隙 `h`；势族逐字复用 `pair_force_over_r`，J1 自动覆盖公式，新增的只有几何装配；
  `WallParams = PotentialParams + radius` 打包传。`h > 0` 是调用方前置条件（h<0 时力反号），
  壁面势开启时越界区间收紧到 `[-½+a, Nz-½-a]`。
- 判据 W8a–W8e（数值见 README）；`wallpotential = none` 时全部旧输出逐位不变。
- `fpd_tool --sed-stats`（复用 `Analysis.cpp` 分块平均）；`--diff-ckpt` 补上 F 数组。
- 端到端（**变密度改造前**，重力为每粒子恒力 -10，等密度惯性）：单粒子 32³ 停在 z≈4.69（独立解
  `z_rest = 4.7189`）；无壁面势对照沉到 `Rz=-0.5027` 后 rc=3 中止。`sed.cfg` 128×64×32、N=95、
  20 万步：前沿 `z_min` 7.68 → 4.95，最终在 4.62–4.75 徘徊，全程最小间隙 1.887；
  `--sed-stats` 在 400 时间单位内找不到分块平台，按设计拒绝给出沉降速度。当时吞吐 116 步/秒。
  日志：`baseline/phase6_*.log`。

### 2026-09-14/15 Phase 8-A：删除三维 FFT 与 z 周期路径

删除 `solve_periodic()`、`cufftPlan3d`、`NS_Config::wall_z`、`boundary_z`、`gravity_compensate` 及全部
`if (wall_z)` 分支；`.fpd` 版本 1 → 2 取代 `FPD_FLAG_WALL_Z`，载入后直接断言上壁 `vz` 为 0。
验收：壁面检查点 5 个时刻逐位相同，各判据输出逐字节相同，重启逐位相同。
排序经验：先删 `make_ns_config` 的 `wall_z=0` 默认值才能穷举调用点；周期路径还活着时先切判据，
用两条臂对照证明切换安全；版本号最后改（否则 `--diff-ckpt` 无法读旧文件做逐位证明）。

覆盖损失：3A `--noise`（dt 偏差定律无法重建）；W5b 子判据「与周期控制同量级」（失去同机基线）；
3B `--equipart`（粒子侧 FDT 空洞，见缺口 1）。

### 2026-09-11 Phase 7-B：z 向无滑移壁面接线

分步提交：修 `--check-poisson` 丢弃返回值 → 抽取 `solve_pressure`（逐位）→ 显式 `v*` 重构（4e-16）→
壁面流体力学（`Wall.h`、棱边 `Nz+1` 层、切向 ghost、噪声 √2）→ 相场按 Loc 截断 → 粒子 z 不折叠 + 越界中止
→ 判据与文档。四个推导写入 `doc/PressurePoisson.md` §14：ghost 代入、棱边多一层、噪声 √2、壁面下能量均分
目标 `⟨Σ|v|²⟩ = kT·dim`、`dim = n_v - size + 1`。当时 W5b 为 γ=2 0.9750±0.0029 vs γ≡1 0.9427±0.0028（8.0σ），
W5'a 的 `|trS-dim|/dim` 按 dt² 下降。

**顺带修复模板盒少一格**：盒取 `[kn-range+1, kn+range]` 在负方向漏层（FACE_X @ Rz=16.5：2140 vs 暴力枚举
2152 点），先于壁面存在，力守恒抓不到；W4s 上下壁镜像判据抓出。修后 `∫φ` 三方向散布 5.7e-6 → 2.3e-7，
与连续积分差 0.008% → 4.7e-7。

### 2026-09-07 Phase 7-A：壁面压力泊松求解器

xy 批量 2D FFT + z 向 Thomas，`build_tridiag_coeffs` 预算前推系数；Neumann 由「壁面法向面不修正」导出；
奇异列 `d_0 -= 1` 定规；归一化 `Nx·Ny`。算子往返判据抓到 `thomas_sweep` z 步长写成 1 的 bug。
spike：`spike_fft2d_batch`（轴序、归一化、inembed）、`spike_tridiag`。日志 `baseline/phase7a_poisson.log`。

### 2026-09-07 Phase 5：粒子间势、`FpdState`、`fpd_check` 拆分

`FpdState` 收拢 4 份分配/映射样板，`require()` 把 present 遗漏变成启动即 abort；势函数 WCA/Morse/LJ126
与 J1–J6（J1 抓到 `s6` 写成三次方、力差 200 倍）；自检拆成 `fpd_check`；3A/3B 统一用 Philox。
spike 推翻「粒子间力在 CPU 算」的计划（`spike_pair_bench`）。发现逐位重启只在支撑域不重叠时成立
（N=8 间距 10.5 跨进程 ~1e-16 差异，potential=none 同样如此）。日志 `baseline/spike_state_and_pair.log`。

### 2026-09-05 Phase 2：配置 + `.fpd` + 逐位重启

字段表三用（解析/覆盖/dump）；C++ 只剩 `.fpd` 一条输入路径；删除 ASCII VTK；Python 侧格式镜像与 VTK 转换。
验收：N=1/N=2 逐位重启、改一字节/截断被校验和拒绝、VTK 往返相对差 0。RNG 实测推翻累计记账方案，改 slot 方案
（`spike_rng_offset/slice/advance`）。顺手修掉：F 输出前漏 `update host`、NaN 早退跳过清理、输出目录不存在
静默丢数据。踩坑：`#pragma acc update` 不能放在 lambda 里。

### 2026-09-05 测试 3A / 3B（周期盒，判据已删除）

**3A 纯流体噪声谱**（32³、η≡1、kT=1、T=600），目标 `⟨|v|²⟩ = 2kT`：

| dt | 0.02 | 0.01 | 0.005 | 0.0025 |
|---|---|---|---|---|
| 偏差 | +6.61% | +3.18% | +1.56% | +0.77% |

`bias% = 334.6·dt − 0.104`，外推截距在误差内为零：噪声标定正确，偏差是显式 Euler 的 O(dt) 项。

**3B 粒子能量均分**（32³、kT=0.25、T=2000），目标 `2kT/M − 2kT/(ρL³)`（k=0 冻结修正减 2 个自由度）：

| 模式 | η_c | dt | 比值 |
|---|---|---|---|
| ghost | 1 | 0.01 | 0.9952 ± 0.0349 |
| frozen（seed 1/2） | 50 | 0.002 | 1.0050 ± 0.0175 / 1.0072 ± 0.0150（seed 2 无平台） |
| moving | 50 | 0.002 | 0.9712 ± 0.0149（平台 b=16 过早，可疑） |

结论：O(dt) 偏差不被 η_c=50 放大（流体式预期 +0.67%，「按扩散数放大 50 倍」预期 +33.5%，被排除约 19σ），
生产可用 `dt=0.002`。`moving` 比 `frozen` 低约 2σ 未定论。日志 `baseline/phase3a_dt_scan.log`、
`phase3b_equipart.log`、`phase3b_constants.txt`。

### 2026-09-05 Phase 0/1：基线、交错网格封装与缺陷修复

Phase 0：git 化、缺陷清单、`baseline/phase0_particle_traj.txt`（修改前 19 万步单粒子轨迹，仅供定性对照）。
spike `spike_template_routine.cpp` 验证 `acc routine seq` 作用于模板函数可行。

Phase 1：新增 `Stencil.h`、`Check.h`、`PhiParams`、`make_ns_config()`（`inv_dt` 与 `W` 的唯一来源）；
`sum_phi` 拆成三套、删除 `phi_grid`；修复 C1–C5、C7–C9、C11（C6 噪声标定由 3A/3B 解决，C10 由工厂保证）：

| 编号 | 缺陷 |
|---|---|
| C1 | 力投影用 x 面心权重，归一化却用体心的 `50·φ` 求和 → 位置相关伪力 |
| C2 | 速度平均用全场 Σφ 除以单粒子归一化，N=1 恰好抵消，骗过 20 万步 |
| C3 | 粘度系数用 `ratio_eta` 而非 `ratio_eta - 1`（粒子内 η=51） |
| C4 | `Vx[n] +=` 前未清零（影响极小） |
| C5 | PBC 用 `fmod`，负数得负值 |
| C6 | `W=10, dt=0.001` 隐含 kT=0.05，噪声几乎关着 |
| C7 | `collapse(4)` 下 atomic 打同一地址，GPU 退化串行 |
| C8 | `present()` 遗漏 `sum_phi` |
| C9 | 旧索引约定的死代码 |
| C10 | `dt` 与 `inv_dt` 人工维护 |
| C11 | cuFFT/cuRAND 返回码未检查 |
