# 失败的尝试：压力分裂法（Dodd & Ferrante 2014）

日期：2026-10-19。实现于提交 `a424dfd`（分裂分支、V6、W5c 分裂行）与 `0f22da7`（`.fpd` v4 压力历史），
已从主路径移除（`.fpd` 回退 v3）。本文件只作记录，原理与数据留给以后参考。

## 动机

`ratio_rho ≠ 1` 时压力方程是变系数 `∇·(a∇p) = b`（a = 1/ρ 在面心），FFT 只能作 PCG 的预条件器，
默认容差下每步约 16 次迭代，`sed.cfg` 吞吐从 134.8 降到 35.9 步/秒。分裂法想把它换成每步一次常系数直解。

## 方法

Dodd & Ferrante, *A fast pressure-correction method for incompressible two-fluid flows*,
J. Comput. Phys. **273**, 416 (2014)。把 `a∇p` 拆成常系数部分与用外推压力近似的变系数部分：

```
ρ0 = min ρ = min(1, ratio_rho)，a0 = 1/ρ0 = max a
p̂  = 2pⁿ − pⁿ⁻¹                                       （外推压力）
∇²pⁿ⁺¹ = ρ0·b + ∇·[(1 − ρ0·a)∇p̂]                       （常系数，直接复用 solve_pressure）
vⁿ⁺¹  = v* − dt·[a0∇pⁿ⁺¹ + (a − a0)∇p̂]
```

散度、梯度、Laplacian 与 PCG 路径同一套离散，所以 `∇·vⁿ⁺¹` 到机器精度；误差转移到动量方程，
来源是 `p̂` 与真压力之差。前 2 步无历史时用 PCG 建立 `pⁿ`、`pⁿ⁻¹`。逐位重启需要把压力历史写进检查点。

## 关键代码（`a424dfd` 的 `src/Stokes.cpp`）

右端的 `∇·[(1−ρ0 a)∇p̂]` 用同一份离散算子，系数换成 `alpha + beta·a`（`apply_laplacian_affine(cfg, 1, -ρ0, ...)`；
`(0, 1)` 时与 `apply_laplacian_var` 逐位相同）。

```cpp
// 压力求解（cg_r 暂存 p̂，cg_s 暂存右端修正项）
for (...) { cg_r[ijk] = 2.0 * p[ijk] - p_prev[ijk]; p_prev[ijk] = p[ijk]; }
apply_laplacian_affine(cfg, 1.0, -rho0, inv_rho_x, inv_rho_y, inv_rho_z, cg_r, cg_s);
for (...) { fft_data[2*ijk] = rho0 * fft_data[2*ijk] + cg_s[ijk]; fft_data[2*ijk+1] = 0.0; }
solve_pressure(cfg, fft_data, tri_w, plan_xy, p, diag);

// 修正步
vx[ijk] = tmp_fx[ijk] - DT * (a0 * (p[ip] - p[ijk]) + (inv_rho_x[ijk] - a0) * (cg_r[ip] - cg_r[ijk]));
vy[ijk] = tmp_fy[ijk] - DT * (a0 * (p[jp] - p[ijk]) + (inv_rho_y[ijk] - a0) * (cg_r[jp] - cg_r[ijk]));
if (!wz_pinned_zface(cfg, k))
    vz[ijk] = tmp_fz[ijk] - DT * (a0 * (p[kp] - p[ijk]) + (inv_rho_z[ijk] - a0) * (cg_r[kp] - cg_r[ijk]));
```

完整实现：`git show a424dfd -- src/Stokes.cpp src/Poisson.cpp src/main.cpp`；
压力历史入档：`git show 0f22da7 -- src/IOBin.cpp src/include/IOBin.h tools/fpd_format.py`。

## 判据结果

| 判据 | 结果 |
|---|---|
| V6：分裂一步后 `max|div v|/max|v|` | 1.3e-15（机器精度） |
| V6：a ≡ 1 时与 PCG 之差 | 3.5e-16 |
| `--pcg-noise`（32³、N=4、ratio_rho=2、kT=0.25；历史取自参考轨迹，最有利情形） | ε_E = 9.3e-2，ε_V = 1.8（门槛 1e-3） |
| W5c（4³、ρ∈[1,2]，300k 步） | 温度 1.0049±0.0035，与 PCG 同 seed 成对差 +1.9e-4 |
| `sed.cfg` 5000 步，同 seed，对 PCG 1e-14 参考 | 粒子位置最大差 8.3e-3（换 seed 差 0.41 的 2%），`<Vz>` 偏 0.55% |
| E2 单粒子（noise_on=0） | 位置差 2.4e-8 |
| 吞吐（`sed.cfg`） | 99.2 步/秒（PCG 默认 35.9，`pcg_tol=1e-6` 53.3） |
| 重启 | v4 存压力历史后 20+30 步与 50 步逐位相同 |

## 为什么失败

- 有热噪声时，右端 b 几乎 100% 是噪声，压力是逐步独立的随机场；相邻两步的压力几乎不相关，
  外推 `2pⁿ − pⁿ⁻¹` 的误差与压力本身同量级（甚至更大）。
- 误差项 `(a − a0)∇(p̂ − p)` 在 ρ = ρ0 的溶剂里为零，只在粒子内部非零，恰好落在粒子速度的平均域上，
  所以单步粒子速度误差比噪声增量本身还大（ε_V = 1.8）。
- 平衡统计（W5c 温度）偏差不大，但单步动力学保真度不合格。

## 可能适用的场合

无噪声（`noise_on=0`）或压力随时间光滑变化的确定性算例：E2 的误差只有 2.4e-8，吞吐约为默认 PCG 的 2.8 倍。
有噪声时，更合适的提速方式是放松 PCG 容差（见 README「压力求解精度」与 `fpd_calib`）。
