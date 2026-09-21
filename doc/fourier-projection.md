# 用二维傅里叶变换求解不可压NS

本页是三部曲的**第二篇**。第一篇把方程推出来了，但是压力$p$无法确定。
它没有自己的演化方程，唯一的作用是强制流体满足 $\nabla\cdot v=0$。
本页要把约束变成可解的方程，再使用二维FFT解决x和y方向。

  * 第一篇：[从 F=ma 推导 Navier–Stokes 方程（FPD 形式）](ns-equation.md)
  * 第三篇：[二维与 z 向离散化](discretization.md)

单位约定同第一篇：格距 $dx\equiv 1$、密度 $\rho\equiv 1$、溶剂粘度 $\eta_\ell\equiv 1$。

## 1 压力是个外行

第一篇推出来的动量方程与约束是

$$\frac{\partial v}{\partial t} + (v\cdot\nabla)v = -\nabla p + \nabla\cdot[\eta(\nabla v + \nabla v^{T})] + f \tag{3}$$

$$\nabla\cdot v = 0 \tag{4}$$

盯着这两行你会发现(3) 里没有 $\partial_t p$，(4) 里也没有。
压力不是状态量，它是为了让不可压条件成立而临时解出来的东西。
每一个时间步都得重新解一个 $p$，使得推进完的速度恰好无散。

## 2 投影法：先按照没有压力计算，算完再投影回无散空间

把一个时间步拆成三步：

$$v^{*} = v^{n} + dt\left(-\nabla\cdot\Pi + f\right) \tag{5}$$

$$\nabla^2 p = \frac{1}{dt}\nabla\cdot v^{*} \tag{6}$$

$$v^{n+1} = v^{*} - dt\,\nabla p \tag{7}$$
其中$\Pi = v \otimes v - \eta(\nabla v+\nabla v^{T})$
分量形式为$\Pi_{ij} = v_i v_j - \eta\left(\frac{\partial v_i}{\partial x_j} + \frac{\partial v_j}{\partial x_i} \right)$

因为不可压条件 \(\nabla\cdot v=0\)，有

$$\nabla\cdot(v\otimes v) = (v\cdot\nabla)v+v(\nabla\cdot v) = (v\cdot\nabla)v$$

所以

$$-\nabla\cdot\Pi = -(v\cdot\nabla)v + \nabla\cdot[\eta(\nabla v+\nabla v^T)]$$

代入式 (5)：

$$v^* = v^n+dt [-(v\cdot\nabla)v + \nabla\cdot\left(\eta(\nabla v+\nabla v^T)\right)+f]$$

(5) 只算对流、粘性、外力，得到一个中间速度 $v^{*}$（一般 $\nabla\cdot v^{*}\neq 0$）；
(6) 解出修正量；(7) 把它扣掉。

**为什么 (6) 一定是对的**：对 (7) 两边取散度，并强制下一时刻无散：

$$\nabla\cdot v^{n+1} = \nabla\cdot v^{*} - dt\,\nabla\cdot(\nabla p) = \nabla\cdot v^{*} - dt\,\nabla^2 p = 0$$

移项就是 (6)。所以泊松方程的右端不是新加的物理，而是「投影必须精确」这个要求的代数后果。

!!! warning ""
    **本项目的关键不变量**：(6) 里的 $\nabla^2$ 必须是「(7) 用的梯度 $\nabla$」与
    「散度检验用的 $\nabla\cdot$」的**精确复合** $div\circ grad$，三处用同一套离散。
    否则 $\nabla\cdot v^{n+1}$ 只到截断误差，投影法的全部意义就没了（离散形式见第三篇）。

## 3 Helmholtz–Hodge 分解，以及为什么这叫「投影」

任何足够光滑的矢量场 $v^{*}$ 都能唯一拆成

$$v^{*} = v + \nabla\chi, \qquad \nabla\cdot v = 0$$

其中 $v$ 无散、$\nabla\chi$ 是纯梯度（无旋）。这两部分在内积意义下**正交**，
证明只要一次分部积分：

$$\langle v, \nabla\chi\rangle = \int_{\Omega}v\cdot\nabla\chi\,d\Omega = \oint_{\partial\Omega}\chi\,(v\cdot\hat{n})\,dS - \int_{\Omega}\chi\,(\nabla\cdot v)\,d\Omega$$

  * 第二项为 0，因为 $\nabla\cdot v = 0$；
  * 第一项也为 0，因为边界上法向速度为零 —— 周期边界是逐位相消，壁面边界是
    $v\cdot\hat{n} = 0$（代码里 `vz` 在两壁被钉死）。

所以整个空间是正交直和，而「从 $v^{*}$ 里取出无散部分」就是**正交投影**，
(7) 式正是它在数值上的实现。

## 4 傅里叶变换解微分方程
**为什么是二维 FFT**：壁面加在 $z$ 方向，$x/y$ 是周期的，
所以 $x/y$ 可以用 FFT 精确求解；$z$ 不一定。
这里先定义一下我使用的傅里叶变换（虽然这不重要，但是还是讲一下）：

$$\hat{f}(k_x, k_y) = \int_{-\infty}^{\infty} \int_{-\infty}^{\infty} f(x, y) e^{-i(k_x x + k_y y)} \, dx dy$$

$$f(x, y) = \frac{1}{(2\pi)^2} \int_{-\infty}^{\infty} \int_{-\infty}^{\infty} \hat{f}(k_x, k_y) e^{i(k_x x + k_y y)} \, dk_x dk_y$$

对周期（或边界上衰减）的函数，分部积分给出求导性质

$$\mathcal{F}\left\{\frac{\partial f}{\partial x}\right\} = ik_x\hat{f}, \qquad \mathcal{F}\left\{\frac{\partial^2 f}{\partial x^2}\right\} = -k_x^2\hat{f}$$

于是拉普拉斯算子退化成纯标量乘法：

$$\mathcal{F}\left\{ \nabla^2 f \right\} = (i k_x)^2 \hat{f} + (i k_y)^2 \hat{f} = -(k_x^2 + k_y^2) \hat{f}$$

泊松方程 (6) 因此在频域里变成**逐模式的代数除法** $\hat{p} = \hat{b}/\lambda$。
这就是谱方法的全部内容：**求导 = 乘 $ik$；解泊松 = 除以本征值**。
对方程$\nabla^2 p = \frac{1}{\Delta t} \nabla \cdot v^*$两侧同时做二维傅里叶变换：


$$-(k_x^2 + k_y^2) \hat{p}(k_x, k_y) = \frac{1}{\Delta t} \mathcal{F}\left\{ \nabla \cdot v^* \right\}(k_x, k_y)$$

定义波数平方标量 $k^2 = k_x^2 + k_y^2$：

* 当 $k^2 \neq 0$ 时，直接通过代数除法解出频域压力：

$$\hat{p}(k_x, k_y) = -\frac{1}{\Delta t (k_x^2 + k_y^2)} \mathcal{F}\left\{ \nabla \cdot v^* \right\}$$


* 当 $k_x = k_y = 0$ 时，对应零频（常数压力分量），物理上压力常数不影响梯度，直接设 $\hat{p}(0, 0) = 0$。

解出 $\hat{p}$ 后，做二维傅里叶逆变换回到实空间，再带入修正步更新速度场

$$p(x, y) = \mathcal{F}^{-1}\left\{ \hat{p}(k_x, k_y) \right\}$$

$$v^{n+1} = v^* - \Delta t \nabla p$$

## 5 连续傅里叶变换离散化

在计算机中，连续变量被采样为有限网格点。为了避免速度和压力解耦引起的“棋盘振荡”，采用 **MAC 交错网格（Staggered Grid）**。

### 5.1 变量空间布置

设网格步长为 $\Delta x = \Delta y = 1$：

* 压力 $p[i, j]$ 存储在**网格中心（胞心）** $(i, j)$
* 速度 $v_x[i, j]$ 存储在**面心（face）** ，在$(i+1/2, j)$，$v_y[i, j]$在$(i, j+1/2)$。

### 5.2 离散微分算子与精确离散拉普拉斯

* **面心梯度（体心 $\to$ 面心）**：

$$(G_x p)[i, j] = p[i+1, j] - p[i, j], \quad (G_y p)[i, j] = p[i, j+1] - p[i, j]$$


* **体心散度（面心 $\to$ 体心）**：

$$(D \mathbf{u})[i, j] = (u[i, j] - u[i-1, j]) + (v[i, j] - v[i, j-1])$$


* **复合二阶拉普拉斯算子 $L = D \circ G$**：

$$(L p)[i, j] = (p[i+1, j] - 2p[i, j] + p[i-1, j]) + (p[i, j+1] - 2p[i, j] + p[i, j-1])$$

代码里右端与修正步是两个分开的 kernel —— 因为 $p$ 的梯度要读相邻线程刚写下的 $p$：

```cpp
// src/Stokes.cpp — §3b 泊松右端 b = (1/dt)·div(v*)
double div_vs = (tmp_fx[ijk] - tmp_fx[IDX(im, j, k)])
              + (tmp_fy[ijk] - tmp_fy[IDX(i, jm, k)])
              + (tmp_fz[ijk] - tmp_fz[IDX(i, j, km)]);

fft_data[ijk * 2]     = cfg.inv_dt * div_vs;   // 注意：乘 inv_dt = 1/dt
fft_data[ijk * 2 + 1] = 0.0;                   // 右端是实场
```

```cpp
// src/Stokes.cpp — §4 求解 / §5 修正步
solve_pressure(cfg, fft_data, tri_w, plan, plan_xy, p, diag);

vx[ijk] = tmp_fx[ijk] - DT * (p[IDX(ip, j, k)] - p[ijk]);
vy[ijk] = tmp_fy[ijk] - DT * (p[IDX(i, jp, k)] - p[ijk]);
if (!wz_pinned_zface(cfg, k)) {
    vz[ijk] = tmp_fz[ijk] - DT * (p[IDX(i, j, kp)] - p[ijk]);
}
```
你可能很好奇为什么临时变量的名字是 `tmp_fx`、`tmp_fy`、`tmp_fz`，而不是 `tmp_vx`、`tmp_vy`、`tmp_vz`。
~~好吧，其实是我随手起的。没有什么原因。要说原因可能是因为要把外力加上去。~~
### 5.3 离散傅里叶变换下的特征值改变

连续情况下算子的本征值为 $-(k_x^2 + k_y^2)$，而在离散有限差分下，离散差分作用于离散傅里叶基函数 $e^{i (\theta_x i + \theta_y j)}$（其中 $\theta_x = \frac{2\pi m}{N_x}, \theta_y = \frac{2\pi n}{N_y}$）：


$$(e^{i \theta_x (i+1)} - 2e^{i \theta_x i} + e^{i \theta_x (i-1)}) = 2(\cos \theta_x - 1) e^{i \theta_x i}$$

因此，离散拉普拉斯算子在频域对应的离散本征值是：


$$\lambda(m, n) = 2\left(\cos \frac{2\pi m}{N_x} - 1\right) + 2\left(\cos \frac{2\pi n}{N_y} - 1\right) = 2\left(\cos \theta_x + \cos \theta_y - 2\right)$$

### 5.4 完整的离散求解流程

每一个时间步循环执行以下步骤：

1. **显式步推进速度**：计算对流和扩散项，得到 $v_x^*, v_y^*$。
2. **计算离散散度（右端项）**：

$$b[i, j] = \frac{1}{\Delta t}\Big( (v_x^*[i, j] - v_x^*[i-1, j]) + (v_y^*[i, j] - v_y^*[i, j-1]) \Big)$$


3. **2D FFT 变换**：

$$\hat{b} = \text{FFT2D}(b)$$


4. **频域代数除法（谱求解）**：

$$\hat{p}[m, n] = \begin{cases} \dfrac{\hat{b}[m, n]}{\lambda(m, n)} = \dfrac{\hat{b}[m, n]}{2(\cos \frac{2\pi m}{N_x} + \cos \frac{2\pi n}{N_y} - 2)}, & (m, n) \neq (0, 0) \\ 0, & (m, n) = (0, 0) \end{cases}$$


5. **2D IFFT 逆变换**：

$$p = \text{IFFT2D}(\hat{p})$$


6. **修正速度场**：

$$u^{n+1}[i, j] = u^*[i, j] - \Delta t \, (p[i+1, j] - p[i, j])$$


$$v^{n+1}[i, j] = v^*[i, j] - \Delta t \, (p[i, j+1] - p[i, j])$$

修正后的速度场散度可以在机器浮点精度内严格满足不可压缩条件。




## 6 cuFFT 的两个静默陷阱

### 6.1 轴序搞反不报错，只是得到一个转置的场

cuFFT 是 row-major，`n[]` 的**最后一维变化最快**。本项目的索引约定是
$IDX(i,j,k) = i + jN_x + kN_xN_y$（$x$ 最快），所以二维变换的尺寸数组必须写$n = \{N_y,\,N_x\}$。
批量 2D 的 $z$-slab 在 `IDX` 布局下天然连续（`idist = Nx*Ny`），**不需要任何重排**。

### 6.2 归一化：3D 除 `size`，2D 批量除 `Nx*Ny`

cuFFT 不做归一化，正向 + 逆向合起来会把结果放大**变换点数**倍：

  * 3D 版是 $N_xN_yN_z$
  * 2D 版是 $N_xN_y$

```cpp
// src/State.cpp — 壁面模式的 xy 批量 2D FFT 计划
int n[2] = { cfg.Ny, cfg.Nx };      // cuFFT row-major，n[] 最后一维最快
CUFFT_CHECK(cufftPlanMany(&plan_xy, 2, n, NULL, 1, cfg.Nx * cfg.Ny,
                          NULL, 1, cfg.Nx * cfg.Ny, CUFFT_Z2Z, cfg.Nz));
```
## 7 公式 ↔ 代码对照

| 公式 | 代码 |
| --- | --- |
| (5) $v^{*} = v^{n} + dt(-\nabla\cdot\Pi + f)$ | `Stokes.cpp` §3a（`tmp_fx/tmp_fy/tmp_fz`） |
| (6) $\nabla^2 p = \nabla\cdot v^{*}/dt$ | `Stokes.cpp` §3b + `Poisson.cpp: solve_pressure()` |
| (7) $v^{n+1} = v^{*} - dt\nabla p$ | `Stokes.cpp` §5 |
| $\lambda(k) = 2(\cos\theta_x+\cos\theta_y+\cos\theta_z-3)$ | `Poisson.cpp: solve_periodic()` 的 `nrm` |
| $\lambda_{xy}(m,n) = 2(\cos\theta_x-1)+2(\cos\theta_y-1)$ | `Poisson.cpp: build_tridiag_coeffs()` 的 `lam` |
| 2D 批量归一化 $N_xN_y$ | `Poisson.cpp: solve_wall()` 的 `nrm` |
| 3D 归一化 $N_xN_yN_z$ | `Poisson.cpp: solve_periodic()` 的 `/(double)size` |
| 轴序 $n=\{N_y,N_x\}$、batch $N_z$ | `State.cpp` 的 `cufftPlanMany` |
| $(0,0)$ 奇异列的相容投影 | `Poisson.cpp: fix_singular_column()` |

## 8 实测判据

测量环境 nvc++ 26.3 + RTX 3060。
~~这些验证都是AI跑的，看个乐子完了。~~

| # | 判据 | 实测 | 阈值 |
| --- | --- | --- | --- |
| P1 | 批量 2D FFT + 手写 z-DFT 对照 3D FFT（轴序） | 1.0e-14（8×6×5）、3.1e-14（32×16×8） | < 1e-13 |
| P2 | 2D 批量 forward + inverse 的放大因子 | 精确等于 $N_xN_y$ | 等于 $N_xN_y$ |
| P3 | `inembed = NULL` 与 `inembed = n` 对照 | 0.000e+00（逐位） | 逐位 |
| P3' | Thomas 对照稠密 Gauss（`--check-tridiag`） | 2.8e-16 | < 1e-12 |
| P4 | 奇异列定规：$\vert x_0\vert$ 与 $\vert A_0x-b\vert$ | 8.2e-15 | < 1e-13 |
| P5 | 壁面算子往返（`--check-poisson`） | 4.3e-16（8×4×6）、2.7e-15（128×64×32） | < 1e-10 |
| P6 | 周期算子往返（`--check-poisson`） | 2.0e-15（128×64×32） | < 1e-10 |
| P6' | 相容性诊断（人为给不相容右端） | diag = 3.811 等于 $\vert\sum_k\hat{b}_k\vert$ | 匹配 |

!!! note ""
    **往返判据为什么最判决性**：`p_ref → 手写 div∘grad → b → solve_pressure → p`，
    输出的 $p$ 必须等于 $p_{ref}$ 加一个常数。手写算子独立于求解器
    （求解器走 FFT + Thomas，手写算子走逐点 ghost 差分），所以这不是自洽性检验。
    它抓到过一个真实 bug：`thomas_sweep` 里 z 向步长写成 1 而不是 `Nx*Ny`，
    往返误差 O(1)。

## 9 下一步

下一章我要解决z向无法使用傅里叶变换，只能离散求解的问题。
[第三篇：二维与 z 向离散化](discretization.md)
