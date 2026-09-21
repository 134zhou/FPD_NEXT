# 二维与 z 向离散化

本页是三部曲的**第三篇**：把前两篇的连续方程与谱解真正落到网格上。
前半是**二维（xy）离散** —— 交错网格、半格偏移、离散 `div`/`grad`、离散本征值与交错混叠；
后半是 **z 向离散** —— 周期本征值、壁面 Neumann 的导出、Thomas 三对角、奇异列定规、
切向 ghost、棱边 Nz+1 层、壁面噪声 √2、离散 Poiseuille 闭式。

  * 第一篇：[从 F=ma 推导 Navier–Stokes 方程（FPD 形式）](ns-equation.md)
  * 第二篇：[用二维傅里叶变换求解不可压 NS](fourier-projection.md)

单位：$dx = 1$、$\rho = 1$、$\eta_\ell = 1$
索引约定$IDX(i,j,k) = i + jN_x + kN_xN_y$，$x$ 变化最快。
## 1 MAC 交错网格

### 1.1 变量都放在哪

| 量 | 位置 | 物理坐标 |
| --- | --- | --- |
| $p$、$\eta$、$\Pi_{xx}$、$\Pi_{yy}$、$\Pi_{zz}$ | 体心 | $(i,\ j,\ k)$ |
| $v_x$、$f_x$ | x 面 | $(i+\frac{1}{2},\ j,\ k)$ |
| $v_y$、$f_y$ | y 面 | $(i,\ j+\frac{1}{2},\ k)$ |
| $v_z$、$f_z$ | z 面 | $(i,\ j,\ k+\frac{1}{2})$ |
| $\eta_{xy}$、$\Pi_{xy}$ | xy 棱 | $(i+\frac{1}{2},\ j+\frac{1}{2},\ k)$ |
| $\eta_{yz}$、$\Pi_{yz}$ | yz 棱 | $(i,\ j+\frac{1}{2},\ k+\frac{1}{2})$ |
| $\eta_{zx}$、$\Pi_{zx}$ | zx 棱 | $(i+\frac{1}{2},\ j,\ k+\frac{1}{2})$ |

**为什么要交错**：简单来说面心的差分会直接落在体心，所以最简单的差分就可以得到精确的结果，不需要计算物理课上教的三点差分。

### 1.2 半格偏移的唯一真值源

各位置相对体心的偏移由 `Loc` 枚举与 `LocOffset` 模板给出：

```cpp
// src/include/Common.h — 索引约定：x 变化最快
#define IDX(i, j, k) ((i) + (j)*Nx + (k)*Nx*Ny)
```

```cpp
// src/include/Stencil.h
// 交错网格上的位置。MAC 布局：
//   CELL              体心   —— p, eta, 对角应力 Π_xx/yy/zz
//   FACE_X/Y/Z        面心   —— vx/vy/vz, fx/fy/fz
//   EDGE_XY/YZ/ZX     棱心   —— etaXY/YZ/ZX, 非对角应力 Π_xy/yz/zx
enum Loc { CELL, FACE_X, FACE_Y, FACE_Z, EDGE_XY, EDGE_YZ, EDGE_ZX };

// 每个位置相对体心的半格偏移
template<Loc L> struct LocOffset;
template<> struct LocOffset<CELL>    { static constexpr double dx=0.0, dy=0.0, dz=0.0; };
template<> struct LocOffset<FACE_X>  { static constexpr double dx=0.5, dy=0.0, dz=0.0; };
template<> struct LocOffset<FACE_Y>  { static constexpr double dx=0.0, dy=0.5, dz=0.0; };
template<> struct LocOffset<FACE_Z>  { static constexpr double dx=0.0, dy=0.0, dz=0.5; };
template<> struct LocOffset<EDGE_XY> { static constexpr double dx=0.5, dy=0.5, dz=0.0; };
template<> struct LocOffset<EDGE_YZ> { static constexpr double dx=0.0, dy=0.5, dz=0.5; };
template<> struct LocOffset<EDGE_ZX> { static constexpr double dx=0.5, dy=0.0, dz=0.5; };
```

`IDX` 与 `cufftPlan3d(Nz, Ny, Nx)` 是配套的：cuFFT 的最后一维变化最快，而这里 $x$ 最快，
所以三维变换的尺寸数组是 $(N_z, N_y, N_x)$。这个搭配**不能动**。
~~我承认这是我写的屎山代码，可以改的但是懒得改了。~~
### 1.3 相场模板盒与球形截断

`stencil_point()`（模板参数 `Loc`）是**唯一**实现半格偏移、球形截断、周期回绕的地方，
`Viscosity` / `Force` / `Velocity` 三个模块共用它。它把局部平铺索引 $l$ 映射到全局下标
$ijk$ 与相场权重 $w$：

```cpp
// src/include/Stencil.h — 模板盒核心（半格偏移 + 球形截断 + 回绕）
const int in = (int)floor(Rnx), jn = (int)floor(Rny), kn = (int)floor(Rnz);
const int ir = li + in - pp.range;
const int jr = lj + jn - pp.range;
const int kr = lk + kn - pp.range;

const double dx = ir - Rnx + LocOffset<L>::dx;
const double dy = jr - Rny + LocOffset<L>::dy;
const double dz = kr - Rnz + LocOffset<L>::dz;

if (dx*dx + dy*dy + dz*dz > pp.range2) { return false; }   // 唯一的截断判据

const int iw = (ir + Nx) % Nx;      // x/y 周期回绕
const int jw = (jr + Ny) % Ny;
...
ijk = (zk == ZEDGE) ? wz_edge_idx(cfg, iw, jw, kr) : IDX(iw, jw, (kr + Nz) % Nz);
w   = order(dx, dy, dz, pp.radius, pp.inv_xi);
```
~~我这里临时变量的命名有点迷，但是懒得改~~
  * **模板盒必须取 $[k_n - range,\ k_n + range]$，宽度 $2\,range + 1$**。
    半径为 $range$ 的球心落在格胞**任意**位置时，能触及的整数层最多有 $2\,range+1$ 个；
    宽 $2\,range$ 的盒会在负方向**静默漏掉一层**。漏掉的权重随亚格点位置与 `Loc` 变化，
    是方向相关的伪偏差。
  * **壁面模式下必须用 `floor` 而不是 `(int)`**：
    粒子中心合法范围是$z\in(-\frac{1}{2},\,N_z-\frac{1}{2})$，$R_z$ 可以为负；
    `(int)(-0.3) = 0` 而`floor(-0.3) = -1`，用 `(int)` 会让整个模板盒偏一格。
    你可能好奇为什么是$z\in(-\frac{1}{2},\,N_z-\frac{1}{2})$，
    那是因为格子的编号是$0$到$N_z-1$，我把壁面的位置各往外推了$\frac{1}{2}$。

```text
               z 坐标                      变量与含义
────────────────────────────────────────────────────────────────
上边界 (Top)   z = Nz - 1/2  ─────── vz[Nz-1] (第 Nz-1 胞元的顶面，即上壁面)
               z = Nz - 1             p[Nz-1]  (第 Nz-1 胞元的中心)
               z = Nz - 3/2  ─────── vz[Nz-2] (内部流体面)
                  ...                  ...
               z = 1/2       ─────── vz[0]    (第 0 胞元的顶面 / 第 1 胞元的底面)
               z = 0                  p[0]     (第 0 胞元的中心)
下边界 (Bottom)z = -1/2      ─────── vz[-1]   (第 0 胞元的底面，即下壁面)
────────────────────────────────────────────────────────────────
```

### 1.4 三套相场的归一化系数要重采样

`Viscosity` 为每个粒子维护**三套**归一化因子 `sum_phix / sum_phiy / sum_phiz`，
分别用 `FACE_X / FACE_Y / FACE_Z` 的权重求和。
因为 $x$ 面、$y$ 面、$z$ 面上的相场是**同一个连续函数在三组错开半格的格点上重新采样**，所以

$$\sum\hat{\phi}_x(r) \neq \sum\hat{\phi}_y(r)$$
当然，相场和粘度场也要每个算一遍。

## 2 z向壁面如何与循环边界结合
### 2.1 我们要解什么方程组
三维拉普拉斯算子是水平和竖直的叠加：

$$\nabla^2 p = L_{xy} p + L_z p$$

方程$\nabla^2 p = b$两边同时做关于水平方向 $(x, y)$ 的傅里叶变换：

$$\mathcal{F}_{xy}\{ (L_{xy} p)_{i,j,k} + (L_z p)_{i,j,k} \} = \mathcal{F}_{xy}\{ b_{i,j,k} \}$$

$L_z$ 只涉及 $k$ 维度的变化，与水平索引 $i, j$ 互相独立，傅里叶变换对 $k$ 来说是线性算子：


$$\mathcal{F}_{xy}\{ p_{i,j,k+1} - 2p_{i,j,k} + p_{i,j,k-1} \} = \hat{p}_{m,n,k+1} - 2\hat{p}_{m,n,k} + \hat{p}_{m,n,k-1}$$

代入水平项的变换结果：


$$\lambda_{xy} \hat{p}_{m,n,k} + \left( \hat{p}_{m,n,k+1} - 2\hat{p}_{m,n,k} + \hat{p}_{m,n,k-1} \right) = \hat{b}_{m,n,k}$$

在固定的某个水平频率 $(m, n)$ 下，省略下标 $m, n$，简记频域变量 $\hat{p}_{m,n,k}$ 为 $p_k$，右端频域项为 $\hat{b}_k$：

$$p_{k-1} + (\lambda_{xy} - 2) p_k + p_{k+1} = \hat{b}_k$$

所以每个水平模式都对应一个 $Nz \times Nz$ 的三对角线性方程组 $A x = b$。
不过我们还没有推导在上下边界会发生什么。
### 2.2 导出三对角

根据修正步，速度更新公式为 $v_z^{n+1} = v_z^* - \Delta t \, G_z p$。
但在物理壁面上，流体不可穿透，壁面法向速度恒为 $0$

* 下壁面：$v_z[-1] = 0$
* 上壁面：$v_z[N_z-1] = 0$
* 内部自由面（$k = 0, \dots, N_z-2$）：执行正常修正

$$v_z[k] = v_z^*[k] - \Delta t (p[k+1] - p[k])$$

**第 $0$ 个胞元（紧贴下壁面）** 在修正步之后的离散散度：

$$(D_z v_z^{n+1})[0] = v_z^{n+1}[0] - v_z^{n+1}[-1]$$

代入：

* $v_z^{n+1}[0] = v_z^*[0] - \Delta t (p[1] - p[0])$（上表面是内部自由面，受压力梯度修正）
* $v_z^{n+1}[-1] = 0$（下表面是刚性壁面，恒等于 $0$）

相减得到：

$$(D_z v_z^{n+1})[0] = v_z^*[0] - \Delta t (p[1] - p[0]) - 0 = v_z^*[0] - \Delta t (p[1] - p[0])$$
同时$(D_zv_z^*)[0] = v_z^*[0]-v_z^*[-1] = v_z^*[0]$
所以$(D_zv_z^{n+1})[0]=(D_zv_z^*)[0]-\Delta t(p[1]-p[0])=0$

因此
$$p[1]-p[0]=\frac1{\Delta t}(D_zv_z^*)[0]$$

**第 $N_z -1$ 个胞元（紧贴上壁面）**


上壁面速度固定为$v_z^{n}[N_z-1]=0$

胞元$N_z-1$的下表面是内部面，所以
$$v_z^{n+1}[N_z-2]=v_z^*[N_z-2]-\Delta t(p[N_z-1]-p[N_z-2])$$

因此

$$
\begin{aligned}
(D_zv_z^{n+1})[N_z-1]
&=v_z^{n+1}[N_z-1]-v_z^{n+1}[N_z-2]\\
&=0-\left[v_z^*[N_z-2]-\Delta t(p[K]-p[N_z-2])\right]\\
&=-v_z^*[N_z-2]+\Delta t(p[K]-p[N_z-2]).
\end{aligned}
$$

因为$v_z^*[N_z-1]=0$，有

$$(D_zv_z^*)[N_z-1]=v_z^*[N_z-1]-v_z^*[N_z-2]=-v_z^*[N_z-2]$$

所以

$$(D_zv_z^{n+1})[N_z-1]=(D_zv_z^*)[N_z-1]+\Delta t(p[N_z-1]-p[N_z-2])\\$$

因此

$$p[N_z-2]-p[N_z-1]=\frac1{\Delta t}(D_zv_z^*)[N_z-1]$$

写成矩阵就是

$$\begin{pmatrix}
\lambda_{xy}-1 & 1 & & & \\
1 & \lambda_{xy}-2 & 1 & & \\
& \ddots & \ddots & \ddots & \\
& & 1 & \lambda_{xy}-2 & 1 \\
& & & 1 & \lambda_{xy}-1
\end{pmatrix}
\begin{pmatrix} p_0 \\ p_1 \\ \vdots \\ p_{N_z-2} \\ p_{N_z-1} \end{pmatrix}
=
\begin{pmatrix} \hat{b}_0 \\ \hat{b}_1 \\ \vdots \\ \hat{b}_{N_z-2} \\ \hat{b}_{N_z-1} \end{pmatrix} \tag{8}$$

### 5.2 Thomas 算法：前推系数可以提前算好

三对角方程 $Ax = b$，次/超对角都是 1。消元的前推系数是

$$w_0 = \frac{1}{d_0}, \qquad w_k = \frac{1}{d_k - w_{k-1}} \quad (k = 1,\ldots,N_z-1) \tag{9}$$

式子里**没有右端 $b$**，所以 $w_k$ 只由网格几何与边界决定
可以在初始化时**算一次**存进`tri_w`（$N_xN_yN_z$ 个 double，生产盒 2 MB），之后直接查表

$$y_0 = b_0w_0, \qquad y_k = \left(b_k - y_{k-1}\right)w_k, \qquad x_{N_z-1} = y_{N_z-1}, \qquad x_k = y_k - w_k x_{k+1}$$

复数右端的实部虚部用**同一组** $w_k$，在同一个 $k$ 循环里同时推进；
整个过程**完全就地**（`y` 覆盖 `b`，`x` 覆盖 `y`），零额外显存。

```cpp
// src/Poisson.cpp — build_tridiag_coeffs：lambda_xy 与 Thomas 前推系数
const double lam = 2.0 * (cos(2.0 * M_PI * i / Nx) - 1.0)
                 + 2.0 * (cos(2.0 * M_PI * j / Ny) - 1.0);
double d0 = lam - 1.0;
if (i == 0 && j == 0) { d0 -= 1.0; }        // (0,0) 奇异列的定规（见 5.3）
double w = 1.0 / d0;
tri_w[IDX(i, j, 0)] = w;
for (int k = 1; k < Nz; k++)
{
    const double dk = (k == Nz - 1) ? (lam - 1.0) : (lam - 2.0);
    w = 1.0 / (dk - w);
    tri_w[IDX(i, j, k)] = w;
}
```

```cpp
// src/Poisson.cpp — thomas_sweep：前代与回代（就地，复右端同步推进）
double prev_re = fft_data[2 * base]     * tri_w[base];
double prev_im = fft_data[2 * base + 1] * tri_w[base];
for (int k = 1; k < Nz; k++)
{
    const int idx = base + k * stride;          // stride = Nx*Ny
    const double w  = tri_w[idx];
    const double re = (fft_data[2 * idx]     - prev_re) * w;
    const double im = (fft_data[2 * idx + 1] - prev_im) * w;
    fft_data[2 * idx]     = re;
    fft_data[2 * idx + 1] = im;
    prev_re = re; prev_im = im;
}
for (int k = Nz - 2; k >= 0; k--)               // 回代 x_k = y_k - w_k * x_{k+1}
{
    const int idx = base + k * stride;
    const double w = tri_w[idx];
    fft_data[2 * idx]     -= w * fft_data[2 * (idx + stride)];
    fft_data[2 * idx + 1] -= w * fft_data[2 * (idx + stride) + 1];
}
```

### 5.3 奇异列 $(0,0)$：不能整列清零

$(m,n) = (0,0)$ 时 $\lambda_{xy} = 0$，该列矩阵 $A_0$ 是纯 Neumann 一维拉普拉斯：
两端对角 $-1$、内部 $-2$、次对角 $1$。每一行元素之和为 0，故

$$A_0\,\mathbf{1} = 0$$

即矩阵奇异，零空间是常向量。可解的条件是相容性

$$\mathbf{1}^{T}\hat{b} = \sum_k \hat{b}_k(0,0) = 0$$

它在精确算术下**自动成立**：$\hat{b}_k(0,0) = \sum_{i,j}(div\,v^{*})_{ijk}$，
其中 xy 部分在周期方向求和为零，z 部分沿 $k$ 望远镜相消：

$$\sum_k\left(v_z^{*}[k] - v_z^{*}[k-1]\right) = v_z^{*}[N_z-1] - v_z^{*}[-1] = 0 - 0 = 0$$

**这正是「壁面法向速度必须在 $v^{*}$ 阶段就钉死」的理由**（见 5.7）。浮点下残差约 1e-13，
不处理会积成 $p$ 的线性漂移。

定规取 $A' = A_0 - e_0e_0^{T}$（即把 $d_0$ 从 $-1$ 改成 $-2$），左乘 $\mathbf{1}^{T}$：

$$\mathbf{1}^{T}A'x = \mathbf{1}^{T}A_0x - x_0 = \left(A_0\mathbf{1}\right)^{T}x - x_0 = -x_0$$

而 $\mathbf{1}^{T}b = 0$，故 $x_0 = 0$；代回得 $A'x = A_0x - x_0e_0 = A_0x = b$。
也就是说 **$x$ 精确满足原方程**，且自动落在 $p_0 = 0$ 的规范上 ——
这不是罚函数、不是近似，没有引入任何误差。

!!! warning ""
    **绝不能把 $(0,0)$ 整列清零**。周期版把单点 $(0,0,0)$ 清零就够了，但壁面版
    这一列的 $z$ 结构是**水平均匀的竖直压力分布**，正是支撑粒子重量的静压。
    代码只做两件事：`build_tridiag_coeffs` 里 `d0 -= 1.0`，
    以及 `fix_singular_column` 减掉均值并记录 $\vert\sum_k\hat{b}_k\vert$ 作为诊断量。

```cpp
// src/Poisson.cpp — fix_singular_column：减均值（相容投影），并记录残差
double sr = 0.0, si = 0.0;
for (int k = 0; k < Nz; k++)
{
    sr += fft_data[2 * IDX(0, 0, k)];
    si += fft_data[2 * IDX(0, 0, k) + 1];
}
diag[0] = sqrt(sr * sr + si * si);          // 诊断量 |sum_k b_hat_k|
const double mr = sr / (double)Nz;
const double mi = si / (double)Nz;
for (int k = 0; k < Nz; k++)
{
    fft_data[2 * IDX(0, 0, k)]     -= mr;
    fft_data[2 * IDX(0, 0, k) + 1] -= mi;
}
```

### 5.4 切向无滑移的壁面通量展开式

壁面上没有 $v_x$ / $v_y$ 自由度（它们在整数 $z$ 层，壁面在半格处 $z = \mp\frac{1}{2}$）。
标准 MAC ghost 是

$$v_x[-1] = -v_x[0], \qquad v_x[N_z] = -v_x[N_z-1] \qquad (v_y\ \text{同理}), \qquad v_z[i,j,-1] = v_z[i,j,N_z-1] = 0$$

把它代入 $\Pi_{zx}$ 的一般表达式，壁面棱边（逻辑 $k=-1$）上得到：

| 项 | 壁面棱边上的值 |
| --- | --- |
| 对流 $v_{sx,z} = (v_x[k]+v_x[k+1])/2$ | $(-v_x[0]+v_x[0])/2 = 0$ |
| 对流 $v_{sz,x} = (v_z[i_p,k]+v_z[i,k])/2$ | $0$（因 $v_z[-1]=0$） |
| 粘性 $v_x[k+1]-v_x[k]$ | $v_x[0]-(-v_x[0]) = 2v_x[0]$ |

**两个对流因子同时为 0 是物理正确的**（壁面不输运动量），$\pm 2v_x$ 正是单侧导数除以半格。
代码不再在内部 kernel 中读取 ghost，而是把代入后的结果直接写进独立壁面 kernel：

```cpp
// src/Stokes.cpp — §2c 两片壁面的 YZ/ZX 棱
pi_ny[eb_bottom] = -etaZX[eb_bottom] * ( 2.0 * vx[bottom]);
pi_ny[eb_top]    = -etaZX[eb_top]    * (-2.0 * vx[top]);
pi_nx[eb_bottom] = -etaYZ[eb_bottom] * ( 2.0 * vy[bottom]);
pi_nx[eb_top]    = -etaYZ[eb_top]    * (-2.0 * vy[top]);
```

因此 §2 分成三个互不混合定义域的 kernel：体心与 XY 棱 $k\in[0,N_z-1]$、
内部 YZ/ZX 棱 $k\in[0,N_z-2]$、两片壁面棱 $k=-1,N_z-1$。

### 5.5 棱边数组为什么必须多一层

$\Pi_{yz}$ / $\Pi_{zx}$ 与 $\eta_{yz}$ / $\eta_{zx}$ 位于 z 面层，逻辑下标 $k$ 表示
$z = k+\frac{1}{2}$。壁面在 $z = \mp\frac{1}{2}$，即逻辑 $k = -1$ 与 $k = N_z-1$：

  * $k = N_z-1$（上壁棱边）**已经在数组里**；
  * $k = -1$（下壁棱边）**在数组之外**。

而这两层承载的 $v_x$ 完全不同（一个是 $-v_x[0]$，一个是 $-v_x[N_z-1]$）。
如果沿用 $N_z$ 层布局，$k=-1$ 会周期回绕到 $N_z-1$，**两者撞同一个槽位**，
于是静默用错壁面剪切。因此棱边数组在壁面模式下开 $N_z+1$ 层：

```cpp
// src/include/Wall.h — 棱边存储映射
static inline int wz_edge_idx(NS_Config cfg, int i, int j, int k)
{
    const int Nx = cfg.Nx, Ny = cfg.Ny;
    return i + j * Nx + (k + 1) * Nx * Ny;
}
```

!!! note ""
    **关键设计**：逻辑 $k$ 的物理含义与数组槽位分开，存储层始终是 $k+1$。
    内部与壁面 kernel 都通过 `wz_edge_idx` 做这一处映射，避免整体重编号。

### 5.6 壁面剪切噪声为什么是 √2

这是做错了不会报错、只会让壁面附近 $kT$ 系统性偏低的一项。

记 $G$ 为「速度到棱边应变率」的算子，$D$ 为「棱边应力到面心力」的算子
（代码里就是 `pi[k] - pi[k-1]`）。由 5.4 的 ghost：

  * 内部棱边的 $G$ 行是 $[-1,\ +1]$；
  * **壁面棱边的 $G$ 行是 $[+2]$**；
  * 而有限体积散度在壁面棱边的系数是 $1$（不是 2）。

于是可以写 $D = -G^{T}W$，其中 $W = diag(w)$，$w_{内部} = 1$、$w_{壁面} = \frac{1}{2}$
（壁面棱边只占**半个控制体** $z\in[-\frac{1}{2},0]$）。耗散算子是

$$A = D\,diag(\eta)\,G = -G^{T}diag(\eta w)G$$

它对称、半负定。噪声力协方差（噪声幅度正比于 $\sqrt{\eta\gamma}$）为

$$\left\langle ff^{T}\right\rangle = D\,diag(\eta\gamma\,2kT/dt)\,D^{T} = \frac{2kT}{dt}G^{T}diag(\eta\gamma w^{2})G$$

离散涨落耗散要求它等于 $-\frac{2kT}{dt}A = \frac{2kT}{dt}G^{T}diag(\eta w)G$，
即 $\gamma w^{2} = w$，所以 $\gamma = 1/w$，壁面棱边 **$\gamma = 2$（幅度 ×$\sqrt{2}$）**。

**独立佐证**：涨落应力的方差反比于控制体体积；壁面棱边控制体是内部的一半，
方差 ×2，幅度 ×$\sqrt{2}$。两条互不依赖的论证给同一个数。

**范围**：只有 `EDGE_YZ` / `EDGE_ZX` 的两片壁面棱边需要它。
$\Pi_{zz}$ 的壁面行是 $[+1]$（因 $v_z[-1]=0$），散度端系数也是 1，故 $w=1$，**不需要修正**；
$\Pi_{xy}$ 与 $\Pi_{xx},\Pi_{yy}$ 全在完整控制体上，同样不需要。

代码的内部棱边 kernel 直接使用 $\gamma=1$；壁面 kernel 同时处理 $k=-1$ 与
$k=N_z-1$，使用 `wall_gamma = noise_gamma1 ? 1.0 : 2.0`。这样不存在把上壁误判成
$k=N_z$ 的 off-by-one 空间。

判据 **W5b** 用 $\gamma\equiv 1$ 对照证明「$\sqrt{2}$ 是被数据选中的，不是被假设的」：
$\gamma=2$ 给 0.9750 ± 0.0029，$\gamma\equiv 1$ 给 0.9427 ± 0.0028，相差 8.0σ。

### 5.7 壁面法向面在 $v^{*}$ 阶段就钉死

```cpp
// src/Stokes.cpp — §3a 显式一步：v* = v + dt*(-div(Pi) + f)
tmp_fx[ijk] = vx[ijk] + DT * (-d_pix + fx[ijk]);
tmp_fy[ijk] = vy[ijk] + DT * (-d_piy + fy[ijk]);
tmp_fz[ijk] = vz[ijk] + DT * (-d_piz + fz[ijk]);

// 壁面法向面【永不修正】：v*z 在两壁钉死为 0。这正是 5.3 的相容性
// 条件 sum_k b_hat_k = 0 的前提。
if (wz_pinned_zface(cfg, k)) { tmp_fz[ijk] = 0.0; }
```

### 5.8 离散 Poiseuille 流：判决性的闭式解

取均匀体力 $f_x = g$、$\eta \equiv 1$、无噪声、稳态。此时 $v_x$ 与 $x,y$ 无关
⇒ 对流恒 0、$div\,v = 0$ ⇒ $p \equiv 0$，只剩粘性与外力平衡。离散方程是

$$v_x[k+1] - 2v_x[k] + v_x[k-1] = -\frac{g}{\eta} \quad (\text{内部各层}), \qquad v_x[1] - 3v_x[0] = -\frac{g}{\eta} \quad (k = 0)$$

第二条来自 ghost $v_x[-1] = -v_x[0]$。解得**精确离散闭式**

$$v_x[k] = \frac{g}{2\eta}\left[\frac{N_z^{2}+1}{4} - \left(k - \frac{N_z-1}{2}\right)^{2}\right] \tag{10}$$

!!! warning ""
    **判据必须用这个离散闭式，不能用连续抛物线** $(g/2)(N_z^2/4 - z^2)$：
    后者带 $O(1/N_z^{2})$ 的假偏差（$N_z=16$ 时 0.39%），足以让人误判 ghost 写错了。
    判据 **W3** 与这条闭式的相对差是 **0.000e+00** —— 漏掉 ghost 的 ×2
    或散度端的系数 1 都会立刻破坏它。

## 6 时间离散与稳定性

时间推进用显式 Euler：5.7 算 $v^{*}$，第二篇 (7) 式修正。显式粘性项要求扩散数满足

$$dt \le \frac{\rho\,dx^{2}}{2d\,\eta_c} \tag{11}$$

本项目 $\eta_c = 50$、$d = 3$、$dx = 1$，故 $dt \le 1/300 \approx 0.0033$。
生产用 $dt = 0.002$，正好在稳定区里。

## 7 公式 ↔ 代码对照

| 公式 / 设计 | 代码 |
| --- | --- |
| $IDX(i,j,k) = i+jN_x+kN_xN_y$ | `Common.h` 的 `IDX` 宏 |
| 半格偏移 $LocOffset$ | `Stencil.h: LocOffset<>` |
| 球形截断 $\vert r-R\vert \le range$ | `Stencil.h: stencil_point()` |
| 模板盒 $2\,range+1$ | `Common.h: make_phi_params()` 的 `n_range` |
| $\eta$ 的体心 / 棱边取值 | `Viscosity.cpp`（`eta/etaXY/etaYZ/etaZX`） |
| $\Pi$ 六个分量 | `Stokes.cpp` §2 |
| $\lambda(k) = 2(\sum\cos-3)$ | `Poisson.cpp: solve_periodic()` |
| $\lambda_{xy} = 2(\cos\theta_x-1)+2(\cos\theta_y-1)$ | `Poisson.cpp: build_tridiag_coeffs()` |
| $w_k = 1/(d_k-w_{k-1})$ | 同上（`tri_w`） |
| Thomas 前代 / 回代 | `Poisson.cpp: thomas_sweep()` |
| $(0,0)$ 列定规 $d_0-1$ | `build_tridiag_coeffs()` + `fix_singular_column()` |
| 切向 ghost 代入后的 $\pm2\eta v$ | `Stokes.cpp` §2c 壁面 kernel |
| 法向 $v_z$ 两壁钉死 | `Wall.h: wz_pinned_zface()` + `Stokes.cpp` §3a |
| 棱边 $N_z+1$ 层映射 | `Wall.h: wz_edge_idx()` / `wz_edge_nz()` |
| 壁面 $\gamma = 2$ | `Stokes.cpp` §2c 的 `wall_gamma` |
| 相场按 `Loc` 的 z 截断 | `Stencil.h: LocZKind<>` 与 `stencil_point()` |
| 离散 Poiseuille 闭式（判据） | `tests/CheckWall.cpp: w3_poiseuille()` |
| 两片壁面剪切通量展开式（判据） | `tests/CheckWall.cpp: w3b_wall_shear_flux()` |
| z 向动量收支恒等式（判据） | `tests/CheckWall.cpp: w6_momentum_budget()` |

## 8 实测判据总表

测量环境：nvc++ 26.3 + RTX 3060（double），见 `PROGRESS.md` 与
`doc/PressurePoisson.md` §14.7。

| # | 判据 | 实测 | 阈值 |
| --- | --- | --- | --- |
| W1 | 壁面模式 $\max\vert div\,v\vert / \max\vert v\vert$，含两个壁面层 | 3.4e-16（16³，20 步） | < 1e-12 |
| W2 | $v_z[:,:,N_z-1]$ 逐位为 0 | 0.000e+00 | 逐位 |
| W3 | 与离散 Poiseuille 闭式的相对差 | 0.000e+00 | < 1e-12 |
| W3b | 两片壁面的 $\Pi_{yz}/\Pi_{zx}$ 与独立展开式 | 待本轮 GPU 复验 | < 1e-13 |
| W4 | 力守恒（远离壁 / 贴下壁 / 贴上壁 / N=2 重叠）+ 均匀流场 $V_i$ | 5.3e-15 / 0.000e+00 | 1e-12 |
| W4s | 上下壁**离散镜像**位置上 `sum_phiz/sum_phix` 逐个数字相等 | 0.000e+00 | 1e-12 |
| W5b | 壁面能量均分，$\gamma\equiv 1$ 对照 | $\gamma=2$: 0.9750 ± 0.0029；$\gamma\equiv 1$: 0.9427 ± 0.0028（差 8.0σ） | > 5σ |
| W6 | z 向动量收支恒等式，逐步 | 1.3e-15 | < 1e-11 |

纯 CPU 的势函数判据（J1–J6）见 `fpd_check --check-potential`；力守恒、亚格点与交错混叠
判据见 `fpd_check --check`；泊松算子往返（P1–P6）见第二篇 §9。

## 9 已知盲区（诚实记录）

  * **算子往返判据对「算子与修正步不匹配」是盲的**。已由 W1 补上：
    $\max\vert div\,v\vert / \max\vert v\vert = 3.4\times 10^{-16}$，
    且**含 $k=0$ 与 $k=N_z-1$ 两个壁面层** —— 那正是盲区所在。
  * **W5b 的绝对量测带约 1.7% 系统偏差**（壁面 0.975、周期控制 0.983）。
    所以判据写成**差分**形式：$\gamma=2$ 必须显著比 $\gamma=1$ 更接近 1，
    且壁面与周期控制的系统偏差同量级。绝对一致需要更大的盒子 + dt 外推，未立项。
  * **模板盒缺陷（已修复）**：窄一号的盒在负方向静默漏一层，$\int\phi$ 带方向相关伪偏差；
    力守恒判据对它完全盲，判决者是暴力枚举与连续球坐标积分。
    修完后 $\int\phi$ 三方向散布从 5.659e-06 降到 2.348e-07。
  * **壁面附近没有粒子–壁面排斥势**：粒子靠初始构型远离两壁（离壁 $\ge range$），
    越界直接报错中止（不 clamp）。
  * **壁面附近各向异性判据会报 FAIL，这是设计行为**：`Analysis.cpp` 的
    `spread_int_phi2 < 1e-4` 与 `spread_subgrid < 0.01` 建立在「相场三方向等价」上；
    粒子靠近壁面时 z 向被截断 ⇒ 各向异性是**预期的**。
    检验截断对称性的正确工具是 W4s（上下壁离散镜像），不是这两个散布判据。
