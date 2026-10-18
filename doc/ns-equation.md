# 从 F=ma 推导 Navier–Stokes 方程（FPD 形式）

本页是胶体悬浮液**流体粒子动力学**模拟器公式推导三部曲的**第一篇**。
目标：从牛顿第二定律出发，一步一步推出代码真正在求解的动量方程，并说明 FPD 相对
教科书 Navier–Stokes（NS）方程的**三处改写**（变粘度、涨落应力、粒子–流体耦合）。

  * 第二篇：[用（二维）傅里叶变换求解不可压 NS](fourier-projection.md)
  * 第三篇：[二维与 z 向离散化](discretization.md)

单位约定（与代码一致）：格距 $dx \equiv 1$、**溶剂**密度 $\rho_\ell \equiv 1$、
溶剂粘度 $\eta_\ell \equiv 1$；粒子/溶剂密度比 `ratio_rho`$\,=\rho_p/\rho_\ell$。
命名约定：**小写 = 流体场**（$v, p, \eta$），
**大写 = 粒子量**（$V, F, R$）。这两个约定贯穿全部源码。

## 1 记号

| 符号 | 含义 | 代码对应 |
| --- | --- | --- |
| $\mathbf{v}(\mathbf{r},t)$ | 流体速度场 | `vx/vy/vz` |
| $p(\mathbf{r},t)$ | 压力 | `p` |
| $\rho(\mathbf r)$ | 质量密度场 $\rho_\ell+(\rho_p-\rho_\ell)\sum_i\phi_i$，$\rho_\ell=1$ | `Mass.cpp` 的面心 $a=1/\rho$ |
| $\sigma$ | Cauchy 应力张量 | `pi_dx/pi_dy/pi_dz/pi_nx/pi_ny/pi_nz` |
| $\Pi$ | 代码里实际构造的动量通量张量 | 同上 |
| $\eta(\mathbf{r})$ | 粘度场 | `eta/etaXY/etaYZ/etaZX` |
| $a,\ \xi$ | 粒子半径、界面厚度 | `radius`, `xi` |
| $\phi_i(\mathbf{r})$ | 粒子 $i$ 的相场 | `order()` |
| $\mathbf{F}_i,\ \mathbf{V}_i$ | 粒子 $i$ 的受力、速度 | `Fx/Fy/Fz`, `Vx/Vy/Vz` |
| $W$ | 噪声幅度系数 $\sqrt{2kT/dt}$ | `NS_Config::W` |

推导中出现的所有量都能在上表找到落点；本页最后一节给出完整的**公式 ↔ 代码**对照表。

## 2 从 F=ma 到连续动量方程

### 2.1 流体微团的加速度：物质导数

取一个**随流体一起运动的微团**（material volume，始终包含同一批流体质量）。
微团中心的轨迹记作 $\mathbf{r}(t)$，其速度就是该处的流场速度：

$$\frac{d r}{dt} = v(r(t), t)$$

对它再求一次时间导数就得到微团的加速度。注意**不能**只对 $t$ 求偏导 ——
因为微团自身在动，$\mathbf{r}$ 也依赖 $t$，必须用链式法则：

$$a = \frac{d^2 r}{dt^2} = \frac{\partial v}{\partial t} + \left(\frac{d r}{dt}\cdot\nabla\right)v$$

即

$$\frac{Dv}{Dt} \equiv \frac{\partial v}{\partial t} + (v\cdot\nabla)v$$

这就是**物质导数**（随体导数）。

!!! note ""

    **要点**：对流项 $(\mathbf{v}\cdot\nabla)\mathbf{v}$ 不是额外假设，
    而是「微团在动」这一事实经过链式法则后的必然产物。它描述流体把自身速度场
    **平流输运**的效果：速度被流体带着走，所以局部加速度里多出一项。
    因为v是速度场，是绑定在位置网格上的。
    如果只是 $\frac{\partial v}{\partial t}$你会发现这一秒的速度与上一秒的速度是不同微团的速度。

### 2.2 微团的牛顿第二定律

微团的质量是 $\rho\,dV$。作用在它上面的力分两类：

  * **体力**：正比于体积，$fdV$。本项目里就是粒子受到的的力被分摊到了流体上（见 §3.3）。
  * **面力**：也就是液体压强。由 Cauchy 应力张量 $\sigma$ 描述：$\sigma\cdot\hat{n}\,dS$，$\hat{n}$ 是外法向。

于是逐字写出 $F=ma$：

$$\rho\,dV \cdot \frac{Dv}{Dt} = \mathbf{f}\,dV + \oint_{\partial\Omega} \sigma\cdot\hat{n}\;dS$$

用高斯散度定理把面积分变成体积分：

$$\oint_{\partial\Omega} \sigma\cdot\hat{n}\;dS = \int_{\Omega} \nabla\cdot\sigma\;dV$$

如果这里的 $\sigma$ 是个矢量会让人感觉更熟悉一点，但是注意这里它是张量。其实展开之后可以推出来是对的。
两边同除 $dV$ 并取微团趋于一点的极限，得到**连续介质的运动方程**：
$$\int_{\Omega} \nabla\cdot\sigma\;dV \approx \nabla\cdot\sigma\;dV$$

$$\rho\left(\frac{\partial v}{\partial t} + (v\cdot\nabla)v\right) = \nabla\cdot\sigma + \mathbf{f} \tag{1}$$

### 2.3 本构关系：牛顿流体

(1) 对各向同性不可压缩流体，应力张量取为：

$$\sigma = -pI + \eta(\nabla v + \nabla v^T)$$

$-pI$是静止流体只承受各向同性的静压 $p$，这个好理解，我都能看懂，就不解释了。

$\eta(\nabla v+\nabla v^T)$：粘滞阻力部分，这是为什么呢？
速度梯度矩阵 $\nabla v$的完整形式为：

$$\nabla v = \begin{pmatrix}  \frac{\partial v_x}{\partial x} & \frac{\partial v_x}{\partial y} & \frac{\partial v_x}{\partial z} \\ \frac{\partial v_y}{\partial x} & \frac{\partial v_y}{\partial y} & \frac{\partial v_y}{\partial z} \\ \frac{\partial v_z}{\partial x} & \frac{\partial v_z}{\partial y} & \frac{\partial v_z}{\partial z} \end{pmatrix}$$

任何方阵都可以严格唯一地拆分为对称矩阵 $E$ 与反对称矩阵 $\Omega$：

$$\nabla v = E + \Omega, \quad \text{其中 } E = \frac{1}{2}\left(\nabla v + \nabla v^T\right), \; \Omega = \frac{1}{2}\left(\nabla v - \nabla v^T\right)$$

$$\Omega
= \frac{1}{2}
\begin{pmatrix}
0 & (\frac{\partial v_x}{\partial y} - \frac{\partial v_y}{\partial x}) & (\frac{\partial v_x}{\partial z} - \frac{\partial v_z}{\partial x}) \\
(\frac{\partial v_y}{\partial x} - \frac{\partial v_x}{\partial y}) & 0 & (\frac{\partial v_y}{\partial z} - \frac{\partial v_z}{\partial y}) \\
(\frac{\partial v_z}{\partial x} - \frac{\partial v_x}{\partial z}) & (\frac{\partial v_z}{\partial y} - \frac{\partial v_y}{\partial z}) & 0
\end{pmatrix}
= \frac{1}{2}
\begin{pmatrix}
0 & -w_z & w_y \\
w_z & 0 & -w_x \\
-w_y & w_x & 0
\end{pmatrix}$$
也就是说：
$$\Omega dx = \frac{1}{2} \omega \times dx$$
也就是说$\Omega$掌管整体选择，而整体旋转是不会产生粘滞阻力的。
讲回到$\eta(\nabla v + \nabla v^T)$
至于它为什么是粘滞阻力，我们可以用一个特例说明。对于一个经典的层流场景，只有$x$方向有速度，存在$\frac{\partial v_x}{y}$
那么
$$\eta(\nabla v + \nabla v^T) =
\eta \begin{pmatrix}
0 & \frac{d v_x}{d y} & 0 \\
\frac{d v_x}{d y} & 0 & 0 \\
0 & 0 & 0
\end{pmatrix}$$
对于$(0,1,0)$面，是经典的$\eta \frac{d v_x}{d y}$，$(0,0,1)$面没有力，这很正常。也许你会惊喜地发现$(1,0,0)$面竟然有力，其实这是对的，不然流体微团就要力矩不平衡了。


### 2.4 不可压约束，以及为什么 $\eta$ 不能提出散度

**不可压条件**

$$\nabla\cdot\mathbf{v} = 0$$
把本构代回 (1)。密度的处理是本页与代码最容易脱节的一处：

- **论文的原始假设**（[SM] 2018 式 (3)）：液体与粒子密度相同，$\rho$ 是常数。
  取 $\rho \equiv 1$ 后得到教科书形式的 (2)：

$$\frac{\partial \mathbf{v}}{\partial t} + (\mathbf{v}\cdot\nabla)\mathbf{v} = -\nabla p + \nabla\cdot\left[\eta\left(\nabla\mathbf{v} + \nabla\mathbf{v}^{T}\right)\right] + \mathbf{f} \tag{2}$$

- **本实现（ratio_rho ≠ 1）**：粒子相场把密度也变成场，
  $\rho(\mathbf r)=\rho_\ell+(\rho_p-\rho_\ell)\sum_i\phi_i$（$\rho_\ell=1$），
  惯性项保留 $\rho$，压力投影因此变成变系数：

$$\rho\left(\partial_t\mathbf v+(\mathbf v\cdot\nabla)\mathbf v\right)
= -\nabla p+\nabla\cdot\left[\eta(\nabla\mathbf v+\nabla\mathbf v^{T})\right]+\mathbf f+(\rho-\rho_\ell)\mathbf g$$

$$\nabla\cdot\left(\frac1\rho\nabla p\right)=\frac1{dt}\nabla\cdot\mathbf v^*,\qquad
\mathbf v=\mathbf v^*-dt\,\frac1\rho\nabla p \tag{2'}$$

重力的均匀项 $\rho_\ell\mathbf g$ 被静压 $\nabla p_h$ 解析减掉，所以体密度只剩
$(\rho-\rho_\ell)\mathbf g$（Boussinesq 拆分，代数上精确）。$\rho\equiv1$ 时 (2')
逐位退化成 (2) 与原常系数泊松。离散与 PCG 见 `doc/PressurePoisson.md` §16。

**这是本页最重要的一步，也是本项目与教科书常粘度 NS 方程的分水岭。**
因为 FPD 的粘度是空间场

$$\eta(\mathbf{r}) = \eta_\ell + (\eta_c - \eta_\ell)\sum_i \phi_i(\mathbf{r})$$

它在界面层内从 $\eta_\ell$ 变到 $\eta_c$（本项目$\eta_c/\eta_\ell = 50$），
所以$\nabla\cdot[\eta(\nabla v+\nabla v^T)]$**不能**化简成 $\eta\nabla^2 v$
代码 `src/Stokes.cpp` 正是因此逐点构造通量张量

不可压条件下的压力 $p$ **没有自己的演化方程**，
唯一的作用是约束 $\nabla\cdot v=0$

## 3 FPD 的三处改写

### 3.1 相场与粘度场

把刚性胶体粒子 $i$（中心 $\mathbf{r}_i$）表示成一个光滑相场（叫它相场是因为它区分了液相和固相）：

$$\phi_i(\mathbf{r}) = \frac{1}{2}\left[\tanh\left(\frac{a - \vert\mathbf{r}-\mathbf{r}_i\vert}{\xi}\right) + 1\right]$$

  * 粒子内部 $\vert\mathbf{r}-\mathbf{r}_i\vert \ll a$：$\phi \to 1$。
  * 远离粒子：$\phi \to 0$。
  * $\xi$ 是界面厚度：$\xi/a \to 0$ 时趋于硬球。

粘度场由各粒子相场叠加得到
$\eta(r) = \eta_\ell + (\eta_c-\eta_\ell)\sum_i\phi_i$。
粒子的「刚性」由大粘度差近似保证：粒子内部粘度是溶剂的 50 倍，内部流场被迅速
动量扩散抹平，宏观上表现得像一个不可变形的流体粒子。

```cpp
// src/include/Common.h — 相场本身（唯一实现）
#pragma acc routine seq
static inline double order(double dx, double dy, double dz,
                           double radius, double inv_xi)
{
    return 0.5*(tanh((radius - sqrt(dx*dx + dy*dy + dz*dz))*inv_xi) + 1.);
}
```

```cpp
// src/Viscosity.cpp — 粘度增量系数
// eta = eta_l + (eta_c - eta_l) * sum(phi)，eta_l = 1  ⇒  系数是 ratio_eta - 1
const double d_eta = pp.ratio_eta - 1.0;
...
#pragma acc atomic update
eta[ijk] += d_eta * w;          // 体心（正应力位置）
```

!!! warning ""

    `d_eta = ratio_eta - 1` 里的「减 1」是$\eta_c-\eta_\ell$

### 3.2 涨落应力：热噪声

**热力学一致性要求**：噪声不是随便加的，必须满足涨落耗散定理（FDT）。
在 FPD 里噪声以**涨落应力张量** $\tilde{\sigma}$ 的形式直接加在流体（含粒子内部）上，
而不是加在粒子质心上，因为后者无法正确重现流体力学记忆效应。
或者这么说：粒子的运动不是独立的，会通过流体介导，使得相邻粒子运动具有相关性。
所以直接给粒子加随机热噪声是不合理的，应该把热噪声加到整个流体场上。

$$\left\langle \tilde{\sigma}_{ij}(r,t)\,\tilde{\sigma}_{mn}(r',t')\right\rangle = 2kT\,\eta(r)\left(\delta_{im}\delta_{jn} + \delta_{in}\delta_{jm} - \frac{2}{3}\delta_{ij}\delta_{mn}\right)\delta(r-r')\,\delta(t-t')$$


这个公式翻译成人话就是
切应力（非对角项）：
$$\left\langle \tilde{\sigma}_{12}^2 \right\rangle = 2kT\eta \times 1$$
正应力（对角项）：
$$\left\langle \tilde{\sigma}_{11}^2 \right\rangle = 2kT\eta \times \frac{4}{3}$$
对角项之间的互相关：
$$\left\langle \tilde{\sigma}_{11}\tilde{\sigma}_{22} \right\rangle = -\frac{2}{3} (2kT\eta)$$
负相关是正常的。毕竟还是无散的，一个方向压缩另一个方向就要拉伸。
我们提取公因式$\sqrt{\frac{2kT}{dt}}$：
```cpp
// src/include/Common.h — W 的唯一真值源（禁止手工维护）
cfg.W = noise_on ? sqrt(2.0 * kT / dt) : 0.0;
```
这里有$dt$是因为离散化之后$\delta(t-t') = \frac{1}{dt}$。
这时候我们发现一个很大的问题，想把噪声加上去还是很麻烦。
因为计算机生成独立的噪声很容易，但是要生成相关噪声就很难，而且不利于并行化。
于是我把代码写成了这样：

```cpp
// src/Stokes.cpp — §2 正应力三分量：对流 - 粘性 - 噪声
pi_dx[ijk]  = vax * vax;                                              // 对流
pi_dx[ijk] -= eta[ijk] * 2.0 * (vx[ijk] - vx[IDX(im, j, k)]);        // 粘性 2*eta*E_xx
pi_dx[ijk] -= sqrt(2.0 * eta[ijk]) * cfg.W * randD[0 * size + ijk];   // 噪声（对角：sqrt(2*eta)）
```

```cpp
// src/Stokes.cpp — §2 切应力（棱边）：对角与棱边幅度相差 sqrt(2)
pi_nz[ijk]  = vsx_y * vsy_x;
pi_nz[ijk] -= etaXY[ijk] * ((vx[IDX(i, jp, k)] - vx[ijk]) + (vy[IDX(ip, j, k)] - vy[ijk]));
pi_nz[ijk] -= sqrt(etaXY[ijk]) * cfg.W * randN[0 * esize + ijk];      // 噪声（棱边：sqrt(eta)）
```
你会发现一个系数是1没有问题，另一个却是2而不是$\frac{4}{3}$，而且还全是独立的。
这是因为之后还有一个重新投影保证无散的步骤。
可以理解为把2中减了一个$\frac{2}{3}$出来，并且减到了互相关头上。
或者更严谨的：
原本是
$$
\begin{pmatrix}
2 & 0 & 0 \\
0 & 2 & 0 \\
0 & 0 & 2
\end{pmatrix}$$
平均迹方差为 $\langle (s_1 + s_2 + s_3)^2 \rangle = 2 + 2 + 2 = 6$。
无散化时$s'_1 = s_1 - \bar{s} = s_1 - \frac{s_1 + s_2 + s_3}{3} = \frac{2}{3}s_1 - \frac{1}{3}s_2 - \frac{1}{3}s_3$
于是
$$\langle {s'_1}^2 \rangle = \left\langle \left( \frac{2}{3}s_1 - \frac{1}{3}s_2 - \frac{1}{3}s_3 \right)^2 \right\rangle = \frac{4}{9}\langle s_1^2 \rangle + \frac{1}{9}\langle s_2^2 \rangle + \frac{1}{9}\langle s_3^2 \rangle = \frac{4}{3}$$

$$\langle s'_1 s'_2 \rangle = \left\langle \left( \frac{2}{3}s_1 - \frac{1}{3}s_2 - \frac{1}{3}s_3 \right) \left( -\frac{1}{3}s_1 + \frac{2}{3}s_2 - \frac{1}{3}s_3 \right) \right\rangle = -\frac{2}{9}\langle s_1^2 \rangle - \frac{2}{9}\langle s_2^2 \rangle + \frac{1}{9}\langle s_3^2 \rangle = -\frac{2}{3}$$

然后就变成了符合公式的
$$\begin{pmatrix}
\frac{4}{3} & -\frac{2}{3} & -\frac{2}{3} \\
-\frac{2}{3} & \frac{4}{3} & -\frac{2}{3} \\
-\frac{2}{3} & -\frac{2}{3} & \frac{4}{3}
\end{pmatrix}$$
### 3.3 力的注入与粒子速度：相场加权平均

粒子 $i$ 受到合力 $\mathbf{F}_i$（粒子间势 + 壁面力；重力已改为体力，不在 $\mathbf F_i$ 里）。
FPD 的做法是按相场权重摊开成体力密度：

$$f(r) = \sum_i \frac{\phi_i(r)\,\mathbf F_i}{\int\phi_i\,dr}$$

**分母是每粒子自己的 $\int\phi_i$，不是全场 $\sum_i\phi_i$。** 本项目的旧文档写成
$\phi F_i/\sum_i\phi_i$，与代码相反（那是被修掉的 C1/C2 类缺陷：分子分母必须用同一个
`Loc`，否则 $\sum_{\rm grid} f \neq \sum_n F_n$，产生位置相关的伪力）。

反过来，粒子速度定义为流场在相场支撑域上的加权平均：

$$V_i = \frac{\int \mathbf{v(r)}\,\phi_i(r)\,dr}{\int\phi_i(r)\,dr}$$

由 [SM] 2018 的界面修正，粒子的质量、修正因子与等效质量是

$$M_i = \int\rho\,\phi_i\,dr,\qquad
\lambda^{T} = \frac{\int\phi_i}{\int\phi_i^2},\qquad
M_{\rm eq} = \frac{\left(\int\phi_i\right)^2}{\int\phi_i^2} = \lambda^{T}M_i,\qquad
M_{\rm eff} = \frac{3}{2}M_i$$

$\lambda^{T}$ 修正「有限厚度界面被当作刚体」引入的偏差，$\xi/a\to 0$ 时 $\lambda^{T}\to 1$。
**代码的力投影没有乘 $\lambda^T$**：均匀体力下局部「力→加速度」的离散比值是
$S_1^2/\int(\phi^2/\rho)\approx ratio\_rho\cdot M_{\rm eq}$，不是 $M_i$。旧表把
$M_{\rm eq}$ 标成 $M_i$（`Analysis.cpp` 的 `mass`），本次已拆开；$\lambda^T$ 缺口只测量、
不修，它影响瞬态与布朗速度，不影响沉降的稳态力平衡。
这些量**必须用与模拟相同的离散求和**计算：

**实测**（`fpd_check --lambda 32 2.0`，纯 CPU，见 `PROGRESS.md`）：

| 常量 | x / y | z |
| --- | --- | --- |
| $\int\phi$ | 170.3061 | 170.3060 |
| $\int\phi^2$ | 100.841 | 100.832 |
| $\lambda^{T}$ | 1.688855 | 1.689010 |
| $M_i=\int\rho\phi$（ratio_rho=2） | 271.1472 | 271.1379 |
| $M_{\rm eq}=S_1^2/S_2$ | 287.6223 | 287.6486 |

$M_{\rm eff}=1.5M_i=406.72$（ratio_rho=2）。三方向 $\sum\phi_\alpha^2$ 散布 9.185e-05，
4³ 亚格点扫描散布 0.0348% —— 后者是整个方法的**系统误差下限**。

```cpp
// src/Force.cpp — 力投影：分子用哪个 Loc，分母就必须用哪个 Loc
if (stencil_point<FACE_X>(pp, cfg, Rx[n], Ry[n], Rz[n], l, ijk, w))
{
    #pragma acc atomic update
    fx[ijk] += Fx[n] * w / sum_phix[n];
}
```

```cpp
// src/Velocity.cpp — 粒子速度 = int(v*phi) / int(phi)，与上面同一个 Loc
if (stencil_point<FACE_X>(pp, cfg, Rx[n], Ry[n], Rz[n], l, ijk, w)) { sx += vx[ijk] * w; }
...
Vx[n] = sx / sum_phix[n];
```

!!! warning ""

    **绝不能用解析球体 $\frac{4}{3}\pi a^3$ 代替离散求和**：
    真实 $\int\phi$ 比它大 **24.1%**，那是扩散界面的曲率项（$+8\pi a\xi^2\pi^2/24$ 量级）。
    误用会让 $\int\phi$、$M_i$、$M_{\rm eq}$ 一起差 24%，直接毁掉验证。

## 3.4 密度对比：从等密度到 $\rho(\mathbf r)$

论文 [SM] 2018 在式 (3) 前写得很清楚："we assume that the density of the liquid is
the same as that of the particles"。本项目原先把这条假设实现成了 $\rho\equiv1$：
流体的惯性密度是 1，粒子的惯性也只来自流体（$M\approx V$），而重力却绕过质量直接
加到每个粒子上。这在本项目要研究的「粒子与溶液密度不同」的体系里是不自洽的。

本实现把它换成

$$\rho(\mathbf r)=1+(ratio\_rho-1)\sum_i\phi_i(\mathbf r),\qquad a(\mathbf r)=\frac1{\rho(\mathbf r)}$$

三件事同时改变，且必须同时改：

1. **惯性**：动量方程左端保留 $\rho$。显式步必须把对流项和力密度分开：
   $$\mathbf v^*=\mathbf v+dt\left[-(\mathbf v\cdot\nabla)\mathbf v
   +\frac1\rho\left(\nabla\cdot\tau+\nabla\cdot\tilde\sigma+\mathbf f\right)\right]$$
   压力投影变成 $\nabla\cdot(a\nabla p)=\frac1{dt}\nabla\cdot v^*$，
   $v=v^*-dt\,a\nabla p$。$a=1/\rho$ 在面心，与力的摊布、$V_i$ 的归一化用同一个
   `Loc`。**$1/\rho$ 只乘力密度（粘性、噪声、外力、压力梯度），不乘对流加速度**：
   $\rho Dv/Dt$ 两边除以 $\rho$ 后对流项在左边、系数是 1。把整个 $\nabla\cdot\Pi$
   都乘 $1/\rho$ 会得到 $\rho\partial_t v+\nabla\cdot(vv)=\dots$，既不是非守恒
   形式也不是守恒形式。这是 `Stokes.cpp` §3a 把总通量散度拆成 `d_adv`/`d_visc`
   的原因；判据是 `--check-poisson` 的 V0。这是 `Mass.cpp` 与 `Poisson.cpp`
   变系数路径的职责。
2. **重力**：改为体力 $\rho\mathbf g$。均匀项 $\rho_\ell\mathbf g$ 被静压
   $\nabla p_h$ 解析减掉（代数上精确；变系数投影无法精确吸收均匀体力，所以必须
   显式减），净体密度是 $(\rho-\rho_\ell)\mathbf g$。粒子受到的净重力自动是
   $(\rho_p-\rho_\ell)\int\phi\,\mathbf g$ —— 浮力由密度差给出。
3. **粒子惯性**：粒子速度仍定义为 $V_i=\int v\phi_i/\int\phi_i$；它的惯性来自
   粒子内部流体的 $\rho_p$，所以密粒子加速更慢、布朗速度更小（FDT 里噪声也随
   $\rho$ 缩放）。$\lambda^T$ 缺口仍在（见 §3.3），本次只测量、不修。

$ratio\_rho=1$ 时 $\rho\equiv1$、$a\equiv1$，变系数路径逐位退化成原常系数求解器
（`--check-poisson` 的 V0–V5、`--check-mass` 的 MA1–MA7）。
**边界**：$\nabla\cdot v=0$ 与「密度场由粒子位置给定」意味着质量守恒的残差是
$\partial_t\rho+v\cdot\nabla\rho=O(\xi/a)$，与 FPD 的界面近似同阶；这是对论文
等密度假设的扩展，不是论文原定理。$\lambda^T$、粒子侧 FDT 仍是独立缺口。

## 4 与 [PRL] 2000 的差异

本项目以 **[SM] 2018** 为准，[PRL] 2000 是简化版：

| 项 | [PRL] 2000 | [SM] 2018 / 本实现 |
| --- | --- | --- |
| 对流项 $(\mathbf{v}\cdot\nabla)\mathbf{v}$ | 忽略 | 保留 |
| 热噪声 | 无 | 涨落应力 + FDT，$W=\sqrt{2kT/dt}$ |
| 力密度归一化 | 未归一化 | 除以 $\int\phi_i$ |
| 转动自由度 | 无 | 无 |
| 界面修正因子 $\lambda^{T}, M_i$ | 无 | **只测量，未进力投影**（[SM] 第 2.1–2.2 节；本实现的缺口） |
| 密度 $\rho$ | 常数（=1） | 场 $\rho(\mathbf r)$，`ratio_rho`（对论文的扩展） |

## 5 公式 ↔ 代码对照

| 公式 | 代码 |
| --- | --- |
| $\frac{D v}{Dt} = \partial_t v + (v\cdot\nabla) v$ | `Stokes.cpp` §2 的 `vax*vax`、`vsx_y*vsy_x` 等对流项 |
| $\Pi = vv - \eta(\nabla v+\nabla v^{T}) - \tilde{\sigma}$ | `Stokes.cpp` §2（`pi_dx/pi_dy/pi_dz/pi_nx/pi_ny/pi_nz`）；单位密度通量 |
| $\phi_i = \frac{1}{2}[\tanh((a-\vert r- r_i\vert)/\xi)+1]$ | `Common.h: order()` |
| $\rho = 1+(ratio\_rho-1)\sum\phi_i$、$a=1/\rho$ | `Mass.cpp: update_mass_fields()`（面心） |
| $\eta = \eta_\ell + (\eta_c-\eta_\ell)\sum\phi_i$ | `Viscosity.cpp`（`d_eta = ratio_eta − 1`） |
| $W = \sqrt{2kT/dt}$ | `Common.h: make_ns_config()` |
| $\langle\tilde{\sigma}\tilde{\sigma}\rangle$ 的 FDT | `Stokes.cpp` §2 的 `sqrt(2*eta)` 与 `sqrt(eta)` 噪声项 |
| $f = \sum_i\phi_i F_i/\int\phi_i$ | `Force.cpp: update_force_field()` |
| $V_i = \int v\phi_i/\int\phi_i$ | `Velocity.cpp: update_particle_velocity()` |
| $\mathbf f_g=(\rho-1)\mathbf g$ | `Mass.cpp: apply_gravity()` |
| $\nabla\cdot(a\nabla p)=b$ + PCG | `Poisson.cpp: apply_laplacian_var` / `solve_pressure_var` |
| $v^*=v+dt\,[-\nabla\cdot(vv)+a(-\nabla\cdot\Pi_{\rm visc}+f)]$、$v=v^*-dt\,a\nabla p$ | `Stokes.cpp` §3a / §5（1/ρ 只乘力密度） |
| V0：对流项不随 $a$ 缩放 | `CheckPoisson.cpp: v0_convection_independent_of_a` |
| $\lambda^{T}, M_i, M_{\rm eq}, M_{\rm eff}$ | `tests/Analysis.cpp: compute_fpd_constants()` |

## 6 下一步

**第二篇** [用（二维）傅里叶变换求解不可压 NS](fourier-projection.md) 说明如何用
投影法把它变成压力泊松方程，再在傅里叶空间里一步解掉。
**第三篇** [二维与 z 向离散化](discretization.md) 说明这些算子在 MAC 交错网格上
到底怎么写、为什么必须这么写。
