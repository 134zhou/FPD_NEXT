# 压力泊松求解器代码阅读记录

日期：2026-09-22

## 当前状态

- 记录基线：`8d6ad6c`（`重构：按离散定义域拆分应力内核`）。
- 工作区已有未提交修改：`src/Poisson.cpp` 把仅作转调的 `solve_wall` 合并为公开入口
  `solve_pressure`。这是接口精简，不改变函数体中的数值运算。
- 本轮工作只整理代码理解，没有重新构建或运行数值判据；上述工作区修改仍应执行
  `cmake --build build -j` 和 `./build/fpd_check --check-poisson`。
- 本文是便于继续阅读代码的说明；离散算子的正式定义仍以
  `doc/PressurePoisson.md` 和 `src/Poisson.cpp` 为准。

## 1. `solve_pressure` 的完整数据流

压力投影前，`Stokes.cpp` 已经计算预测速度 `v*`，然后构造

\[
b=\frac{1}{\Delta t}\,\operatorname{div}_h v^*.
\]

`solve_pressure` 中的数据流是：

```text
物理空间的 b
  -> xy 二维正向 FFT
  -> 修正 (0,0) 奇异列右端的相容性
  -> 每个 (kx,ky) 模式沿 z 做 Thomas 求解
  -> xy 二维逆向 FFT
  -> 除以 Nx*Ny，得到物理空间压力 p
```

这里只有 x、y 做了傅里叶变换。数组下标 `(i,j,k)` 在正向 FFT 后的含义是：

- `i`：x 方向傅里叶模式；
- `j`：y 方向傅里叶模式；
- `k`：实际的 z 网格层，不是 z 向傅里叶模式。

因此 `fft_data[IDX(0,0,k)]` 表示第 k 层在 xy 平面上完全均匀的分量。由于 cuFFT
正变换没有归一化，它等于该层所有 xy 网格值的和，而不是平均值。

## 2. 为什么 `(0,0,k)` 沿 k 求和必须为零

### 2.1 连续形式

对整个计算区域积分：

\[
\int_\Omega b\,dV
=\frac{1}{\Delta t}\int_\Omega \nabla\cdot v^*\,dV
=\frac{1}{\Delta t}\int_{\partial\Omega}v^*\cdot n\,dS.
\]

x、y 是周期边界，两侧通量相互抵消；z 上下是不可穿透壁面，法向速度为零。因此

\[
\int_\Omega b\,dV=0.
\]

### 2.2 离散形式

代码中的离散散度可以简写为

\[
b_{ijk}=(U_{ijk}-U_{i-1,j,k})
       +(V_{ijk}-V_{i,j-1,k})
       +(W_{ijk}-W_{i,j,k-1}).
\]

对所有网格求和时，x、y 周期方向逐项相消；z 方向也望远镜相消，只剩上、下壁面的
法向速度，而两者都为零。因此

\[
\sum_{i,j,k}b_{ijk}=0.
\]

二维 FFT 的零模式满足

\[
\hat b_{00k}=\sum_{i,j}b_{ijk},
\]

所以同一个结论可以写成

\[
\sum_k\hat b_{00k}=0.
\]

要求为零的是整列的和，不是每一个 `b_hat(0,0,k)`。例如 `[2,-1,-1]` 是合法的，
因为它保留了 z 层之间的变化，同时总和为零。

## 3. `fix_singular_column` 做什么

当 `(i,j)=(0,0)` 时，xy 拉普拉斯本征值为零，剩下的是带 Neumann 边界的一维 z
拉普拉斯。它不能产生一个空间常数的右端，因此右端必须满足

\[
\sum_k\hat b_{00k}=0.
\]

精确算术下，这个条件由离散散度和边界条件自动保证；浮点运算会留下很小的求和残差。
`fix_singular_column` 先计算

\[
s=\sum_k\hat b_{00k},
\qquad
\bar b=\frac{s}{N_z},
\]

然后执行

\[
\hat b_{00k}\leftarrow\hat b_{00k}-\bar b.
\]

`diag[0]=|s|` 保存修正前的相容性残差。正常结果应接近机器精度；若它是 `O(1)`，
说明壁面法向速度、散度离散或索引可能有错误，不能把减均值当成真正修复。

这个函数不能清空整个 `(0,0,k)` 列，因为该列随 z 的变化表示水平均匀的竖直压力，
其中包括壁面支撑粒子重量所需的静压。

## 4. 它与后面的压力投影不重复

设 `D` 是离散散度，`G` 是离散梯度：

\[
q=Dv^*.
\]

`fix_singular_column` 只把压力方程的右端投影到可解空间：

\[
\tilde q=q-\operatorname{mean}(q).
\]

它不修改速度。随后求解

\[
DGp=\frac{\tilde q}{\Delta t},
\]

并执行真正的速度投影：

\[
v^{n+1}=v^*-\Delta t\,Gp.
\]

于是

\[
Dv^{n+1}=q-\tilde q=\operatorname{mean}(q)\approx0.
\]

两者的分工是：

- `fix_singular_column`：删除一个理论上应为零、但浮点下可能非零的全局平均分量，
  使压力方程可解；
- 压力投影：用压力梯度消除所有可以由边界条件允许的速度修正所消除的局部散度。

## 5. xy FFT 后的 z 向三对角方程

对每个固定的傅里叶模式 `(i,j)`，定义

\[
\lambda_{xy}(i,j)
=2\left(\cos\frac{2\pi i}{N_x}-1\right)
+2\left(\cos\frac{2\pi j}{N_y}-1\right).
\]

需要沿 z 求解

\[
A\hat p=\hat b,
\]

即

\[
\begin{aligned}
d_0\hat p_0+\hat p_1 &= \hat b_0,\\
\hat p_{k-1}+d_k\hat p_k+\hat p_{k+1} &= \hat b_k,
\quad 1\le k\le N_z-2,\\
\hat p_{N_z-2}+d_{N_z-1}\hat p_{N_z-1} &= \hat b_{N_z-1}.
\end{aligned}
\]

主对角为

\[
d_k=
\begin{cases}
\lambda_{xy}-1,&k=0\text{ 或 }k=N_z-1,\\
\lambda_{xy}-2,&1\le k\le N_z-2.
\end{cases}
\]

上下副对角都等于 1。对于 `(i,j)=(0,0)`，`build_tridiag_coeffs` 额外执行
`d_0 -= 1`，选择 `p_0=0` 的压力规范；这和右端减均值是两个不同操作。

## 6. Thomas 公式与代码变量

前推系数为

\[
w_0=\frac{1}{d_0},
\qquad
w_k=\frac{1}{d_k-w_{k-1}}.
\]

它只依赖网格和 `(i,j,k)`，不依赖每一步变化的右端，因此由
`build_tridiag_coeffs` 在初始化时预算并保存到 `tri_w`。

前代为

\[
y_0=b_0w_0,
\qquad
y_k=(b_k-y_{k-1})w_k.
\]

回代为

\[
x_{N_z-1}=y_{N_z-1},
\qquad
x_k=y_k-w_kx_{k+1}.
\]

最终的 `x_k` 就是压力傅里叶系数 `p_hat(i,j,k)`。

| 数学量 | 代码变量或存储位置 |
|---|---|
| `(i,j)` | 外层并行循环；一个 xy 傅里叶模式 |
| `k` | 内层顺序循环；实际 z 层 |
| `b_k` | 前代前的 `fft_data[2*idx]` / `fft_data[2*idx+1]` |
| `w_k` | `tri_w[idx]` |
| `y_{k-1}` | `prev_re` / `prev_im` |
| `y_k` | `re` / `im`，随后覆盖回 `fft_data` |
| `x_k` | 回代后的 `fft_data`，即最终 `p_hat` |
| `base` | `IDX(i,j,0)=i+j*Nx` |
| `stride` | `Nx*Ny`，相邻 z 层在数组中的距离 |
| `idx` | `base+k*stride=IDX(i,j,k)` |

`fft_data` 被就地复用：

```text
b_hat  --前代覆盖-->  y  --回代覆盖-->  p_hat
```

FFT 结果是复数，但三对角矩阵系数是实数，所以实部和虚部分别解同一个方程；
`prev_re/re` 与 `prev_im/im` 共用同一组 `w_k`。

## 7. CPU 与 GPU 的实际分工

本项目用 `nvc++ -acc` 编译 OpenACC，而不是把整个 `.cpp` 文件指定为 CPU 或 GPU。
普通 C++ 代码运行在 CPU；编译器把 OpenACC 计算区抽取成 GPU kernel，并在原位置生成
CPU 端的启动代码。

| 部分 | 实际执行位置 |
|---|---|
| `build_tridiag_coeffs` 的全部循环 | CPU |
| `solve_pressure` 的控制流程 | CPU |
| `cufftExecZ2Z` API 调用 | CPU 发起 |
| cuFFT 的正、逆变换计算 | GPU |
| `fix_singular_column` 的 `acc serial` 区域 | GPU 单个逻辑线程 |
| `thomas_sweep` 的 `acc parallel loop collapse(2)` | GPU 上 `(i,j)` 并行、每条 k 串行 |
| 逆变换后的归一化循环 | GPU |

初始化时的数据过程是：

```text
CPU: build_tridiag_coeffs 写主机 tri_w
  -> acc_copyin
GPU: 保存 tri_w，供以后每一步 thomas_sweep 使用
```

`present(fft_data,tri_w)` 要求这些主机指针已经在 OpenACC 映射表中拥有对应设备内存。
`host_data use_device(fft_data)` 则取出对应的 GPU 地址，传给在 CPU 端调用的 cuFFT API。

## 8. 下一步验证

当前代码阅读没有提出新的物理算法修改。针对工作区中 `solve_wall` 包装层的精简，后续应：

```bash
cmake --build build -j
ACC_DEVICE_TYPE=nvidia ./build/fpd_check --check-poisson
```

第二项必须在 NVIDIA GPU 环境正常时执行，不能把 CPU 回退结果记作 GPU 验证。
