# CLAUDE.md

给在本仓库工作的 Claude Code 实例的指引。

本文只写**怎么在这里干活**。物理背景、缺陷清单、验证结果见 `README.md`；
进度与决策记录见 `PROGRESS.md`。

## 构建与测试

```bash
cmake -B build && cmake --build build -j     # nvc++ 26.3 + CUDA 13.1 + CMake 4.2.3
./build/main --check                          # 自检：力守恒 / 亚格点不变 / 混叠 / N=2 重叠
./build/main --lambda [L]                     # 常量表 λ^T, M_i（纯 CPU，秒级）
./build/main --noise [L] [dt] [kT] [steps]    # 测试 3A：纯流体噪声谱
./build/main --equipart <ghost|frozen|moving> [L] [dt] [kT] [steps] [seed]   # 测试 3B
./build/main                                  # 生产运行（当前 128x64x32, N=1, 20 万步）
cmake --build build -j 2>&1 | grep -i error   # 编译错误
nvc++ -acc -O3 -Isrc/include -Minfo=accel -c src/X.cpp -o /tmp/x.o   # 看 kernel 生成
```

环境：RTX 3060 12 GB。实测吞吐：**153 步/秒**（128×64×32）、**1128 步/秒**（32³）。
注意 32³ 只有 1/8 的格点却只快 7.4 倍 —— 那个规模已接近**启动延迟主导**。
排期前先按实测数算机时；Phase 4 的 L=192 立方盒是 27 倍网格量。

> ⚠️ **别在后台运行时重建** —— 正在执行的二进制被写会触发 `ETXTBSY`。
> 长跑输出用 `std::endl` 或 `fflush`，否则重定向到文件时会被缓冲，只能盲等。

## 不可违反的约束

### 1. `Stencil.h` 是交错网格的唯一真值源

`src/include/Stencil.h` 的 `stencil_point<Loc L>()` 是**唯一**实现半格偏移约定、
球形截断判据、周期回绕的地方。`Viscosity` / `Force` / `Velocity` 三个模块必须共用它。

**不要在任何地方重新手写模板盒遍历循环。** 缺陷 C1/C2 的成因正是这段循环被复制粘贴了
三份然后各自漂移 —— 力投影用 x-face 的权重、归一化却用 cell-center 的求和。

规则一句话：**做投影用的权重函数、和做归一化的求和，必须是同一个函数在同一组点上的求和。
三个方向各一套，不能共用。**

### 2. 改完必须重跑 `--check`

动过 `Stencil.h` / `Viscosity.cpp` / `Force.cpp` / `Velocity.cpp` 之后**必须**跑
`./build/main --check`，全绿才算数。`stencil_point` 是单点故障，它错了三个模块一起错。

### 3. 两个测试盲区（都真实骗过人）

- **力守恒判据对 `stencil_point` 的内部公式错误是盲的** —— 投影和归一化会一起错，
  `Σf == F` 照样成立。需要独立数值对照，见 `spike/spike_template_routine.cpp`
  （用独立写的 Python 复核过 `∫φ = 170.3051225223`）。
- **单粒子测试对一整类缺陷是盲的** —— C2 在 N=1 时数学上恰好抵消，骗过了 20 万步模拟。
  **任何涉及相场加权的改动，都必须补 N=2 重叠粒子的用例。**
  判决性判据：均匀流场 `v ≡ 1` 下 `V_i` 必须精确等于 1，与重叠无关。

### 4. 旧代码不是真值

`doc/old_code/` 仅供参考。已确认它在四处是错的或有瑕疵（cellList 的 `nnIndex[12]`、
cell 数不足 3、z 向边界 hack、重力硬编码），详见 `README.md`。
「和旧代码对齐」只能作定性 sanity check，**定量真值来自文献公式和自洽性测试**。

## 代码约定

- **命名**：小写 = 流体场（`vx`、`fx`、`eta`），大写 = 粒子量（`Vx`、`Fx`、`Rx`）。贯穿全代码，别破坏
- **注释用中文**，物理术语可混排英文。函数 `snake_case`，文件名 `PascalCase.cpp`，Allman 括号
- **索引** `IDX(i,j,k) = i + j*Nx + k*Nx*Ny`，x 最快。与 `cufftPlan3d(Nz,Ny,Nx)` 匹配，**别动**
- **单位**：格距 `dx ≡ 1`，密度 `ρ ≡ 1`，溶剂粘度 `η_ℓ ≡ 1`
- 风格接近「带 vector 的 C」：裸指针 + OpenACC 手动映射，不用 class/智能指针/异常。
  新代码沿用，别引入现代 C++ 抽象（`Stencil.h` 的模板是唯一例外，且已验证不影响 kernel 生成）
- POD 结构体（`NS_Config`、`PhiParams`、未来的 `PotentialParams`）**按值传进 GPU 函数**。
  `acc routine seq` 里用引用参数会要求对象位于设备内存

## nvc++ / OpenACC 的坑

- **`present()` 遗漏不报错，只给垃圾数据**。nvc++ 对未列出的裸指针走隐式查找，通常「碰巧能跑」。
  每次给函数加数组参数，必须同步检查所有 `present()`
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
- **λ^T、M_i 等常量必须用与模拟相同的离散求和**（`./build/main --lambda`）。
  理由是**自洽性**，不是格点不准 —— 格点其实精确到 5e-6。
  但**绝不能用解析球体 `(4/3)πa³`**：真值比它大 24.1%，那是扩散界面的曲率项
  （`+8πaξ²π²/24`），误用会让 `M_i` 差 24%，直接毁掉验证

## 工作方式

- 每个阶段完成后 git 提交，commit message 用中文
- 实验日志归档到 `baseline/`
- 判据必须是**数值判据**，不是「看起来对」
- 长跑前先小规模测吞吐，别盲目启动几小时的任务
