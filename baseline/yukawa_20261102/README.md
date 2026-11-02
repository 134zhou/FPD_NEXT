# Yukawa / GPU 排序分箱验证

构建：`cmake -S . -B build && cmake --build build --target fpd fpd_check -j4`。
GPU 检查均设置 `ACC_DEVICE_TYPE=nvidia`。

- `neighbor.log`：`fpd_check --check-neighbor`，GPU cell/direct 及 CPU 半矩阵参考、重复建箱逐位比较。
- `bench.log`：`fpd_check --bench-neighbor`，固定近似数密度 1/27，Yukawa sigma=2、kappa=0.2、rcut=8、energy 移位；每路径预热，三批各五次，取批平均的中位数。cell 包含每步 CUDA 排序、建箱、同步和力计算。只有力阶段，不是完整模拟步。
- `config.log`、`potential.log`：配置和势函数测试。
- `check.log`、`wall.log`、`pcg_noise.log`：原有 GPU 回归。
- `restart.py` / `restart.log`：32³、N=1/2、Yukawa direct/cell，6 步连续结果与第 3 步重启结果逐位比较。运行脚本可重新生成初态与存档。

未覆盖的计划项：邻居集合逐对独立枚举、所有退化箱数/随机聚集/截断边界专项、旧版本检查点比较、整步大体系计时及不同截断的性能扫描。当前结果不能替代这些验收。
