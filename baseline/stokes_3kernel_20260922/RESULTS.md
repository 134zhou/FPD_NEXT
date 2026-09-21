# Stokes 应力三 kernel 重构结果

环境：nvc++ 26.3、CMake 4.2.3；起点提交
`f4ae91b37cafbed1ca7baeb87458cba2d73c5b03`。

## 已完成

- 应力装配拆成体心/XY 棱、内部 YZ/ZX 棱、两片壁面三个 OpenACC kernel。
- 壁面 kernel 显式实现上下壁 `∓2ηv`，内部与壁面的噪声因子分别固定为 1 和 2。
- 删除五个 wall helper、`adv_on/noise_skip` 和 W5'a；保留 W5b。
- 新增 W3b，独立检查两片壁面、两个切向分量的通量展开式。
- `cmake --build build -j` 成功。
- `git diff --check` 通过；已删除符号的全树活动代码搜索无残留。
- `fpd_check --check`、`--check-config`、`--check-potential`、`--check-tridiag` 全部通过。

## 尚未验证

2026-09-22 执行 `nvidia-smi` 时无法连接 NVIDIA 驱动。按照仓库约定，不能把 CPU
回退当作 GPU 验证，因此以下项目待 GPU 环境恢复：

- `fpd_check --check-wall`，包括 W3、W3b、W5b
- `fpd_check --check-poisson`
- 无噪声壁面检查点修改前后逐位对照
- 短程吞吐记录

`--check-wall` 与 `--check-poisson` 均在 `src/State.cpp:131` 创建 xy cuFFT plan 时以
`cuFFT 错误 5` 中止。当前结论是“实现完成、成功编译且不依赖 cuFFT 的判据通过”，
不是“壁面 GPU 数值验证通过”。
