# Stokes 应力三 kernel 重构计划

起点提交：`f4ae91b37cafbed1ca7baeb87458cba2d73c5b03`。

目标是让实现直接对应三个离散定义域，而不是追求性能：

1. 体心与 XY 棱：`k = 0 .. Nz-1`；
2. 内部 YZ/ZX 棱：`k = 0 .. Nz-2`；
3. 两片壁面：`k = -1, Nz-1`，显式写代入 ghost 后的 `±2ηv`。

同步删除 `wz_vx/wz_vy/wz_vz`、错误的壁面判定 helper、测试专用
`adv_on/noise_skip` 和依赖它们的 W5'a。保留 W5b 的 `noise_gamma1` 对照，新增
W3b 直接核对上下壁 `Π_yz/Π_zx`。

验收以构建、W3/W3b/W5b、完整 `--check-wall`、`--check`、`--check-poisson` 和
无噪声检查点逐位对照为准；不设置性能提升门槛。
