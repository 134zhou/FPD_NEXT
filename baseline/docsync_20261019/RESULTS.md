# 文档同步复核（2026-10-19）

代码 HEAD `0b9655f`，重新构建 `fpd`/`fpd_check`/`fpd_tool`，`ACC_DEVICE_TYPE=nvidia`。

## 判据复跑

| 模式 | PASS | FAIL | 与归档对照 |
|---|---|---|---|
| `--check-config` | 33 | 0 | 与 `mass_20261017/check_config.log` 逐字节相同 |
| `--check-potential` | 14 | 0 | 逐字节相同 |
| `--check-tridiag` | 5 | 0 | 逐字节相同 |
| `--lambda 32 2.0` | 2 | 0 | 逐字节相同 |
| `--check` | 15 | 0 | 逐字节相同 |
| `--check-poisson` | 20 | 0 | 逐字节相同 |
| `--check-mass` | 11 | 0 | 逐字节相同 |
| `--check-wall` | 31 | 0 | 逐字节相同 |
| `--debug-pcg` | 2 | 0 | 与 `pcg_precision_20261019.log` 的 32³ 段相同 |

## 吞吐（`config/sed.cfg`，128×64×32，N=95，3000 步）

初态用 make_init.py 的随机构型生成到 /tmp（不改动脚本）：

```bash
python3 -c "import sys; sys.path.insert(0,'tools'); import make_init as m; \
m.NX,m.NY,m.NZ=128,64,32; m.POINTS=None; m.OUT='/tmp/fpddoc/sed_init.fpd'; m.main()"
./build/fpd config/sed.cfg --set init_file=/tmp/fpddoc/sed_init.fpd --set out_dir=/tmp/fpddoc/out \
  --set ratio_rho=<r> --set n_steps=3001 --set interval_log=1000 --set interval_ckpt=100000
```

| ratio_rho | 步/秒 | PCG 迭代/步 | 末次相对残差 |
|---|---|---|---|
| 2.0 | 35.9 | 16 | 5.9e-13 |
| 1.0 | 134.8 | — | — |
