# make_init.py 精简结果

基线提交：`5609a3d`。旧命令与新脚本均生成 128×64×32、N=95、step=0 的 v3 文件。

| 判据 | 结果 |
|---|---|
| Python 语法 | `make_init.py`、`fpd_format.py` 通过 `py_compile` |
| 默认沉降构型 | 前后 SHA-256 均为 `c8da174a34bfde187a53b7bda3a70e8b476398384e52590c4417a6efc9b9f6f0` |
| 数组对照 | `fpd_tool --diff-ckpt` 的 15 个数组全部逐位相同 |
| 轴序互操作 | 8×4×2 空体系的 `(3,1,1)` 三分量均为 `3001001` |
| v3 完整性 | 两个文件均由 `fpd_tool` 确认校验和通过 |

脚本 248 行减至 85 行；未修改 `fpd_format.py` 的格式实现或 C++ 数值代码，因此不需要 GPU 判据。
