#!/usr/bin/env python3
"""历史 .fpd v3 -> VTK；参数、并行和增量转换与 fpd2vtk.py 相同。

  $PV/bin/pvpython tools/fpd2vtk_v3.py out/run_*.fpd --config out/run.config.used -o vis/v3/

只替换本进程的格式读取模块，共用现有 VTK 写出实现。
输出 velocity 和 unwrapped，不输出旧文件末尾的 force。
请使用新的输出目录，避免增量转换复用历史 VTK 产物。
"""

import fpd2vtk as converter
import fpd_format_v3

# 放在模块顶层，使 multiprocessing spawn 的 worker 同样使用 v3 布局。
# 不修改 fpd_format 模块；现有 v4 入口在独立进程中仍只接受 v4。
converter.F = fpd_format_v3


if __name__ == "__main__":
    raise SystemExit(converter.main())
