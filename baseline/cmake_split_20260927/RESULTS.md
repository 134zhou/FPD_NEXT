# CMake 按目标构建验收结果

日期：2026-09-27

## 干净构建

在全新的 `/tmp/fpd-cmake-split-20260927-make` 目录执行：

```bash
cmake -S . -B /tmp/fpd-cmake-split-20260927-make
cmake --build /tmp/fpd-cmake-split-20260927-make -j
```

结果：PASS。构建日志只包含 `fpd_core` 和 `fpd`；构建根目录只有
`fpd` 这一个可执行文件，没有编译 `fpd_check` 或 `fpd_tool`。

## 按需目标

```bash
cmake --build /tmp/fpd-cmake-split-20260927-make --target fpd_check -j
cmake --build /tmp/fpd-cmake-split-20260927-make --target fpd_tool -j
```

结果：PASS。两个目标分别构建成功，都复用已构建的 `fpd_core`；
`fpd`、`fpd_check`、`fpd_tool` 均位于构建根目录。

## 纯 CPU 自检

```bash
/tmp/fpd-cmake-split-20260927-make/fpd_check --check-config
/tmp/fpd-cmake-split-20260927-make/fpd_check --check-potential
/tmp/fpd-cmake-split-20260927-make/fpd_check --check-tridiag
```

结果：三项全部 PASS。本次只改构建组织，源文件集、链接库和数值内核
均未改动，因此未重跑 GPU 数值判据或生产模拟。
