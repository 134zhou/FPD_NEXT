# CMake 按目标构建拆分

## 基线

- 根 `CMakeLists.txt` 同时定义 `fpd`、`fpd_check`、`fpd_tool`。
- 三个可执行目标都属于默认 `all`，所以 `cmake --build build -j` 会全部编译。
- 三个目标共用仅含生产代码的 `fpd_core`。

## 实施

1. 根 CMake 保留编译器、OpenACC、CUDA 和输出目录设置。
2. `src/CMakeLists.txt` 定义 `fpd_core` 和 `fpd`。
3. `tests/CMakeLists.txt` 定义 `fpd_check` 和 `fpd_tool`，并通过
   `add_subdirectory(tests EXCLUDE_FROM_ALL)` 排除出默认构建。
4. 保持可执行文件位于 `build/` 根目录，不改调用接口。

## 验收标准

- 全新配置后，默认构建只编译并生成 `fpd`。
- `--target fpd_check` 和 `--target fpd_tool` 能分别独立构建。
- 三个可执行文件仍位于 `build/` 根目录，并且都链接同一个 `fpd_core`。
- 源文件集、编译选项和物理实现不变；本阶段不需要数值轨迹对照。
