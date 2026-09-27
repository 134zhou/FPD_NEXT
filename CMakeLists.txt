cmake_minimum_required(VERSION 3.18)

# 【核心修改】直接硬编码定义编译器路径
# 注意：这必须放在 project() 命令之前
set(CMAKE_CXX_COMPILER "nvc++") 
# 如果系统不识别 nvc++，请替换为绝对路径，例如：
# set(CMAKE_CXX_COMPILER "/opt/nvidia/hpc_sdk/Linux_x86_64/25.5/compilers/bin/nvc++")

project(MyProject LANGUAGES CXX)

# 1. 设置全局编译选项
add_compile_options(-g -O3 -Wall)

# 2. 设置 OpenACC 标志
add_compile_options(-acc)
add_link_options(-acc)

# 3. 查找 CUDA 库
find_package(CUDAToolkit REQUIRED)

# 子目录定义目标，但仍把可执行文件放在 build/ 根目录，
# 保持 ./build/fpd 等既有调用路径不变。
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# 默认 all 只编译生产程序。验证和工具仍保留为可显式构建的目标：
# cmake --build build --target fpd_check -j
# cmake --build build --target fpd_tool -j
add_subdirectory(src)
add_subdirectory(tests EXCLUDE_FROM_ALL)
