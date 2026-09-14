cmake_minimum_required(VERSION 3.18)

# 【核心修改】直接硬编码定义编译器路径
# 注意：这必须放在 project() 命令之前
set(CMAKE_CXX_COMPILER "nvc++") 
# 如果系统不识别 nvc++，请替换为绝对路径，例如：
# set(CMAKE_CXX_COMPILER "/opt/nvidia/hpc_sdk/Linux_x86_64/25.5/compilers/bin/nvc++")

project(MyProject LANGUAGES CXX)

# 1. 设置包含路径
#    只有 src/include 需要在路径上 —— 生产头文件从这里找。
#    tests/ 下的验证代码用 "Analysis.h" 这类相对包含，编译器会先在【本文件所在
#    目录】找，所以 tests/ 的头文件不必进全局包含路径（也就不会被生产代码误引）。
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/include)

# 2. 设置全局编译选项
add_compile_options(-g -O3 -Wall)

# 3. 设置 OpenACC 标志
add_compile_options(-acc)
add_link_options(-acc)

# 4. 查找 CUDA 库
find_package(CUDAToolkit REQUIRED)

# 5. 核心静态库
#    验证程序与生产程序【必须走同一条代码路径】，否则验证结果不能外推到生产。
#    另外 OpenACC 的 TU 编译很慢，拆库避免每个可执行都重编一遍核函数。
#
#    ⚠️ fpd_core 里【只有生产代码】。验证代码（tests/）不进来 —— 否则
#       「验证与生产走同一条路径」这个保证就反过来了：生产会被迫链接测试代码。
add_library(fpd_core STATIC
    src/Stokes.cpp
    src/Force.cpp
    src/Viscosity.cpp
    src/Velocity.cpp
    src/Config.cpp
    src/IOBin.cpp
    src/State.cpp
    src/Potential.cpp
    src/Poisson.cpp
)
target_link_libraries(fpd_core PUBLIC CUDA::cudart CUDA::curand CUDA::cufft)

# 6. 模拟主程序（生产）
add_executable(fpd src/main.cpp)
target_link_libraries(fpd PRIVATE fpd_core)

# ===========================================================================
# tests/ —— 全部验证与测试代码
# ===========================================================================

# 6b. 自检与验证程序 fpd_check
#     Analysis.cpp 只被 fpd_check 用（常量表 + 分块平均），所以它属于这一侧。
add_executable(fpd_check
    tests/Analysis.cpp
    tests/CheckMain.cpp
    tests/CheckStencil.cpp
    tests/CheckPotential.cpp
    tests/CheckPoisson.cpp
    tests/CheckWall.cpp
)
target_link_libraries(fpd_check PRIVATE fpd_core)

# 7. 纯 CPU 的检查工具 fpd_tool（不需要 GPU，可在无卡机器上跑）
add_executable(fpd_tool tests/tool_main.cpp)
target_link_libraries(fpd_tool PRIVATE fpd_core)

# tests/spike/ 不参与构建：一次性的前置风险验证程序，各自用 nvc++/g++ 单独编译，
# 用法见 spike 源码头部与 README「已完成的验证」。