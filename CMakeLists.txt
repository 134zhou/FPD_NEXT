cmake_minimum_required(VERSION 3.18)

# 【核心修改】直接硬编码定义编译器路径
# 注意：这必须放在 project() 命令之前
set(CMAKE_CXX_COMPILER "nvc++") 
# 如果系统不识别 nvc++，请替换为绝对路径，例如：
# set(CMAKE_CXX_COMPILER "/opt/nvidia/hpc_sdk/Linux_x86_64/25.5/compilers/bin/nvc++")

project(MyProject LANGUAGES CXX)

# 1. 设置包含路径
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
add_library(fpd_core STATIC
    src/Stokes.cpp
    src/Force.cpp
    src/Viscosity.cpp
    src/Velocity.cpp
    src/Analysis.cpp
    src/Config.cpp
    src/IOBin.cpp
    src/State.cpp
    src/Potential.cpp
    src/Poisson.cpp
)
target_link_libraries(fpd_core PUBLIC CUDA::cudart CUDA::curand CUDA::cufft)

# 6. 模拟主程序
add_executable(fpd src/main.cpp)
target_link_libraries(fpd PRIVATE fpd_core)

# 6b. 自检与验证程序（从 fpd 拆出；与 fpd 共用 fpd_core，走同一条物理代码路径）
add_executable(fpd_check
    src/CheckMain.cpp
    src/CheckStencil.cpp
    src/CheckNoise.cpp
    src/CheckPotential.cpp
    src/CheckPoisson.cpp
)
target_link_libraries(fpd_check PRIVATE fpd_core)

# 7. 纯 CPU 的检查工具（不需要 GPU，可在无卡机器上跑）
add_executable(fpd_tool src/tool_main.cpp)
target_link_libraries(fpd_tool PRIVATE fpd_core)