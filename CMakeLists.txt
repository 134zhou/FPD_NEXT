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

# 5. 定义源文件列表
set(SOURCES
    src/main.cpp
    src/Stokes.cpp
    src/IO.cpp
    src/Force.cpp
    src/Viscosity.cpp
    src/Velocity.cpp
    src/Analysis.cpp

)

# 6. 生成可执行文件
add_executable(main ${SOURCES})

# 7. 链接 CUDA 相关库
target_link_libraries(main 
    PRIVATE 
    CUDA::cudart 
    CUDA::curand 
    CUDA::cufft
)