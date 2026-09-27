# 验证程序与生产程序必须走同一条代码路径，所以共用 fpd_core。
# fpd_core 只包含生产代码，不得把 tests/ 中的验证实现加进来。
add_library(fpd_core STATIC
    Stokes.cpp
    Force.cpp
    Viscosity.cpp
    Velocity.cpp
    Config.cpp
    IOBin.cpp
    State.cpp
    Potential.cpp
    Poisson.cpp
)
target_include_directories(fpd_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(fpd_core PUBLIC CUDA::cudart CUDA::curand CUDA::cufft)

add_executable(fpd main.cpp)
target_link_libraries(fpd PRIVATE fpd_core)
