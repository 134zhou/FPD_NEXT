# 自检与验证程序。Analysis.cpp 提供常量表和分块平均。
add_executable(fpd_check
    Analysis.cpp
    CheckMain.cpp
    CheckConfig.cpp
    CheckStencil.cpp
    CheckPotential.cpp
    CheckPoisson.cpp
    CheckWall.cpp
)
target_link_libraries(fpd_check PRIVATE fpd_core)

# 纯 CPU 的检查工具。--sed-stats 与 fpd_check 共用 Analysis.cpp，
# 避免在 Python 侧再定义一份分块平均和平台判据。
add_executable(fpd_tool tool_main.cpp Analysis.cpp)
target_link_libraries(fpd_tool PRIVATE fpd_core)

# spike/ 是一次性前置风险验证，不定义 CMake 目标。
