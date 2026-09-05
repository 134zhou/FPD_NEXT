#include "./include/IO.h"

// 辅助函数：保存为 VTK 格式供 ParaView 可视化
void save_fluid_vtk(std::string path, int step, const NS_Config& cfg, const std::vector<double>& vx, const std::vector<double>& vy, const std::vector<double>& vz, const std::vector<double>& p)
{
    std::string filename = path + "/result_" + std::to_string(step) + ".vtk";
    std::ofstream ofs(filename);
    int Nx = cfg.Nx; int Ny = cfg.Ny; int Nz = cfg.Nz;

    ofs << "# vtk DataFile Version 3.0\n";
    ofs << "Stochastic Stokes Flow\n";
    ofs << "ASCII\n";
    ofs << "DATASET STRUCTURED_POINTS\n";
    ofs << "DIMENSIONS " << Nx << " " << Ny << " " << Nz << "\n";
    ofs << "ORIGIN 0 0 0\n";
    ofs << "SPACING 1 1 1\n";
    ofs << "POINT_DATA " << Nx * Ny * Nz << "\n";

    // 保存压力场 (标量)
    ofs << "SCALARS pressure double 1\nLOOKUP_TABLE default\n";
    for (const auto& val : p) ofs << val << "\n";

    // 保存速度场 (矢量)
    ofs << "VECTORS velocity double\n";
    for (int i = 0; i < Nx * Ny * Nz; ++i)
    {
        ofs << vx[i] << " " << vy[i] << " " << vz[i] << "\n";
    }
    ofs.close();
    std::cout << "Saved: " << filename << std::endl;
};

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

void save_particles_vtk
(
    std::string path, int step, int N, 
    const std::vector<double>& Rx, const std::vector<double>& Ry, const std::vector<double>& Rz,
    const std::vector<double>& Vx, const std::vector<double>& Vy, const std::vector<double>& Vz,
    const std::vector<double>& Fx, const std::vector<double>& Fy, const std::vector<double>& Fz
)
{
    std::string filename = path + "/particles_" + std::to_string(step) + ".vtk";
    std::ofstream ofs(filename);

    if (!ofs.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    // VTK 文件头
    ofs << "# vtk DataFile Version 3.0\n";
    ofs << "FPD Particle Data\n";
    ofs << "ASCII\n";
    ofs << "DATASET POLYDATA\n";

    // 1. 写入几何点坐标
    ofs << "POINTS " << N << " double\n";
    for (int i = 0; i < N; ++i)
    {
        ofs << Rx[i] << " " << Ry[i] << " " << Rz[i] << "\n";
    }

    // 2. 写入拓扑结构 (VERTICES)
    // VTK 要求：每个顶点前面要带一个数字表示该单元包含的点数（这里是 1）
    // 第二个参数是 总数据量 = N (个数) + N (每个点的索引)
    ofs << "VERTICES " << N << " " << 2 * N << "\n";
    for (int i = 0; i < N; ++i)
    {
        ofs << "1 " << i << "\n";
    }

    // 3. 写入点属性数据
    ofs << "POINT_DATA " << N << "\n";

    // 保存粒子速度
    ofs << "VECTORS velocity double\n";
    for (int i = 0; i < N; ++i)
    {
        ofs << Vx[i] << " " << Vy[i] << " " << Vz[i] << "\n";
    }

    // 保存粒子受力 (方便分析碰撞或外力)
    ofs << "VECTORS force double\n";
    for (int i = 0; i < N; ++i)
    {
        ofs << Fx[i] << " " << Fy[i] << " " << Fz[i] << "\n";
    }

    ofs.close();
    // 使用 \r 打印可以保持终端整洁，只在同一行更新
    std::cout << "Saved: " << filename << std::endl;
}