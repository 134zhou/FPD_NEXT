#ifndef IO_H
#define IO_H
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <openacc.h>
#include "Common.h"

void save_fluid_vtk(std::string path, int step, const NS_Config& cfg, const std::vector<double>& vx, const std::vector<double>& vy, const std::vector<double>& vz, const std::vector<double>& p);

void save_particles_vtk
(
    std::string path, int step, int N, 
    const std::vector<double>& Rx, const std::vector<double>& Ry, const std::vector<double>& Rz,
    const std::vector<double>& Vx, const std::vector<double>& Vy, const std::vector<double>& Vz,
    const std::vector<double>& Fx, const std::vector<double>& Fy, const std::vector<double>& Fz
);
#endif