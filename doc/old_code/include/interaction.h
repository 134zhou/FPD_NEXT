#ifndef INTERACTION_H
#define INTERACTION_H

void cul_LJ6_cellList3D(double sgm,double eps,double cutoff,cellList3D *cells,std::vector<int> nnIndex,vector3D *R,vector3D *F,vector3D box, int N);

void cul_WCA_cellList3D(double sgm,double eps,cellList3D *cells,int *nnIndex,vector3D *R,vector3D *F,vector3D box, int N);

void cul_Morse_cellList3D(double eps,double rho,double sgm,double cutoff,cellList3D *cells,int *nnIndex,vector3D *R,vector3D *F,vector3D box, int N);

void cul_LJn_cellList3D_notail(double sgm,double eps, int n, double cutoff, cellList3D *cells,int *nnIndex,vector3D *R,vector3D *F,vector3D box,  int N);

void Morse_force(double De,double a,double Re,double cutoff,cellList3D *cells,std::vector<int> nnIndex,vector3D *R,vector3D *F,vector3D box, int N);

#endif
