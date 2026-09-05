#include <stdio.h>
#include <math.h>
#include "../include/vector_mine.h"
#include "../include/cellList.h"

#define IDX3(i, j, k,Nx ,Ny, Nz) ((i) + (j)*Nx + (k)*Nx*Ny)
#define IDX2(i, j,Nx ,Ny) ((i) + (j)*Nx)

#define PRODUCT_3D(R1,R2) (((R1).x)*((R2).x) + ((R1).y)*((R2).y) + ((R1).z)*((R2).z))

	
//计算两个粒子位置的相对向量，并应用周期性边界条件（PBC）
inline void get_relativeVector_PBC(vector3D R1, vector3D R0,vector3D *d01,vector3D box)
{
	int x = (int)(box.x);
	int y = (int)(box.y);
	int z = (int)(box.z);
	d01->x = R1.x - R0.x;
	d01->y = R1.y - R0.y;
	d01->z = R1.z - R0.z;
	if(d01->x >x*0.5) d01->x -= x; else if(d01->x <-x*0.5) d01->x += x;
	if(d01->y >y*0.5) d01->y -= y; else if(d01->y <-y*0.5) d01->y += y;
	if(d01->z >z*0.5) d01->z -= z; else if(d01->z <-z*0.5) d01->z += z;
}



// 3D, morse interaction: u(r) = eps*( exp(rho*(sgm-r))*(exp(rho*(sgm-r)) - 2) )
void cul_Morse_cellList3D(double eps,double rho,double sgm,double cutoff,cellList3D *cells,int *nnIndex,vector3D *R,vector3D *F,vector3D box, int N)
{
	int N_cell = getNumberCells(cells);
	int i,j,nn,n_cell,m_cell;
	double cutoff2 = cutoff*cutoff;
	double amp = 2.*rho*eps;
	
	for(int n=0;n<N;n++)
	{
		F[n].x = 0; 
		F[n].y = 0;
		F[n].z = 0;
	}

    for(n_cell=0;n_cell<N_cell;n_cell++)
	{
		i=cells->HoC[n_cell];
		while(i!=-1)
		{
			vector3D Ri; Ri.x = R[i].x; Ri.y = R[i].y; Ri.z = R[i].z;
			vector3D Fi; Fi.x = 0; Fi.y = 0; Fi.z = 0;
			for(nn=0;nn<14;nn++)
			{
				m_cell = nnIndex[IDX2(n_cell, nn,N_cell ,14)];
				j=cells->HoC[m_cell];
				while(j!=-1)
				{
					if((m_cell==n_cell&&i<j)||m_cell!=n_cell)
					{
						/********************************************/
						//Morse 3D
						vector3D dij;
						get_relativeVector_PBC(R[j],Ri,&dij,box); 
						double Rij2 = PRODUCT_3D(dij,dij);
						if(Rij2<cutoff2)
						{
							double Rij1 = sqrt(Rij2);
							double EXP1 = exp(rho*(sgm-Rij1));
							double f_nrm = amp*(EXP1*(EXP1 - 1.))/Rij1;
							double f_x_ij = -dij.x*f_nrm;
							double f_y_ij = -dij.y*f_nrm;
							double f_z_ij = -dij.z*f_nrm;
							Fi.x +=  f_x_ij;  Fi.y +=  f_y_ij; Fi.z +=  f_z_ij;
							F[j].x += -f_x_ij;  F[j].y += -f_y_ij; F[j].z += -f_z_ij;
						}	
						/********************************************/
					}
					j=cells->LinkedList[j];	
				}
			}
			F[i].x += Fi.x; F[i].y += Fi.y; F[i].z += Fi.z;	
			i=cells->LinkedList[i];	
		}
	}
}
	

// 3D, Lennard-Jones interaction: U(r) = 4*eps*( (sgm/r)^12 -  (sgm/r)^6 ) 
void cul_LJ6_cellList3D(double sgm,double eps,double cutoff,cellList3D *cells,std::vector<int> nnIndex,vector3D *R,vector3D *F,vector3D box, int N){
	int N_cell = getNumberCells(cells);
	int i,j,nn,n_cell,m_cell;
	double cutoff2 = cutoff*cutoff;
	double sgm6 = sgm*sgm*sgm*sgm*sgm*sgm;
	
	for(int n=0;n<N;n++)
	{
		F[n].x = 0; 
		F[n].y = 0;
		F[n].z = 0;
	}
	
    for(n_cell=0;n_cell<N_cell;n_cell++)
	{
		i=cells->HoC[n_cell];
		while(i!=-1)
		{
			vector3D Ri; Ri.x = R[i].x; Ri.y = R[i].y; Ri.z = R[i].z;
			vector3D Fi; Fi.x = 0; Fi.y = 0; Fi.z = 0;
			for(nn=0;nn<14;nn++)
			{
				m_cell = nnIndex[IDX2(n_cell, nn,N_cell ,14)];
				j=cells->HoC[m_cell];
				while(j!=-1)
				{
					if((m_cell==n_cell&&i<j)||m_cell!=n_cell)
					{
						/********************************************/
						// LJ6 3D
						vector3D dij;
						get_relativeVector_PBC(R[j],Ri,&dij,box); 
						double Rij2 = PRODUCT_3D(dij,dij);
						if(Rij2<cutoff2)
						{
							double inv_Rij2 = 1./Rij2;
							double inv_rij6 = sgm6*inv_Rij2*inv_Rij2*inv_Rij2;
							double f_nrm = 24.*eps*inv_Rij2*inv_rij6*(2.*inv_rij6 - 1.);							
							double f_x_ij = -dij.x*f_nrm;
							double f_y_ij = -dij.y*f_nrm;
							double f_z_ij = -dij.z*f_nrm;
							Fi.x +=  f_x_ij;  Fi.y +=  f_y_ij; Fi.z +=  f_z_ij;
							F[j].x += -f_x_ij;  F[j].y += -f_y_ij; F[j].z += -f_z_ij;
						}
						/********************************************/
					}
					j=cells->LinkedList[j];	
				}
			}
			F[i].x += Fi.x; F[i].y += Fi.y; F[i].z += Fi.z;	
			i=cells->LinkedList[i];	
		}
	}
}		

//Morse势
void Morse_force(double De,double a,double Re,double cutoff,cellList3D *cells,std::vector<int> nnIndex,vector3D *R,vector3D *F,vector3D box, int N)
{
	int N_cell = getNumberCells(cells);
	int i,j,nn,n_cell,m_cell;
	double cutoff2 = cutoff*cutoff;

	double r_min = 1e-6;  // 设置一个非常小的最小距离
	
	for(int n=0;n<N;n++)
	{
		F[n].x = 0; 
		F[n].y = 0;
		F[n].z = 0;
	}
	
    for(n_cell=0;n_cell<N_cell;n_cell++)
	{
		i=cells->HoC[n_cell];
		while(i!=-1)
		{
			vector3D Ri; Ri.x = R[i].x; Ri.y = R[i].y; Ri.z = R[i].z;
			vector3D Fi; Fi.x = 0; Fi.y = 0; Fi.z = 0;
			for(nn=0;nn<14;nn++)
			{
				m_cell = nnIndex[IDX2(n_cell, nn,N_cell ,14)];
				j=cells->HoC[m_cell];
				while(j!=-1)
				{
					if((m_cell==n_cell&&i<j)||m_cell!=n_cell)
					{
						/********************************************/
						// LJ6 3D
						vector3D dij;
						get_relativeVector_PBC(R[j],Ri,&dij,box); 
						double Rij2 = PRODUCT_3D(dij,dij);
						
						double r;
						if (Rij2<r_min)
						{
							r = r_min;
						}
						else
						{
							r = sqrt(Rij2);
						}
						
						if(Rij2<cutoff2)
						{
							double f_nrm = -2*De*a*(1-exp(-a*(r-Re)))*exp(-a*(r-Re))/r;
							double f_x_ij = dij.x*f_nrm;
							double f_y_ij = dij.y*f_nrm;
							double f_z_ij = dij.z*f_nrm;
							Fi.x +=  f_x_ij;  Fi.y +=  f_y_ij; Fi.z +=  f_z_ij;
							F[j].x += f_x_ij;  F[j].y += f_y_ij; F[j].z += f_z_ij;
						}
						/********************************************/
					}
					j=cells->LinkedList[j];	
				}
			}
			F[i].x += -Fi.x; F[i].y += -Fi.y; F[i].z += -Fi.z;	
			i=cells->LinkedList[i];	
		}
	}
}

// 3D, WCA interaction: repulsive part of LJ interaction 
void cul_WCA_cellList3D(double sgm,double eps,cellList3D *cells,int *nnIndex,vector3D *R,vector3D *F,vector3D box, int N){
	int N_cell = getNumberCells(cells);
	int i,j,nn,n_cell,m_cell;
	double cutoff2 = 1.12246048*1.12246048*sgm*sgm;
	double sgm6 = sgm*sgm*sgm*sgm*sgm*sgm;
	
	for(int n=0;n<N;n++){
		F[n].x = 0; 
		F[n].y = 0;
		F[n].z = 0;
		}
	
    for(n_cell=0;n_cell<N_cell;n_cell++){
		i=cells->HoC[n_cell];
		while(i!=-1){
			vector3D Ri; Ri.x = R[i].x; Ri.y = R[i].y; Ri.z = R[i].z;
			vector3D Fi; Fi.x = 0; Fi.y = 0; Fi.z = 0;
			for(nn=0;nn<14;nn++){
				m_cell = nnIndex[IDX2(n_cell, nn,N_cell ,14)];
				j=cells->HoC[m_cell];
				while(j!=-1){
					if((m_cell==n_cell&&i<j)||m_cell!=n_cell){
						/********************************************/
						// WCA 3D
						vector3D dij;
						get_relativeVector_PBC(R[j],Ri,&dij,box); 
						double Rij2 = PRODUCT_3D(dij,dij);
						if(Rij2<cutoff2){
							double inv_Rij2 = 1./Rij2;
							double inv_rij6 = sgm6*inv_Rij2*inv_Rij2*inv_Rij2;
							double f_nrm = 24.*eps*inv_Rij2*inv_rij6*(2.*inv_rij6 - 1.);							
							double f_x_ij = -dij.x*f_nrm;
							double f_y_ij = -dij.y*f_nrm;
							double f_z_ij = -dij.z*f_nrm;
							Fi.x +=  f_x_ij;  Fi.y +=  f_y_ij; Fi.z +=  f_z_ij;
							F[j].x += -f_x_ij;  F[j].y += -f_y_ij; F[j].z += -f_z_ij;
							}
						/********************************************/
						}
					j=cells->LinkedList[j];	
					}
				}
			F[i].x += Fi.x; F[i].y += Fi.y; F[i].z += Fi.z;	
			i=cells->LinkedList[i];	
			}
		}
	}	



//3D
void cul_LJn_cellList3D_notail(double sgm,double eps, int n, double cutoff, cellList3D *cells,int *nnIndex,vector3D *R,vector3D *F,vector3D box,  int N){
	int N_cell = getNumberCells(cells);
	int i,j,nn,n_cell,m_cell;
	double cutoff2 = cutoff*cutoff;
	double n4eps = n*4.*eps;

	double c2 = (2.*pow(cutoff, -(2*n+2)) - pow(cutoff, -(n+2)))/(sgm*sgm);

	for(int n=0;n<N;n++){
		F[n].x = 0;
		F[n].y = 0;
		F[n].z = 0;
		}

    for(n_cell=0;n_cell<N_cell;n_cell++){
		i=cells->HoC[n_cell];
		while(i!=-1){
			vector3D Ri; Ri.x = R[i].x; Ri.y = R[i].y; Ri.z = R[i].z;
			vector3D Fi; Fi.x = 0; Fi.y = 0; Fi.z = 0;
			for(nn=0;nn<14;nn++){
				m_cell = nnIndex[IDX2(n_cell, nn,N_cell ,14)];
				j=cells->HoC[m_cell];
				while(j!=-1){
					if((m_cell==n_cell&&i<j)||m_cell!=n_cell){
						/********************************************/
						// LJ6 3D
						vector3D dij;
						get_relativeVector_PBC(R[j],Ri,&dij,box);
						double Rij2 = PRODUCT_3D(dij,dij);
						if(Rij2/(sgm*sgm)<cutoff2){
							double inv_Rij2 = 1./Rij2;
							double inv_rijn = pow(sgm*sqrt(inv_Rij2), n);
							double f_nrm = n4eps*inv_Rij2*( inv_rijn*(2.*inv_rijn-1) - c2*Rij2 );
							double f_x_ij = -dij.x*f_nrm;
							double f_y_ij = -dij.y*f_nrm;
							double f_z_ij = -dij.z*f_nrm;
							Fi.x +=  f_x_ij;  Fi.y +=  f_y_ij; Fi.z +=  f_z_ij;
							F[j].x += -f_x_ij;  F[j].y += -f_y_ij; F[j].z += -f_z_ij;
							}
						/********************************************/
						}
					j=cells->LinkedList[j];
					}
				}
			F[i].x += Fi.x; F[i].y += Fi.y; F[i].z += Fi.z;
			i=cells->LinkedList[i];
			}
		}
	}
