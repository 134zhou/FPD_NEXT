#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vector.h"

void ampMatrix(double **data_in,double amp,double **data_out,vector2D box);
void ampMatrix(vector2D **data_in,double amp,vector2D **data_out,vector2D box);
void ampMatrix(double ***data_in,double amp,double ***data_out,vector3D box);
void ampMatrix(vector3D ***data_in,double amp,vector3D ***data_out,vector3D box);


//3D scalar		
void onesMatrix(double ***data,vector3D box){
	int Nx = (int)(box.x);
	int Ny = (int)(box.y);
	int Nz = (int)(box.z);
	
	#pragma omp parallel for 
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			for(int k=0;k<Nz;k++){
				data[i][j][k] = 1.;
				}
			}
		}	
	}

//3D scalar		
void zerosMatrix(double ***data,vector3D box){
	int Nx = (int)(box.x);
	int Ny = (int)(box.y);
	int Nz = (int)(box.z);
	
	#pragma omp parallel for 
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			for(int k=0;k<Nz;k++){
				data[i][j][k] = 0.;
				}
			}
		}	
	}

//3D scalar		
void zerosMatrix(vector3D ***data,vector3D box){
	int Nx = (int)(box.x);
	int Ny = (int)(box.y);
	int Nz = (int)(box.z);
	
	#pragma omp parallel for 
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			for(int k=0;k<Nz;k++){
				data[i][j][k].x = 0.;
				data[i][j][k].y = 0.;
				data[i][j][k].z = 0.;
				}
			}
		}	
	}


