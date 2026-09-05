#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <vector>

#include "../include/vector_mine.h"
#include "../include/cellList.h"

#define SQR(x) ((x)*(x))

#define IDX3(i, j, k,Nx ,Ny, Nz) ((i) + (j)*Nx + (k)*Nx*Ny)
#define IDX2(i, j,Nx ,Ny) ((i) + (j)*Nx)

static int module(int n,int mo)
{
	n=n%mo;
	
	while (n<0)
	{
		n+=mo;
	}
	
	return n;
}

//3D
cellList3D* getList(vector3D box,double cutoff,int num_particles)
{
	cellList3D *l = new cellList3D[1];
	l->NumberCells_x=(int)(box.x/cutoff);
	l->NumberCells_y=(int)(box.y/cutoff);
	l->NumberCells_z=(int)(box.z/cutoff);
	
	l->CellSize_x=box.x/(double)l->NumberCells_x;
	l->CellSize_y=box.y/(double)l->NumberCells_y;
	l->CellSize_z=box.z/(double)l->NumberCells_z;
	
	l->HoC=(int*)calloc(l->NumberCells_x*l->NumberCells_y*l->NumberCells_z,sizeof(int));
	l->LinkedList=(int*)calloc(num_particles,sizeof(int));
		
	int i;
	int NNN=l->NumberCells_x*l->NumberCells_y*l->NumberCells_z;
	for (i=0;i<NNN;i++)
	{
		(l->HoC)[i]=-1;
	}
	
	return l;
}

int getNumberCells(cellList3D *l)
{
	return l->NumberCells_x*l->NumberCells_y*l->NumberCells_z;
}

void freeList(cellList3D *l)
{
	free(l->HoC);
	free(l->LinkedList);
	free(l);
}

void resetList(cellList3D *l)
{
	int i;
	int nnn=l->NumberCells_x*l->NumberCells_y*l->NumberCells_z;
	
	for (i=0;i<nnn;i++){
		(l->HoC)[i]=-1;
		}
}

void updateList(cellList3D *l,const vector3D *pos,int num)
{
	
	int NNN=l->NumberCells_x*l->NumberCells_y*l->NumberCells_z;
	
	// HoC initialization
	for (int i=0;i<NNN;i++)
		(l->HoC)[i]=-1;
		
	double inv_cellSize_x = 1./(l->CellSize_x);
	double inv_cellSize_y = 1./(l->CellSize_y);
	double inv_cellSize_z = 1./(l->CellSize_z); 
	
	int numCell_x = l->NumberCells_x;
	int numCell_y = l->NumberCells_y;
	int numCell_z = l->NumberCells_z;
	
	// colloids loop
	for (int i=0;i<num;i++)
	{
		int posx = (int)floor(pos[i].x*inv_cellSize_x);
		int posy = (int)floor(pos[i].y*inv_cellSize_y);
		int posz = (int)floor(pos[i].z*inv_cellSize_z);
		
		posx = module(posx,numCell_x);
		posy = module(posy,numCell_y);
		posz = module(posz,numCell_z);
		
		int ncell = posx+(posy+posz*numCell_y)*numCell_x;
		
		(l->LinkedList)[i] = (l->HoC)[ncell];
		(l->HoC)[ncell] = i;
	}
	
}

//number of nearest neighborhood
std::vector<int> get_nnIndex(cellList3D *l)
{
	int NNN = 14;//number of nearest neighborhood
	int N_cell = getNumberCells(l);
	int Nx = (l->NumberCells_x);
	int Ny = (l->NumberCells_y);
	int Nz = (l->NumberCells_z);
	std::vector<int> nnIndex(N_cell*NNN);
	for(int n_cell=0;n_cell<N_cell;n_cell++)
	{
		int i = n_cell%Nx;
		int j = (n_cell%(Nx*Ny)-i)/Nx;
		int k = (n_cell-i-Nx*j)/(Nx*Ny);
		nnIndex[IDX2(n_cell, 0,N_cell ,NNN)] = i + Nx*j + Nx*Ny*k;//x,y,z

		nnIndex[IDX2(n_cell, 1,N_cell ,NNN)] = (i+1+Nx)%Nx + Nx*j + Nx*Ny*k;//x+1 y z
		nnIndex[IDX2(n_cell, 2,N_cell ,NNN)] = i + Nx*((j+1+Ny)%Ny) + Nx*Ny*k;//x y+1 z 
		nnIndex[IDX2(n_cell, 3,N_cell ,NNN)] = i + Nx*j + Nx*Ny*((k+1+Nz)%Nz);//x y z+1 
		
		nnIndex[IDX2(n_cell, 4,N_cell ,NNN)] = (i+1+Nx)%Nx + Nx*((j+1+Ny)%Ny) + Nx*Ny*k;//x+1 y+1 z
		nnIndex[IDX2(n_cell, 5,N_cell ,NNN)] = i + Nx*((j+1+Ny)%Ny) + Nx*Ny*((k+1+Nz)%Nz);//x y+1 z+1
		nnIndex[IDX2(n_cell, 6,N_cell ,NNN)] = (i+1+Nx)%Nx + Nx*j + Nx*Ny*((k+1+Nz)%Nz);//x+1 y z+1
		
		nnIndex[IDX2(n_cell, 7,N_cell ,NNN)] =  (i+1+Nx)%Nx + Nx*((j+1+Ny)%Ny) + Nx*Ny*((k+1+Nz)%Nz);//x+1 y+1 z+1
		
		nnIndex[IDX2(n_cell, 8,N_cell ,NNN)] =  (i+1+Nx)%Nx + Nx*((j+1+Ny)%Ny) + Nx*Ny*((k-1+Nz)%Nz);//x+1 y+1 z-1
		nnIndex[IDX2(n_cell, 9,N_cell ,NNN)] =  (i-1+Nx)%Nx + Nx*((j+1+Ny)%Ny) + Nx*Ny*((k+1+Nz)%Nz);//x-1 y+1 z+1
		nnIndex[IDX2(n_cell, 10,N_cell ,NNN)] = (i+1+Nx)%Nx + Nx*((j-1+Ny)%Ny) + Nx*Ny*((k+1+Nz)%Nz);//x+1 y-1 z+1
		
		nnIndex[IDX2(n_cell, 11,N_cell ,NNN)] = (i+1+Nx)%Nx + Nx*j + Nx*Ny*((k-1+Nz)%Nz);//x+1 y z-1
		nnIndex[IDX2(n_cell, 13,N_cell ,NNN)] = (i-1+Nx)%Nx + Nx*((j+1+Ny)%Ny) + Nx*Ny*k;//x-1 y+1 z
		nnIndex[IDX2(n_cell, 13,N_cell ,NNN)] = i + Nx*((j-1+Ny)%Ny) + Nx*Ny*((k+1+Nz)%Nz);//x y-1 z+1
	}
	return nnIndex;	
}
