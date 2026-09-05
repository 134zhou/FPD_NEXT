#ifndef CELLLIST_H
#define CELLLIST_H
#include <vector>
struct cellList2D{
	int *HoC;                            // Head of Chain for linked list
	int *LinkedList;                     // linked list
	int NumberCells_x;                     // number of cells in one direction
	int NumberCells_y;
	double CellSize_x;                     // cell edge size
	double CellSize_y;
};

struct cellList3D
{
	int *HoC;                            // Head of Chain for linked list
	int *LinkedList;                     // linked list
	int NumberCells_x;                     // number of cells in one direction
	int NumberCells_y;
	int NumberCells_z;
	double CellSize_x;                     // cell edge size
	double CellSize_y;
	double CellSize_z;
};

cellList3D* getList(vector3D box,double cutoff,int num_particles);

int getNumberCells(cellList2D *l);
int getNumberCells(cellList3D *l);

std::vector<int> get_nnIndex(cellList3D *l);

void freeList(cellList3D *l);

void resetList(cellList3D *l);

void updateList(cellList3D *l,const vector3D *pos,int num);

#endif

