#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "../include/vector_mine.h"
#include "../global_definitions.h"
#include "../include/io.h"

#define IDX3(i, j, k,Nx ,Ny, Nz) ((i) + (j)*Nx + (k)*Nx*Ny)
#define IDX2(i, j,Nx ,Ny) ((i) + (j)*Nx)

int getLine(char *line,FILE *pfile)
{
	if (fgets(line, MAX_LINE_LENGTH, pfile) == NULL)
		return 0;
	else
		return strlen(line);
}


vector3D getBoxLength3D(char *input_name){
	FILE *ifile=fopen(input_name, "r");
	
	if (ifile==NULL){
		printf("Error: can't open initial conditions file '%s'\n", input_name);
		exit(1);
	}
	
	vector3D box;
	int state=0;
	char line[MAX_LINE_LENGTH]="";
	
	// HEADER: first line
	getLine(line,ifile);
	
	state=sscanf(line,"%*d %*d %lf %lf %lf\n", &box.x, &box.y, &box.z);
	
	if (state!=3){
		printf("Error while reading box length in '%s' file\n", input_name);
		exit(1);
	}
	
	fclose(ifile);
	
	return box;
}

/*从一个输入文件中读取模拟的初始条件
*格式为"%lld %lf %lf %lf %d\n"，分别是模拟的时间步、盒子的尺寸的xyz、粒子数量
*它读取文件的第一行，并解析其中的时间步、盒子尺寸和粒子数量
*如果文件读取或数据解析失败，函数会输出错误信息并终止程序
*/
void getHeader3D_FPD(char *input_name, steps *time, vector3D *box, int *N)
{
	FILE *ifile=fopen(input_name, "r");
	
	if (ifile==NULL)
	{
		printf("Error: can't open initial conditions file '%s'\n", input_name);
		exit(1);
	}
	
	int state=0;
	
	char line[MAX_LINE_LENGTH]="";
	
	// HEADER: first line
	getLine(line,ifile);

	state=sscanf(line,"%lld %lf %lf %lf %d\n", time, &box->x, &box->y, &box->z, N);
	
	if (state!=5)
	{
		printf("Error while reading header in '%s' file\n", input_name);
		exit(1);
	}
	
	fclose(ifile);	
}


/*从一个输入文件中读取模拟的初始条件
*格式为"%lld %lf %lf %lf %d\n"，分别是模拟的时间步、盒子的尺寸的xyz、粒子数量
*它读取文件的第一行，并解析其中的时间步、盒子尺寸和粒子数量
*如果文件读取或数据解析失败，函数会输出错误信息并终止程序
*/
void readData3D_FPD(char *input_name, steps *time, vector3D *box, vector3D *v, int *N, vector3D *R)
{

	FILE *ifile=fopen(input_name,"r");
	
	if (ifile==NULL)
	{
		printf("Error: can't open initial conditions file '%s'\n", input_name);
		exit(1);
	}
	
	int state=0;
	
	char line[MAX_LINE_LENGTH]="";
	
	// HEADER: first line
	getLine(line,ifile);

	//printf("Read line: %s\n", line);
	
	state=sscanf(line,"%lld %lf %lf %lf %d\n", time, &box->x, &box->y, &box->z,N);
	
	if (state!=5)
	{
		printf("Error while reading header in '%s' file\n", input_name);
		exit(1);
	}
	
	
	for(int n=0;n<*N;n++)
	{
		getLine(line,ifile);

		//printf("Read line: %s\n", line);

		state=sscanf(line,"%lf %lf %lf \n", &(R[n]).x, &(R[n]).y, &(R[n]).z);
		if (state!=3)
		{
			printf("Error while reading the configurations of particles in '%s' file\n", input_name);
			exit(1);
		}	
	}

	int Nx = int((*box).x);
	int Ny = int((*box).y);
	int Nz = int((*box).z);
	
	for(int i=0;i<Nx;i++)
	{
		for(int j=0;j<Ny;j++)
		{
			for(int k=0;k<Nz;k++)
			{
				getLine(line,ifile);
				state=sscanf(line,"%lf %lf %lf\n", &(v[IDX3(i, j, k,Nx ,Ny, Nz)]).x, &(v[IDX3(i, j, k,Nx ,Ny, Nz)]).y, &(v[IDX3(i, j, k,Nx ,Ny, Nz)]).z);
							
				if (state!=3)
				{
					printf("Error while reading velocity field in '%s' file\n", input_name);
					exit(1);
				}
			}
		}
	}	
	fclose(ifile);
}

void readData3D_FPD(char *input_name, steps *time, vector3D *box, vector3D *v, int *N, vector3D *R, double *radius)
{

	FILE *ifile=fopen(input_name, "r");
	
	if (ifile==NULL)
	{
		printf("Error: can't open initial conditions file '%s'\n", input_name);
		exit(1);
	}
	
	int state=0;
	
	char line[MAX_LINE_LENGTH]="";
	
	// HEADER: first line
	getLine(line,ifile);
	
	state=sscanf(line,"%lld %lf %lf %lf %d\n", time, &box->x, &box->y, &box->z, N);
	
	if (state!=5)
	{
		printf("Error while reading header in '%s' file\n", input_name);
		exit(1);
	}
	
	//int n;
	for(int n=0;n<*N;n++)
	{
		getLine(line,ifile);
		
		state=sscanf(line,"%lf %lf %lf %lf\n", &(R[n]).x, &(R[n]).y, &(R[n]).z, &(radius[n]));
		if (state!=4)
		{
			printf("Error while reading the configurations of particles in '%s' file\n", input_name);
			exit(1);
		}
	}
	
	int Nx = int((*box).x);
	int Ny = int((*box).y);
	int Nz = int((*box).z);

	for(int i=0;i<Nx;i++)
	{
		for(int j=0;j<Ny;j++)
		{
			for(int k=0;k<Nz;k++)
			{
				getLine(line,ifile);
				state=sscanf(line,"%lf %lf %lf\n", &(v[IDX3(i, j, k,Nx ,Ny, Nz)]).x, &(v[IDX3(i, j, k,Nx ,Ny, Nz)]).y, &(v[IDX3(i, j, k,Nx ,Ny, Nz)]).z);
							
				if (state!=3)
				{
					printf("Error while reading velocity field in '%s' file\n", input_name);
					exit(1);
				}
			}
		}
	}	
	fclose(ifile);
}


void saveData3D_FPD(char *output_name, steps time, vector3D box, vector3D *v, int N, vector3D *R)
{
	FILE *ofile=fopen(output_name, "w");
	if (ofile==NULL)
	{
		printf("Error: can't open initial conditions file '%s'\n", output_name);
		exit(1);
	}
	
	// HEADER: first line
	fprintf(ofile,"%lld %.9lf %.9lf %.9lf %d \n", time, box.x, box.y, box.z, N);
	
	for(int n=0;n<N;n++)
	{
		fprintf(ofile,"%.9lf %.9lf %.9lf \n", R[n].x, R[n].y, R[n].z);
	}

	int Nx = int(box.x);
	int Ny = int(box.y);
	int Nz = int(box.z);

	for(int i=0;i<Nx;i++)
	{
		for(int j=0;j<Ny;j++)
		{
			for(int k=0;k<Nz;k++)
			{
				fprintf(ofile,"%.9lf %.9lf %.9lf \n", v[IDX3(i, j, k,Nx ,Ny, Nz)].x, v[IDX3(i, j, k,Nx ,Ny, Nz)].y, v[IDX3(i, j, k,Nx ,Ny, Nz)].z);
			}
		}
	}
	fclose(ofile);
}

void saveData3D_FPD(char *output_name, steps time, vector3D box, vector3D *v, int N, vector3D *R, double *radius)
{
	FILE *ofile=fopen(output_name, "w");
	if (ofile==NULL)
	{
		printf("Error: can't open initial conditions file '%s'\n", output_name);
		exit(1);
	}
	
	// HEADER: first line
	fprintf(ofile,"%lld %.9lf %.9lf %.9lf %d \n", time, box.x, box.y, box.z, N);
	
	for(int n=0;n<N;n++)
	{
		fprintf(ofile,"%.9lf %.9lf %.9lf %.9lf\n", R[n].x, R[n].y, R[n].z, radius[n]);
	}

	int Nx = int(box.x);
	int Ny = int(box.y);
	int Nz = int(box.z);
	
	for(int i=0;i<Nx;i++)
	{
		for(int j=0;j<Ny;j++)
		{
			for(int k=0;k<Nz;k++)
			{
				fprintf(ofile,"%.9lf %.9lf %.9lf \n", v[IDX3(i, j, k,Nx ,Ny, Nz)].x, v[IDX3(i, j, k,Nx ,Ny, Nz)].y, v[IDX3(i, j, k,Nx ,Ny, Nz)].z);
			}
		}
	}
	fclose(ofile);
}

void savePositions3D_FPD(char *output_name, steps time, vector3D box, int N, vector3D *R)
{
	FILE *ofile=fopen(output_name,"w");
	if (ofile==NULL)
	{
		printf("Error: can't open initial conditions file '%s'\n", output_name);
		exit(1);
	}
	
	// HEADER: first line
	fprintf(ofile,"%lld %d %.9lf %.9lf %.9lf\n", time, N, box.x, box.y, box.z);
	
	for(int n=0;n<N;n++)
	{
		fprintf(ofile,"%.9lf %.9lf %.9lf \n", R[n].x, R[n].y, R[n].z);
	}
	
	fclose(ofile);
}

void savePositions3D_FPD(char *output_name, steps time, vector3D box, int N, vector3D *R, double *radius)
{
	FILE *ofile=fopen(output_name,"w");
	if (ofile==NULL)
	{
		printf("Error: can't open initial conditions file '%s'\n", output_name);
		exit(1);
	}
	
	// HEADER: first line
	fprintf(ofile,"%lld %d %.9lf %.9lf %.9lf\n", time, N, box.x, box.y, box.z);
	
	for(int n=0;n<N;n++)
	{
		fprintf(ofile,"%.9lf %.9lf %.9lf %.9lf\n", R[n].x, R[n].y, R[n].z, radius[n]);
	}
	
	fclose(ofile);
}


void readPositions3D_FPD(char *input_name, steps *time, vector3D *box, int *N, vector3D *R)
{

	FILE *ifile=fopen(input_name,"r");
	
	if (ifile==NULL)
	{
		printf("Error: can't open initial conditions file '%s'\n", input_name);
		exit(1);
	}
	
	int state=0;
	
	char line[MAX_LINE_LENGTH]="";
	
	// HEADER: first line
	getLine(line,ifile);
	
	state=sscanf(line,"%lld %lf %lf %lf %d\n", time, &box->x, &box->y, &box->z,N);
	
	if (state!=5)
	{
		printf("Error while reading header in '%s' file\n", input_name);
		exit(1);
	}
	
	
	for(int n=0;n<*N;n++)
	{
		getLine(line,ifile);
		state=sscanf(line,"%lf %lf %lf \n", &(R[n]).x, &(R[n]).y, &(R[n]).z);
		if (state!=3)
		{
			printf("Error while reading the configurations of particles in '%s' file\n", input_name);
			exit(1);
		}	
	}
	
	fclose(ifile);
}

