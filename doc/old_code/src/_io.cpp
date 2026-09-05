#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <vector.h>
#include <global_definitions.h>
#include <io.h>


int getLine(char *line,FILE *pfile){
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


void getHeader3D_FPD(char *input_name, steps *time, vector3D *box, int *N)
{
	FILE *ifile=fopen(input_name, "r");
	
	if (ifile==NULL){
		printf("Error: can't open initial conditions file '%s'\n", input_name);
		exit(1);
	}
	
	int state=0;
	
	char line[MAX_LINE_LENGTH]="";
	
	// HEADER: first line
	getLine(line,ifile);

	state=sscanf(line,"%lld %lf %lf %lf %d\n", time, &box->x, &box->y, &box->z, N);
	
	if (state!=5){
		printf("Error while reading header in '%s' file\n", input_name);
		exit(1);
	}
	
	fclose(ifile);	
}



void readData3D_FPD(char *input_name, steps *time, vector3D *box, double ***vx, double ***vy, double ***vz, int *N, vector3D *R){

	FILE *ifile=fopen(input_name,"r");
	
	if (ifile==NULL){
		printf("Error: can't open initial conditions file '%s'\n", input_name);
		exit(1);
	}
	
	int state=0;
	
	char line[MAX_LINE_LENGTH]="";
	
	// HEADER: first line
	getLine(line,ifile);
	
	state=sscanf(line,"%lld %lf %lf %lf %d\n", time, &box->x, &box->y, &box->z,N);
	
	if (state!=5){
		printf("Error while reading header in '%s' file\n", input_name);
		exit(1);
	}
	
	int n;
	for(n=0;n<*N;n++){
		getLine(line,ifile);
		state=sscanf(line,"%lf %lf %lf \n", &(R[n]).x, &(R[n]).y, &(R[n]).z);
		if (state!=3){
			printf("Error while reading the configurations of particles in '%s' file\n", input_name);
			exit(1);
			}	
		}
	
	for(int i=0;i<int((*box).x);i++){
		for(int j=0;j<int((*box).y);j++){
			for(int k=0;k<int((*box).z);k++){
				getLine(line,ifile);
				state=sscanf(line,"%lf %lf %lf\n", &(vx[i][j][k]), &(vy[i][j][k]), &(vz[i][j][k]));
							
				if (state!=3){
					printf("Error while reading velocity field in '%s' file\n", input_name);
					exit(1);
					}
				}
			}
		}	
	fclose(ifile);
}

void readData3D_FPD(char *input_name, steps *time, vector3D *box, double ***vx, double ***vy, double ***vz, int *N, vector3D *R, double *radius){

	FILE *ifile=fopen(input_name, "r");
	
	if (ifile==NULL){
		printf("Error: can't open initial conditions file '%s'\n", input_name);
		exit(1);
	}
	
	int state=0;
	
	char line[MAX_LINE_LENGTH]="";
	
	// HEADER: first line
	getLine(line,ifile);
	
	state=sscanf(line,"%lld %lf %lf %lf %d\n", time, &box->x, &box->y, &box->z, N);
	
	if (state!=5){
		printf("Error while reading header in '%s' file\n", input_name);
		exit(1);
	}
	
	int n;
	for(n=0;n<*N;n++){
		getLine(line,ifile);
		
		state=sscanf(line,"%lf %lf %lf %lf\n", &(R[n]).x, &(R[n]).y, &(R[n]).z, &(radius[n]));
		if (state!=4){
			printf("Error while reading the configurations of particles in '%s' file\n", input_name);
			exit(1);
			}
		}
		
	for(int i=0;i<int((*box).x);i++){
		for(int j=0;j<int((*box).y);j++){
			for(int k=0;k<int((*box).z);k++){
				getLine(line,ifile);
				state=sscanf(line,"%lf %lf %lf\n", &(vx[i][j][k]), &(vy[i][j][k]), &(vz[i][j][k]));
							
				if (state!=3){
					printf("Error while reading velocity field in '%s' file\n", input_name);
					exit(1);
					}
				}
			}
		}	
	fclose(ifile);
}


void saveData3D_FPD(char *output_name, steps time, vector3D box, double ***vx, double ***vy, double ***vz, int N, vector3D *R){
	FILE *ofile=fopen(output_name, "w");
	if (ofile==NULL){
		printf("Error: can't open initial conditions file '%s'\n", output_name);
		exit(1);
	}
	
	// HEADER: first line
	fprintf(ofile,"%lld %.9lf %.9lf %.9lf %d \n", time, box.x, box.y, box.z, N);
	int n;
	for(int n=0;n<N;n++){
		fprintf(ofile,"%.9lf %.9lf %.9lf \n", R[n].x, R[n].y, R[n].z);
		}
	for(int i=0;i<int(box.x);i++){
		for(int j=0;j<int(box.y);j++){
			for(int k=0;k<int(box.z);k++){
				fprintf(ofile,"%.9lf %.9lf %.9lf \n", vx[i][j][k], vy[i][j][k], vz[i][j][k]);
				}
			}
		}
	fclose(ofile);
}

void saveData3D_FPD(char *output_name, steps time, vector3D box, double ***vx, double ***vy, double ***vz, int N, vector3D *R, double *radius){
	FILE *ofile=fopen(output_name, "w");
	if (ofile==NULL){
		printf("Error: can't open initial conditions file '%s'\n", output_name);
		exit(1);
	}
	
	// HEADER: first line
	fprintf(ofile,"%lld %.9lf %.9lf %.9lf %d \n", time, box.x, box.y, box.z, N);
	int n;
	for(int n=0;n<N;n++){
		fprintf(ofile,"%.9lf %.9lf %.9lf %.9lf\n", R[n].x, R[n].y, R[n].z, radius[n]);
		}
	
	for(int i=0;i<int(box.x);i++){
		for(int j=0;j<int(box.y);j++){
			for(int k=0;k<int(box.z);k++){
				fprintf(ofile,"%.9lf %.9lf %.9lf \n", vx[i][j][k], vy[i][j][k], vz[i][j][k]);
				}
			}
		}
	fclose(ofile);
}

void savePositions3D_FPD(char *output_name, steps time, vector3D box, int N, vector3D *R){
	FILE *ofile=fopen(output_name,"w");
	if (ofile==NULL){
		printf("Error: can't open initial conditions file '%s'\n", output_name);
		exit(1);
	}
	
	// HEADER: first line
	fprintf(ofile,"%lld %d %.9lf %.9lf %.9lf\n", time, N, box.x, box.y, box.z);
	int n;
	for(int n=0;n<N;n++){
		fprintf(ofile,"%.9lf %.9lf %.9lf \n", R[n].x, R[n].y, R[n].z);
		}
	
	fclose(ofile);
}

void savePositions3D_FPD(char *output_name, steps time, vector3D box, int N, vector3D *R, double *radius){
	FILE *ofile=fopen(output_name,"w");
	if (ofile==NULL){
		printf("Error: can't open initial conditions file '%s'\n", output_name);
		exit(1);
	}
	
	// HEADER: first line
	fprintf(ofile,"%lld %d %.9lf %.9lf %.9lf\n", time, N, box.x, box.y, box.z);
	int n;
	for(int n=0;n<N;n++){
		fprintf(ofile,"%.9lf %.9lf %.9lf %.9lf\n", R[n].x, R[n].y, R[n].z, radius[n]);
		}
	
	fclose(ofile);
}
