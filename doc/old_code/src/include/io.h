#ifndef INPUT_OUTPUT_H
#define INPUT_OUTPUT_H

#define MAX_LINE_LENGTH 2000

vector3D getBoxLength3D(char *input_name);

void getHeader3D_FPD(char *input_name, steps *time, vector3D *box, int *N);

void readData3D_FPD(char *input_name, steps *time, vector3D *box, vector3D ***v, int *N, vector3D *R);
void readData3D_FPD(char *input_name, steps *time, vector3D *box, vector3D ***v, int *N, vector3D *R, double *radius);

void saveData3D_FPD(char *output_name, steps time, vector3D box, vector3D ***v, int N, vector3D *R);
void saveData3D_FPD(char *output_name, steps time, vector3D box, vector3D ***v, int N, vector3D *R, double *radius);

void savePositions3D_FPD(char *output_name, steps time, vector3D box, int N, vector3D *R);
void savePositions3D_FPD(char *output_name, steps time, vector3D box, int N, vector3D *R, double *radius);

void readPositions3D_FPD(char *input_name, steps *time, vector3D *box, int *N, vector3D *R);

int getLine(char *line, FILE *pfile);

#endif

