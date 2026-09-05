#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
//#include <fftw3.h>
#include <mkl_dfti.h>
#include <iostream>
#include <fstream>
#include "vector.h"
#include "matrix.h"
#include "global_definitions.h"

void printStatus(long status){
	DFTI_DESCRIPTOR_HANDLE desc;
	long class_error = DftiErrorClass(status, DFTI_NO_ERROR);;
	char* error_message;
	if (! class_error) {
		std::cout << "DftiGetValue() fixes the wrong situation and returns the corresponding value" << std::endl;
		error_message = DftiErrorMessage(status);
		std::cout << "error_message = " << error_message << std::endl;
		std::cout << std::endl;
		}
	}


// intel FFT 
void poisson_PBC(double ***data_in, double ***data_out, int Nx, int Ny ,int Nz){
	//std::cout << "b3 " << std::endl;
	double _Complex in[64][64][64];
	//double _Complex out[64][64][64];
	//std::cout << "b3 " << std::endl;
	//double _Complex *** in; Matrix3D(in,Nx,Ny,Nz,double _Complex);
	//double _Complex*** out; Matrix3D(out,Nx,Ny,Nz,double _Complex);
	DFTI_DESCRIPTOR_HANDLE my_desc_handle;
	long status, l[3]; 
	l[0] = Nx; l[1] = Ny; l[2] = Nz;

	for(int j=0;j<l[0];j++){
		for(int k=0;k<l[1];k++){
			for(int s=0;s<l[2];s++){
				in[j][k][s] = data_in[j][k][s];
				}
			}
		}
			
		
	status = DftiCreateDescriptor( &my_desc_handle, DFTI_DOUBLE, DFTI_COMPLEX, 3, l);
	status = DftiCommitDescriptor( my_desc_handle);
	status = DftiComputeForward(my_desc_handle, in);
		 
	const double phase_x= 2.*M_PI/(double)Nx;
	const double phase_y= 2.*M_PI/(double)Ny;
	const double phase_z= 2.*M_PI/(double)Nz;
	
	for(int j=0;j<Nx;j++){
		const double cosXm1 = cos(phase_x*j) - 1.;
		for(int k=0;k<Ny;k++){
			const double cosYm1 = cos(phase_y*k) - 1.;
			for(int s=0;s<Nz;s++){
				const double cosZm1 = cos(phase_z*s) - 1.;
				in[j][k][s] *= 0.5/(cosXm1 + cosYm1 + cosZm1);
				}
			}
		}
	in[0][0][0] = 0;	
	
	status = DftiSetValue(my_desc_handle, DFTI_BACKWARD_SCALE, 1./(Nx*Ny*Nz));
	status = DftiCommitDescriptor( my_desc_handle);
	status = DftiComputeBackward( my_desc_handle, in);
	status = DftiFreeDescriptor( &my_desc_handle);

	for(int j=0;j<l[0];j++){
		for(int k=0;k<l[1];k++){
			for(int s=0;s<l[2];s++){
				data_out[j][k][s] = creal(in[j][k][s]);
				}
			}
		}
		
	//Free3D(in,Nx);
	//Free3D(out,Nx);
	}


/*
// real to complex
void poisson_PBC(double*** B,double*** A, int Nx, int Ny, int Nz){
	int Nz2=Nz/2+1;
	double *fftw_in;
	fftw_complex *fftw_out;
	fftw_plan plan1;
	fftw_in =(double*)malloc(sizeof(double)*Nx*Ny*2*Nz2);
	fftw_out=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Nx*Ny*Nz2);

#pragma omp parallel for
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			for(int k=0;k<Nz;k++){
				fftw_in[i*Ny*Nz+j*Nz+k]=B[i][j][k];
				}
			}
		}

	plan1=fftw_plan_dft_r2c_3d(Nx,Ny,Nz,fftw_in,fftw_out,FFTW_ESTIMATE);
	fftw_execute(plan1);
	fftw_destroy_plan(plan1);

	double dphasex=2.*M_PI/(double)Nx;
	double dphasey=2.*M_PI/(double)Ny;
	double dphasez=2.*M_PI/(double)Nz;
#pragma omp parallel for	
	for(int i=0;i<Nx;i++){
		double cosx=cos(dphasex*(double)i);
		for(int j=0;j<Ny;j++){
			double cosy=cos(dphasey*(double)j);
			for(int k=0;k<=Nz/2;k++){
				if(i==0 && j==0 && k==0){continue;}
				double cosz=cos(dphasez*(double)k);
				double cos2=cosx+cosy+cosz-3.0;
				double norm=0.5/cos2;
				int ijk = (i*Ny+j)*Nz2+k;
				fftw_out[ijk][0]*=norm;
				fftw_out[ijk][1]*=norm;
			}
		}
	}
	fftw_out[0][0]=0.;
	fftw_out[0][1]=0.;
	
	plan1=fftw_plan_dft_c2r_3d(Nx,Ny,Nz,fftw_out,fftw_in,FFTW_ESTIMATE);
	fftw_execute(plan1);
	fftw_destroy_plan(plan1);
	
	double invVol = 1./double(Nx*Ny*Nz);
#pragma omp parallel for
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			for(int k=0;k<Nz;k++){
				A[i][j][k]=invVol*fftw_in[i*Ny*Nz+j*Nz+k];
			}
		}
	}

	free(fftw_in);
	fftw_free(fftw_out);

}

*/

/*
// complex to complex
void poisson_PBC(double ***data_in, double ***data_out, int Nx, int Ny ,int Nz){
	fftw_complex *in, *out;
	fftw_plan pa, pb;

	in = (fftw_complex*)fftw_malloc( Nx * Ny * Nz * sizeof(fftw_complex) );
	out = (fftw_complex*)fftw_malloc( Nx * Ny * Nz * sizeof(fftw_complex) );
	
	
	// cpy data
	#pragma omp parallel for
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			for(int k=0;k<Nz;k++){
			const int ijk = (i*Ny+j)*Nz+k;
			in[ijk][0] = data_in[i][j][k];
			in[ijk][1] = 0;
			}
		}
	}
	
	// FFT
	pa = fftw_plan_dft_3d( Nx, Ny, Nz, in, out, FFTW_BACKWARD, FFTW_ESTIMATE );
	fftw_execute( pa );
	
	// slove poisson eq
	const double phase_x= 2.*M_PI/(double)Nx;
	const double phase_y= 2.*M_PI/(double)Ny;
	const double phase_z= 2.*M_PI/(double)Nz;
	double buf;
	#pragma omp parallel for
	for(int i=0;i<Nx;i++){
		const double cosXm1 = cos(phase_x*i) - 1.;
		for(int j=0;j<Ny;j++){
			const double cosYm1 = cos(phase_y*j) - 1.;
			for(int k=0;k<Nz;k++){
				const double cosZm1 = cos(phase_z*k) - 1.;
				const int ijk = (i*Ny+j)*Nz+k;
				const double buf = 0.5/(cosXm1 + cosYm1 + cosZm1);
				in[ijk][0] = out[ijk][0]*buf;
				in[ijk][1] = out[ijk][1]*buf;
				}
			}
		}
	in[0][0] = 0;	
	in[0][1] = 0;	
	fftw_destroy_plan( pa );
	
	
	// invert FFT
	pb = fftw_plan_dft_3d( Nx, Ny, Nz, in, out, FFTW_FORWARD, FFTW_ESTIMATE );
	fftw_execute( pb );
	const double inv_NxNyNz = 1./double(Nx*Ny*Nz);
	#pragma omp parallel for
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			for(int k=0;k<Nz;k++){
				data_out[i][j][k] = out[(i*Ny+j)*Nz+k][0]*inv_NxNyNz; 
				}
			}
		}
	fftw_destroy_plan( pb );	
	fftw_free( in );
	fftw_free( out );
}
*/
