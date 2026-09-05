#include <stdlib.h>
#include <stdio.h>
#include <fftw3.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include "vector.h"
#include "matrix.h"
#include "global_definitions.h"


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
