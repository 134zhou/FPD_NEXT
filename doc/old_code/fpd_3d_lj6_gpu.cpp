#include <stdio.h>
#include <math.h>
#include <ctime>
//#include <omp.h>

#include <iostream>
#include <fstream>
/*
#include <fftw3.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
*/

#include "global_definitions.h"
#include "vector.h"
#include "matrix.h"
#include "io.h"
//#include "poisson.h"
#include "cellList.h"
#include "interaction.h"

#include <openacc.h>
#include <curand.h>
#include <cufft.h>
//#include "cufft/cufft_c2c.cuh"
//#include "curand/fpd_curand.cuh"

// function to determine the shape of the fluid particles
static inline double order(double dx, double dy, double dz, double radi, double inv_xi)
{
	return 0.5*(tanh((radi - sqrt(dx*dx + dy*dy + dz*dz))*inv_xi) + 1.);
}

// fft
void cufft_c2c_forward(double *d_data, int Nx, int Ny, int Nz, void *stream);
void cufft_c2c_inverse(double *d_data, int Nx, int Ny, int Nz, void *stream);
// randam number generator
curandGenerator_t generator;
void fpd_curand_init(unsigned long int seed);
void fpd_curand_fin(void);
void fpd_curand_normal(float *rand, int n, float mean, float stddev);


int main(int argc,char* argv[])
{
	std::time_t timer; std::time(&timer);
	std::cout.precision(15);
		
		
	if(argc!=6)
	{
		std::cout << argv[0] <<" [init state] [output directory] [time range] [beta] [eps lj6]" << std::endl;
		exit(0);
	}	
	else
	{
		// output setting parametersge
		std::cout << "execute time: " << std::ctime(&timer) << std::endl;
		std::cout << "init state: " << argv[1] << std::endl;
		std::cout << "output directory: " << argv[2] << std::endl;
		std::cout << "time range: " << argv[3] << std::endl; 
		std::cout << "beta=" << argv[4] <<", eps(LJ6)="<< argv[5] << std::endl;
	}
	
	// setting parameters-------------------
	int endTime = atoi(argv[3]);//模拟的结束时间步数
	const double inv_beta = 1./atof(argv[4]);//beta的倒数，beta=1/(k_B*T)
	
	const double inv_dt = 1./DT;//DT在global_definitions.h中定义，这是它的倒数
	const double inv_xi = 1./XI;//XT在global_definitions.h中定义，这是它的倒数

	const double eps = atof(argv[5]);//Lennard-Jones 势能的强度参数
	const double sgm = 2.*RADIUS + 1./inv_xi;//Lennard-Jones 势能的截距距离
	
	const int interval_pos = INTERVAL_POSITIONS*(int)(inv_dt);//每隔多少时间步保存一次粒子的位置数据？
	const int interval_vel = INTERVAL_VELOSITY*(int)(inv_dt);//决定每隔多少时间步保存一次流体的速度数据？
	
	const double ratio_eta_m1 = RATIO_ETA - 1.;//控制流体的粘度变化？
	const double W = sqrt(2.*inv_beta*inv_dt);//控制随机热噪声的强度，用于模拟热波动？
	
	// memory allocation-------------------
	steps time;//模拟的时间步
	vector3D box;//模拟区域的大小
	int N;//粒子数量
	getHeader3D_FPD(argv[1],&time,&box,&N);
	int Nx = (int)(box.x);//= (int)(box.x)
	int Ny = (int)(box.y);//= (int)(box.y)
	int Nz = (int)(box.z);//= (int)(box.z)
	double inv_N = 1./(double)N;
	double inv_boxSize = 1./(double)(Nx*Ny*Nz);
	
	// Rx, Ry, Rz: 分别存储每个粒子的x、y、z坐标。
	// Vx, Vy, Vz: 分别存储每个粒子的x、y、z方向速度。
	// Fx, Fy, Fz: 分别存储每个粒子的x、y、z方向受力。

	double *Rx = new double[N];//每个粒子的x坐标
	double *Ry = new double[N];//每个粒子的y坐标
	double *Rz = new double[N];//每个粒子的z坐标
	double *Vx = new double[N];//每个粒子的x方向速度
	double *Vy = new double[N];//每个粒子的y方向速度
	double *Vz = new double[N];//每个粒子的z方向速度
	double *Fx = new double[N];//每个粒子的x方向受力
	double *Fy = new double[N];//每个粒子的y方向受力
	double *Fz = new double[N];//每个粒子的z方向受力
	
	// fx, fy, fz: 分别存储流体在x、y、z方向的受力场。

	double ***fx; Matrix3D(fx,Nx,Ny,Nz,double);//储流体在x方向的受力场
	double ***fy; Matrix3D(fy,Nx,Ny,Nz,double);//储流体在y方向的受力场
	double ***fz; Matrix3D(fz,Nx,Ny,Nz,double);//储流体在z方向的受力场

	// eta, etaXY, etaYZ, etaZX: 分别存储流体的粘度场及其各方向分量。

	double ***eta; Matrix3D(eta,Nx,Ny,Nz,double);//流体的粘度场
	double ***etaXY; Matrix3D(etaXY,Nx,Ny,Nz,double);//流体的粘度场分量
	double ***etaYZ; Matrix3D(etaYZ,Nx,Ny,Nz,double);//流体的粘度场分量
	double ***etaZX; Matrix3D(etaZX,Nx,Ny,Nz,double);//流体的粘度场分量
	
	// vx, vy, vz: 分别存储流体在x、y、z方向的速度场。
	
	double ***vx; Matrix3D(vx,Nx,Ny,Nz,double);//流体在x方向的速度场
	double ***vy; Matrix3D(vy,Nx,Ny,Nz,double);//流体在y方向的速度场
	double ***vz; Matrix3D(vz,Nx,Ny,Nz,double);//流体在z方向的速度场

	double ***p; Matrix3D(p,Nx,Ny,Nz,double);//p流体的压力场。
	
	//分别存储流体的应力张量场

	double ***restrict pi_dx; Matrix3D(pi_dx,Nx,Ny,Nz,double);
	double ***restrict pi_dy; Matrix3D(pi_dy,Nx,Ny,Nz,double);
	double ***restrict pi_dz; Matrix3D(pi_dz,Nx,Ny,Nz,double);
	double ***pi_nx; Matrix3D(pi_nx,Nx,Ny,Nz,double);
	double ***pi_ny; Matrix3D(pi_ny,Nx,Ny,Nz,double);
	double ***pi_nz; Matrix3D(pi_nz,Nx,Ny,Nz,double);
	
	//分别存储流体的应力张量散度场
	
	double ***div_pix; Matrix3D(div_pix,Nx,Ny,Nz,double);
	double ***div_piy; Matrix3D(div_piy,Nx,Ny,Nz,double);
	double ***div_piz; Matrix3D(div_piz,Nx,Ny,Nz,double);
	
	//用于存储FFT变换的输入数据，大小为2*Nx*Ny*Nz，因为FFT变换需要复数输入，每个复数由两个连续的double表示（实部和虚部）。
	double *fftinput = new double[2*Nx*Ny*Nz];
	
#ifdef _NOISE_
	unsigned long int seed = std::time(NULL);
	std::cout << "NOISE ON, SEED: " << seed << std::endl;
	fpd_curand_init(seed);
	const int randSize = 3*Nx*Ny*Nz;
	float *randD = new float[randSize];
	float *randN = new float[randSize];
#else
	std::cout << "NOISE OFF" << std::endl;
#endif
	
	vector3D ***v; Matrix3D(v,Nx,Ny,Nz,vector3D);//储流体的速度场
	vector3D *R = new vector3D[N];//粒子的位置
	vector3D *V = new vector3D[N];//粒子的速度
	vector3D *F = new vector3D[N];//粒子的受力
	

	// memory allocation: phase field
	const int range = int(2.*(RADIUS+1./inv_xi));//相场函数的作用范围
	const int range_m1 = range-1; 
#ifndef _8BLOCKLOOP_
	const double range2 = range*range;//range的平方
#else
	const double range2 = (double)(RANGE_CORE*RANGE_CORE);
	std::cout << "8BLOCKLOOP ON, RANGE_CORE=" << RANGE_CORE << std::endl; 
#endif
	const int N_range = 2*range;//range的两倍

	double *phiX = new double[N*N_range*N_range*N_range];
	double *phiY = new double[N*N_range*N_range*N_range];
	double *phiZ = new double[N*N_range*N_range*N_range];

	double *sum_phix = new double[N];
	double *sum_phiy = new double[N];
	double *sum_phiz = new double[N];

	// read initial data
	readData3D_FPD(argv[1],&time,&box,v,&N,R);
	for(int n=0;n<N;n++)
	{
		Rx[n] = R[n].x;
		Ry[n] = R[n].y;
		Rz[n] = R[n].z;
		printf("%lf, %lf, %lf\n", Rx[n], Ry[n], Rz[n]);
	}

	for(int i=0;i<Nx;i++)
	{
		for(int j=0;j<Ny;j++)
		{
			for(int k=0;k<Nz;k++)
			{
				vx[i][j][k] = v[i][j][k].x;
				vy[i][j][k] = v[i][j][k].y;
				vz[i][j][k] = v[i][j][k].z;
			}
		}
	}
	

	// cell list
	double cutoff = CUTOFF_POTENTIAL*sgm; 
    cellList3D *cells = getList(box,CELL_SIZE,N);
    int **nnIndex = get_nnIndex(cells);
	//std::cout << "DeBUG8" << std::endl;
    updateList(cells,R,N);   
	
	
	/*
	// gaussian noise
	gsl_rng_env_setup();
	gsl_rng_type *T = (gsl_rng_type *)gsl_rng_default;
	gsl_rng **r = new gsl_rng*[NUMBER_THREADS];
	for(int i=0;i<NUMBER_THREADS;i++){
		r[i] = gsl_rng_alloc(T); 
		gsl_rng_set (r[i], (unsigned long int) std::time(NULL)+i);
		}
	*/
	
	// output setting parameters
	std::cout << "a=" << RADIUS << " xi=" << XI << " ratio_eta=" << RATIO_ETA << std::endl;
	std::cout << "box (" << box.x << " " << box.y << " " << box.z << ") " << std::endl;
	std::cout << "numParitlces=" << N << std::endl; 
	std::cout << "-------------------------------------"<< std::endl; 
	
	
	#pragma acc data create(Vx[:N],Vy[:N],Vz[:N])
	#pragma acc data create(Fx[:N],Fy[:N],Fz[:N])
	#pragma acc data copyin(Rx[:N],Ry[:N],Rz[:N])
	#pragma acc data copyin(vx[:Nx][:Ny][:Nz],vy[:Nx][:Ny][:Nz],vz[:Nx][:Ny][:Nz])
	#pragma acc data create(pi_dx[:Nx][:Ny][:Nz],pi_dy[:Nx][:Ny][:Nz],pi_dz[:Nx][:Ny][:Nz])
	#pragma acc data create(pi_nx[:Nx][:Ny][:Nz],pi_ny[:Nx][:Ny][:Nz],pi_nz[:Nx][:Ny][:Nz])
	#pragma acc data create(div_pix[:Nx][:Ny][:Nz],div_piy[:Nx][:Ny][:Nz],div_piz[:Nx][:Ny][:Nz])
	#pragma acc data create(fx[:Nx][:Ny][:Nz],fy[:Nx][:Ny][:Nz],fz[:Nx][:Ny][:Nz])
	#pragma acc data create(eta[:Nx][:Ny][:Nz],etaXY[:Nx][:Ny][:Nz],etaYZ[:Nx][:Ny][:Nz],etaZX[:Nx][:Ny][:Nz])
	#pragma acc data create(p[:Nx][:Ny][:Nz])
	#pragma acc data create(phiX[:N*N_range*N_range*N_range],phiY[:N*N_range*N_range*N_range],phiZ[:N*N_range*N_range*N_range])
	#pragma acc data create(sum_phix[:N],sum_phiy[:N],sum_phiz[:N])
	#pragma acc data create(fftinput[:2*Nx*Ny*Nz])
#ifdef _NOISE_
	#pragma acc data create(randD[:randSize],randN[:randSize])
#endif
	{// acc data
	
	// time loop
	for(int t=0;t<(int)(endTime*inv_dt);t++)
	{
	//for(int t=0;t<(int)(endTime);t++){
		//std::cout << t << std::endl;
		printf("%d ",t);
		// save 
		if(t%interval_pos==0)
		{
			#pragma acc update host(Rx[:N],Ry[:N],Rz[:N])
			for(int n=0;n<N;n++)
			{
				R[n].x = Rx[n];
				R[n].y = Ry[n];
				R[n].z = Rz[n];
			}
			char filename[200];
			sprintf(filename,"%s/pos_%05d",argv[2],(int)(DT*t)+(int)time);
			savePositions3D_FPD(filename,(int)(DT*t+time),box,N,R);
		}
		
		if(t%interval_vel==0)
		{
			#pragma acc update host(Rx[:N],Ry[:N],Rz[:N])
			for(int n=0;n<N;n++)
			{
				R[n].x = Rx[n];
				R[n].y = Ry[n];
				R[n].z = Rz[n];
			}
			#pragma acc update host(vx[:Nx][:Ny][:Nz],vy[:Nx][:Ny][:Nz],vz[:Nx][:Ny][:Nz])
			for(int i=0;i<Nx;i++)
			{
				for(int j=0;j<Ny;j++)
				{
					for(int k=0;k<Nz;k++)
					{
						v[i][j][k].x = vx[i][j][k];
						v[i][j][k].y = vy[i][j][k];
						v[i][j][k].z = vz[i][j][k];
					}
				}
			}
			char filename[200];
			sprintf(filename,"%s/velocity",argv[2]);
			//sprintf(filename,"%s/vel_%05d",argv[2],(int)(DT*t)+(int)time);
			saveData3D_FPD(filename,(int)(DT*t+time),box,v,N,R);
		}
		
		
// -- a. calculation of phase fields --
#pragma acc parallel loop collapse(2)
		for(int n=0;n<N;n++)
		{
			for(int ijk=0; ijk<N_range*N_range*N_range; ijk++)
			{
				int in = (int)Rx[n]; int jn = (int)Ry[n]; int kn = (int)Rz[n];
				int i = int(ijk/(N_range*N_range));//层同步信号
				int j = int(ijk%(N_range*N_range)/N_range);//列同步信号
				int k = int(ijk%N_range);//行同步信号
				
				//下面的作用是把0到2range-1映射到-range-0.5到range+0.5

				int ir = i+in-range_m1;
				double dx = ir - Rx[n]; double dxh = dx + 0.5;
				int jr = j+jn-range_m1;
				double dy = jr - Ry[n]; double dyh = dy + 0.5;
				int kr = k+kn-range_m1; 
				double dz = kr - Rz[n]; double dzh = dz + 0.5;
				
				int idx = n*(N_range*N_range*N_range)+i*(N_range*N_range)+j*N_range+k;
				if(range2<dx*dx+dy*dy+dz*dz)//半径range外
				{
					phiX[idx] = 0.;
					phiY[idx] = 0.;
					phiZ[idx] = 0.;
				}
				else//半径range内
				{
					phiX[idx] = order(dxh, dy, dz, RADIUS, inv_xi);
					phiY[idx] = order(dx, dyh, dz, RADIUS, inv_xi);
					phiZ[idx] = order(dx, dy, dzh, RADIUS, inv_xi);
				}
			}
		}
		
		
#pragma acc kernels
#pragma acc loop independent
		for(int n=0;n<N;n++)
		{
			double sum_phixn = 0.;
			double sum_phiyn = 0.;
			double sum_phizn = 0.;
			int offset = n*N_range*N_range*N_range;
			for(int i=0;i<N_range;i++)
			{
				for(int j=0;j<N_range;j++)
				{
					for(int k=0;k<N_range;k++)
					{
						int idx = offset+i*(N_range*N_range)+j*N_range+k;
						sum_phixn += phiX[idx];
						sum_phiyn += phiY[idx];
						sum_phizn += phiZ[idx];
					}
				}
			}
			sum_phix[n] = sum_phixn;
			sum_phiy[n] = sum_phiyn;
			sum_phiz[n] = sum_phizn;
		}
		
		
		
		
// -- c. calculation of viscosity fields --
#pragma acc parallel loop collapse(3)
		for(int i=0;i<Nx;i++)
		{
			for(int j=0;j<Ny;j++)
			{
				for(int k=0;k<Nz;k++)
				{
					eta[i][j][k]   = 1.;
					etaXY[i][j][k] = 1.;
					etaYZ[i][j][k] = 1.;
					etaZX[i][j][k] = 1.;
				}
			}
		}
		
		
#ifdef _8BLOCKLOOP_
		for(int b=0;b<8;b++)
		{
			int ib = b/4;
			int jb = (b%4)/2;
			int kb = b%2;
			//printf("b = %d, ib = %d, jb = %d, kb = %d\n", b, ib, jb, kb);
			
			#pragma acc parallel loop collapse(2)
			for(int n=0;n<N;n++)
			{
				for(int ijk=0; ijk<range*range*range; ijk++)
				{
					int in = (int)Rx[n]; int jn = (int)Ry[n]; int kn = (int)Rz[n];
					
					int i = ib*range + int(ijk/(range*range));
					int j = jb*range + int(ijk%(range*range)/range);
					int k = kb*range + int(ijk%range);
					
					int ir = i+in-range_m1;
					int irP = (ir+Nx)%Nx;
					double dx = ir - Rx[n];
					double dxh = dx + 0.5;
					int jr = j+jn-range_m1;
					int jrP = (jr+Ny)%Ny;
					double dy = jr - Ry[n];
					double dyh = dy + 0.5;
					int kr = k+kn-range_m1;
					int krP = (kr+Nz)%Nz;
					double dz = kr - Rz[n];
					double dzh = dz + 0.5;
					
					if(range2<dx*dx+dy*dy+dz*dz) continue;
					eta[irP][jrP][krP] += ratio_eta_m1*order(dx, dy, dz, RADIUS, inv_xi);
					etaXY[irP][jrP][krP] += ratio_eta_m1*order(dxh, dyh, dz, RADIUS, inv_xi);
					etaYZ[irP][jrP][krP] += ratio_eta_m1*order(dx, dyh, dzh, RADIUS, inv_xi);
					etaZX[irP][jrP][krP] += ratio_eta_m1*order(dxh, dy, dzh, RADIUS, inv_xi);
				}
			}
		}
#else // _8BLOCKLOOP_
		for(int n=0;n<N;n++)
		{
#pragma acc parallel loop collapse(3) 
			for(int i=0;i<N_range;i++)
			{
				for(int j=0;j<N_range;j++)
				{
					for(int k=0;k<N_range;k++)
					{
						double Rnx = Rx[n]; double Rny = Ry[n]; double Rnz = Rz[n];
						int in = (int)Rnx; int jn = (int)Rny; int kn = (int)Rnz;
						
						int ir = i+in-range_m1;
						int irP = (ir+Nx)%Nx;
						double dx = ir - Rnx;
						double dxh = dx + 0.5;
						int jr = j+jn-range_m1;
						int jrP = (jr+Ny)%Ny;
						double dy = jr - Rny;
						double dyh = dy + 0.5;
						int kr = k+kn-range_m1;
						int krP = (kr+Nz)%Nz;
						double dz = kr - Rnz;
						double dzh = dz + 0.5;
						
						if(range2<dx*dx+dy*dy+dz*dz) continue;
						eta[irP][jrP][krP] += ratio_eta_m1*order(dx, dy, dz, RADIUS, inv_xi);
						etaXY[irP][jrP][krP] += ratio_eta_m1*order(dxh, dyh, dz, RADIUS, inv_xi);
						etaYZ[irP][jrP][krP] += ratio_eta_m1*order(dx, dyh, dzh, RADIUS, inv_xi);
						etaZX[irP][jrP][krP] += ratio_eta_m1*order(dxh, dy, dzh, RADIUS, inv_xi);
					}
				}
			}
		}
#endif
		
// -- b. calculation of forces acting on particles --
		#pragma acc update host(Rx[:N],Ry[:N],Rz[:N])
		
		for(int n=0;n<N;n++)
		{
			R[n].x = Rx[n];
			R[n].y = Ry[n];
			R[n].z = Rz[n];
		}
		// interparticle forces
		updateList(cells,R,N);
		
		cul_LJ6_cellList3D(sgm,eps,cutoff,cells,nnIndex,R,F,box,N);
		//cul_LJn_cellList3D_notail(sgm,eps,6,cutoff,cells,nnIndex,R,F,box,N);
		//cul_WCA_cellList3D(sgm,eps,cells,nnIndex,R,F,box,N);	

		for(int n=0;n<N;n++)
		{
			Fx[n] = F[n].x;
			Fy[n] = F[n].y;
			Fz[n] = F[n].z;
		}
		
		#pragma acc update device(Fx[:N],Fy[:N],Fz[:N])
		
		
// -- d. calculation of force field acting on fluid --
#pragma acc parallel loop collapse(3)
		for(int i=0;i<Nx;i++)
		{
			for(int j=0;j<Ny;j++)
			{
				for(int k=0;k<Nz;k++)
				{
					fx[i][j][k]=0;
					fy[i][j][k]=0;
					fz[i][j][k]=0;
				}
			}
		}
		
		
#ifdef _8BLOCKLOOP_
		for(int b=0;b<8;b++)
		{
			int ib = b/4;
			int jb = (b%4)/2;
			int kb = b%2;
			
			#pragma acc parallel loop collapse(2)
			for(int n=0;n<N;n++)
			{
				for(int ijk=0; ijk<range*range*range; ijk++)
				{
					
					int in = (int)Rx[n]; int jn = (int)Ry[n]; int kn = (int)Rz[n];
					
					int i = ib*range + int(ijk/(range*range));
					int j = jb*range + int(ijk%(range*range)/range);
					int k = kb*range + int(ijk%range);
					
					int ir = i+in-range_m1;
					int irp = (ir+Nx)%Nx;
					double dx = ir - Rx[n];
					int jr = j+jn-range_m1;
					int jrp = (jr+Ny)%Ny;
					double dy = jr - Ry[n];
					int kr = k+kn-range_m1;
					int krp = (kr+Nz)%Nz;
					double dz = kr - Rz[n];
					
					if(range2<dx*dx+dy*dy+dz*dz) continue;
					
					int idx = n*(N_range*N_range*N_range)+i*(N_range*N_range)+j*N_range+k;
					fx[irp][jrp][krp] += Fx[n]*phiX[idx]/sum_phix[n];
					fy[irp][jrp][krp] += Fy[n]*phiY[idx]/sum_phiy[n]; 
					fz[irp][jrp][krp] += Fz[n]*phiZ[idx]/sum_phiz[n];
				}
			}
		}
#else // _8BLOCKLOOP_
		for(int n=0;n<N;n++)
		{
#pragma acc parallel loop collapse(3)
			for(int i=0;i<N_range;i++)
			{
				for(int j=0;j<N_range;j++)
				{
					for(int k=0;k<N_range;k++)
					{
						double Rnx = Rx[n]; double Rny = Ry[n]; double Rnz = Rz[n];
						int in = (int)Rnx; int jn = (int)Rny; int kn = (int)Rnz;
						
						int ir = i+in-range_m1;
						int irP = (ir+Nx)%Nx;
						double dx = ir - Rnx; 
						int jr = j+jn-range_m1;
						int jrP = (jr+Ny)%Ny;
						double dy = jr - Rny; 
						int kr = k+kn-range_m1;
						int krP = (kr+Nz)%Nz;
						double dz = kr - Rnz; 
						
						int idx = n*(N_range*N_range*N_range)+i*(N_range*N_range)+j*N_range+k;
						if(range2<dx*dx+dy*dy+dz*dz) continue;
						fx[irP][jrP][krP] += Fx[n]*phiX[idx]/sum_phix[n];
						fy[irP][jrP][krP] += Fy[n]*phiY[idx]/sum_phiy[n]; 
						fz[irP][jrP][krP] += Fz[n]*phiZ[idx]/sum_phiz[n]; 
					}
				}
			}
		}
		
#endif		
		
// -- e. time evolution of Navier-Stokes equation --
#pragma acc parallel loop collapse(3)
		for(int i=0;i<Nx;i++)
		{
			for(int j=0;j<Ny;j++)
			{
				for(int k=0;k<Nz;k++)
				{
					int ip=(i+1+Nx)%Nx; int im=(i-1+Nx)%Nx;
					int jp=(j+1+Ny)%Ny; int jm=(j-1+Ny)%Ny;
					int kp=(k+1+Nz)%Nz; int km=(k-1+Nz)%Nz;
					
					double eta_ijk = eta[i][j][k];
					double eta_IJk = etaXY[i][j][k];
					double eta_iJK = etaYZ[i][j][k];
					double eta_IjK = etaZX[i][j][k];
					
					double vx_ijk = vx[i][j][k]; double vy_ijk = vy[i][j][k]; double vz_ijk = vz[i][j][k]; 
					
					// advective term
					pi_dx[i][j][k] = (vx_ijk + vx[im][j][k])*(vx_ijk + vx[im][j][k])*0.25;
					pi_dy[i][j][k] = (vy_ijk + vy[i][jm][k])*(vy_ijk + vy[i][jm][k])*0.25;
					pi_dz[i][j][k] = (vz_ijk + vz[i][j][km])*(vz_ijk + vz[i][j][km])*0.25;
					pi_nz[i][j][k] = (vx_ijk + vx[i][jp][k])*(vy_ijk + vy[ip][j][k])*0.25;
					pi_nx[i][j][k] = (vy_ijk + vy[i][j][kp])*(vz_ijk + vz[i][jp][k])*0.25;
					pi_ny[i][j][k] = (vz_ijk + vz[ip][j][k])*(vx_ijk + vx[i][j][kp])*0.25;
					
					// shear stress term
					pi_dx[i][j][k] -= eta_ijk*2.*(vx_ijk - vx[im][j][k]);
					pi_dy[i][j][k] -= eta_ijk*2.*(vy_ijk - vy[i][jm][k]);
					pi_dz[i][j][k] -= eta_ijk*2.*(vz_ijk - vz[i][j][km]);
					pi_nz[i][j][k] -= eta_IJk*((vx[i][jp][k] - vx_ijk) + (vy[ip][j][k] - vy_ijk));
					pi_nx[i][j][k] -= eta_iJK*((vy[i][j][kp] - vy_ijk) + (vz[i][jp][k] - vz_ijk));
					pi_ny[i][j][k] -= eta_IjK*((vz[ip][j][k] - vz_ijk) + (vx[i][j][kp] - vx_ijk));
					
				}
			}
		}
		
		
		// noise term
#ifdef _NOISE_
#pragma acc host_data use_device(randD,randN)
{
		fpd_curand_normal(randD, randSize, 0.0, 1.0);
		fpd_curand_normal(randN, randSize, 0.0, 1.0);
}
		
#pragma acc parallel loop collapse(3)
		for(int i=0;i<Nx;i++)
		{
			for(int j=0;j<Ny;j++)
			{
				for(int k=0;k<Nz;k++)
				{
					double eta_ijk = eta[i][j][k];
					double eta_IJk = etaXY[i][j][k];
					double eta_iJK = etaYZ[i][j][k];
					double eta_IjK = etaZX[i][j][k];

					//random stress
					double sqrt2_eta_ijk_W = sqrt(2.*eta_ijk)*W;
					int idx = i*(Ny*Nz)+j*Nz+k;
					int Nxyz = Nx*Ny*Nz;
					pi_dx[i][j][k] -= sqrt2_eta_ijk_W*randD[0*Nxyz+idx];
					pi_dy[i][j][k] -= sqrt2_eta_ijk_W*randD[1*Nxyz+idx];
					pi_dz[i][j][k] -= sqrt2_eta_ijk_W*randD[2*Nxyz+idx];
					pi_nz[i][j][k] -= sqrt(eta_IJk)*W*randN[0*Nxyz+idx];
					pi_nx[i][j][k] -= sqrt(eta_iJK)*W*randN[1*Nxyz+idx];
					pi_ny[i][j][k] -= sqrt(eta_IjK)*W*randN[2*Nxyz+idx];
				}
			}
		}
#endif
		
		
		// divergence of pi
#pragma acc parallel loop collapse(3)
		for(int i=0;i<Nx;i++)
		{
			for(int j=0;j<Ny;j++)
			{
				for(int k=0;k<Nz;k++)
				{
					int ip=(i+1+Nx)%Nx; int im=(i-1+Nx)%Nx;
					int jp=(j+1+Ny)%Ny; int jm=(j-1+Ny)%Ny;
					int kp=(k+1+Nz)%Nz; int km=(k-1+Nz)%Nz;
					
					div_pix[i][j][k] = (pi_dx[ip][j][k] - pi_dx[i][j][k])
									 + (pi_nz[i][j][k] - pi_nz[i][jm][k])
									 + (pi_ny[i][j][k] - pi_ny[i][j][km]);
					div_piy[i][j][k] = (pi_nz[i][j][k] - pi_nz[im][j][k])
									 + (pi_dy[i][jp][k] - pi_dy[i][j][k])
									 + (pi_nx[i][j][k] - pi_nx[i][j][km]);
					div_piz[i][j][k] = (pi_ny[i][j][k] - pi_ny[im][j][k])
									 + (pi_nx[i][j][k] - pi_nx[i][jm][k])
									 + (pi_dz[i][j][kp] - pi_dz[i][j][k]);
									 
				}
			}
		}
		
		
		// source term of poisson equation
#pragma acc parallel loop collapse(3)
		for(int i=0;i<Nx;i++)
		{
			for(int j=0;j<Ny;j++)
			{
				for(int k=0;k<Nz;k++)
				{
					
					int im=(i-1+Nx)%Nx;
					int jm=(j-1+Ny)%Ny;
					int km=(k-1+Nz)%Nz;
					int ijk = (i*Ny+j)*Nz+k;
					
					fftinput[ijk*2] = inv_dt*(
						(vx[i][j][k] - vx[im][j][k]) + 
						(vy[i][j][k] - vy[i][jm][k]) + 
						(vz[i][j][k] - vz[i][j][km])
						) + 
						(-1)*( 
						(div_pix[i][j][k] - div_pix[im][j][k]) + 
						(div_piy[i][j][k] - div_piy[i][jm][k]) + 
						(div_piz[i][j][k] - div_piz[i][j][km])
						) + 
						(
						(fx[i][j][k] - fx[im][j][k]) + 
						(fy[i][j][k] - fy[i][jm][k]) + 
						(fz[i][j][k] - fz[i][j][km])
						);	
					fftinput[ijk*2+1] = 0;
				}
			}
		}
		
		
		
		// solve poisson equation
		//--- cufft -----
		
#pragma acc host_data use_device(fftinput)
{
		void *stream = acc_get_cuda_stream(acc_async_sync);
		cufft_c2c_forward(fftinput, Nx, Ny, Nz, stream);
}
		
		
#pragma acc kernels
#pragma acc loop independent
		for(int ijk=0;ijk<Nx*Ny*Nz;ijk++)
		{
			int i = ijk/(Ny*Nz);
			int j = ijk%(Ny*Nz)/Nz;
			int k = ijk%Nz;
			if(i==0 && j==0 && k==0){continue;}
			double phasex=2.*M_PI/(double)Nx;
			double phasey=2.*M_PI/(double)Ny;
			double phasez=2.*M_PI/(double)Nz;
			double cosx=cos(phasex*(double)i);
			double cosy=cos(phasey*(double)j);
			double cosz=cos(phasez*(double)k);
			double nrm = 0.5/(cosx + cosy + cosz - 3.);
			fftinput[ijk*2] *= nrm;
			fftinput[ijk*2+1] *= nrm;
			
		}
		
		
#pragma acc host_data use_device(fftinput)
{
		void *stream = acc_get_cuda_stream(acc_async_sync);
		cufft_c2c_inverse(fftinput, Nx, Ny, Nz, stream);
}
		
#pragma acc parallel loop collapse(3)
		for(int i=0;i<Nx;i++)
		{
			for(int j=0;j<Ny;j++)
			{
				for(int k=0;k<Nz;k++)
				{
					p[i][j][k] = fftinput[2*((i*Ny+j)*Nz+k)]/((double)(Nx*Ny*Nz));
				}
			}
		}
		
		
		// calculation of velosity field in the next time step 
#pragma acc parallel loop collapse(3)
		for(int i=0;i<Nx;i++)
		{
			for(int j=0;j<Ny;j++)
			{
				for(int k=0;k<Nz;k++)
				{
					int ip=(i+1+Nx)%Nx;
					int jp=(j+1+Ny)%Ny;
					int kp=(k+1+Nz)%Nz;
					vx[i][j][k] += DT*(
						-(p[ip][j][k] - p[i][j][k])
						-div_pix[i][j][k]
						+fx[i][j][k]
						);
					vy[i][j][k] += DT*(
						-(p[i][jp][k] - p[i][j][k])
						-div_piy[i][j][k]
						+fy[i][j][k]
						);
					vz[i][j][k] += DT*(
						-(p[i][j][kp] - p[i][j][k])
						-div_piz[i][j][k]
						+fz[i][j][k]
						);
				}
			}
		}
		
		
		
		
// -- f. calculation of velosities of particles --		
		
#pragma acc parallel loop collapse(2)
		for(int n=0;n< N;n++)
		{
			for(int ijk=0; ijk<N_range*N_range*N_range; ijk++)
			{
				int in = (int)Rx[n]; int jn = (int)Ry[n]; int kn = (int)Rz[n];
				
				int i = int(ijk/(N_range*N_range));
				int j = int(ijk%(N_range*N_range)/N_range);
				int k = int(ijk%N_range);
				
				int ir = i+in-range_m1; 
				double dx = ir - Rx[n]; 
				int irp = (ir+Nx)%Nx;
				
				int jr = j+jn-range_m1; 
				double dy = jr - Ry[n]; 
				int jrp = (jr+Ny)%Ny;
					
				int kr = k+kn-range_m1;
				double dz = kr - Rz[n]; 
				int krp = (kr+Nz)%Nz;
						
				// note:this "if" gvies 1/10,000 difference in sum_phi with a=3.2 and xi = 1.
				if(range2<dx*dx+dy*dy+dz*dz) continue;	

				int idx = n*(N_range*N_range*N_range)+i*(N_range*N_range)+j*N_range+k;

				phiX[idx] *= vx[irp][jrp][krp];
				phiY[idx] *= vy[irp][jrp][krp];
				phiZ[idx] *= vz[irp][jrp][krp];
			}
		}
		
		
		
#pragma acc kernels
#pragma acc loop independent
		for(int n=0;n< N;n++)
		{
			double sum_Vx = 0.;
			double sum_Vy = 0.;
			double sum_Vz = 0.;
			int offset = n*N_range*N_range*N_range;
			
			for(int i=0;i<N_range;i++)
			{
				for(int j=0;j<N_range;j++)
				{
					for(int k=0;k<N_range;k++)
					{
						int idx = offset+i*(N_range*N_range)+j*N_range+k;
						sum_Vx += phiX[idx];
						sum_Vy += phiY[idx];
						sum_Vz += phiZ[idx];
					}
				}
			}
			
			Vx[n] = sum_Vx/sum_phix[n];
			Vy[n] = sum_Vy/sum_phiy[n];
			Vz[n] = sum_Vz/sum_phiz[n];
		}
			
		
		
// -- g. time evolution of particles' positions --			
#pragma acc kernels
#pragma acc loop independent
		for(int n=0;n<N;n++)
		{
			Rx[n] += DT*Vx[n]; 			
			Ry[n] += DT*Vy[n];
			Rz[n] += DT*Vz[n];
				
			//PBC	
			if(Rx[n]<0) Rx[n] += box.x; else if(Rx[n]>box.x) Rx[n] -= box.x;			
			if(Ry[n]<0) Ry[n] += box.y; else if(Ry[n]>box.y) Ry[n] -= box.y;
			if(Rz[n]<0) Rz[n] += box.z; else if(Rz[n]>box.z) Rz[n] -= box.z;
		}
		
		
		
		
	} // the end of time loop
	
	
}// acc data
	
	
	// memroy release
	Free1D(R);Free1D(V);Free1D(F);
	Free1D(Rx);Free1D(Ry);Free1D(Rz);
	Free1D(Vx);Free1D(Vy);Free1D(Vz);
	Free1D(Fx);Free1D(Fy);Free1D(Fz);
	Free3D(v);
	
	Free3D(fx);Free3D(fy);Free3D(fz);
	Free1D(phiX);Free1D(phiY);Free1D(phiZ);
	Free3D(eta);Free3D(etaXY);Free3D(etaYZ);Free3D(etaZX);
	Free3D(vx);Free3D(vy);Free3D(vz);
	Free3D(pi_dx);Free3D(pi_dy);Free3D(pi_dz);
	Free3D(pi_nx);Free3D(pi_ny);Free3D(pi_nz);
	Free3D(div_pix);Free3D(div_piy);Free3D(div_piz);
	Free1D(fftinput);
	
#ifdef _NOISE_
	fpd_curand_fin();
	Free1D(randD);
	Free1D(randN);
#endif
		
	std::cout << "end time: " << std::ctime(&timer) << std::endl;
	return 0;
}
	
	
	
	
// cufft
void cufft_c2c_forward(double *d_data, int Nx, int Ny, int Nz, void *stream)
{
    cufftHandle plan;
    cufftPlan3d(&plan, Nx, Ny, Nz, CUFFT_Z2Z);
    cufftSetStream(plan, (cudaStream_t)stream);
    cufftExecZ2Z(plan, (cufftDoubleComplex*)d_data, (cufftDoubleComplex*)d_data,CUFFT_FORWARD);
    cufftDestroy(plan);
}

void cufft_c2c_inverse(double *d_data, int Nx, int Ny, int Nz, void *stream)
{
    cufftHandle plan;
    cufftPlan3d(&plan, Nx, Ny, Nz, CUFFT_Z2Z);
    cufftSetStream(plan, (cudaStream_t)stream);
    cufftExecZ2Z(plan, (cufftDoubleComplex*)d_data, (cufftDoubleComplex*)d_data,CUFFT_INVERSE);
    cufftDestroy(plan);
}

	
// randam number generator
void fpd_curand_init(unsigned long int seed)
{
	curandCreateGenerator(&generator, CURAND_RNG_PSEUDO_XORWOW);
	curandSetPseudoRandomGeneratorSeed(generator, seed);
}

void fpd_curand_fin(void)
{
	curandDestroyGenerator(generator);
}

void fpd_curand_normal(float *rand, int n, float mean, float stddev)
{
	curandGenerateNormal(generator, rand, n, mean, stddev);
	cudaDeviceSynchronize();
}

