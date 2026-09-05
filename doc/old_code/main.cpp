#include <stdio.h>
#include <math.h>
#include <ctime>
#include <chrono>
#include <vector>

#include <iostream>
#include <fstream>

#include "global_definitions.h"
#include "include/vector_mine.h"
#include "include/io.h"
#include "include/cellList.h"
#include "include/interaction.h"

#include <openacc.h>
#include <curand.h>
#include <cufft.h>

#define IDX3(i, j, k,Nx ,Ny, Nz) ((i) + (j)*Nx + (k)*Nx*Ny)//一维转三维
#define IDX2(i, j,Nx ,Ny) ((i) + (j)*Nx)//一维转二维

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
	auto start = std::chrono::high_resolution_clock::now(); // 开始时间
		
		
	if(argc!=6)
	{
		std::cout << argv[0] <<" [init state] [output directory] [time range] [beta] [g]" << std::endl;
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
	const double g = atof(argv[5]);
	int endTime = atoi(argv[3]);//模拟的结束时间步数
	const double inv_beta = 1./atof(argv[4]);//beta的倒数，beta=1/(k_B*T)
	
	const double inv_dt = 1./DT;//DT在global_definitions.h中定义，这是它的倒数
	const double inv_xi = 1./XI;//XT在global_definitions.h中定义，这是它的倒数

	const int interval_pos = INTERVAL_POSITIONS*(int)(inv_dt);//每隔多少时间步保存一次粒子的位置数据？
	//const int interval_pos = 1;//每隔多少时间步保存一次粒子的位置数据？
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

	double *fx = (double*)calloc(Nx*Ny*Nz, sizeof(double));//储流体在x方向的受力场
	double *fy = (double*)calloc(Nx*Ny*Nz, sizeof(double));//储流体在y方向的受力场
	double *fz = (double*)calloc(Nx*Ny*Nz, sizeof(double));//储流体在z方向的受力场

	// eta, etaXY, etaYZ, etaZX: 分别存储流体的粘度场及其各方向分量。

	double *eta = (double*)calloc(Nx*Ny*Nz, sizeof(double));//流体的粘度场分量
	double *etaXY = (double*)calloc(Nx*Ny*Nz, sizeof(double));//流体的粘度场分量
	double *etaYZ = (double*)calloc(Nx*Ny*Nz, sizeof(double));//流体的粘度场分量
	double *etaZX = (double*)calloc(Nx*Ny*Nz, sizeof(double));//流体的粘度场分量
	
	// vx, vy, vz: 分别存储流体在x、y、z方向的速度场。

	double *vx = (double*)calloc(Nx*Ny*Nz, sizeof(double));
	double *vy = (double*)calloc(Nx*Ny*Nz, sizeof(double));
	double *vz = (double*)calloc(Nx*Ny*Nz, sizeof(double));

	vector3D *v = (vector3D*)calloc(Nx * Ny * Nz, sizeof(vector3D));//储流体的速度场
	double *p = (double*)calloc(Nx*Ny*Nz, sizeof(double));//p流体的压力场。
	
	//分别存储流体的应力张量场

	double *pi_dx = (double*)calloc(Nx*Ny*Nz, sizeof(double));
	double *pi_dy = (double*)calloc(Nx*Ny*Nz, sizeof(double));
	double *pi_dz = (double*)calloc(Nx*Ny*Nz, sizeof(double));
	double *pi_nx = (double*)calloc(Nx*Ny*Nz, sizeof(double));
	double *pi_ny = (double*)calloc(Nx*Ny*Nz, sizeof(double));
	double *pi_nz = (double*)calloc(Nx*Ny*Nz, sizeof(double));
	
	//分别存储流体的应力张量散度场

	double *div_pix = (double*)calloc(Nx*Ny*Nz, sizeof(double));
	double *div_piy = (double*)calloc(Nx*Ny*Nz, sizeof(double));
	double *div_piz = (double*)calloc(Nx*Ny*Nz, sizeof(double));
	
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
		//printf("%lf, %lf, %lf\n", Rx[n], Ry[n], Rz[n]);
	}

	for(int i=0;i<Nx;i++)
	{
		for(int j=0;j<Ny;j++)
		{
			for(int k=0;k<Nz;k++)
			{
				vx[IDX3(i, j, k,Nx ,Ny, Nz)] = v[IDX3(i, j, k,Nx ,Ny, Nz)].x;
				vy[IDX3(i, j, k,Nx ,Ny, Nz)] = v[IDX3(i, j, k,Nx ,Ny, Nz)].y;
				vz[IDX3(i, j, k,Nx ,Ny, Nz)] = v[IDX3(i, j, k,Nx ,Ny, Nz)].z;
			}
		}
	}
	

	// cell list
	double cutoff = CUTOFF_POTENTIAL*Re; 
    cellList3D *cells = getList(box,CELL_SIZE,N);
    std::vector<int> nnIndex = get_nnIndex(cells);
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
	std::cout << "g=" << g << std::endl;
	std::cout << "-------------------------------------"<< std::endl; 
	
	
	#pragma acc data create(Vx[:N],Vy[:N],Vz[:N])
	#pragma acc data create(Fx[:N],Fy[:N],Fz[:N])
	#pragma acc data copyin(Rx[:N],Ry[:N],Rz[:N])
	#pragma acc enter data copyin(vx[0:Nx*Ny*Nz],vy[0:Nx*Ny*Nz],vz[0:Nx*Ny*Nz])
	#pragma acc data create(pi_dx[0:Nx*Ny*Nz],pi_dy[0:Nx*Ny*Nz],pi_dz[0:Nx*Ny*Nz])
	#pragma acc data create(pi_nx[0:Nx*Ny*Nz],pi_ny[0:Nx*Ny*Nz],pi_nz[0:Nx*Ny*Nz])
	#pragma acc data create(div_pix[0:Nx*Ny*Nz],div_piy[0:Nx*Ny*Nz],div_piz[0:Nx*Ny*Nz])
	#pragma acc data create(fx[0:Nx*Ny*Nz],fy[0:Nx*Ny*Nz],fz[0:Nx*Ny*Nz])
	#pragma acc data create(eta[0:Nx*Ny*Nz],etaXY[0:Nx*Ny*Nz],etaYZ[0:Nx*Ny*Nz],etaZX[0:Nx*Ny*Nz])
	#pragma acc data create(p[0:Nx*Ny*Nz])
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
		//printf("%d ",t);
		// save 
		if(t%interval_pos==0)
		//if(t%1==0)
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
			//sprintf(filename,"%s/pos_%05d",argv[2],t);
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
			#pragma acc update host(vx[0:Nx*Ny*Nz],vy[0:Nx*Ny*Nz],vz[0:Nx*Ny*Nz])
			for(int i=0;i<Nx;i++)
			{
				for(int j=0;j<Ny;j++)
				{
					for(int k=0;k<Nz;k++)
					{
						v[IDX3(i, j, k,Nx ,Ny, Nz)].x = vx[IDX3(i, j, k,Nx ,Ny, Nz)];
						v[IDX3(i, j, k,Nx ,Ny, Nz)].y = vy[IDX3(i, j, k,Nx ,Ny, Nz)];
						v[IDX3(i, j, k,Nx ,Ny, Nz)].z = vz[IDX3(i, j, k,Nx ,Ny, Nz)];
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
		
		
//#pragma acc kernels
//#pragma acc loop independent
#pragma acc parallel loop collapse(1)
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
					eta[IDX3(i, j, k,Nx ,Ny, Nz)]   = 1.;
					etaXY[IDX3(i, j, k,Nx ,Ny, Nz)] = 1.;
					etaYZ[IDX3(i, j, k,Nx ,Ny, Nz)] = 1.;
					etaZX[IDX3(i, j, k,Nx ,Ny, Nz)] = 1.;
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
					eta[IDX3(irP, jrP, krP,Nx ,Ny, Nz)] += ratio_eta_m1*order(dx, dy, dz, RADIUS, inv_xi);
					etaXY[IDX3(irP, jrP, krP,Nx ,Ny, Nz)] += ratio_eta_m1*order(dxh, dyh, dz, RADIUS, inv_xi);
					etaYZ[IDX3(irP, jrP, krP,Nx ,Ny, Nz)] += ratio_eta_m1*order(dx, dyh, dzh, RADIUS, inv_xi);
					etaZX[IDX3(irP, jrP, krP,Nx ,Ny, Nz)] += ratio_eta_m1*order(dxh, dy, dzh, RADIUS, inv_xi);
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
						eta[IDX3(irP, jrP, krP,Nx ,Ny, Nz)] += ratio_eta_m1*order(dx, dy, dz, RADIUS, inv_xi);
						etaXY[IDX3(irP, jrP, krP,Nx ,Ny, Nz)] += ratio_eta_m1*order(dxh, dyh, dz, RADIUS, inv_xi);
						etaYZ[IDX3(irP, jrP, krP,Nx ,Ny, Nz)] += ratio_eta_m1*order(dx, dyh, dzh, RADIUS, inv_xi);
						etaZX[IDX3(irP, jrP, krP,Nx ,Ny, Nz)] += ratio_eta_m1*order(dxh, dy, dzh, RADIUS, inv_xi);
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
		//cul_LJ6_cellList3D(sgm,eps,cutoff,cells,nnIndex,R,F,box,N);
		Morse_force(De,a,Re,cutoff,cells,nnIndex,R,F,box,N);
		for(int n=0;n<N;n++)
		{
			Fx[n] = F[n].x;
			Fy[n] = F[n].y;
			Fz[n] = F[n].z-10;//手动加个重力  <--------------------------------------------------------------------
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
					fx[IDX3(i, j, k,Nx ,Ny, Nz)]=0;
					fy[IDX3(i, j, k,Nx ,Ny, Nz)]=0;
					fz[IDX3(i, j, k,Nx ,Ny, Nz)]=0;
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
					fx[IDX3(irp, jrp, krp,Nx ,Ny, Nz)] += Fx[n]*phiX[idx]/sum_phix[n];
					fy[IDX3(irp, jrp, krp,Nx ,Ny, Nz)] += Fy[n]*phiY[idx]/sum_phiy[n]; 
					fz[IDX3(irp, jrp, krp,Nx ,Ny, Nz)] += Fz[n]*phiZ[idx]/sum_phiz[n];
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
						fx[IDX3(irP, jrP, krP,Nx ,Ny, Nz)] += Fx[n]*phiX[idx]/sum_phix[n];
						fy[IDX3(irP, jrP, krP,Nx ,Ny, Nz)] += Fy[n]*phiY[idx]/sum_phiy[n]; 
						fz[IDX3(irP, jrP, krP,Nx ,Ny, Nz)] += Fz[n]*phiZ[idx]/sum_phiz[n]; 
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
					
					double eta_ijk = eta[IDX3(i, j, k,Nx ,Ny, Nz)];
					double eta_IJk = etaXY[IDX3(i, j, k,Nx ,Ny, Nz)];
					double eta_iJK = etaYZ[IDX3(i, j, k,Nx ,Ny, Nz)];
					double eta_IjK = etaZX[IDX3(i, j, k,Nx ,Ny, Nz)];
					
					double vx_ijk = vx[IDX3(i, j, k,Nx ,Ny, Nz)]; double vy_ijk = vy[IDX3(i, j, k,Nx ,Ny, Nz)]; double vz_ijk = vz[IDX3(i, j, k,Nx ,Ny, Nz)]; 
					
					// advective term
					pi_dx[IDX3(i, j, k,Nx ,Ny, Nz)] = (vx_ijk + vx[IDX3(im, j, k,Nx ,Ny, Nz)])*(vx_ijk + vx[IDX3(im, j, k,Nx ,Ny, Nz)])*0.25;
					pi_dy[IDX3(i, j, k,Nx ,Ny, Nz)] = (vy_ijk + vy[IDX3(i, jm, k,Nx ,Ny, Nz)])*(vy_ijk + vy[IDX3(i, jm, k,Nx ,Ny, Nz)])*0.25;
					pi_dz[IDX3(i, j, k,Nx ,Ny, Nz)] = (vz_ijk + vz[IDX3(i, j, km,Nx ,Ny, Nz)])*(vz_ijk + vz[IDX3(i, j, km,Nx ,Ny, Nz)])*0.25;
					pi_nz[IDX3(i, j, k,Nx ,Ny, Nz)] = (vx_ijk + vx[IDX3(i, jp, k,Nx ,Ny, Nz)])*(vy_ijk + vy[IDX3(ip, j, k,Nx ,Ny, Nz)])*0.25;
					pi_nx[IDX3(i, j, k,Nx ,Ny, Nz)] = (vy_ijk + vy[IDX3(i, j, kp,Nx ,Ny, Nz)])*(vz_ijk + vz[IDX3(i, jp, k,Nx ,Ny, Nz)])*0.25;
					pi_ny[IDX3(i, j, k,Nx ,Ny, Nz)] = (vz_ijk + vz[IDX3(ip, j, k,Nx ,Ny, Nz)])*(vx_ijk + vx[IDX3(i, j, kp,Nx ,Ny, Nz)])*0.25;
					
					// shear stress term
					pi_dx[IDX3(i, j, k,Nx ,Ny, Nz)] -= eta_ijk*2.*(vx_ijk - vx[IDX3(im, j, k,Nx ,Ny, Nz)]);
					pi_dy[IDX3(i, j, k,Nx ,Ny, Nz)] -= eta_ijk*2.*(vy_ijk - vy[IDX3(i, jm, k,Nx ,Ny, Nz)]);
					pi_dz[IDX3(i, j, k,Nx ,Ny, Nz)] -= eta_ijk*2.*(vz_ijk - vz[IDX3(i, j, km,Nx ,Ny, Nz)]);
					pi_nz[IDX3(i, j, k,Nx ,Ny, Nz)] -= eta_IJk*((vx[IDX3(i, jp, k,Nx ,Ny, Nz)] - vx_ijk) + (vy[IDX3(ip, j, k,Nx ,Ny, Nz)] - vy_ijk));
					pi_nx[IDX3(i, j, k,Nx ,Ny, Nz)] -= eta_iJK*((vy[IDX3(i, j, kp,Nx ,Ny, Nz)] - vy_ijk) + (vz[IDX3(i, jp, k,Nx ,Ny, Nz)] - vz_ijk));
					pi_ny[IDX3(i, j, k,Nx ,Ny, Nz)] -= eta_IjK*((vz[IDX3(ip, j, k,Nx ,Ny, Nz)] - vz_ijk) + (vx[IDX3(i, j, kp,Nx ,Ny, Nz)] - vx_ijk));
					
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
					double eta_ijk = eta[IDX3(i, j, k,Nx ,Ny, Nz)];
					double eta_IJk = etaXY[IDX3(i, j, k,Nx ,Ny, Nz)];
					double eta_iJK = etaYZ[IDX3(i, j, k,Nx ,Ny, Nz)];
					double eta_IjK = etaZX[IDX3(i, j, k,Nx ,Ny, Nz)];

					//random stress
					double sqrt2_eta_ijk_W = sqrt(2.*eta_ijk)*W;
					int idx = i*(Ny*Nz)+j*Nz+k;
					int Nxyz = Nx*Ny*Nz;
					pi_dx[IDX3(i, j, k,Nx ,Ny, Nz)] -= sqrt2_eta_ijk_W*randD[0*Nxyz+idx];
					pi_dy[IDX3(i, j, k,Nx ,Ny, Nz)] -= sqrt2_eta_ijk_W*randD[1*Nxyz+idx];
					pi_dz[IDX3(i, j, k,Nx ,Ny, Nz)] -= sqrt2_eta_ijk_W*randD[2*Nxyz+idx];
					pi_nz[IDX3(i, j, k,Nx ,Ny, Nz)] -= sqrt(eta_IJk)*W*randN[0*Nxyz+idx];
					pi_nx[IDX3(i, j, k,Nx ,Ny, Nz)] -= sqrt(eta_iJK)*W*randN[1*Nxyz+idx];
					pi_ny[IDX3(i, j, k,Nx ,Ny, Nz)] -= sqrt(eta_IjK)*W*randN[2*Nxyz+idx];
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
					
					div_pix[IDX3(i, j, k,Nx ,Ny, Nz)] = (pi_dx[IDX3(ip, j, k,Nx ,Ny, Nz)] - pi_dx[IDX3(i, j, k,Nx ,Ny, Nz)])
														+ (pi_nz[IDX3(i, j, k,Nx ,Ny, Nz)] - pi_nz[IDX3(i, jm, k,Nx ,Ny, Nz)])
														+ (pi_ny[IDX3(i, j, k,Nx ,Ny, Nz)] - pi_ny[IDX3(i, j, km,Nx ,Ny, Nz)]);
					div_piy[IDX3(i, j, k,Nx ,Ny, Nz)] = (pi_nz[IDX3(i, j, k,Nx ,Ny, Nz)] - pi_nz[IDX3(im, j, k,Nx ,Ny, Nz)])
														+ (pi_dy[IDX3(i, jp, k,Nx ,Ny, Nz)] - pi_dy[IDX3(i, j, k,Nx ,Ny, Nz)])
														+ (pi_nx[IDX3(i, j, k,Nx ,Ny, Nz)] - pi_nx[IDX3(i, j, km,Nx ,Ny, Nz)]);
					div_piz[IDX3(i, j, k,Nx ,Ny, Nz)] = (pi_ny[IDX3(i, j, k,Nx ,Ny, Nz)] - pi_ny[IDX3(im, j, k,Nx ,Ny, Nz)])
														+ (pi_nx[IDX3(i, j, k,Nx ,Ny, Nz)] - pi_nx[IDX3(i, jm, k,Nx ,Ny, Nz)])
														+ (pi_dz[IDX3(i, j, kp,Nx ,Ny, Nz)] - pi_dz[IDX3(i, j, k,Nx ,Ny, Nz)]);
									 
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
						(vx[IDX3(i, j, k,Nx ,Ny, Nz)] - vx[IDX3(im, j, k,Nx ,Ny, Nz)]) + 
						(vy[IDX3(i, j, k,Nx ,Ny, Nz)] - vy[IDX3(i, jm, k,Nx ,Ny, Nz)]) + 
						(vz[IDX3(i, j, k,Nx ,Ny, Nz)] - vz[IDX3(i, j, km,Nx ,Ny, Nz)])
						) + 
						(-1)*( 
						(div_pix[IDX3(i, j, k,Nx ,Ny, Nz)] - div_pix[IDX3(im, j, k,Nx ,Ny, Nz)]) + 
						(div_piy[IDX3(i, j, k,Nx ,Ny, Nz)] - div_piy[IDX3(i, jm, k,Nx ,Ny, Nz)]) + 
						(div_piz[IDX3(i, j, k,Nx ,Ny, Nz)] - div_piz[IDX3(i, j, km,Nx ,Ny, Nz)])
						) + 
						(
						(fx[IDX3(i, j, k,Nx ,Ny, Nz)] - fx[IDX3(im, j, k,Nx ,Ny, Nz)]) + 
						(fy[IDX3(i, j, k,Nx ,Ny, Nz)] - fy[IDX3(i, jm, k,Nx ,Ny, Nz)]) + 
						(fz[IDX3(i, j, k,Nx ,Ny, Nz)] - fz[IDX3(i, j, km,Nx ,Ny, Nz)])
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
		
		
//#pragma acc kernels
//#pragma acc loop independent
#pragma acc parallel loop collapse(1)
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
					p[IDX3(i, j, k,Nx ,Ny, Nz)] = fftinput[2*((i*Ny+j)*Nz+k)]/((double)(Nx*Ny*Nz));
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
					vx[IDX3(i, j, k,Nx ,Ny, Nz)] += DT*(
						-(p[IDX3(ip, j, k,Nx ,Ny, Nz)] - p[IDX3(i, j, k,Nx ,Ny, Nz)])
						-div_pix[IDX3(i, j, k,Nx ,Ny, Nz)]
						+fx[IDX3(i, j, k,Nx ,Ny, Nz)]
						);
					vy[IDX3(i, j, k,Nx ,Ny, Nz)] += DT*(
						-(p[IDX3(i, jp, k,Nx ,Ny, Nz)] - p[IDX3(i, j, k,Nx ,Ny, Nz)])
						-div_piy[IDX3(i, j, k,Nx ,Ny, Nz)]
						+fy[IDX3(i, j, k,Nx ,Ny, Nz)]
						);
					vz[IDX3(i, j, k,Nx ,Ny, Nz)] += DT*(
						-(p[IDX3(i, j, kp,Nx ,Ny, Nz)] - p[IDX3(i, j, k,Nx ,Ny, Nz)])
						-div_piz[IDX3(i, j, k,Nx ,Ny, Nz)]
						+fz[IDX3(i, j, k,Nx ,Ny, Nz)]
						);

					// if (k==0 || k==Nz-1 )
					// {
					// 	vz[i][j][k]=0;
					// }
					
				}
			}
		}


		// 边界条件
#pragma acc parallel loop collapse(2)
		for(int i=0;i<Nx;i++)
		{
			for(int j=0;j<Ny;j++)
			{
				vz[IDX3(i, j, 0,Nx ,Ny, Nz)] = fabs(vz[IDX3(i, j, 0,Nx ,Ny, Nz)]);
				vz[IDX3(i, j, Nz-1,Nx ,Ny, Nz)] = -fabs(vz[IDX3(i, j, Nz-1,Nx ,Ny, Nz)]);
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

				phiX[idx] *= vx[IDX3(irp, jrp, krp,Nx ,Ny, Nz)];
				phiY[idx] *= vy[IDX3(irp, jrp, krp,Nx ,Ny, Nz)];
				phiZ[idx] *= vz[IDX3(irp, jrp, krp,Nx ,Ny, Nz)];
			}
		}
		
		
		
//#pragma acc kernels
//#pragma acc loop independent
#pragma acc parallel loop collapse(1)
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
//#pragma acc kernels
//#pragma acc loop independent
#pragma acc parallel loop collapse(1)
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
	// Free1D(R);Free1D(V);Free1D(F);
	// Free1D(Rx);Free1D(Ry);Free1D(Rz);
	// Free1D(Vx);Free1D(Vy);Free1D(Vz);
	// Free1D(Fx);Free1D(Fy);Free1D(Fz);
	// Free3D(v);
	
	#pragma acc exit data delete(R, V, F)
	// free(R); free(V); free(F);
	
	#pragma acc exit data delete(Rx, Ry, Rz)
	// free(Rx); free(Ry); free(Rz);
	
	#pragma acc exit data delete(Vx, Vy, Vz)
	// free(Vx); free(Vy); free(Vz);

	#pragma acc exit data delete(Fx, Fy, Fz)
	// free(Fx); free(Fy); free(Fz);

	#pragma acc exit data delete(v)
	// free(v);

	// Free3D(fx);Free3D(fy);Free3D(fz);
	// Free1D(phiX);Free1D(phiY);Free1D(phiZ);
	// Free3D(eta);Free3D(etaXY);Free3D(etaYZ);Free3D(etaZX);

	#pragma acc exit data delete(fx, fy, fz)
	// free(fx); free(fy); free(fz);
	
	#pragma acc exit data delete(phiX, phiY, phiZ)
	// free(phiX); free(phiY); free(phiZ);

	#pragma acc exit data delete(Fx, Fy, Fz)
	// free(Fx); free(Fy); free(Fz);

	#pragma acc exit data delete(etaXY, etaYZ, etaZX)
	// free(etaXY); free(etaYZ); free(etaZX);

	#pragma acc exit data delete(eta)
	// free(eta);

	// Free3D(pi_dx);Free3D(pi_dy);Free3D(pi_dz);
	// Free3D(pi_nx);Free3D(pi_ny);Free3D(pi_nz);
	// Free3D(div_pix);Free3D(div_piy);Free3D(div_piz);
	// Free1D(fftinput);

	#pragma acc exit data delete(pi_dx, pi_dy, pi_dz)
	// free(pi_dx); free(pi_dy); free(pi_dz);

	#pragma acc exit data delete(pi_nx, pi_ny, pi_nz)
	// free(pi_nx); free(pi_ny); free(pi_nz);

	#pragma acc exit data delete(div_pix, div_piy, div_piz)
	// free(div_pix); free(div_piy); free(div_piz);

	#pragma acc exit data delete(fftinput)
	// free(fftinput);
	
#ifdef _NOISE_
	fpd_curand_fin();
	// Free1D(randD);
	// Free1D(randN);
	#pragma acc exit data delete(randD)
	free(randD);
	#pragma acc exit data delete(randN)
	free(randN);
#endif
		
	std::cout << "end time: " << std::ctime(&timer) << std::endl;

	auto end = std::chrono::high_resolution_clock::now(); // 结束时间
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start); // 转换为毫秒
    std::cout << "Wall-clock time: " << duration.count() << " milliseconds" << std::endl;

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

