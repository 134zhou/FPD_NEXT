#include "State.h"
#include "Neighbor.h"
#include <openacc.h>
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cmath>
#include <vector>

int run_neighbor_check(bool bench)
{
    int fails=0;
    const int sizes[]={1000,10000,30000,100000};
    for(int test=0;test<(bench?4:24);++test)
    {
        int N=bench?sizes[test]:257;
        int L=bench?int(std::ceil(std::cbrt(double(N))*3)):24;
        NS_Config cfg=make_ns_config(L,L,L,0.001,0.0,false);
        FpdState st;
        st.init(cfg,N,ST_PARTICLE,1);
        int side=int(std::ceil(std::cbrt(double(N))));
        for(int i=0;i<N;++i)
        {
            st.Rx[i]=(i%side+0.25)*L/side;
            st.Ry[i]=((i/side)%side+0.35)*L/side;
            st.Rz[i]=(i/(side*side)+0.45)*L/side-0.5;
        }
        st.upload(ST_PARTICLE);
        int type=bench?POT_YUKAWA:1+test%4;
        double rc=bench?8.0:(test%8<4?5.0:14.0);
        PotentialParams p=make_potential_params(type,bench?SHIFT_ENERGY:test/8,1,bench?2:3,1,0.4,2,rc,0.2);
        st.prepare_neighbor(p.rcut);
        auto direct=[&]() { compute_particle_forces(cfg,p,no_wall(),N,st.Rx,st.Ry,st.Rz,st.Fx,st.Fy,st.Fz); acc_wait_all(); };
        auto cell=[&]() { compute_particle_forces_cell(cfg,p,no_wall(),N,st.Rx,st.Ry,st.Rz,st.Fx,st.Fy,st.Fz,st.neighbor); acc_wait_all(); };
        direct(); st.download(ST_PARTICLE);
        std::vector<double> dx(st.Fx,st.Fx+N),dy(st.Fy,st.Fy+N),dz(st.Fz,st.Fz+N);
        cell(); st.download(ST_PARTICLE);
        double error=0,scale=1;
        for(int i=0;i<N;++i)
        {
            error=std::max(error,std::max(std::abs(st.Fx[i]-dx[i]),std::max(std::abs(st.Fy[i]-dy[i]),std::abs(st.Fz[i]-dz[i]))));
            scale=std::max(scale,std::max(std::abs(dx[i]),std::max(std::abs(dy[i]),std::abs(dz[i]))));
        }
        if (!bench)
        {
            std::vector<double> fx(N),fy(N),fz(N);
            compute_particle_forces_cpu(cfg,p,no_wall(),N,st.Rx,st.Ry,st.Rz,fx.data(),fy.data(),fz.data());
            for (int i=0; i<N; ++i)
            {
                error=std::max(error,std::max(std::abs(st.Fx[i]-fx[i]),std::max(std::abs(st.Fy[i]-fy[i]),std::abs(st.Fz[i]-fz[i]))));
            }
        }
        std::vector<double> cx(st.Fx,st.Fx+N),cy(st.Fy,st.Fy+N),cz(st.Fz,st.Fz+N);
        cell(); st.download(ST_PARTICLE);
        bool repeat=true;
        for(int i=0;i<N;++i) repeat=repeat && cx[i]==st.Fx[i] && cy[i]==st.Fy[i] && cz[i]==st.Fz[i];
        if(error/scale>=1e-10 || !repeat) ++fails;
        std::printf("N=%d L=%d type=%d rcut=%.4g error=%.3e repeat=%s",N,L,type,p.rcut,error/scale,repeat?"PASS":"FAIL");
        if(bench)
        {
            double td[3],tc[3];
            for(int batch=0;batch<3;++batch)
            {
                auto start=std::chrono::steady_clock::now();
                for(int k=0;k<5;++k) direct();
                td[batch]=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()/5;
                start=std::chrono::steady_clock::now();
                for(int k=0;k<5;++k) cell();
                tc[batch]=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()/5;
            }
            std::sort(td,td+3);std::sort(tc,tc+3);
            std::printf(" direct_ms=%.4f cell_ms=%.4f speedup=%.3f",td[1]*1000,tc[1]*1000,td[1]/tc[1]);
        }
        std::printf("\n"); std::fflush(stdout);
        st.finish();
    }
    return fails;
}
