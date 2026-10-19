#include "Step.h"

#include "Force.h"
#include "Viscosity.h"
#include "Velocity.h"
#include "Mass.h"
#include "Stokes.h"

StepParams make_step_params(const FpdConfig& effective)
{
    StepParams sp;
    sp.cfg   = make_ns_config(effective);
    sp.pp    = make_phi_params(effective);
    sp.pot   = make_potential_params(effective);
    sp.wall  = make_wall_params(effective);
    sp.ext   = make_external_field(effective);
    sp.d_rho = effective.ratio_rho - 1.0;
    return sp;
}

void step_forces(const StepParams& sp, FpdState& st)
{
    compute_particle_forces(sp.cfg, sp.pot, sp.wall, st.N, st.Rx, st.Ry, st.Rz,
                            st.Fx, st.Fy, st.Fz);
}

void step_fields(const StepParams& sp, FpdState& st)
{
    const NS_Config cfg = sp.cfg;
    const int N = st.N;
    update_viscosity_fields(cfg, sp.pp, N, st.Rx, st.Ry, st.Rz,
                            st.sum_phix, st.sum_phiy, st.sum_phiz,
                            st.eta, st.etaXY, st.etaYZ, st.etaZX);
    // 反密度场只依赖 R；density_on=0 时跳过（数组保持 1.0）。
    if (cfg.density_on)
    {
        update_mass_fields(cfg, sp.pp, sp.d_rho, N, st.Rx, st.Ry, st.Rz,
                           st.inv_rho_x, st.inv_rho_y, st.inv_rho_z);
    }
    // 背景力密度补偿恒为 0：无滑移壁面本身就是真实的动量汇，均匀外场的反冲由
    // 壁面吸收。bg 形参保留（删它会动力路径），只是传 0。
    update_force_field(cfg, sp.pp, N, st.Rx, st.Ry, st.Rz, st.Fx, st.Fy, st.Fz,
                       st.sum_phix, st.sum_phiy, st.sum_phiz, 0.0, 0.0, 0.0,
                       st.fx, st.fy, st.fz);
    // 重力体力必须加在 update_force_field 之后、step_navier_stokes 之前。
    const ExternalField g = sp.ext;
    if (cfg.density_on && (g.gx != 0.0 || g.gy != 0.0 || g.gz != 0.0))
    {
        apply_gravity(cfg, g, st.inv_rho_x, st.inv_rho_y, st.inv_rho_z,
                      st.fx, st.fy, st.fz);
    }
}

int step_fluid(NS_Config cfg, FpdState& st, long step)
{
    return step_navier_stokes(cfg, st.vx, st.vy, st.vz, st.p,
                              st.fx, st.fy, st.fz,
                              st.inv_rho_x, st.inv_rho_y, st.inv_rho_z,
                              st.cg_r, st.cg_s, st.cg_zq,
                              st.eta, st.etaXY, st.etaYZ, st.etaZX,
                              st.pi_dx, st.pi_dy, st.pi_dz, st.pi_nx, st.pi_ny, st.pi_nz,
                              st.fft, st.plan_xy, st.tri_w, st.diag,
                              st.gen, st.randD, st.randN,
                              st.tmp_fx, st.tmp_fy, st.tmp_fz, step);
}

void step_particles(const StepParams& sp, FpdState& st)
{
    const NS_Config cfg = sp.cfg;
    const int N = st.N;
    update_particle_velocity(cfg, sp.pp, N, st.Rx, st.Ry, st.Rz,
                             st.sum_phix, st.sum_phiy, st.sum_phiz,
                             st.vx, st.vy, st.vz, st.Vx, st.Vy, st.Vz);
    update_particle_position(cfg, N, st.Rx, st.Ry, st.Rz,
                             st.Rux, st.Ruy, st.Ruz, st.Vx, st.Vy, st.Vz);
}
