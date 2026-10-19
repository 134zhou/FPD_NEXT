#!/bin/bash
# S4：sed.cfg（128x64x32，N=95，ratio_rho=2，噪声开）同 seed 不同容差 / 分裂法 5000 步
X=/tmp/fpdb2; R=/home/doll/Documents/C++/FPD_NEXT; O=/tmp/fpds4/out
export ACC_DEVICE_TYPE=nvidia
cd $R
common="--set init_file=/tmp/fpddoc/sed_init.fpd --set out_dir=$O --set n_steps=5000 --set interval_ckpt=5000 --set interval_log=1000"
for t in 1e-14 1e-10 1e-8 1e-6 1e-4 1e-3; do
  $X/fpd config/sed.cfg $common --set run_name=tol$t --set pcg_tol=$t > /tmp/fpds4/tol$t.log 2>&1
  echo "tol=$t rc=$? $(grep 'step      4000' /tmp/fpds4/tol$t.log)"
done
$X/fpd config/sed.cfg $common --set run_name=seedB --set pcg_tol=1e-14 --set seed=777 > /tmp/fpds4/seedB.log 2>&1
echo "seedB rc=$? $(grep 'step      4000' /tmp/fpds4/seedB.log)"
$X/fpd config/sed.cfg $common --set run_name=split --set pressure_solver=split > /tmp/fpds4/split.log 2>&1
echo "split rc=$? $(grep 'step      4000' /tmp/fpds4/split.log)"
python3 /tmp/fpds4/cmp.py $O/tol1e-14_0005000.fpd $O/tol1e-10_0005000.fpd $O/tol1e-8_0005000.fpd \
   $O/tol1e-6_0005000.fpd $O/tol1e-4_0005000.fpd $O/tol1e-3_0005000.fpd $O/split_0005000.fpd $O/seedB_0005000.fpd

# 确定性对照（E2 构型：单粒子 32³，无势、noise_on=0，g_z=-0.02），5000 步
cat > /tmp/fpds4/e2.cfg <<'EOC'
dt = 0.002
n_steps = 5000
kT = 0.25
radius = 3.2
xi = 1.0
ratio_eta = 50.0
ratio_rho = 2.0
noise_on = 0
gravity_z = -0.02
init_file = /tmp/fpdw/n1.fpd
out_dir = /tmp/fpds4/out
interval_ckpt = 5000
interval_log = 1000
seed = 1
EOC
for t in 1e-14 1e-10 1e-8 1e-6 1e-4; do
  $X/fpd /tmp/fpds4/e2.cfg --set run_name=e2tol$t --set pcg_tol=$t > /tmp/fpds4/e2tol$t.log 2>&1
done
$X/fpd /tmp/fpds4/e2.cfg --set run_name=e2split --set pressure_solver=split > /tmp/fpds4/e2split.log 2>&1
python3 /tmp/fpds4/cmp.py $O/e2tol1e-14_0005000.fpd $O/e2tol1e-10_0005000.fpd $O/e2tol1e-8_0005000.fpd \
   $O/e2tol1e-6_0005000.fpd $O/e2tol1e-4_0005000.fpd $O/e2split_0005000.fpd
for f in e2tol1e-14 e2tol1e-6 e2split; do echo "$f $(grep 'step      4000' /tmp/fpds4/$f.log)"; done
