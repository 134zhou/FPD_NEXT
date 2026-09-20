#!/usr/bin/env python3
"""配置精简的同机 GPU 基线与回归；从仓库根目录运行。"""
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
BASE = Path(__file__).resolve().parent
os.chdir(ROOT)
sys.path.insert(0, str(ROOT / 'tools'))
import fpd_format as F

phase = sys.argv[1]
assert phase in ('before', 'after')
dest = BASE / phase
bindir = dest if phase == 'before' else ROOT / 'build'
env = dict(os.environ, ACC_DEVICE_TYPE='nvidia')

def run(args, log):
    print('运行', ' '.join(map(str, args)), flush=True)
    with (dest / log).open('w') as out:
        subprocess.run(list(map(str, args)), stdout=out, stderr=subprocess.STDOUT,
                       env=env, check=True)

for flag in ('check', 'check-potential', 'check-poisson', 'check-wall'):
    run([bindir / 'fpd_check', '--' + flag], flag + '.log')

for name, z in [('bulk', 16.5), ('wall', 4.7)]:
    initial = BASE / (name + '_init.fpd')
    if phase == 'before':
        header = F.default_header(32, 32, 32, 1, seed=20260920)
        particles = dict(Rx=[16.3], Ry=[16.7], Rz=[z], Rux=[16.3], Ruy=[16.7], Ruz=[z])
        F.write_fpd(initial, header, dict(vx=None, vy=None, vz=None), particles)
        config = f'''Nx = 32
Ny = 32
Nz = 32
dt = 0.002
n_steps = 200
kT = 0.25
radius = 3.2
xi = 1
ratio_eta = 50
seed = 20260920
init_file = {initial.relative_to(ROOT)}
interval_ckpt = 50
interval_log = 50
save_pressure = 1
'''
        if name == 'wall':
            config += 'wallpotential = wca\nwall_eps = 1\nwall_sigma = 2\ngravity_z = -10\n'
        (BASE / (name + '.cfg')).write_text(config)
    cfg = BASE / (name + '.cfg')
    args = [bindir / 'fpd', cfg, '--set', 'out_dir=' + str(dest), '--set', 'run_name=' + name]
    run(args, name + '.log')
    if phase == 'after':
        for step in (0, 50, 100, 150, 200):
            filename = f'{name}_{step:07d}.fpd'
            run([bindir / 'fpd_tool', '--diff-ckpt', BASE / 'before' / filename, dest / filename],
                f'{name}_diff_{step}.log')
        run(args + ['--set', 'run_name=' + name + '_split', '--set', 'n_steps=100'], name + '_split.log')
        run(args + ['--set', 'run_name=' + name + '_restart', '--set',
                    'init_file=' + str(dest / f'{name}_split_0000100.fpd')], name + '_restart.log')
        run([bindir / 'fpd_tool', '--diff-ckpt', dest / f'{name}_0000200.fpd',
             dest / f'{name}_restart_0000200.fpd'], name + '_restart_diff.log')
print(phase, '全部通过', flush=True)
