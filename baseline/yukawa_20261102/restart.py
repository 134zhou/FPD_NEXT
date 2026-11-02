import sys, pathlib, subprocess, os
sys.path.insert(0, str(pathlib.Path('tools').resolve()))
from fpd_format import default_header, write_fpd
base=pathlib.Path('baseline/yukawa_20261102')
for n in (1,2):
    init=base/f'init{n}.fpd'
    x=[8.,24.][:n]; y=[16.]*n; z=[16.]*n
    write_fpd(init,default_header(32,32,32,n),{},dict(Rx=x,Ry=y,Rz=z,Rux=x,Ruy=y,Ruz=z))
    for mode in ('direct','cell'):
        cfg=base/f'{mode}{n}.cfg'
        cfg.write_text(f'''dt=0.002
n_steps=6
kT=0.01
radius=1
xi=0.5
ratio_eta=2
ratio_rho=1
init_file={init}
out_dir={base}/{mode}{n}
run_name=run
interval_ckpt=3
interval_log=3
seed=123
potential=yukawa
pot_eps=1
pot_sigma=2
pot_kappa=0.2
pot_rcut=20
pair_search={mode}
''')
        env=dict(os.environ,ACC_DEVICE_TYPE='nvidia')
        subprocess.run(['./build/fpd',str(cfg)],env=env,check=True)
        checkpoint=base/f'{mode}{n}/run_0000003.fpd'
        out=base/f'{mode}{n}_restart'
        subprocess.run(['./build/fpd',str(cfg),'--set',f'init_file={checkpoint}','--set',f'out_dir={out}'],env=env,check=True)
        assert (base/f'{mode}{n}/run_0000006.fpd').read_bytes()==(out/'run_0000006.fpd').read_bytes()
        print('RESTART PASS',mode,n,flush=True)
