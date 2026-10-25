#!/usr/bin/env bash
set -eu
export ACC_DEVICE_TYPE=nvidia
archive=baseline/ckpt_v4_20261025
/tmp/fpd-v3-before "$archive/pair_disjoint.cfg" > "$archive/before_pair_disjoint.log" 2>&1
./build/fpd "$archive/pair_disjoint.cfg" --set "init_file=$archive/pair_disjoint_v4.fpd" --set "out_dir=$archive/after" > "$archive/after_pair_disjoint.log" 2>&1
for version in v3 v4; do
    if [ "$version" = v3 ]; then binary=/tmp/fpd-v3-before; else binary=./build/fpd; fi
    "$binary" "$archive/wall.cfg" --set "init_file=$archive/wall_${version}.fpd" --set "out_dir=$archive/zero_$version" --set n_steps=0 > "$archive/zero_$version.log" 2>&1
    set +e
    "$binary" "$archive/wall.cfg" --set "init_file=$archive/wall_${version}.fpd" --set "out_dir=$archive/failed_$version" --set pcg_tol=1e-30 --set n_steps=1 > "$archive/failed_$version.log" 2>&1
    code=$?
    set -e
    test "$code" -eq 3
    echo "PASS $version PCG 未收敛返回 3"
done
