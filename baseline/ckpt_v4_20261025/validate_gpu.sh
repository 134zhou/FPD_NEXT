#!/usr/bin/env bash
# 从仓库根目录执行；当前二进制构建完成后再运行。
set -eu
export ACC_DEVICE_TYPE=nvidia
archive=baseline/ckpt_v4_20261025
for case in wall pair_disjoint; do
    ./build/fpd "$archive/$case.cfg" --set "init_file=$archive/${case}_v4.fpd" --set "out_dir=$archive/after" > "$archive/after_$case.log" 2>&1
    echo "PASS $case 10 步"
done
./build/fpd "$archive/wall.cfg" --set "init_file=$archive/after/wall_0000005.fpd" --set "out_dir=$archive/restart" > "$archive/restart.log" 2>&1
cmp "$archive/after/wall_0000010.fpd" "$archive/restart/wall_0000010.fpd"
echo 'PASS 单粒子 10 vs 5+5 整文件逐位相同'
for check in check-config check check-potential check-wall pcg-noise; do
    ./build/fpd_check "--$check" > "$archive/$check.log" 2>&1
    echo "PASS --$check"
done
