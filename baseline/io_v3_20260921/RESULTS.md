# `.fpd` v3 验证结果

基线 eec6ef8（v2），新格式 v3。RTX 3060，GPU 运行使用 ACC_DEVICE_TYPE=nvidia。

| 判据 | 结果 |
|---|---|
| --check / --check-potential / --check-poisson / --check-wall | 15 / 14 / 4 / 31 PASS；与原基线日志逐字节相同 |
| --check-config | 32 PASS，含 init_file 必填及三维尺寸键在配置和 --set 的拒绝 |
| v2/v3 载荷 | bulk、wall 的 step 0/50 速度与粒子数组逐位相同；v2 的压力数组从 v3 删除 |
| v3 重启 | bulk、wall 的连续 50 步与 25+25 步，整文件逐位相同 |
| 互操作 | Python v3 文件由 C++ 成功读出索引 (3,1,1)=3001001；C++/Python 均明确拒绝 v2；翻转数据字节后校验和拒绝 |
| 力及统计 | 近壁 WCA 的 --verify-forces PASS；--sed-stats 从配置得到 dt/radius/xi |
| VTK | self-test 与 v3→VTK→ParaView 往返通过；时间 t=0.1，速度和粒子数组差 0 |

复查数值载荷：`python3 baseline/io_v3_20260921/validate_payload.py`。
冻结 v2 二进制与检查点、v3 检查点留在本机且被本目录 .gitignore 排除；配置与日志入库。

最终增量核验：生产路径停止回传压力后，--check-wall 日志逐字节相同；bulk/wall 在
step 0/50 的 v3 检查点也与优化前整文件逐位相同。`--check-config` 再次全绿。
