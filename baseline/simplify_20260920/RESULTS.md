# 配置精简验证结果

基线：82cfd69；RTX 3060；GPU 测试显式设置 ACC_DEVICE_TYPE=nvidia。

| 判据 | 结果 |
|---|---|
| --check | 15 PASS；前后输出逐字节相同 |
| --check-potential | 14 PASS；前后输出逐字节相同 |
| --check-poisson | 4 PASS；前后输出逐字节相同 |
| --check-wall | 31 PASS；前后输出逐字节相同 |
| --check-config | 26 PASS，纯 CPU |
| bulk 单粒子 | 5 个检查点前后整文件逐位相同；200 步与 100+100 重启整文件逐位相同 |
| wall 单粒子 | 5 个检查点前后整文件逐位相同；200 步与 100+100 重启整文件逐位相同 |

bulk 初始中心 z=16.5、无壁面势；wall 初始中心 z=4.7、壁面 WCA，
表面间隙 h=2.0 位于力作用区内，gravity_z=-10。两者网格 32³、dt=0.002、seed=20260920。

Config.cpp：674 → 405 行；AGENTS.md：319 → 114 行；新增 CheckConfig.cpp：99 行。
保留基础/势参数必填、基本解析、运行中穿墙、NaN/Inf、初态一致性、文件完整性及 GPU 错误检查。
删除范围/组合/重复项限制、拼写建议、配置预测和 h_rest/z_rest 输出。

执行：`python3 baseline/simplify_20260920/validate.py before`（冻结程序）；
`python3 baseline/simplify_20260920/validate.py after`（当前 build）。
日志与配置已归档；冻结二进制及 .fpd 留在本机，由 .gitignore 排除；before.sha256 记录冻结产物校验和。
本次未修改二进制格式或 Python IO，也未运行可视化转换。
