# .fpd v3 与单一参数来源

基线提交 eec6ef8，工作区干净。始终要求 init_file；头部只存 magic/endian/version、Nx/Ny/Nz/N、step；删除 flags、压力及物理/RNG 元数据，尾部仍 FNV。配置或 --set 出现任一尺寸键时报错，其余参数均读配置。生产单次打开读取；离线工具显式 --config；v2 明确拒绝，不迁移。

先冻结 v2 数值基线，然后改 C++/Python 格式、配置和工具，最后同步文档并验证 GPU 轨迹、断点重启、互操作和工具。旧版二进制、测试产物不入库，日志与结果入 baseline/。
