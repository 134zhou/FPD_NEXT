// spike/spike_state_map_kernel.cpp —— 单独的编译单元，模拟「fpd_core 里的 kernel」
//
// 这里只有 kernel，用 #pragma acc parallel loop present(...)。
// present 表是设备全局的，理论上与「谁建的映射」无关（无论是运行时 API 还是
// enter data pragma）。但 CLAUDE.md 记录 lambda 里的 update 会报
// "data in update host clause was not found on device"，所以本 spike 要实测：
// 用【运行时 API / 成员函数 pragma / 自由函数 pragma】三种方式建的映射，
// 能否被这个【另一个 TU】里的 present() 子句认账。

// 这个 kernel 由 spike_state_map.cpp 的 main 调用。它把 p[0..n) 每个元素 *2，
// 写入一个独立的输出数组 out[0..n)，然后 spike_state_map.cpp 把 out 拉回核对。
extern "C" void kernel_times2(double* p, double* out, int n)
{
    #pragma acc parallel loop present(p, out)
    for (int i = 0; i < n; i++)
    {
        out[i] = p[i] * 2.0;
    }
}
