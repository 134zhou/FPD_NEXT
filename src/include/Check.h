#ifndef CHECK_H
#define CHECK_H

#include <cstdio>
#include <cstdlib>
#include <cufft.h>
#include <curand.h>

// cuFFT / cuRAND 的返回码检查（修 C11：原先全部被忽略）。
// 失败即退出 —— 这类错误没有合理的恢复路径，继续跑只会产生垃圾数据。

#define CUFFT_CHECK(x)                                                        \
    do {                                                                      \
        cufftResult _r = (x);                                                 \
        if (_r != CUFFT_SUCCESS) {                                            \
            fprintf(stderr, "cuFFT 错误 %d  @ %s:%d\n  %s\n",                 \
                    (int)_r, __FILE__, __LINE__, #x);                         \
            exit(EXIT_FAILURE);                                               \
        }                                                                     \
    } while (0)

#define CURAND_CHECK(x)                                                       \
    do {                                                                      \
        curandStatus_t _r = (x);                                              \
        if (_r != CURAND_STATUS_SUCCESS) {                                    \
            fprintf(stderr, "cuRAND 错误 %d  @ %s:%d\n  %s\n",                \
                    (int)_r, __FILE__, __LINE__, #x);                         \
            exit(EXIT_FAILURE);                                               \
        }                                                                     \
    } while (0)

#endif
