/*
###################################################################################
#
# CubeZ
#
# Copyright (C) 2018 Research Institute for Information Technology(RIIT), Kyushu University.
# All rights reserved.
#
###################################################################################
*/

// @ref https://daily.belltail.jp/?p=1520
//      http://koturn.hatenablog.com/entry/2016/07/18/090000

#include "cz.h"
#include <xmmintrin.h>
#include <immintrin.h>

float dot_normal(const float *vec1, const float *vec2, unsigned n)
{
    float sum = 0;
    for(unsigned i = 0; i < n; ++i)
        sum += vec1[i] * vec2[i];
    return sum;
}

float dot_sse(const float *vec1, const float *vec2, unsigned n)
{
    __m128 u = {0};
    for (unsigned i = 0; i < n; i += 4)
    {
        __m128 w = _mm_load_ps(&vec1[i]);
        __m128 x = _mm_load_ps(&vec2[i]);

        x = _mm_mul_ps(w, x);
        u = _mm_add_ps(u, x);
    }
    __attribute__((aligned(16))) float t[4] = {0};
    _mm_store_ps(t, u);
    return t[0] + t[1] + t[2] + t[3];
}

float dot_avx(const float *vec1, const float *vec2, unsigned n)
{
    __m256 u = {0};
    for(unsigned i = 0; i < n; i += 8)
    {
        __m256 w = _mm256_load_ps(&vec1[i]);
        __m256 x = _mm256_load_ps(&vec2[i]);

        x = _mm256_mul_ps(w, x);
        u = _mm256_add_ps(u, x);
    }
    __attribute__((aligned(32))) float t[8] = {0};
    _mm256_store_ps(t, u);
    return t[0] + t[1] + t[2] + t[3] + t[4] + t[5] + t[6] + t[7];
}

float dot_avx_2(const float *vec1, const float *vec2, unsigned n)
{
    __m256 u1 = {0};
    __m256 u2 = {0};
    for(unsigned i = 0; i < n; i += 16)
    {
        __m256 w1 = _mm256_load_ps(&vec1[i]);
        __m256 w2 = _mm256_load_ps(&vec1[i + 8]);
        __m256 x1 = _mm256_load_ps(&vec2[i]);
        __m256 x2 = _mm256_load_ps(&vec2[i + 8]);

        x1 = _mm256_mul_ps(w1, x1);
        x2 = _mm256_mul_ps(w2, x2);
        u1 = _mm256_add_ps(u1, x1);
        u2 = _mm256_add_ps(u2, x2);
    }
    u1 = _mm256_add_ps(u1, u2);

    __attribute__((aligned(32))) float t[8] = {0};
    _mm256_store_ps(t, u1);
    return t[0] + t[1] + t[2] + t[3] + t[4] + t[5] + t[6] + t[7];
}

//FMA版
float dot_avx_fma(const float *vec1, const float *vec2, unsigned n)
{
    __m256 u1 = {0};
    __m256 u2 = {0};
    for(unsigned i = 0; i < n; i += 16)
    {
        __m256 w1 = _mm256_load_ps(&vec1[i]);
        __m256 w2 = _mm256_load_ps(&vec1[i + 8]);
        __m256 x1 = _mm256_load_ps(&vec2[i]);
        __m256 x2 = _mm256_load_ps(&vec2[i + 8]);
        //FMA命令で加算と乗算を行うけど，Haswellアーキテクチャ待ち(´・ω・｀)
        u1 = _mm256_fmadd_ps(w1, x1, u1);
        u2 = _mm256_fmadd_ps(w2, x2, u2);
    }
    u1 = _mm256_add_ps(u1, u2);
    //レジスタから書き戻し
    __attribute__((aligned(32))) float t[8] = {0};
    _mm256_store_ps(t, u1);
    return t[0] + t[1] + t[2] + t[3] + t[4] + t[5] + t[6] + t[7];
}
