#include <stdio.h>
#include <stdlib.h>
#include <cmath>

// std::inner_product, accumulate
#include <numeric>

#include <xmmintrin.h>
#include <immintrin.h>

#include <PerfMonitor.h>

#define PM_NUM_MAX 50
#define TM_LABEL_MAX 50
#define ALIGN_SIZE 64
#define NN 500

int order_of_PM_key=0;

pm_lib::PerfMonitor PM;

// c11のコンパイルオプションが必要
#ifndef _NOAVX512
static constexpr int ALIGN = alignof(__m512);
#endif

/*
#if defined(ENABLE_AVX512)
static constexpr int ALIGN = alignof(__m512);
#elif defined(ENABLE_AVX)
static constexpr int ALIGN = alignof(__m256);
#elif defined(ENABLE_SSE)
static constexpr int ALIGN = alignof(__m128);
#elif defined(ENABLE_NEON)
static constexpr int ALIGN = alignof(float32x4_t);
#else
static constexpr int ALIGN = 8;
#endif  // defined(ENABLE_AVX512)
*/

void set_label(const std::string label, pm_lib::PerfMonitor::Type type)
{
  bool exclusive=true;

  // 登録個数のチェック
  order_of_PM_key++;

  if ( order_of_PM_key > PM_NUM_MAX )
  {
    printf("\tThe number of labels for Performance monitor goes over limit.\n");
    exit(0);
  }

  // 文字数がTM_LABEL_MAX-1を超えるものはカット
  if ( strlen(label.c_str()) > TM_LABEL_MAX-1 )
  {
    printf("\tWarning: Length of timing label must be less than %d\n", TM_LABEL_MAX-1);
  }

  // Performance Monitorへの登録
  PM.setProperties(label, type, exclusive);
}

inline void TIMING_start(const std::string key) {
  PM.start(key);
}

inline void TIMING_stop(const std::string key, double flopPerTask=0.0, int iterationCount=1) {
  PM.stop(key, flopPerTask, (unsigned)iterationCount);
}



float dot_1_normal(const unsigned n, const float *vec1, double& flop)
{
  float sum = 0.0f;
  flop += (double)(n*1.0);
  float e;

  #pragma omp parallel for reduction(+:sum) private(e)
  for(unsigned i = 0; i < n; i++)
    e = vec1[i];
    sum += e * e;
  return sum;
}

float dot_2_normal(const unsigned n, const float *vec1, const float *vec2, double& flop)
{
  float sum = 0.0f;
  flop += (double)(n*1.0);

  #pragma omp parallel for reduction(+:sum)
  for(unsigned i = 0; i < n; i++)
    sum += vec1[i] * vec2[i];
  return sum;
}

float dot_2_fma(const unsigned n, const float *vec1, const float *vec2, double& flop)
{
  float sum = 0.0f;
  flop += (double)(n*1.0);

  for (std::size_t i = 0; i < n; i++) {
    // <cmath>のstd::fma関数を用いると，積和演算がハードウェアのサポートを受けることを期待できる
    sum = std::fma(vec1[i], vec2[i], sum);
  }
  return sum;
}


float dot_1_sse(const unsigned n, const float *vec1, double& flop)
{
  flop += (double)(n*1.0);

  __m128 u = {0};
  for (unsigned i = 0; i < n; i += 4)
  {
      __m128 w = _mm_load_ps(&vec1[i]);
      __m128 x = _mm_load_ps(&vec1[i]);

      x = _mm_mul_ps(w, x);
      u = _mm_add_ps(u, x);
  }
  __attribute__((aligned(16))) float t[4] = {0};
  _mm_store_ps(t, u);
  return t[0] + t[1] + t[2] + t[3];
}

float dot_2_sse(const unsigned n, const float *vec1, const float *vec2, double& flop)
{
  flop += (double)(n*1.0);

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

float dot_2_avx(const unsigned n, const float *vec1, const float *vec2, double& flop)
{
  flop += (double)(n*1.0);

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

float dot_2_avx_u(const unsigned n, const float *vec1, const float *vec2, double& flop)
{
  flop += (double)(n*1.0);

  __m256 u = {0};
  for(unsigned i = 0; i < n; i += 8)
  {
      __m256 w = _mm256_loadu_ps(&vec1[i]);
      __m256 x = _mm256_loadu_ps(&vec2[i]);

      x = _mm256_mul_ps(w, x);
      u = _mm256_add_ps(u, x);
  }
  __attribute__((aligned(32))) float t[8] = {0};
  _mm256_storeu_ps(t, u);
  return t[0] + t[1] + t[2] + t[3] + t[4] + t[5] + t[6] + t[7];
}

float dot_2_avx2(const unsigned n, const float *vec1, const float *vec2, double& flop)
{
  flop += (double)(n*1.0);

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
float dot_2_avx_fma(const unsigned n, const float *vec1, const float *vec2, double& flop)
{
  flop += (double)(n*1.0);

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


#ifndef _NOAVX512
float dot_2_avx512(const unsigned n, const float *vec1, const float *vec2, double& flop)
{
  static constexpr std::size_t INTERVAL = sizeof(__m512) / sizeof(float);

  __m512 sumx16 = {0};

  for (std::size_t i = 0; i < n; i += INTERVAL) {
    __m512 ax16 = _mm512_load_ps(&vec1[i]);
    __m512 bx16 = _mm512_load_ps(&vec2[i]);
    sumx16 = _mm512_add_ps(sumx16, _mm512_mul_ps(ax16, bx16));
  }

  alignas(ALIGN) float s[INTERVAL] = {0};
  _mm512_store_ps(s, sumx16);

  std::size_t offset = n - n % INTERVAL;
  return std::inner_product(
      vec1 + offset,
      vec1 + n,
      vec2 + offset,
      std::accumulate(std::begin(s), std::end(s), 0.0f));
}

float dot_2_avx512_fma(const unsigned n, const float *vec1, const float *vec2, double& flop)
{
  static constexpr std::size_t INTERVAL = sizeof(__m512) / sizeof(float);

  __m512 sumx16 = {0};

  for (std::size_t i = 0; i < n; i += INTERVAL) {
    __m512 ax16 = _mm512_load_ps(&vec1[i]);
    __m512 bx16 = _mm512_load_ps(&vec2[i]);
    sumx16 = _mm512_fmadd_ps(ax16, bx16, sumx16);
  }

  alignas(ALIGN) float s[INTERVAL] = {0};
  _mm512_store_ps(s, sumx16);

  std::size_t offset = n - n % INTERVAL;
  return std::inner_product(
      vec1 + offset,
      vec1 + n,
      vec2 + offset,
      std::accumulate(std::begin(s), std::end(s), 0.0f));
}
#endif

/*
  size_array  データ長
  gc          パディング、アラインメントが崩れた場合の想定

  $ ./test 8192 2
*/
int main(int argc, char *argv[])
{

  if (argc != 3) {
    printf("Usage : ./test size_array gc\n");
    exit(0);
  }
  const int n = atoi(argv[1]);
  const int gc= atoi(argv[2]);


  //const int n = SZ*SZ*SZ;
  printf("size = %d   gc=%d\n", n, gc);

  int numThreads=0;
  double flop = 0.0;

#ifdef _OPENMP
  char* c_env = std::getenv("OMP_NUM_THREADS");
  if (c_env == NULL) {
    omp_set_num_threads(1);	// OMP_NUM_THREADS was not defined. set as 1
  }
  numThreads  = omp_get_max_threads();
#else
  numThreads  = 1;
#endif

  PM.initialize( PM_NUM_MAX );
  PM.setRankInfo( 1 );
  std::string Parallel_str;
  if ( numThreads > 1 ) {
    Parallel_str = "OpenMP";
  }
  else {
    Parallel_str = "Serial";
  }
  PM.setParallelMode(Parallel_str, numThreads, 1);

  {
    using namespace pm_lib;

    set_label("Dot1_normal",  PerfMonitor::CALC);
    set_label("Dot1_SSE" ,    PerfMonitor::CALC);
    set_label("Dot1_AVX" ,    PerfMonitor::CALC);
    set_label("Dot1_AVX2" ,   PerfMonitor::CALC);
    set_label("Dot1_AVX_FMA", PerfMonitor::CALC);

    set_label("Dot2_normal",  PerfMonitor::CALC);
    set_label("Dot2_fma",     PerfMonitor::CALC);
    set_label("Dot2_SSE" ,    PerfMonitor::CALC);
    set_label("Dot2_AVX" ,    PerfMonitor::CALC);
    set_label("Dot2_AVX_U" ,  PerfMonitor::CALC);
    set_label("Dot2_AVX2" ,   PerfMonitor::CALC);
    set_label("Dot2_AVX_FMA", PerfMonitor::CALC);
    set_label("Dot2_AVX512" ,   PerfMonitor::CALC);
    set_label("Dot2_AVX512_FMA", PerfMonitor::CALC);
  }

  // ===========================

  // ALIGN_SIZEでアラインして確保
  float* x = (float*)_mm_malloc( (n+2*gc) * sizeof(float), ALIGN_SIZE);
  float* y = (float*)_mm_malloc( (n+2*gc) * sizeof(float), ALIGN_SIZE);
  float* z = (float*)_mm_malloc( (n+2*gc) * sizeof(float), ALIGN_SIZE);

  // アラインしているか確認
  long int li = (long int)x;
  printf("x   : addrs= %08x align(4=%2u, 8=%2u, 16=%2u, 32=%2d, 64=%2d)\n\n",
         li, li%4, li%8, li%16, li%32, li%64);
  li = (long int)y;
  printf("y   : addrs= %08x align(4=%2u, 8=%2u, 16=%2u, 32=%2d, 64=%2d)\n\n",
                li, li%4, li%8, li%16, li%32, li%64);


  // 全範囲を初期化
  for (unsigned i=0; i<n+2*gc; i++) {
    x[i] = 0.0f;
    y[i] = 0.0f;
    z[i] = 0.0f;
  }

  // ループ対象の該当部分だけに値を入れる
  // (x,y)で総和になる
  for (unsigned i=gc; i<n+gc; i++) {
    x[i] = (float)(i-gc+1);
  }
  for (unsigned i=gc; i<n+gc; i++) {
    y[i] = 1.0f;
  }


  float sum = 0.0, sum2=0.0;
  for (unsigned i=gc; i<n+gc; i++) sum += (float)(i-gc+1);





  for (int m=0; m<NN; m++)
  {
    flop = 0.0;
    TIMING_start("Dot1_normal");
    sum2 = dot_1_normal(n, &x[gc], flop);
    TIMING_stop("Dot1_normal", flop);
  }
  printf("Dot1 normal\t= %10.1f %10.1f %10.1f\n", 0, 0, 0);


  for (int m=0; m<NN; m++)
  {
    flop = 0.0;
    TIMING_start("Dot2_normal");
    sum2 = dot_2_normal(n, &x[gc], &y[gc], flop);
    TIMING_stop("Dot2_normal", flop);
  }
  printf("Dot2 normal\t= %10.1f %10.1f %10.1f\n", sum2, sum, sum2-sum);

/*
  for (int m=0; m<NN; m++)
  {
    flop = 0.0;
    TIMING_start("Dot2_fma");
    sum2 = dot_2_fma(n, x, y, flop);
    TIMING_stop("Dot2_fma", flop);
  }
  printf("Dot2 fma\t= %10.1f %10.1f %10.1f\n", sum2, sum, sum2-sum);


  for (int m=0; m<NN; m++)
  {
    flop = 0.0;
    TIMING_start("Dot2_SSE");
    sum2 = dot_2_sse(n, &x[gc], &y[gc], flop);
    TIMING_stop("Dot2_SSE", flop);
  }
  printf("Dot2 SSE\t= %10.1f %10.1f %10.1f\n", sum2, sum, sum2-sum);
*/

  for (int m=0; m<NN; m++)
  {
    flop = 0.0;
    TIMING_start("Dot2_AVX");
    sum2 = dot_2_avx(n, &x[gc], &y[gc], flop);
    TIMING_stop("Dot2_AVX", flop);
  }
  printf("Dot2 AVX\t= %10.1f %10.1f %10.1f\n", sum2, sum, sum2-sum);

  for (int m=0; m<NN; m++)
  {
    flop = 0.0;
    TIMING_start("Dot2_AVX_U");
    sum2 = dot_2_avx_u(n, &x[gc], &y[gc], flop);
    TIMING_stop("Dot2_AVX_U", flop);
  }
  printf("Dot2 AVX_U\t= %10.1f %10.1f %10.1f\n", sum2, sum, sum2-sum);


  for (int m=0; m<NN; m++)
  {
    flop = 0.0;
    TIMING_start("Dot2_AVX2");
    sum2 = dot_2_avx2(n, &x[gc], &y[gc], flop);
    TIMING_stop("Dot2_AVX2", flop);
  }
  printf("Dot2 AVX2\t= %10.1f %10.1f %10.1f\n", sum2, sum, sum2-sum);


  for (int m=0; m<NN; m++)
  {
    flop = 0.0;
    TIMING_start("Dot2_AVX_FMA");
    sum2 = dot_2_avx_fma(n, &x[gc], &y[gc], flop);
    TIMING_stop("Dot2_AVX_FMA", flop);
  }
  printf("Dot2 AVX_FMA\t= %10.1f %10.1f %10.1f\n", sum2, sum, sum2-sum);

#ifndef _NOAVX512
  for (int m=0; m<NN; m++)
  {
    flop = 0.0;
    TIMING_start("Dot2_AVX512");
    sum2 = dot_2_avx512(n, &x[gc], &y[gc], flop);
    TIMING_stop("Dot2_AVX512", flop);
  }
  printf("Dot2 AVX512\t= %10.1f %10.1f %10.1f\n", sum2, sum, sum2-sum);


  for (int m=0; m<NN; m++)
  {
    flop = 0.0;
    TIMING_start("Dot2_AVX512_FMA");
    sum2 = dot_2_avx512_fma(n, &x[gc], &y[gc], flop);
    TIMING_stop("Dot2_AVX512_FMA", flop);
  }
  printf("Dot2 AVX512_FMA\t= %10.1f %10.1f %10.1f\n", sum2, sum, sum2-sum);
#endif

  _mm_free(x);
  _mm_free(y);
  _mm_free(z);


  // ===================

  char str[100];
  sprintf(str, "SIMD test");
  PM.print(stdout, "hoge", str);
  //PM.printThreads(stdout, 1, 0); // time cost order

  return 0;
}
