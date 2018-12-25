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
#if defined(ENABLE_AVX512)
static constexpr int ALIGN = alignof(__m512);
#elif defined(ENABLE_AVX2)
static constexpr int ALIGN = alignof(__m256);
#elif defined(ENABLE_SSE)
static constexpr int ALIGN = alignof(__m128);
#elif defined(ENABLE_NEON)
static constexpr int ALIGN = alignof(float32x4_t);
#else
static constexpr int ALIGN = 8;
#endif  // defined(ENABLE_AVX512)

void print_m256(__m256 x) {
  printf("%f %f %f %f %f %f %f %f\n",
  x[7], x[6], x[5], x[4],
  x[3], x[2], x[1], x[0]);
}
void print_m128(__m128 x) {
  printf("%f %f %f %f\n",
  x[3], x[2], x[1], x[0]);
}


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


/*
  size_array  データ長
  gc          パディング、アラインメントが崩れた場合の想定

  $ ./test 8192 2
*/
int main(int argc, char *argv[])
{


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
/*
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
  }
*/
  // ===========================

  // ALIGN_SIZEでアラインして確保
  float* x = (float*)_mm_malloc( (32) * sizeof(float), ALIGN_SIZE);


  // 全範囲を初期化
  for (unsigned i=0; i<32; i++) {
    x[i] = (float)i;
  }

  // a = (7, 6, 5, 4, 3, 2, 1, 0) LSB
  __m256 a = _mm256_load_ps(&x[0]);
  print_m256(a);

  // t0 = (6, 5, 4, 7, 2, 1, 0, 3) LSB
  __m256 t0 = _mm256_permute_ps(a, _MM_SHUFFLE(2, 1, 0, 3));
  print_m256(t0);

  // t1 = (2, 1, 0, 3, 0, 0, 0, 0) LSB
  __m256 t1 = _mm256_permute2f128_ps(t0, t0, 0x08);
  print_m256(t1);

  // y  = (6, 5, 4, 3, 2, 1, 0, 0) LSB
  __m256 y = _mm256_blend_ps(t0, t1, 0x11);
  print_m256(y);

  // z = (0, 6, 5, 4, 3, 2, 1, 7)
  __m256i ofst = _mm256_set_epi32(0,6,5,4,3,2,1,7);
  __m256 z = _mm256_permutevar8x32_ps(a, ofst);
  print_m256(z);

  // 7.0
  __m256 ac = _mm256_broadcastss_ps( _mm256_extractf128_ps(z, 0) );
  print_m256(ac);

  // 1.0
  ac = _mm256_broadcastss_ps(
         _mm256_extractf128_ps(
           _mm256_permute_ps(z, _MM_SHUFFLE(3,2,1,1)), 0) );
  print_m256(ac);

  // 2.0
  ac = _mm256_broadcastss_ps(
         _mm256_extractf128_ps(
           _mm256_permute_ps(z, _MM_SHUFFLE(3,2,1,2)), 0) );
  print_m256(ac);

  // 3.0
  ac = _mm256_broadcastss_ps(
         _mm256_extractf128_ps(
           _mm256_permute_ps(z, _MM_SHUFFLE(3,2,1,3)), 0) );
  print_m256(ac);

  // 4.0
  ac = _mm256_broadcastss_ps(
         _mm256_extractf128_ps(
           _mm256_permutevar8x32_ps(z, _mm256_set_epi32(7,6,5,4,3,2,1,4)),
         0) );
  print_m256(ac);

  //__m256 b = _mm256_shuffle_ps(a, a, _MM_SHUFFLE(1,0,2,0));
  //__m256 c = _mm256_broadcastss_ps(b);
  //print_m256(c);






  _mm_free(x);


  // ===================
/*
  FILE* fp;
  fp=fopen("prof.txt", "w");

  char str[100];
  sprintf(str, "SIMD test");
  PM.print(stdout, "hoge", str);
  PM.print(fp, "hoge", str);
  PM.printThreads(fp, 1, 0); // time cost order
  fclose(fp);
*/
  return 0;
}
