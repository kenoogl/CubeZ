/*
###################################################################################
#
# CubeZ
#
# Copyright (C) 2018-2020 Research Institute for Information Technology(RIIT), Kyushu University.
# All rights reserved.
#
###################################################################################
*/

#ifndef _CZ_F_FUNC_H_
#define _CZ_F_FUNC_H_


extern "C" {

// cz_solver.f90

void bc_k_    (int* sz,
               int* g,
               REAL_TYPE* e,
               REAL_TYPE* dh,
               REAL_TYPE* org,
               int* nID);

void jacobi_        (REAL_TYPE* p,
                     int* sz,
                     int* idx,
                     int* g,
                     REAL_TYPE* cf,
                     REAL_TYPE* omg,
                     REAL_TYPE* b,
                     double* res,
                     REAL_TYPE* wk2,
                     double* flop);

void psor_          (REAL_TYPE* p,
                     int* sz,
                     int* idx,
                     int* g,
                     REAL_TYPE* cf,
                     REAL_TYPE* omg,
                     REAL_TYPE* b,
                     double* res,
                     double* flop);

void psor2sma_core_ (REAL_TYPE* p,
                     int* sz,
                     int* idx,
                     int* g,
                     REAL_TYPE* cf,
                     int* ip,
                     int* color,
                     REAL_TYPE* omg,
                     REAL_TYPE* b,
                     double* res,
                     double* flop);

void pcr_rb_ (int* sz,
              int* idx,
              int* g,
              int* pn,
              int* ofst,
              int* color,
              REAL_TYPE* x,
              REAL_TYPE* msk,
              REAL_TYPE* rhs,
              REAL_TYPE* WA,
              REAL_TYPE* WC,
              REAL_TYPE* WD,
              REAL_TYPE* WAA,
              REAL_TYPE* WCC,
              REAL_TYPE* WDD,
              REAL_TYPE* omg,
              double* res,
              double* flop);

void pcr_rb_esa_ (int* sz,
                  int* idx,
                  int* g,
                  int* pn,
                  int* ofst,
                  int* color,
                  int* s,
                  REAL_TYPE* x,
                  REAL_TYPE* msk,
                  REAL_TYPE* rhs,
                  REAL_TYPE* WA,
                  REAL_TYPE* WC,
                  REAL_TYPE* WD,
                  REAL_TYPE* WAA,
                  REAL_TYPE* WCC,
                  REAL_TYPE* WDD,
                  REAL_TYPE* omg,
                  double* res,
                  double* flop);
  
void pcr_(int* sz,
          int* idx,
          int* g,
          int* pn,
          REAL_TYPE* x,
          REAL_TYPE* msk,
          REAL_TYPE* rhs,
          REAL_TYPE* WA,
          REAL_TYPE* WC,
          REAL_TYPE* WD,
          REAL_TYPE* WAA,
          REAL_TYPE* WCC,
          REAL_TYPE* WDD,
          REAL_TYPE* omg,
          double* res,
          double* flop);

void pcr_eda_(int* sz,
           int* idx,
           int* g,
           int* pn,
           REAL_TYPE* x,
           REAL_TYPE* msk,
           REAL_TYPE* rhs,
           REAL_TYPE* WA,
           REAL_TYPE* WC,
           REAL_TYPE* WD,
           REAL_TYPE* omg,
           double* res,
           double* flop);

void pcr_esa_(int* sz,
              int* idx,
              int* g,
              int* pn,
              int* s,
              REAL_TYPE* x,
              REAL_TYPE* msk,
              REAL_TYPE* rhs,
              REAL_TYPE* SA,
              REAL_TYPE* SC,
              REAL_TYPE* SD,
              REAL_TYPE* WA,
              REAL_TYPE* WC,
              REAL_TYPE* WD,
              REAL_TYPE* omg,
              double* res,
              double* flop);
  
void pcr_j_esa_(int* sz,
                int* idx,
                int* g,
                int* pn,
                int* s,
                REAL_TYPE* x,
                REAL_TYPE* msk,
                REAL_TYPE* rhs,
                REAL_TYPE* SA,
                REAL_TYPE* SC,
                REAL_TYPE* SD,
                REAL_TYPE* WA,
                REAL_TYPE* WC,
                REAL_TYPE* WD,
                REAL_TYPE* SRC,
                REAL_TYPE* WRK,
                REAL_TYPE* omg,
                double* res,
                double* flop);
  
  
// cz_maf.f90
void jacobi_maf_   (REAL_TYPE* p,
                    int* sz,
                    int* idx,
                    int* g,
                    REAL_TYPE* X,
                    REAL_TYPE* Y,
                    REAL_TYPE* Z,
                    REAL_TYPE* omg,
                    REAL_TYPE* b,
                    double* res,
                    REAL_TYPE* wk2,
                    REAL_TYPE* tmp,
                    double* flop);
    
void psor_maf_     (REAL_TYPE* p,
                    int* sz,
                    int* idx,
                    int* g,
                    REAL_TYPE* X,
                    REAL_TYPE* Y,
                    REAL_TYPE* Z,
                    REAL_TYPE* omg,
                    REAL_TYPE* b,
                    double* res,
                    double* flop);
    
void psor2sma_core_maf_ (REAL_TYPE* p,
                         int* sz,
                         int* idx,
                         int* g,
                         REAL_TYPE* X,
                         REAL_TYPE* Y,
                         REAL_TYPE* Z,
                         int* ip,
                         int* color,
                         REAL_TYPE* omg,
                         REAL_TYPE* b,
                         double* res,
                         REAL_TYPE* tmp,
                         double* flop);
    
void pcr_rb_maf_(int* sz,
                 int* idx,
                 int* g,
                 int* pn,
                 int* ofst,
                 int* color,
                 REAL_TYPE* x,
                 REAL_TYPE* msk,
                 REAL_TYPE* rhs,
                 REAL_TYPE* XX,
                 REAL_TYPE* YY,
                 REAL_TYPE* ZZ,
                 REAL_TYPE* WA,
                 REAL_TYPE* WC,
                 REAL_TYPE* WD,
                 REAL_TYPE* WAA,
                 REAL_TYPE* WCC,
                 REAL_TYPE* WDD,
                 REAL_TYPE* omg,
                 double* res,
                 REAL_TYPE* tmp,
                 double* flop);
  
void pcr_rb_esa_maf_(int* sz,
                     int* idx,
                     int* g,
                     int* pn,
                     int* ofst,
                     int* color,
                     int* s,
                     REAL_TYPE* x,
                     REAL_TYPE* msk,
                     REAL_TYPE* rhs,
                     REAL_TYPE* XX,
                     REAL_TYPE* YY,
                     REAL_TYPE* ZZ,
                     REAL_TYPE* WA,
                     REAL_TYPE* WC,
                     REAL_TYPE* WD,
                     REAL_TYPE* WAA,
                     REAL_TYPE* WCC,
                     REAL_TYPE* WDD,
                     REAL_TYPE* omg,
                     double* res,
                     REAL_TYPE* tmp,
                     double* flop);
  
void pcr_maf_(int* sz,
              int* idx,
              int* g,
              int* pn,
              REAL_TYPE* x,
              REAL_TYPE* msk,
              REAL_TYPE* rhs,
              REAL_TYPE* XX,
              REAL_TYPE* YY,
              REAL_TYPE* ZZ,
              REAL_TYPE* WA,
              REAL_TYPE* WC,
              REAL_TYPE* WD,
              REAL_TYPE* WAA,
              REAL_TYPE* WCC,
              REAL_TYPE* WDD,
              REAL_TYPE* omg,
              double* res,
              REAL_TYPE* tmp,
              double* flop);

void pcr_eda_maf_(int* sz,
               int* idx,
               int* g,
               int* pn,
               REAL_TYPE* x,
               REAL_TYPE* msk,
               REAL_TYPE* rhs,
               REAL_TYPE* XX,
               REAL_TYPE* YY,
               REAL_TYPE* ZZ,
               REAL_TYPE* WA,
               REAL_TYPE* WC,
               REAL_TYPE* WD,
               REAL_TYPE* omg,
               double* res,
               REAL_TYPE* tmp,
               double* flop);

void pcr_esa_maf_(int* sz,
                  int* idx,
                  int* g,
                  int* pn,
                  int* s,
                  REAL_TYPE* x,
                  REAL_TYPE* msk,
                  REAL_TYPE* rhs,
                  REAL_TYPE* XX,
                  REAL_TYPE* YY,
                  REAL_TYPE* ZZ,
                  REAL_TYPE* A,
                  REAL_TYPE* C,
                  REAL_TYPE* D,
                  REAL_TYPE* WA,
                  REAL_TYPE* WC,
                  REAL_TYPE* WD,
                  REAL_TYPE* omg,
                  double* res,
                  REAL_TYPE* tmp,
                  double* flop);
  
  
// obsolete.f90
void tdma_0_ (int* nx,
              REAL_TYPE* d,
              REAL_TYPE* a,
              REAL_TYPE* b,
              REAL_TYPE* c,
              REAL_TYPE* w);

void tdma_1_ (int* nx,
              REAL_TYPE* d,
              REAL_TYPE* cf,
              REAL_TYPE* w,
              double* flop);

  void tdma_p_(int* nx,
               REAL_TYPE* d,
               REAL_TYPE* a,
               REAL_TYPE* c,
               REAL_TYPE* w);
  
  void tdma_mp_(int* nx,
                int* mp,
                REAL_TYPE* d,
                REAL_TYPE* a,
                REAL_TYPE* c,
                REAL_TYPE* w);
  
// cz_lsor.f90
  /*
void lsor_pcr_kij_(int* sz,
                   int* idx,
                   int* g,
                   int* pn,
                   REAL_TYPE* x,
                   REAL_TYPE* a,
                   REAL_TYPE* c,
                   REAL_TYPE* d,
                   REAL_TYPE* a1,
                   REAL_TYPE* c1,
                   REAL_TYPE* d1,
                   REAL_TYPE* msk,
                   REAL_TYPE* rhs,
                   REAL_TYPE* omg,
                   double* res,
                   double* flop);
  
void lsor_pcr_kij2_(int* sz,
                    int* idx,
                    int* g,
                    int* pn,
                    REAL_TYPE* x,
                    REAL_TYPE* a,
                    REAL_TYPE* c,
                    REAL_TYPE* d,
                    REAL_TYPE* a1,
                    REAL_TYPE* c1,
                    REAL_TYPE* d1,
                    REAL_TYPE* msk,
                    REAL_TYPE* rhs,
                    REAL_TYPE* omg,
                    double* res,
                    double* flop);
  
void lsor_pcr_kij3_(int* sz,
                    int* idx,
                    int* g,
                    int* pn,
                    REAL_TYPE* x,
                    REAL_TYPE* a,
                    REAL_TYPE* c,
                    REAL_TYPE* d,
                    REAL_TYPE* a1,
                    REAL_TYPE* c1,
                    REAL_TYPE* d1,
                    REAL_TYPE* msk,
                    REAL_TYPE* rhs,
                    REAL_TYPE* omg,
                    double* res,
                    double* flop);
  
void lsor_pcr_kij4_(int* sz,
                    int* idx,
                    int* g,
                    int* pn,
                    REAL_TYPE* x,
                    REAL_TYPE* a,
                    REAL_TYPE* c,
                    REAL_TYPE* d,
                    REAL_TYPE* a1,
                    REAL_TYPE* c1,
                    REAL_TYPE* d1,
                    REAL_TYPE* msk,
                    REAL_TYPE* rhs,
                    REAL_TYPE* omg,
                    double* res,
                    double* flop);
  
void lsor_pcr_kij5_(int* sz,
                    int* idx,
                    int* g,
                    int* pn,
                    REAL_TYPE* x,
                    REAL_TYPE* msk,
                    REAL_TYPE* rhs,
                    REAL_TYPE* omg,
                    double* res,
                    double* flop);
  
void lsor_pcr_kij6_(int* sz,
                    int* idx,
                    int* g,
                    int* pn,
                    REAL_TYPE* x,
                    REAL_TYPE* msk,
                    REAL_TYPE* rhs,
                    REAL_TYPE* omg,
                    double* res,
                    double* flop);
  */


// cz_blas.f90

void imask_k_       (REAL_TYPE* x,
                     int* sz,
                     int* idx,
                     int* g);

void blas_clear_    (REAL_TYPE* x,
                     int* sz,
                     int* g);

void blas_copy_     (REAL_TYPE* dst,
                     REAL_TYPE* src,
                     int* sz,
                     int* g);
  
  
void blas_copy_in_   (REAL_TYPE* dst,
                      REAL_TYPE* src,
                      int* sz,
                      int* g);
  
void blas_triad_    (REAL_TYPE* z,
                     REAL_TYPE* x,
                     REAL_TYPE* y,
                     REAL_TYPE* a,
                     int* sz,
                     int* idx,
                     int* g,
                     double* flop);

void blas_dot1_     (REAL_TYPE* r,
                     REAL_TYPE* p,
                     int* sz,
                     int* idx,
                     int* g,
                     double* flop);

void blas_dot2_     (REAL_TYPE* r,
                     REAL_TYPE* p,
                     REAL_TYPE* q,
                     int* sz,
                     int* idx,
                     int* g,
                     double* flop);

void blas_bicg_1_ (REAL_TYPE* p,
                   REAL_TYPE* r,
                   REAL_TYPE* q,
                   REAL_TYPE* beta,
                   REAL_TYPE* omg,
                   int* sz,
                   int* idx,
                   int* g,
                   double* flop);

void blas_bicg_2_   (REAL_TYPE* z,
                     REAL_TYPE* x,
                     REAL_TYPE* y,
                     REAL_TYPE* a,
                     REAL_TYPE* b,
                     int* sz,
                     int* idx,
                     int* g,
                     double* flop);

void blas_calc_ax_  (REAL_TYPE* ap,
                     REAL_TYPE* p,
                     int* sz,
                     int* idx,
                     int* g,
                     REAL_TYPE* cf,
                     double* flop);

void blas_calc_rk_  (REAL_TYPE* r,
                     REAL_TYPE* p,
                     REAL_TYPE* b,
                     int* sz,
                     int* idx,
                     int* g,
                     REAL_TYPE* cf,
                     double* flop);

void calc_rk_maf_  (REAL_TYPE* r,
                    REAL_TYPE* p,
                    REAL_TYPE* b,
                    int* sz,
                    int* idx,
                    int* g,
                    REAL_TYPE* xc,
                    REAL_TYPE* yc,
                    REAL_TYPE* zc,
                    REAL_TYPE* pvt,
                    double* flop);
  
void calc_ax_maf_  (REAL_TYPE* ap,
                    REAL_TYPE* p,
                    int* sz,
                    int* idx,
                    int* g,
                    REAL_TYPE* xc,
                    REAL_TYPE* yc,
                    REAL_TYPE* zc,
                    REAL_TYPE* pvt,
                    double* flop);

void search_pivot_ (REAL_TYPE* pvt,
                    int* sz,
                    int* idx,
                    int* g,
                    REAL_TYPE* X,
                    REAL_TYPE* Y,
                    REAL_TYPE* Z);
  
  
// utility.f90

void fileout_t_ (int* sz,
               int* g,
               REAL_TYPE* s,
               REAL_TYPE* dh,
               REAL_TYPE* org,
               char* fname);

void exact_t_ (int* sz,
               int* g,
               REAL_TYPE* e,
               REAL_TYPE* dh,
               REAL_TYPE* org);

void err_t_   (int* sz,
               int* idx,
               int* g,
               double* d,
               REAL_TYPE* p,
               REAL_TYPE* e,
               int* loc);
}



#endif // _CZ_F_FUNC_H_
