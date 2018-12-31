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

#ifndef _CZ_F_FUNC_H_
#define _CZ_F_FUNC_H_


extern "C" {

// cz_lsolver.f90

void bc_      (int* sz,
               int* g,
               REAL_TYPE* e,
               REAL_TYPE* dh,
               REAL_TYPE* org,
               int* nID);

void bc_k_    (int* sz,
               int* g,
               REAL_TYPE* e,
               REAL_TYPE* dh,
               REAL_TYPE* org,
               int* nID);

void bc_ikj_  (int* sz,
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


// cz_losr.f90
void lsor_lu_rhs_jd_(REAL_TYPE* d,
                   int* sz,
                   int* idx,
                   int* g,
                   int* j,
                   REAL_TYPE* x,
                   REAL_TYPE* rhs,
                   REAL_TYPE* msk,
                   double* flop);

void lsor_relax_d_(REAL_TYPE* d,
                 int* sz,
                 int* idx,
                 int* g,
                 int* j,
                 REAL_TYPE* x,
                 REAL_TYPE* msk,
                 REAL_TYPE* omg,
                 double* res,
                 double* flop);

void lsor_lu_rhs_j_(REAL_TYPE* d,
                   int* sz,
                   int* idx,
                   int* g,
                   int* j,
                   REAL_TYPE* x,
                   REAL_TYPE* rhs,
                   double* flop);

void lsor_lu_rhs_k_(REAL_TYPE* d,
                   int* sz,
                   int* idx,
                   int* g,
                   int* j,
                   int* k,
                   REAL_TYPE* x,
                   REAL_TYPE* d2,
                   REAL_TYPE* msk,
                   double* flop);

void lsor_lu_bc_kst_(REAL_TYPE* d,
                     int* sz,
                     int* idx,
                     int* g,
                     int* j,
                     REAL_TYPE* rhs,
                     REAL_TYPE* msk,
                     double* flop);

void lsor_lu_bc_ked_(REAL_TYPE* d,
                     int* sz,
                     int* idx,
                     int* g,
                     int* j,
                     REAL_TYPE* rhs,
                     REAL_TYPE* msk,
                     double* flop);

void lsor_tdma_f_(REAL_TYPE* d,
                  int* sz,
                  int* idx,
                  int*g,
                  int* j,
                  int* k,
                  REAL_TYPE* a,
                  REAL_TYPE* e,
                  double* flop);

void lsor_tdma_b_(REAL_TYPE* d,
                  int* sz,
                  int* idx,
                  int*g,
                  int* j,
                  REAL_TYPE* w,
                  double* flop);

void lsor_relax_(REAL_TYPE* d,
                 int* sz,
                 int* idx,
                 int* g,
                 int* j,
                 REAL_TYPE* x,
                 REAL_TYPE* msk,
                 REAL_TYPE* omg,
                 double* res,
                 double* flop);

void lsor_lu_fwd_(REAL_TYPE* d,
                   int* sz,
                   int* idx,
                   int* g,
                   int* j,
                   REAL_TYPE* x,
                   REAL_TYPE* a,
                   REAL_TYPE* e,
                   REAL_TYPE* d2,
                   REAL_TYPE* msk,
                   REAL_TYPE* rhs,
                   double* flop);

void lsor_lu_rhs_k4_(REAL_TYPE* d,
                   int* sz,
                   int* idx,
                   int* g,
                   int* j,
                   int* k,
                   REAL_TYPE* x,
                   REAL_TYPE* d2,
                   REAL_TYPE* msk,
                   REAL_TYPE* rhs,
                   double* flop);

void lsor_lu_bc_kst4_(REAL_TYPE* d,
                     int* sz,
                     int* idx,
                     int* g,
                     int* j,
                     REAL_TYPE* rhs,
                     REAL_TYPE* msk,
                     double* flop);

void lsor_lu_bc_ked4_(REAL_TYPE* d,
                     int* sz,
                     int* idx,
                     int* g,
                     int* j,
                     REAL_TYPE* rhs,
                     REAL_TYPE* msk,
                     double* flop);

void lsor_tdma_f4_(REAL_TYPE* d,
                  int* sz,
                  int* idx,
                  int*g,
                  int* j,
                  REAL_TYPE* d2,
                  REAL_TYPE* x,
                  REAL_TYPE* a,
                  REAL_TYPE* e,
                  REAL_TYPE* msk,
                  REAL_TYPE* rhs,
                  double* flop);

void lsor_tdma_fwd_(REAL_TYPE* d,
                  int* sz,
                  int* idx,
                  int*g,
                  int* j,
                  REAL_TYPE* d2,
                  REAL_TYPE* x,
                  REAL_TYPE* a,
                  REAL_TYPE* e,
                  REAL_TYPE* msk,
                  REAL_TYPE* rhs,
                  double* flop);

void lsor_inner_(REAL_TYPE* d,
                  int* sz,
                  int* idx,
                  int* g,
                  int* j,
                  REAL_TYPE* d2,
                  REAL_TYPE* x,
                  REAL_TYPE* a,
                  REAL_TYPE* e,
                  REAL_TYPE* w,
                  REAL_TYPE* msk,
                  REAL_TYPE* rhs,
                  REAL_TYPE* omg,
                  double* res,
                  double* flop);

void lsor_inner_b_(REAL_TYPE* d,
                  int* sz,
                  int* idx,
                  int* g,
                  REAL_TYPE* d2,
                  REAL_TYPE* x,
                  REAL_TYPE* a,
                  REAL_TYPE* e,
                  REAL_TYPE* w,
                  REAL_TYPE* msk,
                  REAL_TYPE* rhs,
                  REAL_TYPE* omg,
                  int* ItrInner,
                  double* rd,
                  double* flop);

void lsor_inner_b4_(REAL_TYPE* d,
                  int* sz,
                  int* idx,
                  int* g,
                  REAL_TYPE* d2,
                  REAL_TYPE* x,
                  REAL_TYPE* a,
                  REAL_TYPE* e,
                  REAL_TYPE* w,
                  REAL_TYPE* msk,
                  REAL_TYPE* rhs,
                  REAL_TYPE* omg,
                  int* ItrInner,
                  double* rd,
                  double* flop);

void lsor_inner_c_(REAL_TYPE* d,
                  int* sz,
                  int* idx,
                  int* g,
                  REAL_TYPE* d2,
                  REAL_TYPE* x,
                  REAL_TYPE* a,
                  REAL_TYPE* e,
                  REAL_TYPE* w,
                  REAL_TYPE* msk,
                  REAL_TYPE* rhs,
                  REAL_TYPE* omg,
                  int* cl_sz,
                  int* nt,
                  int* itrinner,
                  double* rd,
                  double* flop);

void lsor_tdma_b4_(REAL_TYPE* d,
                  int* sz,
                  int* idx,
                  int*g,
                  int* j,
                  REAL_TYPE* w,
                  double* flop);

void lsor_relax4_(REAL_TYPE* d,
                 int* sz,
                 int* idx,
                 int* g,
                 int* j,
                 REAL_TYPE* x,
                 REAL_TYPE* msk,
                 REAL_TYPE* omg,
                 double* res,
                 double* flop);


void tdma_lsor_a_(REAL_TYPE* d,
                int* sz,
                int* idx,
                int* g,
                REAL_TYPE* x,
                REAL_TYPE* w,
                REAL_TYPE* a,
                REAL_TYPE* b,
                REAL_TYPE* c,
                REAL_TYPE* rhs,
                REAL_TYPE* omg,
                double* res,
                double* flop);

void tdma_lsor_b_(REAL_TYPE* d,
                int* sz,
                int* idx,
                int* g,
                REAL_TYPE* x,
                REAL_TYPE* w,
                REAL_TYPE* a,
                REAL_TYPE* c,
                REAL_TYPE* rhs,
                REAL_TYPE* omg,
                double* res,
                double* flop);

void tdma_lsor_d_(REAL_TYPE* d,
                int* sz,
                int* idx,
                int* g,
                REAL_TYPE* x,
                REAL_TYPE* w,
                REAL_TYPE* a,
                REAL_TYPE* b,
                REAL_TYPE* c,
                REAL_TYPE* rhs,
                REAL_TYPE* omg,
                double* res,
                double* flop);

void tdma_lsor_e_(REAL_TYPE* d,
                int* sz,
                int* idx,
                int* g,
                REAL_TYPE* x,
                REAL_TYPE* w,
                REAL_TYPE* a,
                REAL_TYPE* c,
                REAL_TYPE* rhs,
                REAL_TYPE* omg,
                double* res,
                double* flop);

void tdma_lsor_f_(REAL_TYPE* d,
                int* sz,
                int* idx,
                int* g,
                REAL_TYPE* x,
                REAL_TYPE* w,
                REAL_TYPE* e,
                REAL_TYPE* a,
                REAL_TYPE* rhs,
                REAL_TYPE* omg,
                double* res,
                double* flop);

// cz_ljcb.f90

void ms_rhs8_(REAL_TYPE* d,
              int* sz,
              int* idx,
              int* g,
              int* ia,
              int* ja,
              REAL_TYPE* x,
              REAL_TYPE* rhs,
              REAL_TYPE* msk,
              double* flop);

void ms_relax8_(REAL_TYPE* d,
                int* sz,
                int* idx,
                int* g,
                int* ia,
                int* ja,
                REAL_TYPE* x,
                REAL_TYPE* msk,
                REAL_TYPE* omg,
                double* res,
                double* flop);

void ms_bc_(REAL_TYPE* d,
            int* sz,
            int* idx,
            int* g,
            int* ia,
            int* ja,
            REAL_TYPE* rhs,
            REAL_TYPE* msk,
            double* flop);

void ms_tdma_(REAL_TYPE* d,
              int* sz,
              int* idx,
              int*g,
              int* ia,
              int* ja,
              REAL_TYPE* a,
              REAL_TYPE* e,
              REAL_TYPE* w,
              double* flop);

void tdma_ljcb_a_(REAL_TYPE* d,
                int* sz,
                int* idx,
                int* g,
                REAL_TYPE* x,
                REAL_TYPE* w,
                REAL_TYPE* a,
                REAL_TYPE* b,
                REAL_TYPE* c,
                REAL_TYPE* rhs,
                REAL_TYPE* omg,
                double* res,
                double* flop);

void tdma_ljcb_b_(REAL_TYPE* d,
                int* sz,
                int* idx,
                int* g,
                REAL_TYPE* x,
                REAL_TYPE* w,
                REAL_TYPE* a,
                REAL_TYPE* b,
                REAL_TYPE* c,
                REAL_TYPE* rhs,
                REAL_TYPE* m,
                REAL_TYPE* omg,
                double* res,
                double* flop);

void ljcb_f0_(REAL_TYPE* d,
              int* sz,
              int* idx,
              int* g,
              REAL_TYPE* x,
              REAL_TYPE* rhs,
              double* flop);

void ljcb_f04_(REAL_TYPE* d,
              int* sz,
              int* idx,
              int* g,
              REAL_TYPE* x,
              REAL_TYPE* rhs,
              double* flop);

void ljcb_f1_(REAL_TYPE* d,
              int* sz,
              int* idx,
              int* g,
              REAL_TYPE* w,
              REAL_TYPE* b,
              REAL_TYPE* c,
              REAL_TYPE* rhs,
              double* flop);

void ljcb_f2_(REAL_TYPE* d,
              int* sz,
              int* idx,
              int* g,
              REAL_TYPE* w,
              REAL_TYPE* a,
              REAL_TYPE* b,
              REAL_TYPE* c,
              double* flop);

void ljcb_f3_(REAL_TYPE* d,
              int* sz,
              int* idx,
              int* g,
              REAL_TYPE* w,
              double* flop);

void ljcb_f4_(REAL_TYPE* d,
              int* sz,
              int* idx,
              int* g,
              REAL_TYPE* x,
              REAL_TYPE* omg,
              double* res,
              double* flop);

void ljcb_f0t_(REAL_TYPE* d,
              int* sz,
              int* idx,
              int* g,
              REAL_TYPE* x,
              REAL_TYPE* rhs,
              double* flop);

void ljcb_f1t_(REAL_TYPE* d,
              int* sz,
              int* idx,
              int* g,
              REAL_TYPE* w,
              REAL_TYPE* b,
              REAL_TYPE* c,
              REAL_TYPE* rhs,
              double* flop);

void ljcb_f2t_(REAL_TYPE* d,
              int* sz,
              int* idx,
              int* g,
              REAL_TYPE* w,
              REAL_TYPE* a,
              REAL_TYPE* b,
              REAL_TYPE* c,
              double* flop);

void ljcb_f3t_(REAL_TYPE* d,
              int* sz,
              int* idx,
              int* g,
              REAL_TYPE* w,
              double* flop);

void ljcb_f4t_(REAL_TYPE* d,
              int* sz,
              int* idx,
              int* g,
              REAL_TYPE* x,
              REAL_TYPE* omg,
              double* res,
              double* flop);

void ljcb_g1_(REAL_TYPE* d,
              int* sz,
              int* idx,
              int* g,
              REAL_TYPE* rhs,
              double* flop);

void ljcb_g2_(REAL_TYPE* d,
              int* sz,
              int* idx,
              int* g,
              REAL_TYPE* e,
              REAL_TYPE* a,
              double* flop);

void ljcb_g3_(REAL_TYPE* d,
              int* sz,
              int* idx,
              int* g,
              REAL_TYPE* w,
              double* flop);

// cz_blas.f90
void imask_ikj_     (REAL_TYPE* x,
                     int* sz,
                     int* idx,
                     int* g);

void imask_k_       (REAL_TYPE* x,
                     int* sz,
                     int* idx,
                     int* g);

void init_mask_     (REAL_TYPE* x,
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

void blas_triad_    (REAL_TYPE* z,
                     REAL_TYPE* x,
                     REAL_TYPE* y,
                     double* a,
                     int* sz,
                     int* idx,
                     int* g,
                     double* flop);

void blas_dot1_     (double* r,
                     REAL_TYPE* p,
                     int* sz,
                     int* idx,
                     int* g,
                     double* flop);

void blas_dot2_     (double* r,
                     REAL_TYPE* p,
                     REAL_TYPE* q,
                     int* sz,
                     int* idx,
                     int* g,
                     double* flop);

void blas_bicg_1_ (REAL_TYPE* p,
                   REAL_TYPE* r,
                   REAL_TYPE* q,
                   double* beta,
                   double* omg,
                   int* sz,
                   int* idx,
                   int* g,
                   double* flop);

void blas_bicg_2_   (REAL_TYPE* z,
                     REAL_TYPE* x,
                     REAL_TYPE* y,
                     double* a,
                     double* b,
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

void blas_calc_r2_  (double* res,
                     REAL_TYPE* p,
                     REAL_TYPE* b,
                     int* sz,
                     int* idx,
                     int* g,
                     REAL_TYPE* cf,
                     double* flop);

// utility.f90
void fileout_ (int* sz,
               int* g,
               REAL_TYPE* s,
               REAL_TYPE* dh,
               REAL_TYPE* org,
               char* fname);

void exact_   (int* sz,
               int* g,
               REAL_TYPE* e,
               REAL_TYPE* dh,
               REAL_TYPE* org);

void err_     (int* sz,
               int* idx,
               int* g,
               double* d,
               REAL_TYPE* p,
               REAL_TYPE* e,
               int* loc);

void fileout_t_ (int* sz,
               int* g,
               REAL_TYPE* s,
               REAL_TYPE* dh,
               REAL_TYPE* org,
               char* fname);

void fileout_ikj_ (int* sz,
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

void exact_ikj_ (int* sz,
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

void err_ikj_   (int* sz,
               int* idx,
               int* g,
               double* d,
               REAL_TYPE* p,
               REAL_TYPE* e,
               int* loc);
}



#endif // _CZ_F_FUNC_H_
