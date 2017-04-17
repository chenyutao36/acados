/* This function was automatically generated by CasADi */
#ifdef __cplusplus
extern "C" {
#endif

#ifdef CODEGEN_PREFIX
  #define NAMESPACE_CONCAT(NS, ID) _NAMESPACE_CONCAT(NS, ID)
  #define _NAMESPACE_CONCAT(NS, ID) NS ## ID
  #define CASADI_PREFIX(ID) NAMESPACE_CONCAT(CODEGEN_PREFIX, ID)
#else /* CODEGEN_PREFIX */
  #define CASADI_PREFIX(ID) vde_chain_nm2_ ## ID
#endif /* CODEGEN_PREFIX */

#include <math.h>

#ifndef real_t
#define real_t double
#define to_double(x) (double) x
#define to_int(x) (int) x
#endif /* real_t */

/* Pre-c99 compatibility */
#if __STDC_VERSION__ < 199901L
real_t CASADI_PREFIX(fmin)(real_t x, real_t y) { return x<y ? x : y;}
#define fmin(x,y) CASADI_PREFIX(fmin)(x,y)
real_t CASADI_PREFIX(fmax)(real_t x, real_t y) { return x>y ? x : y;}
#define fmax(x,y) CASADI_PREFIX(fmax)(x,y)
#endif

#define PRINTF printf
real_t CASADI_PREFIX(sq)(real_t x) { return x*x;}
#define sq(x) CASADI_PREFIX(sq)(x)

real_t CASADI_PREFIX(sign)(real_t x) { return x<0 ? -1 : x>0 ? 1 : x;}
#define sign(x) CASADI_PREFIX(sign)(x)

static const int CASADI_PREFIX(s0)[] = {6, 1, 0, 6, 0, 1, 2, 3, 4, 5};
#define s0 CASADI_PREFIX(s0)
static const int CASADI_PREFIX(s1)[] = {6, 6, 0, 6, 12, 18, 24, 30, 36, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5};
#define s1 CASADI_PREFIX(s1)
static const int CASADI_PREFIX(s2)[] = {6, 3, 0, 6, 12, 18, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5};
#define s2 CASADI_PREFIX(s2)
static const int CASADI_PREFIX(s3)[] = {3, 1, 0, 3, 0, 1, 2};
#define s3 CASADI_PREFIX(s3)
/* vdeFun */
int vde_chain_nm2(void* mem, const real_t** arg, real_t** res, int* iw, real_t* w) {
    mem = 0; mem += 0; w = 0; w += 0; iw = 0; iw += 0;
  real_t a0=arg[0] ? arg[0][3] : 0;
  if (res[0]!=0) res[0][0]=a0;
  a0=arg[0] ? arg[0][4] : 0;
  if (res[0]!=0) res[0][1]=a0;
  a0=arg[0] ? arg[0][5] : 0;
  if (res[0]!=0) res[0][2]=a0;
  a0=arg[3] ? arg[3][0] : 0;
  if (res[0]!=0) res[0][3]=a0;
  a0=arg[3] ? arg[3][1] : 0;
  if (res[0]!=0) res[0][4]=a0;
  a0=arg[3] ? arg[3][2] : 0;
  if (res[0]!=0) res[0][5]=a0;
  a0=arg[1] ? arg[1][3] : 0;
  if (res[1]!=0) res[1][0]=a0;
  a0=arg[1] ? arg[1][4] : 0;
  if (res[1]!=0) res[1][1]=a0;
  a0=arg[1] ? arg[1][5] : 0;
  if (res[1]!=0) res[1][2]=a0;
  a0=0.;
  if (res[1]!=0) res[1][3]=a0;
  if (res[1]!=0) res[1][4]=a0;
  if (res[1]!=0) res[1][5]=a0;
  real_t a1=arg[1] ? arg[1][9] : 0;
  if (res[1]!=0) res[1][6]=a1;
  a1=arg[1] ? arg[1][10] : 0;
  if (res[1]!=0) res[1][7]=a1;
  a1=arg[1] ? arg[1][11] : 0;
  if (res[1]!=0) res[1][8]=a1;
  if (res[1]!=0) res[1][9]=a0;
  if (res[1]!=0) res[1][10]=a0;
  if (res[1]!=0) res[1][11]=a0;
  a1=arg[1] ? arg[1][15] : 0;
  if (res[1]!=0) res[1][12]=a1;
  a1=arg[1] ? arg[1][16] : 0;
  if (res[1]!=0) res[1][13]=a1;
  a1=arg[1] ? arg[1][17] : 0;
  if (res[1]!=0) res[1][14]=a1;
  if (res[1]!=0) res[1][15]=a0;
  if (res[1]!=0) res[1][16]=a0;
  if (res[1]!=0) res[1][17]=a0;
  a1=arg[1] ? arg[1][21] : 0;
  if (res[1]!=0) res[1][18]=a1;
  a1=arg[1] ? arg[1][22] : 0;
  if (res[1]!=0) res[1][19]=a1;
  a1=arg[1] ? arg[1][23] : 0;
  if (res[1]!=0) res[1][20]=a1;
  if (res[1]!=0) res[1][21]=a0;
  if (res[1]!=0) res[1][22]=a0;
  if (res[1]!=0) res[1][23]=a0;
  a1=arg[1] ? arg[1][27] : 0;
  if (res[1]!=0) res[1][24]=a1;
  a1=arg[1] ? arg[1][28] : 0;
  if (res[1]!=0) res[1][25]=a1;
  a1=arg[1] ? arg[1][29] : 0;
  if (res[1]!=0) res[1][26]=a1;
  if (res[1]!=0) res[1][27]=a0;
  if (res[1]!=0) res[1][28]=a0;
  if (res[1]!=0) res[1][29]=a0;
  a1=arg[1] ? arg[1][33] : 0;
  if (res[1]!=0) res[1][30]=a1;
  a1=arg[1] ? arg[1][34] : 0;
  if (res[1]!=0) res[1][31]=a1;
  a1=arg[1] ? arg[1][35] : 0;
  if (res[1]!=0) res[1][32]=a1;
  if (res[1]!=0) res[1][33]=a0;
  if (res[1]!=0) res[1][34]=a0;
  if (res[1]!=0) res[1][35]=a0;
  a1=arg[2] ? arg[2][3] : 0;
  if (res[2]!=0) res[2][0]=a1;
  a1=arg[2] ? arg[2][4] : 0;
  if (res[2]!=0) res[2][1]=a1;
  a1=arg[2] ? arg[2][5] : 0;
  if (res[2]!=0) res[2][2]=a1;
  a1=1.;
  if (res[2]!=0) res[2][3]=a1;
  if (res[2]!=0) res[2][4]=a0;
  if (res[2]!=0) res[2][5]=a0;
  real_t a2=arg[2] ? arg[2][9] : 0;
  if (res[2]!=0) res[2][6]=a2;
  a2=arg[2] ? arg[2][10] : 0;
  if (res[2]!=0) res[2][7]=a2;
  a2=arg[2] ? arg[2][11] : 0;
  if (res[2]!=0) res[2][8]=a2;
  if (res[2]!=0) res[2][9]=a0;
  if (res[2]!=0) res[2][10]=a1;
  if (res[2]!=0) res[2][11]=a0;
  a2=arg[2] ? arg[2][15] : 0;
  if (res[2]!=0) res[2][12]=a2;
  a2=arg[2] ? arg[2][16] : 0;
  if (res[2]!=0) res[2][13]=a2;
  a2=arg[2] ? arg[2][17] : 0;
  if (res[2]!=0) res[2][14]=a2;
  if (res[2]!=0) res[2][15]=a0;
  if (res[2]!=0) res[2][16]=a0;
  if (res[2]!=0) res[2][17]=a1;
  return 0;
}

int vde_chain_nm2_init(int* n_in, int* n_out, int* n_int, int* n_real) {
  if (n_in) *n_in = 4;
  if (n_out) *n_out = 3;
  if (n_int) *n_int = 0;
  if (n_real) *n_real = 0;
  return 0;
}

int vde_chain_nm2_alloc(void** mem, const int* idata, const double* rdata) {
  if (mem) *mem = 0;
  (void)idata;
  (void)rdata;
  return 0;
}

int vde_chain_nm2_free(void* mem) {
  (void)mem;
  return 0;
}

int vde_chain_nm2_sparsity(int i, int *nrow, int *ncol, const int **colind, const int **row) {
  const int* s;
  switch (i) {
    case 0:
    case 4:
      s = s0; break;
    case 1:
    case 5:
      s = s1; break;
    case 2:
    case 6:
      s = s2; break;
    case 3:
      s = s3; break;
    default:
      return 1;
  }

  if (nrow) *nrow = s[0];
  if (ncol) *ncol = s[1];
  if (colind) *colind = s + 2;
  if (row) *row = s + 3 + s[1];
  return 0;
}

int vde_chain_nm2_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w) {
  if (sz_arg) *sz_arg = 4;
  if (sz_res) *sz_res = 3;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 3;
  return 0;
}


#ifdef __cplusplus
} /* extern "C" */
#endif
