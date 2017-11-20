/* This function was automatically generated by CasADi */
#ifdef __cplusplus
extern "C" {
#endif

#ifdef CODEGEN_PREFIX
  #define NAMESPACE_CONCAT(NS, ID) _NAMESPACE_CONCAT(NS, ID)
  #define _NAMESPACE_CONCAT(NS, ID) NS ## ID
  #define CASADI_PREFIX(ID) NAMESPACE_CONCAT(CODEGEN_PREFIX, ID)
#else /* CODEGEN_PREFIX */
  #define CASADI_PREFIX(ID) impl_jac_u_ ## ID
#endif /* CODEGEN_PREFIX */

#include <math.h>

#ifndef real_t
#define real_t double
#endif /* real_t */

#define to_double(x) (double) x
#define to_int(x) (int) x
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

static const int CASADI_PREFIX(s0)[12] = {8, 1, 0, 8, 0, 1, 2, 3, 4, 5, 6, 7};
#define s0 CASADI_PREFIX(s0)
static const int CASADI_PREFIX(s1)[6] = {2, 1, 0, 2, 0, 1};
#define s1 CASADI_PREFIX(s1)
static const int CASADI_PREFIX(s2)[21] = {8, 2, 0, 8, 16, 0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7};
#define s2 CASADI_PREFIX(s2)
/* impl_jacFun_u */
int impl_jacFun_u(const real_t** arg, real_t** res, int* iw, real_t* w, int mem) {
  real_t a0=0.;
  if (res[0]!=0) res[0][0]=a0;
  if (res[0]!=0) res[0][1]=a0;
  if (res[0]!=0) res[0][2]=a0;
  if (res[0]!=0) res[0][3]=a0;
  real_t a1=-1.;
  if (res[0]!=0) res[0][4]=a1;
  if (res[0]!=0) res[0][5]=a0;
  if (res[0]!=0) res[0][6]=a0;
  real_t a2=arg[0] ? arg[0][6] : 0;
  a2=cos(a2);
  real_t a3=4.7418203070092001e-02;
  a3=(a3*a2);
  a2=arg[0] ? arg[0][2] : 0;
  a3=(a3/a2);
  if (res[0]!=0) res[0][7]=a3;
  if (res[0]!=0) res[0][8]=a0;
  if (res[0]!=0) res[0][9]=a0;
  if (res[0]!=0) res[0][10]=a0;
  if (res[0]!=0) res[0][11]=a0;
  if (res[0]!=0) res[0][12]=a0;
  if (res[0]!=0) res[0][13]=a1;
  if (res[0]!=0) res[0][14]=a0;
  if (res[0]!=0) res[0][15]=a0;
  return 0;
}

void impl_jacFun_u_incref(void) {
}

void impl_jacFun_u_decref(void) {
}

int impl_jacFun_u_n_in(void) { return 3;}

int impl_jacFun_u_n_out(void) { return 1;}

const char* impl_jacFun_u_name_in(int i){
  switch (i) {
  case 0: return "i0";
  case 1: return "i1";
  case 2: return "i2";
  default: return 0;
  }
}

const char* impl_jacFun_u_name_out(int i){
  switch (i) {
  case 0: return "o0";
  default: return 0;
  }
}

const int* impl_jacFun_u_sparsity_in(int i) {
  switch (i) {
  case 0: return s0;
  case 1: return s0;
  case 2: return s1;
  default: return 0;
  }
}

const int* impl_jacFun_u_sparsity_out(int i) {
  switch (i) {
  case 0: return s2;
  default: return 0;
  }
}

int impl_jacFun_u_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w) {
  if (sz_arg) *sz_arg = 3;
  if (sz_res) *sz_res = 1;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 4;
  return 0;
}


#ifdef __cplusplus
} /* extern "C" */
#endif