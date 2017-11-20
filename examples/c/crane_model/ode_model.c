/* This function was automatically generated by CasADi */
#ifdef __cplusplus
extern "C" {
#endif

#ifdef CODEGEN_PREFIX
  #define NAMESPACE_CONCAT(NS, ID) _NAMESPACE_CONCAT(NS, ID)
  #define _NAMESPACE_CONCAT(NS, ID) NS ## ID
  #define CASADI_PREFIX(ID) NAMESPACE_CONCAT(CODEGEN_PREFIX, ID)
#else /* CODEGEN_PREFIX */
  #define CASADI_PREFIX(ID) ode_model_ ## ID
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
/* odeFun */
int odeFun(const real_t** arg, real_t** res, int* iw, real_t* w, int mem) {
  real_t a0=arg[0] ? arg[0][1] : 0;
  if (res[0]!=0) res[0][0]=a0;
  real_t a1=4.7418203070092001e-02;
  real_t a2=arg[0] ? arg[0][4] : 0;
  a2=(a1*a2);
  a0=(a0-a2);
  a2=-7.8182378879940387e+01;
  a2=(a2*a0);
  if (res[0]!=0) res[0][1]=a2;
  a2=arg[0] ? arg[0][3] : 0;
  if (res[0]!=0) res[0][2]=a2;
  a0=3.4087337273385997e-02;
  real_t a3=arg[0] ? arg[0][5] : 0;
  a0=(a0*a3);
  a0=(a2-a0);
  a3=-4.0493711676434543e+01;
  a3=(a3*a0);
  if (res[0]!=0) res[0][3]=a3;
  a3=arg[1] ? arg[1][0] : 0;
  if (res[0]!=0) res[0][4]=a3;
  a0=arg[1] ? arg[1][1] : 0;
  if (res[0]!=0) res[0][5]=a0;
  a0=arg[0] ? arg[0][7] : 0;
  if (res[0]!=0) res[0][6]=a0;
  a1=(a1*a3);
  a3=arg[0] ? arg[0][6] : 0;
  real_t a4=cos(a3);
  a1=(a1*a4);
  a3=sin(a3);
  a4=9.8100000000000005e+00;
  a4=(a4*a3);
  a1=(a1+a4);
  a4=2.;
  a4=(a4*a2);
  a4=(a4*a0);
  a1=(a1+a4);
  a4=arg[0] ? arg[0][2] : 0;
  a1=(a1/a4);
  a1=(-a1);
  if (res[0]!=0) res[0][7]=a1;
  return 0;
}

void odeFun_incref(void) {
}

void odeFun_decref(void) {
}

int odeFun_n_in(void) { return 2;}

int odeFun_n_out(void) { return 1;}

const char* odeFun_name_in(int i){
  switch (i) {
  case 0: return "i0";
  case 1: return "i1";
  default: return 0;
  }
}

const char* odeFun_name_out(int i){
  switch (i) {
  case 0: return "o0";
  default: return 0;
  }
}

const int* odeFun_sparsity_in(int i) {
  switch (i) {
  case 0: return s0;
  case 1: return s1;
  default: return 0;
  }
}

const int* odeFun_sparsity_out(int i) {
  switch (i) {
  case 0: return s0;
  default: return 0;
  }
}

int odeFun_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w) {
  if (sz_arg) *sz_arg = 2;
  if (sz_res) *sz_res = 1;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 5;
  return 0;
}


#ifdef __cplusplus
} /* extern "C" */
#endif
