#ifndef EXAMPLES_CASADI_CRANE_CRANE_MODEL_H_
#define EXAMPLES_CASADI_CRANE_CRANE_MODEL_H_

#include "acados/utils/types.h"

#ifdef __cplusplus
extern "C" {
#endif

int vdeFun(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int adjFun(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int hessFun(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int impl_odeFun(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int impl_jacFun_x(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int impl_jacFun_xdot(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);
int impl_jacFun_u(const real_t** arg, real_t** res, int* iw, real_t* w, int mem);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif