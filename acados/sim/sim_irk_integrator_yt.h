#ifndef ACADOS_SIM_SIM_IRK_INTEGRATOR_YT_H_
#define ACADOS_SIM_SIM_IRK_INTEGRATOR_YT_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "acados/sim/sim_rk_common_yt.h"
#include "acados/sim/sim_common_yt.h"
#include "acados/utils/types.h"

#include "external/blasfeo/include/blasfeo_target.h"
#include "external/blasfeo/include/blasfeo_common.h"

typedef struct {
    struct d_strmat *JG; // jacobian of G over K (nx*ns, nx*ns)
    struct d_strmat *rG; // residuals of G (nx*ns, 1)
    struct d_strmat *K; // internal variables (nx*ns, 1)
    struct d_strmat *JGf; // jacobian of G over x and u (nx*ns, nx+nu);
    struct d_strmat *JKf; // jacobian of K over x and u (nx*ns, nx+nu);
    double *xt0; // temporary states (nx)
    double *kt; // temporary internal variables (nx, 1)
    double *rGt; // temporary residuals of G (nx, 1)
    double *jac_out; // temporary Jacobian of ode (nx, 2*nx+nu)
    double *Jt; // temporary Jacobian of ode (nx, nx)
    double *ode_args; // pointer to ode args
    double *JFK; // jacobian of F over K
    int *ipiv; // index of pivot vector
}sim_irk_memory;

int irk_calculate_memory_size(sim_RK_opts *opts, sim_in *in);

char *assign_irk_memory(sim_RK_opts *opts, sim_in *in, sim_irk_memory **memory, void *raw_memory);

sim_irk_memory *sim_irk_create_memory(sim_RK_opts *opts, sim_in *in);

int sim_irk_yt(const sim_in *in, sim_out *out, void *opts_, void *mem_);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // ACADOS_SIM_SIM_IRK_INTEGRATOR_H_