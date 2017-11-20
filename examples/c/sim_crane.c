#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "acados/sim/sim_common_yt.h"
#include "acados/sim/sim_rk_common_yt.h"
#include "acados/sim/sim_erk_integrator_yt.h"
#include "acados/sim/sim_irk_integrator_yt.h"
#include "acados/sim/sim_casadi_wrapper.h"

#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"

#include "examples/c/crane_model/crane_model.h"

// blasfeo
#include "external/blasfeo/include/blasfeo_target.h"
#include "external/blasfeo/include/blasfeo_common.h"
#include "external/blasfeo/include/blasfeo_d_aux.h"
#include "external/blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "external/blasfeo/include/blasfeo_v_aux_ext_dep.h"
#include "external/blasfeo/include/blasfeo_d_blas.h"

int main() {

    int ii;
    int jj;
    
    int nx = 8;
    int nu = 2;
    int NF = nx + nu; // columns of forward seed
    // u0 = 40.108149413030752 -50.446662212534974

    double T = 0.1;
    int num_stages = 4;
    double *xref;
    xref = (double*)calloc(nx, sizeof(double));
    xref[2] = 0.8;

    double A_rk[] = {0.0, 0.5, 0.0, 0.0,
        0.0, 0.0, 0.5, 0.0,
        0.0, 0.0, 0.0, 1.0,
        0.0, 0.0, 0.0, 0.0};
    double B_rk[] = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};
    // double C_rk[] = {0.0, 0.5, 0.5, 0.0};
    double u_0[] = {40.108149413030752, -50.446662212534974};

    sim_RK_opts *erk_opts = create_sim_RK_opts(num_stages);

    for (ii=0;ii<num_stages*num_stages;ii++){
        erk_opts->A_mat[ii] = A_rk[ii];
    }
    for (ii=0;ii<num_stages;ii++){
        erk_opts->b_vec[ii] = B_rk[ii];
        // rk_opts->c_vec[ii] = C_rk[ii];
    }

    sim_in *in = create_sim_in(nx, nu ,NF);

    in->num_steps = 1;
    in->step = T / in->num_steps;
    in->sens_forw = true;
    in->sens_adj = false;
    in->sens_hess = false;

  //  printf("\n in.step = %8.5f \n ",in->step);

    in->vde = & vdeFun;
    in->adj = & adjFun;
    in->VDE_forw = &vde_fun;
    in->VDE_adj = &adj_fun;
    // in->hess = &hessFun;
    // in->Hess_fun = hess_fun;

    for (ii = 0; ii < nx; ii++) {
        in->x[ii] = xref[ii];
    }

    
    // set controls
    for (ii = 0; ii < nu; ii++) {
        in->u[ii] = u_0[ii];
    }


    printf("in u0 = \n");
    for (ii=0;ii<nu;ii++){
        printf("%8.5f ",in->u[ii]);
    }
    printf("\n in x0 = \n");
    for (ii=0;ii<nx;ii++){
        printf("%8.5f ",in->x[ii]);
    }

    for (ii = 0; ii < nx * NF; ii++)
        in->S_forw[ii] = 0.0;
    for (ii = 0; ii < nx; ii++)
        in->S_forw[ii * (nx + 1)] = 1.0;
        
    for (ii = 0; ii < nx; ii++)
        in->S_adj[ii] = 1.0;

    sim_erk_memory *erk_mem = sim_erk_create_memory(erk_opts, in);

    sim_out *out = create_sim_out(nx, nu, NF);
    
    int flag = sim_erk_yt(in, out, erk_opts, erk_mem);

    double *xn = out->xn;

    printf("\nxn: \n");
    for (ii=0;ii<nx;ii++)
        printf("%8.5f ",xn[ii]);
    printf("\n");

    double *S_forw_out = out->S_forw;  
    printf("\nS_forw_out: \n");
    for (ii=0;ii<nx;ii++){
        for (jj=0;jj<NF;jj++)
            printf("%8.5f ",S_forw_out[jj*nx+ii]);
        printf("\n");
    }
     
    double *S_adj_out = out->S_adj;
    printf("\nS_adj_out: \n");
    for (ii=0;ii<nx+nu;ii++){
        printf("%8.5f ",S_adj_out[ii]);
    }
    printf("\n"); 
 

    printf("\n");
    printf("cpt: %8.4f [ms]\n", out->info->CPUtime);
    printf("AD cpt: %8.4f [ms]\n", out->info->ADtime);

    struct d_strmat sA;
    d_create_strmat(nx, nx+nu, &sA, S_forw_out);

    struct d_strvec sx;
    d_create_strvec(nx, &sx, in->S_adj);

/*    struct d_strvec sz;
    void *mz; 
    v_zeros_align(&mz, d_size_strvec(nx+nu));
    d_create_strvec(nx+nu, &sz, mz);
    dgemv_t_libstr(nx, nx+nu, 1.0, &sA, 0, 0, &sx, 0, 0.0, &sz, 0, &sz, 0);
    
    printf("\nJac times lambdaX:\n");
    d_print_tran_strvec(nx+nu, &sz, 0); */

    free(xref);
    free(erk_opts);
    free(in);
    free(erk_mem);
    free(out);
 //   v_free_align(mz);
    
    return flag;
}