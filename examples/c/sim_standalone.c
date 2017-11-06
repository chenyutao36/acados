#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "acados/sim/sim_common_yt.h"
#include "acados/sim/sim_rk_common_yt.h"
#include "acados/sim/sim_erk_integrator_yt.h"
#include "acados/sim/sim_casadi_wrapper.h"

// #include "include/sim/sim_lifted_irk_integrator.h"
#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"

// #include "examples/c/chain_model/chain_model.h"
#include "examples/c/yutao_model/yutao_model.h"

// blasfeo
#include "external/blasfeo/include/blasfeo_target.h"
#include "external/blasfeo/include/blasfeo_common.h"
#include "external/blasfeo/include/blasfeo_d_aux.h"
#include "external/blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "external/blasfeo/include/blasfeo_v_aux_ext_dep.h"
#include "external/blasfeo/include/blasfeo_d_blas.h"

int main() {

    int ii,jj;
    // int nil;

    // int NMF = 8;
    // int nx = 6 * NMF;
    int nx = 4;
    int nu = 1;
    int NF = nx + nu; // columns of forward seed
    int NA = nx; // number of elements of adjoint seed
    int nX = nx *(1+NF);

    double T = 0.05;
    int num_stages = 4;
    double *xref;
    xref = (double*)calloc(nx, sizeof(double));
    xref[1] = 3.14;

    double A_rk[] = {0.0, 0.5, 0.0, 0.0,
        0.0, 0.0, 0.5, 0.0,
        0.0, 0.0, 0.0, 1.0,
        0.0, 0.0, 0.0, 0.0};
    double B_rk[] = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};
    double C_rk[] = {0.0, 0.5, 0.5, 0.0};

    sim_RK_opts *rk_opts = create_sim_RK_opts(num_stages);
    rk_opts->num_stages = num_stages;

    for (ii=0;ii<num_stages*num_stages;ii++){
        rk_opts->A_mat[ii] = A_rk[ii];
    }
    for (ii=0;ii<num_stages;ii++){
        rk_opts->b_vec[ii] = B_rk[ii];
        rk_opts->c_vec[ii] = C_rk[ii];
    }

    sim_in *in = create_sim_in(nx, nu ,NF, NA);

    in->num_steps = 4;
    in->step = T / in->num_steps;
    in->sens_forw = true;
    in->sens_adj = true;
    in->sens_hess = false;

    in->vde = & vdeFun;
    in->adj = & adjFun;
    // in->vde = &vde_chain_nm9;
    in->VDE_forw = &vde_fun;
    in->VDE_adj = &adj_fun;

    printf("nx=%d nu=%d NF=%d NA=%d nX=%d \n",in->nx, in->nu, in->NF, in->NA, nX);

    // FILE *refStates;
    // refStates = fopen(XN_NM9_FILE, "r");
    for (ii = 0; ii < nx; ii++) {
        // nil = fscanf(refStates, "%lf", &xref[ii]);
        in->x[ii] = xref[ii];
    }
    // if (!nil)
    //     printf("xref read successfully");
    for (ii = 0;ii < nu; ii++){
        in->u[ii] = 0.0;
    }


    printf("x0:\n");
    for (ii = 0; ii < nx; ii++)
        printf("%8.5f ",in->x[ii]);
    printf("\n");

    for (ii = 0; ii < nx * NF; ii++)
        in->S_forw[ii] = 0.0;
    for (ii = 0; ii < nx; ii++)
        in->S_forw[ii * (nx + 1)] = 1.0;

    for (ii = 0; ii < nx; ii++)
        in->S_adj[ii] = 1.0;

    sim_erk_memory *erk_mem = sim_erk_create_memory(rk_opts, in);

    sim_out *out = create_sim_out(nx, nu, NF, NA);
    
    int flag = sim_erk_yt(in, out, rk_opts, erk_mem);

    double *xn = out->xn;

    printf("xn: \n");
    for (ii=0;ii<nx;ii++)
        printf("%8.5f ",xn[ii]);

    printf("\n");
    double *S_forw_out = out->S_forw;
    if (in->sens_forw){      
        printf("S_forw_out: \n");
        for (ii=0;ii<nx;ii++){
            for (jj=0;jj<NF;jj++)
                printf("%8.5f ",S_forw_out[jj*nx+ii]);
            printf("\n");
        }
    }

    printf("\n");
    double *S_adj_out = out->S_adj;
    if (in->sens_adj){
        printf("S_adj_out: \n");
        for (ii=0;ii<nx+nu;ii++){
                printf("%8.5f ",S_adj_out[ii]);
        }
    }

    printf("\n");
    printf("cpt: %8.4f [ms]\n", out->info->CPUtime);
    printf("AD cpt: %8.4f [ms]\n", out->info->ADtime);

    struct d_strmat sA;
    void *mA; 
    v_zeros_align(&mA, d_size_strmat(nx, nx+nu));
    d_create_strmat(nx, nx+nu, &sA, mA);
    d_cvt_mat2strmat(nx, nx+nu, S_forw_out, nx, &sA, 0, 0);

    // d_print_strmat(nx,nx+nu,&sA,0,0);

    struct d_strvec sx;
    void *mx; 
    v_zeros_align(&mx, d_size_strvec(nx));
    d_create_strvec(nx, &sx, mx);
    d_cvt_vec2strvec(nx, in->S_adj, &sx, 0);

    struct d_strvec sz;
    void *mz; 
    v_zeros_align(&mz, d_size_strvec(nx+nu));
    d_create_strvec(nx+nu, &sz, mz);
    dgemv_t_libstr(nx, nx+nu, 1.0, &sA, 0, 0, &sx, 0, 0.0, &sz, 0, &sz, 0);
    
    printf("Jac times lambdaX:\n");
    d_print_tran_strvec(nx+nu, &sz, 0);

    free(xref);
    free(rk_opts);
    free(in);
    free(erk_mem);
    free(out);
    v_free_align(mA);
    v_free_align(mx);
    v_free_align(mz);
    
    return flag;
}