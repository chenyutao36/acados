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

#include "examples/c/chain_model/chain_model.h"

int main() {

    int ii,jj,nil;

    int NMF = 1;
    int nx = 6 * NMF;
    int nu = 3;
    int NF = nx + nu;
    int NA = 0;
    int nX = nx *(1+NF);

    double T = 0.2;
    int num_stages = 4;
    double *xref;
    xref = (double*)calloc(nx, sizeof(double));

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

    printf("A:\n");
    for (ii=0;ii<num_stages;ii++){
        for (jj=0;jj<num_stages;jj++)
            printf("%4.2f ",rk_opts->A_mat[jj*num_stages+ii]);
        printf("\n");
    }
    printf("num_stages: %d", rk_opts->num_stages);
        
    printf("\n");

    printf("b:");
    for (ii=0;ii<num_stages;ii++)
        printf("%4.2f ",rk_opts->b_vec[ii]);
    printf("\n");

    sim_in *in = create_sim_in(nx, nu ,NF, NA);

    in->num_steps = 4;
    in->step = T / in->num_steps;
    in->sens_forw = true;
    in->sens_adj = false;
    in->sens_hess = false;

    in->vde = &vde_chain_nm2;
    in->VDE_forw = &vde_fun;
    // in->jac = &jac_chain_nm2;
    // in->jac_fun = &jac_fun;

    printf("nx=%d nu=%d NF=%d NA=%d nX=%d \n",in->nx, in->nu, in->NF, in->NA, nX);

    FILE *refStates;
    refStates = fopen(XN_NM2_FILE, "r");
    // refStates = fopen("/home/yutaochen/Documents/MATLAB/Packages/acados/examples/c/chain_model/xN_nm2.txt","r");
    for (ii = 0; ii < nx; ii++) {
        nil = fscanf(refStates, "%lf", &xref[ii]);
        in->x[ii] = xref[ii];
    }
    printf("scanf return: %d\n", nil);
    fclose(refStates);
    for (ii = 0;ii < nu; ii++){
        in->u[ii] = 1.0;
    }

    printf("x:");
    for (ii = 0; ii < nx; ii++)
        printf("%4.2f ",in->x[ii]);
    printf("\n");

    for (ii = 0; ii < nx * NF; ii++)
        in->S_forw[ii] = 0.0;
    for (ii = 0; ii < nx; ii++)
        in->S_forw[ii * (nx + 1)] = 1.0;

    double *S_forw_in = in->S_forw;
    printf("S_forw_in: \n");
    for (ii=0;ii<nx;ii++){
        for (jj=0;jj<NF;jj++)
            printf("%4.2f ",S_forw_in[jj*nx+ii]);
        printf("\n");
    }

    sim_erk_memory *erk_mem = sim_erk_create_memory(rk_opts, in);

    sim_out *out = create_sim_out(nx, nu, NF, NA);
    
    int flag = sim_erk_yt(in, out, rk_opts, erk_mem);

    double *xn = out->xn;

    printf("xn: ");
    for (ii=0;ii<nx;ii++)
        printf("%4.2f ",xn[ii]);

    printf("\n");
    if (in->sens_forw){
        double *S_forw_out = out->S_forw;
        printf("S_forw_out: \n");
        for (ii=0;ii<nx;ii++){
            for (jj=0;jj<NF;jj++)
                printf("%4.2f ",S_forw_out[jj*nx+ii]);
            printf("\n");
        }
    }

    printf("\n");
    double *rhs_forw_in = erk_mem->rhs_forw_in;
    printf("rhs_forw_in:\n");
    for (ii=0;ii<nX+nu;ii++){
        printf("%4.2f ",rhs_forw_in[ii]);
    }

    printf("\n");
    printf("cpt: %4.2f [ms]\n", out->info->CPUtime*1000);
    printf("AD cpt: %4.2f [ms]\n", out->info->ADtime*1000);
     
    free(xref);
    free(rk_opts);
    free(in);
    free(erk_mem);
    free(out);
    

    return flag;
}