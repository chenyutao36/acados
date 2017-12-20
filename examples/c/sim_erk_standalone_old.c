#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "acados/sim/sim_common_new.h"
#include "acados/sim/sim_rk_common_new.h"
#include "acados/sim/sim_erk_integrator_new.h"
#include "acados/sim/sim_casadi_wrapper.h"

#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"

#include "examples/c/crane_model/crane_model.h"


#define M_PI 3.14159265358979323846

int main() {

    int ii, jj;
    
    int nx = 4;
    int nu = 1;
    int NF = nx + nu; // columns of forward seed

    double T = 0.05;
    int num_stages = 4;
    double *xref;
    xref = (double*)calloc(nx, sizeof(double));
    xref[1] = M_PI;

    double A_rk[] = {0.0, 0.5, 0.0, 0.0,
        0.0, 0.0, 0.5, 0.0,
        0.0, 0.0, 0.0, 1.0,
        0.0, 0.0, 0.0, 0.0};
    double B_rk[] = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};


    sim_RK_opts *erk_opts = create_sim_RK_opts(num_stages);

    for (ii=0;ii<num_stages*num_stages;ii++){
        erk_opts->A_mat[ii] = A_rk[ii];
    }
    for (ii=0;ii<num_stages;ii++){
        erk_opts->b_vec[ii] = B_rk[ii];
    }

    sim_in *in = create_sim_in(nx, nu ,NF);

    in->num_steps = 4;
    in->step = T / in->num_steps;
    in->sens_forw = false;
    in->sens_adj = false;
    in->sens_hess = false;

    in->vde = & vdeFun;
    in->VDE_forw = &vde_fun;
    in->adj = & adjFun;
    in->VDE_adj = &vde_adj_fun;
    in->hess = &hessFun;
    in->Hess_fun = &vde_hess_fun;

    for (ii = 0; ii < nx; ii++) {
        in->x[ii] = xref[ii];
    }
    for (ii = 0;ii < nu; ii++){
        in->u[ii] = 1.0;
    }

    for (ii = 0; ii < nx * NF; ii++)
        in->S_forw[ii] = 0.0;
    for (ii = 0; ii < nx; ii++)
        in->S_forw[ii * (nx + 1)] = 1.0;

    for (ii = 0; ii < nx; ii++)
        in->S_adj[ii] = 1.0;

    sim_erk_memory *erk_mem = sim_erk_create_memory(erk_opts, in);

    sim_out *out = create_sim_out(nx, nu, NF);
    
    int flag = sim_erk_new(in, out, erk_opts, erk_mem);

    double *xn = out->xn;

    printf("\nxn: \n");
    for (ii=0;ii<nx;ii++)
        printf("%8.5f ",xn[ii]);
    printf("\n");

    double *S_forw_out = out->S_forw;  
    if (in->sens_forw){
        printf("\nS_forw_out: \n");
        for (ii=0;ii<nx;ii++){
            for (jj=0;jj<NF;jj++)
                printf("%8.5f ",S_forw_out[jj*nx+ii]);
            printf("\n");
        }
    }
     
    if (in->sens_adj){
        double *S_adj_out = out->S_adj;
        printf("\nS_adj_out: \n");
        for (ii=0;ii<nx+nu;ii++){
            printf("%8.5f ",S_adj_out[ii]);
        }
        printf("\n"); 
    }

    double zero = 0.0;
    if(in->sens_hess){ 
        double *S_hess_out = out->S_hess;
        printf("\nS_hess_out: \n");
        for (ii=0;ii<NF;ii++){
            for (jj=0;jj<NF;jj++){
                if (jj>ii){
                    printf("%8.5f ",zero);
                }else{
                    printf("%8.5f ",S_hess_out[jj*NF+ii]);
                }
            }
            printf("\n");
        }
    }

    printf("\n");
    printf("cpt: %8.4f [ms]\n", out->info->CPUtime);
    printf("AD cpt: %8.4f [ms]\n", out->info->ADtime);

    free(xref);
    free(erk_opts);
    free(in);
    free(erk_mem);
    free(out);
    
    return flag;
}