#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "acados/sim/sim_common.h"
#include "acados/sim/sim_rk_common.h"
#include "acados/sim/sim_erk_integrator.h"
#include "acados/sim/sim_casadi_wrapper.h"

// #include "include/sim/sim_lifted_irk_integrator.h"
#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"

#include "examples/c/chain_model/chain_model.h"

int main() {
    int_t nil;
    int_t NMF=8;

    int_t NX = 6 * NMF;
    int_t NU = 3;

    double *xref;
    xref = (double*)calloc(NX, sizeof(double));
    FILE *refStates;
    refStates = fopen(XN_NM9_FILE, "r");
    for (int_t i = 0; i < NX; i++) {
        nil = fscanf(refStates, "%lf", &xref[i]);
    }
    fclose(refStates);
    if (!nil)
        printf("xref read successfully");

    real_t T = 0.2;
    sim_in sim_in;
    sim_out sim_out;
    sim_info info;

    sim_RK_opts rk_opts;
    void *sim_work = NULL;

    sim_in.num_steps = 4;
    sim_in.step = T / sim_in.num_steps;
    sim_in.nx = NX;
    sim_in.nu = NU;

    sim_in.sens_forw = true;
    sim_in.sens_adj = false;
    sim_in.sens_hess = false;
    sim_in.num_forw_sens = NX + NU;

    sim_in.vde = &vde_chain_nm9;
    sim_in.VDE_forw = &vde_fun;

    sim_in.x = malloc(sizeof(*sim_in.x) * (NX));
    sim_in.u = malloc(sizeof(*sim_in.u) * (NU));
    sim_in.S_forw = malloc(sizeof(*sim_in.S_forw) * (NX * (NX + NU)));

    for (int_t i = 0; i < NX * (NX + NU); i++)
        sim_in.S_forw[i] = 0.0;
    for (int_t i = 0; i < NX; i++)
        sim_in.S_forw[i * (NX + 1)] = 1.0;

    sim_out.xn = malloc(sizeof(*sim_out.xn) * (NX));
    sim_out.S_forw = malloc(sizeof(*sim_out.S_forw) * (NX * (NX + NU)));
    sim_out.info = &info;

    int_t workspace_size;

    sim_erk_create_arguments(&rk_opts, 4);
    workspace_size =
        sim_erk_calculate_workspace_size(&sim_in, &rk_opts);

    sim_work = (void *)malloc(workspace_size);

    for (int_t j = 0; j < NX; j++)
        sim_in.x[j] = xref[j];
    for (int_t j = 0; j < NU; j++)
        sim_in.u[j] = 0.0;

    sim_erk(&sim_in, &sim_out, &rk_opts, 0, sim_work);

    double *xn = sim_out.xn;
    double *S_forw_out = sim_out.S_forw;

    printf("xn: \n");
    for (int ii=0;ii<NX;ii++)
        printf("%4.2f ",xn[ii]);

    printf("\n");
    printf("S_forw_out: \n");
    for (int ii=0;ii<NX;ii++){
        for (int jj=0;jj<sim_in.num_forw_sens;jj++)
            printf("%4.2f ",S_forw_out[jj*NX+ii]);
        printf("\n");
    }

    printf("\n");
    printf("cpt: %8.4f [ms]\n", sim_out.info->CPUtime*1000);
    printf("AD cpt: %8.4f [ms]\n", sim_out.info->ADtime*1000);

    return 0;

}
