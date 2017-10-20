#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_condensing.h"
#include "acados/ocp_qp/ocp_qp_condensing_common.h"

#include <stdlib.h>

#include "acados/utils/math.h"

#include "blasfeo/include/blasfeo_target.h"
#include "blasfeo/include/blasfeo_common.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_i_aux_ext_dep.h"

#include "hpipm/include/hpipm_d_ocp_qp.h"
#include "hpipm/include/hpipm_d_ocp_qp_sol.h"
#include "hpipm/include/hpipm_d_dense_qp.h"
#include "hpipm/include/hpipm_d_dense_qp_sol.h"
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_blas.h"
#include "hpipm/include/hpipm_d_cond.h"

int ocp_qp_expand(const ocp_qp_in *qp_in, const dense_qp_out *sol, ocp_qp_out *qp_out, void *args_,
    void *memory_) {
        ocp_qp_condensing_memory *memory = (ocp_qp_condensing_memory *)memory_;

        int ii;

        int N = qp_in-> N;
        const int *nb = qp_in->nb;
        const int *ng = qp_in->nc;

        int nvd = sol-> nvd;
        int nbd = sol-> nbd;
        int ngd = sol-> ngd;
        struct d_ocp_qp *qp = memory->qp;
        struct d_ocp_qp_sol *qp_sol = memory->qp_sol;
        struct d_dense_qp_sol *qpd_sol = memory->qpd_sol;
        struct d_cond_qp_ocp2dense_workspace *cond_workspace = memory->cond_workspace;
        double **hlam_lb = memory->hlam_lb;
        double **hlam_ub = memory->hlam_ub;
        double **hlam_lg = memory->hlam_lg;
        double **hlam_ug = memory->hlam_ug;

        double **hx = qp_out->x;
        double **hu = qp_out->u;
        double **hpi = qp_out->pi;
        double **hlam = qp_out->lam;

        for (ii = 0; ii <= N; ii++) {
            hlam_lb[ii] = hlam[ii];
            hlam_ub[ii] = hlam[ii] + nb[ii];
            hlam_lg[ii] = hlam[ii] + 2 * nb[ii];
            hlam_ug[ii] = hlam[ii] + 2 * nb[ii] + ng[ii];
        }

        double *prim_sol = sol->u;
        double *dual_sol = sol->lam;          

        d_cvt_vec2strvec(nvd, prim_sol, qpd_sol->v, 0);
        for (ii=0; ii < 2*nbd+2*ngd; ii++) qpd_sol->lam->pa[ii] = 0.0;
        for (ii=0; ii < nbd; ii++)
            if (dual_sol[ii] >= 0.0)
                qpd_sol->lam->pa[ii] =   dual_sol[ii];
            else
                qpd_sol->lam->pa[nbd+ngd+ii] = - dual_sol[ii];
        for (ii=0; ii < ngd; ii++)
            if (dual_sol[nbd+ii] >= 0.0)
                qpd_sol->lam->pa[nbd+ii] =   dual_sol[nbd+ii];
            else
                qpd_sol->lam->pa[2*nbd+ngd+ii] = - dual_sol[nbd+ii];

        d_expand_sol_dense2ocp(qp, qpd_sol, qp_sol, cond_workspace);

        d_cvt_ocp_qp_sol_to_colmaj(qp, qp_sol, hu, hx, NULL, NULL, hpi,
            hlam_lb, hlam_ub, hlam_lg, hlam_ug, NULL, NULL);

        return 0;
    }

