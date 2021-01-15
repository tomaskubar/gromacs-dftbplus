//#include "gromacs/domdec/collect.h"
//#include "gromacs/domdec/domdec_struct.h"
//#include "gromacs/mdtypes/commrec.h"
//#include "gromacs/mdtypes/state.h"

#include "gromacs/transfer/transfer.h"

/* Current routine, adopted from mdoutf_write_to_trajectory_files() in mdoutf.cpp */    
    
void transfer_collect_x(const t_commrec* cr,
                        t_state*         state_local,
                        t_state*         state_global)
{
    if (DOMAINDECOMP(cr))
    {
        auto globalXRef = MASTER(cr) ? state_global->x : gmx::ArrayRef<gmx::RVec>();
        dd_collect_vec(cr->dd,
                       state_local->ddp_count,
                       state_local->ddp_count_cg_gl,
                       state_local->cg_gl,
                       state_local->x,
                       globalXRef);
    }
}

/* The now obsolete routine for particle decomposition */

#define GMX_LEFT (0)
#define GMX_RIGHT (1)
static void moveit(t_commrec*  cr,
                   int         left,
                   int         right,
                   const char* s,
                   rvec        xx[])
{
    (void) cr;
    (void) left;
    (void) right;
    (void) s;
    (void) xx;
}

void ct_collect_x(t_commrec* cr,
                  t_state*   state_global);
void ct_collect_x(t_commrec* cr,
                  t_state*   state_global)
{
#define MX(xvf) moveit(cr, GMX_LEFT, GMX_RIGHT, #xvf, xvf)

    if (cr->nnodes > 1)
    {
        /* Particle decomposition, collect the data on the master node */
        MX(state_global->x.rvec_array());
    }
}
