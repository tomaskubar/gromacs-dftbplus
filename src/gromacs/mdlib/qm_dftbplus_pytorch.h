/*
 * Hybrid DFTB+ + PyTorch QMMM interface (delta learning)
 */
#ifndef GMX_MDLIB_QM_DFTBPLUS_PYTORCH_H
#define GMX_MDLIB_QM_DFTBPLUS_PYTORCH_H

#include "gromacs/mdlib/qmmm.h"
#include "gromacs/mdlib/qm_dftbplus.h"
#include "gromacs/mdlib/qm_pytorch.h"
#include <vector>

void init_dftbplus_pytorch(QMMM_QMrec* qm, QMMM_rec* qr, const t_inputrec* ir, const t_commrec* cr);

real call_dftbplus_pytorch(QMMM_rec* qr,
                           const t_commrec* cr,
                           QMMM_QMrec* qm,
                           const QMMM_MMrec& mm,
                           rvec f[],
                           rvec fshift[],
                           t_nrnb* nrnb,
                           gmx_wallcycle_t wcycle);
                           
#endif
