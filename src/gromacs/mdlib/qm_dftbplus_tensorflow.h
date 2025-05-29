/*
 * Hybrid DFTB+ + Tensorflow QMMM interface (delta learning)
 */
#ifndef GMX_MDLIB_QM_DFTBPLUS_TENSORFLOW_H
#define GMX_MDLIB_QM_DFTBPLUS_TENSORFLOW_H

#include "gromacs/mdlib/qmmm.h"
#include "gromacs/mdlib/qm_dftbplus.h"
#include "gromacs/mdlib/qm_tensorflow.h"
#include <vector>

void init_dftbplus_tensorflow(QMMM_QMrec* qm, QMMM_rec* qr, const t_inputrec* ir, const t_commrec* cr);

real call_dftbplus_tensorflow(QMMM_rec* qr,
                           const t_commrec* cr,
                           QMMM_QMrec* qm,
                           const QMMM_MMrec& mm,
                           rvec f[],
                           rvec fshift[],
                           t_nrnb* nrnb,
                           gmx_wallcycle_t wcycle);
                           
#endif