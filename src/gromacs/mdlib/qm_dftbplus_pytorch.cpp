/*
 * Hybrid DFTB+ + PyTorch QMMM implementation (delta learning)
 */

#include "config.h"
#include "gromacs/utility/smalloc.h"

// When not built in a configuration with QMMM support, much of this
// code is unreachable by design. Tell clang not to warn about it.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-noreturn"

#if GMX_QMMM_DFTBPLUS_PYTORCH

#include "gromacs/mdlib/qm_dftbplus_pytorch.h"

void init_dftbplus_pytorch(QMMM_QMrec* qm, QMMM_rec* qr, const t_inputrec* ir, const t_commrec* cr)
{
    // Initialize DFTB+
    init_dftbplus(qm, qr, ir, cr);
    // Initialize PyTorch
    init_pytorch(qm);
}

real call_dftbplus_pytorch(QMMM_rec* qr,
                           const t_commrec* cr,
                           QMMM_QMrec* qm,
                           const QMMM_MMrec& mm,
                           rvec f[],
                           rvec fshift[],
                           t_nrnb* nrnb,
                           gmx_wallcycle_t wcycle)
{
    int nQMatoms = qm->nrQMatoms_get();
    // Call DFTB+ for base energy/forces
    real e_dftb = call_dftbplus(qr, cr, qm, mm, f, fshift, nrnb, wcycle);
    // Store DFTB+ forces for later correction
    rvec* f_dftb = nullptr;
    snew(f_dftb, nQMatoms);
    for (int i = 0; i < nQMatoms; ++i)
        for (int d = 0; d < DIM; ++d)
            f_dftb[i][d] = f[i][d];
    // Call PyTorch for delta correction
    real e_delta = call_pytorch(qr, cr, qm, mm, f, fshift, nrnb, wcycle);
    // Add delta correction to DFTB+ forces
    for (int i = 0; i < nQMatoms; ++i)
        for (int d = 0; d < DIM; ++d)
            f[i][d] += f_dftb[i][d];
    sfree(f_dftb);
    // Return total energy (DFTB+ + delta)
    return e_dftb + e_delta;
}

#endif

#pragma GCC diagnostic pop

