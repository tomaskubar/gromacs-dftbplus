/*
 * Hybrid DFTB+ + PyTorch QMMM implementation (delta learning)
 */

#include "config.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/smalloc.h"

// When not built in a configuration with QMMM support, much of this
// code is unreachable by design. Tell clang not to warn about it.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-noreturn"

#if GMX_QMMM_DFTBPLUS_PYTORCH

#include "gromacs/mdlib/qm_dftbplus_pytorch.h"
#include <omp.h>

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
    real e_dftb = 0.0;
    real e_delta = 0.0;
    
    // Allocate memory for both force arrays
    rvec* f_dftb = nullptr;
    rvec* f_delta = nullptr;
    snew(f_dftb, nQMatoms);
    snew(f_delta, nQMatoms);
    
    for (int i=0; i < nQMatoms; i++)
    {
        clear_rvec(f_dftb[i]);
        clear_rvec(f_delta[i]);
    }
    
    #pragma omp parallel
    {
        #pragma omp single
        {
            // Task for DFTB+ calculation
            #pragma omp task shared(e_dftb, f_dftb)
            {
                // Call DFTB+ for base energy/forces (directly storing in f_dftb)
                e_dftb = call_dftbplus(qr, cr, qm, mm, f_dftb, fshift, nrnb, wcycle);
            }
            
            // Task for PyTorch calculation
            #pragma omp task shared(e_delta, f_delta)
            {
                // Call PyTorch for delta correction (using f_delta)
                e_delta = call_pytorch(qr, cr, qm, mm, f_delta, fshift, nrnb, wcycle);
            }
            
            // Wait for both tasks to complete
            #pragma omp taskwait
        }
    }
    
    // Combine the forces
    for (int i = 0; i < nQMatoms; ++i)
        for (int d = 0; d < DIM; ++d)
            f[i][d] = f_dftb[i][d] + f_delta[i][d];
    
    sfree(f_dftb);
    sfree(f_delta);
    
    // Return total energy (DFTB+ + delta)
    return e_dftb + e_delta;
}

#endif

#pragma GCC diagnostic pop

