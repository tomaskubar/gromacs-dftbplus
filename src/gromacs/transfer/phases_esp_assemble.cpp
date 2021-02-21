#include "gromacs/transfer/transfer.h"


#if GMX_MPI

/*
void do_dftb_phase1(charge_transfer_t*                                ct,
                    std::vector<std::unique_ptr<QMMM_rec_transfer> >& dftbplus_phase1,
                    const t_commrec*                                  cr,
                    rvec                                              f[],
                    t_nrnb*                                           nrnb,
                    gmx_wallcycle_t                                   wcycle,
                    MPI_Comm                                          ct_mpi_comm,
                    int                                               ct_mpi_rank,
                    int                                               ct_mpi_size)
{
    (void) ct_mpi_comm;

    printf("do_dftb_phase1 start at rank %d at %f\n", ct_mpi_rank, (double) clock()/CLOCKS_PER_SEC);

    for (int i = 0; i < ct->sites; i++)
    {
        if (i % ct_mpi_size == ct_mpi_rank)
        {
         // printf("Doing residue %d at rank %d\n", ct->site[i].resnr, ct_mpi_rank);
         // run_dftb1(ct, dftb, i); TODO DFTB+
            call_dftbplus_transfer(dftbplus_phase1[i].get(), cr, f, nrnb, wcycle);
        }
    }

    printf("do_dftb_phase1 end at rank %d at %f\n", ct_mpi_rank, (double) clock()/CLOCKS_PER_SEC);
}
*/

void do_esp_only(charge_transfer_t*                                ct,
                 std::vector<std::unique_ptr<QMMM_rec_transfer> >& dftbplus_phase1,
                 const t_commrec*                                  cr,
                 real*                                             q,
                 t_nrnb*                                           nrnb,
                 gmx_wallcycle_t                                   wcycle,
                 MPI_Comm                                          ct_mpi_comm,
                 int                                               ct_mpi_rank,
                 int                                               ct_mpi_size)
{
    (void) ct;
    (void) dftbplus_phase1;
    (void) cr;
    (void) q;
    (void) nrnb;
    (void) wcycle;
    (void) ct_mpi_comm;
    (void) ct_mpi_rank;
    (void) ct_mpi_size;

    printf("ESP not implemented yet!\n");
    exit(-1);
/*
   if (ct->qmmm == 3) {
    for (int i=0; i<ct->sites; i++) {
      if (i % ct_mpi_size == ct_mpi_rank) {
        // set the charges on "QM" atoms
        for (int j=0; j<dftb->phase1[i].nn; j++) {
          dftb->phase1[i].qmat[j] = - q[ct->site[i].atom[j]] + dftb->qzero1[dftb->phase1[i].izp[j]];
        }
        dftb->phase1[i].qmat[0] += - q[ct->site[i].atom[0]+1];
        if (ct->modif_extcharge[i][0] > -1)
          dftb->phase1[i].qmat[0] -= - (ct->is_dna[i] ? EXTCHARGE_SHIFT_DNA : EXTCHARGE_SHIFT_TRP); //muss noch angepasst werden
        if (ct->modif_extcharge[i][1] > -1)
          dftb->phase1[i].qmat[0] -= - (ct->is_dna[i] ? EXTCHARGE_SHIFT_DNA : EXTCHARGE_SHIFT_TRP);
        // calculate the remaining components of PME
        do_pme_for_dftb_part2(ct, dftb, i);
        if (ct_mpi_rank != 0) {
          return_value = MPI_Send(&(dftb->phase1[i].esp), 1, MPI_DOUBLE, 0, 103 + i * 4, ct_mpi_comm);
        }
      }
    }
   }
 */
}

#else

/*
void do_dftb_phase1(charge_transfer_t*                                ct,
                    std::vector<std::unique_ptr<QMMM_rec_transfer> >& dftbplus_phase1,
                    const t_commrec*                                  cr,
                    rvec                                              f[],
                    t_nrnb*                                           nrnb,
                    gmx_wallcycle_t                                   wcycle)
{
    printf("do_dftb_phase1 start at %f\n", (double) clock()/CLOCKS_PER_SEC);

    for (int i = 0; i < ct->sites; i++)
    {
     // run_dftb1(ct, dftb, i); TODO DFTB+
        call_dftbplus_transfer(dftbplus_phase1[i].get(), cr, f, nrnb, wcycle);
    }

    printf("do_dftb_phase1 end   at %f\n", (double) clock()/CLOCKS_PER_SEC);
}
*/

void do_esp_only(charge_transfer_t*                                ct,
                 std::vector<std::unique_ptr<QMMM_rec_transfer> >& dftbplus_phase1,
                 const t_commrec*                                  cr,
                 real*                                             q,
                 t_nrnb*                                           nrnb,
                 gmx_wallcycle_t                                   wcycle)
{
    (void) ct;
    (void) dftbplus_phase1;
    (void) cr;
    (void) q;
    (void) nrnb;
    (void) wcycle;

    printf("ESP not implemented yet!\n");
    exit(-1);
/*
    if (ct->qmmm == 3)
    {
        for (int i=0; i<ct->sites; i++)
        {
            // set the charges on "QM" atoms
            for (int j=0; j<dftb->phase1[i].nn; j++)
            {
                dftb->phase1[i].qmat[j] = - q[ct->atom[i][j]] + dftb->qzero1[dftb->phase1[i].izp[j]];
            }
            dftb->phase1[i].qmat[0] += - q[ct->atom[i][0]+1];
            if (ct->modif_extcharge[i][0] > -1)
                dftb->phase1[i].qmat[0] -= - (ct->is_dna[i] ? EXTCHARGE_SHIFT_DNA : EXTCHARGE_SHIFT_TRP);
            if (ct->modif_extcharge[i][1] > -1)
                dftb->phase1[i].qmat[0] -= - (ct->is_dna[i] ? EXTCHARGE_SHIFT_DNA : EXTCHARGE_SHIFT_TRP);
            // calculate the remaining components of PME
            do_pme_for_dftb_part2(ct, dftb, i);
        }
    }
 */
}
#endif

/*
void do_dftb_phase2(charge_transfer_t*               ct,
                    std::unique_ptr<QMMM_rec_transfer>& dftbplus_phase2,
                    const t_commrec*                 cr,
                    rvec                             f[],
                    t_nrnb*                          nrnb,
                    gmx_wallcycle_t                  wcycle)
{
    (void) ct;

    printf("do_dftb_phase2 start at %f\n", (double) clock()/CLOCKS_PER_SEC);

 // run_dftb2(ct, dftb); TODO DFTB+
    call_dftbplus_transfer(dftbplus_phase2.get(), cr, f, nrnb, wcycle);

    printf("do_dftb_phase2 end   at %f\n", (double) clock()/CLOCKS_PER_SEC);
}
*/

void ct_assemble_hamiltonian(charge_transfer_t *ct) // , dftb_t *dftb)
{
    //  mean coupling for a and b direction in SAM.

    // save old hamilton matrices (they range from t=0 as newest martrix up to n_avg_ham steps into the past)
    // start with filled array from first step for meaningful averaging

    if (ct->jobtype == cteBCJFSSH)
    {
       for (int i=0; i<ct->dim; i++)
           ct->occupation[i] = SQR(ct->tfs_diab[i]) + SQR(ct->tfs_diab[i+ct->dim]);
    }
    else if (ct->jobtype == cteSCCDYNAMIC)
    {
       for (int i=0; i<ct->dim; i++)
           ct->occupation[i] = SQR(ct->wf[i]) + SQR(ct->wf[i+ct->dim]);
    }

    if (ct->first_step)
    {
        for (int t = 0; t < ct->n_avg_ham; t++)
        {
            for (int i = 0; i < ct->dim; i++)
            {
             // ct->hamiltonian_history[i][i][t]  = (ct->is_hole_transfer == 1) ? -dftb->phase2.THamilOrtho[i][i] : dftb->phase2.THamilOrtho[i][i];
                ct->hamiltonian_history[i][i][t]  = (ct->is_hole_transfer == 1) ? -ct->hamiltonian[i][i] : ct->hamiltonian[i][i];
                ct->hamiltonian_history[i][i][t] += ct->fo_shift[i];
                for (int j = i+1; j < ct->dim; j++)
                {
                    // 1.0 or 1.540 scaling factor according to J. Chem. Phys. 140, 104105 (2014)
                    // this is the right way to do it. without changing off-diagonal elements
                    //   the eigenvalues of the CG Ham are not necessarily -eigenvalue of electron transfer.
                    ct->hamiltonian_history[i][j][t] = ct->hamiltonian_history[j][i][t] =
                     // (ct->is_hole_transfer == 1) ? -dftb->phase2.THamilOrtho[i][j] * ct->offdiag_scaling : dftb->phase2.THamilOrtho[i][j] * ct->offdiag_scaling;
                        (ct->is_hole_transfer == 1) ? -ct->hamiltonian[i][j] * ct->offdiag_scaling : ct->hamiltonian[i][j] * ct->offdiag_scaling;
                }
            }
        }
    }
    else
    {
        for (int t = ct->n_avg_ham-1; t > 0; t--)
        {
            for (int i = 0; i < ct->dim; i++)
            {
                for (int j = 0; j < ct->dim; j++)
                {
                    ct->hamiltonian_history[i][j][t] = ct->hamiltonian_history[i][j][t-1];
                }
            }
        }

        for (int i = 0; i < ct->dim; i++)
        {
         // ct->hamiltonian_history[i][i][0]  = (ct->is_hole_transfer == 1) ? -dftb->phase2.THamilOrtho[i][i] : dftb->phase2.THamilOrtho[i][i];
            ct->hamiltonian_history[i][i][0]  = (ct->is_hole_transfer == 1) ? -ct->hamiltonian[i][i] : ct->hamiltonian[i][i];
            ct->hamiltonian_history[i][i][0] += ct->fo_shift[i];
            for (int j = i+1; j < ct->dim; j++)
            {
                ct->hamiltonian_history[i][j][0] = ct->hamiltonian_history[j][i][0] =
                    // 1.0 or 1.540 scaling factor according to J. Chem. Phys. 140, 104105 (2014)
                    // this is the right way to do it. without changing off-diagonal elements
                    // the eigenvalues of the CG Ham are not necessarily -eigenvalue of electron transfer.
                 // (ct->is_hole_transfer == 1) ? -dftb->phase2.THamilOrtho[i][j] * ct->offdiag_scaling :  dftb->phase2.THamilOrtho[i][j] * ct->offdiag_scaling;
                    (ct->is_hole_transfer == 1) ? -ct->hamiltonian[i][j] * ct->offdiag_scaling :  ct->hamiltonian[i][j] * ct->offdiag_scaling;
            }
        }
    }

    //build averaged hamiltonian
    for (int i = 0; i < ct->dim; i++)
    {
        for (int j = 0; j < ct->dim; j++)
        {
            ct->hamiltonian_dftb[i][j] = 0.0;
            for (int t = 0; t < ct->n_avg_ham; t++)
            {
                ct->hamiltonian_dftb[i][j] += ct->hamiltonian_history[i][j][t];
            }
            ct->hamiltonian_dftb[i][j] /= ct->n_avg_ham;
            if (j != i)
            {
                //test: only average diagonal elements. fast oscillations of off-diags can also be caused 
                //  by slow (classical) nuclear oscillation due to complex nodal structure of the orbitals.
                ct->hamiltonian_dftb[i][j] = ct->hamiltonian_history[i][j][0];
            }
            else if (j == i && ct->lambda != 0.0)
            {
                ct->hamiltonian_dftb[i][i] -= ct->occupation[i]*ct->lambda;
            }
        }
    }

    if (ct->freeze_coupling != 0.0)
    {
       for (int i = 0; i < ct->dim - 1; i++)
       {
           for (int j = i + 1; j < ct->dim; j++)
           {
               if (j == i + 1)
               {
                  ct->hamiltonian_dftb[i][j] = ct->freeze_coupling;
               }
               else
               {
                  ct->hamiltonian_dftb[i][j] = 0.0;
               }

               ct->hamiltonian_dftb[j][i] = ct->hamiltonian_dftb[i][j];
           }
       }
    }
}


void get_spectrum(charge_transfer_t *ct) //, dftb_t *dftb) //, FILE *f_ct_spec, FILE *f_ct_spec_evec){
{
 // dftb_phase2_t dftb2 = dftb->phase2;

    // diagonalize FO-Hamiltonian
    for (int i = 0; i < ct->dim; i++)
    {
        ct->evec_spec[i + i * ct->dim] = ct->hamiltonian[i][i];
        for (int j = 0; j < ct->dim; j++)
        {
            if (i != j)
            {
                ct->evec_spec[i + j * ct->dim] = ct->hamiltonian[i][j];
            }
        }
    }
    dsyev(ct->dim, ct->evec_spec, ct->ev_spec, ct->work_spec, 3*ct->dim);

    /* need larger work arrays
    // diagonalize AO-Hamiltonian
    for (int i=0; i<dftb2.norb; i++) {
        ct->evec_spec[i + i * dftb2.norb] = dftb2.hamil[i][i];
        for (int j=0; j<dftb2.norb; j++)
        {
            if (i!=j) ct->evec_spec[i + j * dftb2.norb] = dftb2.hamil[i][j];
        }
    }
    dsyev(dftb2.norb, ct->evec_spec, ct->ev_spec, ct->work_spec, 3*dftb2.norb);

    fprintf(f_ct_spec, "energy eigen values of AO Hamiltonian:");
    for (int i=0; i<dftb2.norb; i++)
        fprintf(f_ct_spec, " %10.6f \n", ct->ev_spec[i] * HARTREE_TO_EV);
    fprintf(f_ct_spec, "\n");
    */
}
