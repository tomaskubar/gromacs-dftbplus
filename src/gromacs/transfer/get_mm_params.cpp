#include "transfer.h"

#if GMX_MPI

void get_MM_params(charge_transfer_t *ct, dftb_t *dftb, MPI_Comm ct_mpi_comm, int ct_mpi_rank, int ct_mpi_size)
{
    int i;
    for (i = 0; i < ct->sites; i++)
    {
        if (i % ct_mpi_size == ct_mpi_rank)
        {
            if (ct->delta_q_mode == 0)
            {
                get_delta_q(dftb, ct, i);
            }
            if (!ct->hamiltonian_type && (ct->do_lambda_i == 2 || ct->do_lambda_i == 3))
            {
                get_internal_forces(dftb, ct, i);
            }
        }
    }
    printf("MM params part 1 end   at %f\n", (double) clock()/CLOCKS_PER_SEC);
    if (!ct->hamiltonian_type && (ct->do_lambda_i == 3 && ct_mpi_rank == 0))
    {
       offdiag_gradient_homo(dftb, dftb->phase2.x, dftb->phase2.grad, ct);
       //for (i=0; i<dftb->phase2.nn; i++)
       //    printf("offdiag grad at atom i %d  QM force %lf \n", i,  -(real) HARTREE_BOHR2MD * fabs(dftb->phase2.grad[i][0])+fabs(dftb->phase2.grad[i][1])+fabs(dftb->phase2.grad[i][2]));
       printf("MM params part 2 end   at %f\n", (double) clock()/CLOCKS_PER_SEC);
    }

    (void) ct_mpi_comm;
    return;
}

#else

void get_MM_params(charge_transfer_t *ct, dftb_t *dftb)
{
    int i;
    for (i = 0; i < ct->sites; i++)
    {
        if (ct->delta_q_mode == 0)
        {
            get_delta_q(dftb, ct, i);
        }
        if (!ct->hamiltonian_type && (ct->do_lambda_i == 2 || ct->do_lambda_i == 3))
        {
            get_internal_forces(dftb, ct, i);
        }
    }

    if (ct->do_lambda_i == 3 && !ct->hamiltonian_type)
    {
       if (((ct->step % JFSSH_OFFDIAG_FORCE == 0 && ct->jobtype == cteJFSSH) || ct->jobtype == cteBCJFSSH) ||
           (ct->step % (JFSSH_OFFDIAG_FORCE/10) == 0 && ct->jobtype == cteSCCDYNAMIC) ||
           (dftb->do_ncv == 1))
       {
           offdiag_gradient_homo(dftb, dftb->phase2.x, dftb->phase2.grad, ct);
           //for (i=0; i<dftb->phase2.nn; i++)
           //    printf("offdiag grad at atom i %d  QM force %lf \n", i,  -(real) HARTREE_BOHR2MD * fabs(dftb->phase2.grad[i][0])+fabs(dftb->phase2.grad[i][1])+fabs(dftb->phase2.grad[i][2]));
           printf("offdiag forces end   at %f\n", (double) clock()/CLOCKS_PER_SEC);
       }
    }

    return;
}

#endif


void get_delta_q(dftb_t *dftb, charge_transfer_t *ct, int i)
{
    int    j, k, lj, m, n, jofn, mj, mo;
    double q, qhelp;
    //calculates atomic charge differences for hole in MO i

//for (i = 0; i < ct->sites; i++)
    for (k = 0; k < ct->site[i].homos; k++)
    {
        mo = ct->site[i].homo[k]-1;
//printf("site %d homo %d\n", i, mo);

        for (m = 0; m < dftb->phase1[i].ndim; m++)
        {
            dftb->phase1[i].qmulli[m] = 0.0;
            for (n = 0; n < dftb->phase1[i].ndim; n++)
            {
//      dftb->phase1[i].qmulli[m] += dftb->phase1[i].overl[m][n] * dftb->phase1[i].occ[mo] * dftb->phase1[i].a[m][mo] * dftb->phase1[i].a[n][mo]; //no sum over occupied MOs, only hole-MO needed for delta_q. Attention: overwrites qmulli
                dftb->phase1[i].qmulli[m] += dftb->phase1[i].overl[m][n] * dftb->phase1[i].a[m][mo] * dftb->phase1[i].a[n][mo]; //no sum over occupied MOs, only hole-MO needed for delta_q. Attention: overwrites qmulli
            }
        }

        for (j = 0; j < dftb->phase1[i].nn; j++)
        {
            q = 0.0;
            // printf("j = %d\n", j);
            for (lj = 0; lj < dftb->lmax[dftb->phase1[i].izp[j]]; lj++)
            {
                jofn = dftb->phase1[i].ind[j] + lj*lj;
                // printf("lj = %d, jofn = %d\n", lj, jofn);
                qhelp = 0.0;
                for (mj = 0; mj <= 2*lj; mj++)
                {
                    // printf("mj = %d\n", mj);
                    qhelp += dftb->phase1[i].qmulli[jofn+mj];
                }
                q += qhelp;
            }
            ct->site[i].delta_q[k][j] = (ct->is_hole_transfer == 1) ? q : -q;
            //printf("mulliken delta_q %f \n", ct->site[i].delta_q[k][j]);
        }
    }
    return;
}


void get_internal_forces(dftb_t *dftb, charge_transfer_t *ct, int site_i)
{
    int           i, j, nn;
 // int           l;
    dftb_phase1_t dftb1;

    dftb1 = dftb->phase1[site_i];
    nn    = dftb1.nn;
    //sign convention: calculate change of gradient in charged system due to additional charge carrier.
    // grad_total = grad_charged-grad_neutral    F= -grad

    for (i = 0; i < nn; i++)
    {
        clear_dvec(dftb1.grad[i]);
    }
    // calculate difference of atomic hamilton shifts (= sum over gamma*delta_charge)
    // no shift due to external field. forces due to electrostatic interacions of charge carrier and MM atoms are handled by Gromacs.
    for (i = 0; i < nn; i++)
    {
        dftb1.shift[i] = 0.0;
        for (j = 0; j < nn; j++)
        {
            dftb1.shift[i] += (dftb1.qmat[j] - dftb->qzero1[dftb1.izp[j]]) * (i > j ? dftb1.gammamat[i][j] : dftb1.gammamat[j][i]);
        }
        //for(l=0; l<ct->site[site_i].homos; l++)
        //dftb1.shift[i] += -ct->site[site_i].delta_q[l][j] * ct->occupation[ct->indFO[site_i]+l] * (i>j ? dftb1.gammamat[i][j] : dftb1.gammamat[j][i]);
//printf("shift %d = %f dQ= %f\n", i, dftb1.shift[i], (dftb1.qmat[i] - dftb->qzero1[dftb1.izp[i]]) );
    }

    //force due to missing electron in HOMO
    for (i = 0; i < nn; i++)
    {
        clear_dvec(dftb1.partgrad[i]);
    }

    usual_gradient_homo(dftb, dftb1.x, dftb1.partgrad, ct, site_i);

    for (i = 0; i < nn; i++)
    {
        copy_dvec(dftb1.partgrad[i], dftb1.grad[i]);
    }

/*
   for (i=0; i<nn; i++)
    for(l=0; l<DIM; l++)
      printf("usual grad at atom i %d  QM force %lf \n", i,  -(real) HARTREE_BOHR2MD * dftb1.partgrad[i][l]);
 */
/*
   //force due to changing atomic charges by extracting electron from HOMO // now we are only calculating DFTB1 forces
   for (i=0; i<nn; i++)
    clear_dvec(dftb1.partgrad[i]);
   gamma_gradient_homo(dftb, dftb1.x, dftb1.partgrad, ct, site_i);
   for (i=0; i<nn; i++)
    dvec_inc(dftb1.grad[i], dftb1.partgrad[i]);

   for (i=0; i<nn; i++)
    for(l=0; l<DIM; l++)
      printf("gamma grad at atom i %d  QM force %lf \n", i,  -(real) HARTREE_BOHR2MD * dftb1.partgrad[i][l]);
 */

/*
   // third force term //this is not needed
   for (i=0; i<nn; i++)
    clear_dvec(dftb1.partgrad[i]);
   //  additional_gradient_homo(dftb, dftb1.x, dftb1.partgrad, ct, site_i);
   // for (i=0; i<nn; i++)
   //    dvec_inc(dftb1.grad[i], dftb1.partgrad[i]);
   for (i=0; i<nn; i++)
    for(l=0; l<DIM; l++)
      printf("additional grad at atom i %d  QM force %lf \n", i,  -(real) HARTREE_BOHR2MD * dftb1.partgrad[i][l]);
 */

    return;
}

