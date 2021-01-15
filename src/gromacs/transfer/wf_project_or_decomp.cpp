#include "gromacs/transfer/transfer.h"

// needs fermi_coef_sum from quantum_auxiliary.cpp

int do_wf_decomp_evec(charge_transfer_t *ct, ct_diis_t *diis, FILE *f)
{
    int    i, j, k, step;
    double old_energy, energy, *ham, dotproduct_squared;

    ham                 = ct->hamiltonian_adiab; // why hamiltonian_adiab aren't we in diab basis dholub July 2018
    diis->n_prev_vector = 1;
    diis->indx          = 0;

    // calculate the charges from the wave function
    for (i = 0; i < ct->dim; i++)
    {
        ct->q_act[i] = SQR(ct->wf[i]);
    }

    // SCC cycle
    for (step = 0, energy = 1.e10;; step++)
    {

        if (step >= MAXITER_DIIS)
        {
            /* DO HERE SOMETHING, BUT DO NOT MODIFY ct->wf !
               for (i=0; i<ct->sites; i++)
               ct->wf[i] = sqrt(ct->q_act[i]);
             */
            fprintf(f, "diagonalization did not converge! ");
            return 1;
        }

        //printf("step %5d:", step);

        //printf(" old q =");
        //for (i=0; i<ct->sites; i++)
        //  printf("%9.5f", ct->q_act[i]);

        // construct the self-consistent Hamiltonian
        for (i = 0; i < ct->dim; i++)
        {
            ham[i + i * ct->dim] = ct->hamiltonian[i][i];
            for (j = 0; j < ct->dim; j++)
            {
                ham[i + i * ct->dim] += ct->hubbard[i][j] * ct->q_act[j];
                if (i != j)
                {
                    ham[i + j * ct->dim] = ct->hamiltonian[i][j];
                }
            }
        }

        // diagonalize this Hamiltonian
        dsyev(ct->dim, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);

        // copy previous charge and the difference to the arrays
        for (i = 0; i < ct->dim; i++)
        {
            diis->prev_q_input[diis->indx][i] = ct->q_act[i];
            diis->prev_q_diff[diis->indx][i]  = -ct->q_act[i]; // here, the new occupation will be added
        }
        // FERMI DISTRIBUTION (electronic temperature)
        // calculate the "Fermi energy", store the Fermi coefficients in diis->fermi_coef
        /* SWITCH IT OFF FOR DECOMPOSITION !
           fermi_lower = ct->ev_adiab[0] - 0.1;
           fermi_upper = ct->ev_adiab[ct->sites - 1] + 0.1;
           while (fermi_upper - fermi_lower > FERMI_CONVERG) {
           fermi_energy = 0.5 * (fermi_lower + fermi_upper);
           if (fermi_coef_sum(fermi_energy, ct->sites, ct->ev_adiab, diis->fermi_coef, ct->fermi_kt) > 1.0)
            fermi_upper = 0.5 * (fermi_lower + fermi_upper);
           else
            fermi_lower = 0.5 * (fermi_lower + fermi_upper);
           }
           fprintf(f, " fermi E = %9.6f, fermi coefs:", fermi_energy);
           for (i=0; i<ct->sites; i++)
           fprintf(f, " %8.5f", diis->fermi_coef[i]);
           for (i=0; i<ct->sites; i++)
           for (j=0; j<ct->sites; j++)
            diis->prev_q_diff[diis->indx][i] += diis->fermi_coef[j] * SQR(ham[j * ct->sites + i]);
         */
        for (i = 0; i < ct->dim; i++)
        {
            diis->prev_q_diff[diis->indx][i] += SQR(ham[i]);
        }

        old_energy = energy;
        /* AGAIN - DO NOT MIX THE ENERGY !
           energy = 0.e0;
           for (i=0; i<ct->sites; i++)
           energy += diis->fermi_coef[i] * ct->ev_adiab[i];
         */
        energy = ct->ev_adiab[0];
        //printf(" energy = %12.9f\n", energy);
        if (fabs(energy - old_energy) < BROYDEN_SCFTOL)
        {
            //printf("Converged!\n");
            /* NOW, WE HAVE THE ADIABATIC STATES IN ham !
               for (i=0; i<ct->sites; i++)
               ct->wf[i] = sqrt(diis->prev_q_input[diis->indx][i] + diis->prev_q_diff[diis->indx][i]);
             */
            break;
        }

        if (diis->n_prev_vector < DIIS_MAX_PREV_VECTORS)
        {
            // SIMPLE MIXING in a couple of first steps
            for (i = 0; i < ct->dim; i++)
            {
                ct->q_act[i] = diis->prev_q_input[diis->indx][i] + DIIS_INIT_MIXING_PARAM * diis->prev_q_diff[diis->indx][i];
            }
            diis->n_prev_vector++;
        }
        else
        {
            // REAL DIIS here
            for (i = 0; i < DIIS_MAX_PREV_VECTORS; i++)
            {
                for (j = 0; j < DIIS_MAX_PREV_VECTORS; j++)
                {
                    diis->aa[i + j * (DIIS_MAX_PREV_VECTORS + 1)] = 0.e0;
                    for (k = 0; k < ct->dim; k++)
                    {
                        diis->aa[i + j * (DIIS_MAX_PREV_VECTORS + 1)] += diis->prev_q_diff[i][k] * diis->prev_q_diff[j][k];
                    }
                }
            }
            for (i = 0; i < DIIS_MAX_PREV_VECTORS; i++)
            {
                diis->aa[i + DIIS_MAX_PREV_VECTORS * (DIIS_MAX_PREV_VECTORS + 1)] = diis->aa[DIIS_MAX_PREV_VECTORS + i * (DIIS_MAX_PREV_VECTORS + 1)] = -1.e0;
            }
            diis->aa[DIIS_MAX_PREV_VECTORS + DIIS_MAX_PREV_VECTORS * (DIIS_MAX_PREV_VECTORS + 1)] = 0.e0;

            for (i = 0; i < DIIS_MAX_PREV_VECTORS; i++)
            {
                diis->bb[i] = 0.e0;
            }
            diis->bb[DIIS_MAX_PREV_VECTORS] = -1.e0;

            dgesv(DIIS_MAX_PREV_VECTORS + 1, diis->aa, diis->bb, diis->ipiv);

            for (j = 0; j < ct->dim; j++)
            {
                ct->q_act[j] = 0.e0;
                for (i = 0; i < DIIS_MAX_PREV_VECTORS; i++)
                {
                    ct->q_act[j] += diis->bb[i] * diis->prev_q_input[i][j];
                }
            }
        }
        diis->indx = (diis->indx + 1) % DIIS_MAX_PREV_VECTORS;
    }

    // print energies of adiabatic states
    //printf("Adiab E: ");
    for (i = 0; i < ct->dim; i++)
    {
        //fprintf(f, "%10.7f", ct->ev_adiab[i]);
        fprintf(f, "%10.7f", ct->ev_adiab[i] - ct->ev_adiab[0]);
    }

    // calculate the expansion coefficients
    fprintf(f, "   WF exp'd in adiab states: ");
    /*
       for (i=0; i<ct->sites; i++) {
       // real part of a_i
       dotproduct_re = 0.e0;
       for (j=0; j<ct->sites; j++) dotproduct_re += ham[i * ct->sites + j] * ct->wf[j];
       // imaginary part of a_i
       dotproduct_im = 0.e0;
       for (j=0; j<ct->sites; j++) dotproduct_im += ham[i * ct->sites + j] * ct->wf[j + ct->sites];
       fprintf(f, "%9.6f", SQR(dotproduct_re)+SQR(dotproduct_im));
       }
     */
    for (i = 0; i < ct->dim; i++)
    {
        // calculate directly the square of a_i
        dotproduct_squared = 0.e0;
        for (j = 0; j < ct->dim; j++)
        {
            for (k = 0; k < ct->dim; k++)
            {
                dotproduct_squared += ham[i * ct->dim + j] * ham[i * ct->dim + k] * (ct->wf[j] * ct->wf[k] + ct->wf[j + ct->dim] * ct->wf[k + ct->dim]);
            }
        }
        fprintf(f, "%9.6f", dotproduct_squared);
    }
    fprintf(f, "\n");

    return 0;
}


int do_wf_decomp_evec_fermi(charge_transfer_t *ct, dftb_broyden_t *broyd, double *fermi_coeff, FILE *f)
{
    int    i, j, k, step;
    double old_energy, energy, *ham, fermi_energy, fermi_upper, fermi_lower, dotproduct_real, dotproduct_imag, q_sum;

    ham = ct->hamiltonian_adiab;

    // calculate the charges from the wave function
    for (i = 0; i < ct->dim; i++)
    {
        ct->q_act[i] = ct->q_old[i] = SQR(ct->wf[i]);
    }

    energy = 1.e10;

    // SCC cycle
    for (step = 0;; step++)
    {

        // too many iterations? break and indicate failure
        if (step >= MAXITER_BROYDEN)
        {
            fprintf(f, "Broyden/Fermi failed to converge!\n");
            return 1;
        }

        // construct the self-consistent Hamiltonian
        for (i = 0; i < ct->dim; i++)
        {
            ham[i + i * ct->dim] = ct->hamiltonian[i][i];
            for (j = 0; j < ct->dim; j++)
            {
                ham[i + i * ct->dim] += ct->hubbard[i][j] * ct->q_act[j];
                if (i != j)
                {
                    ham[i + j * ct->dim] = ct->hamiltonian[i][j];
                }
            }
        }

        // diagonalize this Hamiltonian
        dsyev(ct->dim, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);

        // FERMI DISTRIBUTION + BROYDEN MIXING

        // calculate the "Fermi energy", store the Fermi coefficients in fermi_coeff
        fermi_lower = ct->ev_adiab[0] - 0.1;
        fermi_upper = ct->ev_adiab[ct->dim - 1] + 0.1;
        while (fermi_upper - fermi_lower > FERMI_CONVERG)
        {
            fermi_energy = 0.5 * (fermi_lower + fermi_upper);
            //if (fermi_coef_sum(fermi_energy, ct->dim, ct->ev_adiab, fermi_coeff, ct->fermi_kt) > 1.0)
            if (fermi_coef_sum(fermi_energy, ct->dim, ct->ev_adiab, fermi_coeff, FERMI_KT) > 1.0) // ct->fermi_kt is not read for SCC
            {
                fermi_upper = 0.5 * (fermi_lower + fermi_upper);
            }
            else
            {
                fermi_lower = 0.5 * (fermi_lower + fermi_upper);
            }
        }

        // mix the eigenvectors, with the appropriate coefficients
        for (i = 0; i < ct->dim; i++) // i runs over the sites
        {
            ct->q_act[i] = 0.e0;
            for (j = 0; j < ct->dim; j++) // j runs over the individual eigenvectors
            {
                ct->q_act[i] += fermi_coeff[j] * SQR(ham[j * ct->dim + i]);
            }
        }

        // calculate the energy and check convergence
        old_energy = energy;
        energy     = 0.e0;
        // indices i and j run over the sites
        for (i = 0; i < ct->dim; i++)
        {
            for (j = 0; j < ct->dim; j++)
            {
                // TB Hamiltonian
                // index k runs over the adiabatic states (eigenvectors)
                for (k = 0; k < ct->dim; k++)
                {
                    energy += fermi_coeff[k] * ham[k * ct->dim + i] * ham[k * ct->dim + j] * ct->hamiltonian[i][j];
                }
                // Hubbard / gamma terms
                energy += 0.5 * ct->hubbard[i][j] * ct->q_act[i] * ct->q_act[j];
            }
        }

        if (step > 1 && fabs(energy - old_energy) < BROYDEN_SCFTOL)
        {
            // Broyden mixing has converged
            break;
        }

        // mix the charges - Broyden
        broyden(step, BROYDEN_ALMIX, ct->dim, ct->q_old, ct->q_act, broyd);
        for (i = 0; i < ct->dim; i++)
        {
            ct->q_act[i] = ct->q_old[i];
        }
        // correct numerical errors caused by broyden
        q_sum = 0;
        for (i = 0; i < ct->dim; i++)
        {
            if (ct->q_act[i] < 0.0)
            {
                ct->q_act[i] = 0.0;
            }
            if (ct->q_act[i] > 1.0)
            {
                ct->q_act[i] = 1.0;
            }
            q_sum += ct->q_act[i];
        }
        for (i = 0; i < ct->dim; i++)
        {
            ct->q_act[i] = ct->q_act[i] / q_sum; //normalization
        }
    } /* end SCC cycle */

    // print the energies of "adiabatic states"
    for (i = 0; i < ct->dim; i++)
    {
        fprintf(f, "%10.7f", ct->ev_adiab[i]);
    }

    // calculate the expansion coefficients
    fprintf(f, " WF expansion:");
    for (i = 0; i < ct->dim; i++)
    {
        dotproduct_real = 0.e0;
        for (j = 0; j < ct->dim; j++)
        {
            dotproduct_real += ham[i * ct->dim + j] * ct->wf[j];
        }
        dotproduct_imag = 0.e0;
        for (j = 0; j < ct->dim; j++)
        {
            dotproduct_imag += ham[i * ct->dim + j] * ct->wf[j + ct->dim];
        }
        fprintf(f, "%9.6f", SQR(dotproduct_real)+SQR(dotproduct_imag));
    }
    /*
       for (i=0; i<ct->sites; i++) {
       dotproduct_real = dotproduct_imag = 0.e0;
       for (j=0; j<ct->sites; j++)
        for (k=0; k<ct->sites; k++)
          dotproduct_squared += ham[i * ct->sites + j] * ham[i * ct->sites + k] * (ct->wf[j] * ct->wf[k] + ct->wf[j + ct->sites] * ct->wf[k + ct->sites]);
       fprintf(f, "%9.6f", dotproduct_squared);
       }
     */
    fprintf(f, "\n");

    return 0;
}


void project_wf_on_new_basis_exact(int step, dftb_t *dftb, charge_transfer_t *ct, FILE *f_ct_project_wf, FILE *f_ct_project_wf_ref)
//void project_wf_on_new_basis(int step, dftb_t *dftb, charge_transfer_t *ct, FILE *f_ct_project_wf, FILE *f_ct_project_wf_ref)
{
    int    i, ii, j, k, l, m; // jj, n;
    int    iao, jao, ifo; // offset, jfo;
 // double sum;
    double norm;
// Transform orthogonal eigenvectors (basis of FO Hamiltonian) from non-orthogonal fragment basis to atomic basis //
    for (l = 0; l < ct->dim; l++)
    {
        for (iao = 0; iao < dftb->phase2.norb; iao++)
        {
            dftb->orthogo.evec_ao[iao][l] = 0.0;
        }
    }
    for (l = 0; l < ct->dim; l++)
    {
        k = 0;
        for (i = 0; i < ct->sites; i++)
        {
            for (ii = 0; ii < ct->site[i].homos; ii++)
            {
                ifo = ct->site[i].homo[ii] + dftb->phase2.inf[i] - 1;
                for (iao = 0; iao < dftb->phase2.norb; iao++)
                {
                    dftb->orthogo.evec_ao[iao][l] += dftb->orthogo.sij[k + l * ct->dim] * dftb->phase2.Taf[iao][ifo];
                }
                k++;
            }
        }
    }
    if (ct->first_step)
    {
        printf("initialize state following \n");
        for (l = 0; l < ct->dim; l++)
        {
            for (iao = 0; iao < dftb->phase2.norb; iao++)
            {
                dftb->orthogo.evec_ao_old[iao][l] = dftb->orthogo.evec_ao[iao][l];
                dftb->orthogo.evec_ao_ref[iao][l] = dftb->orthogo.evec_ao[iao][l];
            }
        }
    }
/*
   // calculate overlap for basis function i
   // using combined overlap matrix with different atomic overlaps for intra and inter-fragment calculations
   //no longer needed. now dftb2.overl is already hybrid matrix
   // <evec_ao|evec_ao_old> = evec_ao^T * S_ao * evec_ao_old
   for (iao=0; iao<dftb->phase2.norb; iao++)
   for (jao=0; jao<dftb->phase2.norb; jao++)
   dftb->phase2.overl_hybrid[iao][jao]=dftb->phase2.overl[iao][jao]; // we just need here the inter fragment overlap. intra fragment overlap will be overwritten
   offset=0;
   for (i=0; i<ct->sites; i++){
   for (j=0; j < dftb->phase1[i].norb; j++) {
   for (k=0; k < dftb->phase1[i].norb; k++){
   dftb->phase2.overl_hybrid[offset+j][offset+k] = dftb->phase1[i].overl[j][k];
   dftb->phase2.overl_hybrid[offset+k][offset+j] = dftb->phase1[i].overl[k][j];
   }
   }
   offset+=dftb->phase1[i].norb;
   }
 */
/*
   offset=0;
   for (i=0; i<ct->sites; i++){
   m=0;
   for (j=0; j < dftb->phase1[i].norb; j++) {
   n=0;
   for (k=0; k < dftb->phase1[i].norb; k++){
   dftb->phase2.overl_hybrid[offset+m][offset+n] = dftb->phase1[i].overl[j][k];
   dftb->phase2.overl_hybrid[offset+n][offset+m] = dftb->phase1[i].overl[k][j];
   n++;
   }
   m++;
   }
   offset+=dftb->phase1[i].norb;
   }
 */
/*
   printf("overlap ortho/ortho (hybrid matrix)\n");
   printf( "%10d ", step);
   for (i=0; i<ct->dim; i++) {
   for (m=i; m<ct->dim; m++) printf(" %10.6f ", dftb->overl_test[i][m]);
   }
   printf( "\n");
 */
// overlap with last step //
// <evec_ao|evec_ao_old> = evec_ao^T * S_ao * evec_ao_old
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            dftb->orthogo.overlap[i][j] = 0.0;
            for (iao = 0; iao < dftb->phase2.norb; iao++)
            {
                for (jao = 0; jao < dftb->phase2.norb; jao++)
                {
                    dftb->orthogo.overlap[i][j] += dftb->orthogo.evec_ao[iao][i] * dftb->phase2.overl[iao][jao] * dftb->orthogo.evec_ao_old[jao][j];
                }
            }
        }
    }
// overlap with wf in t=0
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            dftb->orthogo.overlap_ref[i][j] = 0.0;
            for (iao = 0; iao < dftb->phase2.norb; iao++)
            {
                for (jao = 0; jao < dftb->phase2.norb; jao++)
                {
                    dftb->orthogo.overlap_ref[i][j] += dftb->orthogo.evec_ao[iao][i] * dftb->phase2.overl[iao][jao] * dftb->orthogo.evec_ao_ref[jao][j];
                }
            }
        }
    }
///*
    printf("Wave function before proj %d:\n", step);
    for (i = 0; i < ct->dim; i++)
    {
        printf("%10.7f %10.7f\n", ct->wf[i], ct->wf[i + ct->dim]);
    }
//*/
// get new wf by projecting old wf onto new basis. element i is now:
// wf_i(t2) = |fo_i(t2)><fo_i(t2)|wf(t1)>
// wf_i(t2) = sum_j |fo_i(t2)><fo_i(t2)|fo_j(t1)>*c_j
    for (i = 0; i < ct->dim; i++)
    {
        ct->wf_old[i]         = ct->wf[i];
        ct->wf_old[i+ct->dim] = ct->wf[i+ct->dim];
        ct->wf[i]             = 0.0;
        ct->wf[i+ct->dim]     = 0.0;
    }
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            ct->wf[i]         += ct->wf_old[j] * dftb->orthogo.overlap[i][j];         //real part
            ct->wf[i+ct->dim] += ct->wf_old[j+ct->dim] * dftb->orthogo.overlap[i][j]; //imaginary part
        }
    }
// scale new wavefunction (reasons: projection is not complete) //
    norm = 0;
    for (i = 0; i < ct->dim; i++)
    {
        norm += SQR(ct->wf[i]) + SQR(ct->wf[i+ct->dim]);
    }
    for (i = 0; i < 2*ct->dim; i++)
    {
        ct->wf[i] /= sqrt(norm);
    }
///*
    printf("Wave function after proj %d:\n", step);
    for (i = 0; i < ct->dim; i++)
    {
        printf("%10.7f %10.7f\n", ct->wf[i], ct->wf[i + ct->dim]);
    }
//*/
// save wave function for next step //
    for (l = 0; l < ct->dim; l++)
    {
        for (iao = 0; iao < dftb->phase2.norb; iao++)
        {
            dftb->orthogo.evec_ao_old[iao][l] = dftb->orthogo.evec_ao[iao][l];
        }
    }
// OUTPUT INTO FILES//
    fprintf(f_ct_project_wf, "%10d ", step);
    for (i = 0; i < ct->dim; i++)
    {
        for (m = i; m < ct->dim; m++)
        {
            fprintf(f_ct_project_wf, " %10.6f ", dftb->orthogo.overlap[i][m]);
        }
    }
    fprintf(f_ct_project_wf, " %10.6f ", sqrt(norm));
    fprintf(f_ct_project_wf, "\n");
    fprintf(f_ct_project_wf_ref, "%10d ", step);
    for (i = 0; i < ct->dim; i++)
    {
        for (m = i; m < ct->dim; m++)
        {
            fprintf(f_ct_project_wf_ref, " %10.6f ", dftb->orthogo.overlap_ref[i][m]);
        }
    }
    fprintf(f_ct_project_wf_ref, "\n");
    return;
}


//void project_wf_on_new_basis_approx(int step, dftb_t *dftb, charge_transfer_t *ct, FILE *f_ct_project_wf, FILE *f_ct_project_wf_ref )
void project_wf_on_new_basis(int step, dftb_t *dftb, charge_transfer_t *ct, FILE *f_ct_project_wf, FILE *f_ct_project_wf_ref )
{
    int    i, j, k, l, m; // ii, jj, n;
    int    iao, jao, ifo, jfo; // offset;
 // double sum;
    double norm;

    // get new wf by projecting old wf onto new basis. element i will be:
    // wf_i(t2) =       |fo_i(t2)><fo_i(t2)|wf(t1)>
    // wf_i(t2) = sum_j |fo_i(t2)><fo_i(t2)|fo_j(t1)>*c_j

    // exact way would be transformation of orthogonalized basis functions |fo_i> into AO basis. Then calculate overlap with last step as:  <fo_i_ao|fo_i_ao_old> = fo_i_ao^T * S_ao * fo_i_ao_old
    // fast approximate version: use non-orthogonal FOs instead of orthogonalized ones to calculate <fo_i_ao|fo_i_ao_old>. -> is blockdiagonal matrix -> linear scaling
    // furthermore <fo_i_ao|fo_i_ao_old> was already calculated by check_and_invert_orbital_phase().



    // save reference WF at t=0 to get the total change of the basis
    if (ct->first_step)
    {
        for (i = 0; i < ct->sites; i++)
        {
            for (iao = 0; iao < dftb->phase1[i].norb; iao++)
            {
                for (jao = 0; jao < dftb->phase1[i].norb; jao++)
                {
                    dftb->phase1[i].a_ref[iao][jao] = dftb->phase1[i].a[iao][jao];
                }
            }
        }
    }

    // overlap with reference wf in t=0
    for (k = 0; k < ct->sites; k++)
    {
        for (l = 0; l < ct->site[k].homos; l++)
        {
            for (m = 0; m < ct->site[k].homos; m++)
            {
                i = dftb->phase2.ihomo[k]+l;
                j = dftb->phase2.ihomo[k]+m;
                dftb->orthogo.overlap_ref[i][j] = 0.0;
                for (iao = 0; iao < dftb->phase1[k].norb; iao++)
                {
                    for (jao = 0; jao < dftb->phase1[k].norb; jao++)
                    {
                        ifo = ct->site[k].homo[l]-1;
                        jfo = ct->site[k].homo[m]-1;
                        dftb->orthogo.overlap_ref[i][j] += dftb->phase1[k].a[iao][ifo] * dftb->phase1[k].overl[iao][jao] * dftb->phase1[k].a_ref[jao][jfo];
                    }
                }
            }
        }
    }


    // calculate overlap with last step //
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            dftb->orthogo.overlap[i][j] = 0.0;
        }
    }
    for (k = 0; k < ct->sites; k++)
    {
        for (l = 0; l < ct->site[k].homos; l++)
        {
            for (m = 0; m < ct->site[k].homos; m++)
            {
                i = dftb->phase2.ihomo[k]+l;
                j = dftb->phase2.ihomo[k]+m;
                dftb->orthogo.overlap[i][j] = ct->site[k].overlap[l][m];
            }
        }
    }

/*
         printf("Wave function before proj %d:\n", step);
         for (i=0; i<ct->dim; i++)
           printf("%10.7f %10.7f\n", ct->wf[i], ct->wf[i + ct->dim]);
 */


    for (i = 0; i < ct->dim; i++)
    {
        ct->wf_old[i]         = ct->wf[i];
        ct->wf_old[i+ct->dim] = ct->wf[i+ct->dim];
        ct->wf[i]             = 0.0;
        ct->wf[i+ct->dim]     = 0.0;
    }
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            ct->wf[i]         += ct->wf_old[j] * dftb->orthogo.overlap[i][j];         //real part
            ct->wf[i+ct->dim] += ct->wf_old[j+ct->dim] * dftb->orthogo.overlap[i][j]; //imaginary part
        }
    }


    // scale new wavefunction (reasons: projection is not complete) //
    norm = 0;
    for (i = 0; i < ct->dim; i++)
    {
        norm += SQR(ct->wf[i]) + SQR(ct->wf[i+ct->dim]);
    }
    printf("norm changed to %f due to projection. Rescaling wavefunction\n", norm);
    for (i = 0; i < 2*ct->dim; i++)
    {
        ct->wf[i] /= sqrt(norm);
    }


/*
         printf("Wave function  after proj %d:\n", step);
         for (i=0; i<ct->dim; i++)
           printf("%10.7f %10.7f\n", ct->wf[i], ct->wf[i + ct->dim]);
 */



    // OUTPUT INTO FILES//
    fprintf(f_ct_project_wf, "%10d ", step);
    for (i = 0; i < ct->dim; i++)
    {
        for (m = i; m < ct->dim; m++)
        {
            fprintf(f_ct_project_wf, " %10.6f ", dftb->orthogo.overlap[i][m]);
        }
    }
    fprintf(f_ct_project_wf, " %10.6f ", sqrt(norm));
    fprintf(f_ct_project_wf, "\n");

    fprintf(f_ct_project_wf_ref, "%10d ", step);
    for (i = 0; i < ct->dim; i++)
    {
        for (m = i; m < ct->dim; m++)
        {
            fprintf(f_ct_project_wf_ref, " %10.6f ", dftb->orthogo.overlap_ref[i][m]);
        }
    }
    fprintf(f_ct_project_wf_ref, "\n");


    return;
}

