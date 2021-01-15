#include "gromacs/transfer/transfer.h"

extern "C" void dpotrf_(char *, long *, double *, long *, long *);
extern "C" void dpotri_(char *, long *, double *, long *, long *);
extern "C" void dgetrf_(long *, long *, double *, long *, long *, long *);
extern "C" void dgetri_(long *, double *, long *, long *, double *, long *, long *);

double calc_tda(charge_transfer_t *ct)
{
    long   nb, info, lwork = 1000, ipiv[1000];
    int    i, j;
    char   itype = 'U';
    double etun, etun_old, tda;
    double vdb[1000], vba[1000], temp[1000]; //should be allocatable, no. of bridge MOs
    double gb[1000*1000], work[1000];
    double ev_Heff[2], Heff[4];              //2x2
    int    niter, maxiter = 30;

    double Etol = 1.e-7;

    printf("Performing TDA Calculation\n");

    // defining number of bridge units.
    // these are all orbitals besides donor and acceptor HOMO (multiple FOs on bridge molecules possible)
    nb = (long) ct->dim-2;


    // define tunneling energy
    etun = 0.5 * (ct->hamiltonian[0][0] + ct->hamiltonian[ct->dim-1][ct->dim-1]);

    // define donor-bridge coupling
    // define bridge-acceptor coupling
    for (i = 0; i < nb; i++)
    {
        vdb[i] = ct->hamiltonian[0][i+1];
        vba[i] = ct->hamiltonian[ct->dim-1][i+1];
    }


    //start self consistent tuning of Etun
    for (niter = 0; niter < maxiter; niter++)
    {

        // build bridge green function
        for (i = 0; i < nb; i++)
        {
            for (j = 0; j < nb; j++)
            {
                if (i == j)
                {
                    gb[i+j*nb] = etun - ct->hamiltonian[i+1][j+1];
                }
                else
                {
                    gb[i+j*nb] = -1.0 * ct->hamiltonian[i+1][j+1];
                }
            }
        }

        // output
        /*
           for (i=0; i< ct->dim; i++){
           for (j=0; j< ct->dim; j++)
            printf("%lf ", ct->hamiltonian[i][j]*HARTREE_TO_EV);
           printf("\n");
           }
           printf("Etun %lf\n", etun*HARTREE_TO_EV);
           for (i=0; i< nb; i++){
           printf("Vba %lf\n", vba[i]*HARTREE_TO_EV);
           printf("Vdb %lf\n", vdb[i]*HARTREE_TO_EV);
           }

           printf("Gb^-1[eV] =\n");
           for(i=0;i<nb;i++){
           for(j=0;j<nb;j++)
            printf("%lf ", gb[i+j*nb]*HARTREE_TO_EV);
           printf("\n");
           }
         */


        /*
           // inversion of positive definite matrix
           dpotrf_(&itype,&nb,gb,&nb,&info);
           if(info!=0){
           printf("ERROR in calc_tda (factorization)\n");
           exit(-1);
           }
           dpotri_(&itype,&nb,gb,&nb,&info);
           if(info!=0){
           printf("ERROR in calc_tda (inversion)\n");
           exit(-1);
           }
           //  */
        ///*
        //alternative: inversion of general matrix. however, indefinite or positive definite matrix indicates that bridge levels are not correct (compared to D/A)
        dgetrf_(&nb, &nb, gb, &nb, ipiv, &info);
        if ((int) info != 0)
        {
            printf("ERROR in calc_tda (dgetrf) info:%d \n", (int) info);
            exit(-1);
        }
        dgetri_(&nb, gb, &nb, ipiv, work, &lwork, &info);
        if ((int) info != 0)
        {
            printf("ERROR in calc_tda (dgetri) info:%d\n", (int) info);
            exit(-1);
        }
        //*/
        /*
           printf("Gb[eV^-1] =\n");
           for(i=0;i<nb;i++){
           for(j=0;j<nb;j++)
            printf("%lf ", gb[i+j*nb]/HARTREE_TO_EV);
           printf("\n");
           }
         */

        // construct effective two state hamiltonian
        tda       = 0.0;
        Heff[0]   = ct->hamiltonian[0][0];
        Heff[1+2] = ct->hamiltonian[ct->dim-1][ct->dim-1];
        Heff[1]   = ct->hamiltonian[0][ct->dim-1];
        Heff[1+1] = ct->hamiltonian[ct->dim-1][0];
        for (i = 0; i < nb; i++)
        {
            for (j = 0; j < nb; j++)
            {
                Heff[0]   += vdb[i]*gb[i+j*nb]*vdb[j];
                Heff[1+2] += vba[i]*gb[i+j*nb]*vba[j];
                tda       += vdb[i]*gb[i+j*nb]*vba[j];
            }
        }
        Heff[1]   += tda;
        Heff[1+1] += tda;

        // diagonalize effective Hamiltonian
        dsyev(2, Heff, ev_Heff, work, 3*2);

        // get new tunneling energy
        etun = 0.5 * (ev_Heff[0]+ev_Heff[1]);


        printf("niter %d: Direct coupling (HDA)[eV] = %f   Bridge-mediated coupling (TDA)[eV] = %f   Etun[eV] = %f \n",
               niter, ct->hamiltonian[0][ct->dim-1]*HARTREE_TO_EV, tda*HARTREE_TO_EV, etun*HARTREE_TO_EV);

        if (fabs(etun-etun_old) < Etol)
        {
            printf("TDA calculation converged after %d iterations. \n", niter);
            break;
        }
        else if (niter == maxiter-1)
        {
            printf("TDA calculation did not converge after %d iterations. \n", niter);
        }

        etun_old = etun;
    } // end tuning Etun

    return tda;
}


int do_adiabatic(charge_transfer_t *ct, ct_diis_t *diis, FILE *f)
//int do_adiabatic(charge_transfer_t *ct, dftb_broyden_t *broyd, FILE *f)
{
    int    i, j, k, step;
    double old_energy, energy, *ham, almix, fermi_energy = 0, fermi_upper, fermi_lower;

    ham   = ct->hamiltonian_adiab;
    almix = SIMPLE_ALMIX;

    diis->n_prev_vector = 1;
    diis->indx          = 0;

    // calculate the charges from the wave function
    for (i = 0; i < ct->dim; i++)
    {
        ct->q_act[i] = SQR(ct->wf[i]);
    }

    // copy as the first vector

    energy = 1.e10;

    // SCC cycle
    for (step = 0;; step++)
    {

        // too many iterations?
        if (step >= MAXITER_DIIS)
        {
            // copy the solution to the wave function
            for (i = 0; i < ct->dim; i++)
            {
                ct->wf[i] = sqrt(ct->q_act[i]);
            }
            return 1;
        }

        printf("step %5d:", step);
        //fprintf(f, "step %5d:", step);

        // save and print old charges
        printf(" old q =");
        //fprintf(f, " old q =");
        for (i = 0; i < ct->dim; i++)
        {
            //fprintf(f, "%12.9f", ct->q_act[i]);
            printf("%12.9f", ct->q_act[i]);
        }

        //  FOR BROYDEN:
        //  ct->q_old[i] = ct->q_act[i];

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
/*
    // instead of Fermi:
    for (i=0; i<ct->sites; i++)
      diis->prev_q_diff[diis->indx][i] += SQR(ham[i]);
 */
        // FERMI DISTRIBUTION (electronic temperature)
        // calculate the "Fermi energy", store the Fermi coefficients in diis->fermi_coef
        fermi_lower = ct->ev_adiab[0] - 0.1;
        fermi_upper = ct->ev_adiab[ct->dim - 1] + 0.1;
        while (fermi_upper - fermi_lower > FERMI_CONVERG)
        {
            fermi_energy = 0.5 * (fermi_lower + fermi_upper);
            if (fermi_coef_sum(fermi_energy, ct->dim, ct->ev_adiab, diis->fermi_coef, ct->fermi_kt) > 1.0)
            {
                fermi_upper = 0.5 * (fermi_lower + fermi_upper);
            }
            else
            {
                fermi_lower = 0.5 * (fermi_lower + fermi_upper);
            }
        }
        fprintf(f, " fermi E = %9.6f, fermi coefs:", fermi_energy);
        for (i = 0; i < ct->dim; i++)
        {
            fprintf(f, " %8.5f", diis->fermi_coef[i]);
        }
        // mix the eigenvectors, with the appropriate coefficients
        for (i = 0; i < ct->dim; i++)     // i runs over the sites
        {
            for (j = 0; j < ct->dim; j++) // j runs over the individual eigenvectors
            {
                diis->prev_q_diff[diis->indx][i] += diis->fermi_coef[j] * SQR(ham[j * ct->dim + i]);
            }
        }

        // calculate the energy and check convergence
        old_energy = energy;
        energy     = 0.e0;
        /* old - calculate energy by orbitals
           for (i=0; i<ct->sites; i++)
           for (j=0; j<ct->sites; j++)
            energy += ham[i] * ham[j] * ct->hamiltonian[i][j]
         + 0.5 * ct->q_act[i] * ct->q_act[j] * (ct->sic ? (SIC_COEFFICIENT * ct->hubbard[i][j]) : (ct->hubbard[i][j]));
         */
        for (i = 0; i < ct->dim; i++)
        {
            energy += diis->fermi_coef[i] * ct->ev_adiab[i];
        }
        printf(" energy = %12.9f\n", energy);
        //fprintf(f, " energy = %12.9f\n", energy);
        if (fabs(energy - old_energy) < BROYDEN_SCFTOL)
        {
            //fprintf(f, "Converged!\n");
            printf("Converged!\n");
            // copy the solution to the wave function
            for (i = 0; i < ct->dim; i++)
            {
                //ct->wf[i] = sqrt(ct->q_act[i]);
                ct->wf[i] = sqrt(diis->prev_q_input[diis->indx][i] + diis->prev_q_diff[diis->indx][i]);
            }
            //ct->wf[i] = ham[i];
            break;
        }

        // DIIS - Direct Inversion in the Iterative Subspace
        // we need in memory:
        // p        -- matrix (# of bases) x (# of vectors)
        // delta p  -- matrix (# of bases) x (# of vectors)
        // b        -- matrix (# of vectors + 1) x (# of vectors + 1)

        // produce a new vector "ct->q_old", either by simple mixing (in first iterations) or by DIIS (later)
        if (diis->n_prev_vector < DIIS_MAX_PREV_VECTORS)
        {
            // SIMPLE MIXING in a couple of first steps
            for (i = 0; i < ct->dim; i++)
            {
                ct->q_act[i] = diis->prev_q_input[diis->indx][i] + DIIS_INIT_MIXING_PARAM * diis->prev_q_diff[diis->indx][i];
            }
            //ct->q_act[i] = (1.0 - DIIS_INIT_MIXING_PARAM) * ct->q_act[i] + DIIS_INIT_MIXING_PARAM * diis->q_inp_result[i];
            diis->n_prev_vector++;
        }
        else
        {
            // real DIIS here
            // what do the vectors "prev_q_diff" look like?
            //printf("\n DIFF VECTORS\n");
            //for (i=0; i<DIIS_MAX_PREV_VECTORS; i++) {
            //  for (j=0; j<ct->sites; j++)
            //  printf("%14.10f", diis->prev_q_diff[i][j]);
            //printf("\n");
            //}
            // construct the matrix of overlaps
            for (i = 0; i < DIIS_MAX_PREV_VECTORS; i++)
            {
                for (j = 0; j < DIIS_MAX_PREV_VECTORS; j++)
                {
                    diis->aa[i + j * (DIIS_MAX_PREV_VECTORS + 1)] = 0.e0;
                    for (k = 0; k < ct->dim; k++)
                    {
                        diis->aa[i + j * (DIIS_MAX_PREV_VECTORS + 1)] += diis->prev_q_diff[i][k] * diis->prev_q_diff[j][k];
                    }
                    //printf("aa[%d][%d] = %20.15f\n", i, j, diis->aa[i + j * (DIIS_MAX_PREV_VECTORS + 1)]);
                }
            }
            for (i = 0; i < DIIS_MAX_PREV_VECTORS; i++)
            {
                diis->aa[i + DIIS_MAX_PREV_VECTORS * (DIIS_MAX_PREV_VECTORS + 1)] = diis->aa[DIIS_MAX_PREV_VECTORS + i * (DIIS_MAX_PREV_VECTORS + 1)] = -1.e0;
            }
            diis->aa[DIIS_MAX_PREV_VECTORS + DIIS_MAX_PREV_VECTORS * (DIIS_MAX_PREV_VECTORS + 1)] = 0.e0;
            //printf("\n MATRIX AA:\n");
            //for (i=0; i<=DIIS_MAX_PREV_VECTORS; i++) {
            //  for (j=0; j<=DIIS_MAX_PREV_VECTORS; j++)
            //  printf("%14.10f", diis->aa[i + j * (DIIS_MAX_PREV_VECTORS + 1)]);
            //printf("\n");
            //}
            // construct the vector of right hand sides
            for (i = 0; i < DIIS_MAX_PREV_VECTORS; i++)
            {
                diis->bb[i] = 0.e0;
            }
            diis->bb[DIIS_MAX_PREV_VECTORS] = -1.e0;
            //printf(" MATRIX BB:\n");
            //for (i=0; i<=DIIS_MAX_PREV_VECTORS; i++)
            //  printf("%14.10f", diis->bb[i]);
            //printf("\n");
            // solve it with lapack
            //printf(" DGESV: %ld", dgesv(DIIS_MAX_PREV_VECTORS + 1, diis->aa, diis->bb, diis->ipiv));
            dgesv(DIIS_MAX_PREV_VECTORS + 1, diis->aa, diis->bb, diis->ipiv);
            // the column of solution is now in diis->bb -- print it!
            //printf(" DIIS coeffs");
            //for (i=0; i<=DIIS_MAX_PREV_VECTORS; i++)
            //  printf("%9.5f", diis->bb[i]);
            for (j = 0; j < ct->dim; j++)
            {
                ct->q_act[j] = 0.e0;
                for (i = 0; i < DIIS_MAX_PREV_VECTORS; i++)
                {
                    ct->q_act[j] += diis->bb[i] * diis->prev_q_input[i][j];
                }
            }
        }

        // move the counter in a loop
        diis->indx = (diis->indx + 1) % DIIS_MAX_PREV_VECTORS;

/*
    // FERMI DISTRIBUTION + BROYDEN MIXING
    // mix all eigenvectors with appropriate coefficients

    // calculate the "Fermi energy", store the Fermi coefficients in broyd->df
    fermi_lower = ct->ev_adiab[0] - 0.1;
    fermi_upper = ct->ev_adiab[ct->sites - 1] + 0.1;
    while (fermi_upper - fermi_lower > FERMI_CONVERG) {
      fermi_energy = 0.5 * (fermi_lower + fermi_upper);
      if (fermi_coef_sum(fermi_energy, ct->sites, ct->ev_adiab, broyd->df, ct->fermi_kt) > 1.0)
        fermi_upper = 0.5 * (fermi_lower + fermi_upper);
      else
        fermi_lower = 0.5 * (fermi_lower + fermi_upper);
    }
    //fermi_coef2 = 1 / (exp((ct->ev_adiab[1] - ct->ev_adiab[0]) / FERMI_KT) + 1);
    //printf(" fermi coef2 = %7.4f", fermi_coef2);
    fprintf(f, " fermi E = %9.6f, fermi coefs:", fermi_energy);
    for (i=0; i<ct->sites; i++)
      fprintf(f, " %8.5f", broyd->df[i]);

    // mix the eigenvectors, with the appropriate coefficients
    for (i=0; i<ct->sites; i++) { // i runs over the sites
      ct->q_act[i] = 0.e0;
      for (j=0; j<ct->sites; j++) // j runs over the individual eigenvectors
        ct->q_act[i] += broyd->df[j] * SQR(ham[j * ct->sites + i]);
    }

    // mix it - Broyden
    broyden(step, BROYDEN_ALMIX, ct->sites, ct->q_old, ct->q_act, *broyd);
    for (i=0; i<ct->sites; i++)
      ct->q_act[i] = ct->q_old[i];
 */

/*
    // SIMPLE MIXING
    //
    // calculate the charges
    for (i=0; i<ct->sites; i++)
      ct->q_act[i] = SQR(ham[i]);
    // calculate oldness - the content of old vector in the new one
    oldness = 0.e0;
    for (i=0; i<ct->sites; i++)
      oldness += ct->q_act[i] * ct->q_old[i];
    // the orthogonal vector is now: new_ortho = new - oldness * old
    // calculate the new vector as: old + almix(t) * new_ortho
    for (i=0; i<ct->sites; i++)
      ct->q_act[i] = (1.0 - almix) * ct->q_old[i] + almix * (ct->q_act[i] - oldness * ct->q_old[i]);
    // attenuate the mixing coefficient
    almix *= ALMIX_ATTENUATOR;
 */

    }

    return 0;
}


int do_surface_hopping(charge_transfer_t *ct, ct_diis_t *diis, FILE *f)
{
    /* this should be done by the modification of do_adiabatic
     * followed by the testing of scalar product of the new wavefunction - ground and 1st excited state - with the old one
     */
    /* ct->surface = 0 (ground state) or 1 (1st excited state) */
    int    i, j, k, step;
    double old_energy, energy, *ham, s01, s10, nonadiab_coupling;

    ham                 = ct->hamiltonian_adiab;
    diis->n_prev_vector = 1;
    diis->indx          = 0;

    // calculate the charges from the wave function
    for (i = 0; i < ct->sites; i++)
    {
        ct->q_act[i] = SQR(ct->wf[i]);
    }

    // SCC cycle
    for (step = 0, energy = 1.e10;; step++)
    {

        if (step >= MAXITER_DIIS)
        {
            /* DIIS NOT CONVERGED! */
            //printf("DIIS not converged !!!      ");
            //fprintf(f, "DIIS not converged !!!      ");
            fprintf(f, "DIIS   0 ");
            fprintf(f, "             ");
            for (i = 0; i < ct->sites; i++)
            {
                fprintf(f, "        ");
            }
            fprintf(f, "                           DIIS not converged so    hopping not tested... ");
            fprintf(f, " %d", ct->surface);
            for (i = 0; i < ct->sites; i++)
            {
                fprintf(f, "%8.4f", SQR(ct->wf[i]));
            }
            fprintf(f, "\n");

            return 1;
        }

        //printf("step %5d:", step);

        //printf(" old q =");
        //for (i=0; i<ct->sites; i++)
        //printf("%9.5f", ct->q_act[i]);

        // construct the self-consistent Hamiltonian
        for (i = 0; i < ct->sites; i++)
        {
            ham[i + i * ct->sites] = ct->hamiltonian[i][i];
            for (j = 0; j < ct->sites; j++)
            {
                ham[i + i * ct->sites] += ct->hubbard[i][j] * ct->q_act[j];
                if (i != j)
                {
                    ham[i + j * ct->sites] = ct->hamiltonian[i][j];
                }
            }
        }

        /*
           // orbital shift:
           // 1. preliminary diagonalization
           dsyev(ct->sites, ham, ct->ev_adiab, ct->work_adiab, 3*ct->sites);
           // 2. correction to the elements of the Hamiltonian
           //    Delta H = C * delta-Epsilon * C^T
           for (i=0; i<ct->sites; i++)
           ct->work_adiab[i] = ham[i];
           for (i=0; i<ct->sites; i++) {
           ham[i + i * ct->sites] = ct->hamiltonian[i][i] - SFH_ORBITAL_SHIFT * ct->work_adiab[i] * ct->work_adiab[i];
           for (j=0; j<ct->sites; j++) {
            ham[i + i * ct->sites] += ct->hubbard[i][j] * ct->q_act[j];
            if (i!=j) ham[i + j * ct->sites] = ct->hamiltonian[i][j] - SFH_ORBITAL_SHIFT * ct->work_adiab[i] * ct->work_adiab[j];
           }
           }
           // done with orbital shift
         */

        // diagonalize this Hamiltonian
        dsyev(ct->sites, ham, ct->ev_adiab, ct->work_adiab, 3*ct->sites);

        // copy previous charge and the difference to the arrays
        for (i = 0; i < ct->sites; i++)
        {
            diis->prev_q_input[diis->indx][i] = ct->q_act[i];
            diis->prev_q_diff[diis->indx][i]  = SQR(ham[i]) - ct->q_act[i];
        }

        old_energy = energy;
        energy     = ct->ev_adiab[0];
        //printf(" energy = %12.9f\n", energy);
        if (fabs(energy - old_energy) < BROYDEN_SCFTOL)
        {
            //printf("DIIS Converged in %3d steps ", step);
            //fprintf(f, "DIIS Converged in %3d steps ", step);
            fprintf(f, "DIIS %3d ", step);
            break;
        }

        if (diis->n_prev_vector < DIIS_MAX_PREV_VECTORS)
        {
            // SIMPLE MIXING in a couple of first steps
            for (i = 0; i < ct->sites; i++)
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
                    for (k = 0; k < ct->sites; k++)
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

            for (j = 0; j < ct->sites; j++)
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
    //printf(" Eigenvalues:");
    fprintf(f, " Eigenvalues:");
    for (i = 0; i < ct->sites; i++)
    {
        //fprintf(f, "%10.7f", ct->ev_adiab[i]);
        fprintf(f, "%8.5f", ct->ev_adiab[i]);
        //printf("%8.5f", ct->ev_adiab[i] - ct->ev_adiab[0]);
    }
    //printf("    ");
    fprintf(f, "    ");

    // check the energy gap between the ground and the first excited state
    if (ct->ev_adiab[1] - ct->ev_adiab[0] < SFH_LIMIT)
    {
        /* attempt surface hopping
         * a variation of the AS3 method in Fabiano et al., Chem. Phys. 351 (2008) 111-116
         */
        s01 = 0.e0;
        for (i = 0; i < ct->sites; i++)
        {
            s01 += (ct->surface == 0 ? ct->wf[i] : ct->wf_exc[i]) * ham[ct->sites + i];
        }
        s10 = 0.e0;
        for (i = 0; i < ct->sites; i++)
        {
            s10 += (ct->surface == 0 ? ct->wf_exc[i] : ct->wf[i]) * ham[i];
        }
        /* calculate the (quasi-)non-adiabatic coupling */
        nonadiab_coupling = 1 - pow(1 - (fabs(s01) + fabs(s10))*0.5, SFH_EXPONENT);
        /* print them out */
        //printf("wf.ground = %6.3f, wf.excited = %6.3f ", overlap_ground, overlap_excited);
        fprintf(f, "s01 = %7.4f, s10 = %7.4f, coupling = %6.3f ", s01, s10, nonadiab_coupling);

        /* check if hopping should be performed */
        if (nonadiab_coupling > rand()/(double)RAND_MAX)
        {
            fprintf(f, "hopping between states");
            if (ct->surface == 0) // hopping 0->1
            {
                for (i = 0; i < ct->sites; i++)
                {
                    ct->wf[i]     = ham[ct->sites + i];
                    ct->wf_exc[i] = ham[i];
                }
                ct->surface = 1;
            }
            else              // hopping 1->0
            {
                for (i = 0; i < ct->sites; i++)
                {
                    ct->wf[i]     = ham[i];
                    ct->wf_exc[i] = ham[ct->sites + i];
                }
                ct->surface = 0;
            }
        }
        else
        {
            fprintf(f, "no hopping now...     ");
            for (i = 0; i < ct->sites; i++)
            {
                if (ct->surface == 0) // system is remaining in state 0
                {
                    ct->wf[i]     = ham[i];
                    ct->wf_exc[i] = ham[ct->sites + i];
                }
                else            // system is remaining in state 1
                {
                    ct->wf[i]     = ham[ct->sites + i];
                    ct->wf_exc[i] = ham[i];
                }
            }
        }
    }
    else
    {
        //printf("Egap > threshold, no surf hopping ");
        fprintf(f, "                                                hopping not tested... ");
        for (i = 0; i < ct->sites; i++)
        {
            if (ct->surface == 0) // system is remaining in state 0
            {
                ct->wf[i]     = ham[i];
                ct->wf_exc[i] = ham[ct->sites + i];
            }
            else              // system is remaining in state 1
            {
                ct->wf[i]     = ham[ct->sites + i];
                ct->wf_exc[i] = ham[i];
            }
        }
    }
    fprintf(f, " %d", ct->surface);
    for (i = 0; i < ct->sites; i++)
    {
        //printf("%8.4f", SQR(ct->wf[i]));
        fprintf(f, "%8.4f", SQR(ct->wf[i]));
    }
    //printf("\n");
    fprintf(f, "\n");

    return 0;
}


double do_born_oppenheimer(charge_transfer_t *ct, dftb_broyden_t *broyd, double *fermi_coeff, FILE *f)
{
    int    i, j, k, step, ier;
    double old_energy, energy, energy1, energy2, *ham, fermi_energy, fermi_upper, fermi_lower, q_sum, s01, best_overlap;

    ham = ct->hamiltonian_adiab;

    // calculate the charges from the wave function
    for (i = 0; i < ct->dim; i++)
    {
        ct->q_act[i] = ct->q_old[i] = SQR(ct->wf[i]);
    }

    // copy as the first vector

    energy = 1.e10;
    // SCC cycle
    for (step = 0;; step++)
    {

        // too many iterations? break and indicate failure
        if (step >= MAXITER_BROYDEN)
        {
            fprintf(f, "Broyden/Fermi failed to converge!\n");
            return -1;
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
        ier = dsyev(ct->dim, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);
        if ((int) ier)
        {
            printf("\nDSYEV: ier = %d\nEXITING!\n\n", (int) ier);
            exit(-1);
        }


        // FERMI DISTRIBUTION + BROYDEN MIXING
        // mix all eigenvectors with appropriate coefficients

        // calculate the "Fermi energy", store the Fermi coefficients in fermi_coeff
        fermi_lower = ct->ev_adiab[0] - 0.1;
        fermi_upper = ct->ev_adiab[ct->dim - 1] + 0.1;
        while (fermi_upper - fermi_lower > FERMI_CONVERG)
        {
            fermi_energy = 0.5 * (fermi_lower + fermi_upper);
            if (fermi_coef_sum(fermi_energy, ct->dim, ct->ev_adiab, fermi_coeff, ct->fermi_kt) > 1.0)
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
        energy1    = energy2 = 0.0;
        // indices i and j run over the sites
        for (i = 0; i < ct->dim; i++)
        {
            for (j = 0; j < ct->dim; j++)
            {
                // TB Hamiltonian
                // index k runs over the adiabatic states (eigenvectors)
                for (k = 0; k < ct->dim; k++)
                {
                    energy1 += fermi_coeff[k] * ham[k * ct->dim + i] * ham[k * ct->dim + j] * ct->hamiltonian[i][j];
                }
                // Hubbard / gamma terms
                energy2 += 0.5 * ct->hubbard[i][j] * ct->q_act[i] * ct->q_act[j];
            }
        }
        energy = energy1 + energy2;

        if (step > 1 && fabs(energy - old_energy) < BROYDEN_SCFTOL)
        {
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

    fprintf(f, " %2d ", step);
    //fprintf(f, " %2d iters, SCC-Q =", step);
    for (i = 0; i < ct->dim; i++)
    {
        fprintf(f, "%9.6f", ct->q_act[i]);
    }


    // copy the lowest pure adiabatic state to the wave function
    //for (i=0; i<ct->dim; i++)
    //  ct->wf[i] = ham[i];



    // calculate and print the energies of pure adiabatic states (indexed by k)
    //fprintf(f, " E(adiab) =");
    for (k = 0; k < ct->dim; k++)
    {
        energy = 0.0;
        for (i = 0; i < ct->dim; i++)
        {
            for (j = 0; j < ct->dim; j++)
            {
                // TB Hamiltonian
                energy += ham[k * ct->dim + i] * ham[k * ct->dim + j] * ct->hamiltonian[i][j]
                    // Hubbard / gamma terms
                    +  0.5 * ct->hubbard[i][j] * SQR(ham[k * ct->dim + i]) * SQR(ham[k * ct->dim + j]);
            }
        }
        fprintf(f, "%10.7f", energy);
    }




///*
// calc overlap with adiabatic states i //
    for (i = 0; i < ct->dim; i++)
    {
        ct->born_overlap[i] = 0.0;
        for (j = 0; j < ct->dim; j++)
        {
            ct->born_overlap[i] +=   ham[i*ct->dim + j]*ct->wf[j];
        }
    }

    // wf will be adiabatic state that is most similar //
    best_overlap = 0.0;
    for (i = 0; i < ct->dim; i++)
    {
        if (fabs(ct->born_overlap[i]) > fabs(best_overlap))
        {
            best_overlap = ct->born_overlap[i];
            for (j = 0; j < ct->dim; j++)
            {
                ct->wf[j] = ham[i*ct->dim + j];
            }

            printf("best overlap with adiabatic state no. %d. overlap %lf \n", i, ct->born_overlap[i]);
        }
    }
//*/

/*
   for(i=0;i<ct->dim;i++){ //calc overlap with state i
    s01 = 0.e0;
    for (j=0; j<ct->dim; j++)
      s01 += ham[i*ct->dim + j]*ct->wf[j];
   printf("S01 %f\n", s01);
    if(fabs(s01)>0.8){ //threshold
        for (j=0; j<ct->dim; j++)
          ct->wf[j] = ham[i*ct->sites + j];
    break;
    }
   }//end loop over eigenvalues
   if(i==ct->dim){ printf("Eigenfunctions changed too much! Born-Oppenheimer not reasonable!"); exit(-1); }
 */
    return SQR(best_overlap);
}


int do_adiab_fermi(charge_transfer_t *ct, dftb_broyden_t *broyd, double *fermi_coeff, FILE *f)
{
    int    i, j, k, step;
    double old_energy, energy, energy1, energy2, *ham, fermi_energy, fermi_upper, fermi_lower, q_sum;

    ham = ct->hamiltonian_adiab;

    // calculate the charges from the wave function
    for (i = 0; i < ct->dim; i++)
    {
        ct->q_act[i] = ct->q_old[i] = SQR(ct->wf[i]);
    }

    // copy as the first vector

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
        // mix all eigenvectors with appropriate coefficients

        // calculate the "Fermi energy", store the Fermi coefficients in fermi_coeff
        fermi_lower = ct->ev_adiab[0] - 0.1;
        fermi_upper = ct->ev_adiab[ct->dim - 1] + 0.1;
        while (fermi_upper - fermi_lower > FERMI_CONVERG)
        {
            fermi_energy = 0.5 * (fermi_lower + fermi_upper);
            if (fermi_coef_sum(fermi_energy, ct->dim, ct->ev_adiab, fermi_coeff, ct->fermi_kt) > 1.0)
            {
                fermi_upper = 0.5 * (fermi_lower + fermi_upper);
            }
            else
            {
                fermi_lower = 0.5 * (fermi_lower + fermi_upper);
            }
        }
        fprintf(f, " fermi E = %9.6f, coefs", fermi_energy);
        for (i = 0; i < ct->dim; i++)
        {
            fprintf(f, "%8.5f", fermi_coeff[i]);
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

        // print new charges
        fprintf(f, " mixedQ =");
        for (i = 0; i < ct->dim; i++)
        {
            fprintf(f, "%9.6f", ct->q_act[i]);
        }

        // calculate the energy and check convergence
        old_energy = energy;
        energy1    = energy2 = 0.0;
        // indices i and j run over the sites
        for (i = 0; i < ct->dim; i++)
        {
            for (j = 0; j < ct->dim; j++)
            {
                // TB Hamiltonian
                // index k runs over the adiabatic states (eigenvectors)
                for (k = 0; k < ct->dim; k++)
                {
                    energy1 += fermi_coeff[k] * ham[k * ct->dim + i] * ham[k * ct->dim + j] * ct->hamiltonian[i][j];
                }
                // Hubbard / gamma terms
                energy2 += 0.5 * ct->hubbard[i][j] * ct->q_act[i] * ct->q_act[j];
            }
        }
        energy = energy1 + energy2;

        fprintf(f, " E1=%10.7f E2=%10.7f totE=%13.10f", energy1, energy2, energy);

        if (step > 1 && fabs(energy - old_energy) < BROYDEN_SCFTOL)
        {
            fprintf(f, " Converged!\n");
            // copy the solution to the wave function
            for (i = 0; i < ct->dim; i++)
            {
                ct->wf[i] = sqrt(ct->q_act[i]);
            }
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
        // print mixed charges
        fprintf(f, " newQ =");
        for (i = 0; i < ct->dim; i++)
        {
            fprintf(f, "%9.6f", ct->q_act[i]);
        }
        fprintf(f, "\n");

    } /* end SCC cycle */

    return 0;
}

int do_adiab_fermi_onestate(charge_transfer_t *ct, dftb_broyden_t *broyd, double *fermi_coeff, FILE *f)
{
    int    i, j, k, step;
    double old_energy, energy, energy1, energy2, *ham, fermi_energy, fermi_upper, fermi_lower, q_sum;

    ham = ct->hamiltonian_adiab;

    // calculate the charges from the wave function
    for (i = 0; i < ct->dim; i++)
    {
        ct->q_act[i] = ct->q_old[i] = SQR(ct->wf[i]);
    }

    // copy as the first vector

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
        // mix all eigenvectors with appropriate coefficients

        // calculate the "Fermi energy", store the Fermi coefficients in fermi_coeff
        fermi_lower = ct->ev_adiab[0] - 0.1;
        fermi_upper = ct->ev_adiab[ct->dim - 1] + 0.1;
        while (fermi_upper - fermi_lower > FERMI_CONVERG)
        {
            fermi_energy = 0.5 * (fermi_lower + fermi_upper);
            if (fermi_coef_sum(fermi_energy, ct->dim, ct->ev_adiab, fermi_coeff, ct->fermi_kt) > 1.0)
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
        energy1    = energy2 = 0.0;
        // indices i and j run over the sites
        for (i = 0; i < ct->dim; i++)
        {
            for (j = 0; j < ct->dim; j++)
            {
                // TB Hamiltonian
                // index k runs over the adiabatic states (eigenvectors)
                for (k = 0; k < ct->dim; k++)
                {
                    energy1 += fermi_coeff[k] * ham[k * ct->dim + i] * ham[k * ct->dim + j] * ct->hamiltonian[i][j];
                }
                // Hubbard / gamma terms
                energy2 += 0.5 * ct->hubbard[i][j] * ct->q_act[i] * ct->q_act[j];
            }
        }
        energy = energy1 + energy2;

        if (step > 1 && fabs(energy - old_energy) < BROYDEN_SCFTOL)
        {
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

    fprintf(f, " %2d ", step);
    //fprintf(f, " %2d iters, SCC-Q =", step);
    for (i = 0; i < ct->dim; i++)
    {
        fprintf(f, "%9.6f", ct->q_act[i]);
    }

    // copy the lowest pure adiabatic state to the wave function
    for (i = 0; i < ct->dim; i++)
    {
        ct->wf[i] = ham[i];
    }

    // calculate and print the energies of pure adiabatic states (indexed by k)
    //fprintf(f, " E(adiab) =");
    for (k = 0; k < ct->dim; k++)
    {
        energy = 0.0;
        for (i = 0; i < ct->dim; i++)
        {
            for (j = 0; j < ct->dim; j++)
            {
                // TB Hamiltonian
                energy += ham[k * ct->dim + i] * ham[k * ct->dim + j] * ct->hamiltonian[i][j]
                    // Hubbard / gamma terms
                    +  0.5 * ct->hubbard[i][j] * SQR(ham[k * ct->dim + i]) * SQR(ham[k * ct->dim + j]);
            }
        }
        fprintf(f, "%10.7f", energy);
    }
    return 0;
}


// Trajectory surface hopping method using original diabatic Landau-Zener formula W.X Feb.2019

int do_lzsh_diab(int step, charge_transfer_t *ct, dftb_t *dftb, t_state *state_global, t_mdatoms *mdatoms, rvec *x_ct, FILE *f)
{

    int    i, j, k, l, m, iatom, hopsurf, prev_surf, indextmp, counter;
    double *ham, *ekin, *dist, de[1];
    double probmax, random_number, deriv, coupling;
    double ekintot, ekintmp, kai, detmp, disttmp, norm;
    dvec   com, coord, masscoord;


    snew(ekin, ct->dim);
    snew(dist, ct->dim);
//    snew(index, ct->dim);

    ham = ct->tfl_mean_ham;

/*    counter = 0;
    for (i = 0; i < ct->pool_size; i++)
    {
        clear_dvec(com);
        for (j = 0; j < ct->pool_site[i].atoms; j++)
        {
            for (k = 0; k < 3; k++)
            {
                coord[k] =  x_ct[ct->pool_site[i].atom[j]][k] * NM_TO_BOHR;//take here x_ct for wich molecules were made whole again. in gromacs_x all atoms are put in the box.
            }
            dsvmul(mdatoms->massT[ct->pool_site[i].atom[j]], coord, masscoord);
            dvec_inc(com, masscoord);
        }
        dsvmul(ct->adapt_inv_tot_mass, com, ct->pool_site[i].com);
        printf("residue %d  COM[Angstrom]: %f %f %f\n", ct->pool_site[i].resnr, ct->pool_site[i].com[XX] * 0.52, ct->pool_site[i].com[YY] * 0.52, ct->pool_site[i].com[ZZ] * 0.52);
    }

   exit(0);
*/
    for (i = 0; i < SQR(ct->dim); i++)
        ham[i] = *(ct->hamiltonian[0]+i);

// diagonalize the Hamiltonian
    dsyev(ct->dim, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);

// record eigenvectors
    for (i = 0; i < ct->dim; i++)
        for (j = 0; j < ct->dim; j++)
            ct->tfs_vector[i][j] = ham[i * ct->dim + j];

// set up initial surface
    if (ct->tfs_initialization_step)
    {
// determine inital surface using eigenvector as a criterion
       for (i = 0; i < 2*ct->dim; i++)
       {
           ct->tfs_popul[i] = 0.0;
           ct->tfs_diab[i] = 0.0;
       }

       for (i = 0; i < ct->dim; i++)
       {
           ct->tfs_diab[i] = ct->wf[i];
           if (ct->wf[i] == 1.0)
              ct->surface = i;
       }
       ct->tfs_initialization_step = 0;
    }

    if (step > 0)
    {
       probmax = 0.0;
       for (i = 0; i < ct->dim; i++)
       {
// calculate the diabatic energy gaps at two time steps before and previous time step
           de[0] = ct->tfl_mean_ham_full[ct->surface][ct->surface] - ct->tfl_mean_ham_full[i][i];
           de[1] = ct->hamiltonian[ct->surface][ct->surface] - ct->hamiltonian[i][i];
           coupling = (fabs(ct->tfl_mean_ham_full[ct->surface][i]) + fabs(ct->hamiltonian[ct->surface][i]))/2.0;
// examine local energy minimum
           if (de[0]*de[1] < 0.0)
           {
// derivative of diabatic energy gap with respect to time
              deriv = (de[1] - de[0])/ct->rk_timestep;
// Landau-Zener formula (Phy.Rev.A 84, 014701, 2011)
              ct->surf_prob[i] = 1.0 - exp(-2.0*M_PI*SQR(coupling)/fabs(deriv));

              printf("%f \n", ct->surf_prob[i]);
              if (ct->surf_prob[i] > probmax)
              {
                 probmax = ct->surf_prob[i];
                 hopsurf = i;
              }
           }
       }

// compared with a random number to determine the hopping event
       prev_surf = ct->surface;
       random_number = drand48();
       if (probmax > random_number)
       {
// kinetic energy of QM zone
           for (k = 0; k < ct->sites; k++)
               ekin[k] = 0.0;

           ekintot = 0.0;
           for (iatom = 0; iatom < mdatoms->homenr; iatom++)
           {
               for (k = 0; k < ct->sites; k++)
               {
                   for (l = 0; l < ct->site[k].atoms; l++)
                   {
                       if (iatom == ct->site[k].atom[l])
                       {
                          for (m = 0; m < DIM; m++)
                          {
                              ekintot += 0.5*dftb->phase1[k].mass[l]*AMU_TO_AU*SQR(state_global->v[iatom][m]*NM_TO_BOHR/PS_TO_AU);
                              ekin[k] += 0.5*dftb->phase1[k].mass[l]*AMU_TO_AU*SQR(state_global->v[iatom][m]*NM_TO_BOHR/PS_TO_AU);
                          }
                       }
                   }
               }
           }
           printf("Kinetic energy is %f eV\n", ekintot*27.2114);

// frustrated hops
           ekintmp = ekin[hopsurf] + ekin[ct->surface];
           detmp = ct->hamiltonian[hopsurf][hopsurf] - ct->hamiltonian[ct->surface][ct->surface];

           printf("Two-site kinetic energy is %f eV \n", ekintmp*27.2114);
           printf("Potential energy difference is %f eV \n", detmp*27.2114);

           if (1.0 - detmp/ekintmp < 0.0)
           {
              printf("Frustrated hop! \n");
              fprintf(f,"%8.2f %10d %17d %21f %19f Frustrated\n",(step-1)*ct->rk_timestep*1000./PS_TO_AU, ct->surface, hopsurf, probmax, random_number);
           }
           else
           {
              kai = sqrt(1.0 - detmp/ekintmp);

              for (iatom = 0; iatom < mdatoms->homenr; iatom++)
              {
                  for (k = 0; k < ct->sites; k++)
                  {
                      if (k == hopsurf || k == ct->surface)
                      {
                          for (l = 0; l < ct->site[k].atoms; l++)
                          {
                              if (iatom == ct->site[k].atom[l])
                              {
                                 for (m = 0; m < DIM; m++)
                                     state_global->v[iatom][m] *= kai;
                              }
                          }
                      }
                  }
              }

              fprintf(f,"%8.2f %10d %17d %21f %19f \n",(step-1)*ct->rk_timestep*1000./PS_TO_AU, ct->surface, hopsurf, probmax, random_number);

              ct->tfs_diab[ct->surface] = 0.0;
              ct->wf[ct->surface]       = 0.0;

              ct->surface = hopsurf;

              ct->tfs_diab[ct->surface] = 1.0;
              ct->wf[ct->surface]       = 1.0;
           }
       }
    }

// adaptive sites
/*    if (ct->surface != prev_surf)
    {
    counter = 0;
    for (i = 0; i < ct->pool_size; i++)
    {
        clear_dvec(com);
        for (j = 0; j < ct->pool_site[i].atoms; j++)
        {
            for (k = 0; k < 3; k++)
            {
                coord[k] =  x_ct[ct->pool_site[i].atom[j]][k] * NM_TO_BOHR;//take here x_ct for wich molecules were made whole again. in gromacs_x all atoms are put in the box.
            }
            dsvmul(mdatoms->massT[ct->pool_site[i].atom[j]], coord, masscoord);
            dvec_inc(com, masscoord);
        }
        dsvmul(ct->adapt_inv_tot_mass, com, ct->pool_site[i].com);
        printf("residue %d  COM[Angstrom]: %f %f %f\n", ct->pool_site[i].resnr, ct->pool_site[i].com[XX] * 0.52, ct->pool_site[i].com[YY] * 0.52, ct->pool_site[i].com[ZZ] * 0.52);
    }

    for (i = 0; i < ct->dim; i++)
    {
        for (j = i + 1; j < ct->dim; j++)
        {
            if (dist[i] > dist[j])
            {
               disttmp = dist[i];
               dist[i] = dist[j];
               dist[j] = disttmp;

               indextmp = index[i];
               index[i] = index[j];
               index[j] = indextmp;
            }
        }
    }

    }
*/
// store the old hamiltonian
    for (i = 0; i < ct->dim; i++)
        for (j = 0; j < ct->dim; j++)
            ct->tfl_mean_ham_full[i][j] = ct->hamiltonian[i][j];

// calculate adibatic expansion coeffients
    for (i = 0; i < 2*ct->dim; i++)
        ct->tfs_popul[i] = 0.0;

    for (i = 0; i < ct->dim; i++)
        ct->tfs_popul[i] = ct->tfs_vector[i][ct->surface];

    printf("Current Surface: %d\n", ct->surface);

    return 0;
}


// Trajectory surface hopping method using refomulated Landau-Zener formula W.X Nov.2018
// Ref: Phys.Rev A 84, 014701 (2011)

int do_lzsh_adiab(int step, charge_transfer_t *ct, FILE *f)
{

    int    i, j, k, l, hopsurf, *ind;
    double *ham, *initprob, de[2], de_diab[1];
    double probsum, probmax, random_number, deriv2, coupling, deriv;
    double ekintot, ekintmp, detmp;
    double nc, aa, bb, cc, dd, ee;


    ham = ct->tfl_mean_ham;

    snew(initprob, ct->dim);
    snew(ind, ct->dim);

    for (i = 0; i < SQR(ct->dim); i++)
        ham[i] = *(ct->hamiltonian[0]+i);

// diagonalize the Hamiltonian
    dsyev(ct->dim, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);

// record eigenvectors
    for (i = 0; i < ct->dim; i++)
        for (j = 0; j < ct->dim; j++)
            ct->tfs_vector[i][j] = ham[i * ct->dim + j];

// set up initial surface
    if (ct->tfs_initialization_step)
    {
// determine inital surface using eigenvector as a criterion
       probsum = 0.;
       random_number = drand48();
       for (i = 0; i < ct->dim; i++)
       {
           initprob[i] = 0.;
           for (j = 0; j < ct->dim; j++)
               initprob[i] += ct->tfs_vector[i][j]*ct->wf[j];
           initprob[i] = SQR(initprob[i]);

           probsum += initprob[i];
           if (random_number < probsum)
           {
               ct->surface = i;
               printf("\n Initial surface: %10d %12.5f %12.5f \n", ct->surface, probsum, random_number);
               break;
           }
       }

       for (i = 0; i < 2*ct->dim; i++)
       {
           ct->tfs_diab[i] = 0.0;
           ct->wf[i]       = 0.0;
       }

       for (i = 0; i < ct->dim; i++)
       {
           ct->tfs_diab[i] = ct->tfs_vector[ct->surface][i];
           ct->wf[i]       = ct->tfs_vector[ct->surface][i];
       }

       for (i = 0; i < 2*ct->dim; i++)
           ct->tfs_popul[i] = 0.0;

       ct->tfs_popul[ct->surface] = 1.0;

       ct->tfs_initialization_step = 0;
    }

// calculate the adiabatic energy gaps at current time
    if (step > 1)
    {
       probmax = 0.;
       for (i = 0; i < ct->dim; i++)
       {
// calculate the adiabatic energy gaps at two time steps before and previous time step
           de[0] = fabs(ct->ev_adiab_old2[ct->surface] - ct->ev_adiab_old2[i]);
           de[1] = fabs(ct->ev_adiab_old[ct->surface] - ct->ev_adiab_old[i]);
           de[2] = fabs(ct->ev_adiab[ct->surface] - ct->ev_adiab[i]);

           de_diab[0] = ct->ev_adiab[i] - ct->ev_adiab_old2[ct->surface];
           de_diab[1] = ct->ev_adiab[ct->surface] - ct->ev_adiab_old2[i];

           coupling = de[1]/2.0;
// calculate hopping probability
           ct->surf_prob[i] = 0.;
// examine local energy minimum
           if (de[1] < de[0] && de[1] < de[2])
           {
// second derivative of adiabatic energy gap with respect to time
              deriv2 = (de[0] + de[2] - 2.0*de[1])/SQR(ct->rk_timestep);
// Landau-Zener formula (Phy.Rev.A 84, 014701, 2011)
              ct->surf_prob[i] = exp(-M_PI/2.0*sqrt(CUB(de[1])/deriv2));
// derivative of diabatic energy gap with respect to time
//              deriv = fabs(de_diab[1] - de_diab[0])/ct->rk_timestep;
// Landau-Zener formula (Phy.Rev.A 84, 014701, 2011)
//              ct->surf_prob[i] = exp(-2.0*M_PI*SQR(coupling)/fabs(deriv));

              if (ct->surf_prob[i] > probmax)// && fabs(ct->surface - i) < 3)
              {
                 probmax = ct->surf_prob[i];
                 hopsurf = i;
              }
           }
       }

// compared with a random number to determine the hopping event
       random_number = drand48();
       if (probmax > random_number)
       {
          fprintf(f,"%8.2f %10d %17d %21f %19f \n",(step-1)*ct->rk_timestep*1000./PS_TO_AU, ct->surface, hopsurf, probmax, random_number);
          ct->surface = hopsurf;
       }
    }

// calculate diabatic expansion coeffients
    for (i = 0; i < ct->dim; i++)
    {
        ct->tfs_diab[i] = ct->tfs_vector[ct->surface][i];
        ct->wf[i]       = ct->tfs_vector[ct->surface][i];
    }

    for (i = 0; i < ct->dim; i++)
    {
        ct->ev_adiab_old2[i] = ct->ev_adiab_old[i];
        ct->ev_adiab_old[i]  = ct->ev_adiab[i];
    }

// calculate adibatic expansion coeffients
    for (i = 0; i < 2*ct->dim; i++)
        ct->tfs_popul[i] = 0.0;

    ct->tfs_popul[ct->surface] = 1.0;

    printf("Current Surface: %d\n", ct->surface);

    return 0;
}

int do_fssh_diab(int step, charge_transfer_t *ct, dftb_t *dftb, t_state *state_global, t_mdatoms *mdatoms, FILE *f)
{
// Declaration:
//
// surface hopping betwteen different sites
//
//                                          -----------Weiwei Xie on Jan. 2019


    int    i, j, k, l, m, n, prev_surf, hopsurf, iatom, istep, nstep, hop, lamtmp;
    double *ham, *old_ham, **ham_full, *new_ham, *diapop, *diapop_old, *ekin, *ddiapop;
    double tdse_timestep = 1.e-5*PS_TO_AU;
    double cur_weight, rk_timestep0, probsum, random_number;
    double decay_time, curr_surf_fac, ddiapopsum;
    double dot_product, ekintot, ekintmp, de, kai;
    double nc, kai1, kai2, delta, aa, bb, cc, dd, ee;


    old_ham  = ct->tfl_old_ham_adiab;
    new_ham  = ct->hamiltonian_adiab;
    ham      = ct->tfl_mean_ham;
    ham_full = ct->tfl_mean_ham_full;

    snew(diapop_old, ct->dim);
    snew(diapop, ct->dim);
    snew(ddiapop, ct->dim);
    snew(ekin, ct->dim);

// set up initial condition
    if (ct->tfs_initialization_step)
    {
       printf("Input wavefunction is:\n");
       for (j = 0; j < ct->dim; j++)
           printf("%f ", ct->wf[j]);
       printf("\n");

       for (i = 0; i < SQR(ct->dim); i++)
       {
           ct->tfl_old_ham[i] = *(ct->hamiltonian[0]+i);
           ham[i] = ct->tfl_old_ham[i];
       }

// diagonalize this Hamiltonian
       dsyev(ct->dim, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);

       for (i = 0; i < ct->dim; i++)
           for (j = 0; j < ct->dim; j++)
               ct->tfs_vector[i][j] = ham[i * ct->dim + j];

       for (i = 0; i < 2*ct->dim; i++)
       {
           ct->tfs_popul[i] = 0.0;
           ct->tfs_diab[i] = 0.0;
       }

       for (i = 0; i < ct->dim; i++)
       {
           ct->tfs_diab[i] = ct->wf[i];
           if (ct->wf[i] == 1.0)
              ct->surface = i;
       }

       for (i = 0; i < ct->dim; i++)
           ct->tfs_popul[i] = ct->tfs_vector[i][ct->surface];

       ct->tfs_initialization_step = 0;

       return 0;
    }

// calculate the diabatic populations at the previous time step
    for (i = 0; i < ct->dim; i++)
        diapop_old[i]  = SQR(ct->tfs_diab[i]) + SQR(ct->tfs_diab[i+ct->dim]);

// set up the old and new Hamiltonian for linear interpolation
    for (i = 0; i < ct->dim; i++)
        for (j = 0; j < ct->dim; j++)
        {
            old_ham[j+i*ct->dim] = ct->tfl_old_ham[i*ct->dim+j];
            new_ham[j+i*ct->dim] = ct->hamiltonian[i][j];
        }

/*    for (m= 0; m < ct->dim; m++)
    {
        for (n = 0; n < ct->dim; n++)
        {
            ct->fo_overlap[m][n] = 0.0;
            for (i = 0; i < ct->dim; i++)
            {
                for (j = 0; j < ct->dim; j++)
                {
                   for (k = 0; k < dftb->phase2.norb; k++)
                   {
                       for (l = 0; l < dftb->phase2.norb; l++)
                       {
                           ct->fo_overlap[m][n] += dftb->orthogo.sij_old[i+m*ct->dim]*dftb->orthogo.sij[j+n*ct->dim]
                                         *dftb->orthogo.fmo_old[k][i]*dftb->orthogo.fmo[l][j]
                                         *0.5*(dftb->orthogo.sao_old[k][l]+dftb->orthogo.sao[k][l]);
                       }
                   }

                }
            }
        }
    }
    for (m= 0; m < ct->dim; m++)
    {
        for (n = 0; n < ct->dim; n++)
        {
            ct->fo_tdc[m][n] = (ct->fo_overlap[m][n]- ct->fo_overlap[n][m])/2.0/ct->rk_timestep;
//            printf("%d %d %15.5f  %15.5f  %15.5f  \n", m, n, ct->fo_overlap[m][n], ct->fo_overlap[n][m], ct->fo_tdc[m][n]);
        }
    }
*/
    for (i = 0; i < ct->dim; i++)
        ct->surf_prob[i] = 0.0;

// solve time-dependent Schroedinger equation
    nstep = (int)(ct->rk_timestep/tdse_timestep);
    rk_timestep0 = ct->rk_timestep;
    ct->rk_timestep = tdse_timestep;

    for (istep = 1; istep < nstep+1; istep++)
    {
// linear interpolation for Hamiltonian
        cur_weight = (double)istep/(double)nstep;
        for (i = 0; i < ct->dim; i++)
            for (j = 0; j < ct->dim; j++)
            {
                ham_full[i][j] = (1.0 - cur_weight)*ct->tfl_old_ham[i*ct->dim+j] + cur_weight*ct->hamiltonian[i][j];
//                printf("%d  %d  %f  %f \n", i, j, ham_full[i][j], ct->tfl_mean_ham_full[i][j]);
            }
// propagate the TDSE in the diabatic representation
//        RK5th(ct);
        do_rksuite_diab(ct);

/*        for (i = 0; i < SQR(ct->dim); i++)
            ham[i] = (1.0 - cur_weight)*old_ham[i] + cur_weight*new_ham[i];

        for (i = 0; i < ct->dim; i++)
            for (j = 0; j < ct->dim; j++)
                ct->surf_prob[i] += (2.0*ham[i*ct->dim+j]*(ct->tfs_diab[i]*ct->tfs_diab[ct->dim+j]-ct->tfs_diab[j]*ct->tfs_diab[ct->dim+i]))*ct->rk_timestep/diapop_old[ct->surface];
*/
    }
    ct->rk_timestep = rk_timestep0;

// record the adiabatic populations at new time
    ddiapopsum = 0.0;
    for (i = 0; i < ct->dim; i++)
    {
        diapop[i] = SQR(ct->tfs_diab[i])+SQR(ct->tfs_diab[i+ct->dim]);

        ddiapop[i] = diapop[i] - diapop_old[i];

        if (ddiapop[i] < 0.0)
           ddiapopsum += ddiapop[i];
    }

// hopping probability
    for (i = 0; i < ct->dim; i++)
    {
        if (ddiapop[ct->surface] < 0.0 && ddiapop[i] > 0.0)
           ct->surf_prob[i] = ddiapop[i]/diapop_old[ct->surface]*ddiapop[ct->surface]/ddiapopsum;
        else
           ct->surf_prob[i] = 0.0;
    }

/*    for (i = 0; i < ct->dim; i++)
        ct->surf_prob[i] = ddiapop[i]/diapop_old[ct->surface];
*/
    for (i = 0; i < ct->dim; i++)
        printf("%2d  %15.5f  %15.5f \n", i, ddiapop[i]/diapop_old[ct->surface], ct->surf_prob[i]);

// kinetic energy of QM zone
    for (k = 0; k < ct->sites; k++)
        ekin[k] = 0.0;

    ekintot = 0.0;
    for (iatom = 0; iatom < mdatoms->homenr; iatom++)
    {
        for (k = 0; k < ct->sites; k++)
        {
            for (l = 0; l < ct->site[k].atoms; l++)
            {
                if (iatom == ct->site[k].atom[l])
                {
                   for (m = 0; m < DIM; m++)
                   {
                       ekintot += 0.5*dftb->phase1[k].mass[l]*AMU_TO_AU*SQR(state_global->v[iatom][m]*NM_TO_BOHR/PS_TO_AU);
                       ekin[k] += 0.5*dftb->phase1[k].mass[l]*AMU_TO_AU*SQR(state_global->v[iatom][m]*NM_TO_BOHR/PS_TO_AU);
                   }
                }
            }
        }
    }

    printf("Kinetic energy is %f %f eV\n", ekintot*27.2114, BOLTZMANN_HARTREE_KELVIN);

// stochastic algorithm
    probsum = 0.0;
    random_number = drand48();
    for (i = 0; i < ct->dim; i++)
    {
        if (ct->surf_prob[i] > 0.0 && i != ct->surface)
        {
           probsum += ct->surf_prob[i];

           if (random_number < probsum)
           {
              ekintmp = ekin[i] + ekin[ct->surface];
              de = ct->hamiltonian[i][i] - ct->hamiltonian[ct->surface][ct->surface];

              printf("Two-site kinetic energy is %f eV \n", ekintmp*27.2114);
              printf("Potential energy difference is %f eV \n", de*27.2114);

              if (1.0 - de/ekintmp < 0.0)
              {
                 printf("Frustated hop! \n");
                 fprintf(f,"%8.2f %10d %17d %21f %19f Frustrated \n",(step*ct->rk_timestep)*1000./PS_TO_AU, ct->surface, i, ct->surf_prob[i], random_number);
                 break;
              }
              else
              {

                 kai = sqrt(1.0 - de/ekintmp);

                 ekintot -= de;

                 for (iatom = 0; iatom < mdatoms->homenr; iatom++)
                 {
                     for (k = 0; k < ct->sites; k++)
                     {
                         if (k == i || k == ct->surface)
                         {
                             for (l = 0; l < ct->site[k].atoms; l++)
                             {
                                 if (iatom == ct->site[k].atom[l])
                                 {
                                    for (m = 0; m < DIM; m++)
                                        state_global->v[iatom][m] *= kai;
                                 }
                             }
                         }
                     }
                 }
              }

              fprintf(f,"%8.2f %10d %17d %21f %19f \n",(step*ct->rk_timestep)*1000./PS_TO_AU, ct->surface, i, ct->surf_prob[i], random_number);

              ct->surface = i;

              break;
           }
        }
    }

    printf("Current Surface: %d\n", ct->surface);

// decoherence corrections
    dot_product = 0.0;
    for (i = 0; i < ct->dim; i++)
    {
        if (i != ct->surface)
        {
            decay_time = (1.0 + 0.1/ekintot)/fabs(ct->hamiltonian[i][i] - ct->hamiltonian[ct->surface][ct->surface]);

            ct->tfs_diab[i]         *= exp(-ct->rk_timestep/decay_time);
            ct->tfs_diab[i+ct->dim] *= exp(-ct->rk_timestep/decay_time);

            dot_product              += SQR(ct->tfs_diab[i])+SQR(ct->tfs_diab[i+ct->dim]);
        }
    }
    curr_surf_fac = sqrt((1.0 - dot_product)/(SQR(ct->tfs_diab[ct->surface]) + SQR(ct->tfs_diab[ct->dim+ct->surface])));

    ct->tfs_diab[ct->surface]         *= curr_surf_fac;
    ct->tfs_diab[ct->surface+ct->dim] *= curr_surf_fac;

// new state
    for (i = 0; i < ct->dim; i++)
        ct->wf[i] = 0.0;
    ct->wf[ct->surface] = 1.0;

    for (i = 0; i < SQR(ct->dim); i++)
        ham[i] = new_ham[i];

    dsyev(ct->dim, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);

// calculate expansion coefficients in the adiabatic basis
    for (i = 0; i < ct->dim; i++)
        for (j = 0; j < ct->dim; j++)
            ct->tfs_vector[i][j] = ham[i * ct->dim + j];

    for (i = 0; i < ct->dim; i++)
    {
        ct->tfs_popul[i] = 0.0; ct->tfs_popul[i+ct->dim] = 0.0;
        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_popul[i]         += ct->tfs_vector[i][j]*ct->tfs_diab[j];
            ct->tfs_popul[i+ct->dim] += ct->tfs_vector[i][j]*ct->tfs_diab[j+ct->dim];

        }
    }

//store hamiltonian as old_ham
    for (i = 0; i < SQR(ct->dim); i++)
        ct->tfl_old_ham[i] = *(ct->hamiltonian[0]+i);

    return 0;
}

#if GMX_MPI
int do_fssh_johen(int step, charge_transfer_t *ct, dftb_t *dftb, t_state *state_global, t_mdatoms *mdatoms, FILE *f,
                  MPI_Comm ct_mpi_comm, int ct_mpi_rank, int ct_mpi_size)
#else
int do_fssh_johen(int step, charge_transfer_t *ct, dftb_t *dftb, t_state *state_global, t_mdatoms *mdatoms, FILE *f)
#endif
{
// Declaration:
//
// 1. a correction is implemented for decoherence correction-induced spurious long-range charge transfercorrection for the decoherence;
// 2. self-consistent algorithm is applied for avoided crossings;
// 3. state-tracking algorithm is implemented for trivial corrssing detection;
//
//                                          -----------Weiwei Xie on Feb. 2019
// Ref: J. Phys. Chem. Lett 9, 3116 (2018)



    int    i, j, k, l, m, n, tmp, prev_surf, old_surf, hopsurf, new_surf, iatom, istep, nstep, *ind, indmax, *indtmp, *indtmp1, *map, do_lambda_i_tmp;
    double **grad0, *ham, *old_ham, **ham_full, *new_ham, *adiapop_old, *adiapop, *initprob, *probtmp1, *probtmp2, *diapop, *diatmp;
    double tdse_timestep = 1.e-5*PS_TO_AU, diapopmax = 0.999;
    double rk_timestep0, probsum, random_number, cur_weight, random_number1;
    double popul_norm, overlapmax, overlaptmp1, overlaptmp2, adiapoptmp0, diapoptot, diapoptot0;
    double nc, kai1, kai2, kai, delta, aa, bb, cc, dd, ee, de;
    double decay_time, curr_surf_fac, trivial_cross_type;
    double dot_product, ekintot;
    double delta_min;


    old_ham  = ct->tfl_old_ham_adiab;
    new_ham  = ct->hamiltonian_adiab;
    ham      = ct->tfl_mean_ham;
    ham_full = ct->tfl_mean_ham_full;

    snew(diapop, ct->dim);
    snew(diatmp, 2*ct->dim);
    snew(adiapop_old, ct->dim);
    snew(adiapop, ct->dim);
    snew(initprob, ct->dim);
    snew(probtmp1, ct->dim);
    snew(probtmp2, ct->dim);

    snew(ind, ct->dim);
    snew(indtmp, ct->dim);
    snew(indtmp1, ct->dim);
    snew(map, ct->dim);

    if (ct->do_lambda_i == 3)
    {
       snew(grad0, dftb->phase2.nn);
       for (j = 0; j < dftb->phase2.nn; j++)
           snew(grad0[j], 3);
    }

// set up initial surface
    if (ct->tfs_initialization_step)
    {
       printf("Input wavefunction is:\n");
       for (j = 0; j < ct->dim; j++)
           printf("%f ", ct->wf[j]);
       printf("\n");

       for (i = 0; i < SQR(ct->dim); i++)
       {
           ct->tfl_old_ham[i] = *(ct->hamiltonian[0]+i);
           ham[i] = ct->tfl_old_ham[i];
       }

// diagonalize this Hamiltonian
       dsyev(ct->dim, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);

// record eigenvector
       for (i = 0; i < ct->dim; i++)
           for (j = 0; j < ct->dim; j++)
               ct->tfs_vector[i][j] = ham[i * ct->dim + j];

// determine inital surface using eigenvector as a criterion
       probsum = 0.;
       random_number = drand48();
       for (i = 0; i < ct->dim; i++)
       {
           initprob[i] = 0.0;

           for (j = 0; j < ct->dim; j++)
               initprob[i] += ct->tfs_vector[i][j]*ct->wf[j];

           initprob[i] = SQR(initprob[i]);

           probsum += initprob[i];
           if (random_number < probsum)
           {
               ct->surface = i;
               break;
           }
       }

       for (i = 0; i < 2*ct->dim; i++)
       {
           ct->tfs_popul[i] = 0.0; ct->tfs_diab[i]  = 0.0;
       }

       for (i = 0; i < ct->dim; i++)
           ct->tfs_diab[i] = ct->wf[i];

       ct->tfs_popul[ct->surface] = 1.0;

       ct->tfs_initialization_step = 0;

       return 0;
    }

// set up the old and new Hamiltonian for linear interpolation
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            old_ham[i*ct->dim+j] = ct->tfl_old_ham[i*ct->dim+j];
            new_ham[i*ct->dim+j] = ct->hamiltonian[i][j];
        }
    }

// store the vector at preceding time step
/*    for (i = 0; i < ct->dim; i++)
        for (j = 0; j < ct->dim; j++)
            ct->tfs_vector_old2[i][j] = ct->tfs_vector[i][j];
*/
// solve time-dependent Schroedinger equation
    nstep = (int)(ct->rk_timestep/tdse_timestep);
    rk_timestep0 = ct->rk_timestep;
    ct->rk_timestep = tdse_timestep;
    printf("Time step for TDSE: %f fs \n", ct->rk_timestep/PS_TO_AU*1000);

    for (istep = 1; istep < nstep+1; istep++)
    {
// linear interpolation for Hamiltonian
        cur_weight = (double)istep/(double)nstep;
        for (i = 0; i < ct->dim; i++)
            for (j = 0; j < ct->dim; j++)
                ham_full[i][j] = (1.0 - cur_weight)*ct->tfl_old_ham[i*ct->dim+j] + cur_weight*ct->hamiltonian[i][j];
// propagate the TDSE in the diabatic representation
        do_rksuite_diab(ct);
    }
    ct->rk_timestep = rk_timestep0;

// diagonalize this Hamiltonian
    for (i = 0; i < SQR(ct->dim); i++)
        ham[i] = new_ham[i];

    dsyev(ct->dim, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);

// store the eigenvectors at time t-dt and time t
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_vector_old[i][j] = ct->tfs_vector[i][j];
            ct->tfs_vector[i][j]     = ham[i*ct->dim + j];
        }
    }

// adiabatic wave function overlaps between time t-dt and t
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_overlap[i][j] = 0.0;

            for (k = 0; k < ct->dim; k++)
                ct->tfs_overlap[i][j] += ct->tfs_vector_old[i][k]*ct->tfs_vector[j][k];
        }
    }

// reassignment of states according to the maximum wave function overlaps
    for (i = 0; i < ct->dim; i++)
    {
        overlapmax = 0.0;

        for (j = 0; j < ct->dim; j++)
        {
            if (overlapmax < fabs(ct->tfs_overlap[j][i]))
            {
               overlapmax = fabs(ct->tfs_overlap[j][i]);
               indtmp[i] = j;
            }
        }
    }

// reassignment of states according to the min-cost algorithm
//    assignprob(ct, indtmp);

// state tracking
    for (i = 0; i < ct->dim; i++)
    {
        map[i] = 0; ind[i] = -1;
    }

    for (i = 0; i < ct->dim; i++)
    {
        if (fabs(ct->tfs_overlap[indtmp[i]][i]) > 0.9)
        {
           ind[i] = indtmp[i];
           map[ind[i]] = 1;
        }
    }

    for (i = 0; i < ct->dim; i++)
    {
        if (ind[i] == -1)
        {
           for (j = 0; j < ct->dim; j++)
           {
               if (! map[j])
               {
                  ind[i] = j; map[j] = 1;
                  break;
               }
           }
        }
    }

// check signs
    for (i = 0; i < ct->dim; i++)
    {
        j = ind[i];

        if (ct->tfs_overlap[j][i] < 0.0)
        {
            for (k = 0; k < ct->dim; k++)
            {
                ct->tfs_vector[i][k]  *= -1.0;
                ct->tfs_overlap[k][i] *= -1.0;
            }
        }
    }

// reassign the active state
    prev_surf = ct->surface;
    if (ct->surface != ind[ct->surface])
    {
       for (i = 0; i < ct->dim; i++)
       {
           if (ct->surface == ind[i])
           {
              printf("Changed from surface %d to surface %d \n", ct->surface, i);
              ct->surface = i;
              break;
           }
       }
    }

// calculate the adiabatic populations at the previous time step
    for (i = 0; i < ct->dim; i++)
        adiapop_old[i] = SQR(ct->tfs_popul[i]) + SQR(ct->tfs_popul[i+ct->dim]);

// hopping probabilities based on Tully's fewest-switches anzatz
    for (i = 0; i < ct->dim; i++)
    {
        j = ind[i];
//     nonadibatic couplings based on Hammes-Schiffer and Tully, JCP 101, 4657 (1994)
//        nc = (ct->tfs_overlap[prev_surf][i] - ct->tfs_overlap[j][ct->surface])/2.0/ct->rk_timestep;
//     ndibatic couplings based norm-preserving interpolation by Benjamin, JPCL 5, 2351 (2014)
        overlaptmp1 = sqrt(fabs(1.0 - SQR(ct->tfs_overlap[j][i]) - SQR(ct->tfs_overlap[prev_surf][i])));
        overlaptmp2 = (- ct->tfs_overlap[j][ct->surface]*ct->tfs_overlap[j][i] - ct->tfs_overlap[prev_surf][ct->surface]*ct->tfs_overlap[prev_surf][i])/overlaptmp1;

        aa = - sin(acos(ct->tfs_overlap[j][i]) - asin(ct->tfs_overlap[j][ct->surface]))/(acos(ct->tfs_overlap[j][i]) - asin(ct->tfs_overlap[j][ct->surface]));
        bb =   sin(acos(ct->tfs_overlap[j][i]) + asin(ct->tfs_overlap[j][ct->surface]))/(acos(ct->tfs_overlap[j][i]) + asin(ct->tfs_overlap[j][ct->surface]));
        cc =   sin(acos(ct->tfs_overlap[prev_surf][ct->surface]) - asin(ct->tfs_overlap[prev_surf][i]))/(acos(ct->tfs_overlap[prev_surf][ct->surface]) - asin(ct->tfs_overlap[prev_surf][i]));
        dd =   sin(acos(ct->tfs_overlap[prev_surf][ct->surface]) + asin(ct->tfs_overlap[prev_surf][i]))/(acos(ct->tfs_overlap[prev_surf][ct->surface]) + asin(ct->tfs_overlap[prev_surf][i]));
        ee = 2.0*asin(overlaptmp1*(overlaptmp1*overlaptmp2*asin(overlaptmp1) + (sqrt((1.0 - SQR(overlaptmp1))*(1.0 - SQR(overlaptmp2))) - 1.0)*asin(overlaptmp2)))/(SQR(asin(overlaptmp1)) - SQR(asin(overlaptmp2)));

        nc = (acos(ct->tfs_overlap[j][i])*(aa + bb) + asin(ct->tfs_overlap[prev_surf][i])*(cc + dd) + ee)/2.0/ct->rk_timestep;

        probtmp1[i] = (2.0*(ct->tfs_popul[prev_surf]*ct->tfs_popul[j]+ct->tfs_popul[ct->dim+prev_surf]*ct->tfs_popul[ct->dim+j])*nc)*ct->rk_timestep/adiapop_old[prev_surf];

        if (probtmp1[i] < 0.0 || i == ct->surface)
           probtmp1[i] = 0.0;
    }

// calculate expansion coefficients in the adiabatic basis
    for (i = 0; i < ct->dim; i++)
    {
        ct->tfs_popul[i] = 0.0; ct->tfs_popul[i+ct->dim] = 0.0;

        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_popul[i]         += ct->tfs_vector[i][j]*ct->tfs_diab[j];
            ct->tfs_popul[i+ct->dim] += ct->tfs_vector[i][j]*ct->tfs_diab[j+ct->dim];
        }
    }

// record the adiabatic populations at new time
    for (i = 0; i < ct->dim; i++)
        adiapop[i] = SQR(ct->tfs_popul[i])+SQR(ct->tfs_popul[i+ct->dim]);

// hopping probabilities based on populations flux
/*    for (i = 0; i < ct->dim; i++)
    {
        j = ind[i];

        probtmp2[i] = (adiapop[i] - adiapop_old[j])/adiapop_old[prev_surf];

        if (probtmp2[i] < 0.0 || i == ct->surface)
           probtmp2[i] = 0.0;
    }
*/
    new_surf = ct->surface;
    for (i = 0; i < ct->dim; i++)
    {
        if (prev_surf == indtmp[i])
        {
           new_surf = i;
           break;
        }
    }

// SC-FSSH, avoid the numerical instability of the NCs
    if (ct->surface != new_surf)
    {
       probtmp1[new_surf] = (adiapop_old[ct->surface] - adiapop[ct->surface])/adiapop_old[ct->surface];

       for (i = 0; i < ct->dim; i++)
       {
           if (i != new_surf)
           {
              probtmp1[new_surf] -= probtmp1[i];
           }
       }

       printf("SC is actived: %d  %d  %10f  %10f\n", ct->surface, new_surf, (adiapop_old[ct->surface] - adiapop[ct->surface])/adiapop_old[ct->surface], probtmp1[new_surf]);
    }

// define how to do adjustment for velocities
// 0: isotropic rescaling; 1: along NAC; 2: Boltzmann corrections.
    if (ct->jobtype == cteBCJFSSH)
    {
       dftb->do_ncv = 2;
// Boltzmann factor corrected hopping probability
       for (i = ct->surface+1; i < ct->dim; i++)
           probtmp1[i] *= exp(-(ct->ev_adiab[i] - ct->ev_adiab[ct->surface])/ct->fermi_kt);
    }
    else
    {
       dftb->do_ncv = 1;
    }

// print out the surface probabilities
    for (i = 0; i < ct->dim; i++)
        printf("%2d  %2d  %2d  %15.5f \n", i, indtmp[i], ind[i], probtmp1[i]);

// stochastic algorithm
    probsum = 0.0;
    random_number = drand48();
    for (i = 0; i < ct->dim; i++)
    {
        ct->surf_prob[i] = probtmp1[i];

        if (ct->surf_prob[i] > 0.0 && i != ct->surface)
        {
           probsum += ct->surf_prob[i];

           if (random_number < probsum)
           {
              hopsurf = i;

// calculate the nonadibatic coupling vector for velocity adjustment
              if (dftb->do_ncv == 1)
              {

                 ct->ncv_isurf = ct->surface; ct->ncv_jsurf = hopsurf;

                 if (ct->do_lambda_i == 3)
                 {
                    for (iatom = 0; iatom < mdatoms->homenr; iatom++)
                    {
                        n = 0;
                        for (k = 0; k < ct->sites; k++)
                        {
                            for (l = 0; l < ct->site[k].atoms; l++)
                            {
                                if (iatom == ct->site[k].atom[l])
                                {
                                   for (m = 0; m < DIM; m++)
                                   {
                                       grad0[n][m] = dftb->phase2.grad[n][m];
                                   }
                                }
                                n++;
                            }
                        }
                    }
                 }

                 do_lambda_i_tmp = ct->do_lambda_i;
                 ct->do_lambda_i = 3;

                 if(ct->hamiltonian_type == 0)
                 {
#if GMX_MPI
                   get_MM_params(ct, dftb, ct_mpi_comm, ct_mpi_rank, ct_mpi_size);
#else
                   get_MM_params(ct, dftb);
#endif
                 }
                 else if (ct->hamiltonian_type == 1)
                 {
                   get_nn_hamiltonian(ct,state_global, mdatoms, dftb);
                 }


                 ct->do_lambda_i = do_lambda_i_tmp;

                 aa = 0.0; bb = 0.0;
                 cc = ct->ev_adiab[ct->ncv_jsurf] - ct->ev_adiab[ct->ncv_isurf];

                 for (iatom = 0; iatom < mdatoms->homenr; iatom++)
                 {
                     n = 0;
                     for (k = 0; k < ct->sites; k++)
                     {
                         for (l = 0; l < ct->site[k].atoms; l++)
                         {
                             if (iatom == ct->site[k].atom[l])
                             {
                                for (m = 0; m < DIM; m++)
                                {

                                    ct->tfs_ncv[k][l][m] = (dftb->phase1[k].grad[l][m] + dftb->phase2.grad[n][m])/cc *
                                                           (ct->is_hole_transfer ? 1.0 : -1.0);

                                    if (ct->do_lambda_i == 3)
                                    {
                                       dftb->phase2.grad[n][m] = grad0[n][m];
                                    }

                                    aa += SQR(ct->tfs_ncv[k][l][m])/2./dftb->phase1[k].mass[l]/AMU_TO_AU;
                                    bb += ct->tfs_ncv[k][l][m]*state_global->v[iatom][m]*NM_TO_BOHR/PS_TO_AU;

                                }
                             }
                             n++;
                         }
                     }
                 }

// adjust velocity along the direction of nonadiabtic coupling vector
                 delta = SQR(bb) - 4.*aa*cc;
                 if (delta >= 0.0)
                 {
                    kai1 = (- bb + sqrt(delta))/2./aa;
                    kai2 = (- bb - sqrt(delta))/2./aa;

                    kai = (fabs(kai1) < fabs(kai2) ? kai1 : kai2);

                    printf("Hop from %d to %d succeeds! Adjust velocity along NCV with kai = %8.5f \n", ct->surface, hopsurf, kai);
                 }
                 else
                 {
                    kai = - bb/aa;
                    hopsurf = ct->surface;

                    printf("Hop from %d to %d fails! Reverse velocity along NCV with kai = %8.5f \n", ct->surface, i, kai);
                 }

// velocity adjustment
                 for (iatom = 0; iatom < mdatoms->homenr; iatom++)
                 {
                     for (k = 0; k < ct->sites; k++)
                     {
                         for (l = 0; l < ct->site[k].atoms; l++)
                         {
                             if (iatom == ct->site[k].atom[l])
                             {
                                for (m = 0; m < DIM; m++)
                                    state_global->v[iatom][m] += kai*ct->tfs_ncv[k][l][m]/dftb->phase1[k].mass[l]/AMU_TO_AU/NM_TO_BOHR*PS_TO_AU;
                             }
                         }
                     }
                 }
              }
              else if (dftb->do_ncv == 0)
              {
// compute kinetic energy
                 ekintot = 0.0;
                 for (iatom = 0; iatom < mdatoms->homenr; iatom++)
                 {
                     for (k = 0; k < ct->sites; k++)
                     {
                         for (l = 0; l < ct->site[k].atoms; l++)
                         {
                             if (iatom == ct->site[k].atom[l])
                             {
                                for (m = 0; m < DIM; m++)
                                    ekintot += 0.5*dftb->phase1[k].mass[l]*AMU_TO_AU*SQR(state_global->v[iatom][m]*NM_TO_BOHR/PS_TO_AU);
                             }
                         }
                     }
                 }
                 printf("Kinetic energy is %f\n", ekintot);

                 de = ct->ev_adiab[i] - ct->ev_adiab[ct->surface];

                 if (1.0 - de/ekintot < 0.0)
                 {
                    printf("Frustated hop! \n");
                    fprintf(f,"%8.2f %10d %17d %21f %19f Frustrated \n",(step*ct->rk_timestep)*1000./PS_TO_AU, ct->surface, i, ct->surf_prob[i], random_number);
                    break;
                 }
                 else
                 {
                    kai = sqrt(1.0 - de/ekintot);
                    printf("kai is %f \n", kai);

                    ekintot -= de;

                    for (iatom = 0; iatom < mdatoms->homenr; iatom++)
                    {
                        for (k = 0; k < ct->sites; k++)
                        {
                            for (l = 0; l < ct->site[k].atoms; l++)
                            {
                                if (iatom == ct->site[k].atom[l])
                                {
                                   for (m = 0; m < DIM; m++)
                                       state_global->v[iatom][m] *= kai;
                                }
                            }
                        }
                    }
                 }

              }

              if (hopsurf != ct->surface)
                 fprintf(f,"%8.2f %10d %17d %21f %19f \n",(step*ct->rk_timestep)*1000./PS_TO_AU, ct->surface, hopsurf, ct->surf_prob[i], random_number);

              ct->surface = hopsurf;

              break;
           }
        }
    }

    dftb->do_ncv = 0;

    printf("Current Surface: %d\n", ct->surface);

// compute kinetic energy
    ekintot = 0.0;
    for (iatom = 0; iatom < mdatoms->homenr; iatom++)
    {
        for (k = 0; k < ct->sites; k++)
        {
            for (l = 0; l < ct->site[k].atoms; l++)
            {
                if (iatom == ct->site[k].atom[l])
                {
                   for (m = 0; m < DIM; m++)
                       ekintot += 0.5*dftb->phase1[k].mass[l]*AMU_TO_AU*SQR(state_global->v[iatom][m]*NM_TO_BOHR/PS_TO_AU);
                }
            }
        }
    }
    printf("Kinetic energy is %f\n", ekintot);

// decoherence correction
    dot_product = 0.0;
    for (i = 0; i < ct->dim; i++)
    {
        if (i != ct->surface)
        {
            decay_time = (1.0 + 0.1/ekintot)/fabs(ct->ev_adiab[i] - ct->ev_adiab[ct->surface]);

            ct->tfs_popul[i]         *= exp(-ct->rk_timestep/decay_time);
            ct->tfs_popul[i+ct->dim] *= exp(-ct->rk_timestep/decay_time);
            dot_product              += SQR(ct->tfs_popul[i])+SQR(ct->tfs_popul[i+ct->dim]);
        }
    }
    curr_surf_fac = sqrt((1.0 - dot_product)/(SQR(ct->tfs_popul[ct->surface]) + SQR(ct->tfs_popul[ct->dim+ct->surface])));
    ct->tfs_popul[ct->surface]         *= curr_surf_fac;
    ct->tfs_popul[ct->surface+ct->dim] *= curr_surf_fac;

// recalc. diab. from adiab after decoherence corr.
    for (i = 0; i < ct->dim; i++)
    {
        diatmp[i] = ct->tfs_diab[i]; diatmp[i+ct->dim] = ct->tfs_diab[i+ct->dim];

        ct->tfs_diab[i] = 0.0; ct->tfs_diab[i+ct->dim] = 0.0;

        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_diab[i]         += ct->tfs_popul[j]*ct->tfs_vector[j][i];
            ct->tfs_diab[i+ct->dim] += ct->tfs_popul[j+ct->dim]*ct->tfs_vector[j][i];
        }
    }

// johen's corrections
    for (i = 0; i < ct->dim; i++)
        diapop[i] = SQR(ct->tfs_diab[i])+SQR(ct->tfs_diab[i+ct->dim]);

    diapoptot = 0.0;
    for (i = 0; i < ct->dim; i++)
    {
        diapoptot += diapop[i];

        if (diapoptot > diapopmax)
        {
           n = i; diapoptot -= diapop[i];
           break;
        }
    }

    diapoptot0 = 0.0;
    for (i = n; i < ct->dim; i++)
    {
        ct->tfs_diab[i]         = diatmp[i];
        ct->tfs_diab[i+ct->dim] = diatmp[i+ct->dim];

        diapoptot0 += SQR(ct->tfs_diab[i])+SQR(ct->tfs_diab[i+ct->dim]);
    }

    for (i = 0; i < n; i++)
    {
        ct->tfs_diab[i]         *= sqrt((1.0 - diapoptot0)/diapoptot);
        ct->tfs_diab[i+ct->dim] *= sqrt((1.0 - diapoptot0)/diapoptot);
    }

// calculate expansion coefficients in the adiabatic basis
    for (i = 0; i < ct->dim; i++)
    {
        ct->tfs_popul[i] = 0.0; ct->tfs_popul[i+ct->dim] = 0.0;

        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_popul[i]         += ct->tfs_vector[i][j]*ct->tfs_diab[j];
            ct->tfs_popul[i+ct->dim] += ct->tfs_vector[i][j]*ct->tfs_diab[j+ct->dim];
        }
    }

// charge wave function
    for (i = 0; i < ct->dim; i++)
        ct->wf[i] = ct->tfs_vector[ct->surface][i];

//store hamiltonian as old_ham
    for (i = 0; i < SQR(ct->dim); i++)
        ct->tfl_old_ham[i] = *(ct->hamiltonian[0]+i);

    sfree(diapop);
    sfree(diatmp);
    sfree(adiapop_old);
    sfree(adiapop);
    sfree(initprob);
    sfree(probtmp1);
    sfree(probtmp2);

    sfree(ind);
    sfree(indtmp);
    sfree(indtmp1);
    sfree(map);

    if (ct->do_lambda_i == 3)
    {
       for (i = 0; i < dftb->phase2.nn; i++)
           sfree(grad0[i]);
    }

    return 0;
}

#if GMX_MPI
int do_self_consistent_restricted_decoherence_fssh(int step, charge_transfer_t *ct, dftb_t *dftb, t_state *state_global, t_mdatoms *mdatoms, FILE *f,
                                                   MPI_Comm ct_mpi_comm, int ct_mpi_rank, int ct_mpi_size)
#else
int do_self_consistent_restricted_decoherence_fssh(int step, charge_transfer_t *ct, dftb_t *dftb, t_state *state_global, t_mdatoms *mdatoms, FILE *f)
#endif
{
// Declaration:
//
// 1. self-consistent algorithm is implemented for weak avoided crossings;
// 2. restricted decoherence correction only for the adiabatic population larger than a threshold;
//
//                                          -----------Weiwei Xie on Mar. 2019
// Ref: J.Chem.Phys. 148, 104106 (2018)


    int    i, j, k, l, m, n, tmp, prev_surf, old_surf, hopsurf, new_surf, surftmp, iatom, istep, nstep, *ind, indmax, *indtmp, *indtmp1, *map, newline;
    double *ham, *old_ham, **ham_full, *new_ham, *adiapop_old, *adiapop, *initprob, *probtmp1, *probtmp2, *diapop, *diatmp;
    double tdse_timestep = 1.e-5*PS_TO_AU;
    double rk_timestep0, probsum, random_number, cur_weight;
    double popul_norm, overlapmax, overlaptmp1, overlaptmp2, adiapoptmp0, probtmp0, probtmp, diapoptmp;
    double nc, kai1, kai2, kai, delta, aa, bb, cc, dd, ee;
    double decay_time, curr_surf_fac, trivial_cross_type;
    double dot_product, ekintot;
    double delta_min;


    old_ham  = ct->tfl_old_ham_adiab;
    new_ham  = ct->hamiltonian_adiab;
    ham      = ct->tfl_mean_ham;
    ham_full = ct->tfl_mean_ham_full;

    snew(diapop, ct->dim);
    snew(diatmp, 2*ct->dim);
    snew(adiapop_old, ct->dim);
    snew(adiapop, ct->dim);
    snew(initprob, ct->dim);
    snew(probtmp1, ct->dim);
    snew(probtmp2, ct->dim);

    snew(ind, ct->dim);
    snew(indtmp, ct->dim);
    snew(indtmp1, ct->dim);
    snew(map, ct->dim);

// set up initial surface
    if (ct->tfs_initialization_step)
    {
       printf("Input wavefunction is:\n");
       for (j = 0; j < ct->dim; j++)
           printf("%f ", ct->wf[j]);
       printf("\n");

       for (i = 0; i < SQR(ct->dim); i++)
       {
           ct->tfl_old_ham[i] = *(ct->hamiltonian[0]+i);
           ham[i] = ct->tfl_old_ham[i];
       }

// diagonalize this Hamiltonian
       dsyev(ct->dim, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);

// record eigenvector
       for (i = 0; i < ct->dim; i++)
           for (j = 0; j < ct->dim; j++)
               ct->tfs_vector[i][j] = ham[i * ct->dim + j];

// determine inital surface using eigenvector as a criterion
       probsum = 0.;
       random_number = drand48();
       for (i = 0; i < ct->dim; i++)
       {
           initprob[i] = 0.;
           for (j = 0; j < ct->dim; j++)
               initprob[i] += ct->tfs_vector[i][j]*ct->wf[j];
           initprob[i] = SQR(initprob[i]);

           probsum += initprob[i];
           if (random_number < probsum)
           {
               ct->surface = i;
               break;
           }
       }

       for (i = 0; i < 2*ct->dim; i++)
       {
           ct->tfs_popul[i] = 0.0; ct->tfs_diab[i]  = 0.0;
       }

       for (i = 0; i < ct->dim; i++)
           ct->tfs_diab[i] = ct->wf[i];

       ct->tfs_popul[ct->surface] = 1.0;

       ct->tfs_initialization_step = 0;

       return 0;
    }

// set up the old and new Hamiltonian for linear interpolation
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            old_ham[i*ct->dim+j] = ct->tfl_old_ham[i*ct->dim+j];
            new_ham[i*ct->dim+j] = ct->hamiltonian[i][j];
        }
    }

// store the vector at preceding time step
/*    for (i = 0; i < ct->dim; i++)
        for (j = 0; j < ct->dim; j++)
            ct->tfs_vector_old2[i][j] = ct->tfs_vector[i][j];
*/
// solve time-dependent Schroedinger equation
    nstep = (int)(ct->rk_timestep/tdse_timestep);
    rk_timestep0 = ct->rk_timestep;
    ct->rk_timestep = tdse_timestep;
    printf("Time step for TDSE: %f fs \n", ct->rk_timestep/PS_TO_AU*1000);

    for (istep = 1; istep < nstep+1; istep++)
    {
// linear interpolation for Hamiltonian
        cur_weight = (double)istep/(double)nstep;
        for (i = 0; i < ct->dim; i++)
            for (j = 0; j < ct->dim; j++)
                ham_full[i][j] = (1.0 - cur_weight)*ct->tfl_old_ham[i*ct->dim+j] + cur_weight*ct->hamiltonian[i][j];
// propagate the TDSE in the diabatic representation
        do_rksuite_diab(ct);
    }
    ct->rk_timestep = rk_timestep0;

// diagonalize this Hamiltonian
    for (i = 0; i < SQR(ct->dim); i++)
        ham[i] = new_ham[i];

    dsyev(ct->dim, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);

// store the eigenvectors at time t-dt and time t
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_vector_old[i][j] = ct->tfs_vector[i][j];
            ct->tfs_vector[i][j]     = ham[i*ct->dim + j];
        }
    }

// adiabatic wave function overlaps between time t-dt and t
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_overlap[i][j] = 0.0;

            for (k = 0; k < ct->dim; k++)
                ct->tfs_overlap[i][j] += ct->tfs_vector_old[i][k]*ct->tfs_vector[j][k];
        }
    }

// reassignment of states according to the maximum wave function overlaps
/*    for (i = 0; i < ct->dim; i++)
    {
        overlapmax = 0.0;
        for (j = 0; j < ct->dim; j++)
        {
            if (overlapmax < fabs(ct->tfs_overlap[j][i]))
            {
               overlapmax = fabs(ct->tfs_overlap[j][i]);
               indtmp[i] = j;
            }
        }
    }
*/
// check signs
    for (i = 0; i < ct->dim; i++)
    {
        if (ct->tfs_overlap[i][i] < 0.0)
        {
            for (k = 0; k < ct->dim; k++)
            {
                ct->tfs_vector[i][k]  *= -1.0;
                ct->tfs_overlap[k][i] *= -1.0;
            }
        }
    }

// reassignment of states according to the min-cost algorithm
    assignprob(ct, ind);

// calculate the adiabatic populations at the previous time step
    for (i = 0; i < ct->dim; i++)
        adiapop_old[i] = SQR(ct->tfs_popul[i]) + SQR(ct->tfs_popul[i+ct->dim]);

// hopping probabilities based on Tully's fewest-switches anzatz
    for (i = 0; i < ct->dim; i++)
    {
//     nonadibatic couplings based on Hammes-Schiffer and Tully, JCP 101, 4657 (1994)
//        nc = (ct->tfs_overlap[prev_surf][i] - ct->tfs_overlap[j][ct->surface])/2.0/ct->rk_timestep;
//     ndibatic couplings based norm-preserving interpolation by Benjamin, JPCL 5, 2351 (2014)
        overlaptmp1 = sqrt(fabs(1.0 - SQR(ct->tfs_overlap[i][i]) - SQR(ct->tfs_overlap[ct->surface][i])));
        overlaptmp2 = (- ct->tfs_overlap[i][ct->surface]*ct->tfs_overlap[i][i] - ct->tfs_overlap[ct->surface][ct->surface]*ct->tfs_overlap[ct->surface][i])/overlaptmp1;

        aa = - sin(acos(ct->tfs_overlap[i][i]) - asin(ct->tfs_overlap[i][ct->surface]))/(acos(ct->tfs_overlap[i][i]) - asin(ct->tfs_overlap[i][ct->surface]));
        bb =   sin(acos(ct->tfs_overlap[i][i]) + asin(ct->tfs_overlap[i][ct->surface]))/(acos(ct->tfs_overlap[i][i]) + asin(ct->tfs_overlap[i][ct->surface]));
        cc =   sin(acos(ct->tfs_overlap[ct->surface][ct->surface]) - asin(ct->tfs_overlap[ct->surface][i]))/(acos(ct->tfs_overlap[ct->surface][ct->surface]) - asin(ct->tfs_overlap[ct->surface][i]));
        dd =   sin(acos(ct->tfs_overlap[ct->surface][ct->surface]) + asin(ct->tfs_overlap[ct->surface][i]))/(acos(ct->tfs_overlap[ct->surface][ct->surface]) + asin(ct->tfs_overlap[ct->surface][i]));
        ee = 2.0*asin(overlaptmp1*(overlaptmp1*overlaptmp2*asin(overlaptmp1) + (sqrt((1.0 - SQR(overlaptmp1))*(1.0 - SQR(overlaptmp2))) - 1.0)*asin(overlaptmp2)))/(SQR(asin(overlaptmp1)) - SQR(asin(overlaptmp2)));

        nc = (acos(ct->tfs_overlap[i][i])*(aa + bb) + asin(ct->tfs_overlap[ct->surface][i])*(cc + dd) + ee)/2.0/ct->rk_timestep;

        probtmp1[i] = (2.0*(ct->tfs_popul[ct->surface]*ct->tfs_popul[i]+ct->tfs_popul[ct->dim+ct->surface]*ct->tfs_popul[ct->dim+i])*nc)*ct->rk_timestep/adiapop_old[ct->surface];

        if (probtmp1[i] < 0.0 || i == ct->surface)
           probtmp1[i] = 0.0;
    }

// calculate expansion coefficients in the adiabatic basis
    for (i = 0; i < ct->dim; i++)
    {
        ct->tfs_popul[i] = 0.0; ct->tfs_popul[i+ct->dim] = 0.0;

        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_popul[i]         += ct->tfs_vector[i][j]*ct->tfs_diab[j];
            ct->tfs_popul[i+ct->dim] += ct->tfs_vector[i][j]*ct->tfs_diab[j+ct->dim];
        }
    }

// record the adiabatic populations at new time
    for (i = 0; i < ct->dim; i++)
        adiapop[i] = SQR(ct->tfs_popul[i])+SQR(ct->tfs_popul[i+ct->dim]);

// SC-FSSH, avoid the numerical instability of the NCs
    prev_surf = ct->surface;
    for (i = 0; i < ct->dim; i++)
    {
        if (prev_surf == ind[i])
        {
           new_surf = i;
           break;
        }
    }

    newline = 0;
    if (ct->surface != new_surf)
    {
       probtmp1[new_surf] = (adiapop_old[ct->surface] - adiapop[ct->surface])/adiapop_old[ct->surface];

       for (i = 0; i < ct->dim; i++)
       {
           if (i != new_surf)
           {
              probtmp1[new_surf] -= probtmp1[i];
           }
       }

//       printf("SC is actived: %d  %d  %10f  %10f\n", ct->surface, new_surf, (adiapop_old[ct->surface] - adiapop[ct->surface])/adiapop_old[ct->surface], probtmp1[new_surf]);
       fprintf(f,"%8.2f %10d %17d %21f %19f %19f    UnavoidedCrossing \n",(step*ct->rk_timestep)*1000./PS_TO_AU, ct->surface, new_surf,
              probtmp1[new_surf], fabs(ct->tfs_overlap[ct->surface][new_surf]), (ct->ev_adiab[new_surf] - ct->ev_adiab[ct->surface])*HARTREE_TO_EV);
       newline = 1;
    }

    for (i = 0; i < ct->dim; i++)
        printf("%2d  %2d  %15.5f \n", i, ind[i], probtmp1[i]);

// stochastic algorithm
    probsum = 0.0;
    random_number = drand48();
    for (i = 0; i < ct->dim; i++)
    {
        ct->surf_prob[i] = probtmp1[i];

        if (ct->surf_prob[i] >= 0.0 && i != ct->surface)
        {
           probsum += ct->surf_prob[i];

           if (random_number < probsum)
           {
              hopsurf = i;
// calculate the nonadibatic coupling vector for velocity adjustment
              dftb->do_ncv = 1;

              if (dftb->do_ncv && hopsurf != new_surf)
              {
                 ct->ncv_isurf = ct->surface; ct->ncv_jsurf = hopsurf;

#if GMX_MPI
                 get_MM_params(ct, dftb, ct_mpi_comm, ct_mpi_rank, ct_mpi_size);
#else
                 get_MM_params(ct, dftb);
#endif

                 aa = 0.0; bb = 0.0;
                 cc = ct->ev_adiab[hopsurf] - ct->ev_adiab[ct->surface];

                 for (iatom = 0; iatom < mdatoms->homenr; iatom++)
                 {
                     n = 0;
                     for (k = 0; k < ct->sites; k++)
                     {
                         for (l = 0; l < ct->site[k].atoms; l++)
                         {
                             if (iatom == ct->site[k].atom[l])
                             {
                                for (m = 0; m < DIM; m++)
                                {
                                    ct->tfs_ncv[k][l][m] = (real) (dftb->phase1[k].grad[l][m] + dftb->phase2.grad[n][m])/cc;

                                    aa += SQR(ct->tfs_ncv[k][l][m])/2./dftb->phase1[k].mass[l]/AMU_TO_AU;
                                    bb += ct->tfs_ncv[k][l][m]*state_global->v[iatom][m]*NM_TO_BOHR/PS_TO_AU;

                                }
                             }
                             n++;
                         }
                     }
                 }

// adjust velocity along the direction of nonadiabtic coupling vector
                 delta = SQR(bb) - 4.*aa*cc;
                 if (delta >= 0.0)
                 {
                    kai1 = (- bb + sqrt(delta))/2./aa;
                    kai2 = (- bb - sqrt(delta))/2./aa;

                    kai = (fabs(kai1) < fabs(kai2) ? kai1 : kai2);

                    fprintf(f,"Hop from %d to %d succeeds! Adjust velocity along NCV with kai = %8.5f \n", ct->surface, hopsurf, kai);
                 }
                 else
                 {
                    kai = - bb/aa;
                    hopsurf = ct->surface;

                    fprintf(f,"Hop from %d to %d fails! Reverse velocity along NCV with kai = %8.5f \n", ct->surface, i, kai);
                 }

// velocity adjustment
                 for (iatom = 0; iatom < mdatoms->homenr; iatom++)
                 {
                     for (k = 0; k < ct->sites; k++)
                     {
                         for (l = 0; l < ct->site[k].atoms; l++)
                         {
                             if (iatom == ct->site[k].atom[l])
                             {
                                for (m = 0; m < DIM; m++)
                                    state_global->v[iatom][m] += kai*ct->tfs_ncv[k][l][m]/dftb->phase1[k].mass[l]/AMU_TO_AU/NM_TO_BOHR*PS_TO_AU;
                             }
                         }
                     }
                 }

                 dftb->do_ncv = 0;
              }

              if (hopsurf != ct->surface)
              {
                 fprintf(f,"%8.2f %10d %17d %21f %19f %19f \n",(step*ct->rk_timestep)*1000./PS_TO_AU, ct->surface, hopsurf,
                        ct->surf_prob[i], fabs(ct->tfs_overlap[ct->surface][hopsurf]), (ct->ev_adiab[hopsurf] - ct->ev_adiab[ct->surface])*HARTREE_TO_EV);
                 newline = 1;
              }

              ct->surface = hopsurf;

              break;
           }
        }
    }

   if (newline)
      fprintf(f, "\n");

    printf("Current Surface: %d\n", ct->surface);

    ekintot = 0.0;
    for (iatom = 0; iatom < mdatoms->homenr; iatom++)
    {
        for (k = 0; k < ct->sites; k++)
        {
            for (l = 0; l < ct->site[k].atoms; l++)
            {
                if (iatom == ct->site[k].atom[l])
                {
                   for (m = 0; m < DIM; m++)
                       ekintot += 0.5*dftb->phase1[k].mass[l]*AMU_TO_AU*SQR(state_global->v[iatom][m]*NM_TO_BOHR/PS_TO_AU);
                }
            }
        }
    }
    printf("Kinetic energy is %f\n", ekintot);

// restricted decoherence correction
    dot_product = 0.0;
    if (adiapop[ct->surface] > 0.001)
    {
       for (i = 0; i < ct->dim; i++)
       {
           if (i != ct->surface)
           {
               decay_time = (1.0 + 0.1/ekintot)/fabs(ct->ev_adiab[i] - ct->ev_adiab[ct->surface]);

               ct->tfs_popul[i]         *= exp(-ct->rk_timestep/decay_time);
               ct->tfs_popul[i+ct->dim] *= exp(-ct->rk_timestep/decay_time);
               dot_product              += SQR(ct->tfs_popul[i])+SQR(ct->tfs_popul[i+ct->dim]);
           }
       }
       curr_surf_fac = sqrt((1.0 - dot_product)/(SQR(ct->tfs_popul[ct->surface]) + SQR(ct->tfs_popul[ct->dim+ct->surface])));
       ct->tfs_popul[ct->surface]         *= curr_surf_fac;
       ct->tfs_popul[ct->surface+ct->dim] *= curr_surf_fac;

// recalc. diab. from adiab after decoherence corr.
       for (i = 0; i < ct->dim; i++)
       {
           ct->tfs_diab[i] = 0.0; ct->tfs_diab[i+ct->dim] = 0.0;

           for (j = 0; j < ct->dim; j++)
           {
               ct->tfs_diab[i]         += ct->tfs_popul[j]*ct->tfs_vector[j][i];
               ct->tfs_diab[i+ct->dim] += ct->tfs_popul[j+ct->dim]*ct->tfs_vector[j][i];
           }
       }
    }

// charge wave function
    for (i = 0; i < ct->dim; i++)
        ct->wf[i] = ct->tfs_vector[ct->surface][i];

//store hamiltonian as old_ham
    for (i = 0; i < SQR(ct->dim); i++)
        ct->tfl_old_ham[i] = *(ct->hamiltonian[0]+i);

    return 0;
}

#if GMX_MPI
int do_crossing_corrected_fssh(int step, charge_transfer_t *ct, dftb_t *dftb, t_state *state_global, t_mdatoms *mdatoms, FILE *f,
                               MPI_Comm ct_mpi_comm, int ct_mpi_rank, int ct_mpi_size)
#else
int do_crossing_corrected_fssh(int step, charge_transfer_t *ct, dftb_t *dftb, t_state *state_global, t_mdatoms *mdatoms, FILE *f)
#endif
{
// Declaration:
//
// 1. self-consitent algorithm for the weak avoided crossings is used;
// 2. different crossing types are detected and treated differently;
//
//                                          -----------Weiwei Xie on March. 2019
// Ref: J. Phys. Chem. Lett. 915, 4319 (2018)


    int    i, j, k, l, m, n, tmp, prev_surf, old_surf, hopsurf, new_surf, surftmp, iatom, istep, nstep, *ind, indmax, *indtmp, *indtmp1, *map, newline;
    double *ham, *old_ham, **ham_full, *new_ham, *adiapop_old, *adiapop, *initprob, *probtmp1, *probtmp2, *diapop, *diatmp;
    double tdse_timestep = 1.e-5*PS_TO_AU;
    double rk_timestep0, probsum, random_number, cur_weight;
    double popul_norm, overlapmax, overlaptmp1, overlaptmp2, adiapoptmp0, probtmp0, probtmp, diapoptmp;
    double nc, kai1, kai2, kai, delta, aa, bb, cc, dd, ee;
    double decay_time, curr_surf_fac, trivial_cross_type;
    double dot_product, ekintot;
    double delta_min;


    old_ham  = ct->tfl_old_ham_adiab;
    new_ham  = ct->hamiltonian_adiab;
    ham      = ct->tfl_mean_ham;
    ham_full = ct->tfl_mean_ham_full;

    snew(diapop, ct->dim);
    snew(diatmp, 2*ct->dim);
    snew(adiapop_old, ct->dim);
    snew(adiapop, ct->dim);
    snew(initprob, ct->dim);
    snew(probtmp1, ct->dim);
    snew(probtmp2, ct->dim);

    snew(ind, ct->dim);
    snew(indtmp, ct->dim);
    snew(indtmp1, ct->dim);
    snew(map, ct->dim);

// set up initial surface
    if (ct->tfs_initialization_step)
    {
       printf("Input wavefunction is:\n");
       for (j = 0; j < ct->dim; j++)
           printf("%f ", ct->wf[j]);
       printf("\n");

       for (i = 0; i < SQR(ct->dim); i++)
       {
           ct->tfl_old_ham[i] = *(ct->hamiltonian[0]+i);
           ham[i] = ct->tfl_old_ham[i];
       }

// diagonalize this Hamiltonian
       dsyev(ct->dim, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);

// record eigenvector
       for (i = 0; i < ct->dim; i++)
           for (j = 0; j < ct->dim; j++)
               ct->tfs_vector[i][j] = ham[i * ct->dim + j];

// determine inital surface using eigenvector as a criterion
       probsum = 0.;
       random_number = drand48();
       for (i = 0; i < ct->dim; i++)
       {
           initprob[i] = 0.;
           for (j = 0; j < ct->dim; j++)
               initprob[i] += ct->tfs_vector[i][j]*ct->wf[j];
           initprob[i] = SQR(initprob[i]);

           probsum += initprob[i];
           if (random_number < probsum)
           {
               ct->surface = i;
               break;
           }
       }

       for (i = 0; i < 2*ct->dim; i++)
       {
           ct->tfs_popul[i] = 0.0; ct->tfs_diab[i]  = 0.0;
       }

       for (i = 0; i < ct->dim; i++)
           ct->tfs_diab[i] = ct->wf[i];

       ct->tfs_popul[ct->surface] = 1.0;

       ct->tfs_initialization_step = 0;

       return 0;
    }

// set up the old and new Hamiltonian for linear interpolation
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            old_ham[i*ct->dim+j] = ct->tfl_old_ham[i*ct->dim+j];
            new_ham[i*ct->dim+j] = ct->hamiltonian[i][j];
        }
    }

// store the vector at preceding time step
/*    for (i = 0; i < ct->dim; i++)
        for (j = 0; j < ct->dim; j++)
            ct->tfs_vector_old2[i][j] = ct->tfs_vector[i][j];
*/
// solve time-dependent Schroedinger equation
    nstep = (int)(ct->rk_timestep/tdse_timestep);
    rk_timestep0 = ct->rk_timestep;
    ct->rk_timestep = tdse_timestep;
    printf("Time step for TDSE: %f fs \n", ct->rk_timestep/PS_TO_AU*1000);

    for (istep = 1; istep < nstep+1; istep++)
    {
// linear interpolation for Hamiltonian
        cur_weight = (double)istep/(double)nstep;
        for (i = 0; i < ct->dim; i++)
            for (j = 0; j < ct->dim; j++)
                ham_full[i][j] = (1.0 - cur_weight)*ct->tfl_old_ham[i*ct->dim+j] + cur_weight*ct->hamiltonian[i][j];
// propagate the TDSE in the diabatic representation
        do_rksuite_diab(ct);
    }
    ct->rk_timestep = rk_timestep0;

// diagonalize this Hamiltonian
    for (i = 0; i < SQR(ct->dim); i++)
        ham[i] = new_ham[i];

    dsyev(ct->dim, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);

// store the eigenvectors at time t-dt and time t
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_vector_old[i][j] = ct->tfs_vector[i][j];
            ct->tfs_vector[i][j]     = ham[i*ct->dim + j];
        }
    }

// adiabatic wave function overlaps between time t-dt and t
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_overlap[i][j] = 0.0;

            for (k = 0; k < ct->dim; k++)
                ct->tfs_overlap[i][j] += ct->tfs_vector_old[i][k]*ct->tfs_vector[j][k];
        }
    }

// check signs
    for (i = 0; i < ct->dim; i++)
    {
        if (ct->tfs_overlap[i][i] < 0.0)
        {
            for (k = 0; k < ct->dim; k++)
            {
                ct->tfs_vector[i][k]  *= -1.0;
                ct->tfs_overlap[k][i] *= -1.0;
            }
        }
    }

// reassignment of states according to the min-cost algorithm
    assignprob(ct, ind);

// reassignment of states according to the maximum wave function overlaps
/*    for (i = 0; i < ct->dim; i++)
    {
        overlapmax = 0.0;
        for (j = 0; j < ct->dim; j++)
        {
            if (overlapmax < fabs(ct->tfs_overlap[j][i]))
            {
               overlapmax = fabs(ct->tfs_overlap[j][i]);
               ind[i] = j;
            }
        }
    }
*/
// calculate the adiabatic populations at the previous time step
    for (i = 0; i < ct->dim; i++)
        adiapop_old[i] = SQR(ct->tfs_popul[i]) + SQR(ct->tfs_popul[i+ct->dim]);

// hopping probabilities based on Tully's fewest-switches anzatz
    for (i = 0; i < ct->dim; i++)
    {
//     nonadibatic couplings based on Hammes-Schiffer and Tully, JCP 101, 4657 (1994)
//        nc = (ct->tfs_overlap[prev_surf][i] - ct->tfs_overlap[j][ct->surface])/2.0/ct->rk_timestep;
//     ndibatic couplings based norm-preserving interpolation by Benjamin, JPCL 5, 2351 (2014)
        overlaptmp1 = sqrt(fabs(1.0 - SQR(ct->tfs_overlap[i][i]) - SQR(ct->tfs_overlap[ct->surface][i])));
        overlaptmp2 = (- ct->tfs_overlap[i][ct->surface]*ct->tfs_overlap[i][i] - ct->tfs_overlap[ct->surface][ct->surface]*ct->tfs_overlap[ct->surface][i])/overlaptmp1;

        aa = - sin(acos(ct->tfs_overlap[i][i]) - asin(ct->tfs_overlap[i][ct->surface]))/(acos(ct->tfs_overlap[i][i]) - asin(ct->tfs_overlap[i][ct->surface]));
        bb =   sin(acos(ct->tfs_overlap[i][i]) + asin(ct->tfs_overlap[i][ct->surface]))/(acos(ct->tfs_overlap[i][i]) + asin(ct->tfs_overlap[i][ct->surface]));
        cc =   sin(acos(ct->tfs_overlap[ct->surface][ct->surface]) - asin(ct->tfs_overlap[ct->surface][i]))/(acos(ct->tfs_overlap[ct->surface][ct->surface]) - asin(ct->tfs_overlap[ct->surface][i]));
        dd =   sin(acos(ct->tfs_overlap[ct->surface][ct->surface]) + asin(ct->tfs_overlap[ct->surface][i]))/(acos(ct->tfs_overlap[ct->surface][ct->surface]) + asin(ct->tfs_overlap[ct->surface][i]));
        ee = 2.0*asin(overlaptmp1*(overlaptmp1*overlaptmp2*asin(overlaptmp1) + (sqrt((1.0 - SQR(overlaptmp1))*(1.0 - SQR(overlaptmp2))) - 1.0)*asin(overlaptmp2)))/(SQR(asin(overlaptmp1)) - SQR(asin(overlaptmp2)));

        nc = (acos(ct->tfs_overlap[i][i])*(aa + bb) + asin(ct->tfs_overlap[ct->surface][i])*(cc + dd) + ee)/2.0/ct->rk_timestep;

        probtmp1[i] = (2.0*(ct->tfs_popul[ct->surface]*ct->tfs_popul[i]+ct->tfs_popul[ct->dim+ct->surface]*ct->tfs_popul[ct->dim+i])*nc)*ct->rk_timestep/adiapop_old[ct->surface];

        if (probtmp1[i] < 0.0 || i == ct->surface)
           probtmp1[i] = 0.0;
    }

// calculate expansion coefficients in the adiabatic basis
    for (i = 0; i < ct->dim; i++)
    {
        ct->tfs_popul[i] = 0.0; ct->tfs_popul[i+ct->dim] = 0.0;

        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_popul[i]         += ct->tfs_vector[i][j]*ct->tfs_diab[j];
            ct->tfs_popul[i+ct->dim] += ct->tfs_vector[i][j]*ct->tfs_diab[j+ct->dim];
        }
    }

// record the adiabatic populations at new time
    for (i = 0; i < ct->dim; i++)
        adiapop[i] = SQR(ct->tfs_popul[i])+SQR(ct->tfs_popul[i+ct->dim]);

// SC-FSSH, avoid the numerical instability of the NCs
    prev_surf = ct->surface;
    for (i = 0; i < ct->dim; i++)
    {
        if (prev_surf == ind[i])
        {
           new_surf = i;
           break;
        }
    }

    newline = 0;
    if (ct->surface != new_surf)
    {
       probtmp1[new_surf] = (adiapop_old[ct->surface] - adiapop[ct->surface])/adiapop_old[ct->surface];

       for (i = 0; i < ct->dim; i++)
       {
           if (i != new_surf)
           {
              probtmp1[new_surf] -= probtmp1[i];
           }
       }

//       printf("SC is actived: %d  %d  %10f  %10f\n", ct->surface, new_surf, (adiapop_old[ct->surface] - adiapop[ct->surface])/adiapop_old[ct->surface], probtmp1[new_surf]);
       fprintf(f,"%8.2f %10d %17d %21f %19f %19f   UnavoidedCrossing \n",(step*ct->rk_timestep)*1000./PS_TO_AU, ct->surface, new_surf,
              probtmp1[new_surf], fabs(ct->tfs_overlap[ct->surface][new_surf]), (ct->ev_adiab[new_surf] - ct->ev_adiab[ct->surface])*HARTREE_TO_EV);
       newline = 1;
    }

    for (i = 0; i < ct->dim; i++)
        printf("%2d  %2d  %15.5f \n", i, ind[i], probtmp1[i]);

// stochastic algorithm
    probsum = 0.0;
    random_number = drand48();
    for (i = 0; i < ct->dim; i++)
    {
        ct->surf_prob[i] = probtmp1[i];

        if (ct->surf_prob[i] >= 0.0 && i != ct->surface)
        {
           probsum += ct->surf_prob[i];

           if (random_number < probsum)
           {
              hopsurf = i;
// calculate the nonadibatic coupling vector for velocity adjustment
              dftb->do_ncv = 1;

              if (dftb->do_ncv)
              {
                 ct->ncv_isurf = ct->surface; ct->ncv_jsurf = hopsurf;

#if GMX_MPI
                 get_MM_params(ct, dftb, ct_mpi_comm, ct_mpi_rank, ct_mpi_size);
#else
                 get_MM_params(ct, dftb);
#endif

                 aa = 0.0; bb = 0.0;
                 cc = ct->ev_adiab[hopsurf] - ct->ev_adiab[ct->surface];

                 for (iatom = 0; iatom < mdatoms->homenr; iatom++)
                 {
                     n = 0;
                     for (k = 0; k < ct->sites; k++)
                     {
                         for (l = 0; l < ct->site[k].atoms; l++)
                         {
                             if (iatom == ct->site[k].atom[l])
                             {
                                for (m = 0; m < DIM; m++)
                                {
                                    ct->tfs_ncv[k][l][m] = (real) (dftb->phase1[k].grad[l][m] + dftb->phase2.grad[n][m])/cc;

                                    aa += SQR(ct->tfs_ncv[k][l][m])/2./dftb->phase1[k].mass[l]/AMU_TO_AU;
                                    bb += ct->tfs_ncv[k][l][m]*state_global->v[iatom][m]*NM_TO_BOHR/PS_TO_AU;

                                }
                             }
                             n++;
                         }
                     }
                 }

// crossing-corrected FSSH
                 for (j = 0; j < ct->dim; j++)
                 {
                     if (hopsurf == ind[j])
                     {
                        surftmp = j;
                        break;
                     }
                 }

                 delta = SQR(bb) - 4.0*aa*cc;

                 if (ct->surface == new_surf && delta >= 0.0)
                 {
                    kai1 = (- bb + sqrt(delta))/2./aa;
                    kai2 = (- bb - sqrt(delta))/2./aa;

                    kai = (fabs(kai1) < fabs(kai2) ? kai1 : kai2);

                    hopsurf = surftmp;

                    fprintf(f, "Hop from %d to %d succeeds! Adjust velocity along NCV with kai = %8.5f \n", ct->surface, hopsurf, kai);
                 }
                 else if (ct->surface == new_surf && delta < 0.0)
                 {
                    kai = - bb/aa;
                    hopsurf = ct->surface;

                    fprintf(f, "Hop from %d to %d fails! Reverse velocity along NCV with kai = %8.5f \n", ct->surface, i, kai);
                 }
                 else if (ct->surface != new_surf && hopsurf == new_surf)
                 {
                    kai = 0.0;
                    hopsurf = new_surf;

                    fprintf(f, "Hop from %d to %d (new) succeeds! \n", ct->surface, hopsurf);
                 }
                 else if (ct->surface != new_surf && hopsurf != new_surf && delta >= 0.0)
                 {
                    kai1 = (- bb + sqrt(delta))/2./aa;
                    kai2 = (- bb - sqrt(delta))/2./aa;

                    kai = (fabs(kai1) < fabs(kai2) ? kai1 : kai2);

                    hopsurf = surftmp;

                    fprintf(f, "Hop from %d to %d (new) succeeds! Adjust velocity along NCV with kai = %8.5f \n", ct->surface, hopsurf, kai);
                 }
                 else if (ct->surface != new_surf && hopsurf != new_surf && delta < 0.0)
                 {
                    kai = - bb/aa;
                    hopsurf = new_surf;

                    fprintf(f, "Hop from %d to %d (new) fails! Reverse velocity along NCV with kai = %8.5f \n", ct->surface, i, kai);
                 }
              }

              if (hopsurf != ct->surface)
              {
                 fprintf(f,"%8.2f %10d %17d %21f %19f %19f \n",(step*ct->rk_timestep)*1000./PS_TO_AU, ct->surface, hopsurf,
                        ct->surf_prob[i], fabs(ct->tfs_overlap[ct->surface][hopsurf]), (ct->ev_adiab[hopsurf] - ct->ev_adiab[ct->surface])*HARTREE_TO_EV);
                 newline = 1;
              }

              ct->surface = hopsurf;

// velocity adjustment
              if (dftb->do_ncv)
              {
                 if (kai != 0.0)
                 {
                    for (iatom = 0; iatom < mdatoms->homenr; iatom++)
                    {
                        for (k = 0; k < ct->sites; k++)
                        {
                            for (l = 0; l < ct->site[k].atoms; l++)
                            {
                                if (iatom == ct->site[k].atom[l])
                                {
                                   for (m = 0; m < DIM; m++)
                                       state_global->v[iatom][m] += kai*ct->tfs_ncv[k][l][m]/dftb->phase1[k].mass[l]/AMU_TO_AU/NM_TO_BOHR*PS_TO_AU;
                                }
                            }
                        }
                    }
                 }
                 dftb->do_ncv = 0;
              }
              break;
           }
        }
    }

    if (newline)
       fprintf(f,"\n");

    printf("Current Surface: %d\n", ct->surface);

    ekintot = 0.0;
    for (iatom = 0; iatom < mdatoms->homenr; iatom++)
    {
        for (k = 0; k < ct->sites; k++)
        {
            for (l = 0; l < ct->site[k].atoms; l++)
            {
                if (iatom == ct->site[k].atom[l])
                {
                   for (m = 0; m < DIM; m++)
                       ekintot += 0.5*dftb->phase1[k].mass[l]*AMU_TO_AU*SQR(state_global->v[iatom][m]*NM_TO_BOHR/PS_TO_AU);
                }
            }
        }
    }
    printf("Kinetic energy is %f\n", ekintot);

// restricted decoherence correction
//    if (adiapop[ct->surface] > 0.001)
//    {
    dot_product = 0.0;
    for (i = 0; i < ct->dim; i++)
    {
        if (i != ct->surface)
        {
            decay_time = (1.0 + 0.1/ekintot)/fabs(ct->ev_adiab[i] - ct->ev_adiab[ct->surface]);

            ct->tfs_popul[i]         *= exp(-ct->rk_timestep/decay_time);
            ct->tfs_popul[i+ct->dim] *= exp(-ct->rk_timestep/decay_time);
            dot_product              += SQR(ct->tfs_popul[i])+SQR(ct->tfs_popul[i+ct->dim]);
        }
    }
    curr_surf_fac = sqrt((1.0 - dot_product)/(SQR(ct->tfs_popul[ct->surface]) + SQR(ct->tfs_popul[ct->dim+ct->surface])));
    ct->tfs_popul[ct->surface]         *= curr_surf_fac;
    ct->tfs_popul[ct->surface+ct->dim] *= curr_surf_fac;

// recalc. diab. from adiab after decoherence corr.
    for (i = 0; i < ct->dim; i++)
    {
        ct->tfs_diab[i] = 0.0; ct->tfs_diab[i+ct->dim] = 0.0;

        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_diab[i]         += ct->tfs_popul[j]*ct->tfs_vector[j][i];
            ct->tfs_diab[i+ct->dim] += ct->tfs_popul[j+ct->dim]*ct->tfs_vector[j][i];
        }
    }
//    }

// charge wave function
    for (i = 0; i < ct->dim; i++)
        ct->wf[i] = ct->tfs_vector[ct->surface][i];

//store hamiltonian as old_ham
    for (i = 0; i < SQR(ct->dim); i++)
        ct->tfl_old_ham[i] = *(ct->hamiltonian[0]+i);

    return 0;
}

#if GMX_MPI
int do_global_flux_fssh(int step, charge_transfer_t *ct, dftb_t *dftb, t_state *state_global, t_mdatoms *mdatoms, FILE *f,
                        MPI_Comm ct_mpi_comm, int ct_mpi_rank, int ct_mpi_size)
#else
int do_global_flux_fssh(int step, charge_transfer_t *ct, dftb_t *dftb, t_state *state_global, t_mdatoms *mdatoms, FILE *f)
#endif
{
// Declaration:
//
// hopping probabilities are only determined by populations;
// new: add the threshold for the detection of trivial crossings
//
//                                          -----------Weiwei Xie on July. 2019
// Ref: J. Chem. Theory Comput. 10, 9, 3598 (2014)


    int    i, j, k, l, m, n, tmp, prev_surf, old_surf, hopsurf, new_surf, surftmp, iatom, istep, nstep, *ind, newline, *indtmp, *map;
    double *ham, *old_ham, **ham_full, *new_ham, *adiapop_old, *adiapop, *dadiapop,  *initprob, *diapop;
    double tdse_timestep = 1.e-5*PS_TO_AU, threshold = 0.99;
    double rk_timestep0, probsum, random_number, random_number1, cur_weight, dadiapopsum, overlapmax, probmax;
    double nc, kai1, kai2, kai, delta, aa, bb, cc, dd, ee;
    double decay_time, curr_surf_fac;
    double dot_product, ekintot;


    old_ham  = ct->tfl_old_ham_adiab;
    new_ham  = ct->hamiltonian_adiab;
    ham      = ct->tfl_mean_ham;
    ham_full = ct->tfl_mean_ham_full;

    snew(diapop, ct->dim);
    snew(adiapop_old, ct->dim);
    snew(adiapop, ct->dim);
    snew(dadiapop, ct->dim);
    snew(initprob, ct->dim);
    snew(ind, ct->dim);
    snew(map, ct->dim);
    snew(indtmp, ct->dim);

// set up initial surface
    if (ct->tfs_initialization_step)
    {
       printf("Input wavefunction is:\n");
       for (j = 0; j < ct->dim; j++)
           printf("%f ", ct->wf[j]);
       printf("\n");

       for (i = 0; i < SQR(ct->dim); i++)
       {
           ct->tfl_old_ham[i] = *(ct->hamiltonian[0]+i);
           ham[i] = ct->tfl_old_ham[i];
       }

// diagonalize this Hamiltonian
       dsyev(ct->dim, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);

// record eigenvector
       for (i = 0; i < ct->dim; i++)
           for (j = 0; j < ct->dim; j++)
               ct->tfs_vector[i][j] = ham[i * ct->dim + j];

// determine inital surface using eigenvector as a criterion
       probsum = 0.;
       random_number = drand48();
       for (i = 0; i < ct->dim; i++)
       {
           initprob[i] = 0.;
           for (j = 0; j < ct->dim; j++)
               initprob[i] += ct->tfs_vector[i][j]*ct->wf[j];
           initprob[i] = SQR(initprob[i]);

           probsum += initprob[i];
           if (random_number < probsum)
           {
               ct->surface = i;
               break;
           }
       }

       for (i = 0; i < 2*ct->dim; i++)
       {
           ct->tfs_popul[i] = 0.0; ct->tfs_diab[i]  = 0.0;
       }

       for (i = 0; i < ct->dim; i++)
           ct->tfs_diab[i] = ct->wf[i];

       ct->tfs_popul[ct->surface] = 1.0;

       ct->tfs_initialization_step = 0;

       return 0;
    }

// set up the old and new Hamiltonian for linear interpolation
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            old_ham[i*ct->dim+j] = ct->tfl_old_ham[i*ct->dim+j];
            new_ham[i*ct->dim+j] = ct->hamiltonian[i][j];
        }
    }

// solve time-dependent Schroedinger equation
    nstep = (int)(ct->rk_timestep/tdse_timestep);
    rk_timestep0 = ct->rk_timestep;
    ct->rk_timestep = tdse_timestep;
    printf("Time step for TDSE: %f fs \n", ct->rk_timestep/PS_TO_AU*1000);

    for (istep = 1; istep < nstep+1; istep++)
    {
// linear interpolation for Hamiltonian
        cur_weight = (double)istep/(double)nstep;
        for (i = 0; i < ct->dim; i++)
            for (j = 0; j < ct->dim; j++)
                ham_full[i][j] = (1.0 - cur_weight)*ct->tfl_old_ham[i*ct->dim+j] + cur_weight*ct->hamiltonian[i][j];
// propagate the TDSE in the diabatic representation
        do_rksuite_diab(ct);
    }
    ct->rk_timestep = rk_timestep0;

// diagonalize this Hamiltonian
    for (i = 0; i < SQR(ct->dim); i++)
        ham[i] = new_ham[i];

    dsyev(ct->dim, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);

// store the eigenvectors at time t-dt and time t
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_vector_old[i][j] = ct->tfs_vector[i][j];
            ct->tfs_vector[i][j]     = ham[i*ct->dim + j];
        }
    }

// adiabatic wave function overlaps between time t-dt and t
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_overlap[i][j] = 0.0;

            for (k = 0; k < ct->dim; k++)
                ct->tfs_overlap[i][j] += ct->tfs_vector_old[i][k]*ct->tfs_vector[j][k];
        }
    }

// reassignment of states according to the min-cost algorithm
    assignprob(ct, indtmp);

// state tracking
    for (i = 0; i < ct->dim; i++)
    {
        map[i] = 0; ind[i] = -1;
    }

    for (i = 0; i < ct->dim; i++)
    {
        if (fabs(ct->tfs_overlap[indtmp[i]][i]) >= 0.9)
        {
           ind[i] = indtmp[i];

           for (j = 0; j < ct->dim; j++)
           {
               if (ind[i] == j)
               {
                  map[j] = 1;
                  break;
               }
           }
        }
    }

    for (i = 0; i < ct->dim; i++)
    {
        if (ind[i] == -1)
        {
           for (j = 0; j < ct->dim; j++)
           {
               if (! map[j])
               {
                  ind[i] = j; map[j] = 1;
                  break;
               }
           }
        }
    }

// reassign the active state
    prev_surf = ct->surface;

    if (ct->surface != ind[ct->surface])
    {
       for (i = 0; i < ct->dim; i++)
       {
           if (ct->surface == ind[i])
           {
              printf("Changed from surface %d to surface %d \n", ct->surface, i);
              ct->surface = i;
              break;
           }
       }
    }

    for (i = 0; i < ct->dim; i++)
    {
        if (prev_surf == indtmp[i])
        {
           new_surf = i;
           break;
        }
    }

// check signs
    for (i = 0; i < ct->dim; i++)
    {
        j = ind[i];

        if (ct->tfs_overlap[j][i] < 0.0)
        {
            for (k = 0; k < ct->dim; k++)
            {
                ct->tfs_vector[i][k]  *= -1.0;
                ct->tfs_overlap[k][i] *= -1.0;
            }
        }
    }

// calculate the adiabatic populations at the previous time step
    for (i = 0; i < ct->dim; i++)
        adiapop_old[i] = SQR(ct->tfs_popul[i]) + SQR(ct->tfs_popul[i+ct->dim]);

// calculate expansion coefficients in the adiabatic basis
    for (i = 0; i < ct->dim; i++)
    {
        ct->tfs_popul[i] = 0.0; ct->tfs_popul[i+ct->dim] = 0.0;

        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_popul[i]         += ct->tfs_vector[i][j]*ct->tfs_diab[j];
            ct->tfs_popul[i+ct->dim] += ct->tfs_vector[i][j]*ct->tfs_diab[j+ct->dim];
        }
    }

// record the adiabatic populations at new time
    dadiapopsum = 0.0;
    for (i = 0; i < ct->dim; i++)
    {
        j = ind[i];

        adiapop[i] = SQR(ct->tfs_popul[i])+SQR(ct->tfs_popul[i+ct->dim]);

        dadiapop[i] = adiapop[i] - adiapop_old[j];

        if (dadiapop[i] < 0.0)
           dadiapopsum += dadiapop[i];
    }

// hopping probability
    probmax = 0.0;
    for (i = 0; i < ct->dim; i++)
    {
        ct->surf_prob[i] = 0.0;

        if (dadiapop[ct->surface] < 0.0 && dadiapop[i] > 0.0)
           ct->surf_prob[i] = dadiapop[i]/adiapop_old[ind[ct->surface]]*dadiapop[ct->surface]/dadiapopsum;
    }

    for (i = 0; i < ct->dim; i++)
        printf("%2d  %15.5f  %15.5f \n", i, dadiapop[i]/adiapop_old[ct->surface], ct->surf_prob[i]);

    newline = 0;
    if (ct->surface != new_surf)
    {
       fprintf(f,"%8.2f %10d %17d %21f %19f %19f   UnavoidedCrossing \n",(step*ct->rk_timestep)*1000./PS_TO_AU, ct->surface, new_surf,
              ct->surf_prob[new_surf], fabs(ct->tfs_overlap[ct->surface][new_surf]), (ct->ev_adiab[new_surf] - ct->ev_adiab[ct->surface])*HARTREE_TO_EV);
       newline = 1;
    }

// stochastic algorithm
    probsum = 0.0;
    random_number = drand48();
    for (i = 0; i < ct->dim; i++)
    {
        if (ct->surf_prob[i] >= 0.0 && i != ct->surface)
        {
           probsum += ct->surf_prob[i];

           if (random_number < probsum)
           {
              hopsurf = i;

              dftb->do_ncv = 1;

              if (dftb->do_ncv == 1)
              {
                 ct->ncv_isurf = ct->surface; ct->ncv_jsurf = hopsurf;

#if GMX_MPI
                 get_MM_params(ct, dftb, ct_mpi_comm, ct_mpi_rank, ct_mpi_size);
#else
                 get_MM_params(ct, dftb);
#endif

                 aa = 0.0; bb = 0.0;
                 cc = ct->ev_adiab[hopsurf] - ct->ev_adiab[ct->surface];

                 for (iatom = 0; iatom < mdatoms->homenr; iatom++)
                 {
                     n = 0;
                     for (k = 0; k < ct->sites; k++)
                     {
                         for (l = 0; l < ct->site[k].atoms; l++)
                         {
                             if (iatom == ct->site[k].atom[l])
                             {
                                for (m = 0; m < DIM; m++)
                                {
                                    ct->tfs_ncv[k][l][m] = (real) (dftb->phase1[k].grad[l][m] + dftb->phase2.grad[n][m])/cc;

                                    aa += SQR(ct->tfs_ncv[k][l][m])/2./dftb->phase1[k].mass[l]/AMU_TO_AU;
                                    bb += ct->tfs_ncv[k][l][m]*state_global->v[iatom][m]*NM_TO_BOHR/PS_TO_AU;

                                }
                             }
                             n++;
                         }
                     }
                 }

// adjust velocity along the direction of nonadiabtic coupling vector
                 delta = SQR(bb) - 4.*aa*cc;
                 if (delta >= 0.0)
                 {
                    kai1 = (- bb + sqrt(delta))/2./aa;
                    kai2 = (- bb - sqrt(delta))/2./aa;

                    kai = (fabs(kai1) < fabs(kai2) ? kai1 : kai2);

                    fprintf(f, "Hop from %d to %d succeeds! Adjust velocity along NCV with kai = %8.5f \n", ct->surface, hopsurf, kai);
                 }
                 else if (delta < 0.0)
                 {
                    kai = - bb/aa;
                    hopsurf = ct->surface;

                    fprintf(f, "Hop from %d to %d fails! Reverse velocity along NCV with kai = %8.5f \n", ct->surface, i, kai);
                 }

// velocity adjustment
                 for (iatom = 0; iatom < mdatoms->homenr; iatom++)
                 {
                     for (k = 0; k < ct->sites; k++)
                     {
                         for (l = 0; l < ct->site[k].atoms; l++)
                         {
                             if (iatom == ct->site[k].atom[l])
                             {
                                for (m = 0; m < DIM; m++)
                                    state_global->v[iatom][m] += kai*ct->tfs_ncv[k][l][m]/dftb->phase1[k].mass[l]/AMU_TO_AU/NM_TO_BOHR*PS_TO_AU;
                             }
                         }
                     }
                 }
                 dftb->do_ncv = 0;
              }
              else if (dftb->do_ncv == 2)
              {
                 random_number1 = drand48();

                 if (random_number1 > exp(-(ct->ev_adiab[hopsurf] - ct->ev_adiab[ct->surface])/ct->fermi_kt))
                 {
                      hopsurf = ct->surface;
                      fprintf(f,"Hop from %d to %d fails! \n", ct->surface, i);
                 }
                 else
                 {
                    fprintf(f,"Hop from %d to %d succeeds! \n", ct->surface, hopsurf);
                 }
                 dftb->do_ncv = 0;
              }

              if (hopsurf != ct->surface)
              {
                 fprintf(f,"%8.2f %10d %17d %21f %19f %19f \n \n",(step*ct->rk_timestep)*1000./PS_TO_AU, ct->surface, hopsurf,
                        ct->surf_prob[i], fabs(ct->tfs_overlap[ct->surface][hopsurf]), (ct->ev_adiab[hopsurf] - ct->ev_adiab[ct->surface])*HARTREE_TO_EV);
                 newline = 1;
              }

              ct->surface = hopsurf;

              break;
           }
        }
    }

    if (newline)
       fprintf(f,"\n");

    printf("Current Surface: %d\n", ct->surface);

    ekintot = 0.0;
    for (iatom = 0; iatom < mdatoms->homenr; iatom++)
    {
        for (k = 0; k < ct->sites; k++)
        {
            for (l = 0; l < ct->site[k].atoms; l++)
            {
                if (iatom == ct->site[k].atom[l])
                {
                   for (m = 0; m < DIM; m++)
                       ekintot += 0.5*dftb->phase1[k].mass[l]*AMU_TO_AU*SQR(state_global->v[iatom][m]*NM_TO_BOHR/PS_TO_AU);
                }
            }
        }
    }
    printf("Kinetic energy is %f\n", ekintot);

// restricted decoherence correction
    dot_product = 0.0;
    if (adiapop[ct->surface] > 0.001)
    {
    for (i = 0; i < ct->dim; i++)
    {
        if (i != ct->surface)
        {
            decay_time = (1.0 + 0.1/ekintot)/fabs(ct->ev_adiab[i] - ct->ev_adiab[ct->surface]);

            ct->tfs_popul[i]         *= exp(-ct->rk_timestep/decay_time);
            ct->tfs_popul[i+ct->dim] *= exp(-ct->rk_timestep/decay_time);
            dot_product              += SQR(ct->tfs_popul[i])+SQR(ct->tfs_popul[i+ct->dim]);
        }
    }
    curr_surf_fac = sqrt((1.0 - dot_product)/(SQR(ct->tfs_popul[ct->surface]) + SQR(ct->tfs_popul[ct->dim+ct->surface])));
    ct->tfs_popul[ct->surface]         *= curr_surf_fac;
    ct->tfs_popul[ct->surface+ct->dim] *= curr_surf_fac;

// recalc. diab. from adiab after decoherence corr.
    for (i = 0; i < ct->dim; i++)
    {
        ct->tfs_diab[i] = 0.0; ct->tfs_diab[i+ct->dim] = 0.0;

        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_diab[i]         += ct->tfs_popul[j]*ct->tfs_vector[j][i];
            ct->tfs_diab[i+ct->dim] += ct->tfs_popul[j+ct->dim]*ct->tfs_vector[j][i];
        }
    }
    }

// charge wave function
    for (i = 0; i < ct->dim; i++)
        ct->wf[i] = ct->tfs_vector[ct->surface][i];

//store hamiltonian as old_ham
    for (i = 0; i < SQR(ct->dim); i++)
        ct->tfl_old_ham[i] = *(ct->hamiltonian[0]+i);

    return 0;
}


/*int do_global_flux_fssh(int step, charge_transfer_t *ct, dftb_t *dftb, t_state *state_global, t_mdatoms *mdatoms, FILE *f)
{
// Declaration:
//
// hopping probabilities are only determined by populations;
// new: add the threshold for the detection of trivial crossings
//
//                                          -----------Weiwei Xie on July. 2019
// Ref: J. Chem. Theory Comput. 10, 9, 3598 (2014)


    int    i, j, k, l, m, n, tmp, prev_surf, old_surf, hopsurf, new_surf, surftmp, iatom, istep, nstep, *ind;
    double *ham, *old_ham, **ham_full, *new_ham, *adiapop_old, *adiapop, *dadiapop,  *initprob, *diapop;
    double tdse_timestep = 1.e-6*PS_TO_AU, threshold = 0.99;
    double rk_timestep0, probsum, random_number, random_number1, cur_weight, dadiapopsum, overlapmax, probmax;
    double nc, kai1, kai2, kai, delta, aa, bb, cc, dd, ee;
    double decay_time, curr_surf_fac;
    double dot_product, ekintot;


    old_ham  = ct->tfl_old_ham_adiab;
    new_ham  = ct->hamiltonian_adiab;
    ham      = ct->tfl_mean_ham;
    ham_full = ct->tfl_mean_ham_full;

    snew(diapop, ct->dim);
    snew(adiapop_old, ct->dim);
    snew(adiapop, ct->dim);
    snew(dadiapop, ct->dim);
    snew(initprob, ct->dim);
    snew(ind, ct->dim);

// set up initial surface
    if (ct->tfs_initialization_step)
    {
       printf("Input wavefunction is:\n");
       for (j = 0; j < ct->dim; j++)
           printf("%f ", ct->wf[j]);
       printf("\n");

       for (i = 0; i < SQR(ct->dim); i++)
       {
           ct->tfl_old_ham[i] = *(ct->hamiltonian[0]+i);
           ham[i] = ct->tfl_old_ham[i];
       }

// diagonalize this Hamiltonian
       dsyev(ct->dim, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);

// record eigenvector
       for (i = 0; i < ct->dim; i++)
           for (j = 0; j < ct->dim; j++)
               ct->tfs_vector[i][j] = ham[i * ct->dim + j];

// determine inital surface using eigenvector as a criterion
       probsum = 0.;
       random_number = drand48();
       for (i = 0; i < ct->dim; i++)
       {
           initprob[i] = 0.;
           for (j = 0; j < ct->dim; j++)
               initprob[i] += ct->tfs_vector[i][j]*ct->wf[j];
           initprob[i] = SQR(initprob[i]);

           probsum += initprob[i];
           if (random_number < probsum)
           {
               ct->surface = i;
               break;
           }
       }

       for (i = 0; i < 2*ct->dim; i++)
       {
           ct->tfs_popul[i] = 0.0; ct->tfs_diab[i]  = 0.0;
       }

       for (i = 0; i < ct->dim; i++)
           ct->tfs_diab[i] = ct->wf[i];

       ct->tfs_popul[ct->surface] = 1.0;

       ct->tfs_initialization_step = 0;

       return 0;
    }

// set up the old and new Hamiltonian for linear interpolation
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            old_ham[i*ct->dim+j] = ct->tfl_old_ham[i*ct->dim+j];
            new_ham[i*ct->dim+j] = ct->hamiltonian[i][j];
        }
    }


    ekintot = 0.0;
    for (iatom = 0; iatom < mdatoms->homenr; iatom++)
    {
        for (k = 0; k < ct->sites; k++)
        {
            for (l = 0; l < ct->site[k].atoms; l++)
            {
                if (iatom == ct->site[k].atom[l])
                {
                   for (m = 0; m < DIM; m++)
                       ekintot += 0.5*dftb->phase1[k].mass[l]*AMU_TO_AU*SQR(state_global->v[iatom][m]*NM_TO_BOHR/PS_TO_AU);
                }
            }
        }
    }
    printf("Kinetic energy is %f\n", ekintot);

// solve time-dependent Schroedinger equation
    nstep = (int)(ct->rk_timestep/tdse_timestep);
    rk_timestep0 = ct->rk_timestep;
    ct->rk_timestep = tdse_timestep;
    printf("Time step for TDSE: %f fs \n", ct->rk_timestep/PS_TO_AU*1000);

    for (istep = 1; istep < nstep+1; istep++)
    {
// linear interpolation for Hamiltonian
        cur_weight = (double)istep/(double)nstep;
        for (i = 0; i < ct->dim; i++)
            for (j = 0; j < ct->dim; j++)
                ham_full[i][j] = (1.0 - cur_weight)*ct->tfl_old_ham[i*ct->dim+j] + cur_weight*ct->hamiltonian[i][j];
// propagate the TDSE in the diabatic representation
        do_rksuite_diab(ct);

// diagonalize this Hamiltonian
    for (i = 0; i < SQR(ct->dim); i++)
        ham[i] = (1.0-cur_weight)*old_ham[i]+cur_weight*new_ham[i];

    dsyev(ct->dim, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);

// store the eigenvectors at time t-dt and time t
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_vector_old[i][j] = ct->tfs_vector[i][j];
            ct->tfs_vector[i][j]     = ham[i*ct->dim + j];
        }
    }

// adiabatic wave function overlaps between time t-dt and t
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_overlap[i][j] = 0.0;

            for (k = 0; k < ct->dim; k++)
                ct->tfs_overlap[i][j] += ct->tfs_vector_old[i][k]*ct->tfs_vector[j][k];
        }
    }

// check signs
    for (i = 0; i < ct->dim; i++)
    {
        if (ct->tfs_overlap[i][i] < 0.0)
        {
            for (k = 0; k < ct->dim; k++)
            {
                ct->tfs_vector[i][k]  *= -1.0;
                ct->tfs_overlap[k][i] *= -1.0;
            }
        }
    }

// calculate the adiabatic populations at the previous time step
    for (i = 0; i < ct->dim; i++)
        adiapop_old[i] = SQR(ct->tfs_popul[i]) + SQR(ct->tfs_popul[i+ct->dim]);

// calculate expansion coefficients in the adiabatic basis
    for (i = 0; i < ct->dim; i++)
    {
        ct->tfs_popul[i] = 0.0; ct->tfs_popul[i+ct->dim] = 0.0;

        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_popul[i]         += ct->tfs_vector[i][j]*ct->tfs_diab[j];
            ct->tfs_popul[i+ct->dim] += ct->tfs_vector[i][j]*ct->tfs_diab[j+ct->dim];
        }
    }

// record the adiabatic populations at new time
    dadiapopsum = 0.0;
    for (i = 0; i < ct->dim; i++)
    {
        adiapop[i] = SQR(ct->tfs_popul[i])+SQR(ct->tfs_popul[i+ct->dim]);

        dadiapop[i] = adiapop[i] - adiapop_old[i];

        if (dadiapop[i] < 0.0)
           dadiapopsum += dadiapop[i];
    }

// hopping probability
    probmax = 0.0;
    for (i = 0; i < ct->dim; i++)
    {
        ct->surf_prob[i] = 0.0;

        if (dadiapop[ct->surface] < 0.0 && dadiapop[i] > 0.0)
           ct->surf_prob[i] = dadiapop[i]/adiapop_old[ct->surface]*dadiapop[ct->surface]/dadiapopsum;
    }

//    for (i = 0; i < ct->dim; i++)
 //       printf("%2d  %15.5f  %15.5f \n", i, dadiapop[i]/adiapop_old[ct->surface], ct->surf_prob[i]);

    prev_surf = ct->surface;
    for (i = 0; i < ct->dim; i++)
    {
        if (ct->surf_prob[i] > threshold)
           ct->surface = i;
    }

// stochastic algorithm
    probsum = 0.0;
    random_number = drand48();
    if (prev_surf == ct->surface)
    {
    for (i = 0; i < ct->dim; i++)
    {
        if (ct->surf_prob[i] >= 0.0 && i != ct->surface)
        {
           probsum += ct->surf_prob[i];

           if (random_number < probsum)
           {
              hopsurf = i;

              dftb->do_ncv = 2;

              if (dftb->do_ncv == 2)
              {
                 random_number1 = drand48();

                 if (random_number1 > exp(-(ct->ev_adiab[hopsurf] - ct->ev_adiab[ct->surface])/ct->fermi_kt))
                    hopsurf = ct->surface;

                 dftb->do_ncv = 0;
              }

              if (hopsurf != ct->surface)
                 fprintf(f,"%8.2f %10d %17d %21f %19f \n",(istep*ct->rk_timestep+step*rk_timestep0)*1000./PS_TO_AU, ct->surface, hopsurf, ct->surf_prob[i], random_number);

              ct->surface = hopsurf;

              break;
           }
        }
    }
    }

// restricted decoherence correction
    dot_product = 0.0;
    for (i = 0; i < ct->dim; i++)
    {
        if (i != ct->surface)
        {
            decay_time = (1.0 + 0.1/ekintot)/fabs(ct->ev_adiab[i] - ct->ev_adiab[ct->surface]);

            ct->tfs_popul[i]         *= exp(-ct->rk_timestep/decay_time);
            ct->tfs_popul[i+ct->dim] *= exp(-ct->rk_timestep/decay_time);
            dot_product              += SQR(ct->tfs_popul[i])+SQR(ct->tfs_popul[i+ct->dim]);
        }
    }
    curr_surf_fac = sqrt((1.0 - dot_product)/(SQR(ct->tfs_popul[ct->surface]) + SQR(ct->tfs_popul[ct->dim+ct->surface])));
    ct->tfs_popul[ct->surface]         *= curr_surf_fac;
    ct->tfs_popul[ct->surface+ct->dim] *= curr_surf_fac;

// recalc. diab. from adiab after decoherence corr.
    for (i = 0; i < ct->dim; i++)
    {
        ct->tfs_diab[i] = 0.0; ct->tfs_diab[i+ct->dim] = 0.0;

        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_diab[i]         += ct->tfs_popul[j]*ct->tfs_vector[j][i];
            ct->tfs_diab[i+ct->dim] += ct->tfs_popul[j+ct->dim]*ct->tfs_vector[j][i];
        }
    }
    }
    ct->rk_timestep = rk_timestep0;

    printf("Current Surface: %d\n", ct->surface);

// charge wave function
    for (i = 0; i < ct->dim; i++)
        ct->wf[i] = ct->tfs_vector[ct->surface][i];

//store hamiltonian as old_ham
    for (i = 0; i < SQR(ct->dim); i++)
        ct->tfl_old_ham[i] = *(ct->hamiltonian[0]+i);

    return 0;
}
*/

#if GMX_MPI
int do_decoherence_induced_fssh(int step, charge_transfer_t *ct, dftb_t *dftb, t_state *state_global, t_mdatoms *mdatoms, FILE *f,
                                MPI_Comm ct_mpi_comm, int ct_mpi_rank, int ct_mpi_size)
#else
int do_decoherence_induced_fssh(int step, charge_transfer_t *ct, dftb_t *dftb, t_state *state_global, t_mdatoms *mdatoms, FILE *f)
#endif
{
// Declaration:
//
//                                          -----------Weiwei Xie on July. 2019
// Ref: J.Chem.Phys. 137, 22A545 (2012)


    int    i, j, k, l, m, n, tmp, prev_surf, old_surf, hopsurf, new_surf, surftmp, iatom, istep, nstep, *ind, indmax, *indtmp, *indtmp1, *map;
    double *ham, *old_ham, **ham_full, *new_ham, *adiapop_old, *adiapop, *initprob, *probtmp1, *probtmp2, *diapop, *diatmp;
    double tdse_timestep = 1.e-5*PS_TO_AU;
    double rk_timestep0, probsum, random_number, cur_weight;
    double popul_norm, overlapmax, overlaptmp1, overlaptmp2, adiapoptmp0, probtmp0, probtmp, diapoptmp;
    double nc, kai1, kai2, kai, delta, aa, bb, cc, dd, ee;
    double decay_time, curr_surf_fac, trivial_cross_type;
    double dot_product, ekintot;
    double delta_min;


    old_ham  = ct->tfl_old_ham_adiab;
    new_ham  = ct->hamiltonian_adiab;
    ham      = ct->tfl_mean_ham;
    ham_full = ct->tfl_mean_ham_full;

    snew(diapop, ct->dim);
    snew(diatmp, 2*ct->dim);
    snew(adiapop_old, ct->dim);
    snew(adiapop, ct->dim);
    snew(initprob, ct->dim);
    snew(probtmp1, ct->dim);
    snew(probtmp2, ct->dim);

    snew(ind, ct->dim);
    snew(indtmp, ct->dim);
    snew(indtmp1, ct->dim);
    snew(map, ct->dim);

// set up initial surface
    if (ct->tfs_initialization_step)
    {
       printf("Input wavefunction is:\n");
       for (j = 0; j < ct->dim; j++)
           printf("%f ", ct->wf[j]);
       printf("\n");

       for (i = 0; i < SQR(ct->dim); i++)
       {
           ct->tfl_old_ham[i] = *(ct->hamiltonian[0]+i);
           ham[i] = ct->tfl_old_ham[i];
       }

// diagonalize this Hamiltonian
       dsyev(ct->dim, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);

// record eigenvector
       for (i = 0; i < ct->dim; i++)
           for (j = 0; j < ct->dim; j++)
               ct->tfs_vector[i][j] = ham[i * ct->dim + j];

// determine inital surface using eigenvector as a criterion
       probsum = 0.;
       random_number = drand48();
       for (i = 0; i < ct->dim; i++)
       {
           initprob[i] = 0.;
           for (j = 0; j < ct->dim; j++)
               initprob[i] += ct->tfs_vector[i][j]*ct->wf[j];
           initprob[i] = SQR(initprob[i]);

           probsum += initprob[i];
           if (random_number < probsum)
           {
               ct->surface = i;
               break;
           }
       }

       for (i = 0; i < 2*ct->dim; i++)
       {
           ct->tfs_popul[i] = 0.0; ct->tfs_diab[i]  = 0.0;
       }

       for (i = 0; i < ct->dim; i++)
           ct->tfs_diab[i] = ct->wf[i];

       ct->tfs_popul[ct->surface] = 1.0;

       ct->tfs_initialization_step = 0;

       return 0;
    }

// set up the old and new Hamiltonian for linear interpolation
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            old_ham[i*ct->dim+j] = ct->tfl_old_ham[i*ct->dim+j];
            new_ham[i*ct->dim+j] = ct->hamiltonian[i][j];
        }
    }

// store the vector at preceding time step
/*    for (i = 0; i < ct->dim; i++)
        for (j = 0; j < ct->dim; j++)
            ct->tfs_vector_old2[i][j] = ct->tfs_vector[i][j];
*/
// solve time-dependent Schroedinger equation
    nstep = (int)(ct->rk_timestep/tdse_timestep);
    rk_timestep0 = ct->rk_timestep;
    ct->rk_timestep = tdse_timestep;
    printf("Time step for TDSE: %f fs \n", ct->rk_timestep/PS_TO_AU*1000);

    for (istep = 1; istep < nstep+1; istep++)
    {
// linear interpolation for Hamiltonian
        cur_weight = (double)istep/(double)nstep;
        for (i = 0; i < ct->dim; i++)
            for (j = 0; j < ct->dim; j++)
                ham_full[i][j] = (1.0 - cur_weight)*ct->tfl_old_ham[i*ct->dim+j] + cur_weight*ct->hamiltonian[i][j];
// propagate the TDSE in the diabatic representation
        do_rksuite_diab(ct);
    }
    ct->rk_timestep = rk_timestep0;

// diagonalize this Hamiltonian
    for (i = 0; i < SQR(ct->dim); i++)
        ham[i] = new_ham[i];

    dsyev(ct->dim, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);

// store the eigenvectors at time t-dt and time t
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_vector_old[i][j] = ct->tfs_vector[i][j];
            ct->tfs_vector[i][j]     = ham[i*ct->dim + j];
        }
    }

// adiabatic wave function overlaps between time t-dt and t
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_overlap[i][j] = 0.0;

            for (k = 0; k < ct->dim; k++)
                ct->tfs_overlap[i][j] += ct->tfs_vector_old[i][k]*ct->tfs_vector[j][k];
        }
    }

// reassignment of states according to the maximum wave function overlaps
/*    for (i = 0; i < ct->dim; i++)
    {
        overlapmax = 0.0;
        for (j = 0; j < ct->dim; j++)
        {
            if (overlapmax < fabs(ct->tfs_overlap[j][i]))
            {
               overlapmax = fabs(ct->tfs_overlap[j][i]);
               indtmp[i] = j;
            }
        }
    }
*/
// reassignment of states according to the min-cost algorithm
    assignprob(ct, ind);

// calculate the adiabatic populations at the previous time step
    for (i = 0; i < ct->dim; i++)
        adiapop_old[i] = SQR(ct->tfs_popul[i]) + SQR(ct->tfs_popul[i+ct->dim]);

// hopping probabilities based on Tully's fewest-switches anzatz
    for (i = 0; i < ct->dim; i++)
    {
//     nonadibatic couplings based on Hammes-Schiffer and Tully, JCP 101, 4657 (1994)
//        nc = (ct->tfs_overlap[prev_surf][i] - ct->tfs_overlap[j][ct->surface])/2.0/ct->rk_timestep;
//     ndibatic couplings based norm-preserving interpolation by Benjamin, JPCL 5, 2351 (2014)
        overlaptmp1 = sqrt(fabs(1.0 - SQR(ct->tfs_overlap[i][i]) - SQR(ct->tfs_overlap[ct->surface][i])));
        overlaptmp2 = (- ct->tfs_overlap[i][ct->surface]*ct->tfs_overlap[i][i] - ct->tfs_overlap[ct->surface][ct->surface]*ct->tfs_overlap[ct->surface][i])/overlaptmp1;

        aa = - sin(acos(ct->tfs_overlap[i][i]) - asin(ct->tfs_overlap[i][ct->surface]))/(acos(ct->tfs_overlap[i][i]) - asin(ct->tfs_overlap[i][ct->surface]));
        bb =   sin(acos(ct->tfs_overlap[i][i]) + asin(ct->tfs_overlap[i][ct->surface]))/(acos(ct->tfs_overlap[i][i]) + asin(ct->tfs_overlap[i][ct->surface]));
        cc =   sin(acos(ct->tfs_overlap[ct->surface][ct->surface]) - asin(ct->tfs_overlap[ct->surface][i]))/(acos(ct->tfs_overlap[ct->surface][ct->surface]) - asin(ct->tfs_overlap[ct->surface][i]));
        dd =   sin(acos(ct->tfs_overlap[ct->surface][ct->surface]) + asin(ct->tfs_overlap[ct->surface][i]))/(acos(ct->tfs_overlap[ct->surface][ct->surface]) + asin(ct->tfs_overlap[ct->surface][i]));
        ee = 2.0*asin(overlaptmp1*(overlaptmp1*overlaptmp2*asin(overlaptmp1) + (sqrt((1.0 - SQR(overlaptmp1))*(1.0 - SQR(overlaptmp2))) - 1.0)*asin(overlaptmp2)))/(SQR(asin(overlaptmp1)) - SQR(asin(overlaptmp2)));

        nc = (acos(ct->tfs_overlap[i][i])*(aa + bb) + asin(ct->tfs_overlap[ct->surface][i])*(cc + dd) + ee)/2.0/ct->rk_timestep;

        probtmp1[i] = (2.0*(ct->tfs_popul[ct->surface]*ct->tfs_popul[i]+ct->tfs_popul[ct->dim+ct->surface]*ct->tfs_popul[ct->dim+i])*nc)*ct->rk_timestep/adiapop_old[ct->surface];

        if (probtmp1[i] < 0.0 || i == ct->surface)
           probtmp1[i] = 0.0;
    }

// calculate expansion coefficients in the adiabatic basis
    for (i = 0; i < ct->dim; i++)
    {
        ct->tfs_popul[i] = 0.0; ct->tfs_popul[i+ct->dim] = 0.0;

        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_popul[i]         += ct->tfs_vector[i][j]*ct->tfs_diab[j];
            ct->tfs_popul[i+ct->dim] += ct->tfs_vector[i][j]*ct->tfs_diab[j+ct->dim];
        }
    }

// record the adiabatic populations at new time
    for (i = 0; i < ct->dim; i++)
        adiapop[i] = SQR(ct->tfs_popul[i])+SQR(ct->tfs_popul[i+ct->dim]);

// SC-FSSH, avoid the numerical instability of the NCs
    prev_surf = ct->surface;
    for (i = 0; i < ct->dim; i++)
    {
        if (prev_surf == ind[i])
        {
           new_surf = i;
           break;
        }
    }


    if (ct->surface != new_surf)
    {
       probtmp1[new_surf] = (adiapop_old[ct->surface] - adiapop[ct->surface])/adiapop_old[ct->surface];

       for (i = 0; i < ct->dim; i++)
       {
           if (i != new_surf)
           {
              probtmp1[new_surf] -= probtmp1[i];
           }
       }

       printf("SC is actived: %d  %d  %10f  %10f\n", ct->surface, new_surf, (adiapop_old[ct->surface] - adiapop[ct->surface])/adiapop_old[ct->surface], probtmp1[new_surf]);
    }

    for (i = 0; i < ct->dim; i++)
        printf("%2d  %2d  %15.5f \n", i, ind[i], probtmp1[i]);

// stochastic algorithm
    probsum = 0.0;
    random_number = drand48();
    for (i = 0; i < ct->dim; i++)
    {
        ct->surf_prob[i] = probtmp1[i];

        if (ct->surf_prob[i] >= 0.0 && i != ct->surface)
        {
           probsum += ct->surf_prob[i];

           if (random_number < probsum)
           {
              hopsurf = i;
// calculate the nonadibatic coupling vector for velocity adjustment
              dftb->do_ncv = 1;

              if (dftb->do_ncv)
              {
                 ct->ncv_isurf = ct->surface; ct->ncv_jsurf = hopsurf;

#if GMX_MPI
                 get_MM_params(ct, dftb, ct_mpi_comm, ct_mpi_rank, ct_mpi_size);
#else
                 get_MM_params(ct, dftb);
#endif

                 aa = 0.0; bb = 0.0;
                 cc = ct->ev_adiab[hopsurf] - ct->ev_adiab[ct->surface];

                 for (iatom = 0; iatom < mdatoms->homenr; iatom++)
                 {
                     n = 0;
                     for (k = 0; k < ct->sites; k++)
                     {
                         for (l = 0; l < ct->site[k].atoms; l++)
                         {
                             if (iatom == ct->site[k].atom[l])
                             {
                                for (m = 0; m < DIM; m++)
                                {
                                    ct->tfs_ncv[k][l][m] = (real) (dftb->phase1[k].grad[l][m] + dftb->phase2.grad[n][m])/cc;

                                    aa += SQR(ct->tfs_ncv[k][l][m])/2./dftb->phase1[k].mass[l]/AMU_TO_AU;
                                    bb += ct->tfs_ncv[k][l][m]*state_global->v[iatom][m]*NM_TO_BOHR/PS_TO_AU;

                                }
                             }
                             n++;
                         }
                     }
                 }

// adjust velocity along the direction of nonadiabtic coupling vector
                 delta = SQR(bb) - 4.*aa*cc;
                 if (delta >= 0.0)
                 {
                    kai1 = (- bb + sqrt(delta))/2./aa;
                    kai2 = (- bb - sqrt(delta))/2./aa;

                    kai = (fabs(kai1) < fabs(kai2) ? kai1 : kai2);

                    printf("Hop from %d to %d succeeds! Adjust velocity along NCV with kai = %8.5f \n", ct->surface, hopsurf, kai);
                 }
                 else
                 {
                    kai = - bb/aa;
                    hopsurf = ct->surface;

                    printf("Hop from %d to %d fails! Reverse velocity along NCV with kai = %8.5f \n", ct->surface, i, kai);
                 }
              }

              if (hopsurf != ct->surface)
                 fprintf(f,"%8.2f %10d %17d %21f %19f \n",(step*ct->rk_timestep)*1000./PS_TO_AU, ct->surface, hopsurf, ct->surf_prob[i], random_number);

              ct->surface = hopsurf;

// velocity adjustment
              if (dftb->do_ncv)
              {
                 if (kai != 0.0)
                 {
                    for (iatom = 0; iatom < mdatoms->homenr; iatom++)
                    {
                        for (k = 0; k < ct->sites; k++)
                        {
                            for (l = 0; l < ct->site[k].atoms; l++)
                            {
                                if (iatom == ct->site[k].atom[l])
                                {
                                   for (m = 0; m < DIM; m++)
                                       state_global->v[iatom][m] += kai*ct->tfs_ncv[k][l][m]/dftb->phase1[k].mass[l]/AMU_TO_AU/NM_TO_BOHR*PS_TO_AU;
                                }
                            }
                        }
                    }
                 }
                 dftb->do_ncv = 0;
              }
              break;
           }
        }
    }

    printf("Current Surface: %d\n", ct->surface);

    ekintot = 0.0;
    for (iatom = 0; iatom < mdatoms->homenr; iatom++)
    {
        for (k = 0; k < ct->sites; k++)
        {
            for (l = 0; l < ct->site[k].atoms; l++)
            {
                if (iatom == ct->site[k].atom[l])
                {
                   for (m = 0; m < DIM; m++)
                       ekintot += 0.5*dftb->phase1[k].mass[l]*AMU_TO_AU*SQR(state_global->v[iatom][m]*NM_TO_BOHR/PS_TO_AU);
                }
            }
        }
    }
    printf("Kinetic energy is %f\n", ekintot);

// restricted decoherence correction
    dot_product = 0.0;
    if (adiapop[ct->surface] > 0.001)
    {
       for (i = 0; i < ct->dim; i++)
       {
           if (i != ct->surface)
           {
               decay_time = (1.0 + 0.1/ekintot)/fabs(ct->ev_adiab[i] - ct->ev_adiab[ct->surface]);

               ct->tfs_popul[i]         *= exp(-ct->rk_timestep/decay_time);
               ct->tfs_popul[i+ct->dim] *= exp(-ct->rk_timestep/decay_time);
               dot_product              += SQR(ct->tfs_popul[i])+SQR(ct->tfs_popul[i+ct->dim]);
           }
       }
       curr_surf_fac = sqrt((1.0 - dot_product)/(SQR(ct->tfs_popul[ct->surface]) + SQR(ct->tfs_popul[ct->dim+ct->surface])));
       ct->tfs_popul[ct->surface]         *= curr_surf_fac;
       ct->tfs_popul[ct->surface+ct->dim] *= curr_surf_fac;

// recalc. diab. from adiab after decoherence corr.
       for (i = 0; i < ct->dim; i++)
       {
           ct->tfs_diab[i] = 0.0; ct->tfs_diab[i+ct->dim] = 0.0;

           for (j = 0; j < ct->dim; j++)
           {
               ct->tfs_diab[i]         += ct->tfs_popul[j]*ct->tfs_vector[j][i];
               ct->tfs_diab[i+ct->dim] += ct->tfs_popul[j+ct->dim]*ct->tfs_vector[j][i];
           }
       }
    }

// charge wave function
    for (i = 0; i < ct->dim; i++)
        ct->wf[i] = ct->tfs_vector[ct->surface][i];

//store hamiltonian as old_ham
    for (i = 0; i < SQR(ct->dim); i++)
        ct->tfl_old_ham[i] = *(ct->hamiltonian[0]+i);

    return 0;
}


// new sf implemention needs to be added / done May2018

int do_tully_local(int step, charge_transfer_t *ct, dftb_t *dftb, t_state *state_global, t_mdatoms *mdatoms, FILE *f)
{

    //combination of flexible surface hopping (JPCL 2013, 4, 1888-1895) and modified hopping prob. calculation (JPCL 2014, 5, 713-719)  dholub difference found

    int    i, j, k, l;
    double old_energy, energy, energy1, energy2, *ham, *old_ham, *new_ham, **ham_full, fermi_energy, fermi_upper, fermi_lower, *initprob, probsum;
    double random_number, nc;                      // in the interval (0,1)
    double cumulated_probability;              // sum of ct->surf_prob[0..k]
    double popul_norm, prob_b;                         // sum of SQR(ct->tfs_popul) except SQR(ct->tfs_popul[ct->surface])
    double decay_time, current_surface_factor; // in the decoherence algorithm
    double dot_product;
    int    dt_step = 0;

    old_ham  = ct->tfl_old_ham_adiab;
    new_ham  = ct->hamiltonian_adiab;
    ham      = ct->tfl_mean_ham;
    ham_full = ct->tfl_mean_ham_full;

    snew(initprob, ct->dim);
// calculate the charges from the wave function  [JUMP3 turn on?]
//for (i=0; i<ct->dim; i++)
// ct->q_act[i] = ct->q_old[i] = SQR(ct->wf[i]);

    // copy as the first vector
    energy = 1.e10;


    ////////assemble hamiltonian

    k = 0;
    for (i = 0; i < ct->dim; i++)
    {
        if (ct->tfl_is_in_system[i])
        {
            l = 0;
            for (j = 0; j < ct->dim; j++)
            {
                if (ct->tfl_is_in_system[j])
                {
                    new_ham[k+l*ct->tfl_num_of_states] = ct->hamiltonian[i][j];
                    new_ham[l+k*ct->tfl_num_of_states] = ct->hamiltonian[j][i];
                    l++;
                }
            }
            k++;
        }
    }



    if (ct->tfs_initialization_step)
    {
        for (i = 0; i < SQR(ct->dim); i++)
        {
            ct->tfl_old_ham[i] = *(ct->hamiltonian[0]+i);
        }
    }

    k = 0;
    for (i = 0; i < ct->dim; i++)
    {
        if (ct->tfl_is_in_system[i])
        {
            l = 0;
            for (j = 0; j < ct->dim; j++)
            {
                if (ct->tfl_is_in_system[j])
                {
                    old_ham[k+l*ct->tfl_num_of_states] = ct->tfl_old_ham[i*ct->dim+j];
                    old_ham[l+k*ct->tfl_num_of_states] = ct->tfl_old_ham[j*ct->dim+i];
                    l++;
                }
            }
            k++;
        }
    }

    ///print out system
    printf("Current system:\n");
    for (i = 0; i < ct->dim; i++)
    {
        printf("%d ", ct->tfl_is_in_system[i]);
    }
    printf("\n");


    // printf("Hamiltonian:\n"); //debug
    //for(i=0;i<ct->dim;i++) for(j=i;j<ct->dim;j++) printf("%8.5f \n",ct->hamiltonian[i][j]); //debug
    ct->rk_timestep /= 100.0; //dholub stepsize for electronic propagation; tried  50tepps  for CT -> crash  Jun2018
    for (dt_step = 0; dt_step < 101; dt_step++)
    {

        double cur_weight = (double)dt_step/100.0;
        //build hamiltonian
        for (i = 0; i < SQR(ct->tfl_num_of_states); i++)
        {
            //i -> dtstep...
            ham[i] = (1.0-cur_weight)*old_ham[i]+cur_weight*new_ham[i];
        }

        for (i = 0; i < ct->dim; i++)
        {
            for (j = 0; j < ct->dim; j++)
            {
                ham_full[i][j] = (1.0-cur_weight)*ct->tfl_old_ham[i*ct->dim+j]+cur_weight*ct->hamiltonian[i][j];

            }
        }

        // diagonalize this Hamiltonian
        dsyev(ct->tfl_num_of_states, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);

        if (ct->tfs_initialization_step)
        {
            /* let us start with tfs_vector_old == tfs_vector */

            for (i = 0; i < ct->tfl_num_of_states_old; i++)
            {
                k = 0;
                for (j = 0; j < ct->dim; j++)
                {

                    if (ct->tfl_is_in_system[j])
                    {
                        ct->tfs_vector_old[i][j] = ham[i * ct->tfl_num_of_states_old + k]; k++;
                        ct->tfs_vector[i][j] = ct->tfs_vector_old[i][j];
                    }
                    else
                    {
                        ct->tfs_vector_old[i][j] = 0.;
                        ct->tfs_vector[i][j]     = 0.;
                    }
                }
            }
            for (; i < ct->dim; i++)
            {
                for (j = 0; j < ct->dim; j++)
                {
                    ct->tfs_vector_old[i][j] = 0.;
                    ct->tfs_vector[i][j]     = 0.;
                }
            }
// determine inital surface using eigenvector as a criterion
            probsum = 0.;
            random_number = drand48();
            for (i = 0; i < ct->dim; i++)
            {
                initprob[i] = 0.;
                for (j = 0; j < ct->dim; j++)
                {
                    initprob[i] += ct->tfs_vector_old[i][j]*ct->wf[j];
                }
                initprob[i] = SQR(initprob[i]);

                probsum += initprob[i];
                if (random_number < probsum)
                {
                    ct->surface = i;
                    break;
                }
            }
            for (i = 0; i < ct->dim; i++)
            {
                ct->tfs_diab[i] = ct->tfs_vector_old[ct->surface][i];
                ct->wf[i]       = ct->tfs_vector_old[ct->surface][i];
            }
            for (i = 0; i < 2*ct->dim; i++)
            {
                ct->tfs_popul[i] = 0.0;
            }
            ct->tfs_popul[ct->surface] = 1.0;

            ct->tfs_initialization_step = 0;

            ct->rk_timestep *= 100.0;

            return 0;
        }
        else
        {
            /* push tfs_vector to tfs_vector_old */
            for (i = 0; i < ct->dim; i++)
            {
                for (j = 0; j < ct->dim; j++)
                {
                    ct->tfs_vector_old[i][j] = ct->tfs_vector[i][j];
                }
            }
        }
        ////GOT HERE
        /* enter new eigenvectors to tfs_vector */
        for (i = 0; i < ct->tfl_num_of_states; i++)
        {
            k = 0;
            for (j = 0; j < ct->dim; j++)
            {
                if (ct->tfl_is_in_system[j])
                {
                    ct->tfs_vector[i][j] = ham[i * ct->tfl_num_of_states + k]; k++;
                }
                else
                {
                    ct->tfs_vector[i][j] = 0;
                }

            }
        }
        for (; i < ct->dim; i++)
        {
            for (j = 0; j < ct->dim; j++)
            {
                ct->tfs_vector[i][j] = 0;
            }
        }



        /////Attention!!: number of old and new sites may differ. Account for this!
        /* calculate tfs_overlap[i][j] = <tfs_vector_old[i] | tfs_vector[j]> */
        for (i = 0; i < ct->dim; i++)
        {
            for (j = 0; j < ct->dim; j++)
            {
                ct->tfs_overlap[i][j] = 0.;

                for (k = 0; k < ct->dim; k++)
                {
                    ct->tfs_overlap[i][j] += ct->tfs_vector_old[i][k] * ct->tfs_vector[j][k];

                }
            }
        }

        /////GOT HERE
        //////////////////////////////////////check signs
        for (i = 0; i < ct->dim; i++)
        {
            int sim_ind = 0;
            dot_product = 0.0;

            for (k = 0; k < ct->dim; k++)
            {
                if (fabs(ct->tfs_overlap[k][i]) > dot_product)
                {
                    sim_ind = k; dot_product = fabs(ct->tfs_overlap[k][i]);
                }
            }

            if (ct->tfs_overlap[sim_ind][i] < 0.00)
            {
                for (k = 0; k < ct->dim; k++)
                {
                    ct->tfs_vector[i][k] *= -1.0; ct->tfs_overlap[k][i] *= -1.0;
                }
            }
        }

        /////GOT HERE

        ////////store old adiab. population of surface
        double ca_re = ct->tfs_popul[ct->surface];
        double ca_im = ct->tfs_popul[ct->surface+ct->dim];

        //////       Statingpoint after inisialization:   propagate in diabatic representation
        do_rksuite_diab(ct);

        // should work like this, but think again:
        // calc. new popul (remember to store some values from the old one)
        dot_product = 0.0;
        for (i = 0; i < ct->tfl_num_of_states; i++)
        {
            ct->tfs_popul[i] = 0.0; ct->tfs_popul[i+ct->dim] = 0.0;
            for (j = 0; j < ct->dim; j++)
            {
                ct->tfs_popul[i]         += ct->tfs_vector[i][j]*ct->tfs_diab[j];
                ct->tfs_popul[i+ct->dim] += ct->tfs_vector[i][j]*ct->tfs_diab[j+ct->dim];

            }
            dot_product += SQR(ct->tfs_popul[i])+SQR(ct->tfs_popul[i+ct->dim]);
        }
        for (i = 0; i < ct->tfl_num_of_states; i++) //account for limit basis
        {
            ct->tfs_popul[i] /= sqrt(dot_product); ct->tfs_popul[i+ct->dim] /= sqrt(dot_product);
        }

        for (j = 0; j < ct->dim; j++) //diab expansion may have changed due to inexact projection
        {
            ct->tfs_diab[j]         = 0.0;
            ct->tfs_diab[j+ct->dim] = 0.0;
            for (i = 0; i < ct->tfl_num_of_states; i++)
            {
                ct->tfs_diab[j]         += ct->tfs_popul[i]*ct->tfs_vector[i][j];
                ct->tfs_diab[j+ct->dim] += ct->tfs_popul[i+ct->dim]*ct->tfs_vector[i][j];
            }
        }


        int old_surface = ct->surface;

        //////check for trivial crossing:

//	printf("<90");
        if (fabs(ct->tfs_overlap[ct->surface][ct->surface]) < 0.92)  // check  falue 0.85 +- mapping surface
        {
            //|| (SQR(ca_im)+SQR(ca_re))/(SQR(ct->tfs_popul[ct->surface])+SQR(ct->tfs_popul[ct->surface+ct->dim]))>2.0){//1.4 tried too

            //printf("Detected trivial crossing!\n");
            ////using overlap, may use populations instead

            dot_product = 0.0;
            int sim_ind = 0;
            /*   for(k=0;k<ct->dim;k++){
                if(fabs(ct->tfs_overlap[ct->surface][k])>dot_product){
                  sim_ind=k; dot_product=fabs(ct->tfs_overlap[ct->surface][k]);
                }
                }  */
            for (k = 0; k < ct->tfl_num_of_states; k++)
            {
                if ((SQR(ct->tfs_popul[k])+SQR(ct->tfs_popul[k+ct->dim])) > dot_product)
                {
                    dot_product = (SQR(ct->tfs_popul[k])+SQR(ct->tfs_popul[k+ct->dim]));
                    sim_ind     = k;
                }
            }


            printf("Changed surface to %d from %d according to population\n",sim_ind,ct->surface);
            printf("Overlap is: %f\n",ct->tfs_overlap[ct->surface][sim_ind]);
            ct->surface = sim_ind;

        }

        ////////////////


        popul_norm = 0.;
        for (k = 0; k < ct->tfl_num_of_states; k++)
        {
            popul_norm += SQR(ct->tfs_popul[k]) + SQR(ct->tfs_popul[ct->dim + k]);
        }
        // printf("Population norm: %9.6f\n",popul_norm);
        //correct norm. RK not unitary!!!!
        popul_norm = sqrt(popul_norm);
        for (k = 0; k < ct->tfl_num_of_states; k++)
        {
            ct->tfs_popul[k] /= popul_norm; ct->tfs_popul[k+ct->dim] /= popul_norm;
        }

/////////do hopping//////////////////////

        /////calc. overall hopping prob.
        double hop_prob = 1.0-(SQR(ct->tfs_popul[ct->surface])+SQR(ct->tfs_popul[ct->surface+ct->dim]))/(SQR(ca_im)+SQR(ca_re));
        double prob_a   = SQR(ct->tfs_popul[ct->surface])+SQR(ct->tfs_popul[ct->surface+ct->dim]);

        //if(hop_prob>1e-5) printf("Dada: %f",hop_prob); //debug

        ///conventional form
        ct->surf_prob[ct->surface] = 0.;
        double hopping_prob_sum = 0.0;
        for (i = 0; i < ct->tfl_num_of_states; i++)
        {
            if (i != ct->surface)
            {
                //   ct->surf_prob[j]=prob_prefactor[j]*ct->tfs_overlap[ct->surface][j];
                ct->surf_prob[i] = 2.0*(ct->tfs_popul[ct->surface]*ct->tfs_popul[i]+ct->tfs_popul[ct->dim+ct->surface]*ct->tfs_popul[ct->dim+i])*0.5*(ct->tfs_overlap[old_surface][i]-ct->tfs_overlap[i][ct->surface])/prob_a;//ct->tfs_overlap[i][ct->surface];
                if (ct->surf_prob[i] > 0.0)
                {
                    hopping_prob_sum += ct->surf_prob[i];
                }
                if (ct->ev_adiab[i] > ct->ev_adiab[ct->surface])
                {
                    ct->surf_prob[i] *= exp(-(ct->ev_adiab[i]-ct->ev_adiab[ct->surface])/ct->fermi_kt);
                }
            }
        }

        //   if(hop_prob>hopping_prob_sum){
        //   for (j=0; j<ct->tfl_num_of_states; j++){
        //     if (j != ct->surface&&ct->surf_prob[j]>0.0) ct->surf_prob[j]*=hop_prob/hopping_prob_sum;
        //     }
        //  }
        hopping_prob_sum = 0.0;
        for (j = 0; j < ct->tfl_num_of_states; j++)
        {
            if (j != ct->surface && ct->surf_prob[j] > 0.0)
            {
                hopping_prob_sum += ct->surf_prob[j];
            }
        }


        if (hopping_prob_sum > 1.0)
        {
            // printf("Need to rescale probabilities!\n");
            for (j = 0; j < ct->tfl_num_of_states; j++)
            {
                if (j != ct->surface && ct->surf_prob[j] > 0.0)
                {
                    ct->surf_prob[j] /= hopping_prob_sum;
                }
            }
        }
        int prev_surf = ct->surface;
        ///////////////////////////////hopping routine/////////////////////////////////////////
        cumulated_probability = 0.;
        //random_number = rand()/(double)RAND_MAX;
        random_number = drand48();
        for (i = 0; i < ct->tfl_num_of_states; i++)
        {
            if (i != ct->surface && ct->surf_prob[i] > 0.)
            {
                cumulated_probability += ct->surf_prob[i];
                //if(ct->surf_prob[i]>1e-5) //debug
                //printf("Bla: %f\n",ct->surf_prob[i]); //debug
                if (cumulated_probability > random_number)
                {
                    ct->surface = i;
                    break;
                }
            }
        }
        //////////////decoherence correction/////////////////////////////////////////
       /* if (!ct->use_strong_dec)
        {     		reduce this option and use always just the DEC way 16.08.18 dholub	*/
            dot_product = 0.0;
            for (i = 0; i < ct->tfl_num_of_states; i++)
            {
                if (i != ct->surface)
                {
                    double decay_time = TFS_DECAY_CONSTANT/fabs(ct->ev_adiab[i]-ct->ev_adiab[ct->surface]); //dholub changed TFS_DECAY_CONSTANT in charge_trans.h in testsystem without larger effects JUMP
                    ct->tfs_popul[i]         *= exp(-ct->rk_timestep/decay_time);
                    ct->tfs_popul[i+ct->dim] *= exp(-ct->rk_timestep/decay_time);
                    dot_product              += SQR(ct->tfs_popul[i])+SQR(ct->tfs_popul[i+ct->dim]);
                }
            }
            double curr_surf_fac = sqrt((1.0-dot_product)/(SQR(ct->tfs_popul[ct->surface])+SQR(ct->tfs_popul[ct->dim+ct->surface])));
            ct->tfs_popul[ct->surface]         *= curr_surf_fac;
            ct->tfs_popul[ct->surface+ct->dim] *= curr_surf_fac;

            //recalc. diab. from adiab after decoherence corr.
            for (i = 0; i < ct->dim; i++) //bad loop order, accept for now
            {
                ct->tfs_diab[i] = 0.0; ct->tfs_diab[i+ct->dim] = 0.0;
                for (j = 0; j < ct->tfl_num_of_states; j++)
                {
                    ct->tfs_diab[i]         += ct->tfs_popul[j]*ct->tfs_vector[j][i];
                    ct->tfs_diab[i+ct->dim] += ct->tfs_popul[j+ct->dim]*ct->tfs_vector[j][i];
                }
            }

        //}
       /* else
        {
            if (prev_surf != ct->surface)
            {
                for (i = 0; i < ct->tfl_num_of_states; i++)
                {
                    ct->tfs_popul[i]         = 0.0;
                    ct->tfs_popul[i+ct->dim] = 0.0;
                }
                ct->tfs_popul[ct->surface] = 1.0;
            }

            //recalc. diab. from adiab after decoherence corr.

            for (i = 0; i < ct->dim; i++) //bad loop order, accept for now

            {
                ct->tfs_diab[i] = 0.0; ct->tfs_diab[i+ct->dim] = 0.0;
                for (j = 0; j < ct->tfl_num_of_states; j++)
                {
                    ct->tfs_diab[i]         += ct->tfs_popul[j]*ct->tfs_vector[j][i];
                    ct->tfs_diab[i+ct->dim] += ct->tfs_popul[j+ct->dim]*ct->tfs_vector[j][i];
                }
            }

        } */
        //////////////////////////////////////////////////////////////////////////////


    }//End loop over dt_step

    ct->rk_timestep *= 100.0;


    ////////////////Some output:
    printf("TFL Info: ");
    //if(prev_surf!=ct->surface){ printf("Did hopp this step! New surface: %d!\n",ct->surface); }
    //else printf("Did NOT hopp this step!\n");

    //printf("Overlap between active states at prev. and current step: %f %f\n",ct->tfs_overlap[ct->surface][prev_surf],ct->tfs_overlap[old_surface][ct->surface]);

    /////////////////////// tried to comment dholub 1.2019
    ////update system //////////////////////////////////////////////////////////////////7

    //////////determine which sites need to be included
    ct->tfl_num_of_states_old = ct->tfl_num_of_states;
    for (i = 0; i < ct->dim; i++)
    {
        ct->tfl_is_in_system_old[i] = ct->tfl_is_in_system[i];                //copy current to old
    }
    ct->tfl_num_of_states = 0;
    for (i = 0; i < ct->dim; i++)
    {
        //calc. <phi_surface|i>
        double phixi = 0.;
        for (j = 0; j < ct->dim; j++)
        {
            phixi += ct->tfs_vector[ct->surface][j]*ct->hamiltonian[i][j];
        }
        phixi = fabs(phixi);
        //calc. <phi_surface|H|phi_surface>-<i|H|j|i>
        double tfl_en_diff = fabs(ct->ev_adiab[ct->surface]-ct->hamiltonian[i][i]);
        if (phixi/tfl_en_diff > TFL_RC && fabs(phixi) > 0.001)
        {
            ct->tfl_is_in_system[i] = 1; ct->tfl_num_of_states++;
        }
        else
        {
            ct->tfl_is_in_system[i] = 0;
        }
    }
    int tfl_did_sys_change = 0;
    for (i = 0; i < ct->dim; i++)
    {
        if (ct->tfl_is_in_system[i] != ct->tfl_is_in_system_old[i])
        {
            tfl_did_sys_change = 1; break;
        }
    }
    //recalc. Eigenstates and choose appropriate surface
    if (tfl_did_sys_change)
    {

        //rebuild hamiltonian with new set of sites
        k = 0;
        for (i = 0; i < ct->dim; i++)
        {
            if (ct->tfl_is_in_system[i])
            {
                l = 0;
                for (j = 0; j < ct->dim; j++)
                {
                    if (ct->tfl_is_in_system[j])
                    {
                        ham[k+l*ct->tfl_num_of_states] = ct->hamiltonian[i][j];
                        ham[l+k*ct->tfl_num_of_states] = ct->hamiltonian[j][i];
                        l++;
                    }
                }
                k++;
            }
        }

//////////////////////////////////////////

        dsyev(ct->tfl_num_of_states, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);

        printf("System changed! WF before system adaption was: \n");
        double *old_surf_vec = (double*)malloc(sizeof(double)*ct->dim); //vec. to store surf. state.
        for (i = 0; i < ct->dim; i++)
        {
            old_surf_vec[i] = ct->tfs_vector[ct->surface][i];
        }
        for (i = 0; i < ct->dim; i++)
        {
            printf("%f ", old_surf_vec[i]);
        }
        printf("\n");

        for (i = 0; i < ct->tfl_num_of_states; i++)
        {
            k = 0;
            for (j = 0; j < ct->dim; j++)
            {
                if (ct->tfl_is_in_system[j])
                {
                    ct->tfs_vector[i][j] = ham[i * ct->tfl_num_of_states + k]; k++;
                }
                else
                {
                    ct->tfs_vector[i][j] = 0;
                }

            }
        }
        for (; i < ct->dim; i++)
        {
            for (j = 0; j < ct->dim; j++)
            {
                ct->tfs_vector[i][j] = 0;
            }
        }



        /////find which new state is surface by selecting the one with max overlap
        double maximum_overlap = 0.;
        int    max_index       = 0;
        for (i = 0; i < ct->tfl_num_of_states; i++)
        {
            dot_product = 0;
            for (k = 0; k < ct->dim; k++)
            {
                dot_product += old_surf_vec[k]*ct->tfs_vector[i][k];
            }
            //dot_product=fabs(dot_product);
            if (fabs(dot_product) > maximum_overlap)
            {
                maximum_overlap = fabs(dot_product); max_index = i;
            }

        }///////////////
        ct->surface = max_index;
        //if(dot_product<0.) for(j=0;j<ct->dim;j++) ct->tfs_vector[ct->surface][j]*=-1.0; //check sign
        printf("Surface in new system is: %d \n", ct->surface);



        /////update populations
        dot_product = 0.0;
        for (i = 0; i < ct->tfl_num_of_states; i++)
        {
            ct->tfs_popul[i]         = 0.0;
            ct->tfs_popul[i+ct->dim] = 0.0;
            for (j = 0; j < ct->dim; j++)
            {
                ct->tfs_popul[i]         += ct->tfs_vector[i][j]*ct->tfs_diab[j];
                ct->tfs_popul[i+ct->dim] += ct->tfs_vector[i][j]*ct->tfs_diab[j+ct->dim];
            }
            dot_product += SQR(ct->tfs_popul[i])+SQR(ct->tfs_popul[i+ct->dim]);
        }
        for (i = 0; i < ct->tfl_num_of_states; i++) //fix norm as projection may not be exact
        {
            ct->tfs_popul[i]         /= sqrt(dot_product);
            ct->tfs_popul[i+ct->dim] /= sqrt(dot_product);
        }
        for (j = 0; j < ct->dim; j++) //diab expansion may have changed due to inexact projection
        {
            ct->tfs_diab[j]         = 0.0;
            ct->tfs_diab[j+ct->dim] = 0.0;
            for (i = 0; i < ct->tfl_num_of_states; i++)
            {
                ct->tfs_diab[j]         += ct->tfs_popul[i]*ct->tfs_vector[i][j];
                ct->tfs_diab[j+ct->dim] += ct->tfs_popul[i+ct->dim]*ct->tfs_vector[i][j];
            }
        }



        free(old_surf_vec); //clear mem.

        //take care of populations //collapse wave function for now
        // for(i=0;i<2*ct->dim;i++) ct->tfs_popul[i]=0.0;
        //ct->tfs_popul[ct->surface]=1.0;
        ///

    }
    /////end recalc.

////////////////////end update system


    // assign the new wave function to the eigenvector ct->surface

    printf("The new wavefunction is:\n");
    for (i = 0; i < ct->dim; i++)
    {
        ct->wf[i] = ct->tfs_vector[ct->surface][i];
        printf("%f ", ct->wf[i]);
    }
    printf("\n");
/*
    printf("TFL Info: System did ");
    if (tfl_did_sys_change == 1)
    {
        printf("change.\n");
    }
    else
    {
        printf("not change.\n");
    }
*/
//    old end of comment dholub 1.2019 */

    printf("Current Surface: %d\n", ct->surface);
//debug //print eigenvectors
//  printf("Eigenvectors: \n");
// for(i=0;i<ct->dim;i++){
// for(j=0;j<ct->dim;j++) printf("%9.5f ",ct->tfs_vector[i][j]);
// printf("\n");
// }



    ////store hamiltonian as old_ham
    for (i = 0; i < SQR(ct->dim); i++)
    {
        ct->tfl_old_ham[i] = *(ct->hamiltonian[0]+i);
    }

    return 0;
}


int do_persico_diabatic_sfhopping(charge_transfer_t *ct, dftb_broyden_t *broyd, double *fermi_coeff, FILE *f, FILE *f2)
{
    int    i, j, k, step;
    double old_energy, energy, energy1, energy2, *ham, fermi_energy, fermi_upper, fermi_lower, q_sum;
    double prob_factor;
    double random_number;                      // in the interval (0,1)
    double cumulated_probability;              // sum of ct->surf_prob[0..k]
    double popul_norm;                         // sum of SQR(ct->tfs_popul) except SQR(ct->tfs_popul[ct->surface])
    double decay_time, current_surface_factor; // in the decoherence algorithm
    double dot_product;

    ham = ct->hamiltonian_adiab;

    // calculate the charges from the wave function
    for (i = 0; i < ct->dim; i++)
    {
        ct->q_act[i] = ct->q_old[i] = SQR(ct->wf[i]);
    }

    // copy as the first vector

    energy = 1.e10;

    // SCC cycle
    for (step = 0;; step++)
    {
        printf("%d\n", step);
        // too many iterations? break and indicate failure
        if (step >= MAXITER_BROYDEN)
        {
            fprintf(f, "Broyden/Fermi failed to converge!");
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
        // mix all eigenvectors with appropriate coefficients

        // calculate the "Fermi energy", store the Fermi coefficients in fermi_coeff
        fermi_lower = ct->ev_adiab[0] - 0.1;
        fermi_upper = ct->ev_adiab[ct->dim - 1] + 0.1;
        while (fermi_upper - fermi_lower > FERMI_CONVERG)
        {
            fermi_energy = 0.5 * (fermi_lower + fermi_upper);
            if (fermi_coef_sum(fermi_energy, ct->dim, ct->ev_adiab, fermi_coeff, ct->fermi_kt) > 1.0)
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
        energy1    = energy2 = 0.0;
        // indices i and j run over the sites
        for (i = 0; i < ct->dim; i++)
        {
            for (j = 0; j < ct->dim; j++)
            {
                // TB Hamiltonian
                // index k runs over the adiabatic states (eigenvectors)
                for (k = 0; k < ct->dim; k++)
                {
                    energy1 += fermi_coeff[k] * ham[k * ct->dim + i] * ham[k * ct->dim + j] * ct->hamiltonian[i][j];
                }
                // Hubbard / gamma terms
                energy2 += 0.5 * ct->hubbard[i][j] * ct->q_act[i] * ct->q_act[j];
            }
        }
        energy = energy1 + energy2;

        if (step > 1 && fabs(energy - old_energy) < BROYDEN_SCFTOL)
        {
            //fprintf(f, "Broyd %2d ", step);
            fprintf(f, " %2d ", step);
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

    // print energies of adiabatic states - these are pure, but obtained with non-pure (Fermi-mixed) input state
    //fprintf(f, " Eigenvalues:");
    for (i = 0; i < ct->dim; i++)
    {
        fprintf(f, " %8.5f", ct->ev_adiab[i]);
    }
    //fprintf(f, "    ");

    /* surface hopping tried according to Tully JCP 93, 1061 (1990)
     * "fewest switches" algorithm
     * diabatization according to Persico JCP 114, 10608 (2001)
     * approximate decoherence correction adapted from Persico JCP 126, 134114 (2007)
     */

    /* exp_imag_matrix - has to be written yet !!! */

    /* THIS WILL REMAIN! */
    if (ct->tfs_initialization_step)
    {
        /* let us start with tfs_vector_old == tfs_vector */
        for (i = 0; i < ct->dim; i++)
        {
            for (j = 0; j < ct->dim; j++)
            {
                ct->tfs_vector_old[i][j] = ham[i * ct->dim + j];
            }
        }
        for (i = 0; i < ct->dim; i++)
        {
            ct->ev_adiab_old[i] = ct->ev_adiab[i];
        }
    }
    else
    {
        /* push tfs_vector to tfs_vector_old */
        for (i = 0; i < ct->dim; i++)
        {
            for (j = 0; j < ct->dim; j++)
            {
                ct->tfs_vector_old[i][j] = ct->tfs_vector[i][j];
            }
        }
    }

    /* PUSH THE OLD POPULATIONS TO tfs_popul_old */
    for (i = 0; i < ct->dim; i++)
    {
        ct->tfs_popul_old[i]           = ct->tfs_popul[i];
        ct->tfs_popul_old[ct->dim + i] = ct->tfs_popul[ct->dim + i];
    }

    /* THIS WILL REMAIN, TOO
       MAYBE JUST DO NOT PERFORM THE CHECKING FOR SIGN CHANGE! */
    /* enter new eigenvectors to tfs_vector */
    for (i = 0; i < ct->dim; i++)
    {
        /* calculate the dot product with tfs_vector_old[i] */
        //  dot_product = 0.;
        //  for (k=0; k<ct->sites; k++)
        //    dot_product += ct->tfs_vector_old[i][k] * ham[i * ct->sites + k];
        /* for positive dot_product, take the vector as it is
         * for negative dot_product, take (-1) * the vector
         */
        //  for (k=0; k<ct->sites; k++)
        //    ct->tfs_vector[i][k] = dot_product > 0 ? ham[i * ct->sites + k] : - ham[i * ct->sites + k];
        for (k = 0; k < ct->dim; k++)
        {
            ct->tfs_vector[i][k] = ham[i * ct->dim + k];
        }
    }

    /* THIS WILL REMAIN - tfs_overlap[i][j] corresponds to matrix T[i][j] ! */
    /* calculate tfs_overlap[i][j] = <tfs_vector_old[i] | tfs_vector[j]> */
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_overlap[i][j] = 0.;
            for (k = 0; k < ct->dim; k++)
            {
                ct->tfs_overlap[i][j] += ct->tfs_vector_old[i][k] * ct->tfs_vector[j][k];
            }
        }
    }

    /* print out the overlap matrix */
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            fprintf(f2, "%9.5f", ct->tfs_overlap[i][j]);
        }
    }
    fprintf(f2, "\n");

    /* DO THE PERSICO STUFF HERE!
     * BEGINNING FROM THE SECOND STEP OF THE SIMULATION
     */
    if (ct->tfs_initialization_step)
    {
        ct->tfs_initialization_step = 0;
    }
    else
    {

        /* SET UP THE AVERAGE DIABATIC HAMILTONIAN Z = 1/2 * (E(0) + T(dt) * E(dt) * T^t(dt) */
        for (i = 0; i < ct->dim; i++)
        {
            for (j = 0; j < ct->dim; j++)
            {
                if (i == j)
                {
                    ct->per_diab_hamiltonian[i][i] = ct->ev_adiab_old[i] / 2.;
                }
                else
                {
                    ct->per_diab_hamiltonian[i][j] = 0.;
                }
                for (k = 0; k < ct->dim; k++)
                {
                    ct->per_diab_hamiltonian[i][j] += ct->tfs_overlap[i][k] * ct->ev_adiab[k] * ct->tfs_overlap[j][k] / 2.;
                }
            }
        }

        /* CONSTRUCT THE PROPAGATION OPERATOR exp[-i*Z*dt] ! */
        exp_imag_matrix(ct->per_diab_hamiltonian, ct->per_propag_operator, ct->rk_timestep, ct->dim, ct->per_arrays);

        /* CONSTRUCT THE TRANSFORMATION MATRIX U = T^t * exp[-i*Z*dt] */
        for (i = 0; i < ct->dim; i++)
        {
            for (j = 0; j < ct->dim; j++)
            {
                ct->per_transformator[i][j][0] = 0.;
                ct->per_transformator[i][j][1] = 0.;
                for (k = 0; k < ct->dim; k++)
                {
                    ct->per_transformator[i][j][0] += ct->tfs_overlap[k][i] * ct->per_propag_operator[k][j][0];
                    ct->per_transformator[i][j][1] += ct->tfs_overlap[k][i] * ct->per_propag_operator[k][j][1];
                }
            }
        }

        /* OBTAIN THE NEW POPULATIONS */
        for (i = 0; i < ct->dim; i++)
        {
            ct->tfs_popul[i] = ct->tfs_popul[ct->dim + i] = 0.;
            for (j = 0; j < ct->dim; j++)
            {
                ct->tfs_popul[i]             += ct->per_transformator[i][j][0] * ct->tfs_popul_old[j]
                    - ct->per_transformator[i][j][1] * ct->tfs_popul_old[ct->dim + j];
                ct->tfs_popul[ct->dim + i] += ct->per_transformator[i][j][0] * ct->tfs_popul_old[ct->dim + j]
                    + ct->per_transformator[i][j][1] * ct->tfs_popul_old[j];
            }
        }

        /* print out new populations */
        for (k = 0; k < ct->dim; k++)
        {
            fprintf(f, "%9.6f", SQR(ct->tfs_popul[k]) + SQR(ct->tfs_popul[ct->dim + k]));
        }

        /* we are on ct->surface at the moment
         * now, go over all of the other states,
         * and calculate switching probability */

        /* CONTINUE HERE WITH SWITCHING PROBABILITIES! */
        ct->surf_prob[ct->surface] = -1.; /* to skip this surface in the test */
        if (SQR(ct->per_transformator[ct->surface][ct->surface][0]) + SQR(ct->per_transformator[ct->surface][ct->surface][1]) > 0.999999)
        {
            /* the population of current state does not really change
             * numerical problems are to be expected
             * therefore, skip the calculation of probabilities
             * no surface hop will be performed in this step!
             */
            for (j = 0; j < ct->dim; j++)
            {
                if (j != ct->surface)
                {
                    ct->surf_prob[j] = -0.1;
                }
            }
        }
        else
        {
            /* otherwise, there is change in the population
             * and it is safe to calculate surface-hopping probabilities!
             */
            prob_factor = (SQR(ct->tfs_popul[ct->surface]) + SQR(ct->tfs_popul[ct->dim + ct->surface]) - SQR(ct->tfs_popul_old[ct->surface]) - SQR(ct->tfs_popul_old[ct->dim + ct->surface]))
                / ( SQR(ct->tfs_popul[ct->surface]) + SQR(ct->tfs_popul[ct->dim + ct->surface])
                    - ( (ct->per_transformator[ct->surface][ct->surface][0] * ct->tfs_popul_old[ct->surface] - ct->per_transformator[ct->surface][ct->surface][1] * ct->tfs_popul_old[ct->dim + ct->surface])
                        * ct->tfs_popul[ct->surface]
                        + (ct->per_transformator[ct->surface][ct->surface][0] * ct->tfs_popul_old[ct->dim + ct->surface] + ct->per_transformator[ct->surface][ct->surface][1] * ct->tfs_popul_old[ct->surface])
                        * ct->tfs_popul[ct->dim + ct->surface] ) );
            for (j = 0; j < ct->dim; j++)
            {
                if (j != ct->surface)
                {
                    ct->surf_prob[j] = -((ct->per_transformator[ct->surface][j][0] * ct->tfs_popul_old[j] - ct->per_transformator[ct->surface][j][1] * ct->tfs_popul_old[ct->dim + j]) * ct->tfs_popul[ct->surface]
                                         + (ct->per_transformator[ct->surface][j][0] * ct->tfs_popul_old[ct->dim + j] + ct->per_transformator[ct->surface][j][1] * ct->tfs_popul_old[j]) * ct->tfs_popul[ct->dim + ct->surface])
                        * prob_factor;
                }
            }
            /* ct->surf_prob[j] = -2 * (ct->tfs_popul[ct->surface] * ct->tfs_popul[j] + ct->tfs_popul[ct->sites + ct->surface] * ct->tfs_popul[ct->sites + j])
                                  / (SQR(ct->tfs_popul[ct->surface]) + SQR(ct->tfs_popul[ct->sites + ct->surface]))
             * ct->tfs_overlap[j][ct->surface]; */
        }

        /* write out the probabilities */
        for (k = 0; k < ct->dim; k++)
        {
            fprintf(f, "%10.6f", ct->surf_prob[k]);
        }

        /* generate a random number */
        random_number         = rand()/(double)RAND_MAX;
        cumulated_probability = 0.;

        /* and determine the surface to continue on */
        for (j = 0; j < ct->dim; j++)
        {
            if (j != ct->surface && ct->surf_prob[j] > 0.)
            {
                cumulated_probability += ct->surf_prob[j];
                if (cumulated_probability > random_number)
                {
                    ct->surface = j;
                    break;
                }
            }
        }

        /* decoherence correction
         * proposed by Truhlar JCTC 1, 527 (2005)
         * recommended by Persico 126, 134114 (2007)
         * simplified to:
         * tau_km = hbar / |E_k - E_m| * TFS_DECAY_CONSTANT
         * TFS_DECAY_CONSTANT chosen to be 3. (should be quite insensitive, though)
         * /

           popul_norm = 0.;
           for (j=0; j<ct->sites; j++) if (j != ct->surface) {
           decay_time = TFS_DECAY_CONSTANT / fabs(ct->ev_adiab[j] - ct->ev_adiab[ct->surface]);
           ct->tfs_popul[j]             *= exp( - ct->rk_timestep / decay_time);
           ct->tfs_popul[ct->sites + j] *= exp( - ct->rk_timestep / decay_time);
           popul_norm += SQR(ct->tfs_popul[j]) + SQR(ct->tfs_popul[ct->sites + j]);
           }
           current_surface_factor = sqrt((1. - popul_norm) / (SQR(ct->tfs_popul[ct->surface]) + SQR(ct->tfs_popul[ct->sites + ct->surface])));
           ct->tfs_popul[ct->surface]             *= current_surface_factor;
           ct->tfs_popul[ct->sites + ct->surface] *= current_surface_factor;
         */

    } /* end if initialization step else */

    /* assign the new wave function to the eigenvector ct->surface */
    for (i = 0; i < ct->dim; i++)
    {
        ct->wf[i] = ct->tfs_vector[ct->surface][i];
    }

    /* print out current state */
    fprintf(f, " %d", ct->surface);
    for (i = 0; i < ct->dim; i++)
    {
        fprintf(f, "%8.4f", SQR(ct->wf[i]));
    }
    fprintf(f, "\n");

    /* push the current ev_adiab to ev_adiab_old, to be used in the next step */
    for (i = 0; i < ct->dim; i++)
    {
        ct->ev_adiab_old[i] = ct->ev_adiab[i];
    }

    return 0;
}


int do_persico_diabatic_sfhopping_new(charge_transfer_t *ct, dftb_broyden_t *broyd, double *fermi_coeff, FILE *f, FILE *f2, FILE *f3)
{
    int             i, j, k, l, m, n, step;
    double          old_energy, energy, energy1, energy2, *ham, fermi_energy, fermi_upper, fermi_lower;
    double          prob_factor;
    double          random_number;                      // in the interval (0,1)
    double          cumulated_probability;              // sum of ct->surf_prob[0..k]
    double          popul_norm;                         // sum of SQR(ct->tfs_popul) except SQR(ct->tfs_popul[ct->surface])
    double          decay_time, current_surface_factor; // in the decoherence algorithm
    double          dot_product, tot_prob;
    double          dot_prod_r, dot_prod_i;
    /* the names of these variables are taken from Newton-X, file sh_locdiab.f90, subroutine trans_ld */
    int             surf, surf_old; /* originally: istati, istati_old */
    double          popk0, popk1, rekk, rekl, rnum, rden;
    double         *acoer0, *acoei0, *acoer1, *acoei1;
    twodoubles      ac0, ac1, vkk, vkl;
    float           maxval;
    const float     eps_switch = 0.5;//1.e-9; take always state with better overlap. note: overlap should be close to 1 or 0. if offdiag is far from 0 smaller time step is needed
    static double **pop_flow = NULL;
    static int      lwork = 0;
    static double  *hdx = NULL, *ex = NULL, *work = NULL, **rotr = NULL, **roti = NULL, **rotsr = NULL, **rotsi = NULL, **rotr_new = NULL, **roti_new = NULL;
    double          timetau, ld_sk, ld_ck, ld_fac;
    int             ld_istep;
    const int       ld_nstep = 20;

    ham      = ct->hamiltonian_adiab;
    n        = ct->dim;
    surf_old = surf = ct->surface;
    acoer0   = ct->tfs_popul_old;
    acoei0   = ct->tfs_popul_old + n;
    acoer1   = ct->tfs_popul;
    acoei1   = ct->tfs_popul + n;

    if (pop_flow == NULL)
    {
        snew(pop_flow, n);
        snew(pop_flow[0], SQR(n));
        for (i = 1; i < n; i++)
        {
            pop_flow[i] = pop_flow[0] + i * n;
        }
    }

    if (hdx == NULL)
    {
        snew(hdx, SQR(n));
    }

    if (ex == NULL)
    {
        snew(ex, n);
    }

    if (work == NULL)
    {
        lwork = (n > 3) ? SQR(n) : 5;
        snew(work, lwork);
    }

    if (rotr == NULL)
    {
        snew(rotr, n);
        snew(rotr[0], SQR(n));
        for (i = 1; i < n; i++)
        {
            rotr[i]  = rotr[0] + i * n;
        }
        snew(roti, n);
        snew(roti[0], SQR(n));
        for (i = 1; i < n; i++)
        {
            roti[i]  = roti[0] + i * n;
        }
        snew(rotsr, n);
        snew(rotsr[0], SQR(n));
        for (i = 1; i < n; i++)
        {
            rotsr[i] = rotsr[0] + i * n;
        }
        snew(rotsi, n);
        snew(rotsi[0], SQR(n));
        for (i = 1; i < n; i++)
        {
            rotsi[i] = rotsi[0] + i * n;
        }
        snew(rotr_new, n);
        snew(rotr_new[0], SQR(n));
        for (i = 1; i < n; i++)
        {
            rotr_new[i]  = rotr_new[0] + i * n;
        }
        snew(roti_new, n);
        snew(roti_new[0], SQR(n));
        for (i = 1; i < n; i++)
        {
            roti_new[i]  = roti_new[0] + i * n;
        }
    }

    // calculate the charges from the wave function
    // and copy as the first vector
    // NO! set to zero instead - perhaps improve efficiency...
    for (i = 0; i < ct->dim; i++)
    {
        ct->q_act[i] = ct->q_old[i] = 0.0; //alex: test for identical dimers with a large spacial separation
    }
    // ct->q_act[i] = ct->q_old[i] = SQR(ct->wf[i]);
    // NO! take them from the previous iteration

    energy = 1.e10;

    // SCC cycle
    for (step = 0;; step++)
    {

        // too many iterations? break and indicate failure
        if (step >= MAXITER_BROYDEN)
        {
            fprintf(f, "Broyden/Fermi failed to converge!");
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
        // mix all eigenvectors with appropriate coefficients

        // calculate the "Fermi energy", store the Fermi coefficients in fermi_coeff
        fermi_lower = ct->ev_adiab[0] - 0.1;
        fermi_upper = ct->ev_adiab[ct->dim - 1] + 0.1;
        while (fermi_upper - fermi_lower > FERMI_CONVERG)
        {
            fermi_energy = 0.5 * (fermi_lower + fermi_upper);
            if (fermi_coef_sum(fermi_energy, ct->dim, ct->ev_adiab, fermi_coeff, ct->fermi_kt) > 1.0)
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
/*
    for (i=0; i<ct->dim; i++){
    printf("qact: %lf\n", ct->q_act[i]);
    printf("fermicoeff: %lf\n", fermi_coeff[i]);
    for (j=0; j<ct->dim; j++){
        printf("ham %d %d : %lf\n", i,j, ham[j * ct->dim + i]);
    }
    }
 */
        // calculate the energy and check convergence
        old_energy = energy;
        energy1    = energy2 = 0.0;
        // indices i and j run over the sites
        for (i = 0; i < ct->dim; i++)
        {
            for (j = 0; j < ct->dim; j++)
            {
                // TB Hamiltonian
                // index k runs over the adiabatic states (eigenvectors)
                for (k = 0; k < ct->dim; k++)
                {
                    energy1 += fermi_coeff[k] * ham[k * ct->dim + i] * ham[k * ct->dim + j] * ct->hamiltonian[i][j];
                }
                // Hubbard / gamma terms
                energy2 += 0.5 * ct->hubbard[i][j] * ct->q_act[i] * ct->q_act[j];
            }
        }
        energy = energy1 + energy2;

        if (step > 1 && fabs(energy - old_energy) < BROYDEN_SCFTOL)
        {
            //fprintf(f, "Broyd %2d ", step);
            fprintf(f, " %2d ", step);
            break;
        }

        // mix the charges - Broyden
        broyden(step, BROYDEN_ALMIX, ct->dim, ct->q_old, ct->q_act, broyd);
        for (i = 0; i < ct->dim; i++)
        {
            ct->q_act[i] = ct->q_old[i];
        }

    } /* end SCC cycle */

    // print energies of adiabatic states - these are pure, but obtained with non-pure (Fermi-mixed) input state
    //fprintf(f, " Eigenvalues:");
    for (i = 0; i < ct->dim; i++)
    {
        fprintf(f, "%8.5f", ct->ev_adiab[i]);
    }
    //fprintf(f, "    ");

    /* surface hopping tried according to Tully JCP 93, 1061 (1990)
     * "fewest switches" algorithm
     * diabatization according to Persico JCP 114, 10608 (2001)
     * approximate decoherence correction adapted from Persico JCP 126, 134114 (2007)
     */

    /* exp_imag_matrix - has to be written yet !!! */

    /* 1. set up the state vectors */

    /*   a. back up state vectors from the previous step, tfs_vector_old*/

    if (ct->tfs_initialization_step)
    {
        /* let us start with tfs_vector_old == tfs_vector */
        for (i = 0; i < ct->dim; i++)
        {
            for (j = 0; j < ct->dim; j++)
            {
                ct->tfs_vector_old[i][j] = ham[i * ct->dim + j];
            }
        }
        for (i = 0; i < ct->dim; i++)
        {
            ct->ev_adiab_old[i] = ct->ev_adiab[i];
        }

        //ALEX: project input-WF on adiabatic states. i.e. take diabatic state as initial condiation instead of lowest adiabatic state.
        //for (i=0; i<ct->dim; i++) // for every surface
        //  dot_product = 0.;
        //  for (k=0; k<ct->dim; k++) // for every FO
        //    dot_product += ct->tfs_vector_old[i][k] * ct->wf[k];
        //  ct->tfs_popul[i]=dot_product;
    }
    else
    {
        /* push tfs_vector to tfs_vector_old */
        for (i = 0; i < ct->dim; i++)
        {
            for (j = 0; j < ct->dim; j++)
            {
                ct->tfs_vector_old[i][j] = ct->tfs_vector[i][j];
            }
        }
    }

    /*   b. back up populations (complex coefficients) of states from the previous step, tfs_popul_old */

    for (i = 0; i < ct->dim; i++)
    {
        ct->tfs_popul_old[i]             = ct->tfs_popul[i];
        ct->tfs_popul_old[ct->dim + i]   = ct->tfs_popul[ct->dim + i];
    }

    /*  c. save the new state vectors */

    for (i = 0; i < ct->dim; i++)
    {
        /* calculate the dot product with tfs_vector_old[i],
           to check if we shall invert the sign
         */
        dot_product = 0.;
        for (k = 0; k < ct->dim; k++)
        {
            dot_product += ct->tfs_vector_old[i][k] * ham[i * ct->dim + k];
        }
        /* for positive dot_product, take the vector as it is
         * for negative dot_product, take (-1) * the vector
         */
        for (k = 0; k < ct->dim; k++)
        {
            ct->tfs_vector[i][k] = dot_product > 0 ? ham[i * ct->dim + k] : -ham[i * ct->dim + k];
        }
        //ct->tfs_vector[i][k] =  ham[i * ct->dim + k];
    }

    /* 2. calculate the overlap of state vectors
            tfs_overlap[i][j] is the matrix T[i][j] = <tfs_vector_old[i] | tfs_vector[j]>
     */

    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_overlap[i][j] = 0.;
            for (k = 0; k < ct->dim; k++)
            {
                ct->tfs_overlap[i][j] += ct->tfs_vector_old[i][k] * ct->tfs_vector[j][k];
            }
        }
    }

    /* print out the overlap matrix */
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            fprintf(f2, "%9.5f", ct->tfs_overlap[i][j]);
        }
    }
    fprintf(f2, "\n");

    /* print out the state vectors */
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            fprintf(f3, "%9.5f", ct->tfs_vector[i][j]);
        }
    }
    fprintf(f3, "\n");

    /* check if the matrix is unitary - in our case orthogonal (real matrix)
       printf("Testing the orthogonality of T - this should be a unit matrix:\n");
       for (i=0; i<n; i++) {
       for (j=0; j<n; j++) {
        dot_product = 0.;
        for (k=0; k<n; k++)
          dot_product += ct->tfs_overlap[i][k] * ct->tfs_overlap[j][k];
        printf(" %9.6f", dot_product);
       }
       printf("\n");
       }
     */

    /* 3. do the propagation of the wave function / populations
            this is equivalent to the propag_ld function in NewtonX */

    /*    a. do not do it in the first step of the simulation */

    if (ct->tfs_initialization_step)
    {
        ct->tfs_initialization_step = 0;
    }
    else
    {

        /*    b. set up the diabatic hamiltonian at the end of the time step Hdia(dt) = T(dt) * E(dt) * T^t(dt) : Eq. B7 */
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                ct->per_diab_hamiltonian[i][j] = 0.;
                for (k = 0; k < n; k++)
                {
                    ct->per_diab_hamiltonian[i][j] += ct->tfs_overlap[i][k] * ct->ev_adiab[k] * ct->tfs_overlap[j][k];
                }
            }
        }

        /*    c. divide the propagation into smaller elementary steps */
        timetau = ct->rk_timestep / (double) ld_nstep;
        for (ld_istep = 0; ld_istep < ld_nstep; ld_istep++)
        {
            /* matrix Hdx                 = (Hdia(istep*timetau) + Hdia((istep+1)*tau)) / 2
               matrix Hdia(istep*timetau) = E(0) + (2*istep-1)/(2*ld_nstep) * (Hdia(dt) - E(0))
             */
            ld_fac = (double) (2. * ld_istep + 1) / (double) (2 * ld_nstep);
            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    hdx[i + j * n] = timetau * (ld_fac * ct->per_diab_hamiltonian[i][j] + (i == j ? (1. - ld_fac) * ct->ev_adiab_old[i] : 0.));
                }
            }
            /* eigenproblem of exp[-i * Hdx * timetau] */
            dsyev(n, hdx, ex, work, lwork);
            /* construct the rotation matrix (real + imaginary parts) for the elementary step */
            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    rotsr[i][j] = 0.;
                    rotsi[i][j] = 0.;
                }
            }
            for (k = 0; k < n; k++)
            {
                ld_sk = sin(ex[k]);
                ld_ck = cos(ex[k]);
                for (i = 0; i < n; i++)
                {
                    for (j = 0; j < n; j++)
                    {
                        rotsr[i][j] += hdx[i + k * n] * hdx[j + k * n] * ld_ck;
                        rotsi[i][j] -= hdx[i + k * n] * hdx[j + k * n] * ld_sk;
                    }
                }
            }
            /* accumulate the overall rotation matrix (real + imaginary) */
            if (ld_istep == 0)
            {
                for (i = 0; i < n; i++)
                {
                    for (j = 0; j < n; j++)
                    {
                        rotr[i][j] = rotsr[i][j];
                        roti[i][j] = rotsi[i][j];
                    }
                }
            }
            else
            {
                /* tw  = rotsr * rotr
                   tx  = rotsi * roti
                   ty  = rotsr * roti
                   tz  = rotsi * rotr
                   rotr = tw - tx
                   roti = ty + tz
                 */
                for (i = 0; i < n; i++)
                {
                    for (j = 0; j < n; j++)
                    {
                        rotr_new[i][j] = 0.;
                        roti_new[i][j] = 0.;
                        for (k = 0; k < n; k++)
                        {
                            rotr_new[i][j] += rotsr[i][k] * rotr[k][j] - rotsi[i][k] * roti[k][j];
                            roti_new[i][j] += rotsr[i][k] * roti[k][j] + rotsi[i][k] * rotr[k][j];
                        }
                    }
                }
                for (i = 0; i < n; i++)
                {
                    for (j = 0; j < n; j++)
                    {
                        rotr[i][j] = rotr_new[i][j];
                        roti[i][j] = roti_new[i][j];
                    }
                }
            }
            /*    d. update the coefficients in the adiabatic representation */
            for (i = 0; i < n; i++)
            {
                ct->tfs_popul[i]     = 0.;
                ct->tfs_popul[i + n] = 0.;
                for (j = 0; j < n; j++)
                {
                    for (k = 0; k < n; k++)
                    {
                        ct->tfs_popul[i]     += ct->tfs_overlap[j][i] *
                            (ct->tfs_popul_old[k] * rotr[j][k] - ct->tfs_popul_old[k + n] * roti[j][k]);
                        ct->tfs_popul[i + n] += ct->tfs_overlap[j][i] *
                            (ct->tfs_popul_old[k] * roti[j][k] + ct->tfs_popul_old[k + n] * rotr[j][k]);
                    }
                }
            }
        }
        /*
           // print out the total rotation matrix
           printf("===============================================\n");
           printf("Rotation matrix for coefficients in the adiabatic representation:\n");
           for (i=0; i<n; i++) {
           for (j=0; j<n; j++)
            printf("  %9.6f + %9.6fi", rotr[i][j], roti[i][j]);
           printf("\n");
           }
           printf("Testing the unitarity of the rotation matrix - this should be a unit matrix:\n");
           for (i=0; i<n; i++) {
           for (j=0; j<n; j++) {
            dot_prod_r = dot_prod_i = 0.;
            for (k=0; k<n; k++) {
              dot_prod_r += rotr[i][k] * rotr[j][k] + roti[i][k] * roti[j][k];
              dot_prod_i += roti[i][k] * rotr[j][k] - rotr[i][k] * roti[j][k];
           }
            printf("%10.6f%10.6f", dot_prod_r, dot_prod_i);
           }
           printf("\n");
           }
         */
        /*
           printf("TFS_POPUL_NEW");
           for (i=0; i<n; i++)
           printf("%10.6f%10.6f",  ct->tfs_popul[i],  ct->tfs_popul[i + n]);
           printf("\n");
         */

        /* FOR TESTING PURPOSES, COMPARE WITH OUR PREVIOUS IMPLEMENTATION!

           // average diabatic Hamiltonian: Z = 1/2 * (E(0) + T(dt) * E(dt) * T^t(dt)
           for (i=0; i<ct->sites; i++) {
           for (j=0; j<ct->sites; j++) {
            if (i==j) {
              ct->per_diab_hamiltonian[i][i] = ct->ev_adiab_old[i] / 2.;
            } else {
              ct->per_diab_hamiltonian[i][j] = 0.;
            }
            for (k=0; k<ct->sites; k++)
              ct->per_diab_hamiltonian[i][j] += ct->tfs_overlap[i][k] * ct->ev_adiab[k] * ct->tfs_overlap[j][k] / 2.;
           }
           }
           // construct the propagation operator exp[-i*Z*dt] !
           exp_imag_matrix(ct->per_diab_hamiltonian, ct->per_propag_operator, ct->rk_timestep, ct->sites, ct->per_arrays);
           // print out the propagation operator
           printf("===============================================\n");
           printf("Original alternative - the propagation operator:\n");
           for (i=0; i<n; i++) {
           for (j=0; j<n; j++)
            printf("  %9.6f + %9.6fi", ct->per_propag_operator[i][j][0], ct->per_propag_operator[i][j][1]);
           printf("\n");
           }
           printf("===============================================\n");

           printf("Testing the unitarity of the propagation operator:\n");
           for (i=0; i<n; i++) {
           for (j=0; j<n; j++) {
            dot_prod_r = dot_prod_i = 0.;
            for (k=0; k<n; k++) {
              dot_prod_r += ct->per_propag_operator[i][k][0] * ct->per_propag_operator[j][k][0] + ct->per_propag_operator[i][k][1] * ct->per_propag_operator[j][k][1];
              dot_prod_i += ct->per_propag_operator[i][k][1] * ct->per_propag_operator[j][k][0] - ct->per_propag_operator[i][k][0] * ct->per_propag_operator[j][k][1];
           }
            printf("%10.6f%10.6f", dot_prod_r, dot_prod_i);
           }
           printf("\n");
           }
         */

        /* construct the transformation matrix U = T^t * exp[-i*Z*dt]
           for (i=0; i<ct->sites; i++)
           for (j=0; j<ct->sites; j++) {
            ct->per_transformator[i][j][0] = 0.;
            ct->per_transformator[i][j][1] = 0.;
            for (k=0; k<ct->sites; k++) {
              ct->per_transformator[i][j][0] += ct->tfs_overlap[k][i] * ct->per_propag_operator[k][j][0];
              ct->per_transformator[i][j][1] += ct->tfs_overlap[k][i] * ct->per_propag_operator[k][j][1];
            }
           }
         */
        /* then, new populations would be obtained as
           for (i=0; i<n; i++) {
           ct->tfs_popul[i] = ct->tfs_popul[n + i] = 0.;
           for (j=0; j<n; j++) {
            ct->tfs_popul[i]     += ct->per_transformator[i][j][0] * ct->tfs_popul_old[j]
                                  - ct->per_transformator[i][j][1] * ct->tfs_popul_old[n + j];
            ct->tfs_popul[n + i] += ct->per_transformator[i][j][0] * ct->tfs_popul_old[n + j]
         + ct->per_transformator[i][j][1] * ct->tfs_popul_old[j];
           }
           }
         */
        /*
           printf("TFS_POPUL_OLD");
           for (i=0; i<n; i++)
           printf("%10.6f%10.6f",  ct->tfs_popul[i],  ct->tfs_popul[i + n]);
           printf("\n");
         */
        /* END TESTING PURPOSES */



        /* GO ON HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */



        /* print out new populations */
        for (k = 0; k < ct->dim; k++)
        {
            fprintf(f, "%9.6f", SQR(ct->tfs_popul[k]) + SQR(ct->tfs_popul[ct->dim + k]));
        }

        /* we are on ct->surface at the moment
         * now, go over all of the other states,
         * and calculate switching probability */

        /* the following is inspired/taken from Newton-X, file sh_locdiab.f90, subroutine trans_ld */

        /* save the value U_kk */
        //vkk[0] = ct->per_transformator[surf][surf][0];
        //vkk[1] = ct->per_transformator[surf][surf][1];
        vkk[0] = vkk[1] = 0.;
        for (j = 0; j < n; j++)
        {
            vkk[0] += ct->tfs_overlap[j][surf] * rotr[j][surf];
            vkk[1] += ct->tfs_overlap[j][surf] * roti[j][surf];
        }

        //printf("\n CHECK VKK: %9.6f %9.6f VS %9.6f %9.6f\n",
        //  ct->per_transformator[surf][surf][0], ct->per_transformator[surf][surf][1], vkk[0], vkk[1]);

        /* has a switch of states just occured? */
        surf_old = surf;
        if (fabs(ct->tfs_overlap[surf][surf]) < eps_switch)
        {
            printf("state switching has occurred! Overlap is only: %f  \n", ct->tfs_overlap[surf][surf]);
            for (j = 0; j < n; j++)
            {
                ct->surf_prob[j] = 0.;
            }
            maxval = 0.;
            for (j = 0; j < n; j++)
            {
                if (fabs((float) ct->tfs_overlap[j][surf_old]) > maxval)
                {
                    surf   = j;
                    maxval = fabs((float) ct->tfs_overlap[j][surf_old]);
                }
            }
            printf("state %d probably became state %d! Overlap: %f  \n", surf_old, surf, ct->tfs_overlap[surf][surf_old]);
        }

        /* rename / evaluate simple variables */
        ac0[0] = acoer0[surf_old];
        ac0[1] = acoei0[surf_old];
        ac1[0] = acoer1[surf];
        ac1[1] = acoei1[surf];
        popk0  = SQR(ac0[0]) + SQR(ac0[1]);
        popk1  = SQR(ac1[0]) + SQR(ac1[1]);
        /* rekk = U_kk * A_k(0) * A_k^*(dt) */
        rekk = vkk[0] * (ac0[0]*ac1[0] + ac0[1]*ac1[1]) + vkk[1] * (ac0[0]*ac1[1] - ac0[1]*ac1[0]);
        //printf(" CHECK REKK: %9.6f\n", rekk);
        rnum = (popk1 - popk0) / popk0;
        rden = popk1 - rekk;

        /* if population of occupied surface increases,
         * then no hopping!
         */
        if (popk1 > popk0)
        {
            for (j = 0; j < n; j++)
            {
                ct->surf_prob[j] = -1.;
            }
            /* otherwise, calculate the probabilities */
        }
        else
        {
            for (j = 0; j < n; j++)
            {
                /* no hopping from state 'surf' to itself, obviously */
                if (j == surf)
                {
                    ct->surf_prob[j] = -1.;
                    /* calculate the probability for all of the states except 'surf' */
                }
                else
                {
                    //vkl[0] = ct->per_transformator[surf][j][0];
                    //vkl[1] = ct->per_transformator[surf][j][1];
                    vkl[0] = vkl[1] = 0.;
                    for (m = 0; m < n; m++)
                    {
                        vkl[0] += ct->tfs_overlap[m][surf] * rotr[m][j];
                        vkl[1] += ct->tfs_overlap[m][surf] * roti[m][j];
                    }
                    rekl             = vkl[0] * (acoer0[j]*ac1[0] + acoei0[j]*ac1[1]) + vkl[1] * (acoer0[j]*ac1[1] - acoei0[j]*ac1[0]);
                    ct->surf_prob[j] = -rekl * rnum / rden;
                }
            }
        }
        /* probabilities are done at this point */
        for (k = 0; k < n; k++)
        {
            //vkk[0] = ct->per_transformator[k][k][0];
            //vkk[1] = ct->per_transformator[k][k][1];
            vkk[0] = vkk[1] = 0.;
            for (j = 0; j < n; j++)
            {
                vkk[0] += ct->tfs_overlap[j][k] * rotr[j][k];
                vkk[1] += ct->tfs_overlap[j][k] * roti[j][k];
            }
            ac0[0] = acoer0[k];
            ac0[1] = acoei0[k];
            ac1[0] = acoer1[k];
            ac1[1] = acoei1[k];
            popk0  = SQR(ac0[0]) + SQR(ac0[1]);
            popk1  = SQR(ac1[0]) + SQR(ac1[1]);
            /* rekk = U_kk * A_k(0) * A_k^*(dt) */
            rekk = vkk[0] * (ac0[0]*ac1[0] + ac0[1]*ac1[1]) + vkk[1] * (ac0[0]*ac1[1] - ac0[1]*ac1[0]);
            rnum = (popk1 - popk0) / popk0;
            rden = popk1 - rekk;
            for (l = 0; l < n; l++)
            {
                //vkl[0] = ct->per_transformator[k][l][0];
                //vkl[1] = ct->per_transformator[k][l][1];
                vkl[0] = vkl[1] = 0.;
                for (m = 0; m < n; m++)
                {
                    vkl[0] += ct->tfs_overlap[m][k] * rotr[m][l];
                    vkl[1] += ct->tfs_overlap[m][k] * roti[m][l];
                }
                rekl           = vkl[0] * (acoer0[l]*ac1[0] + acoei0[l]*ac1[1]) + vkl[1] * (acoer0[l]*ac1[1] - acoei0[l]*ac1[0]);
                pop_flow[k][l] = -rekl * rnum * popk0 / rden;
                //printf(" %8.5f", pop_flow[k][l]);
            }
        }
        printf("\n");

        tot_prob = 0.;
        for (j = 0; j < n; j++)
        {
            tot_prob += (ct->surf_prob[j] > 0. ? ct->surf_prob[j] : 0.);
        }
        if (tot_prob > 1.)
        {
            for (j = 0; j < n; j++)
            {
                if (ct->surf_prob[j] > 0.)
                {
                    ct->surf_prob[j] /= tot_prob;
                }
            }
        }

        /* END OF PROBABILITIES! */

        /* write out the probabilities */
        for (k = 0; k < ct->dim; k++)
        {
            fprintf(f, "%10.6f", ct->surf_prob[k]);
        }

        /* generate a random number */
        random_number         = rand()/(double)RAND_MAX;
        cumulated_probability = 0.;

        /* and determine the surface to continue on */
        for (j = 0; j < ct->dim; j++)
        {
            if (j != ct->surface && ct->surf_prob[j] > 0.)
            {
                cumulated_probability += ct->surf_prob[j];
                if (cumulated_probability > random_number)
                {
                    ct->surface = j;
                    break;
                }
            }
        }

        /* decoherence correction
         * proposed by Truhlar JCTC 1, 527 (2005)
         * recommended by Persico 126, 134114 (2007)
         * simplified to:
         * tau_km = hbar / |E_k - E_m| * TFS_DECAY_CONSTANT
         * TFS_DECAY_CONSTANT chosen to be 3. (should be quite insensitive, though)
         */

        if (ct->decoherence)
        {
            popul_norm = 0.;
            for (j = 0; j < n; j++)
            {
                if (j != ct->surface)
                {
                    decay_time                    = TFS_DECAY_CONSTANT / fabs(ct->ev_adiab[j] - ct->ev_adiab[ct->surface]);
                    ct->tfs_popul[j]             *= exp( -ct->rk_timestep / decay_time);
                    ct->tfs_popul[ct->dim + j]   *= exp( -ct->rk_timestep / decay_time);
                    popul_norm                   += SQR(ct->tfs_popul[j]) + SQR(ct->tfs_popul[n + j]);
                }
            }
            current_surface_factor                  = sqrt((1. - popul_norm) / (SQR(ct->tfs_popul[ct->surface]) + SQR(ct->tfs_popul[n + ct->surface])));
            ct->tfs_popul[ct->surface]             *= current_surface_factor;
            ct->tfs_popul[ct->dim + ct->surface]   *= current_surface_factor;
        }

    } /* end if initialization step else */

    /* assign the new wave function to the eigenvector ct->surface */
    for (i = 0; i < ct->dim; i++)
    {
        ct->wf[i] = ct->tfs_vector[ct->surface][i];
    }

    /* print out current state */
    fprintf(f, " %d", ct->surface);
    for (i = 0; i < ct->dim; i++)
    {
        fprintf(f, "%8.4f", SQR(ct->wf[i]));
    }
    fprintf(f, "\n");

    /* push the current ev_adiab to ev_adiab_old, to be used in the next step */
    for (i = 0; i < ct->dim; i++)
    {
        ct->ev_adiab_old[i] = ct->ev_adiab[i];
    }

    return 0;
}


int do_prezhdo_sfhopping(charge_transfer_t *ct, dftb_broyden_t *broyd, double *fermi_coeff, FILE *f, FILE *f2, FILE *f3)
{
    int             i, j, k, l, m, n, step;
    double          old_energy, energy, energy1, energy2, *ham, fermi_energy, fermi_upper, fermi_lower;
    double          prob_factor;
    double          random_number;                      // in the interval (0,1)
    double          cumulated_probability;              // sum of ct->surf_prob[0..k]
    double          popul_norm;                         // sum of SQR(ct->tfs_popul) except SQR(ct->tfs_popul[ct->surface])
    double          decay_time, current_surface_factor; // in the decoherence algorithm
    double          dot_product, tot_prob;
    double          dot_prod_r, dot_prod_i;
    /* the names of these variables are taken from Newton-X, file sh_locdiab.f90, subroutine trans_ld */
    int             surf, surf_old; /* originally: istati, istati_old */
    double          popk0, popk1, rekk, rekl, rnum, rden;
    double         *acoer0, *acoei0, *acoer1, *acoei1;
    twodoubles      ac0, ac1, vkk, vkl;
    float           maxval;
    const float     eps_switch = 1.e-9;
    static double **pop_flow = NULL;
    static int      lwork = 0;
    static double  *hdx = NULL, *ex = NULL, *work = NULL, **rotr = NULL, **roti = NULL, **rotsr = NULL, **rotsi = NULL, **rotr_new = NULL, **roti_new = NULL;
    double          timetau, ld_sk, ld_ck, ld_fac;
    int             ld_istep;
    const int       ld_nstep = 20;
    double          diff_flow, gross_flow;
    double          a, b, ski;

    ham      = ct->hamiltonian_adiab;
    n        = ct->dim;
    surf_old = surf = ct->surface;
    acoer0   = ct->tfs_popul_old;
    acoei0   = ct->tfs_popul_old + n;
    acoer1   = ct->tfs_popul;
    acoei1   = ct->tfs_popul + n;

    if (pop_flow == NULL)
    {
        snew(pop_flow, n);
        snew(pop_flow[0], SQR(n));
        for (i = 1; i < n; i++)
        {
            pop_flow[i] = pop_flow[0] + i * n;
        }
    }

    if (hdx == NULL)
    {
        snew(hdx, SQR(n));
    }

    if (ex == NULL)
    {
        snew(ex, n);
    }

    if (work == NULL)
    {
        lwork = (n > 3) ? SQR(n) : 5;
        snew(work, lwork);
    }

    if (rotr == NULL)
    {
        snew(rotr, n);
        snew(rotr[0], SQR(n));
        for (i = 1; i < n; i++)
        {
            rotr[i]  = rotr[0] + i * n;
        }
        snew(roti, n);
        snew(roti[0], SQR(n));
        for (i = 1; i < n; i++)
        {
            roti[i]  = roti[0] + i * n;
        }
        snew(rotsr, n);
        snew(rotsr[0], SQR(n));
        for (i = 1; i < n; i++)
        {
            rotsr[i] = rotsr[0] + i * n;
        }
        snew(rotsi, n);
        snew(rotsi[0], SQR(n));
        for (i = 1; i < n; i++)
        {
            rotsi[i] = rotsi[0] + i * n;
        }
        snew(rotr_new, n);
        snew(rotr_new[0], SQR(n));
        for (i = 1; i < n; i++)
        {
            rotr_new[i]  = rotr_new[0] + i * n;
        }
        snew(roti_new, n);
        snew(roti_new[0], SQR(n));
        for (i = 1; i < n; i++)
        {
            roti_new[i]  = roti_new[0] + i * n;
        }
    }

    // calculate the charges from the wave function
    // and copy as the first vector
    // NO! set to zero instead - perhaps improve efficiency...
    for (i = 0; i < ct->dim; i++)
    {
        ct->q_act[i] = ct->q_old[i] = 0.0; //alex: test for identical dimers with a large spacial separation
    }
    // ct->q_act[i] = ct->q_old[i] = SQR(ct->wf[i]);
    // NO! take them from the previous iteration

    energy = 1.e10;
    // SCC cycle
    for (step = 0;; step++)
    {

        // too many iterations? break and indicate failure
        if (step >= MAXITER_BROYDEN)
        {
            fprintf(f, "Broyden/Fermi failed to converge!");
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
        // mix all eigenvectors with appropriate coefficients

        // calculate the "Fermi energy", store the Fermi coefficients in fermi_coeff
        fermi_lower = ct->ev_adiab[0] - 0.1;
        fermi_upper = ct->ev_adiab[ct->dim - 1] + 0.1;
        while (fermi_upper - fermi_lower > FERMI_CONVERG)
        {
            fermi_energy = 0.5 * (fermi_lower + fermi_upper);
            if (fermi_coef_sum(fermi_energy, ct->dim, ct->ev_adiab, fermi_coeff, ct->fermi_kt) > 1.0)
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
/*
    for (i=0; i<ct->dim; i++){
    printf("qact: %lf\n", ct->q_act[i]);
    printf("fermicoeff: %lf\n", fermi_coeff[i]);
    for (j=0; j<ct->dim; j++){
        printf("ham %d %d : %lf\n", i,j, ham[j * ct->dim + i]);
    }
    }
 */
        // calculate the energy and check convergence
        old_energy = energy;
        energy1    = energy2 = 0.0;
        // indices i and j run over the sites
        for (i = 0; i < ct->dim; i++)
        {
            for (j = 0; j < ct->dim; j++)
            {
                // TB Hamiltonian
                // index k runs over the adiabatic states (eigenvectors)
                for (k = 0; k < ct->dim; k++)
                {
                    energy1 += fermi_coeff[k] * ham[k * ct->dim + i] * ham[k * ct->dim + j] * ct->hamiltonian[i][j];
                }
                // Hubbard / gamma terms
                energy2 += 0.5 * ct->hubbard[i][j] * ct->q_act[i] * ct->q_act[j];
            }
        }
        energy = energy1 + energy2;

        if (step > 1 && fabs(energy - old_energy) < BROYDEN_SCFTOL)
        {
            //fprintf(f, "Broyd %2d ", step);
            fprintf(f, " %2d ", step);
            break;
        }

        // mix the charges - Broyden
        broyden(step, BROYDEN_ALMIX, ct->dim, ct->q_old, ct->q_act, broyd);
        for (i = 0; i < ct->dim; i++)
        {
            ct->q_act[i] = ct->q_old[i];
        }

    } /* end SCC cycle */

    // print energies of adiabatic states - these are pure, but obtained with non-pure (Fermi-mixed) input state
    //fprintf(f, " Eigenvalues:");
    for (i = 0; i < ct->dim; i++)
    {
        fprintf(f, "%8.5f", ct->ev_adiab[i]);
    }
    //fprintf(f, "    ");

    /* surface hopping tried according to Tully JCP 93, 1061 (1990)
     * "fewest switches" algorithm
     * diabatization according to Persico JCP 114, 10608 (2001)
     * approximate decoherence correction adapted from Persico JCP 126, 134114 (2007)
     */

    /* exp_imag_matrix - has to be written yet !!! */

    /* 1. set up the state vectors */

    /*   a. back up state vectors from the previous step, tfs_vector_old*/

    if (ct->tfs_initialization_step)
    {
        /* let us start with tfs_vector_old == tfs_vector */
        for (i = 0; i < ct->dim; i++)
        {
            for (j = 0; j < ct->dim; j++)
            {
                ct->tfs_vector_old[i][j] = ham[i * ct->dim + j];
            }
        }
        for (i = 0; i < ct->dim; i++)
        {
            ct->ev_adiab_old[i] = ct->ev_adiab[i];
        }

        //ALEX: project input-WF on adiabatic states. i.e. take diabatic state as initial condiation instead of lowest adiabatic state.
        //for (i=0; i<ct->dim; i++) // for every surface
        //  dot_product = 0.;
        //  for (k=0; k<ct->dim; k++) // for every FO
        //    dot_product += ct->tfs_vector_old[i][k] * ct->wf[k];
        //  ct->tfs_popul[i]=dot_product;
    }
    else
    {
        /* push tfs_vector to tfs_vector_old */
        for (i = 0; i < ct->dim; i++)
        {
            for (j = 0; j < ct->dim; j++)
            {
                ct->tfs_vector_old[i][j] = ct->tfs_vector[i][j];
            }
        }
    }

    /*   b. back up populations (complex coefficients) of states from the previous step, tfs_popul_old */

    for (i = 0; i < ct->dim; i++)
    {
        ct->tfs_popul_old[i]             = ct->tfs_popul[i];
        ct->tfs_popul_old[ct->dim + i]   = ct->tfs_popul[ct->dim + i];
    }

    /*  c. save the new state vectors */
// IT IS NOT NECESSARY TO CHECK THE SIGN HERE. WE WILL CALCULATE OVERLAP AND THIS INCLUDES SIGN CHECK
    for (i = 0; i < ct->dim; i++)
    {
        /* calculate the dot product with tfs_vector_old[i],
           to check if we shall invert the sign
         */
        dot_product = 0.;
        for (k = 0; k < ct->dim; k++)// WHAT IF THERE IS CROSSING OF ADIABATIC STATES? WE SHOULD CHECK THAT STATE i IS SAME AS IN PREVIOUS STEP
        {
            dot_product += ct->tfs_vector_old[i][k] * ham[i * ct->dim + k];
        }
        /* for positive dot_product, take the vector as it is
         * for negative dot_product, take (-1) * the vector
         */
        for (k = 0; k < ct->dim; k++)
        {
            ct->tfs_vector[i][k] = dot_product > 0 ? ham[i * ct->dim + k] : -ham[i * ct->dim + k];
        }
        //ct->tfs_vector[i][k] = ham[i * ct->dim + k]; //we will check phase later
    }

    /* 2. calculate the overlap of state vectors
            tfs_overlap[i][j] is the matrix T[i][j] = <tfs_vector_old[i] | tfs_vector[j]>
     */

    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_overlap[i][j] = 0.;
            for (k = 0; k < ct->dim; k++)
            {
                ct->tfs_overlap[i][j] += ct->tfs_vector_old[i][k] * ct->tfs_vector[j][k];
            }
        }
    }
/*
   // If T matrix is replaced by <tfs_vector_old[i] | tfs_vector[j]> then there is an additional orthogonalization in newton-X which we are missing here (see Persico JCP 114, 10608 (2001) Appendix B)
   // However, there is also reordering of the states which will be done here
   for (i=0; i<ct->dim; i++)
    for (j=0; j<ct->dim; j++) {
      a=ct->tfs_overlap[i][i]*ct->tfs_overlap[i][i] + ct->tfs_overlap[j][j]*ct->tfs_overlap[j][j];
      b=ct->tfs_overlap[i][j]*ct->tfs_overlap[i][j] + ct->tfs_overlap[i][j]*ct->tfs_overlap[j][i];
      if ( a < b )
      for (k=0; k<ct->dim; k++){
        ski = ct->tfs_overlap[k][i];
        ct->tfs_overlap[k][i] = ct->tfs_overlap[k][j];
        ct->tfs_overlap[k][j] = ski;
      }
    }
 */
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            a = ct->tfs_overlap[i][i]*ct->tfs_overlap[i][i] + ct->tfs_overlap[j][j]*ct->tfs_overlap[j][j];
            b = ct->tfs_overlap[i][j]*ct->tfs_overlap[i][j] + ct->tfs_overlap[i][j]*ct->tfs_overlap[j][i];
            if (a < b)                           //switch occurred between state i and j. we should probably take care of this
            {
                if (ct->tfs_overlap[i][j] < 0.0) //state j has changed sign
                {
                    for (k = 0; k < ct->dim; k++)
                    {
                        ct->tfs_vector[j][k] = -ct->tfs_vector[j][k];
                    }
                }
            }
            else if (i == j && ct->tfs_overlap[i][i] < 0.0)//state stays the same but has changed sign
            {
                for (k = 0; k < ct->dim; k++)
                {
                    ct->tfs_vector[i][k] = -ct->tfs_vector[i][k];
                }
            }
        }
    }



    /* print out the overlap matrix */
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            fprintf(f2, "%9.5f", ct->tfs_overlap[i][j]);
        }
    }
    fprintf(f2, "\n");

    /* print out the state vectors */
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            fprintf(f3, "%9.5f", ct->tfs_vector[i][j]);
        }
    }
    fprintf(f3, "\n");

    /* check if the matrix is unitary - in our case orthogonal (real matrix)
       printf("Testing the orthogonality of T - this should be a unit matrix:\n");
       for (i=0; i<n; i++) {
       for (j=0; j<n; j++) {
        dot_product = 0.;
        for (k=0; k<n; k++)
          dot_product += ct->tfs_overlap[i][k] * ct->tfs_overlap[j][k];
        printf(" %9.6f", dot_product);
       }
       printf("\n");
       }
     */


    /* 3. do the propagation of the wave function / populations
            this is equivalent to the propag_ld function in NewtonX */

    /*    a. do not do it in the first step of the simulation */

    if (ct->tfs_initialization_step)
    {
        ct->tfs_initialization_step = 0;
    }
    else
    {

        /*    b. set up the diabatic hamiltonian at the end of the time step Hdia(dt) = T(dt) * E(dt) * T^t(dt) : Eq. B7 */
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                ct->per_diab_hamiltonian[i][j] = 0.;
                for (k = 0; k < n; k++)
                {
                    ct->per_diab_hamiltonian[i][j] += ct->tfs_overlap[i][k] * ct->ev_adiab[k] * ct->tfs_overlap[j][k];
                }
            }
        }

        /*    c. divide the propagation into smaller elementary steps */
        timetau = ct->rk_timestep / (double) ld_nstep;
        for (ld_istep = 0; ld_istep < ld_nstep; ld_istep++)
        {
            /* matrix Hdx                 = (Hdia(istep*timetau) + Hdia((istep+1)*tau)) / 2
               matrix Hdia(istep*timetau) = E(0) + (2*istep-1)/(2*ld_nstep) * (Hdia(dt) - E(0))
             */
            ld_fac = (double) (2. * ld_istep + 1) / (double) (2 * ld_nstep);
            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    hdx[i + j * n] = timetau * (ld_fac * ct->per_diab_hamiltonian[i][j] + (i == j ? (1. - ld_fac) * ct->ev_adiab_old[i] : 0.));
                }
            }
            /* eigenproblem of exp[-i * Hdx * timetau] */
            dsyev(n, hdx, ex, work, lwork);
            /* construct the rotation matrix (real + imaginary parts) for the elementary step */
            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    rotsr[i][j] = 0.;
                    rotsi[i][j] = 0.;
                }
            }
            for (k = 0; k < n; k++)
            {
                ld_sk = sin(ex[k]);
                ld_ck = cos(ex[k]);
                for (i = 0; i < n; i++)
                {
                    for (j = 0; j < n; j++)
                    {
                        rotsr[i][j] += hdx[i + k * n] * hdx[j + k * n] * ld_ck;
                        rotsi[i][j] -= hdx[i + k * n] * hdx[j + k * n] * ld_sk;
                    }
                }
            }
            /* accumulate the overall rotation matrix (real + imaginary) */
            if (ld_istep == 0)
            {
                for (i = 0; i < n; i++)
                {
                    for (j = 0; j < n; j++)
                    {
                        rotr[i][j] = rotsr[i][j];
                        roti[i][j] = rotsi[i][j];
                    }
                }
            }
            else
            {
                /* tw  = rotsr * rotr
                   tx  = rotsi * roti
                   ty  = rotsr * roti
                   tz  = rotsi * rotr
                   rotr = tw - tx
                   roti = ty + tz
                 */
                for (i = 0; i < n; i++)
                {
                    for (j = 0; j < n; j++)
                    {
                        rotr_new[i][j] = 0.;
                        roti_new[i][j] = 0.;
                        for (k = 0; k < n; k++)
                        {
                            rotr_new[i][j] += rotsr[i][k] * rotr[k][j] - rotsi[i][k] * roti[k][j];
                            roti_new[i][j] += rotsr[i][k] * roti[k][j] + rotsi[i][k] * rotr[k][j];
                        }
                    }
                }
                for (i = 0; i < n; i++)
                {
                    for (j = 0; j < n; j++)
                    {
                        rotr[i][j] = rotr_new[i][j];
                        roti[i][j] = roti_new[i][j];
                    }
                }
            }
            /*    d. update the coefficients in the adiabatic representation */
            for (i = 0; i < n; i++)
            {
                ct->tfs_popul[i]     = 0.;
                ct->tfs_popul[i + n] = 0.;
                for (j = 0; j < n; j++)
                {
                    for (k = 0; k < n; k++)
                    {
                        ct->tfs_popul[i]     += ct->tfs_overlap[j][i] *
                            (ct->tfs_popul_old[k] * rotr[j][k] - ct->tfs_popul_old[k + n] * roti[j][k]);
                        ct->tfs_popul[i + n] += ct->tfs_overlap[j][i] *
                            (ct->tfs_popul_old[k] * roti[j][k] + ct->tfs_popul_old[k + n] * rotr[j][k]);
                    }
                }
            }
        }

        /* print out new populations */
        for (k = 0; k < ct->dim; k++)
        {
            fprintf(f, "%9.6f", SQR(ct->tfs_popul[k]) + SQR(ct->tfs_popul[ct->dim + k]));
        }

        /* we are on ct->surface at the moment
         * now, go over all of the other states,
         * and calculate switching probability */


        /* has a switch of states just occured? */ // I guess this shouldn't happen after I corrected ordering before diabatic hamiltonian is constructed
        surf_old = surf;
        if (fabs(ct->tfs_overlap[surf][surf]) < eps_switch)
        {
            printf("state switching has occurred! Overlap is only: %f  \n", ct->tfs_overlap[surf][surf]);
            for (j = 0; j < n; j++)
            {
                //ct->surf_prob[j] = 0.; why was this here?
                maxval = 0.;
            }
            for (j = 0; j < n; j++)
            {
                if (fabs((float) ct->tfs_overlap[j][surf_old]) > maxval)
                {
                    surf        = j;
                    ct->surface = j; //I guess this was missing
                    maxval      = fabs((float) ct->tfs_overlap[j][surf_old]);
                }
            }
            printf("state %d probably became state %d! Overlap: %f  \n", surf_old, surf, ct->tfs_overlap[surf][surf_old]);
        }



        /* if population of occupied surface increases,
         * then no hopping!
         */
        ac0[0] = acoer0[surf_old];
        ac0[1] = acoei0[surf_old];
        ac1[0] = acoer1[surf];
        ac1[1] = acoei1[surf];
        popk0  = SQR(ac0[0]) + SQR(ac0[1]);
        popk1  = SQR(ac1[0]) + SQR(ac1[1]);
        if (popk1 > popk0)
        {
            for (j = 0; j < n; j++)
            {
                ct->surf_prob[j] = -1.;
            }

            /* otherwise, calculate the probabilities */
        }
        else
        {
            /* calculate the gross flow */
            gross_flow = 0.;
            for (j = 0; j < n; j++)
            {
                diff_flow = SQR(acoer1[j]) + SQR(acoei1[j]) - (SQR(acoer0[j]) + SQR(acoei0[j]));
                if (diff_flow > 0.)
                {
                    gross_flow += diff_flow;
                }
            }
            /* go over the states */
            for (j = 0; j < n; j++)
            {
                /* no hopping from state 'surf' to itself, obviously */
                if (j == surf)
                {
                    ct->surf_prob[j] = -1.;
                    /* calculate the probability for all of the states except 'surf' */
                }
                else
                {
                    diff_flow = SQR(acoer1[j]) + SQR(acoei1[j]) - (SQR(acoer0[j]) + SQR(acoei0[j]));
                    /* no hopping to a state with a decreased population */
                    if (diff_flow < 0.)
                    {
                        ct->surf_prob[j] = -0.1;
                        /* otherwise, calculate probability according to Prezhdo */
                    }
                    else
                    {
                        ct->surf_prob[j] = diff_flow * (popk0 - popk1) / (popk0 * gross_flow);
                    }
                }
            }
        }


        tot_prob = 0.;
        for (j = 0; j < n; j++)
        {
            tot_prob += (ct->surf_prob[j] > 0. ? ct->surf_prob[j] : 0.);
        }
        if (tot_prob > 1.)
        {
            for (j = 0; j < n; j++)
            {
                if (ct->surf_prob[j] > 0.)
                {
                    ct->surf_prob[j] /= tot_prob;
                }
            }
        }

        /* END OF PROBABILITIES! */

        /* write out the probabilities */
        for (k = 0; k < ct->dim; k++)
        {
            fprintf(f, "%10.6f", ct->surf_prob[k]);
        }

        /* generate a random number */
        random_number         = rand()/(double)RAND_MAX; //hopping switched off for testing
        cumulated_probability = 0.;

        /* and determine the surface to continue on */
        for (j = 0; j < ct->dim; j++)
        {
            if (j != ct->surface && ct->surf_prob[j] > 0.)
            {
                cumulated_probability += ct->surf_prob[j];
                if (cumulated_probability > random_number)
                {
                    ct->surface = j;
                    break;
                }
            }
        }

        /* decoherence correction
         * proposed by Truhlar JCTC 1, 527 (2005)
         * recommended by Persico 126, 134114 (2007)
         * simplified to:
         * tau_km = hbar / |E_k - E_m| * TFS_DECAY_CONSTANT
         * TFS_DECAY_CONSTANT chosen to be 3. (should be quite insensitive, though)
         */

        if (ct->decoherence)
        {
            popul_norm = 0.;
            for (j = 0; j < n; j++)
            {
                if (j != ct->surface)
                {
                    decay_time                    = TFS_DECAY_CONSTANT / fabs(ct->ev_adiab[j] - ct->ev_adiab[ct->surface]);
                    ct->tfs_popul[j]             *= exp( -ct->rk_timestep / decay_time);
                    ct->tfs_popul[ct->dim + j]   *= exp( -ct->rk_timestep / decay_time);
                    popul_norm                   += SQR(ct->tfs_popul[j]) + SQR(ct->tfs_popul[n + j]);
                }
            }
            current_surface_factor                  = sqrt((1. - popul_norm) / (SQR(ct->tfs_popul[ct->surface]) + SQR(ct->tfs_popul[n + ct->surface])));
            ct->tfs_popul[ct->surface]             *= current_surface_factor;
            ct->tfs_popul[ct->dim + ct->surface]   *= current_surface_factor;
        }

    } /* end if initialization step else */

    /* assign the new wave function to the eigenvector ct->surface */
    for (i = 0; i < ct->dim; i++)
    {
        ct->wf[i] = ct->tfs_vector[ct->surface][i];
    }

    /* print out current state */
    fprintf(f, " %d", ct->surface);
    for (i = 0; i < ct->dim; i++)
    {
        fprintf(f, "%8.4f", SQR(ct->wf[i]));
    }
    fprintf(f, "\n");

    /* push the current ev_adiab to ev_adiab_old, to be used in the next step */
    for (i = 0; i < ct->dim; i++)
    {
        ct->ev_adiab_old[i] = ct->ev_adiab[i];
    }

    return 0;
}

