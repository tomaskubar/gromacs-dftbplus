/*************************
 * Charge transfer in DNA
 * Tomas Kubar
 * Neural Nets: Mila Kraemer
 *************************/

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "mdatoms.h"

#include "transfer.h"
#include <math.h>

#ifdef GMX_MPI
#define PRINTF(...) if (ct_mpi_rank == 0) printf(__VA_ARGS__)
#else
#define PRINTF(...) printf(__VA_ARGS__)
#endif
#define CUB(x) ((x)*(x)*(x))


/****************************
 * INITIALIZE THE DFTB CODE *
 ****************************/

#ifdef GMX_MPI
void init_dftb(t_mdatoms *mdatoms, dftb_t *dftb, charge_transfer_t *ct, char *slko_path, int ct_mpi_rank)
#else
void init_dftb(t_mdatoms *mdatoms, dftb_t *dftb, charge_transfer_t *ct, char *slko_path)
#endif
{
/* Initializes the dftb calculations. Reads information from Slater-Koster files, copies information that was gathered in init_charge_transfer(), allocates memory */

// PARAMETERS:
// mdatoms     = (in) variable of GROMACS data-structure-type (we need atomic masses)
// dftb        = (out) main data structure for DFTB calculations
// ct          = (in) main data structure with information about the charge transfer calculation
// slko_path   = (in) path to the Slater-Koster files
// ct_mpi_rank = (in) the rank of the process in MPI calculations (each process in MPI calculations can be destinguished by its value of ct_mpi_rank)
//////////////////////////////////////////////////////////////////////////////////////////
    const char *type_symbols = "chnosfklbpi";        // fe is z 16.06.2016 because of one letter relation between element and code
                                                  // K is Br, L is Cl 13Sep2018 Daniel for halogens
    const char *suffix1      = "-c.spl";
    const char *suffix2      = "-uncomp-c.spl"; // now diffrent folder for each parameter set with diefferent r_dens and r_wf

    char        filename[128];
    int         i, j, k, l, izpj, counter;
    double      espin, qzeroh[3], uhubbh[3], mass;
    FILE       *f;

    char       *line = NULL;
    char      * pch;
    size_t      len;

#ifdef GMX_MPI
    printf("DFTB initialization at rank %d\n", ct_mpi_rank);
#else
    printf("DFTB initialization\n");
#endif


    for (i = 0; i < DFTB_MAXTYPES; i++)
    {
        dftb->qzero1[i] = 0.0;
        dftb->uhubb1[i] = 0.0;
        dftb->qzero2[i] = 0.0;
        dftb->uhubb2[i] = 0.0;
        for (j = 0; j < DFTB_MAXTYPES; j++)
        {
            PRINTF("Atomtype pair %d-%d\n", i+1, j+1);
            /* read the tables for DFTB phase 1 - calculation of monomers */
            sprintf(filename, "%s%c%c%s", slko_path, type_symbols[i], type_symbols[j], suffix1);
            f = fopen(filename, "r");
            if (f == NULL)
            {
                PRINTF("Cannot open the parameter file %s, exiting!\n", filename);
                exit(-1);
            }
            //printf("fscanf(f, \"%%lf %%d\", &(dftb->dr1[i][j]), &(dftb->dim1[i][j]))\n");
            //printf("%lf %d\n", dftb->dr1[i][j], dftb->dim1[i][j]);
            if (i == j)
            {
                fscanf(f, "%lf %d %d", &(dftb->dr1[i][j]), &(dftb->dim1[i][j]), &(dftb->lmax[i]));
                //printf("&(dftb->skself1[i][0]), &(dftb->skself1[i][1]), &(dftb->skself1[i][2]), &espin\n");
                fscanf(f, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                       &(dftb->skself1[i][0]), &(dftb->skself1[i][1]), &(dftb->skself1[i][2]), &espin,
                       uhubbh+2, uhubbh+1, uhubbh, qzeroh+2, qzeroh+1, qzeroh);
                //printf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", dftb->skself1[i][0], dftb->skself1[i][1], dftb->skself1[i][2], espin,
                //          uhubbh[2], uhubbh[1], uhubbh[0], qzeroh[2], qzeroh[1], qzeroh[0]);
                dftb->uhubb1[i] = uhubbh[0];
                for (k = 0; k < 3; k++)
                {
                    dftb->qzero1[i] += qzeroh[k];
                }
            }
            else
            {
                fscanf(f, "%lf %d", &(dftb->dr1[i][j]), &(dftb->dim1[i][j]));
            }

            snew(dftb->skhtab1[i][j], dftb->dim1[i][j]);
            snew(dftb->skstab1[i][j], dftb->dim1[i][j]);
            for (k = 0; k < dftb->dim1[i][j]; k++)
            {
                //   printf("for (l=0; l<10; l++) fscanf(f, \"%%lf\", &(dftb->skhtab1[i][j][k][l]))\n");
                for (l = 0; l < 10; l++)
                {
                    fscanf(f, "%lf", &(dftb->skhtab1[i][j][k][l]));
                    // printf("%lf ", dftb->skhtab1[i][j][k][l]);
                }
                //printf("for (l=0; l<10; l++) fscanf(f, \"%%lf\", &(dftb->skstab1[i][j][k][l]))\n");
                for (l = 0; l < 10; l++)
                {
                    fscanf(f, "%lf", &(dftb->skstab1[i][j][k][l]));
                    //printf("%lf ", dftb->skstab1[i][j][k][l]);
                }
            }
            PRINTF("skfile for pair %d-%d, phase 1: %s\n", i+1, j+1, filename);
            fclose(f);

            /* read the tables for DFTB phase 2 - calculation of complex */
            sprintf(filename, "%s%c%c%s", slko_path, type_symbols[i], type_symbols[j], suffix2);
            f = fopen(filename, "r");
            if (f == NULL)
            {
                PRINTF("Cannot open the parameter file %s, exiting!\n", filename);
                exit(-1);
            }

            //in converntional SLKOs there are 3 entries for i=j, in our CT-SLKO lmax is missing. In order to be able to use both formats in phase2, fscan was replaced by getline
            //fscanf(f, "%lf %d", &(dftb->dr2[i][j]), &(dftb->dim2[i][j]));
            getline(&line, &len, f);
            pch              = strtok (line, " ");
            dftb->dr2[i][j]  = atof(pch);
            pch              = strtok (NULL, " ");
            dftb->dim2[i][j] = atoi(pch);

            //printf("%lf %d\n", dftb->dr2[i][j], dftb->dim2[i][j]);
            if (i == j)
            {
                fscanf(f, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                       &(dftb->skself2[i][0]), &(dftb->skself2[i][1]), &(dftb->skself2[i][2]), &espin,
                       uhubbh+2, uhubbh+1, uhubbh, qzeroh+2, qzeroh+1, qzeroh);
                //printf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", dftb->skself2[i][0], dftb->skself2[i][1], dftb->skself2[i][2], espin,
                //          uhubbh[2], uhubbh[1], uhubbh[0], qzeroh[2], qzeroh[1], qzeroh[0]);
                dftb->uhubb2[i] = uhubbh[0];
                for (k = 0; k < 3; k++)
                {
                    dftb->qzero2[i] += qzeroh[k];
                }
            }
            snew(dftb->skhtab2[i][j], dftb->dim2[i][j]);
            snew(dftb->skstab2[i][j], dftb->dim2[i][j]);
            for (k = 0; k < dftb->dim2[i][j]; k++)
            {
                for (l = 0; l < 10; l++)
                {
                    fscanf(f, "%lf", &(dftb->skhtab2[i][j][k][l]));
                    //printf("%lf ", dftb->skhtab2[i][j][k][l]);
                }
                for (l = 0; l < 10; l++)
                {
                    fscanf(f, "%lf", &(dftb->skstab2[i][j][k][l]));
                    //printf("%lf ", dftb->skstab2[i][j][k][l]);
                }
                //printf("\n");
            }
            PRINTF("skfile for pair %d-%d, phase 2: %s\n", i+1, j+1, filename);
            fclose(f);
        }
    }

    /* deal with phase1 and certain parts of phase2 */
    snew(dftb->phase1, ct->sites);
    dftb->phase2.nn   = 0;
    dftb->phase2.ne   = ct->extcharges_cplx;
    dftb->phase2.nel  = 0;
    dftb->phase2.norb = 0;

    /* prepare the broyden structures */
    snew(dftb->broyden, ct->sites);
    for (i = 0; i < ct->sites; i++)
    {
        snew(dftb->broyden[i].f, ct->site[i].atoms);
        snew(dftb->broyden[i].ui, ct->site[i].atoms);
        snew(dftb->broyden[i].vti, ct->site[i].atoms);
        snew(dftb->broyden[i].t1, ct->site[i].atoms);
        snew(dftb->broyden[i].dumvi, ct->site[i].atoms);
        snew(dftb->broyden[i].df, ct->site[i].atoms);
        snew(dftb->broyden[i].vector, ct->site[i].atoms);
        snew(dftb->broyden[i].unit31, ct->site[i].atoms);
        snew(dftb->broyden[i].unit32, ct->site[i].atoms);
    }
    /* determine the numbers of electrons and atom types (in arrays) */
    dftb->phase2.nel = 0;
    for (i = 0; i < ct->sites; i++)
    {
        dftb->phase1[i].nel  = ct->site[i].nel;
        dftb->phase1[i].norb = ct->site[i].norb;
        dftb->phase1[i].nn   = ct->site[i].atoms;
        dftb->phase2.nn     += ct->site[i].atoms;
        dftb->phase2.nel    += dftb->phase1[i].nel;
        dftb->phase2.norb   += dftb->phase1[i].norb;

        /* double arrays */
        snew(dftb->phase1[i].x, ct->site[i].atoms);
        snew(dftb->phase1[i].x_opt, ct->site[i].atoms);
        snew(dftb->phase1[i].xe, ct->site[i].extcharges);
        snew(dftb->phase1[i].mass, ct->site[i].atoms);
        snew(dftb->phase1[i].ze, ct->site[i].extcharges);
        snew(dftb->phase1[i].qmat, ct->site[i].atoms);
        snew(dftb->phase1[i].qmold, ct->site[i].atoms);
        snew(dftb->phase1[i].qmulli, dftb->phase1[i].norb);
        snew(dftb->phase1[i].ev, dftb->phase1[i].norb);
        snew(dftb->phase1[i].occ, dftb->phase1[i].norb);
        snew(dftb->phase1[i].a, dftb->phase1[i].norb);
        snew(dftb->phase1[i].a[0], SQR(dftb->phase1[i].norb));
        for (j = 1; j < dftb->phase1[i].norb; j++)
        {
            dftb->phase1[i].a[j] = dftb->phase1[i].a[0] + j * dftb->phase1[i].norb;
        }

        snew(dftb->phase1[i].a_old, dftb->phase1[i].norb);
        snew(dftb->phase1[i].a_old[0], SQR(dftb->phase1[i].norb));
        for (j = 1; j < dftb->phase1[i].norb; j++)
        {
            dftb->phase1[i].a_old[j] = dftb->phase1[i].a_old[0] + j * dftb->phase1[i].norb;
        }
        snew(dftb->phase1[i].a_ref, dftb->phase1[i].norb);
        snew(dftb->phase1[i].a_ref[0], SQR(dftb->phase1[i].norb));
        for (j = 1; j < dftb->phase1[i].norb; j++)
        {
            dftb->phase1[i].a_ref[j] = dftb->phase1[i].a_ref[0] + j * dftb->phase1[i].norb;
        }

        snew(dftb->phase1[i].b, dftb->phase1[i].norb);
        snew(dftb->phase1[i].b[0], SQR(dftb->phase1[i].norb));
        for (j = 1; j < dftb->phase1[i].norb; j++)
        {
            dftb->phase1[i].b[j] = dftb->phase1[i].b[0] + j * dftb->phase1[i].norb;
        }
        snew(dftb->phase1[i].a_trans, dftb->phase1[i].norb * dftb->phase1[i].norb);
        snew(dftb->phase1[i].b_trans, dftb->phase1[i].norb * dftb->phase1[i].norb);
        snew(dftb->phase1[i].hamil, dftb->phase1[i].norb);
        snew(dftb->phase1[i].hamil[0], SQR(dftb->phase1[i].norb));
        for (j = 1; j < dftb->phase1[i].norb; j++)
        {
            dftb->phase1[i].hamil[j] = dftb->phase1[i].hamil[0] + j * dftb->phase1[i].norb;
        }
        snew(dftb->phase1[i].overl, dftb->phase1[i].norb);
        snew(dftb->phase1[i].overl[0], SQR(dftb->phase1[i].norb));
        for (j = 1; j < dftb->phase1[i].norb; j++)
        {
            dftb->phase1[i].overl[j] = dftb->phase1[i].overl[0] + j * dftb->phase1[i].norb;
        }

        //snew(dftb->phase1[i].overl_diffuse, dftb->phase1[i].norb);
        //  snew(dftb->phase1[i].overl_diffuse[0], SQR(dftb->phase1[i].norb));
        //  for(j = 1; j < dftb->phase1[i].norb; j++)
        //    dftb->phase1[i].overl_diffuse[j] = dftb->phase1[i].overl_diffuse[0] + j * dftb->phase1[i].norb;
        //dftb->phase1[i].gammamat = (double **) malloc(ct->atoms[i] * sizeof(double *));
        //  dftb->phase1[i].gammamat[0] = (double *) malloc(ct->atoms[i] * ct->atoms[i] * sizeof(double));
        snew(dftb->phase1[i].gammamat, ct->site[i].atoms);
        snew(dftb->phase1[i].gammamat[0], SQR(ct->site[i].atoms));
        for (j = 1; j < ct->site[i].atoms; j++)
        {
            dftb->phase1[i].gammamat[j] = dftb->phase1[i].gammamat[0] + j * ct->site[i].atoms;
        }
        snew(dftb->phase1[i].shift, ct->site[i].atoms);
        snew(dftb->phase1[i].dshift, ct->site[i].atoms);
        snew(dftb->phase1[i].shiftE, ct->site[i].atoms);
        snew(dftb->phase1[i].aux, 3 * dftb->phase1[i].norb);

        if (ct->do_lambda_i > 1 || ct->jobtype == cteJFSSH || ct->jobtype == cteBCJFSSH || ct->jobtype == cteNOMOVEMENT || ct->hamiltonian_type)
        {
            snew(dftb->phase1[i].grad, dftb->phase1[i].nn);
            snew(dftb->phase1[i].partgrad, dftb->phase1[i].nn);
        }
        if (ct->qmmm == 3)
        {
            snew(dftb->phase1[i].neighbors_pme, ct->site[i].atoms);
            snew(dftb->phase1[i].neighbor_pme, ct->site[i].atoms);
        }

        /* int arrays */
        snew(dftb->phase1[i].izp, ct->site[i].atoms);
        snew(dftb->phase1[i].ind, ct->site[i].atoms + 1);
        dftb->phase1[i].ind[0] = 0;
        for (j = 0; j < ct->site[i].atoms; j++)
        {
            dftb->phase1[i].izp[j]   = ct->site[i].atomtype[j];
            izpj                     = dftb->phase1[i].izp[j];
            dftb->phase1[i].ind[j+1] = dftb->phase1[i].ind[j] + dftb->lmax[izpj]* dftb->lmax[izpj];
        }
        dftb->phase1[i].ndim = dftb->phase1[i].norb;
        dftb->phase1[i].ne   = ct->site[i].extcharges;

        /* atom masses */
        mass = 0.0;
        for (j = 0; j < ct->site[i].atoms; j++)
        {
            dftb->phase1[i].mass[j] = mdatoms->massT[ct->site[i].atom[j]];
            mass                   += dftb->phase1[i].mass[j];
            printf("site %d  mass %d   %f\n", i, j,  mdatoms->massT[ct->site[i].atom[j]]);
        }
        dftb->phase1[i].inv_tot_mass = 1.0 / mass;
    }

    /* the remainder of phase2 */
    snew(dftb->phase2.x, dftb->phase2.nn);
    snew(dftb->phase2.xe, dftb->phase2.ne);
    snew(dftb->phase2.mass, dftb->phase2.nn);
    snew(dftb->phase2.ze, dftb->phase2.ne);
    snew(dftb->phase2.qmat, dftb->phase2.nn);
    snew(dftb->phase2.shift, dftb->phase2.nn);
    snew(dftb->phase2.shiftE, dftb->phase2.nn);
    snew(dftb->phase2.ind, dftb->phase2.nn + 1);
    snew(dftb->phase2.atind, dftb->phase2.nn + 1);
    snew(dftb->phase2.inf, ct->sites + 1);
    snew(dftb->phase2.ihomo, ct->sites + 1);
    snew(dftb->phase2.izp, dftb->phase2.nn);
    snew(dftb->phase2.hamil, dftb->phase2.norb);
    snew(dftb->phase2.hamil[0], SQR(dftb->phase2.norb));
    for (j = 1; j < dftb->phase2.norb; j++)
    {
        dftb->phase2.hamil[j] = dftb->phase2.hamil[0] + j * dftb->phase2.norb;
    }
    snew(dftb->phase2.overl, dftb->phase2.norb);
    snew(dftb->phase2.overl[0], SQR(dftb->phase2.norb));
    for (j = 1; j < dftb->phase2.norb; j++)
    {
        dftb->phase2.overl[j] = dftb->phase2.overl[0] + j * dftb->phase2.norb;
    }
    snew(dftb->phase2.Taf, dftb->phase2.norb);
    snew(dftb->phase2.Taf[0], SQR(dftb->phase2.norb));
    for (j = 1; j < dftb->phase2.norb; j++)
    {
        dftb->phase2.Taf[j] = dftb->phase2.Taf[0] + j * dftb->phase2.norb;
    }
/*
   snew(dftb->phase2.THamil, dftb->phase2.norb);
    snew(dftb->phase2.THamil[0], SQR(dftb->phase2.norb));
    for(j = 1; j < dftb->phase2.norb; j++)
      dftb->phase2.THamil[j] = dftb->phase2.THamil[0] + j * dftb->phase2.norb;
    snew(dftb->phase2.OverlF, dftb->phase2.norb);
    snew(dftb->phase2.OverlF[0], SQR(dftb->phase2.norb));
    for (j = 1; j < dftb->phase2.norb; j++)
    {
        dftb->phase2.OverlF[j] = dftb->phase2.OverlF[0] + j * dftb->phase2.norb;
    }

 */
    if (ct->do_lambda_i == 3 || ct->jobtype == cteJFSSH || ct->jobtype == cteBCJFSSH || ct->jobtype == cteNOMOVEMENT)
    {
        snew(dftb->phase2.b, dftb->phase2.norb);
        snew(dftb->phase2.b[0], SQR(dftb->phase2.norb));
        for (j = 1; j < dftb->phase2.norb; j++)
        {
            dftb->phase2.b[j] = dftb->phase2.b[0] + j * dftb->phase2.norb;
        }

        snew(dftb->phase2.grad, dftb->phase2.nn);
        snew(dftb->phase2.partgrad, dftb->phase2.nn);

    }
    //snew(dftb->phase2.overl_hybrid, dftb->phase2.norb);
    //  snew(dftb->phase2.overl_hybrid[0], SQR(dftb->phase2.norb));
    //  for(j = 1; j < dftb->phase2.norb; j++)
    //    dftb->phase2.overl_hybrid[j] = dftb->phase2.overl_hybrid[0] + j * dftb->phase2.norb;
    //dftb->phase2.THamilOrtho = (double **) malloc(dftb->phase2.norb * sizeof(double *));
    //  dftb->phase2.THamilOrtho[0] = (double *) malloc(SQR(dftb->phase2.norb) * sizeof(double));

/*
   snew(dftb->phase2.THamilOrtho, dftb->phase2.norb);
    snew(dftb->phase2.THamilOrtho[0], SQR(dftb->phase2.norb));
    for(j = 1; j < dftb->phase2.norb; j++)
      dftb->phase2.THamilOrtho[j] = dftb->phase2.THamilOrtho[0] + j * dftb->phase2.norb;
   snew(dftb->phase2.tij, dftb->phase2.norb);
    snew(dftb->phase2.tij[0], SQR(dftb->phase2.norb));
    for(j = 1; j < dftb->phase2.norb; j++)
      dftb->phase2.tij[j] = dftb->phase2.tij[0] + j * dftb->phase2.norb;
   snew(dftb->phase2.sij, dftb->phase2.norb);
    snew(dftb->phase2.sij[0], SQR(dftb->phase2.norb));
    for(j = 1; j < dftb->phase2.norb; j++)
      dftb->phase2.sij[j] = dftb->phase2.sij[0] + j * dftb->phase2.norb;
 */
//*
    snew(dftb->phase2.THamilOrtho, ct->dim);
    snew(dftb->phase2.THamilOrtho[0], SQR(ct->dim));
    for (j = 1; j < ct->dim; j++)
    {
        dftb->phase2.THamilOrtho[j] = dftb->phase2.THamilOrtho[0] + j * ct->dim;
    }
    snew(dftb->phase2.tij, ct->dim);
    snew(dftb->phase2.tij[0], SQR(ct->dim));
    for (j = 1; j < ct->dim; j++)
    {
        dftb->phase2.tij[j] = dftb->phase2.tij[0] + j * ct->dim;
    }
    snew(dftb->phase2.sij, ct->dim);
    snew(dftb->phase2.sij[0], SQR(ct->dim));
    for (j = 1; j < ct->dim; j++)
    {
        dftb->phase2.sij[j] = dftb->phase2.sij[0] + j * ct->dim;
    }
//*/
    snew(dftb->phase2.gammamat, ct->atoms_cplx);
    snew(dftb->phase2.gammamat[0], SQR(ct->atoms_cplx));
    for (j = 1; j < ct->atoms_cplx; j++)
    {
        dftb->phase2.gammamat[j] = dftb->phase2.gammamat[0] + j * ct->atoms_cplx;
    }
    snew(dftb->nl, ct->atoms_cplx);
    snew(dftb->nl[0], SQR(ct->atoms_cplx));
    for (j = 1; j < ct->atoms_cplx; j++)
    {
        dftb->nl[j] = dftb->nl[0] + j * ct->atoms_cplx;
    }
    snew(dftb->phase2.ev, dftb->phase2.norb);
    snew(dftb->phase2.occ, dftb->phase2.norb);
    snew(dftb->phase2.aux, 3 * dftb->phase2.norb);

    if (ct->qmmm == 3)
    {
        snew(dftb->phase2.neighbors_pme, ct->atoms_cplx);
        snew(dftb->phase2.neighbor_pme, ct->atoms_cplx);
    }

    dftb->phase2.inf[0]   = 0; // where do the orbitals of fragment i start?
    dftb->phase2.ihomo[0] = 0; //index of first homo of site i
    for (i = 0; i < ct->sites; i++)
    {
        dftb->phase2.inf[i+1]   = dftb->phase2.inf[i] + dftb->phase1[i].norb;
        dftb->phase2.ihomo[i+1] = dftb->phase2.ihomo[i] + ct->site[i].homos;
    }

    for (i = 0; i < ct->sites; i++)
    {
        for (j = 0; j < ct->site[i].atoms; j++)
        {
            dftb->phase1[i].qmat[j] = dftb->qzero1[dftb->phase1[i].izp[j]];
        }
    }
    counter               = 0;
    dftb->phase2.ind[0]   = 0;
    dftb->phase2.atind[0] = 0;
    mass                  = 0.0;
    for (i = 0; i < ct->sites; i++)
    {
        dftb->phase2.atind[i+1] = dftb->phase2.atind[i] + ct->site[i].atoms;
        for (j = 0; j < ct->site[i].atoms; j++)
        {
            dftb->phase2.izp[counter] = dftb->phase1[i].izp[j];
            izpj = dftb->phase2.izp[counter];
            dftb->phase2.ind[counter+1] = dftb->phase2.ind[counter] + dftb->lmax[izpj]* dftb->lmax[izpj];
            //dftb->phase2.mass[counter] = mdatoms->massT[ct->site[i].atom[j]]; /* mass */
            dftb->phase2.mass[counter] = dftb->phase1[i].mass[j]; /* mass */
            mass += dftb->phase2.mass[counter];
            counter++;
        }
    }
    dftb->phase2.inv_tot_mass = 1.0 / mass;
    if (counter != dftb->phase2.nn)
    {
        PRINTF("The number of atoms does not match: counter = %d, dftb->phase2.nn = %d !\n", counter, dftb->phase2.nn);
        exit(-1);
    }

    /* auxiliary arrays for function orthogonalize() */
    snew(dftb->orthogo.tij, ct->dim*ct->dim);
    snew(dftb->orthogo.sij, ct->dim*ct->dim);
    snew(dftb->orthogo.evec, ct->dim*ct->dim);
    dftb->orthogo.lwork = 26*ct->dim;
    snew(dftb->orthogo.work, 26*ct->dim*26*ct->dim);
    dftb->orthogo.liwork = 10*ct->dim;
    snew(dftb->orthogo.iwork, 10*ct->dim*10*ct->dim);
    snew(dftb->orthogo.eval, ct->dim);
    snew(dftb->orthogo.issupz, 2*ct->dim);
    /* also includes arrays for time derivative of FMOs - Weiwei X Mar 2019 */
    //  [naos][ct->dim]
    snew(dftb->orthogo.sij_old, ct->dim*ct->dim);

    snew(dftb->orthogo.fmo_old, dftb->phase2.norb);
    snew(dftb->orthogo.fmo_old[0], ct->dim * dftb->phase2.norb);
    for (i=1; i<dftb->phase2.norb; i++)
      dftb->orthogo.fmo_old[i] = dftb->orthogo.fmo_old[0] + i * ct->dim;

    snew(dftb->orthogo.fmo, dftb->phase2.norb);
    snew(dftb->orthogo.fmo[0], ct->dim * dftb->phase2.norb);
    for (i=1; i<dftb->phase2.norb; i++)
      dftb->orthogo.fmo[i] = dftb->orthogo.fmo[0] + i * ct->dim;
    // WRONG? for (i=1; i<ct->dim; i++)
    // WRONG?   dftb->orthogo.fmo[i] = dftb->orthogo.fmo[0] + i * dftb->phase2.norb;
    //  [naos][naos]
/*
    snew(dftb->orthogo.sao_old, dftb->phase2.norb);
    snew(dftb->orthogo.sao_old[0], SQR(dftb->phase2.norb));
    for (i=1; i<dftb->phase2.norb; i++)
      dftb->orthogo.sao_old[i] = dftb->orthogo.sao_old[0] + i * dftb->phase2.norb;

    snew(dftb->orthogo.sao, dftb->phase2.norb);
    snew(dftb->orthogo.sao[0], SQR(dftb->phase2.norb));
    for (i=1; i<dftb->phase2.norb; i++)
      dftb->orthogo.sao[i] = dftb->orthogo.sao[0] + i * dftb->phase2.norb;
*/
    /* arrays for state following */
    /*
       snew(dftb->orthogo.evec_ao, dftb->phase2.norb);
       snew(dftb->orthogo.evec_ao[0], dftb->phase2.norb * ct->dim);
       for(j = 1; j < dftb->phase2.norb; j++)
        dftb->orthogo.evec_ao[j] = dftb->orthogo.evec_ao[0] + j * ct->dim;
       snew(dftb->orthogo.evec_ao_old, dftb->phase2.norb);
       snew(dftb->orthogo.evec_ao_old[0], dftb->phase2.norb * ct->dim);
       for(j = 1; j < dftb->phase2.norb; j++)
        dftb->orthogo.evec_ao_old[j] = dftb->orthogo.evec_ao_old[0] + j * ct->dim;
       snew(dftb->orthogo.evec_ao_ref, dftb->phase2.norb);
       snew(dftb->orthogo.evec_ao_ref[0], dftb->phase2.norb * ct->dim);
       for(j = 1; j < dftb->phase2.norb; j++)
        dftb->orthogo.evec_ao_ref[j] = dftb->orthogo.evec_ao_ref[0] + j * ct->dim;
     */
    snew(dftb->orthogo.overlap, ct->dim);
    snew(dftb->orthogo.overlap[0], SQR(ct->dim));
    for (j = 1; j < ct->dim; j++)
    {
        dftb->orthogo.overlap[j] = dftb->orthogo.overlap[0] + j * ct->dim;
    }
    snew(dftb->orthogo.overlap_ref, ct->dim);
    snew(dftb->orthogo.overlap_ref[0], SQR(ct->dim));
    for (j = 1; j < ct->dim; j++)
    {
        dftb->orthogo.overlap_ref[j] = dftb->orthogo.overlap_ref[0] + j * ct->dim;
    }

    snew(ct->wf_old, 2*ct->dim);

    snew(dftb->overl_test, ct->dim);
    snew(dftb->overl_test[0], SQR(ct->dim));
    for (j = 1; j < ct->dim; j++)
    {
        dftb->overl_test[j] = dftb->overl_test[0] + j * ct->dim;
    }



    /* machine accuracy */
    dftb->racc = 1.0;
    while ((1.0 + dftb->racc) > 1.0)
    {
        dftb->racc /= 2.0;
    }
    dftb->racc *= 2.0;
    dftb->dacc  = 4 * dftb->racc;

#ifdef GMX_MPI
    printf("Completed DFTB initialization at rank %d\n", ct_mpi_rank);
#else
    printf("Completed DFTB initialization\n");
#endif

    return;
}

#ifdef GMX_MPI
void init_dftb_pme(dftb_t *dftb, charge_transfer_t *ct, t_inputrec *ir, int ct_mpi_rank)
#else
void init_dftb_pme(dftb_t *dftb, charge_transfer_t *ct, t_inputrec *ir)
#endif
{
    int i, j, status;
    dftb->ewaldcoeff_pme = calc_ewaldcoeff(ir->rcoulomb, ir->ewald_rtol);
    dftb->rcoulomb_pme   = ir->rcoulomb;
    dftb->nstlist_pme    = ir->nstlist;
    dftb->lastlist_pme   = -1000000; /* to have the list generated before the first PME calculation (in any case) */
    /* for phase 1 -- loop over all sites
     * POSSIBLE MODIFICATION FOR PARALLEL RUNS
     *   -- DO IT ONLY FOR THOSE SITES THAT ARE TO BE CALCULATED ON THIS RESPECTIVE CT_MPI_RANK!
     */
    for (i = 0; i < ct->sites; i++)
    {
        snew(dftb->phase1[i].q_pme, ct->site[i].atoms + ct->site[i].extcharges);
        snew(dftb->phase1[i].x_pme, ct->site[i].atoms + ct->site[i].extcharges);
        snew(dftb->phase1[i].pot, ct->site[i].atoms + ct->site[i].extcharges);
        snew(dftb->phase1[i].pot2, ct->site[i].atoms);
        snew(dftb->phase1[i].pot3, ct->site[i].atoms);
        snew(dftb->phase1[i].pot4, ct->site[i].atoms);
        snew(dftb->phase1[i].nrnb_pme, 1);
        for (j = 0; j < ct->site[i].atoms; j++)
        {
            dftb->phase1[i].qmat[j] = dftb->qzero1[dftb->phase1[i].izp[j]];
        }
        snew(dftb->phase1[i].pmedata, 1);
        status = gmx_pme_init_dftb(dftb->phase1[i].pmedata, ir, ct->site[i].atoms + ct->site[i].extcharges);
    }
    /* for phase 2 -- POSSIBLY ONLY ON CT_MPI_RANK==0
     * ! */
    snew(dftb->phase2.q_pme, ct->atoms_cplx + ct->extcharges_cplx);
    snew(dftb->phase2.x_pme, ct->atoms_cplx + ct->extcharges_cplx);
    snew(dftb->phase2.pot, ct->atoms_cplx + ct->extcharges_cplx);
    snew(dftb->phase2.pot2, ct->atoms_cplx);
    snew(dftb->phase2.pot3, ct->atoms_cplx);
    snew(dftb->phase2.pot4, ct->atoms_cplx);
    snew(dftb->phase2.nrnb_pme, 1);
    snew(dftb->phase2.pmedata, 1);
    status = gmx_pme_init_dftb(dftb->phase2.pmedata, ir, ct->atoms_cplx + ct->extcharges_cplx);
    return;
}


#ifdef GMX_MPI
void do_pme_for_dftb(charge_transfer_t *ct, dftb_t *dftb, int ct_mpi_rank, int ct_mpi_size)
#else
void do_pme_for_dftb(charge_transfer_t *ct, dftb_t *dftb)
#endif
{
    int           i, j, k, nn, ne, status;
    dftb_phase1_t dftb1;
    dftb_phase2_t dftb2;
    real          energy_pme, normbond;
    rvec          bond, x_qm, x_mm;
    dvec          dbond;
    t_pbc         pbc;

    set_pbc(&pbc, epbcXYZ, dftb->box_pme);

    for (i = 0; i < ct->sites; i++)
#ifdef GMX_MPI
    {
        if (i % ct_mpi_size == ct_mpi_rank)
#endif
    {
        dftb1 = dftb->phase1[i];
        nn    = dftb1.nn;
        ne    = dftb1.ne;
        /* do PME for DFTB here!
         * do it with the atomic charges obtained in the previous MD step
         * it is possibly slightly inaccurate,
         * but probably OK as a first guess */

        for (j = 0; j < nn; j++)
        {
            dftb1.pot[j] = 0.0;
        }

        /* PME -- long-range component (reciprocal space) */
        init_nrnb(dftb->phase1[i].nrnb_pme);
        gmx_pme_do_dftb(*(dftb->phase1[i].pmedata), 0, nn+ne, dftb1.x_pme, dftb1.q_pme, dftb->box_pme,
                        dftb->phase1[i].nrnb_pme, dftb->ewaldcoeff_pme, &energy_pme, dftb1.pot);

        /* PME -- corrections */
        for (j = 0; j < nn; j++)
        {
            /* exclude the QM-QM interactions as the shift will be calculated in DFTB for these interactions */
            dftb1.pot2[j] = 0.;
            for (k = 0; k < nn; k++)
            {
                if (j != k)
                {
                    dvec_sub(dftb1.x[j], dftb1.x[k], dbond);
                    dftb1.pot2[j] -= (-dftb1.qmat[k] + dftb->qzero1[dftb1.izp[k]])
                        * gmx_erf(dftb->ewaldcoeff_pme * dnorm(dbond) / NM_TO_BOHR) / dnorm(dbond); // this should be OK
                }
            }
            /* the self-interaction of charge densities */
            dftb1.pot3[j] = -dftb->ewaldcoeff_pme / NM_TO_BOHR * (-dftb1.qmat[j] + dftb->qzero1[dftb1.izp[j]]) / sqrt(M_PI); // this should be OK
        }

        /* PME -- short-range component (real space) using the previously created neighbor list */
        for (j = 0; j < nn; j++)
        {
            dftb1.pot4[j] = 0.;
            x_qm[XX]      = (real) dftb1.x[j][XX] / NM_TO_BOHR;
            x_qm[YY]      = (real) dftb1.x[j][YY] / NM_TO_BOHR;
            x_qm[ZZ]      = (real) dftb1.x[j][ZZ] / NM_TO_BOHR;
            /* loop over the neighbors */
            for (k = 0; k < dftb1.neighbors_pme[j]; k++)
            {
                x_mm[XX] = (real) dftb1.xe[dftb1.neighbor_pme[j][k]][XX] / NM_TO_BOHR;
                x_mm[YY] = (real) dftb1.xe[dftb1.neighbor_pme[j][k]][YY] / NM_TO_BOHR;
                x_mm[ZZ] = (real) dftb1.xe[dftb1.neighbor_pme[j][k]][ZZ] / NM_TO_BOHR;
                /* calculate the distance (considering PBC) */
                status   = pbc_dx_aiuc(&pbc, x_qm, x_mm, bond);
                normbond = norm(bond);
                if (normbond < dftb->rcoulomb_pme)
                {
                    dftb1.pot4[j] += dftb1.ze[dftb1.neighbor_pme[j][k]] / normbond * gmx_erfc(dftb->ewaldcoeff_pme * normbond);
                }
            }
        }

        //printf("Ewald potential kJ/mol/e = external shift\n");
        for (j = 0; j < nn; j++)
        {
            /* convert kJ/mol/e to atomic units of ESP */
            dftb1.pot[j] *= KJMOL_TO_HARTREE;
            if (i == 0 && j == 0)
            {
                printf("%12.7f %12.7f %12.7f %12.7f %12.7f ", dftb1.pot[j], dftb1.pot2[j], dftb1.pot3[j], dftb1.pot4[j] / NM_TO_BOHR, -dftb1.qmat[j] + dftb->qzero1[dftb1.izp[j]]);
            }
            dftb1.pot[j] += dftb1.pot2[j] + dftb1.pot3[j] + dftb1.pot4[j] / NM_TO_BOHR;
            //printf("%12.7f\n", dftb1.pot[j]);
        }

        /* save some an average value in the variable esp */
        dftb->phase1[i].esp = 0.;
        for (j = 0; j < nn; j++)
        {
            dftb->phase1[i].esp += dftb1.pot[j] * dftb1.mass[j];
        }
        dftb->phase1[i].esp *= dftb1.inv_tot_mass / ct->esp_scaling_factor;
        /* the following would be a very approximative alternative...
         * dftb->phase1[i].esp = dftb1.pot[0];
         */

        /* result is in: dftb.phase1[i].pot */
    }
#ifdef GMX_MPI
    }
    if (ct_mpi_rank == 0)
    {
#endif
    /* do it here for the complex */
    /* POSSIBLE OPTIMIZATION -- USE THE DATA OBTAINED FOR THE FRAGMENTS?
     * BUT ATTENTION: NOT ALL OF THE DATA HAVE BEEN CALCULATED ON MPI_RANK==0 !!!
     */
    dftb2 = dftb->phase2;
    nn    = dftb2.nn;
    ne    = dftb2.ne;

    for (j = 0; j < nn; j++)
    {
        dftb2.pot[j] = 0.0;
    }

    /* PME -- long-range component (reciprocal space) */
    init_nrnb(dftb->phase2.nrnb_pme);
    gmx_pme_do_dftb(*(dftb->phase2.pmedata), 0, nn+ne, dftb2.x_pme, dftb2.q_pme, dftb->box_pme,
                    dftb->phase2.nrnb_pme, dftb->ewaldcoeff_pme, &energy_pme, dftb2.pot);

    /* PME -- corrections */
    for (j = 0; j < nn; j++)
    {
        /* exclude the QM-QM interactions as the shift will be calculated in DFTB for these interactions */
        dftb2.pot2[j] = 0.;
        for (k = 0; k < nn; k++)
        {
            if (j != k)
            {
                dvec_sub(dftb2.x[j], dftb2.x[k], dbond);
                dftb2.pot2[j] -= (-dftb2.qmat[k] + dftb->qzero2[dftb2.izp[k]])
                    * gmx_erf(dftb->ewaldcoeff_pme * dnorm(dbond) / NM_TO_BOHR) / dnorm(dbond);
            }
        }
        /* the self-interaction of charge densities */
        dftb2.pot3[j] = -dftb->ewaldcoeff_pme / NM_TO_BOHR * (-dftb2.qmat[j] + dftb->qzero2[dftb2.izp[j]]) / sqrt(M_PI);
    }

    /* PME -- short-range component (real space) using the previously created neighbor list */
    for (j = 0; j < nn; j++)
    {
        dftb2.pot4[j] = 0.;
        x_qm[XX]      = (real) dftb2.x[j][XX] / NM_TO_BOHR;
        x_qm[YY]      = (real) dftb2.x[j][YY] / NM_TO_BOHR;
        x_qm[ZZ]      = (real) dftb2.x[j][ZZ] / NM_TO_BOHR;
        /* loop over the neighbors */
        for (k = 0; k < dftb2.neighbors_pme[j]; k++)
        {
            x_mm[XX] = (real) dftb2.xe[dftb2.neighbor_pme[j][k]][XX] / NM_TO_BOHR;
            x_mm[YY] = (real) dftb2.xe[dftb2.neighbor_pme[j][k]][YY] / NM_TO_BOHR;
            x_mm[ZZ] = (real) dftb2.xe[dftb2.neighbor_pme[j][k]][ZZ] / NM_TO_BOHR;
            /* calculate the distance (considering PBC) */
            status   = pbc_dx_aiuc(&pbc, x_qm, x_mm, bond);
            normbond = norm(bond);
            if (normbond < dftb->rcoulomb_pme)
            {
                dftb2.pot4[j] += dftb2.ze[dftb2.neighbor_pme[j][k]] / normbond * gmx_erfc(dftb->ewaldcoeff_pme * normbond);
            }
        }
    }

    //printf("Ewald potential kJ/mol/e = external shift\n");
    for (j = 0; j < nn; j++)
    {
        /* convert kJ/mol/e to atomic units of ESP */
        dftb2.pot[j] *= KJMOL_TO_HARTREE;
        //printf("%12.7f %12.7f %12.7f %12.7f ", dftb1.pot[j], dftb1.pot2[j], dftb1.pot3[j], dftb1.pot4[j]);
        dftb2.pot[j] += dftb2.pot2[j] + dftb2.pot3[j] + dftb2.pot4[j] / NM_TO_BOHR;
        //printf("%12.7f\n", dftb1.pot[j]);
    }
    /* result is in: dftb.phase2.pot */
#ifdef GMX_MPI
}
#endif
    return;
}

#ifdef GMX_MPI
void do_pme_for_dftb_part1(charge_transfer_t *ct, dftb_t *dftb, int ct_mpi_rank, int ct_mpi_size)
#else
void do_pme_for_dftb_part1(charge_transfer_t *ct, dftb_t *dftb)
#endif
{
    int           i, j, k, nn, ne, status;
    dftb_phase1_t dftb1;
    dftb_phase2_t dftb2;
    real          energy_pme, normbond;
    rvec          bond, x_qm, x_mm;
    dvec          dbond;
    t_pbc         pbc;

    set_pbc(&pbc, epbcXYZ, dftb->box_pme);

    for (i = 0; i < ct->sites; i++)
#ifdef GMX_MPI
    {
        if (i % ct_mpi_size == ct_mpi_rank)
#endif
    {
        dftb1 = dftb->phase1[i];
        nn    = dftb1.nn;
        ne    = dftb1.ne;
        /* do PME for DFTB here! -- "preparatory" part, before actual SCC-DFTB calculation */

        /* PME -- short-range component (real space) using the previously created neighbor list */
        for (j = 0; j < nn; j++)
        {
            dftb1.pot4[j] = 0.;
            x_qm[XX]      = (real) dftb1.x[j][XX] / NM_TO_BOHR;
            x_qm[YY]      = (real) dftb1.x[j][YY] / NM_TO_BOHR;
            x_qm[ZZ]      = (real) dftb1.x[j][ZZ] / NM_TO_BOHR;
            /* loop over the neighbors */
            for (k = 0; k < dftb1.neighbors_pme[j]; k++)
            {
                x_mm[XX] = (real) dftb1.xe[dftb1.neighbor_pme[j][k]][XX] / NM_TO_BOHR;
                x_mm[YY] = (real) dftb1.xe[dftb1.neighbor_pme[j][k]][YY] / NM_TO_BOHR;
                x_mm[ZZ] = (real) dftb1.xe[dftb1.neighbor_pme[j][k]][ZZ] / NM_TO_BOHR;
                /* calculate the distance (considering PBC) */
                status   = pbc_dx_aiuc(&pbc, x_qm, x_mm, bond);
                normbond = norm(bond);
                if (normbond < dftb->rcoulomb_pme)
                {
                    dftb1.pot4[j] += dftb1.ze[dftb1.neighbor_pme[j][k]] / normbond * gmx_erfc(dftb->ewaldcoeff_pme * normbond);
                }
            }
        }
    }
#ifdef GMX_MPI
    }
#endif
    return;
}

void do_pme_for_dftb_part2(charge_transfer_t *ct, dftb_t *dftb, int i)
{
    int           j, k, nn, ne, status;
    dftb_phase1_t dftb1;
    real          energy_pme, normbond;
    rvec          bond, x_qm, x_mm;
    dvec          dbond;
    t_pbc         pbc;

    set_pbc(&pbc, epbcXYZ, dftb->box_pme);

    dftb1 = dftb->phase1[i];
    nn    = dftb1.nn;
    ne    = dftb1.ne;

    /* transfer the QM charges to PME,
     * and reset the ESP array
     */
    for (j = 0; j < nn; j++)
    {
        dftb1.pot[j]   = 0.0;
        dftb1.q_pme[j] = -dftb1.qmat[j] + dftb->qzero1[dftb1.izp[j]];
    }

    /* PME calculation */
    init_nrnb(dftb->phase1[i].nrnb_pme);
    gmx_pme_do_dftb(*(dftb->phase1[i].pmedata), 0, nn+ne, dftb1.x_pme, dftb1.q_pme, dftb->box_pme,
                    dftb->phase1[i].nrnb_pme, dftb->ewaldcoeff_pme, &energy_pme, dftb1.pot);

    /* PME -- corrections */
    for (j = 0; j < nn; j++)
    {
        /* exclude the QM-QM interactions */
        dftb1.pot2[j] = 0.;
        for (k = 0; k < nn; k++)
        {
            if (j != k)
            {
                dvec_sub(dftb1.x[j], dftb1.x[k], dbond);
                dftb1.pot2[j] -= (-dftb1.qmat[k] + dftb->qzero1[dftb1.izp[k]])
                    * gmx_erf(dftb->ewaldcoeff_pme * dnorm(dbond) / NM_TO_BOHR) / dnorm(dbond);
            }
        }
        /* exclude the self-interaction of the charge densities */
        dftb1.pot3[j] = -dftb->ewaldcoeff_pme / NM_TO_BOHR * (-dftb1.qmat[j] + dftb->qzero1[dftb1.izp[j]]) / sqrt(M_PI);
    }

    /* dftb1.pot4 was pre-calculated */

    //printf("Ewald potential kJ/mol/e = external shift\n");
    for (j = 0; j < nn; j++)
    {
        /* convert kJ/mol/e to atomic units of ESP */
        dftb1.pot[j] *= KJMOL_TO_HARTREE;
        //if (i==0 && j==0)
        //  printf("%12.7f %12.7f %12.7f %12.7f %12.7f ", dftb1.pot[j], dftb1.pot2[j], dftb1.pot3[j], dftb1.pot4[j] / NM_TO_BOHR, -dftb1.qmat[j] + dftb->qzero1[dftb1.izp[j]]);
        dftb1.pot[j] += dftb1.pot2[j] + dftb1.pot3[j] + dftb1.pot4[j] / NM_TO_BOHR;
        //if (i==0 && j==0) printf("%12.7f ", dftb1.pot[j]);
    }

    /* save some an average value in the variable esp */
    dftb->phase1[i].esp = 0.;
    for (j = 0; j < nn; j++)
    {
        dftb->phase1[i].esp += dftb1.pot[j] * dftb1.mass[j];
    }
    dftb->phase1[i].esp *= dftb1.inv_tot_mass / ct->esp_scaling_factor;
    /* this was too approximative:
       dftb->phase1[i].esp = dftb1.pot[0];
     */

    return;
}

/* PME for DFTB phase 2 -- this will be done only once per MD step */
void do_pme_for_dftb_phase2(charge_transfer_t *ct, dftb_t *dftb)
{
    int           i, j, k, nn, ne, status;
    dftb_phase2_t dftb2;
    real          energy_pme, normbond;
    rvec          bond, x_qm, x_mm;
    dvec          dbond;
    t_pbc         pbc;

    set_pbc(&pbc, epbcXYZ, dftb->box_pme);

    dftb2 = dftb->phase2;
    nn    = dftb2.nn;
    ne    = dftb2.ne;

    for (j = 0; j < nn; j++)
    {
        dftb2.q_pme[j] = -dftb2.qmat[j] + dftb->qzero2[dftb2.izp[j]];
        dftb2.pot[j]   = 0.0;
    }

    /* PME -- long-range component (reciprocal space) */
    init_nrnb(dftb->phase2.nrnb_pme);
    gmx_pme_do_dftb(*(dftb->phase2.pmedata), 0, nn+ne, dftb2.x_pme, dftb2.q_pme, dftb->box_pme,
                    dftb->phase2.nrnb_pme, dftb->ewaldcoeff_pme, &energy_pme, dftb2.pot);

    /* PME -- corrections */
    for (j = 0; j < nn; j++)
    {
        /* exclude the QM-QM interactions as the shift will be calculated in DFTB for these interactions */
        dftb2.pot2[j] = 0.;
        for (k = 0; k < nn; k++)
        {
            if (j != k)
            {
                dvec_sub(dftb2.x[j], dftb2.x[k], dbond);
                dftb2.pot2[j] -= (-dftb2.qmat[k] + dftb->qzero2[dftb2.izp[k]])
                    * gmx_erf(dftb->ewaldcoeff_pme * dnorm(dbond) / NM_TO_BOHR) / dnorm(dbond);
            }
        }
        /* the self-interaction of charge densities */
        dftb2.pot3[j] = -dftb->ewaldcoeff_pme / NM_TO_BOHR * (-dftb2.qmat[j] + dftb->qzero2[dftb2.izp[j]]) / sqrt(M_PI);
    }

    /* PME -- short-range component (real space) using the previously created neighbor list */
    for (j = 0; j < nn; j++)
    {
        dftb2.pot4[j] = 0.;
        x_qm[XX]      = (real) dftb2.x[j][XX] / NM_TO_BOHR;
        x_qm[YY]      = (real) dftb2.x[j][YY] / NM_TO_BOHR;
        x_qm[ZZ]      = (real) dftb2.x[j][ZZ] / NM_TO_BOHR;
        /* loop over the neighbors */
        for (k = 0; k < dftb2.neighbors_pme[j]; k++)
        {
            x_mm[XX] = (real) dftb2.xe[dftb2.neighbor_pme[j][k]][XX] / NM_TO_BOHR;
            x_mm[YY] = (real) dftb2.xe[dftb2.neighbor_pme[j][k]][YY] / NM_TO_BOHR;
            x_mm[ZZ] = (real) dftb2.xe[dftb2.neighbor_pme[j][k]][ZZ] / NM_TO_BOHR;
            /* calculate the distance (considering PBC) */
            status   = pbc_dx_aiuc(&pbc, x_qm, x_mm, bond);
            normbond = norm(bond);
            if (normbond < dftb->rcoulomb_pme)
            {
                dftb2.pot4[j] += dftb2.ze[dftb2.neighbor_pme[j][k]] / normbond * gmx_erfc(dftb->ewaldcoeff_pme * normbond);
            }
        }
    }

    //printf("Ewald potential kJ/mol/e = external shift\n");
    for (j = 0; j < nn; j++)
    {
        /* convert kJ/mol/e to atomic units of ESP */
        dftb2.pot[j] *= KJMOL_TO_HARTREE;
        //printf("%12.7f %12.7f %12.7f %12.7f ", dftb1.pot[j], dftb1.pot2[j], dftb1.pot3[j], dftb1.pot4[j]);
        dftb2.pot[j] += dftb2.pot2[j] + dftb2.pot3[j] + dftb2.pot4[j] / NM_TO_BOHR;
        //printf("%12.7f\n", dftb1.pot[j]);
    }
    /* result is in: dftb.phase2.pot */
    return;
}

#ifdef GMX_MPI
void do_neighborlist_for_dftb(charge_transfer_t *ct, dftb_t *dftb, rvec *x, int ct_mpi_rank, int ct_mpi_size)
#else
void do_neighborlist_for_dftb(charge_transfer_t *ct, dftb_t *dftb, rvec *x)
#endif
{
    int           i, j, k, nn, ne, counter, status;
    dftb_phase1_t dftb1;
    dftb_phase2_t dftb2;
    rvec          bond;
    t_pbc         pbc;
    real          rcoulomb2;

    set_pbc(&pbc, epbcXYZ, dftb->box_pme);

//printf("pbc %d\n", pbc.ePBCDX);
/*
    switch (pbc.ePBCDX)
    {
        case epbcdxRECTANGULAR:
   printf("tric\n");
        case epbcdxTRICLINIC:
   printf("rect\n");
     }
 */


    //for (i = 0; i < DIM; i++)
    //printf("hbox %f fbox %f mhbox %f \n", pbc.hbox_diag[i], pbc.fbox_diag[i],  pbc.mhbox_diag[i]);
    rcoulomb2 = dftb->rcoulomb_pme * dftb->rcoulomb_pme;

    //printf("  creating PME DFTB neighborlists...");
    for (i = 0; i < ct->sites; i++)
#ifdef GMX_MPI
    {
        if (i % ct_mpi_size == ct_mpi_rank)
#endif
    {
        dftb1 = dftb->phase1[i];
        nn    = dftb1.nn;
        ne    = dftb1.ne;
        /* create a neighborlist for PME calculation for QM/MM DFTB
         * in a naive O(N^2) way:
         * double loop over QM atoms and over MM atoms */

        /* do it for every QM atom */
        for (j = 0; j < nn; j++)
        {
            counter = 0;
            /* check every MM atom */
            for (k = 0; k < ne; k++)
            {
                status = pbc_dx_aiuc(&pbc, x[ct->site[i].atom[j]], x[ct->site[i].extcharge[k]], bond);
                if (norm2(bond) < rcoulomb2)
                {
                    dftb1.neighbor_pme[j][counter] = k;
                    counter++;
                    // printf("Fragment %d, atom %d: neighbor no. %d is extcharge %d\n", i+1, j+1, counter, k+1);
                    if (counter == MAX_PME_NEIGHBORS)
                    {
                        fprintf(stderr, "Too many PME neighbors for atom %d of fragment %d\n  Exiting !!!\n\n", j+1, i+1);
                        exit(-1);
                    }
                }
            }
//      printf("bond %17.8f \n", norm2(bond));
            dftb1.neighbors_pme[j] = counter;
            //fprintf(stderr, "NS for PME/DFTB: atom %d in fragment %d has %d MM neighbors\n", j+1, i+1, counter);
        }
    }
#ifdef GMX_MPI
    }
#endif

#ifdef GMX_MPI
    if (ct_mpi_rank == 0)
    {
#endif
    /* do it here for the complex */
    dftb2 = dftb->phase2;
    nn    = dftb2.nn;
    ne    = dftb2.ne;
    /* do it for every QM atom */
    for (j = 0; j < nn; j++)
    {
        counter = 0;
        /* check every MM atom */
        for (k = 0; k < ne; k++)
        {
            status = pbc_dx_aiuc(&pbc, x[ct->atom_cplx[j]], x[ct->extcharge_cplx[k]], bond);
            if (norm2(bond) < rcoulomb2)
            {
                dftb2.neighbor_pme[j][counter] = k;
                counter++;
                if (counter == MAX_PME_NEIGHBORS)
                {
                    fprintf(stderr, "Too many PME neighbors for atom %d of complex\n  Exiting !!!\n\n", j+1);
                    exit(-1);
                }
            }
        }
        dftb2.neighbors_pme[j] = counter;
        //fprintf(stderr, "NS for PME/DFTB: atom %d in complex has %d MM neighbors\n", j+1, counter);
    }
#ifdef GMX_MPI
}
#endif
    return;
}

