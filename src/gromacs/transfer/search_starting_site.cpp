#include "gromacs/transfer/transfer.h"


void search_starting_site(matrix                              state_box,
                          t_mdatoms*                          mdatoms,
                          const nonbonded_verlet_t*           nbv,
                          const rvec*                         shift_vec,
                          const t_commrec*                    cr,
                          charge_transfer_t*                  ct,
                          std::unique_ptr<QMMM_rec_transfer>* dftbplus_phase1,
                          std::unique_ptr<QMMM_rec_transfer>& dftbplus_phase2,
                          rvec*                               x_ct,
                          const gmx_mtop_t*                   top_global,
                          rvec*                               gromacs_x)
{
    int    original_nsites, i, j, lowest_site, counter;
    double Elowest, Esite;
 // rvec  *dummy_x;
 // matrix dummy_matrix;

    Elowest = 1000.0;
    lowest_site = -1;

    /* perform DFTB calculations on every site to find the lowest energy
       this will determine the starting site */
    printf("Scanning for site with lowest energy.\n scanning site:\n");
    original_nsites = ct->sites;
    ct->sites       = 1;
    for (i = 0; i < ct->pool_size; i++)
    {
        ct->site[0] = ct->pool_site[i];
        //printf("%d resnr %d addr %p\n",i, ct->pool_site[i].resnr, &(ct->pool_site[i].resnr) );
        //printf("%d resnr %d addr %p\n",i, ct->site[0].resnr, &(ct->site[0].resnr));
        prepare_charge_transfer(state_box, mdatoms, true, nbv, shift_vec, cr, ct, dftbplus_phase1, dftbplus_phase2, x_ct);
     // run_dftb1(ct, dftb, 0); TODO DFTB+
        printf("in search_starting_site(): missing DFTB routine!\n");
        exit(-1);
        for (j = 0; j < ct->site[0].homos; j++)
        {
     //     Esite  = dftb->phase1[0].ev[ct->site[0].homo[j]-1]; TODO DFTB+
            Esite *= ct->is_hole_transfer ? -1.0 : 1.0;
            if (Esite < Elowest)
            {
                Elowest     = Esite;
                lowest_site = i;
                copy_dvec(ct->site[0].centerOfMass, ct->coc);
                printf("found residue %d  MO %d  E = %f ha  COM[Angstrom] = %f %f %f\n", ct->pool_site[i].resnr, ct->pool_site[i].homo[j], Esite,  ct->coc[XX] * 0.52, ct->coc[YY] * 0.52, ct->coc[ZZ] * 0.52);
            }
        }
    }
    printf("starting residue is poolsite %d res %d  COC[Angstrom] = %f %f %f\n", lowest_site, ct->pool_site[lowest_site].resnr, ct->coc[XX] * 0.52, ct->coc[YY] * 0.52, ct->coc[ZZ] * 0.52);
    ct->sites = original_nsites;


    /* find nearest neighborst of the starting site */
    for (i = 0; i < ct->sites; i++)
    {
        ct->site[i] = ct->pool_site[i];
    }
    i = 0;
    while (adapt_QMzone(ct, x_ct, mdatoms, top_global, state_box, gromacs_x))
    {
        i++;
    }
    printf("Initially adapted QM zone %d times.\n", i);

    /* set starting wavefunction */
    counter = 0;
    for (i = 0; i < ct->sites; i++)
    {
        for (j = 0; j < ct->site[i].homos; j++)
        {
            if (ct->site[i].resnr == ct->pool_site[lowest_site].resnr && (j == ct->site[i].homos-1))
            {
                ct->wf[counter] = 1.0; ct->wf[counter+ct->dim] = 0.0;
            }
            else
            {
                ct->wf[counter] = 0.0; ct->wf[counter+ct->dim] = 0.0;
            }
            counter++;
        }
    }
    printf("starting wavefunction:\n");
    for (i = 0; i < ct->dim; i++)
    {
        printf("%10.7f %10.7f\n", ct->wf[i], ct->wf[i + ct->dim]);
    }


    return;
}
