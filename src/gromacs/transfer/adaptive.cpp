#include "gromacs/transfer/transfer.h"

int adapt_QMzone(charge_transfer_t *ct, rvec *x_ct, t_mdatoms *mdatoms, const gmx_mtop_t *top_global, matrix state_box, rvec *gromacs_x)
{
    int    i, j, l, k, counter, nearest_inactive = -1, farthest_active = -1, adapted = 0;
    dvec   bond, com, coord, masscoord;
    double best_inactive_dist, best_active_dist, tot_occ, norm;
    rvec   box_center, dx;

    /* get center of mass for every site */
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
        //printf("residue %d  COM[Angstrom]: %f %f %f\n", ct->pool_site[i].resnr, ct->pool_site[i].com[XX] * 0.52, ct->pool_site[i].com[YY] * 0.52, ct->pool_site[i].com[ZZ] * 0.52);
    }


    /* search for nearest inactive site */
    best_inactive_dist = 100000.0;
    for (i = 0; i < ct->pool_size; i++)
    {
        if (!ct->pool_site[i].active)
        {
            dvec_sub(ct->pool_site[i].com, ct->coc, bond);
            if (dnorm(bond) < best_inactive_dist)
            {
                best_inactive_dist = dnorm(bond);
                nearest_inactive   = i;
            }
        }
    }
    //printf("nearest site is %d (res %d) dist %f\n",nearest_inactive,ct->pool_site[nearest_inactive].resnr, best_inactive_dist);

    /* search for farthermost active site */
    best_active_dist = 0.0;
    for (i = 0; i < ct->sites; i++)
    {
        dvec_sub(ct->site[i].com, ct->coc, bond);
        if (dnorm(bond) > best_active_dist)
        {
            best_active_dist = dnorm(bond);
            farthest_active  = i;
        }
    }
    //printf("farthest site is %d (res %d) dist %f\n",farthest_active,ct->site[farthest_active].resnr, best_active_dist);



    /* substitute farthermost active with nearer inactive site */
    if (best_active_dist < best_inactive_dist)
    {
        printf("Did not adapt QM zone. Is already optimal.\n");
    }
    else
    {
        tot_occ = 0.0;
        counter = 0;
        for (i = 0; i < ct->sites; i++)
        {
            for (j = 0; j < ct->site[i].homos; j++)
            {
                if (i == farthest_active)
                {
                    tot_occ += ct->occupation[counter];
                }
                counter++;
            }
        }
        if (tot_occ * ct->sites > 0.1)//more than 10% of the average occupation of the other sites.
        {
            printf("Aborted replacement of residue %d with %d because of non-neglible occupation %f.\n", ct->site[farthest_active].resnr, ct->pool_site[nearest_inactive].resnr, tot_occ);
        }
        else
        {
            printf("Replacing residue %d with %d. Distance to COC reduced form %fnm to %fnm.\n", ct->site[farthest_active].resnr, ct->pool_site[nearest_inactive].resnr, best_active_dist/NM_TO_BOHR, best_inactive_dist/NM_TO_BOHR);
            counter = 0;
            for (i = 0; i < ct->sites; i++)
            {
                for (j = 0; j < ct->site[i].homos; j++)
                {
                    if (i == farthest_active)
                    {
                        ct->wf[counter] = 0.0; ct->wf[counter+ct->dim] = 0.0;
                    }
                    counter++;
                }
            }
            ct->pool_site[nearest_inactive].active = 1;
            for (i = 0; i < ct->pool_size; i++)
            {
                if (ct->site[farthest_active].resnr == ct->pool_site[i].resnr)
                {
                    ct->pool_site[i].active = 0;
                }
            }
            ct->site[farthest_active] = ct->pool_site[nearest_inactive];

            /* center all atoms around COC */
            calc_box_center(ecenterTRIC, state_box, box_center);
            for (i = 0; i < DIM; i++)
            {
                dx[i] = box_center[i] - (real) ct->coc[i]/NM_TO_BOHR;
            }

            for (i = 0; i < top_global->natoms; i++)
            {
                rvec_inc(gromacs_x[i], dx);
                rvec_inc(x_ct[i], dx);
            }
            //printf("coc is %f %f %f  nm.\n ", ct->coc[0]/NM_TO_BOHR,ct->coc[1]/NM_TO_BOHR,ct->coc[2]/NM_TO_BOHR);
            //printf("box_cneter is %f %f %f  nm.\n ", box_center[0],box_center[1], box_center[2]);
            for (i = 0; i < DIM; i++)
            {
                ct->coc[i]       += (double) dx[i] * NM_TO_BOHR;
                ct->coc_start[i] += (double) dx[i] * NM_TO_BOHR;
            }
            printf("displaced all atoms by %f %f %f  nm.\n ", dx[0], dx[1], dx[2]);

            /* rescale the wavefunction */
            norm = 0;
            for (i = 0; i < ct->dim; i++)
            {
                norm += SQR(ct->wf[i]) + SQR(ct->wf[i+ct->dim]);
            }
            //printf("norm has changed to %f due to adaption of the QM zone.\n", norm);
            for (i = 0; i < 2*ct->dim; i++)
            {
                ct->wf[i] /= sqrt(norm);
            }

            /* set the active atoms of the complex */
            counter = 0;
            for (i = 0; i < ct->sites; i++)
            {
                for (j = 0; j < ct->site[i].atoms; j++)
                {
                    ct->atom_cplx[counter]     = ct->site[i].atom[j];
                    ct->atomtype_cplx[counter] = ct->site[i].atomtype[j];
                    counter++;
                }
            }

            /* set the new external charges of the complex */
            if (ct->qmmm > 0)
            {
                for (i = 0; i < top_global->natoms; i++)
                {
                    ct->extcharge_cplx[i] = i;
                }
                for (i = 0; i < ct->sites; i++)
                {
                    k = find_intersection(top_global->natoms, ct->extcharge_cplx, ct->site[i].extcharge, ct->extcharge_cplx); // this should successively reduce the charges in ct->extcharge_cplx.
                    for (j = k+1; j <= top_global->natoms; j++)                                                               // k is index of highest common entry
                    {
                        ct->extcharge_cplx[j] = -1;
                    }
                }
                counter = 0;
                for (l = 0; l < ct->extcharges_cplx; l++)
                {
                    for (i = 0; i < ct->sites; i++)
                    {
                        for (j = 0; j < ct->site[i].bonds; j++)
                        {
                            for (k = 0; k < ct->site[i].addchrs[j]; k++)
                            {
                                if (ct->extcharge_cplx[l] == ct->site[i].extcharge[ ct->site[i].modif_extcharge[j][k] ]) // if one of the extcharges of the complex is the same atom that was modified in the monomer calculation, then also modify it in the complex calculation.
                                {
                                    ct->modif_extcharge_cplx[counter] = l;
                                    counter++;
                                }
                            }
                        }
                    }
                }
            }
            adapted = 1;
        }
    }

    return adapted;
}

