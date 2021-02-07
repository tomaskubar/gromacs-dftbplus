#include "gromacs/transfer/transfer.h"


/******************************************
 * PREREQUISITIES FOR A QM/MM CALCULATION *
 *    TO BE PERFORMED IN EVERY MD STEP    *
 ******************************************/

//do_SCC?????????
void prepare_charge_transfer(const matrix                        state_box,
                             const t_mdatoms*                    mdatoms,
                             const bool                          bNS,
                             const nonbonded_verlet_t*           nbv,
                             const rvec*                         shift_vec,
                             const t_commrec*                    cr,
                             charge_transfer_t*                  ct,
                             std::unique_ptr<QMMM_rec_transfer>* dftbplus_phase1,
                             std::unique_ptr<QMMM_rec_transfer>& dftbplus_phase2,
                             rvec*                               x_ct)
{
 // int           i, j, k, l, counter; // m, n, counter2;
 // ivec          shift, shiftmin;
 // double        bond_length, mindist, curdist, sum;
 // dvec          bond, box, image, com, coord, masscoord, r;
 // dftb_phase1_t dftb1;
 // dftb_phase2_t dftb2;

    // DEBUG OUTPUT
 // printf("Information about system\n");
 // printf("Number of atoms: %d\n", mdatoms->nr);
 // printf("mdatoms->massA = %p, mdatoms->massT = %p, mdatoms->chargeA = %p\n", mdatoms->massA, mdatoms->massT, mdatoms->chargeA);
 // printf("Selected atoms - massT, charge:\n");
 // for (int i=0; i<10; i++)
 //     printf("%d %12.7f %12.7f\n", i, mdatoms->massT[i], mdatoms->chargeA[i]);
 // for (int i=7270; i<7280; i++)
 //     printf("%d %12.7f %12.7f\n", i, mdatoms->massT[i], mdatoms->chargeA[i]);
 // for (int i=8270; i<8280; i++)
 //     printf("%d %12.7f %12.7f\n", i, mdatoms->massT[i], mdatoms->chargeA[i]);

 // printf("prepare_charge_transfer\n");


    // In case the neighborlist has been updated,
    // perform a new search for MM neighbors of QM atoms.
    if (bNS)
	{
        for (int site=0; site<ct->sites; site++)
        {
	        dftbplus_phase1[site]->update_QMMMrec_verlet_ns(cr, nbv, x_ct, mdatoms, state_box);
	        dftbplus_phase1[site]->update_QMMMrec_dftb(cr, shift_vec, x_ct, mdatoms, state_box);
        }
	    dftbplus_phase2->update_QMMMrec_verlet_ns(cr, nbv, x_ct, mdatoms, state_box);
	    dftbplus_phase2->update_QMMMrec_dftb(cr, shift_vec, x_ct, mdatoms, state_box);
	}
    // Update the coordinates in any case.
    for (int site=0; site<ct->sites; site++)
    {
	    dftbplus_phase1[site]->update_QMMM_coord(cr, shift_vec, x_ct, mdatoms, state_box);
        // treat the link hydrogen atoms
        for (int j = 0; j < ct->site[site].bonds; j++)
        {
            dvec bond;
         // rvec_sub(dftbplus_phase1[site]->qm.x[ct->site[site].MMLA[j]], dftb->phase1[i].x[ct->site[site].QMLA[j]], bond);
            dvec_sub(dftbplus_phase1[site]->qm->xQM[ct->site[site].MMLA[j]],
                     dftbplus_phase1[site]->qm->xQM[ct->site[site].QMLA[j]],
                     bond);
            double bond_length = dnorm(bond);
            double new_bond_length;
            switch (dftbplus_phase1[site]->qm->atomicnumberQM_get(ct->site[site].QMLA[j]))
            {
                case  6: new_bond_length = CH_BOND_LENGTH; break;
                case  7: new_bond_length = NH_BOND_LENGTH; break;
                case  8: new_bond_length = OH_BOND_LENGTH; break;
                case 16: new_bond_length = SH_BOND_LENGTH; break;
                case  5: new_bond_length = BH_BOND_LENGTH; break;
                case 15: new_bond_length = PH_BOND_LENGTH; break;
                default: printf("unknown element - QM link atom must be either C, N, O, S, B or P.\n");
                         exit(-1);
            }
            for (int k = 0; k < DIM; k++)
            {
                dftbplus_phase1[site]->qm->xQM[ct->site[site].MMLA[j]][k] = dftbplus_phase1[site]->qm->xQM[ct->site[site].QMLA[j]][k]
                                + bond[k] * new_bond_length / bond_length;
            }
        }
    }
	dftbplus_phase2->update_QMMM_coord(cr, shift_vec, x_ct, mdatoms, state_box);
    // TODO: How to treat the link atoms in PHASE2 ?


    // test - write out the coordinates
 // int counter=1;
 // for (int i=0; i<ct->sites; i++)
 // {
 //  // i=0;
 //     printf("Site %d - %d atoms\n", i+1, dftb->phase1[i].nn);
 //     for (int j=0; j<dftb->phase1[i].nn; j++)
 //     {
 //      // printf("%d %d %12.7f%12.7f%12.7f\n", counter, dftb->phase1[i].izp[j]+1,
 //      //        dftb->phase1[i].x[j][0]*0.52, dftb->phase1[i].x[j][1]*0.52, dftb->phase1[i].x[j][2]*0.52);
 //         printf("C %12.7f%12.7f%12.7f\n",
 //                dftb->phase1[i].x[j][0]/NM_TO_BOHR*10, dftb->phase1[i].x[j][1]/NM_TO_BOHR*10, dftb->phase1[i].x[j][2]/NM_TO_BOHR*10);
 //         counter++;
 //     }
 // }


    // write out total QM zone
    if (ct->first_step)
    {
        char c[3];
        printf("%d \n", dftbplus_phase2->qm->nrQMatoms_get());
        printf("Complex .xyz\n");
        for (int j = 0; j < dftbplus_phase2->qm->nrQMatoms_get(); j++)
        {
            switch (dftbplus_phase2->qm->atomicnumberQM_get(j)) // dftb->phase2.izp[j])
            {
                case  6: strcpy(c,  "C"); break;
                case  1: strcpy(c,  "H"); break;
                case  7: strcpy(c,  "N"); break;
                case  8: strcpy(c,  "O"); break;
                case 16: strcpy(c,  "S"); break;
                case  9: strcpy(c,  "F"); break;
                case 35: strcpy(c, "Br"); break;
                case 17: strcpy(c, "Cl"); break;
                case  5: strcpy(c,  "B"); break;
                case 15: strcpy(c,  "P"); break;
                case 53: strcpy(c,  "I"); break;
                default: strcpy(c,  "X"); break;
            }
            printf("%2s %12.7f%12.7f%12.7f\n", c,
                   dftbplus_phase2->qm->xQM_get(j,0) * 10., dftbplus_phase2->qm->xQM_get(j,1) * 10., dftbplus_phase2->qm->xQM_get(j,2) * 10.);
        }
    }

    /* get the center of mass of every fragment as well as of the complex */
    for (int i = 0; i < ct->sites; i++)
    {
        dvec com, coord, masscoord;
        clear_dvec(com);
        for (int j = 0; j < ct->site[i].atoms; j++)
        {
            coord[XX] = dftbplus_phase1[i]->qm->xQM_get(j,XX);
            coord[YY] = dftbplus_phase1[i]->qm->xQM_get(j,YY);
            coord[ZZ] = dftbplus_phase1[i]->qm->xQM_get(j,ZZ);
            dsvmul(ct->site[i].mass[j], coord, masscoord);
            dvec_inc(com, masscoord);
        }
     // dsvmul(dftb->phase1[i].inv_tot_mass, masscoord, dftb->phase1[i].com); - WRONG, ISN'T IT???
     // dsvmul(dftb->phase1[i].inv_tot_mass, com, dftb->phase1[i].com);
        dsvmul(ct->site[i].inv_tot_mass, com, ct->site[i].centerOfMass);
     // printf("COM base %d: %f %f %f\n", i+1, dftb->phase1[i].com[XX] * 0.52, dftb->phase1[i].com[YY] * 0.52, dftb->phase1[i].com[ZZ] * 0.52);
    }

    // construct the Hubbard matrix (MOVED TO md.c) //

    // read coordinates and magnitudes of the external charges
    // attention - consider the extcharges to be in the nearest periodic image!
    if (ct->qmmm > 0)
    {
        for (int i = 0; i < ct->sites; i++)
        {
            for (int j = 0; j < ct->site[i].extcharges; j++)
            {
                /* COORDINATES ARE HANDLED BY "STANDARD" QM/MM INSTEAD
                // coordinates
                for (int k = 0; k < DIM; k++)
                {
                    dftb->phase1[i].xe[j][k] = NM_TO_BOHR * x_ct[ct->site[i].extcharge[j]][k];
                }
                if (ct->qmmm < 3)
                {
                    //if (0) {  // we will use pbc where molecules are kept whole (do_pbc_mtop , see above)
                                // edit: for MD we will use site energies obtained in phase1 for CT-Hamiltonian.
                                // in phase 2 different sites have different environments (QM vs MM) and therefore energies (very bad!).
                    // when not doing particle--mesh Ewald:
                    // identify the nearest periodic image
                    // attention - doesn't necessarily keep molecules whole
                    // induces artificial dipole !!!
                    double mindist = 1.e10;
                    clear_ivec(shiftmin);
                    for (shift[XX] = -1; shift[XX] <= 1; shift[XX]++)
                    {
                        for (shift[YY] = -1; shift[YY] <= 1; shift[YY]++)
                        {
                            for (shift[ZZ] = -1; shift[ZZ] <= 1; shift[ZZ]++)
                            {
                                image[XX] = dftb->phase1[i].xe[j][XX] + shift[XX] * box[XX];
                                image[YY] = dftb->phase1[i].xe[j][YY] + shift[YY] * box[YY];
                                image[ZZ] = dftb->phase1[i].xe[j][ZZ] + shift[ZZ] * box[ZZ];
                                dvec_sub(image, dftb->phase1[i].com, bond);
                                double curdist = dnorm(bond);
                                if (curdist < mindist)
                                {
                                    mindist = curdist;
                                    copy_ivec(shift, shiftmin);
                                }
                            }
                        }
                    }
                    dftb->phase1[i].xe[j][XX] += shiftmin[XX] * box[XX];
                    dftb->phase1[i].xe[j][YY] += shiftmin[YY] * box[YY];
                    dftb->phase1[i].xe[j][ZZ] += shiftmin[ZZ] * box[ZZ];
                }
                */
                // magnitude, then add remaining charge
             // dftb->phase1[i].ze[j] = mdatoms->chargeA[ct->site[i].extcharge[j]];
                dftbplus_phase1[i]->mm->MMcharges[j] = mdatoms->chargeA[ct->site[i].extcharge[j]];
            }
            for (int k = 0; k < ct->site[i].bonds; k++)
            {
                for (int l = 0; l < ct->site[i].addchrs[k]; l++)
                {
                    if (ct->site[i].modif_extcharge[k][l] > -1)
                    {
                        // the remaining charge (for electro-neutrality) is divided over specified atoms
                     // dftb->phase1[i].ze[ct->site[i].modif_extcharge[k][l]] += ct->site[i].extracharge[k] / ct->site[i].addchrs[k];
                        dftbplus_phase1[i]->mm->MMcharges[ct->site[i].modif_extcharge[k][l]] += ct->site[i].extracharge[k] / ct->site[i].addchrs[k];
                    }
                }
            }
        }
    }

    // test - write out the extcharges for nucleobase 1
 // {
 //     int i=0;
 //     for (int j=0; j<ct->site[i].extcharges; j++) // ct->extcharges[i]; j++)
 //         printf("H %12.7f%12.7f%12.7f%12.7f\n",
 //                dftb->phase1[i].xe[j][XX]*0.52, dftb->phase1[i].xe[j][YY]*0.52, dftb->phase1[i].xe[j][ZZ]*0.52,
 //                dftb->phase1[i].ze[j]);
 //  // printf("atom%d %12.7f%12.7f%12.7f\n", j, x_ct[j][XX]*10, x_ct[j][YY]*10, x_ct[j][ZZ]*10);
 // }
    // end test

    /* center of mass of the complex */
    {
        dvec com, coord, masscoord;
        clear_dvec(com);
        for (int j = 0; j < ct->atoms_cplx; j++)
        {
            coord[XX] = dftbplus_phase2->qm->xQM_get(j,XX);
            coord[YY] = dftbplus_phase2->qm->xQM_get(j,YY);
            coord[ZZ] = dftbplus_phase2->qm->xQM_get(j,ZZ);
            dsvmul(ct->mass_cplx[j], coord, masscoord);
            dvec_inc(com, masscoord);
        }
 //     dsvmul(dftb->phase1[i].inv_tot_mass, masscoord, dftb->phase1[i].com); - WRONG, ISN'T IT???
 //     dsvmul(dftb->phase2.inv_tot_mass, com, dftb->phase2.com);
        dsvmul(ct->inv_tot_mass_cplx, com, ct->centerOfMassCplx);
    }

    /* coordinates and magnitude of external charges for the complex */
    if (ct->qmmm > 0)
    {
        for (int j = 0; j < ct->extcharges_cplx; j++)
        {
            /* COORDINATES ARE HANDLED BY "STANDARD" QM/MM INSTEAD
            for (int k = 0; k < DIM; k++)
            {
                dftb->phase2.xe[j][k] = NM_TO_BOHR * x_ct[ct->extcharge_cplx[j]][k];
            }
            if (ct->qmmm < 3)
            {
                double mindist = 1.e10;
                clear_ivec(shiftmin);
                for (shift[XX] = -1; shift[XX] <= 1; shift[XX]++)
                {
                    for (shift[YY] = -1; shift[YY] <= 1; shift[YY]++)
                    {
                        for (shift[ZZ] = -1; shift[ZZ] <= 1; shift[ZZ]++)
                        {
                            image[XX] = dftb->phase2.xe[j][XX] + shift[XX] * box[XX];
                            image[YY] = dftb->phase2.xe[j][YY] + shift[YY] * box[YY];
                            image[ZZ] = dftb->phase2.xe[j][ZZ] + shift[ZZ] * box[ZZ];
                            dvec_sub(image, dftb->phase2.com, bond);
                            double curdist = dnorm(bond);
                            if (curdist < mindist)
                            {
                                mindist = curdist;
                                copy_ivec(shift, shiftmin);
                            }
                        }
                    }
                }
                dftb->phase2.xe[j][XX] += shiftmin[XX] * box[XX];
                dftb->phase2.xe[j][YY] += shiftmin[YY] * box[YY];
                dftb->phase2.xe[j][ZZ] += shiftmin[ZZ] * box[ZZ];
            }
            */
            // magnitude, then add remaining charge
         // dftb->phase2.ze[j] = mdatoms->chargeA[ct->extcharge_cplx[j]];
            dftbplus_phase2->mm->MMcharges[j] = mdatoms->chargeA[ct->extcharge_cplx[j]];
        }
        int counter = 0;
        for (int i = 0; i < ct->sites; i++)
        {
            for (int j = 0; j < ct->site[i].bonds; j++)
            {
                for (int k = 0; k < ct->site[i].addchrs[j]; k++)
                {
                    if (ct->modif_extcharge_cplx[counter] > -1)
                    {
                        dftbplus_phase2->mm->MMcharges[ct->modif_extcharge_cplx[counter]] += ct->site[i].extracharge[j] / ct->site[i].addchrs[j];
                    }
                    counter++;
                }
            }
        }
    }


    // test - write out the extcharges for complex
 // for (j=0; j<ct->extcharges_cplx; j++)
 //     printf("H %12.7f%12.7f%12.7f   %12.7f\n", dftb->phase2.xe[j][XX]*0.52, dftb->phase2.xe[j][YY]*0.52, dftb->phase2.xe[j][ZZ]*0.52, dftb->phase2.ze[j]);
    // end test

    // begin debug - check of sum of extcharges
    {
        printf("Site   sum of extcharges\n");
        for (int i=0; i<ct->sites; i++)
        {
            double sum = 0.0;
            for (int j=0; j<ct->site[i].extcharges; j++)
                sum += dftbplus_phase1[i]->mm->MMcharges[j];
            printf("%d %12.7f\n", i+1, sum);
        }
        double sum = 0.0;
        for (int j=0; j<ct->extcharges_cplx; j++)
            sum += dftbplus_phase2->mm->MMcharges[j];
        printf("Complex: sum of extcharges = %12.7f\n", sum);
    }
    // end debug

    /* WILL BE PERFORMED BY "STANDARD" QM/MM
    if (ct->qmmm == 3)
    {
        // feed the cordinates into the PME arrays
        for (i = 0; i < ct->sites; i++)
        {
            dftb1 = dftb->phase1[i];
            for (j = 0; j < dftb1.nn; j++)
            {
                dftb1.x_pme[j][0] = (real) dftb1.x[j][0] / NM_TO_BOHR;
                dftb1.x_pme[j][1] = (real) dftb1.x[j][1] / NM_TO_BOHR;
                dftb1.x_pme[j][2] = (real) dftb1.x[j][2] / NM_TO_BOHR;
                // dftb1.q_pme[j]    = 0.; // (real) (-dftb1.qmat[j] + dftb->qzero1[dftb1.izp[j]]);
            }
            for (j = 0; j < dftb1.ne; j++)
            {
                dftb1.x_pme[dftb1.nn + j][0] = (real) dftb1.xe[j][0] / NM_TO_BOHR;
                dftb1.x_pme[dftb1.nn + j][1] = (real) dftb1.xe[j][1] / NM_TO_BOHR;
                dftb1.x_pme[dftb1.nn + j][2] = (real) dftb1.xe[j][2] / NM_TO_BOHR;
                dftb1.q_pme[dftb1.nn + j]    = (real) dftb1.ze[j];
            }
        }
        dftb2 = dftb->phase2;
        for (j = 0; j < dftb2.nn; j++)
        {
            dftb2.x_pme[j][0] = (real) dftb2.x[j][0] / NM_TO_BOHR;
            dftb2.x_pme[j][1] = (real) dftb2.x[j][1] / NM_TO_BOHR;
            dftb2.x_pme[j][2] = (real) dftb2.x[j][2] / NM_TO_BOHR;
            // dftb2.q_pme[j]    = 0.; // (real) (-dftb2.qmat[j] + dftb->qzero2[dftb2.izp[j]]);
        }
        for (j = 0; j < dftb2.ne; j++)
        {
            dftb2.x_pme[dftb2.nn + j][0] = (real) dftb2.xe[j][0] / NM_TO_BOHR;
            dftb2.x_pme[dftb2.nn + j][1] = (real) dftb2.xe[j][1] / NM_TO_BOHR;
            dftb2.x_pme[dftb2.nn + j][2] = (real) dftb2.xe[j][2] / NM_TO_BOHR;
            dftb2.q_pme[dftb2.nn + j]    = (real) dftb2.ze[j];
        }
    }
    */

    // calculate ESP (electro-static potentials) at the sites
    // TODO: check if this is used anywhere, if it is calculated the right way, and if the unit is OK
 // printf("Site   ESP (V)\n");
    if (ct->qmmm == 1 || ct->qmmm == 2)
    {
        for (int i = 0; i < ct->sites; i++)
        {
            double sum = 0.0;
            for (int j = 0; j < ct->site[i].extcharges; j++)
            {
                dvec r;
                dvec_sub(dftbplus_phase1[i]->mm->xMM[j], ct->site[i].centerOfMass, r);
                sum += dftbplus_phase1[i]->mm->MMcharges[j] / dnorm(r);
            }
            ct->site[i].esp = sum / ct->esp_scaling_factor;
         // printf("%d %12.7f\n", i+1, sum * AU_OF_ESP_TO_VOLT);
        }
    }
} // prepare_charge_transfer()

