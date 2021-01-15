#include "gromacs/transfer/transfer.h"


/******************************************
 * PREREQUISITIES FOR A QM/MM CALCULATION *
 *    TO BE PERFORMED IN EVERY MD STEP    *
 ******************************************/
//do_SCC?????????
void prepare_charge_transfer(matrix state_box, t_mdatoms *mdatoms, dftb_t *dftb, charge_transfer_t *ct, rvec *x_ct)
{
    int           i, j, k, l, counter; // m, n, counter2;
    ivec          shift, shiftmin;
    double        bond_length, mindist, curdist, sum;
    dvec          bond, box, image, com, coord, masscoord, r;
    dftb_phase1_t dftb1;
    dftb_phase2_t dftb2;

    /* debug begin
       printf("Information about system\n");
       printf("Number of atoms: %d\n", mdatoms->nr);
       printf("mdatoms->massA = %p, mdatoms->massT = %p, mdatoms->chargeA = %p\n", mdatoms->massA, mdatoms->massT, mdatoms->chargeA);
       printf("Selected atoms - massT, charge:\n");
       for (i=0; i<10; i++)
       printf("%d %12.7f %12.7f\n", i, mdatoms->massT[i], mdatoms->chargeA[i]);
       for (i=7270; i<7280; i++)
       printf("%d %12.7f %12.7f\n", i, mdatoms->massT[i], mdatoms->chargeA[i]);
       for (i=8270; i<8280; i++)
       printf("%d %12.7f %12.7f\n", i, mdatoms->massT[i], mdatoms->chargeA[i]);
       debug end */

    //printf("prepare_charge_transfer\n");

    /* read the box dimensions */
    for (j = 0; j < DIM; j++)
    {
        box[j] = state_box[j][j] * NM_TO_BOHR;
        //printf("BOX     %12.7f %12.7f %12.7f\n", state_box[j][0], state_box[j][1], state_box[j][2]);
    }
    if (ct->qmmm == 3)
    {
        copy_mat(state_box, dftb->box_pme);
        //printf("BOX_PME %12.7f %12.7f %12.7f\n", dftb->box_pme[XX][XX], dftb->box_pme[YY][YY], dftb->box_pme[ZZ][ZZ]);
    }


    /* read coordinates of the quantum system */
    /* conversion from nanometer to bohr */
    counter = 0;
    for (i = 0; i < ct->sites; i++)
    {
        for (j = 0; j < ct->site[i].atoms; j++)
        {
            for (k = 0; k < 3; k++)
            {
                dftb->phase1[i].x[j][k] = x_ct[ct->site[i].atom[j]][k] * NM_TO_BOHR;
            }
        }
        /* the link hydrogen atoms */
        for (j = 0; j < ct->site[i].bonds; j++)
        {
            dvec_sub(dftb->phase1[i].x[ct->site[i].MMLA[j]], dftb->phase1[i].x[ct->site[i].QMLA[j]], bond);
            bond_length = dnorm(bond);
            for (k = 0; k < DIM; k++)
            {
                switch (dftb->phase1[i].izp[ct->site[i].QMLA[j]])
                {
                    case 0: dftb->phase1[i].x[ct->site[i].MMLA[j]][k] = dftb->phase1[i].x[ct->site[i].QMLA[j]][k] + bond[k] / bond_length  * NM_TO_BOHR * CH_BOND_LENGTH / 10; break; // bondlength in Angstrom
                    case 2: dftb->phase1[i].x[ct->site[i].MMLA[j]][k] = dftb->phase1[i].x[ct->site[i].QMLA[j]][k] + bond[k] / bond_length  * NM_TO_BOHR * NH_BOND_LENGTH / 10; break; // bondlength in Angstrom
                    case 3: dftb->phase1[i].x[ct->site[i].MMLA[j]][k] = dftb->phase1[i].x[ct->site[i].QMLA[j]][k] + bond[k] / bond_length  * NM_TO_BOHR * OH_BOND_LENGTH / 10; break; // bondlength in Angstrom
                    case 4: dftb->phase1[i].x[ct->site[i].MMLA[j]][k] = dftb->phase1[i].x[ct->site[i].QMLA[j]][k] + bond[k] / bond_length  * NM_TO_BOHR * SH_BOND_LENGTH / 10; break; // bondlength in Angstrom
                    case 8: dftb->phase1[i].x[ct->site[i].MMLA[j]][k] = dftb->phase1[i].x[ct->site[i].QMLA[j]][k] + bond[k] / bond_length  * NM_TO_BOHR * BH_BOND_LENGTH / 10; break; // bondlength in Angstrom
                    case 9: dftb->phase1[i].x[ct->site[i].MMLA[j]][k] = dftb->phase1[i].x[ct->site[i].QMLA[j]][k] + bond[k] / bond_length  * NM_TO_BOHR * PH_BOND_LENGTH / 10; break; // bondlength in Angstrom
                    default: printf("unknown element - QM link atom must be either C, N, O or S.\n"); exit(-1);
                }
            }
        }
        /* copy to the phase2 - coordinates of the complex */
        for (j = 0; j < ct->site[i].atoms; j++)
        {
            copy_dvec(dftb->phase1[i].x[j], dftb->phase2.x[counter]);
            counter++;
        }
    }

/*
   // test - write out the coordinates
   counter=1;
   for (i=0; i<ct->sites; i++) {
   //i=0;
    printf("Site %d - %d atoms\n", i+1, dftb->phase1[i].nn);
    for (j=0; j<dftb->phase1[i].nn; j++){
      //printf("%d %d %12.7f%12.7f%12.7f\n",counter, dftb->phase1[i].izp[j]+1, dftb->phase1[i].x[j][0]*0.52, dftb->phase1[i].x[j][1]*0.52, dftb->phase1[i].x[j][2]*0.52);
      printf("C %12.7f%12.7f%12.7f\n", dftb->phase1[i].x[j][0]/NM_TO_BOHR*10, dftb->phase1[i].x[j][1]/NM_TO_BOHR*10, dftb->phase1[i].x[j][2]/NM_TO_BOHR*10);
      counter++;
    }
   }
   // */


//write out total QM zone
    if (ct->first_step)
    {
        char c;
        printf("%d \n", dftb->phase2.nn);
        printf("Complex .xyz\n");
        for (j = 0; j < dftb->phase2.nn; j++)
        {
            switch (dftb->phase2.izp[j])
            {
                case 0: c  = 'C'; break;
                case 1: c  = 'H'; break;
                case 2: c  = 'N'; break;
                case 3: c  = 'O'; break;
                case 4: c  = 'S'; break;
                case 5: c  = 'F'; break;
                case 6: c  = 'K'; break;
                case 7: c  = 'L'; break;
                case 8: c  = 'B'; break;
                case 9: c  = 'P'; break;
                case 10: c  = 'I'; break;
                default: c = 'X'; break;
            }
            printf("%c %12.7f%12.7f%12.7f\n", c, dftb->phase2.x[j][0]/NM_TO_BOHR*10, dftb->phase2.x[j][1]/NM_TO_BOHR*10, dftb->phase2.x[j][2]/NM_TO_BOHR*10);
        }
    }

    /* get the center of mass of every fragment as well as of the complex */
    for (i = 0; i < ct->sites; i++)
    {
        clear_dvec(com);
        for (j = 0; j < ct->site[i].atoms; j++)
        {
            coord[XX] = dftb->phase1[i].x[j][XX];
            coord[YY] = dftb->phase1[i].x[j][YY];
            coord[ZZ] = dftb->phase1[i].x[j][ZZ];
            dsvmul(mdatoms->massT[ct->site[i].atom[j]], coord, masscoord);
            dvec_inc(com, masscoord);
        }
        // dsvmul(dftb->phase1[i].inv_tot_mass, masscoord, dftb->phase1[i].com); - WRONG, ISN'T IT???
        dsvmul(dftb->phase1[i].inv_tot_mass, com, dftb->phase1[i].com);
        //printf("COM base %d: %f %f %f\n", i+1, dftb->phase1[i].com[XX] * 0.52, dftb->phase1[i].com[YY] * 0.52, dftb->phase1[i].com[ZZ] * 0.52);
    }

    // construct the Hubbard matrix (MOVED TO md.c) //

    /* read coordinates and magnitudes of the external charges */
    /* attention - consider the extcharges to be in the nearest periodic image! */
    if (ct->qmmm > 0)
    {
        for (i = 0; i < ct->sites; i++)
        {
            for (j = 0; j < ct->site[i].extcharges; j++)
            {
                /* coordinates */
                for (k = 0; k < DIM; k++)
                {
                    dftb->phase1[i].xe[j][k] = NM_TO_BOHR * x_ct[ct->site[i].extcharge[j]][k];
                }
                if (ct->qmmm < 3)
                {
                    //if (0) {  // we will use pbc where molecules are kept whole (do_pbc_mtop , see above) // edit: for MD we will use site energies obtained in phase1 for CT-Hamiltonian. in phase 2 different sites have different environments (QM vs MM) and therefore enrgies (very bad!).
                    /* when not doing particle--mesh Ewald: */
                    /* identify the nearest periodic image */
                    /* attention - doesn't necessarily keep molecules whole
                        induces artificial dipole !!! */
                    mindist = 1.e10;
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
                                curdist = dnorm(bond);
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
                /* magnitude, then add remaining charge */
                dftb->phase1[i].ze[j] = mdatoms->chargeA[ct->site[i].extcharge[j]];
            }
            for (k = 0; k < ct->site[i].bonds; k++)
            {
                for (l = 0; l < ct->site[i].addchrs[k]; l++)
                {
                    if (ct->site[i].modif_extcharge[k][l] > -1)
                    {
                        dftb->phase1[i].ze[ct->site[i].modif_extcharge[k][l]] += ct->site[i].extracharge[k] / ct->site[i].addchrs[k];/* the remaining charge (for electro-neutrality) is divided over specified atoms */
                    }
                }
            }
        }
    }

    /* test - write out the extcharges for nucleobase 1
       i=0;
       for (j=0; j<ct->site[i].extcharges; j++) // ct->extcharges[i]; j++)
        printf("H %12.7f%12.7f%12.7f%12.7f\n", dftb->phase1[i].xe[j][XX]*0.52, dftb->phase1[i].xe[j][YY]*0.52, dftb->phase1[i].xe[j][ZZ]*0.52 , dftb->phase1[i].ze[j]);
        //printf("atom%d %12.7f%12.7f%12.7f\n", j, x_ct[j][XX]*10, x_ct[j][YY]*10, x_ct[j][ZZ]*10);
       //end test */

    /* center of mass of the complex */
    clear_dvec(com);
    for (j = 0; j < ct->atoms_cplx; j++)
    {
        coord[XX] = dftb->phase2.x[j][XX];
        coord[YY] = dftb->phase2.x[j][YY];
        coord[ZZ] = dftb->phase2.x[j][ZZ];
        dsvmul(mdatoms->massT[ct->atom_cplx[j]], coord, masscoord);
        dvec_inc(com, masscoord);
    }
    // dsvmul(dftb->phase1[i].inv_tot_mass, masscoord, dftb->phase1[i].com); - WRONG, ISN'T IT???
    dsvmul(dftb->phase2.inv_tot_mass, com, dftb->phase2.com);

    /* coordinates and magnitude of external charges for the complex */
    if (ct->qmmm > 0)
    {
        for (j = 0; j < ct->extcharges_cplx; j++)
        {
            for (k = 0; k < DIM; k++)
            {
                dftb->phase2.xe[j][k] = NM_TO_BOHR * x_ct[ct->extcharge_cplx[j]][k];
            }
            if (ct->qmmm < 3)
            {
                mindist = 1.e10;
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
                            curdist = dnorm(bond);
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
            /* magnitude, then add remaining charge */
            dftb->phase2.ze[j] = mdatoms->chargeA[ct->extcharge_cplx[j]];
        }
        counter = 0;
        for (i = 0; i < ct->sites; i++)
        {
            for (j = 0; j < ct->site[i].bonds; j++)
            {
                for (k = 0; k < ct->site[i].addchrs[j]; k++)
                {
                    if (ct->modif_extcharge_cplx[counter] > -1)
                    {
                        dftb->phase2.ze[ct->modif_extcharge_cplx[counter]] += ct->site[i].extracharge[j] / ct->site[i].addchrs[j];
                    }
                    counter++;
                }
            }
        }
    }
/*
   // test - write out the extcharges for complex
   for (j=0; j<ct->extcharges_cplx; j++)
      printf("H %12.7f%12.7f%12.7f   %12.7f\n", dftb->phase2.xe[j][XX]*0.52, dftb->phase2.xe[j][YY]*0.52, dftb->phase2.xe[j][ZZ]*0.52, dftb->phase2.ze[j]);
   //end test */
/*
   // begin debug - check of sum of extcharges
   printf("Site   sum of extcharges\n");
   for (i=0; i<ct->sites; i++) {
    sum = 0.0;
    for (j=0; j<ct->site[i].extcharges; j++)
      sum += dftb->phase1[i].ze[j];
    printf("%d %12.7f\n", i+1, sum);
   }
   sum = 0.0;
   for (j=0; j<ct->extcharges_cplx; j++)
    sum += dftb->phase2.ze[j];
   printf("Complex: sum of extcharges = %12.7f\n", sum);
   //   end debug */

    if (ct->qmmm == 3)
    {
        /* feed the cordinates into the PME arrays */
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

    /* calculate ESP (electro-static potentials) at the nucleobases */
    //printf("Site   ESP (V)\n");
    if (ct->qmmm == 1 || ct->qmmm == 2)
    {
        for (i = 0; i < ct->sites; i++)
        {
            sum = 0.0;
            for (j = 0; j < ct->site[i].extcharges; j++)
            {
                dvec_sub(dftb->phase1[i].xe[j], dftb->phase1[i].com, r);
                sum += dftb->phase1[i].ze[j] / dnorm(r);
            }
            dftb->phase1[i].esp = sum / ct->esp_scaling_factor;
            //printf("%d %12.7f\n", i+1, sum * AU_OF_ESP_TO_VOLT);
        }
    }

    return;
}

