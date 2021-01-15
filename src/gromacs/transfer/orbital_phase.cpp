#include "gromacs/transfer/transfer.h"

void check_and_invert_orbital_phase(dftb_phase1_t *dftb1, charge_transfer_t *ct, t_state *state_global, t_mdatoms *mdatoms)
{
    int    iatom, i, nn, j, k, l, m, iao, jao;
    double norm1, norm2, proj;
//  use p-orbital and normal vector of the molecule plane to fix the sign of electronic coupling --W.X
    double vector1[3], vector2[3], normal_vector[3], ***coords;
    int    *orb_index, counter;


    if (ct->first_step && ct->define_orbital_sign)
    {
       snew(orb_index, 3);

       snew(coords, ct->sites);
       for (i = 0; i < ct->sites; i++)
       {
           snew(coords[i], 3);
           for (j = 0; j < 3; j++)
               snew(coords[i][j], DIM);
       }

       for (iatom = 0; iatom < mdatoms->homenr; iatom++)
       {
           for (i = 0; i < ct->sites; i++)
           {
               for (j = 0; j < 3; j++)
               {
                   if (iatom == ct->site[i].atom[ct->atom_index_sign[j]]-1)
                   {
                      for (k = 0; k < DIM; k++)
                          coords[i][j][k] = state_global->x[iatom][k]*NM_TO_BOHR;
                   }
               }
           }
       }
    }


    for (i = 0; i < ct->sites; i++)
    {
        if (ct->first_step)
        {
//  p-orbital scheme to fix the sign of coefficients
// Define normal vector of the molecule plane
           if (ct->define_orbital_sign)
           {
                for (j = 0; j < DIM; j++)
                {
                    vector1[j] = coords[i][0][j] - coords[i][2][j];
                    vector2[j] = coords[i][1][j] - coords[i][2][j];
                }

                normal_vector[0] = vector1[1]*vector2[2] - vector1[2]*vector2[1];
                normal_vector[1] = vector1[2]*vector2[0] - vector1[0]*vector2[2];
                normal_vector[2] = vector1[0]*vector2[1] - vector1[1]*vector2[0];

                counter = 0; nn = 0;

                for (j = 0; j < ct->site[i].atoms; j++)
                {
                    if (j == ct->atom_index_sign[counter]-1)
                    {
                         orb_index[counter] = nn + 1;
                         counter++;

                         if (counter == 3)  break;
                    }

                    switch (ct->site[i].atomtype[j])
                    {
                      case 1:  nn += 1; break; // H
                      case 4 : nn += 9; break; // S
                      case 9 : nn += 9; break; // P
                      default: nn += 4; break; // B, C, O, N, F, Cl, Br, I
                    }
                }

//                       for (k = 0; k < dftb1[i].norb; k++)
//                           printf("%f \n",dftb1[i].a[k][ct->site[i].homo[0]-1]);
//                       exit(0);

//only consider homos
                l = ct->site[i].homo[0]-1;
                counter = 1;

                norm1 = 0.0; norm2 = 0.0; proj = 0.0;

                for (j = 0; j < DIM; j++)
                {
                    norm1 += SQR(normal_vector[j]);
                    norm2 += SQR(dftb1[i].a[orb_index[counter]+j][l]);

                    proj  += normal_vector[j]*dftb1[i].a[orb_index[counter]+j][l];
                }

                proj /= sqrt(norm1)*sqrt(norm2);

                if (acos(proj) > M_PI/2.0)
                {
                   for (k = 0; k < dftb1[i].norb; k++)
                   {
                       dftb1[i].a[k][l] = -dftb1[i].a[k][l];
                   }
                   printf("Change the sign of initial coupling: %f  \n", acos(proj));
                }

            }

            for (j = 0; j < dftb1[i].norb; j++)
            {
                for (k = 0; k < dftb1[i].norb; k++)
                {
                    dftb1[i].a_old[k][j] = dftb1[i].a[k][j];
                }
            }

        }

        // calc overlap with previous step //
        for (j = 0; j < ct->site[i].homos; j++)
        {
            for (k = 0; k < ct->site[i].homos; k++)
            {
                l = ct->site[i].homo[j]-1;
                m = ct->site[i].homo[k]-1;
                ct->site[i].overlap[j][k] = 0.0;
                for (iao = 0; iao < dftb1[i].norb; iao++)
                {
                    for (jao = 0; jao < dftb1[i].norb; jao++)
                    {
                        ct->site[i].overlap[j][k] += dftb1[i].a[iao][l] * dftb1[i].overl[iao][jao] * dftb1[i].a_old[jao][m]; //is later also needed by project_wf_on_new_basis() if FOs are degenerated.
                    }
                }
            }
        }

        // invert sign if needed
        for (j = 0; j < ct->site[i].homos; j++)
        {
            if (ct->site[i].overlap[j][j] < -0.0)
            {
                l = ct->site[i].homo[j]-1;
                for (k = 0; k < dftb1[i].norb; k++)
                {
                    dftb1[i].a[k][l] = -dftb1[i].a[k][l];
                }
            }
            if (ct->site[i].overlap[j][j] > -0.9 && ct->site[i].overlap[j][j] < 0.9)
            {
                printf("warning for site %d: strong change of shape for orbital %d between two steps! overl = %4.3f \n", i, j,  ct->site[i].overlap[j][j]);
                //  exit(-1);
            }
        }

        //calculate correct overlap after inversion of signs
        for (j = 0; j < ct->site[i].homos; j++)
        {
            for (k = 0; k < ct->site[i].homos; k++)
            {
                l = ct->site[i].homo[j]-1;
                m = ct->site[i].homo[k]-1;
                ct->site[i].overlap[j][k] = 0.0;
                for (iao = 0; iao < dftb1[i].norb; iao++)
                {
                    for (jao = 0; jao < dftb1[i].norb; jao++)
                    {
                        ct->site[i].overlap[j][k] += dftb1[i].a[iao][l] * dftb1[i].overl[iao][jao] * dftb1[i].a_old[jao][m];
                    }
                }
            }
        }
        /* update the "old" array */
        for (j = 0; j < dftb1[i].norb; j++)
        {
            for (k = 0; k < dftb1[i].norb; k++)
            {
                dftb1[i].a_old[k][j] = dftb1[i].a[k][j];
            }
        }
    }


    if (ct->first_step && ct->define_orbital_sign)
    {
       sfree(orb_index);

       for (i = 0; i < ct->sites; i++)
           for (j = 0; j < 3; j++)
               sfree(coords[i][j]);
    }

    return;
}


/*
   void check_and_invert_orbital_phase(dftb_phase1_t *dftb1, int sites)
   {
   int i, j, k, l, best_one;
   double dotprod, best_dot, temp;

   for (i=0; i<sites; i++) {
    // look at orbitals in site i
    for (j=0; j<dftb1[i].norb; j++) {
      // j-th orbital
      dotprod = 0.;
      // calculate the dot product of j-th orbital with j-th orbital in the previous step
      for (k=0; k<dftb1[i].norb; k++)
        dotprod += dftb1[i].a_old[k][j] * dftb1[i].a[k][j];
       check the dot product:
 * if negative: invert a[..][j]
 * if too small: look for orbitals with better overlap (maybe the orbital ordering is altered)

   //      if (dotprod >= 0.5) break; // normally this should be the case, so no further calculations are required
      if (dotprod <= -0.5) // invert
        for (k=0; k<dftb1[i].norb; k++)
          dftb1[i].a[k][j] = - dftb1[i].a[k][j];
      if (-0.5 < dotprod && dotprod < 0.5){ //find best overlap
   //        best_dot = dotprod;
        best_one = j;
        for (l=0; l<dftb1[i].norb; l++){
          dotprod = 0.;
          for (k=0; k<dftb1[i].norb; k++)
            dotprod += dftb1[i].a_old[k][l] * dftb1[i].a[k][j];
          if ( fabs(dotprod) > fabs(best_dot) ){
            best_one = l;
            best_dot = dotprod;
          }
        }
         // swap orbital j with best_one
        for (k=0; k<dftb1[i].norb; k++){
          temp = dftb1[i].a[k][j];
          dftb1[i].a[k][j] = ( best_dot > 0 ) ? dftb1[i].a[k][best_one] : -dftb1[i].a[k][best_one]; // invert if needed
          dftb1[i].a[k][best_one] = temp;
          //8888888888888888888888888888888888888888 ich glaube die orbitalenergien aus phase1  müssen nicht vertauscht werden, da sie in phase2 aus den H in atomistischer basis durch transformation (abhängig von den MOs) erhalten werden
        }
        printf("swapped orbital %d and %d of site %d. dotproduct is %f \n", j, best_one, i, best_dot);
    //  }
      // update the "old" array
      for (k=0; k<dftb1[i].norb; k++)
        dftb1[i].a_old[k][j] = dftb1[i].a[k][j];
    }
   }

   return;
   }
 */

