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


void sort_mobasis(dftb_t *dftb, charge_transfer_t *ct, int i)
{
    int           j, k, l, m, n, iao, jao;
    double        dummy;
    dftb_phase1_t dftb1;

    dftb1 = dftb->phase1[i];

    if (ct->first_step)
    {
        for (m = 0; m < dftb1.norb; m++)
        {
            for (n = 0; n < dftb1.norb; n++)
            {
                dftb1.a_old[m][n] = dftb1.a[m][n];
                dftb1.a_ref[m][n] = dftb1.a[m][n];
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
            for (iao = 0; iao < dftb1.norb; iao++)
            {
                for (jao = 0; jao < dftb1.norb; jao++)
                {
                    ct->site[i].overlap[j][k] += dftb1.a_old[iao][l] * dftb1.overl[iao][jao] * dftb1.a[jao][m];
                }
            }
        }
    }
    // swap orbitals //
/*
   for (j = 0; j < ct->site[i].homos; j++)
   for (k = 0; k < ct->site[i].homos; k++){
    l=ct->site[i].homo[j]-1;
    m=ct->site[i].homo[k]-1;
    if (fabs(ct->site[i].overlap[j][k]) > fabs(ct->site[i].overlap[j][j])){
      printf("orbital swap of orb %d and %d! overlap %lf > %lf \n", ct->site[i].homo[j], ct->site[i].homo[k], ct->site[i].overlap[j][k] , ct->site[i].overlap[j][j]);
      for (iao=0; iao<dftb1.norb; iao++){
        dummy = dftb1.a[iao][l]; //save vector in dummy during swap
        dftb1.a[iao][l] = dftb1.a[iao][m];
        dftb1.a[iao][m] = dummy;
      }
      dummy = ct->site[i].overlap[j][k];
      ct->site[i].overlap[j][j] = ct->site[i].overlap[j][k];
      ct->site[i].overlap[j][k] = dummy;
      dummy = dftb1.ev[l];
      dftb1.ev[l] = dftb1.ev[m];
      dftb1.ev[m] = dummy;
    }
   }
 */

    // check orbital phase //
    for (j = 0; j < ct->site[i].homos; j++)
    {
        if (ct->site[i].overlap[j][j] < -0.)
        {
            l = ct->site[i].homo[j]-1;
            for (iao = 0; iao < dftb1.norb; iao++)
            {
                dftb1.a[iao][l] = -dftb1.a[iao][l];
            }
        }
    }

    // recalculate overlap
    for (j = 0; j < ct->site[i].homos; j++) //only sort HOMOs
    {
        for (k = 0; k < ct->site[i].homos; k++)
        {
            l = ct->site[i].homo[j]-1;
            m = ct->site[i].homo[k]-1;
            ct->site[i].overlap[j][k]     = 0.0;
            ct->site[i].overlap_ref[j][k] = 0.0;
            for (iao = 0; iao < dftb1.norb; iao++)
            {
                for (jao = 0; jao < dftb1.norb; jao++)
                {
                    ct->site[i].overlap[j][k]     += dftb1.a_old[iao][l] * dftb1.overl[iao][jao] * dftb1.a[jao][m];
                    ct->site[i].overlap_ref[j][k] += dftb1.a_ref[iao][l] * dftb1.overl[iao][jao] * dftb1.a[jao][m];
                }
            }
            ct->site[i].overlap[k][j]     = ct->site[i].overlap[j][k];
            ct->site[i].overlap_ref[k][j] = ct->site[i].overlap_ref[j][k];
        }
    }

    // print quality of overlap //
    for (j = 0; j < ct->site[i].homos; j++)
    {
        if (ct->site[i].overlap[j][j] > -0.9 && ct->site[i].overlap[j][j] < 0.9)
        {
            printf("warning: strong change of shape for orbital %d between two steps! overl = %4.3f \n", j, ct->site[i].overlap[j][j]);
            //  exit(-1);
        }
    }

    /* update the "old" array */
    for (j = 0; j < dftb1.norb; j++)
    {
        for (k = 0; k < dftb1.norb; k++)
        {
            dftb1.a_old[k][j] = dftb1.a[k][j];
        }
    }

    return;
}
    

void get_neighbor_list(charge_transfer_t *ct, dftb_t *dftb)
{
  int counter1, counter2;
  double dist;
  dvec bond;
  // get neighbor list
  counter1=-1;
  for (int i=0; i<ct->sites; i++)
  for (int ii=0; ii<ct->site[i].atoms; ii++){
    counter1++;
    counter2=-1;
    for (int j=0; j<ct->sites; j++)
    for (int jj=0; jj<ct->site[j].atoms; jj++){
      counter2++;
      if (counter2 > counter1) {break;}
      dvec_sub(dftb->phase1[i].x[ii],dftb->phase1[j].x[jj], bond);
      dist = dnorm(bond);
      if (dist < 20.0){
        dftb->nl[counter1][counter2]=dftb->nl[counter2][counter1]=1;
      }else{
        dftb->nl[counter1][counter2]=dftb->nl[counter2][counter1]=0;
      }
    }
  }
}


int find_union(int size, int array1[], int array2[], int union_array[])
{
    int i = 0, j = 0, k = 0;
    while ((i < size) && (j < size))
    {
        if (array1[i] < array2[j])
        {
            union_array[k] = array1[i];
            i++;
            k++;
        }
        else if (array1[i] > array2[j])
        {
            union_array[k] = array2[j];
            j++;
            k++;
        }
        else
        {
            union_array[k] = array1[i];
            i++;
            j++;
            k++;
        }
    }
    if (i == size)
    {
        while (j < size)
        {
            union_array[k] = array2[j];
            j++;
            k++;
        }
    }
    else
    {
        while (i < size)
        {
            union_array[k] = array1[i];
            i++;
            k++;
        }
    }
    return(k);
}

