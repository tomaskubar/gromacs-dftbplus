#include "gromacs/transfer/transfer.h"

extern "C" void dsyevr_(char *, char *, char *, long *, double *, long *,
                        double *, double *, long *, long *, double *, long *, double *,
                        double *, long *, double *, double *, long *, double *, long *, long *);

void step1(int *step, int nstate, double **cm);
void step2(int *step, int nstate, double **cm, int **star, int *rowc, int *colc);
void step3(int *step, int nstate, int *colc, int **star);
void step4(int *step, int nstate, double **cm, int **star, int *colc, int *rowc, int *path_row_0, int *path_col_0);
void step5(int *step, int nstate, int *path_count, int **star, int *rowc, int *colc, int **path, int path_row_0, int path_col_0);
void step6(int *step, int nstate, double **cm, int *colc, int *rowc);
void step7(int nstate, int **star, int *ind);
void RK5th(charge_transfer_t *ct);
void find_a_zero(int nstate, double **cm, int *rowc, int *colc, int *irow, int *icol);
int star_in_row(int nstate, int irow, int **star);
void find_star_in_row(int nstate, int irow, int **star, int *icol);
void rk_delta_c(double *c0, double *dc, charge_transfer_t *ct);


double fermi_coef_sum(double x, int n, double *a, double *coeffs, double fermi_kt)
{
    int    i;
    double sum;

    sum = 0.0;
    for (i = 0; i < n; i++)
    {
        sum += coeffs[i] = 1.0 / ( exp( (a[i]-x) / fermi_kt ) + 1.0 );
    }

    return sum;
}

/* OLD VERISON
   double fermi_coef_sum(double x, int n, double *a, double *coeffs)
   {
   int i;
   double sum;

   sum = 0.0;
   for (i=0; i<n; i++)
    sum += coeffs[i] = 1.0 / ( exp( (a[i]-x) / FERMI_KT ) + 1.0 );

   return sum;
   }
 */


long exp_imag_matrix(double **in, twodoubles **out, double dt, long n, ct_per_orthogo_t *arrays)
{
    double abstol = 1.0e-8;
 // char   itype = 'U';
    char   jobz = 'V';
    char   range = 'A';
    char   uplo = 'U';
    long   i, j, k, info, eval_no;
    double void_double;
    long   void_long;

    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
        {
            arrays->in[i+j*n] = in[i][j];
        }

/*
   // print the arrays
   printf("Orthogonalize - array sij:\n");
   for (i=0; i<n; i++) {
    for (j=0; j<n; j++) printf("%12.6f", arrays.sij[i+j*n]);
    printf("\n");
   }
   printf("Orthogonalize - array tij:\n");
   for (i=0; i<n; i++) {
    for (j=0; j<n; j++) printf("%12.6f", arrays.tij[i+j*n]);
    printf("\n");
   }
 */

    // Diagonalize the matrix
    // evec - eigenvectors, eval - eigenvalues
    dsyevr_(&jobz, &range, &uplo, &n, arrays->in, &n,
            &void_double, &void_double, &void_long, &void_long, &abstol, &eval_no, arrays->eval,
            arrays->evec, &n, arrays->issupz, arrays->work, &(arrays->lwork), arrays->iwork, &(arrays->liwork), &info);
    //if (info) return info;

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            out[i][j][0] = 0.;
            out[i][j][1] = 0.;
            for (k = 0; k < n; k++)
            {
                out[i][j][0] += arrays->evec[i + n * k] * arrays->evec[j + n * k] * cos(arrays->eval[k] * dt);
                out[i][j][1] -= arrays->evec[i + n * k] * arrays->evec[j + n * k] * sin(arrays->eval[k] * dt);
            }
        }
    }

    return info;
}


void assignprob(charge_transfer_t *ct, int *ind)
{

    int     i, j, k, step, nstate, done, path_row_0, path_col_0, path_count;
 // int     nsmax = 2;
    int     *rowc, *colc, *map, **star, **path;
    double  **cm;

    snew(cm, ct->dim);
    snew(rowc, ct->dim);
    snew(colc, ct->dim);
    snew(star, ct->dim);
    snew(path, 2*ct->dim+1);
    snew(map, ct->dim);

    snew(cm[0], SQR(ct->dim));
    snew(star[0], SQR(ct->dim));
    for (i = 0; i < 2*ct->dim+1; i++)
        snew(path[i], 2);

    for (j = 1; j < ct->dim; j++)
    {
        cm[j] = cm[0] + j*ct->dim;
        star[j] = star[0] + j*ct->dim;
    }

    for (i = 0; i < ct->dim; i++)
    {
        colc[i] = 0; rowc[i] = 0;
        for (j = 0; j < ct->dim; j++)
        {
            cm[i][j] = - SQR(ct->tfs_overlap[j][i]);
            star[i][j] = 0;
        }
    }

/*    for (i = 0; i < ct->dim; i++)
        for (j = 0; j < ct->dim; j++)
            if (j < i-nsmax || j > i+nsmax)
               cm[i][j] = 1000.0;
*/
    for (i = 0; i < ct->dim; i++)
    {
        map[i] = 0; ind[i] = -1;
    }

    nstate = ct->dim;
    step = 1;
    done = 0;
    while (!done)
    {
      switch (step)
      {
         case 1:
              step1(&step, nstate, cm);
              break;
         case 2:
              step2(&step, nstate, cm, star, rowc, colc);
              break;
         case 3:
              step3(&step, nstate, colc, star);
              break;
         case 4:
              step4(&step, nstate, cm, star, colc, rowc, &path_row_0, &path_col_0);
              break;
         case 5:
              step5(&step, nstate, &path_count, star, rowc, colc, path, path_row_0, path_col_0);
              break;
         case 6:
              step6(&step, nstate, cm, colc, rowc);
              break;
         case 7:
              step7(nstate, star, ind);
              done = 1;
              break;
      }
    }

// check the mapping states
    for (i = 0; i < nstate; i++)
    {
        for (j = 0; j < nstate; j++)
        {
            if (ind[i] == j)
            {
               map[j] = 1;
               break;
            }
        }
    }

// reorder the unmapped states
    k = - 1;
    for (i = 0; i < nstate; i++)
    {
        if (ind[i] == -1)
        {
           for (j = k + 1; j < nstate; j++)
           {
               if (! map[j])
               {
                  ind[i] = j; map[j] = 1;
                  printf("State %d is not ordered! \n", j);
                  k = j;
                  break;
               }
           }
        }
    }

   return;
}
//For each row of the cost matrix, find the smallest element and subtract
//it from every element in its row.  When finished, Go to Step 2.
void step1(int *step, int nstate, double **cm)
{
     int    i, j;
     double minval;

     for (i = 0; i < nstate; i++)
     {
         minval = cm[i][0];
         for (j = 0; j < nstate; j++)
             if (minval > cm[i][j])
                minval = cm[i][j];

         for (j = 0; j < nstate; j++)
             cm[i][j] -= minval;
     }

     *step = 2;
}
//Find a zero (Z) in the resulting matrix.  If there is no starred
//zero in its row or column, star Z. Repeat for each element in the
//matrix. Go to Step 3.
void step2(int *step, int nstate, double **cm, int **star, int *rowc, int *colc)
{

     int i, j;


     for (i = 0; i < nstate; i++)
     {
         for (j = 0; j < nstate; j++)
         {
             if (cm[i][j] < 0.000000000000001 && rowc[i] == 0 && colc[j] == 0)
             {
                star[i][j] = 1;
                rowc[i] = 1;
                colc[j] = 1;
             }
         }
    }
     for (i = 0; i < nstate; i++)
     {
         rowc[i] = 0;
         colc[i] = 0;
     }

     *step = 3;
}
//Cover each column containing a starred zero.  If K columns are covered,
//the starred zeros describe a complete set of unique assignments.  In this
//case, Go to DONE, otherwise, Go to Step 4.
void step3(int *step, int nstate, int *colc, int **star)
{
     int i, j, k;

     for (i = 0; i < nstate; i++)
         for (j = 0; j < nstate; j++)
             if (star[i][j] == 1)
                colc[j] = 1;

     k = 0;
     for (j = 0; j < nstate; j++)
         if (colc[j] == 1)
            k += 1;

     if (k >= nstate)
        *step = 7;
     else
        *step = 4;
}
//Find a noncovered zero and prime it.  If there is no starred zero
//in the row containing this primed zero, Go to Step 5.  Otherwise,
//cover this row and uncover the column containing the starred zero.
//Continue in this manner until there are no uncovered zeros left.
//Save the smallest uncovered value and Go to Step 6.
void step4(int *step, int nstate, double **cm, int **star, int *colc, int *rowc, int *path_row_0, int *path_col_0)
{
     int irow = -1, icol = -1, done;

     done = 0;
     while (!done)
     {
         find_a_zero(nstate, cm, rowc, colc, &irow, &icol);
         if (irow == -1)
         {
             done = 1;
             *step = 6;
         }
         else
         {
             star[irow][icol] = 2;
             if (star_in_row(nstate, irow, star))
             {
                 find_star_in_row(nstate, irow, star, &icol);
                 rowc[irow] = 1;
                 colc[icol] = 0;
             }
             else
             {
                 done = 1;
                 *step = 5;
                 *path_row_0 = irow;
                 *path_col_0 = icol;
             }
         }
     }
}
//methods to support step 4
void find_a_zero(int nstate, double **cm, int *rowc, int *colc, int *irow, int *icol)
{
     int i = 0, j, done;

     *irow = -1;
     *icol = -1;
     done = 0;
     while (!done)
     {
         j = 0;
         while (1)
         {
             if (cm[i][j] < 0.000000000000001 && rowc[i] == 0 && colc[j] == 0)
             {
                 *irow = i;
                 *icol = j;
                 done = 1;
             }
             j += 1;
             if (j >= nstate || done)
                 break;
         }
         i += 1;
         if (i >= nstate)
            done = 1;
     }
}

int star_in_row(int nstate, int irow, int **star)
{
    int j, tmp;

    tmp = 0;
    for (j = 0; j < nstate; j++)
        if (star[irow][j] == 1)
            tmp = 1;

    return tmp;
}

void find_star_in_row(int nstate, int irow, int **star, int *icol)
{
    int j;

    *icol = -1;
    for (j = 0; j < nstate; j++)
        if (star[irow][j] == 1)
            *icol = j;
}
//Construct a series of alternating primed and starred zeros as follows.
//Let Z0 represent the uncovered primed zero found in Step 4.  Let Z1 denote
//the starred zero in the column of Z0 (if any). Let Z2 denote the primed zero
//in the row of Z1 (there will always be one).  Continue until the series
//terminates at a primed zero that has no starred zero in its column.
//Unstar each starred zero of the series, star each primed zero of the series,
//erase all primes and uncover every line in the matrix.  Return to Step 3.
void step5(int *step, int nstate, int *path_count, int **star, int *rowc, int *colc, int **path, int path_row_0, int path_col_0)
{
    int i, j, done;
    int irow = -1;
    int icol = -1;

    *path_count = 1;
    path[*path_count - 1][0] = path_row_0;
    path[*path_count - 1][1] = path_col_0;
    done = 0;
    while (!done)
    {
        irow = -1;
        for (i = 0; i < nstate; i++)
            if (star[i][path[*path_count - 1][1]] == 1)
               irow = i;

        if (irow > -1)
        {
            *path_count += 1;
            path[*path_count - 1][0] = irow;
            path[*path_count - 1][1] = path[*path_count - 2][1];
        }
        else
            done = 1;
        if (!done)
        {
           for (j = 0; j < nstate; j++)
               if (star[path[*path_count - 1][0]][j] == 2)
                  icol = j;

            *path_count += 1;
            path[*path_count - 1][0] = path[*path_count - 2][0];
            path[*path_count - 1][1] = icol;
        }
    }

    for (i = 0; i < *path_count; i++)
        if (star[path[i][0]][path[i][1]] == 1)
            star[path[i][0]][path[i][1]] = 0;
        else
            star[path[i][0]][path[i][1]] = 1;

    for (i = 0; i < nstate; i++)
    {
        rowc[i] = 0;
        colc[i] = 0;
    }

    for (i = 0; i < nstate; i++)
        for (j = 0; j < nstate; j++)
            if (star[i][j] == 2)
                star[i][j] = 0;
    *step = 3;
}

//Add the value found in Step 4 to every element of each covered row, and subtract
//it from every element of each uncovered column.  Return to Step 4 without
//altering any stars, primes, or covered lines.
void step6(int *step, int nstate, double **cm, int *colc, int *rowc)
{
    int     i, j;
    double  minval;

    minval = 1e10;
    for (i = 0; i < nstate; i++)
        for (j = 0; j < nstate; j++)
            if (rowc[i] == 0 && colc[j] == 0)
                if (minval > cm[i][j])
                    minval = cm[i][j];

    for (i = 0; i < nstate; i++)
        for (j = 0; j < nstate; j++)
        {
            if (rowc[i] == 1)
                cm[i][j] += minval;
            if (colc[j] == 0)
                cm[i][j] -= minval;
        }
    *step = 4;
}

void step7(int nstate, int **star, int *ind)
{

     int i, j;

     for (i = 0; i < nstate; i++)
         for (j = 0; j < nstate; j++)
             if (star[i][j] == 1) ind[i] = j;

     return;
}

void RK5th(charge_transfer_t *ct)
{
    int i;
    double  *c0, *dc1, *dc2, *dc3, *dc4, *dc5, *dc6;

    snew(c0, SQR(ct->dim));
    snew(dc1, SQR(ct->dim));
    snew(dc2, SQR(ct->dim));
    snew(dc3, SQR(ct->dim));
    snew(dc4, SQR(ct->dim));
    snew(dc5, SQR(ct->dim));
    snew(dc6, SQR(ct->dim));


    for (i = 0; i < SQR(ct->dim); i++)
        c0[i] = ct->tfs_diab[i];

    rk_delta_c(c0, dc1, ct);
    for (i = 0; i < SQR(ct->dim); i++)
        c0[i] = ct->tfs_diab[i] + ct->rk_timestep/4.0*dc1[i];

    rk_delta_c(c0, dc2, ct);
    for (i = 0; i < SQR(ct->dim); i++)
        c0[i] = ct->tfs_diab[i] + ct->rk_timestep/4.0*dc2[i];

    rk_delta_c(c0, dc3, ct);
    for (i = 0; i < SQR(ct->dim); i++)
        c0[i] = ct->tfs_diab[i] + ct->rk_timestep/2.0*dc3[i];

    rk_delta_c(c0, dc4, ct);
    for (i = 0; i < SQR(ct->dim); i++)
        c0[i] = ct->tfs_diab[i] + ct->rk_timestep*3.0/4.0*dc4[i];

    rk_delta_c(c0, dc5, ct);
    for (i = 0; i < SQR(ct->dim); i++)
        c0[i] = ct->tfs_diab[i] + ct->rk_timestep*dc5[i];

    rk_delta_c(c0, dc6, ct);

    for (i = 0; i < SQR(ct->dim); i++)
    {
        dc1[i] += dc6[i];
        dc3[i] += dc5[i];
    }

    for (i = 0; i < SQR(ct->dim); i++)
        ct->tfs_diab[i] += 1.0/90.0*(7.0*dc1[i]+32.0*dc3[i]+12.0*dc4[i])*ct->rk_timestep;
}

void rk_delta_c(double *c0, double *dc, charge_transfer_t *ct)
{

   int i, j;
   double tau;

   for (i = 0; i < ct->dim; i++)
   {
       dc[i] = 0.e0;
       dc[ct->dim+i] = 0.e0;

       for (j = 0; j < ct->dim; j++)
       {
           dc[i]         += ct->tfl_mean_ham_full[i][j]*c0[ct->dim+j];// - ct->fo_tdc[i][j]*c0[j];
           dc[ct->dim+i] -= ct->tfl_mean_ham_full[i][j]*c0[j];// - ct->fo_tdc[i][j]*c0[ct->dim+j];
       }
   }

 // add decoherence corrections to ehrenfest method, weiwei xie, Feb 2019
   if (ct->jobtype == cteCPFSCCDYNAMIC || ct->jobtype == cteDDBSCCDYNAMIC)
   {
      for (i = 0; i < ct->dim; i++)
          ct->occupation[i] = SQR(c0[i]) + SQR(c0[i+ct->dim]);

      for (i = 0; i < ct->dim; i++)
          for (j = 0; j < ct->dim; j++)
              if (i != j)
              {
                 tau = (1.0 + 0.1/ct->neg_imag_pot_rate_constant)/fabs(ct->tfl_mean_ham_full[i][i] - ct->tfl_mean_ham_full[j][j]);
                 dc[i]         += 4.0/tau*ct->occupation[j]*c0[ct->dim+i];
                 dc[ct->dim+i] -= 4.0/tau*ct->occupation[j]*c0[i];
              }
   }
}

