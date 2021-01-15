/* From the old DFTB codebase
 * converted to C by Tomas Kubar.
 */

#include "gromacs/transfer/transfer.h"

extern "C" void dsytrf_(char *, long *, double *, long *, long *, double *, long *, long *);
extern "C" void dsytri_(char *, long *, double *, long *, long *, double *, long *);

static int inverse(int isize, double *b, double *work, long *ipiv)
{
  static long info, imatsz=IMATSZ_BROYDEN, size;
  static char uplo='U';
  size = (long) isize;

  dsytrf_(&uplo, &size, b, &imatsz, ipiv, work, &imatsz, &info);
  // printf("dsytrf info: %ld\n", info);
  if (info) return info;
  dsytri_(&uplo, &size, b, &imatsz, ipiv, work, &info);
  // printf("dsytri info: %ld\n", info);
  return info;
} 

void broyden(int niter, double alpha, int jtop, double *vecin, double *vecout, dftb_broyden_t *arrays)
{
  double w0, amix, dfnorm, fnorm, fac1, fac2, aij, cmj, gmi, wtmp;
  int lastit, iter, i, j, k, lm, ln, ip, ivsiz;
  static double uamix;
  static int ilastit;

  // printf("START BROYDEN - iter %d\n", niter);
  // for (k=0; k<jtop; k++) {
  //   printf("%3d%10.6f%10.6f\n", k+1, vecin[k], vecout[k]);
  // }

  iter = niter;
  if (niter >= MAXITER_BROYDEN)
    iter = iter % MAXITER_BROYDEN;
  /* if (iter == 0)
    return; */

  // ++++++ SET UP THE VECTOR OF THE CURRENT ITERATION FOR MIXING ++++++
  // 
  // FOR THIS METHOD WE HAVE ONLY SAVED INPUT/OUTPUT CHG. DENSITIES,
  for (k=0; k<jtop; k++) {
    arrays->vector[k][0] = vecin[k];
    arrays->vector[k][1] = vecout[k];
  }
  // ++++++ END OF PROGRAM SPECIFIC LOADING OF VECTOR FROM MAIN ++++++++
  // 
  // IVSIZ IS THE LENGTH OF THE VECTOR
 
  ivsiz = jtop;

  // 
  // 
  // *******************  BEGIN BROYDEN'S METHOD  **********************
  // 
  //   WEIGHTING FACTOR FOR THE ZEROTH ITERATION
  w0 = 0.01;

  //     F:  THE DIFFERENCE OF PREVIOUS OUTPUT AND INPUT VECTORS
  // DUMVI:  A DUMMY VECTOR, HERE IT IS THE PREVIOUS INPUT VECTOR
  if (iter == 0) {
    // IF THIS IS THE FIRST ITERATION, THEN LOAD
    // F=VECTOR(OUT)-VECTOR(IN) AND VECTOR(IN)
    // PRINT*,'SIMPLE MIXING THIS ITERATION'
    lastit = 0;
    amix = alpha;
    uamix = amix;
    ilastit = lastit;
    for (k=0; k<ivsiz; k++) {
      arrays->f[k] = arrays->vector[k][1] - arrays->vector[k][0];
      arrays->unit31[k][0] = arrays->f[k];
      arrays->unit31[k][1] = arrays->vector[k][0];
      // SINCE WE ARE ON THE FIRST ITERATION, SIMPLY MIX THE VECTOR.
      arrays->dumvi[k] = arrays->vector[k][0] + amix * arrays->f[k];
    }
  } else {
    // ALL THIS IS PERFORMED IN FURTHER ITERATIONS
    amix = uamix;
    lastit = ilastit;
    // FOR I-TH ITER.,DFNORM IS ( F(I) MINUS F(I-1) ), USED FOR NORMALIZATION
    dfnorm = 0.0;
    fnorm = 0.0;
    for (k=0; k<ivsiz; k++) {
       // ALPHA (OR AMIX) IS SIMPLE MIXING PARAMETER
      arrays->dumvi[k] = arrays->vector[k][0] - arrays->unit31[k][1];
      arrays->df[k] = arrays->vector[k][1] - arrays->vector[k][0] - arrays->unit31[k][0];
      arrays->f[k] = arrays->vector[k][1] - arrays->vector[k][0];
      dfnorm += arrays->df[k] * arrays->df[k];
      fnorm += arrays->f[k] * arrays->f[k];
    }
    dfnorm = sqrt(dfnorm);
    fnorm = sqrt(fnorm);

    fac2 = 1.0 / dfnorm;
    fac1 = amix * fac2;

    for (k=0; k<ivsiz; k++) {
      arrays->ui[k] = fac1 * arrays->df[k] + fac2 * arrays->dumvi[k];
      arrays->vti[k] = fac2 * arrays->df[k];
    }

    // *********** CALCULATION OF COEFFICIENT MATRICES *************
    // ***********    AND THE SUM FOR CORRECTIONS      *************
    // 
    // RECALL: A(I,J) IS A SYMMETRIC MATRIX
    //       : B(I,J) IS THE INVERSE OF [ W0**2 I + A ]
    //
    lastit++;

    // DUMVI IS THE U(OF I) AND T1 IS THE VT(OF I)
    // FROM THE PREVIOUS ITERATIONS
    if (lastit > 1) {
      for (j=0; j<lastit-1; j++) {
        aij = cmj = 0.0;
	for (k=0; k<ivsiz; k++) {
	  arrays->dumvi[k] = arrays->unit32[k][0][j];
	  arrays->t1[k] = arrays->unit32[k][1][j];
	  cmj += arrays->t1[k] * arrays->f[k];
	  aij += arrays->t1[k] * arrays->vti[k];
	}
	arrays->a[lastit-1][j] = arrays->a[j][lastit-1] = aij;
	//printf("a[%d][%d] = %f\n", j, lastit-1, aij);
	arrays->cm[j] = cmj;
      }
    }

    aij = cmj = 0.0;
    for (k=0; k<ivsiz; k++) {
      cmj += arrays->vti[k] * arrays->f[k];
      aij += arrays->vti[k] * arrays->vti[k];
    }
    arrays->a[lastit-1][lastit-1] = aij;
    //printf("a[%d][%d] = %f\n", lastit-1, lastit-1, aij);
    arrays->cm[lastit-1] = cmj;

    for (k=0; k<ivsiz; k++) {
      arrays->unit32[k][0][lastit-1] = arrays->ui[k];
      arrays->unit32[k][1][lastit-1] = arrays->vti[k];
    }

    // THE WEIGHTING FACTORS FOR EACH ITERATION HAVE BEEN CHOSEN
    // EQUAL TO ONE OVER THE R.M.S. ERROR. THIS NEED NOT BE THE CASE.
    if (fnorm > 1.0e-7)
      wtmp = 0.01 / fnorm;
    else
      wtmp = 1.0e5;
    if (wtmp < 1.0)
      wtmp = 1.0;
    arrays->w[lastit-1] = wtmp;
    // WRITE(66,'(''  WEIGHTING SET =  '',E12.6)')WTMP

    // WITH THE CURRENT ITERATIONS F AND VECTOR CALCULATED,
    // WRITE THEM TO UNIT 31 FOR USE LATER.
    uamix = amix;
    ilastit = lastit;
    for (k=0; k<ivsiz; k++) {
      arrays->unit31[k][0] = arrays->f[k];
      arrays->unit31[k][1] = arrays->vector[k][0];
    }

    // SET UP AND CALCULATE BETA MATRIX
    for (lm=0; lm < lastit; lm++) {
      for (ln=0; ln < lastit; ln++)
        arrays->b_lapack[ln * IMATSZ_BROYDEN + lm] 
		= arrays->b[ln][lm] 
		= arrays->a[ln][lm] * arrays->w[ln] * arrays->w[lm];
      arrays->b_lapack[lm * IMATSZ_BROYDEN + lm] 
	      = arrays->b[lm][lm] 
	      = w0*w0 + arrays->a[lm][lm] * arrays->w[lm] * arrays->w[lm];
    }

    // print out the matrix
    // printf("w vector:\n");
    //   for (ln=0; ln<lastit; ln++)
    //     printf("%9.5f", arrays->w[ln]);
    //   printf("\n");
    // printf("beta matrix:\n");
    // for (lm=0; lm<lastit; lm++) {
    //   for (ln=0; ln<lastit; ln++)
    //     printf("%9.5f", arrays->b[lm][ln]);
    //   printf("\n");
    // }

    // INVERT THE MATRIX USING LAPACK, INSTEAD OF CALL INVERSE(D,B,LASTM1)
    if (inverse(lastit, arrays->b_lapack, arrays->work, arrays->ipiv) != 0) {
      printf("Broyden: error in matrix inversion\n");
      exit(-1);
    }

    for (lm=0; lm<lastit; lm++)
      for (ln=0; ln<lastit; ln++)
        arrays->b[lm][ln] = lm < ln ? 
		arrays->b_lapack[ln * IMATSZ_BROYDEN + lm]:
	       	arrays->b_lapack[lm * IMATSZ_BROYDEN + ln];

    // print out the inverse
    // printf("inverse of beta matrix:\n");
    // for (lm=0; lm<lastit; lm++) {
    //   for (ln=0; ln<lastit; ln++)
    //     printf("%9.5f", arrays->b[lm][ln]);
    //   printf("\n");
    // }

    // calculate the vector for the new iteration
    for (k=0; k<ivsiz; k++)
      arrays->dumvi[k] = arrays->vector[k][0] + amix * arrays->f[k];

    for (i=0; i<lastit; i++) {
      for (k=0; k<ivsiz; k++) {
        arrays->ui[k] = arrays->unit32[k][0][i];
        arrays->vti[k] = arrays->unit32[k][1][i];
      }
      gmi = 0.0;
      for (ip=0; ip<lastit; ip++)
        gmi += arrays->cm[ip] * arrays->b[ip][i] * arrays->w[ip];
      for (k=0; k<ivsiz; k++) {
        arrays->dumvi[k] -= gmi * arrays->ui[k] * arrays->w[i];
	// printf("%9.5f%9.5f%9.5f%9.5f\n", arrays->dumvi[k], gmi, arrays->ui[k], arrays->w[i]);
      }
    }
  }
  //     END OF THE CALCULATION OF DUMVI, THE NEW VECTOR
  //
  //**********  THE END OF THE BROYDEN METHOD **************
  //
  //+++++ PROGRAM SPECIFIC CODE OF RELOADING ARRAYS +++++++++
  //
  // NEED TO UNLOAD THE NEW VECTOR INTO THE APPROPRIATE ARRAYS.
  for (k=0; k<jtop; k++)
    vecin[k] = arrays->dumvi[k];
  // printf("END BROYDEN - iter %d\n", niter);
  // for (k=0; k<jtop; k++) {
  //   printf("%3d%10.6f%10.6f\n", k+1, vecin[k], vecout[k]);
  // }
  return;
}
