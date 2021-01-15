#include "gromacs/transfer/transfer.h"

extern "C" void dsyev_(char *, char *, long *, double *, long *, double *, double *, long *, long *);
extern "C" void dgesv_(long *, long *, double *, long *, long *, double *, long *, long *);

long dsyev(int int_n, double *a, double *w, double *work, int int_lwork)
{
    long info, n, lwork;
    char jobz = 'V', uplo = 'U';
    n     = (long) int_n;
    lwork = (long) int_lwork;

    dsyev_(&jobz, &uplo, &n, a, &n, w, work, &lwork, &info);
    if ((int) info)
    {
        printf("\nDSYEV: ier = %d\nEXITING!\n\n", (int) info);
        exit(-1);
    }

    return info;
}

long dgesv(int int_n, double *a, double *b, long ipiv[DIIS_MAX_PREV_VECTORS + 1])
{
    long n, nrhs, lda, ldb, info;
    n    = lda = ldb = (long) int_n;
    nrhs = 1L;

    dgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);

    return info;
}

