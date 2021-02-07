#include "gromacs/transfer/transfer.h"

long negf_const_mat_(double *q, double *w, long *e, long *r)
{
    (void) q;
    (void) w;
    (void) e;
    (void) r;
    return 1l;
}

long negf_jacobi_(double *q, double *w, double *e, long *r)
{
    (void) q;
    (void) w;
    (void) e;
    (void) r;
    return 1l;
}

long negf_construct_nu_(double *q, double *w, double *e, double *r, long *t, double *y, double *u, double *i, double *o, long *p)
{
    (void) q;
    (void) w;
    (void) e;
    (void) r;
    (void) t;
    (void) y;
    (void) u;
    (void) i;
    (void) o;
    (void) p;
    return 1l;
}

long negf_calc_r_(double *q, double *w, double *e, long *r, double *t, double *y, double *u, long *i)
{
    (void) q;
    (void) w;
    (void) e;
    (void) r;
    (void) t;
    (void) y;
    (void) u;
    (void) i;
    return 1l;
}

long negf_create_hi_(double *q, double *w, double *e, double *r, long *t, double *y, double *u, double *i, double *o, long *p, double *a)
{
    (void) q;
    (void) w;
    (void) e;
    (void) r;
    (void) t;
    (void) y;
    (void) u;
    (void) i;
    (void) o;
    (void) p;
    (void) a;
    return 1l;
}

long negf_create_gam_(double *q, double *w, double *e, double *r, double *t, double *y, double *u, double *i, double *o, double *p, long *a, long *s, double *d)
{
    (void) q;
    (void) w;
    (void) e;
    (void) r;
    (void) t;
    (void) y;
    (void) u;
    (void) i;
    (void) o;
    (void) p;
    (void) a;
    (void) s;
    (void) d;
    return 1l;
}

long negf_set_initial_values_(double_complex *q, double *w)
{
    (void) q;
    (void) w;
    return 1l;
}

long negf_allocate_fortran_(long *q, long *w, long *e, long *r)
{
    (void) q;
    (void) w;
    (void) e;
    (void) r;
    return 1l;
}

long negf_rk_init_()
{
    return 1l;
}

void negf_create_h_(double *q)
{
    (void) q;
}

void negf_do_step_(double_complex *q, double_complex *w, double *e)
{
    (void) q;
    (void) w;
    (void) e;
}

void negf_calculate_current_(double_complex *q, double *w, double *e)
{
    (void) q;
    (void) w;
    (void) e;
}

