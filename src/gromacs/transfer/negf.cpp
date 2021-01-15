#include "gromacs/transfer/transfer.h"

int negf_init_arrays(ct_negf_lorentz_t *negf, double *rk_timestep, double *wf)
{
    extern long negf_const_mat_(double *, double *, long *, long *);
    extern long negf_jacobi_(double *, double *, double *, long *);
    extern long negf_construct_nu_(double *, double *, double *, double *, long *, double *, double *, double *, double *, long *);
    extern long negf_calc_r_(double *, double *, double *, long *, double *, double *, double *, long *);
    extern long negf_create_hi_(double *, double *, double *, double *, long *,
                                double *, double *, double *, double *, long *, double *);
    //extern long negf_create_gam_(double_complex *, double_complex *, double_complex *, double_complex *,
    //                             double_complex *, double_complex *, double_complex *, double_complex *,
    //	                       double *, double *, double *, double *, double *, double_complex *, double_complex *,
    //			       double *, double *, double *, double *, double *, double_complex *, double_complex *,
    //			       long *, long *, long *, long *, long *, double *);
    extern long negf_create_gam_(double *, double *, double *, double *, double *,
                                 double *, double *, double *, double *, double *,
                                 long *, long *, double *);
    extern long negf_set_initial_values_(double_complex *, double *);
    extern long negf_allocate_fortran_(long *, long *, long *, long *);
    extern long negf_rk_init_();

    long info;

    (void) rk_timestep;

    /* intialize */
    //ct->negf_arrays->kl[0] = ct->negf_arrays->n_poles_l[0] + ct->negf_arrays->n_lorentz[0];
    //ct->negf_arrays->kr[0] = ct->negf_arrays->n_poles_r[0] + ct->negf_arrays->n_lorentz[0];
    //snew(ct->negf_arrays->hi_left_p, ct->negf_arrays->kl[0]);
    //snew(ct->negf_arrays->hi_left_m, ct->negf_arrays->kl[0]);
    //snew(ct->negf_arrays->hi_right_p, ct->negf_arrays->kr[0]);
    //snew(ct->negf_arrays->hi_right_m, ct->negf_arrays->kr[0]);

    //snew(ct->negf_arrays->gam_greater_l_m, ct->negf_arrays->n[0] * ct->negf_arrays->n[0] *  ct->negf_arrays->kl[0]);
    //snew(ct->negf_arrays->gam_greater_l_p, ct->negf_arrays->n[0] * ct->negf_arrays->n[0] *  ct->negf_arrays->kl[0]);
    //snew(ct->negf_arrays->gam_lesser_l_m,  ct->negf_arrays->n[0] * ct->negf_arrays->n[0] *  ct->negf_arrays->kl[0]);
    //snew(ct->negf_arrays->gam_lesser_l_p,  ct->negf_arrays->n[0] * ct->negf_arrays->n[0] *  ct->negf_arrays->kl[0]);
    //snew(ct->negf_arrays->gam_greater_r_m, ct->negf_arrays->n[0] * ct->negf_arrays->n[0] *  ct->negf_arrays->kr[0]);
    //snew(ct->negf_arrays->gam_greater_r_p, ct->negf_arrays->n[0] * ct->negf_arrays->n[0] *  ct->negf_arrays->kr[0]);
    //snew(ct->negf_arrays->gam_lesser_r_m,  ct->negf_arrays->n[0] * ct->negf_arrays->n[0] *  ct->negf_arrays->kr[0]);
    //snew(ct->negf_arrays->gam_lesser_r_p,  ct->negf_arrays->n[0] * ct->negf_arrays->n[0] *  ct->negf_arrays->kr[0]);

    //snew(ct->negf_arrays->pi_l, ct->negf_arrays->n[0] * ct->negf_arrays->n[0] *  (ct->negf_arrays->kl[0] + 1));
    //snew(ct->negf_arrays->pi_r, ct->negf_arrays->n[0] * ct->negf_arrays->n[0] *  (ct->negf_arrays->kr[0] + 1));

    //snew(ct->negf_arrays->omega_ll1, ct->negf_arrays->n[0] * ct->negf_arrays->n[0] * ct->negf_arrays->n_lorentz[0] * (ct->negf_arrays->n_lorentz[0] + 1));
    //snew(ct->negf_arrays->omega_ll2, ct->negf_arrays->n[0] * ct->negf_arrays->n[0] * ct->negf_arrays->n_lorentz[0] * (ct->negf_arrays->kl[0] + 1 - ct->negf_arrays->n_lorentz[0]));
    //snew(ct->negf_arrays->omega_ll3, ct->negf_arrays->n[0] * ct->negf_arrays->n[0] * (ct->negf_arrays->kl[0] - ct->negf_arrays->n_lorentz[0]) * (ct->negf_arrays->n_lorentz[0] + 1));
    //snew(ct->negf_arrays->omega_lr1, ct->negf_arrays->n[0] * ct->negf_arrays->n[0] * ct->negf_arrays->n_lorentz[0] * (ct->negf_arrays->n_lorentz[0] + 1));
    //snew(ct->negf_arrays->omega_lr2, ct->negf_arrays->n[0] * ct->negf_arrays->n[0] * ct->negf_arrays->n_lorentz[0] * (ct->negf_arrays->kl[0] + 1 - ct->negf_arrays->n_lorentz[0]));
    //snew(ct->negf_arrays->omega_lr3, ct->negf_arrays->n[0] * ct->negf_arrays->n[0] * (ct->negf_arrays->kl[0] - ct->negf_arrays->n_lorentz[0]) * (ct->negf_arrays->n_lorentz[0] + 1));
    //snew(ct->negf_arrays->omega_rl1, ct->negf_arrays->n[0] * ct->negf_arrays->n[0] * ct->negf_arrays->n_lorentz[0] * (ct->negf_arrays->n_lorentz[0] + 1));
    //snew(ct->negf_arrays->omega_rl2, ct->negf_arrays->n[0] * ct->negf_arrays->n[0] * ct->negf_arrays->n_lorentz[0] * (ct->negf_arrays->kr[0] + 1 - ct->negf_arrays->n_lorentz[0]));
    //snew(ct->negf_arrays->omega_rl3, ct->negf_arrays->n[0] * ct->negf_arrays->n[0] * (ct->negf_arrays->kr[0] - ct->negf_arrays->n_lorentz[0]) * (ct->negf_arrays->n_lorentz[0] + 1));
    //snew(ct->negf_arrays->omega_rr1, ct->negf_arrays->n[0] * ct->negf_arrays->n[0] * ct->negf_arrays->n_lorentz[0] * (ct->negf_arrays->n_lorentz[0] + 1));
    //snew(ct->negf_arrays->omega_rr2, ct->negf_arrays->n[0] * ct->negf_arrays->n[0] * ct->negf_arrays->n_lorentz[0] * (ct->negf_arrays->kr[0] + 1 - ct->negf_arrays->n_lorentz[0]));
    //snew(ct->negf_arrays->omega_rr3, ct->negf_arrays->n[0] * ct->negf_arrays->n[0] * (ct->negf_arrays->kr[0] - ct->negf_arrays->n_lorentz[0]) * (ct->negf_arrays->n_lorentz[0] + 1));

    negf->length_rkvec[0] = negf->n[0] * negf->n[0] * (1 + negf->kl[0] + negf->kr[0] +
                                                       4 * negf->n_lorentz[0] * (negf->n_lorentz[0] + negf->n_poles_l[0] + negf->n_poles_r[0]));
    snew(negf->rkvec,  negf->length_rkvec[0]);
    snew(negf->drkvec, negf->length_rkvec[0]);

    negf->beta[0] = 1. / (BOLTZMANN_HARTREE_KELVIN * HARTREE_TO_EV * negf->temp[0]); /* make everything in eV within NEGF !!! */

    snew(negf->eig_val_mat_l,  2 * negf->n_poles_l[0]);
    snew(negf->eig_vect_mat_l, SQR(2 * negf->n_poles_l[0]));
    snew(negf->mat_l,          SQR(2 * negf->n_poles_l[0]));
    snew(negf->nu_l,           negf->n_poles_l[0]);
    snew(negf->r_alpha_l,      2 * negf->n_poles_l[0]);
    snew(negf->r_l,            negf->n_poles_l[0]);
    snew(negf->eig_val_mat_r,  2 * negf->n_poles_r[0]);
    snew(negf->eig_vect_mat_r, SQR(2 * negf->n_poles_r[0]));
    snew(negf->mat_r,          SQR(2 * negf->n_poles_r[0]));
    snew(negf->nu_r,           negf->n_poles_r[0]);
    snew(negf->r_alpha_r,      2 * negf->n_poles_r[0]);
    snew(negf->r_r,            negf->n_poles_r[0]);

    info = negf_allocate_fortran_(negf->n, negf->n_poles_l, negf->n_poles_r, negf->length_rkvec);

    info = negf_const_mat_(negf->mat_l, negf->mat_r, negf->n_poles_l, negf->n_poles_r);
    info = negf_jacobi_(negf->mat_l, negf->eig_val_mat_l, negf->eig_vect_mat_l, negf->n_poles_l);
    info = negf_jacobi_(negf->mat_r, negf->eig_val_mat_r, negf->eig_vect_mat_r, negf->n_poles_r);
    info = negf_calc_r_(negf->r_alpha_l, negf->eig_val_mat_l, negf->eig_vect_mat_l, negf->n_poles_l,
                        negf->r_alpha_r, negf->eig_val_mat_r, negf->eig_vect_mat_r, negf->n_poles_r);
    info = negf_construct_nu_(negf->eig_val_mat_l, negf->nu_l, negf->r_l, negf->r_alpha_l, negf->n_poles_l,
                              negf->eig_val_mat_r, negf->nu_r, negf->r_r, negf->r_alpha_r, negf->n_poles_r);

    //info = negf_create_hi_(negf->hi_left_m,  negf->hi_left_p,  negf->eps_l, negf->w0_l, negf->e_f_left,  negf->nu_l, negf->kl, negf->n_poles_l,
    //                       negf->hi_right_m, negf->hi_right_p, negf->eps_r, negf->w0_r, negf->e_f_right, negf->nu_r, negf->kr, negf->n_poles_r, negf->beta);
    //info = negf_create_gam_(negf->gam_greater_l_m, negf->gam_greater_l_p, negf->gam_lesser_l_m, negf->gam_lesser_l_p,
    //                        negf->gam_greater_r_m, negf->gam_greater_r_p, negf->gam_lesser_r_m, negf->gam_lesser_r_p,
//		          negf->gam_l, negf->w0_l, negf->eps_l, negf->e_f_left,  negf->r_l, negf->hi_left_m,  negf->hi_left_p,
//		          negf->gam_r, negf->w0_r, negf->eps_r, negf->e_f_right, negf->r_r, negf->hi_right_m, negf->hi_right_p,
//		          negf->n, negf->kl, negf->kr, negf->n_poles_l, negf->n_poles_r, negf->beta);
//info = negf_set_initial_values_(negf->rho, negf->pi_l, negf->pi_r,
//                                negf->omega_ll1, negf->omega_ll2, negf->omega_ll3, negf->omega_lr1, negf->omega_lr2, negf->omega_lr3,
//                                negf->omega_rl1, negf->omega_rl2, negf->omega_rl3, negf->omega_rr1, negf->omega_rr2, negf->omega_rr3,
//				  negf->rkvec,
//				  negf->n, negf->n_lorentz, negf->kl, negf->kr, negf->length_rkvec);
    info = negf_create_hi_(negf->eps_l, negf->w0_l, negf->e_f_left,  negf->nu_l, negf->n_poles_l,
                           negf->eps_r, negf->w0_r, negf->e_f_right, negf->nu_r, negf->n_poles_r, negf->beta);
    info = negf_create_gam_(negf->gam_l, negf->w0_l, negf->eps_l, negf->e_f_left,  negf->r_l,
                            negf->gam_r, negf->w0_r, negf->eps_r, negf->e_f_right, negf->r_r,
                            negf->n_poles_l, negf->n_poles_r, negf->beta);
    info = negf_set_initial_values_(negf->rkvec, wf);
    info = negf_rk_init_();

    (void) info;
    return 0;
}

int negf_propagate(charge_transfer_t *ct)
{
    extern void negf_create_h_(double *);
    extern void negf_do_step_(double_complex *, double_complex *, double *);
    extern void negf_calculate_current_(double_complex *, double *, double *);

    negf_create_h_(ct->hamiltonian[0]);
    negf_do_step_(ct->negf_arrays->rkvec, ct->negf_arrays->drkvec, &(ct->rk_timestep));
    negf_calculate_current_(ct->negf_arrays->rkvec, ct->negf_arrays->current, ct->wf);

    return 0;
}

