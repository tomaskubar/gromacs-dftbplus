#include "rksuite.h"


#define TFS_DIV_INTERVAL (15000)
//// Conversion Factors ////
//#define NM_TO_BOHR (18.897259886)
//#define HARTREE_TO_EV (27.211396132)
//#define AU_OF_ESP_TO_VOLT (27.211396132) // was 14.400 when lengths were in angstrom
//#define PS_TO_AU (41341.373) // picosecond to a.u. of time
//#define KJMOL_TO_HARTREE (0.000380880)
//#define AMU_TO_AU (1822.888486)
/* to convert AU to MD units: */
//#define HARTREE2KJ        4.3597482e-21
//#define BOHR2NM           0.0529177249
//#define HARTREE_BOHR2MD   (HARTREE2KJ*AVOGADRO/BOHR2NM)
//#define BOLTZMANN_HARTREE_KELVIN (8.314 / 1000 * KJMOL_TO_HARTREE)

//#define NEG_IMAG_POT (1.e-3)

/* #include "charge_transfer.h"
 *  - cannot be included because of a clash of definitions of bool type
 *    (i)  as standard type in C++
 *    (ii) as typedef int bool in gromacs                                */
//typedef int twointegers[2];
//typedef double twodoubles[2];
//typedef char sixstring[6];
//typedef double dvec[3];
//typedef struct {
//  int type;            /* specific type of this site */
//  int resnr;           /* residue number */
//  int atoms;           /* number of atoms for each residue */
//  int *atom;           /* lists of atoms (atomnumbers) */
//  int *atomtype;       /* atomtypes corresponding to atomnumbers */
//    	       /* C=0; H=1; N=2; O=3; */

//  int bonds;           /* number of bonds cut by QM/MM boundary. the following variables have one array for each bond. */ 
//  int *QMLA;           /* position of QM link atom in list "atom" */
//  int *MMLA;           /* position of MM link atom in list "atom" */
//  int *nochrs;         /* number of additional excluded charges (besides MMLA)*/
//  sixstring **nochr;   /* list of names of excluded atoms */
//  double *extracharge; /* remaining charge for electro-neutrality */
//  int *addchrs;        /* number of atoms over which the remaining charge will be distributed */
//  sixstring **addchr;  /* their names */
//  int **modif_extcharge;/* their atomnumbers */
//  int connections;     /* connections to neighboring fragments */
//  int *QMCA;           /* connection atom belonging to this site */
//  int **QMCN;          /* connection atom used as capping, belonging to a neigboring site */

//  int homos;           /* number of relevant (HO/LU)MOs for this site */
//  int *homo;           /* list of orbital numbers of relevant MOs */
////  double *occupation;  /* occupation of individual MOs of this sites by the excess charge */  
//  double *hubbard;     /* hubbard parameter. different values for holes in different MOs */
//  double *lambda_i;    /* inner-sphere reorganization energy. different values for holes in different MOs */
//  double **delta_q;    /* difference of atomic charges between neutral and cationic site; arrays for holes in different MOs */

//  int extcharges;      /* number of extcharges for this site */
//  int *extcharge;      /* lists of extcharges */
//  double **overlap;    /* overlap of site MOs between two steps */
//  double **overlap_ref;/* overlap of site MOs between actual step and the beginning of the simulation */
//  int nel;             /* total number of electrons */
//  int norb;            /* total number of orbitals */
//  double radius;       /* radius of charged sphere in polarizable continuum */
//  int do_scc;          /* each site may or may not be calculated self consistently */
//  int do_custom_occ;   /* each site may also use a customized occupation, instead of fermi distribution */
//  double *custom_occ;      
//  double *com;         /* center of mass. needed to decide if site should become active or not */
//  int active;          /* switch 0/1 that determines if site is part of the QM calculation or just inactive member of the pool */
//} ct_site_t;

//typedef struct {
//  double *in,
//         *evec,
//         *work,
//         *iwork,
//         *eval,
//         *issupz;
//  long lwork,
//       liwork;
//} ct_per_orthogo_t;

//typedef struct {double dr, di;} double_complex;

//typedef struct {
//  long n[1], n_lorentz[1];
//  double e_f_left[1], e_f_right[1],
//         temp[1];
//  long n_poles_l[1], n_poles_r[1];
//  double t_0[1], t_ned[1], t_step[1],
//         gam_l[1], eps_l[1], w0_l[1],
//         gam_r[1], eps_r[1], w0_r[1],
//         beta[1];
//  long kl[1], kr[1];
//  double_complex *pi_l, *pi_r,
//                 *omega_ll1, *omega_ll2, *omega_ll3, *omega_rr1, *omega_rr2, *omega_rr3,
//                 *omega_lr1, *omega_lr2, *omega_lr3, *omega_rl1, *omega_rl2, *omega_rl3,
//                 *gam_greater_l_m, *gam_greater_l_p, *gam_greater_r_m, *gam_greater_r_p,
//                 *gam_lesser_l_m,  *gam_lesser_l_p,  *gam_lesser_r_m,  *gam_lesser_r_p,
//                 *hi_left_p, *hi_left_m, *hi_right_p, *hi_right_m,
//                 *h,
//                 *rho;
//  long length_rkvec[1];
//  double_complex *rkvec, *drkvec;
//  double *mat_l, *mat_r, *eig_vect_mat_l, *eig_vect_mat_r,
//         *eig_val_mat_l, *eig_val_mat_r, *nu_l, *nu_r, *r_alpha_l, *r_alpha_r, *r_l, *r_r;
//  double current[3];
//} ct_negf_lorentz_t;

//typedef struct {
//  double **U,   //unitary transformation matrix = alignment
//    *sij,
//    *evec,
//    *work,
//    *iwork,
//    *eval,
//    *issupz;
//  long lwork,
//    liwork;
//} align_t;

//typedef struct {
//  /* main data structure of the charge transfer code */
//  int jobtype;         /* cteXXX */
//  int interval;        /* how often should the parameters be calculated? */
//  int qmmm;            /* 1 or 0 */
//  int delta_q_mode;    /* 0 use Mulliken charges, 1 use RESP charges to describe the charge carrier */
//  int sitetypes;       /* number different types of sites e.g. adenine and guanine */
//  ct_site_t *sitetype; /* list of different types of sites */
//  int sites;           /* number of sites */
//  ct_site_t *site;     /* list of the considered nucleobases */   // 		               NOW STRUCT BY ITS OWN
//  int *do_scc;         /* each site may or may not be calculated self consistently */
//  int do_scc_CG;       /* CG Hamiltonian may or may not be calculated self consistently */
//  int dim;             /* dimension of hamilton matrix and so on. (not the same as the number of sites since multiple MOs per site are possible) */
//  int is_hole_transfer; /* electron or hole transfer */
//  double offdiag_scaling; /* 1.540 scaling factor according to J. Chem. Phys. 140, 104105 (2014) */
//  int is_protein;      /* proteins need new residues built by protein_preprocessor. Not implemented yet*/
////  int *homo;           /* which orbital is HOMO? */				//	NOW IN STRUCT SITE
////  int *atoms;          /* number of atoms for each residue */				NOW IN STRUCT SITE
//  int *last_atom;      /* index of the last atom of each residue, in the array of atoms of the complex */
////  int **atom;          /* lists of atoms - one array for each nucleobase */			NOW IN STRUCT SITE
//                       /* first element is C1q while reading coordinates, */
//                       /*               and is then converted to link hydrogen */
////  int **atomtype;      /* similarly as in the previous case */				NOW IN STRUCT SITE
//                       /* C=0; H=1; N=2; O=3; */
//  int atoms_cplx;      /* number of atoms for the complex */
//  int *atom_cplx;      /* list of atoms for the complex */
//                       /* contains C1q in place of link hydrogens */
//  int *atomtype_cplx;  /* similarly as in the previous case */
//                       /* C=0; H=1; N=2; O=3; */
//  int extcharges_cplx; /* number of extcharges for the complex */        
////  int **extcharge;     /* lists of extcharges - one array for each site */                         NOW IN STRUCT SITE
//  int *extcharge_cplx; /* list of extcharges - for the complex */
////  twointegers *modif_extcharge;     /* which extcharges correspond to those of O4a and C2q? will be modified, +0.06080 */        NOW IN STRUCT SITE
//  int modif_extcharges_cplx;        /* atomnumber of charges which will be modified */
//  int *modif_extcharge_cplx;        /* their atomnumber */
////  double **delta_q;    /* difference of charge of neutral nucleobase and cation; array for each nucleobase */			NOW IN STRUCT SITE
//  double *fo_shift;    /* shift that will be applied to the diagonal elements of the FO hamiltonian. can correct the difference between HOMO differences and IP differences */
//  double **hamiltonian;/* CG Hamiltonian, calculated by DFTB, to be used in TDSE integration */
//                       /* n x n matrix, fortran format */
//  double ***hamiltonian_history;/* history over last n_avg_ham steps of CG Hamiltonian, to average fast (nonclassical) oscillations */
//  int n_avg_ham;       /*length of hamiltonian_history */
//  double *ev_spec; 
//  double *evec_spec;   /* used for hamilton-diagonalization to get CG-MO spectrum */
//  double *work_spec;

//  double *hamiltonian_mod; /* the contribution of Hubbard and electrostatics to the site energies - diagonal of hamiltonian[][] */
//  double *hamiltonian_adiab; /* to be used in the adiabatic dynamics */
//  double *tfl_mean_ham;  // new from jjk
//  double **tfl_mean_ham_full; //new jjk
//  double *tfl_old_ham;  //TFL new jjk
//  double *tfl_old_ham_adiab; //new jjk
//  int* indFO;          /* index of site i in Hamilton in FO-basis */
//  double *ev_adiab;    /* to be used in the adiabatic dynamics */
//  double *evec_adiab;  /* to be used in the adiabatic dynamics */
//  double *work_adiab;  /* to be used in the adiabatic dynamics */
//  double survival;     /* magnitude of the remaining excess charge - sum of the occupations */
//  double *occupation;  /* occupation of individual MOs by the excess charge. essentially |wf|^2 ?????? */          
//  double *pop_adia;  /* adiabatic populations*/ 
//  double **hubbard;    /* diag:    Hubbard parameter - chemical hardness of the nucleobase
//                          offdiag: 1/r_ij - Coulomb interaction of parts of hole on two nucleobases */
//  double sic;          /* naive self-interaction correction - scaling factor for the 2nd order term; 1.0 if no self-interaction correction */
//  double esp_scaling_factor; /* divide the ESP (i.e. the DFTB shift) by this value - emulation of electronic polarization; 1.0 of no ESP scaling */
//  int do_lambda_i;     /* shall the inner-sphere reorganization energy be emulated? 1==YES, 0==NO */
//  int do_epol;         /* shall the electronic polarization be calculated? 1==YES, 0==NO */
//  int do_projection;   /* do or do not project the wavefunction at every step onto the new FO basis before propagating */

//  int neg_imag_pot;    /* shall the negative imaginary potential be applied? 0==NO, if > 0 then it means the number of sites to be applied on */
//  int *site_neg_imag_pot; /* list of sites (numbered 1, 2, 3...) to apply the negative imaginary potential */
//  double neg_imag_pot_rate_constant; /* in 1/a.u. of time */
//  double *site_annihilated_occupation;
//  double survival_threshold;
//  int rk_neq;          /* Runge-Kutta: number of equations = 2 * number of sites */
//  double rk_tol;       /* tolerance for Runge-Kutta */
//  int adiabstart;      /* start SCC from lowest adiabatic state */
//  double *wf;          /* CG wavefunction, dimension 2n: real[0] real[1] real[2] ... real[n-1] imag[0] imag[1] imag[2] ... imag[n-1] */
//  double *wf_old;      /* CG wavefunction, before projection on new basis */
//  double *wf_exc;
//  double *q_act;
//  double *q_old;
//  double *dwf;         /* derivative of CG wavefunction */
//  double *rk_ymax;     /* Runge-Kutta: auxiliary array */
//  double *rk_thres;    /* Runge-Kutta: auxiliary array */
//  double *rk_work;     /* Runge-Kutta: auxiliary array */
//  int rk_lenwrk;       /* Runge-Kutta: length of rk_work = 32 * rk_neq */
//  double rk_timestep;  /* time step for TDSE propagation in atomic units of time */
//  int surface;         /* occupied surface for surface hopping */
//  double *surf_overlap, *surf_massey, *surf_prob;
//  double telec;       /* electronic temperature (in Kelvin) for the Fermi distribution */
//  double fermi_kt;     /* kT(elec), for the Fermi distribution */
//  /* variables for Tully's fewest switches - surface hopping */
//  double *tfs_popul;       /* vector of population of adiabatic states */
//  double *tfs_popul_der;   /* time derivative of tfs_popul [1/au_of_time] */
//  double *tfs_diab;       //populations as diabatic wavefct. new  jjk
//  double **tfs_vector;     /* adiabatic states phi_k(t+dt) */
//  double **tfs_vector_old; /* adiab. states in the previous step phi_k(t) */
//  double **tfs_overlap;    /* <phi_k(t) | phi_j(t+dt)> */
//  double ***tfs_ncv;  /* nonadiabatic coupling vector*/
//  double ***v1;  /* horizontal momenta */
//  double ***v2;  /* vertical momenta */
//  int tfs_initialization_step; /* is this the first step of a surface hopping simulation? */
//  //variables for flexible surface hopping
//  int tfl_num_of_states; //current number of sites in sys.                                                                          
//  int tfl_num_of_states_old;  //# of sites at last step                                                                             
//  int *tfl_is_in_system; //...[i]==1 if site i in sys., otherwise ...[i]==0                                                         
//  int *tfl_is_in_system_old;  //at last step                                                                                        
//  double tfl_ca_real;   //coef. of surface before integration (real part)                                                          
//  double tfl_ca_im;    //(imaginary part) 
///* variables for Landau-Zener surface hopping */
///* variables for Persico's local diabatic surface hopping */
//  twodoubles **per_propag_operator; /* exp[-i*Z*dt] */
//  twodoubles **per_transformator;   /* U = T^t * exp[-i*Z*dt] */
//  ct_per_orthogo_t *per_arrays;
//  double *tfs_popul_old;
//  double *ev_adiab_old;
//  double **per_diab_hamiltonian;
//  int decoherence;     /* shall the decoherence correction be applied? 0==NO, 1==YES */
//  ct_negf_lorentz_t *negf_arrays;
//  int *sort_initialization_step; /*(new jjk) is this the first time one has to sort orbitals on site i? */
//  align_t align; //new jjk
//  dvec efield;  // external electric field
//  int first_step; // is first time QM calculations? may differ from "step==0" from gromacs-steps in case of reruns. general purpose variable can be used by different routines (in contrast to tfs_initialization_step)
//  double *born_overlap; //overlap used in cteBORNOPPENHEIMER
//  int *incl_nonad;   //==1 if non-adiab coupling to be calculated, else ==0 //EX new jjk
//  int use_bc_correction; //=1 if Boltzmann correction to be used in Ehrenfest//BC new jjk
//  int use_strong_dec; //=1 if TFL&collapse after hopp new jjk

//  // variables for adaptive QM zone
//  int pool_size;          // number of possible sites (active sites are ct->sites)
//  ct_site_t *pool_site;  // possible site (active ones are ct->site[i])
//  dvec coc;              // center of charge[bohr] (determines which sites become active)
//  dvec coc_start;        // center of charge[bohr] at the beginning of the simulation
//  int opt_QMzone;        // switch 0/1 that derermines if optimal QM zone should be generated from pool of sites 
//  double adapt_inv_tot_mass; // inverse of the total mass of one site

//  double **tfs_overlap2; /* eigenvectors at two time steps before */
//  double **tfs_vector_old2; /* eigenvectors at two time steps before */
//  double *ev_adiab_old2;  // adiabatic energy at two time steps before
//  int ncv_isurf; // index for nonadiabtic coupling vecotr d_{ij}
//  int ncv_jsurf; // index for nonadiabtic coupling vecotr d_{ij}
//  int *ind_alzsh;
// // seed for random number  
//  int rnseed;

//  double **fo_overlap;
//  double **fo_tdc;

//// implicite lambda
//  double lambda;

//// machine learning parameters
//  long   ntrain[2];
//  double sigma[2];
//  double hamiltonian_type;
//  double ***coords;
//  double ***cm_train;
//  double **alpha;
//  double **hamiltonian_ml;
//  double **hamiltonian_dftb;

//// modified eherenfest methods
//  int scctype;  // 0: SCC, 1: coherence penatly functional, 2: decoherence detailed balance
//  double dephasing_rate;

//} charge_transfer_t;

/* calculate the derivative of var and save it in array dvar */
/* read in also the information on the neg imag potential! */
void eqom(double time, double *var, double *dvar, int n, double **hamiltonian, double **hubbard,
          int nip, int *site_nip, double nip_rate_constant)
{
  int i, j;
  double add_to_hamiltonian_ii, occupation_j;
  (void) time;
  
  for (i=0; i<n; i++) {
    add_to_hamiltonian_ii = 0.e0;
    for (j=0; j<n; j++) {
      occupation_j = SQR(var[j]) + SQR(var[n+j]);
      add_to_hamiltonian_ii += hubbard[i][j] * occupation_j;
    }
    dvar[i] = 0.e0;
    dvar[n+i] = 0.e0;
    for (j=0; j<n; j++) {
      /* real component of derivative */
      dvar[i] += (hamiltonian[i][j] + (i==j ? add_to_hamiltonian_ii : 0.e0)) * var[n+j];
      /* imaginary component of derivative */
      dvar[n+i] -= (hamiltonian[i][j] + (i==j ? add_to_hamiltonian_ii : 0.e0)) * var[j];
    }
  }

  /* negative imaginary potential */
  for (i=0; i<nip; i++) {
    /* real component of the derivative */
    dvar[site_nip[i]]   -= nip_rate_constant * var[site_nip[i]];
    /* imaginary component of the derivative */
    dvar[site_nip[i]+n] -= nip_rate_constant * var[site_nip[i]+n];
  }

  return;
}

/* dtto with the negative imaginary potential
void eqom_nip(double time, double *var, double *dvar, int n, double **hamiltonian, double **hubbard,
              int nip, int *site_nip, double nip_rate_constant)
{
  int i, j;
  double add_to_hamiltonian_ii, occupation_j;
  
  for (i=0; i<n; i++) {
    add_to_hamiltonian_ii = 0.e0;
    for (j=0; j<n; j++) {
      occupation_j = SQR(var[j]) + SQR(var[n+j]);
      add_to_hamiltonian_ii += hubbard[i][j] * occupation_j;
    }
    dvar[i] = 0.e0;
    dvar[n+i] = 0.e0;
    for (j=0; j<n; j++) {
      // real component of derivative
      dvar[i] += (hamiltonian[i][j] + (i==j ? add_to_hamiltonian_ii : 0.e0)) * var[n+j];
      // imaginary component of derivative
      dvar[n+i] -= (hamiltonian[i][j] + (i==j ? add_to_hamiltonian_ii : 0.e0)) * var[j];
    }
  }

  // "negative imaginary potential"
  for (i=n-3; i<n; i++) {
    dvar[i] -= NEG_IMAG_POT * var[i];
    dvar[i+n] -= NEG_IMAG_POT * var[i+n];
  }
  // one last site only!
  // dvar[n-1] -= NEG_IMAG_POT * var[n-1];
  // dvar[2*n-1] -= NEG_IMAG_POT * var[2*n-1];

  return;
}
*/

// extern "C" int do_rksuite(charge_transfer_t *ct)
int do_rksuite(charge_transfer_t *ct)
{
  int method = 3;
  const char task[2] = "U";
  bool errass = false;
  bool mesage = false;
  int ifail;
  
double tgot, twant;

  twant = ct->rk_timestep;

  /* initial setup of Runge-Kutta */
  rk_setup(ct->rk_neq, 0.e0, ct->wf, MAX_INTEGRATION_STEP, ct->rk_tol, ct->rk_thres,
           method, task, errass, 0.e0, ct->rk_work, ct->rk_lenwrk, mesage);

  /* do a step */
//  if (ct->neg_imag_pot == 1)
//    rk_ut(eqom_nip, twant, tgot, ct->wf, ct->dwf, ct->rk_ymax, ct->rk_work, ifail, ct->sites, ct->hamiltonian, ct->hubbard);
//  else
    rk_ut(eqom, twant, tgot, ct->wf, ct->dwf, ct->rk_ymax, ct->rk_work, ifail, ct->dim, ct->hamiltonian, ct->hubbard,
          ct->neg_imag_pot, ct->site_neg_imag_pot, ct->neg_imag_pot_rate_constant);
  // on hard or catastrophic errors:
  if (ifail > 4) printf("RKsuite: error %d!\n", ifail);
  // on soft errors:
  while ( (ifail==2 || ifail==3 || ifail==4) && tgot<twant ) {
    ifail = -1;
//    if (ct->neg_imag_pot == 1)
//      rk_ut(eqom_nip, twant, tgot, ct->wf, ct->dwf, ct->rk_ymax, ct->rk_work, ifail, ct->sites, ct->hamiltonian, ct->hubbard);
//    else
      rk_ut(eqom, twant, tgot, ct->wf, ct->dwf, ct->rk_ymax, ct->rk_work, ifail, ct->sites, ct->hamiltonian, ct->hubbard,
            ct->neg_imag_pot, ct->site_neg_imag_pot, ct->neg_imag_pot_rate_constant);
    if (ifail > 4) printf("RKsuite: error %d!\n", ifail);
  }

  return 0;
}

//////////rk_suite for Tully_diab



// extern "C" int do_rksuite_diab(charge_transfer_t *ct)
int do_rksuite_diab(charge_transfer_t *ct)
{ 
  int method = 3;
  const char task[2] = "U";
  bool errass = false;
  bool mesage = false;
  int ifail;
  
double tgot, twant;
 
// printf("check_0"); fflush(stdout);
 
  twant = ct->rk_timestep;
  
//printf("%d,%d,%d,%d,%d,%d,%d,",rk_neq,tfs_diab,MAX_INTEGRATION_STEP,rk_tol,rk_thres,rk_work,rk_lenwrk);
/* initial setup of Runge-Kutta */
  rk_setup(ct->rk_neq, 0.e0, ct->tfs_diab, MAX_INTEGRATION_STEP, ct->rk_tol, ct->rk_thres,
           method, task, errass, 0.e0, ct->rk_work, ct->rk_lenwrk, mesage);


//  printf("check_1"); fflush(stdout);

  /* do a step */
//  if (ct->neg_imag_pot == 1)
//    rk_ut(eqom_nip, twant, tgot, ct->wf, ct->dwf, ct->rk_ymax, ct->rk_work, ifail, ct->sites, ct->hamiltonian, ct->hubbard);
//  else
    rk_ut(eqom, twant, tgot, ct->tfs_diab, ct->dwf, ct->rk_ymax, ct->rk_work, ifail, ct->dim, ct->tfl_mean_ham_full, ct->hubbard,
          ct->neg_imag_pot, ct->site_neg_imag_pot, ct->neg_imag_pot_rate_constant);
//    printf("check_2"); fflush(stdout);


  // on hard or catastrophic errors:
  if (ifail > 4) printf("RKsuite: error %d!\n", ifail);
  // on soft errors:
  while ( (ifail==2 || ifail==3 || ifail==4) && tgot<twant ) {
    ifail = -1;
//    if (ct->neg_imag_pot == 1)
//      rk_ut(eqom_nip, twant, tgot, ct->wf, ct->dwf, ct->rk_ymax, ct->rk_work, ifail, ct->sites, ct->hamiltonian, ct->hubbard);
//    else
      rk_ut(eqom, twant, tgot, ct->tfs_diab, ct->dwf, ct->rk_ymax, ct->rk_work, ifail, ct->sites, ct->tfl_mean_ham_full, ct->hubbard,
            ct->neg_imag_pot, ct->site_neg_imag_pot, ct->neg_imag_pot_rate_constant);
    if (ifail > 4) printf("RKsuite: error %d!\n", ifail);
  }

  return 0;
}




/////////////end rk_suite for Tully_diab new jjk may18





/* calculate the derivative of var
 * to be used in do_rksuite_tfs (Tully's fewest switches)
 * and save it in array dvar
 */
void eqom_tfs(double time, double *var, double *dvar, int n, double *energy, double **overlap, double dt)
{
  int i, j;
  (void) time;

  for (i=0; i<n; i++) {
    dvar[i]   =   var[n+i] * energy[i];
    dvar[n+i] = - var[i] * energy[i];
    for (j=0; j<n; j++) if (i!=j) {
      dvar[i]   -= var[j]   * overlap[i][j] / dt;
      dvar[n+i] -= var[n+j] * overlap[i][j] / dt;
    }
  }
  
  return;
}

// extern "C" int do_rksuite_tfs(charge_transfer_t *ct)
int do_rksuite_tfs(charge_transfer_t *ct)
{
  int method = 3;
  const char task[2] = "U";
  bool errass = false;
  bool mesage = false;
  int ifail;
  
  double tgot;
  /* TWANT IS VARIABLE NOW!
   * twant =  (i+1) * ct->rk_timestep / TFS_DIV_INTERVAL
  twant = ct->rk_timestep;
  */

  int i, j;

  for (i=0; i<ct->sites; i++)
    ct->surf_prob[i] = 0.;
  ct->surf_prob[ct->surface] = -1.;

  /* initial setup of Runge-Kutta */
  rk_setup(ct->rk_neq, 0.e0, ct->tfs_popul, MAX_INTEGRATION_STEP, ct->rk_tol, ct->rk_thres,
           method, task, errass, 0.e0, ct->rk_work, ct->rk_lenwrk, mesage);

  /* do a sequence of multiple steps */
  for (j=0; j<TFS_DIV_INTERVAL; j++) {
    /* this is the step */
    //printf("starting iteration %d\n", j);
    //printf("twant = %f\n", (j+1) * ct->rk_timestep / TFS_DIV_INTERVAL);
    rk_ut_tfs(eqom_tfs, (j+1) * ct->rk_timestep / TFS_DIV_INTERVAL, tgot, ct->tfs_popul, ct->tfs_popul_der, ct->rk_ymax, ct->rk_work, ifail,
          ct->dim, ct->ev_adiab, ct->tfs_overlap, ct->rk_timestep);
    // on hard or catastrophic errors:
    if (ifail > 4) printf("RKsuite: error %d!\n", ifail);
    // on soft errors:
    while ( (ifail==2 || ifail==3 || ifail==4) && tgot < (j+1) * ct->rk_timestep / TFS_DIV_INTERVAL ) {
      ifail = -1;
    //printf("twant = %f\n", (j+1) * ct->rk_timestep / TFS_DIV_INTERVAL);
        rk_ut_tfs(eqom_tfs, (j+1) * ct->rk_timestep / TFS_DIV_INTERVAL, tgot, ct->tfs_popul, ct->tfs_popul_der, ct->rk_ymax, ct->rk_work, ifail,
              ct->sites, ct->ev_adiab, ct->tfs_overlap, ct->rk_timestep);
      if (ifail > 4) printf("RKsuite: error %d!\n", ifail);
    }
    /* end of the step */

    /* the following is for the surface-hopping probability */
    for (i=0; i<ct->sites; i++) if (i != ct->surface)
      ct->surf_prob[i] += ct->tfs_popul[i] * ct->tfs_popul[ct->surface] + ct->tfs_popul[ct->sites + i] * ct->tfs_popul[ct->sites + ct->surface];
  }

  /* normalize the probabilities */
  for (i=0; i<ct->sites; i++) if (i != ct->surface)
    ct->surf_prob[i] /= TFS_DIV_INTERVAL;

  return 0;
}

