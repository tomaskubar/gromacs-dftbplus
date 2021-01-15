#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <ctype.h>
#include <math.h>
#if GMX_MPI
#include <mpi.h>
#endif
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tensorflow/c/c_api.h>

#include "gromacs/domdec/collect.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_network.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/mdsetup.h"
#include "gromacs/domdec/partition.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/checkpointhandler.h"
#include "gromacs/mdlib/compute_io.h"
//#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdlib/resethandler.h"
#include "gromacs/mdlib/update.h"
#include "gromacs/mdlib/vcm.h"
#include "gromacs/mdlib/vsite.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/forcebuffers.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/timing/walltime_accounting.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/index.h"
//#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

#if GMX_MPI
#define PRINTF(...) if (ct_mpi_rank == 0) printf(__VA_ARGS__)
#else
#define PRINTF(...) printf(__VA_ARGS__)
#endif

#define SQR(x) ((x)*(x))
#define CUB(x) ((x)*(x)*(x))
#define QRT(x) ((x)*(x)*(x)*(x))

//////// DFTB/QM parameters ///////
#define DFTB_MAXTYPES (11)     /* number of chemical elements known to DFTB - C,H,N,O,S,F,Br,Cl,B,P,I */
                              /* Br,Cl */
#define MAX_PATH_LENGTH (288) /* max length of the string, which holds the path to the SLKO files */
#define LDIM (9)              /* size of profylactic arrays - 1 s, 3 p and 5 d orbital components = 9 */
#define MAXITER_BROYDEN (80)  /* max number of iterations for self consitent charge calculations with the broyden method */
#define MAXITER_DIIS (200)    /* max number of iterations for self consitent charge calculations with the diis method */
#define IMATSZ_BROYDEN (80)   /* size of inverse matrix that is used in the broyden method*/
#define MAXITER_SCC (70)      /* maybe a relict? SCC iterations are now given directly in the respective function */


///// GAFF parameter used in the construction of link hydrogen atoms /////
   // the length of the bond in Angstrom
#define CH_BOND_LENGTH (1.093)
#define NH_BOND_LENGTH (1.018)
#define OH_BOND_LENGTH (0.974)
#define SH_BOND_LENGTH (1.353)
#define BH_BOND_LENGTH (1.190)
#define PH_BOND_LENGTH (1.420)


//// Conversion Factors ////
#define NM_TO_BOHR (18.897259886)
#define HARTREE_TO_EV (27.211396132)
#define HARTREE_TO_CM (219474.63068)
#define AU_OF_ESP_TO_VOLT (27.211396132) // was 14.400 when lengths were in angstrom
#define PS_TO_AU (41341.373) // picosecond to a.u. of time
#define KJMOL_TO_HARTREE (0.000380880)
#define AMU_TO_AU (1822.888486)
/* to convert AU to MD units: */
//#define HARTREE2KJ        4.3597482e-21
//#define BOHR2NM           0.0529177249
//#define HARTREE_BOHR2MD   (HARTREE2KJ*AVOGADRO/BOHR2NM)
#define BOLTZMANN_HARTREE_KELVIN (RGAS / 1000 * KJMOL_TO_HARTREE)


//// MISC. ////
#define OFFDIAG_FACTOR_HOLE (1.540) // scaling of offdiagonal elements of FO Hamiltonian. See: J. Chem. Phys. 2014, 140, 104105+
#define OFFDIAG_FACTOR_ELEC (1.795) // scaling of offdiagonal elements of FO Hamiltonian. See: Phys. Chem. Chem. Phys. 2015, 17, 14342-14354.
#define EPSILON_OP (2.0)            // quite general value for high frequency permittivity. Needed for implicit solvent polarization.
#define N_OLDVEC (5)                // number of time steps a_old is kept
#define QMCASIZE (50000)            // max number of atoms named connection atoms in the whole system
#define PROTEIN (0)                 // enable protein version (activate protein_preprocessor)
#define GIESEPEPTIDE (0)            // if GIESEPEPTIDE==1, switch to hard-coded implementation of TDA for one specific peptide(Phys. Chem. B 118, 2014, 4261-4272)

#define ADIAB_SCF_CONVERGENCE (1.e-12) // used in the minimization of Hamiltonian within the adiabatic dynamics
#define DIIS_INIT_MIXING_PARAM (0.2)
#define DIIS_MAX_PREV_VECTORS (4)
#define BROYDEN_ALMIX (0.2)
#define BROYDEN_SCFTOL (1.e-9)
#define FERMI_KT (9.5004455e-4) // kT in hartree units, at 300 K
#define FERMI_CONVERG (1.e-12) // 1.e-9 kT in hartree units (at 300 K)
#define SIMPLE_ALMIX (0.01)
#define ALMIX_ATTENUATOR (0.9)
#define MAXLINES (50)      // maximal number of lines in input files
#define MAXWIDTH (10000)    // maximal width of key- and value-field in input files
#define MAXELEMENTS (300)  // maximal elements of in a value-field
#define ELEMENTWIDTH (100) // maximal length of a single element in the value-field
#define MAXSITETYPES (7)   // maximal number of different types of sites. e.g. guanine-adenine-guanine stack would have two different sitetypes

#define NEG_IMAG_POT (1.e-3)  // tau for the negative imaginary potential on the last couple of sites
/* #define SURVIVAL_THRESHOLD (0.1)  // if the survival drops under this value, finish */

#define SFH_ORBITAL_SHIFT (0.3674931) // energy shift of the ground state (down)
                             // 1 eV in hartree units
//#define SFH_LIMIT (0.4042424) // max. energy gap for which surface hopping shall be attempted
#define SFH_LIMIT (0.03674931) // max. energy gap for which surface hopping shall be attempted
                                // 1 eV in hartree units
#define SFH_EXPONENT (0.5) // the parameter alpha in Fabiano et al. Chem. Phys. 351 (2008) 111-116, Eq. 13
                           // 0.5 is their recommended value
#define TFS_DECAY_CONSTANT (3) // used in the decoherence algorithm //dholub changed 3.0 to 5..
                                // within Tully's fewest switches

#define TFL_RC (1.0/40.0)  // (1.0/80) new jjk critical fra5ction for flexible surface hopping approach

#define MAX_PME_NEIGHBORS (1666)

#define JFSSH_OFFDIAG_FORCE (1000)

#define SQR(x) ((x)*(x))




///// Old DNA/TRP/TYR stuff  /////
//#define HUBBARD_ADENINE (0.208)
//#define HUBBARD_GUANINE (0.205)
//#define EXTCHARGE_SHIFT (0.06080)
//#define LAMBDA_I (0.008452341) // inner-sphere reorganization energy of a nucleobase. 0.23 eV in hartree units
//#define SIC_COEFFICIENT (0.2) // self-interaction correction - only 20 % of the QQ interaction involved. Now you have to give the SIC coefficient in the input



/* some simple custom data types */
typedef int twointegers[2];
typedef double tendoubles[10];
typedef double twodoubles[2];
typedef double twodoubles_array[2][MAXITER_BROYDEN];
typedef char sixstring[6]; //gromacs atom names are 5 chars long (+string terminator \0)
typedef int pme_integers[MAX_PME_NEIGHBORS];
typedef struct {double dr, di;} double_complex;


#if (GMX_TENSORFLOW)
typedef struct {
    /* main data structure for NN ML
     * Mila Oct 2020
     */
  TF_Session* session;
  TF_Graph* graph;
  TF_Status* status;
  TF_Output* input;
  TF_Output* output;
  int num_inputs;
  int num_outputs;
  TF_Tensor** input_val;
  TF_Tensor** output_val;
  float x_mean, x_std;
  float y_mean, y_std;
  float grad_std;
} tf_model;
#endif

/* available jobtypes */
enum { cteSCCDYNAMIC, cteCPFSCCDYNAMIC, cteDDBSCCDYNAMIC, cteADIABATIC, cteBORNOPPENHEIMER, cteNONSCCDYNAMIC, ctePARAMETERS, cteADNONSCC, cteNOMOVEMENT, cteSURFACEHOPPING, cteFERMI, cteFERMIADIABATIC, cteFERMISFHOPPING, cteTULLYFEWESTSWITCHES, cteDLZSH,  cteALZSH, cteDFSSH, cteJFSSH, cteBCJFSSH, cteSCRDFSSH, cteCCFSSH, cteDISH, cteGFSH, cteTULLYLOC, cteTULLYDIA, ctePERSICOSFHOPPING, cteNEGFLORENTZ, cteNEGFLORENTZNONSCC, cteESP, cteTDA, ctePREZHDOSFHOPPING, cteNR };

typedef struct {
  /* data structure for non-equilibrium greens function calculation. I guess this was used for calculating transmissions in DNA between electrodes */
  long n[1], n_lorentz[1];
  double e_f_left[1], e_f_right[1],
         temp[1];
  long n_poles_l[1], n_poles_r[1];
  double t_0[1], t_ned[1], t_step[1],
         gam_l[1], eps_l[1], w0_l[1],
         gam_r[1], eps_r[1], w0_r[1],
	 beta[0];
  long kl[1], kr[1];
  double_complex *pi_l, *pi_r,
                 *omega_ll1, *omega_ll2, *omega_ll3, *omega_rr1, *omega_rr2, *omega_rr3,
		 *omega_lr1, *omega_lr2, *omega_lr3, *omega_rl1, *omega_rl2, *omega_rl3,
		 *gam_greater_l_m, *gam_greater_l_p, *gam_greater_r_m, *gam_greater_r_p,
		 *gam_lesser_l_m,  *gam_lesser_l_p,  *gam_lesser_r_m,  *gam_lesser_r_p,
		 *hi_left_p, *hi_left_m, *hi_right_p, *hi_right_m,
		 *h,
		 *rho;
  long length_rkvec[1];
  double_complex *rkvec, *drkvec;
  double *mat_l, *mat_r, *eig_vect_mat_l, *eig_vect_mat_r,
         *eig_val_mat_l, *eig_val_mat_r, *nu_l, *nu_r, *r_alpha_l, *r_alpha_r, *r_l, *r_r;
  double current[3];
} ct_negf_lorentz_t;

typedef struct {
  /* used in persico surface hopping */
  double *in,
         *evec,
         *work,
         *iwork,
         *eval,
         *issupz;
  long lwork,
       liwork;
} ct_per_orthogo_t;


typedef struct {
  /* data structure that stores information about one fragment */
  int type;            /* specific type of this site */
  int resnr;           /* residue number */
  int atoms;           /* number of atoms for each residue */
  int *atom;           /* lists of atoms (atomnumbers) */
  int *atomtype;       /* atomtypes corresponding to atomnumbers */
		       /* C=0; H=1; N=2; O=3; */
  int bonds;           /* number of bonds cut by QM/MM boundary. the following variables have one array for each bond. */
  int *QMLA;           /* position of QM link atom in list "atom" */
  int *MMLA;           /* position of MM link atom in list "atom" */
  int *nochrs;         /* number of additional excluded charges (besides MMLA)*/
  sixstring **nochr;   /* list of names of excluded atoms */
  double *extracharge; /* remaining charge for electro-neutrality */
  int *addchrs;        /* number of atoms over which the remaining charge will be distributed */
  sixstring **addchr;  /* their names */
  int **modif_extcharge;/* their atomnumbers */

  int connections;     /* connections to neighboring fragments */
  int *QMCA;           /* connection atom belonging to this site */
  int **QMCN;          /* connection atom used as capping, belonging to a neigboring site */

  int homos;           /* number of relevant (HO/LU)MOs for this site */
  int *homo;           /* list of orbital numbers of relevant MOs */
//  double *occupation;  /* occupation of individual MOs of this sites by the excess charge */
  double *hubbard;     /* hubbard parameter. different values for holes in different MOs. variation of orbital energy if an other orbital on this site is getting charged is mean value */
  double *lambda_i;    /* inner-sphere reorganization energy. different values for holes in different MOs */
  //
  double **delta_q;    /* difference of atomic charges between neutral and cationic site; arrays for holes in different MOs */


  int extcharges;      /* number of extcharges for this site */
  int *extcharge;      /* lists of extcharges */
  double **overlap;    /* overlap of site MOs between two steps */
  double **overlap_ref;/* overlap of site MOs between actual step and the beginning of the simulation */
  int nel;             /* total number of electrons */
  int norb;            /* total number of orbitals */
  double radius;       /* radius of charged sphere in polarizable continuum */
  int nocc;    //occupied orbs.=nel/2  //EX                                                      |  232   int nocc;
  int nvirt;   //virutal orbs.  //EX#                                                            |  233
  int nlow;    //number of occ. orbs. included in tddft calc.//EX                                |  234
  int *exn;    //list of excitations to use for exciton //EX                                     |  235
  int nexn;    //number of excitations
  //alex code
  int do_scc;          /* each site may or may not be calculated self consistently */
  int do_custom_occ;   /* each site may also use a customized occupation, instead of fermi distribution. This may be used as rough approximation for excited states. */
  double *custom_occ;  /* occupation vector e.g. ...2, 2, 2, 1, 1, 0, 0... for system where one electron from HOMO "gets excited" into LUMO*/
  double *com;         /* center of mass[bohr]. needed to decide if site should become active or not */
  int active;          /* switch 0/1 that determines if site is part of the QM calculation or just inactive member of the pool */
} ct_site_t;

typedef struct {
  double **U,   //unitary transformation matrix = alignment
         *sij,   //new jjk
         *evec,
         *work,
         *iwork,
         *eval,
         *issupz;
  long lwork,
       liwork;  //new jjk
} align_t;


typedef struct {
  /* main data structure of the charge transfer code */
  int jobtype;         /* cteXXX */
  int interval;        /* how often should the parameters be calculated? */
  int qmmm;            /* 1 or 0 */
  int delta_q_mode;    /* 0 use Mulliken charges, 1 use RESP charges to describe the charge carrier */
  int sitetypes;       /* number different types of sites e.g. adenine and guanine */
  ct_site_t *sitetype; /* list of different types of sites */
  int sites;           /* number of sites */
  ct_site_t *site;     /* list of the considered nucleobases */   // 		               NOW STRUCT BY ITS OWN
  int *do_scc;         /* each site may or may not be calculated self consistently */
  int do_scc_CG;       /* CG Hamiltonian may or may not be calculated self consistently */
  int dim;             /* dimension of hamilton matrix and so on. (not the same as the number of sites since multiple MOs per site are possible) */
  int is_hole_transfer; /* electron or hole transfer */
  double offdiag_scaling; /* 1.540 scaling factor according to J. Chem. Phys. 140, 104105 (2014) */
  int is_protein;      /* proteins need new residues built by protein_preprocessor. Not implemented yet*/
//  int *homo;           /* which orbital is HOMO? */				//	NOW IN STRUCT SITE
//  int *atoms;          /* number of atoms for each residue */				NOW IN STRUCT SITE
  int *last_atom;      /* index of the last atom of each residue, in the array of atoms of the complex */
//  int **atom;          /* lists of atoms - one array for each nucleobase */			NOW IN STRUCT SITE
                       /* first element is C1q while reading coordinates, */
                       /*               and is then converted to link hydrogen */
//  int **atomtype;      /* similarly as in the previous case */				NOW IN STRUCT SITE
                       /* C=0; H=1; N=2; O=3; */
  int atoms_cplx;      /* number of atoms for the complex */
  int *atom_cplx;      /* list of atoms for the complex */
                       /* contains C1q in place of link hydrogens */
  int *atomtype_cplx;  /* similarly as in the previous case */
                       /* C=0; H=1; N=2; O=3; */
  int extcharges_cplx; /* number of extcharges for the complex */
//  int **extcharge;     /* lists of extcharges - one array for each site */                         NOW IN STRUCT SITE
  int *extcharge_cplx; /* list of extcharges - for the complex */
//  twointegers *modif_extcharge;     /* which extcharges correspond to those of O4a and C2q? will be modified, +0.06080 */        NOW IN STRUCT SITE
  int modif_extcharges_cplx;        /* atomnumber of charges which will be modified */
  int *modif_extcharge_cplx;        /* their atomnumber */
//  double **delta_q;    /* difference of charge of neutral nucleobase and cation; array for each nucleobase */			NOW IN STRUCT SITE
  double *fo_shift;    /* shift that will be applied to the diagonal elements of the FO hamiltonian. can correct the difference between HOMO differences and IP differences */
  double **hamiltonian;/* CG Hamiltonian, calculated by DFTB, to be used in TDSE integration */
                       /* n x n matrix, fortran format */
  double ***hamiltonian_history;/* history over last n_avg_ham steps of CG Hamiltonian, to average fast (nonclassical) oscillations */
  int n_avg_ham;       /*length of hamiltonian_history */
  double *ev_spec;
  double *evec_spec;   /* used for hamilton-diagonalization to get CG-MO spectrum */
  double *work_spec;

  double *hamiltonian_mod; /* the contribution of Hubbard and electrostatics to the site energies - diagonal of hamiltonian[][] */
  double *hamiltonian_adiab; /* to be used in the adiabatic dynamics */
  double *tfl_mean_ham;  // new from jjk
  double **tfl_mean_ham_full; //new jjk
  double *tfl_old_ham;  //TFL new jjk
  double *tfl_old_ham_adiab; //new jjk
  int* indFO;          /* index of site i in Hamilton in FO-basis */
  double *ev_adiab;    /* to be used in the adiabatic dynamics */
  double *evec_adiab;  /* to be used in the adiabatic dynamics */
  double *work_adiab;  /* to be used in the adiabatic dynamics */
  double survival;     /* magnitude of the remaining excess charge - sum of the occupations */
  double *occupation;  /* occupation of individual MOs by the excess charge. essentially |wf|^2 ?????? */
  double *pop_adia;  /* adiabatic populations*/
  double **hubbard;    /* diag:    Hubbard parameter - chemical hardness of the nucleobase
                          offdiag: 1/r_ij - Coulomb interaction of parts of hole on two nucleobases */
  double sic;          /* naive self-interaction correction - scaling factor for the 2nd order term; 1.0 if no self-interaction correction */
  double esp_scaling_factor; /* divide the ESP (i.e. the DFTB shift) by this value - emulation of electronic polarization; 1.0 of no ESP scaling */
  int do_lambda_i;     /* shall the inner-sphere reorganization energy be emulated? 1==YES, 0==NO */
  int do_epol;         /* shall the electronic polarization be calculated? 1==YES, 0==NO */
  int do_projection;   /* do or do not project the wavefunction at every step onto the new FO basis before propagating */

  int neg_imag_pot;    /* shall the negative imaginary potential be applied? 0==NO, if > 0 then it means the number of sites to be applied on */
  int *site_neg_imag_pot; /* list of sites (numbered 1, 2, 3...) to apply the negative imaginary potential */
  double neg_imag_pot_rate_constant; /* in 1/a.u. of time */
  double *site_annihilated_occupation;
  double survival_threshold;
  int rk_neq;          /* Runge-Kutta: number of equations = 2 * number of sites */
  double rk_tol;       /* tolerance for Runge-Kutta */
  int adiabstart;      /* start SCC from lowest adiabatic state */
  double *wf;          /* CG wavefunction, dimension 2n: real[0] real[1] real[2] ... real[n-1] imag[0] imag[1] imag[2] ... imag[n-1] */
  double *wf_old;      /* CG wavefunction, before projection on new basis */
  double *wf_exc;
  double *q_act;
  double *q_old;
  double *dwf;         /* derivative of CG wavefunction */
  double *rk_ymax;     /* Runge-Kutta: auxiliary array */
  double *rk_thres;    /* Runge-Kutta: auxiliary array */
  double *rk_work;     /* Runge-Kutta: auxiliary array */
  int rk_lenwrk;       /* Runge-Kutta: length of rk_work = 32 * rk_neq */
  double rk_timestep;  /* time step for TDSE propagation in atomic units of time */
  int surface;         /* occupied surface for surface hopping */
  double *surf_overlap, *surf_massey, *surf_prob;
  double telec;       /* electronic temperature (in Kelvin) for the Fermi distribution */
  double fermi_kt;     /* kT(elec), for the Fermi distribution */
  /* variables for Tully's fewest switches - surface hopping */
  double *tfs_popul;       /* vector of population of adiabatic states */
  double *tfs_popul_der;   /* time derivative of tfs_popul [1/au_of_time] */
  double *tfs_diab;       //populations as diabatic wavefct. new  jjk
  double **tfs_vector;     /* adiabatic states phi_k(t+dt) */
  double **tfs_vector_old; /* adiab. states in the previous step phi_k(t) */
  double **tfs_overlap;    /* <phi_k(t) | phi_j(t+dt)> */
  double ***tfs_ncv;  /* nonadiabatic coupling vector*/
  double ***v1;  /* horizontal momenta */
  double ***v2;  /* vertical momenta */
  int tfs_initialization_step; /* is this the first step of a surface hopping simulation? */
  //variables for flexible surface hopping
  int tfl_num_of_states; //current number of sites in sys.
  int tfl_num_of_states_old;  //# of sites at last step
  int *tfl_is_in_system; //...[i]==1 if site i in sys., otherwise ...[i]==0
  int *tfl_is_in_system_old;  //at last step
  double tfl_ca_real;   //coef. of surface before integration (real part)
  double tfl_ca_im;    //(imaginary part)
/* variables for Landau-Zener surface hopping */
/* variables for Persico's local diabatic surface hopping */
  twodoubles **per_propag_operator; /* exp[-i*Z*dt] */
  twodoubles **per_transformator;   /* U = T^t * exp[-i*Z*dt] */
  ct_per_orthogo_t *per_arrays;
  double *tfs_popul_old;
  double *ev_adiab_old;
  double **per_diab_hamiltonian;
  int decoherence;     /* shall the decoherence correction be applied? 0==NO, 1==YES */
  ct_negf_lorentz_t *negf_arrays;
  int *sort_initialization_step; /*(new jjk) is this the first time one has to sort orbitals on site i? */
  align_t align; //new jjk
  dvec efield;  // external electric field
  int first_step; // is first time QM calculations? may differ from "step==0" from gromacs-steps in case of reruns. general purpose variable can be used by different routines (in contrast to tfs_initialization_step)
  double *born_overlap; //overlap used in cteBORNOPPENHEIMER
  int *incl_nonad;   //==1 if non-adiab coupling to be calculated, else ==0 //EX new jjk
  int use_bc_correction; //=1 if Boltzmann correction to be used in Ehrenfest//BC new jjk
  int use_strong_dec; //=1 if TFL&collapse after hopp new jjk

  // variables for adaptive QM zone
  int pool_size;          // number of possible sites (active sites are ct->sites)
  ct_site_t *pool_site;  // possible site (active ones are ct->site[i])
  dvec coc;              // center of charge[bohr] (determines which sites become active)
  dvec coc_start;        // center of charge[bohr] at the beginning of the simulation
  int opt_QMzone;        // switch 0/1 that derermines if optimal QM zone should be generated from pool of sites
  double adapt_inv_tot_mass; // inverse of the total mass of one site

  double **tfs_overlap2; /* eigenvectors at two time steps before */
  double **tfs_vector_old2; /* eigenvectors at two time steps before */
  double *ev_adiab_old2;  // adiabatic energy at two time steps before
  int ncv_isurf; // index for nonadiabtic coupling vecotr d_{ij}
  int ncv_jsurf; // index for nonadiabtic coupling vecotr d_{ij}
  int *ind_alzsh;
 // seed for random number
  int rnseed;

  double **fo_overlap;
  double **fo_tdc;

// implicite lambda
  double lambda;

// general machine learning parameters
  double hamiltonian_type;
  int neighbors_only; //for evaluating the couplings between neighboring sites (1D) only
  int pure_forces; //output force only, no weighing by wf coef. Debugging only
  int force_output; //enforce output to files independent of jobtype.
  double ***coords;
  double **hamiltonian_ml;
  double **hamiltonian_dftb;
  double ***diag_force;
  double ***offdiag_force;

// machine learning parameters for KRR
  long   ntrain[2];
  double sigma[2];
  double ***cm_train;
  double **alpha;

#if (GMX_TENSORFLOW)
// machine learning parameters for NN
  tf_model *site_model;
  tf_model *coupling_model;
#endif

// modified eherenfest methods
  int scctype;  // 0: SCC, 1: coherence penatly functional, 2: decoherence detailed balance
  double dephasing_rate;
// for pca
  double *coords_pca;
  double *coords_pca_tot;
  double *covar_aver;
  double **covar_tot;
  double *eigen_pca;
  double **eigenvec_pca;
  double *dummy_pca;
// coupling sign
  int define_orbital_sign;
  int atom_index_sign[3];
// fft for coupling
  double *fft_vij;

// store the signs of couplings
  double **positive_coupling_factor;
// average coupling
  double freeze_coupling ;
// time step
  int step;

  int **pair_index;
  double ***deriv_offdiag;

} charge_transfer_t;

typedef struct {
  /* data structure for DIIS charge convergence method used to converge coarse grained FO Hamiltonian */
  int n_prev_vector;       // actual # of stored previous vectors
  int n_elem;              // # of elements of vectors
  int indx;
  long ipiv[DIIS_MAX_PREV_VECTORS + 1]; // for dgesv (LAPACK)
  double **prev_q_input;   // stored previous input charges
  double **prev_q_diff;    // stored previous charge differences
  double *aa, *bb;         // auxiliary arrays for LAPACK
  double *q_inp_result;
  double *q_diff;
  double *fermi_coef;
} ct_diis_t;

typedef struct {
  /* data structure for BROYDEN charge convergence method used to converge coarse grained FO Hamiltonian */
  int i_iter;                  // actual iteration
  int n_elem;                  // # of elements of the vectors
  double omega0;               // Jacobi matrix differences
  double alpha;                // mixing parameter
  double minweight, maxweight; // min and max weights :-)
  double weight_fac;           // weighting factor (numerator)
  double *ww;                  // weights for previous iterations
  double *q_diff_last;         // charge difference in last iterations
  double *q_inp_last;          // input charge in last iteration
  double **aa;                 // storage for the matrix 'a'
  double *fifo_df;             // fifo for previous df vectors
  double *fifo_uu;             // fifo for previous u vectors
} ct_broyden_t;

typedef struct {
  /* data structure for BROYDEN charge convergence method used to converge fragment DFTB calculations */
  /* MAXSIZ will be substituted by the (much smaller) actual number of atoms in each nucleobase */
  double *f, //[MAXSIZ],
         *ui, //[MAXSIZ],
         *vti, //[MAXSIZ],
         *t1, //[MAXSIZ],
         // *vector, //[MAXSIZ][2],
         *dumvi, //[MAXSIZ],
         *df, //[MAXSIZ],
         a[IMATSZ_BROYDEN][IMATSZ_BROYDEN],
         b[IMATSZ_BROYDEN][IMATSZ_BROYDEN],
         b_lapack[IMATSZ_BROYDEN * IMATSZ_BROYDEN],
         cm[IMATSZ_BROYDEN],
         w[IMATSZ_BROYDEN],
         // *unit31, //[MAXSIZ][2],
         // *unit32, //[MAXSIZ][2][MAXITER],
         work[IMATSZ_BROYDEN * IMATSZ_BROYDEN];
  long ipiv[IMATSZ_BROYDEN];
  twodoubles *vector,
             *unit31;
  twodoubles_array *unit32;
} dftb_broyden_t;

typedef struct {
/* data structure needed for the orthogonalization of MOs in DFTB calculations of the complex ( run_dftb2() )*/
  double *tij,
         *sij,
         *evec,
         *work,
         *iwork,
         *eval,
         *issupz,
         **overlap,    /* overlap of basis functions between two steps */
         **overlap_ref,/* overlap of basis functions between actual and first step */
         **evec_ao,    /* orthogonal basis functions, expanded in AOs */
         **evec_ao_old,/* orthogonal basis functions of the previous step, expanded in AOs */
         **evec_ao_ref;/* orthogonal basis functions of first step, expanded in AOs */
  long lwork,
       liwork;
  /* this will be used to calculate derivatives of FMOs - Weiwei X Mar 2019 */
  double **fmo_old, // 1st index - AO, 2nd index - FO at preceding time step
         **sao_old, // symmetric matrix of overlap of AOs in the complex at preceding time step
         **fmo, // 1st index - AO, 2nd index - FO
         **sao, // symmetric matrix of overlap of AOs in the complex
         *sij_old;
} dftb_orthogo_t;

typedef struct {
  /* data structure for the DFTB calculations of the individual fragments (called "phase1")*/
  dvec *x, // NNDIM 3 /* coordinates */
       *x_old, //old coordinates, required for non_adiab coupling new jjk
       *grad, // for explicit lambda_i. force on every atom of this site. dim nn
       *partgrad, // for explicit lambda_i. force on every atom of this site. dim nn
       *xe, /* coordinates of external charges, NNDIM 3 */
       com, /* center of mass */
       *x_opt; /* optimized coordinates - static array in FORTRAN FORMAT (3,NNDIM) */
  double *mass, /* masses of atoms */
                /* attention - the link hydrogen is considered to have the mass of C1q */
                /* ... may be corrected in a future release */
         inv_tot_mass, /* inverse of the total mass of the nucleobase */
         esp, /* electro-static potential at the nucleobase */
         *ze, /* magnitudes of external charges, NNDIM */
         *qmat, // NNDIM  /* number of valence electrons on each atom in present SCC iteration */
         *qmold, // NNDIM /* number of valence electrons on each atom in last SCC iteration */
         *qmulli, // MDIM !!! /* help vector needed for mulliken charge calculation */
 //        *d_qmulli, // atomic delta q, calculated via mulliken charges
         au[LDIM][LDIM], /* auxiliary matrix in which the Hamilton matrix between two atoms is read in from the slater koster files */
         bu[LDIM][LDIM], /* auxiliary matrix in which the overlap matrix between two atoms is read in from the slater koster files */
         auh[LDIM][LDIM], //needed for gradient calculation (lambda_i)
         buh[LDIM][LDIM], //needed for gradient calculation (lambda_i)
        auhh[LDIM][LDIM], // new jjk
         buhh[LDIM][LDIM], // new jjk
         **a, // MDIM MDIM /* first this is the charge dependent (DFTB2) Hamilton matrix of the fragment (in AO basis). After diagonalization these are the eigenvectors*/
         **a_old, // MDIM MDIM /* eigenvectors of the last MD step */
         **a_ref, // MDIM MDIM /* used to store the eigenvectors of a reference snapshot (first MD step)*/
         **b, // MDIM MDIM /* overlap matrix */
         *a_trans, //ndim^2 /* essentially same as "a" but as single array, which is needed by LAPACK */
         *b_trans, //ndim^2 /* essentially same as "b" but as single array, which is needed by LAPACK */
         **hamil, // MDIM MDIM
         **overl, // MDIM MDIM
         **overl_diffuse, //overlap of atomic orbitals, calculatet with diffuse basis functions
         **gammamat, // nn nn
         *shift, // NNDIM
         *dshift, // NNDIM delta shift due to deltaQ used for QM calculation of lambda_i
         *shiftE, // NNDIM
         *ev, // MDIM
         *occ, // MDIM
         *aux, // 3*MDIM
         *pot, *pot2, *pot3, *pot4, // for PME - electric potential
         *homo_opt; /* optimized HOMO orbital - static array in FORTRAN FORMAT (MDIM) */
  int nn, /* number of atoms */
      *ind, // nn+1
      ndim, // = ind[nn]
      ne, /* number of external charges */
      *izp, /* list of atom types */
      nel, /* number of electrons */
      norb, /* number of orbitals == ndim!*/
      *neighbors_pme; /* for PME - number of neighbors for each atom */
  pme_integers *neighbor_pme; // for PME - neighbor lists for each atom
  real *q_pme; // for PME - charges
  rvec *x_pme; // for PME - coordinates of atoms and extcharges
  t_nrnb *nrnb_pme; // for PME, whatever it is...
  gmx_pme_t *pmedata; // for PME, various necessary data
} dftb_phase1_t;

typedef struct {
  /* data structure for the DFTB calculations of the whole complex consisting of several fragments  (called "phase2")*/
  dvec *x, // NNDIM 3 /* coordinates */
       *grad, // for explicit lambda_i. force on every atom of this site. dim nn
       *partgrad, // for explicit lambda_i. force on every atom of this site. dim nn
       *xe, /* coordinates of external charges, NNDIM 3 */
       com; /* center of mass */
  double *mass, /* masses of atoms */
         inv_tot_mass, /* inverse of the total mass of the complex */
         *ze, /* magnitudes of external charges, NNDIM */
         *qmat, // NNDIM
         au[LDIM][LDIM],
         bu[LDIM][LDIM],
         **hamil, // MDIM MDIM
         **overl, // MDIM MDIM overlap of atomic orbitals
         **b,
         auh[LDIM][LDIM], //needed for gradient calculation (lambda_i)
         buh[LDIM][LDIM], //needed for gradient calculation (lambda_i)
         **overl_c, // overlap of atomic orbitals with compressed basis (necessary for project_wf)
         **overl_hybrid, // overlap matrix with compressed value for intra-site overlap and uncompressed for inter-site
         **gammamat, // nn nn
         *shift, // NNDIM
         *shiftE, // NNDIM
         *ev, // MDIM
         *occ, // MDIM
         *aux, // 3*MDIM
         **Taf,
         **THamil, // new jjk
         **OverlF,
         **THamilOrtho,
         **tij,
         **sij,
         *pot, *pot2, *pot3, *pot4; // for PME - electric potential
  int nn, /* number of atoms */
      *ind, // nn+1
      *atind, // index of first atom of site i / remaining difference jjk has *at_ind
      ndim, // = ind[nn]
      *inf, // certain orbital index...
      *ihomo,  // ihomo[i] index of first homo of site i
      ne, /* number of external charges */
      *izp, /* list of atom types */
      nel, /* number of electrons */
      norb, /* number of orbitals == ndim!*/
      *neighbors_pme; /* for PME - number of neighbors for each atom */
  pme_integers *neighbor_pme; // for PME - neighbor lists for each atom
  real *q_pme; // for PME - charges
  rvec *x_pme; // for PME - coordinates of atoms and extcharges
  t_nrnb *nrnb_pme; // for PME, whatever it is...
  gmx_pme_t *pmedata; // for PME, various necessary data

  ////tddftb
  double **zeta;  //Zeta(A,B) fct., made to equal gammamat
} dftb_phase2_t;
typedef struct {
  //beware: for matrix no [i][j] dereference provided, use [i*dim+j] or flipped for col. major
  int nn;       //number of atoms   [whole block new in jjk]
  int nocc;     //#occupied orbitals=nel/2 (assume closed ground state!!)
  int nvirt;   //#virtual orbitals
  int nvpo;   //=nvirt+nocc
  int nexn;    //#included excited states on site
  int *exn;    //indices of included excited states
  int exdim;   //N_occ*N_virt, dimension of tddft problem
  int fexdim;  //(N_occ+N_virt)^2
  int nlow;   //number of included occupied orbitals
  int low;     //index of lowest lying state still included in TDDFT calc.
  //double **Cmat;  //matrix whose EV yield excited states in respone theory C*F=Omega^2*F
  double *Fmat;  //matrix whose EV yield excited states and after diagonalization matrix with EV
  double *Fmat_old; //old to determine sign, may for now only store relevant state not whole matrix
  double *ev_i;   //excitation energies
  double *w_i;   //difference of orbital ev for excitation
  double *qA; //Mulliken transition charges of atoms
  double *QAmb; //Many-body transition charges of atoms
  double *gammamat; //gamma matrix onsite SLK file
  double *t;  //Casida assignment rule many-body wf coefficients=X+Y
  ///for forces
  double *XpY;  //X+Y as defined in tddfrt
  double *XmY;  //X-Y
  double *XpY_t; //transposed for faster summation
  double *XmY_t;
  double *ApB;
  double *ApB_inv; //copy of A+B, used for matrix inversion
  double *qAf; //trans. mulliken charges for all pairs of orbitals up to highest virt.
  double *K;  //full K matrix
  double *T;
  double *P;  //difference density matrix
  double *QI;
  double *R;
  double *WI;
  double *dS;
  double *dH;
  double *dG;
  double *Dmn; //D_munu ground state density matrix (includes only included occupied orbitals)
  double *EDmn; //sum_i e_i*c_mui*c_nui for calculation of ground state gradient
  double *Pmn; //P_munu one-particle difference density matrix
  double *Wmn; //W_munu W in orbital basis
  double *XpYmn; //(X+Y)_munu
  int *ind2; //first index of orbtial pairs on atom i
  int norb2; //=SUM(norb(i)^2) over all atoms
  long *ipiv; //array for matrix diagonalization
  double **gammamatx; //gammamat with R_a=R_a+deltax
  double *dummyvec; //nn entries

} tddft_phase1_t;



typedef struct {
  /* main data structure for DFTB related stuff */
  int lmax[DFTB_MAXTYPES];       /* number of shells for each atom type */
  double racc, dacc;             /* machine accuracy */
  //double ****skhtab1, ****skstab1,
  //double (*(skhtab1[10]))[DFTB_MAXTYPES][DFTB_MAXTYPES], (*(skstab1[10]))[DFTB_MAXTYPES][DFTB_MAXTYPES],
  tendoubles *skhtab1[DFTB_MAXTYPES][DFTB_MAXTYPES], *skstab1[DFTB_MAXTYPES][DFTB_MAXTYPES];
  double skself1[DFTB_MAXTYPES][3], dr1[DFTB_MAXTYPES][DFTB_MAXTYPES],
    qzero1[DFTB_MAXTYPES], uhubb1[DFTB_MAXTYPES]; /* SLKO parameters for DFTB phase 1 */
  int dim1[DFTB_MAXTYPES][DFTB_MAXTYPES];
  //double ****skhtab2, ****skstab2,
  //double (*(skhtab2[10]))[DFTB_MAXTYPES][DFTB_MAXTYPES], (*(skstab2[10]))[DFTB_MAXTYPES][DFTB_MAXTYPES],
  tendoubles *skhtab2[DFTB_MAXTYPES][DFTB_MAXTYPES], *skstab2[DFTB_MAXTYPES][DFTB_MAXTYPES];
  double skself2[DFTB_MAXTYPES][3], dr2[DFTB_MAXTYPES][DFTB_MAXTYPES],
    qzero2[DFTB_MAXTYPES], uhubb2[DFTB_MAXTYPES]; /* SLKO parameters for DFTB phase 2 */
  int dim2[DFTB_MAXTYPES][DFTB_MAXTYPES];
  dftb_phase1_t *phase1;
  tddft_phase1_t *tddft_phase1;   //new jjk
  dftb_phase2_t phase2;
  dftb_broyden_t *broyden;
  dftb_orthogo_t orthogo;
  double *tddft_work; //memory for tddft diagonalization new jjk
  matrix box_pme;
  real ewaldcoeff_pme, rcoulomb_pme;
  int nstlist_pme, lastlist_pme;
  double ** overl_test; //just for testing stuff
  int **nl;  //neighbor list matrix of all QM atoms (dim atoms_cplx*atoms_cplx). is there any matrix element between them?
  int do_ncv;  // determine whether the nonadiabatic coupling vectors are calculated
} dftb_t;


#if GMX_MPI
void init_charge_transfer(t_atoms *atoms, const gmx_mtop_t *top_global, t_mdatoms *mdatoms, charge_transfer_t *ct, char *slko_path, t_state *state, int ct_mpi_rank);
//void init_dftb(t_mdatoms *mdatoms, dftb_t *dftb, charge_transfer_t *ct, char *slko_path, int ct_mpi_rank);
//void init_dftb_pme(dftb_t *dftb, charge_transfer_t *ct, t_inputrec *ir, int ct_mpi_rank);
//void do_neighborlist_for_dftb(charge_transfer_t *ct, dftb_t *dftb, rvec *x, int ct_mpi_rank, int ct_mpi_size);
#else
void init_charge_transfer(t_atoms *atoms, const gmx_mtop_t *top_global, t_mdatoms *mdatoms, charge_transfer_t *ct, char *slko_path, t_state *state);
//void init_dftb(t_mdatoms *mdatoms, dftb_t *dftb, charge_transfer_t *ct, char *slko_path);
//void init_dftb_pme(dftb_t *dftb, charge_transfer_t *ct, t_inputrec *ir);
//void do_neighborlist_for_dftb(charge_transfer_t *ct, dftb_t *dftb, rvec *x);
#endif
void init_dftb_stub(dftb_t *dftb, charge_transfer_t *ct);

void prepare_charge_transfer(matrix state_box, t_mdatoms *mdatoms, dftb_t *dftb, charge_transfer_t *ct, rvec *x_ct);
void ct_init_diis(charge_transfer_t *ct, ct_diis_t *diis);
void ct_init_broyden(charge_transfer_t *ct, dftb_broyden_t *broy);

#if GMX_MPI
void do_dftb_phase1(charge_transfer_t *ct, dftb_t *dftb, MPI_Comm ct_mpi_comm, int ct_mpi_rank, int ct_mpi_size);
void get_MM_params(charge_transfer_t *ct, dftb_t *dftb, MPI_Comm ct_mpi_comm, int ct_mpi_rank, int ct_mpi_size);
void do_tddft_phase1(charge_transfer_t *ct, dftb_t *dftb, MPI_Comm ct_mpi_comm, int ct_mpi_rank, int ct_mpi_size);
void prep_exc_state_forces(charge_transfer_t *ct, dftb_t *dftb, MPI_Comm ct_mpi_comm, int ct_mpi_rank, int ct_mpi_size);
void calc_tddft_forces(charge_transfer_t *ct, dftb_t *dftb, MPI_Comm ct_mpi_comm, int ct_mpi_rank, int ct_mpi_size);
void do_esp_only(charge_transfer_t *ct, dftb_t *dftb, real *q, MPI_Comm ct_mpi_comm, int ct_mpi_rank, int ct_mpi_size);
#else
void do_dftb_phase1(charge_transfer_t *ct, dftb_t *dftb);
void get_MM_params(charge_transfer_t *ct, dftb_t *dftb);
void do_tddft_phase1(charge_transfer_t *ct, dftb_t *dftb);// new jjk
void do_esp_only(charge_transfer_t *ct, dftb_t *dftb, real *q);
#endif
void do_dftb_phase2(charge_transfer_t *ct, dftb_t *dftb);
void calc_zeta(charge_transfer_t *ct, dftb_t *dftb);  //for tddft exciton treatment  new jjk

/*
#if GMX_MPI
void do_pme_for_dftb_part1(charge_transfer_t *ct, dftb_t *dftb, int ct_mpi_rank, int ct_mpi_size);
#else
void do_pme_for_dftb_part1(charge_transfer_t *ct, dftb_t *dftb);
#endif
void do_pme_for_dftb_part2(charge_transfer_t *ct, dftb_t *dftb, int i);
void do_pme_for_dftb_phase2(charge_transfer_t *ct, dftb_t *dftb);

int run_dftb1(charge_transfer_t *ct, dftb_t *dftb, int ibase);
int run_esp_only(charge_transfer_t *ct, dftb_t *dftb, int ibase);
int run_dftb2(charge_transfer_t *ct, dftb_t *dftb);
//int run_dftb2_mpi(charge_transfer_t *ct, dftb_t *dftb, MPI_Comm ct_mpi_comm, int ct_mpi_rank, int ct_mpi_size, MPI_Status ct_mpi_status );
*/

void ct_assemble_hamiltonian(charge_transfer_t *ct, dftb_t *dftb);
/* void ct_integrate_tdse(charge_transfer_t *ct, double twant); */

/* do_rksuite.cpp */
int do_rksuite(charge_transfer_t *ct);
int do_rksuite_diab(charge_transfer_t *ct); //new jjk
int do_rksuite_tfs(charge_transfer_t *ct);

int do_fermi_surface_hopping(charge_transfer_t *ct, dftb_broyden_t *broyd, double *fermi_coeff, FILE *f);
//int do_tully_fewest_switches(charge_transfer_t *ct, dftb_broyden_t *broyd, double *fermi_coeff, FILE *f, FILE *f2);
//int do_tully_diab(charge_transfer_t *ct, FILE *f, FILE *f2); // new jjk
int do_surf_hop_dummy(charge_transfer_t *ct, FILE *f, FILE *f2);  //dummy function to play around with different methods new jjk
//int do_persico_diabatic_sfhopping(charge_transfer_t *ct, dftb_broyden_t *broyd, double *fermi_coeff, FILE *f, FILE *f2);
//int do_persico_diabatic_sfhopping_new(charge_transfer_t *ct, dftb_broyden_t *broyd, double *fermi_coeff, FILE *f, FILE *f2, FILE *f3);
long exp_imag_matrix(double **in, twodoubles **out, double dt, long n, ct_per_orthogo_t *arrays);

void check_and_invert_orbital_phase(dftb_phase1_t *dftb1, charge_transfer_t *ct, t_state *state_global, t_mdatoms *mdatoms);

int negf_init_arrays(ct_negf_lorentz_t *negf, double *rk_timestep, double *wf);
int negf_propagate(charge_transfer_t *ct);

//Alex
//int searchkey(int lines, char input[MAXLINES][2][MAXWIDTH], char *key , char value[MAXWIDTH], int required);
t_atoms gmx_mtop_global_atoms(const gmx_mtop_t *mtop);
void get_delta_q(dftb_t *dftb, charge_transfer_t *ct, int i);
void get_internal_forces(dftb_t *dftb, charge_transfer_t *ct, int site_i);
void write_out_MOs(int step,rvec x_ct, t_atoms *ct_atoms, dftb_t *dftb, charge_transfer_t *ct);
void get_spectrum(charge_transfer_t *ct, dftb_t *dftb);
void sort_mobasis(dftb_t *dftb, charge_transfer_t *ct, int i);
void get_alignment_matrix(align_t arrays, double **overlap, long n );      //new jjk
t_atoms* protein_preprocessor(t_atoms *atoms, t_state *state_global);
//void print_time_difference(const char *s, struct timespec start, struct timespec end);
void usual_gradient_homo(dftb_t *dftb, dvec *x, dvec *grad, charge_transfer_t *ct, int site_i); // this will be needed for actual calculation of lambda_i
void offdiag_gradient_homo(dftb_t *dftb, dvec *x, dvec *grad, charge_transfer_t *ct);
void gamma_gradient_homo(dftb_t *dftb, dvec *x, dvec *grad, charge_transfer_t *ct, int site_i);
void additional_gradient_homo(dftb_t *dftb, dvec *x, dvec *grad, charge_transfer_t *ct, int site_i);
void additional_gradient_homo_new(dftb_t *dftb, dvec *x, dvec *grad, charge_transfer_t *ct, int site_i);
//int do_prezhdo_sfhopping(charge_transfer_t *ct, dftb_broyden_t *broyd, double *fermi_coeff, FILE *f, FILE *f2, FILE *f3);
//void adapt_QMzone(charge_transfer_t *ct, rvec *x_ct); //new jjk
void search_starting_site(matrix state_box, t_mdatoms *mdatoms, dftb_t *dftb, charge_transfer_t *ct, rvec *x_ct, char *slko_path, const gmx_mtop_t *top_global, rvec *gromacs_x);
int adapt_QMzone(charge_transfer_t *ct, rvec *x_ct, t_mdatoms *mdatoms, const gmx_mtop_t *top_global, matrix state_box, rvec *gromacs_x);
int find_intersection(int size, int array1[], int array2[], int intersection_array[]);

#if (GMX_TENSORFLOW)
void NoOpDeallocator(void* data, size_t a, void* b);
void initialize_nn_model(tf_model *model, const char *path, const char *hyperpath);
void input_data_to_model(tf_model *model, float data[][3], int natom);
void check_nn_status(tf_model *model);
void get_neighbor_list(charge_transfer_t *ct, dftb_t *dftb);
void get_nn_hamiltonian(charge_transfer_t *ct, t_state *state_global, t_mdatoms *mdatoms, dftb_t *dftb);
#endif

/* communication.cpp */
void transfer_collect_x(const t_commrec* cr,
                        t_state*         state_local,
                        t_state*         state_global);

/* get_mm_params.cpp */
#if GMX_MPI
void get_MM_params(charge_transfer_t *ct, dftb_t *dftb, MPI_Comm ct_mpi_comm, int ct_mpi_rank, int ct_mpi_size);
#else
void get_MM_params(charge_transfer_t *ct, dftb_t *dftb);
#endif

/* lapack.cpp */
long dsyev(int int_n, double *a, double *w, double *work, int int_lwork);
long dgesv(int int_n, double *a, double *b, long ipiv[DIIS_MAX_PREV_VECTORS + 1]);

/* negf* */
long negf_const_mat_(double *, double *, long *, long *);
long negf_jacobi_(double *, double *, double *, long *);
long negf_construct_nu_(double *, double *, double *, double *, long *, double *, double *, double *, double *, long *);
long negf_calc_r_(double *, double *, double *, long *, double *, double *, double *, long *);
long negf_create_hi_(double *, double *, double *, double *, long *, double *, double *, double *, double *, long *, double *);
long negf_create_gam_(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, long *, long *, double *);
long negf_set_initial_values_(double_complex *, double *);
long negf_allocate_fortran_(long *, long *, long *, long *);
long negf_rk_init_();
void negf_create_h_(double *);
void negf_do_step_(double_complex *, double_complex *, double *);
void negf_calculate_current_(double_complex *, double *, double *);

/* quantum_auxiliary.cpp */
double fermi_coef_sum(double x, int n, double *a, double *coeffs, double fermi_kt);
void assignprob(charge_transfer_t *ct, int *ind);

/* quantum_broyden.cpp */
void broyden(int niter, double alpha, int jtop, double *vecin, double *vecout, dftb_broyden_t *arrays);

/* quantum_driver.cpp*/
double calc_tda(charge_transfer_t *ct);
int do_adiabatic(charge_transfer_t *ct, ct_diis_t *diis, FILE *f);
int do_surface_hopping(charge_transfer_t *ct, ct_diis_t *diis, FILE *f);
double do_born_oppenheimer(charge_transfer_t *ct, dftb_broyden_t *broyd, double *fermi_coeff, FILE *f);
int do_adiab_fermi(charge_transfer_t *ct, dftb_broyden_t *broyd, double *fermi_coeff, FILE *f);
int do_adiab_fermi_onestate(charge_transfer_t *ct, dftb_broyden_t *broyd, double *fermi_coeff, FILE *f);
int do_lzsh_diab(int step, charge_transfer_t *ct, dftb_t *dftb, t_state *state_global, t_mdatoms *mdatoms, rvec *x_ct, FILE *f);
int do_lzsh_adiab(int step, charge_transfer_t *ct, FILE *f);
int do_fssh_diab(int step, charge_transfer_t *ct, dftb_t *dftb, t_state *state_global, t_mdatoms *mdatoms, FILE *f);
#if GMX_MPI
int do_fssh_johen(int step, charge_transfer_t *ct, dftb_t *dftb, t_state *state_global, t_mdatoms *mdatoms, FILE *f,
                  MPI_Comm ct_mpi_comm, int ct_mpi_rank, int ct_mpi_size);
int do_self_consistent_restricted_decoherence_fssh(int step, charge_transfer_t *ct, dftb_t *dftb, t_state *state_global, t_mdatoms *mdatoms, FILE *f,
                                                   MPI_Comm ct_mpi_comm, int ct_mpi_rank, int ct_mpi_size);
int do_global_flux_fssh(int step, charge_transfer_t *ct, dftb_t *dftb, t_state *state_global, t_mdatoms *mdatoms, FILE *f,
                        MPI_Comm ct_mpi_comm, int ct_mpi_rank, int ct_mpi_size);
int do_crossing_corrected_fssh(int step, charge_transfer_t *ct, dftb_t *dftb, t_state *state_global, t_mdatoms *mdatoms, FILE *f,
                               MPI_Comm ct_mpi_comm, int ct_mpi_rank, int ct_mpi_size);
int do_decoherence_induced_fssh(int step, charge_transfer_t *ct, dftb_t *dftb, t_state *state_global, t_mdatoms *mdatoms, FILE *f,
                                MPI_Comm ct_mpi_comm, int ct_mpi_rank, int ct_mpi_size);
#else
int do_fssh_johen(int step, charge_transfer_t *ct, dftb_t *dftb, t_state *state_global, t_mdatoms *mdatoms, FILE *f);
int do_self_consistent_restricted_decoherence_fssh(int step, charge_transfer_t *ct, dftb_t *dftb, t_state *state_global, t_mdatoms *mdatoms, FILE *f);
int do_global_flux_fssh(int step, charge_transfer_t *ct, dftb_t *dftb, t_state *state_global, t_mdatoms *mdatoms, FILE *f);
int do_crossing_corrected_fssh(int step, charge_transfer_t *ct, dftb_t *dftb, t_state *state_global, t_mdatoms *mdatoms, FILE *f);
int do_decoherence_induced_fssh(int step, charge_transfer_t *ct, dftb_t *dftb, t_state *state_global, t_mdatoms *mdatoms, FILE *f);
#endif
int do_tully_local(int step, charge_transfer_t *ct, dftb_t *dftb, t_state *state_global, t_mdatoms *mdatoms, FILE *f);
int do_persico_diabatic_sfhopping(charge_transfer_t *ct, dftb_broyden_t *broyd, double *fermi_coeff, FILE *f, FILE *f2);
int do_persico_diabatic_sfhopping_new(charge_transfer_t *ct, dftb_broyden_t *broyd, double *fermi_coeff, FILE *f, FILE *f2, FILE *f3);
int do_prezhdo_sfhopping(charge_transfer_t *ct, dftb_broyden_t *broyd, double *fermi_coeff, FILE *f, FILE *f2, FILE *f3);

/* wf_project_or_decomp.cpp */
int do_wf_decomp_evec(charge_transfer_t *ct, ct_diis_t *diis, FILE *f);
int do_wf_decomp_evec_fermi(charge_transfer_t *ct, dftb_broyden_t *broyd, double *fermi_coeff, FILE *f);
void project_wf_on_new_basis(int step, dftb_t *dftb, charge_transfer_t *ct, FILE *f_ct_project_wf, FILE *f_ct_project_wf_ref);
void project_wf_on_new_basis_exact(int step, dftb_t *dftb, charge_transfer_t *ct, FILE *f_ct_project_wf, FILE *f_ct_project_wf_ref);

