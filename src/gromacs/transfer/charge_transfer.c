/*************************
 * Charge transfer in DNA
 * Tomas Kubar
 * Neural Nets: Mila Kraemer
 *************************/

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <tensorflow/c/c_api.h>
#include "mdatoms.h"

#include "charge_transfer.h"
//#include <Python.h>
#include <math.h>

#ifdef GMX_MPI
#define PRINTF(...) if (ct_mpi_rank == 0) printf(__VA_ARGS__)
#else
#define PRINTF(...) printf(__VA_ARGS__)
#endif
#define CUB(x) ((x)*(x)*(x))

void ct_get_index(int isize[], atom_id *index[])   /* atom_id is typedef int! */
{ /* these were parameters to get_index() - we do not need them here */
    char *fnm = "charge-transfer.ndx";
    /* maybe unnecessary?
       int ngrps = 1;
       char *grpnames;
     */
    char    ***gnames;
    t_blocka  *grps = NULL;
    int        j;

    snew(gnames, 1);

    grps = init_index(fnm, gnames);

    if (grps->nr != 1)
    {
        fprintf(stderr, "\nThe index file %s contains a number of groups different from 1! (%d)\nExiting!\n\n", fnm, grps->nr);
        exit(-1);
    }

    /* ct_read_group(grps,*gnames,isize,index); */

    fprintf(stderr, "   -- group %s has %d atoms\n", **gnames, grps->index[1] - grps->index[0]);

    isize[0] = grps->index[1]-grps->index[0];
    snew(index[0], isize[0]);
    for (j = 0; j < isize[0]; j++)
    {
        index[0][j] = grps->a[grps->index[0]+j];
    }

    return;
}

int ct_atom_in_group(int atom, int *list, int list_size)
{
    int i;

    for (i = 0; i < list_size; i++)
    {
        if (atom == list[i])
        {
            return 1;
        }
    }

    return 0;
}

int negf_init_arrays(ct_negf_lorentz_t *negf, double *rk_timestep, double *wf)
{
    extern long negf_const_mat_(double *, double *, long *, long *);
    extern long negf_jacobi_(double *, double *, double *, long *);
    extern long negf_construct_nu_();
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


int searchkey(int lines, char input[MAXLINES][2][MAXWIDTH], char *key, char value[MAXWIDTH], int required)
{
/* given a keyword, this function checks if it is present in the input. If it is so, the corresponding value will be returned */

// PARAMETERS:
// lines = (in) total number of lines of input
// key   = (in) the key that is searched
// value = (out) the corresponding value to key
// required = (in) is the keyword required or optional for this type of calculation
////////////////////////////////////////////////////////////////////////
    int line;
    for (line = 0; line < lines; line++)
    {
        if (strcmp(input[line][0], key) == 0)
        {
            strcpy(value, input[line][1]);
            return 1; //keyword found
        }
    }
    if (required)
    {
        printf("KEYWORD %s NOT FOUND!\n", key);
        exit(-1);
    }
    else
    {
        return 0; //keyword not found
    }
}


int read_file(char* file_name, char input[MAXLINES][2][MAXWIDTH], char possiblekeys[][2][MAXWIDTH], int nkeys)
{
/* reads the input files that are related to the CT code */

// PARAMETERS:
// file_name = (in) complete path to file
// input     = (out) the key and value pairs that are found in the input files
// possiblekeys = (in) collection of every allowed key
// nkeys     = (in) number of possible keys
///////////////////////////////////////////////////////
    int   i, j;
    int   lines, ch, line, len;
    FILE *f;
    char *s;
    //char input[MAXLINES][2][MAXWIDTH];

    f = fopen(file_name, "r");
    if (f == NULL)
    {
        printf("%s not accessible, exiting!\n", file_name);
        exit(-1);
    }
    printf("Searching in file %s for:\n", file_name);
    for (i = 0; i < nkeys; i++)
    {
        printf("%20s:  %s\n", possiblekeys[i][0], possiblekeys[i][1]);
    }

    /* get length of input file */
    lines = 0;
    while (!feof(f))
    {
        ch = fgetc(f);
        if (ch == '\n')
        {
            lines++;
        }
    }
    //rewind(f);
    fclose(f);
    f = fopen(file_name, "r");


    /* read input file */
    for (line = 0; line < lines; line++)
    {
        //ch=fscanf(f, "%[;a-zA-Z0-9/] = %[-a-zA-Z0-9/{}. ]\n", input[line][0],  input[line][1]);
        ch = fscanf(f, " %[^= \t\r] = %[^=\n\t\r] \n", input[line][0],  input[line][1]);
        if (ch != 2)
        {
            printf("READING ERROR in line %d.\n Use format 'keyword = value'.\n ", line+1); exit(0);
        }

        //remove trailing whitespace
        s = input[line][1] + strlen(input[line][1]);
        while (--s >= input[line][1])
        {
            if (!isspace(*s))
            {
                break;
            }
            *s = 0;
        }

        //convert all to lowercase for better handling (besides file and atom names)
        len = strlen(input[line][0]);
        for (i = 0; i < len; i++)
        {
            input[line][0][i] = tolower((unsigned char)input[line][0][i]);
        }
        if ( (strcmp(input[line][0], "slkopath") && strcmp(input[line][0], "specfiles") && strncmp(input[line][0], "name", 4) ) )
        {
            len = strlen(input[line][1]);
            for (i = 0; i < len; i++)
            {
                input[line][1][i] = tolower((unsigned char)input[line][1][i]);
            }
        }
        // check input file
        j = 1;
        for (i = 0; i < nkeys; i++)
        {
            if (strcmp(input[line][0], possiblekeys[i][0]) == 0)
            {
                j = 0;
            }
        }
        if (j)
        {
            printf("KEYWORD NOT KNOWN: %s\n", input[line][0]);
            exit(-1);
        }
    }
    fclose(f);
    printf("Finished reading file %s\n", file_name);

    /* print copy of the input file */
    printf("---------------------\n");
    printf("INPUT:\n");
    for (line = 0; line < lines; line++)
    {
        printf("%s = %s\n", input[line][0], input[line][1]);
    }
    printf("---------------------\n");
    return lines;
}


int split_string_into_double(char value[MAXWIDTH], int n, double* target, char* name)
{
/* splits a string of whitespace separated numbers and converts them into double data type */

// PARAMETERS:
// value  = (in) is the value string that gets splitted
// n      = (in) expected number of elements
// target = (out) vector containing the individual doubles
// name   = (in) name of the keyword whose corresponding value gets splitted
//////////////////////////////////////////////////////////////////////////////////////////
    int   i;
    char *ptr;

    i   = 0;
    ptr = strtok(value, " ");
    while (ptr != NULL)
    {
        target[i] = atof(ptr);
        ptr       = strtok (NULL, " ");
        i++;
    }
    if (i > n)
    {
        printf("TOO MANY ARGUMENTS FOR %s\n", name); exit(-1);
    }
    else if (i < n)
    {
        printf("TOO FEW ARGUMENTS FOR %s\n", name); exit(-1);
    }
    else
    {
        return 0;
    }
}


int split_string_into_int(char value[MAXWIDTH], int n, int* target, char* name)
{
/* splits a string of whitespace separated numbers and converts them into integer data type */

// PARAMETERS:
// value  = (in) is the value string that gets splitted
// n      = (in) expected number of elements
// target = (out) vector containing the individual inegers
// name   = (in) name of the keyword whose corresponding value gets splitted
//////////////////////////////////////////////////////////////////////////////////////////
    int   i;
    char *ptr;

    i   = 0;
    ptr = strtok(value, " ");
    while (ptr != NULL)
    {
        target[i] = atoi(ptr);
        ptr       = strtok (NULL, " ");
        i++;
    }
    if (i > n)
    {
        printf("TOO MANY ARGUMENTS FOR %s\n", name); exit(-1);
    }
    else if (i < n)
    {
        printf("TOO FEW ARGUMENTS FOR %s\n", name); exit(-1);
    }
    else
    {
        return 0;
    }
}


int split_string_into_string(char value[MAXWIDTH], int n, char target[][ELEMENTWIDTH], char* name)
{
/* splits a string of whitespace separated substrings (names, words,...)*/

// PARAMETERS:
// value  = (in) is the value string that gets splitted
// n      = (in) expected number of elements
// target = (out) vector containing the individual substrings
// name   = (in) name of the keyword whose corresponding value gets splitted
//////////////////////////////////////////////////////////////////////////////////////////
    int   i;
    char *ptr;

    i   = 0;
    ptr = strtok(value, " ");
    while (ptr != NULL)
    {
        strcpy(target[i], ptr);
        ptr = strtok (NULL, " ");
        i++;
    }
    if (i > n)
    {
        printf("TOO MANY ARGUMENTS FOR %s\n", name); exit(-1);
    }
    else if (i < n)
    {
        printf("TOO FEW ARGUMENTS FOR %s\n", name); exit(-1);
    }
    else
    {
        return 0;
    }
}



/***************************************
 * INITIALIZE THE CHARGE TRANSFER CODE *
 ***************************************/

#ifdef GMX_MPI
void init_charge_transfer(t_atoms *atoms, gmx_mtop_t *top_global, t_mdatoms *mdatoms, charge_transfer_t *ct, char *slko_path,  t_state *state, int ct_mpi_rank)
{
#else
void init_charge_transfer(t_atoms *atoms, gmx_mtop_t *top_global, t_mdatoms *mdatoms, charge_transfer_t *ct, char *slko_path, t_state *state)
{
#endif
/* Initializes the charge transfer code. Reads in CT files, sets up the stuff that is fixed variables, allocates memory */

// PARAMETERS:
// atoms       = (in) variable of GROMACS data-structure-type, which we generated outside of this function. Contains info like atomnames and residuenumbers, that we need to distinguish QM from MM atoms.
// top_global  = (in) GROMACS data structure, from which we only need total number of MM atoms
// mdatoms     = (in) variable of GROMACS data-structure-type (we need atomic masses)
// ct          = (out) main data structure with information about the charge transfer calculation
// slko_path   = (out) path to the Slater-Koster files
// state       = (in) GROMACS data structure that holds (among others) the coordinates of the MD system
// ct_mpi_rank = (in) the rank of the process in MPI calculations (each process in MPI calculations can be destinguished by its value of ct_mpi_rank)
//////////////////////////////////////////////////////////////////////////////////////////

    char sitespecdat[MAXSITETYPES][ELEMENTWIDTH], value[MAXWIDTH];
    char possiblekeys1[][2][MAXWIDTH] = {
        {"seed",  " seed for random number generator"},
        {"hamiltonian",  " hamiltonian calculated by dftb or machine learning"},
        {"sitemodelpath",  "{Path} to the NN model for site energy"},
        {"couplingmodelpath",  "{Path} to the NN model for coupling"},
        {"siteparapath",  "{Path} to the hyperparameters for site energy"},
        {"couplingparapath",  "{Path} to the hyperparameters for coupling"},
        {"neighborsonly", "How often to calculate Neighborlist: -1 (once at start), 0 (never, evaluate all interactions), x (every x-th step)"},
        {"pureforces", "Calculate forces without weighing with wf coefficients (for debugging)"},
        {"forcedoutput", "Output CT data to file independently of jobtype"},
        // {"sitetrainpath",  "{Path} to the training data for site energy"},
        // {"couplingtrainpath",  "{Path} to the training data for coupling"},
    	// {"sigmacouplings",  "value of sigma for couplings-ml-model"},
	    // {"sigmasites",  "value of sigma for sites-ml-model"},
    	// {"ntraincouplings",  "number of training data for couplings-ml-model"},
    	// {"ntrainsites",  "number of training data for sites-ml-model"},
        {"frcoupling",  "freeze the electronic coupling to be average coupling"},
        {"slkopath",  "{Path} to directory of the DFTB Slater-Koster files"},
        {"chargecarrier",   "The charge carrier (electron/hole). Will effect sign of the Hamilton matrix and of the charges that are added to the force-field."},
        {"offdiagscaling",   "{yes/no} Scale offdiagonal elements of FO Hamiltonian. See: J. Chem. Phys. 2014, 140, 104105+  and  Phys. Chem. Chem. Phys. 2015, 17, 14342-14354."},
        {"atomindex",   "Atomic indexes to fix initial sign of electronic coupling. "},
        {"extchrmode",   "Treatment of the MM pointcharges. {vacou} is as it says, {qmmm} uses pointcharges with minimum image convention, {pme} is particle-mesh-Ewald treatment"},
        {"espscaling",  "Scales the strength of the electrostatic potential of the environment with 1/espscaling."},
        {"efield", "External electric field vector [V/cm]. X Y and Z direction. Adds shift to the site energies depending on their position."},
        {"nsitetypes", "Number of unique molecules."},
        {"typefiles", "File with specifications for every unique fragment."},
        {"nsites", "Total number of sites that are treated at the QM level."},
        {"zonesize", "Select QM zone as subset of nstites. Default is zonesize=nsites."},
        {"optimizezone", "If {yes/no} construct QM zone of zonesize fragments around site with lowest energy of all nsites fragments."},
        {"sites", "Residue number of all nsites fragments"},
        {"sitetypes", "Type of every site. 1,2,3,etc. corresponding to the order of typefiles."},
        {"sitescc", "{0} Non-self-consistent DFTB1 calculations. {1} Self-consistent DFTB2 calculations. {2} DFTB2 calculations where initial charges are taken from last MD step to accelerater convergence."},
        {"foshift", "Shift that is added to the diagonal elements of the FO hamiltonian [hartree]. This can correct for wrong relative energies due to approximating ionization potentials with HOMO energies"},
        {"jobtype", "{PAR}: calculate Hamiltonmatrix along MD. Also possible  with -rerun option. {SCC}: propagate chargecarrier and nuclei in an Ehrenfest simulation. {TDA}: calculate bridge mediated coupling between FO 1 and FO N with bridge states 2...(N-1)"},
        {"nstqm", "Frequency of the QM calculations for jobs whithout propagation. 1 is every step, 2 is every other step etc."},
        {"tfermi", "Fermi Temperature [K] for certain jobs."},
        {"epol", "Electronic polatization model. {imp} is naive implicit polarization (like born-model) for the charge carrier."},
        {"sic", "Self interaction correction factor. Scales second order terms. See: J. Phys. Chem. B 2010, 114, 11221-11240."},
        {"internalrelax", "Internal relaxation of the sites. {parameter} uses precalculated values. {onsite} relaxes each site according to DFTB forces. {full} additionally calculates inter-molecular forces."},
        {"adiabstart", "{yes/no} Take lowest eigenvector of the FO-Hamiltonian as starting wavefunction."},
        {"wavefunctionreal", "Coefficients of the starting Wavefunction. Real part. Default is vector of zeroes. If non-normalized wavefunction is provided, lowest adiabatic eigenvector will be taken instead."},
        {"wavefunctionim", "Imaginary part of the wavefunction coefficients. Default is vector of zeroes."},
        {"nnegimpot", "Negative imaginary potential. Can be used to annihilate charge."},
        {"negimpotrate", "Rate for charge annihilation."},
        {"negimpotfos", "FOs from which the charge will be annihilated."},
        {"deltaqmode", "Either add precalculated RESP charges to the force-field to describe the charge-carrier {resp} or use internally obtained DFTB Mulliken-charges {mulliken}"},
        {"projection", "{yes/no} do or do not project the charge carrier wavefunction at every step onto the new FO basis before propagating it."}
    };

    char possiblekeys2[][2][MAXWIDTH] = {
        {"natoms", "Number of QM atoms of this site (inclusive capping atoms)."},
        {"nelectrons", "Total number of valence electrons of this site."},
        {"nallorbitals", "Sum of occupied and virtual MOs of this site (note that there are only valence electrons in DFTB)"},
        {"radius", "Radius for spherical approximation of the site. Needed in Born solvation model."},
        {"nqmmmbonds", "Number of bonds that have to be capped."},
        {"nignorechr", "For each bond the number of MM atoms (besides the MM link atom) whose charge will be deleted. Use 'nignorechr=1 2' to delete 1 additional charge for the first and 2 for the second bond"},
        {"nameignorechr", "Atomnames for nignorechr. Use of 'nameignorechr=H1 H3 H4' with nqmmmbonds=2 and 'nignorechr=1 2' will delete H1 for the first bond and H3 and H4 for the second bond"},
        {"naddchr", "For each bond 'totaladdchr' will be distributed over 'naddchr' atoms to restore the integer total charge of the environment."},
        {"nameaddchr", "Atomnames for naddchr. Use of 'nameaddchr=C1 C3 C4' with 'nqmmmbonds=2' , 'naddchr=1 2' and 'totaladdchr=-0.1 0.05' will put -0.1 on C1 and distibute 0.05 over C3 and C4"},
        {"totaladdchr", "Charge for each QM/MM bond to restore the integer total charge of the environment."},
        {"nfragorbs", "Number of molecular orbitals of this site that will be used in the fragment orbital Hamiltonian."},
        {"fragorbs", "The index (starting from 1) of the molecular orbitals."},
        {"hubbard", "Hubbard parameter for each orbital."},
        {"lambda_i", "Lambda_i for each orbital."},
        {"dqresp", "List of RESP charges that will be added to the QM atoms, scaled by the occupation of the HOMO/LUMO. If more than one HOMO is used per site, first the list for the first FO is read then for the second."},
        {"customocc", "list with nallorbitals elements defining a fixed occupation (e.g. 2 2 2 1 1 0 0 0)"}
    };


    int        i, j, k, l, m, n, counter, *counter_array, counter_cplx, QMLAcounter, MMLAcounter, QMCAcounter, QMCApool[QMCASIZE], QMCApoolcounter, environment, environment_cplx, modif_cplx, counter_modif_cplx = 0,
               mm_list_size, *mm_list;
    double     sum, bond_length, bond_length_best, X[3], Y[3], magnitude, mass;
    dvec       bond;
    ct_site_t *site, s;

    char       input1[MAXLINES][2][MAXWIDTH], input2[MAXLINES][2][MAXWIDTH]; // variables to store input files(charge-transfer.dat, and *.spec): several lines, 2 collumns (before and after '='), and string of several chars
    char       dummy[MAXELEMENTS][ELEMENTWIDTH];
    int        lines1, lines2;
//  names of files for machinelearning, W. X + M. K
    char       ml_coupling_model_path[MAX_PATH_LENGTH], ml_site_energy_model_path[MAX_PATH_LENGTH];
    char       ml_site_energy_para_path[MAX_PATH_LENGTH], ml_coupling_para_path[MAX_PATH_LENGTH];
    // char       ml_coupling_train_path[MAX_PATH_LENGTH], ml_site_energy_train_path[MAX_PATH_LENGTH];
    // char       ml_site_energy_para_path[MAX_PATH_LENGTH], ml_coupling_para_path[MAX_PATH_LENGTH];

    ct->first_step = 1;
#ifdef GMX_MPI
    printf("Initializing charge transfer at rank %d\n", ct_mpi_rank);
#else
    printf("Initializing charge transfer\n");
#endif



    /* Get input from file charge-transfer.dat */
    lines1 = read_file("charge-transfer.dat", input1, possiblekeys1, sizeof(possiblekeys1)/sizeof(possiblekeys1[0]));



    /* set default options */

    // either sensible values or "not defined"=-1 for values that have to be specified in the input file and will be checked later.
    //ct->n_avg_ham = 1; // average hamilton over n_avg_ham steps to assimilate fast non-classical vibrations
    //ct->esp_scaling_factor = 1.; // scaling of the electrostatic potential of the environment
    //ct->opt_QMzone=0;
    //ct->neg_imag_pot = 0; //negative imaginary potential to drain the charge at some sites

/* evaluate input */

    /* read in stuff that is always needed */
    /* read seed for random number*/
    if (searchkey(lines1, input1, "seed", value, 0))
       ct->rnseed = atoi(value);

    /* Generation of hamiltonian */
    searchkey(lines1, input1, "hamiltonian", value, 1);
    if (strcmp(value, "dftb") == 0)
    {
        ct->hamiltonian_type = 0;
        PRINTF("Hamiltonian calculated by DFTB \n");
    }
    else if (strcmp(value, "ml") == 0)
    {
        ct->hamiltonian_type = 1;
        PRINTF("Hamiltonian predicted by Machine learning method\n");
    }

    if (searchkey(lines1, input1, "pureforces", value, 0))
    {
      ct->pure_forces = atoi(value);
      PRINTF("Evaluating pure forces only = %d\n", ct->pure_forces);
    }

    if (searchkey(lines1, input1, "forcedoutput", value, 0))
    {
      ct->force_output = atoi(value);
      PRINTF("Printing everything to file = %d\n", ct->force_output);
    }

    if (searchkey(lines1, input1, "neighborsonly", value, 0))
    {
      ct->neighbors_only = atoi(value);
      PRINTF("Calculating NL every %d steps\n", ct->neighbors_only);
    }


    if (ct->hamiltonian_type == 1)
    {

       if (searchkey(lines1, input1, "sitemodelpath", ml_site_energy_model_path, 0))
       {
          PRINTF("Training data for site energy file will be sought in %s\n", ml_site_energy_model_path);
       }
       else
       {
           PRINTF("No path for traning data for site energy ! \n");
           exit(-1);
       }

       if (searchkey(lines1, input1, "couplingmodelpath", ml_coupling_model_path, 0))
       {
          PRINTF("Training train for coupling file will be sought in %s\n", ml_coupling_model_path);
       }
       else
       {
           PRINTF("No path for traning data for couplings ! \n");
           exit(-1);
       }

       if (searchkey(lines1, input1, "siteparapath", ml_site_energy_para_path, 0))
       {
          PRINTF("Hyperparameters for site energy file will be sought in %s\n", ml_site_energy_para_path);
       }
       else
       {
           PRINTF("No hyperparameters for site energy ! \n");
           exit(-1);
       }

       if (searchkey(lines1, input1, "couplingparapath", ml_coupling_para_path, 0))
       {
          PRINTF("Hyperparameters for coupling file will be sought in %s\n", ml_coupling_para_path);
       }
       else
       {
           PRINTF("No hyperparameters for couplings ! \n");
           exit(-1);
       }

       // if (searchkey(lines1, input1, "sigmacouplings", value, 0))
       // {
       //    ct->sigma[1] = atof(value);
       //    PRINTF("Value of sigma for coupling model %f\n", ct->sigma[1]);
       // }
       // else
       // {
       //     PRINTF("No width value for couplings ! \n");
       //     exit(-1);
       // }
       //
       // if (searchkey(lines1, input1, "sigmasites", value, 0))
       // {
       //    ct->sigma[0] = atof(value);
       //    PRINTF("Value of sigma for site model %f\n", ct->sigma[0]);
       // }
       // else
       // {
       //     PRINTF("No width value for site energy ! \n");
       //     exit(-1);
       // }
       // if (searchkey(lines1, input1, "ntrainsites", value, 0))
       // {
       //    ct->ntrain[0] = atof(value);
       //    PRINTF("Number of training data for site model %s\n", value);
       // }
       // else
       // {
       //     PRINTF("Input number of training data for site model ! \n");
       //     exit(-1);
       // }
       // if (searchkey(lines1, input1, "ntraincouplings", value, 0))
       // {
       //    ct->ntrain[1] = atof(value);
       //    PRINTF("Number of training data for coupling model %s\n", value);
       // }
       // else
       // {
       //     PRINTF("Input number of training data for coupling model ! \n");
       //     exit(-1);
       // }
    }

    if (searchkey(lines1, input1, "atomindex", value, 0))
    {
       ct->define_orbital_sign = 1;

       split_string_into_int(value, 3, ct->atom_index_sign, "atomindex");

       PRINTF("Atomic indexes for fixing the initial sign of coupling: %d %d %d\n",
             ct->atom_index_sign[0], ct->atom_index_sign[1], ct->atom_index_sign[2]);
    }
    else
    {
       ct->define_orbital_sign = 0;
       PRINTF("Sign of coupling are randomly initialized !");
    }

    if (searchkey(lines1, input1, "frcoupling", value, 0))
    {
       ct->freeze_coupling  = atof(value)/HARTREE_TO_EV;

       PRINTF("Freeze the electronic coupling to be average coupling: %f eV \n",
             ct->freeze_coupling*HARTREE_TO_EV);
    }
    else
    {
       ct->freeze_coupling = 0.0;
    }

    searchkey(lines1, input1, "slkopath", slko_path, 1);
    PRINTF("SLKO files will be sought in %s\n", slko_path);

    searchkey(lines1, input1, "chargecarrier", value, 1);
    if (strcmp(value, "hole") == 0)
    {
        ct->is_hole_transfer = 1;
        PRINTF("perform hole transfer\n");
    }
    else if (strcmp(value, "electron") == 0)
    {
        ct->is_hole_transfer = 0;
        PRINTF("perform electron transfer\n");
    }
    else
    {
        PRINTF("chargecarrier value not known\n");
        exit(-1);
    }
    if (searchkey(lines1, input1, "offdiagscaling", value, 0))
    {
        if (strcmp(value, "yes") == 0)
        {
            PRINTF("  scaling of off-diagonal elements applied \n");
            switch (ct->is_hole_transfer)
            {
                case 1: ct->offdiag_scaling = OFFDIAG_FACTOR_HOLE; break; // 1.540
                case 0: ct->offdiag_scaling = OFFDIAG_FACTOR_ELEC; break; // 1.795
            }
        }
        else if (strcmp(value, "no") == 0)
        {
            PRINTF("standard hamilton matrix. off-diagonal elements decreased to 1.0 scale\n");
        }
        else
        {
            ct->offdiag_scaling = atof(value);
            PRINTF("standard hamilton matrix. off-diagonal elements  decreased to %4.2f scale\n", ct->offdiag_scaling);
        }
    }
    else
    {
        ct->offdiag_scaling = 1.0;
        PRINTF("standard hamilton matrix. off-diagonal elements unscaled\n");
    }

    if (searchkey(lines1, input1, "extchrmode", value, 0))
    {
        if (strcmp(value, "vacuo") == 0)
        {
            ct->qmmm = 0;
            PRINTF("\"in vacuo\" calculation - no QM/MM\n");
        }
        else if (strcmp(value, "qmmm") == 0)
        {
            ct->qmmm = 1;
            PRINTF("QM/MM calculation - the charge-transfer hamiltonian affected by the electric field\n");
        }
        else if (strcmp(value, "list") == 0) //formally QMG
        {
            ct->qmmm = 2;
            PRINTF("QM/MM calculation - the list of MM atoms will be read from file charge-transfer.ndx\n THIS HAS TO BE TESTED FIRST");
            exit(-1);
            /* read the file "charge-transfer.ndx" here! */
            ct_get_index(&mm_list_size, &mm_list);
            //for (i=0; i<mm_list_size; i++) PRINTF(" %d", mm_list[i]); PRINTF("\n");
        }
        else if (strcmp(value, "pme") == 0)
        {
            ct->qmmm = 3;
            PRINTF("QM/MM calculation with particle--mesh Ewald summation\n");
        }
        else
        {
            PRINTF("Didn't understand treatment of external charges (extcharmode).\n");
            exit(-1);
        }
    }
    else
    {
        ct->qmmm = 0;
        PRINTF("\"in vacuo\" calculation - no QM/MM\n");
    }
    if (searchkey(lines1, input1, "espscaling", value, 0))
    {
        ct->esp_scaling_factor = atof(value);
        PRINTF("the electrostatic interaction with MM atoms will be attenuated by a factor of %f\n", ct->esp_scaling_factor);
    }
    else
    {
        ct->esp_scaling_factor = 1.0; //default
    }
    if (searchkey(lines1, input1, "efield", value, 0))
    {
        split_string_into_double(value, 3, ct->efield, "efield");
        PRINTF("electric field applied: direction = %lf %lf %lf \n  magnitude[V/cm]: %lf \n",
               ct->efield[0]/dnorm(ct->efield), ct->efield[1]/dnorm(ct->efield), ct->efield[2]/dnorm(ct->efield), dnorm(ct->efield));
    }
    else
    {
        for (i = 0; i < 3; i++)
        {
            ct->efield[i] = 0.0;
        }
    }


    /* read in job-specific stuff */

    searchkey(lines1, input1, "jobtype", value, 1);

    ct->decoherence = 0;//default is propagator without decoherence
    if (strcmp(value, "scc") == 0)
    {
        ct->jobtype = cteSCCDYNAMIC;
        PRINTF("Fully coupled electron-ion dynamics\n");
    }
    else if (strcmp(value, "cpfscc") == 0)
    {
        ct->jobtype = cteCPFSCCDYNAMIC;
        PRINTF("Fully coupled electron-ion dynamics with coherence penalty functional\n");
    }
    else if (strcmp(value, "ddbscc") == 0)
    {
        ct->jobtype = cteDDBSCCDYNAMIC;
        PRINTF("Fully coupled electron-ion dynamics with decohernce detailed balance correction\n");
    }
    else if (strcmp(value, "adi") == 0)
    {
        ct->jobtype = cteADIABATIC;
        PRINTF("Adiabatic Born-Oppenheimer (SCF) dynamics of the hole\n");
    }
    else if (strcmp(value, "bod") == 0) //JJK
    {
        ct->jobtype = cteBORNOPPENHEIMER;
        PRINTF("Born-Oppenheimer dynamics with explicit following of the wave function\n");
    }
    else if (strcmp(value, "non") == 0)
    {
        ct->jobtype = cteNONSCCDYNAMIC;
        PRINTF("Uncoupled dynamics of the hole - w/o the polarization of solvent\n");
    }
    else if (strcmp(value, "par") == 0)
    {
        ct->jobtype = ctePARAMETERS;
        PRINTF("Calculation of charge-transfer parameters.\n");
    }
    else if (strcmp(value, "adn") == 0)
    {
        ct->jobtype = cteADNONSCC; // adiabatic non-self-consistent-charges
        PRINTF("Adiabatic Born-Oppenheimer (SCF) dynamics of the hole w/o the polarization of solvent\n");
    }
    else if (strcmp(value, "nom") == 0)
    {
        ct->jobtype = cteNOMOVEMENT;
        PRINTF("Stationary charge, calculation of all contributions to hamiltonian\n");
    }
    else if (strcmp(value, "fad") == 0)
    {
        ct->jobtype = cteFERMIADIABATIC;
        PRINTF("Dynamics of the adiabatic ground state obtained from the Fermi-distribution-based combination\n");
    }
    else if (strcmp(value, "dlzsh") == 0)
    {
        ct->jobtype = cteDLZSH; //lzsh is new one
        PRINTF("Landau-Zener type surface hopping method based on original diabatic LZ formula \n");
    }
    else if (strcmp(value, "alzsh") == 0)
    {
        ct->jobtype = cteALZSH; //lzsh is new one
        PRINTF("Landau-Zener type surface hopping method based on refomulated LZ formula \n");
    }
    else if (strcmp(value, "dfssh") == 0)
    {
        ct->jobtype = cteDFSSH;
        PRINTF("Diabatic fewest switches surface hopping method\n");
        ct->use_strong_dec=1;
    }
    else if (strcmp(value, "jfssh") == 0)
    {
        ct->jobtype = cteJFSSH;
        PRINTF("Adiabatic fewest switches surface hopping method by Johen\n");
        ct->use_strong_dec=1;
    }
    else if (strcmp(value, "bcjfssh") == 0)
    {
      ct->jobtype = cteBCJFSSH;
      PRINTF("Bolzmann corrected adiabatic fewest switches surface hopping method \n");
      ct->use_strong_dec=1;
    }
    else if (strcmp(value, "scrdfssh") == 0)
    {
        ct->jobtype = cteSCRDFSSH;
        PRINTF("Self-consistent restricted-decohenrence fewest switches surface hopping method\n");
	ct->use_strong_dec=1;
    }
    else if (strcmp(value, "gfsh") == 0)
    {
        ct->jobtype = cteGFSH;
        PRINTF("Global flux fewest switches surface hopping method\n");
	ct->use_strong_dec=1;
    }
    else if (strcmp(value, "ccfssh") == 0)
    {
        ct->jobtype = cteCCFSSH;
        PRINTF("Crossing-corrected fewest switches surface hopping method\n");
	ct->use_strong_dec=1;
    }
    else if (strcmp(value, "dish") == 0)
    {
        ct->jobtype = cteDISH;
        PRINTF("Decoherence induced fewest switches surface hopping method\n");
	ct->use_strong_dec=1;
    }
    else if (strcmp(value, "tfl") == 0)
    {
        ct->jobtype = cteTULLYLOC; //tfl is new one
        PRINTF("Tully's fewest switches surface hopping adapted for systems with localized, spatially spread-out adiab. states \n");
	ct->use_strong_dec=1; // ugly Command by Julian needs to be 1 for dholub 15.08.18
    }
    else if ((strcmp(value, "per") == 0 || strcmp(value, "ped") == 0))
    {
        ct->jobtype = ctePERSICOSFHOPPING;
        PRINTF("Persico's locally diabatic surface hopping between adiabatic states from Fermi-dist. combination\n");
        if (strcmp(value, "ped") == 0)
        {
            ct->decoherence = 1;
            PRINTF(" - correction for quantum decoherence switched on!\n");
        }
    }
    else if (strcmp(value, "ngl") == 0)
    {
        ct->jobtype = cteNEGFLORENTZ;
        PRINTF("Calculation of electric current with non-equlibrium Green's function approach + Lorentzian functions\n");
        PRINTF(" - populations of molecules mapped onto MD charges (self-consistent calculation)\n");
    }
    else if (strcmp(value, "ngn") == 0)
    {
        ct->jobtype = cteNEGFLORENTZNONSCC;
        PRINTF("Calculation of electric current with non-equlibrium Green's function approach + Lorentzian functions\n");
        PRINTF(" - no mapping of charges to MD (non-self-consistent calculation)\n");
    }
    else if (strcmp(value, "esp") == 0)
    {
        ct->jobtype = cteESP;
        PRINTF("Calculation of electrostatic potential only\n");
    }
    else if (strcmp(value, "tda") == 0)
    {
        ct->jobtype = cteTDA;
        PRINTF("Calculation of Tunneling matrix elements through bridge.\n");
    }
    else if ((strcmp(value, "gfs") == 0) || strcmp(value, "gfd") == 0)
    {
        ct->jobtype = ctePREZHDOSFHOPPING;
        printf("Global Flux Surface Hopping\n");
        if (strcmp(value, "gfd") == 0)
        {
            ct->decoherence = 1;
            PRINTF(" - correction for quantum decoherence switched on!\n");
        }
    }
    else
    {
        PRINTF("Did not understand jobtype.\n");
        exit(-1);
    }

    if (searchkey(lines1, input1, "nstqm", value, 0))
    {
        ct->interval = atoi(value);
        if (ct->jobtype == ctePARAMETERS || ct->jobtype == cteNOMOVEMENT || ct->jobtype == cteESP || ct->jobtype == cteTDA)
        {
            PRINTF("Performing QM calculation every %d steps.\n", ct->interval);
        }
        else if (ct->interval > 1)
        {
            PRINTF("WARNING: specified nstqm > 1, which will be ignored in the specified jobtype.\n");
        }
    }
    else
    {
        ct->interval = 1;
    }

    if (searchkey(lines1, input1, "nstaverage", value, 0))
    {
        ct->n_avg_ham = atoi(value);
        PRINTF("WARNING: averaging hamilton over %d steps to assimilate fast non-classical vibrations\n This is just for testing, so take care.\n", ct->n_avg_ham);
    }
    else
    {
        ct->n_avg_ham = 1; //default
    }

    if (searchkey(lines1, input1, "tfermi", value, 0))
    {
            ct->telec    = atof(value);
            ct->fermi_kt = BOLTZMANN_HARTREE_KELVIN * ct->telec;

            PRINTF("Temperature is %f K, kT = %f a.u.\n", ct->telec, ct->fermi_kt);
    }


    // search keys that are relevant when an actual charge is in the system //
    if (!(ct->jobtype == ctePARAMETERS || ct->jobtype == cteNEGFLORENTZ || ct->jobtype == cteNEGFLORENTZNONSCC ||
          ct->jobtype == cteESP || ct->jobtype == cteTDA))
    {

        if (searchkey(lines1, input1, "epol", value, 0))
        {
            if (strcmp(value, "imp") == 0)
            {
                ct->do_epol = 1;
                PRINTF("implicit electronic polarization applied (Born model)\n");
            }
            else
            {
                PRINTF("Did not understand electronic polarization model. Use IMP for born-like implicit electronic polarization.\n");
                exit(-1);
            }
            if (ct->jobtype == ctePARAMETERS || ct->jobtype == cteNEGFLORENTZ || ct->jobtype == cteNEGFLORENTZNONSCC ||
                ct->jobtype == cteESP || ct->jobtype == cteTDA)
            {
                PRINTF("WARNING: specified electronic polarization, which makes only sense for jobtypes with actual charge in the system.\n");
            }
        }
        else
        {
            ct->do_epol = 0; //default
        }

        if (searchkey(lines1, input1, "projection", value, 0))
        {
            if (strcmp(value, "yes") == 0)
            {
                ct->do_projection = 1;
                PRINTF("Projecting the charge carrier wavefunction at every step onto the new FO basis before propagating it.\n");
            }
            else if (strcmp(value, "no") == 0)
            {
                PRINTF("No projection applied\n");
            }
            else
            {
                PRINTF("Did not understand projection keyword.\n");
                exit(-1);
            }
            if (ct->jobtype == ctePARAMETERS || ct->jobtype == cteNEGFLORENTZ || ct->jobtype == cteNEGFLORENTZNONSCC ||
                ct->jobtype == cteESP || ct->jobtype == cteTDA)
            {
                PRINTF("WARNING: selected projection of the charge carrier wavefunction, which makes only sense if there is a charge carrier to project.\n");
            }
        }
        else
        {
            ct->do_projection = 0; //default
            PRINTF("No projection applied\n");
        }

        if (searchkey(lines1, input1, "sic", value, 0))
        {
            ct->sic = atof(value);
            PRINTF("Naive self-interaction correction applied, second-order term scaled by factor %f\n", ct->sic);
            if (ct->jobtype == ctePARAMETERS || ct->jobtype == cteNEGFLORENTZ || ct->jobtype == cteNEGFLORENTZNONSCC ||
                ct->jobtype == cteESP || ct->jobtype == cteTDA)
            {
                PRINTF("WARNING: specified self interaction corrction (SIC), which makes only sense for jobtypes with actual charge in the system.\n");
            }
        }
        else
        {
            ct->sic = 0.0; //default
            PRINTF("Omitting second-order terms.\n");
        }
        if (searchkey(lines1, input1, "internalrelax", value, 0))
        {
             ct->lambda = 0.0;

            if (strcmp(value, "no") == 0)
            {
                ct->do_lambda_i = 0;
                PRINTF("No inner-sphere reorganization energy\n");
            }
            else if (strcmp(value, "parameter") == 0) // former L_I
            {
                ct->do_lambda_i = 1;
                PRINTF("Emulation of inner-sphere reorganization energy with precalculated parameter.\n");
            }
            else if (strcmp(value, "onsite") == 0) // former LIQM
            {
                ct->do_lambda_i = 2;
                PRINTF("Emulation of internal relaxation by adding DFTB-QM forces to the force field\n");

            }
            else if (strcmp(value, "full") == 0) // former LQM
            {
                ct->do_lambda_i = 3;
                PRINTF("Emulation of inter- and intra-site relaxation by adding DFTB-QM forces to the force field\n");
            }
            else
            {
                ct->do_lambda_i = 0;
                ct->lambda = atof(value)/HARTREE_TO_EV;
                PRINTF("Emulation of inner-sphere reorganization energy with precalculated parameter of %f eV \n", atof(value));
//                PRINTF("Did not understand relaxation model.\n");
//                exit(-1);
            }
            if (ct->interval != 1)
            {
                PRINTF("Application of internal relaxation makes only sense if QM calculations are performed every step\n");
                exit(-1);
            }
            if (ct->jobtype == ctePARAMETERS || ct->jobtype == cteNEGFLORENTZ || ct->jobtype == cteNEGFLORENTZNONSCC ||
                ct->jobtype == cteESP || ct->jobtype == cteTDA)
            {
                PRINTF("WARNING: specified internal relaxation, which makes only sense for jobtypes with actual charge in the system.\n");
            }
        }
        else
        {
            ct->do_lambda_i = 0; //default
            PRINTF("No inner-sphere reorganization energy\n");
        }

        if (searchkey(lines1, input1, "deltaqmode", value, 0))
        {
            if (strcmp(value, "mulliken") == 0)
            {
                ct->delta_q_mode = 0;
                PRINTF("Representing the charge carrier with Mulliken charges.");
            }
            else if (strcmp(value, "resp") == 0)
            {
                ct->delta_q_mode = 1;
                PRINTF("Representing the charge carrier with RESP charges.");
            }
            else
            {
                PRINTF("UNKOWN OPTION FOR DELTAQMODE");
                exit(-1);
            }
        }
        else
        {
            ct->delta_q_mode = 0; //default
            PRINTF("Representing the charge carrier with Mulliken charges.");
        }
    }


    /* build QM system */

    /* read in site specifications*/
    searchkey(lines1, input1, "nsitetypes", value, 1);
    ct->sitetypes = atoi(value);
    snew(ct->sitetype, ct->sitetypes);
    PRINTF("There are %d different type(s) of sites\n", ct->sitetypes);

    searchkey(lines1, input1, "typefiles", value, 1);
    split_string_into_string(value, ct->sitetypes, sitespecdat, "typefiles");

    for (i = 0; i < ct->sitetypes; i++)
    {
        lines2 = read_file(sitespecdat[i], input2, possiblekeys2, sizeof(possiblekeys2)/sizeof(possiblekeys2[0]));
        searchkey(lines2, input2, "natoms", value, 1);
        ct->sitetype[i].atoms = atoi(value);
        snew(ct->sitetype[i].atom, ct->sitetype[i].atoms);
        snew(ct->sitetype[i].atomtype, ct->sitetype[i].atoms);
        searchkey(lines2, input2, "nelectrons", value, 1);
        ct->sitetype[i].nel = atoi(value);
        searchkey(lines2, input2, "nallorbitals", value, 1);
        ct->sitetype[i].norb = atoi(value);
        if (searchkey(lines2, input2, "radius", value, 0))
        {
            ct->sitetype[i].radius = atof(value);
        }
        else if (ct->do_epol == 1)
        {
            PRINTF("RADIUS OF SITE NEEDED FOR BORN POLARIZATION MODEL\n");
            exit(-1);
        }
        ct->sitetype[i].type = i;
        if (searchkey(lines2, input2, "nqmmmbonds", value, 0))
        {
            ct->sitetype[i].bonds = atoi(value);
            PRINTF("Sitetype %d found to have %d bond(s) between the QM and the MM region \n", i+1, ct->sitetype[i].bonds );

            if (ct->sitetype[i].bonds > 0)
            {
                snew(ct->sitetype[i].QMLA, ct->sitetype[i].bonds);
                snew(ct->sitetype[i].MMLA, ct->sitetype[i].bonds);
                snew(ct->sitetype[i].nochrs, ct->sitetype[i].bonds);
                snew(ct->sitetype[i].addchrs, ct->sitetype[i].bonds);
                snew(ct->sitetype[i].extracharge, ct->sitetype[i].bonds);
                snew(ct->sitetype[i].nochr, ct->sitetype[i].bonds);
                snew(ct->sitetype[i].addchr, ct->sitetype[i].bonds);

                searchkey(lines2, input2, "nignorechr", value, 1);
                split_string_into_int(value, ct->sitetype[i].bonds, ct->sitetype[i].nochrs,  "nignorechr");
                counter = 0;
                for (j = 0; j < ct->sitetype[i].bonds; j++)
                {
                    snew(ct->sitetype[i].nochr[j], ct->sitetype[i].nochrs[j]);
                    counter += ct->sitetype[i].nochrs[j];
                }
                if (counter > 0) //if there are any charges to ignore at all
                {
                    searchkey(lines2, input2, "nameignorechr", value, 1);
                    l = 0;
                    for (j = 0; j < ct->sitetype[i].bonds; j++)
                    {
                        for (k = 0; k < ct->sitetype[i].nochrs[j]; k++)
                        {
                            l++;
                        }
                    }
                    split_string_into_string(value, l, dummy,  "nameignorechr");
                    l = 0;
                    for (j = 0; j < ct->sitetype[i].bonds; j++)
                    {
                        for (k = 0; k < ct->sitetype[i].nochrs[j]; k++)
                        {
                            strcpy(ct->sitetype[i].nochr[j][k], dummy[l]);
                            PRINTF("Ignoring charges on %s\n", ct->sitetype[i].nochr[j][k]);
                            l++;
                        }
                    }
                }

                searchkey(lines2, input2, "naddchr", value, 1);
                split_string_into_int(value, ct->sitetype[i].bonds, ct->sitetype[i].addchrs,  "naddchr");
                counter = 0;
                for (j = 0; j < ct->sitetype[i].bonds; j++)
                {
                    snew(ct->sitetype[i].addchr[j], ct->sitetype[i].addchrs[j]);
                    counter += ct->sitetype[i].addchrs[j];
                }
                if (counter > 0) //if there are any charges to add at all
                {
                    searchkey(lines2, input2, "totaladdchr", value, 1);
                    split_string_into_double(value, ct->sitetype[i].bonds, ct->sitetype[i].extracharge, "totaladdchr");

                    searchkey(lines2, input2, "nameaddchr", value, 1);
                    l = 0;
                    for (j = 0; j < ct->sitetype[i].bonds; j++)
                    {
                        for (k = 0; k < ct->sitetype[i].addchrs[j]; k++)
                        {
                            l++;
                        }
                    }
                    split_string_into_string(value, l, dummy,  "nameaddchr");
                    l = 0;
                    for (j = 0; j < ct->sitetype[i].bonds; j++)
                    {
                        PRINTF("Distributing charge of %f over atoms:\n", ct->sitetype[i].extracharge[j]);
                        for (k = 0; k < ct->sitetype[i].addchrs[j]; k++)
                        {
                            strcpy(ct->sitetype[i].addchr[j][k], dummy[l]);
                            PRINTF(" %s \n", ct->sitetype[i].addchr[j][k]);
                            l++;
                        }
                    }
                }
            }
        }
        else
        {
            ct->sitetype[i].bonds = 0;
        }
        ct->sitetype[i].connections = 0; // not yet fully implemented. can maybe done without reading in
        //PRINTF("Sitetype %d found to have %d connections to other fragments \n", i+1, ct->sitetype[i].connections );
        snew(ct->sitetype[i].QMCA, ct->sitetype[i].connections);
        snew(ct->sitetype[i].QMCN, ct->sitetype[i].connections);
        for (j = 0; j < ct->sitetype[i].connections; j++)
        {
            snew(ct->sitetype[i].QMCN[j], 2);
        }

        searchkey(lines2, input2, "nfragorbs", value, 1);
        ct->sitetype[i].homos = atoi(value);
        PRINTF("Considering %d MOs per site.\n", ct->sitetype[i].homos);
        snew(ct->sitetype[i].homo, ct->sitetype[i].homos);
        snew(ct->sitetype[i].hubbard, ct->sitetype[i].homos);
        snew(ct->sitetype[i].lambda_i, ct->sitetype[i].homos);

        searchkey(lines2, input2, "fragorbs", value, 1);
        split_string_into_int(value, ct->sitetype[i].homos, ct->sitetype[i].homo, "fragorbs");
        if (ct->sic > 0.0)
        {
            searchkey(lines2, input2, "hubbard", value, 1);
            split_string_into_double(value, ct->sitetype[i].homos, ct->sitetype[i].hubbard, "hubbard");
        }
        if (ct->do_lambda_i == 1)
        {
            searchkey(lines2, input2, "lambda_i", value, 1);
            split_string_into_double(value, ct->sitetype[i].homos, ct->sitetype[i].lambda_i, "lambda_i");
        }
        if (ct->delta_q_mode == 1)
        {
            snew(ct->sitetype[i].delta_q, ct->sitetype[i].homos);
            snew(ct->sitetype[i].delta_q[0], ct->sitetype[i].homos* ct->sitetype[i].atoms);
            for (j = 1; j < ct->sitetype[i].homos; j++)
            {
                ct->sitetype[i].delta_q[j] = ct->sitetype[i].delta_q[0]+j*ct->sitetype[i].atoms;
            }

            searchkey(lines2, input2, "dqresp", value, 1);
            l = 0;
            for (j = 0; j < ct->sitetype[i].homos; j++)
            {
                for (k = 0; k < ct->sitetype[i].atoms; k++)
                {
                    l++;
                }
            }
            split_string_into_string(value, l, dummy, "dqresp");
            l = 0;
            for (j = 0; j < ct->sitetype[i].homos; j++)
            {
                for (k = 0; k < ct->sitetype[i].atoms; k++)
                {
                    ct->sitetype[i].delta_q[j][k] = atof(dummy[l]);
                    l++;
                }
            }

        }
        if (searchkey(lines2, input2, "customocc", value, 0))
        {
            ct->sitetype[i].do_custom_occ = 1;
            snew(ct->sitetype[i].custom_occ, ct->sitetype[i].norb);
            split_string_into_double(value, ct->sitetype[i].norb, ct->sitetype[i].custom_occ, "customocc");
            sum = 0.0;
            for (j = 0; j < ct->sitetype[i].norb; j++)
            {
                sum += ct->sitetype[i].custom_occ[j];
            }
            if ((int) sum != ct->sitetype[i].nel)
            {
                PRINTF("WRONG ELECTRON NUMBER IN CUSTOM OCCUPATION\n");
                exit(-1);
            }
        }
        else
        {
            ct->sitetype[i].do_custom_occ = 0;
        }
    }

    /* read in sites*/
    searchkey(lines1, input1, "nsites", value, 1);
    ct->pool_size = atoi(value);
    ct->sites     = ct->pool_size;
    if (searchkey(lines1, input1, "zonesize", value, 0))
    {
        ct->sites = atoi(value);
        if (searchkey(lines1, input1, "optimizezone", value, 0))
        {
            if (strcmp(value, "yes") == 0)
            {
                ct->opt_QMzone = 1;
                PRINTF("Will search pool for energetically best site to start. Overwriting start wavefunction\n");
            }
            else if (strcmp(value, "no") == 0)
            {
                ct->opt_QMzone = 0;
            }
            else
            {
                PRINTF("Did not understand optimizezone option.\n");
                exit(-1);
            }
            if (ct->pool_size <= ct->sites)
            {
                PRINTF("pool of possible sites has to be larger then the QM zone\n");
                exit(-1);
            }
        }
        else
        {
            ct->opt_QMzone = 0; //default
        }
        if (ct->sitetypes > 1 && (ct->pool_size != ct->sites))
        {
            printf("Adaptive QM zone works only if every site is the same.\n");
            exit(-1);
        }
    }

    if (ct->sites == ct->pool_size)
    {
        PRINTF("QM system consists of %d sites:\n", ct->sites);
    }
    else
    {
        PRINTF("%d out of a pool of %d sites will constitute the QM system\n Starting with: ", ct->sites, ct->pool_size);
    }
    snew(ct->site, ct->sites);
    snew(ct->indFO, ct->sites);
    snew(ct->pool_site, ct->pool_size);

    searchkey(lines1, input1, "sites", value, 1);
    split_string_into_string(value, ct->pool_size, dummy, "sites");
    for (i = 0; i < ct->pool_size; i++)
    {
        ct->pool_site[i].resnr = atoi(dummy[i]);
        //ct->pool_site[i].resnr--; //apparently residues are numbered in gromacs starting from 0 but written out as starting from 1  CHANGED IN GROMACS4.6
    }

    searchkey(lines1, input1, "sitetypes", value, 1);
    split_string_into_string(value, ct->pool_size, dummy, "sitetypes");
    for (i = 0; i < ct->pool_size; i++)
    {
        ct->pool_site[i].type = atoi(dummy[i]);
        ct->pool_site[i].type--; //for convinient numbering in charge-transfer.dat from 1 to #_diffrent_sites
    }

    searchkey(lines1, input1, "sitescc", value, 1);
    split_string_into_string(value, ct->pool_size, dummy, "sitescc");
    for (i = 0; i < ct->pool_size; i++)
    {
        ct->pool_site[i].do_scc = atoi(dummy[i]);
        if (ct->do_lambda_i > 1 && ct->pool_site[i].do_scc != 0)
        {
            PRINTF("QM-forces were designed for DFTB1 formalism but you want to use DFTB2\n");
           /*  exit(-1);  dholub 13.08.18 changed the exit and try benchmark of influence of DFTB2 and the derivaties of HOMO*/
        }
    }

    /* build sites according to sitetypes */
    for (i = 0; i < ct->pool_size; i++)
    {
        snew(ct->pool_site[i].nochrs, ct->sitetype[ct->pool_site[i].type].bonds);
        snew(ct->pool_site[i].addchrs, ct->sitetype[ct->pool_site[i].type].bonds);
        snew(ct->pool_site[i].nochr, ct->sitetype[ct->pool_site[i].type].bonds);
        snew(ct->pool_site[i].addchr, ct->sitetype[ct->pool_site[i].type].bonds);
        for (j = 0; j < ct->sitetype[ct->pool_site[i].type].bonds; j++)
        {
            snew(ct->pool_site[i].addchr[j], ct->sitetype[ct->pool_site[i].type].addchrs[j]);
            snew(ct->pool_site[i].nochr[j], ct->sitetype[ct->pool_site[i].type].nochrs[j]);
        }
        snew(ct->pool_site[i].homo, ct->sitetype[ct->pool_site[i].type].homos);
        snew(ct->pool_site[i].lambda_i, ct->sitetype[ct->pool_site[i].type].homos);
        snew(ct->pool_site[i].custom_occ, ct->sitetype[ct->pool_site[i].type].norb);

        l                       = ct->pool_site[i].resnr;              // resnr is parked in l. otherwise the resnr would get lost by copying sitetype to site.
        m                       = ct->pool_site[i].do_scc;             // do_SCC is parked in m. otherwise the resnr would get lost by copying sitetype to site.
        ct->pool_site[i]        = ct->sitetype[ct->pool_site[i].type]; //not sure if copying is that easy. it is
        ct->pool_site[i].resnr  = l;
        ct->pool_site[i].do_scc = m;

        /* stuff that is unique for each site */
        snew(ct->pool_site[i].delta_q, ct->sitetype[ct->pool_site[i].type].homos);
        snew(ct->pool_site[i].delta_q[0], ct->sitetype[ct->pool_site[i].type].homos* ct->sitetype[ct->pool_site[i].type].atoms);
        for (j = 1; j < ct->pool_site[i].homos; j++)
        {
            ct->pool_site[i].delta_q[j] = ct->pool_site[i].delta_q[0]+j*ct->pool_site[i].atoms;
        }
        if (ct->delta_q_mode == 1)
        {
            for (j = 0; j < ct->pool_site[i].homos; j++)
            {
                for (k = 0; k < ct->pool_site[i].atoms; k++)
                {
                    ct->pool_site[i].delta_q[j][k] = ct->sitetype[ct->pool_site[i].type].delta_q[j][k];
                }
            }
        }
        snew(ct->pool_site[i].overlap, ct->pool_site[i].homos);
        snew(ct->pool_site[i].overlap[0], SQR(ct->pool_site[i].homos));
        for (j = 0; j < ct->pool_site[i].homos; j++)
        {
            ct->pool_site[i].overlap[j] = ct->pool_site[i].overlap[0] + j * ct->pool_site[i].homos;
        }
        snew(ct->pool_site[i].overlap_ref, ct->pool_site[i].homos);
        snew(ct->pool_site[i].overlap_ref[0], SQR(ct->pool_site[i].homos));
        for (j = 0; j < ct->pool_site[i].homos; j++)
        {
            ct->pool_site[i].overlap_ref[j] = ct->pool_site[i].overlap_ref[0] + j * ct->pool_site[i].homos;
        }

        snew(ct->pool_site[i].atom, ct->sitetype[ct->pool_site[i].type].atoms);
        snew(ct->pool_site[i].atomtype, ct->sitetype[ct->pool_site[i].type].atoms);
        snew(ct->pool_site[i].QMLA, ct->sitetype[ct->pool_site[i].type].bonds);
        snew(ct->pool_site[i].MMLA, ct->sitetype[ct->pool_site[i].type].bonds);
        snew(ct->pool_site[i].modif_extcharge, ct->sitetype[ct->pool_site[i].type].bonds);
        for (j = 0; j < ct->sitetype[ct->pool_site[i].type].bonds; j++)
        {
            snew(ct->pool_site[i].modif_extcharge[j], ct->sitetype[ct->pool_site[i].type].addchrs[j]);
            for (k = 0; k <  ct->sitetype[ct->pool_site[i].type].addchrs[j]; k++)
            {
                ct->pool_site[i].modif_extcharge[j][k] = -1; // -1 equals undetermined
            }
        }
        snew(ct->pool_site[i].QMCA, ct->sitetype[ct->pool_site[i].type].connections);
        snew(ct->pool_site[i].QMCN, ct->sitetype[ct->pool_site[i].type].connections);
        for (j = 0; j < ct->sitetype[ct->pool_site[i].type].connections; j++)
        {
            snew(ct->pool_site[i].QMCN[j], 2);
        }
        snew(ct->pool_site[i].com, 3);
        //PRINTF("DATA %d %d %d %s %lf %d %d %d %lf %lf \n", ct->pool_site[i].atoms, ct->pool_site[i].bonds,ct->pool_site[i].nochrs[0],ct->pool_site[i].nochr[0][0],ct->pool_site[i].extracharge[0],ct->pool_site[i].addchrs[0],ct->pool_site[i].homos,ct->pool_site[i].homo[0],ct->pool_site[i].hubbard[0],ct->pool_site[i].lambda_i[0]);
    }

    /* get the number of extcharges */
    switch (ct->qmmm)
    {
        case 1:
        case 3:
            ct->extcharges_cplx = top_global->natoms;
            break;
        case 2:
            ct->extcharges_cplx = mm_list_size;
            break;
        default:
            ct->extcharges_cplx = 0;
    }
    for (i = 0; i < ct->pool_size; i++)
    {
        site = &(ct->pool_site[i]);
        switch (ct->qmmm)
        {
            case 1:
            case 3:
                site->extcharges = top_global->natoms - site->atoms;
                for (j = 0; j < site->bonds; j++)
                {
                    site->extcharges -= site->nochrs[j];
                }

                if (i < ct->sites)//complex has only ct->sites sites
                {
                    ct->extcharges_cplx -= site->atoms;
                    for (j = 0; j < site->bonds; j++)
                    {
                        ct->extcharges_cplx -= site->nochrs[j];
                    }
                }
                break;
            case 2:
                /* NUMBER OF EXTCHARGES HERE! */
                printf("check if adaptive QM zone works with extcharge list\n"); exit(-1);
                site->extcharges = mm_list_size;
                /* DO NOT SUBTRACT ANYTHING HERE, YET! */
                break;
            default:
                site->extcharges = 0;
        }
    }
    if (ct->qmmm > 0)
    {
        snew(ct->extcharge_cplx, top_global->natoms);
        for (i = 0; i < ct->pool_size; i++)
        {
            snew(ct->pool_site[i].extcharge, top_global->natoms);// we allocate array a little bit larger than needed (natoms instead of extcharges) and let remaining enrtries blank. This way we can build the intersection of these arrays in order to find the extcharges of the complex
        }
    }

    /* set the first ct->sites of the pool active */
    for (i = 0; i < ct->pool_size; i++)
    {
        ct->pool_site[i].active = i < ct->sites ? 1 : 0;
    }
    for (i = 0; i < ct->sites; i++)
    {
        ct->site[i] = ct->pool_site[i];
    }

    // end build QM system



    /* set arrays regarding the complex */

    ct->dim        = 0;
    ct->atoms_cplx = 0;
    counter        = 0;
    PRINTF("Site   Residue   MO   DO_SCC?\n");
    for (i = 0; i < ct->sites; i++)
    {
        ct->indFO[i]    = ct->dim;
        ct->dim        += ct->site[i].homos;
        ct->atoms_cplx += ct->site[i].atoms;
        for (k = 0; k < ct->site[i].bonds; k++)
        {
            counter += ct->site[i].addchrs[k];
        }
        for (k = 0; k < ct->site[i].homos; k++)
        {
            PRINTF("  %d       %d      %d    %s\n", i+1, ct->site[i].resnr, ct->site[i].homo[k],  (ct->site[i].do_scc == 0) ? "NO" : "YES" );
        }
    }
    snew(ct->atom_cplx, ct->atoms_cplx);
    snew(ct->atomtype_cplx, ct->atoms_cplx);
    ct->modif_extcharges_cplx = counter;
    snew(ct->modif_extcharge_cplx, ct->modif_extcharges_cplx);
    for (i = 0; i < ct->modif_extcharges_cplx; i++)
    {
        ct->modif_extcharge_cplx[i] = -1; // -1 equals undetermined
    }
    snew(ct->hamiltonian, ct->dim);
    snew(ct->hamiltonian[0], SQR(ct->dim));
    for (j = 1; j < ct->dim; j++)
    {
        ct->hamiltonian[j] = ct->hamiltonian[0] + j * ct->dim;
    }

// DFTB hamiltonian
    snew(ct->hamiltonian_dftb, ct->dim);
    snew(ct->hamiltonian_dftb[0], SQR(ct->dim));
    for (j = 1; j < ct->dim; j++)
    {
        ct->hamiltonian_dftb[j] = ct->hamiltonian_dftb[0] + j * ct->dim;
    }

    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            ct->hamiltonian[i][j] = 0.0;
            ct->hamiltonian_dftb[i][j] = 0.0;
        }
    }

// store the factors to force sign of coupling to be positive/negative
    snew(ct->positive_coupling_factor, ct->dim);
    snew(ct->positive_coupling_factor[0], SQR(ct->dim));
    for (j = 1; j < ct->dim; j++)
    {
        ct->positive_coupling_factor[j] = ct->positive_coupling_factor[0] + j * ct->dim;
    }

// machine learning
    if (ct->hamiltonian_type)
    {
// ML hamiltonian
       snew(ct->hamiltonian_ml, ct->dim);
       snew(ct->hamiltonian_ml[0], SQR(ct->dim));
       for (j = 1; j < ct->dim; j++)
       {
           ct->hamiltonian_ml[j] = ct->hamiltonian_ml[0] + j * ct->dim;
       }

    // Neural Network initialization:
       tf_model coupling_model_full;
       tf_model *coupling_model = &coupling_model_full;
       printf("initializing coupling model\n");
       initialize_nn_model(coupling_model, ml_coupling_model_path, ml_coupling_para_path);
       check_nn_status(coupling_model);
       ct->coupling_model = coupling_model;

       tf_model site_model_full;
       tf_model *site_model = &site_model_full;
       printf("initializing site energy model\n");
       initialize_nn_model(site_model, ml_site_energy_model_path, ml_site_energy_para_path);
       check_nn_status(site_model);
       ct->site_model = site_model;
    }

    snew(ct->fo_shift, ct->dim);

    //averaging out fast oscillation of the hamiltonian. H[dim][dim][time]
    snew(ct->hamiltonian_history, ct->dim);
    for (i = 0; i < ct->dim; i++)
    {
        snew(ct->hamiltonian_history[i], ct->dim);
        for (j = 0; j < ct->dim; j++)
        {
            snew(ct->hamiltonian_history[i][j], ct->n_avg_ham);
        }
    }

    /* remaining allocations */
    snew(ct->hamiltonian_mod, ct->dim);
    snew(ct->hamiltonian_adiab, SQR(ct->dim));
    snew(ct->ev_adiab, ct->dim);
    snew(ct->evec_adiab, SQR(ct->dim));
    snew(ct->work_adiab, 3*ct->dim);
    snew(ct->occupation, ct->dim);
    snew(ct->pop_adia, ct->dim);
    snew(ct->ev_spec, ct->dim);
    snew(ct->evec_spec, SQR(ct->dim));
    snew(ct->work_spec, 3*ct->dim);
    snew(ct->hubbard, ct->dim);
    snew(ct->hubbard[0], SQR(ct->dim));
    for (j = 1; j < ct->dim; j++)
    {
        ct->hubbard[j] = ct->hubbard[0] + j * ct->dim;
    }


    /* Runge-Kutta related */
    ct->rk_neq = 2 * ct->dim;
    snew(ct->wf, ct->rk_neq);
    snew(ct->wf_exc, ct->dim);
    snew(ct->dwf, ct->rk_neq);
    snew(ct->rk_ymax, ct->rk_neq);
    ct->rk_tol = 1.e-8;
    snew(ct->rk_thres, ct->rk_neq);
    for (i = 0; i < ct->rk_neq; i++)
    {
        ct->rk_thres[i] = ct->rk_tol;
    }
    ct->rk_lenwrk = 32 * ct->rk_neq;
    snew(ct->rk_work, ct->rk_lenwrk);


    /* shift of hamilton diagonal elements */
    if (searchkey(lines1, input1, "foshift", value, 0))
    {
        PRINTF("Applying shift to Hamiltonian:\n");
        split_string_into_double(value, ct->dim, ct->fo_shift, "foshift");
        for (i = 0; i < ct->dim; i++)
        {
            PRINTF("  %f eV\n", ct->fo_shift[i]*HARTREE_TO_EV);
        }
    }


    /* read the wavefunction */

    if (!(ct->jobtype == ctePARAMETERS || ct->jobtype == cteNEGFLORENTZ || ct->jobtype == cteNEGFLORENTZNONSCC ||
          ct->jobtype == cteESP || ct->jobtype == cteTDA))
    {
        counter = 0;
        if (searchkey(lines1, input1, "adiabstart", value, 0))
        {
            if (strcmp(value, "yes") == 0)
            {
                counter++;
                ct->adiabstart = 1;
                PRINTF("Starting from lowest adiabatic state.\n");
            }
            else if (strcmp(value, "no") == 0)
            {
                ct->adiabstart = 0;
            }
        }
        else
        {
            ct->adiabstart = 0;
        }
        if (searchkey(lines1, input1, "wavefunctionreal", value, 0))
        {
            counter++;
            ct->adiabstart = 0;
            PRINTF("Read the wavefunction:\n");
            split_string_into_double(value, ct->dim, ct->wf, "wavefunctionreal");

            if (searchkey(lines1, input1, "wavefunctionim", value, 0))
            {
                split_string_into_double(value, ct->dim, &ct->wf[ct->dim], "wavefunctionreal");
            }
            ct->survival = 0.0;
            for (i = 0; i < ct->dim; i++)
            {
                PRINTF(" Re_wf[%d] = %7.4f, Im_wf[%d] = %7.4f\n", i+1, ct->wf[i], i+1, ct->wf[i + ct->dim]);
                ct->occupation[i] = SQR(ct->wf[i]) + SQR(ct->wf[ct->dim + i]);
                ct->survival     += ct->occupation[i];
            }
            PRINTF("Sum of occupations = %7.5f\n", ct->survival);
            if (ct->survival < 0.99 || ct->survival > 1.01) //should be normalized, 0.01 tolerance
            {
                PRINTF("WARNING: no normalized starting wave function was specified.\n");
                exit(-1);
            }
        }
        if (counter > 1)
        {
            PRINTF("Providing a starting wavefunction doesn't make sense with option 'adiabstart'.\n");
            exit(-1);
        }
        if (counter < 1)
        {
            PRINTF("Provide either a starting wavefunction via keyword 'wavefunctionreal' (and optional 'wavefunctionim') or alternatively take the lowest adiabatic state with 'adiabstart'.\n");
            exit(-1);
        }
    }

    if (searchkey(lines1, input1, "nnegimpot", value, 0))
    {
        PRINTF("Charge taken out of the several FOs by way of negative imaginary potential,\n");
        PRINTF("NEVER TESTET IN THIS VERSION. CHECK SOURCE CODE BEFORE PROCEEDING!\n");
        exit(-1); //TODO I somehow tried to adapt this feature. However, it was initially considered for one Orbital per site and I'm not sure if it will work with several FOs per site.
        ct->neg_imag_pot = atoi(value);
        snew(ct->site_neg_imag_pot, ct->neg_imag_pot);
        snew(ct->site_annihilated_occupation, ct->neg_imag_pot);

        searchkey(lines1, input1, "negimpotrate", value, 1);
        ct->neg_imag_pot_rate_constant = atof(value);
        PRINTF("   this will be done for %d FOs, with 1/tau = %e au-1 = %f ps-1.\n",
               ct->neg_imag_pot, ct->neg_imag_pot_rate_constant, ct->neg_imag_pot_rate_constant * PS_TO_AU);

        searchkey(lines1, input1, "negimpotfos", value, 1);
        split_string_into_int(value, ct->neg_imag_pot, ct->site_neg_imag_pot, "negimpotfos");
        for (i = 0; i < ct->neg_imag_pot; i++)
        {
            ct->site_neg_imag_pot[i]--;
            ct->site_annihilated_occupation[i] = 0.0;
            PRINTF("   FO index no. %d for NIP occupation removal.\n", ct->site_neg_imag_pot[i]);
        }
        if (!(ct->jobtype == cteSCCDYNAMIC || ct->jobtype == cteCPFSCCDYNAMIC || ct->jobtype == cteDDBSCCDYNAMIC || ct->jobtype == cteNONSCCDYNAMIC ))
        {
            PRINTF("negative imaginary potential is only for SCC NON SFH impelemented.\n");
            exit(-1);
        }
    }
    else
    {
        ct->neg_imag_pot = 0; //default
    }

    /* all read in and allocated at this point*/



    /* set constant hubbard elements */
    counter = 0;
    for (i = 0; i < ct->sites; i++)
    {
        for (j = 0; j < ct->site[i].homos; j++)
        {
            for (k = 0; k < ct->site[i].homos - j; k++)
            {
                ct->hubbard[counter][counter+k] = ct->sic * (ct->site[i].hubbard[j] + ct->site[i].hubbard[k])*0.5; //variation of MO-energy if an other orbital on this site is getting charged.
                if (ct->do_epol == 1)
                {
                    ct->hubbard[counter][counter+k] -= 1.0/(2.0*ct->site[i].radius*NM_TO_BOHR)*(1.0-1.0/EPSILON_OP); //influence of electronic polarization
                }
                ct->hubbard[counter+k][counter] = ct->hubbard[counter][counter+k];
            }
            if (ct->do_lambda_i == 1) //what does that do?? dholub jun18
            {
                ct->hubbard[counter][counter] -= ct->site[i].lambda_i[j];
            }
            counter++;
        }
    }

    /* want to introduce loop with an option 9.9.2016:
        on -> old styl of assigning qm-zone
        off-> new style : shall read a new file:
     */
    /* assign atoms */

    PRINTF("Assigning QM atoms.\n");

    /* build pool of all QM connection atoms */
    QMCApoolcounter = 0;
    for (j = 0; j < atoms->nr; j++)
    {
        if (!strncmp((*(atoms->atomname[j])), "CQMC", 4) ||
            !strncmp((*(atoms->atomname[j])), "NQMC", 4) ||
            !strncmp((*(atoms->atomname[j])), "OQMC", 4) ||
            !strncmp((*(atoms->atomname[j])), "SQMC", 4) ||
            !strncmp((*(atoms->atomname[j])), "KQMC", 4) ||
            !strncmp((*(atoms->atomname[j])), "LQMC", 4) ||
            !strncmp((*(atoms->atomname[j])), "IQMC", 4) ||
            !strncmp((*(atoms->atomname[j])), "BQMC", 4) ||
            !strncmp((*(atoms->atomname[j])), "PQMC", 4) ||
            !strncmp((*(atoms->atomname[j])), "FQMC", 4) ||
            !strncmp((*(atoms->atomname[j])), "CQMT", 4) || //CQMT "branching" atoms
            !strncmp((*(atoms->atomname[j])), "NQMT", 4) )
        {
            QMCApool[QMCApoolcounter] = j; //capping connection atoms are in different residues
            QMCApoolcounter++;
        }
    }
    PRINTF("Total number of connection atoms in the system: %d \n", QMCApoolcounter);


    /////* Assign the atom numbers and types that the sites are composed of */////
    for (i = 0; i < ct->pool_size; i++)
    {
        site        = &(ct->pool_site[i]);
        counter     = 0;
        QMLAcounter = 0;
        MMLAcounter = 0;
        QMCAcounter = 0;
        //PRINTF("Site %d (Residue %d): \n", i+1, site->resnr);
        for (j = 0; j < atoms->nr; j++)
        {
            if (atoms->resinfo[atoms->atom[j].resind].nr == site->resnr) /* atom j is in residue that site i corresponds to */

            {/* find QM link atoms */
                if (!strncmp((*(atoms->atomname[j])), "CQML", 4) || /* QM link atom may be a C, N or O atom. (in most cases only C atoms are recommended) */
                    !strncmp((*(atoms->atomname[j])), "NQML", 4) ||
                    !strncmp((*(atoms->atomname[j])), "OQML", 4) ||
                    !strncmp((*(atoms->atomname[j])), "KQML", 4) ||
                    !strncmp((*(atoms->atomname[j])), "LQML", 4) ||
                    !strncmp((*(atoms->atomname[j])), "IQML", 4) ||
                    !strncmp((*(atoms->atomname[j])), "BQML", 4) ||
                    !strncmp((*(atoms->atomname[j])), "PQML", 4) ||
                    !strncmp((*(atoms->atomname[j])), "FQML", 4) ||
                    !strncmp((*(atoms->atomname[j])), "SQML", 4))
                {
                    site->QMLA[QMLAcounter] = counter;
                    site->atom[counter]     = j;
                    switch ((*(atoms->atomname[j]))[0])
                    {
                        case 'C': site->atomtype[counter] = 0; break;
                        case 'N': site->atomtype[counter] = 2; break;
                        case 'O': site->atomtype[counter] = 3; break;
                        case 'S': site->atomtype[counter] = 4; break;
                        case 'F': site->atomtype[counter] = 5; break;
                        case 'K': site->atomtype[counter] = 6; break;
                        case 'L': site->atomtype[counter] = 7; break;
                        case 'B': site->atomtype[counter] = 8; break;
                        case 'P': site->atomtype[counter] = 9; break;
                        case 'I': site->atomtype[counter] = 10; break;
                        default: PRINTF("Unknown atom type for atom %d (%s), exiting!\n", j, (*(atoms->atomname[j]))); exit(-1);
                    }
                    //PRINTF("%5d (%5s, type %d)\n", site->atom[counter], (*(atoms->atomname[site->atom[counter]])), site->atomtype[counter]+1);
                    counter++;
                    QMLAcounter++;
                }
                else if (!strncmp((*(atoms->atomname[j])), "CQMC", 4) ||
                         !strncmp((*(atoms->atomname[j])), "NQMC", 4) ||
                         !strncmp((*(atoms->atomname[j])), "OQMC", 4) ||
                         !strncmp((*(atoms->atomname[j])), "BQMC", 4) ||
                         !strncmp((*(atoms->atomname[j])), "FQMC", 4) ||
                         !strncmp((*(atoms->atomname[j])), "PQMC", 4) ||
                         !strncmp((*(atoms->atomname[j])), "IQMC", 4) ||
                         !strncmp((*(atoms->atomname[j])), "CQMT", 4) ||
                         !strncmp((*(atoms->atomname[j])), "SQMC", 4) ||
                         !strncmp((*(atoms->atomname[j])), "KQMC", 4) ||
                         !strncmp((*(atoms->atomname[j])), "LQMC", 4) ||
                         !strncmp((*(atoms->atomname[j])), "NQMT", 4) )


                {
                    ct->pool_site[i].QMCA[QMCAcounter] = counter;
                    ct->pool_site[i].atom[counter]     = j;
                    switch ((*(atoms->atomname[j]))[0])
                    {
                        case 'C': site->atomtype[counter] = 0; break;
                        case 'N': site->atomtype[counter] = 2; break;
                        case 'O': site->atomtype[counter] = 3; break;
                        case 'S': site->atomtype[counter] = 4; break;
                        case 'F': site->atomtype[counter] = 5; break;
                        case 'K': site->atomtype[counter] = 6; break;
                        case 'L': site->atomtype[counter] = 7; break;
                        case 'B': site->atomtype[counter] = 8; break;
                        case 'P': site->atomtype[counter] = 9; break;
                        case 'I': site->atomtype[counter] = 10; break;

                        default: PRINTF("Unknown atom type for atom %d (%s), exiting!\n", j, (*(atoms->atomname[j]))); exit(-1);
                    }
                    //PRINTF("%5d (%5s, type %d)\n", site->atom[counter], (*(atoms->atomname[site->atom[counter]])), site->atomtype[counter]+1);
                    counter++;
                    ///// find capping for connection atoms /////
                    bond_length_best = 0.5; // = 0.5nm
                    m                = -1;
                    n                = -1;
                    for (k = 0; k < QMCApoolcounter; k++)
                    {
                        for (l = 0; l < 3; l++)
                        {
                            X[l] = state->x[j][l];
                            Y[l] = state->x[QMCApool[k]][l];
                        }
                        dvec_sub(X, Y, bond);
                        bond_length = dnorm(bond);
                        if (atoms->resinfo[atoms->atom[QMCApool[k]].resind].nr != site->resnr && bond_length < bond_length_best) //find best capping in neighboring residue
                        {
                            m                = QMCApool[k];                                                                      //atomnumber of best QMCA (k) for QMCA (j) is saved in m
                            bond_length_best = bond_length;
                        }
                    }
                    if (m < 0)
                    {
                        printf("error: no connection atom for QMCA no. %d found in a radius of 5 Angstrom \n", j);
                        exit(-1);
                    }
                    else
                    {
                        printf("nearest connection atom %d <- %d\n", j, m);
                    }
                    //cap with best QMCA
                    site->QMCN[QMCAcounter][0] = m;
                    site->QMCN[QMCAcounter][1] = m; // if connection atom is no branching point both neighbors are the same
                    site->atom[counter]        = m;
                    site->atomtype[counter]    = 6; // 6 is pseudo-atom
                    counter++;
                    QMCAcounter++;
                    if (!strncmp((*(atoms->atomname[j])), "CQMT", 4) ||
                        !strncmp((*(atoms->atomname[j])), "NQMT", 4)) /* search also second nearest neighbor */
                    {
                        QMCAcounter--;                                //is the same connection
                        bond_length_best = 0.5;
                        for (k = 0; k < QMCApoolcounter; k++)
                        {
                            if (QMCApool[k] != m)
                            {
                                for (l = 0; l < 3; l++)
                                {
                                    X[l] = state->x[j][l];
                                    Y[l] = state->x[QMCApool[k]][l];
                                }
                                dvec_sub(X, Y, bond);
                                bond_length = dnorm(bond);
                                if (atoms->resinfo[atoms->atom[QMCApool[k]].resind].nr != site->resnr && bond_length < bond_length_best) //find best capping in neighboring residue
                                {
                                    n                = QMCApool[k];                                                                      //atomnumber of best QMCA (k) for QMCA (j) is saved in n
                                    bond_length_best = bond_length;
                                }
                            }
                        }
                        if (n < 0)
                        {
                            printf("error: no connection atom for QMCA no. %d found in a radius of 5 Angstrom \n", j);
                            exit(-1);
                        }
                        else
                        {
                            printf("nearest connection atom %d <- %d\n", j, n);
                        }

                        site->QMCN[QMCAcounter][1] = n;
                        site->atom[counter]        = n;
                        site->atomtype[counter]    = 6; // 6 is pseudo-atom
                        counter++;
                        QMCAcounter++;
                    }
                }

                /* find QM atoms */
                else if (!strncmp((*(atoms->atomname[j])), "CQM", 3) ||
                         !strncmp((*(atoms->atomname[j])), "HQM", 3) ||
                         !strncmp((*(atoms->atomname[j])), "NQM", 3) ||
                         !strncmp((*(atoms->atomname[j])), "OQM", 3) ||
                         !strncmp((*(atoms->atomname[j])), "SQM", 3) ||
                         !strncmp((*(atoms->atomname[j])), "KQM", 3) ||  // Br
                         !strncmp((*(atoms->atomname[j])), "LQM", 3) ||	 // Cl
                                                                       //!strncmp((*(atoms->atomname[j])), "YQM", 3) || // Y was special pseudo atom
                         !strncmp((*(atoms->atomname[j])), "BQM", 3) ||
                         !strncmp((*(atoms->atomname[j])), "FQM", 3) ||
                         !strncmp((*(atoms->atomname[j])), "PQM", 3) ||
                         !strncmp((*(atoms->atomname[j])), "IQM", 3) )


                {
                    site->atom[counter] = j;
                    switch ((*(atoms->atomname[j]))[0])
                    {
                        case 'C':
                            site->atomtype[counter] = 0;
                            break;
                        case 'H':
                            site->atomtype[counter] = 1;
                            break;
                        case 'N':
                            site->atomtype[counter] = 2;
                            break;
                        case 'O':
                            site->atomtype[counter] = 3;
                            break;
                        case 'S':
                            site->atomtype[counter] = 4;
                            break;
                        case 'F':
                            site->atomtype[counter] = 5;
                            break;
                        case 'K':
                            site->atomtype[counter] = 6;
                            break;
                        case 'L':
                            site->atomtype[counter] = 7;
                            break;
                        case 'B':
                            site->atomtype[counter] = 8;
                            break;
                        case 'P':
                            site->atomtype[counter] = 9;
                            break;
                        case 'I':
                            site->atomtype[counter] = 10;
                            break;

                        //case 'Y' : site->atomtype[counter] = 6; break;
                        default: PRINTF("Unknown atom type for atom %d (%s), exiting!\n", j, (*(atoms->atomname[j]))); exit(-1);
                    }
                    //PRINTF("%5d (%5s, type %d)\n", site->atom[counter], (*(atoms->atomname[site->atom[counter]])), site->atomtype[counter]+1);
                    counter++;
                    counter_cplx++;
                }

                /* find MM link atoms */
                else if (!strncmp((*(atoms->atomname[j])), "CMML", 4) || /* MM link atom may be a C, N or O atom. (in most cases only C atoms are recommended) */
                         !strncmp((*(atoms->atomname[j])), "NMML", 4) ||
                         !strncmp((*(atoms->atomname[j])), "OMML", 4) ||
                         !strncmp((*(atoms->atomname[j])), "KMML", 4) ||
                         !strncmp((*(atoms->atomname[j])), "LMML", 4) ||
                         !strncmp((*(atoms->atomname[j])), "BMML", 4) ||
                         !strncmp((*(atoms->atomname[j])), "FMML", 4) ||
                         !strncmp((*(atoms->atomname[j])), "PMML", 4) ||
                         !strncmp((*(atoms->atomname[j])), "IMML", 4) ||
                         !strncmp((*(atoms->atomname[j])), "SMML", 4) ||
                         !strncmp((*(atoms->atomname[j])), "JMML", 4))
                {
                    site->MMLA[MMLAcounter] = counter;
                    site->atom[counter]     = j;
                    site->atomtype[counter] = 1; /* the MM link atom will be substitute by a link hydrogen! */
                    //PRINTF("%5d (%5s, type %d)\n", site->atom[counter], (*(atoms->atomname[site->atom[counter]])), site->atomtype[counter]+1);
                    counter++;
                    counter_cplx++;
                    MMLAcounter++;
                }

                /* ERROR checks */
                if (counter > site->atoms)
                {
                    PRINTF("Site %d found to have %d atoms, which is more than the expected number of %d, exiting!\n", i, counter, site->atoms);
                    exit(-1);
                }
            }
        }
        //PRINTF("\n");
        if (QMLAcounter != MMLAcounter)
        {
            PRINTF("Site %d found to have %d QM link atom(s) but %d MM link atom(s). \n", i, QMLAcounter, MMLAcounter);
            exit(-1);
        }
    }//end atom selection

    /* get atoms of the complex */
    counter = 0;
    for (i = 0; i < ct->sites; i++)
    {
        for (j = 0; j < ct->site[i].atoms; j++)
        {
            ct->atom_cplx[counter]     = ct->site[i].atom[j];
            ct->atomtype_cplx[counter] = ct->site[i].atomtype[j];
            counter++;
        }
    }

    /* find QM/MM caping pairs */
    for (i = 0; i < ct->pool_size; i++)
    {
        for (j = 0; j < ct->pool_site[i].bonds; j++)
        {
            bond_length_best = 0.5; // = 0.5nm
            m                = -1;
            for (k = 0; k < ct->pool_site[i].bonds; k++)
            {
                for (l = 0; l < 3; l++)
                {
                    X[l] = state->x[ct->pool_site[i].atom[ct->pool_site[i].MMLA[k]]][l];
                    Y[l] = state->x[ct->pool_site[i].atom[ct->pool_site[i].QMLA[j]]][l];
                }
                dvec_sub(X, Y, bond);
                bond_length = dnorm(bond);
                if (bond_length < bond_length_best) //find best capping
                {
                    m                = k;           //index of best MMLA (k) for QMLA (j) is saved in m
                    bond_length_best = bond_length;
                }
            }
            if (m < 0)
            {
                PRINTF("error: no MMLA for QMLA no. %d found in a radius of 5 Angstrom \n", j);
                exit(-1);
            }
            //switch best MMLA (m) with MMLA j
            l = ct->pool_site[i].MMLA[j];
            ct->pool_site[i].MMLA[j] = ct->pool_site[i].MMLA[m];
            ct->pool_site[i].MMLA[m] = l;
        }
    }

    mass = 0.0;
    for (j = 0; j < ct->pool_site[0].atoms; j++)
    {
        mass += mdatoms->massT[ct->pool_site[0].atom[j]];
    }
    ct->adapt_inv_tot_mass = 1.0 / mass;



    /////* select the external charges */////
    PRINTF("Assigning MM atoms.\n");
    if (ct->qmmm > 0)
    {
        for (j = 0; j < top_global->natoms; j++)
        {
            ct->extcharge_cplx[j] = -1;
        }
        for (i = 0; i < ct->pool_size; i++)
        {
            for (j = 0; j < top_global->natoms; j++)
            {
                ct->pool_site[i].extcharge[j] = -1;
            }
        }

        snew(counter_array, ct->pool_size);
        for (j = 0; j < top_global->natoms; j++)
        {
            if (ct->qmmm == 1 || ct->qmmm == 3 || ct_atom_in_group(j, mm_list, mm_list_size)) /* either normal QM/MM or group-QM/MM and j is in the group */
            {
                for (i = 0; i < ct->pool_size; i++)
                {
                    site        = &(ct->pool_site[i]);
                    environment = 1;                                             /* default: all atoms are environment */
//PRINTF("j %d i %d atoms.atom[j].resnr %d ct->site[i].resnr %d \n" ,j, i, atoms.atom[j].resnr , ct->site[i].resnr);
                    if (atoms->resinfo[atoms->atom[j].resind].nr == site->resnr) /* in same residue -> further investigations */
                    {
                        for (k = 0; k < site->bonds; k++)
                        {
                            for (l = 0; l < site->addchrs[k]; l++)
                            {
                                if (!strcmp((*(atoms->atomname[j])), site->addchr[k][l])) /* charge will be modified to electro-neutralize */
                                {
                                    site->modif_extcharge[k][l] = counter_array[i];
                                }
                            }
                            for (l = 0; l < site->nochrs[k]; l++)
                            {
                                if (!strcmp((*(atoms->atomname[j])), site->nochr[k][l])) /* charge will be ignored */
                                {
                                    environment = 0;
                                }
                            }
                        }
                        for (k = 0; k < site->atoms; k++)
                        {
                            if (j == site->atom[k]) /* exclude QM atoms */
                            {
                                environment = 0;
                            }
                        }
                    }
                    else /* atoms in neighboring residues are ignored if one connection atom of site i lies in this residue */
                    {
                        for (k = 0; k < site->connections; k++)
                        {
                            if (atoms->resinfo[atoms->atom[j].resind].nr == atoms->resinfo[atoms->atom[site->QMCN[k][0]].resind].nr ||
                                atoms->resinfo[atoms->atom[j].resind].nr == atoms->resinfo[atoms->atom[site->QMCN[k][1]].resind].nr)
                            {
                                environment = 0;
                            }
                        }
                    }

                    if (environment)
                    {
                        site->extcharge[counter_array[i]] = j;
                        counter_array[i]++;
                    }
                }
            }
        }

        /* set extcharges of complex */
        /* and determine which should be modified */
        for (i = 0; i < top_global->natoms; i++)
        {
            ct->extcharge_cplx[i] = i;
        }
        for (i = 0; i < ct->sites; i++)
        {
            k = find_intersection(top_global->natoms, ct->extcharge_cplx, ct->site[i].extcharge, ct->extcharge_cplx); // this should successively reduce the charges in ct->extcharge_cplx
            for (j = k; j < top_global->natoms; j++)                                                                  // k elements are common in both arrays
            {
                ct->extcharge_cplx[j] = -1;
            }
        }
        counter = 0;
        for (l = 0; l < ct->extcharges_cplx; l++)
        {
            for (i = 0; i < ct->sites; i++)
            {
                for (j = 0; j < ct->site[i].bonds; j++)
                {
                    for (k = 0; k < ct->site[i].addchrs[j]; k++)
                    {
                        if (ct->extcharge_cplx[l] == ct->site[i].extcharge[ ct->site[i].modif_extcharge[j][k] ]) // if one of the extcharges of the complex is the same atom that was modified in the monomer calculation, then also modify it in the complex calculation.
                        {
                            ct->modif_extcharge_cplx[counter] = l;
                            counter++;
                        }
                    }
                }
            }
        }

        /* verify the number of ext. charges */
        PRINTF("Number of external charges:\n");
        for (i = 0; i < ct->sites; i++)
        {
            PRINTF("Site %2d: original group %d, (possibly) restricted to %d\n", i+1, ct->site[i].extcharges, counter_array[i]);
            ct->site[i].extcharges = counter_array[i];
            for (j = 0; j < ct->site[i].bonds; j++)
            {
                for (k = 0; k < ct->site[i].addchrs[j]; k++)
                {
                    if (ct->site[i].modif_extcharge[j][k] > -1)
                    {
                        PRINTF("          modified extcharge for bond no. %5d: atom %5d - %s (residue %d)\n", j+1, ct->site[i].extcharge[ct->site[i].modif_extcharge[j][k]]+1,
                               *(atoms->atomname[ct->site[i].extcharge[ct->site[i].modif_extcharge[j][k]]]), atoms->resinfo[atoms->atom[ct->site[i].extcharge[ct->site[i].modif_extcharge[j][k]]].resind].nr+1 );
                    }
                }
            }
        }
        PRINTF("Complex: original group %d, (possibly) restricted to %d\n", ct->extcharges_cplx, counter_cplx);
        //ct->extcharges_cplx = counter_cplx;
        //if (counter_modif_cplx != 2*ct->sites) {
        PRINTF("         number of atoms cutting QM/MM boundary = %d (there are %d sites)\n", counter_modif_cplx, ct->sites);
        //  exit(-1);
        //}
        for (j = 0; j < counter_modif_cplx; j++)
        {
            PRINTF("          modified extcharge no. %5d: atom %5d - %s (residue %d)\n", j+1, ct->extcharge_cplx[ct->modif_extcharge_cplx[j]]+1,
                   *(atoms->atomname[ct->extcharge_cplx[ct->modif_extcharge_cplx[j]]]), atoms->resinfo[atoms->atom[ct->extcharge_cplx[ct->modif_extcharge_cplx[j]]].resind].nr+1);
        }

    }//end QMMM>0


    ///// JOB SPECIFIC PREPARATIONS /////
    /* NEGF data */
    /* TODO: THIS WAS ONLY COMMENTED OUT BECAUSE I HAVE NO IDEA ABOUT THESE CALCULATIONS. MAYBE IMPLEMENT NEGF CALCUlATIONS ALSO IN THE CURRENT CODE
       if (ct->jobtype == cteNEGFLORENTZ || ct->jobtype == cteNEGFLORENTZNONSCC) {
       snew(ct->negf_arrays, 1);
       ct->negf_arrays->n[0] = ct->sites;
       fscanf(f, "%ld\n", ct->negf_arrays->n_lorentz);
       fscanf(f, "%lf %lf\n", ct->negf_arrays->e_f_left, ct->negf_arrays->e_f_right);
       fscanf(f, "%lf %ld %ld\n", ct->negf_arrays->temp, ct->negf_arrays->n_poles_l, ct->negf_arrays->n_poles_r);
       fscanf(f, "%lf %lf %lf\n", ct->negf_arrays->gam_l, ct->negf_arrays->eps_l, ct->negf_arrays->w0_l);
       fscanf(f, "%lf %lf %lf\n", ct->negf_arrays->gam_r, ct->negf_arrays->eps_r, ct->negf_arrays->w0_r);
       PRINTF("Parameters read for the non-eq. Green's functions:\n");
       PRINTF("  n = %ld, n_lorentz = %ld\n", ct->negf_arrays->n[0], ct->negf_arrays->n_lorentz[0]);
       PRINTF("  e_f_left = %lf, e_f_right = %lf\n", ct->negf_arrays->e_f_left[0], ct->negf_arrays->e_f_right[0]);
       PRINTF("  temp = %lf, n_poles_l = %ld, n_poles_r = %ld\n", ct->negf_arrays->temp[0], ct->negf_arrays->n_poles_l[0], ct->negf_arrays->n_poles_r[0]);
       PRINTF("  gam_l = %lf, eps_l = %lf, w0_l = %lf\n", ct->negf_arrays->gam_l[0], ct->negf_arrays->eps_l[0], ct->negf_arrays->w0_l[0]);
       PRINTF("  gam_r = %lf, eps_r = %lf, w0_r = %lf\n", ct->negf_arrays->gam_r[0], ct->negf_arrays->eps_r[0], ct->negf_arrays->w0_r[0]);
       }*/

    /* NEGF initialization including initial density matrix */
    if (ct->jobtype == cteNEGFLORENTZ || ct->jobtype == cteNEGFLORENTZNONSCC)
    {
        PRINTF("Initializing the NEGF calculation\n");
#ifdef GMX_MPI
        if (ct_mpi_rank == 0)
#endif
        negf_init_arrays(ct->negf_arrays, &(ct->rk_timestep), ct->wf);
    }

    /* DO HERE PREPARATIONS FOR BORNOPPENHEIMER! */
    if (ct->jobtype == cteBORNOPPENHEIMER)
    {
        snew(ct->born_overlap, ct->dim);
    }

    /* DO HERE PREPARATIONS FOR TFL! */ //imported by J.J. Kranz Frenkel Version 2018
    if (ct->jobtype == cteSCCDYNAMIC || ct->jobtype == cteCPFSCCDYNAMIC || ct->jobtype == cteDDBSCCDYNAMIC || ct->jobtype == cteDLZSH || ct->jobtype == cteALZSH || ct->jobtype == cteDFSSH || ct->jobtype == cteJFSSH || ct->jobtype == cteBCJFSSH || ct->jobtype == cteSCRDFSSH || ct->jobtype == cteCCFSSH || ct->jobtype == cteGFSH || ct->jobtype == cteDISH   || ct->jobtype == cteSCRDFSSH || ct->jobtype == cteTULLYLOC)
    {
        ct->surface = 0;

        snew(ct->tfs_popul, 2*ct->dim); /* complex array: Re(0), Re(1), ..., Re(n-1), Im(0), Im(1), ..., Im(n-1) */
        snew(ct->tfs_diab, 2*ct->dim);


        /* initial conditions - ground state occupied */
        ct->tfs_popul[0] = 1.; //changed if specific wf choosend as startign point

        for (i = 1; i < 2*ct->dim; i++)
        {
            ct->tfs_popul[i] = 0.;
        }
        snew(ct->tfs_popul_der, 2*ct->dim); /* complex array */
        snew(ct->tfs_vector, ct->dim);      /* tfs_vector[n]: n-th eigenvector of the CG Hamiltonian */
        snew(ct->tfs_vector[0], SQR(ct->dim));
        for (j = 1; j < ct->dim; j++)
        {
            ct->tfs_vector[j] = ct->tfs_vector[0] + j * ct->dim;
        }
        snew(ct->tfs_vector_old, ct->dim); /* tfs_vector_old[n]: tfs_vector[n] in the previous time step */
        snew(ct->tfs_vector_old[0], SQR(ct->dim));
        for (j = 1; j < ct->dim; j++)
        {
            ct->tfs_vector_old[j] = ct->tfs_vector_old[0] + j * ct->dim;
        }
        snew(ct->tfs_vector_old2, ct->dim); /* tfs_vector_old[n]: tfs_vector[n] in the previous time step */
        snew(ct->tfs_vector_old2[0], SQR(ct->dim));
        for (j = 1; j < ct->dim; j++)
        {
            ct->tfs_vector_old2[j] = ct->tfs_vector_old2[0] + j * ct->dim;
        }
        snew(ct->tfs_overlap, ct->dim); /* tfs_overlap[j][k]: <tfs_vector_old[j] | tfs_vector[k]> */
        snew(ct->tfs_overlap[0], SQR(ct->dim));
        for (j = 1; j < ct->dim; j++)
        {
            ct->tfs_overlap[j] = ct->tfs_overlap[0] + j * ct->dim;
        }
        snew(ct->tfs_overlap2, ct->dim); /* tfs_overlap[j][k]: <tfs_vector_old[j] | tfs_vector[k]> */
        snew(ct->tfs_overlap2[0], SQR(ct->dim));
        for (j = 1; j < ct->dim; j++)
        {
            ct->tfs_overlap2[j] = ct->tfs_overlap2[0] + j * ct->dim;
        }
        snew(ct->surf_prob, ct->dim);
        ct->tfs_initialization_step = 1;
        snew(ct->tfl_is_in_system, ct->dim);
        snew(ct->tfl_is_in_system_old, ct->dim);
        for (i = 0; i < ct->dim; i++)
        {
            ct->tfl_is_in_system[i] = 1; ct->tfl_is_in_system_old[i] = 1;
        }
        ct->tfl_num_of_states     = ct->dim;
        ct->tfl_num_of_states_old = ct->dim;
        snew(ct->tfl_mean_ham, SQR(ct->dim));

        snew(ct->tfl_mean_ham_full, ct->dim);
        snew(ct->tfl_mean_ham_full[0], SQR(ct->dim));
        for (j = 1; j < ct->dim; j++)
        {
            ct->tfl_mean_ham_full[j] = ct->tfl_mean_ham_full[0] + j * ct->dim;
        }

// add nonadiabatic coupling vector bewteen pre_surf and hop_surf, W.X, Nov.2018
        snew(ct->tfs_ncv, ct->sites);
        for (i = 0; i < ct->sites; i++)
        {
            snew(ct->tfs_ncv[i], ct->site[i].atoms);

            for (j = 0; j < ct->site[i].atoms; j++)
            {
                snew(ct->tfs_ncv[i][j], DIM);

            }
        }
        snew(ct->v1, ct->sites);
        for (i = 0; i < ct->sites; i++)
        {
            snew(ct->v1[i], ct->site[i].atoms);

            for (j = 0; j < ct->site[i].atoms; j++)
            {
                snew(ct->v1[i][j], DIM);

            }
        }
        snew(ct->v2, ct->sites);
        for (i = 0; i < ct->sites; i++)
        {
            snew(ct->v2[i], ct->site[i].atoms);

            for (j = 0; j < ct->site[i].atoms; j++)
            {
                snew(ct->v2[i][j], DIM);

            }
        }

//  adiabatic energies at previous time step and two time steps before
        snew(ct->ev_adiab_old2, ct->dim);
        snew(ct->ev_adiab_old, ct->dim);
        snew(ct->tfl_old_ham, SQR(ct->dim));
        snew(ct->tfl_old_ham_adiab, SQR(ct->dim));

//        srand48(123);
//        srand((unsigned)time(NULL));
        srand48(ct->rnseed);

        snew(ct->fo_overlap, ct->dim); /* tfs_overlap[j][k]: <tfs_vector_old[j] | tfs_vector[k]> */
        snew(ct->fo_overlap[0], SQR(ct->dim));
        for (j = 1; j < ct->dim; j++)
        {
            ct->fo_overlap[j] = ct->fo_overlap[0] + j * ct->dim;
        }

        snew(ct->fo_tdc, ct->dim); /* tfs_overlap[j][k]: <tfs_vector_old[j] | tfs_vector[k]> */
        snew(ct->fo_tdc[0], SQR(ct->dim));
        for (j = 1; j < ct->dim; j++)
        {
            ct->fo_tdc[j] = ct->fo_tdc[0] + j * ct->dim;
        }
    }

    /* DO HERE PREPARATIONS FOR PREZHDO! */
    if (ct->jobtype == ctePREZHDOSFHOPPING)
    {
        ct->surface = 0;
        snew(ct->tfs_popul, 2*ct->dim); /* complex array: Re(0), Re(1), ..., Re(n-1), Im(0), Im(1), ..., Im(n-1) */
        snew(ct->tfs_popul_old, 2*ct->dim);
        /* initial conditions - ground state occupied */
        ct->tfs_popul[0] = 1.;
        for (i = 1; i < 2*ct->dim; i++)
        {
            ct->tfs_popul[i] = 0.;
        }
        snew(ct->tfs_popul_der, 2*ct->dim); /* complex array */
        snew(ct->tfs_vector, ct->dim);      /* tfs_vector[n]: n-th eigenvector of the CG Hamiltonian */
        snew(ct->tfs_vector[0], SQR(ct->dim));
        for (j = 1; j < ct->dim; j++)
        {
            ct->tfs_vector[j] = ct->tfs_vector[0] + j * ct->dim;
        }
        snew(ct->tfs_vector_old, ct->dim); /* tfs_vector_old[n]: tfs_vector[n] in the previous time step */
        snew(ct->tfs_vector_old[0], SQR(ct->dim));
        for (j = 1; j < ct->dim; j++)
        {
            ct->tfs_vector_old[j] = ct->tfs_vector_old[0] + j * ct->dim;
        }
        snew(ct->tfs_overlap, ct->dim); /* tfs_overlap[j][k]: <tfs_vector_old[j] | tfs_vector[k]> */
        snew(ct->tfs_overlap[0], SQR(ct->dim));
        for (j = 1; j < ct->dim; j++)
        {
            ct->tfs_overlap[j] = ct->tfs_overlap[0] + j * ct->dim;
        }
        snew(ct->per_diab_hamiltonian, ct->dim); /* */
        snew(ct->per_diab_hamiltonian[0], SQR(ct->dim));
        for (j = 1; j < ct->dim; j++)
        {
            ct->per_diab_hamiltonian[j] = ct->per_diab_hamiltonian[0] + j * ct->dim;
        }
        snew(ct->surf_prob, ct->dim);
        snew(ct->ev_adiab_old, ct->dim);
        ct->tfs_initialization_step = 1;
    }

    /* DO HERE PREPARATIONS FOR PERSICO! */
    if (ct->jobtype == ctePERSICOSFHOPPING)
    {
        ct->surface = 0;
        snew(ct->tfs_popul, 2*ct->dim); /* complex array: Re(0), Re(1), ..., Re(n-1), Im(0), Im(1), ..., Im(n-1) */
        snew(ct->tfs_popul_old, 2*ct->dim);
        /* initial conditions - ground state occupied */
        ct->tfs_popul[0] = 1.;
        for (i = 1; i < 2*ct->dim; i++)
        {
            ct->tfs_popul[i] = 0.;
        }
        snew(ct->tfs_vector, ct->dim); /* tfs_vector[n]: n-th eigenvector of the CG Hamiltonian */
        snew(ct->tfs_vector[0], SQR(ct->dim));
        for (j = 1; j < ct->dim; j++)
        {
            ct->tfs_vector[j] = ct->tfs_vector[0] + j * ct->dim;
        }
        snew(ct->tfs_vector_old, ct->dim); /* tfs_vector_old[n]: tfs_vector[n] in the previous time step */
        snew(ct->tfs_vector_old[0], SQR(ct->dim));
        for (j = 1; j < ct->dim; j++)
        {
            ct->tfs_vector_old[j] = ct->tfs_vector_old[0] + j * ct->dim;
        }
        snew(ct->tfs_overlap, ct->dim); /* tfs_overlap[j][k]: <tfs_vector_old[j] | tfs_vector[k]> */
        snew(ct->tfs_overlap[0], SQR(ct->dim));
        for (j = 1; j < ct->dim; j++)
        {
            ct->tfs_overlap[j] = ct->tfs_overlap[0] + j * ct->dim;
        }
        snew(ct->per_propag_operator, ct->dim); /* */
        snew(ct->per_propag_operator[0], SQR(ct->dim));
        for (j = 1; j < ct->dim; j++)
        {
            ct->per_propag_operator[j] = ct->per_propag_operator[0] + j * ct->dim;
        }
        snew(ct->per_transformator, ct->dim); /* */
        snew(ct->per_transformator[0], SQR(ct->dim));
        for (j = 1; j < ct->dim; j++)
        {
            ct->per_transformator[j] = ct->per_transformator[0] + j * ct->dim;
        }
        snew(ct->per_diab_hamiltonian, ct->dim); /* */
        snew(ct->per_diab_hamiltonian[0], SQR(ct->dim));
        for (j = 1; j < ct->dim; j++)
        {
            ct->per_diab_hamiltonian[j] = ct->per_diab_hamiltonian[0] + j * ct->dim;
        }
        snew(ct->surf_prob, ct->dim);
        snew(ct->ev_adiab_old, ct->dim);
        ct->tfs_initialization_step = 1;
        /* auxiliary arrays for the orthogonalizer */
        snew(ct->per_arrays, 1);
        snew(ct->per_arrays->in, ct->dim*ct->dim);
        snew(ct->per_arrays->evec, ct->dim*ct->dim);
        ct->per_arrays->lwork = 26*ct->dim;
        snew(ct->per_arrays->work, 26*ct->dim*26*ct->dim);
        ct->per_arrays->liwork = 10*ct->dim;
        snew(ct->per_arrays->iwork, 10*ct->dim*10*ct->dim);
        snew(ct->per_arrays->eval, ct->dim);
        snew(ct->per_arrays->issupz, 2*ct->dim);
    }

#ifdef GMX_MPI
    printf("Completed charge transfer initialization at rank %d\n", ct_mpi_rank);
#else
    printf("Completed charge transfer initialization\n");
#endif


/*
   // Print out (nearly) all stuff that was read in
   for (k=0; k<ct->pool_size; k++){
    s=ct->pool_site[k];
    PRINTF("%d %d %d %d %d %d %d %d %d %f %d %d\n", s.type, s.resnr, s.atoms, s.bonds, s.connections, s.homos,s.extcharges, s.nel, s.norb, s.radius, s.do_scc, s.active);
    for (i=0; i<s.atoms; i++){
      PRINTF("%d %d \n", s.atom[i], s.atomtype[i]);
    }
    for (i=0; i<s.bonds; i++){
      PRINTF("%d %d %d \n", s.QMLA[i], s.MMLA[i], s.nochrs[i]);
      for (j=0; j<s.nochrs[i]; j++)
        PRINTF("%s\n", s.nochr[i][j]);
      for (j=0; j<s.addchrs[i]; j++)
        PRINTF("%s %d\n", s.addchr[i][j], s.modif_extcharge[i][j]);
    }
    for (i=0; i<s.extcharges; i++)
      PRINTF("%d ", s.extcharge[i]);
    for (i=0; i<s.homos; i++){
      PRINTF("%d %f %f\n", s.homo[i], s.hubbard[i], s.lambda_i[i]);
    }
   }
   PRINTF("%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d \n", ct->jobtype, ct->interval, ct->qmmm, ct->sitetypes, ct->sites, ct->dim, ct->is_hole_transfer, ct->atoms_cplx, ct->extcharges_cplx, ct->modif_extcharges_cplx, ct->n_avg_ham, ct->do_lambda_i, ct->do_epol, ct->decoherence, ct->pool_size, ct->opt_QMzone);
   PRINTF("%f %f %f %f %f    %f %f %f\n", ct->offdiag_scaling, ct->sic, ct->esp_scaling_factor, ct->telec, ct->fermi_kt, ct->efield[0],ct->efield[1],ct->efield[2]);
   for (i=0; i<ct->atoms_cplx; i++)
    PRINTF("%d %d \n", ct->atom_cplx[i], ct->atomtype_cplx[i]);
   for (i=0; i<ct->extcharges_cplx; i++)
    PRINTF("%d ", ct->extcharge_cplx[i]);
   for (i=0; i<ct->modif_extcharges_cplx; i++)
    PRINTF("%d ", ct->modif_extcharge_cplx[i]);
   for (i=0; i< ct->dim; i++)
    PRINTF("%f %f \n", ct->wf[i], ct->wf[i+ct->dim]);

 */

    return;
}

/****************************
 * INITIALIZE THE DFTB CODE *
 ****************************/

#ifdef GMX_MPI
void init_dftb(t_mdatoms *mdatoms, dftb_t *dftb, charge_transfer_t *ct, char *slko_path, int ct_mpi_rank)
#else
void init_dftb(t_mdatoms *mdatoms, dftb_t *dftb, charge_transfer_t *ct, char *slko_path)
#endif
{
/* Initializes the dftb calculations. Reads information from Slater-Koster files, copies information that was gathered in init_charge_transfer(), allocates memory */

// PARAMETERS:
// mdatoms     = (in) variable of GROMACS data-structure-type (we need atomic masses)
// dftb        = (out) main data structure for DFTB calculations
// ct          = (in) main data structure with information about the charge transfer calculation
// slko_path   = (in) path to the Slater-Koster files
// ct_mpi_rank = (in) the rank of the process in MPI calculations (each process in MPI calculations can be destinguished by its value of ct_mpi_rank)
//////////////////////////////////////////////////////////////////////////////////////////
    const char *type_symbols = "chnosfklbpi";        // fe is z 16.06.2016 because of one letter relation between element and code
                                                  // K is Br, L is Cl 13Sep2018 Daniel for halogens
    const char *suffix1      = "-c.spl";
    const char *suffix2      = "-uncomp-c.spl"; // now diffrent folder for each parameter set with diefferent r_dens and r_wf

    char        filename[128];
    int         i, j, k, l, izpj, counter;
    double      espin, qzeroh[3], uhubbh[3], mass;
    FILE       *f;

    char       *line = NULL;
    char      * pch;
    size_t      len;

#ifdef GMX_MPI
    printf("DFTB initialization at rank %d\n", ct_mpi_rank);
#else
    printf("DFTB initialization\n");
#endif


    for (i = 0; i < DFTB_MAXTYPES; i++)
    {
        dftb->qzero1[i] = 0.0;
        dftb->uhubb1[i] = 0.0;
        dftb->qzero2[i] = 0.0;
        dftb->uhubb2[i] = 0.0;
        for (j = 0; j < DFTB_MAXTYPES; j++)
        {
            PRINTF("Atomtype pair %d-%d\n", i+1, j+1);
            /* read the tables for DFTB phase 1 - calculation of monomers */
            sprintf(filename, "%s%c%c%s", slko_path, type_symbols[i], type_symbols[j], suffix1);
            f = fopen(filename, "r");
            if (f == NULL)
            {
                PRINTF("Cannot open the parameter file %s, exiting!\n", filename);
                exit(-1);
            }
            //printf("fscanf(f, \"%%lf %%d\", &(dftb->dr1[i][j]), &(dftb->dim1[i][j]))\n");
            //printf("%lf %d\n", dftb->dr1[i][j], dftb->dim1[i][j]);
            if (i == j)
            {
                fscanf(f, "%lf %d %d", &(dftb->dr1[i][j]), &(dftb->dim1[i][j]), &(dftb->lmax[i]));
                //printf("&(dftb->skself1[i][0]), &(dftb->skself1[i][1]), &(dftb->skself1[i][2]), &espin\n");
                fscanf(f, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                       &(dftb->skself1[i][0]), &(dftb->skself1[i][1]), &(dftb->skself1[i][2]), &espin,
                       uhubbh+2, uhubbh+1, uhubbh, qzeroh+2, qzeroh+1, qzeroh);
                //printf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", dftb->skself1[i][0], dftb->skself1[i][1], dftb->skself1[i][2], espin,
                //          uhubbh[2], uhubbh[1], uhubbh[0], qzeroh[2], qzeroh[1], qzeroh[0]);
                dftb->uhubb1[i] = uhubbh[0];
                for (k = 0; k < 3; k++)
                {
                    dftb->qzero1[i] += qzeroh[k];
                }
            }
            else
            {
                fscanf(f, "%lf %d", &(dftb->dr1[i][j]), &(dftb->dim1[i][j]));
            }

            snew(dftb->skhtab1[i][j], dftb->dim1[i][j]);
            snew(dftb->skstab1[i][j], dftb->dim1[i][j]);
            for (k = 0; k < dftb->dim1[i][j]; k++)
            {
                //   printf("for (l=0; l<10; l++) fscanf(f, \"%%lf\", &(dftb->skhtab1[i][j][k][l]))\n");
                for (l = 0; l < 10; l++)
                {
                    fscanf(f, "%lf", &(dftb->skhtab1[i][j][k][l]));
                    // printf("%lf ", dftb->skhtab1[i][j][k][l]);
                }
                //printf("for (l=0; l<10; l++) fscanf(f, \"%%lf\", &(dftb->skstab1[i][j][k][l]))\n");
                for (l = 0; l < 10; l++)
                {
                    fscanf(f, "%lf", &(dftb->skstab1[i][j][k][l]));
                    //printf("%lf ", dftb->skstab1[i][j][k][l]);
                }
            }
            PRINTF("skfile for pair %d-%d, phase 1: %s\n", i+1, j+1, filename);
            fclose(f);

            /* read the tables for DFTB phase 2 - calculation of complex */
            sprintf(filename, "%s%c%c%s", slko_path, type_symbols[i], type_symbols[j], suffix2);
            f = fopen(filename, "r");
            if (f == NULL)
            {
                PRINTF("Cannot open the parameter file %s, exiting!\n", filename);
                exit(-1);
            }

            //in converntional SLKOs there are 3 entries for i=j, in our CT-SLKO lmax is missing. In order to be able to use both formats in phase2, fscan was replaced by getline
            //fscanf(f, "%lf %d", &(dftb->dr2[i][j]), &(dftb->dim2[i][j]));
            getline(&line, &len, f);
            pch              = strtok (line, " ");
            dftb->dr2[i][j]  = atof(pch);
            pch              = strtok (NULL, " ");
            dftb->dim2[i][j] = atoi(pch);

            //printf("%lf %d\n", dftb->dr2[i][j], dftb->dim2[i][j]);
            if (i == j)
            {
                fscanf(f, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                       &(dftb->skself2[i][0]), &(dftb->skself2[i][1]), &(dftb->skself2[i][2]), &espin,
                       uhubbh+2, uhubbh+1, uhubbh, qzeroh+2, qzeroh+1, qzeroh);
                //printf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", dftb->skself2[i][0], dftb->skself2[i][1], dftb->skself2[i][2], espin,
                //          uhubbh[2], uhubbh[1], uhubbh[0], qzeroh[2], qzeroh[1], qzeroh[0]);
                dftb->uhubb2[i] = uhubbh[0];
                for (k = 0; k < 3; k++)
                {
                    dftb->qzero2[i] += qzeroh[k];
                }
            }
            snew(dftb->skhtab2[i][j], dftb->dim2[i][j]);
            snew(dftb->skstab2[i][j], dftb->dim2[i][j]);
            for (k = 0; k < dftb->dim2[i][j]; k++)
            {
                for (l = 0; l < 10; l++)
                {
                    fscanf(f, "%lf", &(dftb->skhtab2[i][j][k][l]));
                    //printf("%lf ", dftb->skhtab2[i][j][k][l]);
                }
                for (l = 0; l < 10; l++)
                {
                    fscanf(f, "%lf", &(dftb->skstab2[i][j][k][l]));
                    //printf("%lf ", dftb->skstab2[i][j][k][l]);
                }
                //printf("\n");
            }
            PRINTF("skfile for pair %d-%d, phase 2: %s\n", i+1, j+1, filename);
            fclose(f);
        }
    }

    /* deal with phase1 and certain parts of phase2 */
    snew(dftb->phase1, ct->sites);
    dftb->phase2.nn   = 0;
    dftb->phase2.ne   = ct->extcharges_cplx;
    dftb->phase2.nel  = 0;
    dftb->phase2.norb = 0;

    /* prepare the broyden structures */
    snew(dftb->broyden, ct->sites);
    for (i = 0; i < ct->sites; i++)
    {
        snew(dftb->broyden[i].f, ct->site[i].atoms);
        snew(dftb->broyden[i].ui, ct->site[i].atoms);
        snew(dftb->broyden[i].vti, ct->site[i].atoms);
        snew(dftb->broyden[i].t1, ct->site[i].atoms);
        snew(dftb->broyden[i].dumvi, ct->site[i].atoms);
        snew(dftb->broyden[i].df, ct->site[i].atoms);
        snew(dftb->broyden[i].vector, ct->site[i].atoms);
        snew(dftb->broyden[i].unit31, ct->site[i].atoms);
        snew(dftb->broyden[i].unit32, ct->site[i].atoms);
    }
    /* determine the numbers of electrons and atom types (in arrays) */
    dftb->phase2.nel = 0;
    for (i = 0; i < ct->sites; i++)
    {
        dftb->phase1[i].nel  = ct->site[i].nel;
        dftb->phase1[i].norb = ct->site[i].norb;
        dftb->phase1[i].nn   = ct->site[i].atoms;
        dftb->phase2.nn     += ct->site[i].atoms;
        dftb->phase2.nel    += dftb->phase1[i].nel;
        dftb->phase2.norb   += dftb->phase1[i].norb;

        /* double arrays */
        snew(dftb->phase1[i].x, ct->site[i].atoms);
        snew(dftb->phase1[i].x_opt, ct->site[i].atoms);
        snew(dftb->phase1[i].xe, ct->site[i].extcharges);
        snew(dftb->phase1[i].mass, ct->site[i].atoms);
        snew(dftb->phase1[i].ze, ct->site[i].extcharges);
        snew(dftb->phase1[i].qmat, ct->site[i].atoms);
        snew(dftb->phase1[i].qmold, ct->site[i].atoms);
        snew(dftb->phase1[i].qmulli, dftb->phase1[i].norb);
        snew(dftb->phase1[i].ev, dftb->phase1[i].norb);
        snew(dftb->phase1[i].occ, dftb->phase1[i].norb);
        snew(dftb->phase1[i].a, dftb->phase1[i].norb);
        snew(dftb->phase1[i].a[0], SQR(dftb->phase1[i].norb));
        for (j = 1; j < dftb->phase1[i].norb; j++)
        {
            dftb->phase1[i].a[j] = dftb->phase1[i].a[0] + j * dftb->phase1[i].norb;
        }

        snew(dftb->phase1[i].a_old, dftb->phase1[i].norb);
        snew(dftb->phase1[i].a_old[0], SQR(dftb->phase1[i].norb));
        for (j = 1; j < dftb->phase1[i].norb; j++)
        {
            dftb->phase1[i].a_old[j] = dftb->phase1[i].a_old[0] + j * dftb->phase1[i].norb;
        }
        snew(dftb->phase1[i].a_ref, dftb->phase1[i].norb);
        snew(dftb->phase1[i].a_ref[0], SQR(dftb->phase1[i].norb));
        for (j = 1; j < dftb->phase1[i].norb; j++)
        {
            dftb->phase1[i].a_ref[j] = dftb->phase1[i].a_ref[0] + j * dftb->phase1[i].norb;
        }

        snew(dftb->phase1[i].b, dftb->phase1[i].norb);
        snew(dftb->phase1[i].b[0], SQR(dftb->phase1[i].norb));
        for (j = 1; j < dftb->phase1[i].norb; j++)
        {
            dftb->phase1[i].b[j] = dftb->phase1[i].b[0] + j * dftb->phase1[i].norb;
        }
        snew(dftb->phase1[i].a_trans, dftb->phase1[i].norb * dftb->phase1[i].norb);
        snew(dftb->phase1[i].b_trans, dftb->phase1[i].norb * dftb->phase1[i].norb);
        snew(dftb->phase1[i].hamil, dftb->phase1[i].norb);
        snew(dftb->phase1[i].hamil[0], SQR(dftb->phase1[i].norb));
        for (j = 1; j < dftb->phase1[i].norb; j++)
        {
            dftb->phase1[i].hamil[j] = dftb->phase1[i].hamil[0] + j * dftb->phase1[i].norb;
        }
        snew(dftb->phase1[i].overl, dftb->phase1[i].norb);
        snew(dftb->phase1[i].overl[0], SQR(dftb->phase1[i].norb));
        for (j = 1; j < dftb->phase1[i].norb; j++)
        {
            dftb->phase1[i].overl[j] = dftb->phase1[i].overl[0] + j * dftb->phase1[i].norb;
        }

        //snew(dftb->phase1[i].overl_diffuse, dftb->phase1[i].norb);
        //  snew(dftb->phase1[i].overl_diffuse[0], SQR(dftb->phase1[i].norb));
        //  for(j = 1; j < dftb->phase1[i].norb; j++)
        //    dftb->phase1[i].overl_diffuse[j] = dftb->phase1[i].overl_diffuse[0] + j * dftb->phase1[i].norb;
        //dftb->phase1[i].gammamat = (double **) malloc(ct->atoms[i] * sizeof(double *));
        //  dftb->phase1[i].gammamat[0] = (double *) malloc(ct->atoms[i] * ct->atoms[i] * sizeof(double));
        snew(dftb->phase1[i].gammamat, ct->site[i].atoms);
        snew(dftb->phase1[i].gammamat[0], SQR(ct->site[i].atoms));
        for (j = 1; j < ct->site[i].atoms; j++)
        {
            dftb->phase1[i].gammamat[j] = dftb->phase1[i].gammamat[0] + j * ct->site[i].atoms;
        }
        snew(dftb->phase1[i].shift, ct->site[i].atoms);
        snew(dftb->phase1[i].dshift, ct->site[i].atoms);
        snew(dftb->phase1[i].shiftE, ct->site[i].atoms);
        snew(dftb->phase1[i].aux, 3 * dftb->phase1[i].norb);

        if (ct->do_lambda_i > 1 || ct->jobtype == cteJFSSH || ct->jobtype == cteBCJFSSH || ct->jobtype == cteNOMOVEMENT || ct->hamiltonian_type)
        {
            snew(dftb->phase1[i].grad, dftb->phase1[i].nn);
            snew(dftb->phase1[i].partgrad, dftb->phase1[i].nn);
        }
        if (ct->qmmm == 3)
        {
            snew(dftb->phase1[i].neighbors_pme, ct->site[i].atoms);
            snew(dftb->phase1[i].neighbor_pme, ct->site[i].atoms);
        }

        /* int arrays */
        snew(dftb->phase1[i].izp, ct->site[i].atoms);
        snew(dftb->phase1[i].ind, ct->site[i].atoms + 1);
        dftb->phase1[i].ind[0] = 0;
        for (j = 0; j < ct->site[i].atoms; j++)
        {
            dftb->phase1[i].izp[j]   = ct->site[i].atomtype[j];
            izpj                     = dftb->phase1[i].izp[j];
            dftb->phase1[i].ind[j+1] = dftb->phase1[i].ind[j] + dftb->lmax[izpj]* dftb->lmax[izpj];
        }
        dftb->phase1[i].ndim = dftb->phase1[i].norb;
        dftb->phase1[i].ne   = ct->site[i].extcharges;

        /* atom masses */
        mass = 0.0;
        for (j = 0; j < ct->site[i].atoms; j++)
        {
            dftb->phase1[i].mass[j] = mdatoms->massT[ct->site[i].atom[j]];
            mass                   += dftb->phase1[i].mass[j];
            printf("site %d  mass %d   %f\n", i, j,  mdatoms->massT[ct->site[i].atom[j]]);
        }
        dftb->phase1[i].inv_tot_mass = 1.0 / mass;
    }

    /* the remainder of phase2 */
    snew(dftb->phase2.x, dftb->phase2.nn);
    snew(dftb->phase2.xe, dftb->phase2.ne);
    snew(dftb->phase2.mass, dftb->phase2.nn);
    snew(dftb->phase2.ze, dftb->phase2.ne);
    snew(dftb->phase2.qmat, dftb->phase2.nn);
    snew(dftb->phase2.shift, dftb->phase2.nn);
    snew(dftb->phase2.shiftE, dftb->phase2.nn);
    snew(dftb->phase2.ind, dftb->phase2.nn + 1);
    snew(dftb->phase2.atind, dftb->phase2.nn + 1);
    snew(dftb->phase2.inf, ct->sites + 1);
    snew(dftb->phase2.ihomo, ct->sites + 1);
    snew(dftb->phase2.izp, dftb->phase2.nn);
    snew(dftb->phase2.hamil, dftb->phase2.norb);
    snew(dftb->phase2.hamil[0], SQR(dftb->phase2.norb));
    for (j = 1; j < dftb->phase2.norb; j++)
    {
        dftb->phase2.hamil[j] = dftb->phase2.hamil[0] + j * dftb->phase2.norb;
    }
    snew(dftb->phase2.overl, dftb->phase2.norb);
    snew(dftb->phase2.overl[0], SQR(dftb->phase2.norb));
    for (j = 1; j < dftb->phase2.norb; j++)
    {
        dftb->phase2.overl[j] = dftb->phase2.overl[0] + j * dftb->phase2.norb;
    }
    snew(dftb->phase2.Taf, dftb->phase2.norb);
    snew(dftb->phase2.Taf[0], SQR(dftb->phase2.norb));
    for (j = 1; j < dftb->phase2.norb; j++)
    {
        dftb->phase2.Taf[j] = dftb->phase2.Taf[0] + j * dftb->phase2.norb;
    }
/*
   snew(dftb->phase2.THamil, dftb->phase2.norb);
    snew(dftb->phase2.THamil[0], SQR(dftb->phase2.norb));
    for(j = 1; j < dftb->phase2.norb; j++)
      dftb->phase2.THamil[j] = dftb->phase2.THamil[0] + j * dftb->phase2.norb;
    snew(dftb->phase2.OverlF, dftb->phase2.norb);
    snew(dftb->phase2.OverlF[0], SQR(dftb->phase2.norb));
    for (j = 1; j < dftb->phase2.norb; j++)
    {
        dftb->phase2.OverlF[j] = dftb->phase2.OverlF[0] + j * dftb->phase2.norb;
    }

 */
    if (ct->do_lambda_i == 3 || ct->jobtype == cteJFSSH || ct->jobtype == cteBCJFSSH || ct->jobtype == cteNOMOVEMENT)
    {
        snew(dftb->phase2.b, dftb->phase2.norb);
        snew(dftb->phase2.b[0], SQR(dftb->phase2.norb));
        for (j = 1; j < dftb->phase2.norb; j++)
        {
            dftb->phase2.b[j] = dftb->phase2.b[0] + j * dftb->phase2.norb;
        }

        snew(dftb->phase2.grad, dftb->phase2.nn);
        snew(dftb->phase2.partgrad, dftb->phase2.nn);

    }
    //snew(dftb->phase2.overl_hybrid, dftb->phase2.norb);
    //  snew(dftb->phase2.overl_hybrid[0], SQR(dftb->phase2.norb));
    //  for(j = 1; j < dftb->phase2.norb; j++)
    //    dftb->phase2.overl_hybrid[j] = dftb->phase2.overl_hybrid[0] + j * dftb->phase2.norb;
    //dftb->phase2.THamilOrtho = (double **) malloc(dftb->phase2.norb * sizeof(double *));
    //  dftb->phase2.THamilOrtho[0] = (double *) malloc(SQR(dftb->phase2.norb) * sizeof(double));

/*
   snew(dftb->phase2.THamilOrtho, dftb->phase2.norb);
    snew(dftb->phase2.THamilOrtho[0], SQR(dftb->phase2.norb));
    for(j = 1; j < dftb->phase2.norb; j++)
      dftb->phase2.THamilOrtho[j] = dftb->phase2.THamilOrtho[0] + j * dftb->phase2.norb;
   snew(dftb->phase2.tij, dftb->phase2.norb);
    snew(dftb->phase2.tij[0], SQR(dftb->phase2.norb));
    for(j = 1; j < dftb->phase2.norb; j++)
      dftb->phase2.tij[j] = dftb->phase2.tij[0] + j * dftb->phase2.norb;
   snew(dftb->phase2.sij, dftb->phase2.norb);
    snew(dftb->phase2.sij[0], SQR(dftb->phase2.norb));
    for(j = 1; j < dftb->phase2.norb; j++)
      dftb->phase2.sij[j] = dftb->phase2.sij[0] + j * dftb->phase2.norb;
 */
//*
    snew(dftb->phase2.THamilOrtho, ct->dim);
    snew(dftb->phase2.THamilOrtho[0], SQR(ct->dim));
    for (j = 1; j < ct->dim; j++)
    {
        dftb->phase2.THamilOrtho[j] = dftb->phase2.THamilOrtho[0] + j * ct->dim;
    }
    snew(dftb->phase2.tij, ct->dim);
    snew(dftb->phase2.tij[0], SQR(ct->dim));
    for (j = 1; j < ct->dim; j++)
    {
        dftb->phase2.tij[j] = dftb->phase2.tij[0] + j * ct->dim;
    }
    snew(dftb->phase2.sij, ct->dim);
    snew(dftb->phase2.sij[0], SQR(ct->dim));
    for (j = 1; j < ct->dim; j++)
    {
        dftb->phase2.sij[j] = dftb->phase2.sij[0] + j * ct->dim;
    }
//*/
    snew(dftb->phase2.gammamat, ct->atoms_cplx);
    snew(dftb->phase2.gammamat[0], SQR(ct->atoms_cplx));
    for (j = 1; j < ct->atoms_cplx; j++)
    {
        dftb->phase2.gammamat[j] = dftb->phase2.gammamat[0] + j * ct->atoms_cplx;
    }
    snew(dftb->nl, ct->atoms_cplx);
    snew(dftb->nl[0], SQR(ct->atoms_cplx));
    for (j = 1; j < ct->atoms_cplx; j++)
    {
        dftb->nl[j] = dftb->nl[0] + j * ct->atoms_cplx;
    }
    snew(dftb->phase2.ev, dftb->phase2.norb);
    snew(dftb->phase2.occ, dftb->phase2.norb);
    snew(dftb->phase2.aux, 3 * dftb->phase2.norb);

    if (ct->qmmm == 3)
    {
        snew(dftb->phase2.neighbors_pme, ct->atoms_cplx);
        snew(dftb->phase2.neighbor_pme, ct->atoms_cplx);
    }

    dftb->phase2.inf[0]   = 0; // where do the orbitals of fragment i start?
    dftb->phase2.ihomo[0] = 0; //index of first homo of site i
    for (i = 0; i < ct->sites; i++)
    {
        dftb->phase2.inf[i+1]   = dftb->phase2.inf[i] + dftb->phase1[i].norb;
        dftb->phase2.ihomo[i+1] = dftb->phase2.ihomo[i] + ct->site[i].homos;
    }

    for (i = 0; i < ct->sites; i++)
    {
        for (j = 0; j < ct->site[i].atoms; j++)
        {
            dftb->phase1[i].qmat[j] = dftb->qzero1[dftb->phase1[i].izp[j]];
        }
    }
    counter               = 0;
    dftb->phase2.ind[0]   = 0;
    dftb->phase2.atind[0] = 0;
    mass                  = 0.0;
    for (i = 0; i < ct->sites; i++)
    {
        dftb->phase2.atind[i+1] = dftb->phase2.atind[i] + ct->site[i].atoms;
        for (j = 0; j < ct->site[i].atoms; j++)
        {
            dftb->phase2.izp[counter] = dftb->phase1[i].izp[j];
            izpj = dftb->phase2.izp[counter];
            dftb->phase2.ind[counter+1] = dftb->phase2.ind[counter] + dftb->lmax[izpj]* dftb->lmax[izpj];
            //dftb->phase2.mass[counter] = mdatoms->massT[ct->site[i].atom[j]]; /* mass */
            dftb->phase2.mass[counter] = dftb->phase1[i].mass[j]; /* mass */
            mass += dftb->phase2.mass[counter];
            counter++;
        }
    }
    dftb->phase2.inv_tot_mass = 1.0 / mass;
    if (counter != dftb->phase2.nn)
    {
        PRINTF("The number of atoms does not match: counter = %d, dftb->phase2.nn = %d !\n", counter, dftb->phase2.nn);
        exit(-1);
    }

    /* auxiliary arrays for function orthogonalize() */
    snew(dftb->orthogo.tij, ct->dim*ct->dim);
    snew(dftb->orthogo.sij, ct->dim*ct->dim);
    snew(dftb->orthogo.evec, ct->dim*ct->dim);
    dftb->orthogo.lwork = 26*ct->dim;
    snew(dftb->orthogo.work, 26*ct->dim*26*ct->dim);
    dftb->orthogo.liwork = 10*ct->dim;
    snew(dftb->orthogo.iwork, 10*ct->dim*10*ct->dim);
    snew(dftb->orthogo.eval, ct->dim);
    snew(dftb->orthogo.issupz, 2*ct->dim);
    /* also includes arrays for time derivative of FMOs - Weiwei X Mar 2019 */
    //  [naos][ct->dim]
    snew(dftb->orthogo.sij_old, ct->dim*ct->dim);

    snew(dftb->orthogo.fmo_old, dftb->phase2.norb);
    snew(dftb->orthogo.fmo_old[0], ct->dim * dftb->phase2.norb);
    for (i=1; i<dftb->phase2.norb; i++)
      dftb->orthogo.fmo_old[i] = dftb->orthogo.fmo_old[0] + i * ct->dim;

    snew(dftb->orthogo.fmo, dftb->phase2.norb);
    snew(dftb->orthogo.fmo[0], ct->dim * dftb->phase2.norb);
    for (i=1; i<dftb->phase2.norb; i++)
      dftb->orthogo.fmo[i] = dftb->orthogo.fmo[0] + i * ct->dim;
    // WRONG? for (i=1; i<ct->dim; i++)
    // WRONG?   dftb->orthogo.fmo[i] = dftb->orthogo.fmo[0] + i * dftb->phase2.norb;
    //  [naos][naos]
/*
    snew(dftb->orthogo.sao_old, dftb->phase2.norb);
    snew(dftb->orthogo.sao_old[0], SQR(dftb->phase2.norb));
    for (i=1; i<dftb->phase2.norb; i++)
      dftb->orthogo.sao_old[i] = dftb->orthogo.sao_old[0] + i * dftb->phase2.norb;

    snew(dftb->orthogo.sao, dftb->phase2.norb);
    snew(dftb->orthogo.sao[0], SQR(dftb->phase2.norb));
    for (i=1; i<dftb->phase2.norb; i++)
      dftb->orthogo.sao[i] = dftb->orthogo.sao[0] + i * dftb->phase2.norb;
*/
    /* arrays for state following */
    /*
       snew(dftb->orthogo.evec_ao, dftb->phase2.norb);
       snew(dftb->orthogo.evec_ao[0], dftb->phase2.norb * ct->dim);
       for(j = 1; j < dftb->phase2.norb; j++)
        dftb->orthogo.evec_ao[j] = dftb->orthogo.evec_ao[0] + j * ct->dim;
       snew(dftb->orthogo.evec_ao_old, dftb->phase2.norb);
       snew(dftb->orthogo.evec_ao_old[0], dftb->phase2.norb * ct->dim);
       for(j = 1; j < dftb->phase2.norb; j++)
        dftb->orthogo.evec_ao_old[j] = dftb->orthogo.evec_ao_old[0] + j * ct->dim;
       snew(dftb->orthogo.evec_ao_ref, dftb->phase2.norb);
       snew(dftb->orthogo.evec_ao_ref[0], dftb->phase2.norb * ct->dim);
       for(j = 1; j < dftb->phase2.norb; j++)
        dftb->orthogo.evec_ao_ref[j] = dftb->orthogo.evec_ao_ref[0] + j * ct->dim;
     */
    snew(dftb->orthogo.overlap, ct->dim);
    snew(dftb->orthogo.overlap[0], SQR(ct->dim));
    for (j = 1; j < ct->dim; j++)
    {
        dftb->orthogo.overlap[j] = dftb->orthogo.overlap[0] + j * ct->dim;
    }
    snew(dftb->orthogo.overlap_ref, ct->dim);
    snew(dftb->orthogo.overlap_ref[0], SQR(ct->dim));
    for (j = 1; j < ct->dim; j++)
    {
        dftb->orthogo.overlap_ref[j] = dftb->orthogo.overlap_ref[0] + j * ct->dim;
    }

    snew(ct->wf_old, 2*ct->dim);

    snew(dftb->overl_test, ct->dim);
    snew(dftb->overl_test[0], SQR(ct->dim));
    for (j = 1; j < ct->dim; j++)
    {
        dftb->overl_test[j] = dftb->overl_test[0] + j * ct->dim;
    }



    /* machine accuracy */
    dftb->racc = 1.0;
    while ((1.0 + dftb->racc) > 1.0)
    {
        dftb->racc /= 2.0;
    }
    dftb->racc *= 2.0;
    dftb->dacc  = 4 * dftb->racc;

#ifdef GMX_MPI
    printf("Completed DFTB initialization at rank %d\n", ct_mpi_rank);
#else
    printf("Completed DFTB initialization\n");
#endif

    return;
}

#ifdef GMX_MPI
void init_dftb_pme(dftb_t *dftb, charge_transfer_t *ct, t_inputrec *ir, int ct_mpi_rank)
#else
void init_dftb_pme(dftb_t *dftb, charge_transfer_t *ct, t_inputrec *ir)
#endif
{
    int i, j, status;
    dftb->ewaldcoeff_pme = calc_ewaldcoeff(ir->rcoulomb, ir->ewald_rtol);
    dftb->rcoulomb_pme   = ir->rcoulomb;
    dftb->nstlist_pme    = ir->nstlist;
    dftb->lastlist_pme   = -1000000; /* to have the list generated before the first PME calculation (in any case) */
    /* for phase 1 -- loop over all sites
     * POSSIBLE MODIFICATION FOR PARALLEL RUNS
     *   -- DO IT ONLY FOR THOSE SITES THAT ARE TO BE CALCULATED ON THIS RESPECTIVE CT_MPI_RANK!
     */
    for (i = 0; i < ct->sites; i++)
    {
        snew(dftb->phase1[i].q_pme, ct->site[i].atoms + ct->site[i].extcharges);
        snew(dftb->phase1[i].x_pme, ct->site[i].atoms + ct->site[i].extcharges);
        snew(dftb->phase1[i].pot, ct->site[i].atoms + ct->site[i].extcharges);
        snew(dftb->phase1[i].pot2, ct->site[i].atoms);
        snew(dftb->phase1[i].pot3, ct->site[i].atoms);
        snew(dftb->phase1[i].pot4, ct->site[i].atoms);
        snew(dftb->phase1[i].nrnb_pme, 1);
        for (j = 0; j < ct->site[i].atoms; j++)
        {
            dftb->phase1[i].qmat[j] = dftb->qzero1[dftb->phase1[i].izp[j]];
        }
        snew(dftb->phase1[i].pmedata, 1);
        status = gmx_pme_init_dftb(dftb->phase1[i].pmedata, ir, ct->site[i].atoms + ct->site[i].extcharges);
    }
    /* for phase 2 -- POSSIBLY ONLY ON CT_MPI_RANK==0
     * ! */
    snew(dftb->phase2.q_pme, ct->atoms_cplx + ct->extcharges_cplx);
    snew(dftb->phase2.x_pme, ct->atoms_cplx + ct->extcharges_cplx);
    snew(dftb->phase2.pot, ct->atoms_cplx + ct->extcharges_cplx);
    snew(dftb->phase2.pot2, ct->atoms_cplx);
    snew(dftb->phase2.pot3, ct->atoms_cplx);
    snew(dftb->phase2.pot4, ct->atoms_cplx);
    snew(dftb->phase2.nrnb_pme, 1);
    snew(dftb->phase2.pmedata, 1);
    status = gmx_pme_init_dftb(dftb->phase2.pmedata, ir, ct->atoms_cplx + ct->extcharges_cplx);
    return;
}

/****************************
 * INITIALIZE THE DIIS CODE *
 ****************************/

void ct_init_diis(charge_transfer_t *ct, ct_diis_t *diis)
{
    int i;

    diis->n_elem = ct->dim;
    //diis->prev_q_input = (double**) malloc(DIIS_MAX_PREV_VECTORS * sizeof(double *));
    //  diis->prev_q_input[0] = (double*) malloc(diis->n_elem * DIIS_MAX_PREV_VECTORS);
    snew(diis->prev_q_input, DIIS_MAX_PREV_VECTORS);
    snew(diis->prev_q_input[0], diis->n_elem * DIIS_MAX_PREV_VECTORS);
    for (i = 1; i < DIIS_MAX_PREV_VECTORS; i++)
    {
        diis->prev_q_input[i] = diis->prev_q_input[0] + i * diis->n_elem;
    }
    //diis->prev_q_diff = (double**) malloc(DIIS_MAX_PREV_VECTORS * sizeof(double *));
    //  diis->prev_q_diff[0] = (double*) malloc(diis->n_elem * DIIS_MAX_PREV_VECTORS);
    snew(diis->prev_q_diff, DIIS_MAX_PREV_VECTORS);
    snew(diis->prev_q_diff[0], diis->n_elem * DIIS_MAX_PREV_VECTORS);
    for (i = 1; i < DIIS_MAX_PREV_VECTORS; i++)
    {
        diis->prev_q_diff[i] = diis->prev_q_diff[0] + i * diis->n_elem;
    }
    //diis->aa = (double *) malloc(SQR(DIIS_MAX_PREV_VECTORS + 1));
    //diis->bb = (double *) malloc(DIIS_MAX_PREV_VECTORS + 1);
    //diis->q_inp_result = (double *) malloc(ct->sites * sizeof(double));
    //diis->q_diff = (double *) malloc(ct->sites * sizeof(double));
    snew(diis->aa, SQR(DIIS_MAX_PREV_VECTORS + 1));
    snew(diis->bb, DIIS_MAX_PREV_VECTORS + 1);
    //snew(diis->q_inp_result, ct->sites);
    //snew(diis->q_diff, ct->sites);
    snew(diis->fermi_coef, ct->dim);

    return;
}

void ct_init_broyden(charge_transfer_t *ct, dftb_broyden_t *broyd)
{
    snew(broyd->f,      ct->dim);
    snew(broyd->ui,     ct->dim);
    snew(broyd->vti,    ct->dim);
    snew(broyd->t1,     ct->dim);
    snew(broyd->dumvi,  ct->dim);
    snew(broyd->df,     ct->dim);
    snew(broyd->vector, ct->dim);
    snew(broyd->unit31, ct->dim);
    snew(broyd->unit32, ct->dim);

    return;
}

/******************************************
 * PREREQUISITIES FOR A QM/MM CALCULATION *
 *    TO BE PERFORMED IN EVERY MD STEP    *
 ******************************************/
//do_SCC?????????
void prepare_charge_transfer(matrix state_box, t_mdatoms *mdatoms, dftb_t *dftb, charge_transfer_t *ct, rvec *x_ct)
{
    int           i, j, k, l, m, n, counter, counter2;
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

#ifdef GMX_MPI
void do_pme_for_dftb(charge_transfer_t *ct, dftb_t *dftb, int ct_mpi_rank, int ct_mpi_size)
#else
void do_pme_for_dftb(charge_transfer_t *ct, dftb_t *dftb)
#endif
{
    int           i, j, k, nn, ne, status;
    dftb_phase1_t dftb1;
    dftb_phase2_t dftb2;
    real          energy_pme, normbond;
    rvec          bond, x_qm, x_mm;
    dvec          dbond;
    t_pbc         pbc;

    set_pbc(&pbc, epbcXYZ, dftb->box_pme);

    for (i = 0; i < ct->sites; i++)
#ifdef GMX_MPI
    {
        if (i % ct_mpi_size == ct_mpi_rank)
#endif
    {
        dftb1 = dftb->phase1[i];
        nn    = dftb1.nn;
        ne    = dftb1.ne;
        /* do PME for DFTB here!
         * do it with the atomic charges obtained in the previous MD step
         * it is possibly slightly inaccurate,
         * but probably OK as a first guess */

        for (j = 0; j < nn; j++)
        {
            dftb1.pot[j] = 0.0;
        }

        /* PME -- long-range component (reciprocal space) */
        init_nrnb(dftb->phase1[i].nrnb_pme);
        gmx_pme_do_dftb(*(dftb->phase1[i].pmedata), 0, nn+ne, dftb1.x_pme, dftb1.q_pme, dftb->box_pme,
                        dftb->phase1[i].nrnb_pme, dftb->ewaldcoeff_pme, &energy_pme, dftb1.pot);

        /* PME -- corrections */
        for (j = 0; j < nn; j++)
        {
            /* exclude the QM-QM interactions as the shift will be calculated in DFTB for these interactions */
            dftb1.pot2[j] = 0.;
            for (k = 0; k < nn; k++)
            {
                if (j != k)
                {
                    dvec_sub(dftb1.x[j], dftb1.x[k], dbond);
                    dftb1.pot2[j] -= (-dftb1.qmat[k] + dftb->qzero1[dftb1.izp[k]])
                        * gmx_erf(dftb->ewaldcoeff_pme * dnorm(dbond) / NM_TO_BOHR) / dnorm(dbond); // this should be OK
                }
            }
            /* the self-interaction of charge densities */
            dftb1.pot3[j] = -dftb->ewaldcoeff_pme / NM_TO_BOHR * (-dftb1.qmat[j] + dftb->qzero1[dftb1.izp[j]]) / sqrt(M_PI); // this should be OK
        }

        /* PME -- short-range component (real space) using the previously created neighbor list */
        for (j = 0; j < nn; j++)
        {
            dftb1.pot4[j] = 0.;
            x_qm[XX]      = (real) dftb1.x[j][XX] / NM_TO_BOHR;
            x_qm[YY]      = (real) dftb1.x[j][YY] / NM_TO_BOHR;
            x_qm[ZZ]      = (real) dftb1.x[j][ZZ] / NM_TO_BOHR;
            /* loop over the neighbors */
            for (k = 0; k < dftb1.neighbors_pme[j]; k++)
            {
                x_mm[XX] = (real) dftb1.xe[dftb1.neighbor_pme[j][k]][XX] / NM_TO_BOHR;
                x_mm[YY] = (real) dftb1.xe[dftb1.neighbor_pme[j][k]][YY] / NM_TO_BOHR;
                x_mm[ZZ] = (real) dftb1.xe[dftb1.neighbor_pme[j][k]][ZZ] / NM_TO_BOHR;
                /* calculate the distance (considering PBC) */
                status   = pbc_dx_aiuc(&pbc, x_qm, x_mm, bond);
                normbond = norm(bond);
                if (normbond < dftb->rcoulomb_pme)
                {
                    dftb1.pot4[j] += dftb1.ze[dftb1.neighbor_pme[j][k]] / normbond * gmx_erfc(dftb->ewaldcoeff_pme * normbond);
                }
            }
        }

        //printf("Ewald potential kJ/mol/e = external shift\n");
        for (j = 0; j < nn; j++)
        {
            /* convert kJ/mol/e to atomic units of ESP */
            dftb1.pot[j] *= KJMOL_TO_HARTREE;
            if (i == 0 && j == 0)
            {
                printf("%12.7f %12.7f %12.7f %12.7f %12.7f ", dftb1.pot[j], dftb1.pot2[j], dftb1.pot3[j], dftb1.pot4[j] / NM_TO_BOHR, -dftb1.qmat[j] + dftb->qzero1[dftb1.izp[j]]);
            }
            dftb1.pot[j] += dftb1.pot2[j] + dftb1.pot3[j] + dftb1.pot4[j] / NM_TO_BOHR;
            //printf("%12.7f\n", dftb1.pot[j]);
        }

        /* save some an average value in the variable esp */
        dftb->phase1[i].esp = 0.;
        for (j = 0; j < nn; j++)
        {
            dftb->phase1[i].esp += dftb1.pot[j] * dftb1.mass[j];
        }
        dftb->phase1[i].esp *= dftb1.inv_tot_mass / ct->esp_scaling_factor;
        /* the following would be a very approximative alternative...
         * dftb->phase1[i].esp = dftb1.pot[0];
         */

        /* result is in: dftb.phase1[i].pot */
    }
#ifdef GMX_MPI
    }
    if (ct_mpi_rank == 0)
    {
#endif
    /* do it here for the complex */
    /* POSSIBLE OPTIMIZATION -- USE THE DATA OBTAINED FOR THE FRAGMENTS?
     * BUT ATTENTION: NOT ALL OF THE DATA HAVE BEEN CALCULATED ON MPI_RANK==0 !!!
     */
    dftb2 = dftb->phase2;
    nn    = dftb2.nn;
    ne    = dftb2.ne;

    for (j = 0; j < nn; j++)
    {
        dftb2.pot[j] = 0.0;
    }

    /* PME -- long-range component (reciprocal space) */
    init_nrnb(dftb->phase2.nrnb_pme);
    gmx_pme_do_dftb(*(dftb->phase2.pmedata), 0, nn+ne, dftb2.x_pme, dftb2.q_pme, dftb->box_pme,
                    dftb->phase2.nrnb_pme, dftb->ewaldcoeff_pme, &energy_pme, dftb2.pot);

    /* PME -- corrections */
    for (j = 0; j < nn; j++)
    {
        /* exclude the QM-QM interactions as the shift will be calculated in DFTB for these interactions */
        dftb2.pot2[j] = 0.;
        for (k = 0; k < nn; k++)
        {
            if (j != k)
            {
                dvec_sub(dftb2.x[j], dftb2.x[k], dbond);
                dftb2.pot2[j] -= (-dftb2.qmat[k] + dftb->qzero2[dftb2.izp[k]])
                    * gmx_erf(dftb->ewaldcoeff_pme * dnorm(dbond) / NM_TO_BOHR) / dnorm(dbond);
            }
        }
        /* the self-interaction of charge densities */
        dftb2.pot3[j] = -dftb->ewaldcoeff_pme / NM_TO_BOHR * (-dftb2.qmat[j] + dftb->qzero2[dftb2.izp[j]]) / sqrt(M_PI);
    }

    /* PME -- short-range component (real space) using the previously created neighbor list */
    for (j = 0; j < nn; j++)
    {
        dftb2.pot4[j] = 0.;
        x_qm[XX]      = (real) dftb2.x[j][XX] / NM_TO_BOHR;
        x_qm[YY]      = (real) dftb2.x[j][YY] / NM_TO_BOHR;
        x_qm[ZZ]      = (real) dftb2.x[j][ZZ] / NM_TO_BOHR;
        /* loop over the neighbors */
        for (k = 0; k < dftb2.neighbors_pme[j]; k++)
        {
            x_mm[XX] = (real) dftb2.xe[dftb2.neighbor_pme[j][k]][XX] / NM_TO_BOHR;
            x_mm[YY] = (real) dftb2.xe[dftb2.neighbor_pme[j][k]][YY] / NM_TO_BOHR;
            x_mm[ZZ] = (real) dftb2.xe[dftb2.neighbor_pme[j][k]][ZZ] / NM_TO_BOHR;
            /* calculate the distance (considering PBC) */
            status   = pbc_dx_aiuc(&pbc, x_qm, x_mm, bond);
            normbond = norm(bond);
            if (normbond < dftb->rcoulomb_pme)
            {
                dftb2.pot4[j] += dftb2.ze[dftb2.neighbor_pme[j][k]] / normbond * gmx_erfc(dftb->ewaldcoeff_pme * normbond);
            }
        }
    }

    //printf("Ewald potential kJ/mol/e = external shift\n");
    for (j = 0; j < nn; j++)
    {
        /* convert kJ/mol/e to atomic units of ESP */
        dftb2.pot[j] *= KJMOL_TO_HARTREE;
        //printf("%12.7f %12.7f %12.7f %12.7f ", dftb1.pot[j], dftb1.pot2[j], dftb1.pot3[j], dftb1.pot4[j]);
        dftb2.pot[j] += dftb2.pot2[j] + dftb2.pot3[j] + dftb2.pot4[j] / NM_TO_BOHR;
        //printf("%12.7f\n", dftb1.pot[j]);
    }
    /* result is in: dftb.phase2.pot */
#ifdef GMX_MPI
}
#endif
    return;
}

#ifdef GMX_MPI
void do_pme_for_dftb_part1(charge_transfer_t *ct, dftb_t *dftb, int ct_mpi_rank, int ct_mpi_size)
#else
void do_pme_for_dftb_part1(charge_transfer_t *ct, dftb_t *dftb)
#endif
{
    int           i, j, k, nn, ne, status;
    dftb_phase1_t dftb1;
    dftb_phase2_t dftb2;
    real          energy_pme, normbond;
    rvec          bond, x_qm, x_mm;
    dvec          dbond;
    t_pbc         pbc;

    set_pbc(&pbc, epbcXYZ, dftb->box_pme);

    for (i = 0; i < ct->sites; i++)
#ifdef GMX_MPI
    {
        if (i % ct_mpi_size == ct_mpi_rank)
#endif
    {
        dftb1 = dftb->phase1[i];
        nn    = dftb1.nn;
        ne    = dftb1.ne;
        /* do PME for DFTB here! -- "preparatory" part, before actual SCC-DFTB calculation */

        /* PME -- short-range component (real space) using the previously created neighbor list */
        for (j = 0; j < nn; j++)
        {
            dftb1.pot4[j] = 0.;
            x_qm[XX]      = (real) dftb1.x[j][XX] / NM_TO_BOHR;
            x_qm[YY]      = (real) dftb1.x[j][YY] / NM_TO_BOHR;
            x_qm[ZZ]      = (real) dftb1.x[j][ZZ] / NM_TO_BOHR;
            /* loop over the neighbors */
            for (k = 0; k < dftb1.neighbors_pme[j]; k++)
            {
                x_mm[XX] = (real) dftb1.xe[dftb1.neighbor_pme[j][k]][XX] / NM_TO_BOHR;
                x_mm[YY] = (real) dftb1.xe[dftb1.neighbor_pme[j][k]][YY] / NM_TO_BOHR;
                x_mm[ZZ] = (real) dftb1.xe[dftb1.neighbor_pme[j][k]][ZZ] / NM_TO_BOHR;
                /* calculate the distance (considering PBC) */
                status   = pbc_dx_aiuc(&pbc, x_qm, x_mm, bond);
                normbond = norm(bond);
                if (normbond < dftb->rcoulomb_pme)
                {
                    dftb1.pot4[j] += dftb1.ze[dftb1.neighbor_pme[j][k]] / normbond * gmx_erfc(dftb->ewaldcoeff_pme * normbond);
                }
            }
        }
    }
#ifdef GMX_MPI
    }
#endif
    return;
}

void do_pme_for_dftb_part2(charge_transfer_t *ct, dftb_t *dftb, int i)
{
    int           j, k, nn, ne, status;
    dftb_phase1_t dftb1;
    real          energy_pme, normbond;
    rvec          bond, x_qm, x_mm;
    dvec          dbond;
    t_pbc         pbc;

    set_pbc(&pbc, epbcXYZ, dftb->box_pme);

    dftb1 = dftb->phase1[i];
    nn    = dftb1.nn;
    ne    = dftb1.ne;

    /* transfer the QM charges to PME,
     * and reset the ESP array
     */
    for (j = 0; j < nn; j++)
    {
        dftb1.pot[j]   = 0.0;
        dftb1.q_pme[j] = -dftb1.qmat[j] + dftb->qzero1[dftb1.izp[j]];
    }

    /* PME calculation */
    init_nrnb(dftb->phase1[i].nrnb_pme);
    gmx_pme_do_dftb(*(dftb->phase1[i].pmedata), 0, nn+ne, dftb1.x_pme, dftb1.q_pme, dftb->box_pme,
                    dftb->phase1[i].nrnb_pme, dftb->ewaldcoeff_pme, &energy_pme, dftb1.pot);

    /* PME -- corrections */
    for (j = 0; j < nn; j++)
    {
        /* exclude the QM-QM interactions */
        dftb1.pot2[j] = 0.;
        for (k = 0; k < nn; k++)
        {
            if (j != k)
            {
                dvec_sub(dftb1.x[j], dftb1.x[k], dbond);
                dftb1.pot2[j] -= (-dftb1.qmat[k] + dftb->qzero1[dftb1.izp[k]])
                    * gmx_erf(dftb->ewaldcoeff_pme * dnorm(dbond) / NM_TO_BOHR) / dnorm(dbond);
            }
        }
        /* exclude the self-interaction of the charge densities */
        dftb1.pot3[j] = -dftb->ewaldcoeff_pme / NM_TO_BOHR * (-dftb1.qmat[j] + dftb->qzero1[dftb1.izp[j]]) / sqrt(M_PI);
    }

    /* dftb1.pot4 was pre-calculated */

    //printf("Ewald potential kJ/mol/e = external shift\n");
    for (j = 0; j < nn; j++)
    {
        /* convert kJ/mol/e to atomic units of ESP */
        dftb1.pot[j] *= KJMOL_TO_HARTREE;
        //if (i==0 && j==0)
        //  printf("%12.7f %12.7f %12.7f %12.7f %12.7f ", dftb1.pot[j], dftb1.pot2[j], dftb1.pot3[j], dftb1.pot4[j] / NM_TO_BOHR, -dftb1.qmat[j] + dftb->qzero1[dftb1.izp[j]]);
        dftb1.pot[j] += dftb1.pot2[j] + dftb1.pot3[j] + dftb1.pot4[j] / NM_TO_BOHR;
        //if (i==0 && j==0) printf("%12.7f ", dftb1.pot[j]);
    }

    /* save some an average value in the variable esp */
    dftb->phase1[i].esp = 0.;
    for (j = 0; j < nn; j++)
    {
        dftb->phase1[i].esp += dftb1.pot[j] * dftb1.mass[j];
    }
    dftb->phase1[i].esp *= dftb1.inv_tot_mass / ct->esp_scaling_factor;
    /* this was too approximative:
       dftb->phase1[i].esp = dftb1.pot[0];
     */

    return;
}

/* PME for DFTB phase 2 -- this will be done only once per MD step */
void do_pme_for_dftb_phase2(charge_transfer_t *ct, dftb_t *dftb)
{
    int           i, j, k, nn, ne, status;
    dftb_phase2_t dftb2;
    real          energy_pme, normbond;
    rvec          bond, x_qm, x_mm;
    dvec          dbond;
    t_pbc         pbc;

    set_pbc(&pbc, epbcXYZ, dftb->box_pme);

    dftb2 = dftb->phase2;
    nn    = dftb2.nn;
    ne    = dftb2.ne;

    for (j = 0; j < nn; j++)
    {
        dftb2.q_pme[j] = -dftb2.qmat[j] + dftb->qzero2[dftb2.izp[j]];
        dftb2.pot[j]   = 0.0;
    }

    /* PME -- long-range component (reciprocal space) */
    init_nrnb(dftb->phase2.nrnb_pme);
    gmx_pme_do_dftb(*(dftb->phase2.pmedata), 0, nn+ne, dftb2.x_pme, dftb2.q_pme, dftb->box_pme,
                    dftb->phase2.nrnb_pme, dftb->ewaldcoeff_pme, &energy_pme, dftb2.pot);

    /* PME -- corrections */
    for (j = 0; j < nn; j++)
    {
        /* exclude the QM-QM interactions as the shift will be calculated in DFTB for these interactions */
        dftb2.pot2[j] = 0.;
        for (k = 0; k < nn; k++)
        {
            if (j != k)
            {
                dvec_sub(dftb2.x[j], dftb2.x[k], dbond);
                dftb2.pot2[j] -= (-dftb2.qmat[k] + dftb->qzero2[dftb2.izp[k]])
                    * gmx_erf(dftb->ewaldcoeff_pme * dnorm(dbond) / NM_TO_BOHR) / dnorm(dbond);
            }
        }
        /* the self-interaction of charge densities */
        dftb2.pot3[j] = -dftb->ewaldcoeff_pme / NM_TO_BOHR * (-dftb2.qmat[j] + dftb->qzero2[dftb2.izp[j]]) / sqrt(M_PI);
    }

    /* PME -- short-range component (real space) using the previously created neighbor list */
    for (j = 0; j < nn; j++)
    {
        dftb2.pot4[j] = 0.;
        x_qm[XX]      = (real) dftb2.x[j][XX] / NM_TO_BOHR;
        x_qm[YY]      = (real) dftb2.x[j][YY] / NM_TO_BOHR;
        x_qm[ZZ]      = (real) dftb2.x[j][ZZ] / NM_TO_BOHR;
        /* loop over the neighbors */
        for (k = 0; k < dftb2.neighbors_pme[j]; k++)
        {
            x_mm[XX] = (real) dftb2.xe[dftb2.neighbor_pme[j][k]][XX] / NM_TO_BOHR;
            x_mm[YY] = (real) dftb2.xe[dftb2.neighbor_pme[j][k]][YY] / NM_TO_BOHR;
            x_mm[ZZ] = (real) dftb2.xe[dftb2.neighbor_pme[j][k]][ZZ] / NM_TO_BOHR;
            /* calculate the distance (considering PBC) */
            status   = pbc_dx_aiuc(&pbc, x_qm, x_mm, bond);
            normbond = norm(bond);
            if (normbond < dftb->rcoulomb_pme)
            {
                dftb2.pot4[j] += dftb2.ze[dftb2.neighbor_pme[j][k]] / normbond * gmx_erfc(dftb->ewaldcoeff_pme * normbond);
            }
        }
    }

    //printf("Ewald potential kJ/mol/e = external shift\n");
    for (j = 0; j < nn; j++)
    {
        /* convert kJ/mol/e to atomic units of ESP */
        dftb2.pot[j] *= KJMOL_TO_HARTREE;
        //printf("%12.7f %12.7f %12.7f %12.7f ", dftb1.pot[j], dftb1.pot2[j], dftb1.pot3[j], dftb1.pot4[j]);
        dftb2.pot[j] += dftb2.pot2[j] + dftb2.pot3[j] + dftb2.pot4[j] / NM_TO_BOHR;
        //printf("%12.7f\n", dftb1.pot[j]);
    }
    /* result is in: dftb.phase2.pot */
    return;
}

#ifdef GMX_MPI
void do_neighborlist_for_dftb(charge_transfer_t *ct, dftb_t *dftb, rvec *x, int ct_mpi_rank, int ct_mpi_size)
#else
void do_neighborlist_for_dftb(charge_transfer_t *ct, dftb_t *dftb, rvec *x)
#endif
{
    int           i, j, k, nn, ne, counter, status;
    dftb_phase1_t dftb1;
    dftb_phase2_t dftb2;
    rvec          bond;
    t_pbc         pbc;
    real          rcoulomb2;

    set_pbc(&pbc, epbcXYZ, dftb->box_pme);

//printf("pbc %d\n", pbc.ePBCDX);
/*
    switch (pbc.ePBCDX)
    {
        case epbcdxRECTANGULAR:
   printf("tric\n");
        case epbcdxTRICLINIC:
   printf("rect\n");
     }
 */


    //for (i = 0; i < DIM; i++)
    //printf("hbox %f fbox %f mhbox %f \n", pbc.hbox_diag[i], pbc.fbox_diag[i],  pbc.mhbox_diag[i]);
    rcoulomb2 = dftb->rcoulomb_pme * dftb->rcoulomb_pme;

    //printf("  creating PME DFTB neighborlists...");
    for (i = 0; i < ct->sites; i++)
#ifdef GMX_MPI
    {
        if (i % ct_mpi_size == ct_mpi_rank)
#endif
    {
        dftb1 = dftb->phase1[i];
        nn    = dftb1.nn;
        ne    = dftb1.ne;
        /* create a neighborlist for PME calculation for QM/MM DFTB
         * in a naive O(N^2) way:
         * double loop over QM atoms and over MM atoms */

        /* do it for every QM atom */
        for (j = 0; j < nn; j++)
        {
            counter = 0;
            /* check every MM atom */
            for (k = 0; k < ne; k++)
            {
                status = pbc_dx_aiuc(&pbc, x[ct->site[i].atom[j]], x[ct->site[i].extcharge[k]], bond);
                if (norm2(bond) < rcoulomb2)
                {
                    dftb1.neighbor_pme[j][counter] = k;
                    counter++;
                    // printf("Fragment %d, atom %d: neighbor no. %d is extcharge %d\n", i+1, j+1, counter, k+1);
                    if (counter == MAX_PME_NEIGHBORS)
                    {
                        fprintf(stderr, "Too many PME neighbors for atom %d of fragment %d\n  Exiting !!!\n\n", j+1, i+1);
                        exit(-1);
                    }
                }
            }
//      printf("bond %17.8f \n", norm2(bond));
            dftb1.neighbors_pme[j] = counter;
            //fprintf(stderr, "NS for PME/DFTB: atom %d in fragment %d has %d MM neighbors\n", j+1, i+1, counter);
        }
    }
#ifdef GMX_MPI
    }
#endif

#ifdef GMX_MPI
    if (ct_mpi_rank == 0)
    {
#endif
    /* do it here for the complex */
    dftb2 = dftb->phase2;
    nn    = dftb2.nn;
    ne    = dftb2.ne;
    /* do it for every QM atom */
    for (j = 0; j < nn; j++)
    {
        counter = 0;
        /* check every MM atom */
        for (k = 0; k < ne; k++)
        {
            status = pbc_dx_aiuc(&pbc, x[ct->atom_cplx[j]], x[ct->extcharge_cplx[k]], bond);
            if (norm2(bond) < rcoulomb2)
            {
                dftb2.neighbor_pme[j][counter] = k;
                counter++;
                if (counter == MAX_PME_NEIGHBORS)
                {
                    fprintf(stderr, "Too many PME neighbors for atom %d of complex\n  Exiting !!!\n\n", j+1);
                    exit(-1);
                }
            }
        }
        dftb2.neighbors_pme[j] = counter;
        //fprintf(stderr, "NS for PME/DFTB: atom %d in complex has %d MM neighbors\n", j+1, counter);
    }
#ifdef GMX_MPI
}
#endif
    return;
}

#ifdef GMX_MPI
void do_dftb_phase1(charge_transfer_t *ct, dftb_t *dftb, MPI_Comm ct_mpi_comm, int ct_mpi_rank, int ct_mpi_size)
{
    int i, return_value;

    printf("do_dftb_phase1 start at rank %d at %f\n", ct_mpi_rank, (double) clock()/CLOCKS_PER_SEC);

    for (i = 0; i < ct->sites; i++)
    {
        if (i % ct_mpi_size == ct_mpi_rank)
        {
            //printf("Doing residue %d at rank %d\n", ct->site[i].resnr, ct_mpi_rank);
            run_dftb1(ct, dftb, i);
        }
    }


    printf("do_dftb_phase1 end at rank %d at %f\n", ct_mpi_rank, (double) clock()/CLOCKS_PER_SEC);
    return;
}
void get_MM_params(charge_transfer_t *ct, dftb_t *dftb, MPI_Comm ct_mpi_comm, int ct_mpi_rank, int ct_mpi_size)
{
    int i;
    for (i = 0; i < ct->sites; i++)
    {
        if (i % ct_mpi_size == ct_mpi_rank)
        {
            if (ct->delta_q_mode == 0)
            {
                get_delta_q(dftb, ct, i);
            }
            if (!ct->hamiltonian_type && (ct->do_lambda_i == 2 || ct->do_lambda_i == 3))
            {
                get_internal_forces(dftb, ct, i);
            }
        }
    }
    printf("MM params part 1 end   at %f\n", (double) clock()/CLOCKS_PER_SEC);
    if (!ct->hamiltonian_type && (ct->do_lambda_i == 3 && ct_mpi_rank == 0))
    {
       offdiag_gradient_homo(dftb, dftb->phase2.x, dftb->phase2.grad, ct);
       //for (i=0; i<dftb->phase2.nn; i++)
       //    printf("offdiag grad at atom i %d  QM force %lf \n", i,  -(real) HARTREE_BOHR2MD * fabs(dftb->phase2.grad[i][0])+fabs(dftb->phase2.grad[i][1])+fabs(dftb->phase2.grad[i][2]));
       printf("MM params part 2 end   at %f\n", (double) clock()/CLOCKS_PER_SEC);
    }

    return;
}

void do_esp_only(charge_transfer_t *ct, dftb_t *dftb, real *q, MPI_Comm ct_mpi_comm, int ct_mpi_rank, int ct_mpi_size)
{
    int i, j, return_value;

    printf("ESP not implemented yet!\n");
    exit(-1);
/*
   if (ct->qmmm == 3) {
    for (i=0; i<ct->sites; i++) {
      if (i % ct_mpi_size == ct_mpi_rank) {
        // set the charges on "QM" atoms
        for (j=0; j<dftb->phase1[i].nn; j++) {
          dftb->phase1[i].qmat[j] = - q[ct->site[i].atom[j]] + dftb->qzero1[dftb->phase1[i].izp[j]];
        }
        dftb->phase1[i].qmat[0] += - q[ct->site[i].atom[0]+1];
        if (ct->modif_extcharge[i][0] > -1)
          dftb->phase1[i].qmat[0] -= - (ct->is_dna[i] ? EXTCHARGE_SHIFT_DNA : EXTCHARGE_SHIFT_TRP); //muss noch angepasst werden
        if (ct->modif_extcharge[i][1] > -1)
          dftb->phase1[i].qmat[0] -= - (ct->is_dna[i] ? EXTCHARGE_SHIFT_DNA : EXTCHARGE_SHIFT_TRP);
        // calculate the remaining components of PME
        do_pme_for_dftb_part2(ct, dftb, i);
        if (ct_mpi_rank != 0) {
          return_value = MPI_Send(&(dftb->phase1[i].esp), 1, MPI_DOUBLE, 0, 103 + i * 4, ct_mpi_comm);
        }
      }
    }
   }
 */
    return;
}
#else
void do_dftb_phase1(charge_transfer_t *ct, dftb_t *dftb)
{
    int i;

    printf("do_dftb_phase1 start at %f\n", (double) clock()/CLOCKS_PER_SEC);

    for (i = 0; i < ct->sites; i++)
    {
        run_dftb1(ct, dftb, i);
    }

    printf("do_dftb_phase1 end   at %f\n", (double) clock()/CLOCKS_PER_SEC);

    return;
}
void get_MM_params(charge_transfer_t *ct, dftb_t *dftb)
{
    int i;
    for (i = 0; i < ct->sites; i++)
    {
        if (ct->delta_q_mode == 0)
        {
            get_delta_q(dftb, ct, i);
        }
        if (!ct->hamiltonian_type && (ct->do_lambda_i == 2 || ct->do_lambda_i == 3))
        {
            get_internal_forces(dftb, ct, i);
        }
    }

    if (ct->do_lambda_i == 3 && !ct->hamiltonian_type)
    {
       if ((ct->step % JFSSH_OFFDIAG_FORCE == 0 && ct->jobtype == cteJFSSH || ct->jobtype == cteBCJFSSH) ||
           (ct->step % (JFSSH_OFFDIAG_FORCE/10) == 0 && ct->jobtype == cteSCCDYNAMIC) ||
           (dftb->do_ncv == 1))
       {
           offdiag_gradient_homo(dftb, dftb->phase2.x, dftb->phase2.grad, ct);
           //for (i=0; i<dftb->phase2.nn; i++)
           //    printf("offdiag grad at atom i %d  QM force %lf \n", i,  -(real) HARTREE_BOHR2MD * fabs(dftb->phase2.grad[i][0])+fabs(dftb->phase2.grad[i][1])+fabs(dftb->phase2.grad[i][2]));
           printf("offdiag forces end   at %f\n", (double) clock()/CLOCKS_PER_SEC);
    }
    }

    return;
}

void do_esp_only(charge_transfer_t *ct, dftb_t *dftb, real *q)
{
    int i, j;

    printf("ESP not implimented yet!\n");
    exit(-1);
/*
   if (ct->qmmm == 3) {
      for (i=0; i<ct->sites; i++) {
        // set the charges on "QM" atoms
        for (j=0; j<dftb->phase1[i].nn; j++) {
          dftb->phase1[i].qmat[j] = - q[ct->atom[i][j]] + dftb->qzero1[dftb->phase1[i].izp[j]];
        }
        dftb->phase1[i].qmat[0] += - q[ct->atom[i][0]+1];
        if (ct->modif_extcharge[i][0] > -1)
          dftb->phase1[i].qmat[0] -= - (ct->is_dna[i] ? EXTCHARGE_SHIFT_DNA : EXTCHARGE_SHIFT_TRP);
        if (ct->modif_extcharge[i][1] > -1)
          dftb->phase1[i].qmat[0] -= - (ct->is_dna[i] ? EXTCHARGE_SHIFT_DNA : EXTCHARGE_SHIFT_TRP);
        // calculate the remaining components of PME
        do_pme_for_dftb_part2(ct, dftb, i);
      }
   }
 */
    return;
}
#endif

void do_dftb_phase2(charge_transfer_t *ct, dftb_t *dftb)
{
    int i, l;
    printf("do_dftb_phase2 start at %f\n", (double) clock()/CLOCKS_PER_SEC);

    run_dftb2(ct, dftb);

    printf("do_dftb_phase2 end   at %f\n", (double) clock()/CLOCKS_PER_SEC);

    return;
}

t_atoms* protein_preprocessor(t_atoms *atoms, t_state *state_global)
{
    int      i, j, rescounter;
    t_atoms *atoms2, *atoms3;
    /* build copy of ct_atoms */
    //snew(atoms2,1);
    //init_t_atoms(atoms2, atoms->nr, FALSE);
    atoms2 = copy_t_atoms(atoms);

    /* take a different part in the memory for the copied atomnames */
    snew(atoms2->atomname, atoms->nr);
    for (i = 0; i < atoms->nr; i++)
    {
        snew(atoms2->atomname[i], 1);
    }
    for (i = 0; i < atoms->nr; i++)
    {
        snew(atoms2->atomname[i][0], 6);
    }
//  for (i=0; i<atoms->nr; i++)
//  atoms2->atomname[i] = atoms2->atomname[0] + 6*i;
//  for (i=0; i<atoms->nr; i++)
//  sfree(atoms2->atomname[i]);
//  for (i=0; i<atoms->nr; i++)
//  snew(atoms2->atomname[i], atoms->nr);
//  for (i=0; i<atoms->nr; i++)
//  snew(atoms2->atomname[i][0], 6);
//atoms2->atomname[0][0]=&atoms2->atomname[0][0][0];

    printf("pointer %p   %p\n", atoms->atomname[2], atoms2->atomname[2]);
    printf("pointer %p   %p\n", atoms->atomname, atoms2->atomname);
    printf("pointer %p   %p\n", atoms->atomname[0], atoms2->atomname[0]);
    printf("pointer %p   %p\n", *(atoms->atomname[1]), *(atoms2->atomname[1]));
    printf("pointer %p   %p\n", atoms->atomname[1][0], atoms2->atomname[1][0]);
    printf("pointer %p\n", atoms2->atomname[0]);
    printf("pointer %p\n", *(atoms2->atomname[2]));




/*
    for (j=0; j<10; j++){
   printf("1name %s\n",*(atoms->atomname[j]));
   printf("1name %c %c %c %c\n",(atoms->atomname[j])[0][0], (atoms->atomname[j])[0][1],(atoms->atomname[j])[0][2],(atoms->atomname[j])[0][3]);
   printf("1name %c %c %c %c\n",(*(atoms->atomname[j]))[0], (*(atoms->atomname[j]))[1],(*(atoms->atomname[j]))[2],(*(atoms->atomname[j]))[3]);
   }
 */
    rescounter = -1;
    for (i = 0; i < atoms->nres; i++)
    {
        rescounter = rescounter+2;
        for (j = 0; j < atoms->nr; j++)
        {

//printf("resind %d\n",atoms->atom[j].resind);
//printf("1name %c %c %c %c\n",(atoms->atomname[j])[0][0], (atoms->atomname[j])[0][1],(atoms->atomname[j])[0][2],(atoms->atomname[j])[0][3]);
//printf("pointer %p   %p\n",*(atoms->atomname[j]),*(atoms2->atomname[j]));
//printf("1name %c %c %c %c\n",(*(atoms->atomname[j]))[0], (*(atoms->atomname[j]))[1],(*(atoms->atomname[j]))[2],(*(atoms->atomname[j]))[3]);
//printf("1name %s\n",*(atoms->atomname[j]));
            if (atoms->resinfo[atoms->atom[j].resind].nr-1 == i) /* atom j is in residue i. resind is index in resinfo starting from 0 .nr is number of residue starting from 1 */
            {
                if (!strcmp((*(atoms->atomname[j])), "C") ||
                    !strcmp((*(atoms->atomname[j])), "O")) // C and O
                {
                    (*(atoms2->atomname[j])) = (!strcmp((*(atoms->atomname[j])), "C")) ? "CQMB" : "OQM";
                    atoms2->atom[j].resind   = rescounter+1;
                }
                else if (!strcmp((*(atoms->atomname[j])), "CA") ||
                         !strcmp((*(atoms->atomname[j])), "HA"))
                {
                    atoms2->atom[j].resind   = rescounter-1;
                    (*(atoms2->atomname[j])) = (!strcmp((*(atoms->atomname[j])), "C")) ? "CQMB" : "HQM";
                    if (!strcmp((*(atoms->resinfo[atoms->atom[j].resind].name)), "PRO"))
                    {
                        atoms2->atom[j].resind = rescounter+1;
                    }                                                                                                  //merge res 0 1 2 for prolin
                }
                else if (!strcmp((*(atoms->atomname[j])), "N") ||
                         !strcmp((*(atoms->atomname[j])), "H"))
                {
                    (*(atoms2->atomname[j])) = (!strcmp((*(atoms->atomname[j])), "N")) ? "NQMB" : "HQM";
                    atoms2->atom[j].resind   = rescounter-1;
                }
                else
                {
                    printf("sidechain\n");
                    if (!strcmp((*(atoms->atomname[j])), "CG"))
                    {
                        (*(atoms2->atomname[j])) = "CQMB";
                    }
                    else
                    {
//(*(atoms2->atomname[j]))=" ";
                        strncpy((*(atoms2->atomname[j])), (*(atoms->atomname[j])), 1);
                        strcat((*(atoms2->atomname[j])), "QM");
                        // (atoms2->atomname[j])[0][1]='Q'; //crop string. first letter is element type
                        //  (atoms2->atomname[j])[0][2]='M'; //crop string. first letter is element type
                        //  (atoms2->atomname[j])[0][3]='\0'; //crop string. first letter is element type
                        printf("name 2 xname %s  %s\n", (*(atoms->atomname[j])), (*(atoms2->atomname[j])));
                        //strcat((*(atoms2->atomname[j])), "QM");
                    }
                    atoms2->atom[j].resind = rescounter-1;
                    if (!strcmp((*(atoms->resinfo[atoms->atom[j].resind].name)), "PRO"))
                    {
                        atoms2->atom[j].resind = rescounter-3;
                    }                                                                                                  //merge res 0 1 2 for prolin
                }
            }
        }
    }
    atoms2->nres = rescounter;

    for (j = 0; j < atoms->nr; j++)
    {
        printf("final name %s %d\n", *(atoms2->atomname[j]), atoms2->atom[j].resind);
//printf("1name %c %c %c %c\n",(atoms->atomname[j])[0][0], (atoms->atomname[j])[0][1],(atoms->atomname[j])[0][2],(atoms->atomname[j])[0][3]);
//printf("1name %c %c %c %c\n",(*(atoms->atomname[j]))[0], (*(atoms->atomname[j]))[1],(*(atoms->atomname[j]))[2],(*(atoms->atomname[j]))[3]);
    }
    atoms = copy_t_atoms(atoms2);
    for (j = 0; j < atoms->nr; j++)
    {
        printf("final name %s %d\n", *(atoms->atomname[j]), atoms->atom[j].resind);
    }


    //write_sto_conf("preprocessor.gro", ".gro file to check residues",  atoms,  state_global->x[], NULL, NULL, NULL);
    FILE *f_ct_preprocessor = NULL;
    f_ct_preprocessor = fopen("preprocessor.gro", "w");
    fprintf(f_ct_preprocessor, "gro file with new residues\n %d\n", atoms->nr);
    for (j = 0; j < atoms->nr; j++)
    {
        fprintf(f_ct_preprocessor, "%5d%-5.5s%5.5s%5d %lf %lf %lf\n", atoms->atom[j].resind%100000, "NEW", *(atoms->atomname[j]), (j+1)%100000, state_global->x[j][XX], state_global->x[j][YY], state_global->x[j][ZZ]);
    }
    fprintf(f_ct_preprocessor, "10 10 10"); //pseudo-box
    //write_hconf_box(f_ct_preprocessor, -1, state_global->box);
    fclose(f_ct_preprocessor);

    return atoms;
}

// lapack routines
static long dsyev(int int_n, double *a, double *w, double *work, int int_lwork)
{
    extern void dsyev_(char *, char *, long *, double *, long *, double *, double *, long *, long *);
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
static long dgesv(int int_n, double *a, double *b, long ipiv[DIIS_MAX_PREV_VECTORS + 1])
{
    extern void dgesv_(long *, long *, double *, long *, long *, double *, long *, long *);
    long n, nrhs, lda, ldb, info;
    n    = lda = ldb = (long) int_n;
    nrhs = 1L;

    dgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);

    return info;
}


double calc_tda(charge_transfer_t *ct)
{
    extern void dpotrf_(char *, long *, double *, long *, long *);
    extern void dpotri_(char *, long *, double *, long *, long *);
    extern void dgetrf_(long *, long *, double *, long *, long *, long *);
    extern void dgetri_(long *, double *, long *, long *, double *, long *, long *);
    long   nb, info, lwork = 1000, ipiv[1000];
    int    i, j;
    char   itype = 'U';
    double etun, etun_old, tda;
    double vdb[1000], vba[1000], temp[1000]; //should be allocatable, no. of bridge MOs
    double gb[1000*1000], work[1000];
    double ev_Heff[2], Heff[4];              //2x2
    int    niter, maxiter = 30;

    double Etol = 1.e-7;

    printf("Performing TDA Calculation\n");

    // defining number of bridge units.
    // these are all orbitals besides donor and acceptor HOMO (multiple FOs on bridge molecules possible)
    nb = (long) ct->dim-2;


    // define tunneling energy
    etun = 0.5 * (ct->hamiltonian[0][0] + ct->hamiltonian[ct->dim-1][ct->dim-1]);

    // define donor-bridge coupling
    // define bridge-acceptor coupling
    for (i = 0; i < nb; i++)
    {
        vdb[i] = ct->hamiltonian[0][i+1];
        vba[i] = ct->hamiltonian[ct->dim-1][i+1];
    }


    //start self consistent tuning of Etun
    for (niter = 0; niter < maxiter; niter++)
    {

        // build bridge green function
        for (i = 0; i < nb; i++)
        {
            for (j = 0; j < nb; j++)
            {
                if (i == j)
                {
                    gb[i+j*nb] = etun - ct->hamiltonian[i+1][j+1];
                }
                else
                {
                    gb[i+j*nb] = -1.0 * ct->hamiltonian[i+1][j+1];
                }
            }
        }

        // output
        /*
           for (i=0; i< ct->dim; i++){
           for (j=0; j< ct->dim; j++)
            printf("%lf ", ct->hamiltonian[i][j]*HARTREE_TO_EV);
           printf("\n");
           }
           printf("Etun %lf\n", etun*HARTREE_TO_EV);
           for (i=0; i< nb; i++){
           printf("Vba %lf\n", vba[i]*HARTREE_TO_EV);
           printf("Vdb %lf\n", vdb[i]*HARTREE_TO_EV);
           }

           printf("Gb^-1[eV] =\n");
           for(i=0;i<nb;i++){
           for(j=0;j<nb;j++)
            printf("%lf ", gb[i+j*nb]*HARTREE_TO_EV);
           printf("\n");
           }
         */


        /*
           // inversion of positive definite matrix
           dpotrf_(&itype,&nb,gb,&nb,&info);
           if(info!=0){
           printf("ERROR in calc_tda (factorization)\n");
           exit(-1);
           }
           dpotri_(&itype,&nb,gb,&nb,&info);
           if(info!=0){
           printf("ERROR in calc_tda (inversion)\n");
           exit(-1);
           }
           //  */
        ///*
        //alternative: inversion of general matrix. however, indefinite or positive definite matrix indicates that bridge levels are not correct (compared to D/A)
        dgetrf_(&nb, &nb, gb, &nb, ipiv, &info);
        if ((int) info != 0)
        {
            printf("ERROR in calc_tda (dgetrf) info:%d \n", (int) info);
            exit(-1);
        }
        dgetri_(&nb, gb, &nb, ipiv, work, &lwork, &info);
        if ((int) info != 0)
        {
            printf("ERROR in calc_tda (dgetri) info:%d\n", (int) info);
            exit(-1);
        }
        //*/
        /*
           printf("Gb[eV^-1] =\n");
           for(i=0;i<nb;i++){
           for(j=0;j<nb;j++)
            printf("%lf ", gb[i+j*nb]/HARTREE_TO_EV);
           printf("\n");
           }
         */

        // construct effective two state hamiltonian
        tda       = 0.0;
        Heff[0]   = ct->hamiltonian[0][0];
        Heff[1+2] = ct->hamiltonian[ct->dim-1][ct->dim-1];
        Heff[1]   = ct->hamiltonian[0][ct->dim-1];
        Heff[1+1] = ct->hamiltonian[ct->dim-1][0];
        for (i = 0; i < nb; i++)
        {
            for (j = 0; j < nb; j++)
            {
                Heff[0]   += vdb[i]*gb[i+j*nb]*vdb[j];
                Heff[1+2] += vba[i]*gb[i+j*nb]*vba[j];
                tda       += vdb[i]*gb[i+j*nb]*vba[j];
            }
        }
        Heff[1]   += tda;
        Heff[1+1] += tda;

        // diagonalize effective Hamiltonian
        dsyev(2, Heff, ev_Heff, work, 3*2);

        // get new tunneling energy
        etun = 0.5 * (ev_Heff[0]+ev_Heff[1]);


        printf("niter %d: Direct coupling (HDA)[eV] = %f   Bridge-mediated coupling (TDA)[eV] = %f   Etun[eV] = %f \n",
               niter, ct->hamiltonian[0][ct->dim-1]*HARTREE_TO_EV, tda*HARTREE_TO_EV, etun*HARTREE_TO_EV);

        if (fabs(etun-etun_old) < Etol)
        {
            printf("TDA calculation converged after %d iterations. \n", niter);
            break;
        }
        else if (niter == maxiter-1)
        {
            printf("TDA calculation did not converge after %d iterations. \n", niter);
        }

        etun_old = etun;
    } // end tuning Etun

    return tda;
}


void ct_assemble_hamiltonian(charge_transfer_t *ct, dftb_t *dftb)
{
    int i, j, k, l, m, n, t;
//  mean coupling for a and b direction in SAM.

    //save old hamilton matrices (they range from t=0 as newest martrix up to n_avg_ham steps into the past)
    //start with filled array from first step for meaningful averaging

    if (ct->jobtype == cteBCJFSSH)
    {
       for (i=0; i<ct->dim; i++)
           ct->occupation[i] = SQR(ct->tfs_diab[i]) + SQR(ct->tfs_diab[i+ct->dim]);
    }
    else if (ct->jobtype == cteSCCDYNAMIC)
    {
       for (i=0; i<ct->dim; i++)
           ct->occupation[i] = SQR(ct->wf[i]) + SQR(ct->wf[i+ct->dim]);
    }

    if (ct->first_step)
    {
        for (t = 0; t < ct->n_avg_ham; t++)
        {
            for (i = 0; i < ct->dim; i++)
            {
                ct->hamiltonian_history[i][i][t]  = (ct->is_hole_transfer == 1) ? -dftb->phase2.THamilOrtho[i][i] : dftb->phase2.THamilOrtho[i][i];
                ct->hamiltonian_history[i][i][t] += ct->fo_shift[i];
                for (j = i+1; j < ct->dim; j++)
                {
                    ct->hamiltonian_history[i][j][t] = ct->hamiltonian_history[j][i][t] = (ct->is_hole_transfer == 1) ? -dftb->phase2.THamilOrtho[i][j] * ct->offdiag_scaling :  dftb->phase2.THamilOrtho[i][j] * ct->offdiag_scaling;//  1.0 or 1.540 scaling factor according to J. Chem. Phys. 140, 104105 (2014) // this is the right way to do it. without changing off-diagonal elements the eigenvalues of the CG Ham are not necessarily -eigenvalue of electron transfer.
                }
            }
        }
    }
    else
    {
        for (t = ct->n_avg_ham-1; t > 0; t--)
        {
            for (i = 0; i < ct->dim; i++)
            {
                for (j = 0; j < ct->dim; j++)
                {
                    ct->hamiltonian_history[i][j][t] = ct->hamiltonian_history[i][j][t-1];
                }
            }
        }

        for (i = 0; i < ct->dim; i++)
        {
            ct->hamiltonian_history[i][i][0]  = (ct->is_hole_transfer == 1) ? -dftb->phase2.THamilOrtho[i][i] : dftb->phase2.THamilOrtho[i][i];
            ct->hamiltonian_history[i][i][0] += ct->fo_shift[i];
            for (j = i+1; j < ct->dim; j++)
            {
                ct->hamiltonian_history[i][j][0] = ct->hamiltonian_history[j][i][0] = (ct->is_hole_transfer == 1) ? -dftb->phase2.THamilOrtho[i][j] * ct->offdiag_scaling :  dftb->phase2.THamilOrtho[i][j] * ct->offdiag_scaling;//  1.0 or 1.540 scaling factor according to J. Chem. Phys. 140, 104105 (2014) // this is the right way to do it. without changing off-diagonal elements the eigenvalues of the CG Ham are not necessarily -eigenvalue of electron transfer.
            }
        }
    }

    //build averaged hamiltonian
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            ct->hamiltonian_dftb[i][j] = 0.0;
            for (t = 0; t < ct->n_avg_ham; t++)
            {
                ct->hamiltonian_dftb[i][j] += ct->hamiltonian_history[i][j][t];
            }
            ct->hamiltonian_dftb[i][j] /= ct->n_avg_ham;
            if (j != i)
            {
                ct->hamiltonian_dftb[i][j] = ct->hamiltonian_history[i][j][0];//test: only average diagonal elements. fast oscillations of off-diags can also be caused by slow (classical) nuclear oscillation due to complex nodal structure of the orbitals.
            }
            else if (j == i && ct->lambda != 0.0)
            {
                ct->hamiltonian_dftb[i][i] -= ct->occupation[i]*ct->lambda;
            }
        }
    }

    if (ct->freeze_coupling != 0.0)
    {
       for (i = 0; i < ct->dim - 1; i++)
       {
           for (j = i + 1; j < ct->dim; j++)
           {
               if (j == i + 1)
               {
                  ct->hamiltonian_dftb[i][j] = ct->freeze_coupling;
               }
               else
               {
                  ct->hamiltonian_dftb[i][j] = 0.0;
               }

               ct->hamiltonian_dftb[j][i] = ct->hamiltonian_dftb[i][j];
           }
       }
    }

    return;
}

void broyden(int niter, double alpha, int jtop, double *vecin, double *vecout, dftb_broyden_t *arrays);

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

int do_adiabatic(charge_transfer_t *ct, ct_diis_t *diis, FILE *f)
//int do_adiabatic(charge_transfer_t *ct, dftb_broyden_t *broyd, FILE *f)
{
    int    i, j, k, step;
    double old_energy, energy, *ham, almix, fermi_energy = 0, fermi_upper, fermi_lower;

    ham   = ct->hamiltonian_adiab;
    almix = SIMPLE_ALMIX;

    diis->n_prev_vector = 1;
    diis->indx          = 0;

    // calculate the charges from the wave function
    for (i = 0; i < ct->dim; i++)
    {
        ct->q_act[i] = SQR(ct->wf[i]);
    }

    // copy as the first vector

    energy = 1.e10;

    // SCC cycle
    for (step = 0;; step++)
    {

        // too many iterations?
        if (step >= MAXITER_DIIS)
        {
            // copy the solution to the wave function
            for (i = 0; i < ct->dim; i++)
            {
                ct->wf[i] = sqrt(ct->q_act[i]);
            }
            return 1;
        }

        printf("step %5d:", step);
        //fprintf(f, "step %5d:", step);

        // save and print old charges
        printf(" old q =");
        //fprintf(f, " old q =");
        for (i = 0; i < ct->dim; i++)
        {
            //fprintf(f, "%12.9f", ct->q_act[i]);
            printf("%12.9f", ct->q_act[i]);
        }

        //  FOR BROYDEN:
        //  ct->q_old[i] = ct->q_act[i];

        // construct the self-consistent Hamiltonian
        for (i = 0; i < ct->dim; i++)
        {
            ham[i + i * ct->dim] = ct->hamiltonian[i][i];
            for (j = 0; j < ct->dim; j++)
            {
                ham[i + i * ct->dim] += ct->hubbard[i][j] * ct->q_act[j];
                if (i != j)
                {
                    ham[i + j * ct->dim] = ct->hamiltonian[i][j];
                }
            }
        }

        // diagonalize this Hamiltonian
        dsyev(ct->dim, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);

        // copy previous charge and the difference to the arrays
        for (i = 0; i < ct->dim; i++)
        {
            diis->prev_q_input[diis->indx][i] = ct->q_act[i];
            diis->prev_q_diff[diis->indx][i]  = -ct->q_act[i]; // here, the new occupation will be added
        }
/*
    // instead of Fermi:
    for (i=0; i<ct->sites; i++)
      diis->prev_q_diff[diis->indx][i] += SQR(ham[i]);
 */
        // FERMI DISTRIBUTION (electronic temperature)
        // calculate the "Fermi energy", store the Fermi coefficients in diis->fermi_coef
        fermi_lower = ct->ev_adiab[0] - 0.1;
        fermi_upper = ct->ev_adiab[ct->dim - 1] + 0.1;
        while (fermi_upper - fermi_lower > FERMI_CONVERG)
        {
            fermi_energy = 0.5 * (fermi_lower + fermi_upper);
            if (fermi_coef_sum(fermi_energy, ct->dim, ct->ev_adiab, diis->fermi_coef, ct->fermi_kt) > 1.0)
            {
                fermi_upper = 0.5 * (fermi_lower + fermi_upper);
            }
            else
            {
                fermi_lower = 0.5 * (fermi_lower + fermi_upper);
            }
        }
        fprintf(f, " fermi E = %9.6f, fermi coefs:", fermi_energy);
        for (i = 0; i < ct->dim; i++)
        {
            fprintf(f, " %8.5f", diis->fermi_coef[i]);
        }
        // mix the eigenvectors, with the appropriate coefficients
        for (i = 0; i < ct->dim; i++)     // i runs over the sites
        {
            for (j = 0; j < ct->dim; j++) // j runs over the individual eigenvectors
            {
                diis->prev_q_diff[diis->indx][i] += diis->fermi_coef[j] * SQR(ham[j * ct->dim + i]);
            }
        }

        // calculate the energy and check convergence
        old_energy = energy;
        energy     = 0.e0;
        /* old - calculate energy by orbitals
           for (i=0; i<ct->sites; i++)
           for (j=0; j<ct->sites; j++)
            energy += ham[i] * ham[j] * ct->hamiltonian[i][j]
         + 0.5 * ct->q_act[i] * ct->q_act[j] * (ct->sic ? (SIC_COEFFICIENT * ct->hubbard[i][j]) : (ct->hubbard[i][j]));
         */
        for (i = 0; i < ct->dim; i++)
        {
            energy += diis->fermi_coef[i] * ct->ev_adiab[i];
        }
        printf(" energy = %12.9f\n", energy);
        //fprintf(f, " energy = %12.9f\n", energy);
        if (fabs(energy - old_energy) < BROYDEN_SCFTOL)
        {
            //fprintf(f, "Converged!\n");
            printf("Converged!\n");
            // copy the solution to the wave function
            for (i = 0; i < ct->dim; i++)
            {
                //ct->wf[i] = sqrt(ct->q_act[i]);
                ct->wf[i] = sqrt(diis->prev_q_input[diis->indx][i] + diis->prev_q_diff[diis->indx][i]);
            }
            //ct->wf[i] = ham[i];
            break;
        }

        // DIIS - Direct Inversion in the Iterative Subspace
        // we need in memory:
        // p        -- matrix (# of bases) x (# of vectors)
        // delta p  -- matrix (# of bases) x (# of vectors)
        // b        -- matrix (# of vectors + 1) x (# of vectors + 1)

        // produce a new vector "ct->q_old", either by simple mixing (in first iterations) or by DIIS (later)
        if (diis->n_prev_vector < DIIS_MAX_PREV_VECTORS)
        {
            // SIMPLE MIXING in a couple of first steps
            for (i = 0; i < ct->dim; i++)
            {
                ct->q_act[i] = diis->prev_q_input[diis->indx][i] + DIIS_INIT_MIXING_PARAM * diis->prev_q_diff[diis->indx][i];
            }
            //ct->q_act[i] = (1.0 - DIIS_INIT_MIXING_PARAM) * ct->q_act[i] + DIIS_INIT_MIXING_PARAM * diis->q_inp_result[i];
            diis->n_prev_vector++;
        }
        else
        {
            // real DIIS here
            // what do the vectors "prev_q_diff" look like?
            //printf("\n DIFF VECTORS\n");
            //for (i=0; i<DIIS_MAX_PREV_VECTORS; i++) {
            //  for (j=0; j<ct->sites; j++)
            //  printf("%14.10f", diis->prev_q_diff[i][j]);
            //printf("\n");
            //}
            // construct the matrix of overlaps
            for (i = 0; i < DIIS_MAX_PREV_VECTORS; i++)
            {
                for (j = 0; j < DIIS_MAX_PREV_VECTORS; j++)
                {
                    diis->aa[i + j * (DIIS_MAX_PREV_VECTORS + 1)] = 0.e0;
                    for (k = 0; k < ct->dim; k++)
                    {
                        diis->aa[i + j * (DIIS_MAX_PREV_VECTORS + 1)] += diis->prev_q_diff[i][k] * diis->prev_q_diff[j][k];
                    }
                    //printf("aa[%d][%d] = %20.15f\n", i, j, diis->aa[i + j * (DIIS_MAX_PREV_VECTORS + 1)]);
                }
            }
            for (i = 0; i < DIIS_MAX_PREV_VECTORS; i++)
            {
                diis->aa[i + DIIS_MAX_PREV_VECTORS * (DIIS_MAX_PREV_VECTORS + 1)] = diis->aa[DIIS_MAX_PREV_VECTORS + i * (DIIS_MAX_PREV_VECTORS + 1)] = -1.e0;
            }
            diis->aa[DIIS_MAX_PREV_VECTORS + DIIS_MAX_PREV_VECTORS * (DIIS_MAX_PREV_VECTORS + 1)] = 0.e0;
            //printf("\n MATRIX AA:\n");
            //for (i=0; i<=DIIS_MAX_PREV_VECTORS; i++) {
            //  for (j=0; j<=DIIS_MAX_PREV_VECTORS; j++)
            //  printf("%14.10f", diis->aa[i + j * (DIIS_MAX_PREV_VECTORS + 1)]);
            //printf("\n");
            //}
            // construct the vector of right hand sides
            for (i = 0; i < DIIS_MAX_PREV_VECTORS; i++)
            {
                diis->bb[i] = 0.e0;
            }
            diis->bb[DIIS_MAX_PREV_VECTORS] = -1.e0;
            //printf(" MATRIX BB:\n");
            //for (i=0; i<=DIIS_MAX_PREV_VECTORS; i++)
            //  printf("%14.10f", diis->bb[i]);
            //printf("\n");
            // solve it with lapack
            //printf(" DGESV: %ld", dgesv(DIIS_MAX_PREV_VECTORS + 1, diis->aa, diis->bb, diis->ipiv));
            dgesv(DIIS_MAX_PREV_VECTORS + 1, diis->aa, diis->bb, diis->ipiv);
            // the column of solution is now in diis->bb -- print it!
            //printf(" DIIS coeffs");
            //for (i=0; i<=DIIS_MAX_PREV_VECTORS; i++)
            //  printf("%9.5f", diis->bb[i]);
            for (j = 0; j < ct->dim; j++)
            {
                ct->q_act[j] = 0.e0;
                for (i = 0; i < DIIS_MAX_PREV_VECTORS; i++)
                {
                    ct->q_act[j] += diis->bb[i] * diis->prev_q_input[i][j];
                }
            }
        }

        // move the counter in a loop
        diis->indx = (diis->indx + 1) % DIIS_MAX_PREV_VECTORS;

/*
    // FERMI DISTRIBUTION + BROYDEN MIXING
    // mix all eigenvectors with appropriate coefficients

    // calculate the "Fermi energy", store the Fermi coefficients in broyd->df
    fermi_lower = ct->ev_adiab[0] - 0.1;
    fermi_upper = ct->ev_adiab[ct->sites - 1] + 0.1;
    while (fermi_upper - fermi_lower > FERMI_CONVERG) {
      fermi_energy = 0.5 * (fermi_lower + fermi_upper);
      if (fermi_coef_sum(fermi_energy, ct->sites, ct->ev_adiab, broyd->df, ct->fermi_kt) > 1.0)
        fermi_upper = 0.5 * (fermi_lower + fermi_upper);
      else
        fermi_lower = 0.5 * (fermi_lower + fermi_upper);
    }
    //fermi_coef2 = 1 / (exp((ct->ev_adiab[1] - ct->ev_adiab[0]) / FERMI_KT) + 1);
    //printf(" fermi coef2 = %7.4f", fermi_coef2);
    fprintf(f, " fermi E = %9.6f, fermi coefs:", fermi_energy);
    for (i=0; i<ct->sites; i++)
      fprintf(f, " %8.5f", broyd->df[i]);

    // mix the eigenvectors, with the appropriate coefficients
    for (i=0; i<ct->sites; i++) { // i runs over the sites
      ct->q_act[i] = 0.e0;
      for (j=0; j<ct->sites; j++) // j runs over the individual eigenvectors
        ct->q_act[i] += broyd->df[j] * SQR(ham[j * ct->sites + i]);
    }

    // mix it - Broyden
    broyden(step, BROYDEN_ALMIX, ct->sites, ct->q_old, ct->q_act, *broyd);
    for (i=0; i<ct->sites; i++)
      ct->q_act[i] = ct->q_old[i];
 */

/*
    // SIMPLE MIXING
    //
    // calculate the charges
    for (i=0; i<ct->sites; i++)
      ct->q_act[i] = SQR(ham[i]);
    // calculate oldness - the content of old vector in the new one
    oldness = 0.e0;
    for (i=0; i<ct->sites; i++)
      oldness += ct->q_act[i] * ct->q_old[i];
    // the orthogonal vector is now: new_ortho = new - oldness * old
    // calculate the new vector as: old + almix(t) * new_ortho
    for (i=0; i<ct->sites; i++)
      ct->q_act[i] = (1.0 - almix) * ct->q_old[i] + almix * (ct->q_act[i] - oldness * ct->q_old[i]);
    // attenuate the mixing coefficient
    almix *= ALMIX_ATTENUATOR;
 */

    }

    return 0;
}

int do_wf_decomp_evec(charge_transfer_t *ct, ct_diis_t *diis, FILE *f)
{
    int    i, j, k, step;
    double old_energy, energy, *ham, dotproduct_squared;

    ham                 = ct->hamiltonian_adiab; // why hamiltonian_adiab aren't we in diab basis dholub July 2018
    diis->n_prev_vector = 1;
    diis->indx          = 0;

    // calculate the charges from the wave function
    for (i = 0; i < ct->dim; i++)
    {
        ct->q_act[i] = SQR(ct->wf[i]);
    }

    // SCC cycle
    for (step = 0, energy = 1.e10;; step++)
    {

        if (step >= MAXITER_DIIS)
        {
            /* DO HERE SOMETHING, BUT DO NOT MODIFY ct->wf !
               for (i=0; i<ct->sites; i++)
               ct->wf[i] = sqrt(ct->q_act[i]);
             */
            fprintf(f, "diagonalization did not converge! ");
            return 1;
        }

        //printf("step %5d:", step);

        //printf(" old q =");
        //for (i=0; i<ct->sites; i++)
        //  printf("%9.5f", ct->q_act[i]);

        // construct the self-consistent Hamiltonian
        for (i = 0; i < ct->dim; i++)
        {
            ham[i + i * ct->dim] = ct->hamiltonian[i][i];
            for (j = 0; j < ct->dim; j++)
            {
                ham[i + i * ct->dim] += ct->hubbard[i][j] * ct->q_act[j];
                if (i != j)
                {
                    ham[i + j * ct->dim] = ct->hamiltonian[i][j];
                }
            }
        }

        // diagonalize this Hamiltonian
        dsyev(ct->dim, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);

        // copy previous charge and the difference to the arrays
        for (i = 0; i < ct->dim; i++)
        {
            diis->prev_q_input[diis->indx][i] = ct->q_act[i];
            diis->prev_q_diff[diis->indx][i]  = -ct->q_act[i]; // here, the new occupation will be added
        }
        // FERMI DISTRIBUTION (electronic temperature)
        // calculate the "Fermi energy", store the Fermi coefficients in diis->fermi_coef
        /* SWITCH IT OFF FOR DECOMPOSITION !
           fermi_lower = ct->ev_adiab[0] - 0.1;
           fermi_upper = ct->ev_adiab[ct->sites - 1] + 0.1;
           while (fermi_upper - fermi_lower > FERMI_CONVERG) {
           fermi_energy = 0.5 * (fermi_lower + fermi_upper);
           if (fermi_coef_sum(fermi_energy, ct->sites, ct->ev_adiab, diis->fermi_coef, ct->fermi_kt) > 1.0)
            fermi_upper = 0.5 * (fermi_lower + fermi_upper);
           else
            fermi_lower = 0.5 * (fermi_lower + fermi_upper);
           }
           fprintf(f, " fermi E = %9.6f, fermi coefs:", fermi_energy);
           for (i=0; i<ct->sites; i++)
           fprintf(f, " %8.5f", diis->fermi_coef[i]);
           for (i=0; i<ct->sites; i++)
           for (j=0; j<ct->sites; j++)
            diis->prev_q_diff[diis->indx][i] += diis->fermi_coef[j] * SQR(ham[j * ct->sites + i]);
         */
        for (i = 0; i < ct->dim; i++)
        {
            diis->prev_q_diff[diis->indx][i] += SQR(ham[i]);
        }

        old_energy = energy;
        /* AGAIN - DO NOT MIX THE ENERGY !
           energy = 0.e0;
           for (i=0; i<ct->sites; i++)
           energy += diis->fermi_coef[i] * ct->ev_adiab[i];
         */
        energy = ct->ev_adiab[0];
        //printf(" energy = %12.9f\n", energy);
        if (fabs(energy - old_energy) < BROYDEN_SCFTOL)
        {
            //printf("Converged!\n");
            /* NOW, WE HAVE THE ADIABATIC STATES IN ham !
               for (i=0; i<ct->sites; i++)
               ct->wf[i] = sqrt(diis->prev_q_input[diis->indx][i] + diis->prev_q_diff[diis->indx][i]);
             */
            break;
        }

        if (diis->n_prev_vector < DIIS_MAX_PREV_VECTORS)
        {
            // SIMPLE MIXING in a couple of first steps
            for (i = 0; i < ct->dim; i++)
            {
                ct->q_act[i] = diis->prev_q_input[diis->indx][i] + DIIS_INIT_MIXING_PARAM * diis->prev_q_diff[diis->indx][i];
            }
            diis->n_prev_vector++;
        }
        else
        {
            // REAL DIIS here
            for (i = 0; i < DIIS_MAX_PREV_VECTORS; i++)
            {
                for (j = 0; j < DIIS_MAX_PREV_VECTORS; j++)
                {
                    diis->aa[i + j * (DIIS_MAX_PREV_VECTORS + 1)] = 0.e0;
                    for (k = 0; k < ct->dim; k++)
                    {
                        diis->aa[i + j * (DIIS_MAX_PREV_VECTORS + 1)] += diis->prev_q_diff[i][k] * diis->prev_q_diff[j][k];
                    }
                }
            }
            for (i = 0; i < DIIS_MAX_PREV_VECTORS; i++)
            {
                diis->aa[i + DIIS_MAX_PREV_VECTORS * (DIIS_MAX_PREV_VECTORS + 1)] = diis->aa[DIIS_MAX_PREV_VECTORS + i * (DIIS_MAX_PREV_VECTORS + 1)] = -1.e0;
            }
            diis->aa[DIIS_MAX_PREV_VECTORS + DIIS_MAX_PREV_VECTORS * (DIIS_MAX_PREV_VECTORS + 1)] = 0.e0;

            for (i = 0; i < DIIS_MAX_PREV_VECTORS; i++)
            {
                diis->bb[i] = 0.e0;
            }
            diis->bb[DIIS_MAX_PREV_VECTORS] = -1.e0;

            dgesv(DIIS_MAX_PREV_VECTORS + 1, diis->aa, diis->bb, diis->ipiv);

            for (j = 0; j < ct->dim; j++)
            {
                ct->q_act[j] = 0.e0;
                for (i = 0; i < DIIS_MAX_PREV_VECTORS; i++)
                {
                    ct->q_act[j] += diis->bb[i] * diis->prev_q_input[i][j];
                }
            }
        }
        diis->indx = (diis->indx + 1) % DIIS_MAX_PREV_VECTORS;
    }

    // print energies of adiabatic states
    //printf("Adiab E: ");
    for (i = 0; i < ct->dim; i++)
    {
        //fprintf(f, "%10.7f", ct->ev_adiab[i]);
        fprintf(f, "%10.7f", ct->ev_adiab[i] - ct->ev_adiab[0]);
    }

    // calculate the expansion coefficients
    fprintf(f, "   WF exp'd in adiab states: ");
    /*
       for (i=0; i<ct->sites; i++) {
       // real part of a_i
       dotproduct_re = 0.e0;
       for (j=0; j<ct->sites; j++) dotproduct_re += ham[i * ct->sites + j] * ct->wf[j];
       // imaginary part of a_i
       dotproduct_im = 0.e0;
       for (j=0; j<ct->sites; j++) dotproduct_im += ham[i * ct->sites + j] * ct->wf[j + ct->sites];
       fprintf(f, "%9.6f", SQR(dotproduct_re)+SQR(dotproduct_im));
       }
     */
    for (i = 0; i < ct->dim; i++)
    {
        // calculate directly the square of a_i
        dotproduct_squared = 0.e0;
        for (j = 0; j < ct->dim; j++)
        {
            for (k = 0; k < ct->dim; k++)
            {
                dotproduct_squared += ham[i * ct->dim + j] * ham[i * ct->dim + k] * (ct->wf[j] * ct->wf[k] + ct->wf[j + ct->dim] * ct->wf[k + ct->dim]);
            }
        }
        fprintf(f, "%9.6f", dotproduct_squared);
    }
    fprintf(f, "\n");

    return 0;
}

int do_surface_hopping(charge_transfer_t *ct, ct_diis_t *diis, FILE *f)
{
    /* this should be done by the modification of do_adiabatic
     * followed by the testing of scalar product of the new wavefunction - ground and 1st excited state - with the old one
     */
    /* ct->surface = 0 (ground state) or 1 (1st excited state) */
    int    i, j, k, step;
    double old_energy, energy, *ham, s01, s10, nonadiab_coupling;

    ham                 = ct->hamiltonian_adiab;
    diis->n_prev_vector = 1;
    diis->indx          = 0;

    // calculate the charges from the wave function
    for (i = 0; i < ct->sites; i++)
    {
        ct->q_act[i] = SQR(ct->wf[i]);
    }

    // SCC cycle
    for (step = 0, energy = 1.e10;; step++)
    {

        if (step >= MAXITER_DIIS)
        {
            /* DIIS NOT CONVERGED! */
            //printf("DIIS not converged !!!      ");
            //fprintf(f, "DIIS not converged !!!      ");
            fprintf(f, "DIIS   0 ");
            fprintf(f, "             ");
            for (i = 0; i < ct->sites; i++)
            {
                fprintf(f, "        ");
            }
            fprintf(f, "                           DIIS not converged so    hopping not tested... ");
            fprintf(f, " %d", ct->surface);
            for (i = 0; i < ct->sites; i++)
            {
                fprintf(f, "%8.4f", SQR(ct->wf[i]));
            }
            fprintf(f, "\n");

            return 1;
        }

        //printf("step %5d:", step);

        //printf(" old q =");
        //for (i=0; i<ct->sites; i++)
        //printf("%9.5f", ct->q_act[i]);

        // construct the self-consistent Hamiltonian
        for (i = 0; i < ct->sites; i++)
        {
            ham[i + i * ct->sites] = ct->hamiltonian[i][i];
            for (j = 0; j < ct->sites; j++)
            {
                ham[i + i * ct->sites] += ct->hubbard[i][j] * ct->q_act[j];
                if (i != j)
                {
                    ham[i + j * ct->sites] = ct->hamiltonian[i][j];
                }
            }
        }

        /*
           // orbital shift:
           // 1. preliminary diagonalization
           dsyev(ct->sites, ham, ct->ev_adiab, ct->work_adiab, 3*ct->sites);
           // 2. correction to the elements of the Hamiltonian
           //    Delta H = C * delta-Epsilon * C^T
           for (i=0; i<ct->sites; i++)
           ct->work_adiab[i] = ham[i];
           for (i=0; i<ct->sites; i++) {
           ham[i + i * ct->sites] = ct->hamiltonian[i][i] - SFH_ORBITAL_SHIFT * ct->work_adiab[i] * ct->work_adiab[i];
           for (j=0; j<ct->sites; j++) {
            ham[i + i * ct->sites] += ct->hubbard[i][j] * ct->q_act[j];
            if (i!=j) ham[i + j * ct->sites] = ct->hamiltonian[i][j] - SFH_ORBITAL_SHIFT * ct->work_adiab[i] * ct->work_adiab[j];
           }
           }
           // done with orbital shift
         */

        // diagonalize this Hamiltonian
        dsyev(ct->sites, ham, ct->ev_adiab, ct->work_adiab, 3*ct->sites);

        // copy previous charge and the difference to the arrays
        for (i = 0; i < ct->sites; i++)
        {
            diis->prev_q_input[diis->indx][i] = ct->q_act[i];
            diis->prev_q_diff[diis->indx][i]  = SQR(ham[i]) - ct->q_act[i];
        }

        old_energy = energy;
        energy     = ct->ev_adiab[0];
        //printf(" energy = %12.9f\n", energy);
        if (fabs(energy - old_energy) < BROYDEN_SCFTOL)
        {
            //printf("DIIS Converged in %3d steps ", step);
            //fprintf(f, "DIIS Converged in %3d steps ", step);
            fprintf(f, "DIIS %3d ", step);
            break;
        }

        if (diis->n_prev_vector < DIIS_MAX_PREV_VECTORS)
        {
            // SIMPLE MIXING in a couple of first steps
            for (i = 0; i < ct->sites; i++)
            {
                ct->q_act[i] = diis->prev_q_input[diis->indx][i] + DIIS_INIT_MIXING_PARAM * diis->prev_q_diff[diis->indx][i];
            }
            diis->n_prev_vector++;
        }
        else
        {
            // REAL DIIS here
            for (i = 0; i < DIIS_MAX_PREV_VECTORS; i++)
            {
                for (j = 0; j < DIIS_MAX_PREV_VECTORS; j++)
                {
                    diis->aa[i + j * (DIIS_MAX_PREV_VECTORS + 1)] = 0.e0;
                    for (k = 0; k < ct->sites; k++)
                    {
                        diis->aa[i + j * (DIIS_MAX_PREV_VECTORS + 1)] += diis->prev_q_diff[i][k] * diis->prev_q_diff[j][k];
                    }
                }
            }
            for (i = 0; i < DIIS_MAX_PREV_VECTORS; i++)
            {
                diis->aa[i + DIIS_MAX_PREV_VECTORS * (DIIS_MAX_PREV_VECTORS + 1)] = diis->aa[DIIS_MAX_PREV_VECTORS + i * (DIIS_MAX_PREV_VECTORS + 1)] = -1.e0;
            }
            diis->aa[DIIS_MAX_PREV_VECTORS + DIIS_MAX_PREV_VECTORS * (DIIS_MAX_PREV_VECTORS + 1)] = 0.e0;

            for (i = 0; i < DIIS_MAX_PREV_VECTORS; i++)
            {
                diis->bb[i] = 0.e0;
            }
            diis->bb[DIIS_MAX_PREV_VECTORS] = -1.e0;

            dgesv(DIIS_MAX_PREV_VECTORS + 1, diis->aa, diis->bb, diis->ipiv);

            for (j = 0; j < ct->sites; j++)
            {
                ct->q_act[j] = 0.e0;
                for (i = 0; i < DIIS_MAX_PREV_VECTORS; i++)
                {
                    ct->q_act[j] += diis->bb[i] * diis->prev_q_input[i][j];
                }
            }
        }
        diis->indx = (diis->indx + 1) % DIIS_MAX_PREV_VECTORS;
    }

    // print energies of adiabatic states
    //printf(" Eigenvalues:");
    fprintf(f, " Eigenvalues:");
    for (i = 0; i < ct->sites; i++)
    {
        //fprintf(f, "%10.7f", ct->ev_adiab[i]);
        fprintf(f, "%8.5f", ct->ev_adiab[i]);
        //printf("%8.5f", ct->ev_adiab[i] - ct->ev_adiab[0]);
    }
    //printf("    ");
    fprintf(f, "    ");

    // check the energy gap between the ground and the first excited state
    if (ct->ev_adiab[1] - ct->ev_adiab[0] < SFH_LIMIT)
    {
        /* attempt surface hopping
         * a variation of the AS3 method in Fabiano et al., Chem. Phys. 351 (2008) 111-116
         */
        s01 = 0.e0;
        for (i = 0; i < ct->sites; i++)
        {
            s01 += (ct->surface == 0 ? ct->wf[i] : ct->wf_exc[i]) * ham[ct->sites + i];
        }
        s10 = 0.e0;
        for (i = 0; i < ct->sites; i++)
        {
            s10 += (ct->surface == 0 ? ct->wf_exc[i] : ct->wf[i]) * ham[i];
        }
        /* calculate the (quasi-)non-adiabatic coupling */
        nonadiab_coupling = 1 - pow(1 - (fabs(s01) + fabs(s10))*0.5, SFH_EXPONENT);
        /* print them out */
        //printf("wf.ground = %6.3f, wf.excited = %6.3f ", overlap_ground, overlap_excited);
        fprintf(f, "s01 = %7.4f, s10 = %7.4f, coupling = %6.3f ", s01, s10, nonadiab_coupling);

        /* check if hopping should be performed */
        if (nonadiab_coupling > rand()/(double)RAND_MAX)
        {
            fprintf(f, "hopping between states");
            if (ct->surface == 0) // hopping 0->1
            {
                for (i = 0; i < ct->sites; i++)
                {
                    ct->wf[i]     = ham[ct->sites + i];
                    ct->wf_exc[i] = ham[i];
                }
                ct->surface = 1;
            }
            else              // hopping 1->0
            {
                for (i = 0; i < ct->sites; i++)
                {
                    ct->wf[i]     = ham[i];
                    ct->wf_exc[i] = ham[ct->sites + i];
                }
                ct->surface = 0;
            }
        }
        else
        {
            fprintf(f, "no hopping now...     ");
            for (i = 0; i < ct->sites; i++)
            {
                if (ct->surface == 0) // system is remaining in state 0
                {
                    ct->wf[i]     = ham[i];
                    ct->wf_exc[i] = ham[ct->sites + i];
                }
                else            // system is remaining in state 1
                {
                    ct->wf[i]     = ham[ct->sites + i];
                    ct->wf_exc[i] = ham[i];
                }
            }
        }
    }
    else
    {
        //printf("Egap > threshold, no surf hopping ");
        fprintf(f, "                                                hopping not tested... ");
        for (i = 0; i < ct->sites; i++)
        {
            if (ct->surface == 0) // system is remaining in state 0
            {
                ct->wf[i]     = ham[i];
                ct->wf_exc[i] = ham[ct->sites + i];
            }
            else              // system is remaining in state 1
            {
                ct->wf[i]     = ham[ct->sites + i];
                ct->wf_exc[i] = ham[i];
            }
        }
    }
    fprintf(f, " %d", ct->surface);
    for (i = 0; i < ct->sites; i++)
    {
        //printf("%8.4f", SQR(ct->wf[i]));
        fprintf(f, "%8.4f", SQR(ct->wf[i]));
    }
    //printf("\n");
    fprintf(f, "\n");

    return 0;
}







double do_born_oppenheimer(charge_transfer_t *ct, dftb_broyden_t *broyd, double *fermi_coeff, FILE *f)
{
    int    i, j, k, step, ier;
    double old_energy, energy, energy1, energy2, *ham, fermi_energy, fermi_upper, fermi_lower, q_sum, s01, best_overlap;

    ham = ct->hamiltonian_adiab;

    // calculate the charges from the wave function
    for (i = 0; i < ct->dim; i++)
    {
        ct->q_act[i] = ct->q_old[i] = SQR(ct->wf[i]);
    }

    // copy as the first vector

    energy = 1.e10;
    // SCC cycle
    for (step = 0;; step++)
    {

        // too many iterations? break and indicate failure
        if (step >= MAXITER_BROYDEN)
        {
            fprintf(f, "Broyden/Fermi failed to converge!\n");
            return -1;
        }

        // construct the self-consistent Hamiltonian
        for (i = 0; i < ct->dim; i++)
        {
            ham[i + i * ct->dim] = ct->hamiltonian[i][i];
            for (j = 0; j < ct->dim; j++)
            {
                ham[i + i * ct->dim] += ct->hubbard[i][j] * ct->q_act[j];
                if (i != j)
                {
                    ham[i + j * ct->dim] = ct->hamiltonian[i][j];
                }
            }
        }

        // diagonalize this Hamiltonian
        ier = dsyev(ct->dim, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);
        if ((int) ier)
        {
            printf("\nDSYEV: ier = %d\nEXITING!\n\n", (int) ier);
            exit(-1);
        }


        // FERMI DISTRIBUTION + BROYDEN MIXING
        // mix all eigenvectors with appropriate coefficients

        // calculate the "Fermi energy", store the Fermi coefficients in fermi_coeff
        fermi_lower = ct->ev_adiab[0] - 0.1;
        fermi_upper = ct->ev_adiab[ct->dim - 1] + 0.1;
        while (fermi_upper - fermi_lower > FERMI_CONVERG)
        {
            fermi_energy = 0.5 * (fermi_lower + fermi_upper);
            if (fermi_coef_sum(fermi_energy, ct->dim, ct->ev_adiab, fermi_coeff, ct->fermi_kt) > 1.0)
            {
                fermi_upper = 0.5 * (fermi_lower + fermi_upper);
            }
            else
            {
                fermi_lower = 0.5 * (fermi_lower + fermi_upper);
            }
        }

        // mix the eigenvectors, with the appropriate coefficients
        for (i = 0; i < ct->dim; i++) // i runs over the sites
        {
            ct->q_act[i] = 0.e0;
            for (j = 0; j < ct->dim; j++) // j runs over the individual eigenvectors
            {
                ct->q_act[i] += fermi_coeff[j] * SQR(ham[j * ct->dim + i]);
            }
        }

        // calculate the energy and check convergence
        old_energy = energy;
        energy1    = energy2 = 0.0;
        // indices i and j run over the sites
        for (i = 0; i < ct->dim; i++)
        {
            for (j = 0; j < ct->dim; j++)
            {
                // TB Hamiltonian
                // index k runs over the adiabatic states (eigenvectors)
                for (k = 0; k < ct->dim; k++)
                {
                    energy1 += fermi_coeff[k] * ham[k * ct->dim + i] * ham[k * ct->dim + j] * ct->hamiltonian[i][j];
                }
                // Hubbard / gamma terms
                energy2 += 0.5 * ct->hubbard[i][j] * ct->q_act[i] * ct->q_act[j];
            }
        }
        energy = energy1 + energy2;

        if (step > 1 && fabs(energy - old_energy) < BROYDEN_SCFTOL)
        {
            break;
        }

        // mix the charges - Broyden
        broyden(step, BROYDEN_ALMIX, ct->dim, ct->q_old, ct->q_act, broyd);
        for (i = 0; i < ct->dim; i++)
        {
            ct->q_act[i] = ct->q_old[i];
        }
        // correct numerical errors caused by broyden
        q_sum = 0;
        for (i = 0; i < ct->dim; i++)
        {
            if (ct->q_act[i] < 0.0)
            {
                ct->q_act[i] = 0.0;
            }
            if (ct->q_act[i] > 1.0)
            {
                ct->q_act[i] = 1.0;
            }
            q_sum += ct->q_act[i];
        }
        for (i = 0; i < ct->dim; i++)
        {
            ct->q_act[i] = ct->q_act[i] / q_sum; //normalization


        }
    } /* end SCC cycle */

    fprintf(f, " %2d ", step);
    //fprintf(f, " %2d iters, SCC-Q =", step);
    for (i = 0; i < ct->dim; i++)
    {
        fprintf(f, "%9.6f", ct->q_act[i]);
    }


    // copy the lowest pure adiabatic state to the wave function
    //for (i=0; i<ct->dim; i++)
    //  ct->wf[i] = ham[i];



    // calculate and print the energies of pure adiabatic states (indexed by k)
    //fprintf(f, " E(adiab) =");
    for (k = 0; k < ct->dim; k++)
    {
        energy = 0.0;
        for (i = 0; i < ct->dim; i++)
        {
            for (j = 0; j < ct->dim; j++)
            {
                // TB Hamiltonian
                energy += ham[k * ct->dim + i] * ham[k * ct->dim + j] * ct->hamiltonian[i][j]
                    // Hubbard / gamma terms
                    +  0.5 * ct->hubbard[i][j] * SQR(ham[k * ct->dim + i]) * SQR(ham[k * ct->dim + j]);
            }
        }
        fprintf(f, "%10.7f", energy);
    }




///*
// calc overlap with adiabatic states i //
    for (i = 0; i < ct->dim; i++)
    {
        ct->born_overlap[i] = 0.0;
        for (j = 0; j < ct->dim; j++)
        {
            ct->born_overlap[i] +=   ham[i*ct->dim + j]*ct->wf[j];
        }
    }

    // wf will be adiabatic state that is most similar //
    best_overlap = 0.0;
    for (i = 0; i < ct->dim; i++)
    {
        if (fabs(ct->born_overlap[i]) > fabs(best_overlap))
        {
            best_overlap = ct->born_overlap[i];
            for (j = 0; j < ct->dim; j++)
            {
                ct->wf[j] = ham[i*ct->dim + j];
            }

            printf("best overlap with adiabatic state no. %d. overlap %lf \n", i, ct->born_overlap[i]);
        }
    }
//*/

/*
   for(i=0;i<ct->dim;i++){ //calc overlap with state i
    s01 = 0.e0;
    for (j=0; j<ct->dim; j++)
      s01 += ham[i*ct->dim + j]*ct->wf[j];
   printf("S01 %f\n", s01);
    if(fabs(s01)>0.8){ //threshold
        for (j=0; j<ct->dim; j++)
          ct->wf[j] = ham[i*ct->sites + j];
    break;
    }
   }//end loop over eigenvalues
   if(i==ct->dim){ printf("Eigenfunctions changed too much! Born-Oppenheimer not reasonable!"); exit(-1); }
 */
    return SQR(best_overlap);
}













int do_adiab_fermi(charge_transfer_t *ct, dftb_broyden_t *broyd, double *fermi_coeff, FILE *f)
{
    int    i, j, k, step;
    double old_energy, energy, energy1, energy2, *ham, fermi_energy, fermi_upper, fermi_lower, q_sum;

    ham = ct->hamiltonian_adiab;

    // calculate the charges from the wave function
    for (i = 0; i < ct->dim; i++)
    {
        ct->q_act[i] = ct->q_old[i] = SQR(ct->wf[i]);
    }

    // copy as the first vector

    energy = 1.e10;

    // SCC cycle
    for (step = 0;; step++)
    {

        // too many iterations? break and indicate failure
        if (step >= MAXITER_BROYDEN)
        {
            fprintf(f, "Broyden/Fermi failed to converge!\n");
            return 1;
        }

        // construct the self-consistent Hamiltonian
        for (i = 0; i < ct->dim; i++)
        {
            ham[i + i * ct->dim] = ct->hamiltonian[i][i];
            for (j = 0; j < ct->dim; j++)
            {
                ham[i + i * ct->dim] += ct->hubbard[i][j] * ct->q_act[j];
                if (i != j)
                {
                    ham[i + j * ct->dim] = ct->hamiltonian[i][j];
                }
            }
        }

        // diagonalize this Hamiltonian
        dsyev(ct->dim, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);

        // FERMI DISTRIBUTION + BROYDEN MIXING
        // mix all eigenvectors with appropriate coefficients

        // calculate the "Fermi energy", store the Fermi coefficients in fermi_coeff
        fermi_lower = ct->ev_adiab[0] - 0.1;
        fermi_upper = ct->ev_adiab[ct->dim - 1] + 0.1;
        while (fermi_upper - fermi_lower > FERMI_CONVERG)
        {
            fermi_energy = 0.5 * (fermi_lower + fermi_upper);
            if (fermi_coef_sum(fermi_energy, ct->dim, ct->ev_adiab, fermi_coeff, ct->fermi_kt) > 1.0)
            {
                fermi_upper = 0.5 * (fermi_lower + fermi_upper);
            }
            else
            {
                fermi_lower = 0.5 * (fermi_lower + fermi_upper);
            }
        }
        fprintf(f, " fermi E = %9.6f, coefs", fermi_energy);
        for (i = 0; i < ct->dim; i++)
        {
            fprintf(f, "%8.5f", fermi_coeff[i]);
        }

        // mix the eigenvectors, with the appropriate coefficients
        for (i = 0; i < ct->dim; i++) // i runs over the sites
        {
            ct->q_act[i] = 0.e0;
            for (j = 0; j < ct->dim; j++) // j runs over the individual eigenvectors
            {
                ct->q_act[i] += fermi_coeff[j] * SQR(ham[j * ct->dim + i]);
            }
        }

        // print new charges
        fprintf(f, " mixedQ =");
        for (i = 0; i < ct->dim; i++)
        {
            fprintf(f, "%9.6f", ct->q_act[i]);
        }

        // calculate the energy and check convergence
        old_energy = energy;
        energy1    = energy2 = 0.0;
        // indices i and j run over the sites
        for (i = 0; i < ct->dim; i++)
        {
            for (j = 0; j < ct->dim; j++)
            {
                // TB Hamiltonian
                // index k runs over the adiabatic states (eigenvectors)
                for (k = 0; k < ct->dim; k++)
                {
                    energy1 += fermi_coeff[k] * ham[k * ct->dim + i] * ham[k * ct->dim + j] * ct->hamiltonian[i][j];
                }
                // Hubbard / gamma terms
                energy2 += 0.5 * ct->hubbard[i][j] * ct->q_act[i] * ct->q_act[j];
            }
        }
        energy = energy1 + energy2;

        fprintf(f, " E1=%10.7f E2=%10.7f totE=%13.10f", energy1, energy2, energy);

        if (step > 1 && fabs(energy - old_energy) < BROYDEN_SCFTOL)
        {
            fprintf(f, " Converged!\n");
            // copy the solution to the wave function
            for (i = 0; i < ct->dim; i++)
            {
                ct->wf[i] = sqrt(ct->q_act[i]);
            }
            break;
        }

        // mix the charges - Broyden
        broyden(step, BROYDEN_ALMIX, ct->dim, ct->q_old, ct->q_act, broyd);
        for (i = 0; i < ct->dim; i++)
        {
            ct->q_act[i] = ct->q_old[i];
        }
        // correct numerical errors caused by broyden
        q_sum = 0;
        for (i = 0; i < ct->dim; i++)
        {
            if (ct->q_act[i] < 0.0)
            {
                ct->q_act[i] = 0.0;
            }
            if (ct->q_act[i] > 1.0)
            {
                ct->q_act[i] = 1.0;
            }
            q_sum += ct->q_act[i];
        }
        for (i = 0; i < ct->dim; i++)
        {
            ct->q_act[i] = ct->q_act[i] / q_sum; //normalization


        }
        // print mixed charges
        fprintf(f, " newQ =");
        for (i = 0; i < ct->dim; i++)
        {
            fprintf(f, "%9.6f", ct->q_act[i]);
        }
        fprintf(f, "\n");

    } /* end SCC cycle */

    return 0;
}

int do_adiab_fermi_onestate(charge_transfer_t *ct, dftb_broyden_t *broyd, double *fermi_coeff, FILE *f)
{
    int    i, j, k, step;
    double old_energy, energy, energy1, energy2, *ham, fermi_energy, fermi_upper, fermi_lower, q_sum;

    ham = ct->hamiltonian_adiab;

    // calculate the charges from the wave function
    for (i = 0; i < ct->dim; i++)
    {
        ct->q_act[i] = ct->q_old[i] = SQR(ct->wf[i]);
    }

    // copy as the first vector

    energy = 1.e10;
    // SCC cycle
    for (step = 0;; step++)
    {

        // too many iterations? break and indicate failure
        if (step >= MAXITER_BROYDEN)
        {
            fprintf(f, "Broyden/Fermi failed to converge!\n");
            return 1;
        }

        // construct the self-consistent Hamiltonian
        for (i = 0; i < ct->dim; i++)
        {
            ham[i + i * ct->dim] = ct->hamiltonian[i][i];
            for (j = 0; j < ct->dim; j++)
            {
                ham[i + i * ct->dim] += ct->hubbard[i][j] * ct->q_act[j];
                if (i != j)
                {
                    ham[i + j * ct->dim] = ct->hamiltonian[i][j];
                }
            }
        }

        // diagonalize this Hamiltonian
        dsyev(ct->dim, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);

        // FERMI DISTRIBUTION + BROYDEN MIXING
        // mix all eigenvectors with appropriate coefficients

        // calculate the "Fermi energy", store the Fermi coefficients in fermi_coeff
        fermi_lower = ct->ev_adiab[0] - 0.1;
        fermi_upper = ct->ev_adiab[ct->dim - 1] + 0.1;
        while (fermi_upper - fermi_lower > FERMI_CONVERG)
        {
            fermi_energy = 0.5 * (fermi_lower + fermi_upper);
            if (fermi_coef_sum(fermi_energy, ct->dim, ct->ev_adiab, fermi_coeff, ct->fermi_kt) > 1.0)
            {
                fermi_upper = 0.5 * (fermi_lower + fermi_upper);
            }
            else
            {
                fermi_lower = 0.5 * (fermi_lower + fermi_upper);
            }
        }

        // mix the eigenvectors, with the appropriate coefficients
        for (i = 0; i < ct->dim; i++) // i runs over the sites
        {
            ct->q_act[i] = 0.e0;
            for (j = 0; j < ct->dim; j++) // j runs over the individual eigenvectors
            {
                ct->q_act[i] += fermi_coeff[j] * SQR(ham[j * ct->dim + i]);
            }
        }

        // calculate the energy and check convergence
        old_energy = energy;
        energy1    = energy2 = 0.0;
        // indices i and j run over the sites
        for (i = 0; i < ct->dim; i++)
        {
            for (j = 0; j < ct->dim; j++)
            {
                // TB Hamiltonian
                // index k runs over the adiabatic states (eigenvectors)
                for (k = 0; k < ct->dim; k++)
                {
                    energy1 += fermi_coeff[k] * ham[k * ct->dim + i] * ham[k * ct->dim + j] * ct->hamiltonian[i][j];
                }
                // Hubbard / gamma terms
                energy2 += 0.5 * ct->hubbard[i][j] * ct->q_act[i] * ct->q_act[j];
            }
        }
        energy = energy1 + energy2;

        if (step > 1 && fabs(energy - old_energy) < BROYDEN_SCFTOL)
        {
            break;
        }

        // mix the charges - Broyden
        broyden(step, BROYDEN_ALMIX, ct->dim, ct->q_old, ct->q_act, broyd);
        for (i = 0; i < ct->dim; i++)
        {
            ct->q_act[i] = ct->q_old[i];
        }
        // correct numerical errors caused by broyden
        q_sum = 0;
        for (i = 0; i < ct->dim; i++)
        {
            if (ct->q_act[i] < 0.0)
            {
                ct->q_act[i] = 0.0;
            }
            if (ct->q_act[i] > 1.0)
            {
                ct->q_act[i] = 1.0;
            }
            q_sum += ct->q_act[i];
        }
        for (i = 0; i < ct->dim; i++)
        {
            ct->q_act[i] = ct->q_act[i] / q_sum; //normalization


        }
    } /* end SCC cycle */

    fprintf(f, " %2d ", step);
    //fprintf(f, " %2d iters, SCC-Q =", step);
    for (i = 0; i < ct->dim; i++)
    {
        fprintf(f, "%9.6f", ct->q_act[i]);
    }

    // copy the lowest pure adiabatic state to the wave function
    for (i = 0; i < ct->dim; i++)
    {
        ct->wf[i] = ham[i];
    }

    // calculate and print the energies of pure adiabatic states (indexed by k)
    //fprintf(f, " E(adiab) =");
    for (k = 0; k < ct->dim; k++)
    {
        energy = 0.0;
        for (i = 0; i < ct->dim; i++)
        {
            for (j = 0; j < ct->dim; j++)
            {
                // TB Hamiltonian
                energy += ham[k * ct->dim + i] * ham[k * ct->dim + j] * ct->hamiltonian[i][j]
                    // Hubbard / gamma terms
                    +  0.5 * ct->hubbard[i][j] * SQR(ham[k * ct->dim + i]) * SQR(ham[k * ct->dim + j]);
            }
        }
        fprintf(f, "%10.7f", energy);
    }
    return 0;
}

// Trajectory surface hopping method using original diabatic Landau-Zener formula W.X Feb.2019


int do_lzsh_diab(int step, charge_transfer_t *ct, dftb_t *dftb, t_state *state_global, t_mdatoms *mdatoms, rvec *x_ct, FILE *f)
{

    int    i, j, k, l, m, iatom, hopsurf, prev_surf, indextmp, counter;
    double *ham, *ekin, *dist, de[1];
    double probmax, random_number, deriv, coupling;
    double ekintot, ekintmp, kai, detmp, disttmp, norm;
    dvec   com, coord, masscoord;


    snew(ekin, ct->dim);
    snew(dist, ct->dim);
//    snew(index, ct->dim);

    ham = ct->tfl_mean_ham;

/*    counter = 0;
    for (i = 0; i < ct->pool_size; i++)
    {
        clear_dvec(com);
        for (j = 0; j < ct->pool_site[i].atoms; j++)
        {
            for (k = 0; k < 3; k++)
            {
                coord[k] =  x_ct[ct->pool_site[i].atom[j]][k] * NM_TO_BOHR;//take here x_ct for wich molecules were made whole again. in gromacs_x all atoms are put in the box.
            }
            dsvmul(mdatoms->massT[ct->pool_site[i].atom[j]], coord, masscoord);
            dvec_inc(com, masscoord);
        }
        dsvmul(ct->adapt_inv_tot_mass, com, ct->pool_site[i].com);
        printf("residue %d  COM[Angstrom]: %f %f %f\n", ct->pool_site[i].resnr, ct->pool_site[i].com[XX] * 0.52, ct->pool_site[i].com[YY] * 0.52, ct->pool_site[i].com[ZZ] * 0.52);
    }

   exit(0);
*/
    for (i = 0; i < SQR(ct->dim); i++)
        ham[i] = *(ct->hamiltonian[0]+i);

// diagonalize the Hamiltonian
    dsyev(ct->dim, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);

// record eigenvectors
    for (i = 0; i < ct->dim; i++)
        for (j = 0; j < ct->dim; j++)
            ct->tfs_vector[i][j] = ham[i * ct->dim + j];

// set up initial surface
    if (ct->tfs_initialization_step)
    {
// determine inital surface using eigenvector as a criterion
       for (i = 0; i < 2*ct->dim; i++)
       {
           ct->tfs_popul[i] = 0.0;
           ct->tfs_diab[i] = 0.0;
       }

       for (i = 0; i < ct->dim; i++)
       {
           ct->tfs_diab[i] = ct->wf[i];
           if (ct->wf[i] == 1.0)
              ct->surface = i;
       }
       ct->tfs_initialization_step = 0;
    }

    if (step > 0)
    {
       probmax = 0.0;
       for (i = 0; i < ct->dim; i++)
       {
// calculate the diabatic energy gaps at two time steps before and previous time step
           de[0] = ct->tfl_mean_ham_full[ct->surface][ct->surface] - ct->tfl_mean_ham_full[i][i];
           de[1] = ct->hamiltonian[ct->surface][ct->surface] - ct->hamiltonian[i][i];
           coupling = (fabs(ct->tfl_mean_ham_full[ct->surface][i]) + fabs(ct->hamiltonian[ct->surface][i]))/2.0;
// examine local energy minimum
           if (de[0]*de[1] < 0.0)
           {
// derivative of diabatic energy gap with respect to time
              deriv = (de[1] - de[0])/ct->rk_timestep;
// Landau-Zener formula (Phy.Rev.A 84, 014701, 2011)
              ct->surf_prob[i] = 1.0 - exp(-2.0*M_PI*SQR(coupling)/fabs(deriv));

              printf("%f \n", ct->surf_prob[i]);
              if (ct->surf_prob[i] > probmax)
              {
                 probmax = ct->surf_prob[i];
                 hopsurf = i;
              }
           }
       }

// compared with a random number to determine the hopping event
       prev_surf = ct->surface;
       random_number = drand48();
       if (probmax > random_number)
       {
// kinetic energy of QM zone
           for (k = 0; k < ct->sites; k++)
               ekin[k] = 0.0;

           ekintot = 0.0;
           for (iatom = mdatoms->start; iatom < mdatoms->start+mdatoms->homenr; iatom++)
           {
               for (k = 0; k < ct->sites; k++)
               {
                   for (l = 0; l < ct->site[k].atoms; l++)
                   {
                       if (iatom == ct->site[k].atom[l])
                       {
                          for (m = 0; m < DIM; m++)
                          {
                              ekintot += 0.5*dftb->phase1[k].mass[l]*AMU_TO_AU*SQR(state_global->v[iatom][m]*NM_TO_BOHR/PS_TO_AU);
                              ekin[k] += 0.5*dftb->phase1[k].mass[l]*AMU_TO_AU*SQR(state_global->v[iatom][m]*NM_TO_BOHR/PS_TO_AU);
                          }
                       }
                   }
               }
           }
           printf("Kinetic energy is %f eV\n", ekintot*27.2114);

// frustrated hops
           ekintmp = ekin[hopsurf] + ekin[ct->surface];
           detmp = ct->hamiltonian[hopsurf][hopsurf] - ct->hamiltonian[ct->surface][ct->surface];

           printf("Two-site kinetic energy is %f eV \n", ekintmp*27.2114);
           printf("Potential energy difference is %f eV \n", detmp*27.2114);

           if (1.0 - detmp/ekintmp < 0.0)
           {
              printf("Frustrated hop! \n");
              fprintf(f,"%8.2f %10d %17d %21f %19f Frustrated\n",(step-1)*ct->rk_timestep*1000./PS_TO_AU, ct->surface, hopsurf, probmax, random_number);
           }
           else
           {
              kai = sqrt(1.0 - detmp/ekintmp);

              for (iatom = mdatoms->start; iatom < mdatoms->start+mdatoms->homenr; iatom++)
              {
                  for (k = 0; k < ct->sites; k++)
                  {
                      if (k == hopsurf || k == ct->surface)
                      {
                          for (l = 0; l < ct->site[k].atoms; l++)
                          {
                              if (iatom == ct->site[k].atom[l])
                              {
                                 for (m = 0; m < DIM; m++)
                                     state_global->v[iatom][m] *= kai;
                              }
                          }
                      }
                  }
              }

              fprintf(f,"%8.2f %10d %17d %21f %19f \n",(step-1)*ct->rk_timestep*1000./PS_TO_AU, ct->surface, hopsurf, probmax, random_number);

              ct->tfs_diab[ct->surface] = 0.0;
              ct->wf[ct->surface]       = 0.0;

              ct->surface = hopsurf;

              ct->tfs_diab[ct->surface] = 1.0;
              ct->wf[ct->surface]       = 1.0;
           }
       }
    }

// adaptive sites
/*    if (ct->surface != prev_surf)
    {
    counter = 0;
    for (i = 0; i < ct->pool_size; i++)
    {
        clear_dvec(com);
        for (j = 0; j < ct->pool_site[i].atoms; j++)
        {
            for (k = 0; k < 3; k++)
            {
                coord[k] =  x_ct[ct->pool_site[i].atom[j]][k] * NM_TO_BOHR;//take here x_ct for wich molecules were made whole again. in gromacs_x all atoms are put in the box.
            }
            dsvmul(mdatoms->massT[ct->pool_site[i].atom[j]], coord, masscoord);
            dvec_inc(com, masscoord);
        }
        dsvmul(ct->adapt_inv_tot_mass, com, ct->pool_site[i].com);
        printf("residue %d  COM[Angstrom]: %f %f %f\n", ct->pool_site[i].resnr, ct->pool_site[i].com[XX] * 0.52, ct->pool_site[i].com[YY] * 0.52, ct->pool_site[i].com[ZZ] * 0.52);
    }

    for (i = 0; i < ct->dim; i++)
    {
        for (j = i + 1; j < ct->dim; j++)
        {
            if (dist[i] > dist[j])
            {
               disttmp = dist[i];
               dist[i] = dist[j];
               dist[j] = disttmp;

               indextmp = index[i];
               index[i] = index[j];
               index[j] = indextmp;
            }
        }
    }

    }
*/
// store the old hamiltonian
    for (i = 0; i < ct->dim; i++)
        for (j = 0; j < ct->dim; j++)
            ct->tfl_mean_ham_full[i][j] = ct->hamiltonian[i][j];

// calculate adibatic expansion coeffients
    for (i = 0; i < 2*ct->dim; i++)
        ct->tfs_popul[i] = 0.0;

    for (i = 0; i < ct->dim; i++)
        ct->tfs_popul[i] = ct->tfs_vector[i][ct->surface];

    printf("Current Surface: %d\n", ct->surface);

    return 0;
}
// Trajectory surface hopping method using refomulated Landau-Zener formula W.X Nov.2018
// Ref: Phys.Rev A 84, 014701 (2011)

int do_lzsh_adiab(int step, charge_transfer_t *ct, FILE *f)
{

    int    i, j, k, l, hopsurf, *ind;
    double *ham, *initprob, de[2], de_diab[1];
    double probsum, probmax, random_number, deriv2, coupling, deriv;
    double ekintot, ekintmp, detmp;
    double nc, aa, bb, cc, dd, ee;


    ham = ct->tfl_mean_ham;

    snew(initprob, ct->dim);
    snew(ind, ct->dim);

    for (i = 0; i < SQR(ct->dim); i++)
        ham[i] = *(ct->hamiltonian[0]+i);

// diagonalize the Hamiltonian
    dsyev(ct->dim, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);

// record eigenvectors
    for (i = 0; i < ct->dim; i++)
        for (j = 0; j < ct->dim; j++)
            ct->tfs_vector[i][j] = ham[i * ct->dim + j];

// set up initial surface
    if (ct->tfs_initialization_step)
    {
// determine inital surface using eigenvector as a criterion
       probsum = 0.;
       random_number = drand48();
       for (i = 0; i < ct->dim; i++)
       {
           initprob[i] = 0.;
           for (j = 0; j < ct->dim; j++)
               initprob[i] += ct->tfs_vector[i][j]*ct->wf[j];
           initprob[i] = SQR(initprob[i]);

           probsum += initprob[i];
           if (random_number < probsum)
           {
               ct->surface = i;
               printf("\n Initial surface: %10d %12.5f %12.5f \n", ct->surface, probsum, random_number);
               break;
           }
       }

       for (i = 0; i < 2*ct->dim; i++)
       {
           ct->tfs_diab[i] = 0.0;
           ct->wf[i]       = 0.0;
       }

       for (i = 0; i < ct->dim; i++)
       {
           ct->tfs_diab[i] = ct->tfs_vector[ct->surface][i];
           ct->wf[i]       = ct->tfs_vector[ct->surface][i];
       }

       for (i = 0; i < 2*ct->dim; i++)
           ct->tfs_popul[i] = 0.0;

       ct->tfs_popul[ct->surface] = 1.0;

       ct->tfs_initialization_step = 0;
    }

// calculate the adiabatic energy gaps at current time
    if (step > 1)
    {
       probmax = 0.;
       for (i = 0; i < ct->dim; i++)
       {
// calculate the adiabatic energy gaps at two time steps before and previous time step
           de[0] = fabs(ct->ev_adiab_old2[ct->surface] - ct->ev_adiab_old2[i]);
           de[1] = fabs(ct->ev_adiab_old[ct->surface] - ct->ev_adiab_old[i]);
           de[2] = fabs(ct->ev_adiab[ct->surface] - ct->ev_adiab[i]);

           de_diab[0] = ct->ev_adiab[i] - ct->ev_adiab_old2[ct->surface];
           de_diab[1] = ct->ev_adiab[ct->surface] - ct->ev_adiab_old2[i];

           coupling = de[1]/2.0;
// calculate hopping probability
           ct->surf_prob[i] = 0.;
// examine local energy minimum
           if (de[1] < de[0] && de[1] < de[2])
           {
// second derivative of adiabatic energy gap with respect to time
              deriv2 = (de[0] + de[2] - 2.0*de[1])/SQR(ct->rk_timestep);
// Landau-Zener formula (Phy.Rev.A 84, 014701, 2011)
              ct->surf_prob[i] = exp(-M_PI/2.0*sqrt(CUB(de[1])/deriv2));
// derivative of diabatic energy gap with respect to time
//              deriv = fabs(de_diab[1] - de_diab[0])/ct->rk_timestep;
// Landau-Zener formula (Phy.Rev.A 84, 014701, 2011)
//              ct->surf_prob[i] = exp(-2.0*M_PI*SQR(coupling)/fabs(deriv));

              if (ct->surf_prob[i] > probmax)// && fabs(ct->surface - i) < 3)
              {
                 probmax = ct->surf_prob[i];
                 hopsurf = i;
              }
           }
       }

// compared with a random number to determine the hopping event
       random_number = drand48();
       if (probmax > random_number)
       {
          fprintf(f,"%8.2f %10d %17d %21f %19f \n",(step-1)*ct->rk_timestep*1000./PS_TO_AU, ct->surface, hopsurf, probmax, random_number);
          ct->surface = hopsurf;
       }
    }

// calculate diabatic expansion coeffients
    for (i = 0; i < ct->dim; i++)
    {
        ct->tfs_diab[i] = ct->tfs_vector[ct->surface][i];
        ct->wf[i]       = ct->tfs_vector[ct->surface][i];
    }

    for (i = 0; i < ct->dim; i++)
    {
        ct->ev_adiab_old2[i] = ct->ev_adiab_old[i];
        ct->ev_adiab_old[i]  = ct->ev_adiab[i];
    }

// calculate adibatic expansion coeffients
    for (i = 0; i < 2*ct->dim; i++)
        ct->tfs_popul[i] = 0.0;

    ct->tfs_popul[ct->surface] = 1.0;

    printf("Current Surface: %d\n", ct->surface);

    return 0;
}

int do_fssh_diab(int step, charge_transfer_t *ct, dftb_t *dftb, t_state *state_global, t_mdatoms *mdatoms, FILE *f)
{
// Declaration:
//
// surface hopping betwteen different sites
//
//                                          -----------Weiwei Xie on Jan. 2019


    int    i, j, k, l, m, n, prev_surf, hopsurf, iatom, istep, nstep, hop, lamtmp;
    double *ham, *old_ham, **ham_full, *new_ham, *diapop, *diapop_old, *ekin, *ddiapop;
    double tdse_timestep = 1.e-5*PS_TO_AU;
    double cur_weight, rk_timestep0, probsum, random_number;
    double decay_time, curr_surf_fac, ddiapopsum;
    double dot_product, ekintot, ekintmp, de, kai;
    double nc, kai1, kai2, delta, aa, bb, cc, dd, ee;


    old_ham  = ct->tfl_old_ham_adiab;
    new_ham  = ct->hamiltonian_adiab;
    ham      = ct->tfl_mean_ham;
    ham_full = ct->tfl_mean_ham_full;

    snew(diapop_old, ct->dim);
    snew(diapop, ct->dim);
    snew(ddiapop, ct->dim);
    snew(ekin, ct->dim);

// set up initial condition
    if (ct->tfs_initialization_step)
    {
       printf("Input wavefunction is:\n");
       for (j = 0; j < ct->dim; j++)
           printf("%f ", ct->wf[j]);
       printf("\n");

       for (i = 0; i < SQR(ct->dim); i++)
       {
           ct->tfl_old_ham[i] = *(ct->hamiltonian[0]+i);
           ham[i] = ct->tfl_old_ham[i];
       }

// diagonalize this Hamiltonian
       dsyev(ct->dim, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);

       for (i = 0; i < ct->dim; i++)
           for (j = 0; j < ct->dim; j++)
               ct->tfs_vector[i][j] = ham[i * ct->dim + j];

       for (i = 0; i < 2*ct->dim; i++)
       {
           ct->tfs_popul[i] = 0.0;
           ct->tfs_diab[i] = 0.0;
       }

       for (i = 0; i < ct->dim; i++)
       {
           ct->tfs_diab[i] = ct->wf[i];
           if (ct->wf[i] == 1.0)
              ct->surface = i;
       }

       for (i = 0; i < ct->dim; i++)
           ct->tfs_popul[i] = ct->tfs_vector[i][ct->surface];

       ct->tfs_initialization_step = 0;

       return 0;
    }

// calculate the diabatic populations at the previous time step
    for (i = 0; i < ct->dim; i++)
        diapop_old[i]  = SQR(ct->tfs_diab[i]) + SQR(ct->tfs_diab[i+ct->dim]);

// set up the old and new Hamiltonian for linear interpolation
    for (i = 0; i < ct->dim; i++)
        for (j = 0; j < ct->dim; j++)
        {
            old_ham[j+i*ct->dim] = ct->tfl_old_ham[i*ct->dim+j];
            new_ham[j+i*ct->dim] = ct->hamiltonian[i][j];
        }

/*    for (m= 0; m < ct->dim; m++)
    {
        for (n = 0; n < ct->dim; n++)
        {
            ct->fo_overlap[m][n] = 0.0;
            for (i = 0; i < ct->dim; i++)
            {
                for (j = 0; j < ct->dim; j++)
                {
                   for (k = 0; k < dftb->phase2.norb; k++)
                   {
                       for (l = 0; l < dftb->phase2.norb; l++)
                       {
                           ct->fo_overlap[m][n] += dftb->orthogo.sij_old[i+m*ct->dim]*dftb->orthogo.sij[j+n*ct->dim]
                                         *dftb->orthogo.fmo_old[k][i]*dftb->orthogo.fmo[l][j]
                                         *0.5*(dftb->orthogo.sao_old[k][l]+dftb->orthogo.sao[k][l]);
                       }
                   }

                }
            }
        }
    }
    for (m= 0; m < ct->dim; m++)
    {
        for (n = 0; n < ct->dim; n++)
        {
            ct->fo_tdc[m][n] = (ct->fo_overlap[m][n]- ct->fo_overlap[n][m])/2.0/ct->rk_timestep;
//            printf("%d %d %15.5f  %15.5f  %15.5f  \n", m, n, ct->fo_overlap[m][n], ct->fo_overlap[n][m], ct->fo_tdc[m][n]);
        }
    }
*/
    for (i = 0; i < ct->dim; i++)
        ct->surf_prob[i] = 0.0;

// solve time-dependent Schroedinger equation
    nstep = (int)(ct->rk_timestep/tdse_timestep);
    rk_timestep0 = ct->rk_timestep;
    ct->rk_timestep = tdse_timestep;

    for (istep = 1; istep < nstep+1; istep++)
    {
// linear interpolation for Hamiltonian
        cur_weight = (double)istep/(double)nstep;
        for (i = 0; i < ct->dim; i++)
            for (j = 0; j < ct->dim; j++)
            {
                ham_full[i][j] = (1.0 - cur_weight)*ct->tfl_old_ham[i*ct->dim+j] + cur_weight*ct->hamiltonian[i][j];
//                printf("%d  %d  %f  %f \n", i, j, ham_full[i][j], ct->tfl_mean_ham_full[i][j]);
            }
// propagate the TDSE in the diabatic representation
//        RK5th(ct);
        do_rksuite_diab(ct);

/*        for (i = 0; i < SQR(ct->dim); i++)
            ham[i] = (1.0 - cur_weight)*old_ham[i] + cur_weight*new_ham[i];

        for (i = 0; i < ct->dim; i++)
            for (j = 0; j < ct->dim; j++)
                ct->surf_prob[i] += (2.0*ham[i*ct->dim+j]*(ct->tfs_diab[i]*ct->tfs_diab[ct->dim+j]-ct->tfs_diab[j]*ct->tfs_diab[ct->dim+i]))*ct->rk_timestep/diapop_old[ct->surface];
*/
    }
    ct->rk_timestep = rk_timestep0;

// record the adiabatic populations at new time
    ddiapopsum = 0.0;
    for (i = 0; i < ct->dim; i++)
    {
        diapop[i] = SQR(ct->tfs_diab[i])+SQR(ct->tfs_diab[i+ct->dim]);

        ddiapop[i] = diapop[i] - diapop_old[i];

        if (ddiapop[i] < 0.0)
           ddiapopsum += ddiapop[i];
    }

// hopping probability
    for (i = 0; i < ct->dim; i++)
    {
        if (ddiapop[ct->surface] < 0.0 && ddiapop[i] > 0.0)
           ct->surf_prob[i] = ddiapop[i]/diapop_old[ct->surface]*ddiapop[ct->surface]/ddiapopsum;
        else
           ct->surf_prob[i] = 0.0;
    }

/*    for (i = 0; i < ct->dim; i++)
        ct->surf_prob[i] = ddiapop[i]/diapop_old[ct->surface];
*/
    for (i = 0; i < ct->dim; i++)
        printf("%2d  %15.5f  %15.5f \n", i, ddiapop[i]/diapop_old[ct->surface], ct->surf_prob[i]);

// kinetic energy of QM zone
    for (k = 0; k < ct->sites; k++)
        ekin[k] = 0.0;

    ekintot = 0.0;
    for (iatom = mdatoms->start; iatom < mdatoms->start+mdatoms->homenr; iatom++)
    {
        for (k = 0; k < ct->sites; k++)
        {
            for (l = 0; l < ct->site[k].atoms; l++)
            {
                if (iatom == ct->site[k].atom[l])
                {
                   for (m = 0; m < DIM; m++)
                   {
                       ekintot += 0.5*dftb->phase1[k].mass[l]*AMU_TO_AU*SQR(state_global->v[iatom][m]*NM_TO_BOHR/PS_TO_AU);
                       ekin[k] += 0.5*dftb->phase1[k].mass[l]*AMU_TO_AU*SQR(state_global->v[iatom][m]*NM_TO_BOHR/PS_TO_AU);
                   }
                }
            }
        }
    }

    printf("Kinetic energy is %f %f eV\n", ekintot*27.2114, BOLTZMANN_HARTREE_KELVIN);

// stochastic algorithm
    probsum = 0.0;
    random_number = drand48();
    for (i = 0; i < ct->dim; i++)
    {
        if (ct->surf_prob[i] > 0.0 && i != ct->surface)
        {
           probsum += ct->surf_prob[i];

           if (random_number < probsum)
           {
              ekintmp = ekin[i] + ekin[ct->surface];
              de = ct->hamiltonian[i][i] - ct->hamiltonian[ct->surface][ct->surface];

              printf("Two-site kinetic energy is %f eV \n", ekintmp*27.2114);
              printf("Potential energy difference is %f eV \n", de*27.2114);

              if (1.0 - de/ekintmp < 0.0)
              {
                 printf("Frustated hop! \n");
                 fprintf(f,"%8.2f %10d %17d %21f %19f Frustrated \n",(step*ct->rk_timestep)*1000./PS_TO_AU, ct->surface, i, ct->surf_prob[i], random_number);
                 break;
              }
              else
              {

                 kai = sqrt(1.0 - de/ekintmp);

                 ekintot -= de;

                 for (iatom = mdatoms->start; iatom < mdatoms->start+mdatoms->homenr; iatom++)
                 {
                     for (k = 0; k < ct->sites; k++)
                     {
                         if (k == i || k == ct->surface)
                         {
                             for (l = 0; l < ct->site[k].atoms; l++)
                             {
                                 if (iatom == ct->site[k].atom[l])
                                 {
                                    for (m = 0; m < DIM; m++)
                                        state_global->v[iatom][m] *= kai;
                                 }
                             }
                         }
                     }
                 }
              }

              fprintf(f,"%8.2f %10d %17d %21f %19f \n",(step*ct->rk_timestep)*1000./PS_TO_AU, ct->surface, i, ct->surf_prob[i], random_number);

              ct->surface = i;

              break;
           }
        }
    }

    printf("Current Surface: %d  %f \n", ct->surface);

// decoherence corrections
    dot_product = 0.0;
    for (i = 0; i < ct->dim; i++)
    {
        if (i != ct->surface)
        {
            decay_time = (1.0 + 0.1/ekintot)/fabs(ct->hamiltonian[i][i] - ct->hamiltonian[ct->surface][ct->surface]);

            ct->tfs_diab[i]         *= exp(-ct->rk_timestep/decay_time);
            ct->tfs_diab[i+ct->dim] *= exp(-ct->rk_timestep/decay_time);

            dot_product              += SQR(ct->tfs_diab[i])+SQR(ct->tfs_diab[i+ct->dim]);
        }
    }
    curr_surf_fac = sqrt((1.0 - dot_product)/(SQR(ct->tfs_diab[ct->surface]) + SQR(ct->tfs_diab[ct->dim+ct->surface])));

    ct->tfs_diab[ct->surface]         *= curr_surf_fac;
    ct->tfs_diab[ct->surface+ct->dim] *= curr_surf_fac;

// new state
    for (i = 0; i < ct->dim; i++)
        ct->wf[i] = 0.0;
    ct->wf[ct->surface] = 1.0;

    for (i = 0; i < SQR(ct->dim); i++)
        ham[i] = new_ham[i];

    dsyev(ct->dim, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);

// calculate expansion coefficients in the adiabatic basis
    for (i = 0; i < ct->dim; i++)
        for (j = 0; j < ct->dim; j++)
            ct->tfs_vector[i][j] = ham[i * ct->dim + j];

    for (i = 0; i < ct->dim; i++)
    {
        ct->tfs_popul[i] = 0.0; ct->tfs_popul[i+ct->dim] = 0.0;
        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_popul[i]         += ct->tfs_vector[i][j]*ct->tfs_diab[j];
            ct->tfs_popul[i+ct->dim] += ct->tfs_vector[i][j]*ct->tfs_diab[j+ct->dim];

        }
    }

//store hamiltonian as old_ham
    for (i = 0; i < SQR(ct->dim); i++)
        ct->tfl_old_ham[i] = *(ct->hamiltonian[0]+i);

    return 0;
}

int do_fssh_johen(int step, charge_transfer_t *ct, dftb_t *dftb, t_state *state_global, t_mdatoms *mdatoms, FILE *f)
{
// Declaration:
//
// 1. a correction is implemented for decoherence correction-induced spurious long-range charge transfercorrection for the decoherence;
// 2. self-consistent algorithm is applied for avoided crossings;
// 3. state-tracking algorithm is implemented for trivial corrssing detection;
//
//                                          -----------Weiwei Xie on Feb. 2019
// Ref: J. Phys. Chem. Lett 9, 3116 (2018)



    int    i, j, k, l, m, n, tmp, prev_surf, old_surf, hopsurf, new_surf, iatom, istep, nstep, *ind, indmax, *indtmp, *indtmp1, *map, do_lambda_i_tmp;
    double **grad0, *ham, *old_ham, **ham_full, *new_ham, *adiapop_old, *adiapop, *initprob, *probtmp1, *probtmp2, *diapop, *diatmp;
    double tdse_timestep = 1.e-5*PS_TO_AU, diapopmax = 0.999;
    double rk_timestep0, probsum, random_number, cur_weight, random_number1;
    double popul_norm, overlapmax, overlaptmp1, overlaptmp2, adiapoptmp0, diapoptot, diapoptot0;
    double nc, kai1, kai2, kai, delta, aa, bb, cc, dd, ee, de;
    double decay_time, curr_surf_fac, trivial_cross_type;
    double dot_product, ekintot;
    double delta_min;


    old_ham  = ct->tfl_old_ham_adiab;
    new_ham  = ct->hamiltonian_adiab;
    ham      = ct->tfl_mean_ham;
    ham_full = ct->tfl_mean_ham_full;

    snew(diapop, ct->dim);
    snew(diatmp, 2*ct->dim);
    snew(adiapop_old, ct->dim);
    snew(adiapop, ct->dim);
    snew(initprob, ct->dim);
    snew(probtmp1, ct->dim);
    snew(probtmp2, ct->dim);

    snew(ind, ct->dim);
    snew(indtmp, ct->dim);
    snew(indtmp1, ct->dim);
    snew(map, ct->dim);

    if (ct->do_lambda_i == 3)
    {
       snew(grad0, dftb->phase2.nn);
       for (j = 0; j < dftb->phase2.nn; j++)
           snew(grad0[j], 3);
    }

// set up initial surface
    if (ct->tfs_initialization_step)
    {
       printf("Input wavefunction is:\n");
       for (j = 0; j < ct->dim; j++)
           printf("%f ", ct->wf[j]);
       printf("\n");

       for (i = 0; i < SQR(ct->dim); i++)
       {
           ct->tfl_old_ham[i] = *(ct->hamiltonian[0]+i);
           ham[i] = ct->tfl_old_ham[i];
       }

// diagonalize this Hamiltonian
       dsyev(ct->dim, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);

// record eigenvector
       for (i = 0; i < ct->dim; i++)
           for (j = 0; j < ct->dim; j++)
               ct->tfs_vector[i][j] = ham[i * ct->dim + j];

// determine inital surface using eigenvector as a criterion
       probsum = 0.;
       random_number = drand48();
       for (i = 0; i < ct->dim; i++)
       {
           initprob[i] = 0.0;

           for (j = 0; j < ct->dim; j++)
               initprob[i] += ct->tfs_vector[i][j]*ct->wf[j];

           initprob[i] = SQR(initprob[i]);

           probsum += initprob[i];
           if (random_number < probsum)
           {
               ct->surface = i;
               break;
           }
       }

       for (i = 0; i < 2*ct->dim; i++)
       {
           ct->tfs_popul[i] = 0.0; ct->tfs_diab[i]  = 0.0;
       }

       for (i = 0; i < ct->dim; i++)
           ct->tfs_diab[i] = ct->wf[i];

       ct->tfs_popul[ct->surface] = 1.0;

       ct->tfs_initialization_step = 0;

       return 0;
    }

// set up the old and new Hamiltonian for linear interpolation
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            old_ham[i*ct->dim+j] = ct->tfl_old_ham[i*ct->dim+j];
            new_ham[i*ct->dim+j] = ct->hamiltonian[i][j];
        }
    }

// store the vector at preceding time step
/*    for (i = 0; i < ct->dim; i++)
        for (j = 0; j < ct->dim; j++)
            ct->tfs_vector_old2[i][j] = ct->tfs_vector[i][j];
*/
// solve time-dependent Schroedinger equation
    nstep = (int)(ct->rk_timestep/tdse_timestep);
    rk_timestep0 = ct->rk_timestep;
    ct->rk_timestep = tdse_timestep;
    printf("Time step for TDSE: %f fs \n", ct->rk_timestep/PS_TO_AU*1000);

    for (istep = 1; istep < nstep+1; istep++)
    {
// linear interpolation for Hamiltonian
        cur_weight = (double)istep/(double)nstep;
        for (i = 0; i < ct->dim; i++)
            for (j = 0; j < ct->dim; j++)
                ham_full[i][j] = (1.0 - cur_weight)*ct->tfl_old_ham[i*ct->dim+j] + cur_weight*ct->hamiltonian[i][j];
// propagate the TDSE in the diabatic representation
        do_rksuite_diab(ct);
    }
    ct->rk_timestep = rk_timestep0;

// diagonalize this Hamiltonian
    for (i = 0; i < SQR(ct->dim); i++)
        ham[i] = new_ham[i];

    dsyev(ct->dim, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);

// store the eigenvectors at time t-dt and time t
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_vector_old[i][j] = ct->tfs_vector[i][j];
            ct->tfs_vector[i][j]     = ham[i*ct->dim + j];
        }
    }

// adiabatic wave function overlaps between time t-dt and t
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_overlap[i][j] = 0.0;

            for (k = 0; k < ct->dim; k++)
                ct->tfs_overlap[i][j] += ct->tfs_vector_old[i][k]*ct->tfs_vector[j][k];
        }
    }

// reassignment of states according to the maximum wave function overlaps
    for (i = 0; i < ct->dim; i++)
    {
        overlapmax = 0.0;

        for (j = 0; j < ct->dim; j++)
        {
            if (overlapmax < fabs(ct->tfs_overlap[j][i]))
            {
               overlapmax = fabs(ct->tfs_overlap[j][i]);
               indtmp[i] = j;
            }
        }
    }

// reassignment of states according to the min-cost algorithm
//    assignprob(ct, indtmp);

// state tracking
    for (i = 0; i < ct->dim; i++)
    {
        map[i] = 0; ind[i] = -1;
    }

    for (i = 0; i < ct->dim; i++)
    {
        if (fabs(ct->tfs_overlap[indtmp[i]][i]) > 0.9)
        {
           ind[i] = indtmp[i];
           map[ind[i]] = 1;
        }
    }

    for (i = 0; i < ct->dim; i++)
    {
        if (ind[i] == -1)
        {
           for (j = 0; j < ct->dim; j++)
           {
               if (! map[j])
               {
                  ind[i] = j; map[j] = 1;
                  break;
               }
           }
        }
    }

// check signs
    for (i = 0; i < ct->dim; i++)
    {
        j = ind[i];

        if (ct->tfs_overlap[j][i] < 0.0)
        {
            for (k = 0; k < ct->dim; k++)
            {
                ct->tfs_vector[i][k]  *= -1.0;
                ct->tfs_overlap[k][i] *= -1.0;
            }
        }
    }

// reassign the active state
    prev_surf = ct->surface;
    if (ct->surface != ind[ct->surface])
    {
       for (i = 0; i < ct->dim; i++)
       {
           if (ct->surface == ind[i])
           {
              printf("Changed from surface %d to surface %d \n", ct->surface, i);
              ct->surface = i;
              break;
           }
       }
    }

// calculate the adiabatic populations at the previous time step
    for (i = 0; i < ct->dim; i++)
        adiapop_old[i] = SQR(ct->tfs_popul[i]) + SQR(ct->tfs_popul[i+ct->dim]);

// hopping probabilities based on Tully's fewest-switches anzatz
    for (i = 0; i < ct->dim; i++)
    {
        j = ind[i];
//     nonadibatic couplings based on Hammes-Schiffer and Tully, JCP 101, 4657 (1994)
//        nc = (ct->tfs_overlap[prev_surf][i] - ct->tfs_overlap[j][ct->surface])/2.0/ct->rk_timestep;
//     ndibatic couplings based norm-preserving interpolation by Benjamin, JPCL 5, 2351 (2014)
        overlaptmp1 = sqrt(fabs(1.0 - SQR(ct->tfs_overlap[j][i]) - SQR(ct->tfs_overlap[prev_surf][i])));
        overlaptmp2 = (- ct->tfs_overlap[j][ct->surface]*ct->tfs_overlap[j][i] - ct->tfs_overlap[prev_surf][ct->surface]*ct->tfs_overlap[prev_surf][i])/overlaptmp1;

        aa = - sin(acos(ct->tfs_overlap[j][i]) - asin(ct->tfs_overlap[j][ct->surface]))/(acos(ct->tfs_overlap[j][i]) - asin(ct->tfs_overlap[j][ct->surface]));
        bb =   sin(acos(ct->tfs_overlap[j][i]) + asin(ct->tfs_overlap[j][ct->surface]))/(acos(ct->tfs_overlap[j][i]) + asin(ct->tfs_overlap[j][ct->surface]));
        cc =   sin(acos(ct->tfs_overlap[prev_surf][ct->surface]) - asin(ct->tfs_overlap[prev_surf][i]))/(acos(ct->tfs_overlap[prev_surf][ct->surface]) - asin(ct->tfs_overlap[prev_surf][i]));
        dd =   sin(acos(ct->tfs_overlap[prev_surf][ct->surface]) + asin(ct->tfs_overlap[prev_surf][i]))/(acos(ct->tfs_overlap[prev_surf][ct->surface]) + asin(ct->tfs_overlap[prev_surf][i]));
        ee = 2.0*asin(overlaptmp1*(overlaptmp1*overlaptmp2*asin(overlaptmp1) + (sqrt((1.0 - SQR(overlaptmp1))*(1.0 - SQR(overlaptmp2))) - 1.0)*asin(overlaptmp2)))/(SQR(asin(overlaptmp1)) - SQR(asin(overlaptmp2)));

        nc = (acos(ct->tfs_overlap[j][i])*(aa + bb) + asin(ct->tfs_overlap[prev_surf][i])*(cc + dd) + ee)/2.0/ct->rk_timestep;

        probtmp1[i] = (2.0*(ct->tfs_popul[prev_surf]*ct->tfs_popul[j]+ct->tfs_popul[ct->dim+prev_surf]*ct->tfs_popul[ct->dim+j])*nc)*ct->rk_timestep/adiapop_old[prev_surf];

        if (probtmp1[i] < 0.0 || i == ct->surface)
           probtmp1[i] = 0.0;
    }

// calculate expansion coefficients in the adiabatic basis
    for (i = 0; i < ct->dim; i++)
    {
        ct->tfs_popul[i] = 0.0; ct->tfs_popul[i+ct->dim] = 0.0;

        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_popul[i]         += ct->tfs_vector[i][j]*ct->tfs_diab[j];
            ct->tfs_popul[i+ct->dim] += ct->tfs_vector[i][j]*ct->tfs_diab[j+ct->dim];
        }
    }

// record the adiabatic populations at new time
    for (i = 0; i < ct->dim; i++)
        adiapop[i] = SQR(ct->tfs_popul[i])+SQR(ct->tfs_popul[i+ct->dim]);

// hopping probabilities based on populations flux
/*    for (i = 0; i < ct->dim; i++)
    {
        j = ind[i];

        probtmp2[i] = (adiapop[i] - adiapop_old[j])/adiapop_old[prev_surf];

        if (probtmp2[i] < 0.0 || i == ct->surface)
           probtmp2[i] = 0.0;
    }
*/
    new_surf = ct->surface;
    for (i = 0; i < ct->dim; i++)
    {
        if (prev_surf == indtmp[i])
        {
           new_surf = i;
           break;
        }
    }

// SC-FSSH, avoid the numerical instability of the NCs
    if (ct->surface != new_surf)
    {
       probtmp1[new_surf] = (adiapop_old[ct->surface] - adiapop[ct->surface])/adiapop_old[ct->surface];

       for (i = 0; i < ct->dim; i++)
       {
           if (i != new_surf)
           {
              probtmp1[new_surf] -= probtmp1[i];
           }
       }

       printf("SC is actived: %d  %d  %10f  %10f\n", ct->surface, new_surf, (adiapop_old[ct->surface] - adiapop[ct->surface])/adiapop_old[ct->surface], probtmp1[new_surf]);
    }

// define how to do adjustment for velocities
// 0: isotropic rescaling; 1: along NAC; 2: Boltzmann corrections.
    if (ct->jobtype == cteBCJFSSH)
    {
       dftb->do_ncv = 2;
// Boltzmann factor corrected hopping probability
       for (i = ct->surface+1; i < ct->dim; i++)
           probtmp1[i] *= exp(-(ct->ev_adiab[i] - ct->ev_adiab[ct->surface])/ct->fermi_kt);
    }
    else
    {
       dftb->do_ncv = 1;
    }

// print out the surface probabilities
    for (i = 0; i < ct->dim; i++)
        printf("%2d  %2d  %2d  %15.5f \n", i, indtmp[i], ind[i], probtmp1[i]);

// stochastic algorithm
    probsum = 0.0;
    random_number = drand48();
    for (i = 0; i < ct->dim; i++)
    {
        ct->surf_prob[i] = probtmp1[i];

        if (ct->surf_prob[i] > 0.0 && i != ct->surface)
        {
           probsum += ct->surf_prob[i];

           if (random_number < probsum)
           {
              hopsurf = i;

// calculate the nonadibatic coupling vector for velocity adjustment
              if (dftb->do_ncv == 1)
              {

                 ct->ncv_isurf = ct->surface; ct->ncv_jsurf = hopsurf;

                 if (ct->do_lambda_i == 3)
                 {
                    for (iatom = mdatoms->start; iatom < mdatoms->start+mdatoms->homenr; iatom++)
                    {
                        n = 0;
                        for (k = 0; k < ct->sites; k++)
                        {
                            for (l = 0; l < ct->site[k].atoms; l++)
                            {
                                if (iatom == ct->site[k].atom[l])
                                {
                                   for (m = 0; m < DIM; m++)
                                   {
                                       grad0[n][m] = dftb->phase2.grad[n][m];
                                   }
                                }
                                n++;
                            }
                        }
                    }
                 }

                 do_lambda_i_tmp = ct->do_lambda_i;
                 ct->do_lambda_i = 3;

                 if(ct->hamiltonian_type == 0)
                 {
                   get_MM_params(ct, dftb);
                 }
                 else if (ct->hamiltonian_type == 1)
                 {
                   get_nn_hamiltonian(ct,state_global, mdatoms, dftb);
                 }


                 ct->do_lambda_i = do_lambda_i_tmp;

                 aa = 0.0; bb = 0.0;
                 cc = ct->ev_adiab[ct->ncv_jsurf] - ct->ev_adiab[ct->ncv_isurf];

                 for (iatom = mdatoms->start; iatom < mdatoms->start+mdatoms->homenr; iatom++)
                 {
                     n = 0;
                     for (k = 0; k < ct->sites; k++)
                     {
                         for (l = 0; l < ct->site[k].atoms; l++)
                         {
                             if (iatom == ct->site[k].atom[l])
                             {
                                for (m = 0; m < DIM; m++)
                                {

                                    ct->tfs_ncv[k][l][m] = (dftb->phase1[k].grad[l][m] + dftb->phase2.grad[n][m])/cc *
                                                           (ct->is_hole_transfer ? 1.0 : -1.0);

                                    if (ct->do_lambda_i == 3)
                                    {
                                       dftb->phase2.grad[n][m] = grad0[n][m];
                                    }

                                    aa += SQR(ct->tfs_ncv[k][l][m])/2./dftb->phase1[k].mass[l]/AMU_TO_AU;
                                    bb += ct->tfs_ncv[k][l][m]*state_global->v[iatom][m]*NM_TO_BOHR/PS_TO_AU;

                                }
                             }
                             n++;
                         }
                     }
                 }

// adjust velocity along the direction of nonadiabtic coupling vector
                 delta = SQR(bb) - 4.*aa*cc;
                 if (delta >= 0.0)
                 {
                    kai1 = (- bb + sqrt(delta))/2./aa;
                    kai2 = (- bb - sqrt(delta))/2./aa;

                    kai = (fabs(kai1) < fabs(kai2) ? kai1 : kai2);

                    printf("Hop from %d to %d succeeds! Adjust velocity along NCV with kai = %8.5f \n", ct->surface, hopsurf, kai);
                 }
                 else
                 {
                    kai = - bb/aa;
                    hopsurf = ct->surface;

                    printf("Hop from %d to %d fails! Reverse velocity along NCV with kai = %8.5f \n", ct->surface, i, kai);
                 }

// velocity adjustment
                 for (iatom = mdatoms->start; iatom < mdatoms->start+mdatoms->homenr; iatom++)
                 {
                     for (k = 0; k < ct->sites; k++)
                     {
                         for (l = 0; l < ct->site[k].atoms; l++)
                         {
                             if (iatom == ct->site[k].atom[l])
                             {
                                for (m = 0; m < DIM; m++)
                                    state_global->v[iatom][m] += kai*ct->tfs_ncv[k][l][m]/dftb->phase1[k].mass[l]/AMU_TO_AU/NM_TO_BOHR*PS_TO_AU;
                             }
                         }
                     }
                 }
              }
              else if (dftb->do_ncv == 0)
              {
// compute kinetic energy
                 ekintot = 0.0;
                 for (iatom = mdatoms->start; iatom < mdatoms->start+mdatoms->homenr; iatom++)
                 {
                     for (k = 0; k < ct->sites; k++)
                     {
                         for (l = 0; l < ct->site[k].atoms; l++)
                         {
                             if (iatom == ct->site[k].atom[l])
                             {
                                for (m = 0; m < DIM; m++)
                                    ekintot += 0.5*dftb->phase1[k].mass[l]*AMU_TO_AU*SQR(state_global->v[iatom][m]*NM_TO_BOHR/PS_TO_AU);
                             }
                         }
                     }
                 }
                 printf("Kinetic energy is %f\n", ekintot);

                 de = ct->ev_adiab[i] - ct->ev_adiab[ct->surface];

                 if (1.0 - de/ekintot < 0.0)
                 {
                    printf("Frustated hop! \n");
                    fprintf(f,"%8.2f %10d %17d %21f %19f Frustrated \n",(step*ct->rk_timestep)*1000./PS_TO_AU, ct->surface, i, ct->surf_prob[i], random_number);
                    break;
                 }
                 else
                 {
                    kai = sqrt(1.0 - de/ekintot);
                    printf("kai is %f \n", kai);

                    ekintot -= de;

                    for (iatom = mdatoms->start; iatom < mdatoms->start+mdatoms->homenr; iatom++)
                    {
                        for (k = 0; k < ct->sites; k++)
                        {
                            for (l = 0; l < ct->site[k].atoms; l++)
                            {
                                if (iatom == ct->site[k].atom[l])
                                {
                                   for (m = 0; m < DIM; m++)
                                       state_global->v[iatom][m] *= kai;
                                }
                            }
                        }
                    }
                 }

              }

              if (hopsurf != ct->surface)
                 fprintf(f,"%8.2f %10d %17d %21f %19f \n",(step*ct->rk_timestep)*1000./PS_TO_AU, ct->surface, hopsurf, ct->surf_prob[i], random_number);

              ct->surface = hopsurf;

              break;
           }
        }
    }

    dftb->do_ncv = 0;

    printf("Current Surface: %d\n", ct->surface);

// compute kinetic energy
    ekintot = 0.0;
    for (iatom = mdatoms->start; iatom < mdatoms->start+mdatoms->homenr; iatom++)
    {
        for (k = 0; k < ct->sites; k++)
        {
            for (l = 0; l < ct->site[k].atoms; l++)
            {
                if (iatom == ct->site[k].atom[l])
                {
                   for (m = 0; m < DIM; m++)
                       ekintot += 0.5*dftb->phase1[k].mass[l]*AMU_TO_AU*SQR(state_global->v[iatom][m]*NM_TO_BOHR/PS_TO_AU);
                }
            }
        }
    }
    printf("Kinetic energy is %f\n", ekintot);

// decoherence correction
    dot_product = 0.0;
    for (i = 0; i < ct->dim; i++)
    {
        if (i != ct->surface)
        {
            decay_time = (1.0 + 0.1/ekintot)/fabs(ct->ev_adiab[i] - ct->ev_adiab[ct->surface]);

            ct->tfs_popul[i]         *= exp(-ct->rk_timestep/decay_time);
            ct->tfs_popul[i+ct->dim] *= exp(-ct->rk_timestep/decay_time);
            dot_product              += SQR(ct->tfs_popul[i])+SQR(ct->tfs_popul[i+ct->dim]);
        }
    }
    curr_surf_fac = sqrt((1.0 - dot_product)/(SQR(ct->tfs_popul[ct->surface]) + SQR(ct->tfs_popul[ct->dim+ct->surface])));
    ct->tfs_popul[ct->surface]         *= curr_surf_fac;
    ct->tfs_popul[ct->surface+ct->dim] *= curr_surf_fac;

// recalc. diab. from adiab after decoherence corr.
    for (i = 0; i < ct->dim; i++)
    {
        diatmp[i] = ct->tfs_diab[i]; diatmp[i+ct->dim] = ct->tfs_diab[i+ct->dim];

        ct->tfs_diab[i] = 0.0; ct->tfs_diab[i+ct->dim] = 0.0;

        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_diab[i]         += ct->tfs_popul[j]*ct->tfs_vector[j][i];
            ct->tfs_diab[i+ct->dim] += ct->tfs_popul[j+ct->dim]*ct->tfs_vector[j][i];
        }
    }

// johen's corrections
    for (i = 0; i < ct->dim; i++)
        diapop[i] = SQR(ct->tfs_diab[i])+SQR(ct->tfs_diab[i+ct->dim]);

    diapoptot = 0.0;
    for (i = 0; i < ct->dim; i++)
    {
        diapoptot += diapop[i];

        if (diapoptot > diapopmax)
        {
           n = i; diapoptot -= diapop[i];
           break;
        }
    }

    diapoptot0 = 0.0;
    for (i = n; i < ct->dim; i++)
    {
        ct->tfs_diab[i]         = diatmp[i];
        ct->tfs_diab[i+ct->dim] = diatmp[i+ct->dim];

        diapoptot0 += SQR(ct->tfs_diab[i])+SQR(ct->tfs_diab[i+ct->dim]);
    }

    for (i = 0; i < n; i++)
    {
        ct->tfs_diab[i]         *= sqrt((1.0 - diapoptot0)/diapoptot);
        ct->tfs_diab[i+ct->dim] *= sqrt((1.0 - diapoptot0)/diapoptot);
    }

// calculate expansion coefficients in the adiabatic basis
    for (i = 0; i < ct->dim; i++)
    {
        ct->tfs_popul[i] = 0.0; ct->tfs_popul[i+ct->dim] = 0.0;

        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_popul[i]         += ct->tfs_vector[i][j]*ct->tfs_diab[j];
            ct->tfs_popul[i+ct->dim] += ct->tfs_vector[i][j]*ct->tfs_diab[j+ct->dim];
        }
    }

// charge wave function
    for (i = 0; i < ct->dim; i++)
        ct->wf[i] = ct->tfs_vector[ct->surface][i];

//store hamiltonian as old_ham
    for (i = 0; i < SQR(ct->dim); i++)
        ct->tfl_old_ham[i] = *(ct->hamiltonian[0]+i);

    sfree(diapop);
    sfree(diatmp);
    sfree(adiapop_old);
    sfree(adiapop);
    sfree(initprob);
    sfree(probtmp1);
    sfree(probtmp2);

    sfree(ind);
    sfree(indtmp);
    sfree(indtmp1);
    sfree(map);

    if (ct->do_lambda_i == 3)
    {
       for (i = 0; i < dftb->phase2.nn; i++)
           sfree(grad0[i]);
    }

    return 0;
}

int do_self_consistent_restricted_decoherence_fssh(int step, charge_transfer_t *ct, dftb_t *dftb, t_state *state_global, t_mdatoms *mdatoms, FILE *f)
{
// Declaration:
//
// 1. self-consistent algorithm is implemented for weak avoided crossings;
// 2. restricted decoherence correction only for the adiabatic population larger than a threshold;
//
//                                          -----------Weiwei Xie on Mar. 2019
// Ref: J.Chem.Phys. 148, 104106 (2018)


    int    i, j, k, l, m, n, tmp, prev_surf, old_surf, hopsurf, new_surf, surftmp, iatom, istep, nstep, *ind, indmax, *indtmp, *indtmp1, *map, newline;
    double *ham, *old_ham, **ham_full, *new_ham, *adiapop_old, *adiapop, *initprob, *probtmp1, *probtmp2, *diapop, *diatmp;
    double tdse_timestep = 1.e-5*PS_TO_AU;
    double rk_timestep0, probsum, random_number, cur_weight;
    double popul_norm, overlapmax, overlaptmp1, overlaptmp2, adiapoptmp0, probtmp0, probtmp, diapoptmp;
    double nc, kai1, kai2, kai, delta, aa, bb, cc, dd, ee;
    double decay_time, curr_surf_fac, trivial_cross_type;
    double dot_product, ekintot;
    double delta_min;


    old_ham  = ct->tfl_old_ham_adiab;
    new_ham  = ct->hamiltonian_adiab;
    ham      = ct->tfl_mean_ham;
    ham_full = ct->tfl_mean_ham_full;

    snew(diapop, ct->dim);
    snew(diatmp, 2*ct->dim);
    snew(adiapop_old, ct->dim);
    snew(adiapop, ct->dim);
    snew(initprob, ct->dim);
    snew(probtmp1, ct->dim);
    snew(probtmp2, ct->dim);

    snew(ind, ct->dim);
    snew(indtmp, ct->dim);
    snew(indtmp1, ct->dim);
    snew(map, ct->dim);

// set up initial surface
    if (ct->tfs_initialization_step)
    {
       printf("Input wavefunction is:\n");
       for (j = 0; j < ct->dim; j++)
           printf("%f ", ct->wf[j]);
       printf("\n");

       for (i = 0; i < SQR(ct->dim); i++)
       {
           ct->tfl_old_ham[i] = *(ct->hamiltonian[0]+i);
           ham[i] = ct->tfl_old_ham[i];
       }

// diagonalize this Hamiltonian
       dsyev(ct->dim, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);

// record eigenvector
       for (i = 0; i < ct->dim; i++)
           for (j = 0; j < ct->dim; j++)
               ct->tfs_vector[i][j] = ham[i * ct->dim + j];

// determine inital surface using eigenvector as a criterion
       probsum = 0.;
       random_number = drand48();
       for (i = 0; i < ct->dim; i++)
       {
           initprob[i] = 0.;
           for (j = 0; j < ct->dim; j++)
               initprob[i] += ct->tfs_vector[i][j]*ct->wf[j];
           initprob[i] = SQR(initprob[i]);

           probsum += initprob[i];
           if (random_number < probsum)
           {
               ct->surface = i;
               break;
           }
       }

       for (i = 0; i < 2*ct->dim; i++)
       {
           ct->tfs_popul[i] = 0.0; ct->tfs_diab[i]  = 0.0;
       }

       for (i = 0; i < ct->dim; i++)
           ct->tfs_diab[i] = ct->wf[i];

       ct->tfs_popul[ct->surface] = 1.0;

       ct->tfs_initialization_step = 0;

       return 0;
    }

// set up the old and new Hamiltonian for linear interpolation
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            old_ham[i*ct->dim+j] = ct->tfl_old_ham[i*ct->dim+j];
            new_ham[i*ct->dim+j] = ct->hamiltonian[i][j];
        }
    }

// store the vector at preceding time step
/*    for (i = 0; i < ct->dim; i++)
        for (j = 0; j < ct->dim; j++)
            ct->tfs_vector_old2[i][j] = ct->tfs_vector[i][j];
*/
// solve time-dependent Schroedinger equation
    nstep = (int)(ct->rk_timestep/tdse_timestep);
    rk_timestep0 = ct->rk_timestep;
    ct->rk_timestep = tdse_timestep;
    printf("Time step for TDSE: %f fs \n", ct->rk_timestep/PS_TO_AU*1000);

    for (istep = 1; istep < nstep+1; istep++)
    {
// linear interpolation for Hamiltonian
        cur_weight = (double)istep/(double)nstep;
        for (i = 0; i < ct->dim; i++)
            for (j = 0; j < ct->dim; j++)
                ham_full[i][j] = (1.0 - cur_weight)*ct->tfl_old_ham[i*ct->dim+j] + cur_weight*ct->hamiltonian[i][j];
// propagate the TDSE in the diabatic representation
        do_rksuite_diab(ct);
    }
    ct->rk_timestep = rk_timestep0;

// diagonalize this Hamiltonian
    for (i = 0; i < SQR(ct->dim); i++)
        ham[i] = new_ham[i];

    dsyev(ct->dim, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);

// store the eigenvectors at time t-dt and time t
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_vector_old[i][j] = ct->tfs_vector[i][j];
            ct->tfs_vector[i][j]     = ham[i*ct->dim + j];
        }
    }

// adiabatic wave function overlaps between time t-dt and t
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_overlap[i][j] = 0.0;

            for (k = 0; k < ct->dim; k++)
                ct->tfs_overlap[i][j] += ct->tfs_vector_old[i][k]*ct->tfs_vector[j][k];
        }
    }

// reassignment of states according to the maximum wave function overlaps
/*    for (i = 0; i < ct->dim; i++)
    {
        overlapmax = 0.0;
        for (j = 0; j < ct->dim; j++)
        {
            if (overlapmax < fabs(ct->tfs_overlap[j][i]))
            {
               overlapmax = fabs(ct->tfs_overlap[j][i]);
               indtmp[i] = j;
            }
        }
    }
*/
// check signs
    for (i = 0; i < ct->dim; i++)
    {
        if (ct->tfs_overlap[i][i] < 0.0)
        {
            for (k = 0; k < ct->dim; k++)
            {
                ct->tfs_vector[i][k]  *= -1.0;
                ct->tfs_overlap[k][i] *= -1.0;
            }
        }
    }

// reassignment of states according to the min-cost algorithm
    assignprob(ct, ind);

// calculate the adiabatic populations at the previous time step
    for (i = 0; i < ct->dim; i++)
        adiapop_old[i] = SQR(ct->tfs_popul[i]) + SQR(ct->tfs_popul[i+ct->dim]);

// hopping probabilities based on Tully's fewest-switches anzatz
    for (i = 0; i < ct->dim; i++)
    {
//     nonadibatic couplings based on Hammes-Schiffer and Tully, JCP 101, 4657 (1994)
//        nc = (ct->tfs_overlap[prev_surf][i] - ct->tfs_overlap[j][ct->surface])/2.0/ct->rk_timestep;
//     ndibatic couplings based norm-preserving interpolation by Benjamin, JPCL 5, 2351 (2014)
        overlaptmp1 = sqrt(fabs(1.0 - SQR(ct->tfs_overlap[i][i]) - SQR(ct->tfs_overlap[ct->surface][i])));
        overlaptmp2 = (- ct->tfs_overlap[i][ct->surface]*ct->tfs_overlap[i][i] - ct->tfs_overlap[ct->surface][ct->surface]*ct->tfs_overlap[ct->surface][i])/overlaptmp1;

        aa = - sin(acos(ct->tfs_overlap[i][i]) - asin(ct->tfs_overlap[i][ct->surface]))/(acos(ct->tfs_overlap[i][i]) - asin(ct->tfs_overlap[i][ct->surface]));
        bb =   sin(acos(ct->tfs_overlap[i][i]) + asin(ct->tfs_overlap[i][ct->surface]))/(acos(ct->tfs_overlap[i][i]) + asin(ct->tfs_overlap[i][ct->surface]));
        cc =   sin(acos(ct->tfs_overlap[ct->surface][ct->surface]) - asin(ct->tfs_overlap[ct->surface][i]))/(acos(ct->tfs_overlap[ct->surface][ct->surface]) - asin(ct->tfs_overlap[ct->surface][i]));
        dd =   sin(acos(ct->tfs_overlap[ct->surface][ct->surface]) + asin(ct->tfs_overlap[ct->surface][i]))/(acos(ct->tfs_overlap[ct->surface][ct->surface]) + asin(ct->tfs_overlap[ct->surface][i]));
        ee = 2.0*asin(overlaptmp1*(overlaptmp1*overlaptmp2*asin(overlaptmp1) + (sqrt((1.0 - SQR(overlaptmp1))*(1.0 - SQR(overlaptmp2))) - 1.0)*asin(overlaptmp2)))/(SQR(asin(overlaptmp1)) - SQR(asin(overlaptmp2)));

        nc = (acos(ct->tfs_overlap[i][i])*(aa + bb) + asin(ct->tfs_overlap[ct->surface][i])*(cc + dd) + ee)/2.0/ct->rk_timestep;

        probtmp1[i] = (2.0*(ct->tfs_popul[ct->surface]*ct->tfs_popul[i]+ct->tfs_popul[ct->dim+ct->surface]*ct->tfs_popul[ct->dim+i])*nc)*ct->rk_timestep/adiapop_old[ct->surface];

        if (probtmp1[i] < 0.0 || i == ct->surface)
           probtmp1[i] = 0.0;
    }

// calculate expansion coefficients in the adiabatic basis
    for (i = 0; i < ct->dim; i++)
    {
        ct->tfs_popul[i] = 0.0; ct->tfs_popul[i+ct->dim] = 0.0;

        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_popul[i]         += ct->tfs_vector[i][j]*ct->tfs_diab[j];
            ct->tfs_popul[i+ct->dim] += ct->tfs_vector[i][j]*ct->tfs_diab[j+ct->dim];
        }
    }

// record the adiabatic populations at new time
    for (i = 0; i < ct->dim; i++)
        adiapop[i] = SQR(ct->tfs_popul[i])+SQR(ct->tfs_popul[i+ct->dim]);

// SC-FSSH, avoid the numerical instability of the NCs
    prev_surf = ct->surface;
    for (i = 0; i < ct->dim; i++)
    {
        if (prev_surf == ind[i])
        {
           new_surf = i;
           break;
        }
    }

    newline = 0;
    if (ct->surface != new_surf)
    {
       probtmp1[new_surf] = (adiapop_old[ct->surface] - adiapop[ct->surface])/adiapop_old[ct->surface];

       for (i = 0; i < ct->dim; i++)
       {
           if (i != new_surf)
           {
              probtmp1[new_surf] -= probtmp1[i];
           }
       }

//       printf("SC is actived: %d  %d  %10f  %10f\n", ct->surface, new_surf, (adiapop_old[ct->surface] - adiapop[ct->surface])/adiapop_old[ct->surface], probtmp1[new_surf]);
       fprintf(f,"%8.2f %10d %17d %21f %19f %19f    UnavoidedCrossing \n",(step*ct->rk_timestep)*1000./PS_TO_AU, ct->surface, new_surf,
              probtmp1[new_surf], fabs(ct->tfs_overlap[ct->surface][new_surf]), (ct->ev_adiab[new_surf] - ct->ev_adiab[ct->surface])*HARTREE_TO_EV);
       newline = 1;
    }

    for (i = 0; i < ct->dim; i++)
        printf("%2d  %2d  %15.5f \n", i, ind[i], probtmp1[i]);

// stochastic algorithm
    probsum = 0.0;
    random_number = drand48();
    for (i = 0; i < ct->dim; i++)
    {
        ct->surf_prob[i] = probtmp1[i];

        if (ct->surf_prob[i] >= 0.0 && i != ct->surface)
        {
           probsum += ct->surf_prob[i];

           if (random_number < probsum)
           {
              hopsurf = i;
// calculate the nonadibatic coupling vector for velocity adjustment
              dftb->do_ncv = 1;

              if (dftb->do_ncv && hopsurf != new_surf)
              {
                 ct->ncv_isurf = ct->surface; ct->ncv_jsurf = hopsurf;

                 get_MM_params(ct, dftb);

                 aa = 0.0; bb = 0.0;
                 cc = ct->ev_adiab[hopsurf] - ct->ev_adiab[ct->surface];

                 for (iatom = mdatoms->start; iatom < mdatoms->start+mdatoms->homenr; iatom++)
                 {
                     n = 0;
                     for (k = 0; k < ct->sites; k++)
                     {
                         for (l = 0; l < ct->site[k].atoms; l++)
                         {
                             if (iatom == ct->site[k].atom[l])
                             {
                                for (m = 0; m < DIM; m++)
                                {
                                    ct->tfs_ncv[k][l][m] = (real) (dftb->phase1[k].grad[l][m] + dftb->phase2.grad[n][m])/cc;

                                    aa += SQR(ct->tfs_ncv[k][l][m])/2./dftb->phase1[k].mass[l]/AMU_TO_AU;
                                    bb += ct->tfs_ncv[k][l][m]*state_global->v[iatom][m]*NM_TO_BOHR/PS_TO_AU;

                                }
                             }
                             n++;
                         }
                     }
                 }

// adjust velocity along the direction of nonadiabtic coupling vector
                 delta = SQR(bb) - 4.*aa*cc;
                 if (delta >= 0.0)
                 {
                    kai1 = (- bb + sqrt(delta))/2./aa;
                    kai2 = (- bb - sqrt(delta))/2./aa;

                    kai = (fabs(kai1) < fabs(kai2) ? kai1 : kai2);

                    fprintf(f,"Hop from %d to %d succeeds! Adjust velocity along NCV with kai = %8.5f \n", ct->surface, hopsurf, kai);
                 }
                 else
                 {
                    kai = - bb/aa;
                    hopsurf = ct->surface;

                    fprintf(f,"Hop from %d to %d fails! Reverse velocity along NCV with kai = %8.5f \n", ct->surface, i, kai);
                 }

// velocity adjustment
                 for (iatom = mdatoms->start; iatom < mdatoms->start+mdatoms->homenr; iatom++)
                 {
                     for (k = 0; k < ct->sites; k++)
                     {
                         for (l = 0; l < ct->site[k].atoms; l++)
                         {
                             if (iatom == ct->site[k].atom[l])
                             {
                                for (m = 0; m < DIM; m++)
                                    state_global->v[iatom][m] += kai*ct->tfs_ncv[k][l][m]/dftb->phase1[k].mass[l]/AMU_TO_AU/NM_TO_BOHR*PS_TO_AU;
                             }
                         }
                     }
                 }

                 dftb->do_ncv = 0;
              }

              if (hopsurf != ct->surface)
              {
                 fprintf(f,"%8.2f %10d %17d %21f %19f %19f \n",(step*ct->rk_timestep)*1000./PS_TO_AU, ct->surface, hopsurf,
                        ct->surf_prob[i], fabs(ct->tfs_overlap[ct->surface][hopsurf]), (ct->ev_adiab[hopsurf] - ct->ev_adiab[ct->surface])*HARTREE_TO_EV);
                 newline = 1;
              }

              ct->surface = hopsurf;

              break;
           }
        }
    }

   if (newline)
      fprintf(f, "\n");

    printf("Current Surface: %d\n", ct->surface);

    ekintot = 0.0;
    for (iatom = mdatoms->start; iatom < mdatoms->start+mdatoms->homenr; iatom++)
    {
        for (k = 0; k < ct->sites; k++)
        {
            for (l = 0; l < ct->site[k].atoms; l++)
            {
                if (iatom == ct->site[k].atom[l])
                {
                   for (m = 0; m < DIM; m++)
                       ekintot += 0.5*dftb->phase1[k].mass[l]*AMU_TO_AU*SQR(state_global->v[iatom][m]*NM_TO_BOHR/PS_TO_AU);
                }
            }
        }
    }
    printf("Kinetic energy is %f\n", ekintot);

// restricted decoherence correction
    dot_product = 0.0;
    if (adiapop[ct->surface] > 0.001)
    {
       for (i = 0; i < ct->dim; i++)
       {
           if (i != ct->surface)
           {
               decay_time = (1.0 + 0.1/ekintot)/fabs(ct->ev_adiab[i] - ct->ev_adiab[ct->surface]);

               ct->tfs_popul[i]         *= exp(-ct->rk_timestep/decay_time);
               ct->tfs_popul[i+ct->dim] *= exp(-ct->rk_timestep/decay_time);
               dot_product              += SQR(ct->tfs_popul[i])+SQR(ct->tfs_popul[i+ct->dim]);
           }
       }
       curr_surf_fac = sqrt((1.0 - dot_product)/(SQR(ct->tfs_popul[ct->surface]) + SQR(ct->tfs_popul[ct->dim+ct->surface])));
       ct->tfs_popul[ct->surface]         *= curr_surf_fac;
       ct->tfs_popul[ct->surface+ct->dim] *= curr_surf_fac;

// recalc. diab. from adiab after decoherence corr.
       for (i = 0; i < ct->dim; i++)
       {
           ct->tfs_diab[i] = 0.0; ct->tfs_diab[i+ct->dim] = 0.0;

           for (j = 0; j < ct->dim; j++)
           {
               ct->tfs_diab[i]         += ct->tfs_popul[j]*ct->tfs_vector[j][i];
               ct->tfs_diab[i+ct->dim] += ct->tfs_popul[j+ct->dim]*ct->tfs_vector[j][i];
           }
       }
    }

// charge wave function
    for (i = 0; i < ct->dim; i++)
        ct->wf[i] = ct->tfs_vector[ct->surface][i];

//store hamiltonian as old_ham
    for (i = 0; i < SQR(ct->dim); i++)
        ct->tfl_old_ham[i] = *(ct->hamiltonian[0]+i);

    return 0;
}

int do_crossing_corrected_fssh(int step, charge_transfer_t *ct, dftb_t *dftb, t_state *state_global, t_mdatoms *mdatoms, FILE *f)
{
// Declaration:
//
// 1. self-consitent algorithm for the weak avoided crossings is used;
// 2. different crossing types are detected and treated differently;
//
//                                          -----------Weiwei Xie on March. 2019
// Ref: J. Phys. Chem. Lett. 915, 4319 (2018)


    int    i, j, k, l, m, n, tmp, prev_surf, old_surf, hopsurf, new_surf, surftmp, iatom, istep, nstep, *ind, indmax, *indtmp, *indtmp1, *map, newline;
    double *ham, *old_ham, **ham_full, *new_ham, *adiapop_old, *adiapop, *initprob, *probtmp1, *probtmp2, *diapop, *diatmp;
    double tdse_timestep = 1.e-5*PS_TO_AU;
    double rk_timestep0, probsum, random_number, cur_weight;
    double popul_norm, overlapmax, overlaptmp1, overlaptmp2, adiapoptmp0, probtmp0, probtmp, diapoptmp;
    double nc, kai1, kai2, kai, delta, aa, bb, cc, dd, ee;
    double decay_time, curr_surf_fac, trivial_cross_type;
    double dot_product, ekintot;
    double delta_min;


    old_ham  = ct->tfl_old_ham_adiab;
    new_ham  = ct->hamiltonian_adiab;
    ham      = ct->tfl_mean_ham;
    ham_full = ct->tfl_mean_ham_full;

    snew(diapop, ct->dim);
    snew(diatmp, 2*ct->dim);
    snew(adiapop_old, ct->dim);
    snew(adiapop, ct->dim);
    snew(initprob, ct->dim);
    snew(probtmp1, ct->dim);
    snew(probtmp2, ct->dim);

    snew(ind, ct->dim);
    snew(indtmp, ct->dim);
    snew(indtmp1, ct->dim);
    snew(map, ct->dim);

// set up initial surface
    if (ct->tfs_initialization_step)
    {
       printf("Input wavefunction is:\n");
       for (j = 0; j < ct->dim; j++)
           printf("%f ", ct->wf[j]);
       printf("\n");

       for (i = 0; i < SQR(ct->dim); i++)
       {
           ct->tfl_old_ham[i] = *(ct->hamiltonian[0]+i);
           ham[i] = ct->tfl_old_ham[i];
       }

// diagonalize this Hamiltonian
       dsyev(ct->dim, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);

// record eigenvector
       for (i = 0; i < ct->dim; i++)
           for (j = 0; j < ct->dim; j++)
               ct->tfs_vector[i][j] = ham[i * ct->dim + j];

// determine inital surface using eigenvector as a criterion
       probsum = 0.;
       random_number = drand48();
       for (i = 0; i < ct->dim; i++)
       {
           initprob[i] = 0.;
           for (j = 0; j < ct->dim; j++)
               initprob[i] += ct->tfs_vector[i][j]*ct->wf[j];
           initprob[i] = SQR(initprob[i]);

           probsum += initprob[i];
           if (random_number < probsum)
           {
               ct->surface = i;
               break;
           }
       }

       for (i = 0; i < 2*ct->dim; i++)
       {
           ct->tfs_popul[i] = 0.0; ct->tfs_diab[i]  = 0.0;
       }

       for (i = 0; i < ct->dim; i++)
           ct->tfs_diab[i] = ct->wf[i];

       ct->tfs_popul[ct->surface] = 1.0;

       ct->tfs_initialization_step = 0;

       return 0;
    }

// set up the old and new Hamiltonian for linear interpolation
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            old_ham[i*ct->dim+j] = ct->tfl_old_ham[i*ct->dim+j];
            new_ham[i*ct->dim+j] = ct->hamiltonian[i][j];
        }
    }

// store the vector at preceding time step
/*    for (i = 0; i < ct->dim; i++)
        for (j = 0; j < ct->dim; j++)
            ct->tfs_vector_old2[i][j] = ct->tfs_vector[i][j];
*/
// solve time-dependent Schroedinger equation
    nstep = (int)(ct->rk_timestep/tdse_timestep);
    rk_timestep0 = ct->rk_timestep;
    ct->rk_timestep = tdse_timestep;
    printf("Time step for TDSE: %f fs \n", ct->rk_timestep/PS_TO_AU*1000);

    for (istep = 1; istep < nstep+1; istep++)
    {
// linear interpolation for Hamiltonian
        cur_weight = (double)istep/(double)nstep;
        for (i = 0; i < ct->dim; i++)
            for (j = 0; j < ct->dim; j++)
                ham_full[i][j] = (1.0 - cur_weight)*ct->tfl_old_ham[i*ct->dim+j] + cur_weight*ct->hamiltonian[i][j];
// propagate the TDSE in the diabatic representation
        do_rksuite_diab(ct);
    }
    ct->rk_timestep = rk_timestep0;

// diagonalize this Hamiltonian
    for (i = 0; i < SQR(ct->dim); i++)
        ham[i] = new_ham[i];

    dsyev(ct->dim, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);

// store the eigenvectors at time t-dt and time t
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_vector_old[i][j] = ct->tfs_vector[i][j];
            ct->tfs_vector[i][j]     = ham[i*ct->dim + j];
        }
    }

// adiabatic wave function overlaps between time t-dt and t
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_overlap[i][j] = 0.0;

            for (k = 0; k < ct->dim; k++)
                ct->tfs_overlap[i][j] += ct->tfs_vector_old[i][k]*ct->tfs_vector[j][k];
        }
    }

// check signs
    for (i = 0; i < ct->dim; i++)
    {
        if (ct->tfs_overlap[i][i] < 0.0)
        {
            for (k = 0; k < ct->dim; k++)
            {
                ct->tfs_vector[i][k]  *= -1.0;
                ct->tfs_overlap[k][i] *= -1.0;
            }
        }
    }

// reassignment of states according to the min-cost algorithm
    assignprob(ct, ind);

// reassignment of states according to the maximum wave function overlaps
/*    for (i = 0; i < ct->dim; i++)
    {
        overlapmax = 0.0;
        for (j = 0; j < ct->dim; j++)
        {
            if (overlapmax < fabs(ct->tfs_overlap[j][i]))
            {
               overlapmax = fabs(ct->tfs_overlap[j][i]);
               ind[i] = j;
            }
        }
    }
*/
// calculate the adiabatic populations at the previous time step
    for (i = 0; i < ct->dim; i++)
        adiapop_old[i] = SQR(ct->tfs_popul[i]) + SQR(ct->tfs_popul[i+ct->dim]);

// hopping probabilities based on Tully's fewest-switches anzatz
    for (i = 0; i < ct->dim; i++)
    {
//     nonadibatic couplings based on Hammes-Schiffer and Tully, JCP 101, 4657 (1994)
//        nc = (ct->tfs_overlap[prev_surf][i] - ct->tfs_overlap[j][ct->surface])/2.0/ct->rk_timestep;
//     ndibatic couplings based norm-preserving interpolation by Benjamin, JPCL 5, 2351 (2014)
        overlaptmp1 = sqrt(fabs(1.0 - SQR(ct->tfs_overlap[i][i]) - SQR(ct->tfs_overlap[ct->surface][i])));
        overlaptmp2 = (- ct->tfs_overlap[i][ct->surface]*ct->tfs_overlap[i][i] - ct->tfs_overlap[ct->surface][ct->surface]*ct->tfs_overlap[ct->surface][i])/overlaptmp1;

        aa = - sin(acos(ct->tfs_overlap[i][i]) - asin(ct->tfs_overlap[i][ct->surface]))/(acos(ct->tfs_overlap[i][i]) - asin(ct->tfs_overlap[i][ct->surface]));
        bb =   sin(acos(ct->tfs_overlap[i][i]) + asin(ct->tfs_overlap[i][ct->surface]))/(acos(ct->tfs_overlap[i][i]) + asin(ct->tfs_overlap[i][ct->surface]));
        cc =   sin(acos(ct->tfs_overlap[ct->surface][ct->surface]) - asin(ct->tfs_overlap[ct->surface][i]))/(acos(ct->tfs_overlap[ct->surface][ct->surface]) - asin(ct->tfs_overlap[ct->surface][i]));
        dd =   sin(acos(ct->tfs_overlap[ct->surface][ct->surface]) + asin(ct->tfs_overlap[ct->surface][i]))/(acos(ct->tfs_overlap[ct->surface][ct->surface]) + asin(ct->tfs_overlap[ct->surface][i]));
        ee = 2.0*asin(overlaptmp1*(overlaptmp1*overlaptmp2*asin(overlaptmp1) + (sqrt((1.0 - SQR(overlaptmp1))*(1.0 - SQR(overlaptmp2))) - 1.0)*asin(overlaptmp2)))/(SQR(asin(overlaptmp1)) - SQR(asin(overlaptmp2)));

        nc = (acos(ct->tfs_overlap[i][i])*(aa + bb) + asin(ct->tfs_overlap[ct->surface][i])*(cc + dd) + ee)/2.0/ct->rk_timestep;

        probtmp1[i] = (2.0*(ct->tfs_popul[ct->surface]*ct->tfs_popul[i]+ct->tfs_popul[ct->dim+ct->surface]*ct->tfs_popul[ct->dim+i])*nc)*ct->rk_timestep/adiapop_old[ct->surface];

        if (probtmp1[i] < 0.0 || i == ct->surface)
           probtmp1[i] = 0.0;
    }

// calculate expansion coefficients in the adiabatic basis
    for (i = 0; i < ct->dim; i++)
    {
        ct->tfs_popul[i] = 0.0; ct->tfs_popul[i+ct->dim] = 0.0;

        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_popul[i]         += ct->tfs_vector[i][j]*ct->tfs_diab[j];
            ct->tfs_popul[i+ct->dim] += ct->tfs_vector[i][j]*ct->tfs_diab[j+ct->dim];
        }
    }

// record the adiabatic populations at new time
    for (i = 0; i < ct->dim; i++)
        adiapop[i] = SQR(ct->tfs_popul[i])+SQR(ct->tfs_popul[i+ct->dim]);

// SC-FSSH, avoid the numerical instability of the NCs
    prev_surf = ct->surface;
    for (i = 0; i < ct->dim; i++)
    {
        if (prev_surf == ind[i])
        {
           new_surf = i;
           break;
        }
    }

    newline = 0;
    if (ct->surface != new_surf)
    {
       probtmp1[new_surf] = (adiapop_old[ct->surface] - adiapop[ct->surface])/adiapop_old[ct->surface];

       for (i = 0; i < ct->dim; i++)
       {
           if (i != new_surf)
           {
              probtmp1[new_surf] -= probtmp1[i];
           }
       }

//       printf("SC is actived: %d  %d  %10f  %10f\n", ct->surface, new_surf, (adiapop_old[ct->surface] - adiapop[ct->surface])/adiapop_old[ct->surface], probtmp1[new_surf]);
       fprintf(f,"%8.2f %10d %17d %21f %19f %19f   UnavoidedCrossing \n",(step*ct->rk_timestep)*1000./PS_TO_AU, ct->surface, new_surf,
              probtmp1[new_surf], fabs(ct->tfs_overlap[ct->surface][new_surf]), (ct->ev_adiab[new_surf] - ct->ev_adiab[ct->surface])*HARTREE_TO_EV);
       newline = 1;
    }

    for (i = 0; i < ct->dim; i++)
        printf("%2d  %2d  %15.5f \n", i, ind[i], probtmp1[i]);

// stochastic algorithm
    probsum = 0.0;
    random_number = drand48();
    for (i = 0; i < ct->dim; i++)
    {
        ct->surf_prob[i] = probtmp1[i];

        if (ct->surf_prob[i] >= 0.0 && i != ct->surface)
        {
           probsum += ct->surf_prob[i];

           if (random_number < probsum)
           {
              hopsurf = i;
// calculate the nonadibatic coupling vector for velocity adjustment
              dftb->do_ncv = 1;

              if (dftb->do_ncv)
              {
                 ct->ncv_isurf = ct->surface; ct->ncv_jsurf = hopsurf;

                 get_MM_params(ct, dftb);

                 aa = 0.0; bb = 0.0;
                 cc = ct->ev_adiab[hopsurf] - ct->ev_adiab[ct->surface];

                 for (iatom = mdatoms->start; iatom < mdatoms->start+mdatoms->homenr; iatom++)
                 {
                     n = 0;
                     for (k = 0; k < ct->sites; k++)
                     {
                         for (l = 0; l < ct->site[k].atoms; l++)
                         {
                             if (iatom == ct->site[k].atom[l])
                             {
                                for (m = 0; m < DIM; m++)
                                {
                                    ct->tfs_ncv[k][l][m] = (real) (dftb->phase1[k].grad[l][m] + dftb->phase2.grad[n][m])/cc;

                                    aa += SQR(ct->tfs_ncv[k][l][m])/2./dftb->phase1[k].mass[l]/AMU_TO_AU;
                                    bb += ct->tfs_ncv[k][l][m]*state_global->v[iatom][m]*NM_TO_BOHR/PS_TO_AU;

                                }
                             }
                             n++;
                         }
                     }
                 }

// crossing-corrected FSSH
                 for (j = 0; j < ct->dim; j++)
                 {
                     if (hopsurf == ind[j])
                     {
                        surftmp = j;
                        break;
                     }
                 }

                 delta = SQR(bb) - 4.0*aa*cc;

                 if (ct->surface == new_surf && delta >= 0.0)
                 {
                    kai1 = (- bb + sqrt(delta))/2./aa;
                    kai2 = (- bb - sqrt(delta))/2./aa;

                    kai = (fabs(kai1) < fabs(kai2) ? kai1 : kai2);

                    hopsurf = surftmp;

                    fprintf(f, "Hop from %d to %d succeeds! Adjust velocity along NCV with kai = %8.5f \n", ct->surface, hopsurf, kai);
                 }
                 else if (ct->surface == new_surf && delta < 0.0)
                 {
                    kai = - bb/aa;
                    hopsurf = ct->surface;

                    fprintf(f, "Hop from %d to %d fails! Reverse velocity along NCV with kai = %8.5f \n", ct->surface, i, kai);
                 }
                 else if (ct->surface != new_surf && hopsurf == new_surf)
                 {
                    kai = 0.0;
                    hopsurf = new_surf;

                    fprintf(f, "Hop from %d to %d (new) succeeds! \n", ct->surface, hopsurf);
                 }
                 else if (ct->surface != new_surf && hopsurf != new_surf && delta >= 0.0)
                 {
                    kai1 = (- bb + sqrt(delta))/2./aa;
                    kai2 = (- bb - sqrt(delta))/2./aa;

                    kai = (fabs(kai1) < fabs(kai2) ? kai1 : kai2);

                    hopsurf = surftmp;

                    fprintf(f, "Hop from %d to %d (new) succeeds! Adjust velocity along NCV with kai = %8.5f \n", ct->surface, hopsurf, kai);
                 }
                 else if (ct->surface != new_surf && hopsurf != new_surf && delta < 0.0)
                 {
                    kai = - bb/aa;
                    hopsurf = new_surf;

                    fprintf(f, "Hop from %d to %d (new) fails! Reverse velocity along NCV with kai = %8.5f \n", ct->surface, i, kai);
                 }
              }

              if (hopsurf != ct->surface)
              {
                 fprintf(f,"%8.2f %10d %17d %21f %19f %19f \n",(step*ct->rk_timestep)*1000./PS_TO_AU, ct->surface, hopsurf,
                        ct->surf_prob[i], fabs(ct->tfs_overlap[ct->surface][hopsurf]), (ct->ev_adiab[hopsurf] - ct->ev_adiab[ct->surface])*HARTREE_TO_EV);
                 newline = 1;
              }

              ct->surface = hopsurf;

// velocity adjustment
              if (dftb->do_ncv)
              {
                 if (kai != 0.0)
                 {
                    for (iatom = mdatoms->start; iatom < mdatoms->start+mdatoms->homenr; iatom++)
                    {
                        for (k = 0; k < ct->sites; k++)
                        {
                            for (l = 0; l < ct->site[k].atoms; l++)
                            {
                                if (iatom == ct->site[k].atom[l])
                                {
                                   for (m = 0; m < DIM; m++)
                                       state_global->v[iatom][m] += kai*ct->tfs_ncv[k][l][m]/dftb->phase1[k].mass[l]/AMU_TO_AU/NM_TO_BOHR*PS_TO_AU;
                                }
                            }
                        }
                    }
                 }
                 dftb->do_ncv = 0;
              }
              break;
           }
        }
    }

    if (newline)
       fprintf(f,"\n");

    printf("Current Surface: %d\n", ct->surface);

    ekintot = 0.0;
    for (iatom = mdatoms->start; iatom < mdatoms->start+mdatoms->homenr; iatom++)
    {
        for (k = 0; k < ct->sites; k++)
        {
            for (l = 0; l < ct->site[k].atoms; l++)
            {
                if (iatom == ct->site[k].atom[l])
                {
                   for (m = 0; m < DIM; m++)
                       ekintot += 0.5*dftb->phase1[k].mass[l]*AMU_TO_AU*SQR(state_global->v[iatom][m]*NM_TO_BOHR/PS_TO_AU);
                }
            }
        }
    }
    printf("Kinetic energy is %f\n", ekintot);

// restricted decoherence correction
//    if (adiapop[ct->surface] > 0.001)
//    {
    dot_product = 0.0;
    for (i = 0; i < ct->dim; i++)
    {
        if (i != ct->surface)
        {
            decay_time = (1.0 + 0.1/ekintot)/fabs(ct->ev_adiab[i] - ct->ev_adiab[ct->surface]);

            ct->tfs_popul[i]         *= exp(-ct->rk_timestep/decay_time);
            ct->tfs_popul[i+ct->dim] *= exp(-ct->rk_timestep/decay_time);
            dot_product              += SQR(ct->tfs_popul[i])+SQR(ct->tfs_popul[i+ct->dim]);
        }
    }
    curr_surf_fac = sqrt((1.0 - dot_product)/(SQR(ct->tfs_popul[ct->surface]) + SQR(ct->tfs_popul[ct->dim+ct->surface])));
    ct->tfs_popul[ct->surface]         *= curr_surf_fac;
    ct->tfs_popul[ct->surface+ct->dim] *= curr_surf_fac;

// recalc. diab. from adiab after decoherence corr.
    for (i = 0; i < ct->dim; i++)
    {
        ct->tfs_diab[i] = 0.0; ct->tfs_diab[i+ct->dim] = 0.0;

        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_diab[i]         += ct->tfs_popul[j]*ct->tfs_vector[j][i];
            ct->tfs_diab[i+ct->dim] += ct->tfs_popul[j+ct->dim]*ct->tfs_vector[j][i];
        }
    }
//    }

// charge wave function
    for (i = 0; i < ct->dim; i++)
        ct->wf[i] = ct->tfs_vector[ct->surface][i];

//store hamiltonian as old_ham
    for (i = 0; i < SQR(ct->dim); i++)
        ct->tfl_old_ham[i] = *(ct->hamiltonian[0]+i);

    return 0;
}

int do_global_flux_fssh(int step, charge_transfer_t *ct, dftb_t *dftb, t_state *state_global, t_mdatoms *mdatoms, FILE *f)
{
// Declaration:
//
// hopping probabilities are only determined by populations;
// new: add the threshold for the detection of trivial crossings
//
//                                          -----------Weiwei Xie on July. 2019
// Ref: J. Chem. Theory Comput. 10, 9, 3598 (2014)


    int    i, j, k, l, m, n, tmp, prev_surf, old_surf, hopsurf, new_surf, surftmp, iatom, istep, nstep, *ind, newline, *indtmp, *map;
    double *ham, *old_ham, **ham_full, *new_ham, *adiapop_old, *adiapop, *dadiapop,  *initprob, *diapop;
    double tdse_timestep = 1.e-5*PS_TO_AU, threshold = 0.99;
    double rk_timestep0, probsum, random_number, random_number1, cur_weight, dadiapopsum, overlapmax, probmax;
    double nc, kai1, kai2, kai, delta, aa, bb, cc, dd, ee;
    double decay_time, curr_surf_fac;
    double dot_product, ekintot;


    old_ham  = ct->tfl_old_ham_adiab;
    new_ham  = ct->hamiltonian_adiab;
    ham      = ct->tfl_mean_ham;
    ham_full = ct->tfl_mean_ham_full;

    snew(diapop, ct->dim);
    snew(adiapop_old, ct->dim);
    snew(adiapop, ct->dim);
    snew(dadiapop, ct->dim);
    snew(initprob, ct->dim);
    snew(ind, ct->dim);
    snew(map, ct->dim);
    snew(indtmp, ct->dim);

// set up initial surface
    if (ct->tfs_initialization_step)
    {
       printf("Input wavefunction is:\n");
       for (j = 0; j < ct->dim; j++)
           printf("%f ", ct->wf[j]);
       printf("\n");

       for (i = 0; i < SQR(ct->dim); i++)
       {
           ct->tfl_old_ham[i] = *(ct->hamiltonian[0]+i);
           ham[i] = ct->tfl_old_ham[i];
       }

// diagonalize this Hamiltonian
       dsyev(ct->dim, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);

// record eigenvector
       for (i = 0; i < ct->dim; i++)
           for (j = 0; j < ct->dim; j++)
               ct->tfs_vector[i][j] = ham[i * ct->dim + j];

// determine inital surface using eigenvector as a criterion
       probsum = 0.;
       random_number = drand48();
       for (i = 0; i < ct->dim; i++)
       {
           initprob[i] = 0.;
           for (j = 0; j < ct->dim; j++)
               initprob[i] += ct->tfs_vector[i][j]*ct->wf[j];
           initprob[i] = SQR(initprob[i]);

           probsum += initprob[i];
           if (random_number < probsum)
           {
               ct->surface = i;
               break;
           }
       }

       for (i = 0; i < 2*ct->dim; i++)
       {
           ct->tfs_popul[i] = 0.0; ct->tfs_diab[i]  = 0.0;
       }

       for (i = 0; i < ct->dim; i++)
           ct->tfs_diab[i] = ct->wf[i];

       ct->tfs_popul[ct->surface] = 1.0;

       ct->tfs_initialization_step = 0;

       return 0;
    }

// set up the old and new Hamiltonian for linear interpolation
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            old_ham[i*ct->dim+j] = ct->tfl_old_ham[i*ct->dim+j];
            new_ham[i*ct->dim+j] = ct->hamiltonian[i][j];
        }
    }

// solve time-dependent Schroedinger equation
    nstep = (int)(ct->rk_timestep/tdse_timestep);
    rk_timestep0 = ct->rk_timestep;
    ct->rk_timestep = tdse_timestep;
    printf("Time step for TDSE: %f fs \n", ct->rk_timestep/PS_TO_AU*1000);

    for (istep = 1; istep < nstep+1; istep++)
    {
// linear interpolation for Hamiltonian
        cur_weight = (double)istep/(double)nstep;
        for (i = 0; i < ct->dim; i++)
            for (j = 0; j < ct->dim; j++)
                ham_full[i][j] = (1.0 - cur_weight)*ct->tfl_old_ham[i*ct->dim+j] + cur_weight*ct->hamiltonian[i][j];
// propagate the TDSE in the diabatic representation
        do_rksuite_diab(ct);
    }
    ct->rk_timestep = rk_timestep0;

// diagonalize this Hamiltonian
    for (i = 0; i < SQR(ct->dim); i++)
        ham[i] = new_ham[i];

    dsyev(ct->dim, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);

// store the eigenvectors at time t-dt and time t
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_vector_old[i][j] = ct->tfs_vector[i][j];
            ct->tfs_vector[i][j]     = ham[i*ct->dim + j];
        }
    }

// adiabatic wave function overlaps between time t-dt and t
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_overlap[i][j] = 0.0;

            for (k = 0; k < ct->dim; k++)
                ct->tfs_overlap[i][j] += ct->tfs_vector_old[i][k]*ct->tfs_vector[j][k];
        }
    }

// reassignment of states according to the min-cost algorithm
    assignprob(ct, indtmp);

// state tracking
    for (i = 0; i < ct->dim; i++)
    {
        map[i] = 0; ind[i] = -1;
    }

    for (i = 0; i < ct->dim; i++)
    {
        if (fabs(ct->tfs_overlap[indtmp[i]][i]) >= 0.9)
        {
           ind[i] = indtmp[i];

           for (j = 0; j < ct->dim; j++)
           {
               if (ind[i] == j)
               {
                  map[j] = 1;
                  break;
               }
           }
        }
    }

    for (i = 0; i < ct->dim; i++)
    {
        if (ind[i] == -1)
        {
           for (j = 0; j < ct->dim; j++)
           {
               if (! map[j])
               {
                  ind[i] = j; map[j] = 1;
                  break;
               }
           }
        }
    }

// reassign the active state
    prev_surf = ct->surface;

    if (ct->surface != ind[ct->surface])
    {
       for (i = 0; i < ct->dim; i++)
       {
           if (ct->surface == ind[i])
           {
              printf("Changed from surface %d to surface %d \n", ct->surface, i);
              ct->surface = i;
              break;
           }
       }
    }

    for (i = 0; i < ct->dim; i++)
    {
        if (prev_surf == indtmp[i])
        {
           new_surf = i;
           break;
        }
    }

// check signs
    for (i = 0; i < ct->dim; i++)
    {
        j = ind[i];

        if (ct->tfs_overlap[j][i] < 0.0)
        {
            for (k = 0; k < ct->dim; k++)
            {
                ct->tfs_vector[i][k]  *= -1.0;
                ct->tfs_overlap[k][i] *= -1.0;
            }
        }
    }

// calculate the adiabatic populations at the previous time step
    for (i = 0; i < ct->dim; i++)
        adiapop_old[i] = SQR(ct->tfs_popul[i]) + SQR(ct->tfs_popul[i+ct->dim]);

// calculate expansion coefficients in the adiabatic basis
    for (i = 0; i < ct->dim; i++)
    {
        ct->tfs_popul[i] = 0.0; ct->tfs_popul[i+ct->dim] = 0.0;

        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_popul[i]         += ct->tfs_vector[i][j]*ct->tfs_diab[j];
            ct->tfs_popul[i+ct->dim] += ct->tfs_vector[i][j]*ct->tfs_diab[j+ct->dim];
        }
    }

// record the adiabatic populations at new time
    dadiapopsum = 0.0;
    for (i = 0; i < ct->dim; i++)
    {
        j = ind[i];

        adiapop[i] = SQR(ct->tfs_popul[i])+SQR(ct->tfs_popul[i+ct->dim]);

        dadiapop[i] = adiapop[i] - adiapop_old[j];

        if (dadiapop[i] < 0.0)
           dadiapopsum += dadiapop[i];
    }

// hopping probability
    probmax = 0.0;
    for (i = 0; i < ct->dim; i++)
    {
        ct->surf_prob[i] = 0.0;

        if (dadiapop[ct->surface] < 0.0 && dadiapop[i] > 0.0)
           ct->surf_prob[i] = dadiapop[i]/adiapop_old[ind[ct->surface]]*dadiapop[ct->surface]/dadiapopsum;
    }

    for (i = 0; i < ct->dim; i++)
        printf("%2d  %15.5f  %15.5f \n", i, dadiapop[i]/adiapop_old[ct->surface], ct->surf_prob[i]);

    newline = 0;
    if (ct->surface != new_surf)
    {
       fprintf(f,"%8.2f %10d %17d %21f %19f %19f   UnavoidedCrossing \n",(step*ct->rk_timestep)*1000./PS_TO_AU, ct->surface, new_surf,
              ct->surf_prob[new_surf], fabs(ct->tfs_overlap[ct->surface][new_surf]), (ct->ev_adiab[new_surf] - ct->ev_adiab[ct->surface])*HARTREE_TO_EV);
       newline = 1;
    }

// stochastic algorithm
    probsum = 0.0;
    random_number = drand48();
    for (i = 0; i < ct->dim; i++)
    {
        if (ct->surf_prob[i] >= 0.0 && i != ct->surface)
        {
           probsum += ct->surf_prob[i];

           if (random_number < probsum)
           {
              hopsurf = i;

              dftb->do_ncv = 1;

              if (dftb->do_ncv == 1)
              {
                 ct->ncv_isurf = ct->surface; ct->ncv_jsurf = hopsurf;

                 get_MM_params(ct, dftb);

                 aa = 0.0; bb = 0.0;
                 cc = ct->ev_adiab[hopsurf] - ct->ev_adiab[ct->surface];

                 for (iatom = mdatoms->start; iatom < mdatoms->start+mdatoms->homenr; iatom++)
                 {
                     n = 0;
                     for (k = 0; k < ct->sites; k++)
                     {
                         for (l = 0; l < ct->site[k].atoms; l++)
                         {
                             if (iatom == ct->site[k].atom[l])
                             {
                                for (m = 0; m < DIM; m++)
                                {
                                    ct->tfs_ncv[k][l][m] = (real) (dftb->phase1[k].grad[l][m] + dftb->phase2.grad[n][m])/cc;

                                    aa += SQR(ct->tfs_ncv[k][l][m])/2./dftb->phase1[k].mass[l]/AMU_TO_AU;
                                    bb += ct->tfs_ncv[k][l][m]*state_global->v[iatom][m]*NM_TO_BOHR/PS_TO_AU;

                                }
                             }
                             n++;
                         }
                     }
                 }

// adjust velocity along the direction of nonadiabtic coupling vector
                 delta = SQR(bb) - 4.*aa*cc;
                 if (delta >= 0.0)
                 {
                    kai1 = (- bb + sqrt(delta))/2./aa;
                    kai2 = (- bb - sqrt(delta))/2./aa;

                    kai = (fabs(kai1) < fabs(kai2) ? kai1 : kai2);

                    fprintf(f, "Hop from %d to %d succeeds! Adjust velocity along NCV with kai = %8.5f \n", ct->surface, hopsurf, kai);
                 }
                 else if (delta < 0.0)
                 {
                    kai = - bb/aa;
                    hopsurf = ct->surface;

                    fprintf(f, "Hop from %d to %d fails! Reverse velocity along NCV with kai = %8.5f \n", ct->surface, i, kai);
                 }

// velocity adjustment
                 for (iatom = mdatoms->start; iatom < mdatoms->start+mdatoms->homenr; iatom++)
                 {
                     for (k = 0; k < ct->sites; k++)
                     {
                         for (l = 0; l < ct->site[k].atoms; l++)
                         {
                             if (iatom == ct->site[k].atom[l])
                             {
                                for (m = 0; m < DIM; m++)
                                    state_global->v[iatom][m] += kai*ct->tfs_ncv[k][l][m]/dftb->phase1[k].mass[l]/AMU_TO_AU/NM_TO_BOHR*PS_TO_AU;
                             }
                         }
                     }
                 }
                 dftb->do_ncv = 0;
              }
              else if (dftb->do_ncv == 2)
              {
                 random_number1 = drand48();

                 if (random_number1 > exp(-(ct->ev_adiab[hopsurf] - ct->ev_adiab[ct->surface])/ct->fermi_kt))
                 {
                      hopsurf = ct->surface;
                      fprintf(f,"Hop from %d to %d fails! \n", ct->surface, i);
                 }
                 else
                 {
                    fprintf(f,"Hop from %d to %d succeeds! \n", ct->surface, hopsurf);
                 }
                 dftb->do_ncv = 0;
              }

              if (hopsurf != ct->surface)
              {
                 fprintf(f,"%8.2f %10d %17d %21f %19f %19f \n \n",(step*ct->rk_timestep)*1000./PS_TO_AU, ct->surface, hopsurf,
                        ct->surf_prob[i], fabs(ct->tfs_overlap[ct->surface][hopsurf]), (ct->ev_adiab[hopsurf] - ct->ev_adiab[ct->surface])*HARTREE_TO_EV);
                 newline = 1;
              }

              ct->surface = hopsurf;

              break;
           }
        }
    }

    if (newline)
       fprintf(f,"\n");

    printf("Current Surface: %d\n", ct->surface);

    ekintot = 0.0;
    for (iatom = mdatoms->start; iatom < mdatoms->start+mdatoms->homenr; iatom++)
    {
        for (k = 0; k < ct->sites; k++)
        {
            for (l = 0; l < ct->site[k].atoms; l++)
            {
                if (iatom == ct->site[k].atom[l])
                {
                   for (m = 0; m < DIM; m++)
                       ekintot += 0.5*dftb->phase1[k].mass[l]*AMU_TO_AU*SQR(state_global->v[iatom][m]*NM_TO_BOHR/PS_TO_AU);
                }
            }
        }
    }
    printf("Kinetic energy is %f\n", ekintot);

// restricted decoherence correction
    dot_product = 0.0;
    if (adiapop[ct->surface] > 0.001)
    {
    for (i = 0; i < ct->dim; i++)
    {
        if (i != ct->surface)
        {
            decay_time = (1.0 + 0.1/ekintot)/fabs(ct->ev_adiab[i] - ct->ev_adiab[ct->surface]);

            ct->tfs_popul[i]         *= exp(-ct->rk_timestep/decay_time);
            ct->tfs_popul[i+ct->dim] *= exp(-ct->rk_timestep/decay_time);
            dot_product              += SQR(ct->tfs_popul[i])+SQR(ct->tfs_popul[i+ct->dim]);
        }
    }
    curr_surf_fac = sqrt((1.0 - dot_product)/(SQR(ct->tfs_popul[ct->surface]) + SQR(ct->tfs_popul[ct->dim+ct->surface])));
    ct->tfs_popul[ct->surface]         *= curr_surf_fac;
    ct->tfs_popul[ct->surface+ct->dim] *= curr_surf_fac;

// recalc. diab. from adiab after decoherence corr.
    for (i = 0; i < ct->dim; i++)
    {
        ct->tfs_diab[i] = 0.0; ct->tfs_diab[i+ct->dim] = 0.0;

        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_diab[i]         += ct->tfs_popul[j]*ct->tfs_vector[j][i];
            ct->tfs_diab[i+ct->dim] += ct->tfs_popul[j+ct->dim]*ct->tfs_vector[j][i];
        }
    }
    }

// charge wave function
    for (i = 0; i < ct->dim; i++)
        ct->wf[i] = ct->tfs_vector[ct->surface][i];

//store hamiltonian as old_ham
    for (i = 0; i < SQR(ct->dim); i++)
        ct->tfl_old_ham[i] = *(ct->hamiltonian[0]+i);

    return 0;
}


/*int do_global_flux_fssh(int step, charge_transfer_t *ct, dftb_t *dftb, t_state *state_global, t_mdatoms *mdatoms, FILE *f)
{
// Declaration:
//
// hopping probabilities are only determined by populations;
// new: add the threshold for the detection of trivial crossings
//
//                                          -----------Weiwei Xie on July. 2019
// Ref: J. Chem. Theory Comput. 10, 9, 3598 (2014)


    int    i, j, k, l, m, n, tmp, prev_surf, old_surf, hopsurf, new_surf, surftmp, iatom, istep, nstep, *ind;
    double *ham, *old_ham, **ham_full, *new_ham, *adiapop_old, *adiapop, *dadiapop,  *initprob, *diapop;
    double tdse_timestep = 1.e-6*PS_TO_AU, threshold = 0.99;
    double rk_timestep0, probsum, random_number, random_number1, cur_weight, dadiapopsum, overlapmax, probmax;
    double nc, kai1, kai2, kai, delta, aa, bb, cc, dd, ee;
    double decay_time, curr_surf_fac;
    double dot_product, ekintot;


    old_ham  = ct->tfl_old_ham_adiab;
    new_ham  = ct->hamiltonian_adiab;
    ham      = ct->tfl_mean_ham;
    ham_full = ct->tfl_mean_ham_full;

    snew(diapop, ct->dim);
    snew(adiapop_old, ct->dim);
    snew(adiapop, ct->dim);
    snew(dadiapop, ct->dim);
    snew(initprob, ct->dim);
    snew(ind, ct->dim);

// set up initial surface
    if (ct->tfs_initialization_step)
    {
       printf("Input wavefunction is:\n");
       for (j = 0; j < ct->dim; j++)
           printf("%f ", ct->wf[j]);
       printf("\n");

       for (i = 0; i < SQR(ct->dim); i++)
       {
           ct->tfl_old_ham[i] = *(ct->hamiltonian[0]+i);
           ham[i] = ct->tfl_old_ham[i];
       }

// diagonalize this Hamiltonian
       dsyev(ct->dim, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);

// record eigenvector
       for (i = 0; i < ct->dim; i++)
           for (j = 0; j < ct->dim; j++)
               ct->tfs_vector[i][j] = ham[i * ct->dim + j];

// determine inital surface using eigenvector as a criterion
       probsum = 0.;
       random_number = drand48();
       for (i = 0; i < ct->dim; i++)
       {
           initprob[i] = 0.;
           for (j = 0; j < ct->dim; j++)
               initprob[i] += ct->tfs_vector[i][j]*ct->wf[j];
           initprob[i] = SQR(initprob[i]);

           probsum += initprob[i];
           if (random_number < probsum)
           {
               ct->surface = i;
               break;
           }
       }

       for (i = 0; i < 2*ct->dim; i++)
       {
           ct->tfs_popul[i] = 0.0; ct->tfs_diab[i]  = 0.0;
       }

       for (i = 0; i < ct->dim; i++)
           ct->tfs_diab[i] = ct->wf[i];

       ct->tfs_popul[ct->surface] = 1.0;

       ct->tfs_initialization_step = 0;

       return 0;
    }

// set up the old and new Hamiltonian for linear interpolation
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            old_ham[i*ct->dim+j] = ct->tfl_old_ham[i*ct->dim+j];
            new_ham[i*ct->dim+j] = ct->hamiltonian[i][j];
        }
    }


    ekintot = 0.0;
    for (iatom = mdatoms->start; iatom < mdatoms->start+mdatoms->homenr; iatom++)
    {
        for (k = 0; k < ct->sites; k++)
        {
            for (l = 0; l < ct->site[k].atoms; l++)
            {
                if (iatom == ct->site[k].atom[l])
                {
                   for (m = 0; m < DIM; m++)
                       ekintot += 0.5*dftb->phase1[k].mass[l]*AMU_TO_AU*SQR(state_global->v[iatom][m]*NM_TO_BOHR/PS_TO_AU);
                }
            }
        }
    }
    printf("Kinetic energy is %f\n", ekintot);

// solve time-dependent Schroedinger equation
    nstep = (int)(ct->rk_timestep/tdse_timestep);
    rk_timestep0 = ct->rk_timestep;
    ct->rk_timestep = tdse_timestep;
    printf("Time step for TDSE: %f fs \n", ct->rk_timestep/PS_TO_AU*1000);

    for (istep = 1; istep < nstep+1; istep++)
    {
// linear interpolation for Hamiltonian
        cur_weight = (double)istep/(double)nstep;
        for (i = 0; i < ct->dim; i++)
            for (j = 0; j < ct->dim; j++)
                ham_full[i][j] = (1.0 - cur_weight)*ct->tfl_old_ham[i*ct->dim+j] + cur_weight*ct->hamiltonian[i][j];
// propagate the TDSE in the diabatic representation
        do_rksuite_diab(ct);

// diagonalize this Hamiltonian
    for (i = 0; i < SQR(ct->dim); i++)
        ham[i] = (1.0-cur_weight)*old_ham[i]+cur_weight*new_ham[i];

    dsyev(ct->dim, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);

// store the eigenvectors at time t-dt and time t
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_vector_old[i][j] = ct->tfs_vector[i][j];
            ct->tfs_vector[i][j]     = ham[i*ct->dim + j];
        }
    }

// adiabatic wave function overlaps between time t-dt and t
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_overlap[i][j] = 0.0;

            for (k = 0; k < ct->dim; k++)
                ct->tfs_overlap[i][j] += ct->tfs_vector_old[i][k]*ct->tfs_vector[j][k];
        }
    }

// check signs
    for (i = 0; i < ct->dim; i++)
    {
        if (ct->tfs_overlap[i][i] < 0.0)
        {
            for (k = 0; k < ct->dim; k++)
            {
                ct->tfs_vector[i][k]  *= -1.0;
                ct->tfs_overlap[k][i] *= -1.0;
            }
        }
    }

// calculate the adiabatic populations at the previous time step
    for (i = 0; i < ct->dim; i++)
        adiapop_old[i] = SQR(ct->tfs_popul[i]) + SQR(ct->tfs_popul[i+ct->dim]);

// calculate expansion coefficients in the adiabatic basis
    for (i = 0; i < ct->dim; i++)
    {
        ct->tfs_popul[i] = 0.0; ct->tfs_popul[i+ct->dim] = 0.0;

        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_popul[i]         += ct->tfs_vector[i][j]*ct->tfs_diab[j];
            ct->tfs_popul[i+ct->dim] += ct->tfs_vector[i][j]*ct->tfs_diab[j+ct->dim];
        }
    }

// record the adiabatic populations at new time
    dadiapopsum = 0.0;
    for (i = 0; i < ct->dim; i++)
    {
        adiapop[i] = SQR(ct->tfs_popul[i])+SQR(ct->tfs_popul[i+ct->dim]);

        dadiapop[i] = adiapop[i] - adiapop_old[i];

        if (dadiapop[i] < 0.0)
           dadiapopsum += dadiapop[i];
    }

// hopping probability
    probmax = 0.0;
    for (i = 0; i < ct->dim; i++)
    {
        ct->surf_prob[i] = 0.0;

        if (dadiapop[ct->surface] < 0.0 && dadiapop[i] > 0.0)
           ct->surf_prob[i] = dadiapop[i]/adiapop_old[ct->surface]*dadiapop[ct->surface]/dadiapopsum;
    }

//    for (i = 0; i < ct->dim; i++)
 //       printf("%2d  %15.5f  %15.5f \n", i, dadiapop[i]/adiapop_old[ct->surface], ct->surf_prob[i]);

    prev_surf = ct->surface;
    for (i = 0; i < ct->dim; i++)
    {
        if (ct->surf_prob[i] > threshold)
           ct->surface = i;
    }

// stochastic algorithm
    probsum = 0.0;
    random_number = drand48();
    if (prev_surf == ct->surface)
    {
    for (i = 0; i < ct->dim; i++)
    {
        if (ct->surf_prob[i] >= 0.0 && i != ct->surface)
        {
           probsum += ct->surf_prob[i];

           if (random_number < probsum)
           {
              hopsurf = i;

              dftb->do_ncv = 2;

              if (dftb->do_ncv == 2)
              {
                 random_number1 = drand48();

                 if (random_number1 > exp(-(ct->ev_adiab[hopsurf] - ct->ev_adiab[ct->surface])/ct->fermi_kt))
                    hopsurf = ct->surface;

                 dftb->do_ncv = 0;
              }

              if (hopsurf != ct->surface)
                 fprintf(f,"%8.2f %10d %17d %21f %19f \n",(istep*ct->rk_timestep+step*rk_timestep0)*1000./PS_TO_AU, ct->surface, hopsurf, ct->surf_prob[i], random_number);

              ct->surface = hopsurf;

              break;
           }
        }
    }
    }

// restricted decoherence correction
    dot_product = 0.0;
    for (i = 0; i < ct->dim; i++)
    {
        if (i != ct->surface)
        {
            decay_time = (1.0 + 0.1/ekintot)/fabs(ct->ev_adiab[i] - ct->ev_adiab[ct->surface]);

            ct->tfs_popul[i]         *= exp(-ct->rk_timestep/decay_time);
            ct->tfs_popul[i+ct->dim] *= exp(-ct->rk_timestep/decay_time);
            dot_product              += SQR(ct->tfs_popul[i])+SQR(ct->tfs_popul[i+ct->dim]);
        }
    }
    curr_surf_fac = sqrt((1.0 - dot_product)/(SQR(ct->tfs_popul[ct->surface]) + SQR(ct->tfs_popul[ct->dim+ct->surface])));
    ct->tfs_popul[ct->surface]         *= curr_surf_fac;
    ct->tfs_popul[ct->surface+ct->dim] *= curr_surf_fac;

// recalc. diab. from adiab after decoherence corr.
    for (i = 0; i < ct->dim; i++)
    {
        ct->tfs_diab[i] = 0.0; ct->tfs_diab[i+ct->dim] = 0.0;

        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_diab[i]         += ct->tfs_popul[j]*ct->tfs_vector[j][i];
            ct->tfs_diab[i+ct->dim] += ct->tfs_popul[j+ct->dim]*ct->tfs_vector[j][i];
        }
    }
    }
    ct->rk_timestep = rk_timestep0;

    printf("Current Surface: %d\n", ct->surface);

// charge wave function
    for (i = 0; i < ct->dim; i++)
        ct->wf[i] = ct->tfs_vector[ct->surface][i];

//store hamiltonian as old_ham
    for (i = 0; i < SQR(ct->dim); i++)
        ct->tfl_old_ham[i] = *(ct->hamiltonian[0]+i);

    return 0;
}
*/
int do_decoherence_induced_fssh(int step, charge_transfer_t *ct, dftb_t *dftb, t_state *state_global, t_mdatoms *mdatoms, FILE *f)
{
// Declaration:
//
//                                          -----------Weiwei Xie on July. 2019
// Ref: J.Chem.Phys. 137, 22A545 (2012)


    int    i, j, k, l, m, n, tmp, prev_surf, old_surf, hopsurf, new_surf, surftmp, iatom, istep, nstep, *ind, indmax, *indtmp, *indtmp1, *map;
    double *ham, *old_ham, **ham_full, *new_ham, *adiapop_old, *adiapop, *initprob, *probtmp1, *probtmp2, *diapop, *diatmp;
    double tdse_timestep = 1.e-5*PS_TO_AU;
    double rk_timestep0, probsum, random_number, cur_weight;
    double popul_norm, overlapmax, overlaptmp1, overlaptmp2, adiapoptmp0, probtmp0, probtmp, diapoptmp;
    double nc, kai1, kai2, kai, delta, aa, bb, cc, dd, ee;
    double decay_time, curr_surf_fac, trivial_cross_type;
    double dot_product, ekintot;
    double delta_min;


    old_ham  = ct->tfl_old_ham_adiab;
    new_ham  = ct->hamiltonian_adiab;
    ham      = ct->tfl_mean_ham;
    ham_full = ct->tfl_mean_ham_full;

    snew(diapop, ct->dim);
    snew(diatmp, 2*ct->dim);
    snew(adiapop_old, ct->dim);
    snew(adiapop, ct->dim);
    snew(initprob, ct->dim);
    snew(probtmp1, ct->dim);
    snew(probtmp2, ct->dim);

    snew(ind, ct->dim);
    snew(indtmp, ct->dim);
    snew(indtmp1, ct->dim);
    snew(map, ct->dim);

// set up initial surface
    if (ct->tfs_initialization_step)
    {
       printf("Input wavefunction is:\n");
       for (j = 0; j < ct->dim; j++)
           printf("%f ", ct->wf[j]);
       printf("\n");

       for (i = 0; i < SQR(ct->dim); i++)
       {
           ct->tfl_old_ham[i] = *(ct->hamiltonian[0]+i);
           ham[i] = ct->tfl_old_ham[i];
       }

// diagonalize this Hamiltonian
       dsyev(ct->dim, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);

// record eigenvector
       for (i = 0; i < ct->dim; i++)
           for (j = 0; j < ct->dim; j++)
               ct->tfs_vector[i][j] = ham[i * ct->dim + j];

// determine inital surface using eigenvector as a criterion
       probsum = 0.;
       random_number = drand48();
       for (i = 0; i < ct->dim; i++)
       {
           initprob[i] = 0.;
           for (j = 0; j < ct->dim; j++)
               initprob[i] += ct->tfs_vector[i][j]*ct->wf[j];
           initprob[i] = SQR(initprob[i]);

           probsum += initprob[i];
           if (random_number < probsum)
           {
               ct->surface = i;
               break;
           }
       }

       for (i = 0; i < 2*ct->dim; i++)
       {
           ct->tfs_popul[i] = 0.0; ct->tfs_diab[i]  = 0.0;
       }

       for (i = 0; i < ct->dim; i++)
           ct->tfs_diab[i] = ct->wf[i];

       ct->tfs_popul[ct->surface] = 1.0;

       ct->tfs_initialization_step = 0;

       return 0;
    }

// set up the old and new Hamiltonian for linear interpolation
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            old_ham[i*ct->dim+j] = ct->tfl_old_ham[i*ct->dim+j];
            new_ham[i*ct->dim+j] = ct->hamiltonian[i][j];
        }
    }

// store the vector at preceding time step
/*    for (i = 0; i < ct->dim; i++)
        for (j = 0; j < ct->dim; j++)
            ct->tfs_vector_old2[i][j] = ct->tfs_vector[i][j];
*/
// solve time-dependent Schroedinger equation
    nstep = (int)(ct->rk_timestep/tdse_timestep);
    rk_timestep0 = ct->rk_timestep;
    ct->rk_timestep = tdse_timestep;
    printf("Time step for TDSE: %f fs \n", ct->rk_timestep/PS_TO_AU*1000);

    for (istep = 1; istep < nstep+1; istep++)
    {
// linear interpolation for Hamiltonian
        cur_weight = (double)istep/(double)nstep;
        for (i = 0; i < ct->dim; i++)
            for (j = 0; j < ct->dim; j++)
                ham_full[i][j] = (1.0 - cur_weight)*ct->tfl_old_ham[i*ct->dim+j] + cur_weight*ct->hamiltonian[i][j];
// propagate the TDSE in the diabatic representation
        do_rksuite_diab(ct);
    }
    ct->rk_timestep = rk_timestep0;

// diagonalize this Hamiltonian
    for (i = 0; i < SQR(ct->dim); i++)
        ham[i] = new_ham[i];

    dsyev(ct->dim, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);

// store the eigenvectors at time t-dt and time t
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_vector_old[i][j] = ct->tfs_vector[i][j];
            ct->tfs_vector[i][j]     = ham[i*ct->dim + j];
        }
    }

// adiabatic wave function overlaps between time t-dt and t
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_overlap[i][j] = 0.0;

            for (k = 0; k < ct->dim; k++)
                ct->tfs_overlap[i][j] += ct->tfs_vector_old[i][k]*ct->tfs_vector[j][k];
        }
    }

// reassignment of states according to the maximum wave function overlaps
/*    for (i = 0; i < ct->dim; i++)
    {
        overlapmax = 0.0;
        for (j = 0; j < ct->dim; j++)
        {
            if (overlapmax < fabs(ct->tfs_overlap[j][i]))
            {
               overlapmax = fabs(ct->tfs_overlap[j][i]);
               indtmp[i] = j;
            }
        }
    }
*/
// reassignment of states according to the min-cost algorithm
    assignprob(ct, ind);

// calculate the adiabatic populations at the previous time step
    for (i = 0; i < ct->dim; i++)
        adiapop_old[i] = SQR(ct->tfs_popul[i]) + SQR(ct->tfs_popul[i+ct->dim]);

// hopping probabilities based on Tully's fewest-switches anzatz
    for (i = 0; i < ct->dim; i++)
    {
//     nonadibatic couplings based on Hammes-Schiffer and Tully, JCP 101, 4657 (1994)
//        nc = (ct->tfs_overlap[prev_surf][i] - ct->tfs_overlap[j][ct->surface])/2.0/ct->rk_timestep;
//     ndibatic couplings based norm-preserving interpolation by Benjamin, JPCL 5, 2351 (2014)
        overlaptmp1 = sqrt(fabs(1.0 - SQR(ct->tfs_overlap[i][i]) - SQR(ct->tfs_overlap[ct->surface][i])));
        overlaptmp2 = (- ct->tfs_overlap[i][ct->surface]*ct->tfs_overlap[i][i] - ct->tfs_overlap[ct->surface][ct->surface]*ct->tfs_overlap[ct->surface][i])/overlaptmp1;

        aa = - sin(acos(ct->tfs_overlap[i][i]) - asin(ct->tfs_overlap[i][ct->surface]))/(acos(ct->tfs_overlap[i][i]) - asin(ct->tfs_overlap[i][ct->surface]));
        bb =   sin(acos(ct->tfs_overlap[i][i]) + asin(ct->tfs_overlap[i][ct->surface]))/(acos(ct->tfs_overlap[i][i]) + asin(ct->tfs_overlap[i][ct->surface]));
        cc =   sin(acos(ct->tfs_overlap[ct->surface][ct->surface]) - asin(ct->tfs_overlap[ct->surface][i]))/(acos(ct->tfs_overlap[ct->surface][ct->surface]) - asin(ct->tfs_overlap[ct->surface][i]));
        dd =   sin(acos(ct->tfs_overlap[ct->surface][ct->surface]) + asin(ct->tfs_overlap[ct->surface][i]))/(acos(ct->tfs_overlap[ct->surface][ct->surface]) + asin(ct->tfs_overlap[ct->surface][i]));
        ee = 2.0*asin(overlaptmp1*(overlaptmp1*overlaptmp2*asin(overlaptmp1) + (sqrt((1.0 - SQR(overlaptmp1))*(1.0 - SQR(overlaptmp2))) - 1.0)*asin(overlaptmp2)))/(SQR(asin(overlaptmp1)) - SQR(asin(overlaptmp2)));

        nc = (acos(ct->tfs_overlap[i][i])*(aa + bb) + asin(ct->tfs_overlap[ct->surface][i])*(cc + dd) + ee)/2.0/ct->rk_timestep;

        probtmp1[i] = (2.0*(ct->tfs_popul[ct->surface]*ct->tfs_popul[i]+ct->tfs_popul[ct->dim+ct->surface]*ct->tfs_popul[ct->dim+i])*nc)*ct->rk_timestep/adiapop_old[ct->surface];

        if (probtmp1[i] < 0.0 || i == ct->surface)
           probtmp1[i] = 0.0;
    }

// calculate expansion coefficients in the adiabatic basis
    for (i = 0; i < ct->dim; i++)
    {
        ct->tfs_popul[i] = 0.0; ct->tfs_popul[i+ct->dim] = 0.0;

        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_popul[i]         += ct->tfs_vector[i][j]*ct->tfs_diab[j];
            ct->tfs_popul[i+ct->dim] += ct->tfs_vector[i][j]*ct->tfs_diab[j+ct->dim];
        }
    }

// record the adiabatic populations at new time
    for (i = 0; i < ct->dim; i++)
        adiapop[i] = SQR(ct->tfs_popul[i])+SQR(ct->tfs_popul[i+ct->dim]);

// SC-FSSH, avoid the numerical instability of the NCs
    prev_surf = ct->surface;
    for (i = 0; i < ct->dim; i++)
    {
        if (prev_surf == ind[i])
        {
           new_surf = i;
           break;
        }
    }


    if (ct->surface != new_surf)
    {
       probtmp1[new_surf] = (adiapop_old[ct->surface] - adiapop[ct->surface])/adiapop_old[ct->surface];

       for (i = 0; i < ct->dim; i++)
       {
           if (i != new_surf)
           {
              probtmp1[new_surf] -= probtmp1[i];
           }
       }

       printf("SC is actived: %d  %d  %10f  %10f\n", ct->surface, new_surf, (adiapop_old[ct->surface] - adiapop[ct->surface])/adiapop_old[ct->surface], probtmp1[new_surf]);
    }

    for (i = 0; i < ct->dim; i++)
        printf("%2d  %2d  %15.5f \n", i, ind[i], probtmp1[i]);

// stochastic algorithm
    probsum = 0.0;
    random_number = drand48();
    for (i = 0; i < ct->dim; i++)
    {
        ct->surf_prob[i] = probtmp1[i];

        if (ct->surf_prob[i] >= 0.0 && i != ct->surface)
        {
           probsum += ct->surf_prob[i];

           if (random_number < probsum)
           {
              hopsurf = i;
// calculate the nonadibatic coupling vector for velocity adjustment
              dftb->do_ncv = 1;

              if (dftb->do_ncv)
              {
                 ct->ncv_isurf = ct->surface; ct->ncv_jsurf = hopsurf;

                 get_MM_params(ct, dftb);

                 aa = 0.0; bb = 0.0;
                 cc = ct->ev_adiab[hopsurf] - ct->ev_adiab[ct->surface];

                 for (iatom = mdatoms->start; iatom < mdatoms->start+mdatoms->homenr; iatom++)
                 {
                     n = 0;
                     for (k = 0; k < ct->sites; k++)
                     {
                         for (l = 0; l < ct->site[k].atoms; l++)
                         {
                             if (iatom == ct->site[k].atom[l])
                             {
                                for (m = 0; m < DIM; m++)
                                {
                                    ct->tfs_ncv[k][l][m] = (real) (dftb->phase1[k].grad[l][m] + dftb->phase2.grad[n][m])/cc;

                                    aa += SQR(ct->tfs_ncv[k][l][m])/2./dftb->phase1[k].mass[l]/AMU_TO_AU;
                                    bb += ct->tfs_ncv[k][l][m]*state_global->v[iatom][m]*NM_TO_BOHR/PS_TO_AU;

                                }
                             }
                             n++;
                         }
                     }
                 }

// adjust velocity along the direction of nonadiabtic coupling vector
                 delta = SQR(bb) - 4.*aa*cc;
                 if (delta >= 0.0)
                 {
                    kai1 = (- bb + sqrt(delta))/2./aa;
                    kai2 = (- bb - sqrt(delta))/2./aa;

                    kai = (fabs(kai1) < fabs(kai2) ? kai1 : kai2);

                    printf("Hop from %d to %d succeeds! Adjust velocity along NCV with kai = %8.5f \n", ct->surface, hopsurf, kai);
                 }
                 else
                 {
                    kai = - bb/aa;
                    hopsurf = ct->surface;

                    printf("Hop from %d to %d fails! Reverse velocity along NCV with kai = %8.5f \n", ct->surface, i, kai);
                 }
              }

              if (hopsurf != ct->surface)
                 fprintf(f,"%8.2f %10d %17d %21f %19f \n",(step*ct->rk_timestep)*1000./PS_TO_AU, ct->surface, hopsurf, ct->surf_prob[i], random_number);

              ct->surface = hopsurf;

// velocity adjustment
              if (dftb->do_ncv)
              {
                 if (kai != 0.0)
                 {
                    for (iatom = mdatoms->start; iatom < mdatoms->start+mdatoms->homenr; iatom++)
                    {
                        for (k = 0; k < ct->sites; k++)
                        {
                            for (l = 0; l < ct->site[k].atoms; l++)
                            {
                                if (iatom == ct->site[k].atom[l])
                                {
                                   for (m = 0; m < DIM; m++)
                                       state_global->v[iatom][m] += kai*ct->tfs_ncv[k][l][m]/dftb->phase1[k].mass[l]/AMU_TO_AU/NM_TO_BOHR*PS_TO_AU;
                                }
                            }
                        }
                    }
                 }
                 dftb->do_ncv = 0;
              }
              break;
           }
        }
    }

    printf("Current Surface: %d\n", ct->surface);

    ekintot = 0.0;
    for (iatom = mdatoms->start; iatom < mdatoms->start+mdatoms->homenr; iatom++)
    {
        for (k = 0; k < ct->sites; k++)
        {
            for (l = 0; l < ct->site[k].atoms; l++)
            {
                if (iatom == ct->site[k].atom[l])
                {
                   for (m = 0; m < DIM; m++)
                       ekintot += 0.5*dftb->phase1[k].mass[l]*AMU_TO_AU*SQR(state_global->v[iatom][m]*NM_TO_BOHR/PS_TO_AU);
                }
            }
        }
    }
    printf("Kinetic energy is %f\n", ekintot);

// restricted decoherence correction
    dot_product = 0.0;
    if (adiapop[ct->surface] > 0.001)
    {
       for (i = 0; i < ct->dim; i++)
       {
           if (i != ct->surface)
           {
               decay_time = (1.0 + 0.1/ekintot)/fabs(ct->ev_adiab[i] - ct->ev_adiab[ct->surface]);

               ct->tfs_popul[i]         *= exp(-ct->rk_timestep/decay_time);
               ct->tfs_popul[i+ct->dim] *= exp(-ct->rk_timestep/decay_time);
               dot_product              += SQR(ct->tfs_popul[i])+SQR(ct->tfs_popul[i+ct->dim]);
           }
       }
       curr_surf_fac = sqrt((1.0 - dot_product)/(SQR(ct->tfs_popul[ct->surface]) + SQR(ct->tfs_popul[ct->dim+ct->surface])));
       ct->tfs_popul[ct->surface]         *= curr_surf_fac;
       ct->tfs_popul[ct->surface+ct->dim] *= curr_surf_fac;

// recalc. diab. from adiab after decoherence corr.
       for (i = 0; i < ct->dim; i++)
       {
           ct->tfs_diab[i] = 0.0; ct->tfs_diab[i+ct->dim] = 0.0;

           for (j = 0; j < ct->dim; j++)
           {
               ct->tfs_diab[i]         += ct->tfs_popul[j]*ct->tfs_vector[j][i];
               ct->tfs_diab[i+ct->dim] += ct->tfs_popul[j+ct->dim]*ct->tfs_vector[j][i];
           }
       }
    }

// charge wave function
    for (i = 0; i < ct->dim; i++)
        ct->wf[i] = ct->tfs_vector[ct->surface][i];

//store hamiltonian as old_ham
    for (i = 0; i < SQR(ct->dim); i++)
        ct->tfl_old_ham[i] = *(ct->hamiltonian[0]+i);

    return 0;
}


// new sf implemention needs to be added / done May2018

int do_tully_local(int step, charge_transfer_t *ct, dftb_t *dftb, t_state *state_global, t_mdatoms *mdatoms, FILE *f)
{

    //combination of flexible surface hopping (JPCL 2013, 4, 1888-1895) and modified hopping prob. calculation (JPCL 2014, 5, 713-719)  dholub difference found

    int    i, j, k, l;
    double old_energy, energy, energy1, energy2, *ham, *old_ham, *new_ham, **ham_full, fermi_energy, fermi_upper, fermi_lower, *initprob, probsum;
    double random_number, nc;                      // in the interval (0,1)
    double cumulated_probability;              // sum of ct->surf_prob[0..k]
    double popul_norm, prob_b;                         // sum of SQR(ct->tfs_popul) except SQR(ct->tfs_popul[ct->surface])
    double decay_time, current_surface_factor; // in the decoherence algorithm
    double dot_product;
    int    dt_step = 0;

    old_ham  = ct->tfl_old_ham_adiab;
    new_ham  = ct->hamiltonian_adiab;
    ham      = ct->tfl_mean_ham;
    ham_full = ct->tfl_mean_ham_full;

    snew(initprob, ct->dim);
// calculate the charges from the wave function  [JUMP3 turn on?]
//for (i=0; i<ct->dim; i++)
// ct->q_act[i] = ct->q_old[i] = SQR(ct->wf[i]);

    // copy as the first vector
    energy = 1.e10;


    ////////assemble hamiltonian

    k = 0;
    for (i = 0; i < ct->dim; i++)
    {
        if (ct->tfl_is_in_system[i])
        {
            l = 0;
            for (j = 0; j < ct->dim; j++)
            {
                if (ct->tfl_is_in_system[j])
                {
                    new_ham[k+l*ct->tfl_num_of_states] = ct->hamiltonian[i][j];
                    new_ham[l+k*ct->tfl_num_of_states] = ct->hamiltonian[j][i];
                    l++;
                }
            }
            k++;
        }
    }



    if (ct->tfs_initialization_step)
    {
        for (i = 0; i < SQR(ct->dim); i++)
        {
            ct->tfl_old_ham[i] = *(ct->hamiltonian[0]+i);
        }
    }

    k = 0;
    for (i = 0; i < ct->dim; i++)
    {
        if (ct->tfl_is_in_system[i])
        {
            l = 0;
            for (j = 0; j < ct->dim; j++)
            {
                if (ct->tfl_is_in_system[j])
                {
                    old_ham[k+l*ct->tfl_num_of_states] = ct->tfl_old_ham[i*ct->dim+j];
                    old_ham[l+k*ct->tfl_num_of_states] = ct->tfl_old_ham[j*ct->dim+i];
                    l++;
                }
            }
            k++;
        }
    }

    ///print out system
    printf("Current system:\n");
    for (i = 0; i < ct->dim; i++)
    {
        printf("%d ", ct->tfl_is_in_system[i]);
    }
    printf("\n");


    // printf("Hamiltonian:\n"); //debug
    //for(i=0;i<ct->dim;i++) for(j=i;j<ct->dim;j++) printf("%8.5f \n",ct->hamiltonian[i][j]); //debug
    ct->rk_timestep /= 100.0; //dholub stepsize for electronic propagation; tried  50tepps  for CT -> crash  Jun2018
    for (dt_step = 0; dt_step < 101; dt_step++)
    {

        double cur_weight = (double)dt_step/100.0;
        //build hamiltonian
        for (i = 0; i < SQR(ct->tfl_num_of_states); i++)
        {
            //i -> dtstep...
            ham[i] = (1.0-cur_weight)*old_ham[i]+cur_weight*new_ham[i];
        }

        for (i = 0; i < ct->dim; i++)
        {
            for (j = 0; j < ct->dim; j++)
            {
                ham_full[i][j] = (1.0-cur_weight)*ct->tfl_old_ham[i*ct->dim+j]+cur_weight*ct->hamiltonian[i][j];

            }
        }

        // diagonalize this Hamiltonian
        dsyev(ct->tfl_num_of_states, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);

        if (ct->tfs_initialization_step)
        {
            /* let us start with tfs_vector_old == tfs_vector */

            for (i = 0; i < ct->tfl_num_of_states_old; i++)
            {
                k = 0;
                for (j = 0; j < ct->dim; j++)
                {

                    if (ct->tfl_is_in_system[j])
                    {
                        ct->tfs_vector_old[i][j] = ham[i * ct->tfl_num_of_states_old + k]; k++;
                        ct->tfs_vector[i][j] = ct->tfs_vector_old[i][j];
                    }
                    else
                    {
                        ct->tfs_vector_old[i][j] = 0.;
                        ct->tfs_vector[i][j]     = 0.;
                    }
                }
            }
            for (; i < ct->dim; i++)
            {
                for (j = 0; j < ct->dim; j++)
                {
                    ct->tfs_vector_old[i][j] = 0.;
                    ct->tfs_vector[i][j]     = 0.;
                }
            }
// determine inital surface using eigenvector as a criterion
            probsum = 0.;
            random_number = drand48();
            for (i = 0; i < ct->dim; i++)
            {
                initprob[i] = 0.;
                for (j = 0; j < ct->dim; j++)
                {
                    initprob[i] += ct->tfs_vector_old[i][j]*ct->wf[j];
                }
                initprob[i] = SQR(initprob[i]);

                probsum += initprob[i];
                if (random_number < probsum)
                {
                    ct->surface = i;
                    break;
                }
            }
            for (i = 0; i < ct->dim; i++)
            {
                ct->tfs_diab[i] = ct->tfs_vector_old[ct->surface][i];
                ct->wf[i]       = ct->tfs_vector_old[ct->surface][i];
            }
            for (i = 0; i < 2*ct->dim; i++)
            {
                ct->tfs_popul[i] = 0.0;
            }
            ct->tfs_popul[ct->surface] = 1.0;

            ct->tfs_initialization_step = 0;

            ct->rk_timestep *= 100.0;

            return 0;
        }
        else
        {
            /* push tfs_vector to tfs_vector_old */
            for (i = 0; i < ct->dim; i++)
            {
                for (j = 0; j < ct->dim; j++)
                {
                    ct->tfs_vector_old[i][j] = ct->tfs_vector[i][j];
                }
            }
        }
        ////GOT HERE
        /* enter new eigenvectors to tfs_vector */
        for (i = 0; i < ct->tfl_num_of_states; i++)
        {
            k = 0;
            for (j = 0; j < ct->dim; j++)
            {
                if (ct->tfl_is_in_system[j])
                {
                    ct->tfs_vector[i][j] = ham[i * ct->tfl_num_of_states + k]; k++;
                }
                else
                {
                    ct->tfs_vector[i][j] = 0;
                }

            }
        }
        for (; i < ct->dim; i++)
        {
            for (j = 0; j < ct->dim; j++)
            {
                ct->tfs_vector[i][j] = 0;
            }
        }



        /////Attention!!: number of old and new sites may differ. Account for this!
        /* calculate tfs_overlap[i][j] = <tfs_vector_old[i] | tfs_vector[j]> */
        for (i = 0; i < ct->dim; i++)
        {
            for (j = 0; j < ct->dim; j++)
            {
                ct->tfs_overlap[i][j] = 0.;

                for (k = 0; k < ct->dim; k++)
                {
                    ct->tfs_overlap[i][j] += ct->tfs_vector_old[i][k] * ct->tfs_vector[j][k];

                }
            }
        }

        /////GOT HERE
        //////////////////////////////////////check signs
        for (i = 0; i < ct->dim; i++)
        {
            int sim_ind = 0;
            dot_product = 0.0;

            for (k = 0; k < ct->dim; k++)
            {
                if (fabs(ct->tfs_overlap[k][i]) > dot_product)
                {
                    sim_ind = k; dot_product = fabs(ct->tfs_overlap[k][i]);
                }
            }

            if (ct->tfs_overlap[sim_ind][i] < 0.00)
            {
                for (k = 0; k < ct->dim; k++)
                {
                    ct->tfs_vector[i][k] *= -1.0; ct->tfs_overlap[k][i] *= -1.0;
                }
            }
        }

        /////GOT HERE

        ////////store old adiab. population of surface
        double ca_re = ct->tfs_popul[ct->surface];
        double ca_im = ct->tfs_popul[ct->surface+ct->dim];

        //////       Statingpoint after inisialization:   propagate in diabatic representation
        do_rksuite_diab(ct);

        // should work like this, but think again:
        // calc. new popul (remember to store some values from the old one)
        dot_product = 0.0;
        for (i = 0; i < ct->tfl_num_of_states; i++)
        {
            ct->tfs_popul[i] = 0.0; ct->tfs_popul[i+ct->dim] = 0.0;
            for (j = 0; j < ct->dim; j++)
            {
                ct->tfs_popul[i]         += ct->tfs_vector[i][j]*ct->tfs_diab[j];
                ct->tfs_popul[i+ct->dim] += ct->tfs_vector[i][j]*ct->tfs_diab[j+ct->dim];

            }
            dot_product += SQR(ct->tfs_popul[i])+SQR(ct->tfs_popul[i+ct->dim]);
        }
        for (i = 0; i < ct->tfl_num_of_states; i++) //account for limit basis
        {
            ct->tfs_popul[i] /= sqrt(dot_product); ct->tfs_popul[i+ct->dim] /= sqrt(dot_product);
        }

        for (j = 0; j < ct->dim; j++) //diab expansion may have changed due to inexact projection
        {
            ct->tfs_diab[j]         = 0.0;
            ct->tfs_diab[j+ct->dim] = 0.0;
            for (i = 0; i < ct->tfl_num_of_states; i++)
            {
                ct->tfs_diab[j]         += ct->tfs_popul[i]*ct->tfs_vector[i][j];
                ct->tfs_diab[j+ct->dim] += ct->tfs_popul[i+ct->dim]*ct->tfs_vector[i][j];
            }
        }


        int old_surface = ct->surface;

        //////check for trivial crossing:

//	printf("<90");
        if (fabs(ct->tfs_overlap[ct->surface][ct->surface]) < 0.92)  // check  falue 0.85 +- mapping surface
        {
            //|| (SQR(ca_im)+SQR(ca_re))/(SQR(ct->tfs_popul[ct->surface])+SQR(ct->tfs_popul[ct->surface+ct->dim]))>2.0){//1.4 tried too

            //printf("Detected trivial crossing!\n");
            ////using overlap, may use populations instead

            dot_product = 0.0;
            int sim_ind = 0;
            /*   for(k=0;k<ct->dim;k++){
                if(fabs(ct->tfs_overlap[ct->surface][k])>dot_product){
                  sim_ind=k; dot_product=fabs(ct->tfs_overlap[ct->surface][k]);
                }
                }  */
            for (k = 0; k < ct->tfl_num_of_states; k++)
            {
                if ((SQR(ct->tfs_popul[k])+SQR(ct->tfs_popul[k+ct->dim])) > dot_product)
                {
                    dot_product = (SQR(ct->tfs_popul[k])+SQR(ct->tfs_popul[k+ct->dim]));
                    sim_ind     = k;
                }
            }


            printf("Changed surface to %d from %d according to population\n",sim_ind,ct->surface);
            printf("Overlap is: %f\n",ct->tfs_overlap[ct->surface][sim_ind]);
            ct->surface = sim_ind;

        }

        ////////////////


        popul_norm = 0.;
        for (k = 0; k < ct->tfl_num_of_states; k++)
        {
            popul_norm += SQR(ct->tfs_popul[k]) + SQR(ct->tfs_popul[ct->dim + k]);
        }
        // printf("Population norm: %9.6f\n",popul_norm);
        //correct norm. RK not unitary!!!!
        popul_norm = sqrt(popul_norm);
        for (k = 0; k < ct->tfl_num_of_states; k++)
        {
            ct->tfs_popul[k] /= popul_norm; ct->tfs_popul[k+ct->dim] /= popul_norm;
        }

/////////do hopping//////////////////////

        /////calc. overall hopping prob.
        double hop_prob = 1.0-(SQR(ct->tfs_popul[ct->surface])+SQR(ct->tfs_popul[ct->surface+ct->dim]))/(SQR(ca_im)+SQR(ca_re));
        double prob_a   = SQR(ct->tfs_popul[ct->surface])+SQR(ct->tfs_popul[ct->surface+ct->dim]);

        //if(hop_prob>1e-5) printf("Dada: %f",hop_prob); //debug

        ///conventional form
        ct->surf_prob[ct->surface] = 0.;
        double hopping_prob_sum = 0.0;
        for (i = 0; i < ct->tfl_num_of_states; i++)
        {
            if (i != ct->surface)
            {
                //   ct->surf_prob[j]=prob_prefactor[j]*ct->tfs_overlap[ct->surface][j];
                ct->surf_prob[i] = 2.0*(ct->tfs_popul[ct->surface]*ct->tfs_popul[i]+ct->tfs_popul[ct->dim+ct->surface]*ct->tfs_popul[ct->dim+i])*0.5*(ct->tfs_overlap[old_surface][i]-ct->tfs_overlap[i][ct->surface])/prob_a;//ct->tfs_overlap[i][ct->surface];
                if (ct->surf_prob[i] > 0.0)
                {
                    hopping_prob_sum += ct->surf_prob[i];
                }
                if (ct->ev_adiab[i] > ct->ev_adiab[ct->surface])
                {
                    ct->surf_prob[i] *= exp(-(ct->ev_adiab[i]-ct->ev_adiab[ct->surface])/ct->fermi_kt);
                }
            }
        }

        //   if(hop_prob>hopping_prob_sum){
        //   for (j=0; j<ct->tfl_num_of_states; j++){
        //     if (j != ct->surface&&ct->surf_prob[j]>0.0) ct->surf_prob[j]*=hop_prob/hopping_prob_sum;
        //     }
        //  }
        hopping_prob_sum = 0.0;
        for (j = 0; j < ct->tfl_num_of_states; j++)
        {
            if (j != ct->surface && ct->surf_prob[j] > 0.0)
            {
                hopping_prob_sum += ct->surf_prob[j];
            }
        }


        if (hopping_prob_sum > 1.0)
        {
            // printf("Need to rescale probabilities!\n");
            for (j = 0; j < ct->tfl_num_of_states; j++)
            {
                if (j != ct->surface && ct->surf_prob[j] > 0.0)
                {
                    ct->surf_prob[j] /= hopping_prob_sum;
                }
            }
        }
        int prev_surf = ct->surface;
        ///////////////////////////////hopping routine/////////////////////////////////////////
        cumulated_probability = 0.;
        //random_number = rand()/(double)RAND_MAX;
        random_number = drand48();
        for (i = 0; i < ct->tfl_num_of_states; i++)
        {
            if (i != ct->surface && ct->surf_prob[i] > 0.)
            {
                cumulated_probability += ct->surf_prob[i];
                //if(ct->surf_prob[i]>1e-5) //debug
                //printf("Bla: %f\n",ct->surf_prob[i]); //debug
                if (cumulated_probability > random_number)
                {
                    ct->surface = i;
                    break;
                }
            }
        }
        //////////////decoherence correction/////////////////////////////////////////
       /* if (!ct->use_strong_dec)
        {     		reduce this option and use always just the DEC way 16.08.18 dholub	*/
            dot_product = 0.0;
            for (i = 0; i < ct->tfl_num_of_states; i++)
            {
                if (i != ct->surface)
                {
                    double decay_time = TFS_DECAY_CONSTANT/fabs(ct->ev_adiab[i]-ct->ev_adiab[ct->surface]); //dholub changed TFS_DECAY_CONSTANT in charge_trans.h in testsystem without larger effects JUMP
                    ct->tfs_popul[i]         *= exp(-ct->rk_timestep/decay_time);
                    ct->tfs_popul[i+ct->dim] *= exp(-ct->rk_timestep/decay_time);
                    dot_product              += SQR(ct->tfs_popul[i])+SQR(ct->tfs_popul[i+ct->dim]);
                }
            }
            double curr_surf_fac = sqrt((1.0-dot_product)/(SQR(ct->tfs_popul[ct->surface])+SQR(ct->tfs_popul[ct->dim+ct->surface])));
            ct->tfs_popul[ct->surface]         *= curr_surf_fac;
            ct->tfs_popul[ct->surface+ct->dim] *= curr_surf_fac;

            //recalc. diab. from adiab after decoherence corr.
            for (i = 0; i < ct->dim; i++) //bad loop order, accept for now
            {
                ct->tfs_diab[i] = 0.0; ct->tfs_diab[i+ct->dim] = 0.0;
                for (j = 0; j < ct->tfl_num_of_states; j++)
                {
                    ct->tfs_diab[i]         += ct->tfs_popul[j]*ct->tfs_vector[j][i];
                    ct->tfs_diab[i+ct->dim] += ct->tfs_popul[j+ct->dim]*ct->tfs_vector[j][i];
                }
            }

        //}
       /* else
        {
            if (prev_surf != ct->surface)
            {
                for (i = 0; i < ct->tfl_num_of_states; i++)
                {
                    ct->tfs_popul[i]         = 0.0;
                    ct->tfs_popul[i+ct->dim] = 0.0;
                }
                ct->tfs_popul[ct->surface] = 1.0;
            }

            //recalc. diab. from adiab after decoherence corr.

            for (i = 0; i < ct->dim; i++) //bad loop order, accept for now

            {
                ct->tfs_diab[i] = 0.0; ct->tfs_diab[i+ct->dim] = 0.0;
                for (j = 0; j < ct->tfl_num_of_states; j++)
                {
                    ct->tfs_diab[i]         += ct->tfs_popul[j]*ct->tfs_vector[j][i];
                    ct->tfs_diab[i+ct->dim] += ct->tfs_popul[j+ct->dim]*ct->tfs_vector[j][i];
                }
            }

        } */
        //////////////////////////////////////////////////////////////////////////////


    }//End loop over dt_step

    ct->rk_timestep *= 100.0;


    ////////////////Some output:
    printf("TFL Info: ");
    //if(prev_surf!=ct->surface){ printf("Did hopp this step! New surface: %d!\n",ct->surface); }
    //else printf("Did NOT hopp this step!\n");

    //printf("Overlap between active states at prev. and current step: %f %f\n",ct->tfs_overlap[ct->surface][prev_surf],ct->tfs_overlap[old_surface][ct->surface]);

    /////////////////////// tried to comment dholub 1.2019
    ////update system //////////////////////////////////////////////////////////////////7

    //////////determine which sites need to be included
    ct->tfl_num_of_states_old = ct->tfl_num_of_states;
    for (i = 0; i < ct->dim; i++)
    {
        ct->tfl_is_in_system_old[i] = ct->tfl_is_in_system[i];                //copy current to old
    }
    ct->tfl_num_of_states = 0;
    for (i = 0; i < ct->dim; i++)
    {
        //calc. <phi_surface|i>
        double phixi = 0.;
        for (j = 0; j < ct->dim; j++)
        {
            phixi += ct->tfs_vector[ct->surface][j]*ct->hamiltonian[i][j];
        }
        phixi = fabs(phixi);
        //calc. <phi_surface|H|phi_surface>-<i|H|j|i>
        double tfl_en_diff = fabs(ct->ev_adiab[ct->surface]-ct->hamiltonian[i][i]);
        if (phixi/tfl_en_diff > TFL_RC && fabs(phixi) > 0.001)
        {
            ct->tfl_is_in_system[i] = 1; ct->tfl_num_of_states++;
        }
        else
        {
            ct->tfl_is_in_system[i] = 0;
        }
    }
    int tfl_did_sys_change = 0;
    for (i = 0; i < ct->dim; i++)
    {
        if (ct->tfl_is_in_system[i] != ct->tfl_is_in_system_old[i])
        {
            tfl_did_sys_change = 1; break;
        }
    }
    //recalc. Eigenstates and choose appropriate surface
    if (tfl_did_sys_change)
    {

        //rebuild hamiltonian with new set of sites
        k = 0;
        for (i = 0; i < ct->dim; i++)
        {
            if (ct->tfl_is_in_system[i])
            {
                l = 0;
                for (j = 0; j < ct->dim; j++)
                {
                    if (ct->tfl_is_in_system[j])
                    {
                        ham[k+l*ct->tfl_num_of_states] = ct->hamiltonian[i][j];
                        ham[l+k*ct->tfl_num_of_states] = ct->hamiltonian[j][i];
                        l++;
                    }
                }
                k++;
            }
        }

//////////////////////////////////////////

        dsyev(ct->tfl_num_of_states, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);

        printf("System changed! WF before system adaption was: \n");
        double *old_surf_vec = (double*)malloc(sizeof(double)*ct->dim); //vec. to store surf. state.
        for (i = 0; i < ct->dim; i++)
        {
            old_surf_vec[i] = ct->tfs_vector[ct->surface][i];
        }
        for (i = 0; i < ct->dim; i++)
        {
            printf("%f ", old_surf_vec[i]);
        }
        printf("\n");

        for (i = 0; i < ct->tfl_num_of_states; i++)
        {
            k = 0;
            for (j = 0; j < ct->dim; j++)
            {
                if (ct->tfl_is_in_system[j])
                {
                    ct->tfs_vector[i][j] = ham[i * ct->tfl_num_of_states + k]; k++;
                }
                else
                {
                    ct->tfs_vector[i][j] = 0;
                }

            }
        }
        for (; i < ct->dim; i++)
        {
            for (j = 0; j < ct->dim; j++)
            {
                ct->tfs_vector[i][j] = 0;
            }
        }



        /////find which new state is surface by selecting the one with max overlap
        double maximum_overlap = 0.;
        int    max_index       = 0;
        for (i = 0; i < ct->tfl_num_of_states; i++)
        {
            dot_product = 0;
            for (k = 0; k < ct->dim; k++)
            {
                dot_product += old_surf_vec[k]*ct->tfs_vector[i][k];
            }
            //dot_product=fabs(dot_product);
            if (fabs(dot_product) > maximum_overlap)
            {
                maximum_overlap = fabs(dot_product); max_index = i;
            }

        }///////////////
        ct->surface = max_index;
        //if(dot_product<0.) for(j=0;j<ct->dim;j++) ct->tfs_vector[ct->surface][j]*=-1.0; //check sign
        printf("Surface in new system is: %d \n", ct->surface);



        /////update populations
        dot_product = 0.0;
        for (i = 0; i < ct->tfl_num_of_states; i++)
        {
            ct->tfs_popul[i]         = 0.0;
            ct->tfs_popul[i+ct->dim] = 0.0;
            for (j = 0; j < ct->dim; j++)
            {
                ct->tfs_popul[i]         += ct->tfs_vector[i][j]*ct->tfs_diab[j];
                ct->tfs_popul[i+ct->dim] += ct->tfs_vector[i][j]*ct->tfs_diab[j+ct->dim];
            }
            dot_product += SQR(ct->tfs_popul[i])+SQR(ct->tfs_popul[i+ct->dim]);
        }
        for (i = 0; i < ct->tfl_num_of_states; i++) //fix norm as projection may not be exact
        {
            ct->tfs_popul[i]         /= sqrt(dot_product);
            ct->tfs_popul[i+ct->dim] /= sqrt(dot_product);
        }
        for (j = 0; j < ct->dim; j++) //diab expansion may have changed due to inexact projection
        {
            ct->tfs_diab[j]         = 0.0;
            ct->tfs_diab[j+ct->dim] = 0.0;
            for (i = 0; i < ct->tfl_num_of_states; i++)
            {
                ct->tfs_diab[j]         += ct->tfs_popul[i]*ct->tfs_vector[i][j];
                ct->tfs_diab[j+ct->dim] += ct->tfs_popul[i+ct->dim]*ct->tfs_vector[i][j];
            }
        }



        free(old_surf_vec); //clear mem.

        //take care of populations //collapse wave function for now
        // for(i=0;i<2*ct->dim;i++) ct->tfs_popul[i]=0.0;
        //ct->tfs_popul[ct->surface]=1.0;
        ///

    }
    /////end recalc.

////////////////////end update system


    // assign the new wave function to the eigenvector ct->surface

    printf("The new wavefunction is:\n");
    for (i = 0; i < ct->dim; i++)
    {
        ct->wf[i] = ct->tfs_vector[ct->surface][i];
        printf("%f ", ct->wf[i]);
    }
    printf("\n");
/*
    printf("TFL Info: System did ");
    if (tfl_did_sys_change == 1)
    {
        printf("change.\n");
    }
    else
    {
        printf("not change.\n");
    }
*/
//    old end of comment dholub 1.2019 */

    printf("Current Surface: %d\n", ct->surface);
//debug //print eigenvectors
//  printf("Eigenvectors: \n");
// for(i=0;i<ct->dim;i++){
// for(j=0;j<ct->dim;j++) printf("%9.5f ",ct->tfs_vector[i][j]);
// printf("\n");
// }



    ////store hamiltonian as old_ham
    for (i = 0; i < SQR(ct->dim); i++)
    {
        ct->tfl_old_ham[i] = *(ct->hamiltonian[0]+i);
    }

    return 0;
}



int do_persico_diabatic_sfhopping(charge_transfer_t *ct, dftb_broyden_t *broyd, double *fermi_coeff, FILE *f, FILE *f2)
{
    int    i, j, k, step;
    double old_energy, energy, energy1, energy2, *ham, fermi_energy, fermi_upper, fermi_lower, q_sum;
    double prob_factor;
    double random_number;                      // in the interval (0,1)
    double cumulated_probability;              // sum of ct->surf_prob[0..k]
    double popul_norm;                         // sum of SQR(ct->tfs_popul) except SQR(ct->tfs_popul[ct->surface])
    double decay_time, current_surface_factor; // in the decoherence algorithm
    double dot_product;

    ham = ct->hamiltonian_adiab;

    // calculate the charges from the wave function
    for (i = 0; i < ct->dim; i++)
    {
        ct->q_act[i] = ct->q_old[i] = SQR(ct->wf[i]);
    }

    // copy as the first vector

    energy = 1.e10;

    // SCC cycle
    for (step = 0;; step++)
    {
        printf("%d\n", step);
        // too many iterations? break and indicate failure
        if (step >= MAXITER_BROYDEN)
        {
            fprintf(f, "Broyden/Fermi failed to converge!");
            return 1;
        }

        // construct the self-consistent Hamiltonian
        for (i = 0; i < ct->dim; i++)
        {
            ham[i + i * ct->dim] = ct->hamiltonian[i][i];
            for (j = 0; j < ct->dim; j++)
            {
                ham[i + i * ct->dim] += ct->hubbard[i][j] * ct->q_act[j];
                if (i != j)
                {
                    ham[i + j * ct->dim] = ct->hamiltonian[i][j];
                }
            }
        }

        // diagonalize this Hamiltonian
        dsyev(ct->dim, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);

        // FERMI DISTRIBUTION + BROYDEN MIXING
        // mix all eigenvectors with appropriate coefficients

        // calculate the "Fermi energy", store the Fermi coefficients in fermi_coeff
        fermi_lower = ct->ev_adiab[0] - 0.1;
        fermi_upper = ct->ev_adiab[ct->dim - 1] + 0.1;
        while (fermi_upper - fermi_lower > FERMI_CONVERG)
        {
            fermi_energy = 0.5 * (fermi_lower + fermi_upper);
            if (fermi_coef_sum(fermi_energy, ct->dim, ct->ev_adiab, fermi_coeff, ct->fermi_kt) > 1.0)
            {
                fermi_upper = 0.5 * (fermi_lower + fermi_upper);
            }
            else
            {
                fermi_lower = 0.5 * (fermi_lower + fermi_upper);
            }
        }

        // mix the eigenvectors, with the appropriate coefficients
        for (i = 0; i < ct->dim; i++) // i runs over the sites
        {
            ct->q_act[i] = 0.e0;
            for (j = 0; j < ct->dim; j++) // j runs over the individual eigenvectors
            {
                ct->q_act[i] += fermi_coeff[j] * SQR(ham[j * ct->dim + i]);
            }
        }

        // calculate the energy and check convergence
        old_energy = energy;
        energy1    = energy2 = 0.0;
        // indices i and j run over the sites
        for (i = 0; i < ct->dim; i++)
        {
            for (j = 0; j < ct->dim; j++)
            {
                // TB Hamiltonian
                // index k runs over the adiabatic states (eigenvectors)
                for (k = 0; k < ct->dim; k++)
                {
                    energy1 += fermi_coeff[k] * ham[k * ct->dim + i] * ham[k * ct->dim + j] * ct->hamiltonian[i][j];
                }
                // Hubbard / gamma terms
                energy2 += 0.5 * ct->hubbard[i][j] * ct->q_act[i] * ct->q_act[j];
            }
        }
        energy = energy1 + energy2;

        if (step > 1 && fabs(energy - old_energy) < BROYDEN_SCFTOL)
        {
            //fprintf(f, "Broyd %2d ", step);
            fprintf(f, " %2d ", step);
            break;
        }

        // mix the charges - Broyden
        broyden(step, BROYDEN_ALMIX, ct->dim, ct->q_old, ct->q_act, broyd);
        for (i = 0; i < ct->dim; i++)
        {
            ct->q_act[i] = ct->q_old[i];
        }
        // correct numerical errors caused by broyden
        q_sum = 0;
        for (i = 0; i < ct->dim; i++)
        {
            if (ct->q_act[i] < 0.0)
            {
                ct->q_act[i] = 0.0;
            }
            if (ct->q_act[i] > 1.0)
            {
                ct->q_act[i] = 1.0;
            }
            q_sum += ct->q_act[i];
        }
        for (i = 0; i < ct->dim; i++)
        {
            ct->q_act[i] = ct->q_act[i] / q_sum; //normalization


        }
    } /* end SCC cycle */

    // print energies of adiabatic states - these are pure, but obtained with non-pure (Fermi-mixed) input state
    //fprintf(f, " Eigenvalues:");
    for (i = 0; i < ct->dim; i++)
    {
        fprintf(f, " %8.5f", ct->ev_adiab[i]);
    }
    //fprintf(f, "    ");

    /* surface hopping tried according to Tully JCP 93, 1061 (1990)
     * "fewest switches" algorithm
     * diabatization according to Persico JCP 114, 10608 (2001)
     * approximate decoherence correction adapted from Persico JCP 126, 134114 (2007)
     */

    /* exp_imag_matrix - has to be written yet !!! */

    /* THIS WILL REMAIN! */
    if (ct->tfs_initialization_step)
    {
        /* let us start with tfs_vector_old == tfs_vector */
        for (i = 0; i < ct->dim; i++)
        {
            for (j = 0; j < ct->dim; j++)
            {
                ct->tfs_vector_old[i][j] = ham[i * ct->dim + j];
            }
        }
        for (i = 0; i < ct->dim; i++)
        {
            ct->ev_adiab_old[i] = ct->ev_adiab[i];
        }
    }
    else
    {
        /* push tfs_vector to tfs_vector_old */
        for (i = 0; i < ct->dim; i++)
        {
            for (j = 0; j < ct->dim; j++)
            {
                ct->tfs_vector_old[i][j] = ct->tfs_vector[i][j];
            }
        }
    }

    /* PUSH THE OLD POPULATIONS TO tfs_popul_old */
    for (i = 0; i < ct->dim; i++)
    {
        ct->tfs_popul_old[i]           = ct->tfs_popul[i];
        ct->tfs_popul_old[ct->dim + i] = ct->tfs_popul[ct->dim + i];
    }

    /* THIS WILL REMAIN, TOO
       MAYBE JUST DO NOT PERFORM THE CHECKING FOR SIGN CHANGE! */
    /* enter new eigenvectors to tfs_vector */
    for (i = 0; i < ct->dim; i++)
    {
        /* calculate the dot product with tfs_vector_old[i] */
        //  dot_product = 0.;
        //  for (k=0; k<ct->sites; k++)
        //    dot_product += ct->tfs_vector_old[i][k] * ham[i * ct->sites + k];
        /* for positive dot_product, take the vector as it is
         * for negative dot_product, take (-1) * the vector
         */
        //  for (k=0; k<ct->sites; k++)
        //    ct->tfs_vector[i][k] = dot_product > 0 ? ham[i * ct->sites + k] : - ham[i * ct->sites + k];
        for (k = 0; k < ct->dim; k++)
        {
            ct->tfs_vector[i][k] = ham[i * ct->dim + k];
        }
    }

    /* THIS WILL REMAIN - tfs_overlap[i][j] corresponds to matrix T[i][j] ! */
    /* calculate tfs_overlap[i][j] = <tfs_vector_old[i] | tfs_vector[j]> */
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_overlap[i][j] = 0.;
            for (k = 0; k < ct->dim; k++)
            {
                ct->tfs_overlap[i][j] += ct->tfs_vector_old[i][k] * ct->tfs_vector[j][k];
            }
        }
    }

    /* print out the overlap matrix */
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            fprintf(f2, "%9.5f", ct->tfs_overlap[i][j]);
        }
    }
    fprintf(f2, "\n");

    /* DO THE PERSICO STUFF HERE!
     * BEGINNING FROM THE SECOND STEP OF THE SIMULATION
     */
    if (ct->tfs_initialization_step)
    {
        ct->tfs_initialization_step = 0;
    }
    else
    {

        /* SET UP THE AVERAGE DIABATIC HAMILTONIAN Z = 1/2 * (E(0) + T(dt) * E(dt) * T^t(dt) */
        for (i = 0; i < ct->dim; i++)
        {
            for (j = 0; j < ct->dim; j++)
            {
                if (i == j)
                {
                    ct->per_diab_hamiltonian[i][i] = ct->ev_adiab_old[i] / 2.;
                }
                else
                {
                    ct->per_diab_hamiltonian[i][j] = 0.;
                }
                for (k = 0; k < ct->dim; k++)
                {
                    ct->per_diab_hamiltonian[i][j] += ct->tfs_overlap[i][k] * ct->ev_adiab[k] * ct->tfs_overlap[j][k] / 2.;
                }
            }
        }

        /* CONSTRUCT THE PROPAGATION OPERATOR exp[-i*Z*dt] ! */
        exp_imag_matrix(ct->per_diab_hamiltonian, ct->per_propag_operator, ct->rk_timestep, ct->dim, ct->per_arrays);

        /* CONSTRUCT THE TRANSFORMATION MATRIX U = T^t * exp[-i*Z*dt] */
        for (i = 0; i < ct->dim; i++)
        {
            for (j = 0; j < ct->dim; j++)
            {
                ct->per_transformator[i][j][0] = 0.;
                ct->per_transformator[i][j][1] = 0.;
                for (k = 0; k < ct->dim; k++)
                {
                    ct->per_transformator[i][j][0] += ct->tfs_overlap[k][i] * ct->per_propag_operator[k][j][0];
                    ct->per_transformator[i][j][1] += ct->tfs_overlap[k][i] * ct->per_propag_operator[k][j][1];
                }
            }
        }

        /* OBTAIN THE NEW POPULATIONS */
        for (i = 0; i < ct->dim; i++)
        {
            ct->tfs_popul[i] = ct->tfs_popul[ct->dim + i] = 0.;
            for (j = 0; j < ct->dim; j++)
            {
                ct->tfs_popul[i]             += ct->per_transformator[i][j][0] * ct->tfs_popul_old[j]
                    - ct->per_transformator[i][j][1] * ct->tfs_popul_old[ct->dim + j];
                ct->tfs_popul[ct->dim + i] += ct->per_transformator[i][j][0] * ct->tfs_popul_old[ct->dim + j]
                    + ct->per_transformator[i][j][1] * ct->tfs_popul_old[j];
            }
        }

        /* print out new populations */
        for (k = 0; k < ct->dim; k++)
        {
            fprintf(f, "%9.6f", SQR(ct->tfs_popul[k]) + SQR(ct->tfs_popul[ct->dim + k]));
        }

        /* we are on ct->surface at the moment
         * now, go over all of the other states,
         * and calculate switching probability */

        /* CONTINUE HERE WITH SWITCHING PROBABILITIES! */
        ct->surf_prob[ct->surface] = -1.; /* to skip this surface in the test */
        if (SQR(ct->per_transformator[ct->surface][ct->surface][0]) + SQR(ct->per_transformator[ct->surface][ct->surface][1]) > 0.999999)
        {
            /* the population of current state does not really change
             * numerical problems are to be expected
             * therefore, skip the calculation of probabilities
             * no surface hop will be performed in this step!
             */
            for (j = 0; j < ct->dim; j++)
            {
                if (j != ct->surface)
                {
                    ct->surf_prob[j] = -0.1;
                }
            }
        }
        else
        {
            /* otherwise, there is change in the population
             * and it is safe to calculate surface-hopping probabilities!
             */
            prob_factor = (SQR(ct->tfs_popul[ct->surface]) + SQR(ct->tfs_popul[ct->dim + ct->surface]) - SQR(ct->tfs_popul_old[ct->surface]) - SQR(ct->tfs_popul_old[ct->dim + ct->surface]))
                / ( SQR(ct->tfs_popul[ct->surface]) + SQR(ct->tfs_popul[ct->dim + ct->surface])
                    - ( (ct->per_transformator[ct->surface][ct->surface][0] * ct->tfs_popul_old[ct->surface] - ct->per_transformator[ct->surface][ct->surface][1] * ct->tfs_popul_old[ct->dim + ct->surface])
                        * ct->tfs_popul[ct->surface]
                        + (ct->per_transformator[ct->surface][ct->surface][0] * ct->tfs_popul_old[ct->dim + ct->surface] + ct->per_transformator[ct->surface][ct->surface][1] * ct->tfs_popul_old[ct->surface])
                        * ct->tfs_popul[ct->dim + ct->surface] ) );
            for (j = 0; j < ct->dim; j++)
            {
                if (j != ct->surface)
                {
                    ct->surf_prob[j] = -((ct->per_transformator[ct->surface][j][0] * ct->tfs_popul_old[j] - ct->per_transformator[ct->surface][j][1] * ct->tfs_popul_old[ct->dim + j]) * ct->tfs_popul[ct->surface]
                                         + (ct->per_transformator[ct->surface][j][0] * ct->tfs_popul_old[ct->dim + j] + ct->per_transformator[ct->surface][j][1] * ct->tfs_popul_old[j]) * ct->tfs_popul[ct->dim + ct->surface])
                        * prob_factor;
                }
            }
            /* ct->surf_prob[j] = -2 * (ct->tfs_popul[ct->surface] * ct->tfs_popul[j] + ct->tfs_popul[ct->sites + ct->surface] * ct->tfs_popul[ct->sites + j])
                                  / (SQR(ct->tfs_popul[ct->surface]) + SQR(ct->tfs_popul[ct->sites + ct->surface]))
             * ct->tfs_overlap[j][ct->surface]; */
        }

        /* write out the probabilities */
        for (k = 0; k < ct->dim; k++)
        {
            fprintf(f, "%10.6f", ct->surf_prob[k]);
        }

        /* generate a random number */
        random_number         = rand()/(double)RAND_MAX;
        cumulated_probability = 0.;

        /* and determine the surface to continue on */
        for (j = 0; j < ct->dim; j++)
        {
            if (j != ct->surface && ct->surf_prob[j] > 0.)
            {
                cumulated_probability += ct->surf_prob[j];
                if (cumulated_probability > random_number)
                {
                    ct->surface = j;
                    break;
                }
            }
        }

        /* decoherence correction
         * proposed by Truhlar JCTC 1, 527 (2005)
         * recommended by Persico 126, 134114 (2007)
         * simplified to:
         * tau_km = hbar / |E_k - E_m| * TFS_DECAY_CONSTANT
         * TFS_DECAY_CONSTANT chosen to be 3. (should be quite insensitive, though)
         * /

           popul_norm = 0.;
           for (j=0; j<ct->sites; j++) if (j != ct->surface) {
           decay_time = TFS_DECAY_CONSTANT / fabs(ct->ev_adiab[j] - ct->ev_adiab[ct->surface]);
           ct->tfs_popul[j]             *= exp( - ct->rk_timestep / decay_time);
           ct->tfs_popul[ct->sites + j] *= exp( - ct->rk_timestep / decay_time);
           popul_norm += SQR(ct->tfs_popul[j]) + SQR(ct->tfs_popul[ct->sites + j]);
           }
           current_surface_factor = sqrt((1. - popul_norm) / (SQR(ct->tfs_popul[ct->surface]) + SQR(ct->tfs_popul[ct->sites + ct->surface])));
           ct->tfs_popul[ct->surface]             *= current_surface_factor;
           ct->tfs_popul[ct->sites + ct->surface] *= current_surface_factor;
         */

    } /* end if initialization step else */

    /* assign the new wave function to the eigenvector ct->surface */
    for (i = 0; i < ct->dim; i++)
    {
        ct->wf[i] = ct->tfs_vector[ct->surface][i];
    }

    /* print out current state */
    fprintf(f, " %d", ct->surface);
    for (i = 0; i < ct->dim; i++)
    {
        fprintf(f, "%8.4f", SQR(ct->wf[i]));
    }
    fprintf(f, "\n");

    /* push the current ev_adiab to ev_adiab_old, to be used in the next step */
    for (i = 0; i < ct->dim; i++)
    {
        ct->ev_adiab_old[i] = ct->ev_adiab[i];
    }

    return 0;
}

int do_persico_diabatic_sfhopping_new(charge_transfer_t *ct, dftb_broyden_t *broyd, double *fermi_coeff, FILE *f, FILE *f2, FILE *f3)
{
    int             i, j, k, l, m, n, step;
    double          old_energy, energy, energy1, energy2, *ham, fermi_energy, fermi_upper, fermi_lower;
    double          prob_factor;
    double          random_number;                      // in the interval (0,1)
    double          cumulated_probability;              // sum of ct->surf_prob[0..k]
    double          popul_norm;                         // sum of SQR(ct->tfs_popul) except SQR(ct->tfs_popul[ct->surface])
    double          decay_time, current_surface_factor; // in the decoherence algorithm
    double          dot_product, tot_prob;
    double          dot_prod_r, dot_prod_i;
    /* the names of these variables are taken from Newton-X, file sh_locdiab.f90, subroutine trans_ld */
    int             surf, surf_old; /* originally: istati, istati_old */
    double          popk0, popk1, rekk, rekl, rnum, rden;
    double         *acoer0, *acoei0, *acoer1, *acoei1;
    twodoubles      ac0, ac1, vkk, vkl;
    float           maxval;
    const float     eps_switch = 0.5;//1.e-9; take always state with better overlap. note: overlap should be close to 1 or 0. if offdiag is far from 0 smaller time step is needed
    static double **pop_flow = NULL;
    static int      lwork = 0;
    static double  *hdx = NULL, *ex = NULL, *work = NULL, **rotr = NULL, **roti = NULL, **rotsr = NULL, **rotsi = NULL, **rotr_new = NULL, **roti_new = NULL;
    double          timetau, ld_sk, ld_ck, ld_fac;
    int             ld_istep;
    const int       ld_nstep = 20;

    ham      = ct->hamiltonian_adiab;
    n        = ct->dim;
    surf_old = surf = ct->surface;
    acoer0   = ct->tfs_popul_old;
    acoei0   = ct->tfs_popul_old + n;
    acoer1   = ct->tfs_popul;
    acoei1   = ct->tfs_popul + n;

    if (pop_flow == NULL)
    {
        snew(pop_flow, n);
        snew(pop_flow[0], SQR(n));
        for (i = 1; i < n; i++)
        {
            pop_flow[i] = pop_flow[0] + i * n;
        }
    }

    if (hdx == NULL)
    {
        snew(hdx, SQR(n));
    }

    if (ex == NULL)
    {
        snew(ex, n);
    }

    if (work == NULL)
    {
        lwork = (n > 3) ? SQR(n) : 5;
        snew(work, lwork);
    }

    if (rotr == NULL)
    {
        snew(rotr, n);
        snew(rotr[0], SQR(n));
        for (i = 1; i < n; i++)
        {
            rotr[i]  = rotr[0] + i * n;
        }
        snew(roti, n);
        snew(roti[0], SQR(n));
        for (i = 1; i < n; i++)
        {
            roti[i]  = roti[0] + i * n;
        }
        snew(rotsr, n);
        snew(rotsr[0], SQR(n));
        for (i = 1; i < n; i++)
        {
            rotsr[i] = rotsr[0] + i * n;
        }
        snew(rotsi, n);
        snew(rotsi[0], SQR(n));
        for (i = 1; i < n; i++)
        {
            rotsi[i] = rotsi[0] + i * n;
        }
        snew(rotr_new, n);
        snew(rotr_new[0], SQR(n));
        for (i = 1; i < n; i++)
        {
            rotr_new[i]  = rotr_new[0] + i * n;
        }
        snew(roti_new, n);
        snew(roti_new[0], SQR(n));
        for (i = 1; i < n; i++)
        {
            roti_new[i]  = roti_new[0] + i * n;
        }
    }

    // calculate the charges from the wave function
    // and copy as the first vector
    // NO! set to zero instead - perhaps improve efficiency...
    for (i = 0; i < ct->dim; i++)
    {
        ct->q_act[i] = ct->q_old[i] = 0.0; //alex: test for identical dimers with a large spacial separation
    }
    // ct->q_act[i] = ct->q_old[i] = SQR(ct->wf[i]);
    // NO! take them from the previous iteration

    energy = 1.e10;

    // SCC cycle
    for (step = 0;; step++)
    {

        // too many iterations? break and indicate failure
        if (step >= MAXITER_BROYDEN)
        {
            fprintf(f, "Broyden/Fermi failed to converge!");
            return 1;
        }

        // construct the self-consistent Hamiltonian
        for (i = 0; i < ct->dim; i++)
        {
            ham[i + i * ct->dim] = ct->hamiltonian[i][i];
            for (j = 0; j < ct->dim; j++)
            {
                ham[i + i * ct->dim] += ct->hubbard[i][j] * ct->q_act[j];
                if (i != j)
                {
                    ham[i + j * ct->dim] = ct->hamiltonian[i][j];
                }
            }
        }

        // diagonalize this Hamiltonian
        dsyev(ct->dim, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);

        // FERMI DISTRIBUTION + BROYDEN MIXING
        // mix all eigenvectors with appropriate coefficients

        // calculate the "Fermi energy", store the Fermi coefficients in fermi_coeff
        fermi_lower = ct->ev_adiab[0] - 0.1;
        fermi_upper = ct->ev_adiab[ct->dim - 1] + 0.1;
        while (fermi_upper - fermi_lower > FERMI_CONVERG)
        {
            fermi_energy = 0.5 * (fermi_lower + fermi_upper);
            if (fermi_coef_sum(fermi_energy, ct->dim, ct->ev_adiab, fermi_coeff, ct->fermi_kt) > 1.0)
            {
                fermi_upper = 0.5 * (fermi_lower + fermi_upper);
            }
            else
            {
                fermi_lower = 0.5 * (fermi_lower + fermi_upper);
            }
        }

        // mix the eigenvectors, with the appropriate coefficients
        for (i = 0; i < ct->dim; i++) // i runs over the sites
        {
            ct->q_act[i] = 0.e0;
            for (j = 0; j < ct->dim; j++) // j runs over the individual eigenvectors
            {
                ct->q_act[i] += fermi_coeff[j] * SQR(ham[j * ct->dim + i]);
            }
        }
/*
    for (i=0; i<ct->dim; i++){
    printf("qact: %lf\n", ct->q_act[i]);
    printf("fermicoeff: %lf\n", fermi_coeff[i]);
    for (j=0; j<ct->dim; j++){
        printf("ham %d %d : %lf\n", i,j, ham[j * ct->dim + i]);
    }
    }
 */
        // calculate the energy and check convergence
        old_energy = energy;
        energy1    = energy2 = 0.0;
        // indices i and j run over the sites
        for (i = 0; i < ct->dim; i++)
        {
            for (j = 0; j < ct->dim; j++)
            {
                // TB Hamiltonian
                // index k runs over the adiabatic states (eigenvectors)
                for (k = 0; k < ct->dim; k++)
                {
                    energy1 += fermi_coeff[k] * ham[k * ct->dim + i] * ham[k * ct->dim + j] * ct->hamiltonian[i][j];
                }
                // Hubbard / gamma terms
                energy2 += 0.5 * ct->hubbard[i][j] * ct->q_act[i] * ct->q_act[j];
            }
        }
        energy = energy1 + energy2;

        if (step > 1 && fabs(energy - old_energy) < BROYDEN_SCFTOL)
        {
            //fprintf(f, "Broyd %2d ", step);
            fprintf(f, " %2d ", step);
            break;
        }

        // mix the charges - Broyden
        broyden(step, BROYDEN_ALMIX, ct->dim, ct->q_old, ct->q_act, broyd);
        for (i = 0; i < ct->dim; i++)
        {
            ct->q_act[i] = ct->q_old[i];
        }

    } /* end SCC cycle */

    // print energies of adiabatic states - these are pure, but obtained with non-pure (Fermi-mixed) input state
    //fprintf(f, " Eigenvalues:");
    for (i = 0; i < ct->dim; i++)
    {
        fprintf(f, "%8.5f", ct->ev_adiab[i]);
    }
    //fprintf(f, "    ");

    /* surface hopping tried according to Tully JCP 93, 1061 (1990)
     * "fewest switches" algorithm
     * diabatization according to Persico JCP 114, 10608 (2001)
     * approximate decoherence correction adapted from Persico JCP 126, 134114 (2007)
     */

    /* exp_imag_matrix - has to be written yet !!! */

    /* 1. set up the state vectors */

    /*   a. back up state vectors from the previous step, tfs_vector_old*/

    if (ct->tfs_initialization_step)
    {
        /* let us start with tfs_vector_old == tfs_vector */
        for (i = 0; i < ct->dim; i++)
        {
            for (j = 0; j < ct->dim; j++)
            {
                ct->tfs_vector_old[i][j] = ham[i * ct->dim + j];
            }
        }
        for (i = 0; i < ct->dim; i++)
        {
            ct->ev_adiab_old[i] = ct->ev_adiab[i];
        }

        //ALEX: project input-WF on adiabatic states. i.e. take diabatic state as initial condiation instead of lowest adiabatic state.
        //for (i=0; i<ct->dim; i++) // for every surface
        //  dot_product = 0.;
        //  for (k=0; k<ct->dim; k++) // for every FO
        //    dot_product += ct->tfs_vector_old[i][k] * ct->wf[k];
        //  ct->tfs_popul[i]=dot_product;
    }
    else
    {
        /* push tfs_vector to tfs_vector_old */
        for (i = 0; i < ct->dim; i++)
        {
            for (j = 0; j < ct->dim; j++)
            {
                ct->tfs_vector_old[i][j] = ct->tfs_vector[i][j];
            }
        }
    }

    /*   b. back up populations (complex coefficients) of states from the previous step, tfs_popul_old */

    for (i = 0; i < ct->dim; i++)
    {
        ct->tfs_popul_old[i]             = ct->tfs_popul[i];
        ct->tfs_popul_old[ct->dim + i]   = ct->tfs_popul[ct->dim + i];
    }

    /*  c. save the new state vectors */

    for (i = 0; i < ct->dim; i++)
    {
        /* calculate the dot product with tfs_vector_old[i],
           to check if we shall invert the sign
         */
        dot_product = 0.;
        for (k = 0; k < ct->dim; k++)
        {
            dot_product += ct->tfs_vector_old[i][k] * ham[i * ct->dim + k];
        }
        /* for positive dot_product, take the vector as it is
         * for negative dot_product, take (-1) * the vector
         */
        for (k = 0; k < ct->dim; k++)
        {
            ct->tfs_vector[i][k] = dot_product > 0 ? ham[i * ct->dim + k] : -ham[i * ct->dim + k];
        }
        //ct->tfs_vector[i][k] =  ham[i * ct->dim + k];
    }

    /* 2. calculate the overlap of state vectors
            tfs_overlap[i][j] is the matrix T[i][j] = <tfs_vector_old[i] | tfs_vector[j]>
     */

    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_overlap[i][j] = 0.;
            for (k = 0; k < ct->dim; k++)
            {
                ct->tfs_overlap[i][j] += ct->tfs_vector_old[i][k] * ct->tfs_vector[j][k];
            }
        }
    }

    /* print out the overlap matrix */
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            fprintf(f2, "%9.5f", ct->tfs_overlap[i][j]);
        }
    }
    fprintf(f2, "\n");

    /* print out the state vectors */
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            fprintf(f3, "%9.5f", ct->tfs_vector[i][j]);
        }
    }
    fprintf(f3, "\n");

    /* check if the matrix is unitary - in our case orthogonal (real matrix)
       printf("Testing the orthogonality of T - this should be a unit matrix:\n");
       for (i=0; i<n; i++) {
       for (j=0; j<n; j++) {
        dot_product = 0.;
        for (k=0; k<n; k++)
          dot_product += ct->tfs_overlap[i][k] * ct->tfs_overlap[j][k];
        printf(" %9.6f", dot_product);
       }
       printf("\n");
       }
     */

    /* 3. do the propagation of the wave function / populations
            this is equivalent to the propag_ld function in NewtonX */

    /*    a. do not do it in the first step of the simulation */

    if (ct->tfs_initialization_step)
    {
        ct->tfs_initialization_step = 0;
    }
    else
    {

        /*    b. set up the diabatic hamiltonian at the end of the time step Hdia(dt) = T(dt) * E(dt) * T^t(dt) : Eq. B7 */
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                ct->per_diab_hamiltonian[i][j] = 0.;
                for (k = 0; k < n; k++)
                {
                    ct->per_diab_hamiltonian[i][j] += ct->tfs_overlap[i][k] * ct->ev_adiab[k] * ct->tfs_overlap[j][k];
                }
            }
        }

        /*    c. divide the propagation into smaller elementary steps */
        timetau = ct->rk_timestep / (double) ld_nstep;
        for (ld_istep = 0; ld_istep < ld_nstep; ld_istep++)
        {
            /* matrix Hdx                 = (Hdia(istep*timetau) + Hdia((istep+1)*tau)) / 2
               matrix Hdia(istep*timetau) = E(0) + (2*istep-1)/(2*ld_nstep) * (Hdia(dt) - E(0))
             */
            ld_fac = (double) (2. * ld_istep + 1) / (double) (2 * ld_nstep);
            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    hdx[i + j * n] = timetau * (ld_fac * ct->per_diab_hamiltonian[i][j] + (i == j ? (1. - ld_fac) * ct->ev_adiab_old[i] : 0.));
                }
            }
            /* eigenproblem of exp[-i * Hdx * timetau] */
            dsyev(n, hdx, ex, work, lwork);
            /* construct the rotation matrix (real + imaginary parts) for the elementary step */
            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    rotsr[i][j] = 0.;
                    rotsi[i][j] = 0.;
                }
            }
            for (k = 0; k < n; k++)
            {
                ld_sk = sin(ex[k]);
                ld_ck = cos(ex[k]);
                for (i = 0; i < n; i++)
                {
                    for (j = 0; j < n; j++)
                    {
                        rotsr[i][j] += hdx[i + k * n] * hdx[j + k * n] * ld_ck;
                        rotsi[i][j] -= hdx[i + k * n] * hdx[j + k * n] * ld_sk;
                    }
                }
            }
            /* accumulate the overall rotation matrix (real + imaginary) */
            if (ld_istep == 0)
            {
                for (i = 0; i < n; i++)
                {
                    for (j = 0; j < n; j++)
                    {
                        rotr[i][j] = rotsr[i][j];
                        roti[i][j] = rotsi[i][j];
                    }
                }
            }
            else
            {
                /* tw  = rotsr * rotr
                   tx  = rotsi * roti
                   ty  = rotsr * roti
                   tz  = rotsi * rotr
                   rotr = tw - tx
                   roti = ty + tz
                 */
                for (i = 0; i < n; i++)
                {
                    for (j = 0; j < n; j++)
                    {
                        rotr_new[i][j] = 0.;
                        roti_new[i][j] = 0.;
                        for (k = 0; k < n; k++)
                        {
                            rotr_new[i][j] += rotsr[i][k] * rotr[k][j] - rotsi[i][k] * roti[k][j];
                            roti_new[i][j] += rotsr[i][k] * roti[k][j] + rotsi[i][k] * rotr[k][j];
                        }
                    }
                }
                for (i = 0; i < n; i++)
                {
                    for (j = 0; j < n; j++)
                    {
                        rotr[i][j] = rotr_new[i][j];
                        roti[i][j] = roti_new[i][j];
                    }
                }
            }
            /*    d. update the coefficients in the adiabatic representation */
            for (i = 0; i < n; i++)
            {
                ct->tfs_popul[i]     = 0.;
                ct->tfs_popul[i + n] = 0.;
                for (j = 0; j < n; j++)
                {
                    for (k = 0; k < n; k++)
                    {
                        ct->tfs_popul[i]     += ct->tfs_overlap[j][i] *
                            (ct->tfs_popul_old[k] * rotr[j][k] - ct->tfs_popul_old[k + n] * roti[j][k]);
                        ct->tfs_popul[i + n] += ct->tfs_overlap[j][i] *
                            (ct->tfs_popul_old[k] * roti[j][k] + ct->tfs_popul_old[k + n] * rotr[j][k]);
                    }
                }
            }
        }
        /*
           // print out the total rotation matrix
           printf("===============================================\n");
           printf("Rotation matrix for coefficients in the adiabatic representation:\n");
           for (i=0; i<n; i++) {
           for (j=0; j<n; j++)
            printf("  %9.6f + %9.6fi", rotr[i][j], roti[i][j]);
           printf("\n");
           }
           printf("Testing the unitarity of the rotation matrix - this should be a unit matrix:\n");
           for (i=0; i<n; i++) {
           for (j=0; j<n; j++) {
            dot_prod_r = dot_prod_i = 0.;
            for (k=0; k<n; k++) {
              dot_prod_r += rotr[i][k] * rotr[j][k] + roti[i][k] * roti[j][k];
              dot_prod_i += roti[i][k] * rotr[j][k] - rotr[i][k] * roti[j][k];
           }
            printf("%10.6f%10.6f", dot_prod_r, dot_prod_i);
           }
           printf("\n");
           }
         */
        /*
           printf("TFS_POPUL_NEW");
           for (i=0; i<n; i++)
           printf("%10.6f%10.6f",  ct->tfs_popul[i],  ct->tfs_popul[i + n]);
           printf("\n");
         */

        /* FOR TESTING PURPOSES, COMPARE WITH OUR PREVIOUS IMPLEMENTATION!

           // average diabatic Hamiltonian: Z = 1/2 * (E(0) + T(dt) * E(dt) * T^t(dt)
           for (i=0; i<ct->sites; i++) {
           for (j=0; j<ct->sites; j++) {
            if (i==j) {
              ct->per_diab_hamiltonian[i][i] = ct->ev_adiab_old[i] / 2.;
            } else {
              ct->per_diab_hamiltonian[i][j] = 0.;
            }
            for (k=0; k<ct->sites; k++)
              ct->per_diab_hamiltonian[i][j] += ct->tfs_overlap[i][k] * ct->ev_adiab[k] * ct->tfs_overlap[j][k] / 2.;
           }
           }
           // construct the propagation operator exp[-i*Z*dt] !
           exp_imag_matrix(ct->per_diab_hamiltonian, ct->per_propag_operator, ct->rk_timestep, ct->sites, ct->per_arrays);
           // print out the propagation operator
           printf("===============================================\n");
           printf("Original alternative - the propagation operator:\n");
           for (i=0; i<n; i++) {
           for (j=0; j<n; j++)
            printf("  %9.6f + %9.6fi", ct->per_propag_operator[i][j][0], ct->per_propag_operator[i][j][1]);
           printf("\n");
           }
           printf("===============================================\n");

           printf("Testing the unitarity of the propagation operator:\n");
           for (i=0; i<n; i++) {
           for (j=0; j<n; j++) {
            dot_prod_r = dot_prod_i = 0.;
            for (k=0; k<n; k++) {
              dot_prod_r += ct->per_propag_operator[i][k][0] * ct->per_propag_operator[j][k][0] + ct->per_propag_operator[i][k][1] * ct->per_propag_operator[j][k][1];
              dot_prod_i += ct->per_propag_operator[i][k][1] * ct->per_propag_operator[j][k][0] - ct->per_propag_operator[i][k][0] * ct->per_propag_operator[j][k][1];
           }
            printf("%10.6f%10.6f", dot_prod_r, dot_prod_i);
           }
           printf("\n");
           }
         */

        /* construct the transformation matrix U = T^t * exp[-i*Z*dt]
           for (i=0; i<ct->sites; i++)
           for (j=0; j<ct->sites; j++) {
            ct->per_transformator[i][j][0] = 0.;
            ct->per_transformator[i][j][1] = 0.;
            for (k=0; k<ct->sites; k++) {
              ct->per_transformator[i][j][0] += ct->tfs_overlap[k][i] * ct->per_propag_operator[k][j][0];
              ct->per_transformator[i][j][1] += ct->tfs_overlap[k][i] * ct->per_propag_operator[k][j][1];
            }
           }
         */
        /* then, new populations would be obtained as
           for (i=0; i<n; i++) {
           ct->tfs_popul[i] = ct->tfs_popul[n + i] = 0.;
           for (j=0; j<n; j++) {
            ct->tfs_popul[i]     += ct->per_transformator[i][j][0] * ct->tfs_popul_old[j]
                                  - ct->per_transformator[i][j][1] * ct->tfs_popul_old[n + j];
            ct->tfs_popul[n + i] += ct->per_transformator[i][j][0] * ct->tfs_popul_old[n + j]
         + ct->per_transformator[i][j][1] * ct->tfs_popul_old[j];
           }
           }
         */
        /*
           printf("TFS_POPUL_OLD");
           for (i=0; i<n; i++)
           printf("%10.6f%10.6f",  ct->tfs_popul[i],  ct->tfs_popul[i + n]);
           printf("\n");
         */
        /* END TESTING PURPOSES */



        /* GO ON HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */



        /* print out new populations */
        for (k = 0; k < ct->dim; k++)
        {
            fprintf(f, "%9.6f", SQR(ct->tfs_popul[k]) + SQR(ct->tfs_popul[ct->dim + k]));
        }

        /* we are on ct->surface at the moment
         * now, go over all of the other states,
         * and calculate switching probability */

        /* the following is inspired/taken from Newton-X, file sh_locdiab.f90, subroutine trans_ld */

        /* save the value U_kk */
        //vkk[0] = ct->per_transformator[surf][surf][0];
        //vkk[1] = ct->per_transformator[surf][surf][1];
        vkk[0] = vkk[1] = 0.;
        for (j = 0; j < n; j++)
        {
            vkk[0] += ct->tfs_overlap[j][surf] * rotr[j][surf];
            vkk[1] += ct->tfs_overlap[j][surf] * roti[j][surf];
        }

        //printf("\n CHECK VKK: %9.6f %9.6f VS %9.6f %9.6f\n",
        //  ct->per_transformator[surf][surf][0], ct->per_transformator[surf][surf][1], vkk[0], vkk[1]);

        /* has a switch of states just occured? */
        surf_old = surf;
        if (fabs(ct->tfs_overlap[surf][surf]) < eps_switch)
        {
            printf("state switching has occurred! Overlap is only: %f  \n", ct->tfs_overlap[surf][surf]);
            for (j = 0; j < n; j++)
            {
                ct->surf_prob[j] = 0.;
            }
            maxval = 0.;
            for (j = 0; j < n; j++)
            {
                if (fabs((float) ct->tfs_overlap[j][surf_old]) > maxval)
                {
                    surf   = j;
                    maxval = fabs((float) ct->tfs_overlap[j][surf_old]);
                }
            }
            printf("state %d probably became state %d! Overlap: %f  \n", surf_old, surf, ct->tfs_overlap[surf][surf_old]);
        }

        /* rename / evaluate simple variables */
        ac0[0] = acoer0[surf_old];
        ac0[1] = acoei0[surf_old];
        ac1[0] = acoer1[surf];
        ac1[1] = acoei1[surf];
        popk0  = SQR(ac0[0]) + SQR(ac0[1]);
        popk1  = SQR(ac1[0]) + SQR(ac1[1]);
        /* rekk = U_kk * A_k(0) * A_k^*(dt) */
        rekk = vkk[0] * (ac0[0]*ac1[0] + ac0[1]*ac1[1]) + vkk[1] * (ac0[0]*ac1[1] - ac0[1]*ac1[0]);
        //printf(" CHECK REKK: %9.6f\n", rekk);
        rnum = (popk1 - popk0) / popk0;
        rden = popk1 - rekk;

        /* if population of occupied surface increases,
         * then no hopping!
         */
        if (popk1 > popk0)
        {
            for (j = 0; j < n; j++)
            {
                ct->surf_prob[j] = -1.;
            }
            /* otherwise, calculate the probabilities */
        }
        else
        {
            for (j = 0; j < n; j++)
            {
                /* no hopping from state 'surf' to itself, obviously */
                if (j == surf)
                {
                    ct->surf_prob[j] = -1.;
                    /* calculate the probability for all of the states except 'surf' */
                }
                else
                {
                    //vkl[0] = ct->per_transformator[surf][j][0];
                    //vkl[1] = ct->per_transformator[surf][j][1];
                    vkl[0] = vkl[1] = 0.;
                    for (m = 0; m < n; m++)
                    {
                        vkl[0] += ct->tfs_overlap[m][surf] * rotr[m][j];
                        vkl[1] += ct->tfs_overlap[m][surf] * roti[m][j];
                    }
                    rekl             = vkl[0] * (acoer0[j]*ac1[0] + acoei0[j]*ac1[1]) + vkl[1] * (acoer0[j]*ac1[1] - acoei0[j]*ac1[0]);
                    ct->surf_prob[j] = -rekl * rnum / rden;
                }
            }
        }
        /* probabilities are done at this point */
        for (k = 0; k < n; k++)
        {
            //vkk[0] = ct->per_transformator[k][k][0];
            //vkk[1] = ct->per_transformator[k][k][1];
            vkk[0] = vkk[1] = 0.;
            for (j = 0; j < n; j++)
            {
                vkk[0] += ct->tfs_overlap[j][k] * rotr[j][k];
                vkk[1] += ct->tfs_overlap[j][k] * roti[j][k];
            }
            ac0[0] = acoer0[k];
            ac0[1] = acoei0[k];
            ac1[0] = acoer1[k];
            ac1[1] = acoei1[k];
            popk0  = SQR(ac0[0]) + SQR(ac0[1]);
            popk1  = SQR(ac1[0]) + SQR(ac1[1]);
            /* rekk = U_kk * A_k(0) * A_k^*(dt) */
            rekk = vkk[0] * (ac0[0]*ac1[0] + ac0[1]*ac1[1]) + vkk[1] * (ac0[0]*ac1[1] - ac0[1]*ac1[0]);
            rnum = (popk1 - popk0) / popk0;
            rden = popk1 - rekk;
            for (l = 0; l < n; l++)
            {
                //vkl[0] = ct->per_transformator[k][l][0];
                //vkl[1] = ct->per_transformator[k][l][1];
                vkl[0] = vkl[1] = 0.;
                for (m = 0; m < n; m++)
                {
                    vkl[0] += ct->tfs_overlap[m][k] * rotr[m][l];
                    vkl[1] += ct->tfs_overlap[m][k] * roti[m][l];
                }
                rekl           = vkl[0] * (acoer0[l]*ac1[0] + acoei0[l]*ac1[1]) + vkl[1] * (acoer0[l]*ac1[1] - acoei0[l]*ac1[0]);
                pop_flow[k][l] = -rekl * rnum * popk0 / rden;
                //printf(" %8.5f", pop_flow[k][l]);
            }
        }
        printf("\n");

        tot_prob = 0.;
        for (j = 0; j < n; j++)
        {
            tot_prob += (ct->surf_prob[j] > 0. ? ct->surf_prob[j] : 0.);
        }
        if (tot_prob > 1.)
        {
            for (j = 0; j < n; j++)
            {
                if (ct->surf_prob[j] > 0.)
                {
                    ct->surf_prob[j] /= tot_prob;
                }
            }
        }

        /* END OF PROBABILITIES! */

        /* write out the probabilities */
        for (k = 0; k < ct->dim; k++)
        {
            fprintf(f, "%10.6f", ct->surf_prob[k]);
        }

        /* generate a random number */
        random_number         = rand()/(double)RAND_MAX;
        cumulated_probability = 0.;

        /* and determine the surface to continue on */
        for (j = 0; j < ct->dim; j++)
        {
            if (j != ct->surface && ct->surf_prob[j] > 0.)
            {
                cumulated_probability += ct->surf_prob[j];
                if (cumulated_probability > random_number)
                {
                    ct->surface = j;
                    break;
                }
            }
        }

        /* decoherence correction
         * proposed by Truhlar JCTC 1, 527 (2005)
         * recommended by Persico 126, 134114 (2007)
         * simplified to:
         * tau_km = hbar / |E_k - E_m| * TFS_DECAY_CONSTANT
         * TFS_DECAY_CONSTANT chosen to be 3. (should be quite insensitive, though)
         */

        if (ct->decoherence)
        {
            popul_norm = 0.;
            for (j = 0; j < n; j++)
            {
                if (j != ct->surface)
                {
                    decay_time                    = TFS_DECAY_CONSTANT / fabs(ct->ev_adiab[j] - ct->ev_adiab[ct->surface]);
                    ct->tfs_popul[j]             *= exp( -ct->rk_timestep / decay_time);
                    ct->tfs_popul[ct->dim + j]   *= exp( -ct->rk_timestep / decay_time);
                    popul_norm                   += SQR(ct->tfs_popul[j]) + SQR(ct->tfs_popul[n + j]);
                }
            }
            current_surface_factor                  = sqrt((1. - popul_norm) / (SQR(ct->tfs_popul[ct->surface]) + SQR(ct->tfs_popul[n + ct->surface])));
            ct->tfs_popul[ct->surface]             *= current_surface_factor;
            ct->tfs_popul[ct->dim + ct->surface]   *= current_surface_factor;
        }

    } /* end if initialization step else */

    /* assign the new wave function to the eigenvector ct->surface */
    for (i = 0; i < ct->dim; i++)
    {
        ct->wf[i] = ct->tfs_vector[ct->surface][i];
    }

    /* print out current state */
    fprintf(f, " %d", ct->surface);
    for (i = 0; i < ct->dim; i++)
    {
        fprintf(f, "%8.4f", SQR(ct->wf[i]));
    }
    fprintf(f, "\n");

    /* push the current ev_adiab to ev_adiab_old, to be used in the next step */
    for (i = 0; i < ct->dim; i++)
    {
        ct->ev_adiab_old[i] = ct->ev_adiab[i];
    }

    return 0;
}




int do_prezhdo_sfhopping(charge_transfer_t *ct, dftb_broyden_t *broyd, double *fermi_coeff, FILE *f, FILE *f2, FILE *f3)
{
    int             i, j, k, l, m, n, step;
    double          old_energy, energy, energy1, energy2, *ham, fermi_energy, fermi_upper, fermi_lower;
    double          prob_factor;
    double          random_number;                      // in the interval (0,1)
    double          cumulated_probability;              // sum of ct->surf_prob[0..k]
    double          popul_norm;                         // sum of SQR(ct->tfs_popul) except SQR(ct->tfs_popul[ct->surface])
    double          decay_time, current_surface_factor; // in the decoherence algorithm
    double          dot_product, tot_prob;
    double          dot_prod_r, dot_prod_i;
    /* the names of these variables are taken from Newton-X, file sh_locdiab.f90, subroutine trans_ld */
    int             surf, surf_old; /* originally: istati, istati_old */
    double          popk0, popk1, rekk, rekl, rnum, rden;
    double         *acoer0, *acoei0, *acoer1, *acoei1;
    twodoubles      ac0, ac1, vkk, vkl;
    float           maxval;
    const float     eps_switch = 1.e-9;
    static double **pop_flow = NULL;
    static int      lwork = 0;
    static double  *hdx = NULL, *ex = NULL, *work = NULL, **rotr = NULL, **roti = NULL, **rotsr = NULL, **rotsi = NULL, **rotr_new = NULL, **roti_new = NULL;
    double          timetau, ld_sk, ld_ck, ld_fac;
    int             ld_istep;
    const int       ld_nstep = 20;
    double          diff_flow, gross_flow;
    double          a, b, ski;

    ham      = ct->hamiltonian_adiab;
    n        = ct->dim;
    surf_old = surf = ct->surface;
    acoer0   = ct->tfs_popul_old;
    acoei0   = ct->tfs_popul_old + n;
    acoer1   = ct->tfs_popul;
    acoei1   = ct->tfs_popul + n;

    if (pop_flow == NULL)
    {
        snew(pop_flow, n);
        snew(pop_flow[0], SQR(n));
        for (i = 1; i < n; i++)
        {
            pop_flow[i] = pop_flow[0] + i * n;
        }
    }

    if (hdx == NULL)
    {
        snew(hdx, SQR(n));
    }

    if (ex == NULL)
    {
        snew(ex, n);
    }

    if (work == NULL)
    {
        lwork = (n > 3) ? SQR(n) : 5;
        snew(work, lwork);
    }

    if (rotr == NULL)
    {
        snew(rotr, n);
        snew(rotr[0], SQR(n));
        for (i = 1; i < n; i++)
        {
            rotr[i]  = rotr[0] + i * n;
        }
        snew(roti, n);
        snew(roti[0], SQR(n));
        for (i = 1; i < n; i++)
        {
            roti[i]  = roti[0] + i * n;
        }
        snew(rotsr, n);
        snew(rotsr[0], SQR(n));
        for (i = 1; i < n; i++)
        {
            rotsr[i] = rotsr[0] + i * n;
        }
        snew(rotsi, n);
        snew(rotsi[0], SQR(n));
        for (i = 1; i < n; i++)
        {
            rotsi[i] = rotsi[0] + i * n;
        }
        snew(rotr_new, n);
        snew(rotr_new[0], SQR(n));
        for (i = 1; i < n; i++)
        {
            rotr_new[i]  = rotr_new[0] + i * n;
        }
        snew(roti_new, n);
        snew(roti_new[0], SQR(n));
        for (i = 1; i < n; i++)
        {
            roti_new[i]  = roti_new[0] + i * n;
        }
    }

    // calculate the charges from the wave function
    // and copy as the first vector
    // NO! set to zero instead - perhaps improve efficiency...
    for (i = 0; i < ct->dim; i++)
    {
        ct->q_act[i] = ct->q_old[i] = 0.0; //alex: test for identical dimers with a large spacial separation
    }
    // ct->q_act[i] = ct->q_old[i] = SQR(ct->wf[i]);
    // NO! take them from the previous iteration

    energy = 1.e10;
    // SCC cycle
    for (step = 0;; step++)
    {

        // too many iterations? break and indicate failure
        if (step >= MAXITER_BROYDEN)
        {
            fprintf(f, "Broyden/Fermi failed to converge!");
            return 1;
        }

        // construct the self-consistent Hamiltonian
        for (i = 0; i < ct->dim; i++)
        {
            ham[i + i * ct->dim] = ct->hamiltonian[i][i];
            for (j = 0; j < ct->dim; j++)
            {
                ham[i + i * ct->dim] += ct->hubbard[i][j] * ct->q_act[j];
                if (i != j)
                {
                    ham[i + j * ct->dim] = ct->hamiltonian[i][j];
                }
            }
        }

        // diagonalize this Hamiltonian
        dsyev(ct->dim, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);

        // FERMI DISTRIBUTION + BROYDEN MIXING
        // mix all eigenvectors with appropriate coefficients

        // calculate the "Fermi energy", store the Fermi coefficients in fermi_coeff
        fermi_lower = ct->ev_adiab[0] - 0.1;
        fermi_upper = ct->ev_adiab[ct->dim - 1] + 0.1;
        while (fermi_upper - fermi_lower > FERMI_CONVERG)
        {
            fermi_energy = 0.5 * (fermi_lower + fermi_upper);
            if (fermi_coef_sum(fermi_energy, ct->dim, ct->ev_adiab, fermi_coeff, ct->fermi_kt) > 1.0)
            {
                fermi_upper = 0.5 * (fermi_lower + fermi_upper);
            }
            else
            {
                fermi_lower = 0.5 * (fermi_lower + fermi_upper);
            }
        }

        // mix the eigenvectors, with the appropriate coefficients
        for (i = 0; i < ct->dim; i++) // i runs over the sites
        {
            ct->q_act[i] = 0.e0;
            for (j = 0; j < ct->dim; j++) // j runs over the individual eigenvectors
            {
                ct->q_act[i] += fermi_coeff[j] * SQR(ham[j * ct->dim + i]);
            }
        }
/*
    for (i=0; i<ct->dim; i++){
    printf("qact: %lf\n", ct->q_act[i]);
    printf("fermicoeff: %lf\n", fermi_coeff[i]);
    for (j=0; j<ct->dim; j++){
        printf("ham %d %d : %lf\n", i,j, ham[j * ct->dim + i]);
    }
    }
 */
        // calculate the energy and check convergence
        old_energy = energy;
        energy1    = energy2 = 0.0;
        // indices i and j run over the sites
        for (i = 0; i < ct->dim; i++)
        {
            for (j = 0; j < ct->dim; j++)
            {
                // TB Hamiltonian
                // index k runs over the adiabatic states (eigenvectors)
                for (k = 0; k < ct->dim; k++)
                {
                    energy1 += fermi_coeff[k] * ham[k * ct->dim + i] * ham[k * ct->dim + j] * ct->hamiltonian[i][j];
                }
                // Hubbard / gamma terms
                energy2 += 0.5 * ct->hubbard[i][j] * ct->q_act[i] * ct->q_act[j];
            }
        }
        energy = energy1 + energy2;

        if (step > 1 && fabs(energy - old_energy) < BROYDEN_SCFTOL)
        {
            //fprintf(f, "Broyd %2d ", step);
            fprintf(f, " %2d ", step);
            break;
        }

        // mix the charges - Broyden
        broyden(step, BROYDEN_ALMIX, ct->dim, ct->q_old, ct->q_act, broyd);
        for (i = 0; i < ct->dim; i++)
        {
            ct->q_act[i] = ct->q_old[i];
        }

    } /* end SCC cycle */

    // print energies of adiabatic states - these are pure, but obtained with non-pure (Fermi-mixed) input state
    //fprintf(f, " Eigenvalues:");
    for (i = 0; i < ct->dim; i++)
    {
        fprintf(f, "%8.5f", ct->ev_adiab[i]);
    }
    //fprintf(f, "    ");

    /* surface hopping tried according to Tully JCP 93, 1061 (1990)
     * "fewest switches" algorithm
     * diabatization according to Persico JCP 114, 10608 (2001)
     * approximate decoherence correction adapted from Persico JCP 126, 134114 (2007)
     */

    /* exp_imag_matrix - has to be written yet !!! */

    /* 1. set up the state vectors */

    /*   a. back up state vectors from the previous step, tfs_vector_old*/

    if (ct->tfs_initialization_step)
    {
        /* let us start with tfs_vector_old == tfs_vector */
        for (i = 0; i < ct->dim; i++)
        {
            for (j = 0; j < ct->dim; j++)
            {
                ct->tfs_vector_old[i][j] = ham[i * ct->dim + j];
            }
        }
        for (i = 0; i < ct->dim; i++)
        {
            ct->ev_adiab_old[i] = ct->ev_adiab[i];
        }

        //ALEX: project input-WF on adiabatic states. i.e. take diabatic state as initial condiation instead of lowest adiabatic state.
        //for (i=0; i<ct->dim; i++) // for every surface
        //  dot_product = 0.;
        //  for (k=0; k<ct->dim; k++) // for every FO
        //    dot_product += ct->tfs_vector_old[i][k] * ct->wf[k];
        //  ct->tfs_popul[i]=dot_product;
    }
    else
    {
        /* push tfs_vector to tfs_vector_old */
        for (i = 0; i < ct->dim; i++)
        {
            for (j = 0; j < ct->dim; j++)
            {
                ct->tfs_vector_old[i][j] = ct->tfs_vector[i][j];
            }
        }
    }

    /*   b. back up populations (complex coefficients) of states from the previous step, tfs_popul_old */

    for (i = 0; i < ct->dim; i++)
    {
        ct->tfs_popul_old[i]             = ct->tfs_popul[i];
        ct->tfs_popul_old[ct->dim + i]   = ct->tfs_popul[ct->dim + i];
    }

    /*  c. save the new state vectors */
// IT IS NOT NECESSARY TO CHECK THE SIGN HERE. WE WILL CALCULATE OVERLAP AND THIS INCLUDES SIGN CHECK
    for (i = 0; i < ct->dim; i++)
    {
        /* calculate the dot product with tfs_vector_old[i],
           to check if we shall invert the sign
         */
        dot_product = 0.;
        for (k = 0; k < ct->dim; k++)// WHAT IF THERE IS CROSSING OF ADIABATIC STATES? WE SHOULD CHECK THAT STATE i IS SAME AS IN PREVIOUS STEP
        {
            dot_product += ct->tfs_vector_old[i][k] * ham[i * ct->dim + k];
        }
        /* for positive dot_product, take the vector as it is
         * for negative dot_product, take (-1) * the vector
         */
        for (k = 0; k < ct->dim; k++)
        {
            ct->tfs_vector[i][k] = dot_product > 0 ? ham[i * ct->dim + k] : -ham[i * ct->dim + k];
        }
        //ct->tfs_vector[i][k] = ham[i * ct->dim + k]; //we will check phase later
    }

    /* 2. calculate the overlap of state vectors
            tfs_overlap[i][j] is the matrix T[i][j] = <tfs_vector_old[i] | tfs_vector[j]>
     */

    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            ct->tfs_overlap[i][j] = 0.;
            for (k = 0; k < ct->dim; k++)
            {
                ct->tfs_overlap[i][j] += ct->tfs_vector_old[i][k] * ct->tfs_vector[j][k];
            }
        }
    }
/*
   // If T matrix is replaced by <tfs_vector_old[i] | tfs_vector[j]> then there is an additional orthogonalization in newton-X which we are missing here (see Persico JCP 114, 10608 (2001) Appendix B)
   // However, there is also reordering of the states which will be done here
   for (i=0; i<ct->dim; i++)
    for (j=0; j<ct->dim; j++) {
      a=ct->tfs_overlap[i][i]*ct->tfs_overlap[i][i] + ct->tfs_overlap[j][j]*ct->tfs_overlap[j][j];
      b=ct->tfs_overlap[i][j]*ct->tfs_overlap[i][j] + ct->tfs_overlap[i][j]*ct->tfs_overlap[j][i];
      if ( a < b )
      for (k=0; k<ct->dim; k++){
        ski = ct->tfs_overlap[k][i];
        ct->tfs_overlap[k][i] = ct->tfs_overlap[k][j];
        ct->tfs_overlap[k][j] = ski;
      }
    }
 */
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            a = ct->tfs_overlap[i][i]*ct->tfs_overlap[i][i] + ct->tfs_overlap[j][j]*ct->tfs_overlap[j][j];
            b = ct->tfs_overlap[i][j]*ct->tfs_overlap[i][j] + ct->tfs_overlap[i][j]*ct->tfs_overlap[j][i];
            if (a < b)                           //switch occurred between state i and j. we should probably take care of this
            {
                if (ct->tfs_overlap[i][j] < 0.0) //state j has changed sign
                {
                    for (k = 0; k < ct->dim; k++)
                    {
                        ct->tfs_vector[j][k] = -ct->tfs_vector[j][k];
                    }
                }
            }
            else if (i == j && ct->tfs_overlap[i][i] < 0.0)//state stays the same but has changed sign
            {
                for (k = 0; k < ct->dim; k++)
                {
                    ct->tfs_vector[i][k] = -ct->tfs_vector[i][k];
                }
            }
        }
    }



    /* print out the overlap matrix */
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            fprintf(f2, "%9.5f", ct->tfs_overlap[i][j]);
        }
    }
    fprintf(f2, "\n");

    /* print out the state vectors */
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            fprintf(f3, "%9.5f", ct->tfs_vector[i][j]);
        }
    }
    fprintf(f3, "\n");

    /* check if the matrix is unitary - in our case orthogonal (real matrix)
       printf("Testing the orthogonality of T - this should be a unit matrix:\n");
       for (i=0; i<n; i++) {
       for (j=0; j<n; j++) {
        dot_product = 0.;
        for (k=0; k<n; k++)
          dot_product += ct->tfs_overlap[i][k] * ct->tfs_overlap[j][k];
        printf(" %9.6f", dot_product);
       }
       printf("\n");
       }
     */


    /* 3. do the propagation of the wave function / populations
            this is equivalent to the propag_ld function in NewtonX */

    /*    a. do not do it in the first step of the simulation */

    if (ct->tfs_initialization_step)
    {
        ct->tfs_initialization_step = 0;
    }
    else
    {

        /*    b. set up the diabatic hamiltonian at the end of the time step Hdia(dt) = T(dt) * E(dt) * T^t(dt) : Eq. B7 */
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                ct->per_diab_hamiltonian[i][j] = 0.;
                for (k = 0; k < n; k++)
                {
                    ct->per_diab_hamiltonian[i][j] += ct->tfs_overlap[i][k] * ct->ev_adiab[k] * ct->tfs_overlap[j][k];
                }
            }
        }

        /*    c. divide the propagation into smaller elementary steps */
        timetau = ct->rk_timestep / (double) ld_nstep;
        for (ld_istep = 0; ld_istep < ld_nstep; ld_istep++)
        {
            /* matrix Hdx                 = (Hdia(istep*timetau) + Hdia((istep+1)*tau)) / 2
               matrix Hdia(istep*timetau) = E(0) + (2*istep-1)/(2*ld_nstep) * (Hdia(dt) - E(0))
             */
            ld_fac = (double) (2. * ld_istep + 1) / (double) (2 * ld_nstep);
            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    hdx[i + j * n] = timetau * (ld_fac * ct->per_diab_hamiltonian[i][j] + (i == j ? (1. - ld_fac) * ct->ev_adiab_old[i] : 0.));
                }
            }
            /* eigenproblem of exp[-i * Hdx * timetau] */
            dsyev(n, hdx, ex, work, lwork);
            /* construct the rotation matrix (real + imaginary parts) for the elementary step */
            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    rotsr[i][j] = 0.;
                    rotsi[i][j] = 0.;
                }
            }
            for (k = 0; k < n; k++)
            {
                ld_sk = sin(ex[k]);
                ld_ck = cos(ex[k]);
                for (i = 0; i < n; i++)
                {
                    for (j = 0; j < n; j++)
                    {
                        rotsr[i][j] += hdx[i + k * n] * hdx[j + k * n] * ld_ck;
                        rotsi[i][j] -= hdx[i + k * n] * hdx[j + k * n] * ld_sk;
                    }
                }
            }
            /* accumulate the overall rotation matrix (real + imaginary) */
            if (ld_istep == 0)
            {
                for (i = 0; i < n; i++)
                {
                    for (j = 0; j < n; j++)
                    {
                        rotr[i][j] = rotsr[i][j];
                        roti[i][j] = rotsi[i][j];
                    }
                }
            }
            else
            {
                /* tw  = rotsr * rotr
                   tx  = rotsi * roti
                   ty  = rotsr * roti
                   tz  = rotsi * rotr
                   rotr = tw - tx
                   roti = ty + tz
                 */
                for (i = 0; i < n; i++)
                {
                    for (j = 0; j < n; j++)
                    {
                        rotr_new[i][j] = 0.;
                        roti_new[i][j] = 0.;
                        for (k = 0; k < n; k++)
                        {
                            rotr_new[i][j] += rotsr[i][k] * rotr[k][j] - rotsi[i][k] * roti[k][j];
                            roti_new[i][j] += rotsr[i][k] * roti[k][j] + rotsi[i][k] * rotr[k][j];
                        }
                    }
                }
                for (i = 0; i < n; i++)
                {
                    for (j = 0; j < n; j++)
                    {
                        rotr[i][j] = rotr_new[i][j];
                        roti[i][j] = roti_new[i][j];
                    }
                }
            }
            /*    d. update the coefficients in the adiabatic representation */
            for (i = 0; i < n; i++)
            {
                ct->tfs_popul[i]     = 0.;
                ct->tfs_popul[i + n] = 0.;
                for (j = 0; j < n; j++)
                {
                    for (k = 0; k < n; k++)
                    {
                        ct->tfs_popul[i]     += ct->tfs_overlap[j][i] *
                            (ct->tfs_popul_old[k] * rotr[j][k] - ct->tfs_popul_old[k + n] * roti[j][k]);
                        ct->tfs_popul[i + n] += ct->tfs_overlap[j][i] *
                            (ct->tfs_popul_old[k] * roti[j][k] + ct->tfs_popul_old[k + n] * rotr[j][k]);
                    }
                }
            }
        }

        /* print out new populations */
        for (k = 0; k < ct->dim; k++)
        {
            fprintf(f, "%9.6f", SQR(ct->tfs_popul[k]) + SQR(ct->tfs_popul[ct->dim + k]));
        }

        /* we are on ct->surface at the moment
         * now, go over all of the other states,
         * and calculate switching probability */


        /* has a switch of states just occured? */ // I guess this shouldn't happen after I corrected ordering before diabatic hamiltonian is constructed
        surf_old = surf;
        if (fabs(ct->tfs_overlap[surf][surf]) < eps_switch)
        {
            printf("state switching has occurred! Overlap is only: %f  \n", ct->tfs_overlap[surf][surf]);
            for (j = 0; j < n; j++)
            {
                //ct->surf_prob[j] = 0.; why was this here?
                maxval = 0.;
            }
            for (j = 0; j < n; j++)
            {
                if (fabs((float) ct->tfs_overlap[j][surf_old]) > maxval)
                {
                    surf        = j;
                    ct->surface = j; //I guess this was missing
                    maxval      = fabs((float) ct->tfs_overlap[j][surf_old]);
                }
            }
            printf("state %d probably became state %d! Overlap: %f  \n", surf_old, surf, ct->tfs_overlap[surf][surf_old]);
        }



        /* if population of occupied surface increases,
         * then no hopping!
         */
        ac0[0] = acoer0[surf_old];
        ac0[1] = acoei0[surf_old];
        ac1[0] = acoer1[surf];
        ac1[1] = acoei1[surf];
        popk0  = SQR(ac0[0]) + SQR(ac0[1]);
        popk1  = SQR(ac1[0]) + SQR(ac1[1]);
        if (popk1 > popk0)
        {
            for (j = 0; j < n; j++)
            {
                ct->surf_prob[j] = -1.;
            }

            /* otherwise, calculate the probabilities */
        }
        else
        {
            /* calculate the gross flow */
            gross_flow = 0.;
            for (j = 0; j < n; j++)
            {
                diff_flow = SQR(acoer1[j]) + SQR(acoei1[j]) - (SQR(acoer0[j]) + SQR(acoei0[j]));
                if (diff_flow > 0.)
                {
                    gross_flow += diff_flow;
                }
            }
            /* go over the states */
            for (j = 0; j < n; j++)
            {
                /* no hopping from state 'surf' to itself, obviously */
                if (j == surf)
                {
                    ct->surf_prob[j] = -1.;
                    /* calculate the probability for all of the states except 'surf' */
                }
                else
                {
                    diff_flow = SQR(acoer1[j]) + SQR(acoei1[j]) - (SQR(acoer0[j]) + SQR(acoei0[j]));
                    /* no hopping to a state with a decreased population */
                    if (diff_flow < 0.)
                    {
                        ct->surf_prob[j] = -0.1;
                        /* otherwise, calculate probability according to Prezhdo */
                    }
                    else
                    {
                        ct->surf_prob[j] = diff_flow * (popk0 - popk1) / (popk0 * gross_flow);
                    }
                }
            }
        }


        tot_prob = 0.;
        for (j = 0; j < n; j++)
        {
            tot_prob += (ct->surf_prob[j] > 0. ? ct->surf_prob[j] : 0.);
        }
        if (tot_prob > 1.)
        {
            for (j = 0; j < n; j++)
            {
                if (ct->surf_prob[j] > 0.)
                {
                    ct->surf_prob[j] /= tot_prob;
                }
            }
        }

        /* END OF PROBABILITIES! */

        /* write out the probabilities */
        for (k = 0; k < ct->dim; k++)
        {
            fprintf(f, "%10.6f", ct->surf_prob[k]);
        }

        /* generate a random number */
        random_number         = rand()/(double)RAND_MAX; //hopping switched off for testing
        cumulated_probability = 0.;

        /* and determine the surface to continue on */
        for (j = 0; j < ct->dim; j++)
        {
            if (j != ct->surface && ct->surf_prob[j] > 0.)
            {
                cumulated_probability += ct->surf_prob[j];
                if (cumulated_probability > random_number)
                {
                    ct->surface = j;
                    break;
                }
            }
        }

        /* decoherence correction
         * proposed by Truhlar JCTC 1, 527 (2005)
         * recommended by Persico 126, 134114 (2007)
         * simplified to:
         * tau_km = hbar / |E_k - E_m| * TFS_DECAY_CONSTANT
         * TFS_DECAY_CONSTANT chosen to be 3. (should be quite insensitive, though)
         */

        if (ct->decoherence)
        {
            popul_norm = 0.;
            for (j = 0; j < n; j++)
            {
                if (j != ct->surface)
                {
                    decay_time                    = TFS_DECAY_CONSTANT / fabs(ct->ev_adiab[j] - ct->ev_adiab[ct->surface]);
                    ct->tfs_popul[j]             *= exp( -ct->rk_timestep / decay_time);
                    ct->tfs_popul[ct->dim + j]   *= exp( -ct->rk_timestep / decay_time);
                    popul_norm                   += SQR(ct->tfs_popul[j]) + SQR(ct->tfs_popul[n + j]);
                }
            }
            current_surface_factor                  = sqrt((1. - popul_norm) / (SQR(ct->tfs_popul[ct->surface]) + SQR(ct->tfs_popul[n + ct->surface])));
            ct->tfs_popul[ct->surface]             *= current_surface_factor;
            ct->tfs_popul[ct->dim + ct->surface]   *= current_surface_factor;
        }

    } /* end if initialization step else */

    /* assign the new wave function to the eigenvector ct->surface */
    for (i = 0; i < ct->dim; i++)
    {
        ct->wf[i] = ct->tfs_vector[ct->surface][i];
    }

    /* print out current state */
    fprintf(f, " %d", ct->surface);
    for (i = 0; i < ct->dim; i++)
    {
        fprintf(f, "%8.4f", SQR(ct->wf[i]));
    }
    fprintf(f, "\n");

    /* push the current ev_adiab to ev_adiab_old, to be used in the next step */
    for (i = 0; i < ct->dim; i++)
    {
        ct->ev_adiab_old[i] = ct->ev_adiab[i];
    }

    return 0;
}


int do_wf_decomp_evec_fermi(charge_transfer_t *ct, dftb_broyden_t *broyd, double *fermi_coeff, FILE *f)
{
    int    i, j, k, step;
    double old_energy, energy, *ham, fermi_energy, fermi_upper, fermi_lower, dotproduct_real, dotproduct_imag, q_sum;

    ham = ct->hamiltonian_adiab;

    // calculate the charges from the wave function
    for (i = 0; i < ct->dim; i++)
    {
        ct->q_act[i] = ct->q_old[i] = SQR(ct->wf[i]);
    }

    energy = 1.e10;

    // SCC cycle
    for (step = 0;; step++)
    {

        // too many iterations? break and indicate failure
        if (step >= MAXITER_BROYDEN)
        {
            fprintf(f, "Broyden/Fermi failed to converge!\n");
            return 1;
        }

        // construct the self-consistent Hamiltonian
        for (i = 0; i < ct->dim; i++)
        {
            ham[i + i * ct->dim] = ct->hamiltonian[i][i];
            for (j = 0; j < ct->dim; j++)
            {
                ham[i + i * ct->dim] += ct->hubbard[i][j] * ct->q_act[j];
                if (i != j)
                {
                    ham[i + j * ct->dim] = ct->hamiltonian[i][j];
                }
            }
        }

        // diagonalize this Hamiltonian
        dsyev(ct->dim, ham, ct->ev_adiab, ct->work_adiab, 3*ct->dim);

        // FERMI DISTRIBUTION + BROYDEN MIXING

        // calculate the "Fermi energy", store the Fermi coefficients in fermi_coeff
        fermi_lower = ct->ev_adiab[0] - 0.1;
        fermi_upper = ct->ev_adiab[ct->dim - 1] + 0.1;
        while (fermi_upper - fermi_lower > FERMI_CONVERG)
        {
            fermi_energy = 0.5 * (fermi_lower + fermi_upper);
            //if (fermi_coef_sum(fermi_energy, ct->dim, ct->ev_adiab, fermi_coeff, ct->fermi_kt) > 1.0)
            if (fermi_coef_sum(fermi_energy, ct->dim, ct->ev_adiab, fermi_coeff, FERMI_KT) > 1.0) // ct->fermi_kt is not read for SCC
            {
                fermi_upper = 0.5 * (fermi_lower + fermi_upper);
            }
            else
            {
                fermi_lower = 0.5 * (fermi_lower + fermi_upper);
            }
        }

        // mix the eigenvectors, with the appropriate coefficients
        for (i = 0; i < ct->dim; i++) // i runs over the sites
        {
            ct->q_act[i] = 0.e0;
            for (j = 0; j < ct->dim; j++) // j runs over the individual eigenvectors
            {
                ct->q_act[i] += fermi_coeff[j] * SQR(ham[j * ct->dim + i]);
            }
        }

        // calculate the energy and check convergence
        old_energy = energy;
        energy     = 0.e0;
        // indices i and j run over the sites
        for (i = 0; i < ct->dim; i++)
        {
            for (j = 0; j < ct->dim; j++)
            {
                // TB Hamiltonian
                // index k runs over the adiabatic states (eigenvectors)
                for (k = 0; k < ct->dim; k++)
                {
                    energy += fermi_coeff[k] * ham[k * ct->dim + i] * ham[k * ct->dim + j] * ct->hamiltonian[i][j];
                }
                // Hubbard / gamma terms
                energy += 0.5 * ct->hubbard[i][j] * ct->q_act[i] * ct->q_act[j];
            }
        }

        if (step > 1 && fabs(energy - old_energy) < BROYDEN_SCFTOL)
        {
            // Broyden mixing has converged
            break;
        }

        // mix the charges - Broyden
        broyden(step, BROYDEN_ALMIX, ct->dim, ct->q_old, ct->q_act, broyd);
        for (i = 0; i < ct->dim; i++)
        {
            ct->q_act[i] = ct->q_old[i];
        }
        // correct numerical errors caused by broyden
        q_sum = 0;
        for (i = 0; i < ct->dim; i++)
        {
            if (ct->q_act[i] < 0.0)
            {
                ct->q_act[i] = 0.0;
            }
            if (ct->q_act[i] > 1.0)
            {
                ct->q_act[i] = 1.0;
            }
            q_sum += ct->q_act[i];
        }
        for (i = 0; i < ct->dim; i++)
        {
            ct->q_act[i] = ct->q_act[i] / q_sum; //normalization





        }
    } /* end SCC cycle */

    // print the energies of "adiabatic states"
    for (i = 0; i < ct->dim; i++)
    {
        fprintf(f, "%10.7f", ct->ev_adiab[i]);
    }

    // calculate the expansion coefficients
    fprintf(f, " WF expansion:");
    for (i = 0; i < ct->dim; i++)
    {
        dotproduct_real = 0.e0;
        for (j = 0; j < ct->dim; j++)
        {
            dotproduct_real += ham[i * ct->dim + j] * ct->wf[j];
        }
        dotproduct_imag = 0.e0;
        for (j = 0; j < ct->dim; j++)
        {
            dotproduct_imag += ham[i * ct->dim + j] * ct->wf[j + ct->dim];
        }
        fprintf(f, "%9.6f", SQR(dotproduct_real)+SQR(dotproduct_imag));
    }
    /*
       for (i=0; i<ct->sites; i++) {
       dotproduct_real = dotproduct_imag = 0.e0;
       for (j=0; j<ct->sites; j++)
        for (k=0; k<ct->sites; k++)
          dotproduct_squared += ham[i * ct->sites + j] * ham[i * ct->sites + k] * (ct->wf[j] * ct->wf[k] + ct->wf[j + ct->sites] * ct->wf[k + ct->sites]);
       fprintf(f, "%9.6f", dotproduct_squared);
       }
     */
    fprintf(f, "\n");


    return 0;
}

void get_delta_q(dftb_t *dftb, charge_transfer_t *ct, int i)
{
    int    j, k, lj, m, n, jofn, mj, mo;
    double q, qhelp, temp;
    //calculates atomic charge differences for hole in MO i

//for (i = 0; i < ct->sites; i++)
    for (k = 0; k < ct->site[i].homos; k++)
    {
        mo = ct->site[i].homo[k]-1;
//printf("site %d homo %d\n", i, mo);

        for (m = 0; m < dftb->phase1[i].ndim; m++)
        {
            dftb->phase1[i].qmulli[m] = 0.0;
            for (n = 0; n < dftb->phase1[i].ndim; n++)
            {
//      dftb->phase1[i].qmulli[m] += dftb->phase1[i].overl[m][n] * dftb->phase1[i].occ[mo] * dftb->phase1[i].a[m][mo] * dftb->phase1[i].a[n][mo]; //no sum over occupied MOs, only hole-MO needed for delta_q. Attention: overwrites qmulli
                dftb->phase1[i].qmulli[m] += dftb->phase1[i].overl[m][n] * dftb->phase1[i].a[m][mo] * dftb->phase1[i].a[n][mo]; //no sum over occupied MOs, only hole-MO needed for delta_q. Attention: overwrites qmulli
            }
        }

        for (j = 0; j < dftb->phase1[i].nn; j++)
        {
            q = 0.0;
            // printf("j = %d\n", j);
            for (lj = 0; lj < dftb->lmax[dftb->phase1[i].izp[j]]; lj++)
            {
                jofn = dftb->phase1[i].ind[j] + lj*lj;
                // printf("lj = %d, jofn = %d\n", lj, jofn);
                qhelp = 0.0;
                for (mj = 0; mj <= 2*lj; mj++)
                {
                    // printf("mj = %d\n", mj);
                    qhelp += dftb->phase1[i].qmulli[jofn+mj];
                }
                q += qhelp;
            }
            ct->site[i].delta_q[k][j] = (ct->is_hole_transfer == 1) ? q : -q;
            //printf("mulliken delta_q %f \n", ct->site[i].delta_q[k][j]);
        }
    }
    return;
}

void get_internal_forces(dftb_t *dftb, charge_transfer_t *ct, int site_i)
{
    int           i, j, l, nn;
    dftb_phase1_t dftb1;

    dftb1 = dftb->phase1[site_i];
    nn    = dftb1.nn;
    //sign convention: calculate change of gradient in charged system due to additional charge carrier.
    // grad_total = grad_charged-grad_neutral    F= -grad

    for (i = 0; i < nn; i++)
    {
        clear_dvec(dftb1.grad[i]);
    }
    // calculate difference of atomic hamilton shifts (= sum over gamma*delta_charge)
    // no shift due to external field. forces due to electrostatic interacions of charge carrier and MM atoms are handled by Gromacs.
    for (i = 0; i < nn; i++)
    {
        dftb1.shift[i] = 0.0;
        for (j = 0; j < nn; j++)
        {
            dftb1.shift[i] += (dftb1.qmat[j] - dftb->qzero1[dftb1.izp[j]]) * (i > j ? dftb1.gammamat[i][j] : dftb1.gammamat[j][i]);
        }
        //for(l=0; l<ct->site[site_i].homos; l++)
        //dftb1.shift[i] += -ct->site[site_i].delta_q[l][j] * ct->occupation[ct->indFO[site_i]+l] * (i>j ? dftb1.gammamat[i][j] : dftb1.gammamat[j][i]);
//printf("shift %d = %f dQ= %f\n", i, dftb1.shift[i], (dftb1.qmat[i] - dftb->qzero1[dftb1.izp[i]]) );
    }

    //force due to missing electron in HOMO
    for (i = 0; i < nn; i++)
    {
        clear_dvec(dftb1.partgrad[i]);
    }

    usual_gradient_homo(dftb, dftb1.x, dftb1.partgrad, ct, site_i);

    for (i = 0; i < nn; i++)
    {
        copy_dvec(dftb1.partgrad[i], dftb1.grad[i]);
    }

/*
   for (i=0; i<nn; i++)
    for(l=0; l<DIM; l++)
      printf("usual grad at atom i %d  QM force %lf \n", i,  -(real) HARTREE_BOHR2MD * dftb1.partgrad[i][l]);
 */
/*
   //force due to changing atomic charges by extracting electron from HOMO // now we are only calculating DFTB1 forces
   for (i=0; i<nn; i++)
    clear_dvec(dftb1.partgrad[i]);
   gamma_gradient_homo(dftb, dftb1.x, dftb1.partgrad, ct, site_i);
   for (i=0; i<nn; i++)
    dvec_inc(dftb1.grad[i], dftb1.partgrad[i]);

   for (i=0; i<nn; i++)
    for(l=0; l<DIM; l++)
      printf("gamma grad at atom i %d  QM force %lf \n", i,  -(real) HARTREE_BOHR2MD * dftb1.partgrad[i][l]);
 */

/*
   // third force term //this is not needed
   for (i=0; i<nn; i++)
    clear_dvec(dftb1.partgrad[i]);
   //  additional_gradient_homo(dftb, dftb1.x, dftb1.partgrad, ct, site_i);
   // for (i=0; i<nn; i++)
   //    dvec_inc(dftb1.grad[i], dftb1.partgrad[i]);
   for (i=0; i<nn; i++)
    for(l=0; l<DIM; l++)
      printf("additional grad at atom i %d  QM force %lf \n", i,  -(real) HARTREE_BOHR2MD * dftb1.partgrad[i][l]);
 */

    return;
}





void project_wf_on_new_basis_exact(int step, dftb_t *dftb, charge_transfer_t *ct, FILE *f_ct_project_wf, FILE*f_ct_project_wf_ref )
//void project_wf_on_new_basis(int step, dftb_t *dftb, charge_transfer_t *ct, FILE *f_ct_project_wf, FILE*f_ct_project_wf_ref )
{
    int    i, ii, j, jj, k, l, m, n;
    int    offset, iao, jao, ifo, jfo;
    double sum;
    double norm;
// Transform orthogonal eigenvectors (basis of FO Hamiltonian) from non-orthogonal fragment basis to atomic basis //
    for (l = 0; l < ct->dim; l++)
    {
        for (iao = 0; iao < dftb->phase2.norb; iao++)
        {
            dftb->orthogo.evec_ao[iao][l] = 0.0;
        }
    }
    for (l = 0; l < ct->dim; l++)
    {
        k = 0;
        for (i = 0; i < ct->sites; i++)
        {
            for (ii = 0; ii < ct->site[i].homos; ii++)
            {
                ifo = ct->site[i].homo[ii] + dftb->phase2.inf[i] - 1;
                for (iao = 0; iao < dftb->phase2.norb; iao++)
                {
                    dftb->orthogo.evec_ao[iao][l] += dftb->orthogo.sij[k + l * ct->dim] * dftb->phase2.Taf[iao][ifo];
                }
                k++;
            }
        }
    }
    if (ct->first_step)
    {
        printf("initialize state following \n");
        for (l = 0; l < ct->dim; l++)
        {
            for (iao = 0; iao < dftb->phase2.norb; iao++)
            {
                dftb->orthogo.evec_ao_old[iao][l] = dftb->orthogo.evec_ao[iao][l];
                dftb->orthogo.evec_ao_ref[iao][l] = dftb->orthogo.evec_ao[iao][l];
            }
        }
    }
/*
   // calculate overlap for basis function i
   // using combined overlap matrix with different atomic overlaps for intra and inter-fragment calculations
   //no longer needed. now dftb2.overl is already hybrid matrix
   // <evec_ao|evec_ao_old> = evec_ao^T * S_ao * evec_ao_old
   for (iao=0; iao<dftb->phase2.norb; iao++)
   for (jao=0; jao<dftb->phase2.norb; jao++)
   dftb->phase2.overl_hybrid[iao][jao]=dftb->phase2.overl[iao][jao]; // we just need here the inter fragment overlap. intra fragment overlap will be overwritten
   offset=0;
   for (i=0; i<ct->sites; i++){
   for (j=0; j < dftb->phase1[i].norb; j++) {
   for (k=0; k < dftb->phase1[i].norb; k++){
   dftb->phase2.overl_hybrid[offset+j][offset+k] = dftb->phase1[i].overl[j][k];
   dftb->phase2.overl_hybrid[offset+k][offset+j] = dftb->phase1[i].overl[k][j];
   }
   }
   offset+=dftb->phase1[i].norb;
   }
 */
/*
   offset=0;
   for (i=0; i<ct->sites; i++){
   m=0;
   for (j=0; j < dftb->phase1[i].norb; j++) {
   n=0;
   for (k=0; k < dftb->phase1[i].norb; k++){
   dftb->phase2.overl_hybrid[offset+m][offset+n] = dftb->phase1[i].overl[j][k];
   dftb->phase2.overl_hybrid[offset+n][offset+m] = dftb->phase1[i].overl[k][j];
   n++;
   }
   m++;
   }
   offset+=dftb->phase1[i].norb;
   }
 */
/*
   printf("overlap ortho/ortho (hybrid matrix)\n");
   printf( "%10d ", step);
   for (i=0; i<ct->dim; i++) {
   for (m=i; m<ct->dim; m++) printf(" %10.6f ", dftb->overl_test[i][m]);
   }
   printf( "\n");
 */
// overlap with last step //
// <evec_ao|evec_ao_old> = evec_ao^T * S_ao * evec_ao_old
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            dftb->orthogo.overlap[i][j] = 0.0;
            for (iao = 0; iao < dftb->phase2.norb; iao++)
            {
                for (jao = 0; jao < dftb->phase2.norb; jao++)
                {
                    dftb->orthogo.overlap[i][j] += dftb->orthogo.evec_ao[iao][i] * dftb->phase2.overl[iao][jao] * dftb->orthogo.evec_ao_old[jao][j];
                }
            }
        }
    }
// overlap with wf in t=0
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            dftb->orthogo.overlap_ref[i][j] = 0.0;
            for (iao = 0; iao < dftb->phase2.norb; iao++)
            {
                for (jao = 0; jao < dftb->phase2.norb; jao++)
                {
                    dftb->orthogo.overlap_ref[i][j] += dftb->orthogo.evec_ao[iao][i] * dftb->phase2.overl[iao][jao] * dftb->orthogo.evec_ao_ref[jao][j];
                }
            }
        }
    }
///*
    printf("Wave function before proj %d:\n", step);
    for (i = 0; i < ct->dim; i++)
    {
        printf("%10.7f %10.7f\n", ct->wf[i], ct->wf[i + ct->dim]);
    }
//*/
// get new wf by projecting old wf onto new basis. element i is now:
// wf_i(t2) = |fo_i(t2)><fo_i(t2)|wf(t1)>
// wf_i(t2) = sum_j |fo_i(t2)><fo_i(t2)|fo_j(t1)>*c_j
    for (i = 0; i < ct->dim; i++)
    {
        ct->wf_old[i]         = ct->wf[i];
        ct->wf_old[i+ct->dim] = ct->wf[i+ct->dim];
        ct->wf[i]             = 0.0;
        ct->wf[i+ct->dim]     = 0.0;
    }
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            ct->wf[i]         += ct->wf_old[j] * dftb->orthogo.overlap[i][j];         //real part
            ct->wf[i+ct->dim] += ct->wf_old[j+ct->dim] * dftb->orthogo.overlap[i][j]; //imaginary part
        }
    }
// scale new wavefunction (reasons: projection is not complete) //
    norm = 0;
    for (i = 0; i < ct->dim; i++)
    {
        norm += SQR(ct->wf[i]) + SQR(ct->wf[i+ct->dim]);
    }
    for (i = 0; i < 2*ct->dim; i++)
    {
        ct->wf[i] /= sqrt(norm);
    }
///*
    printf("Wave function after proj %d:\n", step);
    for (i = 0; i < ct->dim; i++)
    {
        printf("%10.7f %10.7f\n", ct->wf[i], ct->wf[i + ct->dim]);
    }
//*/
// save wave function for next step //
    for (l = 0; l < ct->dim; l++)
    {
        for (iao = 0; iao < dftb->phase2.norb; iao++)
        {
            dftb->orthogo.evec_ao_old[iao][l] = dftb->orthogo.evec_ao[iao][l];
        }
    }
// OUTPUT INTO FILES//
    fprintf(f_ct_project_wf, "%10d ", step);
    for (i = 0; i < ct->dim; i++)
    {
        for (m = i; m < ct->dim; m++)
        {
            fprintf(f_ct_project_wf, " %10.6f ", dftb->orthogo.overlap[i][m]);
        }
    }
    fprintf(f_ct_project_wf, " %10.6f ", sqrt(norm));
    fprintf(f_ct_project_wf, "\n");
    fprintf(f_ct_project_wf_ref, "%10d ", step);
    for (i = 0; i < ct->dim; i++)
    {
        for (m = i; m < ct->dim; m++)
        {
            fprintf(f_ct_project_wf_ref, " %10.6f ", dftb->orthogo.overlap_ref[i][m]);
        }
    }
    fprintf(f_ct_project_wf_ref, "\n");
    return;
}




//void project_wf_on_new_basis_approx(int step, dftb_t *dftb, charge_transfer_t *ct, FILE *f_ct_project_wf, FILE*f_ct_project_wf_ref )
void project_wf_on_new_basis(int step, dftb_t *dftb, charge_transfer_t *ct, FILE *f_ct_project_wf, FILE*f_ct_project_wf_ref )
{
    int    i, ii, j, jj, k, l, m, n;
    int    offset, iao, jao, ifo, jfo;
    double sum;
    double norm;

    // get new wf by projecting old wf onto new basis. element i will be:
    // wf_i(t2) =       |fo_i(t2)><fo_i(t2)|wf(t1)>
    // wf_i(t2) = sum_j |fo_i(t2)><fo_i(t2)|fo_j(t1)>*c_j

    // exact way would be transformation of orthogonalized basis functions |fo_i> into AO basis. Then calculate overlap with last step as:  <fo_i_ao|fo_i_ao_old> = fo_i_ao^T * S_ao * fo_i_ao_old
    // fast approximate version: use non-orthogonal FOs instead of orthogonalized ones to calculate <fo_i_ao|fo_i_ao_old>. -> is blockdiagonal matrix -> linear scaling
    // furthermore <fo_i_ao|fo_i_ao_old> was already calculated by check_and_invert_orbital_phase().



    // save reference WF at t=0 to get the total change of the basis
    if (ct->first_step)
    {
        for (i = 0; i < ct->sites; i++)
        {
            for (iao = 0; iao < dftb->phase1[i].norb; iao++)
            {
                for (jao = 0; jao < dftb->phase1[i].norb; jao++)
                {
                    dftb->phase1[i].a_ref[iao][jao] = dftb->phase1[i].a[iao][jao];
                }
            }
        }
    }

    // overlap with reference wf in t=0
    for (k = 0; k < ct->sites; k++)
    {
        for (l = 0; l < ct->site[k].homos; l++)
        {
            for (m = 0; m < ct->site[k].homos; m++)
            {
                i = dftb->phase2.ihomo[k]+l;
                j = dftb->phase2.ihomo[k]+m;
                dftb->orthogo.overlap_ref[i][j] = 0.0;
                for (iao = 0; iao < dftb->phase1[k].norb; iao++)
                {
                    for (jao = 0; jao < dftb->phase1[k].norb; jao++)
                    {
                        ifo = ct->site[k].homo[l]-1;
                        jfo = ct->site[k].homo[m]-1;
                        dftb->orthogo.overlap_ref[i][j] += dftb->phase1[k].a[iao][ifo] * dftb->phase1[k].overl[iao][jao] * dftb->phase1[k].a_ref[jao][jfo];
                    }
                }
            }
        }
    }


    // calculate overlap with last step //
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            dftb->orthogo.overlap[i][j] = 0.0;
        }
    }
    for (k = 0; k < ct->sites; k++)
    {
        for (l = 0; l < ct->site[k].homos; l++)
        {
            for (m = 0; m < ct->site[k].homos; m++)
            {
                i = dftb->phase2.ihomo[k]+l;
                j = dftb->phase2.ihomo[k]+m;
                dftb->orthogo.overlap[i][j] = ct->site[k].overlap[l][m];
            }
        }
    }

/*
         printf("Wave function before proj %d:\n", step);
         for (i=0; i<ct->dim; i++)
           printf("%10.7f %10.7f\n", ct->wf[i], ct->wf[i + ct->dim]);
 */


    for (i = 0; i < ct->dim; i++)
    {
        ct->wf_old[i]         = ct->wf[i];
        ct->wf_old[i+ct->dim] = ct->wf[i+ct->dim];
        ct->wf[i]             = 0.0;
        ct->wf[i+ct->dim]     = 0.0;
    }
    for (i = 0; i < ct->dim; i++)
    {
        for (j = 0; j < ct->dim; j++)
        {
            ct->wf[i]         += ct->wf_old[j] * dftb->orthogo.overlap[i][j];         //real part
            ct->wf[i+ct->dim] += ct->wf_old[j+ct->dim] * dftb->orthogo.overlap[i][j]; //imaginary part
        }
    }


    // scale new wavefunction (reasons: projection is not complete) //
    norm = 0;
    for (i = 0; i < ct->dim; i++)
    {
        norm += SQR(ct->wf[i]) + SQR(ct->wf[i+ct->dim]);
    }
    printf("norm changed to %f due to projection. Rescaling wavefunction\n", norm);
    for (i = 0; i < 2*ct->dim; i++)
    {
        ct->wf[i] /= sqrt(norm);
    }


/*
         printf("Wave function  after proj %d:\n", step);
         for (i=0; i<ct->dim; i++)
           printf("%10.7f %10.7f\n", ct->wf[i], ct->wf[i + ct->dim]);
 */



    // OUTPUT INTO FILES//
    fprintf(f_ct_project_wf, "%10d ", step);
    for (i = 0; i < ct->dim; i++)
    {
        for (m = i; m < ct->dim; m++)
        {
            fprintf(f_ct_project_wf, " %10.6f ", dftb->orthogo.overlap[i][m]);
        }
    }
    fprintf(f_ct_project_wf, " %10.6f ", sqrt(norm));
    fprintf(f_ct_project_wf, "\n");

    fprintf(f_ct_project_wf_ref, "%10d ", step);
    for (i = 0; i < ct->dim; i++)
    {
        for (m = i; m < ct->dim; m++)
        {
            fprintf(f_ct_project_wf_ref, " %10.6f ", dftb->orthogo.overlap_ref[i][m]);
        }
    }
    fprintf(f_ct_project_wf_ref, "\n");


    return;
}


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


void search_starting_site(matrix state_box, t_mdatoms *mdatoms, dftb_t *dftb, charge_transfer_t *ct, rvec *x_ct, char *slko_path, gmx_mtop_t *top_global, rvec *gromacs_x)
{
    int    original_nsites, i, j, lowest_site, counter;
    double Elowest, Esite;
    rvec  *dummy_x;
    matrix dummy_matrix;
    Elowest = 1000.0;

    /* perform DFTB calculations on every site to find the lowest energy
       this will determine the starting site */
    printf("Scanning for site with lowest energy.\n scanning site:\n");
    original_nsites = ct->sites;
    ct->sites       = 1;
    for (i = 0; i < ct->pool_size; i++)
    {
        ct->site[0] = ct->pool_site[i];
        //printf("%d resnr %d addr %p\n",i, ct->pool_site[i].resnr, &(ct->pool_site[i].resnr) );
        //printf("%d resnr %d addr %p\n",i, ct->site[0].resnr, &(ct->site[0].resnr));
        prepare_charge_transfer(state_box, mdatoms, dftb, ct, x_ct);
        run_dftb1(ct, dftb, 0);
        for (j = 0; j < ct->site[0].homos; j++)
        {
            Esite  = dftb->phase1[0].ev[ct->site[0].homo[j]-1];
            Esite *= ct->is_hole_transfer ? -1.0 : 1.0;
            if (Esite < Elowest)
            {
                Elowest     = Esite;
                lowest_site = i;
                copy_dvec(dftb->phase1[0].com, ct->coc);
                printf("found residue %d  MO %d  E = %f ha  COM[Angstrom] = %f %f %f\n", ct->pool_site[i].resnr, ct->pool_site[i].homo[j], Esite,  ct->coc[XX] * 0.52, ct->coc[YY] * 0.52, ct->coc[ZZ] * 0.52);
            }
        }
    }
    printf("starting residue is poolsite %d res %d  COC[Angstrom] = %f %f %f\n", lowest_site, ct->pool_site[lowest_site].resnr, ct->coc[XX] * 0.52, ct->coc[YY] * 0.52, ct->coc[ZZ] * 0.52);
    ct->sites = original_nsites;


    /* find nearest neighborst of the starting site */
    for (i = 0; i < ct->sites; i++)
    {
        ct->site[i] = ct->pool_site[i];
    }
    i = 0;
    while (adapt_QMzone(ct, x_ct, mdatoms, top_global, state_box, gromacs_x))
    {
        i++;
    }
    printf("Initially adapted QM zone %d times.\n", i);

    /* set starting wavefunction */
    counter = 0;
    for (i = 0; i < ct->sites; i++)
    {
        for (j = 0; j < ct->site[i].homos; j++)
        {
            if (ct->site[i].resnr == ct->pool_site[lowest_site].resnr && (j == ct->site[i].homos-1))
            {
                ct->wf[counter] = 1.0; ct->wf[counter+ct->dim] = 0.0;
            }
            else
            {
                ct->wf[counter] = 0.0; ct->wf[counter+ct->dim] = 0.0;
            }
            counter++;
        }
    }
    printf("starting wavefunction:\n");
    for (i = 0; i < ct->dim; i++)
    {
        printf("%10.7f %10.7f\n", ct->wf[i], ct->wf[i + ct->dim]);
    }


    return;
}

int adapt_QMzone(charge_transfer_t *ct, rvec *x_ct, t_mdatoms *mdatoms, gmx_mtop_t *top_global, matrix state_box, rvec *gromacs_x)
{
    int    i, j, l, k, counter, nearest_inactive = -1, farthest_active = -1, adapted = 0;
    dvec   bond, com, coord, masscoord;
    double bond_length, best_inactive_dist, best_active_dist, tot_occ, norm;
    rvec   box_center, dx;

    /* get center of mass for every site */
    counter = 0;
    for (i = 0; i < ct->pool_size; i++)
    {
        clear_dvec(com);
        for (j = 0; j < ct->pool_site[i].atoms; j++)
        {
            for (k = 0; k < 3; k++)
            {
                coord[k] =  x_ct[ct->pool_site[i].atom[j]][k] * NM_TO_BOHR;//take here x_ct for wich molecules were made whole again. in gromacs_x all atoms are put in the box.
            }
            dsvmul(mdatoms->massT[ct->pool_site[i].atom[j]], coord, masscoord);
            dvec_inc(com, masscoord);
        }
        dsvmul(ct->adapt_inv_tot_mass, com, ct->pool_site[i].com);
        //printf("residue %d  COM[Angstrom]: %f %f %f\n", ct->pool_site[i].resnr, ct->pool_site[i].com[XX] * 0.52, ct->pool_site[i].com[YY] * 0.52, ct->pool_site[i].com[ZZ] * 0.52);
    }


    /* search for nearest inactive site */
    best_inactive_dist = 100000.0;
    for (i = 0; i < ct->pool_size; i++)
    {
        if (!ct->pool_site[i].active)
        {
            dvec_sub(ct->pool_site[i].com, ct->coc, bond);
            if (dnorm(bond) < best_inactive_dist)
            {
                best_inactive_dist = dnorm(bond);
                nearest_inactive   = i;
            }
        }
    }
    //printf("nearest site is %d (res %d) dist %f\n",nearest_inactive,ct->pool_site[nearest_inactive].resnr, best_inactive_dist);

    /* search for farthermost active site */
    best_active_dist = 0.0;
    for (i = 0; i < ct->sites; i++)
    {
        dvec_sub(ct->site[i].com, ct->coc, bond);
        if (dnorm(bond) > best_active_dist)
        {
            best_active_dist = dnorm(bond);
            farthest_active  = i;
        }
    }
    //printf("farthest site is %d (res %d) dist %f\n",farthest_active,ct->site[farthest_active].resnr, best_active_dist);



    /* substitute farthermost active with nearer inactive site */
    if (best_active_dist < best_inactive_dist)
    {
        printf("Did not adapt QM zone. Is already optimal.\n");
    }
    else
    {
        tot_occ = 0.0;
        counter = 0;
        for (i = 0; i < ct->sites; i++)
        {
            for (j = 0; j < ct->site[i].homos; j++)
            {
                if (i == farthest_active)
                {
                    tot_occ += ct->occupation[counter];
                }
                counter++;
            }
        }
        if (tot_occ * ct->sites > 0.1)//more than 10% of the average occupation of the other sites.
        {
            printf("Aborted replacement of residue %d with %d because of non-neglible occupation %f.\n", ct->site[farthest_active].resnr, ct->pool_site[nearest_inactive].resnr, tot_occ);
        }
        else
        {
            printf("Replacing residue %d with %d. Distance to COC reduced form %fnm to %fnm.\n", ct->site[farthest_active].resnr, ct->pool_site[nearest_inactive].resnr, best_active_dist/NM_TO_BOHR, best_inactive_dist/NM_TO_BOHR);
            counter = 0;
            for (i = 0; i < ct->sites; i++)
            {
                for (j = 0; j < ct->site[i].homos; j++)
                {
                    if (i == farthest_active)
                    {
                        ct->wf[counter] = 0.0; ct->wf[counter+ct->dim] = 0.0;
                    }
                    counter++;
                }
            }
            ct->pool_site[nearest_inactive].active = 1;
            for (i = 0; i < ct->pool_size; i++)
            {
                if (ct->site[farthest_active].resnr == ct->pool_site[i].resnr)
                {
                    ct->pool_site[i].active = 0;
                }
            }
            ct->site[farthest_active] = ct->pool_site[nearest_inactive];

            /* center all atoms around COC */
            calc_box_center(ecenterTRIC, state_box, box_center);
            for (i = 0; i < DIM; i++)
            {
                dx[i] = box_center[i] - (real) ct->coc[i]/NM_TO_BOHR;
            }

            for (i = 0; i < top_global->natoms; i++)
            {
                rvec_inc(gromacs_x[i], dx);
                rvec_inc(x_ct[i], dx);
            }
            //printf("coc is %f %f %f  nm.\n ", ct->coc[0]/NM_TO_BOHR,ct->coc[1]/NM_TO_BOHR,ct->coc[2]/NM_TO_BOHR);
            //printf("box_cneter is %f %f %f  nm.\n ", box_center[0],box_center[1], box_center[2]);
            for (i = 0; i < DIM; i++)
            {
                ct->coc[i]       += (double) dx[i] * NM_TO_BOHR;
                ct->coc_start[i] += (double) dx[i] * NM_TO_BOHR;
            }
            printf("displaced all atoms by %f %f %f  nm.\n ", dx[0], dx[1], dx[2]);

            /* rescale the wavefunction */
            norm = 0;
            for (i = 0; i < ct->dim; i++)
            {
                norm += SQR(ct->wf[i]) + SQR(ct->wf[i+ct->dim]);
            }
            //printf("norm has changed to %f due to adaption of the QM zone.\n", norm);
            for (i = 0; i < 2*ct->dim; i++)
            {
                ct->wf[i] /= sqrt(norm);
            }

            /* set the active atoms of the complex */
            counter = 0;
            for (i = 0; i < ct->sites; i++)
            {
                for (j = 0; j < ct->site[i].atoms; j++)
                {
                    ct->atom_cplx[counter]     = ct->site[i].atom[j];
                    ct->atomtype_cplx[counter] = ct->site[i].atomtype[j];
                    counter++;
                }
            }

            /* set the new external charges of the complex */
            if (ct->qmmm > 0)
            {
                for (i = 0; i < top_global->natoms; i++)
                {
                    ct->extcharge_cplx[i] = i;
                }
                for (i = 0; i < ct->sites; i++)
                {
                    k = find_intersection(top_global->natoms, ct->extcharge_cplx, ct->site[i].extcharge, ct->extcharge_cplx); // this should successively reduce the charges in ct->extcharge_cplx.
                    for (j = k+1; j <= top_global->natoms; j++)                                                               // k is index of highest common entry
                    {
                        ct->extcharge_cplx[j] = -1;
                    }
                }
                counter = 0;
                for (l = 0; l < ct->extcharges_cplx; l++)
                {
                    for (i = 0; i < ct->sites; i++)
                    {
                        for (j = 0; j < ct->site[i].bonds; j++)
                        {
                            for (k = 0; k < ct->site[i].addchrs[j]; k++)
                            {
                                if (ct->extcharge_cplx[l] == ct->site[i].extcharge[ ct->site[i].modif_extcharge[j][k] ]) // if one of the extcharges of the complex is the same atom that was modified in the monomer calculation, then also modify it in the complex calculation.
                                {
                                    ct->modif_extcharge_cplx[counter] = l;
                                    counter++;
                                }
                            }
                        }
                    }
                }
            }
            adapted = 1;
        }
    }

    return adapted;
}


void check_and_invert_orbital_phase(dftb_phase1_t *dftb1, charge_transfer_t *ct, t_state *state_global, t_mdatoms *mdatoms)
{
    int    iatom, i, nn, ii, j, k, l, m, iao, jao;
    double overl, norm1, norm2, proj;
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

       for (iatom = mdatoms->start; iatom < mdatoms->start+mdatoms->homenr; iatom++)
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

void get_spectrum(charge_transfer_t *ct, dftb_t *dftb) //, FILE *f_ct_spec, FILE *f_ct_spec_evec){
{
    int           i, j;
    dftb_phase2_t dftb2;


    dftb2 = dftb->phase2;

    // diagonalize FO-Hamiltonian
    for (i = 0; i < ct->dim; i++)
    {
        ct->evec_spec[i + i * ct->dim] = ct->hamiltonian[i][i];
        for (j = 0; j < ct->dim; j++)
        {
            if (i != j)
            {
                ct->evec_spec[i + j * ct->dim] = ct->hamiltonian[i][j];
            }
        }
    }
    dsyev(ct->dim, ct->evec_spec, ct->ev_spec, ct->work_spec, 3*ct->dim);



/*need larger work arrays
    // diagonalize AO-Hamiltonian
    for (i=0; i<dftb2.norb; i++) {
      ct->evec_spec[i + i * dftb2.norb] = dftb2.hamil[i][i];
      for (j=0; j<dftb2.norb; j++) {
        if (i!=j) ct->evec_spec[i + j * dftb2.norb] = dftb2.hamil[i][j];
      }
    }
    dsyev(dftb2.norb, ct->evec_spec, ct->ev_spec, ct->work_spec, 3*dftb2.norb);

   fprintf(f_ct_spec, "energy eigen values of AO Hamiltonian:");
   for (i=0; i<dftb2.norb; i++)
     fprintf(f_ct_spec, " %10.6f \n", ct->ev_spec[i] * HARTREE_TO_EV);
   fprintf(f_ct_spec, "\n");
 */


    return;
}

void write_out_MOs(int step, rvec x_ct, t_atoms *ct_atoms, dftb_t *dftb, charge_transfer_t *ct)
{
    double amp, phase;
    int    i, j, lj, mj, jofn, jofn_old, mo = 1, counter, k;
    char   c, filename[20];
    FILE  *f_ct_step_pdb = NULL;

    sprintf(filename, "%d", step);
    strcat(filename, ".pdb" );
    //f_ct_step_pdb=fopen(filename,"w");
/*
   phase=1;
   for (j=0; j < dftb->phase2.nn; j++) {
    switch (dftb->phase2.izp[j]) {
      case 0: c = 'C'; break;
      case 1: c = 'H'; break;
      case 2: c = 'N'; break;
      case 3: c = 'O'; break;
    }
    amp = 0.0;
    if (dftb->lmax[dftb->phase2.izp[j]] == 2) {
      lj=1; //only p-orbitals so far
      jofn = dftb->phase2.ind[j] + lj*lj;
      for (mj=0; mj<=2*lj; mj++) {
        amp +=  SQR(dftb->orthogo.evec_ao[jofn+mj][mo]);
      }
      amp= sqrt(amp)*100; //for better reading
      if(j!=0){
        phase = 0;
        for (mj=0; mj<=2*lj; mj++)
          phase +=  dftb->orthogo.evec_ao[jofn+mj][mo]*dftb->orthogo.evec_ao[jofn_old+mj][mo];
        phase=(phase > 0) ? 1 : ((phase < 0) ? -1 : 0);
      }
      jofn_old = jofn;
      fprintf(f_ct_step_pdb,"%6s%5i   %c %s %s%4s    %8.3f%8.3f%8.3f%6s%6.2f\n" ,"ATOM  ",j+1, c ,"res","I","  1 ",dftb->phase2.x[j][0]/NM_TO_BOHR*10,dftb->phase2.x[j][1]/NM_TO_BOHR*10,dftb->phase2.x[j][2]/NM_TO_BOHR*10,"  1.00",phase*amp);//pdb format
    }
    else{
      fprintf(f_ct_step_pdb,"%6s%5i   %c %s %s%4s    %8.3f%8.3f%8.3f%6s%6s\n" \
      ,"ATOM  ",j+1,c,"res","I","  1 ",dftb->phase2.x[j][0]/NM_TO_BOHR*10,dftb->phase2.x[j][1]/NM_TO_BOHR*10,dftb->phase2.x[j][2]/NM_TO_BOHR*10,"  1.00","  0.00");//pdb format
    }
   }
 */
// write interesting data in pbd file as bfactor
    counter = 0;
    for (i = 0; i < ct->sites; i++)
        for (j = 0; j < ct->site[i].atoms; j++)
        {
            ct_atoms->pdbinfo[ct->site[i].atom[j]].bfac = 0.0;
            for (k = 0; k < ct->site[i].homos; k++)
            {
                //ct_atoms->pdbinfo[ct->site[i].atom[j]].bfac =(real) (dftb->phase2.qmat[counter] - dftb->qzero2[dftb->phase2.izp[counter]]);
                ct_atoms->pdbinfo[ct->site[i].atom[j]].bfac += (real) 100*(ct->occupation[ct->indFO[i]+k]);
                counter++;
            }
        }

        write_sto_conf(filename, "written by charge transfer code", ct_atoms, x_ct, NULL, 0, NULL);




        //fclose(f_ct_step_pdb);

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

long exp_imag_matrix(double **in, twodoubles **out, double dt, long n, ct_per_orthogo_t *arrays)
{
    extern void dsyevr_(char *, char *, char *, long *, double *, long *,
                        double *, double *, long *, long *, double *, long *, double *,
                        double *, long *, double *, double *, long *, double *, long *, long *);
    double abstol = 1.0e-8;
    char   itype = 'U';
    char   jobz = 'V';
    char   range = 'A';
    char   uplo = 'U';
    long   i, j, k, m, info, eval_no;
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


            return info;
}

void print_time_difference(char *s, struct timespec start, struct timespec end)
{
    int       sec, nsec;
    long long value = 0ll;

    value = 1000000000ll * ((long long) end.tv_sec - (long long) start.tv_sec) + (long long) (end.tv_nsec - start.tv_nsec);
    printf("%s %12lld\n", s, value);

    return;
}

int find_intersection(int size, int array1[], int array2[], int intersection_array[])
{
    int i = 0, j = 0, k = 0;
    while ((i < size) && (j < size))
    {
        if (array1[i] < array2[j])
        {
            i++;
        }
        else if (array1[i] > array2[j])
        {
            j++;
        }
        else
        {
            intersection_array[k] = array1[i];
            i++;
            j++;
            k++;
        }
    }
    return(k);
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

void assignprob(charge_transfer_t *ct, int *ind)
{

    int     i, j, k, step, nstate, done, path_row_0, path_col_0, path_count;
    int     nsmax = 2;
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

     int i, j, k;


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

int star_in_row(int nstate, int irow, double **star)
{
    int j, tmp;

    tmp = 0;
    for (j = 0; j < nstate; j++)
        if (star[irow][j] == 1)
            tmp = 1;

    return tmp;
}

void find_star_in_row(int nstate, int irow, double **star, int *icol)
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

// KERNEL RIDGE
void get_ml_hamiltonian(charge_transfer_t *ct, t_state *state_global, t_mdatoms *mdatoms)
{
    long    iatom, isite, jsite, count, index, i, j, k, natoms, ncm;
    double  *cm, *atomcharge, **pair_distance_matrix, **coordinates;
    double   gauss_kernel, norm2;

//  store the coordinates of atoms in QM zone
    for (iatom = mdatoms->start; iatom < mdatoms->start+mdatoms->homenr; iatom++)
    {
        for (i = 0; i < ct->sites; i++)
        {
            for (j = 0; j < ct->site[i].atoms; j++)
            {
                if (iatom == ct->site[i].atom[j])
                {
                   for (k = 0; k < DIM; k++)
                       ct->coords[i][j][k] = state_global->x[iatom][k]*10.0;
                }
            }
        }
    }

// calculate hamiltonian by fitting to a linear function by gaussian kernel regression
    natoms = 2*ct->site[0].atoms;
    ncm = natoms*(natoms+1)/2;

    snew(cm, ncm);
    snew(atomcharge, natoms);

    snew(pair_distance_matrix, natoms);
    for (i = 0; i < natoms; i++)
        snew(pair_distance_matrix[i], natoms);

    snew(coordinates, natoms);
    for (i = 0; i < natoms; i++)
        snew(coordinates[i], DIM);

// NOTE ASSUME ALL SITES ARE THE SAME!!
    for (i = 0; i < ct->site[0].atoms; i++)
    {
        switch (ct->site[0].atomtype[i])
        {
          case 0: atomcharge[i] = 6.0; break; // C
          case 1: atomcharge[i] = 1.0; break; // H
          case 2: atomcharge[i] = 7.0; break; // N
          case 3: atomcharge[i] = 8.0; break; // O
        }
        atomcharge[ct->site[0].atoms+i] = atomcharge[i];
    }

// compute the diagonal terms of pair distance matrix
    for (i = 0; i < natoms; i++)
        pair_distance_matrix[i][i] = 0.5 * pow(atomcharge[i], 2.4);

    for (isite = 0; isite < ct->sites; isite++)
        for (jsite = 0; jsite < ct->sites; jsite++)
            ct->hamiltonian_ml[isite][jsite] = 0.0;

// diagonal terms predictions
    natoms = ct->site[0].atoms;
    ncm = natoms*(natoms+1)/2;

    for (isite = 0; isite < ct->sites; isite++)
    {
        count = 0;

        for (i = 0; i < natoms; i++)
        {
            for (j = i + 1; j < natoms; j++)
            {
                norm2 = 0.0;
                for (k = 0; k < DIM; k++)
                    norm2 += SQR(ct->coords[isite][j][k] - ct->coords[isite][i][k]);

// off-diagonal terms pair distance matrix
                pair_distance_matrix[i][j] = atomcharge[i] * atomcharge[j] / sqrt(norm2);
                pair_distance_matrix[j][i] = pair_distance_matrix[i][j];
            }

// coulomb matrix
            for (j = 0; j <= i; j++)
                cm[count+j] = pair_distance_matrix[i][j];

            count += i+1;
        }

// gaussian kernel function
        for (i = 0; i < ct->ntrain[0]; i++)
        {
            norm2 = 0.0;

            for (j = 0; j < ncm; j++)
            {
                norm2 += SQR(cm[j] - ct->cm_train[0][i][j]);
            }

            gauss_kernel = exp(-norm2/2.0/SQR(ct->sigma[0]));

            ct->hamiltonian_ml[isite][isite] += ct->alpha[0][i]*gauss_kernel;
        }

// Make the predictions
        ct->hamiltonian_ml[isite][isite] /= HARTREE_TO_EV;

// implicit reorganization
        if (ct->jobtype == cteNOMOVEMENT)
        {
        ct->hamiltonian_ml[isite][isite] -= ct->lambda*ct->occupation[isite];
    }
        else
        {
            ct->hamiltonian_ml[isite][isite] -= ct->lambda*(SQR(ct->tfs_diab[isite]) + SQR(ct->tfs_diab[isite+ct->dim]));
        }

    }

// off-diagnoal terms
    natoms = 2*ct->site[0].atoms;
    ncm = natoms*(natoms+1)/2;

    for (isite = 0; isite < ct->sites - 1; isite++)
    {
        for (jsite = isite + 1; jsite < ct->sites; jsite++)
        {
            for (i = 0 ; i < ct->site[0].atoms; i++)
            {
                for (k = 0 ; k < DIM; k++)
                {
                    coordinates[i][k] = ct->coords[isite][i][k];
                    coordinates[i+ct->site[0].atoms][k] = ct->coords[jsite][i][k];

                }

            }

            count = 0;
            for (i = 0; i < natoms; i++)
            {
                for (j = i + 1; j < natoms; j++)
                {
                    norm2 = 0.0;

                    for (k = 0; k < DIM; k++)
                        norm2 += SQR(coordinates[j][k] - coordinates[i][k]);

// off-diagonal terms pair distance matrix
                    pair_distance_matrix[i][j] = atomcharge[i] * atomcharge[j] / sqrt(norm2);
                    pair_distance_matrix[j][i] = pair_distance_matrix[i][j];
                }

// coulomb matrix
                for (j = 0; j <= i; j++)
                    cm[count+j] = pair_distance_matrix[i][j];

                count += i+1;
            }

// gaussian kernel function
            for (i = 0; i < ct->ntrain[1]; i++)
            {
                norm2 = 0.0;

                for (j = 0; j < ncm; j++)
                    norm2 += SQR(cm[j] - ct->cm_train[1][i][j]);

                gauss_kernel = exp(-norm2/2.0/SQR(ct->sigma[1]));

                ct->hamiltonian_ml[isite][jsite] += ct->alpha[1][i]*gauss_kernel;
            }

// Make the predictions
            ct->hamiltonian_ml[isite][jsite] /= HARTREE_TO_EV;
            ct->hamiltonian_ml[jsite][isite] = ct->hamiltonian_ml[isite][jsite];
        }
    }

    for (i = 0; i < natoms; i++)
    {
        sfree(pair_distance_matrix[i]);
        sfree(coordinates[i]);
    }

    sfree(cm);
    sfree(atomcharge);
    sfree(pair_distance_matrix);
}


/* MACHINE LEARNING WITH TF NNs -- MK 2020 */

/*TF tensors need this*/
void NoOpDeallocator(void* data, size_t a, void* b) {}

void initialize_nn_model(tf_model *model, const char *path, const char *hyperpath) {
  /* Initialize saved model from folder path into struct. */
  model->graph = TF_NewGraph();
  model->status = TF_NewStatus();
  TF_SessionOptions* SessionOpts = TF_NewSessionOptions();
  TF_Buffer* RunOpts = NULL;
  const char* tags = "serve";
  model->session = TF_LoadSessionFromSavedModel(SessionOpts, RunOpts, path, &tags, 1, model->graph, NULL, model->status);
  printf("session loaded successfully!\n");
  TF_DeleteSessionOptions(SessionOpts);
  TF_DeleteBuffer(RunOpts);

  /* Preparation of Input and output tensors */
  model->num_inputs = 1; //one input tensor
  model->input = malloc(sizeof(TF_Output) * model->num_inputs);
  // serving_default_input_1 is the default name, get from saved_model_cli
  TF_Output t0 = {TF_GraphOperationByName(model->graph, "serving_default_input_1"), 0};
  model->input[0] = t0;

  // we have two output tensors - the coupling and the force
  model->num_outputs = 2;
  model->output = malloc(sizeof(TF_Output) * model->num_outputs);
  // StatefulPartitionedCall is the default name, get from saved_model_cli
  TF_Output t2 = {TF_GraphOperationByName(model->graph, "StatefulPartitionedCall"), 0};
  model->output[0] = t2;
  TF_Output t3 = {TF_GraphOperationByName(model->graph, "StatefulPartitionedCall"), 1};
  model->output[1] = t3;

  model->input_val = (TF_Tensor**)malloc(sizeof(TF_Tensor*)*model->num_inputs);
  model->output_val = (TF_Tensor**)malloc(sizeof(TF_Tensor*)*model->num_outputs);

  FILE *hyperfile = NULL;
  hyperfile = fopen(hyperpath, "r");
  float hypers[5];
  printf("Reading hyperparameter file\n");
  for (int i=0; i<5; i++){
      fscanf(hyperfile, "%f \n", &hypers[i]);
      //printf("%f\n", hypers[i]);
      }
  model->x_mean=hypers[0];
  model->x_std=hypers[1];
  model->y_mean=hypers[2];
  model->y_std=hypers[3];
  model->grad_std=hypers[4];
  printf("xmean %f , xstd %f , ymean %f , ystd %f , grad_std %f \n", model->x_mean, model->x_std, model->y_mean, model->y_std, model->grad_std);
  }

void check_nn_status(tf_model *model){
  // make sure everything worked
  if(TF_GetCode(model->status) == TF_OK)
  {
    printf("TF status OK\n");
  }
  else
  {
    printf("%s", TF_Message(model->status));
  }
}



void input_data_to_model(tf_model *model, float data[][3], int natom){
  /* pass data to model */

  int64_t dims[] = {1,natom,DIM};
  int ndims = 3; //length of above array
  int ndata = sizeof(float)*natom*DIM; // This is tricky, it's number of bytes not number of elements

  /* To get a prediction now we allocate the input tensor with the input data */;

  TF_Tensor* int_tensor = TF_NewTensor(TF_FLOAT, dims, ndims, data, ndata, &NoOpDeallocator, 0);
  // make sure everything worked
  // if (int_tensor != NULL)
  //  {
  //     printf("TF_NewTensor is OK\n");
  //  }
  // else
  //  {
  //   printf("ERROR: Failed TF_NewTensor\n");
  //  }

  // and pass the tensor to the model
  model->input_val[0] = int_tensor;

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


void get_nn_hamiltonian(charge_transfer_t *ct, t_state *state_global, t_mdatoms *mdatoms, dftb_t *dftb)
{
/* NN Model trained with AIMAT group's code for learning energies & forces
 * Note: Model takes coordinates in Angstrom and returns energies in Hartree
 * and forces in Hartree/nm.
 * NOTE: We abuse the dftb struct's fields to pass the gradients to the rest of the code.
 * This avoids having to go through the entire codebase to implement ML force support.
 */
  long    iatom, isite, jsite, count, index, i, j, k, natoms, ncm;
  int atom_qm_idx, idx_i, idx_j;
  long idx_lookup[ct->sites][ct->site[0].atoms];
  int oth_lookup[2*ct->site[0].atoms];
  float coords[ct->sites][ct->site[0].atoms][DIM];
  double tmp_frc, wf_factor;

  //  store the coordinates of atoms in QM zone
  //  and create lookup for indices needed for forces
  for (iatom = mdatoms->start; iatom < mdatoms->start+mdatoms->homenr; iatom++)
  {
    atom_qm_idx = 0;
    for (i = 0; i < ct->sites; i++)
    {
      for (j = 0; j < ct->site[i].atoms; j++)
      {
        if (iatom == ct->site[i].atom[j])
        {
          for (k = 0; k < DIM; k++)
          {
            coords[i][j][k] = (float)state_global->x[iatom][k]*10.0;
          }
          idx_lookup[i][j] = atom_qm_idx;
        }
        atom_qm_idx++;
      }
    }
  }

  // // faster version, but might be buggy for forces:
  // long iatom_global, site, iatom_local;
  // int m;
  //
  // atom_qm_idx = 0;
  // for (site = 0; site < ct->sites; site++)
  // {
  //   for (iatom_local = 0; iatom_local < ct->site[i].atoms; iatom_local++)
  //   {
  //     for (m = 0; m < DIM; m++)
  //     {
  //       iatom_global = ct->site[site].atom[iatom_local];
  //       coords[site][iatom_local][m] = (float)state_global->x[iatom_global][m]*10.0;
  //     }
  //     // construct lookup array for forces later
  //     idx_lookup[site][iatom_local] = atom_qm_idx;
  //     atom_qm_idx++;
  //   }
  // }

  natoms = 2*ct->site[0].atoms;

  //coordinate array for pairs
  float coordinates[natoms][DIM];

  // zero out all terms of the hamiltonian - why? MK
  for (isite = 0; isite < ct->sites; isite++)
      for (jsite = 0; jsite < ct->sites; jsite++)
          ct->hamiltonian_ml[isite][jsite] = 0.0;

  // off-diagonal terms first, because of scaling... MK
  natoms = 2*ct->site[0].atoms;
  // allocate array for forces (reused for every prediction)
  float** offdiag_forces;
  offdiag_forces = malloc(natoms*sizeof(float*));
  tf_model *coupling_model = ct->coupling_model;

  int pairs_done = 0;
  double offdiag_start =  (double)clock()/CLOCKS_PER_SEC;
  printf("calculating offdiag values start at %f\n", offdiag_start);

  double nl_start =  (double)clock()/CLOCKS_PER_SEC;

  if (ct->neighbors_only == 1 || (ct->neighbors_only==-1 && ct->first_step==1) || (ct->neighbors_only>1 && (ct->step % ct->neighbors_only)==0)){
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
        if (dist < 20.0*BOHR2NM){
          dftb->nl[counter1][counter2]=dftb->nl[counter2][counter1]=1;
        }else{
          dftb->nl[counter1][counter2]=dftb->nl[counter2][counter1]=0;
        }
      }
    }
  }

  double nl_end =  (double)clock()/CLOCKS_PER_SEC;
  printf("calculating neighborlist took %f\n", nl_end - nl_start);

  for (isite = 0; isite < ct->sites - 1; isite++)
  {
    for (jsite = isite + 1; jsite < ct->sites; jsite++)
    {
      if (ct->neighbors_only!=0 && (!dftb->nl[isite][jsite]))
      {
        continue;
      }

      pairs_done++;

      for (i = 0 ; i < ct->site[0].atoms; i++)
      {
        for (k = 0 ; k < DIM; k++)
        {
          coordinates[i][k] = (coords[isite][i][k] - coupling_model->x_mean)/coupling_model->x_std;
          coordinates[i+ct->site[0].atoms][k] = (coords[jsite][i][k] - coupling_model->x_mean)/coupling_model->x_std;
         // coordinates[i+ct->site[0].atoms][k] = coords[jsite][i][k]; no scaling
        }
        oth_lookup[i] = idx_lookup[isite][i];
        oth_lookup[i+ct->site[0].atoms] = idx_lookup[jsite][i];
      }

      input_data_to_model(coupling_model, coordinates, natoms);
      // Get the predictions
      TF_SessionRun(coupling_model->session, NULL, coupling_model->input, coupling_model->input_val, coupling_model->num_inputs, coupling_model->output, coupling_model->output_val, coupling_model->num_outputs, NULL, 0,NULL , coupling_model->status);
      //check_nn_status(coupling_model);

      float* buffer = (float*)TF_TensorData(coupling_model->output_val[0]);
      ct->hamiltonian_ml[isite][jsite] = buffer[0] * coupling_model->y_std + coupling_model->y_mean;
//          printf("coupling model predicted %f, that's %f Hartree\n", buffer[0], ct->hamiltonian_ml[isite][jsite]);
      ct->hamiltonian_ml[jsite][isite] = ct->hamiltonian_ml[isite][jsite];

      // we get the forces as 1D array (with natoms*3 elements)
      float* offsets_force = TF_TensorData(coupling_model->output_val[1]);
      // we reverse the Scaling

      for (int i=0; i<DIM*natoms; i++)
      {
        offsets_force[i] *= coupling_model->grad_std;
      }
      // we can reshape into 2D array in-place
      for(int i=0; i<natoms; i++)
      {
        offdiag_forces[i] = offsets_force + i*DIM;
      }

      // now we pass that to dftb->phase2.grad in the order we got earlier:
      for(int l=0; l < natoms; l++)
      {
        idx_i = oth_lookup[l];
        //printf("indices of atoms %d %d\n", idx_i, idx_j);
        for(int m=0; m<DIM; m++)
        {
          //MK tag forces: set wf_factor = 1 for pure gradients w/o wf
          if (ct->pure_forces){
            wf_factor = 1.0;
          }else if(dftb->do_ncv){
            wf_factor = ct->tfs_vector[ct->ncv_isurf][isite]*ct->tfs_vector[ct->ncv_jsurf][jsite];
          }else{
            wf_factor = ct->wf[isite]*ct->wf[jsite] + ct->wf[isite+ct->dim]*ct->wf[jsite+ct->dim];
          }
          tmp_frc = wf_factor * offdiag_forces[l][m];
          dftb->phase2.grad[idx_i][m] += ct->is_hole_transfer ? tmp_frc : -1.0*tmp_frc;
        }
//        printf("frc for atom %d:  %f %f %f\n", idx_i, dftb->phase2.grad[idx_i][0], dftb->phase2.grad[idx_i][1], dftb->phase2.grad[idx_i][2]);
      }
    }
  }

  double offdiag_end = (double)clock()/CLOCKS_PER_SEC;
  double diff = offdiag_end - nl_end;

  printf("offdiagonal predictions took %f, for %d pairs %f per pair\n", diff, pairs_done, diff/(double)pairs_done);

  free(offdiag_forces);

  // set up site prediction
  natoms = ct->site[0].atoms;
  // allocate array for forces (reused for every prediction)
  float** diag_forces;
  diag_forces = malloc(natoms*sizeof(float*));

  tf_model *site_model = ct->site_model;

  //scale according to nn hyperparams
  for (i=0; i < ct->sites; i++)
    {
      for (j=0; j < ct->site[i].atoms; j++)
        {
          for (k=0; k<DIM; k++)
          {
            coords[i][j][k] = (coords[i][j][k] - site_model->x_mean)/site_model->x_std;
//              coords[i][j][k] = coords[i][j][k]; //no scaling
          }
        }
    }

  double diag_start = (double)clock()/CLOCKS_PER_SEC;
  printf("calculating diag values start at %f\n", diag_start);
  for (int isite = 0; isite < ct->sites; isite++)
  {
      input_data_to_model(site_model, coords[isite], natoms);
      // Get the predictions
      TF_SessionRun(site_model->session, NULL, site_model->input, site_model->input_val, site_model->num_inputs, site_model->output, site_model->output_val, site_model->num_outputs, NULL, 0,NULL , site_model->status);
      //check_nn_status(site_model);
      float* buffer = (float*)TF_TensorData(site_model->output_val[0]);
      ct->hamiltonian_ml[isite][isite] = buffer[0] * site_model->y_std + site_model->y_mean;
//      printf("site model predicted %f, that's %f Hartree\n", buffer[0], ct->hamiltonian_ml[isite][isite]);

      // we get the forces as 1D array (with natoms*3 elements)
      float* offsets_force = TF_TensorData(site_model->output_val[1]);
      // we reverse the Scaling
      for (int i=0; i<DIM*natoms; i++){
        offsets_force[i] = offsets_force[i] * site_model->grad_std;
      }
      // we can reshape into 2D array in-place
      for(int i=0; i<natoms; i++)
      {
        diag_forces[i] = offsets_force + i*DIM;
      }
      for(int i=0; i<natoms; i++)
      {
        for(int m=0; m<DIM; m++)
        {
          //MK tag forces: set wf_factor = 1 for pure gradients w/o wf
          if (ct->pure_forces){
            wf_factor = 1.0;
          }else if (dftb->do_ncv){
            wf_factor = ct->tfs_vector[ct->ncv_isurf][ct->indFO[isite]]*ct->tfs_vector[ct->ncv_jsurf][ct->indFO[isite]];
          }else{
            wf_factor = (ct->wf[ct->indFO[isite]]*ct->wf[ct->indFO[isite]] + ct->wf[ct->indFO[isite]+ct->dim]*ct->wf[ct->indFO[isite]+ct->dim]);
          }
          dftb->phase1[isite].grad[i][m] = wf_factor*diag_forces[i][m];
        }
//        printf("diag_frc %f %f %f\n",  dftb->phase1[isite].grad[i][0],  dftb->phase1[isite].grad[i][1], dftb->phase1[isite].grad[i][2]);

      }
      if (ct->jobtype != cteNOMOVEMENT && ct->lambda != 0.0)
      {
      // implicit reorganization
      ct->occupation[isite] = SQR(ct->tfs_diab[isite]) + SQR(ct->tfs_diab[isite+ct->dim]);
      ct->hamiltonian_ml[isite][isite] -= ct->lambda*ct->occupation[isite];
      }
  }

  double diag_end = (double)clock()/CLOCKS_PER_SEC;
  printf("diagonal predictions took %f, %f per site\n", diag_end - diag_start,  (diag_end-diag_start)/(double)ct->sites);

  free(diag_forces);

  printf("hamiltonian calculation done\n");
}
