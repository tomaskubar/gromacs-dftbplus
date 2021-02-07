#include "gromacs/transfer/transfer.h"


void ct_init_diis(charge_transfer_t *ct, ct_diis_t *diis)
{
    diis->n_elem = ct->dim;
 // diis->prev_q_input = (double**) malloc(DIIS_MAX_PREV_VECTORS * sizeof(double *));
 // diis->prev_q_input[0] = (double*) malloc(diis->n_elem * DIIS_MAX_PREV_VECTORS);
    snew(diis->prev_q_input, DIIS_MAX_PREV_VECTORS);
    snew(diis->prev_q_input[0], diis->n_elem * DIIS_MAX_PREV_VECTORS);
    for (int i = 1; i < DIIS_MAX_PREV_VECTORS; i++)
    {
        diis->prev_q_input[i] = diis->prev_q_input[0] + i * diis->n_elem;
    }
 // diis->prev_q_diff = (double**) malloc(DIIS_MAX_PREV_VECTORS * sizeof(double *));
 // diis->prev_q_diff[0] = (double*) malloc(diis->n_elem * DIIS_MAX_PREV_VECTORS);
    snew(diis->prev_q_diff, DIIS_MAX_PREV_VECTORS);
    snew(diis->prev_q_diff[0], diis->n_elem * DIIS_MAX_PREV_VECTORS);
    for (int i = 1; i < DIIS_MAX_PREV_VECTORS; i++)
    {
        diis->prev_q_diff[i] = diis->prev_q_diff[0] + i * diis->n_elem;
    }
 // diis->aa = (double *) malloc(SQR(DIIS_MAX_PREV_VECTORS + 1));
 // diis->bb = (double *) malloc(DIIS_MAX_PREV_VECTORS + 1);
 // diis->q_inp_result = (double *) malloc(ct->sites * sizeof(double));
 // diis->q_diff = (double *) malloc(ct->sites * sizeof(double));
    snew(diis->aa, SQR(DIIS_MAX_PREV_VECTORS + 1));
    snew(diis->bb, DIIS_MAX_PREV_VECTORS + 1);
 // snew(diis->q_inp_result, ct->sites);
 // snew(diis->q_diff, ct->sites);
    snew(diis->fermi_coef, ct->dim);
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
}


void ct_get_index(int isize[], int *index[])
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
}

int ct_atom_in_group(int atom, int *list, int list_size)
{
    for (int i = 0; i < list_size; i++)
    {
        if (atom == list[i])
        {
            return 1;
        }
    }

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
    for (int line = 0; line < lines; line++)
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
    //char input[MAXLINES][2][MAXWIDTH];

    FILE *f = fopen(file_name, "r");
    if (f == NULL)
    {
        printf("%s not accessible, exiting!\n", file_name);
        exit(-1);
    }
    printf("Searching in file %s for:\n", file_name);
    for (int i = 0; i < nkeys; i++)
    {
        printf("%20s:  %s\n", possiblekeys[i][0], possiblekeys[i][1]);
    }

    /* get length of input file */
    int lines = 0;
    while (!feof(f))
    {
        int ch = fgetc(f);
        if (ch == '\n')
        {
            lines++;
        }
    }
    //rewind(f);
    fclose(f);
    f = fopen(file_name, "r");

    /* read input file */
    for (int line = 0; line < lines; line++)
    {
        //ch=fscanf(f, "%[;a-zA-Z0-9/] = %[-a-zA-Z0-9/{}. ]\n", input[line][0], input[line][1]);
        int ch = fscanf(f, " %[^= \t\r] = %[^=\n\t\r] \n", input[line][0], input[line][1]);
        if (ch != 2)
        {
            printf("READING ERROR in line %d.\n Use format 'keyword = value'.\n ", line+1); exit(0);
        }

        //remove trailing whitespace
        char *s = input[line][1] + strlen(input[line][1]);
        while (--s >= input[line][1])
        {
            if (!isspace(*s))
            {
                break;
            }
            *s = 0;
        }

        //convert all to lowercase for better handling (besides file and atom names)
        int len = strlen(input[line][0]);
        for (int i = 0; i < len; i++)
        {
            input[line][0][i] = tolower((unsigned char)input[line][0][i]);
        }
        if ( (strcmp(input[line][0], "slkopath") && strcmp(input[line][0], "specfiles") && strncmp(input[line][0], "name", 4) ) )
        {
            len = strlen(input[line][1]);
            for (int i = 0; i < len; i++)
            {
                input[line][1][i] = tolower((unsigned char)input[line][1][i]);
            }
        }
        // check input file
        int j = 1;
        for (int i = 0; i < nkeys; i++)
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
    for (int line = 0; line < lines; line++)
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

    int i = 0;
    char *ptr = strtok(value, " ");
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

    int i = 0;
    char *ptr = strtok(value, " ");
    while (ptr != NULL)
    {
        target[i] = atoi(ptr);
        ptr       = strtok(NULL, " ");
        i++;
    }
    if (i > n)
    {
        printf("TOO MANY ARGUMENTS FOR %s\n", name);
        exit(-1);
    }
    else if (i < n)
    {
        printf("TOO FEW ARGUMENTS FOR %s\n", name);
        exit(-1);
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

    int i = 0;
    char *ptr = strtok(value, " ");
    while (ptr != NULL)
    {
        strcpy(target[i], ptr);
        ptr = strtok(NULL, " ");
        i++;
    }
    if (i > n)
    {
        printf("TOO MANY ARGUMENTS FOR %s\n", name);
        exit(-1);
    }
    else if (i < n)
    {
        printf("TOO FEW ARGUMENTS FOR %s\n", name);
        exit(-1);
    }
    else
    {
        return 0;
    }
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


/***************************************
 * INITIALIZE THE CHARGE TRANSFER CODE *
 ***************************************/

#if GMX_MPI
void init_charge_transfer(t_atoms           *atoms,
                          const gmx_mtop_t  *top_global,
                          t_mdatoms         *mdatoms,
                          charge_transfer_t *ct,
                          t_state           *state,
                          int                ct_mpi_rank)
#else
void init_charge_transfer(t_atoms           *atoms,
                          const gmx_mtop_t  *top_global,
                          t_mdatoms         *mdatoms,
                          charge_transfer_t *ct,
                          t_state           *state)
#endif
{
/* Initializes the charge transfer code. Reads in CT files, sets up the stuff that is fixed variables, allocates memory */

// PARAMETERS:
// atoms       = (in) variable of GROMACS data-structure-type, which we generated outside of this function. Contains info like atomnames and residuenumbers, that we need to distinguish QM from MM atoms.
// top_global  = (in) GROMACS data structure, from which we only need total number of MM atoms
// mdatoms     = (in) variable of GROMACS data-structure-type (we need atomic masses)
// ct          = (out) main data structure with information about the charge transfer calculation
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

    int        i, j, k, l, m, n, counter, *counter_array, counter_cplx, QMLAcounter, MMLAcounter, QMCAcounter, QMCApool[QMCASIZE], QMCApoolcounter,
               environment, environment_cplx, modif_cplx, counter_modif_cplx = 0, mm_list_size, *mm_list;
    double     sum, bond_length, bond_length_best, X[3], Y[3], magnitude, mass;
    dvec       bond;
    ct_site_t *site, s;

    // variables to store input files(charge-transfer.dat, and *.spec):
    //   several lines, 2 columns (before and after '='), and string of several chars
    char       input1[MAXLINES][2][MAXWIDTH], input2[MAXLINES][2][MAXWIDTH];
    char       dummy[MAXELEMENTS][ELEMENTWIDTH];
    int        lines1, lines2;
    // names of files for machine learning, W. X + M. K
    char       ml_coupling_model_path[MAX_PATH_LENGTH], ml_site_energy_model_path[MAX_PATH_LENGTH];
    char       ml_site_energy_para_path[MAX_PATH_LENGTH], ml_coupling_para_path[MAX_PATH_LENGTH];
 // char       ml_coupling_train_path[MAX_PATH_LENGTH], ml_site_energy_train_path[MAX_PATH_LENGTH];
 // char       ml_site_energy_para_path[MAX_PATH_LENGTH], ml_coupling_para_path[MAX_PATH_LENGTH];

    ct->first_step = 1;
#if GMX_MPI
    printf("Initializing charge transfer at rank %d\n", ct_mpi_rank);
#else
    printf("Initializing charge transfer\n");
#endif

    /* Get input from file charge-transfer.dat */
    lines1 = read_file("charge-transfer.dat", input1, possiblekeys1, sizeof(possiblekeys1)/sizeof(possiblekeys1[0]));

    /* set default options */

    // either sensible values or "not defined"=-1 for values that have to be specified in the input file and will be checked later.
 // ct->n_avg_ham = 1;           // average hamilton over n_avg_ham steps to assimilate fast non-classical vibrations
 // ct->esp_scaling_factor = 1.; // scaling of the electrostatic potential of the environment
 // ct->opt_QMzone=0;
 // ct->neg_imag_pot = 0;        // negative imaginary potential to drain the charge at some sites

    /* evaluate input */

    /* read in stuff that is always needed */
    /* read seed for random number*/
    if (searchkey(lines1, input1, "seed", value, 0))
    {
       ct->rnseed = atoi(value);
    }

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
         // for (i=0; i<mm_list_size; i++) PRINTF(" %d", mm_list[i]); PRINTF("\n");
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
        ct->esp_scaling_factor = 1.0; // default
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

    ct->decoherence = 0; // default is propagator without decoherence
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
    else if (strcmp(value, "bod") == 0) // JJK
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
        ct->jobtype = cteDLZSH; // lzsh is new one
        PRINTF("Landau-Zener type surface hopping method based on original diabatic LZ formula \n");
    }
    else if (strcmp(value, "alzsh") == 0)
    {
        ct->jobtype = cteALZSH; // lzsh is new one
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
        ct->jobtype = cteTULLYLOC; // tfl is new one
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
        ct->n_avg_ham = 1; // default
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
            ct->do_projection = 0; // default
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
            ct->sic = 0.0; // default
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
             // PRINTF("Did not understand relaxation model.\n");
             // exit(-1);
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
            ct->do_lambda_i = 0; // default
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
            ct->delta_q_mode = 0; // default
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
                split_string_into_int(value, ct->sitetype[i].bonds, ct->sitetype[i].addchrs, "naddchr");
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
                    split_string_into_string(value, l, dummy, "nameaddchr");
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
     // PRINTF("Sitetype %d found to have %d connections to other fragments \n", i+1, ct->sitetype[i].connections );
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
            ct->opt_QMzone = 0; // default
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
        // apparently residues are numbered in gromacs starting from 0 but written out as starting from 1
        // CHANGED IN GROMACS 4.6
     // ct->pool_site[i].resnr--;
    }

    searchkey(lines1, input1, "sitetypes", value, 1);
    split_string_into_string(value, ct->pool_size, dummy, "sitetypes");
    for (i = 0; i < ct->pool_size; i++)
    {
        ct->pool_site[i].type = atoi(dummy[i]);
        // for convenient numbering in charge-transfer.dat from 1 to #_different_sites
        ct->pool_site[i].type--;
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
        ct->pool_site[i]        = ct->sitetype[ct->pool_site[i].type]; // not sure if copying is that easy. it is
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
        snew(ct->pool_site[i].mass, ct->sitetype[ct->pool_site[i].type].atoms);
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
     // PRINTF("DATA %d %d %d %s %lf %d %d %d %lf %lf \n",
     //        ct->pool_site[i].atoms,      ct->pool_site[i].bonds, ct->pool_site[i].nochrs[0], ct->pool_site[i].nochr[0][0], ct->pool_site[i].extracharge[0],
     //        ct->pool_site[i].addchrs[0], ct->pool_site[i].homos, ct->pool_site[i].homo[0],   ct->pool_site[i].hubbard[0], ct->pool_site[i].lambda_i[0]);
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
                printf("check if adaptive QM zone works with extcharge list\n");
                exit(-1);
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
            // we allocate array a little bit larger than needed (natoms instead of extcharges) and let remaining enrtries blank.
            // This way we can build the intersection of these arrays in order to find the extcharges of the complex
            snew(ct->pool_site[i].extcharge, top_global->natoms);
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
    snew(ct->mass_cplx, ct->atoms_cplx);
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
#if (GMX_TENSORFLOW)
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
#else
       printf("Compiled without linking to the TensorFlow library.\n");
       printf("Thus, it is impossible to use a neural network model of electronic couplings!\n");
       printf("Exiting...\n");
       exit(-1);
#endif
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
            if (ct->survival < 0.99 || ct->survival > 1.01) // should be normalized, 0.01 tolerance
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
        exit(-1); // TODO: I somehow tried to adapt this feature.
                  // However, it was initially considered for one Orbital per site and
                  // I'm not sure if it will work with several FOs per site.
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
        ct->neg_imag_pot = 0; // default
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
                //variation of MO-energy if an other orbital on this site is getting charged.
                ct->hubbard[counter][counter+k] = ct->sic * (ct->site[i].hubbard[j] + ct->site[i].hubbard[k])*0.5;
                if (ct->do_epol == 1)
                {
                    //influence of electronic polarization
                    ct->hubbard[counter][counter+k] -= 1.0/(2.0*ct->site[i].radius*NM_TO_BOHR)*(1.0-1.0/EPSILON_OP);
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
        on -> old style of assigning qm-zone
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
            !strncmp((*(atoms->atomname[j])), "CQMT", 4) || // CQMT "branching" atoms
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
                 // PRINTF("%5d (%5s, type %d)\n", site->atom[counter], (*(atoms->atomname[site->atom[counter]])), site->atomtype[counter]+1);
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
                 // PRINTF("%5d (%5s, type %d)\n", site->atom[counter], (*(atoms->atomname[site->atom[counter]])), site->atomtype[counter]+1);
                    counter++;
                    ///// find capping for connection atoms /////
                    bond_length_best = 0.5; // = 0.5 nm
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
                        // find best capping in neighboring residue
                        if (atoms->resinfo[atoms->atom[QMCApool[k]].resind].nr != site->resnr && bond_length < bond_length_best)
                        {
                            // atomnumber of best QMCA (k) for QMCA (j) is saved in m
                            m                = QMCApool[k];
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
                        QMCAcounter--;                                // is the same connection
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
                                //find best capping in neighboring residue
                                if (atoms->resinfo[atoms->atom[QMCApool[k]].resind].nr != site->resnr && bond_length < bond_length_best)
                                {
                                    //atomnumber of best QMCA (k) for QMCA (j) is saved in n
                                    n                = QMCApool[k];
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
                         !strncmp((*(atoms->atomname[j])), "KQM", 3) || // Br
                         !strncmp((*(atoms->atomname[j])), "LQM", 3) || // Cl
                      // !strncmp((*(atoms->atomname[j])), "YQM", 3) || // Y was special pseudo atom
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
                     // case 'Y':
                     //     site->atomtype[counter] = 6;
                     //     break;
                        default: PRINTF("Unknown atom type for atom %d (%s), exiting!\n", j, (*(atoms->atomname[j]))); exit(-1);
                    }
                 // PRINTF("%5d (%5s, type %d)\n", site->atom[counter], (*(atoms->atomname[site->atom[counter]])), site->atomtype[counter]+1);
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
                 // PRINTF("%5d (%5s, type %d)\n", site->atom[counter], (*(atoms->atomname[site->atom[counter]])), site->atomtype[counter]+1);
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
     // PRINTF("\n");
        if (QMLAcounter != MMLAcounter)
        {
            PRINTF("Site %d found to have %d QM link atom(s) but %d MM link atom(s). \n", i, QMLAcounter, MMLAcounter);
            exit(-1);
        }
    } //end atom selection

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
            bond_length_best = 0.5; // = 0.5 nm
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
                //find best capping
                if (bond_length < bond_length_best)
                {
                    //index of best MMLA (k) for QMLA (j) is saved in m
                    m                = k;
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
                    environment = 1; // default: all atoms are environment
                 // PRINTF("j %d i %d atoms.atom[j].resnr %d ct->site[i].resnr %d \n" ,j, i, atoms.atom[j].resnr , ct->site[i].resnr);
                    // in the same residue -> further investigations
                    if (atoms->resinfo[atoms->atom[j].resind].nr == site->resnr)
                    {
                        for (k = 0; k < site->bonds; k++)
                        {
                            for (l = 0; l < site->addchrs[k]; l++)
                            {
                                // charge will be modified to electro-neutralize
                                if (!strcmp((*(atoms->atomname[j])), site->addchr[k][l]))
                                {
                                    site->modif_extcharge[k][l] = counter_array[i];
                                }
                            }
                            for (l = 0; l < site->nochrs[k]; l++)
                            {
                                // charge will be ignored
                                if (!strcmp((*(atoms->atomname[j])), site->nochr[k][l]))
                                {
                                    environment = 0;
                                }
                            }
                        }
                        for (k = 0; k < site->atoms; k++)
                        {
                            // exclude QM atoms
                            if (j == site->atom[k])
                            {
                                environment = 0;
                            }
                        }
                    }
                    else // atoms in neighboring residues are ignored if one connection atom of site i lies in this residue
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
            // this should successively reduce the charges in ct->extcharge_cplx
            k = find_intersection(top_global->natoms, ct->extcharge_cplx, ct->site[i].extcharge, ct->extcharge_cplx);
            // k elements are common in both arrays
            for (j = k; j < top_global->natoms; j++)
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
                        // if one of the extcharges of the complex is the same atom that was modified in the monomer calculation,
                        //   then also modify it in the complex calculation.
                        if (ct->extcharge_cplx[l] == ct->site[i].extcharge[ ct->site[i].modif_extcharge[j][k] ])
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
                        PRINTF("          modified extcharge for bond no. %5d: atom %5d - %s (residue %d)\n",
                               j+1, ct->site[i].extcharge[ct->site[i].modif_extcharge[j][k]]+1,
                               *(atoms->atomname[ct->site[i].extcharge[ct->site[i].modif_extcharge[j][k]]]),
                               atoms->resinfo[atoms->atom[ct->site[i].extcharge[ct->site[i].modif_extcharge[j][k]]].resind].nr+1 );
                    }
                }
            }
        }
        PRINTF("Complex: original group %d, (possibly) restricted to %d\n", ct->extcharges_cplx, counter_cplx);
     // ct->extcharges_cplx = counter_cplx;
     // if (counter_modif_cplx != 2*ct->sites)
     // {
            PRINTF("         number of atoms cutting QM/MM boundary = %d (there are %d sites)\n", counter_modif_cplx, ct->sites);
     //     exit(-1);
     // }
        for (j = 0; j < counter_modif_cplx; j++)
        {
            PRINTF("          modified extcharge no. %5d: atom %5d - %s (residue %d)\n",
                   j+1, ct->extcharge_cplx[ct->modif_extcharge_cplx[j]]+1,
                   *(atoms->atomname[ct->extcharge_cplx[ct->modif_extcharge_cplx[j]]]),
                   atoms->resinfo[atoms->atom[ct->extcharge_cplx[ct->modif_extcharge_cplx[j]]].resind].nr+1);
        }

    } // end QMMM>0


    ///// JOB SPECIFIC PREPARATIONS /////

    // NEGF data
    /* TODO: THIS WAS ONLY COMMENTED OUT BECAUSE I HAVE NO IDEA ABOUT THESE CALCULATIONS.
             MAYBE IMPLEMENT NEGF CALCUlATIONS ALSO IN THE CURRENT CODE
    if (ct->jobtype == cteNEGFLORENTZ || ct->jobtype == cteNEGFLORENTZNONSCC)
    {
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
    }

    // NEGF initialization including initial density matrix
    if (ct->jobtype == cteNEGFLORENTZ || ct->jobtype == cteNEGFLORENTZNONSCC)
    {
        PRINTF("Initializing the NEGF calculation\n");
#if GMX_MPI
        if (ct_mpi_rank == 0)
#endif
        negf_init_arrays(ct->negf_arrays, &(ct->rk_timestep), ct->wf);
    }
    */

    // DO HERE PREPARATIONS FOR BORNOPPENHEIMER!
    if (ct->jobtype == cteBORNOPPENHEIMER)
    {
        snew(ct->born_overlap, ct->dim);
    }

    // DO HERE PREPARATIONS FOR TFL! // imported by J.J. Kranz Frenkel Version 2018
    if (ct->jobtype == cteSCCDYNAMIC || ct->jobtype == cteCPFSCCDYNAMIC || ct->jobtype == cteDDBSCCDYNAMIC ||
        ct->jobtype == cteDLZSH || ct->jobtype == cteALZSH || ct->jobtype == cteDFSSH || ct->jobtype == cteJFSSH ||
        ct->jobtype == cteBCJFSSH || ct->jobtype == cteSCRDFSSH || ct->jobtype == cteCCFSSH || ct->jobtype == cteGFSH ||
        ct->jobtype == cteDISH || ct->jobtype == cteSCRDFSSH || ct->jobtype == cteTULLYLOC)
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

        // add nonadiabatic coupling vector between pre_surf and hop_surf, W.X, Nov.2018
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

        // adiabatic energies at previous time step and two time steps before
        snew(ct->ev_adiab_old2, ct->dim);
        snew(ct->ev_adiab_old, ct->dim);
        snew(ct->tfl_old_ham, SQR(ct->dim));
        snew(ct->tfl_old_ham_adiab, SQR(ct->dim));

     // srand48(123);
     // srand((unsigned)time(NULL));
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

#if GMX_MPI
    printf("Completed charge transfer initialization at rank %d\n", ct_mpi_rank);
#else
    printf("Completed charge transfer initialization\n");
#endif

/*
    // Print out (nearly) all stuff that was read in
    for (k=0; k<ct->pool_size; k++)
    {
        s=ct->pool_site[k];
        PRINTF("%d %d %d %d %d %d %d %d %d %f %d %d\n",
               s.type, s.resnr, s.atoms, s.bonds, s.connections, s.homos,
               s.extcharges, s.nel, s.norb, s.radius, s.do_scc, s.active);
        for (i=0; i<s.atoms; i++)
        {
            PRINTF("%d %d \n", s.atom[i], s.atomtype[i]);
        }
        for (i=0; i<s.bonds; i++)
        {
            PRINTF("%d %d %d \n", s.QMLA[i], s.MMLA[i], s.nochrs[i]);
            for (j=0; j<s.nochrs[i]; j++)
                PRINTF("%s\n", s.nochr[i][j]);
            for (j=0; j<s.addchrs[i]; j++)
                PRINTF("%s %d\n", s.addchr[i][j], s.modif_extcharge[i][j]);
        }
        for (i=0; i<s.extcharges; i++)
            PRINTF("%d ", s.extcharge[i]);
        for (i=0; i<s.homos; i++)
            PRINTF("%d %f %f\n", s.homo[i], s.hubbard[i], s.lambda_i[i]);
    }
    PRINTF("%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d \n",
           ct->jobtype, ct->interval, ct->qmmm, ct->sitetypes, ct->sites, ct->dim, ct->is_hole_transfer, ct->atoms_cplx,
           ct->extcharges_cplx, ct->modif_extcharges_cplx, ct->n_avg_ham, ct->do_lambda_i, ct->do_epol, ct->decoherence, ct->pool_size, ct->opt_QMzone);
    PRINTF("%f %f %f %f %f    %f %f %f\n",
           ct->offdiag_scaling, ct->sic, ct->esp_scaling_factor, ct->telec, ct->fermi_kt, ct->efield[0], ct->efield[1], ct->efield[2]);
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

void init_dftb_stub(dftb_t *dftb, charge_transfer_t *ct)
{
    snew(dftb->nl, ct->atoms_cplx);
    snew(dftb->nl[0], SQR(ct->atoms_cplx));
    for (int j = 1; j < ct->atoms_cplx; j++)
    {
        dftb->nl[j] = dftb->nl[0] + j * ct->atoms_cplx;
    }
}

