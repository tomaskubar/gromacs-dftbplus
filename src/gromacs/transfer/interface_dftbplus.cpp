#include "gromacs/transfer/transfer.h"

typedef struct Context {
  bool               pme;
  int                n;
  const t_commrec*   cr;
  QMMM_rec_transfer* qr;
  t_nrnb*            nrnb;
  gmx_wallcycle_t    wcycle;
  real               rcoul;
  real               ewaldcoeff_q;
} Context;


/***********************************************
 * INITIALIZE THE DFTBPLUS QMMM_REC STRUCTURES *
 ***********************************************/

#if (GMX_MPI)
void init_qmmmrec_transfer(std::unique_ptr<QMMM_rec_transfer>* dftbplus_phase1,
                           std::unique_ptr<QMMM_rec_transfer>& dftbplus_phase2,
                           charge_transfer_t*  ct,
                           const t_commrec*    cr,
                           const real          rcoulomb,
                           const real          ewald_rtol,
                           const int           ct_mpi_rank)
#else
void init_qmmmrec_transfer(std::unique_ptr<QMMM_rec_transfer>* dftbplus_phase1,
                           std::unique_ptr<QMMM_rec_transfer>& dftbplus_phase2,
                           charge_transfer_t*  ct,
                           const t_commrec*    cr,
                           const real          rcoulomb,
                           const real          ewald_rtol)
#endif
{
/* Initializes the DFTB+ calculations, copies information that was gathered in init_charge_transfer(), allocates memory */

// PARAMETERS:
// mdatoms         = (in) variable of GROMACS data-structure-type (we need atomic masses)
// dftbplus_phase1 = (out) main QMMMrec data structure for DFTB calculations, phase1
// dftbplus_phase2 = (out) main QMMMrec data structure for DFTB calculations, phase2
// ct              = (in) main data structure with information about the charge transfer calculation
// ct_mpi_rank     = (in) the rank of the process in MPI calculations (each process in MPI calculations can be destinguished by its value of ct_mpi_rank)
//////////////////////////////////////////////////////////////////////////////////////////
 // const char *type_symbols = "chnosfklbpi";        // fe is z 16.06.2016 because of one letter relation between element and code
 //                                               // K is Br, L is Cl 13Sep2018 Daniel for halogens

#if (GMX_MPI)
    printf("DFTB+ initialization at rank %d\n", ct_mpi_rank);
#else
    printf("DFTB+ initialization\n");
#endif

    // New (DFTB+) variables
    char ptrElement[20][3]; // element names may have up to 2 characters (+1 null termination)
    int *ptrSpecies;
    int *atomicNumber;
    int atomicNumberBySpecies[20];
    char chemSymbol[55][3] = { "XX",
          "H",                                      "He",
          "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne",
          "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar",
          "K",  "Ca",
          "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
                      "Ga", "Ge", "As", "Se", "Br", "Kr",
          "Rb", "Sr",
          "Y ", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
                      "In", "Sn", "Sb", "Te", "I",  "Xe" };
 // static DftbPlus  calculator;
    DftbPlusInput    input;
    DftbPlusAtomList atomList;

    // Deal with the atom types for each site
    for (int site = 0; site < ct->sites; site++)
    {
        std::unique_ptr<QMMM_QMrec_transfer>& qm = dftbplus_phase1[site]->qm;
        snew(qm->dftbContext, 1);
        initialize_context_transfer(qm->dftbContext, qm->nrQMatoms_get(), qm->qmmm_variant_get(), dftbplus_phase1[site].get(), rcoulomb, ewald_rtol, cr);

        snew(qm->dpcalc, 1);
     // qm->dpcalc = &calculator;

        /* Initialize the DFTB+ calculator */
        char filename[16];
        sprintf(filename, "dftb1_%d.out", site);
        dftbp_init(qm->dpcalc, filename);
        printf("DFTB+ calculator for phase 1, site %d has been created!\n", site);

        /* Parse the input file and store the input-tree */
        dftbp_get_input_from_file(qm->dpcalc, "dftb1.hsd", &input);
        printf("DFTB+ input has been read!\n");

        /* Pass the list of QM atoms to DFTB+ */
        int nAtom = qm->nrQMatoms_get();
        snew(ptrSpecies, nAtom);
        snew(atomicNumber, nAtom);
        // read the atomic numbers of QM atoms, and set up the lists
        int nSpecies = 0;
        for (int i=0; i<nAtom; i++)
        {
            atomicNumber[i] = qm->atomicnumberQM_get(i);
            // is this a new chemical element?
            bool newElement = true;
            for (int j=0; j<i; j++)
                if (atomicNumber[i] == atomicNumber[j])
                    newElement = false;
            // if it is a new element, introduce it in the list
            if (newElement)
            {
                atomicNumberBySpecies[nSpecies] = atomicNumber[i];
                nSpecies++;
                ptrSpecies[i] = nSpecies; // this numbering will start at 1 (and not 0)
            }
            else
            {
                for (int k=0; k<nSpecies; k++)
                    if (atomicNumber[i] == atomicNumberBySpecies[k])
                        ptrSpecies[i] = k+1; // because numbering starts at 1
            }
        }
        // assemble the list of chemical species
        for (int k=0; k<nSpecies; k++)
            strcpy(ptrElement[k], chemSymbol[atomicNumberBySpecies[k]]);
        // check what is being passed to DFTB+
        printf("This is being passed from Gromacs to DFTB+:\n");
        printf("No. of QM atoms = %d (site %d)\n", nAtom, site);
        printf("No. of chem. species = %d (site %d):", nSpecies, site);
        for (int k=0; k<nSpecies; k++)
            printf(" %d=%s", k+1, ptrElement[k]);
        printf("\n");
        printf("Species by atom:\n");
        for (int i=0; i<nAtom; i++)
            printf("Atom %d is species %d\n", i+1, ptrSpecies[i]);

        // finally, call the DFTB+ routine
        dftbp_get_atom_list(&atomList, &nAtom, &nSpecies, (char *) ptrElement, ptrSpecies);
        printf("DFTB+ has obtained the list of QM atoms!\n");

        sfree(atomicNumber);
        sfree(ptrSpecies);

        /* Set up the calculator by processing the input tree */
        dftbp_process_input(qm->dpcalc, &input, &atomList);
        printf("DFTB+ input has been processed!\n");

        /* Register the callback functions which calculate
         * the external potential and its gradient
         */
        dftbp_register_ext_pot_generator(qm->dpcalc,
                                         qm->dftbContext,
                                         calcqmextpot_transfer,
                                         calcqmextpotgrad_transfer);

    }

    // Deal with the atom types for the complex (phase 2)
    std::unique_ptr<QMMM_QMrec_transfer>& qm = dftbplus_phase2->qm;
    snew(qm->dftbContext, 1);
    initialize_context_transfer(qm->dftbContext, qm->nrQMatoms_get(), qm->qmmm_variant_get(), dftbplus_phase2.get(), rcoulomb, ewald_rtol, cr);

    snew(qm->dpcalc, 1);
 // qm->dpcalc = &calculator;

    /* Initialize the DFTB+ calculator */
    dftbp_init(qm->dpcalc, "dftb2.out");
    printf("DFTB+ calculator for phase 2 has been created!\n");

    /* Parse the input file and store the input-tree */
    dftbp_get_input_from_file(qm->dpcalc, "dftb2.hsd", &input);
    printf("DFTB+ input has been read!\n");

    /* Pass the list of QM atoms to DFTB+ */
    int nAtom = qm->nrQMatoms_get();
    snew(ptrSpecies, nAtom);
    snew(atomicNumber, nAtom);
    // read the atomic numbers of QM atoms, and set up the lists
    int nSpecies = 0;
    for (int i=0; i<nAtom; i++)
    {
        atomicNumber[i] = qm->atomicnumberQM_get(i);
        // is this a new chemical element?
        bool newElement = true;
        for (int j=0; j<i; j++)
            if (atomicNumber[i] == atomicNumber[j])
                newElement = false;
        // if it is a new element, introduce it in the list
        if (newElement)
        {
            atomicNumberBySpecies[nSpecies] = atomicNumber[i];
            nSpecies++;
            ptrSpecies[i] = nSpecies; // this numbering will start at 1 (and not 0)
        }
        else
        {
            for (int k=0; k<nSpecies; k++)
                if (atomicNumber[i] == atomicNumberBySpecies[k])
                    ptrSpecies[i] = k+1; // because numbering starts at 1
        }
    }
    // assemble the list of chemical species
    for (int k=0; k<nSpecies; k++)
        strcpy(ptrElement[k], chemSymbol[atomicNumberBySpecies[k]]);
    // check what is being passed to DFTB+
    printf("This is being passed from Gromacs to DFTB+:\n");
    printf("No. of QM atoms = %d (phase2)\n", nAtom);
    printf("No. of chem. species = %d (phase2):", nSpecies);
    for (int k=0; k<nSpecies; k++)
        printf(" %d=%s", k+1, ptrElement[k]);
    printf("\n");
    printf("Species by atom:\n");
    for (int i=0; i<nAtom; i++)
        printf("Atom %d is species %d\n", i+1, ptrSpecies[i]);

    // finally, call the DFTB+ routine
    dftbp_get_atom_list(&atomList, &nAtom, &nSpecies, (char *) ptrElement, ptrSpecies);
    printf("DFTB+ has obtained the list of QM atoms!\n");

    sfree(atomicNumber);
    sfree(ptrSpecies);

    /* Set up the calculator by processing the input tree */
    dftbp_process_input(qm->dpcalc, &input, &atomList);
    printf("DFTB+ input has been processed!\n");

    /* Register the callback functions which calculate
     * the external potential and its gradient
     */
    dftbp_register_ext_pot_generator(qm->dpcalc,
                                     qm->dftbContext,
                                     calcqmextpot_transfer,
                                     calcqmextpotgrad_transfer);
    // END DFTB+ ATOMTYPE LIST

#if (GMX_MPI)
    printf("Completed DFTB+ initialization at rank %d\n", ct_mpi_rank);
#else
    printf("Completed DFTB+ initialization\n");
#endif

} // init_qmmmrec_transfer()


/*
#if (GMX_MPI)
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

 // printf("pbc %d\n", pbc.ePBCDX);

 // switch (pbc.ePBCDX)
 // {
 //     case epbcdxRECTANGULAR:
 //         printf("tric\n");
 //     case epbcdxTRICLINIC:
 //         printf("rect\n");
 // }

 // for (i = 0; i < DIM; i++)
 // printf("hbox %f fbox %f mhbox %f \n", pbc.hbox_diag[i], pbc.fbox_diag[i],  pbc.mhbox_diag[i]);
    rcoulomb2 = dftb->rcoulomb_pme * dftb->rcoulomb_pme;

    //printf("  creating PME DFTB neighborlists...");
    for (i = 0; i < ct->sites; i++)
#if (GMX_MPI)
    {
        if (i % ct_mpi_size == ct_mpi_rank)
#endif
    {
        dftb1 = dftb->phase1[i];
        nn    = dftb1.nn;
        ne    = dftb1.ne;
        // create a neighborlist for PME calculation for QM/MM DFTB
        //   in a naive O(N^2) way:
        //     double loop over QM atoms and over MM atoms

        // do it for every QM atom
        for (j = 0; j < nn; j++)
        {
            counter = 0;
            // check every MM atom
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
         // printf("bond %17.8f \n", norm2(bond));
            dftb1.neighbors_pme[j] = counter;
         // fprintf(stderr, "NS for PME/DFTB: atom %d in fragment %d has %d MM neighbors\n", j+1, i+1, counter);
        }
    }
#if (GMX_MPI)
    }
#endif

#if (GMX_MPI)
    if (ct_mpi_rank == 0)
    {
#endif
    // do it here for the complex
    dftb2 = dftb->phase2;
    nn    = dftb2.nn;
    ne    = dftb2.ne;
    // do it for every QM atom
    for (j = 0; j < nn; j++)
    {
        counter = 0;
        // check every MM atom
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
      // fprintf(stderr, "NS for PME/DFTB: atom %d in complex has %d MM neighbors\n", j+1, counter);
    }
#if (GMX_MPI)
}
#endif
} // do_neighborlist_for_dftb()
*/
