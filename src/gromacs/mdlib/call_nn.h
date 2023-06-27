/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017,2018, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
//#include "gromacs/mdlib/qmmm.h"

#ifndef GMX_MDLIB_QM_DFTBPLUS_H
#define GMX_MDLIB_QM_DFTBPLUS_H

#include <tensorflow/c/c_api.h>

#define HARTREE_TO_EV     (27.211396132)

typedef struct TF_Model {
    TF_Graph*   Graph;
    TF_Status*  Status;
    TF_Session* Session;
    int         NumInputs;
    TF_Output*  Input;
    int         NumOutputs;
    TF_Output*  Output;
    TF_Tensor** InputValues;
    TF_Tensor** OutputValues;
    int         nAtoms;
    int*        atomicNumbers;
}

void
init_dftbplus(QMMM_QMrec*       qm,
              QMMM_rec*         qr,
           // const t_forcerec* fr,
              const t_inputrec* ir,
              const t_commrec*  cr);
           // gmx_wallcycle_t   wcycle);

real
call_dftbplus(QMMM_rec*         qr,
              const t_commrec*  cr,
              QMMM_QMrec*       qm,
              const QMMM_MMrec& mm,
              rvec              f[],
              rvec              fshift[],
              t_nrnb*           nrnb,
              gmx_wallcycle_t   wcycle);

struct Context;

void initialize_context(Context*          cont,
                        int               nrQMatoms,
                        int               qmmm_variant,
                        QMMM_rec*         qr_in,
                     // const t_forcerec* fr_in,
                        const t_inputrec* ir_in,
                        const t_commrec*  cr_in);
                     // gmx_wallcycle_t   wcycle_in);

void calcQMextPotPME(Context *cont, double *q, double *extpot);

extern "C" void calcqmextpot(void *refptr, double *q, double *extpot);

extern "C" void calcqmextpotgrad(void *refptr, gmx_unused double *q, double *extpotgrad);


/*------------------------------------------------------------------------------------------------*/
/*  DFTB+: general package for performing fast atomistic simulations                              */
/*  Copyright (C) 2006 - 2020  DFTB+ developers group                                             */
/*                                                                                                */
/*  See the LICENSE file for terms of usage and distribution.                                     */
/*------------------------------------------------------------------------------------------------*/
#ifndef __DFTBPLUS_H__
#define __DFTBPLUS_H__

extern "C" {

/**
 * Type containing the list of atoms and species to be passed to DFTB+.
 *
 * Used by DFTB+ as an opaque handler. Do not manipulate the content of this type directly!
 */
typedef struct DftbPlusAtomList {
  void *pDftbPlusAtomList;
} DftbPlusAtomList;


/**
 * Type containing the DFTB+ input tree.
 *
 * Used by DFTB+ as an opaque handler. Do not manipulate the content of this type directly!
 */
typedef struct {
  void *pDftbPlusInput;
} DftbPlusInput;


/**
 * Type containing the DFTB+ calculator.
 *
 * Used by DFTB+ as an opaque handler. Do not manipulate the content of this type directly!
 */
typedef struct DftbPlus {
  void * pDftbPlus;
} DftbPlus;


/**
 * Callback function signature for calculating external population dependant potential.
 *
 * DFTB+ would call it whenever the population has changed and the population dependant external
 * potential must be recalculated.
 *
 * \param refptr[in] Reference pointer. This is the pointer you passed to DFTB+ when the call-back
 *     functions had been registered. You can use it to find the data you want to use to calculate
 *     the external potential.

 * \param dqatom[in] Population difference with respect to reference population (usually the neutral
 *     atom). Shape: [natom]. Note: Population means electrons, so a positive number indicates
 *     electron excess.
 *
 * \param extpotatom[out] Potential at the position of each qm-atom. Shape: [natom]. Note: It
 *     should be the potential as felt by an electron (negative potential value means attraction
 *     for an electron). Unit: Hartree.
 */
typedef void (*ExtPotFunc)(void *refptr, double *dqatom, double *extpotatom);


/**
 * Callback function signature for calculating external population dependant potential gradient.
 *
 * DFTB+ would call it whenever the forces are calculated and the population dependant external
 * potential must be taken into account.
 *
 * \param refptr[in] Reference pointer. This is the pointer you passed to DFTB+ when the call-back
 *     functions had been registered. You can use it to find the data you want to use to calculate
 *     the external potential gradient.
 *
 * \param dqatom[in] Population difference with respect to reference population (usually the
 *     neutral atom). Shape: [natom]. Note: Population means electrons, so a positive number
 *     indicates electron excess.
 *
 * \param extpotatomgrad[out] Potential gradient at the position of each qm-atom. Shape: [natom,
 *     3]. Note: It should be the gradient of the potential as felt by an electron (negative
 *     potential value means attraction for an electron). Unit: Hartree/Bohr.
 */
typedef void (*ExtPotGradFunc)(void *refptr, double *dqatom, double *extpotatomgrad);


/**
 * Returns current version of the DFTB+ API
 *
 * \param[out] major major release number
 *
 * \param[out] minor minor release number
 *
 * \param[out] patch patch of release
 */
void dftbp_api(int *major, int *minor, int *patch);


/**
 * Initializes a DFTB+ calculator.
 *
 * \param[inout] instance  Handler of DFTB+ instance.
 *
 * \param[in] outputfilename  Name of the file, where the DFTB+ screen output should be written.
 *     If you pass NULL, it will be written to standard output. If you pass any other file, it will
 *     be open, and the file will be written there. Pass "/dev/null" to suppress output.
 */
void dftbp_init(DftbPlus *instance, const char *outputfilename);


/**
 * Finalizes a DFTB+ calculator
 *
 *  \param[inout] instance  Handler of the DFTB+ instance.
 */
void dftbp_final(DftbPlus *instance);


/**
 * Obtain no. of atoms and list of species from the MM program
 */
void dftbp_get_atom_list(DftbPlusAtomList *atomListHandler, int *nAtomC, int *nSpeciesC, char *elementC, int *species);


/**
 * Fills up a DFTB+ input tree from a HSD input file.
 *
 * \param[inout] instance Handler of the DFTB+ instance.
 *
 * \param[in] filename Name of the file containing the HSD-input for DFTB+.
 *
 * \param[inout] input Handler containing the input tree parsed from the input file.
 */
void dftbp_get_input_from_file(DftbPlus *instance, const char *filename, DftbPlusInput *input);


/**
 * Sets up the calculator by processing a given input tree.
 *
 * \param[inout] instance Handler of the DFTB+ instance.
 *
 * \param[inout] atomListHandler The structure containing the list of atoms and species.
 *
 * \param[inout] input The tree containing the DFTB+ input. On return, it contains the tree
 *     extended by all the default options set by the parser.
 */
void dftbp_process_input(DftbPlus *instance, DftbPlusInput *input, DftbPlusAtomList *atomListHandler);


/**
 * Sets a constant (charge independent) external potential.
 *
 * Note: It should be the potential as felt by an electron (negative potential value means
 *     attraction for an electron).
 *
 * \param[inout] instance Handler of the DFTB+ instance.
 *
 * \param[in] extpot External potential at the position of each atom. Shape: [natom]. Unit: Hartree.
 *
 * \param[in] extpotgrad Gradient of the external potential at each atom. Shape: [natom, 3].  Unit:
 *     Hartree/Bohr. This parameter is optional, you can pass NULL if you did not ask DFTB+ to
 *     calculate forces.
 */
void dftbp_set_external_potential(DftbPlus *instance, const double *extpot,
                                  const double *extpotgrad);


/**
 * Registers callback functions for population dependant external potential calculation.
 *
 * \param[inout] instance Handler of the DFTB+ instance.
 *
 * \param[in] refptr Arbitrary pointer. DFTB+ will pass back this pointer unaltered when calling
 *     the registered functions. You can typically use it to pass a pointer to the data struct
 *     or class which contains the necessary data for the potential calculation. If your data
 *     in in the global space and you don not need it, pass an arbitrary pointer, e.g. NULL.
 *
 * \param[in] extpot Function pointer to the call-back function which DFTB+ should call, whenever
 *     the population dependant external potential should be calculated.
 *
 * \param[in] extpotgrad Function pointer to the call-back function which DFTB+ should call,
 *     whenever the gradient of the population dependant external potential should be calculated.
 */
void dftbp_register_ext_pot_generator(DftbPlus *instance, const void *refptr, ExtPotFunc extpot,
                                      ExtPotGradFunc extpotgrad);


/**
 * Sets actual atom coordinates.
 *
 * \param[inout] instance Handler of the DFTB+ instance.
 *
 * \param[in] coords Coordinates of the atoms. Shape: [natom, 3]. Unit: Bohr.
 */
void dftbp_set_coords(DftbPlus *instance, const double *coords);


/**
 * Sets actual atom coordinates and lattice vectors.
 *
 * \param[inout] instance Handler of the DFTB+ instance.
 *
 * \param[in] coords Coordinates of the atoms in atomic units. Shape: [natom, 3]. Unit: Bohr.
 *
 * \param[in] latvecs Lattice vectors Shape: [3, 3]. Unit: Bohr.
 */
void dftbp_set_coords_and_lattice_vecs(DftbPlus *instance, const double *coords,
                                       const double *latvecs);


/**
 * Sets actual atom coordinates and lattice vectors.
 *
 * \param[inout] instance Handler of the DFTB+ instance.
 *
 * \param[in] coords Coordinates of the atoms in atomic units. Shape: [natom, 3]. Unit: Bohr.
 *
 * \param[in] latvecs Lattice vectors Shape: [3, 3]. Unit: Bohr.
 *
 * \param[in] origin Coordinate origin Shape: [3]. Unit: Bohr.
 */
void dftbp_set_coords_lattice_origin(DftbPlus *instance, const double *coords,
                                       const double *latvecs, const double *origin);

/**
 * Queries the nr. of atoms in the system.
 *
 * \param[inout] instance Handler of the DFTB+ instance.
 *
 * \return Nr. of atoms
 */
int dftbp_get_nr_atoms(DftbPlus *instance);


/**
 * Queries the energy of the current geometry
 *
 * \param[inout] instance Handler of the DFTB+ instance.
 *
 * \param[out] mermin_energy  Mermin free energy of the current geometry. Unit: Bohr.
 */
void dftbp_get_energy(DftbPlus *instance, double *mermin_energy);


/**
 * Queries the gradients of the current geometry.
 *
 * \param[inout] instance Handler of the DFTB+ instance.
 *
 * \param[out] gradients Gradients (not forces!) on each atom. Shape [natom, 3]. Unit: Hartree/Bohr.
 */
void dftbp_get_gradients(DftbPlus *instance, double *gradients);


/**
 * Queries the stress tensor of the current periodic box
 *
 * \param[inout] instance Handler of the DFTB+ instance.
 *
 * \param[out] stresstensor Stress Tensor for the periodic box. Shape [3, 3]. Unit: Pascals.
 */
void dftbp_get_stress_tensor(DftbPlus *instance, double *stresstensor);


/**
 * Queries the net populations on the atoms.
 *
 * \param[inout] instance Handler of the DFTB+ instance.
 *
 * \param[out] charges Net charges on each atom.  Shape [natom]. Sign convention: Electron has
 *     negative charge, so negative values indicate electron excess.
 */
void dftbp_get_gross_charges(DftbPlus *instance, double *charges);

}

#endif

#endif
