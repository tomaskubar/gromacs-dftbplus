/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017 by the GROMACS development team.
 * Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include "qmmm.h"

#include "config.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>

#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/ewald/ewald_utils.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/nblist.h"
#include "gromacs/nbnxm/grid.h"
#include "gromacs/nbnxm/gridset.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/nbnxm/pairlist.h"
#include "gromacs/nbnxm/pairlistset.h"
#include "gromacs/nbnxm/pairlistsets.h"
#include "gromacs/nbnxm/pairsearch.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/mtop_lookup.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"


// convert the "charge transfer" atomtype to "standard" atomic number
void atomtype_to_atomicnumber(int atomtype, int& atomicnumber, double& atomicmass)
{
    switch(atomtype)
    {
        case  0: atomicnumber =  6; atomicmass =  12.01; break; // C
        case  1: atomicnumber =  1; atomicmass =   1.01; break; // H
        case  2: atomicnumber =  7; atomicmass =  14.01; break; // N
        case  3: atomicnumber =  8; atomicmass =  16.00; break; // O
        case  4: atomicnumber = 16; atomicmass =  32.06; break; // S
        case  5: atomicnumber =  9; atomicmass =  19.00; break; // F
        case  6: atomicnumber = 35; atomicmass =  79.90; break; // K = Br
        case  7: atomicnumber = 17; atomicmass =  35.45; break; // L = Cl
        case  8: atomicnumber =  5; atomicmass =  10.81; break; // B
        case  9: atomicnumber = 15; atomicmass =  30.97; break; // P
        case 10: atomicnumber = 53; atomicmass = 126.90; break; // I
        default: printf("Unknown atomtype %d, exiting!\n", atomtype); exit(-1);
    }
}


// from the "standard" qmmm.cpp
void put_cluster_in_MMlist_verlet(int                             ck, // cluster number
                                  int                             na_ck, // # of atoms in cluster
                                  int                             nrQMatoms,
                                  const int*                      indexQM,
                                        gmx::ArrayRef<const int>  atomIndices,
			                      int*                            shiftMMatom,
                                  // ^ also has a role of "bool* isMMatom"
				                  t_pbc*                          pbc,
				                  const rvec*                     x);


// Update QM and MM coordinates in the QM/MM data structures.
// New version of the function:
// Update the coordinates of the MM atoms on the short-range neighborlist!
// The NBlist needs to have been created previously by either group or Verlet scheme.
//
// THIS WILL GO TO THE "PREPARE" PROCEDURE MAYBE?
//
void QMMM_rec_transfer::update_QMMM_coord(const t_commrec*  cr,
                                          const rvec*       shift_vec,
                                          const rvec        x[],
                                          const t_mdatoms*  md,
                                          const matrix      box)
{
    // Shifts the QM and MM atoms into the central box and
    //   stores the shifted coordinates in the coordinate arrays of QMMMrec.
    // These coordinates are passed on the QM subroutines.
    //
    // Only MM atoms up to the distance fr->rcoulomb from the respective
    //   nearest QM atoms are considered;
    // in case fr->rcoulomb == 0. is detected,
    //   all of the MM atoms are considered.

 // QMMM_QMrec_transfer& qm_ = qm;
 // QMMM_MMrec_transfer& mm_ = mm;
    real rcut = qm->rcoulomb > 0.1 ? qm->rcoulomb : 999999.; // infinity
    std::vector<bool> isCurrentMMatom;
    isCurrentMMatom.resize(mm->nrMMatoms_nbl);

    // shift the QM atoms into the central box
    for (int i = 0; i < qm->nrQMatoms; i++)
    {
        rvec_sub(x[qm->indexQM[i]], shift_vec[qm->shiftQM[i]], qm->xQM[i]);
    }

    // copy box size
    copy_mat(box, qm->box);

    // initialize PBC for MM coordinate manipulation
    t_pbc pbc;
    ivec null_ivec;
    clear_ivec(null_ivec);
    set_pbc_dd(&pbc, pbcType, DOMAINDECOMP(cr) ? cr->dd->numCells : null_ivec, false, box);

 // for (int s = 0; s < pbc.ntric_vec; s++)
 // {
 //     printf("SHIFT[%2d] = %d %d %d\n", s, pbc.tric_shift[s][0], pbc.tric_shift[s][1], pbc.tric_shift[s][2]);
 //  // printf("SHIFT[%2d] = %8.5f %8.5f %8.5f\n", s, pbc.tric_vec[s][0], pbc.tric_vec[s][1], pbc.tric_vec[s][2]);
 // }

    // DECIDE IF WE WANT TO APPLY A CUTOFF ON THE ATOMS FROM THE SR NEIGHBORLIST!

    // DO WE NEED TO RE-ALLOCATE THE ARRAYS TO BE FILLED?
    //   YES!
    //
    // FIRST, IDENTIFY THE MM ATOMS UP TO CUTOFF AT THIS STEP AND COUNT THEM:

    // Among the atoms found as candidates for being MM atoms in neighborsearching,
    // find those that are within electrostatics cut-off.
    // For the cutoff, use the value "rcut"
    int nrMMatoms = 0;
    for (int i = 0; i < mm->nrMMatoms_nbl; i++)
    {
	    isCurrentMMatom[i] = false;
	 // printf("DEBUG_MM TEST %5d %5d %2d", i, mm_.indexMM_nbl[i], mm_.shiftMM_nbl[i]);
	    // loop over all QM atoms here
	    for (int q=0; q<qm->nrQMatoms; q++)
	    {
            rvec bond;
            pbc_dx_aiuc(&pbc, x[qm->indexQM[q]], x[mm->indexMM_nbl[i]], bond);
	     // printf(" %8.5f\n", norm(bond));
	        if (norm(bond) < rcut)
            {
	            isCurrentMMatom[i] = true;
	            nrMMatoms++;
	         // printf("DEBUG_MM %5d %5d %2d %6.4f\n", i, mm_.indexMM_nbl[i], mm_.shiftMM_nbl[i], distance);
                break;
            }
	    }
	 // printf("\n");
    }
 // printf("Number of actual    MM atoms in the current MD step               : %d\n", nrMMatoms);

    // ALLOCATION
    mm->nrMMatoms = nrMMatoms;
    mm->indexMM.resize(nrMMatoms);
    mm->MMcharges.resize(nrMMatoms);
    mm->shiftMM.resize(nrMMatoms);
    mm->xMM.resizeWithPadding(nrMMatoms);

    int index = 0; // runs over the identified MM atoms
    for (int i = 0; i < mm->nrMMatoms_nbl; i++)
    {
	    if (isCurrentMMatom[i])
	    {
	        // Add to list!
	        mm->indexMM[index] = mm->indexMM_nbl[i];

	        // Also add charge
	        mm->MMcharges[index] = md->chargeA[mm->indexMM[index]] * mm->scalefactor;

            // Having obtained the shift at NS time (update_qmmmrec),
            //   merely copy it here to shiftMM[]
	        mm->shiftMM[index] = mm->shiftMM_nbl[i];

	        // one MM atom found => increment index */
	        index++;
	    }
    }

    // also shift the MM atoms into the central box

 //   for (int a=0; a<45; a++)
 //     printf("SHIFT %2d: %7.3f %7.3f %7.3f\n", a,
 //     fr->shift_vec[a][XX], fr->shift_vec[a][YY], fr->shift_vec[a][ZZ]);

    for (int ind = 0; ind < mm->nrMMatoms; ind++)
    {
        rvec_sub(x[mm->indexMM[ind]], shift_vec[mm->shiftMM[ind]], mm->xMM[ind]);
 //     printf("COORD MM %4d %2d\n", mm_.indexMM[ind], mm_.shiftMM[ind]);
    }

    // For DFTB, also update the coordinates of *all* of the MM atoms,
    //   not only those on the short-range neighborlist.
    // Do not shift the MM atoms into the central box!
    //   It might break the calculation of the surface correction in the Ewald sum.
    for (int i = 0; i < mm->nrMMatoms_full; i++)
    {
        copy_rvec(x[mm->indexMM_full[i]], mm->xMM_full[i]);
    }
} // update_QMMM_coord

void QMMM_QMrec_transfer::init_QMrec(const std::vector<int>& qmAtoms,
                                  // const int               QMcharge_in,
                                     const std::vector<int>& atomicnumberQM_in,
                                     const t_inputrec*       ir)
{
    nrQMatoms = qmAtoms.size();
    printf("nrQMatoms = %d\n", nrQMatoms);
    snew(xQM, nrQMatoms);
    snew(indexQM, nrQMatoms);
    snew(shiftQM, nrQMatoms);
    for (int i = 0; i < nrQMatoms; i++)
    {
        indexQM[i] = qmAtoms[i];
    }

    snew(atomicnumberQM, nrQMatoms);
 // nelectrons = 0;
    for (int i = 0; i < nrQMatoms; i++)
    {
        atomicnumberQM[i]  = atomicnumberQM_in[i];
     // nelectrons        += atomicnumberQM[i];
    }

 // QMcharge      = QMcharge_in;
 // nelectrons   -= QMcharge;

    rcoulomb      = ir->rcoulomb;
    ewaldcoeff_q  = calc_ewaldcoeff_q(ir->rcoulomb, ir->ewald_rtol);
    epsilon_r     = ir->epsilon_r;

    snew(pot_qmmm, nrQMatoms);
    snew(pot_qmqm, nrQMatoms);

} // init_QMrec

int QMMM_QMrec_transfer::nrQMatoms_get() const
{
    return nrQMatoms;
}

int QMMM_QMrec_transfer::qmmm_variant_get() const
{
    return qmmm_variant;
}

double QMMM_QMrec_transfer::xQM_get(const int atom, const int coordinate) const
{
    return xQM[atom][coordinate];
}

real QMMM_QMrec_transfer::QMcharges_get(const int atom) const
{
    return QMcharges[atom];
}

void QMMM_QMrec_transfer::QMcharges_set(const int atom, const real value)
{
    QMcharges[atom] = value;
}

double QMMM_QMrec_transfer::pot_qmmm_get(const int atom) const
{
    return pot_qmmm[atom];
}

double QMMM_QMrec_transfer::pot_qmqm_get(const int atom) const
{
    return pot_qmqm[atom];
}

void QMMM_QMrec_transfer::pot_qmmm_set(const int atom, const double value)
{
    pot_qmmm[atom] = value;
}

void QMMM_QMrec_transfer::pot_qmqm_set(const int atom, const double value)
{
    pot_qmqm[atom] = value;
}

int QMMM_QMrec_transfer::atomicnumberQM_get(const int atom)const
{
    return atomicnumberQM[atom];
}

int QMMM_QMrec_transfer::QMcharge_get()const
{
    return QMcharge;
}

/*
int QMMM_QMrec_transfer::nelectrons_get()const
{
    return nelectrons;
}
*/

void QMMM_MMrec_transfer::init_MMrec(real       scalefactor_in,
                                     int        nrMMatoms_full_in,
                                     int        natoms,
                                     int        nrQMatoms,
                                     const int* indexQM,
                                     int*       found_mm_atoms)
{
    scalefactor    = scalefactor_in;
    nrMMatoms_full = nrMMatoms_full_in;
    indexMM_full.resize(nrMMatoms_full); // ???
    xMM_full.resizeWithPadding(nrMMatoms_full); // ???
    MMcharges_full.resize(nrMMatoms_full); // ???
    shiftMM_full.resize(nrMMatoms_full); // ???

    // fill the indexMM_full array
    *found_mm_atoms = 0;
    for (int i=0; i<natoms; i++)
    {
        bool is_mm_atom = true;
        for (int j=0; j<nrQMatoms; j++)
        {
            if (i == indexQM[j])
            {
                 is_mm_atom = false;
            }
        }
        if (is_mm_atom)
        {
            indexMM_full[*found_mm_atoms] = i;
	        (*found_mm_atoms)++;
        }
    }
}

// HOW ARE WE PASSING THE DEFINITION OF THE "QM GROUP" TO THIS PROCEDURE?
// PERHAPS PASS THE ct STRUCTURE AS AN ADDITIONAL ARGUMENT?
QMMM_rec_transfer::QMMM_rec_transfer(const t_commrec*  cr,
                                     const t_inputrec* ir,
                                     const t_forcerec* fr,
                                     const int         nAtoms,
                                     std::vector<int>& qmAtoms,
                                     std::vector<int>& atomicnumberQM_in,
                                     const int         qmmmVariant_in,
                                     const bool        dipCorrection,
                                     const int         phase,
                                     const int         iSite)
{
    // Put the atom numbers of atoms that belong to the QM group
    // into an array that will be copied later to QMMMrec->indexQM[..].
    // Also, it will be used to create an index array QMMMrec->bQMMM[],
    // which contains true/false for QM and MM (the other) atoms.

 // // issue a fatal if the user wants to run with more than one node
 // if (PAR(cr))
 // {
 //     gmx_fatal(FARGS, "QM/MM may not work in parallel due to neighborsearching issues, \
 //           use a single processor instead!\n");
 // }

    // The array bQMMM[] contains true/false for atoms that are QM/not QM.
    // We first set all elements at false.
    // Afterwards we use qm_arr (= MMrec->indexQM) to change
    // the elements corresponding to the QM atoms at true.

    pbcType = fr->pbcType;

 // std::vector<int> qmmmAtoms = qmmmAtomIndices(*ir, *mtop);

    qm = std::make_unique<QMMM_QMrec_transfer>();

    // store QM atoms in the QMrec and initialize
 // qm->init_QMrec(qmAtoms, QMcharge_in, atomicnumberQM_in, ir);
    qm->init_QMrec(qmAtoms, atomicnumberQM_in, ir);

    // print the current layer to allow users to check their input
    fprintf(stderr, "New site with nr of QM atoms %d\n", qm->nrQMatoms_get());

    // MM rec creation
    int nrMMatoms_full_in = nAtoms - (qm->nrQMatoms); // rest of the atoms
    int found_mm_atoms = 0;
    mm = std::make_unique<QMMM_MMrec_transfer>();
    mm->init_MMrec(ir->scalefactor, nrMMatoms_full_in, nAtoms, qm->nrQMatoms, qm->indexQM, &found_mm_atoms); 

    printf ("nAtoms) = %d\n(qr->qm->nrQMatoms) = %d\nmm->nrMMatoms_full = %d\n",
            nAtoms, qm->nrQMatoms_get(), mm->nrMMatoms_full);
    printf ("(found_mm_atoms) = %d\n", found_mm_atoms);

    // Look how the QM/MM electrostatics shall be treated.
    qm->qmmm_variant = qmmmVariant_in;
    switch (qm->qmmm_variant) {
	case eqmmmVACUO: // 0
	                fprintf(stdout, "No electrostatic QM/MM interaction.\nTo change, set environment variable GMX_QMMM_VARIANT.\n");
	                fprintf(stdout, "No electrostatic QM/MM interaction.\n");
	                break;
	case eqmmmPME: // 1
        {
	           if (pbcType != PbcType::Xyz)
               {
	               fprintf(stderr, "PME treatment of QM/MM electrostatics only possible with triclinic periodic system!\n");
	               exit(-1);
	           }
	           fprintf(stdout, "Electrostatic QM/MM interaction calculated with full PME treatment.\n");

               pme.resize(2);
               QMMM_PME_transfer& pme_full   = pme[0];
               QMMM_PME_transfer& pme_qmonly = pme[1];

               // PME data structure for the entire system
               pme_full.x.resizeWithPadding(qm->nrQMatoms + mm->nrMMatoms_full);
               pme_full.q.resize(qm->nrQMatoms + mm->nrMMatoms_full);
               pme_full.f.resizeWithPadding(qm->nrQMatoms + mm->nrMMatoms_full);
               snew(pme_full.pot, qm->nrQMatoms);
               
               // PME data structure for the QM-only system
               pme_qmonly.x.resizeWithPadding(qm->nrQMatoms);
               pme_qmonly.q.resize(qm->nrQMatoms);
               pme_qmonly.f.resizeWithPadding(qm->nrQMatoms);
               snew(pme_qmonly.pot, qm->nrQMatoms);
               
               if (dipCorrection)
               {
                   pme_full.surf_corr_pme   = true;
                   pme_full.epsilon_r       = qm->epsilon_r;
                   pme_qmonly.surf_corr_pme = true;
                   pme_qmonly.epsilon_r     = qm->epsilon_r;
			       fprintf(stdout, "Dipole (surface) correction for QM/MM PME applied ");
			       fprintf(stdout, "with a permittivity of %5.1f.\n", pme_qmonly.epsilon_r);
			       fprintf(stdout, "\nCurrently disabled due to solvent molecules broken across box boundary!\nExiting!\n\n");
                   exit(-1);
               }
               else
               {
                   pme_full.surf_corr_pme   = false;
                   pme_qmonly.surf_corr_pme = false;
			       fprintf(stdout, "No dipole (surface) correction for QM/MM PME, i.e. tin-foil boundary conditions.\n");
               }
	           break;
        }
	case eqmmmSWITCH: // 2
	             fprintf(stdout, "Electrostatic QM/MM interaction calculated with a switched cut-off.\n");
	             break;
	case eqmmmRFIELD: // 3
	             fprintf(stdout, "Electrostatic QM/MM interaction calculated with a reaction-field cut-off.\n");
	             break;
	case eqmmmSHIFT: // 4
	             fprintf(stdout, "Electrostatic QM/MM interaction calculated with a shifted cut-off.\n");
	             break;
	default:
	        fprintf(stderr, "Unrecognized choice for treatment of QM/MM electrostatics.\n");
	        fprintf(stderr, "Set environment variable GMX_QMMM_VARIANT to either 0, 1, 2, 3, or 4.\n");
	        exit(-1);
    }
    snew(qm->QMcharges, qm->nrQMatoms);

    init_dftbplus_transfer(this, ir->rcoulomb, ir->ewald_rtol, cr, phase, iSite);
} // init_QMMMrec_transfer


// Updates the shift and charges of *all of the* MM atoms in QMMMrec.
//   Only with DFTB.
//   (Not nice, should be done in a more elegant way...)
void QMMM_rec_transfer::update_QMMMrec_dftb(const t_commrec*  cr,
                                            const rvec*       shift_vec,
                                            const rvec        x[],
                                            const t_mdatoms*  md,
                                            const matrix      box)
{
    // INHERITED NOTE: is NOT yet working if there are no PBC.
    // Also in ns.c, simple NS needs to be fixed!
    //   As of 2019, ns.c does not exist any longer.

    // copy pointers
 // QMMM_QMrec& qm_ = qm[0]; // in case of normal QMMM, there is only one group
 // QMMM_MMrec& mm_ = mm[0];

    // init_pbc(box); needs to be called first, see pbc.h
    ivec null_ivec;
    clear_ivec(null_ivec);
    t_pbc pbc;
    set_pbc_dd(&pbc, pbcType, DOMAINDECOMP(cr) ? cr->dd->numCells : null_ivec, false, box);

 // printf("There are %d QM atoms, namely:", qm_.nrQMatoms);
 // for (int i=0; i<qm_.nrQMatoms; i++)
 //   printf(" %d", qm_.indexQM[i]);
 // printf("\n");

    // Compute the shift for the MM atoms with respect to QM atom [0].
    // TODO: This looks like a viable first guess, but is that correct?
    // Related to the problem of contributions to virial pressure
    //   in a system treated with particle--mesh Ewald.
    rvec crd;
    rvec_sub(x[qm->indexQM[0]], shift_vec[qm->shiftQM[0]], crd);
    for (int i=0; i<mm->nrMMatoms_full; i++) {
        rvec dx;
        mm->shiftMM_full[i] = pbc_dx_aiuc(&pbc, crd, x[mm->indexMM_full[i]], dx);
    }

 // // previous version of the loop
 // for (i=0; i<mm_.nrMMatoms; i++) {
 //     ivec dx;
 //     current_shift = pbc_dx_aiuc(&pbc, x[qm_.indexQM[0]], x[mm_.indexMM[i]], dx);
 //     crd[0] = IS2X(QMMMlist->shift[i]) + IS2X(qm_i_particles[i].shift);
 //     crd[1] = IS2Y(QMMMlist->shift[i]) + IS2Y(qm_i_particles[i].shift);
 //     crd[2] = IS2Z(QMMMlist->shift[i]) + IS2Z(qm_i_particles[i].shift);
 //     is     = static_cast<int>(XYZ2IS(crd[0], crd[1], crd[2]));
 //     mm_.shiftMM[i] = is;
 // }

    for (int i = 0; i < mm->nrMMatoms_full; i++) // no free energy yet
    {
        mm->MMcharges_full[i] = md->chargeA[mm->indexMM_full[i]] * mm->scalefactor;
    }
} // update_QMMMrec_dftb


// create the SR MM list using the Verlet neighborlist
void QMMM_rec_transfer::update_QMMMrec_verlet_ns(const t_commrec*          cr,
                                                 const nonbonded_verlet_t* nbv,
                                                 const rvec                x[],
                                                 const t_mdatoms*          md,
                                                 const matrix              box)
{
 //  * COMMENTS TO THE FORMER GROUP-SCHEME BASED VERSION OF THIS FUNCTION:
 //  *********************************************************************
 //  * Create/update a number of QMMMrec entries:
 //  * 1) shiftQM -- shifts of the QM atoms
 //  * 2) indexMM -- indices of the MM atoms
 //  * 3) shiftMM -- shifts of the MM atoms
 //  * 4) shifted coordinates of the MM atoms
 //  *       --- NOT THIS ONE, BECAUSE A SEPARATE ROUTINE IS USED!
 //  * (The shifts are used to compute the virial of the QM/MM particles.)
 //  *

 //  * if atom i shall be considered as MM,
 //  *   isMMatom[i] = true
 //  * UPDATE: store the shift in this array,
 //  * and change 'bool' to 'int':
 //  *   shiftMMatom[i] = the value of shift

    std::vector<int> shiftMMatom(md->nr, -1); // ALL ATOMS IN SIMULATION - IS THAT NECESSARY???

    // init PBC
    ivec null_ivec;
    clear_ivec(null_ivec);
    t_pbc pbc;
    set_pbc_dd(&pbc, pbcType, DOMAINDECOMP(cr) ? cr->dd->numCells : null_ivec, false, box);

    // copy pointers
 // QMMM_QMrec&                           qm_  = qm[0];
 // QMMM_MMrec&                           mm_  = mm[0];
    gmx::ArrayRef<const NbnxnPairlistCpu> nbl = nbv->pairlistSets().pairlistSet(gmx::InteractionLocality::Local).cpuLists();
    int                                   nnbl = nbl.ssize();
    const gmx::ArrayRef<const int>        atomIndices = nbv->pairSearch_->gridSet().atomIndices();

    // QM shift array
    // !!! CHECK THIS !!!
    rvec dx;
    qm->shiftQM[0] = XYZ2IS(0, 0, 0);
    for (int i = 1; i < qm->nrQMatoms; i++)
    {
        qm->shiftQM[i] = pbc_dx_aiuc(&pbc, x[qm->indexQM[0]], x[qm->indexQM[i]], dx);
    }
 // for (int i = 0; i < qm->nrQMatoms; i++)
 // {
 //     printf("VERLET QM SHIFT [%d] = %d\n", i, qm->shiftQM[i]);
 // }

    // LOOP OVER THE nnbl NEIGHBORLISTS!
    //   THIS IS NECESSARY WITH MULTITHREADING
    for (int inbl=0; inbl<nnbl; inbl++)
    {
        // loop over CI clusters
        for (unsigned ci=0; ci<nbl[inbl].ci.size(); ci++)
	    {
            // is there a QM atom in this CI cluster?
	        bool qm_atom_in_ci = false;
	        // break the loop if a QM atom has already been found
	        for (int ii=0; ii<nbl[inbl].na_ci && !qm_atom_in_ci; ii++)
	        {
	            // compare to indices of QM atoms
	            for (int iq=0; iq<qm->nrQMatoms && !qm_atom_in_ci; iq++)
		        {
                    const int iIndex = nbl[inbl].na_ci * nbl[inbl].ci[ci].ci + ii;
                    const int iAtom  = atomIndices[iIndex];
                    //  FORMERLY:
		            // const int iAtom  = nbs->a[nbl[inbl].na_ci * nbl[inbl].ci[ci].ci + ii];
		            if (qm->indexQM[iq] == iAtom)
		            {
		                qm_atom_in_ci = true;
		            }
		        }
	        }
	        // get the shift of this CI cluster */
	     // shift = nbl[inbl]->ci[ci].shift & NBNXN_CI_SHIFT;

            // loop over the corresponding CJ clusters
	        for (int cj = nbl[inbl].ci[ci].cj_ind_start; cj < nbl[inbl].ci[ci].cj_ind_end; cj++)
	        {
	            // is there a QM atom in this CJ cluster?
	            bool qm_atom_in_cj = false;
	            // break the loop if a QM atom has already been found
	            for (int jj=0; jj<nbl[inbl].na_cj && !qm_atom_in_cj; jj++)
		        {
	                // compare to indices of QM atoms
	                for (int jq=0; jq<qm->nrQMatoms && !qm_atom_in_cj; jq++)
		            {
                        const int iIndex = nbl[inbl].na_cj * nbl[inbl].cj[cj].cj + jj;
                        const int iAtom  = atomIndices[iIndex];
                        //  FORMERLY:
		                // if (qm->indexQM[jq] == nbs->a[nbl[inbl].na_cj * nbl[inbl].cj[cj].cj + jj])
		                if (qm->indexQM[jq] == iAtom)
			            {
			                qm_atom_in_cj = true;
		                }
		            }
		        }

                // if there is a QM atom in cluster CI,
		        //   then put the non-QM atoms in cluster CJ into the MM list
	            if (qm_atom_in_ci)
		        {
	                put_cluster_in_MMlist_verlet(nbl[inbl].cj[cj].cj, nbl[inbl].na_cj,
		                                        qm->nrQMatoms, qm->indexQM, atomIndices, shiftMMatom.data(), &pbc, x);
	            }

                // if there is a QM atom in cluster CJ,
		        //   then put the non-QM atoms in cluster CI into the MM list
	            if (qm_atom_in_cj)
		        {
	                put_cluster_in_MMlist_verlet(nbl[inbl].ci[ci].ci, nbl[inbl].na_ci,
		                                        qm->nrQMatoms, qm->indexQM, atomIndices, shiftMMatom.data(), &pbc, x);
	            }
	        }
	    }
    }

    // count the MM atoms found in the above search
    int nrMMatoms = 0;
    for (int i=0; i<md->nr; i++) {
        // criterium for MM atom found
        if (shiftMMatom[i] != -1)
	    {
	        nrMMatoms++;
	    }
    }

 // printf("Number of potential MM atoms as found in the Verlet neighbor lists: %d\n", nrMMatoms);

    // allocate space and fill the array with atom numbers
    mm->nrMMatoms_nbl = nrMMatoms;
    mm->indexMM_nbl.resize(nrMMatoms);
    mm->shiftMM_nbl.resize(nrMMatoms);
    // index i runs along isMMatom / shiftMMatom,
    //       j runs along the new indexMM_nbl array

    int count=0;
    for (int atom=0; atom<md->nr; atom++) {
        // criterium for MM atom found
        if (shiftMMatom[atom] != -1)
	    {
	        mm->indexMM_nbl[count] = atom;
	        mm->shiftMM_nbl[count] = shiftMMatom[atom];
	        count++;
	    }
    }
    // check
    if (count != mm->nrMMatoms_nbl) {
        printf("ERROR IN MM ATOM SEARCH -- VERLET BASED SCHEME\n");
	    exit(-1);
    }
} // update_QMMMrec_verlet_ns
/*
void QMMM_rec_transfer::calculate_QMMM(const t_commrec*      cr,
                                       gmx::ForceWithVirial* forceWithVirial,
                                       t_nrnb*               nrnb,
                                       gmx_wallcycle_t       wcycle)
{
 // QMMM_QMrec_transfer* qm_ = qm;

    gmx::RVec *fMM      = forceWithVirial->force_.data();

 // std::vector<rvec> forces;
 // forces.resize(qm->nrQMatoms);
    std::vector<rvec> forces(qm->nrQMatoms);

    call_dftbplus_transfer(this, cr, forces.data(), nrnb, wcycle);

    for (int i = 0; i < qm->nrQMatoms; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            fMM[qm->indexQM[i]][j]        -= forces.data()[i][j];
        }
     // printf("F[%5d] = %8.2f %8.2f %8.2f\n", qm_->indexQM[i], forces[i][0], forces[i][1], forces[i][2]);
    }

    // NO MM FORCES ARE NEEDED -- RIGHT?
} // calculate_QMMM
*/
