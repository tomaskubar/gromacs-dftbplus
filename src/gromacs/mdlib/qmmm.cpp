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
#include "gromacs/mdlib/qm_dftbplus.h"
#include "gromacs/mdlib/qm_gamess.h"
#include "gromacs/mdlib/qm_gaussian.h"
#include "gromacs/mdlib/qm_mopac.h"
#include "gromacs/mdlib/qm_orca.h"
#include "gromacs/mdlib/qm_nn.h"
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

// When not built in a configuration with QMMM support, much of this
// code is unreachable by design. Tell clang not to warn about it.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunreachable-code"
#pragma GCC diagnostic ignored "-Wmissing-noreturn"

void put_cluster_in_MMlist_verlet(int                            ck, // cluster number
                                  int                            na_ck, // # of atoms in cluster
                                  int                            nrQMatoms,
                                  const int*                     indexQM,
                                        gmx::ArrayRef<const int> atomIndices,
			                      int*                           shiftMMatom,
                                  // ^ also has a role of "bool* isMMatom"
				                  t_pbc*                         pbc,
				                  const rvec*                    x);

/*
std::unique_ptr<QMMM_rec>
void init_QMMM_rec(const t_commrec  *cr,
              const gmx_mtop_t *mtop,
              const t_inputrec *ir,
              const t_forcerec *fr,
              const gmx_wallcycle_t gmx_unused wcycle)
{
    return std::make_unique<QMMM_rec>(cr, mtop, ir, fr, wcycle);
}
*/

static real call_QMroutine(const t_commrec*  cr,
                           QMMM_rec*         qr,
                           QMMM_QMrec*       qm,
                           QMMM_MMrec*       mm,
                           rvec              f[],
                           rvec              fshift[],
                           t_nrnb*           nrnb,
                           gmx_wallcycle_t   wcycle)
{
    // Makes a call to the requested QM routine (qm->QMmethod).
    // Note that f is actually the gradient, i.e. -f

    if (GMX_QMMM_MOPAC)
    {
        return call_mopac(*qm, *mm, f, fshift);
    }
    else if (GMX_QMMM_GAMESS)
    {
        return call_gamess(*qm, *mm, f, fshift);
    }
    else if (GMX_QMMM_GAUSSIAN)
    {
        return qm->gaussian.call_gaussian(*qm, *mm, f, fshift);
    }
    else if (GMX_QMMM_ORCA)
    {
        return call_orca(*qm, *mm, f, fshift);
    }
    else if (GMX_QMMM_DFTBPLUS)
    {
        return call_dftbplus(qr, cr, qm, *mm, f, fshift, nrnb, wcycle);
    }
    else if (GMX_QMMM_NN)
    {
        return call_nn(qr, cr, qm, *mm, f, fshift, nrnb, wcycle);
    }
    else
    {
        gmx_fatal(FARGS, "Unknown QM software -- should never land here :-/");
    }
}

// Update QM and MM coordinates in the QM/MM data structures.
// New version of the function:
// Update the coordinates of the MM atoms on the short-range neighborlist!
// The NBlist needs to have been created previously by either group or Verlet scheme.
void QMMM_rec::update_QMMM_coord(const t_commrec*  cr,
                                 rvec*             shift_vec,
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

    QMMM_QMrec& qm_ = qm[0];
    QMMM_MMrec& mm_ = mm[0];
    real rcut = qm_.rcoulomb > 0.1 ? qm_.rcoulomb : 999999.; // infinity
    std::vector<bool> isCurrentMMatom;
    isCurrentMMatom.resize(mm_.nrMMatoms_nbl);

    // shift the QM atoms into the central box
    for (int i = 0; i < qm_.nrQMatoms; i++)
    {
        rvec_sub(x[qm_.indexQM[i]], shift_vec[qm_.shiftQM[i]], qm_.xQM[i]);
    }

    // copy box size
    copy_mat(box, qm_.box);

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
    for (int i = 0; i < mm_.nrMMatoms_nbl; i++)
    {
	    isCurrentMMatom[i] = false;
	 // printf("DEBUG_MM TEST %5d %5d %2d", i, mm_.indexMM_nbl[i], mm_.shiftMM_nbl[i]);
	    // loop over all QM atoms here
	    for (int q=0; q<qm_.nrQMatoms; q++)
	    {
            rvec bond;
            pbc_dx_aiuc(&pbc, x[qm_.indexQM[q]], x[mm_.indexMM_nbl[i]], bond);
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
    mm_.nrMMatoms = nrMMatoms;
    mm_.indexMM.resize(nrMMatoms);
    mm_.MMcharges.resize(nrMMatoms);
    mm_.shiftMM.resize(nrMMatoms);
    mm_.xMM.resizeWithPadding(nrMMatoms);

    int index = 0; // runs over the identified MM atoms
    for (int i = 0; i < mm_.nrMMatoms_nbl; i++)
    {
	    if (isCurrentMMatom[i])
	    {
	        // Add to list!
	        mm_.indexMM[index] = mm_.indexMM_nbl[i];

	        // Also add charge
	        mm_.MMcharges[index] = md->chargeA[mm_.indexMM[index]] * mm_.scalefactor;

            // Having obtained the shift at NS time (update_qmmmrec),
            //   merely copy it here to shiftMM[]
	        mm_.shiftMM[index] = mm_.shiftMM_nbl[i];

	        // one MM atom found => increment index */
	        index++;
	    }
    }

    // also shift the MM atoms into the central box

 //   for (int a=0; a<45; a++)
 //     printf("SHIFT %2d: %7.3f %7.3f %7.3f\n", a,
 //     fr->shift_vec[a][XX], fr->shift_vec[a][YY], fr->shift_vec[a][ZZ]);

    for (int ind = 0; ind < mm_.nrMMatoms; ind++)
    {
        rvec_sub(x[mm_.indexMM[ind]], shift_vec[mm_.shiftMM[ind]], mm_.xMM[ind]);
 //     printf("COORD MM %4d %2d\n", mm_.indexMM[ind], mm_.shiftMM[ind]);
    }

    // For DFTB, also update the coordinates of *all* of the MM atoms,
    //   not only those on the short-range neighborlist.
    // Do not shift the MM atoms into the central box!
    //   It might break the calculation of the surface correction in the Ewald sum.
    if (GMX_QMMM_DFTBPLUS)
    {
        for (int i = 0; i < mm_.nrMMatoms_full; i++)
        {
            copy_rvec(x[mm_.indexMM_full[i]], mm_.xMM_full[i]);
        }
    }
} // update_QMMM_coord

void QMMM_QMrec::init_QMrec(int               grpnr,
                            int               nr,
                            const int*        atomarray,
                            const gmx_mtop_t* mtop,
                            const t_inputrec* ir)
{
    nrQMatoms = nr;
    snew(xQM, nr);
    snew(indexQM, nr);
    snew(shiftQM, nr);
    for (int i = 0; i < nr; i++)
    {
        indexQM[i] = atomarray[i];
    }

    snew(atomicnumberQM, nr);
    int molb = 0;
    for (int i = 0; i < nrQMatoms; i++)
    {
        const t_atom &atom = mtopGetAtomParameters(mtop, indexQM[i], &molb);
        nelectrons        += mtop->atomtypes.atomnumber[atom.type];
        atomicnumberQM[i]  = mtop->atomtypes.atomnumber[atom.type];
    }

    QMcharge      = ir->opts.QMcharge[grpnr];
    multiplicity  = ir->opts.QMmult[grpnr];
    nelectrons   -= ir->opts.QMcharge[grpnr];

    QMmethod      = ir->opts.QMmethod[grpnr];
    QMbasis       = ir->opts.QMbasis[grpnr];

    // hack to prevent gaussian from reinitializing all the time
    gaussian.nQMcpus = 0; // number of CPU's to be used by g01, is set
                          // upon initializing gaussian with init_gaussian()

    rcoulomb      = ir->rcoulomb;
    ewaldcoeff_q  = calc_ewaldcoeff_q(ir->rcoulomb, ir->ewald_rtol);
    epsilon_r     = ir->epsilon_r;

    snew(pot_qmmm, nr);
    snew(pot_qmqm, nr);

    significant_structure = false; // to mark significant deviations
    
} // init_QMrec

int QMMM_QMrec::nrQMatoms_get() const
{
    return nrQMatoms;
}

int QMMM_QMrec::qmmm_variant_get() const
{
    return qmmm_variant;
}

double QMMM_QMrec::xQM_get(const int atom, const int coordinate) const
{
    return xQM[atom][coordinate];
}

real QMMM_QMrec::QMcharges_get(const int atom) const
{
    return QMcharges[atom];
}

void QMMM_QMrec::QMcharges_set(const int atom, const real value)
{
    QMcharges[atom] = value;
}

double QMMM_QMrec::pot_qmmm_get(const int atom) const
{
    return pot_qmmm[atom];
}

double QMMM_QMrec::pot_qmqm_get(const int atom) const
{
    return pot_qmqm[atom];
}

void QMMM_QMrec::pot_qmmm_set(const int atom, const double value)
{
    pot_qmmm[atom] = value;
}

void QMMM_QMrec::pot_qmqm_set(const int atom, const double value)
{
    pot_qmqm[atom] = value;
}

int QMMM_QMrec::atomicnumberQM_get(const int atom)const
{
    return atomicnumberQM[atom];
}

int QMMM_QMrec::QMcharge_get()const
{
    return QMcharge;
}

int QMMM_QMrec::multiplicity_get()const
{
    return multiplicity;
}

int QMMM_QMrec::QMmethod_get()const
{
    return QMmethod;
}

int QMMM_QMrec::QMbasis_get()const
{
    return QMbasis;
}

int QMMM_QMrec::nelectrons_get()const
{
    return nelectrons;
}

int QMMM_QMrec::CASelectrons_get()const
{
    return CASelectrons;
}

int QMMM_QMrec::CASorbitals_get()const
{
    return CASorbitals;
}

void QMMM_MMrec::init_MMrec(real scalefactor_in,
                            int  nrMMatoms_full_in,
                            int  natoms,
                            int  nrQMatoms,
                            const int* indexQM,
                            int* found_mm_atoms)
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

QMMM_rec::QMMM_rec(const t_commrec*                 cr,
                   const gmx_mtop_t*                mtop,
                   const t_inputrec*                ir,
                   const t_forcerec*                fr)
 //                const gmx_wallcycle_t gmx_unused wcycle)
{
#if GMX_QMMM
    // Put the atom numbers of atoms that belong to the QMMM group
    // into an array that will be copied later to QMMMrec->indexQM[..].
    // Also, it will be used to create an index array QMMMrec->bQMMM[],
    // which contains true/false for QM and MM (the other) atoms.

    if (!GMX_QMMM)
    {
        gmx_incons("Compiled without QMMM");
    }

    // issue a fatal if the user wants to run with more than one node
    if (PAR(cr))
    {
        gmx_fatal(FARGS, "QM/MM may not work in parallel due to neighborsearching issues, \
              use a single processor instead!\n");
    }

    // The array bQMMM[] contains true/false for atoms that are QM/not QM.
    // We first set all elements at false.
    // Afterwards we use qm_arr (= MMrec->indexQM) to change
    // the elements corresponding to the QM atoms at true.

    // We take the possibility into account
    // that a user has defined more than one QM group:
    // HOW SHOULD WE PROCEED IN THAT CASE?
    // IT WOULD BE COOL TO BE ABLE TO DO IT!

    // An ugly work-around in case there is only one group.
    // In this case, the whole system is treated as QM.
    // Otherwise, the second group is always the rest of the total system
    //   and is treated as MM.

    // Small problem if there is only QM... so no MM. */
    
    pbcType = fr->pbcType;

    int numQmmmGroups = ir->opts.ngQM;

    if (numQmmmGroups > 1) {
        fprintf(stderr, "\nQM/MM cannot calculate more than 1 group of atoms at the moment\nExiting!\n\n");
        exit(-1);
    }

    // There are numQmmmGroups groups of QM atoms.
    // Previously, multiple QM groups typically meant
    // that the user wanted to do ONIOM.
    // However, maybe it should also be possible to define
    // more than one QM subsystem with independent neighbourlists.
    // Gerrit Groenhof said he would have to think about that...
    // (11-11-2003)

    std::vector<int> qmmmAtoms = qmmmAtomIndices(*ir, *mtop);

    qm.resize(numQmmmGroups);

    // Standard QMMM (no ONIOM).
    // All layers are merged together, so there is one QM subsystem and one MM subsystem.
    // Also, we set the charges to zero in mtop
    //   to prevent the innerloops from doubly counting the electrostatic QM--MM interaction.
    // TODO: Consider doing this in grompp instead.

    // store QM atoms in the QMrec and initialise
    qm[0].init_QMrec(0, qmmmAtoms.size(), qmmmAtoms.data(), mtop, ir);

    // print the current layer to allow users to check their input
    fprintf(stderr, "Layer %d\nnr of QM atoms %d\n", 0, qm[0].nrQMatoms);
    fprintf(stderr, "QMlevel: %s/%s\n\n",
            eQMmethod_names[qm[0].QMmethod], eQMbasis_names[qm[0].QMbasis]);

    // MM rec creation
    int nrMMatoms_full_in = (mtop->natoms)-(qm[0].nrQMatoms); // rest of the atoms
    int found_mm_atoms = 0;
    mm.resize(1);
    QMMM_MMrec& mm_ = mm[0];
    mm_.init_MMrec(ir->scalefactor, nrMMatoms_full_in, mtop->natoms, qm[0].nrQMatoms, qm[0].indexQM, &found_mm_atoms); 

    printf ("(mtop->natoms) = %d\n(qr->qm[0]->nrQMatoms) = %d\nmm->nrMMatoms_full = %d\n",
            (mtop->natoms), (qm[0].nrQMatoms), mm_.nrMMatoms_full);
    printf ("(found_mm_atoms) = %d\n", found_mm_atoms);

    // these variables get updated in the update QMMMrec // ???

    // OLD COMMENT but maybe useful in the future:
    //   With only one layer there is only one initialization needed.
    //   Multilayer is a bit more complicated as it requires
    //   a re-initialization at every step of the simulation.
    //   This is due to the use of COMMON blocks in Fortran QM subroutines.

    if (GMX_QMMM_MOPAC)
    {
        init_mopac(qm[0]);
    }
    else if (GMX_QMMM_GAMESS)
    {
        init_gamess(cr, qm[0], mm_);
    }
    else if (GMX_QMMM_GAUSSIAN)
    {
        qm[0].gaussian.init_gaussian();
    }
    else if (GMX_QMMM_ORCA)
    {
        init_orca(&(qm[0]));
    }
    else if (GMX_QMMM_DFTBPLUS)
    {
        // Look how the QM/MM electrostatics shall be treated.
        // In the future, this could be performed for QM/MM in general,
        //   not only with DFTB+.
        char *env1 = getenv("GMX_QMMM_VARIANT");
        char *env2 = getenv("GMX_QMMM_PME_DIPCOR");
        if (env1 == nullptr)
        {
            qm[0].qmmm_variant = eqmmmVACUO;
		    fprintf(stdout, "No electrostatic QM/MM interaction.\nTo change, set environment variable GMX_QMMM_VARIANT.\n");
        }
        else
        {
            sscanf(env1, "%d", &(qm[0].qmmm_variant));
            switch (qm[0].qmmm_variant) {
		    case eqmmmVACUO: // 0
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
                       pmedata              = &(fr->pmedata);
                       QMMM_PME& pme_full   = pme[0];
                       QMMM_PME& pme_qmonly = pme[1];

                       // PME data structure for the entire system
                       pme_full.x.resizeWithPadding(qm[0].nrQMatoms + mm_.nrMMatoms_full);
                       pme_full.q.resize(qm[0].nrQMatoms + mm_.nrMMatoms_full);
                       pme_full.f.resizeWithPadding(qm[0].nrQMatoms + mm_.nrMMatoms_full);
                       snew(pme_full.pot, qm[0].nrQMatoms);
                       
                       // PME data structure for the QM-only system
                       pme_qmonly.x.resizeWithPadding(qm[0].nrQMatoms);
                       pme_qmonly.q.resize(qm[0].nrQMatoms);
                       pme_qmonly.f.resizeWithPadding(qm[0].nrQMatoms);
                       snew(pme_qmonly.pot, qm[0].nrQMatoms);
                       
                       if (env2 != nullptr)
                       {
                           pme_full.surf_corr_pme   = true;
                           pme_full.epsilon_r       = qm[0].epsilon_r;
                           pme_qmonly.surf_corr_pme = true;
                           pme_qmonly.epsilon_r     = qm[0].epsilon_r;
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
        }
        snew(qm[0].QMcharges, qm[0].nrQMatoms);

        init_dftbplus(&(qm[0]), this, ir, cr); //, wcycle);
    }
    else if (GMX_QMMM_NN)
    {
        // Look how the QM/MM electrostatics shall be treated.
        // In the future, this could be performed for QM/MM in general,
        //   not only with DFTB+.
        char *env1 = getenv("GMX_QMMM_VARIANT");
        char *env2 = getenv("GMX_QMMM_PME_DIPCOR");
        if (env1 == nullptr)
        {
            qm[0].qmmm_variant = eqmmmVACUO;
		    fprintf(stdout, "No electrostatic QM/MM interaction.\nTo change, set environment variable GMX_QMMM_VARIANT.\n");
        }
        else
        {
            sscanf(env1, "%d", &(qm[0].qmmm_variant));
            switch (qm[0].qmmm_variant) {
		    case eqmmmVACUO: // 0
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
                       pmedata              = &(fr->pmedata);
                       QMMM_PME& pme_full   = pme[0];
                       QMMM_PME& pme_qmonly = pme[1];

                       // PME data structure for the entire system
                       pme_full.x.resizeWithPadding(qm[0].nrQMatoms + mm_.nrMMatoms_full);
                       pme_full.q.resize(qm[0].nrQMatoms + mm_.nrMMatoms_full);
                       pme_full.f.resizeWithPadding(qm[0].nrQMatoms + mm_.nrMMatoms_full);
                       snew(pme_full.pot, qm[0].nrQMatoms);
                       
                       // PME data structure for the QM-only system
                       pme_qmonly.x.resizeWithPadding(qm[0].nrQMatoms);
                       pme_qmonly.q.resize(qm[0].nrQMatoms);
                       pme_qmonly.f.resizeWithPadding(qm[0].nrQMatoms);
                       snew(pme_qmonly.pot, qm[0].nrQMatoms);
                       
                       if (env2 != nullptr)
                       {
                           pme_full.surf_corr_pme   = true;
                           pme_full.epsilon_r       = qm[0].epsilon_r;
                           pme_qmonly.surf_corr_pme = true;
                           pme_qmonly.epsilon_r     = qm[0].epsilon_r;
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
            case eqmmmNN: // 5
		                 fprintf(stdout, "Electrostatic QM/MM interaction calculated with the whole box.\n");
		                 break;
		    default:
		            fprintf(stderr, "Unrecognized choice for treatment of QM/MM electrostatics.\n");
		            fprintf(stderr, "Set environment variable GMX_QMMM_VARIANT to either 0, 1, 2, 3, or 4.\n");
	                exit(-1);
		    }
        }
        snew(qm[0].QMcharges, qm[0].nrQMatoms);

        init_nn(&qm[0]);
    }
    else
    {
        gmx_fatal(FARGS, "Unknown QM software -- should never land here :-/");
    }
#else // GMX_QMMM
    gmx_incons("Compiled without QMMM");
    (void) cr;
    (void) mtop;
    (void) ir;
    (void) fr;
#endif
} // init_QMMMrec

std::vector<int> qmmmAtomIndices(const t_inputrec& ir, const gmx_mtop_t& mtop)
{
    const int               numQmmmGroups = ir.opts.ngQM;
    const SimulationGroups& groups        = mtop.groups;
    std::vector<int>        qmmmAtoms;
    for (int i = 0; i < numQmmmGroups; i++)
    {
        for (const AtomProxy atomP : AtomRange(mtop))
        {
            int index = atomP.globalAtomNumber();
            if (getGroupType(groups, SimulationAtomGroupType::QuantumMechanics, index) == i)
            {
                qmmmAtoms.push_back(index);
            }
        }
    }
    return qmmmAtoms;
}

void removeQmmmAtomCharges(gmx_mtop_t* mtop, gmx::ArrayRef<const int> qmmmAtoms)
{
    int molb = 0;
    for (gmx::index i = 0; i < qmmmAtoms.ssize(); i++)
    {
        int indexInMolecule;
        mtopGetMolblockIndex(mtop, qmmmAtoms[i], &molb, nullptr, &indexInMolecule);
        t_atom* atom = &mtop->moltype[mtop->molblock[molb].type].atoms.atom[indexInMolecule];
        atom->q      = 0.0;
        atom->qB     = 0.0;
    }
}

// Updates the shift and charges of *all of the* MM atoms in QMMMrec.
//   Only with DFTB.
//   (Not nice, should be done in a more elegant way...)
void QMMM_rec::update_QMMMrec_dftb(const t_commrec*  cr,
                                   rvec*             shift_vec,
                                   const rvec        x[],
                                   const t_mdatoms*  md,
                                   const matrix      box)
{
    // INHERITED NOTE: is NOT yet working if there are no PBC.
    // Also in ns.c, simple NS needs to be fixed!
    //   As of 2019, ns.c does not exist any longer.

    // copy pointers
    QMMM_QMrec& qm_ = qm[0]; // in case of normal QMMM, there is only one group
    QMMM_MMrec& mm_ = mm[0];

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
    rvec_sub(x[qm_.indexQM[0]], shift_vec[qm_.shiftQM[0]], crd);
    for (int i=0; i<mm_.nrMMatoms_full; i++) {
        rvec dx;
        mm_.shiftMM_full[i] = pbc_dx_aiuc(&pbc, crd, x[mm_.indexMM_full[i]], dx);
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

    for (int i = 0; i < mm_.nrMMatoms_full; i++) // no free energy yet
    {
        mm_.MMcharges_full[i] = md->chargeA[mm_.indexMM_full[i]] * mm_.scalefactor;
    }
} // update_QMMMrec_dftb

// ADD THE NON-QM ATOMS IN THE VERLET CLUSTER ck TO THE LIST OF MM ATOMS
void put_cluster_in_MMlist_verlet(int                            ck, // cluster number
                                  int                            na_ck, // # of atoms in cluster
                                  int                            nrQMatoms,
                                  const int*                     indexQM,
                                  const gmx::ArrayRef<const int> atomIndices,
			                      int*                           shiftMMatom,
                                  // ^ also has a role of "bool* isMMatom"
				                  t_pbc*                         pbc,
				                  const rvec*                    x)
{
 //  * This calculation of shift would be desirable,
 //  * but it does not seem to work properly!
 // ivec crd;
 // crd[XX] = (is_j_cluster ? 1 : -1) * IS2X(shift) + IS2X(qm->shiftQM[qm_atom]);
 // crd[YY] = (is_j_cluster ? 1 : -1) * IS2Y(shift) + IS2Y(qm->shiftQM[qm_atom]);
 // crd[ZZ] = (is_j_cluster ? 1 : -1) * IS2Z(shift) + IS2Z(qm->shiftQM[qm_atom]);
 // int is = IVEC2IS(crd);

    // Loop over the atoms in the cluster ck.
    for (int k=0; k<na_ck; k++)  // NA_CK IS USUALLY 4 (SIMD RELATED)
    {
	    int ck_atom = atomIndices[na_ck * ck + k];
	    if (ck_atom < 0)
	    {
	        // The value of -1 in the Verlet list is for padding purpose only.
	        // It does not correspond to any atom.
	        // Therefore, ignore!
	        continue;
	    }
	    // In the following loop, determine 2 things:
	    // 1: the shift to put the k-th atom (a putative MM atom)
	    //    shortest-distance with respect to the nearest QM atom
	    // 2: whether the k-th atom is a QM atom
	    real dist = 1000.;
	    int sh = -1;
	    bool is_qmatom = false;
	    for (int q=0; q<nrQMatoms; q++)
	    {
	        // 1: the shift -- this calculation looks OK!
	        rvec bond;
	        int sh_t = pbc_dx_aiuc(pbc, x[indexQM[q]], x[ck_atom], bond);
	        if (norm(bond) < dist)
	        {
	            dist = norm(bond);
		        sh = sh_t;
	        }
	        // 2: a QM atom?
	        if (ck_atom == indexQM[q])
	        {
	            is_qmatom = true;
	        }
	    }
	    // If it is not a QM atom, then put it in the list and store the shift.
	    if (!is_qmatom)
	    {
	        // printf("FOUND_MM_ATOM %5d in cluster %4d\n", ck_atom, ck);
	        shiftMMatom[ck_atom] = sh; // true;
	    }
    }
}

// create the SR MM list using the Verlet neighborlist
void QMMM_rec::update_QMMMrec_verlet_ns(const t_commrec*    cr,
                                        nonbonded_verlet_t* nbv,
                                        const rvec          x[],
                                        const t_mdatoms*    md,
                                        const matrix        box)
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
    QMMM_QMrec&                           qm_  = qm[0];
    QMMM_MMrec&                           mm_  = mm[0];
    gmx::ArrayRef<const NbnxnPairlistCpu> nbl = nbv->pairlistSets().pairlistSet(gmx::InteractionLocality::Local).cpuLists();
    int                                   nnbl = nbl.ssize();
    const gmx::ArrayRef<const int>        atomIndices = nbv->pairSearch_->gridSet().atomIndices();

    // QM shift array
    // !!! CHECK THIS !!!
    rvec dx;
    qm_.shiftQM[0] = XYZ2IS(0, 0, 0);
    for (int i = 1; i < qm_.nrQMatoms; i++)
    {
        qm_.shiftQM[i] = pbc_dx_aiuc(&pbc, x[qm_.indexQM[0]], x[qm_.indexQM[i]], dx);
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
	            for (int iq=0; iq<qm_.nrQMatoms && !qm_atom_in_ci; iq++)
		        {
                    const int iIndex = nbl[inbl].na_ci * nbl[inbl].ci[ci].ci + ii;
                    const int iAtom  = atomIndices[iIndex];
                    //  FORMERLY:
		            // const int iAtom  = nbs->a[nbl[inbl].na_ci * nbl[inbl].ci[ci].ci + ii];
		            if (qm_.indexQM[iq] == iAtom)
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
	                for (int jq=0; jq<qm_.nrQMatoms && !qm_atom_in_cj; jq++)
		            {
                        const int iIndex = nbl[inbl].na_cj * nbl[inbl].cj[cj].cj + jj;
                        const int iAtom  = atomIndices[iIndex];
                        //  FORMERLY:
		                // if (qm->indexQM[jq] == nbs->a[nbl[inbl].na_cj * nbl[inbl].cj[cj].cj + jj])
		                if (qm_.indexQM[jq] == iAtom)
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
		                                        qm_.nrQMatoms, qm_.indexQM, atomIndices, shiftMMatom.data(), &pbc, x);
	            }

                // if there is a QM atom in cluster CJ,
		        //   then put the non-QM atoms in cluster CI into the MM list
	            if (qm_atom_in_cj)
		        {
	                put_cluster_in_MMlist_verlet(nbl[inbl].ci[ci].ci, nbl[inbl].na_ci,
		                                        qm_.nrQMatoms, qm_.indexQM, atomIndices, shiftMMatom.data(), &pbc, x);
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
    mm_.nrMMatoms_nbl = nrMMatoms;
    mm_.indexMM_nbl.resize(nrMMatoms);
    mm_.shiftMM_nbl.resize(nrMMatoms);
    // index i runs along isMMatom / shiftMMatom,
    //       j runs along the new indexMM_nbl array

    int count=0;
    for (int atom=0; atom<md->nr; atom++) {
        // criterium for MM atom found
        if (shiftMMatom[atom] != -1)
	    {
	        mm_.indexMM_nbl[count] = atom;
	        mm_.shiftMM_nbl[count] = shiftMMatom[atom];
	        count++;
	    }
    }
    // check
    if (count != mm_.nrMMatoms_nbl) {
        printf("ERROR IN MM ATOM SEARCH -- VERLET BASED SCHEME\n");
	    exit(-1);
    }
} // update_QMMMrec_verlet_ns

real QMMM_rec::calculate_QMMM(const t_commrec*      cr,
                              gmx::ForceWithVirial* forceWithVirial,
                              t_nrnb*               nrnb,
                              gmx_wallcycle_t       wcycle)
{
    if (!GMX_QMMM)
    {
        gmx_incons("Compiled without QMMM");
    }

    real QMener = 0.0;
    // A selection for the QM package depending on which is requested
    // (Gaussian, GAMESS-UK, MOPAC or ORCA) needs to be implemented here.
    // Now it works through defines.
    //   ... Not so nice yet

    QMMM_QMrec* qm_ = &(qm[0]);
    QMMM_MMrec* mm_ = &(mm[0]);

    rvec *forces = nullptr,
         *fshift = nullptr;
 
 // gmx::ArrayRef<gmx::RVec> fMM      = forceWithVirial->force_.data();
                  gmx::RVec *fMM      = forceWithVirial->force_.data();
 // gmx::ArrayRef<gmx::RVec> fshiftMM = forceWithShiftForces->shiftForces();

    if (GMX_QMMM_DFTBPLUS)
    {
        snew(forces, (qm_->nrQMatoms + mm_->nrMMatoms + mm_->nrMMatoms_full));
     // snew(fshift, (qm_->nrQMatoms + mm_->nrMMatoms + mm_->nrMMatoms_full));
    }
    else
    {
        snew(forces, (qm_->nrQMatoms + mm_->nrMMatoms));
     // snew(fshift, (qm_.nrQMatoms + mm_.nrMMatoms));
    }

    QMener = call_QMroutine(cr, this, qm_, mm_, forces, fshift, nrnb, wcycle);

    if (GMX_QMMM_DFTBPLUS)
    {
        for (int i = 0; i < qm_->nrQMatoms; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                fMM[qm_->indexQM[i]][j]        -= forces[i][j];
             // fshiftMM[qm_->shiftQM[i]][j]   += fshift[i][j];
            }
         // printf("F[%5d] = %8.2f %8.2f %8.2f\n", qm_->indexQM[i], forces[i][0], forces[i][1], forces[i][2]);
        }
        for (int i = 0; i < mm_->nrMMatoms; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                fMM[mm_->indexMM[i]][j]        -= forces[qm_->nrQMatoms+i][j];
             // fshiftMM[mm_->shiftMM[i]][j]   += fshift[qm_->nrQMatoms+i][j];
            }
         // if (i<30) if (norm(forces[qm_->nrQMatoms+i]) > 10.)
         //   printf("F_MM[%5d] = %8.2f %8.2f %8.2f\n", mm_->indexMM[i],
         //     forces[qm_->nrQMatoms+i][0], forces[qm_->nrQMatoms+i][1], forces[qm_->nrQMatoms+i][2]);
        }
        for (int i = 0; i < mm_->nrMMatoms_full; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                fMM[mm_->indexMM_full[i]][j]        -= forces[qm_->nrQMatoms+mm_->nrMMatoms+i][j];
             // fshiftMM[mm_->shiftMM_full[i]][j]   += fshift[qm_->nrQMatoms+mm_->nrMMatoms+i][j];
            }
         // if (i<100) if (norm(forces[qm_->nrQMatoms+mm_->nrMMatoms+i]) > 10.)
         //   printf("F_MM_F[%5d] = %8.2f %8.2f %8.2f\n", mm_->indexMM_full[i],
         //   forces[qm_->nrQMatoms+mm_->nrMMatoms+i][0], forces[qm_->nrQMatoms+mm_->nrMMatoms+i][1], forces[qm_->nrQMatoms+mm_->nrMMatoms+i][2]);
        }
    }
    else
    {
        for (int i = 0; i < qm_->nrQMatoms; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                fMM[qm_->indexQM[i]][j]          -= forces[i][j];
             // fshiftMM[qm_->shiftQM[i]][j]     += fshift[i][j];
            }
        }
        for (int i = 0; i < mm_->nrMMatoms; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                fMM[mm_->indexMM[i]][j]      -= forces[qm_->nrQMatoms+i][j];
             // fshiftMM[mm_->shiftMM[i]][j] += fshift[qm_->nrQMatoms+i][j];
            }
        }
    }

    sfree(forces);
 // sfree(fshift);

    return QMener;
} // calculate_QMMM

#pragma GCC diagnostic pop
