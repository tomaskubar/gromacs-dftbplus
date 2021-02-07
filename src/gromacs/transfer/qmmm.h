/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2017,2018 by the GROMACS development team.
 * Copyright (c) 2019,2020, by the GROMACS development team, led by
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
#ifndef GMX_TRANSFER_QMMM_H
#define GMX_TRANSFER_QMMM_H

#include "config.h"

#include <vector>

#include "gromacs/math/paddedvector.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/qmmm.h"
#include "gromacs/mdlib/qm_dftbplus.h"
#include "gromacs/mdlib/tgroup.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/timing/wallcycle.h"

#include "gromacs/transfer/qm_dftbplus.h"

struct t_nrnb;
struct nonbonded_verlet_t;

struct gmx_pme_t;
struct t_pbc;
enum class PbcType : int;

struct DftbPlusTransfer;
//struct Context;
struct charge_transfer_t;

class QMMM_PME_transfer {
private:
public:
  PaddedVector<gmx::RVec> x;
  PaddedVector<gmx::RVec> f;
  std::vector<real>       q;
  matrix                  vir;
  real*                   pot; // electrostatic potential from PME ("external shift" in DFTB),
  gmx_bool                surf_corr_pme; // whether the surface correction shall be considered or not (if not = tin-foil boundary cond.)
  real                    epsilon_r; // if yes, this is the dielectric constant to be considered

  rvec*                   qmgrad; // gradients at QM atoms
  rvec*                   mmgrad; // gradients at external charges
  rvec*                   partmmgrad; // temp. array for components of gradients at external charges
  rvec                    com; // center of mass
  real*                   mass; // masses of the atoms
  real                    inv_tot_mass; // 1 / sum(mass)
  real*                   ze; // magnitudes of external charges, NNDIM
  int                     n; // number of QM atoms
  int                     ne; // number of external charges
  int                     cutoff_qmmm;     // whether a switched cut-off QM/MM calculation shall be done instead of PME
  // output
  int                     output_qm_freq;  // how often (if ever) the QM coordinates shall be written
  int                     output_mm_freq;  // how often (if ever) the MM coordinates shall be written
  int                     output_nbl_freq; // how often (if ever) the short-range MM coordinates shall be written
  // for PME
  double                  rcoulomb_pme; // cut-off
  double                  rlist_pme;    // neighborlist cut-off (for PME, equal to rcoulomb_pme; for switched cut-off, larger)
  int                     nstlist_pme;  // frequency of neighborsearching
  int                     lastlist_pme; // last step when neighborsearching was done
} ;

struct gmx_localtop_t;
struct gmx_mtop_t;
struct t_commrec;
struct t_forcerec;
struct t_inputrec;
struct t_mdatoms;

namespace gmx
{
//class ForceWithShiftForces;
class ForceWithVirial;
}

class QMMM_rec_transfer;
class QMMM_QMrec_transfer;
class QMMM_MMrec_transfer;

class QMMM_QMrec_transfer {
private:
    int     nrQMatoms;      // total nr of QM atoms
    int*    indexQM;        // atom i = atom indexQM[i] in mdrun
    int*    atomicnumberQM; // atomic numbers of QM atoms
    real*   QMcharges;      // atomic charges of QM atoms(ONIOM)
    int*    shiftQM;
    int     QMcharge;       // charge of the QM system
 // int     nelectrons;     // total number of elecs in QM region

    matrix  box;
    int     qmmm_variant;
    real    rcoulomb;
    real    ewaldcoeff_q;
    real    epsilon_r;
    double* pot_qmmm;       // electric potential induced by the MM atoms
    double* pot_qmqm;       // el. pot. induced by the periodic images of the QM atoms

    void init_QMrec(const std::vector<int>& qmAtoms,
                 // const int               QMcharge_in,
                    const std::vector<int>& atomicnumberQM_in,
                    const t_inputrec*       ir);

    friend class QMMM_rec_transfer;
    friend void init_dftbplus_transfer(QMMM_QMrec_transfer* qm,
                                       QMMM_rec_transfer*   qr,
                                       const t_inputrec*    ir,
                                       const t_commrec*     cr);
    friend real call_dftbplus_transfer(QMMM_rec_transfer* qr,
                                       const t_commrec*   cr,
                                       rvec              f[],
                                       t_nrnb*           nrnb,
                                       gmx_wallcycle_t   wcycle);

public:
    DftbPlus        *dpcalc;        // DFTB+ calculator
    Context         *dftbContext;   // some data for DFTB+, referenced to by DFTB through *dpcalc

    rvec*            xQM;            // shifted to center of box

    // output
    int              nrQMatoms_get()const;
    int              qmmm_variant_get()const;
    double           xQM_get(int atom, int coordinate)const;
    real             QMcharges_get(int atom)const;
    double           pot_qmmm_get(int atom)const;
    double           pot_qmqm_get(int atom)const;
    int              atomicnumberQM_get(int atom)const;
    int              QMcharge_get()const;
 // int              nelectrons_get()const;
    // input
    void             QMcharges_set(int atom, real value);
    void             pot_qmmm_set(int atom, double value);
    void             pot_qmqm_set(int atom, double value);
} ;

class QMMM_MMrec_transfer {
private:
public:
    real           scalefactor;

    // There are 3 kinds of MM atom lists.
    // (1) the short-range list that is updated in every step of MD,
    //    and is used in the QM calculation:
    int                     nrMMatoms; // nr of MM atoms
    PaddedVector<gmx::RVec> xMM;       // coordinates shifted to the center of the box
    std::vector<int>        indexMM;   // atom i = atom indexMM[i] in mdrun
    std::vector<real>       MMcharges; // magnitude of MM point charges
    std::vector<int>        shiftMM;

    // (2) the short-range list that is produced
    //    by the (group or Verlet) neighborsearching procedure,
    //    and is updated in every neighborsearching step.
    //    The processing of this list in every step of MD
    //    yields the list under (1).
    //    This list itself is not used in QM calculation directly:
    int              nrMMatoms_nbl;
    std::vector<int> indexMM_nbl;
    std::vector<int> shiftMM_nbl;

    // (3) the list of *all* of the non-QM atoms,
    //    which is static throughout the simulation and never needs to be updated:
    int                     nrMMatoms_full;
    PaddedVector<gmx::RVec> xMM_full;
    std::vector<int>        indexMM_full;
    std::vector<real>       MMcharges_full;
    std::vector<int>        shiftMM_full;

    void init_MMrec(real       scalefactor_in,
                    int        nrMMatoms_full_in,
                    int        natoms,
                    int        nrQMatoms,
                    const int* indexQM,
                    int*       found_mm_atoms);
} ;

class QMMM_rec_transfer {
private:
public:

    std::unique_ptr<QMMM_QMrec_transfer>  qm;     // atoms and run params for each QM group
    std::unique_ptr<QMMM_MMrec_transfer>  mm;     // one MM subsystem
    std::vector<QMMM_PME_transfer>        pme;    // [0] == pme_full, [1] == pme_qmonly
    PbcType                               pbcType;
    struct gmx_pme_t* const*              pmedata;

    QMMM_rec_transfer(const t_commrec*  cr,
                      const t_inputrec* ir,
                      const t_forcerec* fr,
                      const int         nAtoms,
                      std::vector<int>& qmAtoms,
                      std::vector<int>& atomicnumberQM,
                      const int         qmmmVariant_in,
                      const bool        dipCorrection);
    // From topology->atoms.atomname and topology->atoms.atomtype
    //   the atom names and types are read;
    // From inputrec->QMcharge resp. inputrec->QMmult the nelecs are determined
    //   and md->cQMMM gives numbers of the MM and QM atoms

    // update_QMMMrec() fills the MM stuff in QMMMrec.
    // The MM atoms are taken from the neighbourlists of the QM atoms.
    // This routine should be called in every NS step,
    // to process the newly created neighborlists!

    void update_QMMMrec_verlet_ns(const t_commrec*          cr,
                                  const nonbonded_verlet_t* nbv,
                                  const rvec                x[],
                                  const t_mdatoms*          md,
                                  const matrix              box);
    
    void update_QMMMrec_dftb(const t_commrec*  cr,
                             const rvec*       shift_vec,
                             const rvec        x[],
                             const t_mdatoms*  md,
                             const matrix      box);
    
    void update_QMMM_coord(const t_commrec*  cr,
                           const rvec*       shift_vec,
                           const rvec        x[],
                           const t_mdatoms*  md,
                           const matrix      box);
    
    // New routines for QM/MM interactions
    
    void calculate_SR_QM_MM(int   variant,
                            real* pot);
    
    void calculate_LR_QM_MM(const t_commrec*  cr,
                            t_nrnb*           nrnb,
                            gmx_wallcycle_t   wcycle,
                            struct gmx_pme_t* pmedata,
                            real*             pot);
    
    void calculate_complete_QM_QM(const t_commrec*  cr,
                                  t_nrnb*           nrnb,
                                  gmx_wallcycle_t   wcycle,
                                  struct gmx_pme_t* pmedata,
                                  real*             pot);
    
    void gradient_QM_MM(const t_commrec*  cr,
                        t_nrnb*           nrnb,
                        gmx_wallcycle_t   wcycle,
                        struct gmx_pme_t* pmedata,
                        int               variant,
                        rvec*             partgrad,
                        rvec*             MMgrad,
                        rvec*             MMgrad_full);

    real calculate_QMMM(const t_commrec*           cr,
                        gmx::ForceWithVirial*      forceWithVirial,
                              t_nrnb*              nrnb,
                              gmx_wallcycle_t      wcycle);
} ;

/* Auxiliary functions */

void atomtype_to_atomicnumber(int atomtype, int& atomicnumber, double& atomicmass);

/*
void put_cluster_in_MMlist_verlet(int                             ck, // cluster number
                                  int                             na_ck, // # of atoms in cluster
                                  int                             nrQMatoms,
                                  const int*                      indexQM,
                                  const gmx::ArrayRef<const int>& atomIndices,
			                      int*                            shiftMMatom,
                                  // ^ also has a role of "bool* isMMatom"
				                  t_pbc*                          pbc,
				                  const rvec*                     x);
*/

/* DFTB+ function/s */

#if (GMX_MPI)
void init_qmmmrec_transfer(std::unique_ptr<QMMM_rec_transfer>* dftbplus_phase1,
                           std::unique_ptr<QMMM_rec_transfer>& dftbplus_phase2,
                           charge_transfer_t*  ct,
                           const t_commrec*    cr,
                           const real          rcoulomb,
                           const real          ewald_rtol,
                           const int           ct_mpi_rank);
#else
void init_qmmmrec_transfer(std::unique_ptr<QMMM_rec_transfer>* dftbplus_phase1,
                           std::unique_ptr<QMMM_rec_transfer>& dftbplus_phase2,
                           charge_transfer_t*  ct,
                           const t_commrec*    cr,
                           const real          rcoulomb,
                           const real          ewald_rtol);
#endif


/* Functions defined in qm_dftbplus.cpp */

void
init_dftbplus_transfer(QMMM_QMrec_transfer* qm,
                       QMMM_rec_transfer*   qr,
                       const t_inputrec*    ir,
                       const t_commrec*     cr);

real
call_dftbplus_transfer(QMMM_rec_transfer*   qr,
                       const t_commrec*     cr,
                       rvec                 f[],
                       t_nrnb*              nrnb,
                       gmx_wallcycle_t      wcycle);

void initialize_context_transfer(Context*           cont,
                                 int                nrQMatoms,
                                 int                qmmm_variant,
                                 QMMM_rec_transfer* qr_in,
                                 const t_inputrec*  ir_in,
                                 const t_commrec*   cr_in);

/* End of qm_dftbplus.cpp */

#endif
