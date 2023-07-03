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
#ifndef GMX_MDLIB_QMMM_H
#define GMX_MDLIB_QMMM_H

#include "config.h"

#include <vector>

#include "gromacs/math/paddedvector.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/tgroup.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/timing/wallcycle.h"
//#include "gromacs/mdlib/qm_nn.h"
extern "C" {
    #include "tensorflow/c/c_api.h"
}

//#include "gromacs/mdlib/qm_dftbplus.h"
//#include "gromacs/mdlib/qm_gamess.h"
//#include "gromacs/mdlib/qm_gaussian.h"
//#include "gromacs/mdlib/qm_mopac.h"
//#include "gromacs/mdlib/qm_orca.h"

#define GMX_QMMM (GMX_QMMM_MOPAC || GMX_QMMM_GAMESS || GMX_QMMM_GAUSSIAN || GMX_QMMM_ORCA || GMX_QMMM_DFTBPLUS || GMX_QMMM_NN)

struct t_nrnb;
struct nonbonded_verlet_t;

struct gmx_pme_t;
enum class PbcType : int;

struct DftbPlus;
struct Context;

typedef struct {
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
}TFModel;

// THIS STRUCTURE IS TENTATIVE,
// JUST FOR THE BEGINNING.
// MANY THINGS ARE STORED TWICE!
// TO BE CLEANED UP!
// Also, it would be cool to rename the structure.
class QMMM_PME {
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

class QMMM_QMrec;
class QMMM_MMrec;

class QMMM_QMgaussian {
private:
public:
    int      nQMcpus;  // no. of CPUs used for the QM calc.
    int      QMmem;    // memory for the gaussian calc.
    int      accuracy; // convergence criterium (E(-x))
    gmx_bool cpmcscf;  // using cpmcscf(l1003)
    char*    gauss_dir;
    char*    gauss_exe;
    char*    devel_dir;

/*! \brief
 * Initialize gaussian datastructures.
 */
void init_gaussian();

/*! \brief
 * Call gaussian to do qm calculation.
 *
 * \param[in] fr Global forcerec.
 * \param[in] qm QM part of forcerec.
 * \param[in] mm mm part of forcerec.
 * \param[in] f  force vector.
 * \param[in] fshift shift of force vector.
 */
real call_gaussian(const QMMM_QMrec& qm,
                   const QMMM_MMrec& mm,
                   rvec              f[],
                   rvec              fshift[]) const;
} ;

class QMMM_QMrec {
private:
    int     nrQMatoms;      // total nr of QM atoms
    rvec*   xQM;            // shifted to center of box
    int*    indexQM;        // atom i = atom indexQM[i] in mdrun
    int*    atomicnumberQM; // atomic numbers of QM atoms
    real*   QMcharges;      // atomic charges of QM atoms(ONIOM)
    int*    shiftQM;
    int     QMcharge;       // charge of the QM system
    int     multiplicity;   // multipicity (no of unpaired eln)
    int     QMmethod;       // see enums.h for all methods
    int     QMbasis;        // see enums.h for all bases
    int     nelectrons;     // total number of elecs in QM region
    int     CASelectrons;   // electron in active space
    int     CASorbitals;    // orbitals in active space

    matrix  box;
    int     qmmm_variant;
    real    rcoulomb;
    real    ewaldcoeff_q;
    real    epsilon_r;
    double* pot_qmmm;       // electric potential induced by the MM atoms
    double* pot_qmqm;       // el. pot. induced by the periodic images of the QM atoms

    void init_QMrec(int               grpnr,
                    int               nr,
                    const int*        atomarray,
                    const gmx_mtop_t* mtop,
                    const t_inputrec* ir);

    friend class QMMM_rec;
    friend void init_dftbplus(QMMM_QMrec*       qm,
                              const t_forcerec* fr,
                              const t_inputrec* ir,
                              const t_commrec*  cr,
                              gmx_wallcycle_t   wcycle);
    friend real call_dftbplus(const t_forcerec* fr,
                              const t_commrec*  cr,
                              QMMM_QMrec*       qm,
                              QMMM_MMrec*       mm,
                              rvec              f[],
                              rvec              fshift[],
                              t_nrnb*           nrnb,
                              gmx_wallcycle_t   wcycle);
    friend void init_nn(QMMM_QMrec*       qm);
    friend real call_nn(const t_forcerec* fr,
                              const t_commrec*  cr,
                              QMMM_QMrec*       qm,
                              QMMM_MMrec*       mm,
                              rvec              f[],
                              rvec              fshift[],
                              t_nrnb*           nrnb,
                              gmx_wallcycle_t   wcycle);

public:
    DftbPlus        *dpcalc;        // DFTB+ calculator
    Context         *dftbContext;   // some data for DFTB+, referenced to by DFTB through *dpcalc
    TFModel         *model;

    QMMM_QMgaussian  gaussian;

    char            *orca_basename; // basename for I/O with orca
    char            *orca_dir;      // directory for ORCA

    // output
    int              nrQMatoms_get()const;
    int              qmmm_variant_get()const;
    double           xQM_get(int atom, int coordinate)const;
    real             QMcharges_get(int atom)const;
    double           pot_qmmm_get(int atom)const;
    double           pot_qmqm_get(int atom)const;
    int              atomicnumberQM_get(int atom)const;
    int              QMcharge_get()const;
    int              multiplicity_get()const;
    int              QMmethod_get()const;
    int              QMbasis_get()const;
    int              nelectrons_get()const;
    int              CASelectrons_get()const;
    int              CASorbitals_get()const;
    // input
    void             QMcharges_set(int atom, real value);
    void             pot_qmmm_set(int atom, double value);
    void             pot_qmqm_set(int atom, double value);
} ;

class QMMM_MMrec {
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

class QMMM_rec {
private:
public:

 // int                      nrQMlayers; // number of QM groups/layers (total layers +1 (MM))
    std::vector<QMMM_QMrec>  qm;         // atoms and run params for each QM group
    std::vector<QMMM_MMrec>  mm;         // there can only be one MM subsystem !
    std::vector<QMMM_PME>    pme;        // [0] == pme_full, [1] == pme_qmonly
    PbcType                  pbcType;
    struct gmx_pme_t* const* pmedata;

    QMMM_rec(const t_commrec*                 cr,
             const gmx_mtop_t*                mtop,
             const t_inputrec*                ir,
             const t_forcerec*                fr);
          // const gmx_wallcycle_t gmx_unused wcycle);
    // From topology->atoms.atomname and topology->atoms.atomtype
    //   the atom names and types are read;
    // From inputrec->QMcharge resp. inputrec->QMmult the nelecs and multiplicity are determined
    //   and md->cQMMM gives numbers of the MM and QM atoms

    /*
    void update_QMMMrec(const t_commrec*  cr,
                        const t_forcerec* fr,
                        const rvec*       x,
                        const t_mdatoms*  md,
                        const matrix      box);
    */

    // update_QMMMrec() fills the MM stuff in QMMMrec.
    // The MM atoms are taken from the neighbourlists of the QM atoms.
    // In a QMMM run, this routine should be called at every step,
    //   since it updates the MM elements of the t_QMMMrec struct.
    // CORRECTION: only call it in every NS step,
    // to process the newly created neighborlists!

    void update_QMMMrec_verlet_ns(const t_commrec*    cr,
                                  nonbonded_verlet_t* nbv,
                                  const rvec          x[],
                                  const t_mdatoms*    md,
                                  const matrix        box);
    
    void update_QMMMrec_dftb(const t_commrec*  cr,
                             rvec*             shift_vec,
                             const rvec        x[],
                             const t_mdatoms*  md,
                             const matrix      box);
    
    void update_QMMM_coord(const t_commrec*  cr,
                           rvec*             shift_vec,
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

    // QMMM computes the QM forces.
    // This routine makes either function calls to gmx QM routines
    //   (derived from MOPAC7 (semi-emp.) and MPQC (ab initio))
    // or generates input files for an external QM package
    //   (listed in QMMMrec.QMpackage).
    // The binary of the QM package is called by system().
} ;

// for qmmm_variant
enum {eqmmmVACUO,
      eqmmmPME,
      eqmmmSWITCH,
      eqmmmRFIELD,
      eqmmmSHIFT,
      eqmmmNN,
      eqmmmNR};

/*! \brief
 * Return vector of atom indices for atoms in the QMMM region.
 *
 * \param[in] mtop Topology to use for populating array.
 * \param[in] ir   Inputrec used in simulation.
 * \returns Vector of atoms.
 */
std::vector<int> qmmmAtomIndices(const t_inputrec& ir,
                                 const gmx_mtop_t& mtop);

/*! \brief
 * Remove charges from QMMM atoms.
 *
 * \param[in] mtop Topology used for removing atoms.
 * \param[in] qmmmAtoms ArrayRef to vector conatining qmmm atom indices.
 */
void removeQmmmAtomCharges(gmx_mtop_t*              mtop,
                           gmx::ArrayRef<const int> qmmmAtoms);

#endif
