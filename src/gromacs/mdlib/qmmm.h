/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2017,2018,2019, by the GROMACS development team, led by
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

#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/tgroup.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/timing/wallcycle.h"

//#include "gromacs/mdlib/qm_dftbplus.h"

struct t_nrnb;

struct gmx_pme_t;

struct DftbPlus;

/* THIS STRUCTURE IS TENTATIVE,
 * JUST FOR THE BEGINNING.
 * MANY THINGS ARE STORED TWICE!
 * TO BE CLEANED UP!
 * Also, it would be cool to rename the structure.
 */
typedef struct { // TODO CHECK WHERE ALL OF THIS IS ALLOCATED!
  struct gmx_pme_t **pmedata_qmmm;
  struct gmx_pme_t **pmedata_qmqm;
  t_nrnb *nrnb;
  rvec *x;
  real *q;
  rvec *f;
  matrix vir;
  real *pot, *pot_sr; // electrostatic potential from PME ("external shift" in DFTB),
                      // total and short-range component
  gmx_bool surf_corr_pme;   /* whether the surface correction shall be considered or not (if not = tin-foil boundary cond. */
  real epsilon_r; /* if yes, this is the dielectric constant to be considered */

  rvec *qmgrad, /* gradients at QM atoms */
       *mmgrad, /* gradients at external charges */
       *partmmgrad, /* temp. array for components of gradients at external charges */
       com; /* center of mass */
  real *mass, /* masses of the atoms */
        inv_tot_mass, /* 1 / sum(mass) */
       *ze; /* magnitudes of external charges, NNDIM */
  int n, /* number of QM atoms */
      ne; /* number of external charges */
  int cutoff_qmmm;     /* whether a switched cut-off QM/MM calculation shall be done instead of PME */
  /* output */
  int output_qm_freq;  /* how often (if ever) the QM coordinates shall be written */
  int output_mm_freq;  /* how often (if ever) the MM coordinates shall be written */
  int output_nbl_freq; /* how often (if ever) the short-range MM coordinates shall be written */
  /* for PME */
  double rcoulomb_pme; // cut-off
  double rlist_pme; // neighborlist cut-off (for PME, equal to rcoulomb_pme; for switched cut-off, larger)
  int nstlist_pme, lastlist_pme; // frequency of neighborsearching; last step when neighborsearching was done
} t_QMMM_PME;

#define GMX_QMMM (GMX_QMMM_MOPAC || GMX_QMMM_GAMESS || GMX_QMMM_GAUSSIAN || GMX_QMMM_ORCA || GMX_QMMM_DFTBPLUS)

struct gmx_localtop_t;
struct gmx_mtop_t;
struct t_commrec;
struct t_forcerec;
struct t_inputrec;
struct t_mdatoms;
struct t_QMMMrec;

typedef struct {
    int                nrQMatoms;      /* total nr of QM atoms              */
    rvec              *xQM;            /* shifted to center of box          */
    int               *indexQM;        /* atom i = atom indexQM[i] in mdrun */
    int               *atomicnumberQM; /* atomic numbers of QM atoms        */
    real              *QMcharges;      /* atomic charges of QM atoms(ONIOM) */
    int               *shiftQM;
    int                QMcharge;       /* charge of the QM system           */
    int                multiplicity;   /* multipicity (no of unpaired eln)  */
    int                QMmethod;       /* see enums.h for all methods       */
    int                QMbasis;        /* see enums.h for all bases         */
    int                nelectrons;     /* total number of elecs in QM region*/
    /* Gaussian specific stuff */
    int                nQMcpus;        /* no. of CPUs used for the QM calc. */
    int                QMmem;          /* memory for the gaussian calc.     */
    int                accuracy;       /* convergence criterium (E(-x))     */
    gmx_bool           cpmcscf;        /* using cpmcscf(l1003)*/
    char              *gauss_dir;
    char              *gauss_exe;
    char              *devel_dir;
    char              *orca_basename; /* basename for I/O with orca        */
    char              *orca_dir;      /* directory for ORCA                */
    /* Surface hopping stuff */
    gmx_bool           bSH;           /* surface hopping (diabatic only)   */
    real               SAon;          /* at which energy gap the SA starts */
    real               SAoff;         /* at which energy gap the SA stops  */
    int                SAsteps;       /* stepwise switchinng on the SA     */
    int                SAstep;        /* current state of SA               */
    int                CIdim;
    real              *CIvec1;
    real              *CIvec2;
    real              *CIvec1old;
    real              *CIvec2old;
    ivec               SHbasis;
    int                CASelectrons;
    int                CASorbitals;

    matrix             box;
    int                qmmm_variant;
    real               rcoulomb;
    real               ewaldcoeff_q;
    real               epsilon_r;
    DftbPlus          *dpcalc;        /* DFTB+ calculator */
} t_QMrec;

typedef struct {
    real           scalefactor;
 // int            nrMMatoms;   /* nr of MM atoms, updated every step*/
 // rvec          *xMM;         /* shifted to center of box          */
 // int           *indexMM;     /* atom i = atom indexMM[I] in mdrun */
 // real          *MMcharges;   /* MM point charges in std QMMM calc.*/
 // int           *shiftMM;
    /* There are 3 kinds of MM atom lists. */
    /* (1) the short-range list that is updated in every step of MD,
          and is used in the QM calculation: */
    int            nrMMatoms;   /* nr of MM atoms */
    rvec          *xMM;         /* coordinates shifted to the center of the box */
    int           *indexMM;     /* atom i = atom indexMM[I] in mdrun */
    real          *MMcharges;   /* magnitude of MM point charges */
    int           *shiftMM;
    /* (2) the short-range list that is produced
          by the (group or Verlet) neighborsearching procedure,
          and is updated in every neighborsearching step.
          The processing of this list in every step of MD
          yields the list under (1).
          This list itself is not used in QM calculation directly: */
    int            nrMMatoms_nbl;
    int           *indexMM_nbl;
    int           *shiftMM_nbl;
    /* (3) the list of *all* of the non-QM atoms,
          which is static throughout the simulation and never needs to be updated: */
    int            nrMMatoms_full;
    rvec          *xMM_full;
    int           *indexMM_full;
    real          *MMcharges_full;
    int           *shiftMM_full;

   /* the following seem to be unused in the new version of the QM/MM interface */
 // int           *MMatomtype;  /* only important for semi-emp.      */
} t_MMrec;


typedef struct t_QMMMrec {
 // int             QMMMscheme; /* ONIOM (multi-layer) or normal          */
 // int             nrQMlayers; /* number of QM layers (total layers +1 (MM)) */
    t_QMrec       **qm;         /* atoms and run params for each QM group */
    t_MMrec        *mm;         /* there can only be one MM subsystem !   */
    t_QMMM_PME     *pme;
} t_QMMMrec;

// for qmmm_variant
enum {eqmmmVACUO,
      eqmmmPME,
      eqmmmSWITCH,
      eqmmmRFIELD,
      eqmmmSHIFT,
      eqmmmNR};

void atomic_number(int nr, char ***atomtype, int *nucnum);

t_QMMMrec *mk_QMMMrec();
/* allocates memory for QMMMrec */

void init_QMMMrec(const t_commrec  *cr,
                  const gmx_mtop_t *mtop,
                  const t_inputrec *ir,
                  const t_forcerec *fr,
                  const gmx_wallcycle_t wcycle);                

/* init_QMMMrec initializes the QMMM record. From
 * topology->atoms.atomname and topology->atoms.atomtype the atom
 * names and types are read; from inputrec->QMcharge
 * resp. inputrec->QMmult the nelecs and multiplicity are determined
 * and md->cQMMM gives numbers of the MM and QM atoms
 */

/*
void update_QMMMrec(const t_commrec  *cr,
                    const t_forcerec *fr,
                    const rvec       *x,
                    const t_mdatoms  *md,
                    const matrix      box);
*/

/* update_QMMMrec fills the MM stuff in QMMMrec. The MM atoms are
 * taken froom the neighbourlists of the QM atoms. In a QMMM run this
 * routine should be called at every step, since it updates the MM
 * elements of the t_QMMMrec struct.
 * CORRECTION: only call it in every NS step, to process the newly
 * created neighborlists!
 * There are two separate routines, for group and Verlet nblists.
 */

void update_QMMMrec_verlet_ns(const t_commrec      *cr,
                                    t_forcerec     *fr,
                              const rvec            x[],
                              const t_mdatoms      *md,
                              const matrix          box);

void update_QMMMrec_dftb(const t_commrec      *cr,
                               t_forcerec     *fr,
                         const rvec            x[],
                         const t_mdatoms      *md,
                         const matrix          box);

void update_QMMM_coord(const t_commrec  *cr,
                             t_forcerec *fr,
                       const rvec        x[],
                       const t_mdatoms  *md,
                       const matrix      box);

real calculate_QMMM(const t_commrec      *cr,
                    const t_forcerec     *fr,
                          rvec            f[],
                    const gmx_wallcycle_t wcycle);

/* QMMM computes the QM forces. This routine makes either function
 * calls to gmx QM routines (derived from MOPAC7 (semi-emp.) and MPQC
 * (ab initio)) or generates input files for an external QM package
 * (listed in QMMMrec.QMpackage). The binary of the QM package is
 * called by system().
 */

/*! \brief
 * Return vector of atom indices for atoms in the QMMM region.
 *
 * \param[in] mtop Topology to use for populating array.
 * \param[in] ir   Inputrec used in simulation.
 * \returns Vector of atoms.
 */
std::vector<int> qmmmAtomIndices(const t_inputrec &ir, const gmx_mtop_t &mtop);

/*! \brief
 * Remove charges from QMMM atoms.
 *
 * \param[in] mtop Topology used for removing atoms.
 * \param[in] qmmmAtoms ArrayRef to vector conatining qmmm atom indices.
 */
void removeQmmmAtomCharges(gmx_mtop_t *mtop, gmx::ArrayRef<const int> qmmmAtoms);

/* New routines for QM/MM interactions */

void calculate_SR_QM_MM(t_QMMMrec *qr,
                        int variant,
                        real *pot);

void calculate_LR_QM_MM(t_QMMMrec *qr,
                        const t_commrec *cr,
                        gmx_wallcycle_t wcycle,
                        struct gmx_pme_t *pmedata,
                        real *pot);

void calculate_complete_QM_QM(t_QMMMrec *qr,
                              const t_commrec *cr,
                              gmx_wallcycle_t wcycle,
                              struct gmx_pme_t *pmedata,
                              real *pot);

void gradient_QM_MM(t_QMMMrec *qr,
                    const t_commrec *cr,
                    gmx_wallcycle_t wcycle,
                    struct gmx_pme_t *pmedata,
                    int variant,
                    rvec *partgrad,
                    rvec *MMgrad,
                    rvec *MMgrad_full);

#endif
