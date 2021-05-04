/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2017, by the GROMACS development team, led by
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

#include "config.h"

// When not built in a configuration with QMMM support, much of this
// code is unreachable by design. Tell clang not to warn about it.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-noreturn"

#if GMX_QMMM_DFTBPLUS

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gromacs/fileio/confio.h"
#include "gromacs/ewald/pme.h"
//#include "gromacs/ewald/pme-internal.h"
#include "gromacs/ewald/ewald_utils.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/qmmm.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/timing/cyclecounter.h"
#include "gromacs/timing/walltime_accounting.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/pbcutil/pbc.h"

//#include "dftbplus_gromacs.h"
#include "gromacs/mdlib/qm_dftbplus.h"

typedef struct Context {
  bool              pme;
  int               n;
  const t_commrec*  cr;
  QMMM_rec*         qr;
//const t_forcerec* fr;
  t_nrnb*           nrnb;
  gmx_wallcycle_t   wcycle;
  real              rcoul;
  real              ewaldcoeff_q;
} Context;

// extern "C" {void readdftbplusinput();}
// extern "C" {void performdftbpluscalculation(real *e, int *n, real x[], real q[], real extshift[], real f[]);}

/* Fill up the context (status structure) with all the relevant data
 */
void initialize_context(Context*          cont,
                        int               nrQMatoms,
                        int               qmmm_variant,
                        QMMM_rec*         qr_in,
                     // const t_forcerec* fr_in,
                        const t_inputrec* ir_in,
                        const t_commrec*  cr_in)
                     // gmx_wallcycle_t   wcycle_in)
//                      const real        rcoul_in,
//                      const real        ewaldcoeff_q_in)
{
  /* The "cr" and "qr" structures will be initialized
   *   at the start of every MD step,
   *   and will be remembered for all of the SCC iterations.
   * That way, these data need not pass through the DFTB+ program.
   *
   * NOTE: In the new "C++" implementation,
   *       this is called from the constructor QMMM_rec::QMMM_rec(),
   *       so the object qr is not available yet!
   */

  if (qmmm_variant == eqmmmPME) {
      cont->pme = true;
  } else {
      cont->pme = false;
  }
  printf("qmmm_variant = %d\n", qmmm_variant);
  printf("cont->pme = %s\n", cont->pme ? "true" : "false");
//cont->n = fr_in->qr->qm[0].nrQMatoms; // object qr not available yet!
  cont->n = nrQMatoms;
  printf("cont->n = %d\n", cont->n);
  if (cont->pme)
  {
      cont->cr           = cr_in;
      cont->qr           = qr_in;
   // cont->wcycle       = wcycle_in;
      cont->rcoul        = ir_in->rcoulomb;
      cont->ewaldcoeff_q = calc_ewaldcoeff_q(ir_in->rcoulomb, ir_in->ewald_rtol);
      printf("cont->cr = %p\n", cont->cr);
      printf("cont->qr = %p\n", cont->qr);
      printf("cont->qr->pmedata = %p\n", cont->qr->pmedata);
      printf("cont->qr->pmedata* = %p\n", *cont->qr->pmedata);
      printf("cont->rcoul = %f\n", cont->rcoul);
      printf("cont->ewaldcoeff_q = %f\n", cont->ewaldcoeff_q);
  }

  return;
}

/* Calculate the external potential due to periodic images of QM atoms with PME.
 */
void calcQMextPotPME(Context *cont, double *q, double *extpot)
{
//int n = cont->fr->qr->qm[0]->nrQMatoms;
  int n = cont->n;

//printf("calcQMextPotPME: %s\n", setup ? "setup" : "calculation");

  if (cont->pme)
  {
      /* PERFORM THE REAL CALCULATION */
      real *extpot_real;
      snew(extpot_real, n);
      for (int i=0; i<n; i++)
      {
          cont->qr->qm[0].QMcharges_set(i, (real) -q[i]); // check sign TODO
      }
      cont->qr->calculate_complete_QM_QM(cont->cr, cont->nrnb, cont->wcycle, *cont->qr->pmedata, extpot_real);
      for (int i=0; i<n; i++)
      {
          extpot[i] = (double) - extpot_real[i]; // sign OK
      }
      sfree(extpot_real);
  }
  else
  {
      for (int i=0; i<n; i++)
      {
          extpot[i] = 0.;
      }
  }

  return;
}

/* The wrapper for the calculation of external potential
 *   due to the periodic images of QM atoms with PME.
 */
extern "C" void calcqmextpot(void *refptr, double *q, double *extpot)
{
  Context *cont = (Context *) refptr;
  calcQMextPotPME(cont, q, extpot);
  return;
}

/* The wrapper for the calculation of grdient of external potential
 *   due to the periodic images of QM atoms with PME.
 * This is function does not calculate anything,
 *   because it is not needed in this DFTB+Gromacs implementation.
 */
extern "C" void calcqmextpotgrad(void *refptr, gmx_unused double *q, double *extpotgrad)
{
  Context *cont = (Context *) refptr;
  for (int i=0; i<3*cont->n; i++)
      extpotgrad[i] = 0.;
  return;
}

/* DFTBPLUS interface routines */

void init_dftbplus(QMMM_QMrec*       qm,
                   QMMM_rec*         qr,
                // const t_forcerec* fr,
                   const t_inputrec* ir,
                   const t_commrec*  cr)
                // gmx_wallcycle_t   wcycle)
//void init_dftbplus(t_forcerec *fr)
{
    /* perhaps check the geometry first, to see which elements we have? */

    // variables
    int nAtom;
    int nSpecies;
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

    static DftbPlus  calculator;
    DftbPlusInput    input;
    DftbPlusAtomList atomList;

    /* This structure will be passed through DFTB+
     *   into the Gromacs calculator calcQMextPot
     */
  //static Context *cont;
  //snew(cont, 1);

  //initialize_context(cont, qm->nrQMatoms, qm->qmmm_variant, fr, ir, cr, nrnb, wcycle);
    snew(qm->dftbContext, 1);
    initialize_context(qm->dftbContext, qm->nrQMatoms_get(), qm->qmmm_variant_get(), qr, ir, cr); //, wcycle);
  //qm->dftbContext = cont;

    snew(qm->dpcalc, 1);

    /* Fill up the context with all the relevant data! */

    /* Initialize the DFTB+ calculator */
    dftbp_init(&calculator, "dftb_in.out");
    printf("DFTB+ calculator has been created!\n");

    /* Parse the input file and store the input-tree */
    dftbp_get_input_from_file(&calculator, "dftb_in.hsd", &input);
    printf("DFTB+ input has been read!\n");

    /* Pass the list of QM atoms to DFTB+ */
    nAtom = qm->nrQMatoms_get();
    snew(ptrSpecies, nAtom);
    snew(atomicNumber, nAtom);
    // read the atomic numbers of QM atoms, and set up the lists
    nSpecies = 0;
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
    printf("No. of QM atoms = %d\n", nAtom);
    printf("No. of chem. species = %d:", nSpecies);
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
    dftbp_process_input(&calculator, &input, &atomList);
    printf("DFTB+ input has been processed!\n");

    qm->dpcalc = &calculator;

    /* Register the callback functions which calculate
     * the external potential and its gradient
     */
    dftbp_register_ext_pot_generator(qm->dpcalc,
                                   //cont,
                                     qm->dftbContext,
                                     calcqmextpot,
                                     calcqmextpotgrad);

    /* Alternatively - modify the geometry/elements information here,
     * based on what we have in the Gromacs topology!
     */

    return;
} /* init_dftbplus */

real call_dftbplus(QMMM_rec*         qr,
                   const t_commrec*  cr,
                   QMMM_QMrec*       qm,
                   const QMMM_MMrec& mm,
                   rvec              f[],
                   rvec              fshift[],
		           t_nrnb*           nrnb,
                   gmx_wallcycle_t   wcycle)
{
    static int step = 0;
    static FILE *f_q = nullptr;
    static FILE *f_p = nullptr;
    static FILE *f_x_qm = nullptr;
    static FILE *f_x_mm = nullptr;
    static int output_freq_q;
    static int output_freq_p;
    static int output_freq_x_qm;
    static int output_freq_x_mm;

    double QMener;
 // bool lPme = (qm->qmmm_variant == eqmmmPME);
    (void) fshift;

    qm->dftbContext->nrnb = nrnb;
    int n = qm->nrQMatoms_get();

    double *x, *grad, *pot, *potgrad, *q, *xMM; // real instead of rvec, to help pass data to fortran
    real *pot_sr = nullptr, *pot_lr = nullptr;
    rvec *QMgrad = nullptr, *MMgrad = nullptr, *MMgrad_full = nullptr;

    snew(x, 3*n);
    snew(grad, 3*n);
    for (int i=0; i<3*n; i++)
        grad[i] = 0.;
    snew(pot, n);
    snew(potgrad, 3*n); // dummy parameter; not used at this moment
    for (int i=0; i<3*n; i++)
        potgrad[i] = 0.;
    snew(q, n);
    for (int i=0; i<n; i++)
        q[i] = 0.;
    snew(pot_sr, n);
    snew(pot_lr, n);

    snew(QMgrad, qm->nrQMatoms_get());

    if (step == 0) {
        char *env;

        if ((env = getenv("GMX_DFTB_CHARGES")) != nullptr)
        {
            output_freq_q = atoi(env);
            f_q = fopen("qm_dftb_charges.xvg", "a");
            printf("The QM charges will be saved in file qm_dftb_charges.xvg every %d steps.\n", output_freq_q);
        }

        if (qm->qmmm_variant_get() != eqmmmVACUO && (env = getenv("GMX_DFTB_ESP")) != nullptr)
        {
            output_freq_p = atoi(env);
            f_p = fopen("qm_dftb_esp.xvg", "a");
            printf("The MM potential induced on QM atoms will be saved in file qm_dftb_esp.xvg every %d steps.\n", output_freq_p);
        }

        if ((env = getenv("GMX_DFTB_QM_COORD")) != nullptr)
        {
            output_freq_x_qm = atoi(env);
            f_x_qm = fopen("qm_dftb_qm.qxyz", "a");
            printf("The QM coordinates (XYZQ) will be saved in file qm_dftb_qm.qxyz every %d steps.\n", output_freq_x_qm);
        }

        if (qm->qmmm_variant_get() != eqmmmVACUO && (env = getenv("GMX_DFTB_MM_COORD")) != nullptr)
        {
            output_freq_x_mm = atoi(env);
            f_x_mm = fopen("qm_dftb_mm.qxyz", "a");
            printf("The MM coordinates (XYZQ) will be saved in file qm_dftb_mm.qxyz every %d steps.\n", output_freq_x_mm);
        }
    }

    for (int i=0; i<n; i++)
    {
        for (int j=0; j<DIM; j++)
        {
            x[3*i+j] = qm->xQM_get(i,j) / BOHR2NM; // to bohr units for DFTB+
        }
    }
    if (qm->coupledPerturb && mm.nrMMatoms > 0)
    {
        snew(xMM, 3 * mm.nrMMatoms);
        for (int i=0; i<mm.nrMMatoms; i++)
        {
            for (int j=0; j<DIM; j++)
            {
                xMM[3*i+j] = mm.xMM[i][j] / BOHR2NM; // to bohr units for DFTB+
            }
        }
    }

    if (!qm->coupledPerturb)
    {
        /* calculate the QM/MM electrostatics beforehand with Gromacs!
         * with cut-off treatment, this will be the entire QM/MM
         * with PME, this will only include MM atoms,
         *    and the contribution from QM atoms will be added in every SCC iteration
         */
        qr->calculate_SR_QM_MM(qm->qmmm_variant_get(), pot_sr);
        if (qm->qmmm_variant_get() == eqmmmPME)
        {
         // gmx_pme_init_qmmm(&(qr->pme->pmedata), true, fr->pmedata);
            qr->calculate_LR_QM_MM(cr, nrnb, wcycle, *qr->pmedata, pot_lr);
            for (int i=0; i<n; i++)
            {
                pot[i] = (double) - (pot_sr[i] + pot_lr[i]);
            }
        }
        else
        {
            for (int i=0; i<n; i++)
            {
                pot_lr[i] = 0.;
                pot[i] = (double) - pot_sr[i];
            }
        }
        // save the potential in the QMMM_QMrec structure
        for (int j=0; j<n; j++)
        {
            qm->pot_qmmm_set(j, (double) - pot[j] * HARTREE_TO_EV); // in volt units
        }
     // // DEBUG
     // for (int i=0; i<n; i++)
     //     printf("pot_sr[%d] = %9.5f pot_lr[%d] = %9.5f\n", i+1, pot_sr[i], i+1, pot_lr[i]);

        /* Set up the data structures needed for the PME calculation
         *   of QM--imageQM electrostatics.
         * During the iterative SCC calculation, this routine will be called
         *   directly from DFTB+, and will use those data structures.
         */
     //  calcQMextPotPME(nullptr, nullptr, true, lPme, fr, cr, wcycle, fr->rcoulomb, fr->ewaldcoeff_q);
        /* This was already done in the initialization procedure! */
    }

    for (int i=0; i<n; i++)
    {
        q[i] = 0.;
    }

    /* DFTB+ calculation itself */
    wallcycle_start(wcycle, ewcQM);
    dftbp_set_coords(qm->dpcalc, x); // unit OK
    if (qm->coupledPerturb)
    {
        // with solution of coupled-perturbed equations:
        //   pass MM charges to DFTB+ and let it calculate
        //   everything except forces on MM atoms
        dftbp_set_external_charges(qm->dpcalc, &(mm.nrMMatoms), xMM, mm.MMcharges.data());
    }
    else
    {
        // if coupled-perturbed equations are not solved,
        //   use the routines implemented in Gromacs
        dftbp_set_external_potential(qm->dpcalc, pot, potgrad); // unit and sign OK
    }
    dftbp_get_energy(qm->dpcalc, &QMener); // unit OK
    dftbp_get_gross_charges(qm->dpcalc, q);
 // for (int i=0; i<n; i++)
 //     printf("%d %6.3f\n", i+1, q[i]);
    dftbp_get_gradients(qm->dpcalc, grad);
    if (qm->coupledPerturb)
    {
        dftbp_get_charge_derivatives(qm->dpcalc, qm->QMchargeGradients, qm->QMchargeMMgradients);
        for (int i = 0; i < qm->nrQMatoms_get() * 3 * qm->nrQMatoms_get(); i++)
        {
            qm->QMchargeGradients[i] /= BOHR2NM;
        }
        for (int i = 0; i < qm->nrQMatoms_get() * 3 * qm->nrQMchargeMMatoms; i++)
        {
            qm->QMchargeMMgradients[i] /= BOHR2NM;
        }
    }

    wallcycle_stop(wcycle, ewcQM);

    /* Save the gradient on the QM atoms */
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<3; j++)
        {
            QMgrad[i][j] = (real) grad[3*i+j]; // negative of force -- sign OK
        }
    }

 // /* Print the QM pure gradient */
 // for (int i=0; i<n; i++)
 // {
 //     printf("GRAD QM %d: %8.2f %8.2f %8.2f\n", i+1,
 //         QMgrad[i][XX] * HARTREE_BOHR2MD, QMgrad[i][YY] * HARTREE_BOHR2MD, QMgrad[i][ZZ] * HARTREE_BOHR2MD);
 // }

    /* Save the QM charges */
    for (int i=0; i<n; i++)
    {
        qm->QMcharges_set(i, (real) q[i]); // sign OK
     // printf("CHECK CHARGE QM[%d] = %6.3f\n", i+1, qm->QMcharges[i]);
    }

    /* Calculate the QM/MM forces
     *   (these are not covered by DTFB+,
     *    because it does not know the position and charges of individual atoms,
     *    and instead, it only obtains the external potentials induced at QM atoms.
     */
    snew(MMgrad, mm.nrMMatoms);
    if (qm->qmmm_variant_get() == eqmmmPME)
    {
        snew(MMgrad_full, mm.nrMMatoms_full);
    }

    if (!qm->coupledPerturb)
    {
        rvec *partgrad;
        snew(partgrad, qm->nrQMatoms_get());
        qr->gradient_QM_MM(cr, nrnb, wcycle, (qm->qmmm_variant_get() == eqmmmPME ? *qr->pmedata : nullptr),
                       qm->qmmm_variant_get(), partgrad, MMgrad, MMgrad_full);
        for (int i=0; i<n; i++)
        {
            rvec_inc(QMgrad[i], partgrad[i]); // sign OK
         // printf("GRAD QM FULL %d: %8.2f %8.2f %8.2f\n", i+1,
         //     QMgrad[i][XX] * HARTREE_BOHR2MD, QMgrad[i][YY] * HARTREE_BOHR2MD, QMgrad[i][ZZ] * HARTREE_BOHR2MD);
        }
        sfree(partgrad);
    }
    else
    {
        // QM gradient has been covered by DFTB+, so the only remaining thing
        //   is the QM/MM gradient on MM charges here
        qr->gradient_MM_only(cr, nrnb, wcycle, MMgrad);
    }

    /* Put the QMMM forces in the force array and to the fshift.
     * Convert to MD units.
     */
    for (int i = 0; i < qm->nrQMatoms_get(); i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            f[i][j]      = HARTREE_BOHR2MD*QMgrad[i][j];
         // fshift[i][j] = HARTREE_BOHR2MD*QMgrad[i][j];
        }
    }
    for (int i = 0; i < mm.nrMMatoms; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            f[i+qm->nrQMatoms_get()][j]      = HARTREE_BOHR2MD*MMgrad[i][j];
         // fshift[i+qm.nrQMatoms_get()][j] = HARTREE_BOHR2MD*MMgrad[i][j];
        }
    }
    if (qm->qmmm_variant_get() == eqmmmPME)
    {
        for (int i = 0; i < mm.nrMMatoms_full; i++)
        {
            for (int j = 0; j < DIM; j++)
            {
                f[i+qm->nrQMatoms_get()+mm.nrMMatoms][j]      = HARTREE_BOHR2MD*MMgrad_full[i][j];
             // fshift[i+qm.nrQMatoms_get()+mm.nrMMatoms][j] = HARTREE_BOHR2MD*MMgrad_full[i][j];
            }
        }
    }

    if (f_q && step % output_freq_q == 0)
    {
        fprintf(f_q, "%8d", step);
        for (int i=0; i<n; i++)
        {
            fprintf(f_q, " %8.5f", qm->QMcharges_get(i));
        }
        fprintf(f_q, "\n");
    }

    if (f_p && step % output_freq_p == 0)
    {
        fprintf(f_p, "%8d", step);
        for (int i=0; i<n; i++)
        {
            fprintf(f_p, " %8.5f", qm->pot_qmmm_get(i) + qm->pot_qmqm_get(i));
         // if (qm.qmmm_variant == eqmmmPME)
         // {
         //     fprintf(f_p, " %8.5f %8.5f %8.5f", qm.pot_qmmm[i], qm.pot_qmqm[i], qm.pot_qmmm[i] + qm.pot_qmqm[i]);
         // }
         // else
         // {
         //     fprintf(f_p, " %8.5f", qm.pot_qmmm[i]);
         // }
        }
        fprintf(f_p, "\n");
    }

    if (f_x_qm && step % output_freq_x_qm == 0)
    {
        const char periodic_system[37][3]={"XX",
            "h",                               "he",
            "li","be","b", "c", "n", "o", "f", "ne",
            "na","mg","al","si","p", "s", "cl","ar",
            "k", "ca","sc","ti","v", "cr","mn","fe","co",
            "ni","cu","zn","ga","ge","as","se","br","kr"};
        fprintf(f_x_qm, "\nQM coordinates and charges step %d\n", step);
        for (int i=0; i<n; i++) {
            fprintf(f_x_qm, "%-2s %10.5f%10.5f%10.5f %10.7f\n",
                periodic_system[qm->atomicnumberQM_get(i)],
                qm->xQM_get(i, 0) * 10., qm->xQM_get(i, 1) * 10., qm->xQM_get(i, 2) * 10.,
                qm->QMcharges_get(i));
        }
    }

    if (f_x_mm && step % output_freq_x_mm == 0)
    {
        fprintf(f_x_mm, "\nMM coordinates and charges step %d\n", step);
        for (int i=0; i<mm.nrMMatoms; i++) {
            fprintf(f_x_mm, "%10.7f%10.5f%10.5f%10.5f\n",
                mm.MMcharges[i], // CHECK -- HOW TO DO IT WITH PME / WITH CUT-OFF?
                mm.xMM[i][0] * 10., mm.xMM[i][1] * 10., mm.xMM[i][2] * 10.);
        }
    }

    sfree(QMgrad);
    sfree(MMgrad);
    if (qm->qmmm_variant_get() == eqmmmPME)
    {
        sfree(MMgrad_full);
    }

    sfree(x);
    sfree(grad);
    sfree(pot);
    sfree(potgrad);
    sfree(q);
    sfree(pot_sr);
    sfree(pot_lr);

    step++;

    return (real) QMener * HARTREE2KJ * AVOGADRO;
} /* call_dftbplus */

/* end of dftbplus sub routines */

#endif

#pragma GCC diagnostic pop

