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
#include "gromacs/mdlib/qm_dftbplus.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/timing/cyclecounter.h"
#include "gromacs/timing/walltime_accounting.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/pbcutil/pbc.h"

#include "gromacs/transfer/qmmm.h"
//#include "gromacs/transfer/qm_dftbplus.h"

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

/* Fill up the context (status structure) with all the relevant data
 */
void initialize_context_transfer(Context*           cont,
                                 int                nrQMatoms,
                                 int                qmmm_variant,
                                 QMMM_rec_transfer* qr_in,
                                 const real         rcoulomb,
                                 const real         ewald_rtol,     
                                 const t_commrec*   cr_in)
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
  cont->n = nrQMatoms;
  printf("cont->n = %d\n", cont->n);
  if (cont->pme)
  {
      cont->cr           = cr_in;
      cont->qr           = qr_in;
   // cont->wcycle       = wcycle_in;
      cont->rcoul        = rcoulomb;
      cont->ewaldcoeff_q = calc_ewaldcoeff_q(rcoulomb, ewald_rtol);
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
void calcQMextPotPME_transfer(Context *cont, double *q, double *extpot)
{
  int n = cont->n;

  if (cont->pme)
  {
      /* PERFORM THE REAL CALCULATION */
      real *extpot_real;
      snew(extpot_real, n);
      for (int i=0; i<n; i++)
      {
          cont->qr->qm->QMcharges_set(i, (real) -q[i]); // check sign TODO
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
extern "C" void calcqmextpot_transfer(void *refptr, double *q, double *extpot)
{
  Context *cont = (Context *) refptr;
  calcQMextPotPME_transfer(cont, q, extpot);
  return;
}

/* The wrapper for the calculation of grdient of external potential
 *   due to the periodic images of QM atoms with PME.
 * This is function does not calculate anything,
 *   because it is not needed in this DFTB+Gromacs implementation.
 */
extern "C" void calcqmextpotgrad_transfer(void *refptr, double *q, double *extpotgrad)
{
  /* this will never be used */
  (void) refptr;
  (void) q;
  (void) extpotgrad;
}


/* DFTBPLUS interface routines */

void init_dftbplus_transfer(QMMM_rec_transfer*   qr,
                            const real           rcoulomb,
                            const real           ewald_rtol,
                            const t_commrec*     cr)
{
    QMMM_QMrec_transfer* qm = qr->qm.get();
    // variables
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

    // TODO: Replace "calculator" with "qm->dpcalc", which needs to be allocated first
    //       Otherwise, there may be a problem with "calculator" begin static!

    /* This structure will be passed through DFTB+
     *   into the Gromacs calculator calcQMextPot
     */

    /* Fill up the context with all the relevant data! */
    snew(qm->dftbContext, 1);
    initialize_context_transfer(qm->dftbContext, qm->nrQMatoms_get(), qm->qmmm_variant_get(), qr, rcoulomb, ewald_rtol, cr); //, wcycle);

    /* Initialize the DFTB+ calculator */
    dftbp_init(&calculator, "dftb_in.out");
    printf("DFTB+ calculator has been created!\n");

    /* Parse the input file and store the input-tree */
    dftbp_get_input_from_file(&calculator, "dftb_in.hsd", &input);
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

    snew(qm->dpcalc, 1);
    qm->dpcalc = &calculator;

    /* Register the callback functions which calculate
     * the external potential and its gradient
     */
    dftbp_register_ext_pot_generator(qm->dpcalc,
                                     qm->dftbContext,
                                     calcqmextpot_transfer,
                                     calcqmextpotgrad_transfer);

    /* Alternatively - modify the geometry/elements information here,
     * based on what we have in the Gromacs topology!
     */

    return;
} /* init_dftbplus */

real call_dftbplus_transfer(QMMM_rec_transfer*   qr,
                            const t_commrec*     cr,
                            rvec*                f,
		                    t_nrnb*              nrnb,
                            gmx_wallcycle_t      wcycle)
{
    double QMener;

    qr->qm->dftbContext->nrnb = nrnb;
    int n = qr->qm->nrQMatoms_get();

    double *x, *grad, *pot, *potgrad, *q; // real instead of rvec, to help pass data to fortran
    real *pot_sr = nullptr, *pot_lr = nullptr;
    rvec *QMgrad = nullptr;

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

    snew(QMgrad, qr->qm->nrQMatoms_get());

    for (int i=0; i<n; i++)
    {
        for (int j=0; j<DIM; j++)
        {
            x[3*i+j] = qr->qm->xQM_get(i,j) / BOHR2NM; // to bohr units for DFTB+
        }
    }

    /* calculate the QM/MM electrostatics beforehand with Gromacs!
     * with cut-off treatment, this will be the entire QM/MM
     * with PME, this will only include MM atoms,
     *    and the contribution from QM atoms will be added in every SCC iteration
     */
    qr->calculate_SR_QM_MM(qr->qm->qmmm_variant_get(), pot_sr);
    if (qr->qm->qmmm_variant_get() == eqmmmPME)
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
        qr->qm->pot_qmmm_set(j, (double) - pot[j] * HARTREE_TO_EV); // in volt units
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

    for (int i=0; i<n; i++)
    {
        q[i] = 0.;
    }

    /* DFTB+ calculation itself */
    wallcycle_start(wcycle, ewcQM);
    dftbp_set_coords(qr->qm->dpcalc, x); // unit OK
    dftbp_set_external_potential(qr->qm->dpcalc, pot, potgrad); // unit and sign OK
    dftbp_get_energy(qr->qm->dpcalc, &QMener); // unit OK
    dftbp_get_gross_charges(qr->qm->dpcalc, q);
 // for (int i=0; i<n; i++)
 //     printf("%d %6.3f\n", i+1, q[i]);
    dftbp_get_gradients(qr->qm->dpcalc, grad);
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
        qr->qm->QMcharges_set(i, (real) q[i]); // sign OK
     // printf("CHECK CHARGE QM[%d] = %6.3f\n", i+1, qr->qm->QMcharges[i]);
    }

    /* Put the QMMM forces in the force array and to the fshift.
     * Convert to MD units.
     */
    for (int i = 0; i < qr->qm->nrQMatoms_get(); i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            f[i][j]      = HARTREE_BOHR2MD*QMgrad[i][j];
        }
    }


    sfree(QMgrad);

    sfree(x);
    sfree(grad);
    sfree(pot);
    sfree(potgrad);
    sfree(q);
    sfree(pot_sr);
    sfree(pot_lr);

    return (real) QMener * HARTREE2KJ * AVOGADRO;
} /* call_dftbplus */

/* end of dftbplus sub routines */

