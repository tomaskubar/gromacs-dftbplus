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

#if GMX_QMMM_NN

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
#include "gromacs/mdlib/qm_nn.h"
//#include <tensorflow/c/c_api.h>
extern "C" {
    #include "tensorflow/c/c_api.h"
}

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

void init_nn(QMMM_QMrec*       qm)
{
    //********* Read model
    //snew(qm->model, 1);

    qm->model->Graph = TF_NewGraph();
    qm->model->Status = TF_NewStatus();

    TF_SessionOptions* SessionOpts = TF_NewSessionOptions();
    TF_Buffer* RunOpts = NULL;
    
    char* saved_model_dir; 
    if ((saved_model_dir = getenv("GMX_TF_MODEL_PATH")) == nullptr)
    {
        printf("Path to model must be given in the env variable GMX_TF_MODEL_PATH\n");
        exit(-1);
    }
    else {
        printf("Read model from path %s",saved_model_dir);
    }
    const char* tags = "serve"; 
    
    int ntags = 1;
    qm->model->Session = TF_LoadSessionFromSavedModel(SessionOpts, RunOpts, saved_model_dir, &tags, ntags, qm->model->Graph, NULL, qm->model->Status);
    if(TF_GetCode(qm->model->Status) == TF_OK)
    {
        printf("TF_LoadSessionFromSavedModel OK\n");
    }
    else
    {
        printf("%s",TF_Message(qm->model->Status));
        exit(-1);
    }

    //****** Get input tensor
    qm->model->NumInputs = 1;
    qm->model->Input = (TF_Output*)malloc(sizeof(TF_Output) * qm->model->NumInputs);
    TF_Output t0 = {TF_GraphOperationByName(qm->model->Graph, "serving_default_input_1"), 0};
    
    if (t0.oper == NULL) {
        printf("ERROR: Failed TF_GraphOperationByName serving_default_input_1\n");
        exit(-1);
    }
    else {
        printf("TF_GraphOperationByName serving_default_input_1 is OK\n");
    }
    qm->model->Input[0] = t0;
    
    //********* Get Output tensor
    qm->model->NumOutputs = 1;
    qm->model->Output = (TF_Output*)malloc(sizeof(TF_Output) * qm->model->NumOutputs);
    TF_Output t2 = {TF_GraphOperationByName(qm->model->Graph, "StatefulPartitionedCall"), 0};
    
    if(t2.oper == NULL) {
        printf("ERROR: Failed TF_GraphOperationByName StatefulPartitionedCall\n");
        exit(-1);
    }
    else {
      printf("TF_GraphOperationByName StatefulPartitionedCall is OK\n");
    }
    qm->model->Output[0] = t2;

    qm->model->InputValues  = (TF_Tensor**)malloc(sizeof(TF_Tensor*)*qm->model->NumInputs);
    qm->model->OutputValues = (TF_Tensor**)malloc(sizeof(TF_Tensor*)*qm->model->NumOutputs);
    
    TF_DeleteSessionOptions(SessionOpts);
    TF_DeleteBuffer(RunOpts);
    
    /* Pass the list of QM atoms to DFTB+ */
    qm->model->nAtoms = qm->nrQMatoms_get();
    snew(qm->model->atomicNumbers, qm->model->nAtoms);
    
    // read the atomic numbers of QM atoms, and set up the lists
    for (int i=0; i<qm->model->nAtoms; i++)
    {
        qm->model->atomicNumbers[i] = qm->atomicnumberQM_get(i);
    }

    //read charges of QM atoms from file
    char* qm_charges_file; 
    if ((qm_charges_file = getenv("GMX_QMCHARGES_FILE")) == nullptr)
    {
        printf("Path to QM charges must be given in the env variable GMX_QMCHARGES_FILE\n");
        exit(-1);
    }
    else {
        printf("Read charges from path %s",qm_charges_file);
    }

    FILE* file = fopen("testinput.txt", "r");
    if (file == NULL) {
        printf("Error opening charge file.\n");
        exit(-1);
    }
    float q[qm->model->nAtoms];
    int numValues = 0;
    while (fscanf(file, "%f", &q[numValues]) == 1) {
        numValues++;
        if (numValues == qm->model->nAtoms) {
            printf("QM charges read.\n");
            break;
        }
    }
    fclose(file);
    for (int i=0; i<qm->model->nAtoms; i++)
    {
        qm->QMcharges_set(i, (real) q[i]); // sign OK
    }

    return;
} /* init_dftbplus */
void NoOpDeallocator(void* data, size_t a, void* b) {}
real call_nn(QMMM_rec*         qr,
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

    int n = qm->nrQMatoms_get();

    double *x, *grad, *pot, *potgrad, *q; // real instead of rvec, to help pass data to fortran
    real *pot_sr = nullptr, *pot_lr = nullptr;
    rvec *QMgrad = nullptr, *MMgrad = nullptr, *MMgrad_full = nullptr;

    //snew(x, (n,3));
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

    float data[1][n][4];
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<DIM; j++)
        {
            data[0][i][j] = qm->xQM_get(i,j) / BOHR2NM; // to bohr units for DFTB+, Manus Model also expects Bohr
        }
    }

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

    // for (int i=0; i<n; i++)
    // {
    //     q[i] = 0.;
    // }

    //concat x reshaped with pot MANU to data
    for (int i=0; i<n; i++) {
        // for (int j=0; j<3; j++) {
        //     data[0][i][j]=x[i][j];  //x,y,z in Bohr
        // }
        data[0][i][3]=pot[i];   //ESP in Bohr
    }

    /* DFTB+ calculation itself */
    wallcycle_start(wcycle, ewcQM);

    int ndims = 3;
    int64_t dims[] = {1,n,4};
    
    int ndata = sizeof(x); 
    TF_Tensor* int_tensor = TF_NewTensor(TF_FLOAT, dims, ndims, data, ndata, &NoOpDeallocator, 0);
    
    if (int_tensor != NULL){
        printf("TF_NewTensor is OK\n");
    }
    else {
      printf("ERROR: Failed TF_NewTensor\n");
    }
    qm->model->InputValues[0] = int_tensor;

    // Run the Session
    TF_SessionRun(qm->model->Session, NULL, qm->model->Input, qm->model->InputValues, qm->model->NumInputs, qm->model->Output, qm->model->OutputValues, qm->model->NumOutputs, NULL, 0,NULL , qm->model->Status);
    
    if(TF_GetCode(qm->model->Status) == TF_OK) {
      printf("Session is OK\n");
    }
    else {
      printf("%s",TF_Message(qm->model->Status));
    }
    float* predictions = (float*)TF_TensorData(qm->model->OutputValues[0]);

 // dftbp_set_coords(qm->dpcalc, x); // unit OK
 // dftbp_set_external_potential(qm->dpcalc, pot, potgrad); // unit and sign OK
 // dftbp_get_energy(qm->dpcalc, &QMener); // unit OK
 // dftbp_get_gross_charges(qm->dpcalc, q);
 // for (int i=0; i<n; i++)
 //     printf("%d %6.3f\n", i+1, q[i]);
 // dftbp_get_gradients(qm->dpcalc, grad);
    wallcycle_stop(wcycle, ewcQM);
    QMener = predictions[0];
    /* Save the gradient on the QM atoms */
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<3; j++)
        {
            QMgrad[i][j] = (real) - predictions[3*i+j+1]; // negative of force -- sign OK
        }
    }

 // /* Print the QM pure gradient */
 // for (int i=0; i<n; i++)
 // {
 //     printf("GRAD QM %d: %8.2f %8.2f %8.2f\n", i+1,
 //         QMgrad[i][XX] * HARTREE_BOHR2MD, QMgrad[i][YY] * HARTREE_BOHR2MD, QMgrad[i][ZZ] * HARTREE_BOHR2MD);
 // }

    /* Save the QM charges */
    // for (int i=0; i<n; i++)
    // {
    //     qm->QMcharges_set(i, (real) q[i]); // sign OK
    //  // printf("CHECK CHARGE QM[%d] = %6.3f\n", i+1, qm->QMcharges[i]);
    // }

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
} /* call_nn */

/* end of dftbplus sub routines */

#endif

#pragma GCC diagnostic pop

