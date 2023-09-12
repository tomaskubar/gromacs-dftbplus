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

void NoOpDeallocator(void* data, size_t a, void* b);
void prepare_nn_inputs(QMMM_QMrec* qm);

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
    snew(qm->model, 1);
    //printf("Hello from TensorFlow C library version %s\n", TF_Version());
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
        printf("Read model from path %s\n",saved_model_dir);
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
    qm->model->NumInputs = 11;
    qm->model->Input = (TF_Output*)malloc(sizeof(TF_Output) * qm->model->NumInputs);
    
    const char* base_name = "serving_default_args_0";
    char operation_name[100];

    for (int i=0; i<qm->model->NumInputs; i++) {
        if (i == 0) {
            strcpy(operation_name, base_name);
        }
        else {
            snprintf(operation_name, sizeof(operation_name), "%s_%d", base_name, i);
        }
        TF_Output t = {TF_GraphOperationByName(qm->model->Graph, operation_name), 0};

        if(t.oper == NULL) {
            printf("ERROR: Failed TF_GraphOperationByName %s\n", operation_name);
            exit(-1);
        }
        else {
            printf("TF_GraphOperationByName %s is OK\n", operation_name);
            qm->model->Input[i] = t;
        }
    }
    
    //********* Get Output tensor
    qm->model->NumOutputs = 3;
    qm->model->Output = (TF_Output*)malloc(sizeof(TF_Output) * qm->model->NumOutputs);
    TF_Output charge_output = {TF_GraphOperationByName(qm->model->Graph, "StatefulPartitionedCall"), 0};

    if(charge_output.oper == NULL) {
        printf("ERROR: Failed TF_GraphOperationByName StatefulPartitionedCall\n");
        exit(-1);
    }
    else
        printf("TF_GraphOperationByName StatefulPartitionedCall is OK\n");
    
    qm->model->Output[0] = charge_output;

    TF_Output energy_output = {TF_GraphOperationByName(qm->model->Graph, "StatefulPartitionedCall"), 1};

    if(energy_output.oper == NULL) {
        printf("ERROR: Failed TF_GraphOperationByName StatefulPartitionedCall\n");
        exit(-1);
    }
    else {
      printf("TF_GraphOperationByName StatefulPartitionedCall is OK\n");
    }
    qm->model->Output[1] = energy_output;

    TF_Output force_output = {TF_GraphOperationByName(qm->model->Graph, "StatefulPartitionedCall"), 2};

    if(force_output.oper == NULL) {
        printf("ERROR: Failed TF_GraphOperationByName StatefulPartitionedCall\n");
        exit(-1);
    }
    else {
      printf("TF_GraphOperationByName StatefulPartitionedCall is OK\n");
    }
    qm->model->Output[2] = force_output;

    // Allocate data for inputs & outputs
    qm->model->InputValues  = (TF_Tensor**)malloc(sizeof(TF_Tensor*)*qm->model->NumInputs);
    qm->model->OutputValues = (TF_Tensor**)malloc(sizeof(TF_Tensor*)*qm->model->NumOutputs);
    
    TF_DeleteSessionOptions(SessionOpts);
    TF_DeleteBuffer(RunOpts);

    return;
} /* init_nn */

void NoOpDeallocator(void* data, size_t a, void* b) {(void) b; if (a > 0) {(void) data;}} //Nonsense Deallocator to make compiler happy
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

    double *grad, *pot, *potgrad, *q; // real instead of rvec, to help pass data to fortran
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

    snew(QMgrad, n);

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

    /* NN call itself */
    wallcycle_start(wcycle, ewcQM);

    int nAtoms = qm->nrQMatoms_get();

    // node number flat values 0
    int ndims_0 = 1;
    int64_t dims_0[] = {nAtoms};
    int64_t data_0[nAtoms] = {};
    for (int i=0; i<nAtoms; i++)
    {
        data_0[i] = qm->atomicnumberQM_get(i);
        printf("CHECK Atomic numbers %d %ld\n", i, data_0[i]);
    }
    int ndata_0 = sizeof(data_0);
    qm->model->InputValues[0] = TF_NewTensor(TF_INT64, dims_0, ndims_0, data_0, ndata_0, &NoOpDeallocator, data_0);

    // node number row splits 1
    int ndims_1 = 1;
    int64_t dims_1[] = {2};
    int64_t data_1[2] = {0, nAtoms};
    int ndata_1 = sizeof(data_1); 
    qm->model->InputValues[1] = TF_NewTensor(TF_INT64, dims_1, ndims_1, data_1, ndata_1, &NoOpDeallocator, data_1);

    // node coordinates flat values 2
    int ndims_2 = 2;
    int64_t dims_2[] = {nAtoms, 3};
    float data_2[nAtoms][3] = {};
    for (int i=0; i<nAtoms; i++)
    {
        for (int j=0; j<3; j++)
        {
            data_2[i][j] = qm->xQM_get(i,j) / BOHR2NM; // from nm to bohr for Lukas model
        }
        printf("CHECK COORD QM[%d] = %6.3f %6.3f %6.3f\n", i+1, data_2[i][0], data_2[i][1], data_2[i][2]);
    }
    int ndata_2 = sizeof(data_2); 
    qm->model->InputValues[2] = TF_NewTensor(TF_FLOAT, dims_2, ndims_2, data_2, ndata_2, &NoOpDeallocator, data_2);

    // node coordinates row splits 3
    qm->model->InputValues[3] = qm->model->InputValues[1];

    // edge indces flat values 4
    int nEdgeCombinations = nAtoms*(nAtoms-1); // Cheap Combinations nCr(nAtoms 2)
    int ndims_4 = 2;
    int64_t dims_4[] = {nEdgeCombinations, 2};
    int64_t data_4[nEdgeCombinations][2] = {};
    int index = 0;
    for (int i=0; i<nAtoms; i++)
    {   
        for (int j=0; j<nAtoms; j++)
        {
            if (i==j) {
                continue;
            }
            data_4[index][0] = i;
            data_4[index][1] = j;
            index++;
        }
    }
    printf("CHECK Bonds %d %d\n", index, nEdgeCombinations);
    int ndata_4 = sizeof(data_4); 
    qm->model->InputValues[4] = TF_NewTensor(TF_INT64, dims_4, ndims_4, data_4, ndata_4, &NoOpDeallocator, data_4);

    // edge indces row splits 5
    int ndims_5 = 1;
    int64_t dims_5[] = {2};
    int64_t data_5[2] = {0, nEdgeCombinations};
    int ndata_5 = sizeof(data_5); 
    qm->model->InputValues[5] = TF_NewTensor(TF_INT64, dims_5, ndims_5, data_5, ndata_5, &NoOpDeallocator, data_5);

    // angle indces flat values 6
    int nAngleCombinations;
    if (nAtoms > 2) {
        nAngleCombinations = nAtoms*(nAtoms-1)*(nAtoms-2); // Cheap Combinations nAtoms*nCr(nAtoms-1 2)
    }
    else {
        nAngleCombinations = 0;
    }

    int ndims_6 = 2;
    int64_t dims_6[] = {nAngleCombinations, 3};
    int64_t data_6[nAngleCombinations][3] = {};

    index = 0;
    for (int i=0; i<nAtoms; i++)
    {
        for (int j=0; j<nAtoms; j++)
        {
            if (i == j) {
                continue;
            }
            for (int k=0; k<nAtoms; k++)
            {
                if ((i == k) || (j == k))  {
                    continue;
                }
                data_6[index][0] = i;
                data_6[index][1] = j;
                data_6[index][2] = k; 
                index++;
            }
        }
    }
    printf("CHECK Angles %d %d\n", index, nAngleCombinations);
    int ndata_6 = sizeof(data_6); 
    qm->model->InputValues[6] = TF_NewTensor(TF_INT64, dims_6, ndims_6, data_6, ndata_6, &NoOpDeallocator, data_6);

    // angle indces row splits 7
    int ndims_7 = 1;
    int64_t dims_7[] = {2};
    int64_t data_7[2] = {0, nAngleCombinations};
    int ndata_7 = sizeof(data_7); 
    qm->model->InputValues[7] = TF_NewTensor(TF_INT64, dims_7, ndims_7, data_7, ndata_7, &NoOpDeallocator, data_7);

    // total charge 8
    int ndims_8 = 1;
    int64_t dims_8[] = {1};
    float data_8[1] = {(float)qm->QMcharge_get()};
    printf("CHECK Total Charge %f\n", (float)qm->QMcharge_get());
    int ndata_8 = sizeof(data_8);
    qm->model->InputValues[8] = TF_NewTensor(TF_FLOAT, dims_8, ndims_8, data_8, ndata_8, &NoOpDeallocator, data_8);

    // esp flat values 9
    int ndims_9 = 1;
    int64_t dims_9[] = {nAtoms};
    float data_9[nAtoms] = {};
    for (int i=0; i<nAtoms; i++)
    {
        data_9[i] = qm->pot_qmmm_get(i); // in volt units
        printf("CHECK ESP QM[%d] = %6.3f\n", i+1, data_9[i]);
    }
    int ndata_9 = sizeof(data_9); 
    qm->model->InputValues[9] = TF_NewTensor(TF_FLOAT, dims_9, ndims_9, data_9, ndata_9, &NoOpDeallocator, data_9);

    // esp row splits 10
    qm->model->InputValues[10] = qm->model->InputValues[1];

    // Run the Session
    TF_SessionRun(qm->model->Session, NULL,
        qm->model->Input, qm->model->InputValues, qm->model->NumInputs,
        qm->model->Output, qm->model->OutputValues, qm->model->NumOutputs,
        NULL, 0, NULL, qm->model->Status);
    
    if(!(TF_GetCode(qm->model->Status) == TF_OK)) {
      printf("%s",TF_Message(qm->model->Status));
    }
    float* charge_predictions = (float*)TF_TensorData(qm->model->OutputValues[0]); // in e-
    float* energy_predictions = (float*)TF_TensorData(qm->model->OutputValues[1]); // in Hartree
    float* grad_predictions = (float*)TF_TensorData(qm->model->OutputValues[2]); // in Hartree/Bohr

    wallcycle_stop(wcycle, ewcQM);

    /* Save the QM charges */
    for (int i=0; i<n; i++)
    {
        qm->QMcharges_set(i, (real) charge_predictions[i]); // sign OK
        // printf("CHECK CHARGE QM[%d] = %6.3f\n", i+1, qm->QMcharges_get(i));
    }
    
    QMener = energy_predictions[0];
    // printf("CHECK QM Energy = %6.3f\n", QMener);

    /* Save the gradient on the QM atoms */
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<3; j++)
        {
            QMgrad[i][j] = (real) grad_predictions[3*i+j]; // negative of force 
        }
        // printf("CHECK QM Grad[%d] = %6.3f %6.3f %6.3f\n", i+1, QMgrad[i][0], QMgrad[i][1], QMgrad[i][2]);
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

    rvec *partgrad;
    snew(partgrad, qm->nrQMatoms_get());
    qr->gradient_QM_MM(cr, nrnb, wcycle, (qm->qmmm_variant_get() == eqmmmPME ? *qr->pmedata : nullptr),
                   qm->qmmm_variant_get(), partgrad, MMgrad, MMgrad_full);
    // for (int i=0; i<n; i++) // No need to adjust the QMgrad for the MM contribution, handled by the model
    // {
    //     rvec_inc(QMgrad[i], partgrad[i]); // sign OK
    //  //   printf("GRAD QM FULL %d: %8.2f %8.2f %8.2f\n", i+1,
    //  //       QMgrad[i][0] , QMgrad[i][1] , QMgrad[i][2] );
    // }
    sfree(partgrad);
    /* Put the QMMM forces in the force array and to the fshift.
     * Convert to MD units.
     */
    for (int i = 0; i < qm->nrQMatoms_get(); i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            f[i][j]      = -HARTREE_BOHR2MD*QMgrad[i][j];
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
    // for (int i = 0; i < n; i++) {
    // for (int j = 0; j < 3; j++) {
    //     printf("%f ", f[i][j]);
    // }
    // printf("\n");
    // }
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

    sfree(grad);
    sfree(pot);
    sfree(potgrad);
    sfree(q);
    sfree(pot_sr);
    sfree(pot_lr);

    step++;

    return (real) QMener * HARTREE2KJ * AVOGADRO;
} /* call_nn */

// responsible for putting correct NN inputs into qm->model->Input_values
// Couldn't figure out a proper way to keep it allocated as a seperate function
void prepare_nn_inputs(QMMM_QMrec* qm)
{
    int nAtoms = qm->nrQMatoms_get();

    // node number flat values 0
    int ndims_0 = 1;
    int64_t dims_0[] = {nAtoms};
    int64_t data_0[nAtoms] = {};
    for (int i=0; i<nAtoms; i++)
    {
        data_0[i] = qm->atomicnumberQM_get(i);
    }
    int ndata_0 = sizeof(data_0);
    qm->model->InputValues[0] = TF_NewTensor(TF_INT64, dims_0, ndims_0, data_0, ndata_0, &NoOpDeallocator, data_0);

    // node number row splits 1
    int ndims_1 = 1;
    int64_t dims_1[] = {2};
    int64_t data_1[2] = {0, nAtoms};
    int ndata_1 = sizeof(data_1); 
    qm->model->InputValues[1] = TF_NewTensor(TF_INT64, dims_1, ndims_1, data_1, ndata_1, &NoOpDeallocator, data_1);

    // node coordinates flat values 2
    int ndims_2 = 2;
    int64_t dims_2[] = {nAtoms, 3};
    float data_2[nAtoms][3] = {};
    for (int i=0; i<nAtoms; i++)
    {
        for (int j=0; j<3; j++)
        {
            data_2[i][j] = qm->xQM_get(i,j) / BOHR2NM; // to bohr units for Lukas model
        }
    }
    int ndata_2 = sizeof(data_2); 
    qm->model->InputValues[2] = TF_NewTensor(TF_FLOAT, dims_2, ndims_2, data_2, ndata_2, &NoOpDeallocator, data_2);

    // node coordinates row splits 3
    qm->model->InputValues[3] = qm->model->InputValues[1];

    // edge indces flat values 4
    int nEdgeCombinations = nAtoms*(nAtoms-1)/2; // Cheap Binomialcoeffient (nAtoms 2)
    int ndims_4 = 2;
    int64_t dims_4[] = {nEdgeCombinations, 2};
    int64_t data_4[nEdgeCombinations][2] = {};
    int index = 0;
    for (int i=0; i<nAtoms; i++)
    {   
        for (int j=i+1; j<nAtoms; j++)
        {
            data_4[index][0] = i;
            data_4[index][1] = j;
            index++;
        }
    }
    int ndata_4 = sizeof(data_4); 
    qm->model->InputValues[4] = TF_NewTensor(TF_INT64, dims_4, ndims_4, data_4, ndata_4, &NoOpDeallocator, data_4);

    // edge indces row splits 5
    int ndims_5 = 1;
    int64_t dims_5[] = {2};
    int64_t data_5[2] = {0, nEdgeCombinations};
    int ndata_5 = sizeof(data_5); 
    qm->model->InputValues[5] = TF_NewTensor(TF_INT64, dims_5, ndims_5, data_5, ndata_5, &NoOpDeallocator, data_5);

    // angle indces flat values 6
    int nAngleCombinations = nAtoms*(nAtoms-1)*(nAtoms-2)/6; // Cheap Binomialcoeffient (nAtoms 3)
    int ndims_6 = 2;
    int64_t dims_6[] = {nAngleCombinations, 3};
    int64_t data_6[nAngleCombinations][3] = {};

    index = 0;
    for (int i=0; i<nAtoms; i++)
    {
        for (int j=i+1; j<nAtoms; j++)
        {
            for (int k=j+1; k<nAtoms; k++)
            {
                data_6[index][0] = i;
                data_6[index][1] = j;
                data_6[index][2] = k;
                index++;
            }
        }
    }
    int ndata_6 = sizeof(data_6); 
    qm->model->InputValues[6] = TF_NewTensor(TF_INT64, dims_6, ndims_6, data_6, ndata_6, &NoOpDeallocator, data_6);

    // angle indces row splits 7
    int ndims_7 = 1;
    int64_t dims_7[] = {2};
    int64_t data_7[2] = {0, nAngleCombinations};
    int ndata_7 = sizeof(data_7); 
    qm->model->InputValues[7] = TF_NewTensor(TF_INT64, dims_7, ndims_7, data_7, ndata_7, &NoOpDeallocator, data_7);

    // total charge 8
    int ndims_8 = 1;
    int64_t dims_8[] = {1};
    float data_8[1] = {(float)qm->QMcharge_get()};
    int ndata_8 = sizeof(data_8);
    qm->model->InputValues[8] = TF_NewTensor(TF_FLOAT, dims_8, ndims_8, data_8, ndata_8, &NoOpDeallocator, data_8);

    // esp flat values 9
    int ndims_9 = 1;
    int64_t dims_9[] = {nAtoms};
    float data_9[nAtoms] = {};
    for (int i=0; i<nAtoms; i++)
    {
        qm->pot_qmmm_get(i); // in volt units
    }
    int ndata_9 = sizeof(data_9); 
    qm->model->InputValues[9] = TF_NewTensor(TF_FLOAT, dims_9, ndims_9, data_9, ndata_9, &NoOpDeallocator, data_9);

    // esp row splits 10
    qm->model->InputValues[10] = qm->model->InputValues[1];


    // Print content of data for debugging
    // for (int i=0; i<nAtoms; i++)
    // {
    //     printf("%ld\n", data_0[i]); 
    // }
    // for (int i=0; i<2; i++)
    // {
    //     printf("%ld\n", data_1[i]); 
    // }
    // for (int i=0; i<nAtoms; i++)
    // {
    //     printf("%f,%f,%f\n", data_2[i][0], data_2[i][1], data_2[i][1]); 
    // }
    // for (int i=0; i<nEdgeCombinations; i++)
    // {
    //     printf("%ld,%ld\n", data_4[i][0],  data_4[i][1]); 
    // }
    // for (int i=0; i<nAngleCombinations; i++)
    // {
    //     printf("%ld,%ld,%ld\n", data_6[i][0],  data_6[i][1], data_6[i][2]); 
    // }

    return;
} //prepare_nn_inputs

/* end of nn sub routines */

#endif

#pragma GCC diagnostic pop

