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
#include <vector>
#include <string>

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
                        // const real        rcoul_in,
                        // const real        ewaldcoeff_q_in)
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
    // int n = cont->fr->qr->qm[0]->nrQMatoms;
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

void init_nn(QMMM_QMrec* qm)
{
    // Find amount of models, max 10
    char* N_TF_MODELS; // Pointer to env variable for comparison with nullpointer
    int n_models; 
    if ((N_TF_MODELS = getenv("GMX_N_TF_MODELS")) == nullptr)
    {
        printf("Amount of models to use for adaptive sampling must be given via the env variable GMX_N_TF_MODELS\n");
        exit(-1);
    }
    else {
        n_models = std::stoi(N_TF_MODELS);
        if (n_models > 10)
        {
            printf("A maximum of 10 models is supported\n");
            exit(-1);
        }
        printf("Using %i models \n", n_models);
    }

    // Find model prefixes and pathes
    char* model_path_prefix;
    std::vector<std::string> model_paths;
    char model_path[1000];
    if ((model_path_prefix = getenv("GMX_TF_MODEL_PATH_PREFIX")) == nullptr)
    {
        printf("Path prefix to models must be given in the env variable GMX_TF_MODEL_PATH_PREFIX\n");
        exit(-1);
    }
    else {
        printf("Will read models from ");
        for (int i=0; i<n_models; i++) {
            snprintf(model_path, sizeof(model_path), "%s%d", model_path_prefix, i);
            printf("%s ", model_path);
        }
        printf("\n");
    }

    //********* Read models
    for (int model_idx=0; model_idx<n_models; model_idx++) {
        snprintf(model_path, sizeof(model_path), "%s%d", model_path_prefix, model_idx);

        snew(qm->models[model_idx], 1);
        qm->models[model_idx]->Graph = TF_NewGraph();
        qm->models[model_idx]->Status = TF_NewStatus();

        TF_SessionOptions* SessionOpts = TF_NewSessionOptions();
        TF_Buffer* RunOpts = NULL;
        
        const char* tags = "serve";
        
        int ntags = 1;
        qm->models[model_idx]->Session = TF_LoadSessionFromSavedModel(SessionOpts, RunOpts, model_path, &tags, ntags, qm->models[model_idx]->Graph, NULL, qm->models[model_idx]->Status);
        if(TF_GetCode(qm->models[model_idx]->Status) == TF_OK)
        {
            printf("TF_LoadSessionFromSavedModel OK\n");
        }
        else
        {
            printf("%s",TF_Message(qm->models[model_idx]->Status));
            exit(-1);
        }

        //****** Get input tensor
        qm->models[model_idx]->NumInputs = 13;
        qm->models[model_idx]->Input = (TF_Output*)malloc(sizeof(TF_Output) * qm->models[model_idx]->NumInputs);
        
        const char* base_name = "serving_default_args_0";
        char operation_name[100];

        for (int i=0; i<qm->models[model_idx]->NumInputs; i++) {
            if (i == 0) {
                strcpy(operation_name, base_name);
            }
            else {
                snprintf(operation_name, sizeof(operation_name), "%s_%d", base_name, i);
            }
            TF_Output t = {TF_GraphOperationByName(qm->models[model_idx]->Graph, operation_name), 0};

            if(t.oper == NULL) {
                printf("ERROR: Failed TF_GraphOperationByName %s\n", operation_name);
                exit(-1);
            }
            else {
                printf("TF_GraphOperationByName %s is OK\n", operation_name);
                qm->models[model_idx]->Input[i] = t;
            }
        }
        
        //********* Get Output tensor
        qm->models[model_idx]->NumOutputs = 3;
        qm->models[model_idx]->Output = (TF_Output*)malloc(sizeof(TF_Output) * qm->models[model_idx]->NumOutputs);
        TF_Output charge_output = {TF_GraphOperationByName(qm->models[model_idx]->Graph, "StatefulPartitionedCall"), 0};

        if(charge_output.oper == NULL) {
            printf("ERROR: Failed TF_GraphOperationByName StatefulPartitionedCall\n");
            exit(-1);
        }
        else
            printf("TF_GraphOperationByName StatefulPartitionedCall is OK\n");
        
        qm->models[model_idx]->Output[0] = charge_output;

        TF_Output energy_output = {TF_GraphOperationByName(qm->models[model_idx]->Graph, "StatefulPartitionedCall"), 1};

        if(energy_output.oper == NULL) {
            printf("ERROR: Failed TF_GraphOperationByName StatefulPartitionedCall\n");
            exit(-1);
        }
        else {
        printf("TF_GraphOperationByName StatefulPartitionedCall is OK\n");
        }
        qm->models[model_idx]->Output[1] = energy_output;

        TF_Output force_output = {TF_GraphOperationByName(qm->models[model_idx]->Graph, "StatefulPartitionedCall"), 2};

        if(force_output.oper == NULL) {
            printf("ERROR: Failed TF_GraphOperationByName StatefulPartitionedCall\n");
            exit(-1);
        }
        else {
        printf("TF_GraphOperationByName StatefulPartitionedCall is OK\n");
        }
        qm->models[model_idx]->Output[2] = force_output;

        // Allocate data for inputs & outputs
        qm->models[model_idx]->InputValues  = (TF_Tensor**)malloc(sizeof(TF_Tensor*)*qm->models[model_idx]->NumInputs);
        qm->models[model_idx]->OutputValues = (TF_Tensor**)malloc(sizeof(TF_Tensor*)*qm->models[model_idx]->NumOutputs);
        
        TF_DeleteSessionOptions(SessionOpts);
        TF_DeleteBuffer(RunOpts);
    }

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

    static float energy_prediction_std_threshold; // energy threshold deviation
    static float force_prediction_std_threshold; // force threshold deviation
    static int nn_eval_freq; // evaluation frequency of additional models

    double QMener;
    // bool lPme = (qm->qmmm_variant == eqmmmPME);
    (void) fshift;

    int n = qm->nrQMatoms_get();

    double *grad, *pot, *potgrad, *q; // real instead of rvec, to help pass data to fortran
    real *pot_sr = nullptr, *pot_lr = nullptr;
    rvec *QMgrad = nullptr, *MMgrad = nullptr, *MMgrad_full = nullptr;
    rvec *ESPgrad = nullptr, *ESPgrad_full = nullptr;

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

        // Find threshold deviation, energy threshold has priority
        if (((env = getenv("GMX_ENERGY_PREDICTION_STD_THRESHOLD")) == nullptr) && ((env = getenv("GMX_FORCE_PREDICTION_STD_THRESHOLD")) == nullptr))
        {
            printf("Deviation as threshold for significant standard deviation can be given via the env variable GMX_ENERGY_PREDICTION_STD_THRESHOLD or GMX_FORCE_PREDICTION_STD_THRESHOLD to do adaptive sampling\n");
            energy_prediction_std_threshold = -1;
            force_prediction_std_threshold = -1;        
        }
        else if ((env = getenv("GMX_ENERGY_PREDICTION_STD_THRESHOLD")) != nullptr)
        {
            energy_prediction_std_threshold = std::stof(env);
            force_prediction_std_threshold = -1;
            printf("Using %6.3f as standard deviation threshold of the energy prediction for adaptive sampling\n", energy_prediction_std_threshold);
        }
        else if ((env = getenv("GMX_FORCE_PREDICTION_STD_THRESHOLD")) != nullptr)
        {
            force_prediction_std_threshold = std::stof(env);
            energy_prediction_std_threshold = -1;
            printf("Using %6.3f as standard deviation threshold of the force prediction for adaptive sampling\n", force_prediction_std_threshold);
        }
        // Find and set evaluation frequency of additional models(for adaptive sampling as well as safeguard during simulations), default 1
        if ((env = getenv("GMX_NN_EVAL_FREQ")) != nullptr)
        {
            nn_eval_freq = std::stoi(env);
            printf("Using %i as evaluation frequency of additional models\n", nn_eval_freq);
        }
        else
        {
            nn_eval_freq = 1;
            printf("Using %i as default evaluation frequency of additional models\n", nn_eval_freq);
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

    int n_models = std::stoi(getenv("GMX_N_TF_MODELS"));
    int nAtoms = qm->nrQMatoms_get();

    // Decide how many models to use
    int n_active_models;
    if (step % nn_eval_freq == 0)
    {
        n_active_models = n_models;
    }
    else
    {
        n_active_models = 1;
    }
    
    // node number flat values 0
    int ndims_0 = 1;
    int64_t dims_0[] = {nAtoms};
    int64_t data_0[nAtoms] = {};
    for (int i=0; i<nAtoms; i++)
    {
        data_0[i] = qm->atomicnumberQM_get(i);
    }
    int ndata_0 = sizeof(data_0);
    for (int model_idx=0; model_idx<n_active_models; model_idx++) {
        qm->models[model_idx]->InputValues[0] = TF_NewTensor(TF_INT64, dims_0, ndims_0, data_0, ndata_0, &NoOpDeallocator, nullptr);
    }

    // node number row splits 1
    int ndims_1 = 1;
    int64_t dims_1[] = {2};
    int64_t data_1[2] = {0, nAtoms};
    int ndata_1 = sizeof(data_1);
    for (int model_idx=0; model_idx<n_active_models; model_idx++) {
        qm->models[model_idx]->InputValues[1] = TF_NewTensor(TF_INT64, dims_1, ndims_1, data_1, ndata_1, &NoOpDeallocator, nullptr);
    }

    // node coordinates flat values 2
    int ndims_2 = 2;
    int64_t dims_2[] = {nAtoms, 3};
    float data_2[nAtoms][3] = {};
    for (int i=0; i<nAtoms; i++)
    {
        for (int j=0; j<3; j++)
        {
            data_2[i][j] = qm->xQM_get(i,j) / BOHR2NM; // from nm to bohr for NN model
        }
    }
    int ndata_2 = sizeof(data_2);
    for (int model_idx=0; model_idx<n_active_models; model_idx++) {
        qm->models[model_idx]->InputValues[2] = TF_NewTensor(TF_FLOAT, dims_2, ndims_2, data_2, ndata_2, &NoOpDeallocator, nullptr);
    }

    // node coordinates row splits 3
    int ndims_3 = 1;
    int64_t dims_3[] = {2};
    int64_t data_3[2] = {0, nAtoms};
    int ndata_3 = sizeof(data_3);
    for (int model_idx=0; model_idx<n_active_models; model_idx++) {
        qm->models[model_idx]->InputValues[3] = TF_NewTensor(TF_INT64, dims_3, ndims_3, data_3, ndata_3, &NoOpDeallocator, nullptr);
    }

    // edge indices flat values 4
    int nEdgeCombinations = nAtoms*(nAtoms-1); // Cheap Combinations nCr(nAtoms 2)
    int ndims_4 = 2;
    int64_t dims_4[] = {nEdgeCombinations, 2};
    int64_t data_4[nEdgeCombinations][2] = {};
    int edge_index = 0;
    for (int i=0; i<nAtoms; i++)
    {   
        for (int j=0; j<nAtoms; j++)
        {
            if (i==j) {
                continue;
            }
            data_4[edge_index][0] = i;
            data_4[edge_index][1] = j;
            edge_index++;
        }
    }
    
    int ndata_4 = sizeof(data_4);
    for (int model_idx=0; model_idx<n_active_models; model_idx++) {
        qm->models[model_idx]->InputValues[4] = TF_NewTensor(TF_INT64, dims_4, ndims_4, data_4, ndata_4, &NoOpDeallocator, nullptr);
    }

    // edge indces row splits 5
    int ndims_5 = 1;
    int64_t dims_5[] = {2};
    int64_t data_5[2] = {0, nEdgeCombinations};
    int ndata_5 = sizeof(data_5);
    for (int model_idx=0; model_idx<n_active_models; model_idx++) {
        qm->models[model_idx]->InputValues[5] = TF_NewTensor(TF_INT64, dims_5, ndims_5, data_5, ndata_5, &NoOpDeallocator, nullptr);
    }

    // angle indces flat values 6
    int nAngleCombinations;
    nAngleCombinations = nAtoms*(nAtoms-1)*(nAtoms-2); // Cheap Combinations nAtoms*nCr(nAtoms-1 2)

    int ndims_6 = 2;
    int64_t dims_6[] = {nAngleCombinations, 3};
    int64_t data_6[nAngleCombinations][3] = {};

    int angle_index = 0;
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
                data_6[angle_index][0] = i;
                data_6[angle_index][1] = j;
                data_6[angle_index][2] = k; 
                angle_index++;
            }
        }
    }
    
    int ndata_6 = sizeof(data_6);
    for (int model_idx=0; model_idx<n_active_models; model_idx++) {
        qm->models[model_idx]->InputValues[6] = TF_NewTensor(TF_INT64, dims_6, ndims_6, data_6, ndata_6, &NoOpDeallocator, nullptr);
    }

    // angle indces row splits 7
    int ndims_7 = 1;
    int64_t dims_7[] = {2};
    int64_t data_7[2] = {0, nAngleCombinations};
    int ndata_7 = sizeof(data_7); 
    for (int model_idx=0; model_idx<n_active_models; model_idx++) {
        qm->models[model_idx]->InputValues[7] = TF_NewTensor(TF_INT64, dims_7, ndims_7, data_7, ndata_7, &NoOpDeallocator, nullptr);
    }

    // total charge 8
    int ndims_8 = 1;
    int64_t dims_8[] = {1};
    float data_8[1] = {(float)qm->QMcharge_get()};
    int ndata_8 = sizeof(data_8);
    for (int model_idx=0; model_idx<n_active_models; model_idx++) {
        qm->models[model_idx]->InputValues[8] = TF_NewTensor(TF_FLOAT, dims_8, ndims_8, data_8, ndata_8, &NoOpDeallocator, nullptr);
    }

    // esp flat values 9
    float V_to_au = 1/27.211386245988;
    int ndims_9 = 1;
    int64_t dims_9[] = {nAtoms};
    float data_9[nAtoms] = {};
    for (int i=0; i<nAtoms; i++)
    {
        data_9[i] = qm->pot_qmmm_get(i)*V_to_au; // in atomic units
    }
    int ndata_9 = sizeof(data_9); 
    for (int model_idx=0; model_idx<n_active_models; model_idx++) {
        qm->models[model_idx]->InputValues[9] = TF_NewTensor(TF_FLOAT, dims_9, ndims_9, data_9, ndata_9, &NoOpDeallocator, nullptr);
    }

    // esp row splits 10
    int ndims_10 = 1;
    int64_t dims_10[] = {2};
    int64_t data_10[2] = {0, nAtoms};
    int ndata_10 = sizeof(data_1);
    for (int model_idx=0; model_idx<n_active_models; model_idx++) {
        qm->models[model_idx]->InputValues[10] = TF_NewTensor(TF_INT64, dims_10, ndims_10, data_10, ndata_10, &NoOpDeallocator, nullptr);
    }

    // esp_grad flat values 11
    snew(ESPgrad, qm->nrQMatoms_get());
    if (qm->qmmm_variant_get() == eqmmmPME)
    {
        snew(ESPgrad_full, qm->nrQMatoms_get());
    }
    int ndims_11 = 2;
    int64_t dims_11[] = {nAtoms, 3};
    float data_11[nAtoms][3] = {};
    // qr->gradient_ESP(//cr, nrnb, wcycle, nullptr, qm->qmmm_variant_get(), data_11, ESPgrad_full);
    qr->gradient_ESP(qm->qmmm_variant_get(), data_11, ESPgrad_full);
    int ndata_11 = sizeof(data_11);
    for (int model_idx=0; model_idx<n_active_models; model_idx++) {
        qm->models[model_idx]->InputValues[11] = TF_NewTensor(TF_FLOAT, dims_11, ndims_11, data_11, ndata_11, &NoOpDeallocator, nullptr);
    }

    // esp grad row splits 12
    int ndims_12 = 1;
    int64_t dims_12[] = {2};
    int64_t data_12[2] = {0, nAtoms};
    int ndata_12 = sizeof(data_12);
    for (int model_idx=0; model_idx<n_active_models; model_idx++) {
        qm->models[model_idx]->InputValues[12] = TF_NewTensor(TF_INT64, dims_12, ndims_12, data_12, ndata_12, &NoOpDeallocator, nullptr);
    }

    // Run the Session
    float* charge_predictions[n_active_models];
    float* energy_predictions[n_active_models];
    float* grad_predictions[n_active_models];
    for (int model_idx=0; model_idx<n_active_models; model_idx++) {
        TF_SessionRun(qm->models[model_idx]->Session, NULL,
            qm->models[model_idx]->Input, qm->models[model_idx]->InputValues, qm->models[model_idx]->NumInputs,
            qm->models[model_idx]->Output, qm->models[model_idx]->OutputValues, qm->models[model_idx]->NumOutputs,
            NULL, 0, NULL, qm->models[model_idx]->Status);

        if(!(TF_GetCode(qm->models[model_idx]->Status) == TF_OK)) {
        printf("%s",TF_Message(qm->models[model_idx]->Status));
        }

        charge_predictions[model_idx] = (float*)TF_TensorData(qm->models[model_idx]->OutputValues[0]); // in e-
        energy_predictions[model_idx] = (float*)TF_TensorData(qm->models[model_idx]->OutputValues[1]); // in Hartree
        grad_predictions[model_idx] = (float*)TF_TensorData(qm->models[model_idx]->OutputValues[2]); // in Hartree/Bohr
    }

    // mean and std of charges
    float charge_means[nAtoms] = {};
    float charge_stds[nAtoms] = {};
    for (int at_idx=0; at_idx<nAtoms; at_idx++) {
        //calculate and save mean of each model prediction
        float mean = 0.0;
        for (int model_idx=0; model_idx<n_active_models; model_idx++) {
            mean += charge_predictions[model_idx][at_idx];
        }
        mean /= n_active_models;
        charge_means[at_idx] = mean;

        //calculate and save std of each model prediction
        float std = 0.0;
        for (int model_idx=0; model_idx<n_active_models; model_idx++) {
            std += std::pow(charge_predictions[model_idx][at_idx] - mean, 2);
        }
        std /= n_active_models;
        std = std::sqrt(std);
        charge_stds[at_idx] = std;
    }

    // mean and std of energy
    float energy_mean = 0.0;
    for (int model_idx=0; model_idx<n_active_models; model_idx++) {
        energy_mean += energy_predictions[model_idx][0];
    }
    energy_mean /= n_active_models;

    float energy_std = 0.0;
    for (int model_idx=0; model_idx<n_active_models; model_idx++) {
        energy_std += std::pow(energy_predictions[model_idx][0] - energy_mean, 2);
    }
    energy_std /= n_active_models;
    energy_std = std::sqrt(energy_std);
    if ((energy_std > energy_prediction_std_threshold) and (energy_prediction_std_threshold > 0.0))
    {
        qm->significant_structure = true;
    }

    // mean and std of gradients
    float grad_means[nAtoms][3] = {};
    float grad_stds[nAtoms][3] = {};
    for (int at_idx=0; at_idx<nAtoms; at_idx++) {
        for (int dim_idx=0; dim_idx<3; dim_idx++) {
            //calculate and save mean of each model prediction
            float mean = 0.0;
            for (int model_idx=0; model_idx<n_active_models; model_idx++) {
                mean += grad_predictions[model_idx][3*at_idx+dim_idx];
            }
            mean /= n_active_models;
            grad_means[at_idx][dim_idx] = mean;
            
            //calculate and save std of each model prediction
            float std = 0.0;
            for (int model_idx=0; model_idx<n_active_models; model_idx++) {
                std += std::pow(grad_predictions[model_idx][3*at_idx+dim_idx] - mean, 2);
            }
            std /= n_active_models;
            std = std::sqrt(std);

            if ((std > force_prediction_std_threshold) and (force_prediction_std_threshold > 0.0))
            {
                qm->significant_structure = true;
            }

            grad_stds[at_idx][dim_idx] = std;
        }
    }

    wallcycle_stop(wcycle, ewcQM);

    /* Save the QM charges */
    for (int i=0; i<n; i++)
    {
        qm->QMcharges_set(i, (real) charge_predictions[0][i]); // sign OK
        //qm->QMcharges_set(i, charge_means[i]); // sign OK
    }
    
    QMener = energy_predictions[0][0];
    //QMener = energy_mean;

    /* Save the gradient on the QM atoms */
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<3; j++)
        {
            QMgrad[i][j] = (real) grad_predictions[0][3*i+j]; // negative of force
            //QMgrad[i][j] = grad_means[i][j]; // negative of force 
        }
    }

    /* Calculate the MM-part of the QM/MM forces
     *   (these are not covered by the model,
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
            f[i][j]      = HARTREE_BOHR2MD*QMgrad[i][j];
         // fshift[i][j] = HARTREE_BOHR2MD*QMgrad[i][j];
        }
    }
    for (int i = 0; i < mm.nrMMatoms; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            f[i+qm->nrQMatoms_get()][j]      = -HARTREE_BOHR2MD*MMgrad[i][j];
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

    char periodic_system[37][3]={"XX",
        "h",                               "he",
        "li","be","b", "c", "n", "o", "f", "ne",
        "na","mg","al","si","p", "s", "cl","ar",
        "k", "ca","sc","ti","v", "cr","mn","fe","co",
        "ni","cu","zn","ga","ge","as","se","br","kr"};

    if (f_x_qm && step % output_freq_x_qm == 0)
    {
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

    if (qm->significant_structure)
    {   
        // // Print information and save information for debugging
        // FILE *f_input_0 = nullptr;
        // f_input_0 = fopen("input_00.txt", "w");
        // for (int i=0; i<nAtoms; i++)
        // {
        //     printf("CHECK Atomic numbers %d %ld\n", i, data_0[i]);
        //     fprintf(f_input_0, "%ld\n", data_0[i]);
        // }
        // fclose(f_input_0);

        // FILE *f_input_2 = nullptr;
        // f_input_2 = fopen("input_02.txt", "w");
        // for (int i=0; i<nAtoms; i++)
        // {
        //     printf("CHECK COORD QM Bohr[%d] = %6.3f %6.3f %6.3f\n", i+1, data_2[i][0], data_2[i][1], data_2[i][2]);
        //     fprintf(f_input_2, "%6.6f %6.6f %6.6f\n", data_2[i][0], data_2[i][1], data_2[i][2]);
        // }
        // for (int i=0; i<nAtoms; i++)
        // {
        //     printf("CHECK COORD QM MD units[%d] = %6.3f %6.3f %6.3f\n", i+1, BOHR2NM*data_2[i][0], BOHR2NM*data_2[i][1], BOHR2NM*data_2[i][2]);
        // }
        // fclose(f_input_2);

        // FILE *f_input_4 = nullptr;
        // f_input_4 = fopen("input_04.txt", "w");
        // for (int i=0; i<nEdgeCombinations; i++)
        // {
        //     fprintf(f_input_4, "%ld %ld\n", data_4[i][0], data_4[i][1]);
        // }
        // printf("CHECK Bonds %d %d\n", edge_index, nEdgeCombinations);
        // fclose(f_input_4);

        // FILE *f_input_6 = nullptr;
        // f_input_6 = fopen("input_06.txt", "w");
        // for (int i=0; i<nAngleCombinations; i++)
        // {
        //     fprintf(f_input_6, "%ld %ld %ld\n", data_6[i][0], data_6[i][1],  data_6[i][2]);
        // }
        // printf("CHECK Angles %d %d\n", angle_index, nAngleCombinations);
        // fclose(f_input_6);

        // FILE *f_input_8 = nullptr;
        // f_input_8 = fopen("input_08.txt", "w");
        // fprintf(f_input_8, "%1.0f\n", data_8[0]);
        // printf("CHECK Total Charge %1.0f \n", data_8[0]);
        // fclose(f_input_8);

        // FILE *f_input_9 = nullptr;
        // f_input_9 = fopen("input_09.txt", "w");
        // for (int i=0; i<nAtoms; i++)
        // {
        //     printf("CHECK ESP QM[%d] = %6.3f\n", i+1, data_9[i]);
        //     fprintf(f_input_9, "%6.6f\n", data_9[i]);
        // }
        // fclose(f_input_9);

        // FILE *f_input_11 = nullptr;
        // f_input_11 = fopen("input_11.txt", "w");
        // for (int i=0; i<nAtoms; i++)
        // {
        //     printf("CHECK ESP GRAD[%d] = %6.3f %6.3f %6.3f\n", i+1, data_11[i][0], data_11[i][1], data_11[i][2]);
        //     fprintf(f_input_11, "%6.6f %6.6f %6.6f\n", data_11[i][0], data_11[i][1], data_11[i][2]);
        // }
        // fclose(f_input_11);

        // FILE *f_output = nullptr;
        // f_output = fopen("output.txt", "w");
        // for (int i=0; i<nAtoms; i++)
        // {
        //     printf("CHECK CHARGE QM[%d] = %6.3f\n", i+1, qm->QMcharges_get(i));
        //     fprintf(f_output, "%6.6f\n", qm->QMcharges_get(i));
        // }

        // printf("CHECK QM Energy = %6.3f\n", QMener);
        // fprintf(f_output, "%6.6f\n", QMener);

        // for (int i=0; i<n; i++)
        // {
        //     printf("CHECK QM Grad[%d] = %6.3f %6.3f %6.3f\n", i+1, QMgrad[i][0], QMgrad[i][1], QMgrad[i][2]);
        //     fprintf(f_output, "%6.6f %6.6f %6.6f\n", QMgrad[i][0], QMgrad[i][1], QMgrad[i][2]);
        // }

        // fclose(f_output);

        FILE *f_std = nullptr; // file for saving standard deviations of significant structures
        f_std = fopen("qm_mlmm_std.xyz", "a");

        fprintf(f_std, "\nMeans of energy predictions step %d\n", step);
        fprintf(f_std, "%8.4f\n", energy_mean);

        fprintf(f_std, "\nStds of energy predictions step %d\n", step);
        fprintf(f_std, "%8.4f\n", energy_std);

        fprintf(f_std, "\nMeans of force predictions step %d\n", step);
        for (int i=0; i<n; i++) {
            fprintf(f_std, "%-2s %4i %8.4f %8.4f %8.4f\n",
                periodic_system[qm->atomicnumberQM_get(i)], i,
                grad_means[i][0], grad_means[i][1], grad_means[i][2]);
        }
        
        fprintf(f_std, "\nStds of force predictions step %d\n", step);
        for (int i=0; i<n; i++) {
            fprintf(f_std, "%-2s %4i %8.4f %8.4f %8.4f\n",
                periodic_system[qm->atomicnumberQM_get(i)], i,
                grad_stds[i][0], grad_stds[i][1], grad_stds[i][2]);
        }
        fclose(f_std);
    }
    
    sfree(QMgrad);
    sfree(MMgrad);
    if (qm->qmmm_variant_get() == eqmmmPME)
    {
        sfree(MMgrad_full);
    }
    sfree(ESPgrad);
    if (qm->qmmm_variant_get() == eqmmmPME)
    {
        sfree(ESPgrad_full);
    }

    sfree(grad);
    sfree(pot);
    sfree(potgrad);
    sfree(q);
    sfree(pot_sr);
    sfree(pot_lr);

    for (int model_idx=0; model_idx<n_active_models; model_idx++) {
        for (int i=0; i<qm->models[model_idx]->NumInputs; i++) {
            TF_DeleteTensor(qm->models[model_idx]->InputValues[i]);
        }
        for (int i=0; i<qm->models[model_idx]->NumOutputs; i++) {
            TF_DeleteTensor(qm->models[model_idx]->OutputValues[i]);
        }
    }

    step++;

    return (real) QMener * HARTREE2KJ * AVOGADRO;
} /* call_nn */

/* end of nn sub routines */

#endif

#pragma GCC diagnostic pop

