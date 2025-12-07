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

#if GMX_QMMM_TENSORFLOW or GMX_QMMM_DFTBPLUS_TENSORFLOW
#include <json/value.h>
#include <json/reader.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <map>

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

#include "gromacs/mdlib/qm_tensorflow.h"
extern "C" {
    #include "tensorflow/c/c_api.h"
}

/* Tensorflow interface routines */

void init_tensorflow(QMMM_QMrec* qm)
{
    // Find amount of models, max 10
    char* n_models_env; // Pointer to env variable for comparison with nullpointer
    int n_models; 
    if ((n_models_env = getenv("GMX_N_MODELS")) == nullptr)
    {
        printf("INFO: Amount of models to use for adaptive sampling can be given via the env variable GMX_N_MODELS\n");
        n_models = 1;
        printf("Using %i model \n", n_models);
    }
    else {
        n_models = std::stoi(n_models_env);
        if (n_models > 10)
        {
            printf("A maximum of 10 models is supported\n");
            exit(-1);
        }
        printf("Using %i models \n", n_models);
    }
    qm->n_models = n_models;
    
    // Find and set network architecture, default "hdnnp4th"
    std::vector<std::string> implemented_network_architectures = {"hdnnp2nd", "hdnnp4th", "schnet", "painn"};
    std::string network_architecture;
    char* env;
    if ((env = getenv("GMX_NN_ARCHITECTURE")) != nullptr)
    {
        network_architecture = std::string(env);
        if (std::find(implemented_network_architectures.begin(), implemented_network_architectures.end(), network_architecture) == implemented_network_architectures.end())
        {
            printf("The network architecture %s is not implemented, please choose from ", network_architecture.c_str());
            for (long unsigned int i=0; i<implemented_network_architectures.size(); i++)
            {
                printf("%s ", implemented_network_architectures[i].c_str());
            }
            printf("\n");
            exit(-1);
        }
        printf("Using %s as network architecture\n", network_architecture.c_str());
    }
    else
    {
        network_architecture = "hdnnp4th";
        printf("Using %s as default network architecture\n", network_architecture.c_str());
    }

    // Find model prefixes and pathes
    char* model_path_prefix;
    char model_path[1000];
    if ((model_path_prefix = getenv("GMX_MODEL_PATH_PREFIX")) == nullptr)
    {
        printf("Path prefix to models must be given in the env variable GMX_MODEL_PATH_PREFIX\n");
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

    // Architecture differences
    int NumInputs;
    int NumOutputs;
    if (network_architecture == "hdnnp2nd")
    {
        NumInputs = 8;
        NumOutputs = 2; // energy, force
    }
    else if (network_architecture == "hdnnp4th")
    {
        NumInputs = 13;
        NumOutputs = 3; // charge, energy, force
    }
    else if ((network_architecture == "schnet") || (network_architecture == "painn"))
    {
        NumInputs = 6;
        NumOutputs = 2; // energy, force
    }

    //********* Read models
    for (int model_idx=0; model_idx<n_models; model_idx++) {
        snew(qm->models[model_idx], 1); // Allocate memory for model
        printf("Reading model %d\n", model_idx);
        std::strcpy(qm->models[model_idx]->modelArchitecture, network_architecture.c_str()); // Save network architecture for later use

        snprintf(model_path, sizeof(model_path), "%s%d", model_path_prefix, model_idx);

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
        qm->models[model_idx]->NumInputs = NumInputs;
        qm->models[model_idx]->Input = (TF_Output*)malloc(sizeof(TF_Output) * qm->models[model_idx]->NumInputs);
        
        const char* base_name = "serving_default_args_0";
        char operation_name[100];

        for (int i=0; i<qm->models[model_idx]->NumInputs; i++) {
            if (i == 0) {
                strcpy(operation_name, base_name); // First input is named "serving_default_args_0" with ragged tensors by default
            }
            else {
                snprintf(operation_name, sizeof(operation_name), "%s_%d", base_name, i); // Subsequent inputs are named "serving_default_args_0_1", "serving_default_args_0_2", ...
            }
            TF_Output t = {TF_GraphOperationByName(qm->models[model_idx]->Graph, operation_name), 0};

            if(t.oper == NULL) {
                printf("ERROR: Failed TF_GraphOperationByName %s\n", operation_name);
                exit(-1);
            }
            else {
                // printf("TF_GraphOperationByName %s is OK\n", operation_name);
                qm->models[model_idx]->Input[i] = t;
            }
        }
        
        //********* Get Output tensor
        qm->models[model_idx]->NumOutputs = NumOutputs;
        qm->models[model_idx]->Output = (TF_Output*)malloc(sizeof(TF_Output) * qm->models[model_idx]->NumOutputs);

        TF_Output output;
        for (int output_index=0; output_index<qm->models[model_idx]->NumOutputs; output_index++) {
            output = {TF_GraphOperationByName(qm->models[model_idx]->Graph, "StatefulPartitionedCall"), output_index};
            if(output.oper == NULL) {
                printf("ERROR: Failed TF_GraphOperationByName StatefulPartitionedCall %d\n", output_index);
                exit(-1);
            }
            else {
                // printf("TF_GraphOperationByName StatefulPartitionedCall is OK\n");
                qm->models[model_idx]->Output[output_index] = output;
            }
        }

        // Allocate data for inputs & outputs
        qm->models[model_idx]->InputValues  = (TF_Tensor**)malloc(sizeof(TF_Tensor*)*qm->models[model_idx]->NumInputs);
        qm->models[model_idx]->OutputValues = (TF_Tensor**)malloc(sizeof(TF_Tensor*)*qm->models[model_idx]->NumOutputs);
        
        TF_DeleteSessionOptions(SessionOpts);
        TF_DeleteBuffer(RunOpts);
    }

    // Load scaler
    snew (qm->scaler, 1);
    if (((env = getenv("GMX_NN_SCALER")) != nullptr) && (env[0] != '\0'))
    {
        printf("Using scaler %s\n", env);
        
        load_scaler(qm, env);
    }
    else
    {
        printf("No scaler used\n");
        qm->scaler->use_scaler = false;
    }

    // Read r_max from environment variable
    if ((env = getenv("GMX_NN_R_MAX")) != nullptr)
    {
        qm->r_max = atof(env);
        qm->r_max_squared = qm->r_max*qm->r_max;
        printf("Using cutoff distance GMX_NN_R_MAX = %.2f Bohr\n", qm->r_max);
    }
    else
    {
        qm->r_max = 20.0; // Default cutoff distance (in Bohr). Adjust to cutoff
        qm->r_max_squared = qm->r_max*qm->r_max;
        printf("Using default cutoff distance GMX_NN_R_MAX = %.2f Bohr\n", qm->r_max);
    }

    fflush(stdout);
    fflush(stderr);
    return;
} /* init_tensorflow */

void NoOpDeallocator(void* data, size_t a, void* b) {(void)b; if (a > 0) {(void)data;}} //Nonsense Deallocator to make compiler happy
void SFreeDeallocator(void* data, size_t a, void* b) {
    bool is_last_model = (bool)(intptr_t)b;  // Convert pointer value to bool
    if ((a > 0) && data != nullptr && is_last_model) { // Free only for last model
        sfree(data);
    }
} //sfree memory after TF_DeleteTensor
real call_tensorflow(QMMM_rec*         qr,
                const t_commrec*  cr,
                QMMM_QMrec*       qm,
                const QMMM_MMrec& mm,
                rvec              f[],
                rvec              fshift[],
                t_nrnb*           nrnb,
                gmx_wallcycle_t   wcycle)
{
    static int step = 0;
    static int output_freq_extxyz;
    static FILE *f_std = nullptr; // file for saving standard deviations

    static float energy_prediction_std_threshold; // energy threshold deviation
    static float force_prediction_std_threshold; // force threshold deviation
    static int nn_eval_freq; // evaluation frequency of additional models

    float QMenergy;
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

        if (((env = getenv("GMX_NN_EXTXYZ")) != nullptr) || ((env = getenv("GMX_PYTORCH_EXTXYZ")) != nullptr))
        {
            output_freq_extxyz = atoi(env);
            printf("The model inputs/outputs will be saved in file qm_mlmm.extxyz every %d steps.\n", output_freq_extxyz);
        }
        else
        {
            output_freq_extxyz = -1; // default no output
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
        // Open file for standard deviations
        if (qm->n_models > 1)
        {
            f_std = fopen("qm_mlmm_std.xyz", "a");
            printf("The standard deviations of the model predictions will be saved in file qm_mlmm_std.xyz\n");
        }

        fflush(stdout);
        fflush(stderr);
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
        qm->pot_qmmm_set(j, (double) - pot[j] * AU2EV); // in volt units
    }

    /* NN call itself */
    wallcycle_start(wcycle, ewcQM);

    int n_models = qm->n_models;
    int nAtoms = qm->nrQMatoms_get();
    int nMMAtoms = mm.nrMMatoms;

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
    
    // Prepare inputs
    wallcycle_start(wcycle, ewcQM_DATA_PREP);
    float energy_conversion;
    float force_conversion;
    float mm_gradient_conversion;
    if (strcmp(qm->models[0]->modelArchitecture, "hdnnp2nd") == 0)
    {
        prepare_hdnnp2nd_inputs(qr, qm, n_active_models);
        energy_conversion = HARTREE2KJMOL; // Hartree to kJ/mol
        force_conversion = HARTREE_BOHR2MD; // Hartree/Bohr to kJ/(mol*nm)
        mm_gradient_conversion = HARTREE_BOHR2MD; // Hartree/Bohr to kJ/(mol*nm)
    }
    else if (strcmp(qm->models[0]->modelArchitecture, "hdnnp4th") == 0)
    {
        prepare_hdnnp4th_inputs(qr, qm, n_active_models);
        energy_conversion = HARTREE2KJMOL; // Hartree to kJ/mol
        force_conversion = HARTREE_BOHR2MD; // Hartree/Bohr to kJ/(mol*nm)
        mm_gradient_conversion = HARTREE_BOHR2MD; // Hartree/Bohr to kJ/(mol*nm)
    }
    else if ((strcmp(qm->models[0]->modelArchitecture, "schnet") == 0) || (strcmp(qm->models[0]->modelArchitecture, "painn") == 0))
    {
        prepare_schnet_painn_inputs(qr, qm, n_active_models);
        energy_conversion = HARTREE2KJMOL; // Hartree to kJ/mol
        force_conversion = HARTREE_BOHR2MD; // Hartree/Bohr to kJ/(mol*nm)
        mm_gradient_conversion = HARTREE_BOHR2MD; // Hartree/Bohr to kJ/(mol*nm)
    }
    else
    {
        printf("The network architecture %s is not implemented\n", qm->models[0]->modelArchitecture);
        exit(-1);
    }
    
    wallcycle_stop(wcycle, ewcQM_DATA_PREP);

    // Run the Session
    wallcycle_start(wcycle, ewcQM_MODEL_CALL);
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

        if (strcmp(qm->models[model_idx]->modelArchitecture, "hdnnp4th") == 0)
        {
            charge_predictions[model_idx] = (float*)TF_TensorData(qm->models[model_idx]->OutputValues[0]); // in e-
            energy_predictions[model_idx] = (float*)TF_TensorData(qm->models[model_idx]->OutputValues[1]); // in Hartree
            grad_predictions[model_idx] = (float*)TF_TensorData(qm->models[model_idx]->OutputValues[2]); // in Hartree/Bohr, gradients, not forces 
        }
        else if (
        (strcmp(qm->models[model_idx]->modelArchitecture, "hdnnp2nd") == 0) ||
        (strcmp(qm->models[model_idx]->modelArchitecture, "schnet") == 0) || 
        (strcmp(qm->models[model_idx]->modelArchitecture, "painn") == 0)
        )
        {
            float charges[nAtoms] = {};
            for (int at_idx=0; at_idx<nAtoms; at_idx++) {
                charges[at_idx] = 0; // in e-
            }
            charge_predictions[model_idx] = (float*)charges;
            energy_predictions[model_idx] = (float*)TF_TensorData(qm->models[model_idx]->OutputValues[0]); // in Hartree
            grad_predictions[model_idx] = (float*)TF_TensorData(qm->models[model_idx]->OutputValues[1]); // in Hartree/Bohr
        }
        else
        {
            printf("The network architecture %s is not implemented\n", qm->models[model_idx]->modelArchitecture);
            exit(-1);
        }

        // Unscale predictions
        if (qm->scaler->use_scaler)
        {
            inverse_transform(qm, energy_predictions[model_idx], grad_predictions[model_idx]);
        }
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
        //qm->QMcharges_set(i, (real) charge_predictions[0][i]); // Using the first model
        qm->QMcharges_set(i, charge_means[i]); // Using the mean of all models
    }
    
    //QMenergy = energy_predictions[0][0]; // Using the first model
    QMenergy = energy_mean; // Using the mean of all models

    /* Save the gradient on the QM atoms */
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<3; j++)
        {
            //QMgrad[i][j] = (real) grad_predictions[0][3*i+j]; // negative of force // Using the first model
            QMgrad[i][j] = grad_means[i][j]; // negative of force  // Using the mean of all models
        }
    }

    /* Calculate the MM-part of the QM/MM forces
     *   (these are not covered by the model,
     *    because it does not know the position and charges of individual atoms,
     *    and instead, it only obtains the external potentials induced at QM atoms.
     */
    snew(MMgrad, nMMAtoms);
    if (qm->qmmm_variant_get() == eqmmmPME)
    {
        snew(MMgrad_full, mm.nrMMatoms_full);
    }

    for (int j=0; j<nMMAtoms; j++)
    {
        clear_rvec(MMgrad[j]);
    }
    if (qm->qmmm_variant_get() == eqmmmPME)
    {
        for (int j=0; j<mm.nrMMatoms_full; j++)
        {
        clear_rvec(MMgrad_full[j]);
        }
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
            f[i][j]      = force_conversion*QMgrad[i][j];
         // fshift[i][j] = force_conversion*QMgrad[i][j];
        }
    }
    for (int i = 0; i < nMMAtoms; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            f[i+qm->nrQMatoms_get()][j]      = mm_gradient_conversion*MMgrad[i][j];
         // fshift[i+qm->nrQMatoms_get()][j] = mm_gradient_conversion*MMgrad[i][j];
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
                f[i+qm->nrQMatoms_get()+nMMAtoms][j]      = mm_gradient_conversion*MMgrad_full[i][j];
             // fshift[i+qm->nrQMatoms_get()+nMMAtoms][j] = mm_gradient_conversion*MMgrad_full[i][j];
            }
        }
    }

    // Print inputs and outputs of the model every output_freq_extxyz steps
    if ((output_freq_extxyz > 0) && (step % output_freq_extxyz == 0))
    {
        write_inputs_outputs(qm, step);
    }

    
    char periodic_system[37][3]={"XX",
        "H",                               "He",
        "Li","Be","B", "C", "N", "O", "F", "Ne",
        "Na","Mg","Al","Si","P", "S", "Cl","Ar",
        "K", "Ca","Sc","Ti","V", "Cr","Mn","Fe","Co",
        "Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr"
    };


    if (qm->significant_structure)
    {   
        write_inputs_outputs(qm, step);
    }

    // Save the means and stds of the predictions to a file
    if (n_active_models > 1) {    
        // Save the means
        fprintf(f_std, "%d\n", n);
        fprintf(f_std, "Means step %d: %8.4f\n", step, energy_mean);

        for (int i=0; i<n; i++) {
            fprintf(f_std, "%-2s %8.4f %8.4f %8.4f\n",
                periodic_system[qm->atomicnumberQM_get(i)],
                grad_means[i][0], grad_means[i][1], grad_means[i][2]);
        }
        
        // Save the stds
        fprintf(f_std, "%d\n", n);   
        fprintf(f_std, "Stds step %d: %8.4f\n", step, energy_std);
        for (int i=0; i<n; i++) {
            fprintf(f_std, "%-2s %8.4f %8.4f %8.4f\n",
                periodic_system[qm->atomicnumberQM_get(i)],
                grad_stds[i][0], grad_stds[i][1], grad_stds[i][2]);
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

    for (int model_idx=0; model_idx<n_active_models; model_idx++) {
        for (int i=0; i<qm->models[model_idx]->NumInputs; i++) {
            TF_DeleteTensor(qm->models[model_idx]->InputValues[i]);
        }
        for (int i=0; i<qm->models[model_idx]->NumOutputs; i++) {
            TF_DeleteTensor(qm->models[model_idx]->OutputValues[i]);
        }
    }

    step++;

    return (real) QMenergy * energy_conversion; // in kJ/mol
} /* call_tensorflow */

void prepare_hdnnp2nd_inputs(QMMM_rec* qr,
                            QMMM_QMrec* qm,
                            int n_active_models)
{
    typedef int64_t int_triple[3];

    // Prepare base inputs (0-5) using prepare_schnet_painn_inputs
    prepare_schnet_painn_inputs(qr, qm, n_active_models);

    int nAtoms = qm->nrQMatoms_get();
    float r_max = qm->r_max;
    float geometry_conversion = 1.0 / BOHR2NM; // nm to bohr

    // angle indces flat values 6
    // First, count the number of angles within the cutoff distance
    int nAngleCombinations = 0;
    for (int i=0; i<nAtoms; i++)
    {
        for (int j=0; j<nAtoms; j++)
        {
            if (i == j) {
                continue;
            }

            float r_ij = calculate_distance(qm, i, j) * geometry_conversion;
            
            if (r_ij > r_max) {
                continue;
            }
            
            for (int k=j+1; k<nAtoms; k++)
            {
                if (i == k)  {
                    continue;
                }
                
                float r_ik = calculate_distance(qm, i, k) * geometry_conversion;
                
                if (r_ik > r_max) {
                    continue;
                }
                
                float r_jk = calculate_distance(qm, j, k) * geometry_conversion;
                
                if (r_jk > r_max) {
                    continue;
                }
                
                nAngleCombinations += 2;
            }
        }
    }
    
    int ndims_6 = 2;
    int64_t dims_6[] = {nAngleCombinations, 3};
    int_triple* data_6;
    snew(data_6, nAngleCombinations);

    int angle_index = 0;
    for (int i=0; i<nAtoms; i++)
    {
        for (int j=0; j<nAtoms; j++)
        {
            if (i == j) {
                continue;
            }
            
            float r_ij = calculate_distance(qm, i, j) * geometry_conversion;
            
            if (r_ij > r_max) {
                continue;
            }
            
            for (int k=j+1; k<nAtoms; k++)
            {
                if (i == k)  {
                    continue;
                }
                
                float r_ik = calculate_distance(qm, i, k) * geometry_conversion;
                
                if (r_ik > r_max) {
                    continue;
                }
                
                float r_jk = calculate_distance(qm, j, k) * geometry_conversion;
                
                if (r_jk > r_max) {
                    continue;
                }
                
                data_6[angle_index][0] = i;
                data_6[angle_index][1] = j;
                data_6[angle_index][2] = k;
                angle_index++;
                
                data_6[angle_index][0] = i;
                data_6[angle_index][1] = k;
                data_6[angle_index][2] = j; 
                angle_index++;
            }
        }
    }
    int ndata_6 = sizeof(*data_6)*nAngleCombinations;

    // angle indces row splits 7
    int ndims_7 = 1;
    int64_t dims_7[] = {2};
    static int64_t data_7[2] = {0, 0};
    data_7[1] = nAngleCombinations;
    int ndata_7 = sizeof(data_7);

    // Assign inputs to all models in one loop
    for (int model_idx=0; model_idx<n_active_models; model_idx++) {
        bool is_last_model = (model_idx == n_active_models-1);
        
        qm->models[model_idx]->InputValues[6] = TF_NewTensor(TF_INT64, dims_6, ndims_6, data_6, ndata_6, &SFreeDeallocator, (void*)(intptr_t)is_last_model);
        qm->models[model_idx]->InputValues[7] = TF_NewTensor(TF_INT64, dims_7, ndims_7, data_7, ndata_7, &NoOpDeallocator, nullptr);
    }


    // Check input availability
    if (std::strcmp(qm->models[0]->modelArchitecture, "hdnnp2nd") == 0) {
        for (int input_idx=0; input_idx<qm->models[0]->NumInputs; input_idx++) {
            if (qm->models[0]->InputValues[input_idx] == NULL)
            {
                printf("ERROR: Input %d is NULL\n", input_idx);
                exit(-1);
            }
        }
    }
} // end of prepare_hdnnp2nd_inputs

void prepare_hdnnp4th_inputs(QMMM_rec* qr,
                            QMMM_QMrec* qm,
                            int n_active_models)
{
    prepare_hdnnp2nd_inputs(qr, qm, n_active_models); // First prepare the inputs of hdnnp2nd

    int nAtoms = qm->nrQMatoms_get();

    // total charge 8
    int ndims_8 = 1;
    int64_t dims_8[] = {1};
    static float data_8[1] = {(float)qm->QMcharge_get()};
    int ndata_8 = sizeof(data_8);

    // esp flat values 9
    int ndims_9 = 1;
    int64_t dims_9[] = {nAtoms};
    float* data_9;
    snew(data_9, nAtoms);
    for (int i=0; i<nAtoms; i++)
    {
        data_9[i] = qm->pot_qmmm_get(i)*V2AU; // from volt to atomic units
    }
    int ndata_9 = sizeof(*data_9)*nAtoms;

    // esp row splits 10
    int ndims_10 = 1;
    int64_t dims_10[] = {2};
    static int64_t data_10[2] = {0, nAtoms};
    int ndata_10 = sizeof(data_10);

    // esp_grad flat values 11
    rvec *ESPgrad = nullptr, *ESPgrad_full = nullptr;
    snew(ESPgrad, nAtoms); // Not really needed, as I put this data directly into data_11
    if (qm->qmmm_variant_get() == eqmmmPME)
    {
        snew(ESPgrad_full, nAtoms);
    }
    int ndims_11 = 2;
    int64_t dims_11[] = {nAtoms, 3};
    rvec* data_11;
    snew(data_11, nAtoms);
    // qr->gradient_ESP(//cr, nrnb, wcycle, nullptr, qm->qmmm_variant_get(), data_11, ESPgrad_full);
    qr->gradient_ESP(qm->qmmm_variant_get(), data_11, ESPgrad_full);
    int ndata_11 = sizeof(*data_11)*nAtoms;

    // esp grad row splits 12
    int ndims_12 = 1;
    int64_t dims_12[] = {2};
    static int64_t data_12[2] = {0, nAtoms};
    int ndata_12 = sizeof(data_12);

    // Assign additional hdnnp4th inputs to all models in one loop
    for (int model_idx=0; model_idx<n_active_models; model_idx++) {
        bool is_last_model = (model_idx == n_active_models-1);
        
        qm->models[model_idx]->InputValues[8] = TF_NewTensor(TF_FLOAT, dims_8, ndims_8, data_8, ndata_8, &NoOpDeallocator, nullptr);
        qm->models[model_idx]->InputValues[9] = TF_NewTensor(TF_FLOAT, dims_9, ndims_9, data_9, ndata_9, &SFreeDeallocator, (void*)(intptr_t)is_last_model);
        qm->models[model_idx]->InputValues[10] = TF_NewTensor(TF_INT64, dims_10, ndims_10, data_10, ndata_10, &NoOpDeallocator, nullptr);
        qm->models[model_idx]->InputValues[11] = TF_NewTensor(TF_FLOAT, dims_11, ndims_11, data_11, ndata_11, &SFreeDeallocator, (void*)(intptr_t)is_last_model);
        qm->models[model_idx]->InputValues[12] = TF_NewTensor(TF_INT64, dims_12, ndims_12, data_12, ndata_12, &NoOpDeallocator, nullptr);
    }

    sfree(ESPgrad);
    if (qm->qmmm_variant_get() == eqmmmPME)
    {
        sfree(ESPgrad_full);
    }

    // Check input availability
    for (int input_idx=0; input_idx<qm->models[0]->NumInputs; input_idx++) {
        if (qm->models[0]->InputValues[input_idx] == NULL)
        {
            printf("ERROR: Input %d is NULL\n", input_idx);
            exit(-1);
        }
    }
} // end of prepare_hdnnp_inputs


void prepare_schnet_painn_inputs(   QMMM_rec* qr,
                                    QMMM_QMrec* qm,
                                    int n_active_models)
{
    typedef int64_t int_pair[2];

    int nAtoms = qm->nrQMatoms_get();
    float r_max = qm->r_max;
    float geometry_conversion = 1.0 / BOHR2NM; // nm to bohr

    // node number flat values 0
    int ndims_0 = 1;
    int64_t dims_0[] = {nAtoms};
    int64_t* data_0;
    snew(data_0, nAtoms);
    for (int i=0; i<nAtoms; i++)
    {
        data_0[i] = qm->atomicnumberQM_get(i);
    }
    int ndata_0 = nAtoms * sizeof(*data_0);

    // node number row splits 1
    int ndims_1 = 1;
    int64_t dims_1[] = {2};
    static int64_t data_1[2] = {0, nAtoms};
    int ndata_1 = sizeof(data_1);

    // node coordinates flat values 2
    int ndims_2 = 2;
    int64_t dims_2[] = {nAtoms, 3};
    rvec* data_2;
    snew(data_2, nAtoms);
    for (int i=0; i<nAtoms; i++)
    {
        for (int j=0; j<3; j++)
        {
            data_2[i][j] = qm->xQM_get(i,j) * geometry_conversion;
        }
    }
    int ndata_2 = nAtoms * sizeof(*data_2);

    // node coordinates row splits 3
    int ndims_3 = 1;
    int64_t dims_3[] = {2};
    static int64_t data_3[2] = {0, nAtoms};
    int ndata_3 = sizeof(data_3);

    // edge indices flat values 4
    // First, count the number of edges within the cutoff distance (only j>i, then create both directions)
    int nEdgeCombinations = 0;
    for (int i=0; i<nAtoms; i++)
    {
        for (int j=i+1; j<nAtoms; j++)
        {
            float r = calculate_distance(qm, i, j) * geometry_conversion;
            
            if (r > r_max) {
                continue;
            }

            nEdgeCombinations += 2; // both directions
        }
    }
    
    int ndims_4 = 2;
    int64_t dims_4[] = {nEdgeCombinations, 2};
    int_pair* data_4;
    snew(data_4, nEdgeCombinations);
    int edge_index = 0;
    for (int i=0; i<nAtoms; i++)
    {   
        for (int j=i+1; j<nAtoms; j++)
        {
            float r = calculate_distance(qm, i, j) * geometry_conversion;
            
            if (r > r_max) {
                continue;
            }
            
            // Add both directions
            data_4[edge_index][0] = i;
            data_4[edge_index][1] = j;
            edge_index++;
            
            data_4[edge_index][0] = j;
            data_4[edge_index][1] = i;
            edge_index++;
        }
    }
    int ndata_4 = sizeof(*data_4)*nEdgeCombinations;

    // edge indces row splits 5
    int ndims_5 = 1;
    int64_t dims_5[] = {2};
    static int64_t data_5[2] = {0, 0};
    data_5[1] = nEdgeCombinations; // Update of static array after initialization
    int ndata_5 = sizeof(data_5);

    // Assign inputs to all models in one loop
    for (int model_idx=0; model_idx<n_active_models; model_idx++) {
        bool is_last_model = (model_idx == n_active_models-1);
        
        qm->models[model_idx]->InputValues[0] = TF_NewTensor(TF_INT64, dims_0, ndims_0, data_0, ndata_0, &SFreeDeallocator, (void*)(intptr_t)is_last_model);
        qm->models[model_idx]->InputValues[1] = TF_NewTensor(TF_INT64, dims_1, ndims_1, data_1, ndata_1, &NoOpDeallocator, nullptr);
        qm->models[model_idx]->InputValues[2] = TF_NewTensor(TF_FLOAT, dims_2, ndims_2, data_2, ndata_2, &SFreeDeallocator, (void*)(intptr_t)is_last_model);
        qm->models[model_idx]->InputValues[3] = TF_NewTensor(TF_INT64, dims_3, ndims_3, data_3, ndata_3, &NoOpDeallocator, nullptr);
        qm->models[model_idx]->InputValues[4] = TF_NewTensor(TF_INT64, dims_4, ndims_4, data_4, ndata_4, &SFreeDeallocator, (void*)(intptr_t)is_last_model);
        qm->models[model_idx]->InputValues[5] = TF_NewTensor(TF_INT64, dims_5, ndims_5, data_5, ndata_5, &NoOpDeallocator, nullptr);
    }

    // Check input availability
    if (strcmp(qm->models[0]->modelArchitecture, "schnet") == 0 ||
        strcmp(qm->models[0]->modelArchitecture, "painn") == 0) {
        for (int input_idx=0; input_idx<qm->models[0]->NumInputs; input_idx++) {
            if (qm->models[0]->InputValues[input_idx] == NULL)
            {
                printf("ERROR: Input %d is NULL\n", input_idx);
                exit(-1);
            }
        }
    }

    (void)qr; // placeholder to avoid warning
} // end of prepare_schnet_painn_inputs

void write_schnet_painn_inputs(QMMM_QMrec* qm                            
)
{
    // Print information and save information for debugging

    int nAtoms = qm->nrQMatoms_get();
    // Get actual number of edges from tensor dimensions
    int nEdgeCombinations = TF_Dim(qm->models[0]->InputValues[4], 0);

    // node number flat values 0
    FILE* f_input_0 = nullptr;
    f_input_0 = fopen("input_00.txt", "w");
    int64_t* data_0 = (int64_t*)(TF_TensorData(qm->models[0]->InputValues[0]));
    for (int i=0; i<nAtoms; i++)
    {
        printf("CHECK Atomic numbers %d %ld\n", i, data_0[i]);
        fprintf(f_input_0, "%ld\n", data_0[i]);
    }
    fclose(f_input_0);

    // node coordinates flat values 2
    FILE* f_input_2 = nullptr;
    f_input_2 = fopen("input_02.txt", "w");
    float* data_2 = (float*)(TF_TensorData(qm->models[0]->InputValues[2]));
    for (int i=0; i<nAtoms; i++)
    {
        printf("CHECK COORD QM Bohr[%d] = %6.3f %6.3f %6.3f\n", i+1, data_2[3*i+0], data_2[3*i+1], data_2[3*i+2]);
        fprintf(f_input_2, "%6.6f %6.6f %6.6f\n", data_2[3*i+0], data_2[3*i+1], data_2[3*i+2]);
    }
    for (int i=0; i<nAtoms; i++)
    {
        printf("CHECK COORD QM MD units[%d] = %6.3f %6.3f %6.3f\n", i+1, BOHR2NM*data_2[3*i+0], BOHR2NM*data_2[3*i+1], BOHR2NM*data_2[3*i+2]);
    }
    fclose(f_input_2);

    // edge indices flat values 4
    FILE* f_input_4 = nullptr;
    f_input_4 = fopen("input_04.txt", "w");
    int64_t* data_4 = (int64_t*)(TF_TensorData(qm->models[0]->InputValues[4]));
    for (int i=0; i<nEdgeCombinations; i++)
    {
        fprintf(f_input_4, "%ld %ld\n", data_4[2*i+0], data_4[2*i+1]);
    }
    fclose(f_input_4);
} // end of write_schnet_painn_inputs

void write_hdnnp2nd_inputs(QMMM_QMrec* qm)
{
    // Print information and save information for debugging

    // Get actual number of angles from tensor dimensions
    int nAngleCombinations = TF_Dim(qm->models[0]->InputValues[6], 0);

    write_schnet_painn_inputs(qm); // write inputs 0,2,4

    // angle indces flat values 6
    FILE* f_input_6 = nullptr;
    f_input_6 = fopen("input_06.txt", "w");
    int64_t* data_6 = (int64_t*)(TF_TensorData(qm->models[0]->InputValues[6]));
    for (int i=0; i<nAngleCombinations; i++)
    {
        fprintf(f_input_6, "%ld %ld %ld\n", data_6[3*i+0], data_6[3*i+1],  data_6[3*i+2]);
    }
    fclose(f_input_6);
} // end of write_hdnnp2nd_inputs

void write_hdnnp4th_inputs(QMMM_QMrec* qm)
{
    write_hdnnp2nd_inputs(qm); // write inputs 0,2,4,6

    int nAtoms = qm->nrQMatoms_get();

    // total charge 8
    FILE* f_input_8 = nullptr;
    f_input_8 = fopen("input_08.txt", "w");
    float* data_8 = (float*)(TF_TensorData(qm->models[0]->InputValues[8]));
    fprintf(f_input_8, "%1.0f\n", data_8[0]);
    printf("CHECK Total Charge %1.0f \n", data_8[0]);
    fclose(f_input_8);

    // esp flat values 9
    FILE* f_input_9 = nullptr;
    f_input_9 = fopen("input_09.txt", "w");
    float* data_9 = (float*)(TF_TensorData(qm->models[0]->InputValues[9]));
    for (int i=0; i<nAtoms; i++)
    {
        printf("CHECK ESP QM[%d] = %6.3f\n", i+1, data_9[i]);
        fprintf(f_input_9, "%6.6f\n", data_9[i]);
    }
    fclose(f_input_9);

    // esp_grad flat values 11
    FILE* f_input_11 = nullptr;
    f_input_11 = fopen("input_11.txt", "w");
    float* data_11 = (float*)(TF_TensorData(qm->models[0]->InputValues[11]));
    for (int i=0; i<nAtoms; i++)
    {
        printf("CHECK ESP GRAD[%d] = %6.3f %6.3f %6.3f\n", i+1, data_11[3*i+0], data_11[3*i+1], data_11[3*i+2]);
        fprintf(f_input_11, "%6.6f %6.6f %6.6f\n", data_11[3*i+0], data_11[3*i+1], data_11[3*i+2]);
    }
    fclose(f_input_11);
}

void write_base_outputs(QMMM_QMrec* qm)
{
    // Print and save only energy and forces
    int nAtoms = qm->nrQMatoms_get();

    float* energy_predictions = (float*)TF_TensorData(qm->models[0]->OutputValues[0]); // in Hartree
    float* grad_predictions = (float*)TF_TensorData(qm->models[0]->OutputValues[1]); // in Hartree/Bohr

    FILE* f_output = nullptr;
    f_output = fopen("output.txt", "w");

    printf("CHECK QM Energy = %6.3f\n", energy_predictions[0]);
    fprintf(f_output, "%6.6f\n", energy_predictions[0]);

    for (int i=0; i<nAtoms; i++)
    {
        printf("CHECK QM Grad[%d] = %6.3f %6.3f %6.3f\n", i+1, grad_predictions[3*i+0], grad_predictions[3*i+1], grad_predictions[3*i+2]);
        fprintf(f_output, "%6.6f %6.6f %6.6f\n", grad_predictions[3*i+0], grad_predictions[3*i+1], grad_predictions[3*i+2]);
    }
    fclose(f_output);
} // end of write_base_outputs

void write_extended_outputs(QMMM_QMrec* qm)
{
    // Print and save energy, forces, charges
    int nAtoms = qm->nrQMatoms_get();

    float* charge_predictions = (float*)TF_TensorData(qm->models[0]->OutputValues[0]); // in e-
    float* energy_predictions = (float*)TF_TensorData(qm->models[0]->OutputValues[1]); // in Hartree
    float* grad_predictions = (float*)TF_TensorData(qm->models[0]->OutputValues[2]); // in Hartree/Bohr

    FILE* f_output = nullptr;
    f_output = fopen("output.txt", "w");
    for (int i=0; i<nAtoms; i++)
    {
        printf("CHECK CHARGE QM[%d] = %6.3f\n", i+1, charge_predictions[i]);
        fprintf(f_output, "%6.6f\n", charge_predictions[i]);
    }

    printf("CHECK QM Energy = %6.3f\n", energy_predictions[0]);
    fprintf(f_output, "%6.6f\n", energy_predictions[0]);

    for (int i=0; i<nAtoms; i++)
    {
        printf("CHECK QM Grad[%d] = %6.3f %6.3f %6.3f\n", i+1, grad_predictions[3*i+0], grad_predictions[3*i+1], grad_predictions[3*i+2]);
        fprintf(f_output, "%6.6f %6.6f %6.6f\n", grad_predictions[3*i+0], grad_predictions[3*i+1], grad_predictions[3*i+2]);
    }
    fclose(f_output);
} // end of write_extended_outputs

void write_extxyz(
    QMMM_QMrec* qm,
    int step
)
{
    // Write all inputs and outputs in the extended xyz format
    float geometry_conversion = BOHR2A; // Bohr to Angstrom
    float charge_conversion = 1.0; // e- to e-
    float energy_conversion = AU2EV; // Hartree to eV
    float force_conversion = AU2EV_A*-1; // Hartree/Bohr to eV/Ang, inverse sign to get physical forces
    float esp_conversion = AU2EV; // Hartree/e to eV/e=V
    float esp_grad_conversion = AU2EV_A; // Hartree/e/Bohr to eV/e/Ang
    char periodic_system[37][3]={"XX",
        "H",                               "He",
        "Li","Be","B", "C", "N", "O", "F", "Ne",
        "Na","Mg","Al","Si","P", "S", "Cl","Ar",
        "K", "Ca","Sc","Ti","V", "Cr","Mn","Fe","Co",
        "Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr"};

    bool is_hdnnp4th = false;
    if (std::strcmp(qm->models[0]->modelArchitecture, "hdnnp4th") == 0)
    {
        is_hdnnp4th = true;
    }

    int nAtoms = qm->nrQMatoms_get();
    float* data_2 = (float*)(TF_TensorData(qm->models[0]->InputValues[2]));
    float* data_8 = nullptr;
    float* data_9 = nullptr;
    float* data_11 = nullptr;
    float* charge_predictions = nullptr;
    float* energy_predictions = nullptr;
    float* grad_predictions = nullptr;

    if (!is_hdnnp4th)
    {
        energy_predictions = (float*)TF_TensorData(qm->models[0]->OutputValues[0]); // in Hartree
        grad_predictions = (float*)TF_TensorData(qm->models[0]->OutputValues[1]); // in Hartree/Bohr        
    }
    else {
        data_8 = (float*)(TF_TensorData(qm->models[0]->InputValues[8]));
        data_9 = (float*)(TF_TensorData(qm->models[0]->InputValues[9]));
        data_11 = (float*)(TF_TensorData(qm->models[0]->InputValues[11]));
        charge_predictions = (float*)TF_TensorData(qm->models[0]->OutputValues[0]); // in e-
        energy_predictions = (float*)TF_TensorData(qm->models[0]->OutputValues[1]); // in Hartree
        grad_predictions = (float*)TF_TensorData(qm->models[0]->OutputValues[2]); // in Hartree/Bohr
    }

    FILE* f_extxyz = nullptr;
    f_extxyz = fopen("qm_mlmm.extxyz", "a"); 
    fprintf(f_extxyz, "%d\n", nAtoms);
    if (!is_hdnnp4th) {
        fprintf(f_extxyz, "Properties=species:S:1:pos:R:3:gromacs_force:R:3 ");
    }
    else {
        fprintf(f_extxyz, "Properties=species:S:1:pos:R:3:gromacs_force:R:3:gromacs_charge:R:1:esp:R:1:esp_gradient:R:3 ");
        fprintf(f_extxyz, "total_charge=%.1f ", data_8[0]*charge_conversion);
    }
    fprintf(f_extxyz, "gromacs_energy=%.6f ", energy_predictions[0]*energy_conversion);
    fprintf(f_extxyz, "comment=\"Step %d\" ", step);
    fprintf(f_extxyz, "pbc=\"F F F\" \n");
    for (int i=0; i<nAtoms; i++)
    {
        fprintf(f_extxyz, "%-2s ", periodic_system[qm->atomicnumberQM_get(i)]);
        fprintf(f_extxyz, "%8.4f %8.4f %8.4f ", 
            data_2[3*i+0]*geometry_conversion, 
            data_2[3*i+1]*geometry_conversion, 
            data_2[3*i+2]*geometry_conversion
        );
        fprintf(f_extxyz, "%8.4f %8.4f %8.4f ",
            grad_predictions[3*i+0]*force_conversion,
            grad_predictions[3*i+1]*force_conversion,
            grad_predictions[3*i+2]*force_conversion
        );
        if (!is_hdnnp4th) {
            fprintf(f_extxyz, "\n");
            continue;
        }
        // hdnnp4th additional properties
        fprintf(f_extxyz, "%8.4f ", charge_predictions[i]*charge_conversion);
        fprintf(f_extxyz, "%8.4f ", data_9[i]*esp_conversion);
        fprintf(f_extxyz, "%8.4f %8.4f %8.4f",
            data_11[3*i+0]*esp_grad_conversion,
            data_11[3*i+1]*esp_grad_conversion,
            data_11[3*i+2]*esp_grad_conversion
        );
        fprintf(f_extxyz, "\n");
    }
    fclose(f_extxyz);

    // Unscale predictions not necessary, the qm struct got changed with the scaler already due to the pointer
}

void write_inputs_outputs(QMMM_QMrec* qm,
                                int step                               
)
{
    // if (strcmp(qm->models[0]->modelArchitecture, "hdnnp2nd") == 0)
    // {
    //     write_hdnnp2nd_inputs(qm);
    //     write_base_outputs(qm);
    // }
    // else if (strcmp(qm->models[0]->modelArchitecture, "hdnnp4th") == 0)
    // {
    //     write_hdnnp4th_inputs(qm);
    //     write_extended_outputs(qm);
    // }
    // else if ((strcmp(qm->models[0]->modelArchitecture, "schnet") == 0) || (strcmp(qm->models[0]->modelArchitecture, "painn") == 0))
    // {
    //     write_schnet_painn_inputs(qm);
    //     write_base_outputs(qm);
    // }
    // else
    // {
    //     printf("ERROR: Unknown model architecture %s\n", qm->models[0]->modelArchitecture);
    //     exit(-1);
    // }
    write_extxyz(qm, step);
} // end of write_inputs_outputs

void load_scaler(QMMM_QMrec* qm, char* scaler_file)
{
    std::ifstream file(scaler_file);
    Json::Reader reader;
    Json::Value root;
    Json::Value weights;

    // Input checks
    FILE* f_scaler = nullptr;
    f_scaler = fopen(scaler_file, "r");
    if (f_scaler == NULL) {
        printf("ERROR: Could not open scaler file %s\n", scaler_file);
        exit(-1);
    }
    fclose(f_scaler);

    if (!reader.parse(file, root)) {
        std::cout << "Failed to parse JSON\n" << reader.getFormattedErrorMessages();
        exit(-1);
    }

    if (!root.isMember("weights")) {
        std::cout << "No weights in scaler JSON\n";
        exit(-1);
    }
    weights = root["weights"];

    if (!weights.isMember("n_features_in_")) {
        std::cout << "No n_features_in_ in scaler JSON\n";
        exit(-1);
    }
    if (!weights.isMember("scale_")) {
        std::cout << "No scale_ in scaler JSON\n";
        exit(-1);
    }
    if (!weights.isMember("intercept_")) {
        std::cout << "No intercept_ in scaler JSON\n";
        exit(-1);
    }
    if (!weights.isMember("coef_")) {
        std::cout << "No coef_ in scaler JSON\n";
        exit(-1);
    }
    if (!weights["coef_"].isArray()) {
        std::cout << "coef_ is not an array in scaler JSON\n";
        exit(-1);
    }
    if (!weights["_fit_atom_selection"]) {
        std::cout << "No _fit_atom_selection in scaler JSON\n";
        exit(-1);
    }
    if (!weights["_fit_atom_selection"].isArray()) {
        std::cout << "_fit_atom_selection is not an array in scaler JSON\n";
        exit(-1);
    }
    // if (!weights["_fit_atom_selection_mask"]) {
    //     std::cout << "No _fit_atom_selection_mask in scaler JSON\n";
    //     exit(-1);
    // }
    // if (!weights["_fit_atom_selection_mask"].isArray()) {
    //     std::cout << "_fit_atom_selection_mask is not an array in scaler JSON\n";
    //     exit(-1);
    // }

    // Load weights
    qm->scaler->use_scaler = true;

    qm->scaler->n_features_in_ = weights["n_features_in_"].asInt();
    qm->scaler->scale_ = weights["scale_"].asFloat();
    qm->scaler->intercept_ = weights["intercept_"].asFloat();
    float* coef_;
    snew (coef_, qm->scaler->n_features_in_);
    for (int i=0; i<qm->scaler->n_features_in_; i++) {
        coef_[i] = weights["coef_"][i].asFloat();
    }
    qm->scaler->coef_ = coef_;

    int* _fit_atom_selection;
    snew (_fit_atom_selection, qm->scaler->n_features_in_);
    for (int i=0; i<qm->scaler->n_features_in_; i++) {
        _fit_atom_selection[i] = weights["_fit_atom_selection"][i].asInt();
    }
    qm->scaler->_fit_atom_selection = _fit_atom_selection;

    // Didn't use this mask, redundant with _fit_atom_selection
    // bool* _fit_atom_selection_mask;
    // snew (_fit_atom_selection_mask, weights["_fit_atom_selection_mask"].size());
    // for (unsigned int i=0; i<weights["_fit_atom_selection_mask"].size(); i++) {
    //     _fit_atom_selection_mask[i] = weights["_fit_atom_selection_mask"][i].asBool();
    // }
    // qm->scaler->_fit_atom_selection_mask = _fit_atom_selection_mask;

    // // Print loaded weights
    // std::cout << "n_features_in: " << qm->scaler->n_features_in_ << "\n";
    // std::cout << "scale: " << qm->scaler->scale_ << "\n";
    // std::cout << "intercept: " << qm->scaler->intercept_ << "\n";
    // std::cout << "coef_: ";
    // for (int i=0; i<qm->scaler->n_features_in_; i++) {
    //     std::cout << qm->scaler->coef_[i] << " ";
    // }
    // std::cout << "\n";
    // std::cout << "_fit_atom_selection: ";
    // for (int i=0; i<qm->scaler->n_features_in_; i++) {
    //     std::cout << qm->scaler->_fit_atom_selection[i] << " ";
    // }
    // std::cout << "\n";
    // std::cout << "_fit_atom_selection_mask: ";
    // for (unsigned int i=0; i<weights["_fit_atom_selection_mask"].size(); i++) {
    //     std::cout << qm->scaler->_fit_atom_selection_mask[i] << " ";
    // }
    // std::cout << "\n";
} // end of load_scaler

void inverse_transform(QMMM_QMrec* qm, float* energy_predictions, float* grad_predictions)
// Could maybe be done without passing energy_predictions and grad_predictions and grad_predictions
// Because they point to the qm struct property
{
    int n_atoms = qm->nrQMatoms_get();
    int atomic_numbers[n_atoms];
    for (int i=0; i<n_atoms; i++)
    {
        atomic_numbers[i] = qm->atomicnumberQM_get(i);
    }

    std::set<int> unique_elements(atomic_numbers, atomic_numbers + n_atoms); // unique atomic numbers as set of atomic numbers
    int n_unique_elements = (int)unique_elements.size();
    // printf("Unique elements: ");
    // for (std::set<int>::iterator it = unique_elements.begin(); it != unique_elements.end(); ++it) {
    //     printf("%d ", *it);
    // }
    // printf("\n");
    // printf("Number of unique elements: %d\n", n_unique_elements);

    // Ensure, that every atomic number is in the scaler
    for (std::set<int>::iterator it = unique_elements.begin(); it != unique_elements.end(); ++it) {
        if (std::find(qm->scaler->_fit_atom_selection, qm->scaler->_fit_atom_selection + n_unique_elements, *it) == qm->scaler->_fit_atom_selection + n_unique_elements) {
            printf("ERROR: Atomic number %d not in scaler\n", *it);
            exit(-1);
        }
    }

    // Map atomic number to the amount of its occurrences
    std::map<int, int> atomic_number_occurrences;
    for (int i = 0; i < n_atoms; i++) {
        atomic_number_occurrences[atomic_numbers[i]]++;
    }
    // printf("Atomic number occurrences: ");
    // for (std::map<int, int>::iterator it = atomic_number_occurrences.begin(); it != atomic_number_occurrences.end(); ++it) {
    //     printf("%d: %d ", it->first, it->second);
    // }
    // printf("\n");
    
    // Unscale predictions
    *energy_predictions *= qm->scaler->scale_;
    // printf("scale_ %6.3f\n", qm->scaler->scale_);
    *energy_predictions += qm->scaler->intercept_;
    // printf("intercept_ %6.3f\n", qm->scaler->intercept_);
    int atomic_number;
    for (int i=0; i<n_unique_elements; i++) {
        atomic_number = qm->scaler->_fit_atom_selection[i];
        *energy_predictions += qm->scaler->coef_[i] * atomic_number_occurrences[atomic_number];
        // printf("coef_[%d] %6.3f\n", i, qm->scaler->coef_[i]);
        // printf("atomic_number_occurrences %d [%d] %d\n",i, atomic_number, atomic_number_occurrences[atomic_number]);
    }

    for (int i=0; i<n_atoms; i++)
    {
        *(grad_predictions+(3*i+0)) *= qm->scaler->scale_;
        *(grad_predictions+(3*i+1)) *= qm->scaler->scale_;
        *(grad_predictions+(3*i+2)) *= qm->scaler->scale_;
    }

} // end of inverse_transform

/****************************************
 ******   AUXILIARY ROUTINES   **********
 ****************************************/

// adopted from src/gmxlib/pbc.c and qmmm-calculation.cpp
real calculate_distance(QMMM_QMrec* qm, const int x1_idx, const int x2_idx)
{
    rvec bond;
    for(int i=0; i<DIM; i++) {
        bond[i] = qm->xQM_get(x1_idx,i) - qm->xQM_get(x2_idx,i);
    }

    return norm(bond);
}

/* end of NN sub routines */
#endif

#pragma GCC diagnostic pop

