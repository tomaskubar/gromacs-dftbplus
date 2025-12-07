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

#if GMX_QMMM_PYTORCH or GMX_QMMM_DFTBPLUS_PYTORCH
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

#include "gromacs/mdlib/qm_pytorch.h"
#ifdef DIM
#undef DIM
#endif
#include <torch/torch.h>
#include <torch/script.h>
#ifndef DIM
#define DIM 3
#endif

/* Pytorch interface routines */

void init_pytorch(QMMM_QMrec* qm)
{
    // Set device type
    if (!torch::cuda::is_available()) {
        printf("CUDA unavailable, setting device type to CPU.\n");
        qm->device = c10::Device(torch::kCPU);
    } else {
        printf("CUDA found, setting device type to GPU.\n");
        qm->device = c10::Device(torch::kCUDA);
    }

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
    
    // Find and set network architecture, default "maceqeq"
    std::vector<std::string> implemented_network_architectures = {"mace", "maceqeq", "amp"};
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
        network_architecture = "mace";
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
            snprintf(model_path, sizeof(model_path), "%s%d%s", model_path_prefix, i, ".pt");
            printf("%s ", model_path);
        }
        printf("\n");
    }

    //********* Read models
    torch::jit::script::Module model;
    for (int model_idx=0; model_idx<n_models; model_idx++) {
        snew(qm->models[model_idx], 1); // Allocate memory for model
        printf("Reading model %d\n", model_idx);
        std::strcpy(qm->models[model_idx]->modelArchitecture, network_architecture.c_str()); // Save network architecture for later use

        snprintf(model_path, sizeof(model_path), "%s%d%s", model_path_prefix, model_idx, ".pt");

        try {
            // Attempt to load the model
            model = torch::jit::load(model_path, qm->device);
            qm->models[model_idx]->model = new torch::jit::script::Module(model);
            printf("Model %d loaded successfully\n", model_idx);

        } catch (const c10::Error& e) {
            std::cerr << "Error loading the model: " << e.what() << std::endl;
            exit(-1);
        } catch (const std::exception& e) {
            std::cerr << "Standard exception: " << e.what() << std::endl;
            exit(-1);
        } catch (...) {
            std::cerr << "Unknown error occurred while loading the model." << std::endl;
            exit(-1);
        }
    }

    if ((strcmp(qm->models[0]->modelArchitecture, "mace") == 0) ||
        (strcmp(qm->models[0]->modelArchitecture, "maceqeq") == 0))
    {
        // extract r_max from mace model
        qm->r_max = model.attr("r_max").toTensor().item<double>();
        qm->r_max_squared = qm->r_max*qm->r_max;
        printf("Read r_max from model: %6.3f\n", qm->r_max);
        
        // extract atomic numbers from mace model
        torch::Tensor present_atomic_numbers_tensor = model.attr("atomic_numbers").toTensor();
        qm->n_present_atomic_numbers = present_atomic_numbers_tensor.size(0);
        snew (qm->present_atomic_numbers, qm->n_present_atomic_numbers);
        
        for (int i=0; i<qm->n_present_atomic_numbers; i++) {
            qm->present_atomic_numbers[i] = present_atomic_numbers_tensor[i].item<int>();
        }
        printf("Read %d atomic numbers from model\n", qm->n_present_atomic_numbers);
        for (int i=0; i<qm->n_present_atomic_numbers; i++) {
            printf("%d ", qm->present_atomic_numbers[i]);
        }
        printf("\n");
    }
    //TODO: extract default dtype from mace model

    fflush(stdout);
    fflush(stderr);
    return;
} /* init_pytorch */

real call_pytorch(QMMM_rec*       qr,
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
    c10::Dict<std::string, torch::Tensor> input_dict;
    float energy_conversion;
    float force_conversion;
    float mm_gradient_conversion;
    if (strcmp(qm->models[0]->modelArchitecture, "mace") == 0)
    {
        prepare_base_mace_inputs(qm, input_dict);
        energy_conversion = KCAL2KJ; // from kcal/mol to kJ/mol
        force_conversion = KCAL_A2MD; // from kcal/mol/A to kJ/mol/nm
        mm_gradient_conversion = HARTREE_BOHR2MD; // from Hartree/Bohr to kJ/mol/nm
    }
    else if (strcmp(qm->models[0]->modelArchitecture, "maceqeq") == 0)
    {
        prepare_base_mace_inputs(qm, input_dict);
        prepare_maceqeq_inputs(qr, qm, input_dict);
        energy_conversion = EV2KJMOL; // from eV to kJ/mol
        force_conversion = EV_A2MD; // from eV/A to kJ/mol/nm
        mm_gradient_conversion = HARTREE_BOHR2MD; // from Hartree/Bohr to kJ/mol/nm
    }
    else if (strcmp(qm->models[0]->modelArchitecture, "amp") == 0)
    {
        prepare_amp_inputs(qm, mm, input_dict);
        energy_conversion = 1; // from kJ/mol to kJ/mol
        force_conversion = 1/A2NM; // from kJ/mol/A to kJ/mol/nm
        mm_gradient_conversion = 1/A2NM; // from kJ/mol/A to kJ/mol/nm
    }
    else
    {
        printf("The network architecture %s is not implemented\n", qm->models[0]->modelArchitecture);
        exit(-1);
    }

    torch::Tensor mask = torch::zeros(nAtoms, torch::dtype(torch::kBool));
    for (int i=0; i<nAtoms; i++)
    {
        mask[i] = true;
    }

    wallcycle_stop(wcycle, ewcQM_DATA_PREP);

    // Run the Session
    wallcycle_start(wcycle, ewcQM_MODEL_CALL);
    c10::Dict<std::string, torch::Tensor> output_dicts[n_active_models];
    for (int model_idx=0; model_idx<n_active_models; model_idx++) {
        c10::Dict<c10::IValue,c10::IValue> output_dict = qm->models[model_idx]->model->forward({input_dict}).toGenericDict();
        output_dicts[model_idx] = convertDict(qm, output_dict);
    }

    wallcycle_stop(wcycle, ewcQM_MODEL_CALL);

    // Convert tensors to arrays
    torch::Tensor charge_predictions_tensor, energy_predictions_tensor, force_predictions_tensor, mm_gradient_predictions_tensor;
    float charge_predictions[n_active_models][nAtoms];
    float energy_predictions[n_active_models];
    float force_predictions[n_active_models][nAtoms][3];
    float mm_gradient_predictions[n_active_models][nMMAtoms][3];

    for (int model_idx=0; model_idx<n_active_models; model_idx++) {

        if ((strcmp(qm->models[model_idx]->modelArchitecture, "mace") == 0) || (strcmp(qm->models[model_idx]->modelArchitecture, "maceqeq") == 0))
        {
            charge_predictions_tensor = output_dicts[model_idx].at("charges").cpu(); // in e, dummy for mace
            energy_predictions_tensor = output_dicts[model_idx].at("energy").cpu(); // in eV for maceqeq, in kcal/mol for mace(for pascal for now)
            force_predictions_tensor = output_dicts[model_idx].at("forces").cpu(); // in eV/Angstrom for maceqeq, in kcal/mol/A for mace(for pascal for now), actual forces, not gradients
        }
        else if (strcmp(qm->models[model_idx]->modelArchitecture, "amp") == 0)
        {
            charge_predictions_tensor = output_dicts[model_idx].at("charges").cpu(); // dummy for amp
            energy_predictions_tensor = output_dicts[model_idx].at("energy").cpu()[0]; // in kJ
            force_predictions_tensor = output_dicts[model_idx].at("forces").cpu()[0]; // in kJ/nm
            mm_gradient_predictions_tensor = output_dicts[model_idx].at("mm_gradients").cpu()[0]; // in kJ/nm, still gradient
            for (int mm_idx=0; mm_idx<nMMAtoms; mm_idx++) {
                for (int dim_idx=0; dim_idx<DIM; dim_idx++) {
                    mm_gradient_predictions[model_idx][mm_idx][dim_idx] = mm_gradient_predictions_tensor[mm_idx][dim_idx].item<float>();
                }
            }
        }

        for (int at_idx=0; at_idx<nAtoms; at_idx++) {
            charge_predictions[model_idx][at_idx] = charge_predictions_tensor[at_idx].item<float>();
        }
        energy_predictions[model_idx] = energy_predictions_tensor.item<float>();
        for (int at_idx=0; at_idx<nAtoms; at_idx++) {
            for (int dim_idx=0; dim_idx<DIM; dim_idx++) {
                force_predictions[model_idx][at_idx][dim_idx] = force_predictions_tensor[at_idx][dim_idx].item<float>();
            }
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
        energy_mean += energy_predictions[model_idx];
    }
    energy_mean /= n_active_models;

    float energy_std = 0.0;
    for (int model_idx=0; model_idx<n_active_models; model_idx++) {
        energy_std += std::pow(energy_predictions[model_idx] - energy_mean, 2);
    }
    energy_std /= n_active_models;
    energy_std = std::sqrt(energy_std);
    if ((energy_std > energy_prediction_std_threshold) and (energy_prediction_std_threshold > 0.0))
    {
        qm->significant_structure = true;
    }

    // mean and std of gradients
    float force_means[nAtoms][3] = {};
    float force_stds[nAtoms][3] = {};
    for (int at_idx=0; at_idx<nAtoms; at_idx++) {
        for (int dim_idx=0; dim_idx<3; dim_idx++) {
            //calculate and save mean of each model prediction
            float mean = 0.0;
            for (int model_idx=0; model_idx<n_active_models; model_idx++) {
                mean += force_predictions[model_idx][at_idx][dim_idx];
            }
            mean /= n_active_models;
            force_means[at_idx][dim_idx] = mean;
            
            //calculate and save std of each model prediction
            float std = 0.0;
            for (int model_idx=0; model_idx<n_active_models; model_idx++) {
                std += std::pow(force_predictions[model_idx][at_idx][dim_idx] - mean, 2);
            }
            std /= n_active_models;
            std = std::sqrt(std);

            if ((std > force_prediction_std_threshold) and (force_prediction_std_threshold > 0.0))
            {
                qm->significant_structure = true;
            }

            force_stds[at_idx][dim_idx] = std;
        }
    }

    wallcycle_stop(wcycle, ewcQM);

    /* Save the QM charges */
    for (int i=0; i<n; i++)
    {
        //qm->QMcharges_set(i, (real) charge_predictions[0][i]); // Using the first model
        qm->QMcharges_set(i, charge_means[i]); // Using the mean of all models
    }
    
    //QMenergy = energy_predictions[0]; // Using the first model
    QMenergy = energy_mean; // Using the mean of all models

    /* Save the gradient on the QM atoms */
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<3; j++)
        {
            //QMgrad[i][j] = (real) (-1*force_predictions[0][i][j]); // negative of force // Using the first model
            QMgrad[i][j] = -1*force_means[i][j]; // negative of force // Using the mean of all models
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

    for (int j=0; j<mm.nrMMatoms; j++)
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

    if ((strcmp(qm->models[0]->modelArchitecture, "mace") == 0) || (strcmp(qm->models[0]->modelArchitecture, "maceqeq") == 0))
    {
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
    }
    else if (strcmp(qm->models[0]->modelArchitecture, "amp") == 0)
    {
        for (int mm_idx=0; mm_idx<mm.nrMMatoms; mm_idx++)
        { 
            for (int dim_idx=0; dim_idx<DIM; dim_idx++)
            {
                MMgrad[mm_idx][dim_idx] = mm_gradient_predictions[0][mm_idx][dim_idx]; // negative of force // Using the first model
                //MMgrad[mm_idx][dim_idx] = mm_gradient_predictions[0][mm_idx][dim_idx]; // negative of force // Using the mean of all models
            }
            //printf("GRAD MM %d: %8.2f %8.2f %8.2f\n", mm_idx+1, 
            //    MMgrad[mm_idx][0] , MMgrad[mm_idx][1] , MMgrad[mm_idx][2] );
        }
    }

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
    for (int i = 0; i < mm.nrMMatoms; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            f[i+qm->nrQMatoms_get()][j]      = mm_gradient_conversion*MMgrad[i][j];
         // fshift[i+qm.nrQMatoms_get()][j] = mm_gradient_conversion*MMgrad[i][j];
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
                f[i+qm->nrQMatoms_get()+mm.nrMMatoms][j]      = mm_gradient_conversion*MMgrad_full[i][j];
             // fshift[i+qm.nrQMatoms_get()+mm.nrMMatoms][j] = mm_gradient_conversion*MMgrad_full[i][j];
            }
        }
    }

    // Print inputs and outputs of the model every output_freq_extxyz steps
    if ((output_freq_extxyz > 0) && (step % output_freq_extxyz == 0))
    {
        
        if (std::strcmp(qm->models[0]->modelArchitecture, "mace") == 0)
        {
            write_base_mace_inputs_outputs(qm, input_dict, output_dicts[0], step);
        }
        else if (std::strcmp(qm->models[0]->modelArchitecture, "maceqeq") == 0)
        {
            write_maceqeq_inputs_outputs(qm, input_dict, output_dicts[0], step);
        }
        else if (std::strcmp(qm->models[0]->modelArchitecture, "amp") == 0)
        {
            write_amp_inputs_outputs(qm, mm, input_dict, output_dicts[0], step);
        }
        else
        {
            printf("The network architecture %s is not implemented for output\n", qm->models[0]->modelArchitecture);
            exit(-1);
        }
    }


    char periodic_system[37][3]={"XX",
        "H",                               "He",
        "Li","Be","B", "C", "N", "O", "F", "Ne",
        "Na","Mg","Al","Si","P", "S", "Cl","Ar",
        "K", "Ca","Sc","Ti","V", "Cr","Mn","Fe","Co",
        "Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr"
    };


    // if (qm->significant_structure)
    // {   
        // // Print inputs and outputs of the model at the significant structure
        // if (std::strcmp(qm->models[0]->modelArchitecture, "mace") == 0)
        // {
        //     write_base_mace_inputs_outputs(qm, input_dict, output_dicts[0], step);
        // }
        // else if (std::strcmp(qm->models[0]->modelArchitecture, "maceqeq") == 0)
        // {
        //     write_maceqeq_inputs_outputs(qm, input_dict, output_dicts[0], step);
        // }
        // else if (std::strcmp(qm->models[0]->modelArchitecture, "amp") == 0)
        // {
        //     write_amp_inputs_outputs(qm, mm, input_dict, output_dicts[0], step);
        // }
    // } // end significant structure

    // Save the means and stds of the predictions to a file
    if (n_active_models > 1) {
        // Save the means
        fprintf(f_std, "%d\n", n);
        fprintf(f_std, "Means step %d: %8.4f\n", step, energy_mean);

        for (int i=0; i<n; i++) {
            fprintf(f_std, "%-2s %8.4f %8.4f %8.4f\n",
                periodic_system[qm->atomicnumberQM_get(i)], 
                force_means[i][0], force_means[i][1], force_means[i][2]);
            }

        // Save the stds
        fprintf(f_std, "%d\n", n);   
        fprintf(f_std, "Stds step %d: %8.4f\n", step, energy_std);
        for (int i=0; i<n; i++) {
            fprintf(f_std, "%-2s %8.4f %8.4f %8.4f\n",
                periodic_system[qm->atomicnumberQM_get(i)],
                force_stds[i][0], force_stds[i][1], force_stds[i][2]);
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

    return (real) QMenergy * energy_conversion; // in kJ/mol
} /* call_pytorch */


void prepare_base_mace_inputs(QMMM_QMrec* qm,
                              c10::Dict<std::string, torch::Tensor> input_dict)
{
    int nAtoms = qm->nrQMatoms_get();
    torch::ScalarType torch_float_dtype = qm->torch_float_dtype;
    torch::ScalarType torch_int_dtype = qm->torch_int_dtype;
    float geometry_conversion = NM2A; // nm to angstrom

    // one-hot encode atomic numbers
    int64_t atomic_numbers_dims[] = {nAtoms, qm->n_present_atomic_numbers};
    torch::Tensor atomic_numbers = torch::zeros(atomic_numbers_dims, torch_float_dtype);
    for (int i=0; i<nAtoms; i++)
    {
        int atomic_number = qm->atomicnumberQM_get(i);
        for (int j=0; j<qm->n_present_atomic_numbers; j++)
        {
            if (atomic_number == qm->present_atomic_numbers[j])
            {
                atomic_numbers[i][j] = 1.0;
            }
        }
    }

    // coordinates
    int64_t coordinates_dims[] = {nAtoms, 3};
    torch::Tensor coordinates = torch::zeros(coordinates_dims, torch_float_dtype);
    for (int i=0; i<nAtoms; i++)
    {
        for (int j=0; j<3; j++)
        {
            coordinates[i][j] = qm->xQM_get(i,j) * geometry_conversion;
        }
    }

    // box
    torch::Tensor box = torch::zeros({3,3}, torch_float_dtype);
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            box[i][j] = qm->box_get(i,j) * geometry_conversion;
        }
    }

    // ----- edge_index, shift and unit_shift -----
    // count total number of edges
    int n_edges = 0;
    std::vector<int> n_edges_vec(nAtoms, 0);
    #pragma omp parallel for reduction(+:n_edges)
    for (int i=0; i<nAtoms; i++) {
        for (int j=0; j<nAtoms; j++) {
            if (i==j) {
                continue;
            }

            float r = calculate_distance(qm, i, j) * geometry_conversion;

            if (r > qm->r_max) {
                continue;
            }

            n_edges++;
            n_edges_vec[i] += 1;
        }
    }
  
    // make first_edge vector to help with parallelizing following loop, contains the index of the first edge for each atom
    std::vector<int> first_edge(nAtoms);
    first_edge[0] = 0;
    for (int i=0; i<nAtoms-1; i++) {
        first_edge[i+1] = first_edge[i] + n_edges_vec[i];
    }

    torch::Tensor edge_indices = torch::zeros({2, n_edges}, torch_int_dtype);
    torch::Tensor shifts = torch::zeros({n_edges,3}, torch_float_dtype);
    torch::Tensor unit_shifts = torch::zeros({n_edges,3}, torch_float_dtype);

    #pragma omp parallel for
    for (int i=0; i<nAtoms; i++) {
        int k = first_edge[i];
        for (int j=0; j<nAtoms; j++) {
            if (i==j) {
                continue;
            }

            float r = calculate_distance(qm, i, j) * geometry_conversion;
            if (r > qm->r_max) {
                continue;
            }
            
            edge_indices[0][k] = i;
            edge_indices[1][k] = j;
            k++;
        }
    }

    torch::Tensor batch = torch::zeros({nAtoms}, torch_int_dtype); // batch placeholder
    torch::Tensor ptr = torch::empty({2}, torch_int_dtype); // batch property
    torch::Tensor weight = torch::empty({1}, torch_float_dtype); // batch property
    ptr[0] = 0;
    ptr[1] = nAtoms;
    weight[0] = 1.0;

    // transfer data to device
    torch::Device device = qm->device;
    atomic_numbers = atomic_numbers.to(device);
    coordinates = coordinates.to(device);
    edge_indices = edge_indices.to(device);
    box = box.to(device);
    batch = batch.to(device);
    ptr = ptr.to(device);
    shifts = shifts.to(device);
    unit_shifts = unit_shifts.to(device);
    weight = weight.to(device);

    // Prepare input dictionary
    input_dict.insert("node_attrs", atomic_numbers);
    input_dict.insert("positions", coordinates);
    input_dict.insert("edge_index", edge_indices);
    input_dict.insert("cell", box);
    input_dict.insert("batch", batch);
    input_dict.insert("ptr", ptr);
    input_dict.insert("weight", weight);
    input_dict.insert("shifts", shifts);
    input_dict.insert("unit_shifts", unit_shifts);

} // end of base_mace_inputs


void prepare_maceqeq_inputs(QMMM_rec* qr,
                          QMMM_QMrec* qm,
                          c10::Dict<std::string, torch::Tensor> input_dict)
{
    int nAtoms = qm->nrQMatoms_get();
    torch::ScalarType torch_float_dtype = qm->torch_float_dtype;
    // torch::ScalarType torch_int_dtype = qm->torch_int_dtype;

    // total charge
    int64_t total_charge_dims[] = {1};
    torch::Tensor total_charge = torch::zeros(total_charge_dims, torch_float_dtype);
    total_charge[0] = qm->QMcharge_get(); // in e-

    // esp
    int64_t esp_dims[] = {nAtoms};
    torch::Tensor esp = torch::zeros(esp_dims, torch_float_dtype);
    for (int i=0; i<nAtoms; i++)
    {
       esp[i] = qm->pot_qmmm_get(i); // in Volt units, could be changed with *V2AU
    }

    // esp_grad
    rvec *ESPgrad = nullptr, *ESPgrad_full = nullptr;
    snew(ESPgrad, nAtoms);
    if (qm->qmmm_variant_get() == eqmmmPME)
    {
        snew(ESPgrad_full, nAtoms);
    }
    int64_t esp_grad_dims[] = {nAtoms,3};
    torch::Tensor esp_grad = torch::zeros(esp_grad_dims, torch_float_dtype);
    // qr->gradient_ESP(//cr, nrnb, wcycle, nullptr, qm->qmmm_variant_get(), ESPgrad, ESPgrad_full);
    qr->gradient_ESP(qm->qmmm_variant_get(), ESPgrad, ESPgrad_full);
    for (int i=0; i<nAtoms; i++)
    {
        for (int j=0; j<3; j++)
        {
            esp_grad[i][j] = ESPgrad[i][j] *AU2EV/BOHR2NM/10 ; // from H/B to eV/e/Angstrom
        }
    }

    sfree(ESPgrad);
    if (qm->qmmm_variant_get() == eqmmmPME)
    {
        sfree(ESPgrad_full);
    }

    // transfer data to device
    torch::Device device = qm->device;
    total_charge = total_charge.to(device);
    esp = esp.to(device);
    esp_grad = esp_grad.to(device);

    // Prepare input dictionary
    input_dict.insert("total_charge", total_charge);
    input_dict.insert("esp", esp);
    input_dict.insert("esp_gradient", esp_grad);
} // end of prepare_maceqeq_inputs

void prepare_amp_inputs(QMMM_QMrec* qm,
                        QMMM_MMrec mm,
                        c10::Dict<std::string, torch::Tensor> input_dict   
                        )
{
    int nQMAtoms = qm->nrQMatoms_get();
    int nMMAtoms = mm.nrMMatoms;
    torch::ScalarType torch_float_dtype = qm->torch_float_dtype;
    torch::ScalarType torch_int_dtype = qm->torch_int_dtype;
    float geometry_conversion = NM2A; // nm to angstrom

    // atomic numbers
    int64_t atomic_numbers_dims[] = {1, nQMAtoms};
    torch::Tensor atomic_numbers = torch::zeros(atomic_numbers_dims, torch_int_dtype);
    for (int i=0; i<nQMAtoms; i++)
    {
        int atomic_number = qm->atomicnumberQM_get(i);
        atomic_numbers[0][i] = atomic_number;
    }

    // qm_coordinates
    int64_t qm_coordinates_dims[] = {1, nQMAtoms, 3};
    torch::Tensor coordinates = torch::zeros(qm_coordinates_dims, torch_float_dtype);
    for (int i=0; i<nQMAtoms; i++)
    {
        for (int j=0; j<3; j++)
        {
            coordinates[0][i][j] = qm->xQM_get(i,j) * geometry_conversion;
        }
    }

    // mm_charges
    int buffered_nMMAtoms = std::max(1, nMMAtoms); // buffer for 1 MM atom
    int64_t mm_charges_dims[] = {1, buffered_nMMAtoms};
    torch::Tensor mm_charges = torch::zeros(mm_charges_dims, torch_float_dtype);
    for (int i=0; i<nMMAtoms; i++)
    {
        mm_charges[0][i] = mm.MMcharges[i];
    }

    // mm_coordinates
    int64_t mm_coordinates_dims[] = {1, buffered_nMMAtoms, 3};
    torch::Tensor mm_coordinates = torch::zeros(mm_coordinates_dims, torch_float_dtype);
    for (int i=0; i<nMMAtoms; i++)
    {
        for (int j=0; j<3; j++)
        {
            mm_coordinates[0][i][j] = mm.xMM[i][j] * geometry_conversion;
        }
    }

    // transfer data to device
    torch::Device device = qm->device;

    atomic_numbers = atomic_numbers.to(device);
    coordinates = coordinates.to(device);
    mm_charges = mm_charges.to(device);
    mm_coordinates = mm_coordinates.to(device);

    // Prepare input dictionary
    input_dict.insert("qm_charges", atomic_numbers);
    input_dict.insert("qm_coordinates", coordinates);
    input_dict.insert("mm_charges", mm_charges);
    input_dict.insert("mm_coordinates", mm_coordinates);
} // end of prepare_amp_inputs


void write_base_mace_inputs_outputs(QMMM_QMrec* qm,
                                  c10::Dict<std::string, torch::Tensor> input_dict,
                                  c10::Dict<std::string, torch::Tensor> output_dict,
                                  int step)
{
    // Print information and save information for debugging

    int nAtoms = qm->nrQMatoms_get();
    float x, y, z;

    // // atomic numbers
    // FILE* f_atomic_numbers = nullptr;
    // f_atomic_numbers = fopen("atomic_numbers.txt", "w");
    // torch::Tensor atomic_numbers = input_dict.at("node_attrs").cpu();
    // // printf("CHECK Atomic numbers size %ld %ld\n", atomic_numbers.size(0), atomic_numbers.size(1));
    // for (int i=0; i<nAtoms; i++)
    // {
    //     // printf("CHECK Atomic numbers %d", i+1);
    //     for (int j=0; j<qm->n_present_atomic_numbers; j++)
    //     {
    //         // printf(" %6.3f", atomic_numbers[i][j].item<float>());
    //         fprintf(f_atomic_numbers, "%6.3f ", atomic_numbers[i][j].item<float>());
    //     }
    //     // printf("\n");
    //     fprintf(f_atomic_numbers, "\n");
    // }
    // fclose(f_atomic_numbers);

    // node coordinates 
    // FILE* f_coordinates = nullptr;
    // f_coordinates = fopen("coordinates.txt", "w");
    torch::Tensor coordinates_tensor = input_dict.at("positions").cpu();
    // printf("CHECK COORD QM size %ld %ld\n", coordinates_tensor.size(0), coordinates_tensor.size(1));
    float coordinates[nAtoms][3];
    for (int i=0; i<nAtoms; i++)
    {
        x = coordinates_tensor[i][0].item<float>();
        y = coordinates_tensor[i][1].item<float>();
        z = coordinates_tensor[i][2].item<float>();
        coordinates[i][0] = x;
        coordinates[i][1] = y;
        coordinates[i][2] = z;
        // printf("CHECK COORD QM Angstrom[%d] = %6.3f %6.3f %6.3f\n", i+1, x, y, z);
        // fprintf(f_coordinates, "%6.6f %6.6f %6.6f\n", x, y, z);
    }
    // for (int i=0; i<nAtoms; i++)
    // {
    //     x = coordinates[i][0];
    //     y = coordinates[i][1];
    //     z = coordinates[i][2];
    //     printf("CHECK COORD QM MD units[%d] = %6.3f %6.3f %6.3f\n", i+1, A2NM*x, A2NM*y, A2NM*z);
    // }
    // fclose(f_coordinates);

    // // edge indices 
    // FILE* f_edge_indices = nullptr;
    // f_edge_indices = fopen("edge_indices.txt", "w");
    // torch::Tensor edge_indices = input_dict.at("edge_index").cpu();
    // int n_edges = edge_indices.size(1);
    // for (int i=0; i<n_edges; i++)
    // {
    //     //printf("CHECK EDGE INDICES[%d] = %ld %ld\n", i+1, edge_indices[0][i].item<int64_t>(), edge_indices[1][i].item<int64_t>());
    //     fprintf(f_edge_indices, "%ld %ld\n", edge_indices[0][i].item<int64_t>(), edge_indices[1][i].item<int64_t>());
    // }
    // fclose(f_edge_indices);

    // // batch, shifts, unit_shifts, weight
    // FILE* f_batch = nullptr;
    // f_batch = fopen("batch.txt", "w");
    // torch::Tensor batch_tensor = input_dict.at("batch").cpu();
    // int64_t batch[nAtoms];
    // for (int i=0; i<nAtoms; i++)
    // {
    //     batch[i] = batch_tensor[i].item<int64_t>();
    //     //printf("CHECK BATCH[%d] = %ld\n", i+1, batch[i]);
    //     fprintf(f_batch, "%ld\n", batch[i]);
    // }
    // fclose(f_batch);

    // FILE* f_shift = nullptr;
    // f_shift = fopen("shifts.txt", "w");
    // torch::Tensor shifts_tensor = input_dict.at("shifts").cpu();
    // for (int i=0; i<n_edges; i++)
    // {
    //     x = shifts_tensor[i][0].item<float>();
    //     y = shifts_tensor[i][1].item<float>();
    //     z = shifts_tensor[i][2].item<float>();
    //     //printf("CHECK SHIFT[%d] = %6.3f %6.3f %6.3f\n", i+1, x, y, z);
    //     fprintf(f_shift, "%6.6f %6.6f %6.6f\n", x, y, z);
    // }
    // fclose(f_shift);

    // FILE* f_unit_shift = nullptr;
    // f_unit_shift = fopen("unit_shifts.txt", "w");
    // torch::Tensor unit_shifts_tensor = input_dict.at("unit_shifts").cpu();
    // for (int i=0; i<n_edges; i++)
    // {
    //     x = unit_shifts_tensor[i][0].item<float>();
    //     y = unit_shifts_tensor[i][1].item<float>();
    //     z = unit_shifts_tensor[i][2].item<float>();
    //     //printf("CHECK UNIT SHIFT[%d] = %6.3f %6.3f %6.3f\n", i+1, x, y, z);
    //     fprintf(f_unit_shift, "%6.6f %6.6f %6.6f\n", x, y, z);
    // }
    // fclose(f_unit_shift);

    // FILE* f_weight = nullptr;
    // f_weight = fopen("weight.txt", "w");
    // float weight = input_dict.at("weight").cpu().item<float>();
    // //printf("CHECK WEIGHT = %6.3f\n", weight);
    // fprintf(f_weight, "%6.6f\n", weight);
    // fclose(f_weight);

    torch::Tensor box_tensor = input_dict.at("cell").cpu();
    float box[3][3];
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            box[i][j] = box_tensor[i][j].item<float>();
        }
    }

    torch::Tensor energy_predictions_tensor = output_dict.at("energy").cpu() ;
    torch::Tensor force_predictions_tensor = output_dict.at("forces").cpu() ;
    float energy_predictions;
    float force_predictions[nAtoms][3];

    // FILE* f_output = nullptr;
    // f_output = fopen("output.txt", "w");

    energy_predictions = energy_predictions_tensor[0].item<float>();
    //printf("CHECK QM Energy = %6.3f\n", energy_predictions);
    // fprintf(f_output, "%6.6f\n", energy_predictions);

    for (int i=0; i<nAtoms; i++)
    {
        x = force_predictions_tensor[i][0].item<float>();
        y = force_predictions_tensor[i][1].item<float>();
        z = force_predictions_tensor[i][2].item<float>();
        force_predictions[i][0] = x;
        force_predictions[i][1] = y;
        force_predictions[i][2] = z;
        // printf("CHECK QM Grad[%d] = %6.3f %6.3f %6.3f\n", i+1, x, y, z);
        // fprintf(f_output, "%6.6f %6.6f %6.6f\n", x, y, z);
    }
    // fclose(f_output);

    // Only proceed to write extended xyz if the model is mace
    if (std::strcmp(qm->models[0]->modelArchitecture, "mace") != 0)
    {
        return;
    }

    // Write all inputs and outputs in the extended xyz format
    char periodic_system[37][3]={"XX",
        "H",                               "He",
        "Li","Be","B", "C", "N", "O", "F", "Ne",
        "Na","Mg","Al","Si","P", "S", "Cl","Ar",
        "K", "Ca","Sc","Ti","V", "Cr","Mn","Fe","Co",
        "Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr"};
    FILE* f_extxyz = nullptr;
    f_extxyz = fopen("qm_mlmm.extxyz", "a");
    fprintf(f_extxyz, "%d\n", nAtoms);
    fprintf(f_extxyz, "Lattice=\"%.1f %.1f %.1f %.1f %.1f %.1f %.1f %.1f %.1f\" ", box[0][0], box[0][1], box[0][2], box[1][0], box[1][1], box[1][2], box[2][0], box[2][1], box[2][2]);
    fprintf(f_extxyz, "Properties=species:S:1:pos:R:3:gromacs_force:R:3 ");
    fprintf(f_extxyz, "gromacs_energy=%.6f ", energy_predictions);
    fprintf(f_extxyz, "comment=\"Step %d\" ", step);
    fprintf(f_extxyz, "pbc=\"F F F\" \n");
    for (int i=0; i<nAtoms; i++)
    {
        fprintf(f_extxyz, "%-2s ", periodic_system[qm->atomicnumberQM_get(i)]);
        fprintf(f_extxyz, "%8.4f %8.4f %8.4f ", coordinates[i][0], coordinates[i][1], coordinates[i][2]);
        fprintf(f_extxyz, "%8.4f %8.4f %8.4f ", force_predictions[i][0], force_predictions[i][1], force_predictions[i][2]);
        fprintf(f_extxyz, "\n");
    }
    fclose(f_extxyz);
} // end of write_base_mace_inputs_outputs
    

void write_maceqeq_inputs_outputs(QMMM_QMrec* qm,
                                  c10::Dict<std::string, torch::Tensor> input_dict,
                                  c10::Dict<std::string, torch::Tensor> output_dict,
                                  int step)
{
    // Print information and save information for debugging

    int nAtoms = qm->nrQMatoms_get();
    float x, y, z;

    // Call to write_base_mace_inputs_outputs
    // write_base_mace_inputs_outputs(qm, input_dict, output_dict, step);

    // node coordinates, printed in write_base_mace_inputs_outputs
    torch::Tensor coordinates_tensor = input_dict.at("positions").cpu();
    float coordinates[nAtoms][3];
    for (int i=0; i<nAtoms; i++)
    {
        x = coordinates_tensor[i][0].item<float>();
        y = coordinates_tensor[i][1].item<float>();
        z = coordinates_tensor[i][2].item<float>();
        coordinates[i][0] = x;
        coordinates[i][1] = y;
        coordinates[i][2] = z;
    }

    // total charge
    // FILE* f_total_charge = nullptr;
    // f_total_charge = fopen("total_charge.txt", "w");
    double total_charge = input_dict.at("total_charge").cpu().item<double>();
    // fprintf(f_total_charge, "%1.0f\n", total_charge);
    // printf("CHECK Total Charge %1.0f \n", total_charge);
    // fclose(f_total_charge);

    // esp
    // FILE* f_esp = nullptr;
    // f_esp = fopen("esp.txt", "w");
    torch::Tensor esp_tensor = input_dict.at("esp").cpu();
    float esp[nAtoms];
    for (int i=0; i<nAtoms; i++)
    {
        esp[i] = esp_tensor[i].item<float>();
        // printf("CHECK ESP QM[%d] = %6.3f\n", i+1, esp[i]);
        // fprintf(f_esp, "%6.6f\n", esp[i]);
    }
    // fclose(f_esp);

    // esp_grad
    // FILE* f_esp_grad = nullptr;
    // f_esp_grad = fopen("esp_grad.txt", "w");
    torch::Tensor esp_grad_tensor = input_dict.at("esp_gradient").cpu();
    float esp_grad[nAtoms][3];
    for (int i=0; i<nAtoms; i++)
    {
        x = esp_grad_tensor[i][0].item<float>();
        y = esp_grad_tensor[i][1].item<float>();
        z = esp_grad_tensor[i][2].item<float>();
        esp_grad[i][0] = x;
        esp_grad[i][1] = y;
        esp_grad[i][2] = z;
        // printf("CHECK ESP GRAD[%d] = %6.3f %6.3f %6.3f\n", i+1, x, y, z);
        // fprintf(f_esp_grad, "%6.6f %6.6f %6.6f\n", x, y, z);
    }
    // fclose(f_esp_grad);

    torch::Tensor charge_predictions_tensor = output_dict.at("charges").cpu() ; // in e-
    torch::Tensor energy_predictions_tensor = output_dict.at("energy").cpu() ; // in ev
    torch::Tensor force_predictions_tensor = output_dict.at("forces").cpu() ; // in ev/Angstrom
    float charge_predictions[nAtoms];
    float energy_predictions;
    float force_predictions[nAtoms][3];

    // FILE* f_output = nullptr;
    // f_output = fopen("output.txt", "w"); // Overwrites the write_base_mace_inputs_outputs output.txt
    for (int i=0; i<nAtoms; i++)
    {
        charge_predictions[i] = charge_predictions_tensor[i].item<float>();
        // printf("CHECK CHARGE QM[%d] = %6.3f\n", i+1, charge_predictions[i]);
        // fprintf(f_output, "%6.6f\n", charge_predictions[i]);
    }

    energy_predictions = energy_predictions_tensor[0].item<float>(); // printed in write_base_mace_inputs_outputs
    // fprintf(f_output, "%6.6f\n", energy_predictions);

    for (int i=0; i<nAtoms; i++)
    {
        x = force_predictions_tensor[i][0].item<float>();
        y = force_predictions_tensor[i][1].item<float>();
        z = force_predictions_tensor[i][2].item<float>();
        force_predictions[i][0] = x;
        force_predictions[i][1] = y;
        force_predictions[i][2] = z;
        // fprintf(f_output, "%6.6f %6.6f %6.6f\n", x, y, z);
    }
    // fclose(f_output);

    torch::Tensor box_tensor = input_dict.at("cell").cpu();
    float box[3][3];
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            box[i][j] = box_tensor[i][j].item<float>();
        }
    }

    // Write all inputs and outputs in the extended xyz format
    char periodic_system[37][3]={"XX",
        "H",                               "He",
        "Li","Be","B", "C", "N", "O", "F", "Ne",
        "Na","Mg","Al","Si","P", "S", "Cl","Ar",
        "K", "Ca","Sc","Ti","V", "Cr","Mn","Fe","Co",
        "Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr"};
    FILE* f_extxyz = nullptr;
    f_extxyz = fopen("qm_mlmm.extxyz", "a"); // Write to the same file as in write_base_mace_inputs_outputs, TODO: Separate files
    fprintf(f_extxyz, "%d\n", nAtoms);
    fprintf(f_extxyz, "Lattice=\"%.1f %.1f %.1f %.1f %.1f %.1f %.1f %.1f %.1f\" ", box[0][0], box[0][1], box[0][2], box[1][0], box[1][1], box[1][2], box[2][0], box[2][1], box[2][2]);
    fprintf(f_extxyz, "Properties=species:S:1:pos:R:3:gromacs_force:R:3:gromacs_charge:R:1:esp:R:1:esp_gradient:R:3 ");
    fprintf(f_extxyz, "gromacs_energy=%.6f ", energy_predictions);
    fprintf(f_extxyz, "total_charge=%.1f ", total_charge);
    fprintf(f_extxyz, "comment=\"Step %d\" ", step);
    fprintf(f_extxyz, "pbc=\"F F F\" \n");
    for (int i=0; i<nAtoms; i++)
    {
        fprintf(f_extxyz, "%-2s ", periodic_system[qm->atomicnumberQM_get(i)]);
        fprintf(f_extxyz, "%8.4f %8.4f %8.4f ", coordinates[i][0], coordinates[i][1], coordinates[i][2]);
        fprintf(f_extxyz, "%8.4f %8.4f %8.4f ", force_predictions[i][0], force_predictions[i][1], force_predictions[i][2]);
        fprintf(f_extxyz, "%8.4f ", charge_predictions[i]);
        fprintf(f_extxyz, "%8.4f ", esp[i]);
        fprintf(f_extxyz, "%8.4f %8.4f %8.4f", esp_grad[i][0], esp_grad[i][1], esp_grad[i][2]);
        fprintf(f_extxyz, "\n");
    }
    fclose(f_extxyz);
} // end of write_maceqeq_inputs_outputs

void write_amp_inputs_outputs(QMMM_QMrec* qm,
                                QMMM_MMrec mm,
                                c10::Dict<std::string, torch::Tensor> input_dict,
                                c10::Dict<std::string, torch::Tensor> output_dict,
                                int step)
{
    // Print information and save information for debugging

    int nQMAtoms = qm->nrQMatoms_get();
    int nMMAtoms = std::max(1, mm.nrMMatoms); // buffer for 1 MM atom
    // printf("CHECK nQMAtoms = %d\n", nQMAtoms);
    // printf("CHECK nMMAtoms = %d\n", nMMAtoms);
    float x, y, z, mm_charge;

    // // atomic numbers
    // FILE* f_atomic_numbers = nullptr;
    // f_atomic_numbers = fopen("atomic_numbers.txt", "w");
    // torch::Tensor atomic_numbers = input_dict.at("qm_charges").cpu();
    // // printf("CHECK Atomic numbers size %ld\n", atomic_numbers.size(0));
    // for (int i=0; i<nQMAtoms; i++)
    // {
    //     // printf("CHECK Atomic numbers %d: %2.0f", i+1, atomic_numbers[i].item<float>());
    //     fprintf(f_atomic_numbers, "%6.3f ", atomic_numbers[0][i].item<float>());
    //     // printf("\n");
    //     fprintf(f_atomic_numbers, "\n");
    // }
    // fclose(f_atomic_numbers);

    // node coordinates 
    // FILE* f_coordinates = nullptr;
    // f_coordinates = fopen("coordinates.txt", "w");
    torch::Tensor coordinates_tensor = input_dict.at("qm_coordinates").cpu();
    // printf("CHECK COORD QM size %ld %ld\n", coordinates_tensor.size(0), coordinates_tensor.size(1));
    float coordinates[nQMAtoms][3];
    for (int i=0; i<nQMAtoms; i++)
    {
        x = coordinates_tensor[0][i][0].item<float>();
        y = coordinates_tensor[0][i][1].item<float>();
        z = coordinates_tensor[0][i][2].item<float>();
        coordinates[i][0] = x;
        coordinates[i][1] = y;
        coordinates[i][2] = z;
        // printf("CHECK COORD QM Angstrom[%d] = %6.3f %6.3f %6.3f\n", i+1, x, y, z);
        // fprintf(f_coordinates, "%6.6f %6.6f %6.6f\n", x, y, z);
    }
    for (int i=0; i<nQMAtoms; i++)
    {
        x = coordinates[i][0];
        y = coordinates[i][1];
        z = coordinates[i][2];
        // printf("CHECK COORD QM MD units[%d] = %6.3f %6.3f %6.3f\n", i+1, A2NM*x, A2NM*y, A2NM*z);
    }
    // fclose(f_coordinates);

    // mm coordinates
    // FILE* f_mm_coordinates = nullptr;
    // f_mm_coordinates = fopen("mm_coordinates.txt", "w");
    torch::Tensor mm_coordinates_tensor = input_dict.at("mm_coordinates").cpu();
    // printf("CHECK COORD MM size %ld %ld %ld\n", mm_coordinates_tensor.size(0), mm_coordinates_tensor.size(1), mm_coordinates_tensor.size(2));
    float mm_coordinates[nMMAtoms][3];
    for (int i=0; i<nMMAtoms; i++)
    {
        x = mm_coordinates_tensor[0][i][0].item<float>();
        y = mm_coordinates_tensor[0][i][1].item<float>();
        z = mm_coordinates_tensor[0][i][2].item<float>();
        mm_coordinates[i][0] = x;
        mm_coordinates[i][1] = y;
        mm_coordinates[i][2] = z;
        // printf("CHECK COORD MM Angstrom[%d] = %6.3f %6.3f %6.3f\n", i+1, x, y, z);
        // fprintf(f_mm_coordinates, "%6.6f %6.6f %6.6f\n", x, y, z);
    }
    for (int i=0; i<nMMAtoms; i++)
    {
        x = mm_coordinates[i][0];
        y = mm_coordinates[i][1];
        z = mm_coordinates[i][2];
        // printf("CHECK COORD MM MD units[%d] = %6.3f %6.3f %6.3f\n", i+1, A2NM*x, A2NM*y, A2NM*z);
    }
    // fclose(f_mm_coordinates);

    // mm charges
    // FILE* f_mm_charges = nullptr;
    // f_mm_charges = fopen("mm_charges.txt", "w");
    torch::Tensor mm_charges_tensor = input_dict.at("mm_charges").cpu();
    // printf("CHECK MM charges size %ld %ld\n", mm_charges_tensor.size(0), mm_charges_tensor.size(1));
    float mm_charges[nMMAtoms];
    for (int i=0; i<nMMAtoms; i++)
    {
        mm_charges[i] = mm_charges_tensor[0][i].item<float>();
        // printf("CHECK MM charges %d = %6.3f\n", i+1, mm_charges[i]);
        // fprintf(f_mm_charges, "%6.6f\n", mm_charges[i]);
    }
    // fclose(f_mm_charges);

    // mm_charges + coordinates
    FILE* f_pcfile = nullptr;
    f_pcfile = fopen("qm_mlmm.pc", "a");
    fprintf(f_pcfile, "%d\n", nMMAtoms);
    for (int i=0; i<nMMAtoms; i++)
    {
        mm_charge = mm_charges[i];
        x = mm_coordinates[i][0];
        y = mm_coordinates[i][1];
        z = mm_coordinates[i][2];
        fprintf(f_pcfile, "%6.4f %6.8f %6.8f %6.8f\n", mm_charge, x, y, z);
    }
    fclose(f_pcfile);

    // output dictionary
    // FILE* f_output = nullptr;
    FILE* f_mm_output = nullptr;
    // f_output = fopen("output.txt", "w");
    f_mm_output = fopen("qm_mlmm.pcgrad", "a");
    torch::Tensor energy_predictions_tensor = output_dict.at("energy").cpu()[0] ; // in kJ/mol
    torch::Tensor force_predictions_tensor = output_dict.at("forces").cpu()[0] ; // from kJ/mol/Angstrom to H/B
    torch::Tensor dipole_predictions_tensor = output_dict.at("dipoles").cpu()[0] ; // in ev/Angstrom
    torch::Tensor quadrupole_predictions_tensor = output_dict.at("quadrupoles").cpu()[0] ; // in ev/Angstrom^2
    torch::Tensor mm_gradients_tensor = output_dict.at("mm_gradients").cpu()[0] ; // from kJ/mol/Angstrom to H/B
    
    float energy_conversion = 1/EV2KJMOL; // kJ/mol to eV
    float dipole_conversion = 1/0.2081943; // e*Angstrom to Debye
    float quadrupole_conversion = A2BOHR*A2BOHR; // e*Angstrom^2 to e*Bohr^2
    float force_conversion_extxyz = (1/EV2KJMOL) / 1; // kJ/mol/Angstrom to eV/Angstrom
    float force_conversion_pcgrad =  0.00038088 / A2BOHR; // kJ/mol/Angstrom to Hartree/Bohr

    float energy_predictions = energy_predictions_tensor[0].item<float>();
    // printf("CHECK QM Energy = %6.3f\n", energy_predictions);
    // fprintf(f_output, "%6.6f\n", energy_predictions);

    float force_predictions[nQMAtoms][3];
    for (int i=0; i<nQMAtoms; i++)
    {
        x = force_predictions_tensor[i][0].item<float>();
        y = force_predictions_tensor[i][1].item<float>();
        z = force_predictions_tensor[i][2].item<float>();
        force_predictions[i][0] = x;
        force_predictions[i][1] = y;
        force_predictions[i][2] = z;
        // printf("CHECK QM Grad[%d] = %6.3f %6.3f %6.3f\n", i+1, x, y, z);
        // fprintf(f_output, "%6.8f %6.8f %6.8f\n", x, y, z);
    }

    float dipole_predictions[3];
    for (int i=0; i<3; i++)
    {
        dipole_predictions[i] = dipole_predictions_tensor[i].item<float>();
    }
    // printf("CHECK QM Dipole = %6.3f %6.3f %6.3f\n", dipole_predictions[0], dipole_predictions[1], dipole_predictions[2]);
    // fprintf(f_output, "%6.6f %6.6f %6.6f\n", dipole_predictions[0], dipole_predictions[1], dipole_predictions[2]);

    float quadrupole_predictions[3][3];
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<3; j++)
        {
            quadrupole_predictions[i][j] = quadrupole_predictions_tensor[i][j].item<float>();
        }
        // printf("CHECK QM Quadrupole[%d] = %6.3f %6.3f %6.3f\n", i+1, quadrupole_predictions[i][0], quadrupole_predictions[i][1], quadrupole_predictions[i][2]);
        // fprintf(f_output, "%6.6f %6.6f %6.6f\n", quadrupole_predictions[i][0], quadrupole_predictions[i][1], quadrupole_predictions[i][2]);
    }
    // fclose(f_output);

    // mm gradients
    // printf("CHECK MM Gradients size %ld %ld\n", mm_gradients_tensor.size(0), mm_gradients_tensor.size(1));
    fprintf(f_mm_output, "%d\n", nMMAtoms);
    for (int i=0; i<nMMAtoms; i++)
    {
        x = mm_gradients_tensor[i][0].item<float>();
        y = mm_gradients_tensor[i][1].item<float>();
        z = mm_gradients_tensor[i][2].item<float>();
        // printf("CHECK MM Grad[%d] = %6.3f %6.3f %6.3f\n", i+1, x, y, z);
        x = x * force_conversion_pcgrad;
        y = y * force_conversion_pcgrad;
        z = z * force_conversion_pcgrad;
        // printf("CHECK MM Grad[%d] = %6.3f %6.3f %6.3f\n", i+1, x, y, z);
        fprintf(f_mm_output, "%6.8f %6.8f %6.8f\n", x, y, z);
    }
    fclose(f_mm_output);

    // Write all inputs and outputs in the extended xyz format
    char periodic_system[37][3]={"XX",
        "H",                               "He",
        "Li","Be","B", "C", "N", "O", "F", "Ne",
        "Na","Mg","Al","Si","P", "S", "Cl","Ar",
        "K", "Ca","Sc","Ti","V", "Cr","Mn","Fe","Co",
        "Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr"};
    FILE* f_extxyz = nullptr;
    f_extxyz = fopen("qm_mlmm.extxyz", "a");
    fprintf(f_extxyz, "%d\n", nQMAtoms);
    fprintf(f_extxyz, "Properties=species:S:1:pos:R:3:gromacs_force:R:3 ");
    fprintf(f_extxyz, "gromacs_energy=%.6f ", energy_predictions*energy_conversion);
    fprintf(f_extxyz, "gromacs_dipole='%.3f %.3f %.3f' ", 
            dipole_predictions[0] * dipole_conversion, 
            dipole_predictions[1] * dipole_conversion, 
            dipole_predictions[2] * dipole_conversion);
    fprintf(f_extxyz, "gromacs_quadrupole='%.3f %.3f %.3f %.3f %.3f %.3f' ", 
            quadrupole_predictions[0][0] * quadrupole_conversion, 
            quadrupole_predictions[1][1] * quadrupole_conversion, 
            quadrupole_predictions[2][2] * quadrupole_conversion, 
            quadrupole_predictions[0][1] * quadrupole_conversion, 
            quadrupole_predictions[0][2] * quadrupole_conversion, 
            quadrupole_predictions[1][2] * quadrupole_conversion);
    fprintf(f_extxyz, "comment=\"Step %d\"\n", step);
    for (int i=0; i<nQMAtoms; i++)
    {
        fprintf(f_extxyz, "%-2s ", periodic_system[qm->atomicnumberQM_get(i)]);
        fprintf(f_extxyz, "%6.8f %6.8f %6.8f ", coordinates[i][0], coordinates[i][1], coordinates[i][2]);
        fprintf(f_extxyz, "%6.8f %6.8f %6.8f",
            force_predictions[i][0] * force_conversion_extxyz,
            force_predictions[i][1] * force_conversion_extxyz,
            force_predictions[i][2] * force_conversion_extxyz);
        fprintf(f_extxyz, "\n");
    }
    fclose(f_extxyz);
}

c10::Dict<std::string, torch::Tensor> convertDict(QMMM_QMrec* qm, const c10::impl::GenericDict &inputDict) {
    c10::Dict<std::string, torch::Tensor> outputDict;

    for (const auto& item : inputDict) {
        if (item.value().isTensor()) { // to avoid missing values such as stress/virials
            // Convert key to std::string
            std::string key = item.key().toStringRef();

            // Convert value to torch::Tensor
            torch::Tensor value = item.value().toTensor();

            // Insert into the new dictionary
            outputDict.insert(key, value);
        } 
    }

    // Add expected keys, if they are not present in the inputDict
    int nAtoms = qm->nrQMatoms_get();
    if (!outputDict.contains("charges")) {
        outputDict.insert("charges", torch::zeros({nAtoms, 1}, qm->torch_float_dtype));
    }

    return outputDict;
} // end of convertDict

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

