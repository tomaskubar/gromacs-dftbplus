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

#if GMX_QMMM_PYTORCH
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
    // Find amount of models, max 10
    char* n_models_env; // Pointer to env variable for comparison with nullpointer
    int n_models;
    if ((n_models_env = getenv("GMX_N_MODELS")) == nullptr)
    {
        printf("Amount of models to use for adaptive sampling must be given via the env variable GMX_N_MODELS\n");
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
    
    // Find and set network architecture, default "maceqeq"
    std::vector<std::string> implemented_network_architectures = {"maceqeq"};
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
        network_architecture = "maceqeq";
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
            snprintf(model_path, sizeof(model_path), "%s%d%s", model_path_prefix, i, ".model");
            printf("%s ", model_path);
        }
        printf("\n");
    }

    // Architecture differences
    int NumInputs;
    int NumOutputs;
    if (network_architecture == "maceqeq")
    {
        NumInputs = 13;
        NumOutputs = 3; // charge, energy, force
    }

    //********* Read models
    for (int model_idx=0; model_idx<n_models; model_idx++) {
        snew(qm->models[model_idx], 1); // Allocate memory for model
        printf("Reading model %d\n", model_idx);
        std::strcpy(qm->models[model_idx]->modelArchitecture, network_architecture.c_str()); // Save network architecture for later use

        snprintf(model_path, sizeof(model_path), "%s%d%s", model_path_prefix, model_idx, ".model");

        try {
            // Attempt to load the model
            torch::jit::script::Module model = torch::load(model_path);
            
            // Check if the model is not empty
            // if (model) {
                qm->models[model_idx]->model = new torch::jit::script::Module(model);
                std::cout << "Model loaded successfully." << std::endl;
            //} else {
            //    throw std::runtime_error("Model is empty after loading.");
            //}
        } catch (const c10::Error& e) {
            std::cerr << "Error loading the model: " << e.what() << std::endl;
            throw;
        } catch (const std::exception& e) {
            std::cerr << "Standard exception: " << e.what() << std::endl;
            throw;
        } catch (...) {
            std::cerr << "Unknown error occurred while loading the model." << std::endl;
            throw;
        }

        // Prepare input and output tensors
        qm->models[model_idx]->NumInputs = NumInputs;
        qm->models[model_idx]->NumOutputs = NumOutputs;
    }

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

    int n_models = std::stoi(getenv("GMX_N_MODELS"));
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
    
    // Prepare inputs
    c10::Dict<std::string, torch::Tensor> input_dict;
    if (strcmp(qm->models[0]->modelArchitecture, "maceqeq") == 0)
    {
        prepare_base_mace_inputs(qm, input_dict);
        prepare_maceqeq_inputs(qr, qm, input_dict);
    }
    
    // Check input availability outside of preparation functions
    for (const auto& entry : input_dict) {
        const std::string& key = entry.key();
        const torch::Tensor& value = entry.value();
        printf("Key: %s\n", key.c_str());
        printf("Shape: ");
        for (long unsigned int i = 0; i < value.sizes().size(); i++) {
            printf("%ld ", value.sizes()[i]);
        }
        printf("\n");
    }

    torch::Tensor mask = torch::zeros(nAtoms, torch::dtype(torch::kBool));
    for (int i=0; i<nAtoms; i++)
    {
        mask[i] = true;
    }

    // Run the Session
    std::vector<decltype(qm->models[0]->model->forward({input_dict, mask.to(qm->device)}).toGenericDict())> output_dicts;
    float* charge_predictions[n_active_models];
    float* energy_predictions[n_active_models];
    rvec* force_predictions[n_active_models];
    for (int model_idx=0; model_idx<n_active_models; model_idx++) {
        output_dicts[model_idx] = qm->models[model_idx]->model->forward({input_dict, mask.to(qm->device)}).toGenericDict();

        if (strcmp(qm->models[model_idx]->modelArchitecture, "maceqeq") == 0)
        {
            charge_predictions[model_idx] = output_dicts[model_idx].at("charges").toTensor().cpu().data_ptr<float>(); // in e-
            energy_predictions[model_idx] = output_dicts[model_idx].at("total_energy_local").toTensor().cpu().data_ptr<float>(); // in eV
            force_predictions[model_idx] = reinterpret_cast<rvec*>(output_dicts[model_idx].at("forces").toTensor().cpu().data_ptr<float>()) ; // in eV/Angstrom
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
        qm->QMcharges_set(i, (real) charge_predictions[0][i]); // sign OK
        //qm->QMcharges_set(i, charge_means[i]); // sign OK
    }
    
    QMenergy = energy_predictions[0][0];
    //QMenergy = energy_mean;

    /* Save the gradient on the QM atoms */
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<3; j++)
        {
            QMgrad[i][j] = (real) (-1*force_predictions[0][i][j]); // negative of force
            //QMgrad[i][j] = -force_means[i][j]; // negative of force 
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
        if (std::strcmp(qm->models[0]->modelArchitecture, "maceqeq") == 0)
        {
            write_maceqeq_inputs_outputs(qm, input_dict, output_dicts[0]);
        }

        FILE* f_std = nullptr; // file for saving standard deviations of significant structures
        f_std = fopen("qm_mlmm_std.xyz", "a");

        fprintf(f_std, "\nMeans of energy predictions step %d\n", step);
        fprintf(f_std, "%8.4f\n", energy_mean);

        fprintf(f_std, "\nStds of energy predictions step %d\n", step);
        fprintf(f_std, "%8.4f\n", energy_std);

        fprintf(f_std, "\nMeans of force predictions step %d\n", step);
        for (int i=0; i<n; i++) {
            fprintf(f_std, "%-2s %4i %8.4f %8.4f %8.4f\n",
                periodic_system[qm->atomicnumberQM_get(i)], i,
                force_means[i][0], force_means[i][1], force_means[i][2]);
        }
        
        fprintf(f_std, "\nStds of force predictions step %d\n", step);
        for (int i=0; i<n; i++) {
            fprintf(f_std, "%-2s %4i %8.4f %8.4f %8.4f\n",
                periodic_system[qm->atomicnumberQM_get(i)], i,
                force_stds[i][0], force_stds[i][1], force_stds[i][2]);
        }
        fclose(f_std);
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

    return (real) QMenergy * HARTREE2KJ * AVOGADRO;
} /* call_pytorch */


void prepare_base_mace_inputs(QMMM_QMrec* qm,
                              c10::Dict<std::string, torch::Tensor> input_dict)
{
    int nAtoms = qm->nrQMatoms_get();
    torch::ScalarType torch_float_dtype = qm->torch_float_dtype;
    torch::ScalarType torch_int_dtype = qm->torch_int_dtype;

    // atomic numbers
    int64_t atomic_numbers_dims[] = {nAtoms};
    torch::Tensor atomic_numbers = torch::zeros(atomic_numbers_dims, torch_int_dtype);
    for (int i=0; i<nAtoms; i++)
    {
        atomic_numbers[i] = qm->atomicnumberQM_get(i);
    }

    // coordinates
    int64_t coordinates_dims[] = {nAtoms, 3};
    torch::Tensor coordinates = torch::zeros(coordinates_dims, torch_float_dtype);
    for (int i=0; i<nAtoms; i++)
    {
        for (int j=0; j<3; j++)
        {
            coordinates[i][j] = qm->xQM_get(i,j) / BOHR2NM; // from nm to bohr for NN model
        }
    }

    // edge indices
    int nEdgeCombinations = nAtoms*(nAtoms-1); // Cheap Combinations nCr(nAtoms 2)
    int64_t edge_indices_dims[] = {2, nEdgeCombinations};
    torch::Tensor edge_indices = torch::zeros(edge_indices_dims, torch_int_dtype);
    int edge_index = 0;
    for (int i=0; i<nAtoms; i++)
    {   
        for (int j=0; j<nAtoms; j++)
        {
            if (i==j) {
                continue;
            }
            edge_indices[0][edge_index] = i;
            edge_indices[1][edge_index] = j;
            edge_index++;
        }
    }

    // box
    torch::Tensor box = torch::zeros({3,3}, torch_float_dtype);
    box[0][0] = 10.0;
    box[0][1] = 0.0;
    box[0][2] = 0.0;
    box[1][0] = 0.0;
    box[1][1] = 10.0;
    box[1][2] = 0.0;
    box[2][0] = 0.0;
    box[2][1] = 0.0;
    box[2][2] = 10.0;

    torch::Tensor batch = torch::zeros({nAtoms}, torch_int_dtype);
    torch::Tensor energy = torch::empty({1}, torch_float_dtype);
    torch::Tensor forces = torch::empty({nAtoms,3}, torch_float_dtype);
    torch::Tensor ptr = torch::empty({2}, torch_int_dtype);
    torch::Tensor weight = torch::empty({1}, torch_float_dtype);
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
    energy = energy.to(device);
    forces = forces.to(device);
    ptr = ptr.to(device);
    // shifts = shifts.to(device);
    // unit_shifts = unit_shifts.to(device);
    weight = weight.to(device);

    // Prepare input dictionary
    input_dict.insert("node_attrs", atomic_numbers);
    input_dict.insert("positions", coordinates);
    input_dict.insert("edge_index", edge_indices);
    input_dict.insert("cell", box);
    input_dict.insert("batch", batch);
    input_dict.insert("energy", energy);
    input_dict.insert("forces", forces);
    input_dict.insert("ptr", ptr);
    input_dict.insert("weight", weight);

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
    total_charge[0] = (float)qm->QMcharge_get(); // in e-

    // esp
    float V_to_au = 1/27.211386245988;
    int64_t esp_dims[] = {nAtoms};
    torch::Tensor esp = torch::zeros(esp_dims, torch_float_dtype);
    for (int i=0; i<nAtoms; i++)
    {
       esp[i] = qm->pot_qmmm_get(i)*V_to_au; // in atomic units
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
            esp_grad[i][j] = ESPgrad[i][j]; // in atomic units
        }
    }

    sfree(ESPgrad);
    if (qm->qmmm_variant_get() == eqmmmPME)
    {
        sfree(ESPgrad_full);
    }

    // Prepare input dictionary
    input_dict.insert("total_charge", total_charge);
    input_dict.insert("esp", esp);
    input_dict.insert("esp_grad", esp_grad);
} // end of prepare_maceqeq_inputs

void write_maceqeq_inputs_outputs(QMMM_QMrec* qm,
                                  c10::Dict<std::string, torch::Tensor> input_dict,
                                  c10::impl::GenericDict output_dict)
{
    // Print information and save information for debugging

    int nAtoms = qm->nrQMatoms_get();
    int nEdgeCombinations = nAtoms*(nAtoms-1);

    // atomic numbers
    FILE* f_atomic_numbers = nullptr;
    f_atomic_numbers = fopen("atomic_numbers.txt", "w");
    int64_t* atomic_numbers = input_dict.at("node_attrs").data_ptr<int64_t>();
    for (int i=0; i<nAtoms; i++)
    {
        printf("CHECK Atomic numbers %d %ld\n", i, atomic_numbers[i]);
        fprintf(f_atomic_numbers, "%ld\n", atomic_numbers[i]);
    }
    fclose(f_atomic_numbers);

    // node coordinates 
    FILE* f_coordinates = nullptr;
    f_coordinates = fopen("coordinates.txt", "w");
    rvec* coordinates = reinterpret_cast<rvec*>(input_dict.at("positions").data_ptr<float>());
    for (int i=0; i<nAtoms; i++)
    {
        printf("CHECK COORD QM Bohr[%d] = %6.3f %6.3f %6.3f\n", i+1, coordinates[i][0], coordinates[i][1], coordinates[i][2]);
        fprintf(f_coordinates, "%6.6f %6.6f %6.6f\n", coordinates[i][0], coordinates[i][1], coordinates[i][2]);
    }
    for (int i=0; i<nAtoms; i++)
    {
        printf("CHECK COORD QM MD units[%d] = %6.3f %6.3f %6.3f\n", i+1, BOHR2NM*coordinates[i][0], BOHR2NM*coordinates[i][1], BOHR2NM*coordinates[i][2]);
    }
    fclose(f_coordinates);

    // edge indices 
    FILE* f_edge_indices = nullptr;
    f_edge_indices = fopen("edge_indices.txt", "w");
    int64_t* edge_indices = input_dict.at("edge_index").data_ptr<int64_t>();
    for (int i=0; i<nEdgeCombinations; i++)
    {
        fprintf(f_edge_indices, "%ld %ld\n", edge_indices[2*i+0], edge_indices[2*i+1]);
    }
    fclose(f_edge_indices);

    // total charge
    FILE* f_total_charge = nullptr;
    f_total_charge = fopen("total_charge.txt", "w");
    float* total_charge = input_dict.at("total_charge").data_ptr<float>();
    fprintf(f_total_charge, "%1.0f\n", total_charge[0]);
    printf("CHECK Total Charge %1.0f \n", total_charge[0]);
    fclose(f_total_charge);

    // esp
    FILE* f_esp = nullptr;
    f_esp = fopen("esp.txt", "w");
    float* esp = input_dict.at("esp").data_ptr<float>();
    for (int i=0; i<nAtoms; i++)
    {
        printf("CHECK ESP QM[%d] = %6.3f\n", i+1, esp[i]);
        fprintf(f_esp, "%6.6f\n", esp[i]);
    }
    fclose(f_esp);

    // esp_grad
    FILE* f_esp_grad = nullptr;
    f_esp_grad = fopen("esp_grad.txt", "w");
    rvec* esp_grad = reinterpret_cast<rvec*>(input_dict.at("esp_grad").data_ptr<float>());
    for (int i=0; i<nAtoms; i++)
    {
        printf("CHECK ESP GRAD[%d] = %6.3f %6.3f %6.3f\n", i+1, esp_grad[i][0], esp_grad[i][1], esp_grad[i][2]);
        fprintf(f_esp_grad, "%6.6f %6.6f %6.6f\n", esp_grad[i][0], esp_grad[i][1], esp_grad[i][2]);
    }
    fclose(f_esp_grad);


    float* charge_predictions = output_dict.at("charges").toTensor().data_ptr<float>() ; // in e-
    float* energy_predictions = output_dict.at("energy").toTensor().data_ptr<float>() ; // in ev
    rvec* force_predictions = reinterpret_cast<rvec*>(output_dict.at("forces").toTensor().data_ptr<float>()) ; // in ev/Angstrom

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
        printf("CHECK QM Grad[%d] = %6.3f %6.3f %6.3f\n", i+1, force_predictions[i][0], force_predictions[i][1], force_predictions[i][2]);
        fprintf(f_output, "%6.6f %6.6f %6.6f\n", force_predictions[i][0], force_predictions[i][1], force_predictions[i][2]);
    }
    fclose(f_output);
}

/* end of NN sub routines */
#endif

#pragma GCC diagnostic pop

