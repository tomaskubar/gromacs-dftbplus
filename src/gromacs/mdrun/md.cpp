/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2011-2019,2020, by the GROMACS development team, led by
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
/*! \internal \file
 *
 * \brief Implements the integrator for normal molecular dynamics simulations
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_mdrun
 */
#include "gmxpre.h"

#include <cinttypes>
#include <cmath>
#include <cstdio>
#include <cstdlib>
/* TRANSFER */
#include <ctime>
/* END TRANSFER */

#include <algorithm>
#include <memory>
#include <numeric>

#include "gromacs/applied_forces/awh/awh.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/domdec/collect.h"
#include "gromacs/domdec/dlbtiming.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_network.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/gpuhaloexchange.h"
#include "gromacs/domdec/mdsetup.h"
#include "gromacs/domdec/partition.h"
#include "gromacs/essentialdynamics/edsam.h"
#include "gromacs/ewald/pme_load_balancing.h"
#include "gromacs/ewald/pme_pp.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/gpu_utils/device_stream_manager.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/imd/imd.h"
#include "gromacs/listed_forces/listed_forces.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/invertmatrix.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/checkpointhandler.h"
#include "gromacs/mdlib/compute_io.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdlib/coupling.h"
#include "gromacs/mdlib/ebin.h"
#include "gromacs/mdlib/enerdata_utils.h"
#include "gromacs/mdlib/energyoutput.h"
#include "gromacs/mdlib/expanded.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdlib/force_flags.h"
#include "gromacs/mdlib/forcerec.h"
#include "gromacs/mdlib/freeenergyparameters.h"
#include "gromacs/mdlib/md_support.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdlib/mdoutf.h"
#include "gromacs/mdlib/membed.h"
#include "gromacs/mdlib/resethandler.h"
#include "gromacs/mdlib/sighandler.h"
#include "gromacs/mdlib/simulationsignal.h"
#include "gromacs/mdlib/stat.h"
#include "gromacs/mdlib/stophandler.h"
#include "gromacs/mdlib/tgroup.h"
#include "gromacs/mdlib/trajectory_writing.h"
#include "gromacs/mdlib/update.h"
#include "gromacs/mdlib/update_constrain_gpu.h"
#include "gromacs/mdlib/update_vv.h"
#include "gromacs/mdlib/vcm.h"
#include "gromacs/mdlib/vsite.h"
#include "gromacs/mdrunutility/handlerestart.h"
#include "gromacs/mdrunutility/multisim.h"
#include "gromacs/mdrunutility/printtime.h"
#include "gromacs/mdtypes/awh_history.h"
#include "gromacs/mdtypes/awh_params.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/df_history.h"
#include "gromacs/mdtypes/energyhistory.h"
#include "gromacs/mdtypes/fcdata.h"
#include "gromacs/mdtypes/forcebuffers.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/mdrunoptions.h"
#include "gromacs/mdtypes/multipletimestepping.h"
#include "gromacs/mdtypes/observableshistory.h"
#include "gromacs/mdtypes/pullhistory.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/mdtypes/state_propagator_data_gpu.h"
#include "gromacs/modularsimulator/energydata.h"
#include "gromacs/nbnxm/gpu_data_mgmt.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/output.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/swap/swapcoords.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/timing/walltime_accounting.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
/* TRANSFER */
#include "gromacs/transfer/transfer.h"
/* END TRANSFER */
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

#include "legacysimulator.h"
#include "replicaexchange.h"
#include "shellfc.h"

/* PLUMED */
#if (GMX_PLUMED)
#include "../../../Plumed.h"
#include "gromacs/math/units.h"
extern int    plumedswitch;
extern plumed plumedmain;
#endif
/* END PLUMED */

using gmx::SimulationSignaller;

void gmx::LegacySimulator::do_md()
{
    // TODO Historically, the EM and MD "integrators" used different
    // names for the t_inputrec *parameter, but these must have the
    // same name, now that it's a member of a struct. We use this ir
    // alias to avoid a large ripple of nearly useless changes.
    // t_inputrec is being replaced by IMdpOptionsProvider, so this
    // will go away eventually.
    t_inputrec*  ir = inputrec;
    int64_t      step, step_rel;
    double       t, t0 = ir->init_t;
    gmx_bool     bGStatEveryStep, bGStat, bCalcVir, bCalcEnerStep, bCalcEner;
    gmx_bool     bNS = FALSE, bNStList, bStopCM, bFirstStep, bInitStep, bLastStep = FALSE;
    gmx_bool     bDoDHDL = FALSE, bDoFEP = FALSE, bDoExpanded = FALSE;
    gmx_bool     do_ene, do_log, do_verbose;
    gmx_bool     bMasterState;
    unsigned int force_flags;
    tensor force_vir = { { 0 } }, shake_vir = { { 0 } }, total_vir = { { 0 } }, pres = { { 0 } };
    int    i, m;
    rvec   mu_tot;
    matrix pressureCouplingMu, M;
    gmx_repl_ex_t     repl_ex = nullptr;
    gmx_global_stat_t gstat;
    gmx_shellfc_t*    shellfc;
    gmx_bool          bSumEkinhOld, bDoReplEx, bExchanged, bNeedRepartition;
    gmx_bool          bTrotter;
    real              dvdl_constr;
    std::vector<RVec> cbuf;
    matrix            lastbox;
    int               lamnew = 0;
    /* for FEP */
    int       nstfep = 0;
    double    cycles;
    real      saved_conserved_quantity = 0;
    real      last_ekin                = 0;
    t_extmass MassQ;
    char      sbuf[STEPSTRSIZE], sbuf2[STEPSTRSIZE];

    /* PME load balancing data for GPU kernels */
    gmx_bool bPMETune         = FALSE;
    gmx_bool bPMETunePrinting = FALSE;

    bool bInteractiveMDstep = false;

    /* PLUMED */
#if (GMX_PLUMED)
    int plumedNeedsEnergy=0;
    int plumedWantsToStop=0;
    matrix plumed_vir;
#endif
    /* END PLUMED */

    /************
     * TRANSFER *
     ************/

    // TODO: maybe this should be ct_broyden_t? or rename the type dftb_broyden_t?
    dftb_broyden_t *ct_broyden=NULL;
    ct_diis_t *ct_diis=NULL;
    FILE *f_tb_hamiltonian=NULL, *f_tb_hamiltonian_ml=NULL, *f_tb_hamiltonian_hub=NULL, *f_tb_hubbard=NULL, *f_tb_occupation=NULL, *f_ct_esp=NULL,
         *f_ct_shift=NULL, *f_ct_energy=NULL, *f_ct_adiabatic=NULL, *f_ct_exp_adiab=NULL, *f_ct_surfacehopping=NULL, *f_ct_nonadiab_coupling=NULL,
         *f_ct_state_vectors=NULL, *f_ct_orbital=NULL, *f_ct_current=NULL, *f_ct_project_wf=NULL, *f_ct_project_wf_ref=NULL, *f_ct_spec=NULL, *f_ct_spec_evec=NULL,
      // *f_ct_overlap_fo=NULL,*f_ct_overlap_fo_ref=NULL, *f_ct_ortho_ev=NULL, *f_ct_fo_ref=NULL,
         *f_ct_tda=NULL, *f_tb_wfprops=NULL, *f_ct_startwf=NULL, *f_ct_active=NULL, *f_tb_com=NULL, *f_tb_spread=NULL, *f_tb_diapop=NULL,
         *f_tb_adiapop=NULL, *f_tb_occupation_adiab=NULL, *f_tb_diag_force=NULL, *f_tb_offdiag_force=NULL;

#if GMX_MPI
    MPI_Comm ct_mpi_comm;
    MPI_Status ct_mpi_status;
    int *ct_mpi_gathered; //, nprocessed;
    if (DOMAINDECOMP(cr)){
        printf("CT-code has to be tested with Domain decomposition\n");
     // printf("CT-code works only with particle decomposition (mdrun -pd)\n");
        exit(-1);
    }
    int ct_mpi_rank=0, ct_mpi_size=0;
#endif

    static struct timespec time_1, time_end;

 // int pair_index[14][6];
 //
 // for (i=0;i<14;i++)
 //     for (j=0;j<6;j++)
 //         pair_index[i][j] = 0;
 //
 // pair_index[0][0] = 1;
 // pair_index[0][1] = 2;
 // pair_index[0][2] = 3;
 // pair_index[0][3] = 10;
 // pair_index[0][4] = 11;
 //
 // pair_index[1][0] = 3;
 // pair_index[1][1] = 4;
 // pair_index[1][2] = 10;
 // pair_index[1][3] = 11;
 //
 // pair_index[2][0] = 3;
 // pair_index[2][1] = 5;
 // pair_index[2][2] = 7;
 // pair_index[2][3] = 8;
 //
 // pair_index[3][0] = 4;
 // pair_index[3][1] = 5;
 // pair_index[3][2] = 6;
 // pair_index[3][3] = 7;
 // pair_index[3][4] = 8;
 // pair_index[3][5] = 9;
 //
 // pair_index[4][0] = 6;
 // pair_index[4][1] = 8;
 // pair_index[4][2] = 9;
 //
 // pair_index[5][0] = 6;
 // pair_index[5][1] = 12;
 // pair_index[5][2] = 13;
 //
 // pair_index[6][0] = 12;
 // pair_index[6][1] = 13;
 //
 // pair_index[7][0] = 8;
 //
 // pair_index[8][0] = 9;
 //
 // pair_index[10][0] = 11;
 //
 // pair_index[12][0] = 13;

    /****************
     * END TRANSFER *
     ****************/

    /* Domain decomposition could incorrectly miss a bonded
       interaction, but checking for that requires a global
       communication stage, which does not otherwise happen in DD
       code. So we do that alongside the first global energy reduction
       after a new DD is made. These variables handle whether the
       check happens, and the result it returns. */
    bool shouldCheckNumberOfBondedInteractions = false;
    int  totalNumberOfBondedInteractions       = -1;

    SimulationSignals signals;
    // Most global communnication stages don't propagate mdrun
    // signals, and will use this object to achieve that.
    SimulationSignaller nullSignaller(nullptr, nullptr, nullptr, false, false);

    if (!mdrunOptions.writeConfout)
    {
        // This is on by default, and the main known use case for
        // turning it off is for convenience in benchmarking, which is
        // something that should not show up in the general user
        // interface.
        GMX_LOG(mdlog.info)
                .asParagraph()
                .appendText(
                        "The -noconfout functionality is deprecated, and may be removed in a "
                        "future version.");
    }

    /* md-vv uses averaged full step velocities for T-control
       md-vv-avek uses averaged half step velocities for T-control (but full step ekin for P control)
       md uses averaged half step kinetic energies to determine temperature unless defined otherwise by GMX_EKIN_AVE_VEL; */
    bTrotter = (EI_VV(ir->eI)
                && (inputrecNptTrotter(ir) || inputrecNphTrotter(ir) || inputrecNvtTrotter(ir)));

    const bool bRerunMD = false;

    int nstglobalcomm = computeGlobalCommunicationPeriod(mdlog, ir, cr);
    bGStatEveryStep   = (nstglobalcomm == 1);

    const SimulationGroups* groups = &top_global->groups;

    std::unique_ptr<EssentialDynamics> ed = nullptr;
    if (opt2bSet("-ei", nfile, fnm))
    {
        /* Initialize essential dynamics sampling */
        ed = init_edsam(mdlog, opt2fn_null("-ei", nfile, fnm), opt2fn("-eo", nfile, fnm), top_global,
                        ir, cr, constr, state_global, observablesHistory, oenv, startingBehavior);
    }
    else if (observablesHistory->edsamHistory)
    {
        gmx_fatal(FARGS,
                  "The checkpoint is from a run with essential dynamics sampling, "
                  "but the current run did not specify the -ei option. "
                  "Either specify the -ei option to mdrun, or do not use this checkpoint file.");
    }

    int*                fep_state = MASTER(cr) ? &state_global->fep_state : nullptr;
    gmx::ArrayRef<real> lambda    = MASTER(cr) ? state_global->lambda : gmx::ArrayRef<real>();
    initialize_lambdas(fplog, *ir, MASTER(cr), fep_state, lambda);
    Update     upd(*ir, deform);
    const bool doSimulatedAnnealing = initSimulatedAnnealing(ir, &upd);
    const bool useReplicaExchange   = (replExParams.exchangeInterval > 0);

    const t_fcdata& fcdata = *fr->fcdata;

    bool simulationsShareState = false;
    int  nstSignalComm         = nstglobalcomm;
    {
        // TODO This implementation of ensemble orientation restraints is nasty because
        // a user can't just do multi-sim with single-sim orientation restraints.
        bool usingEnsembleRestraints =
                (fcdata.disres->nsystems > 1) || ((ms != nullptr) && (fcdata.orires->nr != 0));
        bool awhUsesMultiSim = (ir->bDoAwh && ir->awhParams->shareBiasMultisim && (ms != nullptr));

        // Replica exchange, ensemble restraints and AWH need all
        // simulations to remain synchronized, so they need
        // checkpoints and stop conditions to act on the same step, so
        // the propagation of such signals must take place between
        // simulations, not just within simulations.
        // TODO: Make algorithm initializers set these flags.
        simulationsShareState = useReplicaExchange || usingEnsembleRestraints || awhUsesMultiSim
#if (GMX_PLUMED)
            // PLUMED hack, if we have multiple sim and plumed we usually want them to be in sync 
            || (plumedswitch && ms != nullptr)
#endif
            ;

        if (simulationsShareState)
        {
            // Inter-simulation signal communication does not need to happen
            // often, so we use a minimum of 200 steps to reduce overhead.
            const int c_minimumInterSimulationSignallingInterval = 200;
            nstSignalComm = ((c_minimumInterSimulationSignallingInterval + nstglobalcomm - 1) / nstglobalcomm)
                            * nstglobalcomm;
        }
    }

    if (startingBehavior != StartingBehavior::RestartWithAppending)
    {
        pleaseCiteCouplingAlgorithms(fplog, *ir);
    }
    gmx_mdoutf* outf =
            init_mdoutf(fplog, nfile, fnm, mdrunOptions, cr, outputProvider, mdModulesNotifier, ir,
                        top_global, oenv, wcycle, startingBehavior, simulationsShareState, ms);
    gmx::EnergyOutput energyOutput(mdoutf_get_fp_ene(outf), top_global, ir, pull_work,
                                   mdoutf_get_fp_dhdl(outf), false, startingBehavior,
                                   simulationsShareState, mdModulesNotifier);

    gstat = global_stat_init(ir);

    const auto& simulationWork     = runScheduleWork->simulationWork;
    const bool  useGpuForPme       = simulationWork.useGpuPme;
    const bool  useGpuForNonbonded = simulationWork.useGpuNonbonded;
    const bool  useGpuForBufferOps = simulationWork.useGpuBufferOps;
    const bool  useGpuForUpdate    = simulationWork.useGpuUpdate;

    /* Check for polarizable models and flexible constraints */
    shellfc = init_shell_flexcon(fplog, top_global, constr ? constr->numFlexibleConstraints() : 0,
                                 ir->nstcalcenergy, DOMAINDECOMP(cr), useGpuForPme);

    {
        double io = compute_io(ir, top_global->natoms, *groups, energyOutput.numEnergyTerms(), 1);
        if ((io > 2000) && MASTER(cr))
        {
            fprintf(stderr, "\nWARNING: This run will generate roughly %.0f Mb of data\n\n", io);
        }
    }

    // Local state only becomes valid now.
    std::unique_ptr<t_state> stateInstance;
    t_state*                 state;

    gmx_localtop_t top(top_global->ffparams);

    auto mdatoms = mdAtoms->mdatoms();

    ForceBuffers f(fr->useMts, ((useGpuForNonbonded && useGpuForBufferOps) || useGpuForUpdate)
                                       ? PinningPolicy::PinnedIfSupported
                                       : PinningPolicy::CannotBePinned);
    if (DOMAINDECOMP(cr))
    {
        stateInstance = std::make_unique<t_state>();
        state         = stateInstance.get();
        dd_init_local_state(cr->dd, state_global, state);

        /* Distribute the charge groups over the nodes from the master node */
        dd_partition_system(fplog, mdlog, ir->init_step, cr, TRUE, 1, state_global, *top_global, ir,
                            imdSession, pull_work, state, &f, mdAtoms, &top, fr, vsite, constr,
                            nrnb, nullptr, FALSE);
        shouldCheckNumberOfBondedInteractions = true;
        upd.setNumAtoms(state->natoms);
    }
    else
    {
        state_change_natoms(state_global, state_global->natoms);
        /* Copy the pointer to the global state */
        state = state_global;

        /* Generate and initialize new topology */
        mdAlgorithmsSetupAtomData(cr, ir, *top_global, &top, fr, &f, mdAtoms, constr, vsite, shellfc);

        upd.setNumAtoms(state->natoms);
    }

    std::unique_ptr<UpdateConstrainGpu> integrator;

    StatePropagatorDataGpu* stateGpu = fr->stateGpu;

    // TODO: the assertions below should be handled by UpdateConstraintsBuilder.
    if (useGpuForUpdate)
    {
        GMX_RELEASE_ASSERT(!DOMAINDECOMP(cr) || ddUsesUpdateGroups(*cr->dd) || constr == nullptr
                                   || constr->numConstraintsTotal() == 0,
                           "Constraints in domain decomposition are only supported with update "
                           "groups if using GPU update.\n");
        GMX_RELEASE_ASSERT(ir->eConstrAlg != econtSHAKE || constr == nullptr
                                   || constr->numConstraintsTotal() == 0,
                           "SHAKE is not supported with GPU update.");
        GMX_RELEASE_ASSERT(useGpuForPme || (useGpuForNonbonded && simulationWork.useGpuBufferOps),
                           "Either PME or short-ranged non-bonded interaction tasks must run on "
                           "the GPU to use GPU update.\n");
        GMX_RELEASE_ASSERT(ir->eI == eiMD,
                           "Only the md integrator is supported with the GPU update.\n");
        GMX_RELEASE_ASSERT(
                ir->etc != etcNOSEHOOVER,
                "Nose-Hoover temperature coupling is not supported with the GPU update.\n");
        GMX_RELEASE_ASSERT(
                ir->epc == epcNO || ir->epc == epcPARRINELLORAHMAN || ir->epc == epcBERENDSEN
                        || ir->epc == epcCRESCALE,
                "Only Parrinello-Rahman, Berendsen, and C-rescale pressure coupling are supported "
                "with the GPU update.\n");
        GMX_RELEASE_ASSERT(!mdatoms->haveVsites,
                           "Virtual sites are not supported with the GPU update.\n");
        GMX_RELEASE_ASSERT(ed == nullptr,
                           "Essential dynamics is not supported with the GPU update.\n");
        GMX_RELEASE_ASSERT(!ir->bPull || !pull_have_constraint(*ir->pull),
                           "Constraints pulling is not supported with the GPU update.\n");
        GMX_RELEASE_ASSERT(fcdata.orires->nr == 0,
                           "Orientation restraints are not supported with the GPU update.\n");
        GMX_RELEASE_ASSERT(
                ir->efep == efepNO
                        || (!haveFepPerturbedMasses(*top_global) && !havePerturbedConstraints(*top_global)),
                "Free energy perturbation of masses and constraints are not supported with the GPU "
                "update.");

        if (constr != nullptr && constr->numConstraintsTotal() > 0)
        {
            GMX_LOG(mdlog.info)
                    .asParagraph()
                    .appendText("Updating coordinates and applying constraints on the GPU.");
        }
        else
        {
            GMX_LOG(mdlog.info).asParagraph().appendText("Updating coordinates on the GPU.");
        }
        GMX_RELEASE_ASSERT(fr->deviceStreamManager != nullptr,
                           "Device stream manager should be initialized in order to use GPU "
                           "update-constraints.");
        GMX_RELEASE_ASSERT(
                fr->deviceStreamManager->streamIsValid(gmx::DeviceStreamType::UpdateAndConstraints),
                "Update stream should be initialized in order to use GPU "
                "update-constraints.");
        integrator = std::make_unique<UpdateConstrainGpu>(
                *ir, *top_global, fr->deviceStreamManager->context(),
                fr->deviceStreamManager->stream(gmx::DeviceStreamType::UpdateAndConstraints),
                stateGpu->xUpdatedOnDevice(), wcycle);

        integrator->setPbc(PbcType::Xyz, state->box);
    }

    if (useGpuForPme || (useGpuForNonbonded && useGpuForBufferOps) || useGpuForUpdate)
    {
        changePinningPolicy(&state->x, PinningPolicy::PinnedIfSupported);
    }
    if (useGpuForUpdate)
    {
        changePinningPolicy(&state->v, PinningPolicy::PinnedIfSupported);
    }

    // NOTE: The global state is no longer used at this point.
    // But state_global is still used as temporary storage space for writing
    // the global state to file and potentially for replica exchange.
    // (Global topology should persist.)

    update_mdatoms(mdatoms, state->lambda[efptMASS]);

    if (ir->bExpanded)
    {
        /* Check nstexpanded here, because the grompp check was broken */
        if (ir->expandedvals->nstexpanded % ir->nstcalcenergy != 0)
        {
            gmx_fatal(FARGS,
                      "With expanded ensemble, nstexpanded should be a multiple of nstcalcenergy");
        }
        init_expanded_ensemble(startingBehavior != StartingBehavior::NewSimulation, ir, state->dfhist);
    }

    if (MASTER(cr))
    {
        EnergyData::initializeEnergyHistory(startingBehavior, observablesHistory, &energyOutput);
    }

    preparePrevStepPullCom(ir, pull_work, mdatoms->massT, state, state_global, cr,
                           startingBehavior != StartingBehavior::NewSimulation);

    // TODO: Remove this by converting AWH into a ForceProvider
    auto awh = prepareAwhModule(fplog, *ir, state_global, cr, ms,
                                startingBehavior != StartingBehavior::NewSimulation,
                                shellfc != nullptr, opt2fn("-awh", nfile, fnm), pull_work);

    if (useReplicaExchange && MASTER(cr))
    {
        repl_ex = init_replica_exchange(fplog, ms, top_global->natoms, ir, replExParams);
    }
    /* PME tuning is only supported in the Verlet scheme, with PME for
     * Coulomb. It is not supported with only LJ PME. */
    bPMETune = (mdrunOptions.tunePme && EEL_PME(fr->ic->eeltype) && !mdrunOptions.reproducible
                && ir->cutoff_scheme != ecutsGROUP);

    pme_load_balancing_t* pme_loadbal = nullptr;
    if (bPMETune)
    {
        pme_loadbal_init(&pme_loadbal, cr, mdlog, *ir, state->box, *fr->ic, *fr->nbv, fr->pmedata,
                         fr->nbv->useGpu());
    }

    if (!ir->bContinuation)
    {
        if (state->flags & (1U << estV))
        {
            auto v = makeArrayRef(state->v);
            /* Set the velocities of vsites, shells and frozen atoms to zero */
            for (i = 0; i < mdatoms->homenr; i++)
            {
                if (mdatoms->ptype[i] == eptVSite || mdatoms->ptype[i] == eptShell)
                {
                    clear_rvec(v[i]);
                }
                else if (mdatoms->cFREEZE)
                {
                    for (m = 0; m < DIM; m++)
                    {
                        if (ir->opts.nFreeze[mdatoms->cFREEZE[i]][m])
                        {
                            v[i][m] = 0;
                        }
                    }
                }
            }
        }

        if (constr)
        {
            /* Constrain the initial coordinates and velocities */
            do_constrain_first(fplog, constr, ir, mdatoms->nr, mdatoms->homenr,
                               state->x.arrayRefWithPadding(), state->v.arrayRefWithPadding(),
                               state->box, state->lambda[efptBONDED]);
        }
        if (vsite)
        {
            /* Construct the virtual sites for the initial configuration */
            vsite->construct(state->x, ir->delta_t, {}, state->box);
        }
    }

    if (ir->efep != efepNO)
    {
        /* Set free energy calculation frequency as the greatest common
         * denominator of nstdhdl and repl_ex_nst. */
        nstfep = ir->fepvals->nstdhdl;
        if (ir->bExpanded)
        {
            nstfep = std::gcd(ir->expandedvals->nstexpanded, nstfep);
        }
        if (useReplicaExchange)
        {
            nstfep = std::gcd(replExParams.exchangeInterval, nstfep);
        }
        if (ir->bDoAwh)
        {
            nstfep = std::gcd(ir->awhParams->nstSampleCoord, nstfep);
        }
    }

    /* Be REALLY careful about what flags you set here. You CANNOT assume
     * this is the first step, since we might be restarting from a checkpoint,
     * and in that case we should not do any modifications to the state.
     */
    bStopCM = (ir->comm_mode != ecmNO && !ir->bContinuation);

    // When restarting from a checkpoint, it can be appropriate to
    // initialize ekind from quantities in the checkpoint. Otherwise,
    // compute_globals must initialize ekind before the simulation
    // starts/restarts. However, only the master rank knows what was
    // found in the checkpoint file, so we have to communicate in
    // order to coordinate the restart.
    //
    // TODO Consider removing this communication if/when checkpoint
    // reading directly follows .tpr reading, because all ranks can
    // agree on hasReadEkinState at that time.
    bool hasReadEkinState = MASTER(cr) ? state_global->ekinstate.hasReadEkinState : false;
    if (PAR(cr))
    {
        gmx_bcast(sizeof(hasReadEkinState), &hasReadEkinState, cr->mpi_comm_mygroup);
    }
    if (hasReadEkinState)
    {
        restore_ekinstate_from_state(cr, ekind, &state_global->ekinstate);
    }

    unsigned int cglo_flags =
            (CGLO_TEMPERATURE | CGLO_GSTAT | (EI_VV(ir->eI) ? CGLO_PRESSURE : 0)
             | (EI_VV(ir->eI) ? CGLO_CONSTRAINT : 0) | (hasReadEkinState ? CGLO_READEKIN : 0));

    bSumEkinhOld = FALSE;

    t_vcm vcm(top_global->groups, *ir);
    reportComRemovalInfo(fplog, vcm);

    /* To minimize communication, compute_globals computes the COM velocity
     * and the kinetic energy for the velocities without COM motion removed.
     * Thus to get the kinetic energy without the COM contribution, we need
     * to call compute_globals twice.
     */
    for (int cgloIteration = 0; cgloIteration < (bStopCM ? 2 : 1); cgloIteration++)
    {
        unsigned int cglo_flags_iteration = cglo_flags;
        if (bStopCM && cgloIteration == 0)
        {
            cglo_flags_iteration |= CGLO_STOPCM;
            cglo_flags_iteration &= ~CGLO_TEMPERATURE;
        }
        compute_globals(gstat, cr, ir, fr, ekind, makeConstArrayRef(state->x),
                        makeConstArrayRef(state->v), state->box, mdatoms, nrnb, &vcm, nullptr,
                        enerd, force_vir, shake_vir, total_vir, pres, constr, &nullSignaller,
                        state->box, &totalNumberOfBondedInteractions, &bSumEkinhOld,
                        cglo_flags_iteration
                                | (shouldCheckNumberOfBondedInteractions ? CGLO_CHECK_NUMBER_OF_BONDED_INTERACTIONS
                                                                         : 0));
        if (cglo_flags_iteration & CGLO_STOPCM)
        {
            /* At initialization, do not pass x with acceleration-correction mode
             * to avoid (incorrect) correction of the initial coordinates.
             */
            auto x = (vcm.mode == ecmLINEAR_ACCELERATION_CORRECTION) ? ArrayRef<RVec>()
                                                                     : makeArrayRef(state->x);
            process_and_stopcm_grp(fplog, &vcm, *mdatoms, x, makeArrayRef(state->v));
            inc_nrnb(nrnb, eNR_STOPCM, mdatoms->homenr);
        }
    }
    checkNumberOfBondedInteractions(mdlog, cr, totalNumberOfBondedInteractions, top_global, &top,
                                    makeConstArrayRef(state->x), state->box,
                                    &shouldCheckNumberOfBondedInteractions);
    if (ir->eI == eiVVAK)
    {
        /* a second call to get the half step temperature initialized as well */
        /* we do the same call as above, but turn the pressure off -- internally to
           compute_globals, this is recognized as a velocity verlet half-step
           kinetic energy calculation.  This minimized excess variables, but
           perhaps loses some logic?*/

        compute_globals(gstat, cr, ir, fr, ekind, makeConstArrayRef(state->x),
                        makeConstArrayRef(state->v), state->box, mdatoms, nrnb, &vcm, nullptr,
                        enerd, force_vir, shake_vir, total_vir, pres, constr, &nullSignaller,
                        state->box, nullptr, &bSumEkinhOld, cglo_flags & ~CGLO_PRESSURE);
    }

    /* Calculate the initial half step temperature, and save the ekinh_old */
    if (startingBehavior == StartingBehavior::NewSimulation)
    {
        for (i = 0; (i < ir->opts.ngtc); i++)
        {
            copy_mat(ekind->tcstat[i].ekinh, ekind->tcstat[i].ekinh_old);
        }
    }

    /* need to make an initiation call to get the Trotter variables set, as well as other constants
       for non-trotter temperature control */
    auto trotter_seq = init_npt_vars(ir, state, &MassQ, bTrotter);

    if (MASTER(cr))
    {
        if (!ir->bContinuation)
        {
            if (constr && ir->eConstrAlg == econtLINCS)
            {
                fprintf(fplog, "RMS relative constraint deviation after constraining: %.2e\n",
                        constr->rmsd());
            }
            if (EI_STATE_VELOCITY(ir->eI))
            {
                real temp = enerd->term[F_TEMP];
                if (ir->eI != eiVV)
                {
                    /* Result of Ekin averaged over velocities of -half
                     * and +half step, while we only have -half step here.
                     */
                    temp *= 2;
                }
                fprintf(fplog, "Initial temperature: %g K\n", temp);
            }
        }

        char tbuf[20];
        fprintf(stderr, "starting mdrun '%s'\n", *(top_global->name));
        if (ir->nsteps >= 0)
        {
            sprintf(tbuf, "%8.1f", (ir->init_step + ir->nsteps) * ir->delta_t);
        }
        else
        {
            sprintf(tbuf, "%s", "infinite");
        }
        if (ir->init_step > 0)
        {
            fprintf(stderr, "%s steps, %s ps (continuing from step %s, %8.1f ps).\n",
                    gmx_step_str(ir->init_step + ir->nsteps, sbuf), tbuf,
                    gmx_step_str(ir->init_step, sbuf2), ir->init_step * ir->delta_t);
        }
        else
        {
            fprintf(stderr, "%s steps, %s ps.\n", gmx_step_str(ir->nsteps, sbuf), tbuf);
        }
        fprintf(fplog, "\n");
    }

    /* PLUMED */
#if (GMX_PLUMED)
    if(plumedswitch){
      /* detect plumed API version */
      int pversion=0;
      plumed_cmd(plumedmain,"getApiVersion",&pversion);
      /* setting kbT is only implemented with api>1) */
      real kbT=ir->opts.ref_t[0]*BOLTZ;
      if(pversion>1) plumed_cmd(plumedmain,"setKbT",&kbT);
      if(pversion>2){
        int res=1;
        if( (startingBehavior != StartingBehavior::NewSimulation) ) plumed_cmd(plumedmain,"setRestart",&res);
      }

    //if(cr->ms && cr->ms->nsim>1) {
      if(ms && ms->numSimulations_>1) {
        if(MASTER(cr)) plumed_cmd(plumedmain,"GREX setMPIIntercomm",&ms->mastersComm_);
        if(PAR(cr)){
          if(DOMAINDECOMP(cr)) {
            plumed_cmd(plumedmain,"GREX setMPIIntracomm",&cr->dd->mpi_comm_all);
          }else{
            plumed_cmd(plumedmain,"GREX setMPIIntracomm",&cr->mpi_comm_mysim);
          }
        }
        plumed_cmd(plumedmain,"GREX init",NULL);
      }
      if(PAR(cr)){
        if(DOMAINDECOMP(cr)) {
          plumed_cmd(plumedmain,"setMPIComm",&cr->dd->mpi_comm_all);
        }
      }
      plumed_cmd(plumedmain,"setNatoms",&top_global->natoms);
      plumed_cmd(plumedmain,"setMDEngine","gromacs");
      plumed_cmd(plumedmain,"setLog",fplog);
      real real_delta_t=ir->delta_t;
      plumed_cmd(plumedmain,"setTimestep",&real_delta_t);
      plumed_cmd(plumedmain,"init",NULL);

      if(PAR(cr)){
        if(DOMAINDECOMP(cr)) {
          int nat_home = dd_numHomeAtoms(*cr->dd);
          plumed_cmd(plumedmain,"setAtomsNlocal",&nat_home);
          plumed_cmd(plumedmain,"setAtomsGatindex",cr->dd->globalAtomIndices.data());
        }
      }
    }
#endif
    /* END PLUMED */

    /************
     * TRANSFER *
     ************/

    /* initialize the simulation of transfer */

#if GMX_MPI
    ct_mpi_comm = MPI_COMM_WORLD;
    ct_mpi_size = cr->nnodes;
    ct_mpi_rank = cr->nodeid;
    if (ct_mpi_size>1 && ct_mpi_rank==0)
    {
        snew(ct_mpi_gathered, ct_mpi_size);
    }
    printf("MPI Initialization on node %d of %d done\n", ct_mpi_rank+1, ct_mpi_size);
#endif

    t_atoms *ct_atoms; //, *ct_atoms2, ct_atoms3;
    snew(ct_atoms, 1);

#if GMX_MPI
    if (ct_mpi_rank == 0)
    {
        *ct_atoms = gmx_mtop_global_atoms(top_global); // expands the topology. now arrays in struct t_atoms contain all atoms.
                                                     // in top_global they are grouped together in moltypes. only possible on rank0.
                                                     // relevant data have to be sent to all other nodes.
        snew(ct_atoms->pdbinfo, ct_atoms->nr); // pdb files can be handy for special output where e.g. charge is written in bfactor
                                             // which can be visualized with VMD

        if (PROTEIN)
        {
            printf("start protein preprocessing\n");
         // snew(ct_atoms2,1);
         // init_t_atoms(ct_atoms2, ct_atoms->nr, FALSE);
         // ct_atoms = protein_preprocessor(ct_atoms, state_global); //changes residues and atomnames for proteins
            printf("protein preprocessing done\n");
            for (j=0; j<ct_atoms->nr; j++)
            {
             // printf("CHECK %s %s\n",*(ct_atoms->atomname[j]), *(ct_atoms->atomname[j]));
             // printf("CHECK %d %d\n",ct_atoms->atom[j].resind, ct_atoms->atom[j].resind );
             // printf("CHECK %s %s\n",*(ct_atoms->resinfo[ct_atoms->atom[j].resind].name), *(ct_atoms->resinfo[ct_atoms->atom[j].resind].name) );
            }
/*
            FILE *f_ct_preprocessor=NULL;
            f_ct_preprocessor = fopen("preprocessor.gro", "w");
            fprintf(f_ct_preprocessor, "gro file with new residues\n %d\n", ct_atoms->nr);
            for (j=0; j<ct_atoms->nr; j++)
                fprintf(f_ct_preprocessor, "%5d%-5.5s%5.5s%5d %lf %lf %lf\n", ct_atoms->atom[j].resind%100000, "NEW",
                        *(ct_atoms->atomname[j]), (j+1)%100000, state_global->x[j][XX], state_global->x[j][YY], state_global->x[j][ZZ]);
            fclose(f_ct_preprocessor);
*/         
         // write_sto_conf("preprocessor.gro", ".gro file to check residues", ct_atoms,  state_global->x, NULL, -1, NULL);
        }

        for (i=1; i < ct_mpi_size; i++)
        {
            MPI_Send(&ct_atoms->nres, 1, MPI_INT, i, 200+i, ct_mpi_comm);
            MPI_Send(&ct_atoms->nr, 1, MPI_INT, i, 300+i, ct_mpi_comm);
            for (j=0; j < ct_atoms->nr; j++)
            {
                MPI_Send(*ct_atoms->atomname[j],  6 , MPI_CHAR, i, (i)*1000000+j, ct_mpi_comm);
            }
            MPI_Send(ct_atoms->atom, ct_atoms->nr * sizeof(ct_atoms->atom[0]), MPI_CHAR, i, 400+i, ct_mpi_comm);
            MPI_Send(ct_atoms->resinfo, ct_atoms->nres * sizeof(ct_atoms->resinfo[0]), MPI_CHAR, i, 500+i, ct_mpi_comm);
        }
    }
    else
    {
        MPI_Recv(&ct_atoms->nr, 1, MPI_INT, 0, 300+ct_mpi_rank, ct_mpi_comm, &ct_mpi_status);
        printf("nr %d\n",ct_atoms->nr);
        snew(ct_atoms->atom, ct_atoms->nr);
        snew(ct_atoms->atomname, ct_atoms->nr);
        for (j=0; j < ct_atoms->nr; j++)
        {
            snew(ct_atoms->atomname[j], 1);
            snew(ct_atoms->atomname[j][0], 6);
        }
        MPI_Recv(&ct_atoms->nres, 1, MPI_INT, 0, 200+ct_mpi_rank, ct_mpi_comm, &ct_mpi_status);
     // printf("nres %d\n",ct_atoms->nres);
        snew(ct_atoms->resinfo, ct_atoms->nres);

        for (j=0; j < ct_atoms->nr; j++)
        {
            MPI_Recv(*ct_atoms->atomname[j],  6 , MPI_CHAR, 0, (ct_mpi_rank)*1000000+j, ct_mpi_comm, &ct_mpi_status);
         // printf("%5s\n",*(ct_atoms->atomname[j]));
        }
        MPI_Recv(ct_atoms->atom, ct_atoms->nr * sizeof(ct_atoms->atom[0]), MPI_CHAR, 0, 400+ct_mpi_rank, ct_mpi_comm, &ct_mpi_status);
        MPI_Recv(ct_atoms->resinfo, ct_atoms->nres * sizeof(ct_atoms->resinfo[0]), MPI_CHAR, 0, 500+ct_mpi_rank , ct_mpi_comm, &ct_mpi_status);
    }
#else
    *ct_atoms = gmx_mtop_global_atoms(top_global);
    if (PROTEIN) {
        *ct_atoms = gmx_mtop_global_atoms(top_global);
        ct_atoms = protein_preprocessor(ct_atoms, state_global); //changes residues for proteins
    }
#endif

    // TRANSFER data structure

    charge_transfer_t *ct=NULL;
    snew(ct, 1);

#if GMX_MPI
    init_charge_transfer(ct_atoms, top_global, mdatoms, ct, state, ct_mpi_rank);
#else
    init_charge_transfer(ct_atoms, top_global, mdatoms, ct, state);
#endif
    ct->rk_timestep = ir->delta_t * PS_TO_AU;
    printf("Time step for quantum mechanics: %f ps = %f a.u.\n", ct->rk_timestep / PS_TO_AU, ct->rk_timestep);

    // data structures for QM/MM and DFTB+
    int qmmmVariant;
    bool dipCorrection = false;
    switch (ct->qmmm) {
        case 0:  qmmmVariant = 0; break; // vacuo, no QM/MM
        case 1:  qmmmVariant = 2; break; // cut-off, apply the switch version of Gromacs
        case 2:  qmmmVariant = 2; break; // dtto, however only considering MM atoms from a specified list
                                         // TODO: maybe what was intended is the *complete* list?
                                         // that might require an additional value of qmmmVariant - but it would be easy
        case 3:  qmmmVariant = 1; break; // PME, currently always without dipole correction
        default: printf("Impossible value of QM/MM variant (%d), exiting.\n", ct->qmmm); exit(-1);
    }

    std::vector<int> qmAtoms, atomicnumberQM;
    std::vector<std::unique_ptr<QMMM_rec_transfer> > dftbplus_phase1;
    dftbplus_phase1.resize(ct->sites);
    for (int i = 0; i < ct->sites; i++)
    {
        // figure out the atomic numbers from atomtypes (that are stored in the CT structure);
        qmAtoms.clear();
        atomicnumberQM.clear();
        double totMass = 0.;
        for (int j = 0; j<ct->site[i].atoms; j++)
        {
            qmAtoms.push_back(ct->site[i].atom[j]);
         // atomicnumberQM.push_back(atomtype_to_atomicnumber(ct->site[i].atomtype[j]));
            int atomicnumber;
            double atomicmass;
            atomtype_to_atomicnumber(ct->site[i].atomtype[j], atomicnumber, atomicmass);
            atomicnumberQM.push_back(atomicnumber);
            ct->site[i].mass[j] = atomicmass;
            totMass += atomicmass;
        }
        dftbplus_phase1[i] = std::make_unique<QMMM_rec_transfer>(cr, ir, fr, state_global->natoms, qmAtoms, atomicnumberQM, qmmmVariant, dipCorrection);
        ct->site[i].inv_tot_mass = 1. / totMass;
    }
    qmAtoms.clear();
    atomicnumberQM.clear();
    double totMassCplx;
    for (int j = 0; j<ct->atoms_cplx; j++)
    {
        qmAtoms.push_back(ct->atom_cplx[j]);
     // atomicnumberQM.push_back(atomtype_to_atomicnumber(ct->atomtype_cplx[j]));
        int atomicnumber;
        double atomicmass;
        atomtype_to_atomicnumber(ct->atomtype_cplx[j], atomicnumber, atomicmass);
        atomicnumberQM.push_back(atomicnumber);
        ct->mass_cplx[j] = atomicmass;
        totMassCplx += atomicmass;
    }
    std::unique_ptr<QMMM_rec_transfer> dftbplus_phase2 =
                          std::make_unique<QMMM_rec_transfer>(cr, ir, fr, state_global->natoms, qmAtoms, atomicnumberQM, qmmmVariant, dipCorrection);
    ct->inv_tot_mass_cplx = 1. / totMassCplx;
#if GMX_MPI
    init_qmmmrec_transfer(dftbplus_phase1.data(), dftbplus_phase2, ct, cr, ir->rcoulomb, ir->ewald_rtol, ct_mpi_rank);
#else
    init_qmmmrec_transfer(dftbplus_phase1.data(), dftbplus_phase2, ct, cr, ir->rcoulomb, ir->ewald_rtol);
#endif

    dftb_t *dftb;
    snew(dftb, 1);
    init_dftb_stub(dftb, ct); // allocate the neighborlist array

    // end QM/MM and DFTB+

    rvec *x_ct;
    snew(x_ct, state_global->natoms);

    if (ct->jobtype == cteBORNOPPENHEIMER) {
      double original_delta_t = ir->delta_t;
      (void) original_delta_t;
    }
/*
    if (ct->qmmm == 3) // this means PME should be used
#if GMX_MPI
        init_dftb_pme(dftb, ct, ir, ct_mpi_rank);
#else
        init_dftb_pme(dftb, ct, ir);
#endif
*/

#if GMX_MPI
    if (ct_mpi_rank == 0) {
#endif

    if (ct->jobtype != cteESP || ct->force_output == 1)
    {
        if (ct->jobtype==cteNOMOVEMENT || ct->force_output == 1)
        {
            f_tb_com = fopen("TB_COM.xvg" , "w"); /* dholub */
            
            if (ct->hamiltonian_type)
            {
                f_tb_hamiltonian_ml = fopen("TB_HAMILTONIAN_ML.xvg", "w");
            }
            else
            {
                f_tb_hamiltonian = fopen("TB_HAMILTONIAN.xvg", "w");
            }
            // MK tag printout
            f_tb_diag_force = fopen("TB_DIAG_FORCE.xvg", "w");
            f_tb_offdiag_force = fopen("TB_OFFDIAG_FORCE.xvg", "w");
        }

     // f_ct_ortho_ev = fopen("CT_ORTHO_EV.xvg", "w");
     // f_ct_overlap_fo = fopen("CT_OVERLAP_FO.xvg", "w");
     // f_ct_overlap_fo_ref = fopen("CT_OVERLAP_FO_REF.xvg", "w");

        if (ct->jobtype==ctePARAMETERS)
        {
            f_ct_spec = fopen("CT_SPEC.xvg", "w");
            f_ct_spec_evec = fopen("CT_SPEC_EVEC.xvg", "w");
        }
        if (ct->do_projection)
        {
            f_ct_project_wf = fopen("CT_PROJECT_WF.xvg", "w");
            f_ct_project_wf_ref = fopen("CT_PROJECT_WF_REF.xvg", "w");
        }
        if (ct->jobtype==cteTDA)
            f_ct_tda = fopen("CT_TDA.xvg", "w");

        }
        if (ct->qmmm > 0)
            f_ct_esp = fopen("CT_ESP.xvg", "w");

        if (ct->jobtype==cteSCCDYNAMIC || ct->jobtype==cteCPFSCCDYNAMIC || ct->jobtype==cteDDBSCCDYNAMIC || ct->jobtype==cteNONSCCDYNAMIC ||
            ct->jobtype==cteADIABATIC || ct->jobtype==cteADNONSCC || ct->jobtype==cteFERMIADIABATIC || ct->jobtype == cteBORNOPPENHEIMER ||
            ct->jobtype==cteDLZSH || ct->jobtype==cteALZSH || ct->jobtype==cteDFSSH || ct->jobtype==cteJFSSH || ct->jobtype==cteBCJFSSH ||
            ct->jobtype==cteSCRDFSSH || ct->jobtype==cteCCFSSH || ct->jobtype==cteGFSH || ct->jobtype==cteDISH  || ct->jobtype==cteSCRDFSSH ||
            ct->jobtype==cteTULLYLOC || ct->jobtype==ctePERSICOSFHOPPING || ct->jobtype==cteNEGFLORENTZ || ct->jobtype==cteNEGFLORENTZNONSCC ||
            ct->jobtype==ctePREZHDOSFHOPPING)
        {
            f_tb_occupation = fopen("TB_OCCUPATION.xvg", "w");
            f_tb_wfprops = fopen("TB_WFPROPS.xvg", "w");
            f_tb_spread = fopen("TB_SPREAD.xvg" , "w");
        }

        if (ct->jobtype==cteJFSSH || ct->jobtype==cteBCJFSSH || ct->jobtype==cteSCRDFSSH || ct->jobtype==cteCCFSSH || ct->jobtype==cteGFSH ||
            ct->jobtype==cteDISH  || ct->jobtype==cteSCRDFSSH || ct->jobtype==cteTULLYLOC)
        {
            f_tb_occupation_adiab = fopen("TB_OCCUPATION_ADIAB.xvg", "w");
            f_tb_diapop = fopen("TB_DIABPOP.xvg", "w");
            f_tb_adiapop = fopen("TB_ADIABPOP.xvg", "w");
        }

        if (ct->jobtype==cteSCCDYNAMIC || ct->jobtype==cteCPFSCCDYNAMIC || ct->jobtype==cteDDBSCCDYNAMIC || ct->jobtype==cteNONSCCDYNAMIC ||
            ct->jobtype==cteADIABATIC || ct->jobtype==cteADNONSCC || ct->jobtype==cteNOMOVEMENT || ct->jobtype==cteFERMIADIABATIC ||
            ct->jobtype == cteBORNOPPENHEIMER  || ct->jobtype==cteDLZSH || ct->jobtype==cteALZSH || ct->jobtype==cteDFSSH || ct->jobtype==cteJFSSH ||
            ct->jobtype==cteBCJFSSH || ct->jobtype==cteSCRDFSSH || ct->jobtype==cteCCFSSH || ct->jobtype==cteGFSH || ct->jobtype==cteDISH ||
            ct->jobtype==cteSCRDFSSH ||  ct->jobtype==cteTULLYLOC || ct->jobtype==ctePERSICOSFHOPPING)
        {
            f_ct_energy = fopen("CT_ENERGY.xvg", "w");
        }

     // if (ct->jobtype==cteSCCDYNAMIC || ct->jobtype==cteCPFSCCDYNAMIC || ct->jobtype==cteDDBSCCDYNAMIC || ct->jobtype==cteADIABATIC ||
     //     ct->jobtype==cteNOMOVEMENT || ct->jobtype==cteFERMIADIABATIC || ct->jobtype == cteBORNOPPENHEIMER  || ct->jobtype==cteDLZSH ||
     //     ct->jobtype==cteALZSH || ct->jobtype==cteDFSSH || ct->jobtype==cteJFSSH || ct->jobtype==cteBCJFSSH || ct->jobtype==cteSCRDFSSH ||
     //     ct->jobtype==cteCCFSSH || ct->jobtype==cteGFSH || ct->jobtype==cteDISH  || ct->jobtype==cteSCRDFSSH ||  ct->jobtype==cteTULLYLOC ||
     //     ct->jobtype==ctePERSICOSFHOPPING || ct->jobtype==cteNEGFLORENTZ)
     // {
     //     f_tb_hamiltonian_hub = fopen("TB_HAMILTONIAN_HUB.xvg", "w");
     //     f_tb_hubbard = fopen("TB_HUBBARD.xvg", "w");
     // }

        if (ct->jobtype==cteSCCDYNAMIC || ct->jobtype==cteCPFSCCDYNAMIC || ct->jobtype==cteDDBSCCDYNAMIC || ct->jobtype==cteNONSCCDYNAMIC ||
            ct->jobtype==cteADIABATIC || ct->jobtype==cteADNONSCC )
        {
            snew(ct_diis, 1);
            ct_init_diis(ct, ct_diis);
            snew(ct->q_act, ct->dim);
            snew(ct->q_old, ct->dim);
        }

        if (ct->jobtype==cteSCCDYNAMIC || ct->jobtype==cteCPFSCCDYNAMIC || ct->jobtype==cteDDBSCCDYNAMIC || ct->jobtype==cteNONSCCDYNAMIC ||
            ct->jobtype==cteFERMIADIABATIC || ct->jobtype == cteBORNOPPENHEIMER || ct->jobtype==cteDLZSH || ct->jobtype==cteALZSH ||
            ct->jobtype==cteDFSSH || ct->jobtype==cteJFSSH || ct->jobtype==cteBCJFSSH || ct->jobtype==cteSCRDFSSH || ct->jobtype==cteCCFSSH ||
            ct->jobtype==cteGFSH || ct->jobtype==cteDISH  || ct->jobtype==cteSCRDFSSH || ct->jobtype==cteTULLYLOC || ct->jobtype==ctePERSICOSFHOPPING ||
            ct->jobtype==cteNOMOVEMENT || ct->jobtype==ctePREZHDOSFHOPPING) // we need here also NOMOVEMENT if no correct wave function is specified
                                                                            // and one has to take the lowest eigenvalue instead
        { 
            snew(ct_broyden, 1);
            ct_init_broyden(ct, ct_broyden);
            snew(ct->q_act, ct->dim);
            snew(ct->q_old, ct->dim);
        }

        if (ct->jobtype==cteADIABATIC || ct->jobtype==cteFERMIADIABATIC || ct->jobtype == cteBORNOPPENHEIMER) {
            f_ct_adiabatic = fopen("CT_ADIABATIC.xvg", "w");
        }

     // if (ct->jobtype==cteSCCDYNAMIC || ct->jobtype==cteCPFSCCDYNAMIC || ct->jobtype==cteDDBSCCDYNAMIC || ct->jobtype==cteNONSCCDYNAMIC) {
     //     f_ct_exp_adiab = fopen("CT_EXPANS_IN_ADIAB_STATES.xvg", "w");
     // }

        if (ct->jobtype==cteDLZSH || ct->jobtype==cteALZSH || ct->jobtype==cteDFSSH || ct->jobtype==cteJFSSH || ct->jobtype==cteBCJFSSH ||
            ct->jobtype==cteSCRDFSSH || ct->jobtype==cteCCFSSH || ct->jobtype==cteGFSH || ct->jobtype==cteDISH  || ct->jobtype==cteSCRDFSSH ||
            ct->jobtype==cteTULLYLOC || ct->jobtype==ctePERSICOSFHOPPING || ct->jobtype==ctePREZHDOSFHOPPING)
        {
            f_ct_surfacehopping = fopen("CT_SURFACE_HOPPING.xvg", "w");
        }

     // if (ct->jobtype==cteDLZSH || ct->jobtype==cteALZSH || ct->jobtype==cteDFSSH || ct->jobtype==cteJFSSH || ct->jobtype==cteBCJFSSH ||
     //     ct->jobtype==cteSCRDFSSH || ct->jobtype==cteCCFSSH || ct->jobtype==cteGFSH || ct->jobtype==cteDISH  || ct->jobtype==cteSCRDFSSH ||
     //     ct->jobtype==cteTULLYLOC  || ct->jobtype==ctePERSICOSFHOPPING || ct->jobtype==ctePREZHDOSFHOPPING)
     // {
     //     f_ct_nonadiab_coupling = fopen("CT_NONADIAB_COUPLING.xvg", "w");
     //     f_ct_state_vectors = fopen("CT_STATE_VECTORS.xvg", "w");
     // }

     // if (ct->jobtype==ctePARAMETERS || ct->jobtype==cteNOMOVEMENT) {
     //     f_ct_orbital = fopen("CT_ORBITAL.xvg", "w");
     // }

        if (ct->jobtype==cteNEGFLORENTZ || ct->jobtype==cteNEGFLORENTZNONSCC) {
          f_ct_current = fopen("CT_CURRENT.xvg", "w");
        }

        if (ct->pool_size > ct->sites){
          f_ct_active = fopen("CT_ACTIVE.xvg", "w");
        }
#if GMX_MPI
    }
#endif

    /****************
     * END TRANSFER *
     ****************/

    walltime_accounting_start_time(walltime_accounting);
    wallcycle_start(wcycle, ewcRUN);
    print_start(fplog, cr, walltime_accounting, "mdrun");

    /***********************************************************
     *
     *             Loop over MD steps
     *
     ************************************************************/

    bFirstStep = TRUE;
    /* Skip the first Nose-Hoover integration when we get the state from tpx */
    bInitStep        = startingBehavior == StartingBehavior::NewSimulation || EI_VV(ir->eI);
    bSumEkinhOld     = FALSE;
    bExchanged       = FALSE;
    bNeedRepartition = FALSE;

    step     = ir->init_step;
    step_rel = 0;

    auto stopHandler = stopHandlerBuilder->getStopHandlerMD(
            compat::not_null<SimulationSignal*>(&signals[eglsSTOPCOND]), simulationsShareState,
            MASTER(cr), ir->nstlist, mdrunOptions.reproducible, nstSignalComm,
            mdrunOptions.maximumHoursToRun, ir->nstlist == 0, fplog, step, bNS, walltime_accounting);

    auto checkpointHandler = std::make_unique<CheckpointHandler>(
            compat::make_not_null<SimulationSignal*>(&signals[eglsCHKPT]), simulationsShareState,
            ir->nstlist == 0, MASTER(cr), mdrunOptions.writeConfout,
            mdrunOptions.checkpointOptions.period);

    const bool resetCountersIsLocal = true;
    auto       resetHandler         = std::make_unique<ResetHandler>(
            compat::make_not_null<SimulationSignal*>(&signals[eglsRESETCOUNTERS]),
            !resetCountersIsLocal, ir->nsteps, MASTER(cr), mdrunOptions.timingOptions.resetHalfway,
            mdrunOptions.maximumHoursToRun, mdlog, wcycle, walltime_accounting);

    const DDBalanceRegionHandler ddBalanceRegionHandler(cr);

    if (MASTER(cr) && isMultiSim(ms) && !useReplicaExchange)
    {
        logInitialMultisimStatus(ms, cr, mdlog, simulationsShareState, ir->nsteps, ir->init_step);
    }

    /* and stop now if we should */
    bLastStep = (bLastStep || (ir->nsteps >= 0 && step_rel > ir->nsteps));
    while (!bLastStep)
    {

        /* Determine if this is a neighbor search step */
        bNStList = (ir->nstlist > 0 && step % ir->nstlist == 0);

        if (bPMETune && bNStList)
        {
            // This has to be here because PME load balancing is called so early.
            // TODO: Move to after all booleans are defined.
            if (useGpuForUpdate && !bFirstStep)
            {
                stateGpu->copyCoordinatesFromGpu(state->x, AtomLocality::Local);
                stateGpu->waitCoordinatesReadyOnHost(AtomLocality::Local);
            }
            /* PME grid + cut-off optimization with GPUs or PME nodes */
            pme_loadbal_do(pme_loadbal, cr, (mdrunOptions.verbose && MASTER(cr)) ? stderr : nullptr,
                           fplog, mdlog, *ir, fr, state->box, state->x, wcycle, step, step_rel,
                           &bPMETunePrinting, simulationWork.useGpuPmePpCommunication);
        }

        wallcycle_start(wcycle, ewcSTEP);

        bLastStep = (step_rel == ir->nsteps);
        t         = t0 + step * ir->delta_t;
        /* TRANSFER */
        // an older version had this:
     // if (ct->jobtype == cteBORNOPPENHEIMER)
     //     t += ir->delta_t;
     // else
     //     t = t0 + step * ir->delta_t;
        printf("\n %ld %ld\n", step_rel, ir->nsteps);
        /* END TRANSFER */

        // TODO Refactor this, so that nstfep does not need a default value of zero
        if (ir->efep != efepNO || ir->bSimTemp)
        {
            /* find and set the current lambdas */
            state->lambda = currentLambdas(step, *(ir->fepvals), state->fep_state);

            bDoDHDL     = do_per_step(step, ir->fepvals->nstdhdl);
            bDoFEP      = ((ir->efep != efepNO) && do_per_step(step, nstfep));
            bDoExpanded = (do_per_step(step, ir->expandedvals->nstexpanded) && (ir->bExpanded)
                           && (!bFirstStep));
        }

        bDoReplEx = (useReplicaExchange && (step > 0) && !bLastStep
                     && do_per_step(step, replExParams.exchangeInterval));

        if (doSimulatedAnnealing)
        {
            update_annealing_target_temp(ir, t, &upd);
        }

        /* Stop Center of Mass motion */
        bStopCM = (ir->comm_mode != ecmNO && do_per_step(step, ir->nstcomm));

        /* Determine whether or not to do Neighbour Searching */
        bNS = (bFirstStep || bNStList || bExchanged || bNeedRepartition);

        /* Note that the stopHandler will cause termination at nstglobalcomm
         * steps. Since this concides with nstcalcenergy, nsttcouple and/or
         * nstpcouple steps, we have computed the half-step kinetic energy
         * of the previous step and can always output energies at the last step.
         */
        bLastStep = bLastStep || stopHandler->stoppingAfterCurrentStep(bNS);

        /* do_log triggers energy and virial calculation. Because this leads
         * to different code paths, forces can be different. Thus for exact
         * continuation we should avoid extra log output.
         * Note that the || bLastStep can result in non-exact continuation
         * beyond the last step. But we don't consider that to be an issue.
         */
        do_log     = (do_per_step(step, ir->nstlog)
                  || (bFirstStep && startingBehavior == StartingBehavior::NewSimulation) || bLastStep);
        do_verbose = mdrunOptions.verbose
                     && (step % mdrunOptions.verboseStepPrintInterval == 0 || bFirstStep || bLastStep);

        if (useGpuForUpdate && !bFirstStep && bNS)
        {
            // Copy velocities from the GPU on search steps to keep a copy on host (device buffers are reinitialized).
            stateGpu->copyVelocitiesFromGpu(state->v, AtomLocality::Local);
            stateGpu->waitVelocitiesReadyOnHost(AtomLocality::Local);
            // Copy coordinate from the GPU when needed at the search step.
            // NOTE: The cases when coordinates needed on CPU for force evaluation are handled in sim_utils.
            // NOTE: If the coordinates are to be written into output file they are also copied separately before the output.
            stateGpu->copyCoordinatesFromGpu(state->x, AtomLocality::Local);
            stateGpu->waitCoordinatesReadyOnHost(AtomLocality::Local);
        }

        if (bNS && !(bFirstStep && ir->bContinuation))
        {
            bMasterState = FALSE;
            /* Correct the new box if it is too skewed */
            if (inputrecDynamicBox(ir))
            {
                if (correct_box(fplog, step, state->box))
                {
                    bMasterState = TRUE;
                    // If update is offloaded, it should be informed about the box size change
                    if (useGpuForUpdate)
                    {
                        integrator->setPbc(PbcType::Xyz, state->box);
                    }
                }
            }
            if (DOMAINDECOMP(cr) && bMasterState)
            {
                dd_collect_state(cr->dd, state, state_global);
            }

            if (DOMAINDECOMP(cr))
            {
                /* Repartition the domain decomposition */
                dd_partition_system(fplog, mdlog, step, cr, bMasterState, nstglobalcomm, state_global,
                                    *top_global, ir, imdSession, pull_work, state, &f, mdAtoms, &top,
                                    fr, vsite, constr, nrnb, wcycle, do_verbose && !bPMETunePrinting);
                shouldCheckNumberOfBondedInteractions = true;
                upd.setNumAtoms(state->natoms);

                /* PLUMED */
#if (GMX_PLUMED)
                if(plumedswitch){
                  int nat_home = dd_numHomeAtoms(*cr->dd);
                  plumed_cmd(plumedmain,"setAtomsNlocal",&nat_home);
                  plumed_cmd(plumedmain,"setAtomsGatindex",cr->dd->globalAtomIndices.data());
                }
#endif
                /* END PLUMED */
            }
        }

        // Allocate or re-size GPU halo exchange object, if necessary
        if (bNS && havePPDomainDecomposition(cr) && simulationWork.useGpuHaloExchange)
        {
            GMX_RELEASE_ASSERT(fr->deviceStreamManager != nullptr,
                               "GPU device manager has to be initialized to use GPU "
                               "version of halo exchange.");
            constructGpuHaloExchange(mdlog, *cr, *fr->deviceStreamManager, wcycle);
        }

        if (MASTER(cr) && do_log)
        {
            gmx::EnergyOutput::printHeader(fplog, step,
                                           t); /* can we improve the information printed here? */
        }

        if (ir->efep != efepNO)
        {
            update_mdatoms(mdatoms, state->lambda[efptMASS]);
        }

        if (bExchanged)
        {
            /* We need the kinetic energy at minus the half step for determining
             * the full step kinetic energy and possibly for T-coupling.*/
            /* This may not be quite working correctly yet . . . . */
            compute_globals(gstat, cr, ir, fr, ekind, makeConstArrayRef(state->x),
                            makeConstArrayRef(state->v), state->box, mdatoms, nrnb, &vcm, wcycle,
                            enerd, nullptr, nullptr, nullptr, nullptr, constr, &nullSignaller,
                            state->box, &totalNumberOfBondedInteractions, &bSumEkinhOld,
                            CGLO_GSTAT | CGLO_TEMPERATURE | CGLO_CHECK_NUMBER_OF_BONDED_INTERACTIONS);
            checkNumberOfBondedInteractions(mdlog, cr, totalNumberOfBondedInteractions, top_global,
                                            &top, makeConstArrayRef(state->x), state->box,
                                            &shouldCheckNumberOfBondedInteractions);
        }
        clear_mat(force_vir);

        /************
         * TRANSFER *
         ************/

        clock_gettime(CLOCK_MONOTONIC, &time_1); //is here to measure the time of al QM stuff in one step

        // do it all if we have a real simulation, or we want to calculate parameters in this step
        // now with NOMOVEMENT just the propagation is stopped.
        // with LIQM internal relaxation is still there every ct->interval step
        // (this if-clause will be referred to as the INITIAL IF)
        if ((ct->jobtype != ctePARAMETERS && ct->jobtype != cteNOMOVEMENT && ct->jobtype != cteESP && ct->jobtype != cteTDA) || step % ct->interval == 0)
     // if ((ct->jobtype != ctePARAMETERS && ct->jobtype != cteESP && ct->jobtype != cteTDA) || step % ct->interval == 0)
        {
            // RESTORE THE ORIGINAL CHARGES ON THE NUCLEOBASES - NO MAPPING OF THE HOLE
            // ... IN ORDER TO AVOID DOUBLE COUNTING OF THE ELECTROSTATIC INTERACTION
            // THIS WILL BE TAKEN INTO ACCOUNT WITHIN THE TIGHT-BINDING HAMILTONIAN!

            for (int i=0; i<ct->sites; i++)
                for (int m=0; m<ct->site[i].atoms; m++)
                    mdatoms->chargeA[ct->site[i].atom[m]] = ct_atoms->atom[ct->site[i].atom[m]].q;    // A.HECK

            /* collect the coordinates from all MPI threads */
            transfer_collect_x(cr, state, state_global);
            if (!PAR(cr) || MASTER(cr))
            {
                for (i=0; i<top_global->natoms; i++)
                    copy_rvec(state_global->x[i], x_ct[i]);
            }

#if GMX_MPI
            if (MASTER(cr))
            {
                do_pbc_mtop(ir->pbcType, state->box, top_global, x_ct);
                // make molecules whole -> no dipole, but now we no longer build pbc for every site
                // edit: for MD we will use site energies obtained in phase1 for CT-Hamiltonian.
                //       in phase 2 different sites have different environments (QM vs MM) and therefore enrgies (very bad!).
                // edit2: now we need again whole molecules to calculate correct COM distance for adaptive QM zones.
                for (int i=1; i<ct_mpi_size; i++)
                {
                    MPI_Send(x_ct, 3 * top_global->natoms, GMX_MPI_REAL, i, 300+i, ct_mpi_comm);
                 // printf("Data of size %d sent to rank %d\n", 3 * top_global->natoms, i);
                }
            }
            else
            {
                int return_value = MPI_Recv(x_ct, 3 * top_global->natoms, GMX_MPI_REAL, 0, 300 + ct_mpi_rank, ct_mpi_comm, &ct_mpi_status);
                if (return_value)
                    printf("Error receiving data on rank %d\n", ct_mpi_rank);
             // else
             //     printf("Data of size %d received on rank %d\n", 3 * top_global->natoms, ct_mpi_rank);
            }
         // for (int i=0; i<10; i++)
         // {
         //     printf("rank%d atom%d %12.7f%12.7f%12.7f\n", ct_mpi_rank, i, x_ct[i][XX]*10, x_ct[i][YY]*10, x_ct[i][ZZ]*10);
         // }
#endif

            // optimize QM zone at the beginning
#if GMX_MPI    
            if (ct->pool_size > ct->sites && ct_mpi_size > 1 )
            {
                printf("Adaptive QM zone works only without MPI\n"); //TODO: implement MPI version (segfaults arise later in adapt_QMzone)
                exit(-1);
            }
#endif
            if(ct->opt_QMzone && ct->first_step)
            {
                transfer_collect_x(cr, state, state_global);
                for (int i=0; i<top_global->natoms; i++)
                copy_rvec(state_global->x[i], x_ct[i]);
                search_starting_site(state->box, mdatoms, fr->nbv.get(), fr->shift_vec, cr, ct, dftbplus_phase1.data(), dftbplus_phase2,
                                     x_ct, top_global, state_global->x.rvec_array());
            }
        
            // TODO: maybe do this AFTER the call to do_force() because the neighbor lists are only ready then?
            printf("start prepare at %f\n", (double) clock()/CLOCKS_PER_SEC);
            prepare_charge_transfer(state->box, mdatoms, bNS, fr->nbv.get(), fr->shift_vec, cr, ct, dftbplus_phase1.data(), dftbplus_phase2, x_ct);
            printf("stop prepare at %f\n", (double) clock()/CLOCKS_PER_SEC);

#if GMX_MPI
            if (ct->hamiltonian_type == 0)
            {
                rvec* force_dftb; // TODO - allocate and check the use
                /* removed -- TODO check if anything similar is needed for DFTB+
                if (ct->qmmm == 3)
                {
                    if (step - dftb->lastlist_pme >= dftb->nstlist_pme)
                    {
                        do_neighborlist_for_dftb(ct, dftb, state->x, ct_mpi_rank, ct_mpi_size);
                        dftb->lastlist_pme = step;
                    }
                    do_pme_for_dftb_part1(ct, dftb, ct_mpi_rank, ct_mpi_size);
                    printf("stop pme-preparation at %f\n", (double) clock()/CLOCKS_PER_SEC);
                }
                */
                if (ct->jobtype != cteESP)
                {
                 // do_dftb_phase1(ct, dftb, ct_mpi_comm, ct_mpi_rank, ct_mpi_size);
                    do_dftb_phase1(ct, dftbplus_phase1, cr, force_dftb, nrnb, wcycle, ct_mpi_comm, ct_mpi_rank, ct_mpi_size);
                }
                else
                {
                    do_esp_only(ct, dftbplus_phase1, cr, mdatoms->chargeA, nrnb, wcycle, ct_mpi_comm, ct_mpi_rank, ct_mpi_size);
                }
                // broadcast the results if all calculations have finished. broadcaster is rank "i % ct_mpi_size" (there was phase1 performed)
                MPI_Barrier(ct_mpi_comm);
                for (i=0; i<ct->sites; i++)
                {
                 // printf("Bcast at rank %d at %f\n", ct_mpi_rank, (double) clock()/CLOCKS_PER_SEC);
                    // needed to build the coarse grained hamilton matrix at main node
                    MPI_Bcast(dftb->phase1[i].a[0], SQR(dftb->phase1[i].norb),          MPI_DOUBLE, i % ct_mpi_size, ct_mpi_comm);
                    // needed in check_and_invert_orbital_phase
                    MPI_Bcast(dftb->phase1[i].overl[0], SQR(dftb->phase1[i].norb),      MPI_DOUBLE, i % ct_mpi_size, ct_mpi_comm);
                    // dtto
                    MPI_Bcast(dftb->phase1[i].qmat, dftb->phase1[i].nn,                 MPI_DOUBLE, i % ct_mpi_size, ct_mpi_comm);
                    // dtto
                    MPI_Bcast(dftb->phase1[i].ev, dftb->phase1[i].norb,                 MPI_DOUBLE, i % ct_mpi_size, ct_mpi_comm);
                 // printf("after Bcast at rank %d at %f\n",i % ct_mpi_size, (double) clock()/CLOCKS_PER_SEC);
                }
                MPI_Barrier(ct_mpi_comm); //wait until broadcasting has finished
          
                if (ct_mpi_rank == 0 && ct->jobtype != cteESP)
                {
                    check_and_invert_orbital_phase(dftb->phase1, ct, state_global, mdatoms);
                    // calc sign on master since we need later also the obtained overlap in project_wf_on_new_basis().
                 // check_and_invert_orbital_phase(dftb->phase1, ct);
                    do_dftb_phase2(ct, dftbplus_phase2, cr, force_dftb, nrnb, wcycle);
                }
            }
#else
            if (ct->hamiltonian_type == 0)
            {
                rvec* force_dftb; // TODO - allocate and check the use
                /* removed -- TODO check if something similar is needed for DFTB+
                if (ct->qmmm == 3)
                {
                     if (step - dftb->lastlist_pme >= dftb->nstlist_pme)
                     {
                         do_neighborlist_for_dftb(ct, dftb, state->x.rvec_array());
                         dftb->lastlist_pme = step;
                     }
                     do_pme_for_dftb_part1(ct, dftb);
                }
                */
                if (ct->jobtype != cteESP)
                {
                    do_dftb_phase1(ct, dftbplus_phase1, cr, force_dftb, nrnb, wcycle);
                    check_and_invert_orbital_phase(dftb->phase1, ct, state_global, mdatoms);
                    // copying of charges and PME for 2nd phase comes here !!!
                    do_dftb_phase2(ct, dftbplus_phase2, cr, force_dftb, nrnb, wcycle);
                }
                else
                {
                    do_esp_only(ct, dftbplus_phase1, cr, mdatoms->chargeA, nrnb, wcycle);
                }
            }
#endif

            ///////////////////////////////////////////
            // write interesting data in PDB file as bfactor in order to combine structre and property.
            /*
            for (int i=0; i<ct->sites; i++)
            for (int j=0; j<ct->site[i].atoms; j++)
            {
                ct_atoms->pdbinfo[ct->site[i].atom[j]].bfac=0;
                for (int k=0; k<ct->site[i].homos; k++){
                    // partial charges
                 // ct_atoms->pdbinfo[ct->site[i].atom[j]].bfac =(real) (dftb->phase1[i].qmat[j] - dftb->qzero1[dftb->phase1[i].izp[j]]);
                    // delta charges
                    ct_atoms->pdbinfo[ct->site[i].atom[j]].bfac +=(real) ct->site[i].delta_q[k][j];
                }
            }
            write_sto_conf("CT.pdb", "written by charge transfer code", ct_atoms, x_ct, NULL, 0, NULL);
            */
            // write_out_MOs(step, x_ct, ct_atoms, dftb, ct);
            //////////////////////////////////////////

#if GMX_MPI
            if (ct_mpi_rank == 0)
            {
#endif

            if (ct->jobtype != cteESP)
            {
                // Calculate Hamiltonian by machine learning
                if (ct->hamiltonian_type)
                {
                 // get_ml_hamiltonian(ct, state_global, mdatoms);
                    for (int i=0; i < dftb->phase2.nn; i++)
                        for (int m=0; m < DIM; m++)
                        {
                            dftb->phase2.grad[i][m] = 0;
                        }
                    get_nn_hamiltonian(ct, state_global, mdatoms, dftb);
          
                    for (int i = 0; i < ct->dim; i++)
                    {
                        for (int j = i; j < ct->dim; j++)
                        {
                         // printf("ML: %f \n",ct->hamiltonian_ml[i][j]);
                            ct->hamiltonian[i][j] = ct->hamiltonian_ml[i][j];
                            ct->hamiltonian[j][i] = ct->hamiltonian[i][j];
                        }
                    }
          
                    printf("sucessfully calculated hamiltonian\n");
                   
                    if (ct->jobtype==cteNOMOVEMENT || ct->force_output == 1)
                    {
                        printf("printing ml hamiltonian\n");
                        fprintf(f_tb_hamiltonian_ml, "%8.2f", t*1000);
                        for (int i = 0; i < ct->dim; i++)
                            for (int m = i; m < ct->dim; m++)
                                fprintf(f_tb_hamiltonian_ml, "%15.8f", ct->hamiltonian[i][m] * HARTREE_TO_EV);
                        fprintf(f_tb_hamiltonian_ml, "\n");
                    }
               }
               else
               {
                   ct_assemble_hamiltonian(ct, dftb); // now including averaging
                   /*
                   for (i=0; i<ct->dim; i++)
                   {
                       ct->hamiltonian[i][i] = (ct->is_hole_transfer==1) ? -dftb->phase2.THamilOrtho[i][i] : dftb->phase2.THamilOrtho[i][i];
                       for (j=i+1; j<ct->dim; j++)
                       {
                           // 1.0 or 1.540 scaling factor according to J. Chem. Phys. 140, 104105 (2014)
                        // ct->hamiltonian[i][j] = ct->hamiltonian[j][i] = dftb->phase2.THamilOrtho[i][j] * ct->offdiag_scaling;
                           ct->hamiltonian[i][j] = ct->hamiltonian[j][i] = (ct->is_hole_transfer==1) ?
                                      -dftb->phase2.THamilOrtho[i][j] * ct->offdiag_scaling:  dftb->phase2.THamilOrtho[i][j] * ct->offdiag_scaling;
                           // 1.0 or 1.540 scaling factor according to J. Chem. Phys. 140, 104105 (2014)
                           // this is the right way to do it. without changing off-diagonal elements
                           //    the eigenvalues of the CG Ham are not necessarily -eigenvalue of electron transfer.
                        }
                    }
                    */
                    if (ct->first_step)
                    {
                        for (int i = 0; i < ct->dim; i++)
                            for (int j = 0; j < ct->dim; j++)
                                ct->positive_coupling_factor[i][j] = 1.0;
          
                     // int na = 12, nb = 15;
                     // for (int i = 0; i < na; i++)
                     //     for (int j = 0; j < nb; j++)
                     //         for (int k = j+1; k < nb; k++)
                     //             ct->positive_coupling_factor[i*nb+j][i*nb+k] = -1.0;
                    }
          
                    for (int i = 0; i < ct->dim; i++)
                    {
                        for (int j = i; j < ct->dim; j++)
                        {
                            ct->hamiltonian[i][j] = ct->hamiltonian_dftb[i][j]*ct->positive_coupling_factor[i][j];
                         // printf("%f\n", ct->hamiltonian[i][j]);
                            ct->hamiltonian[j][i] = ct->hamiltonian[i][j];
                        }
                    }
          
                    if (ct->jobtype==cteNOMOVEMENT || ct->force_output == 1)
                    {
                        fprintf(f_tb_hamiltonian, "%8.2f", t*1000);
                        for (int i = 0; i < ct->dim; i++)
                            for (int m = i; m < ct->dim; m++)
                                fprintf(f_tb_hamiltonian, "%15.8f", ct->hamiltonian[i][m] * HARTREE_TO_EV);
                        fprintf(f_tb_hamiltonian, "\n");
                    }
                } // if Hamiltonian type
            } // if jobtype != ESP

            //overlap ///////////////////////////////////////////////////////

            // write overlap matrix for first site
         // {
         //     int i=0;
         //     fprintf(f_ct_overlap_fo, "%10d ", step);
         //     for (int j = 0; j < ct->site[i].homos; j++)
         //         for (int k = j; k < ct->site[i].homos; k++)
         //             fprintf(f_ct_overlap_fo, " %10.6f ", ct->site[i].overlap[j][k]);
         //     fprintf(f_ct_overlap_fo, "  \n");
         //
         //  // int i=0;
         //     fprintf(f_ct_overlap_fo_ref, "%10d ", step);
         //     for (int j = 0; j < ct->site[i].homos; j++)
         //         for (int k = j; k < ct->site[i].homos; k++)
         //             fprintf(f_ct_overlap_fo_ref, " %10.6f ", ct->site[i].overlap_ref[j][k]);
         //     fprintf(f_ct_overlap_fo_ref, "  \n");
         // }

            //spectrum////////////////////
            if (ct->jobtype==ctePARAMETERS) // && ct->first_step == 1)
            {
                get_spectrum(ct,dftb);
             // fprintf(f_ct_spec, "energy eigen values of FO Hamiltonian:\n");
                fprintf(f_ct_spec, "%10ld", step);
                fprintf(f_ct_spec_evec, "%10ld", step);

                for (int i=0; i<ct->dim; i++)
                    fprintf(f_ct_spec, " %10.6f ", ct->ev_spec[i] * HARTREE_TO_EV);
                fprintf(f_ct_spec, "\n");
               
                for (int i=0; i<ct->dim; i++)
                    for (int j=0; j<ct->dim; j++)
                        fprintf(f_ct_spec_evec, " %10.6f ", ct->evec_spec[j+i*ct->dim] ) ; //prints all evecs in one line
                fprintf(f_ct_spec_evec, "\n");
            }
/////////////////////////////////////////////////////////////////////


            // CONSTRUCT THE HUBBARD MATRIX
            // now here because hubbard matrix depends on delta_q calculated by get_delta_q.

            if (ct->jobtype==cteSCCDYNAMIC || ct->jobtype==cteCPFSCCDYNAMIC || ct->jobtype==cteDDBSCCDYNAMIC || ct->jobtype==cteADIABATIC ||
                ct->jobtype==cteNOMOVEMENT || ct->jobtype == cteFERMIADIABATIC || ct->jobtype == cteBORNOPPENHEIMER || ct->jobtype==cteDLZSH ||
                ct->jobtype==cteALZSH || ct->jobtype==cteDFSSH || ct->jobtype==cteJFSSH || ct->jobtype==cteBCJFSSH || ct->jobtype==cteSCRDFSSH ||
                ct->jobtype==cteCCFSSH || ct->jobtype==cteGFSH || ct->jobtype==cteDISH || ct->jobtype==cteSCRDFSSH || ct->jobtype==cteTULLYLOC ||
                ct->jobtype==ctePERSICOSFHOPPING)
            {
                // diag: Hubbard params of nucleobases (constant - do not calculate)

                // offdiag: 1/r_ij for nucleobases (original code)
             // {
             //     int counter = -1;
             //     for (int i=0; i<ct->sites; i++)
             //         for (int k=0; k<ct->site[i].homos; k++)
             //         {
             //             int counter2=-1;
             //             counter++;
             //             for (j=0; j<ct->sites; j++)
             //             {
             //                 dvec_sub(dftb->phase1[i].com, dftb->phase1[j].com, bond);
             //                 for (int l=0; l<ct->site[j].homos; l++)
             //                 {
             //                     counter2++;
             //                     if (j==i)
             //                         continue;
             //                     // approximation of nucleobases as single charged points
             //                     ct->hubbard[counter][counter2] = ct->hubbard[counter2][counter] = ct->sic / dnorm(bond);
             //                 }
             //             }
             //         }
             // }

                // offdiag: sum_AB[ (dq_iA * dq_jB)/r_AB  (improved code)
                {
                    int counter=-1;
                    for (int i=0; i<ct->sites; i++)
                        for (int k=0; k<ct->site[i].homos; k++)
                        {
                            int counter2=-1;
                            counter++;
                            for (int j=0; j<ct->sites; j++)
                                for (int l=0; l<ct->site[j].homos; l++)
                                {
                                    counter2++;
                                    if (j==i)
                                        continue;
                                    ct->hubbard[counter][counter2]=0;
                                    for (int m=0; m < ct->site[i].atoms; m++)
                                        for (int n=0; n < ct->site[j].atoms; n++)
                                        {
                                            dvec bond;
                                            // TODO: use the QM system coordinates in dftbplus_phase1[i]->qm instead
                                            dvec_sub(dftb->phase1[i].x[m], dftb->phase1[j].x[n], bond);
                                            ct->hubbard[counter][counter2] += ct->sic * ct->site[i].delta_q[k][m] * ct->site[j].delta_q[l][n] / dnorm(bond);
                                        }
                                    if (ct->do_epol == 1)
                                        ct->hubbard[counter][counter2]= (1.0/EPSILON_OP)*ct->hubbard[counter][counter2];
                                }
                        }
                } // block

                // test: calculate diag elements on the fly    THIS IS WRONG (also overwrites fixed hubbard diagonals eg. lambda_i)
             // {
             //     int counter=0;
             //     for (i=0; i<ct->sites; i++)
             //         for (j=0; j<ct->site[i].homos; j++)
             //         {
             //             ct->hubbard[counter][counter] = 0.0;
             //             for (k=0; k<ct->site[i].atoms; k++)
             //                 for (l=k; l<ct->site[i].atoms; l++)
             //                     ct->hubbard[counter][counter] += ct->site[i].delta_q[j][k]*ct->site[i].delta_q[j][l] *
             //                                      (k>l ? dftb->phase1[i].gammamat[k][l] : dftb->phase1[i].gammamat[l][k]);
             //             counter++;
             //         }
             // }

            } // if jobtype
            else
            {
                for (int i=0; i<ct->dim; i++)
                    for (int j=0; j<ct->dim; j++)
                        ct->hubbard[i][j] = 0.e0;
            } // else (if jobtype)

         // printf("\n HUBBARD MATRIX \n");
         // for (int i=0; i< ct->dim; i++)
         // {
         //     for (int j=0; j< ct->dim; j++)
         //         printf("%f ", ct->hubbard[i][j]);
         //     printf("\n");
         // }

            if (ct->do_projection)
            {
                printf("start project at %f\n", (double) clock()/CLOCKS_PER_SEC);
                project_wf_on_new_basis(step, dftb, ct, f_ct_project_wf, f_ct_project_wf_ref); // TODO: change to DFTB+
                printf("stop project at %f\n", (double) clock()/CLOCKS_PER_SEC);
            }

            if (ct->sic > 0.0 && (ct->jobtype==cteSCCDYNAMIC || ct->jobtype==cteCPFSCCDYNAMIC || ct->jobtype==cteDDBSCCDYNAMIC ||
                                  ct->jobtype==cteADIABATIC || ct->jobtype==cteNOMOVEMENT || ct->jobtype==cteFERMIADIABATIC ||
                                  ct->jobtype==cteBORNOPPENHEIMER || ct->jobtype==cteDLZSH || ct->jobtype==cteALZSH || ct->jobtype==cteDFSSH ||
                                  ct->jobtype==cteJFSSH || ct->jobtype==cteBCJFSSH || ct->jobtype==cteSCRDFSSH || ct->jobtype==cteCCFSSH ||
                                  ct->jobtype==cteGFSH || ct->jobtype==cteDISH || ct->jobtype==cteSCRDFSSH ||  ct->jobtype==cteTULLYLOC ||
                                  ct->jobtype==ctePERSICOSFHOPPING || ct->jobtype==cteNEGFLORENTZ) )
            {
             // fprintf(f_tb_hamiltonian_hub, "%10d", step);
             // fprintf(f_tb_hubbard, "%10d", step);
                for (int i=0; i<ct->dim; i++)
                {
                    ct->hamiltonian_mod[i] = 0.0;
                    for (int m=0; m<ct->dim; m++)
                        ct->hamiltonian_mod[i] += ct->hubbard[i][m] * ct->occupation[m];
                 // fprintf(f_tb_hamiltonian_hub, " %10.6f", ct->hamiltonian_mod[i] * HARTREE_TO_EV);
                 // for (int m=i; m<ct->dim; m++)
                 //     fprintf(f_tb_hubbard, " %10.6f", ct->hubbard[i][m]);
                }
             // fprintf(f_tb_hamiltonian_hub, "\n");
             // fprintf(f_tb_hubbard, "\n");
            }

            if (ct->jobtype==cteADIABATIC || ct->jobtype==cteFERMIADIABATIC || ct->jobtype == cteBORNOPPENHEIMER)
            {
                fprintf(f_ct_adiabatic, "%10ld ", step);
            }

            printf("start dynamic at %f \n", (double) clock()/CLOCKS_PER_SEC);
            ct->do_ncv = 0;
            ct->step = step;
            // start propagation only after the running average is stable. steps starts from zero and we want to perform propagation
            //    in the first step if there is no averaging ct->n_avg=1
            if (step >= ct->n_avg_ham-1)
            switch (ct->jobtype)
            {
            case cteSCCDYNAMIC:
            case cteCPFSCCDYNAMIC:
            case cteDDBSCCDYNAMIC:
            case cteNONSCCDYNAMIC:

                // evaluate the amount of annihilated charge due to neg. imag. potentials
                for (int i=0; i<ct->neg_imag_pot; i++)
                    ct->site_annihilated_occupation[i] += 2 * ct->neg_imag_pot_rate_constant * PS_TO_AU * ct->occupation[ct->site_neg_imag_pot[i]] * ir->delta_t;
                // additional option:  take lowest eigenvector as starting function
                if (ct->adiabstart && ct->first_step )
                {
                    f_ct_startwf = fopen ("CT_STARTWF.xvg", "w" );
                    ct->fermi_kt = BOLTZMANN_HARTREE_KELVIN * 300; // T=300K helps converging
                    if (do_adiab_fermi_onestate(ct, ct_broyden, ct_broyden->df , f_ct_startwf) == 1)
                    {
                        printf("ERROR: calculation of lowest eigenvalue did not converge.\n");
                        exit(-1);
                    }
                    fclose(f_ct_startwf);
                }
                // print out the wave function
                if ((ir->nstxout > 0 && step % ir->nstxout == 0) || (ir->nstxout_compressed > 0 && step % ir->nstxout_compressed == 0))
                {
                    printf("Wave function before step %ld:\n", step);
                    for (int i=0; i<ct->dim; i++)
                        printf("%10.7f %10.7f\n", ct->wf[i], ct->wf[i + ct->dim]);
                }

                // SCC dynamics -- perform the integration with RKsuite
                // (the neg-imag-potential capability is included in do_rksuite() )

                if (ct->do_projection)
                {
                    printf("start project at %f\n", (double) clock()/CLOCKS_PER_SEC);
                    project_wf_on_new_basis(step, dftb, ct, f_ct_project_wf, f_ct_project_wf_ref);
                    printf("stop project at %f\n", (double) clock()/CLOCKS_PER_SEC);
                }

                // Weiwei, coherence penalty functional and decoherence detailed balance
                if (ct->jobtype==cteCPFSCCDYNAMIC || ct->jobtype==cteDDBSCCDYNAMIC)
                {
                    for (int i = 0; i < ct->dim; i++)
                        ct->occupation[i] = SQR(ct->wf[i]) + SQR(ct->wf[i+ct->dim]);

                 // if (ct->jobtype==cteDDBSCCDYNAMIC)
                 // {
                 //     for (int i=0; i<ct->dim; i++)
                 //     {
                 //         for (int j=i+1; j<ct->dim; j++)
                 //         {
                 //             double de = ct->hamiltonian[i][i] - ct->hamiltonian[j][j];
                 //
                 //             double hamtmp1 = sqrt(ct->occupation[j]*2.0/(1.0 + exp( de/ct->fermi_kt)));
                 //             double hamtmp2 = sqrt(ct->occupation[i]*2.0/(1.0 + exp(-de/ct->fermi_kt)));
                 //
                 //             ct->hamiltonian[i][j] *= fabs(hamtmp1 - hamtmp2);
                 //             ct->hamiltonian[j][i]  = ct->hamiltonian[i][j];
                 //         }
                 //     }
                 // }

                    // coherence penalty functional correction
                    double ekintot = 0.0;
                    for (int iatom = 0; iatom < mdatoms->homenr; iatom++)
                    {
                        for (int k = 0; k < ct->sites; k++)
                        {
                            for (int l = 0; l < ct->site[k].atoms; l++)
                            {
                                if (iatom == ct->site[k].atom[l])
                                {
                                    for (int m = 0; m < DIM; m++)
                                        ekintot += 0.5 * ct->site[k].mass[l] * AMU_TO_AU * SQR(state_global->v[iatom][m] * NM_TO_BOHR/PS_TO_AU);
                                     // ekintot += 0.5 * dftb->phase1[k].mass[l] * AMU_TO_AU * SQR(state_global->v[iatom][m] * NM_TO_BOHR/PS_TO_AU);
                                }
                            }
                        }
                    }

                 // ct->dephasing_rate = 1.0/20.0*0.0241;

                    for (int i = 0; i < ct->dim; i++)
                    {
                        for (int j = 0; j < ct->dim; j++)
                        {
                            if (j != i)
                            {
                             // double decay_rate = dephasing_rate;
                                double decay_rate = fabs(ct->hamiltonian[i][i] - ct->hamiltonian[j][j])/(1.0 + 0.1/ekintot);
                                ct->hamiltonian[i][j] += ct->occupation[i] * ct->occupation[j] * decay_rate;
                            }
                        }
                    }
                } // coherence penalty

                // solve TDSE
                /*
                if (step == 0)
                {
                    for (int i = 0; i < 2*ct->dim; i++)
                        ct->tfs_diab[i] = ct->wf[i];
                }
                else
                {
                    double tdse_timestep = 1.e-4 * PS_TO_AU;
                    int nstep = (int) (ct->rk_timestep/tdse_timestep);
                    double rk_timestep0 = ct->rk_timestep;
                    ct->rk_timestep = tdse_timestep;

                    for (int istep = 0; istep < nstep; istep++)
                    {
                        // linear interpolation for Hamiltonian
                        double cur_weight = (double) istep / (double) nstep;
                        for (int i = 0; i < ct->dim; i++)
                            for (int j = 0; j < ct->dim; j++)
                                ct->tfl_mean_ham_full[i][j] = (1.0 - cur_weight)*ct->tfl_old_ham[i*ct->dim+j] + cur_weight*ct->hamiltonian[i][j];

                        // propagate the TDSE in the diabatic representation
                        do_rksuite_diab(ct);
                     // RK5th(ct);
                    }
                    ct->rk_timestep = rk_timestep0;
                }

                for (int i = 0; i < 2*ct->dim; i++)
                    ct->wf[i] = ct->tfs_diab[i];

                // store hamiltonian as old_ham
                for (int i = 0; i < SQR(ct->dim); i++)
                    ct->tfl_old_ham[i] = *(ct->hamiltonian[0]+i);
                */

                do_rksuite(ct);

                // the coefficients of decomposition into adiabatic states
             // fprintf(f_ct_exp_adiab, "%10d ", step);
             // do_wf_decomp_evec(ct, ct_diis, f_ct_exp_adiab);

                if (ct->jobtype==cteCPFSCCDYNAMIC || ct->jobtype==cteDDBSCCDYNAMIC)
                {
                    for (int i=0; i<ct->dim; i++)
                        for (int j=i+1; j<ct->dim; j++)
                        {
                            ct->hamiltonian[i][j] = ct->hamiltonian_type == 0 ? ct->hamiltonian_dftb[i][j] : ct->hamiltonian_ml[i][j];
                            ct->hamiltonian[j][i] = ct->hamiltonian[i][j];
                        }
                }

             // do_wf_decomp_evec_fermi(ct, ct_broyden, ct_broyden->df, f_ct_exp_adiab);
                break;

            case cteADIABATIC:
            case cteADNONSCC:

                // adiabatic dynamics -- find the wavefunction that minimizes the energy
             // do_adiabatic(ct, ct_diis, f_ct_adiabatic);
                if (do_adiabatic(ct, ct_diis, f_ct_adiabatic))
                {
                    fprintf(stderr, "Adiabatic - DIIS not converged in step %ld\n", step);
                 // fprintf(f_ct_adiabatic, "Adiabatic - Broyden not converged\n");
                }
                break;

            case cteFERMIADIABATIC:

                if (do_adiab_fermi_onestate(ct, ct_broyden, ct_broyden->df, f_ct_adiabatic))
                {
                     fprintf(stderr, "Fermi distribution based Broyden not converged in step %ld\n", step);
                }
                break;

            case cteBORNOPPENHEIMER:   //JJK

                // scale timestep with overlap^2. i.e. take smaller steps if wf is changing strongly
                do_born_oppenheimer(ct, ct_broyden, ct_broyden->df, f_ct_adiabatic);

                // scale timestep with overlap^2. i.e. take smaller steps if wf is changing strongly
             // ir->delta_t = do_born_oppenheimer(ct, ct_broyden, ct_broyden->df , f_ct_adiabatic);
             // printf("delta_t %f\n", ir->delta_t);

             // if(ir->delta_t < 0)
             // {
             //     fprintf(stderr, "Fermi distribution based Broyden not converged in step %d\n", step);
             // }
             // else
             // {
             //     ir->delta_t *=original_delta_t;
             //     printf("original delta_t %f\n", original_delta_t);
             // }

                // maybe problematic since       t         = t0 + step*ir->delta_t;
             // ct->rk_timestep = ir->delta_t * PS_TO_AU;
                break;

            // Landau-Zener surface hopping method (diabatic LZ formula)
            case cteDLZSH:
                printf("Beginning propagation and hopping routine\n");
                if (step == 0)
                    fprintf(f_ct_surfacehopping,"%8s %16s %16s %16s %20s \n","#Time","InitSurf","HopSurf","Prob","RandNum");
                do_lzsh_diab(step, ct, dftb, state_global, mdatoms, x_ct, f_ct_surfacehopping);
                break;

            // Landau-Zener surface hopping method (adiabatic LZ formula)
            case cteALZSH:
                printf("Beginning propagation and hopping routine\n");
                if (step == 0)
                    fprintf(f_ct_surfacehopping,"%8s %16s %16s %16s %20s \n","#Time","InitSurf","HopSurf","Prob","RandNum");
                do_lzsh_adiab(step, ct, f_ct_surfacehopping);
                break;

            // diabatic FSSH
            case cteDFSSH:
                printf("Beginning propagation and hopping routine\n");
                if (step == 0)
                    fprintf(f_ct_surfacehopping,"%8s %16s %16s %16s %20s \n","#Time","InitSurf","HopSurf","Prob","RandNum");
                do_fssh_diab(step, ct, dftb, state_global, mdatoms, f_ct_surfacehopping);
                break;

            // FSSH by Jochen Blumberger
            case cteJFSSH:
                printf("Beginning propagation and hopping routine\n");
                if (step == 0)
                    fprintf(f_ct_surfacehopping,"%8s %16s %16s %16s %20s \n","#Time","InitSurf","HopSurf","Prob","RandNum");
#if GMX_MPI
                do_fssh_johen(step, ct, dftb, state_global, mdatoms, f_ct_surfacehopping, ct_mpi_comm, ct_mpi_rank, ct_mpi_size);
#else
                do_fssh_johen(step, ct, dftb, state_global, mdatoms, f_ct_surfacehopping);
#endif
                break;

            // Bolzmann corrected JFSSH
            case cteBCJFSSH:
                printf("Beginning propagation and hopping routine\n");
                if (step == 0)
                    fprintf(f_ct_surfacehopping,"%8s %16s %16s %16s %20s \n","#Time","InitSurf","HopSurf","Prob","RandNum");
#if GMX_MPI
                do_fssh_johen(step, ct, dftb, state_global, mdatoms, f_ct_surfacehopping, ct_mpi_comm, ct_mpi_rank, ct_mpi_size);
#else
                do_fssh_johen(step, ct, dftb, state_global, mdatoms, f_ct_surfacehopping);
#endif
                break;

            // self-consistent restricted-decoherence FSSH by Wang
            case cteSCRDFSSH:
                printf("Beginning propagation and hopping routine\n");
                if (step == 0)
                    fprintf(f_ct_surfacehopping,"%8s %16s %16s %16s %20s \n","#Time","InitSurf","HopSurf","Prob","RandNum");
#if GMX_MPI
                do_self_consistent_restricted_decoherence_fssh(step, ct, dftb, state_global, mdatoms, f_ct_surfacehopping, ct_mpi_comm, ct_mpi_rank, ct_mpi_size);
#else
                do_self_consistent_restricted_decoherence_fssh(step, ct, dftb, state_global, mdatoms, f_ct_surfacehopping);
#endif
                break;

            // crossing-corrected FSSH by Wang
            case cteCCFSSH:
                printf("Beginning propagation and hopping routine\n");
                if (step == 0)
                    fprintf(f_ct_surfacehopping,"%8s %16s %16s %16s %20s \n","#Time","InitSurf","HopSurf","Prob","RandNum");
#if GMX_MPI
                do_crossing_corrected_fssh(step, ct, dftb, state_global, mdatoms, f_ct_surfacehopping, ct_mpi_comm, ct_mpi_rank, ct_mpi_size);
#else
                do_crossing_corrected_fssh(step, ct, dftb, state_global, mdatoms, f_ct_surfacehopping);
#endif
                break;

            // global flux FSSH by prezhdo
            case cteGFSH:
                printf("Beginning propagation and hopping routine\n");
                if (step == 0)
                    fprintf(f_ct_surfacehopping,"%8s %16s %16s %16s %20s \n","#Time","InitSurf","HopSurf","Prob","RandNum");
#if GMX_MPI
                do_global_flux_fssh(step, ct, dftb, state_global, mdatoms, f_ct_surfacehopping, ct_mpi_comm, ct_mpi_rank, ct_mpi_size);
#else
                do_global_flux_fssh(step, ct, dftb, state_global, mdatoms, f_ct_surfacehopping);
#endif
                break;

            // decoherence induced FSSH by prezhdo
            case cteDISH:
                printf("Beginning propagation and hopping routine\n");
                if (step == 0)
                    fprintf(f_ct_surfacehopping,"%8s %16s %16s %16s %20s \n","#Time","InitSurf","HopSurf","Prob","RandNum");
#if GMX_MPI
                do_decoherence_induced_fssh(step, ct, dftb, state_global, mdatoms, f_ct_surfacehopping, ct_mpi_comm, ct_mpi_rank, ct_mpi_size);
#else
                do_decoherence_induced_fssh(step, ct, dftb, state_global, mdatoms, f_ct_surfacehopping);
#endif
                break;

            case cteTULLYLOC: //new one of JJK
                printf("Beginning propagation and hopping routine\n");
                if (step == 0)
                    fprintf(f_ct_surfacehopping,"%8s %16s %16s %16s %20s \n","#Time","InitSurf","HopSurf","Prob","RandNum");
             // fprintf(f_ct_surfacehopping, "%10d ", step);   //maybe replace step with time, maybe move into function
             // fprintf(f_ct_nonadiab_coupling, "%10d ", step);
             // do_tully_local(ct, ct_broyden, ct_broyden->df, f_ct_surfacehopping, f_ct_nonadiab_coupling);
                do_tully_local(step, ct, dftb, state_global, mdatoms, f_ct_surfacehopping);
                break;

            case ctePERSICOSFHOPPING:
                fprintf(f_ct_surfacehopping, "%10ld ", step);
             // fprintf(f_ct_nonadiab_coupling, "%10d ", step);
             // fprintf(f_ct_state_vectors, "%10ld ", step);
                printf("%10ld ", step);
                // new corrected version. adapted but not tested yet
                do_persico_diabatic_sfhopping_new(ct, ct_broyden, ct_broyden->df, f_ct_surfacehopping, f_ct_nonadiab_coupling, f_ct_state_vectors);
             // do_persico_diabatic_sfhopping(ct, ct_broyden, ct_broyden->df, f_ct_surfacehopping, f_ct_nonadiab_coupling);
                break;

            case cteNEGFLORENTZ:
            case cteNEGFLORENTZNONSCC:
                negf_propagate(ct);
                fprintf(f_ct_current, "%10ld %12.5e %12.5e %12.5e\n", step, ct->negf_arrays->current[0], ct->negf_arrays->current[1], ct->negf_arrays->current[2]);
                break;

            case ctePARAMETERS:
                break;

            case cteNOMOVEMENT:
                if (ct->first_step)
                {
                    ct->survival=0;
                    for (i=0; i<ct->dim; i++)
                        ct->survival += SQR(ct->wf[i]) + SQR(ct->wf[i+ct->dim]);
                    if (ct->survival < 0.99) // should be normalized, 0.01 tolerance
                    { 
                        printf("WARNING: no correct starting wave function was specified. Lowest eigenvalue will be taken instead.\n");
                        f_ct_startwf = fopen ("CT_STARTWF.xvg", "w");
                        ct->fermi_kt = BOLTZMANN_HARTREE_KELVIN * 300; // T=300K helps converging
                        if (do_adiab_fermi_onestate(ct, ct_broyden, ct_broyden->df, f_ct_startwf) == 1)
                        {
                            printf("ERROR: calculation of lowest eigenvalue did not converge.\n");
                            exit(-1);
                        }
                        fclose(f_ct_startwf);
                    }
                }
                break;

            case cteESP:
                // nothing to do
                break;

            case cteTDA:
                // call the procedure inside the fprintf command (maybe not so nice :-/ )
                fprintf(f_ct_tda, "%10ld %15.12f \n", step, calc_tda(ct) * HARTREE_TO_EV);
                break;

            case ctePREZHDOSFHOPPING:
                printf("%10ld ", step);
                fprintf(f_ct_surfacehopping, "%10ld ", step);
             // fprintf(f_ct_nonadiab_coupling, "%10d ", step);
             // fprintf(f_ct_state_vectors, "%10ld ", step);
                do_prezhdo_sfhopping(ct, ct_broyden, ct_broyden->df, f_ct_surfacehopping, f_ct_nonadiab_coupling, f_ct_state_vectors);
                break;

            } // switch jobtype
            printf("stop dynamic at %f\n", (double) clock()/CLOCKS_PER_SEC);

            // store auxiliary arrays and arrays for time derivative of FMOs at preceding step - Weiwei X Mar 2019
            /*
            for (i=0; i<SQR(ct->dim); i++)
                dftb->orthogo.sij_old[i] = dftb->orthogo.sij[i];

            for (i=0; i<dftb->phase2.norb; i++)
                for (j=0; j<ct->dim; j++)
                    dftb->orthogo.fmo_old[i][j] = dftb->orthogo.fmo[i][j];

            for (i=0; i<dftb->phase2.norb; i++)
                for (j=0; j<dftb->phase2.norb; j++)
                    dftb->orthogo.sao_old[i][j] = dftb->orthogo.sao[i][j];
            */

            if (ct->jobtype==cteADIABATIC || ct->jobtype==cteFERMIADIABATIC || ct->jobtype == cteBORNOPPENHEIMER) {
                fprintf(f_ct_adiabatic, "\n");
            }

            // calculate and print out the occupations  dholub may18 added tullydia
            if (ct->jobtype==cteSCCDYNAMIC || ct->jobtype==cteCPFSCCDYNAMIC || ct->jobtype==cteDDBSCCDYNAMIC || ct->jobtype==cteNONSCCDYNAMIC ||
                ct->jobtype==cteADIABATIC || ct->jobtype==cteADNONSCC || ct->jobtype==cteFERMIADIABATIC || ct->jobtype == cteBORNOPPENHEIMER ||
                ct->jobtype==cteDLZSH || ct->jobtype==cteALZSH || ct->jobtype==cteDFSSH || ct->jobtype==cteJFSSH || ct->jobtype==cteBCJFSSH ||
                ct->jobtype==cteSCRDFSSH || ct->jobtype==cteCCFSSH || ct->jobtype==cteGFSH || ct->jobtype==cteDISH || ct->jobtype==cteSCRDFSSH ||
                ct->jobtype==cteTULLYLOC || ct->jobtype==ctePERSICOSFHOPPING || ct->jobtype==cteNEGFLORENTZ || ct->jobtype==cteNEGFLORENTZNONSCC ||
                ct->jobtype==ctePREZHDOSFHOPPING)
            {
                ct->survival = 0.0;
                fprintf(f_tb_occupation, "%8.2f", t*1000); //1000=converting ps to fs

                for (int i=0; i<ct->dim; i++)
                {
                    // diabatic occupations
                    if (ct->jobtype==cteDLZSH || ct->jobtype==cteDFSSH )
                    {
                        ct->occupation[i] = i == ct->surface ? 1.0 : 0.0;
                    }
                    else
                    {
                        ct->occupation[i] = SQR(ct->wf[i]) + SQR(ct->wf[i+ct->dim]);
                        ct->survival += ct->occupation[i];
                    }
                    fprintf(f_tb_occupation, "%10.5f", ct->occupation[i]);
                }
                fprintf(f_tb_occupation, "%10.5f", ct->survival);

                if (ct->jobtype==cteSCCDYNAMIC || ct->jobtype==cteCPFSCCDYNAMIC || ct->jobtype==cteDDBSCCDYNAMIC || ct->jobtype==cteNONSCCDYNAMIC)
                    for (i=0; i<ct->neg_imag_pot; i++)
                        fprintf(f_tb_occupation, "%10.5f", ct->site_annihilated_occupation[i]);

                fprintf(f_tb_occupation, "\n");
            }

            if (ct->jobtype==cteDLZSH || ct->jobtype==cteALZSH || ct->jobtype==cteDFSSH || ct->jobtype==cteJFSSH || ct->jobtype==cteBCJFSSH ||
                ct->jobtype==cteSCRDFSSH || ct->jobtype==cteCCFSSH || ct->jobtype==cteGFSH || ct->jobtype==cteDISH  || ct->jobtype==cteSCRDFSSH ||
                ct->jobtype==cteTULLYLOC)
            {
                fprintf(f_tb_occupation_adiab, "%8.2f", t*1000); //1000=converting ps to fs

                for (int i=0; i<ct->dim; i++)
                {
                    double pop;
                    // adiabatic population due to adiabatic occupations
                    if (ct->jobtype==cteDLZSH || ct->jobtype==cteDFSSH )
                    {
                        pop = SQR(ct->tfs_vector[i][ct->surface]);
                    }
                    else
                    {
                        pop = i == ct->surface ? 1.0 : 0.0;
                    }
                    fprintf(f_tb_occupation_adiab, "%10.5f", pop);
                }
                fprintf(f_tb_occupation_adiab, "\n");
            }

            if (ct->jobtype==cteJFSSH || ct->jobtype==cteBCJFSSH || ct->jobtype==cteSCRDFSSH || ct->jobtype==cteCCFSSH || ct->jobtype==cteGFSH ||
                ct->jobtype==cteDISH  || ct->jobtype==cteSCRDFSSH || ct->jobtype==cteTULLYLOC)
            {
                fprintf(f_tb_diapop, "%8.2f", t*1000); //1000=converting ps to fs
                fprintf(f_tb_adiapop, "%8.2f", t*1000); //1000=converting ps to fs
               
                for (int i=0; i<ct->dim; i++) {
                    // adiabatic population due to adiabatic occupations
                    fprintf(f_tb_diapop, "%10.5f", SQR(ct->tfs_diab[i])+SQR(ct->tfs_diab[i+ct->dim]));
                    fprintf(f_tb_adiapop, "%10.5f", SQR(ct->tfs_popul[i])+SQR(ct->tfs_popul[i+ct->dim]));
                }
                fprintf(f_tb_diapop, "\n");
                fprintf(f_tb_adiapop, "\n");
            }

            if (ct->jobtype==cteNOMOVEMENT)
            {
                // calculate and print TB_COM
                // print out the center of mass (in angstrom) of each site.
                // needed for following the charge through amorphous materials
                fprintf(f_tb_com, "%8.2f", t*1000);
                for (int i=0; i<ct->sites; i++)
                    // coordinates of COM in Angstrom
                    fprintf(f_tb_com, "%15.5f %15.5f %15.5f   ",
                        10*ct->site[i].centerOfMass[0], 10*ct->site[i].centerOfMass[1], 10*ct->site[i].centerOfMass[2]);
                fprintf(f_tb_com, "\n");
            }

            if (ct->jobtype==cteSCCDYNAMIC || ct->jobtype==cteCPFSCCDYNAMIC || ct->jobtype==cteDDBSCCDYNAMIC || ct->jobtype==cteNONSCCDYNAMIC ||
                ct->jobtype==cteADIABATIC || ct->jobtype==cteADNONSCC || ct->jobtype==cteFERMIADIABATIC || ct->jobtype == cteBORNOPPENHEIMER ||
                ct->jobtype==cteDLZSH || ct->jobtype==cteALZSH || ct->jobtype==cteDFSSH || ct->jobtype==cteJFSSH || ct->jobtype==cteBCJFSSH ||
                ct->jobtype==cteSCRDFSSH || ct->jobtype==cteCCFSSH || ct->jobtype==cteGFSH || ct->jobtype==cteDISH  || ct->jobtype==cteSCRDFSSH ||
                ct->jobtype==cteTULLYLOC || ct->jobtype==ctePERSICOSFHOPPING || ct->jobtype==cteNEGFLORENTZ || ct->jobtype==cteNEGFLORENTZNONSCC ||
                ct->jobtype==ctePREZHDOSFHOPPING)
            {
                // Inverse participation ratio from Wang (J.Chem.Phys. 148, 104106, 2018)
                fprintf(f_tb_spread, "%10.3f", t*1000); //1000=converting ps to fs

                double occup_sqr_sum = 0.;

                for (int i=0; i<ct->dim; i++)
                {
                    if (ct->jobtype==cteDFSSH)
                        ct->occupation[i] = SQR(ct->tfs_diab[i]) + SQR(ct->tfs_diab[i+ct->dim]);
                    else
                        ct->occupation[i] = SQR(ct->wf[i]) + SQR(ct->wf[i+ct->dim]);
                    occup_sqr_sum += SQR(ct->occupation[i]);
                }
                fprintf(f_tb_spread, "%11.5f \n", 1./occup_sqr_sum);

                // IPR-Jochen (J.Phys.Chem.Lett. 9, 3116, 2018)
             // double occup_sqr_sum = 0.;
             // for (int i=0; i<ct->dim; i++)
             // {
             //     occup_sqr_sum += SQR(SQR(ct->tfs_diab[i]) + SQR(ct->tfs_diab[i+ct->dim]));
             // }
             // fprintf(f_tb_spread, "%12.8f", 1./occup_sqr_sum);
            }

            // print out the ESP
            if (ct->qmmm > 0)
            {
                fprintf(f_ct_esp, "%10ld", step);
                for (int i=0; i<ct->sites; i++)
                {
                    fprintf(f_ct_esp, " %9.5f", dftb->phase1[i].esp * AU_OF_ESP_TO_VOLT);
                    // For comparison (debug) purposes,
                    // calculate the SCC-DFTB external shift
                    // (be averaging over the atoms, mass weighted)
                    // Maybe output this instead of ESP!
                 // shiftE = 0.0;
                 // for (j=0; j<dftb->phase1[i].nn; j++)
                 //   shiftE += dftb->phase1[i].mass[j] * dftb->phase1[i].shiftE[j];
                 // shiftE *= dftb->phase1[i].inv_tot_mass;
                 // fprintf(f_ct_shift, " %9.5f", shiftE * AU_OF_ESP_TO_VOLT);
                }
                fprintf(f_ct_esp, "\n");
            }

            // print out center, standard deviation (width) and mean square displacement (needed for mobility) of the charge
            if (ct->jobtype==cteSCCDYNAMIC || ct->jobtype==cteCPFSCCDYNAMIC || ct->jobtype==cteDDBSCCDYNAMIC || ct->jobtype==cteNONSCCDYNAMIC ||
                ct->jobtype==cteADIABATIC || ct->jobtype==cteADNONSCC || ct->jobtype==cteFERMIADIABATIC || ct->jobtype == cteBORNOPPENHEIMER ||
                ct->jobtype==cteDLZSH || ct->jobtype==cteALZSH || ct->jobtype==cteDFSSH || ct->jobtype==cteJFSSH || ct->jobtype==cteBCJFSSH ||
                ct->jobtype==cteSCRDFSSH || ct->jobtype==cteCCFSSH || ct->jobtype==cteGFSH || ct->jobtype==cteDISH  || ct->jobtype==cteSCRDFSSH ||
                ct->jobtype==cteTULLYLOC || ct->jobtype==ctePERSICOSFHOPPING || ct->jobtype==cteNEGFLORENTZ || ct->jobtype==cteNEGFLORENTZNONSCC ||
                ct->jobtype==ctePREZHDOSFHOPPING)
            {
                // COC: <psi|R|psi>                             = sum_i occ_i(t)*R_i(t)             // R is operator that gives COM of site i (R_i)
                // MSD: <psi(t)|(R-R0)^2|psi(t)>                = sum_i occ_i(t)*[R_i(t)-*R0)]^2    // R0 is position (COC) of the WF at t=0
                // MSD-Johen: <psi(t)|R|psi(t)>^2               = sum_i occ_i(t)*[R_i(t)]^2
                // MSD-Wang: <psi_adia(t)|(R-R0)^2|psi_adia(t)> = sum_i vec_i(t)^2*[R_i(t)-*R0)]^2  // R0 is position (COC) of the WF at t=0
                // STD: sqrt(<psi|R^2|psi> - (<psi|R|psi>)^2)   = sqrt(  sum_i [occ_i(t)*R_i^2(t)- (occ_i(t)*R_i(t))^2]  )

                dvec std, msd;

                for (int i=0; i<ct->dim; i++)
                ct->occupation[i] = SQR(ct->wf[i]) + SQR(ct->wf[i+ct->dim]);
                int counter=0;
                clear_dvec(ct->coc);

                for (int i=0; i<ct->sites; i++)
                {
                    for (int j=0; j<ct->site[i].homos; j++){
                       ct->coc[XX]+=ct->occupation[counter]*ct->site[i].centerOfMass[XX];
                       ct->coc[YY]+=ct->occupation[counter]*ct->site[i].centerOfMass[YY];
                       ct->coc[ZZ]+=ct->occupation[counter]*ct->site[i].centerOfMass[ZZ];
                       counter++;
                    }
                }
                if (ct->first_step)
                    copy_dvec(ct->coc, ct->coc_start);
             
                counter=0;
                clear_dvec(std);
                clear_dvec(msd);

                // MSD
                msd[XX] = SQR(ct->coc[XX] - ct->coc_start[XX]);
                msd[YY] = SQR(ct->coc[YY] - ct->coc_start[YY]);
                msd[ZZ] = SQR(ct->coc[ZZ] - ct->coc_start[ZZ]);

                /*
                for (int i=0; i<ct->sites; i++)
                    for (int j=0; j<ct->site[i].homos; j++)
                    {
                        msd[XX]+=ct->occupation[counter]*SQR(dftb->phase1[i].com[XX]-ct->coc_start[XX]);
                        msd[YY]+=ct->occupation[counter]*SQR(dftb->phase1[i].com[YY]-ct->coc_start[YY]);
                        msd[ZZ]+=ct->occupation[counter]*SQR(dftb->phase1[i].com[ZZ]-ct->coc_start[ZZ]);
                        counter++;
                    }
                */
                fprintf(f_tb_wfprops, "%8.2f  %15.5f %15.5f  %15.5f %15.5f %15.5f %15.5f", t*1000,
                                      ct->coc[XX], ct->coc[YY], ct->coc[ZZ], msd[XX], msd[YY], msd[ZZ]); // TODO: note - this is now in nm units

                clear_dvec(msd);

                for (int i=0; i<ct->sites; i++)
                    for (int j=0; j<ct->site[i].homos; j++)
                    {
                        msd[XX]+=ct->occupation[counter]*SQR(ct->site[i].centerOfMass[XX] - ct->coc_start[XX]);
                        msd[YY]+=ct->occupation[counter]*SQR(ct->site[i].centerOfMass[YY] - ct->coc_start[YY]);
                        msd[ZZ]+=ct->occupation[counter]*SQR(ct->site[i].centerOfMass[ZZ] - ct->coc_start[ZZ]);
                        counter++;
                    }
                fprintf(f_tb_wfprops, "%15.5f %15.5f %15.5f \n", msd[XX], msd[YY], msd[ZZ]); // TODO: note - this is now in nm units
            }

#if GMX_MPI
            } // if (ct_mpi_rank == 0)

            if (ct->jobtype==cteSCCDYNAMIC || ct->jobtype==cteCPFSCCDYNAMIC || ct->jobtype==cteDDBSCCDYNAMIC || ct->jobtype==cteADIABATIC ||
                ct->jobtype==cteNONSCCDYNAMIC || ct->jobtype==cteADNONSCC || ct->jobtype==cteFERMIADIABATIC || ct->jobtype == cteBORNOPPENHEIMER ||
                ct->jobtype==cteDLZSH || ct->jobtype==cteALZSH || ct->jobtype==cteDFSSH || ct->jobtype==cteJFSSH || ct->jobtype==cteBCJFSSH ||
                ct->jobtype==cteSCRDFSSH || ct->jobtype==cteCCFSSH || ct->jobtype==cteGFSH || ct->jobtype==cteDISH  || ct->jobtype==cteSCRDFSSH ||
                ct->jobtype==cteTULLYLOC || ct->jobtype==ctePERSICOSFHOPPING || ct->jobtype==cteNEGFLORENTZ)
            {
                // communicate the occupations and wavefunction to the slaves
                if (ct_mpi_rank == 0)
                {
                    for (int i=1; i<ct_mpi_size; i++)
                        MPI_Send(ct->occupation, ct->dim, MPI_DOUBLE, i, 200+i, ct_mpi_comm);
                }
                else
                {
                    MPI_Recv(ct->occupation, ct->dim, MPI_DOUBLE, 0, 200+ct_mpi_rank, ct_mpi_comm, &ct_mpi_status);
                    ct->survival = 0.0;
                    for (int i=0; i<ct->dim; i++)
                        ct->survival += ct->occupation[i];
                }
              
                if (ct_mpi_rank == 0)
                {
                    for (int i=1; i<ct_mpi_size; i++)
                        MPI_Send(ct->wf, 2*ct->dim, MPI_DOUBLE, i, 300+i, ct_mpi_comm);
                }
                else
                {
                    MPI_Recv(ct->wf, 2*ct->dim, MPI_DOUBLE, 0, 300+ct_mpi_rank, ct_mpi_comm, &ct_mpi_status);
                }
            }
#endif

            /*
            // if we have neg. imag. potential and the survival is under the threshold, finish!
            if (ct->jobtype==cteSCCDYNAMIC || ct->jobtype==cteCPFSCCDYNAMIC || ct->jobtype==cteDDBSCCDYNAMIC || ct->jobtype==cteNONSCCDYNAMIC)
                if (ct->survival < ct->survival_threshold) 
                {
#if GMX_MPI
                    if (ct_mpi_rank == 0)
#endif
                    fprintf(stderr, "Survival of %f dropped under the threshold (%f) at step %d, terminating!\n",
                                       ct->survival, ct->survival_threshold, step);
                    bGotTermSignal = TRUE;
                }
            */

        } // INITIAL IF

        if (ct->jobtype==cteSCCDYNAMIC || ct->jobtype==cteCPFSCCDYNAMIC || ct->jobtype==cteDDBSCCDYNAMIC || ct->jobtype==cteADIABATIC ||
            ct->jobtype==cteNOMOVEMENT  || ct->jobtype==cteFERMIADIABATIC  || ct->jobtype == cteBORNOPPENHEIMER || ct->jobtype==cteDLZSH ||
            ct->jobtype==cteALZSH || ct->jobtype==cteDFSSH || ct->jobtype==cteJFSSH || ct->jobtype==cteBCJFSSH || ct->jobtype==cteSCRDFSSH ||
            ct->jobtype==cteCCFSSH || ct->jobtype==cteGFSH || ct->jobtype==cteDISH  || ct->jobtype==cteSCRDFSSH || ct->jobtype==cteTULLYLOC ||
            ct->jobtype==ctePERSICOSFHOPPING || ct->jobtype==cteNEGFLORENTZ)
        {
#if GMX_MPI
            get_MM_params(ct, dftb, ct_mpi_comm, ct_mpi_rank, ct_mpi_size);
            MPI_Barrier(ct_mpi_comm);
            for (i=0; i<ct->sites; i++)
            {
                //needed to build the hubbard matrix at main node
                MPI_Bcast(ct->site[i].delta_q[0], ct->site[i].homos*ct->site[i].atoms, MPI_DOUBLE, i % ct_mpi_size, ct_mpi_comm);
                if ( ct->do_lambda_i >= 2 )
                    //needed to add forces for explicit lambda_i //try to send rvec -> 3*double. have to  test
                    MPI_Bcast(dftb->phase1[i].grad[0], 3*dftb->phase1[i].nn,             MPI_DOUBLE, i % ct_mpi_size, ct_mpi_comm);
                if (ct->qmmm == 3)
                    MPI_Bcast(&(dftb->phase1[i].esp), 1,  MPI_DOUBLE, i % ct_mpi_size, ct_mpi_comm);
            }
            if ( ct->do_lambda_i >= 3 )
                MPI_Bcast(dftb->phase2.grad[0],  3*dftb->phase2.nn, MPI_DOUBLE, 0, ct_mpi_comm );
            MPI_Barrier(ct_mpi_comm);
#else
            get_MM_params(ct, dftb);
#endif

            if (! (ct->do_lambda_i == 0 && ct->jobtype==cteNOMOVEMENT))
            {
                // map the charge carrier on the atomic charges - first check the charges
                int counter=0;
                for (int i=0; i<ct->sites; i++)
                {
                 // printf("Site %d, charges top_global and mdatoms\n", i);
                    for (int m=0; m<ct->site[i].atoms; m++)
                    {
                        mdatoms->chargeA[ct->site[i].atom[m]] = ct_atoms->atom[ct->site[i].atom[m]].q;
                     // printf("q %f\n",mdatoms->chargeA[ct->site[i].atom[m]]);
                    }
                    for (int j=0; j<ct->site[i].homos; j++)
                    {
                        for (int m=0; m<ct->site[i].atoms; m++)
                        {
                            mdatoms->chargeA[ct->site[i].atom[m]] += ct->occupation[counter] * ct->site[i].delta_q[j][m];
                        }
                        counter++;
                    }
                 // for (int m=0; m<ct->site[i].atoms; m++)
                 //     printf("rank%d %3d %8.5f %8.5f\n", ct_mpi_rank, ct->site[i].atom[m], ct_atoms.atom[ct->site[i].atom[m]].q,
                 //                                             mdatoms->chargeA[ct->site[i].atom[m]]);
                }
            }
        }

        /* print out active QM zone and set new QM zone if charge moved to the border */
        if (ct->jobtype==cteSCCDYNAMIC || ct->jobtype==cteCPFSCCDYNAMIC || ct->jobtype==cteDDBSCCDYNAMIC || ct->jobtype==cteNONSCCDYNAMIC)
        {
            if(ct->pool_size > ct->sites)
            {
                fprintf(f_ct_active, "%ld  ", step);
                for (int i=0; i<ct->sites; i++)
                    fprintf(f_ct_active, " %d", ct->site[i].resnr);
                fprintf(f_ct_active, "\n");
                while (adapt_QMzone(ct,x_ct, mdatoms, top_global, state->box, state_global->x.rvec_array()))
                {}
            }
        }

        clock_gettime(CLOCK_MONOTONIC, &time_end);
#if GMX_MPI
        if (ct_mpi_rank==0)
#endif
        print_time_difference("TOTAL QM TIME[ns]:", time_1, time_end);
        ct->first_step = 0;

        /****************
         * END TRANSFER *
         ****************/

        checkpointHandler->decideIfCheckpointingThisStep(bNS, bFirstStep, bLastStep);

        /* Determine the energy and pressure:
         * at nstcalcenergy steps and at energy output steps (set below).
         */
        if (EI_VV(ir->eI) && (!bInitStep))
        {
            bCalcEnerStep = do_per_step(step, ir->nstcalcenergy);
            bCalcVir      = bCalcEnerStep
                       || (ir->epc != epcNO
                           && (do_per_step(step, ir->nstpcouple) || do_per_step(step - 1, ir->nstpcouple)));
        }
        else
        {
            bCalcEnerStep = do_per_step(step, ir->nstcalcenergy);
            bCalcVir = bCalcEnerStep || (ir->epc != epcNO && do_per_step(step, ir->nstpcouple));
        }
        bCalcEner = bCalcEnerStep;

        do_ene = (do_per_step(step, ir->nstenergy) || bLastStep);

        if (do_ene || do_log || bDoReplEx)
        {
            bCalcVir  = TRUE;
            bCalcEner = TRUE;
        }

        /* Do we need global communication ? */
        bGStat = (bCalcVir || bCalcEner || bStopCM || do_per_step(step, nstglobalcomm)
                  || (EI_VV(ir->eI) && inputrecNvtTrotter(ir) && do_per_step(step - 1, nstglobalcomm)));

        force_flags = (GMX_FORCE_STATECHANGED | ((inputrecDynamicBox(ir)) ? GMX_FORCE_DYNAMICBOX : 0)
                       | GMX_FORCE_ALLFORCES | (bCalcVir ? GMX_FORCE_VIRIAL : 0)
                       | (bCalcEner ? GMX_FORCE_ENERGY : 0) | (bDoFEP ? GMX_FORCE_DHDL : 0));
        if (fr->useMts && !do_per_step(step, ir->nstfout))
        {
            force_flags |= GMX_FORCE_DO_NOT_NEED_NORMAL_FORCE;
        }

        if (shellfc)
        {
            /* Now is the time to relax the shells */
            relax_shell_flexcon(fplog, cr, ms, mdrunOptions.verbose, enforcedRotation, step, ir,
                                imdSession, pull_work, bNS, force_flags, &top, constr, enerd,
                                state->natoms, state->x.arrayRefWithPadding(),
                                state->v.arrayRefWithPadding(), state->box, state->lambda,
                                &state->hist, &f.view(), force_vir, mdatoms, nrnb, wcycle, shellfc,
                                fr, runScheduleWork, t, mu_tot, vsite, ddBalanceRegionHandler);
        }
        else
        {
            /* The AWH history need to be saved _before_ doing force calculations where the AWH bias
               is updated (or the AWH update will be performed twice for one step when continuing).
               It would be best to call this update function from do_md_trajectory_writing but that
               would occur after do_force. One would have to divide the update_awh function into one
               function applying the AWH force and one doing the AWH bias update. The update AWH
               bias function could then be called after do_md_trajectory_writing (then containing
               update_awh_history). The checkpointing will in the future probably moved to the start
               of the md loop which will rid of this issue. */
            if (awh && checkpointHandler->isCheckpointingStep() && MASTER(cr))
            {
                awh->updateHistory(state_global->awhHistory.get());
            }

            /* The coordinates (x) are shifted (to get whole molecules)
             * in do_force.
             * This is parallellized as well, and does communication too.
             * Check comments in sim_util.c
             */

            /* PLUMED */
#if (GMX_PLUMED)
            plumedNeedsEnergy=0;
            if(plumedswitch){
              int pversion=0;
              plumed_cmd(plumedmain,"getApiVersion",&pversion);
              long int lstep=step; plumed_cmd(plumedmain,"setStepLong",&lstep);
              plumed_cmd(plumedmain,"setPositions",&state->x[0][0]);
              plumed_cmd(plumedmain,"setMasses",&mdatoms->massT[0]);
              plumed_cmd(plumedmain,"setCharges",&mdatoms->chargeA[0]);
              plumed_cmd(plumedmain,"setBox",&state->box[0][0]);
              plumed_cmd(plumedmain,"prepareCalc",NULL);
              plumed_cmd(plumedmain,"setStopFlag",&plumedWantsToStop);
              int checkp=0; if(checkpointHandler->isCheckpointingStep()) checkp=1;
              if(pversion>3) plumed_cmd(plumedmain,"doCheckPoint",&checkp);
        //    plumed_cmd(plumedmain,"setForces",&f[0][0]);
              auto forceView = &f.view();
              plumed_cmd(plumedmain,"setForces",forceView->forceWithPadding().unpaddedArrayRef().data());
              plumed_cmd(plumedmain,"isEnergyNeeded",&plumedNeedsEnergy);
              if(plumedNeedsEnergy) force_flags |= GMX_FORCE_ENERGY | GMX_FORCE_VIRIAL;
              clear_mat(plumed_vir);
              plumed_cmd(plumedmain,"setVirial",&plumed_vir[0][0]);
            }
#endif
            /* END PLUMED */

         // do_force(fplog, cr, ms, ir, awh.get(), enforcedRotation, imdSession, pull_work, step,
         //          nrnb, wcycle, &top, state->box, state->x.arrayRefWithPadding(), &state->hist,
         //          f.arrayRefWithPadding(), force_vir, mdatoms, enerd, fcd, state->lambda, graph,
         //          fr, runScheduleWork, vsite, mu_tot, t, ed ? ed->getLegacyED() : nullptr,
         //          (bNS ? GMX_FORCE_NS : 0) | force_flags, ddBalanceRegionHandler);
            do_force(fplog, cr, ms, ir, awh.get(), enforcedRotation, imdSession, pull_work, step,
                     nrnb, wcycle, &top, state->box, state->x.arrayRefWithPadding(), &state->hist,
                     &f.view(), force_vir, mdatoms, enerd, state->lambda, fr, runScheduleWork,
                     vsite, mu_tot, t, ed ? ed->getLegacyED() : nullptr,
                     (bNS ? GMX_FORCE_NS : 0) | force_flags, ddBalanceRegionHandler);

            /* PLUMED */
#if (GMX_PLUMED)
            if(plumedswitch){
              if(plumedNeedsEnergy){
                msmul(force_vir,2.0,plumed_vir);
                plumed_cmd(plumedmain,"setEnergy",&enerd->term[F_EPOT]);
                plumed_cmd(plumedmain,"performCalc",NULL);
                msmul(plumed_vir,0.5,force_vir);
              } else {
                msmul(plumed_vir,0.5,plumed_vir);
                m_add(force_vir,plumed_vir,force_vir);
              }
              if(bDoReplEx) plumed_cmd(plumedmain,"GREX savePositions",NULL);
              if(plumedWantsToStop) ir->nsteps=step_rel+1;
            //if(bHREX) plumed_cmd(plumedmain,"GREX cacheLocalUNow",&enerd->term[F_EPOT]);
            }
#endif
            /* END PLUMED */
        }

        /************
         * TRANSFER *
         ************/

        // add QM forces to relax site geometry
        // TODO: this way to add forces is not completely correct.
        // There are differences in simulations on single and multiple nodes.
        // This is, however, already the case with normal GROMACS but now more pronounced.
        // forces should be added somewhere INSIDE do_force(), probably do_force_lowlevel()

        {
            auto forceView = &f.view();
            auto forceData = forceView->forceWithPadding().unpaddedArrayRef().data();

            if (ct->do_lambda_i==2 || ct->do_lambda_i==3)
                for (int i=0; i<mdatoms->homenr; i++) //iterate over home atoms
                {
                    for (int j=0; j<ct->sites; j++)
                        for (int k=0; k<ct->site[j].atoms; k++)
                            if (i==ct->site[j].atom[k])
                            { // atom is at home on this node
                                for (int l=0; l<DIM; l++)
                                {
                                    forceData[i][l] -= (real) HARTREE_BOHR2MD * dftb->phase1[j].grad[k][l];
                                 // printf("atom i %d MM force %lf   QM force %lf \n", i, f[i][l], -(real) HARTREE_BOHR2MD * dftb->phase1[j].grad[k][l]);
                                 // fshift[i][j] = (real) HARTREE_BOHR2MD * dftb1.grad[i][j]; TODO: what is fshift? this is used in QMMM gromacs code
                                }
                            }
                }

            if (ct->do_lambda_i==3)
                for(int i=0; i<mdatoms->homenr; i++)
                {
                    int counter=0;
                    for (int j=0; j<ct->sites; j++)
                        for(int k=0; k<ct->site[j].atoms; k++)
                        {
                            if (i==ct->site[j].atom[k])
                            { //atom is at home on this node
                                for (int l=0; l<DIM; l++)
                                {
                                    forceData[i][l] -= (real) HARTREE_BOHR2MD * dftb->phase2.grad[counter][l];
                                 // printf("atom i %d MM force %lf   QM force %lf \n", i, f[i][l], -(real) HARTREE_BOHR2MD * dftb->phase1[j].grad[k][l]);
                                }
                            }
                            counter++;
                        }
                }
            }

            // MK tag printout
            if (ct->jobtype==cteNOMOVEMENT || ct->force_output==1)
            {
                printf("Printing Forces to File\n");
               
                for (int j=0; j<ct->sites; j++)
                {
                 // fprintf(f_tb_diag_force, "%16.8f\n", ct->hamiltonian[j][j]);
                    for (int k=0; k<ct->site[j].atoms; k++)
                    {
                        long glo_idx = ct->site[j].atom[k];
                        int charge;
                        switch (ct->site[j].atomtype[k])
                        {
                            case 0: charge = 6; break; // C
                            case 1: charge = 1; break; // H
                        }
                     // printf("pos %f %f %f\n", state_global->x[glo_idx][0]*10, state_global->x[glo_idx][1]*10,state_global->x[glo_idx][2]*10);
                     // printf("frc %f %f %f\n", dftb->phase1[j].grad[k][0],dftb->phase1[j].grad[k][1], dftb->phase1[j].grad[k][2]);
                        fprintf(f_tb_diag_force, "%1d%16.8f%16.8f%16.8f%16.8f%16.8f%16.8f\n", charge,
                                state_global->x[glo_idx][0]*10, state_global->x[glo_idx][1]*10,state_global->x[glo_idx][2]*10,
                                dftb->phase1[j].grad[k][0],dftb->phase1[j].grad[k][1],dftb->phase1[j].grad[k][2]);
                    }
               
                    fprintf(f_tb_diag_force, "\n");
                }

             // fprintf(f_tb_offdiag_force, "%16.8f\n",ct->hamiltonian[0][1]);

                for (int i=0; i<mdatoms->homenr;i++)
                {
                    int counter=0;
                    for (int j=0; j<ct->sites; j++)
                    {
                        for (int k=0; k<ct->site[j].atoms; k++)
                        {
                            if (i==ct->site[j].atom[k]) //atom is at home on this node
                            {
                                int charge;
                                switch (ct->site[j].atomtype[k])
                                {
                                    case 0: charge = 6; break; // C
                                    case 1: charge = 1; break; // H
                                }
                             // printf("frc %f %f %f\n", dftb->phase2.grad[counter][0],dftb->phase2.grad[counter][1], dftb->phase2.grad[counter][2]);
                                fprintf(f_tb_offdiag_force, "%1d%16.8f%16.8f%16.8f%16.8f%16.8f%16.8f\n",
                                    charge, state_global->x[i][0]*10, state_global->x[i][1]*10, state_global->x[i][2]*10,
                                    dftb->phase2.grad[counter][0], dftb->phase2.grad[counter][1], dftb->phase2.grad[counter][2]);
                            }
                            counter++;
                        }
                    }
                }
                fprintf(f_tb_diag_force, "\n");
                fprintf(f_tb_offdiag_force, "\n");
            }

        /****************
         * END TRANSFER *
         ****************/

        // VV integrators do not need the following velocity half step
        // if it is the first step after starting from a checkpoint.
        // That is, the half step is needed on all other steps, and
        // also the first step when starting from a .tpr file.
        if (EI_VV(ir->eI))
        {
            integrateVVFirstStep(step, bFirstStep, bInitStep, startingBehavior, nstglobalcomm, ir,
                                 fr, cr, state, mdatoms, fcdata, &MassQ, &vcm, top_global, top, enerd,
                                 ekind, gstat, &last_ekin, bCalcVir, total_vir, shake_vir, force_vir,
                                 pres, M, do_log, do_ene, bCalcEner, bGStat, bStopCM, bTrotter,
                                 bExchanged, &bSumEkinhOld, &shouldCheckNumberOfBondedInteractions,
                                 &saved_conserved_quantity, &f, &upd, constr, &nullSignaller,
                                 trotter_seq, nrnb, mdlog, fplog, wcycle);
        }

        /* ########  END FIRST UPDATE STEP  ############## */
        /* ########  If doing VV, we now have v(dt) ###### */
        if (bDoExpanded)
        {
            /* perform extended ensemble sampling in lambda - we don't
               actually move to the new state before outputting
               statistics, but if performing simulated tempering, we
               do update the velocities and the tau_t. */

            lamnew = ExpandedEnsembleDynamics(fplog, ir, enerd, state, &MassQ, state->fep_state,
                                              state->dfhist, step, state->v.rvec_array(), mdatoms);
            /* history is maintained in state->dfhist, but state_global is what is sent to trajectory and log output */
            if (MASTER(cr))
            {
                copy_df_history(state_global->dfhist, state->dfhist);
            }
        }

        // Copy coordinate from the GPU for the output/checkpointing if the update is offloaded and
        // coordinates have not already been copied for i) search or ii) CPU force tasks.
        if (useGpuForUpdate && !bNS && !runScheduleWork->domainWork.haveCpuLocalForceWork
            && (do_per_step(step, ir->nstxout) || do_per_step(step, ir->nstxout_compressed)
                || checkpointHandler->isCheckpointingStep()))
        {
            stateGpu->copyCoordinatesFromGpu(state->x, AtomLocality::Local);
            stateGpu->waitCoordinatesReadyOnHost(AtomLocality::Local);
        }
        // Copy velocities if needed for the output/checkpointing.
        // NOTE: Copy on the search steps is done at the beginning of the step.
        if (useGpuForUpdate && !bNS
            && (do_per_step(step, ir->nstvout) || checkpointHandler->isCheckpointingStep()))
        {
            stateGpu->copyVelocitiesFromGpu(state->v, AtomLocality::Local);
            stateGpu->waitVelocitiesReadyOnHost(AtomLocality::Local);
        }
        // Copy forces for the output if the forces were reduced on the GPU (not the case on virial steps)
        // and update is offloaded hence forces are kept on the GPU for update and have not been
        // already transferred in do_force().
        // TODO: There should be an improved, explicit mechanism that ensures this copy is only executed
        //       when the forces are ready on the GPU -- the same synchronizer should be used as the one
        //       prior to GPU update.
        // TODO: When the output flags will be included in step workload, this copy can be combined with the
        //       copy call in do_force(...).
        // NOTE: The forces should not be copied here if the vsites are present, since they were modified
        //       on host after the D2H copy in do_force(...).
        if (runScheduleWork->stepWork.useGpuFBufferOps && (simulationWork.useGpuUpdate && !vsite)
            && do_per_step(step, ir->nstfout))
        {
            stateGpu->copyForcesFromGpu(f.view().force(), AtomLocality::Local);
            stateGpu->waitForcesReadyOnHost(AtomLocality::Local);
        }
        /* Now we have the energies and forces corresponding to the
         * coordinates at time t. We must output all of this before
         * the update.
         */
        do_md_trajectory_writing(fplog, cr, nfile, fnm, step, step_rel, t, ir, state, state_global,
                                 observablesHistory, top_global, fr, outf, energyOutput, ekind,
                                 f.view().force(), checkpointHandler->isCheckpointingStep(),
                                 bRerunMD, bLastStep, mdrunOptions.writeConfout, bSumEkinhOld);
        /* Check if IMD step and do IMD communication, if bIMD is TRUE. */
        bInteractiveMDstep = imdSession->run(step, bNS, state->box, state->x.rvec_array(), t);

        /* kludge -- virial is lost with restart for MTTK NPT control. Must reload (saved earlier). */
        if (startingBehavior != StartingBehavior::NewSimulation && bFirstStep
            && (inputrecNptTrotter(ir) || inputrecNphTrotter(ir)))
        {
            copy_mat(state->svir_prev, shake_vir);
            copy_mat(state->fvir_prev, force_vir);
        }

        stopHandler->setSignal();
        resetHandler->setSignal(walltime_accounting);

        if (bGStat || !PAR(cr))
        {
            /* In parallel we only have to check for checkpointing in steps
             * where we do global communication,
             *  otherwise the other nodes don't know.
             */
            checkpointHandler->setSignal(walltime_accounting);
        }

        /* #########   START SECOND UPDATE STEP ################# */

        /* at the start of step, randomize or scale the velocities ((if vv. Restriction of Andersen
           controlled in preprocessing */

        if (ETC_ANDERSEN(ir->etc)) /* keep this outside of update_tcouple because of the extra info required to pass */
        {
            gmx_bool bIfRandomize;
            bIfRandomize = update_randomize_velocities(ir, step, cr, mdatoms, state->v, &upd, constr);
            /* if we have constraints, we have to remove the kinetic energy parallel to the bonds */
            if (constr && bIfRandomize)
            {
                constrain_velocities(constr, do_log, do_ene, step, state, nullptr, false, nullptr);
            }
        }
        /* Box is changed in update() when we do pressure coupling,
         * but we should still use the old box for energy corrections and when
         * writing it to the energy file, so it matches the trajectory files for
         * the same timestep above. Make a copy in a separate array.
         */
        copy_mat(state->box, lastbox);

        dvdl_constr = 0;

        if (!useGpuForUpdate)
        {
            wallcycle_start(wcycle, ewcUPDATE);
        }
        /* UPDATE PRESSURE VARIABLES IN TROTTER FORMULATION WITH CONSTRAINTS */
        if (bTrotter)
        {
            trotter_update(ir, step, ekind, enerd, state, total_vir, mdatoms, &MassQ, trotter_seq, ettTSEQ3);
            /* We can only do Berendsen coupling after we have summed
             * the kinetic energy or virial. Since the happens
             * in global_state after update, we should only do it at
             * step % nstlist = 1 with bGStatEveryStep=FALSE.
             */
        }
        else
        {
            update_tcouple(step, ir, state, ekind, &MassQ, mdatoms);
            update_pcouple_before_coordinates(fplog, step, ir, state, pressureCouplingMu, M, bInitStep);
        }

        /* With leap-frog type integrators we compute the kinetic energy
         * at a whole time step as the average of the half-time step kinetic
         * energies of two subsequent steps. Therefore we need to compute the
         * half step kinetic energy also if we need energies at the next step.
         */
        const bool needHalfStepKineticEnergy =
                (!EI_VV(ir->eI) && (do_per_step(step + 1, nstglobalcomm) || step_rel + 1 == ir->nsteps));

        // Parrinello-Rahman requires the pressure to be availible before the update to compute
        // the velocity scaling matrix. Hence, it runs one step after the nstpcouple step.
        const bool doParrinelloRahman = (ir->epc == epcPARRINELLORAHMAN
                                         && do_per_step(step + ir->nstpcouple - 1, ir->nstpcouple));

        if (EI_VV(ir->eI))
        {
            GMX_ASSERT(!useGpuForUpdate, "GPU update is not supported with VVAK integrator.");

            integrateVVSecondStep(step, ir, fr, cr, state, mdatoms, fcdata, &MassQ, &vcm, pull_work,
                                  enerd, ekind, gstat, &dvdl_constr, bCalcVir, total_vir, shake_vir,
                                  force_vir, pres, M, lastbox, do_log, do_ene, bGStat, &bSumEkinhOld,
                                  &f, &cbuf, &upd, constr, &nullSignaller, trotter_seq, nrnb, wcycle);
        }
        else
        {
            if (useGpuForUpdate)
            {

                wallcycle_stop(wcycle, ewcUPDATE);

                if (bNS && (bFirstStep || DOMAINDECOMP(cr)))
                {
                    integrator->set(stateGpu->getCoordinates(), stateGpu->getVelocities(),
                                    stateGpu->getForces(), top.idef, *mdatoms, ekind->ngtc);

                    // Copy data to the GPU after buffers might have being reinitialized
                    stateGpu->copyVelocitiesToGpu(state->v, AtomLocality::Local);
                    stateGpu->copyCoordinatesToGpu(state->x, AtomLocality::Local);
                }

                if (simulationWork.useGpuPme && !runScheduleWork->simulationWork.useGpuPmePpCommunication
                    && !thisRankHasDuty(cr, DUTY_PME))
                {
                    // The PME forces were recieved to the host, so have to be copied
                    stateGpu->copyForcesToGpu(f.view().force(), AtomLocality::All);
                }
                else if (!runScheduleWork->stepWork.useGpuFBufferOps)
                {
                    // The buffer ops were not offloaded this step, so the forces are on the
                    // host and have to be copied
                    stateGpu->copyForcesToGpu(f.view().force(), AtomLocality::Local);
                }

                const bool doTemperatureScaling =
                        (ir->etc != etcNO && do_per_step(step + ir->nsttcouple - 1, ir->nsttcouple));

                // This applies Leap-Frog, LINCS and SETTLE in succession
                integrator->integrate(
                        stateGpu->getForcesReadyOnDeviceEvent(
                                AtomLocality::Local, runScheduleWork->stepWork.useGpuFBufferOps),
                        ir->delta_t, true, bCalcVir, shake_vir, doTemperatureScaling, ekind->tcstat,
                        doParrinelloRahman, ir->nstpcouple * ir->delta_t, M);

                // Copy velocities D2H after update if:
                // - Globals are computed this step (includes the energy output steps).
                // - Temperature is needed for the next step.
                if (bGStat || needHalfStepKineticEnergy)
                {
                    stateGpu->copyVelocitiesFromGpu(state->v, AtomLocality::Local);
                    stateGpu->waitVelocitiesReadyOnHost(AtomLocality::Local);
                }
            }
            else
            {
                /* With multiple time stepping we need to do an additional normal
                 * update step to obtain the virial, as the actual MTS integration
                 * using an acceleration where the slow forces are multiplied by mtsFactor.
                 * Using that acceleration would result in a virial with the slow
                 * force contribution would be a factor mtsFactor too large.
                 */
                if (fr->useMts && bCalcVir && constr != nullptr)
                {
                    upd.update_for_constraint_virial(*ir, *mdatoms, *state,
                                                     f.view().forceWithPadding(), *ekind);

                    constrain_coordinates(constr, do_log, do_ene, step, state,
                                          upd.xp()->arrayRefWithPadding(), &dvdl_constr, bCalcVir,
                                          shake_vir);
                }

                ArrayRefWithPadding<const RVec> forceCombined =
                        (fr->useMts && step % ir->mtsLevels[1].stepFactor == 0)
                                ? f.view().forceMtsCombinedWithPadding()
                                : f.view().forceWithPadding();
                upd.update_coords(*ir, step, mdatoms, state, forceCombined, fcdata, ekind, M,
                                  etrtPOSITION, cr, constr != nullptr);

                wallcycle_stop(wcycle, ewcUPDATE);

                constrain_coordinates(constr, do_log, do_ene, step, state,
                                      upd.xp()->arrayRefWithPadding(), &dvdl_constr,
                                      bCalcVir && !fr->useMts, shake_vir);

                upd.update_sd_second_half(*ir, step, &dvdl_constr, mdatoms, state, cr, nrnb, wcycle,
                                          constr, do_log, do_ene);
                upd.finish_update(*ir, mdatoms, state, wcycle, constr != nullptr);
            }

            if (ir->bPull && ir->pull->bSetPbcRefToPrevStepCOM)
            {
                updatePrevStepPullCom(pull_work, state);
            }

            enerd->term[F_DVDL_CONSTR] += dvdl_constr;
        }

        if (vsite != nullptr)
        {
            wallcycle_start(wcycle, ewcVSITECONSTR);
            vsite->construct(state->x, ir->delta_t, state->v, state->box);
            wallcycle_stop(wcycle, ewcVSITECONSTR);
        }

        /* ############## IF NOT VV, Calculate globals HERE  ############ */
        /* With Leap-Frog we can skip compute_globals at
         * non-communication steps, but we need to calculate
         * the kinetic energy one step before communication.
         */
        {
            // Organize to do inter-simulation signalling on steps if
            // and when algorithms require it.
            const bool doInterSimSignal = (simulationsShareState && do_per_step(step, nstSignalComm));

            if (bGStat || needHalfStepKineticEnergy || doInterSimSignal)
            {
                // Copy coordinates when needed to stop the CM motion.
                if (useGpuForUpdate && !EI_VV(ir->eI) && bStopCM)
                {
                    stateGpu->copyCoordinatesFromGpu(state->x, AtomLocality::Local);
                    stateGpu->waitCoordinatesReadyOnHost(AtomLocality::Local);
                }
                // Since we're already communicating at this step, we
                // can propagate intra-simulation signals. Note that
                // check_nstglobalcomm has the responsibility for
                // choosing the value of nstglobalcomm that is one way
                // bGStat becomes true, so we can't get into a
                // situation where e.g. checkpointing can't be
                // signalled.
                bool                doIntraSimSignal = true;
                SimulationSignaller signaller(&signals, cr, ms, doInterSimSignal, doIntraSimSignal);

                compute_globals(gstat, cr, ir, fr, ekind, makeConstArrayRef(state->x),
                                makeConstArrayRef(state->v), state->box, mdatoms, nrnb, &vcm,
                                wcycle, enerd, force_vir, shake_vir, total_vir, pres, constr,
                                &signaller, lastbox, &totalNumberOfBondedInteractions, &bSumEkinhOld,
                                (bGStat ? CGLO_GSTAT : 0) | (!EI_VV(ir->eI) && bCalcEner ? CGLO_ENERGY : 0)
                                        | (!EI_VV(ir->eI) && bStopCM ? CGLO_STOPCM : 0)
                                        | (!EI_VV(ir->eI) ? CGLO_TEMPERATURE : 0)
                                        | (!EI_VV(ir->eI) ? CGLO_PRESSURE : 0) | CGLO_CONSTRAINT
                                        | (shouldCheckNumberOfBondedInteractions ? CGLO_CHECK_NUMBER_OF_BONDED_INTERACTIONS
                                                                                 : 0));
                checkNumberOfBondedInteractions(mdlog, cr, totalNumberOfBondedInteractions,
                                                top_global, &top, makeConstArrayRef(state->x),
                                                state->box, &shouldCheckNumberOfBondedInteractions);
                if (!EI_VV(ir->eI) && bStopCM)
                {
                    process_and_stopcm_grp(fplog, &vcm, *mdatoms, makeArrayRef(state->x),
                                           makeArrayRef(state->v));
                    inc_nrnb(nrnb, eNR_STOPCM, mdatoms->homenr);

                    // TODO: The special case of removing CM motion should be dealt more gracefully
                    if (useGpuForUpdate)
                    {
                        stateGpu->copyCoordinatesToGpu(state->x, AtomLocality::Local);
                        // Here we block until the H2D copy completes because event sync with the
                        // force kernels that use the coordinates on the next steps is not implemented
                        // (not because of a race on state->x being modified on the CPU while H2D is in progress).
                        stateGpu->waitCoordinatesCopiedToDevice(AtomLocality::Local);
                        // If the COM removal changed the velocities on the CPU, this has to be accounted for.
                        if (vcm.mode != ecmNO)
                        {
                            stateGpu->copyVelocitiesToGpu(state->v, AtomLocality::Local);
                        }
                    }
                }
            }
        }

        /* #############  END CALC EKIN AND PRESSURE ################# */

        /* Note: this is OK, but there are some numerical precision issues with using the convergence of
           the virial that should probably be addressed eventually. state->veta has better properies,
           but what we actually need entering the new cycle is the new shake_vir value. Ideally, we could
           generate the new shake_vir, but test the veta value for convergence.  This will take some thought. */

        if (ir->efep != efepNO && !EI_VV(ir->eI))
        {
            /* Sum up the foreign energy and dK/dl terms for md and sd.
               Currently done every step so that dH/dl is correct in the .edr */
            accumulateKineticLambdaComponents(enerd, state->lambda, *ir->fepvals);
        }

        update_pcouple_after_coordinates(fplog, step, ir, mdatoms, pres, force_vir, shake_vir,
                                         pressureCouplingMu, state, nrnb, upd.deform(), !useGpuForUpdate);

        const bool doBerendsenPressureCoupling =
                (inputrec->epc == epcBERENDSEN && do_per_step(step, inputrec->nstpcouple));
        const bool doCRescalePressureCoupling =
                (inputrec->epc == epcCRESCALE && do_per_step(step, inputrec->nstpcouple));
        if (useGpuForUpdate
            && (doBerendsenPressureCoupling || doCRescalePressureCoupling || doParrinelloRahman))
        {
            integrator->scaleCoordinates(pressureCouplingMu);
            if (doCRescalePressureCoupling)
            {
                matrix pressureCouplingInvMu;
                gmx::invertBoxMatrix(pressureCouplingMu, pressureCouplingInvMu);
                integrator->scaleVelocities(pressureCouplingInvMu);
            }
            integrator->setPbc(PbcType::Xyz, state->box);
        }

        /* ################# END UPDATE STEP 2 ################# */
        /* #### We now have r(t+dt) and v(t+dt/2)  ############# */

        /* The coordinates (x) were unshifted in update */
        if (!bGStat)
        {
            /* We will not sum ekinh_old,
             * so signal that we still have to do it.
             */
            bSumEkinhOld = TRUE;
        }

        if (bCalcEner)
        {
            /* #########  BEGIN PREPARING EDR OUTPUT  ###########  */

            /* use the directly determined last velocity, not actually the averaged half steps */
            if (bTrotter && ir->eI == eiVV)
            {
                enerd->term[F_EKIN] = last_ekin;
            }
            enerd->term[F_ETOT] = enerd->term[F_EPOT] + enerd->term[F_EKIN];

            if (integratorHasConservedEnergyQuantity(ir))
            {
                if (EI_VV(ir->eI))
                {
                    enerd->term[F_ECONSERVED] = enerd->term[F_ETOT] + saved_conserved_quantity;
                }
                else
                {
                    enerd->term[F_ECONSERVED] = enerd->term[F_ETOT] + NPT_energy(ir, state, &MassQ);
                }
            }
            /* #########  END PREPARING EDR OUTPUT  ###########  */
        }

        /************
         * TRANSFER *
         ************/

        // calculate energies
#if GMX_MPI
        if (ct_mpi_rank == 0) {
#endif
        if (ct->jobtype==cteNONSCCDYNAMIC || ct->jobtype==cteADNONSCC || ct->jobtype==cteNOMOVEMENT || ct->jobtype==cteSCCDYNAMIC ||
            ct->jobtype==cteCPFSCCDYNAMIC || ct->jobtype==cteDDBSCCDYNAMIC || ct->jobtype==cteADIABATIC || ct->jobtype==cteFERMIADIABATIC ||
            ct->jobtype==cteBORNOPPENHEIMER  || ct->jobtype==cteDLZSH || ct->jobtype==cteALZSH || ct->jobtype==cteDFSSH || ct->jobtype==cteJFSSH ||
            ct->jobtype==cteBCJFSSH || ct->jobtype==cteSCRDFSSH || ct->jobtype==cteCCFSSH || ct->jobtype==cteGFSH || ct->jobtype==cteDISH ||
            ct->jobtype==cteSCRDFSSH || ct->jobtype==cteTULLYLOC || ct->jobtype==ctePERSICOSFHOPPING)
        {
            /* calculate the QM energy
            //new jjk dholub may18
             * attention - the correct expression must be used for that!
             * THIS MAY BE WRONG FOR cteFERMI, WHERE WE MIX THE ADIABATIC STATES AND
             *     CONSTRUCT ARTIFICIAL ct->wf !
             */
         // double ct_energy = 0.0, ct_energy1 = 0.0, ct_energy2 = 0.0;
         // for (int i=0; i<ct->dim; i++)
         //     for (int m=0; m<ct->dim; m++)
         //     {
                    /* TB Hamiltonian */
         //         ct_energy1 += ct->hamiltonian[i][m] * (ct->wf[i] * ct->wf[m] + ct->wf[i + ct->dim] * ct->wf[m + ct->dim]);
                    /* Hubbard / gamma terms */
         //         ct_energy2 += 0.5 * ct->hubbard[i][m] * (SQR(ct->wf[i]) + SQR(ct->wf[i + ct->dim])) * (SQR(ct->wf[m]) + SQR(ct->wf[m + ct->dim]));
         //     }
         // ct_energy = ct_energy1 + ct_energy2;
            /* output both MM and QM energy */
            fprintf(f_ct_energy, "%8.2f", t*1000);

            if (ct->jobtype==cteDFSSH || ct->jobtype==cteDLZSH)
            {
                fprintf(f_ct_energy, "%12.5f ", ct->hamiltonian[ct->surface][ct->surface]*27.2114);
                for (int i=0; i < ct->dim; i++)
                    fprintf(f_ct_energy, "%12.5f", ct->hamiltonian[i][i]*27.2114);
            }
            else if (ct->jobtype==cteSCCDYNAMIC || ct->jobtype==cteNONSCCDYNAMIC || ct->jobtype==cteADNONSCC || ct->jobtype==cteCPFSCCDYNAMIC ||
                     ct->jobtype==cteDDBSCCDYNAMIC || ct->jobtype==cteNOMOVEMENT)
            {
                double ct_energy = 0.0;
                for (int i = 0; i < ct->dim; i++)
                    for (int m = 0; m < ct->dim; m++)
                        ct_energy += ct->hamiltonian[i][m] * (ct->wf[i] * ct->wf[m] + ct->wf[i + ct->dim] * ct->wf[m + ct->dim]);
                fprintf(f_ct_energy, "%12.5f", ct_energy*27.2114);
          
             // fprintf(f_ct_energy, "%12.5f %12.5f %12.5f %12.5f %12.5f",
             //         ct->hamiltonian[14][14]*27.2114, ct->hamiltonian[15][15]*27.2114, ct->hamiltonian[16][16]*27.2114,
             //         ct->hamiltonian[14][15]*27.2114, ct->hamiltonian[15][16]*27.2114);
            }
            else
            {
                fprintf(f_ct_energy, "%12.5f ", ct->ev_adiab[ct->surface]*27.2114);
                for (int i=0; i < ct->dim; i++)
                    fprintf(f_ct_energy, "%12.5f", ct->ev_adiab[i]*27.2114);
            }
            fprintf(f_ct_energy, "\n");
         // fprintf(f_ct_energy, "%10d %12.7f %12.7f %12.7f %12.7f %12.7f\n",
         //                      step, ct_energy1, ct_energy2, ct_energy, enerd->term[F_ETOT] * KJMOL_TO_HARTREE,
         //                      enerd->term[F_ETOT] * KJMOL_TO_HARTREE + ct_energy);
        }
#if GMX_MPI
        }
#endif

        /****************
         * END TRANSFER *
         ****************/

        /* Output stuff */
        if (MASTER(cr))
        {
            if (fplog && do_log && bDoExpanded)
            {
                /* only needed if doing expanded ensemble */
                PrintFreeEnergyInfoToFile(fplog, ir->fepvals, ir->expandedvals,
                                          ir->bSimTemp ? ir->simtempvals : nullptr,
                                          state_global->dfhist, state->fep_state, ir->nstlog, step);
            }
            if (bCalcEner)
            {
                energyOutput.addDataAtEnergyStep(
                        bDoDHDL, bCalcEnerStep, t, mdatoms->tmass, enerd, ir->fepvals,
                        ir->expandedvals, lastbox,
                        PTCouplingArrays{ state->boxv, state->nosehoover_xi, state->nosehoover_vxi,
                                          state->nhpres_xi, state->nhpres_vxi },
                        state->fep_state, shake_vir, force_vir, total_vir, pres, ekind, mu_tot, constr);
            }
            else
            {
                energyOutput.recordNonEnergyStep();
            }

            gmx_bool do_dr = do_per_step(step, ir->nstdisreout);
            gmx_bool do_or = do_per_step(step, ir->nstorireout);

            if (doSimulatedAnnealing)
            {
                gmx::EnergyOutput::printAnnealingTemperatures(do_log ? fplog : nullptr, groups,
                                                              &(ir->opts));
            }
            if (do_log || do_ene || do_dr || do_or)
            {
                energyOutput.printStepToEnergyFile(mdoutf_get_fp_ene(outf), do_ene, do_dr, do_or,
                                                   do_log ? fplog : nullptr, step, t,
                                                   fr->fcdata.get(), awh.get());
            }
            if (do_log && ir->bDoAwh && awh->hasFepLambdaDimension())
            {
                const bool isInitialOutput = false;
                printLambdaStateToLog(fplog, state->lambda, isInitialOutput);
            }

            if (ir->bPull)
            {
                pull_print_output(pull_work, step, t);
            }

            if (do_per_step(step, ir->nstlog))
            {
                if (fflush(fplog) != 0)
                {
                    gmx_fatal(FARGS, "Cannot flush logfile - maybe you are out of disk space?");
                }
            }
        }
        if (bDoExpanded)
        {
            /* Have to do this part _after_ outputting the logfile and the edr file */
            /* Gets written into the state at the beginning of next loop*/
            state->fep_state = lamnew;
        }
        else if (ir->bDoAwh && awh->needForeignEnergyDifferences(step))
        {
            state->fep_state = awh->fepLambdaState();
        }
        /* Print the remaining wall clock time for the run */
        if (isMasterSimMasterRank(ms, MASTER(cr)) && (do_verbose || gmx_got_usr_signal()) && !bPMETunePrinting)
        {
            if (shellfc)
            {
                fprintf(stderr, "\n");
            }
            print_time(stderr, walltime_accounting, step, ir, cr);
        }

        /* Ion/water position swapping.
         * Not done in last step since trajectory writing happens before this call
         * in the MD loop and exchanges would be lost anyway. */
        bNeedRepartition = FALSE;
        if ((ir->eSwapCoords != eswapNO) && (step > 0) && !bLastStep && do_per_step(step, ir->swap->nstswap))
        {
            bNeedRepartition =
                    do_swapcoords(cr, step, t, ir, swap, wcycle, as_rvec_array(state->x.data()),
                                  state->box, MASTER(cr) && mdrunOptions.verbose, bRerunMD);

            if (bNeedRepartition && DOMAINDECOMP(cr))
            {
                dd_collect_state(cr->dd, state, state_global);
            }
        }

        /* Replica exchange */
        bExchanged = FALSE;
        if (bDoReplEx)
        {
            bExchanged = replica_exchange(fplog, cr, ms, repl_ex, state_global, enerd, state, step, t);
        }

        if ((bExchanged || bNeedRepartition) && DOMAINDECOMP(cr))
        {
            dd_partition_system(fplog, mdlog, step, cr, TRUE, 1, state_global, *top_global, ir,
                                imdSession, pull_work, state, &f, mdAtoms, &top, fr, vsite, constr,
                                nrnb, wcycle, FALSE);
            shouldCheckNumberOfBondedInteractions = true;
            upd.setNumAtoms(state->natoms);
        }

        bFirstStep = FALSE;
        bInitStep  = FALSE;

        /* #######  SET VARIABLES FOR NEXT ITERATION IF THEY STILL NEED IT ###### */
        /* With all integrators, except VV, we need to retain the pressure
         * at the current step for coupling at the next step.
         */
        if ((state->flags & (1U << estPRES_PREV))
            && (bGStatEveryStep || (ir->nstpcouple > 0 && step % ir->nstpcouple == 0)))
        {
            /* Store the pressure in t_state for pressure coupling
             * at the next MD step.
             */
            copy_mat(pres, state->pres_prev);
        }

        /* #######  END SET VARIABLES FOR NEXT ITERATION ###### */

        if ((membed != nullptr) && (!bLastStep))
        {
            rescale_membed(step_rel, membed, as_rvec_array(state_global->x.data()));
        }

        cycles = wallcycle_stop(wcycle, ewcSTEP);
        if (DOMAINDECOMP(cr) && wcycle)
        {
            dd_cycles_add(cr->dd, cycles, ddCyclStep);
        }

        /* increase the MD step number */
        step++;
        step_rel++;

#if GMX_FAHCORE
        if (MASTER(cr))
        {
            fcReportProgress(ir->nsteps + ir->init_step, step);
        }
#endif

        resetHandler->resetCounters(step, step_rel, mdlog, fplog, cr, fr->nbv.get(), nrnb,
                                    fr->pmedata, pme_loadbal, wcycle, walltime_accounting);

        /* If bIMD is TRUE, the master updates the IMD energy record and sends positions to VMD client */
        imdSession->updateEnergyRecordAndSendPositionsAndEnergies(bInteractiveMDstep, step, bCalcEner);
    }
    /* End of main MD loop */

    /* Closing TNG files can include compressing data. Therefore it is good to do that
     * before stopping the time measurements. */
    mdoutf_tng_close(outf);

    /* Stop measuring walltime */
    walltime_accounting_end_time(walltime_accounting);

    if (!thisRankHasDuty(cr, DUTY_PME))
    {
        /* Tell the PME only node to finish */
        gmx_pme_send_finish(cr);
    }

    if (MASTER(cr))
    {
        if (ir->nstcalcenergy > 0)
        {
            energyOutput.printEnergyConservation(fplog, ir->simulation_part, EI_MD(ir->eI));

            gmx::EnergyOutput::printAnnealingTemperatures(fplog, groups, &(ir->opts));
            energyOutput.printAverages(fplog, groups);
        }
    }
    done_mdoutf(outf);

    if (bPMETune)
    {
        pme_loadbal_done(pme_loadbal, fplog, mdlog, fr->nbv->useGpu());
    }

    done_shellfc(fplog, shellfc, step_rel);

    if (useReplicaExchange && MASTER(cr))
    {
        print_replica_exchange_statistics(fplog, repl_ex);
    }

    walltime_accounting_set_nsteps_done(walltime_accounting, step_rel);

    global_stat_destroy(gstat);

    /************
     * TRANSFER *
     ************/

#if GMX_MPI
    if (ct_mpi_rank == 0) {
#endif
    if (f_ct_energy) fclose(f_ct_energy);
    if (f_ct_esp) fclose(f_ct_esp);
    if (f_ct_shift) fclose(f_ct_shift);
    if (f_tb_hamiltonian) fclose(f_tb_hamiltonian);
    if (f_tb_hamiltonian_ml) fclose(f_tb_hamiltonian_ml);
    if (f_tb_hamiltonian_hub) fclose(f_tb_hamiltonian_hub);
    if (f_tb_hubbard) fclose(f_tb_hubbard);
    if (f_tb_occupation) fclose(f_tb_occupation);
    if (f_tb_diapop) fclose(f_tb_diapop);
    if (f_tb_adiapop) fclose(f_tb_adiapop);
    if (f_tb_occupation_adiab) fclose(f_tb_occupation_adiab);
    if (f_ct_adiabatic) fclose(f_ct_adiabatic);
    if (f_ct_exp_adiab) fclose(f_ct_exp_adiab);
    if (f_ct_surfacehopping) fclose(f_ct_surfacehopping);
    if (f_ct_orbital) fclose(f_ct_orbital);
    if (f_ct_current) fclose(f_ct_current);
    if (f_tb_com) fclose(f_tb_com);  /* dholub */
    if (f_tb_spread) fclose(f_tb_spread);
    if (f_tb_diag_force) fclose(f_tb_diag_force);
    if (f_tb_offdiag_force) fclose(f_tb_offdiag_force);
#if GMX_MPI
    }
#endif

    /****************
     * END TRANSFER *
     ****************/
}
