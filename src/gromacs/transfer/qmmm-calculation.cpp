
#include "gmxpre.h"

#include "qmmm.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmath>

#include <algorithm>

#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/ewald/pme.h"
#include "gromacs/ewald/pme_internal.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/nblist.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/mtop_lookup.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

#define NM2BOHR           (1 / BOHR2NM)
#define HARTREE_TO_EV     (27.211396132)
// #define AU_OF_ESP_TO_VOLT (14.400) -- this is off by an angstrom/bohr factor! use HARTREE_TO_EV instead!
#define HARTREE2KJMOL     (HARTREE2KJ * AVOGADRO)
#define KJMOL2HARTREE     (1 / HARTREE2KJMOL)

#define QMMM_SWITCH       (0.05) // length of the additional switching region beyond cutoff

#define SQR(x) ((x)*(x))
#define CUB(x) ((x)*(x)*(x))
#define QRT(x) ((x)*(x)*(x)*(x))
#define QUI(x) ((x)*(x)*(x)*(x)*(x))
#define HEX(x) ((x)*(x)*(x)*(x)*(x)*(x))
#define OCT(x) ((x)*(x)*(x)*(x)*(x)*(x)*(x)*(x))
#define CHOOSE2(x) ((x)*((x)+1)/2)

#define gmx_erfc(x) (std::erfc(x))
#define gmx_erf(x)  (std::erf(x))

#include<time.h>

/****************************************
 ******   AUXILIARY ROUTINES   **********
 ****************************************/

// adopted from src/gmxlib/pbc.c
static inline void pbc_dx_qmmm(matrix box, const rvec x1, const rvec x2, rvec dx)
{
    int i;

    for(i=0; i<DIM; i++) {
        dx[i] = x1[i] - x2[i];
        if (box != nullptr)
        {
            real length = box[i][i];
            while (dx[i] > length / 2.) {
                dx[i] -= length;
            }
            while (dx[i] < - length / 2.) {
                dx[i] += length;
            }
        }
    }
}

static inline real pbc_dist_qmmm(matrix box, const rvec x1, const rvec x2)
{
    rvec dx;

    for(int i=0; i<DIM; i++) {
        dx[i] = x1[i] - x2[i];

        if (box != nullptr)
        {
            real length = box[i][i];
            while (dx[i] > length / 2.) {
                dx[i] -= length;
            }
            while (dx[i] < - length / 2.) {
                dx[i] += length;
            }
        }
    }
    return norm(dx);
}

/****************************************
 ******   PRE-SCC CALCULATION  **********
 ****************************************/

void QMMM_rec_transfer::calculate_SR_QM_MM(int variant,
                                           real *pot)
{
  real rcoul = qm->rcoulomb;
  real ewaldcoeff_q = qm->ewaldcoeff_q;

  switch (variant) {

  case eqmmmVACUO: // no QM/MM
  {
      for (int j=0; j<qm->nrQMatoms; j++)
      {
          pot[j] = 0.;
      }
      break;
  }

  case eqmmmSWITCH: // QM/MM with switched cut-off
  {
    real r_1 = rcoul;
    real r_d = QMMM_SWITCH;
    real r_c = r_1 + r_d;
  //printf("r_1 = %f, r_d = %f, r_c = %f\n", r_1, r_d, r_c);
    real big_a =   (5 * r_c - 2 * r_1) / (CUB(r_c) * SQR(r_d));
    real big_b = - (4 * r_c - 2 * r_1) / (CUB(r_c) * CUB(r_d));
    real big_c = - 1 / r_c - big_a / 3 * CUB(r_d) - big_b / 4 * QRT(r_d);
    for (int j=0; j<qm->nrQMatoms; j++) { // do it for every QM atom
      pot[j] = 0.;
      int under_r1=0, under_rc=0;
      // add potential from MM atoms
      for (int k=0; k<mm->nrMMatoms; k++) {
        real r = pbc_dist_qmmm(qm->box, qm->xQM[j], mm->xMM[k]);
        if (r < 0.001) { // this may occur on the first step of simulation for link atom(s)
          printf("QM--MM exploding for QM=%d, MM=%d. MM charge is %f\n", j+1, k+1, mm->MMcharges[k]);
          continue;
        }
        if (r < r_1) {
          pot[j] += mm->MMcharges[k] * (1. / r + big_c);
          under_r1++;
          continue;
        }
        if (r < r_c) {
          pot[j] += mm->MMcharges[k] * ( 1. / r + big_a / 3. * CUB(r - r_1) + big_b / 4. * QRT(r - r_1) + big_c);
          under_rc++;
          continue;
        }
      } // for k
    //printf("ATOM %d: %8.5f %4d %4d\n", j+1, pot[j], under_r1, under_rc);
    } // for j
    break;
  }

  case eqmmmRFIELD: // reaction field with epsilon=infinity
  {
    real r_c = rcoul;
    real big_c = 3. / (2. * r_c);
    for (int j=0; j<qm->nrQMatoms; j++) { // do it for every QM atom
      pot[j] = 0.;
      // add potential from MM atoms
      for (int k=0; k<mm->nrMMatoms; k++) {
        real r = pbc_dist_qmmm(qm->box, qm->xQM[j], mm->xMM[k]);
        if (r < 0.001) { // this may occur on the first step of simulation for link atom(s)
          printf("QM--MM exploding for QM=%d, MM=%d. MM charge is %f\n", j+1, k+1, mm->MMcharges[k]);
          continue;
        }
        if (r < r_c) {
          pot[j] += mm->MMcharges[k] * ( 1. / r + SQR(r) / 2. / CUB(r_c) - big_c);
          continue;
        }
      } // for k
    } // for j
    break;
  }

  case eqmmmSHIFT: // shifted cut-off, in a similar spirit as reaction field with epsilon=infinity
  {
    real r_c = rcoul;
    real big_c = 3. / SQR(r_c);
    real big_k = 3. / r_c;
    for (int j=0; j<qm->nrQMatoms; j++) { // do it for every QM atom
      pot[j] = 0.;
      // add potential from MM atoms
      for (int k=0; k<mm->nrMMatoms; k++) {
        real r = pbc_dist_qmmm(qm->box, qm->xQM[j], mm->xMM[k]);
        if (r < 0.001) { // this may occur on the first step of simulation for link atom(s)
          printf("QM--MM exploding for QM=%d, MM=%d. MM charge is %f\n", j+1, k+1, mm->MMcharges[k]);
          continue;
        }
        if (r < r_c) {
          pot[j] += mm->MMcharges[k] * ( 1. / r - SQR(r) / CUB(r_c) + big_c * r - big_k);
          continue;
        }
      } // for k
    } // for j
    break;
  }

  case eqmmmPME: // QM/MM PME preparation -- short-range QM--MM component (real-space)
  {
    // using fr->ewaldcoeff_q, which has the dimension of 1/distance
    for (int j=0; j<qm->nrQMatoms; j++) { // do it for every QM atom
      pot[j] = 0.;
      // add potential from MM atoms
      for (int k=0; k<mm->nrMMatoms; k++) {
        real r = pbc_dist_qmmm(qm->box, qm->xQM[j], mm->xMM[k]);
        if (r < 0.001) { // this may occur on the first step of simulation for link atom(s)
          printf("QM/MM PME QM--MM short range exploding for QM=%d, MM=%d. MM charge is %f\n", j+1, k+1, mm->MMcharges[k]);
        } else {
          pot[j] += mm->MMcharges[k] / r * gmx_erfc(ewaldcoeff_q * r);
        }
      } // for k
    } // for j
    break;
  }

  default: // it should never get this far
    ;
  } // switch variant

  /* Convert the result to atomic units. */
  for (int j=0; j<qm->nrQMatoms; j++)
  {
      pot[j] /= NM2BOHR;
   // printf("SR QM/MM POT in a.u.: %d %8.5f\n", j+1, pot[j]);
  }
} // calculate_SR_QM_MM

/* Calculate the effect of environment with PME,
 * EXCLUDING the periodic images of QM charges at this stage!
 * That contribution will be added within the SCC cycle.
 */
void QMMM_rec_transfer::calculate_LR_QM_MM(const t_commrec *cr,
                                           t_nrnb *nrnb,
                                           gmx_wallcycle_t wcycle,
                                           struct gmx_pme_t *pmedata,
                                           real *pot)
{
  QMMM_PME_transfer& pme_full   = pme[0];
  QMMM_PME_transfer& pme_qmonly = pme[1];
  const int   n  = qm->nrQMatoms;
  const int   ne = mm->nrMMatoms_full;
//const int   ntot = n + ne;

  /* copy the data into PME structures */
  for (int j=0; j<n; j++)
  {
      /* QM atoms -- no charges! */
      pme_full.x[j][XX] = qm->xQM[j][XX];
      pme_full.x[j][YY] = qm->xQM[j][YY];
      pme_full.x[j][ZZ] = qm->xQM[j][ZZ];
      pme_full.q[j]     = 0.;
   // printf("QM %5d %8.5f %8.5f %8.5f\n", j+1, pme_full.x[j][XX], pme_full.x[j][YY], pme_full.x[j][ZZ]);

      pme_qmonly.x[j][XX] = qm->xQM[j][XX];
      pme_qmonly.x[j][YY] = qm->xQM[j][YY];
      pme_qmonly.x[j][ZZ] = qm->xQM[j][ZZ];
      pme_qmonly.q[j]     = 0.;
  }
  for (int j=0; j<ne; j++)
  {
      /* MM atoms -- with charges */
      pme_full.x[n + j][XX] = mm->xMM_full[j][XX];
      pme_full.x[n + j][YY] = mm->xMM_full[j][YY];
      pme_full.x[n + j][ZZ] = mm->xMM_full[j][ZZ];
      pme_full.q[n + j]     = mm->MMcharges_full[j];
   // printf("MM %5d %8.5f %8.5f %8.5f %8.5f\n", j+1, pme->x[n+j][XX], pme->x[n+j][YY], pme->x[n+j][ZZ], pme->q[n + j]);
  }
  
//static struct timespec time1, time2;
//clock_gettime(CLOCK_MONOTONIC, &time1);
  // init_nrnb(pme_full.nrnb); // TODO change to something?
  gmx::StepWorkload stepWork;
  stepWork.computePotentials = true;
  PaddedVector<gmx::RVec> emptyVec;
  gmx_pme_do(pmedata, gmx::makeArrayRef(pme_full.x), gmx::makeArrayRef(emptyVec), pme_full.q.data(), pme_full.q.data(),
             nullptr, nullptr, nullptr, nullptr, qm->box, cr, 0, 0, //pme_full.nrnb->get(),
             nrnb, wcycle, pme_full.vir, pme_full.vir, nullptr, nullptr, 0., 0., nullptr, nullptr,
             stepWork, TRUE, FALSE, n, pme_full.pot);
//clock_gettime(CLOCK_MONOTONIC, &time2);
//print_time_difference("PMETIME 1 ", time1, time2);

  /* Save the potential */
  for (int j=0; j<qm->nrQMatoms; j++)
  {
      pot[j] = pme_full.pot[j] * KJMOL2HARTREE; // conversion OK
  }

  if (pme_full.surf_corr_pme)
  {
      /* optionally evaluate the PME surface correction term.
       * ATTENTION: modified update_QMMM_coord() (qmmm.cpp) is needed here!
       */
       // sum_j q_j vec(x_j)
       rvec qx, sum_qx;
       clear_rvec(sum_qx);
    // rvec subsum_qx;
	   for (int j=0; j<ne; j++) {
           svmul(mm->MMcharges_full[j], mm->xMM_full[j], qx);
           rvec_inc(sum_qx, qx);
        // if (j%3==0) {
        //   printf("MOL %4d DIPOLE %5.1f %5.1f %5.1f\n", j/3,
        //     subsum_qx[XX]/0.020819434, subsum_qx[YY]/0.020819434, subsum_qx[ZZ]/0.020819434);
        //   clear_rvec(subsum_qx);
        // }
        // rvec_inc(subsum_qx, qx);
	   }
	   // contribution to the potential
       std::vector<real> pot_add(n);
	   real vol = qm->box[XX][XX] * qm->box[YY][YY] * qm->box[ZZ][ZZ];
	   for (int j=0; j<n; j++) {
	       pot_add[j] = 4. * M_PI / 3. / vol / pme_full.epsilon_r * iprod(qm->xQM[j], sum_qx) / NM2BOHR;
       }
    // printf("VOL = %8.5f, EPS_R = %8.5f, DIP = %9.5f %9.5f %9.5f\n",
    //         vol, pme_full.epsilon_r, sum_qx[XX], sum_qx[YY], sum_qx[ZZ]);
	   for (int j=0; j<n; j++) {
    //     printf("POT LR [%d] = %8.5f + %8.5f\n", j+1, pot[j], pot_add[j]);
           pot[j] += pot_add[j];
       }
  }
  else
  {
	// for (int j=0; j<n; j++) {
    //     printf("POT LR [%d] = %8.5f\n", j+1, pot[j]);
    // }
  }

//for (int j=0; j<n; j++) {
//    printf("POT LR [%d] = %8.5f\n", j+1, pot[j]);
//}
} // calculate_LR_QM_MM

/****************************************
 *********  IN-SCC CALCULATION  *********
 ****************************************/

void QMMM_rec_transfer::calculate_complete_QM_QM(const t_commrec*  cr,
                                                 t_nrnb*           nrnb,
                                                 gmx_wallcycle_t   wcycle,
                                                 struct gmx_pme_t* pmedata,
                                                 real*             pot)
{
  QMMM_PME_transfer& pme_qmonly   = pme[1];
  int                n            = qm->nrQMatoms;
//real               rcoul        = qm_.rcoulomb;
  real               ewaldcoeff_q = qm->ewaldcoeff_q;

///* (re)allocate */
//srenew(pme->x, qm->nrQMatoms);
//srenew(pme->q, qm->nrQMatoms);

  /* copy the data into PME structures */
  for (int j=0; j<n; j++)
  {
      /* QM atoms with charges */
      pme_qmonly.x[j][XX] = qm->xQM[j][XX];
      pme_qmonly.x[j][YY] = qm->xQM[j][YY];
      pme_qmonly.x[j][ZZ] = qm->xQM[j][ZZ];
      /* Attenuate the periodic images of the QM zone
       * with the same scaling factor that is applied for the MM atoms.
       */
      pme_qmonly.q[j]     = qm->QMcharges[j] * mm->scalefactor;
  }
  
//static struct timespec time1, time2;
//clock_gettime(CLOCK_MONOTONIC, &time1);
  // init_nrnb(pme_qmonly.nrnb); // TODO change to something?
  gmx::StepWorkload stepWork;
  stepWork.computePotentials = true;
  PaddedVector<gmx::RVec> emptyVec;
  int oldNumAtoms = pmedata->atc[0].numAtoms(); // need to resize PME arrays to the number of QM atoms
  pmedata->atc[0].setNumAtoms(n);
  gmx_pme_do(pmedata, pme_qmonly.x, gmx::makeArrayRef(emptyVec), pme_qmonly.q.data(), pme_qmonly.q.data(),
             nullptr, nullptr, nullptr, nullptr, qm->box, cr, 0, 0, // pme_qmonly.nrnb->get(),
             nrnb, wcycle, pme_qmonly.vir, pme_qmonly.vir, nullptr, nullptr, 0., 0., nullptr, nullptr,
             stepWork, TRUE, FALSE, n, pme_qmonly.pot);
  pmedata->atc[0].setNumAtoms(oldNumAtoms); // resize back
//clock_gettime(CLOCK_MONOTONIC, &time2);
//print_time_difference("PMETIME 2 ", time1, time2);

//for (int j=0; j<n; j++) {
//  printf("POT QM QM [%3d] = %12.7f\n", j, pme_qmonly.pot[j]);
//}
  
  /* short-range corrections */
  std::vector<real> pot_corr(n);
  for (int j=0; j<n; j++)
  {
      /* exclude the interaction of atom j with its own charge density */
      pot_corr[j] = - 2. * ewaldcoeff_q * pme_qmonly.q[j] / sqrt(M_PI);
      /* exclude the interactions with the other QM atoms */
      for (int k=0; k<n; k++)
      {
          if (j != k)
          {
              real r = pbc_dist_qmmm(nullptr, pme_qmonly.x[j], pme_qmonly.x[k]);
              pot_corr[j] -= pme_qmonly.q[k] * gmx_erf(ewaldcoeff_q * r) / r;
          }
      }
  }
      
  std::vector<real> pot_surf(n);
  if (pme_qmonly.surf_corr_pme)
  {
      /* optionally evaluate the PME surface correction term.
       * ATTENTION: modified update_QMMM_coord() (qmmm.c) is needed here!
       */
       // sum_j q_j vec(x_j)
       rvec qx, sum_qx;
       clear_rvec(sum_qx);
	   for (int j=0; j<qm->nrQMatoms; j++) {
           svmul(pme_qmonly.q[j], pme_qmonly.x[j], qx);
           rvec_inc(sum_qx, qx);
	   }
	   // contribution to the potential
	   real vol = qm->box[XX][XX] * qm->box[YY][YY] * qm->box[ZZ][ZZ];
	   for (int j=0; j<n; j++) {
	       pot_surf[j] = 4. * M_PI / 3. / vol / pme_qmonly.epsilon_r * iprod(qm->xQM[j], sum_qx);
       }
  }
  else
  {
	   for (int j=0; j<n; j++) {
           pot_surf[j] = 0.;
	   }
  }

  /* return the potential on QM atoms */
  for (int j=0; j<n; j++)
  {
      pot[j] = pme_qmonly.pot[j] * KJMOL2HARTREE + pot_corr[j] * BOHR2NM + pot_surf[j] * BOHR2NM;
   // printf("pot_qm_in_scc[%d] = %12.8f\n", j+1, pot[j]);
   // printf("Ewald atom %d charge %6.3f potential %8.5f (correction %8.5f surfterm %8.5f)\n",
   //         j+1, pme_qmonly.q[j],      pot[j],                 pot_corr[j] * BOHR2NM, pot_surf[j] * BOHR2NM);
  }

  // also, save the potential in the QMMM_QMrec structure
  for (int j=0; j<n; j++)
  {
      qm->pot_qmqm_set(j, static_cast<double>(pot[j]) * HARTREE_TO_EV); // in volt units
  }
} // calculate_complete_QM_QM

