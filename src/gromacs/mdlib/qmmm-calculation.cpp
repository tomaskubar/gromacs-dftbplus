
#include "gmxpre.h"

#include "qmmm.h"

#include "config.h"

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
#define NM_TO_BOHR        (18.897259886)
#define BOHR_TO_NM        (1 / NM_TO_BOHR)
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

void print_time_difference(const char s[],
                           struct timespec start,
                           struct timespec end);

void print_time_difference(const char s[],
                           struct timespec start,
                           struct timespec end)
{
  //int sec, nsec;
  long long value = 1000000000ll * (static_cast<long long>(end.tv_sec)
                                  - static_cast<long long>(start.tv_sec))
                    + static_cast<long long>(end.tv_nsec - start.tv_nsec);
  printf("%s %12lld\n", s, value);
}  

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

static inline rvec* pbc_dist_vec_qmmm(matrix box, const rvec x1, const rvec x2)
{
    rvec* dx;

    for(int i=0; i<DIM; i++) {
        *dx[i] = x1[i] - x2[i];

        if (box != nullptr)
        {
            real length = box[i][i];
            while (*dx[i] > length / 2.) {
                *dx[i] -= length;
            }
            while (*dx[i] < - length / 2.) {
                *dx[i] += length;
            }
        }
    }
    return dx;
}


/****************************************
 ******   PRE-SCC CALCULATION  **********
 ****************************************/

void QMMM_rec::calculate_SR_QM_MM(int variant,
                                  real *pot)
{
  /* as the last argument is expected: fr->ewaldcoeff_q */
  QMMM_QMrec& qm_ = qm[0];
  QMMM_MMrec& mm_ = mm[0];
  real rcoul = qm_.rcoulomb;
  real ewaldcoeff_q = qm_.ewaldcoeff_q;

  switch (variant) {

  case eqmmmVACUO: // no QM/MM
  {
      for (int j=0; j<qm_.nrQMatoms; j++)
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
    for (int j=0; j<qm_.nrQMatoms; j++) { // do it for every QM atom
      pot[j] = 0.;
      int under_r1=0, under_rc=0;
      // add potential from MM atoms
      for (int k=0; k<mm_.nrMMatoms; k++) {
        real r = pbc_dist_qmmm(qm_.box, qm_.xQM[j], mm_.xMM[k]);
        if (r < 0.001) { // this may occur on the first step of simulation for link atom(s)
          printf("QM--MM exploding for QM=%d, MM=%d. MM charge is %f\n", j+1, k+1, mm_.MMcharges[k]);
          continue;
        }
        if (r < r_1) {
          pot[j] += mm_.MMcharges[k] * (1. / r + big_c);
          under_r1++;
          continue;
        }
        if (r < r_c) {
          pot[j] += mm_.MMcharges[k] * ( 1. / r + big_a / 3. * CUB(r - r_1) + big_b / 4. * QRT(r - r_1) + big_c);
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
    for (int j=0; j<qm_.nrQMatoms; j++) { // do it for every QM atom
      pot[j] = 0.;
      // add potential from MM atoms
      for (int k=0; k<mm_.nrMMatoms; k++) {
        real r = pbc_dist_qmmm(qm_.box, qm_.xQM[j], mm_.xMM[k]);
        if (r < 0.001) { // this may occur on the first step of simulation for link atom(s)
          printf("QM--MM exploding for QM=%d, MM=%d. MM charge is %f\n", j+1, k+1, mm_.MMcharges[k]);
          continue;
        }
        if (r < r_c) {
          pot[j] += mm_.MMcharges[k] * ( 1. / r + SQR(r) / 2. / CUB(r_c) - big_c);
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
    for (int j=0; j<qm_.nrQMatoms; j++) { // do it for every QM atom
      pot[j] = 0.;
      // add potential from MM atoms
      for (int k=0; k<mm_.nrMMatoms; k++) {
        real r = pbc_dist_qmmm(qm_.box, qm_.xQM[j], mm_.xMM[k]);
        if (r < 0.001) { // this may occur on the first step of simulation for link atom(s)
          printf("QM--MM exploding for QM=%d, MM=%d. MM charge is %f\n", j+1, k+1, mm_.MMcharges[k]);
          continue;
        }
        if (r < r_c) {
          pot[j] += mm_.MMcharges[k] * ( 1. / r - SQR(r) / CUB(r_c) + big_c * r - big_k);
          continue;
        }
      } // for k
    } // for j
    break;
  }

  case eqmmmPME: // QM/MM PME preparation -- short-range QM--MM component (real-space)
  {
    // using fr->ewaldcoeff_q, which has the dimension of 1/distance
    for (int j=0; j<qm_.nrQMatoms; j++) { // do it for every QM atom
      pot[j] = 0.;
      // add potential from MM atoms
      for (int k=0; k<mm_.nrMMatoms; k++) {
        real r = pbc_dist_qmmm(qm_.box, qm_.xQM[j], mm_.xMM[k]);
        if (r < 0.001) { // this may occur on the first step of simulation for link atom(s)
          printf("QM/MM PME QM--MM short range exploding for QM=%d, MM=%d. MM charge is %f\n", j+1, k+1, mm_.MMcharges[k]);
        } else {
          pot[j] += mm_.MMcharges[k] / r * gmx_erfc(ewaldcoeff_q * r);
        }
      } // for k
    } // for j
    break;
  }

  case eqmmmNN:
  {
    for (int j=0; j<qm_.nrQMatoms; j++) { // do it for every QM atom
      pot[j] = 0.;
      // add potential from MM atoms
      for (int k=0; k<mm_.nrMMatoms; k++) {
        real r = pbc_dist_qmmm(qm_.box, qm_.xQM[j], mm_.xMM[k]);
        if (r < 0.001) { // this may occur on the first step of simulation for link atom(s)
          printf("QM--MM exploding for QM=%d, MM=%d. MM charge is %f\n", j+1, k+1, mm_.MMcharges[k]);
          continue;
        }
        else {
          pot[j] += mm_.MMcharges[k] / r;
        }
      } // for k
    } // for j
    break;
  }

  default: // it should never get this far
    ;
  } // switch variant

  /* Convert the result to atomic units. */
  for (int j=0; j<qm_.nrQMatoms; j++)
  {
      pot[j] /= NM2BOHR;
   // printf("SR QM/MM POT in a.u.: %d %8.5f\n", j+1, pot[j]);
  }
} // calculate_SR_QM_MM

/* Calculate the effect of environment with PME,
 * EXCLUDING the periodic images of QM charges at this stage!
 * That contribution will be added within the SCC cycle.
 */

void QMMM_rec::calculate_SR_NN_MM(int variant,
                                  rvec *dpot,
                                  real *pot)
{
  /* as the last argument is expected: fr->ewaldcoeff_q */
  QMMM_QMrec& qm_ = qm[0];
  QMMM_MMrec& mm_ = mm[0];
  real rcoul = qm_.rcoulomb;
  real ewaldcoeff_q = qm_.ewaldcoeff_q;
  if (variant!=5) {
      printf("Wrong variant for ML/MM");
      exit(-1);
  }
  else {
      for (int j=0; j<qm_.nrQMatoms; j++) { // do it for every QM atom
      pot[j] = 0.;
      for (int l=0; l<3; l++) {
        dpot[j][l] = 0.;
      }
      // add potential from MM atoms
      for (int k=0; k<mm_.nrMMatoms; k++) {
        rvec dr;
        //real r = pbc_dist_qmmm(qm_.box, qm_.xQM[j], mm_.xMM[k]);
        pbc_dx_qmmm(qm_.box, qm_.xQM[j], mm_.xMM[k], dr);
        real r = norm(dr);
        if (r < 0.001) { // this may occur on the first step of simulation for link atom(s)
          printf("QM--MM exploding for QM=%d, MM=%d. MM charge is %f\n", j+1, k+1, mm_.MMcharges[k]);
          continue;
        }
        else {
          pot[j] += mm_.MMcharges[k] / r;
          for (int l=0; l<3; l++) {
              dpot[j][l] -= dr[l]/(r*r);
          }
        }
      } // for k
    } // for j
  }
/* Convert the result to atomic units. */
  for (int j=0; j<qm_.nrQMatoms; j++)
  {
      pot[j] /= NM2BOHR;
      for (int l=0; l<3; l++) {
          dpot[j][l] /= NM2BOHR;
      }
   // printf("SR QM/MM POT in a.u.: %d %8.5f\n", j+1, pot[j]);
  }
} // calculate_SR_QM_MM


void QMMM_rec::calculate_LR_QM_MM(const t_commrec *cr,
                                  t_nrnb *nrnb,
                                  gmx_wallcycle_t wcycle,
                                  struct gmx_pme_t *pmedata,
                                  real *pot)
{
  QMMM_QMrec& qm_      = qm[0];
  QMMM_MMrec& mm_      = mm[0];
  QMMM_PME& pme_full   = pme[0];
  QMMM_PME& pme_qmonly = pme[1];
//QMMM_PME *pme_full   = qr->pme_full;
//QMMM_PME *pme_qmonly = qr->pme_qmonly;
  const int   n  = qm_.nrQMatoms;
  const int   ne = mm_.nrMMatoms_full;
//const int   ntot = n + ne;

  /* copy the data into PME structures */
  for (int j=0; j<n; j++)
  {
      /* QM atoms -- no charges! */
      pme_full.x[j][XX] = qm_.xQM[j][XX];
      pme_full.x[j][YY] = qm_.xQM[j][YY];
      pme_full.x[j][ZZ] = qm_.xQM[j][ZZ];
      pme_full.q[j]     = 0.;
   // printf("QM %5d %8.5f %8.5f %8.5f\n", j+1, pme_full.x[j][XX], pme_full.x[j][YY], pme_full.x[j][ZZ]);

      pme_qmonly.x[j][XX] = qm_.xQM[j][XX];
      pme_qmonly.x[j][YY] = qm_.xQM[j][YY];
      pme_qmonly.x[j][ZZ] = qm_.xQM[j][ZZ];
      pme_qmonly.q[j]     = 0.;
  }
  for (int j=0; j<ne; j++)
  {
      /* MM atoms -- with charges */
      pme_full.x[n + j][XX] = mm_.xMM_full[j][XX];
      pme_full.x[n + j][YY] = mm_.xMM_full[j][YY];
      pme_full.x[n + j][ZZ] = mm_.xMM_full[j][ZZ];
      pme_full.q[n + j]     = mm_.MMcharges_full[j];
   // printf("MM %5d %8.5f %8.5f %8.5f %8.5f\n", j+1, pme->x[n+j][XX], pme->x[n+j][YY], pme->x[n+j][ZZ], pme->q[n + j]);
  }
  
//static struct timespec time1, time2;
//clock_gettime(CLOCK_MONOTONIC, &time1);
  // init_nrnb(pme_full.nrnb); // TODO change to something?
  gmx::StepWorkload stepWork;
  stepWork.computePotentials = true;
  PaddedVector<gmx::RVec> emptyVec;
  gmx_pme_do(pmedata, gmx::makeArrayRef(pme_full.x), gmx::makeArrayRef(emptyVec), pme_full.q.data(), pme_full.q.data(),
             nullptr, nullptr, nullptr, nullptr, qm_.box, cr, 0, 0, //pme_full.nrnb->get(),
             nrnb, wcycle, pme_full.vir, pme_full.vir, nullptr, nullptr, 0., 0., nullptr, nullptr,
             stepWork, TRUE, FALSE, n, pme_full.pot);
//clock_gettime(CLOCK_MONOTONIC, &time2);
//print_time_difference("PMETIME 1 ", time1, time2);

  /* Save the potential */
  for (int j=0; j<qm_.nrQMatoms; j++)
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
           svmul(mm_.MMcharges_full[j], mm_.xMM_full[j], qx);
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
	   real vol = qm_.box[XX][XX] * qm_.box[YY][YY] * qm_.box[ZZ][ZZ];
	   for (int j=0; j<n; j++) {
	       pot_add[j] = 4. * M_PI / 3. / vol / pme_full.epsilon_r * iprod(qm_.xQM[j], sum_qx) / NM2BOHR;
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

/* Calculate the potential induced by the periodic images of QM charges.
 * This needs to be performed in every SCC iteration.
 */
/*
void calculate_complete_QM_QM_ewald(t_QMMMrec *qr,
                              const t_commrec *cr,
                              gmx_wallcycle_t wcycle,
                              struct gmx_pme_t *pmedata,
                              real *pot)
{
  // as the last arguments are expected: fr->ewaldcoeff_q and fr->pmedata
  t_QMrec    *qm           = qr->qm[0];
  t_QMMM_PME *pme          = qr->pme;
  int         n            = qm->nrQMatoms;
//real        rcoul        = qm->rcoulomb;
  real        ewaldcoeff_q = qm->ewaldcoeff_q;

  // copy the data into PME structures
  for (int j=0; j<n; j++)
  {
      // QM atoms with charges
      pme->x[j][XX] = qm->xQM[j][XX];
      pme->x[j][YY] = qm->xQM[j][YY];
      pme->x[j][ZZ] = qm->xQM[j][ZZ];
      // Attenuate the periodic images of the QM zone
      // with the same scaling factor
      // that is applied for the MM atoms.
      pme->q[j]     = qm->QMcharges[j] * qr->mm->scalefactor;
  }
  
  static struct timespec time1, time2;
  clock_gettime(CLOCK_MONOTONIC, &time1);
  init_nrnb(pme->nrnb);

  for (int j=0; j<n; j++) {
    pme->pot[j] = 0.;
  }

  rvec kvec;
  for (int kxi=-pmedata->nkx; kxi<=pmedata->nkx; kxi++) {
    kvec[XX] = (real) kxi * 2. * M_PI / qm->box[XX][XX];
    for (int kyi=-pmedata->nky; kyi<=pmedata->nky; kyi++) {
      kvec[YY] = (real) kyi * 2. * M_PI / qm->box[YY][YY];
      for (int kzi=-pmedata->nkz; kzi<=pmedata->nkz; kzi++) {
        kvec[ZZ] = (real) kzi * 2. * M_PI / qm->box[ZZ][ZZ];

        if (kxi != 0 || kyi != 0 || kzi != 0) {
          real factor = exp(-norm2(kvec) / 4. / SQR(pmedata->ewaldcoeff_q)) / norm2(kvec);
        //printf("KX %2d KY %2d KZ %2d FACTOR %12.7f\n", kxi, kyi, kzi, factor);
          for (int j=0; j<n; j++) {
            for (int k=0; k<n; k++) {
              rvec bond;
              rvec_sub(pme->x[k], pme->x[j], bond);
              real addend = factor * pme->q[k] * cos((double) iprod(kvec, bond));
              pme->pot[j] += addend;
            //printf("J %2d K %2d X %8.5f Y %8.5f Z %8.5f ADDEND %12.7f\n",
            //  j, k, bond[XX], bond[YY], bond[ZZ], addend);
            }
          }
        }

      }
    }
  }

  for (int j=0; j<n; j++) {
    real factor = 4. * M_PI / (qm->box[XX][XX] * qm->box[YY][YY] * qm->box[ZZ][ZZ])
                * ONE_4PI_EPS0 / pmedata->epsilon_r;
  //printf("final factor = %12.7f\n", factor);
    pme->pot[j] *= factor;
  }

  printf("NKX %d NKY %d NKZ %d\n", pmedata->nkx, pmedata->nky, pmedata->nkz);
  for (int j=0; j<n; j++) {
    printf("POT QM QM [%3d] = %12.7f\n", j, pme->pot[j]);
  }

  clock_gettime(CLOCK_MONOTONIC, &time2);
  print_time_difference("EWATIME 2 ", time1, time2);
  
  // short-range corrections
  std::vector<real> pot_corr(n);
  for (int j=0; j<n; j++)
  {
      // exclude the interaction of atom j with its own charge density
      pot_corr[j] = - 2. * ewaldcoeff_q * pme->q[j] / sqrt(M_PI);
      // exclude the interactions with the other QM atoms
      for (int k=0; k<n; k++)
      {
          if (j != k)
          {
              real r = pbc_dist_qmmm(nullptr, pme->x[j], pme->x[k]);
              pot_corr[j] -= pme->q[k] * gmx_erf(ewaldcoeff_q * r) / r;
          }
      }
  }
      
  std::vector<real> pot_surf(n);
  if (pme->surf_corr_pme)
  {
      // optionally evaluate the PME surface correction term.
      // ATTENTION: modified update_QMMM_coord() (qmmm.c) is needed here!
       // sum_j q_j vec(x_j)
       rvec qx, sum_qx;
       clear_rvec(sum_qx);
	   for (int j=0; j<qm->nrQMatoms; j++) {
           svmul(pme->q[j], pme->x[j], qx);
           rvec_inc(sum_qx, qx);
	   }
	   // contribution to the potential
	   real vol = qm->box[XX][XX] * qm->box[YY][YY] * qm->box[ZZ][ZZ];
	   for (int j=0; j<n; j++) {
	       pot_surf[j] = 4. * M_PI / 3. / vol / pme->epsilon_r * iprod(qm->xQM[j], sum_qx);
       }
  }
  else
  {
	   for (int j=0; j<n; j++) {
           pot_surf[j] = 0.;
	   }
  }

  // return the potential on QM atoms
  for (int j=0; j<n; j++)
  {
      pot[j] = pme->pot[j] * KJMOL2HARTREE + pot_corr[j] * BOHR2NM + pot_surf[j] * BOHR2NM;
   // printf("pot_qm_in_scc[%d] = %12.8f\n", j+1, pot[j]);
   // printf("Ewald atom %d charge %6.3f potential %8.5f (correction %8.5f surfterm %8.5f)\n",
   //         j+1, pme->q[j],      pot[j],                 pot_corr[j] * BOHR2NM, pot_surf[j] * BOHR2NM);
  }

  return;
} // calculate_complete_QM_QM_ewald
*/

void QMMM_rec::calculate_complete_QM_QM(const t_commrec*  cr,
                                        t_nrnb*           nrnb,
                                        gmx_wallcycle_t   wcycle,
                                        struct gmx_pme_t* pmedata,
                                        real*             pot)
{
  /* as the last arguments are expected: fr->ewaldcoeff_q and fr->pmedata */
  QMMM_QMrec& qm_          = qm[0];
  QMMM_MMrec& mm_          = mm[0];
//QMMM_PME& pme_full       = pme[0];
  QMMM_PME& pme_qmonly     = pme[1];
//t_QMMM_PME *pme          = qr->pme_qmonly;
  int         n            = qm_.nrQMatoms;
//real        rcoul        = qm_.rcoulomb;
  real        ewaldcoeff_q = qm_.ewaldcoeff_q;

///* (re)allocate */
//srenew(pme->x, qm_.nrQMatoms);
//srenew(pme->q, qm_.nrQMatoms);

  /* copy the data into PME structures */
  for (int j=0; j<n; j++)
  {
      /* QM atoms with charges */
      pme_qmonly.x[j][XX] = qm_.xQM[j][XX];
      pme_qmonly.x[j][YY] = qm_.xQM[j][YY];
      pme_qmonly.x[j][ZZ] = qm_.xQM[j][ZZ];
      /* Attenuate the periodic images of the QM zone
       * with the same scaling factor
       * that is applied for the MM atoms.
       */
      pme_qmonly.q[j]     = qm_.QMcharges[j] * mm_.scalefactor;
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
             nullptr, nullptr, nullptr, nullptr, qm_.box, cr, 0, 0, // pme_qmonly.nrnb->get(),
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
	   for (int j=0; j<qm_.nrQMatoms; j++) {
           svmul(pme_qmonly.q[j], pme_qmonly.x[j], qx);
           rvec_inc(sum_qx, qx);
	   }
	   // contribution to the potential
	   real vol = qm_.box[XX][XX] * qm_.box[YY][YY] * qm_.box[ZZ][ZZ];
	   for (int j=0; j<n; j++) {
	       pot_surf[j] = 4. * M_PI / 3. / vol / pme_qmonly.epsilon_r * iprod(qm_.xQM[j], sum_qx);
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
      qm_.pot_qmqm_set(j, static_cast<double>(pot[j]) * HARTREE_TO_EV); // in volt units
  }
} // calculate_complete_QM_QM

/**********************************
 ***  GRADIENTS       *************
 **********************************/

void QMMM_rec::gradient_QM_MM(const t_commrec*  cr,
                              t_nrnb*           nrnb,
                              gmx_wallcycle_t   wcycle,
                              struct gmx_pme_t* pmedata,
                              int               variant,
                              rvec*             partgrad,
                              rvec*             MMgrad,
                              rvec*             MMgrad_full)
{
  QMMM_QMrec& qm_ = qm[0];
  QMMM_MMrec& mm_ = mm[0];
//t_QMMM_PME *pme_full   = qr->pme_full;
//t_QMMM_PME *pme_qmonly = qr->pme_qmonly;
  real        rcoul = qm_.rcoulomb;
  real        ewaldcoeff_q = qm_.ewaldcoeff_q;
  int         n = qm_.nrQMatoms;
  int         ne = mm_.nrMMatoms;
  int         ne_full = mm_.nrMMatoms_full;
  rvec bond;

  /* all of the contributions to the gradients are calculated in, or immediately converted to,
   * ATOMIC UNITS!
   */

  for (int j=0; j<ne; j++)
  {
      clear_rvec(MMgrad[j]);
  }
  if (variant == eqmmmPME)
  {
    for (int j=0; j<ne_full; j++)
    {
      clear_rvec(MMgrad_full[j]);
    }
  }

  switch (variant)
  {

    case eqmmmSWITCH:
    {
      real r_1 = rcoul;
      real r_d = QMMM_SWITCH;
      real r_c = r_1 + r_d;
      real big_a =   (5 * r_c - 2 * r_1) / (CUB(r_c) * SQR(r_d));
      real big_b = - (4 * r_c - 2 * r_1) / (CUB(r_c) * CUB(r_d));
      for (int j=0; j<n; j++) { // do it for every QM atom
        // add SR potential only from MM atoms in the neighbor list!
        for (int k=0; k<ne; k++) {
          pbc_dx_qmmm(qm_.box, qm_.xQM[j], mm_.xMM[k], bond);
          real r = norm(bond);
          rvec dgr;
          if (r < 0.001)
          { // this may occur on the first step of simulation for link atom(s)
              printf("QM/MM PME QM--MM short range exploding for QM=%d, MM=%d. MM charge is %f\n", j+1, k+1, mm_.MMcharges[k]);
              continue;
          }
          if (r < r_1)
          {
              real fscal = - qm_.QMcharges[j] * mm_.MMcharges[k] / CUB(r) * SQR(BOHR2NM);
              svmul(fscal, bond, dgr);
              //printf("SR: QM %1d -- MM %1d:%12.7f%12.7f%12.7f\n", j+1, k+1, dgr[XX], dgr[YY], dgr[ZZ]);
              rvec_inc(partgrad[j], dgr);
              rvec_dec(MMgrad[k], dgr);
              continue;
          }
          if (r < r_c)
          {
              real fscal = - qm_.QMcharges[j] * mm_.MMcharges[k] / r * (1. / SQR(r)
                           - big_a * SQR(r - r_1) - big_b * CUB(r - r_1)) * SQR(BOHR2NM);
              svmul(fscal, bond, dgr);
              rvec_inc(partgrad[j], dgr);
              rvec_dec(MMgrad[k], dgr);
          	continue;
          }
          // else ... beyond cutoff+switch, nothing to do
        } // for k
      } // for j
      break;
    }

    case eqmmmRFIELD:
    {
      real r_c = rcoul;
      for (int j=0; j<n; j++) { // do it for every QM atom
        // add SR potential only from MM atoms in the neighbor list!
        for (int k=0; k<ne; k++) {
          pbc_dx_qmmm(qm_.box, qm_.xQM[j], mm_.xMM[k], bond);
          real r = norm(bond);
          rvec dgr;
          if (r < 0.001)
          { // this may occur on the first step of simulation for link atom(s)
              printf("QM/MM PME QM--MM short range exploding for QM=%d, MM=%d. MM charge is %f\n", j+1, k+1, mm_.MMcharges[k]);
              continue;
          }
          if (r < r_c)
          {
              real fscal = - qm_.QMcharges[j] * mm_.MMcharges[k] / r * (1. / SQR(r) - r / CUB(r_c)) * SQR(BOHR2NM);
              svmul(fscal, bond, dgr);
              rvec_inc(partgrad[j], dgr);
              rvec_dec(MMgrad[k], dgr);
              continue;
          }
          // else ... beyond cutoff, nothing to do
        } // for k
      } // for j
      break;
    }

    case eqmmmSHIFT:
    {
      real r_c = rcoul;
      real big_c = 3. / SQR(r_c);
      for (int j=0; j<n; j++) { // do it for every QM atom
        // add SR potential only from MM atoms in the neighbor list!
        for (int k=0; k<ne; k++) {
          pbc_dx_qmmm(qm_.box, qm_.xQM[j], mm_.xMM[k], bond);
          real r = norm(bond);
          rvec dgr;
          if (r < 0.001)
          { // this may occur on the first step of simulation for link atom(s)
              printf("QM/MM PME QM--MM short range exploding for QM=%d, MM=%d. MM charge is %f\n", j+1, k+1, mm_.MMcharges[k]);
              continue;
          }
          if (r < r_c)
          {
              real fscal = - qm_.QMcharges[j] * mm_.MMcharges[k] / r * (1. / SQR(r) + 2. * r / CUB(r_c) - big_c) * SQR(BOHR2NM);
              svmul(fscal, bond, dgr);
              rvec_inc(partgrad[j], dgr);
              rvec_dec(MMgrad[k], dgr);
              continue;
          }
        } // for k
      } // for j
      break;
    }

    case eqmmmPME: // QM/MM with PME 
    {
      QMMM_PME& pme_full   = pme[0];
      QMMM_PME& pme_qmonly = pme[1];
    //std::vector<rvec> grad_add(n);
      rvec* grad_add = new rvec [n];

      /* (1) gradient on QM atoms due to MM atoms. */

      // copy the data into PME structures
      for (int j=0; j<n; j++)
      {
          pme_full.x[j][XX] = qm_.xQM[j][XX];
          pme_full.x[j][YY] = qm_.xQM[j][YY];
          pme_full.x[j][ZZ] = qm_.xQM[j][ZZ];
          /* unscaled QM charges.
           * ATTENTION -- the interaction of periodic QM images will be included FULLY,
           * and thus it has to be reduced below, to account for the possibly requested scalefactor.
           */
          pme_full.q[j]     = qm_.QMcharges[j];
      }
      for (int j=0; j<ne_full; j++)
      { 
          pme_full.x[n + j][XX] = mm_.xMM_full[j][XX];
          pme_full.x[n + j][YY] = mm_.xMM_full[j][YY];
          pme_full.x[n + j][ZZ] = mm_.xMM_full[j][ZZ];
          /* the MM charges are already scaled */
          pme_full.q[n + j]     = mm_.MMcharges_full[j];
      }
      // PME -- long-range component
    //static struct timespec time1, time2;
    //clock_gettime(CLOCK_MONOTONIC, &time1);
      // init_nrnb(pme_full.nrnb); // TODO change to something?
      gmx::StepWorkload stepWork;
      stepWork.computeForces = true;
      std::vector<real> emptyVec;
      gmx_pme_do(pmedata, pme_full.x, pme_full.f, pme_full.q.data(), pme_full.q.data(),
                 nullptr, nullptr, nullptr, nullptr, qm_.box, cr, 0, 0, nrnb, // pme_full.nrnb->get(),
                 wcycle, pme_full.vir, pme_full.vir, nullptr, nullptr, 0., 0., nullptr, nullptr,
                 stepWork, TRUE, FALSE, n, nullptr); // emptyVec);
    //clock_gettime(CLOCK_MONOTONIC, &time2);
    //print_time_difference("PMETIME 3 ", time1, time2);
      for (int j=0; j<n; j++)
      {
        for (int m=0; m<DIM; m++)
        {
          grad_add[j][m] = - pme_full.f[j][m] / HARTREE_BOHR2MD; // partgrad is gradient, i.e. the negative of force
        }
      }
   // printf("================================\n");
   // for (int i=0; i<n; i++)
   // {
   //     printf("GRAD QMMM1 %d: %8.5f %8.5f %8.5f\n", i+1, grad_add[i][XX], grad_add[i][YY], grad_add[i][ZZ]);
   // }
      for (int j=0; j<n; j++) {
        rvec_inc(partgrad[j], grad_add[j]);
      }

      // PME corrections -- exclude QM--QM interaction within the basis PBC cell
      for (int j=0; j<n; j++) {
        clear_rvec(grad_add[j]);
      }
      for (int j=0; j<n; j++) {
        // Exclude the QM[j]--QM[k] interactions:
        for (int k=0; k<j; k++) {
          rvec_sub(qm_.xQM[j], qm_.xQM[k], bond);
          real r = norm(bond);
          rvec dgr;
          // negative of gradient -- we want to subtract it from partgrad
          real fscal = qm_.QMcharges[j] * qm_.QMcharges[k] / SQR(r) *
                        (gmx_erf(ewaldcoeff_q * r) / r
                       - M_2_SQRTPI * ewaldcoeff_q * exp(-SQR(ewaldcoeff_q * r))) * SQR(BOHR2NM);
          svmul(fscal, bond, dgr); // vec(dgr) = fscal * vec(bond)
          rvec_inc(grad_add[j], dgr);
          rvec_dec(grad_add[k], dgr);
        }
      }
      for (int j=0; j<n; j++) {
        rvec_inc(partgrad[j], grad_add[j]);
      }
   // printf("================================\n");
   // for (int i=0; i<n; i++)
   // {
   //     printf("GRAD CORR1 %d: %8.5f %8.5f %8.5f\n", i+1, grad_add[i][XX], grad_add[i][YY], grad_add[i][ZZ]);
   // }

      /* (2) gradient on QM atoms due to the periodic images of QM atoms.
       *     This is contained in the above contribution already, however with unscaled QM charges.
       *     The current contribution (2) is meant to correct for this,
       *       therefore it will be scaled with (scalefactor-1) at the end of the calculation!
       */

      // copy the data into PME structures
      for (int j=0; j<n; j++)
      {
          pme_qmonly.x[j][XX] = qm_.xQM[j][XX];
          pme_qmonly.x[j][YY] = qm_.xQM[j][YY];
          pme_qmonly.x[j][ZZ] = qm_.xQM[j][ZZ];
          /* unscaled QM charges; the resulting gradients will be scaled down at the end of the calculation! */
          pme_qmonly.q[j]     = qm_.QMcharges[j];
      }
      // PME -- long-range component
    //clock_gettime(CLOCK_MONOTONIC, &time1);
      // init_nrnb(pme_qmonly.nrnb); // TODO change to something?
      int oldNumAtoms = pmedata->atc[0].numAtoms(); // need to resize PME arrays to the number of QM atoms
      pmedata->atc[0].setNumAtoms(n);
      gmx_pme_do(pmedata, pme_qmonly.x, pme_qmonly.f, pme_qmonly.q.data(), pme_qmonly.q.data(),
                 nullptr, nullptr, nullptr, nullptr, qm_.box, cr, 0, 0, // pme_full.nrnb->get(),
                 nrnb, wcycle, pme_qmonly.vir, pme_qmonly.vir, nullptr, nullptr, 0., 0., nullptr, nullptr,
                 stepWork, TRUE, FALSE, n, nullptr); // emptyVec);
      pmedata->atc[0].setNumAtoms(oldNumAtoms); // resize back
    //clock_gettime(CLOCK_MONOTONIC, &time2);
    //print_time_difference("PMETIME 4 ", time1, time2);
      for (int j=0; j<n; j++)
      {
        for (int m=0; m<DIM; m++)
        {
          grad_add[j][m] = - pme_qmonly.f[j][m] / HARTREE_BOHR2MD * (mm_.scalefactor - 1.);
        }
      }
   // printf("================================\n");
   // for (int i=0; i<n; i++)
   // {
   //     printf("GRAD QMMM2 %d: %8.5f %8.5f %8.5f\n", i+1, grad_add[i][XX], grad_add[i][YY], grad_add[i][ZZ]);
   // }

      // PME corrections -- exclude QM--QM interaction within the basis PBC cell
      for (int j=0; j<n; j++) {
        clear_rvec(grad_add[j]);
      }
      for (int j=0; j<n; j++) {
        // Exclude the QM[j]--QM[k] interactions:
        for (int k=0; k<j; k++) {
          rvec_sub(qm_.xQM[j], qm_.xQM[k], bond);
          real r = norm(bond);
          rvec dgr;
          // negative of gradient -- we want to subtract it from partgrad
          real fscal = qm_.QMcharges[j] * qm_.QMcharges[k] / SQR(r) *
                        (gmx_erf(ewaldcoeff_q * r) / r
                       - M_2_SQRTPI * ewaldcoeff_q * exp(-SQR(ewaldcoeff_q * r))) * SQR(BOHR2NM);
          svmul(fscal, bond, dgr); // vec(dgr) = fscal * vec(bond)
          rvec_inc(grad_add[j], dgr);
          rvec_dec(grad_add[k], dgr);
        }
      }
      for (int j=0; j<n; j++) {
        for (int m=0; m<DIM; m++) {
          grad_add[j][m] *= (mm_.scalefactor - 1.);
        }
        rvec_inc(partgrad[j], grad_add[j]);
      }
   // printf("================================\n");
   // for (int i=0; i<n; i++)
   // {
   //     printf("GRAD CORR2 %d: %8.5f %8.5f %8.5f\n", i+1, grad_add[i][XX], grad_add[i][YY], grad_add[i][ZZ]);
   // }
      // TODO: RE-TEST THIS CONTRIBUTION WITH SCALEFACTOR != 1 !!!

      /* (3) gradient on MM atoms due to QM atoms.
       *     ASK GERRIT IF THIS IS REALLY NOT INCLUDED IN GROMACS MM CALCULATIONS!
       *     But I still think it is not included.
       */

      // copy the data into PME structures
      for (int j=0; j<n; j++)
      {
          pme_full.x[j][XX] = qm_.xQM[j][XX];
          pme_full.x[j][YY] = qm_.xQM[j][YY];
          pme_full.x[j][ZZ] = qm_.xQM[j][ZZ];
          pme_full.q[j]     = qm_.QMcharges[j];
      }
      for (int k=0; k<ne_full; k++)
      {
          pme_full.x[n + k][XX] = mm_.xMM_full[k][XX];
          pme_full.x[n + k][YY] = mm_.xMM_full[k][YY];
          pme_full.x[n + k][ZZ] = mm_.xMM_full[k][ZZ];
          pme_full.q[n + k]     = 0.;
      }
      // PME -- long-range component
    //clock_gettime(CLOCK_MONOTONIC, &time1);
      // init_nrnb(pme_full.nrnb); // TODO change to something?
      gmx_pme_do(pmedata, pme_full.x, pme_full.f, pme_full.q.data(), pme_full.q.data(),
                 nullptr, nullptr, nullptr, nullptr, qm_.box, cr, 0, 0, // pme_full.nrnb->get(),
                 nrnb, wcycle, pme_full.vir, pme_full.vir, nullptr, nullptr, 0., 0., nullptr, nullptr,
                 stepWork, TRUE, TRUE, n, nullptr); // emptyVec);
    //clock_gettime(CLOCK_MONOTONIC, &time2);
    //print_time_difference("PMETIME 5 ", time1, time2);
      for (int j=0; j<ne_full; j++)
      {
          MMgrad_full[j][XX] = - mm_.MMcharges_full[j] * pme_full.f[n + j][XX] / HARTREE_BOHR2MD;
          MMgrad_full[j][YY] = - mm_.MMcharges_full[j] * pme_full.f[n + j][YY] / HARTREE_BOHR2MD;
          MMgrad_full[j][ZZ] = - mm_.MMcharges_full[j] * pme_full.f[n + j][ZZ] / HARTREE_BOHR2MD;
      } // svmul(- mm_.MMcharges_full[j] / HARTREE_BOHR2MD, pme->f[n + j], mm_.grad_full[j]);
   // printf("================================\n");
   // for (int i=0; i<ne_full; i++)
   // {
   //     if (norm(mm_.grad_full[i]) > 0.001)
   //     printf("GRAD MM    %d: %8.5f %8.5f %8.5f\n", i+1, mm_.grad_full[i][XX], mm_.grad_full[i][YY], mm_.grad_full[i][ZZ]);
   // }

   // /* (4) Surface correction term for both QM and MM. */
   // if (pme->surf_corr_pme)
   // {
   //     // copy the data into PME structures
   //     for (int j=0; j<n; j++)
   //     {
   //         pme->x[j][XX] = qm_.xQM[j][XX];
   //         pme->x[j][YY] = qm_.xQM[j][YY];
   //         pme->x[j][ZZ] = qm_.xQM[j][ZZ];
   //         pme->q[j]     = qm_.QMcharges[j];
   //     }
   //     for (int k=0; k<ne_full; k++)
   //     {
   //         pme->x[n + k][XX] = mm_.xMM_full[k][XX];
   //         pme->x[n + k][YY] = mm_.xMM_full[k][YY];
   //         pme->x[n + k][ZZ] = mm_.xMM_full[k][ZZ];
   //         pme->q[n + k]     = mm_.MMcharges[k];
   //     }
   //     rvec qx, sum_qx, dgr;
   //     // sum_j q_j vec(x_j)
   //     clear_rvec(sum_qx);
   //     for (int j=0; j<n+ne_full; j++)
   //     {
   //       svmul(pme->q[j], pme->x[j], qx);
   //       rvec_inc(sum_qx, qx);
   //     }
   //     // is it OK that the QM charges have been scaled down possibly??? (mm_.scalefactor)
   //   //svmul(NM2BOHR, sum_qx, sum_qx); // TODO -- CONVERSION?
   //     // contribution to the potential
   //     real vol = qm_.box[XX][XX] * qm_.box[YY][YY] * qm_.box[ZZ][ZZ]; //  * CUB(NM2BOHR); // TODO -- CONVERSION?
   //     // partgrad is gradient, i.e. negative of force
   //     for (int j=0; j<qm_.nrQMatoms; j++)
   //     {
   //       svmul(4. * M_PI / 3. / vol / pme->epsilon_r * qm_.QMcharges[j], sum_qx, dgr);
   //       rvec_inc(partgrad[j], dgr);
   //     }
   //     // do this correction for MM atoms as well
   //     for (int j=0; j<mm_.nrMMatoms_full; j++) 
   //     {
   //       svmul(4. * M_PI / 3. / vol / pme->epsilon_r * mm_.MMcharges_full[j], sum_qx, dgr);
   //       rvec_inc(mm_.grad_full[j], dgr);
   //     }
   // }

      /* (5) QM--MM short-range gradient on both QM and MM atoms
       *       - only consider MM atoms on the neighbor list!
       */
      for (int j=0; j<n; j++) {
        clear_rvec(grad_add[j]);
        for (int k=0; k<ne; k++) { // do it for every QM atom
          pbc_dx_qmmm(qm_.box, qm_.xQM[j], mm_.xMM[k], bond);
          real r = norm(bond);
          rvec dgr;
          if (r < 0.001)
          { // this may occur on the first step of simulation for link atom(s)
              printf("QM/MM PME QM--MM short range exploding for QM=%d, MM=%d. MM charge is %f\n", j+1, k+1, mm_.MMcharges[k]);
              continue;
          }
          if (r < rcoul)
          {
              real fscal = qm_.QMcharges[j] * mm_.MMcharges[k] / SQR(r) *
                           (- gmx_erfc(ewaldcoeff_q * r) / r
                            - M_2_SQRTPI * ewaldcoeff_q * exp(-SQR(ewaldcoeff_q * r))) * SQR(BOHR2NM);
              svmul(fscal, bond, dgr);
              rvec_inc(grad_add[j], dgr);
              rvec_dec(MMgrad[k], dgr);
          }
        }
      }

   // printf("================================\n");
   // for (int i=0; i<n; i++)
   // {
   //     printf("GRAD QMMM5 %d: %8.5f %8.5f %8.5f\n", i+1, grad_add[i][XX], grad_add[i][YY], grad_add[i][ZZ]);
   // }
      for (int j=0; j<n; j++) {
        rvec_inc(partgrad[j], grad_add[j]);
      }

   // printf("================================\n");
   // for (int i=0; i<n; i++)
   // {
   //     if (norm(mm_.grad[i]) > 0.001)
   //     printf("GRAD MM  5 %d: %8.5f %8.5f %8.5f\n", i+1, mm_.grad[i][XX], mm_.grad[i][YY], mm_.grad[i][ZZ]);
   // }

      delete[] grad_add;
      // end of PME
      break;
    }

    case eqmmmNN:
    {
      for (int j=0; j<n; j++) { // do it for every QM atom
        // add SR potential only from MM atoms in the neighbor list!
        for (int k=0; k<ne; k++) {
          pbc_dx_qmmm(qm_.box, qm_.xQM[j], mm_.xMM[k], bond);
          real r = norm(bond);
          rvec dgr;
          if (r < 0.001)
          { // this may occur on the first step of simulation for link atom(s)
              printf("QM/MM PME QM--MM short range exploding for QM=%d, MM=%d. MM charge is %f\n", j+1, k+1, mm_.MMcharges[k]);
              continue;
          }
          else
          {
              real fscal = - qm_.QMcharges[j] * mm_.MMcharges[k] / CUB(r) * SQR(BOHR2NM);
              svmul(fscal, bond, dgr);
              rvec_inc(partgrad[j], dgr);
              rvec_dec(MMgrad[k], dgr);
              continue;
          }
        } // for k
      } // for j
      break;
    }

    default: // it should never get this far
      ;
  } // switch (cutoff_qmmm)
}

