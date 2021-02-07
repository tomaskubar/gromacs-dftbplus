#include "transfer.h"

// it is all dummy, so let us use an empty function for this
void slkmatrices(int i, int j, double (*xat)[3], double ham[LDIM][LDIM], double over[LDIM][LDIM],
                 int lmax[DFTB_MAXTYPES], int dim[DFTB_MAXTYPES][DFTB_MAXTYPES], double dr[DFTB_MAXTYPES][DFTB_MAXTYPES],
                 int *izp, tendoubles *skstab[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *skhtab[DFTB_MAXTYPES][DFTB_MAXTYPES], double skself[DFTB_MAXTYPES][3]);
void slkmatrices(int i, int j, double (*xat)[3], double ham[LDIM][LDIM], double over[LDIM][LDIM],
                 int lmax[DFTB_MAXTYPES], int dim[DFTB_MAXTYPES][DFTB_MAXTYPES], double dr[DFTB_MAXTYPES][DFTB_MAXTYPES],
                 int *izp, tendoubles *skstab[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *skhtab[DFTB_MAXTYPES][DFTB_MAXTYPES], double skself[DFTB_MAXTYPES][3])
{
  (void) i;
  (void) j;
  (void) xat;
  (void) ham;
  (void) over;
  (void) lmax;
  (void) dim;
  (void) dr;
  (void) izp;
  (void) skstab;
  (void) skhtab;
  (void) skself;
}


void gammaall1(double r, double ui, double uj, double *dcdr);
void gammaall1(double r, double ui, double uj, double *dcdr)
{
//modified gammaall1: no 3rd order and no gamma_h
        const double zero=1.e-4;
        double r2,a,b,a2,a4,b2,b4,z,z2,zr,g,dgdr,dsdr,fab,fba,dfabdr,dfbadr;

        r2 = r*r;
        a  = 3.2 * ui;
        b  = 3.2 * uj;

        if ((a+b < zero) || (r < zero)) {
               *dcdr = 0.;
               //*dcdr3= 0.;
               r    = 99999999.9; /* WHY ??? */
        } else { /* here: 1/r-s */
                if (fabs(a-b) < 1.e-5) {
                        z    = 0.5 * (a+b);
                        z2   = z*z;
                        zr   = z*r;
                        g    = (48. + 33.*zr + 9.*zr*zr + zr*zr*zr) / (48.*r);
                        dgdr = -1./r2 + 3.*z2/16. + z2*zr/24.;
                        dsdr = exp(-zr) * (dgdr - z*g);
                } else {
                        a2  = a*a;
                        a4  = a2*a2;
                        b2  = b*b;
                        b4  = b2*b2;
                        fab = a*b4 / (2.*SQR(a2-b2)) - (b4*b2-3.*a2*b4) / (CUB(a2-b2)*r);
                        fba = b*a4 / (2.*SQR(b2-a2)) - (a4*a2-3.*b2*a4) / (CUB(b2-a2)*r);
                        dfabdr    = (b4*b2 - 3.*a2*b4) / (CUB(a2-b2)*r2);
                        dfbadr    = (a4*a2 - 3.*b2*a4) / (CUB(b2-a2)*r2);
                        dsdr      = exp(-a*r) * (dfabdr-a*fab) + exp(-b*r) * (dfbadr-b*fba);
                }
                *dcdr   = -1./r2 - dsdr;
        } /* end 1/r-s */
        return;
}


void offdiag_gradient_homo(dftb_t *dftb, dvec *x, dvec *grad, charge_transfer_t *ct){
/* calculates how a missing(additional) electron, which is distributed over several HOMO(LUMO) of neighboring fragments gives rise to attracting or repulsive forces,
   depending on bonding or antibonding linear combination of FOs */

// PARAMETERS:
// dftb  = (in) main data structure for DFTB calculations
// x     = (in) coordinates of all atoms of the complex
// grad  = (out) gradient vector for all atoms of the complex
// ct    = (in) main data structure with information about the charge transfer calculation
////////////////////////////////////////////////////////////////////////

  int i, j, k, izpj, izpk, site_i, site_j; 
  int m, n, indj, indk, lj, mj, nu, lk, mk, mu;
  double ocmcc, xhelp, dgrh;
  const double deltax = 0.01;
  dftb_phase1_t *dftb1 = dftb->phase1;
  dftb_phase2_t dftb2 = dftb->phase2;

  for (i=0; i<dftb2.nn; i++)
    clear_dvec(grad[i]);

  for (m=0; m<dftb2.norb; m++)
    for (n=0; n<dftb2.norb; n++)
      dftb2.b[m][n] = 0.e0;

  //transform eigenvectors in orthogonal basis so they can be used for an EXACT AO representation of ct->wf
  //tested and this has nearly no influence on the results since orthogonalization yields linar combination like 99%:1%
/*  for (l=0; l<ct->dim; l++)
    for (iao=0; iao < dftb->phase2.norb; iao++)
      dftb->orthogo.evec_ao[iao][l] = 0.0;
  for (l=0; l<ct->dim; l++){
  k=0;
  for (i=0; i < ct->sites; i++)
  for (ii = 0; ii < ct->site[i].homos; ii++) {
    ifo = ct->site[i].homo[ii] + dftb->phase2.inf[i] - 1;
      for (iao=0; iao < dftb->phase2.norb; iao++){
        dftb->orthogo.evec_ao[iao][l] += dftb->orthogo.sij[k + l * ct->dim] * dftb->phase2.Taf[iao][ifo];
      }
      k++;
  }
  }
*/

  if (ct->do_ncv)
  {
      m=0;
      for (site_i=0; site_i < ct->sites ; site_i++)
      for (i=0; i<ct->site[site_i].homos; i++) {
        n=0;
        for (site_j=0; site_j < site_i ; site_j++)//reduced loop to lower triangular block. we only need this
        for (j=0; j<ct->site[site_j].homos; j++) {
          for (mu=0; mu<dftb1[site_i].norb; mu++)
          for (nu=0; nu<dftb1[site_j].norb; nu++){
          //ocmcc = (a_m^* a_n )*c_mu^n c_nu^n
//            printf("%d  %d \n", ct->ncv_isurf, ct->ncv_jsurf);
              ocmcc = dftb1[site_i].a[mu][ct->site[site_i].homo[i]-1] * dftb1[site_j].a[nu][ct->site[site_j].homo[j]-1];
              if (!ct->pure_forces){
                ocmcc *= ct->tfs_vector[ct->ncv_isurf][m]*ct->tfs_vector[ct->ncv_jsurf][n];
              }
               //dftb2.b[dftb2.inf[site_i]+mu][dftb2.inf[site_j]+nu] = (ct->wf[m]*ct->wf[n]+ ct->wf[m+ct->dim]*ct->wf[n+ct->dim]) * dftb1[site_i].a[mu][ct->site[site_i].homo[i]-1] * dftb1[site_j].a[nu][ct->site[site_j].homo[j]-1];
               //dftb2.b[dftb2.inf[site_i]+mu][dftb2.inf[site_j]+nu] = (ct->wf[m]*ct->wf[n]+ ct->wf[m+ct->dim]*ct->wf[n+ct->dim]) * dftb->orthogo.evec_ao[dftb2.inf[site_i]+mu][m] * dftb->orthogo.evec_ao[dftb2.inf[site_j]+nu][n];
               //dftb2.b[dftb2.inf[site_i]+mu][dftb2.inf[site_j]+nu] = ct->is_hole_transfer ? -1.0* dftb2.b[dftb2.inf[site_i]+mu][dftb2.inf[site_j]+nu] : dftb2.b[dftb2.inf[site_i]+mu][dftb2.inf[site_j]+nu];
               dftb2.b[dftb2.inf[site_i]+mu][dftb2.inf[site_j]+nu] += ct->is_hole_transfer ? -1.0*ocmcc : ocmcc ;
             }
             n++;
           }
           m++;
         }
  }
  else
  {
      m=0;
      for (site_i=0; site_i < ct->sites ; site_i++)
      for (i=0; i<ct->site[site_i].homos; i++) {
        n=0;
        for (site_j=0; site_j < site_i ; site_j++)//reduced loop to lower triangular block. we only need this
        for (j=0; j<ct->site[site_j].homos; j++) {
          for (mu=0; mu<dftb1[site_i].norb; mu++)
          for (nu=0; nu<dftb1[site_j].norb; nu++){
          //ocmcc = (a_m^* a_n )*c_mu^n c_nu^n
//            ocmcc = dftb1[site_i].a[mu][ct->site[site_i].homo[i]-1] * dftb1[site_j].a[nu][ct->site[site_j].homo[j]-1];
            ocmcc =  dftb1[site_i].a[mu][ct->site[site_i].homo[i]-1] * dftb1[site_j].a[nu][ct->site[site_j].homo[j]-1];
            if (!ct->pure_forces){
              ocmcc *= (ct->wf[m]*ct->wf[n] + ct->wf[m+ct->dim]*ct->wf[n+ct->dim]);
            }
            //dftb2.b[dftb2.inf[site_i]+mu][dftb2.inf[site_j]+nu] = (ct->wf[m]*ct->wf[n]+ ct->wf[m+ct->dim]*ct->wf[n+ct->dim]) * dftb1[site_i].a[mu][ct->site[site_i].homo[i]-1] * dftb1[site_j].a[nu][ct->site[site_j].homo[j]-1];
            //dftb2.b[dftb2.inf[site_i]+mu][dftb2.inf[site_j]+nu] = (ct->wf[m]*ct->wf[n]+ ct->wf[m+ct->dim]*ct->wf[n+ct->dim]) * dftb->orthogo.evec_ao[dftb2.inf[site_i]+mu][m] * dftb->orthogo.evec_ao[dftb2.inf[site_j]+nu][n];
            //dftb2.b[dftb2.inf[site_i]+mu][dftb2.inf[site_j]+nu] = ct->is_hole_transfer ? -1.0* dftb2.b[dftb2.inf[site_i]+mu][dftb2.inf[site_j]+nu] : dftb2.b[dftb2.inf[site_i]+mu][dftb2.inf[site_j]+nu];
            dftb2.b[dftb2.inf[site_i]+mu][dftb2.inf[site_j]+nu] += ct->is_hole_transfer ? -1.0*ocmcc : ocmcc ;
          }
          n++;
          }
          m++;
        }
  }
 //printf("FORCE: b-matrix end at %f\n",  (double) clock()/CLOCKS_PER_SEC);

/*
  for (j=0; j<ct->site[site_i].homos; j++) {
    i=ct->site[site_i].homo[j]-1;
    for (m=0; m<dftb1.norb; m++)
      for (n=0; n<m; n++) {
      //ocmcc = dftb1.occ[i] * dftb1.a[m][i] * dftb1.a[n][i];
        ocmcc = ct->is_hole_transfer ? (-1.0 * dftb1.a[m][i] * dftb1.a[n][i] * ct->occupation[ct->indFO[site_i]+j]) : (1.0 * dftb1.a[m][i] * dftb1.a[n][i] * ct->occupation[ct->indFO[site_i]+j]) ;
        //sum over all occupied orbitals changed to sum over frontier orbitals weighted by their occupation.
        dftb2.b[n][m] += ocmcc * dftb1.ev[i];
        dftb2.b[m][n] += ocmcc;
      }
  }
*/
  for (m=0; m<dftb2.norb; m++)
    for (n=0; n<m; n++)
      if (fabs(dftb2.b[m][n]) < dftb->dacc)
        dftb2.b[m][n] = 0.e0;

  // the matrix b is correct here! //
  for (j=0; j<dftb2.nn; j++) { // for every atom that forces act upon
    for (k=0; k<dftb2.nn; k++) { // for every atom acting on the studied atom j
      // neighborsearching disabled, maybe temporarily! TODO
      if (k != j) { // && dftb->nl[j][k]) { // no force on self and check if k is neighbor of j
        indj = dftb2.ind[j];
        izpj = dftb2.izp[j];
        indk = dftb2.ind[k];
        izpk = dftb2.izp[k];
        for (i=0; i<3; i++) { // XX, YY and ZZ
          // derivative of the slko matrices
          xhelp = dftb2.x[j][i];
          x[j][i] = xhelp + deltax;
          //if(site_i == site_j){ //confined basis
          if(0 == 1){ //confined basis
            slkmatrices(k, j, dftb2.x, dftb2.au,  dftb2.bu,  dftb->lmax, dftb->dim1, dftb->dr1, dftb2.izp, dftb->skstab1, dftb->skhtab1, dftb->skself1);
            x[j][i] = xhelp - deltax;
            slkmatrices(k, j, dftb2.x, dftb2.auh, dftb2.buh, dftb->lmax, dftb->dim1, dftb->dr1, dftb2.izp, dftb->skstab1, dftb->skhtab1, dftb->skself1);
            x[j][i] = xhelp;
          }else{ // unconfined basis
            slkmatrices(k, j, dftb2.x, dftb2.au,  dftb2.bu,  dftb->lmax, dftb->dim2, dftb->dr2, dftb2.izp, dftb->skstab2, dftb->skhtab2, dftb->skself2);
            x[j][i] = xhelp - deltax;
            slkmatrices(k, j, dftb2.x, dftb2.auh, dftb2.buh, dftb->lmax, dftb->dim2, dftb->dr2, dftb2.izp, dftb->skstab2, dftb->skhtab2, dftb->skself2);
            x[j][i] = xhelp;
          }

          // use summation over angular momentum and magnetic quantum numbers
          // because shift is actually defined for the orbitals
          for (lj=1; lj<=dftb->lmax[izpj]; lj++)
            for (mj=1; mj<2*lj; mj++) {
              n  = SQR(lj-1) + mj - 1;
              nu = n + indj;
              for (lk=1; lk<=dftb->lmax[izpk]; lk++)
                for (mk=1; mk<2*lk; mk++) {
                  m  = SQR(lk-1) + mk - 1;
                  mu = m + indk;
                  // dgrh = 2 * ( d H_{k,j}^{m,n}/ d R_{j} )
                  dgrh = (dftb2.au[m][n] - dftb2.auh[m][n]) / deltax;

                  // only lower triangle contains sum_i n(i) c_{mu,i} c_{nu,i}
                  // only upper triangle contains sum_i epsilon_i n(i) c_{mu,i} c_{nu,i} //
                  if (mu > nu) {
                    dgrh *= dftb2.b[mu][nu];
                  } else {
                    dgrh *= dftb2.b[nu][mu];
                  }
// Xie, March 11th
                  dgrh *= ct->offdiag_scaling;

                  grad[j][i] += dgrh; //+ dgrs + dgr;
                }
            }
        }
      }
    }
//printf("offdiag grad on atom %d sum %f ",j, fabs(grad[j][0])+fabs(grad[j][1])+fabs(grad[j][2]));
  }
  return;
}


void usual_gradient_homo(dftb_t *dftb, dvec *x, dvec *grad, charge_transfer_t *ct, int site_i){
/* calculates the fraction of the HOMO(LUMO) orbital to the total gradient of the neutral fragment
   this can be seen as approximative deltaForce term between the neutral reference system and the charged one. */

// PARAMETERS:
// dftb   = (in) main data structure for DFTB calculations
// x      = (in) coordinates of all atoms of this fragment
// grad   = (out) gradient vector for all atoms of this fragment
// ct     = (in) main data structure with information about the charge transfer calculation
// site_i = (in) index of the fragment
////////////////////////////////////////////////////////////////////////
  int i, j, k, izpj, izpk;
  int m, n, indj, indk, lj, mj, nu, lk, mk, mu;
  double ocmcc, xhelp, dgrh, dgrs, dgr;

  const double deltax = 0.01;
  dftb_phase1_t dftb1 = dftb->phase1[site_i];

  dgr = 0.;
//dgr3 = 0.;

  for (m=0; m<dftb1.norb; m++)
    for (n=0; n<dftb1.norb; n++)
      dftb1.b[m][n] = 0.e0;

  //for (i=0; i<dftb1.norb; i++) {
  //  if (dftb1.occ[i] < dftb->dacc)
  //    break;
  /*for (j=0; j<ct->site[site_i].homos; j++) {
    i=ct->site[site_i].homo[j]-1;
    for (m=0; m<dftb1.norb; m++)
      for (n=0; n<m; n++) {
      //ocmcc = dftb1.occ[i] * dftb1.a[m][i] * dftb1.a[n][i];
        ocmcc = ct->is_hole_transfer ? (-1.0 * dftb1.a[m][i] * dftb1.a[n][i] * ct->occupation[ct->indFO[site_i]+j]) : (1.0 * dftb1.a[m][i] * dftb1.a[n][i] * ct->occupation[ct->indFO[site_i]+j]) ;
        //sum over all occupied orbitals changed to sum over frontier orbitals weighted by their occupation.
        dftb1.b[n][m] += ocmcc * dftb1.ev[i];
        dftb1.b[m][n] += ocmcc;
      }
  }
  for (m=0; m<dftb1.norb; m++)
    for (n=0; n<m; n++)
      if (fabs(dftb1.b[m][n]) < dftb->dacc)
        dftb1.b[m][n] = 0.e0;
  */
  if (ct->do_ncv)
  {
  //new b matrix including off-diagonal elements originating from two orbitals on the same site.
     for (i=0; i<ct->site[site_i].homos; i++)
     for (j=0; j<ct->site[site_i].homos; j++) {
     //for (j=i; j<ct->site[site_i].homos; j++) {
       m=ct->site[site_i].homo[i]-1;
       n=ct->site[site_i].homo[j]-1;
       for (mu=0; mu<dftb1.norb; mu++)
       for (nu=0; nu<mu; nu++) {
          //ocmcc = (a_m^* a_n)*c_mu^n c_nu^n
          ocmcc = dftb1.a[mu][m] * dftb1.a[nu][n];
          if (!ct->pure_forces){
            ocmcc *= ct->tfs_vector[ct->ncv_isurf][ct->indFO[site_i]+i]*ct->tfs_vector[ct->ncv_jsurf][ct->indFO[site_i]+j];
          }
           //sum over all occupied orbitals changed to sum over frontier orbitals weighted by their occupation.
           //dftb1.b[nu][mu] += ct->is_hole_transfer ? -1.0*ocmcc* dftb1.ev[m] : ocmcc* dftb1.ev[m] ;
           //dftb1.b[mu][nu] += ct->is_hole_transfer ? -1.0*ocmcc : ocmcc;
           dftb1.b[nu][mu] += m == n ? ocmcc* dftb1.ev[m] : 0.0 ; //norm constraint
           dftb1.b[mu][nu] +=          ocmcc;
       }
     }
  }
  else
  {
  //new b matrix including off-diagonal elements originating from two orbitals on the same site.
     for (i=0; i<ct->site[site_i].homos; i++)
     for (j=0; j<ct->site[site_i].homos; j++) {
     //for (j=i; j<ct->site[site_i].homos; j++) {
       m=ct->site[site_i].homo[i]-1;
       n=ct->site[site_i].homo[j]-1;
       for (mu=0; mu<dftb1.norb; mu++)
       for (nu=0; nu<mu; nu++) {
          //ocmcc = (a_m^* a_n)*c_mu^n c_nu^n
          ocmcc = dftb1.a[mu][m] * dftb1.a[nu][n];
          if (!ct->pure_forces){
            ocmcc *= (ct->wf[ct->indFO[site_i]+i]*ct->wf[ct->indFO[site_i]+j] + ct->wf[ct->indFO[site_i]+i+ct->dim]*ct->wf[ct->indFO[site_i]+j+ct->dim]);
          }
           //sum over all occupied orbitals changed to sum over frontier orbitals weighted by their occupation.
           //dftb1.b[nu][mu] += ct->is_hole_transfer ? -1.0*ocmcc* dftb1.ev[m] : ocmcc* dftb1.ev[m] ;
           //dftb1.b[mu][nu] += ct->is_hole_transfer ? -1.0*ocmcc : ocmcc;
           dftb1.b[nu][mu] += m == n ? ocmcc* dftb1.ev[m] : 0.0 ; //norm constraint
           dftb1.b[mu][nu] +=          ocmcc;
       }
     }
  }
  for (mu=0; mu<dftb1.norb; mu++)
    for (nu=0; nu<dftb1.norb; nu++)
      if (fabs(dftb1.b[mu][nu]) < dftb->dacc){
        dftb1.b[mu][nu] = 0.e0;
      }else if (ct->is_hole_transfer){
         dftb1.b[mu][nu] *= -1.0 ;
      }

  // the matrix b is correct here! //

  for (j=0; j<dftb1.nn; j++) { // for every atom that forces act upon
    indj = dftb1.ind[j];
    izpj = dftb1.izp[j];
    for (k=0; k<dftb1.nn; k++) if (k != j) { // for every atom acting on the studied atom j
      indk = dftb1.ind[k];
      izpk = dftb1.izp[k];

      for (i=0; i<3; i++) { // XX, YY and ZZ
        // derivative of the slko matrices
        xhelp = dftb1.x[j][i];
	x[j][i] = xhelp + deltax;
	slkmatrices(k, j, dftb1.x, dftb1.au,  dftb1.bu,  dftb->lmax, dftb->dim1, dftb->dr1, dftb1.izp, dftb->skstab1, dftb->skhtab1, dftb->skself1);
	x[j][i] = xhelp - deltax;
	slkmatrices(k, j, dftb1.x, dftb1.auh, dftb1.buh, dftb->lmax, dftb->dim1, dftb->dr1, dftb1.izp, dftb->skstab1, dftb->skhtab1, dftb->skself1);
	x[j][i] = xhelp;
	// use summation over angular momentum and magnetic quantum numbers
        // because shift is actually defined for the orbitals
        for (lj=1; lj<=dftb->lmax[izpj]; lj++)
          for (mj=1; mj<2*lj; mj++) {
            n  = SQR(lj-1) + mj - 1;
            nu = n + indj;
            for (lk=1; lk<=dftb->lmax[izpk]; lk++)
              for (mk=1; mk<2*lk; mk++) {
                m  = SQR(lk-1) + mk - 1;
                mu = m + indk;
                // dgrh = 2 * ( d H_{k,j}^{m,n}/ d R_{j} )
                dgrh = (dftb1.au[m][n] - dftb1.auh[m][n]) / deltax;
                // dgrs = - 2 * ( d S_{k,j}^{m,n}/ d R_{j} )
                dgrs = -(dftb1.bu[m][n] - dftb1.buh[m][n]) / deltax;
                // dgr = ( d S_{k,j}^{m,n}/ d R_{j} ) * (shift(k)+shift(j))
                // dgr3 =  ( d S_{k,j)^{m,n}/ d R_{j} ) * (2*shift3(k)+shift3A(k))/3.0
                dgr = -0.5 * dgrs * (dftb1.shift[k] + dftb1.shift[j]);
                // NOTE THAT WITH CDKO, shiftE2 would have to come in here as well !!!

                //if (dftb->sccmode == 3)
                  //dgr3 = -0.5 * dgrs * (2.*dftb1.shift3[k] + dftb1.shift3a[k] + 2.*dftb1.shift3[j] + dftb1.shift3a[j])/3.;

                // only lower triangle contains sum_i n(i) c_{mu,i} c_{nu,i}
                // only upper triangle contains sum_i epsilon_i n(i) c_{mu,i} c_{nu,i} //
                if (mu > nu) {
                  dgrh *= dftb1.b[mu][nu];
                  dgrs *= dftb1.b[nu][mu];
                  dgr  *= dftb1.b[mu][nu];
                  //if (dftb->sccmode == 3) dgr3 *= dftb1.b[mu][nu];
                } else {
                  dgrh *= dftb1.b[nu][mu];
                  dgrs *= dftb1.b[mu][nu];
                  dgr  *= dftb1.b[nu][mu];
                  //if (dftb->sccmode == 3) dgr3 *= dftb1.b[nu][mu];
                }
                //grad[j][i] += dgrh + dgrs + dgr;
                grad[j][i] += dgrh + dgrs; //+ dgr; use DFTB1 forces
//printf("usualgrad on atom %d sum %f H %f   norm %f    2nd %f \n",j, grad[j][i], dgrh, dgrs, dgr);
                //if (dftb->sccmode == 3) grad[j][i] += dgr3;
              }
          }
      }

    }
  }
  return;
}


void gamma_gradient_homo(dftb_t *dftb, dvec *x, dvec *grad, charge_transfer_t *ct, int site_i){
//calculates the fraction of the HOMO orbital to the total grad on the neutral molecule
//total grad is grad of charged system minus grad of neutral one
        int ix,k,j,l;
        dvec tmp, r;
        double bond,dgdrkj,dgdrjk,qdiffj,qdiffk;

        dftb_phase1_t dftb1 = dftb->phase1[site_i];
        int nn = dftb1.nn;
        int *izp = dftb1.izp;
        double *uhubb = dftb->uhubb1;
        //double *uhder = 0.0; //dummy. no hubbard derivative
        //int *izpxh = dftb1.izpxh; // no gamma_h
        double *qmat = dftb1.qmat;
        double *qzero = dftb->qzero1;

        for (k=0; k<nn; k++)
                clear_dvec(grad[k]);
        /*  get dgamma/dr and dGamma/dr   (r=|R_j-R_k|)
         *  change with respect to original DFTB3:
         *  get dr/dR_{kx} and build gammagradient contribution
         *  right away in this cycle!
         */
	//loop for charged system
        for (k=0; k<nn; k++) {
                clear_dvec(tmp);
		qdiffk =0;
                for(l=0; l<ct->site[site_i].homos; l++)
                  qdiffk -= ct->site[site_i].delta_q[l][k] * ct->occupation[ct->indFO[site_i]+l]; //qdiff is electron number
                qdiffk = qmat[k] - qzero[izp[k]] + (ct->is_hole_transfer ? -qdiffk : qdiffk );
                for (j=0; j<nn; j++) if (j != k) {
                        dvec_sub(x[k], x[j], r);
                        bond = dnorm(r);
                        gammaall1(bond, uhubb[izp[k]], uhubb[izp[j]], &dgdrkj);
                        // fprintf(stderr, "gammaall1(%d,%d): dgdr = %7.5f, dgdr3 = %7.5f\n", k+1, j+1, dgdrkj, dgdr3kj);
                        gammaall1(bond, uhubb[izp[j]], uhubb[izp[k]], &dgdrjk);
                        // fprintf(stderr, "gammaall1(%d,%d): dgdr = %7.5f, dgdr3 = %7.5f\n", j+1, k+1, dgdrjk, dgdr3jk);
                        /*dgdr[k][j]  = dgdrkj;
                          dgdr3[k][j] = dgdr3kj; */

                        qdiffj =0;
                        for(l=0; l<ct->site[site_i].homos; l++)
                          qdiffj -= ct->site[site_i].delta_q[l][j] * ct->occupation[ct->indFO[site_i]+l]; //qdiff is electron number
                        qdiffj = qmat[j] - qzero[izp[j]] + (ct->is_hole_transfer ? -qdiffj : qdiffj );
                        for (ix=0; ix<3; ix++) {
                                tmp[ix] += qdiffj * (dgdrjk + dgdrkj) / bond * r[ix];
                        }
                }

                        for (ix=0; ix<3; ix++)
                                grad[k][ix] = qdiffk * 0.5 * tmp[ix];
        }
	//loop for neutral system
        for (k=0; k<nn; k++) {
                clear_dvec(tmp);
                qdiffk = qmat[k] - qzero[izp[k]];
                //qdiffk is deltaQ of neutral system
                for (j=0; j<nn; j++) if (j != k) {
                        dvec_sub(x[k], x[j], r);
                        bond = dnorm(r);
                        gammaall1(bond, uhubb[izp[k]], uhubb[izp[j]], &dgdrkj);
                        // fprintf(stderr, "gammaall1(%d,%d): dgdr = %7.5f, dgdr3 = %7.5f\n", k+1, j+1, dgdrkj, dgdr3kj);
                        gammaall1(bond, uhubb[izp[j]], uhubb[izp[k]], &dgdrjk);
                        // fprintf(stderr, "gammaall1(%d,%d): dgdr = %7.5f, dgdr3 = %7.5f\n", j+1, k+1, dgdrjk, dgdr3jk);
                        /*dgdr[k][j]  = dgdrkj;
                          dgdr3[k][j] = dgdr3kj; */

                        qdiffj = qmat[j] - qzero[izp[j]];
                        for (ix=0; ix<3; ix++) {
                                tmp[ix] += qdiffj * (dgdrjk + dgdrkj) / bond * r[ix];
                        }
                }

                        for (ix=0; ix<3; ix++)
                                grad[k][ix] -= qdiffk * 0.5 * tmp[ix];
        }
        return;
}


void additional_gradient_homo(dftb_t *dftb, dvec *x, dvec *grad, charge_transfer_t *ct, int site_i)
{
//calculates the fraction of the HOMO orbital to the total grad on the neutral molecule
  int i, j, k,l, izpj, izpk;
  int m, n, indj, indk, lj, mj, nu, lk, mk, mu;
  double ocmcc, xhelp, dgrs, dgr;

  const double deltax = 0.01;
  dftb_phase1_t dftb1 = dftb->phase1[site_i];
  int nn;
  nn= dftb1.nn;

  dgr = 0.;

  // calculate change in atomic hamilton shift (= sum over gamma* deltacharge)
  for (i=0; i<nn; i++) {
    dftb1.dshift[i] = 0.0;
    dftb1.shift[i] = 0.0;
    for (j=0; j<nn; j++)
      for(l=0; l<ct->site[site_i].homos; l++){
//printf(" orb %d at %d dQ =  %f occ %f\n", l,j ,  ct->site[site_i].delta_q[l][j], ct->occupation[ct->indFO[site_i]+l] );
        dftb1.dshift[i] += ct->site[site_i].delta_q[l][j] * ct->occupation[ct->indFO[site_i]+l] * (i>j ? dftb1.gammamat[i][j] : dftb1.gammamat[j][i]); //minus because delta_q is actual charge compared to qdiff counting electrons
//        dftb1.shift[i] += (dftb1.qmat[j] - dftb->qzero1[dftb1.izp[j]]) * (i>j ? dftb1.gammamat[i][j] : dftb1.gammamat[j][i]);

      }
printf("dshift %d = %f dQ= %f\n", i, dftb1.dshift[i],  ct->site[site_i].delta_q[0][i] );
//printf("shift %d = %f dQ= %f\n", i, dftb1.shift[i],  (dftb1.qmat[i] - dftb->qzero1[dftb1.izp[i]]) );

  }

  for (m=0; m<dftb1.norb; m++)
    for (n=0; n<dftb1.norb; n++)
      dftb1.b[m][n] = 0.e0;

  for (i=0; i<dftb1.norb; i++) {
    if (dftb1.occ[i] < dftb->dacc)
      break;
    for (m=0; m<dftb1.norb; m++)
      for (n=0; n<m; n++) {
        ocmcc = dftb1.occ[i] * dftb1.a[m][i] * dftb1.a[n][i];
        //ocmcc = 1.0          * dftb1.a[m][i] * dftb1.a[n][i] * ct->occupation[ct->indFO[site_i]+j];
        //sum over all occupied orbitals changed to sum over frontier orbitals weighted by their occupation /////////////////////////////////////////////////////////
        dftb1.b[n][m] += ocmcc ; //* dftb1.ev[i];
        dftb1.b[m][n] += ocmcc;
      }
  }
  for (m=0; m<dftb1.norb; m++)
    for (n=0; n<m; n++)
      if (fabs(dftb1.b[m][n]) < dftb->dacc)
        dftb1.b[m][n] = 0.e0;

  // the matrix b is correct here! //

  for (j=0; j<dftb1.nn; j++) { // for every atom that forces act upon
    indj = dftb1.ind[j];
    izpj = dftb1.izp[j];
    for (k=0; k<dftb1.nn; k++) if (k != j) { // for every atom acting on the studied atom j
      indk = dftb1.ind[k];
      izpk = dftb1.izp[k];

      for (i=0; i<3; i++) { // XX, YY and ZZ
        // derivative of the slko matrices
        xhelp = dftb1.x[j][i];
	x[j][i] = xhelp + deltax;
	slkmatrices(k, j, dftb1.x, dftb1.au,  dftb1.bu,  dftb->lmax, dftb->dim1, dftb->dr1, dftb1.izp, dftb->skstab1, dftb->skhtab1, dftb->skself1);
	x[j][i] = xhelp - deltax;
	slkmatrices(k, j, dftb1.x, dftb1.auh, dftb1.buh, dftb->lmax, dftb->dim1, dftb->dr1, dftb1.izp, dftb->skstab1, dftb->skhtab1, dftb->skself1);
	x[j][i] = xhelp;
	// use summation over angular momentum and magnetic quantum numbers
        // because shift is actually defined for the orbitals
        for (lj=1; lj<=dftb->lmax[izpj]; lj++)
          for (mj=1; mj<2*lj; mj++) {
            n  = SQR(lj-1) + mj - 1;
            nu = n + indj;
            for (lk=1; lk<=dftb->lmax[izpk]; lk++)
              for (mk=1; mk<2*lk; mk++) {
                m  = SQR(lk-1) + mk - 1;
                mu = m + indk;
                // dgrh = 2 * ( d H_{k,j}^{m,n}/ d R_{j} )
                //dgrh = (dftb1.au[m][n] - dftb1.auh[m][n]) / deltax;
                // dgrs = - 2 * ( d S_{k,j}^{m,n}/ d R_{j} )
                dgrs = -(dftb1.bu[m][n] - dftb1.buh[m][n]) / deltax;
                // dgr = ( d S_{k,j}^{m,n}/ d R_{j} ) * (shift(k)+shift(j))
                //dgr = -0.5 * dgrs * (dftb1.dshift[k] + dftb1.dshift[j]);
                dgr = -0.5 * dgrs * (dftb1.dshift[k] + dftb1.dshift[j]);


                // only lower triangle contains sum_i n(i) c_{mu,i} c_{nu,i}
                // only upper triangle contains sum_i epsilon_i n(i) c_{mu,i} c_{nu,i} //
                if (mu > nu) {
                  dgrs *= dftb1.b[nu][mu];
                  dgr  *= dftb1.b[mu][nu];
                } else {
                  //dgrh *= dftb1.b[nu][mu];
                  dgrs *= dftb1.b[mu][nu];
                  dgr  *= dftb1.b[nu][mu];
                }
                grad[j][i] +=  dgr ;
//printf(" addigrad on atom %d grad %15.12f  current norm %15.12f   current 2nd %20.18f \n",j, grad[j][i], dgrs, dgr);
              }
          }
      }
    }
  }
  return;
}


void additional_gradient_homo_new(dftb_t *dftb, dvec *x, dvec *grad, charge_transfer_t *ct, int site_i)
{
//calculates the fraction of the HOMO orbital to the total grad on the neutral molecule
  int i, j, k,l, izpj, izpk;
  int m, n, indj, indk, lj, mj, nu, lk, mk, mu;
  double ocmcc, xhelp, dgrs, dgr;

  const double deltax = 0.01;
  dftb_phase1_t dftb1 = dftb->phase1[site_i];
  int nn;
  nn= dftb1.nn;

  dgr = 0.;

  // calculate change in atomic hamilton shift (= sum over gamma* deltacharge)
  for (i=0; i<nn; i++) {
    dftb1.dshift[i] = 0.0;
    for (j=0; j<dftb1.nn; j++)
      for(l=0; l<ct->site[site_i].homos; l++)
        dftb1.dshift[i] -= ct->site[site_i].delta_q[l][j] * ct->occupation[ct->indFO[site_i]+l] * (i>j ? dftb1.gammamat[i][j] : dftb1.gammamat[j][i]); //minus because delta_q is actual charge compared to qdiff counting electrons
printf("shift %d = %f\n", i, dftb1.dshift[i]);
  }

  for (m=0; m<dftb1.norb; m++)
    for (n=0; n<dftb1.norb; n++)
      dftb1.b[m][n] = 0.e0;

  for (i=0; i<dftb1.norb; i++) {
    if (dftb1.occ[i] < dftb->dacc)
      break;
    for (m=0; m<dftb1.norb; m++)
      for (n=0; n<m; n++) {
        ocmcc = dftb1.occ[i] * dftb1.a[m][i] * dftb1.a[n][i];
        //ocmcc = 1.0          * dftb1.a[m][i] * dftb1.a[n][i] * ct->occupation[ct->indFO[site_i]+j];
        //sum over all occupied orbitals changed to sum over frontier orbitals weighted by their occupation /////////////////////////////////////////////////////////
        dftb1.b[n][m] += ocmcc; //* dftb1.ev[i];
        dftb1.b[m][n] += ocmcc;
      }
  }
  for (m=0; m<dftb1.norb; m++)
    for (n=0; n<m; n++)
      if (fabs(dftb1.b[m][n]) < dftb->dacc)
        dftb1.b[m][n] = 0.e0;

  // the matrix b is correct here! //

  for (l=0; l<dftb1.nn; l++)  // for every atom that forces act upon
  for (j=0; j<dftb1.nn; j++) {
    indj = dftb1.ind[j];
    izpj = dftb1.izp[j];
    for (k=0; k<dftb1.nn; k++)   // only then there is a derivative
      if (k == l ){
      indk = dftb1.ind[k];
      izpk = dftb1.izp[k];
      for (i=0; i<3; i++) { // XX, YY and ZZ
        // derivative of the slko matrices
        xhelp = dftb1.x[j][i];
        x[k][i] = xhelp + deltax;
        slkmatrices(k, j, dftb1.x, dftb1.au,  dftb1.bu,  dftb->lmax, dftb->dim1, dftb->dr1, dftb1.izp, dftb->skstab1, dftb->skhtab1, dftb->skself1);
        x[k][i] = xhelp - deltax;
        slkmatrices(k, j, dftb1.x, dftb1.auh, dftb1.buh, dftb->lmax, dftb->dim1, dftb->dr1, dftb1.izp, dftb->skstab1, dftb->skhtab1, dftb->skself1);
        x[k][i] = xhelp;
        // use summation over angular momentum and magnetic quantum numbers
        // because shift is actually defined for the orbitals
        for (lj=1; lj<=dftb->lmax[izpj]; lj++)
          for (mj=1; mj<2*lj; mj++) {
            n  = SQR(lj-1) + mj - 1;
            nu = n + indj;
            for (lk=1; lk<=dftb->lmax[izpk]; lk++)
              for (mk=1; mk<2*lk; mk++) {
                m  = SQR(lk-1) + mk - 1;
                mu = m + indk;
                // dgrh = 2 * ( d H_{k,j}^{m,n}/ d R_{j} )
                //dgrh = (dftb1.au[m][n] - dftb1.auh[m][n]) / deltax;
                // dgrs = - 2 * ( d S_{k,j}^{m,n}/ d R_{j} )
                dgrs = -(dftb1.bu[m][n] - dftb1.buh[m][n]) / deltax;
                // dgr = ( d S_{k,j}^{m,n}/ d R_{j} ) * (shift(k)+shift(j))
                //dgr = -0.5 * dgrs * (dftb1.dshift[k] + dftb1.dshift[j]);
                dgr = -0.5 * dgrs * (dftb1.dshift[j]);


                // only lower triangle contains sum_i n(i) c_{mu,i} c_{nu,i}
                // only upper triangle contains sum_i epsilon_i n(i) c_{mu,i} c_{nu,i} //
                if (mu > nu) {
                  dgrs *= dftb1.b[nu][mu];
                  dgr  *= dftb1.b[mu][nu];
                } else {
                  //dgrh *= dftb1.b[nu][mu];
                  dgrs *= dftb1.b[mu][nu];
                  dgr  *= dftb1.b[nu][mu];
                }
                grad[l][i] +=  dgr;
printf("current additional grad on atom %d =  %f %f %f \n",l, grad[l][0],grad[l][1],grad[l][2]);
              }
          }
      }
    }
    else if (j == l) {
indk = dftb1.ind[k];
      izpk = dftb1.izp[k];
      for (i=0; i<3; i++) { // XX, YY and ZZ
        // derivative of the slko matrices
        xhelp = dftb1.x[j][i];
        x[j][i] = xhelp + deltax;
        slkmatrices(k, j, dftb1.x, dftb1.au,  dftb1.bu,  dftb->lmax, dftb->dim1, dftb->dr1, dftb1.izp, dftb->skstab1, dftb->skhtab1, dftb->skself1);
        x[j][i] = xhelp - deltax;
        slkmatrices(k, j, dftb1.x, dftb1.auh, dftb1.buh, dftb->lmax, dftb->dim1, dftb->dr1, dftb1.izp, dftb->skstab1, dftb->skhtab1, dftb->skself1);
        x[j][i] = xhelp;
        // use summation over angular momentum and magnetic quantum numbers
        // because shift is actually defined for the orbitals
        for (lj=1; lj<=dftb->lmax[izpj]; lj++)
          for (mj=1; mj<2*lj; mj++) {
            n  = SQR(lj-1) + mj - 1;
            nu = n + indj;
            for (lk=1; lk<=dftb->lmax[izpk]; lk++)
              for (mk=1; mk<2*lk; mk++) {
                m  = SQR(lk-1) + mk - 1;
                mu = m + indk;
                // dgrh = 2 * ( d H_{k,j}^{m,n}/ d R_{j} )
                //dgrh = (dftb1.au[m][n] - dftb1.auh[m][n]) / deltax;
                // dgrs = - 2 * ( d S_{k,j}^{m,n}/ d R_{j} )
                dgrs = -(dftb1.bu[m][n] - dftb1.buh[m][n]) / deltax;
                // dgr = ( d S_{k,j}^{m,n}/ d R_{j} ) * (shift(k)+shift(j))
                //dgr = -0.5 * dgrs * (dftb1.dshift[k] + dftb1.dshift[j]);
                dgr = -0.5 * dgrs * (dftb1.dshift[j]);


                // only lower triangle contains sum_i n(i) c_{mu,i} c_{nu,i}
                // only upper triangle contains sum_i epsilon_i n(i) c_{mu,i} c_{nu,i} //
                if (mu > nu) {
                  dgrs *= dftb1.b[nu][mu];
                  dgr  *= dftb1.b[mu][nu];
                } else {
                  //dgrh *= dftb1.b[nu][mu];
                  dgrs *= dftb1.b[mu][nu];
                  dgr  *= dftb1.b[nu][mu];
                }
                grad[l][i] +=  dgr;
printf("current additional grad on atom %d =  %f %f %f \n",l, grad[l][0],grad[l][1],grad[l][2]);
              }
          }
      }
    }
  }
  return;
}
