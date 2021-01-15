#include "gromacs/transfer/transfer.h"

void write_out_MOs(int step, rvec *x_ct, t_atoms *ct_atoms, dftb_t *dftb, charge_transfer_t *ct);
void write_out_MOs(int step, rvec *x_ct, t_atoms *ct_atoms, dftb_t *dftb, charge_transfer_t *ct)
{
 // double amp, phase;
    int    i, j, counter, k;
 // int    lj, mj, jofn, jofn_old, mo = 1;
 // char   c;
    char   filename[20];
 // FILE  *f_ct_step_pdb = NULL;

    (void) dftb;

    sprintf(filename, "%d", step);
    strcat(filename, ".pdb" );
    //f_ct_step_pdb=fopen(filename,"w");
/*
   phase=1;
   for (j=0; j < dftb->phase2.nn; j++) {
    switch (dftb->phase2.izp[j]) {
      case 0: c = 'C'; break;
      case 1: c = 'H'; break;
      case 2: c = 'N'; break;
      case 3: c = 'O'; break;
    }
    amp = 0.0;
    if (dftb->lmax[dftb->phase2.izp[j]] == 2) {
      lj=1; //only p-orbitals so far
      jofn = dftb->phase2.ind[j] + lj*lj;
      for (mj=0; mj<=2*lj; mj++) {
        amp +=  SQR(dftb->orthogo.evec_ao[jofn+mj][mo]);
      }
      amp= sqrt(amp)*100; //for better reading
      if(j!=0){
        phase = 0;
        for (mj=0; mj<=2*lj; mj++)
          phase +=  dftb->orthogo.evec_ao[jofn+mj][mo]*dftb->orthogo.evec_ao[jofn_old+mj][mo];
        phase=(phase > 0) ? 1 : ((phase < 0) ? -1 : 0);
      }
      jofn_old = jofn;
      fprintf(f_ct_step_pdb,"%6s%5i   %c %s %s%4s    %8.3f%8.3f%8.3f%6s%6.2f\n" ,"ATOM  ",j+1, c ,"res","I","  1 ",dftb->phase2.x[j][0]/NM_TO_BOHR*10,dftb->phase2.x[j][1]/NM_TO_BOHR*10,dftb->phase2.x[j][2]/NM_TO_BOHR*10,"  1.00",phase*amp);//pdb format
    }
    else{
      fprintf(f_ct_step_pdb,"%6s%5i   %c %s %s%4s    %8.3f%8.3f%8.3f%6s%6s\n" \
      ,"ATOM  ",j+1,c,"res","I","  1 ",dftb->phase2.x[j][0]/NM_TO_BOHR*10,dftb->phase2.x[j][1]/NM_TO_BOHR*10,dftb->phase2.x[j][2]/NM_TO_BOHR*10,"  1.00","  0.00");//pdb format
    }
   }
 */
// write interesting data in pbd file as bfactor
    counter = 0;
    for (i = 0; i < ct->sites; i++)
        for (j = 0; j < ct->site[i].atoms; j++)
        {
            ct_atoms->pdbinfo[ct->site[i].atom[j]].bfac = 0.0;
            for (k = 0; k < ct->site[i].homos; k++)
            {
                //ct_atoms->pdbinfo[ct->site[i].atom[j]].bfac =(real) (dftb->phase2.qmat[counter] - dftb->qzero2[dftb->phase2.izp[counter]]);
                ct_atoms->pdbinfo[ct->site[i].atom[j]].bfac += (real) 100*(ct->occupation[ct->indFO[i]+k]);
                counter++;
            }
        }

    write_sto_conf(filename, "written by charge transfer code", ct_atoms, x_ct, NULL, PbcType::No, NULL);

    //fclose(f_ct_step_pdb);

    return;
}


/* This is already implemented in the standard QM/MM:
void print_time_difference(const char *s, struct timespec start, struct timespec end)
{
 // int       sec, nsec;
    long long value = 0ll;

    value = 1000000000ll * ((long long) end.tv_sec - (long long) start.tv_sec) + (long long) (end.tv_nsec - start.tv_nsec);
    printf("%s %12lld\n", s, value);

    return;
}
*/

