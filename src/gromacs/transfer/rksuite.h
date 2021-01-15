#ifndef RKSUITE_H
#define RKSUITE_H

//#include <cstdio>
//#include <cmath>
//#include <cstring>
//#include <cstdlib>
#include <cfloat>
#include "gromacs/transfer/transfer.h"

#define MAX_INTEGRATION_STEP (84.) /* 2 fs in atomic units */

#define max(a,b) ((a) >= (b) ? (a) : (b))
#define min(a,b) ((a) <= (b) ? (a) : (b))
//#define SQR(a) ((a)*(a))

void eqom(double time, double *var, double *dvar, int n, double **hamiltonian, double **hubbard,
          int nip, int *site_nip, double nip_rate_constant);
int do_rksuite(charge_transfer_t *ct);
int do_rksuite_diab(charge_transfer_t *ct);
void eqom_tfs(double time, double *var, double *dvar, int n, double *energy, double **overlap, double dt);
int do_rksuite_tfs(charge_transfer_t *ct);

void rk_setup(int neq, double tstart, double ystart[], double tend,
                    double tol, double thres[], int method, const char* task,
            bool errass, double hstart, double work[], int lenwrk,
            bool mesage);
void rk_ut(void (*f)(double, double*, double*, int, double**, double**, int, int*, double), double twant, double& tgot,
     double ygot[], double ypgot[], double ymax[], double work[], int& uflag, int sites, double **hamiltonian, double **hubbard,
     int nip, int *site_nip, double nip_rate_constant);
void rk_ut_tfs(void (*f)(double, double*, double*, int, double*, double**, double), double twant, double& tgot,
     double ygot[], double ypgot[], double ymax[], double work[], int& uflag, int sites, double *energy, double **overlap, double dt);
void stat(int& totfcn, int& stpcst, double& waste, int& stpsok, double& hnext);
void glberr(double rmserr[], double& errmax, double& terrmx, double work[]);
void rk_ct(void (*f)(double, double*, double*, int, double**, double**, int, int*, double), double& tnow, double ynow[],
     double ypnow[], double work[],int& cflag, int sites, double **hamiltonian, double **hubbard,
     int nip, int *site_nip, double nip_rate_constant);
void rk_ct_tfs(void (*f)(double, double*, double*, int, double*, double**, double), double& tnow, double ynow[],
     double ypnow[], double work[],int& cflag, int sites, double *energy, double **overlap, double dt);
void intrp(double twant, const char reqest[], int nwant, double ywant[], double ypwant[],
     void (*f)(double, double*, double*, int, double**, double**, int, int*, double), double work[], double wrkint[],
     int lenint, int sites, double **hamiltonian, double **hubbard, int nip, int *site_nip, double nip_rate_constant);
void intrp_tfs(double twant, const char reqest[], int nwant, double ywant[], double ypwant[],
     void (*f)(double, double*, double*, int, double*, double**, double), double work[], double wrkint[],
     int lenint, int sites, double *energy, double **overlap, double dt);
void reset(double tendnu);
void mconst(int method);
void envirn(int& outch, double& mcheps, double& dwarf);
void rkconst(int method, int& vecstg, bool& reqstg, int& lintpl);
void rkmsg(int ier, const char* srname, int nrec, int& flag);
void rksit(bool ask, const char* srname, int& state);
void step(void (*f)(double, double*, double*, int, double**, double**, int, int*, double), int neq, double tnow,
	double* y, double* yp, double stages /* [neq][] */ [], double tol, double& htry,
	double* weight, double* ynew, double* errest, double& err, bool main,
	double hmin, double* thres, bool& phase2, int sites, double **hamiltonian, double **hubbard,
        int nip, int *site_nip, double nip_rate_constant);
void stepa(double tnow, double y[], double yp[], double tstg, double ystg[],
	double ypstg[], double& htry, double weight[], bool& cutbak);
void step_tfs(void (*f)(double, double*, double*, int, double*, double**, double), int neq, double tnow,
	double* y, double* yp, double stages /* [neq][] */ [], double tol, double& htry,
	double* weight, double* ynew, double* errest, double& err, bool main,
	double hmin, double* thres, bool& phase2, int sites, double *energy, double **overlap, double dt);
void stepb(int neq, double y[], double yp[], double h, double ynew[],
	double stages/*[neq][]*/[], double thres[], double& err, bool main,
	double weight[]);
void stiff(void (*f)(double, double*, double*, int, double**, double**, int, int*, double), double havg, int& jflstp,
	bool toomch, int maxfcn, double work[], int& ier, int& nrec, int sites, double **hamiltonian, double **hubbard,
        int nip, int *site_nip, double nip_rate_constant);
void stiffa(void (*f)(double, double*, double*, int, double**, double**, int, int*, double),double x, double y[],
	double hnow, double havg, double xend, int maxfcn, double wt[],
	double fxy[], double v0[], bool& unsure, bool& stif, double v1[],
	double v2[], double v3[], double vtemp[], int sites, double **hamiltonian, double **hubbard,
        int nip, int *site_nip, double nip_rate_constant);
void stiff_tfs(void (*f)(double, double*, double*, int, double*, double**, double), double havg, int& jflstp,
	bool toomch, int maxfcn, double work[], int& ier, int& nrec, int sites, double *energy, double **overlap, double dt);
void stiffa_tfs(void (*f)(double, double*, double*, int, double*, double**, double),double x, double y[],
	double hnow, double havg, double xend, int maxfcn, double wt[],
	double fxy[], double v0[], bool& unsure, bool& stif, double v1[],
	double v2[], double v3[], double vtemp[], int sites, double *energy, double **overlap, double dt);
void stiffb(double v1v1, double v0v1, double v0v0, double& rold,
	double& rho, double root1[], double root2[], bool& rootre);
void stiffc(double alpha, double beta, double r1[], double r2[]);
void stiffd(double v[], double havg, double x, double y[],
	void (*f)(double, double[], double[], int, double**, double**, int, int*, double), double fxy[], double wt[],
	double scale, double vdotv, double z[], double& zdotz, double vtemp[], int sites, double **hamiltonian, double **hubbard,
        int nip, int *site_nip, double nip_rate_constant);
void stiffd_tfs(double v[], double havg, double x, double y[],
	void (*f)(double, double[], double[], int, double*, double**, double), double fxy[], double wt[],
	double scale, double vdotv, double z[], double& zdotz, double vtemp[], int sites, double *energy, double **overlap, double dt);
double dotprd(double u[], double v[], double wt[], int neq);
void softfl(bool ask, bool* on);
void chkfl(bool ask, bool& error);
void evali(double y[], double yp[], double p /* [nwant][] */ [], double twant,
	const char reqest[], int nwant, double ywant[], double ypwant[]);
void formi(void (*f)(double, double*, double*, int, double**, double**, int, int*, double), int neq, int nwant, double y[], double yp[],
     double yold[], double ypold[], double stages/*[neq][]*/[], bool calstg,
     double xstage[], double p[], int sites, double **hamiltonian, double **hubbard,
     int nip, int *site_nip, double nip_rate_constant);
void truerr(void (*f)(double, double*, double*, int, double**, double**, int, int*, double), int neq, double y[],
     double tol, double weight[], double zy[], double zyp[], double zerror[],
     double zynew[], double zerres[], double zstage /* [neq][] */ [], int& ier, int sites, double **hamiltonian, double **hubbard,
     int nip, int *site_nip, double nip_rate_constant);
void formi_tfs(void (*f)(double, double*, double*, int, double*, double**, double), int neq, int nwant, double y[], double yp[],
     double yold[], double ypold[], double stages/*[neq][]*/[], bool calstg,
     double xstage[], double p[], int sites, double *energy, double **overlap, double dt);
void truerr_tfs(void (*f)(double, double*, double*, int, double*, double**, double), int neq, double y[],
     double tol, double weight[], double zy[], double zyp[], double zerror[],
     double zynew[], double zerres[], double zstage /* [neq][] */ [], int& ier, int sites, double *energy, double **overlap, double dt);
#endif

