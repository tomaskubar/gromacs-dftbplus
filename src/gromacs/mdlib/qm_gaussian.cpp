/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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

#include "qm_gaussian.h"

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cmath>

#include "gromacs/fileio/confio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdlib/forcerec.h"
#include "gromacs/mdlib/qmmm.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

// When not built in a configuration with QMMM support, much of this
// code is unreachable by design. Tell clang not to warn about it.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-noreturn"

/* TODO: this should be made thread-safe */

/* Gaussian interface routines */

//void QMMM_QMgaussian::init_gaussian(QMMM_QMrec& qm)
void QMMM_QMgaussian::init_gaussian(QMMM_QMsurfaceHopping &surfaceHopping,
                                    int QMbasis)
{
    ivec
        basissets[eQMbasisNR] = {{0, 3, 0},
                                 {0, 3, 0}, /*added for double sto-3g entry in names.c*/
                                 {5, 0, 0},
                                 {5, 0, 1},
                                 {5, 0, 11},
                                 {5, 6, 0},
                                 {1, 6, 0},
                                 {1, 6, 1},
                                 {1, 6, 11},
                                 {4, 6, 0}};
    char
       *buf = nullptr;
    int
        i;

    if (!GMX_QMMM_GAUSSIAN)
    {
        gmx_fatal(FARGS, "Cannot call GAUSSIAN unless linked against it. Use cmake -DGMX_QMMM_PROGRAM=GAUSSIAN, and ensure that linking will work correctly.");
    }

    /* using the ivec above to convert the basis read form the mdp file
     * in a human readable format into some numbers for the gaussian
     * route. This is necessary as we are using non standard routes to
     * do SH.
     */

    /* per layer we make a new subdir for integral file, checkpoint
     * files and such. These dirs are stored in the QMrec for
     * convenience
     */


    if (!nQMcpus) /* this we do only once per layer
                       * as we call g01 externally
                       */

    {
        for (i = 0; i < DIM; i++)
        {
            surfaceHopping.SHbasis[i] = basissets[QMbasis][i];
        }

        /* init gradually switching on of the SA */
        surfaceHopping.SAstep = 0;
        /* we read the number of cpus and environment from the environment
         * if set.
         */
        buf = getenv("GMX_QM_GAUSSIAN_NCPUS");
        if (buf)
        {
            sscanf(buf, "%d", &nQMcpus);
        }
        else
        {
            nQMcpus = 1;
        }
        fprintf(stderr, "number of CPUs for gaussian = %d\n", nQMcpus);
        buf = getenv("GMX_QM_GAUSSIAN_MEMORY");
        if (buf)
        {
            sscanf(buf, "%d", &QMmem);
        }
        else
        {
            QMmem = 50000000;
        }
        fprintf(stderr, "memory for gaussian = %d\n", QMmem);
        buf = getenv("GMX_QM_ACCURACY");
        if (buf)
        {
            sscanf(buf, "%d", &accuracy);
        }
        else
        {
            accuracy = 8;
        }
        fprintf(stderr, "accuracy in l510 = %d\n", accuracy);

        buf = getenv("GMX_QM_CPMCSCF");
        if (buf)
        {
            sscanf(buf, "%d", &i);
            cpmcscf = (i != 0);
        }
        else
        {
            cpmcscf = FALSE;
        }
        if (cpmcscf)
        {
            fprintf(stderr, "using cp-mcscf in l1003\n");
        }
        else
        {
            fprintf(stderr, "NOT using cp-mcscf in l1003\n");
        }
        buf = getenv("GMX_QM_SA_STEP");
        if (buf)
        {
            sscanf(buf, "%d", &surfaceHopping.SAstep);
        }
        else
        {
            /* init gradually switching on of the SA */
            surfaceHopping.SAstep = 0;
        }
        /* we read the number of cpus and environment from the environment
         * if set.
         */
        fprintf(stderr, "Level of SA at start = %d\n", surfaceHopping.SAstep);
        /* gaussian settings on the system */
        buf = getenv("GMX_QM_GAUSS_DIR");

        if (buf)
        {
            gauss_dir = gmx_strdup(buf);
        }
        else
        {
            gmx_fatal(FARGS, "no $GMX_QM_GAUSS_DIR, check gaussian manual\n");
        }

        buf = getenv("GMX_QM_GAUSS_EXE");
        if (buf)
        {
            gauss_exe = gmx_strdup(buf);
        }
        else
        {
            gmx_fatal(FARGS, "no $GMX_QM_GAUSS_EXE set, check gaussian manual\n");
        }
        buf = getenv("GMX_QM_MODIFIED_LINKS_DIR");
        if (buf)
        {
            devel_dir = gmx_strdup(buf);
        }
        else
        {
            gmx_fatal(FARGS, "no $GMX_QM_MODIFIED_LINKS_DIR, this is were the modified links reside.\n");
        }

        /*  if(fr->bRF){*/
        /* reactionfield, file is needed using gaussian */
        /*    rffile=fopen("rf.dat","w");*/
        /*   fprintf(rffile,"%f %f\n",fr->epsilon_r,fr->rcoulomb/BOHR2NM);*/
        /* fclose(rffile);*/
        /*  }*/
    }
    fprintf(stderr, "gaussian initialised...\n");
}


static void write_gaussian_SH_input(int step, gmx_bool swap,
                                    gmx_unused const t_forcerec *fr, QMMM_QMrec& qm, QMMM_MMrec& mm)
{
    int
        i;
    gmx_bool
        bSA;
    FILE
       *out;
 // t_QMMMrec
 //    *QMMMrec;
 // QMMMrec = fr->qr;
    bSA     = (qm.surfaceHopping.SAstep > 0);

    out = fopen("input.com", "w");
    /* write the route */
    fprintf(out, "%s", "%scr=input\n");
    fprintf(out, "%s", "%rwf=input\n");
    fprintf(out, "%s", "%int=input\n");
    fprintf(out, "%s", "%d2e=input\n");
/*  if(step)
 *   fprintf(out,"%s","%nosave\n");
 */
    fprintf(out, "%s", "%chk=input\n");
    fprintf(out, "%s%d\n", "%mem=", qm.gaussian.QMmem);
    fprintf(out, "%s%3d\n", "%nprocshare=", qm.gaussian.nQMcpus);

    /* use the versions of
     * l301 that computes the interaction between MM and QM atoms.
     * l510 that can punch the CI coefficients
     * l701 that can do gradients on MM atoms
     */

    /* local version */
    fprintf(out, "%s%s%s",
            "%subst l510 ",
            qm.gaussian.devel_dir,
            "/l510\n");
    fprintf(out, "%s%s%s",
            "%subst l301 ",
            qm.gaussian.devel_dir,
            "/l301\n");
    fprintf(out, "%s%s%s",
            "%subst l701 ",
            qm.gaussian.devel_dir,
            "/l701\n");

    fprintf(out, "%s%s%s",
            "%subst l1003 ",
            qm.gaussian.devel_dir,
            "/l1003\n");
    fprintf(out, "%s%s%s",
            "%subst l9999 ",
            qm.gaussian.devel_dir,
            "/l9999\n");
    /* print the nonstandard route
     */
    fprintf(out, "%s",
            "#P nonstd\n 1/18=10,20=1,38=1/1;\n");
    fprintf(out, "%s",
            " 2/9=110,15=1,17=6,18=5,40=1/2;\n");
    if (mm.nrMMatoms)
    {
        fprintf(out,
                " 3/5=%d,6=%d,7=%d,25=1,32=1,43=1,94=-2/1,2,3;\n",
                qm.surfaceHopping.SHbasis[0],
                qm.surfaceHopping.SHbasis[1],
                qm.surfaceHopping.SHbasis[2]); /*basisset stuff */
    }
    else
    {
        fprintf(out,
                " 3/5=%d,6=%d,7=%d,25=1,32=1,43=0,94=-2/1,2,3;\n",
                qm.surfaceHopping.SHbasis[0],
                qm.surfaceHopping.SHbasis[1],
                qm.surfaceHopping.SHbasis[2]); /*basisset stuff */
    }
    /* development */
    if (step+1) /* fetch initial guess from check point file */
    {           /* hack, to alyays read from chk file!!!!! */
        fprintf(out, "%s%d,%s%d%s", " 4/5=1,7=6,17=",
                qm.surfaceHopping.CASelectrons,
                "18=", qm.surfaceHopping.CASorbitals, "/1,5;\n");
    }
    else /* generate the first checkpoint file */
    {
        fprintf(out, "%s%d,%s%d%s", " 4/5=0,7=6,17=",
                qm.surfaceHopping.CASelectrons,
                "18=", qm.surfaceHopping.CASorbitals, "/1,5;\n");
    }
    /* the rest of the input depends on where the system is on the PES
     */
    if (swap && bSA)             /* make a slide to the other surface */
    {
        if (qm.surfaceHopping.CASorbitals > 6) /* use direct and no full diag */
        {
            fprintf(out, " 5/5=2,16=-2,17=10000000,28=2,32=2,38=6,97=100/10;\n");
        }
        else
        {
            if (qm.gaussian.cpmcscf)
            {
                fprintf(out, " 5/5=2,6=%d,17=31000200,28=2,32=2,38=6,97=100/10;\n",
                        qm.gaussian.accuracy);
                if (mm.nrMMatoms > 0)
                {
                    fprintf(out, " 7/7=1,16=-2,30=1/1;\n");
                }
                fprintf(out, " 11/31=1,42=1,45=1/1;\n");
                fprintf(out, " 10/6=1,10=700006,28=2,29=1,31=1,97=100/3;\n");
                fprintf(out, " 7/30=1/16;\n 99/10=4/99;\n");
            }
            else
            {
                fprintf(out, " 5/5=2,6=%d,17=11000000,28=2,32=2,38=6,97=100/10;\n",
                        qm.gaussian.accuracy);
                fprintf(out, " 7/7=1,16=-2,30=1/1,2,3,16;\n 99/10=4/99;\n");
            }
        }
    }
    else if (bSA)                /* do a "state-averaged" CAS calculation */
    {
        if (qm.surfaceHopping.CASorbitals > 6) /* no full diag */
        {
            fprintf(out, " 5/5=2,16=-2,17=10000000,28=2,32=2,38=6/10;\n");
        }
        else
        {
            if (qm.gaussian.cpmcscf)
            {
                fprintf(out, " 5/5=2,6=%d,17=31000200,28=2,32=2,38=6/10;\n",
                        qm.gaussian.accuracy);
                if (mm.nrMMatoms > 0)
                {
                    fprintf(out, " 7/7=1,16=-2,30=1/1;\n");
                }
                fprintf(out, " 11/31=1,42=1,45=1/1;\n");
                fprintf(out, " 10/6=1,10=700006,28=2,29=1,31=1/3;\n");
                fprintf(out, " 7/30=1/16;\n 99/10=4/99;\n");
            }
            else
            {
                fprintf(out, " 5/5=2,6=%d,17=11000000,28=2,32=2,38=6/10;\n",
                        qm.gaussian.accuracy);
                fprintf(out, " 7/7=1,16=-2,30=1/1,2,3,16;\n 99/10=4/99;\n");
            }
        }
    }
    else if (swap) /* do a "swapped" CAS calculation */
    {
        if (qm.surfaceHopping.CASorbitals > 6)
        {
            fprintf(out, " 5/5=2,16=-2,17=0,28=2,32=2,38=6,97=100/10;\n");
        }
        else
        {
            fprintf(out, " 5/5=2,6=%d,17=1000000,28=2,32=2,38=6,97=100/10;\n",
                    qm.gaussian.accuracy);
        }
        fprintf(out, " 7/7=1,16=-2,30=1/1,2,3,16;\n 99/10=4/99;\n");
    }
    else /* do a "normal" CAS calculation */
    {
        if (qm.surfaceHopping.CASorbitals > 6)
        {
            fprintf(out, " 5/5=2,16=-2,17=0,28=2,32=2,38=6/10;\n");
        }
        else
        {
            fprintf(out, " 5/5=2,6=%d,17=1000000,28=2,32=2,38=6/10;\n",
                    qm.gaussian.accuracy);
        }
        fprintf(out, " 7/7=1,16=-2,30=1/1,2,3,16;\n 99/10=4/99;\n");
    }
    fprintf(out, "\ninput-file generated by gromacs\n\n");
    fprintf(out, "%2d%2d\n", qm.QMcharge_get(), qm.multiplicity_get());
    for (i = 0; i < qm.nrQMatoms_get(); i++)
    {
        fprintf(out, "%3d %10.7f  %10.7f  %10.7f\n",
                qm.atomicnumberQM_get(i),
                qm.xQM_get(i,XX)/BOHR2NM,
                qm.xQM_get(i,YY)/BOHR2NM,
                qm.xQM_get(i,ZZ)/BOHR2NM);
    }
    /* MM point charge data */
    if (mm.nrMMatoms)
    {
        fprintf(out, "\n");
        for (i = 0; i < mm.nrMMatoms; i++)
        {
            fprintf(out, "%10.7f  %10.7f  %10.7f %8.4f\n",
                    mm.xMM[i][XX]/BOHR2NM,
                    mm.xMM[i][YY]/BOHR2NM,
                    mm.xMM[i][ZZ]/BOHR2NM,
                    mm.MMcharges[i]);
        }
    }
    if (bSA) /* put the SA coefficients at the end of the file */
    {
        fprintf(out, "\n%10.8f %10.8f\n",
                qm.surfaceHopping.SAstep*0.5/qm.surfaceHopping.SAsteps,
                1-qm.surfaceHopping.SAstep*0.5/qm.surfaceHopping.SAsteps);
        fprintf(stderr, "State Averaging level = %d/%d\n", qm.surfaceHopping.SAstep, qm.surfaceHopping.SAsteps);
    }
    fprintf(out, "\n");
    fclose(out);
}  /* write_gaussian_SH_input */

static void write_gaussian_input(int step, gmx_unused const t_forcerec *fr, QMMM_QMrec& qm, QMMM_MMrec& mm)
{
    int
        i;
 // t_QMMMrec
 //    *QMMMrec;
 // QMMMrec = fr->qr;
    FILE
       *out;

    out     = fopen("input.com", "w");
    /* write the route */

    if (qm.QMmethod_get() >= eQMmethodRHF)
    {
        fprintf(out, "%s",
                "%chk=input\n");
    }
    else
    {
        fprintf(out, "%s",
                "%chk=se\n");
    }
    if (qm.gaussian.nQMcpus > 1)
    {
        fprintf(out, "%s%3d\n",
                "%nprocshare=", qm.gaussian.nQMcpus);
    }
    fprintf(out, "%s%d\n",
            "%mem=", qm.gaussian.QMmem);
    fprintf(out, "%s%s%s",
            "%subst l701 ", qm.gaussian.devel_dir, "/l701\n");
    fprintf(out, "%s%s%s",
            "%subst l301 ", qm.gaussian.devel_dir, "/l301\n");
    fprintf(out, "%s%s%s",
            "%subst l9999 ", qm.gaussian.devel_dir, "/l9999\n");
    if (step)
    {
        fprintf(out, "%s",
                "#T ");
    }
    else
    {
        fprintf(out, "%s",
                "#P ");
    }
    if (qm.QMmethod_get() == eQMmethodB3LYPLAN)
    {
        fprintf(out, " %s",
                "B3LYP/GEN Pseudo=Read");
    }
    else
    {
        fprintf(out, " %s",
                eQMmethod_names[qm.QMmethod_get()]);

        if (qm.QMmethod_get() >= eQMmethodRHF)
        {
            if (qm.QMmethod_get() == eQMmethodCASSCF)
            {
                /* in case of cas, how many electrons and orbitals do we need?
                 */
                fprintf(out, "(%d,%d)",
                        qm.surfaceHopping.CASelectrons, qm.surfaceHopping.CASorbitals);
            }
            fprintf(out, "/%s",
                    eQMbasis_names[qm.QMbasis_get()]);
        }
    }
    if (mm.nrMMatoms)
    {
        fprintf(out, " %s",
                "Charge ");
    }
    if (step || qm.QMmethod_get() == eQMmethodCASSCF)
    {
        /* fetch guess from checkpoint file, always for CASSCF */
        fprintf(out, "%s", " guess=read");
    }
    fprintf(out, "\nNosymm units=bohr\n");

    fprintf(out, "FORCE Punch=(Derivatives) ");
    fprintf(out, "iop(3/33=1)\n\n");
    fprintf(out, "input-file generated by gromacs\n\n");
    fprintf(out, "%2d%2d\n", qm.QMcharge_get(), qm.multiplicity_get());
    for (i = 0; i < qm.nrQMatoms_get(); i++)
    {
        fprintf(out, "%3d %10.7f  %10.7f  %10.7f\n",
                qm.atomicnumberQM_get(i),
                qm.xQM_get(i,XX)/BOHR2NM,
                qm.xQM_get(i,YY)/BOHR2NM,
                qm.xQM_get(i,ZZ)/BOHR2NM);
    }

    /* Pseudo Potential and ECP are included here if selected (MEthod suffix LAN) */
    if (qm.QMmethod_get() == eQMmethodB3LYPLAN)
    {
        fprintf(out, "\n");
        for (i = 0; i < qm.nrQMatoms_get(); i++)
        {
            if (qm.atomicnumberQM_get(i) < 21)
            {
                fprintf(out, "%d ", i+1);
            }
        }
        fprintf(out, "\n%s\n****\n", eQMbasis_names[qm.QMbasis_get()]);

        for (i = 0; i < qm.nrQMatoms_get(); i++)
        {
            if (qm.atomicnumberQM_get(i) > 21)
            {
                fprintf(out, "%d ", i+1);
            }
        }
        fprintf(out, "\n%s\n****\n\n", "lanl2dz");

        for (i = 0; i < qm.nrQMatoms_get(); i++)
        {
            if (qm.atomicnumberQM_get(i) > 21)
            {
                fprintf(out, "%d ", i+1);
            }
        }
        fprintf(out, "\n%s\n", "lanl2dz");
    }



    /* MM point charge data */
    if (mm.nrMMatoms)
    {
        fprintf(stderr, "nr mm atoms in gaussian.c = %d\n", mm.nrMMatoms);
        fprintf(out, "\n");
        for (i = 0; i < mm.nrMMatoms; i++)
        {
            fprintf(out, "%10.7f  %10.7f  %10.7f %8.4f\n",
                    mm.xMM[i][XX]/BOHR2NM,
                    mm.xMM[i][YY]/BOHR2NM,
                    mm.xMM[i][ZZ]/BOHR2NM,
                    mm.MMcharges[i]);
        }
    }
    fprintf(out, "\n");


    fclose(out);

}  /* write_gaussian_input */

static real read_gaussian_output(rvec QMgrad[], rvec MMgrad[], QMMM_QMrec& qm, QMMM_MMrec& mm)
{
    int
        i;
    char
        buf[300];
    real
        QMener;
    FILE
       *in;

    in = fopen("fort.7", "r");

    /* (There was additional content in the file in case
     *    of QM optimizations / transition state search,
     *    which was removed.
     */
    /* the next line is the energy and in the case of CAS, the energy
     * difference between the two states.
     */
    if (nullptr == fgets(buf, 300, in))
    {
        gmx_fatal(FARGS, "Error reading Gaussian output");
    }

#if GMX_DOUBLE
    sscanf(buf, "%lf\n", &QMener);
#else
    sscanf(buf, "%f\n", &QMener);
#endif
    /* next lines contain the gradients of the QM atoms */
    for (i = 0; i < qm.nrQMatoms_get(); i++)
    {
        if (nullptr == fgets(buf, 300, in))
        {
            gmx_fatal(FARGS, "Error reading Gaussian output");
        }
#if GMX_DOUBLE
        sscanf(buf, "%lf %lf %lf\n",
               &QMgrad[i][XX],
               &QMgrad[i][YY],
               &QMgrad[i][ZZ]);
#else
        sscanf(buf, "%f %f %f\n",
               &QMgrad[i][XX],
               &QMgrad[i][YY],
               &QMgrad[i][ZZ]);
#endif
    }
    /* the next lines are the gradients of the MM atoms */
    if (qm.QMmethod_get() >= eQMmethodRHF)
    {
        for (i = 0; i < mm.nrMMatoms; i++)
        {
            if (nullptr == fgets(buf, 300, in))
            {
                gmx_fatal(FARGS, "Error reading Gaussian output");
            }
#if GMX_DOUBLE
            sscanf(buf, "%lf %lf %lf\n",
                   &MMgrad[i][XX],
                   &MMgrad[i][YY],
                   &MMgrad[i][ZZ]);
#else
            sscanf(buf, "%f %f %f\n",
                   &MMgrad[i][XX],
                   &MMgrad[i][YY],
                   &MMgrad[i][ZZ]);
#endif
        }
    }
    fclose(in);
    return(QMener);
}

static real read_gaussian_SH_output(rvec QMgrad[], rvec MMgrad[], int step, QMMM_QMrec& qm, QMMM_MMrec& mm)
{
    int
        i;
    char
        buf[300];
    real
        QMener, DeltaE;
    FILE
       *in;

    in = fopen("fort.7", "r");
    /* first line is the energy and in the case of CAS, the energy
     * difference between the two states.
     */
    if (nullptr == fgets(buf, 300, in))
    {
        gmx_fatal(FARGS, "Error reading Gaussian output");
    }

#if GMX_DOUBLE
    sscanf(buf, "%lf %lf\n", &QMener, &DeltaE);
#else
    sscanf(buf, "%f %f\n",  &QMener, &DeltaE);
#endif

    /* switch on/off the State Averaging */

    if (DeltaE > qm.surfaceHopping.SAoff)
    {
        if (qm.surfaceHopping.SAstep > 0)
        {
            qm.surfaceHopping.SAstep--;
        }
    }
    else if (DeltaE < qm.surfaceHopping.SAon || (qm.surfaceHopping.SAstep > 0))
    {
        if (qm.surfaceHopping.SAstep < qm.surfaceHopping.SAsteps)
        {
            qm.surfaceHopping.SAstep++;
        }
    }

    /* for debugging: */
    fprintf(stderr, "Gap = %5f,SA = %3d\n", DeltaE, static_cast<int>(qm.surfaceHopping.SAstep > 0));
    /* next lines contain the gradients of the QM atoms */
    for (i = 0; i < qm.nrQMatoms_get(); i++)
    {
        if (nullptr == fgets(buf, 300, in))
        {
            gmx_fatal(FARGS, "Error reading Gaussian output");
        }

#if GMX_DOUBLE
        sscanf(buf, "%lf %lf %lf\n",
               &QMgrad[i][XX],
               &QMgrad[i][YY],
               &QMgrad[i][ZZ]);
#else
        sscanf(buf, "%f %f %f\n",
               &QMgrad[i][XX],
               &QMgrad[i][YY],
               &QMgrad[i][ZZ]);
#endif
    }
    /* the next lines, are the gradients of the MM atoms */

    for (i = 0; i < mm.nrMMatoms; i++)
    {
        if (nullptr == fgets(buf, 300, in))
        {
            gmx_fatal(FARGS, "Error reading Gaussian output");
        }
#if GMX_DOUBLE
        sscanf(buf, "%lf %lf %lf\n",
               &MMgrad[i][XX],
               &MMgrad[i][YY],
               &MMgrad[i][ZZ]);
#else
        sscanf(buf, "%f %f %f\n",
               &MMgrad[i][XX],
               &MMgrad[i][YY],
               &MMgrad[i][ZZ]);
#endif
    }

    /* the next line contains the two CI eigenvector elements */
    if (nullptr == fgets(buf, 300, in))
    {
        gmx_fatal(FARGS, "Error reading Gaussian output");
    }
    if (!step)
    {
        sscanf(buf, "%d", &qm.surfaceHopping.CIdim);
        snew(qm.surfaceHopping.CIvec1, qm.surfaceHopping.CIdim);
        snew(qm.surfaceHopping.CIvec1old, qm.surfaceHopping.CIdim);
        snew(qm.surfaceHopping.CIvec2, qm.surfaceHopping.CIdim);
        snew(qm.surfaceHopping.CIvec2old, qm.surfaceHopping.CIdim);
    }
    else
    {
        /* before reading in the new current CI vectors, copy the current
         * CI vector into the old one.
         */
        for (i = 0; i < qm.surfaceHopping.CIdim; i++)
        {
            qm.surfaceHopping.CIvec1old[i] = qm.surfaceHopping.CIvec1[i];
            qm.surfaceHopping.CIvec2old[i] = qm.surfaceHopping.CIvec2[i];
        }
    }
    /* first vector */
    for (i = 0; i < qm.surfaceHopping.CIdim; i++)
    {
        if (nullptr == fgets(buf, 300, in))
        {
            gmx_fatal(FARGS, "Error reading Gaussian output");
        }
#if GMX_DOUBLE
        sscanf(buf, "%lf\n", &qm.surfaceHopping.CIvec1[i]);
#else
        sscanf(buf, "%f\n", &qm.surfaceHopping.CIvec1[i]);
#endif
    }
    /* second vector */
    for (i = 0; i < qm.surfaceHopping.CIdim; i++)
    {
        if (nullptr == fgets(buf, 300, in))
        {
            gmx_fatal(FARGS, "Error reading Gaussian output");
        }
#if GMX_DOUBLE
        sscanf(buf, "%lf\n", &qm.surfaceHopping.CIvec2[i]);
#else
        sscanf(buf, "%f\n", &qm.surfaceHopping.CIvec2[i]);
#endif
    }
    fclose(in);
    return(QMener);
}

static real inproduct(const real *a, const real *b, int n)
{
    int
        i;
    real
        dot = 0.0;

    /* computes the inner product between two vectors (a.b), both of
     * which have length n.
     */
    for (i = 0; i < n; i++)
    {
        dot += a[i]*b[i];
    }
    return(dot);
}

static int hop(int step, QMMM_QMrec& qm)
{
    int
        swap = 0;
    real
        d11 = 0.0, d12 = 0.0, d21 = 0.0, d22 = 0.0;

    /* calculates the inproduct between the current Ci vector and the
     * previous CI vector. A diabatic hop will be made if d12 and d21
     * are much bigger than d11 and d22. In that case hop returns true,
     * otherwise it returns false.
     */
    if (step) /* only go on if more than one step has been done */
    {
        d11 = inproduct(qm.surfaceHopping.CIvec1, qm.surfaceHopping.CIvec1old, qm.surfaceHopping.CIdim);
        d12 = inproduct(qm.surfaceHopping.CIvec1, qm.surfaceHopping.CIvec2old, qm.surfaceHopping.CIdim);
        d21 = inproduct(qm.surfaceHopping.CIvec2, qm.surfaceHopping.CIvec1old, qm.surfaceHopping.CIdim);
        d22 = inproduct(qm.surfaceHopping.CIvec2, qm.surfaceHopping.CIvec2old, qm.surfaceHopping.CIdim);
    }
    fprintf(stderr, "-------------------\n");
    fprintf(stderr, "d11 = %13.8f\n", d11);
    fprintf(stderr, "d12 = %13.8f\n", d12);
    fprintf(stderr, "d21 = %13.8f\n", d21);
    fprintf(stderr, "d22 = %13.8f\n", d22);
    fprintf(stderr, "-------------------\n");

    if ((fabs(d12) > 0.5) && (fabs(d21) > 0.5))
    {
        swap = 1;
    }

    return(swap);
}

static void do_gaussian(int step, char *exe)
{
    char
        buf[STRLEN];

    /* make the call to the gaussian binary through system()
     * The location of the binary will be picked up from the
     * environment using getenv().
     */
    if (step) /* hack to prevent long inputfiles */
    {
        sprintf(buf, "%s < %s > %s",
                exe,
                "input.com",
                "input.log");
    }
    else
    {
        sprintf(buf, "%s < %s > %s",
                exe,
                "input.com",
                "input.log");
    }
    fprintf(stderr, "Calling '%s'\n", buf);
    if (system(buf) != 0)
    {
        gmx_fatal(FARGS, "Call to '%s' failed\n", buf);
    }
}

real QMMM_QMgaussian::call_gaussian(const t_forcerec *fr, QMMM_QMrec& qm, QMMM_MMrec& mm, rvec f[], rvec fshift[])
{
    /* normal gaussian jobs */
    static int
        step = 0;
    int
        i, j;
    real
        QMener = 0.0;
    rvec
       *QMgrad, *MMgrad;
    char
       *exe;

    snew(exe, 30);
    sprintf(exe, "%s/%s", qm.gaussian.gauss_dir, qm.gaussian.gauss_exe);
    snew(QMgrad, qm.nrQMatoms_get());
    snew(MMgrad, mm.nrMMatoms);

    write_gaussian_input(step, fr, qm, mm);
    do_gaussian(step, exe);
    QMener = read_gaussian_output(QMgrad, MMgrad, qm, mm);
    /* put the QMMM forces in the force array and to the fshift
     */
    for (i = 0; i < qm.nrQMatoms_get(); i++)
    {
        for (j = 0; j < DIM; j++)
        {
            f[i][j]      = HARTREE_BOHR2MD*QMgrad[i][j];
            fshift[i][j] = HARTREE_BOHR2MD*QMgrad[i][j];
        }
    }
    for (i = 0; i < mm.nrMMatoms; i++)
    {
        for (j = 0; j < DIM; j++)
        {
            f[i+qm.nrQMatoms_get()][j]      = HARTREE_BOHR2MD*MMgrad[i][j];
            fshift[i+qm.nrQMatoms_get()][j] = HARTREE_BOHR2MD*MMgrad[i][j];
        }
    }
    QMener = QMener*HARTREE2KJ*AVOGADRO;
    step++;
    free(exe);
    return(QMener);

} /* call_gaussian */

real QMMM_QMgaussian::call_gaussian_SH(const t_forcerec *fr, QMMM_QMrec &qm, QMMM_MMrec &mm, rvec f[], rvec fshift[])
{
    /* a gaussian call routine intended for doing diabatic surface
     * "sliding". See the manual for the theoretical background of this
     * TSH method.
     */
    static int
        step = 0;
    int
        state, i, j;
    real
        QMener = 0.0;
    static  gmx_bool
        swapped = FALSE; /* handle for identifying the current PES */
    gmx_bool
        swap = FALSE;    /* the actual swap */
    rvec
       *QMgrad, *MMgrad;
    char
       *buf;
    char
       *exe;

    snew(exe, 30);
    sprintf(exe, "%s/%s", qm.gaussian.gauss_dir, qm.gaussian.gauss_exe);
    /* hack to do ground state simulations */
    if (!step)
    {
        snew(buf, 20);
        buf = getenv("GMX_QM_GROUND_STATE");
        if (buf)
        {
            sscanf(buf, "%d", &state);
        }
        else
        {
            state = 2;
        }
        if (state == 1)
        {
            swapped = TRUE;
        }
    }
    /* end of hack */


    /* copy the QMMMrec pointer */
    snew(QMgrad, qm.nrQMatoms_get());
    snew(MMgrad, mm.nrMMatoms);
    /* at step 0 there should be no SA */
    /*  if(!step)
     * qr->bSA=FALSE;*/
    /* temporray set to step + 1, since there is a chk start */
    write_gaussian_SH_input(step, swapped, fr, qm, mm);

    do_gaussian(step, exe);
    QMener = read_gaussian_SH_output(QMgrad, MMgrad, step, qm, mm);

    /* check for a surface hop. Only possible if we were already state
     * averaging.
     */
    if (qm.surfaceHopping.SAstep > 0)
    {
        if (!swapped)
        {
            swap    = ((step != 0) && (hop(step, qm) != 0));
            swapped = swap;
        }
        else /* already on the other surface, so check if we go back */
        {
            swap    = ((step != 0) && (hop(step, qm) != 0));
            swapped = !swap; /* so swapped shoud be false again */
        }
        if (swap)            /* change surface, so do another call */
        {
            write_gaussian_SH_input(step, swapped, fr, qm, mm);
            do_gaussian(step, exe);
            QMener = read_gaussian_SH_output(QMgrad, MMgrad, step, qm, mm);
        }
    }
    /* add the QMMM forces to the gmx force array and fshift
     */
    for (i = 0; i < qm.nrQMatoms_get(); i++)
    {
        for (j = 0; j < DIM; j++)
        {
            f[i][j]      = HARTREE_BOHR2MD*QMgrad[i][j];
            fshift[i][j] = HARTREE_BOHR2MD*QMgrad[i][j];
        }
    }
    for (i = 0; i < mm.nrMMatoms; i++)
    {
        for (j = 0; j < DIM; j++)
        {
            f[i+qm.nrQMatoms_get()][j]      = HARTREE_BOHR2MD*MMgrad[i][j];
            fshift[i+qm.nrQMatoms_get()][j] = HARTREE_BOHR2MD*MMgrad[i][j];
        }
    }
    QMener = QMener*HARTREE2KJ*AVOGADRO;
    fprintf(stderr, "step %5d, SA = %5d, swap = %5d\n",
            step, static_cast<int>(qm.surfaceHopping.SAstep > 0), static_cast<int>(swapped));
    step++;
    free(exe);
    return(QMener);

} /* call_gaussian_SH */

#pragma GCC diagnostic pop
