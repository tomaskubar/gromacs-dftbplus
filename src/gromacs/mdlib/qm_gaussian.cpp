/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017 by the GROMACS development team.
 * Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
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
#include "gromacs/mdlib/qmmm.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

// When not built in a configuration with QMMM support, much of this
// code is unreachable by design. Tell clang not to warn about it.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-noreturn"

// TODO: this should be made thread-safe 

// Gaussian interface routines

void QMMM_QMgaussian::init_gaussian()
{
    char* buf = nullptr;
    int   i;

    if (!GMX_QMMM_GAUSSIAN)
    {
        gmx_fatal(FARGS,
                  "Cannot call GAUSSIAN unless linked against it. Use cmake "
                  "-DGMX_QMMM_PROGRAM=GAUSSIAN, and ensure that linking will work correctly.");
    }

    // Per layer, we make a new subdir for integral file, checkpoint files and such.
    // These dirs are stored in the QMrec for convenience.

    if (!nQMcpus) // this we do only once per layer as we call g01 externally
    {
        // Read the number of cpus and environment from the environment if set.
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

        // Gaussian settings on the system

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
            gmx_fatal(FARGS,
                      "no $GMX_QM_MODIFIED_LINKS_DIR, this is were the modified links reside.\n");
        }

        //  // reactionfield, file is needed using gaussian 
        //  if (fr->bRF)
        //  {
        //    rffile = fopen("rf.dat","w");
        //    fprintf(rffile,"%f %f\n",fr->epsilon_r,fr->rcoulomb/BOHR2NM);
        //    fclose(rffile);
        //  }
    }
    fprintf(stderr, "gaussian initialised...\n");
}


static void write_gaussian_input(int step, QMMM_QMrec& qm, QMMM_MMrec& mm)
{
    FILE* out = fopen("input.com", "w");

    // write the route
    if (qm.QMmethod_get() >= eQMmethodRHF)
    {
        fprintf(out, "%s", "%chk=input\n");
    }
    else
    {
        fprintf(out, "%s", "%chk=se\n");
    }
    if (qm.gaussian.nQMcpus > 1)
    {
        fprintf(out, "%s%3d\n", "%nprocshare=", qm.gaussian.nQMcpus);
    }
    fprintf(out, "%s%d\n", "%mem=", qm.gaussian.QMmem);
    fprintf(out, "%s%s%s", "%subst l701 ", qm.gaussian.devel_dir, "/l701\n");
    fprintf(out, "%s%s%s", "%subst l301 ", qm.gaussian.devel_dir, "/l301\n");
    fprintf(out, "%s%s%s", "%subst l9999 ", qm.gaussian.devel_dir, "/l9999\n");
    if (step)
    {
        fprintf(out, "%s", "#T ");
    }
    else
    {
        fprintf(out, "%s", "#P ");
    }
    if (qm.QMmethod_get() == eQMmethodB3LYPLAN)
    {
        fprintf(out, " %s", "B3LYP/GEN Pseudo=Read");
    }
    else
    {
        fprintf(out, " %s", eQMmethod_names[qm.QMmethod_get()]);

        if (qm.QMmethod_get() >= eQMmethodRHF)
        {
            if (qm.QMmethod_get() == eQMmethodCASSCF)
            {
                // in case of cas, how many electrons and orbitals do we need?
                fprintf(out, "(%d,%d)", qm.CASelectrons_get(), qm.CASorbitals_get());
            }
            fprintf(out, "/%s", eQMbasis_names[qm.QMbasis_get()]);
        }
    }
    if (mm.nrMMatoms)
    {
        fprintf(out, " %s", "Charge ");
    }
    if (step || qm.QMmethod_get() == eQMmethodCASSCF)
    {
        // fetch guess from checkpoint file, always for CASSCF
        fprintf(out, "%s", " guess=read");
    }
    fprintf(out, "\nNosymm units=bohr\n");

    fprintf(out, "FORCE Punch=(Derivatives) ");
    fprintf(out, "iop(3/33=1)\n\n");
    fprintf(out, "input-file generated by gromacs\n\n");
    fprintf(out, "%2d%2d\n", qm.QMcharge_get(), qm.multiplicity_get());
    for (int i = 0; i < qm.nrQMatoms_get(); i++)
    {
        fprintf(out, "%3d %10.7f  %10.7f  %10.7f\n", qm.atomicnumberQM_get(i),
                qm.xQM_get(i, XX) / BOHR2NM, qm.xQM_get(i, YY) / BOHR2NM, qm.xQM_get(i, ZZ) / BOHR2NM);
    }

    // Pseudo Potential and ECP are included here if selected (MEthod suffix LAN)
    if (qm.QMmethod_get() == eQMmethodB3LYPLAN)
    {
        fprintf(out, "\n");
        for (int i = 0; i < qm.nrQMatoms_get(); i++)
        {
            if (qm.atomicnumberQM_get(i) < 21)
            {
                fprintf(out, "%d ", i + 1);
            }
        }
        fprintf(out, "\n%s\n****\n", eQMbasis_names[qm.QMbasis_get()]);

        for (int i = 0; i < qm.nrQMatoms_get(); i++)
        {
            if (qm.atomicnumberQM_get(i) > 21)
            {
                fprintf(out, "%d ", i + 1);
            }
        }
        fprintf(out, "\n%s\n****\n\n", "lanl2dz");

        for (int i = 0; i < qm.nrQMatoms_get(); i++)
        {
            if (qm.atomicnumberQM_get(i) > 21)
            {
                fprintf(out, "%d ", i + 1);
            }
        }
        fprintf(out, "\n%s\n", "lanl2dz");
    }

    // MM point charge data
    if (mm.nrMMatoms)
    {
        fprintf(stderr, "nr mm atoms in gaussian.c = %d\n", mm.nrMMatoms);
        fprintf(out, "\n");
        for (int i = 0; i < mm.nrMMatoms; i++)
        {
            fprintf(out, "%10.7f  %10.7f  %10.7f %8.4f\n",
                    mm.xMM[i][XX]/BOHR2NM, mm.xMM[i][YY]/BOHR2NM, mm.xMM[i][ZZ]/BOHR2NM, mm.MMcharges[i]);
        }
    }
    fprintf(out, "\n");

    fclose(out);

} // write_gaussian_input

static real read_gaussian_output(rvec QMgrad[], rvec MMgrad[], QMMM_QMrec& qm, QMMM_MMrec& mm)
{
    char  buf[300];
    real  QMener;
    FILE* in = fopen("fort.7", "r");

    // (There was additional content in the file in case
    //    of QM optimizations / transition state search,
    //    which was removed.
    //
    // the next line is the energy and in the case of CAS, the energy
    // difference between the two states.
    if (nullptr == fgets(buf, 300, in))
    {
        gmx_fatal(FARGS, "Error reading Gaussian output");
    }

#if GMX_DOUBLE
    sscanf(buf, "%lf\n", &QMener);
#else
    sscanf(buf, "%f\n", &QMener);
#endif
    // next lines contain the gradients of the QM atoms
    for (int i = 0; i < qm.nrQMatoms_get(); i++)
    {
        if (nullptr == fgets(buf, 300, in))
        {
            gmx_fatal(FARGS, "Error reading Gaussian output");
        }
#if GMX_DOUBLE
        sscanf(buf, "%lf %lf %lf\n", &QMgrad[i][XX], &QMgrad[i][YY], &QMgrad[i][ZZ]);
#else
        sscanf(buf, "%f %f %f\n", &QMgrad[i][XX], &QMgrad[i][YY], &QMgrad[i][ZZ]);
#endif
    }
    // the next lines are the gradients of the MM atoms
    if (qm.QMmethod_get() >= eQMmethodRHF)
    {
        for (int i = 0; i < mm.nrMMatoms; i++)
        {
            if (nullptr == fgets(buf, 300, in))
            {
                gmx_fatal(FARGS, "Error reading Gaussian output");
            }
#if GMX_DOUBLE
            sscanf(buf, "%lf %lf %lf\n", &MMgrad[i][XX], &MMgrad[i][YY], &MMgrad[i][ZZ]);
#else
            sscanf(buf, "%f %f %f\n", &MMgrad[i][XX], &MMgrad[i][YY], &MMgrad[i][ZZ]);
#endif
        }
    }
    fclose(in);
    return (QMener);
}

static void do_gaussian(int step, char* exe)
{
    char buf[STRLEN];

    // make the call to the gaussian binary through system()
    // The location of the binary will be picked up from the
    // environment using getenv().
    if (step) // hack to prevent long inputfiles
    {
        sprintf(buf, "%s < %s > %s", exe, "input.com", "input.log");
    }
    else
    {
        sprintf(buf, "%s < %s > %s", exe, "input.com", "input.log");
    }
    fprintf(stderr, "Calling '%s'\n", buf);
    if (system(buf) != 0)
    {
        gmx_fatal(FARGS, "Call to '%s' failed\n", buf);
    }
}

real QMMM_QMgaussian::call_gaussian(QMMM_QMrec& qm, QMMM_MMrec& mm, rvec f[], rvec fshift[])
{
    static int step = 0;
    real       QMener = 0.0;
    rvec      *QMgrad, *MMgrad;
    char*      exe;

    snew(exe, 30);
    sprintf(exe, "%s/%s", qm.gaussian.gauss_dir, qm.gaussian.gauss_exe);
    snew(QMgrad, qm.nrQMatoms_get());
    snew(MMgrad, mm.nrMMatoms);

    write_gaussian_input(step, qm, mm);
    do_gaussian(step, exe);
    QMener = read_gaussian_output(QMgrad, MMgrad, qm, mm);

    // put the QMMM forces in the force array and to the fshift
    for (int i = 0; i < qm.nrQMatoms_get(); i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            f[i][j]      = HARTREE_BOHR2MD * QMgrad[i][j];
            fshift[i][j] = HARTREE_BOHR2MD * QMgrad[i][j];
        }
    }
    for (int i = 0; i < mm.nrMMatoms; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            f[i + qm.nrQMatoms_get()][j]      = HARTREE_BOHR2MD * MMgrad[i][j];
            fshift[i + qm.nrQMatoms_get()][j] = HARTREE_BOHR2MD * MMgrad[i][j];
        }
    }
    QMener = QMener * HARTREE2KJ * AVOGADRO;
    step++;
    free(exe);
    return (QMener);

} // call_gaussian

#pragma GCC diagnostic pop
