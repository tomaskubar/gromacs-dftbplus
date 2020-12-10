/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2017,2018 by the GROMACS development team.
 * Copyright (c) 2019,2020, by the GROMACS development team, led by
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

#include "qm_orca.h"

#include "config.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "gromacs/fileio/confio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/qmmm.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

// When not built in a configuration with QMMM support, much of this
// code is unreachable by design. Tell clang not to warn about it.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-noreturn"

/* ORCA interface routines */

void init_orca(QMMM_QMrec* qm)
{
    char* buf;
    snew(buf, 200);

    if (!GMX_QMMM_ORCA)
    {
        gmx_fatal(FARGS,
                  "Cannot call ORCA unless linked against it. Use cmake -DGMX_QMMM_PROGRAM=ORCA, "
                  "and ensure that linking will work correctly.");
    }

    /* ORCA settings on the system */
    buf = getenv("GMX_QM_ORCA_BASENAME");
    if (buf)
    {
        snew(qm->orca_basename, 200);
        sscanf(buf, "%s", qm->orca_basename);
    }
    else
    {
        gmx_fatal(FARGS, "$GMX_QM_ORCA_BASENAME is not set\n");
    }

    /* ORCA directory on the system */
    snew(buf, 200);
    buf = getenv("GMX_ORCA_PATH");

    if (buf)
    {
        snew(qm->orca_dir, 200);
        sscanf(buf, "%s", qm->orca_dir);
    }
    else
    {
        gmx_fatal(FARGS, "$GMX_ORCA_PATH not set, check manual\n");
    }

    fprintf(stderr, "Setting ORCA path to: %s...\n", qm->orca_dir);
    fprintf(stderr, "ORCA initialised...\n\n");
    /* since we append the output to the BASENAME.out file,
       we should delete an existent old out-file here. */
    sprintf(buf, "%s.out", qm->orca_basename);
    remove(buf);
}


static void write_orca_input(const QMMM_QMrec& qm, const QMMM_MMrec& mm)
{
    int   i;
    FILE *pcFile, *addInputFile;
    char *buf, *orcaInput, *addInputFilename, *pcFilename;

    /* write the first part of the input-file */
    snew(orcaInput, 200);
    sprintf(orcaInput, "%s.inp", qm.orca_basename);
    FILE* out = fopen(orcaInput, "w");

    snew(addInputFilename, 200);
    sprintf(addInputFilename, "%s.ORCAINFO", qm.orca_basename);
    addInputFile = fopen(addInputFilename, "r");

    fprintf(out, "#input-file generated by GROMACS\n");

    fprintf(out, "!EnGrad TightSCF\n");

    /* here we include the insertion of the additional orca-input */
    snew(buf, 200);
    if (addInputFile != nullptr)
    {
        while (!feof(addInputFile))
        {
            if (fgets(buf, 200, addInputFile) != nullptr)
            {
                fputs(buf, out);
            }
        }
    }
    else
    {
        gmx_fatal(FARGS, "No information on the calculation given in %s\n", addInputFilename);
    }

    fclose(addInputFile);

    /* write charge and multiplicity */
    fprintf(out, "*xyz %2d%2d\n", qm.QMcharge_get(), qm.multiplicity_get());

    /* write the QM coordinates */
    for (i = 0; i < qm.nrQMatoms_get(); i++)
    {
        int atomNr;
        if (qm.atomicnumberQM_get(i) == 0)
        {
            atomNr = 1;
        }
        else
        {
            atomNr = qm.atomicnumberQM_get(i);
        }
        fprintf(out, "%3d %10.7f  %10.7f  %10.7f\n", atomNr,
                qm.xQM_get(i, XX) / 0.1, qm.xQM_get(i, YY) / 0.1, qm.xQM_get(i, ZZ) / 0.1);
    }
    fprintf(out, "*\n");

    /* write the MM point charge data */
    if (mm.nrMMatoms)
    {
        /* name of the point charge file */
        snew(pcFilename, 200);
        sprintf(pcFilename, "%s.pc", qm.orca_basename);
        fprintf(out, "%s%s%s\n", "%pointcharges \"", pcFilename, "\"");
        pcFile = fopen(pcFilename, "w");
        fprintf(pcFile, "%d\n", mm.nrMMatoms);
        for (i = 0; i < mm.nrMMatoms; i++)
        {
            fprintf(pcFile, "%8.4f %10.7f  %10.7f  %10.7f\n", mm.MMcharges[i],
                    mm.xMM[i][XX] / 0.1, mm.xMM[i][YY] / 0.1, mm.xMM[i][ZZ] / 0.1);
        }
        fprintf(pcFile, "\n");
        fclose(pcFile);
    }
    fprintf(out, "\n");

    fclose(out);
} /* write_orca_input */

static real read_orca_output(rvec QMgrad[], rvec MMgrad[], const QMMM_QMrec& qm, const QMMM_MMrec& mm)
{
    char  buf[300], orca_pcgradFilename[300], orca_engradFilename[300];
    real  QMener;
    FILE *pcgrad, *engrad;

    /* the energy and gradients for the QM part are stored in the engrad file
     * and the gradients for the point charges are stored in the pc file.
     */
    sprintf(orca_engradFilename, "%s.engrad", qm.orca_basename);
    engrad = fopen(orca_engradFilename, "r");
    /* we read the energy and the gradient for the qm-atoms from the engrad file
     */
    /* we can skip the first seven lines
     */
    for (int j = 0; j < 7; j++)
    {
        if (fgets(buf, 300, engrad) == nullptr)
        {
            gmx_fatal(FARGS, "Unexpected end of ORCA output");
        }
    }
    /* now comes the energy
     */
    if (fgets(buf, 300, engrad) == nullptr)
    {
        gmx_fatal(FARGS, "Unexpected end of ORCA output");
    }
#if GMX_DOUBLE
    sscanf(buf, "%lf\n", &QMener);
#else
    sscanf(buf, "%f\n", &QMener);
#endif
    /* we can skip the next three lines
     */
    for (int j = 0; j < 3; j++)
    {
        if (fgets(buf, 300, engrad) == nullptr)
        {
            gmx_fatal(FARGS, "Unexpected end of ORCA output");
        }
    }
    /* next lines contain the gradients of the QM atoms
     * now comes the gradient, one value per line:
     * (atom1 x \n atom1 y \n atom1 z \n atom2 x ...
     */

    for (int i = 0; i < 3 * qm.nrQMatoms_get(); i++)
    {
        int k = i / 3;
        if (fgets(buf, 300, engrad) == nullptr)
        {
            gmx_fatal(FARGS, "Unexpected end of ORCA output");
        }
#if GMX_DOUBLE
        if (i % 3 == 0)
        {
            sscanf(buf, "%lf\n", &QMgrad[k][XX]);
        }
        else if (i % 3 == 1)
        {
            sscanf(buf, "%lf\n", &QMgrad[k][YY]);
        }
        else if (i % 3 == 2)
        {
            sscanf(buf, "%lf\n", &QMgrad[k][ZZ]);
        }
#else
        if (i % 3 == 0)
        {
            sscanf(buf, "%f\n", &QMgrad[k][XX]);
        }
        else if (i % 3 == 1)
        {
            sscanf(buf, "%f\n", &QMgrad[k][YY]);
        }
        else if (i % 3 == 2)
        {
            sscanf(buf, "%f\n", &QMgrad[k][ZZ]);
        }
#endif
    }
    fclose(engrad);
    /* write the MM point charge data
     */
    if (mm.nrMMatoms)
    {
        sprintf(orca_pcgradFilename, "%s.pcgrad", qm.orca_basename);
        pcgrad = fopen(orca_pcgradFilename, "r");

        /* we read the gradient for the mm-atoms from the pcgrad file
         */
        /* we can skip the first line
         */
        if (fgets(buf, 300, pcgrad) == nullptr)
        {
            gmx_fatal(FARGS, "Unexpected end of ORCA output");
        }
        for (int i = 0; i < mm.nrMMatoms; i++)
        {
            if (fgets(buf, 300, pcgrad) == nullptr)
            {
                gmx_fatal(FARGS, "Unexpected end of ORCA output");
            }
#if GMX_DOUBLE
            sscanf(buf, "%lf%lf%lf\n", &MMgrad[i][XX], &MMgrad[i][YY], &MMgrad[i][ZZ]);
#else
            sscanf(buf, "%f%f%f\n", &MMgrad[i][XX], &MMgrad[i][YY], &MMgrad[i][ZZ]);
#endif
        }
        fclose(pcgrad);
    }
    return (QMener);
}

static void do_orca(char* orca_dir, char* basename)
{

    /* make the call to the orca binary through system()
     * The location of the binary is set through the
     * environment.
     */
    char buf[100];
    sprintf(buf, "%s/%s %s.inp >> %s.out", orca_dir, "orca", basename, basename);
    fprintf(stderr, "Calling '%s'\n", buf);
    if (system(buf) != 0)
    {
        gmx_fatal(FARGS, "Call to '%s' failed\n", buf);
    }
}

real call_orca(const QMMM_QMrec& qm, const QMMM_MMrec& mm, rvec f[], rvec fshift[])
{
    /* normal orca jobs */
    static int step = 0;
    int        i, j;
    real       QMener;
    rvec *     QMgrad, *MMgrad;
    char*      exe;

    snew(exe, 30);
    sprintf(exe, "%s", "orca");
    snew(QMgrad, qm.nrQMatoms_get());
    snew(MMgrad, mm.nrMMatoms);

    write_orca_input(qm, mm);
    do_orca(qm.orca_dir, qm.orca_basename);
    QMener = read_orca_output(QMgrad, MMgrad, qm, mm);
    /* put the QMMM forces in the force array and to the fshift
     */
    for (i = 0; i < qm.nrQMatoms_get(); i++)
    {
        for (j = 0; j < DIM; j++)
        {
            f[i][j]      = HARTREE_BOHR2MD * QMgrad[i][j];
            fshift[i][j] = HARTREE_BOHR2MD * QMgrad[i][j];
        }
    }
    for (i = 0; i < mm.nrMMatoms; i++)
    {
        for (j = 0; j < DIM; j++)
        {
            f[i + qm.nrQMatoms_get()][j]      = HARTREE_BOHR2MD * MMgrad[i][j];
            fshift[i + qm.nrQMatoms_get()][j] = HARTREE_BOHR2MD * MMgrad[i][j];
        }
    }
    QMener = QMener * HARTREE2KJ * AVOGADRO;
    step++;
    free(exe);
    return (QMener);
} /* call_orca */

/* end of orca sub routines */

#pragma GCC diagnostic pop
