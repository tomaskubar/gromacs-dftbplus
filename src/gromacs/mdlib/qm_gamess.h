/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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
#ifndef GMX_MDLIB_QMGAMESS_H
#define GMX_MDLIB_QMGAMESS_H

#include "gromacs/mdlib/qmmm.h"

/*! \brief
 * Initialize gamess datastructures.
 *
 * \param[in] cr Commrec datastructure.
 * \param[in] qm QM forcerec.
 * \param[in] mm MM part of forcerec.
 */
//void init_gamess(const t_commrec *cr, t_QMrec *qm, t_MMrec *mm);
void init_gamess(const t_commrec* cr, const QMMM_QMrec& qm, const QMMM_MMrec& mm);

/*! \brief
 * Run calculation with Gamess.
 *
 * \param[in] qm QM part of forcerec.
 * \param[in] mm MM part of forcerec.
 * \param[in] f  Force vector.
 * \param[in] fshift Force shift vector.
 */
//real call_gamess(const t_QMrec *qm, const t_MMrec *mm,
real call_gamess(const QMMM_QMrec& qm, const QMMM_MMrec& mm, rvec f[], rvec fshift[]);


#endif
