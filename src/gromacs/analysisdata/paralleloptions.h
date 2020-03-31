/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2019, by the GROMACS development team, led by
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
/*! \libinternal \file
 * \brief
 * Declares gmx::AnalysisDataParallelOptions.
 *
 * \if internal
 * Implementation of this class is currently in datastorage.cpp.
 * \endif
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_analysisdata
 */
#ifndef GMX_ANALYSISDATA_PARALLELOPTIONS_H
#define GMX_ANALYSISDATA_PARALLELOPTIONS_H

namespace gmx
{

/*! \libinternal \brief
 * Parallelization options for analysis data objects.
 *
 * Methods in this class do not throw.
 *
 * \inlibraryapi
 * \ingroup module_analysisdata
 */
class AnalysisDataParallelOptions
{
public:
    //! Constructs options for serial execution.
    AnalysisDataParallelOptions();
    /*! \brief
     * Constructs options for parallel execution with given number of
     * concurrent frames.
     *
     * \param[in] parallelizationFactor
     *      Number of frames that may be constructed concurrently.
     *      Must be >= 1.
     */
    explicit AnalysisDataParallelOptions(int parallelizationFactor);

    //! Returns the number of frames that may be constructed concurrently.
    int parallelizationFactor() const { return parallelizationFactor_; }

private:
    int parallelizationFactor_;
};

} // namespace gmx

#endif
