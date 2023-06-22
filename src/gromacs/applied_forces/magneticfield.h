/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2015- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \libinternal \file
 * \brief
 * Declares module creation function for applied magneitc fields
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Tomáš Kubař <tomas.kubar@kit.edu>
 * \ingroup module_applied_forces
 * \inlibraryapi
 */
#ifndef GMX_APPLIED_FORCES_MAGNETICFIELD_H
#define GMX_APPLIED_FORCES_MAGNETICFIELD_H

#include <memory>

namespace gmx
{

class IMDModule;

/*! \brief
 * Creates a module for an external magnetic field.
 *
 * The returned class describes the magneitc field that can be applied
 * to all charges in a simulation. The field is a vector with 3 components:
 *     B = (Bx, By, Bz)
 *
 * force is charge * (velocity x B)
 * its unit is kJ mol^-1 nm^-1 = e * (nm ps^-1 * (kJ mol^-1 nm^-2 ps^-1 / e))
 * internal Gromacs unit of magnetic field is therefore
 * 1 kJ mol^-1 nm^-2 ps^-1 / e = 10^9 / (Faraday const.) T = 10364.27 T
 *
 * WARNING -- FOR E FIELD BUT MAY BE RELEVANT HERE AS WELL:
 * There can be problems with the virial.
 * Since the field is not self-consistent this is unavoidable.
 * For neutral molecules the virial is correct within this approximation.
 * For neutral systems with many charged molecules the error is small.
 * But for systems with a net charge or a few charged molecules
 * the error can be significant when the field is high.
 * Solution: implement a self-consistent electric field into PME.
 */
std::unique_ptr<IMDModule> createMagneticFieldModule();

} // namespace gmx

#endif
