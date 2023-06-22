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
/*! \internal \file
 * \brief
 * Defines data structure and utilities for magnetic fields
 *
 * \author Tomáš Kubař <tomas.kubar@kit.edu>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "magneticfield.h"

#include <cmath>

#include <memory>
#include <string>

#include "gromacs/commandline/filenm.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/units.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/mdtypes/imdoutputprovider.h"
#include "gromacs/mdtypes/imdpoptionprovider.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/ioptionscontainerwithsections.h"
#include "gromacs/options/optionsection.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/keyvaluetreetransform.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/strconvert.h"

namespace gmx
{

namespace
{

/*! \internal
 * \brief Describes an external magnetic field
 *
 * Class that implements a force to be evaluated in mdrun.
 * The magnetic field is constant.
 *
 * Supports operations to read and form corresponding .mdp contents.
 */
class MagneticField final : public IMDModule, public IMdpOptionProvider, public IMDOutputProvider, public IForceProvider
{
public:
    // From IMDModule
    IMdpOptionProvider* mdpOptionProvider() override { return this; }
    IMDOutputProvider*  outputProvider() override { return this; }
    void                initForceProviders(ForceProviders* forceProviders) override
    {
        if (isActive())
        {
            forceProviders->addForceProvider(this);
        }
    }

    // From IMdpOptionProvider
    void initMdpTransform(IKeyValueTreeTransformRules* transform) override;
    void initMdpOptions(IOptionsContainerWithSections* options) override;
    void buildMdpOutput(KeyValueTreeObjectBuilder* builder) const override;

    // From IMDOutputProvider
    void initOutput(FILE* fplog, int nfile, const t_filenm fnm[], bool bAppendFiles, const gmx_output_env_t* oenv) override;
    void finishOutput() override;

    // From IForceProvider
    //! \copydoc IForceProvider::calculateForces()
    void calculateForces(const ForceProviderInput& forceProviderInput,
                         ForceProviderOutput*      forceProviderOutput) override;

    void subscribeToSimulationSetupNotifications(MDModulesNotifiers* /* notifiers */) override {}
    void subscribeToPreProcessingNotifications(MDModulesNotifiers* /* notifiers */) override {}

private:
    //! Return whether or not to apply a field
    bool isActive() const;

    /*! \brief Return the field strength
     *
     * \param[in] dim The spatial direction
     * \return The field strength in Gromacs units (10^9/66485 T)
     */
    real field(int dim) const;

    //! The components of the magnetic field  (10^9 / 96485 T)
    real b0x_ = 0;
    real b0y_ = 0;
    real b0z_ = 0;
};

//! Converts dynamic parameters from new mdp format to (B0_x, B0_y, B0_z).
void convertParameters(gmx::KeyValueTreeObjectBuilder* builder, const std::string& value)
{
    const std::vector<std::string> sxt = splitString(value);
    if (sxt.empty())
    {
        return;
    }
    if (sxt.size() != 3)
    {
        GMX_THROW(InvalidInputError("Please specify 3 components (x, y, z) for the magnetic field"));
    }
    builder->addValue<real>("magnetic-field-x", fromString<real>(sxt[0]));
    builder->addValue<real>("magnetic-field-y", fromString<real>(sxt[1]));
    builder->addValue<real>("magnetic-field-z", fromString<real>(sxt[2]));
}

void MagneticField::initMdpTransform(IKeyValueTreeTransformRules* rules)
{
    rules->addRule().from<std::string>("/magnetic-field-x").toObject("/magnetic-field/x").transformWith(&convertParameters);
    rules->addRule().from<std::string>("/magnetic-field-y").toObject("/magnetic-field/y").transformWith(&convertParameters);
    rules->addRule().from<std::string>("/magnetic-field-z").toObject("/magnetic-field/z").transformWith(&convertParameters);
}

void MagneticField::initMdpOptions(IOptionsContainerWithSections* options) //, const char* sectionName)
{
    auto section = options->addSection(OptionSection("magnetic-field"));
    section.addOption(RealOption("-x").store(&b0x_));
    section.addOption(RealOption("-y").store(&b0y_));
    section.addOption(RealOption("-z").store(&b0z_));
}


void MagneticField::buildMdpOutput(KeyValueTreeObjectBuilder* builder) const
{
    std::string comment =
            R"(; Magnetic fields
; Format is: three real variables, magnetic-field-x magentic-field-y
; magnetic-field-z (unit T).)";
    builder->addValue<std::string>("comment-magnetic-field", comment);
    builder->addValue<real>("magnetic-field-x",  b0x_);
    builder->addValue<real>("magnetic-field-y",  b0y_);
    builder->addValue<real>("magnetic-field-z",  b0z_);
}


void MagneticField::initOutput(FILE* fplog, int nfile, const t_filenm fnm[], bool bAppendFiles, const gmx_output_env_t* oenv)
{
    if (isActive())
    {
        // TODO: replace with any relevant reference?
        please_cite(fplog, "Caleman2008a");
        // keep the compiler silent
        (void) nfile;
        (void) fnm;
        (void) bAppendFiles;
        (void) oenv;
    }
}

void MagneticField::finishOutput()
{
}

real MagneticField::field(int dim) const
{
    switch (dim) {
        case XX: return b0x_;
        case YY: return b0y_;
        case ZZ: return b0z_;
    }
}

bool MagneticField::isActive() const
{
    return (b0x_ != 0 || b0y_ != 0 || b0z_ != 0);
}

void MagneticField::calculateForces(const ForceProviderInput& forceProviderInput,
                                    ForceProviderOutput*      forceProviderOutput)
{
    if (isActive())
    {
        // NOTE: The non-conservative magnetic field does not have a virial
        ArrayRef<RVec> f = forceProviderOutput->forceWithVirial_.force_;

        auto chargeA = forceProviderInput.chargeA_;
        auto v       = forceProviderInput.v_;
        const rvec fieldStrength = {
            static_cast<real>(gmx::c_faraday * gmx::c_nano) * field(XX),
            static_cast<real>(gmx::c_faraday * gmx::c_nano) * field(YY),
            static_cast<real>(gmx::c_faraday * gmx::c_nano) * field(ZZ)
        };

        for (int m = 0; (m < DIM); m++)
        {
            // TODO: Check parallellism
            for (int i = 0; i < forceProviderInput.homenr_; ++i)
            {
                // NOTE: Not correct with perturbed charges
                f[i][m] += chargeA[i] * ( v[i][(m+1)%3] * fieldStrength[(m+2)%3]
                                         -v[i][(m+2)%3] * fieldStrength[(m+1)%3]);
            }
        }
    }
}

} // namespace

std::unique_ptr<IMDModule> createMagneticFieldModule()
{
    return std::make_unique<MagneticField>();
}

} // namespace gmx
