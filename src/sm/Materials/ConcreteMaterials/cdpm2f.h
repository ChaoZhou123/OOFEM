/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef CDPM2F_h
#define CDPM2F_h

#include "sm/Materials/ConcreteMaterials/concretedpm2.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "cltypes.h"
#include "sm/Materials/structuralms.h"
#include "sm/Materials/isolinearelasticmaterial.h"
#include "gausspoint.h"
#include "mathfem.h"

///@name Input fields for CDPM2F
//@{
#define _IFT_CDPM2F_Name "cdpm2f"
#define _IFT_CDPM2F_lf "lf"
//@}

namespace oofem {
/**
 * This class implements the material status associated to CDPM2F.
 * @author Peter Grassl, Dimitrios Xenos
 */
class CDPM2FStatus : public ConcreteDPM2Status
{
public:


protected:
  
public:
    /// Constructor
    CDPM2FStatus(GaussPoint *gp);

    const char *giveClassName() const override { return "CDPM2FStatus"; }

};

//   ********************************
//   *** CLASS CONCRETEDPM2   ***
//   ********************************

/**
 * This class is an extension of concretedpm2. The fibre bridging stress according to the model proposed by Lin and Li (1997) is added to the standard CDPM2 model.
 *
 * @author Peter Grassl, Chao Zhou
 */
class CDPM2F : public ConcreteDPM2
{
public:

protected:
  double lf;
  
public:
    /// Constructor
    CDPM2F(int n, Domain *d);

    void initializeFrom(InputRecord &ir) override;

    const char *giveClassName() const override { return "CDPM2F"; }
    const char *giveInputRecordName() const override { return _IFT_CDPM2F_Name; }

    //void giveRealStressVector_1d(FloatArray &answer, GaussPoint *gp, const FloatArray &totalStrain, TimeStep *tStep) override;

    FloatArrayF< 6 >giveRealStressVector_3d(const FloatArrayF< 6 > &strain, GaussPoint *gp, TimeStep *tStep) const override;

    bool hasMaterialModeCapability(MaterialMode mode) const override;

protected:
    MaterialStatus *CreateStatus(GaussPoint *gp) const override;
};
} //end namespace oofem
#endif
