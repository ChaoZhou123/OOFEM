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

#ifndef CDPM2F2_h
#define CDPM2F2_h

#include "sm/Materials/ConcreteMaterials/concretedpm2.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "cltypes.h"
#include "sm/Materials/structuralms.h"
#include "sm/Materials/isolinearelasticmaterial.h"
#include "gausspoint.h"
#include "mathfem.h"

//Chao: Here we need to add the input parameters which are in your cdpm2f

#define _IFT_CDPM2F2_Name "cdpm2f2"
#define _IFT_CDPM2F2_Lf "lf"
#define _IFT_CDPM2F2_Vf "vf"
#define _IFT_CDPM2F2_Df "df"
#define _IFT_CDPM2F2_Tau0 "tau0"
#define _IFT_CDPM2F2_Beta "beta"
#define _IFT_CDPM2F2_f "f"
#define _IFT_CDPM2F2_Ef "ef"
#define _IFT_CDPM2F2_Sm "sm"
#define _IFT_CDPM2F2_iterType "ctype"
#define _IFT_CDPM2F2_crackType "cracktype"
#define _IFT_CDPM2F2_bondType "bondtype"

//@}

namespace oofem {
/**
 * This class implements the material status associated to CDPM2F2.
 * @author Peter Grassl, Dimitrios Xenos
 */
class CDPM2F2Status : public ConcreteDPM2Status
{
public:
    /// Constructor
    CDPM2F2Status(GaussPoint *gp);
};


//   ********************************
//   *** CLASS CDPM2F2   ***
//   ********************************

/**
 * CDPM2F is an extension of CDPM2 to strain hardening. It uses the fibre model in "Lin and Li (1997): Crack bridging in fibre reinforced cementitious composites with slip-hardening interfaces, Journal of the Mechanics and Physics of Solids, Vol. 45, No. 5, pp. 763-787. The stress-displacement relation is transformed into a stress-strain response by phenomenologically introducing a law for the evolution of the crack spacing. The resulting stress-cracking strain relation is then incorporated into CDPM2 by adjusting the tensile damage evolution.
 *
 * @author Chao Zhou, Peter Grassl
 */
class CDPM2F2 : public ConcreteDPM2
{
public:

protected:

  //Input parameters
  ///length of fibre
    double lf = 0.;
  ///volume of fibre
    double vf = 0.;
  ///diameter of fibre
    double df = 0.;
  ///Young's modulus of fibre
    double ef = 0.;
  ///shear strength of fibre matrix interface
    double tau0 = 0.;
  ///hardening parameter of fibre matrix interface
    double beta = 0.;
  ///snabbing factor
    double f = 0.;
  ///crack spacing
   double sm = 0.;
  /// type of iterative approach to calulate damage varialble
  /// 0 = bisection method, 1-Newton method 
  int ctype = 0;
  /// type of debonding relationship
  /// 0 = debonding according to Lin and Li (with small adjustments to achieve smooth transition)
  /// 1 = phenomenological debonding to avoid issues with hardening
  int bondType = 0;

  ///type of cracking strain assumption
  /// 0 = linear law
  /// 1 = constant
  int crackType = 0;
  

public:
    /// Constructor
    CDPM2F2(int n, Domain *d);

    void initializeFrom(InputRecord &ir) override;

    const char *giveClassName() const override { return "CDPM2F2"; }
    const char *giveInputRecordName() const override { return _IFT_CDPM2F2_Name; }

    /// Compute damage parameter in tension.
    virtual double computeDamageParamTension(double equivStrain, double kappaOne, double kappaTwo, double le, double omegaOld, double rateFactor) const;
};
} //end namespace oofem
#endif
