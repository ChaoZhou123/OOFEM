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
#define _IFT_CDPM2F2_Vf0 "vf0"
#define _IFT_CDPM2F2_Df "df"
#define _IFT_CDPM2F2_Tau0 "tau0"
#define _IFT_CDPM2F2_Beta "beta"
#define _IFT_CDPM2F2_f "f"
#define _IFT_CDPM2F2_Ef "ef"
#define _IFT_CDPM2F2_Sm "sm"
#define _IFT_CDPM2F2_Alpha "alpha"
#define _IFT_CDPM2F2_convergenceType "ctype"
#define _IFT_CDPM2F2_deltarelation "drelation"
#define _IFT_CDPM2F2_fibredebondingtype "fdtype"
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
 * Info about the model
 *
 * @author Chao Zhou, Peter Grassl
 */
class CDPM2F2 : public ConcreteDPM2
{
public:

protected:

    //Chao: These are the varialbes that you had in cdpm2f

    double lf = 0.;
    double vf0 = 0.;
    double vf = 0.;
    double vfm = 0.;
    double z = 0.;
    double df = 0.;
    double ef = 0.;
    double tau0 = 0.;
    double beta = 0.;
    double f = 0.;
    double eta = 0.;
    double g = 0.;
    double s0 = 0.;
    double omega = 0.;
    double k = 0.;
    double lamda = 0.;
    double deltaStar = 0.;
    double c = 0.;
    double deltaCu = 0.;
    double sm = 0.;
    double alpha = 0.;
    double alphamin = 0.;
    int ctype = 0;
    int drelation = 1;
    int fdtype = 0;
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
