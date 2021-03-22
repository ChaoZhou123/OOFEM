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
 *               Copyright (C) 1993 - 2021   Borek Patzak
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

#include "cdpm2f.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "sm/Materials/structuralms.h"
#include "gausspoint.h"
#include "intarray.h"
#include "datastream.h"
#include "contextioerr.h"
#include "timestep.h"
#include "sm/Materials/structuralmaterial.h"
#include "sm/Materials/isolinearelasticmaterial.h"
#include "sm/CrossSections/structuralcrosssection.h"
#include "mathfem.h"
#include "classfactory.h"
#include <limits>

namespace oofem {
REGISTER_Material(CDPM2F);

CDPM2FStatus::CDPM2FStatus(GaussPoint *gp) :
    ConcreteDPM2Status(gp)
{

}

//   ********************************
//   *** CLASS CDPM2F ***
//   ********************************

CDPM2F::CDPM2F(int n, Domain *d) :
    ConcreteDPM2(n, d)
{}


bool
CDPM2F::hasMaterialModeCapability(MaterialMode mode) const
//
// returns whether receiver supports given mode
//
{
    return mode == _3dMat;
}


void
CDPM2F::initializeFrom(InputRecord &ir)
{
    // call the corresponding service for the linear elastic material
  ConcreteDPM2::initializeFrom(ir);

  /*Peter: Chao, you need to add here calls to the input macro to read in the input parameters for the fibre model.
    Below is an example. The variable should be part of the class so that you can use it anywhere in the object.
   */
  
  IR_GIVE_FIELD(ir, this->lf, _IFT_CDPM2F_lf);
  printf("fibre length = %e\n", lf);
  
}

FloatArrayF< 6 >
CDPM2F::giveRealStressVector_3d(const FloatArrayF< 6 > &fullStrainVector, GaussPoint *gp, TimeStep *tStep) const
{

  //Call ConcreteDPM2 to compute the stress in CDPM2
  auto stressConcrete = ConcreteDPM2::giveRealStressVector_3d(fullStrainVector,gp,tStep);
  
  //Peter: Add the bridging stress due to fibres
  FloatArrayF< 6 > stressFibres;
  
  FloatArrayF< 6 >  stress = stressConcrete+stressFibres;
  
  return stress;

}

MaterialStatus *
CDPM2F::CreateStatus(GaussPoint *gp) const
{
    return new  CDPM2FStatus(gp);
}
  
} //end of namespace
