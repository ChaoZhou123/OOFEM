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

FloatMatrixF< 6, 6 >
CDPM2F::give3dMaterialStiffnessMatrix(MatResponseMode mode,
                                            GaussPoint *gp,
                                            TimeStep *tStep) const
{
  return this->linearElasticMaterial.give3dMaterialStiffnessMatrix(mode, gp, tStep);
}

void
CDPM2F::initializeFrom(InputRecord &ir)
{
    // call the corresponding service for the linear elastic material
    ConcreteDPM2::initializeFrom( ir );

    IR_GIVE_FIELD( ir, lf, _IFT_CDPM2F_Lf );
    printf( "fibre length = %e\n", lf );

    IR_GIVE_FIELD( ir, vf, _IFT_CDPM2F_Vf );
    printf( "fibre volume fraction = %e\n", vf );

    IR_GIVE_FIELD( ir, df, _IFT_CDPM2F_Df );
    printf( "fibre diameter = %e\n", df );

    IR_GIVE_FIELD( ir, ef, _IFT_CDPM2F_Ef );
    printf( "fibre Young's modulus = %e\n", ef );

    IR_GIVE_FIELD( ir, tau0, _IFT_CDPM2F_Tau0 );
    printf( "initail bond strength = %e\n", tau0 );

    IR_GIVE_FIELD( ir, beta, _IFT_CDPM2F_Beta );
    printf( "slip hardening coefficient = %e\n", beta );

    IR_GIVE_FIELD( ir, f, _IFT_CDPM2F_f );
    printf( "snubbing coefficient = %e\n", f );

    IR_GIVE_FIELD( ir, sm, _IFT_CDPM2F_Sm );
    printf( "minimum crack spacing = %e\n", sm );

    double em;
    IR_GIVE_FIELD( ir, em, _IFT_IsotropicLinearElasticMaterial_e );

    this->eta = this->ef * this->vf / ( em * ( 1. - this->vf ) );

    this->g = 2. * ( 1. + exp( M_PI * this->f / 2. ) ) / ( 4. + this->f * this->f );

    this->s0 = 1. / 2. * g * this->tau0 * this->vf * ( 1. + this->eta ) * this->lf / this->df;

    this->omega = sqrt( 4. * ( 1. + this->eta ) * this->beta * this->tau0 / this->ef );

    this->k = this->omega * this->lf / ( 2. * this->df );

    this->lamda = cosh( this->k ) - 1.0;

    this->delta_star = ( 2. * this->df ) / beta * this->lamda;

    this->c = this->beta * this->lf / ( 2. * this->df );

    this->delta_cu = 0.5 * this->lf * ( this->c - 2. ) / ( 3. * this->c );
}

FloatArrayF< 6 >
  CDPM2F::giveRealStressVector_3d(const FloatArrayF< 6 > &fullStrainVector, GaussPoint *gp, TimeStep *tStep) const

{
  //CDPM2F status
  auto status = static_cast<CDPM2FStatus*>(this->giveStatus(gp));
  
  //Call ConcreteDPM2 to compute the stress in CDPM2
  auto stressConcrete = ConcreteDPM2::giveRealStressVector_3d(fullStrainVector,gp,tStep);

  auto tempDamageTension = status->giveTempDamageTension();

  FloatArrayF< 6 > stressFibres;

  auto tempStateFlag = status->giveTempStateFlag();

  if (tempStateFlag == CDPM2FStatus::ConcreteDPM2_Elastic || tempStateFlag == CDPM2FStatus::ConcreteDPM2_Plastic){
    printf("No crack yet. Fibres should not be active\n");
  }
  else if(tempDamageTension >0. && (tempStateFlag == CDPM2FStatus::ConcreteDPM2_Damage || tempStateFlag == CDPM2FStatus::ConcreteDPM2_PlasticDamage)){//Calculate the cracking strain and fibre stress only if there is tensile damage. This will fix problems with zero length. In CDPM2, the length is only calculated if there is damage.

    //auto tempKappaDOne = status->giveTempKappaDOne();
    //auto tempKappaDTwo = status->giveTempKappaDTwo();

    auto tempStrain = status->giveTempReducedStrain();

    auto tempPlasticStrain = status->giveTempPlasticStrain();

    FloatArrayF< 6 > crackingStrain;
    crackingStrain = tempPlasticStrain + tempDamageTension*(tempStrain-tempPlasticStrain);

    //Calculate principal values of cracking strain
    //Determine the principal values of the strain
    auto tmp = computePrincipalValDir( from_voigt_strain(crackingStrain) );  ///@todo CHECK
    auto principalCrackingStrain = tmp.first;
    auto strainPrincipalDir = tmp.second;

    //Peter: Why do you calculate the norm of the cracking strain. This would give you a positive value even if the cracking strain is negative. Better to do it component by component?  
    //  double ps = norm(crackingStrain);

    FloatArrayF< 3 > delta;
    
    double Le = status->giveLe();

    double D = sm*delta_cu/(Le-sm);

  
    double e_cu = delta_cu/sm;

    double e_star =  ((delta_star*delta_star)*(Le-sm)/(sm*delta_cu)+delta_star)/Le;

    double e_ul = (0.5*lf-delta_cu)/Le+delta_cu/sm;
  
  
    for (int i =1; i<=3; i++){ 
      if( principalCrackingStrain.at(i) > 0 && principalCrackingStrain.at(i) <= e_cu ){
	delta.at(i) = sqrt(Le*principalCrackingStrain.at(i)*D+D*D/2.)-D/2.;
      }
      else if( principalCrackingStrain.at(i) > e_cu && principalCrackingStrain.at(i) <= e_ul ){
	delta.at(i) = Le*(principalCrackingStrain.at(i)-delta_cu/sm)+delta_cu;
      }
    }

  FloatArrayF<6> principalStressFibres; //Only use first 3 components, but need to have 6 to be able to convert it later
    
    for (int i = 1; i<=3; i++){
        if ( principalCrackingStrain.at( i ) > 0 && principalCrackingStrain.at( i ) <= e_star )
	  principalStressFibres.at( i ) =( 2. / k * ( ( 1. - acosh( 1. + lamda * delta.at( i ) / delta_star ) / k ) * sqrt( pow( ( 1. + lamda * delta.at( i ) / delta_star ), 2. ) - 1. ) + ( lamda * delta.at( i ) ) / ( k * delta_star ) ) * s0);
        else if ( principalCrackingStrain.at( i ) > e_star && principalCrackingStrain.at( i ) <= e_ul )
	  principalStressFibres.at( i ) =  ( 1. + beta * delta.at( i ) / df ) * ( pow( (1.-2. * delta.at( i ) / lf), 2. )) * s0;
        else if ( principalCrackingStrain.at( i ) > e_ul || principalCrackingStrain.at( i ) <= 0 )
            principalStressFibres.at( i ) = 0;
    }

    stressFibres = transformStressVectorTo(strainPrincipalDir, principalStressFibres, 1);    

    
  }//end of loading
  else{
    printf("Here we need to implement the unloading. For this we need to have a measure of irreversible and reversible displacements.\n");
  }
  
  FloatArrayF< 6 >  stress = stressConcrete+stressFibres;
  status->letTempStrainVectorBe(fullStrainVector);
  status->letTempStressVectorBe(stress);
  return stress;
}

MaterialStatus *
CDPM2F::CreateStatus(GaussPoint *gp) const
{
    return new  CDPM2FStatus(gp);
}
} //end of namespace
