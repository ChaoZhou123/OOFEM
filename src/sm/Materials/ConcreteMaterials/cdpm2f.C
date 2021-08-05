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
  
    IR_GIVE_FIELD(ir, lf, _IFT_CDPM2F_Lf);
    printf("fibre length = %e\n", lf);
    
  
    IR_GIVE_FIELD(ir, vf, _IFT_CDPM2F_Vf);
    printf("fibre volume fraction = %e\n", vf);

    IR_GIVE_FIELD(ir, df, _IFT_CDPM2F_Df);
    printf("fibre diameter = %e\n", df);

     IR_GIVE_FIELD(ir, ef, _IFT_CDPM2F_Ef);
    printf("fibre Young's modulus = %e\n", ef);

    IR_GIVE_FIELD(ir, tau0, _IFT_CDPM2F_Tau0);
    printf("initail bond strength = %e\n", tau0);

    IR_GIVE_FIELD(ir, beta, _IFT_CDPM2F_Beta);
    printf("slip hardening coefficient = %e\n", beta);

    IR_GIVE_FIELD(ir, f, _IFT_CDPM2F_f);
    printf("snubbing coefficient = %e\n", f);

    IR_GIVE_FIELD(ir, sm, _IFT_CDPM2F_Sm);
    printf("snubbing coefficient = %e\n", sm);

    double em;
    IR_GIVE_FIELD(ir, em, _IFT_IsotropicLinearElasticMaterial_e);

    this->eta = this->ef * this->vf / (em * ( 1. - this->vf ) );

    this->g = 2. * ( 1. + exp(M_PI * this->f / 2.) ) / ( 4. + this->f *this->f );

    this->s0 = 1./2. * g * this->tau0 * this->vf * (1. + this->eta)* this->lf / this->df;
      
    this->omega = sqrt(4. *(1. + this->eta) * this->beta * this->tau0 / this->ef); 

    this->k = this->omega * this->lf / (2. * this->df);  

    this->lamda = cosh(this->k)-1.0;

    this->delta_star = (2. * this->df) / beta * this->lamda;  

    this-> c = this->beta * this->lf / (2. * this->df);

    this-> delta_cu = 0.5 * this->lf * (this->c - 2.) / (3. * this->c);




  FloatArrayF< 6 >
  CDPM2F::giveRealStressVector_3d(const FloatArrayF< 6 > &fullStrainVector, GaussPoint *gp, TimeStep *tStep) const

{
  //CDPM2F status
  auto status = static_cast<CDPM2FStatus*>(this->giveStatus(gp));
  
  //Call ConcreteDPM2 to compute the stress in CDPM2
  auto stressConcrete = ConcreteDPM2::giveRealStressVector_3d(fullStrainVector,gp,tStep);

  //Peter: Add the bridging stress due to fibres
  FloatArrayF< 6 > stressFibres;

  //Peter: Calculate crackingStrain only when there is tensile damage.
  //The definition of the cracking has also a component which is present when there is no damage during hardening plasticity. I would suggest to switch on the calculation of the fibre stress only when there is tensile damage and use the part of the cracking strain in the post-peak. In CDPM2, damage is calculated by history variables which are related to plastic and elastic strains. I looked it up in the paper once more. There is kappadt1 which is norm of the rate of plastic strain, but only for the post peak (damage has started). Also, kappadt2 is related to the equivalent strain which in uniaxial stress states is something like the elastic strain in axial direction. See (44) and (45) in CDPM2 article. The only issue is that there is a ductility measure, which is used to take into account multiaxial stresses. We do not really are interested in it. For uniaxial tension, the ductility measure should be one so o that is not a problem. The suggestion is then to calculate the fibre stress as a double first. Then, come up with a method to add it to the stress somehow.


  auto tempDamageTension = status->giveTempDamageTension();
  
  if(tempDamageTension >0.){//Calculate the cracking strain and fibre stress only if there is tensile damage. This will fix problems with zero length. In CDPM2, the length is only calculated if there is damage.

    auto tempKappaDOne = status->giveTempKappaDOne();
    auto tempKappaDTwo = status->giveTempKappaDTwo();
    
    // FloatArrayF< 6 > crackingStrain;

    //  auto tempStrain = status->giveTempReducedStrain();

    //  auto tempPlasticStrain = status->giveTempPlasticStrain();
 
    double ps =  tempKappaDOne + tempDamageTension*tempKappaDTwo;

    //Peter: If possible, use variable names which tell reader what it is. 
    //  double ps = norm(crackingStrain);

    FloatArrayF< 6 > delta;
    

    double Le = status->giveLe();

    double D = sm*delta_cu/(Le-sm);

  
    double e_cu = delta_cu/sm;

    double e_star =  ((delta_star*delta_star)*(Le-sm)/(sm*delta_cu)+delta_star)/Le;

    double e_ul = (0.5*lf-delta_cu)/Le+delta_cu/sm;
  
  
 
    if(ps <= e_cu&&ps > 0){
      //  for (int i =0; i<3; i++){
    
      //Peter: I replaced with ps. Only a double.
      delta = sqrt(Le*ps*D+D*D/2.)-D/2.;
      //  }
    }
    else if(ps>e_cu&&ps<=e_ul){
      //    for (int i =0; i<3; i++){
      //Peter: I replaced with ps. Only a double.
      delta = Le*(ps-delta_cu/sm)+delta_cu;
      //  }
    
    }

    
    if (ps <= e_star&&ps > 0){
      // for (int i = 0; i<3; i++){
      stressFibres = 2./k*((1.-acosh(1.+lamda*delta/delta_star)/k)*sqrt(pow((1.+lamda*delta/delta_star),2.)-1.)+(lamda*delta)/(k*delta_star))*s0;
      //  }
    }
    else if (ps > e_star&&ps<=e_ul){
      //    for (int i = 0; i<3; i++){
      stressFibres = ((1.+beta*delta/df)*(1.-pow(2.*delta/lf,2.)))*s0;
      //    }
    }
    //  else if (ps>e_ul&&ps<=0){ //This condition does not make sense? How can it be below or equal zero but at the same time larger than e_ul? Do you mean || for or?
    else if (ps>e_ul&&ps<=0){
      //    for (int i = 0; i<3; i++){
      stressFibres = 0;
      //    }
    }
  }

  //Now add the stress fibres to to the concrete stress. However, you cannot do this to all the stresses because this would be wrong.


  
  FloatArrayF< 6 >  stress = stressConcrete+stressFibres;

  //Chao's original part

    //Peter: Add the bridging stress due to fibres
  FloatArrayF< 6 > stressFibres;
  
  FloatArrayF< 6 > delta;
  auto status = static_cast<CDPM2FStatus*>(this->giveStatus(gp));

  double Le = status->giveLe();

  double D = sm*delta_cu/(Le-sm);

  
    double e_cu = delta_cu/sm;

    double e_star =  ((delta_star*delta_star)*(Le-sm)/(sm*delta_cu)+delta_star)/Le;

    double e_ul = (0.5*lf-delta_cu)/Le+delta_cu/sm;
  
  

  FloatArrayF< 6 > delta;
  auto status = static_cast<CDPM2FStatus*>(this->giveStatus(gp));

  double Le = status->giveLe();

  double D = sm*delta_cu/(Le-sm);

  
    double e_cu = delta_cu/sm;

    double e_star =  ((delta_star*delta_star)*(Le-sm)/(sm*delta_cu)+delta_star)/Le;

    double e_ul = (0.5*lf-delta_cu)/Le+delta_cu/sm;
  

auto tempStrain = status->giveTempReducedStrain();

auto tempPlasticStrain = status->giveTempPlasticStrain();

auto tempDamageTension = status->giveTempDamageTension();

FloatArrayF< 6 > crackingStrain;

 crackingStrain = tempPlasticStrain + tempDamageTension*(tempStrain-tempPlasticStrain);
 double ps = norm(crackingStrain); 
  if(ps <= e_cu&&ps > 0){
  for (int i =0; i<3; i++){
  delta[i] = sqrt(Le*crackingStrain[i]*D+D*D/2.)-D/2.;
  }
  }
  else if(ps>e_cu&&ps<=e_ul){
    for (int i =0; i<3; i++){
      delta[i] = Le*(crackingStrain[i]-delta_cu/sm)+delta_cu;
  }
    

  if (ps <= e_star&&ps > 0){
    for (int i = 0; i<3; i++){
      stressFibres[i] = 2./k*((1.-acosh(1.+lamda*delta[i]/delta_star)/k)*sqrt(pow((1.+lamda*delta[i]/delta_star),2.)-1.)+(lamda*delta[i])/(k*delta_star))*s0;
  }
  }
    else if (ps > e_star&&ps<=e_ul){
    for (int i = 0; i<3; i++){
      stressFibres[i] = ((1.+beta*delta[i]/df)*(1.-pow(2.*delta[i]/lf,2.)))*s0;
    }
    }
    else if (ps>e_ul&&ps<=0){
      for (int i = 0; i<3; i++){
      stressFibres[i] = 0;
    }
    }
    
  FloatArrayF< 6 >  stress = stressConcrete+stressFibres;
  
  return stress;
  }

}

MaterialStatus *
CDPM2F::CreateStatus(GaussPoint *gp) const
{
    return new  CDPM2FStatus(gp);
}
} //end of namespace
