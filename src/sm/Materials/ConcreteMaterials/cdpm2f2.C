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

#include "cdpm2f2.h"
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
REGISTER_Material(CDPM2F2);

CDPM2F2Status::CDPM2F2Status(GaussPoint *gp) :
    ConcreteDPM2Status(gp)
{}


//   ********************************
//   *** CLASS CDPM2F2 ***
//   ********************************

#define IDM_ITERATION_LIMIT 1.e-8

CDPM2F2::CDPM2F2(int n, Domain *d) :
    ConcreteDPM2(n, d)
{}


void CDPM2F2::initializeFrom(InputRecord &ir)
{
    // call the corresponding service for the linear elastic material
    ConcreteDPM2::initializeFrom(ir);

    // Chao: your input

    IR_GIVE_FIELD(ir, lf, _IFT_CDPM2F2_Lf);
    if(lf == 0){
      OOFEM_ERROR("Lf is zero. No fibre contribution. Use CDPM2 instead.\n"); 
    }

    IR_GIVE_FIELD(ir, vf, _IFT_CDPM2F2_Vf);
    if(vf == 0){
      OOFEM_ERROR("Vf is zero. No fibre contribution. Use CDPM2 instead.\n"); 
    }

    IR_GIVE_FIELD(ir, df, _IFT_CDPM2F2_Df);
    if(df == 0){
      OOFEM_ERROR("df is zero. No fibre contribution. Use CDPM2 instead.\n"); 
    }

    IR_GIVE_FIELD(ir, ef, _IFT_CDPM2F2_Ef);
    if(ef == 0){
      OOFEM_ERROR("Ef is zero. No fibre contribution. Use CDPM2 instead.\n"); 
    }

    IR_GIVE_FIELD(ir, tau0, _IFT_CDPM2F2_Tau0);

    IR_GIVE_FIELD(ir, beta, _IFT_CDPM2F2_Beta);

    IR_GIVE_FIELD(ir, f, _IFT_CDPM2F2_f);

    IR_GIVE_FIELD(ir, sm, _IFT_CDPM2F2_Sm);

    // 0-Bisectional method; 1-Newton method;
    this->ctype = 0;    
    IR_GIVE_OPTIONAL_FIELD(ir, this->ctype, _IFT_CDPM2F2_iterType);
    if ( this->ctype > 1 ) {
        throw ValueInputException(ir, _IFT_CDPM2F2_iterType, "iteration type not implemented.");
    }

    // 0-LinLi; 1-phenomenological;
    this->bondType = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, this->bondType, _IFT_CDPM2F2_bondType);
    if ( this->bondType > 1 ) {
        throw ValueInputException(ir, _IFT_CDPM2F2_bondType, "bond type not implemented");
    }

    // 0-linear; 1-constant;
    this->crackType = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, crackType, _IFT_CDPM2F2_crackType);
    if ( this->crackType > 1 ) {
        throw ValueInputException(ir, _IFT_CDPM2F2_crackType, "crack type not implemented");
    }


    // @TODO: We need to find a way to check hardeningModulus Hp based on the fibre input and adjust it so that it is large enough to keep damage positive. Initial inclination is critical. We should then check here the value of Hp and adjust if necessary.
}


double
CDPM2F2::computeDamageParamTension(double equivStrain, double kappaOne, double kappaTwo, double le, double damageOld, double rateFactor) const
{

      //Parameters calculated here for convinience
    double eta = this->ef * this->vf / ( this->eM * ( 1. - this->vf ) );

    double g = 2. * ( 1. + exp(M_PI * this->f / 2.) ) / ( 4. + this->f * this->f );

    double s0 = 1. / 2. * g * this->tau0 * this->vf * ( 1. + eta ) * this->lf / this->df;

    double omega = sqrt(4. * ( 1. + eta ) * this->beta * this->tau0 / this->ef); //not damage

    double k = omega * this->lf / ( 2. * this->df );

    double lambda = cosh(k) - 1.;

    double deltaStar = ( 2. * this->df ) / beta * lambda;

    double c = this->beta * this->lf / ( 2. * this->df );

    double deltaCu = 0.5 * this->lf * ( c - 2. ) / ( 3. * c );


  
    double damage = 0.;
    double fibre    = 0.;
    double concrete = 0.;

    // So that damage does not turn out to be negative if function is entered for equivstrains smaller thatn e0.
    double ftTemp = this->ft * ( 1. - yieldTolDamage );

    double wfMod    = this->wf;
    double wfOneMod = this->wfOne;

    if ( this->strengthRateType > 0 ) {
        if ( this->energyRateType == 0 ) {
            wfMod /= pow(rateFactor, 2.);
            wfOneMod /= pow(rateFactor, 2.);
        } else if ( this->energyRateType == 1 ) {
            wfMod /= rateFactor;
            wfOneMod /= rateFactor;
        }
    }


    double D = 0., eCu = 0., eStar = 0., eUl = 0., residual = 0., a = 0., eCr = 0., delta = 0.;
    int nite = 0;
    if ( equivStrain > e0 * ( 1. - yieldTolDamage ) ) { // Check if damage should start
        D   = this->sm * deltaCu / ( le - this->sm );
        eCu = deltaCu / this->sm;

        //eStar = this->deltaStar / this->sm;
        eStar = deltaStar * deltaCu / ( le * deltaCu - deltaStar * ( le - this->sm ) );//linear strain

        // eStar =  ((deltaStar*deltaStar)*(le-sm)/(sm*deltaCu)+deltaStar)/le;

        eUl = ( 0.5 * lf - deltaCu ) / le + deltaCu / this->sm;
        double deltaStarS = deltaStar * 2. / lf;
        double dd1        = 2. / k * ( ( 1 - 1 / k * acosh(1 + lambda) ) * sqrt(pow( ( 1 + lambda ), 2 ) - 1.) + lambda / k ) * s0;
        double dd2        = ( 1. + beta * lf * deltaStarS / ( 2. * df ) ) * pow( ( 1 - deltaStarS ), 2 ) * s0;
        double dd         = ( dd2 - dd1 ) / s0;
        double ddd1       = ( 1. / k * ( 1. - 1. / k * acosh(1. + lambda) ) * 2. * lambda * ( 1. + lambda ) / ( deltaStarS * sqrt(pow( ( 1. + lambda ), 2 ) - 1.) ) ) * s0;
        double ddd2       = ( beta * lf / ( 2. * df ) * ( 1. - 3. * deltaStarS ) - 2. ) * ( 1. - deltaStarS ) * s0;
        double bbb2       = ddd2 * 2. / lf;
        double ddd        = ( ddd2 - ddd1 ) / s0;
        double aa         = 2. * ( ddd * deltaStarS - dd ) / pow(deltaStarS, 2);
        double b          = ddd - aa * deltaStarS;
        double fibre_star = ( 2. / k * ( ( 1. - acosh(1. + lambda * deltaStar / deltaStar) / k ) * sqrt(pow( ( 1. + lambda * deltaStar / deltaStar ), 2. ) - 1.) + ( lambda * deltaStar ) / ( k * deltaStar ) ) + ( aa * pow( ( deltaStar * 2. / lf ), 2 ) / 2. ) + ( b * deltaStar * 2. / lf ) ) * s0;
        double ddd0       =  fibre_star   / deltaStar;

        damage    = 1.; // initial guess
        residual = 0.;
        a        = 0.;
        nite     = 0;
        if ( ctype == 0 ) {
            do {
                nite++;
                if ( nite == 10000 ) {
                    OOFEM_ERROR("In computeDamageTension: bisection method not converged");
                }

                eCr   = kappaOne + damage * kappaTwo;
                delta = 0.;
                if ( eCr > 0 && eCr <= eCu ) {
                    // delta = sqrt( le * eCr * D + D * D / 4. ) - D / 2.;
                    //delta = eCr * sm;
                    delta = eCr * le / ( ( ( le - sm ) / deltaCu ) * eCr + 1. );
                } else if ( eCr > eCu && eCr <= eUl ) {
                    delta = le * ( eCr - deltaCu / sm ) + deltaCu;
                }


                if ( eCr > 0 && eCr <= eStar ) {
                    fibre = ( 2. / k * ( ( 1. - acosh(1. + lambda * delta / deltaStar) / k ) * sqrt(pow( ( 1. + lambda * delta / deltaStar ), 2. ) - 1.) + ( lambda * delta ) / ( k * deltaStar ) ) + ( aa * pow( ( delta * 2. / lf ), 2 ) / 2. ) + ( b * delta * 2. / lf ) ) * s0;
                    //fibre = fibre_star*delta/deltaStar+(bbb2-ddd0)* pow(delta,2)/deltaStar+(ddd0-bbb2)*delta;
                } else if ( eCr > eStar && eCr <= eUl ) {
                    fibre = ( 1. + beta * delta / df ) * ( pow( ( 1. - 2. * delta / lf ), 2. ) ) * s0;
                } else if ( eCr > eUl || eCr <= 0 ) {
                    fibre = 0;
                }

                if ( this->softeningType == 0 ) {
                    // Todo: implement the linear law
                } else if ( this->softeningType == 1 ) {
                    // Todo: Implement the bilinear law
                } else if ( this->softeningType == 2 ) {
                    concrete = ( 1 - vf ) * ftTemp * exp(-delta / wfMod);
                }

                residual = ( 1. - damage ) * this->eM * equivStrain - fibre - concrete;
                // printf("damage = %e, residual = %e, fibre = %e, concrete = %e\n",damage, residual, fibre, concrete);

                if ( residual < 0 ) {
                    damage = ( damage + a ) / 2.;
                } else if ( residual > 0 ) {
                    double c = damage * 2. - a;
                    a        = damage;
                    damage    = ( c + a ) / 2.;
                }
            } while ( fabs(residual / this->ft) > 0.000001 );
        } else if ( ctype == 1 ) {
            do {
                nite++;
                int newtonIte = 10000;
                if ( nite > newtonIte ) {
                    OOFEM_ERROR("algorithm not converging");
                }
                eCr                    = kappaOne + damage * kappaTwo;
                double Ddelta          = 0.;
                double dResidualDDamage = 0.;
                double dfibre          = 0.;
                double dconcrte        = 0.;

                if ( eCr > 0 && eCr <= eCu ) {
		  delta = eCr * le / ( ( le - sm ) / deltaCu * eCr + 1 );
		  Ddelta = le / pow( ( ( le - sm ) * eCr / deltaCu + 1 ), 2 );
                } else if ( eCr > eCu && eCr <= eUl ) {
		  delta  = le * ( eCr - deltaCu / sm ) + deltaCu;
		  Ddelta = le;
                }

                if ( eCr > 0 && eCr <= eStar ) {
                    fibre = ( 2. / k * ( ( 1. - acosh(1. + lambda * delta / deltaStar) / k ) * sqrt(pow( ( 1. + lambda * delta / deltaStar ), 2. ) - 1.) + ( lambda * delta ) / ( k * deltaStar ) ) + ( aa * pow( ( delta * 2. / lf ), 2 ) / 2. ) + ( b * delta * 2. / lf ) ) * s0;
                    dfibre = ( s0 * ( 1. / k * ( 1. - 1 / k * acosh(1 + lambda * delta / deltaStar) ) * 2. * lambda * ( 1 + lambda * delta / deltaStar ) / ( deltaStar * sqrt(pow( ( 1. + lambda * delta / deltaStar ), 2 ) - 1.) ) + aa * 2. * delta / lf + b ) ) * 2. / lf * Ddelta * kappaTwo;
                } else if ( eCr > eStar && eCr <= eUl ) {
                    fibre  = ( 1. + beta * delta / df ) * ( pow( ( 1. - 2. * delta / lf ), 2. ) ) * s0;
                    dfibre = ( beta * lf / ( 2. * df ) * ( 1. - 6. * delta / lf ) - 2. ) * ( 1. - 2. * delta / lf ) * s0 * 2. / lf * Ddelta * kappaTwo;
                } else if ( eCr > eUl || eCr <= 0 ) {
                    fibre  = 0.000000000000000000000000000000001;
                    dfibre = 0;
                }
                if ( this->softeningType == 0 ) {
                    // Todo: implement the linear law
                } else if ( this->softeningType == 1 ) {
                    // Todo: Implement the bilinear law
                } else if ( this->softeningType == 2 ) {
                    concrete = ( 1 - vf ) * ftTemp * exp(-delta / wfMod);
                    dconcrte = -ftTemp * Ddelta * kappaTwo / wfMod * ( 1 - vf ) * exp(-delta / wfMod);
                }
                residual        = ( 1. - damage ) * this->eM * equivStrain - fibre - concrete;
                dResidualDDamage = -this->eM * equivStrain - dfibre - dconcrte;
                damage -= residual / dResidualDDamage;
                //printf("damage=%e,residual=%e\n,",damage,residual);
            } while ( fabs(residual / this->ft) >= 1.e-8 );
            //printf("\n");
        }
    }   else {
        damage = 0.;
    }


    if ( damage > 1. ) {
        damage = 1.;
    }

    if ( damage < 0. || damage < damageOld ) {
        damage = damageOld;
    }

    printf("concrete = %e and fibre = %e\n", concrete, fibre);
    
    return damage;
}

// end of namespace
}
