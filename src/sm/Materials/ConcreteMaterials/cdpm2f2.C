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


void
CDPM2F2::initializeFrom(InputRecord &ir)
{
    // call the corresponding service for the linear elastic material
    ConcreteDPM2::initializeFrom(ir);

    //Chao: your input

    IR_GIVE_FIELD(ir, lf, _IFT_CDPM2F2_Lf);

    IR_GIVE_FIELD(ir, vf, _IFT_CDPM2F2_Vf);

    IR_GIVE_FIELD(ir, df, _IFT_CDPM2F2_Df);

    IR_GIVE_FIELD(ir, ef, _IFT_CDPM2F2_Ef);

    IR_GIVE_FIELD(ir, tau0, _IFT_CDPM2F2_Tau0);

    IR_GIVE_FIELD(ir, beta, _IFT_CDPM2F2_Beta);

    IR_GIVE_FIELD(ir, f, _IFT_CDPM2F2_f);

    IR_GIVE_FIELD(ir, sm, _IFT_CDPM2F2_Sm);

    double em;
    IR_GIVE_FIELD(ir, em, _IFT_IsotropicLinearElasticMaterial_e);

    this->eta = this->ef * this->vf / ( em * ( 1. - this->vf ) );

    this->g = 2. * ( 1. + exp(M_PI * this->f / 2.) ) / ( 4. + this->f * this->f );

    this->s0 = 1. / 2. * g * this->tau0 * this->vf * ( 1. + this->eta ) * this->lf / this->df;

    this->omega = sqrt(4. * ( 1. + this->eta ) * this->beta * this->tau0 / this->ef);

    this->k = this->omega * this->lf / ( 2. * this->df );

    this->lamda = cosh(this->k) - 1.0;

    this->deltaStar = ( 2. * this->df ) / beta * this->lamda;

    this->c = this->beta * this->lf / ( 2. * this->df );

    this->deltaCu = 0.5 * this->lf * ( this->c - 2. ) / ( 3. * this->c );

    //Todo: We need to find a way to check hardeningModulus Hp based on the fibre input and adjust it so that it is large enough to keep damage positive. Initial inclination is critical. We should then check here the value of Hp and adjust if necessary.
}



double
CDPM2F2::computeDamageParamTension(double equivStrain, double kappaOne, double kappaTwo, double le, double omegaOld, double rateFactor) const
{
    double omega = 0.;

    double fibre = 0.;
    double concrete = 0.;

    //So that damage does not turn out to be negative if function is entered for equivstrains smaller thatn e0.
    double ftTemp = this->ft * ( 1. - yieldTolDamage );

    double wfMod = this->wf;
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
    if ( equivStrain > e0 * ( 1. - yieldTolDamage ) ) {//Check if damage should start
        D = this->sm * this->deltaCu / ( le - this->sm );
        eCu = this->deltaCu / this->sm;

        eStar = this->deltaStar / this->sm;

        eUl = ( 0.5 * lf - this->deltaCu ) / le + this->deltaCu / this->sm;
        omega           = 1.;    // initial guess
        residual = 0.;
        a = 0.;
        nite = 0;
        do{
            nite++;
            if ( nite == 100 ) {
                OOFEM_ERROR("In computeDamageTension: bisection method not converged");
            }

            eCr  = kappaOne + omega * kappaTwo;
            delta = 0.;
            if ( eCr > 0 && eCr <= eCu ) {
                //               double delta = sqrt( le * eCr * D + D * D / 4. ) - D / 2.;
                delta = eCr * sm;
            } else if ( eCr > eCu && eCr <= eUl ) {
                delta = le * ( eCr - deltaCu / sm ) + deltaCu;
            }

            //Chao: I have put initialisation to start of function
            // double fibre = 0.;
            if ( eCr > 0 && eCr <= eStar ) {
                fibre = ( 2. / k * ( ( 1. - acosh(1. + lamda * delta / deltaStar) / k ) * sqrt(pow( ( 1. + lamda * delta / deltaStar ), 2.) - 1.) + ( lamda * delta ) / ( k * deltaStar ) ) * s0 );
            } else if ( eCr > eStar && eCr <= eUl ) {
                fibre = ( 1. + beta * delta / df ) * ( pow( ( 1. - 2. * delta / lf ), 2.) ) * s0;
            } else if ( eCr > eUl || eCr <= 0 ) {
                fibre = 0;
            }

            if ( this->softeningType == 0 ) {
                //Todo: implement the linear law
            } else if ( this->softeningType == 1 )      {
                //Todo: Implement the bilinear law
            } else if ( this->softeningType == 2 )      {
                concrete = ftTemp * exp(-le * ( omega * kappaTwo + kappaOne ) / wfMod);
            }

            residual = ( 1. - omega ) * this->eM * equivStrain - concrete - fibre;

            if ( residual < 0 ) {
                omega = ( omega + a ) / 2.;
            } else if ( residual > 0 )       {
                double c = omega * 2. - a;
                a = ( c + a ) / 2.;
                omega = ( c + a ) / 2.;
            }
        }while( fabs(residual / this->ft) > 0.000001 );
    } else   {
        omega = 0.;
    }


    if ( omega > 1. ) {
        omega = 1.;
    }

    if ( omega < 0. || omega < omegaOld ) {
        omega = omegaOld;
    }

    return omega;
}
} //end of namespace
