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
{

}


//   ********************************
//   *** CLASS CONCRETE DAMAGE PLASTICITY MODEL 2 ***
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
    
    IR_GIVE_FIELD( ir, lf, _IFT_CDPM2F2_Lf );
    printf( "fibre length = %e\n", lf );

    IR_GIVE_FIELD( ir, vf, _IFT_CDPM2F2_Vf );
    printf( "fibre volume fraction = %e\n", vf );

    IR_GIVE_FIELD( ir, df, _IFT_CDPM2F2_Df );
    printf( "fibre diameter = %e\n", df );

    IR_GIVE_FIELD( ir, ef, _IFT_CDPM2F2_Ef );
    printf( "fibre Young's modulus = %e\n", ef );

    IR_GIVE_FIELD( ir, tau0, _IFT_CDPM2F2_Tau0 );
    printf( "initial bond strength = %e\n", tau0 );

    IR_GIVE_FIELD( ir, beta, _IFT_CDPM2F2_Beta );
    printf( "slip hardening coefficient = %e\n", beta );

    IR_GIVE_FIELD( ir, f, _IFT_CDPM2F2_f );
    printf( "snubbing coefficient = %e\n", f );

    IR_GIVE_FIELD( ir, sm, _IFT_CDPM2F2_Sm );
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
    printf( "delta_cu = %e\n", delta_cu );
    
}



double
CDPM2F2::computeDamageParamTension(double equivStrain, double kappaOne, double kappaTwo, double le, double omegaOld, double rateFactor) const
{

  //Chao: This is the function that you would need to edit to obtain your fibre response. The hardening part should be modelled as a stress-inelastic strain curve (not crack band approach). The  softening part should be formulated as a stress-crack opening curve.
  
    double omega = 0.;

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

    double help;
    if ( equivStrain > e0 * ( 1. - yieldTolDamage ) ) {
        if ( softeningType == 0 ) { //linear
            omega = ( this->eM * equivStrain * wfMod - ftTemp * wfMod + ftTemp * kappaOne * le ) /
                    ( this->eM * equivStrain * wfMod - ftTemp * le * kappaTwo );
        } else if ( softeningType == 1 ) { //bilinear: Calculate damage parameter for both parts of bilinear curve  and check which fulfils limits.
            omega = ( this->eM * equivStrain * wfOneMod - ftTemp * wfOneMod - ( this->ftOne - ftTemp ) * kappaOne * le ) /
                    ( this->eM * equivStrain * wfOneMod + ( this->ftOne - ftTemp ) * le * kappaTwo );
            help = le * kappaOne + le * omega * kappaTwo;

            if ( help >= 0. && help < wfOneMod ) {
                return omega;
            }

            omega = ( this->eM * equivStrain * ( wfMod - wfOneMod ) - this->ftOne * ( wfMod - wfOneMod ) +
                      this->ftOne * kappaOne * le  - this->ftOne * wfOneMod ) /
                    ( this->eM * equivStrain * ( wfMod - wfOneMod )  - this->ftOne * le * kappaTwo );
            help = le * kappaOne + le * omega * kappaTwo;

            if ( help > wfOneMod && help < wfMod ) {
                return omega;
            }
        } else if ( softeningType == 2 ) { //exponential: Iterative solution
            omega = 1.; //Initial guess
            double residual = 0.;
            double dResidualDOmega = 0.;
            int nite = 0;

            do {
                nite++;

                residual  = ( 1 - omega ) * this->eM * equivStrain - ftTemp * exp(-le * ( omega * kappaTwo + kappaOne ) / wfMod);
                dResidualDOmega = -this->eM * equivStrain + ftTemp * le * kappaTwo / wfMod * exp(-le * ( omega * kappaTwo + kappaOne ) / wfMod);

                omega -= residual / dResidualDOmega;
                if ( nite > newtonIter ) {
                    OOFEM_ERROR("algorithm not converging");
                }
            } while ( fabs(residual / ftTemp) >= 1.e-8 );
        }
    } else {
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
