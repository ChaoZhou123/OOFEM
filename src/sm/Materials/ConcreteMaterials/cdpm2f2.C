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

        IR_GIVE_FIELD(ir, vf0, _IFT_CDPM2F2_Vf0);

        IR_GIVE_FIELD(ir, df, _IFT_CDPM2F2_Df);

        IR_GIVE_FIELD(ir, ef, _IFT_CDPM2F2_Ef);

        IR_GIVE_FIELD(ir, tau0, _IFT_CDPM2F2_Tau0);

        IR_GIVE_FIELD(ir, beta, _IFT_CDPM2F2_Beta);

        IR_GIVE_FIELD(ir, f, _IFT_CDPM2F2_f);

        IR_GIVE_FIELD(ir, sm, _IFT_CDPM2F2_Sm);

        IR_GIVE_FIELD(ir, alpha, _IFT_CDPM2F2_Alpha);

        // 0-Bisectional method; 1-Newton method;
        IR_GIVE_FIELD(ir, ctype, _IFT_CDPM2F2_convergenceType);

        if ( ctype > 1 ) {
            throw ValueInputException(ir, _IFT_CDPM2F2_convergenceType, "convergence type not implemented");
        }

        //0-constant relation; 1-linear relation; 2-sigmoid relation;
        IR_GIVE_FIELD(ir, drelation, _IFT_CDPM2F2_deltarelation);

        if ( drelation > 2 ) {
            throw ValueInputException(ir, _IFT_CDPM2F2_deltarelation, "delta relation not implemented");
        }


        IR_GIVE_FIELD(ir, fdtype, _IFT_CDPM2F2_fibredebondingtype);

        if ( fdtype > 2 ) {
            throw ValueInputException(ir, _IFT_CDPM2F2_fibredebondingtype, "fibre debonding type not implemented");
        }
        double em;
        IR_GIVE_FIELD(ir, em, _IFT_IsotropicLinearElasticMaterial_e);

        this->vf = ( 1. + log(alpha) ) * this->vf0;
        //this->vf = vf0;



        this->eta = this->ef * this->vf / ( em * ( 1. - this->vf ) );

        this->g = 2. * ( 1. + exp(M_PI * this->f / 2.) ) / ( 4. + this->f * this->f );

        this->s0 = 1. / 2. * g * this->tau0 * this->vf * ( 1. + this->eta ) * this->lf / this->df;

        this->omega = sqrt(4. * ( 1. + this->eta ) * this->beta * this->tau0 / this->ef);

        this->k = this->omega * this->lf / ( 2. * this->df );

        this->lamda = cosh(this->k) - 1.0;

        this->deltaStar = ( 2. * this->df ) / beta * this->lamda;

        this->c = this->beta * this->lf / ( 2. * this->df );

        this->deltaCu = ( 0.5 * this->lf * ( this->c - 2. ) / ( 3. * this->c ) );

        double ftTemp = this->ft * ( 1. - yieldTolDamage );

	this->z = ftTemp*27.*beta*c / (4.*g*tau0*pow((c+1.),3.)); 
	
	this->vfm= ( -( this->em + this->em * this->z ) + sqrt(pow(( this-> em + this->em * z ),2.) + 4.* ( this->ef - this->em ) * this->em * this->z ) ) / ( 2. * ( this->ef - this->em ) );  
	 
        this->alphamin = exp( ( vfm - vf0 ) / vf0);

        this->ap = ( 1. + beta * deltaStar / df ) * pow( ( 1. - deltaStar * 2. / lf ), 2. ) - 2. * lamda / ( pow(k, 2.) );

        if ( alpha < alphamin ) {
            OOFEM_ERROR("alpha should be larger than alphamin %e\n", alphamin);
            printf("alphamin= %e\n", alphamin);
        }
	
    }


    double
    CDPM2F2::computeDamageParamTension(double equivStrain, double kappaOne, double kappaTwo, double le, double omegaOld, double rateFactor) const
    {
        double omega = 0.;

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


        double eCu = 0., eCu2 = 0., eStar = 0., eUl = 0., residual = 0., a = 0., eCr = 0.00, delta = 0., e0 = 0., gammac = 0., p = 0., q = 0., A = 0., B = 0., C = 0, D = 0, x = 0., y = 0., z = 0., o = 0., he_min = 0., he_max = 0., Z = 0., Y = 0.;
        e0 = kappaOne + ftTemp / eM;
        gammac = ( 1 - alpha ) * e0 * sm / ( deltaCu * ( 1 - alphamin ) ) + ( ( alpha - alphamin ) / ( 1 - alphamin ) );
        eCu = this->deltaCu * gammac / this->sm;
        eCu2 = this->deltaCu  / ( this->sm );
        p = 1. - ( 2. * deltaCu / lf );
        q = 1. + ( beta * deltaCu / df );

        A          =   ( eCu2 ) / ( q * pow(p, 2) ) * ( gammac * le - sm ) * 4. * beta / ( pow(lf, 2) * df * le );
        B          =   ( eCu2 ) / ( q * pow(p, 2) ) * ( gammac * le - sm ) * 4. * ( df - beta * lf ) / ( df * pow(lf, 2) * le );
        C          =   ( eCu2 ) / ( q * pow(p, 2) ) * ( gammac * le - sm ) * ( beta * lf - 4. * df ) / ( df * lf * le ) + 1. / le;
        D          =   ( eCu2 ) / ( q * pow(p, 2) ) * ( gammac * le - sm ) / le;
        // p=(gammac*le-sm)*deltaCu/sm * 1./((1.+beta*deltaCu/df)* pow((1.-2*deltaCu/lf),2));
        //A=4.*beta*p/(lf*df*le);
        //B=p*(4.*df-4.*beta*lf)/(df* pow(lf,2)*le);
        //C= 1/le+(lf*beta-4*df)*p/(df*lf*le);
        //D = p/le;
        he_min = ( sm + 3. * q * pow(p, 2) * df * pow(lf, 2) * beta / ( eCu2 * ( pow( ( beta * lf + 2. * df ), 2 ) ) ) ) / gammac;
        if ( le > he_min ) {
            OOFEM_ERROR("element size should  not be larger than he_max %e\n", he_min);
            printf("he_max= %e\n", he_min);
        }
        he_max = lf * sm / ( 2 * gammac * deltaCu );
        if ( le > he_max ) {
            OOFEM_ERROR("element size should not be larger than he_max %e\n", he_max);
            printf("he_max= %e\n", he_max);
        }
        //double gammad = gammac;

        double gammad = gammac + ( 1 - gammac ) * ( eCu - eCr ) / eCu;
        int nite = 0;
        if ( equivStrain > e0 * ( 1. - yieldTolDamage ) ) { // Check if damage should start
            if ( drelation == 0 ) {
                eStar = this->deltaStar * gammad / this->sm;
            }

            if ( drelation == 1 ) {
                eStar = ( this->deltaStar * deltaCu  / ( le * deltaCu - deltaStar * gammad * ( le - this->sm ) ) ) * gammad;
            }
            //eStar = this->deltaStar / this->sm;
            //eStar = this->deltaStar*deltaCu / (le*deltaCu-deltaStar*(le-this->sm));//linear strain
            // eStar =  ((deltaStar*deltaStar)*(le-sm)/(sm*deltaCu)+deltaStar)/le;
            if ( drelation == 2 ) {
                eStar = ( -0.1 * eCu * log(1 - deltaStar / ( deltaCu  ) * ( 1 - exp(-10.) ) ) ) * gammad;
            }

            eUl = 0.5;
            double deltaStarS = deltaStar * 2. / lf;
            double dd1        = 2. / k * ( ( 1 - 1 / k * acosh(1 + lamda) ) * sqrt(pow( ( 1 + lamda ), 2 ) - 1.) + lamda / k ) * s0;
            double dd2        = ( 1. + beta * lf * deltaStarS / ( 2. * df ) ) * pow( ( 1 - deltaStarS ), 2 ) * s0;
            double dd         = ( dd2 - dd1 ) / s0;
            double ddd1       = ( 1. / k * ( 1. - 1. / k * acosh(1. + lamda) ) * 2. * lamda * ( 1. + lamda ) / ( deltaStarS * sqrt(pow( ( 1. + lamda ), 2 ) - 1.) ) ) * s0;
            double ddd2       = ( beta * lf / ( 2. * df ) * ( 1. - 3. * deltaStarS ) - 2. ) * ( 1. - deltaStarS ) * s0;
            double bbb2       = ddd2 * 2. / lf;
            double ddd        = ( ddd2 - ddd1 ) / s0;
            double aa         = 2. * ( ddd * deltaStarS - dd ) / pow(deltaStarS, 2);
            double b          = ddd - aa * deltaStarS;
            double fibre_star = ( 2. / k * ( ( 1. - acosh(1. + lamda * deltaStar / deltaStar) / k ) * sqrt(pow( ( 1. + lamda * deltaStar / deltaStar ), 2. ) - 1.) + ( lamda * deltaStar ) / ( k * deltaStar ) ) + ( aa * pow( ( deltaStar * 2. / lf ), 2 ) / 2. ) + ( b * deltaStar * 2. / lf ) ) * s0;
            double ddd0       =  fibre_star   / deltaStar;





            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            omega    = 1.;// initial guess
            residual = 0.;
            a        = 0.;
            nite     = 0;
            if ( ctype == 0 ) {
                do {
                    nite++;
                    if ( nite == 100000 ) {
                        OOFEM_ERROR("In computeDamageTension: bisection method not converged");
                    }

                    eCr   = kappaOne + omega * kappaTwo;
                    delta = 0.;
                    if ( eCr >= 0 && eCr <= eCu ) {
                        // delta = sqrt( le * eCr * D + D * D / 4. ) - D / 2.;
                        if ( drelation == 0 ) {
                            delta = ( eCr * sm ) / gammad;//constant delta relation
                        }
                        if ( drelation == 1 ) {
                            delta = ( eCr / gammad * le / ( ( le - sm ) / deltaCu   * eCr + 1. ) );//linear delta relation
                        }
                        if ( drelation == 2 ) {
                            delta = ( deltaCu * ( 1 - exp(-eCr / ( 0.1 * eCu ) ) ) / ( 1 - exp(-( eCu ) / ( 0.1 * eCu ) ) ) );//sigmoid delta relation
                        }
                    } else if ( eCr > eCu && eCr <= eUl ) {
                        //OOFEM_ERROR("In computeDamageTension: bisection method not converged");
                        //delta = le * ( eCr - eCu ) + deltaCu;
                        x = -B / ( 3. * A );
                        o = pow( ( pow( ( ( 27. * pow(A, 2) * ( D - eCr ) - 9. * A * B * C + 2 * pow(B, 3) ) / ( 54. * pow(A, 3) ) ), 2 ) + pow( ( ( 3 * A * C - pow(B, 2) ) / ( 9. * pow(A, 2) ) ), 3 ) ), 1. / 2 );
                        Y = ( 0. - ( 27. * pow(A, 2) * ( D - eCr ) - 9. * A * B * C + 2 * pow(B, 3) ) / ( 54. * pow(A, 3) ) ) + o;
                        Z = ( ( 0. - ( 27. * pow(A, 2) * ( D - eCr ) - 9. * A * B * C + 2 * pow(B, 3) ) / ( 54. * pow(A, 3) ) ) - o );
                        if ( Y <= 0 && Z <= 0 ) {
                            y = -pow(-Y, 1. / 3);
                            z = -pow(-Z, 1. / 3);
                        } else if ( Y > 0 && Z <= 0 )    {
                            y = pow(Y, 1. / 3);
                            z = -pow(-Z, 1. / 3);
                        } else if ( Y <= 0 && Z > 0 )    {
                            y = -pow(-Y, 1. / 3);
                            z = pow(Z, 1. / 3);
                        } else if ( Y > 0 && Z > 0 )    {
                            y = pow(Y, 1. / 3);
                            z = pow(Z, 1. / 3);
                        }

                        delta = x  + y + z;
                    }



                    // Chao: I have put initialisation to start of function
                    //  double fibre = 0.;

                    if ( eCr >= 0 && eCr <= eStar ) {
                        if ( fdtype == 0 ) {
                            fibre = ( 2. / k * ( ( 1. - acosh(1. + lamda * delta / deltaStar) / k ) * sqrt(pow( ( 1. + lamda * delta / deltaStar ), 2. ) - 1.) + ( lamda * delta ) / ( k * deltaStar ) ) + ( aa * pow( ( delta * 2. / lf ), 2 ) / 2. ) + ( b * delta * 2. / lf ) ) * s0;
                        }
                        if ( fdtype == 1 ) {
                            fibre = fibre_star * delta / deltaStar + ( bbb2 - ddd0 ) * pow(delta, 2) / deltaStar + ( ddd0 - bbb2 ) * delta;
                        }
                        if ( fdtype == 2 ) {
                            fibre = ( 2. / k * ( ( 1. - acosh(1. + lamda * delta / deltaStar) / k ) * sqrt(pow( ( 1. + lamda * delta / deltaStar ), 2. ) - 1.) + ( lamda * delta ) / ( k * deltaStar ) ) + ap * delta / deltaStar ) * s0;
                        }
                    } else if ( eCr > eStar && eCr <= eUl ) {
                        fibre = ( 1. + beta * delta / df ) * ( pow( ( 1. - 2. * delta / lf ), 2. ) ) * s0;
                        if ( fibre <= 0. ) {
                            fibre = 0.;
                            //printf("omega = %e, residual = %e, fibre = %e, concrete = %e\n",omega, residual, fibre, concrete);
                        }
                    } else if ( eCr > eUl || eCr < 0 ) {
                        fibre = 0.;
                    }

                    if ( this->softeningType == 0 ) {
                        // Todo: implement the linear law
                    } else if ( this->softeningType == 1 ) {
                        // Todo: Implement the bilinear law
                    } else if ( this->softeningType == 2 ) {
                        concrete = ( 1 - vf ) * ftTemp * exp(-delta / wfMod);
                    }

                    residual = ( 1. - omega ) * this->eM * equivStrain - fibre - concrete;
                    //printf("e_cu = %e\n",eCu);
                    if ( residual < 0 ) {
                        omega = ( omega + a ) / 2.;
                    } else if ( residual > 0 ) {
                        double c = omega * 2. - a;
                        a        = omega;
                        omega    = ( c + a ) / 2.;
                    }
                } while ( fabs(residual / this->ft) > 0.000001 );
            } else if ( ctype == 1 ) {
                do {
                    nite++;
                    int newtonIte = 10000;
                    if ( nite > newtonIte ) {
                        OOFEM_ERROR("algorithm not converging");
                    }
                    eCr                    = kappaOne + omega * kappaTwo;
                    double Ddelta          = 0.;
                    double dResidualDOmega = 0.;
                    double dfibre          = 0.;
                    double dconcrte        = 0.;

                    if ( eCr > 0 && eCr < eCu ) {
                        // delta = sqrt( le * eCr * D + D * D / 4. ) - D / 2.;
                        // delta = eCr * sm;
                        if ( drelation == 0 ) {
                            delta = eCr * sm;//constant delta relation
                            Ddelta = sm;
                        }
                        if ( drelation == 1 ) {
                            delta = eCr * le / ( ( ( le - sm ) / deltaCu ) * eCr + 1. );//linear delta relation
                            Ddelta = le / pow( ( ( le - sm ) * eCr / deltaCu + 1 ), 2 );
                        }
                        if ( drelation == 2 ) {
                            delta = deltaCu * ( 1 - exp(-eCr / ( 0.1 * eCu ) ) ) / ( 1 - exp(-eCu / ( 0.1 * eCu ) ) );  //sigmoid delta relation
                        }
                    } else if ( eCr >= eCu && eCr <= eUl ) {
                        delta  = le * ( eCr - eCu ) + deltaCu / alpha;
                        Ddelta = le;
                    }

                    // Chao: I have put initialisation to start of function
                    //  double fibre = 0.;
                    if ( eCr > 0 && eCr <= eStar ) {
                        if ( fdtype == 0 ) {
                            fibre = ( 2. / k * ( ( 1. - acosh(1. + lamda * delta / deltaStar) / k ) * sqrt(pow( ( 1. + lamda * delta / deltaStar ), 2. ) - 1.) + ( lamda * delta ) / ( k * deltaStar ) ) + ( aa * pow( ( delta * 2. / lf ), 2 ) / 2. ) + ( b * delta * 2. / lf ) ) * s0;
                            dfibre = ( s0 * ( 1. / k * ( 1. - 1 / k * acosh(1 + lamda * delta / deltaStar) ) * 2. * lamda * ( 1 + lamda * delta / deltaStar ) / ( deltaStar * sqrt(pow( ( 1. + lamda * delta / deltaStar ), 2) - 1.) ) + aa * 2. * delta / lf + b ) ) * 2. / lf * Ddelta * kappaTwo;
                        }
                        if ( fdtype == 1 ) {
                            fibre = fibre_star * delta / deltaStar + ( bbb2 - ddd0 ) * pow(delta, 2) / deltaStar + ( ddd0 - bbb2 ) * delta;
                            dfibre = ( fibre_star / deltaStar + 2 * ( bbb2 - ddd0 ) * delta / deltaStar + ( ddd0 - bbb2 ) ) * Ddelta * kappaTwo;
                        }
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
                    residual        = ( 1. - omega ) * this->eM * equivStrain - fibre - concrete;
                    dResidualDOmega = -this->eM * equivStrain - dfibre - dconcrte;
                    omega -= residual / dResidualDOmega;
                    //printf("omega=%e,residual=%e\n,",omega,residual);
                } while ( fabs(residual / this->ft) >= 1.e-8 );
                //printf("\n");
            }
        }   else {
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

    // end of namespace
}
