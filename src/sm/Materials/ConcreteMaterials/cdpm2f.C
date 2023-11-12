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
    {}


    //   ********************************
    //   *** CLASS CDPM2F ***
    //   ********************************

#define IDM_ITERATION_LIMIT 1.e-8

    CDPM2F::CDPM2F(int n, Domain *d) :
        ConcreteDPM2(n, d)
    {}


    void CDPM2F::initializeFrom(InputRecord &ir)
    {
        // call the corresponding service for the linear elastic material
        ConcreteDPM2::initializeFrom(ir);

        // Chao: your input

        IR_GIVE_FIELD(ir, lf, _IFT_CDPM2F_Lf);

	
        IR_GIVE_FIELD(ir, vf0, _IFT_CDPM2F_Vf0);

	
        IR_GIVE_FIELD(ir, df, _IFT_CDPM2F_Df);

	
        IR_GIVE_FIELD(ir, ef, _IFT_CDPM2F_Ef);

	//Default 1 MPa
        IR_GIVE_FIELD(ir, tau0, _IFT_CDPM2F_Tau0);

	//Default 0.015
        IR_GIVE_FIELD(ir, beta, _IFT_CDPM2F_Beta);

	//default 0.8
        IR_GIVE_FIELD(ir, f, _IFT_CDPM2F_f);

	//default not available.  Theroetical value very small 
        IR_GIVE_FIELD(ir, sm, _IFT_CDPM2F_Sm);

	//default 1 
        IR_GIVE_FIELD(ir, alpha, _IFT_CDPM2F_Alpha);

	//double em;
        IR_GIVE_FIELD(ir, this->em, _IFT_IsotropicLinearElasticMaterial_e);

	//default 0.3
	IR_GIVE_OPTIONAL_FIELD(ir, this->t, _IFT_CDPM2F_sigmoidRatio);
	

	//Introduce reduction in vf due to dispersion
        this->vf = ( 1. + log(this->alpha) ) * this->vf0;

	//OK
        this->eta = this->ef * this->vf / ( this->em * ( 1. - this->vf ) );


	//OK
        this->g = 2. * ( 1. + exp(M_PI * this->f / 2.) ) / ( 4. + pow(this->f,2.) );

	//OK
        this->s0 = 1. / 2. * this->g * this->tau0 * this->vf * ( 1. + this->eta ) * this->lf / this->df;

	//OK
        this->omega = sqrt(4. * ( 1. + this->eta ) * this->beta * this->tau0 / this->ef);

	//OK
        this->k = this->omega * this->lf / ( 2. * this->df );

	//OK
        this->lambda = cosh(this->k) - 1.;

	//OK
        this->deltaStar = ( 2. * this->df ) / this-> beta * this->lambda;


	//OK
	this->c = this->beta * this->lf / ( 2. * this->df );

	//softening starts at end of debonding
	if(this->c <= 6*this->lambda + 2){	  
	  this->deltaCu = this->lf*this->lambda/this->c;
	}
	else{//softening starts later
	  this->deltaCu = ( 0.5 * this->lf * ( this->c - 2. ) / ( 3. * this->c ) );
	}


	//softening starts at end of debonding
	if(this->c <= 6*this->lambda + 2){	  
	  this->stressCu = (1.+2.*this->lambda)*(2.*lambda/this->c-1.);
	}
	else{
	  this->stressCu = 4.*(c+1.)/27./pow(c,3.);
	}
	
		
	//Calculate alphamin to check alpha input.
	
        double ftTemp = this->ft * ( 1. - this->yieldTolDamage );

	//Chao: What is this? Where is this described?
	// z= ftTemp/(sigma_cu/(Vfm*eta))
        this->z = ftTemp / ( 0.5 * g * tau0 * lf / df * ( ( 1. + beta * deltaCu / df ) * ( pow( ( 1. - 2. * deltaCu / lf ), 2. ) ) ) );

	//Chao: again, not clear where this is from
	//Needs to be derived
        this->vfm = ( -( this->em - z ) + sqrt(pow( ( this->em - z ), 2 ) + 4. * ( this->ef - this->em ) * this->z * this->em) ) / ( 2. * ( this->ef - this->em ) );

	//Chao: Where are the equations defined?
        this->alphaMin = exp( ( this->vfm - this->vf0 ) / this->vf0);

	this->deltaUl = this->lf/2;
	
        if ( alpha < alphaMin ) {
            OOFEM_ERROR("alpha should be larger than alphamin %e\n", alphaMin);
            printf("alphamin= %e\n", alphaMin);
        }

    }


    double CDPM2F::computeFibreStress(double delta) const
    {
      double fibreStress = 0.;


      if ( delta >= 0 && delta <= this->deltaStar ) {
	fibreStress = ( 2. / this->k * ( ( 1. - acosh(1. + this->lambda * delta / this->deltaStar) / this->k ) * sqrt(pow( ( 1. + this->lambda * delta / this->deltaStar ), 2. ) - 1.) + ( this->lambda * delta ) / ( this->k * this->deltaStar ) ) + this->ap * delta / this->deltaStar ) * this->s0;                     
      }
      else if ( delta > this->deltaStar && delta <= this->deltaUl ) {
	fibreStress = ( 1. + beta * delta / df ) * ( pow( ( 1. - 2. * delta / lf ), 2. ) ) * s0;
	
	if ( fibreStress <= 0. ) {
	  fibreStress = 0.;
	}
	
      } else if ( delta > this->deltaUl || delta < 0 ) {
	fibreStress = 0.;
      }

      return fibreStress;  
    }


    double CDPM2F::computeMatrixStress(double delta) const
    {
      double matrixStress = 0.;
      double ftTemp = this->ft * ( 1. - this->yieldTolDamage );

      
      if ( this->softeningType == 2 ) {
	matrixStress = ( 1 - this->vf ) * ftTemp * exp(- delta / this->wf);
      }
      else{
	   OOFEM_ERROR("concrete softening must be exponential (stype = 2)");
      }
      
      return matrixStress;
      
    }
    
    double CDPM2F::computeCrackOpening(double crackingStrain, const double le) const
    {
      double delta = 0., gammaR0 = 0.;

      double residual = 1;
      double dResidualDDelta = 0.;
      
      double ftTemp = this->ft * ( 1. - this->yieldTolDamage );

      //Why would gammaC be a function of kappaOne?
      double gammaCu = ( 1. - this->alpha ) * ftTemp/this->eM * this->sm / ( deltaCu * ( 1. - this->alphaMin ) ) + ( ( this->alpha - this->alphaMin ) / ( 1. - this->alphaMin ) );

      //Our standard ecu
      double eCu = this->deltaCu * gammaCu / this->sm;
      double eUl = this->lf/2./le;
      
      double deltaCuUnloading = deltaCu*(gammaCu*le-sm)/(le-sm);

      if ( crackingStrain >= 0 && crackingStrain <= eCu ) {//pre-preak

	delta = deltaCu * ( 1. - exp(- crackingStrain / ( 0.1 ) ) ) / ( 1. - exp( - ( eCu ) / ( 0.1 ) ) );//sigmoid delta relation

      } else if ( crackingStrain > eCu && crackingStrain <= eUl ) {
	//initial guess of delta
	delta = this->deltaCu;
	
	//Two cases: 1) Unloading occurs in the element. 2) No unloading occurs.
	if(le>gammaCu*sm){
	  //Case 1: Unloading in element occurs.
	  //Apply Newton method to solve third order equation. We could also use analytical soltioon
	  delta = 0.;
	  do{
	  	  
	  residual = 1./le*(delta + (le/this->sm - 1)*deltaCuUnloading/this->stressCu*this->s0*(1. + this->beta*delta/this->lf)*pow(1.-2.*delta/lf,2.));

	  dResidualDDelta =
	    1./le + 1./le*(le/sm-1.)*deltaCuUnloading/stressCu*this->s0*this->beta/lf*
	    pow(1.-2.*delta/lf,2.) + 1./le*(le/sm-1.)*deltaCuUnloading/stressCu*
	    this->s0*(1.+beta*delta/lf)*2.*(1.-2.*delta/lf)*2./lf;

	  delta -= residual/dResidualDDelta;

	  }
	  while(fabs(residual)/eUl>1.e-6);
	  
	}
	else{
	  //Element so small so that no unloading occurs in element.

	  gammaR0 = gammaCu*le/sm;
	  
	  delta = 1./(4.*(gammaR0-1))*
	    sqrt(pow(lf,2.)*pow(gammaR0,2.)-4.*lf*deltaCu*gammaR0 -
		 8.*crackingStrain*le*lf*gammaR0 + 4*pow(deltaCu,2.) -
		 16*(1.-gammaR0)*crackingStrain*le*deltaCu) - 2.*deltaCu + lf*gammaR0;
	}	  
      }
      /* 	else{ */
	
	
      /* 	//delta = le * ( eCr - eCu ) + deltaCu; */
      /* 	x = -B / ( 3. * A ); */
      /* 	o = pow( ( pow( ( ( 27. * pow(A, 2) * ( D - eCr ) - 9. * A * B * C + 2 * pow(B, 3) ) / ( 54. * pow(A, 3) ) ), 2 ) + pow( ( ( 3 * A * C - pow(B, 2) ) / ( 9. * pow(A, 2) ) ), 3 ) ), 1. / 2 ); */
      /* 	Y = ( 0. - ( 27. * pow(A, 2) * ( D - eCr ) - 9. * A * B * C + 2 * pow(B, 3) ) / ( 54. * pow(A, 3) ) ) + o; */
      /* 	Z = ( ( 0. - ( 27. * pow(A, 2) * ( D - eCr ) - 9. * A * B * C + 2 * pow(B, 3) ) / ( 54. * pow(A, 3) ) ) - o ); */
      /* 	if ( Y <= 0 && Z <= 0 ) { */
      /* 	  y = -pow(-Y, 1. / 3); */
      /* 	  z = -pow(-Z, 1. / 3); */
      /* 	} else if ( Y > 0 && Z <= 0 )    { */
      /* 	  y = pow(Y, 1. / 3); */
      /* 	  z = -pow(-Z, 1. / 3); */
      /* 	} else if ( Y <= 0 && Z > 0 )    { */
      /* 	  y = -pow(-Y, 1. / 3); */
      /* 	  z = pow(Z, 1. / 3); */
      /* 	} else if ( Y > 0 && Z > 0 )    { */
      /* 	  y = pow(Y, 1. / 3); */
      /* 	  z = pow(Z, 1. / 3); */
      /* 	} */
	
      /* 	delta = x  + y + z; */
      /* } */
      return delta;      
    }
    
    
    double CDPM2F::computeStressResidual(double equivStrain, double omega, double kappaOne, double kappaTwo, double le) const
    {
      //calculate cracking strain
      double  crackingStrain  = kappaOne + omega * kappaTwo;
      
      double delta = computeCrackOpening(crackingStrain, le);

      double fibreStress = computeFibreStress(delta);
      
      double matrixStress = computeMatrixStress(delta);
      
      double residual = ( 1. - omega ) * this->eM * equivStrain - fibreStress - matrixStress;		   
      
      return residual;		      
    }
    
		   
    double
    CDPM2F::computeDamageParamTension(double equivStrain, double kappaOne, double kappaTwo, double le, double omegaOld, double rateFactor) const
    {

      //Check first if element length is small enough. This needs to be done here because le is caluclated from the principal directions at the onset of cracking.
      if(le > (7.*c-2.)*this->sm/(3.*c-6.)){
	OOFEM_ERROR("element size should be less than %e. Your element size is le = %e\n", (7.*c-2.)*this->sm/(3.*c-6.), le);
      }


      double omega = 0.;
      double fibre    = 0.;
      double concrete = 0.;

      // So that damage does not turn out to be negative if function is entered for equivstrains smaller than e0.
      double ftTemp = this->ft * ( 1. - yieldTolDamage );

      

      /* double eCu = 0., eCu2 = 0., eStar = 0., eUl = 0., residual = 0., a = 0., eCr = 0.00, delta = 0., e0 = 0., gammac = 0., p = 0., q = 0., A = 0., B = 0., C = 0, D = 0, x = 0., y = 0., z = 0., o = 0., Z = 0., Y = 0.; */

      /* 	//Chao: What is this? Should this not be ftTemp/eM? */
      /*   e0 = ftTemp / this->eM; */

      /* 	//Why would gammaC be a function of kappaOne? */
      /* 	gammac = ( 1 - alpha ) * e0 * sm / ( deltaCu * ( 1 - alphamin ) ) + ( ( alpha - alphamin ) / ( 1 - alphamin ) ); */

      /* 	//Our standard ecu */
      /* 	eCu = this->deltaCu * gammac / this->sm; */

      /* 	//for cubic law only */
      /*   eCu2 = this->deltaCu  / ( this->sm ); */


      /* 	p = 1. - ( 2. * deltaCu / lf ); */
      /*   q = 1. + ( beta * deltaCu / df ); */

      /*   A          =   ( eCu2 ) / ( q * pow(p, 2) ) * ( gammac * le - sm ) * 4. * beta / ( pow(lf, 2) * df * le ); */
      /*   B          =   ( eCu2 ) / ( q * pow(p, 2) ) * ( gammac * le - sm ) * 4. * ( df - beta * lf ) / ( df * pow(lf, 2) * le ); */
      /*   C          =   ( eCu2 ) / ( q * pow(p, 2) ) * ( gammac * le - sm ) * ( beta * lf - 4. * df ) / ( df * lf * le ) + 1. / le; */
      /*   D          =   ( eCu2 ) / ( q * pow(p, 2) ) * ( gammac * le - sm ) / le; */
        // p=(gammac*le-sm)*deltaCu/sm * 1./((1.+beta*deltaCu/df)* pow((1.-2*deltaCu/lf),2));
        //A=4.*beta*p/(lf*df*le);
        //B=p*(4.*df-4.*beta*lf)/(df* pow(lf,2)*le);
        //C= 1/le+(lf*beta-4*df)*p/(df*lf*le);
        //D = p/le;

        //double gammad = gammac;


	//Chao: According to the article he>sm and he < (7c-2)sm/(3c-6). Why did you use different expressions. Are the paper descriptions not correct? Please let us discuss.
	//Should sm be here the real sm or gamma*sm? Probably the real sm?

	/* he_max = ; */
	/* he_min = this->sm; */

	
	/* double  he_min = 0., he_max = 0. ; */

	/* double eCu2 = this->deltaCu  / this->sm ; */


	
        /* he_min = ( sm + 3. * q * pow(p, 2) * df * pow(lf, 2) * beta / ( eCu2 * ( pow( ( beta * lf + 2. * df ), 2 ) ) ) ) / gammac; */


	
        /* if ( le > he_min ) { */
        /*     OOFEM_ERROR("element size should  not be larger than he_max %e\n", he_min); */
        /*     printf("he_max= %e\n", he_min); */
        /* } */
        /* he_max = lf * sm / ( 2 * gammac * deltaCu ); */

        /* if ( le > he_max ) { */
        /*     OOFEM_ERROR("element size should not be larger than he_max %e\n", he_max); */
        /*     printf("he_max= %e\n", he_max); */
        /* } */


	/*The main idea of the function is to determine the damage variable, which provides equilibrium of forces. On the constitutive level this means that we have to find the damage variable, which gives a cracking strain which is used to calculate the fibre stress. The procedure is carried out in two steps.
	  1) Firstly, the crack opening is determined from the cracking strain. There are two regions which are separated, namely distributed and localised cracking.
	  2) Secondly, use this crack opening to calculate the fibre stress. Since these expressions are highly nonlinear and implicit, we need to use a bisection method for solving for the damage variable.	 
	  */

	/* //Chao: What is gammad?  */
        /* double gammad = gammac + ( 1 - gammac ) * ( eCu - eCr ) / eCu; */

        int nite = 0;
        if ( equivStrain > e0 * ( 1. - yieldTolDamage ) ) { // Check if damage should start

	  /* eStar = ( -0.1 * eCu * log(1 - deltaStar / ( deltaCu  ) * ( 1 - exp(-10.) ) ) ) * gammad; */

	  
	    
            /* eUl = 0.5; */
            /* double deltaStarS = deltaStar * 2. / lf; */
            /* double dd1        = 2. / k * ( ( 1 - 1 / k * acosh(1 + lambda) ) * sqrt(pow( ( 1 + lambda ), 2 ) - 1.) + lambda / k ) * s0; */
            /* double dd2        = ( 1. + beta * lf * deltaStarS / ( 2. * df ) ) * pow( ( 1 - deltaStarS ), 2 ) * s0; */
            /* double dd         = ( dd2 - dd1 ) / s0; */
            /* double ddd1       = ( 1. / k * ( 1. - 1. / k * acosh(1. + lambda) ) * 2. * lambda * ( 1. + lambda ) / ( deltaStarS * sqrt(pow( ( 1. + lambda ), 2 ) - 1.) ) ) * s0; */
            /* double ddd2       = ( beta * lf / ( 2. * df ) * ( 1. - 3. * deltaStarS ) - 2. ) * ( 1. - deltaStarS ) * s0; */
            /* double bbb2       = ddd2 * 2. / lf; */
            /* double ddd        = ( ddd2 - ddd1 ) / s0; */
            /* double aa         = 2. * ( ddd * deltaStarS - dd ) / pow(deltaStarS, 2); */
            /* double b          = ddd - aa * deltaStarS; */
            /* double fibre_star = ( 2. / k * ( ( 1. - acosh(1. + lambda * deltaStar / deltaStar) / k ) * sqrt(pow( ( 1. + lambda * deltaStar / deltaStar ), 2. ) - 1.) + ( lambda * deltaStar ) / ( k * deltaStar ) ) + ( aa * pow( ( deltaStar * 2. / lf ), 2 ) / 2. ) + ( b * deltaStar * 2. / lf ) ) * s0; */
            /* double ddd0       =  fibre_star   / deltaStar; */





            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            /* omega    = 1.;// initial guess */
            /* residual = 0.; */
            /* a        = 0.; */
            nite     = 0;


	    
	    double residual = 0.;
	    
	    double omegaOne = 0.;
	    double boundOne = computeStressResidual(equivStrain, omegaOne, kappaOne, kappaTwo, le);
	    
	    double omegaTwo = 1.;
	    double boundTwo = computeStressResidual(equivStrain, omegaTwo, kappaOne, kappaTwo, le);
	    
	    if( (boundOne < 0 && boundTwo <0) || (boundOne > 0 && boundTwo >0)){
	      OOFEM_ERROR("Bisection method will not work because solution is not bracketed. The two residuals are %e and %e\n", boundOne, boundTwo);
	    }
	    
	    do {
	      nite++;
	      if ( nite == 100000 ) {
		OOFEM_ERROR("In computeDamageTension: bisection method not converged");
	      }

	      omega = (omegaOne+omegaTwo)/2.;
	      residual = computeStressResidual(equivStrain, omega, kappaOne, kappaTwo, le);

	      if((residual <0 && boundOne <0) || (residual >0 && boundOne >0)){
		omegaOne = omega;
		boundOne = residual;
	      }
	      else if((residual < 0 && boundTwo < 0) || (residual >0 && boundTwo >0)){
		omegaTwo = omega;
		boundTwo = residual;
	      }
	    }
	    while (fabs(residual / this->ft) > 1.e-6);	      
	      
	      
		    /* //Define cracking strain in 1D tension */
                    /* crackingStrain  = kappaOne + omega * kappaTwo; */
		    /* computeResidual(crackingStrain); */
		    

		    /* delta = 0.; */


		    /* if ( eCr >= 0 && eCr <= eCu ) {//pre-preak */
		    /*   delta = ( deltaCu * ( 1 - exp(-eCr / ( 0.1 * eCu ) ) ) / ( 1 - exp(-( eCu ) / ( 0.1 * eCu ) ) ) );//sigmoid delta relation */
                    /* } else if ( eCr > eCu && eCr <= eUl ) { */
                    /*     //OOFEM_ERROR("In computeDamageTension: bisection method not converged"); */
                    /*     //delta = le * ( eCr - eCu ) + deltaCu; */
                    /*     x = -B / ( 3. * A ); */
                    /*     o = pow( ( pow( ( ( 27. * pow(A, 2) * ( D - eCr ) - 9. * A * B * C + 2 * pow(B, 3) ) / ( 54. * pow(A, 3) ) ), 2 ) + pow( ( ( 3 * A * C - pow(B, 2) ) / ( 9. * pow(A, 2) ) ), 3 ) ), 1. / 2 ); */
                    /*     Y = ( 0. - ( 27. * pow(A, 2) * ( D - eCr ) - 9. * A * B * C + 2 * pow(B, 3) ) / ( 54. * pow(A, 3) ) ) + o; */
                    /*     Z = ( ( 0. - ( 27. * pow(A, 2) * ( D - eCr ) - 9. * A * B * C + 2 * pow(B, 3) ) / ( 54. * pow(A, 3) ) ) - o ); */
                    /*     if ( Y <= 0 && Z <= 0 ) { */
                    /*         y = -pow(-Y, 1. / 3); */
                    /*         z = -pow(-Z, 1. / 3); */
                    /*     } else if ( Y > 0 && Z <= 0 )    { */
                    /*         y = pow(Y, 1. / 3); */
                    /*         z = -pow(-Z, 1. / 3); */
                    /*     } else if ( Y <= 0 && Z > 0 )    { */
                    /*         y = -pow(-Y, 1. / 3); */
                    /*         z = pow(Z, 1. / 3); */
                    /*     } else if ( Y > 0 && Z > 0 )    { */
                    /*         y = pow(Y, 1. / 3); */
                    /*         z = pow(Z, 1. / 3); */
                    /*     } */

                    /*     delta = x  + y + z; */
                    /* } */

	    

		   /*  //This is the jump to fix transition from debonding to pull out. */
		   /*  this->ap = ( 1. + beta * deltaStar / df ) * pow( ( 1. - deltaStar * 2 / lf ), 2 ) - 2 * lambda / ( pow(k, 2) ); */

		   /* double fibreStress = computeFibreStress(epscr); */

		   /* double concreteStress = computeMatrixStress(epscr); */


		   /* residual = ( 1. - omega ) * this->eM * equivStrain - fibreStress - concreteStress;		    */

                    // Chao: I have put initialisation to start of function
                    //  double fibre = 0.;

                    /* if ( eCr >= 0 && eCr <= eStar ) { */
                    /*     if ( fdtype == 0 ) { */
                    /*         fibre = ( 2. / k * ( ( 1. - acosh(1. + lambda * delta / deltaStar) / k ) * sqrt(pow( ( 1. + lambda * delta / deltaStar ), 2. ) - 1.) + ( lambda * delta ) / ( k * deltaStar ) ) + ( aa * pow( ( delta * 2. / lf ), 2 ) / 2. ) + ( b * delta * 2. / lf ) ) * s0; */
                    /*     } */
                    /*     if ( fdtype == 1 ) { */
                    /*         fibre = fibre_star * delta / deltaStar + ( bbb2 - ddd0 ) * pow(delta, 2) / deltaStar + ( ddd0 - bbb2 ) * delta; */
                    /*     } */
                    /*     if ( fdtype == 2 ) { */
                    /*         fibre = ( 2. / k * ( ( 1. - acosh(1. + lambda * delta / deltaStar) / k ) * sqrt(pow( ( 1. + lambda * delta / deltaStar ), 2. ) - 1.) + ( lambda * delta ) / ( k * deltaStar ) ) + ap * delta / deltaStar ) * s0; */
                    /*     } */
                    /* } else if ( eCr > eStar && eCr <= eUl ) { */
                    /*     fibre = ( 1. + beta * delta / df ) * ( pow( ( 1. - 2. * delta / lf ), 2. ) ) * s0; */
                    /*     if ( fibre <= 0. ) { */
                    /*         fibre = 0.; */
                    /*         //printf("omega = %e, residual = %e, fibre = %e, concrete = %e\n",omega, residual, fibre, concrete); */
                    /*     } */
                    /* } else if ( eCr > eUl || eCr < 0 ) { */
                    /*     fibre = 0.; */
                    /* } */

                    /* if ( this->softeningType == 2 ) { */
                    /*     concrete = ( 1 - vf ) * ftTemp * exp(-delta / wfMod); */
                    /* } */
		    /* else{ */
		    /*   OOFEM_ERROR("concrete sotening must be exponential (stype = 2)"); */
		    /* } */
		    

                    
                    /* if ( residual < 0 ) { */
                    /*     omega = ( omega + a ) / 2.; */
                    /* } else if ( residual > 0 ) { */
                    /*     double c = omega * 2. - a; */
                    /*     a        = omega; */
                /*     /\*     omega    = ( c + a ) / 2.; *\/ */
                /*     } */
                /* } while ( fabs(residual / this->ft) > 0.000001 ); */

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
