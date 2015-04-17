#include "config.h"
#ifdef HYDRO

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "auxiliary.h"
#include "tree.h"

#define gammin  	(1.01)
#define gammax  	(10.0)
#define eps			(1e-6)
#define maxit		(50)
#define diffusion	(0.1)
#define dviscmax	(0.1)
#define drhomax		(0.2)
#define small_R     1.0e-20

/* factor = cell_volume_inverse change: 1 or 0.125 */
/* v[0-6]               0:rho 1:Pressure 234:velocity[012] 5:gamma[i] 6:Gamma */
/* backup_hvars[0-5]    0:rho 1:energy   234:momentum[012] 5:internal_energy */
/* f:flux[0-6+advected] 0:rho 123:momentum[012] 4:energy   5:internal_energy 6:velocity */
/* solver for interface between cell 1 and 2 */
/* stl/r[0-4] "l/r 0:rhohalf 1:velocityhalf 2:Pressurehalf 3:Gamma_polytropic 4:\gamma_{*S}" */
/* sta[0-3]   "avg 0:rhohalf 1:velocityhalf 2:Pressurehalf 3:\hat{\gamma}_{*S/F}" */
/* ref = rhoflux_rate/dx (speed)*/
        
extern int smooth_density_gradients;
int step;

void riemann( double stl[5], double str[5], double sta[4]) {
    /*see Colella and Glaz 1985*/
/*     F_2=p^{nu+1}_* */
/*     F_1=p^{nu}_* */
/*     F_0=p^{nu-1}_* */
/*     F[r,l]0 = initial st[r,l] */
/*     [0] == left [1] == right */
        int i;
        int indd, ind_r;
        double rho0[2], p0[2], u0[2], bgam0[2], rgam0[2], gam0[2] ;
	double q[2] ;
	double x2[2], x3[2], a[2], b[2], c[2] ;
	double u_0[2], u_1[2], p_0, p_1, xx[2], w2[2] ;
	double p_2, dev, uf,rhof,gamf, rho_s, u_s, p_s, bgam_s, gam_s;
	double a_s, b_s, c_s, w_s, rho, xx4, xx5;
	double a2, a3, fs;
        
	int iter;

	/* initial guess for secant iteration */
        rho0[0]  = stl[0];
	u0[0]	 = stl[1];
	p0[0]    = stl[2];
	bgam0[0]  = stl[3];
        gam0[0]   = stl[4];
        
        rho0[1]  = str[0];
	u0[1]	 = str[1];
	p0[1]    = str[2];
	bgam0[1]  = str[3];
        gam0[1]   = str[4];
        
        p_0     = 0.5 * ( p0[0] + p0[1] );
        
        for(i=0;i<2;i++){
            rgam0[i]  = 1.0 / bgam0[i];
            x2[i]    = 2.0 * gam0[i] * rgam0[i] - 1.0;
            x3[i]    = 0.5 * ( gam0[i] - 1.0 ) * x2[i];
            a[i]     = gam0[i] - x3[i];
            b[i]     = x3[i] * p0[i];
            c[i]     = x2[i] * p0[i];
            q[i]    = sqrt( p0[i] * rho0[i] * ((1.0+bgam0[i]) * p_0/p0[i] + bgam0[i] - 1.0 ));
            u_0[i]   = u0[i];
        }
	
	p_1	= max( small_R, (((u_0[0]-u_0[1])*q[0]*0.707106781 + p0[0])*q[1] + p0[1]*q[0]) / (q[0] + q[1]) ); //? like eq 19 q=Z
        
	/* Riemann solver - secant iterations for pressure */
        dev = eps+eps;
	iter = 0;
	while ( iter <= maxit && dev > eps ) {
            for(i=0;i<2;i++){
		/* inner iteration loop for ws eq 34 */
		xx[i] = ( a[i] * p_1 + b[i] ) / ( p_1 + c[i] ); 
		w2[i] = 1.0/sqrt(max(small_R, xx[i] * rho0[i] * (p_1 + p0[i]))); // w2=1/(W^nu_S) 
		/* calculating velocity eq 18 top*/
                if(i==0){u_1[i] = u0[i] + ( p0[i] - p_1 ) * w2[i];}
                if(i==1){u_1[i] = u0[i] - ( p0[i] - p_1 ) * w2[i];}
            }
	    /* iterating pressure eq 18 bottom*/
            p_2 = max( small_R, 1.0000001 * p_1 - ( u_1[1] - u_1[0] )
                       * fabs( p_1 - p_0 )
                       / (  fabs( u_1[1] - u_0[1] )
                           +fabs( u_1[0] - u_0[0] )
                            +small_R ) ); 
            dev = fabs( p_2 - p_1 ) / ( p_2 + p_1 );
            p_0 = p_1;
            p_1 = p_2;
            u_0[0] = u_1[0];
            u_0[1] = u_1[1];
            iter++;
	}
        
	if ( dev > eps ) {
            cart_error("Riemann solver did not converge!\ndev = %e, p_0 = %e, p_1 = %e\nLeft: %e %e %e %e %e\nRight: %e %e %e %e %e",
                       dev, p_0, p_1,
                       stl[0], stl[1], stl[2], stl[3], stl[4],
                       str[0], str[1], str[2], str[3], str[4] );
	}
        
	uf	= 0.5 * ( u_0[0] + u_0[1] );
        if( uf >=0 ){ ind_r = 0; } else { ind_r = 1;}
        rho_s      = rho0[ind_r];
        u_s	   = u0[ind_r];
        p_s	   = p0[ind_r];
        bgam_s     = bgam0[ind_r];
        gam_s      = gam0[ind_r];
        a_s	   = a[ind_r];
        b_s	   = b[ind_r];
        c_s	   = c[ind_r];
	w_s	= ( a_s * p_1 + b_s ) / ( p_1 + c_s );
	w_s	= max( small_R, w_s * rho_s * ( p_1 + p_s ) );
	
	rhof	= max( small_R, rho_s / ( 1.0 - rho_s * ( p_1 - p_s ) / w_s ) ); //eq 15
	gamf	= gam_s + 2.0 * ( gam_s - 1.0 ) * ( 1.0 - gam_s / bgam_s ) //eq 31 
			/ ( p_1 + p_s ) * ( p_1 - p_s );

 	indd	= sign(1.0,-uf) ;
        
	if ( p_1  >  p_s ) { //p_star > p_side : shock
            a2 = indd * u_s + sqrt( w_s ) / rho_s;
            a3 = 1.0e-10;
	} else {             //p_star < p_side : rarefaction
            xx4 = indd * u_s + sqrt( bgam_s * p_s/rho_s );
            xx5 = indd * uf + sqrt( 0.5 * ( bgam0[0]+bgam0[1] ) * p_1 / rhof );
            a2 = xx5; 
            a3 = xx4 - xx5; 
	}

	fs = -a2/a3;
        if ( fs < 0.0 ){ 
	    fs = 0.0;        //using star states
        }else if ( fs > 1.0 ){ 
	    fs = 1.0;        //using side states 
	}	             //interpolate between
        
        sta[0] = rhof   + fs * ( rho_s - rhof );
        sta[1] = uf     + fs * ( u_s - uf );
        sta[2] = p_1    + fs * ( p_s - p_1 );
        sta[3] = gamf   + fs * ( gam_s - gamf );
        
}       

void slope_limit(double center_value[4],double c[2],double dv0f[2],double dv[2]){
	double dv0, dv1, dv2, dv0a, dv1a, dv2a;
        double dv11, dv20, dlq0, dlq1;
		dv0 = center_value[1] - center_value[0];
		dv1 = center_value[2] - center_value[1];
		dv2 = center_value[3] - center_value[2];

		dv0a = 2.0*fabs(dv0);
		dv1a = 2.0*fabs(dv1);
		dv2a = 2.0*fabs(dv2);

		dv11 = center_value[2] - center_value[0];
		dv20 = center_value[3] - center_value[1];

		dlq0 = ( dv1*dv0 < 0.0 ) ? 0.0 : min( dv0a, dv1a );
		dlq1 = ( dv1*dv2 < 0.0 ) ? 0.0 : min( dv1a, dv2a );

		dv0f[0] = sign( min( 0.5*fabs(dv11), dlq0 ), dv11 ) * c[0];
		dv0f[1] = sign( min( 0.5*fabs(dv20), dlq1 ), dv20 ) * c[1];

                dv[0] = sign( min( min( fabs(dv0f[0]), dv1a ), dv0a ), dv11 );
                dv[1] = sign( min( min( fabs(dv0f[1]), dv2a ), dv1a ), dv20 );
}
#ifdef GRAVITY_IN_RIEMANN
void fluxh( double dtx, double dtx2, double v[num_hydro_vars-1][4], double g[4], double c[2], double f[num_hydro_vars-1])
#else
void fluxh( double dtx, double dtx2, double v[num_hydro_vars-1][4], double c[2], double f[num_hydro_vars-1] )
#endif
{
    /*see Colella and Glaz 1985*/
	int j;
	double rhor_l, u_l, a_l, C_l, C2_l, cp_l, cm_l, x_l, rhow_l;
	double pw_l, uw_l, vw_l, ww_l, gamw_l, b_l, x0_l, rho0_l;
	double v0_l, w0_l, gam0_l, b0_l, v_l, w_l, p_l, gam_l;
	double rhor_r, u_r, a_r, C_r, C2_r, cp_r, cm_r, x_r, rhow_r;
	double pw_r, uw_r, vw_r, ww_r, gamw_r, b_r, x0_r, rho0_r;
	double v0_r, w0_r, gam0_r, b0_r, v_r, w_r, p_r, gam_r;
	double rho, vu, pre, gam, xup_r, xup_l, vv, vw, fmass, predtx,pvudtx;
	double evudtx, vudtx;
	double fu, fv, fw, fe;
	double fl, fr;
	double dv[6][2];
	double dv0[2];
	double dv_other[2];
	double dg[2];
	double stl[5], str[5], sta[4];

	/* slopes */
#ifdef GRAVITY_IN_RIEMANN 
        //g[] adds velocity to the riemann problem. 
         slope_limit(g,c,dv0,dg); 
#endif
	for ( j = 0; j < 6; j++ ) {
            slope_limit(v[j],c,dv0,dv[j]);
        }

	/* left states */
	rhor_l	= 1.0/v[0][1];
	u_l	= v[2][1];
	a_l	= sqrt( v[6][1] * v[1][1] * rhor_l ); //sound speed
	C_l	= a_l * v[0][1]; //lagrangian soundspeed 
	C2_l	= C_l*C_l;
	cp_l	= u_l + a_l;  //characteristic lambda+
	cm_l	= u_l - a_l;  //characteristic lambda-

	/* interpolated states eq 44 */
	x_l	= 0.5 * (1.0 - dtx * max( (float)cp_l, 0.0 )); 
 	rhow_l	= max( (float)(v[0][1] + x_l * dv[0][0]), small_R ); 
        pw_l	= max( (float)(v[1][1] + x_l * dv[1][0]), small_R );
	uw_l	=              v[2][1] + x_l * dv[2][0];
	vw_l	=              v[3][1] + x_l * dv[3][0];
	ww_l	=              v[4][1] + x_l * dv[4][0];
	gamw_l	=              v[5][1] + x_l * dv[5][0];
        
	x0_l	= 0.5 - dtx2 * v[2][1];
	rho0_l	= v[0][1] + x0_l * dv[0][0];
	v0_l 	= v[3][1] + x0_l * dv[3][0];
	w0_l	= v[4][1] + x0_l * dv[4][0];
	gam0_l	= v[5][1] + x0_l * dv[5][0];
        
	/* characteristics eq 45 */
	b_l	= ( cm_l < 0.0 ) ? 0.0 : - dtx2 * rhor_l * ( dv[2][0] - dv[1][0] / C_l  );  //beta_L
	if ( v[2][1] > 0 ) {
            b0_l  = - a_l * dtx2 * (  dv[1][0] / C2_l -  dv[0][0] / ( rhow_l * rho0_l ) );
                                                  
            v_l	  = v0_l;
            w_l	  = w0_l;
            gam_l = ( gam0_l + 2.0 * ( 1.0 - v[5][1] / v[6][1] )
			        	       * ( v[5][1] - 1.0 )
				               * ( p_l - v[1][1] )
                                               / ( p_l + v[1][1] ) );
        }else{
            b0_l  = 0;
            v_l	  = vw_l;
            w_l	  = ww_l;
            gam_l = gamw_l;
        }

	/* states reconstructed from characteristics at x=i+1/2, t=n+1/2  eq 43/46*/
	stl[0]	= max( small_R, (float)(rhow_l / ( 1.0 - ( b0_l + b_l ) * rhow_l ) ) );
#ifdef GRAVITY_IN_RIEMANN 
 	stl[1]	= uw_l - b_l * C_l + (g[1] + x_l*dg[0]);  
#else 
	stl[1]	= uw_l - b_l * C_l;
#endif
	p_l	= max( small_R, (pw_l + b_l * C2_l) );
 	stl[2]	= p_l;
	stl[3]	= v[6][1];
	stl[4]	= max( gammin, min( gammax, (float)(gam_l) ) );

	/* right states */
	rhor_r  = 1.0/v[0][2];
	u_r     = v[2][2];
	a_r     = sqrt( v[6][2] * v[1][2] * rhor_r );
	C_r     = a_r * v[0][2];
	C2_r    = C_r*C_r;
	cp_r    = u_r + a_r;
	cm_r    = u_r - a_r;
        
	/* interpolated states eq 44 */
	x_r     = -(float)0.5*(1.0 + dtx * min( (float)cm_r, 0.0 ));
	rhow_r  = max( small_R, (float)(v[0][2] + x_r*dv[0][1]));
	pw_r    = max( small_R, (float)(v[1][2] + x_r*dv[1][1]));
	uw_r    = v[2][2] + x_r * dv[2][1];
	vw_r    = v[3][2] + x_r * dv[3][1];
	ww_r    = v[4][2] + x_r * dv[4][1];
	gamw_r  = v[5][2] + x_r * dv[5][1];

	x0_r    = -(float)0.5*(1.0 + dtx * v[2][2]);
	rho0_r  = v[0][2] + x0_r * dv[0][1];
	v0_r    = v[3][2] + x0_r * dv[3][1];
	w0_r    = v[4][2] + x0_r * dv[4][1];
	gam0_r  = v[5][2] + x0_r * dv[5][1];

	/* characteristics eq 45 */
	b_r	= ( cp_r >= 0.0 ) ? 0.0 : - dtx2 * rhor_r * ( dv[2][1] + dv[1][1] / C_r );
        if( v[2][2] < 0 ){
            b0_r    = a_r * dtx2 * (   dv[1][1] / C2_r - dv[0][1] / ( rhow_r * rho0_r ) );
					
            v_r  = v0_r ; 
            w_r  = w0_r ; 
            gam_r   = ( gam0_r + 2.0 * ( 1.0 - v[5][2] / v[6][2] )
                                     * ( v[5][2] - 1.0 )
                                     * ( p_r - v[1][2] )
                                     / ( p_r + v[1][2] ) );
        }else{
            b0_r    = 0;
            v_r     = vw_r; 
            w_r     = ww_r; 
            gam_r   = gamw_r;
        }

	/* states reconstructed from characteristics at x=i+1/2, t=n+1/2  eq 43/46*/
	str[0]  = max( small_R,(float)(rhow_r / ( 1.0 - ( b0_r + b_r ) * rhow_r ) ) );
#ifdef GRAVITY_IN_RIEMANN
	str[1]	= uw_r + b_r * C_r + (g[2] + x_r*dg[1]); 
#else 
	str[1]	= uw_r + b_r * C_r;
#endif 
//	if(b_r != 0){ cart_debug("=====================;br>0 b_rC=%e,u=%e",b_r * C_r,uw_r);}
	p_r     = max( small_R,(float)(pw_r + b_r * C2_r) );
	str[2]  = p_r;
	str[3]  = v[6][2];
	str[4]  = max( gammin, min( gammax, (float)(gam_r) ) );

	/* call to riemann function */
	riemann( stl, str, sta );

	/* compute fluxes of hydro variables */
	rho	= sta[0];
	vu	= sta[1];
	pre	= sta[2];
	gam	= sta[3];

	if ( vu < 0.0 ) {
		vv 	= v_r;
		vw 	= w_r;
	        xup_r   = 1.0;
		xup_l   = 0.0;
	} else {
		vv 	= v_l;
		vw	= w_l;
		xup_r	= 0.0;
		xup_l	= 1.0;
	}

	vudtx	= vu * dtx;
	fmass	= rho*vudtx;
	predtx	= pre * dtx;
        pvudtx  = vu * predtx;
	fu	= fmass * vu + predtx;
	fv	= fmass * vv;
	fw	= fmass * vw;
	evudtx	= pvudtx / ( gam - 1.0 );
  	fe	= evudtx + pvudtx + 0.5 * fmass * ( vu*vu + vv*vv + vw*vw ); 
        
	f[0]	= fmass;
	f[1]	= fu;
	f[2]	= fv;
	f[3]	= fw;
	f[4]	= fe;
	f[5]	= evudtx;
	f[6]	= vu;

	j=7;
#ifdef ELECTRON_ION_NONEQUILIBRIUM 
        slope_limit(v[j],c,dv0,dv_other);  
	fl = v[j][1] + ( 0.5 - dtx2 * vu ) * dv0[0] ;
	fr = v[j][2] - ( 0.5 + dtx2 * vu ) * dv0[1] ;
	f[j] = vudtx * ( fl * xup_l + fr * xup_r );
	j++;
#endif /* ELECTRON_ION_NONEQUILIBRIUM */
#ifdef TURBULENT_ENERGY
#error is this what you want for the pdv+advection? energy is an advected variable
	//reconstruct left and right state
	//x_l	  =  ( 0.5 - dtx2 * cp_l);
	//x_r     = -( 0.5 + dtx2 * cm_r)
	//v[2][2] + x_r * dv[2][1];
	//v[2][1] + x_l * dv[2][0];
        slope_limit(v[j],c,dv0,dv_other);  
	fl = v[j][1] + ( 0.5 - dtx2 * vu ) * dv0[0] ;  //state reconstruction on left
	fr = v[j][2] - ( 0.5 + dtx2 * vu ) * dv0[1] ;  //state reconstruction on right vu tells which direction to get flux from
	f[j] = vudtx * ( fl * xup_l + fr * xup_r ) + vudtx * (gamma_turb[ig]-1)/2.*(fl+fr) ;    
	j++;
	//E^{n+1}_i =  E^n_i +  vudtx_{i-0.5}*E_{i-0.5} -vudtx_{i+0.5}*E_{i+0.5} )  <-- above you're looking at fluxes across 1 interface, but if you shift to cell-centered and add 2 fluxes
	//P* is the same on each side it is the average of the reconstructions or just:
	// f[j] = vudtx*() + vu*Pa <--- BUT pre is not the total hydro pressure it is the cell_centered pressure
	// just say it is the average of the edge reconstructed pressure at time n so as to maintain monotonicity
#endif /* TURBULENT_ENERGY */

	/* compute fluxes of advected species */
	for ( j = num_hydro_vars-num_chem_species-1; j < num_hydro_vars - 1; j++ ) {
                slope_limit(v[j],c,dv0,dv_other);  
                     
		fl = v[j][1] + ( 0.5 - dtx2 * vu ) * dv0[0]; 
		fr = v[j][2] - ( 0.5 + dtx2 * vu ) * dv0[1]; 

		/* advected variables in v[j] must be in dimensionless units (var/density) */
		f[j] = fmass * ( fl * xup_l + fr * xup_r );
	}
}

void lapidus( double dtx2, int L1, int R1, int sweep_direction, int j3, int j4, int j5, 
		double v[num_hydro_vars-1][4], double f[num_hydro_vars-1] ) {
	int i, j;
	double diffk;
	double diff;
	double gvisc;
	double dvisc;
	int neighborsL1[num_neighbors];
	int neighborsR1[num_neighbors];
	int iCh00, iCh01, iCh10, iCh11;
	double v00, v01, v10, v11;
	double xx;
	int j0, j1;

	diffk = dtx2*diffusion;

	gvisc = 2.0 * ( v[2][1] - v[2][2] );

	/* Compute neighbors */
	cell_all_neighbors( L1, neighborsL1 );
	cell_all_neighbors( R1, neighborsR1 );

	/* velocity difference in directions orthogonal to sweep direction */
	for ( i = 0; i < nDim; i++ ) {
		j0 = 2*i;
		j1 = j0+1;

		if ( j0 != sweep_direction && j1 != sweep_direction ) {
			iCh00 = neighborsL1[j0];
			iCh01 = neighborsL1[j1];
			iCh10 = neighborsR1[j0];
			iCh11 = neighborsR1[j1];

			v00 = cell_momentum(iCh00,i) / cell_gas_density(iCh00);
			v01 = cell_momentum(iCh01,i) / cell_gas_density(iCh01);
			v10 = cell_momentum(iCh10,i) / cell_gas_density(iCh10);
			v11 = cell_momentum(iCh11,i) / cell_gas_density(iCh11);

			gvisc += v10 - v11 + v00 - v01;
		}
	}
		
	diff = diffk * max( 0.0, gvisc );

	if(smooth_density_gradients)
	  {
	    xx = drhomax * max( v[0][1], v[0][2] ) / ( min( v[0][1], v[0][2] ) );
	    dvisc = max( 0.0, dviscmax * ( xx - 1.0 ) / ( xx + 1.0 ) );
	    diff = max( diff, dvisc );
	  }

	f[0] += diff * ( cell_gas_density(L1) - cell_gas_density(R1) );
	f[1] += diff * ( cell_momentum(L1,j3) - cell_momentum(R1,j3) );
	f[2] += diff * ( cell_momentum(L1,j4) - cell_momentum(R1,j4) );
	f[3] += diff * ( cell_momentum(L1,j5) - cell_momentum(R1,j5) );
	f[4] += diff * ( cell_gas_energy(L1) - cell_gas_energy(R1) );
	f[5] += diff * ( cell_gas_internal_energy(L1) - cell_gas_internal_energy(R1) );

#ifdef ELECTRON_ION_NONEQUILIBRIUM
	f[7] += diff * ( cell_electron_internal_energy(L1) - cell_electron_internal_energy(R1) );
#endif /* ELECTRON_ION_NONEQUILIBRIUM */

	for ( j = 0; j < num_chem_species; j++ ) {
		f[j+num_hydro_vars-num_chem_species-1] += diff * ( cell_advected_variable(L1,j) - cell_advected_variable(R1,j) );
	}
}
#endif /* HYDRO */
