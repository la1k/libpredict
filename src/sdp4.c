#define _XOPEN_SOURCE 600
#include <math.h>
#include <stdbool.h>
#include "sdp4.h"
#include "defs.h"
#include "unsorted.h"

/// Entry points of deep()
#define DPInit		0
#define DPSecular	1
#define DPPeriodic	2

void sdp4_init(struct _sdp4 *m)
{
	m->initialized = 0;
	m->lunarTermsDone = 0;
	m->resonanceFlag = 0;
	m->synchronousFlag = 0;
	m->loopFlag = 0;
	m->epochRestartFlag = 0;

}

void sdp4_predict(struct _sdp4 *m, double tsince, predict_tle_t * tle, double pos[3], double vel[3])
{

	int i;
	double a, axn, ayn, aynl, beta, betal, capu, cos2u, cosepw, cosik,
	cosnok, cosu, cosuk, ecose, elsq, epw, esine, pl, theta4, rdot,
	rdotk, rfdot, rfdotk, rk, sin2u, sinepw, sinik, sinnok, sinu,
	sinuk, tempe, templ, tsq, u, uk, ux, uy, uz, vx, vy, vz, xl,
	xlt, xmam, xmdf, xmx, xmy, xnoddf, xll, a1, a3ovk2, ao, c2,
	coef, coef1, x1m5th, xhdot1, del1, r, delo, eeta, eta, etasq,
	perigee, psisq, tsi, qoms24, s4, pinvsq, temp, tempa, temp1,
	temp2, temp3, temp4, temp5, temp6;

	/* Initialization (if flag not set) */
	if (!m->initialized) {
		
		//Set initialized flag:
		m->initialized = true;

		/* Recover original mean motion (xnodp) and   */
		/* semimajor axis (aodp) from input elements. */
	  
		a1=pow(xke/tle->xno,tothrd);
		m->deep_arg.cosio=cos(tle->xincl);
		m->deep_arg.theta2=m->deep_arg.cosio*m->deep_arg.cosio;
		m->x3thm1=3*m->deep_arg.theta2-1;
		m->deep_arg.eosq=tle->eo*tle->eo;
		m->deep_arg.betao2=1-m->deep_arg.eosq;
		m->deep_arg.betao=sqrt(m->deep_arg.betao2);
		del1=1.5*ck2*m->x3thm1/(a1*a1*m->deep_arg.betao*m->deep_arg.betao2);
		ao=a1*(1-del1*(0.5*tothrd+del1*(1+134/81*del1)));
		delo=1.5*ck2*m->x3thm1/(ao*ao*m->deep_arg.betao*m->deep_arg.betao2);
		m->deep_arg.xnodp=tle->xno/(1+delo);
		m->deep_arg.aodp=ao/(1-delo);

		/* For perigee below 156 km, the values */
		/* of s and qoms2t are altered.         */
	  
		s4=s;
		qoms24=qoms2t;
		perigee=(m->deep_arg.aodp*(1-tle->eo)-ae)*xkmper;
	  
		if (perigee<156.0)
		{
			if (perigee<=98.0)
				s4=20.0;
			else
				s4=perigee-78.0;
	
			qoms24=pow((120-s4)*ae/xkmper,4);
			s4=s4/xkmper+ae;
		}

		pinvsq=1/(m->deep_arg.aodp*m->deep_arg.aodp*m->deep_arg.betao2*m->deep_arg.betao2);
		m->deep_arg.sing=sin(tle->omegao);
		m->deep_arg.cosg=cos(tle->omegao);
		tsi=1/(m->deep_arg.aodp-s4);
		eta=m->deep_arg.aodp*tle->eo*tsi;
		etasq=eta*eta;
		eeta=tle->eo*eta;
		psisq=fabs(1-etasq);
		coef=qoms24*pow(tsi,4);
		coef1=coef/pow(psisq,3.5);
		c2=coef1*m->deep_arg.xnodp*(m->deep_arg.aodp*(1+1.5*etasq+eeta*(4+etasq))+0.75*ck2*tsi/psisq*m->x3thm1*(8+3*etasq*(8+etasq)));
		m->c1=tle->bstar*c2;
		m->deep_arg.sinio=sin(tle->xincl);
		a3ovk2=-xj3/ck2*pow(ae,3);
		m->x1mth2=1-m->deep_arg.theta2;
		m->c4=2*m->deep_arg.xnodp*coef1*m->deep_arg.aodp*m->deep_arg.betao2*(eta*(2+0.5*etasq)+tle->eo*(0.5+2*etasq)-2*ck2*tsi/(m->deep_arg.aodp*psisq)*(-3*m->x3thm1*(1-2*eeta+etasq*(1.5-0.5*eeta))+0.75*m->x1mth2*(2*etasq-eeta*(1+etasq))*cos(2*tle->omegao)));
		theta4=m->deep_arg.theta2*m->deep_arg.theta2;
		temp1=3*ck2*pinvsq*m->deep_arg.xnodp;
		temp2=temp1*ck2*pinvsq;
		temp3=1.25*ck4*pinvsq*pinvsq*m->deep_arg.xnodp;
		m->deep_arg.xmdot=m->deep_arg.xnodp+0.5*temp1*m->deep_arg.betao*m->x3thm1+0.0625*temp2*m->deep_arg.betao*(13-78*m->deep_arg.theta2+137*theta4);
		x1m5th=1-5*m->deep_arg.theta2;
		m->deep_arg.omgdot=-0.5*temp1*x1m5th+0.0625*temp2*(7-114*m->deep_arg.theta2+395*theta4)+temp3*(3-36*m->deep_arg.theta2+49*theta4);
		xhdot1=-temp1*m->deep_arg.cosio;
		m->deep_arg.xnodot=xhdot1+(0.5*temp2*(4-19*m->deep_arg.theta2)+2*temp3*(3-7*m->deep_arg.theta2))*m->deep_arg.cosio;
		m->xnodcf=3.5*m->deep_arg.betao2*xhdot1*m->c1;
		m->t2cof=1.5*m->c1;
		m->xlcof=0.125*a3ovk2*m->deep_arg.sinio*(3+5*m->deep_arg.cosio)/(1+m->deep_arg.cosio);
		m->aycof=0.25*a3ovk2*m->deep_arg.sinio;
		m->x7thm1=7*m->deep_arg.theta2-1;

		/* initialize Deep() */
		sdp4_deep(m, DPInit,tle,&m->deep_arg);
	}
		
	/* Update for secular gravity and atmospheric drag */
	xmdf=tle->xmo+m->deep_arg.xmdot*tsince;
	m->deep_arg.omgadf=tle->omegao+m->deep_arg.omgdot*tsince;
	xnoddf=tle->xnodeo+m->deep_arg.xnodot*tsince;
	tsq=tsince*tsince;
	m->deep_arg.xnode=xnoddf+m->xnodcf*tsq;
	tempa=1-m->c1*tsince;
	tempe=tle->bstar*m->c4*tsince;
	templ=m->t2cof*tsq;
	m->deep_arg.xn=m->deep_arg.xnodp;

	/* Update for deep-space secular effects */
	m->deep_arg.xll=xmdf;
	m->deep_arg.t=tsince;

	sdp4_deep(m, DPSecular, tle, &m->deep_arg);

	xmdf=m->deep_arg.xll;
	a=pow(xke/m->deep_arg.xn,tothrd)*tempa*tempa;
	m->deep_arg.em=m->deep_arg.em-tempe;
	xmam=xmdf+m->deep_arg.xnodp*templ;

	/* Update for deep-space periodic effects */
	m->deep_arg.xll=xmam;

	sdp4_deep(m, DPPeriodic,tle,&m->deep_arg);

	xmam=m->deep_arg.xll;
	xl=xmam+m->deep_arg.omgadf+m->deep_arg.xnode;
	beta=sqrt(1-m->deep_arg.em*m->deep_arg.em);
	m->deep_arg.xn=xke/pow(a,1.5);

	/* Long period periodics */
	axn=m->deep_arg.em*cos(m->deep_arg.omgadf);
	temp=1/(a*beta*beta);
	xll=temp*m->xlcof*axn;
	aynl=temp*m->aycof;
	xlt=xl+xll;
	ayn=m->deep_arg.em*sin(m->deep_arg.omgadf)+aynl;

	/* Solve Kepler's Equation */
	capu=FMod2p(xlt-m->deep_arg.xnode);
	temp2=capu;
	i=0;

	do
	{
		sinepw=sin(temp2);
		cosepw=cos(temp2);
		temp3=axn*sinepw;
		temp4=ayn*cosepw;
		temp5=axn*cosepw;
		temp6=ayn*sinepw;
		epw=(capu-temp4+temp3-temp2)/(1-temp5-temp6)+temp2;
	  
		if (fabs(epw-temp2)<=e6a)
			break;

		temp2=epw;
	  
	} while (i++<10);
		
	/* Short period preliminary quantities */
	ecose=temp5+temp6;
	esine=temp3-temp4;
	elsq=axn*axn+ayn*ayn;
	temp=1-elsq;
	pl=a*temp;
	r=a*(1-ecose);
	temp1=1/r;
	rdot=xke*sqrt(a)*esine*temp1;
	rfdot=xke*sqrt(pl)*temp1;
	temp2=a*temp1;
	betal=sqrt(temp);
	temp3=1/(1+betal);
	cosu=temp2*(cosepw-axn+ayn*esine*temp3);
	sinu=temp2*(sinepw-ayn-axn*esine*temp3);
	u=AcTan(sinu,cosu);
	sin2u=2*sinu*cosu;
	cos2u=2*cosu*cosu-1;
	temp=1/pl;
	temp1=ck2*temp;
	temp2=temp1*temp;

	/* Update for short periodics */
	rk=r*(1-1.5*temp2*betal*m->x3thm1)+0.5*temp1*m->x1mth2*cos2u;
	uk=u-0.25*temp2*m->x7thm1*sin2u;
	m->xnodek=m->deep_arg.xnode+1.5*temp2*m->deep_arg.cosio*sin2u;
	m->xinck=m->deep_arg.xinc+1.5*temp2*m->deep_arg.cosio*m->deep_arg.sinio*cos2u;
	rdotk=rdot-m->deep_arg.xn*temp1*m->x1mth2*sin2u;
	rfdotk=rfdot+m->deep_arg.xn*temp1*(m->x1mth2*cos2u+1.5*m->x3thm1);

	/* Orientation vectors */
	sinuk=sin(uk);
	cosuk=cos(uk);
	sinik=sin(m->xinck);
	cosik=cos(m->xinck);
	sinnok=sin(m->xnodek);
	cosnok=cos(m->xnodek);
	xmx=-sinnok*cosik;
	xmy=cosnok*cosik;
	ux=xmx*sinuk+cosnok*cosuk;
	uy=xmy*sinuk+sinnok*cosuk;
	uz=sinik*sinuk;
	vx=xmx*cosuk-cosnok*sinuk;
	vy=xmy*cosuk-sinnok*sinuk;
	vz=sinik*cosuk;

	/* Position and velocity */
	pos[0] = rk*ux;
	pos[1] = rk*uy;
	pos[2] = rk*uz;
	vel[0] = rdotk*ux+rfdotk*vx;
	vel[1] = rdotk*uy+rfdotk*vy;
	vel[2] = rdotk*uz+rfdotk*vz;

	/* Phase in radians */
	m->phase=xlt-m->deep_arg.xnode-m->deep_arg.omgadf+twopi;
    
	if (m->phase<0.0)
		m->phase+=twopi;

	m->phase=FMod2p(m->phase);
}

/**
 * Calculates the Greenwich Mean Sidereal Time
 * for an epoch specified in the format used in the NORAD two-line 
 * element sets. 
 * It has been adapted for dates beyond the year 1999.
 * Reference:  The 1992 Astronomical Almanac, page B6. 
 * Modification to support Y2K. Valid 1957 through 2056.
 *
 * \param epoch TLE epoch
 * \param deep_arg Deep arg
 * \copyright GPLv2+
 **/
double ThetaG(double epoch, deep_arg_t *deep_arg)
{
	double year, day, UT, jd, TU, GMST, ThetaG;

	/* Modification to support Y2K */
	/* Valid 1957 through 2056     */

	day=modf(epoch*1E-3,&year)*1E3;

	if (year<57)
		year+=2000;
	else
		year+=1900;

	UT=modf(day,&day);
	jd=Julian_Date_of_Year(year)+day;
	TU=(jd-2451545.0)/36525;
	GMST=24110.54841+TU*(8640184.812866+TU*(0.093104-TU*6.2E-6));
	GMST=Modulus(GMST+secday*omega_E*UT,secday);
	ThetaG = 2*M_PI*GMST/secday;
	deep_arg->ds50=jd-2433281.5+UT;
	ThetaG=FMod2p(6.3003880987*deep_arg->ds50+1.72944494);

	return ThetaG;
}

void sdp4_deep(struct _sdp4 *m, int ientry, predict_tle_t * tle, deep_arg_t * deep_arg)
{
	/* This function is used by SDP4 to add lunar and solar */
	/* perturbation effects to deep-space orbit objects.    */

	double a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, ainv2, alfdp, aqnv,
	sgh, sini2, sinis, sinok, sh, si, sil, day, betdp, dalf, bfact, c,
	cc, cosis, cosok, cosq, ctem, f322, zx, zy, dbet, dls, eoc, eq, f2,
	f220, f221, f3, f311, f321, xnoh, f330, f441, f442, f522, f523,
	f542, f543, g200, g201, g211, pgh, ph, s1, s2, s3, s4, s5, s6, s7,
	se, sel, ses, xls, g300, g310, g322, g410, g422, g520, g521, g532,
	g533, gam, sinq, sinzf, sis, sl, sll, sls, stem, temp, temp1, x1,
	x2, x2li, x2omi, x3, x4, x5, x6, x7, x8, xl, xldot, xmao, xnddt,
	xndot, xno2, xnodce, xnoi, xomi, xpidot, z1, z11, z12, z13, z2,
	  z21, z22, z23, z3, z31, z32, z33, ze, zf, zm, /*zmo,*/ zn, zsing,
	zsinh, zsini, zcosg, zcosh, zcosi, delt=0, ft=0;

	switch (ientry)
	{
		case DPInit:  /* Entrance for deep space initialization */
		m->thgr=ThetaG(tle->epoch,deep_arg);
		eq=tle->eo;
		m->xnq=deep_arg->xnodp;
		aqnv=1/deep_arg->aodp;
		m->xqncl=tle->xincl;
		xmao=tle->xmo;
		xpidot=deep_arg->omgdot+deep_arg->xnodot;
		sinq=sin(tle->xnodeo);
		cosq=cos(tle->xnodeo);
		m->omegaq=tle->omegao;

		/* Initialize lunar solar terms */
		day=deep_arg->ds50+18261.5;  /* Days since 1900 Jan 0.5 */
	  
		if (day!=m->preep)
		{
			m->preep=day;
			xnodce=4.5236020-9.2422029E-4*day;
			stem=sin(xnodce);
			ctem=cos(xnodce);
			m->zcosil=0.91375164-0.03568096*ctem;
			m->zsinil=sqrt(1-m->zcosil*m->zcosil);
			m->zsinhl=0.089683511*stem/m->zsinil;
			m->zcoshl=sqrt(1-m->zsinhl*m->zsinhl);
			c=4.7199672+0.22997150*day;
			gam=5.8351514+0.0019443680*day;
			m->zmol=FMod2p(c-gam);
			zx=0.39785416*stem/m->zsinil;
			zy=m->zcoshl*ctem+0.91744867*m->zsinhl*stem;
			zx=AcTan(zx,zy);
			zx=gam+zx-xnodce;
			m->zcosgl=cos(zx);
			m->zsingl=sin(zx);
			m->zmos=6.2565837+0.017201977*day;
			m->zmos=FMod2p(m->zmos);
		    }

		  /* Do solar terms */
		  m->savtsn=1E20;
		  zcosg=zcosgs;
		  zsing=zsings;
		  zcosi=zcosis;
		  zsini=zsinis;
		  zcosh=cosq;
		  zsinh= sinq;
		  cc=c1ss;
		  zn=zns;
		  ze=zes;
		  /* zmo=m->zmos; */
		  xnoi=1/m->xnq;

		  /* Loop breaks when Solar terms are done a second */
		  /* time, after Lunar terms are initialized        */
	  
		for (;;)
		{
			/* Solar terms done again after Lunar terms are done */
			a1=zcosg*zcosh+zsing*zcosi*zsinh;
			a3=-zsing*zcosh+zcosg*zcosi*zsinh;
			a7=-zcosg*zsinh+zsing*zcosi*zcosh;
			a8=zsing*zsini;
			a9=zsing*zsinh+zcosg*zcosi*zcosh;
			a10=zcosg*zsini;
			a2=deep_arg->cosio*a7+deep_arg->sinio*a8;
			a4=deep_arg->cosio*a9+deep_arg->sinio*a10;
			a5=-deep_arg->sinio*a7+deep_arg->cosio*a8;
			a6=-deep_arg->sinio*a9+deep_arg->cosio*a10;
			x1=a1*deep_arg->cosg+a2*deep_arg->sing;
			x2=a3*deep_arg->cosg+a4*deep_arg->sing;
			x3=-a1*deep_arg->sing+a2*deep_arg->cosg;
			x4=-a3*deep_arg->sing+a4*deep_arg->cosg;
			x5=a5*deep_arg->sing;
			x6=a6*deep_arg->sing;
			x7=a5*deep_arg->cosg;
			x8=a6*deep_arg->cosg;
			z31=12*x1*x1-3*x3*x3;
			z32=24*x1*x2-6*x3*x4;
			z33=12*x2*x2-3*x4*x4;
			z1=3*(a1*a1+a2*a2)+z31*deep_arg->eosq;
			z2=6*(a1*a3+a2*a4)+z32*deep_arg->eosq;
			z3=3*(a3*a3+a4*a4)+z33*deep_arg->eosq;
			z11=-6*a1*a5+deep_arg->eosq*(-24*x1*x7-6*x3*x5);
			z12=-6*(a1*a6+a3*a5)+deep_arg->eosq*(-24*(x2*x7+x1*x8)-6*(x3*x6+x4*x5));
			z13=-6*a3*a6+deep_arg->eosq*(-24*x2*x8-6*x4*x6);
			z21=6*a2*a5+deep_arg->eosq*(24*x1*x5-6*x3*x7);
			z22=6*(a4*a5+a2*a6)+deep_arg->eosq*(24*(x2*x5+x1*x6)-6*(x4*x7+x3*x8));
			z23=6*a4*a6+deep_arg->eosq*(24*x2*x6-6*x4*x8);
			z1=z1+z1+deep_arg->betao2*z31;
			z2=z2+z2+deep_arg->betao2*z32;
			z3=z3+z3+deep_arg->betao2*z33;
			s3=cc*xnoi;
			s2=-0.5*s3/deep_arg->betao;
			s4=s3*deep_arg->betao;
			s1=-15*eq*s4;
			s5=x1*x3+x2*x4;
			s6=x2*x3+x1*x4;
			s7=x2*x4-x1*x3;
			se=s1*zn*s5;
			si=s2*zn*(z11+z13);
			sl=-zn*s3*(z1+z3-14-6*deep_arg->eosq);
			sgh=s4*zn*(z31+z33-6);
			sh=-zn*s2*(z21+z23);
		
			if (m->xqncl<5.2359877E-2)
				sh=0;
		    
			m->ee2=2*s1*s6;
			m->e3=2*s1*s7;
			m->xi2=2*s2*z12;
			m->xi3=2*s2*(z13-z11);
			m->xl2=-2*s3*z2;
			m->xl3=-2*s3*(z3-z1);
			m->xl4=-2*s3*(-21-9*deep_arg->eosq)*ze;
			m->xgh2=2*s4*z32;
			m->xgh3=2*s4*(z33-z31);
			m->xgh4=-18*s4*ze;
			m->xh2=-2*s2*z22;
			m->xh3=-2*s2*(z23-z21);

			//Skip lunar terms?
			if (m->lunarTermsDone) {
				break;
			}

			/* Do lunar terms */
			m->sse=se;
			m->ssi=si;
			m->ssl=sl;
			m->ssh=sh/deep_arg->sinio;
			m->ssg=sgh-deep_arg->cosio*m->ssh;
			m->se2=m->ee2;
			m->si2=m->xi2;
			m->sl2=m->xl2;
			m->sgh2=m->xgh2;
			m->sh2=m->xh2;
			m->se3=m->e3;
			m->si3=m->xi3;
			m->sl3=m->xl3;
			m->sgh3=m->xgh3;
			m->sh3=m->xh3;
			m->sl4=m->xl4;
			m->sgh4=m->xgh4;
			zcosg=m->zcosgl;
			zsing=m->zsingl;
			zcosi=m->zcosil;
			zsini=m->zsinil;
			zcosh=m->zcoshl*cosq+m->zsinhl*sinq;
			zsinh=sinq*m->zcoshl-cosq*m->zsinhl;
			zn=znl;
			cc=c1l;
			ze=zel;
			/* zmo=m->zmol; */
			//Set lunarTermsDone flag:
			m->lunarTermsDone = true;
		}

		m->sse=m->sse+se;
		m->ssi=m->ssi+si;
		m->ssl=m->ssl+sl;
		m->ssg=m->ssg+sgh-deep_arg->cosio/deep_arg->sinio*sh;
		m->ssh=m->ssh+sh/deep_arg->sinio;

		/* Geopotential resonance initialization for 12 hour orbits */
		m->resonanceFlag = 0;
		m->synchronousFlag = 0;

		if (!((m->xnq<0.0052359877) && (m->xnq>0.0034906585)))
		{
			if ((m->xnq<0.00826) || (m->xnq>0.00924))
			    return;
	
			if (eq<0.5)
			    return;
	
			m->resonanceFlag = 1;
			eoc=eq*deep_arg->eosq;
			g201=-0.306-(eq-0.64)*0.440;
		
			if (eq<=0.65)
			{
				g211=3.616-13.247*eq+16.290*deep_arg->eosq;
				g310=-19.302+117.390*eq-228.419*deep_arg->eosq+156.591*eoc;
				g322=-18.9068+109.7927*eq-214.6334*deep_arg->eosq+146.5816*eoc;
				g410=-41.122+242.694*eq-471.094*deep_arg->eosq+313.953*eoc;
				g422=-146.407+841.880*eq-1629.014*deep_arg->eosq+1083.435 * eoc;
				g520=-532.114+3017.977*eq-5740*deep_arg->eosq+3708.276*eoc;
			}
		
			else
			{
				g211=-72.099+331.819*eq-508.738*deep_arg->eosq+266.724*eoc;
				g310=-346.844+1582.851*eq-2415.925*deep_arg->eosq+1246.113*eoc;
				g322=-342.585+1554.908*eq-2366.899*deep_arg->eosq+1215.972*eoc;
				g410=-1052.797+4758.686*eq-7193.992*deep_arg->eosq+3651.957*eoc;
				g422=-3581.69+16178.11*eq-24462.77*deep_arg->eosq+12422.52*eoc;
		      
				if (eq<=0.715)
					g520=1464.74-4664.75*eq+3763.64*deep_arg->eosq;
			  
				else
					g520=-5149.66+29936.92*eq-54087.36*deep_arg->eosq+31324.56*eoc;
			}

			if (eq<0.7)
			{
				g533=-919.2277+4988.61*eq-9064.77*deep_arg->eosq+5542.21*eoc;
				g521=-822.71072+4568.6173*eq-8491.4146*deep_arg->eosq+5337.524*eoc;
				g532=-853.666+4690.25*eq-8624.77*deep_arg->eosq+5341.4*eoc;
			}
		
			else
			{
				g533=-37995.78+161616.52*eq-229838.2*deep_arg->eosq+109377.94*eoc;
				g521 =-51752.104+218913.95*eq-309468.16*deep_arg->eosq+146349.42*eoc;
				g532 =-40023.88+170470.89*eq-242699.48*deep_arg->eosq+115605.82*eoc;
			}

			sini2=deep_arg->sinio*deep_arg->sinio;
			f220=0.75*(1+2*deep_arg->cosio+deep_arg->theta2);
			f221=1.5*sini2;
			f321=1.875*deep_arg->sinio*(1-2*deep_arg->cosio-3*deep_arg->theta2);
			f322=-1.875*deep_arg->sinio*(1+2*deep_arg->cosio-3*deep_arg->theta2);
			f441=35*sini2*f220;
			f442=39.3750*sini2*sini2;
			f522=9.84375*deep_arg->sinio*(sini2*(1-2*deep_arg->cosio-5*deep_arg->theta2)+0.33333333*(-2+4*deep_arg->cosio+6*deep_arg->theta2));
			f523=deep_arg->sinio*(4.92187512*sini2*(-2-4*deep_arg->cosio+10*deep_arg->theta2)+6.56250012*(1+2*deep_arg->cosio-3*deep_arg->theta2));
			f542=29.53125*deep_arg->sinio*(2-8*deep_arg->cosio+deep_arg->theta2*(-12+8*deep_arg->cosio+10*deep_arg->theta2));
			f543=29.53125*deep_arg->sinio*(-2-8*deep_arg->cosio+deep_arg->theta2*(12+8*deep_arg->cosio-10*deep_arg->theta2));
			xno2=m->xnq*m->xnq;
			ainv2=aqnv*aqnv;
			temp1=3*xno2*ainv2;
			temp=temp1*root22;
			m->d2201=temp*f220*g201;
			m->d2211=temp*f221*g211;
			temp1=temp1*aqnv;
			temp=temp1*root32;
			m->d3210=temp*f321*g310;
			m->d3222=temp*f322*g322;
			temp1=temp1*aqnv;
			temp=2*temp1*root44;
			m->d4410=temp*f441*g410;
			m->d4422=temp*f442*g422;
			temp1=temp1*aqnv;
			temp=temp1*root52;
			m->d5220=temp*f522*g520;
			m->d5232=temp*f523*g532;
			temp=2*temp1*root54;
			m->d5421=temp*f542*g521;
			m->d5433=temp*f543*g533;
			m->xlamo=xmao+tle->xnodeo+tle->xnodeo-m->thgr-m->thgr;
			bfact=deep_arg->xmdot+deep_arg->xnodot+deep_arg->xnodot-thdt-thdt;
			bfact=bfact+m->ssl+m->ssh+m->ssh;
		}
	
		else
		{
			m->resonanceFlag = 1;
			m->synchronousFlag = 1;
	
			/* Synchronous resonance terms initialization */
			g200=1+deep_arg->eosq*(-2.5+0.8125*deep_arg->eosq);
			g310=1+2*deep_arg->eosq;
			g300=1+deep_arg->eosq*(-6+6.60937*deep_arg->eosq);
			f220=0.75*(1+deep_arg->cosio)*(1+deep_arg->cosio);
			f311=0.9375*deep_arg->sinio*deep_arg->sinio*(1+3*deep_arg->cosio)-0.75*(1+deep_arg->cosio);
			f330=1+deep_arg->cosio;
			f330=1.875*f330*f330*f330;
			m->del1=3*m->xnq*m->xnq*aqnv*aqnv;
			m->del2=2*m->del1*f220*g200*q22;
			m->del3=3*m->del1*f330*g300*q33*aqnv;
			m->del1=m->del1*f311*g310*q31*aqnv;
			m->fasx2=0.13130908;
			m->fasx4=2.8843198;
			m->fasx6=0.37448087;
			m->xlamo=xmao+tle->xnodeo+tle->omegao-m->thgr;
			bfact=deep_arg->xmdot+xpidot-thdt;
			bfact=bfact+m->ssl+m->ssg+m->ssh;
		}

		m->xfact=bfact-m->xnq;

		/* Initialize integrator */
		m->xli=m->xlamo;
		m->xni=m->xnq;
		m->atime=0;
		m->stepp=720;
		m->stepn=-720;
		m->step2=259200;

		return;

		case DPSecular:  /* Entrance for deep space secular effects */

		deep_arg->xll=deep_arg->xll+m->ssl*deep_arg->t;
		deep_arg->omgadf=deep_arg->omgadf+m->ssg*deep_arg->t;
		deep_arg->xnode=deep_arg->xnode+m->ssh*deep_arg->t;
		deep_arg->em=tle->eo+m->sse*deep_arg->t;
		deep_arg->xinc=tle->xincl+m->ssi*deep_arg->t;
	  
		if (deep_arg->xinc<0)
		{
			deep_arg->xinc=-deep_arg->xinc;
			deep_arg->xnode=deep_arg->xnode+pi;
			deep_arg->omgadf=deep_arg->omgadf-pi;
		}
	
		if (!m->resonanceFlag) {
			return;
		}

		do
		{
			if ((m->atime==0) || ((deep_arg->t>=0) && (m->atime<0)) || ((deep_arg->t<0) && (m->atime>=0)))
			{
				/* Epoch restart */

				if (deep_arg->t>=0)
					delt=m->stepp;
				else
					delt=m->stepn;

				m->atime=0;
				m->xni=m->xnq;
				m->xli=m->xlamo;
			}

			else
			{
				if (fabs(deep_arg->t)>=fabs(m->atime))
				{
					if (deep_arg->t>0)
						delt=m->stepp;
					else
						delt=m->stepn;
				}
			}
	    
			do
			{
				if (fabs(deep_arg->t-m->atime)>=m->stepp)
				{
					m->loopFlag = 1;
					m->epochRestartFlag = 0;
				}
		
				else
				{
					ft=deep_arg->t-m->atime;
					m->loopFlag = 0;
				}

				if (fabs(deep_arg->t)<fabs(m->atime))
				{
					if (deep_arg->t>=0)
						delt=m->stepn;
					else
						delt=m->stepp;

					m->loopFlag = 1;
					m->epochRestartFlag = 1;
				}

				/* Dot terms calculated */
				if (m->synchronousFlag) {
					xndot=m->del1*sin(m->xli-m->fasx2)+m->del2*sin(2*(m->xli-m->fasx4))+m->del3*sin(3*(m->xli-m->fasx6));
					xnddt=m->del1*cos(m->xli-m->fasx2)+2*m->del2*cos(2*(m->xli-m->fasx4))+3*m->del3*cos(3*(m->xli-m->fasx6));
				}
		
				else
				{
					xomi=m->omegaq+deep_arg->omgdot*m->atime;
					x2omi=xomi+xomi;
					x2li=m->xli+m->xli;
					xndot=m->d2201*sin(x2omi+m->xli-g22)+m->d2211*sin(m->xli-g22)+m->d3210*sin(xomi+m->xli-g32)+m->d3222*sin(-xomi+m->xli-g32)+m->d4410*sin(x2omi+x2li-g44)+m->d4422*sin(x2li-g44)+m->d5220*sin(xomi+m->xli-g52)+m->d5232*sin(-xomi+m->xli-g52)+m->d5421*sin(xomi+x2li-g54)+m->d5433*sin(-xomi+x2li-g54);
					xnddt=m->d2201*cos(x2omi+m->xli-g22)+m->d2211*cos(m->xli-g22)+m->d3210*cos(xomi+m->xli-g32)+m->d3222*cos(-xomi+m->xli-g32)+m->d5220*cos(xomi+m->xli-g52)+m->d5232*cos(-xomi+m->xli-g52)+2*(m->d4410*cos(x2omi+x2li-g44)+m->d4422*cos(x2li-g44)+m->d5421*cos(xomi+x2li-g54)+m->d5433*cos(-xomi+x2li-g54));
				}

				xldot=m->xni+m->xfact;
				xnddt=xnddt*xldot;

				if (m->loopFlag) {
					m->xli=m->xli+xldot*delt+xndot*m->step2;
					m->xni=m->xni+xndot*delt+xnddt*m->step2;
					m->atime=m->atime+delt;
				}
			} while (m->loopFlag && !m->epochRestartFlag);
		} while (m->loopFlag && m->epochRestartFlag);

		deep_arg->xn=m->xni+xndot*ft+xnddt*ft*ft*0.5;
		xl=m->xli+xldot*ft+xndot*ft*ft*0.5;
		temp=-deep_arg->xnode+m->thgr+deep_arg->t*thdt;

		if (!m->synchronousFlag) {
			deep_arg->xll=xl+temp+temp;
		}else{
			deep_arg->xll=xl-deep_arg->omgadf+temp;
		}

		return;

		case DPPeriodic:	 /* Entrance for lunar-solar periodics */
		sinis=sin(deep_arg->xinc);
		cosis=cos(deep_arg->xinc);

		if (fabs(m->savtsn-deep_arg->t)>=30)
		{
			m->savtsn=deep_arg->t;
			zm=m->zmos+zns*deep_arg->t;
			zf=zm+2*zes*sin(zm);
			sinzf=sin(zf);
			f2=0.5*sinzf*sinzf-0.25;
			f3=-0.5*sinzf*cos(zf);
			ses=m->se2*f2+m->se3*f3;
			sis=m->si2*f2+m->si3*f3;
			sls=m->sl2*f2+m->sl3*f3+m->sl4*sinzf;
			m->sghs=m->sgh2*f2+m->sgh3*f3+m->sgh4*sinzf;
			m->shs=m->sh2*f2+m->sh3*f3;
			zm=m->zmol+znl*deep_arg->t;
			zf=zm+2*zel*sin(zm);
			sinzf=sin(zf);
			f2=0.5*sinzf*sinzf-0.25;
			f3=-0.5*sinzf*cos(zf);
			sel=m->ee2*f2+m->e3*f3;
			sil=m->xi2*f2+m->xi3*f3;
			sll=m->xl2*f2+m->xl3*f3+m->xl4*sinzf;
			m->sghl=m->xgh2*f2+m->xgh3*f3+m->xgh4*sinzf;
			m->sh1=m->xh2*f2+m->xh3*f3;
			m->pe=ses+sel;
			m->pinc=sis+sil;
			m->pl=sls+sll;
		}

		pgh=m->sghs+m->sghl;
		ph=m->shs+m->sh1;
		deep_arg->xinc=deep_arg->xinc+m->pinc;
		deep_arg->em=deep_arg->em+m->pe;

		if (m->xqncl>=0.2)
		{
			/* Apply periodics directly */
			ph=ph/deep_arg->sinio;
			pgh=pgh-deep_arg->cosio*ph;
			deep_arg->omgadf=deep_arg->omgadf+pgh;
			deep_arg->xnode=deep_arg->xnode+ph;
			deep_arg->xll=deep_arg->xll+m->pl;
		}
	
		else
		{
			/* Apply periodics with Lyddane modification */
			sinok=sin(deep_arg->xnode);
			cosok=cos(deep_arg->xnode);
			alfdp=sinis*sinok;
			betdp=sinis*cosok;
			dalf=ph*cosok+m->pinc*cosis*sinok;
			dbet=-ph*sinok+m->pinc*cosis*cosok;
			alfdp=alfdp+dalf;
			betdp=betdp+dbet;
			deep_arg->xnode=FMod2p(deep_arg->xnode);
			xls=deep_arg->xll+deep_arg->omgadf+cosis*deep_arg->xnode;
			dls=m->pl+pgh-m->pinc*deep_arg->xnode*sinis;
			xls=xls+dls;
			xnoh=deep_arg->xnode;
			deep_arg->xnode=AcTan(alfdp,betdp);

			/* This is a patch to Lyddane modification */
			/* suggested by Rob Matson. */
		
			if (fabs(xnoh-deep_arg->xnode)>pi)
			{
			      if (deep_arg->xnode<xnoh)
				  deep_arg->xnode+=twopi;
			      else
				  deep_arg->xnode-=twopi;
			}

			deep_arg->xll=deep_arg->xll+m->pl;
			deep_arg->omgadf=xls-deep_arg->xll-cos(deep_arg->xinc)*deep_arg->xnode;
		}
		return;
	}
}

