#include "orbit.h"
#include "defs.h"
#include "vec3.h"

#define EPHEMERIS_SGP4	0
#define EPHEMERIS_SDP4	1
#define EPHEMERIS_SGP8	2
#define EPHEMERIS_SDP8	3

/**
 * \brief Prepares internal data for use in prediction
 * functions.
 *
 * This function is a combination of the original PreCalc() and
 * select_ephemeris() functions.
 **/
int orbit_init(struct orbit *m, const char *tle[])
{
	//Parse TLE
	double tempnum;
	m->catnum = atol(SubString(tle[0],2,6));
	m->setnum = atol(SubString(tle[0],64,67));
	m->year = atoi(SubString(tle[0],18,19));
	m->refepoch = atof(SubString(tle[0],20,31));
	m->incl = atof(SubString(tle[1],8,15));
	m->raan = atof(SubString(tle[1],17,24));
	m->eccn = 1.0e-07*atof(SubString(tle[1],26,32));
	m->argper = atof(SubString(tle[1],34,41));
	m->meanan = atof(SubString(tle[1],43,50));
	m->meanmo = atof(SubString(tle[1],52,62));
	m->drag  = atof(SubString(tle[0],33,42));
	tempnum=1.0e-5*atof(SubString(tle[0],44,49));
	m->nddot6 = tempnum/pow(10.0,(tle[0][51]-'0'));
	tempnum=1.0e-5*atof(SubString(tle[0],53,58));
	m->bstar = tempnum/pow(10.0,(tle[0][60]-'0'));
	m->orbitnum = atof(SubString(tle[1],63,67));
	
	/* Fill in tle structure: */
	double temp = twopi/xmnpda/xmnpda;
	m->tle.catnr = m->catnum;
	m->tle.epoch = 1000.0*m->year + m->refepoch;
	m->tle.xndt2o = m->drag*temp;
	m->tle.xndd6o = m->nddot6*temp/xmnpda;
	m->tle.bstar = m->bstar / ae;
	m->tle.xincl = m->incl * M_PI / 180.0;
	m->tle.xnodeo = m->raan * M_PI / 180.0;
	m->tle.eo = m->eccn;
	m->tle.omegao = m->argper * M_PI / 180.0;
	m->tle.xmo = m->meanan * M_PI / 180.0;
	m->tle.xno = m->meanmo*temp*xmnpda;
	m->tle.revnum = m->orbitnum;
	
	/* Period > 225 minutes is deep space */
	double ao, xnodp, dd1, dd2, delo, a1, del1, r1;
	dd1=(xke/m->tle.xno);
	dd2=tothrd;
	a1=pow(dd1,dd2);
	r1=cos(m->tle.xincl);
	dd1=(1.0-m->tle.eo*m->tle.eo);
	temp=ck2*1.5f*(r1*r1*3.0-1.0)/pow(dd1,1.5);
	del1=temp/(a1*a1);
	ao=a1*(1.0-del1*(tothrd*.5+del1*(del1*1.654320987654321+1.0)));
	delo=temp/(ao*ao);
	xnodp=m->tle.xno/(delo+1.0);

	/* Select a deep-space/near-earth ephemeris */
	if (twopi/xnodp/xmnpda >= 0.15625) {
		m->ephemeris = EPHEMERIS_SDP4;
		sdp4_init(&m->sdp4);
	}else {
		m->ephemeris = EPHEMERIS_SGP4;
		sgp4_init(&m->sgp4);
	}
	return 0;
}

/**
 * \brief Is this satellite considered geostationary?
 *
 * This function returns true if the satellite is geostationary.
 **/
int orbit_is_geostationary(const struct orbit *m)
{
	if (fabs(m->meanmo-1.0027)<0.0002) {
		return 1;
	}else {
		return 0;
	}
}

double orbit_apogee(const struct orbit *m)
{
	double sma = 331.25*exp(log(1440.0/m->meanmo)*(2.0/3.0));
	return sma*(1.0+m->eccn)-xkmper;
}
		
double orbit_perigee(const struct orbit *m)
{
	double xno = m->meanmo*twopi/xmnpda;
	double a1=pow(xke/xno,tothrd);
	double cosio=cos(m->incl*M_PI/180.0);
	double theta2=cosio*cosio;
	double x3thm1=3*theta2-1.0;
	double eosq=m->eccn*m->eccn;
	double betao2=1.0-eosq;
	double betao=sqrt(betao2);
	double del1=1.5*ck2*x3thm1/(a1*a1*betao*betao2);
	double ao=a1*(1.0-del1*(0.5*tothrd+del1*(1.0+134.0/81.0*del1)));
	double delo=1.5*ck2*x3thm1/(ao*ao*betao*betao2);
	double aodp=ao/(1.0-delo);

	return (aodp*(1-m->eccn)-ae)*xkmper;
}

int orbit_aos_happens(const struct orbit *m, double latitude)
{
	/* This function returns a 1 if the satellite pointed to by
	   "x" can ever rise above the horizon of the ground station. */

	double lin, sma, apogee;

	if (m->meanmo==0.0)
		return false;
	else
	{
		lin = m->incl;

		if (lin >= 90.0) lin = 180.0-lin;

		sma = 331.25*exp(log(1440.0/m->meanmo)*(2.0/3.0));
		apogee = sma*(1.0+m->eccn)-xkmper;

		if ((acos(xkmper/(apogee+xkmper))+(lin*M_PI/180.0)) > fabs(latitude*M_PI/180.0))
			return true;
		else
			return false;
	}
}

/* This is the stuff we need to do repetitively while tracking. */
/* This is the old Calc() function. */
int orbit_predict(struct orbit *m, double utc)
{
	/* Set time to now if now time is provided: */
	if (utc == 0) utc = CurrentDaynum();
	
	/* Satellite position and velocity vectors */
	vec3_set(m->position, 0, 0, 0);
	vec3_set(m->velocity, 0, 0, 0);

	/* Solar ECI position vector  */
	double solar_vector[3];
	vec3_set(solar_vector, 0, 0, 0);

	/* Solar observed azi and ele vector  */
	//vector_t solar_set;

	m->time = utc;
	double julTime = utc + 2444238.5;

	/* Convert satellite's epoch time to Julian  */
	/* and calculate time since epoch in minutes */

	double jul_epoch = Julian_Date_of_Epoch(m->tle.epoch);
	double tsince = (julTime - jul_epoch)*xmnpda;
	//double age=m->time-jul_epoch;

	/* Call NORAD routines according to deep-space flag. */
	switch (m->ephemeris) {
		case EPHEMERIS_SDP4:
			sdp4_predict(&m->sdp4, tsince, &m->tle, m->position, m->velocity);
			break;
		case EPHEMERIS_SGP4:
			sgp4_predict(&m->sgp4, tsince, &m->tle, m->position, m->velocity);
			break;
		default:
			//TODO: Panic!
			break;
	}

	/* TODO: Remove? Scale position and velocity vectors to km and km/sec */
	Convert_Sat_State(m->position, m->velocity);

	/* Calculate satellite Lat North, Lon East and Alt. */
	geodetic_t sat_geodetic;
	Calculate_LatLonAlt(utc, m->position, &sat_geodetic);

	m->latitude = sat_geodetic.lat;
	m->longitude = sat_geodetic.lon;
	m->altitude = sat_geodetic.alt;

	/* TODO: ?? Calculate squint angle */
	//if (calc_squint) squint=(acos(-(ax*rx+ay*ry+az*rz)/obs_set.z))/deg2rad;

	/* Calculate solar position and satellite eclipse depth. */
	/* Also set or clear the satellite eclipsed flag accordingly. */
//	Calculate_Solar_Position(julTime, &solar_vector);
	//Calculate_Obs(julTime, &solar_vector, &zero_vector, &obs_geodetic, &solar_set);
//	double eclipseDepth;
//	m->eclipsed = Sat_Eclipsed(&pos, &solar_vector, &eclipseDepth);

/*
	fk=12756.33*acos(xkmper/(xkmper+sat_alt));
	fm=fk/1.609344;
	rv=(long)floor((tle.xno*xmnpda/twopi+age*tle.bstar*ae)*age+tle.xmo/twopi)+tle.revnum;
	sun_azi=Degrees(solar_set.x); 
	sun_ele=Degrees(solar_set.y);
	irk=(long)rint(sat_range);
	isplat=(int)rint(sat_lat);
	isplong=(int)rint(360.0-sat_lon);
	iaz=(int)rint(sat_azi);
	iel=(int)rint(sat_ele);
	ma256=(int)rint(256.0*(phase/twopi));

	if (!eclipsed())
	{
		if (sun_ele<=-12.0 && rint(sat_ele)>=0.0)
			findsun='+';
		else
			findsun='*';
	}
	else
		findsun=' ';
*/

	return 0;
}
