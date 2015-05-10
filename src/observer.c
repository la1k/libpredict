#include <predict/observer.h>
#include "unsorted.h"
#include <stdlib.h>
#include <string.h>
#include "defs.h"

observer_t *observer_create(const char *name, double lat, double lon, double alt)
{
	// Allocate memory
	observer_t *obs = (observer_t*)malloc(sizeof(observer_t));
	if (obs == NULL) return NULL;

	snprintf(obs->name, 128, name);
	obs->latitude = lat;
	obs->longitude = lon;
	obs->altitude = alt;

	return obs;
}

void observer_destroy(observer_t *obs)
{
	if (obs != NULL) {
		free(obs);
	}
}

/**
 * \brief Calculates range, azimut, elevation and relative velocity.
 *
 * Calculated range, azimut, elevation and relative velocity from the
 * given observer position.
 **/
void observer_find_orbit(const observer_t *observer, const orbit_t *orbit, struct observation *obs)
{
	if (obs == NULL) return;
	
	double julTime = orbit->time + 2444238.5;

	double sin_lat, cos_lat, sin_theta, cos_theta, azim, top_s, top_e, top_z;
	double obs_pos[3], obs_vel[3], range[3], rgvel[3];

	geodetic_t geodetic;
	geodetic.lat = observer->latitude;
	geodetic.lon = observer->longitude;
	geodetic.alt = observer->altitude / 1000.0;
	geodetic.theta = 0.0;

	Calculate_User_PosVel(julTime, &geodetic, obs_pos, obs_vel);

	//Find QTH -> SAT vector:
	vec3_sub(orbit->position, obs_pos, range);

	//Find QTH -> SAT relative velocity:
	vec3_sub(orbit->velocity, obs_vel, rgvel);

	//Direct distance: 
	double rangeLength = vec3_length(range);

	sin_lat = sin(geodetic.lat);
	cos_lat = cos(geodetic.lat);
	sin_theta = sin(geodetic.theta);
	cos_theta = cos(geodetic.theta);
	top_s = sin_lat*cos_theta*range[0]+sin_lat*sin_theta*range[1]-cos_lat*range[2];
	top_e = -sin_theta*range[0]+cos_theta*range[1];
	top_z = cos_lat*cos_theta*range[0]+cos_lat*sin_theta*range[1]+sin_lat*range[2];
	azim = atan(-top_e/top_s); /* Azimuth */

	if (top_s > 0.0) azim = azim + M_PI;
	if (azim < 0.0) azim = azim + 2*M_PI;

	obs->azimut = azim;
	obs->elevation = ArcSin(top_z/rangeLength);
	obs->range = rangeLength;
	obs->rangeDot = vec3_dot(range, rgvel) / rangeLength;

	/* Corrections for atmospheric refraction */
	/* Reference:  Astronomical Algorithms by Jean Meeus, pp. 101-104    */
	/* Correction is meaningless when apparent elevation is below horizon */
	obs->correctedElevation = obs->elevation+Radians((1.02/tan(Radians(Degrees(obs->elevation)+10.3/(Degrees(obs->elevation)+5.11))))/60);

	//Above horizon?
	obs->visible = (obs->elevation >= 0);

}

void observer_find_sun(const observer_t *observer, double time, struct observation *obs)
{
	
	double mjd, year, T, M, L, e, C, O, Lsa, nu, R, eps;

	double jul_utc = time+2444238.5;
	mjd=jul_utc-2415020.0;
	year=1900+mjd/365.25;
	T=(mjd+Delta_ET(year)/secday)/36525.0;
	M=Radians(Modulus(358.47583+Modulus(35999.04975*T,360.0)-(0.000150+0.0000033*T)*Sqr(T),360.0));
	L=Radians(Modulus(279.69668+Modulus(36000.76892*T,360.0)+0.0003025*Sqr(T),360.0));
	e=0.01675104-(0.0000418+0.000000126*T)*T;
	C=Radians((1.919460-(0.004789+0.000014*T)*T)*sin(M)+(0.020094-0.000100*T)*sin(2*M)+0.000293*sin(3*M));
	O=Radians(Modulus(259.18-1934.142*T,360.0));
	Lsa=Modulus(L+C-Radians(0.00569-0.00479*sin(O)), 2*M_PI);
	nu=Modulus(M+C, 2*M_PI);
	R=1.0000002*(1.0-Sqr(e))/(1.0+e*cos(nu));
	eps=Radians(23.452294-(0.0130125+(0.00000164-0.000000503*T)*T)*T+0.00256*cos(O));
	R=AU*R;

	double solar_vector[3];
	solar_vector[0] = R*cos(Lsa);
	solar_vector[1] = R*sin(Lsa)*cos(eps);
	solar_vector[2] = R*sin(Lsa)*sin(eps);

	
	/* Zero vector for initializations */
	double zero_vector[3] = {0,0,0};

	/* Solar observed azi and ele vector  */
	vector_t solar_set;

	/* Solar right ascension and declination vector */
	vector_t solar_rad;

	/* Solar lat, long, alt vector */
	geodetic_t solar_latlonalt;

	geodetic_t geodetic;
	geodetic.lat = observer->latitude;
	geodetic.lon = observer->longitude;
	geodetic.alt = observer->altitude / 1000.0;
	geodetic.theta = 0.0;
	Calculate_Obs(jul_utc, solar_vector, zero_vector, &geodetic, &solar_set);
	
	double sun_azi = solar_set.x; 
	double sun_ele = solar_set.y;

	double sun_range = 1.0+((solar_set.z-AU)/AU);
	double sun_range_rate = 1000.0*solar_set.w;

	Calculate_LatLonAlt(jul_utc, solar_vector, &solar_latlonalt);

	double sun_lat = Degrees(solar_latlonalt.lat);
	double sun_lon = 360.0-Degrees(solar_latlonalt.lon);

	Calculate_RADec(jul_utc, solar_vector, zero_vector, &geodetic, &solar_rad);

	double sun_ra = solar_rad.x ;
	double sun_dec = solar_rad.y;

	obs->time = time;
	obs->azimut = sun_azi;
	obs->elevation = sun_ele;
	obs->visible = (obs->elevation > 0);
	obs->range = sun_range;
	obs->rangeDot = sun_range_rate;
}


/* This function determines the position of the moon, including the azimuth and elevation headings, relative to the latitude
	   and longitude of the tracking station.  This code was derived
	   from a Javascript implementation of the Meeus method for
	   determining the exact position of the Moon found at:
	   http://www.geocities.com/s_perona/ingles/poslun.htm. */
void observer_find_moon(const observer_t *observer, double time, struct observation *obs)
{
	
	double	jd, ss, t, t1, t2, t3, d, ff, l1, m, m1, ex, om, l,
		b, w1, w2, bt, p, lm, h, ra, dec, z, ob, n, e, el,
		az, teg, th, mm, dv;

	jd = time + 2444238.5;

	t=(jd-2415020.0)/36525.0;
	t2=t*t;
	t3=t2*t;
	l1=270.434164+481267.8831*t-0.001133*t2+0.0000019*t3;
	m=358.475833+35999.0498*t-0.00015*t2-0.0000033*t3;
	m1=296.104608+477198.8491*t+0.009192*t2+0.0000144*t3;
	d=350.737486+445267.1142*t-0.001436*t2+0.0000019*t3;
	ff=11.250889+483202.0251*t-0.003211*t2-0.0000003*t3;
	om=259.183275-1934.142*t+0.002078*t2+0.0000022*t3;
	om = om * M_PI/180.0;
	
	/* Additive terms */

	l1=l1+0.000233*sin((51.2+20.2*t)*M_PI/180.0);
	ss=0.003964*sin((346.56+132.87*t-0.0091731*t2)*M_PI/180.0);
	l1=l1+ss+0.001964*sin(om);
	m=m-0.001778*sin((51.2+20.2*t)*M_PI/180.0);
	m1=m1+0.000817*sin((51.2+20.2*t)*M_PI/180.0);
	m1=m1+ss+0.002541*sin(om);
	d=d+0.002011*sin((51.2+20.2*t)*M_PI/180.0);
	d=d+ss+0.001964*sin(om);
	ff=ff+ss-0.024691*sin(om);
	ff=ff-0.004328*sin(om+(275.05-2.3*t)*M_PI/180.0);
	ex=1.0-0.002495*t-0.00000752*t2;
	om=om*M_PI/180.0;

	l1=PrimeAngle(l1);
	m=PrimeAngle(m);
	m1=PrimeAngle(m1);
	d=PrimeAngle(d);
	ff=PrimeAngle(ff);
	om=PrimeAngle(om);

	m=m*M_PI/180.0;
	m1=m1*M_PI/180.0;
	d=d*M_PI/180.0;
	ff=ff*M_PI/180.0;

	/* Ecliptic Longitude */

	l=l1+6.28875*sin(m1)+1.274018*sin(2.0*d-m1)+0.658309*sin(2.0*d);
	l=l+0.213616*sin(2.0*m1)-ex*0.185596*sin(m)-0.114336*sin(2.0*ff);
	l=l+0.058793*sin(2.0*d-2.0*m1)+ex*0.057212*sin(2.0*d-m-m1)+0.05332*sin(2.0*d+m1);
	l=l+ex*0.045874*sin(2.0*d-m)+ex*0.041024*sin(m1-m)-0.034718*sin(d);
	l=l-ex*0.030465*sin(m+m1)+0.015326*sin(2.0*d-2.0*ff)-0.012528*sin(2.0*ff+m1);
	
	l=l-0.01098*sin(2.0*ff-m1)+0.010674*sin(4.0*d-m1)+0.010034*sin(3.0*m1);
	l=l+0.008548*sin(4.0*d-2.0*m1)-ex*0.00791*sin(m-m1+2.0*d)-ex*0.006783*sin(2.0*d+m);
	
	l=l+0.005162*sin(m1-d)+ex*0.005*sin(m+d)+ex*0.004049*sin(m1-m+2.0*d);
	l=l+0.003996*sin(2.0*m1+2.0*d)+0.003862*sin(4.0*d)+0.003665*sin(2.0*d-3.0*m1);

	l=l+ex*0.002695*sin(2.0*m1-m)+0.002602*sin(m1-2.0*ff-2.0*d)+ex*0.002396*sin(2.0*d-m-2.0*m1);

	l=l-0.002349*sin(m1+d)+ex*ex*0.002249*sin(2.0*d-2.0*m)-ex*0.002125*sin(2.0*m1+m);

	l=l-ex*ex*0.002079*sin(2.0*m)+ex*ex*0.002059*sin(2.0*d-m1-2.0*m)-0.001773*sin(m1+2.0*d-2.0*ff);

	l=l+ex*0.00122*sin(4.0*d-m-m1)-0.00111*sin(2.0*m1+2.0*ff)+0.000892*sin(m1-3.0*d);

	l=l-ex*0.000811*sin(m+m1+2.0*d)+ex*0.000761*sin(4.0*d-m-2.0*m1)+ex*ex*.000717*sin(m1-2.0*m);

	l=l+ex*ex*0.000704*sin(m1-2.0*m-2.0*d)+ex*0.000693*sin(m-2.0*m1+2.0*d)+ex*0.000598*sin(2.0*d-m-2.0*ff)+0.00055*sin(m1+4.0*d);

	l=l+0.000538*sin(4.0*m1)+ex*0.000521*sin(4.0*d-m)+0.000486*sin(2.0*m1-d);

	l=l-0.001595*sin(2.0*ff+2.0*d);

	/* Ecliptic latitude */

	b=5.128189*sin(ff)+0.280606*sin(m1+ff)+0.277693*sin(m1-ff)+0.173238*sin(2.0*d-ff);
	b=b+0.055413*sin(2.0*d+ff-m1)+0.046272*sin(2.0*d-ff-m1)+0.032573*sin(2.0*d+ff);

	b=b+0.017198*sin(2.0*m1+ff)+9.266999e-03*sin(2.0*d+m1-ff)+0.008823*sin(2.0*m1-ff);
	b=b+ex*0.008247*sin(2.0*d-m-ff)+0.004323*sin(2.0*d-ff-2.0*m1)+0.0042*sin(2.0*d+ff+m1);

	b=b+ex*0.003372*sin(ff-m-2.0*d)+ex*0.002472*sin(2.0*d+ff-m-m1)+ex*0.002222*sin(2.0*d+ff-m);

	b=b+0.002072*sin(2.0*d-ff-m-m1)+ex*0.001877*sin(ff-m+m1)+0.001828*sin(4.0*d-ff-m1);

	b=b-ex*0.001803*sin(ff+m)-0.00175*sin(3.0*ff)+ex*0.00157*sin(m1-m-ff)-0.001487*sin(ff+d)-ex*0.001481*sin(ff+m+m1)+ex*0.001417*sin(ff-m-m1)+ex*0.00135*sin(ff-m)+0.00133*sin(ff-d);

	b=b+0.001106*sin(ff+3.0*m1)+0.00102*sin(4.0*d-ff)+0.000833*sin(ff+4.0*d-m1);

	b=b+0.000781*sin(m1-3.0*ff)+0.00067*sin(ff+4.0*d-2.0*m1)+0.000606*sin(2.0*d-3.0*ff);

	b=b+0.000597*sin(2.0*d+2.0*m1-ff)+ex*0.000492*sin(2.0*d+m1-m-ff)+0.00045*sin(2.0*m1-ff-2.0*d);

	b=b+0.000439*sin(3.0*m1-ff)+0.000423*sin(ff+2.0*d+2.0*m1)+0.000422*sin(2.0*d-ff-3.0*m1);

	b=b-ex*0.000367*sin(m+ff+2.0*d-m1)-ex*0.000353*sin(m+ff+2.0*d)+0.000331*sin(ff+4.0*d);

	b=b+ex*0.000317*sin(2.0*d+ff-m+m1)+ex*ex*0.000306*sin(2.0*d-2.0*m-ff)-0.000283*sin(m1+3.0*ff);

	w1=0.0004664*cos(om*M_PI/180.0);
	w2=0.0000754*cos((om+275.05-2.3*t)*M_PI/180.0);
	bt=b*(1.0-w1-w2);

	/* Parallax calculations */

	p=0.950724+0.051818*cos(m1)+0.009531*cos(2.0*d-m1)+0.007843*cos(2.0*d)+0.002824*cos(2.0*m1)+0.000857*cos(2.0*d+m1)+ex*0.000533*cos(2.0*d-m)+ex*0.000401*cos(2.0*d-m-m1);

	p=p+0.000173*cos(3.0*m1)+0.000167*cos(4.0*d-m1)-ex*0.000111*cos(m)+0.000103*cos(4.0*d-2.0*m1)-0.000084*cos(2.0*m1-2.0*d)-ex*0.000083*cos(2.0*d+m)+0.000079*cos(2.0*d+2.0*m1);

	p=p+0.000072*cos(4.0*d)+ex*0.000064*cos(2.0*d-m+m1)-ex*0.000063*cos(2.0*d+m-m1);

	p=p+ex*0.000041*cos(m+d)+ex*0.000035*cos(2.0*m1-m)-0.000033*cos(3.0*m1-2.0*d);

	p=p-0.00003*cos(m1+d)-0.000029*cos(2.0*ff-2.0*d)-ex*0.000029*cos(2.0*m1+m);

	p=p+ex*ex*0.000026*cos(2.0*d-2.0*m)-0.000023*cos(2.0*ff-2.0*d+m1)+ex*0.000019*cos(4.0*d-m-m1);

	b=bt*M_PI/180.0;
	lm=l*M_PI/180.0;
	double moon_dx=3.0/(M_PI*p);

	/* Semi-diameter calculation */
	/* sem=10800.0*asin(0.272488*p*M_PI/180.0)/pi; */
	/* Convert ecliptic coordinates to equatorial coordinates */

	z=(jd-2415020.5)/365.2422;
	ob=23.452294-(0.46845*z+5.9e-07*z*z)/3600.0;
	ob=ob*M_PI/180.0;
	dec=asin(sin(b)*cos(ob)+cos(b)*sin(ob)*sin(lm));
	ra=acos(cos(b)*cos(lm)/cos(dec));
	
	if (lm > M_PI)
		ra = 2*M_PI - ra;

	/* ra = right ascension */
	/* dec = declination */

	n = observer->latitude;    /* North latitude of tracking station */
	e = observer->longitude;  /* East longitude of tracking station */

	/* Find siderial time in radians */

	t=(jd-2451545.0)/36525.0;
	teg=280.46061837+360.98564736629*(jd-2451545.0)+(0.000387933*t-t*t/38710000.0)*t;

	while (teg>360.0)
		teg-=360.0;

	th = FixAngle(teg*M_PI/180.0 + e);
	h=th-ra;

	az=atan2(sin(h),cos(h)*sin(n)-tan(dec)*cos(n))+M_PI;
	el=asin(sin(n)*sin(dec)+cos(n)*cos(dec)*cos(h));

	double moon_az=az;
	double moon_el=el;

	/* Radial velocity approximation.  This code was derived
	   from "Amateur Radio Software", by John Morris, GM4ANB,
	   published by the RSGB in 1985. */

	mm=FixAngle(1.319238+time*0.228027135);  /* mean moon position */
	t2=0.10976;
	t1=mm+t2*sin(mm);
	dv=0.01255*moon_dx*moon_dx*sin(t1)*(1.0+t2*cos(mm));
	dv=dv*4449.0;
	t1=6378.0;
	t2=384401.0;
	t3=t1*t2*(cos(dec)*cos(n)*sin(h));
	t3=t3/sqrt(t2*t2-t2*t1*sin(el));

	//double moon_dv=dv+t3*0.0753125;
	//double moon_dec=dec/M_PI/180.0;
	double moon_ra=ra/M_PI/180.0;
	double moon_gha=teg-moon_ra;

	if (moon_gha<0.0) moon_gha+=360.0;

	obs->azimut = moon_az;
	obs->elevation = moon_el;

}


#define ELEVATION_ZERO_TOLERANCE 0.3 //threshold for fine-tuning of AOS/LOS
#define DAYNUM_MINUTE 1.0/(24*60) //number of days corresponding to a minute
double observer_get_next_aos(const observer_t *observer, orbit_t *orbit, double start_utc)
{
	double ret_aos_time = 0;
	double curr_time = start_utc;
	struct observation obs;
	double time_step = 0;
	
	orbit_predict(orbit, curr_time);
	observer_find_orbit(observer, orbit, &obs);

	//check whether AOS can happen after specified start time
	if (orbit_aos_happens(orbit, observer->latitude) && !orbit_is_geostationary(orbit) && !orbit_decayed(orbit, curr_time))
	{
		//TODO: Time steps have been found in FindAOS/LOS(). 
		//Might be based on some pre-existing source, root-finding techniques
		//or something. Find them, and improve readability of the code and so that
		//the mathematical stability of the iteration can be checked. 
		//Bisection method, Brent's algorithm? Given a coherent root finding algorithm, 
		//can rather have one function for iterating the orbit and then let get_next_aos/los 
		//specify bounding intervals for the root finding. 

		//skip the rest of the pass if the satellite is currently in range, since we want the _next_ AOS. 
		if (obs.elevation > 0.0)
		{
			curr_time = observer_get_next_los(observer, orbit, curr_time);
			curr_time += DAYNUM_MINUTE*20; //skip 20 minutes. LOS might still be within the elevation threshold. (rough quickfix from predict) 
			orbit_predict(orbit, curr_time);
			observer_find_orbit(observer, orbit, &obs);
		}

		//iteration until the orbit is roughly in range again, before the satellite pass
		while (obs.elevation*180.0/M_PI < -1.0)
		{
			time_step = 0.00035*(obs.elevation*180.0/M_PI*((orbit->altitude/8400.0)+0.46)-2.0); 
			curr_time -= time_step;
			orbit_predict(orbit, curr_time);
			observer_find_orbit(observer, orbit, &obs);
		}

		//fine tune the results until the elevation is within a low enough threshold
		while (fabs(obs.elevation*180/M_PI) > ELEVATION_ZERO_TOLERANCE)
		{
			time_step = obs.elevation*180.0/M_PI*sqrt(orbit->altitude)/530000.0;
			curr_time -= time_step;
			orbit_predict(orbit, curr_time);
			observer_find_orbit(observer, orbit, &obs);
		}

		ret_aos_time = curr_time;
	}
	return ret_aos_time;
}

double observer_get_next_los(const observer_t *observer, orbit_t *orbit, double start_utc)
{
	double ret_los_time = 0;
	double curr_time = start_utc;
	struct observation obs;
	double time_step = 0;

	orbit_predict(orbit, curr_time);
	observer_find_orbit(observer, orbit, &obs);

	//check whether AOS/LOS can happen after specified start time
	if (orbit_aos_happens(orbit, observer->latitude) && !orbit_is_geostationary(orbit) && !orbit_decayed(orbit, curr_time))
	{
		//iterate until next satellite pass
		if (obs.elevation < 0.0)
		{
			curr_time = observer_get_next_aos(observer, orbit, curr_time);
			orbit_predict(orbit, curr_time);
			observer_find_orbit(observer, orbit, &obs);
		}

		//step through the pass
		do 
		{
			time_step = cos(obs.elevation - 1.0)*sqrt(orbit->altitude)/25000.0; 
			curr_time += time_step;
			orbit_predict(orbit, curr_time);
			observer_find_orbit(observer, orbit, &obs);
		} 
		while (obs.elevation >= 0.0);
		
		//fine tune to elevation threshold
		do 
		{
			time_step = obs.elevation*180.0/M_PI*sqrt(orbit->altitude)/502500.0;
			curr_time += time_step;
			orbit_predict(orbit, curr_time);
			observer_find_orbit(observer, orbit, &obs);
		}
		while (fabs(obs.elevation*180.0/M_PI) > ELEVATION_ZERO_TOLERANCE);

		ret_los_time = curr_time;
	}
	return ret_los_time;

}

