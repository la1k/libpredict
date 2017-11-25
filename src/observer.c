#include <predict/predict.h>
#include "unsorted.h"
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "sun.h"

void observer_calculate(const predict_observer_t *observer, double time, const double pos[3], const double vel[3], struct predict_observation *result);

predict_observer_t *predict_create_observer(const char *name, double lat, double lon, double alt)
{
	// Allocate memory
	predict_observer_t *obs = (predict_observer_t*)malloc(sizeof(predict_observer_t));
	if (obs == NULL) return NULL;

	strncpy(obs->name, name, 128);
	obs->name[127] = '\0';
	obs->latitude = lat;
	obs->longitude = lon;
	obs->altitude = alt;

	return obs;
}

void predict_destroy_observer(predict_observer_t *obs)
{
	if (obs != NULL) {
		free(obs);
	}
}


#define VISIBILITY_SUN_ELE_UPPER_THRESH -12.0
#define VISIBILITY_ORBIT_ELE_LOWER_THRESH 0.0
/**
 * \brief Calculates range, azimuth, elevation and relative velocity.
 *
 * Calculated range, azimuth, elevation and relative velocity from the
 * given observer position.
 **/
void predict_observe_orbit(const predict_observer_t *observer, const struct predict_orbit *orbit, struct predict_observation *obs)
{
	if (obs == NULL) return;
	
	double julTime = orbit->time + 2444238.5;

	observer_calculate(observer, julTime, orbit->position, orbit->velocity, obs);

	// Calculate visibility status of the orbit: Orbit is visible if sun elevation is low enough and the orbit is above the horizon, but still in sunlight.
	obs->visible = false;
	struct predict_observation sun_obs;
	predict_observe_sun(observer, orbit->time, &sun_obs);
	if (!(orbit->eclipsed) && (sun_obs.elevation*180.0/M_PI < VISIBILITY_SUN_ELE_UPPER_THRESH) && (obs->elevation*180.0/M_PI > VISIBILITY_ORBIT_ELE_LOWER_THRESH)) {
		obs->visible = true;
	}
	obs->time = orbit->time;
}

void observer_calculate(const predict_observer_t *observer, double time, const double pos[3], const double vel[3], struct predict_observation *result)
{
	
		/* The procedures Calculate_Obs and Calculate_RADec calculate         */
	/* the *topocentric* coordinates of the object with ECI position,     */
	/* {pos}, and velocity, {vel}, from location {geodetic} at {time}.    */
	/* The {obs_set} returned for Calculate_Obs consists of azimuth,      */
	/* elevation, range, and range rate (in that order) with units of     */
	/* radians, radians, kilometers, and kilometers/second, respectively. */
	/* The WGS '72 geoid is used and the effect of atmospheric refraction */
	/* (under standard temperature and pressure) is incorporated into the */
	/* elevation calculation; the effect of atmospheric refraction on     */
	/* range and range rate has not yet been quantified.                  */

	/* The {obs_set} for Calculate_RADec consists of right ascension and  */
	/* declination (in that order) in radians.  Again, calculations are   */
	/* based on *topocentric* position using the WGS '72 geoid and        */
	/* incorporating atmospheric refraction.                              */


	double obs_pos[3];
	double obs_vel[3];
	double range[3];
	double rgvel[3];
	
	geodetic_t geodetic;
	geodetic.lat = observer->latitude;
	geodetic.lon = observer->longitude;
	geodetic.alt = observer->altitude / 1000.0;
	geodetic.theta = 0.0;
	Calculate_User_PosVel(time, &geodetic, obs_pos, obs_vel);

	vec3_sub(pos, obs_pos, range);
	vec3_sub(vel, obs_vel, rgvel);
	
	double range_length = vec3_length(range);
	double range_rate_length = vec3_dot(range, rgvel) / range_length;

	double theta_dot = 2*M_PI*omega_E/secday;
	double sin_lat = sin(geodetic.lat);
	double cos_lat = cos(geodetic.lat);
	double sin_theta = sin(geodetic.theta);
	double cos_theta = cos(geodetic.theta);
	
	double top_s = sin_lat*cos_theta*range[0] + sin_lat*sin_theta*range[1] - cos_lat*range[2];
	double top_e = -sin_theta*range[0] + cos_theta*range[1];
	double top_z = cos_lat*cos_theta*range[0] + cos_lat*sin_theta*range[1] + sin_lat*range[2];


	double top_s_dot = sin_lat*(cos_theta*rgvel[0] - sin_theta*range[0]*theta_dot) + 
						sin_lat*(sin_theta*rgvel[1] + cos_theta*range[1]*theta_dot) -
						cos_lat*rgvel[2];
	double top_e_dot = - (sin_theta*rgvel[0] + cos_theta*range[0]*theta_dot) + 
						(cos_theta*rgvel[1] - sin_theta*range[1]*theta_dot);

	double top_z_dot = cos_lat * ( cos_theta*(rgvel[0] + range[1]*theta_dot) + 
								sin_theta*(rgvel[1] - range[0]*theta_dot) ) +
								sin_lat*rgvel[2];
	
	// Azimut
	double y = -top_e / top_s;
	double az = atan(-top_e / top_s);

	if (top_s > 0.0) az = az + M_PI;
	if (az < 0.0) az = az + 2*M_PI;

	// Azimut rate
	double y_dot = - (top_e_dot*top_s - top_s_dot*top_e) / (top_s*top_s);
	double az_dot = y_dot / (1 + y*y);

	// Elevation
	double x = top_z / range_length;
	double el = ArcSin(x);

	// Elevation rate
	double x_dot = (top_z_dot*range_length - range_rate_length*top_z) / (range_length * range_length);
	double el_dot = x_dot / sqrt( 1 - x*x );
	
	result->azimuth = az;
	result->azimuth_rate = az_dot;
	result->elevation = el;
	result->elevation_rate = el_dot;
	result->range = range_length;
	result->range_rate = range_rate_length; 
	result->range_x = range[0];
	result->range_y = range[1];
	result->range_z = range[2];

}

void predict_observe_sun(const predict_observer_t *observer, double time, struct predict_observation *obs)
{
	
	// Find sun position
	double solar_vector[3];
	sun_predict(time, solar_vector);

	/* Zero vector for initializations */
	double zero_vector[3] = {0,0,0};

	/* Solar observed azimuth and elevation vector  */
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
	
	double jul_utc = time + 2444238.5;
	Calculate_Obs(jul_utc, solar_vector, zero_vector, &geodetic, &solar_set);
	
	double sun_azi = solar_set.x; 
	double sun_ele = solar_set.y;

	double sun_range = 1.0+((solar_set.z-AU)/AU);
	double sun_range_rate = 1000.0*solar_set.w;

	Calculate_LatLonAlt(jul_utc, solar_vector, &solar_latlonalt);

	/*
	double sun_lat = Degrees(solar_latlonalt.lat);
	double sun_lon = 360.0-Degrees(solar_latlonalt.lon);
	*/

	Calculate_RADec(jul_utc, solar_vector, zero_vector, &geodetic, &solar_rad);

	/*
	double sun_ra = solar_rad.x ;
	double sun_dec = solar_rad.y;
	*/
	
	obs->time = time;
	obs->azimuth = sun_azi;
	obs->elevation = sun_ele;
	obs->range = sun_range;
	obs->range_rate = sun_range_rate;
}


/* This function determines the position of the moon, including the azimuth and elevation headings, relative to the latitude
	   and longitude of the tracking station.  This code was derived
	   from a Javascript implementation of the Meeus method for
	   determining the exact position of the Moon found at:
	   http://www.geocities.com/s_perona/ingles/poslun.htm. */
void predict_observe_moon(const predict_observer_t *observer, double time, struct predict_observation *obs)
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

	obs->azimuth = moon_az;
	obs->elevation = moon_el;

}

#define ELEVATION_ZERO_TOLERANCE 0.3 //threshold for fine-tuning of AOS/LOS
#define DAYNUM_MINUTE 1.0/(24*60) //number of days corresponding to a minute
struct predict_observation predict_next_aos(const predict_observer_t *observer, const predict_orbital_elements_t *orbital_elements, double start_utc)
{
	double curr_time = start_utc;
	struct predict_observation obs;
	double time_step = 0;

	struct predict_orbit orbit;
	predict_orbit(orbital_elements, &orbit, curr_time);
	predict_observe_orbit(observer, &orbit, &obs);

	//check whether AOS can happen after specified start time
	if (predict_aos_happens(orbital_elements, observer->latitude) && !predict_is_geosynchronous(orbital_elements) && !orbit.decayed) {
		//TODO: Time steps have been found in FindAOS/LOS().
		//Might be based on some pre-existing source, root-finding techniques
		//or something. Find them, and improve readability of the code and so that
		//the mathematical stability of the iteration can be checked.
		//Bisection method, Brent's algorithm? Given a coherent root finding algorithm,
		//can rather have one function for iterating the orbit and then let get_next_aos/los
		//specify bounding intervals for the root finding.

		//skip the rest of the pass if the satellite is currently in range, since we want the _next_ AOS.
		if (obs.elevation > 0.0) {
			struct predict_observation los = predict_next_los(observer, orbital_elements, curr_time);
			curr_time = los.time;
			curr_time += DAYNUM_MINUTE*20; //skip 20 minutes. LOS might still be within the elevation threshold. (rough quickfix from predict)
			predict_orbit(orbital_elements, &orbit, curr_time);
			predict_observe_orbit(observer, &orbit, &obs);
		}

		//iteration until the orbit is roughly in range again, before the satellite pass
		while ((obs.elevation*180.0/M_PI < -1.0) || (obs.elevation_rate < 0)) {
			time_step = 0.00035*(obs.elevation*180.0/M_PI*((orbit.altitude/8400.0)+0.46)-2.0);
			curr_time -= time_step;
			predict_orbit(orbital_elements, &orbit, curr_time);
			predict_observe_orbit(observer, &orbit, &obs);
		}

		//fine tune the results until the elevation is within a low enough threshold
		while (fabs(obs.elevation*180/M_PI) > ELEVATION_ZERO_TOLERANCE) {
			time_step = obs.elevation*180.0/M_PI*sqrt(orbit.altitude)/530000.0;
			curr_time -= time_step;
			predict_orbit(orbital_elements, &orbit, curr_time);
			predict_observe_orbit(observer, &orbit, &obs);
		}
	}
	return obs;
}

/**
 * Pass stepping direction used for pass stepping function below.
 **/
enum step_pass_direction{POSITIVE_DIRECTION, NEGATIVE_DIRECTION};

/**
 * Rough stepping through a pass. Uses weird time steps from Predict.
 *
 * \param observer Ground station
 * \param orbital_elements Orbital elements of satellite
 * \param curr_time Time from which to start stepping
 * \param direction Either POSITIVE_DIRECTION (step from current time to pass end) or NEGATIVE_DIRECTION (step from current time to start of pass). In case of the former, the pass will be stepped until either elevation is negative or the derivative of the elevation is negative
 * \return Time for when we have stepped out of the pass
 * \copyright GPLv2+
 **/
double step_pass(const predict_observer_t *observer, const predict_orbital_elements_t *orbital_elements, double curr_time, enum step_pass_direction direction) {
	struct predict_orbit orbit;
	struct predict_observation obs;
	do {
		predict_orbit(orbital_elements, &orbit, curr_time);
		predict_observe_orbit(observer, &orbit, &obs);

		//weird time stepping from Predict, but which magically works
		double time_step = cos(obs.elevation - 1.0)*sqrt(orbit.altitude)/25000.0;
		if (((direction == POSITIVE_DIRECTION) && time_step < 0) || ((direction == NEGATIVE_DIRECTION) && time_step > 0)) {
			time_step = -time_step;
		}

		curr_time += time_step;
	} while ((obs.elevation >= 0) || ((direction == POSITIVE_DIRECTION) && (obs.elevation_rate > 0.0)));
	return curr_time;
}

struct predict_observation predict_next_los(const predict_observer_t *observer, const predict_orbital_elements_t *orbital_elements, double start_utc)
{
	double curr_time = start_utc;
	struct predict_observation obs;
	double time_step = 0;

	struct predict_orbit orbit;
	predict_orbit(orbital_elements, &orbit, curr_time);
	predict_observe_orbit(observer, &orbit, &obs);

	//check whether AOS/LOS can happen after specified start time
	if (predict_aos_happens(orbital_elements, observer->latitude) && !predict_is_geosynchronous(orbital_elements) && !orbit.decayed) {
		//iteration algorithm from Predict, see comments in predict_next_aos().

		//iterate until next satellite pass
		if (obs.elevation < 0.0) {
			struct predict_observation aos = predict_next_aos(observer, orbital_elements, curr_time);
			curr_time = aos.time;
			predict_orbit(orbital_elements, &orbit, curr_time);
			predict_observe_orbit(observer, &orbit, &obs);
		}

		//step through the pass
		curr_time = step_pass(observer, orbital_elements, curr_time, POSITIVE_DIRECTION);

		//fine tune to elevation threshold
		do {
			time_step = obs.elevation*180.0/M_PI*sqrt(orbit.altitude)/502500.0;
			curr_time += time_step;
			predict_orbit(orbital_elements, &orbit, curr_time);
			predict_observe_orbit(observer, &orbit, &obs);
		} while (fabs(obs.elevation*180.0/M_PI) > ELEVATION_ZERO_TOLERANCE);
	}
	return obs;
}

/**
 * Convenience function for calculation of derivative of elevation at specific time.
 *
 * \param observer Ground station
 * \param orbital_elements Orbital elements for satellite
 * \param time Time
 * \return Derivative of elevation at input time
 **/
double elevation_derivative(const predict_observer_t *observer, const predict_orbital_elements_t *orbital_elements, double time)
{
	struct predict_orbit orbit;
	struct predict_observation observation;
	predict_orbit(orbital_elements, &orbit, time);
	predict_observe_orbit(observer, &orbit, &observation);
	return observation.elevation_rate;
}

#include <float.h>

//Threshold used for comparing lower and upper brackets in find_max_elevation
#define TIME_THRESHOLD FLT_EPSILON

//Maximum number of iterations in find_max_elevation
#define MAX_ITERATIONS 10000

/**
 * Find maximum elevation bracketed by input lower and upper time.
 *
 * \param observer Ground station
 * \param orbital_elements Orbital elements of satellite
 * \param lower_time Lower time bracket
 * \param upper_time Upper time bracket
 * \return Observation of satellite for maximum elevation between lower and upper time brackets
 **/
struct predict_observation find_max_elevation(const predict_observer_t *observer, const predict_orbital_elements_t *orbital_elements, double lower_time, double upper_time)
{
	double max_ele_time_candidate = (upper_time + lower_time)/2.0;
	int iteration = 0;
	while ((fabs(lower_time - upper_time) > TIME_THRESHOLD) && (iteration < MAX_ITERATIONS)) {
		max_ele_time_candidate = (upper_time + lower_time)/2.0;

		//calculate derivatives for lower, upper and candidate
		double candidate_deriv = elevation_derivative(observer, orbital_elements, max_ele_time_candidate);
		double lower_deriv = elevation_derivative(observer, orbital_elements, lower_time);
		double upper_deriv = elevation_derivative(observer, orbital_elements, upper_time);

		//check whether derivative has changed sign
		if (candidate_deriv*lower_deriv < 0) {
			upper_time = max_ele_time_candidate;
		} else if (candidate_deriv*upper_deriv < 0) {
			lower_time = max_ele_time_candidate;
		} else {
			break;
		}
		iteration++;
	}

	//prepare output
	struct predict_orbit orbit;
	struct predict_observation observation;
	predict_orbit(orbital_elements, &orbit, max_ele_time_candidate);
	predict_observe_orbit(observer, &orbit, &observation);
	return observation;
}

struct predict_observation predict_at_max_elevation(const predict_observer_t *observer, const predict_orbital_elements_t *orbital_elements, predict_julian_date_t start_time)
{
	struct predict_observation ret_observation = {0};

	if (predict_is_geosynchronous(orbital_elements)) {
		return ret_observation;
	}

	struct predict_orbit orbit;
	struct predict_observation observation;
	predict_orbit(orbital_elements, &orbit, start_time);
	if (orbit.decayed) {
		return ret_observation;
	}

	predict_observe_orbit(observer, &orbit, &observation);

	//bracket the solution by approximate times for AOS and LOS of the current or next pass
	double lower_time = start_time;
	double upper_time = start_time;
	if (observation.elevation < 0) {
		struct predict_observation aos = predict_next_aos(observer, orbital_elements, start_time);
		lower_time = aos.time;
	} else {
		lower_time = step_pass(observer, orbital_elements, lower_time, NEGATIVE_DIRECTION);
	}
	struct predict_observation los = predict_next_los(observer, orbital_elements, lower_time);
	upper_time = los.time;

	//assume that we can only have two potential local maxima along the elevation curve, and be content with that. For most cases, we will only have one, unless it is a satellite in deep space with long passes and weird effects.

	//bracket by AOS/LOS
	struct predict_observation candidate_center = find_max_elevation(observer, orbital_elements, lower_time, upper_time);

	//bracket by a combination of the found candidate above and either AOS or LOS (will thus cover solutions within [aos, candidate] and [candidate, los])
	struct predict_observation candidate_lower = find_max_elevation(observer, orbital_elements, candidate_center.time - TIME_THRESHOLD, upper_time);
	struct predict_observation candidate_upper = find_max_elevation(observer, orbital_elements, lower_time, candidate_center.time + TIME_THRESHOLD);

	//return the best candidate
	if ((candidate_center.elevation > candidate_lower.elevation) && (candidate_center.elevation > candidate_upper.elevation)) {
		return candidate_center;
	} else if (candidate_lower.elevation > candidate_upper.elevation) {
		return candidate_lower;
	} else {
		return candidate_upper;
	}
}

double predict_doppler_shift(const struct predict_observation *obs, double frequency)
{
	double sat_range_rate = obs->range_rate*1000.0; //convert to m/s
	return -frequency*sat_range_rate/SPEED_OF_LIGHT; //assumes that sat_range <<<<< speed of light, which is very ok
}
