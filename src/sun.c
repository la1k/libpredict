#include "sun.h"
#include "unsorted.h"
#include "defs.h"

/**
 * The function Delta_ET has been added to allow calculations on the position of the sun.  It provides the difference between UT (approximately the same as UTC) and ET (now referred to as TDT). This function is based on a least squares fit of data from 1950 to 1991 and will need to be updated periodically. Values determined using data from 1950-1991 in the 1990 Astronomical Almanac.  See DELTA_ET.WQ1 for details.
 *
 * \copyright GPLv2+
 **/
double Delta_ET(double year)
{
	double delta_et;

	delta_et=26.465+0.747622*(year-1950)+1.886913*sin(2*M_PI*(year-1975)/33);

	return delta_et;
}

/**
 * Returns angle in radians from argument in degrees.
 *
 * \copyright GPLv2+
 **/
double Radians(double arg)
{
	/* Returns angle in radians from argument in degrees */
	return (arg*M_PI/180.0);
}

/**
 * Returns angle in degrees from argument in radians.
 *
 * \copyright GPLv2+
 **/
double Degrees(double arg)
{
	/* Returns angle in degrees from argument in radians */
	return (arg*180.0/M_PI);
}

void sun_predict(double time, double position[3])
{
	double JD = time + 2444238.5; // julian date from predict time
	double T = (JD - 2451545.0) / 36525; // julian centuries - Meeus (25.1)
	double M = (357.52911 + 35999.05029*T-0.0001537*pow(T,2))/180*M_PI; // mean anomaly - Meeus (25.3)
	double L_0 = (280.46646 + 36000.76983*T+0.0003032*pow(T,2))/180*M_PI; // mean longitude - Meeus (25.2)
	double e = 0.016708634 - 0.000042037*T - 0.0000001267*pow(T,2); // eccentricyit of Earths orbit - Meeus (25.4)
	double C = ((1.914602 - 0.004817*T - 0.000014*pow(T,2)) * sin(M) + (0.019993 - 0.000101 * T) * sin(2*M) + 0.000289 * sin(3*M))/180*M_PI; // Sun's equation of center - Meeus p. 164
	double L = L_0 + C; // Sun's true longitude
	double nu = M+C;	// Sun's true anomaly
	double R = 1.000001018*(1-pow(e,2)) / (1 + e * cos(nu)); // Sun's radius vector - Meeus (25.5)
	// result here is the true geometric longitude referred to the mean equinox of the date
	R = ASTRONOMICAL_UNIT_KM*R; // R in AU to R in km
	const double eps = 23.43928 / 180 * M_PI; // obliquity of the ecliptic at J2000.0
	// TODO: check if we need to calc eps as mean obliquiy of the ecliptic and not just using the constant there

	// rectangular coordinates of the sun reffered to the mean equator and equinox of the date - Meeus (26.1)
	position[0] = R*cos(L);
	position[1] = R*sin(L)*cos(eps);
	position[2] = R*sin(L)*sin(eps);
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

	geodetic_t geodetic;
	geodetic.lat = observer->latitude;
	geodetic.lon = observer->longitude;
	geodetic.alt = observer->altitude / 1000.0;
	geodetic.theta = 0.0;

	double jul_utc = time + JULIAN_TIME_DIFF;
	Calculate_Obs(jul_utc, solar_vector, zero_vector, &geodetic, &solar_set);

	double sun_azi = solar_set.x;
	double sun_ele = solar_set.y;

	double sun_range = 1.0+((solar_set.z-ASTRONOMICAL_UNIT_KM)/ASTRONOMICAL_UNIT_KM);
	double sun_range_rate = 1000.0*solar_set.w;

	obs->time = time;
	obs->azimuth = sun_azi;
	obs->elevation = sun_ele;
	obs->range = sun_range;
	obs->range_rate = sun_range_rate;
}

/**
 * Calculate RA and dec for the sun.
 *
 * \param time Time
 * \param ra Right ascension
 * \param dec Declination
 * \copyright GPLv2+
 **/
void predict_sun_ra_dec(predict_julian_date_t time, double *ra, double *dec)
{
	//predict absolute position of the sun
	double solar_vector[3];
	sun_predict(time, solar_vector);

	//prepare for radec calculation
	double jul_utc = time + JULIAN_TIME_DIFF;
	double zero_vector[3] = {0,0,0};
	vector_t solar_rad;

	//for some reason, RADec requires QTH coordinates, though
	//the properties to be calculated are observer-independent.
	//Pick some coordinates, will be correct anyway.
	geodetic_t geodetic;
	geodetic.lat = 10;
	geodetic.lon = 10;
	geodetic.alt = 10 / 1000.0;
	geodetic.theta = 0.0;

	//calculate right ascension/declination
	Calculate_RADec(jul_utc, solar_vector, zero_vector, &geodetic, &solar_rad);
	*ra = solar_rad.x;
	*dec = solar_rad.y;
}

double predict_sun_ra(predict_julian_date_t time)
{
	double ra, dec;
	predict_sun_ra_dec(time, &ra, &dec);
	return ra;
}

double predict_sun_declination(predict_julian_date_t time)
{
	double ra, dec;
	predict_sun_ra_dec(time, &ra, &dec);
	return dec;
}

double predict_sun_gha(predict_julian_date_t time)
{
	//predict absolute position of sun
	double solar_vector[3];
	sun_predict(time, solar_vector);

	//convert to lat/lon/alt
	geodetic_t solar_latlonalt;
	Calculate_LatLonAlt(time, solar_vector, &solar_latlonalt);

	//return longitude as the GHA
	double sun_lon = 360.0-Degrees(solar_latlonalt.lon);
	return sun_lon*M_PI/180.0;
}
