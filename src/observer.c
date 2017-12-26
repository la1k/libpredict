#include <predict/predict.h>
#include "unsorted.h"
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "sun.h"
#include "pass_utils.h"

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


/**
 * \brief Calculates range, azimuth, elevation and relative velocity.
 *
 * Calculated range, azimuth, elevation and relative velocity from the
 * given observer position.
 **/
void predict_observe_orbit(const predict_observer_t *observer, const struct predict_position *orbit, struct predict_observation *obs)
{
	if (obs == NULL) return;
	
	double julTime = orbit->time + JULIAN_TIME_DIFF;

	observer_calculate(observer, julTime, orbit->position, orbit->velocity, obs);

	// Calculate visibility status of the orbit: Orbit is visible if sun elevation is low enough and the orbit is above the horizon, but still in sunlight.
	obs->visible = false;
	struct predict_observation sun_obs;
	predict_observe_sun(observer, orbit->time, &sun_obs);
	if (!(orbit->eclipsed) && (sun_obs.elevation*180.0/M_PI < NAUTICAL_TWILIGHT_SUN_ELEVATION) && (obs->elevation*180.0/M_PI > 0)) {
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

	double theta_dot = 2*M_PI*EARTH_ROTATIONS_PER_SIDERIAL_DAY/SECONDS_PER_DAY;
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
	double el = asin_(x);

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

struct predict_observation predict_at_max_elevation(const predict_observer_t *observer, const predict_orbital_elements_t *orbital_elements, predict_julian_date_t start_time)
{
	struct predict_observation ret_observation = {0};

	if (predict_is_geosynchronous(orbital_elements)) {
		return ret_observation;
	}

	struct predict_position orbit;
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
		lower_time = step_pass(observer, orbital_elements, lower_time, NEGATIVE_DIRECTION, STEP_OUT_OF_PASS);
	}
	struct predict_observation los = predict_next_los(observer, orbital_elements, lower_time);
	upper_time = los.time;

	//assume that we can only have two potential local maxima along the elevation curve, and be content with that. For most cases, we will only have one, unless it is a satellite in deep space with long passes and weird effects.

	//bracket by AOS/LOS
	struct predict_observation candidate_center = find_elevation_derivative_root(observer, orbital_elements, lower_time, upper_time);

	//bracket by a combination of the found candidate above and either AOS or LOS (will thus cover solutions within [aos, candidate] and [candidate, los])
	struct predict_observation candidate_lower = find_elevation_derivative_root(observer, orbital_elements, candidate_center.time - MAXELE_TIME_EQUALITY_THRESHOLD, upper_time);
	struct predict_observation candidate_upper = find_elevation_derivative_root(observer, orbital_elements, lower_time, candidate_center.time + MAXELE_TIME_EQUALITY_THRESHOLD);

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

struct predict_observation predict_next_aos(const predict_observer_t *observer, const predict_orbital_elements_t *orbital_elements, predict_julian_date_t start_utc)
{
	return next_aos_los(observer, orbital_elements, start_utc, AOS_SEARCH);
}


struct predict_observation predict_next_los(const predict_observer_t *observer, const predict_orbital_elements_t *orbital_elements, predict_julian_date_t start_utc)
{
	return next_aos_los(observer, orbital_elements, start_utc, LOS_SEARCH);
}
