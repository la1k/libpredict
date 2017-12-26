#include "pass_utils.h"
#include <math.h>
#include <float.h>
#include "defs.h"

void observe_orbit_at(const predict_observer_t *observer, const predict_orbital_elements_t *orbital_elements, predict_julian_date_t curr_time, struct predict_observation *observation)
{
	struct predict_position position;
	predict_orbit(orbital_elements, &position, curr_time);
	predict_observe_orbit(observer, &position, observation);
}

predict_julian_date_t step_pass(const predict_observer_t *observer, const predict_orbital_elements_t *orbital_elements, predict_julian_date_t curr_time, enum step_pass_direction direction, enum step_pass_type type)
{
	struct predict_observation observation;
	observe_orbit_at(observer, orbital_elements, curr_time, &observation);

	int elevation_sign = 1;
	if (type == STEP_UNTIL_PASS) {
		elevation_sign = 1;
	} else if (type == STEP_OUT_OF_PASS) {
		elevation_sign = -1;
	}

	if (observation.elevation*elevation_sign > 0) {
		return curr_time;
	}

	double mean_period = 1.0/(orbital_elements->mean_motion);

	//if the mean period is larger than the earth rotational period, then the periods in the elevation curve will be
	//governed by earth's rotation
	mean_period = fmin(mean_period, 1.0/EARTH_ROTATIONS_PER_SIDERIAL_DAY);

	double step = fabs(mean_period/60.0);
	if (direction == NEGATIVE_DIRECTION) {
		step = -step;
	}

	bool found_pass = false;

	struct predict_observation obs;
//	struct predict_observation prev_obs = observation;
//	predict_julian_date_t prev_time = curr_time;

	while (!found_pass) {
		curr_time += step;
		observe_orbit_at(observer, orbital_elements, curr_time, &obs);
		observation = obs;

/*		if (prev_obs.elevation_rate * obs.elevation_rate < 0) {
			//there is a local minima or maxima between obs_1 and obs_2
			//find approximate location by extrapolating from the two points using straight lines
			double b_1 = prev_obs.elevation - prev_obs.elevation_rate*prev_time*SECONDS_PER_DAY;
			double b_2 = obs.elevation - obs.elevation_rate*curr_time*SECONDS_PER_DAY;

			//candidate at approximate location of extremum
			predict_julian_date_t candidate_time = (b_2 - b_1)/(prev_obs.elevation_rate - obs.elevation_rate)/SECONDS_PER_DAY;
			observe_orbit_at(observer, orbital_elements, candidate_time, &observation);

			//adjust current time point with respect to position of the extremum
			curr_time = candidate_time + step;
			observe_orbit_at(observer, orbital_elements, curr_time, &obs);
		}
*/
		if (observation.elevation*elevation_sign > 0) {
			found_pass = true;
		}
//		prev_obs = obs;
//		prev_time = curr_time;
	}
	return observation.time;
}

bool aoslos_is_feasible(const predict_observer_t *observer, const predict_orbital_elements_t *orbital_elements, predict_julian_date_t start_utc)
{
	struct predict_position orbit;
	predict_orbit(orbital_elements, &orbit, start_utc);
	return !predict_is_geosynchronous(orbital_elements) && predict_aos_happens(orbital_elements, observer->latitude) && !orbit.decayed;
}

struct predict_observation next_aos_los(const predict_observer_t *observer, const predict_orbital_elements_t *orbital_elements, predict_julian_date_t start_utc, enum aos_los_search search_type)
{
	struct predict_observation observation = {0};
	if (!aoslos_is_feasible(observer, orbital_elements, start_utc)) {
		return observation;
	}

	observe_orbit_at(observer, orbital_elements, start_utc, &observation);

	//find lower and upper brackets for solution
	predict_julian_date_t lower_bracket;
	predict_julian_date_t upper_bracket;

	if (search_type == AOS_SEARCH) {
		if (observation.elevation > 0) {
			lower_bracket = step_pass(observer, orbital_elements, start_utc, POSITIVE_DIRECTION, STEP_OUT_OF_PASS);
		} else {
			lower_bracket = start_utc;
		}
		upper_bracket = step_pass(observer, orbital_elements, lower_bracket, POSITIVE_DIRECTION, STEP_UNTIL_PASS);
	} else {
		if (observation.elevation > 0) {
			lower_bracket = start_utc;
		} else {
			lower_bracket = step_pass(observer, orbital_elements, start_utc, POSITIVE_DIRECTION, STEP_UNTIL_PASS);
		}
		upper_bracket = step_pass(observer, orbital_elements, lower_bracket, POSITIVE_DIRECTION, STEP_OUT_OF_PASS);
	}

	const double AOSLOS_TIME_EQUALITY_THRESHOLD = FLT_EPSILON;
	const int AOSLOS_MAX_NUM_ITERATIONS = 1000;
	int iteration = 0;

	predict_julian_date_t time_candidate;
	struct predict_observation candidate, lower, upper;

	while ((fabs(lower_bracket - upper_bracket) > AOSLOS_TIME_EQUALITY_THRESHOLD) && (iteration < AOSLOS_MAX_NUM_ITERATIONS)) {
		time_candidate = (upper_bracket + lower_bracket)/2.0;

		observe_orbit_at(observer, orbital_elements, time_candidate, &candidate);
		observe_orbit_at(observer, orbital_elements, lower_bracket, &lower);
		observe_orbit_at(observer, orbital_elements, upper_bracket, &upper);

		if (candidate.elevation*lower.elevation < 0) {
			upper_bracket = time_candidate;
		} else if (candidate.elevation*upper.elevation < 0) {
			lower_bracket = time_candidate;
		} else {
			break;
		}
		iteration++;
	}

	observe_orbit_at(observer, orbital_elements, (upper_bracket + lower_bracket)/2.0, &candidate);
	return candidate;
}

predict_julian_date_t bisection_method(predict_orbital_elements_t *orbital_elements, predict_observer_t *observer, predict_julian_date_t lower_bracket, predict_julian_date_t upper_bracket)
{
	const double AOSLOS_TIME_EQUALITY_THRESHOLD = FLT_EPSILON;
	const int AOSLOS_MAX_NUM_ITERATIONS = 1000;
	struct predict_observation lower, upper, candidate;
	predict_julian_date_t time_candidate;
	int iteration = 0;
	while ((fabs(lower_bracket - upper_bracket) > AOSLOS_TIME_EQUALITY_THRESHOLD) && (iteration < AOSLOS_MAX_NUM_ITERATIONS)) {
		time_candidate = (upper_bracket + lower_bracket)/2.0;

		observe_orbit_at(observer, orbital_elements, time_candidate, &candidate);
		observe_orbit_at(observer, orbital_elements, lower_bracket, &lower);
		observe_orbit_at(observer, orbital_elements, upper_bracket, &upper);

		if (candidate.elevation*lower.elevation < 0) {
			upper_bracket = time_candidate;
		} else if (candidate.elevation*upper.elevation < 0) {
			lower_bracket = time_candidate;
		} else {
			break;
		}
		iteration++;
	}
	return time_candidate;
}
