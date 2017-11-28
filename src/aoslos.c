#include "aoslos.h"

double step_pass(const predict_observer_t *observer, const predict_orbital_elements_t *orbital_elements, predict_julian_date_t curr_time, enum step_pass_direction direction)
{
	return 0;
}

struct predict_observation predict_next_aos(const predict_observer_t *observer, const predict_orbital_elements_t *orbital_elements, predict_julian_date_t start_utc)
{
	struct predict_observation ret_obs = {0};
	return ret_obs;
}


struct predict_observation predict_next_los(const predict_observer_t *observer, const predict_orbital_elements_t *orbital_elements, predict_julian_date_t start_utc)
{
	struct predict_observation ret_obs = {0};
	return ret_obs;
}
