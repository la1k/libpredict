#ifndef AOSLOS_H_DEFINED
#define AOSLOS_H_DEFINED

#include <predict/predict.h>

/**
 * Pass stepping direction used for pass stepping function below.
 **/
enum step_pass_direction{POSITIVE_DIRECTION, NEGATIVE_DIRECTION};

enum step_pass_type{STEP_UNTIL_PASS, STEP_OUT_OF_PASS};

predict_julian_date_t step_pass(const predict_observer_t *observer, const predict_orbital_elements_t *orbital_elements, predict_julian_date_t curr_time, enum step_pass_direction direction, enum step_pass_type type);

enum aos_los_search {
	AOS_SEARCH,
	LOS_SEARCH
};

struct predict_observation next_aos_los(const predict_observer_t *observer, const predict_orbital_elements_t *orbital_elements, predict_julian_date_t start_utc, enum aos_los_search search_type);

#endif
