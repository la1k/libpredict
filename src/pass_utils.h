#ifndef AOSLOS_H_DEFINED
#define AOSLOS_H_DEFINED

#include <predict/predict.h>

/**
 * Pass stepping direction used for pass stepping function below.
 **/
enum step_pass_direction{POSITIVE_DIRECTION, NEGATIVE_DIRECTION};

/**
 * Rough stepping through a pass.
 *
 * \param observer Ground station
 * \param orbital_elements Orbital elements of satellite
 * \param curr_time Time from which to start stepping
 * \param direction Either POSITIVE_DIRECTION (step from current time to pass end) or NEGATIVE_DIRECTION (step from current time to start of pass). In case of the former, the pass will be stepped until either elevation is negative or the derivative of the elevation is negative
 * \return Time for when we have stepped out of the pass
 **/
predict_julian_date_t step_pass(const predict_observer_t *observer, const predict_orbital_elements_t *orbital_elements, predict_julian_date_t curr_time, enum step_pass_direction direction);

#endif
