#ifndef _PREDICT_H_
#define _PREDICT_H_

/* The representation of time used by libpredict: The number of days since 31Dec79 00:00:00 UTC. */
typedef double predict_julian_date_t;

#include <predict/orbit.h>
#include <predict/observer.h>
#include <time.h>

/* Convert time_t to julian date. */
predict_julian_date_t predict_get_julian_date_from_time(time_t time);

/* Convert julian date back to a time_t. */
time_t predict_get_time_from_julian_date(predict_julian_date_t date);

#endif //_PREDICT_H_
