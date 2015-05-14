#ifndef _JULIANDATE_H
#define _JULIANDATE_H

#include <time.h>

/* The representation of time used by libpredict: The number of days since 31Dec79 00:00:00 UTC. */
typedef double predict_julian_date_t;

/* Convert time_t to julian date in UTC. */
predict_julian_date_t predict_get_julian_date_from_time(time_t time);

/* Convert julian date in UTC back to a time_t. */
time_t predict_get_time_from_julian_date(predict_julian_date_t date);


#endif
