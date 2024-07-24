#include <predict/predict.h>

#include <stdio.h>
#include <time.h>

#include "defs.h"

predict_julian_date_t predict_to_julian(time_t input_time)
{
	//get number of seconds since 1979-12-31 00:00:00 UTC, convert to days
	double seconds = difftime(input_time, JULIAN_START_DAY);
	return seconds/SECONDS_PER_DAY;
}

time_t predict_from_julian(predict_julian_date_t date)
{
	return date * SECONDS_PER_DAY + JULIAN_START_DAY;
}
