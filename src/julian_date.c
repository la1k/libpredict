#include <predict/predict.h>
#include <stdio.h>

//helper function for getting the julian day start date as time_t (a bit difficult to define this as a const...)
time_t get_julian_start_day()
{
	struct tm start_time;
	start_time.tm_sec = 0;
	start_time.tm_min = 0;
	start_time.tm_hour = 0;
	start_time.tm_mday = 31;
	start_time.tm_mon = 11;
	start_time.tm_year = 1979-1900;
	time_t ret_time = timegm(&start_time);
	return ret_time;
}

#define NUM_SECONDS_IN_DAY (24.0*60.0*60.0)
predict_julian_date_t predict_get_julian_date_from_time(time_t input_time)
{
	//get number of seconds since 1979-12-31 00:00:00 UTC, convert to days
	double seconds = difftime(input_time, get_julian_start_day());
	return seconds/NUM_SECONDS_IN_DAY;
}

time_t predict_get_time_from_julian_date(predict_julian_date_t date)
{
	double seconds_since = date*NUM_SECONDS_IN_DAY;
	time_t ret_time = get_julian_start_day();
	
	//add number of seconds since julian start day to the julian start day, get current time_t
	struct tm timeinfo;
	gmtime_r(&ret_time, &timeinfo); 
	timeinfo.tm_sec += seconds_since;
	ret_time = timegm(&timeinfo);
	return ret_time;
}
