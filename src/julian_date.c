#define _POSIX_C_SOURCE 1
#include <predict/predict.h>
#include <predict/juliandate.h>
#include <stdio.h>


#define SECONDS_IN_HOUR (60*60)
/* assume timeinfo_utc describes time in UTC, return time_t value */
time_t mktime_utc(const struct tm* timeinfo_utc)
{
	time_t curr_time = time(NULL);
	int timezone_diff = 0; //deviation of the current timezone from UTC in seconds

	//get UTC time, interpret resulting tm as a localtime
	struct tm timeinfo_gmt;
	gmtime_r(&curr_time, &timeinfo_gmt);
	time_t time_gmt = mktime(&timeinfo_gmt);

	//get localtime, interpret resulting tm as localtime
	struct tm timeinfo_local;
	localtime_r(&curr_time, &timeinfo_local);
	time_t time_local = mktime(&timeinfo_local);

	//find the time difference between the two interpretations
	timezone_diff += difftime(time_local, time_gmt);

	//add daylight saving time, if any
	if (timeinfo_local.tm_isdst)
	{
		timezone_diff += SECONDS_IN_HOUR;
	}

	//hack for preventing mktime from assuming localtime: subtract timezone difference from the input struct.
	struct tm ret_timeinfo;
	ret_timeinfo.tm_sec = timeinfo_utc->tm_sec - timezone_diff; 
	ret_timeinfo.tm_min = timeinfo_utc->tm_min;
	ret_timeinfo.tm_hour = timeinfo_utc->tm_hour;
	ret_timeinfo.tm_mday = timeinfo_utc->tm_mday;
	ret_timeinfo.tm_mon = timeinfo_utc->tm_mon;
	ret_timeinfo.tm_year = timeinfo_utc->tm_year;
	return mktime(&ret_timeinfo);
}

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
	return mktime_utc(&start_time);
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
	ret_time = mktime_utc(&timeinfo);
	return ret_time;
}
