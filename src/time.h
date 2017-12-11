#ifndef _PREDICT_TIME_H_
#define _PREDICT_TIME_H_

//TODO: implement strong typedefs for different times based on double

// Julian days since 1979-12-31 00:00:00 UTC (historical reasons - old predict)
typedef double predict_time_t;

typedef double predict_julian_day_t;
typedef double predict_julian_ephimeris_day_t;

// Julian millenia since J2000.0
typedef double predict_julian_centuries_t;

// Julian millenia since J2000.0
typedef double predict_julian_millenia_t;

predict_julian_day_t predict_time_to_julian_day(predict_time_t t);

predict_julian_ephimeris_day_t predict_julian_day_to_ephimeris_day(predict_julian_day_t JD);

predict_julian_centuries_t predict_julian_day_to_centuries(predict_julian_day_t JD);

predict_julian_centuries_t predict_julian_ephimeris_day_to_centuries(predict_julian_ephimeris_day_t JDE);

#endif
