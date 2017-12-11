#include <math.h>

#include "time.h"
#include <stdio.h> // debug

predict_julian_day_t predict_time_to_julian_day(predict_time_t t) {
  return t + 2444238.5;
}

predict_julian_ephimeris_day_t predict_julian_day_to_ephimeris_day(predict_julian_day_t JD) {
  // src: https://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html
  // delta_t estimation for years 2005 to 2050
  double y = (JD - 2451545.0) / 365; // years since 2000
  double delta_t = 62.92 + 0.32217 * y + 0.005589 * pow(y,2); // in seconds
  //printf("delta t %f\n", delta_t);
  return JD + delta_t/(60*60*24);
}

predict_julian_centuries_t predict_julian_day_to_centuries(predict_julian_day_t JD) {
  return (JD - 2451545.0) / 36525;
}

predict_julian_centuries_t predict_julian_ephimeris_day_to_centuries(predict_julian_ephimeris_day_t JDE) {
  return (JDE - 2451545.0) / 36525;
}
