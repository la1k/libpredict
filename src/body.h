#ifndef _BODY_H_
#define _BODY_H_

#include "coordinates.h"
#include "time.h"

typedef enum {
  SUN       = 0,
  MERCURY   = 1,
  VENUS     = 2,
  EARTH     = 3,
  MARS      = 4,
  JUPITER   = 5,
  SATURN    = 6,
  URANUS    = 7,
  NEPTUNE   = 8,
  PLUTO     = 9
} predict_body_t;

#define N_COEFFS 7

typedef struct {
  const double a[N_COEFFS];       // semi-major axis [au]
  const double lambda[N_COEFFS];  // mean longitude []
  const double e[N_COEFFS];       // eccentricity
  const double pi[N_COEFFS];      // longitude of the perihelion
  const double i[N_COEFFS];       // inclination
  const double Omega[N_COEFFS];   // longitude of the ascending node
} predict_body_keplarian_elements_t;

void predict_body_calc_pos(predict_julian_ephimeris_day_t JDE, predict_body_t body, predict_pos_t *pos);

#endif
