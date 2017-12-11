#ifndef _PREDICT_COORDINATES_H_
#define _PREDICT_COORDINATES_H_

#include <predict/predict.h>
#include "time.h"

typedef enum {
  CS_CARTESIAN,
  CS_SPHERICAL
} predict_cs_type_t;

typedef enum {
  CS_OBSERVER,
  CS_EARTH,
  CS_SUN
} predict_cs_origin_t;

typedef enum {
  CS_HORIZONTAL,
  CS_EQUATORIAL,
  CS_ECLIPTIC,
  CS_GALACTIC
} predict_cs_fund_plane_t;

typedef struct {
  double coords[3];
  predict_cs_type_t type;
  predict_cs_origin_t origin;
  predict_cs_fund_plane_t plane;
} predict_pos_t;

void predict_helio_to_geo(predict_pos_t *pos_body, predict_pos_t *pos_earth);

void predict_cart_to_sph(predict_pos_t *position);

void predict_ecliptic_to_equatorial(predict_julian_ephimeris_day_t JDE, predict_pos_t *position);

void predict_equatorial_to_horizontal(predict_julian_day_t JD, predict_pos_t *position, const predict_observer_t *observer);

#endif
