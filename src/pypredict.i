%module pypredict
%{
#include "predict/predict.h"
%}

/* typedefs */
typedef double predict_julian_date_t;

typedef struct {
int satellite_number;
long element_number;
char designator[10];
int epoch_year;
double epoch_day;
double inclination;
double right_ascension;
double eccentricity;
double argument_of_perigee;
double mean_anomaly;
double mean_motion;
double derivative_mean_motion;
double second_derivative_mean_motion;
double bstar_drag_term;
int revolutions_at_epoch;
} predict_orbital_elements_t;

typedef struct {
char name[128];
predict_julian_date_t time;
double position[3];
double velocity[3];
double latitude;
double longitude;
double altitude;
double footprint;
int eclipsed;
double eclipse_depth;
double phase;
long revolutions;
enum predict_ephemeris ephemeris;
predict_orbital_elements_t orbital_elements;
void *ephemeris_data;
} predict_orbit_t;

typedef struct {
char name[128];
double latitude;
double longitude;
double altitude;
} predict_observer_t;

struct predict_observation {
predict_julian_date_t time;
double azimuth;
double azimuth_rate;
double elevation;
double elevation_rate;
double range;
double range_x, range_y, range_z;
double range_rate;
bool visible;
};

typedef long time_t;

/* functions */
char *predict_version_string();
time_t time(time_t *time);
predict_julian_date_t predict_to_julian(time_t time);
time_t predict_from_julian(predict_julian_date_t date);
predict_orbit_t *predict_create_orbit(predict_orbital_elements_t orbital_elements);
void predict_destroy_orbit(predict_orbit_t *orbit);
int predict_orbit(predict_orbit_t *x, predict_julian_date_t time);
double predict_apogee(const predict_orbit_t *x);
double predict_perigee(const predict_orbit_t *x);
predict_observer_t *predict_create_observer(const char *name, double lat, double lon, double alt);
void predict_destroy_observer(predict_observer_t *obs);
void predict_observe_orbit(const predict_observer_t *observer, const predict_orbit_t *orbit, struct predict_observation *obs);
predict_julian_date_t predict_next_aos(const predict_observer_t *observer, predict_orbit_t *orbit, predict_julian_date_t start_time);
predict_julian_date_t predict_next_los(const predict_observer_t *observer, predict_orbit_t *orbit, predict_julian_date_t start_time);
double predict_doppler_shift(const predict_observer_t *observer, const predict_orbit_t *orbit, double downlink_frequency);
double predict_squint_angle(const predict_observer_t *observer, const predict_orbit_t *orbit, double alon, double alat);
