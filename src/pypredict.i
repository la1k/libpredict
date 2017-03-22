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

struct predict_orbit {
predict_julian_date_t time;
bool decayed;
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
double inclination;
double right_ascension;
double argument_of_perigee;
};

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
    enum predict_ephemeris ephemeris;
    void *ephemeris_data;
} predict_orbital_elements_t;

/* functions */
char *predict_version_string();
time_t time(time_t *time);
predict_julian_date_t predict_to_julian(time_t time);
time_t predict_from_julian(predict_julian_date_t date);
void predict_destroy_orbital_elements(predict_orbital_elements_t *orbital_elements);
int predict_orbit(const predict_orbital_elements_t *orbital_elements, struct predict_orbit *x, predict_julian_date_t time);
double predict_apogee(const predict_orbital_elements_t *x);
double predict_perigee(const predict_orbital_elements_t *x);
predict_orbital_elements_t* predict_parse_tle(const char *tle_line_1, const char *tle_line_2);
predict_observer_t *predict_create_observer(const char *name, double lat, double lon, double alt);
void predict_destroy_observer(predict_observer_t *obs);
void predict_observe_orbit(const predict_observer_t *observer, const struct predict_orbit *orbit, struct predict_observation *obs);
predict_julian_date_t predict_next_aos(const predict_observer_t *observer, const predict_orbital_elements_t *orbital_elements, predict_julian_date_t start_time);
predict_julian_date_t predict_next_los(const predict_observer_t *observer, const predict_orbital_elements_t *orbital_elements, predict_julian_date_t start_time);
double predict_doppler_shift(const predict_observer_t *observer, const struct predict_orbit *orbit, double downlink_frequency);
double predict_squint_angle(const predict_observer_t *observer, const struct predict_orbit *orbit, double alon, double alat);
