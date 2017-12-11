#include <stdio.h>
#include <math.h>

#include "body.h"
#include "coordinates.h"
#include "time.h"

void predict_helio_to_geo(predict_pos_t *pos_body, predict_pos_t *pos_earth) {
  if (pos_body->origin != CS_SUN) {
    printf("WARNING: Converting non-heliocentric position!");
  }
  pos_body->coords[0] -= pos_earth->coords[0];
  pos_body->coords[1] -= pos_earth->coords[1];
  // pos_body->coords[2] -= pos_earth->coords[2];
  pos_body->origin = CS_EARTH;
}

void predict_cart_to_sph(predict_pos_t *position) {
  if (position->type != CS_CARTESIAN) {
    printf("WARNING: Converting non-cartesian position!");
  }
  // TODO: find out which spherical convention is used in Meeus
  double r = sqrt(pow(position->coords[0],2) + pow(position->coords[1],2) + pow(position->coords[2],2));
  double phi = atan2(position->coords[1], position->coords[0]);
  double theta = atan2(position->coords[2], sqrt(pow(position->coords[0],2) + pow(position->coords[1],2)));

  position->coords[0] = r;
  position->coords[1] = theta;
  position->coords[2] = phi;
  position->type = CS_SPHERICAL;
}

void predict_ecliptic_to_equatorial(predict_julian_ephimeris_day_t JDE, predict_pos_t *position) {
  if (position->plane != CS_ECLIPTIC) {
    printf("WARNING: Converting non-ecliptic position!");
  }
  double T = predict_julian_ephimeris_day_to_centuries(JDE);
  // TODO: should be defined as a constant somewhere
  //const double epsilon_j2000 = 23.4392911 / 180 * M_PI; // obliquity of the ecliptic
  double epsilon_0 = (23.43929 + 0.01300417 * T + 1.666667e-7 * pow(T,2) +  5.027778e-7 * pow(T,3))/180*M_PI; // mean obliquity of the ecliptic
  //const double epsilon_j2000 = 23.43928 / 180 * M_PI;

  double Omega = (125.04452 - 1934.136261*T + 0.0020708 * pow(T,2) + pow(T,3) / 450000)/180*M_PI;
  double L = (280.4665 + 36000.7698*T)/180*M_PI; // mean longitude Sun
  double L_prime = (218.3165 + 481267.8813 * T)/180*M_PI; // mean longitude Moon

  double delta_epsilon = (9.2 * cos(Omega) + 0.57 * cos(2*L) + 0.1 * cos(2*L_prime) - 0.09 * cos(2*Omega))/M_PI/648000; // Meeus ch 22

  double epsilon = epsilon_0 + delta_epsilon;

  double lambda = position->coords[2]; // ecliptic longitude
  double beta = position->coords[1]; // ecliptic latitude

  // src: Meeus (13.3)
  double alpha = atan2(sin(lambda)*cos(epsilon)-tan(beta)*sin(epsilon),cos(lambda));
  // src: Meeus (13.4)
  double delta = asin(sin(beta)*cos(epsilon)+cos(beta)*sin(epsilon)*sin(lambda));

  position->coords[2] = alpha; // right acension
  position->coords[1] = delta; // declination
  position->plane = CS_EQUATORIAL;
}

void predict_equatorial_to_horizontal(predict_julian_day_t JD, predict_pos_t *position, const predict_observer_t *observer) {
  double alpha = position->coords[2];
  double delta = position->coords[1];

  /* Find siderial time in radians */
  // src: Meeus (11.1)
  double T=predict_julian_day_to_centuries(JD);
  // Greenwich Mean Siderial Time - Meeus (11.4)
  double theta_0 = (280.46061837+360.98564736629*(JD-2451545.0)+0.000387933*pow(T,2)-pow(T,3)/38710000.0)/180.0*M_PI;

  // Nutation corrections
  T = predict_julian_ephimeris_day_to_centuries(predict_julian_day_to_ephimeris_day(JD));
  //const double epsilon_j2000 = 23.4392911 / 180 * M_PI; // obliquity of the ecliptic
  double epsilon_0 = (23.43929 + 0.01300417 * T + 1.666667e-7 * pow(T,2) +  5.027778e-7 * pow(T,3))/180*M_PI; // mean obliquity of the ecliptic
  double Omega = (125.04452 - 1934.136261*T + 0.0020708 * pow(T,2) + pow(T,3) / 450000)/180*M_PI;
  double L = (280.4665 + 36000.7698*T)/180*M_PI; // mean longitude Sun
  double L_prime = (218.3165 + 481267.8813 * T)/180*M_PI; // mean longitude Moon

  double delta_psi = (-17.20*sin(Omega) - 1.32*sin(2*L) - 0.23*sin(2*L_prime) + 0.21 * sin(2*Omega))/M_PI/648000; // nutation in longitude
  double delta_epsilon = (9.2 * cos(Omega) + 0.57 * cos(2*L) + 0.1 * cos(2*L_prime) - 0.09 * cos(2*Omega))/M_PI/648000; // nutation in obliquity
  double epsilon = epsilon_0 + delta_epsilon;


  double theta = theta_0 + delta_psi*cos(epsilon); // Apparent Greenwich Siderial Time

  double H = theta + observer->longitude - alpha; // local hour angle (considering L positively eastwards)

  // Azimuth - Meeus (12.5)
  double A = atan2(sin(H),cos(H)*sin(observer->latitude)-tan(delta)*cos(observer->latitude))+M_PI;
  // Altitude - Meeus (12.6)
  double h = asin(sin(observer->latitude)*sin(delta)+cos(observer->latitude)*cos(delta)*cos(H));

  // TODO: correct distance of observation
  position->coords[2] = A;
  position->coords[1] = h;
}
