#include <math.h>
#include <stdio.h> // debug
#include "body.h"
#include "coordinates.h"
#include "time.h"

// Mean elements of the planets referred to the mean dynamical ecliptic and equinox of date
// Coefficients from Listings 5.9.*
// Simon et al. - Numerical expressions for precession formulae and mean elements for the Moon and the planets
// a, lambda, e, pi, i, Omega
// Memory Allocation for 7 coefficients per variable
// If definition is shorter, intrinsic zeros for the rest of the coefficients
const predict_body_keplarian_elements_t KEPLERIAN_ELEMENTS_VENUS = {
  {0.7233298200},
  {181.97980085, 2106691666.31989, 111.65021, 0.05368, -0.23516, -0.00179, 0.0020},
  {0.0067719164, -0.0004776521, 98127e-10, 4639e-10, 123e-10, -3e-10},
  {131.56370300, 50477.47081, -387.42545, -20.44048, -0.95948, 0.00044, 0.00020},
  {3.39466189, 36.13261, -0.31523, -0.02525, 0.00085, -0.00008},
  {76.67992019, 32437.57636, 146.22586, -0.33446, 0.23007, -0.00088, 0.00009}
};

const predict_body_keplarian_elements_t KEPLERIAN_ELEMENTS_EARTH = {
  {1.0000010178},
  {100.46645683, 1296027711.03429, 109.15809, 0.07207, -0.23530, -0.00180, 0.00020},
  {0.0167086342, -0.0004203654, -0.0000126734, 1444e-10, -2e-10, 3e-10},
  {102.93734808, 61900.55290, 164.47797, -0.06365, -0.12090, 0.00298, 0.00020},
  {0.},
  {0.}
};

const predict_body_keplarian_elements_t KEPLERIAN_ELEMENTS_JUPITER = {
  {5.2026032092, 19132e-10, -39e-10, -60e-10, -10e-10, 1e-10},
  {34.35151874, 109306899.89453, 80.38700, 0.13327, -0.18850, 0.00411, -0.00014},
  {0.0484979255, 0.0016322542, -0.0000471366, -20063e-10, 1018e-10, -21e-10, 1e-10},
  {14.33120687, 58054.86625, 370.95016, -16.07110, 0.51186, -0.02268, 0.00004},
  {1.30326698, -197.87442, 1.67744, -0.00838, -0.00735, 0.00085, 0.00004},
  {100.46440702, 36755.18747, 145.13295, 1.45556, -0.59609, -0.04324, 0.00175}
};

// TODO: add additional bodies if necessary


const predict_body_keplarian_elements_t *KEPLERIAN_ELEMENTS[10] = {
  NULL,
  NULL, //&KEPLERIAN_ELEMENTS_MERCURY,
  &KEPLERIAN_ELEMENTS_VENUS,
  &KEPLERIAN_ELEMENTS_EARTH,
  NULL, //&KEPLERIAN_ELEMENTS_MARS,
  &KEPLERIAN_ELEMENTS_JUPITER,
  NULL, //&KEPLERIAN_ELEMENTS_SATURN,
  NULL, //&KEPLERIAN_ELEMENTS_URANUS,
  NULL, //&KEPLERIAN_ELEMENTS_NEPTUNE,
  NULL, //&KEPLERIAN_ELEMENTS_PLUTO
};

void predict_body_calc_pos(predict_julian_ephimeris_day_t JDE, predict_body_t body, predict_pos_t *pos) {
  double t = predict_julian_ephimeris_day_to_centuries(JDE)/10.; // millenia past J2000.0
  // calculate Keplerian elements - constants
  double a =      KEPLERIAN_ELEMENTS[body]->a[0];
  double lambda = KEPLERIAN_ELEMENTS[body]->lambda[0] * 3600;
  double e =      KEPLERIAN_ELEMENTS[body]->e[0];
  double pi =     KEPLERIAN_ELEMENTS[body]->pi[0] * 3600;
  double i =      KEPLERIAN_ELEMENTS[body]->i[0] * 3600;
  double Omega =  KEPLERIAN_ELEMENTS[body]->Omega[0] * 3600;

  // calculate Keplerian elements - time series
  int n;
  double t_pow;
  for (n=1; n<N_COEFFS; n++) {
    if (n > 1) {
      t_pow = pow(t,n);
    } else {
      t_pow = t;
    }
    a       += KEPLERIAN_ELEMENTS[body]->a[n]       * t_pow;
    lambda  += KEPLERIAN_ELEMENTS[body]->lambda[n]  * t_pow; // time coefficients in arcsec
    e       += KEPLERIAN_ELEMENTS[body]->e[n]       * t_pow;
    pi      += KEPLERIAN_ELEMENTS[body]->pi[n]      * t_pow; // time coefficients in arcsec
    i       += KEPLERIAN_ELEMENTS[body]->i[n]       * t_pow; // time coefficients in arcsec
    Omega   += KEPLERIAN_ELEMENTS[body]->Omega[n]   * t_pow; // time coefficients in arcsec
  }

  // Keplerian elements from arcsec to radians
  lambda  = lambda / 3600 / 180*M_PI;
  pi      = pi     / 3600 / 180*M_PI;
  i       = i      / 3600 / 180*M_PI;
  Omega   = Omega  / 3600 / 180*M_PI;

  double omega = pi - Omega; // argument of perihelion
  double M = lambda - pi; // mean anomaly
  // for 3000 BC to 3000 AD + bT^2 + c * cos(f*T) + s * sin(f*T)

  // Solving Kepler's Equation
  // Meeus (30.7)
  double E = M;
  double delta_E = 1;
  double delta_M = 0;

  while (fabs(delta_E) > 1e-10) { // 1e-6 tolerance sufficient for mean ephimerides stuff
    delta_M = M + e * sin(E) - E;
    delta_E = delta_M / (1 - e * cos(E));
    E += delta_E;
  }

  double x_prime = a * (cos(E) - e);
  double y_prime = a * sqrt(1-pow(e,2)) * sin(E);

  double x_ecl = (cos(omega)*cos(Omega) - sin(omega)*sin(Omega)*cos(i)) * x_prime + (-sin(omega)*cos(Omega) - cos(omega)*sin(Omega)*cos(i))* y_prime;
  double y_ecl = (cos(omega)*sin(Omega) + sin(omega)*cos(Omega)* cos(i)) * x_prime + (-sin(omega)*sin(Omega) + cos(omega)*cos(Omega)*cos(i))* y_prime;
  double z_ecl = (sin(omega) * sin(i)) * x_prime + (cos(omega) * sin(i)) * y_prime;

  pos->coords[0] = x_ecl;
  pos->coords[1] = y_ecl;
  pos->coords[2] = z_ecl;
  pos->type = CS_CARTESIAN;
  pos->origin = CS_SUN;
  pos->plane = CS_ECLIPTIC;

}
