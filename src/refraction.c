#include <predict/refraction.h>
#include <math.h>

/* Corrections for atmospheric refraction */
/* Reference:  Astronomical Algorithms by Jean Meeus, pp. 101-104    */
/* 			   http://en.wikipedia.org/wiki/Atmospheric_refraction */

#define A (1.02*M_PI/180.0)
#define B (10.3*M_PI/180.0)
#define C (5.11*M_PI/180.0)

/*!
 * \brief Calculate refraction angle.
 *
 * This function assumes atmospheric pressure of 101.0kPa and temperature 10deg celsius.
 *
 * \param el True elevation angle (rad).
 *
 * \return Refraction angle (rad).
 */
double predict_refraction(double el)
{
	return A / tan( el + B / ( el + C ) );
}

/*!
 * \brief Calculate refraction angle.
 *
 * Corrects for different atmospheric pressure and temperature.
 *
 * \param el True elevation angle in rads.
 * \param pressure Atmospheric pressure in kPa.
 * \param temp Temperature in deg celsius.
 *
 * \return Refraction angle (rad).
 */
double predict_refraction_ext(double el, double pressure, double temp)
{
	double x = 283*pressure / (101 * (273 + temp));
	return x * predict_refraction(el);
}

/*!
 * \brief Calculate refraction angle from apparent elevation.
 *
 * This function assumes atmospheric pressure of 101.0kPa and temperature 10deg celsius.
 *
 * \param apparent_el Apparent elevation angle (rad).
 *
 * \return Refraction angle (rad).
 */
double predict_refraction_from_apparent(double apparent_el)
{
	return 1.0 / tan( apparent_el + 7.31*M_PI/180.0 / (apparent_el + 4.4*M_PI/180.0));
}

/*!
 * \brief Calculate refraction angle from apparent elevation.
 *
 * Corrects for different atmospheric pressure and temperature.
 *
 * \param apparent_el Apparent elevation angle (rad).
 * \param pressure Atmospheric pressure in kPa.
 * \param temp Temperature in deg celsius.
 *
 * \return Refraction angle (rad).
 */
double predict_refraction_from_apparent_ext(double apparent_el, double pressure, double temp)
{
	double x = 283*pressure / (101 * (273 + temp));
	return x / tan( apparent_el + 7.31*M_PI/180.0 / (apparent_el + 4.4*M_PI/180.0));
}


/*!
 * \brief Calculate refraction rate of change.
 *
 * \param el True elevation angle (rad).
 * \param el_rate Rate of change of true elevation angle (rad/s).
 *
 * \return Refraction rate of change (rad/s).
 */
double predict_refraction_rate(double el, double el_rate)
{
	double u0 = el + C;
	double u1 = sin(el + B / (el + C));
	return A * el_rate * (B / (u0*u0) - 1.0) / (u1*u1);
}

/*!
 * \brief Calculate refraction rate of change.
 *
 * Corrects for different atmospheric pressure and temerature.
 *
 * \param el True elevation angle (rad).
 * \param el_rate Rate of change of true elevation angle (rad/s).
 * \param pressure Atmospheric pressure in kPa.
 * \param temp Temperature in deg celsius.
 *
 * \return Apparent elevation (rad).
 */
double predict_refraction_rate_ext(double el, double el_rate, double pressure, double temp)
{
	double x = 283*pressure / (101 * (273 + temp));
	return x * predict_refraction_rate(el, el_rate);
}

/*!
 * \brief Calculate apparent elevation from true elevation.
 *
 * \param el True elevation angle (rad).
 *
 * \return Apparent elevation (rad).
 */
double predict_apparent_elevation(double el)
{
	double apparent = el + predict_refraction(el);
	if (apparent >= 0.0) return apparent;
	else return el;
}

/*!
 * \brief Calculate apparent elevation from true elevation.
 *
 * Corrects for different atmospheric pressures and temperatures.
 *
 * \param el True elevation angle (rad).
 * \param pressure Atmospheric pressure (kPa).
 * \param temp Temperature (deg C).
 *
 * \return Apparent elevation (rad).
 */
double predict_apparent_elevation_ext(double el, double pressure, double temp)
{
	double apparent = el + predict_refraction_ext(el, pressure, temp);
	if (apparent >= 0.0) return apparent;
	else return el;
}

/*!
 * \brief Calculate apparent elevation rate.
 *
 * \param el True elevation angle (rad).
 * \param el_rate Rate of change of true elevation angle (rad/s).
 *
 * \return Rate of change of apparent elevation (rad/s).
 */
double predict_apparent_elevation_rate(double el, double el_rate)
{
	return el_rate * (1 + predict_refraction_rate(el, el_rate));
}

/*!
 * \brief Calculate apparent elevation rate.
 *
 * Corrects for different atmospheric pressures and temperatures.
 *
 * \param el True elevation angle (rad).
 * \param el_rate Rate of change of true elevation angle (rad/s).
 * \param pressure Atmospheric pressure (kPa).
 * \param temp Temperature (deg C).
 *
 * \return Rate of change of apparent elevation (rad/s).
 */
double predict_apparent_elevation_rate_ext(double el, double el_rate, double pressure, double temp)
{
	return el_rate * (1 + predict_refraction_rate_ext(el, el_rate, pressure, temp));
}

