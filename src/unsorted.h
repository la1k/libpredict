#ifndef UNSORTED_H_
#define UNSORTED_H_

#include <predict/predict.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>

/**
 * Set three-element vector to specified components.
 *
 * \param v Output vector
 * \param x x-component
 * \param y y-component
 * \param z z-component
 **/
void vec3_set(double v[3], double x, double y, double z);

/**
 * Get length of vector.
 *
 * \param v Input vector
 * \return Euclidean length
 **/
double vec3_length(const double v[3]);

/**
 * Multiply vector by scalar.
 *
 * \param v Input vector. Overwritten by result
 * \param a Scalar to be multiplied into vector
 * \param r Resulting vector
 **/
void vec3_mul_scalar(const double v[3], double a, double r[3]);

/**
 * Subtract a vector 2 from vector 1.
 *
 * \param v1 Vector 1
 * \param v2 Vector 2
 * \param r Resulting vector
 **/
void vec3_sub(const double v1[3], const double v2[3], double *r);

/**
 * Dot product between two vectors.
 *
 * \param v Vector
 * \param u Another vector
 * \return Dot product
 **/
double vec3_dot(const double v[3], const double u[3]);

/**
 * Geodetic position structure used by SGP4/SDP4 code.
 **/
typedef struct	{
		   double lat, lon, alt, theta;
		}  geodetic_t;

/**
 * General three-dimensional vector structure used by SGP4/SDP4 code.
 **/
typedef struct	{
		   double x, y, z, w;
		}  vector_t;

/**
 * This function returns a substring based on the starting
 * and ending positions provided. Trims whitespaces.
 *
 * \param input_string Full input string
 * \param buffer_length Length of output_buffer
 * \param output_buffer Returned substring
 * \param start Start position
 * \param end End position
 * \return Pointer to output_buffer
 * \copyright GPLv2+
 **/
char *SubString(const char *input_string, int buffer_length, char *output_buffer, int start, int end);

/**
 * Returns sign of a double.
 *
 * \copyright GPLv2+
 **/
int Sign(double arg);

/**
 * Returns square of a double.
 *
 * \copyright GPLv2+
 **/
double Sqr(double arg);

/**
 * Returns cube of a double.
 *
 * \copyright GPLv2+
 **/
double Cube(double arg);

/**
 * This function reduces angles greater than.
 *
 * \copyright GPLv2+
 **/
double FixAngle(double x);

/**
 * Returns angle in radians from argument in degrees.
 *
 * \copyright GPLv2+
 **/
double Radians(double arg);

/**
 * Returns angle in degrees from argument in radians.
 *
 * \copyright GPLv2+
 **/
double Degrees(double arg);

/**
 * Calculates scalar magnitude of a vector_t argument.
 *
 * \copyright GPLv2+
 **/
void Magnitude(vector_t *v);

/**
 * Adds vectors v1 and v2 together to produce v3.
 *
 * \copyright GPLv2+
 **/
void Vec_Add(vector_t *v1, vector_t *v2, vector_t *v3);

/**
 * Subtracts vector v2 from v1 to produce v3.
 *
 * \copyright GPLv2+
 **/
void Vec_Sub(vector_t *v1, vector_t *v2, vector_t *v3);

/**
 * Multiplies the vector v1 by the scalar k to produce the vector v2.
 *
 * \copyright GPLv2+
 **/
void Scalar_Multiply(double k, vector_t *v1, vector_t *v2);

/**
 * Multiplies the vector v1 by the scalar k.
 *
 * \copyright GPLv2+
 **/
void Scale_Vector(double k, vector_t *v);

/**
 * Returns the dot product of two vectors.
 *
 * \copyright GPLv2+
 **/
double Dot(vector_t *v1, vector_t *v2);

/**
 * Calculates the angle between vectors v1 and v2.
 *
 * \copyright GPLv2+
 **/
double Angle(vector_t *v1, vector_t *v2);

/**
 * Produces cross product of v1 and v2, and returns in v3.
 *
 * \copyright GPLv2+
 **/
void Cross(vector_t *v1, vector_t *v2 ,vector_t *v3);

/**
 * Normalizes a vector.
 *
 * \copyright GPLv2+
 **/
void Normalize(vector_t *v);

/**
 * Returns mod 2PI of argument.
 *
 * \copyright GPLv2+
 **/
double FMod2p(double x);

/* predict's old date/time management functions. */

/**
 * Converts the satellite's position and velocity vectors from normalized values to km and km/sec.
 *
 * \copyright GPLv2+
 **/
void Convert_Sat_State(double pos[3], double vel[3]);

/**
 * The function Julian_Date_of_Year calculates the Julian Date of Day 0.0 of {year}. This function is used to calculate the Julian Date of any date by using Julian_Date_of_Year, DOY, and Fraction_of_Day. Astronomical Formulae for Calculators, Jean Meeus, pages 23-25. Calculate Julian Date of 0.0 Jan year.
 *
 * \copyright GPLv2+
 **/
double Julian_Date_of_Year(double year);

/**
 * The function Julian_Date_of_Epoch returns the Julian Date of an epoch specified in the format used in the NORAD two-line element sets. It has been modified to support dates beyond the year 1999 assuming that two-digit years in the range 00-56 correspond to 2000-2056. Until the two-line element set format is changed, it is only valid for dates through 2056 December 31. Modification to support Y2K. Valid 1957 through 2056.
 *
 * \copyright GPLv2+
 **/
double Julian_Date_of_Epoch(double epoch);

/**
 * The function Delta_ET has been added to allow calculations on the position of the sun.  It provides the difference between UT (approximately the same as UTC) and ET (now referred to as TDT). This function is based on a least squares fit of data from 1950 to 1991 and will need to be updated periodically. Values determined using data from 1950-1991 in the 1990 Astronomical Almanac.  See DELTA_ET.WQ1 for details.
 *
 * \copyright GPLv2+
 **/
double Delta_ET(double year);

/**
 * Reference:  The 1992 Astronomical Almanac, page B6.
 *
 * \copyright GPLv2+
 **/
double ThetaG_JD(double jd);

/**
 * Calculates the day number from m/d/y. Needed for orbit_decay.
 *
 * \copyright GPLv2+
 **/
long DayNum(int month, int day, int year);

/**
 * Procedure Calculate_LatLonAlt will calculate the geodetic position of an object given its ECI position pos and time. It is intended to be used to determine the ground track of a satellite.  The calculations  assume the earth to be an oblate spheroid as defined in WGS '72. Reference:  The 1992 Astronomical Almanac, page K12.
 *
 * \copyright GPLv2+
 **/
void Calculate_LatLonAlt(double time, const double pos[3], geodetic_t *geodetic);

/**
 * The procedures Calculate_Obs and Calculate_RADec calculate
 * the *topocentric* coordinates of the object with ECI position,
 * {pos}, and velocity, {vel}, from location {geodetic} at {time}.
 * The {obs_set} returned for Calculate_Obs consists of azimuth,
 * elevation, range, and range rate (in that order) with units of
 * radians, radians, kilometers, and kilometers/second, respectively.
 * The WGS '72 geoid is used and the effect of atmospheric refraction
 * (under standard temperature and pressure) is incorporated into the
 * elevation calculation; the effect of atmospheric refraction on
 * range and range rate has not yet been quantified.
 * The {obs_set} for Calculate_RADec consists of right ascension and
 * declination (in that order) in radians.  Again, calculations are
 * based on *topocentric* position using the WGS '72 geoid and
 * incorporating atmospheric refraction.
 *
 * \copyright GPLv2+
 **/
void Calculate_Obs(double time,  const double pos[3], const double vel[3], geodetic_t *geodetic, vector_t *obs_set);

/**
 * This function is used in the FindMoon() function.
 *
 * \copyright GPLv2+
 **/
double PrimeAngle(double x);

/**
 *
 * \copyright GPLv2+
 **/
void Calculate_User_PosVel(double time, geodetic_t *geodetic, double obs_pos[3], double obs_vel[3]);


#endif
