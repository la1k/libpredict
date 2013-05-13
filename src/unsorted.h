#ifndef UNSORTED_H_
#define UNSORTED_H_

#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <curses.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <pthread.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <unistd.h>
#include <fcntl.h>
#include <termios.h>

//#include <TLE.h>
#include "vec3.h"

/* Constants used by SGP4/SDP4 code */
#define	km2mi		0.621371		/* km to miles */
#define deg2rad		1.745329251994330E-2	/* Degrees to radians */
#define pi		3.14159265358979323846	/* Pi */
#define pio2		1.57079632679489656	/* Pi/2 */
#define x3pio2		4.71238898038468967	/* 3*Pi/2 */
#define twopi		6.28318530717958623	/* 2*Pi  */
#define e6a		1.0E-6
#define tothrd		6.6666666666666666E-1	/* 2/3 */
#define xj2		1.0826158E-3		/* J2 Harmonic (WGS '72) */
#define xj3		-2.53881E-6		/* J3 Harmonic (WGS '72) */   
#define xj4		-1.65597E-6		/* J4 Harmonic (WGS '72) */
#define xke		7.43669161E-2
#define xkmper		6.378137E3		/* WGS 84 Earth radius km */
#define xmnpda		1.44E3			/* Minutes per day */
#define ae		1.0
#define ck2		5.413079E-4
#define ck4		6.209887E-7
#define f		3.35281066474748E-3	/* Flattening factor */
#define ge		3.986008E5 	/* Earth gravitational constant (WGS '72) */
#define s		1.012229
#define qoms2t		1.880279E-09
#define secday		8.6400E4	/* Seconds per day */
#define omega_E		1.00273790934	/* Earth rotations/siderial day */
#define omega_ER	6.3003879	/* Earth rotations, rads/siderial day */
#define zns		1.19459E-5
#define c1ss		2.9864797E-6
#define zes		1.675E-2
#define znl		1.5835218E-4
#define c1l		4.7968065E-7
#define zel		5.490E-2
#define zcosis		9.1744867E-1
#define zsinis		3.9785416E-1
#define zsings		-9.8088458E-1
#define zcosgs		1.945905E-1
#define zcoshs		1
#define zsinhs		0
#define q22		1.7891679E-6
#define q31		2.1460748E-6
#define q33		2.2123015E-7
#define g22		5.7686396
#define g32		9.5240898E-1
#define g44		1.8014998
#define g52		1.0508330
#define g54		4.4108898
#define root22		1.7891679E-6
#define root32		3.7393792E-7
#define root44		7.3636953E-9
#define root52		1.1428639E-7
#define root54		2.1765803E-9
#define thdt		4.3752691E-3
#define rho		1.5696615E-1
#define mfactor		7.292115E-5
#define sr		6.96000E5	/* Solar radius - km (IAU 76) */
#define AU		1.49597870691E8	/* Astronomical unit - km (IAU 76) */

//Seems to be used as a container for processed tle data from
//the original TLE.
typedef struct	{
	double epoch;
	double xndt2o;
	double xndd6o;
	double bstar;
	double xincl;
	double xnodeo;
	double eo;
	double omegao;
	double xmo;
	double xno;
 	int catnr;
	int elset;
	int revnum;
}tle_t; 


/* Geodetic position structure used by SGP4/SDP4 code. */

typedef struct	{
		   double lat, lon, alt, theta;
		}  geodetic_t;

/* General three-dimensional vector structure used by SGP4/SDP4 code. */

typedef struct	{
		   double x, y, z, w;
		}  vector_t;

/* Common arguments between deep-space functions used by SGP4/SDP4 code. */

typedef struct	{
		   	   /* Used by dpinit part of Deep() */
		   double  eosq, sinio, cosio, betao, aodp, theta2,
			   sing, cosg, betao2, xmdot, omgdot, xnodot, xnodp;
	   
			   /* Used by dpsec and dpper parts of Deep() */
		   double  xll, omgadf, xnode, em, xinc, xn, t;
    
		 	   /* Used by thetg and Deep() */
		   double  ds50;
		}  deep_arg_t;

char *SubString(const char *string, int start, int end);

/* Returns sign of a double */
int Sign(double arg);

/* Returns square of a double */
double Sqr(double arg);

/* Returns cube of a double */
double Cube(double arg);

/* This function reduces angles greater than*/
double FixAngle(double x);

/* Returns angle in radians from argument in degrees */
double Radians(double arg);

/* Returns angle in degrees from argument in radians */
double Degrees(double arg);

/* Returns the arcsine of the argument */
double ArcSin(double arg);

/* Returns arccosine of argument */
double ArcCos(double arg);

/* Calculates scalar magnitude of a vector_t argument */
void Magnitude(vector_t *v);

/* Adds vectors v1 and v2 together to produce v3 */
void Vec_Add(vector_t *v1, vector_t *v2, vector_t *v3);

	/* Subtracts vector v2 from v1 to produce v3 */
void Vec_Sub(vector_t *v1, vector_t *v2, vector_t *v3);

/* Multiplies the vector v1 by the scalar k to produce the vector v2 */
void Scalar_Multiply(double k, vector_t *v1, vector_t *v2);

/* Multiplies the vector v1 by the scalar k */
void Scale_Vector(double k, vector_t *v);

/* Returns the dot product of two vectors */
double Dot(vector_t *v1, vector_t *v2);

/* Calculates the angle between vectors v1 and v2 */
double Angle(vector_t *v1, vector_t *v2);

/* Produces cross product of v1 and v2, and returns in v3 */
void Cross(vector_t *v1, vector_t *v2 ,vector_t *v3);

/* Normalizes a vector */
void Normalize(vector_t *v);

/* Four-quadrant arctan function */
double AcTan(double sinx, double cosx);

/* Returns mod 2PI of argument */
double FMod2p(double x);

/* Returns arg1 mod arg2 */
double Modulus(double arg1, double arg2);

/* Returns fractional part of double argument */
double Frac(double arg);

/* Returns argument rounded up to nearest integer */
int Round(double arg);

/* Returns the floor integer of a double arguement, as double */
double Int(double arg);

/* Converts the satellite's position and velocity  */
/* vectors from normalized values to km and km/sec */ 
void Convert_Sat_State(double pos[3], double vel[3]);

/* The function Julian_Date_of_Year calculates the Julian Date  */
/* of Day 0.0 of {year}. This function is used to calculate the */
/* Julian Date of any date by using Julian_Date_of_Year, DOY,   */
/* and Fraction_of_Day. */
/* Astronomical Formulae for Calculators, Jean Meeus, */
/* pages 23-25. Calculate Julian Date of 0.0 Jan year */
double Julian_Date_of_Year(double year);

/* The function Julian_Date_of_Epoch returns the Julian Date of     */
/* an epoch specified in the format used in the NORAD two-line      */
/* element sets. It has been modified to support dates beyond       */
/* the year 1999 assuming that two-digit years in the range 00-56   */
/* correspond to 2000-2056. Until the two-line element set format   */
/* is changed, it is only valid for dates through 2056 December 31. */
/* Modification to support Y2K */
/* Valid 1957 through 2056     */
double Julian_Date_of_Epoch(double epoch);

/* The function DOY calculates the day of the year for the specified */
/* date. The calculation uses the rules for the Gregorian calendar   */
/* and is valid from the inception of that calendar system.          */
int DOY (int yr, int mo, int dy);

/* Fraction_of_Day calculates the fraction of */
/* a day passed at the specified input time.  */
double Fraction_of_Day(int hr, int mi, double se);

/* The function Julian_Date converts a standard calendar   */
/* date and time to a Julian Date. The procedure Date_Time */
/* performs the inverse of this function. */
double Julian_Date(struct tm *cdate);

/* The function Date_Time() converts a Julian Date to
standard calendar date and time. The function
Julian_Date() performs the inverse of this function. */
void Date_Time(double julian_date, struct tm *cdate);

	/* The function Delta_ET has been added to allow calculations on   */
	/* the position of the sun.  It provides the difference between UT */
	/* (approximately the same as UTC) and ET (now referred to as TDT).*/
	/* This function is based on a least squares fit of data from 1950 */
	/* to 1991 and will need to be updated periodically. */

	/* Values determined using data from 1950-1991 in the 1990 
	Astronomical Almanac.  See DELTA_ET.WQ1 for details. */
double Delta_ET(double year);

	/* The function ThetaG calculates the Greenwich Mean Sidereal Time */
	/* for an epoch specified in the format used in the NORAD two-line */
	/* element sets. It has now been adapted for dates beyond the year */
	/* 1999, as described above. The function ThetaG_JD provides the   */
	/* same calculation except that it is based on an input in the     */
	/* form of a Julian Date. */

	/* Reference:  The 1992 Astronomical Almanac, page B6. */
	/* Modification to support Y2K */
	/* Valid 1957 through 2056     */
double ThetaG(double epoch, deep_arg_t *deep_arg);

/* Reference:  The 1992 Astronomical Almanac, page B6. */
double ThetaG_JD(double jd);

double CurrentDaynum();

/* Procedure Calculate_LatLonAlt will calculate the geodetic  */
/* position of an object given its ECI position pos and time. */
/* It is intended to be used to determine the ground track of */
/* a satellite.  The calculations  assume the earth to be an  */
/* oblate spheroid as defined in WGS '72.                     */
/* Reference:  The 1992 Astronomical Almanac, page K12. */
void Calculate_LatLonAlt(double time, const double pos[3], geodetic_t *geodetic);

	
/* The procedures Calculate_Obs and Calculate_RADec calculate         */
/* the *topocentric* coordinates of the object with ECI position,     */
/* {pos}, and velocity, {vel}, from location {geodetic} at {time}.    */
/* The {obs_set} returned for Calculate_Obs consists of azimuth,      */
/* elevation, range, and range rate (in that order) with units of     */
/* radians, radians, kilometers, and kilometers/second, respectively. */
/* The WGS '72 geoid is used and the effect of atmospheric refraction */
/* (under standard temperature and pressure) is incorporated into the */
/* elevation calculation; the effect of atmospheric refraction on     */
/* range and range rate has not yet been quantified.                  */
/* The {obs_set} for Calculate_RADec consists of right ascension and  */
/* declination (in that order) in radians.  Again, calculations are   */
/* based on *topocentric* position using the WGS '72 geoid and        */
/* incorporating atmospheric refraction.                              */
void Calculate_Obs(double time,  const double pos[3], const double vel[3], geodetic_t *geodetic, vector_t *obs_set);

/* Reference:  Methods of Orbit Determination by  */
/*             Pedro Ramon Escobal, pp. 401-402   */
void Calculate_RADec(double time, const double pos[3], const double vel[3], geodetic_t *geodetic, vector_t *obs_set);

/* This function is used in the FindMoon() function. */
double PrimeAngle(double x);

void Calculate_User_PosVel(double time, geodetic_t *geodetic, double obs_pos[3], double obs_vel[3]);

/* Calculates solar position vector */
void Calculate_Solar_Position(double time, double solar_vector[3]);

/* Calculates satellite's eclipse status and depth */
bool Sat_Eclipsed(vector_t *pos, vector_t *sol, double *depth);

#endif
