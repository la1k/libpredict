#include "sun.h"
#include "unsorted.h"
#include "defs.h"

void sun_predict(double time, double position[3])
{
	double jul_utc = time + 2444238.5;
	double mjd = jul_utc - 2415020.0;
	double year = 1900 + mjd / 365.25;
	double T = (mjd + Delta_ET(year) / secday) / 36525.0;
	double M = Radians(Modulus(358.47583+Modulus(35999.04975*T,360.0)-(0.000150+0.0000033*T)*Sqr(T),360.0));
	double L = Radians(Modulus(279.69668+Modulus(36000.76892*T,360.0)+0.0003025*Sqr(T),360.0));
	double e = 0.01675104-(0.0000418+0.000000126*T)*T;
	double C = Radians((1.919460-(0.004789+0.000014*T)*T)*sin(M)+(0.020094-0.000100*T)*sin(2*M)+0.000293*sin(3*M));
	double O = Radians(Modulus(259.18-1934.142*T,360.0));
	double Lsa = Modulus(L+C-Radians(0.00569-0.00479*sin(O)), 2*M_PI);
	double nu = Modulus(M+C, 2*M_PI);
	double R = 1.0000002*(1.0-Sqr(e))/(1.0+e*cos(nu));
	double eps = Radians(23.452294-(0.0130125+(0.00000164-0.000000503*T)*T)*T+0.00256*cos(O));
	R = AU*R;

	position[0] = R*cos(Lsa);
	position[1] = R*sin(Lsa)*cos(eps);
	position[2] = R*sin(Lsa)*sin(eps);
}
