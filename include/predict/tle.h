#ifndef _TLE_H_
#define _TLE_H_

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

#endif
