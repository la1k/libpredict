#include <stdio.h>
#include <math.h>
#include <unistd.h>

#include <predict/predict.h>
#include <predict/observer.h>

double observer_next_sunset(const observer_t *observer, double time)
{
	struct observation sun;
	
	// Scan for first elevation > 0 (t0)
	double t0 = time;
	observer_find_sun(observer, t0, &sun);
	while (sun.elevation < 0) {
		t0 += 1.0/24.0/3600.0;
		observer_find_sun(observer, t0, &sun);
	}

	//Scan for elevation < 0 (t1);
	double t1 = t0 + 1.0/24.0/3600.0;
	observer_find_sun(observer, t1, &sun);
	while (sun.elevation > 0) {
		t1 += 1.0/24.0/3600.0;
		observer_find_sun(observer, t1, &sun);
	}
	
//	int i = 0;
	while ( fabs(t0-t1) > 0.01/24/3600 ) {

		time = (t0 + t1) / 2.0;
		observer_find_sun(observer, time, &sun);

		if (sun.elevation > 0) t0 = time;
		else t1 = time;

//		printf("iteration %i: elev=%f\n", i++, sun.elevation*180.0/M_PI);
	}

	return time;
}

double observer_next_sunrise(const observer_t *observer, double time)
{
	struct observation sun;
	
	// Scan for first elevation < 0 (t0)
	double t0 = time;
	observer_find_sun(observer, t0, &sun);
	while (sun.elevation > 0) {
		t0 += 1.0/24.0/3600.0;
		observer_find_sun(observer, t0, &sun);
	}

	//Scan for elevation > 0 (t1);
	double t1 = t0 + 1.0/24.0/3600.0;
	observer_find_sun(observer, t1, &sun);
	while (sun.elevation < 0) {
		t1 += 1.0/24.0/3600.0;
		observer_find_sun(observer, t1, &sun);
	}
	
//	int i = 0;
	while ( fabs(t0-t1) > 0.01/24/3600 ) {

		time = (t0 + t1) / 2.0;
		observer_find_sun(observer, time, &sun);

		if (sun.elevation < 0) t0 = time;
		else t1 = time;

	//	printf("iteration %i: elev=%f\n", i++, sun.elevation*180.0/M_PI);
	}

	return time;
}

int main(int argc, char **argv)
{

	// Create observer object
	observer_t *obs = observer_create("Me", 59.95*M_PI/180.0, 10.75*M_PI/180.0, 0);
	if (!obs) {
		fprintf(stderr, "Failed to initialize observer!");
		exit(1);
	}

	double sunset = observer_next_sunset(obs, CurrentDaynum());

	// Convert to hour, minute, seconds
	double timeto = (sunset - CurrentDaynum())*24*3600;
	int h = timeto / 3600;
	int m = (timeto-h*3600) / 60;
	int s = ((int)timeto)%60;
	
	time_t t = (time_t)(86400.0 * (sunset + 3651.0));
	printf("Next sunset in %02i:%02i:%02i at UTC %s", h, m, s, asctime(gmtime(&t)));
	
	double sunrise = observer_next_sunrise(obs, CurrentDaynum());

	// Convert to hour, minute, seconds
	timeto = (sunrise - CurrentDaynum())*24*3600;
	h = timeto / 3600;
	m = (timeto-h*3600) / 60;
	s = ((int)timeto)%60;
	
	t = (time_t)(86400.0 * (sunrise + 3651.0));
	printf("Next sunrise in %02i:%02i:%02i at UTC %s", h, m, s, asctime(gmtime(&t)));
	
	// Free memory
	observer_destroy(obs);

	return 0;
}
