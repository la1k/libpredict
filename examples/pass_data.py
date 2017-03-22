#!/usr/bin/env python3

import datetime
import numpy as np
import pypredict

# OSCAR 7 obtained at 2017-02-27
TLE = ( "1 07530U 74089B   17058.02491442 -.00000048  00000-0 -22049-4 0  9995",
        "2 07530 101.6163  28.9438 0011984 174.4353 227.0960 12.53625643935054")

OBSERVER_LATITUDE = 63.9
OBSERVER_LONGITUDE = 10.9
OBSERVER_HEIGHT = 10

if __name__ == "__main__":
    t = int(datetime.datetime(2017, 1, 28, 11, 30).timestamp())
    orb_ele = pypredict.predict_parse_tle(TLE[0], TLE[1])
    obs = pypredict.predict_create_observer("LA1K", np.radians(OBSERVER_LATITUDE), np.radians(OBSERVER_LONGITUDE), OBSERVER_HEIGHT)
    start_time = pypredict.predict_to_julian(t)
    
    aos_time = pypredict.predict_next_aos(obs, orb_ele, start_time)
    los_time = pypredict.predict_next_los(obs, orb_ele, start_time)

    aos_time = pypredict.predict_from_julian(aos_time)
    los_time = pypredict.predict_from_julian(los_time)

    aos_time = datetime.datetime.fromtimestamp(aos_time)
    los_time = datetime.datetime.fromtimestamp(los_time)

    print(aos_time)
    print(los_time)

