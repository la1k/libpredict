#!/usr/bin/env python3

import datetime
import numpy as np
import pypredict
import pytz

# OSCAR 7 obtained at 2017-02-27
TLE = ( "1 07530U 74089B   17058.02491442 -.00000048  00000-0 -22049-4 0  9995",
        "2 07530 101.6163  28.9438 0011984 174.4353 227.0960 12.53625643935054")

OBSERVER_LATITUDE = 63.9
OBSERVER_LONGITUDE = 10.9
OBSERVER_HEIGHT = 10

if __name__ == "__main__":
    t = datetime.datetime(2017, 2, 28, 11, 30, tzinfo=pytz.utc)
    print('search start:', t)
    t = int(t.timestamp())
    print('search start unix timestamp:', t)
    orb_ele = pypredict.predict_parse_tle(TLE[0], TLE[1])
    obs = pypredict.predict_create_observer("LA1K", np.radians(OBSERVER_LATITUDE), np.radians(OBSERVER_LONGITUDE), OBSERVER_HEIGHT)
    start_time = pypredict.predict_to_julian(t)
    print('search start julian timestamp:', start_time)
    
    aos_time = pypredict.predict_next_aos(obs, orb_ele, start_time)
    los_time = pypredict.predict_next_los(obs, orb_ele, start_time)
    print('julian timestamps: ', aos_time, los_time)

    aos_time = pypredict.predict_from_julian(aos_time)
    los_time = pypredict.predict_from_julian(los_time)
    print('unix timestamps: ', aos_time, los_time)

    aos_time = datetime.datetime.utcfromtimestamp(aos_time)
    los_time = datetime.datetime.utcfromtimestamp(los_time)

    print('datetime objects: ',aos_time, los_time)
    print('reference from c implementation: 12:56:07 og 13:18:28')

