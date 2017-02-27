#!/usr/bin/env python

import sys

if len(sys.argv) < 2:
    print("Usage: plot_pass <input_data_file>")
    sys.exit()

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline

#######################
# Plot pass on a map. #
#######################

qth_lat = 63.422;
qth_lon = 10.39;

#prepare world map approximately limited to the pass
plt.figure(1);
map_width = 8000000;
map_height = 10000000;
basemap = Basemap(width=map_width,height=map_height,\
                    resolution='l',projection='gnom',\
                                        lat_0=qth_lat,lon_0=qth_lon)
basemap.drawcoastlines(linewidth=0.25)
basemap.drawcountries(linewidth=0.25)
basemap.drawmapboundary(fill_color='aqua')
basemap.fillcontinents(color='#cc9966',lake_color='#99ffff')

#load satellite data
sat_data = np.loadtxt(sys.argv[1]);
sat_data[sat_data[:,3] > 180, 3] -= 360;

sat_visible = sat_data[sat_data[:,2] > 0]; #nonnegative elevation

#plot satellite track on map
basemap.plot(sat_data[:,1], sat_data[:,0], linewidth=6.0, color='black', latlon=True, label='Elevation < 0')

#indicate in satellite track where the satellite is observable from the QTH coordinates
basemap.plot(sat_visible[:,1], sat_visible[:,0], linewidth=6.0, color='green', latlon=True, label='Elevation > 0')

elevation_thresh_10 = sat_visible[:,2] > 10;
basemap.plot(sat_visible[elevation_thresh_10,1], sat_visible[elevation_thresh_10,0], linewidth=4.0, color='yellow', latlon=True, label='Elevation > 10')

elevation_thresh_50 = sat_visible[:,2] > 50;
basemap.plot(sat_visible[elevation_thresh_50,1], sat_visible[elevation_thresh_50,0], linewidth=2.0, color='red', latlon=True, label='Elevation > 50')

#plot QTH point
qth_x,qth_y = basemap(qth_lon, qth_lat);
basemap.plot(qth_lon, qth_lat, 'bo', latlon=True);
plt.text(qth_x, qth_y, 'QTH',fontsize=14,
                            ha='right',va='bottom',color='black')
plt.legend(fontsize=10)

#save map
plt.savefig("map.png", bbox_inches='tight');

#########################
# Plot pass properties. #
#########################

f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)

#plot elevation and azimuth in simple plot
x = np.arange(0, len(sat_data[:,1]));
x = x/60.0; #convert to minutes
ax1.plot(x, sat_data[:,2], label="Elevation");
ax1.plot(x, sat_data[:,3], label="Azimuth");
ax1.plot(x, np.repeat(0, len(x)), label="Horizon", linestyle="dashed");
legendsize=10
ax1.legend(loc='best', fontsize=legendsize);
ax1.set_ylabel("Degrees");

#plot elevation rate
ax2.plot(x, sat_data[:,5], label="Elevation rate");

#plot numerical derivative estimated from an interpolating spline
f = InterpolatedUnivariateSpline(x, sat_data[:,2], k=1)
dfdx = f.derivative()
dydx = dfdx(x)
ax2.plot(x, sat_data[0,5]/dydx[0]*dydx, label="Scaled numerical derivative");
ax2.legend(loc='best', fontsize=legendsize);
ax2.set_ylabel("Degrees");
ax2.set_ylabel("Degrees");

#plot doppler shift
ax3.plot(x, np.repeat(1000, len(x)), label="Downlink frequency");
ax3.plot(x, sat_data[:,4], label="Doppler-shifted downlink frequency");
ax3.legend(loc='best', fontsize=legendsize);
ax3.set_ylabel("Frequency (MHz)");
tickrange=np.arange(999.98, 1000.025, 0.01);
ax3.set_yticks(tickrange)
ax3.set_yticklabels(map(str, tickrange));
plt.xlabel("Timestep (min)");
plt.savefig("pass.png");

#plot azimuth and elevation in a polar diagram
f, ax = plt.subplots(1, subplot_kw=dict(projection='polar'))
x = np.arange(0, len(sat_visible[:,2]));
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)
ax.plot(sat_visible[:,3]*np.pi/180.0, 90 - sat_visible[:,2]);
ax.set_yticks(range(0, 90, 10))                   # Define the yticks
ax.set_yticklabels(map(str, range(90, 0, -10)))   # Change the labels
plt.savefig("azimuth.png", bbox_inches='tight');
