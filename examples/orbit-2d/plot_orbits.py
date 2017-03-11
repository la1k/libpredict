#!/usr/bin/env python

import sys

if len(sys.argv) < 2:
    print("Usage: plot_orbits <input_data_files>")
    sys.exit()

import matplotlib.pyplot as plt
import numpy as np
import re
from mpl_toolkits.mplot3d import Axes3D

#2d plots
f_2d, (ax1, ax2) = plt.subplots(2, figsize=(15,15));#, sharex=True, sharey=True)

#3d plot
f_3d = plt.figure();
ax_3d = f_3d.add_subplot(111, projection='3d')

#regexp for extracting satellite name from filename
prog = re.compile("sat_track_(.*)\.dat");

for i in range(1, len(sys.argv)):
    sat_data = np.loadtxt(sys.argv[i]);
    satname = prog.match(sys.argv[i]).group(1);

    #orbits projected onto xy and yz planes
    ax1.plot(sat_data[:,0], sat_data[:,1], label=satname);
    ax2.plot(sat_data[:,1], sat_data[:,2], label=satname);

    #orbits in 3d plot
    ax_3d.plot(sat_data[:,0], sat_data[:,1], sat_data[:,2], label=satname);

#set aspect ratios of 2d plot
ax1.set_aspect(1);
ax2.set_aspect(1);

#set fixed limits of 2d plot
range_start = -100000;
range_y_end = 100000;
range_x_end = 200000;

ax1.set_xlim(range_start, range_x_end);
ax2.set_xlim(range_start, range_x_end);

ax1.set_ylim(range_start, range_y_end);
ax2.set_ylim(range_start, range_y_end);
ax1.legend(loc='best', fontsize=10);
ax2.legend(loc='best', fontsize=10);

f_2d.savefig("2d_plot.png", bbox_inches='tight');
f_3d.savefig("3d_plot.png", bbox_inches='tight');


