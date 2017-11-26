#
# Regular cron jobs for the libpredict2 package
#
0 4	* * *	root	[ -x /usr/bin/libpredict2_maintenance ] && /usr/bin/libpredict2_maintenance
