#!/bin/bash

function generate_satellite_testcase(){
	tle_file="$1"
	qth_file="$2"
	satellite_name="$3"
	start_time="$4"
	tot_secs="$5"
	testcase_filename="$6"

	#parse tle and qth information into testcase file
	echo "[tle]" > $testcase_filename
	grep -A 2 $satellite_name $tle_file | tail -2 >> $testcase_filename

	echo "" >> $testcase_filename
	echo "[qth]" >> $testcase_filename
	echo "lat=$(sed '2!d' $qth_file | sed -rn 's/\s//p')" >> $testcase_filename
	echo "lon=$(sed '3!d' $qth_file | sed -rn 's/\s//p')" >> $testcase_filename
	echo "alt=$(sed '4!d' $qth_file | sed -rn 's/\s//p')" >> $testcase_filename
	echo "freq=0" >> $testcase_filename

	#predict orbit
	echo "" >> $testcase_filename
	echo "[data]" >> $testcase_filename
	faketime "$start_time" predict -t $tle_file -q $qth_file -s &
	steps=$(seq 1 $tot_secs)
	for timestemp in $steps; do
		sleep 1
		predict_response=($(./predict_client "GET_SAT $satellite_name"))
		time=$(./predict_client "GET_TIME")
		doppler_shift=($(./predict_client "GET_DOPPLER $satellite_name"))
		satname=${predict_response[0]}
		lon=${predict_response[1]}
		lat=${predict_response[2]}
		az=${predict_response[3]}
		el=${predict_response[4]}
		aostime=${predict_response[5]}
		footprint=${predict_response[6]}
		range=${predict_response[7]}
		alt=${predict_response[8]}
		vel=${predict_response[9]}
		orbitnum=${predict_response[10]}
		visibility=${predict_response[11]}
		phase=${predict_response[12]}
		eclipse_depth=${predict_response[13]}
		squint=${predict_response[14]}
		echo "$time $lat $lon $alt $az $el $doppler_shift" >> $testcase_filename
	done
	killall predict
	sleep 1
}

function generate_sun_testcase(){
	start_time="$1"
	tot_secs="$2"
	qth_file="$3"
	testcase_filename="$4"

	#parse qth information into testcase file
	echo "[qth]" > $testcase_filename
	echo "lat=$(sed '2!d' $qth_file | sed -rn 's/\s//p')" >> $testcase_filename
	echo "lon=$(sed '3!d' $qth_file | sed -rn 's/\s//p')" >> $testcase_filename
	echo "alt=$(sed '4!d' $qth_file | sed -rn 's/\s//p')" >> $testcase_filename

	faketime "$start_time" predict -q $qth_file -s &

	echo "" >> $testcase_filename
	echo "[data]" >> $testcase_filename
	steps=$(seq 1 $tot_secs)
	for timestemp in $steps; do
		sleep 1
		predict_response=($(./predict_client "GET_SUN"))
		time=$(./predict_client "GET_TIME")
		lon=${predict_response[0]}
		lat=${predict_response[1]}
		az=${predict_response[2]}
		el=${predict_response[3]}
		alt=${predict_response[4]}
		echo $time $lat $lon $alt $az $el >> $testcase_filename
	done
	killall predict
	sleep 1
}

function generate_moon_testcase(){
	start_time="$1"
	tot_secs="$2"
	qth_file="$3"
	testcase_filename="$4"

	#parse qth information into testcase file
	echo "[qth]" > $testcase_filename
	echo "lat=$(sed '2!d' $qth_file | sed -rn 's/\s//p')" >> $testcase_filename
	echo "lon=$(sed '3!d' $qth_file | sed -rn 's/\s//p')" >> $testcase_filename
	echo "alt=$(sed '4!d' $qth_file | sed -rn 's/\s//p')" >> $testcase_filename

	faketime "$start_time" predict -q $qth_file -s &

	echo "" >> $testcase_filename
	echo "[data]" >> $testcase_filename
	steps=$(seq 1 $tot_secs)
	for timestemp in $steps; do
		sleep 1
		predict_response=($(./predict_client "GET_MOON"))
		time=$(./predict_client "GET_TIME")
		lon=${predict_response[0]}
		lat=${predict_response[1]}
		az=${predict_response[2]}
		el=${predict_response[3]}
		alt=${predict_response[4]}
		echo $time $lat $lon $alt $az $el >> $testcase_filename
	done
	killall predict
	sleep 1
}


killall predict
generate_satellite_testcase "predict.tle" "predict.qth" "OSCAR-7" "2015-09-20 19:33" "5" "testcase"
generate_sun_testcase "2015-09-20 19:33" "5" "predict.qth" sun_testcase.dat
generate_moon_testcase "2015-09-20 19:33" "5" "predict.qth" moon_testcase.dat
