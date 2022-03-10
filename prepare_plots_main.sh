#!/bin/bash
. copy_data.sh


if [ "$#" != "2" ]; then
 	echo "usage: ./prepare_plots_main.sh <experiments folder> <plots folder>"
 	exit 65
fi

expfolder=$1
plotsfolder=$2

link_array=("4" "12" "40" "120" "200")
rtt_array=("5" "10" "20" "50" "100")

testcases_to_plot=("s1d0s1d0" "s1dhs1dh")
testcases_to_plot_tb=("s1d0s1d0" "s1dhs1dh" "s5d0s5d0")
testcases_to_plot_extra=("s1d0s1d0" "s2d0s2d0" "s3d0s3d0" "s4d0s4d0" "s5d0s5d0" "s6d0s6d0" "s7d0s7d0" "s8d0s8d0" "s9d0s9d0" "s10d0s10d0"
                         "s0d0s10d0" "s1d0s9d0" "s2d0s8d0" "s3d0s7d0" "s4d0s6d0" "s6d0s4d0" "s7d0s3d0" "s8d0s2d0" "s9d0s1d0" "s10d0s0d0")

testcase_id=""

function copy_data_m() {
	for aqm in "${aqm_array[@]}"; do
		aqmname="m_${aqm}"
		copy_data_line_m
	done
}

function copy_data_mrtt2() {
	for aqm in "${aqm_array[@]}"; do
		aqmname=mr2_${aqm}
		copy_data_line_mr2
	done
}
function copy_data_m_extra() {
	declare -i index=0
	for aqm in "${aqm_array[@]}"; do

		aqmname="m_${aqm}"
		copy_data_link40_rtt10_extra
		index=$index+1
	done
}

function copy_data_overload() {
	declare -i index=0
	for aqm in "${aqm_array[@]}"; do
		aqmname="o_${aqm}"
		copy_data_line_o
	done
}


declare -a aqm_array=( "pcdualpi2" "ecpie" "ecfqcodel" )

targetfolder="${plotsfolder}/plots_er"
mainfolder="${expfolder}/res_er"
mkdir ${targetfolder}/data
copy_data_m

declare -a aqm_array=( "pcdualpi2" "ecpie" "ecfqcodel" "prdualpi2" "erpie" "erfqcodel")
targetfolder="${plotsfolder}/plots_mr"
mainfolder="${expfolder}/res_mr"
mkdir ${targetfolder}/data
copy_data_mrtt2

targetfolder="${plotsfolder}/plots_extra"
mainfolder="${expfolder}/res_tb"
mkdir ${targetfolder}/data
copy_data_m_extra

aqm_array=( "pcdualpi2" )
targetfolder="${plotsfolder}/plots_overload"
mainfolder="${expfolder}/res_overload"
mkdir ${targetfolder}/data
copy_data_overload

declare -a aqm_array=( "ppdualpi2" "eefqcodel1ms" "eefqcodel5ms" )
testcases_to_plot=("s1d0s0d0")
targetfolder="${plotsfolder}/plots_er"
mainfolder="${expfolder}/res_er"
copy_data_m
