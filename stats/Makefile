# run this like this:
# CPATH=/path/to/aqmt/common make

CPP=g++
AR=ar

all: calc_mix calc_qpd compl pcapstats

calc_mix: calc_mix.cpp Makefile
	$(CPP) calc_mix.cpp -std=c++11 -O3 -o $@
calc_qpd: calc_queue_packets_drops.cpp Makefile
	$(CPP) calc_queue_packets_drops.cpp -std=c++11 -O3 -o $@
compl: compl.cpp Makefile
	$(CPP) compl.cpp -O3 -o $@
pcapstats: pcapstats.cpp Makefile
	$(CPP) pcapstats.cpp -std=c++11 -lpcap -O3 -o $@

clean:
	rm -rf calc_mix calc_qpd compl pcapstats *.a *.o
