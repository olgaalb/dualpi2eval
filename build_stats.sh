#!/bin/bash

(cd stats && make)
ln -s stats/calc_mix stats/calc_qpd stats/compl stats/pcapstats .
