#!/bin/bash

(cd exp_source

for ccaqm in "pcdualpi2" "ecpie" "ecfqcodel"; do for resfolder in "res_er" "res_mr" "res_tb" "res_ccdf"; do
	tar -xvf ${ccaqm}_${resfolder}.tar.gz -C .
done; done

for ccaqm in "prdualpi2" "erpie" "erfqcodel"; do for resfolder in "res_mr" "res_tb"; do
        tar -I 'gzip -9' -cf ${ccaqm}_${resfolder}.tar.gz ${resfolder}/*${ccaqm}*
done; done

tar -xvf pcdualpi2_res_overload.tar.gz -C .
tar -xvf ppdualpi2_res_er.tar.gz -C .)
