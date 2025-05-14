#!/bin/bash
# Set proj_root to the path leading to mpif_article
proj_root=/path/to/mpif_article
wksp_rt=mpif_example

models=(cohort iota many)
particles=(500 5000 10000)
fitr_types=(mpif pif)

for mod in "${models[@]}"; do
	for np in "${particles[@]}"; do
		for fitr in "${fitr_types[@]}"; do
			out_dir=$wksp_rt/$mod/np_$np/$fitr
			sbatch --export=mod=$mod,np=$np,fitr=$fitr,proj_root=$proj_root,wksp_rt=$wksp_rt,out_dir=$out_dir\
				--output=$proj_root/$out_dir/slurm-%j.out $proj_root/$wksp_rt/fit_models.sbat
		done
	done
done
