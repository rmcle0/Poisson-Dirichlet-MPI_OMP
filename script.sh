#!/bin/bash

GRID_VALS=(80 160)
MPI_NUMS=(1 2 4)
THREAD_NUMS=(1 2 4 8)

OUT_DIR="./TABLE_1"
ERR_DIR="./ERRORS"

EXE="./mpi_openmp"

mkdir -p "$OUT_DIR"
mkdir -p "$ERR_DIR"

for M in ${GRID_VALS[@]}; do
	for mpi_proc in ${MPI_NUMS[@]}; do
		for threads in ${THREAD_NUMS[@]}; do
			if [[ $M -eq 80 ]]
			then
				N=90
			fi

			if [[ $M -eq 160 ]]
			then
				N=180
			fi

			OUT_FILE="$OUT_DIR/para_${M}_${N}_${mpi_proc}_${threads}.txt"
			ERR_FILE="$ERR_DIR/err_para_${M}_${N}_${mpi_proc}_${threads}.txt"
			
			cores=2

			if [[ $mpi_proc -eq 8 ]]
			then
				cores=1
			fi
			
			bsub  -x -n ${mpi_proc} -W 00:30 -q short -oo "$OUT_FILE" -eo "$ERR_FILE"  -R "span[hosts=1] affinity[core(${cores})]" mpiexec ${EXE} "$M" "$N"  "$threads"
			
			#mpisubmit.pl -p ${mpi_proc}  -t ${threads} --stdout "$OUT_FILE" --stderr "$ERR_FILE"  mpi_openmp -- "$M" "$N" "$threads"
			
			echo "Launched ${EXE}, M=$M, N=$N, delta=1e$delta_power, on $threads threads"

		done
	done
done


#mpisubmit.pl -p 4  -t 4 --stdout TEST --stderr ERROR  mpi_openmp -- 160 180 4
bsub -n 2 -W 00:30 -q short -oo "TEST" -eo "ERROR" -R "span[hosts=1]"  "./mpi_openmp" 80 90 2