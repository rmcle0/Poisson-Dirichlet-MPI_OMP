all: 	clean  consequential openmp


mpi_openmp: main.cpp matrix_filling.cpp matrix_filling.h
	mpicxx -fopenmp -O3 -std=c++11 -g3  main.cpp matrix_filling.cpp  -o mpi_openmp





init: 
	mkdir -p archive
	mkdir -p TABLE_1

launch_mpi: mpi
	mpisubmit.pl --stdout "mpi_40_40_1" -p 1 ./mpi -- 40 40 1
	mpisubmit.pl --stdout "mpi_40_40_2" -p 2 ./mpi -- 40 40 1
	mpisubmit.pl --stdout "mpi_40_40_4" -p 4 ./mpi -- 40 40 1




mpi: main.cpp matrix_filling.cpp matrix_filling.h
	mpicxx -O3 -std=c++11 -g3  main.cpp matrix_filling.cpp  -o mpi 




launch_mpiomp: mpi_openmp
	mpisubmit.pl --stdout "mpiomp_40_40_4_1" -p 4 -t 1 ./mpi_openmp -- 40 40 1
	mpisubmit.pl --stdout "mpiomp_40_40_4_2" -p 4 -t 2 ./mpi_openmp -- 40 40 2
	mpisubmit.pl --stdout "mpiomp_40_40_4_4" -p 4 -t 4 ./mpi_openmp -- 40 40 4
	mpisubmit.pl --stdout "mpiomp_40_40_4_8" -p 4 -t 8 ./mpi_openmp -- 40 40 8







rmpi: mpi
	mpirun ./mpi 40 40 1

dbg:
	#--log-file="dbg.log"
	mpirun -n 4 valgrind  ./mpi 10 10 1

consequential: matrix_filling.cpp matrix_filling.h main.cpp
	g++ -std=c++11 -g3  main.cpp matrix_filling.cpp  -o consequential


openmp:  matrix_filling.cpp matrix_filling.h main.cpp
	g++ -g3 -fopenmp -std=c++11 -g3  main.cpp matrix_filling.cpp  -o openmp


comparator: comparator.cpp
	g++ -std=c++11  comparator.cpp -o comparator

launch_compare_tasks:  consequential openmp 
	bsub -n 1 -W 00:30 -q short  -R "span[hosts=1]" OMP_NUM_THREADS="1" ./consequential "40" "40" "1"
	bsub -n 1 -W 00:30 -q short  -R "span[hosts=1]" OMP_NUM_THREADS="1" ./openmp "40" "40" "1"
	bsub -n 1 -W 00:30 -q short  -R "span[hosts=1]" OMP_NUM_THREADS="4" ./openmp "40" "40" "4"
	bsub -n 1 -W 00:30 -q short  -R "span[hosts=1]" OMP_NUM_THREADS="16" ./openmp "40" "40" "16"

compare: comparator
	./comparator ./TABLE_1/CONS_40_40_1  ./TABLE_1/OMP_40_40_1 
	./comparator  ./TABLE_1/CONS_40_40_1  ./TABLE_1/OMP_40_40_4
	./comparator  ./TABLE_1/CONS_40_40_1  ./TABLE_1/OMP_40_40_16 
	

omp:
	./task1_openmp

archive: init
	mv *out *err archive

send:
	scp -i /home/sdev/.ssh/id_rsa_hpc  *.cpp *.h *lsf Makefile script.sh edu-cmc-skmodel24-607-01@polus.hpc.cs.msu.ru:/home_edu/edu-cmc-skmodel24-607/edu-cmc-skmodel24-607-01/task1


lsf:
	scp -i /home/sdev/.ssh/id_rsa_hpc  my_task.lsf edu-cmc-skmodel24-607-01@polus.hpc.cs.msu.ru:/home_edu/edu-cmc-skmodel24-607/edu-cmc-skmodel24-607-01/task1


script:
	scp -i /home/sdev/.ssh/id_rsa_hpc  script.sh edu-cmc-skmodel24-607-01@polus.hpc.cs.msu.ru:/home_edu/edu-cmc-skmodel24-607/edu-cmc-skmodel24-607-01/task1



clean:
	rm -f task1_openmp task1


