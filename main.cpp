#include <mpi.h>

#include <chrono>

#include "matrix_filling.h"

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

//  try {
    // Get mesh params from the command line
    int M = std::stoi(argv[1]);
    int N = std::stoi(argv[2]);
    int THREADS = std::stoi(argv[3]);

    FinDiffSolution solution(M, N, THREADS, MPI_COMM_WORLD);

    int size, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    solution.SetMPISize(size);
    solution.SetMPIRank(rank);

    std::chrono::steady_clock::time_point begin =
        std::chrono::steady_clock::now();

    solution.solve();

    MPI_Barrier(MPI_COMM_WORLD);

    std::chrono::steady_clock::time_point end =
        std::chrono::steady_clock::now();
    std::cout << "Time difference = " << M << " " << N << " " << std::to_string(size) << " " << std::to_string(THREADS) << " "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                       begin)
                     .count()
              << "[ms]" << std::endl;

    //solution.dump_Xij();
//  } catch (...) {
//    std::cerr << "exception was caught, exiting" << std::endl;
//  };


  MPI_Finalize();
}