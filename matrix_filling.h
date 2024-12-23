#include <iostream>
#include <array>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <cassert>
#include <unordered_map>
#include <map>
#include <deque>
#include <cstdio>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <chrono>
#include <string>

#include <omp.h>
#include <mpi.h>

using namespace std;

#define EPS 0.01
#define NANTHROW(x)                              \
    if (std::isnan(x))                           \
    {                                            \
        throw std::runtime_error("nan occured"); \
    }
    

#define PRINTFOVER
#define DEBUGOVER


#define DELTA   0.00000003
#define ITMAX   10000









//Checked
struct point{
    double x;
    double y;

    point(double x, double y) : x(x), y(y) {         };

    struct point operator +(struct point& p){
        return point(x + p.x, y + p.y);
    };

    struct point operator -(struct point& p){
        return point(x - p.x, y - p.y);
    };

    struct point operator /(double t){
        return point(x / t, y / t);
    };

    double length(){
        return sqrt(abs(x * x + y * y));
    }

    //avoid costly root calculations
    double manhattan_length(){
        return abs(x) + abs(y);
    }
};

typedef struct point Point;




template<class T>
class Matrix {
  private:
    std::vector<T> data;

    int M = 0, N = 0; //real size of allocated (domain size)

 public:
  Matrix(Matrix& a) = delete;

  Matrix(Matrix&& a){
    data = std::move(a.data);
    M = a.M;
    N = a.N;
  }

  Matrix& operator=(Matrix&& other)
    {
    M = other.M;
    N = other.N;

    data = std::move(other.data);
    return *this;
    }




  Matrix(const int &rows, const int &cols, const T& val ) {
    M = rows; 
    N = cols;
    Reset(rows, cols, val);
  }

  void Reset(const int &rows, const int &cols, const T& val) {
    M = rows;
    N = cols;

    data.resize(M * N, val);
  }

  T At(const int &row, const int &col) const {
    return data[row * N + col];
  }

  T& At(const int &row, const int &col) {
    return data[row * N + col];
  }

  inline T& At_WithShadow(const int &row, const int &col) {
    //i th in original
    // is (i+1) in shadow
    return data[ (row + 1) * N   + (col + 1) + 1];
  }

  int GetNumRows() const {
    return M;
  }

  int GetNumColumns() const {
    return N;
  }

    void show(){
    for(int j = N - 1; j >= 0; j--){
        for(int i = 0; i < M; i++){
            printf("%10.5f  ", data[i * N + j]);
            }   
        std::cout << std::endl;
        }
  }

  
};



//For five-point equation
class SparseRows{
 
 public:
  enum Direction{Center, Up, Down, Left, Right, Q_koeff_in_one_eq};

  SparseRows(int num_equations, int domain_size_x, int domain_size_y):domain_size_x(domain_size_x), domain_size_y(domain_size_y){
    sparse_equations.resize(num_equations * Q_koeff_in_one_eq, 0);
    sparse_equations.resize(num_equations * Q_koeff_in_one_eq, 0);
  }

  void set(int x, int y, Direction dir, double num){
    int eq_num = x * domain_size_y + y;
    sparse_equations[eq_num * Q_koeff_in_one_eq + dir] = num;
  }

  inline double get(int x, int y, Direction dir){
    int eq_num = x * domain_size_y + y;
    return sparse_equations[eq_num * Q_koeff_in_one_eq + dir];
  }

  
  //multiplies given (i, j) equation with the 
  inline  double multiply(int i, int j, Matrix<double>& w){
    double res = 0;
    
    res += get(i, j, Center) * w.At_WithShadow(i, j) + 
           get(i, j, Right) * w.At_WithShadow(i + 1, j) + 
           get(i, j, Left) * w.At_WithShadow(i - 1, j) + 
           get(i, j, Up) * w.At_WithShadow(i, j + 1) +
             get(i, j, Down) * w.At_WithShadow(i, j - 1);  

    return res;
  }

  
    std::vector<double> sparse_equations;

    int domain_size_x;
    int domain_size_y;
};



class FinDiffSolution
{
public:
    bool verbose = true;

    double mesh_eps = 0.01;

    int M = 10;
    int N = 10;
    int THREADS = 1;

    double A1 = -1.5, B1 = 2, A2 = -1.5, B2 = 2;
    double h1 = (B1 - A1) / M;
    double h2 = (B2 - A2) / N;

    
    //int ind(int i, int j) {return i * N + M;};
    
    //1..M-1 , 1..N-1
    int map_to_line(int i, int j) {
        return (j - 1) * (M - 1) + (i - 1);
    };


    Matrix<double> X;


    //MPI
    int domain_coord_along_x = 0, domain_coord_along_y = 0;
    int base_domain_size_x = 0, base_domain_size_y = 0;
    int domains_quantity_along_x = 0, domains_quantity_along_y = 0;

    int ind_of_smaller_domains_X = 0;
    int ind_of_smaller_domains_Y = 0;

    int current_domain_size_x, current_domain_size_y;
    int current_domain_size;

    int mpi_rank = 0;
    int mpi_size = 0;
    
    int active_proc = 0;



    MPI_Comm world_comm;
    MPI_Comm cart_comm;
    int cart_rank;

    int domain_coords[2];



public:
    FinDiffSolution(int M, int N, int THREADS, MPI_Comm world_comm):M(M), N(N), THREADS(THREADS), X(M+1,N+1, 0), 
      world_comm(world_comm)
    {

        double h1 = (B1 - A1) / M;
        double h2 = (B2 - A2) / N;
    
        double hmax = std::max(h1, h2);
        mesh_eps = hmax * hmax; 


        MPI_Comm_size(world_comm, &mpi_size); 

        std::cout << "[!]MPI SIZE = " << mpi_size << std::endl;  

        omp_set_num_threads(THREADS);
    }

    ~FinDiffSolution(){
        
    }

    void SetMPISize(int size){
        mpi_size = size;
    }

    void SetMPIRank(int rank){
        mpi_rank = rank;
    }


    //check whether point belongs to area D
    bool D(const Point& p);

    double x(double i);
    double y(double j);

    //length of the segment part between points which belongs to D
    double l_ij(Point p1, Point p2);

    double a_ij(int i, int j);

    double b_ij(int i, int j);

    double F_ij(int i, int j);


    void exchange(Matrix<double>& exchangeable);
    Matrix<double> iteration_solver(SparseRows &A, Matrix<double> &B);

    void SpecifyDomain();
    bool Split(int m, int n, int n_proc);



    void show_Dij();
    void show_aij();
    void show_bij();
    void show_Fij();

    void show_Xij();
    void dump_Xij();

    void solve();

    bool approx(double x1, double x2, int cnt);


    int dom_i_min, dom_i_max, dom_j_min, dom_j_max;






    //_________________________________
    //Statistics of MPI usage
    int bytes_sent = 0;
    int total_bytes_sent = 0;

    void reset_statistics(){
      bytes_sent = 0;
    }

    void global_statistics(){
      MPI_Allreduce(&bytes_sent, &total_bytes_sent,   1, MPI_INT, MPI_SUM, cart_comm );
    }
    
};
