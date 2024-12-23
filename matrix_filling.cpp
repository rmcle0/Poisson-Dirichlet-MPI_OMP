#include "matrix_filling.h"



bool FinDiffSolution::D(const Point& p){
    if( abs(p.x) < 1 && abs(p.y) < 1){
        if( 0 < p.x && 0 < p.y){
            return false;
        }
        return true;
    }

    return false;
}

double FinDiffSolution::x(double i){
    if(i < 0 || i > M)
        throw std::runtime_error("x: invalid index");

    return A1 + h1 * i;
}

double FinDiffSolution::y(double j){
    if(j < 0 || j > N)
        throw std::runtime_error("y: invalid index");

    return A2 + h2 * j;
}

// Searching for the length inside D between 2 points
double FinDiffSolution::l_ij(Point p1, Point p2)
{
    bool swap_flag = false;

    if (D(p1) && D(p2))
        return (p1 - p2).length();
    else if (!D(p1) && !D(p2))
        return 0;
    else // search procedure
    {
        Point pi1 = p1;
        Point pi2 = p2;

        // For convinience, let p1 be that D(p1) be true
        if (D(p2))
        {
            swap(p1, p2);
            swap_flag = true;
        }

        Point pb = (p1 + p2) / 2; // border point where contour border lies

        while (true)
        {
            pb = (p1 + p2) / 2;

            if (D(pb))
                p1 = pb;
            else
                p2 = pb;

            if ((pb - p1).manhattan_length() < EPS &&
                (pb - p2).manhattan_length() < EPS)
                break;
        }

        if (!swap_flag) // p1 was initially in D
            return (pb - pi1).length();
        else
            return (pb - pi2).length();
    }

    return 0;
}



//from y_(j-1/2) to y_(j+1/2), x is fixed = i-1/2
double FinDiffSolution::a_ij(int i, int j){
    double y_jm05 = y((double)j - 0.5);
    double y_jp05 = y((double)j + 0.5);
    double x_im05 = x((double)i - 0.5);

    Point p1(x_im05, y_jm05);
    Point p2(x_im05, y_jp05);

    double l = l_ij(p1, p2);

    return l / h2  + (1.0 - l / h2 ) / mesh_eps;
}


double FinDiffSolution::b_ij(int i, int j){
    double x_im05 = x((double)i - 0.5);
    double x_ip05 = x((double)i + 0.5);

    double y_jm05 = y((double)j - 0.5);



    Point p1(x_im05, y_jm05);
    Point p2(x_ip05, y_jm05);

    double l = l_ij(p1, p2);

    return l / h1  + (1.0 - l / h1 ) / mesh_eps;
}




void FinDiffSolution::show_Dij(){
    for(int j = N; j >= 0; j--){
    for(int i = 0; i <= M; i++){
           printf("%2d  ", D( Point(A1 + i * h1 , A2 + j * h2)) );
        }
    //cout << endl;
   
    printf("\n");
    }
}


void FinDiffSolution::show_aij(){

    for(int j = N - 1; j > 0; j--){
    for(int i = 1; i < M; i++){
           printf("%10.6f  ", a_ij(i, j));
        }   
    printf("\n");
    }
}

void FinDiffSolution::show_bij(){
    for(int j = N - 1; j > 0; j--){
    for(int i = 1; i < M; i++){
           printf("%10.6f  ", b_ij(i, j));
        }   
    printf("\n");
    }
}

void FinDiffSolution::show_Fij(){
    for(int j = N - 1; j > 0; j--){
    for(int i = 1; i < M; i++){
           printf("%10.6f  ", F_ij(i, j));
        }   
    printf("\n");
    }
}



double FinDiffSolution::F_ij(int i, int j){
    double x_im05 = x((double)i - 0.5);
    double x_ip05 = x((double)i + 0.5);

    double y_jm05 = y((double)j - 0.5);
    double y_jp05 = y((double)j + 0.5);


    std::deque<Point> rc;
    rc.push_back({x_im05, y_jm05});  
    rc.push_back({x_ip05, y_jm05});
    rc.push_back({x_ip05, y_jp05});
    rc.push_back({x_im05, y_jp05});


    if(D(rc[0]) && D(rc[1]) && D(rc[2]) && D(rc[3])){
        //return h1 * h2; ERROR!

        return 1; 
    }
    else if(!D(rc[0]) && !D(rc[1]) && !D(rc[2]) && !D(rc[3])){
        return 0;
    }

    
    while( ! (D(rc[0]) && !D(rc[1])) ){
        rc.push_back(rc[0]);
        rc.pop_front();
    }

    double a = (rc[0] - rc[1]).length();
    double b = (rc[1] - rc[2]).length();

    double S = 0;

    //now 0 is in D, 1 is not. Only 3 options  are left
    if(D(rc[3])){
        if(D(rc[2])){
            S = a * b - 
            0.5 * (a - l_ij(rc[0], rc[1]) ) * (b - l_ij(rc[1], rc[2]) );
        }
        else 
            {
            S = a * b -  0.5 * (2 * a - l_ij(rc[0], rc[1]) - l_ij(rc[2], rc[3]) )   * b;
            }
    }
    else {
        S = /*a * b -*/ //ERROR 
            0.5 *  l_ij(rc[0], rc[3]) * l_ij(rc[0], rc[1]) ;
    }


    return S / (h1 * h2) * 1; //* f(rc[0]);
    
}






void FinDiffSolution::exchange(Matrix<double>& w_k){

    //if(mpi_size == 1)
    //    return;

    static MPI_Request r[8];
    static MPI_Status s[8];

    static vector<double> w_up(current_domain_size_x, 0), w_down(current_domain_size_x, 0), w_left(current_domain_size_y, 0), w_right(current_domain_size_y, 0);
    static vector<double> w_up2(current_domain_size_x, 0), w_down2(current_domain_size_x, 0), w_left2(current_domain_size_y, 0), w_right2(current_domain_size_y, 0);



    for (int j = 0; j < current_domain_size_y; j++) {
      w_left[j] = w_k.At_WithShadow(0, j);
    }

    for (int j = 0; j < current_domain_size_y; j++) {
      w_right[j] = w_k.At_WithShadow(current_domain_size_x - 1, j);
    }

    for (int i = 0; i < current_domain_size_x; i++) {
      w_down[i] = w_k.At_WithShadow(i, 0);
    }

    for (int i = 0; i < current_domain_size_x; i++) {
      w_up[i] = w_k.At_WithShadow(i, current_domain_size_y - 1);
    }




    int source, destination;
    MPI_Cart_shift(cart_comm, 1, 1,  &source, &destination );
        MPI_Irecv( w_down2.data(), current_domain_size_x, MPI_DOUBLE, source, 1, cart_comm, &r[0]);
        MPI_Isend( w_up.data(), current_domain_size_x, MPI_DOUBLE, destination, 1, cart_comm, &r[1]);

     bytes_sent += current_domain_size_x * sizeof(double);

    MPI_Cart_shift(cart_comm, 0, 1,  &source, &destination );
        MPI_Irecv( w_left2.data(), current_domain_size_y, MPI_DOUBLE, source, 2, cart_comm, &r[2]);
        MPI_Isend( w_right.data(), current_domain_size_y, MPI_DOUBLE, destination, 2, cart_comm, &r[3]);

     bytes_sent += current_domain_size_y * sizeof(double);

    MPI_Cart_shift(cart_comm, 1, -1,  &source, &destination );
        MPI_Irecv( w_up2.data(), current_domain_size_x, MPI_DOUBLE, source, 3, cart_comm, &r[4]);
        MPI_Isend( w_down.data(), current_domain_size_x, MPI_DOUBLE, destination, 3, cart_comm, &r[5]);

     bytes_sent += current_domain_size_x * sizeof(double);

    MPI_Cart_shift(cart_comm, 0, -1,  &source, &destination );
        MPI_Irecv( w_right2.data(), current_domain_size_y, MPI_DOUBLE, source, 4, cart_comm, &r[6]);
        MPI_Isend( w_left.data(), current_domain_size_y, MPI_DOUBLE, destination, 4, cart_comm, &r[7]);


     bytes_sent += current_domain_size_y * sizeof(double);
     

        MPI_Waitall(8, r, s);

        MPI_Barrier(cart_comm);



            for (int j = 0; j < current_domain_size_y; j++) 
                w_k.At_WithShadow(-1, j) = w_left2[j];
        

            for (int j = 0; j < current_domain_size_y; j++) 
                w_k.At_WithShadow(current_domain_size_x, j) = w_right2[j];
        
            for (int i = 0; i < current_domain_size_x; i++)
                w_k.At_WithShadow(i,  -1) = w_down2[i];
        

            for (int i = 0; i < current_domain_size_x; i++)
                w_k.At_WithShadow(i, current_domain_size_y) = w_up2[i];
        



}



Matrix<double> FinDiffSolution::iteration_solver(SparseRows &A, Matrix<double> &B)
{
    Matrix<double> w_k(current_domain_size_x + 2 , current_domain_size_y + 2, 0);  //current_domain_size_x
    Matrix<double> w_kp1(current_domain_size_x + 2 , current_domain_size_y + 2, 0);

    Matrix<double> Aw_k(current_domain_size_x + 2 , current_domain_size_y + 2, 0);

    Matrix<double> r_k(current_domain_size_x + 2 , current_domain_size_y + 2, 0);
    Matrix<double> Ar_k(current_domain_size_x + 2 , current_domain_size_y + 2, 0);



    

    double max_difference = 0;
    double global_max_difference = __INT16_MAX__;
    int iteration = 0;


       std::chrono::steady_clock::time_point begin =
        std::chrono::steady_clock::now();


    
    while ( global_max_difference > DELTA && iteration < ITMAX )
    {
        MPI_Barrier(cart_comm);

        //Sent bytes set to zero
        reset_statistics();

        
        //checking for deadlocks
        if(cart_rank == 0 && iteration % 100 == 0){
            //std::cout << "it = " << iteration << " " <<  global_max_difference << std::endl;
            std::chrono::steady_clock::time_point end =
                std::chrono::steady_clock::now();
            std::cout << "it " << iteration << " Time difference = "  << " "
                    << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                            begin)
                            .count()
                    << "[ms]" << std::endl;


        }


        // calculating A * w_k (w_k is known for all processes)
        // for into domain only
        #pragma omp parallel for
        for (int i = 0; i < current_domain_size_x; i++)
          for (int j = 0; j < current_domain_size_y; j++) {
            Aw_k.At_WithShadow(i, j) = A.multiply(i, j, w_k);

            r_k.At_WithShadow(i, j) = Aw_k.At_WithShadow(i, j) - B.At(i, j);
          }

        //Exchange r_k
        exchange(r_k);
            

        


        // Calculate Tau_k+1
        double numerator = 0, denominator = 0;

        #pragma omp parallel for reduction (+:denominator, numerator)
        for (int i = 0; i < current_domain_size_x; i++)
          for (int j = 0; j < current_domain_size_y; j++) {
            Ar_k.At_WithShadow(i, j) = A.multiply(i, j, r_k);

            numerator += r_k.At_WithShadow(i, j) * r_k.At_WithShadow(i, j);
            denominator += Ar_k.At_WithShadow(i, j) * r_k.At_WithShadow(i, j);
          }

        

        double global_numerator, global_denominator;

        MPI_Allreduce(&numerator, &global_numerator,  1, MPI_DOUBLE, MPI_SUM, cart_comm );
        MPI_Allreduce(&denominator, &global_denominator,  1, MPI_DOUBLE, MPI_SUM, cart_comm );


        double tau = global_numerator / global_denominator;

        // find max:
        max_difference = 0;

        #pragma omp parallel for reduction(max:max_difference)
        for (int i = 0; i < current_domain_size_x; i++)
          for (int j = 0; j < current_domain_size_y; j++)
            {  
            //w_kp1.At_WithShadow(i,j) = w_k.At_WithShadow(i,j) - tau * r_k.At_WithShadow(i,j);
            w_k.At_WithShadow(i,j) -= tau * r_k.At_WithShadow(i,j);

            max_difference = std::max( std::abs( tau * r_k.At_WithShadow(i,j)) , max_difference);

            //max_difference = std::max( std::abs( w_kp1.At_WithShadow(i,j) - w_k.At_WithShadow(i,j)) , max_difference);
            }

        MPI_Allreduce(&max_difference, &global_max_difference,   1, MPI_DOUBLE, MPI_MAX, cart_comm );

        //swap(w_k, w_kp1);

        exchange(w_k);

        global_statistics();

        iteration++;
    }


    std::cout << "bytes sent " << bytes_sent << std::endl;

    if(cart_rank == 0){
    std::cout << "res_eps = " << global_max_difference << std::endl;
    std::cout << "iteration = " <<iteration << std::endl;
    

    std::cout << "last iteration sends " << total_bytes_sent << std::endl;
    }

    return w_k; 
}









/* Fills class data members in an attempt to split the mesh */
bool FinDiffSolution::Split(int m, int n, int n_proc){
    int& a = base_domain_size_x;
    int& b = base_domain_size_y;
    

    for(a = 1; a <= m; a++){
		for(b = 1; b <= n; b++){
			//condition on  a,b size relations
			if( ((double)std::max(a,b)) / ((double)std::min(a,b)) > (2 + 0.001) )
				continue;



            int domains_along_x = m / a;
            int domains_along_y = n / b;
			
			active_proc =  domains_along_x * domains_along_y;
			
			if(  active_proc < ((double)n_proc) * 0.9 || active_proc > n_proc )
				continue;	

            //trying to distribute the remnants between the domains
            //let the first some domains be bigger than those with bigger coordinates
            
            // Example where impossible to distribute: 39, size 10, module 10, quantity of domains is 3
            if(m % a > domains_along_x || n % b > domains_along_y)
                continue;

            //first m%a domains are bigger
            //smaller begin with exactly m%a index
            ind_of_smaller_domains_X = (m % a);
            ind_of_smaller_domains_Y = (n % b);

            if (verbose) {
              std::cout << "ind of smallerX " << ind_of_smaller_domains_X
                        << std::endl;
              std::cout << "ind of smallerY " << ind_of_smaller_domains_Y
                        << std::endl;
            }

            domains_quantity_along_x = (m / a);
            domains_quantity_along_y = (n / b);


			return true;
		}
	}
	
	return false;
}



void FinDiffSolution::SpecifyDomain(){
    current_domain_size_x = base_domain_size_x;
    current_domain_size_y = base_domain_size_y;

    if(domain_coords[0] <  ind_of_smaller_domains_X)
        current_domain_size_x++;

    if(domain_coords[1] < ind_of_smaller_domains_Y)
        current_domain_size_y++;

    current_domain_size = current_domain_size_x * current_domain_size_y;




    int large_domains_x = min(domain_coords[0],   ind_of_smaller_domains_X);     
    int large_domains_y = min(domain_coords[1],   ind_of_smaller_domains_Y);

    dom_i_min = max(large_domains_x * (base_domain_size_x + 1) + 
                    (domain_coords[0] - large_domains_x) * (base_domain_size_x )  + 1, 1);
    dom_i_max =  min( dom_i_min + current_domain_size_x , M);
    
    dom_j_min = max(large_domains_y * (base_domain_size_y + 1) + 
                    (domain_coords[1] - large_domains_y) * (base_domain_size_y ) + 1, 1);
    dom_j_max =  min( dom_j_min + current_domain_size_y, N);
}





void FinDiffSolution::solve()
{
    //Create the partitioning of П internal equations to the domains
    // each domain corresponds to a set of equations for points in this domain
    if(!Split(M - 1, N - 1, mpi_size)){
        std::cerr << "impossible to split " << std::endl;
        return;
    }

    std::cout << "here" << std::endl;

    printf(" Was trying to split %d %d  to %d processes \n", M - 1, N - 1, mpi_size);
    printf(" domain base size %d %d \n", base_domain_size_x, base_domain_size_y);
    printf(" active_proc %d \n", active_proc);



    //Creating a cart communicator
    int dim[2] = {domains_quantity_along_x, domains_quantity_along_y};
    int period[2]={false, false};
    int reorder = true; 
    
    
    std::cout << "trying to split to " << domains_quantity_along_x << " " << domains_quantity_along_y << std::endl;

    MPI_Cart_create(world_comm, 2, dim, period,  reorder, &cart_comm);

    if(cart_comm == MPI_COMM_NULL)
        return;

    MPI_Cart_coords(cart_comm, mpi_rank, 2,  domain_coords);
    MPI_Cart_rank(cart_comm, domain_coords, &cart_rank);


    SpecifyDomain();




    // Creating right side Fij, and left side of koefficients for wij
    Matrix<double> Fij((current_domain_size_x), (current_domain_size_y), 0);

    SparseRows A(current_domain_size, current_domain_size_x, current_domain_size_y);


    // writing down equations from the linear system
    for( int i = dom_i_min; i < dom_i_max; i++ )
        for(int j = dom_j_min; j <  dom_j_max; j++ )
        {
            int iw = i - dom_i_min;
            int jw = j - dom_j_min;

            A.set(iw, jw, SparseRows::Center, 1.0 / h1 * a_ij(i + 1, j) / h1 + 
            1.0 / h1 * a_ij(i, j) / h1  + 
            1.0 / h2 * b_ij(i, j + 1) / h2 + 
            1.0 / h2 * b_ij(i, j) / h2 );


            if(i + 1 < M)
                 A.set(iw, jw, SparseRows::Right, -1.0 / h1 * a_ij(i + 1, j) / h1);
            
            if(i - 1 > 0)
                A.set(iw, jw, SparseRows::Left, -1.0 / h1 * a_ij(i, j) / h1);

            if(j + 1 < N)
                A.set(iw, jw, SparseRows::Up, -1.0 / h2 * b_ij(i, j + 1) / h2);

            if(j - 1 > 0)
                A.set(iw, jw, SparseRows::Down, -1.0 / h2 * b_ij(i, j) / h2);
        }


    // filling right side
    for( int i = dom_i_min; i < dom_i_max; i++ )
        for(int j = dom_j_min; j <  dom_j_max; j++ )
        {
        Fij.At(i - dom_i_min, j - dom_j_min) = F_ij(i, j);
        }



    // TRY TO CALL Solver.
    X = iteration_solver(A, Fij);

    std::cout << "here" << std::endl;

    //dump_Xij();

    return;
    



    //Forming equal-sized send arrays
    /*double *sendarray = new double[domain_size_x * domain_size_y];

    for(int i = 0; i < domain_size_x; i++){
        for(int j = 0; j < domain_size_y; j++){
            try{
            sendarray[j * domain_size_x + i] = X.At(domain_coord_along_x * domain_size_x + i,
                                                    domain_coord_along_y * domain_size_y + j );
            }
            catch(...){
                sendarray[j * domain_size_x + i] = 0;
            }
        }
    }

    double *recv = nullptr;
    if(mpi_rank == 0)
        recv = new double[active_proc * domain_size_x * domain_size_y];


    MPI_Gather(sendarray, domain_size_x * domain_size_y, MPI_DOUBLE, recv, domain_size_x * domain_size_y, MPI_DOUBLE, 0, cart_comm);
    

    if(mpi_rank != 0)
        return;


    




   
    if(!mpi_rank){
        for(int i = 0; i < M + 1; i++)
            for(int j = 0; j < N + 1; j++){
                int coords[2] = {i / domain_size_x, j / domain_size_y};
                int rank;

                MPI_Cart_rank(cart_comm, coords, &rank);

                int in_dom_x = i - coords[0] * domain_size_x;
                int in_dom_y = j - coords[1] * domain_size_y;

                int offset = rank * domain_size_x * domain_size_y;

                X.At(i, j) = recv[offset + in_dom_y * domain_size_x + in_dom_x];
            }
    }
*/


    

    /*if(mpi_rank != 0)
        return;*/

    

    /*for (int i = 1; i <= M - 1; i++)
        for (int j = 1; j <= N - 1; j++){
            SparseRow &eq = A.At(i,j);

            for (int i1 = 1; i1 <= M - 1; i1++)
                for (int j1 = 1; j1 <= N - 1; j1++){
                    if(eq.find({i1,j1}) != eq.end()){
                        std::cout << eq[{i1,j1}] << " ";
                    }
                    else 
                        cout << "0" << " ";

                }
            
            cout << std::endl;
        }*/
}

bool FinDiffSolution::approx(double x1, double x2, int cnt){
    return abs(x1 - x2) < 0.05;
}



void FinDiffSolution::show_Xij(){
    for(int j = current_domain_size_y; j >= 1; j--){
    for(int i = 1; i < current_domain_size_x; i++){
           printf("%10.5f  ", X.At(i, j));
        }   
    printf("\n");
    }
}

void FinDiffSolution::dump_Xij()
{
    std::ofstream ofs;

    #ifdef _OPENMP
    std::string filename = std::string("TABLE_1/H_") + std::to_string(M) + "_" + std::to_string(N) + "_"
                         + std::to_string(domain_coords[0]) + "_" + std::to_string(domain_coords[1]) + "_" + std::to_string(THREADS);
    #else 
    std::string filename = std::string("TABLE_1/CONS_") + std::to_string(M) + "_" + std::to_string(N) + "_" + 
        std::to_string(THREADS) + "_" + std::to_string(mpi_rank);
    #endif 

    ofs.open(filename, std::ofstream::out | std::ofstream::trunc);

    for(int j = current_domain_size_y; j >= 1; j--){
    for(int i = 1; i <= current_domain_size_x; i++){
            ofs << X.At(i, j) << " ";
            }   
        ofs << endl;
        }

    ofs.close();
}