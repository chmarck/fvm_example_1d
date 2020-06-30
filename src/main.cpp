#include <iostream>
#include <chrono>
#include <eigen3/Eigen/Dense>

/*
Refer to:
An Introduction to Computational Fluid Dynamics THE FINITE VOLUME METHOD
H.K. Versteeg and W. Malalasekera
Example 4.2 pg. 122
http://ftp.demec.ufpr.br/disciplinas/TM702/Versteeg_Malalasekera_2ed.pdf

*/


void build_matrix(Eigen::MatrixXd& M, Eigen::VectorXd& b, unsigned int N, double dx, double k, double A, double q, double T_A, double T_B) {
     //Interior nodes
     double Su = q*A*dx;
     double aW = k*A/dx;
     double aE = k*A/dx;
     double aP = aW + aE;
     
     for (unsigned int i=1; i<N-1; i++) {
          M(i, i-1) = -aW;
          M(i, i) = aP;
          M(i, i+1) = -aE;

          b(i) = Su;
     }

     //Left boundary
     aE = k*A/dx;
     aP = aE + 2*k*A/dx;
     
     M(0,0) = aP;
     M(0,1) = -aE;

     b(0) = q*A*dx + 2*k*A*T_A/dx;

     //Right boundary
     aW = k*A/dx;
     aP = aW +2*k*A/dx;

     M(N-1,N-1) = aP;
     M(N-1,N-2) = -aW;

     b(N-1) = q*A*dx + 2*k*A*T_B/dx;   
}

void solve(Eigen::MatrixXd& M, Eigen::VectorXd& b, Eigen::VectorXd& result) {
     std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
     result = M.colPivHouseholderQr().solve(b);
     std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
     std::chrono::duration<double, std::nano> solver_time = end - start;

     std::cout << "Time to solve: "
          << std::chrono::duration_cast<std::chrono::microseconds>(solver_time).count()
          << "s" << std::endl;
}

int main(int argc, char* argv[]){
     unsigned int N = 2E3;//Number of nodes
     double L = 0.02;    //Domain length, m
     double q = 5E6;     //Heat generation, W/m^3
     double k = 0.5;     //Thermal conductivity, W/m.K
     double dx = L / N;  //Step size
     double A = dx;

     double T_A = 100; //Left boundary temperature
     double T_B = 500; //Right boundary temperature

     Eigen::MatrixXd M = Eigen::MatrixXd::Zero(N,N);
     Eigen::VectorXd b = Eigen::VectorXd::Zero(N);
     Eigen::VectorXd result = Eigen::VectorXd::Zero(N);

     build_matrix(M, b, N, dx, k, A, q, T_A, T_B);

     solve(M,b, result);
     
}