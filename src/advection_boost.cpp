#include <iostream>
#include <math.h> // sqrt, exp, pow
#include <string> // stoi, stod
#include <fstream> // file writing
#include <chrono> // std::chrono


#include <omp.h>


void writeToFile(std::string const& fileName, double** const& result, int const& NT, int const& N) {
    std::ofstream myfile;

    myfile.open(fileName);

    myfile << "[" << NT << "], ";

    myfile << "[";
    for (int i = 0; i < N; i++) {
        if (i == 0) {
            myfile << "[";
        } else {
            myfile << ", [";
        }
        for (int j = 0; j < N; j++) {
            if (j == N - 1) {
                myfile << result[i][j];
            } else {
                myfile << result[i][j] << ", ";
            }
        }
        myfile << "]";
    }
    myfile << "]";
    myfile.close();
}


double** createMatrix(int const& N) {
    double** matrix = new double*[N];

    for (int i = 0; i < N; i++) {
        matrix[i] = new double[N];

        for (int j = 0; j < N; j++) {
            matrix[i][j] = 0;
        }
    }

    return matrix;
}

void applyBoundaryConditions(double** C_n, double** C_nplusone, int const& N, double const& delta_x, double const& delta_t, double const& u, double const& v) {
    int i = 0, j = 0;
    double up;
    double right;
    double down;
    double left;
    #pragma omp parallel for default(none) private(i, j, up, down, left, right) shared(C_n, C_nplusone, N, delta_x, delta_t, u, v) schedule(guided) 
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            if (i != 0 && i != N-1) {
                up = C_n[i-1][j];
                down = C_n[i+1][j];
                if (j != 0 && j != N-1) {
                    left = C_n[i][j-1];
                    right = C_n[i][j+1];
                }

                if (j == 0) {
                    left = C_n[i][N-1];
                    right = C_n[i][j+1];
                } else if (j == N-1) {
                    left = C_n[i][j-1];
                    right = C_n[i][0];
                }
            }

            // apply periodic boundary conditions 
            if (i == 0) {
                up = C_n[N-1][j];
                down = C_n[i+1][j];

                if (j != 0 && j != N-1) {
                    left = C_n[i][j-1];
                    right = C_n[i][j+1];
                }

                if (j == 0) {
                    left = C_n[i][N-1];
                    right = C_n[i][j+1];
                } else if (j == N-1) {
                    right = C_n[i][0];
                    left = C_n[i][j-1];
                }
            } else if (i == N-1) {
                up = C_n[i-1][j];
                down = C_n[0][j];

                if (j != 0 && j != N-1) {
                    left = C_n[i][j-1];
                    right = C_n[i][j+1];
                }

                if (j == 0) {
                    left = C_n[i][N-1];
                    right = C_n[i][j+1];
                } else if (j == N-1) {
                    right = C_n[i][0];
                    left = C_n[i][j-1]; 
                }
            }

            // update result vector
            C_nplusone[i][j] = 0.25*(up + down + left + right) - ((delta_t/(2*delta_x)) * (u*(down - up) + v*(right - left)));
        }
    }
}

void advection(int const& N, int const& NT, double const& L, double const& T, double const& u, double const& v) {
    double** C_n = createMatrix(N);
    double** C_nplusone = createMatrix(N);

    double delta_x = L/N;
    double delta_t = T/NT;

    if (delta_t > delta_x / sqrt(2*(u*u + v*v))) {
        // check for Courant stability condition
        return;
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            // initialise as Gaussian
            C_n[i][j] = exp (-(pow(((i*L/N)-(L/2.0)), 2)/(2*pow(L/4.0, 2)) + pow(((j*L/N)-(L/2.0)), 2)/(2*pow(L/4.0, 2))));
        }
    }
    
    applyBoundaryConditions(C_n, C_nplusone, N, delta_x, delta_t, u, v);

    // copy forecast to current grid
    std::swap(C_n, C_nplusone);

    writeToFile("../out/out_0.txt", C_n, 0, N);

    for (int n = 0; n < NT; n++) {
        std::cout << "n: " << n << std::endl;
	double t1 = omp_get_wtime();
        
        applyBoundaryConditions(C_n, C_nplusone, N, delta_x, delta_t, u, v);

        // copy forecast to current grid
        std::swap(C_n, C_nplusone);

	double t2 = omp_get_wtime();
        std::cout << "iteration " << n << " (s) -> " << t2-t1 << std::endl;

        if (n == (NT/2) - 1) {
            writeToFile("../out/out_NTby2.txt", C_n, NT/2, N);
        } else if (n == NT - 1) {
            writeToFile("../out/out_NT.txt", C_n, NT, N);
        }
    }
}





int main(int argc, char** argv) {
    if (argc != 7) {
        std::cout << "Invalid Command Line Arguments. Please check and try again." << std::endl;
        std::cout << "Enter simulation arguments in this format: ./advection N NT L T u v" << std::endl;
        return 0;
    }

    int N;
    int NT;
    double L;
    double T;
    double u;
    double v;
    int numThreads = 12;
    omp_set_num_threads(numThreads);
    try
    {
        N = std::stoi(argv[1]);
        NT = std::stoi(argv[2]);
        L = std::stod(argv[3]);
        T = std::stod(argv[4]);
        u = std::stod(argv[5]);
        v = std::stod(argv[6]);
    }
    catch (std::exception const& e)
    {
        std::cout << "Invalid Command Line Arguments. Please check and try again." << std::endl;
        std::cout << "Error was caused by ~ " << e.what() << std::endl;
        return 0;
    }

    // Print parameters and estimated memory usage
    std::cout << "Estimated memory usage: " << N*N*sizeof(double) << " bytes" << std::endl;
    std::cout << "Simulating with parameters..." << std::endl;
    std::cout << "N: " << N << std::endl;
    std::cout << "NT: " << NT << std::endl;
    std::cout << "L: " << L << std::endl;
    std::cout << "T: " << T << std::endl;
    std::cout << "u: " << u << std::endl;
    std::cout << "v: " << v << std::endl;

    // Start simulation
    double t1 = omp_get_wtime();
    advection(N, NT, L, T, u, v);
    double t2 = omp_get_wtime();

    std::cout << "Simulation complete. Check out_*.txt files for output." << std::endl;
    std::cout << "time(s): " << t2-t1 << std::endl;
}
