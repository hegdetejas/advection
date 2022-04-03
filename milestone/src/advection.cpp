#include <iostream>
#include <vector>
#include <math.h> // sqrt, exp, pow
#include <string> // stoi, stod
#include <fstream> // file writing
#include <chrono> // std::chrono


void writeToFile(std::string const& fileName, std::vector<std::vector<double>> const& result, int const& NT) {
    std::ofstream myfile;

    myfile.open(fileName);

    myfile << "[" << NT << "], ";

    myfile << "[";
    for (int i = 0; i < result.size(); i++) {
        if (i == 0) {
            myfile << "[";
        } else {
            myfile << ", [";
        }
        for (int j = 0; j < result[i].size(); j++) {
            if (j == result[i].size() - 1) {
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



std::vector<std::vector<double>> advection(int const& N, int const& NT, double const& L, double const& T, double const& u, double const& v) {
    std::vector<std::vector<double>> C_n(N, std::vector<double>(N, 0.0)), C_nplusone(N, std::vector<double>(N, 0.0));

    double delta_x = L/N;
    double delta_t = T/NT;

    if (delta_t > delta_x / sqrt(2*(u*u + v*v))) {
        // check for Courant stability condition
        return std::vector<std::vector<double>>{};
    }

    for (int i = 0; i < C_n.size(); i++) {
        for (int j = 0; j < C_n[i].size(); j++) {
            // initialise as Gaussian
            C_n[i][j] = exp (-(pow((i-(L/2)), 2)/(2*pow(L/4, 2)) + pow((j-(L/2)), 2)/(2*pow(L/4, 2))));
        }
    }

    writeToFile("../out/out_0.txt", C_n, 0);

    for (int n = 0; n < NT; n++) {
        for (double i = 0; i < N; i++) {
            for (double j = 0; j < N; j++) {
                double left;
                double right;
                double up;
                double down;

                // apply periodic boundary conditions 
                if (i == 0) {
                    up = C_n[N-1][j];
                    down = C_n[i+1][j];
                    if (j == 0) {
                        left = C_n[i][N-1];
                        right = C_n[i][j+1];
                    } else if (j == N-1) {
                        right = C_n[i][0];
                        left = C_n[i][j-1];
                    } else {
                        left = C_n[i][j-1];
                        right = C_n[i][j+1];
                    }
                } else if (i == N-1) {
                    up = C_n[i-1][j];
                    down = C_n[0][j];
                    if (j == 0) {
                        left = C_n[i][N-1];
                        right = C_n[i][j+1];
                    } else if (j == N-1) {
                        right = C_n[i][0];
                        left = C_n[i][j-1]; 
                    } else {
                        left = C_n[i][j-1];
                        right = C_n[i][j+1];
                    }
                } else {
                    up = C_n[i-1][j];
                    down = C_n[i+1][j];
                    if (j == 0) {
                        left = C_n[i][N-1];
                        right = C_n[i][j+1];
                    } else if (j == N-1) {
                        left = C_n[i][j-1];
                        right = C_n[i][0];
                    } else {
                        left = C_n[i][j-1];
                        right = C_n[i][j+1];
                    }
                }

                // update result vector
                C_nplusone[i][j] = 0.25*(up + down + left + right) - ((delta_t/(2*delta_x)) * (u*(down - up) + v*(right - left)));
            }
        }
        // copy forecast to current grid
        std::swap(C_n, C_nplusone);

        if (n == (NT/2) - 1) {
            writeToFile("../out/out_NTby2.txt", C_n, NT/2);
        } else if (n == NT - 1) {
            writeToFile("../out/out_NT.txt", C_n, NT);
        }
    }
    return C_n;
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
    auto start = std::chrono::steady_clock::now();
    advection(N, NT, L, T, u, v);
    auto end = std::chrono::steady_clock::now();

    std::cout << "Simulation complete. Check out_*.txt files for output." << std::endl;
    std::cout << "Simulation took " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;
}