#include <iostream>
#include <vector>
#include "Simplex.hpp"
// #include "RevisedSimplexMPFI.hpp"
#include "RevisedSimplexPFI.hpp"

void compareSolvers(SolverAlgorithm &solver1, SolverAlgorithm &solver2, const std::vector<double> &coefficients);

int main()
{

    std::ifstream input("lp.txt");
    std::ofstream outfile("benchmark_results.txt", std::ios::out);

    if (!input.is_open())
    {
        std::cerr << "Failed to open input file" << std::endl;
        return 1;
    }

    if (!outfile.is_open())
    {
        std::cerr << "Failed to open output file" << std::endl;
        return 1;
    }

    int trial = 0;

    for (string line; getline(input, line); ++trial)
    {

        tuple<vector<Entry>, unsigned, unsigned> myData = parseCoefficientMatrix(line);

        Simplex simplexSolver;
        // initialize the coefficient matrix sparse
        simplexSolver.matrix = get<0>(myData);
        unsigned numRules = get<1>(myData);
        unsigned numVariables = get<2>(myData);
        simplexSolver.m = numRules;
        // initialize the b vector
        vector<double> b(numRules, 1.0);
        simplexSolver.B = b;
        // initialize the c vector
        vector<double> c;
        c.resize(numVariables + numRules); 
        fill(c.begin(), c.begin() + numVariables, -1.0); // Fill the first n elements with -1
        fill(c.begin() + numVariables, c.end(), 0.0);    // Fill the remaining m elements with 0
        simplexSolver.C = c;

        auto start = std::chrono::high_resolution_clock::now();

        auto z = simplexSolver.solve(numVariables, 30);

        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

        // Output the result to console and file
        printf("Trial %d: Optimal value of %5.3f has been reached. Time taken: %d microseconds\n", trial + 1, z, duration.count());
        outfile << "Trial " << trial + 1 << ": Optimal value of " << z << " has been reached. Time taken: " << duration.count() << " microseconds.\n";
    }

    input.close();
    outfile.close();

    return 0;
}