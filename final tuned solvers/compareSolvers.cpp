#include <iostream>
#include <vector>
#include "Simplex.hpp"
#include "RevisedSimplexMPFI.hpp"
#include "RevisedSimplexPFI.hpp"

int main()
{

    std::ifstream input("lp.txt");
    std::ofstream outfileTableau("benchmark_results.txt", std::ios::out);
    std::ofstream outfilePFI("benchmark_results_PFI.txt", std::ios::out);
    std::ofstream outfileMPFI("benchmark_results_MPFI.txt", std::ios::out);

    if (!input.is_open())
    {
        std::cerr << "Failed to open input file" << std::endl;
        return 1;
    }

    if (!outfileTableau.is_open())
    {
        std::cerr << "Failed to open output file" << std::endl;
        return 1;
    }
    if (!outfilePFI.is_open())
    {
        std::cerr << "Failed to open output file" << std::endl;
        return 1;
    }
    if (!outfileMPFI.is_open())
    {
        std::cerr << "Failed to open output file" << std::endl;
        return 1;
    }

    int trial = 0;
    const int MAX_REPS = 100; // Maximum number of repetitions
    int reps;                 // Actual number of repetitions

    for (string line; getline(input, line); ++trial)
    {

        tuple<vector<Entry>, unsigned, unsigned> myData = parseCoefficientMatrix(line);

        unsigned numRules = get<1>(myData);
        unsigned numVariables = get<2>(myData);
        vector<double> c;
        c.resize(numVariables + numRules);
        fill(c.begin(), c.begin() + numVariables, -1.0); // Fill the first n elements with -1
        fill(c.begin() + numVariables, c.end(), 0.0);    // Fill the remaining m elements with 0

        // initialize the b vector
        vector<double> b(numRules, 1.0);
        //---------------------------------------------------------------------------------------------
        Simplex simplexSolver;
        simplexSolver.matrix = get<0>(myData);
        simplexSolver.m = numRules;
        simplexSolver.B = b;
        simplexSolver.C = c;

        auto totalDuration = 0;
        double zTableau;
        for (reps = 0; reps < MAX_REPS; ++reps)
        {
            simplexSolver.B = b;
            simplexSolver.C = c;
            auto start = std::chrono::high_resolution_clock::now();
            zTableau = simplexSolver.solve(numVariables, ~0u);
            auto stop = std::chrono::high_resolution_clock::now();
            totalDuration += std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
        }
        auto avgDuration = (double)totalDuration / reps;
        //---------------------------------------------------------------------------------------------
        RevisedSimplexMPFI revisedSimplexMPFISolver;
        revisedSimplexMPFISolver.matrix = get<0>(myData);
        revisedSimplexMPFISolver.m = numRules;

        auto totalDurationMPFI = 0;
        double zMPFI;
        for (reps = 0; reps < MAX_REPS; ++reps)
        {
            auto startMPFI = std::chrono::high_resolution_clock::now();
            zMPFI = revisedSimplexMPFISolver.solveSimplex(numVariables, ~0u);
            auto stopMPFI = std::chrono::high_resolution_clock::now();
            totalDurationMPFI += std::chrono::duration_cast<std::chrono::microseconds>(stopMPFI - startMPFI).count();
        }

        auto avgDurationMPFI = (double)totalDurationMPFI / reps;
        //---------------------------------------------------------------------------------------------
        RevisedSimplexPFI revisedSimplexSolver;
        revisedSimplexSolver.matrix = get<0>(myData);
        revisedSimplexSolver.m = numRules;

        auto totalDurationPFI = 0;
        double zPFI;
        for (reps = 0; reps < MAX_REPS; ++reps)
        {
            auto startPFI = std::chrono::high_resolution_clock::now();
            zPFI = revisedSimplexSolver.solve(numVariables, ~0u);
            auto stopPFI = std::chrono::high_resolution_clock::now();
            totalDurationPFI += std::chrono::duration_cast<std::chrono::microseconds>(stopPFI - startPFI).count();
        }

        auto avgDurationPFI = (double)totalDurationPFI / reps;
        //---------------------------------------------------------------------------------------------

        //---------------------------------------------------------------------------------------------

        // Output the result to  file

        outfileTableau << "Trial " << trial + 1 << ": Optimal value of " << zTableau
                       << " has been reached. Time taken: " << avgDuration << " microseconds. NumIter: " << simplexSolver.numberStepsLastLP << " iterations.\n";
        outfilePFI << "Trial " << trial + 1 << ": Optimal value of " << zPFI
                   << " has been reached. Time taken: " << avgDurationPFI << " microseconds. NumIter: " << revisedSimplexSolver.numberStepsLastLP << " iterations.\n";
        outfileMPFI << "Trial " << trial + 1 << ": Optimal value of " << zMPFI
                    << " has been reached. Time taken: " << avgDurationMPFI << " microseconds. NumIter: " << revisedSimplexMPFISolver.lastStepCount << " iterations.\n";
    }

    input.close();
    outfileTableau.close();
    outfilePFI.close();
    outfileMPFI.close();

    return 0;
}