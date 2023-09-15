#pragma once

#include "SolverAlgorithm.hpp"

class Simplex : public SolverAlgorithm
{
public:
    std::vector<double> denseMatrix;
    double solve(unsigned n, unsigned stepLimit) override
    {
        if (n == 0)
            return 0.0;
        if (m == 0)
            return numeric_limits<double>::infinity();

        // initial feasibility is guaranteed since b is all ones
        denseMatrix = prepareDenseMatrix(matrix, m, n);

        double z = 0.0;
        unsigned step = 0;

        for (; step < stepLimit; step++)
        {
            vector<unsigned> enteringVars = findPivotColumnCandidates();

            if (enteringVars.size() == 0)
            {
                printf("\nNo entering var. Optimal value of %5.3f has been reached.\n", z);
                return z;
            }
            for (unsigned pivotColumn : enteringVars)
            {
                unsigned pivotRow = findPivotRowOrCheckUnboundedness(pivotColumn, n);
                if (pivotRow == -1)
                {
                    continue;
                }
                double pivotValue = denseMatrix[pivotRow * (n + m) + pivotColumn];
                if (pivotValue > epsilon2)
                {
                    doPivotting(pivotRow, pivotColumn, z, n);
                    break;
                }
            }
        }
        return z;
    }

    void doPivotting(int pivotRow, int pivotColumn, double &z, int n)
    {
        // Calculate the pivot value
        double pivotValue = denseMatrix[pivotRow * (n + m) + pivotColumn];

        // Update the objective function value
        z -= C[pivotColumn] * (B[pivotRow] / pivotValue);

        // Normalize the pivot row
        for (int i = 0; i < (n + m); ++i)
        {
            denseMatrix[pivotRow * (n + m) + i] /= pivotValue;
        }
        B[pivotRow] /= pivotValue;

        // Update the other rows in A and B
        for (int i = 0; i < m; ++i)
        {
            if (i != pivotRow)
            {
                double multiplier = denseMatrix[i * (n + m) + pivotColumn];
                for (int j = 0; j < (n + m); ++j)
                {
                    denseMatrix[i * (n + m) + j] -= multiplier * denseMatrix[pivotRow * (n + m) + j];
                }
                B[i] -= multiplier * B[pivotRow];
            }
        }

        // Update the objective function coefficients
        double multiplier = C[pivotColumn];
        for (int i = 0; i < (n + m); ++i)
        {
            C[i] -= multiplier * denseMatrix[pivotRow * (n + m) + i];
        }
    }

    std::vector<double> prepareDenseMatrix(const std::vector<Entry> &sparseMatrix, int m, int n)
    {
        // Initialize dense vector with zeros
        std::vector<double> denseVector(m * (m + n), 0.0);

        // Fill in the non-zero values from the sparse matrix
        for (const auto &entry : sparseMatrix)
        {
            int index = entry.rule * n + entry.variable;
            denseVector[index] = entry.coef;
        }
        // Add slack variables
        for (int i = 0; i < m; ++i)
        {
            int index = i * (n + m) + n + i;
            denseVector[index] = 1.0;
        }

        return denseVector;
    }
    vector<unsigned> findPivotColumnCandidates()
    {
        vector<unsigned> negativeIndices;
        for (int i = 0; i < C.size(); ++i)
        {
            if (C[i] < 0)
            {
                negativeIndices.push_back(i);
            }
        }
        return negativeIndices;
    }

    int findPivotRowOrCheckUnboundedness(int pivotColumn, unsigned n)
    {
        double value;
        int pivotRow = 0;
        int nonpositiveCount = 0;
        double min_denominator = 1.0; // Initialize to 1 to avoid division by zero
        double min_numerator = numeric_limits<double>::max();
        ;

        for (int i = 0; i < m; i++)
        {
            value = denseMatrix[i * (n + m) + pivotColumn]; // Adjusted for 1D representation
            double b_value = B[i];

            if (value > 0)
            {
                if (value * min_denominator < b_value * min_numerator)
                {
                    min_numerator = value;
                    min_denominator = b_value;
                    pivotRow = i;
                }
            }
            else
            {
                nonpositiveCount++;
            }
        }

        if (nonpositiveCount == m)
        {
            return -1;
        }

        return pivotRow;
    }
};