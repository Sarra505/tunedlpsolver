#pragma once

#include "SolverAlgorithm.hpp"

class RevisedSimplexPFI : public SolverAlgorithm
{
public:
    // a structure that represents an eta matrix
    // containing the column that's different from an identity matrix (ranges from 1 to m - 1)
    // and the values in that column
    struct eta
    {
        size_t col;
        vector<double> values;
    };
    double solve(unsigned n, unsigned stepLimit) override
    // Solve the packing LP using simplex
    {
        if (n == 0)
            return 0.0;
        if (m == 0)
            return numeric_limits<double>::infinity();

        prepareCCSFormat(n);

        nonBasic.resize(n);
        basic.resize(m);

        // Initialize the nonbasic variables' labels to be 0 through (n - 1)
        for (size_t i = 0; i < n; ++i)
        {
            nonBasic[i] = i;
        }
        // Initialize the basic variables to be n through m+n-1 and their values to be the right-hand-side vector b which is 1
        // initial feasible solution is (1,....,1) with z=0
        for (size_t i = 0; i < m; ++i)
        {
            basic[i].label = i + n;
            basic[i].value = 1.0;
        }

        // problem is initially feasible since b >= 0 for our dataset

        unsigned step = 0;
        bool foundSolution = false;

        /***************variables needed during each iteration****************/
        // An array of eta matrices representing previous pivots
        // the most recent pivot should be at the beginning
        vector<eta> pivots{};
        vector<pair<variable, unsigned>> enteringVarCandidates;

        // initial value of objective function
        double z = 0.0;
        /********************************************************************************/

        for (; step < stepLimit; step++)
        {
            // 1. find the entering variable:
            enteringVarCandidates = findPivotColumnCandidates(pivots, n);
            if (enteringVarCandidates.size() == 0)
            {
                numberStepsLastLP = step+1;
                return z;
            }
            sort(enteringVarCandidates.begin(), enteringVarCandidates.end(), [](pair<variable, unsigned> a, pair<variable, unsigned> b)
                 { return (a.first.value > b.first.value); });
            
            unsigned candidateIndex = 0;
            
            vector<double> d(m);
            unsigned enteringVarLabel;
            unsigned enteringColIndex;
            pair<variable, unsigned> pivotRow;
            unsigned leavingRowIndex;
            variable leavingVar;

            // 2. find the leaving variable
            for (candidateIndex = 0; candidateIndex < enteringVarCandidates.size(); candidateIndex++)
            {

                enteringVarLabel = enteringVarCandidates[candidateIndex].first.label;
                enteringColIndex = enteringVarCandidates[candidateIndex].second;

                pivotRow = findPivotRow(enteringVarLabel, pivots, n, d);
                leavingRowIndex = pivotRow.second;
                leavingVar = pivotRow.first;

                if (leavingRowIndex == -1)
                {
                    printf("\nThe given LP is unbounded.\n");
                    return numeric_limits<double>::infinity();
                }

                // check the diagonal element in the eta column
                // to see if the current choice of entering variable has to be rejected
                if (d[leavingRowIndex] > epsilon2)
                {
                    break;
                }      
            }
            if(candidateIndex == enteringVarCandidates.size()){
                //here we check if we saw all candidates for entering variables and none of them provide
                // a diagonal entry larger than epsilon2
                numberStepsLastLP = step+1;
                return z;
            }
            // 3. update the current basic feasible solution
            // set entering variable at t = min Ratio
            // other basic variables = x*_B - td
            // replace leaving column of B by the entering column , or add d as an eta matrix

            variable enteringVar = {enteringVarLabel, leavingVar.value};
            basic[leavingRowIndex] = enteringVar;

            for (size_t row = 0; row < m; ++row)
            {
                if (row != leavingRowIndex)
                {
                    basic[row].value = max(basic[row].value - d[row] * leavingVar.value, 0.0);
                }
            }

            // push a new eta matrix onto the vector
            eta pivot = {leavingRowIndex, d};
            pivots.push_back(pivot);

            nonBasic[enteringColIndex] = leavingVar.label;
            z = 0;
            for (size_t i = 0; i < basic.size(); i++)
            {
                z += basic[i].label < n ? basic[i].value : 0.0;
            }
        }
    }
    pair<variable, unsigned> findPivotRow(unsigned pivotColumn, vector<eta> pivots, int n, vector<double> &d)
    {
        // initialise d to be a the entering column

        for (unsigned j = 0; j < m; j++)
        {
            unsigned col = pivotColumn;
            if (col < n)
            {
                const double *__restrict cs = getCoefs(col);
                auto rs = getRules(col);
                auto nbElements = colOffsets[col + 1] - colOffsets[col];

                for (unsigned k = 0; k < nbElements; k++)
                    d[rs[k]] = cs[k];
            }
            else
            {
                d[col - n] = 1.0;
            }
        }

        // Go through eta matrices from pivot 1 to pivot k
        for (auto iter = pivots.cbegin(); iter != pivots.cend(); ++iter)
        {
            eta pivot = *iter;
            size_t rowToChange = pivot.col;
            double dOriginal = d[rowToChange];

            d[rowToChange] = dOriginal / pivot.values[rowToChange];

            for (size_t row = 0; row < d.size(); ++row)
            {
                if (row != rowToChange)
                {
                    d[row] = d[row] - pivot.values[row] * d[rowToChange];
                }
            }
        }
        // compute t for each b[i].value / d[i]
        // where d[i] > 0
        // choose the corresponding i for the smallest ratio
        // as the leaving variable
        double value;
        unsigned leavingVarLabel = -1;
        unsigned leavingRowIndex = -1;

        double ratio;
        double minRatio = numeric_limits<double>::max();

        for (size_t row = 0; row < d.size(); ++row)
        {
            if (d[row] > 0.0)
            {
                ratio = basic[row].value / d[row];
                if (ratio < minRatio)
                {
                    leavingVarLabel = basic[row].label;
                    leavingRowIndex = row;
                    minRatio = ratio;
                }
            }
        }
        variable leavingVar = {leavingVarLabel, minRatio};
        return make_pair(leavingVar, leavingRowIndex);
    }
    vector<pair<variable, unsigned>> findPivotColumnCandidates(vector<eta> pivots, int n)
    {
        // we solve the system yB = C_b using eta matrices
        // initialise y to C_b

        vector<double> y(m);

        for (size_t i = 0; i < m; i++)
        {
            y[i] = basic[i].label < n ? 1.0 : 0.0;
        }
        // solving y using a succession of BTRAN operations
        for (auto rIter = pivots.crbegin(); rIter != pivots.crend(); ++rIter)
        {
            eta pivot = *rIter;
            size_t colToChange = pivot.col;
            double yOriginal = y[colToChange];

            for (size_t row = 0; row < pivot.values.size(); ++row)
            {
                if (row != colToChange)
                {
                    yOriginal -= pivot.values[row] * y[row];
                }
            }

            double yNew = yOriginal / pivot.values[colToChange];
            y[colToChange] = yNew;
        }
        vector<pair<variable, unsigned>> rs;

        for (unsigned i = 0; i < n; i++)
        {
            unsigned col = nonBasic[i];
            double r;
            if (col < n)
            {
                const double *__restrict cs = getCoefs(col);
                auto rs = getRules(col);
                auto nbElements = colOffsets[col + 1] - colOffsets[col];

                r = 1.0;
                for (unsigned k = 0; k < nbElements; k++)
                    r -= y[rs[k]] * cs[k];
            }
            else
            {
                r = -y[col - n];
            }
            if (r > epsilon1)
            {
                variable v = {col, r};
                rs.push_back(make_pair(v, i));
            }
        }
        return rs;
    }
};
