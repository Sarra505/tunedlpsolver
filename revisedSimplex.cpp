// during a revised simplex iteration:

// 2. find the leaving variable
// calculate d with d=B^(-1)a with a being the entering column => solve the system Bd = a
// compute the largest ratio x*_B/d
// the index of that t is the leaving variable

// 3. update the current basic feasible solution
// set entering variable at t
// other basic variables = x*_B - td
// replace leaving column of B by the entering column

#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>
#include <span>
#include <string>
#include <vector>

class LPSolver
{
public:
   // a structure that represents a basic variable
   // containing its label (ranges from 0 to m + n - 1) and its current value
   struct variable
   {
      size_t label;
      float value;
   };
   // a structure that represents an eta matrix
   // containing the column that's different from an identity matrix (ranges from 1 to m - 1)
   // and the values in that column
   struct eta
   {
      size_t col;
      vector<double> values;
   };
   // The sparse coefficient matrix M
   std::vector<Entry> matrix;
   /// The number of rules
   unsigned m = 0;
   // indices of non basic variables
   std::vector<unsigned> nonBasic;
   // indices and current values of basic variables
   std::vector<variable> basic;

   float solveSimplex(unsigned n, unsigned stepLimit = ~0u);
   void solveEtaFTRAN(eta matrix, vector<float> *b);
   void solveEtaBTRAN(eta matrix, vector<float> *b);
   std::vector<Entry> getAllRow(unsigned int label);
   std::vector<Entry> getAllColumn(unsigned int label);
};
/// An entry in the sparse coefficient matrix
struct Entry
{
   /// The coefficient
   float coef;
   /// The corresponding rule/row
   unsigned rule;
   /// The corresponding variable/column
   unsigned variable;
};
struct eta
{
   size_t col;
   vector<double> values;
};
using namespace std;
float LPSolver::solveSimplex(unsigned n, unsigned stepLimit)
// Solve the given LP using simplex
{
   if (n == 0)
      return 0.0;
   if (m == 0)
      return numeric_limits<float>::infinity();
   // initialise c vector which is (1,1,...,1,0,...,0)
   float objFuncCoeff[n + m];
   size_t i;
   for (i = 0; i < n; i++)
   {
      objFuncCoeff[i] = 1.0;
   }
   for (i = n; i < m; i++)
   {
      objFuncCoeff[i] = 0.0f;
   }
   float objFuncCoeff[n + m];

   // Initialize the nonbasic variables' labels to be 0 through (n - 1)
   for (size_t i = 0; i < n; ++i)
   {
      nonBasic[i] = i;
   }
   // Initialize the basic variables to be n through m and their values to be the right-hand-side vector b which is 1
   for (size_t i = 0; i < n; ++i)
   {
      basic[i].label = i;
      basic[i].value = 1.0f;
   }

   // problem is initially feasible since b >= 0
   // An array of eta matrices representing previous pivots
   // the most recent pivot should be at the beginning
   vector<eta> pivots{};

   // initial value of objective function
   float z = 0;

   unsigned step = 0;
   bool foundSolution = false;
   for (; step < stepLimit; step++)
   {
      // 1. find the entering variable: solve the system yB = C_b
      // if none is found the problem is optimal

      // we solve the system yB = C_b using eta matrices
      // initialise y to C_b

      vector<float> y(m);

      for (size_t i = 0; i < m; i++)
      {
         y[i] = objFuncCoeff[basic[i].label];
      }
      // solving y using a succession of BTRAN operations
      for (auto rIter = pivots.crbegin(); rIter != pivots.crend(); ++rIter)
      {
         eta pivot = *rIter;
         solveEtaBTRAN(pivot, &y);
      }
      // now calculate C_n -yA_n and then choose A_j column with
      // first positive coefficient --> pricing

      float r;
      unsigned int nonBasicColumnIndex;
      unsigned int enteringIndex = -1;
      vector<Entry> nonBasicColumn;
      for (size_t i = 0; i < n; i++)
      {
         nonBasicColumnIndex = nonBasic[i];
         nonBasicColumn = getAllColumn(nonBasicColumnIndex);
         r = objFuncCoeff[nonBasicColumnIndex];
         for (size_t i = 0; i < nonBasicColumn.size(); i++)
         {
            r -= nonBasicColumn[i].coef * y[nonBasicColumn[i].rule];
         }
         if (r > 0)
         {
            enteringIndex = nonBasicColumnIndex;
            break;
         }
      }
      // if no entering index has been selected, no cadidates for entering variables, the optimal is reached
      if (enteringIndex == -1)
      {
         printf("\nNo entering var. Optimal value of %5.3f has been reached.\n", z);
         return 0;
      }

      // we solve the system Bd = a using eta matrices
      // initialise d to be a the entering column

      vector<float> d(m);

      for (size_t i = 0; i < nonBasicColumn.size(); i++)
      {
         d[nonBasicColumn[i].rule] = nonBasicColumn[i].coef;
         // the rest of indices should be filled with 0.0 !!!
      }
      // solving d using a succession of FTRAN operations
      for (auto rIter = pivots.crbegin(); rIter != pivots.crend(); ++rIter)
      {
         eta pivot = *rIter;
         solveEtaFTRAN(pivot, &d);
      }
      unsigned int leavingLabel = -1;
      unsigned int leavingRow = -1;
      float smallest_t = -1;
      // compute t for each b[i].value / d[i]
      // where d[i] > 0
      // choose the corresponding i for the smallest ratio
      // as the leaving variable

      // initialize smallest_t to be the first ratio where
      // the coefficient of the entering variable in that row is negative
      for (size_t row = 0; row < d.size(); ++row)
      {
         if (d[row] > 0.0)
         {
            leavingLabel = basic[row].label;
            leavingRow = row;
            smallest_t = basic[row].value / d[row];
         }
      }

      // if no ratio is computed, then the LP is unbounded
      if (leavingLabel == -1)
      {
         printf("\nThe given LP is unbounded. The family of solutions is:\n");
         return 0;
      }
   }

   if (!foundSolution)
      return -1.0f;

   float result = 0.0f;

   return result;
}

// function to solve Ex = b with E being an eta matrix
void LPSolver::solveEtaFTRAN(eta matrix, vector<float> *b)
{
   size_t pivot = matrix.col;
   (*b)[pivot] = (*b)[pivot] / matrix.values[pivot];

   for (size_t i = 0; i < m; i++)
   {
      if (i != pivot)
      {
         (*b)[i] -= (*b)[pivot] * matrix.values[i];
      }
   }
}
// function to solve xE = b with E being an eta matrix
void LPSolver::solveEtaBTRAN(eta pivotMatrix, vector<float> *b)
{
   size_t colToChange = pivotMatrix.col;
   float bp = (*b)[colToChange];

   for (size_t i = 0; i < m; ++i)
   {
      if (i != colToChange)
      {
         bp -= pivotMatrix.values[i] * (*b)[i];
      }
   }

   float bNew = bp / pivotMatrix.values[colToChange];
   (*b)[colToChange] = bNew;
}
// function that retrieves all entries within a certain row/column
std::vector<Entry> LPSolver::getAllRow(unsigned int row)
{
   std::vector<Entry> result;
   for (size_t i = 0; i < matrix.size(); i++)
   {
      if (matrix[i].rule == row)
      {
         result.push_back(matrix[i]);
      }
   }
   return result;
}
std::vector<Entry> LPSolver::getAllColumn(unsigned int column)
{
   std::vector<Entry> result;
   for (size_t i = 0; i < matrix.size(); i++)
   {
      if (matrix[i].variable == column)
      {
         result.push_back(matrix[i]);
      }
   }
   return result;
}