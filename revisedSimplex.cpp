// during a revised simplex iteration:

// 1. find the entering variable: solve the system yB = C_b and then calculate C_n -yA_n and then choose A_j column with
// first positive coefficient
// if none is found the problem is optimal

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
   float solveSimplex(unsigned n, unsigned stepLimit = ~0u);
   vector<float> solveEtaFTRAN(eta matrix, vector<float> b);
   vector<float> solveEtaBTRAN(eta matrix, vector<float> b);
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

   unsigned step = 0;
   bool foundSolution = false;
   for (; step < stepLimit; step++)
   {
   }

   if (!foundSolution)
      return -1.0f;

   float result = 0.0f;

   return result;
}
// function to solve yB = c_b for y

// LU factroization of a sparsely-represented matrix B which is mxm

void factorizeLU()
{
}

// function to solve Ex = b with E being an eta matrix
vector<float> LPSolver::solveEtaFTRAN(eta matrix, vector<float> b)
{
   size_t pivot = matrix.col;
   b[pivot] = b[pivot] / matrix.values[pivot];

   for (size_t i = 0; i < m; i++)
   {
      if (i != pivot)
      {
         b[i] -= b[pivot] * matrix.values[i];
      }
   }
}
// function to solve xE = b with E being an eta matrix
vector<float> LPSolver::solveEtaBTRAN(eta pivotMatrix, vector<float> b)
{
   size_t colToChange = pivotMatrix.col;
   double bp = b[colToChange];

   for (size_t i = 0; i < m; ++i)
   {
      if (i != colToChange)
      {
         bp -= pivotMatrix.values[i] * b[i];
      }
   }

   double bNew = bp / pivotMatrix.values[colToChange];
   b[colToChange] = bNew;
}
