// during a revised simplex iteration:

// 1. find the entering variable: solve the system yB = C_b and then calculate C_n -yA_n and then choose A_j column with 
//first positive coefficient
//if none is found the problem is optimal




// 2. find the leaving variable
// calculate d with d=B^(-1)a with a being the entering column => solve the system Bd = a
// compute the largest ratio x*_B/d 
//the index of that t is the leaving variable


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

class LPSolver {
   public:
    /// The number of rules
   unsigned m = 0;
   float solveSimplex(unsigned n, unsigned stepLimit = ~0u);

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
   for (; step < stepLimit; step++) {
      
   }

   if (!foundSolution)
      return -1.0f;

   float result = 0.0f;

   return result;
}