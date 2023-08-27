#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>
#include <span>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

struct variable {
    size_t label;
    double value;
};

struct Entry
{
   /// The coefficient
   double coef;
   /// The corresponding rule/row
   unsigned rule;
   /// The corresponding variable/column
   unsigned variable;
};
class LPSolver
{
public:
   // a structure that represents a basic variable
   // containing its label (ranges from 0 to m + n - 1) and its current value
   struct variable
   {
      size_t label;
      double value;
   };

   // a structure that represents an eta matrix
   // containing the column that's different from an identity matrix (ranges from 1 to m - 1)
   // and the values in that column
   struct eta
   {
      size_t col;
      std::vector<double> values;
   };
   // The sparse coefficient matrix M
   std::vector<Entry> matrix;

   // the matrix is also stored in CCS format to facilitate column-based operations
   /// coefficients of non-zero elements in M
   std::vector<double> coefs;
   /// row indices of non-zero elements in M
   std::vector<unsigned> rules;
   /// colOffsets[i] tells you where the i-th column starts in the rules and coefs arrays
   std::vector<unsigned> colOffsets;

   /// The number of rules
   unsigned m = 0;
   // indices of non basic variables
   std::vector<unsigned> nonBasic;
   // indices and current values of basic variables
   std::vector<variable> basic;

   LPSolver(std::vector<Entry> coefMatrix, unsigned numRules)
   {
      matrix = coefMatrix;
      m = numRules;
   }
   /// Prepare for computation
   void prepareCCSFormat(unsigned n);
   double solveSimplex(unsigned n, unsigned stepLimit = ~0u);
   auto getCoefs(unsigned i);
   auto getRules(unsigned i);
};

using namespace std;
// zero tolerances
static const double epsilon1 = 0.00001;
static const double epsilon2 = 0.00000001;

auto LPSolver::getCoefs(unsigned i) { return coefs.data() + colOffsets[i]; };
auto LPSolver::getRules(unsigned i) { return rules.data() + colOffsets[i]; };
void LPSolver::prepareCCSFormat(unsigned n)
// Prepare for computation
{
   // Build optimized column representation for fast operations
   colOffsets.clear();
   colOffsets.resize(n + 1);
   rules.resize(matrix.size());
   coefs.resize(matrix.size());

   sort(matrix.begin(), matrix.end(), [](auto &a, auto &b)
        { return tie(a.variable, a.rule) < tie(b.variable, b.rule); });

   for (size_t i = 0; i < matrix.size(); i++)
   {
      assert(matrix[i].rule < m);
      assert(matrix[i].variable < n);
      rules[i] = matrix[i].rule;
      coefs[i] = matrix[i].coef;
      colOffsets[matrix[i].variable + 1]++;
   }

   for (unsigned i = 1; i <= n; i++)
      colOffsets[i] += colOffsets[i - 1];
}
double LPSolver::solveSimplex(unsigned n, unsigned stepLimit)
// Solve the given LP using simplex
{
   if (n == 0)
      return 0.0;
   if (m == 0)
      return numeric_limits<double>::infinity();

   prepareCCSFormat(n);
   // initialise c vector which is (1,1,...,1,0,...,0)
   double objFuncCoef[n + m];
   size_t i;
   for (i = 0; i < n; i++)
   {
      objFuncCoef[i] = 1.0;
   }
   for (i = n; i < m + n; i++)
   {
      objFuncCoef[i] = 0.0;
   }
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
   // initial value of objective function
   double z = 0.0;
   /********************************************************************************/

   for (; step < stepLimit; step++)
   {
      // 1. find the entering variable: solve the system yB = C_b
      // if none is found the problem is optimal

      // we solve the system yB = C_b using eta matrices
      // initialise y to C_b

      vector<double> y(m);

      for (size_t i = 0; i < m; i++)
      {
         y[i] = objFuncCoef[basic[i].label];
      }
      // solving y using a succession of BTRAN operations (in-place)
      for (auto rIter = pivots.crbegin(); rIter != pivots.crend(); ++rIter)
      {
         eta pivot = *rIter;
         size_t colToChange = pivot.col;
            double yOriginal = y[colToChange];
            
            for (size_t row = 0; row < pivot.values.size(); ++row) {
                if (row != colToChange) {
                    yOriginal -= pivot.values[row] * y[row];
                }
            }
            
            double yNew = yOriginal / pivot.values[colToChange];
            y[colToChange] = yNew;
      }
      // now calculate C_n -yA_n and then choose A_j column with
      // first positive coefficient --> pricing

      vector<variable> rs;
      unsigned int enteringLabel = -1;
      unsigned int enteringIndex = -1;
      double bestR = 1e-5;

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
            rs.push_back(v);
                
            if (r > bestR) {
               bestR = r;
               enteringLabel = col;
               enteringIndex = i;
            }
            // break;
         }
      }
      // if no entering index has been selected, no cadidates for entering variables, the optimal is reached
      if (enteringLabel == -1)
      {
         printf("\nNo entering var. Optimal value of %5.3f has been reached.\n", z);
         foundSolution = true;
         break;
      }
      sort(rs.begin(), rs.end(), [](variable a, variable b) {
            return (a.value > b.value);
        }
      );

      // 2. find the leaving variable
      // calculate d with d=B^(-1)a with a being the entering column => solve the system Bd = a
      // compute the largest ratio x*_B/d
      // the index of that t is the leaving variable
      // we solve the system Bd = a using eta matrices

      size_t enteringVariable_index = 0;
        
        // compute the column d in Anbar
        // for the entering variable
        // using eta matrices (Bd = a)
        
        size_t leavingLabel;
        size_t leavingRow;
        double smallest_t;
        vector<double> d(m);
        
        while (true) {
            
            leavingLabel = -1;
            leavingRow = -1;
            smallest_t = numeric_limits<float>::infinity();
            
            if (enteringVariable_index > 0) {
                printf("\n\nRechoosing entering variable since the diagonal element in the eta column is close to zero.\n");
            }
            
            if (enteringVariable_index < rs.size()) {
                enteringLabel = rs[enteringVariable_index].label;
                
                if (enteringVariable_index > 0) {
                    printf("Entering variable is x%lu \n", enteringLabel + 1);
                }
            } else {
                printf("\nNo entering var. Optimal value of %5.3f has been reached.\n", z);
                foundSolution = true;
                return 0;
            }
            
            // initialise d to be a the entering column
            

            for (unsigned j = 0; j < m; j++)
            {
               unsigned col = enteringLabel;
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
            for (auto iter = pivots.cbegin(); iter != pivots.cend(); ++iter) {
                eta pivot = *iter;
                size_t rowToChange = pivot.col;
                double dOriginal = d[rowToChange];
                
                d[rowToChange] = dOriginal / pivot.values[rowToChange];
                
                for (size_t row = 0; row < d.size(); ++row) {
                    if (row != rowToChange) {
                        d[row] = d[row] - pivot.values[row] * d[rowToChange];
                    }
                }
            }
            
            // print out d (abarj)
            printf("d = ");
            
            for (auto iter = d.cbegin(); iter != d.cend(); ++iter) {
                printf("%5.3f ", *iter);
            }
            
            printf("\n");
            
            // compute t for each b[i].value / d[i]
            // where d[i] > 0
            // choose the corresponding i for the smallest ratio
            // as the leaving variable
            for (size_t row = 0; row < d.size(); ++row) {
                if (d[row] > 0.0) {
                    leavingLabel = basic[row].label;
                    leavingRow = row;
                    smallest_t = basic[row].value / d[row];
                }
            }
            
            // if no ratio is computed, then the LP is unbounded
            if (leavingLabel == -1) {
                printf("\nThe given LP is unbounded. The family of solutions is:\n");
                return 0;
            }
            
            
            for (size_t row = 0; row < d.size(); ++row) {
                if (d[row] < 0.0) {
                    continue;
                }
                
                double t_row = basic[row].value / d[row];
                
                if (t_row >= 0.0) {
                    printf("x%lu %5.3f ", basic[row].label + 1, t_row);
                }
                
                if (t_row < smallest_t) {
                    leavingLabel = basic[row].label;
                    leavingRow = row;
                    smallest_t = t_row;
                }
            }
            
            // check the diagonal element in the eta column
            // to see if the current choice of entering variable has to be rejected
            if (d[leavingRow] > epsilon2) {
                printf("\nLeaving variable is x%lu\n", leavingLabel + 1);
                break;
            } else {
                enteringVariable_index++;
                continue;
            }
        }

      // 3. update the current basic feasible solution
      // set entering variable at t
      // other basic variables = x*_B - td
      // replace leaving column of B by the entering column , or add d as an eta matrix
      variable enteringVar = {enteringLabel, smallest_t};
      basic[leavingRow] = enteringVar;

      for (size_t row = 0; row < m; ++row)
      {
         if (row != leavingRow)
         {
            basic[row].value = max(basic[row].value - d[row] * smallest_t, 0.0);
         }
      }

      // push a new eta matrix onto the vector
      eta pivot = {leavingRow, d};
      pivots.push_back(pivot);

      nonBasic[enteringIndex] = leavingLabel;

      // recalculate z TODO: change thisss!!!! to  a more efficient implem
      z = 0;
      for (size_t i = 0; i < basic.size(); i++)
      {
         z += basic[i].label < n ? basic[i].value : 0.0;
      }
   }

   if (!foundSolution)
      return -1.0;

   return z;
}
bool CompareMaxColumn(const Entry &_a, const Entry &_b)
{
   return _a.variable < _b.variable;
}
// a function for sorting variables in a vector
// into descending order
bool mComparator(variable v1, variable v2) {
    return v1.value > v2.value;
}

std::tuple<std::vector<Entry>, unsigned, unsigned> parseLPProblem(const std::string &lpString)
{
   std::vector<Entry> coefficientMatrix;

   // columns = number of structural variables + number of slack variables
   // number of structural variables = maximu column index + 1
   // number of slack variables = number of rules
   // rows = number of rules + 1
   // number of structural variables = max column index + 1
   std::istringstream iss(lpString);

   int numRules;
   iss >> numRules;

   Entry entry;

   // initialize A
   for (int i = 0; i < numRules; i++)
   {
      int numEntries;
      iss >> numEntries;

      for (int j = 0; j < numEntries; j++)
      {
         int columnNumber;
         float coefficient;      
         iss >> columnNumber >> coefficient;
         entry.coef = coefficient;
         entry.rule = i;
         entry.variable = columnNumber;

         // Add the entry to the matrix
         coefficientMatrix.push_back(entry);
      }
   }
   auto numVariables = *max_element(coefficientMatrix.begin(), coefficientMatrix.end(), CompareMaxColumn);

   return std::make_tuple(coefficientMatrix, numRules, numVariables.variable + 1);
}

int main()
{

   std::ifstream input("lp.txt");

   /*for (std::string line; getline(input, line);){
       myData = parseLPProblem(line);
   }*/

   std::string line;
   getline(input, line);
   std::tuple<std::vector<Entry>, unsigned, unsigned> myData = parseLPProblem(line);

   LPSolver simplex(get<0>(myData), get<1>(myData));

   double z = simplex.solveSimplex(get<2>(myData), 1000);
   printf("\nOptimal value of %5.3f has been reached.\n", z);

   return 0;
}
