#include <iostream>
#include <cmath>
#include <vector>
#include <sstream>
#include <iomanip>
#include <string>
#include <fstream>
#include <algorithm>
#include <cassert>
// zero tolerances
static const double epsilon1 = 0.00001;
static const double epsilon2 = 0.00000001; // TODO: find a way to filter out nan and -inf
// In IEEE 754, the most common binary representation of floating-point numbers, 
//the positive infinity is the value with all bits of the exponent set and all bits of the fraction cleared. 
//he NaN in C++ can occur in the given three conditions which are when the log of a negative number or 
//zero is taken when we divide a number by zero, and the square root of the negative number.
using namespace std;

struct lpProblemData
{
    std::vector<std::vector<double>> A;
    std::vector<double> B;
    std::vector<double> C;
};
class Simplex
{

private:
    int rows, cols;
    // stores coefficients of all the variables
    std::vector<std::vector<double>> A;  //TODO: change A to compact sparse representation
    // stores constants of constraints
    std::vector<double> B;
    // stores the coefficients of the objective function
    std::vector<double> C;

    double maximum;

    bool isUnbounded;
    bool degeneracy;

public:
    Simplex(std::vector<std::vector<double>> matrix, std::vector<double> b, std::vector<double> c)
    {
        maximum = 0;
        isUnbounded = false;
        rows = matrix.size();
        cols = matrix[0].size();
        A.resize(rows, vector<double>(cols, 0));
        B.resize(b.size());
        C.resize(c.size());

        for (int i = 0; i < rows; i++)
        { // pass A[][] values to the matrix
            for (int j = 0; j < cols; j++)
            {
                A[i][j] = matrix[i][j];
            }
        }

        for (int i = 0; i < c.size(); i++)
        { // pass c[] values to the B vector
            C[i] = c[i];
        }
        for (int i = 0; i < b.size(); i++)
        { // pass b[] values to the B vector
            B[i] = b[i];
        }
    }

    bool simplexAlgorithmCalculation()
    {
        // check whether the table is optimal,if optimal no need to process further
        if (checkOptimality() == true)
        {
            return true;
        }

        // find the column which has the pivot.The least coefficient of the objective function(C array).
        int pivotColumn = findPivotColumn();

        if (isUnbounded == true)
        {
            cout << "Error unbounded" << endl;
            return true;
        }

        // find the row with the pivot value.The least value item's row in the B array
        int pivotRow = findPivotRow(pivotColumn);

        // form the next table according to the pivot value
        doPivotting(pivotRow, pivotColumn);

        if(degeneracy == true){
            cout << "Degeneracy occured" << endl;
            return true;
        }

        return false;
    }

    bool checkOptimality()
    {
        // if the table has further negative constraints,then it is not optimal
        bool isOptimal = true;
        // as long as we still have negative coefficients, the optimum is not found
        for (int i = 0; i < C.size(); i++)
        {
            double value = C[i];
            if (value < 0)
            {
                isOptimal = false;
                break;
            }
        }
        // if all the constraints are positive now,the table is optimal
        return isOptimal;
    }
    //TODO: refactor after changing matrix A representation
    void doPivotting(int pivotRow, int pivotColumn)
    {

        double pivotValue = A[pivotRow][pivotColumn]; // gets the pivot value
        if(pivotValue == 0) {
            degeneracy = true;
            return;
            
        }
        /*Yes, a row with ratio 0 would be the pivot row. A zero ratio signals a condition called degeneracy. G
        In geometric terms, there are more constraints binding at the current vertex of the feasible region than the number needed to define a point. 
        (You need one binding constraint for each decision variable, excluding slack and surplus variables.)
        When the simplex method arrives at a degenerate corner, cycling can occur. There are a variety of adjustments you 
        can make to the basic simplex algorithm to avoid being stuck in a cycle. See the relevant Wikipedia page for a brief discussion.*/




        double pivotRowVals[cols]; // the column with the pivot

        double pivotColVals[rows]; // the row with the pivot

        double rowNew[cols]; // the row after processing the pivot value

        maximum = maximum - (C[pivotColumn] * (B[pivotRow] / pivotValue)); // set the maximum step by step
        // get the row that has the pivot value
        for (int i = 0; i < cols; i++)
        {
            pivotRowVals[i] = A[pivotRow][i];
        }
        // get the column that has the pivot value
        for (int j = 0; j < rows; j++)
        {
            pivotColVals[j] = A[j][pivotColumn];
        }

        // set the row values that has the pivot value divided by the pivot value and put into new row
        for (int k = 0; k < cols; k++)
        {
            rowNew[k] = pivotRowVals[k] / pivotValue;
        }

        B[pivotRow] = B[pivotRow] / pivotValue;

        // process the other coefficients in the A array by subtracting
        for (int m = 0; m < rows; m++)
        {
            // ignore the pivot row as we already calculated that
            if (m != pivotRow)
            {
                for (int p = 0; p < cols; p++)
                {
                    double multiplyValue = pivotColVals[m];
                    A[m][p] = A[m][p] - (multiplyValue * rowNew[p]);
                    // C[p] = C[p] - (multiplyValue*C[pivotRow]);
                    // B[i] = B[i] - (multiplyValue*B[pivotRow]);
                }
            }
        }

        // process the values of the B array
        for (int i = 0; i < B.size(); i++)
        {
            if (i != pivotRow)
            {

                double multiplyValue = pivotColVals[i];
                B[i] = B[i] - (multiplyValue * B[pivotRow]);
            }
        }
        // the least coefficient of the constraints of the objective function
        double multiplyValue = C[pivotColumn];
        // process the C array
        for (int i = 0; i < C.size(); i++)
        {
            C[i] = C[i] - (multiplyValue * rowNew[i]);
        }

        // replacing the pivot row in the new calculated A array
        for (int i = 0; i < cols; i++)
        {
            A[pivotRow][i] = rowNew[i];
        }
    }

    // find the least coefficients of constraints in the objective function's position
    int findPivotColumn()
    {

        int location = 0;
        double minm = C[0];

        for (int i = 1; i < C.size(); i++)
        {
            if (C[i] < minm)
            {
                minm = C[i];
                location = i;
            }
        }

        return location;
    }

    // find the row with the pivot value. conduct a ratio test only on strictly positive items
    int findPivotRow(int pivotColumn)
    {
        double min = 99999999;
        double value;
        int pivotRow = 0;
        int nonpositiveCount = 0;

        for (int i = 0; i < rows; i++)
        {
            value = A[i][pivotColumn];
            

            if (value > 0 && (value / B[i]) < min)
            {
                min = (value / B[i]);
                pivotRow = i;
            }
            else if (value <= 0)
            {
                nonpositiveCount++;
            }
        }
        if (nonpositiveCount == rows)
        {
            isUnbounded = true;
        } 
        return pivotRow;
    }

    void DisplaySimplex()
    {
        // Display the coefficient matrix
        for (const auto &ruleEntries : A)
        {
            for (const auto &coefficient : ruleEntries)
            {
                std::cout << coefficient << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
        for (const auto &element : B)
        {
            std::cout << element << " ";
        }
        std::cout << std::endl;
        std::cout << std::endl;
        for (const auto &element : C)
        {
            std::cout << element << " ";
        }
        std::cout << std::endl;
        std::cout << std::endl;   

    }
    void CalculateSimplex()
    {
        bool end = false;

        cout << "initial array(Not optimal)" << endl;
        DisplaySimplex();

        cout << " " << endl;
        cout << "final array(Optimal solution)" << endl;

        while (!end)
        {

            bool result = simplexAlgorithmCalculation();
            DisplaySimplex();
            cout << " " << endl;

            end = result;
        }
        cout << "Answers for the Constraints of variables" << endl;

        for (int i = 0; i < A.size(); i++)
        { // every basic column has the values, get it form B array
            int count0 = 0;
            int index = 0;
            for (int j = 0; j < rows; j++)
            {
                if (A[j][i] == 0.0)
                {
                    count0 += 1;
                }
                else if (A[j][i] == 1)
                {
                    index = j;
                }
            }

            if (count0 == rows - 1)
            {

                cout << "variable" << index + 1 << ": " << B[index] << endl; // every basic column has the values, get it form B array
            }
            else
            {
                cout << "variable" << index + 1 << ": " << 0 << endl;
            }
        }

        cout << "" << endl;
        cout << "maximum value: " << maximum << endl; // print the maximum values
    }
};
struct lpProblemData parseLPProblem(const std::string &lpString)
{
    std::vector<std::vector<double>> coefficientMatrix;
    std::vector<double> C;
    std::vector<double> B;

    // columns = number of structural variables + number of slack variables
    // number of structural variables = maximu column index + 1
    // number of slack variables = number of rules
    // rows = number of rules + 1
    // number of structural variables = max column index + 1
    std::istringstream iss(lpString);
    std::istringstream issCopy(lpString);

    int numRules;
    iss >> numRules;

    int maxCol = std::numeric_limits<int>::min();
    double value;
    int firstItem;
    issCopy >> firstItem;
    while (issCopy >> value)
    {
        if (value > maxCol)
        {
            maxCol = value;
        }
    }
    int cols = maxCol + 1 + numRules;
    int rows = numRules;
    // initialize B array
    B.resize(rows);
    for (int i = 0; i < B.size(); i++)
    {
        B[i] = 1;
    }
    // initialize C array
    C.resize(cols);
    for (int i = 0; i < maxCol + 1; i++)
    {
        C[i] = -1;
    }
    for (int i = maxCol + 1; i < C.size(); i++)
    {
        C[i] = 0;
    }

    // initialize A
    coefficientMatrix.resize(rows, vector<double>(cols, 0));
    for (int i = 0; i < numRules; i++)
    {
        int numEntries;
        iss >> numEntries;

        std::vector<double> ruleEntries;
        for (int j = 0; j < numEntries; j++)
        {
            int columnNumber;
            double coefficient;

            // Read the column number and coefficient as strings
            std::string columnNumberStr, coefficientStr;
            iss >> columnNumberStr >> coefficientStr;

            // Convert the strings to the respective types
            columnNumber = std::stoi(columnNumberStr);
            coefficient = std::stod(coefficientStr);

            // Add the coefficient to the rule entries
            coefficientMatrix[i][columnNumber] = coefficient;
        }
        coefficientMatrix[i][maxCol + 1 + i] = 1;
    }
    lpProblemData myData;
    myData.A = coefficientMatrix;
    myData.B = B;
    myData.C = C;

    return myData;
}

int main()
{

    std::ifstream input("lp.txt");
    lpProblemData myData;

    /*for (std::string line; getline(input, line);){
        myData = parseLPProblem(line);
    }*/

    std::string line;
    getline(input, line);
    myData = parseLPProblem(line);

    // hear the make the class parameters with A[m][n] vector b[] vector and c[] vector
    Simplex simplex(myData.A, myData.B, myData.C);

    simplex.CalculateSimplex();

    return 0;
}