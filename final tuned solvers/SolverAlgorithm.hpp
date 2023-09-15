#pragma once
#include "common.hpp"
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
#include <chrono>
#include <vector>
using namespace std;
// zero tolerances
static const double epsilon1 = 0.00001;
static const double epsilon2 = 0.000000001;
static const double epsilon3 = 1e-50;


class SolverAlgorithm
{
public:
    // The sparse coefficient matrix M
    std::vector<Entry> matrix;

    std::vector<double> B;
    std::vector<double> C;
    // the matrix is also stored in CCS format to facilitate column-based operations
    /// coefficients of non-zero elements in M
    std::vector<double> values;
    /// row indices of non-zero elements in M
    std::vector<unsigned> rows;
    /// colOffsets[i] tells you where the i-th column starts in the rules and coefs arrays
    std::vector<unsigned> colOffsets;
    /// The number of rules
    unsigned m = 0;
    // indices of non basic variables
    std::vector<unsigned> nonBasic;
    // indices and current values of basic variables
    std::vector<variable> basic;
    virtual double solve(unsigned n, unsigned stepLimit) = 0;
    void prepareCCSFormat(unsigned n);
    auto getCoefs(unsigned i);
    auto getRules(unsigned i);

};
auto SolverAlgorithm::getCoefs(unsigned i) { return values.data() + colOffsets[i]; };

auto SolverAlgorithm::getRules(unsigned i) { return rows.data() + colOffsets[i]; };

bool CompareMaxColumn(const Entry &_a, const Entry &_b)
{
   return _a.variable < _b.variable;
}

tuple<vector<Entry>, unsigned, unsigned> parseCoefficientMatrix(const string &lpString)
{
   vector<Entry> coefficientMatrix;

   // columns = number of decision variables + number of slack variables
   // number of decision variables = maximu column index + 1
   // number of slack variables = number of rules
   // rows = number of rules + 1
   // number of decision variables = max column index + 1
   
   istringstream iss(lpString);

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

   return make_tuple(coefficientMatrix, numRules, numVariables.variable + 1);
}

void SolverAlgorithm::prepareCCSFormat(unsigned n)
{
   // Build compressed column representation to faciliate fast column operations
   colOffsets.clear();
   colOffsets.resize(n + 1);
   rows.resize(matrix.size());
   values.resize(matrix.size());

   sort(matrix.begin(), matrix.end(), [](auto &a, auto &b)
        { return tie(a.variable, a.rule) < tie(b.variable, b.rule); });

   for (size_t i = 0; i < matrix.size(); i++)
   {
      assert(matrix[i].rule < m);
      assert(matrix[i].variable < n);
      rows[i] = matrix[i].rule;
      values[i] = matrix[i].coef;
      colOffsets[matrix[i].variable + 1]++;
   }

   for (unsigned i = 1; i <= n; i++)
      colOffsets[i] += colOffsets[i - 1];
}