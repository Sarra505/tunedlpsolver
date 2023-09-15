#pragma once

struct variable
{
   unsigned label;
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