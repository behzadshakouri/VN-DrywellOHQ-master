#pragma once
#include "fieldgenerator.h"
#include <string>
#include <vector>

// Helper to generate, analyze, and save a stochastic field
void generateAndAnalyzeField(FieldGenerator &gen,
                             const std::string &csvPath,
                             const std::string &paramName,
                             double correlationLength,
                             const std::vector<double> &testDistances,
                             const std::string &outputPrefix);
