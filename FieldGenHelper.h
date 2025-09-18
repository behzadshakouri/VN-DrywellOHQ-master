#pragma once
#include "fieldgenerator.h"
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>

// Helper to generate, analyze, and save a stochastic field
inline void generateAndAnalyzeField(FieldGenerator &gen,
                                    const std::string &csvPath,
                                    const std::string &paramName,
                                    double correlationLength,
                                    const std::vector<double> &testDistances,
                                    const std::string &outputPrefix)
{
    std::string fieldName = paramName + "_normal";

    std::cout << "Reading data from '" << csvPath << "'" << std::endl;

    // Generate field
    gen.generateFieldFromMeasurementData(csvPath, paramName, fieldName, correlationLength);

    // Save CDFs
    gen.getMeasuredCDFs().write(outputPrefix + "Params_CDFs.txt");

    // Print stats
    std::cout << "Mean: " << gen.mean(fieldName) << std::endl;
    std::cout << "Standard Deviation: " << gen.standardDeviation(fieldName) << std::endl;
    std::cout << "Field size: " << 1000 << " points" << std::endl;
    std::cout << "Grid spacing: " << gen.getDx() << " m" << std::endl;
    std::cout << "Total length: " << 1000 * gen.getDx() << " m" << std::endl;

    // Auto-correlation analysis
    std::cout << "\n=== AUTO-CORRELATION ANALYSIS for " << paramName << " ===" << std::endl;
    std::cout << "Distance (m)  |  Correlation  |  Theoretical*" << std::endl;
    std::cout << "------------- | ------------- | -------------" << std::endl;

    for (double dist : testDistances) {
        double empirical = gen.autoCorrelation(fieldName, dist);
        double theoretical = std::exp(-dist / correlationLength);
        std::cout << std::setw(10) << dist << "   |  "
                  << std::setw(10) << empirical << "   |  "
                  << std::setw(10) << theoretical << std::endl;
    }

    // Save fields
    std::string normalFile = outputPrefix + fieldName + "_field.csv";
    std::string realFile   = outputPrefix + paramName + "_field.csv";

    if (gen.saveFieldToCSV(fieldName, normalFile))
        std::cout << "Saved " << normalFile << " successfully!\n";
    else
        std::cout << "Failed to save " << normalFile << "\n";

    if (gen.saveFieldToCSV(paramName, realFile))
        std::cout << "Saved " << realFile << " successfully!\n";
    else
        std::cout << "Failed to save " << realFile << "\n";
}
