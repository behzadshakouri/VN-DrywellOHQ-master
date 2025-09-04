#ifndef FIELDGENERATOR_H
#define FIELDGENERATOR_H

#include <map>
#include <string>
#include <vector>
#include <utility>  // for std::pair
#include <gsl/gsl_rng.h>
#include <Matrix_arma.h>
#include <Vector_arma.h>
#include "TimeSeries.h"

// Correlation matrix structure for kriging
    struct correl_mat_vec
{
    CMatrix_arma M_22;  // Correlation matrix between determined points
    CVector_arma V_21;               // Correlation vector between target and determined points
    CVector_arma V_RHS;              // Values at determined points
};

struct PropertyPoint {
    // Track which parameters have been determined (for conditional generation)
    std::map<std::string, bool> determined;

    // Default constructor
    PropertyPoint() = default;

    // Check if a parameter has been determined
    bool isDetermined(const std::string& param) const;

    // Set determination status for a parameter
    void setDetermined(const std::string& param, bool isDet = true);

    // Reset all determination flags
    void reset();
};

class FieldGenerator
{
public:
    // Constructors
    FieldGenerator();
    FieldGenerator(int seed);
    FieldGenerator(unsigned int numPoints, int seed = 0);

    // Destructor
    ~FieldGenerator();

    // Copy constructor and assignment operator (needed for GSL cleanup)
    FieldGenerator(const FieldGenerator& other);
    FieldGenerator& operator=(const FieldGenerator& other);

    // Normal score field generation with conditional simulation
    void generateNormalScoreField(const std::string& fieldName,
                                  double correlationLength,
                                  const std::vector<std::pair<int, double>>& knownPoints = {});

    // Save field to CSV file
    bool saveFieldToCSV(const std::string& fieldName, const std::string& filename) const;

    // Getters and setters
    double getDx() const;
    void setDx(double gridSpacing);

    // Statistical functions
    double mean(const std::string& fieldName) const;
    double standardDeviation(const std::string& fieldName) const;
    double meanLog(const std::string& fieldName) const;
    double standardDeviationLog(const std::string& fieldName) const;
    double autoCorrelation(const std::string& fieldName, double distance) const;
    double minimum(const std::string& fieldName) const;
    double maximum(const std::string& fieldName) const;

    // Interpolation
    double interpolate(const std::string& fieldName, double x) const;


    // Distribution transformations
    void normalToUniform(const std::string& sourceField, const std::string& targetField);
    void normalToCDF(const std::string& sourceField, const std::string& targetField, const TimeSeries<double> &CDF);
    void exponentialTransform(const std::string& sourceField, const std::string& targetField, double a, double b);

private:
    // Field properties
    std::vector<PropertyPoint> points;
    std::map<std::string, std::vector<double>> parameters;

    // Spatial properties
    double dx = 1.0;  // Grid spacing
    double correlation_length_scale = 10.0;  // Correlation length

    // GSL random number generator
    const gsl_rng_type* rng_type;
    gsl_rng* rng;

    // Helper methods for RNG
    void initializeRNG(int seed);
    unsigned long generateSeed(int seed) const;
    void cleanup();

    // Helper methods for conditional generation
    std::vector<int> getDeterminedIndices(const std::string& fieldName) const;
    correl_mat_vec buildCorrelationMatrix(const std::string& fieldName, int targetIndex,
                                          double correlationLength) const;
    void buildCorrelationVectors(correl_mat_vec& correlationData,
                                 const std::vector<int>& determinedIndices,
                                 const std::string& fieldName, int targetIndex,
                                 double correlationLength) const;
    void buildCorrelationMatrix(correl_mat_vec& correlationData,
                                const std::vector<int>& determinedIndices,
                                double correlationLength) const;
    CVector_arma solveKrigingSystem(const correl_mat_vec& correlationData) const;
    void initializeKnownPoints(const std::string& fieldName,
                               const std::vector<std::pair<int, double>>& knownPoints);
    void initializeFirstPoint(const std::string& fieldName);
    unsigned int selectRandomUndeterminedPoint(const std::string& fieldName) const;
    bool allPointsDetermined(const std::string& fieldName) const;
    void generateConditionalValue(const std::string& fieldName, unsigned int targetIndex,
                                  double correlationLength);
    double computeCorrelation(int i, int j, double correlationLength) const;

    // Validation methods
    void validateFieldGeneration() const;
    void validateTargetIndex(unsigned int targetIndex, const std::string& fieldName) const;
    bool hasAnyDeterminedPoints(const std::string& fieldName) const;
};

#endif // FIELDGENERATOR_H
