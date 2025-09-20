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
#include "TimeSeriesSet.h"


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
    FieldGenerator(const FieldGenerator& other) = delete;
    FieldGenerator& operator=(const FieldGenerator& other) = delete;

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

    // Interpolate field value at any x position using linear interpolation
    double interpolateAt(double x, const std::string& quantity = "K_sat_normal_score") const;

    // Interpolate multiple points at once
    std::vector<double> interpolateAt(const std::vector<double>& xPositions,
                                      const std::string& quantity = "K_sat_normal_score") const;

    // Get the x-coordinate for a given grid index
    double indexToPosition(unsigned int index) const;

    // Get the grid index closest to a given x position
    unsigned int positionToIndex(double x) const;


    // Find the closest grid node index to a given x position
    unsigned int findClosestNode(double x) const;

    // Set known values from a TimeSeries where t=x position and value=field value
    void setKnownPointsFromTimeSeries(const std::string& fieldName,
                                      const TimeSeries<double>& timeSeries,
                                      double correlationLength);

    // Write field to ASCII format file
    bool writeFieldToCSV(const std::string& fieldName, const std::string& filename,
                         int precision = 6, bool includeHeader = true) const;

    // Generate field from measurement data file with full workflow
    bool generateFieldFromMeasurementData(const std::string& dataFilePath,
                                          const std::string& parameterName,
                                          const std::string& fieldName,
                                          double correlationLength,
                                          bool verbose = true);

    // Getter for measured CDFs
    const TimeSeriesSet<double>& getMeasuredCDFs() const;

    // Check if CDFs are available
    bool hasMeasuredCDFs() const;

    // Clear stored CDFs
    void clearMeasuredCDFs();

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
    void validateFieldForInterpolation(const std::string& quantity) const;


    // Member variable to store measured CDFs
    TimeSeriesSet<double> measured_CDFs_;
    bool has_measured_CDFs_ = false;
};

#endif // FIELDGENERATOR_H
