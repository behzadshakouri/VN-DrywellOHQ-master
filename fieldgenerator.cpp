#include "fieldgenerator.h"
#ifndef _WINDOWS
#include <sys/time.h>
#endif
#include <ctime>
#include <stdexcept>
#include <gsl/gsl_randist.h>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <gsl/gsl_randist.h>
#include <cmath>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <gsl/gsl_cdf.h>
#include "TimeSeriesSet.h"


// Default constructor
FieldGenerator::FieldGenerator() : rng_type(gsl_rng_taus), rng(nullptr) {
    rng = gsl_rng_alloc(rng_type);
    initializeRNG(0);  // Use time-based seed
}

// Constructor with seed
FieldGenerator::FieldGenerator(int seed) : rng_type(gsl_rng_taus), rng(nullptr) {
    rng = gsl_rng_alloc(rng_type);
    initializeRNG(seed);
}

// Constructor with number of points and seed
FieldGenerator::FieldGenerator(unsigned int numPoints, int seed)
    : rng_type(gsl_rng_taus), rng(nullptr) {
    rng = gsl_rng_alloc(rng_type);
    points.resize(numPoints);
    initializeRNG(seed);
}

// Destructor
FieldGenerator::~FieldGenerator() {
    cleanup();
}

// Copy constructor
/*
FieldGenerator::FieldGenerator(const FieldGenerator& other)
    : points(other.points), parameters(other.parameters),
    dx(other.dx), correlation_length_scale(other.correlation_length_scale),
    rng_type(other.rng_type), rng(nullptr) {

    rng = gsl_rng_alloc(rng_type);
    gsl_rng_memcpy(rng, other.rng);  // Copy RNG state
}

// Assignment operator
FieldGenerator& FieldGenerator::operator=(const FieldGenerator& other) {
    if (this != &other) {
        // Clean up existing RNG
        cleanup();

        // Copy data
        points = other.points;
        parameters = other.parameters;
        dx = other.dx;
        correlation_length_scale = other.correlation_length_scale;
        rng_type = other.rng_type;

        // Initialize new RNG and copy state
        rng = gsl_rng_alloc(rng_type);
        gsl_rng_memcpy(rng, other.rng);
    }
    return *this;
}
*/

// Helper method to initialize RNG
void FieldGenerator::initializeRNG(int seed) {
    if (rng == nullptr) {
        throw std::runtime_error("RNG not allocated");
    }

    unsigned long actualSeed = generateSeed(seed);
    gsl_rng_set(rng, actualSeed);
}

// Generate seed based on time and input seed
unsigned long FieldGenerator::generateSeed(int seed) const {
#ifndef _WINDOWS
    struct timeval tv;
    if (gettimeofday(&tv, nullptr) != 0) {
        // Fallback if gettimeofday fails
        return static_cast<unsigned long>(time(nullptr)) + seed;
    }
    return static_cast<unsigned long>(tv.tv_sec) + tv.tv_usec + seed;
#else
    return static_cast<unsigned long>(time(nullptr)) + seed;
#endif
}

// Cleanup GSL resources
void FieldGenerator::cleanup() {
    if (rng != nullptr) {
        gsl_rng_free(rng);
        rng = nullptr;
    }
}

bool FieldGenerator::saveFieldToCSV(const std::string& fieldName, const std::string& filename) const {
    // Check if field exists
    auto it = parameters.find(fieldName);
    if (it == parameters.end()) {
        return false; // Field doesn't exist
    }

    const std::vector<double>& fieldData = it->second;

    // Open file for writing
    std::ofstream file(filename);
    if (!file.is_open()) {
        return false; // Failed to open file
    }

    // Set precision for floating point output
    file << std::fixed << std::setprecision(8);

    // Write header
    file << "z," << fieldName << "\n";

    // Write data
    for (size_t i = 0; i < fieldData.size(); ++i) {
        file << i*dx << "," << fieldData[i] << "\n";
    }

    file.close();
    return true;
}

// Getter for dx
double FieldGenerator::getDx() const {
    return dx;
}

// Setter for dx with validation
void FieldGenerator::setDx(double gridSpacing) {
    if (gridSpacing <= 0.0) {
        throw std::invalid_argument("Grid spacing (dx) must be positive, got: " +
                                    std::to_string(gridSpacing));
    }
    dx = gridSpacing;
}

bool PropertyPoint::isDetermined(const std::string& param) const {
    auto it = determined.find(param);
    return (it != determined.end()) ? it->second : false;
}

void PropertyPoint::setDetermined(const std::string& param, bool isDet) {
    determined[param] = isDet;
}

void PropertyPoint::reset() {
    determined.clear();
}

void FieldGenerator::generateNormalScoreField(const std::string& fieldName,
                                              double correlationLength,
                                              const std::vector<std::pair<int, double>>& knownPoints) {
    validateFieldGeneration();

    if (points.empty()) {
        throw std::invalid_argument("Field is empty - no points to generate");
    }

    // Initialize or resize parameter vector
    if (parameters.find(fieldName) == parameters.end()) {
        parameters[fieldName] = std::vector<double>(points.size(), 0.0);
    } else if (parameters[fieldName].size() != points.size()) {
        parameters[fieldName].resize(points.size(), 0.0);
    }

    // Reset determination flags for this field
    for (auto& point : points) {
        point.setDetermined(fieldName, false);
    }

    // Initialize known points if provided
    if (!knownPoints.empty()) {
        initializeKnownPoints(fieldName, knownPoints);
    }

    // If no known points, initialize first point randomly
    if (!hasAnyDeterminedPoints(fieldName)) {
        initializeFirstPoint(fieldName);
    }

    // Sequential simulation: generate remaining points conditionally
    while (!allPointsDetermined(fieldName)) {
        const unsigned int targetIndex = selectRandomUndeterminedPoint(fieldName);
        generateConditionalValue(fieldName, targetIndex, correlationLength);
    }
}

void FieldGenerator::initializeKnownPoints(const std::string& fieldName,
                                           const std::vector<std::pair<int, double>>& knownPoints) {
    for (const auto& knownPoint : knownPoints) {
        const int index = knownPoint.first;
        const double value = knownPoint.second;

        // Validate index
        if (index < 0 || index >= static_cast<int>(points.size())) {
            throw std::out_of_range("Known point index " + std::to_string(index) +
                                    " is out of bounds [0, " + std::to_string(points.size()) + ")");
        }

        // Set value and mark as determined
        parameters[fieldName][index] = value;
        points[index].setDetermined(fieldName, true);
    }
}

void FieldGenerator::initializeFirstPoint(const std::string& fieldName) {
    if (points.empty()) {
        throw std::runtime_error("Cannot initialize first point: field is empty");
    }

    // Generate random index
    const unsigned int startIndex = static_cast<unsigned int>(
        gsl_rng_uniform(rng) * points.size());

    // Generate unconditional Gaussian value
    const double unconditionalValue = gsl_ran_ugaussian(rng);

    // Set value and mark as determined
    parameters[fieldName][startIndex] = unconditionalValue;
    points[startIndex].setDetermined(fieldName, true);
}

unsigned int FieldGenerator::selectRandomUndeterminedPoint(const std::string& fieldName) const {
    std::vector<unsigned int> undeterminedIndices;
    undeterminedIndices.reserve(points.size());

    // Collect all undetermined point indices
    for (unsigned int i = 0; i < points.size(); ++i) {
        if (!points[i].isDetermined(fieldName)) {
            undeterminedIndices.push_back(i);
        }
    }

    if (undeterminedIndices.empty()) {
        throw std::runtime_error("No undetermined points available for selection");
    }

    // Select random undetermined point
    const size_t randomIndex = static_cast<size_t>(
        gsl_rng_uniform(rng) * undeterminedIndices.size());

    return undeterminedIndices[randomIndex];
}

bool FieldGenerator::allPointsDetermined(const std::string& fieldName) const {
    for (const auto& point : points) {
        if (!point.isDetermined(fieldName)) {
            return false;
        }
    }
    return true;
}

void FieldGenerator::generateConditionalValue(const std::string& fieldName,
                                              unsigned int targetIndex,
                                              double correlationLength) {
    validateTargetIndex(targetIndex, fieldName);

    // Build correlation matrix system for kriging
    const correl_mat_vec correlationData = buildCorrelationMatrix(fieldName, targetIndex, correlationLength);

    // Solve kriging system to get weights
    CVector_arma krigingWeights = solveKrigingSystem(correlationData);

    // Compute conditional mean: weights^T * observations
    double conditionalMean = dotproduct(krigingWeights,correlationData.V_RHS);

    // Compute conditional variance: 1 - weights^T * V_21
    double conditionalVariance = 1.0 - dotproduct(krigingWeights,correlationData.V_21);


    // Ensure variance is non-negative (numerical stability)
    conditionalVariance = std::max(0.0, conditionalVariance);

    // Generate conditional value
    const double standardGaussian = gsl_ran_ugaussian(rng);
    const double conditionalValue = conditionalMean + standardGaussian * std::sqrt(conditionalVariance);

    // Set value and mark as determined
    parameters[fieldName][targetIndex] = conditionalValue;
    points[targetIndex].setDetermined(fieldName, true);
}

std::vector<int> FieldGenerator::getDeterminedIndices(const std::string& fieldName) const {
    std::vector<int> determinedIndices;
    for (int i = 0; i < static_cast<int>(points.size()); ++i) {
        if (points[i].isDetermined(fieldName)) {
            determinedIndices.push_back(i);
        }
    }
    return determinedIndices;
}

correl_mat_vec FieldGenerator::buildCorrelationMatrix(const std::string& fieldName,
                                                      int targetIndex,
                                                      double correlationLength) const {
    // Get determined points
    const std::vector<int> determinedIndices = getDeterminedIndices(fieldName);
    const int numDetermined = static_cast<int>(determinedIndices.size());

    if (numDetermined == 0) {
        throw std::runtime_error("No determined points available for correlation matrix construction");
    }

    // Initialize correlation matrix and vectors
    correl_mat_vec correlationData;
    correlationData.M_22 = CMatrix_arma(numDetermined, numDetermined);
    correlationData.V_21 = CVector_arma(numDetermined);
    correlationData.V_RHS = CVector_arma(numDetermined);

    // Build correlation vectors and matrix
    buildCorrelationVectors(correlationData, determinedIndices, fieldName, targetIndex, correlationLength);
    buildCorrelationMatrix(correlationData, determinedIndices, correlationLength);

    return correlationData;
}

void FieldGenerator::buildCorrelationVectors(correl_mat_vec& correlationData,
                                             const std::vector<int>& determinedIndices,
                                             const std::string& fieldName,
                                             int targetIndex,
                                             double correlationLength) const {
    const int numDetermined = static_cast<int>(determinedIndices.size());

    for (int i = 0; i < numDetermined; ++i) {
        const int determinedIndex = determinedIndices[i];

        // V_21: correlation between target point and determined points
        correlationData.V_21[i] = computeCorrelation(targetIndex, determinedIndex, correlationLength);

        // V_RHS: observed values at determined points
        correlationData.V_RHS[i] = parameters.at(fieldName)[determinedIndex];
    }
}

void FieldGenerator::buildCorrelationMatrix(correl_mat_vec& correlationData,
                                            const std::vector<int>& determinedIndices,
                                            double correlationLength) const {
    const int numDetermined = static_cast<int>(determinedIndices.size());

    // Build M_22: correlation matrix between determined points
    for (int i = 0; i < numDetermined; ++i) {
        for (int j = 0; j < numDetermined; ++j) {
            correlationData.M_22(i,j) = computeCorrelation(determinedIndices[i], determinedIndices[j], correlationLength);
        }
    }
}

CVector_arma FieldGenerator::solveKrigingSystem(const correl_mat_vec& correlationData) const {

    CVector_arma weights;

    // Simple solution: M_22 * weights = V_21
    // For now, use Gaussian elimination (you can replace with more robust solver)

    // Copy M_22 and V_21 for solving
    weights = CMatrix_arma(inv(correlationData.M_22))*correlationData.V_21;

    return weights;
}

double FieldGenerator::computeCorrelation(int i, int j, double correlationLength) const {
    if (correlationLength <= 0.0) {
        throw std::invalid_argument("Correlation length must be positive");
    }

    const double distance = std::abs(i - j) * dx;
    return std::exp(-distance / correlationLength);
}

void FieldGenerator::validateFieldGeneration() const {
    if (rng == nullptr) {
        throw std::runtime_error("Random number generator not initialized");
    }

    if (dx <= 0.0) {
        throw std::invalid_argument("Grid spacing (dx) must be positive");
    }
}

void FieldGenerator::validateTargetIndex(unsigned int targetIndex, const std::string& fieldName) const {
    if (targetIndex >= points.size()) {
        throw std::out_of_range("Target index " + std::to_string(targetIndex) +
                                " is out of bounds (size: " + std::to_string(points.size()) + ")");
    }

    if (rng == nullptr) {
        throw std::runtime_error("Random number generator not initialized");
    }
}

// Helper method to check if any points are determined
bool FieldGenerator::hasAnyDeterminedPoints(const std::string& fieldName) const {
    for (const auto& point : points) {
        if (point.isDetermined(fieldName)) {
            return true;
        }
    }
    return false;
}

// Calculate mean of field values
double FieldGenerator::mean(const std::string& fieldName) const {
    auto it = parameters.find(fieldName);
    if (it == parameters.end()) {
        throw std::invalid_argument("Field '" + fieldName + "' does not exist");
    }

    const std::vector<double>& fieldData = it->second;
    if (fieldData.empty()) {
        throw std::runtime_error("Field '" + fieldName + "' is empty");
    }

    double sum = std::accumulate(fieldData.begin(), fieldData.end(), 0.0);
    return sum / fieldData.size();
}

// Calculate standard deviation of field values
double FieldGenerator::standardDeviation(const std::string& fieldName) const {
    auto it = parameters.find(fieldName);
    if (it == parameters.end()) {
        throw std::invalid_argument("Field '" + fieldName + "' does not exist");
    }

    const std::vector<double>& fieldData = it->second;
    if (fieldData.empty()) {
        throw std::runtime_error("Field '" + fieldName + "' is empty");
    }

    if (fieldData.size() == 1) {
        return 0.0;  // Standard deviation of single value is 0
    }

    double fieldMean = mean(fieldName);
    double sumSquaredDiffs = 0.0;

    for (double value : fieldData) {
        double diff = value - fieldMean;
        sumSquaredDiffs += diff * diff;
    }

    return std::sqrt(sumSquaredDiffs / (fieldData.size() - 1));  // Sample standard deviation
}

// Calculate mean of log-transformed field values
double FieldGenerator::meanLog(const std::string& fieldName) const {
    auto it = parameters.find(fieldName);
    if (it == parameters.end()) {
        throw std::invalid_argument("Field '" + fieldName + "' does not exist");
    }

    const std::vector<double>& fieldData = it->second;
    if (fieldData.empty()) {
        throw std::runtime_error("Field '" + fieldName + "' is empty");
    }

    double sumLogs = 0.0;
    for (double value : fieldData) {
        if (value <= 0.0) {
            throw std::domain_error("Cannot compute log of non-positive value: " +
                                    std::to_string(value) + " in field '" + fieldName + "'");
        }
        sumLogs += std::log(value);
    }

    return sumLogs / fieldData.size();
}

// Calculate standard deviation of log-transformed field values
double FieldGenerator::standardDeviationLog(const std::string& fieldName) const {
    auto it = parameters.find(fieldName);
    if (it == parameters.end()) {
        throw std::invalid_argument("Field '" + fieldName + "' does not exist");
    }

    const std::vector<double>& fieldData = it->second;
    if (fieldData.empty()) {
        throw std::runtime_error("Field '" + fieldName + "' is empty");
    }

    if (fieldData.size() == 1) {
        return 0.0;  // Standard deviation of single value is 0
    }

    // Validate all values are positive
    for (double value : fieldData) {
        if (value <= 0.0) {
            throw std::domain_error("Cannot compute log of non-positive value: " +
                                    std::to_string(value) + " in field '" + fieldName + "'");
        }
    }

    double logMean = meanLog(fieldName);
    double sumSquaredDiffs = 0.0;

    for (double value : fieldData) {
        double logValue = std::log(value);
        double diff = logValue - logMean;
        sumSquaredDiffs += diff * diff;
    }

    return std::sqrt(sumSquaredDiffs / (fieldData.size() - 1));  // Sample standard deviation
}

// Calculate auto-correlation at a given distance using interpolation
double FieldGenerator::autoCorrelation(const std::string& fieldName, double distance) const {
    auto it = parameters.find(fieldName);
    if (it == parameters.end()) {
        throw std::invalid_argument("Field '" + fieldName + "' does not exist");
    }

    const std::vector<double>& fieldData = it->second;
    if (fieldData.empty()) {
        throw std::runtime_error("Field '" + fieldName + "' is empty");
    }

    if (fieldData.size() < 2) {
        return 1.0;  // Single point field has correlation 1
    }

    if (std::abs(distance) < 1e-12) {
        return 1.0;  // Auto-correlation at zero distance is always 1
    }

    // Convert distance to lag (continuous)
    double continuousLag = std::abs(distance) / dx;

    // Check if distance is too large
    if (continuousLag >= static_cast<double>(fieldData.size() - 1)) {
        return 0.0;  // Beyond field extent, assume zero correlation
    }

    const int n = static_cast<int>(fieldData.size());

    // Calculate field mean and variance once
    double fieldMean = mean(fieldName);
    double fieldVariance = 0.0;

    for (double value : fieldData) {
        double diff = value - fieldMean;
        fieldVariance += diff * diff;
    }
    fieldVariance /= n;  // Population variance

    if (fieldVariance < 1e-12) {
        return 1.0;  // Constant field has correlation 1
    }

    // Find integer lag bounds
    int lowerLag = static_cast<int>(std::floor(continuousLag));
    int upperLag = lowerLag + 1;

    // Calculate auto-correlation at both integer lags
    auto calculateCorrelationAtLag = [&](int lag) -> double {
        if (lag == 0) return 1.0;
        if (lag >= n) return 0.0;

        int validPairs = n - lag;
        if (validPairs <= 0) return 0.0;

        double autoCovariance = 0.0;
        for (int i = 0; i < validPairs; ++i) {
            int j = i + lag;
            double val1 = fieldData[i] - fieldMean;
            double val2 = fieldData[j] - fieldMean;
            autoCovariance += val1 * val2;
        }

        autoCovariance /= validPairs;
        return autoCovariance / fieldVariance;
    };

    // Get correlations at integer lag bounds
    double corrLower = calculateCorrelationAtLag(lowerLag);
    double corrUpper = calculateCorrelationAtLag(upperLag);

    // Linear interpolation between the two correlations
    double weight = continuousLag - lowerLag;  // Weight for upper lag (0 to 1)
    double interpolatedCorrelation = (1.0 - weight) * corrLower + weight * corrUpper;

    return interpolatedCorrelation;
}
// Piecewise linear interpolation to estimate field value at any x position
double FieldGenerator::interpolate(const std::string& fieldName, double x) const {
    auto it = parameters.find(fieldName);
    if (it == parameters.end()) {
        throw std::invalid_argument("Field '" + fieldName + "' does not exist");
    }

    const std::vector<double>& fieldData = it->second;
    if (fieldData.empty()) {
        throw std::runtime_error("Field '" + fieldName + "' is empty");
    }

    if (fieldData.size() == 1) {
        return fieldData[0];  // Single point - return that value
    }

    // Convert position to grid index (continuous)
    double continuousIndex = x / dx;

    // Handle boundary cases
    if (continuousIndex <= 0.0) {
        return fieldData[0];  // Extrapolate with first value
    }

    if (continuousIndex >= static_cast<double>(fieldData.size() - 1)) {
        return fieldData.back();  // Extrapolate with last value
    }

    // Find the two surrounding grid points
    int leftIndex = static_cast<int>(std::floor(continuousIndex));
    int rightIndex = leftIndex + 1;

    // Ensure indices are valid (should be guaranteed by boundary checks above)
    leftIndex = std::max(0, std::min(leftIndex, static_cast<int>(fieldData.size() - 1)));
    rightIndex = std::max(0, std::min(rightIndex, static_cast<int>(fieldData.size() - 1)));

    // Get values at surrounding points
    double leftValue = fieldData[leftIndex];
    double rightValue = fieldData[rightIndex];

    // Calculate interpolation weight
    double weight = continuousIndex - leftIndex;  // Weight for right point (0 to 1)

    // Linear interpolation: value = (1-w) * left + w * right
    double interpolatedValue = (1.0 - weight) * leftValue + weight * rightValue;

    return interpolatedValue;
}

// Transform normal scores to uniform distribution [0,1]
void FieldGenerator::normalToUniform(const std::string& sourceField, const std::string& targetField) {
    // Check if source field exists
    auto sourceIt = parameters.find(sourceField);
    if (sourceIt == parameters.end()) {
        throw std::invalid_argument("Source field '" + sourceField + "' does not exist");
    }

    const std::vector<double>& sourceData = sourceIt->second;
    if (sourceData.empty()) {
        throw std::runtime_error("Source field '" + sourceField + "' is empty");
    }

    // Create or resize target field
    if (parameters.find(targetField) == parameters.end()) {
        parameters[targetField] = std::vector<double>(sourceData.size());
    } else if (parameters[targetField].size() != sourceData.size()) {
        parameters[targetField].resize(sourceData.size());
    }

    std::vector<double>& targetData = parameters[targetField];

    // Transform each normal score to uniform using standard normal CDF
    for (size_t i = 0; i < sourceData.size(); ++i) {
        // Use GSL's standard normal CDF: Φ(z) maps N(0,1) to U(0,1)
        targetData[i] = gsl_cdf_ugaussian_P(sourceData[i]);

        // Ensure values are strictly within [0,1] (handle numerical precision)
        targetData[i] = std::max(0.0, std::min(1.0, targetData[i]));
    }

    // Mark all points as determined for the target field
    for (size_t i = 0; i < points.size(); ++i) {
        points[i].setDetermined(targetField, true);
    }
}

void FieldGenerator::normalToCDF(const std::string& sourceField, const std::string& targetField, const TimeSeries<double> &CDF)
{
    // Check if source field exists
    auto sourceIt = parameters.find(sourceField);
    if (sourceIt == parameters.end()) {
        throw std::invalid_argument("Source field '" + sourceField + "' does not exist");
    }

    const std::vector<double>& sourceData = sourceIt->second;
    if (sourceData.empty()) {
        throw std::runtime_error("Source field '" + sourceField + "' is empty");
    }

    // Create or resize target field
    if (parameters.find(targetField) == parameters.end()) {
        parameters[targetField] = std::vector<double>(sourceData.size());
    } else if (parameters[targetField].size() != sourceData.size()) {
        parameters[targetField].resize(sourceData.size());
    }

    std::vector<double>& targetData = parameters[targetField];

    // Transform each normal score to uniform using standard normal CDF
    for (size_t i = 0; i < sourceData.size(); ++i) {
        // Use GSL's standard normal CDF: Φ(z) maps N(0,1) to U(0,1)
        targetData[i] = gsl_cdf_ugaussian_P(sourceData[i]);

        // Ensure values are strictly within [0,1] (handle numerical precision)
        targetData[i] = CDF.inverse_CDF(std::max(0.0, std::min(1.0, targetData[i])));
    }

    // Mark all points as determined for the target field
    for (size_t i = 0; i < points.size(); ++i) {
        points[i].setDetermined(targetField, true);
    }
}

void FieldGenerator::normalToCDF(const std::string& sourceField, const std::string& targetField, const double &mean_log, const double &std_log)
{
    // Check if source field exists
    auto sourceIt = parameters.find(sourceField);
    if (sourceIt == parameters.end()) {
        throw std::invalid_argument("Source field '" + sourceField + "' does not exist");
    }

    const std::vector<double>& sourceData = sourceIt->second;
    if (sourceData.empty()) {
        throw std::runtime_error("Source field '" + sourceField + "' is empty");
    }

    // Create or resize target field
    if (parameters.find(targetField) == parameters.end()) {
        parameters[targetField] = std::vector<double>(sourceData.size());
    } else if (parameters[targetField].size() != sourceData.size()) {
        parameters[targetField].resize(sourceData.size());
    }

    std::vector<double>& targetData = parameters[targetField];

    // Transform each normal score to uniform using standard normal CDF
    for (size_t i = 0; i < sourceData.size(); ++i) {
        // map to lognormal distribution
        targetData[i] = exp(sourceData[i]*std_log + mean_log);
    }

    // Mark all points as determined for the target field
    for (size_t i = 0; i < points.size(); ++i) {
        points[i].setDetermined(targetField, true);
    }
}


// Transform field using exponential function: target = exp(a*source + b)
void FieldGenerator::exponentialTransform(const std::string& sourceField, const std::string& targetField,
                                          double a, double b) {
    // Check if source field exists
    auto sourceIt = parameters.find(sourceField);
    if (sourceIt == parameters.end()) {
        throw std::invalid_argument("Source field '" + sourceField + "' does not exist");
    }

    const std::vector<double>& sourceData = sourceIt->second;
    if (sourceData.empty()) {
        throw std::runtime_error("Source field '" + sourceField + "' is empty");
    }

    // Create or resize target field
    if (parameters.find(targetField) == parameters.end()) {
        parameters[targetField] = std::vector<double>(sourceData.size());
    } else if (parameters[targetField].size() != sourceData.size()) {
        parameters[targetField].resize(sourceData.size());
    }

    std::vector<double>& targetData = parameters[targetField];

    // Apply exponential transformation: target[i] = exp(a * source[i] + b)
    for (size_t i = 0; i < sourceData.size(); ++i) {
        double exponent = a * sourceData[i] + b;

        // Check for potential overflow before computing exp
        if (exponent > 700.0) {  // exp(700) ≈ 1e304, near double max
            throw std::runtime_error("Exponential transformation would cause overflow at index " +
                                     std::to_string(i) + ": exp(" + std::to_string(exponent) + ")");
        }

        targetData[i] = std::exp(exponent);

        // Additional check for resulting infinity or NaN
        if (!std::isfinite(targetData[i])) {
            throw std::runtime_error("Exponential transformation resulted in non-finite value at index " +
                                     std::to_string(i) + ": exp(" + std::to_string(exponent) + ") = " +
                                     std::to_string(targetData[i]));
        }
    }

    // Mark all points as determined for the target field
    for (size_t i = 0; i < points.size(); ++i) {
        points[i].setDetermined(targetField, true);
    }
}

// Find minimum value in field
double FieldGenerator::minimum(const std::string& fieldName) const {
    auto it = parameters.find(fieldName);
    if (it == parameters.end()) {
        throw std::invalid_argument("Field '" + fieldName + "' does not exist");
    }

    const std::vector<double>& fieldData = it->second;
    if (fieldData.empty()) {
        throw std::runtime_error("Field '" + fieldName + "' is empty");
    }

    // Use std::min_element to find minimum value
    auto minIt = std::min_element(fieldData.begin(), fieldData.end());
    return *minIt;
}

// Find maximum value in field
double FieldGenerator::maximum(const std::string& fieldName) const {
    auto it = parameters.find(fieldName);
    if (it == parameters.end()) {
        throw std::invalid_argument("Field '" + fieldName + "' does not exist");
    }

    const std::vector<double>& fieldData = it->second;
    if (fieldData.empty()) {
        throw std::runtime_error("Field '" + fieldName + "' is empty");
    }

    // Use std::max_element to find maximum value
    auto maxIt = std::max_element(fieldData.begin(), fieldData.end());
    return *maxIt;
}

// Interpolate field value at any x position using linear interpolation
double FieldGenerator::interpolateAt(double x, const std::string& quantity) const {
    validateFieldForInterpolation(quantity);

    auto it = parameters.find(quantity);
    if (it == parameters.end()) {
        throw std::invalid_argument("Field '" + quantity + "' does not exist");
    }

    const std::vector<double>& fieldData = it->second;
    if (fieldData.empty()) {
        throw std::runtime_error("Field '" + quantity + "' is empty");
    }

    if (fieldData.size() == 1) {
        return fieldData[0];
    }

    // Find the grid position
    const double gridPosition = x / dx;

    // Handle boundary cases
    if (gridPosition <= 0.0) {
        return fieldData[0];
    }
    if (gridPosition >= static_cast<double>(fieldData.size() - 1)) {
        return fieldData.back();
    }

    // Find surrounding grid points
    const unsigned int lowerIndex = static_cast<unsigned int>(std::floor(gridPosition));
    const unsigned int upperIndex = lowerIndex + 1;

    // Linear interpolation weights
    const double fraction = gridPosition - static_cast<double>(lowerIndex);
    const double lowerValue = fieldData[lowerIndex];
    const double upperValue = fieldData[upperIndex];

    return lowerValue + fraction * (upperValue - lowerValue);
}

// Interpolate multiple points at once
std::vector<double> FieldGenerator::interpolateAt(const std::vector<double>& xPositions,
                                                  const std::string& quantity) const {
    std::vector<double> results;
    results.reserve(xPositions.size());

    for (double x : xPositions) {
        results.push_back(interpolateAt(x, quantity));
    }

    return results;
}

// Get the x-coordinate for a given grid index
double FieldGenerator::indexToPosition(unsigned int index) const {
    if (index >= points.size()) {
        throw std::out_of_range("Index out of bounds for position calculation");
    }
    return static_cast<double>(index) * dx;
}

// Get the grid index closest to a given x position
unsigned int FieldGenerator::positionToIndex(double x) const {
    if (x < 0.0) {
        return 0;
    }

    const double indexDouble = x / dx;
    const unsigned int index = static_cast<unsigned int>(std::round(indexDouble));

    return std::min(index, static_cast<unsigned int>(points.size() - 1));
}

// Validation method for interpolation (add this to the private section)
void FieldGenerator::validateFieldForInterpolation(const std::string& quantity) const {
    if (points.empty()) {
        throw std::runtime_error("Cannot interpolate: field is empty");
    }

    if (dx <= 0.0) {
        throw std::invalid_argument("Grid spacing (dx) must be positive for interpolation");
    }

    auto it = parameters.find(quantity);
    if (it == parameters.end()) {
        throw std::invalid_argument("Field '" + quantity + "' does not exist");
    }
}

// Implementation:
unsigned int FieldGenerator::findClosestNode(double x) const {
    if (points.empty()) {
        throw std::runtime_error("Cannot find closest node: field is empty");
    }

    if (dx <= 0.0) {
        throw std::invalid_argument("Grid spacing (dx) must be positive");
    }

    // Convert x position to continuous grid index
    const double continuousIndex = x / dx;

    // Round to nearest integer index
    int nearestIndex = static_cast<int>(std::round(continuousIndex));

    // Clamp to valid range [0, points.size()-1]
    nearestIndex = std::max(0, std::min(nearestIndex, static_cast<int>(points.size() - 1)));

    return static_cast<unsigned int>(nearestIndex);
}

void FieldGenerator::setKnownPointsFromTimeSeries(const std::string& fieldName,
                                                  const TimeSeries<double>& timeSeries,
                                                  double correlationLength) {
    // Validate inputs
    if (points.empty()) {
        throw std::runtime_error("Cannot set known points: field is empty");
    }

    if (timeSeries.size() == 0) {
        throw std::invalid_argument("TimeSeries is empty");
    }

    if (dx <= 0.0) {
        throw std::invalid_argument("Grid spacing (dx) must be positive");
    }

    // Create or resize parameter field
    if (parameters.find(fieldName) == parameters.end()) {
        parameters[fieldName] = std::vector<double>(points.size(), 0.0);
    } else if (parameters[fieldName].size() != points.size()) {
        parameters[fieldName].resize(points.size(), 0.0);
    }

    // Reset determination flags for this field
    for (auto& point : points) {
        point.setDetermined(fieldName, false);
    }

    // Convert TimeSeries data to known points vector
    std::vector<std::pair<int, double>> knownPoints;
    knownPoints.reserve(timeSeries.size());

    for (int i = 0; i < timeSeries.size(); ++i) {
        double xPosition = timeSeries.getTime(i);  // t represents x position
        double value = timeSeries.getValue(i);      // represents field value

        // Find closest grid node to this x position
        unsigned int closestIndex = findClosestNode(xPosition);

        // Add to known points (avoiding duplicates by using the closest node)
        knownPoints.emplace_back(static_cast<int>(closestIndex), value);
    }

    // Remove duplicate indices (keep the last value if multiple TimeSeries points map to same node)
    std::map<int, double> uniqueKnownPoints;
    for (const auto& knownPoint : knownPoints) {
        uniqueKnownPoints[knownPoint.first] = knownPoint.second;
    }

    // Convert back to vector format
    std::vector<std::pair<int, double>> finalKnownPoints;
    finalKnownPoints.reserve(uniqueKnownPoints.size());
    for (const auto& entry : uniqueKnownPoints) {
        finalKnownPoints.emplace_back(entry.first, entry.second);
    }

    // Generate the field using the known points
    generateNormalScoreField(fieldName, correlationLength, finalKnownPoints);
}

bool FieldGenerator::writeFieldToCSV(const std::string& fieldName, const std::string& filename,
                                     int precision, bool includeHeader) const {
    // Check if field exists
    auto it = parameters.find(fieldName);
    if (it == parameters.end()) {
        return false; // Field doesn't exist
    }

    const std::vector<double>& fieldData = it->second;

    // Check if field has data
    if (fieldData.empty()) {
        return false; // No data to write
    }

    // Check if field size matches points size
    if (fieldData.size() != points.size()) {
        return false; // Inconsistent data
    }

    // Open file for writing
    std::ofstream file(filename);
    if (!file.is_open()) {
        return false; // Failed to open file
    }

    // Set precision for floating point output
    file << std::fixed << std::setprecision(precision);

    // Write CSV header if requested
    if (includeHeader) {
        file << "x," << fieldName << std::endl;
    }

    // Write data in CSV format
    for (size_t i = 0; i < fieldData.size(); ++i) {
        double xPosition = static_cast<double>(i) * dx;
        file << xPosition << "," << fieldData[i] << std::endl;
    }

    file.close();
    return true;
}

const TimeSeriesSet<double>& FieldGenerator::getMeasuredCDFs() const {
    if (!has_measured_CDFs_) {
        throw std::runtime_error("No measured CDFs available. Call generateFieldFromMeasurementData first.");
    }
    return measured_CDFs_;
}

bool FieldGenerator::hasMeasuredCDFs() const {
    return has_measured_CDFs_;
}

void FieldGenerator::clearMeasuredCDFs() {
    measured_CDFs_.clear();
    has_measured_CDFs_ = false;
}

bool FieldGenerator::generateFieldFromMeasurementData(const std::string& dataFilePath,
                                                      const std::string& parameterName,
                                                      const std::string& fieldName,
                                                      double correlationLength,
                                                      bool verbose) {
    try {
        if (verbose) {
            std::cout << "Reading data from '" << dataFilePath << "'" << std::endl;
        }

        // Step 1: Read measurement data
        TimeSeriesSet<double> measured_data;
        if (!measured_data.read(dataFilePath)) {
            if (verbose) {
                std::cerr << "Error: Failed to read data from " << dataFilePath << std::endl;
            }
            return false;
        }

        // Step 2: Check if parameter exists in the data
        if (!measured_data.Contains(parameterName)) {
            if (verbose) {
                std::cerr << "Error: Parameter '" << parameterName << "' not found in data file" << std::endl;
            }
            return false;
        }

        // Step 3: Generate cumulative distribution functions and store them
        CVector mean_log_values;
        CVector std_log_values;
        if (pdfmod_ == pdfmode::nonparametric)
        {   measured_CDFs_ = measured_data.GetCummulativeDistribution();
            has_measured_CDFs_ = true;
        }
        else
        {
            mean_log_values = measured_data.Log().mean();
            std_log_values = measured_data.Log().standardDeviation();
            measured_CDFs_ = TimeSeriesSet<double>::LogNormalCDF(mean_log_values, std_log_values,100);
            measured_CDFs_.setSeriesNames(measured_data.getSeriesNames());
            has_measured_CDFs_ = true;
        }

        // Step 4: Convert to normal scores
        TimeSeriesSet<double> normal_scores = measured_data.ConverttoNormalScore();

        // Step 5: Check if the converted parameter exists
        if (!normal_scores.Contains(parameterName)) {
            if (verbose) {
                std::cerr << "Error: Failed to convert parameter '" << parameterName << "' to normal scores" << std::endl;
            }
            return false;
        }

        // Step 6: Set known points and generate field
        setKnownPointsFromTimeSeries(fieldName, normal_scores[parameterName], correlationLength);

        // Step 7: Transform normal scores back to real values using measured CDF
        if (!measured_CDFs_.Contains(parameterName)) {
            if (verbose) {
                std::cerr << "Error: CDF for parameter '" << parameterName << "' not found" << std::endl;
            }
            return false;
        }

        normalToCDF(fieldName, parameterName, measured_CDFs_[parameterName]);

        if (verbose) {
            std::cout << "Successfully generated field '" << fieldName << "' from parameter '"
                      << parameterName << "'" << std::endl;
            std::cout << "Transformed to real values in field '" << parameterName << "'" << std::endl;
        }

        return true;

    } catch (const std::exception& e) {
        if (verbose) {
            std::cerr << "Error in generateFieldFromMeasurementData: " << e.what() << std::endl;
        }
        return false;
    }
}
