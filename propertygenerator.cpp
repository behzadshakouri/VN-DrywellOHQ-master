#include "propertygenerator.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_cdf.h"
#ifndef  _WINDOWS
#include <sys/time.h>
#endif // ! _WINDOWS


PropertyGenerator::PropertyGenerator(int seed):std::vector<Propfull>()
{
#ifndef _WINDOWS
    struct timeval tv;
    gettimeofday(&tv, 0);
    unsigned long mySeed = tv.tv_sec + tv.tv_usec+seed;
#else
    long int mySeed = static_cast<long int>(time(NULL)+seed);
#endif // !_WINDOWS
gsl_rng_set(r, mySeed);
}


PropertyGenerator::PropertyGenerator(unsigned int n,int seed):std::vector<Propfull>(n)
{
    initializeRNG(seed);

}


void PropertyGenerator::initializeRNG(int seed)
{
    // Initialize RNG type and allocate
    initializeRNG(seed);
}

PropertyGenerator::~PropertyGenerator()
{
    if (r != nullptr) {
        gsl_rng_free(r);
        r = nullptr;
    }
}

unsigned long PropertyGenerator::generateSeed(int seed) const
{
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

correl_mat_vec PropertyGenerator::get_correll_matrix_vec(int i)
{
    correl_mat_vec Correl_Matrix_Vector;
    int num_determined = GetNumberOfPointsDetermined();
    std::vector<int> determined = Determined();
#ifdef _arma
    Correl_Matrix_Vector.M_22 = CMatrix_arma(num_determined);
    Correl_Matrix_Vector.V_21 = CVector_arma(num_determined);
    Correl_Matrix_Vector.V_RHS = CVector_arma(num_determined);
#else
    Correl_Matrix_Vector.M_22 = CMatrix(num_determined);
    Correl_Matrix_Vector.V_21 = CVector(num_determined);
    Correl_Matrix_Vector.V_RHS = CVector(num_determined);
#endif //  arma
    for (int ii = 0; ii < num_determined; ii++)
    {
        Correl_Matrix_Vector.V_21[ii] = exp(-fabs(determined[ii]-i)*dx/correlation_length_scale);
        Correl_Matrix_Vector.V_RHS[ii] = at(determined[ii]).normal_scores.K_sat;
        for (int jj = 0; jj < num_determined; jj++)
        {
#ifdef _arma
            Correl_Matrix_Vector.M_22(ii,jj) = exp(-fabs(determined[ii]-determined[jj])*dx/correlation_length_scale);
#else
            Correl_Matrix_Vector.M_22[ii][jj] = exp(-fabs(determined[ii]-determined[jj])*dx/correlation_length_scale);
#endif // arma
        }
    }
    return Correl_Matrix_Vector;
}

correl_mat_vec PropertyGenerator::get_correll_matrix_vec(int i)
{
    return buildCorrelationMatrix(i);
}

double PropertyGenerator::computeCorrelation(int i, int j) const
{
    if (correlation_length_scale <= 0.0) {
        throw std::invalid_argument("Correlation length scale must be positive");
    }

    const double distance = std::abs(i - j) * dx;
    return std::exp(-distance / correlation_length_scale);
}

void PropertyGenerator::validateCorrelationInputs(int targetIndex) const
{
    if (targetIndex < 0 || targetIndex >= static_cast<int>(size())) {
        throw std::out_of_range("Target index out of bounds");
    }

    if (dx <= 0.0) {
        throw std::invalid_argument("Grid spacing (dx) must be positive");
    }

    if (correlation_length_scale <= 0.0) {
        throw std::invalid_argument("Correlation length scale must be positive");
    }

    if (GetNumberOfPointsDetermined() == 0) {
        throw std::runtime_error("No points have been determined yet");
    }
}

void PropertyGenerator::assign_Gaussian_scores(unsigned int targetIndex)
{
    validateTargetIndex(targetIndex);

    if (at(targetIndex).k_det) {
        return; // Already determined, nothing to do
    }

    // Build correlation matrix and vectors for conditional simulation
    const correl_mat_vec correlationData = buildCorrelationMatrix(static_cast<int>(targetIndex));

    // Compute conditional Gaussian parameters (kriging)
    const ConditionalGaussianParams params = computeConditionalGaussianParams(correlationData);

    // Generate random value from conditional distribution
    const double standardGaussian = gsl_ran_ugaussian(r);
    const double conditionalValue = params.mean + standardGaussian * std::sqrt(params.variance);

    // Update the target point
    at(targetIndex).normal_scores.K_sat = conditionalValue;
    at(targetIndex).k_det = true;
}

PropertyGenerator::ConditionalGaussianParams
PropertyGenerator::computeConditionalGaussianParams(const correl_mat_vec& correlationData) const
{
    ConditionalGaussianParams params;

    try {
        // Solve the kriging system: M_22 * weights = V_21
#ifdef _arma
        const CMatrix_arma matrixInverse = inv(correlationData.M_22);
#else
        const CMatrix matrixInverse = Invert(correlationData.M_22);
#endif

        // Compute kriging weights
        const CVector_arma krigingWeights = matrixInverse * correlationData.V_21;

        // Conditional mean: weights^T * observations
        params.mean = dotproduct(krigingWeights, correlationData.V_RHS);

        // Conditional variance: 1 - weights^T * V_21
        params.variance = 1.0 - dotproduct(krigingWeights, correlationData.V_21);

        // Ensure variance is non-negative (numerical stability)
        if (params.variance < 0.0) {
            params.variance = 0.0;
        }

    } catch (const std::exception& e) {
        throw std::runtime_error(
            std::string("Failed to compute conditional Gaussian parameters: ") + e.what());
    }

    return params;
}

void PropertyGenerator::validateTargetIndex(unsigned int targetIndex) const
{
    if (targetIndex >= size()) {
        throw std::out_of_range(
            "Target index " + std::to_string(targetIndex) +
            " is out of bounds (size: " + std::to_string(size()) + ")");
    }

    if (GetNumberOfPointsDetermined() == 0) {
        throw std::runtime_error(
            "Cannot assign conditional K_sat value: no points have been determined yet");
    }

    if (r == nullptr) {
        throw std::runtime_error("Random number generator not initialized");
    }
}


void PropertyGenerator::generateField()
{
    validateFieldGeneration();

    if (size() == 0) {
        return; // Nothing to generate
    }

    // Initialize the first point with unconditional Gaussian value
    initializeFirstPoint();

    // Sequential simulation: generate remaining points conditionally
    while (!allPointsDetermined()) {
        const unsigned int targetIndex = selectRandomUndeterminedPoint();
        assign_Gaussian_scores(targetIndex);
    }
}

void PropertyGenerator::initializeFirstPoint(InitializationStrategy strategy)
{
    if (size() == 0) {
        throw std::runtime_error("Cannot initialize first point: field is empty");
    }

    switch (strategy) {
    case InitializationStrategy::RANDOM_LOCATION:
        initializeAtRandomLocation();
        break;
    case InitializationStrategy::CENTER_POINT:
        initializeAtCenter();
        break;
    case InitializationStrategy::MULTIPLE_SEEDS:
        initializeMultipleSeeds();
        break;
    default:
        throw std::invalid_argument("Unknown initialization strategy");
    }
}

void PropertyGenerator::initializeAtRandomLocation()
{
    const unsigned int startIndex = generateRandomIndex();
    const double unconditionalValue = gsl_ran_ugaussian(r);

    setPointValue(startIndex, unconditionalValue);
}

void PropertyGenerator::initializeAtCenter()
{
    const unsigned int centerIndex = static_cast<unsigned int>(size() / 2);
    const double unconditionalValue = gsl_ran_ugaussian(r);

    setPointValue(centerIndex, unconditionalValue);
}

void PropertyGenerator::initializeMultipleSeeds(unsigned int numSeeds)
{
    if (numSeeds == 0) {
        throw std::invalid_argument("Number of seeds must be positive");
    }

    // Limit seeds to available points
    const unsigned int maxSeeds = std::min(numSeeds, static_cast<unsigned int>(size()));

    // Generate well-distributed seed points
    std::vector<unsigned int> seedIndices;
    if (maxSeeds == 1) {
        seedIndices.push_back(generateRandomIndex());
    } else {
        // Distribute seeds evenly across the domain
        const double spacing = static_cast<double>(size()) / maxSeeds;
        for (unsigned int i = 0; i < maxSeeds; ++i) {
            const unsigned int baseIndex = static_cast<unsigned int>(i * spacing);
            // Add small random perturbation
            const unsigned int perturbation = static_cast<unsigned int>(
                gsl_rng_uniform(r) * std::min(spacing * 0.3, 10.0));
            const unsigned int seedIndex = std::min(baseIndex + perturbation,
                                                    static_cast<unsigned int>(size() - 1));
            seedIndices.push_back(seedIndex);
        }
    }

    // Initialize all seed points
    for (unsigned int index : seedIndices) {
        const double unconditionalValue = gsl_ran_ugaussian(r);
        setPointValue(index, unconditionalValue);
    }
}

unsigned int PropertyGenerator::generateRandomIndex() const
{
    if (size() == 0) {
        throw std::runtime_error("Cannot generate random index: field is empty");
    }

    // Generate random index in range [0, size-1]
    // Using proper scaling to avoid bias
    const double uniformRandom = gsl_rng_uniform(r);  // [0, 1)
    const unsigned int index = static_cast<unsigned int>(uniformRandom * size());

    // Ensure we don't exceed bounds due to floating point precision
    return std::min(index, static_cast<unsigned int>(size() - 1));
}

void PropertyGenerator::setPointValue(unsigned int index, double value)
{
    validatePointForInitialization(index);

    at(index).normal_scores.K_sat = value;
    at(index).k_det = true;
}

void PropertyGenerator::validatePointForInitialization(unsigned int index) const
{
    if (index >= size()) {
        throw std::out_of_range(
            "Point index " + std::to_string(index) +
            " is out of bounds (size: " + std::to_string(size()) + ")");
    }

    if (at(index).k_det) {
        throw std::runtime_error(
            "Point at index " + std::to_string(index) + " is already determined");
    }
}

// Overloaded version for backward compatibility and convenience
void PropertyGenerator::initializeFirstPoint()
{
    initializeFirstPoint(InitializationStrategy::RANDOM_LOCATION);
}

// Additional utility method to check if any points are initialized
bool PropertyGenerator::hasInitializedPoints() const
{
    for (unsigned int i = 0; i < size(); ++i) {
        if (at(i).k_det) {
            return true;
        }
    }
    return false;
}

// Method to reset all points (useful for re-running simulations)
void PropertyGenerator::resetAllPoints()
{
    for (unsigned int i = 0; i < size(); ++i) {
        at(i).k_det = false;
        at(i).normal_scores.K_sat = 0.0;
        // Could also reset other properties if needed
    }
}

void PropertyGenerator::initializeFirstPoint()
{
    // Select random starting point
    const unsigned int startIndex = static_cast<unsigned int>(
        gsl_rng_uniform(r) * size());

    // Generate unconditional Gaussian value
    const double unconditionalValue = gsl_ran_ugaussian(r);

    // Set the first point
    at(startIndex).normal_scores.K_sat = unconditionalValue;
    at(startIndex).k_det = true;
}

unsigned int PropertyGenerator::selectRandomUndeterminedPoint() const
{
    std::vector<unsigned int> undeterminedIndices;
    undeterminedIndices.reserve(size());

    // Collect all undetermined point indices
    for (unsigned int i = 0; i < size(); ++i) {
        if (!at(i).k_det) {
            undeterminedIndices.push_back(i);
        }
    }

    if (undeterminedIndices.empty()) {
        throw std::runtime_error("No undetermined points available for selection");
    }

    // Select random undetermined point
    const size_t randomIndex = static_cast<size_t>(
        gsl_rng_uniform(r) * undeterminedIndices.size());

    return undeterminedIndices[randomIndex];
}

bool PropertyGenerator::allPointsDetermined() const
{
    for (unsigned int i = 0; i < size(); ++i) {
        if (!at(i).k_det) {
            return false;
        }
    }
    return true;
}

void PropertyGenerator::validateFieldGeneration() const
{
    if (r == nullptr) {
        throw std::runtime_error("Random number generator not initialized");
    }

    if (dx <= 0.0) {
        throw std::invalid_argument("Grid spacing (dx) must be positive");
    }

    if (correlation_length_scale <= 0.0) {
        throw std::invalid_argument("Correlation length scale must be positive");
    }
}

// Backward compatibility wrapper
void PropertyGenerator::assign_gaussian_score()
{
    generateField();
}

// Alternative implementation with progress tracking (optional enhancement):
void PropertyGenerator::generateKSatFieldWithProgress(std::function<void(unsigned int, unsigned int)> progressCallback)
{
    validateFieldGeneration();

    if (size() == 0) {
        return;
    }

    initializeFirstPoint();

    unsigned int pointsGenerated = 1;
    const unsigned int totalPoints = static_cast<unsigned int>(size());

    if (progressCallback) {
        progressCallback(pointsGenerated, totalPoints);
    }

    while (!allPointsDetermined()) {
        const unsigned int targetIndex = selectRandomUndeterminedPoint();
        assignKSatGaussian(targetIndex);

        ++pointsGenerated;
        if (progressCallback) {
            progressCallback(pointsGenerated, totalPoints);
        }
    }
}

void PropertyGenerator::assign_K_gauss()
{
    unsigned int n_filled = 0;
    srand(time(NULL));
    int i = gsl_rng_uniform(r)*(size()-1) + 0.5;
    at(i).normal_scores.K_sat = gsl_ran_ugaussian(r);
    at(i).k_det = true;
    n_filled++;
    while (n_filled<size())
    {
        i = gsl_rng_uniform(r)*(size()-1) + 0.5;
        if (!at(i).k_det)
        {
            assign_K_gauss(i);
            n_filled++;
        }
    }

}


unsigned int PropertyGenerator::GetNumberOfPointsDetermined() const
{
    int num_determined = 0;
    for (unsigned int i=0; i<size(); i++)
        if (at(i).k_det)
            num_determined++;
    return num_determined;
}

std::vector<int> PropertyGenerator::Determined()
{
    std::vector<int> determined;
    for (unsigned int i=0; i<size(); i++)
        if (at(i).k_det)
            determined.push_back(i);
    return determined;
}


void PropertyGenerator::Normalize_Ksat_normal_scores(const double &new_mean, const double &new_std)
{
    double m = mean("K_sat_normal_score");
    double s= std("K_sat_normal_score");
    for (unsigned int i=0; i<size(); i++)
    {
        at(i).normal_scores.K_sat = (at(i).normal_scores.K_sat - m)/s*new_std + new_mean;
    }
}

double PropertyGenerator::mean(const std::string &quan, bool log) const
{
    CVector V(vals(quan));
    if (!log)
        return V.mean();
    else
        return exp(V.Log().mean());
}

double PropertyGenerator::std(const std::string &quan, bool log) const
{
    CVector V(vals(quan));
    if (!log)
        return V.stdev();
    else
        return V.Log().stdev();
}

void PropertyGenerator::Normalize(const std::string &quan,const double &denominator)
{
    for (int i=0; i<size(); i++)
    {
        SetVal(quan,i,val(quan,i)/denominator);
    }
}

bool PropertyGenerator::write(const std::string &quan, const std::string &filename) const
{
    CVector V(vals(quan));
    V.writetofile(filename);
    return true;
}


std::vector<double> PropertyGenerator::vals(const std::string &quan) const
{
    std::vector<double> out;
    for (unsigned int i=0; i<size(); i++)
    {
        if (quan == "K_sat_normal_score")
            out.push_back(at(i).normal_scores.K_sat);
        else if (quan == "alpha_normal_score")
            out.push_back(at(i).normal_scores.alpha);
        else if (quan == "n_normal_score")
            out.push_back(at(i).normal_scores.n);
        if (quan == "K_sat")
            out.push_back(at(i).realvalues.K_sat);
        else if (quan == "alpha")
            out.push_back(at(i).realvalues.alpha);
        else if (quan == "n")
            out.push_back(at(i).realvalues.n);

    }

    return out;
}

double PropertyGenerator::val(const std::string &quan, int i) const
{

    if (quan == "K_sat_normal_score")
        return at(i).normal_scores.K_sat;
    else if (quan == "alpha_normal_score")
        return at(i).normal_scores.alpha;
    else if (quan == "n_normal_score")
        return at(i).normal_scores.n;
    if (quan == "K_sat")
        return at(i).realvalues.K_sat;
    else if (quan == "alpha")
        return at(i).realvalues.alpha;
    else if (quan == "n")
        return at(i).realvalues.n;
    else
        return -999;


}

bool PropertyGenerator::SetVal(const std::string &quan, int i, const double &value)
{

    if (quan == "K_sat_normal_score")
        at(i).normal_scores.K_sat = value;
    else if (quan == "alpha_normal_score")
        at(i).normal_scores.alpha = value;
    else if (quan == "n_normal_score")
        at(i).normal_scores.n = value;
    if (quan == "K_sat")
        at(i).realvalues.K_sat = value;
    else if (quan == "alpha")
        at(i).realvalues.alpha = value;
    else if (quan == "n")
        at(i).realvalues.n = value;
    else
        return false;

    return true;


}

void PropertyGenerator::SetMarginalDistribution(const std::string &quan, const CTimeSeries<double> series)
{
    marginal_distributions[quan] = series;
}
CTimeSeries<double> PropertyGenerator::MarginalDistribution(const std::string &quan)
{
    if (marginal_distributions.count(quan)==1)
    {
        return marginal_distributions[quan];
    }
    else
        return CTimeSeries<double>();
}

void PropertyGenerator::PopulateRealValue(const std::string &quan, const std::string &quanfrom)
{
    //cout<<"PopulateRealValue-1"<<std::endl;
    CTimeSeries<double> LogTransformed = marginal_distributions[quan].LogTransformX();
    //cout<<"LogTransformed"<<std::endl;
    CTimeSeries<double> inverseCDF = LogTransformed.inverse_cumulative_uniform(100);
    //cout<<"PopulateRealValue-2"<<std::endl;
    for (unsigned int i=0; i<size(); i++)
    {
        //cout<<i<<std::endl;
        double score = gsl_cdf_ugaussian_P(val(quanfrom,i));
        double value = exp(inverseCDF.interpol(score));
        SetVal(quan,i,value);
    }
}

void PropertyGenerator::SetCorr(params p, const double &value)
{
    if (p==params::alpha)
        K_sat_alpha_correlation = value;
    else if (p==params::n)
        K_sat_n_correlation = value;
}

void PropertyGenerator::Populate_Alpha_n_normal_scores(params p)
{
    for (unsigned int i=0; i<size(); i++)
    {
        if (p == params::alpha)
            at(i).normal_scores.alpha = K_sat_alpha_correlation + gsl_ran_ugaussian(r)*sqrt(1-pow(K_sat_alpha_correlation,2));
        if (p== params::n)
            at(i).normal_scores.n = K_sat_n_correlation + gsl_ran_ugaussian(r)*sqrt(1-pow(K_sat_n_correlation,2));
    }
}
