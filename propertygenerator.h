#ifndef PROPERTYGENERATOR_H
#define PROPERTYGENERATOR_H

#include <vector>
#include <Matrix.h>
#include <Vector.h>
#include <gsl/gsl_rng.h>
#include "BTC.h"

//using namespace std;

struct correl_mat_vec
{
#ifdef  _arma
    CMatrix_arma M_22;
    CVector_arma V_21;
    CVector_arma V_RHS;
#else
    CMatrix M_22;
    CVector V_21;
    CVector V_RHS;
#endif //  arma
};

struct ival
{
    int i;
    double val;
};

struct Prop
{
    double K_sat;
    double alpha;
    double n;
};

struct Propfull
{
    Prop realvalues;
    Prop normal_scores;
    bool k_det = false;
};

using KnownPoint = std::pair<unsigned int, double>;  // (index, value)
using KnownPoints = std::vector<KnownPoint>;

enum class params {K_sat, alpha, n};

class PropertyGenerator: public std::vector<Propfull>
{

public:
    enum class InitializationStrategy {
        RANDOM_LOCATION,
        CENTER_POINT,
        MULTIPLE_SEEDS
    };

    PropertyGenerator(int seed=0);
    PropertyGenerator(unsigned int n, int seed);
    ~PropertyGenerator();
    double dx;
    double correlation_length_scale;
    void assign_gaussian_score();
    bool write(const std::string &quan, const std::string &filename) const;
    void Normalize_Normal_Scores(const double &mean, const double &std);
    double mean(const std::string &quan, bool log=false) const;
    double std(const std::string &quan, bool log=false) const;
    std::vector<double> vals(const std::string &quan) const;
    void generateField();
    void SetMarginalDistribution(const std::string &quan, const CTimeSeries<double> series);
    CTimeSeries<double> MarginalDistribution(const std::string &quan);
    void PopulateRealValue(const std::string &quan, const std::string &quanfrom);
    double val(const std::string &quan, int i) const;
    bool SetVal(const std::string &quan, int i, const double &value);
    void SetCorr(params, const double &value);
    void Populate_Correlated_Normal_Scores(params p);
    void Normalize(const std::string &quan, const double &denominator);
    void initializeFirstPoint(const KnownPoints& knownPoints);
    void generateKSatField(const KnownPoints& knownPoints);
    void setKnownPoints(const KnownPoints& knownPoints);
    std::vector<unsigned int> getKnownPointIndices() const;
    KnownPoints getKnownPoints() const;
    bool hasInitializedPoints() const;
    void resetAllPoints();
private:
    void initializeRNG(int seed);
    unsigned long generateSeed(int seed) const;
    CMatrix K_alpha_n_corr_matrix;
    double K_sat_normal_score_mean;
    double K_sat_normal_score_std;
    void assign_Gaussian_scores(unsigned int targetIndex);
    struct ConditionalGaussianParams {
        double mean;
        double variance;
    };
    ConditionalGaussianParams computeConditionalGaussianParams(const correl_mat_vec& correlationData) const;
    void validateTargetIndex(unsigned int targetIndex) const;

    correl_mat_vec buildCorrelationMatrix(int targetIndex) const;
    void initializeFirstPoint(InitializationStrategy strategy = InitializationStrategy::RANDOM_LOCATION);
    void initializeAtRandomLocation();
    void initializeAtCenter();
    void initializeMultipleSeeds(unsigned int numSeeds = 3);
    unsigned int generateRandomIndex() const;
    void setPointValue(unsigned int index, double value);
    void validatePointForInitialization(unsigned int index) const;
    double computeCorrelation(int i, int j) const;
    void validateCorrelationInputs(int targetIndex) const;
    void initializeFirstPoint();
    unsigned int selectRandomUndeterminedPoint() const;
    bool allPointsDetermined() const;
    void validateFieldGeneration() const;
    void generateKSatFieldWithProgress(std::function<void(unsigned int, unsigned int)> progressCallback = nullptr);
    correl_mat_vec get_correll_matrix_vec(int i);
    std::vector<ival> get_top_n(const std::vector<ival> &vec);
    std::vector<ival> get_closest_val_dets(unsigned int i);
    unsigned int GetNumberOfPointsDetermined() const;
    std::vector<int> Determined();
    std::map<std::string,CTimeSeries<double>> marginal_distributions;
    double K_sat_alpha_correlation = 1;
    double K_sat_n_correlation = 1;
    const gsl_rng_type * T;
    gsl_rng * r = gsl_rng_alloc (gsl_rng_taus);

};

#endif // PROPERTYGENERATOR_H
