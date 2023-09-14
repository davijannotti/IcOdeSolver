#ifndef H_ODESOLVER
#define H_ODESOLVER

#include <cmath>
#include <ctime>
#include <ostream>
#include <sstream>
#include <unordered_map> 
#include <fstream> 
#include <numeric>
#include <boost/numeric/odeint.hpp>
#include <boost/algorithm/string.hpp>
//#include "matplotlibcpp.h"

using namespace std;
using namespace boost::algorithm;
using namespace boost::numeric::odeint;
//namespace plt = matplotlibcpp;

typedef struct TODE {    
    std::vector<double> u; //Variables's vector     
    std::function<void (const std::vector<double> &u , std::vector<double> &dudt , const double /* t */)> odeSystem;
    TODE(){}
    TODE(std::function<void (const std::vector<double> &u , std::vector<double> &dudt , const double /* t */)> s){
        odeSystem = s;
    }
} ODE; 

void odesystem(const std::vector<double> &u , std::vector<double> &dudt , const double /* t */) ;

std::vector<double> advanceStep(double t, double dt, std::vector<double> u);

void solve(ODE *ode, double tfinal, double dt, std::vector<std::string> varNames, 
                std::vector<double> u0, std::string fname); 

std::vector< std::vector<double> > readCSV_to_MultidimensionalArray(std::string fname, bool readFirstLine);

void createFile(std::string name, std::vector<std::string> varNames); 

void save(std::ofstream &fp, double t, std::vector<double> values); 

extern double k_i_tke,
    pi_v,
    k_v1,
    k_v2,
    k_v3,
    beta_Ap,
    c_ap1,
    c_ap2,
    delta_Apm,
    beta_apm,
    alpha_th,
    alpha_tk,
    beta_the,
    beta_th,
    pi_th,
    delta_th,
    beta_tk,
    pi_tk,
    delta_tk,
    alpha_B,
    pi_B1,
    pi_B2,
    beta_ps,
    beta_pl,
    beta_Bm,
    delta_S,
    delta_L,
    gamma_bm,
    pi_AS,
    pi_AL,
    delta_am,
    delta_ag,
    pi_c_apm,
    pi_c_i,
    pi_c_tke,
    delta_c,
    V0,
    Ap0,
    Tkn0,
    Thn0,
    B0,
    k_bm1,
    k_bm2,
    lambda_th,
    lambda_the,
    react_thmi,
    react_thm,
    beta_thn,
    delta_thm,
    p_igg,
    p_igm,
    k_thmi,
    alpha_ap,
    delta_the,
    beta_thm;

#endif 