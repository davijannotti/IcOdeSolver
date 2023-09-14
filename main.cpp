#include <cmath>
#include <ctime>
#include <iterator>
#include <map>
#include <cstdlib>
#include <fstream>
#include <numeric>
#include <iostream>
#include <bits/stdc++.h>
#include <boost/numeric/odeint.hpp>
#include "./de/DifferentialEvolution.h"

#include "OdeSolver.h"

using namespace std;
using namespace de;
using namespace boost::numeric::odeint;

class TestDE : public de::IOptimizable{
private:
    double tfinal;
    double dt;
    std::vector<std::vector<double>> vData;
    std::vector<std::vector<double>> il6Data;
    std::vector<std::vector<double>> antibodyData;

public:

    TestDE(double tf, double deltat) : tfinal(tf), dt(deltat) {}
    void init(){
        vData = readCSV_to_MultidimensionalArray("data/viremia_data.csv", false);
    }
    double EvaluateCost(std::vector<double> inputs) const override{
        int s = 18;
        std::vector<double> u;
        u.reserve(s);
        u.resize(s);
        V0 = 80;
        Ap0 = 1e06;
        Thn0 = 1e06;
        Tkn0 = 5e05;
        B0 = 2.5e05;

        //{"V", "Ap", "Apm", "I", "Thn", "The", "Tkn", "Tke", "B", "Ps", "Pl", "Bm", "IgM", "IgG", "C"}
        u[0] = V0;
        u[1] = Ap0;
        u[2] = 0;
        u[3] = 0;
        u[4] = Thn0;
        u[5] = 0;
        u[6] = Tkn0;
        u[7] = 0;
        u[8] = B0;
        u[9] = 0;
        u[10] = 0;
        u[11] = 0;
        u[12] = 0;
        u[13] = 0;
        u[14] = 0;
        u[15] = 0;
        u[16] = 0;
        u[17] = 0;
    
        pi_v = inputs[0];
        beta_thm = inputs[1];
        beta_the = inputs[2];
        lambda_the = inputs[3];
        k_i_tke = inputs[4];
        //k_thmi = inputs[4];
        

        double vError = 0, vSum = 0;
        int indexV = 0;

        for (double t = 0; t <= tfinal; t += dt){

            if (indexV < vData.size() && abs(t - vData[indexV][0]) < dt){
                double V = vData[indexV][1];
                double vOde = u[0];
                vError += (vOde - V) * (vOde - V);
                vSum += V * V;

                indexV++;
            }

            u = advanceStep(t, dt, u);
            if (u[0] < 0)
                break;
        }

        /*if (vError > 0)
            vError = sqrt(vError/vSum);*/

        return vError;
    }

    unsigned int NumberOfParameters() const override{
        return 5;
    }

    std::vector<Constraints> GetConstraints() const override{
        std::vector<Constraints> constr(NumberOfParameters());
        constr[0] = Constraints(0.1, 2, true);
        constr[1] = Constraints(1e-06, 1e-02, true);
        constr[2] = Constraints(1e-06, 1e-02, true);
        constr[3] = Constraints(1e-06, 1e-02, true);
        constr[4] = Constraints(1e-05, 1e-02, true);
        //constr[5] = Constraints(1e-05, 1e-02, true);

        // Gravar em um arquivo as restricoes usadas
        return constr;
    }

    static bool terminationCondition(const DifferentialEvolution &de){
        if (de.GetBestCost() <= 0.1)
            return true;
        return false;
    }
};

std::vector<double> ajuste(int tfinal, double dt){
    TestDE deInstance(tfinal, dt);
    deInstance.init();
    int populationSize = 100, maxIterations = 30;
    de::DifferentialEvolution de(deInstance, populationSize,
                                 std::time(nullptr), true, TestDE::terminationCondition);
    std::pair<std::vector<double>, std::vector<double>> costs = de.Optimize(maxIterations, true);

    return de.GetBestAgent();
}



int main(){
    int tfinal = 800;
    double dt = 0.0001;
    // 1.20963 0.03049 0.12209 0.00383 0.96523 52.00420 4.01182 177.29352
    // 1.07266 0.04945 0.03945 0.00481 0.97041 164.07535 2.73649 137.60471
    std::vector<double> best = ajuste(tfinal, dt);

    pi_v = best[0];
    beta_thm = best[1];
    beta_the = best[2];
    lambda_the = best[3];
    k_i_tke = best[4];
    //k_thmi = best[4];
    

    ofstream fp("best_params.txt");
    fp << best[0] << std::endl;
    fp << best[1] << std::endl;
    fp << best[2] << std::endl;
    fp << best[3] << std::endl;
    fp << best[4] << std::endl;
    // fp << best[5] << std::endl;
    // fp << best[5] << std::endl;   
    // fp << best[6] << std::endl;
    // fp << best[7] << std::endl;
    // fp << best[8] << std::endl;
    // fp << best[9] << std::endl;
    fp.close();

    ODE ode(odesystem);
    std::vector<std::string> populations = {"V", "Ap", "Apm", "I", "Thn", "The", "Tkn",
                                            "Tke", "B", "Ps", "Pl", "Bm", "IgM", "IgG", "C", "Thmi", "Thm"};
    std::string fname = "ode_output.csv";

    V0 = 80;
    Ap0 = 1e06;
    Thn0 = 1e06;
    Tkn0 = 5e05;
    B0 = 2.5e05;

    solve(&ode, tfinal, dt, populations, {V0, Ap0, 0, 0, Thn0, 0, Tkn0, 0, B0, 0, 0, 0, 0, 0, 0, 0, 0}, fname);

    return 0;
}