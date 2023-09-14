#include "OdeSolver.h"

double k_i_tke,
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

void odesystem(const std::vector<double> &u, std::vector<double> &dudt, const double t){
    //pi_v = 1.9989,//1.47,
    //k_i_tke = 0.0875702,
    k_v1 = 0.00024600, // 0.00982,
    k_v2 = 0.000061,
    k_v3 = 0.0645,
    beta_Ap = 1, // 0.179,
    c_ap1 = 8, //0.968483,//8,
    c_ap2 = 8080000, //1e08, // 8080000
    delta_Apm = 0.04,
    beta_apm = 0.00047707,  
    alpha_th = 0.000217,
    alpha_tk = 0.0217,
    pi_th = 1e-08,
    delta_th = 0.3,
    beta_tk = 0.0000143,
    pi_tk = 1e-08,
    delta_tk = 0.3,
    alpha_B = 8,
    pi_B1 = 0.0000898,
    pi_B2 = 1.27e-08,
    beta_ps = 6e-06,
    beta_pl = 5e-06,
    beta_Bm = 1e-06,
    delta_S = 2.5,
    delta_L = 0.35,
    gamma_bm = 9.75e-04,
    //pi_AS = 0.087,
    //pi_AL = 0.001,
    delta_am = 0.07,
    delta_ag = 0.07,
    pi_c_apm = 0.000313313,
    pi_c_i = 0,//0.000945578,
    pi_c_tke = 2.22783e-05, 
    delta_c = 17.8634, 
    k_bm1 = 0.01,
    k_bm2 = 25000,
    beta_th = 0.000018,
    lambda_th = beta_th,
    //lambda_the = 10 *lambda_th,
    react_thmi = 0,
    react_thm = 0,
    delta_thm = 1e-03,
    beta_thn = 9.98996e-06,
    //beta_the = 6.74371e-06,
    p_igg = 0.01,
    p_igm = 0.05,
    k_thmi = 0,
    alpha_ap = 0.5, 
    delta_the = 2*delta_th;
    //beta_thm = 9.99043e-06;//beta_the/10; //beta_the;

    double V = u[0], Ap = u[1], Apm = u[2], I = u[3], Thn = u[4], The = u[5],
           Tkn = u[6], Tke = u[7], B = u[8], Ps = u[9], Pl = u[10], Bm = u[11],
           IgM = u[12], IgG = u[13], C = u[14], Thmi = u[15], Thm = u[16];
    
    dudt[0] = pi_v * (I + Thmi) - k_v1 * V * IgM - k_v2 * V * IgG; //K_v2 para IgG = 2*k_v2 para IgM
    dudt[1] = (alpha_ap) * (Ap0 - Ap) - Ap * (c_ap1 * (V) / (c_ap2 + V));
    dudt[2] = Ap * (c_ap1 * (V) / (c_ap2 + V)) - delta_Apm * Apm;
    dudt[3] = beta_thn * Thn * V + beta_the * The * V + beta_thm * Thm * V - k_i_tke * I * Tke;
    dudt[4] = alpha_th * (Thn0 - Thn) - lambda_th * Apm * Thn - beta_thn * Thn * V;
    dudt[5] = lambda_th * Apm * Thn + pi_th * Apm * The - delta_the * The - lambda_the * The - beta_the * The * V + react_thm * V * Thm;
    dudt[6] = (alpha_tk) * (Tkn0 - Tkn) - beta_tk * Apm * Tkn;
    dudt[7] = beta_tk * Apm * Tkn + pi_tk * Apm * Tke - delta_tk * Tke;
    dudt[8] = alpha_B * (B0 - B) + pi_B1 * V * B + pi_B2 * The * B - beta_ps * Apm * B -
              beta_pl * The * B - beta_Bm * The * B;
    dudt[9] = beta_ps * Apm * B - delta_S * Ps;
    dudt[10] = beta_pl * The * B - delta_L * Pl + gamma_bm * Bm;
    dudt[11] = beta_Bm * The * B + k_bm1 * Bm * (1 - Bm / (k_bm2)) - gamma_bm * Bm;
    dudt[12] = p_igm * Ps - delta_am * IgM;
    dudt[13] = p_igg * Pl - delta_ag * IgG;
    dudt[14] = pi_c_apm * V * Apm + pi_c_i * I + pi_c_tke * V * Tke - delta_c * C;
    dudt[15] = beta_thm * Thm * V - k_thmi * Thmi * Tke;
    dudt[16] = lambda_the * The - delta_thm * Thm - beta_thm * Thm * V - react_thm * V * Thm;
}

std::vector<double> advanceStep(double t, double dt, std::vector<double> u){
    runge_kutta_cash_karp54<std::vector<double>> stepper;
    stepper.do_step(odesystem, u, t, dt);
    return u;
}

void solve(ODE *ode, double tfinal, double dt, std::vector<std::string> varNames,
           std::vector<double> u0, std::string fname){
    double t = 0;
    ode->u = u0;
    createFile(fname, varNames);

    std::ofstream fp;
    fp.open(fname, ios::app);
    save(fp, 0, ode->u);

    runge_kutta_cash_karp54<std::vector<double>> stepper;
    auto c_stepper = make_controlled(1.E-08, 1.E-08, stepper);

    for (t = dt; t <= tfinal; t += dt){
        c_stepper.stepper().do_step(ode->odeSystem, ode->u, t, dt);
        if (((int)t % 10) <=  1e-06)
            save(fp, t, ode->u);
    }
    fp.close();
}

std::vector<std::vector<double>> readCSV_to_MultidimensionalArray(std::string fname, bool readFirstLine){
    std::ifstream f(fname);
    std::string line, val;                  /* string for line & value */
    std::vector<std::vector<double>> array; /* vector of vector<double>  */
    if (readFirstLine == false)
        std::getline(f, line);

    while (std::getline(f, line)){                                    /* read each line */
        std::vector<double> v;           /* row vector v */
        std::stringstream s(line);       /* stringstream line */
        while (getline(s, val, ','))     /* get each value (',' delimited) */
            v.push_back(std::stod(val)); /* add to row vector */
        array.push_back(v);              /* add row vector to array */
    }

    return array;
}

void createFile(std::string name, std::vector<std::string> varNames){
    std::ofstream fp;
    fp.open(name);
    fp << "t,";
    for (std::string n : varNames)
        fp << n << ",";
    fp << std::endl;
    fp.close();
}

void save(std::ofstream &fp, double t, std::vector<double> values){
    fp << t;
    int cont = 0;
    for (double v : values){
        /*if (cont == 12 || cont == 13){
            if (v <= 1)
                fp << "," << 0;
            else
                fp << "," << log2(v);
        }
        else*/

        fp << "," << (v);
        cont++;
    }
    fp << endl;
}
