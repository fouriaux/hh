/*
 *                Solving HH model with Crank-Nicholson method
 *
 *     dv      1   /                                                        \
 *     ---  = ----*| I - g_Na*m³*h(v-E_Na) - G_K*n⁴(v - E_K) - gl(v - El)   |
 *     dt      Cm  \                                                        /
 *
 *     dn
 *     ---  =  alpha_n(v)*(1-n) - beta_n(v)*(n)
 *     dt
 *
 *     idem for m and h
 *
 *
 *     Crank Nicholson solution:
 *     We first interleave time steps for v and all m,n,h computations.
 *
 *
 *
 */

#include <iostream>
#include <math.h>
#include <cstdlib>

static const double PERIOD = 0.1;

static const double E_K  = -72.14;   // Resting potential of Potatium channel (mV)
static const double E_Na = 55.17;    // Resting potential of Sodium channel   (mV)
static const double E_l  = -49.42;   // Resting potential of passive membrane (mV)
static const double g_K  = 0.36;     //  ion channels conductances (mS/cm²)
static const double g_Na = 1.2;      //  ion channels conductances (mS/cm²)
static const double g_l  = 0.003;    //  ion channels conductances (mS/cm²)

static double Cm         = 0.01;     //  Membrane capacitance (µF/cm²)
static double I          = 0.0;

static double v_prev     = -60.0;    //  Potential at t = n
static double v          = -60.0;    //  Potential at t = n+1

static double n          = 0.0;      // state variable at t = n+½
static double m          = 0.0;      // state variable at t = n+½
static double h          = 0.0;      // state variable at t = n+½
 
static double n_prev     = 0.0;      // state variable at t = n-½
static double m_prev     = 0.0;      // state variable at t = n-½
static double h_prev     = 0.0;      // state variable at t = n-½

static double delta_t = 0.005; // timestep size (ms) 
static double delta_t_inv = 1.0/ delta_t;
static const double end_of_times = 200.0;

double get_alpha_n (double v) {
    return (0.01*(v + 50.0))/(1 - exp(-(v + 50.0)/10.0));
}

double get_beta_n  (double v) {
    return 0.125*exp(-(v + 60.0)/80.0);
}

double get_alpha_m (double v) {
    return (0.1*(v+35.0))/(1- exp(-(v+35.0)/10.0));
}

double get_beta_m (double v) {
    return 4.0*exp(-0.0556*(v+60.0));
}

double get_alpha_h (double v) {
    return 0.07*exp(-0.05*(v+60.0));
}

double get_beta_h (double v) {
    return 1/(1+exp(-0.1*(v+30.0)));
}
/*
 * Some Crank Nicholson theory
 *
 * solving dv/dt = f(t,v) with v(0) = v₀ (-60mV in this particular case)
 * Crank-Nicholson Method is
 * 
 *                   f(v₊₁) + f(v)
 *    v₊₁ = v + Δt * -------------
 *                         2
 *
 * In our case we will consider 2 interleaved intervals 
 * to solve "independently" channel gates and potential equations.
 * 
 * v will be solve with initial t'₀ = t₀+½
 * When computing gates equations on interval t = [n, n₊₁], 
 *    v is "constant" with value determined in middle of [n, n₊₁] at t₊½.
 * Conversely when computing v we will consider channels constant 
 *    on interval t = [n₋½;n₊½]
 *
 *
 */


void crank_nicholson (double I) {

    // solve potential
    double G_Na = g_Na*m*m*m*h;
    double G_K  = g_K*n*n*n*n;
//    G_K = 0; G_Na = 0;
    double new_v     =  (2*delta_t*(g_l*E_l+G_Na*E_Na+G_K*E_K+I)) - v*(delta_t*(g_l+G_Na+G_K) - 2*Cm);
    new_v           /=  (2*Cm + delta_t*(g_l+G_Na+G_K));
    v = new_v;
    
    // compute channel states
     double gamma_n       = (get_alpha_n(v) + get_beta_n(v))/2.0;
     double new_n =  get_alpha_n(v)/(delta_t_inv + gamma_n); 
     new_n -= n *(gamma_n - delta_t_inv)/(gamma_n + delta_t_inv);
     n = new_n;
     
     double  gamma_m       = (get_alpha_m(v) + get_beta_m(v))/2.0;
     double new_m =  get_alpha_m(v)/(delta_t_inv + gamma_m); 
     new_m       -= m *(gamma_m - delta_t_inv)/(gamma_m + delta_t_inv);
     m = new_m;
     
     double gamma_h       = (get_alpha_h(v) + get_beta_h(v))/2.0;
     double new_h =  get_alpha_h(v)/(delta_t_inv + gamma_h); 
     new_h       -= h *(gamma_h - delta_t_inv)/(gamma_h + delta_t_inv);
     h = new_h;
}


double zero_current (double t) {
    return 0.0;
}

double constant_current(double t) {
    return 0.1;
}

double pulse (double t) {
    double value    = sin(t/PERIOD * 2.0 * M_PI);
    return  (value >= 0.0 && value <= 0.5) ? 1.0:0.0;
}

int read_data (const char* input_file) {
 return 0;
}

void printHeader () {
    std::cout << "time(ms),method,potential(mV),n,m,h,g,I(mAmp)" << std::endl;
}

void print (double t) {
   std::cout << t << "," << "CrankNicholson," 
             << v << ","
             << n << ","
             << m << ","
             << h << ","
             << m*m*m*h*g_Na + n*n*n*n*g_K << ","
             << I << std::endl;
}

void init () {
    v = -60.0;
    v_prev = -60.0;
    I = 0;
    for (double init_t = 0; init_t < 20.0; init_t+=delta_t) {
//        print(init_t);
        crank_nicholson(I);
    }
}

int main (int argc, char** argv) {
    double t = 0.0;
    if(argc == 2) {
        delta_t = atof(argv[1]);
        delta_t_inv = 1.0/ delta_t;
    }
    printHeader ();
    init();
    for (t = 0.0; t < end_of_times; t+=delta_t) {
        I = constant_current(t);
        //I = 0;
        crank_nicholson(I);
        print (t);
    }
    return 0;
}
