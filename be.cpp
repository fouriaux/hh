/*
 *                Solving HH model with Backward Euler method
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
 */

#include <iostream>
#include <math.h>

static const double PERIOD = 0.1;

static const double E_k  = -72.14;   // Resting potential of Potatium channel (mV)
static const double E_Na = 55.17;    // Resting potential of Sodium channel   (mV)
static const double E_l  = -49.42;   // Resting potential of passive membrane (mV)
static       double g_K  = 0.36;     //  ion channels conductances (mS/cm²)
static       double g_Na = 1.2;      //  ion channels conductances (mS/cm²)
static       double g_l  = 0.003;    //  ion channels conductances (mS/cm²)

static double Cm      = 0.01;        //  Membrane capacitance (µF/cm²)
static double I       = 0.0;

static double v_prev  = -60.0;       //  Potential at t = n
static double v       = -60.0;       //  Potential at t = n+1
static double n       = 0.0;         // state variable at t = n+½
static double m       = 0.0;         // state variable at t = n+½
static double h       = 0.0;         // state variable at t = n+½
 
static double n_prev  = 0.0;         // state variable at t = n-½
static double m_prev  = 0.0;         // state variable at t = n-½
static double h_prev  = 0.0;         // state variable at t = n-½

static const double delta_t = 0.05; // timestep size (ms) 
static const double delta_t_inv = 1.0/ delta_t;
static const double end_of_times = 25.0;
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
 * Some backward Euler theory
 *
 * solving dv/dt = f(t,v) with v(0) = v₀ (-60mV in this particular case)
 * Backward Euler Method is
 * 
 *    v₊₁ = v + Δt*f(t₊₁,v₊₁)
 * →  v₊₁ - v - Δt*f(t₊₁, v₊₁) = 0
 *
 * In our case we will consider 2 interleaved intervals 
 * to solve "independently" channel gates and potential equations.
 * 
 * v will be solve with initial t'₀ = t₀+½
 * then when computing gates equations we consider v to be constant on t = [n, n₊₁]
 * conversely when computing v we will consider channels   constant on t = [n₋½;n₊½]
 * 
 * 
 */




void backward_euler (double I) {

    // solve potential
    double G_Na = g_Na*m*m*m*h;
    double G_K  = g_K*n*n*n*n;
    double new_v = 2*delta_t * ( I + G_Na*E_Na + G_k*E_k + g_l*E_l);
    //     new_v = ------------------------------------------------
    new_v       /= (2*c_m + (G_Na + G_k + G_l)*delta_t);
    new_v += v;
    v_prev = v;
    v = new_v;
    
    // compute channel states
     double gamma_n       = (get_alpha_n(v) + get_beta_n(v))/2.0;
     double new_n =  get_alpha_n(v)*1.0/(delta_t_inv + gamma_n); 
     new_n -= n_prev *(gamma_n - delta_t_inv)/(gamma_n + delta_t_inv);
     n_prev = n;
     n = new_n;
     
     double  gamma_m       = (get_alpha_m(v) + get_beta_m(v))/2.0;
     double new_m =  get_alpha_m(v)*1.0/(delta_t_inv + gamma_m); 
     new_m -= n_prev *(gamma - delta_t_inv)/(gamma + delta_t_inv);
     m_prev = m;
     m = new_m; 
     
     double gamma_h       = (get_alpha_h(v) + get_beta_h(v))/2.0;
     double new_h =  get_alpha_h(v)*1.0/(delta_t_inv + gamma_h); 
     double new_h -= n_prev *(gamma_h - delta_t_inv)/(gamma_h + delta_t_inv);
     h_prev = h;
     h = new_h;
}

void init () {
    v = -60.0;
    v_prev = -60.0;
    I = 0;
    for (double init_t = 0; init_t < 20.0; init_t+=delta_t) {
        backward_euler(I);
    }
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
   std::cout << t << "," << "ForwardEuler," 
             << v_euler << ","
             << n << ","
             << m << ","
             << h << ","
             << m*m*m*h*g_Na + n*n*n*n*g_K << ","
             << I << std::endl;
}

int main (int argc, char** argv) {
    double t = 0.0;
    printHeader ();
    init(v);
    for (t = 0.0; t < end_of_times; t+=delta_t) {
        I = constant_current(t);
        backward_euler(I);
        print (t);
    }
    return 0;
}
