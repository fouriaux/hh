/*
 *                                        CABLE EQUATION
 *                           Solving HH model using Forward Euler Method
 *
 *     dVi      1    /                                                                              \
 *     ---  = ---- * | Ii - g_Na*m³*h(Vi-E_Na) - G_K*n⁴(Vi - E_K) - gl(Vi - El) - Gij * (Vi - Vj)   |
 *     dt      Cm    \                                                                              /
 *
 *     => this is the cable equation for everywhere exept at branching points
 *
 *             
 *          \j \  / k/
 *           \  \/  /
 *            \/__\/.....Center Have 0 Conductance because 0 surface
 *             |  |
 *             |  |
 *             |l |
 *
 *
 *     => We first compute U₁ on every surounding compartments with U₋₁ on center
 *     => Then center is computed algebraicaly.
 * 
 *        Ii + GijVj + GikVk + GilVl
 *  Vi = ----------------------------
 *               Gij + Gik + Gil
 *
 *
 *     dn
 *     ---  =  alpha_n(v)*(1-n) - beta_n(v)*(n)
 *     dt
 *
 *     idem for m and h
 */

#include <iostream>
#include <math.h>
#include <cstdlib>
static const double PERIOD = 0.1;

static const double E_k  = -72.14;   // (mV)
static const double E_Na = 55.17;    // (mV)
static const double E_l  = -49.42;   // (mV)
static       double g_K  = 0.36;     //  ion channels conductances (mS/cm²)
static       double g_Na = 1.2;      //  ion channels conductances (mS/cm²)
static       double g_l  = 0.003;    //  ion channels conductances (mS/cm²)

static double* v;                    //  Potential at t computed by Forward Euler method
static double* Cm;                   //  Membrane capacitance (µF/cm²)
static double* I;
static double* n;
static double* m;
static double* h;

static double delta_t = 0.025; // timestep size (ms) 
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

void forward_euler () {
    for (int i = 0 ; i < nb_comp; i++) {
        // solve potential
        if (is_real_node[i]) {
            double dv = (1/Cm) * ( I[i] -
                               g_Na*pow(m[i],3)*h[i]*(v[i]- E_Na) -
                               g_K*pow(n[i],4)*(v[i]-E_k) -
                               g_l*(v[i]-E_l) -
                               G[i-1]*(v[i] - v[i-1]) -
                               G[i+1]*(v[i] - v[i+1]));
            v[i] += dv*delta_t;
    
            // compute channel states
            double dn = get_alpha_n(v[i])*(1-n[i]) - get_beta_n(v[i])*n[i];
            n[i] += dn * delta_t;
        
            double dm = get_alpha_m(v[i])*(1-m[i]) - get_beta_m(v[i])*m[i];
            m[i] += dm * delta_t;
    
            double dh = get_alpha_h(v[i])*(1-h[i]) - get_beta_h(v[i])*h[i];
            h[i] += dh * delta_t;
        }
    }
    for (int i = 0; i < nb_virtual; i++) {
        int  idx    = virtual_nodes[i];
        int* childs = get_child(idx);
        v[idx] = (I[idx] + Gi[idx-1] * v[idx-1] + Gi[childs[0]] * Vi[childs[0]] +Gi[childs[1]] * Vi[childs[1]]) /
                        (Gi[idx] + Gi[childs[0]] + Gi[childs[1]]);
    }
}

void init (double& v) {
    v = -60.0;
    I = 0;
    for (double init_t = 0; init_t < 20.0; init_t+=delta_t) {
        forward_euler(I, v_euler);
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
    if(argc == 2) {
        delta_t = atof(argv[1]);
    }
    printHeader ();
    init(v_euler);
    for (t = 0.0; t < end_of_times; t+=delta_t) {
        I = constant_current(t);
        forward_euler(I, v_euler);
        print (t);
    }
    return 0;
}
