//=== README ===================================//{{{
//·
//· Description:
//·    An implementation of
//·    M J Shelley, L Tao's
//·    Efficient and accurate time-stepping schemes for iaf neuronal networks
//·    published in Journal of Computational Neuroscience, 11, 111, 2001
//·
//· gcc -std=gnu11   %.c  -o nif   -g -lm -lGL -lGLU -lglut -Wall
//· run: nif [num_nif] [RK_order] [ifAdj] [num_t_click] [m]
//·
//· This is my 3rd try on object+functional programing.
//· mogeiwang@gmail.com
//· 5.24-6.12, 2015
//· Shanghai, China
//·
//==============================================//}}}

#include <assert.h>//{{{ //checked
#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <locale.h>
#include <math.h> // tgmath.h
#include <setjmp.h>
#include <signal.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h> // unix.h :-)
#include <GL/freeglut.h> // openGL
#include <sys/stat.h> // unix files
#include <sys/types.h>

#define BEGIN {
#define END }
#define integer long int
#define boolean integer
#define real double
#define elif else if
#define true ((integer)(1==1))   // gcc: 1=true, 0=false
#define false ((integer)(!true))
#define EQ(f1, f2)          ( fabs((real)(f1)-(real)(f2))<(1e-13))

#define max(a,b)  ( ((a)>=(b)) ? (a):(b) )
#define min(a,b)  ( ((a)<=(b)) ? (a):(b) )
#define IdentifierLiteral(T) #T // return "T"
#define IdentifierJoint(T,x) T##x // return Tx
#define  print        printf
#define fprint       fprintf
#define pprint       fprintf
#define sprint       sprintf
#define  println()    print("\n")
#define fprintln(f)  fprint(f, "\n")
#define pprintln(p)  pprint(p, "\n")
#define sprintln(s)  sprint(s, "\n")
#define  printab()    print("\t")
#define fprintab(f)  fprint(f, "\t")
#define pprintab(p)  pprint(p, "\t")
#define sprintab(s)  sprint(s, "\t")
#define print_integer(pi)   print("\t%ld", pi)
#define print_real(pr)      print("\t%.15lf", pr)
#define print_boolean(pb)   if (EQ(pb,1)) print("\ttrue"); else print("\tfalse")
#define DebugWrong(s)           { print("\n[LINE %d] [FILE %s] [WRONG]: %s\n",__LINE__,__FILE__, s);}
#define DebugCheck(c, s)        { if (c) print("\n[LINE %d] [FILE %s] [CHECK]: %s\n",__LINE__,__FILE__, s);}
#define DebugW_T_F(s)           { print("\n[LINE %d] [FILE %s] [W-T-F]: %s\n",__LINE__,__FILE__, s);}
#define DebugPrint(s)           { print("\n[LINE %d] [FILE %s] [Print]: %s\n",__LINE__,__FILE__, s);}
#define DebugError(s, t)        { print("\n[LINE %d] [FILE %s] [ERROR]: %s\nreturning...\n",__LINE__,__FILE__, s); return t;}
#define Debug0_Div(divisor)     { if (EQ(divisor,0)) print("\n[LINE %d] [FILE %s] [0_Div]: ATTEMPT TO USE 0 AS DENOMINATOR!\n",__LINE__,__FILE__);}
#define RandomSeed(s)            { if (s==0) srand(time(NULL)); else srand((int)s);}
#define string_to_integer(s,t)  t=strtol(s, NULL, 10)  // use this macro in this way:   string_to_integer(string, integer);
#define string_to_real(s,t)     t=strtod(s, NULL, 10)
#define integer_to_string(s, t) sprint(t, "%ld", s) // integer_to_string(integer, string); Reserve the memory for string first!
#define real_to_string(s, t)    sprint(t, "%lf", s)
#define BALL(x)                    typedef struct {x} // functions are threads, structs are balls ;-)
#define StartFromXOS   int main(int argc, char *argv[]) { if (!initialization(argc, argv)) DebugError("program cannot be initializated!", false);
#define EndWithHonor   return deconstruction(); }
#define array_size(a)  (sizeof((a)) / sizeof((a[0])))
#define iSwap(x,y)     { integer tmp; tmp=x; x=y; y=tmp; }
#define rSwap(x,y)     { real    tmp; tmp=x; x=y; y=tmp; }
//}}}

// --- includes and macros {{{
#define t_end                (1L) //(4L)
#define num_t_click_default  (1250000L)
#define RK_order_default     (4L)
#define ifAdj_default        (true)
#define m_default            (5L)
#define num_nif_default      (200L)


integer num_nif    = num_nif_default; // how many time clicks / time steps
integer num_t_click= num_t_click_default; // how many time clicks / time steps
integer RK_order   = RK_order_default; // RK2 or RK4
integer m          = m_default;        // m affect the synaptic induced conductance change, see eq(3) in the paper.
boolean ifAdj      = ifAdj_default;    // use the linear interpolation or not.
real    t_step     = 1.0 * t_end / num_t_click_default; // time step = t_end / time clicks, (t_start = 0); assigned in initialization()
char    myDatahubFilename[256]; // for saving filename. temp var.

const real g_leak      = 50.0; // --- for single neurons
const real V_leak      = 0.0; // the V parameters
const real V_threshold = 1.0;
const real V_E         = 14/3.0;
const real V_I         = -2/3.0;

const real Tau_e       = 0.6e-3; // used in G_sicc
const real Tau_i       = 1e-3; // used in G_sicc
const real Ge0_intensity = 25.0; // m到底用了没用？？
const real Gi0_intensity = 0.0;
const real G_angleField  = 0; //2 * M_PI;

const real coupleSpatial_s    = 1.0;
const real coupleSpatial_mu   = 0.0;
const real coupleSpatial_delta= 1.0;
//const char *myDatahubFilenameHead =  "/home/mw/neuronSci/data/vol"; // "neuron_voltage";
const char *myDatahubFilenameHead =  "vol"; // "neuron_voltage";
const char *myDatahubFilenameTail = ".txt";
//}}}

// --- nif system {{{
BALL (
    integer     ID;

    real        v;
    real        t;
    real        v_prev;
    real        t_prev;

    real        ls_t_spike[ t_end*300 ]; // t_end is in second. never seen more than 50 spikes in a second.
    integer     lenLs_t_spike;

    real        Ge0_anglePosition; // to be set to  (ID/num_nif)*G_angleField
    real        Gi0_anglePosition;
    real        Ge1_next; // 2nd factor of g_e. see Eq.(5)
    real        Gi1_next;
    real        ge0_next; // 1st term of g_e(t), i.e., outter drive. see Eq.(5)
    real        gi0_next;
    real        ge1_next; // 2nd term of g_e(t), i.e., neighbors' drive. see Eq.(5)
    real        gi1_next;

    real        alpha_curr; // alpha 0
    real        alpha_halfNext; // alpha 0.5
    real        alpha_next; // alpha 1
    real        beta_curr; // beta 0
    real        beta_halfNext; // beta 0.5
    real        beta_next; // beta 1
    FILE        *datahub;

    real        Ge1_curr; // really necessary???
    real        Gi1_curr;
    real        Ge1_halfNext;
    real        Gi1_halfNext; // really necessary???

    real        ge0_curr; // the previous value of g_e0, for computing alpha and beta
    real        gi0_curr;
    real        ge0_halfNext; // the previous value of g_e0, for computing alpha and beta
    real        gi0_halfNext;

    real        ge1_curr; // they are temporary variable, but it should be kept..
    real        gi1_curr;
    real        ge1_halfNext; // they are temporary variable, but it should be kept..
    real        gi1_halfNext;
) NIF;
NIF nif[num_nif_default] = {{ 0, 0.0l, 0.0l, 0.0l, 0.0l, {0.0l}, 0, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, NULL, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l, 0.0l }};

real mat_e[num_nif_default][num_nif_default] = {{0.0l}}; //coupling matrix
real mat_i[num_nif_default][num_nif_default] = {{0.0l}};

// parameter pool
BALL (
    integer     gl_RK_order;
    integer     nif_num;
    real        gl_g_leak;
    real        gl_V_leak;
    real        gl_V_E;
    real        gl_V_I;
    real        gl_V_threshold;
    real        gl_t_step;

    real        ge0_next;
    real        ge1_next;
    real        gi0_next;
    real        gi1_next;

    real        Ge1_curr; // really necessary???
    real        Gi1_curr;
    real        Ge1_halfNext;
    real        Gi1_halfNext; // really necessary???

    real        ge0_curr; // the previous value of g_e0, for computing alpha and beta
    real        gi0_curr;
    real        ge0_halfNext; // the previous value of g_e0, for computing alpha and beta
    real        gi0_halfNext;

    real        ge1_curr; // they are temporary variable, but it should be kept..
    real        gi1_curr;
    real        ge1_halfNext; // they are temporary variable, but it should be kept..
    real        gi1_halfNext;

    real        alpha_curr;
    real        alpha_halfNext;
    real        alpha_next;
    real        beta_curr;
    real        beta_halfNext;
    real        beta_next;

    real        t_prev;
    real        t;
    real        v_prev;
    real        v;

    real        last_t_spike;
) paraPool;
paraPool pp;
//}}}

// function declarations {{{
integer connection_distance(integer i, integer j, integer num_nif_ring);
real connection_intensity(integer i, integer j, integer num_nif_ring);
real RK2_stepper(real(*df)(real, real, paraPool), real t0, real x0, real ts, paraPool pp);
real RK4_stepper(real(*df)(real, real, paraPool), real t0, real x0, real ts, paraPool pp);
real dynamics(real t, real v, paraPool pp);
real adj_t_spikeDelta(paraPool pp);
real NewtonMethod(real(*eq)(real, paraPool), real(*eqDerivative)(real, paraPool), real t0, paraPool pp);
real TimeCubicInterpolation(real t_spikeDelta, paraPool pp);
real TimeCubicInterpolationDerivative(real t_spikeDelta, paraPool pp);
real adj_v(paraPool pp);
real DichotomySearch(real(*eq)(real, paraPool), real x0, real x1, integer deepth, paraPool pp);
real VoltageCubicInterpolation(real v_n, paraPool pp);
real G_sicc(real t, real tau, real m) { return pow(1.0*t/tau, m)*exp(-1.0*t/tau); };// G_spike_induced_conductance_change
//}}}

// init and decons --- part of main() //{{{
boolean init_para(int num_userCommand, char * userCommands[])//{{{
{
    if (num_userCommand == 1)
    {// use default
        print("using program default parameters, i.e.,\n");
        return true;
        // userCommands[0] is the filename of the running command.
    }
    if (num_userCommand >= 2)
    {// use userCommands: num_nif, RK_order, ifAdj, num_t_click, m
        string_to_integer(userCommands[1], num_nif);
        if (num_nif>num_nif_default || num_nif<1)  num_nif = (num_nif%num_nif_default)+1; //(1-100 neurons)
    }
    if (num_userCommand >= 3)
    {
        string_to_integer(userCommands[2], RK_order);
        RK_order = (((RK_order-2)/2)%2)*2+2; // RK2 or RK4. (map any integer to 2,4, and keep 2->2, 4->4 unchanged)
    }
    if (num_userCommand >= 4)
    {
        string_to_integer(userCommands[3], ifAdj);
        ifAdj = ifAdj%2; //use the linear interpolation or not. In gcc, 1 is true, and 0 is false.
    }
    if (num_userCommand >= 5)
    {
        string_to_integer(userCommands[4], num_t_click); // how many time clicks (number of time steps)
         t_step = 1.0*t_end/num_t_click;
    }
    if (num_userCommand >= 6)
    {
        string_to_integer(userCommands[5], m); // m
    }
    if (num_userCommand >= 7)
    {
        print("Usage:  nif  [num_nif] [RK_order] [ifAdj] [num_t_click] [m]\n");
        print("\t num_nif is number of neurons (1-100)\n");
        print("\t RK_order is 2 or 4\n");
        print("\t ifAdj is 1 or 0\n");
        print("\t num_t_click is an integer\n");
        print("\t m is an integer usually  1 <= m <= 5 \n");
        return false;
    }
    print("using user set parameters, i.e.,\n");
    return true;
}//}}}

boolean init_nif()//{{{
{
    integer i, j;
    for (i=0; i<num_nif; i++)
    {
        nif[i].ID = i;

        nif[i].v_prev = (1.0*rand())/RAND_MAX; // or V_leak?
        nif[i].t_prev = -t_step;
        nif[i].v  = nif[i].v_prev;
        nif[i].t  = 0; // always starts from time 0
        print("Init voltage of neuron[%ld]: %lf\n", i, nif[i].v);

        memset(nif[i].ls_t_spike, 0.0l, sizeof(nif[i].ls_t_spike));
        nif[i].lenLs_t_spike = 0;

         //[buffer overflow] may occur!!! vvv
        nif[i].datahub = NULL;
        sprint(myDatahubFilename, "%s%ld%s", myDatahubFilenameHead, i, myDatahubFilenameTail);
        nif[i].datahub = fopen(myDatahubFilename, "w");
        if (nif[i].datahub == NULL) DebugWrong("data file cannot open\n");

        if (nif[i].datahub) fprint(nif[i].datahub, "      %3.15lf", nif[i].t);
        if (nif[i].datahub) fprint(nif[i].datahub, "      %3.15lf", nif[i].v);
        if (nif[i].datahub) fprintln(nif[i].datahub);

        nif[i].Ge0_anglePosition = G_angleField * i / num_nif;
        nif[i].Gi0_anglePosition = G_angleField * i / num_nif + M_PI;

        // compute g_e0 and g_i0:
        nif[i].ge0_curr     = Ge0_intensity * sin(nif[i].t + nif[i].Ge0_anglePosition); // used to compute v_{t+1}
        nif[i].gi0_curr     = Gi0_intensity * sin(nif[i].t + nif[i].Gi0_anglePosition);
        nif[i].ge0_halfNext = Ge0_intensity * sin(nif[i].t + 0.5*t_step + nif[i].Ge0_anglePosition);
        nif[i].gi0_halfNext = Gi0_intensity * sin(nif[i].t + 0.5*t_step + nif[i].Gi0_anglePosition);
        nif[i].ge0_next     = Ge0_intensity * sin(nif[i].t + 1.0*t_step + nif[i].Ge0_anglePosition);
        nif[i].gi0_next     = Gi0_intensity * sin(nif[i].t + 1.0*t_step + nif[i].Gi0_anglePosition);

        // 1: compute G_e and G_i
        nif[i].Ge1_curr     = 0.0l;
        nif[i].Gi1_curr     = 0.0l;
        nif[i].Ge1_halfNext = 0.0l;
        nif[i].Gi1_halfNext = 0.0l;
        nif[i].Ge1_next     = 0.0l;
        nif[i].Gi1_next     = 0.0l;
        for (j=0; j<nif[i].lenLs_t_spike; j++)
        {
            nif[i].Ge1_curr     +=  G_sicc(nif[i].t - nif[i].ls_t_spike[j],     Tau_e, m);
            nif[i].Gi1_curr     +=  G_sicc(nif[i].t - nif[i].ls_t_spike[j],     Tau_i, m);
            nif[i].Ge1_halfNext +=  G_sicc(nif[i].t + 0.5*t_step - nif[i].ls_t_spike[j], Tau_e, m);
            nif[i].Gi1_halfNext +=  G_sicc(nif[i].t + 0.5*t_step - nif[i].ls_t_spike[j], Tau_i, m);
            nif[i].Ge1_next     +=  G_sicc(nif[i].t + 1.0*t_step - nif[i].ls_t_spike[j], Tau_e, m);
            nif[i].Gi1_next     +=  G_sicc(nif[i].t + 1.0*t_step - nif[i].ls_t_spike[j], Tau_i, m);
        }
    }

    for (i=0; i<num_nif; i++) // 2: compute g_e1 and g_i1
    {
        nif[i].ge1_curr     = 0.0l;
        nif[i].gi1_curr     = 0.0l;
        nif[i].ge1_halfNext = 0.0l;
        nif[i].gi1_halfNext = 0.0l;
        nif[i].ge1_next     = 0.0l;
        nif[i].gi1_next     = 0.0l;
        for (j = 0; j < num_nif; j++)
        {
            nif[i].ge1_curr     += mat_e[i][j] * nif[j].Ge1_curr    ;
            nif[i].gi1_curr     += mat_i[i][j] * nif[j].Gi1_curr    ;
            nif[i].ge1_halfNext += mat_e[i][j] * nif[j].Ge1_halfNext;
            nif[i].gi1_halfNext += mat_i[i][j] * nif[j].Gi1_halfNext;
            nif[i].ge1_next     += mat_e[i][j] * nif[j].Ge1_next    ;
            nif[i].gi1_next     += mat_i[i][j] * nif[j].Gi1_next    ;
        }
    }

    for (i=0; i<num_nif; i++)
    {
        nif[i].alpha_curr    =g_leak + (nif[i].ge0_curr     + nif[i].ge1_curr)    + (nif[i].gi0_curr     + nif[i].gi1_curr);
        nif[i].alpha_halfNext=g_leak + (nif[i].ge0_halfNext + nif[i].ge1_halfNext)+ (nif[i].gi0_halfNext + nif[i].gi1_halfNext);
        nif[i].alpha_next    =g_leak + (nif[i].ge0_next     + nif[i].ge1_next)    + (nif[i].gi0_next     + nif[i].gi1_next);
        nif[i].beta_curr     =g_leak*V_leak + (nif[i].ge0_curr + nif[i].ge1_curr)*V_E + (nif[i].gi0_curr + nif[i].gi1_curr)*V_I;
        nif[i].beta_halfNext =g_leak*V_leak + (nif[i].ge0_halfNext + nif[i].ge1_halfNext)*V_E + (nif[i].gi0_halfNext + nif[i].gi1_halfNext)*V_I;
        nif[i].beta_next     =g_leak*V_leak + (nif[i].ge0_next + nif[i].ge1_next)*V_E + (nif[i].gi0_next+ nif[i].gi1_next) *V_I;
    }
    return true;
}//}}}

boolean init_mat()//{{{
{
    integer i;
    for (i=0; i<num_nif; i++) // init the coupling matrix
    {
        integer j;
        for (j=0; j<num_nif; j++)
        {
            real ci = connection_intensity(i, j, num_nif);
            mat_e[i][j] = ci;
        }
    }
//    mat_i[i][j] = ci; // this matrix is kept to be 0s
    return true;
}//}}}

boolean initialization(int num_userCommand, char * userCommands[])//{{{
{
    boolean b_ret;
    RandomSeed(0);
    print("\nInitialization starts:\n");

    b_ret = init_para(num_userCommand, userCommands);
    if (!b_ret) return false;

    b_ret = init_nif();
    if (!b_ret) return false;

    b_ret = init_mat();
    if (!b_ret) return false;

    print("There are %ld neurons:\n", num_nif);
    print("\t integer num_t_click = %ld\n", num_t_click);
    print("\t integer RK_order = %ld\n", RK_order);
    print("\t integer ifAdj = %ld\n",  ifAdj);
    print("\t integer m = %ld\n", m);
    print("\t integer t_end = %ld\n", t_end);
    print("\t real t_step = %.15lf\n", t_step);
    println();

    print("For each neuron:\n");
    print("\t real g_leak     = %.15lf \n", g_leak);
    print("\t real V_leak     = %.15lf \n", V_leak);
    print("\t real V_threshold= %.15lf \n", V_threshold);
    print("\t real V_E        = %.15lf \n", V_E);
    print("\t real V_I        = %.15lf \n", V_I);
    print("\t real Tau_e      = %.15lf \n", Tau_e);
    print("\t real Tau_i      = %.15lf \n", Tau_i);
    print("\t real Ge0_intensity=%.15lf \n",Ge0_intensity);
    print("\t real i0_intensity=%.15lf \n",Gi0_intensity);
    print("\t real G_angleField =%.15lf \n",G_angleField);
    print("\t real coupleSpatial_s    = %.15lf \n",coupleSpatial_s);
    print("\t real coupleSpatial_mu   = %.15lf \n",coupleSpatial_mu);
    print("\t real coupleSpatial_delta= %.15lf \n",coupleSpatial_delta);
    println();

    if (num_nif > 6) return true; // do not print matrix, when there are too many neurons.
    print("The E coupling matrix is\n");
    integer i, j;
    for (i=0; i<num_nif; i++)
    {
        for (j=0; j<num_nif; j++)
        {
            print_real(mat_e[i][j]);
        }
        println();
    }
    println();

    print("The i coupling matrix is\n");
    for (i=0; i<num_nif; i++)
    {
        for (j=0; j<num_nif; j++)
        {
            print_real(mat_i[i][j]);
        }
        println();
    }
    println();

    return true;
}//}}}

boolean ppUpdate(integer i)//{{{
{
    pp.gl_RK_order = RK_order;
    pp.nif_num = nif[i].ID;
    pp.gl_g_leak = g_leak;
    pp.gl_V_leak = V_leak;
    pp.gl_V_E    = V_E;
    pp.gl_V_I    = V_I;
    pp.gl_V_threshold = V_threshold;
    pp.gl_t_step = t_step;

    pp.ge0_next          = nif[i].ge0_next;
    pp.ge1_next          = nif[i].ge1_next;
    pp.gi0_next          = nif[i].gi0_next;
    pp.gi1_next          = nif[i].gi1_next;

    pp.Ge1_curr          = nif[i].Ge1_curr;
    pp.Gi1_curr          = nif[i].Gi1_curr;
    pp.Ge1_halfNext      = nif[i].Ge1_halfNext;
    pp.Gi1_halfNext      = nif[i].Gi1_halfNext;

    pp.ge0_curr          = nif[i].ge0_curr;
    pp.gi0_curr          = nif[i].gi0_curr;
    pp.ge0_halfNext      = nif[i].ge0_halfNext;
    pp.gi0_halfNext      = nif[i].gi0_halfNext;

    pp.ge1_curr          = nif[i].ge1_curr;
    pp.gi1_curr          = nif[i].gi1_curr;
    pp.ge1_halfNext      = nif[i].ge1_halfNext;
    pp.gi1_halfNext      = nif[i].gi1_halfNext;

    pp.t_prev            = nif[i].t_prev;
    pp.t                 = nif[i].t;
    pp.v_prev            = nif[i].v_prev;
    pp.v                 = nif[i].v;
    pp.last_t_spike      = 0.0l;
    if (nif[i].lenLs_t_spike > 0)   pp.last_t_spike = nif[i].ls_t_spike[nif[i].lenLs_t_spike-1];

    pp.alpha_curr        = nif[i].alpha_curr;
    pp.alpha_halfNext    = nif[i].alpha_halfNext;
    pp.alpha_next        = nif[i].alpha_next;
    pp.beta_curr         = nif[i].beta_curr;
    pp.beta_halfNext     = nif[i].beta_halfNext;
    pp.beta_next         = nif[i].beta_next;
    return true;
}//}}}

boolean deconstruction()//{{{
{
    integer i;
    for (i = 0; i < num_nif; i++) if (nif[i].datahub) fclose(nif[i].datahub);
    for (i = 0; i < num_nif; i++) print("\nNeuron #%ld @ Time %3.15lf : Vol. %3.15lf", i, nif[i].t, nif[i].v); //time and vol at the last step
    print("\nDeconstruction finished.\n");
    return true;
} //}}}
//}}}

StartFromXOS
    integer  ijk, iii, jjj;
    for (ijk = 1; ijk <= num_t_click; ijk ++)//{{{
    { // main loop --- --- ... ...
//         print("Time %d / %ld\n", ijk, num_t_click);

        for (iii = 0; iii < num_nif; iii ++)//{{{
        { // loop of the current neuron
            nif[iii].v_prev = nif[iii].v;
            ppUpdate(iii);
            if   (RK_order == 2)    nif[iii].v = RK2_stepper(dynamics, nif[iii].t, nif[iii].v, t_step, pp);
            elif (RK_order == 4)    nif[iii].v = RK4_stepper(dynamics, nif[iii].t, nif[iii].v, t_step, pp);
            else                    DebugError("only RK2 or RK4 are allowed!!",false);
            nif[iii].t_prev = nif[iii].t;
            nif[iii].t += t_step;

            if (nif[iii].v >= V_threshold)
            { // adjust time and voltage now //{{{
                if (ifAdj)
                {
                    ppUpdate(iii);
                    real t_spikeDelta = adj_t_spikeDelta(pp);
                    nif[iii].ls_t_spike[nif[iii].lenLs_t_spike] = nif[iii].t_prev + t_spikeDelta;
                    nif[iii].lenLs_t_spike += 1;
                    ppUpdate(iii);
                    nif[iii].v = adj_v(pp);
                }
               else
                {
                    nif[iii].ls_t_spike[nif[iii].lenLs_t_spike] = nif[iii].t;
                    nif[iii].lenLs_t_spike += 1;
                    nif[iii].v = nif[iii].v - (V_threshold - V_leak);
                }
            }//}}}

            if (nif[iii].datahub)
            {
                fprint(nif[iii].datahub, "      %3.15lf", nif[iii].t);
                fprint(nif[iii].datahub, "      %3.15lf", nif[iii].v);
                fprintln(nif[iii].datahub);
             }
        } // the current neuron finished//}}}

        // process the couping here...  //{{{
        //// 1) neuron computes G_e, // captial G_e is a factor of g_e1
        //// 2) main() computes g_e1 with the use of G_e and e_couple_matrix,
        //// 3) all drive information are obtained (g_e0 can be computed according to the current time),
        //      and dynamics() can be called to compute the next v again.

        // compute g_e0 and g_i0:
        for (iii = 0; iii < num_nif; iii++)
        {
            nif[iii].ge0_curr     = Ge0_intensity * sin(nif[iii].t + nif[iii].Ge0_anglePosition);
            nif[iii].gi0_curr     = Gi0_intensity * sin(nif[iii].t + nif[iii].Gi0_anglePosition);
            nif[iii].ge0_halfNext = Ge0_intensity * sin(nif[iii].t + 0.5*t_step + nif[iii].Ge0_anglePosition);
            nif[iii].gi0_halfNext = Gi0_intensity * sin(nif[iii].t + 0.5*t_step + nif[iii].Gi0_anglePosition);
            nif[iii].ge0_next     = Ge0_intensity * sin(nif[iii].t + 1.0*t_step + nif[iii].Ge0_anglePosition);
            nif[iii].gi0_next     = Gi0_intensity * sin(nif[iii].t + 1.0*t_step + nif[iii].Gi0_anglePosition);
        }

        // compute g_e1 and g_i1 (1: compute G_e and G_i; 2: compute g_e1 and g_i1):
        // 1: compute G_e and G_i
        for (iii = 0; iii < num_nif; iii++)
        {
            nif[iii].Ge1_curr     = 0.0l;
            nif[iii].Gi1_curr     = 0.0l;
            nif[iii].Ge1_halfNext = 0.0l;
            nif[iii].Gi1_halfNext = 0.0l;
            nif[iii].Ge1_next     = 0.0l;
            nif[iii].Gi1_next     = 0.0l;
            for (jjj=0; jjj<nif[iii].lenLs_t_spike; jjj++)
            {
                nif[iii].Ge1_curr     +=  G_sicc(nif[iii].t - nif[iii].ls_t_spike[jjj],     Tau_e, m);
                nif[iii].Gi1_curr     +=  G_sicc(nif[iii].t - nif[iii].ls_t_spike[jjj],     Tau_i, m);
                nif[iii].Ge1_halfNext +=  G_sicc(nif[iii].t + 0.5*t_step - nif[iii].ls_t_spike[jjj], Tau_e, m);
                nif[iii].Gi1_halfNext +=  G_sicc(nif[iii].t + 0.5*t_step - nif[iii].ls_t_spike[jjj], Tau_i, m);
                nif[iii].Ge1_next     +=  G_sicc(nif[iii].t + 1.0*t_step - nif[iii].ls_t_spike[jjj], Tau_e, m);
                nif[iii].Gi1_next     +=  G_sicc(nif[iii].t + 1.0*t_step - nif[iii].ls_t_spike[jjj], Tau_i, m);
            }
        }

        // 2: compute g_e1 and g_i1
        for (iii = 0; iii < num_nif; iii++)
        { // Compute C = A B ;; A:m*p; B:p*n; C:m*n
            nif[iii].ge1_curr     = 0.0l;
            nif[iii].gi1_curr     = 0.0l;
            nif[iii].ge1_halfNext = 0.0l;
            nif[iii].gi1_halfNext = 0.0l;
            nif[iii].ge1_next     = 0.0l;
            nif[iii].gi1_next     = 0.0l;
            for (jjj = 0; jjj < num_nif; jjj++)
            {
                nif[iii].ge1_curr     += mat_e[iii][jjj] * nif[jjj].Ge1_curr     ;
                nif[iii].gi1_curr     += mat_i[iii][jjj] * nif[jjj].Gi1_curr     ;
                nif[iii].ge1_halfNext += mat_e[iii][jjj] * nif[jjj].Ge1_halfNext ;
                nif[iii].gi1_halfNext += mat_i[iii][jjj] * nif[jjj].Gi1_halfNext ;
                nif[iii].ge1_next     += mat_e[iii][jjj] * nif[jjj].Ge1_next     ;
                nif[iii].gi1_next     += mat_i[iii][jjj] * nif[jjj].Gi1_next     ;
            }
        }

        for (iii = 0; iii < num_nif; iii++)
        {
            nif[iii].alpha_curr    =g_leak + (nif[iii].ge0_curr     + nif[iii].ge1_curr)    + (nif[iii].gi0_curr     + nif[iii].gi1_curr);
            nif[iii].alpha_halfNext=g_leak + (nif[iii].ge0_halfNext + nif[iii].ge1_halfNext)+ (nif[iii].gi0_halfNext + nif[iii].gi1_halfNext);
            nif[iii].alpha_next    =g_leak + (nif[iii].ge0_next     + nif[iii].ge1_next)    + (nif[iii].gi0_next     + nif[iii].gi1_next);
            nif[iii].beta_curr     =g_leak*V_leak + (nif[iii].ge0_curr + nif[iii].ge1_curr)*V_E+ (nif[iii].gi0_curr + nif[iii].gi1_curr)*V_I;
            nif[iii].beta_halfNext =g_leak*V_leak + (nif[iii].ge0_halfNext + nif[iii].ge1_halfNext)*V_E+ (nif[iii].gi0_halfNext + nif[iii].gi1_halfNext)*V_I;
            nif[iii].beta_next     =g_leak*V_leak + (nif[iii].ge0_next + nif[iii].ge1_next)*V_E+ (nif[iii].gi0_next + nif[iii].gi1_next) *V_I;
        }
        //}}}

    }// main loop ended   }}}
EndWithHonor

integer connection_distance(integer i, integer j, integer num_nif_ring)//{{{ {{{
{
    // This function computes how long is the distance of two nodes on a same ring.
    // i and j are node numbers, num_nif_ring are how many nodes on the ring.
    // the node numbers should be labeled as 0, 1, 2, ..., num_nif_ring-1...
    if (i<0 || i>=num_nif_ring || j<0 || j>=num_nif_ring) return false;
    return min( labs(i-j),  num_nif_ring + min(i,j) - max(i,j) );
}//}}}

real connection_intensity(integer i, integer j, integer num_nif_ring)//{{{
{
    // This function computes the connection intensity between two nodes on a ring.
    // the intensity damps following the Gauss's low with an extra strength parameter
    // i and j are node numbers, num_nif_ring are how many nodes on the ring.
    // mu:mean deviation (usually 0); delta:standard deviation; s:strength
    integer distance = 0;
    distance = connection_distance(i, j, num_nif_ring);
    if (distance == 0) return 0.0; // a neuron does not affect itself
    real intensity = 0.0l;
    Debug0_Div(coupleSpatial_delta);
    intensity = coupleSpatial_s*exp(-0.5*pow((distance-coupleSpatial_mu)/coupleSpatial_delta, 2))/sqrt(2*M_PI*coupleSpatial_delta);
    if (intensity > 1.0) intensity = 1.0;
    return intensity;
}//}}}

real RK2_stepper(real(*df)(real, real, paraPool), real t0, real x0, real ts, paraPool pp)//{{{
{
    real K1, K2;
    K1 = df(t0     , x0        , pp);
    K2 = df(t0 + ts, x0 + ts*K1, pp);
    return x0 + ts * (K1 + K2) / 2.0;
}//}}}

real RK4_stepper(real(*df)(real, real, paraPool), real t0, real x0, real ts, paraPool pp)//{{{
{
    real K1, K2, K3, K4;
    K1 = df(t0         , x0            , pp);
    K2 = df(t0 + ts/2.0, x0 + ts*K1/2.0, pp);
    K3 = df(t0 + ts/2.0, x0 + ts*K2/2.0, pp);
    K4 = df(t0 + ts    , x0 + ts*K3    , pp);
    return x0 + ts * (K1 + 2*K2 + 2*K3 + K4) / 6.0;
}//}}}

real dynamics(real t, real v, paraPool pp) //{{{
{ // the membrane potential dynamics, i.e., Eq.(4)
    // dv/dt = -g_leak(v-V_leak)   -g_e(t)(v-V_E)   -g_i(t)(v-V_I)
    // do not use the form of g_e(t), g_i(t); because it will refer to global vars
    if   EQ(t, pp.t)                  return -pp.gl_g_leak * (v-pp.gl_V_leak) - (pp.ge0_curr    +pp.ge1_curr)     * (v-pp.gl_V_E)- (pp.gi0_curr    +pp.gi1_curr) * (v-pp.gl_V_I);
    elif EQ(t, pp.t+0.5*pp.gl_t_step) return -pp.gl_g_leak * (v-pp.gl_V_leak) - (pp.ge0_halfNext+pp.ge1_halfNext) * (v-pp.gl_V_E)- (pp.gi0_halfNext+pp.gi1_halfNext) * (v-pp.gl_V_I);
    elif EQ(t, pp.t + pp.gl_t_step)   return -pp.gl_g_leak * (v-pp.gl_V_leak) - (pp.ge0_next    +pp.ge1_next)     * (v-pp.gl_V_E)- (pp.gi0_next    +pp.gi1_next) * (v-pp.gl_V_I);
    else DebugError("a strange step (length)",false);
} //}}}

real adj_t_spikeDelta(paraPool pp)//{{{
{ // linera interpolation codes for RK2
    real t_spikeDelta;
    if (pp.gl_V_threshold < pp.v_prev) DebugWrong("the last voltage > V_threshold!");
    if (pp.v < pp.v_prev) DebugWrong("the current voltage < the last one!");
    t_spikeDelta = pp.gl_t_step * (pp.gl_V_threshold - pp.v_prev)/(pp.v - pp.v_prev);
    if   (pp.gl_RK_order == 2) return t_spikeDelta;// vvv RK4 use NewtonMethod() and search from t_ret; DichotomySearch() also works well.
    elif (pp.gl_RK_order == 4) return NewtonMethod(TimeCubicInterpolation, TimeCubicInterpolationDerivative, t_spikeDelta, pp);
    else DebugError("only RK2 or RK4 are allowed!!", false);
}//}}}

real NewtonMethod(real(*eq)(real, paraPool), real(*eqDerivative)(real, paraPool), real t0, paraPool pp)//{{{
// find a root t_spikeDelta (t_spikeDelta = t_spike - t_n) (near t0) for eq
// eq is the equation to solve,
// eqDerivative is eq's Derivative,
// t0 is the point from where we begin to search.
{
    integer iii;
    real t_iteration = t0;
    for (iii=0; iii<10; iii++)       t_iteration -= eq(t_iteration, pp)/eqDerivative(t_iteration, pp);
    return t_iteration;
}//}}}

real TimeCubicInterpolation(real t_spikeDelta, paraPool pp)//{{{
// eq(15) with left hand v(t) moved to right and set to V_threshold
// this is used in NewtonMethod() to adjust t_spike.
// ---- can also used with DichotomySearch() to adjust t_spike.
{
    real v_n = pp.v_prev;
    real v_nPlus1 = pp.v;
    real vPrime_n = pp.beta_curr - pp.alpha_curr*v_n;
    real vPrime_nPlus1 = pp.beta_next - pp.alpha_next*v_nPlus1;
    return v_n + vPrime_n*t_spikeDelta - pp.gl_V_threshold
        + ( 3*(v_nPlus1-v_n) - pp.gl_t_step * (2*vPrime_n + vPrime_nPlus1))*pow(t_spikeDelta/pp.gl_t_step, 2)
        + (-2*(v_nPlus1-v_n) + pp.gl_t_step * (  vPrime_n + vPrime_nPlus1))*pow(t_spikeDelta/pp.gl_t_step, 3);
}//}}}

real TimeCubicInterpolationDerivative(real t_spikeDelta, paraPool pp)//{{{
// Derivation of eq(15), i.e., eq(15)'t, (t is t_spikeDelta here)
// this is used with NewtonMethod() to adjust t_spike.
{
    real v_n = pp.v_prev;
    real v_nPlus1 = pp.v;
    real vPrime_n = pp.beta_curr - pp.alpha_curr*v_n;
    real vPrime_nPlus1 = pp.beta_next - pp.alpha_next*v_nPlus1;
    return vPrime_n
        + ( 3*(v_nPlus1 - v_n) - pp.gl_t_step * (2*vPrime_n + vPrime_nPlus1))*2*    t_spikeDelta    /pow(pp.gl_t_step, 2)
        + (-2*(v_nPlus1 - v_n) + pp.gl_t_step * (  vPrime_n + vPrime_nPlus1))*3*pow(t_spikeDelta, 2)/pow(pp.gl_t_step, 3);
}//}}}

real adj_v(paraPool pp)//{{{
{
    if   (pp.gl_RK_order != 2 && pp.gl_RK_order != 4) DebugError("only RK2 or RK4 are allowed!!", false);

    real v_nPlus1, v_n, t_spikeDelta, K1, K2, v_Delta;
    t_spikeDelta = pp.last_t_spike - pp.t_prev;
    v_n = (2.0 * pp.gl_V_leak - t_spikeDelta * ( pp.beta_next +  pp.beta_next -  pp.beta_curr * pp.alpha_next * pp.gl_t_step) ) //V-leak ...
                    / (2.0 - t_spikeDelta * (   pp.alpha_next + pp.alpha_next - pp.alpha_curr * pp.alpha_next * pp.gl_t_step) ); ////???? WRONG...
    K1 = pp.beta_curr - pp.alpha_curr *  v_n;
    K2 = pp.beta_next - pp.alpha_next * (v_n + K1 * pp.gl_t_step);
    v_Delta = pp.gl_t_step * ( K1 + K2 ) * 0.5;
    v_nPlus1 = v_n + v_Delta;

    if (pp.gl_RK_order == 2) return v_nPlus1;

    real K3, K4, v_spike;

    v_spike = VoltageCubicInterpolation(0, pp);
    v_n = pp.gl_V_leak - (pp.v_prev - pp.gl_V_leak) * (v_spike - pp.gl_V_leak)/(pp.gl_V_threshold - v_spike);
    K1 = pp.beta_curr     - pp.alpha_curr     * (v_n + 0                      );
    K2 = pp.beta_halfNext - pp.alpha_halfNext * (v_n + 0.5 * K1 * pp.gl_t_step);
    K3 = pp.beta_halfNext - pp.alpha_halfNext * (v_n + 0.5 * K2 * pp.gl_t_step);
    K4 = pp.beta_next     - pp.alpha_next     * (v_n +       K3 * pp.gl_t_step);
    v_Delta = pp.gl_t_step * (K1 + 2*K2 + 2*K3 + K4)/6.0;

    v_nPlus1 = v_n + v_Delta;
    DebugCheck(v_nPlus1<0, "triple adjusted volume less than 0...");
    return v_nPlus1;
} //}}}

real VoltageCubicInterpolation(real v_n, paraPool pp)//{{{
{
    // eq(15)=0 with left hand v(t) moved to right and set to V_leak
    real t_n= pp.t_prev;
    real K1 = pp.beta_curr     - pp.alpha_curr     * (v_n + 0                      );
    real K2 = pp.beta_halfNext - pp.alpha_halfNext * (v_n + 0.5 * K1 * pp.gl_t_step);
    real K3 = pp.beta_halfNext - pp.alpha_halfNext * (v_n + 0.5 * K2 * pp.gl_t_step);
    real K4 = pp.beta_next     - pp.alpha_next     * (v_n +       K3 * pp.gl_t_step);
    real v_Delta  = pp.gl_t_step * (K1 + 2*K2 + 2*K3 + K4) / 6.0;
    real vPrime_n = pp.beta_curr - pp.alpha_curr*v_n;
    real vPrime_nPlus1 = pp.beta_next - pp.alpha_next*(v_n + v_Delta);
    real t_spikeDelta  = pp.last_t_spike - t_n;
    // NOTE: according to eq(14), we have v_nPlus1-v_n = (K1+2*K2+2*K3+k4)*t_step/6;
    // by substituting (v_t-v_n) into eq(15), then we get:
    return v_n + vPrime_n*t_spikeDelta - pp.gl_V_leak
        + ( 3*v_Delta - pp.gl_t_step * (2*vPrime_n + vPrime_nPlus1))*pow(t_spikeDelta/pp.gl_t_step, 2)
        + (-2*v_Delta + pp.gl_t_step * (  vPrime_n + vPrime_nPlus1))*pow(t_spikeDelta/pp.gl_t_step, 3);
}//}}}





