//!/usr/bin/env c
//-*- coding:utf-8 -*-

//=== README ===================================/*{{{*/
//·
//· Author: Mogei Wang
//·
//· MogeiWang@GMail.com
//·
//· Filename: stack.c
//·
//· COPYRIGHT 2015
//·
//· Description:
//· gcc -std=gnu99 stack.c -o stack -g -lm -lgsl -lgslcblas
//·
//· There may be small errors in DichotomySearch(), chech again!!!
//·
//==============================================/*}}}*/

// --- includes and macros {{{
#include <stdio.h>
#include <math.h> //tgmath.h
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <GL/freeglut.h>

#define BEGIN {
#define END }
#define integer long int
#define boolean integer
#define real double
#define elif else if
#define true (1==1)   // gcc: 1=true, 0=false
#define false (!true)
#define IdentifierLiteral(T) #T // return "T"
#define IdentifierJoint(T,x) T##x // return Tx
#define print printf
#define fprint fprintf
#define pprint fprintf
#define DebugWrong(s)   {printf("\n[LINE %d] [FILe %s] [WRONG]: %s\n",__LINE__,__FILE__, s);}
#define DebugW_T_F(s)   {printf("\n[LINE %d] [FILe %s] [W-T-F]: %s\n",__LINE__,__FILE__, s);}
#define DebugError(s, t)  {printf("\n[LINE %d] [FILe %s] [ERROR]: %s\nreturning...\n",__LINE__,__FILE__, s); return t;}
#define Debug0_Div(divisor) {if(divisor==0) printf("\n[LINE %d] [FILe %s] [ERROR]: ATTEMPT TO USE 0 AS DENOMINATOR!\n",__LINE__,__FILE__);}
#define StartFromXOS  int main(int argc, char * argv[]){ if (!initialization(argc, argv)) ErrorReturn("program cannot be initializated!", false);
#define EndWithHonor  return deconstruction(); }
//}}}

// --- system {{{
#define t_end (1.0l)
#define t_click_number_default (1024000L)
#define RK_order_default (4L)
#define adjust_mode_default (true)

integer t_click_number = t_click_number_default; // how many time clicks (number of time steps)
integer RK_order=RK_order_default; // RK2 or RK4. both of them use the linear interpolation, so both errors should be of the order Dt^2.
boolean adjust_mode = adjust_mode_default; //true; // use the linear interpolation or not.

real t_step=t_end/t_click_number_default; // compute time step according to the end time and number of time clicks (start time = 0)
real t_list[1+t_click_number_default]; // from 0 to current time
real v_list[1+t_click_number_default]; // this array is pretty large, usually do not use this large...
real t_spike_list[t_click_number_default]; // some memory is wasted here :-)
integer len_t_list, len_v_list, len_t_spike_list;

real g_leak = 50.0; // --- for single neurons
real V_leak = 0.0; // the V parameters
real V_E = 14/3.0;
real V_I = -2/3.0;
real V_threshold = 1.0;

FILE * plotter;
FILE * datahub;
char * myPlotterFilename;
char * myDatahubFilename;
//}}}

// functions {{{
real RK2_stepper(real(*df)(real,real),real x0,real t0,real ts);
real RK4_stepper(real(*df)(real,real),real x0,real t0,real ts);
real neuron(real t, real v);
real alpha(int n);
real beta(int n);
real adjust_t_spike_delta();
real adjust_v_n();
real DichotomySearch(real(*eq)(real), real x0, real x1, int deepth);
real VoltageCubicInterpolation(real v_n);
real TimeCubicInterpolation(real t_spike_delta);
real TimeCubicInterpolationDerivative(real t_spike_delta);
real NewtonRaphsonMethod(real(*eq)(real), real(*eqDerivative)(real), real t0);
// }}}

boolean initialization(int user_command_number, char * user_commands[])/*{{{*/
{
    print("\nInitialization starts:\n");
    if (user_command_number <= 1)/*{{{*/
    {// use default
        print("using program default parameters, i.e.,\n");
        print("integer t_click_number = %ld\n", t_click_number_default);
        print("integer RK_order = %ld\n", RK_order_default);
        print("integer adjust_mode = %d\n",  adjust_mode_default);
    }
    elif (user_command_number != 4)
    {
        print("Usage: singleNeuron RK_order adjust_mode t_click_number\n");
        print("\t RK_order is 2 or 4\n");
        print("\t adjust_mode is 1 or 0\n");
        print("\t t_click_number is an integer\n");
        return false;
    }
    else
    {// use user_commands: RK_order, adjust_mode, t_click_number
        RK_order = (((atol(user_commands[1])-2)/2)%2)*2+2; // RK2 or RK4. (map any integer to 2,4, and keep 2->2, 4->4 unchanged)
        adjust_mode = atol(user_commands[2]) %2; //use the linear interpolation or not. In gcc, 1 is true, and 0 is false.
        t_click_number = atol(user_commands[3]); // how many time clicks (number of time steps)
        // user_commands[0] is the filename of the running command.

        print("using user set parameters, i.e.,\n");
        print("integer t_click_number = %ld\n", t_click_number);
        print("integer RK_order = %ld\n", RK_order);
        print("integer adjust_mode = %ld\n",  adjust_mode);
        t_step=t_end/t_click_number; // compute time step according to the end time and number of time clicks (start time = 0)
    }/*}}}*/
    { // init the lists/*{{{*/
        for (int i=0; i<=t_click_number; i++)
        {
            t_list[i]=0;
            v_list[i]=0;
        }
        v_list[0]=V_leak;
        len_v_list=1;
        t_list[0]=0;
        len_t_list=1;
    }/*}}}*/
    { // init the gnuplot things/*{{{*/
        plotter = popen("gnuplot","w");
        if (plotter==NULL) ErrorReturn("gnuplot pipe cannot open\n", false);
        pprint(plotter,"set terminal post enh \n");
        myPlotterFilename="neuron_voltage.pdf";
        myDatahubFilename="neuron_voltage.txt";
        datahub = fopen(myDatahubFilename, "w+");
        if (datahub==NULL) ErrorReturn("data file cannot open\n", false);
    }/*}}}*/
    return true;
}/*}}}*/

boolean deconstruction()/*{{{*/
{
    pprint(plotter, "set output\n");
    pprint(plotter, "set term wxt\n");
    pprint(plotter, "\n");
    pclose(plotter);
    fclose(datahub);
    print("\nDeconstruction finished.\n");
    return true;
} /*}}}*/

StartFromXOS
    for (int i=1; i<=t_click_number; i++) // main loop --- --- ... .../*{{{*/
    {
        if   (RK_order==2) v_list[len_v_list] = RK2_stepper(neuron, v_list[len_v_list-1], t_list[len_t_list-1], t_step); // not matter RK2 or RK4,
        elif (RK_order==4) v_list[len_v_list] = RK4_stepper(neuron, v_list[len_v_list-1], t_list[len_t_list-1], t_step); // linear interpolation is used.
        else ErrorReturn("only RK2 or RK4 are allowed!!",-1.0l);
        len_v_list ++;
        t_list[len_t_list] = t_list[len_t_list-1]+t_step; // the previous error was put this before v!!!
        len_t_list ++;                                    // this (t) and v can be exchanged, but not to change v to t_list[len_v_list-2]!!!
        if (v_list[len_v_list-1] >= V_threshold)
        { // adjust time and voltage now /*{{{*/
            if (adjust_mode)
            {
                real t_spike_delta = adjust_t_spike_delta();
                t_spike_list[len_t_spike_list] = t_list[len_t_list-2] + t_spike_delta;
                len_t_spike_list ++;
                v_list[len_v_list-1] = adjust_v_n();
            }
            else
            {
                t_spike_list[len_t_spike_list] = t_list[len_t_list-1];
                len_t_spike_list ++;
                v_list[len_v_list-1] -= (V_threshold - V_leak);
            }
        }/*}}}*/
    }/*}}}*/
    { // print the spike list/*{{{*/
        print("There are %ld spkikes, and \n", len_t_spike_list);
        print("The spkike time list is: ");
        for (int i=0; i<len_t_spike_list; i++) print("%3.15lf ", t_spike_list[i]);
        print("\n\n");
    }/*}}}*/
    { // plot/*{{{*/
        if (len_t_list != len_t_list) ErrorReturn("length of time and v is not the same.", false);
        print("saving v_list to datahub...\n");
        for (int i=0; i<len_t_list; i++)
        {
            fprint(datahub, "      %3.15lf", t_list[i]);
            fprint(datahub, "      %3.15lf", v_list[i]);
            fprint(datahub, "\n");
        }
        { // print the voltage at the last click.
            print("\n");
            print("Time: %3.15lf", t_list[len_t_list-1]);
            print(" Vol: %3.15lf", v_list[len_t_list-1]);
            print("\n");
        }
        print("done.\n");
        pprint(plotter, "set output '%s' \n", myPlotterFilename);
        pprint(plotter, "plot \"%s\" using 1:2 with l lc 1  lw 1 title \"\"\n", myDatahubFilename);

        print("\n");
    }/*}}}*/
EndWithHonor

real RK2_stepper(real(*df)(real,real),real x0,real t0,real ts)
{/*{{{*/
    real K1, K2;
    K1 = df(t0     , x0        );
    K2 = df(t0 + ts, x0 + ts*K1);
    return x0 + ts * (K1 + K2) / 2.0;
}/*}}}*/

real RK4_stepper(real(*df)(real,real),real x0,real t0,real ts)
{/*{{{*/
    real K1, K2, K3, K4;
    K1 = df(t0         , x0            );
    K2 = df(t0 + ts/2.0, x0 + ts*K1/2.0);
    K3 = df(t0 + ts/2.0, x0 + ts*K2/2.0);
    K4 = df(t0 + ts    , x0 + ts*K3    );
    return x0 + ts * (K1 + 2*K2 + 2*K3 + K4) / 6.0;
}/*}}}*/

real neuron(real t, real v)
{/*{{{*/
    return -g_leak * (v - V_leak) - 25.0 * sin(t) * (v - V_E); //- g_i(t) * (v - V_I)
}/*}}}*/

real alpha(int n)
{/*{{{*/
    if (n==0) return g_leak + 25.0 * sin(t_list[len_t_list-2]);
    if (n==1) return g_leak + 25.0 * sin(t_list[len_t_list-1]);
    if (n==2) return g_leak + 25.0 * sin(t_list[len_t_list-2] + 0.5*t_step); // this returns alphs(1/2) !!!;
}/*}}}*/

real beta(int n)
{/*{{{*/
    if (n==0) return g_leak * V_leak + 25.0 * sin(t_list[len_t_list-2]) * V_E;
    if (n==1) return g_leak * V_leak + 25.0 * sin(t_list[len_t_list-1]) * V_E;
    if (n==2) return g_leak * V_leak + 25.0 * sin(t_list[len_t_list-2] + 0.5*t_step) * V_E; // this returns beta(1/2) !!!
}/*}}}*/

real adjust_t_spike_delta()
{/*{{{*/
    // linera interpolation codes for RK2
    real t_spike_delta;
    if (V_threshold < v_list[len_v_list-2]) ErrorPrint("the last voltage > V_threshold!");
    if (v_list[len_v_list-1] < v_list[len_v_list-2])
    {
        ErrorPrint("the current voltage < the last one!");
        print("\t the current voltage is %lf ; while the last one %3.15lf!", v_list[len_v_list-1], v_list[len_v_list-2]);
    }
    t_spike_delta = t_step * (V_threshold - v_list[len_v_list-2]) / (v_list[len_v_list-1] - v_list[len_v_list-2]);
    if   (RK_order == 2) return t_spike_delta;
    elif (RK_order == 4) // use NewtonRaphsonMethod() and search from t_ret; DichotomySearch() also works well, but need to be more careful.
                         return NewtonRaphsonMethod(TimeCubicInterpolation, TimeCubicInterpolationDerivative, t_spike_delta);
    else ErrorReturn("only RK2 or RK4 are allowed!!",-1.0l);
}/*}}}*/

real adjust_v_n()
{/*{{{*/
    if (RK_order == 4)
    {
        real v_n, K1, K2, K3, K4, v_Delta, v_nPlus1;
        v_n = DichotomySearch(VoltageCubicInterpolation, V_leak-0.5, V_leak, 0); //print("Careful, using DichotomySearch() to find v_n: %lf ", v_n);
        K1 = beta(0) - alpha(0) * (v_n + 0                );
        K2 = beta(2) - alpha(2) * (v_n + 0.5 * K1 * t_step);
        K3 = beta(2) - alpha(2) * (v_n + 0.5 * K2 * t_step);
        K4 = beta(1) - alpha(1) * (v_n +       K3 * t_step);
        v_Delta = t_step * (K1 + 2*K2 + 2*K3 + K4)/6;
        v_nPlus1 = v_n + v_Delta;
        return v_nPlus1;
    }
    elif (RK_order != 2) ErrorReturn("only RK2 or RK4 are allowed!!",-1.0l);
    // the following are linera interpolation codes for RK2
    real v_nPlus1, v_n, t_spike_delta, K1, K2, t_n;
    t_n = t_list[len_t_list-2];
    t_spike_delta = t_spike_list[len_t_spike_list-1] - t_n;
    v_n = (2 * V_leak - t_spike_delta * (  beta(0) + beta(1) - alpha(1) * beta(0) * t_step) )
                 / (2 - t_spike_delta * ( alpha(0) + alpha(1) - alpha(1) * alpha(0) * t_step) );
    K1 = beta(0) - alpha(0) * v_n;
    K2 = beta(1) - alpha(1) * (v_n + K1 * t_step);
    v_nPlus1 = v_n + t_step * ( K1 + K2 ) / 2.0;
    return v_nPlus1;
} /*}}}*/

real VoltageCubicInterpolation(real v_n)
// eq(15)=0 with left hand v(t) moved to right and set to V_leak/*{{{*/
// this is used with DichotomySearch() to compute v_n,
// once v_n is got, v_nPlus1 can be obtained by eq(14)...
{
    real t_n = t_list[len_t_list-2];
    real K1 = beta(0) - alpha(0) * (v_n + 0                );
    real K2 = beta(2) - alpha(2) * (v_n + 0.5 * K1 * t_step);
    real K3 = beta(2) - alpha(2) * (v_n + 0.5 * K2 * t_step);
    real K4 = beta(1) - alpha(1) * (v_n +       K3 * t_step);
    real v_Delta = t_step * (K1 + 2*K2 + 2*K3 + K4) / 6;
    real v_prime_n = beta(0) - alpha(0)*v_n;
    real v_prime_nPlus1 = beta(1) - alpha(1)*(v_n+v_Delta);
    real t_spike_delta = t_spike_list[len_t_spike_list-1] - t_n;
    // NOTE: according to eq(14), we have v_nPlus1-v_n = (K1+2*K2+2*K3+k4)*t_step/6;
    // by substituting (v_nPlus1-v_n) into eq(15), then we get:
    return v_n + v_prime_n*t_spike_delta - V_leak
        + ( 3*v_Delta - t_step * (2*v_prime_n + v_prime_nPlus1))*pow(t_spike_delta/t_step, 2)
        + (-2*v_Delta + t_step * (  v_prime_n + v_prime_nPlus1))*pow(t_spike_delta/t_step, 3);
}/*}}}*/

real TimeCubicInterpolation(real t_spike_delta)
// eq(15) with left hand v(t) moved to right and set to V_threshold/*{{{*/
// this is used with DichotomySearch() to adjust t_spike.
{
    real v_n = v_list[len_v_list-2];
    real v_nPlus1 = v_list[len_v_list-1];
    real v_prime_n = beta(0) - alpha(0)*v_n;
    real v_prime_nPlus1 = beta(1) - alpha(1)*v_nPlus1;
    real t_n = t_list[len_t_list-2];
    return v_n + v_prime_n*t_spike_delta - V_threshold
        + ( 3*(v_nPlus1-v_n) - t_step * (2*v_prime_n + v_prime_nPlus1))*pow(t_spike_delta/t_step, 2)
        + (-2*(v_nPlus1-v_n) + t_step * (  v_prime_n + v_prime_nPlus1))*pow(t_spike_delta/t_step, 3);
}/*}}}*/

real TimeCubicInterpolationDerivative(real t_spike_delta)
// Derivation of eq(15), i.e., eq(15)'t, (t is t_spike_delta here) /*{{{*/
// this is used with NewtonRaphsonMethod() to adjust t_spike.
{
    real v_n = v_list[len_v_list-2];
    real v_nPlus1 = v_list[len_v_list-1];
    real v_prime_n = beta(0) - alpha(0)*v_n;
    real v_prime_nPlus1 = beta(1) - alpha(1)*v_nPlus1;
    real t_n = t_list[len_t_list-2];
    return v_prime_n
        + ( 3*(v_nPlus1 - v_n) - t_step * (2*v_prime_n + v_prime_nPlus1))*2*t_spike_delta/pow(t_step, 2)
        + (-2*(v_nPlus1 - v_n) + t_step * (  v_prime_n + v_prime_nPlus1))*3*pow(t_spike_delta, 2)/pow(t_step, 3);
}/*}}}*/

real DichotomySearch(real(*eq)(real), real x0, real x1, int deepth)
// find a root x for eq, x0<x<x1 OR x1<x<x0/*{{{*/
// eq is the equation to solve,
// x0 and x1 are two ends of the region the root is in.
// deepth is how many iterations the DichotomySearch() has been called.
// NOTE, eq must be monotonous on (x0, x1) or (x1, x0)! otherwise this function does not work!!
// NOTE, ensure eq(x0)*eq(x1)<0 stands at initial.
{
    real mid = (x0+x1)/2.0;
    real y_at_x0 =eq(x0);
    real y_at_x1 =eq(x1);
    real y_at_mid=eq(mid);
    if   (y_at_x0  == 0) return x0;
    elif (y_at_x1  == 0) return x1;
    elif (y_at_mid == 0) return mid;
    elif (fabs(x0-x1) <= 1e-10) return mid;
    elif (y_at_x0*y_at_mid < 0) DichotomySearch(eq, x0, mid, deepth+1);
    elif (y_at_x1*y_at_mid < 0) DichotomySearch(eq, x1, mid, deepth+1);
    elif (deepth>=50)
    {
        ErrorPrint("DetectZeroDenominator cannot find the root");
        return mid;
    }
    else
    { //        if (eq(x0)*eq(x1) > 0)
        ErrorPrint("You see this because your  eq(x0)*eq(x1)>0, or eq is not monotonous on the region between x0 and x1.");
        return mid;
    }
}/*}}}*/

real NewtonRaphsonMethod(real(*eq)(real), real(*eqDerivative)(real), real t0)
// find a root t_spike_delta (t_spike_delta = t_spike - t_n) (near t0) for eq /*{{{*/
// eq is the equation to solve,
// eqDerivative is eq's Derivative,
// x0 is the point from where we begin to search.
{
  real t_spike_delta = t0;
  for (int i=0; i<5; i++) t_spike_delta -= TimeCubicInterpolation(t_spike_delta)/TimeCubicInterpolationDerivative(t_spike_delta);
  return t_spike_delta;
}/*}}}*/

