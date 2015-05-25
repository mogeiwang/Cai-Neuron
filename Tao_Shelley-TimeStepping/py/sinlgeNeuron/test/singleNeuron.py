# !/usr/bin/env python
# -*- coding:utf-8 -*-

# ==============================================#{{{
# ·
# · Author: Mogei Wang
# ·
# · MogeiWang@GMail.com
# ·
# · Filename: singleNeuron.py
# ·
# · COPYRIGHT 2015
# ·
# · Description:
# · This is an rewritten version of singleNeuron.c,
# ·      and the latter which has been debuged.T_T
# ==============================================#}}}


# --- system --- parameters {{{
# t_click_number_default = 1024 # these parameters will be set by test program
# RK_order_default = 4
# adjust_mode_default = True

t_end = 1.0
t_click_number = t_click_number_default # how many time clicks (number of time steps)
RK_order = RK_order_default # RK2 or RK4. both of them use the interpolation, so both errors should be of the order Dt^2.
adjust_mode = adjust_mode_default #True # use the interpolation or not.

t_step = t_end/t_click_number # compute time step according to the end time and number of time clicks (start time = 0)
t_list = [] # len[1+t_click_number] # from 0 to current time
v_list = [] # len[1+t_click_number] # this array is pretty large, usually do not use this large...
t_spike_list = [] #len[t_click_number] # some memory is wasted here :-)

g_leak = 50.0 # --- for single neurons
V_leak = 0.0 # the V parameters
V_E = 14/3.0
V_I = -2/3.0
V_threshold = 1.0
#}}}

# functions ----------------------------------------------------
def initialization():#{{{
    print("\nInitialization starts:\n")
    v_list.append(V_leak)
    t_list.append(0)
    return True

def deconstruction():
    print("\nDeconstruction finished.\n")
    return True

def RK2_stepper(df, x0, t0, ts):
    K1 = df(t0     , x0        )
    K2 = df(t0 + ts, x0 + ts*K1)
    return x0 + ts*( K1+K2 )/2.0

def RK4_stepper(df, x0, t0, ts):
    K1 = df(t0         , x0            )
    K2 = df(t0 + ts/2.0, x0 + ts*K1/2.0)
    K3 = df(t0 + ts/2.0, x0 + ts*K2/2.0)
    K4 = df(t0 + ts    , x0 + ts*K3    )
    return x0 + ts*( K1 + 2*K2 + 2*K3 + K4 )/ 6.0

def dynamics(t, v):
    return -g_leak * (v - V_leak) - 25.0 * sin(t) * (v - V_E) #- g_i(t) * (v - V_I)

def alpha(n):
    if (n==0): return g_leak + 25.0 * sin(t_list[-2])
    if (n==1): return g_leak + 25.0 * sin(t_list[-1])
    if (n==2): return g_leak + 25.0 * sin(t_list[-2] + 0.5*t_step) # this returns alphs(1/2) !!!
    else:
        print("ERROR parameter while calling beta")
        return False

def beta(n):
    if (n==0): return g_leak * V_leak + 25.0 * sin(t_list[-2]) * V_E
    if (n==1): return g_leak * V_leak + 25.0 * sin(t_list[-1]) * V_E
    if (n==2): return g_leak * V_leak + 25.0 * sin(t_list[-2] + 0.5*t_step) * V_E # this returns beta(1/2) !!!
    else:
        print("ERROR parameter while calling beta")
        return False

def adjust_t_spike_delta():
    # linera interpolation codes for RK2
    t_spike_delta = 0.0
    if ( LS(V_threshold, v_list[-2]) ): # all the float > < >= <= checked
        print("the last voltage > V_threshold!")
    if ( LS(v_list[-1], v_list[-2]) ):
        print("the current voltage < the last one!", v_list[-1], v_list[-2])
    t_spike_delta = t_step * (V_threshold - v_list[-2]) / (v_list[-1] - v_list[-2])
    if   (RK_order == 2):
        return t_spike_delta
    elif (RK_order == 4): # use NewtonRaphsonMethod() and search from t_ret; DichotomySearch() is also ok, but need to be more carefully
        return NewtonRaphsonMethod(TimeCubicInterpolation, TimeCubicInterpolationDerivative, t_spike_delta)
    else:
        print("only RK2 or RK4 are allowed!!")
        return False

def adjust_v_n():
    if (RK_order == 4):
        v_n = DichotomySearch(VoltageCubicInterpolation, V_leak-0.5, V_leak, 0) #Careful using DichotomySearch() to find v_n
#         v_n = float(v_n)
        K1 = beta(0) - alpha(0) * (v_n + 0.0              )
        K2 = beta(2) - alpha(2) * (v_n + 0.5 * K1 * t_step)
        K3 = beta(2) - alpha(2) * (v_n + 0.5 * K2 * t_step)
        K4 = beta(1) - alpha(1) * (v_n +       K3 * t_step)
        v_Delta = t_step*( K1 + 2*K2 + 2*K3 + K4 )/6.0
        v_nPlus1 = v_n + v_Delta
        return v_nPlus1
    elif (RK_order == 2):
        # the following are linera interpolation codes for RK2
        t_n = t_list[-2]
        t_spike_delta = t_spike_list[-1] - t_n
        v_n = (2 * V_leak - t_spike_delta * (  beta(0) +  beta(1) - alpha(1) *  beta(0) * t_step) ) \
                     / (2 - t_spike_delta * ( alpha(0) + alpha(1) - alpha(1) * alpha(0) * t_step) )
        K1 = beta(0) - alpha(0) * v_n
        K2 = beta(1) - alpha(1) * (v_n + K1 * t_step)
        v_nPlus1 = v_n + t_step*( K1 + K2 )/2.0
        return v_nPlus1
    else:
        print ("only RK2 or RK4 are allowed!!")
        return False

def VoltageCubicInterpolation(v_n):
# eq(15)=0 with left hand v(t) moved to right and set to V_leak
# this is used with DichotomySearch() to compute v_n,
# once v_n is got, v_nPlus1 can be obtained by eq(14)...
    t_n = t_list[-2]
    K1 = beta(0) - alpha(0) * (v_n + 0                )
    K2 = beta(2) - alpha(2) * (v_n + 0.5 * K1 * t_step)
    K3 = beta(2) - alpha(2) * (v_n + 0.5 * K2 * t_step)
    K4 = beta(1) - alpha(1) * (v_n +       K3 * t_step)
    v_Delta = t_step*( K1 + 2*K2 + 2*K3 + K4 )/6.0
    v_prime_n = beta(0) - alpha(0)*v_n
    v_prime_nPlus1 = beta(1) - alpha(1)*(v_n+v_Delta)
    t_spike_delta = t_spike_list[-1] - t_n
    # NOTE: according to eq(14), we have v_nPlus1-v_n = (K1+2*K2+2*K3+k4)*t_step/6
    # by substituting (v_nPlus1-v_n) into eq(15), then we get:
    # check all the linebreak problem
    return v_n + v_prime_n*t_spike_delta - V_leak \
        + ( 3*v_Delta - t_step * (2*v_prime_n + v_prime_nPlus1))*pow(t_spike_delta/t_step, 2) \
        + (-2*v_Delta + t_step * (  v_prime_n + v_prime_nPlus1))*pow(t_spike_delta/t_step, 3)

def TimeCubicInterpolation(t_spike_delta):
# eq(15) with left hand v(t) moved to right and set to V_threshold
# this is used with DichotomySearch() to adjust t_spike.
    v_n = v_list[-2]
    v_nPlus1 = v_list[-1]
    v_prime_n = beta(0) - alpha(0)*v_n
    v_prime_nPlus1 = beta(1) - alpha(1)*v_nPlus1
    t_n = t_list[-2]
    return v_n + v_prime_n*t_spike_delta - V_threshold \
        + ( 3*(v_nPlus1-v_n) - t_step * (2*v_prime_n + v_prime_nPlus1))*pow(t_spike_delta/t_step, 2) \
        + (-2*(v_nPlus1-v_n) + t_step * (  v_prime_n + v_prime_nPlus1))*pow(t_spike_delta/t_step, 3)

def TimeCubicInterpolationDerivative(t_spike_delta):
# Derivation of eq(15), i.e., eq(15)'t, (t is t_spike_delta here)
# this is used with NewtonRaphsonMethod() to adjust t_spike.
    v_n = v_list[-2]
    v_nPlus1 = v_list[-1]
    v_prime_n = beta(0) - alpha(0)*v_n
    v_prime_nPlus1 = beta(1) - alpha(1)*v_nPlus1
    t_n = t_list[-2]
    return v_prime_n \
        + ( 3*(v_nPlus1 - v_n) - t_step * (2*v_prime_n + v_prime_nPlus1))*2*t_spike_delta/pow(t_step, 2) \
        + (-2*(v_nPlus1 - v_n) + t_step * (  v_prime_n + v_prime_nPlus1))*3*pow(t_spike_delta, 2)/pow(t_step, 3)

def DichotomySearch(eq, x0, x1, deepth):
# find a root x for eq, x0<x<x1 OR x1<x<x0
# eq is the equation to solve,
# x0 and x1 are two ends of the region the root is in.
# deepth is how many iterations the DichotomySearch() has been called.
# NOTE, eq must be monotonous on (x0, x1) or (x1, x0)! otherwise this function does not work!!
# NOTE, ensure eq(x0)*eq(x1)<0 stands at initial.
    mid = (x0+x1)/2.0
    y_at_x0 =eq(x0)
    y_at_x1 =eq(x1)
    y_at_mid=eq(mid)
    if (deepth>=50):
        print ("DetectZeroDenominator cannot find the root")
        return mid
    if ( LE( fabs(x0-x1), 1e-10 ) ):
        return mid
    if   (y_at_x0  == 0):
        return x0
    elif (y_at_x1  == 0):
        return x1
    elif (y_at_mid == 0):
        return mid
    elif ( LS(y_at_x0*y_at_mid, 0) ):
        return DichotomySearch(eq, x0, mid, deepth+1)
    elif ( LS(y_at_x1*y_at_mid, 0) ):
        return DichotomySearch(eq, x1, mid, deepth+1)
    else:
        print("DichotomySearch() found eq(x0)*eq(x1)>0 at deepth: %5d, x0: %3.15f, x1: %3.15f,"%(deepth,x0,x1))
        print("\tIs the region between x1 and x0 right? Is the eq monotonous on the region?")
        return mid

def NewtonRaphsonMethod(eq, eqDerivative, t0):
# find a root t_spike_delta (t_spike_delta = t_spike - t_n) (near t0) for eq
# eq is the equation to solve,
# eqDerivative is eq's Derivative,
# x0 is the pofrom where we begin to search.
  t_spike_delta = t0
  for i in range(5):
      t_spike_delta -= TimeCubicInterpolation(t_spike_delta)/TimeCubicInterpolationDerivative(t_spike_delta)
  return t_spike_delta
#}}}

# ___main___ ---------------------------------------------------
initialization()
print("Runing")
for i in range(1, t_click_number+1): # main loop --- --- ... ...
    if   (RK_order==2): v_list.append( RK2_stepper(dynamics, v_list[-1], t_list[-1], t_step) )
    elif (RK_order==4): v_list.append( RK4_stepper(dynamics, v_list[-1], t_list[-1], t_step) )
    else: print("only RK2 or RK4 are allowed!!")
    t_list.append( t_list[-1] + t_step ) # the previous error was put this before v!!!
    if GE(v_list[-1], V_threshold):
        if adjust_mode:
            t_spike_delta = 0.0
            t_spike_delta = adjust_t_spike_delta()
            t_spike_list.append( t_list[-2] + t_spike_delta )
            v_list[-1] = adjust_v_n()
        else:
            t_spike_list.append( t_list[-1] )
            v_list[-1] -= (V_threshold - V_leak)

print("There are %ld spkikes, "%len(t_spike_list))
if (len(t_list) != len(t_list)):
    print("ERROR: length of time and v is not the same.")
else:
    print("   at the last click,", "TIME:",t_list[-1], "; VOLTAGE: ",v_list[-1])
deconstruction()
