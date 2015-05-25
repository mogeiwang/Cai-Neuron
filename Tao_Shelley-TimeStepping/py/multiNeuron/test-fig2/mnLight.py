# !/usr/bin/env python
# -*- coding:utf-8 -*-

# ==============================================#{{{
# ·
# · Author: Mogei Wang
# ·
# · MogeiWang@GMail.com
# ·
# · Filename: multiNeurons.py
# ·
# · COPYRIGHT 2015
# ·
# · Description:
# · An implementation of
# · M J Shelley, L Tao's
# · Efficient and accurate time-stepping schemes for iaf neuronal networks
# · published in Journal of Computational Neuroscience, 11, 111, 2001
# ·
# · This is based on singleNeuron.py and singleNeuron.c.
# ·
# · Mogei Wang
# · 2015-5-6 --- 2015-5-21
# · in Shanghai
# ==============================================#}}}

# - imports ... -----------------------------------------------------#{{{
from __future__ import division
from __future__ import print_function
# from SimPy.Simulation import *
# from networkx import *
from pylab import *
from matplotlib.pyplot import * # static image
from matplotlib.animation import * # dynamic gif
from pygraphviz import *
from random import *
from math import *
from time import *
from decimal import *
import warnings
import pdb
import sys
import os
#
getcontext().prec = 100
warnings.simplefilter('ignore', DeprecationWarning)
# pdb.set_trace()
#}}}

# --- my mini function lib#{{{
def random_range(x):
    # return a random list, with each element being one of [0..x-1]
    a=[i for i in range(x)]
    shuffle(a)
    return a

def floatcmp(x, y):
    # returns -1, 0, 1, when x <, ==, > y
    try: return cmp(float(x), float(y))
    except ValueError: return cmp(x, y)

def EQ(x,y): # equal with
    if floatcmp(x,y) == 0: return True
    return False

def GT(x,y): # greater than
    if floatcmp(x,y) >  0: return True
    return False

def LS(x,y): # less than
    if floatcmp(x,y) <  0: return True
    return False

def GE(x,y): # greater or equal
    if floatcmp(x,y) >= 0: return True
    return False

def LE(x,y): # less or equal
    if floatcmp(x,y) <= 0: return True
    return False

def RK2_stepper(df, x0, t0, ts):
    # runge_kutta2 method.
    # this is used to move the system a step forward.
    # df is the dynamics function,
    # double x is input state,
    # double t0 is init time,
    # double ts is time step length,
    # returns the next state.
    K1 = df(t0     , x0        )
    K2 = df(t0 + ts, x0 + ts*K1)
    return x0 + ts*( K1+K2 )/2.0

def RK4_stepper(df, x0, t0, ts):
    # runge_kutta4 method.
    K1 = df(t0         , x0            )
    K2 = df(t0 + ts/2.0, x0 + ts*K1/2.0)
    K3 = df(t0 + ts/2.0, x0 + ts*K2/2.0)
    K4 = df(t0 + ts    , x0 + ts*K3    )
    return x0 + ts*( K1 + 2*K2 + 2*K3 + K4 )/ 6.0

def G_t_minus_tSpike(t=0.0, tau=0.0, m=5): # Eq.(3)
    # G(t−t_{spike}) describes synaptic induced conductance change,
    # both the excitatory and inhibitory postsynaptic conductance.
    # This function will be called by neuron class to compute G_e() and G_i().
    # Note that, when calling this function, t should minus a t_spike firstly.
    if LE(t, 0.0): return 0.0
    if LE(tau, 0.0): return 0.0
    return pow(t/tau, m) * exp(-t/tau)

def DichotomySearch(eq, x0, x1, deepth):
    # find a root x for eq, x0<x<x1 OR x1<x<x0
    # eq is the equation to solve,
    # x0 and x1 are two ends of the region the root is in.
    # deepth is how many iterations the DichotomySearch() has been called.
    # NOTE, eq must be monotonous on (x0, x1) or (x1, x0)! otherwise this function does not work!!
    # NOTE, ensure eq(x0)*eq(x1)<0 stands at initial.
    mid = (x0+x1)/2.0
    y_at_x0 = eq(x0)
    y_at_x1 = eq(x1)
    y_at_mid=eq(mid)
    if ( deepth >= 50 ): return mid
    if LE( fabs(x0-x1), 1e-10 ): return mid
    if   EQ( y_at_x0 , 0.0 ): return x0
    elif EQ( y_at_x1 , 0.0 ): return x1
    elif EQ( y_at_mid, 0.0 ): return mid
    elif ( LS(y_at_x0 * y_at_mid, 0.0) ): return DichotomySearch(eq, x0, mid, deepth+1)
    elif ( LS(y_at_x1 * y_at_mid, 0.0) ): return DichotomySearch(eq, x1, mid, deepth+1)
    else: return mid
#   ^^^ - print("DichotomySearch() found eq(x0)*eq(x1)>0 at deepth: %5d, x0: %3.15f, x1: %3.15f,"%(deepth,x0,x1))
#     | - print("\tIs the region between x1 and x0 right? Is the eq monotonous on the region?")

def NewtonRaphsonMethod(eq, eqDerivative, t0):
    # find a root t_spike_delta (t_spike_delta = t_spike - t_n) (near t0) for eq
    # eq is the equation to solve,
    # eqDerivative is eq's Derivative,
    # x0 is the pofrom where we begin to search.
    t_spike_delta = t0
    for i in range(5): # fixed to 5 iterations. The paper said, 2 or 3 iterations are ok.
        t_spike_delta -= eq(t_spike_delta)/eqDerivative(t_spike_delta)
    return t_spike_delta

def connection_distance(i, j, neuron_number_ring):
    # This function computes how long is the distance of two nodes on a same ring.
    # i and j are node numbers, neuron_number_ring are how many nodes on the ring.
    # the node numbers should be labeled as 0, 1, 2, ..., neuron_number_ring-1...
    if (i<0 or i>=neuron_number_ring or j<0 or j>=neuron_number_ring): return -1
    return min( abs(i-j),  neuron_number_ring + min(i,j) - max(i,j) )

def connection_intensity(i, j, neuron_number_ring, s=1.0, mu=0.0, delta=1.0):
    # This function computes the connection intensity between two nodes on a ring.
    # the intensity damps following the Gauss's low with an extra strength parameter
    # i and j are node numbers, neuron_number_ring are how many nodes on the ring.
    # mu:mean deviation (usually 0); delta:standard deviation; s:strength
    distance = connection_distance(i, j, neuron_number_ring)
    if LS(delta, 1e-5): return False
    if (distance == 0): intensity = 0.0 # avoid being divided by 0.
    else: intensity = s * exp( -0.5 * ( (distance - mu) / delta )**2 ) / sqrt( 2 * pi * delta )
    if (intensity > 1.0): intensity = 1.0
    return intensity
#}}}

# --- system --- parameters {{{
# t_click_number_default = 1000
# RK_order_default = 4
# adjust_mode_default = True

neuron_number = 5 # how many neurons are there in the system
t_end_default = 1.0 # will compute time step according to the end time and number of time clicks (start time = 0)

# potential: rest at 0 & threshold is 1
g_leak_default = 50.0 # --- for single neurons
V_leak_default = 0.0 # the V parameters
V_E_default = 14/3.0
V_I_default = -2/3.0
V_threshold_default = 1.0

m_default = 5 # m affects synaptic-induced conductance change, see Eq.(3) and G_t_minus_tSpike()
tau_e_default = 0.6e-3
tau_i_default = 1.0e-3

g_e0_intensity_default = 25.0
g_i0_intensity_default = 0.0
angle_field_default    = 2*pi
#}}}

class neuron_iaf_class:#{{{
    def __init__(self, ID=0):#{{{
        self.ID = ID # neuron ID
        self.adjust_mode = adjust_mode_default # use the interpolation or not.
        self.RK_order = RK_order_default # RK2 or RK4.
        self.t_click_number = t_click_number_default # how many time clicks (not used yet)
        self.t_end = t_end_default # end time, not yet used. maybe later in run_over()
        self.t_step = self.t_end/self.t_click_number # compute time step

        # the below parameters are shifted: x -> (x+70)/15
        # dv/dt =                           Eq.(4) in the paper
        #   -g_leak(v-V_leak)   -g_e(t)(v-V_E)   -g_i(t)(v-V_I)
        self.g_leak = g_leak_default
        self.V_leak = V_leak_default
        self.V_E = V_E_default
        self.V_I = V_I_default
        self.V_threshold = V_threshold_default

        self.m = m_default # m affects synaptic-induced conductance change @G_t_minus_tSpike()
        self.tau_e = tau_e_default # used in G_t_minus_tSpike
        self.tau_i = tau_i_default # used in G_t_minus_tSpike
        self.angle_field    = angle_field_default
        self.g_e0_intensity = g_e0_intensity_default
        self.g_i0_intensity = g_i0_intensity_default
        self.g_e0_angle_delta = self.angle_field * self.ID / neuron_number
        self.g_i0_angle_delta = self.angle_field * self.ID / neuron_number + pi
#         self.g_e0 = 0.0 # g_e1 driven by outter
#         self.g_i0 = 0.0 # g_i1 driven by outter

        self.g_e1 = 0.0 # 2nd term of g_e(t), i.e., neighbors' drive. see Eq.(5):
        self.g_i1 = 0.0 # they are temporary variable, but it should be kept...
        self.G_e  = 0.0 # 2nd factor of g_e. see Eq.(5):
        self.G_i  = 0.0 # 2nd factor of g_i. see Eq.(5):

        self.t_spike_list = [] #len would be t_click_number
        self.t_list = [0] # len would be 1+t_click_number, always starts from 0
        self.v_list = [self.V_leak] # len would be 1+t_click_number
        self.data_file_name = "/home/mw/neuronSci/data/neuronData%i.txt"%self.ID # the files
        self.data_file_handle = open(self.data_file_name,'w') # file handle, w+ or r+
#}}}

    def g_e0(self, t):#{{{
        # g_e0() first term of g_e(t), i.e., outter drive system. see Eq.(5):
        return self.g_e0_intensity * sin(t + self.g_e0_angle_delta)

    def g_i0(self, t):
        # g_i0() first term of g_i(t), i.e., outter drive system. see Eq.(5):
        return self.g_i0_intensity * sin(t + self.g_i0_angle_delta)

    def g_e(self, t):
        # g_e() computes Eq.(5), both outter drive and neighbor neurons' drive.
        # in g_e(), the outter drive g_e0() is computed directly,
        # and g_e1(), i.e., neighbor neurons' drive p4, is in fact computed in main(),
        # although g_e() is usually called by the neuron class
        # -- once dynamic is called, g_e() and g_i() will also be called...
        return self.g_e0(t) + self.g_e1

    def g_i(self, t):
        # g_i() computes Eq.(5), it is similar with g_e()
        return self.g_i0(t) + self.g_i1
#}}}

    def dynamics(self, t, v):#{{{
        # the membrane potential dynamics, i.e., Eq.(4)
        # dv/dt = -g_leak(v-V_leak)   -g_e(t)(v-V_E)   -g_i(t)(v-V_I)
        return -self.g_leak * (v - self.V_leak) - self.g_e(t) * (v - self.V_E) - self.g_i(t) * (v - self.V_I)

    # 1) neuron computes G_e and sends it to main(),  // captial G_e is a factor of g_e1
    # 2) main() computes g_e1 and sends it back to neuron,
    # 3) all drive information are stored in class, dynamics() can be called.
    # NOTE: the neuron dynamics() will call g_e() to compute all the drives...
#}}}


    def alpha(self, n=0): #{{{
        if EQ(n, 0): #alpha0:
            return self.g_leak + self.g_e(self.t_list[-2]) + self.g_i(self.t_list[-2])
        if EQ(n, 1): #alpha1:
            return self.g_leak + self.g_e(self.t_list[-1]) + self.g_i(self.t_list[-1])
        if EQ(n,0.5): #alpha(1/2)
            return self.g_leak + self.g_e(self.t_list[-2] + 0.5*self.t_step) + self.g_i(self.t_list[-2] + 0.5*self.t_step)
        else: return False

    def beta(self, n=0):
        if EQ(n, 0): #beta0:
                return self.g_leak * self.V_leak + self.g_e(self.t_list[-2]) * self.V_E + self.g_i(self.t_list[-2]) * self.V_I
        if EQ(n, 1): #beta1:
                return self.g_leak * self.V_leak + self.g_e(self.t_list[-1]) * self.V_E + self.g_i(self.t_list[-1]) * self.V_I
        if EQ(n,0.5): #beta(1/2)
                return self.g_leak * self.V_leak + self.g_e(self.t_list[-2] + 0.5*self.t_step) * self.V_E + self.g_i(self.t_list[-2] + 0.5*self.t_step) * self.V_I
        else: return False
#}}}

    def adjust_t_spike_delta(self):#{{{
        # linera interpolation codes for RK2
        t_spike_delta = 0.0
        if ( LS(self.V_threshold, self.v_list[-2]) ):
            print("the last voltage > V_threshold!")
        if ( LS(self.v_list[-1], self.v_list[-2]) ):
            print("the current voltage < the last one!", self.v_list[-1], self.v_list[-2])
        t_spike_delta = self.t_step * (self.V_threshold - self.v_list[-2]) / (self.v_list[-1] - self.v_list[-2])
        if   (self.RK_order == 2):
            return t_spike_delta
        elif (self.RK_order == 4): # use NewtonRaphsonMethod() and search from t_ret; DichotomySearch() is also ok, but need to be more carefully
            return NewtonRaphsonMethod(self.TimeCubicInterpolation, self.TimeCubicInterpolationDerivative, t_spike_delta)
        else:
            print("only RK2 or RK4 are allowed!!")
            return False

    def TimeCubicInterpolation(self, t_spike_delta):
        # eq(15) with left hand v(t) moved to right and set to V_threshold
        # this is used with DichotomySearch() to adjust t_spike.
        v_n = self.v_list[-2]
        v_nPlus1 = self.v_list[-1]
        v_prime_n = self.beta(0) - self.alpha(0)*v_n
        v_prime_nPlus1 = self.beta(1) - self.alpha(1)*v_nPlus1
        t_n = self.t_list[-2]
        return v_n + v_prime_n*t_spike_delta - self.V_threshold \
            + ( 3*(v_nPlus1-v_n) - self.t_step * (2*v_prime_n + v_prime_nPlus1))*pow(t_spike_delta/self.t_step, 2) \
            + (-2*(v_nPlus1-v_n) + self.t_step * (  v_prime_n + v_prime_nPlus1))*pow(t_spike_delta/self.t_step, 3)

    def TimeCubicInterpolationDerivative(self, t_spike_delta):
        # Derivation of eq(15), i.e., eq(15)'t, (t is t_spike_delta here)
        # this is used with NewtonRaphsonMethod() to adjust t_spike.
        v_n = self.v_list[-2]
        v_nPlus1 = self.v_list[-1]
        v_prime_n = self.beta(0) - self.alpha(0)*v_n
        v_prime_nPlus1 = self.beta(1) - self.alpha(1)*v_nPlus1
        t_n = self.t_list[-2]
        return v_prime_n \
            + ( 3*(v_nPlus1 - v_n) - self.t_step * (2*v_prime_n + v_prime_nPlus1))*2*t_spike_delta/pow(self.t_step, 2) \
            + (-2*(v_nPlus1 - v_n) + self.t_step * (  v_prime_n + v_prime_nPlus1))*3*pow(t_spike_delta, 2)/pow(self.t_step, 3)
#}}}

    def adjust_v_n(self): #{{{
        if (self.RK_order == 4):
            v_n = DichotomySearch(self.VoltageCubicInterpolation, self.V_leak-0.5, self.V_leak, 0) #Careful! using DichotomySearch() to find v_n
            K1 = self.beta(0)   - self.alpha(0)   * (v_n + 0.0                   )
            K2 = self.beta(0.5) - self.alpha(0.5) * (v_n + 0.5 * K1 * self.t_step)
            K3 = self.beta(0.5) - self.alpha(0.5) * (v_n + 0.5 * K2 * self.t_step)
            K4 = self.beta(1)   - self.alpha(1)   * (v_n +       K3 * self.t_step)
            v_Delta = self.t_step*( K1 + 2*K2 + 2*K3 + K4 )/6.0
            v_nPlus1 = v_n + v_Delta
            return v_nPlus1
        elif (self.RK_order == 2):
            # the following are linera interpolation codes for RK2
            t_n = self.t_list[-2]
            t_spike_delta = self.t_spike_list[-1] - t_n
            v_n = (2 * self.V_leak - t_spike_delta * (  self.beta(0) +  self.beta(1) - self.alpha(1) *  self.beta(0) * self.t_step) ) \
                            / (2 - t_spike_delta * ( self.alpha(0) + self.alpha(1) - self.alpha(1) * self.alpha(0) * self.t_step) )
            K1 = self.beta(0) - self.alpha(0) *  v_n
            K2 = self.beta(1) - self.alpha(1) * (v_n + K1 * self.t_step)
            v_nPlus1 = v_n + self.t_step*( K1 + K2 )/2.0
            return v_nPlus1
        else:
            print ("only RK2 or RK4 are allowed!!")
            return False

    def VoltageCubicInterpolation(self, v_n):
        # eq(15)=0 with left hand v(t) moved to right and set to V_leak
        # this is used with DichotomySearch() to compute v_n,
        # once v_n is got, v_nPlus1 can be obtained by eq(14)...
        t_n = self.t_list[-2]
        K1 = self.beta(0)   - self.alpha(0)   * (v_n + 0                     )
        K2 = self.beta(0.5) - self.alpha(0.5) * (v_n + 0.5 * K1 * self.t_step)
        K3 = self.beta(0.5) - self.alpha(0.5) * (v_n + 0.5 * K2 * self.t_step)
        K4 = self.beta(1)   - self.alpha(1)   * (v_n +       K3 * self.t_step)
        v_Delta = self.t_step*( K1 + 2*K2 + 2*K3 + K4 )/6.0
        v_prime_n      = self.beta(0) - self.alpha(0)* v_n
        v_prime_nPlus1 = self.beta(1) - self.alpha(1)*(v_n+v_Delta)
        t_spike_delta = self.t_spike_list[-1] - t_n
        # NOTE: according to eq(14), we have v_nPlus1-v_n = (K1+2*K2+2*K3+k4)*t_step/6
        # by substituting (v_nPlus1-v_n) into eq(15), then we get:
        return v_n + v_prime_n*t_spike_delta - self.V_leak \
            + ( 3*v_Delta - self.t_step * (2*v_prime_n + v_prime_nPlus1))*pow(t_spike_delta/self.t_step, 2) \
            + (-2*v_Delta + self.t_step * (  v_prime_n + v_prime_nPlus1))*pow(t_spike_delta/self.t_step, 3)
#}}}

    def run_step(self): #{{{
        if   (self.RK_order == 2): self.v_list.append( RK2_stepper( self.dynamics, self.v_list[-1], self.t_list[-1], self.t_step ) )
        elif (self.RK_order == 4): self.v_list.append( RK4_stepper( self.dynamics, self.v_list[-1], self.t_list[-1], self.t_step ) )
        else: return False
        self.t_list.append( self.t_list[-1] + self.t_step ) # the previous BUG was put this before updating v_list!!!

        if GE(self.v_list[-1], self.V_threshold) and LS(self.v_list[-2], self.V_threshold):
            if self.adjust_mode:
                self.t_spike_delta = self.adjust_t_spike_delta()
                self.t_spike_list.append( self.t_list[-2] + self.t_spike_delta )
                self.v_list[-1] = self.adjust_v_n()
            else:
                self.t_spike_list.append( self.t_list[-1] )
                self.v_list[-1] -= (self.V_threshold - self.V_leak)
        return (self.t_list[-1], self.v_list[-1])

        # computes g_e1 and g_i1, i.e., drive as a neighbors'
        if (len(self.t_spike_list) > 0):
            t = self.t_list[-1]
            self.G_e = 0.0
            self.G_i = 0.0
            for tt in self.t_spike_list:
                    self.G_e += G_t_minus_tSpike(t - tt, self.tau_e, self.m)
                    self.G_i += G_t_minus_tSpike(t - tt, self.tau_i, self.m)
#}}}
#}}}

# __main__ --------------------------------------------#{{{
# pdb.set_trace()

e_couple_matrix=[[connection_intensity(i, j, neuron_number)  for i in range(neuron_number)] for j in range(neuron_number)]
i_couple_matrix=[[connection_intensity(i, j, neuron_number)  for i in range(neuron_number)] for j in range(neuron_number)]

neuron_iaf = [] # instantiate all the iaf neurons, in order 0,1,..,n-1
for i in range(neuron_number): neuron_iaf.append(neuron_iaf_class(i))

current_G_e  = zeros(neuron_number) # get from neuron_iaf[i]
current_G_i  = zeros(neuron_number)
current_g_e1 = zeros(neuron_number) # compute in main(), and put to neuron_iaf[j]
current_g_i1 = zeros(neuron_number)

# system evolutions.
t_current = 0.0
t_step_default = t_end_default/t_click_number_default

for tttt in range(t_click_number_default):
#     print(t_current)
    for i in range(neuron_number): neuron_iaf[i].run_step() # NOT in order
    t_current += t_step_default

    for i in range(neuron_number): current_G_e[i] = neuron_iaf[i].G_e
    for i in range(neuron_number): current_G_i[i] = neuron_iaf[i].G_i

    # calculate... matrix*vector
    current_g_e1 = dot(e_couple_matrix, current_G_e)
    current_g_i1 = dot(i_couple_matrix, current_G_i)

    # evolution...
    for i in range(neuron_number): neuron_iaf[i].g_e1 = current_g_e1[i]
    for i in range(neuron_number): neuron_iaf[i].g_i1 = current_g_i1[i]
#}}}

# figure()#{{{
# subplot(4,1,1)
# plot(neuron_iaf[0].t_list, neuron_iaf[0].v_list, label='$v^1$')
# xlim([0, t_end_default])
# ylabel("voltage of neuron 1")
# subplot(4,1,2)
# plot(neuron_iaf[1].t_list, neuron_iaf[1].v_list, label='$v^2$')
# xlim([0, t_end_default])
# ylabel("voltage of neuron 2")
# subplot(4,1,3)
# plot(neuron_iaf[2].t_list, neuron_iaf[2].v_list, label='$v^3$')
# xlim([0, t_end_default])
# ylabel("voltage of neuron 3")
# subplot(4,1,4)
# plot(neuron_iaf[3].t_list, neuron_iaf[3].v_list, label='$v^4$')
# xlim([0, t_end_default])
# ylabel("voltage of neuron 4")
# xlabel("time ($S$)")
# # title("")
# savefig('coupled_Neurons.pdf', dpi=600)

print("The last click:")
for i in range(neuron_number): print("\t T%d:"%i,neuron_iaf[i].t_list[-1], "VOL%d:"%i, neuron_iaf[i].v_list[-1])#}}}
