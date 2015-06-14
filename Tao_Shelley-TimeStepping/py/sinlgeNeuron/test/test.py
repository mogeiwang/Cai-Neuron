# !/usr/bin/env python
# -*- coding:utf-8 -*-


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

getcontext().prec = 100

# --- my mini function lib
def random_range(x):
    # return a random list, with each element being one of [0..x-1]
    a=[i for i in range(x)]
    shuffle(a)
    return a

def floatcmp(x, y):
    # returns -1, 0, 1, when x <, ==, > y
    try:
        return cmp(float(x), float(y))
    except ValueError:
        return cmp(x, y)

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
#}}}

t_click_number_default = 1024000#{{{
RK_order_default = 4
adjust_mode_default = True
print("canonical")
execfile("singleNeuron.py")#}}}


t_click_number_default = 1000#{{{
RK_order_default = 2
adjust_mode_default = False
execfile("singleNeuron.py")
t_click_number_default = 2000
RK_order_default = 2
adjust_mode_default = False
execfile("singleNeuron.py")
t_click_number_default = 4000
RK_order_default = 2
adjust_mode_default = False
execfile("singleNeuron.py")
t_click_number_default = 8000
RK_order_default = 2
adjust_mode_default = False
execfile("singleNeuron.py")
t_click_number_default = 16000
RK_order_default = 2
adjust_mode_default = False
execfile("singleNeuron.py")
t_click_number_default = 32000
RK_order_default = 2
adjust_mode_default = False
execfile("singleNeuron.py")
t_click_number_default = 64000
RK_order_default = 2
adjust_mode_default = False
execfile("singleNeuron.py")
t_click_number_default = 128000
RK_order_default = 2
adjust_mode_default = False
execfile("singleNeuron.py")#}}}

t_click_number_default = 1000#{{{
RK_order_default = 2
adjust_mode_default = True
execfile("singleNeuron.py")
t_click_number_default = 2000
RK_order_default = 2
adjust_mode_default = True
execfile("singleNeuron.py")
t_click_number_default = 4000
RK_order_default = 2
adjust_mode_default = True
execfile("singleNeuron.py")
t_click_number_default = 8000
RK_order_default = 2
adjust_mode_default = True
execfile("singleNeuron.py")
t_click_number_default = 16000
RK_order_default = 2
adjust_mode_default = True
execfile("singleNeuron.py")
t_click_number_default = 32000
RK_order_default = 2
adjust_mode_default = True
execfile("singleNeuron.py")
t_click_number_default = 64000
RK_order_default = 2
adjust_mode_default = True
execfile("singleNeuron.py")
t_click_number_default = 128000
RK_order_default = 2
adjust_mode_default = True
execfile("singleNeuron.py")#}}}


t_click_number_default = 1000#{{{
RK_order_default = 4
adjust_mode_default = False
execfile("singleNeuron.py")
t_click_number_default = 2000
RK_order_default = 4
adjust_mode_default = False
execfile("singleNeuron.py")
t_click_number_default = 4000
RK_order_default = 4
adjust_mode_default = False
execfile("singleNeuron.py")
t_click_number_default = 8000
RK_order_default = 4
adjust_mode_default = False
execfile("singleNeuron.py")
t_click_number_default = 16000
RK_order_default = 4
adjust_mode_default = False
execfile("singleNeuron.py")
t_click_number_default = 32000
RK_order_default = 4
adjust_mode_default = False
execfile("singleNeuron.py")
t_click_number_default = 64000
RK_order_default = 4
adjust_mode_default = False
execfile("singleNeuron.py")
t_click_number_default = 128000
RK_order_default = 4
adjust_mode_default = False
execfile("singleNeuron.py")#}}}

t_click_number_default = 1000#{{{
RK_order_default = 4
adjust_mode_default = True
execfile("singleNeuron.py")
t_click_number_default = 2000
RK_order_default = 4
adjust_mode_default = True
execfile("singleNeuron.py")
t_click_number_default = 4000
RK_order_default = 4
adjust_mode_default = True
execfile("singleNeuron.py")
t_click_number_default = 8000
RK_order_default = 4
adjust_mode_default = True
execfile("singleNeuron.py")
t_click_number_default = 16000
RK_order_default = 4
adjust_mode_default = True
execfile("singleNeuron.py")
t_click_number_default = 32000
RK_order_default = 4
adjust_mode_default = True
execfile("singleNeuron.py")
t_click_number_default = 64000
RK_order_default = 4
adjust_mode_default = True
execfile("singleNeuron.py")
t_click_number_default = 128000
RK_order_default = 4
adjust_mode_default = True
execfile("singleNeuron.py")#}}}


