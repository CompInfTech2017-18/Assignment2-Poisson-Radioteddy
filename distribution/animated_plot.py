#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  3 01:20:16 2018

@author: radioteddy
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def evolution():
    U = np.loadtxt('/home/radioteddy/scripts/poisson/distribution/potential_'+str(i+1)+'.dat')
 