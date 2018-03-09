#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  9 15:13:35 2018

@author: radioteddy
"""

import numpy as np
#initial parameters
X = 100
Y = 100
w = 60
d = 20
h = 40
l = int((Y-w)/2)
#initial distribution
U = np.zeros((X, Y))
U[h, l:l+w] = 100
U[h+d-1, l:l+w] = -100
np.savetxt('capacitor.dat', U, fmt='%.0f')
x = np.nonzero(U)[0]
y = np.nonzero(U)[1]
x = np.unique(x)
y = np.unique(y)
y1 = y[len(y)-1]
print(y1)