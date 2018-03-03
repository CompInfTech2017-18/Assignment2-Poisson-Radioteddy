# -*- coding: utf-8 -*-
"""
Created on Fri Mar  2 21:24:18 2018

@author: Элвис
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
#import os

#dir_path = os.path.abspath('Poisson.py')
#print(dir_path)

class solution:
    def __init__(self, Ox, Oy, eps, cut): #cut is fourier clipping
        self.err = eps
        self.X = Ox
        self.Y = Oy
        self.init = np.zeros((Ox, Oy))
        self.cut = cut        
        
    def Fourier(self):
        x = np.linspace(0, self.X, self.X)
        y = np.linspace(0, self.Y, self.Y)
        xgrid, ygrid = np.meshgrid(x, y)
        x = xgrid
        y = ygrid
        a = 100
        U0 = 100
        U = 0
        n = 0
        while n <= self.cut:
            U += 4*U0/(np.pi*(2*n+1))*( np.cosh(x*(2*n+1)*np.pi/a)-1/np.tanh((2*n+1)*np.pi)*np.sinh(x*(2*n+1)*np.pi/a) )* np.sin(y*(2*n+1)*np.pi/a) 
            n += 1
        return U
        
    def Jacobi(self):
        self.U = self.init
        self.U[:, 0] = 100
        self.U[self.X-1, 1:] = 0
        self.U[0, 1:] = 0
        self.U[:, self.Y-1] = 0       
        for i in range(5000):
            self.init = self.U
            self.U[1:self.X-1, 1:self.Y-1] = (self.init[0:-2, 1:-1]+self.init[2:, 1:-1]+self.init[1:-1, 2:]+self.init[1:-1, 0:-2])/4
#            np.savetxt('/home/radioteddy/scripts/poisson/distribution/potential_'+str(i+1)+'.dat', self.U, fmt='%.5f') #intermediate distribution
        return self.U
    
    def Gauss_Seidel(self):
        U_init = 0
        U_final = 1
        self.U = self.init
        self.U[:, 0] = 100
        self.U[self.X-1, 1:] = 0
        self.U[0, 1:] = 0
        self.U[:, self.Y-1] = 0
        while (np.abs(U_init - U_final) > self.err):
            U_init = np.abs(np.trace(self.U))
            self.U[1:self.X-1, 1:self.Y-1] = (self.U[0:-2, 1:-1]+self.U[2:, 1:-1]+self.U[1:-1, 2:]+self.U[1:-1, 0:-2])/4 
            U_final = np.abs(np.trace(self.U))
        return self.U
            
    def plotter(self, U, method):
        x = np.linspace(0, self.X, self.X)
        y = np.linspace(0, self.Y, self.Y)
        self.xgrid, self.ygrid = np.meshgrid(x, y)
        
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        face = ax.plot_surface(self.xgrid, self.ygrid, U, cmap=cm.plasma, linewidth=0, antialiased=False)        
        ax.set_zlim(0, 100)
        ax.zaxis.set_major_locator(LinearLocator(6))
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        fig.colorbar(face, shrink=0.5, aspect=10)
        ax.set_title(method)
        plt.show()
        
    def run(self, method):
        if method == 'Jacobi':
            self.plotter(self.Jacobi(), 'Jacobi')
        elif method == 'Fourier':
            self.plotter(self.Fourier(), 'Fourier')
        elif method == 'Gauss-Seidel':
            self.plotter(self.Gauss_Seidel(), 'Gauss-Seidel')
    
test = solution(100, 100, 0.001, 100)
test.run('Fourier')
test.run('Jacobi')
test.run('Gauss-Seidel')

        
        