# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 14:11:29 2016

@author: Bossa

Algoritmo 6.3.2: Primal-Dual Path Following

"""
import numpy as np
from scipy import linalg

class PrimalDualPF(object):
    '''Implementação do método Primal Dual Path Following'''
    def __init__(self, A, b, c, s, x, w):
        '''Inicia as variáveis.'''
        self.A = np.matrix(A)
        (m,n) = self.A.shape
        # Vamos deixar todos os vetores em pé
        self.b = np.matrix(b).reshape([n,1])
        self.c = np.matrix(c).reshape([m,1])
        self.s = np.matrix(s).reshape([n,1])
        self.x = np.matrix(x).reshape([n,1])
        self.w = np.matrix(w).reshape(-1).transpose()
        # constantes
        self.sigma = 0.9
        self.mu = np.vdot(self.x,self.s)/n
        self.update()
        
    def update(self):
        (m,n) = self.A.shape
        #Calculando a matriz Jacobiana
        J1 = np.concatenate( 
        (np.zeros([n,n]), self.A.transpose(), np.eye(n)),axis=1)
        J2 = np.concatenate(
        (self.A, np.zeros([m,m]), np.zeros([m,n])),axis=1)
        J3 = np.concatenate(
        (np.diag(self.s), np.zeros([n,m]), np.diag(self.x)),axis=1)
        self.jacobian = np.concatenate((J1,J2,J3))
        # Calculando o lado direito
        self.RHS = np.concatenate((
        self.c - self.A.transpose()*self.w - self.s,
        self.b - self.A*self.x,
        np.diag(self.s)*np.diag(self.x)*np.ones([n,1]) 
                        -self.sigma*self.mu*np.ones([n,1])
        ))
    def step(self):
        newsol = linalg.solve(self.jacobian,self.RHS)
        alpha = 
        
        

A = [[ 1, -2, -1,  0],
     [-1, -1,  0, -1]]
b = [-4,-4]
c = [-1,-3]

