# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 14:11:29 2016

@author: Bossa

Algoritmo 6.3.2: Primal-Dual Path Following

Powered by Spider 2

"""
import numpy as np 

np.set_printoptions(precision=3)

class PrimalDualPF(object):
    '''
    Implementação do método Primal Dual Path Following para o problema 
    de programação linear 
        min  c*x
        s.a. Ax = b
             x >= 0   
    
    Parâmetros
    ----------
    A: matriz mxn com as variáveis de folga já inseridas.
    b: vetor nx1 
    c: vetor de custo mx1
    x: solução inicial do problema primal 
    w: solução inicial do problema dual
    s: folga dual associada a solução w
    '''
    def __init__(self, A, b, c, x, w, s):
        self.A = np.array(A) 
        # vetores serão arrays deitados
        self.b = np.array(b)
        # completando o vetor c com zeros
        self.c = np.zeros_like(A[0])
        self.c[0:len(c)] = np.array(c)
        self.x = np.array(x)
        self.w = np.array(w)
        self.s = np.array(s) 
        (m,n) = self.A.shape
        self.mu = np.dot(self.x,self.s)/n
        self.gamma = 1e-2
        
    def update(self):
        (m,n) = self.A.shape
        self.mu = np.dot(self.x,self.s)/n
        #Calculando a matriz Jacobiana
        J1 = np.concatenate( 
        (np.zeros((n,n)), np.transpose(self.A), np.eye(n)),axis=1)
        J2 = np.concatenate(
        (self.A, np.zeros((m,m)), np.zeros((m,n))),axis=1)
        J3 = np.concatenate(
        (np.diag(self.s), np.zeros((n,m)), np.diag(self.x)),axis=1)
        self.jacobian = np.concatenate((J1,J2,J3))
        # Calculando o lado direito, que será um array deitado
        D1 = self.c - self.A.T.dot(self.w) - self.s
        D2 = self.b - self.A.dot(self.x)
        D3 = -1*np.multiply(self.s,self.x)*np.ones(n) + \
               self.sigma*self.mu*np.ones(n)
        self.RHS = np.concatenate((D1,D2,D3))
#        RHS.write(str(np.asarray(self.RHS)))
#        RHS.write('\n')
#        LHS.write(str(np.asarray(self.jacobian)))
#        LHS.write('\n')
        
    def step(self):
        (m,n) = self.A.shape 
        solucao = np.linalg.solve(self.jacobian,self.RHS) 
        # quebramsos a solucao do sistema linear em 
        # solucao = [px,pw,ps]
        px = solucao[0:n]
        pw = solucao[n:m+n]
        ps = solucao[m+n:m+2*n]
        alpha = 1.0
        backtracking = True
        while backtracking:
            # sejam 
            # (x_k+1,w_k+1,s_k+1) = (x_k,w_k,s_k) +alpha(px_k,pw_k,ps_k)
            novox = self.x + alpha*px
            novow = self.w + alpha*pw
            novos = self.s + alpha*ps 
            const = self.gamma*self.mu
            if validacao(novox,const) and validacao(novos,const):
                # se (x_k+1,s_k+1) > gamma*mu, tudo bem
                backtracking = False
            else:
                # senão, diminuímos o tamanho do passo
                alpha = 0.9*alpha
        self.x = novox
        self.w = novow
        self.s = novos
        
    def solve(self, TOL=1e-3,MAXSTEPS=1000, SIGMA=.2):
        '''Itera os passos de resolução até que a precisão TOL seja alcançada'''
        self.sigma = SIGMA
        self.caminho = [[],[]]
        for k in range(MAXSTEPS):
            if self.mu < TOL:
                break
            else:
                for k in range(2):
                    self.caminho[k].append(self.x[k])
                self.update()
                self.step()
            
    def solution(self):
        return (self.x,self.w,self.s)
            
def validacao(vector,constante):
    '''Verifica se cada entrada de um dado vetor é maior que a constante.'''
    test = True
    for x in vector:
        if x <= constante:
            test = False
            break
    return test
    


A = [[ 1, -2, -1,  0],
     [-1, -1,  0, -1]]
b = [-4,-4]
c = [-1,-3]
x = [1, 1, 3, 2]
s = [1, 2, 1, 3]
w = [1, 3]

#for k in range(1,10):
#    problema = PrimalDualPF(A,b,c,x,w,s)
#    problema.solve(SIGMA=1/k,MAXSTEPS=200)
#    print('sigma:', 1/k,';','mu:', problema.mu) 
import matplotlib.pyplot as plt
 
TOTAL_ITERATIONS = 25
for i in range(1,TOTAL_ITERATIONS):
    sigma = i/TOTAL_ITERATIONS
    problema = PrimalDualPF(A,b,c,x,w,s) 
    problema.solve(SIGMA=sigma)
    plt.plot(problema.caminho[0],problema.caminho[1],'b',alpha=1-sigma,
             hold=True)

plt.show() 
