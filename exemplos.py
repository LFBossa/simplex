# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 14:54:10 2016

@author: Bossa
""" 
from metodosimplex import SimplexPrimal, SimplexDual
# Problema com base completa
#A = [[ 1,  1,  2,  1,  0,  0],
#     [ 1,  1, -1,  0,  1,  0],
#     [ -1,  1,  1,  0,  0,  1]]
#b = [9, 2, 4] 

## Problema que não temos uma base completa
#A = [[ 1,  1, -1,  0,  0],
#     [-1,  1,  0, -1,  0],
#     [ 0,  1,  0,  0,  1]]
#b = [2, 1, 3]
#c = [1, -2]
 
## lista 2 problema 4
#A = [[1,2,1,0], 
#     [2,1,0,1]]
#b = [16,12]
#c = [-2,-5]


##lista 2 problema 7
#A = [[1,2,1,0,0],
#     [2,-1,0,1,0],
#     [5,3,0,0,1]]
#b = [6,4,15]
#c = [-5,-4]

##lista 2 problema 8
#A = [[3,4,0],
#     [2,-1,1]]
#b = [12,12]
#c = [1,-2]

##lista 2 problema 9
#A = [[-1, 1, 3,-3,1,0,0,0],
#     [-2, 2,-3, 3,0,1,0,0],
#     [ 2,-2, 1,-1,0,0,1,0],
#     [ 4,-4,-1, 1,0,0,0,1]]
#b = [3,6,8,16]
#c = [3,-3,-1,1] 

##apostila, exemplo regiao viavel vazia, pg 43
#A = [[1,1,1,0],
#     [2,3,0,-1]]
#b = [4,18]
#c = [-3,4] 

##Exemplo que criei de conjunto ilimitado
#A = [[-1,1,1,0],
#     [1,-2,0,1]]
#b = [4,4]
#c = [-1,1] 

# Exemplo redundância, apostila pg 44
#A = [[1, 1, 1, 0],
#    [-1, 1, 2, 0],
#     [0, 2, 3, 0],
#     [0, 0, 1, 1]]
#b = [6, 4, 10, 2]
#c = [-1,2,-3]



#problema = SimplexPrimal(A,b,c,output='minimal,file=bossa.txt')
#problema.resolver()

# Exemplo do Dual

#A = [[-1,0,-3,1,0],
#     [0,-1,-2,0,1]]
#b = [-3,-5]
#c = [4,6,18]
#
#problema = SimplexDual(A,b,c)
#problema.resolver()

# Lista 3
#A = [[-1,-2,-3,-1,1,0],
#     [-2, 1,-1, 3,0,1]]
#b = [-2,-3]
#c = [2,3,5,6]


# Exemplo Apostila, pg 73
A = [[-1,-2,-1,1,0],
     [-2, 1,-3,0,1]]
b = [-3,-4]
c = [2,3,4]
problema = SimplexDual(A,b,c)
problema.resolver() 