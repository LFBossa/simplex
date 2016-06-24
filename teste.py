# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 23:21:23 2016

@author: Bossa
"""

import numpy as np
from scipy import linalg
AA = [[1, 2, 3],
      [0, 3, 4],
      [0, 0, 5]]
bb = [6, 7, 5]
A = np.matrix(AA)
b = np.array(bb)
sol = linalg.solve(A,b)
print(sol)

print(type(A),A.shape,type(b),b.shape)
b=b.reshape((3,1))
print(b)