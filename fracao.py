# -*- coding: utf-8 -*-
"""
Created on Sun May 22 19:19:35 2016

@author: Bossa
"""

from fractions import Fraction

class Fracao(Fraction):
    def __str__(self):
        return str(self.numerator)+'/'+str(self.denominator)

if __name__ == "__main__":
    A = Fracao("3/4")
    print(A)