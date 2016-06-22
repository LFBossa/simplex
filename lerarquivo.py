import fractions as fr
import numpy as np
from metodosimplex import Simplex


class Fracao(fr.Fraction):
    def __str__(self):
        if self.denominator == 1:
            txt = self.numerator
        else:
            txt = str(self.numerator)+'/'+str(self.denominator)
        return txt

class Matriz(np.matrix):
    def __srt__(self):
        (m,n) = self.shape
        txt = ''
        for i in range(m):
            
        
# le um arquivo e retorna as matrizes A, b e c 
# do problema de otimzacão min c*x s.a. Ax <= b
class SimplexFromFile(object):
    def __init__(self,filename):
        self.A = []
        self.b = []
        self.c = []
        sinais  = []
        l = 0
        with open(filename, 'r') as arquivo: 
            for linha in arquivo: 
                if linha[-1]=='\n': #caso a quebra de linha seja colada com o caractere
                    linha = linha[:-1]
                    #separa a linha por espaços, e ignora mais do que 1 espaço seguido
                linhaprepronta = [block for block in linha.split(' ') if block!='']
                if l==0:
                    self.c = [Fracao(fracao) for fracao in linhaprepronta]
                else:
                    linhafrac = [Fracao(fracao) for fracao in linhaprepronta[0:-2]]
                    self.b.append(Fracao(linhaprepronta[-1]))
                    self.A.append(linhafrac)
                    sinais.append(linhaprepronta[-2])
                l+=1
        ## Vamos adicionar as variáveis de folga/excesso        
        m = len(self.A)
        n = len(self.A[0])
        VarFolga  = [[Fracao(0) for x in range(m)] for y in range(m)]
        VarArt = [[] for y in range(m)]
        # retorna uma matriz de zeros, mxm
        self.base = [False for i in range(m)] #um vetor com os caras que estarao na base
        for k in range(m):
            if sinais[k]=='<=': # variável de folga
                VarFolga[k][k] =  Fracao(1)
                self.base[k] = n+k
            elif sinais[k]=='>=': # variável de excesso
                VarFolga[k][k] = Fracao(-1)
            else:                # gambiarra 
                VarFolga[k][k] = Fracao(1)
        for k in range(m):  
            if self.base[k]==False:  # remove as colunas de zeros
                for i in range(m):
                    VarArt[i].extend([VarFolga[i][k]])
                    del(VarFolga[i][k])
        # colando as varáveis de folga na matriz A
        for k in range(m):
            self.A[k].extend(VarFolga[k]) 
        if len(VarFolga[0]) == m:
            self.tembase = True
        else:
            self.tembase = False
    def problema(self):
        return Simplex(self.A,self.b,self.c)


fromfile = SimplexFromFile('matrix.txt')
print(fromfile.A)
#problema = fromfile.problema()
#problema.show(True)
#problema.fase1()
#problema.fase2()
            
    


    