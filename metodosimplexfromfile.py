import numpy as np
np.set_printoptions(precision=3)

class Simplex(object): 
    ''' Classe que contém todo o procedimento de implementação do método simplex.
    @author: lfbossa@gmail.com
    @date: 2016
    '''
    def __init__(self,A,b,c):
        self.A = np.matrix(A)
        self.b = np.matrix(b)
        self.c = np.matrix(c)  
    
    def set_base(self,base):
        '''Define quem são os caras da base.'''
        self.base = base
        
    def solve(self):
        '''Supondo que temos uma identidade nas últimas m colunas, 
        e que 0 é solução inicial, executa o simplex.'''
        (m,n) = self.A.shape
        self.tableau = np.zeros((m+1,n+1))
        self.tableau[0:m,0:n] = self.A     # Colocamos A, b e c no tableau 
        self.tableau[0:m,n] = self.b       # (lembrando sempre que os indices começam em 0)
        self.tableau[m,0:self.c.shape[1]] = -self.c
        self.base = [n-m+i for i in range(m)]
        if self.show:
            print(self.tableau)
        self.run()
        
    def fase1(self):
        '''Resolve o problema auxiliar min 1^Tx_a para encontrar uma primeira 
        solução viável'''
        (m,n) = self.A.shape
        # cria um tableau de zeros
        self.tableau = np.zeros((m+1,n+m+1))
        self.tableau[0:m,0:n] = self.A       # Colocamos A, b e c no tableau 
        self.tableau[0:m,n+m] = self.b       # (lembrando sempre que os indices começam em 0)
        self.tableau[0:m,n:n+m] = np.eye(m)  # adicionamos as variáveis artificiais 
        self.base = [i for i in range(n,m+n)]
        if self.show:
            print('iniciando fase 1 \n\n')
        (m,n) = self.A.shape
        self.tableau[m,n:m+n] = -np.ones((1,m))
        if self.show:
            print(self.tableau)
        print("pivotando nas variáveis artificiais")
        for i in range(m):
            self.tableau[m,:] = self.tableau[m,:] + self.tableau[i,:]
        print(self.tableau)
        self.run()
        deletar = []
        artificialnabase = False
        # começamos supondo que nenhuma variável aritificial está na base
        for k in range(n,m+n):
            if k in self.base:
                artificialnabase = True
                # se alguma das variáveis artificiais pertencer a base
                i = self.tableau[0:m,k].argmax()
                if self.tableau[i,-1] != 0:
                    # e seu valor for diferente de zero, 
                    # então não temos soluções viáveis
                    raise SemSolucoesViaveis()
                else:
                    # caso contrário, temos redundância e vamos eliminar 
                    # a linha a qual adicionamos a variável artificial
                    deletar.append(i-n)
        if artificialnabase:
            # deletar as linhas 
            self.tableau = np.delete(self.tableau, deletar,axis=0)
        pass     
    
    def fase2(self):
        '''Método para rodar a fase 2 do simplex.'''
        if self.show:
            print('iniciando fase 2 \n\n')
        (m,n) = self.A.shape
        deletar = [i for i in range(n,n+m)]
        self.tableau = np.delete(self.tableau, deletar,axis=1)
        self.tableau[m,0:self.c.shape[1]] = -self.c
        if self.show:
            print('vetor de custo relativo adicionado.')
            print(self.tableau)
        for j in self.base:
            i = self.tableau[0:m,j].argmax()
            self.tableau[m,:] = self.tableau[m,:] - self.tableau[m,j]*self.tableau[i,:]
        if self.show:
            print('custo relativo dos caras da base zerados.')
            print(self.tableau)
        self.run()
        pass
    
    def solucao(self):
        '''Retorna o vetor solução do problema.'''
        (m,n) = self.A.shape
        x = np.zeros((m+n,1))
        for i in range(m):
            x[self.base[i],0] = self.tableau[i,-1]
        return x
    
    def run(self):
        '''Roda a parte do pivotamento usado em casa fase do simplex.'''
        (m,n) = self.tableau.shape 
        status = 0
        maximo = self.tableau[m-1,0:n-1].max()
        while  maximo > 0:
            quemEntra = self.tableau[m-1,0:n-1].argmax()   # indice do maior elemento da ultima linha do tableau
            quemSai = self.quem_sai_da_base(quemEntra)     
            self.base[quemSai] = quemEntra                # atualiza a base
            if self.show:
                print("pivotando na linha {} e coluna {}".format(quemSai+1, quemEntra +1))
            self.pivotar(quemSai,quemEntra)               # pivota na linha e coluna indicadas
            maximo = self.tableau[m-1,0:n-1].max()        # procura um novo cara pra entrar na base 
        pass
    
    def quem_sai_da_base(self,k):        
        '''Dado que x_k entra na base, essa função determina quem sai da base.'''
        (m,n) = self.tableau.shape
        j = -1 
        I = [i for i in range(m-1) if self.tableau[i,k] > 0] 
        if len(I) == 0:
            # se y_k <= 0, o problema é ilimitado
            raise ProblemaIlimitado()
        # I-1 é o conjunto dos indices i tais que yik > 0 
        razoes = np.array([self.tableau[i,j]/self.tableau[i,k] for i in I]) 
        # calcula as razoes b_i/y_ik para y_ik > 0
        razaominima = razoes.min()         
        # encontra o valor da razao mínima
        I = [I[i] for i in range(razoes.size) if razoes[i] == razaominima]  
        # conjunto I0 dos indices r tais que br/yrk é igual a razao minima
        # se I0  tiver mais de  1 elemento, usamos a validação lexicográfica
        while len(I) > 1:
            j = j + 1 # calcula os conjuntos I_j, para j <= n 
            razoes = np.array([self.tableau[i,j]/self.tableau[i,k] for i in I])
            razaominima = razoes.min()
            I = [I[i] for i in range(razoes.size) if razoes[i] == razaominima] 
            print(razoes, I)
            if j == n-1:
                # caso I_j contenha mais do que 1 elemento para todo j de 1 a n, 
                # temos que o tableau tem pelo menos duas linhas l.d. 
                raise LinhasLD(I) 
        return I[0] 
    
    def pivotar(self,i,j): 
        '''Pivota o tableau na linha i e coluna j.'''
        (m,n) = self.tableau.shape 
        self.tableau[i,:] = self.tableau[i,:]/self.tableau[i,j] # Divide a linha i pelo pivo 
        for k in range(m):
            if k != i:
                self.tableau[k,:] = self.tableau[k,:] - self.tableau[k,j]*self.tableau[i,:]  
                # Faz as eliminações nas outras linhas
        if self.show == True:
            print(self.tableau,'\n')
        pass
    def show(self,booleano):
        '''Método que diz se iremos mostrar o tableau a cada iteração ou não.'''
        self.show = booleano
        pass
    
    
    
# Fim da Classe Simplex 


class LinhasLD(Exception):
    def __init__(self, linhas):
        self.linhas = linhas
        self.texto = ''
        self.txt()
    def txt(self):
        K = len(self.linhas)
        lnhs =  ''
        for i in range(K):
            if i < K-1:
                lnhs =  lnhs + '{}, '
            else:
                lnhs = lnhs + 'e {} '
        self.texto = 'As linhas ' + lnhs + 'sao linearmente dependentes.'
        self.texto = self.texto.format(*self.linhas)
    def __str__(self):
        return self.texto


class SemSolucoesViaveis(Exception):
    def __str__(self):
        return 'O conjunto viável é vazio.'

    
class ProblemaIlimitado(Exception):
    def __str__(self):
        return 'O -gradiente aponta em uma direção ilimitada do conjunto.'    
        
    