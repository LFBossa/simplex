# -*- coding: utf-8 -*-
import numpy as np
np.set_printoptions(precision=3)

class SimplexPrimal():
    '''
SimplexPrimal
=============
Classe que contém todo a implementação do método simplex primal para resolver
o problema de programação linear na forma padrão
::
    min  <c,x>
    s.a. Ax  = b
          x >= 0

Parâmetros
----------
A : m x n array
    Matriz das restrições lineares, com as variáveis de folga já inclusas.
b : m array
    Lado direito das restrições lineares, que deve necessariamente ser
    não-negativo.
c : p array
    Gradiente da função objetivo. Se p <= n, o algoritmo completa c com zeros.
output : string
    String contendo as configurações de onde será feito o output do problema.
    Aceita as seguintes opções, separadas por vírgula.

    * `file=filename.txt`: imprime o output no arquivo filename.txt

    * `latex=filename.tex`: imprime o output em formato LaTeX no arquivo filename.tex

    * `screen`: imprime o output na tela

    * `minimal`: reduz o tamanho do output, não imprimindo todos os passos do pivotamento



Autor
-----

lfbossa@gmail.com, 2016
    '''
    def __init__(self,A,b,c, output='screen,minimal'):
        self.A = np.array(A)
        self.b = np.array(b)
        self.c = np.array(c)
        self.output = output
        self.check()
        # checa se b >= 0

    def check(self):
        zeros = np.zeros_like(self.b)
        posicoes, = np.where(self.b < zeros)
        # retorna as posições em que b < 0
        if len(posicoes) > 0: # Se alguma entrada de b for negativa
            raise ValueError("O vetor b deve ser maior que zero.")


    def resolver(self):
        '''Analisa a matriz A e determina se temos ou não uma base completa
        nas últimas m colunas.
        Se sim, executa o método `jatembase`.
        Se não, executa a `fase1` e `fase2`.
        '''
        self.preparar_output('abrir')
        self._output('','definirproblema')
        (m,n) = self.A.shape
        # Para cada vetor da base canônica, vamos verificar se e onde
        # ele se encontra na matriz A
        self.base = m*[False]
        # Iniciamos com uma base "vazia"
        for j in range(m):
            ej = basecanonica(j,m)
            ej = ej.reshape((m,1)) # deixamos o vetor em pé
            vetortruefalse = np.all(self.A==ej,axis=0)
            # vetor 1xn tal que vetortruefalse[j] = True se e somente se
            # A[:,j] = basecanonica[j,m]
            posicao, = np.where(vetortruefalse)
            # posicao é um vetor que contém os indices das colunas para
            # as quais vetortruefalse é True
            if len(posicao) > 0:
                # i.e., se existir alguma coluna da identidade em A
                coluna = posicao[0] # pegamos o primeiro indice
                self.base[j] = coluna
        # Se não tivermos uma base completa, rodamos a `fase1` e `fase2`
        if False in self.base:
            msg = "Problema não tem base, vamos para a fase 1."
            self._output(msg,'[m]resolver')
            try:
                self.fase1()
                self.fase2()
                self._output("Solução x =",'[m]solucao')
                self._output(self.solucao(),'[m]solucao')
            except Exception as e:
                self._output(str(e),'[m]error')
                self.preparar_output('fechar')
        else:
            msg = "Problema tem base, vamos começar a pivotar."
            self._output(msg,'[m]resolver')
            try:
                self.jatembase()
                self._output("Solução x =",'[m]solucao')
                self._output(self.solucao(),'[m]solucao')
            except Exception as e:
                self._output(str(e),'[m]error')
                self.preparar_output('fechar')
        try:
            self.preparar_output('fechar')
        except:
            pass

    def jatembase(self):
        '''Supondo que temos uma identidade nas últimas m colunas,
        e que 0 é solução inicial, executa o simplex.'''
        (m,n) = self.A.shape
        # Iniciamos com um tableau de zeros
        self.tableau = np.zeros((m+1,n+1))
        # Colocamos A, b e c no tableau
        self.tableau[0:m,0:n] = self.A
        self.tableau[0:m,n] = self.b
        # (lembrando sempre que os indices começam em 0)
        # Agora, colocamos c no tableau, completando com zeros
        self.tableau[m,0:len(self.c)] = -self.c
        self._output('Tableau','[m]tableau')
        self._output(self.tableau,'[m]jatembase')
        self.run()

    def fase1(self):
        '''Resolve o problema auxiliar min 1^Tx_a para encontrar uma primeira
        solução viável'''
        (m,n) = self.A.shape
        a = self.base.count(False)
        # quantidade de variaveis artificiais que adicionaremos
        self.artificiais = []
        self.tableau = np.zeros((m+1,n+a+1))
        # cria um tableau de zeros
        self.tableau[0:m,0:n] = self.A
        # Colocamos A, b  no tableau
        self.tableau[0:m,n+a] = self.b
        # (lembrando sempre que os indices começam em 0)
        for j in range(m):
            if self.base[j] == False:
                # para cada entrada False da base, adicionamos o vetor da
                # base canonica faltante
                self.tableau[0:m,n+j] = basecanonica(j,m)
                self.artificiais.append(n+j)
                self.base[j] = n+j
        self._output('Iniciando a fase 1.','[m]fase1')
        self.tableau[m,n:n+a] = -np.ones(a)
        self._output(self.tableau,'[m]custorelarti')
        self._output("Pivotando nas variáveis artificiais.",'[m]inipivotart')
        for j in self.artificiais:
            i = self.base.index(j)
            self.tableau[m,:] = self.tableau[m,:] + self.tableau[i,:]
            self._output(self.tableau,'pivotart')
        self._output(self.tableau,'[m]fimpivotart')
        self.run()
        # começamos supondo que nenhuma variável aritificial está na base
        self.deletar = []
        for k in self.base:
            if k in self.artificiais:
                # se alguma das variáveis artificiais pertencer a base
                i = self.base.index(k)
                if self.tableau[i,-1] != 0:
                    # e seu valor for diferente de zero,
                    # então não temos soluções viáveis
                    raise SemSolucoesViaveis()
                else:
                    # caso contrário, temos que tentar tirar a variável
                    # artificial da base
                    self.retirar_artificial_da_base(k)
        # se o método acima não conseguir remover a variável artificial
        # da base, significa que temos uma linha de zeros, que pode ser
        # eliminada
        if len(self.deletar)>0:
            self.tableau = np.delete(self.tableau, self.deletar,axis=0)

    def fase2(self):
        '''Método para rodar a fase 2 do simplex.'''
        self._output('Iniciando a fase 2.','[m]fase2')
        (M,N) = self.tableau.shape
        # aqui eu me ferrei muito até aprender que poderia ter acontecido de
        # termos eliminado linhas do tableau devido a variáveis artificiais
        # na base no final da fase 1.
        m = M - 1 #esse é o índice correto
        # primeiro, apagamos as coluas das variáveis artificiais
        self.tableau = np.delete(self.tableau, self.artificiais,axis=1)
        self.tableau[m,0:len(self.c)] = -self.c
        self._output('Vetor de custo relativo adicionado.','[m]fase2vetor')
        self._output(self.tableau,'[m]fase2vetor')
        # Pivotamos para zerar os custos relativos dos caras da base
        for j in self.base:
            i = self.tableau[0:m,j].argmax()
            self.tableau[m,:] = self.tableau[m,:] - self.tableau[m,j]*self.tableau[i,:]
            self._output(self.tableau,'zerarcustrel')
        self._output('Custo relativo dos elementos da base zerados.','zerarcustrel')
        self._output(self.tableau,'[m]custorelativozerado')
        self.run()

    def run(self):
        '''Roda a parte do pivotamento usado em cada fase do simplex.'''
        (M,N) = self.tableau.shape
        # usamos M,N maiúsuculos para o tamanho do tableau
        maximo = self.tableau[M-1,0:N-1].max()
        # enquanto o  maior elemento da ultima linha for maior que zero,
        # continua pivotando
        while  maximo > 0:
            quemEntra = self.tableau[M-1,0:N-1].argmax()
            # indice da coluna do maior elemento da ultima linha do tableau
            quemSai = self.quem_sai_da_base(quemEntra)
            # decide tem sai da base
            self.base[quemSai] = quemEntra
            # atualiza a base
            msg = "Pivotando na linha {} e coluna {}."
            self._output(msg.format(quemSai+1, quemEntra+1),'[m]entraesai')
            self.pivotar(quemSai,quemEntra)
            # pivota na linha e coluna indicadas
            self._output(self.tableau,'[m]run')
            maximo = self.tableau[M-1,0:N-1].max()
            # procura um novo cara pra entrar na base

    def quem_sai_da_base(self,k):
        '''Dado que x_k entra na base, essa função determina quem sai da base.'''
        (M,N) = self.tableau.shape
        # usando M,N maiúsculos para o tamanho do tableau para não confundir
        # com m,n tamanho da matriz A
        j = -1
        I = [i for i in range(M-1) if self.tableau[i,k] > 0]
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
            if j == N-1:
                # caso I_j contenha mais do que 1 elemento para todo j de 1 a n,
                # temos que o tableau tem pelo menos duas linhas l.d.
                raise LinhasLD(I)
        return I[0]

    def pivotar(self,i,j):
        '''Pivota o tableau na linha i e coluna j.'''
        (M,N) = self.tableau.shape
        self.tableau[i,:] = self.tableau[i,:]/self.tableau[i,j]
        # Divide a linha i pelo pivo
        self._output(self.tableau,'linhaipelopivo')
        for k in range(M):
            if k != i:
                self.tableau[k,:] = self.tableau[k,:] - self.tableau[k,j]*self.tableau[i,:]
                # Faz as eliminações nas outras linhas
                self._output(self.tableau, 'pivotando')

    def retirar_artificial_da_base(self,coluna):
        '''Tenta retirar a variável artificial da base. Se não conseguir,
        inclui ela na lista para eliminação.'''
        linha = self.base.index(coluna)
        (m,n) = self.A.shape
        naobasicas = [x for x in range(n) if x not in self.base]
        # guarda os índices das colunas das variáveis não-básicas
        # não-artificiais
        colunas = [j for j in naobasicas if self.tableau[linha,j] > 0]
        if len(colunas) > 0:
            # se existir alguém maior que zero para podermos pivotar
            self.base[linha] = colunas[0]
            msg = "Retirando a variável artifical {} da base."
            self._output(msg.format(coluna+1),'[m]retbase')
            self.pivotar(linha,colunas[0])
            self._output(self.tableau,'[m]retbase')
        else:
            # se não for possível retirar a variável artificial da base
            # diminuimos o tamanho da base pois teremos que eliminar linhas
            # do tableau
            self.base.remove(coluna)
            self.deletar.append(linha)

    def solucao(self):
        '''Retorna o vetor solução do problema.'''
        (M,N) = self.tableau.shape
        (m,n) = self.A.shape
        x = np.zeros(n)
        for i in range(M-1):
            x[self.base[i]] = self.tableau[i,-1]
        return x

    # Aqui é só coisa pra deixar o output bonitinho e converter pra latex
    def _output(self,thing,source):
        '''Imprime strings ou arrays nos devidos locais definidos em output.'''
        if source=='definirproblema':
            # vamos declarar as variáveis em jogo
            if 'screen' in self.output:
                print('A = \n',self.A,'\n')
                print('b = \n',self.b,'\n')
                print('c = \n',self.c,'\n')

            if 'latex' in self.output:
                texto = '\\begin{eqnarray*} \\min & '+ bmatrix(self.c)  +'x \\\\'
                texto = texto + '\\text{s. a.} &' + bmatrix(self.A) + 'x = '
                bempe = bmatrix(self.b.reshape((len(self.b),1)))
                texto = texto + bempe + '\end{eqnarray*}\n'
                self.latex.write(texto)
            if 'file' in self.output:
                self.file.write('A = \n'+str(self.A)+'\n')
                self.file.write('b = \n'+str(self.b)+'\n')
                self.file.write('c = \n'+str(self.c)+'\n')
        if 'minimal' in self.output:
            # Se a opção minimal for passada, temos menos output
            if 'screen' in self.output and '[m]' in source:
                print(thing)
            if 'latex' in self.output  and '[m]' in source:
                if type(thing) == str:
                    self.latex.write('\\par ' + thing + '\n')
                if type(thing) == np.ndarray:
                    matriz = bmatrix(thing)
                    self.latex.write('\\['+ matriz +'\\]\n')
            if 'file' in self.output  and '[m]' in source:
                self.file.write(str(thing)+'\n')
        else:
            # se não, imprimimos todas as informações
            if 'screen' in self.output:
                print(thing)
            if 'latex' in self.output:
                if type(thing) == str:
                    self.latex.write('\\par ' + thing + '\n')
                if type(thing) == np.ndarray:
                    matriz = bmatrix(thing)
                    self.latex.write('\\['+ matriz +'\\]\n')
            if 'file' in self.output:
                self.file.write(str(thing)+'\n')

    def preparar_output(self,mode):
        if mode=='abrir':
            if 'latex' in self.output:
                preamble = """\\documentclass[a4paper,12pt]{article}
\\usepackage[utf8]{inputenc}
\\usepackage[brazil]{babel}
\\usepackage{amsmath,amssymb,amsfonts}
\\begin{document}
"""

                start = self.output.index('latex=') + len('latex=')
                end = self.output.index('.tex') + len('.tex')
                filename = self.output[start:end]
                self.latex = open(filename,'w',encoding='utf-8')
                self.latex.write(preamble)
            if 'file' in self.output:
                start = self.output.index('file=') + len('file=')
                end = self.output.index('.txt') + len('.txt')
                filename = self.output[start:end]
                self.file = open(filename,'w',encoding='utf-8')
        if mode=='fechar':
            if 'latex' in self.output:
                self.latex.write('\\end{document}')
                self.latex.close()
            if 'file' in self.output:
                self.file.close()



# Fim da Classe Simplex

def basecanonica(j,n):
    '''Retorna o vetor nx1 e_j, da base canonica de R^n.
    PS: índices começam em 0.'''
    ej = np.zeros(n)
    ej[j] = 1
    return ej

def bmatrix(a):
    """Returns a LaTeX bmatrix

    :a: numpy array
    :returns: LaTeX bmatrix as a string
    """
    if len(a.shape) > 2:
        raise ValueError('bmatrix can at most display two dimensions')
    lines = str(a).replace('[', '').replace(']', '').splitlines()
    rv = [r'\begin{bmatrix}']
    rv += ['  ' + ' & '.join(l.split()) + r'\\' for l in lines]
    rv +=  [r'\end{bmatrix}']
    return '\n'.join(rv)


class LinhasLD(Exception):
    def __init__(self, linhas):
        self.linhas = linhas
        self.word = 'linhas'
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
        self.texto = 'As '+self.word+' ' + lnhs + 'sao linearmente dependentes.'
        self.texto = self.texto.format(*self.linhas)
    def __str__(self):
        return self.texto


class ColunasLD(LinhasLD):
    def __init__(self,colunas):
        self.linhas = colunas
        self.word = 'colunas'
        self.texto = ''
        self.txt()


class SemSolucoesViaveis(Exception):
    def __str__(self):
        return "O conjunto viável é vazio."


class ProblemaIlimitado(Exception):
    def __str__(self):
        return "O -gradiente aponta em uma direção ilimitada do conjunto."


class SimplexDual(SimplexPrimal):
    '''
SimplexDual
===========

Classe que contém todo a implementação do método simplex dual para resolver
o problema de programação linear na forma padrão
::
    min  <c,x>
    s.a. Ax  = b
          x >= 0

Parâmetros
----------
A : m x n array
    Matriz das restrições lineares, com as variáveis de folga já inclusas.
    Assumimos que a matriz A já contém a identidade m x m em suas últimas
    m colunas.
b : m array
    Lado direito das restrições lineares.
c : p array
    Gradiente da função objetivo, que deve necessariamente ser não-negativo.
    Se p <= n, o algoritmo completa c com zeros.
output : string
    String contendo as configurações de onde será feito o output do problema.
    Aceita as seguintes opções, separadas por vírgula.

    * `file=filename.txt`: imprime o output no arquivo filename.txt

    * `latex=filename.tex`: imprime o output em formato LaTeX no arquivo filename.tex

    * `screen`: imprime o output na tela

    * `minimal`: reduz o tamanho do output, não imprimindo todos os passos do pivotamento



Autor
-----

lfbossa@gmail.com, 2016
    '''

    def check(self):
        zeros = np.zeros_like(self.c)
        posicoes, = np.where(self.c < zeros)
        # retorna as posições em que b < 0
        if len(posicoes) > 0: # Se alguma entrada de c for negativa
            raise ValueError("O vetor c deve ser maior que zero.")
        (m,n) = self.A.shape


    def resolver(self):
        '''Supõe que as últimas m colunas de A formam uma base'''
        (m,n) = self.A.shape
        self.base = [j for j in range(n-m,n)]
        self.preparar_output('abrir')
        self._output('','definirproblema')
        try:
                self.jatembase()
                self._output("Solução x =",'[m]solucao')
                self._output(self.solucao(),'[m]solucao')
        except Exception as e:
                self._output(str(e),'[m]error')
                self.preparar_output('fechar')
        try:
            self.preparar_output('fechar')
        except:
            pass

    def run(self):
        '''Roda a parte do pivotamento usado em cada fase do simplex.'''
        (M,N) = self.tableau.shape
        # usamos M,N maiúsuculos para o tamanho do tableau
        minimo = self.tableau[0:M-1,N-1].min()
        # enquanto o  menor elemento do lado direito for menor que zero,
        # continua pivotando
        while  minimo < 0:
            quemSai = self.tableau[0:M-1,N-1].argmin()
            # indice  da linha do menor elemento do lado direito do tableau
            quemEntra = self.quem_entra_na_base(quemSai)
            # decide tem entra da base
            self.base[quemSai] = quemEntra
            # atualiza a base
            msg = "Pivotando na linha {} e coluna {}."
            self._output(msg.format(quemSai+1, quemEntra+1),'[m]entraesai')
            self.pivotar(quemSai,quemEntra)
            # pivota na linha e coluna indicadas
            self._output(self.tableau,'[m]run')
            minimo = self.tableau[0:M-1,N-1].min()
            # procura um novo cara pra sair na base

    def quem_entra_na_base(self,k):
        '''Dado que x_k sai da base, essa função determina quem
        entra na base.'''
        (M,N) = self.tableau.shape
        # usando M,N maiúsculos para o tamanho do tableau para não confundir
        # com m,n tamanho da matriz A
        i = -1
        I = [j for j in range(N-1) if self.tableau[k,j] < 0]
        # I-1 é o conjunto dos indices j tais que y_kj < 0
        if len(I) == 0:
            # se y_k >= 0, o problema é ilimitado
            raise ProblemaIlimitado()
        razoes = np.array([self.tableau[i,j]/self.tableau[k,j] for j in I])
        # calcula as razoes z_j-c_j/y_kj para y_kj < 0
        razaominima = razoes.min()
        # encontra o valor da razao mínima
        I = [I[j] for j in range(razoes.size) if razoes[j] == razaominima]
        # conjunto I0 dos indices r tais que z_r-c_r/y_kr é igual a razao minima
        # se I0  tiver mais de  1 elemento, usamos a validação lexicográfica
        while len(I) > 1:
            i = i + 1 # calcula os conjuntos I_i, para i <= m
            razoes = np.array([self.tableau[i,j]/self.tableau[k,j] for j in I])
            razaominima = razoes.min()
            I = [I[j] for j in range(razoes.size) if razoes[j] == razaominima]
            if i == M-1:
                # caso I_i contenha mais do que 1 elemento para todo i de 1 a m,
                # temos que o tableau tem pelo menos duas colunas l.d.
                raise ColunasLD(I)
        return I[0]
