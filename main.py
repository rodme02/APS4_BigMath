from funcoesTermosol import *
from calculo_matrizes import *
from solucoes import *

import numpy as np

[nn,N,nm,Inc,nc,F,nr,R] = importa('entrada.xls')

lista_KE = calcula_KE(Inc, N)
KG = calcula_KG(lista_KE, Inc)

R = np.array(R).flatten().astype(int)

KG_com_restricoes = CondicoesContorno_KG(KG, R)
F_com_restricoes = CondicoesContorno_F(F, R)

#Tolerancias usadas em ambas foi recomendada pelo professor
solucao_jacobi, erro_max_jacobi, iteracoes_jacobi = jacobi(100, 10e-10, KG_com_restricoes, F_com_restricoes)

# Vamos usar o gauss_seidel porque o resultado dele possui menos iterações e porque o erro máximo é menor
solucao_gauss, erro_max_gauss, iteracoes_gauss = gauss_seidel(100, 10e-10, KG_com_restricoes, F_com_restricoes)

# Deformações nodais
deslocamentos= np.zeros((nm*2, 1))
j = 0
for i in range(nm*2):
    if i not in R:
        deslocamentos[i] = solucao_gauss[j]
        j += 1

# Reações de apoio
PG = KG.dot(deslocamentos)
R1x = PG[0][0]          # Reação de apoio no nó 1 na direção x
R2x = PG[2][0]          # Reação de apoio no nó 2 na direção x
R2y = PG[3][0]          # Reação de apoio no nó 2 na direção y
Reacoes = np.array([R1x, R2x, R2y])

# Tensões nos elementos
lista_deformacoes = []
lista_tensoes = []
lista_forcas_internas = []
lista_cossenos, lista_senos, lista_L, lista_E, lista_A, lista_n1, lista_n2 = CalculaParametros(Inc, N) # Calcula os parâmetros dos elementos

for i in range(nm):
        n1 = lista_n1[i]
        n2 = lista_n2[i]

        x1, y1 = N[0, n1], N[1, n1]
        x2, y2 = N[0, n2], N[1, n2]
        
        L = lista_L[i]

        c = lista_cossenos[i]
        s = lista_senos[i]

        u1 = deslocamentos[2*n1][0]     # Deslocamento na direção x do nó 1
        v1 = deslocamentos[2*n1+1][0]   # Deslocamento na direção y do nó 1
        u2 = deslocamentos[2*n2][0]     # Deslocamento na direção x do nó 2
        v2 = deslocamentos[2*n2+1][0]   # Deslocamento na direção y do nó 2

        deformacao = (1/L)*np.array([-c, -s, c, s]).dot(np.array([u1, v1, u2, v2])) # Deformação no elemento
        tensao = Inc[i, 2]*deformacao                   # Tensão no elemento
        forca_interna = tensao*Inc[i, 3]                # Força interna no elemento

        lista_deformacoes.append(deformacao)
        lista_tensoes.append(tensao)
        lista_forcas_internas.append(forca_interna)

geraSaida("saida.txt",Reacoes, deslocamentos, lista_deformacoes, lista_forcas_internas, lista_tensoes)
