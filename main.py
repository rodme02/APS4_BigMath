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
deformacao_nodal = np.zeros((nm*2, 1))
j = 0
for i in range(nm*2):
    if i not in R:
        deformacao_nodal[i] = solucao_gauss[j]
        j += 1

# Reações de apoio
PG = KG.dot(deformacao_nodal)
R1x = PG[0][0]
R2x = PG[2][0]
R2y = PG[3][0]
print("R1x = ", R1x)
print("R2x = ", R2x)
print("R2y = ", R2y)

