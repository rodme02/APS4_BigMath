from funcoesTermosol import *
from calculo_matrizes import *
from solucoes import *

import numpy as np

[nn,N,nm,Inc,nc,F,nr,R] = importa('entrada.xls')

lista_KE = calcula_KE(Inc, N)
KG = calcula_KG(lista_KE, Inc)

KG_com_restricoes = CondicoesContorno_KG(KG, R)
F_com_restricoes = CondicoesContorno_F(F, R)
#print(F)
#print(KG)
#print(R)
#print(nr)
#print(KG_com_restricoes)
#print(F_com_restricoes)

#Tolerancias usadas em ambas foi recomendada pelo professor
jacobi_feito = jacobi(100, 10e-10, KG_com_restricoes, F_com_restricoes)
#print(jacobi_feito)

# Vamos usar o gauss_seidel porque o resultado dele possui menos iterações e porque 
gauss_seidel_feito = gauss_seidel(100, 10e-10, KG_com_restricoes, F_com_restricoes)
#print(gauss_seidel_feito)

deslocamentos =[]
for i in gauss_seidel_feito[0]:
    deslocamentos.append(i[0])

print(deslocamentos)

#controle = np.linalg.solve(KG_com_restricoes, F_com_restricoes)
#print(controle)