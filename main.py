from funcoesTermosol import *
from calculo_matrizes import *
from solucoes import *

import numpy as np

[nn,N,nm,Inc,nc,F,nr,R] = importa('entrada.xls')

lista_KE = calcula_KE(Inc, N)
KG = calcula_KG(lista_KE, Inc)

KG_com_restricoes = CondicoesContorno_KG(KG, R)
# print(F)
print(KG)
print(R)
print(nr)
print(KG_com_restricoes)

