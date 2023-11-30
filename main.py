from funcoesTermosol import *
from calculo_matrizes import *
from solucoes import *

import numpy as np

[nn,N,nm,Inc,nc,F,nr,R] = importa('validacao.xls')

KG = calcula_KG(Inc, N)

R = R.flatten().astype(int)

KG_com_restricoes, F_com_restricoes = reduz_matrizes(KG, F, R)

deslocamentos = sistema_de_equacoes(KG_com_restricoes, F_com_restricoes, nn, R)

reacoes_de_apoio = reacao_de_apoio(deslocamentos, KG, R)

lista_deformacoes, lista_tensoes, lista_forcas_internas = deformacao_tensao_forca_interna(Inc, N, deslocamentos)

geraSaida("saida.txt", reacoes_de_apoio, deslocamentos, lista_deformacoes, lista_forcas_internas, lista_tensoes)

plota(N, Inc)

N_novo = np.array(N)
for i in range(len(N_novo[0])):
        N_novo[0, i] += deslocamentos[i*2]
        N_novo[1, i] += deslocamentos[i*2+1]

plota(N_novo, Inc)