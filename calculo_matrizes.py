from funcoesTermosol import *
import numpy as np
from math import *
from solucoes import *


def calcula_KG(Inc,N): 
    nn = N.shape[1]
    elementos = Inc.shape[0]
    lista_KE = []
    KG = np.zeros((nn*2,nn*2))

    for i in range(elementos):
        E = Inc[i,2]
        A = Inc[i,3]
        n1 = int(Inc[i,0])
        n2 = int(Inc[i,1])
        
        x1 = N[0,n1-1]
        y1 = N[1,n1-1]
        x2 = N[0,n2-1]
        y2 = N[1,n2-1]
        
        L = sqrt((x2-x1)**2+(y2-y1)**2)
        c = (x2-x1)/L
        s = (y2-y1)/L
        
        KE = ((E*A)/L)*np.array([  [c**2,c*s,-c**2,-c*s],
                                    [c*s,s**2,-c*s,-s**2],
                                    [-c**2,-c*s,c**2,c*s],
                                    [-c*s,-s**2,c*s,s**2]])
        
        lista_KE.append(KE)
        
        KG[int(Inc[i,0])*2-2:int(Inc[i,0])*2,int(Inc[i,0])*2-2:int(Inc[i,0])*2] += KE[0:2,0:2]
        KG[int(Inc[i,0])*2-2:int(Inc[i,0])*2,int(Inc[i,1])*2-2:int(Inc[i,1])*2] += KE[0:2,2:4]
        KG[int(Inc[i,1])*2-2:int(Inc[i,1])*2,int(Inc[i,0])*2-2:int(Inc[i,0])*2] += KE[2:4,0:2]
        KG[int(Inc[i,1])*2-2:int(Inc[i,1])*2,int(Inc[i,1])*2-2:int(Inc[i,1])*2] += KE[2:4,2:4]

    return KG

def reduz_matrizes(KG, F, R):
    KG_reduzido = np.delete(KG, R, axis=0)
    KG_reduzido = np.delete(KG_reduzido, R, axis=1)
    F_reduzido = np.delete(F, R, axis=0)

    return (KG_reduzido, F_reduzido)

def sistema_de_equacoes(KG_reduzido, F_reduzido, nn, R):
    U = gauss_seidel(1000, 1e-100, KG_reduzido, F_reduzido)

    deslocamento = np.zeros((nn*2,1))

    j = 0
    for i in range(nn*2):
        if i not in R:
            deslocamento[i] = U[j]
            j += 1
    return deslocamento

def reacao_de_apoio(deslocamento, kg, R):
    PG = kg.dot(deslocamento)
    
    reacoes_de_apoio = PG[R]
    reacoes_de_apoio = np.array([reacoes_de_apoio]).T
    return reacoes_de_apoio

def deformacao_tensao_forca_interna(Inc, N, deslocamentos):
    lista_deformacoes = []
    lista_tensoes = []
    lista_forcas_internas = []
    for i in range(Inc.shape[0]):
        E = Inc[i,2]
        A = Inc[i,3]
        n1 = int(Inc[i,0]) - 1
        n2 = int(Inc[i,1]) - 1
        
        x1 = N[0,n1]
        y1 = N[1,n1]
        x2 = N[0,n2]
        y2 = N[1,n2]
        
        L = sqrt((x2-x1)**2+(y2-y1)**2)
        c = (x2-x1)/L
        s = (y2-y1)/L

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

    return np.array([lista_deformacoes]).T, np.array([lista_tensoes]).T, np.array([lista_forcas_internas]).T