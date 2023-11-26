from funcoesTermosol import *
import numpy as np


def CalculaParametros(Inc, N):
    lista_cossenos = []
    lista_senos = []
    lista_L = []
    lista_E = []
    lista_A = []
    lista_n1 = []
    lista_n2 = []

    for i in range(Inc.shape[0]):
        # Nós do elemento
        n1 = int(Inc[i, 0]) - 1
        n2 = int(Inc[i, 1]) - 1

        # Módulo de Young do elemento
        E = Inc[i, 2]

        # Área da seção transversal do elemento
        A = Inc[i, 3]

        # Coordenadas dos nós do elemento
        x1, y1 = N[0, n1], N[1, n1]
        x2, y2 = N[0, n2], N[1, n2]
        
        # Comprimento do elemento
        L = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)

        # Cosseno e seno do ângulo de inclinação do elemento
        c = (x2 - x1) / L
        s = (y2 - y1) / L
        lista_L.append(L)
        lista_cossenos.append(c)
        lista_senos.append(s)
        lista_E.append(E)
        lista_A.append(A)
        lista_n1.append(n1)
        lista_n2.append(n2)


    return lista_cossenos, lista_senos, lista_L, lista_E, lista_A, lista_n1, lista_n2

# Função para verificar se uma matriz é simétrica
def is_symmetric(arr):
    arr = np.array(arr)
    return np.array_equal(arr, arr.T)

def calcula_KE(Inc, N):
    lista_cossenos, lista_senos, lista_L, lista_E, lista_A, n1, n2 = CalculaParametros(Inc, N)
    lista_KE = []
    
    # Loop para calcular a matriz de rigidez local de cada elemento
    for i in range(Inc.shape[0]):
        # Nós do elemento
        c = lista_cossenos[i]
        s = lista_senos[i]
        L = lista_L[i]
        E = lista_E[i]
        A = lista_A[i]

        # Matriz de rigidez local do elemento
        KE = (E * A / L) * np.array(   [[c**2, c*s, -c**2, -c*s],
                                        [c*s, s**2, -c*s, -s**2],
                                        [-c**2, -c*s, c**2, c*s],
                                        [-c*s, -s**2, c*s, s**2]])

        lista_KE.append(KE)

    return lista_KE

def calcula_KG(lista_KE, Inc, nn):
  # Inicializa a matriz de rigidez global com dimensão 2n x 2n
  KG = np.zeros((2*nn, 2*nn))

  # Loop para adicionar as matrizes de rigidez locais à matriz de rigidez global
  for i in range(len(lista_KE)):
    KE = lista_KE[i]

    # Nós do elemento
    n1 = int(Inc[i, 0]) - 1
    n2 = int(Inc[i, 1]) - 1
    
    # Índices dos nós do elemento na matriz de rigidez global
    i1 = 2*n1
    i2 = 2*n1 + 1
    i3 = 2*n2
    i4 = 2*n2 + 1
    
    # Adiciona as matrizes de rigidez locais à matriz de rigidez global
    KG[i1:i2+1, i1:i2+1] += KE[0:2, 0:2]
    KG[i1:i2+1, i3:i4+1] += KE[0:2, 2:4]
    KG[i3:i4+1, i1:i2+1] += KE[2:4, 0:2]
    KG[i3:i4+1, i3:i4+1] += KE[2:4, 2:4]
  
  return KG

def CondicoesContorno_KG(KG, R):
    KG = np.delete(KG, R, axis=0)
    KG = np.delete(KG, R, axis=1)
    return KG

def CondicoesContorno_F(F, R): 
    F = np.delete(F, R, axis=0)
    return F
    
