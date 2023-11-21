from funcoesTermosol import *
import numpy as np


def calcula_KE(E, A, N, nn):
    """
    Calcula a matriz de rigidez de um elemento de barra
    """
    x = []
    y = []
    L = []
    c = []
    s = []
    KE = []

    #calcula as coordenadas dos nos
    for i in range(nn):
        x.append(N[0,i])
        y.append(N[1,i])

    print(f"X: {x} \n Y: {y}")

    for i in range(len(x)):
        #calcula o comprimento do elemento
        if i != nn-1:
            L.append(np.sqrt((x[i+1]-x[i])**2+(y[i+1]-y[i])**2))
        else:
            L.append(np.sqrt((x[0]-x[i])**2+(y[0]-y[i])**2))
    
    print(f"L: {L} \n")

    #calcula seno e cosseno
    for i in range(len(L)):
        if i != nn-1:
            c.append((x[i+1]-x[i])/L[i])
            s.append((y[i+1]-y[i])/L[i])
        else:
            c.append((x[0]-x[i])/L[i])
            s.append((y[0]-y[i])/L[i])

    print(f"C: {c} \n S: {s} \n")

    # Calcula a matriz de rigidez
    for i in range(nn):
        KE.append((E*A/L[i])*np.array([[c[i]**2, c[i]*s[i], -c[i]**2, -c[i]*s[i]],
                                     [c[i]*s[i], s[i]**2, -c[i]*s[i], -s[i]**2],
                                     [-c[i]**2, -c[i]*s[i], c[i]**2, c[i]*s[i]],
                                     [-c[i]*s[i], -s[i]**2, c[i]*s[i], s[i]**2]]))
        print(f"\n K{i+1}: {KE[i]} \n")
        print(f"K{i+1} simetrica: {is_symmetric(KE[i])} \n")
        print("--------------------------------------------------")

    return KE


def is_symmetric(arr):
    # Convert the array to a NumPy array for easy manipulation
    arr = np.array(arr)
    
    # Check if the array is symmetric by comparing with its transpose
    return np.array_equal(arr, arr.T)


def calcula_K(nm, nn, Inc, N, E, A):
    """
    Calcula a matriz de rigidez global
    """
    # Inicializa a matriz de rigidez global
    K = np.zeros((nn*2,nn*2))
    
    # Loop sobre os membros
    for e in range(nm):
        # Obtem os nos do elemento
        no1 = int(Inc[e,0])
        no2 = int(Inc[e,1])
        
        # Obtem as coordenadas dos nos
        x1 = N[0,no1-1]
        y1 = N[1,no1-1]
        x2 = N[0,no2-1]
        y2 = N[1,no2-1]
        
        # Calcula o comprimento do elemento
        L = np.sqrt((x2-x1)**2+(y2-y1)**2)
        
        # Calcula o cosseno e o seno
        c = (x2-x1)/L
        s = (y2-y1)/L
        
        # Calcula a matriz de rigidez do elemento
        KE = calcula_KE(E, A, L, c, s)
        
        # Obtem os graus de liberdade do elemento
        GDL = np.array([no1*2-1, no1*2, no2*2-1, no2*2])
        
        # Adiciona a contribuicao do elemento a matriz de rigidez global
        for i in range(4):
            for j in range(4):
                K[GDL[i]-1,GDL[j]-1] = K[GDL[i]-1,GDL[j]-1] + KE[i,j]
    
    return K

def calcula_comprimento(nm, Inc, N):
    """
    Calcula o comprimento de cada elemento
    """
    # Inicializa o vetor de comprimentos
    L = np.zeros((nm,1))
    
    # Loop sobre os membros
    for e in range(nm):
        # Obtem os nos do elemento
        no1 = int(Inc[e,0])
        no2 = int(Inc[e,1])
        
        # Obtem as coordenadas dos nos
        x1 = N[0,no1-1]
        y1 = N[1,no1-1]
        x2 = N[0,no2-1]
        y2 = N[1,no2-1]
        
        # Calcula o comprimento do elemento
        L[e] = np.sqrt((x2-x1)**2+(y2-y1)**2)
    
    return L

def calcula_theta(nm, Inc, N):
    """
    Calcula o angulo de cada elemento
    """
    # Inicializa o vetor de angulos
    theta = np.zeros((nm,1))
    
    # Loop sobre os membros
    for e in range(nm):
        # Obtem os nos do elemento
        no1 = int(Inc[e,0])
        no2 = int(Inc[e,1])
        
        # Obtem as coordenadas dos nos
        x1 = N[0,no1-1]
        y1 = N[1,no1-1]
        x2 = N[0,no2-1]
        y2 = N[1,no2-1]
        
        # Calcula o angulo do elemento
        theta[e] = np.arctan((y2-y1)/(x2-x1))
    
    return theta


# print('Lendo o arquivo de entrada')
[nn,N,nm,Inc,nc,F,nr,R] = importa('entrada.xls')

calcula_KE(Inc[0,2], Inc[0,3], N, nn )