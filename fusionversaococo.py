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



def assemble_global_stiffness_matrix(KE_matrices, nn):
    """
    Assembles the global stiffness matrix from a list of individual element stiffness matrices.
    """
    # Total number of nodes
    total_nodes = nn * len(KE_matrices)

    # Initialize the global stiffness matrix
    K_global = np.zeros((total_nodes, total_nodes))

    for i, KE in enumerate(KE_matrices):
        # Determine the size of the element stiffness matrix
        size_i, size_j = KE.shape

        # Assemble the element stiffness matrix into the global stiffness matrix
        for j in range(size_i):
            i = i % nn
            print(f"i: {i} \n")
            for k in range(size_j):
                # Assemble each term into the corresponding location in the global matrix
                j = j % nn
                print(f"j: {j} \n")
                K_global[nn * i + j, nn * i + k] += KE[j, k]

    return K_global

# print('Lendo o arquivo de entrada')
[nn,N,nm,Inc,nc,F,nr,R] = importa('entrada.xls')

KE = calcula_KE(Inc[0,2], Inc[0,3], N, nn )
K = assemble_global_stiffness_matrix(KE, nn)