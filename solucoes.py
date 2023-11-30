import numpy as np

def jacobi(ite, tol, K, F):
    x1_old, x2_old, x3_old = 0, 0, 0
    iteracoes = 0
    erro_max = np.inf
    for i in range(ite):
        x1 = (F[0] - K[0,1]*x2_old - K[0,2]*x3_old)/K[0,0]
        x2 = (F[1] - K[1,0]*x1_old - K[1,2]*x3_old)/K[1,1]
        x3 = (F[2] - K[2,0]*x1_old - K[2,1]*x2_old)/K[2,2]
        
        erro_x1 = np.abs((x1 - x1_old) / x1) if x1 != 0 else np.abs(x1 - x1_old)
        erro_x2 = np.abs((x2 - x2_old) / x2) if x2 != 0 else np.abs(x2 - x2_old)
        erro_x3 = np.abs((x3 - x3_old) / x3) if x3 != 0 else np.abs(x3 - x3_old)
        
        erros = np.array([erro_x1, erro_x2, erro_x3])
        erro_max = np.max(erros)
        
        if erro_max <= tol:
            break
        
        x1_old = x1
        x2_old = x2
        x3_old = x3
        iteracoes = i + 1

    solucao = np.array([x1, x2, x3])
    return solucao, erro_max, iteracoes
 
def gauss_seidel(ite, tol, K, F):
    n = len(F)
    U = np.zeros(n)
    Ei = np.inf
    cont_ite = 0

    for i in range(ite):
        U_old = np.copy(U)
        erros = []
        for i in range(n):
            s1 = sum(K[i][j] * U[j] for j in range(i))
            s2 = sum(K[i][j] * U_old[j] for j in range(i + 1, n))
            U[i] = (F[i] - s1 - s2) / K[i][i]
            erro = abs((U[i] - U_old[i]) / U[i]) if U[i] != 0 else abs(U[i] - U_old[i])
            erros.append(erro)

        cont_ite += 1

        Ei = max(erros)

        if Ei < tol:
            break


    return U
