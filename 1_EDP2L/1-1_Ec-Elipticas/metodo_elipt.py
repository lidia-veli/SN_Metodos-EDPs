   
import matplotlib.pyplot as plt
import numpy as np


def gauss_seidel_elipt_poisson(a,b,c,d, N,M, h,k, w, f_lapl):
    '''f_lapl = u_xx + u_yy'''

    # Gauss-Seidel para la solución de la EDP
    for _ in range(100):
        for i in range(1, N):
            for j in range(1, M):
                w[i][j] = (k**2 *(w[i+1][j] + w[i-1][j]) + h**2 *(w[i][j+1] + w[i][j-1]) - h**2 * k**2 * f_lapl(i,j)) / (2*(k**2 + h**2))


    # Convertir 'w' a un array de NumPy para la visualización
    w_np = np.array(w)

    # malla para las coordenadas x e y
    x = np.linspace(a, b, N+1)
    y = np.linspace(c, d, M+1)
    X, Y = np.meshgrid(x, y)

    # Visualización en 3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(X, Y, w_np, cmap='viridis', edgecolor='none')
    ax.set_xlabel('eje X')
    ax.set_ylabel('eje Y')
    ax.set_zlabel('U(x, y)')
    ax.set_title('Solución de la EDP con Gauss-Seidel')

    # Añadir una barra de colores que mapea los valores a colores
    fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.show()



def gauss_seidel_elipt_helmholtz(a,b,c,d, N,M, h,k, w, lamb, func):
    '''el codigo es igual pero sale un -h**2*k**2*l**2 en el denominador de la ecuacion de relajacion ecuacion de Helmholtz
    u_xx + u_yy + lambda*u = func(x,y)
    donde f(x,y) = 0
    u(0,y) = u(x,0) = u(x, 1) = 0
    u(1, y) = 1
    vamos a hacer 3 casos
    1. lambda = 1
    2. lambda = 300
    3. lambda = 1000
    '''

        
    #recorremos los puntos interiores de la malla
    for _ in range(100): # iteramos 100 veces. de momento, luego pondremos condiciones de parada
        for i in range(1, N):
            for j in range(1, M):
                w[i][j] = (k**2 * (w[i+1][j] + w[i-1][j]) + h**2 * (w[i][j+1] + w[i][j-1]) - h**2 * k**2 * func(i,j)) / (2*(h**2 + k**2) - (h*k*lamb)**2) 


    # Convertir 'w' a un array de NumPy para la visualización
    w_np = np.array(w)

    # malla para las coordenadas x e ys
    x = np.linspace(a, b, N+1)
    y = np.linspace(c, d, M+1)
    X, Y = np.meshgrid(x, y)

    # Grafica la función como una superficie en 3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, w_np, cmap='viridis')
    ax.set_xlabel('eje X')
    ax.set_ylabel('eje Y')
    ax.set_zlabel('U(x, y)')
    ax.set_title('Solución de la EDP con Gauss-Seidel')

    # Muestra la gráfica
    plt.show()

