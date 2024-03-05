import matplotlib.pyplot as plt
import numpy as np


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np


a = float(input('a: '))
b = float(input('b: '))
c = float(input('c: '))
d = float(input('d: '))
N = int(input('N: '))
M = int(input('M: '))


h = (b-a) / N
k = (d-c) / M

w = np.zeros((N,M))
v = float(input('Ingrese velocidad: '))
p = v*k/h

def f(x):
    return 0

def g(x):
    return 0

for j in range(1,M):
    w[j][0]=0
    w[j][N]=0

for i in range(1,N):
    w[0][i]=f(h*i)
    w[1][i]=w[0][i]+k*g(h*i)

for j in range(0,M):
    for i in range(0,N):
        w[j+1][i]=2*(1-p**2)*w[j][i]+(p**2)(w[j][i+1]+w[j][i-1])-w[j-1][i]

#definir los puntos x,y,z para la superficie
x=np.linspace(a,b,N+1)
y=np.linspace(c,d,M+1)
X,Y=np.meshgrid(x,y)

#graficar la superficie
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X,Y,w,cmap='viridis', edgecolor='none')
ax.set_xlabel('X')
ax.set_ylabel('Y ')
ax.set_zlabel('Z')
ax.set_title('Solución de la EDP con Ondas')
plt.show()

# a=0 b=5 c=0 d=10 N=40 M=400 v=0.5




"""
    '''func = u_xx + u_yy'''

    # Gauss-Seidel para la solución de la EDP
    for _ in range(100):
        for i in range(1, N):
            for j in range(1, M):
                w[i][j] = (k**2 *(w[i+1][j] + w[i-1][j]) + h**2 *(w[i][j+1] + w[i][j-1]) - h**2 * k**2 * func(i,j)) / (2*(k**2 + h**2))


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

"""
