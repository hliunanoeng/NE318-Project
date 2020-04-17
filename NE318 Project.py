from __future__ import division
from mpl_toolkits.mplot3d import Axes3D
import sympy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import sys
sys.setrecursionlimit(10000)
# all necessary modules imported from above

#the tolerance value is set to be less than 1% of the maximum peak in the surface plot
def Project(m,n,t,type,x_num=50,y_num=50,tolerance=0.001):
    #save the input
    M=m
    N=n
    tt=t
    #symbolize the parameters
    [t,x1,x2,m,n]=sympy.symbols(["t","x1","x2","m","n"])
    #divide the final equation into four parts, integrate seperately, then recombine in the end
    T=0
    f1=sympy.cos(m*(sympy.pi)*x1)*((x1**4)-(2*(x1**3))+(x1**2)+1)
    f2=sympy.sin(n*(sympy.pi)*x2)*((x2**2)-x2)
    f3=(sympy.cos(m*(sympy.pi)*x1))**2
    f4=(sympy.sin(n*(sympy.pi)*x2))**2
    f1=sympy.integrate(f1,(x1,0,1))
    f2=sympy.integrate(f2,(x2,0,1))
    f3=sympy.integrate(f3,(x1,0,1))
    f4=sympy.integrate(f4,(x2,0,1))
    A_mn=(f1*f2)/(f3*f4)
    #for every different mn pair, calculate the sum of T equations
    for i in range(0,M):
        for j in range(1,N):
            T_ij=A_mn.subs({m:i,n:j})*(sympy.exp(-1*((i**2)+(j**2))*(((sympy.pi)**2)*t)))*(sympy.cos(i*(sympy.pi)*x1))*(sympy.sin(j*(sympy.pi)*x2))
            T=T+T_ij
    T=T*(-1) #the surface plot has a negative minimum value; it is my personal preference to flip the plot, making all values positive and creating a local maximum peak
    T_test=T.subs({t:0})
    T_exact=(-1)*((x1**4)-(2*(x1**3))+(x1**2)+1)*((x2**2)-x2)

    #make a grid to compare the real and approximated values
    x=np.linspace(0,1,x_num)
    y=np.linspace(0,1,y_num)
    i_matrix=np.zeros(shape=(x_num,y_num))
    matrix_exact=np.zeros(shape=(x_num,y_num))

    #calculate the real and approximated value at each grid point
    for ix in range(x_num):
        for iy in range(y_num):
            i_matrix[ix,iy]=T_test.subs({x1:x[ix],x2:y[iy]}).evalf()
            matrix_exact[ix,iy]=T_exact.subs({x1:x[ix],x2:y[iy]}).evalf()

    #the average error equals to the total sum of absolute values of the difference, then divided by the number of grid points
    error=(((np.absolute(i_matrix-matrix_exact))).sum()/(x_num * y_num))

    #only proceed if the error value is smaller than the tolerance,
    if (error>tolerance):
        raise ValueError("Tolerance test failed.")
        #sys.exit()
    else:
        print("Tolerance test passed.")

    #surface plot
    if type=="surface":
        T=T.subs({t:tt})
        sympy.plotting.plot3d(T,(x1,0,1),(x2,0,1),xlabel="x1",ylabel="x2",title="t ="+str(tt))
    #contour plot
    elif type=="contour":
        T=T.subs({t:tt})
        i_matrix=np.zeros(shape=(x_num,y_num))
        for ix in range(x_num):
            for iy in range(y_num):
                i_matrix[ix,iy]=T.subs({x1:x[ix],x2:y[iy]}).evalf()
        fig=plt.figure()
        ax=fig.add_subplot(111)
        ax.contourf(x,y,i_matrix)
        ax.set_ylabel('x2')
        ax.set_xlabel('x1')
        ax.set_title("t ="+str(tt))
        plt.show()
    #streamline plot
    elif type=="streamline":
        T=T.subs({t:tt})
        u = sympy.lambdify((x1, x2), T.diff(x2), 'numpy')
        v = sympy.lambdify((x1, x2), T.diff(x1), 'numpy')
        Y, X =  np.ogrid[0:1:100j, 0:1:100j]
        fig, ax = plt.subplots(figsize=(5, 5))
        ax.streamplot(X, Y, u(X, Y), v(X, Y))
        ax.set_ylabel('X2')
        ax.set_xlabel('X1')
        ax.set_title("t ="+str(tt))
        plt.show()
    #arrow plot
    elif type=="arrow":
        T=T.subs({t:tt})
        u = sympy.lambdify((x1, x2), T.diff(x2), 'numpy')
        v = sympy.lambdify((x1, x2), T.diff(x1), 'numpy')
        x = np.linspace(0,1,10)
        X,Y =np.meshgrid(x,x)
        plt.quiver(X, Y, u(X, Y), v(X, Y))
        plt.show()
    else:
        raise ValueError("Input graph type not recognizable, please enter a value in all lower case.")

#Project(10,10,0,"surface")
Project(10,10,1,"surface")
#Project(10,10,0,"contour")
#Project(10,10,1,"contour")
#Project(10,10,0,"streamline")
#Project(10,10,1,"streamline")
#Project(10,10,0,"arrow")
#Project(10,10,1,"arrow")
