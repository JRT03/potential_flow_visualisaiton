import numpy as np
import matplotlib.pyplot as plt
plt.style.use('dark_background')

x_array = np.linspace(-10,10,1000)
y_array = np.linspace(-5,5,500)

X_mesh = np.zeros([len(x_array),len(y_array)])
Y_mesh = np.zeros([len(x_array),len(y_array)])

for i in range(len(x_array)):
    for j in range(len(y_array)):
        X_mesh[i][j] = x_array[i]
        Y_mesh[i][j] = y_array[j]        

U = 2
s = 1

phi = U*(Y_mesh + s*np.arccos(X_mesh/np.sqrt(X_mesh**2 + Y_mesh**2)))
C_p = - ( ((2*s)/(X_mesh**2 + Y_mesh**2))*((X_mesh)/(np.sqrt(X_mesh**2 + Y_mesh**2))) + ((s**2)/(X_mesh**2 + Y_mesh**2)) )

phi_mask = np.ones(phi.shape)
for i in range(phi_mask.shape[0]):
    for j in range(phi_mask.shape[1]):
        if np.abs(phi[i][j]) < np.pi*s*U:
            phi_mask[i][j] = 0

C_p = C_p * phi_mask

contour_list = []
for i in range(10,30,2):
    contour_list.append(np.pi*U*s*(i/10))

plt.contour(X_mesh,Y_mesh,phi,contour_list)
plt.contour(X_mesh,-Y_mesh,phi,contour_list)
plt.title('Flow around Rankine Half Body')
plt.xlabel('x')
plt.ylabel('y')
plt.show()


contour_levels = np.linspace(-1, 1, 100)  
plt.contour(X_mesh,Y_mesh,phi,contour_list,colors='black')
plt.contourf(X_mesh, Y_mesh, C_p, levels=contour_levels)
#plt.ylim(0,4*s)
plt.colorbar()
plt.title('Pressure Field Around Rankine Half Body')
plt.xlabel('x')
plt.ylabel('y')
plt.show()



## from equation of stremlines, for a given y find the x that satisfys the equation then plot these sets of points (x that is closest to it)
## problems will be run into long y = 0 ! there will be multiple points with same value could plot all points with same lowst val, still problem at x= y= 0 where theta undefined