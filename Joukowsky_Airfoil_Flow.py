import numpy as np
import matplotlib.pyplot as plt
plt.style.use('dark_background')


'''Defining flow parameters'''
U = 2.0
alpha = (np.pi/180)*0
a = 2.0
c = 0.1
f = 0.1
Gam = -4*np.pi*U*a*(f*np.cos(alpha) + (1+c)*np.sin(alpha))

'''Generating complex plane'''
x_array = np.linspace(-6*a, 6*a, 600)
y_array = np.linspace(-5*a, 5*a, 500)
X_mesh, Y_mesh = np.meshgrid(x_array, 1j * y_array)
Z_mesh = X_mesh + Y_mesh

'''Apply conformal map to create second plane'''
Zeta_mesh = Z_mesh + (a**2)/(Z_mesh)

'''Define the complex potential'''
F_pot = U * (Z_mesh + a*c - a*f*1j) * (np.cos(alpha) - 1j * np.sin(alpha)) + U * ((a**2)*(f**2 + (1+c)**2))/(Z_mesh + a*c - a*f*1j) * (np.cos(alpha) + 1j * np.sin(alpha)) - (1j*Gam*np.log(Z_mesh + a*c - a*f*1j))/(2*np.pi)


'''Mask out unwanted points inside circular region'''
phi_mask = np.ones(F_pot.shape)
for i in range(phi_mask.shape[0]):
    for j in range(phi_mask.shape[1]):
        if (np.real(Z_mesh[i][j] + a*c))**2 + (np.imag(Z_mesh[i][j])-a*f)**2 < ((a**2)*(f**2 + (1+c)**2)):
            phi_mask[i][j] = 0

F_pot = phi_mask * F_pot


contour_list = np.linspace(-15,20,30)

'''Create airfoil shape'''
theta_range = np.linspace(-np.pi,np.pi,200)
circle_z = np.zeros(len(theta_range),dtype=np.complex_)
for i in range(len(theta_range)):
    circle_z[i] = np.sqrt((a**2)*(f**2 + (1+c)**2))*np.cos(theta_range[i]) - a*c + 1j*(np.sqrt((a**2)*(f**2 + (1+c)**2))*np.sin(theta_range[i]) + a*f)
circle_zeta = circle_z + (a**2)/(circle_z)

'''Apply roatation'''
for i in range(Zeta_mesh.shape[0]):
    for j in range(Zeta_mesh.shape[1]):
        Zeta_mesh[i][j] = (np.real(Zeta_mesh[i][j])*np.cos(alpha) + np.imag(Zeta_mesh[i][j])*np.sin(alpha)) + 1j*(-np.real(Zeta_mesh[i][j])*np.sin(alpha) + np.imag(Zeta_mesh[i][j])*np.cos(alpha))

for i in range(circle_zeta.shape[0]):
    circle_zeta[i] = (np.real(circle_zeta[i])*np.cos(alpha) + np.imag(circle_zeta[i])*np.sin(alpha)) + 1j*(-np.real(circle_zeta[i])*np.sin(alpha) + np.imag(circle_zeta[i])*np.cos(alpha))

'''Plot flow'''
plt.contour(np.real(Zeta_mesh), np.imag(Zeta_mesh), np.imag(F_pot), contour_list)
plt.plot(np.real(circle_zeta),np.imag(circle_zeta),color='orange')
plt.xlim(-5*a,5*a)
plt.ylim(-3*a,3*a)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Flow Around Airfoil')
plt.savefig('Airfoil_flow.jpg',dpi=300)
plt.show()
