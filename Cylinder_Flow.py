import numpy as np
import matplotlib.pyplot as plt
plt.style.use('dark_background')

U = 2.0
alpha = 0
a = 2.0
Gamma_list = [0, -10, -30, -50]  


fig, axs = plt.subplots(2, 2, figsize=(12, 10))
axs = axs.flatten()

for i, Gam in enumerate(Gamma_list):
    x_array = np.linspace(-10*a, 10*a, 1000)
    y_array = np.linspace(-5*a, 5*a, 500)

    X_mesh, Y_mesh = np.meshgrid(x_array, 1j * y_array)
    Z_mesh = X_mesh + Y_mesh

    F_pot = U * Z_mesh * (np.cos(alpha) - 1j * np.sin(alpha)) + U * (a**2 / Z_mesh) * (np.cos(alpha) + 1j * np.sin(alpha)) - (1j*Gam*np.log(Z_mesh))/(2*np.pi)

    phi_mask = np.ones(F_pot.shape)
    for j in range(phi_mask.shape[0]):
        for k in range(phi_mask.shape[1]):
            if (np.real(Z_mesh[j][k])**2 + np.imag(Z_mesh[j][k])**2) < a**2:
                phi_mask[j][k] = 0

    F_pot = phi_mask * F_pot

    contour_list = np.linspace(-9*a, 20*a, 30)

    x_circle = np.linspace(-a, a, 300)
    y_circle_p = np.sqrt(a**2 - x_circle**2)
    y_circle_n = -np.sqrt(a**2 - x_circle**2)


    axs[i].contour(np.real(Z_mesh), np.imag(Z_mesh), np.imag(F_pot), contour_list)
    axs[i].plot(x_circle, y_circle_p, color='orange')
    axs[i].plot(x_circle, y_circle_n, color='orange')
    axs[i].set_xlabel('x')
    axs[i].set_ylabel('y')
    axs[i].set_title('Circulation: ' + str(Gam))


plt.tight_layout()
plt.show()
