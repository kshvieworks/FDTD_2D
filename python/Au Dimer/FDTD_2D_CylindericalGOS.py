import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import time
from matplotlib.ticker import MultipleLocator
from matplotlib.colors import SymLogNorm
import multiprocessing as mp
import pathlib
import atexit



def animate_func4(k):
    im4.set_data(np.real(Render_Ez[k, :, :]))
    t4.set_text(f"Time: {round(dt * k * 1E15, 3)} [1E-15 s]")

def animate_func5(k):
    im5.set_data(np.real(Render_Ebz[k, :, :]))
    t5.set_text(f"Time: {round(dt * k * 1E15, 3)} [1E-15 s]")

def animate_func6(k):
    im6.set_data(np.real(Render_Ez[k, :, :] - Render_Ebz[k, :, :]))
    t6.set_text(f"Time: {round(dt * k * 1E15, 3)} [1E-15 s]")

def handle_exit():
    global fig4, ax4, im4, t4, fig5, ax5, im5, t5, fig6, ax6, im6, t6
    print("Stopped")
    fig4 = plt.figure(figsize=(6.4, 4.8))
    ax4 = fig4.subplots()
    im4 = ax4.imshow(np.real(Render_Ez[0, :, :]), cmap='bwr', origin='lower', aspect=1, vmin=-1/4, vmax=1/4, alpha=1,
               extent=[-Render_Ez.shape[2] * dx * nunit / 2, Render_Ez.shape[2] * dx * nunit / 2,
                       -Render_Ez.shape[1] * dy * nunit / 2, Render_Ez.shape[1] * dy * nunit / 2])
    cbar4 = fig4.colorbar(im4, ax=ax4)
    cbar4.ax.tick_params(labelsize=15)
    ax4.set_title(f1, fontsize=20)
    ax4.set_xlabel(xlabel, fontsize=15)
    ax4.set_ylabel(ylabel, fontsize=15)
    ax4.text(0, 3 * Render_Ez.shape[1] * dy * nunit / 7, f"Wavelength: {wl * 1E9} [nm]")
    t4 = ax4.text(3, 2 * Render_Ez.shape[1] * dy * nunit / 5, f"Time: {0 * dt} 1E-15 [s]", fontsize=8)

    ax4.imshow(e > e0, alpha=0.2, cmap='Purples', origin='lower', aspect=1, vmin=0, vmax=1,
               extent=[-Render_Ez.shape[2] * dx * nunit / 2, Render_Ez.shape[2] * dx * nunit / 2,
                       -Render_Ez.shape[1] * dy * nunit / 2, Render_Ez.shape[1] * dy * nunit / 2])
    ax4.imshow(e_sigma > 0, alpha=0.1, cmap='Greys', origin='lower', aspect=1, vmin=0, vmax=1,
               extent=[-Render_Ez.shape[2] * dx * nunit / 2, Render_Ez.shape[2] * dx * nunit / 2,
                       -Render_Ez.shape[1] * dy * nunit / 2, Render_Ez.shape[1] * dy * nunit / 2])

    fig5 = plt.figure(figsize=(6.4, 4.8))
    ax5 = fig5.subplots()
    im5 = ax5.imshow(np.real(Render_Ez[0, :, :]), cmap='bwr', origin='lower', aspect=1, vmin=-1/4, vmax=1/4, alpha=1,
               extent=[-Render_Ez.shape[2] * dx * nunit / 2, Render_Ez.shape[2] * dx * nunit / 2,
                       -Render_Ez.shape[1] * dy * nunit / 2, Render_Ez.shape[1] * dy * nunit / 2])
    cbar5 = fig5.colorbar(im5, ax=ax5)
    cbar5.ax.tick_params(labelsize=15)
    ax5.set_title(f2, fontsize=20)
    ax5.set_xlabel(xlabel, fontsize=15)
    ax5.set_ylabel(ylabel, fontsize=15)
    ax5.text(0, 3 * Render_Ez.shape[1] * dy * nunit / 7, f"Wavelength: {wl * 1E9} [nm]")
    t5 = ax5.text(3, 2 * Render_Ez.shape[1] * dy * nunit / 5, f"Time: {0 * dt} 1E-15 [s]", fontsize=8)

    ax5.imshow(e > e0, alpha=0.2, cmap='Purples', origin='lower', aspect=1, vmin=0, vmax=1,
               extent=[-Render_Ez.shape[2] * dx * nunit / 2, Render_Ez.shape[2] * dx * nunit / 2,
                       -Render_Ez.shape[1] * dy * nunit / 2, Render_Ez.shape[1] * dy * nunit / 2])
    ax5.imshow(e_sigma > 0, alpha=0.1, cmap='Greys', origin='lower', aspect=1, vmin=0, vmax=1,
               extent=[-Render_Ez.shape[2] * dx * nunit / 2, Render_Ez.shape[2] * dx * nunit / 2,
                       -Render_Ez.shape[1] * dy * nunit / 2, Render_Ez.shape[1] * dy * nunit / 2])

    fig6 = plt.figure(figsize=(6.4, 4.8))
    ax6 = fig6.subplots()
    im6 = ax6.imshow(np.real(Render_Ez[0, :, :]), cmap='bwr', origin='lower', aspect=1, vmin=-1/4, vmax=1/4, alpha=1,
               extent=[-Render_Ez.shape[2] * dx * nunit / 2, Render_Ez.shape[2] * dx * nunit / 2,
                       -Render_Ez.shape[1] * dy * nunit / 2, Render_Ez.shape[1] * dy * nunit / 2])
    cbar6 = fig6.colorbar(im6, ax=ax6)
    cbar6.ax.tick_params(labelsize=15)
    ax6.set_title(f3, fontsize=20)
    ax6.set_xlabel(xlabel, fontsize=15)
    ax6.set_ylabel(ylabel, fontsize=15)
    ax6.text(0, 3 * Render_Ez.shape[1] * dy * nunit / 7, f"Wavelength: {wl * 1E9} [nm]")
    t6 = ax6.text(3, 2 * Render_Ez.shape[1] * dy * nunit / 5, f"Time: {0 * dt} 1E-15 [s]", fontsize=8)

    ax6.imshow(e > e0, alpha=0.2, cmap='Purples', origin='lower', aspect=1, vmin=0, vmax=1,
               extent=[-Render_Ez.shape[2] * dx * nunit / 2, Render_Ez.shape[2] * dx * nunit / 2,
                       -Render_Ez.shape[1] * dy * nunit / 2, Render_Ez.shape[1] * dy * nunit / 2])
    ax6.imshow(e_sigma > 0, alpha=0.1, cmap='Greys', origin='lower', aspect=1, vmin=0, vmax=1,
               extent=[-Render_Ez.shape[2] * dx * nunit / 2, Render_Ez.shape[2] * dx * nunit / 2,
                       -Render_Ez.shape[1] * dy * nunit / 2, Render_Ez.shape[1] * dy * nunit / 2])

    anim4 = FuncAnimation(fig4, animate_func4, frames = t, interval = 100)
    anim4.save('Total Field.gif', fps=10)

    anim5 = FuncAnimation(fig5, animate_func5, frames = t, interval = 100)
    anim5.save('Background Field.gif', fps=10)

    anim6 = FuncAnimation(fig6, animate_func6, frames = t, interval = 100)
    anim6.save('Scattered Field.gif', fps=10)

    print("Rendering Finished")


atexit.register(handle_exit)



##########################################FDTD Settings###################################################

# Render_Ez = np.empty((0,0,0))
# Render_Ebz = np.empty((0,0,0))

# Define Figure Settings
fs = (12.8, 9.6)
funit = '$\mu$m'
nunit = 1E6
f1 = f"Total Field, |E$_T$|"
f2 = f"Background Field, |E$_B$|"
f3 = f"Scattering Field, |E$_S$|"
xlabel = f"x Location [{funit}]"
ylabel = f"y Location [{funit}]"
fd = pathlib.Path(__file__).parent.resolve()

# Define Constants
e0 = 8.854E-12 # [F/m]
m0 = 4*np.pi*1E-7 # [H/m]
c = 1/np.sqrt(e0*m0) # [m/s]

[e_r, e_sigma] = [1, 0]
[m_r, m_sigma] = [1, 0]
wl = 545E-9 # [m]

# Define Spatial System
[Lx, Ly] = [2E-5, 2E-5] # [m]
[Nx, Ny] = [501, 501] # Number of Mesh Point
[x, y] = [np.linspace(0, Lx, Nx), np.linspace(0, Ly, Ny)]
[dx, dy] = [x[1], y[1]]

# Define Temporal System
CFL = 1
dt = CFL/np.sqrt(np.power(c/dx, 2) + np.power(c/dy, 2))
v = c/wl

# Sigma for PML Structure
PML_Nx = 100
PML_Ny = 100

PML_Lx = PML_Nx * dx
sigma_x = np.ones((Ny, PML_Nx)) * (e0 / (2*dt)) * np.power((dx/PML_Lx) * np.linspace(1, PML_Nx, PML_Nx)[np.newaxis], 3)

PML_Ly = PML_Ny * dy
sigma_y = np.ones((PML_Ny, Nx)) * (e0 / (2*dt)) * np.power((dy/PML_Ly) * np.linspace(1, PML_Ny, PML_Ny)[np.newaxis].T, 3)

sigma_xy = np.ones((PML_Ny, PML_Nx)) * (e0 / (2*dt)) * (np.power((dx/PML_Lx) * np.linspace(1, PML_Nx, PML_Nx), 3))
sigma_xy = sigma_xy + (e0 / (2*dt))*np.power((dy/PML_Ly) * np.linspace(1, PML_Ny, PML_Ny)[np.newaxis].T, 3)

sigma_y = np.append(np.flip(sigma_xy, axis=1), sigma_y, axis=1)
sigma_y = np.append(sigma_y, sigma_xy, axis=1)

e_sigma = np.ones((Ny, Nx), dtype=complex) * e_sigma
e_sigma = np.append(e_sigma, sigma_x, axis=1)
e_sigma = np.append(np.flip(sigma_x, axis=1), e_sigma, axis=1)
e_sigma = np.append(np.flip(sigma_y, axis=0), e_sigma, axis=0)
e_sigma = np.append(e_sigma, sigma_y, axis=0)

e = np.ones((Ny+2*PML_Ny, Nx+2*PML_Nx), dtype=complex) * e_r * e0
m = np.ones((Ny+2*PML_Ny, Nx+2*PML_Nx), dtype=complex) * m_r * m0

m_sigma = e_sigma * m / e

########################################Structure Settings###################################################

n = 2.2
k = 0
material_e = np.power(n, 2) - np.power(k, 2) + 1j*n*k

cr = 1.75E-6
[cx1, cy1] = [0, 0]
x_grid = np.arange(0, Nx + 2*PML_Nx)
y_grid = np.arange(0, Ny + 2*PML_Ny)

eb = e.copy()
mb = m.copy()
e[((x_grid[np.newaxis, :] - (e.shape[0]/2 + cx1/dx))**2) + ((y_grid[:, np.newaxis] - (e.shape[1]/2 + cy1/dy))**2) < (cr/dx)**2] = material_e * e0

#############################################################################################################
# Coefficients
A = (2*m-m_sigma*dt)/(2*m+m_sigma*dt)
B = (2*dt)/(2*m+m_sigma*dt)
C = (2*e-e_sigma*dt)/(2*e+e_sigma*dt)
D = (2*dt)/(2*e+e_sigma*dt)

# Background Coefficients
bA = (2*mb-m_sigma*dt)/(2*mb+m_sigma*dt)
bB = (2*dt)/(2*mb+m_sigma*dt)
bC = (2*eb-e_sigma*dt)/(2*eb+e_sigma*dt)
bD = (2*dt)/(2*eb+e_sigma*dt)

# Initialize Field
Hx = np.zeros((Ny+2*PML_Ny, Nx+2*PML_Nx), dtype=complex)
Hy = np.zeros((Ny+2*PML_Ny, Nx+2*PML_Nx), dtype=complex)
Ez = np.zeros((Ny+2*PML_Ny, Nx+2*PML_Nx), dtype=complex)

# Background Initial Field
Hbx = np.zeros((Ny+2*PML_Ny, Nx+2*PML_Nx), dtype=complex)
Hby = np.zeros((Ny+2*PML_Ny, Nx+2*PML_Nx), dtype=complex)
Ebz = np.zeros((Ny+2*PML_Ny, Nx+2*PML_Nx), dtype=complex)

########################################### Open Figure Window ################################################
# Total Field
fig1 = plt.figure(figsize=fs)
ax1 = fig1.subplots()

cx1 = ax1.imshow(np.real(Ez), cmap='bwr', origin='lower', aspect=1, vmin=-1/4, vmax=1/4, alpha=1,
               extent=[-Ez.shape[1] * dx * nunit / 2, Ez.shape[1] * dx * nunit / 2,
                       -Ez.shape[0] * dy * nunit / 2, Ez.shape[0] * dy * nunit / 2])
cbar1 = fig1.colorbar(cx1, ax=ax1)
cbar1.ax.tick_params(labelsize=15)
ax1.set_title(f1, fontsize=30)
ax1.set_xlabel(xlabel, fontsize=30)
ax1.set_ylabel(ylabel, fontsize=30)
ax1.text(0, 3*Ez.shape[0] * dy * nunit / 7, f"Wavelength: {wl*1E9} [nm]")
t1 = ax1.text(0, 2*Ez.shape[0] * dy * nunit / 5, f"Time: {0*dt} 1E-15 [s]")

ax1.imshow(e > e0, alpha = 0.2, cmap='Purples', origin='lower', aspect=1, vmin=0, vmax=1,
           extent=[-Ez.shape[1] * dx * nunit / 2, Ez.shape[1] * dx * nunit / 2,
                   -Ez.shape[0] * dy * nunit / 2, Ez.shape[0] * dy * nunit / 2])

ax1.imshow(e_sigma > 0, alpha=0.1, cmap='Greys', origin='lower', aspect=1, vmin=0, vmax=1,
           extent=[-Ez.shape[1] * dx * nunit / 2, Ez.shape[1] * dx * nunit / 2,
                   -Ez.shape[0] * dy * nunit / 2, Ez.shape[0] * dy * nunit / 2])
plt.show(block=False)

# Background Field
fig2, ax2 = plt.subplots(figsize=fs)
cx2 = ax2.imshow(np.real(Ebz), cmap='bwr', origin='lower', aspect=1, vmin=-1/4, vmax=1/4,
               extent=[-Ez.shape[1] * dx * nunit / 2, Ez.shape[1] * dx * nunit / 2,
                       -Ez.shape[0] * dy * nunit / 2, Ez.shape[0] * dy * nunit / 2])
cbar2 = fig2.colorbar(cx2, ax=ax2)
cbar2.ax.tick_params(labelsize=15)
ax2.set_title(f2, fontsize=30)
ax2.set_xlabel(xlabel, fontsize=30)
ax2.set_ylabel(ylabel, fontsize=30)
ax2.text(0, 3*Ez.shape[0] * dy * nunit / 7, f"Wavelength: {wl*1E9}[nm]")
t2 = ax2.text(0, 2*Ez.shape[0] * dy * nunit / 5, f"Time: {0*dt} 1E-15 [s]")


ax2.imshow(e > e0, alpha = 0.2, cmap='Purples', origin='lower', aspect=1, vmin=0, vmax=1,
           extent=[-Ez.shape[1] * dx * nunit / 2, Ez.shape[1] * dx * nunit / 2,
                   -Ez.shape[0] * dy * nunit / 2, Ez.shape[0] * dy * nunit / 2])

ax2.imshow(e_sigma > 0, alpha=0.1, cmap='Greys', origin='lower', aspect=1, vmin=0, vmax=1,
           extent=[-Ez.shape[1] * dx * nunit / 2, Ez.shape[1] * dx * nunit / 2,
                   -Ez.shape[0] * dy * nunit / 2, Ez.shape[0] * dy * nunit / 2])

plt.show(block=False)

# Scattering Field
fig3, ax3 = plt.subplots(figsize=fs)
cx3 = ax3.imshow(np.real(Ez - Ebz), cmap='bwr', origin='lower', aspect=1, vmin=-1/4, vmax=1/4,
               extent=[-Ez.shape[1] * dx * nunit / 2, Ez.shape[1] * dx * nunit / 2,
                       -Ez.shape[0] * dy * nunit / 2, Ez.shape[0] * dy * nunit / 2])
cbar3 = fig1.colorbar(cx3, ax=ax3)
cbar3.ax.tick_params(labelsize=15)
ax3.set_title(f3, fontsize=30)
ax3.set_xlabel(xlabel, fontsize=30)
ax3.set_ylabel(ylabel, fontsize=30)
ax3.text(0, 3*Ez.shape[0] * dy * nunit / 7, f"Wavelength: {wl*1E9}[nm]")
t3 = ax3.text(0, 2*Ez.shape[0] * dy * nunit / 5, f"Time: {0*dt} 1E-15 [s]")

ax3.imshow(e > e0, alpha = 0.2, cmap='Purples', origin='lower', aspect=1, vmin=0, vmax=1,
           extent=[-Ez.shape[1] * dx * nunit / 2, Ez.shape[1] * dx * nunit / 2,
                   -Ez.shape[0] * dy * nunit / 2, Ez.shape[0] * dy * nunit / 2])

ax3.imshow(e_sigma > 0, alpha=0.1, cmap='Greys', origin='lower', aspect=1, vmin=0, vmax=1,
           extent=[-Ez.shape[1] * dx * nunit / 2, Ez.shape[1] * dx * nunit / 2,
                   -Ez.shape[0] * dy * nunit / 2, Ez.shape[0] * dy * nunit / 2])

plt.show(block=False)

t = 0

Render_Ez = np.zeros((1, Ez.shape[0], Ez.shape[1]), dtype=complex)
Render_Ebz = np.zeros((1, Ez.shape[0], Ez.shape[1]), dtype=complex)


# Update Equations
while True:

# Total Field
    Hx[:-1, :-1] = A[:-1, :-1] * Hx[:-1, :-1] - B[:-1, :-1] * (Ez[1:, :-1] - Ez[:-1, :-1]) / dy
    Hy[:-1, :-1] = A[:-1, :-1] * Hy[:-1, :-1] + B[:-1, :-1] * (Ez[:-1, 1:] - Ez[:-1, :-1]) / dx
    Ez[1:-1, 1:-1] = C[1:-1, 1:-1] * Ez[1:-1, 1:-1] + D[1:-1, 1:-1] * ((Hy[1:-1, 1:-1] - Hy[1:-1, :-2])/dx - (Hx[1:-1, 1:-1] - Hx[:-2, 1:-1])/dy)
    Ez[-1, :] = Ez[-2, :]
    Ez[:, -1] = Ez[:, -2]

# Background Field
    Hbx[:-1, :-1] = bA[:-1, :-1] * Hbx[:-1, :-1] - bB[:-1, :-1] * (Ebz[1:, :-1] - Ebz[:-1, :-1]) / dy
    Hby[:-1, :-1] = bA[:-1, :-1] * Hby[:-1, :-1] + bB[:-1, :-1] * (Ebz[:-1, 1:] - Ebz[:-1, :-1]) / dx
    Ebz[1:-1, 1:-1] = bC[1:-1, 1:-1] * Ebz[1:-1, 1:-1] + bD[1:-1, 1:-1] * (
                      (Hby[1:-1, 1:-1] - Hby[1:-1, :-2]) / dx - (Hbx[1:-1, 1:-1] - Hbx[:-2, 1:-1]) / dy)
    Ebz[-1, :] = Ebz[-2, :]
    Ebz[:, -1] = Ebz[:, -2]

# Define Source
    Ebz[int(Ebz.shape[0]/2), int(Ebz.shape[1]/4)] = np.sin(2*np.pi*v*dt*t)
    Ez[int(Ez.shape[0] / 2), int(Ez.shape[1] / 4)] = np.sin(2 * np.pi * v * dt * t)

    t = t + 1

    cx1.set_data(np.real(Ez))
    cx2.set_data(np.real(Ebz))
    cx3.set_data(np.real(Ez - Ebz))
    t1.set_text(f"Time: {round(dt*t*1E15, 3)} [1E-15 s]")
    t2.set_text(f"Time: {round(dt*t*1E15, 3)} [1E-15 s]")
    t3.set_text(f"Time: {round(dt*t*1E15, 3)} [1E-15 s]")

    plt.draw()
    plt.pause(0.01)



    Render_Ez = np.append(Render_Ez, Ez[np.newaxis], axis=0)
    Render_Ebz = np.append(Render_Ebz, Ebz[np.newaxis], axis=0)
