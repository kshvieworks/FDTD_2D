import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.ticker import MultipleLocator
from matplotlib.colors import SymLogNorm
import time
import multiprocessing as mp
import pathlib
import threading


'''''''''''''''''''''''''''''''''''''''''Executable Example'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
# import numpy as np
# import multiprocessing as mp
# import testFDTD as TM
# l = TM.FDTD2D(1, 4.48-2.53j, 532E-9, 40E-9, 2E-7, 0, 150, np.arange(200,600))
# l.DefinePML()
# l.DefineSphere()
# if __name__ == "__main__":
#     p = mp.Pool(12)
#     p.map(l.RunFDTD())
#     p.close()
#     p.join()
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

'''''''''''''''''''''''''''''''''''''''''Physical Constants'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
e0 = 8.854E-12 # [F/m]
m0 = 4*np.pi*1E-7 # [H/m]
c = 1/np.sqrt(e0*m0) # [m/s]
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

'''''''''''''''''''''''''''''''''''''''''Define Spatial & Temporal System'''''''''''''''''''''''''''''''''''''''''''''

[Lx, Ly] = [10E-6, 10E-6] # [m]
[Nx, Ny] = [401, 401] # Number of Mesh Point
[x, y] = [np.linspace(0, Lx, Nx), np.linspace(0, Ly, Ny)]
[dx, dy] = [x[1], y[1]]

CFL = 1
dt = CFL/np.sqrt(np.power(c/dx, 2) + np.power(c/dy, 2))

PML_Nx = 100
PML_Ny = 100
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

'''''''''''''''''''''''''''''''''''''''''Figure Settings'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
FS = (12.8, 9.6)
FUNIT = '$\mu$m'
NUNIT = 1E6
F1 = f"Total Field, |E$_T$|"
F2 = f"Background Field, E$_B$"
F3 = f"Scattering Field, |E$_S$|"
F4 = f"Background Field in LogScale, E$_B$"
XLABEL = f"x Location [{FUNIT}]"
YLABEL = f"y Location [{FUNIT}]"
FD = pathlib.Path(__file__).parent.resolve()
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

class FDTD2D:
    def __init__(self, e_m, e_p, wl, a, x, y, srcx, srcy):

        e_m = np.real(e_m) - 1j*np.imag(e_m)
        e_p = np.real(e_p) - 1j*np.imag(e_p)

        [self.e_m, self.e_p, self.wl, self.a, self.x, self.y, self.srcx, self.srcy] = [e_m, e_p, wl, a, x, y, srcx, srcy]

        self.c = []
        # for k in range(self.x.__len__()):
        #     self.c.append((self.x[k], self.y[k]))

        self.n_m = np.sqrt((np.abs(self.e_m)+np.real(self.e_m))/2) + 1j*np.sqrt((np.abs(self.e_m)-np.real(self.e_m))/2)
        self.n_p = np.sqrt((np.abs(self.e_p)+np.real(self.e_p))/2) + 1j*np.sqrt((np.abs(self.e_p)-np.real(self.e_p))/2)

        self.k0 = 2 * np.pi / self.wl
        self.k_m = self.n_m * self.k0
        self.k_p = self.n_p * self.k0
        self.v = c / self.wl
        self.w = 2 * np.pi * self.v

        self.mu_m = 1 # Medium
        self.mu_p = 1 # Particle

        self.sigma_e_m = 0
        self.sigma_e_p = 0
        self.sigma_mu_m = 0
        self.sigma_mu_p = 0

        self.e = e0 * np.real(self.e_m) * np.ones((Ny, Nx), dtype=complex)
        self.sigma_e = self.sigma_e_m * np.ones((Ny, Nx), dtype=complex)
        self.mu = m0 * np.real(self.mu_m) * np.ones((Ny, Nx), dtype=complex)
        self.sigma_mu = self.sigma_mu_m * np.ones((Ny, Nx), dtype=complex)

    def DefinePML(self):

        if self.sigma_e.shape[1] > Nx:
            return

        PML_Lx = PML_Nx * dx
        PML_Ly = PML_Ny * dy

        sigma_x = np.ones((Ny, PML_Nx)) * (e0/(2*dt)) * np.power((dx/PML_Lx) * np.linspace(1, PML_Nx, PML_Nx)[np.newaxis], 3)
        sigma_y = np.ones((PML_Ny, Nx)) * (e0/(2*dt)) * np.power((dy/PML_Ly) * np.linspace(1, PML_Ny, PML_Ny)[np.newaxis].T, 3)
        sigma_xy = np.ones((PML_Ny, PML_Nx)) * (e0 / (2*dt)) * (np.power((dx/PML_Lx) * np.linspace(1, PML_Nx, PML_Nx), 3))
        sigma_xy = sigma_xy + (e0 / (2*dt))*np.power((dy/PML_Ly) * np.linspace(1, PML_Ny, PML_Ny)[np.newaxis].T, 3)

        sigma_y = np.append(np.flip(sigma_xy, axis=1), sigma_y, axis=1)
        sigma_y = np.append(sigma_y, sigma_xy, axis=1)

        x_template = np.ones((sigma_x.shape[0], sigma_x.shape[1]), dtype=complex)
        y_template = np.ones((sigma_y.shape[0], sigma_y.shape[1]), dtype=complex)

        self.sigma_e = self.AppendStructure(self.sigma_e, sigma_x, sigma_y)
        self.sigma_mu = self.AppendStructure(self.sigma_mu, sigma_x*self.mu_m*m0/(self.e_m*e0), sigma_y*self.mu_m*m0/(self.e_m*e0))
        self.e = self.AppendStructure(self.e, e0*self.e_m*x_template, e0*self.e_m*y_template)
        self.mu = self.AppendStructure(self.mu, m0*self.mu_m*x_template, m0*self.mu_m*y_template)

        self.e_bg = self.e.copy()
        self.mu_bg = self.mu.copy()
        self.sigma_e_bg = self.sigma_e.copy()
        self.sigma_mu_bg = self.sigma_mu.copy()


    def AppendStructure(self, data, x, y):

        data = np.append(data, x, axis=1)
        data = np.append(np.flip(x, axis=1), data, axis=1)
        data = np.append(np.flip(y, axis=0), data, axis=0)
        data = np.append(data, y, axis=0)

        return data

    def DefineSphere(self):

        x_grid = np.arange(0, self.e.shape[1])
        y_grid = np.arange(0, self.e.shape[0])

        for (cx, cy) in self.c:
            [cx_mesh, cy_mesh, a_mesh] = [cx/dx, cy/dy, self.a/((dx+dy)/2)]

            x_logic = x_grid[np.newaxis, :] - (x_grid.max()/2 + cx_mesh)
            y_logic = y_grid[:, np.newaxis] - (y_grid.max()/2 + cy_mesh)

            self.e = self.SphereParams(self.e, x_logic, y_logic, a_mesh, e0 * np.real(self.e_p))
            self.sigma_e = self.SphereParams(self.sigma_e, x_logic, y_logic, a_mesh, e0 * np.imag(self.e_p) * self.w)
            self.mu = self.SphereParams(self.mu, x_logic, y_logic, a_mesh, m0 * np.real(self.mu_p))
            self.sigma_mu = self.SphereParams(self.sigma_mu, x_logic, y_logic, a_mesh, m0 * np.imag(self.mu_p) * self.w)

    def SphereParams(self, data, x, y, a, sdata):

        data[(np.power(x, 2) + np.power(y, 2) <= np.power(a, 2))] = sdata

        return data

    def RunFDTD(self):

        [A, B, C, D] = self.UpdateCoeff(self.e, self.mu, self.sigma_e, self.sigma_mu)
        [bA, bB, bC, bD] = self.UpdateCoeff(self.e_bg, self.mu_bg, self.sigma_e_bg, self.sigma_mu_bg)
        self.InitializeField(self.e.shape[1], self.e.shape[0])

        t = 0
        FigureParam = self.FigureOpen(self.Hz, self.Hzb, True, False, False, True)

        while True:

            # print(self.Hz[self.srcy, self.srcx])
            # Total Field
            self.Ex[:-1, :-1] = A[:-1, :-1] * self.Ex[:-1, :-1] + B[:-1, :-1] * (self.Hz[1:, :-1] - self.Hz[:-1, :-1]) / dy
            self.Ey[:-1, :-1] = A[:-1, :-1] * self.Ey[:-1, :-1] - B[:-1, :-1] * (self.Hz[:-1, 1:] - self.Hz[:-1, :-1]) / dx
            self.Hz[1:-1, 1:-1] = C[1:-1, 1:-1] * self.Hz[1:-1, 1:-1] + \
                                  D[1:-1, 1:-1] * (-((self.Ey[1:-1, 1:-1] - self.Ey[1:-1, :-2]) / dx) + ((self.Ex[1:-1, 1:-1] - self.Ex[:-2, 1:-1]) / dy))
            self.Hz[-1, :] = self.Hz[-2, :]
            self.Hz[:, -1] = self.Hz[:, -2]

            # Background Field
            self.Exb[:-1, :-1] = bA[:-1, :-1] * self.Exb[:-1, :-1] + bB[:-1, :-1] * (self.Hzb[1:, :-1] - self.Hzb[:-1, :-1]) / dy
            self.Eyb[:-1, :-1] = bA[:-1, :-1] * self.Eyb[:-1, :-1] - bB[:-1, :-1] * (self.Hzb[:-1, 1:] - self.Hzb[:-1, :-1]) / dx
            self.Hzb[1:-1, 1:-1] = bC[1:-1, 1:-1] * self.Hzb[1:-1, 1:-1] + \
                                   bD[1:-1, 1:-1] * (-((self.Eyb[1:-1, 1:-1] - self.Eyb[1:-1, :-2]) / dx) + ((self.Exb[1:-1, 1:-1] - self.Exb[:-2, 1:-1]) / dy))
            self.Hzb[-1, :] = self.Hzb[-2, :]
            self.Hzb[:, -1] = self.Hzb[:, -2]

            self.Source(self.Hz, self.Hzb, self.srcx, self.srcy, t, 'Gaussian', 1, 36.15)

            t = t+1

            self.FigureUpdate(FigureParam, t)

    def UpdateCoeff(self, e, mu, sigma_e, sigma_mu):

        A = (2*e - sigma_e*dt) / (2*e + sigma_e*dt)
        B = 2*dt / (2*e + sigma_e*dt)
        C = (2*mu - sigma_mu*dt) / (2*mu + sigma_mu*dt)
        D = 2*dt / (2*mu + sigma_mu*dt)

        return [A, B, C, D]

    def InitializeField(self, lx, ly):

        self.Ex = np.zeros((ly, lx), dtype=complex)
        self.Ey = np.zeros((ly, lx), dtype=complex)
        self.Hz = np.zeros((ly, lx), dtype=complex)

        self.Exb = np.zeros((ly, lx), dtype=complex)
        self.Eyb = np.zeros((ly, lx), dtype=complex)
        self.Hzb = np.zeros((ly, lx), dtype=complex)

    def Source(self, E, Eb, x, y, t, sourcetype, z=1, theta=90):

        if sourcetype == 'Point':
            E[y, x] = np.sin(self.w*dt*t)
            Eb[y, x] = np.sin(self.w*dt*t)

        if sourcetype == 'Plane':
            if 2*self.w*dt*t<np.pi:
                E[y, x] = np.sin(2*self.w*dt*t)
                Eb[y, x] = np.sin(2*self.w*dt*t)
            else:
                E[y, x] = 0
                Eb[y, x] = 0

        if sourcetype == 'Gaussian':
                E[y, x] = np.exp(-np.power((PML_Ny+y-y[int(y.size/2)])/(z*np.tan(np.deg2rad(theta))), 2))*np.sin(2*self.w*dt*t)
                Eb[y, x] = np.exp(-np.power((PML_Ny+y-y[int(y.size/2)])/(z*np.tan(np.deg2rad(theta))), 2))*np.sin(2*self.w*dt*t)


    def FigureOpen(self, E, Eb, f1=True, f2=True, f3=True, f4=True):

        returnvalue = []

        if f1:
            fig1 = plt.figure(figsize=FS)
            ax1 = fig1.subplots()
            self.cx1 = ax1.imshow(np.abs(E), cmap='bwr', origin='lower', aspect=1, vmin=-1, vmax=1, alpha=1,
                             extent=[-E.shape[1] * dx * NUNIT / 2, E.shape[1] * dx * NUNIT / 2,
                                     -E.shape[0] * dy * NUNIT / 2, E.shape[0] * dy * NUNIT / 2])
            self.t1 = ax1.text(0, 2 * E.shape[0] * dy * NUNIT / 5, f"Time: {0 * dt} 1E-15 [s]")

            self.FigureSettings(E, fig1, ax1, self.cx1, F1)

            returnvalue.append((self.cx1, self.t1, E))

        if f2:
            fig2 = plt.figure(figsize=FS)
            ax2 = fig2.subplots()
            self.cx2 = ax2.imshow(np.real(Eb), cmap='bwr', origin='lower', aspect=1, vmin=-1, vmax=1, alpha=1,
                             extent=[-E.shape[1] * dx * NUNIT / 2, E.shape[1] * dx * NUNIT / 2,
                                     -E.shape[0] * dy * NUNIT / 2, E.shape[0] * dy * NUNIT / 2])
            self.t2 = ax2.text(0, 2 * E.shape[0] * dy * NUNIT / 5, f"Time: {0 * dt} 1E-15 [s]")

            self.FigureSettings(Eb, fig2, ax2, self.cx2, F2)

            returnvalue.append((self.cx2, self.t2, Eb))


        if f3:
            fig3 = plt.figure(figsize=FS)
            ax3 = fig3.subplots()
            self.cx3 = ax3.imshow(np.abs(E - Eb), cmap='bwr', origin='lower', aspect=1, vmin=-1E-3, vmax=1E-3, alpha=1,
                                  extent=[-E.shape[1] * dx * NUNIT / 2, E.shape[1] * dx * NUNIT / 2,
                                          -E.shape[0] * dy * NUNIT / 2, E.shape[0] * dy * NUNIT / 2])
            self.t3 = ax3.text(0, 2 * E.shape[0] * dy * NUNIT / 5, f"Time: {0 * dt} 1E-15 [s]")

            self.FigureSettings(E, fig3, ax3, self.cx3, F3)

            returnvalue.append((self.cx3, self.t3, E - Eb))

        if f4:
            fig4 = plt.figure(figsize=FS)
            ax4 = fig4.subplots()
            self.cx4 = ax4.imshow(np.real(Eb), cmap='bwr', origin='lower', aspect=1, alpha=1,
                                  extent=[-E.shape[1] * dx * NUNIT / 2, E.shape[1] * dx * NUNIT / 2,
                                  -E.shape[0] * dy * NUNIT / 2, E.shape[0] * dy * NUNIT / 2],
                                  norm = SymLogNorm(linthresh=np.power(10, 30), vmin=-1, vmax=1, base=10))
            self.t4 = ax4.text(0, 2 * E.shape[0] * dy * NUNIT / 5, f"Time: {0 * dt} 1E-15 [s]")

            self.FigureSettings(E, fig4, ax4, self.cx4, F4)

            returnvalue.append((self.cx4, self.t4, Eb))

        return returnvalue


    def FigureSettings(self, E, fig, ax, cx, ftitle):
        cbar = fig.colorbar(cx, ax=ax)
        cbar.ax.tick_params(labelsize=15)
        ax.set_title(ftitle, fontsize=30)
        ax.set_xlabel(XLABEL, fontsize=30)
        ax.set_ylabel(YLABEL, fontsize=30)
        ax.text(0, 3 * E.shape[0] * dy * NUNIT / 7, f"Wavelength: {self.wl * 1E9} [nm]")

        ax.imshow(self.e == e0 * np.real(self.e_p), alpha=0.2, cmap='Purples', origin='lower', aspect=1, vmin=0, vmax=1,
                  extent=[-E.shape[1] * dx * NUNIT / 2, E.shape[1] * dx * NUNIT / 2,
                          -E.shape[0] * dy * NUNIT / 2, E.shape[0] * dy * NUNIT / 2])

        ax.imshow(self.sigma_e > 0, alpha=0.1, cmap='Greys', origin='lower', aspect=1, vmin=0, vmax=1,
                  extent=[-E.shape[1] * dx * NUNIT / 2, E.shape[1] * dx * NUNIT / 2,
                          -E.shape[0] * dy * NUNIT / 2, E.shape[0] * dy * NUNIT / 2])
        plt.show(block=False)

    def FigureUpdate(self, fp, t):

        for (cx, tx, E) in fp:
            cx.set_data(np.real(E))
            tx.set_text(f"Time: {round(dt*t*1E15, 3)} [1E-15 s]")

        plt.draw()
        plt.pause(0.0001)



# if __name__ == "__main__":
#     l = FDTD2D(1, -5.0687 - 1j * 2.7106, 532E-9, 15E-9, 0, 0, 45, 90)
#     l.DefinePML()
#     l.DefineSphere()
#     p = mp.Process(target=l.RunFDTD)
#     p.start()
#     p.join()
