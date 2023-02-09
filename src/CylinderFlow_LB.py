import numpy as np
from matplotlib import pyplot

plot_every = 50
def distance(x1, y1, x2, y2):
    return np.sqrt((x2-x1)**2 + (y2-y1)**2) # Enforce Euclidian Distance Metric on the space


def main():
    Nx = 4000
    Ny = 100
    tau = .53
    Nt = 3000

    # lattice weights and speeds
    NL = 9
    cxs = np.array([0, 0, 1, 1, 1, 0, -1, -1, -1])
    cys = np.array([0, 1, 1, 0, -1, -1, -1, 0, 1])
    weights = np.array([4/9, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36])

    # I.C.'s
    F = np.ones((Ny, Nx, NL)) + 0.01*np.random.randn(Ny, Nx, NL) # np.ones is more stable than np.zeros
    F[:, :, 3] = 2.3 # feel free to mess with this value

    cylinder = np.full((Ny, Nx), False)

    for y in range(0, Ny):
        for x in range(0, Nx):
            if(distance(Nx//4, Ny//2, x, y)<13): # Set this equation to whatever you want to edit the boundary of the cylinder object
                cylinder[y][x] = True # The way we have our axes set up, y comes before x

    # Main loop
    for it in range(Nt):
        print(it)

        F[:, -1, [6, 7, 8]] = F[:, -2, [6, 7, 8]] # Apply a modified Naumann B.C. so that wall does not reflect particles
        F[:, 0, [2, 3, 4]] = F[:, 1, [2, 3, 4]] # Apply a modified Naumann B.C. so that wall does not reflect particles
        
        for i, cx, cy in zip(range(NL), cxs, cys):
            F[:, :, i] = np.roll(F[:, :, i], cx, axis = 1)
            F[:, :, 1] = np.roll(F[:, :, i], cy, axis = 0)

        bndryF = F[cylinder, :]
        bndryF = bndryF[:, [0, 5, 6, 7, 8, 1, 2, 3, 4]]

        rho = np.sum(F, 2)
        ux = np.sum(F * cxs, 2) / rho
        uy = np.sum(F * cys, 2) / rho

        F[cylinder, :] = bndryF
        ux[cylinder] = 0
        uy[cylinder] = 0

        # collision
        Feq = np.zeros(F.shape)
        for i, cx, cy, w in zip(range(NL), cxs, cys, weights):
            Feq[:, :, i] = rho * w * (1 + 3 * (cx*ux + cy*uy) + 9*(cx*ux + cy*uy)**2 / 2 - 3*(ux**2 + uy**2)/2)

        F = F + -(1.0/tau)*(F - Feq)

        if(it%plot_every == 0):
            pyplot.imshow(np.sqrt(ux**2 + uy**2))
            pyplot.pause(0.01)
            pyplot.cla()

if __name__ == "__main__":
    main()