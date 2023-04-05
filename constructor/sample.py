import numpy as np
from bondset import Bondtype
from dpd_files import print_bonds, print_coord, print_fixed
from mol import Dendron
from periodic_box import Box
from scipy import interpolate

if __name__ == '__main__':
    main_chain = 50
    n = 3
    g = 2
    m = 2
    length_bond = (1/3)**(1/3)
    N = (2**(g+1) - 1) * n
    lmbd = [1.571, 1.846, 2.379, 3.164,	4.294] #for intepolation
    x_g = [0, 1, 2, 3, 4] #for intepolation
    lmbd_g = interpolate.InterpolatedUnivariateSpline(x_g, lmbd)
    #h_xy = (g+1) * n * length_bond * 2 + 1
    h_xy = 2 * (4*N**3/(3*np.pi*lmbd_g(g)**2))**0.25 * length_bond + 2.0
    
    box = Box(x=h_xy, y=h_xy, z=main_chain*length_bond)
    x = np.array([])
    y = np.array([])
    z = np.array([])
    b_type = np.array([], dtype=int)
    bonds: Bondtype = []
    fixed_list = []
    dendron = Dendron(n=n, g=g)
    k = 1
    for i in range(main_chain):
        if i%m == 0:
            fixed_list.append(k)
            b_type = np.hstack(
                [b_type, np.array([1] + [2] * (dendron.num_beads - 1))])
            z0 = 0.5 * (-box.z + length_bond) + i * length_bond
            d = {0: (0., 0., z0)}
            (x_t, y_t, z_t) = dendron.get_coords(box=box, fixed_coords=d)
            x = np.hstack([x, x_t])
            y = np.hstack([y, y_t])
            z = np.hstack([z, z_t])
            b = [(bond[0] + k, bond[1] + k) for bond in dendron.bonds]
            bonds += b
            k += dendron.num_beads
        else:
            fixed_list.append(k)
            b_type = np.hstack([b_type, np.array([1])])
            z0 = 0.5 * (-box.z + length_bond) + i * length_bond
            x = np.hstack([x, np.array([0.])])
            y = np.hstack([y, np.array([0.])])
            z = np.hstack([z, np.array([z0])])
            k += 1
    n_solvent = int(box.volume * 3 - k)
    x = np.hstack([x, np.random.uniform(-box.x*0.5, box.x*0.5, n_solvent)])
    y = np.hstack([y, np.random.uniform(-box.y*0.5, box.y*0.5, n_solvent)])
    z = np.hstack([z, np.random.uniform(-box.z*0.5, box.z*0.5, n_solvent)])
    b_type = np.hstack([b_type, np.array([3] * n_solvent)])
    print_coord(box, x, y, z, b_type)
    print_bonds(box=box, num_atoms=len(x), bonds=bonds)
    print_fixed(fixed_list)
