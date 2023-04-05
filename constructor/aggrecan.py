import numpy as np
from bondset import Bondtype
from dpd_files import print_bonds, print_coord, print_fixed
from mol import Aggrecan
from periodic_box import Box


if __name__ == '__main__':
    N1 = 30
    N2 = 30
    n1 = 0
    n2 = 20
    m1 = 1
    m2 = 1
    num_chain = 64
    length_bond = (1/3)**(1/3)

    box = Box(x=39.222, y=39.222, z=39.222)
    x = np.array([])
    y = np.array([])
    z = np.array([])
    fixed_list = []
    b_type = np.array([], dtype=int)
    bonds: Bondtype = []

    aggrecan = Aggrecan(N1, N2, n1, n2, m1, m2)
    aggrecan_btypes = np.array(aggrecan.btypes)

    k = 1
    d_x = box.x / 8
    d_y = box.y / 8
    zt = - box.z * 0.5 + length_bond * 0.5
    for i in range(8):
        xt = - box.x * 0.5 + d_x * 0.5 + d_x * i
        for j in range(8):
            yt = - box.y * 0.5 + d_y * 0.5 + d_y * j
            b_type = np.hstack([b_type, aggrecan_btypes])
            (x_t, y_t, z_t) = aggrecan.get_coords(box=box,
                                                  bond_length=length_bond,
                                                  fixed_coords={
                                                      0: (xt, yt, zt)},
                                                  periodic=False)
            x = np.hstack([x, x_t])
            y = np.hstack([y, y_t])
            z = np.hstack([z, z_t])
            b = [(bond[0] + k, bond[1] + k) for bond in aggrecan.bonds]
            bonds += b
            fixed_list.append(k)
            k += aggrecan.M
    n_solvent = int(box.volume * 3 - k)
    x = np.hstack([x, np.random.uniform(-box.x*0.5, box.x*0.5, n_solvent)])
    y = np.hstack([y, np.random.uniform(-box.y*0.5, box.y*0.5, n_solvent)])
    z = np.hstack([z, np.random.uniform(-box.z*0.5, box.z*0.5, n_solvent)])
    b_type = np.hstack([b_type, np.array([3] * n_solvent)])
    print_coord(box, x, y, z, b_type)
    print_bonds(box=box, num_atoms=len(x), bonds=bonds)
    print_fixed(fixed_list)
