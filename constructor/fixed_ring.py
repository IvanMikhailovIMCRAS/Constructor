import numpy as np
from bondset import Bondtype
from dpd_files import print_bonds, print_coord, print_fixed
from mol import Dendron, Chain
from periodic_box import Box

class Ring:
  
    def __init__(self, num_beads, length_bond):
        self.num_beads = num_beads
        self.lenth_bond = length_bond
        self.radius = length_bond / (2*np.sin(np.pi/num_beads))

    def get_bonds(self):
        bonds = []
        for i in range(1, self.num_beads):
            bonds.append([i, i+1])
        bonds.append([1, self.num_beads])
        return bonds

    def get_coords(self):
        x = []
        y = []
        z = []
        for i in range(self.num_beads):
            alfa = 2*i*np.pi/self.num_beads
            x.append(self.radius*np.cos(alfa))
            y.append(self.radius*np.sin(alfa))
            z.append(0.0)
        return x, y, z
    
    def get_angles(self):
        angles = []
        for i in range(1, self.num_beads-1):
            angles.append([i, i+1, i+2])
        angles.append([self.num_beads-1, self.num_beads, 1])
        angles.append([self.num_beads, 1, 2])
        return angles 

if __name__ == '__main__':
    main_chain = 100
    bufer = 0
    n = 3
    g = 3
    m = 1
    length_bond = (1/3)**(1/3)
    h_xyz = main_chain * length_bond / np.pi + 11.0
    box = Box(x=h_xyz, y=h_xyz, z=h_xyz)
    x = np.array([])
    y = np.array([])
    z = np.array([])
    b_type = np.array([], dtype=int)
    bonds: Bondtype = []
    fixed_list = []
    ring = Ring(num_beads=main_chain, length_bond=length_bond)
    dendron = Dendron(n=n, g=g)
    (x_0, y_0, z_0) = ring.get_coords()
    k = 1
    for i in range(main_chain):
        fixed_list.append(k)
        if i%m == 0:
            b_type = np.hstack(
                [b_type, np.array([1] + [2] * (dendron.num_beads - 1))])
            d = {0: (x_0[i], y_0[i], z_0[i])}
            (x_t, y_t, z_t) = dendron.get_coords(box=box, fixed_coords=d)
            x = np.hstack([x, x_t])
            y = np.hstack([y, y_t])
            z = np.hstack([z, z_t])
            b = [(bond[0] + k, bond[1] + k) for bond in dendron.bonds]
            bonds += b
            k += dendron.num_beads
            bonds += [(k-dendron.num_beads, k)]
    bonds.pop()
    bonds += [(1, k-dendron.num_beads)]
    n_solvent = int(box.volume * 3 - k)
    x = np.hstack([x, np.random.uniform(-box.x*0.5, box.x*0.5, n_solvent)])
    y = np.hstack([y, np.random.uniform(-box.y*0.5, box.y*0.5, n_solvent)])
    z = np.hstack([z, np.random.uniform(-box.z*0.5, box.z*0.5, n_solvent)])
    b_type = np.hstack([b_type, np.array([3] * n_solvent)])
    print_coord(box, x, y, z, b_type)
    print_bonds(box=box, num_atoms=len(x), bonds=bonds)
    print_fixed(fixed_list)