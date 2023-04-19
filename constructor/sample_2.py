import numpy as np
from bondset import Bondtype
from dpd_files import print_bonds, print_coord, print_ent_file
from mol import Brush
from periodic_box import Box
import shutil
import os
def main(m, pd, n, q, n_end_ch, l_end_ch):
    num_ch = 1000 #number of chains
    

    # # m = 3 # spacer length
    # pd = 75/m #polymerization degree
    # # n = 4 # side chain length
    # q = 1 #number of side chains grafting into one point
    # n_end_ch = 2 #number of ending chains
    # l_end_ch = 4 # length of ending chains 
    box_size = (num_ch*(pd * m + pd * n * q + n_end_ch * l_end_ch) / 3) ** (1 / 3)
    
    box = Box(x=box_size, y=box_size, z=box_size)
    x = np.array([])
    y = np.array([])
    z = np.array([])
    b_type = np.array([], dtype=int)
    bonds: Bondtype = []
    fixed_list = []
    b_type = []
    brush = Brush(pd=pd, m=m, n=n, q=q, n_end_ch=n_end_ch, l_end_ch=l_end_ch)
    cur = 1
    for i in range(num_ch):
        # cr = ( 0.5 *  box_size * np.random.rand(), 0.5 *  box_size * np.random.rand(), 0.5 *  box_size * np.random.rand())
        (x_t, y_t, z_t) = brush.get_coords(box=box)
        for tp in brush.types:
            b_type.append(tp)
        b = [(bond[0] + cur, bond[1] + cur) for bond in brush.bonds]
        bonds += b
        cur += brush.num_beads 
        x = np.hstack([x, x_t])
        y = np.hstack([y, y_t])
        z = np.hstack([z, z_t])
       
    print_coord(box, x, y, z, b_type)
    print_bonds(box=box, num_atoms=len(x), bonds=bonds)
    print_ent_file(x=x, y=y, z=z, b_type=b_type)

if __name__ == '__main__':
    total_length = 50
    if not os.path.isdir('OUTPUT'):
        os.mkdir('OUTPUT')
    os.chdir('OUTPUT')
    for l_end_ch in range(1,18):
        if not os.path.isdir(f"lenf{l_end_ch}_pd{total_length-l_end_ch}"):
            os.mkdir(f"lenf{l_end_ch}_pd{total_length-l_end_ch}")
        os.chdir(f"lenf{l_end_ch}_pd{total_length-l_end_ch}")
        if not os.path.exists("dpd.exe"):
            shutil.copy('../../dpd.exe','dpd.exe')
        if not os.path.exists("FIELD"):
            shutil.copy('../../FIELD','FIELD')
        if not os.path.exists("CONTR"):
            shutil.copy('../../CONTR','CONTR')
        if not os.path.exists("COORD"):    
            main(m=1, pd=total_length-l_end_ch, n=0, q=1, n_end_ch=1, l_end_ch=l_end_ch) 
        os.chdir('../')   # for n in [4, 8, 12, 16, 20]:
    #     if not os.path.isdir(f'n={n}'):
    #         os.mkdir(f'n={n}')
    #     os.chdir(f'n={n}')
           
    #     for m in [1, 3, 5, 15]:
    #         if not os.path.isdir(f'n={n}'):
    #             os.mkdir(f'n={n}')
    #         os.chdir(f'n={n}')
    #         main(m, n)
    #         os.chdir('../')
    #     os.chdir('../')     
            
