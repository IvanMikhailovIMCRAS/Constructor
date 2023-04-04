from bondset import Bondtype

from constructor import MolGraph


class Chain(MolGraph):
    """
    Chain molecular structure inherited from MolGraph

    Args:
        MolGraph (_type_): see constructor.MolGraph
    Initial attribute:
        self.n (int): chain length
    """

    def __init__(self, n: int) -> None:
        if n < 1:
            raise ValueError(f'(n = {n}) set is invalid ')
        self.n = n
        bonds: Bondtype = []
        for i in range(self.n-1):
            bonds.append((i, i+1))
        super().__init__(bonds)


class Aggrecan(MolGraph):
    """
    Aggrecan molecular structure inherited from MolGraph

    Args:
        MolGraph (_type_): see constructor.MolGraph
    Initial attribute:
        self.N1 (int): backbone lenth (root block)
        self.N2 (int): backbone lenth (head block)
        self.n1 (int): side chain length (root block)
        self.n2 (int): side chain length (head block)
        self.m1 (int): spacer length (root block)
        self.m2 (int): spacer length (head block)
    """

    def __init__(self, N1: int, N2: int, n1: int, n2: int, m1: int, m2: int) -> None:
        if m1 < 1 or m2 < 1:
            raise ValueError(f'(m1 = {m1}) and (m2 = {m2}) sets are invalid ')
        if N1 % m1 != 0 or N2 % m2 != 0:
            raise ValueError(
                'N1 and / or N2 must be evenly divisible by m1 and / or m2')
        if N1 <= 1 or N2 <= 0:
            raise ValueError('values must be N1 >= 1 or N2 >= 0')
        if n1 <= 0 or n2 <= 0:
            raise ValueError('values must be n1 >= 0 or n2 >= 0 ')

        self.N1 = N1
        self.N2 = N2
        self.n1 = n1
        self.n2 = n2
        self.m1 = m1
        self.m2 = m2
        self.M = N1 + N2 + N1 * n1 // m1 + N2 * n2 // m2  # Polymerization degree

        bonds: Bondtype = []
        self.btypes = []

        self.btypes.append(1)
        for i in range(self.N1-1):
            bonds.append((i, i+1))
            self.btypes.append(1)
        current_N = self.N1-1

        for i in range(self.N2):
            bonds.append((i+current_N, i+current_N+1))
            self.btypes.append(2)
        current_N += self.N2

        for i in range(N1 // m1):
            current_N += 1
            self.btypes.append(4)
            nc = m1 // 2 + i * m1
            bonds.append((nc, current_N))
            for j in range(1, n1):
                current_N += 1
                self.btypes.append(4)
                bonds.append((current_N-1, current_N))

        for i in range(N2 // m2):
            current_N += 1
            self.btypes.append(5)
            nc = N1 + m2 // 2 + i * m2
            bonds.append((nc, current_N))
            for j in range(1, n2):
                current_N += 1
                self.btypes.append(5)
                bonds.append((current_N-1, current_N))

        super().__init__(bonds, sort=False)


class Dendron(MolGraph):
    """
    Dendrit molecular structure inherited from MolGraph

    Args:
        MolGraph (_type_): see constructor.MolGraph
    Initial attributes:
        self.n (int): spacer length
        self.g (int): number of generations in dendron
    """

    def __init__(self, n: int, g: int = 0, q: int = 2) -> None:
        if q <= 1 or g < 0 or n < 1:
            raise ValueError(f'(n = {n}, g = {g}, q = {q}) set is invalid ')
        self.n = n
        self.g = g
        self.q = q
        bonds: Bondtype = []
        self.num_spacers: int = (self.q**(self.g) - 1) // (self.q - 1)
        for i in range(self.n):
            bonds.append((i, i+1))
        id0 = 0
        current_id = self.n
        for s in range(self.num_spacers):
            id0 += self.n
            for i in range(self.q):
                for j in range(self.n):
                    current_id += 1
                    if j == 0:
                        bonds.append((id0, current_id))
                    else:
                        bonds.append((current_id-1, current_id))
        super().__init__(bonds, sort=False)
        # super().__init__(bonds)


class Brush(MolGraph):
    """Molecular brush structure

    Args:
        MolGraph (_type_): see constructor.MolGraph
    Initial attributes:
        self.pd (int): polymerization degree
        self.m (int): spacer length
        self.n (int): side chain length
        self.q (int): number of side chains grafting into one point
        self.n_end_ch (int): number of ending chains
        self.l_end_ch (int): length of ending chains
    """

    def __init__(self, pd: int, m: int, n: int, l_end_ch: int, q: int = 1, n_end_ch: int = 0) -> None:
        if pd < 1 or m < 1 or l_end_ch < 0:
            raise ValueError(
                f'(pd = {pd}, m = {m}, l_end_ch = {l_end_ch}) set is invalid - the main core is empty')
        if n < 0 or q < 0:
            raise ValueError(
                f'(n = {n}, q = {q}) set is invalid - length and number of side chains must be positive')
        if n_end_ch > 2 or n_end_ch < 0:
            raise ValueError(
                f'(n_end_ch ={n_end_ch}) set is invalid - n_end_ch should be 0 or 1 or 2')
        self.pd = pd
        self.m = m
        self.n = n
        self.q = q
        self.n_end_ch = n_end_ch
        self.l_end_ch = l_end_ch
        self.types = list()
        bonds: Bondtype = []
        self.types = [1] * self.l_end_ch + [2] * \
            self.pd * self.m + [1] * self.l_end_ch

        for i in range(self.n_end_ch * self.l_end_ch + self.pd * self.m - 1):
            bonds.append((i, i+1))
            if (i < self.l_end_ch) or (i >= self.l_end_ch + self.pd * self.m):
                self.types[i] = 1
                if i == self.n_end_ch * self.l_end_ch + self.pd * self.m - 2:
                    self.types[i + 1] = 1
            else:
                self.types[i] = 2

        n_curent = self.pd * self.m - 1 + self.n_end_ch * self.l_end_ch
        for i in range(1, self.pd + 1):
            if self.n_end_ch == 0:
                num = (self.m * (i-1)) + (self.m // 2 + 1) - 1
            else:
                num = (self.m * (i-1)) + (self.m // 2 + 1) - 1 + self.l_end_ch
            for t in range(self.q):
                for j in range(self.n):
                    n_curent += 1
                    if j == 0:
                        bonds.append((num, n_curent))
                    else:
                        bonds.append((n_curent-1, n_curent))

        super().__init__(bonds, sort=False)
        # for element in bonds:
        #     for bond in element:
        #         if bond not in self.types:
        #             self.types[bond] = 2
        self.types = self.types + \
            [3] * (self.num_beads - self.n_end_ch *
                   self.l_end_ch - self.pd*self.m)


if __name__ == '__main__':
    x = Aggrecan(N1=10, N2=10, n1=2, n2=3, m1=5, m2=1)
    print(x.bonds)
    print(x.M)
    print(x.btypes)
    print(len(x.btypes))
