import numpy as np
from math_func import rotate_structure, rotate_struc_euler, vec_ang
from compare_files import compare
import math


def coord_from_system(chain_atoms, coord_places):
    chain_coords = []
    for line in chain_atoms:
        at_coord = []
        for i in coord_places:
            at_coord.append(line[i])
        chain_coords.append(at_coord)
    return np.array(chain_coords)


def copy_coord(chain_atoms, coord_places, coords):
    new_system = []
    for count, line in enumerate(chain_atoms):
        new_line = line.copy()
        for j, i in enumerate(coord_places):
            new_line[i] = coords[count][j]
        new_system.append(new_line)
    return new_system


class Polymer:
    def __init__(self):
        self.chain_coords = []
        self.chain_reses = []
        self.chain_atoms = []
        self.system = []
        self.coord_places = None
        self.second_chain_atoms = []
        self.second_chain_coords = None
        self.two_chains_atoms = []
        self.two_chains_coords = None
        self.distance_c1c2 = np.array([1.5, 8, 0])
        self.distance1 = np.array([0, 17, 0])
        self.distance2 = np.array([10, 0, 0])

    def append(self, res):
        if len(self.chain_reses) == 0:
            self.chain_reses.append(res)

        else:
            end = self.chain_reses[-1].edj_atoms_coords[1]
            begin = res.edj_atoms_coords[0]
            last_line = self.chain_reses[-1].strucs[-1]
            res.connect(point=((begin - end) + [0, 0, 2.226]), res_num=last_line[4], at_num=last_line[1])
            self.chain_reses.append(res)

    def create_one_chain(self):
        for res in self.chain_reses:
            for atom in res.strucs:
                self.chain_atoms.append(atom)
        self.coord_places = self.chain_reses[0].coord_place
        for line in self.chain_atoms:
            at_coord = []
            for i in self.coord_places:
                at_coord.append(line[i])
            self.chain_coords.append(at_coord)
        self.chain_coords = coord_from_system(coord_places=self.coord_places, chain_atoms=self.chain_atoms)
        self.chain_coords = np.array(self.chain_coords)
        # for line in self.chain_atoms:
        #     self.system.append(line.copy())

    def create_second_chain(self):
        self.second_chain_coords = rotate_struc_euler(self.chain_coords, anglez=45)
        self.second_chain_coords += self.distance_c1c2
        self.second_chain_atoms = copy_coord(chain_atoms=self.chain_atoms, coord_places=self.coord_places,
                                            coords=self.second_chain_coords)

    def move_chain(self, chain, move, level=False):
        if chain == 0 or level:
            self.chain_coords = coord_from_system(self.chain_atoms, coord_places=self.coord_places) + move
            self.chain_atoms = copy_coord(chain_atoms=self.chain_atoms, coord_places=self.coord_places, coords=self.chain_coords)
        if chain == 1 or level:
            self.second_chain_coords = coord_from_system(self.second_chain_atoms, coord_places=self.coord_places) + move
            self.second_chain_atoms = copy_coord(chain_atoms=self.second_chain_atoms, coord_places=self.coord_places,
                       coords=self.second_chain_coords)

    def multiply(self, n, m):
        for j in range(m):
            for i in range(n):
                if i % 2 == 0:
                    for line in self.chain_atoms:
                        self.system.append(line.copy())
                else:
                    for line in self.second_chain_atoms:
                        self.system.append(line.copy())
                self.move_chain(i % 2, self.distance1)
            self.move_chain(1, move=self.distance2 - self.distance1 * math.ceil(n/2), level=True)

    def write(self, file):
        with open(file, 'w') as file:
            # for res in self.system_reses:
            for count, atom in enumerate(self.system):  # res.strucs:
                if atom[2] == 'O5' and atom[3] == 'HYn' and atom[4] == 1:
                    print(atom[6])
                atom[1] = count + 1
                line = '{:4s}{:7d}  {:4s}{:8s}{:<5d}{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}\n'.format(*atom)
                file.write(line)
            file.write('TER\nEND\n')


class Residue:
    def __init__(self, input_file, coord_place=(5, 6, 7), prop_num=(1, 4, 8, 9), change_orientation=False,
                 orient_atom_names=None, edj_at=None, distance=(0, 0, 2.226)):
        self.input_file = input_file
        self.struc = []
        self.coord_place = coord_place
        self.prop_num = prop_num
        self.read()
        self.orient_atom_names = orient_atom_names
        self.begining = []
        self.end = []
        self.system_coords = []
        self.orient_atom1 = [-1, -1, -1]
        self.edj_atoms = edj_at
        self.edj_atoms_position = [0, 0]
        self.edj_atoms_coords = [0, 0]
        self.__create_orient()
        self.system_coords = np.array(self.system_coords)
        self.struc_vec1 = self.system_coords[self.orient_atom1[0]] - self.system_coords[self.orient_atom1[-1]]

        self.in_centre()

        if change_orientation:
            self.change_orientation()
        self.__refrash_coord()
        self.distance = distance

    def __create_orient(self):
        for count, line in enumerate(self.struc):
            at_coord = []
            if line[2] in self.orient_atom_names:
                if line[2] == 'C4':
                    self.orient_atom1[0] = count
                if line[2] == 'C2':
                    self.orient_atom1[1] = count
                if line[2] == 'O5':
                    self.orient_atom1[2] = count
            if line[2] == self.edj_atoms[0]:
                self.edj_atoms_position[0] = count
            if line[2] == self.edj_atoms[1]:
                self.edj_atoms_position[1] = count
            for i in self.coord_place:
                at_coord.append(line[i])
            self.system_coords.append(at_coord)

    def __refrash_coord(self):
        for i in range(len(self.system_coords)):
            for count, coord in enumerate(self.coord_place):
                self.struc[i][coord] = self.system_coords[i][count]
        self.edj_atoms_coords[0] = self.system_coords[self.edj_atoms_position[0]]
        self.edj_atoms_coords[1] = self.system_coords[self.edj_atoms_position[1]]

    def read(self):
        with open(self.input_file) as file:
            for line in file:
                atom = line.split()
                for i in self.coord_place:
                    atom[i] = float(atom[i])
                for i in self.prop_num:
                    atom[i] = int(round(float(atom[i])))
                self.struc.append(atom)

    def write(self, file):
        with open(file, 'w') as file:
            for atom in self.strucs:
                line = '{:4s}{:7d}  {:4s}{:8s}{:<5d}{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}\n'.format(*atom)
                file.write(line)
            file.write('TER\nEND\n')

    def in_centre(self):
        pass
        self.system_coords -= self.system_coords[self.orient_atom1[-1]]
        self.system_coords = rotate_structure(structure=self.system_coords, str_vector=self.struc_vec1,
                                              oriented_vector=(0, 0, 1))
        self.system_coords[:, -1] -= self.system_coords[self.orient_atom1[1]][-1]
        self.system_coords = rotate_struc_euler(structure=self.system_coords,
                                                anglez=vec_ang(self.system_coords[self.orient_atom1[1]], (0, 1, 0),
                                                               degree=True))

    def change_orientation(self):
        self.system_coords = rotate_struc_euler(structure=self.system_coords, anglez=180)

    def connect(self, point, res_num, at_num):
        self.system_coords -= point
        for line in self.struc:
            line[1] += at_num
            line[4] += res_num
        self.__refrash_coord()

    def replace(self):
        pass

    @property
    def dist(self):
        return self.distance

    @property
    def edje_atoms(self):
        return self.edj_atoms_coords

    @property
    def strucs(self):
        return self.struc.copy()

    @property
    def coords(self):
        return self.system_coords


class CHT0(Residue):
    def __init__(self, change_orientation=0, orient_atom_names=None, edj_atoms=None):
        super().__init__(input_file='files\\HYn.pdb', change_orientation=change_orientation,
                         orient_atom_names=orient_atom_names, edj_at=edj_atoms)


class CHT(Residue):
    def __init__(self, change_orientation=0, orient_atom_names=None, edj_atoms=None):
        super().__init__(input_file='files\\BYn.pdb', change_orientation=change_orientation,
                         orient_atom_names=orient_atom_names, edj_at=edj_atoms)


class CHTN(Residue):
    def __init__(self, change_orientation=0, orient_atom_names=None, edj_atoms=None):
        super().__init__(input_file='files\\TYn.pdb', change_orientation=change_orientation,
                         orient_atom_names=orient_atom_names, edj_at=edj_atoms)


or_at = {'C2', 'O5', 'C4'}
edj_at = ('C4', 'C1')
polymer = Polymer()
polymer.append(CHT0(change_orientation=False, orient_atom_names=or_at, edj_atoms=edj_at))
for i in range(18):
    polymer.append(CHT(change_orientation=(i % 2 == 0), orient_atom_names=or_at, edj_atoms=edj_at))
polymer.append(CHTN(change_orientation=(8 % 2 == 0), orient_atom_names=or_at, edj_atoms=edj_at))
polymer.create_one_chain()
polymer.create_second_chain()
polymer.multiply(4, 4)
polymer.write('cht_20mon_16chains.pdb')
# compare('test2.pdb', 'test1.pdb')
