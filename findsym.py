# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 10:16:11 2016
modified on Mon Dec 18, 2017
@author: nenian

Simple script to convert VASP poscar inorganic 
structure files to cifs with symmetry information 
using Harold Stokes ISOTROPY for Linux 
(https://stokes.byu.edu/iso/isolinux.php)
"""

from pymatgen.io.vasp import Poscar
import os, time, string
import argparse

def findsym(args):
    POSCAR = args.parent
    poscar = Poscar.from_file(POSCAR)

    struct = poscar.structure

    atoms = dict([[str(sym),(i+1)] for i, sym in enumerate(struct.symbol_set)])
    associ = dict([[let, atom] for atom in  atoms for let in string.ascii_uppercase[(atoms[atom]-1)] ])
    lat_symprec = '1d-5'
    atom_symprec = '1d-3'
    with open('findsym_input.in', 'w') as fid:
        fid.write("{}\n".format(struct.composition.formula))
        fid.write("{}\ttolerance for lattice parameters a,b,c (angstroms)\n".format(lat_symprec))
        fid.write("{}\ttolerance for atomic positions (angstroms)\n".format(atom_symprec))
        fid.write("{}\ttolerance for magnetic moments\n".format(atom_symprec))        
        fid.write("2\tform of lattice parameters:\
        to be entered as lengths and angles\n")
        fid.write("{} {} {} {} {} {}\n".format(struct.lattice.a, struct.lattice.b, struct.lattice.c,\
        struct.lattice.alpha,struct.lattice.beta, struct.lattice.gamma))
        fid.write("2\tform of vectors defining unit cell\n")
        fid.write("P\n")
        fid.write("{}\n".format(struct.num_sites))
        for specie in struct.species:
        
            fid.write("{} ".format(atoms[str(specie.symbol)]))
        fid.write("\n")
        for vec in struct.frac_coords:
            fid.write(" ".join(str(p) for p in vec))
            fid.write("\n")

    operation = "findsym < findsym_input.in > findsym_input.out"
    os.system(operation)
    time.sleep(1) #"Failsafe!!"
    if args.output == "":
        if POSCAR.endswith(".vasp"):
            filename = POSCAR.replace(".vasp", ".cif")
            cif = open(filename, 'w')
        else:
            filename = "{}.cif".format(POSCAR)
            cif = open(filename, 'w')
    else:
        cif = open(args.output, 'w')
    with open("findsym_input.out") as sym_data:
        for line in sym_data:
            if "# CIF file" in line.strip():
                for line in sym_data:
                    cif.write(line)
                    if "_atom_site_occupancy" in line.strip():
                        for line in sym_data:
                            split_line = line.split()
                            for i in range(len(split_line)):
                                if i == 0:
                                    cif.write("{} ".format(split_line[0].replace(split_line[0][0],associ[split_line[1]] )))
                                elif i == 1:
                                
                                    cif.write("{} ".format(associ[split_line[1]]))
                                else:
                                    cif.write("{} ".format(split_line[i]))
                            cif.write("\n")
                
                
    cif.close()

    # Remove excess
    os.system("rm findsym_input.in findsym_input.out findsym.log")

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--parent", help="parent file in vasp POSCAR format",
                   required=True, metavar="FILE")
parser.add_argument("-o", "--output", help="name you want to give output cif",
                   required=False, metavar='STRING', default="", type=str)
findsym(parser.parse_args())
