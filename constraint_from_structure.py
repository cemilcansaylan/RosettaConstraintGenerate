#!/usr/bin/env python
# coding: utf-8

"""
Created on Wed Sep 15 22:27:56 2021

@author: cemil can saylan
"""

from Bio.PDB import *
from Bio import SeqIO
import numpy

from absl import app
from absl import flags

flags.DEFINE_string('structure', None, 'PDB Structure')
flags.DEFINE_string('output', 'constraint.cst', 'constraint file name')
flags.DEFINE_boolean('printfasta', False, 'Print sequence of structure')

FLAGS = flags.FLAGS
flags.mark_flag_as_required("structure")

RADIUS=20

def PrintFasta(structure):
    f = open("input.fasta", "w")
    for record in SeqIO.parse(structure, "pdb-atom"):
        f.write(str(record.seq))
    f.close()

def AtomPairNearest(structure):
    """ 
    Only for CA 
    """
    
    res_list = Selection.unfold_entities(structure, "R")
    atoms  = Selection.unfold_entities(structure, 'A')
    
    ns = NeighborSearch(atoms)
    Nearest_Atom = []
    
    for res in res_list:
        near = ns.search(res["CA"].coord, RADIUS)
        for a in near:
            if a.id == "CA":
                diff = res["CA"].coord - a.coord
                dist = numpy.sqrt(numpy.sum(diff * diff))
                if dist == 0.0:
                    continue
                Nearest_Atom.append([res["CA"].id, 
                                     res["CA"].get_parent().id[1],
                                     a.id, 
                                     a.get_parent().id[1], 
                                     dist])
    return Nearest_Atom

def Dih_Angle(structure):
    """
    Get the phi, psi information from the structure
    """
    dihedrals = []
    angles = []
    structure.atom_to_internal_coordinates()
    for r in structure.get_residues():
        if r.internal_coord:
            angles.append(r.internal_coord.get_angle("N:CA:C"))
            dihedrals.append([r.internal_coord.get_angle("phi"),
                             r.internal_coord.get_angle("psi")])
    return dihedrals, angles

def PrintConstraint(pdb, filename):
    """     
    psi N-Cα-C-N
    phi C-N-Cα-C
    """
    res_list = Selection.unfold_entities(pdb, "R")
    
    Nearest_Atom = AtomPairNearest(pdb)
    Dihedrals, Angles = Dih_Angle(pdb)
    
    f = open(filename, "w")
    for res_num in range(1,len(res_list)+1):
        if res_num == 1:
            f.write(f"Dihedral N {res_num} CA {res_num} C {res_num} N {res_num+1} LINEAR_PENALTY {Dihedrals[res_num-1][1]*0.01745329252} 0 0 10.0\n")
            f.write(f"Angle N {res_num} CA {res_num} C {res_num} LINEAR_PENALTY {Angles[res_num-1]} 0 0 10.0\n")
        elif res_num == len(res_list):
            f.write(f"Dihedral C {res_num-1} N {res_num} CA {res_num} C {res_num} LINEAR_PENALTY {Dihedrals[res_num-1][0]*0.01745329252} 0 0 10.0\n")
            f.write(f"Angle N {res_num} CA {res_num} C {res_num} LINEAR_PENALTY {Angles[res_num-1]} 0 0 10.0\n")
        else:
            f.write(f"Dihedral N {res_num} CA {res_num} C {res_num} N {res_num+1} LINEAR_PENALTY {Dihedrals[res_num-1][1]*0.01745329252} 0 0 10.0\n")
            f.write(f"Dihedral C {res_num-1} N {res_num} CA {res_num} C {res_num} LINEAR_PENALTY {Dihedrals[res_num-1][0]*0.01745329252} 0 0 10.0\n")
            f.write(f"Angle N {res_num} CA {res_num} C {res_num} LINEAR_PENALTY {Angles[res_num-1]} 0 0 10.0\n")

            
    for atompair in Nearest_Atom:
        f.write(f"AtomPair CA {atompair[1]} CA {atompair[3]} LINEAR_PENALTY {atompair[4]} 0 0 10.0\n")
    
    f.close()
  
def main(argv):
    parser = PDBParser()
    pdb_file = parser.get_structure("demo", FLAGS.structure)
    
    PrintConstraint(pdb_file, FLAGS.output)
    
    if FLAGS.printfasta:
        PrintFasta(FLAGS.structure)
 
if __name__ == '__main__':
    app.run(main)

