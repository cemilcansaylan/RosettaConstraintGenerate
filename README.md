# RosettaConstraintGenerate

Generate Rosetta constraints (Angle (N,Ca,C), Dihedral (phi, psi), AtomPair (between Ca) from given structure

#### Required packages:
```
Biopython
numpy
absl
```

#### Install Required Packages:
```bash
pip install biopython
pip install numpy
pip install absl-py
```
#### Usage:
```
python3 constraint_from_structure.py --structure demo.pdb --distance 20
```
