# trpzip2

trpzip2 peptide in explicit solvent with counterions

## Prerequisites

```bash
conda create -n openmm631 python==3.5 openmm==6.3.1 progressbar2
source activate openmm631
```

## Manifest

* `1le1.pdb` - NMR structure
* `prepare-for-fah.py` - script to solvate, minimize, parameterize, and equilibrate with OpenMM
* `initial.pdb` - initial solvated structure
* `minimized.pdb` - minimized structure
* `equilibrated.pdb` - PDB file of equilibrated structure
* `system.xml` - serialized System
* `integrator.xml` - serialized integrator
* `state.xml` - serialized equilibrated state

## Sample output

```
[chodera@lilac:fah-setup]$ python prepare-for-fah.py 
Loading 1le1.pdb
Loading forcefield: ['amber14/protein.ff14SB.xml', 'amber14/tip3p.xml']
Adding solvent...
System has 12552 atoms
Writing initial solvated system to solvated.pdb
Creating OpenMM System...
Adding barostat...
Serializing System to system.xml
Serializing integrator to integrator.xml
Minimizing energy...
  initial : -17323.150 kcal/mol
  final   : -52095.049 kcal/mol
Equilibrating...
100% (5000 of 5000) |#################################################################################################################################################################################################################| Elapsed Time: 0:25:10 Time:  0:25:10    Equilibration took 1512.053 s for 5.000 ns ( 285.704 ns/day)
  final   : -35865.146 kcal/mol
Serializing state to state.xml
```
