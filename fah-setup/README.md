# trpzip2

trpzip2 peptide in explicit solvent with counterions

## Prerequisites

```bash
conda install --yes progressbar2 openmm==6.3.1
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
