import time
import progressbar
from simtk import openmm, unit
from simtk.openmm import app

#
# Global simulation parameters
#

water_model = 'tip3p'
solvent_padding = 15.0 * unit.angstroms
ionic_strength = 39 * unit.millimolar

ffxml_filenames = ['amber96.xml', 'tip3p.xml']

pressure = 1.0 * unit.atmospheres
#temperature = simulation_temperature * unit.kelvin
collision_rate = 1.0 / unit.picoseconds
timestep = 2.0 * unit.femtoseconds
nsteps = 500 # 1 ps
niterations = 5000 # 5 ns

solvated_pdb_filename = 'solvated.pdb'

# Read in the NMR model
pdb_filename = '1le1.pdb'
print('Loading %s' % pdb_filename)
pdb = app.PDBFile(pdb_filename)

print("Loading forcefield: %s" % ffxml_filenames)
forcefield = app.ForceField(*ffxml_filenames)

# Solvate
print('Adding solvent...')
modeller = app.Modeller(pdb.topology, pdb.positions)
modeller.addSolvent(forcefield, model=water_model, padding=solvent_padding, ionicStrength=ionic_strength)
#print('System has %d atoms' % modeller.topology.getNumAtoms())

# Write initial model
print('Writing initial solvated system to %s' % solvated_pdb_filename)
with open(solvated_pdb_filename, 'w') as outfile:
    app.PDBFile.writeFile(modeller.topology, modeller.positions, file=outfile, keepIds=True)


