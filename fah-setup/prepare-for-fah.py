#!/bin/env python

import progressbar
from simtk import openmm, unit
from simtk.openmm import app

#
# Global simulation parameters
#

water_model = 'tip3p'
solvent_padding = 15.0 * unit.angstroms
ionic_strength = 200 * unit.millimolar
hydrogen_mass = 2.0 * unit.amu

ffxml_filenames = ['amber14/protein.ff14SB.xml', 'amber14/tip3p.xml']

pressure = 1.0 * unit.atmospheres
temperature = 425 * unit.kelvin
collision_rate = 1.0 / unit.picoseconds
timestep = 2.0 * unit.femtoseconds
nsteps = 500 # 1 ps
niterations = 1000 # 1 ns

solvated_pdb_filename = 'solvated.pdb'
minimized_pdb_filename = 'minimized.pdb'
equilibrated_pdb_filename = 'equilibrated.pdb'

system_xml_filename = 'system.xml'
integrator_xml_filename = 'integrator.xml'
state_xml_filename = 'state.xml'

# Read in the NMR model
pdb_filename = '1le1.pdb'
print('Loading %s' % pdb_filename)
pdb = app.PDBFile(pdb_filename)

# Use Amber 14SB (which has parameters for C-terminal -NH2)
print("Loading forcefield: %s" % ffxml_filenames)
forcefield = app.ForceField(*ffxml_filenames)

# Solvate
print('Adding solvent...')
modeller = app.Modeller(pdb.topology, pdb.positions)
modeller.addSolvent(forcefield, model=water_model, padding=solvent_padding, ionicStrength=ionic_strength)
print('System has %d atoms' % modeller.topology.getNumAtoms())

# Write initial model
print('Writing initial solvated system to %s' % solvated_pdb_filename)
with open(solvated_pdb_filename, 'w') as outfile:
    app.PDBFile.writeFile(modeller.topology, modeller.positions, file=outfile, keepIds=True)

# Create the system
print('Creating OpenMM System...')
system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.PME, constraints=app.HBonds, removeCMMotion=False, hydrogenMass=hydrogen_mass)

# Add a barostat
print('Adding barostat...')
barostat = openmm.MonteCarloBarostat(pressure, temperature)
system.addForce(barostat)

# Serialize system
print('Serializing System to %s' % system_xml_filename)
with open(system_xml_filename, 'w') as outfile:
    xml = openmm.XmlSerializer.serialize(system)
    outfile.write(xml)

# Serialize integrator
print('Serializing integrator to %s' % integrator_xml_filename)
integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)
with open(integrator_xml_filename, 'w') as outfile:
    xml = openmm.XmlSerializer.serialize(integrator)
    outfile.write(xml)

# Minimize
print('Minimizing energy...')
context = openmm.Context(system, integrator)
context.setPositions(modeller.positions)
print('  initial : %8.3f kcal/mol' % (context.getState(getEnergy=True).getPotentialEnergy()/unit.kilocalories_per_mole))
openmm.LocalEnergyMinimizer.minimize(context)
print('  final   : %8.3f kcal/mol' % (context.getState(getEnergy=True).getPotentialEnergy()/unit.kilocalories_per_mole))
with open(minimized_pdb_filename, 'w') as outfile:
    app.PDBFile.writeFile(modeller.topology, context.getState(getPositions=True,enforcePeriodicBox=True).getPositions(), file=outfile, keepIds=True)

# Equilibrate
print('Equilibrating...')
for iteration in progressbar.progressbar(range(niterations)):
    integrator.step(nsteps)
with open(equilibrated_pdb_filename, 'w') as outfile:
    app.PDBFile.writeFile(modeller.topology, context.getState(getPositions=True,enforcePeriodicBox=True).getPositions(), file=outfile, keepIds=True)
print('  final   : %8.3f kcal/mol' % (context.getState(getEnergy=True).getPotentialEnergy()/unit.kilocalories_per_mole))

# Serialize state
print('Serializing state to %s' % state_xml_filename)
state = context.getState(getPositions=True, getVelocities=True, getEnergy=True, getForces=True)
with open(state_xml_filename, 'w') as outfile:
    xml = openmm.XmlSerializer.serialize(state)
    outfile.write(xml)
