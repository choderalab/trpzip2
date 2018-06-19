import time
import progressbar
from simtk import openmm, unit
from simtk.openmm import app
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('simulation_temperature', type=int)
args = parser.parse_args()
simulation_temperature = args.simulation_temperature

#
# Global simulation parameters
#

water_model = 'tip3p'
solvent_padding = 15.0 * unit.angstroms
ionic_strength = 39 * unit.millimolar

ffxml_filenames = ['amber96.xml', 'tip3p.xml']

pressure = 1.0 * unit.atmospheres
temperature = simulation_temperature * unit.kelvin
collision_rate = 1.0 / unit.picoseconds
timestep = 2.0 * unit.femtoseconds
nsteps = 500 # 1 ps
niterations = 5000 # 5 ns

solvated_pdb_filename = 'solvated.pdb'
minimized_pdb_filename = '%s/minimized.pdb' % simulation_temperature
equilibrated_pdb_filename = '%s/equilibrated.pdb' % simulation_temperature

system_xml_filename = '%s/system.xml' % simulation_temperature
integrator_xml_filename = '%s/integrator.xml' % simulation_temperature
state_xml_filename = '%s/state.xml' % simulation_temperature

# Read in the solvated model
print('Loading %s' % solvated_pdb_filename)
pdb = app.PDBFile(solvated_pdb_filename)

print("Loading forcefield: %s" % ffxml_filenames)
forcefield = app.ForceField(*ffxml_filenames)

modeller = app.Modeller(pdb.topology, pdb.positions)

# Create the system
print('Creating OpenMM System...')
system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.PME, constraints=app.HBonds, removeCMMotion=False)

# Add a barostat
print('Adding barostat...')
barostat = openmm.MonteCarloBarostat(pressure, temperature)
system.addForce(barostat)

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
initial_time = time.time()
for iteration in progressbar.progressbar(range(niterations)):
    integrator.step(nsteps)
elapsed_time = (time.time() - initial_time) * unit.seconds
simulation_time = niterations * nsteps * timestep
print('    Equilibration took %.3f s for %.3f ns (%8.3f ns/day)' % (elapsed_time / unit.seconds, simulation_time / unit.nanoseconds, simulation_time / elapsed_time * unit.day / unit.nanoseconds))
with open(equilibrated_pdb_filename, 'w') as outfile:
    app.PDBFile.writeFile(modeller.topology, context.getState(getPositions=True,enforcePeriodicBox=True).getPositions(), file=outfile, keepIds=True)
print('  final   : %8.3f kcal/mol' % (context.getState(getEnergy=True).getPotentialEnergy()/unit.kilocalories_per_mole))

# Serialize state
print('Serializing state to %s' % state_xml_filename)
state = context.getState(getPositions=True, getVelocities=True, getEnergy=True, getForces=True)
with open(state_xml_filename, 'w') as outfile:
    xml = openmm.XmlSerializer.serialize(state)
    outfile.write(xml)

# Serialize system
print('Serializing System to %s' % system_xml_filename)
system.setDefaultPeriodicBoxVectors(*state.getPeriodicBoxVectors())
with open(system_xml_filename, 'w') as outfile:
    xml = openmm.XmlSerializer.serialize(system)
    outfile.write(xml)

