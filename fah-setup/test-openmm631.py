#!/bin/env python

from simtk import openmm, unit

system_xml_filename = 'system.xml'
integrator_xml_filename = 'integrator.xml'
state_xml_filename = 'state.xml'

print('Deserializing...')
with open(system_xml_filename, 'r') as infile:
    system = openmm.XmlSerializer.deserialize(infile.read())

with open(integrator_xml_filename, 'r') as infile:
    integrator = openmm.XmlSerializer.deserialize(infile.read())

with open(state_xml_filename, 'r') as infile:
    state = openmm.XmlSerializer.deserialize(infile.read())

print('Creating Context...')
context = openmm.Context(system, integrator)
context.setPeriodicBoxVectors(*state.getPeriodicBoxVectors())
context.setPositions(state.getPositions())
print('Stepping...')
integrator.step(500)
print('Done.')
