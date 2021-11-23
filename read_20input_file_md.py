"""Demonstrates molecular dynamics with constant energy."""
# Import Module
import os
from asap3 import Trajectory
from ase.lattice.cubic import FaceCenteredCubic
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase import units
from ase.io import read
import matplotlib.pyplot as plt
import numpy as np

# Use Asap for a huge performance increase if it is installed
use_asap = True

if use_asap:
    from asap3 import EMT
    size = 10
else:
    from ase.calculators.emt import EMT
    size = 3

def calcenergy(atoms):
    """Calculate the potential, kinetic and total energy per atom."""
    epot = atoms.get_potential_energy() / len(atoms)
    ekin = atoms.get_kinetic_energy() / len(atoms)
    Tinst = ekin / (1.5 * units.kB)
    etot = epot + ekin
    return (epot, ekin, Tinst, etot)

def readinputfile(file_name):
    # read input cif file
    atoms = read(file_name)
    return atoms

def run_md(atoms):
    # Set up a crystal
    #atoms = FaceCenteredCubic(directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    #                          symbol="Cu",
    #                          size=(size, size, size),
    #                          pbc=True)

    # Describe the interatomic interactions with the Effective Medium Theory
    atoms.calc = EMT()

    # Set the momenta corresponding to T=300K
    MaxwellBoltzmannDistribution(atoms, temperature_K=300)

    # We want to run MD with constant energy using the VelocityVerlet algorithm.
    dyn = VelocityVerlet(atoms, 5 * units.fs)  # 5 fs time step.

    #traj = Trajectory("CuNi_mp-1225687.traj", "w", atoms)
    #dyn.attach(traj.write, interval=10)

    def printenergy(a=atoms):  # store a reference to atoms in the definition.
        """Function to print the potential, kinetic and total energy."""
        res = calcenergy(a)
        print('Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  '
              'Etot = %.3feV' % res)

    # Now run the dynamics
    dyn.attach(printenergy, interval=100)
    printenergy()
    res = calcenergy(atoms)
    dyn.run(200)
    return res

if __name__ == "__main__":
    # Folder Path
    path = "./test_samples"

    # Change the directory
    os.chdir(path)

    # define the array to store the data
    res_20 = []

    # iterate through all file
    for file in os.listdir():
        # Check whether file is in cif format or not
        if file.endswith(".cif"):
            atoms = readinputfile(file)
            res_20.append(run_md(atoms))

    res_array = np.array(res_20)
    #print(res_array[:, 1])
    #scatter section
    plt.scatter(res_array[:, 0], res_array[:, 1], alpha=0.5)

    plt.ylabel('Kinetic_E')
    plt.xlabel('Potential_E')

    plt.show()