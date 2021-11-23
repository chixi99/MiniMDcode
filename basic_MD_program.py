"""Demonstrates molecular dynamics with constant energy."""

from asap3 import Trajectory
from ase.lattice.cubic import FaceCenteredCubic
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase import units

def init():
    """Initialization of MD program"""
    sumv=0
    sumv2=0
    for i in range(npart):
        x[i] = lattice_pos(i)    # Place the particles on a lattice
        v[i] = (random() - 0.5)  # Give random velocities
        sumv += m[i]*v[i]        # Velocity center of mass
        sumv2 += m[i]*v[i]**2    # Kinetic energy
        summ += m[i]             # Total mass
    sumv = sumv/(npart*summ)     # Velocity center of mass
    sumv2 = sumv2/npart          # Mean squared velocity
    fs = sqrt(3*kB*temp/sumv2)   # Scale factor of the velocities

    for i in range(npart):       # set desired kinetic energy and set
        v[i] = (v[i] - sumv)*fs  # velocity center of mass to zero
        xm[i] = x[i] - v[i]*dt   # position previous time step

def lattice_pos(i):
    """Return the coordinates of ideal lattice position i"""


def force(f, en):
    """Determine the force and energy """
    en=0
    f[:] = 0

    for i in range(npart-1):
      for j in range(i+1, npart):
        xr = x[i] - x[j]
        xr -= box*round(xr/box)
        r2=xr**2
        if r2 < rc2:
          r2i = sigma[i][j]/r2
          r6i = r2i**3
          ff = 48*r2i*r6i*(r6i-0.5)
          f[i] += ff*xr
          f[j] -= ff*xr
          ecut = 4*(1/rc2**6 - 1/rc2**3) # Potential at cutoff
          en += 4*r6i*(r6i-1)-ecut # Pot energy for all pairs
    en = en*epsilon
    return f, en

def integrate(f, en):
"""Integrate equations of motion """
    sumv = 0
    sumv2 = 0
    for i in range(npart):                # Verlet algorithm
        xx = 2*x[i] - xm[i] + dt**2*f[i]  # New position
        vi = (xx-xm[i])/(2*dt)            # New velocity
        sumv += vi                        # Velocity center of mass
        sumv2 += vi**2                    # Total kinetic energy
        xm[i]=x[i]                        # Update positions prev time step
        x[i]=xx                           # Update positions current time step
    temp = sumv2/(3*npart)                # Instantaneous temperature
    etot = (en+0.5*sumv2)/npart           # Add kinetic energy to total energy

def sample():
"""Calculate quantities"""
    # Instantaneous potential energy calculated each step in force()
    epot = en/npart
    # Instantaneous total kinetic energy
    ekin = 0
    for i in range(npart):
        ekin += m[i]*v[i]**2 ekin = sumv2/npart
    print('Energy per atom: Epot = %.3feV Ekin = %.3feV (T=%3.0fK) '
          'Etot = %.3feV' % (epot, ekin, ekin / (1.5 * units.kB), epot + ekin))


"""
Molecular dynamics program using Verlet integrator
"""
init()
t=0
while t < tmax:
# Initalization
# MD loop
# Determine forces on all atoms
# Integrate equations of motion
# Calculate averages, e.g., pressure, temp.
f, en = force()
integrate(f, en)
sample()
t = t + dt # <=====================================================
