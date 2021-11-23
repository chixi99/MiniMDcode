#!/usr/bin/env python3

from read_tmqm import read_tmqm_properties
import matplotlib.pyplot as plt

tmqm_properties = read_tmqm_properties('tmQM_y.csv','tmQM_X.xyz')

print("Available properties:",tmqm_properties.dtype.names)

plt.scatter(tmqm_properties["Electronic_E"], tmqm_properties["Dispersion_E"], alpha=0.5)

plt.ylabel('Dispersion_E')
plt.xlabel('Electronic_E');

plt.show()
