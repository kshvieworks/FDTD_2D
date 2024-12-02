import numpy as np
import matplotlib.pyplot as plt

import meep as mp
from meep.materials import SiO2
from IPython.display import Video

resolution = 100  # pixels/um

dpml = 1.0
pml_layers = [mp.PML(thickness=dpml)]

r = 1.0     # radius of cylinder
dair = 2.0  # air padding thickness

s = 2*(dpml+dair+r)
cell_size = mp.Vector3(s,s)

wvl = 1.0
fcen = 1/wvl

# is_integrated=True necessary for any planewave source extending into PML
sources = [mp.Source(mp.ContinuousSource(wavelength=wvl, width=20),
                     center=mp.Vector3(-0.5*s+dpml),
                     size=mp.Vector3(0,s),
                     component=mp.Ez)]

symmetries = [mp.Mirror(mp.Y)]

geometry = [mp.Cylinder(material=SiO2,
                        center=mp.Vector3(),
                        radius=r,
                        height=mp.inf)]

sim = mp.Simulation(resolution=resolution,
                    cell_size=cell_size,
                    boundary_layers=pml_layers,
                    sources=sources,
                    k_point=mp.Vector3(),
                    symmetries=symmetries,
                    geometry=geometry)

f = plt.figure(dpi=150)
sim.plot2D(ax = f.gca())
plt.show()

f = plt.figure(dpi=150)
Animate = mp.Animate2D(sim, fields=mp.Ez, f=f, realtime=False, normalize=True)
sim.run(mp.at_every(0.5,Animate),until=100)
plt.close()

filename = "media/nanosphere.mp4"
fps = 10
Animate.to_mp4(fps,filename)
Video(filename)