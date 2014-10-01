from pytriqs.applications.realevol.realevol import rmesh
from numpy import linspace

def print_mesh(m):
    for n,v in enumerate(m): print n,v

m1 = rmesh(0,1,11)

print len(m1)
print m1[5]
print_mesh(m1)
