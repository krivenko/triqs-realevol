from pytriqs.applications.realevol.meshes import umesh
from numpy import linspace

def print_mesh(m):
    for n,v in enumerate(m): print n,v

m1 = umesh(0,1,11)

print len(m1)
print m1[5]
print_mesh(m1)
