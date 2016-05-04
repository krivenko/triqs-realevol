from pytriqs.gf.local import *
from realevol.gf_retime import *
import numpy as np

m1 = MeshReTime(-1,1,50)
m2 = MeshReTime(-2,2,100)

m = MeshReTime2(m1,m2)
assert m.size_of_components() == [len(m1), len(m2)]
assert len(m) == len(m1)*len(m2)
assert m.components() == (m1,m2)

g1 = GfReTime2(window = (-3,3), n_points = 60, indices = ["a","b"], name = "g1")
assert g1.name == "g1"
assert g1.mesh == MeshReTime2(MeshReTime(-3,3,60), MeshReTime(-3,3,60))
assert g1.target_shape == [2,2]
assert (g1.data == np.zeros((60,60,2,2))).all()

g2 = GfReTime2(m, indices = ["a","b","c"], name = "g2")
assert g2.name == "g2"
assert g2.mesh == MeshReTime2(m1, m2)
assert g2.target_shape == [3,3]
assert (g2.data == np.zeros((50,100,3,3))).all()
