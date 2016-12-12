from pytriqs.gf.local._imp_tools import get_indices_in_dict

def init(mesh = None, shape = None, name = 'g', **kwargs):
    """
    """
    if mesh is None:
        from pytriqs.gf.local import MeshReTime
        from gf_retime import MeshReTimeReTime
        window = kwargs.pop('window')
        t_min = window[0]
        t_max = window[1]
        n_max = kwargs.pop('n_points',10000)
        m = MeshReTime(t_min, t_max, n_max)
        mesh = MeshReTimeReTime(m,m)

    indices_pack = get_indices_in_dict(kwargs)
    if shape is None:
      assert indices_pack, "No shape, no indices !"
      indicesL, indicesR = indices_pack
      shape = len(indicesL), len(indicesR)
    if kwargs: raise ValueError, "GfReTime_x_ReTime: Unused parameters %s were passed." % kwargs.keys()

    return (mesh, shape, indices_pack, name), {}


