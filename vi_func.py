# ##### BEGIN GPL LICENSE BLOCK #####
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software Foundation,
#  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# ##### END GPL LICENSE BLOCK #####

import bpy, os, sys, inspect, multiprocessing, mathutils, bmesh, datetime, colorsys, blf, bpy_extras, math
from subprocess import Popen
from numpy import array, digitize, amax, amin, average, clip, char, int8, int16, frombuffer, uint8, multiply, float32, zeros
from math import sin, cos, asin, acos, pi, ceil, log10
from math import e as expo
from mathutils import Vector, Matrix
from mathutils.bvhtree import BVHTree
from xml.dom import minidom
from bpy.props import IntProperty, StringProperty, EnumProperty, FloatProperty, BoolProperty, FloatVectorProperty
checked_groups_names_list = []
materials_from_group = set()
dtdf = datetime.date.fromordinal

try:
    from netgen import occ
    from netgen.meshing import MeshingStep
    ng = 1
except Exception:
    ng = 0


def li_calcob(ob, li):
    ovp = ob.vi_params

    if not ob.data.materials:
        ovp.vi_type_string = ''
    else:
        ovp.vi_type_string = 'LiVi Calc' if [face.index for face in ob.data.polygons if face.material_index < len(ob.material_slots) and ob.material_slots[face.material_index].material and ob.material_slots[face.material_index].material and ob.material_slots[face.material_index].material.vi_params.mattype == '1'] else ''

    return ovp.vi_type_string == 'LiVi Calc'


def ec_update(self, context):
    self['ecentries'] = []

    if not self.ee.updated:
        self.ee.update()

    if self.embodiedclass in self.ee.propdict.keys():
        if self.embodiedtype not in self.ee.propdict[self.embodiedclass]:
            self.embodiedtype = list(self.ee.propdict[self.embodiedclass].keys())[0]

        if self.embodiedmat not in self.ee.propdict[self.embodiedclass][self.embodiedtype]:
            self.embodiedmat = list(self.ee.propdict[self.embodiedclass][self.embodiedtype])[0]

        self['ecdict'] = self.ee.propdict[self.embodiedclass][self.embodiedtype][self.embodiedmat]
        self['ecentries'] = [(k, self.ee.propdict[self.embodiedclass][self.embodiedtype][self.embodiedmat][k]) for k in self['ecdict'].keys()]


def ret_datab(fname, r_w):
    # Need the following to initiate a scene held database
    db = ''

    if bpy.context.scene.vi_params.get('viparams'):
        db = bpy.context.scene.vi_params['viparams'].get('datab')

    addondir = os.path.basename(os.path.dirname(os.path.abspath(__file__)))
    addonfolder = os.path.dirname(os.path.abspath(__file__))
    vi_prefs = bpy.context.preferences.addons['{}'.format(addondir)].preferences

    if vi_prefs and vi_prefs.datab and os.path.isdir(vi_prefs.datab):
        if os.path.isfile(os.path.join(vi_prefs.datab, fname)) or r_w == 'w':
            datab = os.path.join(vi_prefs.datab, fname)
        else:
            datab = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'EPFiles', fname)

    elif db == os.path.join(addonfolder, 'EPFiles', fname):
        datab = db

    else:
        datab = os.path.join(addonfolder, 'EPFiles', fname)
        logentry('No database directory set but the scene has one built-in. Using the scene directory')

    return datab


def ret_plt():
    try:
        import matplotlib
        matplotlib.use('qtagg', force=True)
        from matplotlib import pyplot as plt
        plt.figure()
        return plt

    except Exception as e:
        logentry('Matplotlib error: {}'.format(e))
        return 0


def ret_mcm():
    try:
        import matplotlib.cm as mcm
        return mcm
    except Exception as e:
        logentry('Matplotlib error: {}'.format(e))
        return 0


def create_coll(c, name):
    if bpy.data.collections.get(name):
        coll = bpy.data.collections[name]
    else:
        coll = bpy.data.collections.new(name)
        c.scene.collection.children.link(coll)

    for lcc in c.view_layer.layer_collection.children:
        if lcc.name == name:
            c.view_layer.active_layer_collection = lcc

    return coll


def create_empty_coll(c, name):
    coll = create_coll(c, name)

    for o in coll.objects:
        if name == 'LiVi Results' and o.vi_params.vi_type_string == 'LiVi Res':
            bpy.data.objects.remove(o)

    return coll


def move_to_coll(context, coll, o):
    if o.parent:
        o.parent = None
    collection = create_coll(context, coll)

    if o.name not in collection.objects:
        collection.objects.link(o)
        for c in bpy.data.collections:
            if c.name != coll and o.name in c.objects:
                c.objects.unlink(o)
        if o.name in context.scene.collection.objects:
            context.scene.collection.objects.unlink(o)


def clear_coll(c, coll):
    for o in coll.objects:
        if coll.name == 'LiVi Results' and o.vi_params.vi_type_string != 'LiVi Res':
            c.scene.collection.link(o)
            coll.objects.unlink(o)
        elif coll.name == 'LiVi Results':
            coll.objects.unlink(o)
            bpy.data.objects.remove(o)

        elif coll.name == 'FloVi Mesh' and o.vi_params.vi_type_string != 'FloVi Mesh':
            c.scene.collection.objects.link(o)
            coll.objects.unlink(o)
        elif coll.name == 'FloVi Mesh':
            coll.objects.unlink(o)
            bpy.data.objects.remove(o)


def rm_coll(c, colls):
    for coll in colls:
        print(coll.name)
        for obj in coll.objects:
            bpy.data.objects.remove(obj, do_unlink=True)

        bpy.data.collections.remove(coll)


def ret_coll_bb(coll):
    for oi, ob in enumerate(coll.objects):
        if ob.type == 'MESH':
            gbbs = [ob.matrix_world@mathutils.Vector(bb) for bb in ob.bound_box]
            lbbs = list(zip(*gbbs))

            if not oi:
                (minx, miny, minz) = [min(lbb) for lbb in lbbs]
                (maxx, maxy, maxz) = [max(lbb) for lbb in lbbs]

            else:
                mins = [min(lbb) for lbb in lbbs]
                maxs = [max(lbb) for lbb in lbbs]
                minx = mins[0] if mins[0] < minx else minx
                miny = mins[1] if mins[1] < miny else miny
                minz = mins[2] if mins[2] < minz else minz
                maxx = maxs[0] if maxs[0] > maxx else maxx
                maxy = maxs[1] if maxs[1] > maxy else maxy
                maxz = maxs[2] if maxs[2] > maxz else maxz

    minx, miny, minz, maxx, maxy, maxz = minx - 0.1, miny - 0.1, minz - 0.1, maxx + 0.1, maxy + 0.1, maxz + 0.1
    bb_mesh = bpy.data.meshes.new("bb_mesh")
    bb_v_cos = [[minx, miny, minz], [maxx, miny, minz], [maxx, maxy, minz], [minx, maxy, minz],
                [minx, miny, maxz], [maxx, miny, maxz], [maxx, maxy, maxz], [minx, maxy, maxz]]
    bb_f_is = [[0, 3, 2, 1], [0, 1, 5, 4], [0, 4, 7, 3], [2, 3, 7, 6], [1, 2, 6, 5], [4, 5, 6, 7]]
    bb_mesh.from_pydata(bb_v_cos, [], bb_f_is)
    return bb_mesh


CIE_X = (1.299000e-04, 2.321000e-04, 4.149000e-04, 7.416000e-04, 1.368000e-03,
         2.236000e-03, 4.243000e-03, 7.650000e-03, 1.431000e-02, 2.319000e-02,
         4.351000e-02, 7.763000e-02, 1.343800e-01, 2.147700e-01, 2.839000e-01,
         3.285000e-01, 3.482800e-01, 3.480600e-01, 3.362000e-01, 3.187000e-01,
         2.908000e-01, 2.511000e-01, 1.953600e-01, 1.421000e-01, 9.564000e-02,
         5.795001e-02, 3.201000e-02, 1.470000e-02, 4.900000e-03, 2.400000e-03,
         9.300000e-03, 2.910000e-02, 6.327000e-02, 1.096000e-01, 1.655000e-01,
         2.257499e-01, 2.904000e-01, 3.597000e-01, 4.334499e-01, 5.120501e-01,
         5.945000e-01, 6.784000e-01, 7.621000e-01, 8.425000e-01, 9.163000e-01,
         9.786000e-01, 1.026300e+00, 1.056700e+00, 1.062200e+00, 1.045600e+00,
         1.002600e+00, 9.384000e-01, 8.544499e-01, 7.514000e-01, 6.424000e-01,
         5.419000e-01, 4.479000e-01, 3.608000e-01, 2.835000e-01, 2.187000e-01,
         1.649000e-01, 1.212000e-01, 8.740000e-02, 6.360000e-02, 4.677000e-02,
         3.290000e-02, 2.270000e-02, 1.584000e-02, 1.135916e-02, 8.110916e-03,
         5.790346e-03, 4.106457e-03, 2.899327e-03, 2.049190e-03, 1.439971e-03,
         9.999493e-04, 6.900786e-04, 4.760213e-04, 3.323011e-04, 2.348261e-04,
         1.661505e-04, 1.174130e-04, 8.307527e-05, 5.870652e-05, 4.150994e-05,
         2.935326e-05, 2.067383e-05, 1.455977e-05, 1.025398e-05, 7.221456e-06,
         5.085868e-06, 3.581652e-06, 2.522525e-06, 1.776509e-06, 1.251141e-06)

CIE_Y = (3.917000e-06, 6.965000e-06, 1.239000e-05, 2.202000e-05, 3.900000e-05,
         6.400000e-05, 1.200000e-04, 2.170000e-04, 3.960000e-04, 6.400000e-04,
         1.210000e-03, 2.180000e-03, 4.000000e-03, 7.300000e-03, 1.160000e-02,
         1.684000e-02, 2.300000e-02, 2.980000e-02, 3.800000e-02, 4.800000e-02,
         6.000000e-02, 7.390000e-02, 9.098000e-02, 1.126000e-01, 1.390200e-01,
         1.693000e-01, 2.080200e-01, 2.586000e-01, 3.230000e-01, 4.073000e-01,
         5.030000e-01, 6.082000e-01, 7.100000e-01, 7.932000e-01, 8.620000e-01,
         9.148501e-01, 9.540000e-01, 9.803000e-01, 9.949501e-01, 1.000000e+00,
         9.950000e-01, 9.786000e-01, 9.520000e-01, 9.154000e-01, 8.700000e-01,
         8.163000e-01, 7.570000e-01, 6.949000e-01, 6.310000e-01, 5.668000e-01,
         5.030000e-01, 4.412000e-01, 3.810000e-01, 3.210000e-01, 2.650000e-01,
         2.170000e-01, 1.750000e-01, 1.382000e-01, 1.070000e-01, 8.160000e-02,
         6.100000e-02, 4.458000e-02, 3.200000e-02, 2.320000e-02, 1.700000e-02,
         1.192000e-02, 8.210000e-03, 5.723000e-03, 4.102000e-03, 2.929000e-03,
         2.091000e-03, 1.484000e-03, 1.047000e-03, 7.400000e-04, 5.200000e-04,
         3.611000e-04, 2.492000e-04, 1.719000e-04, 1.200000e-04, 8.480000e-05,
         6.000000e-05, 4.240000e-05, 3.000000e-05, 2.120000e-05, 1.499000e-05,
         1.060000e-05, 7.465700e-06, 5.257800e-06, 3.702900e-06, 2.607800e-06,
         1.836600e-06, 1.293400e-06, 9.109300e-07, 6.415300e-07, 4.518100e-07)

CIE_Z = (6.061000e-04, 1.086000e-03, 1.946000e-03, 3.486000e-03, 6.450001e-03,
         1.054999e-02, 2.005001e-02, 3.621000e-02, 6.785001e-02, 1.102000e-01,
         2.074000e-01, 3.713000e-01, 6.456000e-01, 1.039050e+00, 1.385600e+00,
         1.622960e+00, 1.747060e+00, 1.782600e+00, 1.772110e+00, 1.744100e+00,
         1.669200e+00, 1.528100e+00, 1.287640e+00, 1.041900e+00, 8.129501e-01,
         6.162000e-01, 4.651800e-01, 3.533000e-01, 2.720000e-01, 2.123000e-01,
         1.582000e-01, 1.117000e-01, 7.824999e-02, 5.725001e-02, 4.216000e-02,
         2.984000e-02, 2.030000e-02, 1.340000e-02, 8.749999e-03, 5.749999e-03,
         3.900000e-03, 2.749999e-03, 2.100000e-03, 1.800000e-03, 1.650001e-03,
         1.400000e-03, 1.100000e-03, 1.000000e-03, 8.000000e-04, 6.000000e-04,
         3.400000e-04, 2.400000e-04, 1.900000e-04, 1.000000e-04, 4.999999e-05,
         3.000000e-05, 2.000000e-05, 1.000000e-05, 0.000000e+00, 0.000000e+00,
         0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
         0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
         0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
         0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
         0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
         0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
         0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00)


XYZtosRGB = (3.2404542, -1.5371385, -0.4985314,
             -0.9692660, 1.8760108, 0.0415560,
             0.0556434, -0.2040259, 1.0572252)


def planck(w, t):
    wm = w * 1e-9
    c1 = 3.7402e-16
    c2 = 1.43848e-2
    return (c1 * wm**-5.0) / (expo**(c2 / (wm * t)) - 1)


def ct2RGB(t):
    xyz = [0.0, 0.0, 0.0]
    rgb = [0.0, 0.0, 0.0]
    step = 5  # unit: nanometer
    startWavelength = 360
    endWavelength = 830
    i = 0

    for w in range(startWavelength, endWavelength + step, step):
        Irr = planck(w, t)
        xyz[0] += Irr * CIE_X[i]
        xyz[1] += Irr * CIE_Y[i]
        xyz[2] += Irr * CIE_Z[i]
        i += 1
    # mxyz = max(xyz)
    mxyz = xyz[0] + xyz[1] + xyz[2]
    xyz[0] /= mxyz
    xyz[1] /= mxyz
    xyz[2] /= mxyz
    rgb[0] = max(0, xyz[0] * XYZtosRGB[0] + xyz[1] * XYZtosRGB[1] + xyz[2] * XYZtosRGB[2])
    rgb[1] = max(0, xyz[0] * XYZtosRGB[3] + xyz[1] * XYZtosRGB[4] + xyz[2] * XYZtosRGB[5])
    rgb[2] = max(0, xyz[0] * XYZtosRGB[6] + xyz[1] * XYZtosRGB[7] + xyz[2] * XYZtosRGB[8])
    maxrgb = max(rgb)

    if maxrgb > 0.0:
        rgb[0] /= maxrgb
        rgb[1] /= maxrgb
        rgb[2] /= maxrgb
    return rgb


def retcols(cmap, levels):
    try:
        rgbas = [cmap(int(i * 255 / (levels - 1))) for i in range(levels)]
    except Exception:
        hs = [0.75 - 0.75 * (i / (levels - 1)) for i in range(levels)]
        rgbas = [(*colorsys.hsv_to_rgb(h, 1.0, 1.0), 1.0) for h in hs]
    return rgbas


def move_obs(s_coll, coll, t_string):
    for ob in coll.objects:
        if ob.vi_params.vi_type_string != t_string:
            s_coll.objects.link(ob)
            coll.objects.unlink(ob)


def cmap(svp):
    cols = retcols(ret_mcm().get_cmap(svp.vi_leg_col), svp.vi_leg_levels)
    cols = [[col[0], col[1], col[2], 1 - svp.vi_disp_trans] for col in cols]

    for i in range(svp.vi_leg_levels):
        matname = '{}#{}'.format('vi-suite', i)

        if not bpy.data.materials.get(matname):
            bpy.data.materials.new(matname)
            bpy.data.materials[matname].specular_intensity = 0
            bpy.data.materials[matname].specular_color = (0, 0, 0)

        bpy.data.materials[matname].diffuse_color = cols[i][0:4]
        bpy.data.materials[matname].use_nodes = True
        nodes = bpy.data.materials[matname].node_tree.nodes
        links = bpy.data.materials[matname].node_tree.links

        for node in nodes:
            nodes.remove(node)

        node_output = nodes.new(type='ShaderNodeOutputMaterial')
        node_output.location = 400, 0

        if svp.vi_disp_trans > 0:
            node_material = nodes.new(type='ShaderNodeBsdfTransparent')
            node_material.inputs[0].default_value[3] = 1 - svp.vi_disp_trans

        elif svp.vi_disp_mat:
            node_material = nodes.new(type='ShaderNodeEmission')
            node_material.inputs[1].default_value = svp.vi_disp_ems
            node_lp = nodes.new(type='ShaderNodeLightPath')
            node_mix = nodes.new(type='ShaderNodeMixShader')
            node_trans = nodes.new(type='ShaderNodeBsdfTransparent')
            links.new(node_lp.outputs[0], node_mix.inputs[0])
            links.new(node_material.outputs[0], node_mix.inputs[2])
            links.new(node_trans.outputs[0], node_mix.inputs[1])
            links.new(node_mix.outputs[0], node_output.inputs[0])
        else:
            # create diffuse node
            node_material = nodes.new(type='ShaderNodeBsdfDiffuse')
            node_material.inputs[1].default_value = 0.5

        node_material.inputs[0].default_value = (cols[i][0:4])  # green RGBA
        node_material.location = 0, 0

        # create output node
        if 1 - svp.vi_disp_trans == 1 and svp.vi_disp_mat:
            links.new(node_mix.outputs[0], node_output.inputs[0])
        else:
            links.new(node_material.outputs[0], node_output.inputs[0])


def regresults(scene, frames, simnode, res):
    for i, f in enumerate(frames):
        simnode['maxres'][str(f)] = amax(res[i])
        simnode['minres'][str(f)] = amin(res[i])
        simnode['avres'][str(f)] = average(res[i])

    scene.vi_leg_max, scene.vi_leg_min = max(simnode['maxres'].values()), min(simnode['minres'].values())


def clearlayers(bm, ltype):
    if ltype in ('a', 'f'):
        while bm.faces.layers.float:
            bm.faces.layers.float.remove(bm.faces.layers.float[0])
        while bm.verts.layers.float:
            bm.verts.layers.float.remove(bm.verts.layers.float[0])
    if ltype in ('a', 's'):
        while bm.faces.layers.string:
            bm.faces.layers.string.remove(bm.faces.layers.string[0])
        while bm.verts.layers.string:
            bm.verts.layers.string.remove(bm.verts.layers.string[0])
    if ltype in ('a', 'i'):
        while bm.faces.layers.int:
            bm.faces.layers.int.remove(bm.faces.layers.int[0])
        while bm.verts.layers.string:
            bm.verts.layers.int.remove(bm.verts.layers.int[0])


def retvpvloc(context):
    return bpy_extras.view3d_utils.region_2d_to_origin_3d(context.region, context.space_data.region_3d,
                                                          (context.region.width / 2.0, context.region.height / 2.0))


def rettree(scene, obs, ignore):
    bmob = bmesh.new()
    dp = bpy.context.evaluated_depsgraph_get()

    for soi, so in enumerate(obs):
        btemp = bpy.data.meshes.new("temp")
        bmtemp = bmesh.new()
        bmtemp.from_object(so, dp)
        bmtemp.transform(so.matrix_world)
        bmtemp.normal_update()
        delfaces = [face for face in bmtemp.faces if face.material_index >= len(so.data.materials) or so.data.materials[face.material_index].vi_params.mattype == ignore]
        bmesh.ops.delete(bmtemp, geom=delfaces, context='FACES')
        bmtemp.to_mesh(btemp)
        bmob.from_mesh(btemp)
        bpy.data.meshes.remove(btemp)

    tree = BVHTree.FromBMesh(bmob)
    bmob.free()
    bmtemp.free()
    return tree


def ret_empty_menu(self, context):
    e_names = [(o.name, o.name, 'Name of empty') for o in bpy.data.objects if o.type == 'EMPTY']

    if not e_names:
        e_names = [('None', 'None', 'None')]

    return e_names


class progressfile():
    def __init__(self, folder, starttime, calcsteps):
        self.starttime = starttime
        self.calcsteps = calcsteps
        self.curres = 0
        self.folder = folder
        self.pfile = os.path.join(folder, 'viprogress')

        with open(self.pfile, 'w') as pfile:
            pfile.write('STARTING')

    def check(self, curres):
        if curres == 'CANCELLED':
            with open(self.pfile, 'w') as pfile:
                pfile.write('CANCELLED')

        with open(self.pfile, 'r') as pfile:
            if 'CANCELLED' in pfile.read():
                return 'CANCELLED'

        with open(self.pfile, 'w') as pfile:
            if curres:
                if curres > self.curres:
                    self.ldt = datetime.datetime.now()
                    self.curres = curres

                dt = (self.ldt - self.starttime) / (curres / self.calcsteps) - (datetime.datetime.now() - self.starttime)
                pfile.write('{} {}'.format(int(100 * curres / self.calcsteps), datetime.timedelta(seconds=dt.seconds if dt.total_seconds() > 0 else 0)))
            else:
                pfile.write('0 Initialising')


class fvprogressfile():
    def __init__(self, folder):
        self.pfile = os.path.join(folder, 'viprogress')

        with open(self.pfile, 'w') as pfile:
            pfile.write('STARTING')

    def check(self, curres):
        if curres == 'CANCELLED':
            if 'CANCELLED' in self.pfile.read():
                return 'CANCELLED'

        with open(self.pfile, 'w') as pfile:
            if curres:
                pfile.write(curres)
            else:
                pfile.write('0 Initialising')


def qtprogressbar(file, pdll_path, calctype):
    #addon_path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    qttext = "# -*- coding: " + sys.getfilesystemencoding() + " -*-\n\
import os, sys\n\
if sys.platform == 'win32':\n\
    os.add_dll_directory(r'" + pdll_path + "')\n\
from PySide6.QtWidgets import (QApplication, QMainWindow, QPushButton, QVBoxLayout, QWidget, QProgressBar, QLabel)\n\
from PySide6.QtCore import QProcess, QTimer, Qt\n\
\n\
class MainWindow(QMainWindow):\n\
    def __init__(self):\n\
        super().__init__()\n\
        self.setWindowTitle('" + calctype + "')\n\
        self.setFixedWidth(350)\n\
        self.btn = QPushButton('Cancel')\n\
        #self.connect(self.update)\n\
        self.btn.pressed.connect(self.stop_process)\n\
        #self.text = QPlainTextEdit()\n\
        self.progress = QProgressBar()\n\
        self.progress.setTextVisible(False)\n\
        self.label = QLabel('Initialsing')\n\
        font = self.label.font()\n\
        font.setPointSize(12)\n\
        self.label.setFont(font)\n\
        self.label.setAlignment(Qt.AlignmentFlag.AlignHCenter)\n\
        self.progress.setStyleSheet('QProgressBar {border: 1px solid #1DACD6; border-radius: 2px; background-color: #E0E0E0;}'\n\
                                    'QProgressBar::chunk {background-color: #1DACD6; width: 13px; margin: 0.5px;}')\n\
\n\
        self.progress.setRange(0, 100)\n\
        self.timer = QTimer()\n\
        self.timer.timeout.connect(self.p_update)\n\
        self.timer.start(1000)\n\
\n\
        #self.text.setReadOnly(False)\n\
\n\
        l = QVBoxLayout()\n\
        l.addWidget(self.label)\n\
        l.addWidget(self.progress)\n\
        l.addWidget(self.btn)\n\
        \n\
        #l.addWidget(self.text)\n\
\n\
        w = QWidget()\n\
        w.setLayout(l)\n\
\n\
        self.setCentralWidget(w)\n\
    \n\
    def p_update(self):\n\
        with open(r'" + file + "', 'r') as pffile:\n\
            try:    (percent, tr) = pffile.readlines()[0].split()\n\
            except: percent, tr = 0, 'Not known'\n\
\n\
        self.progress.setValue(int(percent))\n\
        self.label.setText(f'{percent}% Complete - Time remaining: {tr}')\n\
        \n\
    def stop_process(self):\n\
        with open(r'" + file + "', 'w') as pffile:\n\
            pffile.write('CANCELLED')\n\
        self.close()\n\
\n\
#QApplication.shutdown()\n\
#app = QApplication.instance()\n\
try:\n\
    app = QApplication(sys.argv)\n\
except:\n\
    app = QApplication.instance()\n\
\n\
w = MainWindow()\n\
w.show()\n\
\n\
app.exec()\n\
QApplication.shutdown(app)"

    with open(file + "qt.py", 'w') as qtfile:
        qtfile.write(qttext)
    return Popen([sys.executable, file + "qt.py"])


def progressbar(file, calctype):
    addonpath = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    kivytext = "# -*- coding: " + sys.getfilesystemencoding() + " -*-\n\
import os, sys\n\
sys.path.append(os.path.join(r'" + addonpath + "', 'Python', sys.platform))\n\
from kivy.app import App \n\
from kivy.clock import Clock \n\
from kivy.uix.progressbar import ProgressBar\n\
from kivy.uix.boxlayout import BoxLayout\n\
from kivy.uix.button import Button\n\
from kivy.uix.label import Label\n\
from kivy.config import Config\n\
Config.set('graphics', 'width', '500')\n\
Config.set('graphics', 'height', '200')\n\
Config.set('kivy', 'log_level', 'warning')\n\
\n\
class CancelButton(Button):\n\
    def on_touch_down(self, touch):\n\
        if 'button' in touch.profile:\n\
            if self.collide_point(*touch.pos):\n\
                with open(r'" + file + "', 'w') as pffile:\n\
                    pffile.write('CANCELLED')\n\
                App.get_running_app().stop()\n\
        else:\n\
            return\n\
    def on_open(self, widget, parent):\n\
        self.focus = True\n\
\n\
class Calculating(App):\n\
    bl = BoxLayout(orientation='vertical')\n\
    rpb = ProgressBar()\n\
    label = Label(text=' 0% Complete', font_size=20)\n\
    button = CancelButton(text='Cancel', font_size=20)\n\
    bl.add_widget(rpb)\n\
    bl.add_widget(label)\n\
    bl.add_widget(button)\n\
\n\
    def build(self):\n\
        self.title = 'Calculating " + calctype + "'\n\
        refresh_time = 1\n\
        Clock.schedule_interval(self.timer, refresh_time)\n\
        return self.bl\n\
\n\
    def timer(self, dt):\n\
        with open(r'" + file + "', 'r') as pffile:\n\
            try:    (percent, tr) = pffile.readlines()[0].split()\n\
            except: percent, tr = 0, 'Not known'\n\
        self.rpb.value = int(percent)\n\
        self.label.text = '{}% Complete - Time remaining: {}'.format(percent, tr)\n\
\n\
if __name__ == '__main__':\n\
    Calculating().run()\n"

    with open(file + ".py", 'w') as kivyfile:
        kivyfile.write(kivytext)
    return Popen([sys.executable, file + ".py"])


def cancel_window(file, pdll_path, calc_type):
    qttext = "# -*- coding: " + sys.getfilesystemencoding() + " -*-\n\
import os, sys\n\
if sys.platform == 'win32':\n\
    os.add_dll_directory(r'" + pdll_path + "')\n\
from PySide6.QtWidgets import (QApplication, QMainWindow, QPushButton, QPlainTextEdit, QVBoxLayout, QWidget)\n\
from PySide6.QtCore import QProcess, QTimer\n\
\n\
class MainWindow(QMainWindow):\n\
    def __init__(self):\n\
        super().__init__()\n\
        self.setWindowTitle('{}')\n\
        self.setFixedWidth(300)\n\
        self.btn = QPushButton('Cancel')\n\
        self.btn.pressed.connect(self.stop_process)\n\
        l = QVBoxLayout()\n\
        l.addWidget(self.btn)\n\
        w = QWidget()\n\
        w.setLayout(l)\n\
        self.setCentralWidget(w)\n\
\n\
    def stop_process(self):\n\
        self.close()\n\
\n\
try:\n\
    app = QApplication(sys.argv)\n\
except:\n\
    app = QApplication.instance()\n\
\n\
w = MainWindow()\n\
w.show()\n\
\n\
app.exec()\n\
QApplication.shutdown(app)".format(calc_type)

    with open(file + ".py", 'w') as qtfile:
        qtfile.write(qttext)
    return Popen([sys.executable, file + ".py"])


def qtfvprogress(file, pdll_path, et, residuals, frame):
    #addonpath = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    qttext = "# -*- coding: " + sys.getfilesystemencoding() + " -*-\n\
import os, sys\n\
if sys.platform == 'win32':\n\
    os.add_dll_directory(r'" + pdll_path + "')\n\
from PySide6.QtWidgets import (QApplication, QMainWindow, QPushButton, QPlainTextEdit, QVBoxLayout, QHBoxLayout, QWidget, QProgressBar, QLabel)\n\
from PySide6.QtCore import QProcess, QTimer, Qt\n\
r_dict = {'epsilon': 'ep', 'Time': 'Ti', 'k': 'k', 'Ux': 'Ux', 'Uy': 'Uy', 'Uz': 'Uz', 'p': 'p', 'p_rgh': 'pr', 'e': 'en', 'T': 'T'}\n\
\n\
class MainWindow(QMainWindow):\n\
    def __init__(self):\n\
        super().__init__()\n\
        self.rpbs = []\n\
        self.nums = []\n\
        self.setWindowTitle('OpenFOAM Residuals')\n\
        self.setFixedWidth(370)\n\
        self.btn = QPushButton('Cancel')\n\
        self.btn.pressed.connect(self.stop_process)\n\
        l = QVBoxLayout()\n\
        residual = QHBoxLayout()\n\
        progress = QProgressBar()\n\
        progress.setTextVisible(False)\n\
        label = QLabel('Ti')\n\
        label.setFixedWidth(20)\n\
        num = QLabel('0')\n\
        num.setFixedWidth(60)\n\
        num.setAlignment(Qt.AlignmentFlag.AlignLeft)\n\
        self.nums.append(num)\n\
        progress.setRange(0, " + str(et) + ")\n\
        progress.setStyleSheet('QProgressBar {border: 1px solid #1DACD6; border-radius: 2px; background-color: #E0E0E0;}'\n\
                                    'QProgressBar::chunk {background-color: #1DACD6; width: 5px; margin: 0.5px;}')\n\
        self.rpbs.append(progress)\n\
        residual.addWidget(label)\n\
        residual.addWidget(progress)\n\
        residual.addWidget(num)\n\
        l.addLayout(residual)\n\
\n\
        for r in " + residuals + ":\n\
            residual = QHBoxLayout()\n\
            label = QLabel(r_dict[r])\n\
            label.setFixedWidth(20)\n\
            num = QLabel('1.00000')\n\
            num.setFixedWidth(60)\n\
            num.setAlignment(Qt.AlignmentFlag.AlignLeft)\n\
            self.nums.append(num)\n\
            residual.addWidget(label)\n\
            progress = QProgressBar()\n\
            progress.setTextVisible(False)\n\
            progress.setRange(0, 100)\n\
            progress.setStyleSheet('QProgressBar {border: 1px solid #1DACD6; border-radius: 2px; background-color: #E0E0E0; margin-left: 0em; margin-right: 0em;}'\n\
                                        'QProgressBar::chunk {background-color: #1DACD6; width: 5px; margin: 0.5px;}')\n\
            self.rpbs.append(progress)\n\
            residual.addWidget(progress)\n\
            residual.addWidget(num)\n\
            l.addLayout(residual)\n\
\n\
        l.addWidget(self.btn)\n\
        w = QWidget()\n\
        w.setLayout(l)\n\
        self.setCentralWidget(w)\n\
        self.timer = QTimer()\n\
        self.timer.timeout.connect(self.p_update)\n\
        self.timer.start(1000)\n\
\n\
    def p_update(self):\n\
        with open(r'" + file + "', 'r') as pffile:\n\
            for ri, r in enumerate(pffile.readlines()):\n\
                if r.split()[0] in " + residuals + ":\n\
                    try:\n\
                        li = " + residuals + ".index(r.split()[0])\n\
                        self.rpbs[li+1].setValue(100 * abs(float(r.split()[1]))**0.5)\n\
                        self.nums[li+1].setText('{:.5f}'.format(abs(float(r.split()[1]))))\n\
                    except Exception as e:\n\
                        print(e)\n\
                elif r.split()[0] == 'Time':\n\
                    self.rpbs[0].setValue(float(r.split()[1]))\n\
                    self.nums[0].setText('{:.1f}'.format(float(r.split()[1])))\n\
\n\
    def stop_process(self):\n\
        self.close()\n\
\n\
try:\n\
    app = QApplication(sys.argv)\n\
except:\n\
    app = QApplication.instance()\n\
\n\
w = MainWindow()\n\
w.show()\n\
\n\
app.exec()\n\
QApplication.shutdown(app)"

    with open(file + "qt.py", 'w') as qtfile:
        qtfile.write(qttext)

    return Popen([sys.executable, file + "qt.py"])


def logentry(text):
    log = bpy.data.texts.new('vi-suite-log') if 'vi-suite-log' not in bpy.data.texts else bpy.data.texts['vi-suite-log']
    log.write('')
    log.write('{}: {}\n'.format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), text))


def chunks(li, n):
    for v in range(0, len(li), n):
        yield li[v:v + n]


# This function can be used to modify results with a driver function
def ret_res_vals(svp, reslist):
    if svp.vi_res_process == '2' and svp.script_file:
        try:
            if svp.vi_leg_levels == len(bpy.app.driver_namespace['restext']()) + 1:
                return bpy.app.driver_namespace['resmod'](reslist)
            else:
                logentry('Set legend levels to one more than the length of the restext return list')
                return reslist
        except Exception as e:
            logentry('User script error {}. Check console'.format(e))
            return reslist

    elif svp.vi_res_process == '1' and svp.vi_res_mod:
        try:
            return array([eval('{}{}'.format(r, svp.vi_res_mod)) for r in reslist])
        except Exception:
            return reslist
    else:
        return reslist


def lividisplay(self, scene):
    svp = scene.vi_params
    svp.vi_view_update = 1

    if self.id_data.vi_params.vi_type_string == 'LiVi Res':
        frames = range(svp['liparams']['fs'], svp['liparams']['fe'] + 1)
        ll = svp.vi_leg_levels
        increment = 1 / ll

        if len(frames) > 1:
            if not self.id_data.data.animation_data:
                self.id_data.data.animation_data_create()

            self.id_data.data.animation_data.action = bpy.data.actions.new(name="LiVi {} MI".format(self.name))
            fis = [str(face.index) for face in self.id_data.data.polygons]
            lms = {fi: self.id_data.data.animation_data.action.fcurves.new(data_path='polygons[{}].material_index'.format(fi)) for fi in fis}

            for fi in fis:
                lms[fi].keyframe_points.add(len(frames))

        for f, frame in enumerate(frames):
            bm = bmesh.new()
            bm.from_mesh(self.id_data.data)
            geom = bm.verts if svp['liparams']['cp'] == '1' else bm.faces
            sf = str(frame)

            if svp.li_disp_menu == 'aga1v':
                res_name = f'aga{svp.vi_views}v{frame}'
            elif svp.li_disp_menu == 'ago1v':
                res_name = f'ago{svp.vi_views}v{frame}'
            elif svp.li_disp_menu == 'rt':
                res_name = f'{svp.au_sources}_rt{frame}'
            else:
                res_name = f'{svp.li_disp_menu}{frame}'

            if geom.layers.float.get(res_name):
                livires = geom.layers.float[res_name]
                res = geom.layers.float[res_name]
                oreslist = [g[livires] for g in geom]
                self['omax'][sf], self['omin'][sf], self['oave'][sf] = max(oreslist), min(oreslist), sum(oreslist) / len(oreslist)

                if svp.vi_res_process == '2' and svp.script_file:
                    try:
                        vals = ret_res_vals(svp, array([f[livires] for f in bm.faces])) if svp['liparams']['cp'] == '0' else \
                            ret_res_vals(svp, array([(sum([vert[livires] for vert in f.verts])) / len(f.verts) for f in bm.faces]))
                        nmatis = array(vals).astype(int8)
                    except Exception:
                        nmatis = zeros(len(bm.faces))

                else:
                    smaxres, sminres = ret_res_vals(svp, [max(svp['liparams']['maxres'].values()), min(svp['liparams']['minres'].values())])[:2]

                    if smaxres > sminres:
                        vals = (array(ret_res_vals(svp, array([f[livires] for f in bm.faces]))) - sminres) / (smaxres - sminres) if svp['liparams']['cp'] == '0' else \
                            (array(ret_res_vals(svp, array([(sum([vert[livires] for vert in f.verts])) / len(f.verts) for f in bm.faces]))) - sminres) / (smaxres - sminres)
                    else:
                        vals = array(ret_res_vals(svp, array([max(svp['liparams']['maxres'].values()) for x in range(len(bm.faces))])))

                    if livires != res:
                        for g in geom:
                            g[res] = g[livires]

                    if svp['liparams']['unit'] == 'SVF (%)X':
                        nmatis = [(0, ll - 1)[v == 1] for v in vals]
                    else:
                        bins = array([increment * i for i in range(ll)])
                        nmatis = clip(digitize(vals, bins, right=True) - 1, 0, ll - 1, out=None)

                bm.to_mesh(self.id_data.data)
                bm.free()

                if len(frames) == 1:
                    self.id_data.data.polygons.foreach_set('material_index', nmatis)
                elif len(frames) > 1:
                    for fii, fi in enumerate(fis):
                        lms[fi].keyframe_points[f].co = frame, nmatis[fii]


# def lividisplay(self, scene):
#     svp = scene.vi_params
#     svp.vi_view_update = 1
#     anim_name = "LiVi {} MI".format(self.name)
#     disp_menu = svp.li_disp_menu
#     params = 'liparams'
#     type_strings = ('LiVi Res', 'LiVi_Calc')
#
#     if self.id_data.vi_params.vi_type_string == type_strings[0]:
#         frames = range(svp[params]['fs'], svp[params]['fe'] + 1)
#         ll = svp.vi_leg_levels
#         increment = 1 / ll
#
#         if len(frames) > 1:
#             if not self.id_data.data.animation_data:
#                 self.id_data.data.animation_data_create()
#
#             self.id_data.data.animation_data.action = bpy.data.actions.new(name=anim_name)
#             fis = [str(face.index) for face in self.id_data.data.polygons]
#             lms = {fi: self.id_data.data.animation_data.action.fcurves.new(data_path='polygons[{}].material_index'.format(fi)) for fi in fis}
#
#             for fi in fis:
#                 lms[fi].keyframe_points.add(len(frames))
#
#         for f, frame in enumerate(frames):
#             bm = bmesh.new()
#             bm.from_mesh(self.id_data.data)
#             geom = bm.verts if svp[params]['cp'] == '1' else bm.faces
#             sf = str(frame)
#
#             if disp_menu == 'aga1v':
#                 res_name = f'aga{svp.vi_views}v{frame}'
#             elif disp_menu == 'ago1v':
#                 res_name = f'ago{svp.vi_views}v{frame}'
#             elif disp_menu == 'rt':
#                 res_name = f'{svp.au_sources}_rt{frame}'
#             else:
#                 res_name = f'{disp_menu}{frame}'
#
#             if geom.layers.float.get(res_name):
#                 vires = geom.layers.float[res_name]
#                 res = geom.layers.float[res_name]
#                 oreslist = [g[vires] for g in geom]
#                 self['omax'][sf], self['omin'][sf], self['oave'][sf] = max(oreslist), min(oreslist), sum(oreslist) / len(oreslist)
#
#                 if svp.vi_res_process == '2' and svp.script_file:
#                     try:
#                         vals = ret_res_vals(svp, array([f[vires] for f in bm.faces])) if svp[params]['cp'] == '0' else ret_res_vals(svp, array([(sum([vert[vires] for vert in f.verts])) / len(f.verts) for f in bm.faces]))
#                         nmatis = array(vals).astype(int8)
#                     except Exception:
#                         nmatis = zeros(len(bm.faces))
#
#                 else:
#                     smaxres, sminres = ret_res_vals(svp, [max(svp[params]['maxres'].values()), min(svp[params]['minres'].values())])[:2]
#
#                     if smaxres > sminres:
#                         vals = (array(ret_res_vals(svp, array([f[vires] for f in bm.faces]))) - sminres) / (smaxres - sminres) if svp[params]['cp'] == '0' else \
#                             (array(ret_res_vals(svp, array([(sum([vert[vires] for vert in f.verts])) / len(f.verts) for f in bm.faces]))) - sminres) / (smaxres - sminres)
#                     else:
#                         vals = array(ret_res_vals(svp, array([max(svp[params]['maxres'].values()) for x in range(len(bm.faces))])))
#
#                     if vires != res:
#                         for g in geom:
#                             g[res] = g[vires]
#
#                     if svp[params]['unit'] == 'SVF (%)X':
#                         nmatis = [(0, ll - 1)[v == 1] for v in vals]
#                     else:
#                         bins = array([increment * i for i in range(ll)])
#                         nmatis = clip(digitize(vals, bins, right=True) - 1, 0, ll - 1, out=None)
#
#                 bm.to_mesh(self.id_data.data)
#                 bm.free()
#
#                 if len(frames) == 1:
#                     self.id_data.data.polygons.foreach_set('material_index', nmatis)
#
#                 elif len(frames) > 1:
#                     for fii, fi in enumerate(fis):
#                         lms[fi].keyframe_points[f].co = frame, nmatis[fii]


def ret_vp_loc(context):
    cspr = context.space_data.region_3d
    cr = context.region
    return bpy_extras.view3d_utils.region_2d_to_origin_3d(cr, cspr, (cr.width / 2.0, cr.height / 2.0))


def viparams(op, scene):
    svp = scene.vi_params
    bdfp = bpy.data.filepath

    if not bdfp:
        op.report({'ERROR'}, "The Blender file has not been saved. Save the Blender file before exporting")
        return 'Save file'

    #isascii = lambda s: len(s) == len(s.encode())

    if len(bdfp) != len(bdfp.encode()):
        op.report({'WARNING'}, "The directory path or Blender filename has non-ascii characters in it. Photon mapping may not work")

    fd, fn = os.path.dirname(bpy.data.filepath), os.path.splitext(os.path.basename(bpy.data.filepath))[0]
    if not os.path.isdir(os.path.join(fd, fn)):
        os.makedirs(os.path.join(fd, fn))
    if not os.path.isdir(os.path.join(fd, fn, 'obj')):
        os.makedirs(os.path.join(fd, fn, 'obj'))
    if not os.path.isdir(os.path.join(fd, fn, 'bsdfs')):
        os.makedirs(os.path.join(fd, fn, 'bsdfs'))
    if not os.path.isdir(os.path.join(fd, fn, 'images')):
        os.makedirs(os.path.join(fd, fn, 'images'))
    if not os.path.isdir(os.path.join(fd, fn, 'lights')):
        os.makedirs(os.path.join(fd, fn, 'lights'))
    if not os.path.isdir(os.path.join(fd, fn, 'textures')):
        os.makedirs(os.path.join(fd, fn, 'textures'))
    if not os.path.isdir(os.path.join(fd, fn, 'Openfoam')):
        os.makedirs(os.path.join(fd, fn, 'Openfoam'))

    nd = os.path.join(fd, fn)
    fb, ofb, lfb, tfb, offb, idf = os.path.join(nd, fn), os.path.join(nd, 'obj'), os.path.join(nd, 'lights'), os.path.join(nd, 'textures'), os.path.join(nd, 'Openfoam'), os.path.join(nd, 'in.idf')

    if not svp.get('viparams'):
        svp['viparams'] = {}

    svp['viparams']['cat'] = ('cat ', 'type ')[str(sys.platform) == 'win32']
    svp['viparams']['nproc'] = str(multiprocessing.cpu_count())
    svp['viparams']['wnproc'] = str(multiprocessing.cpu_count()) if str(sys.platform) != 'win32' else '1'
    svp['viparams']['addonpath'] = os.path.dirname(os.path.abspath(__file__))
    svp['viparams']['filepath'] = bpy.data.filepath
    svp['viparams']['filename'] = fn
    svp['viparams']['filedir'] = fd
    svp['viparams']['newdir'] = nd
    svp['viparams']['filebase'] = fb
    svp['viparams']['drivers'] = []

    if not svp.get('spparams'):
        svp['spparams'] = {}

    if not svp.get('liparams'):
        svp['liparams'] = {}
        svp['liparams']['disp_count'] = 0

    svp['liparams']['objfilebase'] = ofb
    svp['liparams']['lightfilebase'] = lfb
    svp['liparams']['texfilebase'] = tfb
    svp['liparams']['livig'] = []
    svp['liparams']['livic'] = []
    svp['liparams']['livir'] = []
    svp['liparams']['shadc'] = []

    if not svp.get('enparams'):
        svp['enparams'] = {}

    svp['enparams']['idf_file'] = idf
    svp['enparams']['epversion'] = '25.1'

    if not svp.get('flparams'):
        svp['flparams'] = {}

    svp['flparams']['offilebase'] = offb


def nodestate(self, opstate):
    if self['exportstate'] != opstate:
        self.exported = False
        if self.bl_label[0] != '*':
            self.bl_label = '*' + self.bl_label
    else:
        self.exported = True
        if self.bl_label[0] == '*':
            self.bl_label = self.bl_label[1:-1]


def face_centre(ob, obresnum, f):
    if obresnum:
        vsum = mathutils.Vector((0, 0, 0))
        for v in f.vertices:
            vsum = ob.active_shape_key.data[v].co + vsum
        return (vsum / len(f.vertices))
    else:
        return (f.center)


def v_pos(ob, v):
    return (ob.active_shape_key.data[v].co if ob.name in bpy.context.scene.vi_params['liparams']['livir'] else ob.data.vertices[v].co)


def newrow(layout, s1, root, s2):
    row = layout.row()
    row.label(text=s1)
    row.prop(root, s2)


def newrow2(row, s1, root, s2):
    row.label(text=s1)
    row.prop(root, s2)


def retobj(name, fr, node, scene):
    if node.animmenu == "Geometry":
        return (os.path.join(scene['liparams']['objfilebase'], "{}-{}.obj".format(name.replace(" ", "_"), fr)))
    else:
        return (os.path.join(scene['liparams']['objfilebase'], "{}-{}.obj".format(name.replace(" ", "_"), bpy.context.scene.frame_start)))


def retelaarea(node):
    inlinks = [sock.links[0] for sock in node.inputs if sock.bl_idname in ('So_En_Net_SSFlow', 'So_En_Net_SFlow') and sock.links]
    outlinks = [sock.links[:] for sock in node.outputs if sock.bl_idname in ('So_En_Net_SSFlow', 'So_En_Net_SFlow') and sock.links]
    inosocks = [link.from_socket for link in inlinks if inlinks and link.from_socket.node.get('zone') and link.from_socket.node.zone in [o.name for o in bpy.data.objects]]
    outosocks = [link.to_socket for x in outlinks for link in x if link.to_socket.node.get('zone') and link.to_socket.node.zone in [o.name for o in bpy.data.objects]]

    if outosocks or inosocks:
        elaarea = max([facearea(bpy.data.objects[sock.node.zone], bpy.data.objects[sock.node.zone].data.polygons[int(sock.sn)]) for sock in outosocks + inosocks])
        node["_RNA_UI"] = {"ela": {"max": elaarea, "min": 0.0001}}


def objmode():
    if bpy.context.active_object and bpy.context.active_object.type == 'MESH' and not bpy.context.active_object.hide_viewport:
        bpy.ops.object.mode_set(mode='OBJECT')


def retmesh(name, fr, node, scene):
    if node.animmenu in ("Geometry", "Material"):
        return (os.path.join(scene['liparams']['objfilebase'], '{}-{}.mesh'.format(name.replace(" ", "_"), fr)))
    else:
        return (os.path.join(scene['liparams']['objfilebase'], '{}-{}.mesh'.format(name.replace(" ", "_"), bpy.context.scene.frame_start)))


def nodeinputs(node):
    try:
        ins = [i for i in node.inputs if not i.hide]
        if ins and not all([i.links for i in ins]):
            return 0
        elif ins and any([i.links[0].from_node.use_custom_color for i in ins if i.links]):
            return 0
        else:
            inodes = [i.links[0].from_node for i in ins if i.links[0].from_node.inputs]
            for inode in inodes:
                iins = [i for i in inode.inputs if not i.hide]
                if iins and not all([i.is_linked for i in iins]):
                    return 0
                elif iins and not all([i.links[0].from_node.use_custom_color for i in iins if i.is_linked]):
                    return 0
        return 1
    except Exception as e:
        print(e)


def retmat(fr, node, scene):
    if node.animmenu == "Material":
        return ("{}-{}.rad".format(scene['viparams']['filebase'], fr))
    else:
        return ("{}-{}.rad".format(scene['viparams']['filebase'], scene.frame_start))


def retsky(fr, node, scene):
    if node.animmenu == "Time":
        return ("{}-{}.sky".format(scene['viparams']['filebase'], fr))
    else:
        return ("{}-{}.sky".format(scene['viparams']['filebase'], scene.frame_start))


def nodeexported(self):
    self.exported = 0


def negneg(x):
    x = 0 if float(x) < 0 else x
    return float(x)


def clearanim(scene, obs):
    for o in obs:
        selobj(scene, o)
        o.animation_data_clear()
        o.data.animation_data_clear()

        while o.data.shape_keys:
            bpy.context.object.active_shape_key_index = 0
            bpy.ops.object.shape_key_remove(all=True)


def clearfiles(filebase):
    if os.path.isdir(filebase):
        fileList = os.listdir(filebase)

        for fileName in fileList:
            try:
                os.remove(os.path.join(filebase, fileName))
            except Exception:
                pass


def clearscene(context, op):
    scene = context.scene
    svp = scene.vi_params
    svp['viparams']['vidisp'] = ''

    if context.mode == 'EDIT_MESH':
        bpy.ops.object.mode_set(mode='OBJECT')

    if context.view_layer.layer_collection.children.get('LiVi Results'):
        context.view_layer.layer_collection.children['LiVi Results'].exclude = 0

        for ob in context.view_layer.layer_collection.children['LiVi Results'].collection.objects:
            if ob.vi_params.vi_type_string == 'LiVi Res':
                delobj(context.view_layer, ob)

        context.view_layer.layer_collection.children['LiVi Results'].exclude = 1

    for ob in [ob for ob in scene.objects if ob.type == 'MESH' and ob.visible_get()]:
        if ob.vi_params.vi_type_string not in ('LiVi Calc', 'LiVi Res') and 'cindex' in [a.name for a in ob.data.attributes]:
            if 'export' in op.name or 'simulation' in op.name:
                for r in ob.data.attributes:
                    if not r.is_required and r not in ('UVMap', 'sharp_face'):
                        ob.data.attributes.remove(r)

    for mesh in bpy.data.meshes:
        if mesh.users == 0:
            bpy.data.meshes.remove(mesh)

    for lamp in bpy.data.lights:
        if lamp.users == 0:
            bpy.data.lights.remove(lamp)

    for oldgeo in bpy.data.objects:
        if oldgeo.users == 0:
            bpy.data.objects.remove(oldgeo)

    for sk in bpy.data.shape_keys:
        if sk.users == 0:
            for keys in sk.keys():
                keys.animation_data_clear()

    scene.vi_params['liparams']['livir'] = []


def rtupdate(self, context):
    try:
        rl = self.links[0].from_node['reslists']
        zri = set([(zr[1], zr[1], 'Plot {}'.format(zr[1])) for zr in rl if zr[0] == self.framemenu])
        return zri
    except Exception:
        return []


def iprop(iname, idesc, imin, imax, idef):
    return (IntProperty(name=iname, description=idesc, min=imin, max=imax, default=idef))


def eprop(eitems, ename, edesc, edef):
    return (EnumProperty(items=eitems, name=ename, description=edesc, default=edef))


def bprop(bname, bdesc, bdef):
    return (BoolProperty(name=bname, description=bdesc, default=bdef))


def sprop(sname, sdesc, smaxlen, sdef):
    return (StringProperty(name=sname, description=sdesc, maxlen=smaxlen, default=sdef))


def fprop(fname, fdesc, fmin, fmax, fdef):
    return (FloatProperty(name=fname, description=fdesc, min=fmin, max=fmax, default=fdef))


def fvprop(fvsize, fvname, fvattr, fvdef, fvsub, fvmin, fvmax):
    return (FloatVectorProperty(size=fvsize, name=fvname, attr=fvattr, default=fvdef, subtype=fvsub, min=fvmin, max=fvmax))


def niprop(iname, idesc, imin, imax, idef):
    return (IntProperty(name=iname, description=idesc, min=imin, max=imax, default=idef, update=nodeexported))


def neprop(eitems, ename, edesc, edef):
    return (EnumProperty(items=eitems, name=ename, description=edesc, default=edef, update=nodeexported))


def nbprop(bname, bdesc, bdef):
    return (BoolProperty(name=bname, description=bdesc, default=bdef, update=nodeexported))


def nsprop(sname, sdesc, smaxlen, sdef):
    return (StringProperty(name=sname, description=sdesc, maxlen=smaxlen, default=sdef, update=nodeexported))


def nfprop(fname, fdesc, fmin, fmax, fdef):
    return (FloatProperty(name=fname, description=fdesc, min=fmin, max=fmax, default=fdef, update=nodeexported))


def nfvprop(fvname, fvattr, fvdef, fvsub):
    return (FloatVectorProperty(name=fvname, attr=fvattr, default=fvdef, subtype=fvsub, update=nodeexported))


def vertarea(mesh, vert):
    area = 0
    faces = [face for face in vert.link_faces]

    if hasattr(mesh.verts, "ensure_lookup_table"):
        mesh.verts.ensure_lookup_table()

    if len(faces) > 1:
        for f, face in enumerate(faces):
            ovs, oes = [], []
            fvs = [le.verts[(0, 1)[le.verts[0] == vert]] for le in vert.link_edges]
            ofaces = [oface for oface in faces if len([v for v in oface.verts if v in face.verts]) == 2]

            for oface in ofaces:
                oes.append([e for e in face.edges if e in oface.edges])
                ovs.append([i for i in face.verts if i in oface.verts])

            if len(ovs) == 1:
                sedgevs = (vert.index, [v.index for v in fvs if v not in ovs][0])
                sedgemp = mathutils.Vector([((mesh.verts[sedgevs[0]].co)[i] + (mesh.verts[sedgevs[1]].co)[i]) / 2 for i in range(3)])

                if not mathutils.geometry.intersect_line_line(face.calc_center_median(), ofaces[0].calc_center_median(), ovs[0][0].co, ovs[0][1].co):
                    return 0
                eps = [mathutils.geometry.intersect_line_line(face.calc_center_median(), ofaces[0].calc_center_median(), ovs[0][0].co, ovs[0][1].co)[1]] + [sedgemp]

            elif len(ovs) == 2:
                if None in [mathutils.geometry.intersect_line_line(face.calc_center_median(), ofaces[i].calc_center_median(), ovs[i][0].co, ovs[i][1].co) for i in range(2)]:
                    return 0
                eps = [mathutils.geometry.intersect_line_line(face.calc_center_median(), ofaces[i].calc_center_median(), ovs[i][0].co, ovs[i][1].co)[1] for i in range(2)]
            else:
                return 0

            area += mathutils.geometry.area_tri(vert.co, *eps) + mathutils.geometry.area_tri(face.calc_center_median(), *eps)

    elif len(faces) == 1:
        eps = [(ev.verts[0].co + ev.verts[1].co) / 2 for ev in vert.link_edges]
        eangle = (vert.link_edges[0].verts[0].co - vert.link_edges[0].verts[1].co).angle(vert.link_edges[1].verts[0].co - vert.link_edges[1].verts[1].co)

        if eangle:
            area = mathutils.geometry.area_tri(vert.co, *eps) + mathutils.geometry.area_tri(faces[0].calc_center_median(), *eps) * 2 * pi / eangle
        else:
            return 0

    return area


def facearea(obj, face):
    omw = obj.matrix_world
    vs = [omw @ mathutils.Vector(face.center)] + [omw @ obj.data.vertices[v].co for v in face.vertices] + [omw @ obj.data.vertices[face.vertices[0]].co]
    return (vsarea(obj, vs))


def vsarea(obj, vs):
    if len(vs) == 5:
        cross = mathutils.Vector.cross(vs[3] - vs[1], vs[3] - vs[2])
        return (0.5 * (cross[0]**2 + cross[1]**2 + cross[2]**2)**0.5)
    else:
        i, area = 0, 0
        while i < len(vs) - 2:
            cross = mathutils.Vector.cross(vs[0] - vs[1 + i], vs[0] - vs[2 + i])
            area += 0.5 * (cross[0]**2 + cross[1]**2 + cross[2]**2)**0.5
            i += 1
        return (area)


def wind_rose(wro, maxws, wrsvg, wrtype, colors):
    zp = 0
    bm = bmesh.new()
    wro.location = (0, 0, 0)
    svg = minidom.parse(wrsvg)
    pos_strings = [path.getAttribute('d') for path in svg.getElementsByTagName('path')]
    style_strings = [path.getAttribute('style').split(';') for path in svg.getElementsByTagName('path')]
    dimen = [eval(path.getAttribute('height').strip('pt')) for path in svg.getElementsByTagName('svg')][0]
    scale = 0.04 * dimen
    svg.unlink()
    sposnew = [[(eval(ss.split()[ss.index('M') + 1]) - dimen / 2) * 0.1, (eval(ss.split()[ss.index('M') + 2]) - dimen / 2) * -0.1, 0.05] for ss in pos_strings]
    lposnew = [[[(eval(ss.split()[li + 1]) - dimen / 2) * 0.1, (eval(ss.split()[li + 2]) - dimen / 2) * -0.1, 0.05] for li in [si for si, s in enumerate(ss.split()) if s == 'L']] for ss in pos_strings]

    for stsi, sts in enumerate(style_strings):
        if ('fill: #' in sts[0] or 'fill:#' in sts[0]) and sts[0][-6:] != 'ffffff':
            hexcol, col = sts[0][-7:], sts[0][-6:]
            fillrgb = colors.hex2color(hexcol)

            if 'wr-{}'.format(col) not in [mat.name for mat in bpy.data.materials]:
                bpy.data.materials.new('wr-{}'.format(col))

            bpy.data.materials['wr-{}'.format(col)].diffuse_color = [c for c in fillrgb] + [1]

            if 'wr-{}'.format(col) not in [mat.name for mat in wro.data.materials]:
                bpy.ops.object.material_slot_add()
                wro.material_slots[-1].material = bpy.data.materials['wr-{}'.format(col)]

            if wrtype in ('2', '3', '4'):
                vs = [bm.verts.new(pos) for pos in [sposnew[stsi]] + lposnew[stsi][:-1]]
            else:
                vs = [bm.verts.new(pos) for pos in [sposnew[stsi]] + lposnew[stsi]]

            vs.reverse()

            if len(vs) > 2:
                nf = bm.faces.new(vs[::])
                nf.material_index = wro.data.materials[:].index(wro.data.materials['wr-{}'.format(col)])

                if wrtype in ('2', '3', '4'):
                    zp += 0.0005 * scale

                    for vert in nf.verts:
                        vert.co[2] = zp

    if 'wr-000000' not in [mat.name for mat in bpy.data.materials]:
        bpy.data.materials.new('wr-000000')
        bpy.data.materials['wr-000000'].diffuse_color = (0, 0, 0, 1)

    if 'wr-000000' not in [mat.name for mat in wro.data.materials]:
        bpy.ops.object.material_slot_add()
        wro.material_slots[-1].material = bpy.data.materials['wr-000000']

    bmesh.ops.remove_doubles(bm, verts=vs, dist=scale * 0.01)

    if wrtype in ('0', '1', '3', '4'):
        thick = scale * 0.005 if wrtype == '4' else scale * 0.0005
        faces = bmesh.ops.inset_individual(bm, faces=bm.faces, thickness=thick, use_even_offset=True)['faces']

        if wrtype == '4':
            [bm.faces.remove(f) for f in bm.faces if f not in faces]
        else:
            bi = wro.data.materials[:].index(wro.data.materials['wr-000000'])

            for face in faces:
                face.material_index = bi

    bm.to_mesh(wro.data)
    bm.free()

    bpy.ops.mesh.primitive_circle_add(vertices=132, radius=scale * 1.2, fill_type='NGON', align='WORLD', enter_editmode=False, location=(0, 0, -0.01))
    wrbo = bpy.context.active_object

    if 'wr-base'not in [mat.name for mat in bpy.data.materials]:
        bpy.data.materials.new('wr-base')
        bpy.data.materials['wr-base'].diffuse_color = (1, 1, 1, 1)

    bpy.ops.object.material_slot_add()
    wrbo.material_slots[-1].material = bpy.data.materials['wr-base']
    return ((wrbo, wro), scale)


def compass(loc, scale, platmat, basemat, greymat):
    bpy.ops.wm.append(filepath="sp.blend", directory=os.path.join(os.path.dirname(os.path.realpath(__file__)),
                      'Images/sp.blend', 'Object'), filename="SPathMesh", autoselect=True)

    coo = bpy.data.objects['SPathMesh']
    bpy.context.view_layer.objects.active = coo
    coo.material_slots[0].material = basemat
    coo.material_slots[1].material = platmat
    coo.material_slots[2].material = greymat
    return coo


def spathrange(mats):
    sprme = bpy.data.meshes.new("SPRange")
    spro = bpy.data.objects.new('SPRrange', sprme)
    bpy.context.scene.collection.objects.link(spro)
    bpy.context.view_layer.objects.active = spro
    spro.location = (0, 0, 0)
    bm = bmesh.new()
    params = ((177, 0.05, 0), (80, 0.1, 1), (355, 0.15, 2)) if bpy.context.scene.latitude >= 0 else ((355, 0.05, 0), (80, 0.1, 1), (177, 0.15, 2))

    for param in params:
        bpy.ops.object.material_slot_add()
        spro.material_slots[-1].material = mats[param[2]]
        morn = solarRiseSet(param[0], 0, bpy.context.scene.latitude, bpy.context.scene.longitude, 'morn')
        eve = solarRiseSet(param[0], 0, bpy.context.scene.latitude, bpy.context.scene.longitude, 'eve')

        if morn or param[2] == 0:
            if morn:
                mornevediff = eve - morn if bpy.context.scene.latitude >= 0 else 360 - eve + morn
            else:
                mornevediff = 360

            startset = morn if bpy.context.scene.latitude >= 0 else eve
            angrange = [startset + a * 0.0125 * mornevediff for a in range(0, 81)]
            bm.verts.new().co = (95 * sin(angrange[0] * pi / 180), 95 * cos(angrange[0] * pi / 180), param[1])

            for a in angrange[1:-1]:
                bm.verts.new().co = (92 * sin(a * pi / 180), 92 * cos(a * pi / 180), param[1])

            bm.verts.new().co = (95 * sin(angrange[len(angrange) - 1] * pi / 180), 95 * cos(angrange[len(angrange) - 1] * pi / 180), param[1])
            angrange.reverse()

            for b in angrange[1:-1]:
                bm.verts.new().co = (98 * sin(b * pi / 180), 98 * cos(b * pi / 180), param[1])

            bm.faces.new(bm.verts[-160:])
            bm.faces.ensure_lookup_table()
            bm.faces[-1].material_index = param[2]

    bm.to_mesh(sprme)
    bm.free()
    return spro


def windnum(maxws, loc, scale, wr):
    txts = []
    matrot = Matrix.Rotation(-pi * 0.05, 4, 'Z')
    direc = Vector((0, 1, 0))

    for i in range(2, 6):
        bpy.ops.object.text_add(align='WORLD', enter_editmode=False, location=((i**2) / 25) * scale * (matrot @ direc))
        txt = bpy.context.active_object
        txt.data.body, txt.scale, txt.location[2] = '{:.1f}'.format((i**2) * maxws / 25), (scale * 0.05, scale * 0.05, scale * 0.05), scale * 0.01
        bpy.ops.object.convert(target='MESH')
        bpy.ops.object.material_slot_add()
        txt.material_slots[-1].material = bpy.data.materials['wr-000000']
        txts.append(txt)

    joinobj(bpy.context.view_layer, txts + [wr]).name = 'Wind Rose'
    bpy.context.active_object.visible_shadow = False
    bpy.context.active_object.display.show_shadows = False
    bpy.context.active_object['rpe'] = 'Wind_Plane'


def wind_compass(loc, scale, wro, mat):
    txts = []
    come = bpy.data.meshes.new("Compass")
    coo = bpy.data.objects.new('Compass', come)
    coo.location = loc
    bpy.context.scene.collection.objects.link(coo)
    bpy.context.view_layer.objects.active = coo
    bpy.ops.object.material_slot_add()
    coo.material_slots[-1].material = mat
    bm = bmesh.new()
    matrot = Matrix.Rotation(pi * 0.25, 4, 'Z')

    for i in range(1, 11):
        bmesh.ops.create_circle(bm, cap_ends=False, radius=scale * ((i**2) / 10) * 0.1, segments=132, matrix=Matrix.Rotation(pi / 64, 4, 'Z') @ Matrix.Translation((0, 0, 0)))
    for edge in bm.edges:
        edge.select_set(False) if edge.index % 3 or edge.index > 1187 else edge.select_set(True)

    bmesh.ops.delete(bm, geom=[edge for edge in bm.edges if edge.select], context='FACES_ONLY')
    newgeo = bmesh.ops.extrude_edge_only(bm, edges=bm.edges, use_select_history=False)

    for v, vert in enumerate(newgeo['geom'][:1320]):
        vert.co = vert.co - (vert.co - coo.location).normalized() * scale * (0.0025, 0.005)[v > 1187]
        vert.co[2] = 0

    bmesh.ops.create_circle(bm, cap_ends=True, radius=scale * 0.0025, segments=8, matrix=Matrix.Rotation(-pi / 8, 4, 'Z') @ Matrix.Translation((0, 0, 0)))
    matrot = Matrix.Rotation(pi * 0.25, 4, 'Z')
    degmatrot = Matrix.Rotation(pi * 0.125, 4, 'Z')
    tmatrot = Matrix.Rotation(0, 4, 'Z')
    direc = Vector((0, 1, 0))

    for i, edge in enumerate(bm.edges[-8:]):
        verts = bmesh.ops.extrude_edge_only(bm, edges=[edge], use_select_history=False)['geom'][:2]

        for vert in verts:
            vert.co += 1.0 * scale * (tmatrot @ direc)
            vert.co[2] = 0

        bpy.ops.object.text_add(align='WORLD', enter_editmode=False, location=Vector(loc) + scale * 1.13 * (tmatrot @ direc), rotation=tmatrot.to_euler())
        txt = bpy.context.active_object
        txt.scale, txt.data.body, txt.data.align_x, txt.data.align_y, txt.location[2] = (scale * 0.075, scale * 0.075, scale * 0.075), ('N', 'NW', 'W', 'SW', 'S', 'SE', 'E', 'NE')[i], 'CENTER', 'CENTER', txt.location[2]
        bpy.ops.object.convert(target='MESH')
        bpy.ops.object.material_slot_add()
        txt.material_slots[-1].material = mat
        txts.append(txt)
        tmatrot = tmatrot @ matrot

    tmatrot = Matrix.Rotation(0, 4, 'Z')

    for d in range(16):
        bpy.ops.object.text_add(align='WORLD', enter_editmode=False, location=Vector(loc) + scale * 1.04 * (tmatrot @ direc), rotation=tmatrot.to_euler())
        txt = bpy.context.active_object
        txt.scale, txt.data.body, txt.data.align_x, txt.data.align_y, txt.location[2] = (scale * 0.05, scale * 0.05, scale * 0.05), (u'0\u00B0', u'337.5\u00B0', u'315\u00B0', u'292.5\u00B0', u'270\u00B0', u'247.5\u00B0', u'225\u00B0', u'202.5\u00B0', u'180\u00B0', u'157.5\u00B0', u'135\u00B0', u'112.5\u00B0', u'90\u00B0', u'67.5\u00B0', u'45\u00B0', u'22.5\u00B0')[d], 'CENTER', 'CENTER', txt.location[2]
        bpy.ops.object.convert(target='MESH')
        bpy.ops.object.material_slot_add()
        txt.material_slots[-1].material = mat
        txts.append(txt)
        tmatrot = tmatrot @ degmatrot

    bm.to_mesh(come)
    bm.free()
    return joinobj(bpy.context.view_layer, txts + [coo] + [wro])


def rgb2h(rgb):
    return colorsys.rgb_to_hsv(rgb[0] / 255.0, rgb[1] / 255.0, rgb[2] / 255.0)[0]


def mp2im(fig, imname):
    fig.canvas.draw()
    ipwidth, ipheight = fig.canvas.width(), fig.canvas.height()

    if imname not in [im.name for im in bpy.data.images]:
        bpy.ops.image.new(name=imname, width=ipwidth, height=ipheight, color=(1, 1, 1, 1), alpha=True, generated_type='BLANK', float=False, use_stereo_3d=False)
        im = bpy.data.images[imname]

    else:
        im = bpy.data.images[imname]
        im.buffers_free()

        if im.size[:] != (ipwidth, ipheight):
            im.scale(ipwidth, ipheight)

    rgba = +frombuffer(fig.canvas.buffer_rgba(), dtype=uint8)
    rgba.shape = (ipheight, ipwidth, 4)
    rgba = rgba[::-1].ravel()
    rgba = multiply(rgba, 0.00392, dtype=float32)
    im.pixels.foreach_set(rgba)


def livisimacc(simnode):
    context = simnode.inputs['Context in'].links[0].from_node['Options']['Context']
    return (simnode.csimacc if context in ('Compliance', 'CBDM') else simnode.simacc)


def radial2xy(c, theta, phi, w, h):
    return c[0] + theta * sin(math.pi * phi / 180) * w, c[1] + theta * math.cos(math.pi * phi / 180) * h


def drawfont(text, fi, lencrit, height, x1, y1):
    blf.position(fi, x1, height - y1 - lencrit * 26, 0)
    blf.draw(fi, text)


def xy2radial(c, pos, w, h):
    dx, dy = pos[0] - c[0], pos[1] - c[1]
    hypo = (((dx / w)**2 + (dy / h)**2)**0.5)
    at = math.atan((dy / h) / (dx / w))

    if dx == 0:
        azi = 0 if dy >= 0 else math.pi
    elif dx > 0:
        at = math.atan((dy / h) / (dx / w))
        azi = math.pi * 0.5 - at
    elif dx < 0:
        at = math.atan((dy / h) / (dx / w))
        azi = math.pi * 1.5 - at
    return hypo, azi


def bres(scene, o):
    bm = bmesh.new()
    bm.from_mesh(o.data)

    if scene['liparams']['cp'] == '1':
        rtlayer = bm.verts.layers.int['cindex']
        reslayer = bm.verts.layers.float['res{}'.format(scene.frame_current)]
        res = [v[reslayer] for v in bm.verts if v[rtlayer] > 0]
    elif scene['liparams']['cp'] == '0':
        rtlayer = bm.faces.layers.int['cindex']
        reslayer = bm.faces.layers.float['res{}'.format(scene.frame_current)]
        res = [f[reslayer] for f in bm.faces if f[rtlayer] > 0]
    
    bm.free()
    return res


def framerange(scene, anim):
    if anim == 'Static':
        return (range(scene.frame_current, scene.frame_current + 1))
    else:
        return (range(scene.frame_start, scene['liparams']['fe'] + 1))


def frameindex(scene, anim):
    if anim == 'Static':
        return (range(0, 1))
    else:
        return (range(0, scene.frame_end - scene.frame_start + 1))


def ret_camera_menu(self, context):
    cameras = [o for o in bpy.data.objects if o.type == 'CAMERA']
    if cameras:
        return [(o.name, o.name, 'Name of the section plane') for o in cameras]
    else:
        return [('None', 'None', 'None')]


def retrobjs(obs):
    robjs = []
    if bpy.data.collections.get('LiVi Res'):
        for o in bpy.data.collections['LiVi Res'].objects:
            robjs.append(o)
    if bpy.data.collections.get('EnVi Geometry'):
        for o in bpy.data.collections['EnVi Geometry'].objects:
            robjs.append(o)
    if bpy.data.collections.get('FloVi Mesh'):
        for o in bpy.data.collections['FloVi Mesh'].objects:
            robjs.append(o)
    if bpy.data.collections.get('SunPath'):
        for o in bpy.data.collections['SunPath'].objects:
            robjs.append(o)
    return robjs


def retobjs(otypes):
    scene = bpy.context.scene
    validobs = [o for o in scene.objects if o.visible_get() and '/' not in o.name and o not in retrobjs(scene.objects)]

    for o in scene.objects:
        if '/' in o.name:
            logentry('Object {} has a "/" in the name and will not be exported'.format(o.name))
        if o not in validobs and o.vi_params.vi_type_string != "LiVi Res":
            o.vi_params.vi_type_string = ''

    if otypes == 'livig':
        return [o for o in validobs if o.type == 'MESH' and o.data.polygons and any(o.data.materials) and not (o.parent and os.path.isfile(o.vi_params.ies_name))
                and o.vi_params.vi_type not in ('4', '5')
                and o.vi_params.vi_type_string != 'LiVi Res'
                and o.get('VIType') not in ('SPathMesh', 'SunMesh', 'Wind_Plane', 'SkyMesh')]
    elif otypes == 'livigeno':
        return [o for o in validobs if o.type == 'MESH' and o.data.polygons and o.data.materials and not any([m.vi_params.livi_sense for m in o.data.materials])]
    elif otypes == 'livigengeosel':
        return [o for o in validobs if o.type == 'MESH' and o.data.polygons and o.select is True and o.data.materials and not any([m.vi_params.livi_sense for m in o.data.materials])]
    elif otypes == 'livil':
        return [o for o in validobs if o.type == 'LIGHT' or o.vi_params.vi_type == '4']
    elif otypes == 'livic':
        return [o for o in validobs if o.type == 'MESH' and o.data.polygons and li_calcob(o, 'livi') and o.vi_params.vi_type_string != 'LiVi Res']
    elif otypes == 'livir':
        return [o for o in validobs if o.type == 'MESH' and o.data.polygons and True in [m.vi_params.livi_sense for m in o.data.materials] and o.vi_params.vi_type_string != 'LiVi Calc']
    elif otypes == 'envig':
        return [o for o in scene.objects if o.type == 'MESH' and o.hide is False]
    elif otypes == 'ssc':
        return [o for o in validobs if o.type == 'MESH' and o.data.polygons and li_calcob(o, 'livi') and o.vi_params.vi_type_string != 'LiVi Res']
    elif otypes == 'selected':
        return [o for o in [bpy.context.active_object] if o.type == 'MESH']


def radmesh(scene, obs, export_op):
    for o in obs:
        for mat in o.data.materials:
            if mat['radentry'] and mat['radentry'].split(' ')[1] in ('light', 'mirror', 'antimatter') or mat.pport:
                export_op.report({'INFO'}, o.name + " has an antimatter, photon port, emission or mirror material. Basic export routine used with no modifiers.")
                o.vi_params['merr'] = 1
        selobj(scene, o)

        if not o.get('merr'):
            o.vi_params['merr'] = 0


def viewdesc(context):
    region = context.region
    width, height = region.width, region.height
    mid_x, mid_y = width / 2, height / 2
    return (mid_x, mid_y, width, height)


def skfpos(o, frame, vis):
    vcos = [o.matrix_world * o.data.shape_keys.key_blocks[str(frame)].data[v].co for v in vis]
    maxx = max([vco[0] for vco in vcos])
    minx = min([vco[0] for vco in vcos])
    maxy = max([vco[1] for vco in vcos])
    miny = min([vco[1] for vco in vcos])
    maxz = max([vco[2] for vco in vcos])
    minz = min([vco[2] for vco in vcos])
    return mathutils.Vector(((maxx + minx) * 0.5, (maxy + miny) * 0.5, (maxz + minz) * 0.5))


def selmesh(sel):
    if bpy.context.active_object.mode != 'EDIT':
        bpy.ops.object.mode_set(mode='EDIT')
    if sel == 'selenm':
        bpy.ops.mesh.select_mode(type="EDGE")
        bpy.ops.mesh.select_all(action='DESELECT')
        bpy.ops.mesh.select_non_manifold()
    elif sel == 'desel':
        bpy.ops.mesh.select_all(action='DESELECT')
    elif sel in ('delf', 'rd'):
        if sel == 'delf':
            bpy.ops.mesh.delete(type='FACE')
        bpy.ops.mesh.select_all(action='SELECT')
        bpy.ops.mesh.remove_doubles()
        bpy.ops.mesh.select_all(action='DESELECT')
    elif sel == 'dele':
        bpy.ops.mesh.delete(type='EDGE')
    elif sel == 'delv':
        bpy.ops.mesh.delete(type='VERT')
    elif sel == 'mc':
        bpy.ops.mesh.select_all(action='SELECT')
        bpy.ops.mesh.vert_connect_concave()
        bpy.ops.mesh.select_all(action='DESELECT')
    elif sel == 'mp':
        bpy.ops.mesh.select_all(action='SELECT')
        bpy.ops.mesh.vert_connect_nonplanar(angle_limit=0.001)
        bpy.ops.mesh.select_all(action='DESELECT')
    elif sel in ('SELECT', 'INVERT', 'PASS'):
        if sel in ('SELECT', 'INVERT'):
            bpy.ops.mesh.select_all(action=sel)
        bpy.ops.object.vertex_group_assign()
    bpy.ops.object.mode_set(mode='OBJECT')


def retdp(mres, dp):
    try:
        dp = 0 if ceil(log10(100 / mres)) < 0 or not dp else ceil(log10(100 / mres))
    except Exception:
        dp = 0
    return dp


def draw_index_distance(posis, res, fontsize, fontcol, shadcol, distances):
    if distances.size:
        try:
            dp = 0 if max(res) > 100 else 2 - int(log10(max(res)))
            dpi = bpy.context.preferences.system.dpi
            nres = char.mod('%.{}f'.format(dp), res)
            fsdist = (fontsize / distances).astype(int16)
            xposis = posis[0::2]
            yposis = posis[1::2]
            alldata = zip(nres, fsdist, xposis, yposis)
            ysize = int(0.5 * blf.dimensions(0, nres[0])[1])

            for ad in alldata:
                blf.size(0, ad[1] * dpi / 72)
                blf.position(0, ad[2] - int(0.5 * blf.dimensions(0, ad[0])[0]), ad[3] - ysize, 0)
                blf.draw(0, ad[0])

        except Exception as e:
            print('Drawing index error: ', e)


def draw_index(posis, res, dists, fontsize, fontcol, shadcol, **kwargs):
    if not kwargs.get('hour'):
        nres = ['{}'.format(format(r, '.{}f'.format(retdp(max(res), 0)))) for ri, r in enumerate(res)]
    else:
        nres = res

    for ri, nr in enumerate(nres):
        blf.size(0, int(0.25 * fontsize + 0.25 * fontsize * (max(dists) - dists[ri]) / (max(dists) - min(dists))) * 150 / 72)
        blf.position(0, posis[ri][0] - int(0.5 * blf.dimensions(0, nr)[0]), posis[ri][1] - int(0.5 * blf.dimensions(0, nr)[1]), 0.0)
        blf.draw(0, nr)

    blf.disable(0, 4)


def draw_time(pos, time, fontsize, fontcol, shadcol):
    blf.size(0, int(fontsize))
    blf.position(0, pos[0], pos[1] - blf.dimensions(0, time)[1] * 0.5, 0)
    blf.draw(0, time)
    blf.disable(0, 4)


def blf_props(scene, width, height):
    svp = scene.vi_params
    blf.enable(0, 2)
    blf.clipping(0, 0, 0, width, height)

    if svp.vi_display_rp_sh:
        blf.enable(0, 4)
        blf.shadow(0, 3, *svp.vi_display_rp_fsh)

    blf.color(0, *svp.vi_display_rp_fc)
    blf.size(0, svp.vi_display_rp_fs * int(width / 20) / 72)


def blf_unprops():
    blf.disable(0, 2)
    blf.disable(0, 4)
    blf.color(0, 0, 0, 1)


def edgelen(ob, edge):
    omw = ob.matrix_world
    vdiff = omw * (ob.data.vertices[edge.vertices[0]].co - ob.data.vertices[edge.vertices[1]].co)
    mathutils.Vector(vdiff).length


def sunpath1(self, context):
    sunpath(bpy.context)


def sunpath2():
    sunpath(bpy.context)


def sunpath(context):
    scene = context.scene
    svp = scene.vi_params
    suns = [ob for ob in scene.objects if ob.parent and ob.type == 'LIGHT' and ob.data.type == 'SUN' and ob.parent.get('VIType') == "SPathMesh"]

    if svp.get('spparams') and svp['spparams'].get('suns') and svp['spparams']['suns'] == '0':
        if suns:
            alt, azi, beta, phi = solarPosition(svp.sp_sd, svp.sp_sh, svp.latitude, svp.longitude)
            suns[0].location.z = 100 * sin(beta)
            suns[0].location.x = -(100**2 - (suns[0].location.z)**2)**0.5 * sin(phi)
            suns[0].location.y = -(100**2 - (suns[0].location.z)**2)**0.5 * cos(phi)
            suns[0].data.energy = svp.sp_sun_strength
            suns[0].data.angle = svp.sp_sun_angle
            suns[0].rotation_euler = pi * 0.5 - beta, 0, -phi
            sky_vec = suns[0].location.normalized()
            sky_vec.rotate(suns[0].parent.matrix_world.to_euler())
            scene.display.light_direction = (sky_vec[0], sky_vec[2], -sky_vec[1])

            if scene.render.engine == 'CYCLES':
                if scene.world.node_tree:
                    for stnode in [no for no in scene.world.node_tree.nodes if no.bl_label == 'Sky Texture']:
                        stnode.sun_direction = -sin(phi), -cos(phi), sin(beta)
                        for bnode in [no for no in scene.world.node_tree.nodes if no.bl_label == 'Background']:
                            bnode.inputs[1].default_value = 1.5 + sin(beta) * 0.5

                if suns[0].data.node_tree:
                    for blnode in [node for node in suns[0].data.node_tree.nodes if node.bl_label == 'Blackbody']:
                        blnode.inputs[0].default_value = 3000 + 2500 * sin(beta)**0.5 if beta > 0 else 2500
                    for emnode in [node for node in suns[0].data.node_tree.nodes if node.bl_label == 'Emission']:
                        emnode.inputs[1].default_value = 10 * sin(beta)**0.5 if beta > 0 else 0

            suns[0]['solhour'], suns[0]['solday'] = svp.sp_sh, svp.sp_sd
            suns[0].hide_viewport = True if alt <= 0 else False

            if suns[0].children:
                suns[0].children[0].hide_viewport = True if alt <= 0 else False

            return

    elif svp.get('spparams') and svp['spparams'].get('suns') and svp['spparams']['suns'] == '1':
        all_alts = [solarPosition(d, svp.sp_sh, svp.latitude, svp.longitude)[0] for d in (20, 50, 80, 110, 140, 171, 201, 231, 261, 292, 323, 354)]
        valid_suns = len([aa for aa in all_alts if aa > 0])

        for d, day in enumerate((20, 50, 80, 110, 140, 171, 201, 231, 261, 292, 323, 354)):
            alt, azi, beta, phi = solarPosition(day, svp.sp_sh, svp.latitude, svp.longitude)
            suns[d].location.z = 100 * sin(beta)
            suns[d].location.x = -(100**2 - (suns[d].location.z)**2)**0.5 * sin(phi)
            suns[d].location.y = -(100**2 - (suns[d].location.z)**2)**0.5 * cos(phi)
            suns[d].rotation_euler = pi * 0.5 - beta, 0, -phi
            suns[d].hide_viewport = True if alt <= 0 else False

            if suns[d].children:
                suns[d].children[0].hide_viewport = True if alt <= 0 else False
            if alt > 0:
                suns[d].data.energy = 1.5 * svp.sp_sun_strength / (sin(beta) * valid_suns)

                if scene.render.engine == 'CYCLES':
                    if suns[d].data.node_tree:
                        for emnode in [node for node in suns[d].data.node_tree.nodes if node.bl_label == 'Emission']:
                            emnode.inputs[1].default_value = 1.5 * svp.sp_sun_strength / (sin(beta) * valid_suns)

                suns[d].data.angle = svp.sp_sun_angle

    elif svp.get('spparams') and svp['spparams'].get('suns') and svp['spparams']['suns'] == '2':
        all_alts = [solarPosition(svp.sp_sd, h, svp.latitude, svp.longitude)[0] for h in range(24)]
        valid_suns = len([aa for aa in all_alts if aa > 0])

        for h in range(24):
            alt, azi, beta, phi = solarPosition(svp.sp_sd, h, svp.latitude, svp.longitude)
            if alt < 0:
                suns[h].hide_viewport = True
                if suns[h].children:
                    suns[h].children[0].hide_viewport = True
            else:
                sun_strength = 1.5 * svp.sp_sun_strength / (max(sin(beta), 0.2) * valid_suns)
                suns[h].hide_viewport = False

                if suns[h].children:
                    suns[h].children[0].hide_viewport = False

                suns[h].location.z = 100 * sin(beta)
                suns[h].location.x = -(100**2 - (suns[h].location.z)**2)**0.5 * sin(phi)
                suns[h].location.y = -(100**2 - (suns[h].location.z)**2)**0.5 * cos(phi)
                suns[h].rotation_euler = pi * 0.5 - beta, 0, -phi
                suns[h].data.energy = sun_strength

                if scene.render.engine == 'CYCLES':
                    if suns[h].data.node_tree:
                        for emnode in [node for node in suns[h].data.node_tree.nodes if node.bl_label == 'Emission']:
                            emnode.inputs[1].default_value = sun_strength

                suns[h].data.angle = math.pi * svp.sp_sun_angle / 180


def epentry(header, params, paramvs):
    return '{}\n'.format(header + (',', '')[header == '']) + '\n'.join([('    ', '')[header == ''] + '{:{width}}! - {}'.format(str(pv[0]) + (',', ';')[pv[1] == params[-1]],
                                                                                                                               pv[1], width=80 + (0, 4)[header == '']) for pv in zip(paramvs, params)]) + ('\n\n', '')[header == '']


def epwlatilongi(scene, node):
    with open(node.weather, "r") as epwfile:
        fl = epwfile.readline()
        latitude, longitude = float(fl.split(",")[6]), float(fl.split(",")[7])
    return latitude, longitude


# Compute solar position (altitude and azimuth in degrees) based on day of year (doy; integer), local solar time (lst; decimal hours), latitude (lat; decimal degrees), and longitude (lon; decimal degrees).
def solarPosition(doy, lst, lat, lon):
    # Set the local standard time meridian (lsm) (integer degrees of arc)
    lsm = round(lon / 15, 0) * 15
    # Approximation for equation of time (et) (minutes) comes from the Wikipedia article on Equation of Time
    b = 2 * pi * (doy - 81) / 364
    et = 9.87 * sin(2 * b) - 7.53 * cos(b) - 1.5 * sin(b)
    # The following formulas adapted from the 2005 ASHRAE Fundamentals, pp. 31.13-31.16
    # Conversion multipliers
    degToRad = 2 * pi / 360
    radToDeg = 1 / degToRad
    # Apparent solar time (ast)
    if lon > lsm:
        ast = lst + et / 60 - (lsm - lon) / 15
    else:
        ast = lst + et / 60 + (lsm - lon) / 15
    # Solar declination (delta) (radians)
    delta = degToRad * 23.45 * sin(2 * pi * (284 + doy) / 365)
    # Hour angle (h) (radians)
    h = degToRad * 15 * (ast - 12)
    # Local latitude (ll) (radians)
    ll = degToRad * lat
    # Solar altitude (beta) (radians)
    beta = asin(cos(ll) * cos(delta) * cos(h) + sin(ll) * sin(delta))
    # Solar azimuth phi (radians)
    phi = acos((sin(beta) * sin(ll) - sin(delta)) / (cos(beta) * cos(ll)))
    # Convert altitude and azimuth from radians to degrees, since the Spatial Analyst's Hillshade function inputs solar angles in degrees
    altitude = radToDeg * beta
    phi = 2 * pi - phi if ast <= 12 or ast >= 24 else phi
    azimuth = radToDeg * phi
    return ([altitude, azimuth, beta, phi])


def solarRiseSet(doy, beta, lat, lon, riseset):
    degToRad = 2 * pi / 360
    radToDeg = 1 / degToRad
    delta = degToRad * 23.45 * sin(2 * pi * (284 + doy) / 365)
    ll = degToRad * lat

    try:
        phi = acos((sin(beta) * sin(ll) - sin(delta)) / (cos(beta) * cos(ll)))
        phi = pi - phi if riseset == 'morn' else pi + phi
    except Exception:
        phi = 0

    return (phi * radToDeg)


def skframe(pp, scene, oblist, params):
    svp = scene.vi_params
    for frame in range(svp[params]['fs'], svp[params]['fe'] + 1):
        scene.frame_set(frame)
        for o in [o for o in oblist if o.data.shape_keys]:
            for shape in o.data.shape_keys.key_blocks:
                if shape.name.isdigit():
                    shape.value = shape.name == str(frame)
                    shape.keyframe_insert("value")


def gentarget(tarnode, result):
    if tarnode.stat == '0':
        res = sum(result) / len(result)
    elif tarnode.stat == '1':
        res = max(result)
    elif tarnode.stat == '2':
        res = min(result)
    elif tarnode.stat == '3':
        res = sum(result)

    if tarnode.value > res and tarnode.ab == '0':
        return (1)
    elif tarnode.value < res and tarnode.ab == '1':
        return (1)
    else:
        return (0)


def selobj(vl, geo):
    vl = bpy.context.view_layer
    if vl.objects.active and vl.objects.active.hide_viewport == 'False':
        bpy.ops.object.mode_set(mode='OBJECT')
    for ob in vl.objects:
        bpy.context.view_layer.objects.active
        ob.select_set(1) if ob == geo else ob.select_set(0)
    try:
        vl.objects.active = geo
    except Exception:
        pass


def actselobj(vl, act, geos):
    if vl.objects.active and vl.objects.active.hide_viewport == 'False':
        bpy.ops.object.mode_set(mode='OBJECT')
    for ob in vl.objects:
        bpy.context.view_layer.objects.active
        ob.select_set(1) if ob.name in geos else ob.select_set(0)
    act.select_set(1)
    vl.objects.active = act


def selobs(vl, geos):
    if vl.objects.active and vl.objects.active.hide_viewport == 'False':
        bpy.ops.object.mode_set(mode='OBJECT')
    for ob in vl.objects:
        ob.select_set(1) if ob.name in geos else ob.select_set(0)
    vl.objects.active = vl.objects[geos[0]]


def delobj(vl, delgeo):
    selobj(vl, delgeo)
    bpy.ops.object.delete(use_global=True)


def joinobj(vl, obs):
    bpy.ops.object.select_all(action='DESELECT')

    for o in obs:
        o.select_set(state=True)

    vl.objects.active = obs[-1]
    bpy.ops.object.join()
    return bpy.context.active_object


def nodeid(node):
    for ng in bpy.data.node_groups:
        if node in ng.nodes[:]:
            return node.name + '@' + ng.name


def nodecolour(node, prob):
    (node.use_custom_color, node.color) = (1, (0.8, 0.3, 0.3)) if prob else (0, (0.8, 0.3, 0.3))

    if prob:
        node.hide = False

    return not prob


def remlink(ng, links):
    for link in links:
        ng.links.remove(link)


def sockhide(node, lsocknames):
    try:
        for ins in [insock for insock in node.inputs if insock.name in lsocknames]:
            node.outputs[ins.name].hide = True if ins.links else False
        for outs in [outsock for outsock in node.outputs if outsock.name in lsocknames]:
            node.inputs[outs.name].hide = True if outs.links else False
    except Exception as e:
        print('sockhide', e)


def socklink(sock, ng):
    if ng in [g.name for g in bpy.data.node_groups]:
        try:
            valid1 = sock.valid if not sock.get('valid') else sock['valid']

            for link in sock.links:
                valid2 = link.to_socket.valid if not link.to_socket.get('valid') else link.to_socket['valid']
                valset = set(valid1) & set(valid2)

                if not valset or len(valset) < min((len(valid1), len(valid2))):
                    bpy.data.node_groups[ng].links.remove(link)

        except Exception:
            if sock.links:
                bpy.data.node_groups[ng].links.remove(sock.links[-1])


def socklink2(sock, ng):
    try:
        valid1 = sock.ret_valid(sock.node)

        for link in sock.links:
            valid2 = link.to_socket.ret_valid(link.to_socket.node)
            valset = set(valid1) & set(valid2)

            if not valset:
                ng.links.remove(link)

    except Exception:
        if sock.links:
            ng.links.remove(sock.links[-1])


def uvsocklink(sock, ng):
    try:
        uv1 = sock.uvalue
        for link in sock.links:
            uv2 = link.to_socket.uvalue
            if uv1 != uv2:
                bpy.data.node_groups[ng].links.remove(link)
    except Exception:
        pass


def uvsocklink2(sock, ng):
    try:
        uv1 = sock.uvalue
        for link in sock.links:
            uv2 = link.to_socket.uvalue
            if uv1 != uv2:
                ng.links.remove(link)
    except Exception:
        pass


def rettimes(ts, fs, us):
    tot = range(min(len(ts), len(fs), len(us)))
    fstrings, ustrings, tstrings = [[] for t in tot], [[] for t in tot], ['Through: {}/{}'.format(dtdf(ts[t]).month, dtdf(ts[t]).day) for t in tot]

    for t in tot:
        fstrings[t] = ['For: ' + ''.join(f.strip()) for f in fs[t].split(' ') if f.strip(' ') != '']

        for uf, ufor in enumerate(us[t].split(';')):
            ustrings[t].append([])

            for ut, utime in enumerate(ufor.split(',')):
                ustrings[t][uf].append(['Until: ' + ','.join([u.strip() for u in utime.split(' ') if u.strip(' ')])])

    return (tstrings, fstrings, ustrings)


def retdates(sdoy, edoy, y):
    (y1, y2) = (y, y) if edoy >= sdoy else (y - 1, y)
    sdate = datetime.datetime(y1, 1, 1) + datetime.timedelta(sdoy - 1)
    edate = datetime.datetime(y2, 1, 1) + datetime.timedelta(edoy - 1)
    return (sdate, edate)


def ret_param(param, val):
    if isinstance(param, float):
        return float(val)
    elif isinstance(param, int):
        return int(val)
    else:
        return str(val)


def sunposh(context, suns):
    scene = context.scene
    sps = [solarPosition(scene.solday, i, scene.latitude, scene.longitude) for i in range(0, 23)]
    spsvalid = [sp[0] > 0 for sp in sps]

    if sum(spsvalid) > len(suns):
        for i in range(sum(spsvalid) - len(suns)):
            bpy.ops.object.lamp_add(type="SUN")
            suns.append(context.active_object)
    elif sum(spsvalid) < len(suns):
        [scene.objects.unlink(sun) for sun in suns[sum(spsvalid)]]

    for sun in suns:
        pass


def sunapply(scene, sun, values, solposs, frames, sdist):
    sun.data.animation_data_clear()
    sun.animation_data_clear()
    action = bpy.data.actions.get("EnVi Sun")

    if action:
        bpy.data.actions.remove(action, do_unlink=True)

    
    action = bpy.data.actions.new(name="EnVi Sun")
    sun.animation_data_create().action = action
    # es_slot = sun.animation_data.action.slots.new(id_type='OBJECT', name=sun.name)
    # es_layer = sun.animation_data.action.layers.new("Layer")
    # strip = es_layer.strips.new(type='KEYFRAME')
    # channelbag = strip.channelbag(es_slot, ensure=True)

    sunposx = action.fcurve_ensure_for_datablock(sun, "location", index=0)
    sunposy = action.fcurve_ensure_for_datablock(sun, "location", index=1)
    sunposz = action.fcurve_ensure_for_datablock(sun, "location", index=2)
    sunposx.keyframe_points.add(len(frames))
    sunposy.keyframe_points.add(len(frames))
    sunposz.keyframe_points.add(len(frames))
    sunrotx = action.fcurve_ensure_for_datablock(sun, "rotation_euler", index=0)
    sunroty = action.fcurve_ensure_for_datablock(sun, "rotation_euler", index=1)
    sunrotz = action.fcurve_ensure_for_datablock(sun, "rotation_euler", index=2)
    sunrotx.keyframe_points.add(len(frames))
    sunroty.keyframe_points.add(len(frames))
    sunrotz.keyframe_points.add(len(frames))
    sunenergy = action.fcurve_ensure_for_datablock(sun, "energy")
    sunenergy.keyframe_points.add(len(frames))

# This is an attempt to use low level routines for node value animation but it don't work.
    if sun.data.node_tree:
        sun.data.node_tree.animation_data_clear()
        action = bpy.data.actions.get("EnVi Sun Node")

        if action:
            bpy.data.actions.remove(action, do_unlink=True)

        sun.data.node_tree.animation_data_create()
        sun.data.node_tree.animation_data.action = bpy.data.actions.new(name="EnVi Sun Node")
        emnodes = [emnode for emnode in sun.data.node_tree.nodes if emnode.bl_label == 'Emission']
        bbnodes = [bbnode for bbnode in sun.data.node_tree.nodes if bbnode.bl_label == 'Blackbody']

        for emnode in emnodes:
            em1 = sun.data.node_tree.animation_data.action.fcurves.new(data_path='nodes["{}"].inputs[1].default_value'.format(emnode.name))
            em1.keyframe_points.add(len(frames))

        for bbnode in bbnodes:
            bb1 = sun.data.node_tree.animation_data.action.fcurves.new(data_path='nodes["{}"].inputs[0].default_value'.format(bbnode.name))
            bb1.keyframe_points.add(len(frames))

    if scene.world.node_tree:
        scene.world.node_tree.animation_data_clear()
        action = bpy.data.actions.get("EnVi World Node")

        if action:
            bpy.data.actions.remove(action, do_unlink=True)

        scene.world.node_tree.animation_data_create()
        scene.world.node_tree.animation_data.action = bpy.data.actions.new(name="EnVi World Node")
        stnodes = [stnode for stnode in scene.world.node_tree.nodes if stnode.bl_label == 'Sky Texture']

        for stnode in stnodes:
            st1x = scene.world.node_tree.animation_data.action.fcurves.new(data_path='nodes["{}"].sun_direction'.format(stnode.name), index=0)
            st1y = scene.world.node_tree.animation_data.action.fcurves.new(data_path='nodes["{}"].sun_direction'.format(stnode.name), index=1)
            st1z = scene.world.node_tree.animation_data.action.fcurves.new(data_path='nodes["{}"].sun_direction'.format(stnode.name), index=2)
            st1x.keyframe_points.add(len(frames))
            st1y.keyframe_points.add(len(frames))
            st1z.keyframe_points.add(len(frames))

    for f, frame in enumerate(frames):
        (sun.data.shadow_soft_size, sun.data.energy) = values[f][:2]
        sunz = sdist * sin(solposs[f][2])
        sunx = -(sdist**2 - (sunz)**2)**0.5 * sin(solposs[f][3])
        suny = -(sdist**2 - (sunz)**2)**0.5 * cos(solposs[f][3])
        sunpos = [sunx, suny, sunz]
        sunrot = [(pi / 2) - solposs[f][2], 0, -solposs[f][3]]
        scene.display.light_direction = (-sin(solposs[f][3]) * cos(solposs[f][2]), sin(solposs[f][2]), cos(solposs[f][3]) * cos(solposs[f][2]))

        if scene.render.engine == 'CYCLES' and scene.world.node_tree:
            if 'Sky Texture' in [no.bl_label for no in scene.world.node_tree.nodes]:
                skydir = -sin(solposs[f][3]), -cos(solposs[f][3]), sin(solposs[f][2])
                st1x.keyframe_points[f].co = frame, skydir[0]
                st1y.keyframe_points[f].co = frame, skydir[1]
                st1z.keyframe_points[f].co = frame, skydir[2]

        if scene.render.engine == 'CYCLES' and sun.data.node_tree:
            for emnode in emnodes:
                em1.keyframe_points[f].co = frame, values[f][1]
        if sun.data.node_tree:
            for bbnode in bbnodes:
                bb1.keyframe_points[f].co = frame, retsunct(solposs[f][2])

        sunposx.keyframe_points[f].co = frame, sunpos[0]
        sunposy.keyframe_points[f].co = frame, sunpos[1]
        sunposz.keyframe_points[f].co = frame, sunpos[2]
        sunrotx.keyframe_points[f].co = frame, sunrot[0]
        sunroty.keyframe_points[f].co = frame, sunrot[1]
        sunrotz.keyframe_points[f].co = frame, sunrot[2]
        sunenergy.keyframe_points[f].co = frame, values[f][1]

    sun.data.cycles.use_multiple_importance_sampling = True


def retsunct(beta):
    return 2500 + 3000 * sin(beta)**0.5 if beta > 0 else 2500


def spfc(self):
    scene = bpy.context.scene
    svp = scene.vi_params

    if not svp['viparams'].get('newframe'):
        svp['viparams']['newframe'] = 1
    else:
        svp['viparams']['newframe'] = 0
        scene.frame_set(scene.frame_current)

    if svp['viparams']['resnode'] == 'VI Sun Path':
        spoblist = {ob.get('VIType'): ob for ob in scene.objects if ob.get('VIType') in ('Sun', 'SPathMesh')}
        beta, phi = solarPosition(scene.sp_sd, scene.sp_sh, scene.latitude, scene.longitude)[2:]

        if not scene.world.use_nodes:
            scene.world.use_nodes = True

        nt = bpy.data.worlds[0].node_tree

        if nt and nt.nodes.get('Sky Texture'):
            scene.world.node_tree.nodes['Sky Texture'].sun_direction = -sin(phi), -cos(phi), sin(beta)

        for ob in scene.objects:
            if ob.get('VIType') == 'Sun':
                ob.rotation_euler = pi * 0.5 - beta, 0, -phi
                ob.location.z = spoblist['SPathMesh'].location.z + 100 * sin(beta)
                ob.location.x = spoblist['SPathMesh'].location.x - (100**2 - (spoblist['Sun'].location.z - spoblist['SPathMesh'].location.z)**2)**0.5 * sin(phi)
                ob.location.y = spoblist['SPathMesh'].location.y - (100**2 - (spoblist['Sun'].location.z - spoblist['SPathMesh'].location.z)**2)**0.5 * cos(phi)

                if ob.data.node_tree:
                    for blnode in [blnode for blnode in ob.data.node_tree.nodes if blnode.bl_label == 'Blackbody']:
                        blnode.inputs[0].default_value = 2500 + 3000 * sin(beta)**0.5
                    for emnode in [emnode for emnode in ob.data.node_tree.nodes if emnode.bl_label == 'Emission']:
                        emnode.inputs[1].default_value = 10 * sin(beta)

            elif ob.get('VIType') == 'SkyMesh':
                ont = ob.data.materials['SkyMesh'].node_tree
                if ont and ont.nodes.get('Sky Texture'):
                    ont.nodes['Sky Texture'].sun_direction = sin(phi), -cos(phi), sin(beta)

            elif ob.get('VIType') == 'SunMesh':
                ob.location = (0, 0, 0)

                if ob.data.materials[0].node_tree:
                    for smblnode in [smblnode for smblnode in ob.data.materials[0].node_tree.nodes if ob.data.materials and smblnode.bl_label == 'Blackbody']:
                        smblnode.inputs[0].default_value = 2500 + 3000 * sin(beta)**0.5
    else:
        return


def bm_to_stl(bm, stl_path):
    with open(stl_path, 'w') as stlfile:
        stlfile.write('solid\n')

        for face in bm.faces:
            stlfile.write('facet normal {0[0]:.6f} {0[1]:.6f} {0[2]:.6f}\nouter loop\n'.format(face.normal.normalized()))

            for vert in face.verts:
                stlfile.write('vertex {0[0]:.6f} {0[1]:.6f} {0[2]:.6f}\n'.format(vert.co))

            stlfile.write('endloop\nendfacet\n')

        stlfile.write('endsolid\n')

    bm.free()


def meshes_to_solids(context, coll):
    dp = context.evaluated_depsgraph_get()
    scene = context.scene
    svp = scene.vi_params
    c_geos = [ob for ob in coll.objects if ob.visible_get() and ob.type == "MESH"]
    g_geos = []
    g_names = []
    fi = 1
    bb_mesh = ret_coll_bb(coll)
    bm = bmesh.new()
    bm.from_mesh(bb_mesh)
    lbm = len(bm.faces)
    faces = []

    for fi, face in enumerate(bm.faces[:lbm]):
        edges = [occ.Segment(occ.gp_Pnt(tuple(loop.vert.co)), occ.gp_Pnt(tuple(loop.link_loop_next.vert.co))) for loop in face.loops]
        wire = occ.Wire(edges)
        f = occ.Face(wire)
        f.name = None
        faces.append(f)
        fi += 1

    bm.free()
    d_geo = occ.OCCGeometry(occ.Compound(faces))
    d_geo.Heal(tolerance=0.01)

    for oi, ob in enumerate(c_geos):
        fi = 0
        faces = []
        bm = bmesh.new()
        bm.from_object(ob, dp)
        bmesh.ops.triangulate(bm, faces=bm.faces)
        bm.transform(ob.matrix_world)
        lbm = len(bm.faces)

        for fi, face in enumerate(bm.faces[:lbm]):
            edges = [occ.Segment(occ.gp_Pnt(tuple(round(co, 4) for co in loop.vert.co)), occ.gp_Pnt(tuple(round(co, 4) for co in loop.link_loop_next.vert.co))) for loop in face.loops]
            wire = occ.Wire(edges)
            f = occ.Face(wire)

            if len(f.edges) > 2 and ob.material_slots and ob.material_slots[face.material_index].material:
                matname = ob.material_slots[face.material_index].material.name
                f.name = matname
                f.mat(matname)
                f.bc(matname)
                f.layer = oi

            faces.append(f)
            fi += 1

        g_geo = occ.OCCGeometry(occ.Compound(faces))
        fns = [face.name for face in g_geo.faces]
        fcs = [face.center for face in g_geo.faces]
        print(f'Healing {ob.name}')
        g_geo.Heal(tolerance=0.001)

        if None in set([face.name for face in g_geo.shape.faces]):
            for fi, face in enumerate(g_geo.shape.faces):
                if face.name is None:
                    for fci, fc in enumerate(fcs):
                        if (mathutils.Vector(face.center) - mathutils.Vector(fc)).length < 0.0001:
                            face.name = fns[fci]
                            break
                    if face.name is None:
                        face.name = fns[fi]

        g_geos.append(g_geo)
        g_names.append(ob.name)
        bm.free()

    solids = []

    for gi, g_geo in enumerate(g_geos):
        if len(g_geo.shape.SubShapes(occ.SOLID)) == len(g_geo.shape.SubShapes(occ.SHELL)):
            for g_geo_solid in g_geo.shape.SubShapes(occ.SOLID):
                g_geo_solid.name = g_names[gi]
                solids.append(g_geo_solid)

        else:
            for g_shell in g_geo.shape.SubShapes(occ.SHELL):
                g_shell_solid = occ.OCCGeometry(g_shell)
                fns = [face.name for face in g_shell_solid.faces]
                fcs = [face.center for face in g_shell_solid.faces]
                g_shell_solid.Heal(tolerance=0.0001)

                if None in set([face.name for face in g_shell_solid.shape.faces]):
                    for fi, face in enumerate(g_shell_solid.shape.faces):
                        if face.name is None:
                            for fci, fc in enumerate(fcs):
                                if (mathutils.Vector(face.center) - mathutils.Vector(fc)).length < 0.0001:
                                    face.name = fns[fci]
                                    break
                            if face.name is None:
                                face.name = fns[fi]

                for gs, g_solid in enumerate(g_shell_solid.shape.SubShapes(occ.SOLID)):
                    solids.append(g_solid)
                    print('shell to solid', g_names[gi])

                if not g_shell_solid.shape.SubShapes(occ.SOLID):
                    g_shell_solid.shape.WriteStep(os.path.join(svp['viparams']['newdir'], 'temp.step'))
                    g_shell_solid = occ.OCCGeometry(os.path.join(svp['viparams']['newdir'], 'temp.step'))
                    fns = [face.name for face in g_shell_solid.faces]
                    fcs = [face.center for face in g_shell_solid.faces]
                    g_shell_solid.Heal(tolerance=0.0001)

                    if None in set([face.name for face in g_shell_solid.shape.faces]):
                        for fi, face in enumerate(g_shell_solid.shape.faces):
                            if face.name is None:
                                for fci, fc in enumerate(fcs):
                                    if (mathutils.Vector(face.center) - mathutils.Vector(fc)).length < 0.0001:
                                        face.name = fns[fci]
                                        break
                                if face.name is None:
                                    face.name = fns[fi]

                    for gs, g_solid in enumerate(g_shell_solid.shape.SubShapes(occ.SOLID)):
                        solids.append(g_solid)

    solids = [occ.Fuse(solids)]

    for si, solid in enumerate(solids):
        d_geo = d_geo.shape - solid if not si else d_geo - solid
        fns = [face.name for face in d_geo.faces]
        fcs = [face.center for face in d_geo.faces]

        if None in set([face.name for face in d_geo.faces]):
            for fi, face in enumerate(d_geo.faces):
                if face.name is None:
                    for fci, fc in enumerate(fcs):
                        if (mathutils.Vector(face.center) - mathutils.Vector(fc)).length < 0.0001:
                            face.name = fns[fci]
                            break

                    if face.name is None:
                        face.name = fns[fi]

    return d_geo.solids[1:]

def solid_to_mesh(svp, solid, si):
    if any([len(face.wires) > 1 for face in solid.faces]):
        fis, verts, mat_list, vno = [], [], [], 0
        solid.MakeTriangulation()
        solid.WriteBrep(os.path.join(svp['viparams']['newdir'], 'temp.brep'), withTriangles=1)

        with open(os.path.join(svp['viparams']['newdir'], 'temp.brep')) as brep_file:
            brep_lines = brep_file.readlines()
        
        for li, line in enumerate(brep_lines):
            if line.split() and line.split()[0] == 'Triangulations':
                t_lines = brep_lines[li + 1: li + 1 + 2 * int(line.split()[1])]
                break
        
        for li, line in enumerate(t_lines):
            if not li % 2:
                verts.append(t_lines[li + 1].split()[:int(line.split()[0]) * 3])
                fis.append([int(f) - 1 + vno for f in t_lines[li + 1].split()[-int(line.split()[1]) * 3:]])
                vno += int(line.split()[0])

                for m in range(int(line.split()[1])):
                    mat_list.append(solid.faces[int(li/2)].name)

        fis = [int(f) for fi in fis for f in fi]
        step_fs = [[fis[i], fis[i + 1], fis[i + 2]] for i in range(len(fis)) if not i % 3]
        vps = [float(v) for p in verts for v in p]
        step_vs = [[vps[i], vps[i + 1], vps[i + 2]] for i in range(len(vps)) if not i % 3]
        
    else:   
        all_vps = []
        mat_list = [face.name for face in solid.faces]

        if mat_list and len(set(mat_list)) > 1:
            for fi, face in enumerate(solid.faces):
                vps = [tuple(v.p) for v in face.wires[0].vertices]

                if not len(vps) % 2:
                    if vps[0] not in (vps[-2], vps[-1]):
                        vps[0], vps[1] = vps[1], vps[0]

                    for vi in [vi for vi in range(1, len(vps) - 1) if vi % 2]:
                        if vps[vi] != vps[vi + 1]:
                            vps[vi + 1], vps[vi + 2] = vps[vi + 2], vps[vi + 1]

                    if not all([vps[vi] == vps[vi+1] for vi in range(len(vps) -1) if vi % 2]):
                        evps = [[tuple([round(p, 5) for p in e.start]), tuple([round(p, 5) for p in e.end])] for e in face.edges]

                        for evi in range(len(evps) - 1):
                            if evps[evi] not in evps[evi+1] and evps[evi] not in evps[evi+1]:
                                for vi in range(evi + 2, len(evps)):
                                    if evps[evi][0] in evps[vi] or evps[evi][1] in evps[vi]:
                                        evps[evi + 1], evps[vi] = evps[vi], evps[evi + 1]
                                        break

                        for evi in range(len(evps) - 1):
                            if evps[evi][0] in evps[evi + 1]:
                                evps[evi][0], evps[evi][1] = evps[evi][1], evps[evi][0]
                            if evps[evi + 1][1] in evps[evi]:
                                evps[evi + 1][0], evps[evi + 1][1] = evps[evi + 1][1], evps[evi + 1][0]

                        vps = [x for xs in evps for x in xs]

                    all_vps.append([vps[vi] for vi in range(len(vps)) if not vi % 2])

            fis = [len(step_v) for step_v in all_vps]
            step_vs = [x for xs in all_vps for x in xs ]
            step_fs = []
            f = 0

            for fi in fis:
                step_fs.append([i for i in range(f, f + int(fi))])
                f += int(fi)

        elif len(set(mat_list)) < 2:
            logentry('Error: Solids create a volume with only one material type. Create or duplicate an additional material on the zone boundary')
    
    if f'Zone{si}' not in [mesh.name for mesh in bpy.data.meshes]:
        mesh = bpy.data.meshes.new(f'Zone{si}')
    else:
        mesh = bpy.data.meshes[f'Zone{si}']
        bm = bmesh.new()
        bm.to_mesh(mesh)
        bm.free()

    mesh.from_pydata(step_vs, [], step_fs, shade_flat=True)
    bm = bmesh.new()
    bm.from_mesh(mesh)
    bmesh.ops.remove_doubles(bm, verts=bm.verts, dist = 0.00001)
    manifold = all([e.is_manifold for e in bm.edges])

    if not manifold:
        bmesh.ops.recalc_face_normals(bm, faces=bm.faces)

    manifold = all([e.is_manifold for e in bm.edges])

    for mi, mat in enumerate(set(mat_list)):
        if mat:
            if len(mesh.materials) < mi + 1:
                mesh.materials.append(bpy.data.materials[mat])
            else:
                mesh.materials[mi] = bpy.data.materials[mat]

    for poly in bm.faces:
        poly.material_index = list(set(mat_list)).index(mat_list[poly.index])

    bmesh.ops.recalc_face_normals(bm, faces=bm.faces)
    bmesh.ops.dissolve_limit(bm, angle_limit=0.1, use_dissolve_boundaries=False, verts=bm.verts, edges=bm.edges, delimit={'NORMAL', 'MATERIAL'})
    bm.to_mesh(mesh)
    bm.free()
    return manifold, mesh

def ob_to_stl(self, dp, stl_path):
    o = self.id_data

    if o.visible_get():
        bm = bmesh.new()
        bm.from_object(o, dp)
        bm.transform(o.matrix_world)
        bmesh.ops.triangulate(bm, faces=bm.faces, quad_method='BEAUTY', ngon_method='BEAUTY')

        with open(stl_path, 'w') as stlfile:
            stlfile.write('solid\n')

            for face in bm.faces:
                stlfile.write('facet normal {0[0]:.6f} {0[1]:.6f} {0[2]:.6f}\nouter loop\n'.format(face.normal.normalized()))

                for vert in face.verts:
                    stlfile.write('vertex {0[0]:.6f} {0[1]:.6f} {0[2]:.6f}\n'.format(vert.co))

                stlfile.write('endloop\nendfacet\n')

            stlfile.write('endsolid\n')

        bm.free()


def find_materials_in_groupinstances(empty):
    if empty.instance_collection.name in checked_groups_names_list:
        return None

    for obj in bpy.data.collections[empty.instance_collection.name].objects:
        if obj.instance_type == 'COLLECTION' and obj.instance_collection is not None and obj.type == 'EMPTY':
            return find_materials_in_groupinstances(obj)
        elif obj.type == "MESH":
            for slot in obj.material_slots:
                if slot.material:
                    materials_from_group.add(slot.material)

    checked_groups_names_list.append(empty.instance_collection.name)  # or no empty mat in group
    return None


def material_on_sel_obj(mat):
    obj = bpy.context.active_object

    for slot in obj.material_slots:
        if slot.material == mat:
            return True

    return False


def get_materials():
    materials = []

    for mat in bpy.data.materials:
        conditions = [mat.vi_params.envi_nodes]

        if all(conditions):
            materials.append(mat)

    additional_mats = set()
    checked_groups_names_list.clear()

    for obj in bpy.context.selected_objects:
        if obj.instance_type == 'COLLECTION' and obj.instance_collection is not None and obj.type == 'EMPTY':
            find_materials_in_groupinstances(obj)
            additional_mats = additional_mats | materials_from_group
            materials_from_group.clear()

    all_mats = list(set(materials) | additional_mats)
    all_mats = sorted(all_mats, key=lambda x: x.name.lower())
    return all_mats


# def au_pb(file):
#     kivytext = "# -*- coding: "+sys.getfilesystemencoding()+" -*-\n\
# from kivy.app import App \n\
# from kivy.clock import Clock \n\
# from kivy.uix.progressbar import ProgressBar\n\
# from kivy.uix.boxlayout import BoxLayout\n\
# from kivy.uix.button import Button\n\
# from kivy.uix.label import Label\n\
# from kivy.config import Config\n\
# Config.set('graphics', 'width', '500')\n\
# Config.set('graphics', 'height', '200')\n\
# Config.set('kivy', 'log_level', 'warning')\n\
# \n\
# class CancelButton(Button):\n\
#     def on_touch_down(self, touch):\n\
#         if 'button' in touch.profile:\n\
#             if self.collide_point(*touch.pos):\n\
#                 App.get_running_app().stop()\n\
#         else:\n\
#             return\n\
#     def on_open(self, widget, parent):\n\
#         self.focus = True\n\
# \n\
# class Calculating(App):\n\
#     bl = BoxLayout(orientation='vertical')\n\
#     rpb = ProgressBar()\n\
#     button = CancelButton(text='Cancel', font_size=20)\n\
#     bl.add_widget(rpb)\n\
#     bl.add_widget(button)\n\
# \n\
#     def build(self):\n\
#         self.title = 'Calculating RTs'\n\
#         refresh_time = 1\n\
#         return self.bl\n\
# \n\
# if __name__ == '__main__':\n\
#     Calculating().run()\n"

#     with open(file+".py", 'w') as kivyfile:
#         kivyfile.write(kivytext)
#     return Popen([sys.executable, file+".py"])

# def fvprogressbar(file, et, residuals, frame):
#     addonpath = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
#     kivytext = "# -*- coding: "+sys.getfilesystemencoding()+" -*-\n\
# import os, sys\n\
# sys.path.append(os.path.join(r'"+addonpath+"', 'Python', sys.platform))\n\
# from kivy.app import App\n\
# from kivy.clock import Clock\n\
# from kivy.uix.progressbar import ProgressBar\n\
# from kivy.uix.boxlayout import BoxLayout\n\
# from kivy.uix.gridlayout import GridLayout\n\
# from kivy.uix.button import Button\n\
# from kivy.uix.label import Label\n\
# from kivy.config import Config\n\
# Config.set('graphics', 'width', '500')\n\
# Config.set('graphics', 'height', '250')\n\
# \n\
# class CancelButton(Button):\n\
#     def on_touch_down(self, touch):\n\
#         if 'button' in touch.profile:\n\
#             if self.collide_point(*touch.pos):\n\
#                 App.get_running_app().stop()\n\
#         else:\n\
#             return\n\
#     def on_open(self, widget, parent):\n\
#         self.focus = True\n\
# \n\
# class Calculating(App):\n\
#     rpbs, labels, nums, oftime = [], [], [], '0'\n\
#     bl = BoxLayout(orientation='vertical')\n\
#     gl  = GridLayout(cols=3, height = 250)\n\
#     t = Label(text='Time:', font_size=20, size_hint=(0.2, .2))\n\
#     tpb = ProgressBar(max = "+str(et)+")\n\
#     tt = Label(text=oftime, font_size=20, size_hint=(0.2, .2))\n\
#     gl.add_widget(t)\n\
#     gl.add_widget(tpb)\n\
#     gl.add_widget(tt)\n\
# \n\
#     for r in "+residuals+":\n\
#         rpb = ProgressBar(max = 1)\n\
#         rpbs.append(rpb)\n\
#         label = Label(text=r, font_size=20, size_hint=(0.2, .2))\n\
#         num = Label(text='1', font_size=20, size_hint=(0.2, .2))\n\
#         labels.append(r)\n\
#         nums.append(num)\n\
#         gl.add_widget(label)\n\
#         gl.add_widget(rpb)\n\
#         gl.add_widget(num)\n\
#     bl.add_widget(gl)\n\
#     button = CancelButton(text='Cancel', font_size=20, size_hint=(1, .2))\n\
#     bl.add_widget(button)\n\
# \n\
#     def build(self):\n\
#         self.title = 'OpenFOAM Residuals - "+str(frame)+"'\n\
#         refresh_time = 1\n\
#         Clock.schedule_interval(self.timer, refresh_time)\n\
#         return self.bl\n\
# \n\
#     def timer(self, dt):\n\
#         with open(r'"+file+"', 'r') as pffile:\n\
#             for ri, r in enumerate(pffile.readlines()):\n\
#                 if r.split()[0] in ('Time', 'Ux', 'Uy', 'Uz', 'p', 'k', 'epsilon', 'p_rgh', 'e'):\n\
#                     try:\n\
#                         if r.split()[0] == 'Time':\n\
#                             self.tpb.value = float(r.split()[1])\n\
#                             self.tt.text = '{:.1f}'.format(float(r.split()[1]))\n\
#                         else:\n\
#                             li = self.labels.index(r.split()[0])\n\
#                             self.rpbs[li].value = abs(float(r.split()[1]))**0.5\n\
#                             self.nums[li].text = '{:.5f}'.format(abs(float(r.split()[1])))\n\
#                     except Exception as e:\n\
#                         pass\n\
# \n\
# if __name__ == '__main__':\n\
#     Calculating().run()"

#     with open(file+".py", 'w') as kivyfile:
#         kivyfile.write(kivytext)

#     return Popen([sys.executable, file+".py"])

#    txts = []
#    come = bpy.data.meshes.new("Compass")
#    coo = bpy.data.objects.new('Compass', come)
#   coo.location = loc
#    bpy.context.scene.collection.objects.link(coo)
#    bpy.context.view_layer.objects.active = coo
#    bpy.ops.object.material_slot_add()

#    bm = bmesh.new()
#    matrot = Matrix.Rotation(pi*0.25, 4, 'Z')
#    bmesh.ops.create_circle(bm, cap_ends=True, radius=100, segments=132,  matrix=Matrix.Rotation(0, 4, 'Z')@Matrix.Translation((0, 0, 0)))

#    newgeo = bmesh.ops.extrude_edge_only(bm, edges = bm.edges, use_select_history=False)

#    for face in [f for f in newgeo['geom'] if isinstance(f, bmesh.types.BMFace)]:
#        face.material_index = 1
#    newverts0 = []
#    for v, vert in enumerate([v for v in newgeo['geom'] if isinstance(v, bmesh.types.BMVert)]):
#        vert.co = vert.co + (vert.co - coo.location).normalized() * scale * 0.0025
#        vert.co[2] = 0
#        newverts0.append(vert)
#
#    newgeo = bmesh.ops.extrude_edge_only(bm, edges = [e for e in newgeo['geom'] if isinstance(e, bmesh.types.BMEdge) and e.calc_length() > 0.05], use_select_history=False)
#
#    for face in [f for f in newgeo['geom'] if isinstance(f, bmesh.types.BMFace)]:
#        face.material_index = 0
#
#    newverts = []
#    for v, vert in enumerate([v for v in newgeo['geom'] if isinstance(v, bmesh.types.BMVert)]):
#        vert.co = vert.co + (vert.co - coo.location).normalized() * scale * 0.05
#        vert.co[2] = 0
#        newverts.append(vert)
#
#    for v in newverts:
#        if abs(v.co[1]) < 0.01:
#            v.co[0] += v.co[0] * 0.025
#        elif abs(v.co[0]) < 0.01:
#            v.co[1] += v.co[1] * 0.025
#
#    newgeo = bmesh.ops.extrude_edge_only(bm, edges = [e for e in newgeo['geom'] if isinstance(e, bmesh.types.BMEdge) and e.verts[0] in newverts and e.verts[1] in newverts], use_select_history=False)
#
#    for face in [f for f in newgeo['geom'] if isinstance(f, bmesh.types.BMFace)]:
#        face.material_index = 1
#
#    for v, vert in enumerate([v for v in newgeo['geom'] if isinstance(v, bmesh.types.BMVert)]):
#        vert.co = vert.co + (vert.co - coo.location).normalized() * scale * 0.0025
#        vert.co[2] = 0
#
#    bmesh.ops.dissolve_edges(bm, edges = [e for e in bm.edges if e.verts[0] in (newverts0 + newverts) and e.verts[1] in (newverts0 + newverts) and len(e.link_faces) == 2 and abs(e.link_faces[0].calc_area() - e.link_faces[1].calc_area()) < 1.0],
#                                          use_verts = True, use_face_split = False)
#
#    #matrot = Matrix.Rotation(pi*0.25, 4, 'Z')
#    degmatrot = Matrix.Rotation(pi*0.125, 4, 'Z')
#    tmatrot = Matrix.Rotation(0, 4, 'Z')
#    direc = Vector((0, 1, 0))
#    tmatrot = Matrix.Rotation(0, 4, 'Z')
#    f_sizes = (0.08, 0.06, 0.06, 0.06, 0.08, 0.06, 0.06, 0.06, 0.08, 0.06, 0.06, 0.06, 0.08, 0.06, 0.06, 0.06)
#    f_texts = ('N', u'337.5\u00B0', u'315\u00B0', u'292.5\u00B0', 'W', u'247.5\u00B0', u'225\u00B0', u'202.5\u00B0', 'S', u'157.5\u00B0', u'135\u00B0', u'112.5\u00B0', 'E', u'67.5\u00B0', u'45\u00B0', u'22.5\u00B0')
#    f_texts = ('N', 'NNW', 'NW', 'WNW', 'W', 'WSW', 'SW', 'SSW', 'S', 'SSE', 'SE', 'ESE', 'E', 'ENE', 'NE', 'NNE')
#    f_texts = ('N', u'337.5\u00B0', 'NW', u'292.5\u00B0', 'W', u'247.5\u00B0', 'SW', u'202.5\u00B0', 'S', u'157.5\u00B0', 'SE', u'112.5\u00B0', 'E', u'67.5\u00B0', 'NE', u'22.5\u00B0')
#
#    for d in range(16):
#        bpy.ops.object.text_add(align='WORLD', enter_editmode=False, location=Vector(loc) + scale*1.0005*(tmatrot@direc), rotation=tmatrot.to_euler())
#        txt = bpy.context.active_object
#        txt.data.font = bpy.data.fonts.load(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'Fonts/NotoSans-Regular.ttf'))
#        txt.scale, txt.data.body, txt.data.align_x, txt.data.align_y, txt.location[2]  = (scale*f_sizes[d], scale*f_sizes[d], scale*f_sizes[d]), f_texts[d], 'CENTER', 'BOTTOM', 0.05
#        bpy.ops.object.convert(target='MESH')
#        bpy.ops.object.material_slot_add()
#        txt.material_slots[-1].material = basemat
#        txts.append(txt)
#        tmatrot = tmatrot @ degmatrot
#
#    bm.to_mesh(come)
#    bm.free()

# def set_legend(ax):
#    l = ax.legend(borderaxespad = -4)
#    plt.setp(l.get_texts(), fontsize=8)

# def wr_axes(plt):
#    fig = plt.figure(figsize=(8, 8), dpi=150, facecolor='w', edgecolor='w')
#    rect = [0.1, 0.1, 0.8, 0.8]
#    ax = WindroseAxes(fig, rect, facecolor='w')
#    fig.add_axes(ax)
#    return(fig, ax)
# def drawloop(x1, y1, x2, y2):
#     bgl.glLineWidth(1)
#     bgl.glColor4f(0.0, 0.0, 0.0, 1.0)
#     bgl.glBegin(bgl.GL_LINE_LOOP)
#     bgl.glVertex2i(x1, y2)
#     bgl.glVertex2i(x2, y2)
#     bgl.glVertex2i(x2, y1)
#     bgl.glVertex2i(x1, y1)
#     bgl.glEnd()


# def drawsquare(c, w, h, col):
#     vxs = (c[0] + 0.5 * w, c[0] + 0.5 * w, c[0] - 0.5 * w, c[0] - 0.5 * w)
#     vys = (c[1] - 0.5 * h, c[1] + 0.5 * h, c[1] + 0.5 * h, c[1] - 0.5 * h)

#     if col:
#         bgl.glColor4f(*col)
#         z = 0.1
#         bgl.glBegin(bgl.GL_POLYGON)
#     else:
#         z = 0.05
#         bgl.glLineWidth(1)
#         bgl.glColor4f(0, 0, 0, 1)
#         bgl.glBegin(bgl.GL_LINE_LOOP)
#     for v in range(4):
#         bgl.glVertex3f(vxs[v], vys[v], z)
#     bgl.glLineWidth(1)
#     bgl.glColor4f(0, 0, 0, 1)
#     bgl.glEnd()

    # bm.transform(o.matrix_world.inverted())
    # bmesh.ops.recalc_face_normals(bm, faces = bm.faces)
    # bm.normal_update()

    # tris = bm.calc_loop_triangles()

    # with open(stl_path, 'w') as stlfile:
    #     stlfile.write('solid\n')

    #     for tri in tris:
    #         # Below can correct post-processing normal but creates incorrect STLs
    #         blender_normal = (o.matrix_world.inverted_safe().transposed().to_3x3() @ tri[0].face.normal).normalized()
    #         face_normal = (o.rotation_euler.to_matrix() @ tri[0].face.normal).to_3d().normalized()
    #         stlfile.write('facet normal {0[0]:.3f} {0[1]:.3f} {0[2]:.3f}\nouter loop\n'.format((o.matrix_world.inverted_safe().transposed().to_3x3() @ tri[0].face.normal).normalized()))
    #         # stlfile.write('facet normal {0[0]:.3f} {0[1]:.3f} {0[2]:.3f}\nouter loop\n'.format(tri[0].face.normal.normalized()))

    #         if face_normal.angle(blender_normal) < 0.1 or o.vi_params.vi_type in ('2', '3'):
    #             for t in tri:
    #                 stlfile.write('vertex {0[0]:.5f} {0[1]:.5f} {0[2]:.5f}\n'.format(o.matrix_world@t.vert.co))
    #         else:
    #             for t in tri[::-1]:
    #                 stlfile.write('vertex {0[0]:.5f} {0[1]:.5f} {0[2]:.5f}\n'.format(o.matrix_world@t.vert.co))

    #         stlfile.write('endloop\nendfacet\n')
    #     stlfile.write('endsolid\n')
    # bm.transform(o.matrix_world)
# def drawpoly(x1, y1, x2, y2, r, g, b, a):
#     bgl.glLineWidth(1)
#     bgl.glColor4f(r, g, b, a)
#     bgl.glBegin(bgl.GL_POLYGON)
#     bgl.glVertex2i(x1, y2)
#     bgl.glVertex2i(x2, y2)
#     bgl.glVertex2i(x2, y1)
#     bgl.glVertex2i(x1, y1)
#     bgl.glEnd()
#     bgl.glColor4f(0.0, 0.0, 0.0, 1.0)


# def drawtri(posx, posy, l, d, hscale, radius):
#     r, g, b = colorsys.hsv_to_rgb(0.75 - l * 0.75, 1.0, 1.0)
#     a = 0.9
#     bgl.glEnable(bgl.GL_BLEND)
#     bgl.glBegin(bgl.GL_POLYGON)
#     bgl.glColor4f(r, g, b, a)
#     bgl.glVertex2f(posx - l * 0.5 * hscale * (radius - 20)*sin(d*pi/180), posy - l * 0.5 * hscale * (radius - 20)*cos(d*pi/180))
#     bgl.glVertex2f(posx + hscale * (l**0.5) * (radius/4 - 5)*cos(d*pi/180), posy - hscale * (l**0.5) * (radius/4 - 5)*sin(d*pi/180))
#     bgl.glVertex2f(posx + l**0.5 * hscale * (radius - 20)*sin(d*pi/180), posy + l**0.5 * hscale * (radius - 20)*cos(d*pi/180))
#     bgl.glVertex2f(posx - hscale * (l**0.5) * (radius/4 - 5)*cos(d*pi/180), posy + hscale * (l**0.5) * (radius/4 - 5)*sin(d*pi/180))
#     bgl.glEnd()
#     bgl.glDisable(bgl.GL_BLEND)


# def drawcircle(center, radius, resolution, fill, a, r, g, b):
#     bgl.glColor4f(r, g, b, a)
#     bgl.glEnable(bgl.GL_LINE_SMOOTH)
#     bgl.glEnable(bgl.GL_BLEND)
#     bgl.glBlendFunc(bgl.GL_SRC_ALPHA, bgl.GL_ONE_MINUS_SRC_ALPHA)
#     bgl.glHint(bgl.GL_LINE_SMOOTH_HINT, bgl.GL_NICEST)
#     bgl.glLineWidth(1.5)

#     if fill:
#         bgl.glBegin(bgl.GL_POLYGON)
#     else:
#         bgl.glBegin(bgl.GL_LINE_STRIP)

#     for i in range(resolution+1):
#         vec = Vector((cos(i/resolution*2*pi), sin(i/resolution*2*pi)))
#         v = vec * radius + center
#         bgl.glVertex2f(v.x, v.y)
#     bgl.glEnd()


# def drawbsdfcircle(centre, radius, resolution, fill, col, w, h, z, lw):
#     bgl.glEnable(bgl.GL_BLEND)
#     bgl.glEnable(bgl.GL_LINE_SMOOTH)
#     bgl.glBlendFunc(bgl.GL_SRC_ALPHA, bgl.GL_ONE_MINUS_SRC_ALPHA)
#     bgl.glHint(bgl.GL_LINE_SMOOTH_HINT, bgl.GL_NICEST)
#     bgl.glLineWidth(lw)
#     if not fill:
#         if col:
#             bgl.glColor4f(*col)

#         bgl.glBegin(bgl.GL_LINE_LOOP)
#     else:
#         bgl.glColor4f(*col)
#         bgl.glLineWidth(2.5)
#         bgl.glBegin(bgl.GL_POLYGON)
#     for p in range(0, resolution):
#         bgl.glVertex3f(centre[0] + radius * math.sin(math.pi * p/180) * w, centre[1] + radius * math.cos(math.pi * p/180) * h, z)

#     bgl.glEnd()
#     bgl.glDisable(bgl.GL_BLEND)


# def drawwedge(c, phis, rs, col, w, h):
#     bgl.glEnable(bgl.GL_BLEND)
#     bgl.glEnable(bgl.GL_LINE_SMOOTH)
#     bgl.glBlendFunc(bgl.GL_SRC_ALPHA, bgl.GL_ONE_MINUS_SRC_ALPHA)
#     bgl.glHint(bgl.GL_LINE_SMOOTH_HINT, bgl.GL_FASTEST)
#     (z, lw, col) = (0.1, 3, col) if col else (0.05, 1.5, [0, 0, 0, 0.25])
#     bgl.glColor4f(*col)
#     bgl.glLineWidth(lw)
#     bgl.glBegin(bgl.GL_LINE_LOOP)
#     for p in range(phis[0], phis[1] + 1):
#         bgl.glVertex3f(*radial2xy(c, rs[0], p, w, h), z)
#     for p in range(phis[1], phis[0] - 1, -1):
#         bgl.glVertex3f(*radial2xy(c, rs[1], p, w, h), z)
#     bgl.glLineWidth(1)

#     bgl.glEnd()
#     bgl.glDisable(bgl.GL_BLEND)
