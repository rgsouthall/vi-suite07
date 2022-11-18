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

import bpy, bmesh, os, mathutils
from .vi_func import selobj
from math import cos, sin, pi
from mathutils import Vector

ofheader = r'''/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM:    The Open Source CFD Toolbox        |
|  \\    /   O peration     | Version:     8                                  |
|   \\  /    A nd           | Web:         www.OpenFOAM.org                   |
|    \\/     M anipulation  | Created by:  FloVi (part of the VI-Suite)       |
\*---------------------------------------------------------------------------*/''' + '\n\n'


def fileheader(o):
    return '''FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      {};
}}
// * *

'''.format(o)


flovi_p_bounds = {'sf': {'0': ('zeroGradient', 'fixedValue', 'fixedMean', 'fixedMeanOutletInlet', 'freestreamPressure', 'totalPressure',
                               'inletOutlet', 'pressureInletOutletVelocity', 'PressureInletVelocity', 'calculated'),
                         '1': ('zeroGradient', 'calculated'), '2': ['None']},
                  'bsf': {'0': ('calculated'), '1': ('calculated'), '2': ('None',)}}

flovi_u_bounds = {'simpleFoam': {'0': ('zeroGradient', 'fixedValue', 'inletOutlet', 'freestream', 'pressureInletOutletVelocity', 'atmBoundaryLayerInletVelocity', 'slip'), '1': ('noSlip', 'fixedValue', 'slip'), '2': ['None']},
                'buoyantSimpleFoam': {'0': ('zeroGradient', 'fixedValue', 'inletOutlet', 'freestream', 'pressureInletOutletVelocity', 'atmBoundaryLayerInletVelocity', 'slip'), '1': ('noSlip', 'fixedValue', 'slip'), '2': ['None']}}

flovi_nut_bounds = {'simpleFoam': {'0': ['calculated'], '1': ['nutkWallFunction'], '2': ['None']},
                    'buoyantSimpleFoam': {'0': ['calculated'], '1': ['nutkWallFunction'], '2': ['None']}}

flovi_nutilda_bounds = {'simpleFoam': {'0': ('zeroGradient', 'fixedValue'), '1': ('zeroGradient', 'fixedValue'), '2': ['None']},
                        'buoyantSimpleFoam': {'0': ('zeroGradient', 'fixedValue'), '1': ('zeroGradient', 'fixedValue'), '2': ['None']}}

flovi_k_bounds = {'simpleFoam': {'0': ('fixedValue', 'inletOutlet'), '1': ['kqRWallFunction'], '2': ['None']},
                  'buoyantSimpleFoam': {'0': ('fixedValue', 'inletOutlet', 'turbulentIntensityKineticEnergyInlet'), '1': ['kqRWallFunction'], '2': ['None']}}

flovi_epsilon_bounds = {'simpleFoam': {'0': ('fixedValue', 'inletOutlet'), '1': ['epsilonWallFunction'], '2': ('None')},
                        'buoyantSimpleFoam': {'0': ('fixedValue', 'inletOutlet'), '1': ['epsilonWallFunction'], '2': ['None']}}

flovi_omega_bounds = {'simpleFoam': {'0': ('zeroGradient', 'fixedValue'), '1': ['omegaWallFunction'], '2': ['None']},
                      'buoyantSimpleFoam': {'0': ('zeroGradient', 'fixedValue'), '1': ['omegaWallFunction'], '2': ['None']}}

flovi_t_bounds = {'buoyantSimpleFoam': {'0': ('zeroGradient', 'fixedValue', 'inletOutlet'), '1': ('zeroGradient', 'fixedValue'), '2': ['None']}}

flovi_prgh_bounds = {'buoyantSimpleFoam': {'0': ('fixedFluxPressure', 'fixedValue', 'prghPressure', 'prghTotalPressure', 'prghEntrainmentPressure',
                                                 'fixedMean', 'fixedMeanOutletInlet'),
                                           '1': ('fixedFluxPressure', 'fixedValue', 'totalPressure'),
                                           '2': ['None']}}

flovi_a_bounds = {'buoyantSimpleFoam': {'0': ('calculated',),
                                        '1': ('compressible::alphatWallFunction', 'compressible::alphatJayatillekeWallFunction'),
                                        '2': ['None']}}
flovi_rad_bounds = {'buoyantSimpleFoam': {'0': ('MarshakRadiation',),
                                          '1': ('MarshakRadiation',),
                                          '2': ['None']}}


def ret_fvbp_menu(mat, context):
    return [('{}'.format(b), '{}'.format(b), '{} boundary type'.format(b)) for b in flovi_p_bounds[context.scene.vi_params['flparams']['solver_type']][mat.flovi_bmb_type]]


def ret_fvbu_menu(mat, context):
    return [('{}'.format(b), '{}'.format(b), '{} boundary type'.format(b)) for b in flovi_u_bounds[context.scene.vi_params['flparams']['solver']][mat.flovi_bmb_type]]


def ret_fvbnut_menu(mat, context):
    return [('{}'.format(b), '{}'.format(b), '{} boundary type'.format(b)) for b in flovi_nut_bounds[context.scene.vi_params['flparams']['solver']][mat.flovi_bmb_type]]


def ret_fvbnutilda_menu(mat, context):
    return [('{}'.format(b), '{}'.format(b), '{} boundary type'.format(b)) for b in flovi_nutilda_bounds[context.scene.vi_params['flparams']['solver']][mat.flovi_bmb_type]]


def ret_fvbk_menu(mat, context):
    return [('{}'.format(b), '{}'.format(b), '{} boundary type'.format(b)) for b in flovi_k_bounds[context.scene.vi_params['flparams']['solver']][mat.flovi_bmb_type]]


def ret_fvbepsilon_menu(mat, context):
    return [('{}'.format(b), '{}'.format(b), '{} boundary type'.format(b)) for b in flovi_epsilon_bounds[context.scene.vi_params['flparams']['solver']][mat.flovi_bmb_type]]


def ret_fvbomega_menu(mat, context):
    return [('{}'.format(b), '{}'.format(b), '{} boundary type'.format(b)) for b in flovi_omega_bounds[context.scene.vi_params['flparams']['solver']][mat.flovi_bmb_type]]


def ret_fvbt_menu(mat, context):
    return [('{}'.format(b), '{}'.format(b), '{} boundary type'.format(b)) for b in flovi_t_bounds[context.scene.vi_params['flparams']['solver']][mat.flovi_bmb_type]]

def ret_fvba_menu(mat, context):
    return [('{}'.format(b), '{}'.format(b), '{} boundary type'.format(b)) for b in flovi_a_bounds[context.scene.vi_params['flparams']['solver']][mat.flovi_bmb_type]]


def ret_fvbprgh_menu(mat, context):
    return [('{}'.format(b), '{}'.format(b), '{} boundary type'.format(b)) for b in flovi_prgh_bounds[context.scene.vi_params['flparams']['solver']][mat.flovi_bmb_type]]

def ret_fvrad_menu(mat, context):
    return [('{}'.format(b), '{}'.format(b), '{} boundary type'.format(b)) for b in flovi_rad_bounds[context.scene.vi_params['flparams']['solver']][mat.flovi_bmb_type]]


def write_header(func):
    def wrapper(o, expnode):
        return ofheader + func(o, expnode)
    return wrapper


def fventry(func):
    return '    {\n'


def write_ffile(cla, loc, obj):
    location = 'location    "{}";\n    '.format(loc) if loc else ''
    return 'FoamFile\n  {{\n    version   2.0;\n    format    ascii;\n    {}class     {};\n    object    {};\n  }}\n\n'.format(location, cla, obj)


def write_fvdict(text, fvdict):
    for d in fvdict:
        if isinstance(fvdict[d], str):
            text += '{} {};\n\n'.format(d, fvdict[d])
        elif isinstance(fvdict[d], dict):
            text += '{}\n{{\n'.format(d)
            for sd in fvdict[d]:
                if isinstance(fvdict[d][sd], str):
                    text += '  {} {};\n'.format(sd, fvdict[d][sd])
                elif isinstance(fvdict[d][sd], dict):
                    text += '  {}\n  {{\n'.format(sd)
                    for ssd in fvdict[d][sd]:
                        if isinstance(fvdict[d][sd][ssd], str):
                            text += '    {} {};\n'.format(ssd, fvdict[d][sd][ssd])
                        elif isinstance(fvdict[d][sd][ssd], dict):
                            text += '  {}\n  {{\n'.format(ssd)
                            for sssd in fvdict[d][sd][ssd]:
                                if isinstance(fvdict[d][sd][ssd][sssd], str):
                                    text += '    {} {};\n'.format(sssd, fvdict[d][sd][ssd][sssd])
                            text += '  }\n'
                    text += '  }\n'
            text += '}\n'
    return text


def fvboundwrite(o):
    boundary = ''

    for mat in o.data.materials:
        boundary += "  {}\n  {{\n    type {};\n    faces\n    (\n".format(mat.name, ("patch", "wall", "symmetry", "empty")[int(mat.flovi_bmb_type)])#;\n\n"
        faces = [face for face in o.data.polygons if o.data.materials[face.material_index] == mat]

        for face in faces:
            boundary += "      ("+" ".join([str(v) for v in face.vertices])+")\n"

        boundary += "    );\n  }\n"
    boundary += ");\n\nmergePatchPairs\n(\n);"
    return boundary


def write_bound(o, m, ns, nf):
    t = ("patch", "wall", "symmetry", "empty")[int(m.vi_params.flovi_bmb_type)]

    return '''   {0}_{1}
    {{
        type            {2};
        nFaces          {3};
        startFace       {4};
    }}\n'''.format(o.name, m.name, t, nf, ns)


# def flovi_bm_update(self, context):
#     scene = context.scene
#     svp = scene.vi_params
#     ovp = context.object.vi_params
#     svp = scene.vi_params
#    svp['flparams']['solver'] = ovp.flovi_solver
#    svp['flparams']['turbulence'] = ovp.flovi_turb


# # @writeheader
# def fvbmwrite(o, expnode):
#     bm = bmesh.new()
#     tempmesh = o.to_mesh(scene=bpy.context.scene, apply_modifiers=True, settings='PREVIEW')
#     bm.from_mesh(tempmesh)
#     bm.verts.ensure_lookup_table()
#     bm.transform(o.matrix_world)
#     bpy.data.meshes.remove(tempmesh)
#     [xs, ys, zs] = [[v.co[i] for v in bm.verts] for i in range(3)]
#     bm.transform(mathutils.Matrix.Translation((-min(xs), -min(ys), -min(zs))))
#     o['flovi_translate'] = (-min(xs), -min(ys), -min(zs))
#     lengths = [mathutils.Vector(v.co).length for v in bm.verts]
#     vert0 = bm.verts[lengths.index(min(lengths))]
#     angles = [mathutils.Vector(v.co).angle(mathutils.Vector((0, 0, 1))) for v in bm.verts if v != vert0]
# #    vert0 = [v for v in bm.verts if v.co[:] == (min(xs), min(ys), min(zs))][0]
#     vert4 = bm.verts[angles.index(min(angles)) + 1]
# #    print(vert0.index, vert4.index)
# #    vert4 = [v for v in bm.verts if (v.co[0], v.co[1]) == (vert0.co[0], vert0.co[1]) and v.co[2] != vert0.co[2]][0]

#     for face in bm.faces:
#         if vert0 in face.verts and vert4 not in face.verts:
#             vis = [vert.index for vert in face.verts][::-1]
#             vertorder1 = vis[vis.index(vert0.index):] + vis[:vis.index(vert0.index)]
# #            vertorder1 = [vertorder1[0], vertorder1[3], vertorder1[2], vertorder1[1]]
#         if vert4 in face.verts and vert0 not in face.verts:
#             vis = [vert.index for vert in face.verts]
#             vertorder2 = vis[vis.index(vert4.index):] + vis[:vis.index(vert4.index)]

#     vertorder = ''.join(['{} '.format(v) for v in vertorder1 + vertorder2])

# #    omw, bmovs = o.matrix_world, [vert for vert in o.data.vertices]
# #    xvec, yvec, zvec = (omw*bmovs[3].co - omw*bmovs[0].co).normalized(), (omw*bmovs[2].co - omw*bmovs[3].co).normalized(), (omw*bmovs[4].co - omw*bmovs[0].co).normalized()
# #    ofvpos = [[(omw*bmov.co - omw*bmovs[0].co)*vec for vec in (xvec, yvec, zvec)] for bmov in bmovs]
# #    bmdict = "vertices\n(\n" + "\n".join(["  ({0:.3f} {1:.3f} {2:.3f})" .format(*ofvpo) for ofvpo in ofvpos]) +"\n);\n\n"
#     bmdict = "vertices\n(\n" + "\n".join(["  ({0[0]:.3f} {0[1]:.3f} {0[2]:.3f})" .format(v.co) for v in bm.verts]) +"\n);\n\n"
#     bmdict += "blocks\n(\n  hex ({}) ({} {} {}) simpleGrading ({} {} {})\n);\n\n".format(vertorder, expnode.bm_xres, expnode.bm_yres, expnode.bm_zres, expnode.bm_xgrad, expnode.bm_ygrad, expnode.bm_zgrad)
#     bmdict += "edges\n(\n);\n\nboundary\n(\n"
#     bmdict += fvboundwrite(o)
#     bm.free()
#     return ofheader + write_ffile('dictionary', '', 'blockMeshDict') + bmdict


# def fvblbmgen(mats, ffile, vfile, bfile, meshtype):
#     scene = bpy.context.scene
#     matfacedict = {mat.name: [0, 0] for mat in mats}

#     for line in bfile.readlines():
#         if line.strip() in matfacedict:
#             mat = line.strip()
#         elif line.strip() in [o.name for o in bpy.data.objects]:
#             mat = bpy.data.objects[line.strip()].data.materials[0].name
#         if 'nFaces' in line:
#             matfacedict[mat][1] = int(line.split()[1].strip(';'))
#         if 'startFace' in line:
#             matfacedict[mat][0] = int(line.split()[1].strip(';'))
#     bobs = [ob for ob in scene.objects if ob.get('VIType') and ob['VIType'] == 'FloViMesh']

#     if bobs:
#         o = bobs[0]
#         selobj(scene, o)
#         while o.data.materials:
#             bpy.ops.object.material_slot_remove()
#     else:
#         bpy.ops.object.add(type='MESH', layers=(False, False, False, False, False, False,
#                                                 False, False, False, False, False, False,
#                                                 False, False, False, False, False, False, False, True))
#         o = bpy.context.object
#         o['VIType'] = 'FloViMesh'

#     o.name = meshtype
#     for mat in mats:
#         if mat.name not in o.data.materials:
#             bpy.ops.object.material_slot_add()
#             o.material_slots[-1].material = mat

#     matnamedict = {mat.name: m for m, mat in enumerate(o.data.materials)}
#     bm = bmesh.new()

#     for line in [line for line in vfile.readlines() if line[0] == '(' and len(line.split(' ')) == 3]:
#         bm.verts.new([float(vpos) for vpos in line[1:-2].split(' ')])

#     if hasattr(bm.verts, "ensure_lookup_table"):
#         bm.verts.ensure_lookup_table()

#     for l, line in enumerate([line for line in ffile.readlines() if '(' in line and line[0].isdigit() and len(line.split(' ')) == int(line[0])]):
#         newf = bm.faces.new([bm.verts[int(fv)] for fv in line[2:-2].split(' ')])

#         for facerange in matfacedict.items():
#             if l in range(facerange[1][0], facerange[1][0] + facerange[1][1]):
#                 newf.material_index = matnamedict[facerange[0]]

#     bm.transform(o.matrix_world.inverted())
#     bm.to_mesh(o.data)
#     bm.free()


def fvmat(self, svp, mn, bound, frame):
    begin = '\n  {}\n  {{\n    type    '.format(mn)
    end = ';\n  }\n'

    if bound == 'p':
        val = 'uniform {}'.format(self.flovi_bmbp_val) if not self.flovi_p_field else '$internalField'
        pdict = {'0': self.flovi_bmbp_subtype, '1': self.flovi_bmbp_subtype, '2': 'symmetry', '3': 'empty'}
        ptdict = {'zeroGradient': 'zeroGradient',
                  'fixedValue': 'fixedValue;\n    value    {}'.format(val),
                  'calculated': 'calculated;\n    value    $internalField',
                  'freestreamPressure': 'freestreamPressure;\n  freestreamValue {}'.format(val),
                  'fixedMeanOutletInlet': 'fixedMeanOutletInlet;\n    meanValue   {0};\n    value {0}'.format(val),
                  'fixedMean': 'fixedMean;\n    meanValue   {0};\n    value {0}'.format(val),
                  'totalPressure': 'totalPressure;\n    p0      uniform {};\n    gamma    {};\n    value    {}'.format(self.flovi_bmbp_p0val, self.flovi_bmbp_gamma, val),
                  'symmetry': 'symmetry', 'empty': 'empty'}
        entry = ptdict[pdict[self.flovi_bmb_type]]

    elif bound == 'U':
        if self.flovi_u_type == '0':
            v_vec = self.flovi_bmbu_val
            val = 'uniform ({:.4f} {:.4f} {:.4f})'.format(*v_vec) if not self.flovi_u_field else '$internalField'
        else:
            v_vec = Vector((self.flovi_u_speed * sin(pi*self.flovi_u_azi/180), self.flovi_u_speed * cos(pi*self.flovi_u_azi/180), 0))
            val = 'uniform ({:.4f} {:.4f} {:.4f})'.format(*v_vec) if not self.flovi_u_field else '$internalField'

        if self.flovi_bmbu_subtype == 'atmBoundaryLayerInletVelocity':
            fd_val = self.flovi_u_fdir if not self.flovi_u_field else svp['flparams'][str(frame)]['Udir']
            s_val = self.flovi_u_ref if not self.flovi_u_field else svp['flparams'][str(frame)]['Uspeed']
        else:
            fd_val, s_val = (0, 0, 0), 0

        Udict = {'0': self.flovi_bmbu_subtype, '1': self.flovi_bmbu_subtype, '2': 'symmetry', '3': 'empty'}
        Utdict = {'fixedValue': 'fixedValue;\n    value    {}'.format(val), 'slip': 'slip', 'noSlip': 'noSlip',
                  'inletOutlet': 'inletOutlet;\n    inletValue    $internalField;\n    value    {}'.format(val),
                  'pressureInletOutletVelocity': 'pressureInletOutletVelocity;\n    value    {}'.format(val),
                  'zeroGradient': 'zeroGradient', 'symmetry': 'symmetry',
                  'freestream': 'freestream;\n    freestreamValue    $internalField',
                  'calculated': 'calculated;\n    value    $internalField',
                  'atmBoundaryLayerInletVelocity': 'atmBoundaryLayerInletVelocity;\n    Uref {0:.3f};\n    Zref {1:.3f};\n    zDir ({2[0]:.3f} {2[1]:.3f} {2[2]:.3f});\n    flowDir    ({3[0]:.3f} {3[1]:.3f} {3[2]:.3f});\n    z0 uniform {4:.3f};\n    zGround uniform {5:.3f};\n    d uniform {6:.3f}\n    value {7}'.format(s_val,
                                                                                                                                                                                                                                 self.flovi_u_zref,
                                                                                                                                                                                                                                 self.flovi_u_zdir,
                                                                                                                                                                                                                                 fd_val,
                                                                                                                                                                                                                                 self.flovi_u_z0,
                                                                                                                                                                                                                                 self.flovi_u_zground,
                                                                                                                                                                                                                                 self.flovi_u_d,
                                                                                                                                                                                                                                 val),
                  'empty': 'empty'}

        entry = Utdict[Udict[self.flovi_bmb_type]]

    elif bound == 'nut':
        ndict = {'0': self.flovi_bmbnut_subtype, '1': self.flovi_bmbnut_subtype, '2': 'symmetry', '3': 'empty'}
        ntdict = {'nutkWallFunction': 'nutkWallFunction;\n    value    $internalField',
                  'nutUSpaldingWallFunction': 'nutUSpaldingWallFunction;\n    value    $internalField',
                  'calculated': 'calculated;\n    value    $internalField',
                  'inletOutlet': 'inletOutlet;\n    inletValue    $internalField\n    value    $internalField',
                  'symmetry': 'symmetry', 'empty': 'empty'}
        entry = ntdict[ndict[self.flovi_bmb_type]]

    elif bound == 'k':
        val = '{:.4f}'.format(self.flovi_k_val) if not self.flovi_k_field else '$internalField'
        ival = '{:.4f}'.format(self.flovi_k_intensity) if not self.flovi_k_field else '$internalField'
        kdict = {'0': self.flovi_k_subtype, '1': self.flovi_k_subtype, '2': 'symmetry', '3': 'empty'}
        ktdict = {'fixedValue': 'fixedValue;\n    value    $internalField',
                  'kqRWallFunction': 'kqRWallFunction;\n    value    $internalField',
                  'inletOutlet': 'inletOutlet;\n    inletValue    $internalField;\n    value    $internalField',
                  'calculated': 'calculated;\n    value    $internalField',
                  'symmetry': 'symmetry',
                  'empty': 'empty',
                  'turbulentIntensityKineticEnergyInlet': 'turbulentIntensityKineticEnergyInlet;\n    intensity       {};\n    value      {}'.format(ival, val)}
        entry = ktdict[kdict[self.flovi_bmb_type]]

    elif bound == 't':
        val = 'uniform {}'.format(self.flovi_bmbt_val) if not self.flovi_t_field else '$internalField'
        ival = 'uniform {}'.format(self.flovi_bmbti_val) if not self.flovi_t_field else '$internalField'
        tdict = {'0': self.flovi_bmbt_subtype, '1': self.flovi_bmbt_subtype, '2': 'symmetry', '3': 'empty'}
        ttdict = {'zeroGradient': 'zeroGradient',
                  'fixedValue': 'fixedValue;\n    value    {}'.format(val),
                  'inletOutlet': 'inletOutlet;\n    inletValue    {}\n    value    {}'.format(ival, val),
                  'calculated': 'calculated;\n    value    $internalField',
                  'symmetry': 'symmetry',
                  'empty': 'empty'}
        entry = ttdict[tdict[self.flovi_bmb_type]]

    elif bound == 'p_rgh':
        val = 'uniform {:.4f}'.format(self.flovi_prgh_val) if not self.flovi_prgh_field else '$internalField'
        p0val = 'uniform {:.4f}'.format(self.flovi_prgh_p0) if not self.flovi_prgh_field else '$internalField'
        prghdict = {'0': self.flovi_prgh_subtype, '1': self.flovi_prgh_subtype, '2': 'symmetry', '3': 'empty'}
        prghtdict = {'fixedFluxPressure': 'fixedFluxPressure;\n    value    {}'.format(val),
                     'fixedValue': 'fixedValue;\n    value    {}'.format(val),
                     'prghEntrainmentPressure': 'prghEntrainmentPressure;\n    p0       $internalField'.format(p0val),
                     'calculated': 'calculated;\n    value    $internalField',
                     'prghPressure': 'prghPressure;\n    p    $internalField;\n    value    $internalField',
                     'freestreamPressure': 'freestreamPressure',
                     'prghTotalPressure': 'prghTotalPressure;\n    p0      {};\n    value    {}'.format(p0val, val),
                     'fixedMeanOutletInlet': 'fixedMeanOutletInlet;\n    meanValue   {0};\n    value {0}'.format(val),
                     'fixedMean': 'fixedMean;\n    meanValue   {0};\n    value {0}'.format(val),
                     'symmetry': 'symmetry', 'empty': 'empty'}
        entry = prghtdict[prghdict[self.flovi_bmb_type]]

    elif bound == 'a':
        val = 'uniform {}'.format(self.flovi_a_val) if not self.flovi_a_field else '$internalField'
        tdict = {'0': self.flovi_a_subtype, '1': self.flovi_a_subtype, '2': 'symmetry', '3': 'empty'}
        ttdict = {'compressible::alphatWallFunction': 'compressible::alphatWallFunction;\n    value     $internalField',
                  'compressible::alphatJayatillekeWallFunction': 'compressible::alphatJayatillekeWallFunction;\n    Prt    0.85;\n    value     $internalField',
                  'fixedValue': 'fixedValue;\n    value    {}'.format(val),
                  'inletOutlet': 'inletOutlet;\n    inletValue    $internalField\n    value    $internalField',
                  'calculated': 'calculated;\n    value    $internalField', 'symmetry': 'symmetry', 'empty': 'empty'}
        entry = ttdict[tdict[self.flovi_bmb_type]]

    elif bound == 'e':
        edict = {'0': self.flovi_bmbe_subtype, '1': self.flovi_bmbe_subtype, '2': 'symmetry', '3': 'empty'}
        etdict = {'symmetry': 'symmetry',
                  'empty': 'empty',
                  'inletOutlet': 'inletOutlet;\n    inletValue    $internalField;\n    value    $internalField',
                  'fixedValue': 'fixedValue;\n    value    $internalField',
                  'epsilonWallFunction': 'epsilonWallFunction;\n    value    $internalField',
                  'calculated': 'calculated;\n    value    $internalField'}
        entry = etdict[edict[self.flovi_bmb_type]]

    elif bound == 'o':
        odict = {'0': self.flovi_bmbo_subtype, '1': self.flovi_bmbo_subtype, '2': 'symmetry', '3': 'empty'}
        otdict = {'symmetry': 'symmetry',
                  'empty': 'empty',
                  'inletOutlet': 'inletOutlet;\n    inletValue    $internalField\n    value    $internalField',
                  'zeroGradient': 'zeroGradient',
                  'omegaWallFunction': 'omegaWallFunction;\n    value    $internalField',
                  'fixedValue': 'fixedValue;\n    value    $internalField'}
        entry = otdict[odict[self.flovi_bmb_type]]

    elif bound == 'nutilda':
        ntdict = {'0': self.flovi_bmbnutilda_subtype, '1': self.flovi_bmbnutilda_subtype, '2': 'symmetry', '3': 'empty'}
        nttdict = {'fixedValue': 'fixedValue;\n    value    $internalField', 'inletOutlet': 'inletOutlet;\n    inletValue    $internalField\n    value    $internalField', 'empty': 'empty',
                   'zeroGradient': 'zeroGradient', 'freestream': 'freestream\n    freeStreamValue  $internalField\n', 'symmetry': 'symmetry'}
        entry = nttdict[ntdict[self.flovi_bmb_type]]

    elif bound == 'G':
        raddict = {'0': self.flovi_rad_subtype, '1': self.flovi_rad_subtype, '2': 'symmetry', '3': 'empty'}
        radtdict = {'MarshakRadiation': 'MarshakRadiation;\n    emissivityMode    {};\n    emissivity    uniform {};\n    value    uniform {}'.format(self.flovi_rad_em, self.flovi_rad_e, self.flovi_rad_val), 'symmetry': 'symmetry'}
        entry = radtdict[raddict[self.flovi_bmb_type]]

    return begin + entry + end


def fvvarwrite(scene, obs, node):
    '''Turbulence modelling: k and epsilon required for kEpsilon, k and omega required for kOmega, nutilda required for SpalartAllmaras, nut required for all
        Buoyancy modelling: T'''
    svp = scene.vi_params

    for frame in range(svp['flparams']['start_frame'], svp['flparams']['end_frame'] + 1):
        scene.frame_set(frame)
        frame_of0fb = os.path.join(svp['flparams']['offilebase'], str(frame), '0')

        if not os.path.isdir(frame_of0fb):
            os.makedirs(frame_of0fb)

        uval = node.uval if node.uval_type == '0' else Vector((sin(node.uval_azi * pi/180) * node.umag, cos(node.uval_azi * pi/180) * node.umag, 0))
        svp['flparams'][str(frame)] = {}
        svp['flparams'][str(frame)]['Udir'] = uval.normalized()
        svp['flparams'][str(frame)]['Uspeed'] = uval.length

        if not node.buoyancy:  # or (node.buoyancy and not node.buossinesq):
            pentry = "dimensions [{} {} {} {} 0 0 0];\ninternalField   uniform {};\n\nboundaryField\n{{\n".format('0', '2', '-2', '0',
                                                                                                                '{}'.format(node.pnormval))
        else:
            pentry = "dimensions [{} {} {} {} 0 0 0];\ninternalField   uniform {};\n\nboundaryField\n{{\n".format('1', '-1', '-2', '0',
                                                                                                                '{}'.format(node.pabsval))

        (Uentry, nutildaentry, nutentry, kentry, eentry, oentry, tentry, p_rghentry, aentry, Gentry) = ["dimensions [{} {} {} {} 0 0 0];\ninternalField   uniform {};\n\nboundaryField\n{{\n".format(*var) for var in (
                                                                                    ('0', '1', '-1', '0', '({:.4f} {:.4f} {:.4f})'.format(*uval)),
                                                                                    ('0', '2', '-1', '0', '{:.4f}'.format(node.nutildaval)),
                                                                                    ('0', '2', '-1', '0', '{:.4f}'.format(node.nutval)),
                                                                                    ('0', '2', '-2', '0', '{:.4f}'.format(node.kval)),
                                                                                    ('0', '2', '-3', '0', '{:.4f}'.format(node.epval)),
                                                                                    ('0', '0', '-1', '0', '{:.4f}'.format(node.oval)),
                                                                                    ('0', '0', '0', '1', '{:.4f}'.format(node.tval)),
                                                                                    ('1', '-1', '-2', '0', '{:.4f}'.format(node.p_rghval)),
                                                                                    ('1', '-1', '-1', '0', '{:.4f}'.format(node.aval)),
                                                                                    ('1', '0', '-3', '0', '{:.4f}'.format(node.Gval)))]
        for o in obs:
            for mat in o.data.materials:
                mvp = mat.vi_params
                matname = '{}_{}'.format(o.name, mat.name)

                if mvp.mattype == '2':
                    pentry += mvp.flovi_mat(svp, matname, 'p', frame)
                    Uentry += mvp.flovi_mat(svp, matname, 'U', frame)

    #                if node.turbulence != 'laminar':
                    if node.turbulence != 'laminar':
                        nutentry += mvp.flovi_mat(svp, matname, 'nut', frame)
                        if node.turbulence == 'SpalartAllmaras':
                            nutildaentry += mvp.flovi_mat(svp, matname, 'nutilda', frame)
                        elif node.turbulence == 'kEpsilon':
                            kentry += mvp.flovi_mat(svp, matname, 'k', frame)
                            eentry += mvp.flovi_mat(svp, matname, 'e', frame)
                        elif node.turbulence == 'kOmega':
                            kentry += mvp.flovi_mat(svp, matname, 'k', frame)
                            oentry += mvp.flovi_mat(svp, matname, 'o', frame)
                    if node.buoyancy:
                        tentry += mvp.flovi_mat(svp, matname, 't', frame)
                        p_rghentry += mvp.flovi_mat(svp, matname, 'p_rgh', frame)
                        aentry += mvp.flovi_mat(svp, matname, 'a', frame)

                        if node.radiation:
                            Gentry += mvp.flovi_mat(svp, matname, 'G', frame)

        pentry += '}'
        Uentry += '}'
        nutentry += '}'
        nutildaentry += '}'
        kentry += '}'
        eentry += '}'
        oentry += '}'
        tentry += '}'
        p_rghentry += '}'
        aentry += '}'
        Gentry += '}'

        with open(os.path.join(frame_of0fb, 'p'), 'w') as pfile:
            pfile.write(ofheader + write_ffile('volScalarField', '', 'p') + pentry)
        with open(os.path.join(frame_of0fb, 'U'), 'w') as Ufile:
            Ufile.write(ofheader + write_ffile('volVectorField', '', 'U') + Uentry)

        if node.turbulence != 'laminar':
            with open(os.path.join(frame_of0fb, 'nut'), 'w') as nutfile:
                nutfile.write(ofheader + write_ffile('volScalarField', '', 'nut') + nutentry)

            if node.turbulence == 'SpalartAllmaras':
                with open(os.path.join(frame_of0fb, 'nuTilda'), 'w') as nutildafile:
                    nutildafile.write(ofheader + write_ffile('volScalarField', '', 'nut') + nutildaentry)

            if node.turbulence == 'kEpsilon':
                with open(os.path.join(frame_of0fb, 'k'), 'w') as kfile:
                    kfile.write(ofheader + write_ffile('volScalarField', '', 'k') + kentry)
                with open(os.path.join(frame_of0fb, 'epsilon'), 'w') as efile:
                    efile.write(ofheader + write_ffile('volScalarField', '', 'epsilon') + eentry)

            if node.turbulence == 'kOmega':
                with open(os.path.join(frame_of0fb, 'k'), 'w') as kfile:
                    kfile.write(ofheader + write_ffile('volScalarField', '', 'k') + kentry)
                with open(os.path.join(frame_of0fb, 'omega'), 'w') as ofile:
                    ofile.write(ofheader + write_ffile('volScalarField', '', 'omega') + oentry)

        if node.buoyancy:
            with open(os.path.join(frame_of0fb, 'T'), 'w') as tfile:
                tfile.write(ofheader + write_ffile('volScalarField', '', 'T') + tentry)
            with open(os.path.join(frame_of0fb, 'alphat'), 'w') as afile:
                afile.write(ofheader + write_ffile('volScalarField', '', 'alphat') + aentry)
            with open(os.path.join(frame_of0fb, 'p_rgh'), 'w') as prghfile:
                prghfile.write(ofheader + write_ffile('volScalarField', '', 'p_rgh') + p_rghentry)

            if node.radiation:
                with open(os.path.join(frame_of0fb, 'G'), 'w') as Gfile:
                    Gfile.write(ofheader + write_ffile('volScalarField', '', 'G') + Gentry)


def fvmattype(mat, var):
    if mat.flovi_bmb_type == '0':
        matbptype = ['zeroGradient'][int(mat.flovi_bmwp_type)]
        matbUtype = ['fixedValue'][int(mat.flovi_bmwu_type)]
    elif mat.flovi_bmb_type in ('1', '2'):
        matbptype = ['freestreamPressure'][int(mat.flovi_bmiop_type)]
        matbUtype = ['fixedValue'][int(mat.flovi_bmiou_type)]
    elif mat.flovi_bmb_type == '3':
        matbptype = 'empty'
        matbUtype = 'empty'


def fvcdwrite(svp, dp, solver, st, dt, et):
    pw = 0 if solver == 'icoFoam' else 0
    ps, ss, bs = [], [], []
    htext = ofheader + write_ffile('dictionary', 'system', 'controlDict')
    cdict = {'application': solver, 'startFrom': 'startTime', 'startTime': '{}'.format(st), 'stopAt': 'endTime',
             'endTime': '{}'.format(et), 'deltaT': '{:.5f}'.format(dt), 'writeControl': 'timeStep', 'writeInterval': '10',
             'purgeWrite': '{}'.format(pw), 'writeFormat': 'ascii', 'writePrecision': '6', 'writeCompression': 'off',
             'timeFormat': 'general', 'timePrecision': '6', 'runTimeModifiable': 'true', 'functions': {}, 'libs': '("libatmosphericModels.so")'}

    for o in bpy.data.objects:
        ovp = o.vi_params
        if o.type == 'MESH' and ovp.vi_type == '2':
            dom = o
        elif o.type == 'EMPTY' and ovp.flovi_probe:
            ps.append(o)
        elif o.type == 'MESH' and ovp.vi_type == '6':
            for frame in range(svp['flparams']['start_frame'], svp['flparams']['end_frame'] + 1):
                if not os.path.isdir(os.path.join(svp['flparams']['offilebase'], str(frame), 'constant', 'triSurface')):
                    os.makedirs(os.path.join(svp['flparams']['offilebase'], str(frame), 'constant', 'triSurface'))
                    ovp.write_stl(dp, os.path.join(svp['flparams']['offilebase'], str(frame), 'constant', 'triSurface', '{}.stl'.format(o.name)))
            ss.append(o.name)

        if o.type == 'MESH' and ovp.vi_type in ('2', '3') and any([m.vi_params.flovi_probe for m in o.data.materials]):
            for mat in o.data.materials:
                if mat.vi_params.flovi_probe:
                    bs.append('{}_{}'.format(o.name, mat.name))

    if ps:
        svp['flparams']['probes'] = [p.name for p in ps]
        probe_vars = 'p U T'

        for p in ps:
            cdict['functions'][p.name] = {'libs': '("libsampling.so")', 'type': 'probes', 'name': '{}'.format(p.name), 'writeControl': 'timeStep',
                                          'writeInterval': '1', 'fields': '({0})'.format(probe_vars),
                                          'probeLocations\n(\n{}\n)'.format('   ({0[0]} {0[1]} {0[2]})'.format(p.location)): ''}

        probe_text = '''functions
{{
    probes
    {{
        libs            ("libsampling.so");
        type            probes;
        name            {2};
        writeControl    timeStep;
        writeInterval   1;
        fields          ({0});
        probeLocations
        (
            {1}
        );
    }}
}}'''.format(probe_vars, ''.join([' ({0[0]} {0[1]} {0[2]})\n'.format(p.location) for p in ps]), ','.join(['{}'.format(p.name) for p in ps]))

    else:
        probe_text = ''
        bpy.context.scene.vi_params['flparams']['probes'] = []

    bpy.context.scene.vi_params['flparams']['s_probes'] = ss
    bpy.context.scene.vi_params['flparams']['b_probes'] = bs
        # for o in ss:
        #     cdict['functions'][o.name] = {'name': o.name, 'type': 'surfaceFieldValue', 'libs': '("libfieldFunctionObjects.so")',
        #         'writeControl': 'timeStep', 'writeInterval': '1', 'writeFields': 'false', 'log': 'true', 'operation': 'average',
        #         'fields': '(U)', 'operation': 'areaNormalIntegrate', 'regionType': 'sampledSurface', 'surfaceFormat': 'stl',
        #         'sampledSurfaceDict': {'type': 'triSurfaceMesh', 'surface': '{}.stl'.format(o.name), 'source': 'cells', 'interpolate': 'true'}}

    if bs:
        for b in bs:
            cdict['functions'][b] = {'type': 'surfaceFieldValue', 'libs': '("libfieldFunctionObjects.so")', 'writeControl': 'writeTime',
                                     'writeFields': 'true', 'surfaceFormat': 'raw', 'regionType': 'patch', 'name': '{}'.format(b),
                                     'operation': 'none', 'fields    (      p    )': ''}
    return write_fvdict(htext, cdict)
    return 'FoamFile\n{\n  version     2.0;\n  format      ascii;\n  class       dictionary;\n  location    "system";\n  object      controlDict;\n}\n\n' + \
           'application     {};\nstartFrom       startTime;\nstartTime       {};\nstopAt          endTime;\nendTime         {:.5f};\n'.format(solver, st, et) +\
           'deltaT          {:.5f};\nwriteControl    timeStep;\nwriteInterval   {};\npurgeWrite      {};\nwriteFormat     ascii;\nwritePrecision  6;\n'.format(dt, 1, pw) +\
           'writeCompression off;\ntimeFormat      general;\ntimePrecision   6;\nrunTimeModifiable true;\n\n' + probe_text


def fvprefwrite(node):
    htext = ofheader + write_ffile('uniformDimensionedScalarField', 'constant', 'pRef')
    pdict = {'dimensions': '[1 -1 -2 0 0 0 0]', 'value': '{}'.format(node.pabsval)}
    text = write_fvdict(htext, pdict)
    return text


def fvsolwrite(node, features):
    # turb = 'k' in residuals or 'nuTilda' in residuals
    # buoy = 'p_rgh' in residuals
    # print(residuals, 'p_rgh' in residuals)
    # buoss = 'e' in residuals
    # rad = 'G' in residuals
    soldict = {}
    soldict['SIMPLE'] = {}
    soldict['SIMPLE']['residualControl'] = {}
    soldict['SIMPLE']['residualControl']['U'] = '{:.5f}'.format(node.uresid)
    soldict['SIMPLE']['nNonOrthogonalCorrectors'] = '0'

    if node.p_ref != '0':
        soldict['SIMPLE']['pRefValue'] = '{}'.format(node.p_ref_val)

        if node.p_ref:
            soldict['SIMPLE']['pRefPoint'] = '({0[0]} {0[1]} {0[2]})'.format(bpy.data.objects[node.p_ref_point].location)
        else:
            soldict['SIMPLE']['pRefCell'] = '0'

    soldict['relaxationFactors' ] = {}
    soldict['relaxationFactors']['equations'] = {}
    soldict['relaxationFactors']['equations']['U'] = '0.7'
    soldict['relaxationFactors']['fields'] = {}
    soldict['solvers'] = {}
    soldict['solvers']['p'] = {'solver': 'GAMG', 'smoother': 'GaussSeidel', 'tolerance': '1e-6', 'relTol': '0.1'}
    soldict['solvers']['"U|k|epsilon|omega|nuTilda"'] = {'solver': 'smoothSolver', 'smoother': 'symGaussSeidel', 'tolerance': '1e-6', 'relTol': '0.1'}

#     if not features['turb']:
#         soldict['solvers']['U'] = {'solver': 'GAMG', 'smoother': 'GaussSeidel', 'tolerance': '1e-6', 'relTol': '0.1'}
# # #        soldict['SIMPLE']['nNonOrthogonalCorrectors'] = '0'
# #         # soldict['solvers']['p'] = {'solver': 'PCG', 'preconditioner': 'DIC', 'tolerance': '1e-6', 'relTol': '0.05'}


    if features['turb']:
    #    soldict['solvers']['"U|k|epsilon|omega|nuTilda"'] = {'solver': 'smoothSolver', 'smoother': 'symGaussSeidel', 'tolerance': '1e-6', 'relTol': '0.1'}
        soldict['SIMPLE']['residualControl']['"(k|epsilon|omega|nuTilda)"'] = '{:.5f}'.format(node.keoresid)
        soldict['relaxationFactors']['equations']['"(k|epsilon|omega|nuTilda)"'] = '0.7'

    if not features['buoy']:
        soldict['SIMPLE']['residualControl']['p'] = '{:.5f}'.format(node.presid)
        soldict['solvers']['U'] = {'solver': 'smoothSolver', 'smoother': 'symGaussSeidel', 'tolerance': '1e-6', 'relTol': '0.05'}
#        soldict['solvers']['p'] = {'solver': 'PCG', 'preconditioner': 'DIC', 'tolerance': '1e-6', 'relTol': '0.05'}
        soldict['relaxationFactors']['fields'] = {'p': '0.3'}
        soldict['potentialFlow'] = {'nNonOrthogonalCorrectors': '10'}

    else:
        soldict['SIMPLE']['residualControl']['p_rgh'] = '{:.5f}'.format(node.presid)
        soldict['solvers'] = {}
        soldict['relaxationFactors']['fields'] = {'p_rgh': '0.3'}

        if not features['buoss']:
            soldict['solvers']['p_rgh'] = {'solver': 'GAMG', 'smoother': 'DICGaussSeidel', 'tolerance': '1e-6', 'relTol': '0.01'}
            soldict['solvers']['"(U|h|k|epsilon|omega|nuTilda)"'] = {'solver': 'PBiCGStab', 'preconditioner': 'DILU', 'tolerance': '1e-07', 'relTol': '0.01'}
            soldict['SIMPLE']['residualControl']['h'] = '{:.5f}'.format(node.enresid)
            soldict['relaxationFactors']['equations']['rho'] = '1'
            soldict['relaxationFactors']['fields']['p_rgh'] = '0.3'
            soldict['relaxationFactors']['equations']['U'] = '0.7'
            soldict['relaxationFactors']['equations']['h'] = '0.7'
        else:
            soldict['solvers']['p_rgh'] = {'solver': 'PCG', 'preconditioner': 'DIC', 'tolerance': '1e-06', 'relTol': '0.01'}
            soldict['solvers']['"(U|e|k|epsilon|omega|nuTilda)"'] = {'solver': 'PBiCGStab', 'preconditioner': 'DILU', 'tolerance': '1e-07', 'relTol': '0.1'}
            soldict['SIMPLE']['residualControl']['e'] = '{:.5f}'.format(node.enresid)
            soldict['relaxationFactors']['fields'] = {'p_rgh': '0.7'}
            soldict['relaxationFactors']['equations']['e'] = '0.1'

    if features['rad']:
        soldict['solvers']['G'] = {'$p_rgh': '', 'tolerance' : '1e-05', 'relTol': '0.1'}
        soldict['SIMPLE']['residualControl']['G'] = '1e-3'
        soldict['relaxationFactors']['equations']['G'] = '0.7'


    htext = ofheader + write_ffile('dictionary', 'system', 'fvSolution')
    # soldict = {'lsf': {'solvers': {'p': {'solver': 'PCG', 'preconditioner': 'DIC', 'tolerance': '1e-6', 'relTol': '0.05'},
    #                             'pFinal': {'$p': '', 'relTol': '0'},
    #                             'U': {'solver': 'smoothSolver', 'smoother': 'symGaussSeidel', 'tolerance': '1e-6', 'relTol': '0.05'}},
    #                     'SIMPLE': {'nNonOrthogonalCorrectors': '0', 'pRefCell': '0', 'pRefValue': '0',
    #                            'residualControl': {'p': '{:.5f}'.format(node.presid), 'U': '{:.5f}'.format(node.uresid), '"(k|epsilon|omega)"': '{:.5f}'.format(node.keoresid)}},
    #                            'potentialFlow': {'nNonOrthogonalCorrectors': '10'},
    #                            'relaxationFactors': {'fields': {'p': '0.3'}, 'equations': {'U': '0.7', '"(k|epsilon|omega)"': '0.7'}}},
    #             'sf': {'solvers': {'p': {'solver': 'GAMG', 'smoother': 'GaussSeidel', 'tolerance': '1e-6', 'relTol': '0.1'},
    #                                'U': {'solver': 'smoothSolver', 'smoother': 'symGaussSeidel', 'tolerance': '1e-6', 'relTol': '0.05'},
    #                                 '"k|epsilon|omega"': {'solver': 'smoothSolver', 'smoother': 'symGaussSeidel', 'tolerance': '1e-6', 'relTol': '0.1'}},
    #                     'SIMPLE': {'nNonOrthogonalCorrectors': '0', 'pRefCell': '0', 'pRefValue': '0',
    #                            'residualControl': {'p': '{:.5f}'.format(node.presid), 'U': '{:.5f}'.format(node.uresid), '"(k|epsilon|omega)"': '{:.5f}'.format(node.keoresid)}},
    #                            'potentialFlow': {'nNonOrthogonalCorrectors': '10'},
    #                            'relaxationFactors': {'fields': {'p': '0.3'}, 'equations': {'U': '0.7', '"(k|epsilon|omega)"': '0.7'}}},
    #             'bsf': {'solvers': {'p_rgh': {'solver': 'GAMG', 'smoother': 'DICGaussSeidel', 'tolerance': '1e-6', 'relTol': '0.01'},
    #                                 '"(U|h|k|epsilon|omega)"': {'solver': 'PBiCGStab', 'preconditioner': 'DILU', 'tolerance': '1e-07', 'relTol': '0.01'}},
    #                     'SIMPLE': {'momentumPredictor': 'no', 'nNonOrthogonalCorrectors': '0', 'pRefCell': '0', 'pRefValue': '0',
    #                     'residualControl': {'p_rgh': '{:.5f}'.format(node.presid), 'U': '{:.5f}'.format(node.uresid), 'h': '{:.5f}'.format(node.enresid), '"(k|epsilon|omega)"': '{:.5f}'.format(node.keoresid)}},
    #                     'relaxationFactors': {'rho': '1', 'p_rgh': '0.7', 'U': '0.3', 'h': '0.7', '"(k|epsilon|omega)"': '0.3'}},
    #             'bbsf': {'solvers': {'p_rgh': {'solver': 'PCG', 'preconditioner': 'DIC', 'tolerance': '1e-06', 'relTol': '0.01'},
    #                                 '"(U|e|k|epsilon|omega)"': {'solver': 'PBiCGStab', 'preconditioner': 'DILU', 'tolerance': '1e-07', 'relTol': '0.1'}},
    #                     'SIMPLE': {'nNonOrthogonalCorrectors': '0',
    #                     'residualControl': {'p_rgh': '{:.5f}'.format(node.presid), 'U': '{:.5f}'.format(node.uresid), 'e': '{:.5f}'.format(node.enresid), '"(k|epsilon|omega)"': '{:.5f}'.format(node.keoresid)}},
    #                     'relaxationFactors': {'fields': {'p_rgh': '0.7'}, 'equations':{'U': '0.2', 'e': '0.1', '"(k|epsilon|R|omega)"': '0.7'}}}}

    # if node.buoyancy and node.radiation:
    #     soldict[solver]['solvers']['G'] = {'$p_rgh': '', 'tolerance' : '1e-05', 'relTol': '0.1'}
    #     soldict[solver]['solvers']['p_rgh'] = {'solver': 'PCG', 'preconditioner': 'DIC', 'tolerance' : '1e-06', 'relTol': '0.01'}
    #     soldict[solver]['SIMPLE']['residualControl']['G'] = '1e-3'
    #     soldict[solver]['relaxationFactors']= {'fields': {'rho': '1.0', 'p_rgh': '0.7'}, 'equations': {'U': '0.2', 'h': '0.2', '"(k|epsilon|omega|R)"': '0.5', 'G': '0.7'}}

    # return write_fvdict(htext, soldict[solver])
    return write_fvdict(htext, soldict)


def fvtppwrite(node, features):
    htext = ofheader + write_ffile('dictionary', 'constant', 'thermophysicalProperties')

    if not features['buoss']:
        tppdict = {'thermoType': {'type': 'heRhoThermo', 'mixture': 'pureMixture',
                                  'transport': 'const', 'thermo': 'hConst', 'equationOfState': 'perfectGas',
                                  'specie': 'specie', 'energy': 'sensibleEnthalpy'},
                   'mixture': {'specie': {'molWeight': '28.96'},
                               'thermodynamics': {'Cp': '1004.4', 'Hf': '0'},
                               'transport': {'mu': '1e-05', 'Pr': '0.7'}}}

    else:
        tppdict = {'thermoType': {'type': 'heRhoThermo', 'mixture': 'pureMixture',
                    'transport': 'const', 'thermo': 'eConst', 'equationOfState': 'Boussinesq',
                    'specie': 'specie', 'energy': 'sensibleInternalEnergy'},
                    'mixture': {'specie': {'molWeight': '28.96'},
                    'equationOfState': {'rho0': '1', 'T0': '300', 'beta': '3e-03'},
                    'thermodynamics': {'Cv': '712', 'Hf': '0'},
                    'transport': {'mu': '1e-05', 'Pr': '0.7'}}}
    return write_fvdict(htext, tppdict)


def fvmtwrite(node, features):
    mtdict = {}
    mtdict['RAS'] = {}

    if not features['turb']:
        mtdict['simulationType'] = 'laminar'
        mtdict['RAS']['model'] = 'laminar'
        mtdict['RAS']['turbulence'] = 'off'
        mtdict['RAS']['printCoeffs'] = 'off'

    else:
        mtdict['simulationType'] = 'RAS'
        mtdict['RAS']['model'] = node.turbulence
        mtdict['RAS']['turbulence'] = 'on'
        mtdict['RAS']['printCoeffs'] = 'on'

    htext = ofheader + write_ffile('dictionary', 'constant', ('RASProperties', 'momentumTransport')[features['buoss']])
    # mtdict = {'bbsf': {'simulationType': 'RAS', 'RAS': {'model': '{}'.format(node.turbulence), 'turbulence': 'on', 'printCoeffs': 'on'}},
    #             'sf': {'simulationType': 'RAS', 'RAS': {'model': '{}'.format(node.turbulence), 'turbulence': 'on', 'printCoeffs': 'on'}},
    #             'bsf':{'simulationType': 'RAS', 'RAS': {'model': '{}'.format(node.turbulence), 'turbulence': 'on', 'printCoeffs': 'on'}},
    #             'lsf':{'simulationType': 'laminar', 'RAS': {'model': '{}'.format(node.turbulence), 'turbulence': 'off', 'printCoeffs': 'off'}}}
    return write_fvdict(htext, mtdict)

    # {'simulationType': 'RAS', 'p_rgh': {'solver': 'PCG', 'preconditioner': 'DIC', 'tolerance': '1e-8', 'relTol': '0.01'},
    #                  '"(U|e|k|omega|epsilon)"': {'solver': 'PBiCGStab', 'preconditioner':'DILU;', 'tolerance': '1e-7', 'relTol': '0.1'}}
    # text = ofheader + write_ffile('dictionary', 'constant', ('momentumTransport', 'RASProperties')[solver == 'bsf']) + '\nsimulationType RAS;\n\n'

    # for sol in ras[solver]:
    #     if sol != 'names':
    #         text += sol + '\n  {\n' + ';\n'.join(['    {} {}'.format(s, ras[solver][sol][si]) for si, s in enumerate(ras[solver]['names'])]) + ';\n}'
    # return text


def fvtpwrite():
    htext = ofheader + write_ffile('dictionary', 'constant', 'transportProperties')
    tpdict = {'transportModel': 'Newtonian', 'nu': '[0 2 -1 0 0 0 0] 1e-05'}
    # tpdict = {'lsf': {'transportModel': 'Newtonian',
    #           'nu': '[0 2 -1 0 0 0 0] 1e-05'},
    #           'sf': {'transportModel': 'Newtonian',
    #           'nu': '[0 2 -1 0 0 0 0] 1e-05'}}
    return write_fvdict(htext, tpdict)


def fvrpwrite(node):
    htext = ofheader + write_ffile('dictionary', 'constant', 'radiationProperties')
    raddict = {'0': {'radiation': 'on', 'radiationModel': 'P1',
                'solverFreq': '1', 'absorptionEmissionModel': 'constant',
                'constantCoeffs': {'absorptivity': '0.5', 'emissivity': '0.5', 'E': '0'},
                'scatterModel': 'none', 'sootModel': 'none'},
                '1': {'radiation': 'on', 'radiationModel': 'fvDOM',
                'fvDOMCoeffs': {'nPhi': '0.5', 'nTheta': '0.5', 'tolerance': '1e-3', 'maxIter': '10'},
                'solverFreq':'10', 'absorptionEmissionModel':'constant',
                'constantCoeffs': {'absorptivity': '0.5', 'emissivity': '0.5', 'E': '0'},
                'scatterModel': 'none', 'sootModel': 'none'}}
    if node.solar:
        raddict[node.radmodel]['SolarLoadCoeffs'] = {'sunDirectionModel': 'sunDirConstant',
        'sunDirection': '({0[0]} {0[1]} {0[2]})'.format(-bpy.data.objects[node.sun].location), 'gridUp': '(0 0 1)', 'gridEast': '(1 0 0)',
        'sunLoadModel': 'sunLoadFairWeatherConditions', 'skyCloudCoverFraction': '0', 'A': '1088', 'B': '0.205', 'groundReflectivity': '0.2',
        'C': '0.134', 'useReflectedRays': 'true', 'reflecting': {'nPhi': '10', 'nTheta': '10'}, 'absorptionEmissionModel': 'none',
	    'scatterModel': 'none', 'sootModel': 'none'}

#    soldict = {'SolarLoadCoeffs': {'sunDirectionModel': 'sunDirConstant', 'sunDirection': '({0[0]} {0[1]} {0[2]})'.format(-bpy.data.lights[node.sun].location)}}

    return (write_fvdict(htext, raddict[node.radmodel]))
    # raddict = p1dict if node.radmodel == '0' else fvdomdict


    # for d in raddict:
    #     if isinstance(raddict[d], str):
    #         text += '{} {};\n\n'.format(d,raddict[d])
    #     elif isinstance(raddict[d], dict):
    #         text += '{}\n{{\n'.format(d)
    #         for sd in raddict[d]:
    #             if isinstance(raddict[d][sd], str):
    #                 text += '    {} {};\n'.format(sd, raddict[d][sd])
    #         text += '}\n'
    # return text


def fvschwrite(node, features):
    scdict = {}
    scdict['ddtSchemes'] = {'default': 'steadyState'}
    scdict['gradSchemes'] = {}
    scdict['gradSchemes']['default'] =  'Gauss linear'
    scdict['divSchemes'] = {}
    scdict['divSchemes']['default'] = 'none'
    scdict['laplacianSchemes'] = {}
    scdict['interpolationSchemes'] = {}
    scdict['interpolationSchemes']['default'] = 'linear'
    scdict['snGradSchemes'] = {}
    scdict['snGradSchemes']['default'] = 'corrected'

    if not features['turb']:
#        scdict['gradSchemes']['limited'] = 'cellLimited Gauss linear 1'
        scdict['divSchemes']['div((nuEff*dev2(T(grad(U)))))'] =  'Gauss linear'
        scdict['divSchemes']['div(phi,U)'] =  'bounded Gauss linearUpwind limited'
        scdict['laplacianSchemes']['default'] = 'Gauss linear corrected'

    else:
#        scdict['gradSchemes']['limited'] = 'cellLimited Gauss linear 1'
#        scdict['divSchemes']['turbulence'] = 'bounded Gauss limitedLinear 1'
        scdict['divSchemes']['div((nuEff*dev2(T(grad(U)))))'] = 'Gauss linear'
        scdict['divSchemes']['div(phi,U)'] = 'bounded Gauss linearUpwind limited'

        if features['turb'] == 'kE':
            scdict['divSchemes']['div(phi,k)'] = 'bounded Gauss limitedLinear 1'
            scdict['divSchemes']['div(phi,epsilon)'] = 'bounded Gauss limitedLinear 1'
#            scdict['gradSchemes']['grad(k)'] = '$limited'
#            scdict['gradSchemes']['grad(epsilon)'] = '$limited'
        elif features['turb'] == 'kO':
            scdict['divSchemes']['div(phi,k)'] = 'bounded Gauss limitedLinear 1'
            scdict['divSchemes']['div(phi,omega)'] = 'bounded Gauss limitedLinear 1'
#            scdict['gradSchemes']['grad(k)'] = '$limited'
#            scdict['gradSchemes']['grad(omega)'] = '$limited'
        elif features['turb'] == 'sA':
            scdict['divSchemes']['div(phi,nuTilda)'] = 'bounded Gauss limitedLinear 1'
            scdict['gradSchemes']['grad(nuTilda)'] = 'bounded Gauss limitedLinear 1'

        scdict['laplacianSchemes']['default'] = 'Gauss linear corrected'


#    if not features['buoy']:
#        scdict['gradSchemes']['grad(U)'] = '$limited'

    if features['buoy']:
        scdict['divSchemes']['div(((rho*nuEff)*dev2(T(grad(U)))))'] = 'Gauss linear'

        if not features['buoss']:
#            scdict['divSchemes']['turbulence'] = 'bounded Gauss upwind'
#            scdict['divSchemes']['div(phi,U)'] = 'bounded Gauss limitedLinear 0.2'
            scdict['divSchemes']['div(phi,K)'] = 'bounded Gauss limitedLinear 0.2'
            scdict['divSchemes']['div(phi,h)'] = 'bounded Gauss limitedLinear 0.2'

        else:
            scdict['divSchemes']['turbulence'] = 'bounded Gauss limitedLinear 1'
            scdict['divSchemes']['div(phi,e)'] = 'bounded Gauss upwind'
#            scdict['divSchemes']['div(phi,U)'] = 'bounded Gauss upwind',
            scdict['divSchemes']['div(phi,Ekp)'] = 'bounded Gauss linear'
#            scdict['laplacianSchemes']['default'] = 'Gauss linear corrected'

    htext = ofheader + write_ffile('dictionary', 'system', 'fvSchemes')
    return write_fvdict(htext, scdict)


def fvgwrite():
    htext = ofheader + write_ffile('uniformDimensionedVectorField', '"constant"', 'g')
    gdict = {'dimensions': '[0 1 -2 0 0 0 0]', 'value': '(0 0 -9.81)'}
    return write_fvdict(htext, gdict)


def fvshmlayers(oname, node):
    surfdict = {"0": (("firstLayerThickness", node.frlayer), ("thickness", node.olayer)),
                "1": (("firstLayerThickness", node.frlayer), ("expansionRatio", node.expansion)),
                "2": (("finalLayerThickness", node.fnlayer), ("expansionRatio", node.expansion)),
                "3": (("finalLayerThickness", node.fnlayer), ("thickness", node.olayer)),
                "4": (("thickness", node.olayer), ("expansionRatio", node.expansion))}


    return 'addLayersControls\n{{\n  relativeSizes true;\n  layers\n  {{\n    "{}.*"\n    {{\n      nSurfaceLayers {};\n    }}\n  }}\n\n'.format(oname, node.layers)
    '  expansionRatio 1.0;\n  finalLayerThickness 0.3;\n  minThickness 0.1;\n  nGrow 0;\n  featureAngle 60;\n  slipFeatureAngle 30;\n  nRelaxIter 3;\n  nSmoothSurfaceNormals 1;\n  nSmoothNormals 3;\n' + \
    '  nSmoothThickness 10;\n  maxFaceThicknessRatio 0.5;\n  maxThicknessToMedialRatio 0.3;\n  minMedianAxisAngle 90;\n  nBufferCellsNoExtrude 0;\n  nLayerIter 50;\n}\n\n'


def fvshmwrite(node, fvos, bmo, **kwargs):
    surfdict = {"0": ("firstLayerThickness", node.frlayer, "thickness", node.olayer),
                "1": ("firstLayerThickness", node.frlayer, "expansionRatio", node.expansion),
                "2": ("finalLayerThickness", node.fnlayer, "expansionRatio", node.expansion),
                "3": ("finalLayerThickness", node.fnlayer, "thickness", node.olayer),
                "4": ("thickness", node.olayer, "expansionRatio", node.expansion)}

#    layersurf = '({}|{})'.format(kwargs['ground'][0].name, fvos[0].name) if kwargs and kwargs['ground'] else fvos[0].name
    ofheader = 'FoamFile\n{\n    version     2.0;\n    format      ascii;\n    class       dictionary;\n    object      snappyHexMeshDict;\n}\n\n'
    ofheader += 'castellatedMesh    {};\nsnap    {};\naddLayers    {};\ndebug    {};\n\n'.format('true', 'true', ('false', 'true')[node.layers], 0)

    ofheader += 'geometry\n{\n'

    for o in fvos:
        ofheader += '    {0}\n    {{\n        type triSurfaceMesh;\n        file "{0}.obj";\n    \n}}'.format(o.name)

    ofheader += '};\n\n'
    ofheader += 'castellatedMeshControls\n{{\n  maxLocalCells {};\n  maxGlobalCells {};\n  minRefinementCells {};\n  maxLoadUnbalance 0.10;\n  nCellsBetweenLevels {};\n\n'.format(node.lcells, node.gcells, int(node.gcells/100), node.ncellsbl)
    ofheader += '  features\n  (\n'

    for o in fvos:
        ofheader += '    {{\n      file "{}.eMesh";\n      level {};\n    }}\n\n'.format(o.name, o.flovi_fl)

    ofheader += ');\n\n'
    ofheader += '  refinementSurfaces\n  {\n'

    for o in fvos:
        ofheader += '    {}\n    {{\n      level ({} {});\n    }}\n\n  '.format(o.name, o.flovi_slmin, o.flovi_slmax)

    ofheader += '};\n\n'
    ofheader += '  resolveFeatureAngle 30;\n  refinementRegions\n  {}\n\n'
    ofheader += '  locationInMesh ({0[0]:} {0[1]} {0[2]});\n  allowFreeStandingZoneFaces true;\n}}\n\n'.format(mathutils.Matrix.Translation(bmo['flovi_translate']) * bpy.data.objects[node.empties].location)
    ofheader += 'snapControls\n{\n  nSmoothPatch 3;\n  tolerance 2.0;\n  nSolveIter 30;\n  nRelaxIter 5;\n  nFeatureSnapIter 10;\n  implicitFeatureSnap false;\n  explicitFeatureSnap true;\n  multiRegionFeatureSnap false;\n}\n\n'
    ofheader += 'addLayersControls\n{\n  relativeSizes true;\n  layers\n  {\n'

    for o in fvos:
        ofheader += '"{}.*"\n    {{\n      nSurfaceLayers {};\n    }}\n'.format(o.name, o.flovi_sl)

    ofheader += '}}\n\n'.format(o.name, node.layers)
    ofheader += '  {0[0]} {0[1]};\n  {0[2]} {0[3]};\n  minThickness 0.1;\n  nGrow 0;\n  featureAngle 60;\n  slipFeatureAngle 30;\n  nRelaxIter 5;\n  nSmoothSurfaceNormals 1;\n  nSmoothNormals 3;\n'.format(surfdict[node.layerspec][:]) + \
                '  nSmoothThickness 10;\n  maxFaceThicknessRatio 0.5;\n  maxThicknessToMedialRatio 0.3;\n  minMedianAxisAngle 90;\n  nBufferCellsNoExtrude 0;\n  nLayerIter 50;\n}\n\n'
    ofheader += 'meshQualityControls\n{\n  #include "meshQualityDict"\n  nSmoothScale 4;\n  errorReduction 0.75;\n}\n\n'
    ofheader += 'writeFlags\n(\n  scalarLevels\n  layerSets\n  layerFields\n);\n\nmergeTolerance 1e-6;\n'
    return ofheader


def fvdcpwrite(p):
    body = 'numberOfSubdomains {0};\n\nmethod          simple;\n\nsimpleCoeffs\n{{\n    n               ({0} 1 1);\n    delta           0.001;\n}}\n\nhierarchicalCoeffs\n{{\n    n               (1 1 1);\n    delta           0.001;\n    order           xyz;\n}}\n\nmanualCoeffs\n{{\n    dataFile        "";\n}}\ndistributed     no;\nroots           ( );'.format(p)
    return ofheader + write_ffile("dictionary", "system", "decomposeParDict") + body


def fvmqwrite():
    ofheader = 'FoamFile\n{\n  version     2.0;\n  format      ascii;\n  class       dictionary;\n  object      meshQualityDict;\n}\n\n'
    ofheader += '#include "$WM_PROJECT_DIR/etc/caseDicts/mesh/generation/meshQualityDict"'
    return ofheader


def fvsfewrite(fvos):
    ofheader = 'FoamFile\n{\n  version     2.0;\n  format      ascii;\n  class       dictionary;\n  object      surfaceFeatureExtractDict;\n}\n\n'
    for o in fvos:
        ofheader += '{}.obj\n{{\n  extractionMethod    extractFromSurface;\n\n  extractFromSurfaceCoeffs\n  {{\n    includedAngle   150;\n  }}\n\n    writeObj\n    yes;\n}}\n'.format(o.name)
    return ofheader


def fvobjwrite(scene, fvos, bmo):
    objheader = '# FloVi obj exporter\n'
    svp = scene.vi_params

    for frame in range(svp['flparams']['start_frame'], svp['flparams']['end_frame'] + 1):
        for o in fvos:
            with open(os.path.join(scene['flparams']['offilebase'], str(frame), 'constant', 'trisurface', '{}.obj'.format(o.name)), 'w') as objfile:
                bm = bmesh.new()
                tempmesh = o.to_mesh(scene=scene, apply_modifiers=True, settings='PREVIEW')
                bm.from_mesh(tempmesh)
                bm.transform(o.matrix_world)
                bm.transform(mathutils.Matrix.Translation(bmo['flovi_translate']))
                bpy.data.meshes.remove(tempmesh)
                vcos = ''.join(['v {0[0]} {0[1]} {0[2]}\n'.format(v.co) for v in bm.verts])
                objfile.write(objheader+vcos)

                for m, mat in enumerate(o.data.materials):
                    objfile.write('g {}\n'.format(mat.name) + ''.join(['f {} {} {}\n'.format(*[v.index + 1 for v in f.verts]) for f in bmesh.ops.triangulate(bm, faces=bm.faces)['faces'] if f.material_index == m]))

                objfile.write('#{}'.format(len(bm.faces)))
                bm.free()


def oftomesh(ofb, vl, fomats, st, ns, nf):
    mesh = bpy.data.meshes.new("mesh")
    vcoords = []
    findices = []
    fi = []
    fn = 0
    prevline = ''

    with open(os.path.join(ofb, st, 'polyMesh', 'points'), 'r') as mfile:
        for line in mfile.readlines():
            if '(' in line and ')' in line:
                vcoords.append([float(x) for x in line.split('(')[1].split(')')[0].split()])

    with open(os.path.join(ofb, st, 'polyMesh', 'faces'), 'r') as mfile:
        for line in mfile.readlines():
            if line:
                if fn:
                    try:
                        fi.append(int(line))
                        fn -= 1
                    except:
                        pass

                if not fn and fi:
                    findices.append(fi)
                    fi = []
                    fn = 0
                elif '(' in line and ')' in line:
                    findices.append([int(x) for x in line.split('(')[1].split(')')[0].split()])

                else:
                    try:
                        if prevline == '\n' and int(line) < 100:
                            fn = int(line)
                    except:
                        fn = 0
            prevline = line

    mesh.from_pydata(vcoords, [], findices)
    o = bpy.data.objects.new('Mesh', mesh)
    o.vi_params.vi_type_string = 'FloVi Mesh'
    o.show_wire = True
    bpy.context.view_layer.active_layer_collection.collection.objects.link(o)
    selobj(vl, o)

    for mat in fomats:
        bpy.ops.object.material_slot_add()
        o.material_slots[-1].material = mat

    bpy.ops.object.material_slot_add()
    o.material_slots[-1].material = bpy.data.materials.new("Volume") if 'Volume' not in [m.name for m in bpy.data.materials] else bpy.data.materials["Volume"]

    for face in o.data.polygons:
        mi = 0
        for ni, n in enumerate(ns):
            if face.index >= n and face.index <= n + nf[ni]:
                face.material_index = ni
                mi = 1
        if not mi:
            face.material_index = len(ns)

# scdict = {'lsf': {'ddtSchemes': {'default': 'steadyState'},
    #                   'gradSchemes': {'default': 'Gauss linear',
    #                                   'limited': 'cellLimited Gauss linear 1',
    #                                   'grad(U)': '$limited'},
    #                   'divSchemes': {'default': 'none',
    #                                  'div((nuEff*dev2(T(grad(U)))))': 'Gauss linear',
    #                                  'div(phi,U)': 'bounded Gauss linearUpwind limited'},
    #                'laplacianSchemes': {'default': 'Gauss linear corrected'},
    #                'interpolationSchemes': {'default': 'linear'},
    #                'snGradSchemes': {'default': 'corrected'}},
    #             'sf': {'ddtSchemes': {'default': 'steadyState'},
    #                 'gradSchemes': {'default': 'Gauss linear', 'limited': 'cellLimited Gauss linear 1', 'grad(U)': '$limited'},
    #                 'divSchemes': {'default': 'none',
    #                     'turbulence': 'bounded Gauss limitedLinear 1',
    #                     'div((nuEff*dev2(T(grad(U)))))': 'Gauss linear', 'div(phi,U)': 'bounded Gauss linearUpwind limited'},
    #                'laplacianSchemes': {'default': 'Gauss linear corrected'},
    #                'interpolationSchemes': {'default': 'linear'},
    #                'snGradSchemes': {'default': 'corrected'}},
    #             'bsf': {'ddtSchemes': {'default': 'steadyState'},
    #                 'gradSchemes': {'default': 'Gauss linear'},
    #                 'divSchemes': {'default': 'none',
    #                     'turbulence': 'bounded Gauss upwind',
    #                     'div(((rho*nuEff)*dev2(T(grad(U)))))': 'Gauss linear',
    #                     'div(phi,U)': 'bounded Gauss limitedLinear 0.2', 'div(phi,K)': 'bounded Gauss limitedLinear 0.2',
    #                     'div(phi,h)': 'bounded Gauss limitedLinear 0.2'},
    #                 'laplacianSchemes': {'default': 'Gauss linear uncorrected'},
    #                 'interpolationSchemes': {'default': 'linear'},
    #                 'snGradSchemes': {'default': 'corrected'}},
    #             'bbsf': {'ddtSchemes': {'default': 'steadyState'},
    #                 'gradSchemes': {'default': 'Gauss linear'},
    #                 'divSchemes': {'default': 'none',
    #                     'turbulence': 'bounded Gauss limitedLinear 1',
    #                     'div(phi,e)': 'bounded Gauss upwind',
    #                     'div(phi,U)': 'bounded Gauss upwind',
    #                     'div(((rho*nuEff)*dev2(T(grad(U)))))': 'Gauss linear', 'div(phi,Ekp)': 'bounded Gauss linear'},
    #                 'laplacianSchemes': {'default': 'Gauss linear corrected'},
    #                 'interpolationSchemes': {'default': 'linear'},
    #                 'snGradSchemes': {'default': 'corrected'}}}

    # if node.turbulence == 'kEpsilon':
    #     scdict[solver]['divSchemes']['div(phi,k)'] = '$turbulence'
    #     scdict[solver]['divSchemes']['div(phi,epsilon)'] = '$turbulence'
    #     if solver == 'sf':
    #         scdict[solver]['gradSchemes']['grad(k)'] = '$limited'
    #         scdict[solver]['gradSchemes']['grad(epsilon)'] = '$limited'
    # elif node.turbulence == 'kOmega':
    #     scdict[solver]['divSchemes']['div(phi,k)'] = '$turbulence'
    #     scdict[solver]['divSchemes']['div(phi,omega)'] = '$turbulence'
    #     if solver == 'sf':
    #         scdict[solver]['gradSchemes']['grad(k)'] = '$limited'
    #         scdict[solver]['gradSchemes']['grad(epsilon)'] = '$limited'
    # elif node.turbulence == 'SpalartAllmaras':
    #     scdict['divSchemes']['div(phi,nuTilda)'] = '$turbulence'
    #     if solver == 'sf':
    #         scdict[solver]['gradSchemes']['grad(nuTilda)'] = '$limited'

    # return write_fvdict(htext, scdict[solver])