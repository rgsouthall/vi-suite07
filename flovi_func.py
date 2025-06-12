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
from math import cos, sin, pi
from mathutils import Vector
from .vi_func import selobj
from .vi_dicts import flovi_b_dict, flovi_p_dict, flovi_u_dict, flovi_nut_dict, flovi_k_dict, flovi_epsilon_dict, flovi_t_dict, flovi_a_dict, flovi_prgh_dict, flovi_rad_dict, flovi_id_dict

ofheader = r'''/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM:    The Open Source CFD Toolbox        |
|  \\    /   O peration     | Version:     10                                 |
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


flovi_p_bounds = {'sf': {'0': ('zeroGradient', 'fixedValue', 'freestreamPressure', 'totalPressure', 'calculated', 'uniformFixedValue'),
                         '1': ('calculated', 'zeroGradient', 'fixedFluxPressure'),
                         '2': ('symmetry', ),
                         '3': ('empty', )},
                  'bf': {'0': ('calculated',),
                         '1': ('calculated',),
                         '2': ('symmetry',),
                         '3': ('empty',)}}

flovi_u_bounds = {'sf': {'0': ('zeroGradient', 'fixedValue', 'inletOutlet', 'outletInlet', 'freestream', 'pressureInletOutletVelocity', 'atmBoundaryLayerInletVelocity', 'freestreamVelocity'),
                         '1': ('noSlip', 'fixedValue', 'slip'),
                         '2': ('symmetry', ),
                         '3': ('empty',)},
                  'bf': {'0': ('zeroGradient', 'fixedValue', 'inletOutlet', 'outletInlet', 'freestream', 'pressureInletOutletVelocity', 'atmBoundaryLayerInletVelocity', 'freestreamVelocity'),
                         '1': ('noSlip', 'fixedValue', 'slip'),
                         '2': ('symmetry', ),
                         '3': ('empty',)}}

flovi_nut_bounds = {'sf': {'0': ['calculated'], '1': ['nutkWallFunction'], '2': ('symmetry', ), '3': ('empty',)},
                    'bf': {'0': ['calculated'], '1': ['nutkWallFunction'], '2': ('symmetry', ), '3': ('empty',)}}

flovi_nutilda_bounds = {'sf': {'0': ('zeroGradient', 'fixedValue'), '1': ('zeroGradient', 'fixedValue'), '2': ('symmetry', ), '3': ('empty',)},
                        'bf': {'0': ('zeroGradient', 'fixedValue'), '1': ('zeroGradient', 'fixedValue'), '2': ('symmetry', ), '3': ('empty',)}}

flovi_k_bounds = {'sf': {'0': ('fixedValue', 'inletOutlet'), '1': ['kqRWallFunction', 'zeroGradient'], '2': ('symmetry', ), '3': ('empty',)},
                  'bf': {'0': ('fixedValue', 'inletOutlet', 'turbulentIntensityKineticEnergyInlet'), '1': ['kqRWallFunction'], '2': ('symmetry', ), '3': ('empty',)}}

flovi_epsilon_bounds = {'sf': {'0': ('fixedValue', 'inletOutlet'), '1': ['epsilonWallFunction', 'zeroGradient'], '2': ('symmetry', ), '3': ('empty',)},
                        'bf': {'0': ('fixedValue', 'inletOutlet', 'turbulentMixingLengthDissipationRateInlet'), '1': ['epsilonWallFunction'], '2': ('symmetry',), '3': ('empty',)}}

flovi_omega_bounds = {'sf': {'0': ('zeroGradient', 'fixedValue'), '1': ['omegaWallFunction'], '2': ('symmetry', ), '3': ('empty',)},
                      'bf': {'0': ('zeroGradient', 'fixedValue'), '1': ['omegaWallFunction'], '2': ('symmetry', ), '3': ('empty',)}}

flovi_t_bounds = {'sf': {'0': ('None', ), '1': ['None'], '2': ('None', ), '3': ('None',)},
                  'bf': {'0': ('zeroGradient', 'fixedValue', 'inletOutlet'), '1': ('zeroGradient', 'fixedValue'), '2': ('symmetry', ), '3': ('empty',)}}

flovi_prgh_bounds = {'sf': {'0': ('None', ), '1': ('None', ), '2': ('None', ), '3': ('None',)},
                     'bf': {'0': ('totalPressure', 'fixedFluxPressure', 'fixedValue', 'prghPressure', 'prghTotalPressure', 'prghEntrainmentPressure'),
                            '1': ('fixedFluxPressure', 'fixedValue', 'totalPressure'),
                            '2': ('symmetry', ), '3': ('empty',)}}

flovi_a_bounds = {'sf': {'0': ('None', ), '1': ['None'], '2': ('None', ), '3': ('None',)},
                  'bf': {'0': ('calculated', 'inletOutlet'), '1': ('compressible::alphatWallFunction', 'compressible::alphatJayatillekeWallFunction'), '2': ('symmetry', ), '3': ('empty',)}}

flovi_rad_bounds = {'sf': {'0': ('None', ), '1': ('None', ), '2': ('None', ), '3': ('None',)},
                    'bf': {'0': ('MarshakRadiation', 'zeroGradient'), '1': ('MarshakRadiation',), '2': ('symmetry', ), '3': ('empty',)}}

flovi_id_bounds = {'sf': {'0': ('None', ), '1': ('None', ), '2': ('None', ), '3': ('None',)},
                   'bf': {'0': ('greyDiffusiveRadiation',), '1': ('greyDiffusiveRadiation',), '2': ('symmetry',), '3': ('empty',)}}


def ret_fvb_menu(mat, context):
    svp = context.scene.vi_params

    if svp.get('flparams') and svp['flparams'].get('scenario') and svp['flparams']['scenario'] != '4':
        return [('{}'.format(b), '{}'.format(b), '{} boundary type'.format(b)) for b in flovi_p_dict[svp['flparams']['scenario']]]
    else:
        return [("0", "Patch", "Inlet/outlet"), ("1", "Wall", "Wall boundary"), ("2", "Symmetry", "Symmetry plane boundary"), ("3", "Empty", "Empty boundary")]


def ret_fvbp_menu(mat, context):
    svp = context.scene.vi_params

    if svp.get('flparams') and svp['flparams'].get('scenario') and svp['flparams']['scenario'] != '4':
        return [('{}'.format(b), '{}'.format(b), '{} boundary type'.format(b)) for b in flovi_p_dict[svp['flparams']['scenario']][mat.flovi_bmb_type]]
    else:
        return [('{}'.format(b), '{}'.format(b), '{} boundary type'.format(b)) for b in flovi_p_bounds[svp['flparams']['solver_type']][mat.flovi_bmb_type]]


def ret_fvbu_menu(mat, context):
    svp = context.scene.vi_params

    if svp.get('flparams') and svp['flparams'].get('scenario') and svp['flparams']['scenario'] != '4':
        return [('{}'.format(b), '{}'.format(b), '{} boundary type'.format(b)) for b in flovi_u_dict[svp['flparams']['scenario']][mat.flovi_bmb_type]]
    else:
        return [('{}'.format(b), '{}'.format(b), '{} boundary type'.format(b)) for b in flovi_u_bounds[svp['flparams']['solver_type']][mat.flovi_bmb_type]]


def ret_fvbnut_menu(mat, context):
    svp = context.scene.vi_params

    if svp.get('flparams') and svp['flparams'].get('scenario') and svp['flparams']['scenario'] != '4':
        return [('{}'.format(b), '{}'.format(b), '{} boundary type'.format(b)) for b in flovi_nut_dict[svp['flparams']['scenario']][mat.flovi_bmb_type]]
    else:
        return [('{}'.format(b), '{}'.format(b), '{} boundary type'.format(b)) for b in flovi_nut_bounds[svp['flparams']['solver_type']][mat.flovi_bmb_type]]


def ret_fvbk_menu(mat, context):
    svp = context.scene.vi_params

    if svp.get('flparams') and svp['flparams'].get('scenario') and svp['flparams']['scenario'] != '4':
        return [('{}'.format(b), '{}'.format(b), '{} boundary type'.format(b)) for b in flovi_k_dict[svp['flparams']['scenario']][mat.flovi_bmb_type]]
    else:
        return [('{}'.format(b), '{}'.format(b), '{} boundary type'.format(b)) for b in flovi_k_bounds[svp['flparams']['solver_type']][mat.flovi_bmb_type]]


def ret_fvbepsilon_menu(mat, context):
    svp = context.scene.vi_params

    if svp.get('flparams') and svp['flparams'].get('scenario') and svp['flparams']['scenario'] != '4':
        return [('{}'.format(b), '{}'.format(b), '{} boundary type'.format(b)) for b in flovi_epsilon_dict[svp['flparams']['scenario']][mat.flovi_bmb_type]]
    else:
        return [('{}'.format(b), '{}'.format(b), '{} boundary type'.format(b)) for b in flovi_epsilon_bounds[svp['flparams']['solver_type']][mat.flovi_bmb_type]]


def ret_fvbomega_menu(mat, context):
    return [('{}'.format(b), '{}'.format(b), '{} boundary type'.format(b)) for b in flovi_omega_bounds[context.scene.vi_params['flparams']['solver_type']][mat.flovi_bmb_type]]


def ret_fvbt_menu(mat, context):
    svp = context.scene.vi_params

    if svp.get('flparams') and svp['flparams'].get('scenario') and svp['flparams']['scenario'] in ('2', '3'):
        return [('{}'.format(b), '{}'.format(b), '{} boundary type'.format(b)) for b in flovi_t_dict[svp['flparams']['scenario']][mat.flovi_bmb_type]]
    elif svp.get('flparams') and svp['flparams'].get('scenario') and svp['flparams']['scenario'] == '4':
        return [('{}'.format(b), '{}'.format(b), '{} boundary type'.format(b)) for b in flovi_t_bounds[svp['flparams']['solver_type']][mat.flovi_bmb_type]]
    else:
        return [('None', 'None', 'None')]


def ret_fvba_menu(mat, context):
    svp = context.scene.vi_params

    if svp.get('flparams') and svp['flparams'].get('scenario') and svp['flparams']['scenario'] in ('2', '3'):
        return [('{}'.format(b), '{}'.format(b), '{} boundary type'.format(b)) for b in flovi_a_dict[svp['flparams']['scenario']][mat.flovi_bmb_type]]
    elif svp.get('flparams') and svp['flparams'].get('scenario') and svp['flparams']['scenario'] == '4':
        return [('{}'.format(b), '{}'.format(b), '{} boundary type'.format(b)) for b in flovi_a_bounds[svp['flparams']['solver_type']][mat.flovi_bmb_type]]
    else:
        return [('None', 'None', 'None')]


def ret_fvbprgh_menu(mat, context):
    svp = context.scene.vi_params

    if svp.get('flparams') and svp['flparams'].get('scenario') and svp['flparams']['scenario'] in ('2', '3'):
        return [('{}'.format(b), '{}'.format(b), '{} boundary type'.format(b)) for b in flovi_prgh_dict[svp['flparams']['scenario']][mat.flovi_bmb_type]]
    elif svp.get('flparams') and svp['flparams'].get('scenario') and svp['flparams']['scenario'] == '4':
        return [('{}'.format(b), '{}'.format(b), '{} boundary type'.format(b)) for b in flovi_prgh_bounds[svp['flparams']['solver_type']][mat.flovi_bmb_type]]
    else:
        return [('None', 'None', 'None')]


def ret_fvrad_menu(mat, context):
    svp = context.scene.vi_params

    if svp.get('flparams') and svp['flparams'].get('scenario') and svp['flparams']['scenario'] in ('2', '3'):
        return [('{}'.format(b), '{}'.format(b), '{} boundary type'.format(b)) for b in flovi_rad_dict[svp['flparams']['scenario']][mat.flovi_bmb_type]]

    elif svp.get('flparams') and svp['flparams'].get('scenario') and svp['flparams']['scenario'] == '4':
        return [('{}'.format(b), '{}'.format(b), '{} boundary type'.format(b)) for b in flovi_rad_bounds[svp['flparams']['solver_type']][mat.flovi_bmb_type]]

    else:
        return [('None', 'None', 'None')]

def ret_fvi_menu(mat, context):
    svp = context.scene.vi_params

    if svp.get('flparams') and svp['flparams'].get('scenario') and svp['flparams']['scenario'] in ('2', '3'):
        return [('{}'.format(b), '{}'.format(b), '{} boundary type'.format(b)) for b in flovi_id_dict[svp['flparams']['scenario']][mat.flovi_bmb_type]]
    elif svp.get('flparams') and svp['flparams'].get('scenario') and svp['flparams']['scenario'] == '4':
        return [('{}'.format(b), '{}'.format(b), '{} boundary type'.format(b)) for b in flovi_id_bounds[svp['flparams']['solver_type']][mat.flovi_bmb_type]]
    else:
        return [('None', 'None', 'None')]

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
        boundary += "  {}\n  {{\n    type {};\n    faces\n    (\n".format(mat.name, ("patch", "wall", "symmetry", "empty")[int(mat.flovi_bmb_type)])
        faces = [face for face in o.data.polygons if o.data.materials[face.material_index] == mat]

        for face in faces:
            boundary += "      ("+" ".join([str(v) for v in face.vertices])+")\n"

        boundary += "    );\n  }\n"
    boundary += ");\n\nmergePatchPairs\n(\n);"
    return boundary


def write_bound(mn, m, ns, nf):
    t = flovi_b_dict[m.vi_params.flovi_bmb_type]

    return '''   {0}
    {{
        type            {1};
        nFaces          {2};
        startFace       {3};
    }}\n'''.format(mn, t, nf, ns)


def fvmat(self, svp, mn, bound, frame):
    begin = '\n  {}\n  {{\n    type    '.format(mn)
    end = ';\n  }\n'
    entry = ''

    if bound == 'p':
        val = 'uniform {:.4f}'.format(self.flovi_bmbp_val) if not self.flovi_p_field else '$internalField'
        ptdict = {'zeroGradient': 'zeroGradient',
                  'fixedValue': 'fixedValue;\n    value    {}'.format(val),
                  'calculated': 'calculated;\n    value    $internalField',
                  'uniformFixedValue': 'uniformFixedValue;\n    uniformValue    constant {}'.format(self.flovi_bmbp_val),
                  'freestreamPressure': 'freestreamPressure;\n  freestreamValue {}'.format(val),
                  'fixedMeanOutletInlet': 'fixedMeanOutletInlet;\n    meanValue   {0};\n    value {0}'.format(val),
                  'fixedMean': 'fixedMean;\n    meanValue   {0};\n    value {0}'.format(val),
                  'totalPressure': 'totalPressure;\n    p0      uniform {};\n    value    {}'.format(self.flovi_bmbp_p0val, val),
                  'symmetry': 'symmetry', 'empty': 'empty'}

        entry = ptdict[self.flovi_bmbp_subtype]

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

        Utdict = {'fixedValue': 'fixedValue;\n    value    {}'.format(val), 'slip': 'slip', 'noSlip': 'noSlip',
                  'uniformFixedValue': 'uniformFixedValue;\n    value    {0};\n    uniformValue    {0}'.format(val),
                  'inletOutlet': 'inletOutlet;\n    inletValue    {0};\n    value    {0}'.format(val),
                  'outletInlet': 'outletInlet;\n    outletValue    {0};\n    value    {0}'.format(val),
                  'pressureInletOutletVelocity': 'pressureInletOutletVelocity;\n    value    {}'.format(val),
                  'zeroGradient': 'zeroGradient', 'symmetry': 'symmetry',
                  'freestream': 'freestream;\n    freestreamValue    $internalField',
                  'freestreamVelocity': 'freestreamVelocity;\n    freestreamValue    $internalField',
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

        if self.flovi_bmbu_subtype in Utdict:
            entry = Utdict[self.flovi_bmbu_subtype]

    elif bound == 'nut':
        ntdict = {'nutkWallFunction': 'nutkWallFunction;\n    value    $internalField',
                  'nutUSpaldingWallFunction': 'nutUSpaldingWallFunction;\n    value    $internalField',
                  'calculated': 'calculated;\n    value    $internalField',
                  'inletOutlet': 'inletOutlet;\n    inletValue    $internalField\n    value    $internalField',
                  'symmetry': 'symmetry', 'empty': 'empty'}

        entry = ntdict[self.flovi_bmbnut_subtype]

    elif bound == 'k':
        val = 'uniform {:.4f}'.format(self.flovi_k_val) if not self.flovi_k_field else '$internalField'
        ival = '{:.4f}'.format(self.flovi_k_intensity)
        ktdict = {'zeroGradient': 'zeroGradient',
                  'fixedValue': 'fixedValue;\n    value    {}'.format(val),
                  'kqRWallFunction': 'kqRWallFunction;\n    value    {}'.format(val),
                  'inletOutlet': 'inletOutlet;\n    inletValue    {}'.format(val),
                  'calculated': 'calculated;\n    value    {}'.format(val),
                  'symmetry': 'symmetry',
                  'empty': 'empty',
                  'turbulentIntensityKineticEnergyInlet': 'turbulentIntensityKineticEnergyInlet;\n    intensity       {};\n    value      {}'.format(ival, val)}

        entry = ktdict[self.flovi_k_subtype]

    elif bound == 't':
        val = 'uniform {}'.format(self.flovi_bmbt_val) if not self.flovi_t_field else '$internalField'
        ival = 'uniform {}'.format(self.flovi_bmbti_val) if not self.flovi_t_field else '$internalField'
        ttdict = {'zeroGradient': 'zeroGradient',
                  'fixedValue': 'fixedValue;\n    value    {}'.format(val),
                  'inletOutlet': 'inletOutlet;\n    inletValue    {}\n    value    {}'.format(ival, val),
                  'calculated': 'calculated;\n    value    $internalField',
                  'symmetry': 'symmetry',
                  'empty': 'empty'}
        entry = ttdict[self.flovi_bmbt_subtype]

    elif bound == 'p_rgh':
        val = 'uniform {:.4f}'.format(self.flovi_prgh_val) if not self.flovi_prgh_field else '$internalField'
        p0val = 'uniform {:.4f}'.format(self.flovi_prgh_p0) if not self.flovi_prgh_field else '$internalField'
        prghtdict = {'zeroGradient': 'zeroGradient',
                     'totalPressure': 'totalPressure;\n    rho  rho;\n    gamma   1;\n    p0    uniform 0;\n    value    {}'.format(val),
                     'inletOutlet': 'inletOutlet;\n    inletValue    $internalField;\n    value    $internalField',
                     'fixedFluxPressure': 'fixedFluxPressure;\n    value    {}'.format(val),
                     'fixedValue': 'fixedValue;\n    value    {}'.format(val),
                     'prghEntrainmentPressure': 'prghEntrainmentPressure;\n    p0       $internalField'.format(p0val),
                     'calculated': 'calculated;\n    value    $internalField',
                     'prghPressure': 'prghPressure;\n    p    $internalField;\n    value    $internalField',
                     'freestreamPressure': 'freestreamPressure',
                     'prghTotalPressure': 'prghTotalPressure;\n    p0      {};\n    value    {}'.format(p0val, val),
                     'fixedMeanOutletInlet': 'fixedMeanOutletInlet;\n    meanValue   {0};\n    value {0}'.format(val),
                     'fixedMean': 'fixedMean;\n    meanValue   {0};\n    value {0}'.format(val),
                     'symmetry': 'symmetry', 'empty': 'empty'}
        entry = prghtdict[self.flovi_prgh_subtype]

    elif bound == 'a':
        val = 'uniform {}'.format(self.flovi_a_val) if not self.flovi_a_field else '$internalField'
        ttdict = {'compressible::alphatWallFunction': 'compressible::alphatWallFunction;\n    value     $internalField',
                  'compressible::alphatJayatillekeWallFunction': 'compressible::alphatJayatillekeWallFunction;\n    Prt    0.85;\n    value     $internalField',
                  'fixedValue': 'fixedValue;\n    value    {}'.format(val),
                  'inletOutlet': 'inletOutlet;\n    inletValue    $internalField\n    value    $internalField',
                  'calculated': 'calculated;\n    value    $internalField', 'symmetry': 'symmetry', 'empty': 'empty'}

        entry = ttdict[self.flovi_a_subtype]

    elif bound == 'e':
        val = 'uniform {:.4f}'.format(self.flovi_bmbe_val) if not self.flovi_e_field else '$internalField'
        etdict = {'zeroGradient': 'zeroGradient',
                  'symmetry': 'symmetry',
                  'empty': 'empty',
                  'inletOutlet': 'inletOutlet;\n    inletValue    {}'.format(val),
                  'fixedValue': 'fixedValue;\n    value    {}'.format(val),
                  'epsilonWallFunction': 'epsilonWallFunction;\n    value    {}'.format(val),
                  'calculated': 'calculated;\n    value    {}'.format(val),
                  'turbulentMixingLengthDissipationRateInlet': 'turbulentMixingLengthDissipationRateInlet;\n    mixingLength  {:.5f};\n    value    {}'.format(self.flovi_eml_val, val)}

        entry = etdict[self.flovi_bmbe_subtype]

    elif bound == 'G':
        radtdict = {'MarshakRadiation': 'MarshakRadiation;\n    emissivityMode    {};\n    emissivity    uniform {};\n    value    uniform {}'.format(self.flovi_rad_em, self.flovi_rad_e, self.flovi_rad_val),
                    'zeroGradient': 'zeroGradient',
                    'symmetry': 'symmetry'}
        entry = radtdict[self.flovi_rad_subtype]

    elif bound == 'IDefault':
        idict = {'greyDiffusiveRadiation': 'greyDiffusiveRadiation;\n   emissivityMode  lookup;\n   emissivity      uniform 1.0;\n  value           uniform 0'}
        entry = idict[self.flovi_i_subtype]

    if entry:
        return begin + entry + end
    else:
        return ''


def fvvarwrite(scene, obs, node):
    '''Turbulence modelling: k and epsilon required for kEpsilon, k and omega required for kOmega, nutilda required for SpalartAllmaras, nut required for all
        Buoyancy modelling: T'''
    svp = scene.vi_params
    b_dict = {}

    for frame in range(svp['flparams']['start_frame'], svp['flparams']['end_frame'] + 1):
        scene.frame_set(frame)
        frame_of0fb = os.path.join(svp['flparams']['offilebase'], str(frame), '0')

        if not os.path.isdir(frame_of0fb):
            os.makedirs(frame_of0fb)
        else:
            for f in os.listdir(frame_of0fb):
                try:
                    os.remove(os.path.join(frame_of0fb, f))
                except Exception:
                    pass

        uval = node.uval if node.uval_type == '0' else Vector((sin(node.uval_azi * pi/180) * node.umag, cos(node.uval_azi * pi/180) * node.umag, 0))
        svp['flparams'][str(frame)] = {}
        svp['flparams'][str(frame)]['Udir'] = uval.normalized()
        svp['flparams'][str(frame)]['Uspeed'] = uval.length

        if not node.buoyancy:
            pentry = "dimensions [{} {} {} {} 0 0 0];\ninternalField   uniform {};\n\nboundaryField\n{{\n".format('0', '2', '-2', '0',
                                                                                                                '{}'.format(node.pnormval))
        else:
            pentry = "dimensions [{} {} {} {} 0 0 0];\ninternalField   uniform {};\n\nboundaryField\n{{\n".format('1', '-1', '-2', '0',
                                                                                                                '{}'.format(node.pabsval))

        (Uentry, nutildaentry, nutentry, kentry, eentry, oentry, tentry, p_rghentry, aentry, Gentry, Ientry) = ["dimensions [{} {} {} {} 0 0 0];\ninternalField   uniform {};\n\nboundaryField\n{{\n".format(*var) for var in (
                                                                                    ('0', '1', '-1', '0', '({:.4f} {:.4f} {:.4f})'.format(*uval)),
                                                                                    ('0', '2', '-1', '0', '{:.4f}'.format(node.nutildaval)),
                                                                                    ('0', '2', '-1', '0', '{:.4f}'.format(node.nutval)),
                                                                                    ('0', '2', '-2', '0', '{:.4f}'.format(node.kval)),
                                                                                    ('0', '2', '-3', '0', '{:.4f}'.format(node.epval)),
                                                                                    ('0', '0', '-1', '0', '{:.4f}'.format(node.oval)),
                                                                                    ('0', '0', '0', '1', '{:.4f}'.format(node.tval)),
                                                                                    ('1', '-1', '-2', '0', '{:.4f}'.format(node.p_rghval)),
                                                                                    ('1', '-1', '-1', '0', '{:.4f}'.format(node.aval)),
                                                                                    ('1', '0', '-3', '0', '{:.4f}'.format(node.Gval)),
                                                                                    ('1', '0', '-3', '0', '{:.4f}'.format(node.Gval)))]
        for o in obs:
            for mat in o.data.materials:
                mvp = mat.vi_params
                matname = '{}_{}'.format(o.name, mat.name)
                b_dict[matname] = {}

                if mvp.mattype == '2':
                    pentry += mvp.flovi_mat(svp, matname, 'p', frame)
                    b_dict[matname]['p'] = mvp.flovi_mat(svp, matname, 'p', frame)
                    Uentry += mvp.flovi_mat(svp, matname, 'U', frame)
                    b_dict[matname]['U'] = mvp.flovi_mat(svp, matname, 'U', frame)
                    nutentry += mvp.flovi_mat(svp, matname, 'nut', frame)
                    b_dict[matname]['nut'] = mvp.flovi_mat(svp, matname, 'nut', frame)
                    kentry += mvp.flovi_mat(svp, matname, 'k', frame)
                    b_dict[matname]['k'] = mvp.flovi_mat(svp, matname, 'k', frame)
                    eentry += mvp.flovi_mat(svp, matname, 'e', frame)
                    b_dict[matname]['e'] = mvp.flovi_mat(svp, matname, 'e', frame)

                    if node.buoyancy:
                        tentry += mvp.flovi_mat(svp, matname, 't', frame)
                        b_dict[matname]['t'] = mvp.flovi_mat(svp, matname, 't', frame)
                        p_rghentry += mvp.flovi_mat(svp, matname, 'p_rgh', frame)
                        b_dict[matname]['p_rgh'] = mvp.flovi_mat(svp, matname, 'p_rgh', frame)
                        aentry += mvp.flovi_mat(svp, matname, 'a', frame)
                        b_dict[matname]['a'] = mvp.flovi_mat(svp, matname, 'a', frame)

                        if node.radiation:
                            Gentry += mvp.flovi_mat(svp, matname, 'G', frame)
                            b_dict[matname]['G'] = mvp.flovi_mat(svp, matname, 'G', frame)
                            if node.radmodel == '1':
                                Ientry += mvp.flovi_mat(svp, matname, 'IDefault', frame)
                                b_dict[matname]['IDefault'] = mvp.flovi_mat(svp, matname, 'IDefault', frame)

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

        with open(os.path.join(frame_of0fb, 'nut'), 'w') as nutfile:
            nutfile.write(ofheader + write_ffile('volScalarField', '', 'nut') + nutentry)

        with open(os.path.join(frame_of0fb, 'k'), 'w') as kfile:
            kfile.write(ofheader + write_ffile('volScalarField', '', 'k') + kentry)
        with open(os.path.join(frame_of0fb, 'epsilon'), 'w') as efile:
            efile.write(ofheader + write_ffile('volScalarField', '', 'epsilon') + eentry)

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
                if node.radmodel == '1':
                    with open(os.path.join(frame_of0fb, 'IDefault'), 'w') as IDfile:
                        IDfile.write(ofheader + write_ffile('volScalarField', '', 'IDefault') + Ientry)
    return b_dict

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


def fvcdwrite(svp, node, dp):
    ps, ss, bs, htcs = [], [], [], []
    solver = svp['flparams']['solver']
    htext = ofheader + write_ffile('dictionary', 'system', 'controlDict')
    cdict = {'application': solver, 'startFrom': 'startTime', 'startTime': '0', 'stopAt': 'endTime',
             'endTime': f'{node.etime}', 'deltaT': f'{node.dtime:.5f}', 'writeControl': 'timeStep', 'writeInterval': f'{node.w_int}',
             'purgeWrite': '{}'.format(0), 'writeFormat': 'ascii', 'writePrecision': '6', 'writeCompression': 'off',
             'timeFormat': 'general', 'timePrecision': '6', 'runTimeModifiable': 'true', 'functions': {}, 'libs': '("libatmosphericModels.so")'}

    if node.comfort:
        cdict['functions']['comfort'] = {'libs': '("libfieldFunctionObjects.so")', 'type': 'comfort', 'clothing': f'{node.clo*0.155:.2f}',
                                         'metabolicRate': f'{node.met:.2f}', 'relHumidity': f'{node.rh * 0.01:.2f}', 
                                         'writeControl': 'writeTime', 'executeControl': 'writeTime'}
    if node.age:
        cdict['functions']['age'] = {'libs': '("libfieldFunctionObjects.so")', 'type': 'age', 'diffusion': 'on', 'writeControl': 'writeTime', 'executeControl': 'writeTime'}

    for o in bpy.data.objects:
        ovp = o.vi_params

        if o.type == 'MESH' and ovp.vi_type == '2':
            dom = o
        elif o.type == 'EMPTY' and ovp.flovi_probe:
            ps.append(o)
        elif o.type == 'EMPTY' and o.name == node.p_ref_point and node.p_ref != '0':
            ps.append(o)
        elif o.type == 'MESH' and ovp.vi_type == '6' and o.visible_get():
            for frame in range(svp['flparams']['start_frame'], svp['flparams']['end_frame'] + 1):
                if not os.path.isdir(os.path.join(svp['flparams']['offilebase'], str(frame), 'constant', 'triSurface')):
                    os.makedirs(os.path.join(svp['flparams']['offilebase'], str(frame), 'constant', 'triSurface'))
                ovp.write_stl(dp, os.path.join(svp['flparams']['offilebase'], str(frame), 'constant', 'triSurface', '{}.stl'.format(o.name)))
            
            ss.append(o.name)

        if o.type == 'MESH' and ovp.vi_type in ('2', '3') and any([m.vi_params.flovi_probe for m in o.data.materials]):
            for mat in o.data.materials:
                if mat.vi_params.flovi_probe:
                    bs.append('{}_{}'.format(o.name, mat.name))

        if o.type == 'MESH' and ovp.vi_type in ('2', '3') and any([m.vi_params.flovi_htc for m in o.data.materials]):
            for mat in o.data.materials:
                if mat.vi_params.flovi_htc and mat.vi_params.flovi_bmb_type in ('Solid', 'Wall'):
                    htcs.append('{}_{}'.format(o.name, mat.name))

    if ps:
        svp['flparams']['probes'] = [p.name.replace(" ", "_") for p in ps]
        probe_vars = 'p U T k epsilon'

        for p in ps:
            cdict['functions'][p.name.replace(" ", "_")] = {'libs': '("libsampling.so")', 'type': 'probes', 'name': '{}'.format(p.name.replace(" ", "_")), 'writeControl': 'timeStep',
                                          'writeInterval': f'{node.w_int}', 'fields': '({0})'.format(probe_vars),
                                          'probeLocations\n(\n{}\n)'.format('   ({0[0]} {0[1]} {0[2]})'.format(p.location)): ''}

        probe_text = '''functions
{{
    probes
    {{
        libs            ("libsampling.so");
        type            probes;
        name            {2};
        writeControl    timeStep;
        writeInterval   '{3}';
        fields          ({0});
        probeLocations
        (
            {1}
        );
    }}
}}'''.format(probe_vars, ''.join([' ({0[0]} {0[1]} {0[2]})\n'.format(p.location) for p in ps]), ','.join(['{}'.format(p.name.replace(" ", "_")) for p in ps]), node.w_int)

    else:
        probe_text = ''
        bpy.context.scene.vi_params['flparams']['probes'] = []

    bpy.context.scene.vi_params['flparams']['s_probes'] = ss
    bpy.context.scene.vi_params['flparams']['b_probes'] = bs

    if bs:
        for b in bs:
            cdict['functions'][b] = {'type': 'surfaceFieldValue', 'libs': '("libfieldFunctionObjects.so")', 'writeControl': 'timeStep',
                                     'writeInterval': f'{node.w_int}', 'writeFields': 'true', 'surfaceFormat': 'raw', 'regionType': 'patch', 'name': '{}'.format(b),
                                     'operation': 'areaAverage', 'fields    (p U)': ''}
            cdict['functions'][b+'_vf'] = {'type': 'surfaceFieldValue', 'libs': '("libfieldFunctionObjects.so")', 'writeControl': 'timeStep',
                                     'writeInterval': f'{node.w_int}', 'writeFields': 'true', 'surfaceFormat': 'raw', 'regionType': 'patch', 'name': '{}'.format(b),
                                     'operation': 'areaNormalIntegrate', 'fields    (U)': ''}
    if htcs:
        cdict['functions']['htc'] = {'type': 'wallHeatTransferCoeff', 'libs': '("libfieldFunctionObjects.so")', 'model': 'kappaEff', 'patches' : '({})'.format(' '.join(list(set(htcs)))), 'writeControl': 'timeStep',
                                     'writeInterval': f'{node.w_int}', 'rho': '1.225', 'Cp': '1005', 'Pr': '0.707', 'Prt': '0.9'}

    return write_fvdict(htext, cdict)


def fvprefwrite(node):
    htext = ofheader + write_ffile('uniformDimensionedScalarField', 'constant', 'pRef')
    pdict = {'dimensions': '[1 -1 -2 0 0 0 0]', 'value': '{}'.format(node.pabsval)}
    text = write_fvdict(htext, pdict)
    return text


def fvsolwrite(svp, node):
    sol_dict = {}
    sol_dict['relaxationFactors'] = {}
    sol_dict['relaxationFactors']['equations'] = {}
    sol_dict['relaxationFactors']['fields'] = {}
    sol_dict['solvers'] = {}

    if svp['flparams']['solver_type'] == 'bf':
        sol_sol_dict = sol_dict['PIMPLE'] = {}
    else:
        sol_sol_dict = sol_dict['SIMPLE'] = {}

    sol_sol_dict['residualControl'] = {}
    sol_sol_dict['residualControl']['U'] = '{:.5f}'.format(node.uresid)

    if svp['flparams']['solver_type'] == 'sf':
        sol_dict['solvers']['p'] = {'solver': 'GAMG', 'smoother': 'GaussSeidel', 'tolerance': '1e-6', 'relTol': '0.1'}
        sol_dict['solvers']['"U|k|epsilon"'] = {'solver': 'smoothSolver', 'smoother': 'symGaussSeidel', 'tolerance': '1e-6', 'relTol': '0.1'}
        sol_sol_dict['residualControl']['"(k|epsilon)"'] = '{:.5f}'.format(node.keoresid)
        sol_sol_dict['residualControl']['p'] = '{:.5f}'.format(node.presid)
        sol_sol_dict['residualControl']['U'] = '{:.5f}'.format(node.uresid)
        sol_sol_dict['nNonOrthogonalCorrectors'] = '1'
        sol_dict['potentialFlow'] = {'nNonOrthogonalCorrectors': '10'}
        sol_dict['relaxationFactors']['equations']['U'] = '0.7'
        sol_dict['relaxationFactors']['equations']['"(k|epsilon)"'] = '0.7'
        sol_dict['relaxationFactors']['fields'] = {'p': '0.3'}

        if node.p_ref != '0':
            sol_sol_dict['pRefValue'] = '{}'.format(node.p_ref_val)

            if node.p_ref:
                sol_sol_dict['pRefPoint'] = '({0[0]} {0[1]} {0[2]})'.format(bpy.data.objects[node.p_ref_point].location)
            else:
                sol_sol_dict['pRefCell'] = '0'

    if svp['flparams']['solver_type'] == 'bf':
        sol_sol_dict['nNonOrthogonalCorrectors'] = '1'
        sol_dict['relaxationFactors']['equations']['U'] = '0.2'
        sol_dict['relaxationFactors']['equations']['e'] = '0.2'
        sol_dict['relaxationFactors']['equations']['"(k|epsilon|R)"'] = '0.7'
        sol_dict['relaxationFactors']['fields'] = {'p_rgh': '0.7'}
        sol_sol_dict['residualControl']['"(k|epsilon)"'] = '{:.5f}'.format(node.keoresid)
        sol_sol_dict['residualControl']['e'] = '{:.5f}'.format(node.enresid)
        sol_sol_dict['residualControl']['p_rgh'] = '{:.5f}'.format(node.presid)

        if node.p_ref != '0':
            sol_sol_dict['pRefValue'] = '{}'.format(node.p_ref_val)

            if node.p_ref:
                sol_sol_dict['pRefPoint'] = '({0[0]} {0[1]} {0[2]})'.format(bpy.data.objects[node.p_ref_point].location)
            else:
                sol_sol_dict['pRefCell'] = '0'

        sol_dict['solvers']['p_rgh'] = {'solver': 'PCG', 'preconditioner': 'DIC', 'tolerance': '1e-8', 'relTol': '0.01'}
        sol_dict['solvers']['"U|e|k|epsilon"'] = {'solver': 'PBiCGStab', 'preconditioner': 'DILU', 'tolerance': '1e-07', 'relTol': '0.1'}

        if svp['flparams']['features']['rad']:
            if svp['flparams']['radmodel'] == '0':
                sol_dict['solvers']['G'] = {'$p_rgh': '', 'tolerance': '1e-05', 'relTol': '0.1'}
                sol_sol_dict['residualControl']['G'] = '1e-3'
                sol_dict['relaxationFactors']['equations']['G'] = '0.7'
            else:
                sol_dict['solvers']['"Ii.*"'] = {'solver': 'GAMG', 'tolerance': '1e-4', 'relTol': '0', 'smoother': 'symGaussSeidel', 'maxIter': '5', 'nPostSweeps': '1'}
                sol_sol_dict['residualControl']['ILambda.*'] = '1e-3'
                sol_dict['relaxationFactors']['equations']['ILambda.*'] = '0.7'
    
    if node.age:
        sol_dict['solvers']['age'] = {'$U': '', 'relTol': '0.001'}
    
    htext = ofheader + write_ffile('dictionary', 'system', 'fvSolution')
    return write_fvdict(htext, sol_dict)


def fvtppwrite(svp):
    htext = ofheader + write_ffile('dictionary', 'constant', 'physicalProperties')

    if svp['flparams']['solver_type'] == 'sf':
        tppdict = {'thermoType': {'type': 'heRhoThermo', 'mixture': 'pureMixture',
                                  'transport': 'const', 'thermo': 'hConst', 'equationOfState': 'perfectGas',
                                  'specie': 'specie', 'energy': 'sensibleEnthalpy'},
                   'mixture': {'specie': {'molWeight': '28.96'},
                               'thermodynamics': {'Cp': '1004.4', 'Hf': '0'},
                               'transport': {'mu': '1e-05', 'Pr': '0.7'}}}

    elif svp['flparams']['solver_type'] == 'bf':
        tppdict = {'thermoType': {'type': 'heRhoThermo', 'mixture': 'pureMixture',
                   'transport': 'const', 'thermo': 'eConst', 'equationOfState': 'Boussinesq',
                   'specie': 'specie', 'energy': 'sensibleInternalEnergy'},
                   'mixture': {'specie': {'molWeight': '28.96'},
                   'equationOfState': {'rho0': '1', 'T0': '300', 'beta': '3e-03'},
                   'thermodynamics': {'Cv': '772', 'Hf': '0'},
                   'transport': {'mu': '1e-05', 'Pr': '0.7'}}}

    return write_fvdict(htext, tppdict)


def fvmtwrite(node, features):
    mtdict = {}
    mtdict['RAS'] = {}
    mtdict['simulationType'] = 'RAS'
    mtdict['RAS']['model'] = 'kEpsilon'
    mtdict['RAS']['turbulence'] = 'on'
    mtdict['RAS']['printCoeffs'] = 'on'
    htext = ofheader + write_ffile('dictionary', 'constant', 'momentumTransport')
    return write_fvdict(htext, mtdict)


def fvtpwrite():
    htext = ofheader + write_ffile('dictionary', 'constant', 'physicalProperties')
    tpdict = {'transportModel': 'Newtonian', 'nu': '[0 2 -1 0 0 0 0] 1.5e-05'}
    return write_fvdict(htext, tpdict)


def fvrpwrite(node):
    htext = ofheader + write_ffile('dictionary', 'constant', 'radiationProperties')
    raddict = {'0': {'radiation': 'on', 'radiationModel': 'P1',
                'solverFreq': '1', 'absorptionEmissionModel': 'constant',
                'constantCoeffs': {'absorptivity': '0.5', 'emissivity': '0.5', 'E': '0'},
                'scatterModel': 'none', 'sootModel': 'none'},
                '1': {'radiation': 'on', 'radiationModel': 'fvDOM',
                'fvDOMCoeffs': {'nPhi': '3', 'nTheta': '5', 'tolerance': '1e-3', 'maxIter': '10'},
                'solverFreq':'10', 'absorptionEmissionModel':'constant',
                'constantCoeffs': {'absorptivity': '0.5', 'emissivity': '0.5', 'E': '0'},
                'scatterModel': 'none', 'sootModel': 'none', 'useSolarLoad': ('false', 'true')[node.radiation and node.radmodel == '1' and node.solar]}}

    if node.solar:
        raddict[node.radmodel]['SolarLoadCoeffs'] = {'sunDirectionModel': 'sunDirConstant', 'sunDirection': '({0[0]} {0[1]} {0[2]})'.format(-bpy.data.objects[node.sun].location),
                                                     'gridUp': '(0 0 1)', 'gridEast': '(1 0 0)', 'sunLoadModel': 'sunLoadFairWeatherConditions', 'skyCloudCoverFraction': '0',
                                                     'A': '1088', 'B': '0.205', 'groundReflectivity': '0.2', 'C': '0.134', 'useReflectedRays': 'true', 'reflecting': {'nPhi': '10', 'nTheta': '10'},
                                                     'absorptionEmissionModel': 'none', 'scatterModel': 'none', 'sootModel': 'none'}

    return (write_fvdict(htext, raddict[node.radmodel]))

def fvmodwrite(node):
    htext = ofheader + write_ffile('dictionary', 'constant', 'fvModels')
    moddict = {'0': {'radiation': {'type': 'radiation', 'libs': '("libradiationModels.so")'}}}
    return write_fvdict(htext, moddict['0'])


def fvschwrite(svp, node):
    scdict = {}
    scdict['ddtSchemes'] = {'default': 'steadyState'}
    scdict['gradSchemes'] = {}
    scdict['divSchemes'] = {}
    scdict['divSchemes']['default'] = 'none'
    scdict['laplacianSchemes'] = {}
    scdict['interpolationSchemes'] = {}
    scdict['snGradSchemes'] = {}

    if svp['flparams']['solver_type'] == 'sf':
        scdict['ddtSchemes'] = {'default': 'steadyState'}
        scdict['gradSchemes']['default'] = 'Gauss linear'
        scdict['gradSchemes']['limited'] = 'cellLimited Gauss linear 1'
        scdict['gradSchemes']['grad(U)'] = '$limited'
        scdict['gradSchemes']['grad(k)'] = '$limited'
        scdict['gradSchemes']['grad(epsilon)'] = '$limited'
        scdict['divSchemes']['div(phi,U)'] = 'bounded Gauss linearUpwind limited'
        scdict['divSchemes']['turbulence'] = 'bounded Gauss limitedLinear 1'
        scdict['divSchemes']['div(phi,k)'] = '$turbulence'
        scdict['divSchemes']['div(phi,epsilon)'] = '$turbulence'
        scdict['divSchemes']['div((nuEff*dev2(T(grad(U)))))'] = 'Gauss linear'
        scdict['laplacianSchemes']['default'] = 'Gauss linear corrected'
        scdict['interpolationSchemes']['default'] = 'linear'
        scdict['snGradSchemes']['default'] = 'corrected'

    elif svp['flparams']['solver_type'] == 'bf':
        scdict['ddtSchemes'] = {'default': 'steadyState'}
        scdict['gradSchemes']['default'] = 'Gauss linear'
        scdict['divSchemes']['div(phi,U)'] = 'bounded Gauss upwind'
        scdict['divSchemes']['div(phi,e)'] = 'bounded Gauss upwind'
        scdict['divSchemes']['div(phi,K)'] = 'Gauss linear'
        scdict['divSchemes']['div(phi,(p|rho))'] = 'Gauss linear'
        scdict['divSchemes']['div(phi,k)'] = 'bounded Gauss upwind'
        scdict['divSchemes']['div(phi,epsilon)'] = 'bounded Gauss upwind'
        scdict['divSchemes']['div(((rho*nuEff)*dev2(T(grad(U)))))'] = 'Gauss linear'

        if svp['flparams']['features']['rad'] and svp['flparams']['radmodel'] == '1':
            scdict['divSchemes']['div(Ji,Ii_h)'] = 'bounded Gauss linearUpwind grad(Ii_h)'

        scdict['laplacianSchemes']['default'] = 'Gauss linear orthogonal'
        scdict['snGradSchemes']['default'] = 'orthogonal'
        scdict['interpolationSchemes']['default'] = 'linear'

    if node.age:
        scdict['divSchemes']['div(phi,age)'] = 'bounded Gauss upwind'

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


def oftomesh(ofb, vl, fomats, st, ns, nf, bo):
    vcoords = []
    findices = []
    fi = []
    fn, p, f = 0, 0, 0
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
                        
                    except Exception:
                        pass

                if not fn and fi:
                    f += 1

                    if f >= ns[0] or not bo:
                        findices.append(fi)

                    fi = []
                    fn = 0
                                       
                elif '(' in line and ')' in line:
                    if f >= ns[0] or not bo:
                        findices.append([int(x) for x in line.split('(')[1].split(')')[0].split()])

                    f += 1

                else:
                    try:
                        if prevline == '\n' and int(line) < 100:
                            fn = int(line)
                    except Exception:
                        fn = 0
            
            prevline = line

    mesh = bpy.data.meshes.new("mesh")
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
            if not bo:
                if face.index >= n and face.index <= n + nf[ni]:
                    face.material_index = ni
                    mi = 1
            else:
                if face.index + 1 > n - ns[0] and face.index < sum(nf[:ni + 1]):
                    face.material_index = ni
                    mi = 1

        if not mi and not bo:
            face.material_index = len(ns)


# def ret_fvbnutilda_menu(mat, context):
#     if context.scene.vi_params.get('flparams') and context.scene.vi_params['flparams'].get('scenario') and context.scene.vi_params['flparams']['scenario'] != '4':
#         # print(mat.name, flovi_nut_dict[context.scene.vi_params['flparams']['scenario']][mat.flovi_bmb_type].keys())
#         return [('{}'.format(b), '{}'.format(b), '{} boundary type'.format(b)) for b in flovi_nut_dict[context.scene.vi_params['flparams']['scenario']][mat.flovi_bmb_type]]
#     else:
#         return [('{}'.format(b), '{}'.format(b), '{} boundary type'.format(b)) for b in flovi_nutilda_bounds[context.scene.vi_params['flparams']['solver_type']][mat.flovi_bmb_type]]

# def fvidwrite(node):
#     htext = ofheader + write_ffile('ascii', 'volScalerField', 'IDefault')
#     iddict = {'0': 'dimensions': '[1 0 -3 0 0 0 0]', 'internalField': 'uniform 0', 'boundaryField': {'".*"': {'type': 'greyDiffusiveRadiation', 'emissivityMode': 'lookup', 'emissivity', 'uniform 1.0', 'value': 'uniform 0'}}}
#     return write_fvdict(htext, iddict['0'])