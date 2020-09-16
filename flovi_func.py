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

flovi_p_bounds = {'if': {'0': ('zeroGradient', 'fixedValue', 'totalPressure'), '1': ('zeroGradient', 'fixedValue'), '2': ['None']}, 
                'sf': {'0': ('zeroGradient', 'fixedValue', 'freestreamPressure', 'totalPressure'), '1': ['zeroGradient'], '2': ['None']},
                'bsf': {'0': ('calculated', 'calculated', 'calculated', 'calculated'), '1': ('calculated', 'calculated'), '2': ('None',)},
                'bbsf': {'0': ('calculated', 'calculated', 'calculated', 'calculated'), '1': ('calculated', 'calculated'), '2': ('None',)}}

flovi_u_bounds = {'icoFoam': {'0': ('zeroGradient','noSlip', 'fixedValue'), '1': ('zeroGradient', 'noSlip', 'fixedValue'), '2': ['None']}, 
                'simpleFoam': {'0': ('zeroGradient', 'fixedValue', 'inletOutlet', 'freestream', 'pressureInletOutletVelocity', 'slip'), '1': ('noSlip', 'fixedValue', 'slip'), '2': ['None']},
                'buoyantSimpleFoam': {'0': ('zeroGradient', 'fixedValue', 'inletOutlet', 'freestream', 'pressureInletOutletVelocity', 'slip'), '1': ('noSlip', 'fixedValue', 'slip'), '2': ['None']}}

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

flovi_prgh_bounds = {'buoyantSimpleFoam': {'0': ('fixedFluxPressure', 'prghTotalHydrostaticPressure', 'fixedValue', 'prghPressure'), '1': ('fixedFluxPressure', 'fixedValue'), '2': ['None']}}

flovi_a_bounds = {'buoyantSimpleFoam': {'0': ('calculated',), '1': ('compressible::alphatJayatillekeWallFunction',), '2': ['None']}}
flovi_rad_bounds = {'buoyantSimpleFoam': {'0': ('MarshakRadiation',), '1': ('MarshakRadiation',), '2': ['None']}}

#flovi_p_dimens


#ico_p_bounds = {'p': ('zeroGradient', 'fixedValue'), 
#              'U': ('zeroGradient','noSlip', 'fixedValue')}
#ico_w_bounds = {'p': ('zeroGradient', 'fixedValue'), 
#              'U': ('zeroGradient', 'noSlip', 'fixedValue')}
#sim_p_bounds = {'p': ('zeroGradient', 'fixedValue', 'freestreamPressure'), 
#              'U': ('zeroGradient', 'fixedValue', 'inletOutlet', 'freestream', 'pressureInletOutletVelocity', 'slip'),
#              'nut':('nutkWallFunction', 'calculated'),
#              'nuTilda': ('zeroGradient', 'fixedValue'),
#              'k': ('fixedValue', 'inletOutlet'),
#              'epsilon': ('fixedValue', 'inletOutlet'),
#              'omega': ('zeroGradient', 'fixedValue')}
#sim_w_bounds = {'p': ['zeroGradient'], 
#              'U': ('noSlip', 'fixedValue', 'slip'),
#              'nut':('nutkWallFunction'),
#              'nuTilda': ('zeroGradient', 'fixedValue'),
#              'k': ('kqRWallFunction'),
#              'epsilon': ('epsilonWallFunction'),
#              'omega': ('omegaWallFunction')}
#bsim_p_bounds = {'T': ('zeroGradient', 'fixedValue', 'inletOutlet'),
#                'p_rgh': ('fixedFluxPressure', 'prghTotalHydrostaticPressure'),
#                'alphat': ('calculated', 'compressible::alphatWallFunction')}
#bsim_w_bounds = {'T': ('zeroGradient', 'fixedValue'),
#                'p_rgh': ('fixedFluxPressure', 'fixedValue'),
#                'alphat': ('calculated', 'compressible::alphatWallFunction')}
#rbsimbounddict = {'G': 'MarshakRadiation'}
#fvrbsimbounddict = {'IDefault': ('greyDiffusiveRadiation', 'calculated')}

#bound_dict = {'icoFoam': (ico_p_bounds, ico_w_bounds), 'simpleFoam': (sim_p_bounds, sim_w_bounds), 'boussinesc': (bsim_p_bounds, bsim_w_bounds)}

#def ret_fvb_menu(mat, context):
#    return [('{}'.format(b), '{}'.format(b), '{} boundary type'.format(b)) for b in bound_dict[context.scene['flparams']['solver']][int(mat.flovi_bmb_type)]]
# 

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
    location = 'location    "{}";\n'.format(loc) if loc else ''
    return 'FoamFile\n  {{\n    version   2.0;\n    format    ascii;\n    {}    class     {};\n    object    {};\n  }}\n\n'.format(location, cla, obj)

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

def flovi_bm_update(self, context):
    scene = context.scene
    svp = scene.vi_params
    ovp = context.object.vi_params
    svp = scene.vi_params
#    svp['flparams']['solver'] = ovp.flovi_solver
#    svp['flparams']['turbulence'] = ovp.flovi_turb

#@writeheader    
def fvbmwrite(o, expnode):
    bm = bmesh.new()
    tempmesh = o.to_mesh(scene = bpy.context.scene, apply_modifiers = True, settings = 'PREVIEW')
    bm.from_mesh(tempmesh)
    bm.verts.ensure_lookup_table()
    bm.transform(o.matrix_world)
    bpy.data.meshes.remove(tempmesh)   
    [xs, ys, zs] = [[v.co[i] for v in bm.verts] for i in range(3)]
    bm.transform(mathutils.Matrix.Translation((-min(xs), -min(ys), -min(zs))))
    o['flovi_translate'] = (-min(xs), -min(ys), -min(zs))
    lengths = [mathutils.Vector(v.co).length for v in bm.verts]
    vert0 = bm.verts[lengths.index(min(lengths))]
    angles = [mathutils.Vector(v.co).angle(mathutils.Vector((0, 0, 1))) for v in bm.verts if v != vert0]
#    vert0 = [v for v in bm.verts if v.co[:] == (min(xs), min(ys), min(zs))][0]
    vert4 = bm.verts[angles.index(min(angles)) + 1]
#    print(vert0.index, vert4.index)
#    vert4 = [v for v in bm.verts if (v.co[0], v.co[1]) == (vert0.co[0], vert0.co[1]) and v.co[2] != vert0.co[2]][0]

    for face in bm.faces:
        if vert0 in face.verts and vert4 not in face.verts:
            vis = [vert.index for vert in face.verts][::-1]
            vertorder1 = vis[vis.index(vert0.index):] + vis[:vis.index(vert0.index)]
#            vertorder1 = [vertorder1[0], vertorder1[3], vertorder1[2], vertorder1[1]]
        if vert4 in face.verts and vert0 not in face.verts:
            vis = [vert.index for vert in face.verts]
            vertorder2 = vis[vis.index(vert4.index):] + vis[:vis.index(vert4.index)]
        
    vertorder = ''.join(['{} '.format(v) for v in vertorder1 + vertorder2])

#    omw, bmovs = o.matrix_world, [vert for vert in o.data.vertices]
#    xvec, yvec, zvec = (omw*bmovs[3].co - omw*bmovs[0].co).normalized(), (omw*bmovs[2].co - omw*bmovs[3].co).normalized(), (omw*bmovs[4].co - omw*bmovs[0].co).normalized() 
#    ofvpos = [[(omw*bmov.co - omw*bmovs[0].co)*vec for vec in (xvec, yvec, zvec)] for bmov in bmovs]
#    bmdict = "vertices\n(\n" + "\n".join(["  ({0:.3f} {1:.3f} {2:.3f})" .format(*ofvpo) for ofvpo in ofvpos]) +"\n);\n\n"
    bmdict = "vertices\n(\n" + "\n".join(["  ({0[0]:.3f} {0[1]:.3f} {0[2]:.3f})" .format(v.co) for v in bm.verts]) +"\n);\n\n"
    bmdict += "blocks\n(\n  hex ({}) ({} {} {}) simpleGrading ({} {} {})\n);\n\n".format(vertorder, expnode.bm_xres, expnode.bm_yres, expnode.bm_zres, expnode.bm_xgrad, expnode.bm_ygrad, expnode.bm_zgrad) 
    bmdict += "edges\n(\n);\n\nboundary\n(\n"  
    bmdict += fvboundwrite(o)
    bm.free()
    return ofheader + write_ffile('dictionary', '', 'blockMeshDict') + bmdict
    
def fvblbmgen(mats, ffile, vfile, bfile, meshtype):
    scene = bpy.context.scene
    matfacedict = {mat.name:[0, 0] for mat in mats}
    
    for line in bfile.readlines():
        if line.strip() in matfacedict:
            mat = line.strip()
        elif line.strip() in [o.name for o in bpy.data.objects]:
            mat = bpy.data.objects[line.strip()].data.materials[0].name
        if 'nFaces' in line:
            matfacedict[mat][1] = int(line.split()[1].strip(';'))
        if 'startFace' in line:
            matfacedict[mat][0] = int(line.split()[1].strip(';'))
    bobs = [ob for ob in scene.objects if ob.get('VIType') and ob['VIType'] == 'FloViMesh']
    
    if bobs:
        o = bobs[0]
        selobj(scene, o)
        while o.data.materials:
            bpy.ops.object.material_slot_remove()
    else:
        bpy.ops.object.add(type='MESH', layers=(False, False, False, False, False, False, 
                                                False, False, False, False, False, False, 
                                                False, False, False, False, False, False, False, True))
        o = bpy.context.object
        o['VIType'] = 'FloViMesh'
    
    o.name = meshtype
    for mat in mats:
        if mat.name not in o.data.materials:
            bpy.ops.object.material_slot_add()
            o.material_slots[-1].material = mat 
    
    matnamedict = {mat.name: m for  m, mat in enumerate(o.data.materials)}    
    bm = bmesh.new()

    for line in [line for line in vfile.readlines() if line[0] == '(' and len(line.split(' ')) == 3]:
        bm.verts.new([float(vpos) for vpos in line[1:-2].split(' ')])

    if hasattr(bm.verts, "ensure_lookup_table"):
        bm.verts.ensure_lookup_table()

    for l, line in enumerate([line for line in ffile.readlines() if '(' in line and line[0].isdigit() and len(line.split(' ')) == int(line[0])]):
        newf = bm.faces.new([bm.verts[int(fv)] for fv in line[2:-2].split(' ')])
        
        for facerange in matfacedict.items():
            if l in range(facerange[1][0], facerange[1][0] + facerange[1][1]):
                newf.material_index = matnamedict[facerange[0]]

    bm.transform(o.matrix_world.inverted())
    bm.to_mesh(o.data)
    bm.free()

#def fvbmr(scene, o):
#    points = ''.join(['({} {} {})\n'.format(o.matrix_world * v.co) for v in o.data.verts]) + ')'
#    with open(os.path.join(scene['flparams']['ofcpfilebase'], 'points'), 'r') as pfile:
#        pfile.write(ofheader + write_ffile('vectorField', '"constant/polyMesh"', 'points') + points)
#    faces = ''.join(['({} {} {} {})\n'.format(f.vertices) for f in o.data.faces]) + ')'
#    with open(os.path.join(scene['flparams']['ofcpfilebase'], 'faces'), 'r') as ffile:
#        ffile.write(ofheader + write_ffile('vectorField', '"constant/polyMesh"', 'points') + faces)


def fvmat(self, mn, bound):
#    fvname = on.replace(" ", "_") + self.name.replace(" ", "_") 
    begin = '\n  {}\n  {{\n    type    '.format(mn)  
    end = ';\n  }\n'
    
    if bound == 'p':
        val = 'uniform {}'.format(self.flovi_bmbp_val) if not self.flovi_p_field else '$internalField'
        pdict = {'0': self.flovi_bmbp_subtype, '1': self.flovi_bmbp_subtype, '2': 'symmetry', '3': 'empty'}
        ptdict = {'zeroGradient': 'zeroGradient', 'fixedValue': 'fixedValue;\n    value    {}'.format(val), 
                'calculated': 'calculated;\n    value    $internalField', 
                'freestreamPressure': 'freestreamPressure', 
                'totalPressure': 'totalPressure;\n    p0      uniform {};\n    gamma    {};\n    value    {}'.format(self.flovi_bmbp_p0val, self.flovi_bmbp_gamma, val), 'symmetry': 'symmetry', 'empty': 'empty'}
#        if pdict[self.flovi_bmb_type] == 'zeroGradient':
        entry = ptdict[pdict[self.flovi_bmb_type]]            
#        return begin + entry + end 
    
    elif bound == 'U':
        val = 'uniform ({} {} {})'.format(*self.flovi_bmbu_val) if not self.flovi_u_field else '$internalField'
        Udict = {'0': self.flovi_bmbu_subtype, '1': self.flovi_bmbu_subtype, '2': 'symmetry', '3': 'empty'}
        Utdict = {'fixedValue': 'fixedValue;\n    value    {}'.format(val), 'slip': 'slip', 'noSlip': 'noSlip', 'inletOutlet': 'inletOutlet;\n    inletValue    $internalField;\n    value    $internalField',
                  'pressureInletOutletVelocity': 'pressureInletOutletVelocity;\n    value    $internalField', 'zeroGradient': 'zeroGradient', 'symmetry': 'symmetry', 
                  'freestream': 'freestream;\n    freestreamValue    $internalField','calculated': 'calculated;\n    value    $internalField', 'empty': 'empty'}
        entry = Utdict[Udict[self.flovi_bmb_type]]            
#        return begin + entry + end
        
    elif bound == 'nut':
        ndict = {'0': self.flovi_bmbnut_subtype, '1': self.flovi_bmbnut_subtype, '2': 'symmetry', '3': 'empty'}
        ntdict = {'nutkWallFunction': 'nutkWallFunction;\n    value    $internalField', 'nutUSpaldingWallFunction': 'nutUSpaldingWallFunction;\n    value    $internalField', 
        'calculated': 'calculated;\n    value    $internalField', 'inletOutlet': 'inletOutlet;\n    inletValue    $internalField\n    value    $internalField',  'symmetry': 'symmetry','empty': 'empty'}
        entry = ntdict[ndict[self.flovi_bmb_type]]            
#        return begin + entry + end

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
#        return begin + entry + end

    elif bound == 't':
        val = 'uniform {}'.format(self.flovi_bmbt_val) if not self.flovi_t_field else '$internalField'
        ival = 'uniform {}'.format(self.flovi_bmbti_val) if not self.flovi_t_field else '$internalField'
        tdict = {'0': self.flovi_bmbt_subtype, '1': self.flovi_bmbt_subtype, '2': 'symmetry', '3': 'empty'}
        ttdict = {'zeroGradient': 'zeroGradient', 'fixedValue': 'fixedValue;\n    value    {}'.format(val), 'inletOutlet': 'inletOutlet;\n    inletValue    {}\n    value    {}'.format(ival, val),
        'calculated': 'calculated;\n    value    $internalField', 'symmetry': 'symmetry', 'empty': 'empty'}
        entry = ttdict[tdict[self.flovi_bmb_type]]  
        
    elif bound == 'p_rgh':
        val = 'value {:.4f}'.format(self.flovi_prgh_val) if not self.flovi_prgh_field else '$internalField'
        prghdict = {'0': self.flovi_prgh_subtype, '1': self.flovi_prgh_subtype, '2': 'symmetry', '3': 'empty'}
        prghtdict = {'fixedFluxPressure': 'fixedFluxPressure;\n    value    {}'.format(val), 'fixedValue': 'fixedValue;\n    value    {}'.format(val), 
                     'prghTotalHydrostaticPressure': 'prghTotalHydrostaticPressure;\n    p0              $internalField;\n     gamma           1;\n     value    {}'.format(val), 
                     'fixedValue': 'fixedValue;\n    value    {}'.format(val),
                     'calculated': 'calculated;\n    value    $internalField', 
                     'prghPressure': 'prghPressure;\n    p    $internalField;\n    value    $internalField', 
                     'symmetry': 'symmetry', 'empty': 'empty'}
        entry = prghtdict[prghdict[self.flovi_bmb_type]] 

    elif bound == 'a':
        val = 'uniform {}'.format(self.flovi_a_val) if not self.flovi_a_field else '$internalField'
        tdict = {'0': self.flovi_a_subtype, '1': self.flovi_a_subtype, '2': 'symmetry', '3': 'empty'}
        ttdict = {'compressible::alphatJayatillekeWallFunction': 'compressible::alphatJayatillekeWallFunction;\n    Prt    0.85;\n    value           $internalField;', 
                  'fixedValue': 'fixedValue;\n    value    {}'.format(val), 
                  'inletOutlet': 'inletOutlet;\n    inletValue    $internalField\n    value    $internalField',
                  'calculated': 'calculated;\n    value    $internalField', 'symmetry': 'symmetry', 'empty': 'empty'}
        entry = ttdict[tdict[self.flovi_bmb_type]] 
        
    elif bound == 'e':
        edict = {'0': self.flovi_bmbe_subtype, '1': self.flovi_bmbe_subtype, '2': 'symmetry', '3': 'empty'}
        etdict = {'symmetry': 'symmetry', 'empty': 'empty', 'inletOutlet': 'inletOutlet;\n    inletValue    $internalField;\n    value    $internalField', 'fixedValue': 'fixedValue;\n    value    $internalField', 
                  'epsilonWallFunction': 'epsilonWallFunction;\n    value    $internalField', 'calculated': 'calculated;\n    value    $internalField', 'symmetry': 'symmetry', 'empty': 'empty'}
        entry = etdict[edict[self.flovi_bmb_type]]            
#        return begin + entry + end
        
    elif bound == 'o':
        odict = {'0': self.flovi_bmbo_subtype, '1': self.flovi_bmbo_subtype, '2': 'symmetry', '3': 'empty'}
        otdict = {'symmetry': 'symmetry', 'empty': 'empty', 'inletOutlet': 'inletOutlet;\n    inletValue    $internalField\n    value    $internalField', 'zeroGradient': 'zeroGradient', 
                  'omegaWallFunction': 'omegaWallFunction;\n    value    $internalField', 'fixedValue': 'fixedValue;\n    value    $internalField'}
        entry = otdict[odict[self.flovi_bmb_type]]            
#        return begin + entry + end
        
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
    if not node.buoyancy:# or (node.buoyancy and not node.buossinesq):
        pentry = "dimensions [{} {} {} {} 0 0 0];\ninternalField   uniform {};\n\nboundaryField\n{{\n".format('0', '2', '-2', '0', '{}'.format(node.pnormval))
    else:
        pentry = "dimensions [{} {} {} {} 0 0 0];\ninternalField   uniform {};\n\nboundaryField\n{{\n".format('1', '-1', '-2', '0', '{}'.format(node.pabsval))
        
    (Uentry, nutildaentry, nutentry, kentry, eentry, oentry, tentry, p_rghentry, aentry, Gentry) = ["dimensions [{} {} {} {} 0 0 0];\ninternalField   uniform {};\n\nboundaryField\n{{\n".format(*var) for var in ( 
                                                                                ('0', '1', '-1', '0', '({:.4f} {:.4f} {:.4f})'.format(*node.uval)), 
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
        ovp = o.vi_params
        
        for mat in o.data.materials: 
            mvp = mat.vi_params
            matname = '{}_{}'.format(o.name, mat.name)

            if mvp.mattype == '2':
                pentry += mvp.flovi_mat(matname, 'p')
                Uentry += mvp.flovi_mat(matname, 'U')
                if node.solver != 'icoFoam':
                    if node.turbulence != 'laminar':
                        nutentry += mvp.flovi_mat(matname, 'nut')                    
                        if node.turbulence ==  'SpalartAllmaras':
                            nutildaentry += mvp.flovi_mat(matname, 'nutilda')
                        elif node.turbulence ==  'kEpsilon':
                            kentry += mvp.flovi_mat(matname, 'k')
                            eentry += mvp.flovi_mat(matname, 'e')
                        elif node.turbulence ==  'kOmega':
                            kentry += mvp.flovi_mat(matname, 'k')
                            oentry += mvp.flovi_mat(matname, 'o')
                    if node.buoyancy:
                        tentry += mvp.flovi_mat(matname, 't')
                        p_rghentry += mvp.flovi_mat(matname, 'p_rgh')
                        aentry += mvp.flovi_mat(matname, 'a')
                        if node.radiation:
                            Gentry += mvp.flovi_mat(matname, 'G')

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
    
    with open(os.path.join(svp['flparams']['of0filebase'], 'p'), 'w') as pfile:
        pfile.write(ofheader + write_ffile('volScalarField', '', 'p') + pentry)
    with open(os.path.join(svp['flparams']['of0filebase'], 'U'), 'w') as Ufile:
        Ufile.write(ofheader + write_ffile('volVectorField', '', 'U') + Uentry)
        
    if node.solver != 'icoFoam':
        with open(os.path.join(svp['flparams']['of0filebase'], 'nut'), 'w') as nutfile:
            nutfile.write(ofheader + write_ffile('volScalarField', '', 'nut') + nutentry)
        if node.turbulence == 'SpalartAllmaras':
            with open(os.path.join(svp['flparams']['of0filebase'], 'nuTilda'), 'w') as nutildafile:
                nutildafile.write(ofheader + write_ffile('volScalarField', '', 'nut') + nutildaentry)
        if node.turbulence == 'kEpsilon':
            with open(os.path.join(svp['flparams']['of0filebase'], 'k'), 'w') as kfile:
                kfile.write(ofheader + write_ffile('volScalarField', '', 'k') + kentry)
            with open(os.path.join(svp['flparams']['of0filebase'], 'epsilon'), 'w') as efile:
                efile.write(ofheader + write_ffile('volScalarField', '', 'epsilon') + eentry)
        if node.turbulence == 'kOmega':
            with open(os.path.join(svp['flparams']['of0filebase'], 'k'), 'w') as kfile:
                kfile.write(ofheader + write_ffile('volScalarField', '', 'k') + kentry)
            with open(os.path.join(svp['flparams']['of0filebase'], 'omega'), 'w') as ofile:
                ofile.write(ofheader + write_ffile('volScalarField', '', 'omega') + oentry)
        if node.buoyancy:
            with open(os.path.join(svp['flparams']['of0filebase'], 'T'), 'w') as tfile:
                tfile.write(ofheader + write_ffile('volScalarField', '', 'T') + tentry)
            with open(os.path.join(svp['flparams']['of0filebase'], 'alphat'), 'w') as afile:
                afile.write(ofheader + write_ffile('volScalarField', '', 'alphat') + aentry)
            with open(os.path.join(svp['flparams']['of0filebase'], 'p_rgh'), 'w') as prghfile:
                prghfile.write(ofheader + write_ffile('volScalarField', '', 'p_rgh') + p_rghentry)  
            if node.radiation:
                with open(os.path.join(svp['flparams']['of0filebase'], 'G'), 'w') as Gfile:
                    Gfile.write(ofheader + write_ffile('volScalarField', '', 'G') + Gentry) 
#            with open(os.path.join(svp['flparams']['of0filebase'], 'g'), 'w') as gfile:
#                gfile.write(fvgwrite())

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
    
def fvcdwrite(solver, st, dt, et):
    pw = 0 if solver == 'icoFoam' else 1
    ps = []
    

    for o in bpy.data.objects:
        if o.type == 'MESH' and o.vi_params.vi_type == '2':
            dom = o
        if o.type == 'EMPTY' and o.vi_params.flovi_probe:
            ps.append(o)
    if ps:
        bpy.context.scene.vi_params['flparams']['probes'] = [p.name for p in ps]
        probe_vars = 'p U T'          
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
}}'''.format(probe_vars, ''.join(['( {0[0]} {0[1]} {0[2]} )\n'.format(p.location) for p in ps]), ','.join(['{}'.format(p.name) for p in ps]))

    else:
        probe_text = ''
        bpy.context.scene.vi_params['flparams']['probes'] = []
    return 'FoamFile\n{\n  version     2.0;\n  format      ascii;\n  class       dictionary;\n  location    "system";\n  object      controlDict;\n}\n\n' + \
            'application     {};\nstartFrom       startTime;\nstartTime       {};\nstopAt          endTime;\nendTime         {};\n'.format(solver, st, et, dt)+\
            'deltaT          {};\nwriteControl    timeStep;\nwriteInterval   {};\npurgeWrite      {};\nwriteFormat     ascii;\nwritePrecision  6;\n'.format(dt, 1, pw)+\
            'writeCompression off;\ntimeFormat      general;\ntimePrecision   6;\nrunTimeModifiable true;\n\n' + probe_text



def fvsolwrite(node, solver):
    basedict = {'solvers': {}, }
 
    if node.transience == '0' and node.turbulence != 'laminar':
        if not node.buoyancy:
            soldict = {'solvers': {'p': {'solver': 'GAMG', 'smoother': 'GaussSeidel', 'tolerance': '1e-6', 'relTol': '0.1'}},
                       'SIMPLE': {'nNonOrthogonalCorrectors': '0', 'pRefCell': '0', 'pRefValue': '0', 
                               'residualControl': {'p': '1e-4', 'U': '1e-4'}},
                       'potentialFlow': {'nNonOrthogonalCorrectors': '10'},
                       'relaxationFactors': {'fields': {'p': '0.3'}, 'equations': {'U': '0.7'}}}
            if node.turbulence == 'kEpsilon':
                soldict['solvers']['"(U|k|epsilon)"'] = {'solver': 'smoothSolver', 'smoother': 'symGaussSeidel', 'tolerance': '1e-6', 'relTol': '0.1'}
                soldict['SIMPLE']['residualControl']['"k|epsilon"'] = '1e-3' 
                soldict['relaxationFactors']['equations']['"k|epsilon"'] = '1e-3' 
            if node.turbulence == 'kOmega':
                soldict['solvers']['"(U|k|omega)"'] = {'solver': 'smoothSolver', 'smoother': 'symGaussSeidel', 'tolerance': '1e-6', 'relTol': '0.1'}
                soldict['SIMPLE']['residualControl']['"k|omega"'] = '1e-3'
                soldict['relaxationFactors']['equations']['"k|omega"'] = '1e-3' 
            if node.turbulence == 'SpalartAllmaras': 
                soldict['solvers']['U'] = {'solver': 'smoothSolver', 'smoother': 'GaussSeidel', 'nSweeps': '2', 'tolerance': '1e-08', 'relTol': '0.1'}
                soldict['solvers']['nuTilda'] = {'solver': 'smoothSolver', 'smoother': 'GaussSeidel', 'nSweeps': '2', 'tolerance': '1e-6', 'relTol': '0.1'}
                soldict['SIMPLE']['residualControl']['nuTilda'] = '1e-3'
                soldict['relaxationFactors']['equations']['nuTilda'] = '1e-3' 
        else:
            soldict = {'solvers': {}, 
                    'SIMPLE': {'nNonOrthogonalCorrectors': '0', 'pRefCell': '0', 'pRefValue': '{}'.format(node.pabsval), 
                               'residualControl': {'U': '1e-3'}},
                    'relaxationFactors': {'fields': {}, 'equations': {'U': '0.2'}}}
           
            if node.buossinesq:
                soldict['solvers']['p_rgh'] = {'solver': 'PCG', 'preconditioner': 'DIC', 'tolerance': '1e-06', 'relTol': '0.01'}
                soldict['solvers']['"(U|e)"'] = {'solver': 'PBiCGStab', 'preconditioner': 'DILU', 'tolerance': '1e-05', 'relTol': '0.1'}
                soldict['SIMPLE']['residualControl']['p_rgh'] = '1e-3'
                soldict['SIMPLE']['residualControl']['e'] = '1e-3'
                soldict['relaxationFactors']['fields']['p_rgh'] = '0.7'
                soldict['relaxationFactors']['equations']['e'] = '0.2'  
                
                
                if node.turbulence == 'kEpsilon':
                    soldict['solvers']['"(k|epsilon)"'] = {'solver': 'PBiCGStab', 'preconditioner': 'DILU', 'tolerance': '1e-05', 'relTol': '0.1'}
                    soldict['SIMPLE']['residualControl']['"(k|epsilon)"'] = '1e-3'
                    soldict['relaxationFactors']['equations']['"(k|epsilon)"'] = '0.7'
                    
                elif node.turbulence == 'kOmega':
                    soldict['solvers']['"(k|omega)"'] = {'solver': 'PBiCGStab', 'preconditioner': 'DILU', 'tolerance': '1e-05', 'relTol': '0.1'}
                    soldict['SIMPLE']['residualControl']['"(k|omega)"'] = '1e-3'
                elif node.turbulence == 'SpalartAllmaras':
                    soldict['solvers']['nuTilda'] = {'solver': 'PBiCGStab', 'preconditioner': 'DILU', 'tolerance': '1e-05', 'relTol': '0.1'}
                    soldict['SIMPLE']['residualControl']['nuTilda'] = '1e-3'
            else:
                soldict['solvers']['"(U|h|k|epsilon|omega|nutilda)"'] = {'solver': 'PBiCGStab', 'preconditioner': 'DILU', 'tolerance': '1e-05', 'relTol': '0.1'}
                soldict['SIMPLE']['residualControl']['h'] = '1e-3'
                soldict['relaxationFactors']['fields']['rho'] = '1'            
            if node.radiation:
                soldict['solvers']['"G.*"'] = {'$p_rgh': '', 'tolerance': '1e-05', 'relTol': '0.1'}
                
                if node.radmodel == '0':
                    soldict['SIMPLE']['residualControl']['G'] = '1e-3'
                    soldict['relaxationFactors']['equations']['nuTilda'] = '1e-3' 

                elif node.radmodel == '1':
                    soldict['solvers']['"G.*"'] = {'$p_rgh': '', 'tolerance': '1e-05', 'relTol': '0.1'}

    htext = ofheader + write_ffile('dictionary', 'system', 'fvSolution') 
    ntext = write_fvdict(htext, soldict)
    
    # for d in soldict:
    #     if isinstance(soldict[d], str):
    #         ntext += '{} {};\n\n'.format(d, soldict[d]) 
    #     elif isinstance(soldict[d], dict):
    #         ntext += '{}\n{{\n'.format(d)
    #         for sd in soldict[d]:
    #             if isinstance(soldict[d][sd], str):
    #                 ntext += '  {} {};\n'.format(sd, soldict[d][sd]) 
    #             elif isinstance(soldict[d][sd], dict):
    #                 ntext += '  {}\n  {{\n'.format(sd)
    #                 for ssd in soldict[d][sd]:
    #                     if isinstance(soldict[d][sd][ssd], str):
    #                         ntext += '    {} {};\n'.format(ssd, soldict[d][sd][ssd])
    #                 ntext += '  }\n'
    #         ntext += '}\n'
    print(ntext)
        
    solvers = {'sf': {'names': ('solver', 'smoother', 'tolerance', 'relTol'), 'p': ('GAMG', 'GaussSeidel', '1e-6', '0.1'), '"(U|k|omega|epsilon)"': ('smoothSolver', 'symGaussSeidel', '1e-6', '0.1')},
     'bbsf': {'names': ('solver', 'preconditioner', 'tolerance', 'relTol'), 'p_rgh': ('PCG', 'DIC', '1e-8', '0.01'), '"(U|e|k|omega|epsilon)"': ('PBiCGStab', 'DILU;', '1e-7', '0.1')},
     'bsf':{}}
    

    resids = {'sf': {'names': ('p', 'U', '"(k|epsilon|omega)"'), 'residualControl': (node.presid, node.uresid, node.keoresid)},
              'bbsf': {'names': ('p_rgh', 'U', 'e', '"(k|epsilon|omega)"'), 'residualControl': (node.presid, node.uresid, '1e-3', node.keoresid)},
         'bsf':{}}
    
    pot = {'sf': {'names': ('nNonOrthogonalCorrectors', ), 'potentialFlow': ('10',)},
         'bsf':{}}
    
    relax_f = {'bbsf': {'names': ('p_rgh',), 'fields': ('0.7',)},
              'sf': {'names': ('p', ), 'fields': ('0.3',)},
         'bsf':{}}
    
    relax_e = {'sf': {'names': ('U', '"(k|epsilon)"'), 'equations': ('0.7', '0.7')},
              'bbsf': {'names': ('U', 'e', '"(k|epsilon|omega)"'), 'equations': ('0.2', '0.2', '0.7')},
         'bsf':{}}
    
    text = ofheader + write_ffile('dictionary', 'system', 'fvSolution') + 'solvers\n{\n  '
    
    for sol in solvers[solver]:
        if sol != 'names':
            text += sol + '\n  {\n' + ';\n'.join(['    {} {}'.format(s, solvers[solver][sol][si]) for si, s in enumerate(solvers[solver]['names'])]) + ';\n  }\n'
    text += '}\n\n'
    text = ofheader + write_ffile('dictionary', 'system', 'fvSolution') + 'solvers\n{\n  ' + ntext
    text += 'SIMPLE\n{\n  '
    
    if solver in ('bbsf', 'sf'):
        text += '''nNonOrthogonalCorrectors {};\n  pRefCell    0;\n  pRefValue   0;\n\n  '''.format(('0', '1')[solver == 'bbsf'])
       
    for sol in resids[solver]:
        if sol != 'names':
            text += sol + '\n  {\n' + ';\n'.join(['    {} {}'.format(s, resids[solver][sol][si]) for si, s in enumerate(resids[solver]['names'])]) + ';\n  }\n'
    text += '}\n\n'
    
    if solver == 'sf':
        for sol in pot[solver]:
            if sol != 'names':
                text += sol + '\n  {\n' + ';\n'.join(['    {} {}'.format(s, pot[solver][sol][si]) for si, s in enumerate(pot[solver]['names'])]) + ';\n  }\n'
        text += '\n\n'
    
    text += 'relaxationFactors\n{\n  '
    
    for sol in relax_f[solver]:
        if sol != 'names':
            text += sol + '\n  {\n' + ';\n'.join(['    {} {}'.format(s, relax_f[solver][sol][si]) for si, s in enumerate(relax_f[solver]['names'])]) + ';\n  }\n'
    
    for sol in relax_e[solver]:
        if sol != 'names':
            text += sol + '\n  {\n' + ';\n'.join(['    {} {}'.format(s, relax_e[solver][sol][si]) for si, s in enumerate(relax_e[solver]['names'])]) + ';\n  }\n'
    text += '}'
    return ntext

# def fventry(n, ps, vs)
#     return '  {}\n
def fvtppwrite(node, solver): 
    thermo = {'bbsf': {'names': ('type', 'mixture', 'transport', 'thermo', 'equationOfState', 'specie', 'energy'), 
                       'thermoType': ('heRhoThermo', 'pureMixture', 'const', 'eConst', 'Boussinesq', 'specie', 'sensibleInternalEnergy')},
              'sf': {'names': ('solver', 'preconditioner', 'tolerance', 'relTol'), 
                     'p_rgh': ('PCG', 'DIC', '1e-8', '0.01'), '"(U|e|k|omega|epsilon)"': ('PBiCGStab', 'DILU;', '1e-7', '0.1')},
              'bsf':{'names': ('type', 'mixture', 'transport', 'thermo', 'equationOfState', 'specie', 'energy'), 
                       'thermoType': ('heRhoThermo', 'pureMixture', 'const', 'hConst', 'perfectGas', 'specie', 'sensibleEnthalpy')}}
    
    specie = {'bbsf': {'names': ('molWeight',), 
                       'specie': ('28.96',)},
              'sf': {'names': ('solver', 'preconditioner', 'tolerance', 'relTol'), 
                     'p_rgh': ('PCG', 'DIC', '1e-8', '0.01'), '"(U|e|k|omega|epsilon)"': ('PBiCGStab', 'DILU;', '1e-7', '0.1')},
              'bsf':{'names': ('molWeight',), 
                       'specie': ('28.96',)}}

    eos = {'bbsf': {'names': ('rho0', 'T0', 'beta'), 
                    'equationOfState': ('1', '300', '3e-03')},
              'sf': {'names': ('solver', 'preconditioner', 'tolerance', 'relTol'), 
                     'p_rgh': ('PCG', 'DIC', '1e-8', '0.01'), '"(U|e|k|omega|epsilon)"': ('PBiCGStab', 'DILU;', '1e-7', '0.1')},
              'bsf':{}}
    thermod = {'bbsf': {'names': ('Cv', 'Hf'), 
                'thermodynamics': ('712', '0')},
              'sf': {'names': ('solver', 'preconditioner', 'tolerance', 'relTol'), 
                     'p_rgh': ('PCG', 'DIC', '1e-8', '0.01'), '"(U|e|k|omega|epsilon)"': ('PBiCGStab', 'DILU;', '1e-7', '0.1')},
              'bsf':{'names': ('Cp', 'Hf'), 
                'thermodynamics': ('1004.4', '0')}}
    trans = {'bbsf': {'names': ('mu', 'Pr'), 
                      'transport': ('1e-05', '0.7')},
              'sf': {'names': ('solver', 'preconditioner', 'tolerance', 'relTol'), 
                     'p_rgh': ('PCG', 'DIC', '1e-8', '0.01'), '"(U|e|k|omega|epsilon)"': ('PBiCGStab', 'DILU;', '1e-7', '0.1')},
              'bsf':{'names': ('mu', 'Pr'), 
                      'transport': ('1e-05', '0.7')}}
    text = ofheader + write_ffile('dictionary', 'constant', 'thermophysicalProperties') + 'thermoType\n{\n'
    
    for sol in thermo[solver]:
        if sol != 'names':
            text += ';\n'.join(['  {} {}'.format(s, thermo[solver][sol][si]) for si, s in enumerate(thermo[solver]['names'])]) + ';\n }\n'
    text += '\n\n mixture\n{\n  '
    
    for sol in specie[solver]:
        if sol != 'names':
            text += sol + '\n  {\n' + ';\n'.join(['    {} {}'.format(s, specie[solver][sol][si]) for si, s in enumerate(specie[solver]['names'])]) + ';\n  }\n  '
    for sol in eos[solver]:
        if sol != 'names':
            text += sol + '\n  {\n' + ';\n'.join(['    {} {}'.format(s, eos[solver][sol][si]) for si, s in enumerate(eos[solver]['names'])]) + ';\n  }\n  '
    for sol in thermod[solver]:
        if sol != 'names':
            text += sol + '\n  {\n' + ';\n'.join(['    {} {}'.format(s, thermod[solver][sol][si]) for si, s in enumerate(thermod[solver]['names'])]) + ';\n  }\n  '        
    for sol in trans[solver]:
        if sol != 'names':
            text += sol + '\n  {\n' + ';\n'.join(['    {} {}'.format(s, trans[solver][sol][si]) for si, s in enumerate(trans[solver]['names'])]) + ';\n  }\n'   
    text +='}'
    return text
    
    # ofheader = 'FoamFile\n{\n    version     2.0;\n    format      ascii;\n    class       dictionary;\n    location    "constant";\n    object      transportProperties;\n}\n\n'
    # if node.transience == '0' and node.turbulence == 'laminar':
    #     return ofheader + 'nu              nu [ 0 2 -1 0 0 0 0 ] 0.01;\n'
    # elif node.transience == '0' and node.turbulence != 'laminar':
    #     return ofheader + 'transportModel  Newtonian;\n\nrho             rho [ 1 -3 0 0 0 0 0 ] 1;\n\nnu              nu [ 0 2 -1 0 0 0 0 ] 1e-05;\n\n' + \
    #     'CrossPowerLawCoeffs\n{\n    nu0             nu0 [ 0 2 -1 0 0 0 0 ] 1e-06;\n    nuInf           nuInf [ 0 2 -1 0 0 0 0 ] 1e-06;\n    m               m [ 0 0 1 0 0 0 0 ] 1;\n' + \
    #     '    n               n [ 0 0 0 0 0 0 0 ] 1;\n}\n\n' + \
    #     'BirdCarreauCoeffs\n{\n    nu0             nu0 [ 0 2 -1 0 0 0 0 ] 1e-06;\n    nuInf           nuInf [ 0 2 -1 0 0 0 0 ] 1e-06;\n' + \
    #     '    k               k [ 0 0 1 0 0 0 0 ] 0;\n    n               n [ 0 0 0 0 0 0 0 ] 1;\n}'
    # elif node.transience == '1' and node.turbulence != 'laminar':   
    #     return ''

def fvmtwrite(node, solver):
    ras = {'bbsf': {'names': ('model', 'turbulence', 'printCoeffs'), 
                    'RAS': ('kEpsilon', 'on', 'on')},
              'sf': {'names': ('solver', 'preconditioner', 'tolerance', 'relTol'), 
                     'p_rgh': ('PCG', 'DIC', '1e-8', '0.01'), '"(U|e|k|omega|epsilon)"': ('PBiCGStab', 'DILU;', '1e-7', '0.1')},
              'bsf':{'names': ('model', 'turbulence', 'printCoeffs'), 
                    'RAS': ('kEpsilon', 'on', 'on')}}
    text = ofheader + write_ffile('dictionary', 'constant', 'momentumTransport') + '\nsimulationType RAS;\n\n'

    for sol in ras[solver]:
        if sol != 'names':
            text += sol + '\n  {\n' + ';\n'.join(['    {} {}'.format(s, ras[solver][sol][si]) for si, s in enumerate(ras[solver]['names'])]) + ';\n}' 
    return text
    
def fvtpwrite(node, solver):    
    ofheader = 'FoamFile\n{\n    version     2.0;\n    format      ascii;\n    class       dictionary;\n    location    "constant";\n    object      transportProperties;\n}\n\n'
    if node.transience == '0' and node.turbulence == 'laminar':
        return ofheader + 'nu              nu [ 0 2 -1 0 0 0 0 ] 0.01;\n'
    elif node.transience == '0' and node.turbulence != 'laminar':
        return ofheader + 'transportModel  Newtonian;\n\nrho             rho [ 1 -3 0 0 0 0 0 ] 1;\n\nnu              nu [ 0 2 -1 0 0 0 0 ] 1e-05;\n\n' + \
        'CrossPowerLawCoeffs\n{\n    nu0             nu0 [ 0 2 -1 0 0 0 0 ] 1e-06;\n    nuInf           nuInf [ 0 2 -1 0 0 0 0 ] 1e-06;\n    m               m [ 0 0 1 0 0 0 0 ] 1;\n' + \
        '    n               n [ 0 0 0 0 0 0 0 ] 1;\n}\n\n' + \
        'BirdCarreauCoeffs\n{\n    nu0             nu0 [ 0 2 -1 0 0 0 0 ] 1e-06;\n    nuInf           nuInf [ 0 2 -1 0 0 0 0 ] 1e-06;\n' + \
        '    k               k [ 0 0 1 0 0 0 0 ] 0;\n    n               n [ 0 0 0 0 0 0 0 ] 1;\n}'
    elif node.transience == '1' and node.turbulence != 'laminar':   
        return ''

def fvrpwrite(node, solver):
    p1dict = {'radiation': 'on', 'radiationModel': 'P1', 'solverFreq': '1', 'absorptionEmissionModel': 'constant',
              'constantCoeffs': {'absorptivity': '0.5', 'emissivity': '0.5', 'E': '0'}, 'scatterModel': 'none', 'sootModel': 'none'}
    fvdomdict = {'radiation': 'on', 'radiationModel': 'fvDOM',
              'fvDOMCoeffs': {'nPhi': '0.5', 'nTheta': '0.5', 'tolerance': '1e-3', 'maxIter': '10'}, 'solverFreq':'10', 'absorptionEmissionModel':'constant', 
              'constantCoeffs': {'absorptivity': '0.5', 'emissivity': '0.5', 'E': '0'}, 'scatterModel': 'none', 'sootModel': 'none'}
    
    raddict = p1dict if node.radmodel == '0' else fvdomdict
    text = ofheader + write_ffile('dictionary', 'constant', 'radiationProperties')
    
    for d in raddict:
        if isinstance(raddict[d], str):
            text += '{} {};\n\n'.format(d,raddict[d]) 
        elif isinstance(raddict[d], dict):
            text += '{}\n{{\n'.format(d)
            for sd in raddict[d]:
                if isinstance(raddict[d][sd], str):
                    text += '    {} {};\n'.format(sd, raddict[d][sd]) 
            text += '}\n'
    return text
    
def fvschwrite(node, solver):  
    print(solver)
    if node.transience == '0' and node.turbulence != 'laminar':
        schdict = {'ddtSchemes': {'default': 'steadyState'}, 
                   'divSchemes': {'default': 'none'},
                   'interpolationSchemes': {'default': 'linear'},
                   'wallDist': {'method': 'meshWave'}}
        if node.turbulence == 'kEpsilon':
            if node.buoyancy:
                if not node.buossinesq:
                    schdict['divSchemes']['div(phi,k)'] = 'bounded Gauss limitedLinear 0.2'
                    schdict['divSchemes']['div(phi,epsilon)'] = 'bounded Gauss limitedLinear 0.2'
                else:
                    schdict['divSchemes']['div(phi,k)'] = 'bounded Gauss upwind'
                    schdict['divSchemes']['div(phi,epsilon)'] = 'bounded Gauss upwind'

        elif node.turbulence == 'kOmega': 
            if node.buoyancy:
                if not node.buossinesq:
                    schdict['divSchemes']['div(phi,k)'] = 'bounded Gauss limitedLinear 0.2'
                    schdict['divSchemes']['div(phi,omega)'] = 'bounded Gauss limitedLinear 0.2'
            else:
                schdict['divSchemes']['div(phi,k)'] = 'bounded Gauss upwind'
                schdict['divSchemes']['div(phi,omega)'] = 'bounded Gauss upwind'

        elif node.turbulence == 'SpalartAllmaras':
            schdict['divSchemes']['div(phi,nuTilda)'] = 'bounded Gauss upwind'

        if node.buoyancy:
            if not node.buossinesq:
                schdict['gradSchemes'] = {'default': 'Gauss linear'}
                schdict['divSchemes']['div(((rho*nuEff)*dev2(T(grad(U)))))'] = 'Gauss linear'
                schdict['divSchemes']['div(phi,Ekp)'] = 'bounded Gauss linear'
                schdict['divSchemes']['div(phi,U)'] = 'bounded Gauss limitedLinear 0.2'
                schdict['divSchemes']['div(phi,K)'] = 'bounded Gauss limitedLinear 0.2'
                schdict['divSchemes']['div(phi,h)'] = 'bounded Gauss limitedLinear 0.2'
                schdict['divSchemes']['laplacianSchemes'] = {'default': 'Gauss linear uncorrected'}
                schdict['divSchemes']['snGradSchemes'] = {'default': 'uncorrected'}
            else:
                schdict['gradSchemes']['default'] = 'Gauss linear'
                schdict['divSchemes']['div(phi,e)'] = 'bounded Gauss upwind'
                schdict['divSchemes']['div(phi,U)'] = 'bounded Gauss upwind'
                schdict['divSchemes']['div(((rho*nuEff)*dev2(T(grad(U)))))'] = 'Gauss linear'
                schdict['laplacianSchemes'] = {'default': 'Gauss linear corrected'}
                schdict['snGradSchemes'] = {'default': 'corrected'}
        else:
            schdict['gradSchemes'] = {'default': 'Gauss linear', 'limited': 'cellLimited Gauss linear 1', 'grad(U)': '$limited'},
            schdict['divSchemes']['div((nuEff*dev2(T(grad(U)))))'] = 'Gauss linear'

    htext = ofheader + write_ffile('dictionary', 'system', 'fvSchemes') 
    ntext = write_fvdict(htext, schdict)  
          
    ddt = {'bbsf': {'names': ('default',), 
                    'ddtSchemes': ('steadyState',)},
              'sf': {'names': ('default',), 
                     'ddtSchemes': ('steadyState',)},
              'bsf':{'names': ('default',), 
                    'ddtSchemes': ('steadyState',)}}
    grad = {'bbsf': {'names': ('default',), 
                    'gradSchemes': ('Gauss linear',)},
              'sf': {'names': ('default', 'limited', 'grad(U)', 'grad(k)', 'grad(epsilon)'), 
                     'gradSchemes': ('Gauss linear', 'cellLimited Gauss linear 1', '$limited', '$limited', '$limited')},
              'bsf':{'names': ('default',), 
                    'gradSchemes': ('Gauss linear',)}}
    div = {'bbsf': {'names': ('default', 'div(phi,U)', 'div(phi,e)', 'div(phi,k)', 'div(phi,epsilon)', 'div(phi,Ekp)', 'div(((rho*nuEff)*dev2(T(grad(U)))))'), 
                    'divSchemes': ('none', 'bounded Gauss upwind', 'bounded Gauss upwind', 'bounded Gauss upwind', 'bounded Gauss upwind', 'bounded Gauss linear', 'Gauss linear')},
              'sf': {'names': ('default', 'div(phi,U)', 'turbulence', 'div(phi,k)', 'div(phi,epsilon)', 'div((nuEff*dev2(T(grad(U)))))'), 
                     'divSchemes': ('none', 'bounded Gauss linearUpwind limited', 'bounded Gauss limitedLinear 1', '$turbulence', '$turbulence', 'Gauss linear')},
              'bsf':{}}
    lap = {'bbsf': {'names': ('default',), 
                    'laplacianSchemes': ('Gauss linear limited corrected 0.33',)},
              'sf': {'names': ('default',), 
                     'laplacianSchemes': ('Gauss linear corrected',)},
              'bsf':{'names': ('default',), 
                     'laplacianSchemes': ('Gauss linear uncorrected',)}}
    terp = {'bbsf': {'names': ('default',), 
                    'interpolationSchemes': ('linear',)},
              'sf': {'names': ('default',), 
                     'interpolationSchemes': ('linear',)},
              'bsf':{'names': ('default',), 
                    'interpolationSchemes': ('linear',)}}
    sng = {'bbsf': {'names': ('default',), 
                    'snGradSchemes': ('limited corrected 0.33',)},
              'sf': {'names': ('default',), 
                     'snGradSchemes': ('corrected',)},
              'bsf':{'names': ('default',), 
                     'snGradSchemes': ('uncorrected',)}}

    text = ofheader + write_ffile('dictionary', 'system', 'fvSchemes')
    for sol in ddt[solver]:
        if sol != 'names':
            text += sol + '\n{\n' + ';\n'.join(['  {} {}'.format(s, ddt[solver][sol][si]) for si, s in enumerate(ddt[solver]['names'])]) + ';\n}' 
    text += '\n'
    for sol in grad[solver]:
        if sol != 'names':
            text += sol + '\n{\n' + ';\n'.join(['  {} {}'.format(s, grad[solver][sol][si]) for si, s in enumerate(grad[solver]['names'])]) + ';\n}'    
    text += '\n'
    for sol in div[solver]:
        if sol != 'names':
            text += sol + '\n{\n' + ';\n'.join(['  {} {}'.format(s, div[solver][sol][si]) for si, s in enumerate(div[solver]['names'])]) + ';\n}'
    text += '\n'
    for sol in lap[solver]:
        if sol != 'names':
            text += sol + '\n{\n' + ';\n'.join(['  {} {}'.format(s, lap[solver][sol][si]) for si, s in enumerate(lap[solver]['names'])]) + ';\n}'
    text += '\n'
    for sol in terp[solver]:
        if sol != 'names':
            text += sol + '\n{\n' + ';\n'.join(['  {} {}'.format(s, terp[solver][sol][si]) for si, s in enumerate(terp[solver]['names'])]) + ';\n}'
    text += '\n'
    for sol in sng[solver]:
        if sol != 'names':
            text += sol + '\n{\n' + ';\n'.join(['  {} {}'.format(s, sng[solver][sol][si]) for si, s in enumerate(sng[solver]['names'])]) + ';\n}'
    text += '\n'
    return ntext       
    # ofheader = 'FoamFile\n{\n  version     2.0;\n  format      ascii;\n  class       dictionary;\n  location    "system";\n  object    fvSchemes;\n}\n\n'
    # if node.solver == 'icoFoam':
    #     return ofheader + 'ddtSchemes\n{\n  default         Euler;\n}\n\ngradSchemes\n{\n  default         Gauss linear;\n  grad(p)         Gauss linear;\n}\n\n' + \
    #         'divSchemes\n{\n  default         none;\n  div(phi,U)      Gauss linear;\n}\n\nlaplacianSchemes\n{\n  default         Gauss linear orthogonal;\n}\n\n' + \
    #         'interpolationSchemes\n{\n  default         linear;\n}\n\n' + \
    #         'snGradSchemes{  default         orthogonal;}\n\nfluxRequired{  default         no;  p;\n}'
    # elif not node.buoyancy:
    #     ofheader += 'ddtSchemes\n{\n    default         steadyState;\n}\n\ngradSchemes\n{\n    default         Gauss linear;\n}\n\ndivSchemes\n{\n    '
    #     if node.turbulence == 'laminar':
    #         ofheader += 'default         none;\n    div(phi,U)   bounded Gauss upwind;\n    div(phi,k)      bounded Gauss upwind;\n    div(phi,epsilon)  bounded Gauss upwind;\n    div((nuEff*dev2(T(grad(U))))) Gauss linear;\n}\n\n'

    #     elif node.turbulence == 'kEpsilon':
    #         ofheader += 'default         none;\n    div(phi,U)   bounded Gauss upwind;\n    div(phi,k)      bounded Gauss upwind;\n    div(phi,epsilon)  bounded Gauss upwind;\n    div((nuEff*dev2(T(grad(U))))) Gauss linear;\n}\n\n'
    #     elif node.turbulence == 'kOmega':
    #         ofheader += 'default         none;\n    div(phi,U)   bounded Gauss upwind;\n    div(phi,k)      bounded Gauss upwind;\n    div(phi,omega)  bounded Gauss upwind;\n    div((nuEff*dev2(T(grad(U))))) Gauss linear;\n}\n\n'
    #     elif node.turbulence == 'SpalartAllmaras':
    #         ofheader += 'default         none;\n    div(phi,U)   bounded Gauss linearUpwind grad(U);\n    div(phi,nuTilda)      bounded Gauss linearUpwind grad(nuTilda);\n    div((nuEff*dev2(T(grad(U))))) Gauss linear;\n}\n\n'
    #     ofheader += 'laplacianSchemes\n{\n    default         Gauss linear corrected;\n}\n\n' + \
    #     'interpolationSchemes\n{\n    default         linear;\n}\n\nsnGradSchemes\n{\n    default         corrected;\n}\n\n' + \
    #     'fluxRequired\n{\n    default         no;\n    p               ;\n}\n\nwallDist\n{\n    method meshWave;\n}\n\n'
    # elif node.buoyancy:
    #     ofheader += 'ddtSchemes\n{\n    default         steadyState;\n}\n\ngradSchemes\n{\n    default         Gauss linear;\n}\n\ndivSchemes\n{\n    '
    #     if node.turbulence == 'laminar':
    #         ofheader += 'default         none;\n    div(phi,U)   bounded Gauss upwind;\n    div(phi,k)      bounded Gauss upwind;\n    div(phi,epsilon)  bounded Gauss upwind;\n    div((nuEff*dev2(T(grad(U))))) Gauss linear;\n}\n\n'

    #     elif node.turbulence == 'kEpsilon':
    #         ofheader += '''default         none;
    # div(phi,U)      bounded Gauss limitedLinear 0.2;
    # div(phi,K)      bounded Gauss limitedLinear 0.2;
    # div(phi,h)      bounded Gauss limitedLinear 0.2;
    # div(phi,k)      bounded Gauss limitedLinear 0.2;
    # div(phi,epsilon) bounded Gauss limitedLinear 0.2;
    # div(phi,omega) bounded Gauss limitedLinear 0.2;
    #     div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;\n}\n\n'''
    
    #     elif node.turbulence == 'kOmega':
    #         ofheader += 'default         none;\n    div(phi,U)   bounded Gauss upwind;\n    div(phi,k)      bounded Gauss upwind;\n    div(phi,omega)  bounded Gauss upwind;\n    div((nuEff*dev2(T(grad(U))))) Gauss linear;\n}\n\n'
    #     elif node.turbulence == 'SpalartAllmaras':
    #         ofheader += 'default         none;\n    div(phi,U)   bounded Gauss linearUpwind grad(U);\n    div(phi,nuTilda)      bounded Gauss linearUpwind grad(nuTilda);\n    div((nuEff*dev2(T(grad(U))))) Gauss linear;\n}\n\n'
    #     ofheader += 'laplacianSchemes\n{\n    default         Gauss linear corrected;\n}\n\n' + \
    #     'interpolationSchemes\n{\n    default         linear;\n}\n\nsnGradSchemes\n{\n    default         corrected;\n}\n\n' + \
    #     'fluxRequired\n{\n    default         no;\n    p               ;\n}\n\nwallDist\n{\n    method meshWave;\n}\n\n'
        
    # return ofheader


    
def fvraswrite(turb):
    ofheader = 'FoamFile\n{\n    version     2.0;\n    format      ascii;\n    class       dictionary;\n    location    "constant";\n    object      turbulenceProperties;\n}\n\n'
    if turb == 'laminar':
        ofheader += 'simulationType laminar;\n\n'
    else:
        ofheader += 'simulationType RAS;\n\n'
        ofheader += 'RAS\n{{\nRASModel        {};\n\nturbulence      on;\n\nprintCoeffs     on;\n}}\n\n'.format(turb)
    return ofheader

def fvgwrite():
    ofheader = '''/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  8
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       uniformDimensionedVectorField;
    location    "constant";
    object      g;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 1 -2 0 0 0 0];
value           (0 0 -9.81);


// ************************************************************************* //'''
    return ofheader



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
    
    layersurf = '({}|{})'.format(kwargs['ground'][0].name, fvos[0].name) if kwargs and kwargs['ground'] else fvos[0].name 
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
    ofheader +='  refinementSurfaces\n  {\n'

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

#
#numberOfSubdomains 16;
#
#method          simple;
#
#simpleCoeffs
#{
#    n               (4 4 1);
#    delta           0.001;
#}
#
#hierarchicalCoeffs
#{
#    n               (1 1 1);
#    delta           0.001;
#    order           xyz;
#}
#
#manualCoeffs
#{
#    dataFile        "";
#}
#
#distributed     no;
#
#roots           ( );
#
#
#// ************************************************************************* //

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
#    bmomw, bmovs = bmo.matrix_world, [vert for vert in bmo.data.vertices]
    for o in fvos:
        with open(os.path.join(scene['flparams']['ofctsfilebase'], '{}.obj'.format(o.name)), 'w') as objfile:
            bm = bmesh.new()
            tempmesh = o.to_mesh(scene = bpy.context.scene, apply_modifiers = True, settings = 'PREVIEW')
            bm.from_mesh(tempmesh)
            bm.transform(o.matrix_world)
            bm.transform(mathutils.Matrix.Translation(bmo['flovi_translate']))
            bpy.data.meshes.remove(tempmesh)
#            omw, ovs = o.matrix_world, [vert for vert in o.data.vertices]
#            xvec, yvec, zvec = (bmomw*bmovs[3].co - bmomw*bmovs[0].co).normalized(), (bmomw*bmovs[2].co - bmomw*bmovs[3].co).normalized(), (bmomw*bmovs[4].co - bmomw*bmovs[0].co).normalized() 
#            ofvpos = [[(omw*ov.co - bmomw*bmovs[0].co)*vec for vec in (xvec, yvec, zvec)] for ov in ovs]
#            bm = bmesh.new()
#            bm.from_mesh(o.data)
#            vcos = ''.join(['v {} {} {}\n'.format(*ofvpo) for ofvpo in ofvpos])    
            vcos =  ''.join(['v {0[0]} {0[1]} {0[2]}\n'.format(v.co) for v in bm.verts])    
            objfile.write(objheader+vcos)
            for m, mat in enumerate(o.data.materials):
                objfile.write('g {}\n'.format(mat.name) + ''.join(['f {} {} {}\n'.format(*[v.index + 1 for v in f.verts]) for f in bmesh.ops.triangulate(bm, faces = bm.faces)['faces'] if f.material_index == m]))
            objfile.write('#{}'.format(len(bm.faces)))
            bm.free()
            
# def fvsolwrite(node):
#     header = ofheader + fileheader('fvSolution')
#     if not node.buoyancy:        
#         text = header +'''solvers
# {
#     p
#     {
#         solver          GAMG;
#         smoother        GaussSeidel;
#         tolerance       1e-6;
#         relTol          0.1;
#     }

#     "(U|k|omega|epsilon)"
#     {
#         solver          smoothSolver;
#         smoother        symGaussSeidel;
#         tolerance       1e-6;
#         relTol          0.1;
#     }
# }

# '''
#     else:
#          text = header +'''solvers
# {    
#     p_rgh
#     {
#         solver          PCG;
#         preconditioner  DIC;
#         tolerance       1e-8;
#         relTol          0.01;
#     }

#     "(U|h|k|epsilon)"
#     {
#         solver          PBiCGStab;
#         preconditioner  DILU;
#         tolerance       1e-7;
#         relTol          0.1;
#     }
# }

# '''
#     if node.turbulence == 'laminar' and node.transience == '0':
#         text += 'PISO\n{\n  nCorrectors     2;\n  nNonOrthogonalCorrectors 0;\n  pRefCell        0;\n  pRefValue       0;\n}\n\n' + \
#         'solvers\n{\n    p\n    {\n        solver          GAMG;\n        tolerance       1e-06;\n        relTol          0.1;\n        smoother        GaussSeidel;\n' + \
#         '        nPreSweeps      0;\n        nPostSweeps     2;\n        cacheAgglomeration true;\n        nCellsInCoarsestLevel 10;\n        agglomerator    faceAreaPair;\n'+ \
#         '        mergeLevels     1;\n    }\n\npFinal\n{\n    $p;\n    relTol 0;\n}\n\n    U\n    {\n        solver          smoothSolver;\n        smoother        GaussSeidel;\n        nSweeps         2;\n' + \
#         '        tolerance       1e-08;\n        relTol          0.1;\n    }\n\n    nuTilda\n    {\n        solver          smoothSolver;\n        smoother        GaussSeidel;\n' + \
#         '        nSweeps         2;\n        tolerance       1e-08;\n        relTol          0.1;\n    }\n}\n\n'
    
#     elif node.turbulence != 'laminar' and node.transience == '0':           
#         text += '''SIMPLE
# {{
#     residualControl
#     {{
#         p               {:.6f};
#         U               {:.6f};
#         "(k|omega|epsilon)" {:.6f};
#     }}
#     nNonOrthogonalCorrectors 0;
#     pRefCell        0;
#     pRefValue       0;

# }}

# '''.format(node.presid, node.uresid, node.keoresid)
        
#         text += '''potentialFlow
# {
#     nNonOrthogonalCorrectors 10;
# }

# '''

#         text += '''relaxationFactors
# {
#     fields
#     {
#         p               0.3;
#     }
#     equations
#     {
#         U               0.7;
#         "(k|omega|epsilon).*" 0.7;
#     }
# }'''
# #        if node.turbulence == 'kEpsilon':
# #            ofheader += 'relaxationFactors\n{\n    fields\n    {\n        p               0.3;\n    }\n    equations\n    {\n' + \
# #            '        U               0.7;\n        k               0.7;\n        epsilon           0.7;\n    }\n}\n\n'
# #        elif node.turbulence == 'kOmega':
# #            ofheader += 'relaxationFactors\n{\n    fields\n    {\n        p               0.3;\n    }\n    equations\n    {\n' + \
# #            '        U               0.7;\n        k               0.7;\n        omega           0.7;\n    }\n}\n\n'
# #        elif node.turbulence == 'SpalartAllmaras':
# #            ofheader += 'relaxationFactors\n{\n    fields\n    {\n        p               0.3;\n    }\n    equations\n    {\n' + \
# #            '        U               0.7;\n        k               0.7;\n        nuTilda           0.7;\n    }\n}\n\n'
#     return text

# def fvtphwrite(buoss):
#     ofheader = '''FoamFile
# {{
#     version     2.0;
#     format      ascii;
#     class       dictionary;
#     location    "constant";
#     object      thermophysicalProperties;
# }}
# // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

# thermoType
# {{
#     type            heRhoThermo;
#     mixture         pureMixture;
#     transport       const;
#     thermo          {}Const;
#     equationOfState {};
#     specie          specie;
#     energy          {};
# }}

# {}

# mixture
# {{
#     specie
#     {{
#         molWeight       28.9;
#     }}
#     thermodynamics
#     {{
#         C{}              {};
#         Hf              0;
#     }}
#     transport
#     {{
#         mu              1.8e-05;
#         Pr              0.7;
#     }}
#     {}
# }}'''.format(('h', 'e')[buoss], ('perfectGas', 'Boussinesq')[buoss], 
#             ('sensibleEnthalpy', 'sensibleInternalEnergy')[buoss],
#             ('pRef            100000;', '')[buoss],
#             ('p', 'v')[buoss],
#             ('1000', '712')[buoss], ('', 'equationOfState\n    {\n        rho0            1;\n        T0              300;\n        beta            3e-03;\n    }\n')[buoss])
#     print(buoss)
#     return ofheader  