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

import bpy, mathutils, colorsys, os
from collections import OrderedDict
from numpy import arange, array
from numpy import sum as nsum
from .vi_func import logentry, facearea
from .vi_dicts import rnu, arnu, zresdict, envdict, enresdict, presdict, lresdict


def ret_areas(o):
    ovp = o.vi_params
    omats = [ms.material for ms in o.material_slots]
    ofa = 0

    for face in o.data.polygons:
        if omats[face.material_index]:
            if omats[face.material_index].vi_params.envi_nodes:
                for node in omats[face.material_index].vi_params.envi_nodes.nodes:
                    if node.bl_idname == 'No_En_Mat_Con' and node.active:
                        if node.envi_con_type == 'Floor':
                            ofa += facearea(o, face)

            omats[face.material_index].vi_params['enparams']['area'] += facearea(o, face)

    return ofa


def get_mat(node, ee):
    for material in bpy.data.materials:
        if node.id_data == material.vi_params.envi_nodes and material.vi_params.envi_export == ee:
            break
    return material


def get_con_node(mvp):
    if mvp.get('envi_nodes'):
        ecnodes = [n for n in mvp.envi_nodes.nodes if n.bl_idname == 'No_En_Mat_Con']
        ecanodes = [n for n in ecnodes if n.active]

        if not ecanodes:
            if not ecnodes[0].active:
                ecnodes[0].active = True
            return ecnodes[0]
        else:
            return ecanodes[0]
    else:
        return None


def get_con_node2(mat):
    mvp = mat.vi_params

    if mvp.get('envi_nodes'):
        ecnodes = [n for n in mvp.envi_nodes.nodes if n.bl_idname == 'No_En_Mat_Con']
        ecanodes = [n for n in ecnodes if n.active]

        if not ecanodes:
            if not ecnodes[0].active:
                ecnodes[0].active = True
            return ecnodes[0]
        else:
            return ecanodes[0]
    else:
        return None


def get_zone_node(coll, enng):
    for node in enng.nodes:
        if node.bl_idname == 'No_En_Net_Zone' and node.zone == coll.name:
            return node


def boundpoly(obj, emnode, poly, enng):
    mat = obj.material_slots[poly.material_index].material

    if emnode.envi_con_con == 'Zone':
        nodes = [node for node in enng.nodes if hasattr(node, 'zone') and node.zone == obj.name]

        for node in nodes:
            insock = node.inputs['{}_{}_b'.format(mat.name, poly.index)]
            outsock = node.outputs['{}_{}_b'.format(mat.name, poly.index)]

            if insock.links:
                bobj = bpy.data.objects[insock.links[0].from_node.zone]
                bpoly = bobj.data.polygons[int(insock.links[0].from_socket.name.split('_')[-2])]
                bmat = bobj.material_slots[bpoly.material_index].material

                if emnode.ret_uv() != get_con_node(bmat.vi_params).ret_uv():
                    logentry('U-values of the paired boundary surfaces {0} and {1} do not match. {1} construction takes precedence'.format(mat.name+'_'+str(poly.index),
                                                                                                                                           insock.links[0].to_node.zone+'_'+str(bpoly.index)))
                    return (('', '', '', ''))
                else:
                    return (("Surface", insock.links[0].from_node.zone+'_'+str(bpoly.index), "NoSun", "NoWind"))

            elif outsock.links:
                bobj = bpy.data.objects[outsock.links[0].to_node.zone]
                bpoly = bobj.data.polygons[int(outsock.links[0].to_socket.name.split('_')[-2])]
                bmat = bobj.data.materials[bpoly.material_index]

                if emnode.ret_uv() != get_con_node(bmat.vi_params).ret_uv():
                    logentry('U-values of the paired boundary surfaces {0} and {1} do not match. {0} construction takes precedence'.format(mat.name+'_'+str(poly.index),
                                                                                                                                           outsock.links[0].to_node.zone+'_'+str(bpoly.index)))
                    return (("Zone", bobj.name, "NoSun", "NoWind"))
                else:
                    return (("Surface", outsock.links[0].to_node.zone+'_'+str(bpoly.index), "NoSun", "NoWind"))

            else:
                return (("Adiabatic", "", "NoSun", "NoWind"))

    elif emnode.envi_con_con == 'Thermal mass':
        return (("Adiabatic", "", "NoSun", "NoWind"))
    elif emnode.envi_con_con == 'Ground':
        return (("Ground", "", "NoSun", "NoWind"))
    else:
        return (("Outdoors", "", "SunExposed", "WindExposed"))


def retenresdict(scene):
    return {'Temp': ('Temperature (degC)', scene.en_temp_max, scene.en_temp_min, u"\u00b0C"), 'Hum': ('Humidity (%)', scene.en_hum_max, scene.en_hum_min, '%'),
            'CO2': ('CO2 (ppm)', scene.en_co2_max, scene.en_co2_min, 'ppm'), 'Heat': ('Heating (W)', scene.en_heat_max, scene.en_heat_min, 'W'),
            'Cool': ('Cooling (W)', scene.en_cool_max, scene.en_cool_min, 'W'),
            'PMV': ('PMV', scene.en_pmv_max, scene.en_pmv_min, 'PMV'), 'PPD': ('PPD (%)', scene.en_ppd_max, scene.en_ppd_min, 'PPD'),
            'SHG': ('Solar gain (W)', scene.en_ppd_max, scene.en_ppd_min, 'SHG'),
            'MaxHeat': ('Max heating (W)', scene.en_maxheat_max, scene.en_maxheat_min, 'W'), 'MaxTemp': ('Max temp (C)', scene.en_maxtemp_max, scene.en_maxtemp_min, u"\u00b0C"),
            'HRheat': ('HR heating (W)', scene.en_hrheat_max, scene.en_hrheat_min, 'hrW')}


def resnameunits():
    return [bpy.props.BoolProperty(name=rnu[str(rnum)][0], description=rnu[str(rnum)][1], default=False) for rnum in range(len(rnu))]


def aresnameunits():
    return [bpy.props.BoolProperty(name=arnu[str(arnum)][0], description=arnu[str(arnum)][1], default=False) for arnum in range(len(arnu))]


def enresprops(disp):
    return {'0': (0, "restt{}".format(disp), "resh{}".format(disp), 0, "restwh{}".format(disp), "restwc{}".format(disp), 0,
                  "ressah{}".format(disp), "reshrhw{}".format(disp), 0, "ressac{}".format(disp), "reswsg{}".format(disp), 0,
                  "resfhb{}".format(disp), "resoeg{}".format(disp)),
            '1': (0, "rescpp{}".format(disp), "rescpm{}".format(disp), 0, 'resmrt{}'.format(disp), 'resocc{}'.format(disp)),
            '2': (0, "resim{}".format(disp), "resiach{}".format(disp), 0, "resco2{}".format(disp), "resihl{}".format(disp)),
            '3': (0, "resl12ms{}".format(disp), "reslof{}".format(disp), 0, "resldp{}".format(disp)),
            '4': (0, "respve{}".format(disp), "respvw{}".format(disp), 0, "respveff{}".format(disp), "respvt{}".format(disp))}


def recalculate_text(scene):
    resdict = {'Temp': ('envi_temp', u'\u00b0C'), 'Hum': ('envi_hum', '%'), 'CO2': ('envi_co2', 'ppm'), 'Heat': ('envi_heat', 'hW'), 'Cool': ('envi_cool', 'cW'),
               'PPD': ('envi_ppd', 'PPD'), 'PMV': ('envi_pmv', 'PMV'), 'SHG': ('envi_shg', 'SHG'), 'HRheat': ('envi_hrheat', 'hrW'),
               'MaxTemp': ('envi_maxtemp', u'Max\u00b0C'), 'MaxHeat': ('envi_maxheat', 'MaxW')}
    resstring = retenvires(scene)

    for res in resdict:
        for o in [o for o in bpy.data.objects if o.get('VIType') and o['VIType'] == resdict[res][0] and o.children]:
            txt = o.children[0]
            sf = scene.frame_current if scene.frame_current <= scene.frame_end else scene.frame_end
            txt.data.body = ("{:.1f}", "{:.0f}")[res in ('MaxHeat', 'Heat', 'Cool', 'SHG', 'CO2', 'HRheat')].format(o[resstring][res][sf]) + resdict[res][1]


def retenvires(scene):
    if scene.en_disp_type == '0':
        if scene['enparams']['fs'] == scene['enparams']['fe']:
            resstring = 'envires'
        else:
            resstring = 'envires{}'.format(bpy.data.node_groups[scene['viparams']['resnode'].split('@')[1]].nodes[scene['viparams']['resnode'].split('@')[0]]['AStart'])
    else:
        resstring = 'envires'
    return resstring


def envilres(scene, resnode):
    for rd in resnode['resdict']:
        if resnode['resdict'][rd][0][:4] == 'WIN-':
            baseob = [o for o in bpy.data.objects if o.name.upper() == resnode['resdict'][rd][0][7:][:-2]][0]
            basefacecent = baseob.matrix_world * baseob.data.polygons[int(resnode['resdict'][rd][0][4:].split('_')[-1])].center

            if scene.envi_flink:
                posobs = [o for o in bpy.data.objects if o.vi_type == '0' and o.layers[0]]
                dists = [(o.location - basefacecent).length for o in posobs]
                resob = posobs[dists.index(min(dists))]
                if not resob.get('envires'):
                    resob['envires'] = {}
            else:
                resob = baseob

            if resob.data.shape_keys and resnode['resdict'][rd][1] == 'Opening Factor':
                resob['envires']['LOF'] = resnode['allresdict'][rd]

                for frame in range(scene.frame_start, scene.frame_end + 1):
                    scene.frame_set(frame)
                    resob.data.shape_keys.key_blocks[1].value = resob['envires']['LOF'][frame]
                    resob.data.shape_keys.key_blocks[1].keyframe_insert(data_path='value', frame=frame)

            if resob.data.shape_keys and resnode['resdict'][rd][1] == 'Linkage Flow in':
                bpy.ops.mesh.primitive_cone_add()
                fcone = bpy.context.active_object
                fcone.rotation_euler = resob.rotation_euler if scene.envi_flink else mathutils.angle(fcone.matrix_world * fcone.data.polygons[-1].normal,
                                                                                                     resob.matrix_word * resob.data.polygons[int(resnode['resdict'][rd][0].split('_')[-1])].normal)
                fcone.parent = resob
                fcone['envires'] = {}
                fi = resnode['allresdict'][rd]

                for frd in resnode['resdict']:
                    if resnode['resdict'][frd][0] == resnode['resdict'][rd][0] and resnode['resdict'][frd][1] == 'Linkage Flow out':
                        fo = resnode['allresdict'][frd]
                fcone['envires']['flow'] = [float(fival) - float(foval) for fival, foval in zip(fi, fo)]

                for frame in range(scene.frame_start, scene.frame_end + 1):
                    scene.frame_set(frame)
                    fcone.rotation_euler = fcone.rotation_euler.to_matrix().inverted().to_euler()
                    fcone.scale = [10*float(fcone['envires']['flow'][frame]) for i in range(3)]
                    fcone.keyframe_insert(data_path='scale', frame=frame)
                    fcone.keyframe_insert(data_path='rotation_euler', frame=frame)


def envizres(scene, eresobs, resnode, restype):
    rl = list(resnode['reslists'])
    zrl = list(zip(*rl))
    resdict = retenresdict(scene)

    if scene.en_disp_type == '0':
        resstart = 24 * (resnode['Start'] - resnode.dsdoy)
        frames = [str(scene['enparams']['fs'])] if scene['enparams']['fs'] == scene['enparams']['fe'] else [str(f) for f in range(scene['enparams']['fs'], scene['enparams']['fe'] + 1)]
        resend = resstart + 24 * (1 + resnode['End'] - resnode['Start'])
        resstrings = ['envires'] if scene['enparams']['fs'] == scene['enparams']['fe'] else ['envires{}'.format(f) for f in range(scene['enparams']['fs'], scene['enparams']['fe'] + 1)]

    elif scene.en_disp_type == '1':
        resstart = resnode['AStart']
        frames = ['All']
        resend = resnode['AEnd'] + 1
        resstrings = ['envires']

    maxval = max([[max(float(r) for r in zrl[4][ri].split())][0] for ri, r in enumerate(zrl[3]) if r == resdict[restype][0] and zrl[1][ri] == 'Zone temporal'])
    minval = min([[min(float(r) for r in zrl[4][ri].split())][0] for ri, r in enumerate(zrl[3]) if r == resdict[restype][0] and zrl[1][ri] == 'Zone temporal'])

    for eo in eresobs:
        o = bpy.data.objects[eo[3:]]
        opos = o.matrix_world * mathutils.Vector([sum(ops)/8 for ops in zip(*o.bound_box)])

        if not any([oc['VIType'] == 'envi_{}'.format(restype.lower()) for oc in o.children if oc.get('VIType')]):
            if scene.en_disp == '1':
                bpy.ops.mesh.primitive_plane_add()
            elif scene.en_disp == '0':
                bpy.ops.mesh.primitive_circle_add(fill_type='NGON')

            ores = bpy.context.active_object
            ores['VIType'] = 'envi_{}'.format(restype.lower())

            for rs, resstring in enumerate(resstrings):
                valstring = [r[4].split()[resstart:resend] for r in rl if r[0] == frames[rs] and r[2] == eo.upper() and r[3] == resdict[restype][0]]
                vals = [float(v) for v in valstring[0]]

                if not ores.get(resstring):
                    ores[resstring] = {}
                    ores[resstring][restype] = vals

            bpy.ops.object.editmode_toggle()
            bpy.ops.mesh.extrude_region_move(MESH_OT_extrude_region={"mirror": False}, TRANSFORM_OT_translate={"value": (0, 0, 1), "constraint_axis": (False, False, True),
                                                                                                               "constraint_orientation": 'NORMAL', "mirror": False,
                                                                                                               "proportional": 'DISABLED', "proportional_edit_falloff": 'SMOOTH',
                                                                                                               "proportional_size": 1, "snap": False, "snap_target": 'CLOSEST',
                                                                                                               "snap_point": (0, 0, 0), "snap_align": False, "snap_normal": (0, 0, 0),
                                                                                                               "gpencil_strokes": False, "texture_space": False, "remove_on_cancel": False,
                                                                                                               "release_confirm": False})
            bpy.ops.object.editmode_toggle()
            ores.scale, ores.parent = (0.25, 0.25, 0.25), o
            ores.location = o.matrix_world.inverted() * opos
            bpy.ops.object.material_slot_add()
            mat = bpy.data.materials.new(name='{}_{}'.format(o.name, restype.lower()))
            ores.material_slots[0].material = mat
            bpy.ops.object.text_add(radius=1, view_align=False, enter_editmode=False, layers=(True, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False))
            txt = bpy.context.active_object
            bpy.context.object.data.extrude = 0.005
            bpy.ops.object.material_slot_add()
            txt.parent = ores
            txt.location, txt.scale = (0, 0, 0), (ores.scale[0]*2, ores.scale[1]*2, 1)
            txt.data.align_x, txt.data.align_y = 'CENTER', 'CENTER'
            txt.name = '{}_{}_text'.format(o.name, restype)
            tmat = bpy.data.materials.new(name='{}'.format(txt.name))
            tmat.diffuse_color = (0, 0, 0)
            txt.material_slots[0].material = tmat
        else:
            ores = [o for o in o.children if o.get(resstrings[0]) and restype in o[resstrings[0]]][0]
            mat = ores.material_slots[0].material

            for rs, resstring in enumerate(resstrings):
                valstring = [r[4].split()[resstart:resend] for r in rl if r[0] == frames[rs] and r[2] == eo.upper() and r[3] == resdict[restype][0]]
                vals = [float(v) for v in valstring[0]] if valstring else []
                ores[resstring][restype] = vals

            if not ores.children:
                bpy.ops.object.text_add(radius=1, view_align=False, enter_editmode=False, layers=(True, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False))
                txt = bpy.context.active_object
                bpy.context.object.data.extrude = 0.005
                bpy.ops.object.material_slot_add()
                txt.parent, txt.location, txt.scale, txt.name = ores, (0, 0, 0), (ores.scale[0], ores.scale[1], 1), '{}_{}_text'.format(o.name, restype)
                txt.data.align_x, txt.data.align_y = 'CENTER', 'CENTER'
                tmat = bpy.data.materials.new(name='{}'.format(txt.name))
                tmat.diffuse_color = (0, 0, 0)
                txt.material_slots[0].material = tmat
            else:
                txt = ores.children[0]

        txt.data.body = "{:.1f}{}".format(ores[resstring][restype][0], resdict[restype][3]) if restype not in ('SHG', 'CO2') else "{:.0f}{}".format(ores[resstring][restype][0], resdict[restype][2])

        if maxval - minval:
            scalevel = [(vals[frame] - minval)/(maxval - minval) for frame in range(0, len(vals))] if maxval - minval else [0] * len(vals)
            colval = [colorsys.hsv_to_rgb(0.667 * (maxval - vals[vi])/(maxval - minval), 1, 1) for vi in range(len(vals))]
        else:
            scalevel = colval = [0] * len(vals)

        sv = [(sv, 0.1)[sv <= 0.1] for sv in scalevel]
        cv = [(((0, 1)[vals[c] >= maxval], 0, (0, 1)[vals[c] <= minval]), cv)[minval < vals[c] < maxval] for c, cv in enumerate(colval)]
        ores.animation_data_clear()
        ores.animation_data_create()
        ores['max'], ores['min'], ores['cmap'], ores['days'], ores['hours'] = maxval, minval, scene.vi_leg_col, o['days'], o['hours']
        ores.animation_data.action = bpy.data.actions.new(name="EnVi Zone")
        oresz = ores.animation_data.action.fcurves.new(data_path="scale", index=2)
        oresz.keyframe_points.add(len(sv))
        mat.animation_data_clear()
        mat.animation_data_create()
        mat.animation_data.action = bpy.data.actions.new(name="EnVi Zone Material")
        mdcr = mat.animation_data.action.fcurves.new(data_path="diffuse_color", index=0)
        mdcg = mat.animation_data.action.fcurves.new(data_path="diffuse_color", index=1)
        mdcb = mat.animation_data.action.fcurves.new(data_path="diffuse_color", index=2)
        mdcr.keyframe_points.add(len(sv))
        mdcg.keyframe_points.add(len(sv))
        mdcb.keyframe_points.add(len(sv))
        txt.animation_data_clear()
        txt.animation_data_create()
        txt.animation_data.action = bpy.data.actions.new(name="EnVi Zone Text")
        txtl = txt.animation_data.action.fcurves.new(data_path="location", index=2)
        txtl.keyframe_points.add(len(sv))

        for frame in range(len(sv)):
            oresz.keyframe_points[frame].co = frame, sv[frame]
            mdcr.keyframe_points[frame].co = frame, cv[frame][0]
            mdcg.keyframe_points[frame].co = frame, cv[frame][1]
            mdcb.keyframe_points[frame].co = frame, cv[frame][2]
            txtl.keyframe_points[frame].co = frame, 1


def epentry(header, params, paramvs):
    return '{}\n'.format(header+(',', '')[header == ''])+'\n'.join([('    ', '')[header == '']+'{:{width}}! - {}'.format(str(pv[0])+(',', ';')[pv[1] == params[-1]],
                                                                    pv[1], width=80 + (0, 4)[header == '']) for pv in zip(paramvs, params)]) + ('\n\n', '')[header == '']


def epschedwrite(name, stype, ts, fs, us):
    params = ['Name', 'Schedule Type Limits Name']
    paramvs = [name, stype]

    for t in range(len(ts)):
        params.append('Field {}'.format(len(params)-2))
        paramvs .append(ts[t])

        for f in range(len(fs[t])):
            params.append('Field {}'.format(len(params)-2))
            paramvs.append(fs[t][f])

            for u in range(len(us[t][f])):
                params.append('Field {}'.format(len(params)-2))
                paramvs.append(us[t][f][u][0])

    return epentry('Schedule:Compact', params, paramvs)


def enunits(self, context):
    try:
        resstring = retenvires(context.scene)
        return [(k, k, 'Display {}'.format(k)) for k in sorted(context.active_object[resstring].keys())]
    except:
        return [('', '', '')]


def enpunits(self, context):
    try:
        resstring = retenvires(context.scene)
        return [(k, k, 'Display {}'.format(k)) for k in context.active_object[resstring].keys()]
    except:
        return []


def enparametric(self, context):
    try:
        resnode = bpy.data.node_groups[context.scene['viparams']['resnode'].split('@')[1]].nodes[context.scene['viparams']['resnode'].split('@')[0]]
        rl = resnode['reslists']
        zrl = list(zip(*rl))

        if len(set(zrl[0])) > 1:
            return [("0", "Static", "Static results"), ("1", "Parametric", "Parametric results")]
        else:
            return [("0", "Static", "Static results")]
    except Exception:
        return [("0", "Static", "Static results")]


def retrmenus(innode, node, axis, zrl):
    ftype = [(frame, frame, "Plot "+frame) for frame in list(OrderedDict.fromkeys(zrl[0])) if frame != 'All']
    frame = 'All' if node.parametricmenu == '1' and len(ftype) > 1 else zrl[0][0]

    invalids = ['Zone temporal'] if frame == 'All' else ['Embodied carbon', 'Zone spatial']

    if axis == 'X-axis':
        rtypes = list(OrderedDict.fromkeys([zrl[1][ri] for ri, r in enumerate(zrl[1]) if zrl[0][ri] == frame and zrl[1][ri] not in invalids]))
    else:
        rtypes = list(OrderedDict.fromkeys([zrl[1][ri] for ri, r in enumerate(zrl[1]) if zrl[0][ri] == frame and zrl[1][ri] not in ['Time'] + invalids]))

    rtype = [(metric, metric, "Plot " + metric) for metric in rtypes]
    rtype = rtype if rtype else [('None', 'None', 'None')]
    ctype = [(metric, metric, "Plot " + metric) for m, metric in enumerate(zrl[3]) if zrl[1][m] == 'Climate' and zrl[0][m] == frame]
    ctype = ctype if ctype else [('None', 'None', 'None')]
    ztypes = list(OrderedDict.fromkeys([metric for m, metric in enumerate(zrl[2]) if zrl[1][m] == ('Zone temporal', 'Zone spatial')[frame == 'All'] and zrl[0][m] == frame]))

    # if frame == 'All':
    #     ztypes += list(OrderedDict.fromkeys([metric for m, metric in enumerate(zrl[2]) if zrl[1][m] == 'Zone spatial' and zrl[0][m] == frame]))
    #     ztypes += list(OrderedDict.fromkeys([metric for m, metric in enumerate(zrl[2]) if zrl[1][m] == 'Embodied carbon' and zrl[0][m] == frame]))

    # if frame == 'All':
    #     ztypes += list(OrderedDict.fromkeys([metric for m, metric in enumerate(zrl[2]) if zrl[1][m] == 'Embodied carbon' and zrl[0][m] == frame]))
    #     print(ztypes)
    ztype = [(metric, metric, "Plot " + metric) for metric in ztypes]
    ztype = ztype if ztype else [('None', 'None', 'None')]

    # try:
    #     zrtype = [(zr[3], zr[3], 'Plot {}'.format(zr[3])) for zr in rl if zr[2] == node.zonemenu and zr[0] == node.framemenu] if node.parametricmenu == '0' else [(zr[3], zr[3], 'Plot {}'.format(zr[3])) for zr in rl if zr[2] == node.zonemenu and zr[0] == 'All']
    # except Exception as e:
    #     zrtype = [('None', 'None', 'None')]

    ptypes = list(OrderedDict.fromkeys([metric for m, metric in enumerate(zrl[2]) if zrl[1][m] == 'Position' and zrl[0][m] == frame]))
    ptype = [(metric, metric, "Plot " + metric) for metric in ptypes]
    ptype = ptype if ptype else [('None', 'None', 'None')]
    prtypes = list(OrderedDict.fromkeys([metric for m, metric in enumerate(zrl[3]) if zrl[1][m] == 'Position' and zrl[0][m] == frame]))
    prtype = [(metric, metric, "Plot " + metric) for metric in prtypes]
    prtype = prtype if prtype else [('None', 'None', 'None')]
    camtypes = list(OrderedDict.fromkeys([metric for m, metric in enumerate(zrl[2]) if zrl[1][m] == 'Camera' and zrl[0][m] == frame]))
    camtype = [(metric, metric, "Plot " + metric) for metric in camtypes]
    camtype = camtype if camtype else [('None', 'None', 'None')]
    camrtypes = list(OrderedDict.fromkeys([metric for m, metric in enumerate(zrl[3]) if zrl[1][m] == 'Camera' and zrl[0][m] == frame]))
    camrtype = [(metric, metric, "Plot " + metric) for metric in camrtypes]
    ltypes = list(OrderedDict.fromkeys([metric for m, metric in enumerate(zrl[2]) if zrl[1][m] == 'Linkage' and zrl[0][m] == frame]))
    ltype = [(metric, metric, "Plot " + metric) for metric in ltypes]
    ltype = ltype if ltype else [('None', 'None', 'None')]
    lrtypes = list(OrderedDict.fromkeys([metric for m, metric in enumerate(zrl[3]) if zrl[1][m] == 'Linkage' and zrl[0][m] == frame]))
    lrtype = [(metric, metric, "Plot " + metric) for metric in lrtypes]
    lrtype = lrtype if lrtype else [('None', 'None', 'None')]
    entypes = list(OrderedDict.fromkeys([metric for m, metric in enumerate(zrl[2]) if zrl[1][m] == 'External' and zrl[0][m] == frame]))
    entype = [(metric, metric, "Plot " + metric) for metric in entypes]
    entype = entype if entype else [('None', 'None', 'None')]
    enrtypes = list(OrderedDict.fromkeys([metric for m, metric in enumerate(zrl[3]) if zrl[1][m] == 'External' and zrl[0][m] == frame]))
    enrtype = [(metric, metric, "Plot " + metric) for metric in enrtypes]
    enrtype = enrtype if enrtype else [('None', 'None', 'None')]
    powtypes = list(OrderedDict.fromkeys([metric for m, metric in enumerate(zrl[2]) if zrl[1][m] == 'Power' and zrl[0][m] == frame]))
    powtype = [(metric, metric, "Plot " + metric) for metric in powtypes]
    powtype = powtype if powtype else [('None', 'None', 'None')]
    powrtypes = list(OrderedDict.fromkeys([metric for m, metric in enumerate(zrl[3]) if zrl[1][m] == 'Power' and zrl[0][m] == frame]))
    powrtype = [(metric, metric, "Plot " + metric) for metric in powrtypes]
    powrtype = powrtype if powrtype else [('None', 'None', 'None')]
    probetypes = list(OrderedDict.fromkeys([metric for m, metric in enumerate(zrl[2]) if zrl[1][m] == 'Probe' and zrl[0][m] == frame]))
    probetype = [(metric, metric, "Plot " + metric) for metric in probetypes]
    probetype = probetype if probetype else [('None', 'None', 'None')]
    probertypes = list(OrderedDict.fromkeys([metric for m, metric in enumerate(zrl[3]) if zrl[1][m] == 'Probe' and zrl[0][m] == frame]))
    probertype = [(metric, metric, "Plot " + metric) for metric in probertypes]
    probertype = probertype if probertype else [('None', 'None', 'None')]
    ectypes = list(OrderedDict.fromkeys([metric for m, metric in enumerate(zrl[2]) if zrl[1][m] == 'Embodied carbon' and zrl[0][m] == frame]))
    ectype = [(metric, metric, "Plot " + metric) for metric in ectypes]
    ectype = ectype if ectype else [('None', 'None', 'None')]
    ecrtypes = list(OrderedDict.fromkeys([metric for m, metric in enumerate(zrl[3]) if zrl[1][m] == 'Embodied carbon' and zrl[0][m] == frame]))
    ecrtype = [(metric, metric, "Plot " + metric) for metric in ecrtypes]
    ecrtype = ecrtype if ecrtype else [('None', 'None', 'None')]
    fmenu = bpy.props.EnumProperty(items=ftype, name="", description="Frame number", default=ftype[0][0])
    rtypemenu = bpy.props.EnumProperty(items=rtype, name="", description="Result types", default=rtype[0][0])
    statmenu = bpy.props.EnumProperty(items=[('Average', 'Average', 'Average Value'), ('Maximum', 'Maximum', 'Maximum Value'), ('Minimum', 'Minimum', 'Minimum Value'), ('Sum', 'Sum', 'Sum Value')],
                                      name="", description="Zone result", default='Average')
#    valid = ['Vi Results']
    climmenu = bpy.props.EnumProperty(items=ctype, name="", description="Climate type", default=ctype[0][0]) if ctype else ''
    zonemenu = bpy.props.EnumProperty(items=ztype, name="", description="Zone", default=ztype[0][0]) if ztype else ''
    zonermenu = bpy.props.EnumProperty(items=zrupdate, name="", description="Zone result")  # if ztype else ''
    linkmenu = bpy.props.EnumProperty(items=ltype, name="", description="Flow linkage", default=ltype[0][0]) if ltype else ''
    linkrmenu = bpy.props.EnumProperty(items=lrtype, name="", description="Flow linkage result", default=lrtype[0][0]) if ltype else ''
    enmenu = bpy.props.EnumProperty(items=entype, name="", description="External node", default=entype[0][0]) if entype else ''
    enrmenu = bpy.props.EnumProperty(items=enrtype, name="", description="External node result", default=enrtype[0][0]) if entype else ''
    posmenu = bpy.props.EnumProperty(items=ptype, name="", description="Position", default=ptype[0][0]) if ptype else ''
    posrmenu = bpy.props.EnumProperty(items=prtype, name="", description="Position result", default=prtype[0][0]) if ptypes else ''
    cammenu = bpy.props.EnumProperty(items=camtype, name="", description="Camera", default=camtype[0][0]) if camtype else ''
    camrmenu = bpy.props.EnumProperty(items=camrtype, name="", description="Camera result", default=camrtype[0][0]) if camtypes else ''
    powmenu = bpy.props.EnumProperty(items=powtype, name="", description="Power", default=powtype[0][0]) if powtype else ''
    powrmenu = bpy.props.EnumProperty(items=powrtype, name="", description="Power result", default=powrtype[0][0]) if powrtype else ''
    probemenu = bpy.props.EnumProperty(items=probetype, name="", description="Probe", default=probetype[0][0]) if probetype else ''
    probermenu = bpy.props.EnumProperty(items=probertype, name="", description="Probe result", default=probertype[0][0]) if probertype else ''
    ecmenu = bpy.props.EnumProperty(items=ectype, name="", description="EC", default=ectype[0][0]) if ectype else ''
    ecrmenu = bpy.props.EnumProperty(items=ecrupdate, name="", description="EC result") # if ecrtype else ''
    multmenu = bpy.props.FloatProperty(name="", description="Result multiplication factor", min=-10000, max=10000, default=1)
    return (fmenu, rtypemenu, climmenu, zonemenu, zonermenu, linkmenu, linkrmenu, enmenu, enrmenu, posmenu, posrmenu,
            cammenu, camrmenu, powmenu, powrmenu, probemenu, probermenu, ecmenu, ecrmenu, multmenu, statmenu)


def processh(lines, znlist):
    hdict, li = {}, 0

    for li, line in enumerate(lines):
        linesplit = line.strip('\n').split(',')

        if len(linesplit) > 3:
            if linesplit[2] == 'Day of Simulation[]':
                hdict[linesplit[0]] = ['Time']
            elif linesplit[3] in envdict:
                hdict[linesplit[0]] = ['Climate',  '', envdict[linesplit[3]]]
            elif linesplit[3] in zresdict and retzonename(linesplit[2]) in znlist:
                hdict[linesplit[0]] = ['Zone temporal',  retzonename(linesplit[2]),  zresdict[linesplit[3]]]
            elif linesplit[3] in enresdict and 'ExtNode' in linesplit[2]:
                hdict[linesplit[0]] = ['External',  linesplit[2],  enresdict[linesplit[3]]]
            elif linesplit[3] in lresdict:
                hdict[linesplit[0]] = ['Linkage',  linesplit[2],  lresdict[linesplit[3]]]
            elif linesplit[3] in presdict:
                hdict[linesplit[0]] = ['Power',  linesplit[2],  presdict[linesplit[3]]]

        if line == 'End of Data Dictionary\n':
            break

    return hdict,  li + 1


def retzonename(zn):
    if zn[-10:] == '_OCCUPANCY':
        return zn[:-10]
    elif zn[-4:] == '_AIR':
        return zn[:-4]
    else:
        return zn


def checkenvierrors(file, sim_op):
    efile = file.read()
    if '** Severe  **' in efile:
        sim_op.report({'ERROR'}, "There is a fatal error in the EnVi model, check the error file in Blender's text editor")


def processf(pro_op, node, con_node):
    scene = bpy.context.scene
    svp = scene.vi_params
    reslists, areslists = con_node['reslists'], []
    frames = range(svp['enparams']['fs'], svp['enparams']['fe'] + 1) if node.bl_label == 'EnVi Simulation' else [scene.frame_current]

    for frame in frames:
        pvps = []
        pow_dict = {}
        en_dict = {}
        ef_dict = {}
        node['envires{}'.format(frame)] = {}
        resfileloc = os.path.join(svp['viparams']['newdir'], '{}{}out.eso'.format(pro_op.resname, frame)) if node.bl_label == 'EnVi Simulation' else node.resfilename

        with open(resfileloc, 'r') as resfile:
            lines = resfile.readlines()
            hdict, lstart = processh(lines, [coll.name.upper() for coll in bpy.data.collections['EnVi Geometry'].children])
            splitlines = [li.strip('\n').split(',') for li in lines[lstart:-2]]
            bdict = {li: ' '.join(['{:.5f}'.format(float(sl[1])) for sl in splitlines if sl[0] == li]) for li in hdict}

            for k in sorted(hdict.keys(), key=int):
                if hdict[k] == ['Time']:
                    reslists.append([str(frame), 'Time', '', 'Month', ' '.join([sl[2] for sl in splitlines if sl[0] == k])])
                    reslists.append([str(frame), 'Time', '', 'Day', ' '.join([sl[3] for sl in splitlines if sl[0] == k])])
                    reslists.append([str(frame), 'Time', '', 'Hour', ' '.join([sl[5] for sl in splitlines if sl[0] == k])])
                    reslists.append([str(frame), 'Time', '', 'DOS', ' '.join([sl[1] for sl in splitlines if sl[0] == k])])
                else:
                    reslists.append([str(frame)] + hdict[k] + [bdict[k]])

            for zn in [coll.name.upper() for coll in bpy.data.collections['EnVi Geometry'].children]:
                for k in sorted(hdict.keys(), key=int):
                    if hdict[k][0] == 'Zone temporal' and hdict[k][1] == zn and hdict[k][2] == 'Ventilation heat (W)':
                        vhs = [float(sl[1]) for sl in splitlines if sl[0] == k]
                    if hdict[k][0] == 'Zone temporal' and hdict[k][1] == zn and hdict[k][2] == 'Heating (W)':
                        hs = [float(sl[1]) for sl in splitlines if sl[0] == k]

                try:
                    reslists.append([str(frame), 'Zone temporal', zn, 'Ventilation heating contribution (W)', ' '.join([str(vhs[i]) if hs[i] > 0.01 else '0' for i in range(len(hs))])])
                    reslists.append([str(frame), 'Zone temporal', zn, 'Fabric heating contribution (W)', ' '.join([str(hs[i] - vhs[i]) if hs[i] > 0.01 else '0' for i in range(len(hs))])])
                except Exception:
                    pass

            for r in reslists:
                if r[0] == str(frame) and r[1] == 'Power' and r[3] == 'PV power (W)':
                    pvps.append([float(r) for r in r[4].split()])

            if pvps:
                reslists.append([str(frame), 'Power', 'All', 'PV power (W)', ' '.join(['{:.3f}'.format(sum(pvpz)) for pvpz in zip(*pvps)])])


    rls = reslists
    zrls = list(zip(*rls))

    try:
        svp['enparams']['lmetrics'] = list(set([zr for zri, zr in enumerate(zrls[3]) if zrls[1][zri] == 'Linkage' and zrls[0][zri] == str(node["AStart"])]))
    except Exception:
        svp['enparams']['lmetrics'] = []

    try:
        svp['enparams']['zmetrics'] = list(set([zr for zri, zr in enumerate(zrls[3]) if zrls[1][zri] == 'Zone temporal' and zrls[0][zri] == str(node["AStart"])]))
    except Exception:
        svp['enparams']['zmetrics'] = []

    zonerls = [zonerl for zonerl in rls if zonerl[1] == 'Zone temporal']
    zzonerls = list(zip(*zonerls))

    # if ecs:
    #     pass

    if len(frames) > 1:
        areslists = []
        areslists.append(['All', 'Frames', '', 'Frames', ' '.join([str(f) for f in frames])])

        if zzonerls:
            temps = [(zrls[2][zi], [float(t) for t in zrls[4][zi].split()]) for zi, z in enumerate(zrls[1]) if z == 'Zone temporal' and zrls[3][zi] == 'Temperature (degC)']
            hums = [(zrls[2][zi], [float(t) for t in zrls[4][zi].split()]) for zi, z in enumerate(zrls[1]) if z == 'Zone temporal' and zrls[3][zi] == 'Humidity (%)']
            heats = [(zrls[2][zi], [float(t) for t in zrls[4][zi].split()]) for zi, z in enumerate(zrls[1]) if z == 'Zone temporal' and zrls[3][zi] == 'Heating (W)']
            cools = [(zrls[2][zi], [float(t) for t in zrls[4][zi].split()]) for zi, z in enumerate(zrls[1]) if z == 'Zone temporal' and zrls[3][zi] == 'Cooling (W)']
            aheats = [(zrls[2][zi], [float(t) for t in zrls[4][zi].split()]) for zi, z in enumerate(zrls[1]) if z == 'Zone temporal' and zrls[3][zi] == 'Air Heating (W)']
            acools = [(zrls[2][zi], [float(t) for t in zrls[4][zi].split()]) for zi, z in enumerate(zrls[1]) if z == 'Zone temporal' and zrls[3][zi] == 'Air Cooling (W)']
            co2s = [(zrls[2][zi], [float(t) for t in zrls[4][zi].split()]) for zi, z in enumerate(zrls[1]) if z == 'Zone temporal' and zrls[3][zi] == 'CO2 (ppm)']
            comfppds = [(zrls[2][zi], [float(t) for t in zrls[4][zi].split()]) for zi, z in enumerate(zrls[1]) if z == 'Zone temporal' and zrls[3][zi] == 'PPD']
            comfpmvs = [(zrls[2][zi], [float(t) for t in zrls[4][zi].split()]) for zi, z in enumerate(zrls[1]) if z == 'Zone temporal' and zrls[3][zi] == 'PMV']
            shgs = [(zrls[2][zi], [float(t) for t in zrls[4][zi].split()]) for zi, z in enumerate(zrls[1]) if z == 'Zone temporal' and zrls[3][zi] == 'Solar gain (W)']

            for zn in set(zzonerls[2]):
                fas = [bpy.data.collections[zn].vi_params['enparams']['floorarea'][str(f)] for f in frames]
                allfas = all([fa > 0 for fa in fas])

                try:
                    if temps:
                        areslists.append(['All', 'Zone spatial', zn, 'Max temp (C)', ' '.join([str(max(t[1])) for t in temps if t[0] == zn])])
                        areslists.append(['All', 'Zone spatial', zn, 'Min temp (C)', ' '.join([str(min(t[1])) for t in temps if t[0] == zn])])
                        areslists.append(['All', 'Zone spatial', zn, 'Avg temp (C)', ' '.join([str(sum(t[1])/len(t[1])) for t in temps if t[0] == zn])])
                    if hums:
                        areslists.append(['All', 'Zone spatial', zn, 'Max humidity (C)', ' '.join([str(max(h[1])) for h in hums if h[0] == zn])])
                        areslists.append(['All', 'Zone spatial', zn, 'Min humidity (C)', ' '.join([str(min(h[1])) for h in hums if h[0] == zn])])
                        areslists.append(['All', 'Zone spatial', zn, 'Avg humidity (C)', ' '.join([str(sum(h[1])/len(h[1])) for h in hums if h[0] == zn])])
                    if heats:
                        areslists.append(['All', 'Zone spatial', zn, 'Max heating (W)', ' '.join([str(max(h[1])) for h in heats if h[0] == zn])])
                        areslists.append(['All', 'Zone spatial', zn, 'Min heating (W)', ' '.join([str(min(h[1])) for h in heats if h[0] == zn])])
                        areslists.append(['All', 'Zone spatial', zn, 'Avg heating (W)', ' '.join([str(sum(h[1])/len(h[1])) for h in heats if h[0] == zn])])
                        areslists.append(['All', 'Zone spatial', zn, 'Total heating (kWh)', ' '.join([str(sum(h[1])*0.001) for h in heats if h[0] == zn])])

                        if allfas:
                            areslists.append(['All', 'Zone spatial', zn, 'Total heating (kWh/m2)', ' '.join([str(sum(h[1])*0.001/fas[hi]) for hi, h in enumerate([h for h in heats if h[0] == zn])])])

                    if cools:
                        areslists.append(['All', 'Zone spatial', zn, 'Max cooling (W)', ' '.join([str(max(h[1])) for h in cools if h[0] == zn])])
                        areslists.append(['All', 'Zone spatial', zn, 'Min cooling (W)', ' '.join([str(min(h[1])) for h in cools if h[0] == zn])])
                        areslists.append(['All', 'Zone spatial', zn, 'Avg cooling (W)', ' '.join([str(sum(h[1])/len(h[1])) for h in cools if h[0] == zn])])
                        areslists.append(['All', 'Zone spatial', zn, 'Total cooling (kWh)', ' '.join([str(sum(h[1])*0.001) for h in cools if h[0] == zn])])

                        if allfas:
                            areslists.append(['All', 'Zone spatial', zn, 'Total cooling (kWh/m2)', ' '.join([str(sum(c[1])*0.001/fas[ci]) for ci, c in enumerate([c for c in cools if c[0] == zn])])])

                    if aheats:
                        areslists.append(['All', 'Zone spatial', zn, 'Max air heating (W)', ' '.join([str(max(h[1])) for h in aheats if h[0] == zn])])
                        areslists.append(['All', 'Zone spatial', zn, 'Min air heating (W)', ' '.join([str(min(h[1])) for h in aheats if h[0] == zn])])
                        areslists.append(['All', 'Zone spatial', zn, 'Avg air heating (W)', ' '.join([str(sum(h[1])/len(h[1])) for h in aheats if h[0] == zn])])
                        areslists.append(['All', 'Zone spatial', zn, 'Total air heating (kWh)', ' '.join([str(sum(h[1])*0.001) for h in aheats if h[0] == zn])])

                        if allfas:
                            areslists.append(['All', 'Zone spatial', zn, 'Total air heating (kWh/m2)', ' '.join([str(sum(h[1])*0.001/fas[hi]) for hi, h in enumerate([h for h in aheats if h[0] == zn])])])

                    if acools:
                        areslists.append(['All', 'Zone spatial', zn, 'Max air cool (W)', ' '.join([str(max(h[1])) for h in acools if h[0] == zn])])
                        areslists.append(['All', 'Zone spatial', zn, 'Min air cool (W)', ' '.join([str(min(h[1])) for h in acools if h[0] == zn])])
                        areslists.append(['All', 'Zone spatial', zn, 'Avg air cool (W)', ' '.join([str(sum(h[1])/len(h[1])) for h in acools if h[0] == zn])])
                        areslists.append(['All', 'Zone spatial', zn, 'Air cooling (kWh)', ' '.join([str(sum(h[1])*0.001) for h in acools if h[0] == zn])])

                        if allfas:
                            areslists.append(['All', 'Zone spatial', zn, 'Total air cooling (kWh/m2)', ' '.join([str(sum(c[1])*0.001/fas[ci]) for ci, c in enumerate([c for c in acools if c[0] == zn])])])

                    if co2s:
                        areslists.append(['All', 'Zone spatial', zn, 'Max CO2 (ppm)', ' '.join([str(max(t[1])) for t in co2s if t[0] == zn])])
                        areslists.append(['All', 'Zone spatial', zn, 'Min CO2 (ppm)', ' '.join([str(min(t[1])) for t in co2s if t[0] == zn])])
                        areslists.append(['All', 'Zone spatial', zn, 'Avg CO2 (ppm)', ' '.join([str(sum(t[1])/len(t[1])) for t in co2s if t[0] == zn])])

                    if comfppds:
                        areslists.append(['All', 'Zone spatial', zn, 'Max PPD', ' '.join([str(max(t[1])) for t in comfppds if t[0] == zn])])
                        areslists.append(['All', 'Zone spatial', zn, 'Min PPD', ' '.join([str(min(t[1])) for t in comfppds if t[0] == zn])])
                        areslists.append(['All', 'Zone spatial', zn, 'Avg PPD', ' '.join([str(sum(t[1])/len(t[1])) for t in comfppds if t[0] == zn])])

                    if comfpmvs:
                        areslists.append(['All', 'Zone spatial', zn, 'Max PMV', ' '.join([str(max(t[1])) for t in comfpmvs if t[0] == zn])])
                        areslists.append(['All', 'Zone spatial', zn, 'Min PMV', ' '.join([str(min(t[1])) for t in comfpmvs if t[0] == zn])])
                        areslists.append(['All', 'Zone spatial', zn, 'Avg PMV', ' '.join([str(sum(t[1])/len(t[1])) for t in comfpmvs if t[0] == zn])])

                    if shgs:
                        areslists.append(['All', 'Zone spatial', zn, 'Max SHG (W)', ' '.join([str(max(t[1])) for t in shgs if t[0] == zn])])
                        areslists.append(['All', 'Zone spatial', zn, 'Min SHG (W)', ' '.join([str(min(t[1])) for t in shgs if t[0] == zn])])
                        areslists.append(['All', 'Zone spatial', zn, 'Avg SHG (W)', ' '.join([str(sum(t[1])/len(t[1])) for t in shgs if t[0] == zn])])
                        areslists.append(['All', 'Zone spatial', zn, 'Total SHG (kWh)', ' '.join([str(sum(t[1])*0.001) for t in shgs if t[0] == zn])])

                        if allfas:
                            areslists.append(['All', 'Zone spatial', zn, 'Total SHG (kWh/m2)', ' '.join([str(sum(s[1])*0.001/fas[si]) for si, s in enumerate([s for s in shgs if s[0] == zn])])])

                except Exception as e:
                    print(e)
                    pro_op.report({'ERROR'}, "There are no zone results to plot. Make sure you have selected valid metrics to calculate and try re-exporting/simulating")
                    return

            if heats and cools:
                try:
                    conds = [sum(x) for x in zip(*[[sum(h[1])*0.001 for h in heats if h[0] == zn], [sum(h[1])*0.001 for h in cools if h[0] == zn]])]
                    areslists.append(['All', 'Zone spatial', zn, 'Total conditioning (kWh)', ' '.join([str(cond) for cond in conds])])

                    if allfas:
                        areslists.append(['All', 'Zone spatial', zn, 'Total conditioning (kWh/m2)', ' '.join([str(cond/fas[ci]) for ci, cond in enumerate([c for c in conds if c[0] == zn])])])
                except:
                    pass
            if aheats and acools:
                try:
                    aconds = [sum(x) for x in zip(*[[sum(h[1])*0.001 for h in aheats if h[0] == zn], [sum(h[1])*0.001 for h in acools if h[0] == zn]])]
                    areslists.append(['All', 'Zone spatial', zn, 'Total air conditioning (kWh)', ' '.join([str(cond) for cond in conds])])

                    if allfas:
                        areslists.append(['All', 'Zone spatial', zn, 'Total air conditioing (kWh/m2)', ' '.join([str(acond/fas[ai]) for ai, acond in enumerate([a for a in aconds if a[0] == zn])])])
                except:
                    pass

        powrls = [powrl for powrl in rls if powrl[1] == 'Power']
        zpowrls = list(zip(*powrls))
        pows = [(zrls[2][zi], [float(t) for t in zrls[4][zi].split()]) for zi, z in enumerate(zrls[1]) if z == 'Power' and zrls[3][zi] == 'PV power (W)']

        if pows:
            ap_dict = {}
            apz_dict = {}

            for pn in set(zpowrls[2]):
                ap_dict[pn] = []

                for pl in pows:
                    if pl[0] == pn:
                        ap_dict[pn].append(sum(pl[1]))

                for ap in ap_dict:
                    areslists.append(['All', 'Power', ap, 'Power (kWh)',
                                      ' '.join([str(p * 0.001) for p in ap_dict[ap]])])

            for pzn in set(['_'.join(ap.split('_')[:-1]) for ap in ap_dict]):
                ap_lists = [ap_dict[ap] for ap in ap_dict if '_'.join(ap.split('_')[:-1]) == pzn]
                apz_dict[pzn] = nsum(array(ap) * 0.001 for ap in ap_lists)

                areslists.append(['All', 'Power', pzn, 'Total power (kWh)',
                                  ' '.join([str(p) for p in apz_dict[pzn]])])

        node['envires'] = {'Invalid object': []}
    else:
        node['envires'] = node['envires{}'.format(frames[0])]

    node['reslists'] = reslists + areslists

    if node.outputs['Results out'].links:
        node.outputs['Results out'].links[0].to_node.update()


def zrupdate(self, context):
    try:
        rl = self.links[0].from_node['reslists']
        zri = [(zr[3], zr[3], 'Plot {}'.format(zr[3])) for zr in rl if zr[1] == self.rtypemenu and zr[2] == self.zonemenu and zr[0] == self.framemenu] if self.node.parametricmenu == '0' else [(zr[3], zr[3], 'Plot {}'.format(zr[3])) for zr in rl if zr[0] == 'All' and zr[1] == self.rtypemenu and zr[2] == self.zonemenu]
        return zri
    except Exception as e:
        print(e)
        return []

def ecrupdate(self, context):
    try:
        rl = self.links[0].from_node['reslists']
        zri = [(zr[3], zr[3], 'Plot {}'.format(zr[3])) for zr in rl if zr[1] == self.rtypemenu and zr[2] == self.ecmenu and zr[0] == self.framemenu] if self.node.parametricmenu == '0' else [(zr[3], zr[3], 'Plot {}'.format(zr[3])) for zr in rl if zr[0] == 'All' and zr[1] == self.rtypemenu and zr[2] == self.ecmenu]
        return zri
    except Exception as e:
        print(e)
        return []

def retmenu(dnode, axis, mtype):
    if mtype == 'Climate':
        return ['', dnode.inputs[axis].climmenu]
    if mtype in ('Zone spatial', 'Zone temporal'):
        return [dnode.inputs[axis].zonemenu, dnode.inputs[axis].zonermenu]
    elif mtype == 'Linkage':
        return [dnode.inputs[axis].linkmenu, dnode.inputs[axis].linkrmenu]
    elif mtype == 'External':
        return [dnode.inputs[axis].enmenu, dnode.inputs[axis].enrmenu]
    elif mtype == 'Chimney':
        return [dnode.inputs[axis].chimmenu, dnode.inputs[axis].chimrmenu]
    elif mtype == 'Position':
        return [dnode.inputs[axis].posmenu, dnode.inputs[axis].posrmenu]
    elif mtype == 'Camera':
        return [dnode.inputs[axis].cammenu, dnode.inputs[axis].camrmenu]
    elif mtype == 'Frames':
        return ['', 'Frames']
    elif mtype == 'Power':
        return [dnode.inputs[axis].powmenu, dnode.inputs[axis].powrmenu]
    elif mtype == 'Probe':
        return [dnode.inputs[axis].probemenu, dnode.inputs[axis].probermenu]
    if mtype == 'Embodied carbon':
        return [dnode.inputs[axis].ecmenu, dnode.inputs[axis].ecrmenu]

def write_ec(scene, frames, coll, reslists):
    for frame in frames:
        scene.frame_set(frame)
        mat_dict = {}
        zone_dict = {}
        fa = coll.vi_params['enparams']['floorarea'][str(frame)]

        for chil in coll.children:
            try:
                chil_fa = chil.vi_params['enparams']['floorarea'][str(frame)]

                for ob in chil.objects:
                    if chil.name not in zone_dict:
                        zone_dict[chil.name] = {}
                        zone_dict[chil.name]['ec'] = 0
                        zone_dict[chil.name]['ecy'] = 0
                        zone_dict[chil.name]['area'] = 0

                    for mat in ob.data.materials:
                        con_node = get_con_node(mat.vi_params)

                        if mat.name not in mat_dict:
                            mat_dict[mat.name] = {}
                            mat_dict[mat.name]['area'] = 0
                            mat_dict[mat.name]['ec'] = 0
                            mat_dict[mat.name]['ecy'] = 0

                        for poly in ob.data.polygons:
                            if ob.material_slots[poly.material_index].material == mat:
                                mat_dict[mat.name]['area'] += poly.area
                                zone_dict[chil.name]['area'] += poly.area
                                mat_ec = con_node.ret_ec()

                                if mat_ec[0] !='N/A':
                                    mat_dict[mat.name]['ec'] += float(mat_ec[0]) * poly.area
                                    zone_dict[chil.name]['ec'] += float(mat_ec[0]) * poly.area
                                    mat_dict[mat.name]['ecy'] += float(mat_ec[1]) * poly.area
                                    zone_dict[chil.name]['ecy'] += float(mat_ec[1]) * poly.area
            except:
                pass

        for zone in zone_dict:
            reslists.append([str(frame), 'Embodied carbon', zone, 'Zone EC (kgCO2e/y)', '{:.3f}'.format(zone_dict[zone]['ecy'])])
            reslists.append([str(frame), 'Embodied carbon', zone, 'Zone EC (kgCO2e/m2/y)', '{:.3f}'.format(zone_dict[zone]['ecy']/chil_fa)])
            reslists.append([str(frame), 'Embodied carbon', zone, 'Surface area (m2)', '{:.3f}'.format(zone_dict[zone]['area'])])

        reslists.append([str(frame), 'Embodied carbon', 'All', 'Total EC (kgCO2e/y)', '{:.3f}'.format(sum([zone_dict[zone]['ecy'] for zone in zone_dict]))])
        reslists.append([str(frame), 'Embodied carbon', 'All', 'Total surface area (m2)', '{:.3f}'.format(sum([mat_dict[mat]['area'] for mat in mat_dict]))])

        if fa:
            reslists.append([str(frame), 'Embodied carbon', 'All', 'Total EC (kgCO2e/m2/y)', '{:.3f}'.format(sum([zone_dict[zone]['ecy'] for zone in zone_dict])/fa)])

        for mat in mat_dict:
            reslists.append([str(frame), 'Embodied carbon', mat, 'Surface area (m2)', '{:.3f}'.format(mat_dict[mat]['area'])])
            reslists.append([str(frame), 'Embodied carbon', mat, 'Surface EC (kgCO2e/y)', '{:.3f}'.format(mat_dict[mat]['ecy'])])

            if fa:
                reslists.append([str(frame), 'Embodied carbon', mat, 'Surface EC (kgCO2e/m2/y)', '{:.3f}'.format(mat_dict[mat]['ecy']/fa)])

    if len(frames) > 1:
        for zone in zone_dict:
            reslists.append(['All', 'Embodied carbon', zone, 'Zone EC (kgCO2e/y)', ' '.join([ec[4] for ec in reslists if ec[2] == zone and ec[3] == 'Zone EC (kgCO2e/y)'])])
        for mat in mat_dict:
            reslists.append(['All', 'Embodied carbon', mat, 'Surface EC (kgCO2e/y)', ' '.join([ec[4] for ec in reslists if ec[2] == mat and ec[3] == 'Surface EC (kgCO2e/y)'])])
    scene.frame_set(frames[0])
    return (reslists)


# def sunposenvi(scene, sun, dirsol, difsol, mdata, ddata, hdata):
#     frames = range(scene.frame_start, scene.frame_end)
#     times = [datetime.datetime(2015, mdata[hi], ddata[hi], h - 1, 0) for hi, h in enumerate(hdata)]
#     solposs = [solarPosition(time.timetuple()[7], time.hour + (time.minute)*0.016666, scene.latitude, scene.longitude) for time in times]
#     beamvals = [0.01 * d for d in dirsol]
#     skyvals =  [1 + 0.01 * d for d in difsol]
#     sizevals = [beamvals[t]/skyvals[t] for t in range(len(times))]
#     values = list(zip(sizevals, beamvals, skyvals))
#     sunapply(scene, sun, values, solposs, frames)

        #             if hdict[k][0] == 'Power' and hdict[k][2] == 'PV power (W)':
        #                 pzn = '_'.join(hdict[k][1].split('_')[:-1])
        #                 if pzn not in pow_dict:
        #                     pow_dict[pzn] = array([float(p) for p in bdict[k].split(' ')])
        #                 else:
        #                     pow_dict[pzn] += array([float(p) for p in bdict[k].split(' ')])

        #             elif hdict[k][0] == 'Power' and hdict[k][2] == 'PV energy (J)':
        #                 pzn = '_'.join(hdict[k][1].split('_')[:-1])
        #                 if pzn not in en_dict:
        #                     en_dict[pzn] = array([float(p) for p in bdict[k].split(' ')])
        #                 else:
        #                     en_dict[pzn] += array([float(p) for p in bdict[k].split(' ')])

        #             elif hdict[k][0] == 'Power' and hdict[k][2] == 'PV efficiency (%)':
        #                 pzn = '_'.join(hdict[k][1].split('_')[:-1])
        #                 if pzn not in ef_dict:
        #                     ef_dict[pzn] = array([float(p) for p in bdict[k].split(' ')])
        #                 else:
        #                     ef_dict[pzn] += array([float(p) for p in bdict[k].split(' ')])

        #             elif hdict[k][0] == 'Power' and hdict[k][2] == 'PV efficiency (%)':
        #                 pzn = '_'.join(hdict[k][1].split('_')[:-1])
        #                 if pzn not in ef_dict:
        #                     ef_dict[pzn] = array([float(p) for p in bdict[k].split(' ')])
        #                 else:
        #                     ef_dict[pzn] += array([float(p) for p in bdict[k].split(' ')])

        #     for pd in pow_dict:
        #         reslists.append([str(frame), 'Power', pd, 'Total power (W)', ' '.join([str(p) for p in pow_dict[pd]])])
        #     for ed in en_dict:
        #         reslists.append([str(frame), 'Power', ed, 'Total energy (J)', ' '.join([str(e) for e in en_dict[ed]])])
        #     for ef in ef_dict:
        #         fno = len([hd for hd in hdict if '_'.join(hdict[k][1].split('_')[:-1]) == ef and hdict[k][2] == 'PV efficiency (%)'])
        #         reslists.append([str(frame), 'Power', ef, 'Average efficiency (%)', ' '.join([str(e) for e in ef_dict[ef]/fno])])
        # print(pow_dict)


# def retdata(dnode, axis, mtype, resdict, frame):
#     if mtype == 'Climate':
#         return resdict[frame][mtype][dnode.inputs[axis].climmenu]
#     if mtype == ('Zone spatial', 'Zone temporal', 'Embodied carbon'):
#         return resdict[frame][mtype][dnode.inputs[axis].zonemenu][dnode.inputs[axis].zonermenu]
#     elif mtype == 'Linkage':
#         return resdict[frame][mtype][dnode.inputs[axis].linkmenu][dnode.inputs[axis].linkrmenu]
#     elif mtype == 'External':
#         return resdict[frame][mtype][dnode.inputs[axis].enmenu][dnode.inputs[axis].enrmenu]
#     elif mtype == 'Chimney':
#         return resdict[frame][mtype][dnode.inputs[axis].chimmenu][dnode.inputs[axis].chimrmenu]
#     elif mtype == 'Position':
#         return resdict[frame][mtype][dnode.inputs[axis].posmenu][dnode.inputs[axis].posrmenu]
#     elif mtype == 'Camera':
#         return resdict[frame][mtype][dnode.inputs[axis].cammenu][dnode.inputs[axis].camrmenu]
#     elif mtype == 'Power':
#         return resdict[frame][mtype][dnode.inputs[axis].powmenu][dnode.inputs[axis].powrmenu]
#     elif mtype == 'Probe':
#         return resdict[frame][mtype][dnode.inputs[axis].probemenu][dnode.inputs[axis].probermenu]
