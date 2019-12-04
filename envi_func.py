import bpy, mathutils, colorsys, os
from collections import OrderedDict
from numpy import arange, array
from numpy import sum as nsum
from .vi_func import logentry 

def get_mat(node, ee):
    for material in bpy.data.materials:
        if node.id_data == material.vi_params.envi_nodes and material.vi_params.envi_export == ee:
            break
    return material

def get_con_node(mvp):
    if mvp.get('envi_nodes'):
#        print('names', mvp.envi_nodes.name, print(mvp.id_data.name))
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
    print('mat', mat.name)
    if mvp.get('envi_nodes'):
        print('names100', print(mat), mvp.envi_nodes.name, 'hello')
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
                bmat = bobj.data.materials[bpoly.material_index]
                
                if emnode.resist != get_con_node(bmat.vi_params).resist:
                    logentry('U-values of the paired boundary surfaces {0} and {1} do not match. {1} construction takes precedence'.format(mat.name+'_'+str(poly.index), insock.links[0].to_node.zone+'_'+str(bpoly.index)))
                    return(('', '', '', ''))
                else:
                    return(("Surface", insock.links[0].from_node.zone+'_'+str(bpoly.index), "NoSun", "NoWind"))
                
            elif outsock.links:
                bobj = bpy.data.objects[outsock.links[0].to_node.zone]
                bpoly = bobj.data.polygons[int(outsock.links[0].to_socket.name.split('_')[-2])]
                bmat = bobj.data.materials[bpoly.material_index]
                
                if emnode.resist != get_con_node(bmat.vi_params).resist: 
                    logentry('U-values of the paired boundary surfaces {0} and {1} do not match. {0} construction takes precedence'.format(mat.name+'_'+str(poly.index), outsock.links[0].to_node.zone+'_'+str(bpoly.index)))
                    return(("Zone", bobj, "NoSun", "NoWind"))
                else:
                    return(("Surface", outsock.links[0].to_node.zone+'_'+str(bpoly.index), "NoSun", "NoWind"))
                
            else:
                return(("Adiabatic", "", "NoSun", "NoWind"))

    elif emnode.envi_con_con == 'Thermal mass':
        return(("Adiabatic", "", "NoSun", "NoWind"))
    elif emnode.envi_con_con == 'Ground':
        return(("Ground", "", "NoSun", "NoWind"))
    else:
        return(("Outdoors", "", "SunExposed", "WindExposed"))

def retenresdict(scene):
    return {'Temp': ('Temperature (degC)', scene.en_temp_max, scene.en_temp_min, u"\u00b0C"), 'Hum': ('Humidity (%)', scene.en_hum_max, scene.en_hum_min, '%'),
           'CO2': ('CO2 (ppm)', scene.en_co2_max, scene.en_co2_min, 'ppm'), 'Heat': ('Heating (W)', scene.en_heat_max, scene.en_heat_min, 'W'), 'Cool': ('Cooling (W)', scene.en_cool_max, scene.en_cool_min, 'W'),
            'PMV': ('PMV', scene.en_pmv_max, scene.en_pmv_min, 'PMV'), 'PPD': ('PPD (%)', scene.en_ppd_max, scene.en_ppd_min, 'PPD'), 'SHG': ('Solar gain (W)', scene.en_ppd_max, scene.en_ppd_min, 'SHG'),
            'MaxHeat': ('Max heating (W)', scene.en_maxheat_max, scene.en_maxheat_min, 'W'), 'MaxTemp': ('Max temp (C)', scene.en_maxtemp_max, scene.en_maxtemp_min, u"\u00b0C"),
            'HRheat': ('HR heating (W)', scene.en_hrheat_max, scene.en_hrheat_min, 'hrW')}

def resnameunits():
    rnu = {'0': ("Air", "Ambient air metrics"),'1': ("Wind Speed", "Ambient Wind Speed (m/s)"), '2': ("Wind Direction", "Ambient Wind Direction (degrees from North)"),
                '3': ("Humidity", "Ambient Humidity"),'4': ("Solar", 'Ambient solar metrics'), '5': ("Temperature", "Zone Temperature"), '6': ("Humidity", "Zone Humidity"),
                '7': ("Heating Watts", "Zone Heating Requirement (Watts)"), '8': ("Cooling Watts", "Zone Cooling Requirement (Watts)"),
                '9': ("Solar Gain", "Window Solar Gain (Watts)"), '10': ("PPD", "Percentage Proportion Dissatisfied"), '11': ("PMV", "Predicted Mean Vote"),
                '12': ("Ventilation (l/s)", "Zone Ventilation rate (l/s)"), '13': (u'Ventilation (m\u00b3/h)', u'Zone Ventilation rate (m\u00b3/h)'),
                '14': (u'Infiltration (m\u00b3/hr)',  u'Zone Infiltration (m\u00b3/hr)'), '15': ('Infiltration (ACH)', 'Zone Infiltration rate (ACH)'), '16': ('CO2 (ppm)', 'Zone CO2 concentration (ppm)'),
                '17': ("Heat loss (W)", "Ventilation Heat Loss (W)"), '18': (u'Flow (m\u00b3/s)', u'Linkage flow (m\u00b3/s)'), '19': ('Opening factor', 'Linkage Opening Factor'),
                '20': ("MRT (K)", "Mean Radiant Temperature (K)"), '21': ('Occupancy', 'Occupancy count'), '22': ("Humidity", "Zone Humidity"),
                '23': ("Fabric HB (W)", "Fabric convective heat balance"), '24': ("Air Heating", "Zone air heating"), '25': ("Air Cooling", "Zone air cooling"),
                '26': ("HR Heating", "Heat recovery heating (W)"), '27': ("Volume flow", "Thermal chimney volume flow rate (m3/2)"), '28': ("Mass flow", "Thermal chmimney mass flow rate (kg/s"),
                '29': ("Out temp.", "Thermal chimney outlet temperature (C)"), '30': ("Heat loss", "Thermal chimney heat loss (W)"), '31': ("Heat gain", "Thermal chimney heat gain (W)"),
                '32': ("Volume", "Thermal chimnwey volume (m3)"), '33': ("Mass", "Thermal chimney mass (kg)"), '34': ('delta P', 'Linkage Pressure Differential (Pa)'),
                '35': ('Equipment', 'Other equipment heat gains (W)'), '36': ('PV Energy', 'PV energy (J)'), '37': ('PV Power', 'PV power (W)'),
                '38': ('PV Temp.', 'PV Temperature (C)'), '39': ('PV Eff.', 'PV efficiency (%)')}

    return [bpy.props.BoolProperty(name = rnu[str(rnum)][0], description = rnu[str(rnum)][1], default = False) for rnum in range(len(rnu))]

def aresnameunits():
    rnu = {'0': (u"Max temp (\u2103)", "Maximum zone temperature"), '1': (u"Min temp (\u2103)", "Minimum zone temperature"), '2': (u"Avg temp (\u2103)", "Average zone temperature"), 
                '3': ("Max heating (W)", "Max Zone heating"), '4': ("Min heating (W)", "Min Zone heating"), '5': ("Avg heating (W)", "Avg Zone heating"), 
                '6': ("Total heating (kWh)", "Total zone heating"), '7': (u"Total heating (kWh/m\u00b2)", "Total zone heating per floor area"),
                '8': ("Max cooling (W)", "Max Zone cooling"), '9': ("Min cooling (W)", "Min Zone cooling"), '10': ("Avg cooling (W)", "Avg Zone colling"), 
                '11': ("Total cooling (kWh)", "Total zone cooling"), '12': (u"Total cooling (kWh/m\u00b2)", "Total zone cooling per floor area"), 
                '13': ("Max CO2 (ppm)", u"Maximum zone CO\u2082 level"), '14': ("Avg CO2 (ppm)", u"Average zone CO\u2082 level"), '15': ("Min CO2 (ppm)", u"Minimum zone CO\u2082 level"),
                '16': (u"Max flow in (m\u00b3/s)", u"Maximum linkage flow level"), '17': (u"Min flow in (m\u00b3/s)", u"Minimum linkage flow level"), '18': (u"Avg flow in (m\u00b3/s)", u"Average linkage flow level"),
                '19': ('Max SHG (W)', 'Maximum Solar Heat Gain'), '20': ('Min SHG (W)', 'Minimum Solar Heat Gain'), '21': ('Avg SHG (W)', 'Average Solar Heat Gain'),
                '22': ('Total SHG (kWh)', 'Total solar heat gain'),'23': ('Total SHG (kWh/m2)', 'Total Ssolar heat gain per floor area')}
    return [bpy.props.BoolProperty(name = rnu[str(rnum)][0], description = rnu[str(rnum)][1], default = False) for rnum in range(len(rnu))]

def enresprops(disp):
    return {'0': (0, "restt{}".format(disp), "resh{}".format(disp), 0, "restwh{}".format(disp), "restwc{}".format(disp), 0, 
                  "ressah{}".format(disp), "reshrhw{}".format(disp), 0, "ressac{}".format(disp), "reswsg{}".format(disp), 0, 
                  "resfhb{}".format(disp), "resoeg{}".format(disp)),
            '1': (0, "rescpp{}".format(disp), "rescpm{}".format(disp), 0, 'resmrt{}'.format(disp), 'resocc{}'.format(disp)), 
            '2': (0, "resim{}".format(disp), "resiach{}".format(disp), 0, "resco2{}".format(disp), "resihl{}".format(disp)), 
            '3': (0, "resl12ms{}".format(disp), "reslof{}".format(disp), 0, "resldp{}".format(disp)), 
            '4':(0, "restcvf{}".format(disp), "restcmf{}".format(disp), 0, "restcot{}".format(disp), "restchl{}".format(disp),
                 0, "restchg{}".format(disp), "restcv{}".format(disp), 0, "restcm{}".format(disp)),
            '5':(0, "respve{}".format(disp), "respvw{}".format(disp), 0, "respveff{}".format(disp), "respvt{}".format(disp))}

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
                    resob.data.shape_keys.key_blocks[1].keyframe_insert(data_path = 'value', frame = frame)
            
            if resob.data.shape_keys and resnode['resdict'][rd][1] == 'Linkage Flow in':
                bpy.ops.mesh.primitive_cone_add()
                fcone = bpy.context.active_object
                fcone.rotation_euler = resob.rotation_euler if scene.envi_flink else mathutils.angle(fcone.matrix_world * fcone.data.polygons[-1].normal, resob.matrix_word * resob.data.polygons[int(resnode['resdict'][rd][0].split('_')[-1])].normal)
                fcone.parent = resob
                fcone['envires'] = {}
                fi = resnode['allresdict'][rd]
                
                for frd in resnode['resdict']:
                    if resnode['resdict'][frd][0] == resnode['resdict'][rd][0] and resnode['resdict'][frd][1] == 'Linkage Flow out':
                        fo = resnode['allresdict'][frd]
                fcone['envires']['flow'] = [float(fival) - float(foval) for fival, foval in zip(fi,fo)]
                
                for frame in range(scene.frame_start, scene.frame_end + 1):
                    scene.frame_set(frame)
                    fcone.rotation_euler = fcone.rotation_euler.to_matrix().inverted().to_euler()
                    fcone.scale = [10*float(fcone['envires']['flow'][frame]) for i in range(3)]
                    fcone.keyframe_insert(data_path = 'scale', frame = frame)
                    fcone.keyframe_insert(data_path = 'rotation_euler', frame = frame)

def envizres(scene, eresobs, resnode, restype):
    rl = resnode['reslists']
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
    
    maxval = max([[max(float(r) for r in zrl[4][ri].split())][0] for ri, r in enumerate(zrl[3]) if r == resdict[restype][0] and zrl[1][ri] == 'Zone']) 
    minval = min([[min(float(r) for r in zrl[4][ri].split())][0] for ri, r in enumerate(zrl[3]) if r == resdict[restype][0] and zrl[1][ri] == 'Zone'])

    for eo in eresobs:
        o = bpy.data.objects[eo[3:]]        
        opos = o.matrix_world * mathutils.Vector([sum(ops)/8 for ops in zip(*o.bound_box)])

        if not any([oc['VIType'] == 'envi_{}'.format(restype.lower()) for oc in o.children if oc.get('VIType')]):    
            if scene.en_disp == '1':
                bpy.ops.mesh.primitive_plane_add()  
            elif scene.en_disp == '0':
                bpy.ops.mesh.primitive_circle_add(fill_type = 'NGON')   
            ores = bpy.context.active_object
            ores['VIType'] = 'envi_{}'.format(restype.lower())

            for rs, resstring in enumerate(resstrings):
                valstring = [r[4].split()[resstart:resend] for r in rl if r[0] == frames[rs] and r[2] == eo.upper() and r[3] == resdict[restype][0]]
                vals = [float(v) for v in valstring[0]]
                if not ores.get(resstring):
                    ores[resstring] = {}
                    ores[resstring][restype] = vals

            bpy.ops.object.editmode_toggle()
            bpy.ops.mesh.extrude_region_move(MESH_OT_extrude_region={"mirror":False}, TRANSFORM_OT_translate={"value":(0, 0, 1), "constraint_axis":(False, False, True), "constraint_orientation":'NORMAL', "mirror":False, "proportional":'DISABLED', "proportional_edit_falloff":'SMOOTH', "proportional_size":1, "snap":False, "snap_target":'CLOSEST', "snap_point":(0, 0, 0), "snap_align":False, "snap_normal":(0, 0, 0), "gpencil_strokes":False, "texture_space":False, "remove_on_cancel":False, "release_confirm":False})
            bpy.ops.object.editmode_toggle()
            ores.scale, ores.parent = (0.25, 0.25, 0.25), o
            ores.location = o.matrix_world.inverted() * opos
            bpy.ops.object.material_slot_add()
            mat = bpy.data.materials.new(name = '{}_{}'.format(o.name, restype.lower()))
            ores.material_slots[0].material = mat 
            bpy.ops.object.text_add(radius=1, view_align=False, enter_editmode=False, layers=(True, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False))
            txt = bpy.context.active_object
            bpy.context.object.data.extrude = 0.005
            bpy.ops.object.material_slot_add()
            txt.parent = ores
            txt.location, txt.scale = (0,0,0), (ores.scale[0]*2, ores.scale[1]*2, 1)
            txt.data.align_x,txt.data.align_y  = 'CENTER', 'CENTER'
            txt.name = '{}_{}_text'.format(o.name, restype)
            tmat = bpy.data.materials.new(name = '{}'.format(txt.name))
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
                txt.parent, txt.location, txt.scale, txt.name = ores, (0,0,0), (ores.scale[0], ores.scale[1], 1), '{}_{}_text'.format(o.name, restype)
                txt.data.align_x,txt.data.align_y  = 'CENTER', 'CENTER'
                tmat = bpy.data.materials.new(name = '{}'.format(txt.name))
                tmat.diffuse_color = (0, 0, 0)
                txt.material_slots[0].material = tmat
            else:
                txt = ores.children[0]
        txt.data.body = "{:.1f}{}".format(ores[resstring][restype][0], resdict[restype][3]) if restype not in ('SHG', 'CO2') else "{:.0f}{}".format(ores[resstring][restype][0], resdict[restype][2])
      
        if maxval - minval:
            scalevel =  [(vals[frame] - minval)/(maxval - minval) for frame in range(0, len(vals))] if maxval - minval else [0] * len(vals)            
            colval = [colorsys.hsv_to_rgb(0.667 * (maxval - vals[vi])/(maxval - minval), 1, 1) for vi in range(len(vals))]
        else:
            scalevel = colval = [0] * len(vals)

        sv = [(sv, 0.1)[sv <= 0.1] for sv in scalevel]    
        cv = [(((0, 1)[vals[c] >= maxval], 0, (0, 1)[vals[c] <= minval]), cv)[minval < vals[c] < maxval] for c, cv in enumerate(colval)]    
        ores.animation_data_clear()
        ores.animation_data_create()
        ores['max'], ores['min'], ores['cmap'], ores['days'], ores['hours'] = maxval, minval, scene.vi_leg_col, o['days'], o['hours']
        ores.animation_data.action = bpy.data.actions.new(name="EnVi Zone")
        oresz = ores.animation_data.action.fcurves.new(data_path="scale", index = 2)
        oresz.keyframe_points.add(len(sv))
        mat.animation_data_clear()
        mat.animation_data_create()
        mat.animation_data.action = bpy.data.actions.new(name="EnVi Zone Material")
        mdcr = mat.animation_data.action.fcurves.new(data_path="diffuse_color", index = 0)
        mdcg = mat.animation_data.action.fcurves.new(data_path="diffuse_color", index = 1)
        mdcb = mat.animation_data.action.fcurves.new(data_path="diffuse_color", index = 2)
        mdcr.keyframe_points.add(len(sv))
        mdcg.keyframe_points.add(len(sv))
        mdcb.keyframe_points.add(len(sv))
        txt.animation_data_clear()
        txt.animation_data_create()
        txt.animation_data.action = bpy.data.actions.new(name="EnVi Zone Text")
        txtl = txt.animation_data.action.fcurves.new(data_path="location", index = 2)
        txtl.keyframe_points.add(len(sv))

        for frame in range(len(sv)):
            oresz.keyframe_points[frame].co = frame, sv[frame]
            mdcr.keyframe_points[frame].co = frame, cv[frame][0]
            mdcg.keyframe_points[frame].co = frame, cv[frame][1]
            mdcb.keyframe_points[frame].co = frame, cv[frame][2]
            txtl.keyframe_points[frame].co = frame, 1

def epentry(header, params, paramvs):
    return '{}\n'.format(header+(',', '')[header == ''])+'\n'.join([('    ', '')[header == '']+'{:{width}}! - {}'.format(str(pv[0])+(',', ';')[pv[1] == params[-1]], pv[1], width = 80 + (0, 4)[header == '']) for pv in zip(paramvs, params)]) + ('\n\n', '')[header == '']

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
    except:
        return [("0", "Static", "Static results")]

#def powres():
#    pass
def retrmenus(innode, node, axis): 
    rl = innode['reslists']
    zrl = list(zip(*rl))
    ftype = [(frame, frame, "Plot "+frame) for frame in list(OrderedDict.fromkeys(zrl[0])) if frame != 'All']        
    frame = 'All' if node.parametricmenu == '1' and len(ftype) > 1 else zrl[0][0]
    if axis == 'X-axis':
        rtypes = list(OrderedDict.fromkeys([zrl[1][ri] for ri, r in enumerate(zrl[1]) if zrl[0][ri] == frame]))
    else:
        rtypes = list(OrderedDict.fromkeys([zrl[1][ri] for ri, r in enumerate(zrl[1]) if zrl[0][ri] == frame and zrl[1][ri] != 'Time']))
    rtype = [(metric, metric, "Plot " + metric) for metric in rtypes]
    ctype = [(metric, metric, "Plot " + metric) for m, metric in enumerate(zrl[3]) if zrl[1][m] == 'Climate' and zrl[0][m] == frame]
    ztypes = list(OrderedDict.fromkeys([metric for m, metric in enumerate(zrl[2]) if zrl[1][m] == 'Zone' and zrl[0][m] == frame]))
    ztype = [(metric, metric, "Plot " + metric) for metric in ztypes]
    ptypes = list(OrderedDict.fromkeys([metric for m, metric in enumerate(zrl[2]) if zrl[1][m] == 'Position' and zrl[0][m] == frame]))
    ptype = [(metric, metric, "Plot " + metric) for metric in ptypes]
    prtypes = list(OrderedDict.fromkeys([metric for m, metric in enumerate(zrl[3]) if zrl[1][m] == 'Position' and zrl[0][m] == frame]))
    prtype = [(metric, metric, "Plot " + metric) for metric in prtypes]
    camtypes = list(OrderedDict.fromkeys([metric for m, metric in enumerate(zrl[2]) if zrl[1][m] == 'Camera' and zrl[0][m] == frame]))
    camtype = [(metric, metric, "Plot " + metric) for metric in camtypes]
    camrtypes = list(OrderedDict.fromkeys([metric for m, metric in enumerate(zrl[3]) if zrl[1][m] == 'Camera' and zrl[0][m] == frame]))
    camrtype = [(metric, metric, "Plot " + metric) for metric in camrtypes]
    ltypes = list(OrderedDict.fromkeys([metric for m, metric in enumerate(zrl[2]) if zrl[1][m] == 'Linkage' and zrl[0][m] == frame]))
    ltype = [(metric, metric, "Plot " + metric) for metric in ltypes]
    lrtypes = list(OrderedDict.fromkeys([metric for m, metric in enumerate(zrl[3]) if zrl[1][m] == 'Linkage' and zrl[0][m] == frame]))
    lrtype = [(metric, metric, "Plot " + metric) for metric in lrtypes]
    entypes = list(OrderedDict.fromkeys([metric for m, metric in enumerate(zrl[2]) if zrl[1][m] == 'External' and zrl[0][m] == frame]))
    entype = [(metric, metric, "Plot " + metric) for metric in entypes]
    enrtypes = list(OrderedDict.fromkeys([metric for m, metric in enumerate(zrl[3]) if zrl[1][m] == 'External' and zrl[0][m] == frame]))       
    enrtype = [(metric, metric, "Plot " + metric) for metric in enrtypes]    
    chimtypes = list(OrderedDict.fromkeys([metric for m, metric in enumerate(zrl[2]) if zrl[1][m] == 'Chimney' and zrl[0][m] == frame]))
    chimtype = [(metric, metric, "Plot " + metric) for metric in chimtypes]
    chimrtypes = list(OrderedDict.fromkeys([metric for m, metric in enumerate(zrl[3]) if zrl[1][m] == 'Chimney' and zrl[0][m] == frame]))       
    chimrtype = [(metric, metric, "Plot " + metric) for metric in chimrtypes]   
    powtypes = list(OrderedDict.fromkeys([metric for m, metric in enumerate(zrl[2]) if zrl[1][m] == 'Power' and zrl[0][m] == frame]))
    powtype = [(metric, metric, "Plot " + metric) for metric in powtypes]
    powrtypes = list(OrderedDict.fromkeys([metric for m, metric in enumerate(zrl[3]) if zrl[1][m] == 'Power' and zrl[0][m] == frame]))
    powrtype = [(metric, metric, "Plot " + metric) for metric in powrtypes]
    fmenu = bpy.props.EnumProperty(items=ftype, name="", description="Frame number", default = ftype[0][0])
    rtypemenu = bpy.props.EnumProperty(items=rtype, name="", description="Result types", default = rtype[0][0])
    statmenu = bpy.props.EnumProperty(items=[('Average', 'Average', 'Average Value'), ('Maximum', 'Maximum', 'Maximum Value'), ('Minimum', 'Minimum', 'Minimum Value')], name="", description="Zone result", default = 'Average')
    valid = ['Vi Results']    
    climmenu = bpy.props.EnumProperty(items=ctype, name="", description="Climate type", default = ctype[0][0]) if ctype else ''     
    zonemenu = bpy.props.EnumProperty(items=ztype, name="", description="Zone", default = ztype[0][0]) if ztype else ''
    zonermenu = bpy.props.EnumProperty(items=zrupdate, name="", description="Flow linkage result")# if ztype else ''
    linkmenu = bpy.props.EnumProperty(items=ltype, name="", description="Flow linkage result", default = ltype[0][0]) if ltype else ''
    linkrmenu = bpy.props.EnumProperty(items=lrtype, name="", description="Flow linkage result", default = lrtype[0][0]) if ltype else ''
    enmenu = bpy.props.EnumProperty(items=entype, name="", description="External node result", default = entype[0][0]) if entype else ''
    enrmenu = bpy.props.EnumProperty(items=enrtype, name="", description="External node result", default = enrtype[0][0]) if entype else ''
    chimmenu = bpy.props.EnumProperty(items=chimtype, name="", description="External node result", default = chimtype[0][0]) if chimtype else ''
    chimrmenu = bpy.props.EnumProperty(items=chimrtype, name="", description="External node result", default = chimrtype[0][0]) if chimtype else ''
    posmenu =  bpy.props.EnumProperty(items=ptype, name="", description="Position result", default = ptype[0][0]) if ptype else ''
    posrmenu = bpy.props.EnumProperty(items=prtype, name="", description="Position result", default = prtype[0][0]) if ptypes else ''
    cammenu =  bpy.props.EnumProperty(items=camtype, name="", description="Camera result", default = camtype[0][0]) if camtype else ''
    camrmenu = bpy.props.EnumProperty(items=camrtype, name="", description="Camera result", default = camrtype[0][0]) if camtypes else ''
    powmenu = bpy.props.EnumProperty(items=powtype, name="", description="Power result", default = powtype[0][0]) if powtype else ''
    powrmenu = bpy.props.EnumProperty(items=powrtype, name="", description="Power result", default = powrtype[0][0]) if powrtype else ''
    multfactor = bpy.props.FloatProperty(name = "", description = "Result multiplication factor", min = -10000, max = 10000, default = 1)
    
    return (valid, fmenu, statmenu, rtypemenu, climmenu, zonemenu, zonermenu, linkmenu, linkrmenu, enmenu, enrmenu, chimmenu, chimrmenu, posmenu, posrmenu, cammenu, camrmenu, powmenu, powrmenu, multfactor)

def processh(lines, znlist):    
    envdict = {'Site Outdoor Air Drybulb Temperature [C] !Hourly': "Temperature (degC)",
               'Site Outdoor Air Relative Humidity [%] !Hourly': 'Humidity (%)',
                'Site Wind Direction [deg] !Hourly': 'Wind Direction (deg)',
                'Site Wind Speed [m/s] !Hourly': 'Wind Speed (m/s)',
                'Site Diffuse Solar Radiation Rate per Area [W/m2] !Hourly': "Diffuse Solar (W/m^2)",
                'Site Direct Solar Radiation Rate per Area [W/m2] !Hourly': "Direct Solar (W/m^2)"}
    zresdict = {'Zone Air Temperature [C] !Hourly': "Temperature (degC)",
                'Zone Air Relative Humidity [%] !Hourly': 'Humidity (%)',
                'Zone Air System Sensible Heating Rate [W] !Hourly': 'Heating (W)',
                'Zone Air System Sensible Cooling Rate [W] !Hourly': 'Cooling (W)',
                'Zone Ideal Loads Supply Air Sensible Heating Rate [W] !Hourly': 'Air heating (W)',
                'Zone Ideal Loads Heat Recovery Sensible Heating Rate [W] !Hourly': 'HR heating (W)',
                'Zone Ideal Loads Supply Air Sensible Cooling Rate [W] !Hourly': 'Air cooling (W)',
                'Zone Windows Total Transmitted Solar Radiation Rate [W] !Hourly': 'Solar gain (W)',
                'AFN Zone Infiltration Volume [m3] !Hourly': 'Infiltration (m3/hr)',
                'AFN Zone Infiltration Air Change Rate [ach] !Hourly': 'Infiltration (ACH)',
                'Zone Infiltration Current Density Volume [m3] !Hourly': 'Infiltration (m3/hr)',
                'Zone Infiltration Air Change Rate [ach] !Hourly': 'Infiltration (ACH)',
                'Zone Mean Air Temperature [C] ! Hourly': 'Mean Temperature (degC)', 
                'Zone Thermal Comfort Fanger Model PPD [%] !Hourly' :'PPD (%)',
                'Zone Thermal Comfort Fanger Model PMV [] !Hourly' :'PMV',               
                'AFN Node CO2 Concentration [ppm] !Hourly': 'CO2 (ppm)',
                'Zone Air CO2 Concentration [ppm] !Hourly': 'CO2 (ppm)',
                'Zone Mean Radiant Temperature [C] !Hourly': 'MRT', 
                'Zone People Occupant Count [] !Hourly': 'Occupancy', 
                'Zone Air Heat Balance Surface Convection Rate [W] !Hourly': 'Heat balance (W)',
                'Zone Thermal Chimney Current Density Air Volume Flow Rate [m3/s] !Hourly': 'Volume flow (m3/s)', 
                'Zone Thermal Chimney Mass Flow Rate [kg/s] !Hourly': 'Mass flow (kg/s)',
                'Zone Thermal Chimney Outlet Temperature [C] !Hourly': 'Outlet temperature (C)',
                'Zone Thermal Chimney Heat Loss Energy [J] !Hourly': 'TC heat loss (J)',
                'Zone Thermal Chimney Heat Gain Energy [J] !Hourly': 'TC heat gain (J)',
                'Zone Thermal Chimney Volume [m3] !Hourly': 'TC VOLUME (m3)',
                'Zone Thermal Chimney Mass [kg] !Hourly':'TC mass(kg)',
                'Zone Other Equipment Total Heating Rate [W] !Hourly': 'Equipment (W)'}
    enresdict = {'AFN Node CO2 Concentration [ppm] !Hourly': 'CO2 (ppm)',
                 'AFN Node Wind Pressure [Pa] !Hourly': 'WP (PA)'}
    lresdict = {'AFN Linkage Node 1 to Node 2 Volume Flow Rate [m3/s] !Hourly': 'Linkage Flow out',
                'AFN Linkage Node 2 to Node 1 Volume Flow Rate [m3/s] !Hourly': 'Linkage Flow in',
                'AFN Surface Venting Window or Door Opening Factor [] !Hourly': 'Opening Factor',
                'AFN Linkage Node 1 to Node 2 Pressure Difference [Pa] !Hourly': 'delta P (Pa)'}
    presdict = {'Generator Produced DC Electric Energy [J] !Hourly': 'PV Energy (J)',
                'Generator Produced DC Electric Power [W] !Hourly': 'PV Power (W)',
                'Generator PV Array Efficiency [] !Hourly': 'PV Efficiency (%)',
                'Generator PV Cell Temperature [C] !Hourly': 'PV Temperature (C)'}
    hdict = {}
    
    for l, line in enumerate(lines):
        linesplit = line.strip('\n').split(',')
        if len(linesplit) > 3:
            if linesplit[2] == 'Day of Simulation[]':
                hdict[linesplit[0]] = ['Time'] 
            elif linesplit[3] in envdict:
                hdict[linesplit[0]] = ['Climate',  '', envdict[linesplit[3]]]  
            elif linesplit[3] in zresdict and retzonename(linesplit[2]) in znlist:
                hdict[linesplit[0]] = ['Zone',  retzonename(linesplit[2]),  zresdict[linesplit[3]]]
            elif linesplit[3] in enresdict and 'ExtNode' in linesplit[2]:
                hdict[linesplit[0]] = ['External',  linesplit[2],  enresdict[linesplit[3]]]
            elif linesplit[3] in lresdict:
                hdict[linesplit[0]] = ['Linkage',  linesplit[2],  lresdict[linesplit[3]]]
            elif linesplit[3] in presdict:  
                hdict[linesplit[0]] = ['Power',  linesplit[2],  presdict[linesplit[3]]]
            print(hdict)
        if line == 'End of Data Dictionary\n':
            break
    return hdict,  l + 1
    
def retzonename(zn):
    if  zn[-10:] == '_OCCUPANCY':
        return zn.strip('_OCCUPANCY')
    elif zn[-4:] == '_AIR':
        return zn.strip('_AIR')
    else:
        return zn

def checkenvierrors(file, sim_op):
    efile = file.read()
    if '** Severe  **' in efile:
        sim_op.report({'ERROR'}, "There is a fatal error in the EnVi model, check the error file in Blender's text editor")
        
def processf(pro_op, node):
    scene = bpy.context.scene
    svp = scene.vi_params
    reslists, areslists = [], []
    frames = range(svp['enparams']['fs'], svp['enparams']['fe'] + 1) if node.bl_label == 'EnVi Simulation' else [scene.frame_current] 
    
    for frame in frames:
        node['envires{}'.format(frame)] = {}
        resfileloc = os.path.join(svp['viparams']['newdir'], '{}{}out.eso'.format(pro_op.resname, frame)) if node.bl_label == 'EnVi Simulation' else node.resfilename

        with open(resfileloc, 'r') as resfile:
            lines = resfile.readlines()
            hdict, lstart = processh(lines, [coll.name.upper() for coll in bpy.data.collections['EnVi Geometry'].children])  
            splitlines = [l.strip('\n').split(',') for l in lines[lstart:-2]]
            bdict = {li: ' '.join(['{:.3f}'.format(float(sl[1])) for sl in splitlines if sl[0] == li]) for li in hdict}
  
            for k in sorted(hdict.keys(), key=int):
                if hdict[k] == ['Time']:
                    reslists.append([str(frame), 'Time', '', 'Month', ' '.join([sl[2] for sl in splitlines if sl[0] == k])])
                    reslists.append([str(frame), 'Time', '', 'Day', ' '.join([sl[3] for sl in splitlines if sl[0] == k])])                    
                    reslists.append([str(frame), 'Time', '', 'Hour', ' '.join([sl[5] for sl in splitlines if sl[0] == k])])
                    reslists.append([str(frame), 'Time', '', 'DOS', ' '.join([sl[1] for sl in splitlines if sl[0] == k])])
                else:
                    reslists.append([str(frame)] + hdict[k] + [bdict[k]])
        
    rls = reslists
    zrls = list(zip(*rls))
    svp['enparams']['lmetrics'] = list(set([zr for zri, zr in enumerate(zrls[3]) if zrls[1][zri] == 'Linkage' and zrls[0][zri] == str(node["AStart"])]))
    svp['enparams']['zmetrics'] = list(set([zr for zri, zr in enumerate(zrls[3]) if zrls[1][zri] == 'Zone' and zrls[0][zri] == str(node["AStart"])]))

    for frame in frames:
        zonerls = [zonerl for zonerl in rls if zonerl[1] == 'Zone' and zonerl[0] == str(frame)]
        zzonerls = list(zip(*zonerls))

        try:
            for resname in set(zzonerls[3]):
                hczres = [zres[4].split() for zres in zonerls if zres[0] == str(frame) and zres[3] == resname]
                node['envires{}'.format(frame)][resname] = nsum(array([array([float(zr) for zr in hc]) for hc in hczres]), axis = 0)
            node['hours'] = arange(1, 25, dtype = float)
            node['days'] = arange(node.dsdoy, node.dedoy + 1, dtype = float) 
        except:
            pro_op.report({'ERROR'}, "There are no results to plot. Make sure you have selected valid metrics to calculate and try re-exporting/simulating")
        
        for o in scene.objects:
            if 'EN_' + o.name.upper() in zrls[2]:
                envires = {}
                oress = [[zrls[3][z], zrls[4][z]]  for z, zr in enumerate(zrls[0]) if zr == str(frame) and zrls[2][z] == 'EN_' + o.name.upper()]
                for ores in oress:                    
                    envires[ores[0]] = array([float(val) for val in ores[1].split()])
                if frame == frames[0]:
                    o['hours'] = arange(1, 25, dtype = float)
                    o['days'] = arange(node.dsdoy, node.dedoy + 1, dtype = float)
                o['envires{}'.format((frame, '')[len(frames) == 1])] = envires
                        
    if len(frames) > 1:  
        areslists = []
        areslists.append(['All', 'Frames', '', 'Frames', ' '.join([str(f) for f in frames])])
        temps = [(zrls[2][zi], [float(t) for t in zrls[4][zi].split()]) for zi, z in enumerate(zrls[1]) if z == 'Zone' and zrls[3][zi] == 'Temperature (degC)']
        hums = [(zrls[2][zi], [float(t) for t in zrls[4][zi].split()]) for zi, z in enumerate(zrls[1]) if z == 'Zone' and zrls[3][zi] == 'Humidity (%)']
        heats = [(zrls[2][zi], [float(t) for t in zrls[4][zi].split()]) for zi, z in enumerate(zrls[1]) if z == 'Zone' and zrls[3][zi] == 'Heating (W)']
        cools = [(zrls[2][zi], [float(t) for t in zrls[4][zi].split()]) for zi, z in enumerate(zrls[1]) if z == 'Zone' and zrls[3][zi] == 'Cooling (W)']
        aheats = [(zrls[2][zi], [float(t) for t in zrls[4][zi].split()]) for zi, z in enumerate(zrls[1]) if z == 'Zone' and zrls[3][zi] == 'Air Heating (W)']
        acools = [(zrls[2][zi], [float(t) for t in zrls[4][zi].split()]) for zi, z in enumerate(zrls[1]) if z == 'Zone' and zrls[3][zi] == 'Air Cooling (W)']
        co2s = [(zrls[2][zi], [float(t) for t in zrls[4][zi].split()]) for zi, z in enumerate(zrls[1]) if z == 'Zone' and zrls[3][zi] == 'CO2 (ppm)']
        comfppds = [(zrls[2][zi], [float(t) for t in zrls[4][zi].split()]) for zi, z in enumerate(zrls[1]) if z == 'Zone' and zrls[3][zi] == 'PPD']
        comfpmvs = [(zrls[2][zi], [float(t) for t in zrls[4][zi].split()]) for zi, z in enumerate(zrls[1]) if z == 'Zone' and zrls[3][zi] == 'PMV']
        shgs = [(zrls[2][zi], [float(t) for t in zrls[4][zi].split()]) for zi, z in enumerate(zrls[1]) if z == 'Zone' and zrls[3][zi] == 'Solar gain (W)']
               
        for zn in set(zzonerls[2]):
            if temps:
                areslists.append(['All', 'Zone', zn, 'Max temp (C)', ' '.join([str(max(t[1])) for t in temps if t[0] == zn])])
                areslists.append(['All', 'Zone', zn, 'Min temp (C)', ' '.join([str(min(t[1])) for t in temps if t[0] == zn])])
                areslists.append(['All', 'Zone', zn, 'Avg temp (C)', ' '.join([str(sum(t[1])/len(t[1])) for t in temps if t[0] == zn])])
            if hums:
                areslists.append(['All', 'Zone', zn, 'Max humidity (C)', ' '.join([str(max(h[1])) for h in hums if h[0] == zn])])
                areslists.append(['All', 'Zone', zn, 'Min humidity (C)', ' '.join([str(min(h[1])) for h in hums if h[0] == zn])])
                areslists.append(['All', 'Zone', zn, 'Avg humidity (C)', ' '.join([str(sum(h[1])/len(h[1])) for h in hums if h[0] == zn])])
            if heats:
                areslists.append(['All', 'Zone', zn, 'Max heating (W)', ' '.join([str(max(h[1])) for h in heats if h[0] == zn])])
                areslists.append(['All', 'Zone', zn, 'Min heating (W)', ' '.join([str(min(h[1])) for h in heats if h[0] == zn])])
                areslists.append(['All', 'Zone', zn, 'Avg heating (W)', ' '.join([str(sum(h[1])/len(h[1])) for h in heats if h[0] == zn])])
                areslists.append(['All', 'Zone', zn, 'Total heating (kWh)', ' '.join([str(sum(h[1])*0.001) for h in heats if h[0] == zn])])
                areslists.append(['All', 'Zone', zn, 'Total heating (kWh/m2)', ' '.join([str(sum(h[1])*0.001/[o for o in bpy.data.objects if o.name.upper() == zn][0]['floorarea']) for h in heats if h[0] == zn])])
            if cools:
                areslists.append(['All', 'Zone', zn, 'Max cooling (W)', ' '.join([str(max(h[1])) for h in cools if h[0] == zn])])
                areslists.append(['All', 'Zone', zn, 'Min cooling (W)', ' '.join([str(min(h[1])) for h in cools if h[0] == zn])])
                areslists.append(['All', 'Zone', zn, 'Avg cooling (W)', ' '.join([str(sum(h[1])/len(h[1])) for h in cools if h[0] == zn])])
                areslists.append(['All', 'Zone', zn, 'Total cooling (kWh)', ' '.join([str(sum(h[1])*0.001) for h in cools if h[0] == zn])])
                areslists.append(['All', 'Zone', zn, 'Total cooling (kWh/m2)', ' '.join([str(sum(h[1])*0.001/[o for o in bpy.data.objects if o.name.upper() == zn][0]['floorarea']) for h in cools if h[0] == zn])])
            if aheats:
                areslists.append(['All', 'Zone', zn, 'Max air heating (W)', ' '.join([str(max(h[1])) for h in aheats if h[0] == zn])])
                areslists.append(['All', 'Zone', zn, 'Min air heating (W)', ' '.join([str(min(h[1])) for h in aheats if h[0] == zn])])
                areslists.append(['All', 'Zone', zn, 'Avg air heating (W)', ' '.join([str(sum(h[1])/len(h[1])) for h in aheats if h[0] == zn])])
                areslists.append(['All', 'Zone', zn, 'Total air heating (kWh)', ' '.join([str(sum(h[1])*0.001) for h in aheats if h[0] == zn])])
                areslists.append(['All', 'Zone', zn, 'Total air heating (kWh/m2)', ' '.join([str(sum(h[1])*0.001/[o for o in bpy.data.objects if o.name.upper() == zn][0]['floorarea']) for h in aheats if h[0] == zn])])
            if acools:
                areslists.append(['All', 'Zone', zn, 'Max air cool (W)', ' '.join([str(max(h[1])) for h in acools if h[0] == zn])])
                areslists.append(['All', 'Zone', zn, 'Min air cool (W)', ' '.join([str(min(h[1])) for h in acools if h[0] == zn])])
                areslists.append(['All', 'Zone', zn, 'Avg air cool (W)', ' '.join([str(sum(h[1])/len(h[1])) for h in acools if h[0] == zn])])
                areslists.append(['All', 'Zone', zn, 'Air cooling (kWh)', ' '.join([str(sum(h[1])*0.001) for h in acools if h[0] == zn])])
                areslists.append(['All', 'Zone', zn, 'Air cooling (kWh/m2)', ' '.join([str(sum(h[1])*0.001/[o for o in bpy.data.objects if o.name.upper() == zn][0]['floorarea']) for h in acools if h[0] == zn])])
            if co2s:
                areslists.append(['All', 'Zone', zn, 'Max CO2 (ppm)', ' '.join([str(max(t[1])) for t in co2s if t[0] == zn])])
                areslists.append(['All', 'Zone', zn, 'Min CO2 (ppm)', ' '.join([str(min(t[1])) for t in co2s if t[0] == zn])])
                areslists.append(['All', 'Zone', zn, 'Avg CO2 (ppm)', ' '.join([str(sum(t[1])/len(t[1])) for t in co2s if t[0] == zn])])
            if comfppds:
                areslists.append(['All', 'Zone', zn, 'Max PPD', ' '.join([str(max(t[1])) for t in comfppds if t[0] == zn])])
                areslists.append(['All', 'Zone', zn, 'Min PPD', ' '.join([str(min(t[1])) for t in comfppds if t[0] == zn])])
                areslists.append(['All', 'Zone', zn, 'Avg PPD', ' '.join([str(sum(t[1])/len(t[1])) for t in comfppds if t[0] == zn])])
            if comfpmvs:
                areslists.append(['All', 'Zone', zn, 'Max PMV', ' '.join([str(max(t[1])) for t in comfpmvs if t[0] == zn])])
                areslists.append(['All', 'Zone', zn, 'Min PMV', ' '.join([str(min(t[1])) for t in comfpmvs if t[0] == zn])])
                areslists.append(['All', 'Zone', zn, 'Avg PMV', ' '.join([str(sum(t[1])/len(t[1])) for t in comfpmvs if t[0] == zn])])
            if shgs:
                areslists.append(['All', 'Zone', zn, 'Max SHG (W)', ' '.join([str(max(t[1])) for t in shgs if t[0] == zn])])
                areslists.append(['All', 'Zone', zn, 'Min SHG (W)', ' '.join([str(min(t[1])) for t in shgs if t[0] == zn])])
                areslists.append(['All', 'Zone', zn, 'Avg SHG (W)', ' '.join([str(sum(t[1])/len(t[1])) for t in shgs if t[0] == zn])])
                areslists.append(['All', 'Zone', zn, 'Total SHG (kWh)', ' '.join([str(sum(t[1])*0.001) for t in shgs if t[0] == zn])])
                areslists.append(['All', 'Zone', zn, 'Total SHG (kWh/m2)', ' '.join([str(sum(t[1])*0.001/[o for o in bpy.data.objects if o.name.upper() == zn][0]['floorarea']) for t in shgs if t[0] == zn])])
            if heats and cools:
                try:
                    conds = [sum(x) for x in zip(*[[sum(h[1])*0.001 for h in heats if h[0] == zn], [sum(h[1])*0.001 for h in cools if h[0] == zn]])]
                    areslists.append(['All', 'Zone', zn, 'Total conditioning (kWh)', ' '.join([str(c) for c in conds])])#join([str(sum(h[1])*0.001) for h in heats if h[0] == zn])])
                    areslists.append(['All', 'Zone', zn, 'Total conditioing (kWh/m2)', ' '.join([str(c/[o for o in bpy.data.objects if o.name.upper() == zn][0]['floorarea']) for c in conds])])
                except:
                    pass                
            if aheats and acools:
                try:
                    aconds = [sum(x) for x in zip(*[[sum(h[1])*0.001 for h in aheats if h[0] == zn], [sum(h[1])*0.001 for h in acools if h[0] == zn]])]
                    areslists.append(['All', 'Zone', zn, 'Total air conditioning (kWh)', ' '.join([str(c) for c in conds])])#join([str(sum(h[1])*0.001) for h in heats if h[0] == zn])])
                    areslists.append(['All', 'Zone', zn, 'Total air conditioing (kWh/m2)', ' '.join([str(c/[o for o in bpy.data.objects if o.name.upper() == zn][0]['floorarea']) for c in aconds])])
                except:
                    pass   
                
        for o in scene.objects:
            o['envires'] = {}
            for arl in areslists:
                if arl[1] == 'Zone' and 'EN_' + o.name.upper() == arl[2]:
                    o['envires'][arl[3]] = [float(val) for val in arl[4].split()]
                        
        node['envires'] = {'Invalid object': []}
    else:
        node['envires'] = node['envires{}'.format(frames[0])]  
                                                                                          
    node['reslists'] = reslists + areslists
    
    if node.outputs['Results out'].links:
       node.outputs['Results out'].links[0].to_node.update() 

def zrupdate(self, context):
    try: 
        rl = self.links[0].from_node['reslists']
        zri = [(zr[3], zr[3], 'Plot {}'.format(zr[3])) for zr in rl if zr[2] == self.zonemenu and zr[0] == self.framemenu] if self.node.parametricmenu == '0' else [(zr[3], zr[3], 'Plot {}'.format(zr[3])) for zr in rl if zr[2] == self.zonemenu and zr[0] == 'All']
        return zri
    except:
        return []

def retmenu(dnode, axis, mtype):
    if mtype == 'Climate':
        return ['', dnode.inputs[axis].climmenu]
    if mtype == 'Zone':
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
        
def retdata(dnode, axis, mtype, resdict, frame):
    if mtype == 'Climate':
        return resdict[frame][mtype][dnode.inputs[axis].climmenu]
    if mtype == 'Zone':
        return resdict[frame][mtype][dnode.inputs[axis].zonemenu][dnode.inputs[axis].zonermenu]
    elif mtype == 'Linkage':
        return resdict[frame][mtype][dnode.inputs[axis].linkmenu][dnode.inputs[axis].linkrmenu]
    elif mtype == 'External':
        return resdict[frame][mtype][dnode.inputs[axis].enmenu][dnode.inputs[axis].enrmenu]
    elif mtype == 'Chimney':
        return resdict[frame][mtype][dnode.inputs[axis].chimmenu][dnode.inputs[axis].chimrmenu]
    elif mtype == 'Position':
        return resdict[frame][mtype][dnode.inputs[axis].posmenu][dnode.inputs[axis].posrmenu]
    elif mtype == 'Camera':
        return resdict[frame][mtype][dnode.inputs[axis].cammenu][dnode.inputs[axis].camrmenu]
    elif mtype == 'Power':
        return resdict[frame][mtype][dnode.inputs[axis].powmenu][dnode.inputs[axis].powrmenu]    