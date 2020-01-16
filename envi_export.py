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

import bpy, os, itertools, subprocess, datetime, shutil, mathutils, bmesh
from .vi_func import selobj, facearea, selmesh, create_coll, selobs, logentry
from .envi_func import epentry, epschedwrite, get_con_node, boundpoly
from .envi_mat import retuval

dtdf = datetime.date.fromordinal
caidict = {"0": "", "1": "Simple", "2": "Detailed", "3": "TrombeWall", "4": "AdaptiveConvectionAlgorithm"}
caodict = {"0": "", "1": "SimpleCombined", "2": "TARP", "3": "DOE-2", "4": "MoWiTT", "5": "AdaptiveConvectionAlgorithm"}

def enpolymatexport(exp_op, node, locnode, em, ec):
    scene = bpy.context.scene 
    svp = scene.vi_params
    svp['viparams']['hvactemplate'] = 0
    depsgraph = bpy.context.evaluated_depsgraph_get()
    
    for frame in range(svp['enparams']['fs'], svp['enparams']['fe'] + 1):
        pvs = []
        gen = 0
#        scene.update()
        scene.frame_set(frame)
        geo_colls = bpy.data.collections['EnVi Geometry'].children
        en_idf = open(os.path.join(svp['viparams']['newdir'], 'in{}.idf'.format(frame)), 'w')
        enng = [ng for ng in bpy.data.node_groups if ng.bl_label == 'EnVi Network'][0]
        badnodes = [node for node in enng.nodes if node.use_custom_color]
        
        for node in badnodes:
            node.hide = 0
            exp_op.report({'ERROR'}, 'Bad {} node in the EnVi network. Delete the node if not needed or make valid connections'.format(node.name))
            return
#        if any([node.bl_idname in 'No_En_Net_SSFlow'])
        enng['enviparams']['afn'] = 1
        en_idf.write("!- Blender -> EnergyPlus\n!- Using the EnVi export scripts\n!- Author: Ryan Southall\n!- Date: {}\n\nVERSION,{};\n\n".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M"), svp['enparams']['epversion']))    
        params = ('Name', 'North Axis (deg)', 'Terrain', 'Loads Convergence Tolerance Value', 'Temperature Convergence Tolerance Value (deltaC)',
                  'Solar Distribution', 'Maximum Number of Warmup Days(from MLC TCM)')
        paramvs = ((node.loc, 'Default')[not node.loc], '0.00', ("City", "Urban", "Suburbs", "Country", "Ocean")[int(node.terrain)], '0.004', '0.4', 'FullInteriorAndExteriorWithReflections', '15')
        en_idf.write(epentry('Building', params, paramvs))
        params = ('Time Step in Hours', 'Algorithm', 'Algorithm', 'Default frequency of calculation', 'no zone sizing, system sizing, plant sizing, no design day, use weather file')
        paramvs = ('Timestep, {}'.format(node.timesteps), 'SurfaceConvectionAlgorithm:Inside, TARP', 'SurfaceConvectionAlgorithm:Outside, TARP',
                   'ShadowCalculation, AverageOverDaysInFrequency, 10', 'SimulationControl, No,No,No,No,Yes')
    
        for ppair in zip(params, paramvs):
            en_idf.write(epentry('', [ppair[0]], [ppair[1]]) + ('', '\n\n')[ppair[0] == params[-1]])

        en_idf.write('HeatBalanceAlgorithm, ConductionTransferFunction;\n\n')
   
        params = ('Name', 'Begin Month', 'Begin Day of Month', 'Begin Year', 'End Month', 'End Day of Month', 'End Year', 'Day of Week for Start Day', 'Use Weather File Holidays and Special Days', 'Use Weather File Daylight Saving Period',\
        'Apply Weekend Holiday Rule', 'Use Weather File Rain Indicators', 'Use Weather File Snow Indicators')
        paramvs = ((node.loc, 'Default')[not node.loc], node.sdate.month, node.sdate.day, '', node.edate.month, node.edate.day, '', "", "Yes", "Yes", "No", "Yes", "Yes")
        en_idf.write(epentry('RunPeriod', params, paramvs))    

        en_idf.write("!-   ===========  ALL OBJECTS IN CLASS: MATERIAL & CONSTRUCTIONS ===========\n\n")

        gen = 0        

        for mat in bpy.data.materials:   
            
            mvp = mat.vi_params            
            if  mvp.envi_nodes and mvp.envi_nodes.nodes and mvp.envi_export == True:  
                print(mat.name)
                for emnode in mvp.envi_nodes.nodes:
                    if emnode.bl_idname == 'No_En_Mat_Con' and emnode.active:
                        if emnode.envi_con_type == 'Window':    
#                            if emnode.envi_simple_glazing:
#                                em.sg_write(en_idf, mat.name+'_sg', emnode.envi_sg_uv, emnode.envi_sg_shgc, emnode.envi_sg_vt)
#                                ec.con_write(en_idf, emnode.envi_con_type, mat.name, mat.name+'_sg', mat.name, [mat.name+'_sg'])
#                            else:                            
                            en_idf.write(emnode.ep_write(mat.name))                        
                        else:                            
                            if emnode.envi_con_type not in ('None', 'Shading', 'Aperture'):
                                en_idf.write(emnode.ep_write(mat.name))
                        if emnode.inputs['PV'].links:
                            gen = 1
                            pvs.append(emnode)    
        em.namedict = {}
        em.thickdict = {}
    
        en_idf.write("!-   ===========  ALL OBJECTS IN CLASS: ZONES ===========\n\n")
        
        
        for coll in geo_colls:
#        for obj in [obj for obj in bpy.context.scene.objects if obj.layers[1] == True and obj.envi_type in ('0', '2')]:
            if coll.objects[0].vi_params.envi_type in ('0', '2'):
                params = ('Name', 'Direction of Relative North (deg)', 'X Origin (m)', 'Y Origin (m)', 'Z Origin (m)', 'Type', 'Multiplier', 'Ceiling Height (m)', 'Volume (m3)',
                          'Floor Area (m2)', 'Zone Inside Convection Algorithm', 'Zone Outside Convection Algorithm', 'Part of Total Floor Area')
                paramvs = (coll.name, 0, 0, 0, 0, 1, 1, 'autocalculate', '{:.1f}'.format(coll.objects[0]['volume']), 'autocalculate', caidict[coll.objects[0].vi_params.envi_ica], caodict[coll.objects[0].vi_params.envi_oca], 'Yes')
                en_idf.write(epentry('Zone', params, paramvs))
        
        params = ('Starting Vertex Position', 'Vertex Entry Direction', 'Coordinate System')
        paramvs = ('UpperRightCorner', 'Counterclockwise', 'World')
        en_idf.write(epentry('GlobalGeometryRules', params, paramvs))
    
        en_idf.write("!-   ===========  ALL OBJECTS IN CLASS: SURFACE DEFINITIONS ===========\n\n")
        
        wfrparams = ['Name', 'Surface Type', 'Construction Name', 'Zone Name', 'Outside Boundary Condition', 'Outside Boundary Condition Object', 'Sun Exposure', 'Wind Exposure', 'View Factor to Ground', 'Number of Vertices']
        onames = [o.name for coll in geo_colls for o in coll.objects if o.vi_params.envi_type == '0'] 
        tcnames = [o.name for coll in geo_colls for o in coll.objects if o.vi_params.envi_type == '2']
        zonenames = [c.name for c in geo_colls]
        zonenodes = [n for n in enng.nodes if hasattr(n, 'zone') and n.zone in zonenames]
        for zn in zonenodes:
            zn.update()
        tcnodes = [n for n in enng.nodes if hasattr(n, 'zone') and n.zone in tcnames]
        gens = []
        
        for coll in geo_colls:
#            zonenames = [o.name for o in coll.objects if o.vi_params.envi_type == '0']
#            tcnames = [o.name for o in coll.objects if o.vi_params.envi_type == '2']
#            zonenodes = [n for n in enng.nodes if hasattr(n, 'zone') and n.zone in zonenames]
#            for zn in zonenodes:
#                zn.update()
#            tcnodes = [n for n in enng.nodes if hasattr(n, 'zone') and n.zone in tcnames]
#            gens = []
            for obj in coll.objects:
                mats = [ms.material for ms in obj.material_slots]
#                obj = scene.objects[zn]
#            for obj in [obj for obj in coll.objects if obj.type == 'MESH' and obj.vi_type == '1']:
#                me = obj.to_mesh()
                me = obj.evaluated_get(depsgraph).to_mesh()
                bm = bmesh.new()
                bm.from_mesh(me)
                bm.transform(obj.matrix_world)
                obj.to_mesh_clear()
    
                for face in [f for f in bm.faces if mats[f.material_index].vi_params.envi_nodes]:
                    mat = mats[face.material_index]
                    for emnode in mat.vi_params.envi_nodes.nodes:
                        if emnode.bl_idname == 'No_En_Mat_Con' and emnode.active:
                            vcos = [v.co for v in face.verts]
                            (obc, obco, se, we) = boundpoly(obj, emnode, face, enng)
                            
                            if obc:
                                if emnode.envi_con_type in ('Wall', "Floor", "Roof"):
                                    if emnode.envi_con_makeup != "2":
                                        params = list(wfrparams) + ["X,Y,Z ==> Vertex {} (m)".format(v.index) for v in face.verts]                     
                                        paramvs = ['{}_{}'.format(obj.name, face.index), emnode.envi_con_type, mat.name, coll.name, obc, obco, se, we, 'autocalculate', len(face.verts)]+ ["  {0[0]:.4f}, {0[1]:.4f}, {0[2]:.4f}".format(vco) for vco in vcos]
                                        en_idf.write(epentry('BuildingSurface:Detailed', params, paramvs))
                                                                            
                                    if emnode.inputs['PV'].links:
                                        pvnode = emnode.inputs['PV'].links[0].from_node
                                        en_idf.write(pvnode.ep_write('{}_{}'.format(obj.name, face.index)))
                        #                en_idf.write(emnode.pv_ep_write('{}_{}'.format(obj.name, face.index)))
                                        gens.append('{}_{}-pv'.format(obj.name, face.index))
                                        
                                elif emnode.envi_con_type in ('Door', 'Window')  and emnode.envi_con_makeup != "2":
                                    if len(face.verts) > 4:
                                        exp_op.report({'ERROR'}, 'Window/door in {} has more than 4 vertices'.format(obj.name))
                                        
                                    xav, yav, zav = mathutils.Vector(face.calc_center_median())
                                    params = list(wfrparams) + ["X,Y,Z ==> Vertex {} (m)".format(v.index) for v in face.verts]
                                    paramvs = ['{}_{}'.format(obj.name, face.index), 'Wall', '{}-frame'.format(mat.name), coll.name, obc, obco, se, we, 'autocalculate', len(face.verts)] + ["  {0[0]:.4f}, {0[1]:.4f}, {0[2]:.4f}".format(vco) for vco in vcos]
                                    en_idf.write(epentry('BuildingSurface:Detailed', params, paramvs))    
                                    obound = ('win-', 'door-')[emnode.envi_con_type == 'Door']+obco if obco else obco
                                    params = ['Name', 'Surface Type', 'Construction Name', 'Building Surface Name', 'Outside Boundary Condition Object', 'View Factor to Ground', 'Frame and Divider Name', 'Multiplier', 'Number of Vertices'] + \
                                    ["X,Y,Z ==> Vertex {} (m)".format(v.index) for v in face.verts]
                                    
                                    if emnode.fclass in ('0', '2'):
                                        paramvs = [('win-', 'door-')[emnode.envi_con_type == 'Door']+'{}_{}'.format(obj.name, face.index), emnode.envi_con_type, mat.name, '{}_{}'.format(obj.name, face.index), obound, 'autocalculate', '', '1', len(face.verts)] + \
                                        ["  {0[0]:.4f}, {0[1]:.4f}, {0[2]:.4f}".format((xav+(vco[0]-xav)*(1 - emnode.farea * 0.01), yav+(vco[1]-yav)*(1 - emnode.farea * 0.01), zav+(vco[2]-zav)*(1 - emnode.farea * 0.01))) for vco in vcos]
                                    else:
                                        paramvs = [('win-', 'door-')[mat.envi_con_type == 'Door']+'{}_{}'.format(obj.name, face.index), emnode.envi_con_type, mat.name, '{}_{}'.format(obj.name, face.index), obound, 'autocalculate', '{}-fad'.format(mat.name), '1', len(face.verts)] + \
                                        ["  {0[0]:.4f}, {0[1]:.4f}, {0[2]:.4f}".format((vco[0] + (1, -1)[vco[0] - xav > 0]*(0.001+emnode.fw, 0)[abs(vco[0] - xav) < 0.0001], vco[1] + (1, -1)[vco[1] - yav > 0]*(0.001+emnode.fw, 0)[abs(vco[1] - yav) < 0.0001], vco[2] + (1, -1)[vco[2] - zav > 0]*(0.001+emnode.fw, 0)[abs(vco[2] - zav) < 0.0001])) for vco in vcos]
                                    
                                    en_idf.write(epentry('FenestrationSurface:Detailed', params, paramvs))
                    
                                elif emnode.envi_con_type == 'Shading' and obj.vi_params.envi_type != '1':
                                    params = ['Name', 'Transmittance Schedule Name', 'Number of Vertices'] + ['X,Y,Z ==> Vertex {} (m)'.format(v.index) for v in face.verts]
                                    paramvs = ['{}_{}'.format(obj.name, face.index), '', len(face.verts)] + ['{0[0]:.4f}, {0[1]:.4f}, {0[2]:.4f}'.format(vco) for vco in vcos]
                                    en_idf.write(epentry('Shading:{}:Detailed'.format('Building'), params, paramvs))
                                    
                                    if emnode.inputs['PV'].links:
                                        pvnode = emnode.inputs['PV'].links[0].from_node
                                        en_idf.write(pvnode.ep_write('{}_{}'.format(obj.name, face.index), face.calc_area()))
                        #                en_idf.write(emnode.pv_ep_write('{}_{}'.format(obj.name, face.index)))
                                        gens.append('{}_{}-pv'.format(obj.name, face.index))
                                        
                            if obj.vi_params.envi_type == '1':
    #                            for face in bm.faces:
                                params = ['Name', 'Transmittance Schedule Name', 'Number of Vertices'] + ['X,Y,Z ==> Vertex {} (m)'.format(v.index) for v in face.verts]
                                paramvs = ['{}_{}'.format(obj.name, face.index), '', len(face.verts)] + ['{0[0]:.4f}, {0[1]:.4f}, {0[2]:.4f}'.format(vco) for vco in vcos]
                                en_idf.write(epentry('Shading:{}:Detailed'.format('Site'), params, paramvs))
                                
                                if emnode.inputs['PV'].links:
                                    pvnode = emnode.inputs['PV'].links[0].from_node
                                    en_idf.write(pvnode.ep_write('{}_{}'.format(obj.name, face.index), face.calc_area()))
                    #                en_idf.write(emnode.pv_ep_write('{}_{}'.format(obj.name, face.index)))
                                    gens.append('{}_{}-pv'.format(obj.name, face.index))
                bm.free()
    
        en_idf.write("\n!-   ===========  ALL OBJECTS IN CLASS: SCHEDULES ===========\n\n")
        params = ('Name', 'Lower Limit Value', 'Upper Limit Value', 'Numeric Type', 'Unit Type')
        paramvs = ("Temperature", -60, 200, "CONTINUOUS", "Temperature")
        en_idf.write(epentry('ScheduleTypeLimits', params, paramvs))
        params = ('Name', 'Lower Limit Value', 'Upper Limit Value', 'Numeric Type')
        paramvs = ("Control Type", 0, 4, "DISCRETE")
        en_idf.write(epentry('ScheduleTypeLimits', params, paramvs))
        params = ('Name', 'Lower Limit Value', 'Upper Limit Value', 'Numeric Type')
        paramvs = ("Fraction", 0, 1, "CONTINUOUS")
        en_idf.write(epentry('ScheduleTypeLimits', params, paramvs))
        params = ['Name']
        paramvs = ["Any Number"]
        en_idf.write(epentry('ScheduleTypeLimits', params, paramvs))
        en_idf.write(epschedwrite('Default outdoor CO2 levels 400 ppm', 'Any number', ['Through: 12/31'], [['For: Alldays']], [[[['Until: 24:00,{}'.format('400')]]]]))
    
        for zn in zonenodes:
            for schedtype in ('VASchedule', 'TSPSchedule', 'HVAC', 'Occupancy', 'Equipment', 'Infiltration'):
                if schedtype == 'HVAC' and zn.inputs[schedtype].links:
                    en_idf.write(zn.inputs[schedtype].links[0].from_node.eptcwrite(zn.zone))
                    try:
                        en_idf.write(zn.inputs[schedtype].links[0].from_node.inputs['Schedule'].links[0].from_node.epwrite(zn.zone+'_hvacsched', 'Fraction'))                            
                    except:
                        en_idf.write(epschedwrite(zn.zone + '_hvacsched', 'Fraction', ['Through: 12/31'], [['For: Alldays']], [[[['Until: 24:00, 1']]]]))
    
                    hsdict = {'HSchedule': '_htspsched', 'CSchedule': '_ctspsched'}
                    tvaldict = {'HSchedule': zn.inputs[schedtype].links[0].from_node.envi_htsp, 'CSchedule': zn.inputs[schedtype].links[0].from_node.envi_ctsp}
                    for sschedtype in hsdict: 
                        if zn.inputs[schedtype].links[0].from_node.inputs[sschedtype].links:
                            en_idf.write(zn.inputs[schedtype].links[0].from_node.inputs[sschedtype].links[0].from_node.epwrite(zn.zone+hsdict[sschedtype], 'Temperature'))                            
                        else:
                            en_idf.write(epschedwrite(zn.zone + hsdict[sschedtype], 'Temperature', ['Through: 12/31'], [['For: Alldays']], [[[['Until: 24:00,{}'.format(tvaldict[sschedtype])]]]]))
    
                elif schedtype == 'Occupancy' and zn.inputs[schedtype].links:
                    osdict = {'OSchedule': '_occsched', 'ASchedule': '_actsched', 'WSchedule': '_wesched', 'VSchedule': '_avsched', 'CSchedule': '_closched'}
                    ovaldict = {'OSchedule': 1, 'ASchedule': zn.inputs[schedtype].links[0].from_node.envi_occwatts, 
                                'WSchedule': zn.inputs[schedtype].links[0].from_node.envi_weff, 'VSchedule': zn.inputs[schedtype].links[0].from_node.envi_airv, 
                                'CSchedule': zn.inputs[schedtype].links[0].from_node.envi_cloth}
                    for sschedtype in osdict:
                        svariant = 'Fraction' if sschedtype == 'OSchedule' else 'Any Number'
                        if zn.inputs[schedtype].links[0].from_node.inputs[sschedtype].links:
                            en_idf.write(zn.inputs[schedtype].links[0].from_node.inputs[sschedtype].links[0].from_node.epwrite(zn.zone + osdict[sschedtype], svariant))
                        else:
                            en_idf.write(epschedwrite(zn.zone + osdict[sschedtype], svariant, ['Through: 12/31'], [['For: Alldays']], [[[['Until: 24:00,{:.3f}'.format(ovaldict[sschedtype])]]]]))
    
                elif schedtype == 'Equipment' and zn.inputs[schedtype].links:
                    if not zn.inputs[schedtype].links[0].from_node.inputs['Schedule'].links:
                        en_idf.write(epschedwrite(zn.zone + '_eqsched', 'Fraction', ['Through: 12/31'], [['For: Alldays']], [[[['Until: 24:00,1']]]]))
                    else:
                        en_idf.write(zn.inputs[schedtype].links[0].from_node.inputs['Schedule'].links[0].from_node.epwrite(zn.zone+'_eqsched', 'Fraction'))
                elif schedtype == 'Infiltration' and zn.inputs[schedtype].links:
                    if not zn.inputs[schedtype].links[0].from_node.inputs['Schedule'].links:
                        en_idf.write(epschedwrite(zn.zone + '_infsched', 'Fraction', ['Through: 12/31'], [['For: Alldays']], [[[['Until: 24:00,{}'.format(1)]]]]))
                    else:
                        en_idf.write(zn.inputs[schedtype].links[0].from_node.inputs['Schedule'].links[0].from_node.epwrite(zn.zone+'_infsched', 'Fraction'))
                elif schedtype == 'VASchedule' and zn.inputs[schedtype].links:
                    en_idf.write(zn.inputs[schedtype].links[0].from_node.epwrite(zn.zone+'_vasched', 'Fraction'))
    
                elif schedtype == 'TSPSchedule' and zn.inputs[schedtype].links:
                    en_idf.write(zn.inputs[schedtype].links[0].from_node.epwrite(zn.zone+'_tspsched', 'Temperature'))
                    
        ssafnodes = [enode for enode in enng.nodes if enode.bl_idname == 'EnViSSFlow']

        for zn in ssafnodes:
            for schedtype in ('VASchedule', 'TSPSchedule'):
                if schedtype == 'VASchedule' and zn.inputs[schedtype].links:
                    en_idf.write(zn.inputs[schedtype].links[0].from_node.epwrite('{}_vasched'.format(zn.name), 'Fraction'))
    
                elif schedtype == 'TSPSchedule' and zn.inputs[schedtype].links:
                    en_idf.write(zn.inputs[schedtype].links[0].from_node.epwrite('{}_tspsched'.format(zn.name), 'Temperature'))
                
        en_idf.write("\n!-   ===========  ALL OBJECTS IN CLASS: GENERATORS ===========\n\n")
        
        if gen:
            params = ('Name', 'Generator List Name', 'Generator Operation Scheme Type', 
                      'Generator Demand Limit Scheme Purchased Electric Demand Limit {W}',
                      'Generator Track Schedule Name Scheme Schedule Name',
                      'Generator Track Meter Scheme Meter Name', 'Electrical Buss Type',
                      'Inverter Name')
            paramvs = 'Electric Load Center', 'PVs', 'Baseload', '0', '', '', 'DirectCurrentWithInverter', 'Simple Ideal Inverter'
            en_idf.write(epentry("ElectricLoadCenter:Distribution", params, paramvs))
            
            params = ('Name', 'Availability Schedule Name', 'Zone Name', 'Radiative Fraction', 'Inverter Efficiency')
            paramvs = ('Simple Ideal Inverter', '', '', 0.0, 1.0)
            en_idf.write(epentry("ElectricLoadCenter:Inverter:Simple", params, paramvs))
            
            
            params = ["Name"] + [item for i, pv in enumerate(gens) for item in ['Generator {} Name'.format(i), 'Generator {} Object Type'.format(i), 
                     'Generator {} Rated Electric Power Output (W)'.format(i), 'Generator {} Availability Schedule Name'.format(i),
                     'Generator {} Rated Thermal to Electrical Power Ratio'.format(i)]] 

            paramvs = ['PVs'] + [item for pv in gens for item in [pv, 'Generator:Photovoltaic', '20000', '', '']] 
            en_idf.write(epentry("ElectricLoadCenter:Generators", params, paramvs))
                        
#            for pv in pvs:
#                en_idf.write(pv.pv_ep_write())


        en_idf.write("\n!-   ===========  ALL OBJECTS IN CLASS: THERMOSTSTATS ===========\n\n")
        for zn in zonenodes:
            for hvaclink in zn.inputs['HVAC'].links:
                en_idf.write(hvaclink.from_node.eptspwrite(zn.zone))
    
        en_idf.write("\n!-   ===========  ALL OBJECTS IN CLASS: EQUIPMENT ===========\n\n")
        for zn in zonenodes:
            for hvaclink in zn.inputs['HVAC'].links:
                hvaczone = hvaclink.from_node
                if not hvaczone.envi_hvact:
                    en_idf.write(zn.inputs['HVAC'].links[0].from_node.epewrite(zn.zone))
        
        en_idf.write("\n!-   ===========  ALL OBJECTS IN CLASS: HVAC ===========\n\n")
        for zn in zonenodes:
            for hvaclink in zn.inputs['HVAC'].links:
                hvacnode = hvaclink.from_node
                if hvacnode.envi_hvact:
                    en_idf.write(hvacnode.hvactwrite(zn.zone))
                else:
                    en_idf.write(hvacnode.ephwrite(zn.zone))
    
        en_idf.write("\n!-   ===========  ALL OBJECTS IN CLASS: OCCUPANCY ===========\n\n")
        for zn in zonenodes:
            for occlink in zn.inputs['Occupancy'].links:
                en_idf.write(occlink.from_node.epwrite(zn.zone))

        en_idf.write("\n!-   ===========  ALL OBJECTS IN CLASS: OTHER EQUIPMENT ===========\n\n")
        for zn in zonenodes:
            for eqlink in zn.inputs['Equipment'].links:
                en_idf.write(eqlink.from_node.oewrite(zn.zone))
       
        en_idf.write("\n!-   ===========  ALL OBJECTS IN CLASS: CONTAMINANTS ===========\n\n")
        zacb = 0
        for zn in zonenodes:
            if not zacb:
                for occlink in zn.inputs['Occupancy'].links:
                    if occlink.from_node.envi_co2 and occlink.from_node.envi_comfort:
                        params = ('Carbon Dioxide Concentration', 'Outdoor Carbon Dioxide Schedule Name', 'Generic Contaminant Concentration', 'Outdoor Generic Contaminant Schedule Name')
                        paramvs = ('Yes', 'Default outdoor CO2 levels 400 ppm', 'No', '')
                        en_idf.write(epentry('ZoneAirContaminantBalance', params, paramvs))
                        zacb = 1
                        break
    
        en_idf.write("\n!-   ===========  ALL OBJECTS IN CLASS: INFILTRATION ===========\n\n")
        for zn in zonenodes:
            for inflink in zn.inputs['Infiltration'].links:
                en_idf.write(inflink.from_node.epwrite(zn.zone))
                
        en_idf.write("\n!-   ===========  ALL OBJECTS IN CLASS: TH ===========\n\n")
        for zn in tcnodes:
            if zn.bl_idname == 'EnViTC':
                en_idf.write(zn.epwrite())
    
        en_idf.write("\n!-   ===========  ALL OBJECTS IN CLASS: AIRFLOW NETWORK ===========\n\n")
        
        if enng and enng['enviparams']['afn']:
            writeafn(exp_op, en_idf, enng)
    
        en_idf.write("!-   ===========  ALL OBJECTS IN CLASS: EMS ===========\n\n")   
        emsprognodes = [pn for pn in enng.nodes if pn.bl_idname == 'No_En_Net_Prog' and not pn.use_custom_color and pn.text_file]
        for prognode in emsprognodes:
            en_idf.write(prognode.epwrite())
        
        en_idf.write("!-   ===========  ALL OBJECTS IN CLASS: REPORT VARIABLE ===========\n\n")
        epentrydict = {"Output:Variable,*,Zone Air Temperature,hourly;\n": node.restt, "Output:Variable,*,Zone Other Equipment Total Heating Rate,hourly;\n": node.resoeg,
                       "Output:Variable,*,Zone Air System Sensible Heating Rate,hourly;\n": node.restwh, "Output:Variable,*,Zone Air System Sensible Cooling Rate,hourly;\n": node.restwc,
                       "Output:Variable,*,Zone Ideal Loads Supply Air Sensible Heating Rate, hourly;\n": node.ressah, "Output:Variable,*,Zone Ideal Loads Heat Recovery Sensible Heating Rate, hourly;\n": node.reshrhw, 
                       "Output:Variable,*,Zone Ideal Loads Supply Air Sensible Cooling Rate,hourly;\n": node.ressac,
                       "Output:Variable,*,Zone Thermal Comfort Fanger Model PMV,hourly;\n": node.rescpm, "Output:Variable,*,Zone Thermal Comfort Fanger Model PPD,hourly;\n": node.rescpp, "Output:Variable,*,AFN Zone Infiltration Volume, hourly;\n": node.resim and enng['enviparams']['afn'],
                       "Output:Variable,*,AFN Zone Infiltration Air Change Rate, hourly;\n": node.resiach and enng['enviparams']['afn'], "Output:Variable,*,Zone Infiltration Current Density Volume,hourly;\n": node.resim and not enng['enviparams']['afn'],
                       "Output:Variable,*,Zone Infiltration Air Change Rate, hourly;\n": node.resiach and not enng['enviparams']['afn'], "Output:Variable,*,Zone Windows Total Transmitted Solar Radiation Rate,hourly;\n": node.reswsg,
                       "Output:Variable,*,AFN Node CO2 Concentration,hourly;\n": node.resco2 and enng['enviparams']['afn'], "Output:Variable,*,Zone Air CO2 Concentration,hourly;\n": node.resco2 and not enng['enviparams']['afn'],
                       "Output:Variable,*,Zone Mean Radiant Temperature,hourly;\n": node.resmrt, "Output:Variable,*,Zone People Occupant Count,hourly;\n": node.resocc,
                       "Output:Variable,*,Zone Air Relative Humidity,hourly;\n": node.resh, "Output:Variable,*,Zone Air Heat Balance Surface Convection Rate, hourly;\n": node.resfhb,
                       "Output:Variable,*,Zone Thermal Chimney Current Density Air Volume Flow Rate,hourly;\n": node.restcvf,
                       "Output:Variable,*,Zone Thermal Chimney Mass Flow Rate,hourly;\n": node.restcmf, "Output:Variable,*,Zone Thermal Chimney Outlet Temperature,hourly;\n": node.restcot,
                       "Output:Variable,*,Zone Thermal Chimney Heat Loss Energy,hourly;\n": node.restchl,"Output:Variable,*,Zone Thermal Chimney Heat Gain Energy,hourly;\n": node.restchg,
                       "Output:Variable,*,Zone Thermal Chimney Volume,hourly;\n": node.restcv, "Output:Variable,*,Zone Thermal Chimney Mass,hourly;\n": node.restcm,
                       "Output:Variable,*,Generator Produced DC Electric Energy,hourly;\n": node.respve, "Output:Variable,*,Generator Produced DC Electric Power,hourly;\n": node.respvw,
                       "Output:Variable,*,Generator PV Array Efficiency,hourly;\n": node.respveff, "Output:Variable,*,Generator PV Cell Temperature,hourly;\n": node.respvt}
        
        for amb in ("Output:Variable,*,Site Outdoor Air Drybulb Temperature,Hourly;\n", "Output:Variable,*,Site Wind Speed,Hourly;\n", "Output:Variable,*,Site Wind Direction,Hourly;\n",
                    "Output:Variable,*,Site Outdoor Air Relative Humidity,hourly;\n", "Output:Variable,*,Site Direct Solar Radiation Rate per Area,hourly;\n", "Output:Variable,*,Site Diffuse Solar Radiation Rate per Area,hourly;\n"):
            en_idf.write(amb)            
        
        for ep in epentrydict:
            if epentrydict[ep]:
                en_idf.write(ep)
        
        if enng['enviparams']['afn']:
            if node.resl12ms:
                for cnode in [cnode for cnode in enng.nodes if cnode.bl_idname == 'No_En_Net_SFlow']:
                    for sno in cnode['sname']:
                        en_idf.write("Output:Variable,{0},AFN Linkage Node 1 to Node 2 Volume Flow Rate,hourly;\nOutput:Variable,{0},AFN Linkage Node 2 to Node 1 Volume Flow Rate,hourly;\n".format(sno))
                        en_idf.write("Output:Variable,{0},AFN Linkage Node 1 to Node 2 Pressure Difference,hourly;\n".format(sno))
                for snode in [snode for snode in enng.nodes if snode.bl_idname == 'No_En_Net_SSFlow']:
                    for sno in snode['sname']:
                        en_idf.write("Output:Variable,{0},AFN Linkage Node 1 to Node 2 Volume Flow Rate,hourly;\nOutput:Variable,{0},AFN Linkage Node 2 to Node 1 Volume Flow Rate,hourly;\n".format(sno))
                        en_idf.write("Output:Variable,{0},AFN Linkage Node 1 to Node 2 Pressure Difference,hourly;\n".format(sno))
            if node.reslof == True:
                for snode in [snode for snode in enng.nodes if snode.bl_idname == 'No_En_Net_SSFlow']:
                    if snode.linkmenu in ('SO', 'DO', 'HO'):
                        for sno in snode['sname']:
                            en_idf.write("Output:Variable,{},AFN Surface Venting Window or Door Opening Factor,hourly;\n".format(sno))
    
        en_idf.write("Output:Table:SummaryReports,\
        AllSummary;              !- Report 1 Name")
        en_idf.close()
        
        if svp['enparams'].get('hvactemplate'):
            os.chdir(svp['viparams']['newdir'])
            ehtempcmd = "ExpandObjects {}".format(os.path.join(svp['viparams']['newdir'], 'in.idf'))
            subprocess.call(ehtempcmd.split())
            shutil.copyfile(os.path.join(svp['viparams']['newdir'], 'expanded.idf'), os.path.join(svp['viparams']['newdir'], 'in.idf')) 
        
        if 'in{}.idf'.format(frame) not in [im.name for im in bpy.data.texts]:
            bpy.data.texts.load(os.path.join(svp['viparams']['newdir'], 'in{}.idf'.format(frame)))
        else:
            bpy.data.texts['in{}.idf'.format(frame)].filepath = os.path.join(svp['viparams']['newdir'], 'in{}.idf'.format(frame))

def pregeo(context, op):
    scene = context.scene
    svp = scene.vi_params
    depsgraph = bpy.context.evaluated_depsgraph_get()

    if context.active_object and context.active_object.mode == 'EDIT':
        bpy.ops.object.editmode_toggle()
    eg = create_coll(context, 'EnVi Geometry')

    for chil in eg.children:
        [bpy.data.objects.remove(o, do_unlink = True, do_id_user=True, do_ui_user=True) for o in chil.objects]
        eg.children.unlink(chil)
        bpy.data.collections.remove(chil)
#    depsgraph = bpy.context.evaluated_depsgraph_get()
    for mesh in bpy.data.meshes:
        if mesh.users == 0:
            bpy.data.meshes.remove(mesh)
            
    for material in bpy.data.materials:
        material.vi_params.envi_export = 0
        material.vi_params['enparams'] = {'area': 0}
        
        if material.users == 0:
            bpy.data.materials.remove(material)
                
    for c in [c for c in bpy.data.collections if c.name != 'EnVi Geometry' and c.name not in [c.name for c in bpy.data.collections['EnVi Geometry'].children]]:
        c.vi_params.envi_zone = 1 if any([o.vi_params.vi_type == '1' for o in c.objects]) else 0
        
        if c.vi_params.envi_zone:
            bpy.data.collections['EnVi Geometry'].children.link(bpy.data.collections.new('EN_{}'.format(c.name.upper())))
            for o in c.objects:
                if o.type == 'MESH' and o.vi_params.envi_type in ('0', '1'):
                    if [f for f in o.data.polygons if o.material_slots and o.material_slots[f.material_index].material and o.material_slots[f.material_index].material.vi_params.envi_nodes]:
                        selobj(context.view_layer, o)
                        bpy.ops.object.duplicate(linked=False)
                        no = context.active_object.copy()  
                        no.location += context.node.geo_offset                       
                        bpy.ops.object.delete()
                        no.name = 'en_{}'.format(c.name) 
                        bpy.data.collections['EN_{}'.format(c.name.upper())].objects.link(no)
            
    for chil in  bpy.data.collections['EnVi Geometry'].children:  
        for o in chil.objects:
            bm = bmesh.new()
            bm.from_mesh(o.evaluated_get(depsgraph).to_mesh())
            o.to_mesh_clear()
            bmesh.ops.remove_doubles(bm, verts = bm.verts, dist = 0.001)                
            bmesh.ops.delete(bm, geom = [e for e in bm.edges if not e.link_faces] + [v for v in bm.verts if not v.link_faces], context = 'VERTS')
            bmesh.ops.delete(bm, geom = [f for f in bm.faces if f.calc_area() < 0.001], context = 'FACES') 
            bmesh.ops.delete(bm, geom = [f for f in bm.faces if o.material_slots[f.material_index].material == None], context = 'FACES') 
            bmesh.ops.delete(bm, geom = [f for f in bm.faces if not o.material_slots[f.material_index].material.vi_params.envi_nodes], context = 'FACES')
            bmesh.ops.delete(bm, geom = [f for f in bm.faces if get_con_node(o.material_slots[f.material_index].material.vi_params).envi_con_type == 'None'], context = 'FACES')
            
            if o.vi_params.envi_type != '1':
                bm.faces.layers.string.new('oname')
                fo = bm.faces.layers.string['oname']
            
                for face in bm.faces:
                    face[fo] = o.name.encode()
                
            bm.to_mesh(o.data)
        
        selobs(context.view_layer, [o.name for o in chil.objects])
        bpy.ops.object.join()
        new_ob = bpy.context.active_object
        new_ob.name = '{}'.format(chil.name)
        bm = bmesh.new()
        bm.from_mesh(new_ob.evaluated_get(depsgraph).to_mesh())
        bm.transform(new_ob.matrix_world)
        new_ob['auto_volume'] = bm.calc_volume()
        bm.free()
        
#        if new_ob.vi_params.envi_type == '1':
#            if 'en_shading' not in [m.name for m in bpy.data.materials]:
#                bpy.data.materials.new('en_shading')
#                bpy.data.materials['en_shading'].use_nodes = 1
#            if not new_ob.material_slots:
#                bpy.ops.object.material_slot_add()                
#            else:
#                while len(new_ob.material_slots) > 1:
#                    bpy.ops.object.material_slot_remove()
#            shadmat = bpy.data.materials['en_shading']
#            new_ob.material_slots[0].material = shadmat
#            new_ob.material_slots[0].material.diffuse_color = (1, 0, 0, 1)
#            
#            if not bpy.data.materials['en_shading'].vi_params.envi_nodes:
#                bpy.ops.material.envi_node('INVOKE_DEFAULT')
#                bpy.data.materials['en_shading'].vi_params.envi_nodes.nodes[0].envi_con_type = 'Shading' 
                        
            

#            shadmat = bpy.data.materials['en_shading']
##            emnode = get_con_node(shadmatmat)
##            shadmat.envi_nodes['envi_con_type']= 'Shading'
#            new_ob.material_slots[0].material = shadmat
#            new_ob.material_slots[0].material.diffuse_color = (1, 0, 0)
        
        if new_ob.vi_params.envi_type == '1':
            if not new_ob.material_slots:
                if 'en_shading' not in [m.name for m in bpy.data.materials]:
                    bpy.data.materials.new('en_shading')
                    bpy.ops.object.material_slot_add()
                    bpy.ops.material.new()
#                    elif len(new_ob.material_slots) > 1:
#                        while len(new_ob.material_slots) > 1:
#                            bpy.ops.object.material_slot_remove()
                    new_ob.material_slots[0].material = bpy.data.materials['en_shading']
                    new_ob.material_slots[0].material.diffuse_color = (1, 0, 0, 1)
            
    if not [ng for ng in bpy.data.node_groups if ng.bl_label == 'EnVi Network']:
        bpy.ops.node.new_node_tree(type='EnViN', name ="EnVi Network") 
        for screen in bpy.data.screens:
            for area in [area for area in screen.areas if area.type == 'NODE_EDITOR' and area.spaces[0].tree_type == 'ViN']:
                area.spaces[0].node_tree = context.node.id_data
#    
    enng = [ng for ng in bpy.data.node_groups if ng.bl_label == 'EnVi Network'][0]
    enng.use_fake_user = True
    enng['enviparams'] = {'wpca': 0, 'wpcn': 0, 'crref': 0, 'afn': 0, 'pcm':0}
#    print([node for node in enng.nodes if hasattr(node, 'zone') and node.zone not in [c.name for c in eg.children]])
    [enng.nodes.remove(node) for node in enng.nodes if hasattr(node, 'zone') and node.zone not in [c.name for c in eg.children]]
#    [enng.nodes.remove(node) for node in enng.nodes if hasattr(node, 'zone') and node.bl_idname == 'EnViZone' and bpy.data.collections[node.zone].vi_params.envi_zone]
#    [enng.nodes.remove(node) for node in enng.nodes if hasattr(node, 'zone') and node.bl_idname == 'No_En_Net_TC' and bpy.data.collections[node.zone].vi_params.envi_zone == 1]            
#
    dcdict = {'Wall':(1, 1, 1, 1), 'Partition':(1, 1, 0, 1), 'Window':(0, 1, 1, 1), 'Roof':(0, 1, 0, 1), 'Ceiling':(1, 1, 0, 1), 'Floor':(0.44,0.185,0.07, 1), 'Shading':(1, 0, 0, 1)}
    pva = 0
    linklist = []
    for coll in eg.children:
        cvp = coll.vi_params
        cvp['enparams'] = {}
        cvp['enparams']['floorarea'] = 0.0
        
        for oi, obj in enumerate(coll.objects):
            ovp = obj.vi_params
            
            for k in ovp.keys():
                if k not in ('envi_type', 'vi_type', 'envi_oca', 'envi_ica'):
                    del ovp[k]
            ovp['enparams'] = {}
            omats = [ms.material for ms in obj.material_slots]
    
            if ovp.envi_type in ('0', '2') and not omats:
                op.report({'WARNING'}, 'Object {} is specified as a thermal zone but has no materials'.format(obj.name))
            elif None in omats:
                op.report({'WARNING'}, 'Object {} has an empty material slot'.format(obj.name))
            else:
                try:
                    ofa = 0
                    for face in obj.data.polygons:
                        for node in omats[face.material_index].vi_params.envi_nodes.nodes:
                            if node.bl_idname == 'No_En_Mat_Con' and node.active:
                                if node.envi_con_type =='Floor':
                                    ofa += facearea(obj, face)
                        omats[face.material_index].vi_params['enparams']['area'] += facearea(obj, face)
#                            elif node.inputs['PV'].links:
#                                if node.inputs['PV'].links[0].from_node.inputs['PV Generator'].links:
#                                    
                                
                    ovp['enparams']["floorarea"] = ofa
                    cvp['enparams']["floorarea"] += ofa
                except Exception as e:
                    logentry('There is a problem with the area calculation for object {}: {}'.format(obj.name, e))
                    ovp['enparams']["floorarea"] = 1
                    
                if ovp.envi_type in ('0', '2'):                        
                    ezdict = {'0': 'No_En_Net_Zone', '2': 'No_En_Net_TC'}            
                
        
        
#                for mats in omats:
#                    if 'en_'+mats.name not in [mat.name for mat in bpy.data.materials]:
#                        mats.copy().name = 'en_'+mats.name
##                if '8' in (mats.envi_type_lo, mats.envi_type_l1, mats.envi_type_l2, mats.envi_type_l3, mats.envi_type_l4):
##                    enng['enviparams']['pcm'] = 1
#    
#            selobj(context.view_layer, obj)
#            en_obj = context.view_layer.objects.active
##            bpy.ops.object.duplicate()    
##            bpy.ops.object.transform_apply(location=False, rotation=True, scale=True)
##            en_obj = context.view_layer.objects.active
##            bpy.ops.object.origin_set(type='ORIGIN_GEOMETRY')
#            selmesh('desel')
#            enomats = [enom for enom in en_obj.data.materials if enom and enom.vi_params.envi_nodes and enom.vi_params.envi_nodes.nodes and not any([n.use_custom_color for n in enom.vi_params.envi_nodes.nodes])]
##            print('e', [enom for enom in en_obj.data.materials if enom and enom.envi_nodes and enom.envi_nodes.nodes])
#            en_obj.name, en_obj.data.name = 'en_'+obj.name, en_obj.data.name
#            mis = [f.material_index for f in no.data.polygons]
            done_mats = []
            for s, sm in enumerate(obj.material_slots): 
                if sm.material and sm.material not in done_mats:
                    done_mats.append(sm.material)
                    mat = sm.material
                    mvp = mat.vi_params
                    
                    if mvp.envi_nodes and mvp.envi_nodes.name != mat.name:
                        mvp.envi_nodes = mvp.envi_nodes.copy()
                        mvp.envi_nodes.name = mat.name
                    
                    if mvp.envi_nodes and mvp.envi_nodes.users > 1:
                       mvp.envi_nodes.user_clear() 
                       mvp.envi_nodes = bpy.data.node_groups[mat.name]
                       
                    emnode = get_con_node(mvp)
                    
                    if not emnode:
                        op.report({'WARNING'}, 'The {} material has no node tree. This material has not been exported.'.format(mat.name))
                    elif any([n.use_custom_color for n in emnode.ret_nodes()]):
                        op.report({'WARNING'}, 'There is a red node in the {} material node tree. This material has not been exported.'.format(mat.name))
                    else:                       
                        emnode.ret_uv()
                        mct = 'Partition' if emnode.envi_con_con == 'Zone' else emnode.envi_con_type
    ##                    mat.envi_nodes['enmatparams']['boundary'] = emnode.envi_boundary
    ##                    mat.envi_nodes['enmatparams']['airflow'] = emnode.af_surface
    ##                    mat.envi_nodes['enmatparams']['tm'] = emnode.envi_thermalmass
    #                    
    #                    if s in mis:
                        mvp.envi_export = True    

                        if emnode.envi_con_type in dcdict:
                            mat.diffuse_color = dcdict[mct] 
                        if emnode.inputs['PV'].links:
                            mat.diffuse_color = (1, 1, 0, 1)
                 
##                if mct not in ('None', 'Shading', 'Aperture', 'Window'):
##                    print('mct', mct)
##                retuval(mat)
#                
#            for poly in no.data.polygons:
#                poly.select = False
#                mat = no.data.materials[poly.material_index]
#                emnode = get_con_node(mat.vi_params)
#                
#                try:
#                    if poly.area < 0.001:
#                        poly.select = True 
#                    if not mat.vi_params.envi_nodes:                        
#                        poly.select = True
#                    elif emnode.envi_con_type == 'None':
#                        poly.select = True
#                    elif any([n.use_custom_color for n in mat.vi_params.envi_nodes.nodes]):
#                        poly.select = True 
#                except Exception as e:
#                    print(e)
##                    logentry("Problem with EnVi export of {}".format(en_obj.name))
#                    poly.select = True

#            selmesh('delf')  
#            for edge in no.data.edges:
#                if edge.is_loose:
#                    edge.select = True
#                    for vi in edge.vertices:
#                        no.data.vertices[vi].select = True          
#            selmesh('delv')
##            selmesh('dele')        
#            if not oi:
#                bm_total = bmesh.new()
#                bm_total.faces.layers.string.new('oname')
#                fo = bm_total.faces.layers.string['oname']
#                om = obj.matrix_world
#            tmesh = en_obj.evaluated_get(depsgraph).to_mesh()
##            tmesh = en_obj.data.copy()
#            tmesh.transform(en_obj.matrix_world)
#            bm_total.from_mesh(tmesh)
#            en_obj.to_mesh_clear()
#            
#            for face in bm_total.faces[-len(tmesh.polygons):]:
#                face[fo] = obj.name.encode()
            
            
                
#            bm = bmesh.new()
#            bm.from_mesh(en_obj.data)
#            bm.transform(en_obj.matrix_world)
        
#        en_obj['auto_volume'] = bm_total.calc_volume()
#            print(bm.calc_volume())
        
#
#            obj['children'] = en_obj.name

        for link in enng.links:
            if link.from_socket.bl_idname in ('So_En_Net_Bound', 'So_En_Net_SFlow', 'So_En_Net_SSFlow'):
                linklist.append([link.from_socket.node.name, link.from_socket.name, link.to_socket.node.name, link.to_socket.name])
                enng.links.remove(link)
        if context.view_layer.objects[coll.name].vi_params.envi_type != '1':
            if coll.name not in [node.zone for node in enng.nodes if hasattr(node, 'zone')]:
                enng.nodes.new(type = ezdict[context.view_layer.objects[coll.name].vi_params.envi_type]).zone = coll.name
            else:
                for node in enng.nodes:
                    if hasattr(node, 'zone') and node.zone == coll.name:
                        node.zupdate(bpy.context)
                    if hasattr(node, 'emszone') and node.emszone == coll.name:
                        node.zupdate(bpy.context) 
                    if hasattr(node, 'zone') and node.zone == coll.name:
                        node.uvsockupdate()    
#        for node in enng.nodes:
#            
#                
        for node in enng.nodes:
#            if hasattr(node, 'zone') and node.zone == coll.name:
#                node.uvsockupdate()
    #    eg['enparams']['pva'] = sum([mat.vi_params['enparams']['area'] for mat in bpy.data.materials if get_con_node(mat).links])
    #                
    ##            print(enomats, [get_con_node(mat).name for mat in enomats])    
    #                if node.envi_afsurface:
            if [sock.bl_idname in ('So_En_Net_SFlow', 'So_En_Net_SFlow') for sock in node.inputs]:
                enng['enviparams']['afn'] = 1
#                
            if 'No_En_Net_ACon' not in [node.bl_idname for node in enng.nodes]:
                enng.nodes.new(type = 'No_En_Net_ACon')         
                enng.use_fake_user = 1
#                
##            scene.layers[0], scene.layers[1] = True, False
#            selobj(context.view_layer, obj)
##            scene.objects.active = obj
#        
        
#            en_obj.layers[1], en_obj.layers[0] = True, False
          
    for ll in linklist:
        try:
            enng.links.new(enng.nodes[ll[0]].outputs[ll[1]], enng.nodes[ll[2]].inputs[ll[3]])
        except Exception as e:
            print('Link', e)
            
def writeafn(exp_op, en_idf, enng):
    if [enode for enode in enng.nodes if enode.bl_idname == 'No_En_Net_ACon'] and not [enode for enode in enng.nodes if enode.bl_idname == 'No_En_Net_Zone']:
        [enng.nodes.remove(enode) for enode in enng.nodes if enode.bl_idname == 'No_En_Net_ACon']

    for connode in [enode for enode in enng.nodes if enode.bl_idname == 'No_En_Net_ACon']:
         en_idf.write(connode.epwrite(exp_op, enng))        

    for crnode in [enode for enode in enng.nodes if enode.bl_idname == 'EnViCrRef']:
        en_idf.write(crnode.epwrite())
        enng['enviparams']['crref'] = 1

    extnodes = [enode for enode in enng.nodes if enode.bl_idname == 'No_En_Net_Ext']
    zonenodes = [enode for enode in enng.nodes if enode.bl_idname == 'No_En_Net_Zone']
    ssafnodes = [enode for enode in enng.nodes if enode.bl_idname == 'No_En_Net_SSFlow']
    safnodes = [enode for enode in enng.nodes if enode.bl_idname == 'No_En_Net_SFlow']

    if enng['enviparams']['wpca'] == 1:
        for extnode in extnodes:
            en_idf.write(extnode.epwrite(enng))
    for enode in zonenodes:
        en_idf.write(enode.epwrite())
    for enode in ssafnodes + safnodes:
        en_idf.write(enode.epwrite(exp_op, enng))

