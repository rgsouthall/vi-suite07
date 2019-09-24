import bpy
from collections import OrderedDict
from .vi_func import newrow, retdates, logentry
#from .envi_mat import envi_materials, envi_constructions

#envi_mats = envi_materials()
#envi_cons = envi_constructions()

class VI_PT_3D(bpy.types.Panel):
    '''VI-Suite 3D view panel'''
    bl_label = "VI Display"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "VI-Suite"
    
    def draw(self, context):
        scene = context.scene
        cao = context.active_object
        layout = self.layout

        try:
            if cao and cao.active_material and cao.active_material.get('bsdf') and cao.active_material['bsdf']['type'] == ' ' and cao.vi_type == '5' and scene['viparams'].get('vidisp'):                
                if scene['viparams']['vidisp'] != 'bsdf_panel':
                    row = layout.row()
                    row.operator("view3d.bsdf_display", text="BSDF Display") 
                else:
                    newrow(layout, 'BSDF max:', scene, "vi_bsdfleg_max")
                    newrow(layout, 'BSDF min:', scene, "vi_bsdfleg_min")
                    newrow(layout, 'BSDF scale:', scene, "vi_bsdfleg_scale")
                    newrow(layout, 'BSDF colour:', scene, "vi_leg_col")
        
        except Exception as e:
            logentry("Problem with BSDF panel display: {}".format(e))

        if scene.get('viparams') and scene['viparams'].get('vidisp'): 
            if scene['viparams']['vidisp'] == 'wr' and 'Wind_Plane' in [o['VIType'] for o in bpy.data.objects if o.get('VIType')]:
                row = layout.row()
                row.operator('view3d.wrdisplay', text = 'Wind Metrics')#('INVOKE_DEFAULT'')
                
            elif scene['viparams']['vidisp'] == 'wrpanel' and scene.vi_display:
                newrow(layout, 'Wind metric:', scene, 'wind_type')
                newrow(layout, 'Colour:', scene, 'vi_leg_col')
                
            elif scene['viparams']['vidisp'] == 'sp' and scene.vi_display:
                newrow(layout, "Latitude:", scene.vi_params, 'latitude')
                newrow(layout, "Longitude:", scene.vi_params, 'longitude')

                (sdate, edate) = retdates(scene.solday, 365, 2015)
                    
                time_disps = ((("Day of year: {}/{}".format(sdate.day, sdate.month), "sp_sd"), ("Time of day:", "sp_sh")), [("Time of day:", "sp_sh")], [("Day of year: {}/{}".format(sdate.day, sdate.month), "sp_sd")])
                
                for i in time_disps[int(scene['spparams']['suns'])]:
                    newrow(layout, i[0], scene.vi_params, i[1])
                
                for i in (("Sun strength:", "sp_sun_strength"), ("Sun angle:", "sp_sun_angle")):
                        newrow(layout, i[0], scene.vi_params, i[1])
                
                newrow(layout, "Line width:", scene.vi_params, 'sp_line_width')
                newrow(layout, "Solstice colour:", scene.vi_params, 'sp_season_main') 
                newrow(layout, "Hour main colour:", scene.vi_params, 'sp_hour_main')            
                newrow(layout, "Hour dash colour:", scene.vi_params, 'sp_hour_dash')  
                newrow(layout, "Hour dash ratio:", scene.vi_params, 'sp_hour_dash_ratio')
                newrow(layout, "Hour dash density:", scene.vi_params, 'sp_hour_dash_density')
                newrow(layout, "Sun size:", scene.vi_params, 'sp_sun_size')
                newrow(layout, "Sun colour:", scene.vi_params, 'sp_sun_colour')
                newrow(layout, "Globe colour:", scene.vi_params, 'sp_globe_colour')
                
                time_disps = ((("Display time:", "sp_td"), ("Display hours:", "sp_hd")), [("Display hours:", "sp_hd")], [("Display hours:", "sp_hd")])
                
                for i in time_disps[int(scene['spparams']['suns'])]:
                    newrow(layout, i[0], scene.vi_params, i[1])
                
                if (scene['spparams']['suns'] == '0' and (scene.vi_params.sp_td or scene.vi_params.sp_hd)) or scene.vi_params.sp_hd:
                    for i in (("Font size:", "display_rp_fs"), ("Font colour:", "display_rp_fc"), ("Font shadow:", "display_rp_sh")):
                        newrow(layout, i[0], scene.vi_params, i[1])
                    if scene.vi_params.display_rp_sh:
                        newrow(layout, "Shadow colour:", scene.vi_params, "display_rp_fsh")
                
            elif scene['viparams']['vidisp'] in ('svf', 'ss', 'li', 'lc'):
                row = layout.row()
                row.prop(scene, "vi_disp_3d")                 
                row = layout.row()
                
                if scene['viparams']['vidisp'] == 'svf':
                    row.operator("view3d.svfdisplay", text="Sky View Display")
                elif scene['viparams']['vidisp'] == 'ss':
                    row.operator("view3d.ssdisplay", text="Shadow Display")
                else:
                    row.operator("view3d.livibasicdisplay", text="Radiance Display")

            elif scene['viparams']['vidisp'] in ('sspanel', 'lipanel', 'lcpanel', 'svfpanel') and [o for o in bpy.data.objects if o.lires] and scene.vi_display:
                row = layout.row()
                row.prop(context.space_data, "show_only_render")

                if not scene.ss_disp_panel:
                    if scene['viparams']['visimcontext'] == 'LiVi CBDM':
                        if scene['liparams']['unit'] in ('DA (%)', 'sDA (%)', 'UDI-f (%)', 'UDI-s (%)', 'UDI-a (%)', 'UDI-e (%)', 'ASE (hrs)', 'Min lux', 'Max lux', 'Avg lux'):
                            newrow(layout, 'Result type:', scene, "li_disp_da")
                        elif scene['liparams']['unit'] in ('Mlxh', u'kWh/m\u00b2 (f)', u'kWh/m\u00b2 (v)', 'kWh (f)', 'kWh (v)'):
                            newrow(layout, 'Result type:', scene, "li_disp_exp")
                        elif scene['liparams']['unit'] in ('kWh', 'kWh/m2'):
                            newrow(layout, 'Result type:', scene, "li_disp_irrad")

                    elif scene['viparams']['visimcontext'] == 'LiVi Compliance': 
                        if scene['liparams']['unit'] in ('sDA (%)', 'ASE (hrs)'):
                            newrow(layout, 'Metric:', scene, 'li_disp_sda')
                        else:
                            newrow(layout, 'Metric:', scene, 'li_disp_sv')
                            
                    elif scene['viparams']['visimcontext'] == 'LiVi Basic':
                        newrow(layout, 'Metric:', scene, 'li_disp_basic')
                    
                    newrow(layout, 'Legend unit:', scene, "vi_leg_unit")
                    newrow(layout, 'Processing:', scene, "vi_res_process") 

                    if scene.vi_res_process == '1':                    
                        newrow(layout, 'Modifier:', scene, "vi_res_mod")
                    elif scene.vi_res_process == '2':
                        layout.prop_search(scene, 'script_file', bpy.data, 'texts', text='File', icon='TEXT')
                       
                    newrow(layout, 'Legend max:', scene, "vi_leg_max")
                    newrow(layout, 'Legend min:', scene, "vi_leg_min")
                    newrow(layout, 'Legend scale:', scene, "vi_leg_scale")
                    newrow(layout, 'Legend colour:', scene, "vi_leg_col")
                    newrow(layout, 'Legend levels:', scene, "vi_leg_levels")
                    newrow(layout, 'Emitter materials:', scene, "vi_disp_mat")
                    
                    if scene.vi_disp_mat:
                        newrow(layout, 'Emitter strength:', scene, "vi_disp_ems")
                    
                    if scene['liparams']['unit'] in ('DA (%)', 'sDA (%)', 'UDI-f (%)', 'UDI-s (%)', 'UDI-a (%)', 'UDI-e (%)', 'ASE (hrs)', 'Max lux', 'Avg lux', 'Min lux', 'kWh', 'kWh/m2'):
                        newrow(layout, 'Scatter max:', scene, "vi_scatter_max")
                        newrow(layout, 'Scatter min:', scene, "vi_scatter_min")
                        
                if cao and cao.type == 'MESH':
                    newrow(layout, 'Draw wire:', scene, 'vi_disp_wire')                    
                
                if int(context.scene.vi_disp_3d) == 1:
                    newrow(layout, "3D Level", scene, "vi_disp_3dlevel")                        
                
                newrow(layout, "Transparency", scene, "vi_disp_trans")

                if context.mode != "EDIT":
                    row = layout.row()
                    row.label(text="{:-<48}".format("Point visualisation "))
                    propdict = OrderedDict([('Enable', "vi_display_rp"), ("Selected only:", "vi_display_sel_only"), ("Visible only:", "vi_display_vis_only"), ("Font size:", "vi_display_rp_fs"), ("Font colour:", "vi_display_rp_fc"), ("Font shadow:", "vi_display_rp_sh"), ("Shadow colour:", "vi_display_rp_fsh"), ("Position offset:", "vi_display_rp_off")])
                    for prop in propdict.items():
                        newrow(layout, prop[0], scene, prop[1])
                    row = layout.row()
                    row.label(text="{:-<60}".format(""))
 
            elif scene['viparams']['vidisp'] in ('en', 'enpanel'):
                fs, fe = scene['enparams']['fs'], scene['enparams']['fe']
                sedt = scene.en_disp_type
                resnode = bpy.data.node_groups[scene['viparams']['resnode'].split('@')[1]].nodes[scene['viparams']['resnode'].split('@')[0]]

                if sedt == '1':
                    zresdict = {}
                    lmetrics = []
                    vresdict = {"Max Flow in": "resazlmaxf_disp", "Min Flow in": "resazlminf_disp", "Avg Flow in": "resazlavef_disp"} 
                else: 
                    lmetrics, zmetrics = scene['enparams']['lmetrics'], scene['enparams']['zmetrics']
                    zresdict = {"Temperature (degC)": "reszt_disp", 'Humidity (%)': 'reszh_disp', 'Heating (W)': 'reszhw_disp', 'Cooling (W)': 'reszcw_disp', 
                                'CO2 (ppm)': 'reszco_disp', 'PMV': 'reszpmv_disp', 'PPD (%)': 'reszppd_disp', 'Solar gain (W)': 'reszsg_disp', 
                                'Air heating (W)': 'reszahw_disp', 'Air cooling (W)': 'reszacw_disp', 'HR heating (W)': 'reshrhw_disp'}
                    vresdict = {"Opening Factor": "reszof_disp", "Linkage Flow in": "reszlf_disp"}  

                if scene['viparams']['vidisp'] == 'en': 
                    newrow(layout, 'Static/Parametric', scene, 'en_disp_type')
                    if sedt == '1':
                        row = layout.row()               
                        row.prop(resnode, '["AStart"]')
                        row.prop(resnode, '["AEnd"]')
                    else:  
                        if fe > fs:                        
                            newrow(layout, 'Frame:', resnode, '["AStart"]')

                        row = layout.row() 
                        row.label(text = 'Start/End day:')
                        row.prop(resnode, '["Start"]')
                        row.prop(resnode, '["End"]')
                        row = layout.row() 
                        row.label(text = 'Ambient')
                        row = layout.row() 
                        row.prop(scene, 'resaa_disp')
                        row.prop(scene, 'resas_disp')
                        
                        for ri, rzname in enumerate(zmetrics):
                            if ri == 0:                    
                                row = layout.row()
                                row.label(text = 'Zone')                    
                            if not ri%2:
                                row = layout.row()  
                            if rzname in zresdict:
                                row.prop(scene, zresdict[rzname])
                        
                        for ri, rname in enumerate(lmetrics):
                            if ri == 0:                    
                                row = layout.row()
                                row.label(text = 'Ventilation')                    
                            if not ri%2:
                                row = layout.row()                            
                            if rname in vresdict:
                                row.prop(scene, vresdict[rname])  
                        if lmetrics:    
                            newrow(layout, 'Link to object', scene, 'envi_flink')  
                        
                        row = layout.row() 
   
                    if sedt == '0':
                        row.operator("view3d.endisplay", text="EnVi Display")
                    elif sedt == '1':
                        row.operator("view3d.enpdisplay", text="EnVi Display")
                        
            if scene['viparams']['vidisp'] == 'enpanel':                                
                if sedt == '0':
                    newrow(layout, 'Display unit:', scene, 'en_disp_unit')  
                    newrow(layout, 'Bar colour:', scene, "vi_leg_col")

                    if fe > fs:
                        newrow(layout, 'Parametric frame:', resnode, '["AStart"]')

                    envimenudict = {'Temperature (degC)': ('en_temp_min', 'en_temp_max'), 'Humidity (%)' : ('en_hum_min', 'en_hum_max'), 'Heating (W)': ('en_heat_min', 'en_heat_max'),
                                'Cooling (W)': ('en_cool_min', 'en_cool_max'), 'Solar gain (W)': ('en_shg_min', 'en_shg_max'), 'CO2 (ppm)': ('en_co2_min', 'en_co2_max'),
                                'PMV': ('en_pmv_min', 'en_pmv_max'), 'PPD (%)': ('en_ppd_min', 'en_ppd_max'), 'Air heating (W)': ('en_aheat_min', 'en_aheat_max'), 
                                'Air cooling (W)': ('en_acool_min', 'en_acool_max'), 'HR heating (W)': ('en_hrheat_min', 'en_hrheat_max'), 'Heat balance (W)': ('en_heatb_min', 'en_heatb_max'),
                                'Occupancy': ('en_occ_min', 'en_occ_max'), 'Infiltration (ACH)': ('en_iach_min', 'en_iach_max'), 'Infiltration (m3/s)': ('en_im3s_min', 'en_im3s_max'),
                                'Equipment (W)': ('en_eq_min', 'en_eq_max')}
               
                    for envirt in envimenudict:
                        if envirt in zmetrics:
                            row = layout.row()
                            row.label(envirt)
                            row.prop(scene, envimenudict[envirt][0])
                            row.prop(scene, envimenudict[envirt][1])
                
                elif sedt == '1':
                    newrow(layout, 'Display unit:', scene, 'en_disp_punit')  
                    newrow(layout, 'Legend colour:', scene, "vi_leg_col")
                    row = layout.row()
                    row.label('Bar chart range:')
                    row.prop(scene, 'bar_min')
                    row.prop(scene, 'bar_max')
                                            
            if scene.vi_display:            
                newrow(layout, 'Display active', scene, 'vi_display')
        
            
class VIMatPanel(bpy.types.Panel):
    bl_label = "VI-Suite Material"
    bl_space_type = "PROPERTIES"
    bl_region_type = "WINDOW"
    bl_context = "material"

    @classmethod
    def poll(cls, context):
        return context.material

    def draw(self, context):
        cm, scene = context.material, context.scene
        layout = self.layout
        newrow(layout, 'Material type', cm, "mattype")
        if cm.mattype == '0':
            rmmenu(layout, cm)
            if not cm.envi_nodes or (cm.envi_nodes.name != cm.name and cm.envi_nodes.name in [m.name for m in bpy.data.materials]):# in bpy.data.node_groups:
                row = layout.row()
                row.operator("material.envi_node", text = "Create EnVi Nodes")
        
        elif cm.mattype == '1':  
            if scene.get('viparams') and scene['viparams'].get('viexpcontext') and scene['viparams']['viexpcontext'] == 'LiVi Compliance':
                connode = bpy.data.node_groups[scene['viparams']['connode'].split('@')[1]].nodes[scene['viparams']['connode'].split('@')[0]]
                coptions = connode['Options']

                if coptions['canalysis'] == '0':
                    if coptions['bambuild'] == '2':
                        newrow(layout, "Space type:", cm, 'hspacemenu')
                    elif coptions['bambuild'] == '3':
                        newrow(layout, "Space type:", cm, 'brspacemenu')
                        if cm.brspacemenu == '2':
                            row = layout.row()
                            row.prop(cm, 'gl_roof')
                    elif coptions['bambuild'] == '4':
                        newrow(layout, "Space type:", cm, 'respacemenu')
                elif coptions['canalysis'] == '1':
                    newrow(layout, "Space type:", cm, 'crspacemenu')
                elif coptions['canalysis'] == '2':
                    if coptions['bambuild'] == '2':
                        newrow(layout, "Space type:", cm, 'hspacemenu')
                    if coptions['bambuild'] == '3':
                        newrow(layout, "Space type:", cm, 'brspacemenu')
#                    elif coptions['canalysis'] == '3':
#                        newrow(layout, "Space type:", cm, 'lespacemenu')                   
            rmmenu(layout, cm)
        
        elif cm.mattype == '2':
            fvsimnode = bpy.data.node_groups[scene['viparams']['fvsimnode'].split('@')[1]].nodes[scene['viparams']['fvsimnode'].split('@')[0]] if scene.get('viparams') and 'fvsimnode' in scene['viparams'] else 0
            newrow(layout, "Type:", cm, "flovi_bmb_type")
            if fvsimnode:
                context.scene['flparams']['solver'] = fvsimnode.solver
#            newrow(layout, "Type:", cm, "flovi_bmb_subtype")
                if cm.flovi_bmb_type in ('0', '1'):
                    newrow(layout, "p type:", cm, "flovi_bmbp_subtype")
                    
                    if cm.flovi_bmbp_subtype in ('fixedValue', 'totalPressure'):
                        if cm.flovi_bmbp_subtype == 'totalPressure':
                            newrow(layout, "p0 value:", cm, "flovi_bmbp_p0val")
                            newrow(layout, "Gamma value:", cm, "flovi_bmbp_gamma")
                        newrow(layout, "Field value:", cm, "flovi_p_field")
                        
                        if not cm.flovi_p_field:
                            newrow(layout, "Pressure value:", cm, "flovi_bmbp_val")                        
                            
                    newrow(layout, "U type:", cm, "flovi_bmbu_subtype")
                    if cm.flovi_bmbu_subtype in ('fixedValue', 'pressureInletOutletVelocity'):
                        newrow(layout, "Field value:", cm, "flovi_u_field")
                        if not cm.flovi_u_field:
                            newrow(layout, "Velocity value:", cm, "flovi_bmbu_val")
                            
                    if fvsimnode.solver in ('simpleFoam', 'buoyantSimpleFoam', 'buoyantBoussinesqSimpleFoam') and fvsimnode.turbulence != '0':    
                        newrow(layout, "Nut type:", cm, "flovi_bmbnut_subtype")
                        if cm.flovi_bmbnut_subtype == 'fixedValue':
                            newrow(layout, "Nut field:", cm, "flovi_nut_field")
                            if not cm.flovi_u_field:
                                newrow(layout, "Nut value:", cm, "flovi_bmbnut_val")
                        if fvsimnode.turbulence == 'kEpsilon':
                            newrow(layout, "k type:", cm, "flovi_bmbk_subtype")
                            if cm.flovi_bmbk_subtype == 'fixedValue':
                                newrow(layout, "K field:", cm, "flovi_k_field")
                                if not cm.flovi_k_field:
                                    newrow(layout, "K value:", cm, "flovi_bmbk_val")
                            newrow(layout, "Epsilon type:", cm, "flovi_bmbe_subtype")
                            if cm.flovi_bmbe_subtype == 'fixedValue':
                                newrow(layout, "Epsilon field:", cm, "flovi_e_field")
                                if not cm.flovi_e_field:
                                    newrow(layout, "Epsilon value:", cm, "flovi_bmbe_val")
                        elif fvsimnode.turbulence == 'kOmega':
                            newrow(layout, "k type:", cm, "flovi_bmbk_subtype")
                            if cm.flovi_bmbk_subtype == 'fixedValue':
                                newrow(layout, "k field:", cm, "flovi_k_field")
                                if not cm.flovi_k_field:
                                    newrow(layout, "k value:", cm, "flovi_bmbk_val")
                            newrow(layout, "Omega type:", cm, "flovi_bmbo_subtype")
                            if cm.flovi_bmbo_subtype == 'fixedValue':
                                newrow(layout, "Omega field:", cm, "flovi_o_field")
                                if not cm.flovi_o_field:
                                    newrow(layout, "Omega value:", cm, "flovi_bmbo_val")
                        
                        elif fvsimnode.turbulence == 'SpalartAllmaras':
                            newrow(layout, "Nutilda type:", cm, "flovi_bmbnutilda_subtype")
                            if cm.flovi_bmbnutilda_subtype == 'fixedValue':
                                newrow(layout, "Nutilda field:", cm, "flovi_nutilda_field")
                                if not cm.flovi_nutilda_field:
                                    newrow(layout, "Nutilda value:", cm, "flovi_bmbnutilda_val")
                        
                    if fvsimnode.solver in ('buoyantSimpleFoam', 'buoyantBoussinesqSimpleFoam'):  
                        newrow(layout, "T type:", cm, "flovi_bmbt_subtype")
                        if cm.flovi_bmbt_subtype == 'fixedValue':
                            newrow(layout, "T field:", cm, "flovi_t_field")
                            if not cm.flovi_t_field:
                                newrow(layout, "T value:", cm, "flovi_bmbt_val")
                            
#                    newrow(layout, "Pressure type:", cm, "flovi_bmwp_type")
#                    if cm.flovi_bmwp_type == 'fixedValue':
#                        newrow(layout, "Pressure value:", cm, "flovi_b_sval")
#                        
#                    newrow(layout, "Velocity type:", cm, "flovi_bmwu_type")
#                    newrow(layout, "Field value:", cm, "flovi_u_field")
#                    if not cm.flovi_u_field:
#                        newrow(layout, 'Velocity:', cm, 'flovi_b_vval')
##                split = layout.split()
##                col = split.column(align=True)
##                col.label(text="Velocity:")
##                col.prop(cm, "flovi_bmu_x")
##                col.prop(cm, "flovi_bmu_y")
##                col.prop(cm, "flovi_bmu_z")
#                
#                if fvsimnode and fvsimnode.solver != 'icoFoam':
#                    if fvsimnode.buoyancy or fvsimnode.radiation:
#                       newrow(layout, "Temperature type:", cm, "flovi_bmot_type")
#                       if cm.flovi_bmot_type == 'fixedValue':
#                           newrow(layout, "Temperature value:", cm, "flovi_temp")
#                    newrow(layout, "nut type:", cm, "flovi_bmwnut_type")
#                    if fvsimnode.turbulence == 'SpalartAllmaras':                        
#                        newrow(layout, "nuTilda type:", cm, "flovi_bmwnutilda_type")
#                    elif fvsimnode.turbulence == 'kEpsilon':
#                        newrow(layout, "k type:", cm, "flovi_bmwk_type")
#                        newrow(layout, "Epsilon type:", cm, "flovi_bmwe_type")
#                    elif fvsimnode.turbulence == 'komega':
#                        newrow(layout, "k type:", cm, "flovi_bmwk_type")
#                        newrow(layout, "Omega type:", cm, "flovi_bmwe_type")
#
##                newrow(layout, "nuTilda:", cm, "flovi_bmnutilda")
##                split = layout.split()
##                col = split.column(align=True)
##                col.label(text="nuTilda:")
##                col.prop(cm, "flovi_bmnut")
##                col.prop(cm, "flovi_bmwnut_y")
##                col.prop(cm, "flovi_bmwnut_z")
#            elif cm.flovi_bmb_type == '1':
#                newrow(layout, "Pressure sub-type:", cm, "flovi_bmip_type")
#                if cm.flovi_bmip_type == 'fixedValue':
#                    newrow(layout, "Pressure value:", cm, "flovi_b_sval")
#                newrow(layout, "Velocity sub-type:", cm, "flovi_bmiu_type")
#                newrow(layout, "Field value:", cm, "flovi_u_field")
#                if not cm.flovi_u_field:
#                    newrow(layout, 'Velocity:', cm, 'flovi_b_vval')
#                if fvsimnode and fvsimnode.solver != 'icoFoam':
#                    newrow(layout, "nut type:", cm, "flovi_bminut_type")
#                    if fvsimnode.turbulence == 'SpalartAllmaras':                        
#                        newrow(layout, "nuTilda type:", cm, "flovi_bminutilda_type")
#                    elif fvsimnode.turbulence == 'kEpsilon':
#                        newrow(layout, "k type:", cm, "flovi_bmik_type")
#                        newrow(layout, "Epsilon type:", cm, "flovi_bmie_type")
#                    elif fvsimnode.turbulence == 'kOmega':
#                        newrow(layout, "k type:", cm, "flovi_bmik_type")
#                        newrow(layout, "Omega type:", cm, "flovi_bmio_type")
#            
#            elif cm.flovi_bmb_type == '2':
#                newrow(layout, "Pressure sub-type:", cm, "flovi_bmop_type")
#                if cm.flovi_bmop_type == 'fixedValue':
#                    newrow(layout, "Pressure value:", cm, "flovi_b_sval")
#                newrow(layout, "Velocity sub-type:", cm, "flovi_bmou_type")
#                newrow(layout, "Field value:", cm, "flovi_u_field")
#                if not cm.flovi_u_field:
#                    newrow(layout, 'Velocity:', cm, 'flovi_b_vval')
#                if fvsimnode and fvsimnode.solver != 'icoFoam':
#                    newrow(layout, "nut type:", cm, "flovi_bmonut_type")
#                    if fvsimnode.turbulence == 'SpalartAllmaras':                        
#                        newrow(layout, "nuTilda type:", cm, "flovi_bmonutilda_type")
#                    elif fvsimnode.turbulence == 'kEpsilon':
#                        newrow(layout, "k type:", cm, "flovi_bmok_type")
#                        newrow(layout, "Epsilon type:", cm, "flovi_bmoe_type")
#                    elif fvsimnode.turbulence == 'kOmega':
#                        newrow(layout, "k type:", cm, "flovi_bmok_type")
#                        newrow(layout, "Omega type:", cm, "flovi_bmoo_type")
                
class VIObPanel(bpy.types.Panel):
    bl_label = "VI-Suite Object Definition"
    bl_space_type = "PROPERTIES"
    bl_region_type = "WINDOW"
    bl_context = "data"

    @classmethod
    def poll(cls, context):
        if context.object and context.object.type in ('LAMP', 'MESH'):
            return True

    def draw(self, context):
        obj = context.active_object
        layout = self.layout

        if obj.type == 'MESH':
            row = layout.row()
            row.prop(obj, "vi_type")
            if obj.vi_type == '1':
                row = layout.row()
                row.prop(obj, "envi_type")
                if obj.envi_type == '0':
                    newrow(layout, 'Inside convection:', obj, "envi_ica")
                    newrow(layout, 'Outside convection:', obj, "envi_oca")
                    
            if obj.vi_type == '3':
                newrow(layout, 'Feature level:', obj, "flovi_fl")
                newrow(layout, 'Surface max level:', obj, "flovi_slmax")
                newrow(layout, 'Surface min level:', obj, "flovi_slmin")
                newrow(layout, 'Surface layers:', obj, "flovi_slmin")
                
        if (obj.type == 'LAMP' and obj.data.type != 'SUN') or obj.vi_type == '4':
            newrow(layout, 'IES file:', obj, "ies_name")
            newrow(layout, 'IES Dimension:', obj, "ies_unit")
            newrow(layout, 'IES Strength:', obj, "ies_strength")
            newrow(layout, 'IES Colour:', obj, "ies_colmenu")
            if obj.ies_colmenu == '0':
                newrow(layout, 'IES RGB:', obj, "ies_rgb")
            else:
                newrow(layout, 'IES Temperature:', obj, "ies_ct")

        elif obj.vi_type == '5':                
            newrow(layout, 'Direction:', obj, 'li_bsdf_direc')
            newrow(layout, 'Proxy:', obj, 'li_bsdf_proxy')
            newrow(layout, 'Klems/Tensor:', obj, 'li_bsdf_tensor')
            if obj.li_bsdf_tensor != ' ':
                newrow(layout, 'resolution:', obj, 'li_bsdf_res')
                newrow(layout, 'Samples:', obj, 'li_bsdf_tsamp')
            else:
                newrow(layout, 'Samples:', obj, 'li_bsdf_ksamp')
            newrow(layout, 'RC params:', obj, 'li_bsdf_rcparam')
            
            if any([obj.data.materials[i].radmatmenu == '8' for i in [f.material_index for f in obj.data.polygons]]):
                row = layout.row()
                row.operator("object.gen_bsdf", text="Generate BSDF")

def rmmenu(layout, cm):
    row = layout.row()
    row.label('LiVi Radiance type:')
    row.prop(cm, 'radmatmenu')
    row = layout.row()

    for prop in cm.radmatdict[cm.radmatmenu]:
        if prop:
             row.prop(cm, prop)
        else:
            row = layout.row()
            
    if cm.radmatmenu == '8':
        newrow(layout, 'Proxy depth:', cm, 'li_bsdf_proxy_depth')
        row = layout.row()
        row.operator("material.load_bsdf", text="Load BSDF")
    elif cm.radmatmenu == '9':
        layout.prop_search(cm, 'radfile', bpy.data, 'texts', text='File', icon='TEXT')
    if cm.get('bsdf'):
        row.operator("material.del_bsdf", text="Delete BSDF")
        row = layout.row()
        row.operator("material.save_bsdf", text="Save BSDF")
    if cm.radmatmenu in ('1', '2', '3', '7'):
        newrow(layout, 'Photon port:', cm, 'pport')
    if cm.radmatmenu in ('0', '1', '2', '3', '6'):
        newrow(layout, 'Textured:', cm, 'radtex')
        if cm.radtex:
            newrow(layout, 'Normal map:', cm, 'radnorm')
            if cm.radnorm:
                newrow(layout, 'Strength:', cm, 'ns')
                newrow(layout, 'Image green vector:', cm, 'nu')
                newrow(layout, 'Image red vector:', cm, 'nside')

    row = layout.row()
    row.label("-----------------------------------------")
    
class MESH_Gridify_Panel(bpy.types.Panel):
    bl_label = "VI Gridify"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_context = "mesh_edit"
    bl_category = "VI-Suite"
      
    def draw(self, context):         
#        scene = context.scene
        layout = self.layout
#        newrow(layout, 'Rotation:', scene, 'vi_gridify_rot')
#        newrow(layout, 'Size 1:', scene, 'vi_gridify_us')
#        newrow(layout, 'Size 2:', scene, 'vi_gridify_as')
        row = layout.row()
        row.operator("object.vi_gridify2", text="Grid the object")