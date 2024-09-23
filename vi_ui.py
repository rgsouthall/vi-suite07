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

import bpy
from collections import OrderedDict
from .vi_func import newrow, retdates, logentry, get_materials


class VI_PT_3D(bpy.types.Panel):
    '''VI-Suite 3D view panel'''
    bl_label = "VI Display"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "VI-Suite"

    def draw(self, context):
        scene = context.scene
        svp = scene.vi_params
        cao = context.active_object
        layout = self.layout

        if cao:
            covp = cao.vi_params

        if cao and cao.active_material and cao.active_material.vi_params.get('bsdf'):
            if cao.active_material.vi_params['bsdf'].get('type') and cao.active_material.vi_params['bsdf']['type'] == 'LBNL/Klems Full' and not svp.vi_display:
                row = layout.row()
                row.operator("view3d.bsdf_display", text="BSDF Display")

            if svp['viparams'].get('vidisp') and svp['viparams']['vidisp'] == 'bsdf_panel':
                newrow(layout, 'Direction:', svp, "vi_bsdf_direc")
                newrow(layout, 'BSDF max:', svp, "vi_bsdfleg_max")
                newrow(layout, 'BSDF min:', svp, "vi_bsdfleg_min")
                newrow(layout, 'BSDF scale:', svp, "vi_bsdfleg_scale")
                newrow(layout, 'BSDF colour:', svp, "vi_leg_col")
                newrow(layout, 'BSDF font:', svp, "vi_bsdf_font")

        if svp.get('viparams') and svp['viparams'].get('vidisp'):
            if not svp.vi_display and svp['viparams']['vidisp'] == 'wr' and 'Wind_Plane' in [o.vi_params['VIType'] for o in bpy.data.objects if o.vi_params.get('VIType')]:
                row = layout.row()
                row.operator('view3d.wrdisplay', text='Wind Metrics')

            elif svp['viparams']['vidisp'] == 'sp' and svp.vi_display:
                row = layout.row()
                row.prop(context.space_data.shading, "light")
                newrow(layout, "Latitude:", svp, 'latitude')
                newrow(layout, "Longitude:", svp, 'longitude')
                (sdate, edate) = retdates(svp.sp_sd, 365, 2015)
                time_disps = ((("Day of year: {}/{}".format(sdate.day, sdate.month), "sp_sd"),
                               ("Time of day: {}:{}".format(int(svp.sp_sh), int((svp.sp_sh*60) % 60)), "sp_sh")),
                              [("Time of day: {}:{}".format(int(svp.sp_sh), int((svp.sp_sh*60)) % 60), "sp_sh")],
                              [("Day of year: {}/{}".format(sdate.day, sdate.month), "sp_sd")])

                for i in time_disps[int(svp['spparams']['suns'])]:
                    newrow(layout, i[0], svp, i[1])

                for i in (("Sun strength:", "sp_sun_strength"), ("Sun angle:", "sp_sun_angle")):
                    newrow(layout, i[0], svp, i[1])

                newrow(layout, "Line width:", svp, 'sp_line_width')
                newrow(layout, "Solstice colour:", svp, 'sp_season_main')
                newrow(layout, "Hour main colour:", svp, 'sp_hour_main')
                newrow(layout, "Hour dash colour:", svp, 'sp_hour_dash')
                newrow(layout, "Hour dash ratio:", svp, 'sp_hour_dash_ratio')
                newrow(layout, "Hour dash density:", svp, 'sp_hour_dash_density')
                newrow(layout, "Sun size:", svp, 'sp_sun_size')
                newrow(layout, "Sun colour:", svp, 'sp_sun_colour')
                newrow(layout, "Globe colour:", svp, 'sp_globe_colour')

                time_disps = ((("Display time:", "sp_td"), ("Display hours:", "sp_hd")), [("Display hours:", "sp_hd")], [("Display hours:", "sp_hd")])

                for i in time_disps[int(svp['spparams']['suns'])]:
                    newrow(layout, i[0], svp, i[1])

                if (svp['spparams']['suns'] == '0' and (svp.sp_td or svp.sp_hd)) or svp.sp_hd:
                    for i in (("Font size:", "vi_display_rp_fs"), ("Font colour:", "vi_display_rp_fc"), ("Font shadow:", "vi_display_rp_sh")):
                        newrow(layout, i[0], svp, i[1])

                    if svp.vi_display_rp_sh:
                        newrow(layout, "Shadow colour:", svp, "vi_display_rp_fsh")

            elif svp['viparams']['vidisp'] in ('svf', 'ss', 'li', 'lc', 'rt'):
                if not svp.vi_display:
                    newrow(layout, "Process", svp, "vi_disp_process")
                    row = layout.row()

                    if svp['viparams']['vidisp'] == 'svf':
                        row.operator("view3d.svfdisplay", text="Sky View Display")
                    elif svp['viparams']['vidisp'] == 'ss':
                        row.operator("view3d.ssdisplay", text="Shadow Display")
                    elif svp['viparams']['vidisp'] == 'li':
                        row.operator("view3d.libd", text="Radiance Display")
                    elif svp['viparams']['vidisp'] == 'rt':
                        row.operator("view3d.rtdisplay", text="Acoustics Display")

                elif [o for o in bpy.data.objects if o.vi_params.vi_type_string == 'LiVi Res']:
                    if not svp.ss_disp_panel:
                        newrow(layout, 'Result type:', svp, "li_disp_menu")

                        if svp.li_disp_menu in ('ago1v', 'aga1v'):
                            newrow(layout, 'GA view:', svp, "vi_views")

                        elif svp['viparams']['vidisp'] == 'rt':
                            newrow(layout, 'Source:', svp, "au_sources")

                        newrow(layout, 'Legend unit:', svp, "vi_leg_unit")
                        newrow(layout, 'Processing:', svp, "vi_res_process")

                    if svp.vi_res_process == '1':
                        newrow(layout, 'Modifier:', svp, "vi_res_mod")
                    elif svp.vi_res_process == '2':
                        layout.prop_search(svp, 'script_file', bpy.data, 'texts', text='File', icon='TEXT')

                    newrow(layout, 'Frame:', svp, "vi_frames")
                    newrow(layout, 'Legend max:', svp, "vi_leg_max")
                    newrow(layout, 'Legend min:', svp, "vi_leg_min")
                    newrow(layout, 'Legend scale:', svp, "vi_leg_scale")
                    newrow(layout, 'Legend colour:', svp, "vi_leg_col")
                    newrow(layout, 'Legend levels:', svp, "vi_leg_levels")
                    newrow(layout, 'Emitter materials:', svp, "vi_disp_mat")

                    if svp.vi_disp_mat:
                        newrow(layout, 'Emitter strength:', svp, "vi_disp_ems")

                    if cao and cao.type == 'MESH':
                        newrow(layout, 'Draw wire:', svp, 'vi_disp_wire')

                    if svp.vi_disp_process == "1":
                        newrow(layout, "3D Level:", svp, "vi_disp_3dlevel")

                    newrow(layout, "Transparency:", svp, "vi_disp_trans")

                    if svp.vi_disp_process == "3":
                        newrow(layout, "Arrow size:", svp, "vi_arrow_size")
                    
                    if context.mode != "EDIT" and svp.vi_disp_process != "2":
                        row = layout.row()
                        row.label(text="{:-<48}".format("Point visualisation "))
                        newrow(layout, 'Enable:', svp, 'vi_display_rp')

                        if svp.vi_display_rp:
                            propdict = OrderedDict([("Selected only:", "vi_display_sel_only"), ("Visible only:", "vi_display_vis_only"),
                                                    ("Font size:", "vi_display_rp_fs"), ("Font colour:", "vi_display_rp_fc"), ("Font shadow:", "vi_display_rp_sh"),
                                                    ("Shadow colour:", "vi_display_rp_fsh"), ("Position offset:", "vi_display_rp_off")])

                            for prop in propdict.items():
                                if svp.vi_display_rp_sh or prop[0] != "Shadow colour:":
                                    newrow(layout, prop[0], svp, prop[1])

                        row = layout.row()
                        row.label(text="{:-<60}".format(""))

                    elif context.mode != "EDIT" and svp.vi_disp_process == "2":
                        newrow(layout, 'Placement', svp, "vi_disp_pos")

            if svp.vi_display:
                newrow(layout, 'Display active', svp, 'vi_display')


class VI_PT_Mat(bpy.types.Panel):
    bl_label = "VI-Suite Material"
    bl_space_type = "PROPERTIES"
    bl_region_type = "WINDOW"
    bl_context = "material"

    def draw(self, context):
        cm, scene = context.material, context.scene

        if cm:
            svp = scene.vi_params
            mvp = cm.vi_params
            layout = self.layout
            newrow(layout, 'Material type:', mvp, "mattype")

            if mvp.mattype in ('0', '1'):
                rmmenu(layout, cm)

                if mvp.mattype == '0':
                    if not mvp.envi_nodes or (mvp.envi_nodes.name != cm.name and mvp.envi_nodes.name in [m.name for m in bpy.data.materials]):
                        row = layout.row()
                        row.operator("material.envi_node", text="Create EnVi Nodes")

            elif mvp.mattype == '2':
                if not context.scene.vi_params.get('flparams') or not context.scene.vi_params['flparams'].get('scenario'):
                    row = layout.row()
                    row.label(text='Export FloVi case')
                else:
                    newrow(layout, "Netgen max cell size:", mvp, "flovi_ng_max")
                    newrow(layout, "Type:", mvp, "flovi_bmb_type")

                    if svp.get('flparams') and svp['flparams'].get('solver_type'):
                        newrow(layout, "Pressure type:", mvp, "flovi_bmbp_subtype")

                        if mvp.flovi_bmbp_subtype in ('fixedValue', 'totalPressure'):
                            if mvp.flovi_bmbp_subtype == 'totalPressure':
                                newrow(layout, "p0 value:", mvp, "flovi_bmbp_p0val")

                            newrow(layout, "Pressure field value:", mvp, "flovi_p_field")

                            if not mvp.flovi_p_field:
                                newrow(layout, "Pressure value:", mvp, "flovi_bmbp_val")

                        newrow(layout, "U type:", mvp, "flovi_bmbu_subtype")

                        if mvp.flovi_bmbu_subtype in ('fixedValue', 'pressureInletOutletVelocity', 'inletOutlet', 'atmBoundaryLayerInletVelocity', 'outletInlet', 'freestream'):
                            newrow(layout, "U field value:", mvp, "flovi_u_field")

                            if mvp.flovi_bmbu_subtype == 'atmBoundaryLayerInletVelocity':
                                newrow(layout, "Reference height:", mvp, "flovi_u_zref")
                                newrow(layout, "Up vector:", mvp, "flovi_u_zdir")
                                newrow(layout, "Roughness:", mvp, "flovi_u_z0")
                                newrow(layout, "Ground height:", mvp, "flovi_u_zground")
                                newrow(layout, "Displacement:", mvp, "flovi_u_d")

                                if not mvp.flovi_u_field:
                                    newrow(layout, "Reference speed:", mvp, "flovi_u_uref")
                                    newrow(layout, "Flow direction:", mvp, "flovi_u_fdir")

                            if not mvp.flovi_u_field:
                                newrow(layout, "U direction:", mvp, "flovi_u_type")

                                if mvp.flovi_u_type == '0':
                                    newrow(layout, "U value:", mvp, "flovi_bmbu_val")
                                else:
                                    newrow(layout, "Azimuth:", mvp, "flovi_u_azi")
                                    newrow(layout, "Speed:", mvp, "flovi_u_speed")

                        if svp.get('flparams') and svp['flparams'].get('params'):
                            newrow(layout, "Nut type:", mvp, "flovi_bmbnut_subtype")

                            if mvp.flovi_bmbnut_subtype == 'fixedValue':
                                newrow(layout, "Nut field value:", mvp, "flovi_nut_field")

                                if not mvp.flovi_u_field:
                                    newrow(layout, "Nut value:", mvp, "flovi_bmbnut_val")

                            if 'k' in svp['flparams']['params']:
                                newrow(layout, "k type:", mvp, "flovi_k_subtype")

                                if mvp.flovi_k_subtype == 'fixedValue':
                                    newrow(layout, "k field value:", mvp, "flovi_k_field")

                                    if not mvp.flovi_k_field:
                                        newrow(layout, "k value:", mvp, "flovi_k_val")

                                if mvp.flovi_k_subtype == 'inletOutlet':
                                    newrow(layout, "k field:", mvp, "flovi_k_field")

                                    if not mvp.flovi_k_field:
                                        newrow(layout, "k inlet value:", mvp, "flovi_k_val")

                                elif mvp.flovi_k_subtype == 'turbulentIntensityKineticEnergyInlet':
                                    newrow(layout, "k intensity:", mvp, "flovi_k_intensity")
                                    newrow(layout, "k field value:", mvp, "flovi_k_field")

                                    if not mvp.flovi_k_field:
                                        newrow(layout, "k value:", mvp, "flovi_k_val")

                                newrow(layout, "Epsilon type:", mvp, "flovi_bmbe_subtype")

                                if mvp.flovi_bmbe_subtype == 'fixedValue':
                                    newrow(layout, "Epsilon field value:", mvp, "flovi_e_field")

                                    if not mvp.flovi_e_field:
                                        newrow(layout, "Epsilon value:", mvp, "flovi_bmbe_val")

                                elif mvp.flovi_bmbe_subtype == 'inletOutlet':
                                    newrow(layout, "Epsilon field value:", mvp, "flovi_e_field")

                                    if not mvp.flovi_e_field:
                                        newrow(layout, "Epsilon inlet value:", mvp, "flovi_bmbe_val")

                                elif mvp.flovi_bmbe_subtype == 'turbulentMixingLengthDissipationRateInlet':
                                    newrow(layout, "Epsilon ml:", mvp, "flovi_eml_val")

                            if svp['flparams']['solver_type'] == 'bf':
                                newrow(layout, "T type:", mvp, "flovi_bmbt_subtype")

                                if mvp.flovi_bmbt_subtype == 'fixedValue':
                                    newrow(layout, "T field:", mvp, "flovi_t_field")

                                    if not mvp.flovi_t_field:
                                        newrow(layout, "T value:", mvp, "flovi_bmbt_val")

                                elif mvp.flovi_bmbt_subtype == 'inletOutlet':
                                    newrow(layout, "T field:", mvp, "flovi_t_field")

                                    if not mvp.flovi_t_field:
                                        newrow(layout, "T inlet value:", mvp, "flovi_bmbti_val")
                                        newrow(layout, "T value:", mvp, "flovi_bmbt_val")

                                newrow(layout, "p_rgh type:", mvp, "flovi_prgh_subtype")
                                newrow(layout, "p_rgh field:", mvp, "flovi_prgh_field")

                                if mvp.flovi_prgh_subtype == 'totalPressure':
                                    newrow(layout, "p_rgh gamma:", mvp, "flovi_prgh_gamma")

                                    if not mvp.flovi_prgh_field:
                                        newrow(layout, "p0:", mvp, "flovi_prgh_p0")
                                        newrow(layout, "p_rgh value:", mvp, "flovi_prgh_val")
                                else:
                                    if not mvp.flovi_prgh_field:
                                        newrow(layout, "p_rgh p:", mvp, "flovi_prgh_p")
                                        newrow(layout, "p_rgh value:", mvp, "flovi_prgh_val")

                                newrow(layout, "alphat type:", mvp, "flovi_a_subtype")

                                if svp['flparams']['features']['rad']:
                                    newrow(layout, "Rad type:", mvp, "flovi_rad_subtype")
                                    if mvp.flovi_rad_subtype == 'MarshakRadiation':
                                        newrow(layout, "Emissivity mode:", mvp, "flovi_rad_em")
                                        newrow(layout, "Emissivity value:", mvp, "flovi_rad_e")
                                        newrow(layout, "Radiation value:", mvp, "flovi_rad_val")

                                    if svp['flparams']['radmodel'] == '1':
                                        newrow(layout, "I type:", mvp, "flovi_i_subtype")
                                        newrow(layout, "I mode:", mvp, "flovi_i_em")
                                        newrow(layout, "I emissivity:", mvp, "flovi_i_e")
                                        newrow(layout, "I value:", mvp, "flovi_i_val")

                    newrow(layout, "Probe:", mvp, "flovi_probe")

                    if mvp.flovi_bmb_type in ('Wall', 'Solid'):
                        newrow(layout, "HTC:", mvp, "flovi_htc")
                
            elif mvp.mattype == '3':
                if not mvp.am.updated:
                    mvp.am.update()

                newrow(layout, "AuVi absorption class:", mvp, "auvi_abs_class")

                if mvp.auvi_abs_class == '0':
                    newrow(layout, "Absorption type:", mvp, "auvi_type_abs")
                    newrow(layout, "Absorption material:", mvp, "auvi_mat_abs")

                elif mvp.auvi_abs_class == '1': 
                    newrow(layout, "Flat response:", mvp, "auvi_abs_flat")

                    if mvp.auvi_abs_flat:
                        newrow(layout, "Absorption coefficient:", mvp, "auvi_o1_abs")
                    else:
                        row = layout.row()
                        row.label(text="Absorption coefficients:")
                        newrow(layout, "125Hz", mvp, "auvi_o1_abs")
                        newrow(layout, "250Hz", mvp, "auvi_o2_abs")
                        newrow(layout, "500Hz", mvp, "auvi_o3_abs")
                        newrow(layout, "1000Hz", mvp, "auvi_o4_abs")
                        newrow(layout, "2000Hz", mvp, "auvi_o5_abs")
                        newrow(layout, "4000Hz", mvp, "auvi_o6_abs")
                        newrow(layout, "8000Hz", mvp, "auvi_o7_abs")
                    
                newrow(layout, "AuVi scatter class:", mvp, "auvi_scatt_class")
                
                if mvp.auvi_scatt_class == '0':
                    newrow(layout, "Scatter type:", mvp, "auvi_type_scatt")
                    newrow(layout, "Scatter material", mvp, "auvi_mat_scatt")
                else:
                    newrow(layout, "Flat response", mvp, "auvi_scatt_flat")

                    if mvp.auvi_scatt_flat:
                        newrow(layout, "Scatter coefficient:", mvp, "auvi_o1_scatt")
                    else:
                        row = layout.row()
                        row.label(text="Scatter coefficients")
                        newrow(layout, "125Hz", mvp, "auvi_o1_scatt")
                        newrow(layout, "250Hz", mvp, "auvi_o2_scatt")
                        newrow(layout, "500Hz", mvp, "auvi_o3_scatt")
                        newrow(layout, "1000Hz", mvp, "auvi_o4_scatt")
                        newrow(layout, "2000Hz", mvp, "auvi_o5_scatt")
                        newrow(layout, "4000Hz", mvp, "auvi_o6_scatt")
                        newrow(layout, "8000Hz", mvp, "auvi_o7_scatt")


class VI_PT_Ob(bpy.types.Panel):
    bl_label = "VI-Suite Object"
    bl_space_type = "PROPERTIES"
    bl_region_type = "WINDOW"
    bl_context = "object"

    @classmethod
    def poll(cls, context):
        if context.object and context.object.type in ('LIGHT', 'MESH', 'EMPTY'):
            return True

    def draw(self, context):
        obj = context.object
        ovp = obj.vi_params
        layout = self.layout

        if obj.type == 'MESH':
            row = layout.row()
            row.prop(ovp, 'vi_type')

            if ovp.vi_type == '0':
                row = layout.row()
                row.label(text='-- Octree generation --')
                newrow(layout, 'Triangulate:', ovp, 'triangulate')
                newrow(layout, 'Mesh:', ovp, 'mesh')
                
                if context.scene.vi_params.get('viparams'):
                    row = layout.row()
                    row.operator('object.vi_genoct', text="Generate Octree")

            elif ovp.vi_type == '2':
                pass

            newrow(layout, 'Embodied:', ovp, 'embodied')

            if ovp.embodied:
                newrow(layout, 'Embodied class:', ovp, 'embodiedclass')

                if ovp.embodiedclass != 'Custom':
                    newrow(layout, 'Embodied type:', ovp, 'embodiedtype')
                    newrow(layout, 'Embodied material:', ovp, 'embodiedmat')
                    newrow(layout, 'Embodied modifier:', ovp, 'ec_amount_mod')
                    newrow(layout, 'Service life:', ovp, 'ec_life')

                    if ovp.get('ecentries') and [ec[1] for ec in ovp['ecentries'] if ec[0] == 'unit'][0] == 'each':
                        newrow(layout, "Object represents:", ovp, "ec_rep")
                        
                        if ovp.ec_rep == '1':
                            newrow(layout, "No. of items:", ovp, "ec_items")

                    if ovp.get('ecentries'):
                        for ec in ovp['ecentries']:
                            row = layout.row()
                            row.label(text='{}: {}'.format(ec[0], ec[1]))

                    row = layout.row()
                    row.operator("object.ec_edit", text="Edit")

                else:                   
                    newrow(layout, "Embodied id:", ovp, "ec_id")
                    newrow(layout, "New class:", ovp, "ec_class")
                    newrow(layout, "Embodied type:", ovp, "ec_type")
                    newrow(layout, "Embodied modules:", ovp, "ec_mod")
                    newrow(layout, "Embodied name:", ovp, "ec_name")
                    newrow(layout, "Declared unit:", ovp, "ec_unit")
                    newrow(layout, "Amount of DU:", ovp, "ec_amount")
                    newrow(layout, "Modifier per amount:", ovp, "ec_amount_mod")
                    newrow(layout, "GWP per amount:", ovp, "ec_du")

                    if ovp.ec_unit not in ('kg', 'm3', 'tonnes'):
                        newrow(layout, "Embodied weight:", ovp, "ec_weight")

                    if ovp.ec_unit not in ('each'):
                        newrow(layout, "Embodied density:", ovp, "ec_density")

                    newrow(layout, "Service life:", ovp, "ec_life")

                    if ovp.get('ecentries') and [ec[1] for ec in ovp['ecentries'] if ec[0] == 'unit'][0] == 'each':
                        newrow(layout, "Object represents:", ovp, "ec_rep")

                        if ovp.ec_rep == '1':
                            newrow(layout, "No. of items:", ovp, "ec_items")

                    if all((ovp.ec_id, ovp.ec_type, ovp.ec_class, ovp.ec_name, ovp.ec_du, ovp.ec_density, ovp.ec_mod)):
                        row = layout.row()
                        row.operator("object.ec_save", text="Save")

        if (obj.type == 'LIGHT' and obj.data.type != 'SUN') or ovp.vi_type == '4':
            newrow(layout, 'IES file:', ovp, "ies_name")
            newrow(layout, 'IES Dimension:', ovp, "ies_unit")
            newrow(layout, 'IES Strength:', ovp, "ies_strength")
            newrow(layout, 'IES Colour:', ovp, "ies_colmenu")

            if ovp.ies_colmenu == '0':
                newrow(layout, 'IES RGB:', ovp, "ies_rgb")
            else:
                newrow(layout, 'IES Temperature:', ovp, "ies_ct")

        elif ovp.vi_type == '5':
            if any([obj.material_slots[i].material.vi_params.radmatmenu == '8' for i in [f.material_index for f in obj.data.polygons if f.material_index <= len(obj.material_slots)]]):
                newrow(layout, 'Direction:', ovp, 'li_bsdf_direc')
                newrow(layout, 'Klems/Tensor:', ovp, 'li_bsdf_tensor')

                if ovp.li_bsdf_tensor != ' ':
                    newrow(layout, 'resolution:', ovp, 'li_bsdf_res')
                    newrow(layout, 'Samples:', ovp, 'li_bsdf_tsamp')
                else:
                    newrow(layout, 'Samples:', ovp, 'li_bsdf_ksamp')
                newrow(layout, 'RC params:', ovp, 'li_bsdf_rcparam')

                if not ovp.bsdf_running:
                    row = layout.row()
                    row.operator("object.gen_bsdf", text="Generate BSDF")
            else:
                row = layout.row()
                row.label(text='No BSDF material applied')

        if obj.type == 'EMPTY' or (ovp.vi_type == '3' and obj.type == 'MESH' and len(obj.data.polygons) == 1):
            newrow(layout, 'CFD probe:', ovp, 'flovi_probe')
            
            if obj.type == 'EMPTY':
                newrow(layout, 'AuVi type:', ovp, 'auvi_sl')

def rmmenu(layout, cm):
    mvp = cm.vi_params
    row = layout.row()
    row.label(text='LiVi Radiance type:')
    row.prop(mvp, 'radmatmenu')
    row = layout.row()

    for prop in mvp.radmatdict[mvp.radmatmenu]:
        if prop:
            row.prop(mvp, prop)
        else:
            row = layout.row()

    if mvp.radmatmenu == '8':
        newrow(layout, 'Proxy depth:', mvp, 'li_bsdf_proxy_depth')
        newrow(layout, 'Up vector:', mvp, 'li_bsdf_up')
        row = layout.row()
        row.operator("material.load_bsdf", text="Load BSDF")
    
    elif mvp.radmatmenu == '9':
        layout.prop_search(mvp, 'radfile', bpy.data, 'texts', text='File', icon='TEXT')
    
    if mvp.get('bsdf'):
        row.operator("material.del_bsdf", text="Delete BSDF")
        row = layout.row()
        row.operator("material.save_bsdf", text="Save BSDF")
    
    if mvp.radmatmenu in ('1', '2', '3', '7'):
        newrow(layout, 'Photon port:', mvp, 'pport')
    
    if mvp.mattype == '0' and mvp.radmatmenu in ('0', '1', '2', '3', '6'):
        newrow(layout, 'Textured:', mvp, 'radtex')

        if mvp.radtex:
            newrow(layout, 'Dirt noise:', mvp, 'li_dirt')

            if mvp.li_dirt:
                newrow(layout, 'Dirt spacing:', mvp, 'li_dirt_spacing')
                newrow(layout, 'Dirt level:', mvp, 'li_dirt_level')

            newrow(layout, 'Texture image:', mvp, 'li_tex')

            if not mvp.li_am:
                newrow(layout, 'Normal map:', mvp, 'li_norm')

            if not mvp.li_norm:
                newrow(layout, 'Blend mask:', mvp, 'li_am')

                if mvp.li_am:
                    newrow(layout, 'Background image:', mvp, 'li_tex_black')

            else:
                newrow(layout, 'Strength:', mvp, 'li_norm_strength')
                newrow(layout, 'Image green vector:', mvp, 'nu')
                newrow(layout, 'Image red vector:', mvp, 'nside')

    row = layout.row()
    row.label(text="-----------------------------------------")


class VI_PT_Col(bpy.types.Panel):
    bl_label = "VI-Suite Collection"
    bl_space_type = "PROPERTIES"
    bl_region_type = "WINDOW"
    bl_context = "collection"

    @classmethod
    def poll(cls, context):
        if context.collection:
            return True

    def draw(self, context):
        coll = context.collection
        cvp = coll.vi_params
        layout = self.layout
        newrow(layout, 'EnVi zone:', context.collection.vi_params, "envi_zone")

        try:
            fa = '{:.2f}'.format(context.collection.vi_params['enparams']['floorarea'][str(context.scene.frame_current)])
        except Exception:
            fa = 'N/A'

        row = layout.row()
        row.label(text='EnergyPlus Metrics:')
        row = layout.row()
        row.label(text='Floor area (m2): {}'.format(fa))
        newrow(layout, 'Embodied:', cvp, 'embodied')

        if cvp.embodied:
            newrow(layout, 'Embodied class:', cvp, 'embodiedclass')

            if cvp.embodiedclass != 'Custom':
                newrow(layout, 'Embodied type:', cvp, 'embodiedtype')
                newrow(layout, 'Embodied material:', cvp, 'embodiedmat')
                newrow(layout, 'Embodied modifer:', cvp, 'ec_amount_mod')
                newrow(layout, 'Service life:', cvp, 'ec_life')

                if cvp.get('ecentries'):
                    for ec in cvp['ecentries']:
                        row = layout.row()
                        row.label(text='{}: {}'.format(ec[0], ec[1]))

                row = layout.row()
                row.operator("object.ec_edit", text="Edit")

            else:
                newrow(layout, "Embodied id:", cvp, "ec_id")
                newrow(layout, "New embodied class:", cvp, "ec_class")
                newrow(layout, "Embodied type:", cvp, "ec_type")
                newrow(layout, "Embodied modules:", cvp, "ec_mod")
                newrow(layout, "Embodied name:", cvp, "ec_name")
                newrow(layout, "Declared unit:", cvp, "ec_unit")
                newrow(layout, "Amount of DU:", cvp, "ec_amount")
                newrow(layout, "GWP per amount:", cvp, "ec_du")

                if cvp.ec_unit not in ('kg', 'm3', 'tonnes'):
                    newrow(layout, "Embodied weight:", cvp, "ec_weight")

                newrow(layout, "Embodied density:", cvp, "ec_density")
                newrow(layout, "Lifespan:", cvp, "ec_life")

                if all((cvp.ec_id, cvp.ec_type, cvp.ec_class, cvp.ec_name, cvp.ec_du, cvp.ec_density, cvp.ec_mod)):
                    row = layout.row()
                    row.operator("object.ec_save", text="Save")

class VI_PT_Gridify(bpy.types.Panel):
    bl_label = "VI Gridify"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_context = "mesh_edit"
    bl_category = "VI-Suite"

    def draw(self, context):
        layout = self.layout
        row = layout.row()
        row.operator("object.vi_gridify", text="Grid the object")


class TREE_PT_vi(bpy.types.Panel):
    bl_label = "VI-Suite"
    bl_space_type = "NODE_EDITOR"
    bl_region_type = "UI"
    bl_category = "VI-Suite"

    @classmethod
    def poll(cls, context):
        return context.space_data.tree_type in ('ViN', 'EnViN', 'EnViMatN')

    def draw(self, context):
        layout = self.layout
        col = layout.column(align=True)
        visuite_groups = [g for g in bpy.data.node_groups if g.bl_idname == 'ViN']

        for g in visuite_groups:
            emboss = False

            if len(context.space_data.path) > 0:
                emboss = context.space_data.path[-1].node_tree.name == g.name

            op = col.operator('tree.goto_group', text=g.name, emboss=emboss, icon='NODETREE')
            op.tree_type = "ViN"
            op.tree = g.name

        col.separator()
        col.separator()
        col.separator()


class TREE_PT_envim(bpy.types.Panel):
    bl_label = "EnVi Materials"
    bl_space_type = "NODE_EDITOR"
    bl_region_type = "UI"
    bl_category = "VI-Suite"

    @classmethod
    def poll(cls, context):
        return context.space_data.tree_type in ('ViN', 'EnViN', 'EnViMatN')

    def draw(self, context):
        layout = self.layout
        materials = get_materials()
        col = layout.column(align=True)

        for mat in materials:
            name = mat.name

            try:
                icon_val = layout.icon(mat)
            except:
                icon_val = 1

            if mat.users:
                op = col.operator('tree.goto_mat',
                                  text=name,
                                  emboss=(mat == context.space_data.id),
                                  icon_value=icon_val)
                op.mat = name
            else:
                row = col.row(align=True)
                op = row.operator('tree.goto_mat',
                                  text=name,
                                  emboss=(mat == context.space_data.id),
                                  icon_value=icon_val)
                op.mat = name
                op = row.operator('tree.goto_mat',
                                  text="",
                                  emboss=(mat == context.space_data.id),
                                  icon='ORPHAN_DATA')
                op.mat = name

        if not materials:
            col.label(text="No EnVi Materials")
        else:
            row = layout.row()
            row.operator('envi_node.remove', text='Clear unused', icon='NODETREE')


class TREE_PT_envin(bpy.types.Panel):
    bl_label = "EnVi Networks"
    bl_space_type = "NODE_EDITOR"
    bl_region_type = "UI"
    bl_category = "VI-Suite"

    @classmethod
    def poll(cls, context):
        return context.space_data.tree_type in ('ViN', 'EnViN', 'EnViMatN')

    def draw(self, context):
        layout = self.layout
        col = layout.column(align=True)
        envin_groups = [g for g in bpy.data.node_groups if g.bl_idname == 'EnViN']

        for g in envin_groups:
            emboss = False

            if len(context.space_data.path) > 0:
                emboss = context.space_data.path[-1].node_tree.name == g.name

            op = col.operator('tree.goto_group', text=g.name, emboss=emboss, icon='NODETREE')
            op.tree_type = "EnViN"
            op.tree = g.name

        if not envin_groups:
            col.label(text="No EnVi Network")

        col.separator()

