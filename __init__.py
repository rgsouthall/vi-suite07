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

bl_info = {
    "name": "VI-Suite",
    "author": "Ryan Southall",
    "version": (0, 6, 0),
    "blender": (2, 80, 0),
    "api":"",
    "location": "Node Editor & 3D View > Properties Panel",
    "description": "Radiance/EnergyPlus exporter and results visualiser",
    "warning": "This is a beta script. Some functionality is buggy",
    "wiki_url": "",
    "tracker_url": "",
    "category": "Import-Export"}

if "bpy" in locals():
    import imp
    imp.reload(vi_node)
    imp.reload(vi_operators)
    imp.reload(vi_ui)
    imp.reload(vi_func)
#    imp.reload(envi_mat)
else:
    from .vi_node import vinode_categories, ViNetwork, ViLoc, ViLocSock, ViSPNode
#    from .envi_mat import envi_materials, envi_constructions, envi_layero, envi_layer1, envi_layer2, envi_layer3, envi_layer4, envi_layerotype, envi_layer1type, envi_layer2type, envi_layer3type, envi_layer4type, envi_con_list
    from .vi_func import iprop, bprop, eprop, fprop, sprop, fvprop, sunpath1, radmat, radbsdf, retsv, cmap
    from .vi_func import rtpoints, lhcalcapply, udidacalcapply, compcalcapply, basiccalcapply, lividisplay, setscenelivivals, py_path
#    from .envi_func import enunits, enpunits, enparametric, resnameunits, aresnameunits
#    from .flovi_func import fvmat, ret_fvbp_menu, ret_fvbu_menu, ret_fvbnut_menu, ret_fvbnutilda_menu, ret_fvbk_menu, ret_fvbepsilon_menu, ret_fvbomega_menu, ret_fvbt_menu, ret_fvba_menu, ret_fvbprgh_menu
#    from .vi_display import setcols
    from .vi_operators import NODE_OT_SunPath, VIEW3D_OT_SPNumDisplay
    from .vi_operators import *
    from .vi_ui import VI_PT_3D
    from .vi_ui import *

import sys, os, inspect, bpy, nodeitems_utils, bmesh, math, mathutils
from bpy.app.handlers import persistent
from numpy import array, digitize, logspace, multiply
from numpy import log10 as nlog10
from bpy.props import StringProperty, EnumProperty, IntProperty, FloatVectorProperty, FloatProperty
from bpy.types import AddonPreferences

evsep = {'linux': ':', 'darwin': ':', 'win32': ';'}
platpath = {'linux': ':', 'darwin': ':', 'win32': ';'}
addonpath = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

if sys.platform in ('darwin', 'win32'):
    if 'PYTHONPATH' not in os.environ:
        os.environ['PYTHONPATH'] =  os.path.join(addonpath, 'Python')
    elif os.path.join(addonpath, 'Python') not in os.environ['PYTHONPATH']:
        os.environ['PYTHONPATH'] =  os.environ['PYTHONPATH'] + evsep[str(sys.platform)] + os.path.join(addonpath, 'Python')

def return_preferences():
    return bpy.context.preferences.addons[__name__].preferences

def abspath(self, context):
    if self.radbin != bpy.path.abspath(self.radbin):
        self.radbin = bpy.path.abspath(self.radbin)
    if self.radlib != bpy.path.abspath(self.radlib):
        self.radlib = bpy.path.abspath(self.radlib)
    if self.epbin != bpy.path.abspath(self.epbin):
        self.epbin = bpy.path.abspath(self.epbin)
    if self.epweath != bpy.path.abspath(self.epweath):
        self.epweath = bpy.path.abspath(self.epweath)
    if self.ofbin != bpy.path.abspath(self.ofbin):
        self.ofbin = bpy.path.abspath(self.ofbin)
    if self.oflib != bpy.path.abspath(self.oflib):
        self.oflib = bpy.path.abspath(self.oflib)  
    if self.ofetc != bpy.path.abspath(self.ofetc):
        self.ofetc = bpy.path.abspath(self.ofetc)
        
class VIPreferences(AddonPreferences):
    bl_idname = __name__    
    radbin: StringProperty(name = '', description = 'Radiance binary directory location', default = '', subtype='DIR_PATH', update=abspath)
    radlib: StringProperty(name = '', description = 'Radiance library directory location', default = '', subtype='DIR_PATH', update=abspath)
    epbin: StringProperty(name = '', description = 'EnergyPlus binary directory location', default = '', subtype='DIR_PATH', update=abspath)
    epweath: StringProperty(name = '', description = 'EnergyPlus weather directory location', default = '', subtype='DIR_PATH', update=abspath)
    ofbin: StringProperty(name = '', description = 'OpenFOAM binary directory location', default = '', subtype='DIR_PATH', update=abspath)
    oflib: StringProperty(name = '', description = 'OpenFOAM library directory location', default = '', subtype='DIR_PATH', update=abspath)
    ofetc: StringProperty(name = '', description = 'OpenFOAM letc directory location', default = '', subtype='DIR_PATH', update=abspath)
    ui_dict = {"Radiance bin directory:": 'radbin', "Radiance lib directory:": 'radlib', "EnergyPlus bin directory:": 'epbin',
               "EnergyPlus weather directory:": 'epweath', 'OpenFOAM bin directory': 'ofbin', 'OpenFOAM lib directory': 'oflib', 'OpenFOAM etc directory': 'ofetc'}

    def draw(self, context):
        layout = self.layout 
        
        for entry in self.ui_dict:
            row = layout.row()
            row.label(text=entry)   
            row.prop(self, self.ui_dict[entry])

class VI_Params(bpy.types.PropertyGroup):   
    sp_hour_dash: fvprop(4, "",'Main colour of the hour lines', [1.0, 1.0, 1.0, 1.0], 'COLOR', 0, 1) 
    sp_hour_main: fvprop(4, "",'Dash colour of the hour lines', [1.0, 1.0, 1.0, 1.0], 'COLOR', 0, 1)
    sp_season_main: fvprop(4, "",'Main colour of the season lines', [1.0, 1.0, 1.0, 1.0], 'COLOR', 0, 1)
    sp_season_dash: fvprop(4, "",'Dash colour of the season lines', [1.0, 1.0, 1.0, 1.0], 'COLOR', 0, 1)
    sp_sun_colour: fvprop(4, "",'Sun colour', [0.0, 1.0, 1.0, 1.0], 'COLOR', 0, 1)
    sp_sun_size: iprop("",'Sun size', 1, 50, 10)
    sp_season_dash_ratio: fprop("", "Ratio of line to dash of season lines", 0, 5, 0)
    sp_hour_dash_ratio: fprop("", "Ratio of line to dash of hour lines", -1, 1, -0.5)
    sp_hour_dash_density: fprop("", "Ratio of line to dash of hour lines", 0, 5, 0)
    sp_line_width: iprop("", "Sun path line width", 0, 50, 4)
    latitude: fprop("", "Site decimal latitude (N is positive)", -89.99, 89.99, 52.0)
    longitude: fprop("", "Site decimal longitude (E is positive)", -180, 180, 0.0)
    sp_suns: EnumProperty(items = [('0', 'Single', 'Single sun'), ('1', 'Monthly', 'Monthly sun for chosen time'), ('2', 'Hourly', 'Hourly sun for chosen date')], name = '', description = 'Sunpath sun type', default = '0', update=sunpath1)
    sp_sst: FloatProperty(name = "", description = "Sun strength", min = 0, max = 100, default = 0.1, update=sunpath1)
    sp_ssi: FloatProperty(name = "", description = "Sun size", min = 0, max = 1, default = 0.01, update=sunpath1)
    sp_sd: IntProperty(name = "", description = "Day of year", min = 1, max = 365, default = 1, update=sunpath1)
    sp_sh: FloatProperty(name = "", description = "Time of day", subtype='TIME', unit='TIME', min = 0, max = 24, default = 12, update=sunpath1)
    sp_hd: bprop("", "",0)
    sp_up: bprop("", "",0)
    sp_td: bprop("", "",0)
    vi_li_disp_panel: iprop("Display Panel", "Shows the Display Panel", -1, 2, 0)
    
@persistent
def update_chart_node(dummy):
    try:
        for ng in [ng for ng in bpy.data.node_groups if ng.bl_idname == 'ViN']:
            [node.update() for node in ng.nodes if node.bl_label == 'VI Chart']
    except Exception as e:
        print('Chart node update failure:', e)

@persistent        
def update_dir(dummy):
    if bpy.context.scene.get('viparams'):
        fp = bpy.data.filepath
        bpy.context.scene['viparams']['newdir'] = os.path.join(os.path.dirname(fp), os.path.splitext(os.path.basename(fp))[0])
               
@persistent
def display_off(dummy):
    if bpy.context.scene.get('viparams') and bpy.context.scene['viparams'].get('vidisp'):
        
        ifdict = {'sspanel': 'ss', 'lipanel': 'li', 'enpanel': 'en', 'bsdf_panel': 'bsdf'}
        if bpy.context.scene['viparams']['vidisp'] in ifdict:
            bpy.context.scene['viparams']['vidisp'] = ifdict[bpy.context.scene['viparams']['vidisp']]
        bpy.context.scene.vi_display = 0
        
@persistent
def mesh_index(dummy):
    try:
        cao = bpy.context.active_object
    
        if cao and cao.layers[1] and cao.mode == 'EDIT':
            if not bpy.app.debug:
                bpy.app.debug = True
        elif bpy.app.debug:
            bpy.app.debug = False
    except:
        pass
            
@persistent
def select_nodetree(dummy): 
    for space in getViEditorSpaces():
        vings = [ng for ng in bpy.data.node_groups if ng.bl_idname == 'ViN']
        if vings:
            space.node_tree = vings[0]

    for space in getEnViEditorSpaces():
        envings = [ng for ng in bpy.data.node_groups if ng.bl_idname == 'EnViN']
        if envings:
            space.node_tree = envings[0]
            
    for space in getEnViMaterialSpaces():
        try:
            if space.node_tree != bpy.context.active_object.active_material.envi_nodes:
                envings = [ng for ng in bpy.data.node_groups if ng.bl_idname == 'EnViMatN' and ng == bpy.context.active_object.active_material.envi_nodes]
                if envings:
                    space.node_tree = envings[0]
        except Exception as e:
            pass
        
#bpy.app.handlers.scene_update_post.append(select_nodetree)
        
def getViEditorSpaces():
    if bpy.context.screen:
        return [area.spaces.active for area in bpy.context.screen.areas if area and area.type == "NODE_EDITOR" and area.spaces.active.tree_type == "ViN" and not area.spaces.active.edit_tree]
    else:
        return []
        
def getEnViEditorSpaces():
    if bpy.context.screen:
        return [area.spaces.active for area in bpy.context.screen.areas if area and area.type == "NODE_EDITOR" and area.spaces.active.tree_type == "EnViN" and not area.spaces.active.edit_tree]
    else:
        return []

def getEnViMaterialSpaces():
    if bpy.context.screen:        
        return [area.spaces.active for area in bpy.context.screen.areas if area and area.type == "NODE_EDITOR" and area.spaces.active.tree_type == "EnViMatN"]
    else:
        return []
        
#bpy.app.handlers.scene_update_post.append(select_nodetree)
#bpy.app.handlers.scene_update_post.append(mesh_index)
            
epversion = "8-9-0"
#envi_mats, envi_cons, conlayers = envi_materials(), envi_constructions(), 5

def path_update():
    vi_prefs = bpy.context.preferences.addons[__name__].preferences
    epdir = vi_prefs.epbin if vi_prefs and vi_prefs.epbin and os.path.isdir(vi_prefs.epbin) else os.path.join('{}'.format(addonpath), 'EPFiles')
    radldir = vi_prefs.radlib if vi_prefs and os.path.isdir(vi_prefs.radlib) else os.path.join('{}'.format(addonpath), 'Radfiles', 'lib')
    radbdir = vi_prefs.radbin if vi_prefs and os.path.isdir(vi_prefs.radbin) else os.path.join('{}'.format(addonpath), 'Radfiles', 'bin') 
    ofbdir = vi_prefs.ofbin if vi_prefs and os.path.isdir(vi_prefs.ofbin) else os.path.join('{}'.format(addonpath), 'OFFiles', 'bin') 
    ofldir = vi_prefs.oflib if vi_prefs and os.path.isdir(vi_prefs.oflib) else os.path.join('{}'.format(addonpath), 'OFFiles', 'lib')
    ofedir = vi_prefs.ofetc if vi_prefs and os.path.isdir(vi_prefs.ofetc) else os.path.join('{}'.format(addonpath), 'OFFiles')
    os.environ["PATH"] += "{0}{1}".format(evsep[str(sys.platform)], os.path.dirname(bpy.app.binary_path))

    if not os.environ.get('RAYPATH'):# or radldir not in os.environ['RAYPATH'] or radbdir not in os.environ['PATH']  or epdir not in os.environ['PATH']:
        if vi_prefs and os.path.isdir(vi_prefs.radlib):
            os.environ["RAYPATH"] = '{0}{1}{2}'.format(radldir, evsep[str(sys.platform)], os.path.join(addonpath, 'Radfiles', 'lib'))
        else:
            os.environ["RAYPATH"] = radldir
           
        os.environ["PATH"] = os.environ["PATH"] + "{0}{1}{0}{2}{0}{3}".format(evsep[str(sys.platform)], radbdir, epdir, ofbdir)    
        os.environ["LD_LIBRARY_PATH"] = os.environ["LD_LIBRARY_PATH"] + "{0}{1}".format(evsep[str(sys.platform)], ofldir) if os.environ.get("LD_LIBRARY_PATH") else "{0}{1}".format(evsep[str(sys.platform)], ofldir)
        os.environ["WM_PROJECT_DIR"] = ofedir
        
def colupdate(self, context):
    cmap(self)

def eupdate(self, context):
    scene = context.scene
    maxo, mino = scene.vi_leg_max, scene.vi_leg_min
    odiff = scene.vi_leg_max - scene.vi_leg_min

    if context.active_object.mode == 'EDIT':
        return
    if odiff:      
        for frame in range(scene['liparams']['fs'], scene['liparams']['fe'] + 1):
            for o in [obj for obj in bpy.data.objects if obj.lires == 1 and obj.data.shape_keys and str(frame) in [sk.name for sk in obj.data.shape_keys.key_blocks]]:  
                bm = bmesh.new()
                bm.from_mesh(o.data)  
                bm.transform(o.matrix_world)            
                skb = bm.verts.layers.shape['Basis']
                skf = bm.verts.layers.shape[str(frame)]
                
                if str(frame) in o['omax']:
                    if bm.faces.layers.float.get('res{}'.format(frame)):
                        extrude = bm.faces.layers.int['extrude']
                        res = bm.faces.layers.float['res{}'.format(frame)] #if context.scene['cp'] == '0' else bm.verts.layers.float['res{}'.format(frame)]                
                        faces = [f for f in bm.faces if f[extrude]]
                        fnorms = array([f.normal.normalized() for f in faces]).T
                        fres = array([f[res] for f in faces])
                        extrudes = (0.1 * scene.vi_disp_3dlevel * (nlog10(maxo * (fres + 1 - mino)/odiff)) * fnorms).T if scene.vi_leg_scale == '1' else \
                            multiply(fnorms, scene.vi_disp_3dlevel * ((fres - mino)/odiff)).T

                        for f, face in enumerate(faces):
                            for v in face.verts:
                                v[skf] = v[skb] + mathutils.Vector(extrudes[f])
                    
                    elif bm.verts.layers.float.get('res{}'.format(frame)):
                        res = bm.verts.layers.float['res{}'.format(frame)]
                        vnorms = array([v.normal.normalized() for v in bm.verts]).T
                        vres = array([v[res] for v in bm.verts])
                        extrudes = multiply(vnorms, scene.vi_disp_3dlevel * ((vres-mino)/odiff)).T if scene.vi_leg_scale == '0' else \
                            [0.1 * scene.vi_disp_3dlevel * (math.log10(maxo * (v[res] + 1 - mino)/odiff)) * v.normal.normalized() for v in bm.verts]  
                        for v, vert in enumerate(bm.verts):
                            vert[skf] = vert[skb] + mathutils.Vector(extrudes[v])

                bm.transform(o.matrix_world.inverted())
                bm.to_mesh(o.data)
                bm.free()

def tupdate(self, context):
    for o in [o for o in context.scene.objects if o.type == 'MESH'  and 'lightarray' not in o.name and o.hide == False and o.layers[context.scene.active_layer] == True and o.get('lires')]:
        o.show_transparent = 1
    for mat in [bpy.data.materials['{}#{}'.format('vi-suite', index)] for index in range(1, context.scene.vi_leg_levels + 1)]:
        mat.use_transparency, mat.transparency_method, mat.alpha = 1, 'MASK', context.scene.vi_disp_trans
    cmap(self)
        
def wupdate(self, context):
    o = context.active_object
    if o and o.type == 'MESH':
        (o.show_wire, o.show_all_edges) = (1, 1) if context.scene.vi_disp_wire else (0, 0)

def legupdate(self, context):
    scene = context.scene
    frames = range(scene['liparams']['fs'], scene['liparams']['fe'] + 1)
    obs = [o for o in scene.objects if o.get('lires')]
    increment = 1/scene.vi_leg_levels
    
    if scene.vi_leg_scale == '0':
        bins = array([increment * i for i in range(1, scene.vi_leg_levels)])
        
    elif scene.vi_leg_scale == '1':
        slices = logspace(0, 2, scene.vi_leg_levels + 1, True)
        bins = array([(slices[i] - increment * (scene.vi_leg_levels - i))/100 for i in range(scene.vi_leg_levels + 1)])
        bins = array([1 - math.log10(i)/math.log10(scene.vi_leg_levels + 1) for i in range(1, scene.vi_leg_levels + 2)][::-1])
        bins = bins[1:-1]
    
    for o in obs:
        bm = bmesh.new()
        bm.from_mesh(o.data)
        cmap(self)
        
        if len(o.material_slots) != scene.vi_leg_levels + 1:
            for matname in ['{}#{}'.format('vi-suite', i) for i in range(0, scene.vi_leg_levels + 1)]:
                if bpy.data.materials[matname] not in o.data.materials[:]:
                    bpy.ops.object.material_slot_add()
                    o.material_slots[-1].material = bpy.data.materials[matname]
            while len(o.material_slots) > scene.vi_leg_levels + 1:
                    bpy.ops.object.material_slot_remove()
                    
        for f, frame in enumerate(frames):
            if bm.faces.layers.float.get('res{}'.format(frame)):
                livires = bm.faces.layers.float['res{}'.format(frame)] 
                ovals = array([f[livires] for f in bm.faces])
            elif bm.verts.layers.float.get('res{}'.format(frame)):
                livires = bm.verts.layers.float['res{}'.format(frame)] 
                ovals = array([sum([vert[livires] for vert in f.verts])/len(f.verts) for f in bm.faces])
            
            if scene.vi_leg_max > scene.vi_leg_min:
                vals = ovals - scene.vi_leg_min
                vals = vals/(scene.vi_leg_max - scene.vi_leg_min)
            else:
                vals = array([scene.vi_leg_max for f in bm.faces])
                        
            nmatis = digitize(vals, bins) + 1

            if len(frames) == 1:                
                o.data.polygons.foreach_set('material_index', nmatis)
                o.data.update()

            elif len(frames) > 1:
                for fi, fc in enumerate(o.data.animation_data.action.fcurves):
                    fc.keyframe_points[f].co = frame, nmatis[fi]
        bm.free()
    scene.frame_set(scene.frame_current)
    
    
def liviresupdate(self, context):
    setscenelivivals(context.scene)
    for o in [o for o in bpy.data.objects if o.lires]:
        o.lividisplay(context.scene)  
    eupdate(self, context)

def script_update(self, context):
    if context.scene.vi_res_process == '2':
        script = bpy.data.texts[context.scene.script_file]
        exec(script.as_string())
        
def flovi_levels(self, context):
    if self.flovi_slmin > self.flovi_slmax:
       self.flovi_slmin -= 1 


classes = (VIPreferences, ViNetwork, ViLoc, ViLocSock, ViSPNode, NODE_OT_SunPath, VIEW3D_OT_SPNumDisplay, VI_PT_3D, VI_Params)


#def register():
#    for cl in classes:
#        bpy.utils.register_class(cl)
#
#    bpy.types.Scene.slicer_settings = bpy.props.PointerProperty(type=Slicer_Settings)
#
#def unregister():
#    bpy.types.Scene.slicer_settings
#    
#    for cl in classes:
#        bpy.utils.unregister_class(cl)                      

def register():
    for cl in classes:
        bpy.utils.register_class(cl)
  
    Object, Scene, Material = bpy.types.Object, bpy.types.Scene, bpy.types.Material
    bpy.types.Scene.vi_params = bpy.props.PointerProperty(type = VI_Params)
    
    


# VI-Suite object definitions
    Object.vi_type = eprop([("0", "None", "Not a VI-Suite zone"), ("1", "EnVi Zone", "Designates an EnVi Thermal zone"), 
                            ("2", "CFD Domain", "Specifies an OpenFoam BlockMesh"), ("3", "CFD Geometry", "Specifies an OpenFoam geometry"),
                            ("4", "Light Array", "Specifies a LiVi lighting array"), ("5", "Complex Fenestration", "Specifies complex fenestration for BSDF generation")], "", "Specify the type of VI-Suite zone", "0")

# LiVi object properties
    Object.livi_merr = bprop("LiVi simple mesh export", "Boolean for simple mesh export", False)
    Object.ies_name = bpy.props.StringProperty(name="", description="Name of the IES file", default="", subtype="FILE_PATH")
    Object.ies_strength = fprop("", "Strength of IES lamp", 0, 1, 1)
    Object.ies_unit = eprop([("m", "Meters", ""), ("c", "Centimeters", ""), ("f", "Feet", ""), ("i", "Inches", "")], "", "Specify the IES file measurement unit", "m")
    Object.ies_colmenu = eprop([("0", "RGB", ""), ("1", "Temperature", "")], "", "Specify the IES colour type", "0")
    Object.ies_rgb = fvprop(3, "",'IES Colour', [1.0, 1.0, 1.0], 'COLOR', 0, 1)
    Object.ies_ct = iprop("", "Colour temperature in Kelven", 0, 12000, 4700)
    (Object.licalc, Object.lires, Object.limerr, Object.manip, Object.bsdf_proxy) = [bprop("", "", False)] * 5
    Object.compcalcapply = compcalcapply    
    Object.basiccalcapply = basiccalcapply 
    Object.rtpoints = rtpoints
    Object.udidacalcapply = udidacalcapply
    Object.lividisplay = lividisplay
    Object.lhcalcapply = lhcalcapply
    Object.li_bsdf_direc = EnumProperty(items = [('+b -f', 'Backwards', 'Backwards BSDF'), ('+f -b', 'Forwards', 'Forwards BSDF'), ('+b +f', 'Bi-directional', 'Bi-directional BSDF')], name = '', description = 'BSDF direction', default = '+b -f')
    Object.li_bsdf_proxy = bprop("", "Include proxy geometry in the BSDF", False)
    Object.li_bsdf_tensor = EnumProperty(items = [(' ', 'Klems', 'Uniform Klems sample'), ('-t3', 'Symmentric', 'Symmetric Tensor BSDF'), ('-t4', 'Assymmetric', 'Asymmetric Tensor BSDF')], name = '', description = 'BSDF tensor', default = ' ')
    Object.li_bsdf_res = EnumProperty(items = [('1', '2x2', '2x2 sampling resolution'), ('2', '4x4', '4x4 sampling resolution'), ('3', '8x8', '8x8 sampling resolution'), ('4', '16x16', '16x16 sampling resolution'), ('5', '32x32', '32x32 sampling resolution'), ('6', '64x64', '64x64 sampling resolution'), ('7', '128x128', '128x128 sampling resolution')], name = '', description = 'BSDF resolution', default = '4')
    Object.li_bsdf_tsamp = IntProperty(name = '', description = 'Tensor samples', min = 1, max = 20, default = 4)
    Object.li_bsdf_ksamp = IntProperty(name = '', description = 'Klem samples', min = 1, default = 200)
    Object.li_bsdf_rcparam = sprop("", "rcontrib parameters", 1024, "")
    Object.radbsdf = radbsdf
    Object.retsv = retsv

# EnVi zone definitions
    Object.envi_type = eprop([("0", "Thermal", "Thermal Zone"), ("1", "Shading", "Shading Object"), ("2", "Chimney", "Thermal Chimney Object")], "EnVi object type", "Specify the EnVi object type", "0")
    Object.envi_oca = eprop([("0", "Default", "Use the system wide convection algorithm"), ("1", "Simple", "Use the simple convection algorithm"), ("2", "TARP", "Use the detailed convection algorithm"), ("3", "DOE-2", "Use the Trombe wall convection algorithm"), ("4", "MoWitt", "Use the adaptive convection algorithm"), ("5", "Adaptive", "Use the adaptive convection algorithm")], "", "Specify the EnVi zone outside convection algorithm", "0")
    Object.envi_ica = eprop([("0", "Default", "Use the system wide convection algorithm"), ("1", "Simple", "Use the simple convection algorithm"), ("2", "Detailed", "Use the detailed convection algorithm"), ("3", "Trombe", "Use the Trombe wall convection algorithm"), ("4", "Adaptive", "Use the adaptive convection algorithm")], "", "Specify the EnVi zone inside convection algorithm", "0")

# FloVi object definitions
    Object.flovi_fl = IntProperty(name = '', description = 'SnappyHexMesh object features levels', min = 1, max = 20, default = 4) 
    Object.flovi_slmax = IntProperty(name = '', description = 'SnappyHexMesh surface maximum levels', min = 1, max = 20, default = 4, update=flovi_levels)   
    Object.flovi_slmin = IntProperty(name = '', description = 'SnappyHexMesh surface minimum levels', min = 1, max = 20, default = 3, update=flovi_levels)     
    Object.flovi_sl = IntProperty(name = '', description = 'SnappyHexMesh surface minimum levels', min = 0, max = 20, default = 3) 
# Vi_suite material definitions
    Material.mattype = eprop([("0", "Geometry", "Geometry"), ("1", 'Light sensor', "LiVi sensing material".format(u'\u00b3')), ("2", "FloVi boundary", 'FloVi blockmesh boundary')], "", "VI-Suite material type", "0")
                                 
# LiVi material definitions                              
    Material.radmat = radmat
    Material.radmatdict = {'0': ['radcolour', 0, 'radrough', 'radspec'], '1': ['radcolour'], '2': ['radcolour', 0, 'radior'], '3': ['radcolour', 0, 'radspec', 'radrough', 0, 'radtrans',  'radtranspec'], '4': ['radcolour'], 
    '5': ['radcolmenu', 0, 'radcolour', 0, 'radct',  0, 'radintensity'], '6': ['radcolour', 0, 'radrough', 'radspec'], '7': [], '8': [], '9': []}
    Material.pport = bprop("", "Flag to signify whether the material represents a Photon Port", False)
    Material.radtex = bprop("", "Flag to signify whether the material has a texture associated with it", False)
    Material.radnorm = bprop("", "Flag to signify whether the material has a normal map associated with it", False)
    Material.ns = fprop("", "Strength of normal effect", 0, 5, 1)
    Material.nu = fvprop(3, '', 'Image up vector', [0, 0, 1], 'VELOCITY', -1, 1)
    Material.nside = fvprop(3, '', 'Image side vector', [-1, 0, 0], 'VELOCITY', -1, 1)
    radtypes = [('0', 'Plastic', 'Plastic Radiance material'), ('1', 'Glass', 'Glass Radiance material'), ('2', 'Dielectric', 'Dialectric Radiance material'),
                ('3', 'Translucent', 'Translucent Radiance material'), ('4', 'Mirror', 'Mirror Radiance material'), ('5', 'Light', 'Emission Radiance material'),
                ('6', 'Metal', 'Metal Radiance material'), ('7', 'Anti-matter', 'Antimatter Radiance material'), ('8', 'BSDF', 'BSDF Radiance material'), ('9', 'Custom', 'Custom Radiance material')]
    Material.radmatmenu = eprop(radtypes, "", "Type of Radiance material", '0')
    Material.radcolour = fvprop(3, "Material Colour",'Material Colour', [0.8, 0.8, 0.8], 'COLOR', 0, 1)
    Material.radcolmenu = eprop([("0", "RGB", "Specify colour temperature"), ("1", "Temperature", "Specify colour temperature")], "Colour type:", "Specify the colour input", "0")
    Material.radrough = fprop("Roughness", "Material roughness", 0, 1, 0.1)
    Material.radspec = fprop("Specularity", "Material specular reflection", 0, 1, 0.0)
    Material.radtrans = fprop("Transmission", "Material diffuse transmission", 0, 1, 0.1)
    Material.radtranspec  = fprop("Trans spec", "Material specular transmission", 0, 1, 0.1)
    Material.radior  = fprop("IOR", "Material index of refractionn", 0, 5, 1.5)
    Material.radct = iprop("Temperature (K)", "Colour temperature in Kelven", 0, 12000, 4700)
    Material.radintensity = fprop("Intensity", u"Material radiance (W/sr/m\u00b2)", 0, 100, 1)   
    Material.radfile = sprop("", "Radiance file material description", 1024, "")
    Material.vi_shadow = bprop("VI Shadow", "Flag to signify whether the material represents a VI Shadow sensing surface", False)
    Material.livi_sense = bprop("LiVi Sensor", "Flag to signify whether the material represents a LiVi sensing surface", False)
    Material.livi_compliance = bprop("LiVi Compliance Surface", "Flag to signify whether the material represents a LiVi compliance surface", False)
    Material.gl_roof = bprop("Glazed Roof", "Flag to signify whether the communal area has a glazed roof", False)
    hspacetype = [('0', 'Public/Staff', 'Public/Staff area'), ('1', 'Patient', 'Patient area')]
    rspacetype = [('0', "Kitchen", "Kitchen space"), ('1', "Living/Dining/Study", "Living/Dining/Study area"), ('2', "Communal", "Non-residential or communal area")]
    respacetype = [('0', "Sales", "Sales space"), ('1', "Occupied", "Occupied space")]
    lespacetype = [('0', "Healthcare", "Healthcare space"), ('1', "Other", "Other space")]
    
    Material.hspacemenu = eprop(hspacetype, "", "Type of healthcare space", '0')
    Material.brspacemenu = eprop(rspacetype, "", "Type of residential space", '0')
    Material.crspacemenu = eprop(rspacetype[:2], "", "Type of residential space", '0')
    Material.respacemenu = eprop(respacetype, "", "Type of retail space", '0')
    Material.lespacemenu = eprop(lespacetype, "", "Type of space", '0')
    Material.BSDF = bprop("", "Flag to signify a BSDF material", False)


    
# Scene parameters
    Scene.latitude = bpy.props.FloatProperty(name = "Latitude", description = "Site decimal latitude (N is positive)", min = -89.99, max = 89.99, default = 52.0)
    Scene.longitude = bpy.props.FloatProperty(name = "Longitude", description = "Site decimal longitude (E is positive)", min = -180, max = 180, default = 0.0)
    Scene.wind_type = eprop([("0", "Speed", "Wind Speed (m/s)"), ("1", "Direction", "Wind Direction (deg. from North)")], "", "Wind metric", "0")
    Scene.vipath = sprop("VI Path", "Path to files included with the VI-Suite ", 1024, addonpath)        
    Scene.suns = EnumProperty(items = [('0', 'Single', 'Single sun'), ('1', 'Monthly', 'Monthly sun for chosen time'), ('2', 'Hourly', 'Hourly sun for chosen date')], name = '', description = 'Sunpath sun type', default = '0', update=sunpath1)
    Scene.sunsstrength = bpy.props.FloatProperty(name = "", description = "Sun strength", min = 0, max = 100, default = 0.1, update=sunpath1)
    Scene.sunssize = bpy.props.FloatProperty(name = "", description = "Sun size", min = 0, max = 1, default = 0.01, update=sunpath1)
    Scene.solday = IntProperty(name = "", description = "Day of year", min = 1, max = 365, default = 1, update=sunpath1)
    Scene.solhour = bpy.props.FloatProperty(name = "", description = "Time of day", subtype='TIME', unit='TIME', min = 0, max = 24, default = 12, update=sunpath1)
    (Scene.hourdisp, Scene.spupdate, Scene.timedisp) = [bprop("", "",0)] * 3
    Scene.li_disp_panel = iprop("Display Panel", "Shows the Display Panel", -1, 2, 0)
    Scene.li_disp_count = iprop("", "", 0, 1000, 0)
    Scene.vi_disp_3d = bprop("VI 3D display", "Boolean for 3D results display",  False)
    Scene.vi_disp_3dlevel = bpy.props.FloatProperty(name = "", description = "Level of 3D result plane extrusion", min = 0, max = 500, default = 0, update = eupdate)
    Scene.ss_disp_panel = iprop("Display Panel", "Shows the Display Panel", -1, 2, 0)
    (Scene.lic_disp_panel, Scene.vi_display, Scene.sp_disp_panel, Scene.wr_disp_panel, Scene.ss_leg_display, Scene.en_disp_panel, Scene.li_compliance, Scene.vi_display_rp, Scene.vi_leg_display, 
     Scene.vi_display_sel_only, Scene.vi_display_vis_only) = [bprop("", "", False)] * 11
    Scene.vi_leg_max = bpy.props.FloatProperty(name = "", description = "Legend maximum", min = 0, max = 1000000, default = 1000, update=legupdate)
    Scene.vi_leg_min = bpy.props.FloatProperty(name = "", description = "Legend minimum", min = 0, max = 1000000, default = 0, update=legupdate)
    Scene.vi_scatter_max = bpy.props.FloatProperty(name = "", description = "Scatter maximum", min = 0, max = 1000000, default = 1000, update=legupdate)
    Scene.vi_scatter_min = bpy.props.FloatProperty(name = "", description = "Scatter minimum", min = 0, max = 1000000, default = 0, update=legupdate)
    Scene.vi_leg_scale = EnumProperty(items = [('0', 'Linear', 'Linear scale'), ('1', 'Log', 'Logarithmic scale')], name = "", description = "Legend scale", default = '0', update=legupdate)    
    Scene.vi_leg_col = EnumProperty(items = [('rainbow', 'Rainbow', 'Rainbow colour scale'), ('gray', 'Grey', 'Grey colour scale'), ('hot', 'Hot', 'Hot colour scale'),
                                             ('CMRmap', 'CMR', 'CMR colour scale'), ('jet', 'Jet', 'Jet colour scale'), ('plasma', 'Plasma', 'Plasma colour scale'), 
                                             ('hsv', 'HSV', 'HSV colour scale'), ('viridis', 'Viridis', 'Viridis colour scale')], 
                                            name = "", description = "Legend scale", default = 'rainbow', update=colupdate)
    Scene.vi_res_mod = sprop("", "Result modifier", 1024, "")
#    Scene.vi_res_py = bprop("", "Boolean for Python function modification of results",  False)
    Scene.vi_res_process = EnumProperty(items = [("0", "None", ""), ("1", "Modifier", ""), ("2", "Script", "")], name = "", description = "Specify the type of data processing", default = "0", update = script_update)
    Scene.script_file = bpy.props.StringProperty(description="Text file to show", update = script_update)
    Scene.vi_leg_unit = sprop("", "Legend unit", 1024, "")
    Scene.vi_leg_levels = IntProperty(name = "", description = "Day of year", min = 2, max = 100, default = 20, update=legupdate)
    Scene.vi_bsdfleg_max = bpy.props.FloatProperty(name = "", description = "Legend maximum", min = 0, max = 1000000, default = 100)
    Scene.vi_bsdfleg_min = bpy.props.FloatProperty(name = "", description = "Legend minimum", min = 0, max = 1000000, default = 0)
    Scene.vi_bsdfleg_scale = EnumProperty(items = [('0', 'Linear', 'Linear scale'), ('1', 'Log', 'Logarithmic scale')], name = "", description = "Legend scale", default = '0')    
    Scene.vi_gridify_rot = fprop("deg", "Rotation around face normal", 0.0, 360, 0.0)
    Scene.vi_gridify_us = fprop("m", "Up direction size", 0.01, 10, 0.6)
    Scene.vi_gridify_as = fprop("m", "Side direction size", 0.01, 10, 0.6)


    Scene.vi_display_rp_fs = iprop("", "Point result font size", 4, 24, 24)
    Scene.vi_display_rp_fc = fvprop(4, "", "Font colour", [0.0, 0.0, 0.0, 1.0], 'COLOR', 0, 1)
    Scene.vi_display_rp_sh = bprop("", "Toggle for font shadow display",  False)
    Scene.vi_display_rp_fsh = fvprop(4, "", "Font shadow", [0.0, 0.0, 0.0, 1.0], 'COLOR', 0, 1)
    Scene.vi_display_rp_off = fprop("", "Surface offset for number display", 0, 5, 0.001)
    Scene.vi_disp_trans = bpy.props.FloatProperty(name = "", description = "Sensing material transparency", min = 0, max = 1, default = 1, update = tupdate)
    Scene.vi_disp_wire = bpy.props.BoolProperty(name = "", description = "Draw wire frame", default = 0, update=wupdate)
    Scene.vi_disp_mat = bpy.props.BoolProperty(name = "", description = "Turn on/off result material emission", default = 0, update=colupdate)
    Scene.vi_disp_ems = bpy.props.FloatProperty(name = "", description = "Emissive strength", default = 1, min = 0, update=colupdate)
    Scene.li_disp_sv = EnumProperty(items = [("0", "Daylight Factor", "Display Daylight factor"),("1", "Sky view", "Display the Sky View")], name = "", description = "Compliance data type", default = "0", update = liviresupdate)
    Scene.li_disp_sda = EnumProperty(items = [("0", "sDA (%)", "Display spatial Daylight Autonomy"), ("1", "ASE (hrs)", "Display the Annual Solar Exposure")], name = "", description = "Compliance data type", default = "0", update = liviresupdate)
    Scene.li_disp_wr = EnumProperty(items = [("0", "Wind Speed", "Wind speed (m/s)"),("1", "Wind Direction", "Wind direction (deg from North)")], name = "", description = "Compliance data type", default = "0", update = liviresupdate)
 #   Scene.li_disp_lh = EnumProperty(items = [("0", "Mluxhours", "Display mega luxhours"), ("1", "Visible Irradiance", "Display visible irradiance"), ("1", "Full Irradiance", "Display full irradiance")], name = "", description = "Exposure data type", default = "0", update = liviresupdate)
#    Scene.li_projname = sprop("", "Name of the building project", 1024, '')
#    Scene.li_assorg = sprop("", "Name of the assessing organisation", 1024, '')
#    Scene.li_assind = sprop("", "Name of the assessing individual", 1024, '')
#    Scene.li_jobno = sprop("", "Project job number", 1024, '')
    Scene.li_disp_basic = EnumProperty(items = [("0", "Illuminance", "Display Illuminance values"), ("1", "Visible Irradiance", "Display Irradiance values"), ("2", "Full Irradiance", "Display Irradiance values"), ("3", "DF", "Display Daylight factor values")], name = "", description = "Basic metric selection", default = "0", update = liviresupdate)
    Scene.li_disp_da = EnumProperty(items = [("0", "DA", "Daylight Autonomy"), ("1", "sDA", "Spatial Daylight Autonomy"), ("2", "UDILow", "Spatial Daylight Autonomy"), ("3", "UDISup", "Spatial Daylight Autonomy"), 
                                             ("4", "UDIAuto", "Spatial Daylight Autonomy"), ("5", "UDIHigh", "Spatial Daylight Autonomy"), ("6", "ASE", "Annual sunlight exposure"), ("7", "Max lux", "Maximum lux level"), 
                                             ("8", "Avg Lux", "Average lux level"), ("9", "Min lux", "Minimum lux level")], name = "", description = "Result selection", default = "0", update = liviresupdate)
    Scene.li_disp_exp = EnumProperty(items = [("0", "LuxHours", "Display LuhHours values"), ("1", "Full Irradiance", "Display full spectrum radiation exposure values"), ("2", "Visible Irradiance", "Display visible spectrum radiation exposure values"),
                                              ("3", "Full Irradiance Density", "Display full spectrum radiation exposure values"), ("4", "Visible Irradiance Density", "Display visible spectrum radiation exposure values")], name = "", description = "Result selection", default = "0", update = liviresupdate)
    Scene.li_disp_irrad = EnumProperty(items = [("0", "kWh", "Display kWh values"), ("1", "kWh/m2", "Display kWh/m2 values")], name = "", description = "Result selection", default = "0", update = liviresupdate)
         
    Scene.envi_flink = bprop("", "Associate flow results with the nearest object", False)

    nodeitems_utils.register_node_categories("Vi Nodes", vinode_categories)

    
    if update_chart_node not in bpy.app.handlers.load_post:
        bpy.app.handlers.load_post.append(update_chart_node)
        
    if display_off not in bpy.app.handlers.load_post:
        bpy.app.handlers.load_post.append(display_off)
    
    if mesh_index not in bpy.app.handlers.load_post:
        bpy.app.handlers.load_post.append(mesh_index)
        
    if update_dir not in bpy.app.handlers.load_post:
        bpy.app.handlers.load_post.append(update_dir)
        
    path_update()

def unregister():
#    bpy.types.Scene.slicer_settings
    
    for cl in classes:
        bpy.utils.unregister_class(cl)
        
    nodeitems_utils.unregister_node_categories("Vi Nodes")
    nodeitems_utils.unregister_node_categories("EnVi Nodes")
    nodeitems_utils.unregister_node_categories("EnVi Mat Nodes")
        
#def unregister():
#    bpy.utils.unregister_module(__name__)
#    nodeitems_utils.unregister_node_categories("Vi Nodes")
#    nodeitems_utils.unregister_node_categories("EnVi Nodes")
#    nodeitems_utils.unregister_node_categories("EnVi Mat Nodes")
#
#    if update_chart_node in bpy.app.handlers.load_post:
#        bpy.app.handlers.load_post.remove(update_chart_node)
#
#    if display_off in bpy.app.handlers.load_post:
#        bpy.app.handlers.load_post.remove(display_off)
#
#    if mesh_index in bpy.app.handlers.load_post:
#        bpy.app.handlers.load_post.remove(mesh_index)
#        
#    if update_dir in bpy.app.handlers.load_post:
#        bpy.app.handlers.load_post.remove(update_dir)
#if __name__ == "__main__":
#    register()