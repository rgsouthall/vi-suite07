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

import bpy, datetime, mathutils, os, bmesh, shutil, sys, math, shlex, gpu, bgl
import numpy
from numpy import arange, histogram, array, int8, float16
import bpy_extras.io_utils as io_utils
from subprocess import Popen, PIPE, call
from collections import OrderedDict
from datetime import datetime as dt
from math import cos, sin, pi, ceil, tan, radians
from time import sleep
from mathutils import Euler, Vector
from gpu_extras.batch import batch_for_shader
#from multiprocessing import Pool
from .livi_export import radgexport, spfc, createoconv, createradfile
from .livi_calc  import li_calc
from .vi_display import li_display, linumdisplay, spnumdisplay, en_air, wr_legend, wr_disp, wr_scatter, wr_table, ss_disp, ss_legend, svf_disp, svf_legend, basic_legend, basic_table, basic_disp, ss_scatter, en_disp, en_pdisp, en_scatter, en_table, en_barchart, comp_table, comp_disp, leed_scatter, cbdm_disp, cbdm_scatter, envals, bsdf, bsdf_disp#, en_barchart, li3D_legend
from .envi_export import enpolymatexport, pregeo
from .envi_mat import envi_materials, envi_constructions
from .vi_func import selobj, delobj, joinobj, livisimacc, solarPosition, clearscene, clearfiles, viparams, objmode, nodecolour, cmap, wind_rose, compass, windnum, leg_min_max
from .flovi_func import fvcdwrite, fvbmwrite, fvblbmgen, fvvarwrite, fvsolwrite, fvschwrite, fvtppwrite, fvraswrite, fvshmwrite, fvmqwrite, fvsfewrite, fvobjwrite, fvdcpwrite
from .vi_func import retobjs, rettree, retpmap, progressbar, spathrange, progressfile, ret_plt
from .vi_func import chunks, xy2radial, logentry, sunpath, radpoints, sunposenvi, clearlayers, fvprogressbar, fvprogressfile
from .envi_func import processf, retenvires, envizres, envilres, recalculate_text
from .vi_chart import chart_disp

try:    
    import matplotlib
    matplotlib.use('qt5agg', warn = False, force = True)
    import matplotlib.cm as mcm
    import matplotlib.colors as mcolors
    mp = 1    
except Exception as e:
#    logentry('Matplotlib error: {}'.format(e))    
    mp = 0

if mp:
    plt = ret_plt()
    if plt:
        from .windrose import WindroseAxes

try:
    import psutil
    psu = 1
except: 
    psu = 0    

envi_mats = envi_materials()
envi_cons = envi_constructions()

rvuerrdict = {'view up parallel to view direction': "Camera cannot point directly upwards", 
              ' x11': "No X11 display server found. You may need to install XQuartz", 
              'source center': "A light source has concave faces. Use mesh - cleanup - split concave faces"}
pmerrdict = {'fatal - too many prepasses, no global photons stored\n': "Too many prepasses have occurred. Make sure light sources can see your geometry",
             'fatal - too many prepasses, no global photons stored, no caustic photons stored\n': "Too many prepasses have occurred. Turn off caustic photons and encompass the scene",
               'fatal - zero flux from light sources\n': "No light flux, make sure there is a light source and that photon port normals point inwards",
               'fatal - no light sources in distribPhotons\n': "No light sources. Photon mapping does not work with HDR skies",
               'fatal - no valid photon ports found\n': 'Make sure photon ports are valid', 
               'fatal - failed photon distribution\n': 'Do the lights see enough geometry?'}

class NODE_OT_LiGExport(bpy.types.Operator):
    bl_idname = "node.ligexport"
    bl_label = "LiVi geometry export"
    nodeid = bpy.props.StringProperty()

    def invoke(self, context, event):
        scene = context.scene
        if viparams(self, scene):
            return {'CANCELLED'}
        scene['viparams']['vidisp'] = ''
        scene['viparams']['viexpcontext'] = 'LiVi Geometry'
        objmode()
        clearfiles(scene['liparams']['objfilebase'])
        clearfiles(scene['liparams']['lightfilebase'])
        node = bpy.data.node_groups[self.nodeid.split('@')[1]].nodes[self.nodeid.split('@')[0]]
        node.preexport(scene)
        radgexport(self, node)
        node.postexport(scene)
        return {'FINISHED'}
        
class OBJECT_GenBSDF(bpy.types.Operator):
    bl_idname = "object.gen_bsdf"
    bl_label = "Gen BSDF"
    bl_description = "Generate a BSDF for the current selected object"
    bl_register = True
    bl_undo = False
    
    def modal(self, context, event):
        if self.bsdfrun.poll() is None:
            if self.pfile.check(0) == 'CANCELLED':                   
                self.bsdfrun.kill()   
                
                if psu:
                    for proc in psutil.process_iter():
                        if 'rcontrib' in proc.name():
                            proc.kill()  
                else:
                    self.report({'ERROR'}, 'psutil not found. Kill rcontrib processes manually')
                        
                return {'CANCELLED'}
            else:
                return{'PASS_THROUGH'}
        else:
            if self.kivyrun.poll() is None:
                self.kivyrun.kill() 
            
            if self.o.li_bsdf_proxy:
                self.pkgbsdfrun = Popen(shlex.split("pkgBSDF -s {}".format(os.path.join(context.scene['viparams']['newdir'], 'bsdfs', '{}.xml'.format(self.mat.name)))), stdin = PIPE, stdout = PIPE)
                self.mat['radentry'] = ''.join([line.decode() for line in self.pkgbsdfrun.stdout])
                print(self.mat['radentry'])
            with open(os.path.join(context.scene['viparams']['newdir'], 'bsdfs', '{}.xml'.format(self.mat.name)), 'r') as bsdffile:
                
                self.mat['bsdf']['xml'] = bsdffile.read()

            context.scene['viparams']['vidisp'] = 'bsdf'
            self.mat['bsdf']['type'] = self.o.li_bsdf_tensor
            return {'FINISHED'}
    
    def execute(self, context):
        scene = context.scene
        self.o = context.active_object
        scene.objects.active = None
        
        if viparams(self, scene):
            return {'CANCELLED'}
        bsdfmats = [mat for mat in self.o.data.materials if mat.radmatmenu == '8']
        
        if bsdfmats:
            self.mat = bsdfmats[0]
            self.mat['bsdf'] = {} 
        else:
            self.report({'ERROR'}, '{} does not have a BSDF material attached'.format(self.o.name))
        
        tm = self.o.to_mesh(scene = scene, apply_modifiers = True, settings = 'PREVIEW')
        bm = bmesh.new()    
        bm.from_mesh(tm) 
        bpy.data.meshes.remove(tm)
        bm.transform(self.o.matrix_world)
        bm.normal_update()
        bsdffaces = [face for face in bm.faces if self.o.data.materials[face.material_index].radmatmenu == '8']    
        
        if bsdffaces:
            fvec = bsdffaces[0].normal
            self.mat['bsdf']['normal'] = '{0[0]:.4f} {0[1]:.4f} {0[2]:.4f}'.format(fvec)
        else:
            self.report({'ERROR'}, '{} does not have a BSDF material associated with any faces'.format(self.o.name))
            return
        
        self.pfile = progressfile(scene['viparams']['newdir'], datetime.datetime.now(), 100)
        self.kivyrun = progressbar(os.path.join(scene['viparams']['newdir'], 'viprogress'), 'BSDF')
        zvec, xvec = mathutils.Vector((0, 0, 1)), mathutils.Vector((1, 0, 0))
        svec = mathutils.Vector.cross(fvec, zvec)
        bm.faces.ensure_lookup_table()
        bsdfrotz = mathutils.Matrix.Rotation(mathutils.Vector.angle(fvec, zvec), 4, svec)
        bm.transform(bsdfrotz)
        bsdfrotx = mathutils.Matrix.Rotation(math.pi + mathutils.Vector.angle_signed(mathutils.Vector(xvec[:2]), mathutils.Vector(svec[:2])), 4, zvec)#mathutils.Vector.cross(svec, xvec))
        bm.transform(bsdfrotx)
        vposis = list(zip(*[v.co[:] for v in bm.verts]))
        (maxx, maxy, maxz) = [max(p) for p in vposis]
        (minx, miny, minz) = [min(p) for p in vposis]
        bsdftrans = mathutils.Matrix.Translation(mathutils.Vector((-(maxx + minx)/2, -(maxy + miny)/2, -maxz)))
        bm.transform(bsdftrans)
        mradfile = ''.join([m.radmat(scene) for m in self.o.data.materials if m.radmatmenu != '8'])                  
        gradfile = radpoints(self.o, [face for face in bm.faces if self.o.data.materials and face.material_index < len(self.o.data.materials) and self.o.data.materials[face.material_index].radmatmenu != '8'], 0)
        bm.free()  
        bsdfsamp = self.o.li_bsdf_ksamp if self.o.li_bsdf_tensor == ' ' else 2**(int(self.o.li_bsdf_res) * 2) * int(self.o.li_bsdf_tsamp) 
        gbcmd = "genBSDF +geom meter -r '{}' {} {} -c {} {} -n {}".format(self.o.li_bsdf_rcparam,  self.o.li_bsdf_tensor, (self.o.li_bsdf_res, ' ')[self.o.li_bsdf_tensor == ' '], bsdfsamp, self.o.li_bsdf_direc, scene['viparams']['nproc'])

        with open(os.path.join(scene['viparams']['newdir'], 'bsdfs', '{}_mg'.format(self.mat.name)), 'w') as mgfile:
            mgfile.write(mradfile+gradfile)

        with open(os.path.join(scene['viparams']['newdir'], 'bsdfs', '{}_mg'.format(self.mat.name)), 'r') as mgfile: 
            with open(os.path.join(scene['viparams']['newdir'], 'bsdfs', '{}.xml'.format(self.mat.name)), 'w') as bsdffile:
#            self.bsdfrun = Popen(shlex.split(gbcmd), stdin = mgfile, stdout = PIPE)
                self.bsdfrun = Popen(shlex.split(gbcmd), stdin = mgfile, stdout = bsdffile)

        wm = context.window_manager
        self._timer = wm.event_timer_add(1, context.window)
        wm.modal_handler_add(self)        
        return {'RUNNING_MODAL'}
        
class MATERIAL_LoadBSDF(bpy.types.Operator, io_utils.ImportHelper):
    bl_idname = "material.load_bsdf"
    bl_label = "Select BSDF file"
    filename_ext = ".XML;.xml;"
    filter_glob = bpy.props.StringProperty(default="*.XML;*.xml;", options={'HIDDEN'})
    filepath = bpy.props.StringProperty(subtype='FILE_PATH', options={'HIDDEN', 'SKIP_SAVE'})
    
    def draw(self,context):
        layout = self.layout
        row = layout.row()
        row.label(text="Import BSDF XML file with the file browser", icon='WORLD_DATA')
        row = layout.row()

    def execute(self, context):
        context.material['bsdf'] = {}
        if " " in self.filepath:
            self.report({'ERROR'}, "There is a space either in the filename or its directory location. Remove this space and retry opening the file.")
            return {'CANCELLED'}
        else:
            with open(self.filepath, 'r') as bsdffile:
                context.material['bsdf']['xml'] = bsdffile.read()
            return {'FINISHED'}

    def invoke(self,context,event):
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}
    
class MATERIAL_DelBSDF(bpy.types.Operator):
    bl_idname = "material.del_bsdf"
    bl_label = "Del BSDF"
    bl_description = "Delete a BSDF for the current selected object"
    bl_register = True
    bl_undo = False
    
    def execute(self, context):
        del context.material['bsdf']
        return {'FINISHED'}
        
class MATERIAL_SaveBSDF(bpy.types.Operator):
    bl_idname = "material.save_bsdf"
    bl_label = "Save BSDF"
    bl_description = "Save a BSDF for the current selected object"
    bl_register = True

    filename_ext = ".XML;.xml;"
    filter_glob = bpy.props.StringProperty(default="*.XML;*.xml;", options={'HIDDEN'})
    filepath = bpy.props.StringProperty(subtype='FILE_PATH', options={'HIDDEN', 'SKIP_SAVE'})
    
    def draw(self,context):
        layout = self.layout
        row = layout.row()
        row.label(text="Save BSDF XML file with the file browser", icon='WORLD_DATA')

    def execute(self, context):
        with open(self.filepath, 'w') as bsdfsave:
            bsdfsave.write(context.material['bsdf']['xml'])
        return {'FINISHED'}

    def invoke(self,context,event):
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}
        
class VIEW3D_OT_BSDF_Disp(bpy.types.Operator):
    bl_idname = "view3d.bsdf_display"
    bl_label = "BSDF display"
    bl_description = "Display BSDF"
    bl_register = True
    bl_undo = False
        
    def modal(self, context, event):
        if event.type != 'INBETWEEN_MOUSEMOVE' and context.region and context.area.type == 'VIEW_3D' and context.region.type == 'WINDOW':  
            if context.scene['viparams']['vidisp'] != 'bsdf_panel':
                self.remove(context)
                context.scene['viparams']['vidisp'] = self.olddisp
                return {'CANCELLED'}

            mx, my = event.mouse_region_x, event.mouse_region_y
            if self.bsdf.spos[0] < mx < self.bsdf.epos[0] and self.bsdf.spos[1] < my < self.bsdf.epos[1]:
                self.bsdf.hl = (0, 1, 1, 1)  
                
                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.bsdfpress = 1
                        self.bsdfmove = 0
                        return {'RUNNING_MODAL'}
                    elif event.value == 'RELEASE':
                        if not self.bsdfmove:
                            self.bsdf.expand = 0 if self.bsdf.expand else 1
                        self.bsdfpress = 0
                        self.bsdfmove = 0
                        context.area.tag_redraw()
                        return {'RUNNING_MODAL'}
                        
                elif event.type == 'ESC':
                    self.remove(context)
                    context.scene['viparams']['vidisp'] = self.olddisp
                    return {'CANCELLED'}                   
                elif self.bsdfpress and event.type == 'MOUSEMOVE':
                     self.bsdfmove = 1
                     self.bsdfpress = 0
                            
            elif abs(self.bsdf.lepos[0] - mx) < 10 and abs(self.bsdf.lspos[1] - my) < 10:
                self.bsdf.hl = (0, 1, 1, 1) 
                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.bsdf.resize = 1
                    if self.bsdf.resize and event.value == 'RELEASE':
                        self.bsdf.resize = 0
                    return {'RUNNING_MODAL'}  
            
            elif all((self.bsdf.expand, self.bsdf.lspos[0] + 0.45 * self.bsdf.xdiff < mx < self.bsdf.lspos[0] + 0.8 * self.bsdf.xdiff, self.bsdf.lspos[1] + 0.06 * self.bsdf.ydiff < my < self.bsdf.lepos[1] - 5)):
                if event.type == 'LEFTMOUSE' and event.value == 'PRESS':
                    self.bsdf.plt.show()
            
            else:
                for butrange in self.bsdf.buttons:
                    if self.bsdf.buttons[butrange][0] - 0.0075 * self.bsdf.xdiff < mx < self.bsdf.buttons[butrange][0] + 0.0075 * self.bsdf.xdiff and self.bsdf.buttons[butrange][1] - 0.01 * self.bsdf.ydiff < my < self.bsdf.buttons[butrange][1] + 0.01 * self.bsdf.ydiff:
                        if event.type == 'LEFTMOUSE' and event.value == 'PRESS' and self.bsdf.expand:
                            if butrange in ('Front', 'Back'):
                                self.bsdf.dir_select = butrange
                            elif butrange in ('Visible', 'Solar', 'Discrete'):
                                self.bsdf.rad_select = butrange
                            elif butrange in ('Transmission', 'Reflection'):
                                self.bsdf.type_select = butrange
                            self.bsdf.plot(context)

                self.bsdf.hl = (1, 1, 1, 1)
                                
            if event.type == 'MOUSEMOVE':                
                if self.bsdfmove:
                    self.bsdf.pos = [mx, my]
                    context.area.tag_redraw()
                    return {'RUNNING_MODAL'}
                if self.bsdf.resize:
                    self.bsdf.lepos[0], self.bsdf.lspos[1] = mx, my
            
            if self.bsdf.expand and self.bsdf.lspos[0] < mx < self.bsdf.lepos[0] and self.bsdf.lspos[1] < my < self.bsdf.lepos[1]:
                theta, phi = xy2radial(self.bsdf.centre, (mx, my), self.bsdf.pw, self.bsdf.ph)
                phi = math.atan2(-my + self.bsdf.centre[1], mx - self.bsdf.centre[0]) + math.pi

                if theta < self.bsdf.radii[-1]:
                    for ri, r in enumerate(self.bsdf.radii):
                        if theta < r:
                            break

                    upperangles = [p * 2 * math.pi/self.bsdf.phis[ri] + math.pi/self.bsdf.phis[ri]  for p in range(int(self.bsdf.phis[ri]))]
                    uai = 0

                    if ri > 0:
                        for uai, ua in enumerate(upperangles): 
                            if phi > upperangles[-1]:
                                uai = 0
                                break
                            if phi < ua:
                                break

                    self.bsdf.patch_hl = sum(self.bsdf.phis[0:ri]) + uai
                    if event.type in ('LEFTMOUSE', 'RIGHTMOUSE')  and event.value == 'PRESS':                        
                        self.bsdf.num_disp = 1 if event.type == 'RIGHTMOUSE' else 0    
                        self.bsdf.patch_select = sum(self.bsdf.phis[0:ri]) + uai
                        self.bsdf.plot(context)
                        context.area.tag_redraw()
                        return {'RUNNING_MODAL'}
                        
                else:
                    self.bsdf.patch_hl = None
                    
            if self.bsdf.expand and any((self.bsdf.leg_max != context.scene.vi_bsdfleg_max, self.bsdf.leg_min != context.scene.vi_bsdfleg_min, self.bsdf.col != context.scene.vi_leg_col, self.bsdf.scale_select != context.scene.vi_bsdfleg_scale)):
                self.bsdf.col = context.scene.vi_leg_col
                self.bsdf.leg_max = context.scene.vi_bsdfleg_max
                self.bsdf.leg_min = context.scene.vi_bsdfleg_min
                self.bsdf.scale_select = context.scene.vi_bsdfleg_scale
                self.bsdf.plot(context)
            
            context.area.tag_redraw()
        
        return {'PASS_THROUGH'}
                
    def invoke(self, context, event):
        cao = context.active_object
        if cao and cao.active_material.get('bsdf') and cao.active_material['bsdf']['xml'] and cao.active_material['bsdf']['type'] == ' ':
            width, height = context.region.width, context.region.height
            self.bsdf = bsdf([160, height - 40], width, height, 'bsdf.png', 750, 400)
            self.bsdf.update(context)
            self.bsdfpress, self.bsdfmove, self.bsdfresize = 0, 0, 0
            self._handle_bsdf_disp = bpy.types.SpaceView3D.draw_handler_add(bsdf_disp, (self, context), 'WINDOW', 'POST_PIXEL')
            self.olddisp = context.scene['viparams']['vidisp']
            context.window_manager.modal_handler_add(self)
            context.scene['viparams']['vidisp'] = 'bsdf_panel'
            context.area.tag_redraw()            
            return {'RUNNING_MODAL'}
        else:
            self.report({'ERROR'},"Selected material contains no BSDF information or contains the wrong BSDF type (only Klems is supported)")
            return {'CANCELLED'}
            
    def remove(self, context):
        self.bsdf.plt.close()
        bpy.types.SpaceView3D.draw_handler_remove(self._handle_bsdf_disp, 'WINDOW')
        context.scene['viparams']['vidisp'] = 'bsdf'
        bpy.data.images.remove(self.bsdf.gimage)
        context.area.tag_redraw()
        
class NODE_OT_FileSelect(bpy.types.Operator, io_utils.ImportHelper):
    bl_idname = "node.fileselect"
    bl_label = "Select file"
    filename = ""
    bl_register = True
    bl_undo = True

    def draw(self,context):
        layout = self.layout
        row = layout.row()
        row.label(text="Import {} file with the file browser".format(self.filename), icon='WORLD_DATA')
        row = layout.row()

    def execute(self, context):
        node = bpy.data.node_groups[self.nodeid.split('@')[1]].nodes[self.nodeid.split('@')[0]]
        if self.filepath.split(".")[-1] in self.fextlist:
            if self.nodeprop == 'epwname':
                node.epwname = self.filepath
            elif self.nodeprop == 'hdrname':
                node.hdrname = self.filepath
            elif self.nodeprop == 'skyname':
                node.skyname = self.filepath
            elif self.nodeprop == 'mtxname':
                node.mtxname = self.filepath
            elif self.nodeprop == 'resfilename':
                node.resfilename = self.filepath
            elif self.nodeprop == 'idffilename':
                node.idffilename = self.filepath
        if " " in self.filepath:
            self.report({'ERROR'}, "There is a space either in the filename or its directory location. Remove this space and retry opening the file.")
        return {'FINISHED'}

    def invoke(self,context,event):
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}

class NODE_OT_HdrSelect(NODE_OT_FileSelect):
    bl_idname = "node.hdrselect"
    bl_label = "Select HDR/VEC file"
    bl_description = "Select the HDR sky image or vector file"
    filename_ext = ".HDR;.hdr;"
    filter_glob = bpy.props.StringProperty(default="*.HDR;*.hdr;", options={'HIDDEN'})
    nodeprop = 'hdrname'
    filepath = bpy.props.StringProperty(subtype='FILE_PATH', options={'HIDDEN', 'SKIP_SAVE'})
    fextlist = ("HDR", "hdr")
    nodeid = bpy.props.StringProperty()

class NODE_OT_SkySelect(NODE_OT_FileSelect):
    bl_idname = "node.skyselect"
    bl_label = "Select RAD file"
    bl_description = "Select the Radiance sky file"
    filename_ext = ".rad;.RAD;"
    filter_glob = bpy.props.StringProperty(default="*.RAD;*.rad;", options={'HIDDEN'})
    nodeprop = 'skyname'
    filepath = bpy.props.StringProperty(subtype='FILE_PATH', options={'HIDDEN', 'SKIP_SAVE'})
    fextlist = ("RAD", "rad")
    nodeid = bpy.props.StringProperty()

class NODE_OT_MtxSelect(NODE_OT_FileSelect):
    bl_idname = "node.mtxselect"
    bl_label = "Select MTX file"
    bl_description = "Select the matrix file"
    filename_ext = ".MTX;.mtx;"
    filter_glob = bpy.props.StringProperty(default="*.MTX;*.mtx;", options={'HIDDEN'})
    nodeprop = 'mtxname'
    filepath = bpy.props.StringProperty(subtype='FILE_PATH', options={'HIDDEN', 'SKIP_SAVE'})
    fextlist = ("MTX", "mtx")
    nodeid = bpy.props.StringProperty()

#class NODE_OT_EpwSelect(bpy.types.Operator, io_utils.ImportHelper):
#    bl_idname = "node.epwselect"
#    bl_label = "Select EPW file"
#    bl_description = "Select the EnergyPlus weather file"
#    filename_ext = ".HDR;.hdr;.epw;.EPW;"
#    filter_glob = bpy.props.StringProperty(default="*.HDR;*.hdr;*.epw;*.EPW;", options={'HIDDEN'})
#    nodeprop = 'epwname'
#    filepath = bpy.props.StringProperty(subtype='FILE_PATH', options={'HIDDEN', 'SKIP_SAVE'})
#    fextlist = ("epw", "EPW", "HDR", "hdr")
#    nodeid = bpy.props.StringProperty()

class NODE_OT_LiExport(bpy.types.Operator, io_utils.ExportHelper):
    bl_idname = "node.liexport"
    bl_label = "LiVi context export"
    bl_description = "Export the scene to the Radiance file format"
    bl_register = True
    bl_undo = False
    nodeid = bpy.props.StringProperty()

    def invoke(self, context, event):
#        print(dir(context.node))
        scene = context.scene
        if viparams(self, scene):
            return {'CANCELLED'}
#        node = bpy.data.node_groups[self.nodeid.split('@')[1]].nodes[self.nodeid.split('@')[0]]
        node = context.node
        scene['viparams']['vidisp'] = ''
        scene['viparams']['viexpcontext'] = 'LiVi {}'.format(node.contextmenu)
        scene['viparams']['connode'] = self.nodeid
        objmode()
        node.preexport()
        if not node.export(scene, self):
            node.postexport()
        return {'FINISHED'}

class NODE_OT_RadPreview(bpy.types.Operator, io_utils.ExportHelper):
    bl_idname = "node.radpreview"
    bl_label = "LiVi preview"
    bl_description = "Prevew the scene with Radiance"
    bl_register = True
    bl_undo = False
    nodeid = bpy.props.StringProperty()
    
    def modal(self, context, event):
        if event.type == 'TIMER':
            if self.rvurun.poll() is not None: # If finished
                for line in self.rvurun.stderr:
                    logentry(line)
                    for rvuerr in rvuerrdict:
                        if rvuerr in line.decode():
                            self.report({'ERROR'}, rvuerrdict[rvuerr])
                            return {'CANCELLED'}

                self.simnode.run = 0
                return {'FINISHED'}
            else:           
                return {'PASS_THROUGH'}
        else:
            return {'PASS_THROUGH'}
        

    def invoke(self, context, event):
        scene = context.scene
        
        if viparams(self, scene):
            return {'CANCELLED'}
        
        objmode()
        self.simnode, frame = bpy.data.node_groups[self.nodeid.split('@')[1]].nodes[self.nodeid.split('@')[0]], scene.frame_current
        self.simnode.presim()
        scene['liparams']['fs'] = min([c['fs'] for c in (self.simnode['goptions'], self.simnode['coptions'])])
        scene['liparams']['fe'] = max([c['fe'] for c in (self.simnode['goptions'], self.simnode['coptions'])])

        if frame not in range(scene['liparams']['fs'], scene['liparams']['fe'] + 1):
            self.report({'ERROR'}, "Current frame is not within the exported frame range")
            return {'CANCELLED'}
            
        cam = scene.camera
        
        if cam:
            curres = 0.1
            createradfile(scene, frame, self, self.simnode)
            createoconv(scene, frame, self, self.simnode)
            cang = '180 -vth ' if self.simnode['coptions']['Context'] == 'Basic' and self.simnode['coptions']['Type'] == '1' else cam.data.angle_x*180/pi
            vv = 180 if self.simnode['coptions']['Context'] == 'Basic' and self.simnode['coptions']['Type'] == '1' else cam.data.angle_y*180/pi
            vd = (0.001, 0, -1*cam.matrix_world[2][2]) if (round(-1*cam.matrix_world[0][2], 3), round(-1*cam.matrix_world[1][2], 3)) == (0.0, 0.0) else [-1*cam.matrix_world[i][2] for i in range(3)]
            
            if self.simnode.pmap:
                self.pfile = progressfile(scene['viparams']['newdir'], datetime.datetime.now(), 100)
                self.kivyrun = progressbar(os.path.join(scene['viparams']['newdir'], 'viprogress'), 'Photon Map')
                amentry, pportentry, cpentry, cpfileentry = retpmap(self.simnode, frame, scene)
                open('{}.pmapmon'.format(scene['viparams']['filebase']), 'w')
                pmcmd = 'mkpmap -t 20 -e {1}.pmapmon -fo+ -bv+ -apD 0.001 {0} -apg {1}-{2}.gpm {3} {4} {5} {1}-{2}.oct'.format(pportentry, scene['viparams']['filebase'], frame, self.simnode.pmapgno, cpentry, amentry)
                logentry('Photon map command: {}'.format(pmcmd))
                pmrun = Popen(pmcmd.split(), stderr = PIPE, stdout = PIPE)

                while pmrun.poll() is None:   
                    sleep(10)
                    with open('{}.pmapmon'.format(scene['viparams']['filebase']), 'r') as vip:
                        for line in vip.readlines()[::-1]:
                            if '%' in line:
                                curres = float(line.split()[6][:-2])
                                break
                                
                    if self.pfile.check(curres) == 'CANCELLED': 
                        pmrun.kill()                                   
                        return {'CANCELLED'}
                
                if self.kivyrun.poll() is None:
                    self.kivyrun.kill()
                        
                with open('{}.pmapmon'.format(scene['viparams']['filebase']), 'r') as pmapfile:
                    for line in pmapfile.readlines():
                        if line in pmerrdict:
                            logentry(line)
                            self.report({'ERROR'}, pmerrdict[line])
                            return {'CANCELLED'}
                                        
                rvucmd = "rvu -w {11} -ap {8} 50 {9} -n {0} -vv {1:.3f} -vh {2} -vd {3[0]:.3f} {3[1]:.3f} {3[2]:.3f} -vp {4[0]:.3f} {4[1]:.3f} {4[2]:.3f} -vu {10[0]:.3f} {10[1]:.3f} {10[2]:.3f} {5} {6}-{7}.oct".format(scene['viparams']['wnproc'], vv, cang, vd, cam.location, self.simnode['radparams'], scene['viparams']['filebase'], scene.frame_current, '{}-{}.gpm'.format(scene['viparams']['filebase'], frame), cpfileentry, cam.matrix_world.to_quaternion() * mathutils.Vector((0, 1, 0)), ('', '-i')[self.simnode.illu])
                
            else:
                rvucmd = "rvu -w {9} -n {0} -vv {1} -vh {2} -vd {3[0]:.3f} {3[1]:.3f} {3[2]:.3f} -vp {4[0]:.3f} {4[1]:.3f} {4[2]:.3f} -vu {8[0]:.3f} {8[1]:.3f} {8[2]:.3f} {5} {6}-{7}.oct".format(scene['viparams']['wnproc'], vv, cang, vd, cam.location, self.simnode['radparams'], scene['viparams']['filebase'], scene.frame_current, cam.matrix_world.to_quaternion() * mathutils.Vector((0, 1, 0)), ('', '-i')[self.simnode.illu])

            logentry('Rvu command: {}'.format(rvucmd))
            self.rvurun = Popen(rvucmd.split(), stdout = PIPE, stderr = PIPE)
            self.simnode.run = 1
            wm = context.window_manager
            self._timer = wm.event_timer_add(5, context.window)
            wm.modal_handler_add(self)
            return {'RUNNING_MODAL'}

        else:
            self.report({'ERROR'}, "There is no camera in the scene. Radiance preview will not work")
            return {'CANCELLED'}

class NODE_OT_RadImage(bpy.types.Operator):
    bl_idname = "node.radimage"
    bl_label = "LiVi Image"
    bl_description = "Generate an image with Rpict"
    bl_register = True
    bl_undo = False
    nodeid = bpy.props.StringProperty()

    def modal(self, context, event):
        if self.pfile.check(self.percent) == 'CANCELLED':                                    
            return {self.terminate()}
        
        while sum([pm.poll() is None for pm in self.pmruns]) < self.processors and self.p < self.frames:
            if self.pmaps[self.p]:
                self.pmruns.append(Popen(self.pmcmds[self.p].split(), stderr = PIPE))
            self.p += 1

        if all([pm.poll() is not None for pm in self.pmruns]) and sum(self.pmaps) == self.p and not self.rpruns:  
            self.pmfin = 1
            
            if len(self.pmruns):
                self.kivyrun.kill()

        if self.pmfin:   
            if len(self.rpruns) == 0:
                self.pfile = progressfile(self.folder, datetime.datetime.now(), 100)
                if len(self.pmruns):
                    self.kivyrun = progressbar(os.path.join(self.folder, 'viprogress'), 'Radiance Image')
            
            if self.mp:
                while self.xindex < self.processes and sum([rp.poll() is None for rp in self.rpruns]) < self.processors and self.frame <= self.fe:
                    echo = Popen(['echo', '{}'.format(self.xindex), '0'], stdout = PIPE)
                    self.rpruns.append(Popen(self.rpiececmds[self.frame - self.fs].split(), stdin = echo.stdout, stderr = PIPE))
                    if self.xindex == 0:
                        if os.path.isfile("{}-{}.hdr".format(os.path.join(self.folder, 'images', self.basename), self.frame)):
                            os.remove("{}-{}.hdr".format(os.path.join(self.folder, 'images', self.basename), self.frame))
                        logentry('rpiece command: {}'.format(self.rpiececmds[self.frame - self.fs]))
                    self.xindex += 1
                
                if self.xindex == self.processes:   
                    self.frame += 1                
                    self.xindex = 0
                    self.percent = (self.frame - self.fs) * 100
                    self.images.append(os.path.join(self.folder, 'images', '{}-{}.hdr'.format(self.basename, self.frameold)))
                    self.frameold = self.frame
            else:
                while sum([rp.poll() is None for rp in self.rpruns]) == 0 and len(self.rpruns) < self.frames:
                    with open("{}-{}.hdr".format(os.path.join(self.folder, 'images', self.basename), self.frame), 'w') as imfile:
                        self.rpruns.append(Popen(self.rpictcmds[self.frame - self.fs].split(), stdout=imfile, stderr = PIPE))
                if [rp.poll() for rp in self.rpruns][self.frame - self.fs] is not None:
                    self.images.append(os.path.join(self.folder, 'images', '{}-{}.hdr'.format(self.basename, self.frame)))
                    self.frame += 1
                                        
        if event.type == 'TIMER':            
            f = self.frame if self.frame <= self.fe else self.fe
            if self.pmfin and not self.rpruns:
                for line in self.pmruns[0].stderr:
                    logentry('Photon mapping error: {}'.format(line.decode()))
                    
                    for pmerr in pmerrdict:
                        if pmerr in line.decode():
                            self.report({'ERROR'}, pmerrdict[pmerr])
                            self.kivyrun.kill()
                            self.simnode.run = 0
                            return {'CANCELLED'}
                                
            elif os.path.isfile(self.pmfile) and not self.rpruns:
                if sum(self.pmaps):
                    self.percent = 0
                    for vip in [open('{}-{}'.format(self.pmfile, frame), 'r') for frame in range(self.fs, self.fe + 1)]:
                        for line in vip.readlines()[::-1]:
                            if '% after' in line:
                                self.percent += [float(ls[:-2]) for ls in line.split() if '%' in ls][0]/sum(self.pmaps)
                                break
                            elif line in pmerrdict:
                                logentry(line)
                                self.report({'ERROR'}, pmerrdict[line])
                                return {'CANCELLED'}

            if self.pmfin and self.rpruns and all([rp.poll() is not None for rp in self.rpruns]):
                for line in self.rpruns[0].stderr:
                    logentry('Rpict error: {}'.format(line.decode()))
                    
                    for rvuerr in rvuerrdict:
                        if rvuerr in line.decode():
                            self.report({'ERROR'}, rvuerrdict[rvuerr])
                            self.kivyrun.kill()
                            self.simnode.run = 0
                            return {'CANCELLED'}
                        
                self.imupdate(f)
                return {self.terminate()}

            elif self.pmfin and self.mp:
                if self.percent != 100 * sum([r.poll() is not None for r in self.rpruns])/(self.processes * self.frames):
                    self.percent = 100 * sum([r.poll() is not None for r in self.rpruns])/(self.processes * self.frames)
                    self.imupdate(self.fs + int(sum([rp.poll() is not None for rp in self.rpruns])/self.processes))
            
            elif self.pmfin and os.path.isfile(self.rpictfile):
                lines = [line for line in open(self.rpictfile, 'r') if '% after' in line][::-1]                
                if lines:
                    for lineentry in lines[0].split():
                        if '%' in lineentry and self.percent != (float(lineentry.strip('%')) + (f - self.fs) * 100)/self.frames:
                            self.percent = (float(lineentry.strip('%')) + (f - self.fs) * 100)/self.frames
                            self.imupdate(f)
                
            return {'PASS_THROUGH'}
        else:
            return {'PASS_THROUGH'}
    
    def imupdate(self, f):
        if 'liviimage' not in bpy.data.images:
            im = bpy.data.images.load("{}-{}.hdr".format(os.path.join(self.folder, 'images', self.basename), f))
            im.name = 'liviimage'
        bpy.data.images['liviimage'].filepath = "{}-{}.hdr".format(os.path.join(self.folder, 'images', self.basename), f)
        bpy.data.images['liviimage'].reload()

        for area in bpy.context.screen.areas:
            if area.type =='IMAGE_EDITOR':
                area.tag_redraw()
        
    def terminate(self):
        nodecolour(self.simnode, 0)
        self.kivyrun.kill() 
        self.simnode.run = 0
        for pm in self.pmruns:
            if pm.poll() is None:
                pm.kill()
        for rp in self.rpruns:
            if rp.poll() is None:                         
                rp.kill()
       
        self.simnode.postsim(self.images)
        if os.path.isfile(self.rpictfile):
            os.remove(self.rpictfile)
        return 'FINISHED'

    def execute(self, context):        
        scene = context.scene
        self.xindex, self.p = 0, 0
        self.cam = scene.camera
        simnode = bpy.data.node_groups[self.nodeid.split('@')[1]].nodes[self.nodeid.split('@')[0]]
        self.fs, self.fe = simnode.retframes()
        self.simnode = simnode
        
        if simnode.camera and bpy.data.cameras.get(simnode.camera):
            self.percent = 0
            self.reslists, self.images = [], []
            self.res = []
            self.rpictfile = os.path.join(scene['viparams']['newdir'], 'rpictprogress')
            self.pmfile = os.path.join(scene['viparams']['newdir'], 'pmprogress')
            simnode.presim()
            scene['liparams']['fs'], scene['liparams']['fe'] =  simnode.retframes()
            self.frames = self.fe - self.fs + 1
            self.frame = self.fs
            self.frameold = self.frame
            self.rpruns, self.pmruns = [], []   
            self.processors = simnode['Processors']
            self.processes = simnode.processes
            self.radparams = simnode['radparams']
            self.viewparams = simnode['viewparams']
            self.pmparams = simnode['pmparams']
            self.pmaps = simnode['pmaps']
            self.pmapgnos = simnode['pmapgnos']
            self.pmapcnos = simnode['pmapcnos']
            self.folder = scene['viparams']['newdir']
            self.fb = scene['viparams']['filebase']
            self.mp = 0 if sys.platform == 'win32' else simnode.mp 
            self.basename = simnode['basename']
            
            for frame in range(self.fs, self.fe + 1):
                createradfile(scene, frame, self, simnode)
                createoconv(scene, frame, self, simnode)
                with open('{}-{}'.format(self.pmfile, frame), 'w'):
                    pass
                
            scene.frame_set(scene['liparams']['fs'])
            self.pmcmds = ['mkpmap -t 10 -e {6} -bv+ +fo -apD 0.001 {0} -apg {1}-{2}.gpm {3} {4} {5} {1}-{2}.oct'.format(self.pmparams[str(frame)]['pportentry'], scene['viparams']['filebase'], frame, self.pmapgnos[str(frame)], self.pmparams[str(frame)]['cpentry'], self.pmparams[str(frame)]['amentry'], '{}-{}'.format(self.pmfile, frame)) for frame in range(self.fs, self.fe + 1)]                   
            self.rppmcmds = [('', ' -ap {} {}'.format('{}-{}.gpm'.format(scene['viparams']['filebase'], frame), self.pmparams[str(frame)]['cpfileentry']))[self.pmaps[frame - self.fs]] for frame in range(self.fs, self.fe + 1)]
            self.rpictcmds = ["rpict -t 10 -e {} ".format(self.rpictfile) + ' '.join(['{0[0]} {0[1]}'.format(i) for i in self.viewparams[str(frame)].items()]) + self.rppmcmds[frame - self.fs] + self.radparams + "{0}-{1}.oct".format(scene['viparams']['filebase'], frame, os.path.join(scene['viparams']['newdir'], 'images', self.basename)) for frame in range(self.fs, self.fe + 1)]
            self.rpiececmds = ["rpiece -t 10 -e {} ".format(self.rpictfile) + ' '.join(['{0[0]} {0[1]}'.format(i) for i in self.viewparams[str(frame)].items()]) + self.rppmcmds[frame - self.fs] + self.radparams + "-o {2}-{1}.hdr {0}-{1}.oct".format(scene['viparams']['filebase'], frame, os.path.join(scene['viparams']['newdir'], 'images', self.basename)) for frame in range(self.fs, self.fe + 1)]
            self.starttime = datetime.datetime.now()
            self.pfile = progressfile(self.folder, datetime.datetime.now(), 100)
            (self.pmfin, flag) = (0, 'Photon Maps') if sum(self.pmaps) else (1, 'Radiance Images')
            self.kivyrun = progressbar(os.path.join(self.folder, 'viprogress'), flag)

            if os.path.isfile("{}-{}.hdr".format(os.path.join(scene['viparams']['newdir'], 'images', self.basename), self.frame)):
               os.remove("{}-{}.hdr".format(os.path.join(scene['viparams']['newdir'], 'images', self.basename), self.frame))
                
            wm = context.window_manager
            self._timer = wm.event_timer_add(2, context.window)
            wm.modal_handler_add(self)                
            return {'RUNNING_MODAL'}
        else:
            self.report({'ERROR'}, "There is no camera in the scene or selected in the node. Create one for rpict image creation")
            return {'FINISHED'}

class NODE_OT_LiFC(bpy.types.Operator):            
    bl_idname = "node.livifc"
    bl_label = "LiVi False Colour Image"
    bl_description = "False colour an image with falsecolor"
    bl_register = True
    bl_undo = False
    nodeid = bpy.props.StringProperty()

    def execute(self, context):
        fcnode = bpy.data.node_groups[self.nodeid.split('@')[1]].nodes[self.nodeid.split('@')[0]] 
        fcnode.presim()
        imnode = fcnode.inputs['Image'].links[0].from_node 
        lmax = '-s {}'.format(fcnode.lmax) if fcnode.lmax else '-s a'
        scaling = '' if fcnode.nscale == '0' else '-log {}'.format(fcnode.decades) 
        mult = '-m {}'.format(fcnode.unitmult[fcnode.unit]) 
        legend = '-l {} -lw {} -lh {} {} {} {}'.format(fcnode.unitdict[fcnode.unit], fcnode.lw, fcnode.lh, lmax, scaling, mult) if fcnode.legend else ''
        bands = '-cb' if fcnode.bands else ''
        contour = '-cl {}'.format(bands) if fcnode.contour else ''
        divisions = '-n {}'.format(fcnode.divisions) if fcnode.divisions != 8 else ''
        
        for i, im in enumerate(imnode['images']): 
            fcim = os.path.join(context.scene['viparams']['newdir'], 'images', '{}-{}.hdr'.format(fcnode['basename'], i + context.scene['liparams']['fs']))
            ofile = bpy.path.abspath(fcnode.ofile) if os.path.isfile(bpy.path.abspath(fcnode.ofile)) and fcnode.overlay else bpy.path.abspath(im)
                        
            with open(fcim, 'w') as fcfile:
                if sys.platform == 'win32':
                    temp_file = os.path.join(context.scene['viparams']['newdir'], 'images', 'temp.hdr')
                    
                    with open(temp_file, 'w') as tfile:
                        Popen('pcond -e {} {}'.format(fcnode.disp, os.path.abspath(im)).split(), stdout = tfile)
                    
                    poverlay = '-p {}'.format(os.path.join(context.scene['viparams']['newdir'], 'images', 'temp.hdr')) if fcnode.contour and fcnode.overlay else ''
                    fccmd = 'falsecolor -i {} {} -pal {} {} {} {}'.format(os.path.abspath(im), poverlay, fcnode.coldict[fcnode.colour], legend, contour, divisions) 
                    fcrun = Popen(fccmd.split(), stdout=fcfile, stderr = PIPE) 
                    os.remove(temp_file)
                else:
                    poverlay = '-p <(pcond -e {0} {1})' .format(fcnode.disp, ofile) if fcnode.contour and fcnode.overlay else ''
                    fccmd = "bash -c 'falsecolor -i {} {} -pal {} {} {} {}'".format(bpy.path.abspath(im), poverlay, fcnode.coldict[fcnode.colour], legend, contour, divisions) 
                    fcrun = Popen(shlex.split(fccmd), stdout=fcfile, stderr = PIPE)

                logentry(fccmd)
               
                for line in fcrun.stderr:
                    logentry('Falsecolour error: {}'.format(line))
                    
            if fcim not in [i.filepath for i in bpy.data.images]:
                bpy.data.images.load(fcim)
            else:
                for i in bpy.data.images:
                    if bpy.path.abspath(i.filepath) == fcim:
                        i.reload()
                        [area.tag_redraw() for area in bpy.context.screen.areas if area and area.type == "IMAGE_EDITOR" and area.spaces[0].image == i]               
        fcnode.postsim()                              
        return {'FINISHED'}
    
class NODE_OT_LiGl(bpy.types.Operator):            
    bl_idname = "node.liviglare"
    bl_label = "LiVi Glare Node"
    bl_description = "Glare analysis node"
    bl_register = True
    bl_undo = False
    nodeid = bpy.props.StringProperty()
    
    def execute(self, context):
        scene = context.scene
        res = []
        reslists = []
        glnode = bpy.data.node_groups[self.nodeid.split('@')[1]].nodes[self.nodeid.split('@')[0]] 
        imnode = glnode.inputs[0].links[0].from_node
        glnode.presim()
        
        for i, im in enumerate(imnode['images']):
            glfile = os.path.join(scene['viparams']['newdir'], 'images', '{}-{}.hdr'.format(glnode['hdrname'], i + scene['liparams']['fs']))
            egcmd = 'evalglare {} -c {}'.format(('-u {0[0]} {0[1]} {0[2]}'.format(glnode.gc), '')[glnode.rand], glfile)

            with open(im, 'r') as hdrfile:
                egrun = Popen(egcmd.split(), stdin = hdrfile, stdout = PIPE)

            time = datetime.datetime(2014, 1, 1, imnode['coptions']['shour'], 0) + datetime.timedelta(imnode['coptions']['sdoy'] - 1) if imnode['coptions']['anim'] == '0' else \
                datetime.datetime(2014, 1, 1, int(imnode['coptions']['shour']), int(60*(imnode['coptions']['shour'] - int(imnode['coptions']['shour'])))) + datetime.timedelta(imnode['coptions']['sdoy'] - 1) + datetime.timedelta(hours = int(imnode['coptions']['interval']*i), seconds = int(60*(imnode['coptions']['interval']*i - int(imnode['coptions']['interval']*i))))
            
            with open(os.path.join(scene['viparams']['newdir'], 'images', "temp.glare"), "w") as glaretf:
                for line in egrun.stdout:
                    if line.decode().split(",")[0] == 'dgp':                            
                        glaretext = line.decode().replace(',', ' ').replace("#INF", "").split(' ')
                        res = [float(x) for x in glaretext[6:12]]
                        glaretf.write("{0:0>2d}/{1:0>2d} {2:0>2d}:{3:0>2d}\ndgp: {4:.2f}\ndgi: {5:.2f}\nugr: {6:.2f}\nvcp: {7:.2f}\ncgi: {8:.2f}\nLv: {9:.0f}\n".format(time.day, time.month, time.hour, time.minute, *res))
                        res.append(res)
                        reslists += [[str(i + scene['liparams']['fs']), 'Camera', 'Camera', 'DGP', '{0[0]}'.format(res)], [str(i + scene['liparams']['fs']), 'Camera', 'Camera', 'DGI', '{0[1]}'.format(res)], [str(i + scene['liparams']['fs']), 'Camera', 'Camera' 'UGR', '{0[2]}'.format(res)], [str(i + scene['liparams']['fs']), 'Camera', 'Camera', 'VCP', '{0[3]}'.format(res)], [str(i + scene['liparams']['fs']), 'Camera', 'Camera', 'CGI', '{[4]}'.format(res)], [str(i + scene['liparams']['fs']), 'Camera', 'Camera', 'LV', '{[5]}'.format(res)]]
            
            pcondcmd = "pcond -h+ -u 300 {}.hdr".format(os.path.join(context.scene['viparams']['newdir'], 'images', 'glare-'+str(i + scene['liparams']['fs'])))

            with open('{}.temphdr'.format(os.path.join(scene['viparams']['newdir'], 'images', 'glare')), 'w') as temphdr:
                Popen(pcondcmd.split(), stdout = temphdr).communicate()

            catcmd = "{0} {1}.glare".format(scene['viparams']['cat'], os.path.join(scene['viparams']['newdir'], 'images', 'temp'))
            catrun = Popen(catcmd, stdout = PIPE, shell = True)
            psigncmd = "psign -h {} -cb 0 0 0 -cf 1 1 1".format(int(0.04 * imnode.y))
            psignrun = Popen(psigncmd.split(), stdin = catrun.stdout, stdout = PIPE)
            pcompcmd = "pcompos {0}.temphdr 0 0 - {1} {2}".format(os.path.join(scene['viparams']['newdir'], 'images', 'glare'), imnode.x, imnode.y*550/800)

            with open("{}.hdr".format(os.path.join(scene['viparams']['newdir'], 'images', 'glare-'+str(i + scene['liparams']['fs']))), 'w') as ghdr:
                Popen(pcompcmd.split(), stdin = psignrun.stdout, stdout = ghdr).communicate()

            os.remove(os.path.join(scene['viparams']['newdir'], 'images', 'glare.temphdr'.format(i + scene['liparams']['fs'])))                               
        return {'FINISHED'}
        
class NODE_OT_LiViCalc(bpy.types.Operator):
    bl_idname = "node.livicalc"
    bl_label = "LiVi simulation"
    nodeid = bpy.props.StringProperty()
    bl_register = True
    bl_undo = False

    def invoke(self, context, event):
        scene = context.scene
        if viparams(self, scene):
            return {'CANCELLED'}
                    
        objmode()
        clearscene(scene, self)
        simnode = bpy.data.node_groups[self.nodeid.split('@')[1]].nodes[self.nodeid.split('@')[0]]
        simnode.presim()
        contextdict = {'Basic': 'LiVi Basic', 'Compliance': 'LiVi Compliance', 'CBDM': 'LiVi CBDM'}        
        
        # Set scene parameters
        scene['viparams']['visimcontext'] = contextdict[simnode['coptions']['Context']]
        scene['liparams']['fs'] = min((simnode['coptions']['fs'], simnode['goptions']['fs'])) 
        scene['liparams']['fe'] = max((simnode['coptions']['fe'], simnode['goptions']['fe'])) 
        scene['liparams']['cp'] = simnode['goptions']['cp']
        scene['liparams']['unit'] = simnode['coptions']['unit']
        scene['liparams']['type'] = simnode['coptions']['Type']
        scene['viparams']['vidisp'] = ''
        scene.frame_start, scene.frame_end = scene['liparams']['fs'], scene['liparams']['fe']
        
        simnode.sim(scene)

        for frame in range(scene['liparams']['fs'], scene['liparams']['fe'] + 1):
            if createradfile(scene, frame, self, simnode) == 'CANCELLED' or createoconv(scene, frame, self, simnode) == 'CANCELLED':
                return {'CANCELLED'}
        
        calcout = li_calc(self, simnode, livisimacc(simnode))

        if calcout == 'CANCELLED':
            return {'CANCELLED'}
        else:
            simnode['reslists'] = calcout
        if simnode['coptions']['Context'] != 'CBDM' and simnode['coptions']['Context'] != '1':
            scene.vi_display = 1

        scene['viparams']['vidisp'] = 'li'
        scene['viparams']['resnode'] = simnode.name
        scene['viparams']['restree'] = self.nodeid.split('@')[1]
        simnode.postsim()
        self.report({'INFO'},"Simulation is finished")
        return {'FINISHED'}
        
#class NODE_OT_LiVIGlare(bpy.types.Operator):
#    bl_idname = "node.liviglare"
#    bl_label = "LiVi glare"
#    bl_description = "Create a glare fisheye image from the Blender camera perspective"
#    bl_register = True
#    bl_undo = False
#    nodeid = bpy.props.StringProperty()
#
#    def modal(self, context, event):
#        if event.type == 'TIMER':
#            if self.egrun.poll() is not None: # If finished
#                if self.frame > self.scene['liparams']['fe']:
#                    self.reslists += [['All', 'Frames', '', 'Frames', ' '.join([str(f) for f in range(self.scene['liparams']['fs'], self.scene['liparams']['fe'] + 1)])]] + [['All', 'Camera', self.cam.name, ('DGP', 'DGI', 'UGR', 'VCP', 'CGI', 'LV')[ri], ' '.join([str(res) for res in r])] for ri, r in enumerate(zip(*self.res))]
#                    self.simnode['reslists'] = self.reslists
#                    self.simnode['frames'] = [f for f in range(self.scene['liparams']['fs'], self.scene['liparams']['fe'] + 1)]
#                    return {self.terminate()}
#
#                elif self.frame > self.frameold:
#                    self.percent = (self.frame - self.scene['liparams']['fs']) * 100
#                    self.frameold = self.frame
#                    os.remove(self.rpictfile)
#                    if self.simnode.pmap:
#                        amentry, pportentry, cpentry, cpfileentry = retpmap(self.simnode, self.frame, self.scene)
#                        pmcmd = ('mkpmap -bv+ +fo -apD 0.001 {0} -apg {1}-{2}.gpm {3} {4} {5} {1}-{2}.oct'.format(pportentry, self.scene['viparams']['filebase'], self.frame, self.simnode.pmapgno, cpentry, amentry))                   
#                        pmrun = Popen(pmcmd.split(), stderr = PIPE)
#                        for line in pmrun.stderr: 
#                            logentry('Photon map error: ', line)#        draw_image(self, self.ydiff * 0.1)
#                            if line.decode() in pmerrdict:
#                                self.report({'ERROR'}, pmerrdict[line.decode()])
#                                return {'CANCELLED'}
#                        rpictcmd = "rpict -w -e {7} -t 10 -vth -vh 180 -vv 180 -x 800 -y 800 -vd {0[0][2]:.3f} {0[1][2]} {0[2][2]} -vp {1[0]} {1[1]} {1[2]} {2} -ap {5} 50 {6} {3}-{4}.oct".format(-1*self.cam.matrix_world, self.cam.location, self.simnode['radparams'], self.scene['viparams']['filebase'], self.frame, '{}-{}.gpm'.format(self.scene['viparams']['filebase'], self.frame), cpfileentry, self.rpictfile)
#                    else:
#                        rpictcmd = "rpict -w -e {5} -t 10 -vth -vh 180 -vv 180 -x 800 -y 800 -vd {0[0][2]} {0[1][2]} {0[2][2]} -vp {1[0]} {1[1]} {1[2]} {2} {3}-{4}.oct".format(-1*self.cam.matrix_world, self.cam.location, self.simnode['radparams'], self.scene['viparams']['filebase'], self.frame, self.rpictfile)
#                    self.rprun = Popen(rpictcmd.split(), stdout = PIPE)                    
#                    self.egcmd = 'evalglare {} -c {}'.format(('-u 1 0 0', '')[sys.platform == 'win32'], os.path.join(self.scene['viparams']['newdir'], 'glare{}.hdr'.format(self.frame)))                    
#                    self.egrun = Popen(self.egcmd.split(), stdin = self.rprun.stdout, stdout = PIPE)
#                    return {'RUNNING_MODAL'}
#
#                time = datetime.datetime(2014, 1, 1, self.simnode['coptions']['shour'], 0) + datetime.timedelta(self.simnode['coptions']['sdoy'] - 1) if self.simnode['coptions']['anim'] == '0' else \
#                    datetime.datetime(2014, 1, 1, int(self.simnode['coptions']['shour']), int(60*(self.simnode['coptions']['shour'] - int(self.simnode['coptions']['shour'])))) + datetime.timedelta(self.simnode['coptions']['sdoy'] - 1) + datetime.timedelta(hours = int(self.simnode['coptions']['interval']*(self.frame-self.scene['liparams']['fs'])), seconds = int(60*(self.simnode['coptions']['interval']*(self.frame-self.scene['liparams']['fs']) - int(self.simnode['coptions']['interval']*(self.frame-self.scene['liparams']['fs'])))))
#                with open(self.scene['viparams']['filebase']+".glare", "w") as glaretf:
#                    for line in self.egrun.stdout:
#                        if line.decode().split(",")[0] == 'dgp':                            
#                            glaretext = line.decode().replace(',', ' ').replace("#INF", "").split(' ')
#                            res = [float(x) for x in glaretext[6:12]]
#                            glaretf.write("{0:0>2d}/{1:0>2d} {2:0>2d}:{3:0>2d}\ndgp: {4:.2f}\ndgi: {5:.2f}\nugr: {6:.2f}\nvcp: {7:.2f}\ncgi: {8:.2f}\nLv: {9:.0f}\n".format(time.day, time.month, time.hour, time.minute, *res))
#                            self.res.append(res)
#                            self.reslists += [[str(self.frame), 'Camera', self.cam.name, 'DGP', '{0[0]}'.format(res)], [str(self.frame), 'Camera', self.cam.name, 'DGI', '{0[1]}'.format(res)], [str(self.frame), 'Camera', self.cam.name, 'UGR', '{0[2]}'.format(res)], [str(self.frame), 'Camera', self.cam.name, 'VCP', '{0[3]}'.format(res)], [str(self.frame), 'Camera', self.cam.name, 'CGI', '{[4]}'.format(res)], [str(self.frame), 'Camera', self.cam.name, 'LV', '{[5]}'.format(res)]]
#                
#                pcondcmd = "pcond -u 300 {0}.hdr".format(os.path.join(self.scene['viparams']['newdir'], 'glare'+str(self.frame)))
#                with open('{}.temphdr'.format(os.path.join(self.scene['viparams']['newdir'], 'glare'+str(self.frame))), 'w') as temphdr:
#                    Popen(pcondcmd.split(), stdout = temphdr).communicate()
#                catcmd = "{0} {1}.glare".format(self.scene['viparams']['cat'], self.scene['viparams']['filebase'])
#                catrun = Popen(catcmd, stdout = PIPE, shell = True)
#                psigncmd = "psign -h 32 -cb 0 0 0 -cf 40 40 40"
#                psignrun = Popen(psigncmd.split(), stdin = catrun.stdout, stdout = PIPE)
#                pcompcmd = "pcompos {0}.temphdr 0 0 - 800 550".format(os.path.join(self.scene['viparams']['newdir'], 'glare'+str(self.frame)))
#                with open("{}.hdr".format(os.path.join(self.scene['viparams']['newdir'], 'glare'+str(self.frame))), 'w') as ghdr:
#                    Popen(pcompcmd.split(), stdin = psignrun.stdout, stdout = ghdr).communicate()
#                os.remove(os.path.join(self.scene['viparams']['newdir'], 'glare{}.temphdr'.format(self.frame)))
#
#                if 'glare{}.hdr'.format(self.frame) in bpy.data.images:
#                    bpy.data.images['glare{}.hdr'.format(self.frame)].filepath = os.path.join(self.scene['viparams']['newdir'], 'glare{}.hdr'.format(self.frame))
#                    bpy.data.images['glare{}.hdr'.format(self.frame)].reload()
#                else:
#                    bpy.data.images.load(os.path.join(self.scene['viparams']['newdir'], 'glare{}.hdr'.format(self.frame)))
#                self.frame += 1
#                return {'RUNNING_MODAL'}
#            else:
#                with open(self.rpictfile) as rpictfile:
#                    for line in rpictfile.readlines()[::-1]:
#                        if '%' in line:
#                            for lineentry in line.split():
#                                if '%' in lineentry:
#                                    self.percent = (float(lineentry.strip('%')) + (self.frame - self.scene['liparams']['fs']) * 100)/self.frames
#                            break
#     
#                if self.percent:
#                    if self.pfile.check(self.percent) == 'CANCELLED':                                    
#                        return {self.terminate()}
#                
#                return {'PASS_THROUGH'}
#        else:
#            return {'PASS_THROUGH'}
#            
#    def terminate(self):
#        nodecolour(self.simnode, 0)
#        self.kivyrun.kill() 
#        self.simnode.run = 0
#
#        if self.egrun.poll() == None:                          
#            self.egrun.kill()
#            
#        self.rprun.kill()        
#        self.simnode.postsim()
#        return 'FINISHED'
#
#    def execute(self, context):        
#        self.scene = bpy.context.scene
#        self.cam = self.scene.camera
#        
#        if self.cam:
#            self.percent = 0
#            self.reslists = []
#            self.res = []
#            self.rpictfile = os.path.join(self.scene['viparams']['newdir'], 'rpictprogress')
#            self.simnode = bpy.data.node_groups[self.nodeid.split('@')[1]].nodes[self.nodeid.split('@')[0]]
#            self.simnode.presim()
#            self.simnode.run = 1
#            nodecolour(self.simnode, 1)
#            self.scene['liparams']['fs'] = min([c['fs'] for c in (self.simnode['goptions'], self.simnode['coptions'])])
#            self.scene['liparams']['fe'] = max([c['fe'] for c in (self.simnode['goptions'], self.simnode['coptions'])])
#            self.frames = self.scene['liparams']['fe'] - self.scene['liparams']['fs'] + 1
#            self.frame = self.scene['liparams']['fs']
#            self.frameold = self.frame
#            
#            for frame in range(self.scene['liparams']['fs'], self.scene['liparams']['fe'] + 1):
#                createradfile(self.scene, frame, self, self.simnode)
#                createoconv(self.scene, frame, self, self.simnode)
#                
#            if self.simnode.pmap:
#                amentry, pportentry, cpentry, cpfileentry = retpmap(self.simnode, self.frame, self.scene)
#                pmcmd = ('mkpmap -bv+ +fo -apD 0.001 {0} -apg {1}-{2}.gpm {3} {4} {5} {1}-{2}.oct'.format(pportentry, self.scene['viparams']['filebase'], self.frame, self.simnode.pmapgno, cpentry, amentry))                   
#                pmrun = Popen(pmcmd.split(), stderr = PIPE)
#                for line in pmrun.stderr: 
#                    logentry(line)
#                    if line.decode() in pmerrdict:
#                        self.report({'ERROR'}, pmerrdict[line.decode()])
#                        return {'FINISHED'}
#                rpictcmd = "rpict -w -e {7} -t 1 -vth -vh 180 -vv 180 -x 800 -y 800 -vd {0[0][2]:.3f} {0[1][2]} {0[2][2]} -vp {1[0]} {1[1]} {1[2]} {2} -ap {5} 50 {6} {3}-{4}.oct".format(-1*self.cam.matrix_world, self.cam.location, self.simnode['radparams'], self.scene['viparams']['filebase'], self.frame, '{}-{}.gpm'.format(self.scene['viparams']['filebase'], self.frame), cpfileentry, self.rpictfile)
#            else:
#                rpictcmd = "rpict -w -vth -vh 180 -e {5} -t 1 -vv 180 -x 800 -y 800 -vd {0[0][2]:.3f} {0[1][2]} {0[2][2]} -vp {1[0]} {1[1]} {1[2]} {2} {3}-{4}.oct".format(-1*self.cam.matrix_world, self.cam.location, self.simnode['radparams'], self.scene['viparams']['filebase'], self.frame, self.rpictfile)
#                rpiececmd = "rpiece -X 1 -Y {6} -o {} -vth -vh 180 -e {5} -t 1 -vv 180 -x 800 -y 800 -vd {0[0][2]:.3f} {0[1][2]} {0[2][2]} -vp {1[0]} {1[1]} {1[2]} {2} {3}-{4}.oct".format(-1*self.cam.matrix_world, self.cam.location, self.simnode['radparams'], self.scene['viparams']['filebase'], self.frame, self.rpictfile)
#            self.starttime = datetime.datetime.now()
#            self.pfile = progressfile(self.scene['viparams']['newdir'], datetime.datetime.now(), 100)
#            self.kivyrun = progressbar(os.path.join(self.scene['viparams']['newdir'], 'viprogress'), 'Glare')
#            self.rprun = Popen(rpictcmd.split(), stdout=PIPE, stderr = PIPE)
#            egcmd = "evalglare {} -c {}".format(('-u 1 0 0', '')[sys.platform == 'win32'], os.path.join(self.scene['viparams']['newdir'], 'glare{}.hdr'.format(self.frame)))
#            self.egrun = Popen(egcmd.split(), stdin = self.rprun.stdout, stdout=PIPE, stderr = PIPE)
#            wm = context.window_manager
#            self._timer = wm.event_timer_add(10, context.window)
#            wm.modal_handler_add(self)
#            return {'RUNNING_MODAL'}
#        else:
#            self.report({'ERROR'}, "There is no camera in the scene. Create one for glare analysis")
#            return {'FINISHED'}

class IES_Select(bpy.types.Operator, io_utils.ImportHelper):
    bl_idname = "livi.ies_select"
    bl_label = "Select IES file"
    bl_description = "Select the lamp IES file"
    filename = ""
    filename_ext = ".ies; .IES"
    filter_glob = bpy.props.StringProperty(default="*.ies; *.IES", options={'HIDDEN'})
    bl_register = True
    bl_undo = True

    def draw(self,context):
        layout = self.layout
        row = layout.row()
        row.label(text="Open an IES File with the file browser", icon='WORLD_DATA')

    def execute(self, context):
        lamp = bpy.context.active_object
        lamp['ies_name'] = self.filepath
        return {'FINISHED'}

    def invoke(self,context,event):
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}

class NODE_OT_ESOSelect(NODE_OT_FileSelect):
    bl_idname = "node.esoselect"
    bl_label = "Select EnVi results file"
    bl_description = "Select the EnVi results file to process"
    filename_ext = ".eso"
    filter_glob = bpy.props.StringProperty(default="*.eso", options={'HIDDEN'})
    nodeprop = 'resfilename'
    filepath = bpy.props.StringProperty(subtype='FILE_PATH', options={'HIDDEN', 'SKIP_SAVE'})
    fextlist = ("eso")
    nodeid = bpy.props.StringProperty()

class NODE_OT_IDFSelect(NODE_OT_FileSelect):
    bl_idname = "node.idfselect"
    bl_label = "Select EnergyPlus input file"
    bl_description = "Select the EnVi input file to process"
    filename_ext = ".idf"
    filter_glob = bpy.props.StringProperty(default="*.idf", options={'HIDDEN'})
    nodeprop = 'idffilename'
    filepath = bpy.props.StringProperty(subtype='FILE_PATH', options={'HIDDEN', 'SKIP_SAVE'})
    fextlist = ("idf")
    nodeid = bpy.props.StringProperty()
    
class NODE_OT_ASCImport(bpy.types.Operator, io_utils.ImportHelper):
    bl_idname = "node.ascimport"
    bl_label = "Select ESRI Grid file"
    bl_description = "Select the ESRI Grid file to process"
    filename = ""
    filename_ext = ".asc"
    filter_glob = bpy.props.StringProperty(default="*.asc", options={'HIDDEN'})
    bl_register = True
    bl_undo = False
    nodeid = bpy.props.StringProperty()

    def draw(self,context):
        layout = self.layout
        row = layout.row()
        row.label(text="Open an asc file with the file browser", icon='WORLD_DATA')

    def execute(self, context):
        node = bpy.data.node_groups[self.nodeid.split('@')[1]].nodes[self.nodeid.split('@')[0]]
        startxs, startys, vlen = [], [], 0
        ascfiles = [self.filepath] if node.single else [os.path.join(os.path.dirname(os.path.realpath(self.filepath)), file) for file in os.listdir(os.path.dirname(os.path.realpath(self.filepath))) if file.endswith('.asc')]
        obs = []
        headerdict = {'ncols': 0, 'nrows': 0, 'xllcorner': 0, 'yllcorner': 0, 'cellsize': 0, 'NODATA_value': 0}

        for file in ascfiles:
            basename = file.split(os.sep)[-1].split('.')[0]
            me = bpy.data.meshes.new("{} mesh".format(basename))
            bm = bmesh.new()
            l = 0
            with open(file, 'r') as ascfile:
                lines = ascfile.readlines()
                
                while len(lines[l].split()) == 2:
                    if lines[l].split()[0] in headerdict:
                        headerdict[lines[l].split()[0]] = eval(lines[l].split()[1])
                    l += 1
   
                vlen = headerdict['nrows'] * headerdict['ncols']                   
                startxs.append(headerdict['xllcorner'])
                startys.append(headerdict['yllcorner'])                
                x, y = 0, headerdict['nrows']
                
                for l, line in enumerate(lines[l:]):   
                    for zval in line.split():
                        [bm.verts.new((x * headerdict['cellsize'], y * headerdict['cellsize'], float(zval)))]                         
                        x += 1
                    x = 0
                    y -=1
            
            bm.verts.ensure_lookup_table()
            faces = [(i+1, i, i+headerdict['ncols'], i+headerdict['ncols'] + 1) for i in range(0, vlen - headerdict['ncols']) if (i+1)%headerdict['ncols']]
            [bm.faces.new([bm.verts[fv] for fv in face]) for face in faces]
            
            if node.clear_nodata == '1':
                bmesh.ops.delete(bm, geom = [v for v in bm.verts if v.co[2] == headerdict['NODATA_value']], context = 1)
            
            elif node.clear_nodata == '0':
                for v in bm.verts:
                    if v.co[2] == headerdict['NODATA_value']:
                        v.co[2] = 0
                        
            bm.to_mesh(me)
            bm.free()
            ob = bpy.data.objects.new(basename, me)
            bpy.context.scene.objects.link(ob)
            obs.append(ob)

        minstartx,  minstarty = min(startxs), min(startys)

        for o, ob in enumerate(obs):
            ob.location = (startxs[o] - minstartx, startys[o] - minstarty, 0)
            
        return {'FINISHED'}
        
    def invoke(self,context,event):
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}
    
class NODE_OT_CSVExport(bpy.types.Operator, io_utils.ExportHelper):
    bl_idname = "node.csvexport"
    bl_label = "Export a CSV file"
    bl_description = "Select the CSV file to export"
    filename = "results"
    filename_ext = ".csv"
    filter_glob = bpy.props.StringProperty(default="*.csv", options={'HIDDEN'})
    bl_register = True
    bl_undo = True
    nodeid = bpy.props.StringProperty()

    def draw(self,context):
        layout = self.layout
        row = layout.row()
        row.label(text="Specify the CSV export file with the file browser", icon='WORLD_DATA')

    def execute(self, context):
        node = bpy.data.node_groups[self.nodeid.split('@')[1]].nodes[self.nodeid.split('@')[0]]
        resstring = ''
        resnode = node.inputs['Results in'].links[0].from_node
        rl = resnode['reslists']
        zrl = list(zip(*rl))

        if len(set(zrl[0])) > 1 and node.animated:
            resstring = ''.join(['{} {},'.format(r[2], r[3]) for r in rl if r[0] == 'All']) + '\n'
            metriclist = list(zip(*[r.split() for ri, r in enumerate(zrl[4]) if zrl[0][ri] == 'All']))

        else:
            resstring = ''.join(['{} {} {},'.format(r[0], r[2], r[3]) for r in rl if r[0] != 'All']) + '\n'
            metriclist = list(zip(*[r.split() for ri, r in enumerate(zrl[4]) if zrl[0][ri] != 'All']))

        for ml in metriclist:
            resstring += ''.join(['{},'.format(m) for m in ml]) + '\n'

        resstring += '\n'

        with open(self.filepath, 'w') as csvfile:
            csvfile.write(resstring)
        return {'FINISHED'}

    def invoke(self,context,event):
        if self.filepath.split('.')[-1] not in ('csv', 'CSV'):
            self.filepath = os.path.join(context.scene['viparams']['newdir'], context.scene['viparams']['filebase'] + '.csv')            
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}

class NODE_OT_PLYExport(bpy.types.Operator):
    bl_idname = "node.plyexport"
    bl_label = "Export a Ply file for I-Simpa acoustic analysis"
    bl_description = "Export the Ply file"
    bl_register = True
    bl_undo = False
    
    def execute(self, context):
        scene = context.scene
        vsum, fsum, msum, vtext, ftext, mtext = 0, 0, 0, '', '', ''
        
        for ob in [o for o in bpy.data.objects if o.select and o.type == 'MESH']:
            bm = bmesh.new()
            tm = ob.to_mesh(scene = scene, apply_modifiers = True, settings = 'PREVIEW') 
            bm.from_mesh(tm)
            bpy.data.meshes.remove(tm)            
            bmesh.ops.connect_verts(bm, verts = bm.verts[:], faces_exclude = [], check_degenerate = 1)
            bmesh.ops.triangulate(bm, faces=bm.faces[:], quad_method=3, ngon_method=0)
            bm.transform(ob.matrix_world)
            vsum += len(bm.verts)
            fsum += len(bm.faces)
            msum += len(ob.data.materials)
            vtext += "\n".join(['{0[0]:.4f} {0[1]:.4f} {0[2]:.4f}'.format(v.co) for v in bm.verts]) + "\n" 
            ftext += "\n".join([' '.join([str(len(f.verts))] + [str(fv.index) for fv in f.verts] + [str(f.material_index)]) for f in bm.faces]) + "\n"
            mtext += "\n".join([' '.join([str(len(m.name))] + [str(ord(c)) for c in m.name]) for m in ob.data.materials]) + "\n"
            bm.free()
        
        if vtext:
            with open(bpy.path.abspath(context.node.ofile), 'w') as ply_file:
                ply_file.write('ply\nformat ascii 1.0\nelement vertex {}\nproperty float x\nproperty float y\nproperty float z\nelement face {}\nproperty list uchar int vertex_indices\nproperty int layer_id\nelement layer {}\nproperty list uchar uchar layer_name\nend_header\n'.format(vsum, fsum, msum))
                ply_file.write(vtext + ftext + mtext)
            return {'FINISHED'}
        
        else:
            self.report({'ERROR'}, 'No valid geometry was selected')
            return {'CANCELLED'}
        
class NODE_OT_TGPolyExport(bpy.types.Operator):
    bl_idname = "node.polyexport"
    bl_label = "Create a poly file for tetgen meshing"
    bl_description = "Create the poly file"
    bl_register = True
    bl_undo = False
    
    def execute(self, context):
        scene = context.scene
        vsum, fsum, msum, vtext, ftext, mtext = 0, 0, 0, '', '', ''
        
        for ob in [o for o in bpy.data.objects if o.select and o.type == 'MESH']:
            bm = bmesh.new()
            tm = ob.to_mesh(scene = scene, apply_modifiers = True, settings = 'PREVIEW') 
            bm.from_mesh(tm)
            bpy.data.meshes.remove(tm)  
            bm.transform(ob.matrix_world)
#            bmesh.ops.connect_verts(bm, verts = bm.verts[:], faces_exclude = [], check_degenerate = 1)
#            bmesh.ops.triangulate(bm, faces=bm.faces[:], quad_method=3, ngon_method=0)
            bmesh.ops.connect_verts_nonplanar(bm, angle_limit = 0.01, faces = bm.faces)
            
            vsum += len(bm.verts)
            fsum += len(bm.faces)
            msum += len(ob.data.materials)
            vtext += "\n".join(['{0} {1[0]:.4f} {1[1]:.4f} {1[2]:.4f}'.format(vi + 1, v.co) for vi, v in enumerate(bm.verts)]) + "\n" 
            ftext += "\n".join(["1 0 {}\n".format(f.index) + ' '.join([str(len(f.verts))] + [str(fv.index + 1) for fv in f.verts] + [str(f.material_index + 1)]) for f in bm.faces]) + "\n"
            mtext += "\n".join([' '.join([str(len(m.name))] + [str(ord(c)) for c in m.name]) for m in ob.data.materials]) + "\n"
            bm.free()
        
        if vtext:
            with open(bpy.path.abspath(context.node.ofile), 'w') as poly_file:
                poly_file.write('#Part 1 - node list\n{} 3 0 0\n'.format(vsum))
                poly_file.write(vtext)
                poly_file.write('#Part 2 - facet list\n{} 1\n'.format(fsum))
                poly_file.write(ftext)
    
            return {'FINISHED'}
        
        else:
            self.report({'ERROR'}, 'No valid geometry was selected')
            return {'CANCELLED'}
        
class NODE_OT_TextUpdate(bpy.types.Operator):
    bl_idname = "node.textupdate"
    bl_label = "Update a text file"
    bl_description = "Update a text file"

    nodeid = bpy.props.StringProperty()

    def execute(self, context):
        tenode = bpy.data.node_groups[self.nodeid.split('@')[1]].nodes[self.nodeid.split('@')[0]]
        tenode.textupdate(tenode['bt'])
        return {'FINISHED'}

class NODE_OT_TextExport(bpy.types.Operator, io_utils.ExportHelper):
    bl_idname = "node.textexport"
    bl_label = "Export a text file"
    bl_description = "Select the text file to export"
    filename = ""
    filename_ext = ".txt"
    filter_glob = bpy.props.StringProperty(default="*.txt", options={'HIDDEN'})
    bl_register = True
    bl_undo = True
    nodeid = bpy.props.StringProperty()

    def draw(self,context):
        layout = self.layout
        row = layout.row()
        row.label(text="Specify the Text export file with the file browser", icon='WORLD_DATA')

    def execute(self, context):
        hostnode = bpy.data.node_groups[self.nodeid.split('@')[1]].nodes[self.nodeid.split('@')[0]]
        textsocket = bpy.data.node_groups[self.nodeid.split('@')[1]].nodes[self.nodeid.split('@')[0]].inputs['Text in'].links[0].from_socket
        resstring = '\n'.join(textsocket['Text'])
        with open(self.filepath, 'w') as textfile:
            textfile.write(resstring)
        if hostnode.etoggle:
            if self.filepath not in [im.filepath for im in bpy.data.texts]:
                bpy.data.texts.load(self.filepath)

            imname = [im.name for im in bpy.data.texts if im.filepath == self.filepath][0]
            text = bpy.data.texts[imname]
            for area in bpy.context.screen.areas:
                if area.type == 'TEXT_EDITOR':
                    area.spaces.active.text = text
                    ctx = bpy.context.copy()
                    ctx['edit_text'] = text
                    ctx['area'] = area
                    ctx['region'] = area.regions[-1]
                    bpy.ops.text.resolve_conflict(ctx, resolution = 'RELOAD')

        return {'FINISHED'}

    def invoke(self,context,event):
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}

class NODE_OT_EnGExport(bpy.types.Operator):
    bl_idname = "node.engexport"
    bl_label = "VI-Suite export"
    bl_context = "scene"
    nodeid = bpy.props.StringProperty()

    def invoke(self, context, event):
        scene = context.scene
        if viparams(self, scene):
            return {'CANCELLED'}
        scene['viparams']['vidisp'] = ''
        node = bpy.data.node_groups[self.nodeid.split('@')[1]].nodes[self.nodeid.split('@')[0]]
        node.preexport(scene)
        pregeo(self)
        node.postexport()
        return {'FINISHED'}

class MAT_EnVi_Node(bpy.types.Operator):
    bl_idname = "material.envi_node"
    bl_label = "EnVi Material export"
    nodeid = bpy.props.StringProperty()

    def invoke(self, context, event):
        if not context.material.envi_nodes:
            bpy.ops.node.new_node_tree(type='EnViMatN', name = context.material.name) 
            context.material.envi_nodes = bpy.data.node_groups[context.material.name]
            context.material.envi_nodes.nodes.new('EnViCon')
            context.material.envi_nodes['envi_con_type'] = 'None'
            context.material.envi_nodes.nodes[0].active = True
            context.material.envi_nodes['enmatparams'] = {'airflow': 0, 'boundary': 0, 'tm': 0}

        elif context.material.name != context.material.envi_nodes.name and context.material.envi_nodes.name in [m.name for m in bpy.data.materials]:
            context.material.envi_nodes = context.material.envi_nodes.copy()
            context.material.envi_nodes.name = context.material.name
        return {'FINISHED'}
    
class NODE_EnVi_UV(bpy.types.Operator):
    bl_idname = "node.envi_uv"
    bl_label = "EnVi Material U-Value Calculation"

    def execute(self, context):
        resists = []
        node = context.node
        lsock = node.inputs['Outer layer']
                
        while lsock.links:
            resists.append(lsock.links[0].from_node.ret_resist())
            lsock = lsock.links[0].from_node.inputs['Layer']

        node.uv = '{:.3f}'.format(1/(sum(resists) + 0.12 + 0.08))
        return {'FINISHED'}
    
class NODE_OT_EnExport(bpy.types.Operator, io_utils.ExportHelper):
    bl_idname = "node.enexport"
    bl_label = "Export"
    bl_description = "Export the scene to the EnergyPlus file format"
    bl_register = True
    bl_undo = False
    nodeid = bpy.props.StringProperty()

    def invoke(self, context, event):
        scene = context.scene
        
        if viparams(self, scene):
            return {'CANCELLED'}
        
        scene['viparams']['vidisp'] = ''
        node = bpy.data.node_groups[self.nodeid.split('@')[1]].nodes[self.nodeid.split('@')[0]]
        (scene['enparams']['fs'], scene['enparams']['fe']) = (node.fs, node.fe) if node.animated else (scene.frame_current, scene.frame_current)
        locnode = node.inputs['Location in'].links[0].from_node
        
        if not os.path.isfile(locnode.weather):
            self.report({'ERROR'}, 'Location node weather file is not valid')
            node.use_custom_color = 1
            return {'CANCELLED'}
        
        node.preexport(scene)
        
        for frame in range(node.fs, node.fe + 1):
            scene.frame_set(frame)
            shutil.copyfile(locnode.weather, os.path.join(scene['viparams']['newdir'], "in{}.epw".format(frame)))
        scene.frame_set(node.fs)

        if bpy.context.active_object and not bpy.context.active_object.hide:
            if bpy.context.active_object.type == 'MESH':
                bpy.ops.object.mode_set(mode = 'OBJECT')

        enpolymatexport(self, node, locnode, envi_mats, envi_cons)
        node.bl_label = node.bl_label[1:] if node.bl_label[0] == '*' else node.bl_label
        node.exported, node.outputs['Context out'].hide = True, False
        node.postexport()
        return {'FINISHED'}
    
class NODE_OT_EnSim(bpy.types.Operator):
    bl_idname = "node.ensim"
    bl_label = "Simulate"
    bl_description = "Run EnergyPlus"
    bl_register = True
    bl_undo = False
    nodeid = bpy.props.StringProperty()

    def modal(self, context, event):
        if self.pfile.check(self.percent) == 'CANCELLED':                                    
            return {self.terminate('CANCELLED', context)}
        
        while sum([esim.poll() is None for esim in self.esimruns]) < self.processors and self.e < self.lenframes:
            self.esimruns.append(Popen(self.esimcmds[self.e].split(), stderr = PIPE))
            self.e += 1
    
        if event.type == 'TIMER':
            if len(self.esimruns) > 1:
                self.percent = 100 * sum([esim.poll() is not None for esim in self.esimruns])/self.lenframes 
            else:
                try:
                    with open(os.path.join(self.nd, '{}{}out.eso'.format(self.resname, self.frame)), 'r') as resfile:
                        for resline in [line for line in resfile.readlines()[::-1] if line.split(',')[0] == '2' and len(line.split(',')) == 9]:
                            self.percent = 100 * int(resline.split(',')[1])/(self.simnode.dedoy - self.simnode.dsdoy)
                            break
                except:
                    pass
#                    logentry('There was an error in the EnVi simulation. Check the error log in the text editor')
#                    return {self.terminate('CANCELLED', context)}
                
            if all([esim.poll() is not None for esim in self.esimruns]) and self.e == self.lenframes:
                for fname in [fname for fname in os.listdir('.') if fname.split(".")[0] == self.simnode.resname]:
                    os.remove(os.path.join(self.nd, fname))

                for f in range(self.frame, self.frame + self.e):
                    nfns = [fname for fname in os.listdir('.') if fname.split(".")[0] == "{}{}out".format(self.resname, f)]
                    
                    for fname in nfns:
                        os.rename(os.path.join(self.nd, fname), os.path.join(self.nd,fname.replace("eplusout", self.simnode.resname)))
                      
                    efilename = "{}{}out.err".format(self.resname, f)
                    
                    if efilename not in [im.name for im in bpy.data.texts]:
                        bpy.data.texts.load(os.path.join(self.nd, efilename))
                    else:
                        bpy.data.texts[efilename].filepath = os.path.join(self.nd, efilename)
                    
                    if '** Severe  **' in bpy.data.texts[efilename]:
                        self.report({'ERROR'}, "Fatal error reported in the {} file. Check the file in Blender's text editor".format(efilename))
                        return {self.terminate('CANCELLED', context)}
        
                    if 'EnergyPlus Terminated--Error(s) Detected' in self.esimruns[f - self.frame].stderr.read().decode() or not [f for f in nfns if f.split(".")[1] == "eso"] or self.simnode.run == 0:
                        errtext = "There is no results file. Check you have selected results outputs and that there are no errors in the .err file in the Blender text editor." if not [f for f in nfns if f.split(".")[1] == "eso"] else "There was an error in the input IDF file. Check the *.err file in Blender's text editor."
                        self.report({'ERROR'}, errtext)
                        return {self.terminate('CANCELLED', context)}
                                
                self.report({'INFO'}, "Calculation is finished.") 
                return {self.terminate('FINISHED', context)}
            else:
                return {'PASS_THROUGH'} 
        else:
            return {'PASS_THROUGH'}   
        
    def invoke(self, context, event):
        scene = context.scene
         
        if viparams(self, scene):
            return {'CANCELLED'}
        
        if shutil.which('energyplus') is None:
            self.report({'ERROR'}, "Energyplus binary is not executable")
            return {'CANCELLED'}
        
        self.frame = scene['enparams']['fs']
        self.frames = range(scene['enparams']['fs'], scene['enparams']['fe'] + 1)
        self.lenframes = len(self.frames)             
        context.scene['viparams']['visimcontext'] = 'EnVi'
        self.pfile = progressfile(scene['viparams']['newdir'], datetime.datetime.now(), 100)
        self.kivyrun = progressbar(os.path.join(scene['viparams']['newdir'], 'viprogress'), 'EnergyPlus Results')
        wm = context.window_manager
        self._timer = wm.event_timer_add(1, context.window)
        wm.modal_handler_add(self)
        self.simnode = bpy.data.node_groups[self.nodeid.split('@')[1]].nodes[self.nodeid.split('@')[0]]
        self.connode = self.simnode.inputs[0].links[0].from_node.name
        self.simnode.presim(context)        
        self.expand = "-x" if scene['viparams'].get('hvactemplate') else ""
        self.resname = (self.simnode.resname, 'eplus')[self.simnode.resname == '']
        os.chdir(scene['viparams']['newdir'])
        self.esimcmds = ["energyplus {0} -w in{1}.epw -p {2} in{1}.idf".format(self.expand, frame, ('{}{}'.format(self.resname, frame))) for frame in self.frames] 
        self.esimruns = []
        self.simnode.run = 1
        self.processors = self.simnode.processors if self.simnode.mp else 1
        self.percent = 0
        self.e = 0
        self.nd = scene['viparams']['newdir']
        return {'RUNNING_MODAL'}
    
    def terminate(self, condition, context):
        self.simnode.postsim(self, condition)
        if condition == 'FINISHED':
            context.scene['viparams']['resnode'] = self.nodeid
            context.scene['viparams']['connode'] = '{}@{}'.format(self.connode, self.nodeid.split('@')[1])
            context.scene['viparams']['vidisp'] = 'en'

        self.kivyrun.kill() 
             
        for es in self.esimruns:
            if es.poll() is None:
                es.kill()
                
        return condition
    
class VIEW3D_OT_EnDisplay(bpy.types.Operator):
    bl_idname = "view3d.endisplay"
    bl_label = "EnVi display"
    bl_description = "Display the EnVi results"
    bl_options = {'REGISTER'}
#    bl_undo = True
    _handle = None
    disp =  bpy.props.IntProperty(default = 1)
    
    @classmethod
    def poll(cls, context):
        return context.area.type   == 'VIEW_3D' and \
               context.region.type == 'WINDOW'

    def modal(self, context, event):
        scene = context.scene
        if event.type != 'INBETWEEN_MOUSEMOVE' and context.region and context.area.type == 'VIEW_3D' and context.region.type == 'WINDOW':  
            if context.scene.vi_display == 0 or context.scene['viparams']['vidisp'] != 'enpanel':
                context.scene['viparams']['vidisp'] = 'en'
                bpy.types.SpaceView3D.draw_handler_remove(self._handle_en_disp, 'WINDOW')
                
                try:
                    bpy.types.SpaceView3D.draw_handler_remove(self._handle_air, 'WINDOW')
                except:
                    pass
                
                for o in [o for o in scene.objects if o.get('VIType') and o['VIType'] in ('envi_temp', 'envi_hum', 'envi_heat', 'envi_cool', 'envi_co2', 'envi_shg', 'envi_ppd', 'envi_pmv', 'envi_aheat', 'envi_acool', 'envi_hrheat')]:
                    for oc in o.children:                        
                        [scene.objects.unlink(oc) for oc in o.children]
                        bpy.data.objects.remove(oc)                    
                    scene.objects.unlink(o)
                    bpy.data.objects.remove(o)

                context.area.tag_redraw()
                return {'CANCELLED'}

            mx, my, redraw = event.mouse_region_x, event.mouse_region_y, 0
            
            if self.dhscatter.spos[0] < mx < self.dhscatter.epos[0] and self.dhscatter.spos[1] < my < self.dhscatter.epos[1]:
                if self.dhscatter.hl != (0, 1, 1, 1):  
                    self.dhscatter.hl = (0, 1, 1, 1)
                    redraw = 1  
                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.dhscatter.press = 1
                        self.dhscatter.move = 0
                        return {'RUNNING_MODAL'}
                    elif event.value == 'RELEASE':
                        if not self.dhscatter.move:
                            self.dhscatter.expand = 0 if self.dhscatter.expand else 1
                        self.dhscatter.press = 0
                        self.dhscatter.move = 0
                        context.area.tag_redraw()
                        return {'RUNNING_MODAL'}
                
                elif event.type == 'ESC':
                    bpy.data.images.remove(bpy.data.images[self.dhscatter.gimage])
                    self.dhscatter.plt.close()
                    bpy.types.SpaceView3D.draw_handler_remove(self._handle_en_disp, 'WINDOW')
                    context.area.tag_redraw()
                    return {'CANCELLED'}
                    
                elif self.dhscatter.press and event.type == 'MOUSEMOVE':
                     self.dhscatter.move = 1
                     self.dhscatter.press = 0
        
            elif self.dhscatter.expand and self.dhscatter.lspos[0] < mx < self.dhscatter.lepos[0] and self.dhscatter.lspos[1] < my < self.dhscatter.lepos[1] and abs(self.dhscatter.lepos[0] - mx) > 20 and abs(self.dhscatter.lspos[1] - my) > 20: 
                self.dhscatter.hl = (1, 1, 1, 1)
                if event.type == 'LEFTMOUSE' and event.value == 'PRESS' and self.dhscatter.expand and self.dhscatter.lspos[0] < mx < self.dhscatter.lepos[0] and self.dhscatter.lspos[1] < my < self.dhscatter.lspos[1] + 0.9 * self.dhscatter.ydiff:
                    self.dhscatter.show_plot()
                    context.area.tag_redraw()
                    return {'RUNNING_MODAL'}    
                    
            elif self.dhscatter.hl != (1, 1, 1, 1):
                self.dhscatter.hl = (1, 1, 1, 1)
                redraw = 1
                
            if self.table.spos[0] < mx < self.table.epos[0] and self.table.spos[1] < my < self.table.epos[1]: 
                if self.table.hl != (0, 1, 1, 1):  
                    self.table.hl = (0, 1, 1, 1)
                    redraw = 1
                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.table.press = 1
                        self.table.move = 0
                        return {'RUNNING_MODAL'}
                    elif event.value == 'RELEASE':
                        if not self.table.move:
                            self.table.expand = 0 if self.table.expand else 1
                        self.table.press = 0
                        self.table.move = 0
                        context.area.tag_redraw()
                        return {'RUNNING_MODAL'}
                
                elif event.type == 'ESC':
                    bpy.data.images.remove(self.table.gimage)
                    self.table.plt.close()
                    bpy.types.SpaceView3D.draw_handler_remove(self._handle_en_disp, 'WINDOW')
                    context.area.tag_redraw()
                    return {'CANCELLED'}
                    
                elif self.table.press and event.type == 'MOUSEMOVE':
                     self.table.move = 1
                     self.table.press = 0
            
            elif self.table.hl != (1, 1, 1, 1):
                self.table.hl = (1, 1, 1, 1)
                redraw = 1
                
            if abs(self.dhscatter.lepos[0] - mx) < 20 and abs(self.dhscatter.lspos[1] - my) < 20 and self.dhscatter.expand:
                self.dhscatter.hl = (0, 1, 1, 1) 
                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.dhscatter.resize = 1
                    if self.dhscatter.resize and event.value == 'RELEASE':
                        self.dhscatter.resize = 0
                    return {'RUNNING_MODAL'}

            if abs(self.table.lepos[0] - mx) < 20 and abs(self.table.lspos[1] - my) < 20 and self.table.expand:
                self.table.hl = (0, 1, 1, 1) 
                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.table.resize = 1
                    if self.table.resize and event.value == 'RELEASE':
                        self.table.resize = 0
                    return {'RUNNING_MODAL'}
    
            if event.type == 'MOUSEMOVE':                
                if self.dhscatter.move:
                    self.dhscatter.pos = [mx, my]
                    redraw = 1
                if self.dhscatter.resize:
                    self.dhscatter.lepos[0], self.dhscatter.lspos[1] = mx, my
                    redraw = 1
                if self.table.move:
                    self.table.pos = [mx, my]
                    redraw = 1
                if self.table.resize:
                    self.table.lepos[0], self.table.lspos[1] = mx, my
                    redraw = 1
    
            if self.dhscatter.unit != scene.en_disp_unit or self.dhscatter.cao != context.active_object or \
                self.dhscatter.col != scene.vi_leg_col or self.dhscatter.resstring != retenvires(scene) or \
                self.dhscatter.minmax != envals(scene.en_disp_unit, scene, [0, 100]):
                print(self.dhscatter.unit, scene.en_disp_unit, self.dhscatter.cao, context.active_object, self.dhscatter.col, scene.vi_leg_col, self.dhscatter.resstring, retenvires(scene), self.dhscatter.minmax, envals(scene.en_disp_unit, scene, [0, 100]))
                self.dhscatter.update(context)
                self.table.update(context)

            if redraw:
                context.area.tag_redraw()
                
            return {'PASS_THROUGH'}
        else:
            return {'PASS_THROUGH'}

    def execute(self, context):
        self.i = 0
        scene = context.scene
        scene.en_frame = scene.frame_current
#        (resnode, restree) = scene['viparams']['resnode'].split('@')
        resnode = bpy.data.node_groups[scene['viparams']['resnode'].split('@')[1]].nodes[scene['viparams']['resnode'].split('@')[0]]
        zrl = list(zip(*resnode['reslists']))
        eresobs = {o.name: o.name.upper() for o in bpy.data.objects if o.name.upper() in zrl[2]}
        resstart, resend = 24 * (resnode['Start'] - 1), 24 * (resnode['End']) - 1
        scene.frame_start, scene.frame_end = 0, len(zrl[4][0].split()) - 1
        
        if scene.resas_disp:
            suns = [o for o in bpy.data.objects if o.type == 'LAMP' and o.data.type == 'SUN']
            if not suns:
                bpy.ops.object.lamp_add(type='SUN')
                sun = bpy.context.object
            else:
                sun = suns[0]
            
            for mi, metric in enumerate(zrl[3]):
                if metric == 'Direct Solar (W/m^2)':
                    dirsol = [float(ds) for ds in zrl[4][mi].split()[resstart:resend]]
                elif metric == 'Diffuse Solar (W/m^2)':
                    difsol = [float(ds) for ds in zrl[4][mi].split()[resstart:resend]]
                elif metric == 'Month':
                    mdata = [int(m) for m in zrl[4][mi].split()[resstart:resend]]
                elif metric == 'Day':
                    ddata = [int(d) for d in zrl[4][mi].split()[resstart:resend]]
                elif metric == 'Hour':
                    hdata = [int(h) for h in zrl[4][mi].split()[resstart:resend]]

            sunposenvi(scene, sun, dirsol, difsol, mdata, ddata, hdata)

        if scene.resaa_disp:
            for mi, metric in enumerate(zrl[3]):
                if metric == 'Temperature (degC)' and zrl[1][mi] == 'Climate':
                    temp = [float(ds) for ds in zrl[4][mi].split()[24 * resnode['Start']:24 * resnode['End'] + 1]]
                elif metric == 'Wind Speed (m/s)' and zrl[1][mi] == 'Climate':
                    ws = [float(ds) for ds in zrl[4][mi].split()[24 * resnode['Start']:24 * resnode['End'] + 1]]
                elif metric == 'Wind Direction (deg)' and zrl[1][mi] == 'Climate':
                    wd = [float(m) for m in zrl[4][mi].split()[24 * resnode['Start']:24 * resnode['End'] + 1]]
                elif metric == 'Humidity (%)' and zrl[1][mi] == 'Climate':
                    hu = [float(d) for d in zrl[4][mi].split()[24 * resnode['Start']:24 * resnode['End'] + 1]]
            
            self._handle_air = bpy.types.SpaceView3D.draw_handler_add(en_air, (self, context, temp, ws, wd, hu), 'WINDOW', 'POST_PIXEL')        
        zmetrics = set([zr for zri, zr in enumerate(zrl[3]) if zrl[1][zri] == 'Zone'  and zrl[0][zri] != 'All'])
        
        if scene.reszt_disp and 'Temperature (degC)' in zmetrics:
            envizres(scene, eresobs, resnode, 'Temp')
        if scene.reszsg_disp and  'Solar gain (W)' in zmetrics:
            envizres(scene, eresobs, resnode, 'SHG')
        if scene.reszh_disp and 'Humidity (%)' in zmetrics:
            envizres(scene, eresobs, resnode, 'Hum')
        if scene.reszco_disp and 'CO2 (ppm)' in zmetrics:
            envizres(scene, eresobs, resnode, 'CO2')
        if scene.reszhw_disp and 'Heating (W)' in zmetrics:
            envizres(scene, eresobs, resnode, 'Heat')
        if scene.reszhw_disp and 'Cooling (W)' in zmetrics:
            envizres(scene, eresobs, resnode, 'Cool')
        if scene.reszpmv_disp and 'PMV' in zmetrics:
            envizres(scene, eresobs, resnode, 'PMV')
        if scene.reszppd_disp and 'PPD (%)' in zmetrics:
            envizres(scene, eresobs, resnode, 'PPD')
        if scene.reshrhw_disp and 'HR heating (W)' in zmetrics:
            envizres(scene, eresobs, resnode, 'HRheat')
        if scene.reszof_disp:
            envilres(scene, resnode)
        if scene.reszlf_disp:
            envilres(scene, resnode)

        scene.frame_set(scene.frame_start)
        bpy.app.handlers.frame_change_pre.clear()
        bpy.app.handlers.frame_change_pre.append(recalculate_text)
        self.dhscatter = en_scatter([160, context.region.height - 40], context.region.width, context.region.height, 'scat.png', 600, 400)
        self.dhscatter.update(context)
        self.table = en_table([240, context.region.height - 40], context.region.width, context.region.height, 'table.png', 600, 150)
        self.table.update(context)           
        self._handle_en_disp = bpy.types.SpaceView3D.draw_handler_add(en_disp, (self, context, resnode), 'WINDOW', 'POST_PIXEL')
        scene['viparams']['vidisp'] = 'enpanel'
        scene.vi_display = True
        context.window_manager.modal_handler_add(self)
#        scene.update()
        return {'RUNNING_MODAL'}

class VIEW3D_OT_EnPDisplay(bpy.types.Operator):
    bl_idname = "view3d.enpdisplay"
    bl_label = "EnVi parametric display"
    bl_description = "Display the parametric EnVi results"
    bl_options = {'REGISTER'}
#    bl_undo = False
    _handle = None
    disp =  bpy.props.IntProperty(default = 1)
    
    @classmethod
    def poll(cls, context):
        return context.area.type  == 'VIEW_3D' and \
               context.region.type == 'WINDOW'
    
    def modal(self, context, event):
        redraw = 0
        scene = context.scene

        if event.type != 'INBETWEEN_MOUSEMOVE':   
            if scene.vi_display == 0 or scene['viparams']['vidisp'] != 'enpanel':
                scene['viparams']['vidisp'] = 'en'
                bpy.types.SpaceView3D.draw_handler_remove(self._handle_en_pdisp, 'WINDOW')
                for o in [o for o in scene.objects if o.get('VIType') and o['VIType'] in ('envi_maxtemp', 'envi_maxhum', 'envi_maxheat', 'envi_maxcool', 'envi_maxco2', 'envi_maxshg', 'envi_maxppd', 'envi_maxpmv')]:
                    for oc in o.children:                        
                        [scene.objects.unlink(oc) for oc in o.children]
                        bpy.data.objects.remove(oc)                    
                    scene.objects.unlink(o)
                    bpy.data.objects.remove(o)    
                context.area.tag_redraw()
                return {'CANCELLED'}

            mx, my = event.mouse_region_x, event.mouse_region_y 
            
            if self.barchart.spos[0] < mx < self.barchart.epos[0] and self.barchart.spos[1] < my < self.barchart.epos[1]:
                if self.barchart.hl != (0, 1, 1, 1):
                    self.barchart.hl = (0, 1, 1, 1) 
                    redraw = 1
                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.barchart.press = 1
                        self.barchart.move = 0
                        return {'RUNNING_MODAL'}
                    elif event.value == 'RELEASE':
                        if not self.barchart.move:
                            self.barchart.expand = 0 if self.barchart.expand else 1
                        self.barchart.press = 0
                        self.barchart.move = 0
                        context.area.tag_redraw()
                        return {'RUNNING_MODAL'}
                
                elif event.type == 'ESC':
                    bpy.data.images.remove(self.barchart.gimage)
                    self.barchart.plt.close()
                    bpy.types.SpaceView3D.draw_handler_remove(self._handle_en_disp, 'WINDOW')
                    context.area.tag_redraw()
                    return {'CANCELLED'}
                    
                elif self.barchart.press and event.type == 'MOUSEMOVE':
                     self.barchart.move = 1
                     self.barchart.press = 0
        
            elif self.barchart.lspos[0] < mx < self.barchart.lepos[0] and self.barchart.lspos[1] < my < self.barchart.lepos[1] and abs(self.barchart.lepos[0] - mx) > 20 and abs(self.barchart.lspos[1] - my) > 20:
                if self.barchart.expand: 
                    if self.barchart.hl != (0, 1, 1, 1):
                        self.barchart.hl = (0, 1, 1, 1)
                        redraw = 1
                    if event.type == 'LEFTMOUSE' and event.value == 'PRESS' and self.barchart.expand and self.barchart.lspos[0] < mx < self.barchart.lepos[0] and self.barchart.lspos[1] < my < self.barchart.lspos[1] + 0.9 * self.barchart.ydiff:
                        self.barchart.show_plot()
                        context.area.tag_redraw()
                        return {'RUNNING_MODAL'}
                        
            elif abs(self.barchart.lepos[0] - mx) < 20 and abs(self.barchart.lspos[1] - my) < 20 and self.barchart.expand:
                if self.barchart.hl != (0, 1, 1, 1):
                    self.barchart.hl = (0, 1, 1, 1) 
                    redraw = 1
                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.barchart.resize = 1
                    if self.barchart.resize and event.value == 'RELEASE':
                        self.barchart.resize = 0
                    return {'RUNNING_MODAL'}
                                       
            else:
                if self.barchart.hl != (1, 1, 1, 1):
                    self.barchart.hl = (1, 1, 1, 1)
                    redraw = 1
                
            if self.table.spos[0] < mx < self.table.epos[0] and self.table.spos[1] < my < self.table.epos[1]:
                if self.table.hl != (0, 1, 1, 1):
                    self.table.hl = (0, 1, 1, 1)  
                    redraw = 1
                    
                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.table.press = 1
                        self.table.move = 0
                        return {'RUNNING_MODAL'}
                    elif event.value == 'RELEASE':
                        if not self.table.move:
                            self.table.expand = 0 if self.table.expand else 1
                        self.table.press = 0
                        self.table.move = 0
                        context.area.tag_redraw()
                        return {'RUNNING_MODAL'}
                
                elif event.type == 'ESC':
                    bpy.data.images.remove(self.table.gimage)
                    self.table.plt.close()
                    bpy.types.SpaceView3D.draw_handler_remove(self._handle_en_disp, 'WINDOW')
                    context.area.tag_redraw()
                    return {'CANCELLED'}
                    
                elif self.table.press and event.type == 'MOUSEMOVE':
                     self.table.move = 1
                     self.table.press = 0

            elif abs(self.table.lepos[0] - mx) < 20 and abs(self.table.lspos[1] - my) < 20 and self.table.expand:
                if self.table.hl != (0, 1, 1, 1):
                    self.table.hl = (0, 1, 1, 1)
                    redraw = 1
                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.table.resize = 1
                        return {'RUNNING_MODAL'}
                    if self.table.resize and event.value == 'RELEASE':
                        self.table.resize = 0
                        return {'RUNNING_MODAL'}
                                    
            else:
                if self.table.hl != (1, 1, 1, 1):
                    self.table.hl = (1, 1, 1, 1)
                    redraw = 1
                    
            if event.type == 'MOUSEMOVE':                
                if self.barchart.move:
                    self.barchart.pos = [mx, my]
                    redraw = 1
                           
                if self.barchart.resize:
                    self.barchart.lepos[0], self.barchart.lspos[1] = mx, my
                    redraw = 1
               
                if self.table.move:
                    self.table.pos = [mx, my]
                    redraw = 1
               
                if self.table.resize:
                    self.table.lepos[0], self.table.lspos[1] = mx, my
                    redraw = 1
            
            if self.barchart.unit != scene.en_disp_punit or self.barchart.cao != context.active_object or \
                self.barchart.resstring != retenvires(scene) or self.barchart.col != scene.vi_leg_col or self.barchart.minmax != (scene.bar_min, scene.bar_max):
                self.barchart.update(context)
                self.table.update(context)
                redraw = 1

            if redraw:
                context.area.tag_redraw()
                
            return {'PASS_THROUGH'}
        else:
            return {'PASS_THROUGH'}
                
    def execute(self, context):
        scene = context.scene
        scene.en_frame = scene.frame_current
        resnode = bpy.data.node_groups[scene['viparams']['resnode'].split('@')[1]].nodes[scene['viparams']['resnode'].split('@')[0]]
        zrl = list(zip(*resnode['reslists']))
        eresobs = {o.name: o.name.upper() for o in bpy.data.objects if o.name.upper() in zrl[2]}
        scene.frame_start, scene.frame_end = scene['enparams']['fs'], scene['enparams']['fe']                
        zmetrics = set([zr for zri, zr in enumerate(zrl[3]) if zrl[1][zri] == 'Zone' and zrl[0][zri] == 'All'])

        if scene.resazmaxt_disp and 'Max temp (C)' in zmetrics:
            envizres(scene, eresobs, resnode, 'MaxTemp')
        if scene.resazavet_disp and 'Avg temp (C)' in zmetrics:
            envizres(scene, eresobs, resnode, 'AveTemp')
        if scene.resazmint_disp and 'Min temp (C)' in zmetrics:
            envizres(scene, eresobs, resnode, 'MinTemp')
        if scene.resazmaxhw_disp and 'Max heating (W)' in zmetrics:
            envizres(scene, eresobs, resnode, 'MaxHeat')
        if scene.resazavehw_disp and 'Avg heating (W)' in zmetrics:
            envizres(scene, eresobs, resnode, 'AveHeat')
        if scene.resazminhw_disp and 'Min heating (W)' in zmetrics:
            envizres(scene, eresobs, resnode, 'MinHeat')
        if scene.reszof_disp:
            envilres(scene, resnode)
        if scene.reszlf_disp:
            envilres(scene, resnode)

        scene.frame_set(scene.frame_start)
        bpy.app.handlers.frame_change_pre.clear()
        bpy.app.handlers.frame_change_pre.append(recalculate_text)
        scene['viparams']['vidisp'] = 'enpanel'
        scene.vi_display = True
        context.window_manager.modal_handler_add(self)
        self.barchart = en_barchart([160, context.region.height - 40], context.region.width, context.region.height, 'stats.png', 600, 400)
        self.barchart.update(context)
        self.table = en_table([240, context.region.height - 40], context.region.width, context.region.height, 'table.png', 600, 150)
        self.table.update(context)
        self._handle_en_pdisp = bpy.types.SpaceView3D.draw_handler_add(en_pdisp, (self, context, resnode), 'WINDOW', 'POST_PIXEL')
        return {'RUNNING_MODAL'}

class NODE_OT_Chart(bpy.types.Operator, io_utils.ExportHelper):
    bl_idname = "node.chart"
    bl_label = "Chart"
    bl_description = "Create a 2D graph from the results file"
    bl_register = True
    bl_undo = True
    nodeid = bpy.props.StringProperty()

    def invoke(self, context, event):
        node = bpy.data.node_groups[self.nodeid.split('@')[1]].nodes[self.nodeid.split('@')[0]]
        innodes = list(OrderedDict.fromkeys([inputs.links[0].from_node for inputs in node.inputs if inputs.links]))
        rl = innodes[0]['reslists']
        zrl = list(zip(*rl))
        year = innodes[0]['year']
        
        if node.inputs['X-axis'].framemenu not in zrl[0]:
            self.report({'ERROR'},"There are no results in the results file. Check the results.err file in Blender's text editor")
            return {'CANCELLED'}
        
        if not mp:
            self.report({'ERROR'},"Matplotlib cannot be found by the Python installation used by Blender")
            return {'CANCELLED'}
        plt.clf()
        Sdate = dt.fromordinal(dt(year, 1, 1).toordinal() + node['Start'] - 1)# + datetime.timedelta(hours = node.dsh - 1)
        Edate = dt.fromordinal(dt(year, 1, 1).toordinal() + node['End'] - 1)# + datetime.timedelta(hours = node.deh - 1)
        chart_disp(self, plt, node, innodes, Sdate, Edate)
        return {'FINISHED'}

class NODE_OT_FileProcess(bpy.types.Operator, io_utils.ExportHelper):
    bl_idname = "node.fileprocess"
    bl_label = "Process"
    bl_description = "Process EnergyPlus results file"
    bl_register = True
    bl_undo = True
    nodeid = bpy.props.StringProperty()

    def invoke(self, context, event):
        node = bpy.data.node_groups[self.nodeid.split('@')[1]].nodes[self.nodeid.split('@')[0]]
        self.resname = node.filebase
        processf(self, context.scene, node)
        node.export()
        return {'FINISHED'}

class NODE_OT_SunPath(bpy.types.Operator):
    bl_idname = "node.sunpath"
    bl_label = "Sun Path"
    bl_description = "Create a Sun Path"
    bl_register = True
    bl_undo = True
#    nodeid = bpy.props.StringProperty()

    def invoke(self, context, event):
        scene = context.scene
        if viparams(self, scene):
            self.report({'ERROR'},"Save the Blender file before continuing")
            return {'CANCELLED'}
        
        try:
            spcoll = bpy.data.collections['SunPath']
        except:
            spcoll = bpy.data.collections.new('SunPath')
            bpy.context.scene.collection.children.link(spcoll)
            
        solringnum, sd, numpos = 0, 100, {}
        node = context.node
        node.export()
        scene['viparams']['resnode'], scene['viparams']['restree'] = node.name, node.id_data.name
        scene.cursor.location = (0.0, 0.0, 0.0)
        suns = [ob for ob in scene.objects if ob.type == 'LIGHT' and ob.data.type == 'SUN']
        sunmeshes = [sunmesh for sunmesh in scene.objects if sunmesh.get('VIType') == "SunMesh"]

        for sm in sunmeshes:
            bpy.data.objects.remove(sm, do_unlink=True, do_id_user=True, do_ui_user=True)
#            spcoll.objects.unlink(sm)
#            delobj(context.view_layer, sm)
            
        print([sunmesh for sunmesh in bpy.data.objects if sunmesh.get('VIType') == "SunMesh"])

        requiredsuns = {'0': 1, '1': 12, '2': 24}[node.suns]

        matdict = {'SolEquoRings': (1, 0, 0, 1), 'HourRings': (1, 1, 0, 1), 'SPBase': (1, 1, 1, 1), 'Sun': (1, 1, 1, 1), 'PathDash': (1, 1, 1, 1),
                   'SumAng': (1, 0, 0, 1), 'EquAng': (0, 1, 0, 1), 'WinAng': (0, 0, 1, 1)}
        
        for mat in [mat for mat in matdict if mat not in bpy.data.materials]:
            bpy.data.materials.new(mat)
            bpy.data.materials[mat].diffuse_color = matdict[mat]
#            bpy.data.materials[mat].use_shadeless = 1
            bpy.data.materials[mat].use_nodes = True
            nodes = bpy.data.materials[mat].node_tree.nodes

            for n in nodes:
                nodes.remove(n)
            
            if mat == 'PathDash':
                bpy.data.materials[mat].diffuse_color[3] = 0
                node_material = nodes.new(type='ShaderNodeBsdfTransparent')
            else:
                node_material = nodes.new(type='ShaderNodeEmission')
                node_material.inputs[1].default_value = 1.0

            node_material.inputs[0].default_value = matdict[mat]
            node_material.location = 0,0
            node_output = nodes.new(type='ShaderNodeOutputMaterial')   
            node_output.location = 400,0            
            links = bpy.data.materials[mat].node_tree.links
            links.new(node_material.outputs[0], node_output.inputs[0])
                            
        if suns:
            for sun in suns[requiredsuns:]: 
                bpy.data.objects.remove(sun, do_unlink=True, do_id_user=True, do_ui_user=True)
#                spcoll.objects.unlink(sun)
#                delobj(context.view_layer, sun)
#            [bpy.data.objects.remove(sun) for sun in suns[requiredsuns:]]
            suns = [ob for ob in context.scene.objects if ob.type == 'LIGHT' and ob.data.type == 'SUN']            
            [sun.animation_data_clear() for sun in suns]

        if not suns or len(suns) < requiredsuns: 
            for rs in range(requiredsuns - len(suns)):
                bpy.ops.object.light_add(type='SUN', radius=1, view_align=False, location=(0, 0, 0))

#                bpy.ops.object.lamp_add(type = "SUN")
                suns.append(context.active_object)
       
        if scene.render.engine == 'CYCLES' and scene.world.get('node_tree') and 'Sky Texture' in [no.bl_label for no in scene.world.node_tree.nodes]:
            scene.world.node_tree.animation_data_clear()    
        
        if bpy.context.active_object and not bpy.context.active_object.hide_viewport:
            if bpy.context.active_object.type == 'MESH':
                bpy.ops.object.mode_set(mode = 'OBJECT')
        
        for ob in context.scene.objects:
            if ob.get('VIType') == "SPathMesh": 
#                selobj(context.view_layer, ob)
                bpy.data.objects.remove(ob, do_unlink=True, do_id_user=True, do_ui_user=True)
#                spcoll.objects.unlink(ob)
#                delobj(context.view_layer, ob)
#                bpy.ops.object.delete(use_global=True)

#                context.scene.objects.unlink(ob)
#                ob.name = 'oldspathmesh'

#        if "SkyMesh" not in [ob.get('VIType') for ob in context.scene.objects]:
#            bpy.data.materials.new('SkyMesh')
#            bpy.ops.mesh.primitive_uv_sphere_add(segments=32, ring_count=16, radius=52.5)
#            smesh = context.active_object
#            smesh.location, smesh.rotation_euler[0], smesh.cycles_visibility.shadow, smesh.name, smesh['VIType']  = (0,0,0), pi, False, "SkyMesh", "SkyMesh"
#            bpy.ops.object.material_slot_add()
#            smesh.material_slots[0].material = bpy.data.materials['SkyMesh']
#            bpy.ops.object.shade_smooth()
#            smesh.hide_viewport, smesh.hide_render = True, True
#        else:
#            smesh =  [ob for ob in context.scene.objects if ob.get('VIType') and ob['VIType'] == "SkyMesh"][0]
          
            
            
        bpy.ops.object.add(type = "MESH")
        spathob = context.active_object
        if spathob.name not in spcoll.objects:
            spcoll.objects.link(spathob)
            scene.collection.objects.unlink(spathob)
#        scene.collection.objects.unlink(spathob)
        spathob.location, spathob.name,  spathob['VIType'], spathmesh = (0, 0, 0), "SPathMesh", "SPathMesh", spathob.data
#        smesh.parent = spathob
        
        for s, sun in enumerate(suns):
            if sun.name not in spcoll.objects:
                spcoll.objects.link(sun)
                scene.collection.objects.unlink(sun)
                
            sun.data.shadow_soft_size = 0.01            
            sun['VIType'] = 'Sun'
            sun['solhour'], sun['solday'] = scene.solhour, scene.solday
            sun.name = sun.data.name ='Sun{}'.format(s)
            bpy.ops.mesh.primitive_uv_sphere_add(segments=12, ring_count=12, radius=0.5)
            sunob = context.active_object
            
            if sunob.name not in spcoll.objects:
                spcoll.objects.link(sunob)
                scene.collection.objects.unlink(sunob)
#            scene.collection.objects.unlink(sunob)
            sunob.location, sunob.cycles_visibility.shadow, sunob.name, sunob['VIType'] = (0, 0, 0), 0, "SunMesh{}".format(s), "SunMesh"
            sunob.cycles_visibility.diffuse, sunob.cycles_visibility.shadow, sunob.cycles_visibility.glossy, sunob.cycles_visibility.transmission, sunob.cycles_visibility.scatter = [False] * 5

            if len(sunob.material_slots) == 0:
                 bpy.ops.object.material_slot_add()
                 sunob.material_slots[0].material = bpy.data.materials['Sun']
                 
            sun.parent = spathob
            sunob.parent = sun
        
        bm = bmesh.new()
        bm.from_mesh(spathmesh)

        for doy in range(0, 365, 2):
            for hour in range(1, 25):
                ([solalt, solazi]) = solarPosition(doy, hour, scene.latitude, scene.longitude)[2:]
                bm.verts.new().co = [-(sd-(sd-(sd*cos(solalt))))*sin(solazi), -(sd-(sd-(sd*cos(solalt))))*cos(solazi), sd*sin(solalt)]
        
        if hasattr(bm.verts, "ensure_lookup_table"):
            bm.verts.ensure_lookup_table()
        for v in range(24, len(bm.verts)):
            bm.edges.new((bm.verts[v], bm.verts[v - 24]))
        if v in range(8568, 8761):
            bm.edges.new((bm.verts[v], bm.verts[v - 8568]))

        for doy in (79, 172, 355):
            for hour in range(1, 241):
                ([solalt, solazi]) = solarPosition(doy, hour*0.1, scene.latitude, scene.longitude)[2:]
                vcoord = [-(sd-(sd-(sd*cos(solalt))))*sin(solazi), -(sd-(sd-(sd*cos(solalt))))*cos(solazi), sd*sin(solalt)]
                bm.verts.new().co = vcoord
                if hasattr(bm.verts, "ensure_lookup_table"):
                    bm.verts.ensure_lookup_table()
                if bm.verts[-1].co.z >= 0 and doy in (172, 355) and not hour%10:
                    numpos['{}-{}'.format(doy, int(hour*0.1))] = vcoord
                if hour != 1:
                    bm.edges.new((bm.verts[-2], bm.verts[-1]))
                    solringnum += 1
                if hour == 240:
                    bm.edges.new((bm.verts[-240], bm.verts[-1]))
                    solringnum += 1
        
        bm.to_mesh(spathmesh)
        bm.free()
        selobj(context.view_layer, spathob)
        bpy.ops.object.convert(target='CURVE')
        spathob.data.bevel_depth, spathob.data.bevel_resolution = node.th, node.res
        bpy.context.object.data.fill_mode = 'FULL'
        bpy.ops.object.convert(target='MESH')
        
        bpy.ops.object.material_slot_add()
        spathob.material_slots[0].material, spathob['numpos'] = bpy.data.materials['HourRings'], numpos
        bpy.ops.object.material_slot_add()
        spathob.material_slots[1].material = bpy.data.materials['PathDash']

        for face in spathob.data.polygons:
            face.material_index = 0 if not int(face.index/16)%2 else 1
                
        for vert in spathob.data.vertices[0:16 * (solringnum + 3)]:
            vert.select = True

        bpy.ops.object.material_slot_add()
        spathob.material_slots[-1].material = bpy.data.materials['SolEquoRings']
        spathob.active_material_index = 2
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.select_mode(type="VERT")
        bpy.ops.object.material_slot_assign()
        bpy.ops.mesh.select_all(action='SELECT')
        bpy.ops.mesh.bisect(plane_co=(0.0, 0.0, 0.0), plane_no=(0.0, 0.0, 1.0), use_fill=True, clear_inner=True, clear_outer=False)
        bpy.ops.object.mode_set(mode='OBJECT')
        bpy.ops.object.select_all(action='DESELECT')
        compassos = compass((0,0,0.01), sd, spathob, bpy.data.materials['SPBase'])
        spro = spathrange([bpy.data.materials['SumAng'], bpy.data.materials['EquAng'], bpy.data.materials['WinAng']])
        joinobj(context.view_layer, [compassos] + [spro] + [spathob])

#        for ob in (spathob, smesh):
        spathob.cycles_visibility.diffuse, spathob.cycles_visibility.shadow, spathob.cycles_visibility.glossy, spathob.cycles_visibility.transmission, spathob.cycles_visibility.scatter = [False] * 5
        spathob.show_transparent = True

        if spfc not in bpy.app.handlers.frame_change_post:
            bpy.app.handlers.frame_change_post.append(spfc)

        scene['viparams']['vidisp'] = 'sp'
        scene['spparams']['suns'] = node.suns
        context.scene['viparams']['visimcontext'] = 'SunPath'
        bpy.ops.view3d.spnumdisplay('INVOKE_DEFAULT')
        sunpath(scene)
        return {'FINISHED'}

class VIEW3D_OT_SPNumDisplayold(bpy.types.Operator):
    '''Display results legend and stats in the 3D View'''
    bl_idname = "view3d.spnumdisplay"
    bl_label = "Point numbers"
    bl_description = "Display the times and solstices on the sunpath"
    bl_register = True
    bl_undo = False

    def modal(self, context, event):
        scene = context.scene
        if context.area:
            context.area.tag_redraw()
        if scene.vi_display == 0 or scene['viparams']['vidisp'] != 'sp':
            bpy.types.SpaceView3D.draw_handler_remove(self._handle_spnum, 'WINDOW')
            [bpy.data.objects.remove(o, do_unlink=True, do_id_user=True, do_ui_user=True) for o in bpy.data.objects if o.get('VIType') and o['VIType'] in ('SunMesh', 'SkyMesh')]
            return {'CANCELLED'}
        return {'PASS_THROUGH'}

    def invoke(self, context, event):
        scene = context.scene
        simnode = bpy.data.node_groups[scene['viparams']['restree']].nodes[scene['viparams']['resnode']]
        
        if simnode.suns != '0':
            scene.timedisp = 0
            
        self._handle_spnum = bpy.types.SpaceView3D.draw_handler_add(spnumdisplay, (self, context, simnode), 'WINDOW', 'POST_PIXEL')
        context.window_manager.modal_handler_add(self)
        scene.vi_display = 1
        return {'RUNNING_MODAL'}
    
class VIEW3D_OT_SPNumDisplay(bpy.types.Operator):
    '''Display results legend and stats in the 3D View'''
    bl_idname = "view3d.spnumdisplay"
    bl_label = "Point numbers"
    bl_description = "Display the times and solstices on the sunpath"
    bl_register = True
    bl_undo = False
    
   
    def create_batch(self, scene, node):
        vertex_shader = '''
            uniform mat4 viewProjectionMatrix;
            uniform vec4 color1;
            uniform vec4 color2;
            in vec3 position;
            in float arcLength;
            in uint line_break;
            
            out vec4 v_color1;
            out vec4 v_color2;
            out float v_ArcLength;
            out float zpos;
            flat out uint lb;
            
            void main()
            {
                v_color1 = color1;
                v_color2 = color2;
                v_ArcLength = arcLength;
                gl_Position = viewProjectionMatrix * vec4(position, 1.0f);
                zpos = vec3(position)[2];
                lb = line_break;
            }
        '''
#        vertex_shader = '''
#            uniform mat4 viewProjectionMatrix;
#            uniform vec4 color1;
#            uniform vec4 color2;
#            in vec3 position;
#            in float arcLength;
#            
#            out vec4 v_color1;
#            out vec4 v_color2;
#            out float v_ArcLength;
#            void main()
#            {
#                v_color1 = color1;
#                v_color2 = color2;
#                v_ArcLength = arcLength;
#                gl_Position = viewProjectionMatrix * vec4(position, 1.0f);
#            }
#        '''
        
#        fragment_shader = '''
#            uniform float u_Scale;
#            in vec4 v_color1;
#            in vec4 v_color2;
#            in float v_ArcLength;
#            void main()
#            {
#                if (step(sin(v_ArcLength * u_Scale), -0.5) == 1) {gl_FragColor = v_color1;} else {gl_FragColor = v_color2;};
#            }
#        '''
        fragment_shader = '''
            uniform float u_Scale;
            in float zpos;
            in vec4 v_color1;
            in vec4 v_color2;
            in float v_ArcLength;
            flat in uint lb;
 
            void main()
            {
                if (zpos < 0) {discard;}; 
                if (lb == uint(1)) {discard;}; 
                if (step(sin(v_ArcLength * u_Scale), -0.5) == 1) {gl_FragColor = v_color1;} else {gl_FragColor = v_color2;};
            }
        '''
        breaks, coords, sd, d = [], [], 100, 0
        
        
        for hour in range(1, 25):
            for doy in range(0, 365, 2):

                ([solalt, solazi]) = solarPosition(doy, hour, scene.latitude, scene.longitude)[2:]
                coords.append(Vector([-(sd-(sd-(sd*cos(solalt))))*sin(solazi), -(sd-(sd-(sd*cos(solalt))))*cos(solazi), sd*sin(solalt)]))
                if d%183 == 0:
                    breaks.append(1)
                else:
                    breaks.append(0)
                d += 1
#        print(breaks, len(breaks), len(coords))    
        line_lengths = [0]
        
        for a, b in zip(coords[:-1], coords[1:]):
            line_lengths.append(line_lengths[-1] + (a - b).length)        
        
        self.shader = gpu.types.GPUShader(vertex_shader, fragment_shader)        
        self.batch = batch_for_shader(self.shader, 'LINE_STRIP', {"position": coords, "arcLength": line_lengths, "line_break": breaks})
        print('hello')
        
    def draw_sp(self, op, context):
        # Draw lines
        bgl.glEnable(bgl.GL_BLEND)
        bgl.glEnable(bgl.GL_LINE_SMOOTH)
        bgl.glLineWidth(context.scene.vi_params.sp_line_width)
        self.shader.bind()
        matrix = bpy.context.region_data.perspective_matrix
#        print(dir(self.shader))
#        self.shader.uniform_vector_float("position", coords)
        self.shader.uniform_float("viewProjectionMatrix", matrix)
        self.shader.uniform_float("color1", context.scene.vi_params.sp_hour_main)
        self.shader.uniform_float("color2", context.scene.vi_params.sp_hour_dash)
        self.shader.uniform_float("u_Scale", 1)
        self.batch.draw(self.shader)
        
        bgl.glDisable(bgl.GL_LINE_SMOOTH)
        bgl.glDisable(bgl.GL_BLEND)

    def modal(self, context, event):
        scene = context.scene
        
        
        if context.area:
            context.area.tag_redraw()
        if scene.vi_display == 0 or scene['viparams']['vidisp'] != 'sp':
            bpy.types.SpaceView3D.draw_handler_remove(self._handle_spnum, 'WINDOW')
            [bpy.data.objects.remove(o, do_unlink=True, do_id_user=True, do_ui_user=True) for o in bpy.data.objects if o.get('VIType') and o['VIType'] in ('SunMesh', 'SkyMesh')]
            return {'CANCELLED'}
        return {'PASS_THROUGH'}

    def invoke(self, context, event):
        
        scene = context.scene
        node = context.node
        self.create_batch(scene, node)
        
        self.draw_handle_spnum = bpy.types.SpaceView3D.draw_handler_add(self.draw_sp, (self, context), "WINDOW", "POST_VIEW")

#        self.draw_handle_2d = bpy.types.SpaceView3D.draw_handler_add(
#            self.draw_callback_2d, args, "WINDOW", "POST_PIXEL")

#        self.draw_event = context.window_manager.event_timer_add(0.1, window=context.window)
        
#        simnode = bpy.data.node_groups[scene['viparams']['restree']].nodes[scene['viparams']['resnode']]
        
#        if simnode.suns != '0':
#            scene.timedisp = 0
            
#        self._handle_spnum = bpy.types.SpaceView3D.draw_handler_add(self.draw, (self, context, node), 'WINDOW', 'POST_PIXEL')
        context.window_manager.modal_handler_add(self)
        scene.vi_display = 1
        return {'RUNNING_MODAL'}
        
class NODE_OT_WindRose(bpy.types.Operator):
    bl_idname = "node.windrose"
    bl_label = "Wind Rose"
    bl_description = "Create a Wind Rose"
    bl_register = True
    bl_undo = True
    nodeid = bpy.props.StringProperty()

    def invoke(self, context, event):
        scene = context.scene
        simnode = bpy.data.node_groups[self.nodeid.split('@')[1]].nodes[self.nodeid.split('@')[0]]
        plt = ret_plt()
        
        if viparams(self, scene):
            return {'CANCELLED'}
        if not plt:
            self.report({'ERROR'},"There is something wrong with your matplotlib installation")
            return {'FINISHED'}

        simnode.export()
        locnode = simnode.inputs['Location in'].links[0].from_node
        scene['viparams']['resnode'], scene['viparams']['restree'] = simnode.name, self.nodeid.split('@')[1]
        scene['viparams']['vidisp'], scene.vi_display = 'wr', 1
        context.scene['viparams']['visimcontext'] = 'Wind'
        rl = locnode['reslists']
        cdoys = [float(c) for c in [r[4].split() for r in rl if r[0] == '0' and r[1] == 'Time' and r[2] == '' and r[3] == 'DOS'][0]]
        cwd = [float(c) for c in [r[4].split() for r in rl if r[0] == '0' and r[1] == 'Climate' and r[2] == '' and r[3] == 'Wind Direction (deg)'][0]]
        cws = [float(c) for c in [r[4].split() for r in rl if r[0] == '0' and r[1] == 'Climate' and r[2] == '' and r[3] == 'Wind Speed (m/s)'][0]]        
        doys = list(range(simnode.sdoy, simnode.edoy + 1)) if simnode.edoy > simnode.sdoy else list(range(1, simnode.edoy + 1)) + list(range(simnode.sdoy, 366))
        awd = array([wd for di, wd in enumerate(cwd) if cdoys[di] in doys])
        aws = array([ws for di, ws in enumerate(cws) if cdoys[di] in doys])
        validdata = numpy.where(awd > 0) if max(cwd) == 360 else numpy.where(awd > -1)
        vawd = awd[validdata]
        vaws = aws[validdata]
        simnode['maxres'], simnode['minres'], simnode['avres'] = max(cws), min(cws), sum(cws)/len(cws)
        fig = plt.figure(figsize=(8, 8), dpi=150, facecolor='w', edgecolor='w')
        rect = [0.1, 0.1, 0.8, 0.8]
        ax = WindroseAxes(fig, rect, facecolor='w')
        fig.add_axes(ax)
#        (fig, ax) = wr_axes(plt)
        sbinvals = arange(0,int(ceil(max(cws))),2)
        dbinvals = arange(-11.25,372.25,22.5)
        dfreq = histogram(cwd, bins=dbinvals)[0]
        adfreq = histogram(vawd, bins=dbinvals)[0]
        dfreq[0] = dfreq[0] + dfreq[-1]
        dfreq = dfreq[:-1]
        
        if simnode.wrtype == '0':
            ax.bar(vawd, vaws, bins=sbinvals, normed=True, opening=0.8, edgecolor='white', cmap=mcm.get_cmap(scene.vi_leg_col))
        elif simnode.wrtype == '1':
            ax.box(vawd, vaws, bins=sbinvals, normed=True, cmap=mcm.get_cmap(scene.vi_leg_col))
        elif simnode.wrtype in ('2', '3', '4'):
            ax.contourf(vawd, vaws, bins=sbinvals, normed=True, cmap=mcm.get_cmap(scene.vi_leg_col))
        
        if simnode.max_freq == '1':
            ax.set_rmax(simnode.max_freq_val)
            
        plt.savefig(scene['viparams']['newdir']+'/disp_wind.svg')
        (wro, scale) = wind_rose(simnode['maxres'], scene['viparams']['newdir']+'/disp_wind.svg', simnode.wrtype, mcolors)
        
        wro['maxres'], wro['minres'], wro['avres'], wro['nbins'], wro['VIType'] = max(aws), min(aws), sum(aws)/len(aws), len(sbinvals), 'Wind_Plane'
        simnode['maxfreq'] = 100*numpy.max(adfreq)/len(vawd) if simnode.max_freq == '0' else simnode.max_freq_val
        windnum(simnode['maxfreq'], (0,0,0), scale, compass((0,0,0), scale, wro, wro.data.materials['wr-000000']))
        plt.close()
        wro['table'] = array([["", 'Minimum', 'Average', 'Maximum'], ['Speed (m/s)', wro['minres'], '{:.1f}'.format(wro['avres']), wro['maxres']], ['Direction (\u00B0)', min(awd), '{:.1f}'.format(sum(awd)/len(awd)), max(awd)]])
        wro['ws'] = aws.reshape(len(doys), 24).T.tolist()
        wro['wd'] = awd.reshape(len(doys), 24).T.tolist()
        wro['days'] = array(doys, dtype = float)
        wro['hours'] = arange(1, 25, dtype = float)        
        wro['maxfreq'] = 100*numpy.max(adfreq)/len(vawd)
        simnode['nbins'] = len(sbinvals)        
        simnode['ws'] = array(cws).reshape(365, 24).T.tolist()
        simnode['wd'] = array(cwd).reshape(365, 24).T.tolist()        
        simnode['days'] = arange(1, 366, dtype = float)
        simnode['hours'] = arange(1, 25, dtype = float)        
        return {'FINISHED'}

class VIEW3D_OT_WRDisplay(bpy.types.Operator):
    '''Display results legend and stats in the 3D View'''
    bl_idname = "view3d.wrdisplay"
    bl_label = "Wind rose display"
    bl_description = "Display wind metrics"
    bl_register = True
    bl_undo = False

    def modal(self, context, event): 
        if context.scene.vi_display == 0 or context.scene['viparams']['vidisp'] != 'wrpanel' or 'Wind_Plane' not in [o['VIType'] for o in bpy.data.objects if o.get('VIType')]:
                bpy.types.SpaceView3D.draw_handler_remove(self._handle_wr_disp, 'WINDOW')
                context.area.tag_redraw()
                return {'CANCELLED'}           

        if event.type != 'INBETWEEN_MOUSEMOVE' and context.region and context.area.type == 'VIEW_3D' and context.region.type == 'WINDOW':            
            mx, my = event.mouse_region_x, event.mouse_region_y 
            
            # Legend routine 
            
            if self.legend.spos[0] < mx < self.legend.epos[0] and self.legend.spos[1] < my < self.legend.epos[1]:
                self.legend.hl = (0, 1, 1, 1)  
                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.legend.press = 1
                        self.legend.move = 0
                        return {'RUNNING_MODAL'}
                    elif event.value == 'RELEASE':
                        if not self.legend.move:
                            self.legend.expand = 0 if self.legend.expand else 1
                        self.legend.press = 0
                        self.legend.move = 0
                        context.area.tag_redraw()
                        return {'RUNNING_MODAL'}
                
                elif event.type == 'ESC':
                    bpy.types.SpaceView3D.draw_handler_remove(self._handle_wr_disp, 'WINDOW')
                    context.area.tag_redraw()
                    return {'CANCELLED'}
                    
                elif self.legend.press and event.type == 'MOUSEMOVE':
                     self.legend.move = 1
                     self.legend.press = 0
            
            elif abs(self.legend.lepos[0] - mx) < 10 and abs(self.legend.lspos[1] - my) < 10:
                self.legend.hl = (0, 1, 1, 1) 
                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.legend.resize = 1
                    if self.legend.resize and event.value == 'RELEASE':
                        self.legend.resize = 0
                    return {'RUNNING_MODAL'}
                    
            else:
                self.legend.hl = (1, 1, 1, 1)
                
            # Scatter routine
                
            if self.dhscatter.spos[0] < mx < self.dhscatter.epos[0] and self.dhscatter.spos[1] < my < self.dhscatter.epos[1]:
                self.dhscatter.hl = (0, 1, 1, 1)  
                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.dhscatter.press = 1
                        self.dhscatter.move = 0
                        return {'RUNNING_MODAL'}
                    elif event.value == 'RELEASE':
                        if not self.dhscatter.move:
                            self.dhscatter.expand = 0 if self.dhscatter.expand else 1
                        self.dhscatter.press = 0
                        self.dhscatter.move = 0
                        context.area.tag_redraw()
                        return {'RUNNING_MODAL'}
                
                elif event.type == 'ESC':
                    bpy.data.images.remove(self.dhscatter.gimage)
                    self.dhscatter.plt.close()
                    bpy.types.SpaceView3D.draw_handler_remove(self._handle_wr_disp, 'WINDOW')
                    context.area.tag_redraw()
                    return {'CANCELLED'}
                    
                elif self.dhscatter.press and event.type == 'MOUSEMOVE':
                     self.dhscatter.move = 1
                     self.dhscatter.press = 0
        
            elif self.dhscatter.lspos[0] < mx < self.dhscatter.lepos[0] and self.dhscatter.lspos[1] < my < self.dhscatter.lepos[1] and abs(self.dhscatter.lepos[0] - mx) > 20 and abs(self.dhscatter.lspos[1] - my) > 20:
                if self.dhscatter.expand: 
                    self.dhscatter.hl = (1, 1, 1, 1)
                    if event.type == 'LEFTMOUSE' and event.value == 'PRESS' and self.dhscatter.expand and self.dhscatter.lspos[0] < mx < self.dhscatter.lepos[0] and self.dhscatter.lspos[1] < my < self.dhscatter.lspos[1] + 0.9 * self.dhscatter.ydiff:
                        self.dhscatter.show_plot()
                    context.area.tag_redraw()
                    return {'RUNNING_MODAL'}
                   
            else:
                self.dhscatter.hl = (1, 1, 1, 1)
                                    
            # Table routine
                     
            if self.table.spos[0] < mx < self.table.epos[0] and self.table.spos[1] < my < self.table.epos[1]:
                self.table.hl = (0, 1, 1, 1)  
                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.table.press = 1
                        self.table.move = 0
                        return {'RUNNING_MODAL'}
                    elif event.value == 'RELEASE':
                        if not self.table.move:
                            self.table.expand = 0 if self.table.expand else 1
                        self.table.press = 0
                        self.table.move = 0
                        context.area.tag_redraw()
                        return {'RUNNING_MODAL'}
                
                elif event.type == 'ESC':
                    bpy.types.SpaceView3D.draw_handler_remove(self._handle_wr_disp, 'WINDOW')
                    context.area.tag_redraw()
                    return {'CANCELLED'}
                    
                elif self.table.press and event.type == 'MOUSEMOVE':
                     self.table.move = 1
                     self.table.press = 0
                     
            else:
                self.table.hl = (1, 1, 1, 1)
                     
            # Resize routines
            
            if abs(self.legend.lepos[0] - mx) < 20 and abs(self.legend.lspos[1] - my) < 20:
                self.legend.hl = (0, 1, 1, 1) 
                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.legend.resize = 1
                    if self.legend.resize and event.value == 'RELEASE':
                        self.legend.resize = 0
                    return {'RUNNING_MODAL'}
                    
            elif abs(self.dhscatter.lepos[0] - mx) < 20 and abs(self.dhscatter.lspos[1] - my) < 20:
                self.dhscatter.hl = (0, 1, 1, 1) 
                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.dhscatter.resize = 1
                    if self.dhscatter.resize and event.value == 'RELEASE':
                        self.dhscatter.resize = 0
                    return {'RUNNING_MODAL'}
                    
            elif abs(self.table.lepos[0] - mx) < 20 and abs(self.table.lspos[1] - my) < 20:
                self.table.hl = (0, 1, 1, 1) 
                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.table.resize = 1
                    if self.table.resize and event.value == 'RELEASE':
                        self.table.resize = 0
                    return {'RUNNING_MODAL'}
            
            # Move routines
                     
            if event.type == 'MOUSEMOVE':                
                if self.legend.move:
                    self.legend.pos = [mx, my]
                if self.legend.resize:
                    self.legend.lepos[0], self.legend.lspos[1] = mx, my
                if self.dhscatter.move:
                    self.dhscatter.pos = [mx, my]
                if self.dhscatter.resize:
                    self.dhscatter.lepos[0], self.dhscatter.lspos[1] = mx, my
                if self.table.move:
                    self.table.pos = [mx, my]
                if self.table.resize:
                    self.table.lepos[0], self.table.lspos[1] = mx, my
                                                
        # Object update routines 
        
            if self.legend.cao != context.active_object:
                self.legend.update(context)
            
            if self.dhscatter.cao != context.active_object or self.dhscatter.unit != context.scene.wind_type or context.scene.vi_leg_col != self.dhscatter.col:
                self.dhscatter.update(context)
                
            if self.table.cao != context.active_object:
                self.table.update(context)
            
            context.area.tag_redraw()
        return {'PASS_THROUGH'}

    def invoke(self, context, event):
        context.scene.vi_display = 1
        context.scene['viparams']['vidisp'] = 'wrpanel'
        simnode = bpy.data.node_groups[context.scene['viparams']['restree']].nodes[context.scene['viparams']['resnode']]
        self.legend = wr_legend([80, context.region.height - 40], context.region.width, context.region.height, 'legend.png', 150, 350)
        self.dhscatter = wr_scatter([160, context.region.height - 40], context.region.width, context.region.height, 'scat.png', 600, 400)
        self.table = wr_table([240, context.region.height - 40], context.region.width, context.region.height, 'table.png', 600, 150)       
        self.legend.update(context)
        self.dhscatter.update(context)
        self.table.update(context)
        self._handle_wr_disp = bpy.types.SpaceView3D.draw_handler_add(wr_disp, (self, context, simnode), 'WINDOW', 'POST_PIXEL')
        context.window_manager.modal_handler_add(self)
        return {'RUNNING_MODAL'}
    
class NODE_OT_SVF(bpy.types.Operator):
    bl_idname = "node.svf"
    bl_label = "Sky View Factor"
    bl_description = "Undertake a sky view factor study"
    bl_register = True
    bl_undo = False
    nodeid = bpy.props.StringProperty()

    def invoke(self, context, event):
        scene = context.scene  
        scene.vi_display = 0
        
        if viparams(self, scene):            
            return {'CANCELLED'}

        shadobs = retobjs('livig')
        if not shadobs:
            self.report({'ERROR'},"No shading objects have a material attached.")
            return {'CANCELLED'}
            
        scene['liparams']['shadc'] = [ob.name for ob in retobjs('ssc')]
        if not scene['liparams']['shadc']:
            self.report({'ERROR'},"No objects have a light sensor material attached.")
            return {'CANCELLED'}

        scene['viparams']['restree'] = self.nodeid.split('@')[1]
        clearscene(scene, self)
        simnode = bpy.data.node_groups[self.nodeid.split('@')[1]].nodes[self.nodeid.split('@')[0]]
        scene['viparams']['visimcontext'] = 'SVF'

        if not scene.get('liparams'):
           scene['liparams'] = {}
           
        scene['liparams']['cp'], scene['liparams']['unit'], scene['liparams']['type'] = simnode.cpoint, '% Sunlit', 'VI Shadow'
        simnode.preexport()
        (scene['liparams']['fs'], scene['liparams']['fe']) = (scene.frame_current, scene.frame_current) if simnode.animmenu == 'Static' else (simnode.startframe, simnode.endframe)      
        scene['viparams']['resnode'], simnode['Animation'] = simnode.name, simnode.animmenu
        (scmaxres, scminres, scavres) = [[x] * (scene['liparams']['fe'] - scene['liparams']['fs'] + 1) for x in (0, 1, 0)]        
        frange = range(scene['liparams']['fs'], scene['liparams']['fe'] + 1)
        x, y, z = [], [], []
        
        if simnode.skypatches == '0':
            alts = (6, 18, 30, 42, 54, 66, 78, 90)
            azis = (30, 30, 24, 24, 18, 12, 6, 1)

        elif simnode.skypatches == '1':
            alts = [(rrow+0.5)*90/(2*7+0.5) for rrow in range(0, 15)]
            azis = (60, 60, 60, 60, 48, 48, 48, 48, 36, 36, 24, 24, 12, 12, 1)
        
        elif simnode.skypatches == '2':
            alts = [(rrow+0.5)*90/(4*7+0.5) for rrow in range(0, 29)]
            azis = (120, 120, 120, 120, 120, 120, 120, 120, 96, 96, 96, 96, 96, 96, 96, 96, 72, 72, 72, 72, 48, 48, 48, 48, 24, 24, 24, 24, 1)

        for a, azi in enumerate(azis):
            for az in arange(0, 360, 360/azi):
                x.append(sin(az * pi/180) * cos(alts[a] * pi/180))
                y.append(cos(az * pi/180) * cos(alts[a] * pi/180))
                z.append(sin(alts[a] * pi/180))   
                    
        valdirecs = [v for v in zip(x, y, z)]
        lvaldirecs = len(valdirecs)
        calcsteps = len(frange) * sum(len([f for f in o.data.polygons if o.data.materials[f.material_index].mattype == '1']) for o in [scene.objects[on] for on in scene['liparams']['shadc']])
        curres, reslists = 0, []
        pfile = progressfile(scene['viparams']['newdir'], datetime.datetime.now(), calcsteps)
        kivyrun = progressbar(os.path.join(scene['viparams']['newdir'], 'viprogress'), 'Sky View')
        
        for oi, o in enumerate([scene.objects[on] for on in scene['liparams']['shadc']]):
            for k in o.keys():
                del o[k]
                
            if any([s < 0 for s in o.scale]):
                logentry('Negative scaling on calculation object {}. Results may not be as expected'.format(o.name))
                self.report({'WARNING'}, 'Negative scaling on calculation object {}. Results may not be as expected'.format(o.name))

            o['omin'], o['omax'], o['oave'] = {}, {}, {}
            bm = bmesh.new()
            bm.from_mesh(o.data)
            clearlayers(bm, 'a')
            bm.transform(o.matrix_world)
            geom = bm.faces if simnode.cpoint == '0' else bm.verts
            geom.layers.int.new('cindex')
            cindex = geom.layers.int['cindex']
            [geom.layers.float.new('res{}'.format(fi)) for fi in frange]
            avres, minres, maxres, g = [], [], [], 0
            
            if simnode.cpoint == '0':
                gpoints = [f for f in geom if o.data.materials[f.material_index].mattype == '1']
            elif simnode.cpoint == '1':
                gpoints = [v for v in geom if any([o.data.materials[f.material_index].mattype == '1' for f in v.link_faces])]

            for g, gp in enumerate(gpoints):
                gp[cindex] = g + 1

            for frame in frange: 
                g, oshadres = 0, array([])                
                scene.frame_set(frame)
                shadtree = rettree(scene, shadobs, ('', '2')[simnode.signore])
                shadres = geom.layers.float['res{}'.format(frame)]
                                    
                if gpoints:
                    posis = [gp.calc_center_bounds() + gp.normal.normalized() * simnode.offset for gp in gpoints] if simnode.cpoint == '0' else [gp.co + gp.normal.normalized() * simnode.offset for gp in gpoints]
                    allpoints = numpy.zeros((len(gpoints), lvaldirecs), dtype=int8)
                    
                    for chunk in chunks(gpoints, int(scene['viparams']['nproc']) * 200):
                        for gp in chunk:
#                           Attempt to multi-process but Pool does not work with class instances
#                            p = Pool(4) 
#                            pointres = array(p.starmap(shadtree.ray_cast, [(posis[g], direc) for direc in direcs]), dtype = int8)
                            pointres = array([(0, 1)[shadtree.ray_cast(posis[g], direc)[3] == None] for direc in valdirecs], dtype = int8)
                            numpy.place(allpoints[g], pointres == 1, pointres)
                            gp[shadres] = ((numpy.sum(pointres)/lvaldirecs)).astype(float16)
                            g += 1

                        curres += len(chunk)
                        if pfile.check(curres) == 'CANCELLED':
                            return {'CANCELLED'}
              
                    shadres = [gp[shadres] for gp in gpoints]
                    o['omin']['res{}'.format(frame)], o['omax']['res{}'.format(frame)], o['oave']['res{}'.format(frame)] = min(shadres), max(shadres), sum(shadres)/len(shadres)
                    reslists.append([str(frame), 'Zone', o.name, 'X', ' '.join(['{:.3f}'.format(p[0]) for p in posis])])
                    reslists.append([str(frame), 'Zone', o.name, 'Y', ' '.join(['{:.3f}'.format(p[1]) for p in posis])])
                    reslists.append([str(frame), 'Zone', o.name, 'Z', ' '.join(['{:.3f}'.format(p[2]) for p in posis])])
                    reslists.append([str(frame), 'Zone', o.name, 'SVF', ' '.join(['{:.3f}'.format(sr) for sr in oshadres])])
                    avres.append(o['oave']['res{}'.format(frame)])
                    minres.append(o['omin']['res{}'.format(frame)])
                    maxres.append(o['omax']['res{}'.format(frame)])

            reslists.append(['All', 'Frames', '', 'Frames', ' '.join(['{}'.format(f) for f in frange])])
            reslists.append(['All', 'Zone', o.name, 'Minimum', ' '.join(['{:.3f}'.format(mr) for mr in minres])])
            reslists.append(['All', 'Zone', o.name, 'Average', ' '.join(['{:.3f}'.format(mr) for mr in avres])])
            reslists.append(['All', 'Zone', o.name, 'Maximum', ' '.join(['{:.3f}'.format(mr) for mr in maxres])])
            bm.transform(o.matrix_world.inverted())
            bm.to_mesh(o.data)
            bm.free()

        scene.vi_leg_max, scene.vi_leg_min = 1, 0

        if kivyrun.poll() is None:
            kivyrun.kill()
        
        scene.frame_start, scene.frame_end = scene['liparams']['fs'], scene['liparams']['fe']
        scene['viparams']['vidisp'] = 'svf'
        simnode['reslists'] = reslists
        simnode['frames'] = [f for f in frange]
        simnode.postexport(scene)
        return {'FINISHED'}

class NODE_OT_Shadow(bpy.types.Operator):
    bl_idname = "node.shad"
    bl_label = "Shadow Map"
    bl_description = "Undertake a shadow mapping"
    bl_register = True
    bl_undo = False
    nodeid = bpy.props.StringProperty()

    def invoke(self, context, event):
        scene = context.scene     
        scene.vi_display = 0
        
        if viparams(self, scene):            
            return {'CANCELLED'}

        shadobs = retobjs('livig')
        if not shadobs:
            self.report({'ERROR'},"No shading objects have a material attached.")
            return {'CANCELLED'}
            
        scene['liparams']['shadc'] = [ob.name for ob in retobjs('ssc')]
        if not scene['liparams']['shadc']:
            self.report({'ERROR'},"No objects have a VI Shadow material attached.")
            return {'CANCELLED'}

        scene['viparams']['restree'] = self.nodeid.split('@')[1]
        
        clearscene(scene, self)
        simnode = bpy.data.node_groups[self.nodeid.split('@')[1]].nodes[self.nodeid.split('@')[0]]
        scene['viparams']['visimcontext'] = 'Shadow'
        if not scene.get('liparams'):
           scene['liparams'] = {}
        scene['liparams']['cp'], scene['liparams']['unit'], scene['liparams']['type'] = simnode.cpoint, '% Sunlit', 'VI Shadow'
        simnode.preexport()
        (scene['liparams']['fs'], scene['liparams']['fe']) = (scene.frame_current, scene.frame_current) if simnode.animmenu == 'Static' else (simnode.startframe, simnode.endframe)
        cmap(scene)

        if simnode.starthour > simnode.endhour:
            self.report({'ERROR'},"End hour is before start hour.")
            return{'CANCELLED'}
        
        scene['viparams']['resnode'], simnode['Animation'] = simnode.name, simnode.animmenu
        (scmaxres, scminres, scavres) = [[x] * (scene['liparams']['fe'] - scene['liparams']['fs'] + 1) for x in (0, 100, 0)]
        
        frange = range(scene['liparams']['fs'], scene['liparams']['fe'] + 1)
        time = datetime.datetime(2014, simnode.sdate.month, simnode.sdate.day, simnode.starthour - 1)
        y =  2014 if simnode.edoy >= simnode.sdoy else 2014 + 1
        endtime = datetime.datetime(y, simnode.edate.month, simnode.edate.day, simnode.endhour - 1)
        interval = datetime.timedelta(hours = 1/simnode.interval)
        
        times = [time + interval*t for t in range(int((endtime - time)/interval) + simnode.interval) if simnode.starthour - 1 <= (time + interval*t).hour <= simnode.endhour  - 1]
        sps = [solarPosition(t.timetuple().tm_yday, t.hour+t.minute/60, scene.latitude, scene.longitude)[2:] for t in times]
        valmask = array([sp[0] > 0 for sp in sps], dtype = int8)
        direcs = [mathutils.Vector((-sin(sp[1]), -cos(sp[1]), tan(sp[0]))) for sp in sps]  
        valdirecs = [mathutils.Vector((-sin(sp[1]), -cos(sp[1]), tan(sp[0]))) for sp in sps if sp[0] > 0]  
        lvaldirecs = len(valdirecs)
        ilvaldirecs = 1/lvaldirecs
        calcsteps = len(frange) * sum(len([f for f in o.data.polygons if o.data.materials[f.material_index].mattype == '1']) for o in [scene.objects[on] for on in scene['liparams']['shadc']])
        curres, reslists = 0, []
        pfile = progressfile(scene['viparams']['newdir'], datetime.datetime.now(), calcsteps)
        kivyrun = progressbar(os.path.join(scene['viparams']['newdir'], 'viprogress'), 'Shadow Map')
        logentry('Conducting shadow map calculation with {} samples per hour for {} total hours and {} available sun hours'.format(simnode.interval, int(len(direcs)/simnode.interval), lvaldirecs))
        # Below uses the new string literal formatting in python 3.6
#        logentry(f'Conducting shadow map calculation with {simnode.interval} samples per hour for {int(len(direcs)/simnode.interval)} total hours and {lvaldirecs} available sun hours')
        
        for oi, o in enumerate([scene.objects[on] for on in scene['liparams']['shadc']]):
            for k in o.keys():
                del o[k]
                
            if any([s < 0 for s in o.scale]):
                logentry('Negative scaling on calculation object {}. Results may not be as expected'.format(o.name))
                self.report({'WARNING'}, 'Negative scaling on calculation object {}. Results may not be as expected'.format(o.name))

            o['omin'], o['omax'], o['oave'] = {}, {}, {}
            
            if simnode.sdoy <= simnode.edoy:
                o['days'] = arange(simnode.sdoy, simnode.edoy + 1, dtype = float)
            else:
                o['days'] = arange(simnode.sdoy, simnode.edoy + 1, dtype = float)
                
            o['hours'] = arange(simnode.starthour, simnode.endhour + 1, 1/simnode.interval, dtype = float)
            bm = bmesh.new()
            bm.from_mesh(o.data)
            clearlayers(bm, 'a')
            bm.transform(o.matrix_world)
            geom = bm.faces if simnode.cpoint == '0' else bm.verts
            geom.layers.int.new('cindex')
            cindex = geom.layers.int['cindex']
            [geom.layers.float.new('res{}'.format(fi)) for fi in frange]
            [geom.layers.float.new('hourres{}'.format(fi)) for fi in frange]
            avres, minres, maxres, g = [], [], [], 0
            
            if simnode.cpoint == '0':
                gpoints = [f for f in geom if o.data.materials[f.material_index].mattype == '1']
            elif simnode.cpoint == '1':
                gpoints = [v for v in geom if any([o.data.materials[f.material_index].mattype == '1' for f in v.link_faces])]

            for g, gp in enumerate(gpoints):
                gp[cindex] = g + 1
            
            for frame in frange: 
                g, oshadres = 0, array([])                
                scene.frame_set(frame)
                shadtree = rettree(scene, shadobs, ('', '2')[simnode.signore])
                shadres = geom.layers.float['res{}'.format(frame)]
                                  
                if gpoints:
                    posis = [gp.calc_center_bounds() + gp.normal.normalized() * simnode.offset for gp in gpoints] if simnode.cpoint == '0' else [gp.co + gp.normal.normalized() * simnode.offset for gp in gpoints]
                    allpoints = numpy.zeros((len(gpoints), len(direcs)), dtype=int8)
                    
                    for chunk in chunks(gpoints, int(scene['viparams']['nproc']) * 200):
                        for gp in chunk:
#                           Attempt to multi-process but Pool does not work with class instances
#                            p = Pool(4) 
#                            pointres = array(p.starmap(shadtree.ray_cast, [(posis[g], direc) for direc in direcs]), dtype = int8)
                            pointres = array([(0, 1)[shadtree.ray_cast(posis[g], direc)[3] == None] for direc in valdirecs], dtype = int8)
                            numpy.place(allpoints[g], valmask == 1, pointres)
                            gp[shadres] = (100 * (numpy.sum(pointres) * ilvaldirecs)).astype(float16)
                            g += 1

                        curres += len(chunk)
                        if pfile.check(curres) == 'CANCELLED':
                            return {'CANCELLED'}
                    
                    ap = numpy.average(allpoints, axis=0)                
                    shadres = [gp[shadres] for gp in gpoints]
                    o['dhres{}'.format(frame)] = array(100 * ap).reshape(len(o['days']), len(o['hours'])).T.tolist()
                    o['omin']['res{}'.format(frame)], o['omax']['res{}'.format(frame)], o['oave']['res{}'.format(frame)] = min(shadres), max(shadres), sum(shadres)/len(shadres)
                    reslists.append([str(frame), 'Zone', o.name, 'X', ' '.join(['{:.3f}'.format(p[0]) for p in posis])])
                    reslists.append([str(frame), 'Zone', o.name, 'Y', ' '.join(['{:.3f}'.format(p[1]) for p in posis])])
                    reslists.append([str(frame), 'Zone', o.name, 'Z', ' '.join(['{:.3f}'.format(p[2]) for p in posis])])
                    reslists.append([str(frame), 'Zone', o.name, 'Sunlit %', ' '.join(['{:.3f}'.format(sr) for sr in oshadres])])
                    avres.append(o['oave']['res{}'.format(frame)])
                    minres.append(o['omin']['res{}'.format(frame)])
                    maxres.append(o['omax']['res{}'.format(frame)])
            
            reslists.append(['All', 'Frames', '', 'Frames', ' '.join(['{}'.format(f) for f in frange])])
            reslists.append(['All', 'Zone', o.name, 'Minimum', ' '.join(['{:.3f}'.format(mr) for mr in minres])])
            reslists.append(['All', 'Zone', o.name, 'Average', ' '.join(['{:.3f}'.format(mr) for mr in avres])])
            reslists.append(['All', 'Zone', o.name, 'Maximum', ' '.join(['{:.3f}'.format(mr) for mr in maxres])])
            
            bm.transform(o.matrix_world.inverted())
            bm.to_mesh(o.data)
            bm.free()

        scene.vi_leg_max, scene.vi_leg_min = 100, 0

        if kivyrun.poll() is None:
            kivyrun.kill()
        
        scene.frame_start, scene.frame_end = scene['liparams']['fs'], scene['liparams']['fe']
        simnode['reslists'] = reslists
        simnode['frames'] = [f for f in frange]
        simnode['year'] = 2015
        simnode.postexport(scene)
        scene['viparams']['vidisp'] = 'ss'
        return {'FINISHED'}
    
class VIEW3D_OT_SVFDisplay(bpy.types.Operator):
    '''Display results legend and stats in the 3D View'''
    bl_idname = "view3d.svfdisplay"
    bl_label = "Shadow study metric display"
    bl_description = "Display shadow study metrics"
    bl_register = True
    bl_undo = False

    def modal(self, context, event):         
        redraw = 0  
        
        if self.scene.vi_display == 0 or self.scene['viparams']['vidisp'] != 'svfpanel' or not [o.lires for o in bpy.data.objects]:
            bpy.types.SpaceView3D.draw_handler_remove(self._handle_svf_disp, 'WINDOW')
            bpy.types.SpaceView3D.draw_handler_remove(self._handle_pointres, 'WINDOW')
            context.area.tag_redraw()
            self.scene['viparams']['vidisp'] = 'svf'
            return {'CANCELLED'}
            
        if event.type != 'INBETWEEN_MOUSEMOVE' and context.region and context.area.type == 'VIEW_3D' and context.region.type == 'WINDOW':  
            mx, my = event.mouse_region_x, event.mouse_region_y 
            
            if any((self.scene.vi_leg_levels != self.legend.levels, self.scene.vi_leg_col != self.legend.col, self.scene.vi_leg_scale != self.legend.scale, (self.legend.minres, self.legend.maxres) != leg_min_max(self.scene))):               
                self.legend.update(context)
                redraw = 1
            
            # Legend routine 
            
            if self.legend.spos[0] < mx < self.legend.epos[0] and self.legend.spos[1] < my < self.legend.epos[1]:
                if self.legend.hl != (0, 1, 1, 1):  
                    self.legend.hl = (0, 1, 1, 1) 
                    redraw = 1
                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.legend.press = 1
                        self.legend.move = 0
                        return {'RUNNING_MODAL'}
                    elif event.value == 'RELEASE':
                        if not self.legend.move:
                            self.legend.expand = 0 if self.legend.expand else 1
                        self.legend.press = 0
                        self.legend.move = 0
                        context.area.tag_redraw()
                        return {'RUNNING_MODAL'}
                
                elif event.type == 'ESC':
                    bpy.types.SpaceView3D.draw_handler_remove(self._handle_ss_disp, 'WINDOW')
                    context.area.tag_redraw()
                    return {'CANCELLED'}
                    
                elif self.legend.press and event.type == 'MOUSEMOVE':
                     self.legend.move = 1
                     self.legend.press = 0
            
            elif abs(self.legend.lepos[0] - mx) < 10 and abs(self.legend.lspos[1] - my) < 10:
                self.legend.hl = (0, 1, 1, 1) 
                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.legend.resize = 1
                    if self.legend.resize and event.value == 'RELEASE':
                        self.legend.resize = 0
                    return {'RUNNING_MODAL'}
                    
            else:
                if self.legend.hl != (1, 1, 1, 1):
                    self.legend.hl = (1, 1, 1, 1)
                    redraw = 1
                                             
            # Resize routines
            
            if abs(self.legend.lepos[0] - mx) < 20 and abs(self.legend.lspos[1] - my) < 20:
                self.legend.hl = (0, 1, 1, 1) 
                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.legend.resize = 1
                    if self.legend.resize and event.value == 'RELEASE':
                        self.legend.resize = 0
                    return {'RUNNING_MODAL'}
                        
            # Move routines
                     
            if event.type == 'MOUSEMOVE':                
                if self.legend.move:
                    self.legend.pos = [mx, my]
                    redraw = 1
                if self.legend.resize:
                    self.legend.lepos[0], self.legend.lspos[1] = mx, my
                    redraw = 1
            
            if redraw:
                context.area.tag_redraw()

        return {'PASS_THROUGH'}

    def invoke(self, context, event):
        self.scene = context.scene
        try:
            bpy.types.SpaceView3D.draw_handler_remove(self._handle_pointres, 'WINDOW')
            bpy.types.SpaceView3D.draw_handler_remove(self._handle_svf_disp, 'WINDOW')
        except:
            pass
        clearscene(self.scene, self)
        self.scene['viparams']['vidisp'] = 'svfpanel'
        self.simnode = bpy.data.node_groups[self.scene['viparams']['restree']].nodes[self.scene['viparams']['resnode']]
        li_display(self, self.simnode)
        self.scene.vi_disp_wire, self.scene.vi_display = 1, 1
        lnd = linumdisplay(self, context, self.simnode)
        self._handle_pointres = bpy.types.SpaceView3D.draw_handler_add(lnd.draw, (context, ), 'WINDOW', 'POST_PIXEL')
        self.legend = svf_legend([80, context.region.height - 40], context.region.width, context.region.height, 'legend.png', 150, 600)
        self.legend.update(context)
        self._handle_svf_disp = bpy.types.SpaceView3D.draw_handler_add(svf_disp, (self, context, self.simnode), 'WINDOW', 'POST_PIXEL')
        context.window_manager.modal_handler_add(self)
        return {'RUNNING_MODAL'}
        
class VIEW3D_OT_SSDisplay(bpy.types.Operator):
    '''Display results legend and stats in the 3D View'''
    bl_idname = "view3d.ssdisplay"
    bl_label = "Shadow study metric display"
    bl_description = "Display shadow study metrics"
    bl_register = True
    bl_undo = False

    def modal(self, context, event):  
        redraw = 0 
        
        if self.scene.vi_display == 0 or self.scene['viparams']['vidisp'] != 'sspanel' or not [o.lires for o in bpy.data.objects]:
            bpy.types.SpaceView3D.draw_handler_remove(self._handle_ss_disp, 'WINDOW')
            bpy.types.SpaceView3D.draw_handler_remove(self._handle_pointres, 'WINDOW')
            context.area.tag_redraw()
            self.scene['viparams']['vidisp'] = 'ss'
            return {'CANCELLED'}        
        
        if event.type != 'INBETWEEN_MOUSEMOVE' and context.region and context.area.type == 'VIEW_3D' and context.region.type == 'WINDOW':                    
            mx, my = event.mouse_region_x, event.mouse_region_y 
            
            if any((self.scene.vi_leg_levels != self.legend.levels, self.scene.vi_leg_col != self.legend.col, self.scene.vi_leg_scale != self.legend.scale, (self.legend.minres, self.legend.maxres) != leg_min_max(self.scene))):               
                self.legend.update(context)                
                redraw = 1
                 
            # Legend routine 
            
            if self.legend.spos[0] < mx < self.legend.epos[0] and self.legend.spos[1] < my < self.legend.epos[1]:
                if self.legend.hl != (0, 1, 1, 1):  
                    self.legend.hl = (0, 1, 1, 1) 
                    redraw = 1
                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.legend.press = 1
                        self.legend.move = 0
                        return {'RUNNING_MODAL'}
                    elif event.value == 'RELEASE':
                        if not self.legend.move:
                            self.legend.expand = 0 if self.legend.expand else 1
                        self.legend.press = 0
                        self.legend.move = 0
                        context.area.tag_redraw()
                        return {'RUNNING_MODAL'}
                
                elif event.type == 'ESC':
                    bpy.types.SpaceView3D.draw_handler_remove(self._handle_ss_disp, 'WINDOW')
                    context.area.tag_redraw()
                    return {'CANCELLED'}
                    
                elif self.legend.press and event.type == 'MOUSEMOVE':
                     self.legend.move = 1
                     self.legend.press = 0
            
            elif abs(self.legend.lepos[0] - mx) < 10 and abs(self.legend.lspos[1] - my) < 10:
                self.legend.hl = (0, 1, 1, 1) 
                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.legend.resize = 1
                    if self.legend.resize and event.value == 'RELEASE':
                        self.legend.resize = 0
                    return {'RUNNING_MODAL'}
                    
            else:
                if self.legend.hl != (1, 1, 1, 1):
                    self.legend.hl = (1, 1, 1, 1)
                    redraw = 1
            
            # Scatter routine
                
            if self.dhscatter.spos[0] < mx < self.dhscatter.epos[0] and self.dhscatter.spos[1] < my < self.dhscatter.epos[1]:
                if self.dhscatter.hl != (0, 1, 1, 1):  
                    self.dhscatter.hl = (0, 1, 1, 1)
                    redraw = 1
                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.dhscatter.press = 1
                        self.dhscatter.move = 0
                        return {'RUNNING_MODAL'}
                    elif event.value == 'RELEASE':
                        if not self.dhscatter.move:
                            self.dhscatter.expand = 0 if self.dhscatter.expand else 1
                        self.dhscatter.press = 0
                        self.dhscatter.move = 0
                        context.area.tag_redraw()
                        return {'RUNNING_MODAL'}
                
                elif event.type == 'ESC':
                    bpy.data.images.remove(self.dhscatter.gimage)
                    self.dhscatter.plt.close()
                    bpy.types.SpaceView3D.draw_handler_remove(self._handle_wr_disp, 'WINDOW')
                    context.area.tag_redraw()
                    return {'CANCELLED'}
                    
                elif self.dhscatter.press and event.type == 'MOUSEMOVE':
                     self.dhscatter.move = 1
                     self.dhscatter.press = 0
        
            elif self.dhscatter.lspos[0] < mx < self.dhscatter.lepos[0] and self.dhscatter.lspos[1] < my < self.dhscatter.lepos[1] and abs(self.dhscatter.lepos[0] - mx) > 20 and abs(self.dhscatter.lspos[1] - my) > 20:
                if self.dhscatter.expand: 
                    self.dhscatter.hl = (0, 1, 1, 1)
                    if event.type == 'LEFTMOUSE' and event.value == 'PRESS' and self.dhscatter.expand and self.dhscatter.lspos[0] < mx < self.dhscatter.lepos[0] and self.dhscatter.lspos[1] < my < self.dhscatter.lspos[1] + 0.9 * self.dhscatter.ydiff:
                        self.dhscatter.show_plot()

                    context.area.tag_redraw()
                    return {'RUNNING_MODAL'}
                   
            else:
                if self.dhscatter.hl != (1, 1, 1, 1):
                    self.dhscatter.hl = (1, 1, 1, 1)
                    redraw = 1

            # Update routine
                
            if self.dhscatter.frame != self.scene.frame_current or self.dhscatter.cao != context.active_object or self.dhscatter.col != self.scene.vi_leg_col:
                self.dhscatter.update(context)
                redraw = 1
                         
            # Resize routines
            
            if abs(self.legend.lepos[0] - mx) < 20 and abs(self.legend.lspos[1] - my) < 20:
                self.legend.hl = (0, 1, 1, 1) 
                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.legend.resize = 1
                    if self.legend.resize and event.value == 'RELEASE':
                        self.legend.resize = 0
                    return {'RUNNING_MODAL'}
            
            if abs(self.dhscatter.lepos[0] - mx) < 20 and abs(self.dhscatter.lspos[1] - my) < 20:
                self.dhscatter.hl = (0, 1, 1, 1) 
                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.dhscatter.resize = 1
                    if self.dhscatter.resize and event.value == 'RELEASE':
                        self.dhscatter.resize = 0
                    return {'RUNNING_MODAL'}
            
            # Move routines
                     
            if event.type == 'MOUSEMOVE':                
                if self.legend.move:
                    self.legend.pos = [mx, my]
                    redraw = 1
                if self.legend.resize:
                    self.legend.lepos[0], self.legend.lspos[1] = mx, my
                    redraw = 1
                if self.dhscatter.move:
                    self.dhscatter.pos = [mx, my]
                    redraw = 1
                if self.dhscatter.resize:
                    self.dhscatter.lepos[0], self.dhscatter.lspos[1] = mx, my
                    redraw = 1
            
            if redraw:
                context.area.tag_redraw()

        return {'PASS_THROUGH'}

    def invoke(self, context, event):
        self.scene = context.scene
        
        try:
            bpy.types.SpaceView3D.draw_handler_remove(self._handle_pointres, 'WINDOW')
            bpy.types.SpaceView3D.draw_handler_remove(self._handle_ss_disp, 'WINDOW')
        except:
            pass
        
        clearscene(self.scene, self)
        self.scene.vi_display = 1
        self.scene['viparams']['vidisp'] = 'sspanel'
        self.simnode = bpy.data.node_groups[self.scene['viparams']['restree']].nodes[self.scene['viparams']['resnode']]
        li_display(self, self.simnode)
        self.scene.vi_disp_wire, self.scene.vi_display = 1, 1
        lnd = linumdisplay(self, context, self.simnode)
        self._handle_pointres = bpy.types.SpaceView3D.draw_handler_add(lnd.draw, (context, ), 'WINDOW', 'POST_PIXEL')
        self.legend = ss_legend([80, context.region.height - 40], context.region.width, context.region.height, 'legend.png', 150, 600)
        self.dhscatter = ss_scatter([160, context.region.height - 40], context.region.width, context.region.height, 'scat.png', 600, 400)
        self.legend.update(context)
        self.dhscatter.update(context)
        self._handle_ss_disp = bpy.types.SpaceView3D.draw_handler_add(ss_disp, (self, context, self.simnode), 'WINDOW', 'POST_PIXEL')
        context.window_manager.modal_handler_add(self)
        return {'RUNNING_MODAL'}
        
class VIEW3D_OT_LiViBasicDisplay(bpy.types.Operator):
    '''Display results legend and stats in the 3D View'''
    bl_idname = "view3d.livibasicdisplay"
    bl_label = "LiVi basic metric display"
    bl_description = "Display basic lighting metrics"
    bl_register = True
    bl_undo = False

    def modal(self, context, event): 
        redraw = 0 
        
        if self.scene.vi_display == 0 or context.scene['viparams']['vidisp'] != 'lipanel' or not any([o.lires for o in bpy.data.objects]):
            self.scene.vi_display = 0
            bpy.types.SpaceView3D.draw_handler_remove(self._handle_disp, 'WINDOW')
            bpy.types.SpaceView3D.draw_handler_remove(self._handle_pointres, 'WINDOW')
            self.scene['viparams']['vidisp'] = 'li'
            context.area.tag_redraw()
            return {'CANCELLED'}

        if event.type != 'INBETWEEN_MOUSEMOVE' and context.region and context.area.type == 'VIEW_3D' and context.region.type == 'WINDOW':            
            mx, my = event.mouse_region_x, event.mouse_region_y 
            
            if any((self.scene.vi_leg_levels != self.legend.levels, self.scene.vi_leg_col != self.legend.col, self.scene.vi_leg_scale != self.legend.scale, (self.legend.minres, self.legend.maxres) != leg_min_max(self.scene))):               
                self.legend.update(context)
                redraw = 1
            
            # Legend routine 
            
            if self.legend.spos[0] < mx < self.legend.epos[0] and self.legend.spos[1] < my < self.legend.epos[1]:
                if self.legend.hl != (0, 1, 1, 1):
                    self.legend.hl = (0, 1, 1, 1)
                    redraw = 1  
                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.legend.press = 1
                        self.legend.move = 0
                        return {'RUNNING_MODAL'}
                    elif event.value == 'RELEASE':
                        if not self.legend.move:
                            self.legend.expand = 0 if self.legend.expand else 1
                        self.legend.press = 0
                        self.legend.move = 0
                        context.area.tag_redraw()
                        return {'RUNNING_MODAL'}
                
                elif event.type == 'ESC':
                    bpy.types.SpaceView3D.draw_handler_remove(self._handle_disp, 'WINDOW')
                    bpy.types.SpaceView3D.draw_handler_remove(self._handle_pointres, 'WINDOW')
                    self.scene['viparams']['vidisp'] = 'li'
                    context.area.tag_redraw()
                    return {'CANCELLED'}
                    
                elif self.legend.press and event.type == 'MOUSEMOVE':
                     self.legend.move = 1
                     self.legend.press = 0
            
            elif abs(self.legend.lepos[0] - mx) < 10 and abs(self.legend.lspos[1] - my) < 10:
                if self.legend.hl != (0, 1, 1, 1):
                    self.legend.hl = (0, 1, 1, 1)
                    redraw = 1
                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.legend.resize = 1
                    if self.legend.resize and event.value == 'RELEASE':
                        self.legend.resize = 0
                    return {'RUNNING_MODAL'}
                    
            elif self.legend.hl != (1, 1, 1, 1):
                self.legend.hl = (1, 1, 1, 1)
                redraw = 1

            # Table routine
            
            if self.frame != context.scene.frame_current or self.table.unit != context.scene['liparams']['unit'] or self.table.cao != context.active_object:
                self.table.update(context)                
                redraw = 1
            
            if self.table.spos[0] < mx < self.table.epos[0] and self.table.spos[1] < my < self.table.epos[1]:
                if self.table.hl != (0, 1, 1, 1):
                    self.table.hl = (0, 1, 1, 1)
                    redraw = 1  
                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.table.press = 1
                        self.table.move = 0
                        return {'RUNNING_MODAL'}
                    elif event.value == 'RELEASE':
                        if not self.table.move:
                            self.table.expand = 0 if self.table.expand else 1
                        self.table.press = 0
                        self.table.move = 0
                        context.area.tag_redraw()
                        return {'RUNNING_MODAL'}
                
                elif event.type == 'ESC':
                    bpy.types.SpaceView3D.draw_handler_remove(self._handle_disp, 'WINDOW')
                    bpy.types.SpaceView3D.draw_handler_remove(self._handle_pointres, 'WINDOW')
                    context.scene['viparams']['vidisp'] = 'li'
                    context.area.tag_redraw()
                    return {'CANCELLED'}
                    
                elif self.table.press and event.type == 'MOUSEMOVE':
                     self.table.move = 1
                     self.table.press = 0                     
            elif self.table.hl != (1, 1, 1, 1):
                self.table.hl = (1, 1, 1, 1)
                redraw = 1
                
            if context.scene['viparams']['visimcontext'] == 'LiVi Compliance':
                if self.frame != context.scene.frame_current:
                    self.tablecomp.update(context)
                    redraw = 1
                if self.tablecomp.unit != context.scene['liparams']['unit']:
                    self.tablecomp.update(context)
                    self.tablecomp.unit = context.scene['liparams']['unit']
                    redraw = 1
                if self.tablecomp.cao != context.active_object:
                    self.tablecomp.update(context)
                    redraw = 1
                
                if self.tablecomp.spos[0] < mx < self.tablecomp.epos[0] and self.tablecomp.spos[1] < my < self.tablecomp.epos[1]:
                    if self.tablecomp.hl != (0, 1, 1, 1):
                        self.tablecomp.hl = (0, 1, 1, 1)
                        redraw = 1  
                    if event.type == 'LEFTMOUSE':
                        if event.value == 'PRESS':
                            self.tablecomp.press = 1
                            self.tablecomp.move = 0
                            return {'RUNNING_MODAL'}
                        elif event.value == 'RELEASE':
                            if not self.tablecomp.move:
                                self.tablecomp.expand = 0 if self.tablecomp.expand else 1
                            self.tablecomp.press = 0
                            self.tablecomp.move = 0
                            context.area.tag_redraw()
                            return {'RUNNING_MODAL'}
                    
                    elif event.type == 'ESC':
                        bpy.types.SpaceView3D.draw_handler_remove(self._handle_disp, 'WINDOW')
                        bpy.types.SpaceView3D.draw_handler_remove(self._handle_pointres, 'WINDOW')
                        context.scene['viparams']['vidisp'] = 'li'
                        context.area.tag_redraw()
                        return {'CANCELLED'}
                        
                    elif self.tablecomp.press and event.type == 'MOUSEMOVE':
                         self.tablecomp.move = 1
                         self.tablecomp.press = 0                     
                elif self.tablecomp.hl != (1, 1, 1, 1):
                    self.tablecomp.hl = (1, 1, 1, 1)
                    redraw = 1
                
            if context.scene['liparams']['unit'] in ('ASE (hrs)', 'sDA (%)', 'DA (%)', 'UDI-f (%)', 'UDI-s (%)', 'UDI-e (%)', 'UDI-a (%)', 'Max lux', 'Min lux', 'Avg lux', 'kWh', 'kWh/m2'):
                if self.dhscatter.frame != context.scene.frame_current:
                    self.dhscatter.update(context)
                    redraw = 1
                if self.dhscatter.unit != context.scene['liparams']['unit']:
                    self.dhscatter.update(context)
                    redraw = 1
                if self.dhscatter.cao != context.active_object:
                    self.dhscatter.update(context)
                    redraw = 1
                if self.dhscatter.col != context.scene.vi_leg_col:
                    self.dhscatter.update(context)
                    redraw = 1
                if context.scene['liparams']['unit'] in ('Max lux', 'Min lux', 'Avg lux', 'kWh', 'kWh/m2'):
                    if (self.dhscatter.vmin, self.dhscatter.vmax) != (context.scene.vi_scatter_min, context.scene.vi_scatter_max):
                       self.dhscatter.update(context) 
                       redraw = 1
                                        
                if self.dhscatter.spos[0] < mx < self.dhscatter.epos[0] and self.dhscatter.spos[1] < my < self.dhscatter.epos[1]:
                    if self.dhscatter.hl != (0, 1, 1, 1):
                        self.dhscatter.hl = (0, 1, 1, 1)
                        redraw = 1 
                    if event.type == 'LEFTMOUSE':
                        if event.value == 'PRESS':
                            self.dhscatter.press = 1
                            self.dhscatter.move = 0
                            return {'RUNNING_MODAL'}
                        elif event.value == 'RELEASE':
                            if not self.dhscatter.move:
                                self.dhscatter.expand = 0 if self.dhscatter.expand else 1
                            self.dhscatter.press = 0
                            self.dhscatter.move = 0
                            context.area.tag_redraw()
                            return {'RUNNING_MODAL'}
                    
                    elif event.type == 'ESC':
                        bpy.types.SpaceView3D.draw_handler_remove(self._handle_disp, 'WINDOW')
                        bpy.types.SpaceView3D.draw_handler_remove(self._handle_pointres, 'WINDOW')
                        context.scene['viparams']['vidisp'] = 'li'
                        context.area.tag_redraw()
                        return {'CANCELLED'}
                        
                    elif self.dhscatter.press and event.type == 'MOUSEMOVE':
                         self.dhscatter.move = 1
                         self.dhscatter.press = 0   
                                            
                else:
                    if self.dhscatter.hl != (1, 1, 1, 1):
                        self.dhscatter.hl = (1, 1, 1, 1)
                        redraw = 1
                    if self.dhscatter.lspos[0] < mx < self.dhscatter.lepos[0] and self.dhscatter.lspos[1] < my < self.dhscatter.lepos[1] and abs(self.dhscatter.lepos[0] - mx) > 20 and abs(self.dhscatter.lspos[1] - my) > 20:
                        if self.dhscatter.expand: 
                            self.dhscatter.hl = (1, 1, 1, 1)
                            if event.type == 'LEFTMOUSE' and event.value == 'PRESS' and self.dhscatter.expand and self.dhscatter.lspos[0] < mx < self.dhscatter.lepos[0] and self.dhscatter.lspos[1] < my < self.dhscatter.lspos[1] + 0.9 * self.dhscatter.ydiff:
                                self.dhscatter.show_plot()
                                                         
            # Resize routines
            
            if abs(self.legend.lepos[0] - mx) < 20 and abs(self.legend.lspos[1] - my) < 20:
                self.legend.hl = (0, 1, 1, 1) 
                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.legend.resize = 1
                    if self.legend.resize and event.value == 'RELEASE':
                        self.legend.resize = 0
                    return {'RUNNING_MODAL'}
                    
            elif abs(self.table.lepos[0] - mx) < 20 and abs(self.table.lspos[1] - my) < 20:
                self.table.hl = (0, 1, 1, 1) 
                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.table.resize = 1
                    if self.table.resize and event.value == 'RELEASE':
                        self.table.resize = 0
                    return {'RUNNING_MODAL'}
            
            elif context.scene['viparams']['visimcontext'] == 'LiVi Compliance' and abs(self.tablecomp.lepos[0] - mx) < 20 and abs(self.tablecomp.lspos[1] - my) < 20:
                self.tablecomp.hl = (0, 1, 1, 1) 
                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.tablecomp.resize = 1
                    if self.tablecomp.resize and event.value == 'RELEASE':
                        self.tablecomp.resize = 0
                    return {'RUNNING_MODAL'}

            elif context.scene['liparams']['unit'] in ('ASE (hrs)', 'sDA (%)', 'DA (%)', 'UDI-s (%)', 'UDI-e (%)', 'UDI-f (%)', 'UDI-a (%)', 'Max lux', 'Min lux', 'Avg lux', 'kWh', 'kWh/m2') and abs(self.dhscatter.lepos[0] - mx) < 20 and abs(self.dhscatter.lspos[1] - my) < 20:
                self.dhscatter.hl = (0, 1, 1, 1) 
                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.dhscatter.resize = 1
                    if self.dhscatter.resize and event.value == 'RELEASE':
                        self.dhscatter.resize = 0
                    return {'RUNNING_MODAL'}
            # Move routines
                     
            if event.type == 'MOUSEMOVE':                
                if self.legend.move:
                    self.legend.pos = [mx, my]
                    redraw = 1
                if self.legend.resize:
                    self.legend.lepos[0], self.legend.lspos[1] = mx, my
                    redraw = 1
                if self.table.move:
                    self.table.pos = [mx, my]
                    redraw = 1
                if self.table.resize:
                    self.table.lepos[0], self.table.lspos[1] = mx, my
                    redraw = 1
                if context.scene['viparams']['visimcontext'] == 'LiVi Compliance':
                    if self.tablecomp.move:
                        self.tablecomp.pos = [mx, my]
                        redraw = 1
                    if self.tablecomp.resize:
                        self.tablecomp.lepos[0], self.tablecomp.lspos[1] = mx, my
                        redraw = 1
                try:
                    if self.dhscatter.move:
                        self.dhscatter.pos = [mx, my]
                        redraw = 1
                    if self.dhscatter.resize:
                        self.dhscatter.lepos[0], self.dhscatter.lspos[1] = mx, my
                        redraw = 1
                except:
                    pass
                                
            if redraw:
                context.area.tag_redraw()
                self.frame = context.scene.frame_current
                
        return {'PASS_THROUGH'}

    def invoke(self, context, event):
        self.scene = context.scene
        clearscene(self.scene, self)
        self.scene.vi_display = 1
        self.scene['viparams']['vidisp'] = 'lipanel'
        self.simnode = bpy.data.node_groups[self.scene['viparams']['restree']].nodes[self.scene['viparams']['resnode']]
        self.frame = self.scene.frame_current
        
        if li_display(self, self.simnode) == 'CANCELLED':
            return {'CANCELLED'}
        
        self.scene.vi_disp_wire, self.scene.vi_display = 1, 1
        lnd = linumdisplay(self, context, self.simnode)
        self._handle_pointres = bpy.types.SpaceView3D.draw_handler_add(lnd.draw, (context, ), 'WINDOW', 'POST_PIXEL')
        self.legend = basic_legend([80, context.region.height - 40], context.region.width, context.region.height, 'legend.png', 150, 600)
        self.legend.update(context)
#        self.dhscatter = wr_scatter([160, context.region.height - 40], context.region.width, context.region.height, 'stats.png', 600, 400)
#        if scene['viparams']['visimcontext'] == 'LiVi Basic':
        self.table = basic_table([240, context.region.height - 40], context.region.width, context.region.height, 'table.png', 600, 100)  
        self.table.update(context)
        
        if self.scene['viparams']['visimcontext'] == 'LiVi Compliance':
            self.tablecomp = comp_table([300, context.region.height - 40], context.region.width, context.region.height, 'compliance.png', 600, 200)
            self.tablecomp.update(context)
            if self.simnode['coptions']['canalysis'] == '3':
                self.dhscatter = leed_scatter([160, context.region.height - 40], context.region.width, context.region.height, 'scat.png', 600, 400)
                self.dhscatter.update(context)        
            self._handle_disp = bpy.types.SpaceView3D.draw_handler_add(comp_disp, (self, context, self.simnode), 'WINDOW', 'POST_PIXEL')

#        self.dhscatter.update(context)
        
#        self._handle_spnum = bpy.types.SpaceView3D.draw_handler_add(viwr_legend, (self, context, simnode), 'WINDOW', 'POST_PIXEL')
        elif self.scene['viparams']['visimcontext'] == 'LiVi Basic':
            self._handle_disp = bpy.types.SpaceView3D.draw_handler_add(basic_disp, (self, context, self.simnode), 'WINDOW', 'POST_PIXEL')
#        if self.scene['viparams']['visimcontext'] == 'LiVi Compliance':
#            self._handle_disp = bpy.types.SpaceView3D.draw_handler_add(comp_disp, (self, context, self.simnode), 'WINDOW', 'POST_PIXEL')
        elif self.scene['viparams']['visimcontext'] == 'LiVi CBDM':
            if self.simnode['coptions']['cbanalysis'] != '0':
                self.dhscatter = cbdm_scatter([160, context.region.height - 40], context.region.width, context.region.height, 'scat.png', 600, 400)
                self.dhscatter.update(context)
            self._handle_disp = bpy.types.SpaceView3D.draw_handler_add(cbdm_disp, (self, context, self.simnode), 'WINDOW', 'POST_PIXEL')
        context.window_manager.modal_handler_add(self)
        return {'RUNNING_MODAL'}

# Openfoam operators

class NODE_OT_Blockmesh(bpy.types.Operator):
    bl_idname = "node.blockmesh"
    bl_label = "Blockmesh export"
    bl_description = "Export an Openfoam blockmesh"
    bl_register = True
    bl_undo = False
    nodeid = bpy.props.StringProperty()

    def execute(self, context):
        scene = context.scene
        expnode = context.node if context.node.bl_label == "FloVi BlockMesh" else context.node.inputs[0].links[0].from_node
        bmos = [o for o in scene.objects if o.vi_type == '2']
        
        if viparams(self, scene):
            return {'CANCELLED'} 
        
        if len(bmos) != 1:
            self.report({'ERROR'},"One and only one object with the CFD Domain property is allowed")
            return {'CANCELLED'}
        elif [f.material_index for f in bmos[0].data.polygons if f.material_index + 1 > len(bmos[0].data.materials)]:
            self.report({'ERROR'},"Not every domain face has a material attached")
            logentry("Not every face has a material attached")
            return {'CANCELLED'}
        with open(os.path.join(scene['flparams']['ofsfilebase'], 'controlDict'), 'w') as cdfile:
            cdfile.write(fvcdwrite("simpleFoam", 0.005, 5))
        with open(os.path.join(scene['flparams']['ofsfilebase'], 'fvSolution'), 'w') as fvsolfile:
            fvsolfile.write(fvsolwrite(expnode))
        with open(os.path.join(scene['flparams']['ofsfilebase'], 'fvSchemes'), 'w') as fvschfile:
            fvschfile.write(fvschwrite(expnode))
        with open(os.path.join(scene['flparams']['ofcpfilebase'], 'blockMeshDict'), 'w') as bmfile:
            bmfile.write(fvbmwrite(bmos[0], expnode))

        call(("blockMesh", "-case", "{}".format(scene['flparams']['offilebase'])))
        fvblbmgen(bmos[0].data.materials, open(os.path.join(scene['flparams']['ofcpfilebase'], 'faces'), 'r'), open(os.path.join(scene['flparams']['ofcpfilebase'], 'points'), 'r'), open(os.path.join(scene['flparams']['ofcpfilebase'], 'boundary'), 'r'), 'blockMesh')
        expnode.export()
        return {'FINISHED'}

class NODE_OT_Snappymesh(bpy.types.Operator):
    bl_idname = "node.snappy"
    bl_label = "SnappyHexMesh export"
    bl_description = "Export an Openfoam snappyhexmesh"
    bl_register = True
    bl_undo = True
    nodeid = bpy.props.StringProperty()

    def execute(self, context):
        bpy.ops.node.blockmesh()        
        scene, mats = context.scene, []

        for dirname in os.listdir(scene['flparams']['offilebase']):
            if os.path.isdir(os.path.join(scene['flparams']['offilebase'], dirname)) and dirname not in ('0', 'constant', 'system'):
                shutil.rmtree(os.path.join(scene['flparams']['offilebase'], dirname))
        for fname in os.listdir(scene['flparams']['ofcpfilebase']):
            if os.path.isfile(os.path.join(scene['flparams']['ofcpfilebase'], fname)) and fname in ('cellLevel', 'pointLevel', 'surfaceIndex', 'level0Edge', 'refinementHistory'):
                os.remove(os.path.join(scene['flparams']['ofcpfilebase'], fname))

        expnode = context.node
#        expnode = bpy.data.node_groups[self.nodeid.split('@')[1]].nodes[self.nodeid.split('@')[0]]
        fvos = [o for o in scene.objects if o.vi_type == '3']
        
        if fvos:
            bmos = [o for o in scene.objects if o.vi_type == '2']
#                bpy.ops.export_mesh.stl(filepath=os.path.join(scene['flparams']['ofctsfilebase'], '{}.obj'.format(o.name)), check_existing=False, filter_glob="*.stl", axis_forward='Y', axis_up='Z', global_scale=1.0, use_scene_unit=True, ascii=False, use_mesh_modifiers=True)
            fvobjwrite(scene, fvos, bmos[0])
#            bpy.ops.export_scene.obj(check_existing=True, filepath=os.path.join(scene['flparams']['ofctsfilebase'], '{}.obj'.format(fvos[0].name)), axis_forward='Y', axis_up='Z', filter_glob="*.obj;*.mtl", use_selection=True, use_animation=False, use_mesh_modifiers=True, use_edges=True, use_smooth_groups=False, use_smooth_groups_bitflags=False, use_normals=False, use_uvs=True, use_materials=True, use_triangles=True, use_nurbs=False, use_vertex_groups=False, use_blen_objects=True, group_by_object=False, group_by_material=True, keep_vertex_order=False, global_scale=1.0, path_mode='AUTO')
            gmats = [mat for mat in fvos[0].data.materials if mat.flovi_ground]
#            if gmats:
            with open(os.path.join(scene['flparams']['ofsfilebase'], 'snappyHexMeshDict'), 'w') as shmfile:
                shmfile.write(fvshmwrite(expnode, fvos, bmos[0], ground = gmats))
            with open(os.path.join(scene['flparams']['ofsfilebase'], 'meshQualityDict'), 'w') as mqfile:
                mqfile.write(fvmqwrite())
            with open(os.path.join(scene['flparams']['ofsfilebase'], 'surfaceFeatureExtractDict'), 'w') as sfefile:
                sfefile.write(fvsfewrite(fvos))
        call(('surfaceFeatureExtract', "-case", "{}".format(scene['flparams']['offilebase'])))
        call(('snappyHexMesh', "-overwrite", "-case", "{}".format(scene['flparams']['offilebase'])))
        
        for mat in fvos[0].data.materials:
#            mat.name = '{}_{}'.format(fvos[0].name, mat.name)
            mats.append(mat)
        for mat in [o for o in scene.objects if o.vi_type == '2'][0].data.materials:
            mats.append(mat)
        fvblbmgen(mats, open(os.path.join(scene['flparams']['ofcpfilebase'], 'faces'), 'r'), open(os.path.join(scene['flparams']['ofcpfilebase'], 'points'), 'r'), open(os.path.join(scene['flparams']['ofcpfilebase'], 'boundary'), 'r'), 'hexMesh')

        expnode.export()
        return {'FINISHED'}
    
class NODE_OT_FVExport(bpy.types.Operator):
    bl_idname = "node.fvexport"
    bl_label = "FloVi Export"
    bl_description = "Export an OpenFOAM case"
    bl_register = True
    bl_undo = True
    nodeid = bpy.props.StringProperty()
    
    def execute(self, context):
        scene = context.scene
        expnode = bpy.data.node_groups[self.nodeid.split('@')[1]].nodes[self.nodeid.split('@')[0]]
        self.residuals  = expnode.preexport(scene)
#        self.fpfile = os.path.join(scene['viparams']['newdir'], 'floviprogress')
#        self.pfile = fvprogressfile(scene['viparams']['newdir'])
#        self.pfile = progressfile(scene['viparams']['newdir'], datetime.datetime.now(), 100)
#        self.pfile = progressfile(self.folder, datetime.datetime.now(), 100)
#        self.kivyrun = fvprogressbar(os.path.join(scene['viparams']['newdir'], 'viprogress'), str(self.residuals))
        bmos = [o for o in scene.objects if o.vi_type in ('2', '3')]
        
       
        with open(os.path.join(scene['flparams']['ofsfilebase'], 'controlDict'), 'w') as cdfile:
            cdfile.write(fvcdwrite(expnode.solver, expnode.dt, expnode.et))
        fvvarwrite(scene, bmos, expnode)
        with open(os.path.join(scene['flparams']['ofsfilebase'], 'fvSolution'), 'w') as fvsolfile:
            fvsolfile.write(fvsolwrite(expnode))
        with open(os.path.join(scene['flparams']['ofsfilebase'], 'fvSchemes'), 'w') as fvschfile:
            fvschfile.write(fvschwrite(expnode))
        with open(os.path.join(scene['flparams']['ofcfilebase'], 'transportProperties'), 'w') as fvtppfile:
            fvtppfile.write(fvtppwrite(expnode.solver))
        if expnode.solver != 'icoFoam':
            with open(os.path.join(scene['flparams']['ofcfilebase'], 'turbulenceProperties'), 'w') as fvrasfile:
                fvrasfile.write(fvraswrite(expnode.turbulence))
            if expnode.buoyancy:
                with open(os.path.join(scene['flparams']['ofcfilebase'], 'g'), 'w') as fvgrafile:
                    fvgrafile.write(fvgrawrite())
            if expnode.radiation:
                with open(os.path.join(scene['flparams']['ofcfilebase'], 'radiationProperties'), 'w') as fvradfile:
                    fvradfile.write(fvradrite())
        
#        with open(self.fpfile, 'w') as fvprogress:
#            if expnode.processes > 1:
#                with open(os.path.join(scene['flparams']['ofsfilebase'], 'decomposeParDict'), 'w') as fvdcpfile:
#                    fvdcpfile.write(fvdcpwrite(expnode))
#                call(("decomposePar -force -case {}".format(scene['flparams']['offilebase']).split()))
#                self.run = Popen('mpirun -np {} {} -parallel -case {}'.format(expnode.processes, expnode.solver, scene['flparams']['offilebase']).split(), stdout = fvprogress)
#            else:
#                self.run = Popen((expnode.solver, "-case", "{}".format(scene['flparams']['offilebase'])), stdout = fvprogress)
#        
        self.convergence = expnode.convergence
        self.econvergence = expnode.econvergence
        self.pv = expnode.pv
            
        expnode.postexport()
        
        return {'FINISHED'}

class NODE_OT_FVSolve(bpy.types.Operator):
    bl_idname = "node.fvsolve"
    bl_label = "FloVi simulation"
    bl_description = "Solve an OpenFOAM case"
    bl_register = True
    bl_undo = True
    nodeid = bpy.props.StringProperty()
    
    def modal(self, context, event):
        if self.run.poll() is None and self.kivyrun.poll() is None:
            with open(self.fpfile, 'r') as fpfile:
                lines = fpfile.readlines()[::-1]
                residict = {}

                for line in lines:
                    if 'Solving for' in line:
                        residict[line.split()[3][:-1]] = float(line.split()[7][:-1])
                    if len(residict) == len(self.residuals):
                        break

                for var in residict:
                    if var == 'e':
                        residict[var] -= self.econvergence
                    else:
                        residict[var] -= self.convergence

                if residict:
                    self.pfile.check("\n".join(['{0[0]} {0[1]}'.format(i) for i in residict.items()]))
            return {'PASS_THROUGH'}
        elif self.run.poll() is None and self.kivyrun.poll() is not None:
            self.run.kill()
            return {'CANCELLED'}
        else:
            self.kivyrun.kill()
            if self.processes > 1:
                Popen(("reconstructPar", "-case", "{}".format(context.scene['flparams']['offilebase'])))
                
            Popen(("postProcess", "-func", "writeCellCentres", "-case", "{}".format(context.scene['flparams']['offilebase'])))
#            if self.pv:
#                Popen(("paraFoam", "-case", "{}".format(context.scene['flparams']['offilebase'])))
            return {'FINISHED'}

    def invoke(self, context, event):
        wm = context.window_manager
        scene = context.scene
        simnode = bpy.data.node_groups[self.nodeid.split('@')[1]].nodes[self.nodeid.split('@')[0]]
        (self.convergence, self.econvergence, self.residuals, self.processes, self.solver)  = simnode.presim()
        self.fpfile = os.path.join(scene['viparams']['newdir'], 'floviprogress')
        self.pfile = fvprogressfile(scene['viparams']['newdir'])
#        self.pfile = progressfile(scene['viparams']['newdir'], datetime.datetime.now(), 100)
#        self.pfile = progressfile(self.folder, datetime.datetime.now(), 100)
        self.kivyrun = fvprogressbar(os.path.join(scene['viparams']['newdir'], 'viprogress'), str(self.residuals))

        with open(self.fpfile, 'w') as fvprogress:
            if self.processes > 1:
                with open(os.path.join(scene['flparams']['ofsfilebase'], 'decomposeParDict'), 'w') as fvdcpfile:
                    fvdcpfile.write(fvdcpwrite(self.processes))
                call(("decomposePar -force -case {}".format(scene['flparams']['offilebase']).split()))
                self.run = Popen('mpirun --oversubscribe -np {} {} -parallel -case {}'.format(self.processes, self.solver, scene['flparams']['offilebase']).split(), stdout = fvprogress)
            else:
                self.run = Popen((self.solver, "-case", "{}".format(scene['flparams']['offilebase'])), stdout = fvprogress)
        

#        self.pv = simnode.pv
            
#        simnode.postsim()
        self._timer = wm.event_timer_add(5, context.window)
        wm.modal_handler_add(self)        
        return {'RUNNING_MODAL'}
    
    def terminate(self, scene):
        self.run.kill()

class OBJECT_OT_VIGridify(bpy.types.Operator):
    ''''''
    bl_idname = "object.vi_gridify"
    bl_label = "VI Gridify"
     
    def modal(self, context, event):
        scene = context.scene
        if self.rotate != scene.vi_gridify_rot or self.us != scene.vi_gridify_us or self.acs != context.scene.vi_gridify_as or self.ft:
            self.bmnew = self.bm.copy()
            self.bmnew.transform(self.o.matrix_world)
            self.ft = 0
            self.us = context.scene.vi_gridify_us
            self.acs = context.scene.vi_gridify_as
            self.rotate = context.scene.vi_gridify_rot
            self.bmnew.faces.ensure_lookup_table()
            self.bmnew.verts.ensure_lookup_table()
            self.upv = self.bmnew.faces[0].calc_tangent_edge_pair().copy().normalized()
            self.norm = self.bmnew.faces[0].normal.copy()
            self.acv = self.upv.copy()
            eul = Euler(radians(-90) * self.norm, 'XYZ')
            self.acv.rotate(eul)
            rotation = Euler(radians(self.rotate) * self.norm, 'XYZ')
            self.upv.rotate(rotation)
            self.acv.rotate(rotation)
            vertdots = [Vector.dot(self.upv, vert.co) for vert in self.bmnew.verts]
            vertdots2 = [Vector.dot(self.acv, vert.co) for vert in self.bmnew.verts]
            svpos = self.bmnew.verts[vertdots.index(min(vertdots))].co
            svpos2 = self.bmnew.verts[vertdots2.index(min(vertdots2))].co
            res1, res2, ngs1, ngs2, gs1, gs2 = 1, 1, self.us, self.acs, self.us, self.acs
            vs = self.bmnew.verts[:]
            es = self.bmnew.edges[:]
            fs = [f for f in self.bmnew.faces[:]]
            gs = vs + es + fs
              
            while res1:
                res = bmesh.ops.bisect_plane(self.bmnew, geom = gs, dist = 0.001, plane_co = svpos + ngs1 * self.upv, plane_no = self.upv, use_snap_center = 0, clear_outer = 0, clear_inner = 0)
                res1 = res['geom_cut']
                gs = self.bmnew.verts[:] + self.bmnew.edges[:] + [v for v in res['geom'] if isinstance(v, bmesh.types.BMFace)]
                ngs1 += gs1
        
            while res2:
                res = bmesh.ops.bisect_plane(self.bmnew, geom = gs, dist = 0.001, plane_co = svpos2 + ngs2 * self.acv, plane_no = self.acv, use_snap_center = 0, clear_outer = 0, clear_inner = 0)
                res2 = res['geom_cut']
                gs = self.bmnew.verts[:] + self.bmnew.edges[:] + [v for v in res['geom'] if isinstance(v, bmesh.types.BMFace)]
                ngs2 += gs2
            
            self.bmnew.transform(self.o.matrix_world.inverted())
            self.bmnew.to_mesh(self.o.data)
            self.bmnew.free()
            context.area.tag_redraw()
            return {'RUNNING_MODAL'}

        elif event.type == 'ESC':  
            self.bm.to_mesh(self.o.data)
            self.bm.free()
            self.bmnew.free()
            context.area.tag_redraw()
            return {'CANCELLED'}

        elif event.ctrl and event.type == 'RET':
            self.bmnew.free()
            self.bm.free()
            return {'FINISHED'}
            
        else:
            return {'PASS_THROUGH'}
     
    def invoke(self, context, event):
        scene = context.scene
        self.o = bpy.context.active_object
        self.bm = bmesh.new()
        tm = self.o.to_mesh(scene = scene, apply_modifiers = True, settings = 'PREVIEW')
        self.bm.from_mesh(tm)
        bpy.data.meshes.remove(tm)
        self.ft = 1
        self.rotate = scene.vi_gridify_rot
        self.us = scene.vi_gridify_us
        self.acs = scene.vi_gridify_as
        context.window_manager.modal_handler_add(self)
        return {'RUNNING_MODAL'}

class OBJECT_OT_VIGridify2(bpy.types.Operator):
    ''''''
    bl_idname = "object.vi_gridify2"
    bl_label = "VI Gridify"
    bl_options = {"REGISTER", 'UNDO'}
    
    rotate =  bpy.props.FloatProperty(name = 'Rotation', default = 0, min = 0, max = 360) 
    us =  bpy.props.FloatProperty(default = 0.6, min = 0.01) 
    acs =  bpy.props.FloatProperty(default = 0.6, min = 0.01) 
    
    @classmethod
    def poll(cls, context):
        obj = context.active_object
        return (obj and obj.type == 'MESH')
    
    def execute(self, context):
        self.o = bpy.context.active_object
        mesh = bmesh.from_edit_mesh(self.o.data)
        mesh.faces.ensure_lookup_table()
        mesh.verts.ensure_lookup_table()
        self.upv = mesh.faces[0].calc_tangent_edge_pair().copy().normalized()
        self.norm = mesh.faces[0].normal.copy()
        self.acv = self.upv.copy()
        eul = Euler(radians(-90) * self.norm, 'XYZ')
        self.acv.rotate(eul)
        rotation = Euler(radians(self.rotate) * self.norm, 'XYZ')
        self.upv.rotate(rotation)
        self.acv.rotate(rotation)
        vertdots = [Vector.dot(self.upv, vert.co) for vert in mesh.verts]
        vertdots2 = [Vector.dot(self.acv, vert.co) for vert in mesh.verts]
        svpos = mesh.verts[vertdots.index(min(vertdots))].co
        svpos2 = mesh.verts[vertdots2.index(min(vertdots2))].co
        res1, res2, ngs1, ngs2, gs1, gs2 = 1, 1, self.us, self.acs, self.us, self.acs
        vs = mesh.verts[:]
        es = mesh.edges[:]
        fs = [f for f in mesh.faces[:] if f.select]
        gs = vs + es + fs
          
        while res1:
            res = bmesh.ops.bisect_plane(mesh, geom = gs, dist = 0.001, plane_co = svpos + ngs1 * self.upv, plane_no = self.upv, use_snap_center = 0, clear_outer = 0, clear_inner = 0)
            res1 = res['geom_cut']
            gs = mesh.verts[:] + mesh.edges[:] + [v for v in res['geom'] if isinstance(v, bmesh.types.BMFace)]
            ngs1 += gs1
    
        while res2:
            res = bmesh.ops.bisect_plane(mesh, geom = gs, dist = 0.001, plane_co = svpos2 + ngs2 * self.acv, plane_no = self.acv, use_snap_center = 0, clear_outer = 0, clear_inner = 0)
            res2 = res['geom_cut']
            gs = mesh.verts[:] + mesh.edges[:] + [v for v in res['geom'] if isinstance(v, bmesh.types.BMFace)]
            ngs2 += gs2
        bmesh.update_edit_mesh(self.o.data)
        return {'FINISHED'}