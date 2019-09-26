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

import bpy, datetime, mathutils, os, bmesh, shutil, sys, math, shlex, gpu, bgl, blf
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
from bpy_extras import view3d_utils
#from multiprocessing import Pool
#from .livi_export import radgexport, spfc, createoconv, createradfile
#from .livi_calc  import li_calc
#from .vi_display import li_display, linumdisplay, spnumdisplay, en_air, wr_legend, wr_disp, wr_scatter, wr_table, ss_disp, ss_legend, svf_disp, svf_legend, basic_legend, basic_table, basic_disp, ss_scatter, en_disp, en_pdisp, en_scatter, en_table, en_barchart, comp_table, comp_disp, leed_scatter, cbdm_disp, cbdm_scatter, envals, bsdf, bsdf_disp#, en_barchart, li3D_legend
#from .envi_export import enpolymatexport, pregeo
#from .envi_mat import envi_materials, envi_constructions
from .vi_func import selobj, joinobj, solarPosition, viparams, compass, spfc, solarRiseSet
#from .flovi_func import fvcdwrite, fvbmwrite, fvblbmgen, fvvarwrite, fvsolwrite, fvschwrite, fvtppwrite, fvraswrite, fvshmwrite, fvmqwrite, fvsfewrite, fvobjwrite, fvdcpwrite
from .vi_func import spathrange, ret_plt, blf_props, draw_index, viewdesc, ret_vp_loc, draw_time
from .vi_func import sunpath
from .vi_display import spnumdisplay
#from .envi_func import processf, retenvires, envizres, envilres, recalculate_text
#from .vi_chart import chart_disp

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
#    if plt:
#        from .windrose import WindroseAxes

try:
    import psutil
    psu = 1
except: 
    psu = 0    


rvuerrdict = {'view up parallel to view direction': "Camera cannot point directly upwards", 
              ' x11': "No X11 display server found. You may need to install XQuartz", 
              'source center': "A light source has concave faces. Use mesh - cleanup - split concave faces"}
pmerrdict = {'fatal - too many prepasses, no global photons stored\n': "Too many prepasses have occurred. Make sure light sources can see your geometry",
             'fatal - too many prepasses, no global photons stored, no caustic photons stored\n': "Too many prepasses have occurred. Turn off caustic photons and encompass the scene",
               'fatal - zero flux from light sources\n': "No light flux, make sure there is a light source and that photon port normals point inwards",
               'fatal - no light sources in distribPhotons\n': "No light sources. Photon mapping does not work with HDR skies",
               'fatal - no valid photon ports found\n': 'Make sure photon ports are valid', 
               'fatal - failed photon distribution\n': 'Do the lights see enough geometry?'}
    
class NODE_OT_SunPath(bpy.types.Operator):
    bl_idname = "node.sunpath"
    bl_label = "Sun Path"
    bl_description = "Create a Sun Path"
    bl_register = True
    bl_undo = False
    
    def ret_coords(self, scene, node):
        breaks, coords, sd, d, line_lengths, sumcoords, wincoords = [], [], 100, 0, [0], [], []
                
        for hour in range(1, 25):
            for doy in range(0, 365, 2):
                ([solalt, solazi]) = solarPosition(doy, hour, scene.vi_params.latitude, scene.vi_params.longitude)[2:]
                coord = Vector([-(sd-(sd-(sd*cos(solalt))))*sin(solazi), -(sd-(sd-(sd*cos(solalt))))*cos(solazi), sd*sin(solalt)])
                coords.append(coord)
                if d%183 == 0:
                    breaks.append(1)
                else:
                    breaks.append(0)
                d += 1

        for doy in (79, 172, 355):
            for hour in range(1, 241):
                ([solalt, solazi]) = solarPosition(doy, hour*0.1, scene.vi_params.latitude, scene.vi_params.longitude)[2:]
                coord = Vector([-(sd-(sd-(sd*cos(solalt))))*sin(solazi), -(sd-(sd-(sd*cos(solalt))))*cos(solazi), sd*sin(solalt)])
                coords.append(coord)
                breaks.append(2)
                
                if doy == 172:
                    if hour in (120, 240):
                        sumcoords.append(coord)
                elif doy == 355:
                    if hour in (120, 240):
                        wincoords.append(coord)
                    
        self.summid = (sumcoords[0]+sumcoords[1])/2
        self.winmid = (wincoords[0]+wincoords[1])/2               
        self.sumnorm = mathutils.Matrix().Rotation(math.pi/2, 4, 'X')@mathutils.Vector([0] + list((sumcoords[0]-sumcoords[1])[1:])).normalized()   
        self.winnorm = mathutils.Matrix().Rotation(math.pi/2, 4, 'X')@mathutils.Vector([0] + list((wincoords[0]-wincoords[1])[1:])).normalized()    
                
        for a, b in zip(coords[:-1], coords[1:]):
            line_lengths.append(line_lengths[-1] + (a - b).length)
            
        return(coords, line_lengths, breaks)
        
    def create_batch(self, scene, node):
        sp_vertex_shader = '''
            uniform mat4 viewProjectionMatrix;
            uniform mat4 sp_matrix;
            uniform vec4 colour1;
            uniform vec4 colour2;
            uniform vec4 colour3;
            in vec3 position;
            in float arcLength;
            in uint line_break;
            
            out vec4 v_colour1;
            out vec4 v_colour2;
            out vec4 v_colour3;
            out float v_ArcLength;
            out float zpos;
            flat out uint lb;
            
            void main()
            {
                v_colour1 = colour1;
                v_colour2 = colour2;
                v_colour3 = colour3;
                v_ArcLength = arcLength;
                gl_Position = viewProjectionMatrix * sp_matrix * vec4(position, 1.0f);
                zpos = vec3(position)[2];
                lb = line_break;
            }
        '''
        
        sp_fragment_shader = '''
            uniform float dash_density;
            uniform float dash_ratio;
            in float zpos;
            in vec4 v_colour1;
            in vec4 v_colour2;
            in vec4 v_colour3;
            in float v_ArcLength;
            flat in uint lb;
            out vec4 FragColour;
 
            void main()
            {
                if (zpos < 0) {discard;}
                else if (lb == uint(1)) {discard;}
                else if (sin(v_ArcLength * dash_density) > dash_ratio) {FragColour = v_colour1;} else {FragColour = v_colour2;}
                if (lb == uint(2)) {FragColour = v_colour3;}
            }
        '''
        
        sun_vertex_shader = '''
            uniform mat4 viewProjectionMatrix;
            uniform mat4 sp_matrix;
            in vec3 position;
            
            void main()
                {
                    gl_Position = viewProjectionMatrix * sp_matrix * vec4(position, 1.0f);
                    gl_Position[2] -= 0.001;
                }
            '''
        sun_fragment_shader = '''
            uniform vec4 sun_colour;
            out vec4 FragColour;
            
            void main()
                {
                    vec2 pos = gl_PointCoord - vec2(0.5);
                    if (length(pos) < 0.4) {FragColour = sun_colour;}
                    if (length(pos) <= 0.5) {
//                            sun_colour[3] = (0.5 - length(pos)) * 10;
                            FragColour = sun_colour;
                            FragColour[3] = (0.5 - length(pos)) * 10;
                        }
                    if (length(pos) > 0.5) {discard;}
                    
                }            
            '''
        sun2_vertex_shader = '''
            uniform mat4 viewProjectionMatrix;
            uniform mat4 sp_matrix;
            in vec3 position;
            in vec3 sun_position;
//            uniform float sun_radius;
//            out vec4 sun_position;
//            out mat4 sp_matrix;
//            out mat4 viewProjectionMatrix;
            out float sun_alpha;
            out vec4 sun_Position;
            out vec3 sp;
            
            void main()
                {
                    gl_Position = viewProjectionMatrix * sp_matrix * vec4(position, 1.0f);
                    sun_Position = viewProjectionMatrix * sp_matrix * vec4(sun_position, 1.0f); 
//                    sun_Position = ftransform(vec4(sun_position, 1.0f));
//                    sun_alpha = length(gl_Position.xy - sun_Position.xy) * 0.1;
//                    sp = sun_position;
                }
            '''
        sun2_fragment_shader = '''
            uniform vec4 sun_colour;
            uniform vec4 viewport;
            out vec4 FragColour;
            in vec4 sun_Position;
            
            void main()
                {
                    vec3 ndc = sun_Position.xyz/sun_Position.w;
                    vec2 vp_coord = ndc.xy * 0.5 + 0.5;
                    vec2 vpp_coord = vp_coord * viewport.zw;
                    float pos = length(gl_FragCoord.xy - vpp_coord);
                    float radius = sun_Position.z * 100;
                    if (pos > radius) 
                        {discard;}
                    if (pos >= radius - 10) 
                        {
                            FragColour = sun_colour;
                            FragColour[3] = 0.1*(radius - 10 - pos);
                        }
                    if (pos < radius -10) 
                        {FragColour = sun_colour;}
                }
           
            '''
            
        globe_vertex_shader = '''
            uniform mat4 viewProjectionMatrix;
            uniform mat4 sp_matrix;
            in vec3 position;
            
            void main()
                {
                    gl_Position = viewProjectionMatrix * sp_matrix * vec4(position, 1.0f);
                }
            '''
        globe_fragment_shader = '''
            uniform vec4 colour;
            out vec4 FragColour;
            
            void main()
                {
                    FragColour = colour;
                }
           
            '''
            
        range_vertex_shader = '''
            uniform mat4 viewProjectionMatrix;
            uniform mat4 sp_matrix;
            in vec3 position;
            in vec3 colour;
            out vec3 tri_colour;
            
            void main()
                {
                    gl_Position = viewProjectionMatrix * sp_matrix * vec4(position, 1.0f);
                    tri_colour = colour;
                }
            '''
        range_fragment_shader = '''
            in vec3 tri_colour;
//            float colours[9] = float[9](1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
//            vec4 green = vec4(0.0, 1.0, 0.0, 1.0);
//            vec4 blue = vec4(0.0, 0.0, 1.0, 1.0);
//            vec3 colours = vec3(red, green, blue);
            out vec4 FragColour;
            
            void main()
                {
//                    FragColour = vec4(colours[tri_colour*3], colours[tri_colour*3 + 1], colours[tri_colour*3 + 2], 1.0);
                        FragColour = vec4(tri_colour, 1.0);
                }
           
            '''
            
        self.sp_shader = gpu.types.GPUShader(sp_vertex_shader, sp_fragment_shader) 
        self.sun_shader = gpu.types.GPUShader(sun_vertex_shader, sun_fragment_shader)  
        self.globe_shader = gpu.types.GPUShader(globe_vertex_shader, globe_fragment_shader) 
        self.range_shader = gpu.types.GPUShader(range_vertex_shader, range_fragment_shader)
        (coords, line_lengths, breaks) = self.ret_coords(scene, node)
        sun_pos = [so.location[:] for so in scene.objects if so.type == 'LIGHT' and so.data.type == 'SUN' and not so.hide_viewport]
        sun_v_coords, sun_f_indices = self.ret_sun_geometry(scene.vi_params.sp_sun_size, self.suns)
        globe_v_coords, globe_f_indices = self.ret_globe_geometry(self.latitude, self.longitude)
        range_v_coords, range_f_indices, range_col_indices = self.ret_range_geometry(self.latitude, self.longitude)
        self.sp_batch = batch_for_shader(self.sp_shader, 'LINE_STRIP', {"position": coords, "arcLength": line_lengths, "line_break": breaks})
        self.sun_batch = batch_for_shader(self.sun_shader, 'POINTS', {"position": sun_pos})
        self.globe_batch = batch_for_shader(self.globe_shader, 'TRIS', {"position": globe_v_coords}, indices=globe_f_indices)
        self.range_batch = batch_for_shader(self.range_shader, 'TRIS', {"position": range_v_coords, "colour": range_col_indices})
        
    def draw_sp(self, op, context, node):
#        context.scene.vi_params.latitude = context.scene.vi_params.latitude
        # Draw lines
        bgl.glEnable(bgl.GL_DEPTH_TEST)
        bgl.glDepthFunc(bgl.GL_LESS)
        bgl.glDepthMask(bgl.GL_FALSE)
        bgl.glEnable(bgl.GL_BLEND)
        
        
#        bgl.glBlendFunc(bgl.GL_SRC_ALPHA, bgl.GL_ONE_MINUS_SRC_ALPHA)
#        bgl.glBlendFunc(bgl.GL_SRC_ALPHA, bgl.GL_ONE_MINUS_SRC_ALPHA)
#        bgl.glBlendFunc(bgl.GL_SRC_ALPHA,bgl.GL_SRC_ALPHA)
#        bgl.glBlendFuncSeparate(bgl.GL_SRC_ALPHA, bgl.GL_ONE_MINUS_SRC_ALPHA, bgl.GL_SRC_ALPHA, bgl.GL_DST_ALPHA )
#        bgl.glEnable(bgl.GL_MULTISAMPLE)
#        bgl.glEnable(bgl.GL_LINE_SMOOTH)
#        bgl.glEnable(bgl.GL_CULL_FACE)
#        bgl.glCullFace(bgl.GL_BACK)
#        bgl.glEnable(bgl.GL_POLYGON_SMOOTH)
#        bgl.glHint(bgl.GL_POLYGON_SMOOTH_HINT, bgl.GL_NICEST)
        bgl.glLineWidth(context.scene.vi_params.sp_line_width)
        bgl.glPointSize(context.scene.vi_params.sp_sun_size)
        
        self.sp_shader.bind()
        matrix = bpy.context.region_data.perspective_matrix
        sp_matrix = context.scene.objects['SPathMesh'].matrix_world
        sun_pos = [so.location[:] for so in context.scene.objects if so.type == 'LIGHT' and so.data.type == 'SUN' and not so.hide_viewport]
        self.sp_shader.uniform_float("viewProjectionMatrix", matrix)
        self.sp_shader.uniform_float("sp_matrix", sp_matrix)
        self.sp_shader.uniform_float("colour1", context.scene.vi_params.sp_hour_dash)
        self.sp_shader.uniform_float("colour2", context.scene.vi_params.sp_hour_main)
        self.sp_shader.uniform_float("colour3", context.scene.vi_params.sp_season_main)
        self.sp_shader.uniform_float("dash_ratio", context.scene.vi_params.sp_hour_dash_ratio)
        self.sp_shader.uniform_float("dash_density", context.scene.vi_params.sp_hour_dash_density) 
        self.sun_shader.bind()
        self.sun_shader.uniform_float("viewProjectionMatrix", matrix)
        self.sun_shader.uniform_float("sp_matrix", sp_matrix)
        self.sun_shader.uniform_float("sun_colour", context.scene.vi_params.sp_sun_colour) 
        self.globe_shader.bind()
        self.globe_shader.uniform_float("viewProjectionMatrix", matrix)
        self.globe_shader.uniform_float("sp_matrix", sp_matrix)
        self.globe_shader.uniform_float("colour", context.scene.vi_params.sp_globe_colour) 
        self.range_shader.bind()
        self.range_shader.uniform_float("viewProjectionMatrix", matrix)
        self.range_shader.uniform_float("sp_matrix", sp_matrix)
        
        if self.latitude != context.scene.vi_params.latitude or self.longitude != context.scene.vi_params.longitude or \
            self.sd != context.scene.vi_params.sp_sd or self.sh != context.scene.vi_params.sp_sh or self.ss != context.scene.vi_params.sp_sun_size:
            (coords, line_lengths, breaks) = self.ret_coords(context.scene, node)        
            self.sp_batch = batch_for_shader(self.sp_shader, 'LINE_STRIP', {"position": coords, "arcLength": line_lengths, "line_break": breaks})
            sun_pos = [so.location[:] for so in context.scene.objects if so.type == 'LIGHT' and so.data.type == 'SUN' and not so.hide_viewport]
            self.sun_batch = batch_for_shader(self.sun_shader, 'POINTS', {"position": sun_pos})
            globe_v_coords, globe_f_indices = self.ret_globe_geometry(self.latitude, self.longitude)            
            self.globe_batch = batch_for_shader(self.globe_shader, 'TRIS', {"position": globe_v_coords}, indices=globe_f_indices)
            range_v_coords, range_f_indices, range_col_indices = self.ret_range_geometry(self.latitude, self.longitude)
            self.range_batch = batch_for_shader(self.range_shader, 'TRIS', {"position": range_v_coords, "colour": range_col_indices})#, indices=range_f_indices)
            self.latitude = context.scene.vi_params.latitude
            self.longitude = context.scene.vi_params.longitude
            self.sd = context.scene.vi_params.sp_sd
            self.sh = context.scene.vi_params.sp_sh
            self.ss = context.scene.vi_params.sp_sun_size
            
        self.range_batch.draw(self.range_shader)    
        self.globe_batch.draw(self.globe_shader)
        bgl.glEnable(bgl.GL_LINE_SMOOTH)
        bgl.glEnable(bgl.GL_MULTISAMPLE)
        self.sun_batch.draw(self.sun_shader)
        self.sp_batch.draw(self.sp_shader)
        bgl.glDisable(bgl.GL_MULTISAMPLE)
        bgl.glDisable(bgl.GL_LINE_SMOOTH)
        
        bgl.glDisable(bgl.GL_BLEND)
        bgl.glClear(bgl.GL_DEPTH_BUFFER_BIT)
        bgl.glDisable(bgl.GL_DEPTH_TEST) 
        bgl.glDepthMask(bgl.GL_TRUE)
        
        
#        bgl.glEnable(bgl.GL_LINE_SMOOTH)
#        bgl.glEnable(bgl.GL_MULTISAMPLE)
#        bgl.glEnable(bgl.GL_POLYGON_SMOOTH)
#        bgl.glHint(bgl.GL_POLYGON_SMOOTH_HINT, bgl.GL_NICEST)
#        bgl.glDisable(bgl.GL_CULL_FACE)
        
#        bgl.glDisable(bgl.GL_LINE_SMOOTH)
#        bgl.glEnable(bgl.GL_BLEND)
        
#        bgl.glDisable(bgl.GL_BLEND)
#        bgl.glDisable(bgl.GL_POLYGON_SMOOTH)
        bgl.glPointSize(1)
        
    def draw_spnum(self, op, context):
        scene = context.scene
    
        if bpy.data.objects.get('SPathMesh'):
            spob = bpy.data.objects['SPathMesh'] 
            ob_mat = spob.matrix_world
            mid_x, mid_y, width, height = viewdesc(context)
            vl = ret_vp_loc(context)
            blf_props(scene, width, height)
            
            if scene.vi_params.sp_hd:
                coords, scene, sd = {}, context.scene, 100
                dists, hs, pos = [], [], []
                
                for doy in (172, 355):
                    for hour in range(24):
                        ([solalt, solazi]) = solarPosition(doy, hour, scene.vi_params.latitude, scene.vi_params.longitude)[2:]
                        if solalt > 0:
                            coords['{}-{}'.format(str(doy), str(hour))] = (Vector([-(sd-(sd-(sd*cos(solalt))))*sin(solazi), 
                                                                                   -(sd-(sd-(sd*cos(solalt))))*cos(solazi), 
                                                                                   sd*sin(solalt)])) 

                for co in coords:
                    dists.append((vl - coords[co]).length) 
                    hs.append(int(co.split('-')[1]))
                    pos.append(view3d_utils.location_3d_to_region_2d(context.region, 
                                                                     context.region_data, 
                                                                     ob_mat@coords[co]))
                
                if pos:
                    draw_index(pos, hs, dists, scene.vi_params.display_rp_fs, scene.vi_params.display_rp_fc, scene.vi_params.display_rp_fsh)
                    
            if [ob.get('VIType') == 'Sun' for ob in bpy.data.objects] and scene.vi_params['spparams']['suns'] == '0':
                sobs = [ob for ob in bpy.data.objects if ob.get('VIType') == 'Sun']
                
                if sobs and scene.vi_params.sp_td:
                    sunloc = ob_mat@sobs[0].location
                    solpos = view3d_utils.location_3d_to_region_2d(context.region, context.region_data, sunloc)
                    
                    try:
                        if 0 < solpos[0] < width and 0 < solpos[1] < height and not scene.ray_cast(context.view_layer, sobs[0].location + 0.05 * (vl - sunloc), vl - sunloc)[0]:
                            soltime = datetime.datetime.fromordinal(scene.vi_params.sp_sd)
                            soltime += datetime.timedelta(hours = scene.vi_params.sp_sh)
                            sre = sobs[0].rotation_euler
                            blf_props(scene, width, height)
                            sol_text = soltime.strftime('     %d %b %X') + ' alt: {:.1f} azi: {:.1f}'.format(90 - sre[0]*180/pi, (180, -180)[sre[2] < -pi] - sre[2]*180/pi)
                            draw_time(solpos, sol_text, scene.vi_params.display_rp_fs, 
                                      scene.vi_params.display_rp_fc, scene.vi_params.display_rp_fsh)
                            
                    except Exception as e:
                        print(e)
            blf.disable(0, 4)
        else:
            return
        
    def modal(self, context, event):
        scene = context.scene
       
        if context.area:
            context.area.tag_redraw()
            
        if scene.vi_display == 0 or scene.vi_params['viparams']['vidisp'] != 'sp':
            bpy.types.SpaceView3D.draw_handler_remove(self.draw_handle_sp, "WINDOW")
            bpy.types.SpaceView3D.draw_handler_remove(self.draw_handle_spnum, 'WINDOW')
            [bpy.data.objects.remove(o, do_unlink=True, do_id_user=True, do_ui_user=True) for o in bpy.data.objects if o.get('VIType') and o['VIType'] in ('SunMesh', 'SkyMesh')]
            return {'CANCELLED'}
        return {'PASS_THROUGH'}
    
    def ret_sun_geometry(self, dia, suns):
        sun_v_coords, sun_f_indices = [], []
        
        for sun in suns:
            sunbm = bmesh.new()
            bmesh.ops.create_uvsphere(sunbm, u_segments = 12, v_segments = 12, diameter = dia, matrix = sun.matrix_world, calc_uvs = 0)
            bmesh.ops.triangulate(sunbm, faces = sunbm.faces, quad_method = 'BEAUTY', ngon_method = 'BEAUTY')
            sun_v_coords += [v.co[:] for v in sunbm.verts]
            sun_f_indices += [[v.index for v in face.verts] for face in sunbm.faces]
            sunbm.free()
        return sun_v_coords, sun_f_indices

    def ret_globe_geometry(self, lat, long):
        midalt = solarPosition(79, 12, lat, long)[2]
        globebm = bmesh.new()
        altrot = mathutils.Matrix().Rotation(midalt, 4, 'X')
        bmesh.ops.create_uvsphere(globebm, u_segments = 48, v_segments = 48, diameter = 100, 
                                  matrix = altrot, calc_uvs = 0)        
        bmesh.ops.bisect_plane(globebm, geom = globebm.verts[:] + globebm.edges[:] + globebm.faces[:], dist = 0.01, 
                               plane_co = (0, 0, 0), plane_no = (0, 0, 1), use_snap_center = False, clear_outer = False, clear_inner = True)
        bmesh.ops.bisect_plane(globebm, geom = globebm.verts[:] + globebm.edges[:] + globebm.faces[:], dist = 0.1, 
                               plane_co = self.winmid, plane_no = self.winnorm, use_snap_center = False, 
                               clear_outer = True, clear_inner = False)
        bmesh.ops.bisect_plane(globebm, geom = globebm.verts[:] + globebm.edges[:] + globebm.faces[:], dist = 0.1, 
                               plane_co = self.summid, plane_no = self.sumnorm, use_snap_center = False, 
                               clear_outer = False, clear_inner = True)
        bmesh.ops.triangulate(globebm, faces = globebm.faces, quad_method = 'BEAUTY', ngon_method = 'BEAUTY')
        globe_v_coords = [v.co[:] for v in globebm.verts]
        globe_f_indices = [[v.index for v in face.verts] for face in globebm.faces]
        globebm.free()
        return globe_v_coords, globe_f_indices
    
    def ret_range_geometry(self, lat, long):
        bm = bmesh.new()
        params = ((177, 0.2, 0), (80, 0.25, 1), (355, 0.3, 2)) if lat >= 0 else ((355, 0.2, 0), (80, 0.25, 1), (177, 0.3, 2))
        v1s = range(1, 159)
        v2s = range(2, 160)
        v3s = range(159, 0, -1)
        
        for param in params:
            morn = solarRiseSet(param[0], 0, lat, long, 'morn')
            eve = solarRiseSet(param[0], 0, lat, long, 'eve')
            if morn or param[2] == 0:
                if morn:
                    mornevediff = eve - morn if lat >= 0 else 360 - eve + morn
                else:
                    mornevediff = 360# if bpy.context.scene.latitude >= 0 else 360
                
                startset = morn if lat >= 0 else eve
                angrange = [startset + a * 0.0125 * mornevediff for a in range (0, 81)]
                bm.verts.new().co = (97*sin(angrange[0]*pi/180), 97*cos(angrange[0]*pi/180), param[1])
            
                for a in angrange[1:-1]:
                    bm.verts.new().co = (95*sin(a*pi/180), 95*cos(a*pi/180), param[1])
                
                bm.verts.new().co = (97*sin(angrange[len(angrange) - 1]*pi/180), 97*cos(angrange[len(angrange) - 1]*pi/180), param[1])
                angrange.reverse()
        
                for b in angrange[1:-1]:
                    bm.verts.new().co = (99*sin(b*pi/180), 99*cos(b*pi/180), param[1])

                bm.verts.ensure_lookup_table()
                face = bm.faces.new((bm.verts[0 + 160 * param[2]], bm.verts[1 + 160 * param[2]], bm.verts[159 + 160 * param[2]]))
                face.material_index = param[2]
                   
                for i in range(79):
                    j = i + 80
                    face = bm.faces.new((bm.verts[v1s[i] + 160 * param[2]], 
                                         bm.verts[v2s[i] + 160 * param[2]], bm.verts[v3s[i] + 160 * param[2]]))
                    face.material_index = param[2]

                    if j < 158:
                        face = bm.faces.new((bm.verts[v1s[j] + 160 * param[2]], 
                                             bm.verts[v2s[j] + 160 * param[2]], bm.verts[v3s[j] + 160 * param[2]]))
                        face.material_index = param[2]
                
        range_v_coords = [[v.co[:] for v in face.verts] for face in bm.faces]
        range_v_coords = [item for sublist in range_v_coords for item in sublist]
        range_f_indices = [[v.index for v in face.verts] for face in bm.faces]        
        range_col_indices = [[((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))[face.material_index] for v in face.verts] for face in bm.faces]
        range_col_indices = [item for sublist in range_col_indices for item in sublist]
        bm.free()
        return (range_v_coords, range_f_indices, range_col_indices)
        
    def invoke(self, context, event):        
        scene = context.scene
        scene.vi_params['viparams'] = {}
        scene.vi_params['spparams'] = {}
#        viparams(self, scene)
#            self.report({'ERROR'},"Save the Blender file before continuing")
#            return {'CANCELLED'}
        
        try:
            spcoll = bpy.data.collections['SunPath']
        except:
            spcoll = bpy.data.collections.new('SunPath')
            context.scene.collection.children.link(spcoll)
            
        for lcc in context.view_layer.layer_collection.children:
            if lcc.name == 'SunPath':
                context.view_layer.active_layer_collection = lcc
            
        sd = 100
        node = context.node
        # Set the node colour
        node.export()
        
        scene.vi_params['viparams']['resnode'], scene.vi_params['viparams']['restree'] = node.name, node.id_data.name
        scene.cursor.location = (0.0, 0.0, 0.0)
        suns = [ob for ob in scene.objects if ob.type == 'LIGHT' and ob.data.type == 'SUN' and ob.parent.get('VIType') == "SPathMesh" ]
        requiredsuns = {'0': 1, '1': 12, '2': 24}[node.suns]
        matdict = {'SPBase': (0, 0, 0, 1), 'SPPlat': (1, 1, 1, 1)}
        
        for mat in [mat for mat in matdict if mat not in bpy.data.materials]:
            bpy.data.materials.new(mat)
            bpy.data.materials[mat].diffuse_color = matdict[mat][:4]
            bpy.data.materials[mat].use_nodes = True
            nodes = bpy.data.materials[mat].node_tree.nodes

            for n in nodes:
                nodes.remove(n)
            if mat == 'SPPlat':
                node_material = nodes.new(type='ShaderNodeBsdfDiffuse')
                node_material.inputs[0].default_value = matdict[mat]
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

            suns = [ob for ob in context.scene.objects if ob.type == 'LIGHT' and ob.data.type == 'SUN' and ob.parent.get('VIType') == "SPathMesh"]            
            [sun.animation_data_clear() for sun in suns]

        if not suns or len(suns) < requiredsuns: 
            for rs in range(requiredsuns - len(suns)):
                bpy.ops.object.light_add(type='SUN', radius=1, align='WORLD', location=(0, 0, 0))
                suns.append(context.active_object)
       
        if scene.render.engine == 'CYCLES' and scene.world.get('node_tree') and 'Sky Texture' in [no.bl_label for no in scene.world.node_tree.nodes]:
            scene.world.node_tree.animation_data_clear()    
        
        if bpy.context.active_object and not bpy.context.active_object.hide_viewport:
            if bpy.context.active_object.type == 'MESH':
                bpy.ops.object.mode_set(mode = 'OBJECT')
        
        if any(ob.get('VIType') == "SPathMesh" for ob in context.scene.objects):
            spathob = [ob for ob in context.scene.objects if ob.get('VIType') == "SPathMesh"][0]
        else:
            bpy.ops.object.add(type = "MESH")
            spathob = context.active_object
        
            if spathob.name not in spcoll.objects:
                spcoll.objects.link(spathob)
                if spathob.name in scene.collection.objects:
                    scene.collection.objects.unlink(spathob)
    
            spathob.location, spathob.name,  spathob['VIType'] = (0, 0, 0), "SPathMesh", "SPathMesh"
            selobj(context.view_layer, spathob)
            compassos = compass((0,0,0.01), sd, spathob, bpy.data.materials['SPPlat'], bpy.data.materials['SPBase'])
            joinobj(context.view_layer, [compassos] + [spathob])
            spathob.cycles_visibility.diffuse, spathob.cycles_visibility.shadow, spathob.cycles_visibility.glossy, spathob.cycles_visibility.transmission, spathob.cycles_visibility.scatter = [False] * 5
            spathob.show_transparent = True
        
        for s, sun in enumerate(suns):
            if sun.name not in spcoll.objects:
                spcoll.objects.link(sun)
                scene.collection.objects.unlink(sun)
                
            sun.data.shadow_soft_size = 0.01            
            sun['VIType'] = 'Sun'
            sun['solhour'], sun['solday'] = scene.solhour, scene.solday
            sun.name = sun.data.name ='Sun{}'.format(s)
            sun.parent = spathob



        if spfc not in bpy.app.handlers.frame_change_post:
            bpy.app.handlers.frame_change_post.append(spfc)

        scene.vi_params['viparams']['vidisp'] = 'sp'
        scene.vi_params['spparams']['suns'] = node.suns
        scene.vi_params['viparams']['visimcontext'] = 'SunPath'
        sunpath(scene)
        node = context.node
        self.suns = [sun for sun in scene.objects if sun.type == "LIGHT" and sun.data.type == 'SUN']
        self.sp = scene.objects['SPathMesh']
        self.latitude = scene.vi_params.latitude
        self.longitude = scene.vi_params.longitude
        self.sd = scene.vi_params.sp_sd
        self.sh = scene.vi_params.sp_sh
        self.ss = scene.vi_params.sp_sun_size
        self.create_batch(scene, node)
        self.draw_handle_sp = bpy.types.SpaceView3D.draw_handler_add(self.draw_sp, (self, context, node), "WINDOW", "POST_VIEW")
        self.draw_handle_spnum = bpy.types.SpaceView3D.draw_handler_add(self.draw_spnum, (self, context), 'WINDOW', 'POST_PIXEL')
        context.window_manager.modal_handler_add(self)
        scene.vi_display = 1
        return {'RUNNING_MODAL'}
        
