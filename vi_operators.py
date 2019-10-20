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
from .vi_func import selobj, joinobj, solarPosition, viparams, compass, wind_compass, spfc, solarRiseSet
#from .flovi_func import fvcdwrite, fvbmwrite, fvblbmgen, fvvarwrite, fvsolwrite, fvschwrite, fvtppwrite, fvraswrite, fvshmwrite, fvmqwrite, fvsfewrite, fvobjwrite, fvdcpwrite
from .vi_func import ret_plt, blf_props, draw_index, viewdesc, ret_vp_loc, draw_time, logentry, rettree
from .vi_func import sunpath, windnum, wind_rose, create_coll, retobjs, progressfile, progressbar, chunks, clearlayers, clearscene
from .vi_display import wr_legend, wr_scatter, wr_table, wr_disp
#from .envi_func import processf, retenvires, envizres, envilres, recalculate_text
#from .vi_chart import chart_disp
from .vi_misc import wr_legend2, results_bar

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
            out vec4 FragColour;
            
            void main()
                {
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
            
        if scene.vi_params.vi_display == 0 or scene.vi_params['viparams']['vidisp'] != 'sp':
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
        spcoll = create_coll(context, 'SunPath')            
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
        scene.vi_params.vi_display = 1
        return {'RUNNING_MODAL'}
        
class NODE_OT_WindRose(bpy.types.Operator):
    bl_idname = "node.windrose"
    bl_label = "Wind Rose"
    bl_description = "Create a Wind Rose"
    bl_register = True
    bl_undo = True
#    nodeid = bpy.props.StringProperty()

    def invoke(self, context, event):
        scene = context.scene
        simnode = context.node
        wrcoll = create_coll(context, 'WindRoses')
        plt = ret_plt()
        
        if viparams(self, scene):
            return {'CANCELLED'}
        if not plt:
            self.report({'ERROR'},"There is something wrong with your matplotlib installation")
            return {'FINISHED'}

        simnode.export()
        locnode = simnode.inputs['Location in'].links[0].from_node
        scene.vi_params['viparams']['resnode'], scene.vi_params['viparams']['restree'] = simnode.name, simnode.id_data.name
        scene.vi_params['viparams']['vidisp'], scene.vi_params.vi_display = 'wr', 0
        context.scene.vi_params['viparams']['visimcontext'] = 'Wind'
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
        dfreq = histogram(awd, bins=dbinvals)[0]
        adfreq = histogram(cwd, bins=dbinvals)[0]
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
            
        plt.savefig(scene.vi_params['viparams']['newdir']+'/disp_wind.svg')
#        (wro, scale) = wind_rose(simnode['maxres'], scene['viparams']['newdir']+'/disp_wind.svg', simnode.wrtype, mcolors)
        wrme = bpy.data.meshes.new("Wind_rose")   
        wro = bpy.data.objects.new('Wind_rose', wrme) 
        
#        context.view_layer.objects.active = wro
        
        if wro.name not in wrcoll.objects:
            wrcoll.objects.link(wro)
            if wro.name in scene.collection.objects:
                scene.collection.objects.unlink(wro)
                
        selobj(context.view_layer, wro)       
        (wro, scale) = wind_rose(wro, simnode['maxres'], scene.vi_params['viparams']['newdir']+'/disp_wind.svg', simnode.wrtype, mcolors)
        wro = joinobj(context.view_layer, wro)        
        wro['maxres'], wro['minres'], wro['avres'], wro['nbins'], wro['VIType'] = max(aws), min(aws), sum(aws)/len(aws), len(sbinvals), 'Wind_Plane'
        simnode['maxfreq'] = 100*numpy.max(adfreq)/len(cwd)
        simnode['maxfreq'] = 100*numpy.max(adfreq)/len(vawd) if simnode.max_freq == '0' else simnode.max_freq_val

        windnum(simnode['maxfreq'], (0,0,0), scale, wind_compass((0,0,0), scale, wro, wro.data.materials['wr-000000']))
        plt.close()
        wro['table'] = array([["", 'Minimum', 'Average', 'Maximum'], ['Speed (m/s)', wro['minres'], '{:.1f}'.format(wro['avres']), wro['maxres']], ['Direction (\u00B0)', min(awd), '{:.1f}'.format(sum(awd)/len(awd)), max(awd)]])
        wro['ws'] = aws.reshape(len(doys), 24).T.tolist()
        wro['wd'] = awd.reshape(len(doys), 24).T.tolist()
        wro['days'] = array(doys, dtype = float)
        wro['hours'] = arange(1, 25, dtype = float)        
        wro['maxfreq'] = 100*numpy.max(dfreq)/len(awd)
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
        if context.scene.vi_params.vi_display == 0 or context.scene.vi_params['viparams']['vidisp'] != 'wrpanel' or 'Wind_Plane' not in [o['VIType'] for o in bpy.data.objects if o.get('VIType')]:
            bpy.types.SpaceView3D.draw_handler_remove(self._handle_wr_disp, 'WINDOW')
            context.area.tag_redraw()
            return {'CANCELLED'}           

        if event.type != 'INBETWEEN_MOUSEMOVE' and context.region and context.area.type == 'VIEW_3D' and context.region.type == 'WINDOW':            
            mx, my = event.mouse_region_x, event.mouse_region_y 
            
            # Legend routine 
            
            if self.legend.spos[0] < mx < self.legend.epos[0] and self.legend.spos[1] < my < self.legend.epos[1]:
                self.legend.hl = (0.5, 0.5, 0.5, 0.5)  
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
                self.legend.hl = (0.5, 0.5, 0.5, 0.5) 
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
        context.scene.vi_params.vi_display = 1
        context.scene.vi_params['viparams']['vidisp'] = 'wrpanel'
        simnode = bpy.data.node_groups[context.scene.vi_params['viparams']['restree']].nodes[context.scene.vi_params['viparams']['resnode']]
        self.legend = wr_legend([80, context.region.height - 40], context.region.width, context.region.height, 'legend.png', 150, 350)
        self.dhscatter = wr_scatter([160, context.region.height - 40], context.region.width, context.region.height, 'scat.png', 600, 400)
        self.table = wr_table([240, context.region.height - 40], context.region.width, context.region.height, 'table.png', 600, 150)       
        self.legend.update(context)
        self.dhscatter.update(context)
        self.table.update(context)
        self._handle_wr_disp = bpy.types.SpaceView3D.draw_handler_add(wr_disp, (self, context, simnode), 'WINDOW', 'POST_PIXEL')
        context.window_manager.modal_handler_add(self)
        return {'RUNNING_MODAL'}
    
    
class VIEW3D_OT_WRDisplay2(bpy.types.Operator):
    bl_idname = "view3d.wrdisplay2"
    bl_label = "Wind rose number display"
    bl_description = "Project the windrose numbers on to the viewport"
    bl_register = True
    bl_undo = False
    
    def invoke(self, context, event):   
        area = context.area
        context.scene.vi_params.vi_display = 1
        context.scene.vi_params['viparams']['vidisp'] = 'wr'
        self.results_bar = results_bar(('legend.png', 'table_new.png', 'scatter_new.png'), 300, area)
#        self.leg_image = bpy.data.images['legend.png']
        self.legend = wr_legend2(context, 'Speed (m/s)', [305, area.height - 80], area.width, area.height, 'legend.png', 125, 300)
#        self.table = wr_table([240, area.height - 40], area.width, area.height, 'table_new.png', 600, 150)
#        self.dhscatter = wr_scatter([160, area.height - 40], area.width, area.height, 'scatter_new.png', 600, 400)
        self.legend.update(context)
        self.height = area.height
#        self.create_batch(context.area)
        self.draw_handle_wrnum = bpy.types.SpaceView3D.draw_handler_add(self.draw_wrnum, (context, ), 'WINDOW', 'POST_PIXEL')
        self.leg_show = 1
        area.tag_redraw()
        context.window_manager.modal_handler_add(self)
        return {'RUNNING_MODAL'}
    
    def modal(self, context, event):
        scene = context.scene
        ah, aw = context.area.height, context.area.width
        redraw = 0
           
        if scene.vi_params.vi_display == 0 or scene.vi_params['viparams']['vidisp'] != 'wr':
            bpy.types.SpaceView3D.draw_handler_remove(self.draw_handle_wrnum, 'WINDOW')
            return {'CANCELLED'}

        if event.type != 'INBETWEEN_MOUSEMOVE' and context.region and context.area.type == 'VIEW_3D' and context.region.type == 'WINDOW':            
            mx, my = event.mouse_region_x, event.mouse_region_y 
            
            # Legend routine 
            if self.legend.ispos[0] < mx < self.legend.iepos[0] and self.legend.ah - 80 < my < self.legend.ah - 40:
                self.legend.hl = (0.8, 0.8, 0.8, 0.8) 
                redraw = 1
                if event.type == 'LEFTMOUSE':
                    if event.value == 'RELEASE':
                        self.legend.expand = 0 if self.legend.expand else 1
            
            elif abs(self.legend.lspos[0] - mx) < 10 and abs(self.legend.lepos[1] - my) < 10:
                self.legend.hl = (0.8, 0.8, 0.8, 0.8) 
                redraw = 1   
                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.legend.move = 1
                        self.legend.draw(ah, aw)
                        context.area.tag_redraw()
                        return {'RUNNING_MODAL'}
                    elif self.legend.move and event.value == 'RELEASE':
                        self.legend.move = 0                        
                        return {'RUNNING_MODAL'}
                
                elif event.type == 'ESC':
                    bpy.types.SpaceView3D.draw_handler_remove(self._handle_wr_disp, 'WINDOW')
                    redraw = 1
                    return {'CANCELLED'}
                                
            elif abs(self.legend.lepos[0] - mx) < 10 and abs(self.legend.lspos[1] - my) < 10:
                self.legend.hl = (0.8, 0.8, 0.8, 0.8) 
                context.area.tag_redraw()
                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.legend.resize = 1
                        self.legend.draw(ah, aw)
                        context.area.tag_redraw()
                    elif self.legend.resize and event.value == 'RELEASE':
                        self.legend.resize = 0
                    return {'RUNNING_MODAL'}
                    
            elif self.legend.hl == (0.8, 0.8, 0.8, 0.8):                 
                self.legend.hl = (1, 1, 1, 1)
                redraw = 1
                                
            # Move routines
                     
            if event.type == 'MOUSEMOVE':                
                if self.legend.move:
                    self.legend.lspos[0], self.legend.lepos[1] = mx, my
                    self.legend.draw(ah, aw)
                    context.area.tag_redraw() 
                elif self.legend.resize:
                    self.legend.lepos[0], self.legend.lspos[1] = mx, my
                    self.legend.draw(ah, aw)
                    context.area.tag_redraw() 
                
            if redraw:
                context.area.tag_redraw()        
        return {'PASS_THROUGH'}
    
    def ret_bar_geometry(self, area): 
#        bm = bmesh.new()
#        xi = (0, 1, 1, 0)
#        yi = (0, 0, 1, 1)
#        
#        for i in range(4):
#            pos = (150 + xi[i] * 100, 25 + yi[i] * 50)
#            bm.verts.new((pos[0], pos[1], 0))
#            
#        bm.faces.new(bm.verts)
#        bmesh.ops.bevel(bm, geom = bm.verts[:] + bm.faces[:], offset = 15, offset_type = 'OFFSET', segments = 10, profile = 0.5, vertex_only = True)
#        
#        bm.verts.ensure_lookup_table()
#        bm.faces.ensure_lookup_table()
#        bar_v_coords = [[v.co[:2] for v in edge.verts] for edge in bm.edges]
        
#        bar_v_coords = [item for sublist in bar_v_coords for item in sublist]
#        barl_v_coords = [bm.faces[0].loops[i].vert.co[:2] for i in range(len(bm.verts))]
#        e = bm.edges[0]
#        bar_v_coords = [e.verts[0].co[:2]]
#        for edge in bm.edges[1:]:
#            if e.verts[1] in edge.verts:
#                bar_v_coords.append()
#                edge = e
                
            
        bar_v_coords = ((300, area.height - 85), (450, area.height - 85),
                         (450, area.height - 35), (300, area.height - 35))
#        bmesh.ops.triangulate(bm, faces = bm.faces, quad_method = 'BEAUTY', ngon_method= 'BEAUTY')
#        bar_v_coords = [[v.co[:2] for v in face.verts] for face in bm.faces]
#        barf_v_coords = [item for sublist in bar_v_coords for item in sublist]
#        bar_f_indices = [[v.index for v in face.verts] for face in bm.faces]
#        bm.free()
        return (bar_v_coords, ((0, 1, 2), (2, 3, 0)))
        
    def create_batch(self, area):
        self.wrbl_shader = gpu.shader.from_builtin('2D_UNIFORM_COLOR')
        self.wrbf_shader = gpu.shader.from_builtin('2D_UNIFORM_COLOR')
        self.leg_shader = gpu.shader.from_builtin('2D_IMAGE')
        self.tab_shader = gpu.shader.from_builtin('2D_IMAGE')
        self.hm_shader = gpu.shader.from_builtin('2D_IMAGE')
        (v_coords, f_indices) = self.ret_bar_geometry(area)
        self.wrbl_batch = batch_for_shader(self.wrbl_shader, 'LINE_LOOP', {"pos": v_coords})
        self.wrbf_batch = batch_for_shader(self.wrbf_shader, 'TRIS', {"pos": v_coords}, indices = f_indices)
        leg_pos = ((305, area.height - 80), (345, area.height - 80),(345, area.height - 40), (305, area.height - 40))
        self.leg_batch = batch_for_shader(self.leg_shader, 'TRI_FAN', {"pos": leg_pos, "texCoord": ((0, 0), (1, 0), (1, 1), (0, 1)),},)
        tab_pos = ((355, area.height - 80), (395, area.height - 80),(395, area.height - 40), (355, area.height - 40))
        self.tab_batch = batch_for_shader(self.tab_shader, 'TRI_FAN', {"pos": tab_pos, "texCoord": ((0, 0), (1, 0), (1, 1), (0, 1)),},)
        hm_pos = ((405, area.height - 80), (445, area.height - 80),(445, area.height - 40), (405, area.height - 40))
        self.hm_batch = batch_for_shader(self.hm_shader, 'TRI_FAN', {"pos": hm_pos, "texCoord": ((0, 0), (1, 0), (1, 1), (0, 1)),},)
        
    def draw_wrnum(self, context):
        area = context.area
        ah = area.height
                
        self.results_bar.draw(ah, 3)
        self.legend.draw(ah, area.width)
        
        
#        if self.leg_show:
#            self.legend.drawopen(context) 

#        width, height = area.width, context.region.height
#        self.legend.draw(context, width, height)
#        self.dhscatter.draw(context, width, height)
#        self.table.draw(context, width, height)

class NODE_OT_SVF(bpy.types.Operator):
    bl_idname = "node.svf"
    bl_label = "Sky View Factor"
    bl_description = "Undertake a sky view factor study"
    bl_register = True
    bl_undo = False
#    nodeid = bpy.props.StringProperty()

    def invoke(self, context, event):
        scene = context.scene  
        scene.vi_params.vi_display = 0
        
        if viparams(self, scene):            
            return {'CANCELLED'}

        shadobs = retobjs('livig')
        if not shadobs:
            self.report({'ERROR'},"No shading objects have a material attached.")
            return {'CANCELLED'}
            
        scene.vi_params['liparams']['shadc'] = [ob.name for ob in retobjs('ssc')]
        if not scene.vi_params['liparams']['shadc']:
            self.report({'ERROR'},"No objects have a light sensor material attached.")
            return {'CANCELLED'}

        scene['viparams']['restree'] = self.nodeid.split('@')[1]
        clearscene(scene, self)
        simnode = context.node
        scene.vi_params['viparams']['visimcontext'] = 'SVF'

        if not scene.vi_params.get('liparams'):
           scene['liparams'] = {}
           
        scene.vi_params['liparams']['cp'], scene.vi_params['liparams']['unit'], scene.vi_params['liparams']['type'] = simnode.cpoint, '% Sunlit', 'VI Shadow'
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