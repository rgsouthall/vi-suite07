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
from .livi_export import radgexport, spfc, createoconv, createradfile
from .livi_calc  import li_calc
from .envi_export import enpolymatexport, pregeo
#from .envi_mat import envi_materials, envi_constructions

from .vi_func import selobj, joinobj, solarPosition, viparams, compass, wind_compass, spfc, solarRiseSet, livisimacc, retpmap

#from .flovi_func import fvcdwrite, fvbmwrite, fvblbmgen, fvvarwrite, fvsolwrite, fvschwrite, fvtppwrite, fvraswrite, fvshmwrite, fvmqwrite, fvsfewrite, fvobjwrite, fvdcpwrite
from .vi_func import ret_plt, blf_props, draw_index, viewdesc, ret_vp_loc, draw_time, logentry, rettree, cmap

from .vi_func import sunpath, windnum, wind_rose, create_coll, retobjs, progressfile, progressbar
from .vi_func import chunks, clearlayers, clearscene, clearfiles, objmode
#from .vi_display import wr_legend, wr_scatter, wr_table, wr_disp
#from .envi_func import processf, retenvires, envizres, envilres, recalculate_text
#from .vi_chart import chart_disp

from .vi_misc import wr_legend, results_bar, wr_table, wr_scatter, svf_legend
from .vi_misc import li_display, linumdisplay, ss_legend, ss_scatter, livi_legend#, spnumdisplay, en_air, wr_legend, wr_disp, wr_scatter, wr_table, ss_disp, ss_legend, svf_disp, svf_legend, basic_legend, basic_table, basic_disp, ss_scatter, en_disp, en_pdisp, en_scatter, en_table, en_barchart, comp_table, comp_disp, leed_scatter, cbdm_disp, cbdm_scatter, envals, bsdf, bsdf_disp#, en_barchart, li3D_legend

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
        svp = scene.vi_params
    
        if bpy.data.objects.get('SPathMesh'):
            spob = bpy.data.objects['SPathMesh'] 
            ob_mat = spob.matrix_world
            mid_x, mid_y, width, height = viewdesc(context)
            vl = ret_vp_loc(context)
            blf_props(scene, width, height)
            
            if svp.sp_hd:
                coords, scene, sd = {}, context.scene, 100
                dists, hs, pos = [], [], []
                
                for doy in (172, 355):
                    for hour in range(24):
                        ([solalt, solazi]) = solarPosition(doy, hour, svp.latitude, svp.longitude)[2:]
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
                    draw_index(pos, hs, dists, svp.vi_display_rp_fs, svp.vi_display_rp_fc, svp.vi_display_rp_fsh)
                    
            if [ob.get('VIType') == 'Sun' for ob in bpy.data.objects] and svp['spparams']['suns'] == '0':
                sobs = [ob for ob in bpy.data.objects if ob.get('VIType') == 'Sun']
                
                if sobs and svp.sp_td:
                    sunloc = ob_mat@sobs[0].location
                    solpos = view3d_utils.location_3d_to_region_2d(context.region, context.region_data, sunloc)
                    
                    try:
                        if 0 < solpos[0] < width and 0 < solpos[1] < height and not scene.ray_cast(context.view_layer, sobs[0].location + 0.05 * (vl - sunloc), vl - sunloc)[0]:
                            soltime = datetime.datetime.fromordinal(svp.sp_sd)
                            soltime += datetime.timedelta(hours = svp.sp_sh)
                            sre = sobs[0].rotation_euler
                            blf_props(scene, width, height)
                            sol_text = soltime.strftime('     %d %b %X') + ' alt: {:.1f} azi: {:.1f}'.format(90 - sre[0]*180/pi, (180, -180)[sre[2] < -pi] - sre[2]*180/pi)
                            draw_time(solpos, sol_text, svp.vi_display_rp_fs, 
                                      svp.vi_display_rp_fc, svp.vi_display_rp_fsh)
                            
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
            [bpy.data.objects.remove(o, do_unlink=True, do_id_user=True, do_ui_user=True) for o in bpy.data.objects if o.vi_params.get('VIType') and o.vi_params['VIType'] in ('SunMesh', 'SkyMesh')]
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
        ovp = wro.vi_params
        ovp['maxres'], ovp['minres'], ovp['avres'], ovp['nbins'], ovp['VIType'] = max(aws), min(aws), sum(aws)/len(aws), len(sbinvals), 'Wind_Plane'
        simnode['maxfreq'] = 100*numpy.max(adfreq)/len(cwd)
        simnode['maxfreq'] = 100*numpy.max(adfreq)/len(vawd) if simnode.max_freq == '0' else simnode.max_freq_val

        windnum(simnode['maxfreq'], (0,0,0), scale, wind_compass((0,0,0), scale, wro, wro.data.materials['wr-000000']))
        plt.close()
        ovp['table'] = array([["", 'Minimum', 'Average', 'Maximum'], 
                             ['Speed (m/s)', ovp['minres'], '{:.1f}'.format(ovp['avres']), ovp['maxres']], 
                             ['Direction (\u00B0)', min(awd), '{:.1f}'.format(sum(awd)/len(awd)), max(awd)]])
        ovp['ws'] = aws.reshape(len(doys), 24).T.tolist()
        ovp['wd'] = awd.reshape(len(doys), 24).T.tolist()
        ovp['days'] = array(doys, dtype = float)
        ovp['hours'] = arange(1, 25, dtype = float)        
        ovp['maxfreq'] = 100*numpy.max(dfreq)/len(awd)
        simnode['nbins'] = len(sbinvals)  
        simnode['ws'] = array(cws).reshape(365, 24).T.tolist()
        simnode['wd'] = array(cwd).reshape(365, 24).T.tolist()        
        simnode['days'] = arange(1, 366, dtype = float)
        simnode['hours'] = arange(1, 25, dtype = float)        
        return {'FINISHED'}
    
class VIEW3D_OT_WRDisplay(bpy.types.Operator):
    bl_idname = "view3d.wrdisplay"
    bl_label = "Wind rose number display"
    bl_description = "Project the windrose numbers on to the viewport"
    bl_register = True
    bl_undo = False
    
    def invoke(self, context, event):   
        area = context.area
        svp = context.scene.vi_params
        svp.vi_display = 1
        svp['viparams']['vidisp'] = 'wr'
        self.results_bar = results_bar(('legend.png', 'table_new.png', 'scatter_new.png'), 300, area)
        self.legend = wr_legend(context, 'Speed (m/s)', [305, area.height - 80], area.width, area.height, 125, 300)
        self.table = wr_table(context, [355, area.height - 80], area.width, area.height, 400, 60) 
        self.dhscatter = wr_scatter(context, [405, area.height - 80], area.width, area.height, 600, 200)
        self.height = area.height
        self.draw_handle_wrnum = bpy.types.SpaceView3D.draw_handler_add(self.draw_wrnum, (context, ), 'WINDOW', 'POST_PIXEL')
        self.cao = context.active_object
        area.tag_redraw()
        context.window_manager.modal_handler_add(self)
        return {'RUNNING_MODAL'}
    
    def modal(self, context, event):
        scene = context.scene
        svp = scene.vi_params
        ah, aw = context.area.height, context.area.width
        redraw = 0
           
        if svp.vi_display == 0 or svp['viparams']['vidisp'] != 'wr' or event.type == 'ESC':
            svp.vi_display = 0
            bpy.types.SpaceView3D.draw_handler_remove(self.draw_handle_wrnum, 'WINDOW')
            context.area.tag_redraw()
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
            
            elif self.table.ispos[0] < mx < self.table.iepos[0] and self.table.ah - 80 < my < self.table.ah - 40:
                self.table.hl = (0.8, 0.8, 0.8, 0.8) 
                redraw = 1
                if event.type == 'LEFTMOUSE':
                    if event.value == 'RELEASE':
                        self.table.expand = 0 if self.table.expand else 1
                        
            elif self.dhscatter.ispos[0] < mx < self.dhscatter.iepos[0] and self.dhscatter.ah - 80 < my < self.dhscatter.ah - 40:
                self.dhscatter.hl = (0.8, 0.8, 0.8, 0.8) 
                redraw = 1
                if event.type == 'LEFTMOUSE':
                    if event.value == 'RELEASE':
                        self.dhscatter.expand = 0 if self.dhscatter.expand else 1
                        
            elif self.dhscatter.expand and self.dhscatter.lspos[0] + 0.1 * self.dhscatter.xdiff < mx < self.dhscatter.lepos[0] - 0.1 * self.dhscatter.xdiff and self.dhscatter.lspos[1] + 0.1 * self.dhscatter.ydiff  < my < self.dhscatter.lepos[1] - 0.1 * self.dhscatter.ydiff:
                self.dhscatter.hl = (0.8, 0.8, 0.8, 0.8) 
                redraw = 1
                if event.type == 'LEFTMOUSE' and event.value == 'RELEASE':
                    self.dhscatter.show_plot(context)
                        
            elif self.legend.expand and abs(self.legend.lspos[0] - mx) < 10 and abs(self.legend.lepos[1] - my) < 10:
                self.legend.hl = (0.8, 0.8, 0.8, 0.8) 
                redraw = 1   
                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.legend.move = 1
                        self.legend.draw(ah, aw)
                        context.area.tag_redraw()
#                        return {'RUNNING_MODAL'}
                    elif self.legend.move and event.value == 'RELEASE':
                        self.legend.move = 0                        
                    return {'RUNNING_MODAL'}
            
            elif self.table.expand and abs(self.table.lspos[0] - mx) < 10 and abs(self.table.lepos[1] - my) < 10:
                self.table.hl = (0.8, 0.8, 0.8, 0.8) 
                redraw = 1   
                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.table.move = 1
                        self.table.draw(ah, aw)
                        context.area.tag_redraw()
#                        return {'RUNNING_MODAL'}
                    elif self.table.move and event.value == 'RELEASE':
                        self.table.move = 0                        
                    return {'RUNNING_MODAL'}
            
            elif self.dhscatter.expand and abs(self.dhscatter.lspos[0] - mx) < 10 and abs(self.dhscatter.lepos[1] - my) < 10:
                self.dhscatter.hl = (0.8, 0.8, 0.8, 0.8) 
                redraw = 1   
                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.dhscatter.move = 1
                        self.dhscatter.draw(ah, aw)
                        context.area.tag_redraw()
#                        return {'RUNNING_MODAL'}
                    elif self.dhscatter.move and event.value == 'RELEASE':
                        self.dhscatter.move = 0                        
                    return {'RUNNING_MODAL'}
                    
            elif self.legend.expand and abs(self.legend.lepos[0] - mx) < 10 and abs(self.legend.lspos[1] - my) < 10:
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
            
            elif self.table.expand and abs(self.table.lepos[0] - mx) < 10 and abs(self.table.lspos[1] - my) < 10:
                self.table.hl = (0.8, 0.8, 0.8, 0.8) 
                context.area.tag_redraw()
                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.table.resize = 1
                        self.table.draw(ah, aw)
                        context.area.tag_redraw()
                    elif self.table.resize and event.value == 'RELEASE':
                        self.table.resize = 0
                    return {'RUNNING_MODAL'}
                
            elif self.dhscatter.expand and abs(self.dhscatter.lepos[0] - mx) < 10 and abs(self.dhscatter.lspos[1] - my) < 10:
                self.dhscatter.hl = (0.8, 0.8, 0.8, 0.8) 
                context.area.tag_redraw()
                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.dhscatter.resize = 1
                        self.dhscatter.draw(ah, aw)
                        context.area.tag_redraw()
                    elif self.dhscatter.resize and event.value == 'RELEASE':
                        self.dhscatter.resize = 0
                    return {'RUNNING_MODAL'}
                
            elif self.legend.hl == (0.8, 0.8, 0.8, 0.8):                 
                self.legend.hl = (1, 1, 1, 1)
                redraw = 1
            
            elif self.table.hl == (0.8, 0.8, 0.8, 0.8):                 
                self.table.hl = (1, 1, 1, 1)
                redraw = 1
                
            elif self.dhscatter.hl == (0.8, 0.8, 0.8, 0.8):                 
                self.dhscatter.hl = (1, 1, 1, 1)
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
                elif self.table.move:
                    self.table.lspos[0], self.table.lepos[1] = mx, my
                    self.table.draw(ah, aw)
                    context.area.tag_redraw() 
                elif self.table.resize:
                    self.table.lepos[0], self.table.lspos[1] = mx, my
                    self.table.draw(ah, aw)
                    context.area.tag_redraw()
                elif self.dhscatter.resize:
                    self.dhscatter.lepos[0], self.dhscatter.lspos[1] = mx, my
                    self.dhscatter.draw(ah, aw)
                    context.area.tag_redraw() 
                elif self.dhscatter.move:
                    self.dhscatter.lspos[0], self.dhscatter.lepos[1] = mx, my
                    self.dhscatter.draw(ah, aw)
                    context.area.tag_redraw() 
        
            if self.cao != context.active_object:
                self.legend.update(context)
                self.legend.draw(ah, aw)
                self.table.update(context)
                self.table.draw(ah, aw)
                self.dhscatter.update(context)
                self.dhscatter.draw(ah, aw)
                context.area.tag_redraw() 
                self.cao = context.active_object
                
            if svp.vi_wr_refresh:
                self.dhscatter.update(context)
                self.dhscatter.draw(ah, aw)
                svp.vi_wr_refresh = 0
                context.area.tag_redraw()
            
            if redraw:
                context.area.tag_redraw()        
        return {'PASS_THROUGH'}
    
    def draw_wrnum(self, context):
        area = context.area
        ah = area.height                
        self.results_bar.draw(ah, 3)
        self.legend.draw(ah, area.width)
        self.table.draw(ah, area.width)
        self.dhscatter.draw(ah, area.width)
                
class NODE_OT_SVF(bpy.types.Operator):
    bl_idname = "node.svf"
    bl_label = "Sky View Factor"
    bl_description = "Undertake a sky view factor study"
    bl_register = True
    bl_undo = False
#    nodeid = bpy.props.StringProperty()

    def invoke(self, context, event):
        scene = context.scene 
        svp = scene.vi_params
        svp.vi_display = 0
        
        if viparams(self, scene):            
            return {'CANCELLED'}

        shadobs = retobjs('livig')
        if not shadobs:
            self.report({'ERROR'},"No shading objects have a material attached.")
            return {'CANCELLED'}
            
        svp['liparams']['shadc'] = [ob.name for ob in retobjs('ssc')]
        if not svp['liparams']['shadc']:
            self.report({'ERROR'},"No objects have a light sensor material attached.")
            return {'CANCELLED'}
        simnode = context.node
        svp['viparams']['restree'] = simnode.id_data.name
        clearscene(scene, self)
        
        svp['viparams']['visimcontext'] = 'SVF'

        if not svp.get('liparams'):
           svp['liparams'] = {}
           
        svp['liparams']['cp'], svp['liparams']['unit'], svp['liparams']['type'] = simnode.cpoint, '% Sunlit', 'VI Shadow'
        simnode.preexport()
        (svp['liparams']['fs'], svp['liparams']['fe']) = (scene.frame_current, scene.frame_current) if simnode.animmenu == 'Static' else (simnode.startframe, simnode.endframe)      
        svp['viparams']['resnode'], simnode['Animation'] = simnode.name, simnode.animmenu
        (scmaxres, scminres, scavres) = [[x] * (svp['liparams']['fe'] - svp['liparams']['fs'] + 1) for x in (0, 1, 0)]        
        frange = range(svp['liparams']['fs'], svp['liparams']['fe'] + 1)
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
        calcsteps = len(frange) * sum(len([f for f in o.data.polygons if o.data.materials[f.material_index].vi_params.mattype == '1']) for o in [scene.objects[on] for on in svp['liparams']['shadc']])
        curres, reslists = 0, []
        pfile = progressfile(svp['viparams']['newdir'], datetime.datetime.now(), calcsteps)
        kivyrun = progressbar(os.path.join(svp['viparams']['newdir'], 'viprogress'), 'Sky View')
        
        for oi, o in enumerate([scene.objects[on] for on in svp['liparams']['shadc']]):
            ovp = o.vi_params
            for k in ovp.keys():
                del ovp[k]
                
            if any([s < 0 for s in o.scale]):
                logentry('Negative scaling on calculation object {}. Results may not be as expected'.format(o.name))
                self.report({'WARNING'}, 'Negative scaling on calculation object {}. Results may not be as expected'.format(o.name))

            ovp['omin'], ovp['omax'], ovp['oave'] = {}, {}, {}
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
                gpoints = [f for f in geom if o.data.materials[f.material_index].vi_params.mattype == '1']
            elif simnode.cpoint == '1':
                gpoints = [v for v in geom if any([o.data.materials[f.material_index].vi_params.mattype == '1' for f in v.link_faces])]

            for g, gp in enumerate(gpoints):
                gp[cindex] = g + 1

            for frame in frange: 
                g, oshadres = 0, array([])                
                scene.frame_set(frame)
                shadtree = rettree(scene, shadobs, ('', '2')[simnode.signore])
                shadres = geom.layers.float['res{}'.format(frame)]
                                    
                if gpoints:
                    posis = [gp.calc_center_bounds() + gp.normal.normalized() * simnode.offset for gp in gpoints] if simnode.cpoint == '0' else [gp.co + gp.normal.normalized() * simnode.offset for gp in gpoints]
                 
                    for chunk in chunks(gpoints, int(scene['viparams']['nproc']) * 200):
                        for gp in chunk:
                            pointres = array([(0, 1)[shadtree.ray_cast(posis[g], direc)[3] == None] for direc in valdirecs], dtype = int8)
                            gp[shadres] = ((numpy.sum(pointres)/lvaldirecs)).astype(float16)
                            g += 1

                        curres += len(chunk)
                        if pfile.check(curres) == 'CANCELLED':
                            return {'CANCELLED'}
              
                    shadres = [gp[shadres] for gp in gpoints]
                    ovp['omin']['res{}'.format(frame)], ovp['omax']['res{}'.format(frame)], ovp['oave']['res{}'.format(frame)] = min(shadres), max(shadres), sum(shadres)/len(shadres)
                    reslists.append([str(frame), 'Zone', o.name, 'X', ' '.join(['{:.3f}'.format(p[0]) for p in posis])])
                    reslists.append([str(frame), 'Zone', o.name, 'Y', ' '.join(['{:.3f}'.format(p[1]) for p in posis])])
                    reslists.append([str(frame), 'Zone', o.name, 'Z', ' '.join(['{:.3f}'.format(p[2]) for p in posis])])
                    reslists.append([str(frame), 'Zone', o.name, 'SVF', ' '.join(['{:.3f}'.format(sr) for sr in oshadres])])
                    avres.append(ovp['oave']['res{}'.format(frame)])
                    minres.append(ovp['omin']['res{}'.format(frame)])
                    maxres.append(ovp['omax']['res{}'.format(frame)])

            reslists.append(['All', 'Frames', '', 'Frames', ' '.join(['{}'.format(f) for f in frange])])
            reslists.append(['All', 'Zone', o.name, 'Minimum', ' '.join(['{:.3f}'.format(mr) for mr in minres])])
            reslists.append(['All', 'Zone', o.name, 'Average', ' '.join(['{:.3f}'.format(mr) for mr in avres])])
            reslists.append(['All', 'Zone', o.name, 'Maximum', ' '.join(['{:.3f}'.format(mr) for mr in maxres])])
            bm.transform(o.matrix_world.inverted())
            bm.to_mesh(o.data)
            bm.free()

        svp.vi_leg_max, svp.vi_leg_min = 1, 0

        if kivyrun.poll() is None:
            kivyrun.kill()
        
        scene.frame_start, scene.frame_end = svp['liparams']['fs'], svp['liparams']['fe']
        svp['viparams']['vidisp'] = 'svf'
        simnode['reslists'] = reslists
        simnode['frames'] = [f for f in frange]
        simnode.postexport(scene)
        return {'FINISHED'} 
  
class VIEW3D_OT_SVFDisplay(bpy.types.Operator):
    '''Display results legend and stats in the 3D View'''
    bl_idname = "view3d.svfdisplay"
    bl_label = "Shadow study metric display"
    bl_description = "Display shadow study metrics"
    bl_register = True
    bl_undo = False
    
    def invoke(self, context, event):
        area = context.area
        svp = context.scene.vi_params
        try:
            bpy.types.SpaceView3D.draw_handler_remove(self.draw_handle_svfnum, 'WINDOW')
#            bpy.types.SpaceView3D.draw_handler_remove(self._handle_ss_disp, 'WINDOW')
        except:
            pass
        svp.vi_display = 1
        svp['viparams']['vidisp'] = 'svf'
        self.simnode = bpy.data.node_groups[svp['viparams']['restree']].nodes[svp['viparams']['resnode']]
        li_display(self, self.simnode)
        self.results_bar = results_bar(('legend.png',), 300, area)
        self.legend = svf_legend(context, 'Sky View (%)', [305, area.height - 80], area.width, area.height, 125, 400)
        self.legend_num = linumdisplay(self, context)
        self.height = area.height
        self.draw_handle_svfnum = bpy.types.SpaceView3D.draw_handler_add(self.draw_svfnum, (context, ), 'WINDOW', 'POST_PIXEL')
        self.cao = context.active_object
        area.tag_redraw()
        context.window_manager.modal_handler_add(self)
        return {'RUNNING_MODAL'}
    
    def draw_svfnum(self, context):
        area = context.area
        ah = area.height                
        self.results_bar.draw(ah)
        self.legend.draw(context)
        self.legend_num.draw(context)
        
    def modal(self, context, event):    
        scene = context.scene
        svp = scene.vi_params
        redraw = 0
           
        if svp.vi_display == 0 or svp['viparams']['vidisp'] != 'svf' or event.type == 'ESC':
            svp.vi_display = 0
            bpy.types.SpaceView3D.draw_handler_remove(self.draw_handle_ssvfnum, 'WINDOW')
            context.area.tag_redraw()
            return {'CANCELLED'}

        if event.type != 'INBETWEEN_MOUSEMOVE' and context.region and context.area.type == 'VIEW_3D' and context.region.type == 'WINDOW':            
            mx, my = event.mouse_region_x, event.mouse_region_y 

#           for ui_element in (self.ui_elements):
                
            # Legend routine 
            if self.legend.ispos[0] < mx < self.legend.iepos[0] and self.legend.ah - 80 < my < self.legend.ah - 40:
                self.legend.hl = (0.8, 0.8, 0.8, 0.8) 
                redraw = 1
                if event.type == 'LEFTMOUSE':
                    if event.value == 'RELEASE':
                        self.legend.expand = 0 if self.legend.expand else 1
                        
            elif self.legend.expand and abs(self.legend.lspos[0] - mx) < 10 and abs(self.legend.lepos[1] - my) < 10:
                self.legend.hl = (0.8, 0.8, 0.8, 0.8) 
                redraw = 1   
                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.legend.move = 1
                        self.legend.draw(context)
                        context.area.tag_redraw()
                    elif self.legend.move and event.value == 'RELEASE':
                        self.legend.move = 0                        
                    return {'RUNNING_MODAL'}
              
            elif self.legend.expand and abs(self.legend.lepos[0] - mx) < 10 and abs(self.legend.lspos[1] - my) < 10:
                self.legend.hl = (0.8, 0.8, 0.8, 0.8) 
                context.area.tag_redraw()
                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.legend.resize = 1
                        self.legend.draw(context)
                        context.area.tag_redraw()
                    elif self.legend.resize and event.value == 'RELEASE':
                        self.legend.resize = 0
                    return {'RUNNING_MODAL'}  
                
            elif self.legend.hl == (0.8, 0.8, 0.8, 0.8):                 
                self.legend.hl = (1, 1, 1, 1)
                redraw = 1
                
            if event.type == 'MOUSEMOVE':                
                if self.legend.move:
                    self.legend.lspos[0], self.legend.lepos[1] = mx, my
                    self.legend.draw(context)
                    context.area.tag_redraw() 
                elif self.legend.resize:
                    self.legend.lepos[0], self.legend.lspos[1] = mx, my
                    self.legend.draw(context)
                    context.area.tag_redraw()
    
        if redraw:
            context.area.tag_redraw()
        return {'PASS_THROUGH'}
    
class NODE_OT_Shadow(bpy.types.Operator):
    bl_idname = "node.shad"
    bl_label = "Shadow Map"
    bl_description = "Undertake a shadow mapping"
    bl_register = True
    bl_undo = False

    def invoke(self, context, event):
        scene = context.scene 
        svp = scene.vi_params
        svp.vi_display = 0
        
       
        
        if viparams(self, scene):            
            return {'CANCELLED'}

        shadobs = retobjs('livig')
        if not shadobs:
            self.report({'ERROR'},"No shading objects have a material attached.")
            return {'CANCELLED'}
            
        svp['liparams']['shadc'] = [ob.name for ob in retobjs('ssc')]
        if not svp['liparams']['shadc']:
            self.report({'ERROR'},"No objects have a VI Shadow material attached.")
            return {'CANCELLED'}
        simnode = context.node
        svp['viparams']['restree'] = simnode.id_data.name
        
        clearscene(scene, self)
        
        svp['viparams']['visimcontext'] = 'Shadow'
        if not svp.get('liparams'):
           svp['liparams'] = {}
        svp['liparams']['cp'], svp['liparams']['unit'], svp['liparams']['type'] = simnode.cpoint, '% Sunlit', 'VI Shadow'
        simnode.preexport()
        (svp['liparams']['fs'], svp['liparams']['fe']) = (scene.frame_current, scene.frame_current) if simnode.animmenu == 'Static' else (simnode.startframe, simnode.endframe)
        cmap(svp)

        if simnode.starthour > simnode.endhour:
            self.report({'ERROR'},"End hour is before start hour.")
            return{'CANCELLED'}
        
        svp['viparams']['resnode'], simnode['Animation'] = simnode.name, simnode.animmenu
        (scmaxres, scminres, scavres) = [[x] * (svp['liparams']['fe'] - svp['liparams']['fs'] + 1) for x in (0, 100, 0)]
        
        frange = range(svp['liparams']['fs'], svp['liparams']['fe'] + 1)
        time = datetime.datetime(2014, simnode.sdate.month, simnode.sdate.day, simnode.starthour - 1)
        y =  2014 if simnode.edoy >= simnode.sdoy else 2014 + 1
        endtime = datetime.datetime(y, simnode.edate.month, simnode.edate.day, simnode.endhour - 1)
        interval = datetime.timedelta(hours = 1/simnode.interval)
        
        times = [time + interval*t for t in range(int((endtime - time)/interval) + simnode.interval) if simnode.starthour - 1 <= (time + interval*t).hour <= simnode.endhour  - 1]
        sps = [solarPosition(t.timetuple().tm_yday, t.hour+t.minute/60, svp.latitude, svp.longitude)[2:] for t in times]
        valmask = array([sp[0] > 0 for sp in sps], dtype = int8)
        direcs = [mathutils.Vector((-sin(sp[1]), -cos(sp[1]), tan(sp[0]))) for sp in sps]  
        valdirecs = [mathutils.Vector((-sin(sp[1]), -cos(sp[1]), tan(sp[0]))) for sp in sps if sp[0] > 0]  
        lvaldirecs = len(valdirecs)
        ilvaldirecs = 1/lvaldirecs
        calcsteps = len(frange) * sum(len([f for f in o.data.polygons if o.data.materials[f.material_index].vi_params.mattype == '1']) for o in [scene.objects[on] for on in svp['liparams']['shadc']])
        curres, reslists = 0, []
        pfile = progressfile(svp['viparams']['newdir'], datetime.datetime.now(), calcsteps)
        kivyrun = progressbar(os.path.join(scene.vi_params['viparams']['newdir'], 'viprogress'), 'Shadow Map')
        logentry(f'Conducting shadow map calculation with {simnode.interval} samples per hour for {int(len(direcs)/simnode.interval)} total hours and {lvaldirecs} available sun hours')
        
        for oi, o in enumerate([scene.objects[on] for on in svp['liparams']['shadc']]):
            ovp = o.vi_params
            for k in ovp.keys():
                del ovp[k]
                
            if any([s < 0 for s in o.scale]):
                logentry('Negative scaling on calculation object {}. Results may not be as expected'.format(o.name))
                self.report({'WARNING'}, 'Negative scaling on calculation object {}. Results may not be as expected'.format(o.name))

            ovp['omin'], ovp['omax'], ovp['oave'] = {}, {}, {}
            
            if simnode.sdoy <= simnode.edoy:
                ovp['days'] = arange(simnode.sdoy, simnode.edoy + 1, dtype = float)
            else:
                ovp['days'] = arange(simnode.sdoy, simnode.edoy + 1, dtype = float)
                
            ovp['hours'] = arange(simnode.starthour, simnode.endhour + 1, 1/simnode.interval, dtype = float)
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
                gpoints = [f for f in geom if o.data.materials[f.material_index].vi_params.mattype == '1']
            elif simnode.cpoint == '1':
                gpoints = [v for v in geom if any([o.data.materials[f.material_index].vi_params.mattype == '1' for f in v.link_faces])]

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

                            pointres = array([(0, 1)[shadtree.ray_cast(posis[g], direc)[3] == None] for direc in valdirecs], dtype = int8)
                            numpy.place(allpoints[g], valmask == 1, pointres)
                            gp[shadres] = (100 * (numpy.sum(pointres) * ilvaldirecs)).astype(float16)
                            g += 1

                        curres += len(chunk)
                        if pfile.check(curres) == 'CANCELLED':
                            return {'CANCELLED'}
                    
                    ap = numpy.average(allpoints, axis=0)                
                    shadres = [gp[shadres] for gp in gpoints]
                    ovp['dhres{}'.format(frame)] = array(100 * ap).reshape(len(ovp['days']), len(ovp['hours'])).T.tolist()
                    ovp['omin']['res{}'.format(frame)], ovp['omax']['res{}'.format(frame)], ovp['oave']['res{}'.format(frame)] = min(shadres), max(shadres), sum(shadres)/len(shadres)
                    reslists.append([str(frame), 'Zone', o.name, 'X', ' '.join(['{:.3f}'.format(p[0]) for p in posis])])
                    reslists.append([str(frame), 'Zone', o.name, 'Y', ' '.join(['{:.3f}'.format(p[1]) for p in posis])])
                    reslists.append([str(frame), 'Zone', o.name, 'Z', ' '.join(['{:.3f}'.format(p[2]) for p in posis])])
                    reslists.append([str(frame), 'Zone', o.name, 'Sunlit %', ' '.join(['{:.3f}'.format(sr) for sr in oshadres])])
                    avres.append(ovp['oave']['res{}'.format(frame)])
                    minres.append(ovp['omin']['res{}'.format(frame)])
                    maxres.append(ovp['omax']['res{}'.format(frame)])
            
            reslists.append(['All', 'Frames', '', 'Frames', ' '.join(['{}'.format(f) for f in frange])])
            reslists.append(['All', 'Zone', o.name, 'Minimum', ' '.join(['{:.3f}'.format(mr) for mr in minres])])
            reslists.append(['All', 'Zone', o.name, 'Average', ' '.join(['{:.3f}'.format(mr) for mr in avres])])
            reslists.append(['All', 'Zone', o.name, 'Maximum', ' '.join(['{:.3f}'.format(mr) for mr in maxres])])
            
            bm.transform(o.matrix_world.inverted())
            bm.to_mesh(o.data)
            bm.free()

        svp.vi_leg_max, svp.vi_leg_min = 100, 0

        if kivyrun.poll() is None:
            kivyrun.kill()

        scene.frame_start, scene.frame_end = svp['liparams']['fs'], svp['liparams']['fe']
        simnode['reslists'] = reslists
        simnode['frames'] = [f for f in frange]
        simnode['year'] = 2015
        simnode.postexport(scene)
        svp['viparams']['vidisp'] = 'ss'
        return {'FINISHED'}

class VIEW3D_OT_SSDisplay(bpy.types.Operator):
    '''Display results legend and stats in the 3D View'''
    bl_idname = "view3d.ssdisplay"
    bl_label = "Shadow study metric display"
    bl_description = "Display shadow study metrics"
    bl_register = True
    bl_undo = False

    def modal(self, context, event): 
        scene = context.scene
        svp = scene.vi_params
        redraw = 0 

        if svp.vi_display == 0 or svp['viparams']['vidisp'] != 'ss' or not [o for o in bpy.data.objects if o.name in svp['liparams']['shadc']]:
            bpy.types.SpaceView3D.draw_handler_remove(self.draw_handle_ssnum, 'WINDOW')
            context.area.tag_redraw()
            return {'CANCELLED'}        
        
        if event.type != 'INBETWEEN_MOUSEMOVE' and context.region and context.area.type == 'VIEW_3D' and context.region.type == 'WINDOW':                    
            mx, my = event.mouse_region_x, event.mouse_region_y 
            
#            if any((svp.vi_leg_levels != self.legend.levels, svp.vi_leg_col != self.legend.col, svp.vi_leg_scale != self.legend.scale, (self.legend.minres, self.legend.maxres) != leg_min_max(self.scene))):               
#                self.legend.update(context)                
#                redraw = 1
                 
            # Legend routine 
            
            if self.legend.ispos[0] < mx < self.legend.iepos[0] and self.legend.ah - 80 < my < self.legend.ah - 40:
                self.legend.hl = (0.8, 0.8, 0.8, 0.8) 
                redraw = 1
                if event.type == 'LEFTMOUSE':
                    if event.value == 'RELEASE':
                        self.legend.expand = 0 if self.legend.expand else 1
                        
            elif self.legend.expand and abs(self.legend.lspos[0] - mx) < 10 and abs(self.legend.lepos[1] - my) < 10:
                self.legend.hl = (0.8, 0.8, 0.8, 0.8) 
                redraw = 1   
                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.legend.move = 1
                        self.legend.draw(context)
                        context.area.tag_redraw()
                    elif self.legend.move and event.value == 'RELEASE':
                        self.legend.move = 0                        
                    return {'RUNNING_MODAL'}
              
            elif self.legend.expand and abs(self.legend.lepos[0] - mx) < 10 and abs(self.legend.lspos[1] - my) < 10:
                self.legend.hl = (0.8, 0.8, 0.8, 0.8) 
                context.area.tag_redraw()
                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.legend.resize = 1
                        self.legend.draw(context)
                        context.area.tag_redraw()
                    elif self.legend.resize and event.value == 'RELEASE':
                        self.legend.resize = 0
                    return {'RUNNING_MODAL'}  
                
            elif self.legend.hl == (0.8, 0.8, 0.8, 0.8):                 
                self.legend.hl = (1, 1, 1, 1)
                redraw = 1
                
            if event.type == 'MOUSEMOVE':                
                if self.legend.move:
                    self.legend.lspos[0], self.legend.lepos[1] = mx, my
                    self.legend.draw(context)
                    context.area.tag_redraw() 
                elif self.legend.resize:
                    self.legend.lepos[0], self.legend.lspos[1] = mx, my
                    self.legend.draw(context)
                    context.area.tag_redraw()
            
            if redraw:
                context.area.tag_redraw()

        return {'PASS_THROUGH'}

    def invoke(self, context, event):
        self.scene = context.scene
        area = context.area
        svp = context.scene.vi_params
        svp.vi_display = 1
        
        
        try:
            bpy.types.SpaceView3D.draw_handler_remove(self.draw_handle_ssnum, 'WINDOW')
#            bpy.types.SpaceView3D.draw_handler_remove(self._handle_ss_disp, 'WINDOW')
        except:
            pass
        
        clearscene(self.scene, self)
        svp['viparams']['vidisp'] = 'ss' 
        self.simnode = bpy.data.node_groups[svp['viparams']['restree']].nodes[svp['viparams']['resnode']]
        li_display(self, self.simnode)
        self.results_bar = results_bar(('legend.png', 'scatter.png'), 300, area)
        self.legend = ss_legend(context, 'Sunlit (%)', [305, area.height - 80], area.width, area.height, 125, 400)
        self.num_display = linumdisplay(self, context)
        self.dhscatter = ss_scatter(context, [355, area.height - 80], area.width, area.height, 600, 200)
        svp.vi_disp_wire = 1
        self.draw_handle_ssnum = bpy.types.SpaceView3D.draw_handler_add(self.draw_ssnum, (context, ), 'WINDOW', 'POST_PIXEL')        
#        self.legend.update(context)
#        self.dhscatter.update(context)
#        self._handle_ss_disp = bpy.types.SpaceView3D.draw_handler_add(ss_disp, (self, context, self.simnode), 'WINDOW', 'POST_PIXEL')
        context.window_manager.modal_handler_add(self)
        return {'RUNNING_MODAL'}

    def draw_ssnum(self, context):
        area = context.area
        ah = area.height                
        self.results_bar.draw(ah)
        self.legend.draw(context)
        self.dhscatter.draw(context, area.width)
        self.num_display.draw(context)
        
class NODE_OT_Li_Geo(bpy.types.Operator):
    bl_idname = "node.ligexport"
    bl_label = "LiVi geometry export"
#    nodeid = bpy.props.StringProperty()

    def invoke(self, context, event):
        scene = context.scene
        svp = scene.vi_params
        if viparams(self, scene):
            return {'CANCELLED'}
        scene['viparams']['vidisp'] = ''
        scene['viparams']['viexpcontext'] = 'LiVi Geometry'
        objmode()
        clearfiles(svp['liparams']['objfilebase'])
        clearfiles(svp['liparams']['lightfilebase'])
        node = context.node
        node.preexport(scene)
        radgexport(self, node)
        node.postexport(scene)
        return {'FINISHED'}
    
class NODE_OT_Li_Con(bpy.types.Operator, io_utils.ExportHelper):
    bl_idname = "node.liexport"
    bl_label = "LiVi context export"
    bl_description = "Export the scene to the Radiance file format"
    bl_register = True
    bl_undo = False

    def invoke(self, context, event):
        scene = context.scene
        self.svp = scene.vi_params
        if viparams(self, scene):
            return {'CANCELLED'}
        node = context.node
        self.svp['viparams']['vidisp'] = ''
        self.svp['viparams']['viexpcontext'] = 'LiVi {}'.format(node.contextmenu)
        self.svp['viparams']['connode'] = '{}@{}'.format(node.name, node.id_data.name)
        objmode()
        node.preexport()
        if not node.export(scene, self):
            node.postexport()
        return {'FINISHED'}
    
#    def modal(self, context, event):
#        if self.svp['visimcontext'] != 'LiVi' or not [o for o in context.scene.objects if o.get('VIType') and o['VIType'] == 'Sun']:
#            return {'FINISHED'}
#        else:
#            # This could be used to se light direction
#            fc = scene.frame_start if scene.frame_current > scene.frame_end else scene.frame_current
#            scene.display.light_direction = (-sin(solposs[fc][3]) * cos(solposs[fc][2]), sin(solposs[fc][2]),  cos(solposs[fc][3]) * cos(solposs[fc][2])) 


class NODE_OT_Li_Pre(bpy.types.Operator, io_utils.ExportHelper):
    bl_idname = "node.radpreview"
    bl_label = "LiVi preview"
    bl_description = "Prevew the scene with Radiance"
    bl_register = True
    bl_undo = False
#    nodeid = bpy.props.StringProperty()
    
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
        svp = scene.vi_params
        if viparams(self, scene):
            return {'CANCELLED'}
        
        objmode()
        self.simnode = context.node
        frame = scene.frame_current
        self.simnode.presim()
        svp['liparams']['fs'] = min([c['fs'] for c in (self.simnode['goptions'], self.simnode['coptions'])])
        svp['liparams']['fe'] = max([c['fe'] for c in (self.simnode['goptions'], self.simnode['coptions'])])

        if frame not in range(svp['liparams']['fs'], svp['liparams']['fe'] + 1):
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
                pmcmd = 'mkpmap -t 20 -e {1}.pmapmon -fo+ -bv+ -apD 0.001 {0} -apg {1}-{2}.gpm {3} {4} {5} {1}-{2}.oct'.format(pportentry, svp['viparams']['filebase'], frame, self.simnode.pmapgno, cpentry, amentry)
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
                                        
                rvucmd = "rvu -w {11} -ap {8} 50 {9} -n {0} -vv {1:.3f} -vh {2} -vd {3[0]:.3f} {3[1]:.3f} {3[2]:.3f} -vp {4[0]:.3f} {4[1]:.3f} {4[2]:.3f} -vu {10[0]:.3f} {10[1]:.3f} {10[2]:.3f} {5} {6}-{7}.oct".format(svp['viparams']['wnproc'], 
                                 vv, cang, vd, cam.location, self.simnode['radparams'], svp['viparams']['filebase'], scene.frame_current, '{}-{}.gpm'.format(scene['viparams']['filebase'], frame), cpfileentry, cam.matrix_world.to_quaternion()@ mathutils.Vector((0, 1, 0)), ('', '-i')[self.simnode.illu])
                
            else:
                rvucmd = "rvu -w {9} -n {0} -vv {1} -vh {2} -vd {3[0]:.3f} {3[1]:.3f} {3[2]:.3f} -vp {4[0]:.3f} {4[1]:.3f} {4[2]:.3f} -vu {8[0]:.3f} {8[1]:.3f} {8[2]:.3f} {5} {6}-{7}.oct".format(svp['viparams']['wnproc'], 
                                 vv, cang, vd, cam.location, self.simnode['radparams'], svp['viparams']['filebase'], scene.frame_current, cam.matrix_world.to_quaternion()@ mathutils.Vector((0, 1, 0)), ('', '-i')[self.simnode.illu])

            logentry('Rvu command: {}'.format(rvucmd))
            self.rvurun = Popen(rvucmd.split(), stdout = PIPE, stderr = PIPE)
            self.simnode.run = 1
            wm = context.window_manager
            self._timer = wm.event_timer_add(5, window = context.window)
            wm.modal_handler_add(self)
            return {'RUNNING_MODAL'}

        else:
            self.report({'ERROR'}, "There is no camera in the scene. Radiance preview will not work")
            return {'CANCELLED'}

class NODE_OT_Li_Sim(bpy.types.Operator):
    bl_idname = "node.livicalc"
    bl_label = "LiVi simulation"
    bl_register = True
    bl_undo = False

    def invoke(self, context, event):
        scene = context.scene
        svp = scene.vi_params
        svp.vi_display = 0
        if viparams(self, scene):
            return {'CANCELLED'}
                    
        objmode()
        clearscene(scene, self)
        simnode = context.node
        simnode.presim()
        contextdict = {'Basic': 'LiVi Basic', 'Compliance': 'LiVi Compliance', 'CBDM': 'LiVi CBDM'}        
        
        # Set scene parameters
        svp['viparams']['visimcontext'] = contextdict[simnode['coptions']['Context']]
        svp['liparams']['fs'] = min((simnode['coptions']['fs'], simnode['goptions']['fs'])) 
        svp['liparams']['fe'] = max((simnode['coptions']['fe'], simnode['goptions']['fe'])) 
        svp['liparams']['cp'] = simnode['goptions']['cp']
        svp['liparams']['unit'] = simnode['coptions']['unit']
#        svp['liparams']['metric'] = simnode['coptions']['metric']
        svp['liparams']['type'] = simnode['coptions']['Type']
        scene.frame_start, scene.frame_end = svp['liparams']['fs'], svp['liparams']['fe']
        
        simnode.sim(scene)

        for frame in range(svp['liparams']['fs'], svp['liparams']['fe'] + 1):
            if createradfile(scene, frame, self, simnode) == 'CANCELLED' or createoconv(scene, frame, self, simnode) == 'CANCELLED':
                return {'CANCELLED'}
        
        calcout = li_calc(self, simnode, livisimacc(simnode))

        if calcout == 'CANCELLED':
            return {'CANCELLED'}
        else:
            simnode['reslists'] = calcout
#        if simnode['coptions']['Context'] != 'CBDM' and simnode['coptions']['Context'] != '1':
#            svp.vi_display = 1

        svp['viparams']['vidisp'] = 'li'
        svp['viparams']['resnode'] = simnode.name
        svp['viparams']['restree'] = simnode.id_data.name
        simnode.postsim()
        self.report({'INFO'},"Simulation is finished")
        return {'FINISHED'}
    
class NODE_OT_Li_Im(bpy.types.Operator):
    bl_idname = "node.radimage"
    bl_label = "LiVi Image"
    bl_description = "Generate an image with Rpict"
    bl_register = True
    bl_undo = False

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
        self.kivyrun.kill() 

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
        svp = scene.vi_params
        self.xindex, self.p = 0, 0
        self.cam = scene.camera
        simnode = context.node
        self.fs, self.fe = simnode.retframes()
        self.simnode = simnode
        
        if simnode.camera and bpy.data.cameras.get(simnode.camera):
            self.percent = 0
            self.reslists, self.images = [], []
            self.res = []
            self.rpictfile = os.path.join(scene['viparams']['newdir'], 'rpictprogress')
            self.pmfile = os.path.join(scene['viparams']['newdir'], 'pmprogress')
            simnode.presim()
            svp['liparams']['fs'], svp['liparams']['fe'] =  simnode.retframes()
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
                
            scene.frame_set(svp['liparams']['fs'])
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
            self._timer = wm.event_timer_add(2, window = context.window)
            wm.modal_handler_add(self)                
            return {'RUNNING_MODAL'}
        else:
            self.report({'ERROR'}, "There is no camera in the scene or selected in the node. Create one for rpict image creation")
            return {'FINISHED'}

class NODE_OT_Li_Gl(bpy.types.Operator):            
    bl_idname = "node.liviglare"
    bl_label = "LiVi Glare Node"
    bl_description = "Glare analysis node"
    bl_register = True
    bl_undo = False
    
    def execute(self, context):
        scene = context.scene
        svp = scene.vi_params
        res = []
        reslists = []
        glnode = context.node
        imnode = glnode.inputs[0].links[0].from_node
        glnode.presim()
        
        for i, im in enumerate(imnode['images']):
            glfile = os.path.join(svp['viparams']['newdir'], 'images', '{}-{}.hdr'.format(glnode['hdrname'], i + svp['liparams']['fs']))
            egcmd = 'evalglare {} -c {}'.format(('-u {0[0]} {0[1]} {0[2]}'.format(glnode.gc), '')[glnode.rand], glfile)

            with open(im, 'r') as hdrfile:
                egrun = Popen(egcmd.split(), stdin = hdrfile, stdout = PIPE)

            time = datetime.datetime(2014, 1, 1, imnode['coptions']['shour'], 0) + datetime.timedelta(imnode['coptions']['sdoy'] - 1) if imnode['coptions']['anim'] == '0' else \
                datetime.datetime(2014, 1, 1, int(imnode['coptions']['shour']), int(60*(imnode['coptions']['shour'] - int(imnode['coptions']['shour'])))) + datetime.timedelta(imnode['coptions']['sdoy'] - 1) + datetime.timedelta(hours = int(imnode['coptions']['interval']*i), 
                                  seconds = int(60*(imnode['coptions']['interval']*i - int(imnode['coptions']['interval']*i))))
            
            with open(os.path.join(svp['viparams']['newdir'], 'images', "temp.glare"), "w") as glaretf:
                for line in egrun.stdout:
                    if line.decode().split(",")[0] == 'dgp':                            
                        glaretext = line.decode().replace(',', ' ').replace("#INF", "").split(' ')
                        res = [float(x) for x in glaretext[6:12]]
                        glaretf.write("{0:0>2d}/{1:0>2d} {2:0>2d}:{3:0>2d}\ndgp: {4:.2f}\ndgi: {5:.2f}\nugr: {6:.2f}\nvcp: {7:.2f}\ncgi: {8:.2f}\nLv: {9:.0f}\n".format(time.day, time.month, time.hour, time.minute, *res))
                        res.append(res)
                        reslists += [[str(i + svp['liparams']['fs']), 'Camera', 'Camera', 'DGP', '{0[0]}'.format(res)], 
                                      [str(i + svp['liparams']['fs']), 'Camera', 'Camera', 'DGI', '{0[1]}'.format(res)], 
                                      [str(i + svp['liparams']['fs']), 'Camera', 'Camera' 'UGR', '{0[2]}'.format(res)], 
                                      [str(i + svp['liparams']['fs']), 'Camera', 'Camera', 'VCP', '{0[3]}'.format(res)], 
                                      [str(i + svp['liparams']['fs']), 'Camera', 'Camera', 'CGI', '{[4]}'.format(res)], 
                                      [str(i + svp['liparams']['fs']), 'Camera', 'Camera', 'LV', '{[5]}'.format(res)]]
            
            pcondcmd = "pcond -h+ -u 300 {}.hdr".format(os.path.join(svp['viparams']['newdir'], 'images', '{}-{}'.format(glnode['hdrname'], str(i + svp['liparams']['fs']))))

            with open('{}.temphdr'.format(os.path.join(svp['viparams']['newdir'], 'images', 'glare')), 'w') as temphdr:
                Popen(pcondcmd.split(), stdout = temphdr).communicate()

            catcmd = "{0} {1}.glare".format(svp['viparams']['cat'], os.path.join(svp['viparams']['newdir'], 'images', 'temp'))
            catrun = Popen(catcmd, stdout = PIPE, shell = True)
            psigncmd = "psign -h {} -cb 0 0 0 -cf 1 1 1".format(int(0.04 * imnode.y))
            psignrun = Popen(psigncmd.split(), stdin = catrun.stdout, stdout = PIPE)
            pcompcmd = "pcompos {0}.temphdr 0 0 - {1} {2}".format(os.path.join(svp['viparams']['newdir'], 'images', 'glare'), imnode.x, imnode.y*550/800)

            with open("{}.hdr".format(os.path.join(svp['viparams']['newdir'], 'images', '{}-{}'.format(glnode['hdrname'], str(i + svp['liparams']['fs'])))), 'w') as ghdr:
                Popen(pcompcmd.split(), stdin = psignrun.stdout, stdout = ghdr).communicate()

            os.remove(os.path.join(svp['viparams']['newdir'], 'images', 'glare.temphdr'.format(i + svp['liparams']['fs'])))                               
        return {'FINISHED'}

class NODE_OT_Li_Fc(bpy.types.Operator):            
    bl_idname = "node.livifc"
    bl_label = "LiVi False Colour Image"
    bl_description = "False colour an image with falsecolor"
    bl_register = True
    bl_undo = False

    def execute(self, context):
        fcnode = context.node
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
    
class VIEW3D_OT_Li_BD(bpy.types.Operator):
    '''Display results legend and stats in the 3D View'''
    bl_idname = "view3d.libd"
    bl_label = "LiVi basic metric display"
    bl_description = "Display basic lighting metrics"
    bl_register = True
    bl_undo = False
    
    def modal(self, context, event): 
        scene = context.scene
        svp = scene.vi_params
        redraw = 0 
        
        if svp.vi_display == 0 or svp['viparams']['vidisp'] != 'li' or not [o for o in context.scene.objects if o.name in svp['liparams']['livir']]:
            bpy.types.SpaceView3D.draw_handler_remove(self.draw_handle_linum, 'WINDOW')
            context.area.tag_redraw()
            return {'CANCELLED'}        
        
        if event.type != 'INBETWEEN_MOUSEMOVE' and context.region and context.area.type == 'VIEW_3D' and context.region.type == 'WINDOW':                    
            mx, my = event.mouse_region_x, event.mouse_region_y 
            
#            if any((svp.vi_leg_levels != self.legend.levels, svp.vi_leg_col != self.legend.col, svp.vi_leg_scale != self.legend.scale, (self.legend.minres, self.legend.maxres) != leg_min_max(self.scene))):               
#                self.legend.update(context)                
#                redraw = 1
                 
            # Legend routine 
            
            if self.legend.ispos[0] < mx < self.legend.iepos[0] and self.legend.ah - 80 < my < self.legend.ah - 40:
                self.legend.hl = (0.8, 0.8, 0.8, 0.8) 
                redraw = 1
                if event.type == 'LEFTMOUSE':
                    if event.value == 'RELEASE':
                        self.legend.expand = 0 if self.legend.expand else 1
                        
            elif self.legend.expand and abs(self.legend.lspos[0] - mx) < 10 and abs(self.legend.lepos[1] - my) < 10:
                self.legend.hl = (0.8, 0.8, 0.8, 0.8) 
                redraw = 1   
                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.legend.move = 1
                        self.legend.draw(context)
                        context.area.tag_redraw()
                    elif self.legend.move and event.value == 'RELEASE':
                        self.legend.move = 0                        
                    return {'RUNNING_MODAL'}
              
            elif self.legend.expand and abs(self.legend.lepos[0] - mx) < 10 and abs(self.legend.lspos[1] - my) < 10:
                self.legend.hl = (0.8, 0.8, 0.8, 0.8) 
                context.area.tag_redraw()
                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.legend.resize = 1
                        self.legend.draw(context)
                        context.area.tag_redraw()
                    elif self.legend.resize and event.value == 'RELEASE':
                        self.legend.resize = 0
                    return {'RUNNING_MODAL'}  
                
            elif self.legend.hl == (0.8, 0.8, 0.8, 0.8):                 
                self.legend.hl = (1, 1, 1, 1)
                redraw = 1
                
            if event.type == 'MOUSEMOVE':                
                if self.legend.move:
                    self.legend.lspos[0], self.legend.lepos[1] = mx, my
                    self.legend.draw(context)
                    context.area.tag_redraw() 
                elif self.legend.resize:
                    self.legend.lepos[0], self.legend.lspos[1] = mx, my
                    self.legend.draw(context)
                    context.area.tag_redraw()
            
            if redraw:
                context.area.tag_redraw()

        return {'PASS_THROUGH'}
    
#    def modal(self, context, event): 
#        redraw = 0 
#        
#        if self.scene.vi_display == 0 or context.scene['viparams']['vidisp'] != 'lipanel' or not any([o.lires for o in bpy.data.objects]):
#            self.scene.vi_display = 0
#            bpy.types.SpaceView3D.draw_handler_remove(self._handle_disp, 'WINDOW')
#            bpy.types.SpaceView3D.draw_handler_remove(self._handle_pointres, 'WINDOW')
#            self.scene['viparams']['vidisp'] = 'li'
#            context.area.tag_redraw()
#            return {'CANCELLED'}
#
#        if event.type != 'INBETWEEN_MOUSEMOVE' and context.region and context.area.type == 'VIEW_3D' and context.region.type == 'WINDOW':            
#            mx, my = event.mouse_region_x, event.mouse_region_y 
#            
#            if any((self.scene.vi_leg_levels != self.legend.levels, self.scene.vi_leg_col != self.legend.col, self.scene.vi_leg_scale != self.legend.scale, (self.legend.minres, self.legend.maxres) != leg_min_max(self.scene))):               
#                self.legend.update(context)
#                redraw = 1
#            
#            # Legend routine 
#            
#            if self.legend.spos[0] < mx < self.legend.epos[0] and self.legend.spos[1] < my < self.legend.epos[1]:
#                if self.legend.hl != (0, 1, 1, 1):
#                    self.legend.hl = (0, 1, 1, 1)
#                    redraw = 1  
#                if event.type == 'LEFTMOUSE':
#                    if event.value == 'PRESS':
#                        self.legend.press = 1
#                        self.legend.move = 0
#                        return {'RUNNING_MODAL'}
#                    elif event.value == 'RELEASE':
#                        if not self.legend.move:
#                            self.legend.expand = 0 if self.legend.expand else 1
#                        self.legend.press = 0
#                        self.legend.move = 0
#                        context.area.tag_redraw()
#                        return {'RUNNING_MODAL'}
#                
#                elif event.type == 'ESC':
#                    bpy.types.SpaceView3D.draw_handler_remove(self._handle_disp, 'WINDOW')
#                    bpy.types.SpaceView3D.draw_handler_remove(self._handle_pointres, 'WINDOW')
#                    self.scene['viparams']['vidisp'] = 'li'
#                    context.area.tag_redraw()
#                    return {'CANCELLED'}
#                    
#                elif self.legend.press and event.type == 'MOUSEMOVE':
#                     self.legend.move = 1
#                     self.legend.press = 0
#            
#            elif abs(self.legend.lepos[0] - mx) < 10 and abs(self.legend.lspos[1] - my) < 10:
#                if self.legend.hl != (0, 1, 1, 1):
#                    self.legend.hl = (0, 1, 1, 1)
#                    redraw = 1
#                if event.type == 'LEFTMOUSE':
#                    if event.value == 'PRESS':
#                        self.legend.resize = 1
#                    if self.legend.resize and event.value == 'RELEASE':
#                        self.legend.resize = 0
#                    return {'RUNNING_MODAL'}
#                    
#            elif self.legend.hl != (1, 1, 1, 1):
#                self.legend.hl = (1, 1, 1, 1)
#                redraw = 1

            # Table routine
            
#            if self.frame != context.scene.frame_current or self.table.unit != context.scene['liparams']['unit'] or self.table.cao != context.active_object:
#                self.table.update(context)                
#                redraw = 1
#            
#            if self.table.spos[0] < mx < self.table.epos[0] and self.table.spos[1] < my < self.table.epos[1]:
#                if self.table.hl != (0, 1, 1, 1):
#                    self.table.hl = (0, 1, 1, 1)
#                    redraw = 1  
#                if event.type == 'LEFTMOUSE':
#                    if event.value == 'PRESS':
#                        self.table.press = 1
#                        self.table.move = 0
#                        return {'RUNNING_MODAL'}
#                    elif event.value == 'RELEASE':
#                        if not self.table.move:
#                            self.table.expand = 0 if self.table.expand else 1
#                        self.table.press = 0
#                        self.table.move = 0
#                        context.area.tag_redraw()
#                        return {'RUNNING_MODAL'}
                
#                elif event.type == 'ESC':
#                    bpy.types.SpaceView3D.draw_handler_remove(self._handle_disp, 'WINDOW')
#                    bpy.types.SpaceView3D.draw_handler_remove(self._handle_pointres, 'WINDOW')
#                    context.scene['viparams']['vidisp'] = 'li'
#                    context.area.tag_redraw()
#                    return {'CANCELLED'}
#                    
#                elif self.table.press and event.type == 'MOUSEMOVE':
#                     self.table.move = 1
#                     self.table.press = 0                     
#            elif self.table.hl != (1, 1, 1, 1):
#                self.table.hl = (1, 1, 1, 1)
#                redraw = 1
#                
#            if context.scene['viparams']['visimcontext'] == 'LiVi Compliance':
#                if self.frame != context.scene.frame_current:
#                    self.tablecomp.update(context)
#                    redraw = 1
#                if self.tablecomp.unit != context.scene['liparams']['unit']:
#                    self.tablecomp.update(context)
#                    self.tablecomp.unit = context.scene['liparams']['unit']
#                    redraw = 1
#                if self.tablecomp.cao != context.active_object:
#                    self.tablecomp.update(context)
#                    redraw = 1
#                
#                if self.tablecomp.spos[0] < mx < self.tablecomp.epos[0] and self.tablecomp.spos[1] < my < self.tablecomp.epos[1]:
#                    if self.tablecomp.hl != (0, 1, 1, 1):
#                        self.tablecomp.hl = (0, 1, 1, 1)
#                        redraw = 1  
#                    if event.type == 'LEFTMOUSE':
#                        if event.value == 'PRESS':
#                            self.tablecomp.press = 1
#                            self.tablecomp.move = 0
#                            return {'RUNNING_MODAL'}
#                        elif event.value == 'RELEASE':
#                            if not self.tablecomp.move:
#                                self.tablecomp.expand = 0 if self.tablecomp.expand else 1
#                            self.tablecomp.press = 0
#                            self.tablecomp.move = 0
#                            context.area.tag_redraw()
#                            return {'RUNNING_MODAL'}
#                    
#                    elif event.type == 'ESC':
#                        bpy.types.SpaceView3D.draw_handler_remove(self._handle_disp, 'WINDOW')
#                        bpy.types.SpaceView3D.draw_handler_remove(self._handle_pointres, 'WINDOW')
#                        context.scene['viparams']['vidisp'] = 'li'
#                        context.area.tag_redraw()
#                        return {'CANCELLED'}
#                        
#                    elif self.tablecomp.press and event.type == 'MOUSEMOVE':
#                         self.tablecomp.move = 1
#                         self.tablecomp.press = 0                     
#                elif self.tablecomp.hl != (1, 1, 1, 1):
#                    self.tablecomp.hl = (1, 1, 1, 1)
#                    redraw = 1
                
#            if context.scene['liparams']['unit'] in ('ASE (hrs)', 'sDA (%)', 'DA (%)', 'UDI-f (%)', 'UDI-s (%)', 'UDI-e (%)', 'UDI-a (%)', 'Max lux', 'Min lux', 'Avg lux', 'kWh', 'kWh/m2'):
#                if self.dhscatter.frame != context.scene.frame_current:
#                    self.dhscatter.update(context)
#                    redraw = 1
#                if self.dhscatter.unit != context.scene['liparams']['unit']:
#                    self.dhscatter.update(context)
#                    redraw = 1
#                if self.dhscatter.cao != context.active_object:
#                    self.dhscatter.update(context)
#                    redraw = 1
#                if self.dhscatter.col != context.scene.vi_leg_col:
#                    self.dhscatter.update(context)
#                    redraw = 1
#                if context.scene['liparams']['unit'] in ('Max lux', 'Min lux', 'Avg lux', 'kWh', 'kWh/m2'):
#                    if (self.dhscatter.vmin, self.dhscatter.vmax) != (context.scene.vi_scatter_min, context.scene.vi_scatter_max):
#                       self.dhscatter.update(context) 
#                       redraw = 1
#                                        
#                if self.dhscatter.spos[0] < mx < self.dhscatter.epos[0] and self.dhscatter.spos[1] < my < self.dhscatter.epos[1]:
#                    if self.dhscatter.hl != (0, 1, 1, 1):
#                        self.dhscatter.hl = (0, 1, 1, 1)
#                        redraw = 1 
#                    if event.type == 'LEFTMOUSE':
#                        if event.value == 'PRESS':
#                            self.dhscatter.press = 1
#                            self.dhscatter.move = 0
#                            return {'RUNNING_MODAL'}
#                        elif event.value == 'RELEASE':
#                            if not self.dhscatter.move:
#                                self.dhscatter.expand = 0 if self.dhscatter.expand else 1
#                            self.dhscatter.press = 0
#                            self.dhscatter.move = 0
#                            context.area.tag_redraw()
#                            return {'RUNNING_MODAL'}
                    
#                    elif event.type == 'ESC':
#                        bpy.types.SpaceView3D.draw_handler_remove(self.draw_handle_ssnum, 'WINDOW')
##                        bpy.types.SpaceView3D.draw_handler_remove(self._handle_pointres, 'WINDOW')
#                        context.scene['viparams']['vidisp'] = 'li'
#                        context.area.tag_redraw()
#                        return {'CANCELLED'}
#                        
#                    elif self.dhscatter.press and event.type == 'MOUSEMOVE':
#                         self.dhscatter.move = 1
#                         self.dhscatter.press = 0   
#                                            
#                else:
#                    if self.dhscatter.hl != (1, 1, 1, 1):
#                        self.dhscatter.hl = (1, 1, 1, 1)
#                        redraw = 1
#                    if self.dhscatter.lspos[0] < mx < self.dhscatter.lepos[0] and self.dhscatter.lspos[1] < my < self.dhscatter.lepos[1] and abs(self.dhscatter.lepos[0] - mx) > 20 and abs(self.dhscatter.lspos[1] - my) > 20:
#                        if self.dhscatter.expand: 
#                            self.dhscatter.hl = (1, 1, 1, 1)
#                            if event.type == 'LEFTMOUSE' and event.value == 'PRESS' and self.dhscatter.expand and self.dhscatter.lspos[0] < mx < self.dhscatter.lepos[0] and self.dhscatter.lspos[1] < my < self.dhscatter.lspos[1] + 0.9 * self.dhscatter.ydiff:
#                                self.dhscatter.show_plot()
#                                                         
#            # Resize routines
#            
#            if abs(self.legend.lepos[0] - mx) < 20 and abs(self.legend.lspos[1] - my) < 20:
#                self.legend.hl = (0, 1, 1, 1) 
#                if event.type == 'LEFTMOUSE':
#                    if event.value == 'PRESS':
#                        self.legend.resize = 1
#                    if self.legend.resize and event.value == 'RELEASE':
#                        self.legend.resize = 0
#                    return {'RUNNING_MODAL'}
                    
#            elif abs(self.table.lepos[0] - mx) < 20 and abs(self.table.lspos[1] - my) < 20:
#                self.table.hl = (0, 1, 1, 1) 
#                if event.type == 'LEFTMOUSE':
#                    if event.value == 'PRESS':
#                        self.table.resize = 1
#                    if self.table.resize and event.value == 'RELEASE':
#                        self.table.resize = 0
#                    return {'RUNNING_MODAL'}
#            
#            elif context.scene['viparams']['visimcontext'] == 'LiVi Compliance' and abs(self.tablecomp.lepos[0] - mx) < 20 and abs(self.tablecomp.lspos[1] - my) < 20:
#                self.tablecomp.hl = (0, 1, 1, 1) 
#                if event.type == 'LEFTMOUSE':
#                    if event.value == 'PRESS':
#                        self.tablecomp.resize = 1
#                    if self.tablecomp.resize and event.value == 'RELEASE':
#                        self.tablecomp.resize = 0
#                    return {'RUNNING_MODAL'}
#
#            elif context.scene['liparams']['unit'] in ('ASE (hrs)', 'sDA (%)', 'DA (%)', 'UDI-s (%)', 'UDI-e (%)', 'UDI-f (%)', 'UDI-a (%)', 'Max lux', 'Min lux', 'Avg lux', 'kWh', 'kWh/m2') and abs(self.dhscatter.lepos[0] - mx) < 20 and abs(self.dhscatter.lspos[1] - my) < 20:
#                self.dhscatter.hl = (0, 1, 1, 1) 
#                if event.type == 'LEFTMOUSE':
#                    if event.value == 'PRESS':
#                        self.dhscatter.resize = 1
#                    if self.dhscatter.resize and event.value == 'RELEASE':
#                        self.dhscatter.resize = 0
#                    return {'RUNNING_MODAL'}
            # Move routines
                     
#            if event.type == 'MOUSEMOVE':                
#                if self.legend.move:
#                    self.legend.pos = [mx, my]
#                    redraw = 1
#                if self.legend.resize:
#                    self.legend.lepos[0], self.legend.lspos[1] = mx, my
#                    redraw = 1
#                if self.table.move:
#                    self.table.pos = [mx, my]
#                    redraw = 1
#                if self.table.resize:
#                    self.table.lepos[0], self.table.lspos[1] = mx, my
#                    redraw = 1
#                if context.scene['viparams']['visimcontext'] == 'LiVi Compliance':
#                    if self.tablecomp.move:
#                        self.tablecomp.pos = [mx, my]
#                        redraw = 1
#                    if self.tablecomp.resize:
#                        self.tablecomp.lepos[0], self.tablecomp.lspos[1] = mx, my
#                        redraw = 1
#                try:
#                    if self.dhscatter.move:
#                        self.dhscatter.pos = [mx, my]
#                        redraw = 1
#                    if self.dhscatter.resize:
#                        self.dhscatter.lepos[0], self.dhscatter.lspos[1] = mx, my
#                        redraw = 1
#                except:
#                    pass
#                                
#            if redraw:
#                context.area.tag_redraw()
#                self.frame = context.scene.frame_current
#                
#        return {'PASS_THROUGH'}

    def invoke(self, context, event):
        area = context.area
        self.scene = context.scene
        svp = context.scene.vi_params
        svp.vi_display, svp.vi_disp_wire = 1, 1        
        clearscene(self.scene, self)
        svp['viparams']['vidisp'] = 'li' 
        self.simnode = bpy.data.node_groups[svp['viparams']['restree']].nodes[svp['viparams']['resnode']]
        
        self.images = ['legend.png']#, 'table_new.png']
#        if svp['viparams']['visimcontext'] == 'LiVi Compliance':
#            self.images.append('compliance.png')
#            if self.simnode['coptions']['canalysis'] == '3':
#                self.images.append('scatter.png')
#        elif svp['viparams']['visimcontext'] == 'LiVi CBDM':
#            self.images.append('scatter.png')
        self.results_bar = results_bar(self.images, 300, area)
        
#        self.scene['viparams']['vidisp'] = 'lipanel'
        self.frame = self.scene.frame_current
        
        if li_display(self, self.simnode) == 'CANCELLED':
            return {'CANCELLED'}

        lnd = linumdisplay(self, context)
   #     self._handle_pointres = bpy.types.SpaceView3D.draw_handler_add(lnd.draw, (context, ), 'WINDOW', 'POST_PIXEL')
        self.legend = livi_legend(context, svp['liparams']['unit'], [305, area.height - 80], area.width, area.height, 125, 400)

#        self.dhscatter = wr_scatter([160, context.region.height - 40], context.region.width, context.region.height, 'stats.png', 600, 400)

#        if svp['viparams']['visimcontext'] == 'LiVi Basic':
#            self.table = basic_table([240, context.region.height - 40], context.region.width, context.region.height, 'table.png', 600, 100)  
#
#        if svp['viparams']['visimcontext'] == 'LiVi Compliance':
#            self.table_comp = comp_table([300, context.region.height - 40], context.region.width, context.region.height, 'compliance.png', 600, 200)
#
#            if self.simnode['coptions']['canalysis'] == '3':
#                self.dhscatter = leed_scatter([160, context.region.height - 40], context.region.width, context.region.height, 'scat.png', 600, 400)
#                self.dhscatter.update(context)        
#            self._handle_disp = bpy.types.SpaceView3D.draw_handler_add(comp_disp, (self, context, self.simnode), 'WINDOW', 'POST_PIXEL')

#        self.dhscatter.update(context)
        
#        self._handle_spnum = bpy.types.SpaceView3D.draw_handler_add(viwr_legend, (self, context, simnode), 'WINDOW', 'POST_PIXEL')
#        elif self.scene['viparams']['visimcontext'] == 'LiVi Basic':
#            self._handle_disp = bpy.types.SpaceView3D.draw_handler_add(basic_disp, (self, context, self.simnode), 'WINDOW', 'POST_PIXEL')
#        if self.scene['viparams']['visimcontext'] == 'LiVi Compliance':
#            self._handle_disp = bpy.types.SpaceView3D.draw_handler_add(comp_disp, (self, context, self.simnode), 'WINDOW', 'POST_PIXEL')
#        elif self.scene['viparams']['visimcontext'] == 'LiVi CBDM':
#            if self.simnode['coptions']['cbanalysis'] != '0':
#                self.dhscatter = cbdm_scatter([160, context.region.height - 40], context.region.width, context.region.height, 'scat.png', 600, 400)
#                self.dhscatter.update(context)
#            self._handle_disp = bpy.types.SpaceView3D.draw_handler_add(cbdm_disp, (self, context, self.simnode), 'WINDOW', 'POST_PIXEL')
        self.draw_handle_linum = bpy.types.SpaceView3D.draw_handler_add(self.draw_linum, (context, ), 'WINDOW', 'POST_PIXEL')        
#     
        context.window_manager.modal_handler_add(self)
        return {'RUNNING_MODAL'}

    def draw_linum(self, context):
        area = context.area
        ah = area.height                
        self.results_bar.draw(ah)
        self.legend.draw(context)
#        self.dhscatter.draw(context, area.width)
#        self.num_display.draw(context)

        
class MAT_EnVi_Node(bpy.types.Operator):
    bl_idname = "material.envi_node"
    bl_label = "EnVi Material export"
    nodeid: bpy.props.StringProperty()

    def invoke(self, context, event):
        cm = context.material
        mvp = cm.vi_params
        if not mvp.envi_nodes:
            bpy.ops.node.new_node_tree(type='EnViMatN', name = cm.name) 
            mvp.envi_nodes = bpy.data.node_groups[cm.name]
            mvp.envi_nodes.nodes.new('No_En_Mat_Con')
            mvp.envi_nodes['envi_con_type'] = 'None'
            mvp.envi_nodes.nodes[0].active = True
            mvp.envi_nodes['enmatparams'] = {'airflow': 0, 'boundary': 0, 'tm': 0}

        elif cm.name != mvp.envi_nodes.name and mvp.envi_nodes.name in [m.name for m in bpy.data.materials]:
            mvp.envi_nodes = mvp.envi_nodes.copy()
            mvp.envi_nodes.name = cm.name
        return {'FINISHED'}

class NODE_OT_En_Geo(bpy.types.Operator):
    bl_idname = "node.engexport"
    bl_label = "VI-Suite export"
    bl_context = "scene"

    def invoke(self, context, event):
        scene = context.scene
        svp = scene.vi_params
        
        if viparams(self, scene):
            return {'CANCELLED'}
        
        svp['viparams']['vidisp'] = ''
        node = context.node
        node.preexport(scene)
        pregeo(context, self)
        node.postexport()
        return {'FINISHED'}
    
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