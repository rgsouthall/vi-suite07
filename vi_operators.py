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
#from .livi_export import radgexport, spfc, createoconv, createradfile
#from .livi_calc  import li_calc
#from .vi_display import li_display, linumdisplay, spnumdisplay, en_air, wr_legend, wr_disp, wr_scatter, wr_table, ss_disp, ss_legend, svf_disp, svf_legend, basic_legend, basic_table, basic_disp, ss_scatter, en_disp, en_pdisp, en_scatter, en_table, en_barchart, comp_table, comp_disp, leed_scatter, cbdm_disp, cbdm_scatter, envals, bsdf, bsdf_disp#, en_barchart, li3D_legend
#from .envi_export import enpolymatexport, pregeo
#from .envi_mat import envi_materials, envi_constructions
from .vi_func import selobj, joinobj, solarPosition, viparams, compass, spfc
#from .flovi_func import fvcdwrite, fvbmwrite, fvblbmgen, fvvarwrite, fvsolwrite, fvschwrite, fvtppwrite, fvraswrite, fvshmwrite, fvmqwrite, fvsfewrite, fvobjwrite, fvdcpwrite
from .vi_func import spathrange, ret_plt
from .vi_func import sunpath
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
    bl_undo = True
#    nodeid = bpy.props.StringProperty()

    def invoke(self, context, event):
        scene = context.scene
#        print(dir(scene))
        if viparams(self, scene):
            self.report({'ERROR'},"Save the Blender file before continuing")
            return {'CANCELLED'}
        
        try:
            spcoll = bpy.data.collections['SunPath']
        except:
            spcoll = bpy.data.collections.new('SunPath')
            context.scene.collection.children.link(spcoll)
            
        for lcc in context.view_layer.layer_collection.children:
            if lcc.name == 'SunPath':
                context.view_layer.active_layer_collection = lcc
            
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
            


        requiredsuns = {'0': 1, '1': 12, '2': 24}[node.suns]

        matdict = {'SolEquoRings': (1, 0, 0, 1), 'HourRings': (1, 1, 0, 1), 'SPBase': (1, 1, 1, 1), 'Sun': (1, 1, 1, 1), 'PathDash': (1, 1, 1, 1),
                   'SumAng': (1, 0, 0, 1), 'EquAng': (0, 1, 0, 1), 'WinAng': (0, 0, 1, 1)}
        
        for mat in [mat for mat in matdict if mat not in bpy.data.materials]:
            bpy.data.materials.new(mat)
            bpy.data.materials[mat].diffuse_color = matdict[mat][:4]
#            bpy.data.materials[mat].use_shadeless = 1
            bpy.data.materials[mat].use_nodes = True
            nodes = bpy.data.materials[mat].node_tree.nodes

            for n in nodes:
                nodes.remove(n)
            
            if mat == 'PathDash':
#                bpy.data.materials[mat].diffuse_color[3] = 0
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
                bpy.ops.object.light_add(type='SUN', radius=1, align='WORLD', location=(0, 0, 0))

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
            if spathob.name in scene.collection.objects:
                scene.collection.objects.unlink(spathob)
#        scene.collection.objects.unlink(spathob)
        spathob.location, spathob.name,  spathob['VIType'] = (0, 0, 0), "SPathMesh", "SPathMesh"
#        smesh.parent = spathob
        
        for s, sun in enumerate(suns):
            if sun.name not in spcoll.objects:
                spcoll.objects.link(sun)
                scene.collection.objects.unlink(sun)
                
            sun.data.shadow_soft_size = 0.01            
            sun['VIType'] = 'Sun'
            sun['solhour'], sun['solday'] = scene.solhour, scene.solday
            sun.name = sun.data.name ='Sun{}'.format(s)
#            bpy.ops.mesh.primitive_uv_sphere_add(segments=12, ring_count=12, radius=0.5)
#            sunob = context.active_object
            
#            if sunob.name not in spcoll.objects:
#                spcoll.objects.link(sunob)
#                scene.collection.objects.unlink(sunob)
#            scene.collection.objects.unlink(sunob)
#            sunob.location, sunob.cycles_visibility.shadow, sunob.name, sunob['VIType'] = (0, 0, 0), 0, "SunMesh{}".format(s), "SunMesh"
#            sunob.cycles_visibility.diffuse, sunob.cycles_visibility.shadow, sunob.cycles_visibility.glossy, sunob.cycles_visibility.transmission, sunob.cycles_visibility.scatter = [False] * 5

#            if len(sunob.material_slots) == 0:
#                 bpy.ops.object.material_slot_add()
#                 sunob.material_slots[0].material = bpy.data.materials['Sun']
#                 
            sun.parent = spathob
#            sunob.parent = sun
        
#        bm = bmesh.new()
#        bm.from_mesh(spathmesh)

#        for doy in range(0, 365, 2):
#            for hour in range(1, 25):
#                ([solalt, solazi]) = solarPosition(doy, hour, scene.latitude, scene.longitude)[2:]
#                bm.verts.new().co = [-(sd-(sd-(sd*cos(solalt))))*sin(solazi), -(sd-(sd-(sd*cos(solalt))))*cos(solazi), sd*sin(solalt)]
#        
#        if hasattr(bm.verts, "ensure_lookup_table"):
#            bm.verts.ensure_lookup_table()
#        for v in range(24, len(bm.verts)):
#            bm.edges.new((bm.verts[v], bm.verts[v - 24]))
#        if v in range(8568, 8761):
#            bm.edges.new((bm.verts[v], bm.verts[v - 8568]))
#
#        for doy in (79, 172, 355):
#            for hour in range(1, 241):
#                ([solalt, solazi]) = solarPosition(doy, hour*0.1, scene.latitude, scene.longitude)[2:]
#                vcoord = [-(sd-(sd-(sd*cos(solalt))))*sin(solazi), -(sd-(sd-(sd*cos(solalt))))*cos(solazi), sd*sin(solalt)]
#                bm.verts.new().co = vcoord
#                if hasattr(bm.verts, "ensure_lookup_table"):
#                    bm.verts.ensure_lookup_table()
#                if bm.verts[-1].co.z >= 0 and doy in (172, 355) and not hour%10:
#                    numpos['{}-{}'.format(doy, int(hour*0.1))] = vcoord
#                if hour != 1:
#                    bm.edges.new((bm.verts[-2], bm.verts[-1]))
#                    solringnum += 1
#                if hour == 240:
#                    bm.edges.new((bm.verts[-240], bm.verts[-1]))
#                    solringnum += 1
        
#        bm.to_mesh(spathmesh)
#        bm.free()
        selobj(context.view_layer, spathob)
#        bpy.ops.object.convert(target='CURVE')
#        spathob.data.bevel_depth, spathob.data.bevel_resolution = node.th, node.res
#        bpy.context.object.data.fill_mode = 'FULL'
#        bpy.ops.object.convert(target='MESH')
        
#        bpy.ops.object.material_slot_add()
#        spathob.material_slots[0].material, spathob['numpos'] = bpy.data.materials['HourRings'], numpos
#        bpy.ops.object.material_slot_add()
#        spathob.material_slots[1].material = bpy.data.materials['PathDash']
#
#        for face in spathob.data.polygons:
#            face.material_index = 0 if not int(face.index/16)%2 else 1
#                
#        for vert in spathob.data.vertices[0:16 * (solringnum + 3)]:
#            vert.select = True
#
#        bpy.ops.object.material_slot_add()
#        spathob.material_slots[-1].material = bpy.data.materials['SolEquoRings']
#        spathob.active_material_index = 2
#        bpy.ops.object.mode_set(mode='EDIT')
#        bpy.ops.mesh.select_mode(type="VERT")
#        bpy.ops.object.material_slot_assign()
#        bpy.ops.mesh.select_all(action='SELECT')
#        bpy.ops.mesh.bisect(plane_co=(0.0, 0.0, 0.0), plane_no=(0.0, 0.0, 1.0), use_fill=True, clear_inner=True, clear_outer=False)
#        bpy.ops.object.mode_set(mode='OBJECT')
#        bpy.ops.object.select_all(action='DESELECT')
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
    
    def ret_coords(self, scene, node):
        breaks, coords, sd, d, line_lengths = [], [], 100, 0, [0]
                
        for hour in range(1, 25):
            for doy in range(0, 365, 2):
                ([solalt, solazi]) = solarPosition(doy, hour, scene.vi_params.latitude, scene.vi_params.longitude)[2:]
                coords.append(Vector([-(sd-(sd-(sd*cos(solalt))))*sin(solazi), -(sd-(sd-(sd*cos(solalt))))*cos(solazi), sd*sin(solalt)]))
                if d%183 == 0:
                    breaks.append(1)
                else:
                    breaks.append(0)
                d += 1

        for doy in (79, 172, 355):
            for hour in range(1, 241):
                ([solalt, solazi]) = solarPosition(doy, hour*0.1, scene.vi_params.latitude, scene.vi_params.longitude)[2:]
                coords.append(Vector([-(sd-(sd-(sd*cos(solalt))))*sin(solazi), -(sd-(sd-(sd*cos(solalt))))*cos(solazi), sd*sin(solalt)]))
                breaks.append(2)
                    
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
#            }
#        uniform sampler2D tex0;
#        uniform float border; // 0.01
#        uniform float circle_radius; // 0.5
#        uniform vec4 sun_colour; 
#        uniform vec2 circle_center; // vec2(0.5, 0.5)    
#        void main(void)
#        {
#          vec2 uv = gl_TexCoord[0].xy;
#          
#          vec4 bkg_color = texture2D(tex0,uv * vec2(1.0, -1.0));
#        
#          // Offset uv with the center of the circle.
#          uv -= circle_center;
#          
#          float dist =  sqrt(dot(uv, uv));
#          if (dist < (circle_radius-border)) {gl_FragColor = sun_colour;}
#          else {discard;}
#        }
#        '''
        
#        sp_vertex_shader = '''
#            in float longitude;
#            in float latitude;
#            uniform mat4 basematrix;
#            
#            void
#        '''
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

        
#        sun_vertex_shader = '''
#        '''
#        sun_geometry_shader = '''
#        '''
#        sun_fragment_shader = '''
#        '''
        
#        breaks, coords, sd, d = [], [], 100, 0
#        
#        
#        for hour in range(1, 25):
#            for doy in range(0, 365, 2):
#                ([solalt, solazi]) = solarPosition(doy, hour, scene.latitude, scene.longitude)[2:]
#                coords.append(Vector([-(sd-(sd-(sd*cos(solalt))))*sin(solazi), -(sd-(sd-(sd*cos(solalt))))*cos(solazi), sd*sin(solalt)]))
#                if d%183 == 0:
#                    breaks.append(1)
#                else:
#                    breaks.append(0)
#                d += 1
#
##        blen = len(breaks)
#                
#        for doy in (79, 172, 355):
#            for hour in range(1, 241):
#                ([solalt, solazi]) = solarPosition(doy, hour*0.1, scene.latitude, scene.longitude)[2:]
#                coords.append(Vector([-(sd-(sd-(sd*cos(solalt))))*sin(solazi), -(sd-(sd-(sd*cos(solalt))))*cos(solazi), sd*sin(solalt)]))
#                breaks.append(2)
#            
#        line_lengths = [0]
#        
#        for a, b in zip(coords[:-1], coords[1:]):
#            line_lengths.append(line_lengths[-1] + (a - b).length)        
        
        self.sp_shader = gpu.types.GPUShader(sp_vertex_shader, sp_fragment_shader) 
        self.sun_shader = gpu.types.GPUShader(sun_vertex_shader, sun_fragment_shader)  
#        self.sun2_shader = gpu.types.GPUShader(sun2_vertex_shader, sun2_fragment_shader) 
        (coords, line_lengths, breaks) = self.ret_coords(scene, node)
#        print(coords, [Vector([so for so in scene.objects if so.type == 'LIGHT' and so.data.type == 'SUN'][0].location[:])])
        sun_pos = [so.location[:] for so in scene.objects if so.type == 'LIGHT' and so.data.type == 'SUN' and not so.hide_viewport]
#        sun_pos = [item for so in scene.objects if so.type == 'LIGHT' and so.data.type == 'SUN' and not so.hide_viewport for item in so.location]

        
        
        sun_v_coords, sun_f_indices = self.ret_sun_geometry(scene.vi_params.sp_sun_size, self.suns)
        self.sp_batch = batch_for_shader(self.sp_shader, 'LINE_STRIP', {"position": coords, "arcLength": line_lengths, "line_break": breaks})
        self.sun_batch = batch_for_shader(self.sun_shader, 'POINTS', {"position": sun_pos})
#        self.sun2_batch = batch_for_shader(self.sun2_shader, 'TRIS', {"position": sun_v_coords, "sun_position": sun_pos}, indices=sun_f_indices)
        
    def draw_sp(self, op, context, node):
        # Draw lines
        bgl.glEnable(bgl.GL_BLEND)
        bgl.glEnable(bgl.GL_LINE_SMOOTH)
        bgl.glEnable(bgl.GL_CULL_FACE)
        bgl.glCullFace(bgl.GL_BACK)
#        bgl.glEnable(bgl.GL_POLYGON_SMOOTH)
        bgl.glLineWidth(context.scene.vi_params.sp_line_width)
        bgl.glPointSize(context.scene.vi_params.sp_sun_size)
        
        self.sp_shader.bind()
        matrix = bpy.context.region_data.perspective_matrix
        sp_matrix = context.scene.objects['SPathMesh'].matrix_world
#        sun_pos = [[l for l in so.location] for so in context.scene.objects if so.type == 'LIGHT' and so.data.type == 'SUN' and not so.hide_viewport]
#        [s.location for s in so.location for so in context.scene.objects if so.type == 'LIGHT' and so.data.type == 'SUN' and not so.hide_viewport]
#        sun_pos = [item for so in context.scene.objects if so.type == 'LIGHT' and so.data.type == 'SUN' and not so.hide_viewport for item in so.location]
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
#        self.sun2_shader.bind()
#        self.sun2_shader.uniform_float("viewProjectionMatrix", matrix)
#        self.sun2_shader.uniform_float("invviewProjectionMatrix", matrix.inverted())
#        self.sun2_shader.uniform_float("sp_matrix", sp_matrix)
#        self.sun2_shader.uniform_float("invsp_matrix", sp_matrix.inverted())
#        self.sun2_shader.uniform_float("sun_colour", context.scene.vi_params.sp_sun_colour)

#        self.sun2_shader.uniform_float("sun_position", sun_pos)
        
        viewport = (0, 265, 1468, 746)        
#        self.sun2_shader.uniform_float("viewport", viewport)
        
        if self.latitude != context.scene.vi_params.latitude or self.longitude != context.scene.vi_params.longitude or \
            self.sd != context.scene.vi_params.sp_sd or self.sh != context.scene.vi_params.sp_sh or self.ss != context.scene.vi_params.sp_sun_size:
            (coords, line_lengths, breaks) = self.ret_coords(context.scene, node)        
            self.sp_batch = batch_for_shader(self.sp_shader, 'LINE_STRIP', {"position": coords, "arcLength": line_lengths, "line_break": breaks})
#            sun_pos = [Vector([so for so in context.scene.objects if so.type == 'LIGHT' and so.data.type == 'SUN'][0].location[:])]
            self.sun_batch = batch_for_shader(self.sun_shader, 'POINTS', {"position": sun_pos})
            sun_v_coords, sun_f_indices = self.ret_sun_geometry(context.scene.vi_params.sp_sun_size, self.suns)
#            self.sun2_batch = batch_for_shader(self.sun2_shader, 'TRIS', {"position": sun_v_coords}, indices=sun_f_indices)
            
#            self.sun2_batch = batch_for_shader(self.sun2_shader, 'TRIS', {"sun_position": sun_pos})
            self.latitude = context.scene.vi_params.latitude
            self.longitude = context.scene.vi_params.longitude
            self.sd = context.scene.vi_params.sp_sd
            self.sh = context.scene.vi_params.sp_sh
            self.ss = context.scene.vi_params.sp_sun_size
        self.sun_batch.draw(self.sun_shader)
        self.sp_batch.draw(self.sp_shader)
#        self.sun2_batch.draw(self.sun2_shader)
        bgl.glDisable(bgl.GL_LINE_SMOOTH)
        bgl.glDisable(bgl.GL_BLEND)
        bgl.glDisable(bgl.GL_CULL_FACE)
#        bgl.glDisable(bgl.GL_POLYGON_SMOOTH)
        bgl.glPointSize(1)

    def modal(self, context, event):
        scene = context.scene
       
        if context.area:
            context.area.tag_redraw()
        if scene.vi_display == 0 or scene['viparams']['vidisp'] != 'sp':
            bpy.types.SpaceView3D.draw_handler_remove(self._handle_spnum, 'WINDOW')
            [bpy.data.objects.remove(o, do_unlink=True, do_id_user=True, do_ui_user=True) for o in bpy.data.objects if o.get('VIType') and o['VIType'] in ('SunMesh', 'SkyMesh')]
            return {'CANCELLED'}
        return {'PASS_THROUGH'}
    
    def ret_sun_geometry(self, dia, suns):
        sun_v_coords, sun_f_indices = [], []
        for sun in suns:
            sunbm = bmesh.new()
            bmesh.ops.create_uvsphere(sunbm, u_segments = 12, v_segments = 12, diameter = dia, matrix = sun.matrix_world, calc_uvs = 0)
            bmesh.ops.triangulate(sunbm, faces = sunbm.faces, quad_method = 'BEAUTY', ngon_method = 'BEAUTY')
#            sun_coords.append([(v.co for v in face.verts) for face in sunbm.faces])
            sun_v_coords += [v.co[:] for v in sunbm.verts]
            
#            self.sun_f_indices = [v.index for v in f.verts for f in sunbm.faces]
#            self.sun_f_indices = [item.index for sublist in sunbm.faces for item in sublist.verts]
            sun_f_indices += [[v.index for v in face.verts] for face in sunbm.faces]
            sunbm.free()
        return sun_v_coords, sun_f_indices

    def invoke(self, context, event):        
        scene = context.scene
        node = context.node
        self.suns = [sun for sun in scene.objects if sun.type == "LIGHT" and sun.data.type == 'SUN']
        self.sp = scene.objects['SPathMesh']
        self.latitude = scene.vi_params.latitude
        self.longitude = scene.vi_params.longitude
        self.sd = scene.vi_params.sp_sd
        self.sh = scene.vi_params.sp_sh
        self.ss = scene.vi_params.sp_sun_size
        self.create_batch(scene, node)
        self.draw_handle_spnum = bpy.types.SpaceView3D.draw_handler_add(self.draw_sp, (self, context, node), "WINDOW", "POST_VIEW")

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
        
