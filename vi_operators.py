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

import bpy, datetime, mathutils, os, bmesh, shutil, sys, math, shlex, itertools
import subprocess
import numpy
from numpy import arange, histogram, array, int8, float16, empty, uint8, transpose, where, delete, ndarray
from bpy_extras.io_utils import ExportHelper, ImportHelper
from subprocess import Popen, PIPE, call
from collections import OrderedDict
from datetime import datetime as dt
from math import cos, sin, pi, ceil, tan, radians
from time import sleep
from mathutils import Euler, Vector, Matrix
from xml.dom.minidom import parse, parseString
from .livi_export import radgexport, createoconv, createradfile, gen_octree, radpoints
#from .livi_calc  import li_calc
from .envi_export import enpolymatexport, pregeo
from .envi_mat import envi_materials, envi_constructions, envi_embodied
from .vi_func import selobj, joinobj, solarPosition, viparams, wind_compass, livisimacc
from .flovi_func import ofheader, fvcdwrite, fvbmwrite, fvblbmgen, fvvarwrite, fvsolwrite, fvschwrite, fvtpwrite, fvmtwrite
from .flovi_func import  fvshmwrite, fvmqwrite, fvsfewrite, fvobjwrite, fvdcpwrite, write_ffile, write_bound, fvtppwrite, fvgwrite, fvrpwrite, fvprefwrite, oftomesh
from .vi_func import ret_plt, logentry, rettree, cmap, fvprogressfile, fvprogressbar
from .vi_func import windnum, wind_rose, create_coll, create_empty_coll, move_to_coll, retobjs, progressfile, progressbar
from .vi_func import chunks, clearlayers, clearscene, clearfiles, objmode, clear_coll, actselobj
from .livi_func import retpmap
from .vi_chart import chart_disp, hmchart_disp
from .vi_dicts import rvuerrdict, pmerrdict
from PyQt5.QtGui import QImage, QColor

try:
    import netgen
    from netgen.meshing import MeshingParameters, FaceDescriptor, Element2D, Mesh
    from netgen.stl import STLGeometry
    from pyngcore import SetNumThreads, TaskManager

except Exception as e:
    print(e)

try:    
    import matplotlib
    matplotlib.use('qt5agg', force = True)
    import matplotlib.cm as mcm
    import matplotlib.colors as mcolors
    mp = 1    
except Exception as e:
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

class NODE_OT_ASCImport(bpy.types.Operator, ImportHelper):
    bl_idname = "node.ascimport"
    bl_label = "Select ESRI Grid file"
    bl_description = "Select the ESRI Grid file to process"
    filename = ""
    filename_ext = ".asc"
    filter_glob: bpy.props.StringProperty(default="*.asc", options={'HIDDEN'})
    bl_register = True
    bl_undo = False
    
    def draw(self,context):
        layout = self.layout
        row = layout.row()
        row.label(text="Open an asc file with the file browser", icon='WORLD_DATA')

    def execute(self, context):
        scene = context.scene
        asccoll = create_coll(context, 'Terrain')
        startxs, startys, vlen = [], [], 0
        ascfiles = [self.filepath] if self.node.single else [os.path.join(os.path.dirname(os.path.realpath(self.filepath)), file) for file in os.listdir(os.path.dirname(os.path.realpath(self.filepath))) if file.endswith('.asc')]
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
            
            if self.node.clear_nodata == '1':
                bmesh.ops.delete(bm, geom = [v for v in bm.verts if v.co[2] == headerdict['NODATA_value']], context = 1)
            
            elif self.node.clear_nodata == '0':
                for v in bm.verts:
                    if v.co[2] == headerdict['NODATA_value']:
                        v.co[2] = 0
                        
            bm.to_mesh(me)
            bm.free()
            ob = bpy.data.objects.new(basename, me)

            if ob.name not in asccoll.objects:
                asccoll.objects.link(ob)
                if ob.name in scene.collection.objects:
                    scene.collection.objects.unlink(ob)

            obs.append(ob)

        minstartx,  minstarty = min(startxs), min(startys)

        for o, ob in enumerate(obs):
            ob.location = (startxs[o] - minstartx, startys[o] - minstarty, 0)
            
        return {'FINISHED'}
        
    def invoke(self,context,event):
        self.node = context.node
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}
    
class NODE_OT_TextUpdate(bpy.types.Operator):
    bl_idname = "node.textupdate"
    bl_label = "Update a text file"
    bl_description = "Update a text file"

    def execute(self, context):
        tenode = context.node
        tenode.textupdate(tenode['bt'])
        return {'FINISHED'}
            
class NODE_OT_WindRose(bpy.types.Operator):
    bl_idname = "node.windrose"
    bl_label = "Wind Rose"
    bl_description = "Create a Wind Rose"
    bl_register = True
    bl_undo = True

    def invoke(self, context, event):
        scene = context.scene
        svp = scene.vi_params
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
        svp['viparams']['resnode'], svp['viparams']['restree'] = simnode.name, simnode.id_data.name
        svp['viparams']['vidisp'], svp.vi_display = 'wr', 0
        svp['viparams']['visimcontext'] = 'Wind'
        rl = locnode['reslists']
        cdoys = [float(c) for c in [r[4].split() for r in rl if r[0] == '0' and r[1] == 'Time' and r[2] == '' and r[3] == 'DOS'][0]]
        cwd = [float(c) for c in [r[4].split() for r in rl if r[0] == '0' and r[1] == 'Climate' and r[2] == '' and r[3] == 'Wind Direction (deg)'][0]]

        if simnode.temp:
            cd = [float(c) for c in [r[4].split() for r in rl if r[0] == '0' and r[1] == 'Climate' and r[2] == '' and r[3] == 'Temperature (degC)'][0]]  
        else:
            cd = [float(c) for c in [r[4].split() for r in rl if r[0] == '0' and r[1] == 'Climate' and r[2] == '' and r[3] == 'Wind Speed (m/s)'][0]] 

        doys = list(range(simnode.sdoy, simnode.edoy + 1)) if simnode.edoy > simnode.sdoy else list(range(1, simnode.edoy + 1)) + list(range(simnode.sdoy, 366))
        awd = array([wd for di, wd in enumerate(cwd) if cdoys[di] in doys])
        ad = array([d for di, d in enumerate(cd) if cdoys[di] in doys])
        validdata = where(awd > 0) if max(cwd) == 360 else where(awd > -1)
        vawd = awd[validdata]
        vad = ad[validdata]
        simnode['maxres'], simnode['minres'], simnode['avres'] = max(cd), min(cd), sum(cd)/len(cd)    
        maxf = simnode.max_freq_val if simnode.max_freq == '1' else 0
        sbinvals = arange(0,int(ceil(max(vad))),2)
        dbinvals = arange(-11.25,372.25,22.5)
        dfreq = histogram(awd, bins=dbinvals)[0]
        adfreq = histogram(cd, bins=dbinvals)[0]
        dfreq[0] = dfreq[0] + dfreq[-1]
        dfreq = dfreq[:-1]  
        fig = plt.figure(figsize=(8, 8), dpi=150, facecolor='w', edgecolor='w')
        rect = [0.1, 0.1, 0.8, 0.8]
        ax = WindroseAxes(fig, rect, facecolor='w')
        fig.add_axes(ax)      
        
        if simnode.wrtype == '0':
            ax.bar(vawd, vad, bins=sbinvals, normed=True, opening=0.8, edgecolor='white', cmap=mcm.get_cmap(svp.vi_scatt_col))
        elif simnode.wrtype == '1':
            ax.box(vawd, vad, bins=sbinvals, normed=True, cmap=mcm.get_cmap(svp.vi_scatt_col))
        elif simnode.wrtype in ('2', '3', '4'):
            ax.contourf(vawd, vad, bins=sbinvals, normed=True, cmap=mcm.get_cmap(svp.vi_scatt_col))
        
        if simnode.max_freq == '1':
            ax.set_rmax(simnode.max_freq_val)
        else:
            ax.set_rmax(100*numpy.max(dfreq)/len(awd) + 0.5)
 
        plt.savefig(svp['viparams']['newdir']+'/disp_wind.svg') 
        wrme = bpy.data.meshes.new("Wind_rose")   
        wro = bpy.data.objects.new('Wind_rose', wrme) 
        
        if wro.name not in wrcoll.objects:
            wrcoll.objects.link(wro)
            if wro.name in scene.collection.objects:
                scene.collection.objects.unlink(wro)
                
        selobj(context.view_layer, wro)    
            
        (wro, scale) = wind_rose(wro, (simnode['maxres'], simnode.max_freq_val)[simnode.max_freq == '1'], svp['viparams']['newdir']+'/disp_wind.svg', simnode.wrtype, mcolors)
        
        wro = joinobj(context.view_layer, wro)  
        ovp = wro.vi_params
        ovp['maxres'], ovp['minres'], ovp['avres'], ovp['nbins'], ovp['VIType'] = max(ad), min(ad), sum(ad)/len(ad), len(sbinvals), 'Wind_Plane'
        simnode['maxfreq'] = 100*numpy.max(dfreq)/len(awd)
        windnum((100*numpy.max(dfreq)/len(awd) + 0.5, simnode.max_freq_val)[simnode.max_freq == '1'], (0,0,0), scale, wind_compass((0,0,0), scale, wro, wro.data.materials['wr-000000']))        
        plt.close()        
        ovp['table'] = array([["", 'Minimum', 'Average', 'Maximum'], 
                             [('Speed (m/s)', 'Temperature (C)')[simnode.temp], ovp['minres'], '{:.1f}'.format(ovp['avres']), ovp['maxres']], 
                             ['Direction (\u00B0)', min(awd), '{:.1f}'.format(sum(awd)/len(awd)), max(awd)]])
        ovp['d'] = ad.reshape(len(doys), 24).T.tolist()
        ovp['wd'] = awd.reshape(len(doys), 24).T.tolist()
        ovp['days'] = array(doys, dtype = float)
        ovp['hours'] = arange(1, 25, dtype = float)        
        ovp['maxfreq'] = 100*numpy.max(dfreq)/len(awd)
        simnode['nbins'] = len(sbinvals)  
        simnode['d'] = array(cd).reshape(365, 24).T.tolist()
        simnode['wd'] = array(cwd).reshape(365, 24).T.tolist()        
        simnode['days'] = arange(1, 366, dtype = float)
        simnode['hours'] = arange(1, 25, dtype = float)               
        return {'FINISHED'}
    
class NODE_OT_SVF(bpy.types.Operator):
    bl_idname = "node.svf"
    bl_label = "Sky View Factor"
    bl_description = "Undertake a sky view factor study"
    bl_register = True
    bl_undo = False

    def invoke(self, context, event):
        scene = context.scene 
        svp = scene.vi_params
        svp.vi_display = 0
        dp = context.evaluated_depsgraph_get()
        
        if viparams(self, scene):            
            return {'CANCELLED'}
        
        shadobs = retobjs('livig')

        if not shadobs:
            self.report({'ERROR'},"No shading objects with a material attached.")
            return {'CANCELLED'}
            
        simnode = context.node
        svp['viparams']['restree'] = simnode.id_data.name
        clearscene(context, self)

        for o in scene.objects:
            o.vi_params.vi_type_string = ''
            
        calcobs = retobjs('ssc')

        if not calcobs:
            self.report({'ERROR'},"No objects have a light sensor material attached.")
            return {'CANCELLED'}

        svp['viparams']['visimcontext'] = 'SVF'

        if not svp.get('liparams'):
           svp['liparams'] = {}
           
        svp['liparams']['cp'], svp['liparams']['unit'], svp['liparams']['type'] = simnode.cpoint, 'SVF (%)', 'VI Shadow'
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
        calcsteps = len(frange) * sum(len([f for f in o.data.polygons if o.data.materials[f.material_index].vi_params.mattype == '1']) for o in calcobs)
        curres, reslists = 0, []
        pfile = progressfile(svp['viparams']['newdir'], datetime.datetime.now(), calcsteps)
        kivyrun = progressbar(os.path.join(svp['viparams']['newdir'], 'viprogress'), 'Sky View')

        for o in calcobs:
            ovp = o.vi_params
            for k in [k for k in ovp.keys()]:
                del ovp[k]
                
            if any([s < 0 for s in o.scale]):
                logentry('Negative scaling on calculation object {}. Results may not be as expected'.format(o.name))
                self.report({'WARNING'}, 'Negative scaling on calculation object {}. Results may not be as expected'.format(o.name))

            ovp['omin'], ovp['omax'], ovp['oave'] = {}, {}, {}
            bm = bmesh.new()
            bm.from_mesh(o.to_mesh())
            o.to_mesh_clear()
            clearlayers(bm, 'a')
            bm.transform(o.matrix_world)
            bm.normal_update()
            geom = bm.faces if simnode.cpoint == '0' else bm.verts
            geom.layers.int.new('cindex')
            cindex = geom.layers.int['cindex']
            [geom.layers.float.new('svf{}'.format(fi)) for fi in frange]
            avres, minres, maxres, g = [], [], [], 0
            
            if simnode.cpoint == '0':
                gpoints = [f for f in geom if o.data.materials[f.material_index].vi_params.mattype == '1']
            elif simnode.cpoint == '1':
                gpoints = [v for v in geom if any([o.data.materials[f.material_index].vi_params.mattype == '1' for f in v.link_faces])]

            for g, gp in enumerate(gpoints):
                gp[cindex] = g + 1

            for frame in frange: 
                g = 0               
                scene.frame_set(frame)
                shadtree = rettree(scene, shadobs, ('', '2')[simnode.signore])
                shadres = geom.layers.float['svf{}'.format(frame)]
                                    
                if gpoints:
                    posis = [gp.calc_center_median() + gp.normal.normalized() * simnode.offset for gp in gpoints] if simnode.cpoint == '0' else [gp.co + gp.normal.normalized() * simnode.offset for gp in gpoints]
                 
                    for chunk in chunks(gpoints, int(svp['viparams']['nproc']) * 200):
                        for gp in chunk:
                            pointres = array([(0, 1)[shadtree.ray_cast(posis[g], direc)[3] == None] for direc in valdirecs], dtype = int8)
                            gp[shadres] = (100*(numpy.sum(pointres)/lvaldirecs)).astype(int8)
                            g += 1

                        curres += len(chunk)

                        if pfile.check(curres) == 'CANCELLED':
                            return {'CANCELLED'}
              
                    shadres = [gp[shadres] for gp in gpoints]
                    ovp['omin']['svf{}'.format(frame)], ovp['omax']['svf{}'.format(frame)], ovp['oave']['svf{}'.format(frame)] = min(shadres), max(shadres), sum(shadres)/len(shadres)
                    reslists.append([str(frame), 'Zone', o.name, 'X', ' '.join(['{:.3f}'.format(p[0]) for p in posis])])
                    reslists.append([str(frame), 'Zone', o.name, 'Y', ' '.join(['{:.3f}'.format(p[1]) for p in posis])])
                    reslists.append([str(frame), 'Zone', o.name, 'Z', ' '.join(['{:.3f}'.format(p[2]) for p in posis])])
                    reslists.append([str(frame), 'Zone', o.name, 'SVF', ' '.join(['{:.3f}'.format(sr) for sr in shadres])])
                    avres.append(ovp['oave']['svf{}'.format(frame)])
                    minres.append(ovp['omin']['svf{}'.format(frame)])
                    maxres.append(ovp['omax']['svf{}'.format(frame)])

            reslists.append(['All', 'Frames', '', 'Frames', ' '.join(['{}'.format(f) for f in frange])])
            reslists.append(['All', 'Zone', o.name, 'Minimum', ' '.join(['{:.3f}'.format(mr) for mr in minres])])
            reslists.append(['All', 'Zone', o.name, 'Average', ' '.join(['{:.3f}'.format(mr) for mr in avres])])
            reslists.append(['All', 'Zone', o.name, 'Maximum', ' '.join(['{:.3f}'.format(mr) for mr in maxres])])
            bm.transform(o.matrix_world.inverted())
            bm.to_mesh(o.data)
            bm.free()
            o.vi_params.vi_type_string = 'LiVi Calc'

        svp.vi_leg_max, svp.vi_leg_min = 100, 0

        if kivyrun.poll() is None:
            kivyrun.kill()
        
        scene.frame_start, scene.frame_end = svp['liparams']['fs'], svp['liparams']['fe']
        svp['viparams']['vidisp'] = 'svf'
        simnode['reslists'] = reslists
#        simnode['year'] = 2018
        simnode['frames'] = [f for f in frange]
        simnode.postexport(scene)
        return {'FINISHED'} 
      
class NODE_OT_Shadow(bpy.types.Operator):
    bl_idname = "node.shad"
    bl_label = "Shadow Map"
    bl_description = "Undertake a shadow mapping"
    bl_register = True
    bl_undo = False

    def invoke(self, context, event):
        scene = context.scene 
        svp = scene.vi_params
        dp = bpy.context.evaluated_depsgraph_get()
        svp.vi_display = 0

        if viparams(self, scene):            
            return {'CANCELLED'}

        for o in scene.objects:
            o.vi_params.vi_type_string = ''
            
        shadobs = retobjs('livig')
        
        if not shadobs:
            self.report({'ERROR'},"No shading objects or none with a material attached.")
            return {'CANCELLED'}
            
        calcobs = retobjs('ssc')
        
        if not calcobs:
            self.report({'ERROR'},"No objects have a light sensor material attached.")
            return {'CANCELLED'}
        
        simnode = context.node
        svp['viparams']['restree'] = simnode.id_data.name
        clearscene(context, self)        
        svp['viparams']['visimcontext'] = 'Shadow'

        if not svp.get('liparams'):
           svp['liparams'] = {}

        svp['liparams']['cp'], svp['liparams']['unit'], svp['liparams']['type'] = simnode.cpoint, 'Sunlit time (%)', 'VI Shadow'
        simnode.preexport()
        (svp['liparams']['fs'], svp['liparams']['fe']) = (scene.frame_current, scene.frame_current) if simnode.animmenu == 'Static' else (simnode.startframe, simnode.endframe)
        cmap(svp)

        if simnode.starthour > simnode.endhour:
            self.report({'ERROR'},"End hour is before start hour.")
            return{'CANCELLED'}
        
        svp['viparams']['resnode'], simnode['Animation'] = simnode.name, simnode.animmenu
        (scmaxres, scminres, scavres) = [[x] * (svp['liparams']['fe'] - svp['liparams']['fs'] + 1) for x in (0, 100, 0)]
        nt = datetime.datetime.now()
        frange = range(svp['liparams']['fs'], svp['liparams']['fe'] + 1)
        time = datetime.datetime(2018, simnode.sdate.month, simnode.sdate.day, simnode.starthour - 1)
        y = 2018 if simnode.edoy >= simnode.sdoy else 2019
        endtime = datetime.datetime(y, simnode.edate.month, simnode.edate.day, simnode.endhour - 1)
        interval = datetime.timedelta(hours = 1/simnode.interval)        
        times = [time + interval*t for t in range(int((endtime - time)/interval) + simnode.interval) if simnode.starthour - 1 <= (time + interval*t).hour <= simnode.endhour  - 1]
        sps = array([solarPosition(t.timetuple().tm_yday, t.hour+t.minute/60, svp.latitude, svp.longitude)[2:] for t in times])
        valmask = array([sp[0] > 0 for sp in sps], dtype = int8)
        direcs = array([(-sin(sp[1]), -cos(sp[1]), tan(sp[0])) for sp in sps])  
        valdirecs = [mathutils.Vector((-sin(sp[1]), -cos(sp[1]), tan(sp[0]))) for sp in sps if sp[0] > 0]  
        lvaldirecs = len(valdirecs)
        ilvaldirecs = 1/lvaldirecs
        calcsteps = len(frange) * sum(len([f for f in o.data.polygons if o.data.materials[f.material_index].vi_params.mattype == '1']) for o in calcobs)
        curres, reslists = 0, []
        pfile = progressfile(svp['viparams']['newdir'], datetime.datetime.now(), calcsteps)
        kivyrun = progressbar(os.path.join(scene.vi_params['viparams']['newdir'], 'viprogress'), 'Shadow Map')
        logentry(f'Conducting shadow map calculation with {simnode.interval} samples per hour for {int(len(direcs)/simnode.interval)} total hours and {lvaldirecs} available sun hours')
        
        for oi, o in enumerate(calcobs):
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
            bm.from_mesh(o.to_mesh())
            o.to_mesh_clear()
            clearlayers(bm, 'a')
            bm.transform(o.matrix_world)
            bm.normal_update()
            geom = bm.faces if simnode.cpoint == '0' else bm.verts
            geom.layers.int.new('cindex')
            cindex = geom.layers.int['cindex']
            [geom.layers.float.new('sm{}'.format(fi)) for fi in frange]
            [geom.layers.float.new('hourres{}'.format(fi)) for fi in frange]
            avres, minres, maxres, g = [], [], [], 0
            
            if simnode.cpoint == '0':
                gpoints = [f for f in geom if o.data.materials[f.material_index].vi_params.mattype == '1']
            elif simnode.cpoint == '1':
                gpoints = [v for v in geom if any([o.data.materials[f.material_index].vi_params.mattype == '1' for f in v.link_faces])]

            for g, gp in enumerate(gpoints):
                gp[cindex] = g + 1
            
            for frame in frange: 
                g = 0               
                scene.frame_set(frame)
                shadtree = rettree(scene, shadobs, ('', '2')[simnode.signore])
                shadres = geom.layers.float['sm{}'.format(frame)]
                                  
                if gpoints:
                    posis = [gp.calc_center_median() + gp.normal.normalized() * simnode.offset for gp in gpoints] if simnode.cpoint == '0' else [gp.co + gp.normal.normalized() * simnode.offset for gp in gpoints]
                    allpoints = numpy.zeros((len(gpoints), len(direcs)), dtype=int8)
                    
                    for chunk in chunks(gpoints, int(svp['viparams']['nproc']) * 200):
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
                    ovp['ss{}'.format(frame)] = array(100 * ap).reshape(len(ovp['days']), len(ovp['hours'])).T.tolist()
                    ovp['omin']['sm{}'.format(frame)] = min(shadres)
                    ovp['omax']['sm{}'.format(frame)] = max(shadres)
                    ovp['oave']['sm{}'.format(frame)] = sum(shadres)/len(shadres)
                    reslists.append([str(frame), 'Zone', o.name, 'X', ' '.join(['{:.3f}'.format(p[0]) for p in posis])])
                    reslists.append([str(frame), 'Zone', o.name, 'Y', ' '.join(['{:.3f}'.format(p[1]) for p in posis])])
                    reslists.append([str(frame), 'Zone', o.name, 'Z', ' '.join(['{:.3f}'.format(p[2]) for p in posis])])
                    reslists.append([str(frame), 'Zone', o.name, 'Sunlit %', ' '.join(['{:.3f}'.format(sr) for sr in shadres])])
                    avres.append(ovp['oave']['sm{}'.format(frame)])
                    minres.append(ovp['omin']['sm{}'.format(frame)])
                    maxres.append(ovp['omax']['sm{}'.format(frame)])
            
            reslists.append(['All', 'Frames', '', 'Frames', ' '.join(['{}'.format(f) for f in frange])])
            reslists.append(['All', 'Zone', o.name, 'Minimum', ' '.join(['{:.3f}'.format(mr) for mr in minres])])
            reslists.append(['All', 'Zone', o.name, 'Average', ' '.join(['{:.3f}'.format(mr) for mr in avres])])
            reslists.append(['All', 'Zone', o.name, 'Maximum', ' '.join(['{:.3f}'.format(mr) for mr in maxres])])
            
            bm.transform(o.matrix_world.inverted())
            bm.to_mesh(o.data)
            bm.free()
#            o.vi_params.licalc = 1
            o.vi_params.vi_type_string = 'LiVi Calc'

        svp.vi_leg_max, svp.vi_leg_min = 100, 0

        if kivyrun.poll() is None:
            kivyrun.kill()

        scene.frame_start, scene.frame_end = svp['liparams']['fs'], svp['liparams']['fe']
        simnode['reslists'] = reslists
        simnode['frames'] = [f for f in frange]
#        simnode['year'] = 2015
        simnode.postexport(scene)
        svp['viparams']['vidisp'] = 'ss'
        return {'FINISHED'}

class NODE_OT_Li_Geo(bpy.types.Operator):
    bl_idname = "node.ligexport"
    bl_label = "LiVi geometry export"

    def execute(self, context):
        scene = context.scene
        svp = scene.vi_params
        svp.vi_display = 0

        if viparams(self, scene):
            return {'CANCELLED'}

        svp['viparams']['vidisp'] = ''
        svp['viparams']['viexpcontext'] = 'LiVi Geometry'
        objmode()
        clearfiles(svp['liparams']['objfilebase'])
        clearfiles(svp['liparams']['lightfilebase'])
        node = context.node
        node.preexport(scene)
        radgexport(self, node)
        node.postexport(scene)
        return {'FINISHED'}
    
class NODE_OT_Li_Con(bpy.types.Operator, ExportHelper):
    bl_idname = "node.liexport"
    bl_label = "LiVi context export"
    bl_description = "Export the scene to the Radiance file format"
    bl_register = True
    bl_undo = False

    def invoke(self, context, event):
        scene = context.scene
        self.svp = scene.vi_params
        self.svp.vi_display = 0
        
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

class NODE_OT_FileSelect(bpy.types.Operator, ImportHelper):
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
        if self.filepath.split(".")[-1] in self.fextlist:
            if self.nodeprop == 'epwname':
                self.node.epwname = self.filepath
            elif self.nodeprop == 'hdrname':
                self.node.hdrname = self.filepath
            elif self.nodeprop == 'skyname':
                self.node.skyname = self.filepath
            elif self.nodeprop == 'mtxname':
                self.node.mtxname = self.filepath
            elif self.nodeprop == 'resfilename':
                self.node.resfilename = self.filepath
            elif self.nodeprop == 'idffilename':
                self.node.idffilename = self.filepath
        if " " in self.filepath:
            self.report({'ERROR'}, "There is a space either in the filename or its directory location. Remove this space and retry opening the file.")
        return {'FINISHED'}

    def invoke(self,context,event):
        self.node = context.node
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}
    
class NODE_OT_HdrSelect(NODE_OT_FileSelect):
    bl_idname = "node.hdrselect"
    bl_label = "Select HDR/VEC file"
    bl_description = "Select the HDR sky image or vector file"
    filename_ext = ".HDR;.hdr;"
    filter_glob: bpy.props.StringProperty(default="*.HDR;*.hdr;", options={'HIDDEN'})
    nodeprop = 'hdrname'
    filepath:bpy.props.StringProperty(subtype='FILE_PATH', options={'HIDDEN', 'SKIP_SAVE'})
    fextlist = ("HDR", "hdr")
    
# BSDF Operators
class OBJECT_OT_Li_GBSDF(bpy.types.Operator):
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
                self.o.vi_params.bsdf_running = 0        
                return {'CANCELLED'}
            else:
                return{'PASS_THROUGH'}
        else:
            self.o.vi_params.bsdf_running = 0
            filepath = os.path.join(context.scene.vi_params['viparams']['newdir'], 'bsdfs', '{}.xml'.format(self.mat.name))
            
            if self.kivyrun.poll() is None:
                self.kivyrun.kill() 
            
            with open(filepath, 'r') as bsdffile:               
                self.mat.vi_params['bsdf']['xml'] = bsdffile.read()
                bsdf = parseString(self.mat.vi_params['bsdf']['xml'])
                self.mat.vi_params['bsdf']['direcs'] = [path.firstChild.data for path in bsdf.getElementsByTagName('WavelengthDataDirection')]
                self.mat.vi_params['bsdf']['type'] = [path.firstChild.data for path in bsdf.getElementsByTagName('AngleBasisName')][0]
                self.mat.vi_params['bsdf']['filepath'] = filepath
                self.mat.vi_params['bsdf']['proxied'] = self.o.vi_params.li_bsdf_proxy

            context.scene.vi_params['viparams']['vidisp'] = 'bsdf'
            return {'FINISHED'}
    
    def execute(self, context):
        scene = context.scene
        svp = scene.vi_params
        dp = bpy.context.evaluated_depsgraph_get()
        vl = context.view_layer
        self.o = context.object
        ovp = self.o.vi_params
        vl.objects.active = None
        
        if viparams(self, scene):
            return {'CANCELLED'}
        bsdfmats = [mat for mat in self.o.data.materials if mat.vi_params.radmatmenu == '8']
        
        if bsdfmats:
            self.mat = bsdfmats[0]
            mvp = self.mat.vi_params
            mvp['bsdf'] = {} 
        else:
            self.report({'ERROR'}, '{} does not have a BSDF material attached'.format(self.o.name))

        bm = bmesh.new()    
        bm.from_object(self.o, dp)
        bm.transform(self.o.matrix_world)
        bm.normal_update()
        bsdffaces = [face for face in bm.faces if self.o.material_slots[face.material_index].material.vi_params.radmatmenu == '8']    
        
        if bsdffaces:
            fvec = bsdffaces[0].normal
    #        mvp['bsdf']['up'] = '{0[0]:.4f} {0[1]:.4f} {0[2]:.4f}'.format(fvec)
        else:
            self.report({'ERROR'}, '{} does not have a BSDF material associated with any faces'.format(self.o.name))
            return
        
        self.pfile = progressfile(svp['viparams']['newdir'], datetime.datetime.now(), 100)
        self.kivyrun = progressbar(os.path.join(svp['viparams']['newdir'], 'viprogress'), 'BSDF')
        zvec, xvec, yvec = mvp.li_bsdf_up, mathutils.Vector((1, 0, 0)), mathutils.Vector((0, 1, 0))
        svec = mathutils.Vector.cross(fvec, zvec)
        bsdfrotz = mathutils.Matrix.Rotation(mathutils.Vector.angle(fvec, zvec), 4, svec)
        bm.transform(bsdfrotz)
        bm.normal_update()
        zvec.rotate(bsdfrotz.to_euler())

        try:
            bsdfrotz = mathutils.Matrix.Rotation(mathutils.Vector.angle(zvec, yvec), 4, mathutils.Vector((0, 0, -1)))
            bm.transform(bsdfrotz)
            bm.normal_update()
            zvec.rotate(bsdfrotz.to_euler())
        except Exception as e:
            logentry('Transform unneccesary')

#        mvp['bsdf']['up'] = '{0[0]:.4f} {0[1]:.4f} {0[2]:.4f}'.format(mvp.li_bsdf_up)
        vposis = list(zip(*[v.co[:] for v in bm.verts]))
        (maxx, maxy, maxz) = [max(p) for p in vposis]
        (minx, miny, minz) = [min(p) for p in vposis]
        bsdftrans = mathutils.Matrix.Translation(mathutils.Vector((-(maxx + minx)/2, -(maxy + miny)/2, -maxz)))
        bm.transform(bsdftrans)
        mradfile = ''.join([m.vi_params.radmat(scene) for m in self.o.data.materials if m.vi_params.radmatmenu != '8'])                  
        gradfile = radpoints(self.o, [face for face in bm.faces if self.o.material_slots and face.material_index < len(self.o.material_slots) and self.o.material_slots[face.material_index].material.vi_params.radmatmenu != '8'], 0)
        bm.free()  
        bsdfsamp = ovp.li_bsdf_ksamp if ovp.li_bsdf_tensor == ' ' else 2**(int(ovp.li_bsdf_res) * 2) * int(ovp.li_bsdf_tsamp)    
#        gbcmd = "genBSDF -geom {} -r '{}' {} {} -c {} {} -n {}".format(ovp.li_bsdf_dimen,  ovp.li_bsdf_rcparam,  ovp.li_bsdf_tensor, (ovp.li_bsdf_res, ' ')[ovp.li_bsdf_tensor == ' '], bsdfsamp, ovp.li_bsdf_direc, svp['viparams']['nproc'])
        # Adding MGF geometry does not work (black inner face)
        gbcmd = "genBSDF +geom {} -r '{}' {} {} -c {} {} -n {}".format(ovp.li_bsdf_dimen,  ovp.li_bsdf_rcparam,  ovp.li_bsdf_tensor, (ovp.li_bsdf_res, ' ')[ovp.li_bsdf_tensor == ' '], bsdfsamp, ovp.li_bsdf_direc, svp['viparams']['nproc'])
        logentry('genBSDF running with command: {}'.format(gbcmd))
        
        with open(os.path.join(svp['viparams']['newdir'], 'bsdfs', '{}_mg'.format(self.mat.name)), 'w') as mgfile:
            mgfile.write(mradfile+gradfile)

        with open(os.path.join(svp['viparams']['newdir'], 'bsdfs', '{}_mg'.format(self.mat.name)), 'r') as mgfile: 
            with open(os.path.join(svp['viparams']['newdir'], 'bsdfs', '{}.xml'.format(self.mat.name)), 'w') as bsdffile:
                self.bsdfrun = Popen(shlex.split(gbcmd), stdin = mgfile, stdout = bsdffile)

        mvp['bsdf']['type'] = 'LBNL/Klems Full' if  ovp.li_bsdf_tensor == ' ' else 'Tensor'      
        vl.objects.active = self.o
        wm = context.window_manager
        self._timer = wm.event_timer_add(1, window = context.window)
        wm.modal_handler_add(self) 
        ovp.bsdf_running = 1       
        return {'RUNNING_MODAL'}
        
class MATERIAL_OT_Li_LBSDF(bpy.types.Operator, ImportHelper):
    bl_idname = "material.load_bsdf"
    bl_label = "Select BSDF file"
    filename_ext = ".XML;.xml;"
    filter_glob: bpy.props.StringProperty(default="*.XML;*.xml;", options={'HIDDEN'})
    filepath: bpy.props.StringProperty(subtype='FILE_PATH', options={'HIDDEN', 'SKIP_SAVE'})
    
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
                context.material['bsdf']['filepath'] = self.filepath
            return {'FINISHED'}

    def invoke(self,context,event):
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}
    
class MATERIAL_OT_Li_DBSDF(bpy.types.Operator):
    bl_idname = "material.del_bsdf"
    bl_label = "Del BSDF"
    bl_description = "Delete a BSDF for the current selected object"
    bl_register = True
    bl_undo = False
    
    def execute(self, context):
        del context.material['bsdf']
        return {'FINISHED'}
        
class MATERIAL_OT_Li_SBSDF(bpy.types.Operator, ExportHelper):
    bl_idname = "material.save_bsdf"
    bl_label = "Save BSDF"
    bl_description = "Save a BSDF for the current selected object"
    bl_register = True
    filename = "material"
    filename_ext = ".xml"
    filter_glob: bpy.props.StringProperty(default="*.XML;*.xml;", options={'HIDDEN'})
    filepath: bpy.props.StringProperty(subtype='FILE_PATH', options={'HIDDEN', 'SKIP_SAVE'})
    
    def draw(self,context):
        layout = self.layout
        row = layout.row()
        row.label(text="Save BSDF XML file with the file browser", icon='WORLD_DATA')

    def execute(self, context):        
        with open(self.filepath, 'w') as bsdfsave:
            bsdfsave.write(context.material.vi_params['bsdf']['xml'])
        return {'FINISHED'}

    def invoke(self,context,event):
        self.filepath= '{}.xml'.format(context.material.name)
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}
        
class NODE_OT_Li_Pre(bpy.types.Operator, ExportHelper):
    bl_idname = "node.radpreview"
    bl_label = "LiVi preview"
    bl_description = "Prevew the scene with Radiance"
    bl_register = True
    bl_undo = False
    
    def modal(self, context, event):
        if event.type == 'TIMER':
            if self.rvurun.poll() is not None: # If finished
                self.simnode.run = 0
                
                for line in self.rvurun.stderr:
                    if  b'fatal IO error' not in line and b'events remaining' not in line and b'Broken pipe' not in line and b'explicit kill' not in line:
                        logentry(line)
                    for rvuerr in rvuerrdict:
                        if rvuerr in line.decode():
                            self.report({'ERROR'}, rvuerrdict[rvuerr])
                            return {'CANCELLED'}
                
                return {'FINISHED'}
            else:           
                return {'PASS_THROUGH'}
        else:
            return {'PASS_THROUGH'}
        

    def invoke(self, context, event):
        scene = context.scene
        svp = scene.vi_params
        svp.vi_display = 0

        if viparams(self, scene):
            return {'CANCELLED'}
        
        objmode()
        self.simnode = context.node
        cam = bpy.data.objects[self.simnode.camera.lstrip()]

        if not cam:
            self.report({'ERROR'}, "There is no camera in the scene. Radiance preview will not work")
            return {'CANCELLED'}
        else:
            frame = scene.frame_current
            self.simnode.presim()
            self.pmfile = os.path.join(svp['viparams']['newdir'], 'pmprogress')
            svp['liparams']['fs'] = min([c['fs'] for c in (self.simnode['goptions'], self.simnode['coptions'])])
            svp['liparams']['fe'] = max([c['fe'] for c in (self.simnode['goptions'], self.simnode['coptions'])])

            if frame not in range(svp['liparams']['fs'], svp['liparams']['fe'] + 1):
                frame = svp['liparams']['fe'] if frame > svp['liparams']['fe'] else frame
                frame = svp['liparams']['fs'] if frame < svp['liparams']['fs'] else frame
                scene.frame_set(frame)
                self.report({'WARNING'}, "Current frame is not within the exported frame range and has been adjusted")

            curres = 0.1
            createradfile(scene, frame, self, self.simnode)
            createoconv(scene, frame, self, self.simnode)
            cang = '180 -vth ' if self.simnode['coptions']['Context'] == 'Basic' and self.simnode['coptions']['Type'] == '1' else cam.data.angle_x*180/pi
            vv = 180 if self.simnode['coptions']['Context'] == 'Basic' and self.simnode['coptions']['Type'] == '1' else cam.data.angle_y*180/pi
            vd = (0.001, 0, -1*cam.matrix_world[2][2]) if (round(-1*cam.matrix_world[0][2], 3), round(-1*cam.matrix_world[1][2], 3)) == (0.0, 0.0) else [-1*cam.matrix_world[i][2] for i in range(3)]
            
            if self.simnode.pmap:
                self.pfile = progressfile(svp['viparams']['newdir'], datetime.datetime.now(), 100)
                self.kivyrun = progressbar(os.path.join(svp['viparams']['newdir'], 'viprogress'), 'Photon Map')
                amentry, pportentry, cpentry, cpfileentry = retpmap(self.simnode, frame, scene)
                open('{}-{}'.format(self.pmfile, frame), 'w')
                pmcmd = 'mkpmap {8} -t 2 -e "{1}" {6} -fo+ -bv{9} -apD 0.1 {0} -apg "{7}-{2}.gpm" {3} {4} {5} "{7}-{2}.oct"'.format(pportentry, '{}-{}'.format(self.pmfile, frame), frame, 
                        self.simnode.pmapgno, cpentry, amentry, ('-n {}'.format(svp['viparams']['wnproc']), '')[sys.platform == 'win32'], svp['viparams']['filebase'], self.simnode.pmapoptions, ('-', '+')[self.simnode.bfv])
                logentry('Photon map command: {}'.format(pmcmd))
                os.chdir(svp['viparams']['newdir'])
                pmrun = Popen(shlex.split(pmcmd), stderr = PIPE, stdout = PIPE)

                for line in pmrun.stderr:
                    logentry('Photon mapping error: {}'.format(line.decode()))
                
                while pmrun.poll() is None:   
                    sleep(5)
                    with open('{}-{}'.format(self.pmfile, frame), 'r') as vip:
                        for line in vip.readlines()[::-1]:
                            if '%' in line:
                                for entry in line.split():
                                    if '%' in entry:
                                        curres = float(entry[:-2])
                                        break
                            break
                                
                    if self.pfile.check(curres) == 'CANCELLED': 
                        pmrun.kill()    
                        return {'CANCELLED'}
                
                if self.kivyrun.poll() is None:
                    self.kivyrun.kill()
                        
                with open('{}-{}'.format(self.pmfile, frame), 'r') as pmapfile:
                    for line in pmapfile.readlines():
                        if line in pmerrdict:
                            logentry(line)
                            self.report({'ERROR'}, pmerrdict[line])
                            return {'CANCELLED'}
 
                if self.simnode.pmappreview:
                    create_empty_coll(context, 'LiVi Results')
                    gpmbm = bmesh.new()

                    for l, line in enumerate(Popen(shlex.split('pmapdump -a -c 0 0 1 {0}-{1}.gpm'.format(svp['viparams']['filebase'], frame)), stdout = PIPE, stderr = PIPE).stdout):
                        dl = line.decode().split()
                        matrix = Matrix.Translation(Vector([float(x) for x in dl[:3]]))
                        bmesh.ops.create_icosphere(gpmbm, subdivisions = 2, diameter = 0.05, matrix = matrix, calc_uvs = False)
                        
                        if l > self.simnode.pmapvno:
                            break

                    gpmmesh = bpy.data.meshes.new("GPM_Mesh")
                    gpmbm.to_mesh(gpmmesh)
                    gpmbm.free()
                    gpmobj = bpy.data.objects.new("GlobalPM", gpmmesh)
                    gpmobj.vi_params.vi_type_string = 'LiVi Res'
                    scene.collection.objects.link(gpmobj)
                    move_to_coll(bpy.context, 'LiVi Results', gpmobj)

                    if cpentry:
                        cpmbm = bmesh.new()

                        for l, line in enumerate(Popen(shlex.split('pmapdump -a -c 0 0 1 {0}-{1}.cpm'.format(svp['viparams']['filebase'], frame)), stdout = PIPE, stderr = PIPE).stdout):
                            dl = line.decode().split()
                            matrix = Matrix.Translation(Vector([float(x) for x in dl[:3]]))
                            bmesh.ops.create_icosphere(cpmbm, subdivisions = 2, diameter = 0.05, matrix = matrix, calc_uvs = False)
                            
                            if l > self.simnode.pmapvno:
                                break

                        cpmmesh = bpy.data.meshes.new("CPM_Mesh")
                        cpmbm.to_mesh(cpmmesh)
                        cpmbm.free()
                        cpmobj = bpy.data.objects.new("CausticPM", cpmmesh)
                        cpmobj.vi_params.vi_type_string = 'LiVi Res'
                        scene.collection.objects.link(cpmobj)
                        move_to_coll(bpy.context, 'LiVi Results', cpmobj)

                    # with open("{0}-{1}pmd.oct".format(svp['viparams']['filebase'], frame), 'wb') as octfile:
                    #     occmd = 'oconv -i "{0}-{1}.oct" "!pmapdump {0}-{1}.gpm{2}"'.format(svp['viparams']['filebase'], frame, (' {}-{}.cpm'.format(svp['viparams']['filebase'], frame), '')[not self.simnode.pmapcno])
                    #     logentry('Running pmapdump: {}'.format(occmd))
                    #     ocrun = Popen(shlex.split(occmd), stdout = octfile, stderr = PIPE)
                    #     ocrun.wait()
                    # for line in ocrun.stderr:
                    #     logentry('Pmap preview error: {}'.format(line.decode()))

                    # rvucmd = 'rvu -w {9} -n {0} -vv {1:.3f} -vh {2:.3f} -vd {3[0]:.3f} {3[1]:.3f} {3[2]:.3f} -vp {4[0]:.3f} {4[1]:.3f} {4[2]:.3f} -vu {8[0]:.3f} {8[1]:.3f} {8[2]:.3f} {5} "{6}-{7}pmd.oct"'.format(svp['viparams']['wnproc'], 
                    #              vv, cang, vd, cam.location, self.simnode['rvuparams'], svp['viparams']['filebase'], scene.frame_current, cam.matrix_world.to_quaternion()@mathutils.Vector((0, 1, 0)), ('', '-i')[self.simnode.illu])

#                else:
                rvucmd = 'rvu -w {11} -ap "{8}" 50 {9} -n {0} -vv {1:.3f} -vh {2:.3f} -vd {3[0]:.3f} {3[1]:.3f} {3[2]:.3f} -vp {4[0]:.3f} {4[1]:.3f} {4[2]:.3f} -vu {10[0]:.3f} {10[1]:.3f} {10[2]:.3f} {5} "{6}-{7}.oct"'.format(svp['viparams']['wnproc'], 
                                 vv, cang, vd, cam.location, self.simnode['rvuparams'], svp['viparams']['filebase'], scene.frame_current, '{}-{}.gpm'.format(svp['viparams']['filebase'], frame), cpfileentry, cam.matrix_world.to_quaternion()@mathutils.Vector((0, 1, 0)), ('', '-i')[self.simnode.illu])

            else:
                rvucmd = 'rvu -w {9} -n {0} -vv {1:.3f} -vh {2:.3f} -vd {3[0]:.3f} {3[1]:.3f} {3[2]:.3f} -vp {4[0]:.3f} {4[1]:.3f} {4[2]:.3f} -vu {8[0]:.3f} {8[1]:.3f} {8[2]:.3f} {5} "{6}-{7}.oct"'.format(svp['viparams']['wnproc'], 
                                 vv, cang, vd, cam.location, self.simnode['rvuparams'], svp['viparams']['filebase'], scene.frame_current, cam.matrix_world.to_quaternion()@ mathutils.Vector((0, 1, 0)), ('', '-i')[self.simnode.illu])

            logentry('Rvu command: {}'.format(rvucmd))
            self.rvurun = Popen(shlex.split(rvucmd), stdout = PIPE, stderr = PIPE)
            context.node.run = 1
            wm = context.window_manager
            self._timer = wm.event_timer_add(1, window = context.window)
            wm.modal_handler_add(self)
            self.simnode.hide = 0
            return {'RUNNING_MODAL'}

class NODE_OT_Li_Sim(bpy.types.Operator):
    bl_idname = "node.livicalc"
    bl_label = "LiVi simulation"
    bl_register = True
    bl_undo = False

    def modal(self, context, event):
#        if self.kivyrun.poll() is not None:
        self.simnode.postsim(self.reslists)
        self.report({'INFO'}, "Simulation is finished")
        return {'FINISHED'}

    def invoke(self, context, event):
        scene = context.scene
        vl = context.view_layer
        frame = scene.frame_current
        svp = scene.vi_params
        svp.vi_display = 0

        if viparams(self, scene):
            return {'CANCELLED'}
                    
        objmode()
        clearscene(context, self)
        self.simnode = context.node
        self.simnode.presim()
        contextdict = {'Basic': 'LiVi Basic', 'CBDM': 'LiVi CBDM'}        
        
        # Set scene parameters
        svp['viparams']['visimcontext'] = contextdict[self.simnode['coptions']['Context']]
        svp['liparams']['fs'] = min((self.simnode['coptions']['fs'], self.simnode['goptions']['fs'])) 
        svp['liparams']['fe'] = max((self.simnode['coptions']['fe'], self.simnode['goptions']['fe'])) 
        svp['liparams']['cp'] = self.simnode['goptions']['cp']
        svp['liparams']['unit'] = self.simnode['coptions']['unit']
        svp['liparams']['type'] =self. simnode['coptions']['Type']
        scene.frame_start, scene.frame_end = svp['liparams']['fs'], svp['liparams']['fe']       
        self.simnode.sim(scene)
        pfs, epfs, curres = [], [], 0
        rtcmds, rccmds = [], []
        scontext = self.simnode['coptions']['Context']
        subcontext = self.simnode['coptions']['Type']
        patches = self.simnode['coptions']['cbdm_res']
        svp['liparams']['maxres'], svp['liparams']['minres'], svp['liparams']['avres'] = {}, {}, {}

        if frame not in range(svp['liparams']['fs'], svp['liparams']['fe'] + 1):
            frame = svp['liparams']['fe'] if frame > svp['liparams']['fe'] else frame
            frame = svp['liparams']['fs'] if frame < svp['liparams']['fs'] else frame
            scene.frame_set(frame)
            self.report({'WARNING'}, "Current frame is not within the exported frame range and has been adjusted")
        
        frames = range(svp['liparams']['fs'], svp['liparams']['fe'] + 1)
        
        for f, frame in enumerate(frames):
            if createradfile(scene, frame, self, self.simnode) == 'CANCELLED' or createoconv(scene, frame, self, self.simnode) == 'CANCELLED':
                return {'CANCELLED'}
        
            if self.simnode.pmap:
                pmappfile = open(os.path.join(svp['viparams']['newdir'], 'viprogress'), 'w')
                pmappfile.close()
                pfile = progressfile(svp['viparams']['newdir'], datetime.datetime.now(), 100)
                self.kivyrun = progressbar(os.path.join(svp['viparams']['newdir'], 'viprogress'), 'Photon map')
                errdict = {'fatal - too many prepasses, no global photons stored\n': "Too many prepasses have occurred. Make sure light sources can see your geometry",
                            'fatal - too many prepasses, no global photons stored, no caustic photons stored\n': "Too many prepasses have occurred. Turn off caustic photons and encompass the scene",
                            'fatal - zero flux from light sources\n': "No light flux, make sure there is a light source and that photon port normals point inwards",
                            'fatal - no light sources\n': "No light sources. Photon mapping does not work with HDR skies",
                            'fatal - no valid photon ports found\n': 'Re-export the geometry'}
                amentry, pportentry, cpentry, cpfileentry = retpmap(self.simnode, frame, scene)
                open('{}.pmapmon'.format(svp['viparams']['filebase']), 'w')

                if scontext == 'Basic' or (scontext == 'CBDM' and subcontext == '0'):
                    pmcmd = 'mkpmap -n {6} -t 10 -e "{1}.pmapmon" -fo+ -bv+ -apD 0.001 {0} -apg "{1}-{2}.gpm" {3} {4} {5} "{1}-{2}.oct"'.format(pportentry, svp['viparams']['filebase'], frame, self.simnode.pmapgno, cpentry, amentry, svp['viparams']['wnproc'])
                else:
                    pmcmd = 'mkpmap -n {3} -t 10 -e "{1}.pmapmon" -fo+ -bv+ -apC "{1}.cpm" {0} "{1}-{2}.oct"'.format(self.simnode.pmapgno, svp['viparams']['filebase'], frame, svp['viparams']['wnproc'])
            
                logentry('Generating photon map: {}'.format(pmcmd))
                pmrun = Popen(shlex.split(pmcmd), stderr = PIPE, stdout = PIPE)
                
                while pmrun.poll() is None:   
                    sleep(10)
                    with open('{}.pmapmon'.format(svp['viparams']['filebase']), 'r') as vip:
                        for line in vip.readlines()[::-1]:
                            if '%' in line:
                                curres = float(line.split()[6][:-2])/len(frames)
                                break
                                    
                    if pfile.check(curres) == 'CANCELLED': 
                        pmrun.kill()                                   
                        return {'CANCELLED'}
                
                if self.kivyrun.poll() is None:
                    self.kivyrun.kill()
                        
                with open('{}.pmapmon'.format(svp['viparams']['filebase']), 'r') as pmapfile:
                    pmlines = pmapfile.readlines()
                    if pmlines:
                        for line in pmlines:
                            if line in errdict:
                                self.report({'ERROR'}, errdict[line])
                                return {'CANCELLED'}
                            if 'fatal - ' in line:
                                self.report({'ERROR'}, line)
                                return {'CANCELLED'}
                    else:
                        self.report({'ERROR'}, 'There is a problem with pmap generation. Check there are no non-ascii characters in the project directory file path')
                        return {'CANCELLED'}
                
            if scontext == 'Basic' or (scontext == 'CBDM' and subcontext == '0'):# or (context == 'Compliance' and int(subcontext) < 3):
                if os.path.isfile("{}-{}.af".format(svp['viparams']['filebase'], frame)):
                    os.remove("{}-{}.af".format(svp['viparams']['filebase'], frame))
                if self.simnode.pmap:
                    rtcmds.append('rtrace -n {0} -w {1} -ap "{2}-{3}.gpm" 50 {4} -faa -h -ov -I "{2}-{3}.oct"'.format(svp['viparams']['nproc'], self.simnode['radparams'], svp['viparams']['filebase'], frame, cpfileentry)) #+" | tee "+lexport.newdir+lexport.fold+self.simlistn[int(lexport.metric)]+"-"+str(frame)+".res"
                else:
                    rtcmds.append('rtrace -n {0} -w {1} -faa -h -ov -I "{2}-{3}.oct"'.format(svp['viparams']['nproc'], self.simnode['radparams'], svp['viparams']['filebase'], frame)) #+" | tee "+lexport.newdir+lexport.fold+self.simlistn[int(lexport.metric)]+"-"+str(frame)+".res"
            else:
                if self.simnode.pmap:
                    rccmds.append('rcontrib -w  -h -I -fo -ap {2}.cpm -bn {4} {0} -n {1} -f tregenza.cal -b tbin -m sky_glow "{2}-{3}.oct"'.format(self.simnode['radparams'], svp['viparams']['nproc'], svp['viparams']['filebase'], frame, patches))
                else:   
                    rccmds.append('rcontrib -w  -h -I -fo -bn {} {} -n {} -f tregenza.cal -b tbin -m sky_glow "{}-{}.oct"'.format(patches, self.simnode['radparams'], svp['viparams']['nproc'], svp['viparams']['filebase'], frame))

        try:
            tpoints = [o.vi_params['rtpnum'] for o in bpy.data.objects if o.name in svp['liparams']['livic']]
        except:
            self.report({'ERROR'}, 'Re-export the LiVi geometry')
            return {'CANCELLED'}

        calcsteps = sum(tpoints) * len(frames)
        pfile = progressfile(svp['viparams']['newdir'], datetime.datetime.now(), calcsteps)
        self.kivyrun = progressbar(os.path.join(svp['viparams']['newdir'], 'viprogress'), 'Lighting')
        self.reslists = []
        obs = [o for o in bpy.data.objects if o.name in svp['liparams']['livic']]

        for oi, o in enumerate(obs):
            ovp = o.vi_params
            curres = sum(tpoints[:oi] * len(frames))
            selobj(vl, o)
            ovp['omax'], ovp['omin'], ovp['oave']  = {}, {}, {}
            
            if scontext == 'Basic':
                bccout = ovp.basiccalcapply(scene, frames, rtcmds, self.simnode, curres, pfile)
                if bccout == 'CANCELLED':
                    if self.kivyrun.poll() is None:
                       self.kivyrun.kill()
                    return {'CANCELLED'}
                else:
                    self.reslists += bccout
                    
            elif scontext == 'CBDM' and subcontext == '0':
                lhout = ovp.lhcalcapply(scene, frames, rtcmds, self.simnode, curres, pfile)
                if lhout  == 'CANCELLED':
                    if self.kivyrun.poll() is None:
                        self.kivyrun.kill()
                    return {'CANCELLED'}
                else:
                    self.reslists += lhout
            
            elif (scontext == 'CBDM' and subcontext in ('1', '2')):# or (context == 'Compliance' and subcontext == '3'):
                cbdmout = ovp.udidacalcapply(scene, frames, rccmds, self.simnode, curres, pfile)
                if cbdmout == 'CANCELLED':
                    if self.kivyrun.poll() is None:
                        self.kivyrun.kill()
                    return {'CANCELLED'}
                else:
                    self.reslists += cbdmout
            
        if self.kivyrun.poll() is None:
            self.kivyrun.kill()
#        calcout = li_calc(self, simnode, livisimacc(simnode))
        # context.evaluated_depsgraph_get()
        # simnode = context.node

        # if calcout == 'CANCELLED':
        #     self.report({'ERROR'},"Simulation was cancelled. See log file")
        #     return {'CANCELLED'}
        # else:
        #     try:
        #         simnode['reslists'] = calcout
        #     except:
        #         self.report({'ERROR'}, "Previous instance of rvu still running?")
        #         logentry('ID problem')
        #         return {'CANCELLED'}

        svp['viparams']['vidisp'] = 'li'
        svp['viparams']['resnode'] = self.simnode.name
        svp['viparams']['restree'] = self.simnode.id_data.name   
        context.window_manager.modal_handler_add(self)    
        return {'RUNNING_MODAL'}
    
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
                self.pmruns.append(Popen(shlex.split(self.pmcmds[self.p]), stderr = PIPE))
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
                    self.rpruns.append(Popen(shlex.split(self.rpiececmds[self.frame - self.fs]), stdin = echo.stdout, stderr = PIPE))

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
                while sum([rp.poll() is None for rp in self.rpruns]) < self.np and len(self.rpruns) < self.frames:
                    f = self.fs + len(self.rpruns)
                    
                    with open("{}-{}.hdr".format(os.path.join(self.folder, 'images', self.basename), f), 'w') as imfile:
                        self.rpruns.append(Popen(shlex.split(self.rpictcmds[len(self.rpruns)]), stdout=imfile, stderr = PIPE))

                try:
                    if [rp.poll() for rp in self.rpruns][self.frame - self.fs] is not None:
                        self.images.append(os.path.join(self.folder, 'images', '{}-{}.hdr'.format(self.basename, self.frame)))
                        if self.frame < self.fe:
                            self.frame += 1
                except:
                    print('Frame passing')
                                        
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
                            newpercent = (float(lineentry.strip('%')) * sum([r.poll() is None for r in self.rpruns]) + 100 * sum([r.poll() is not None for r in self.rpruns]))/self.frames
                            
                            if self.percent != newpercent:
                                self.percent = newpercent
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
        
        if simnode.camera and bpy.data.cameras.get(simnode.camera.lstrip()):
            self.percent = 0
            self.reslists, self.images = [], []
            self.res = []
            self.rpictfile = os.path.join(svp['viparams']['newdir'], 'rpictprogress')
            self.pmfile = os.path.join(svp['viparams']['newdir'], 'pmprogress')
            simnode.presim()
            svp['liparams']['fs'], svp['liparams']['fe'] =  simnode.retframes()
            self.frames = self.fe - self.fs + 1
            self.frame = self.fs
            self.frameold = self.frame
            self.rpruns, self.pmruns = [], []   
            self.processors = simnode.processors
            self.processes = simnode.processes
            self.radparams = simnode['rpictparams']
            self.viewparams = simnode['viewparams']
            self.pmparams = simnode['pmparams']
            self.pmaps = simnode['pmaps']
            self.pmapgnos = simnode['pmapgnos']
            self.pmapcnos = simnode['pmapcnos']
            self.folder = svp['viparams']['newdir']
            self.fb = svp['viparams']['filebase']
            self.np = int(svp['viparams']['nproc'])
            self.mp = 0 if sys.platform == 'win32' else simnode.mp 
            self.basename = simnode['basename']
            
            for frame in range(self.fs, self.fe + 1):
                createradfile(scene, frame, self, simnode)
                createoconv(scene, frame, self, simnode)
                with open('{}-{}'.format(self.pmfile, frame), 'w'):
                    pass
                
            scene.frame_set(svp['liparams']['fs'])
            self.pmcmds = ['mkpmap {7} -t 10 -e "{6}" -bv+ +fo -apD 0.001 {0} -apg "{1}-{2}.gpm" {3} {4} {5} "{1}-{2}.oct"'.format(self.pmparams[str(frame)]['pportentry'], 
                            svp['viparams']['filebase'], frame, self.pmapgnos[str(frame)], self.pmparams[str(frame)]['cpentry'], self.pmparams[str(frame)]['amentry'], 
                            '{}-{}'.format(self.pmfile, frame),  ('-n {}'.format(svp['viparams']['wnproc']), '')[sys.platform == 'win32']) for frame in range(self.fs, self.fe + 1)]                   

            self.rppmcmds = [('', ' -ap "{}" {}'.format('{}-{}.gpm'.format(svp['viparams']['filebase'], frame), self.pmparams[str(frame)]['cpfileentry']))[self.pmaps[frame - self.fs]] for frame in range(self.fs, self.fe + 1)]
            self.rpictcmds = ['rpict -t 10 -e "{}" '.format(self.rpictfile) + ' '.join(['{0[0]} {0[1]}'.format(i) for i in self.viewparams[str(frame)].items()]) + self.rppmcmds[frame - self.fs] + self.radparams + '"{0}-{1}.oct"'.format(svp['viparams']['filebase'], frame) for frame in range(self.fs, self.fe + 1)]
            self.rpiececmds = ['rpiece -t 10 -e "{}" '.format(self.rpictfile) + ' '.join(['{0[0]} {0[1]}'.format(i) for i in self.viewparams[str(frame)].items()]) + self.rppmcmds[frame - self.fs] + self.radparams + '-o "{2}-{1}.hdr" "{0}-{1}.oct"'.format(svp['viparams']['filebase'], frame, os.path.join(svp['viparams']['newdir'], 'images', self.basename)) for frame in range(self.fs, self.fe + 1)]
            self.starttime = datetime.datetime.now()
            self.pfile = progressfile(self.folder, datetime.datetime.now(), 100)
            (self.pmfin, flag) = (0, 'Photon Maps') if sum(self.pmaps) else (1, 'Radiance Images')
            self.kivyrun = progressbar(os.path.join(self.folder, 'viprogress'), flag)

            if os.path.isfile("{}-{}.hdr".format(os.path.join(svp['viparams']['newdir'], 'images', self.basename), self.frame)):
               os.remove("{}-{}.hdr".format(os.path.join(svp['viparams']['newdir'], 'images', self.basename), self.frame))
                
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
        imnode = glnode.inputs['Image'].links[0].from_node if glnode.inputs['Image'].links else glnode
        # imnode = glnode.inputs[0].links[0].from_node
        glnode.presim()
        
        for i, im in enumerate(imnode['images']):
            glfile = os.path.join(svp['viparams']['newdir'], 'images', '{}-{}.hdr'.format(glnode['hdrname'], i + svp['liparams']['fs']))
            egcmd = 'evalglare {} -c {}'.format(('-u {0[0]} {0[1]} {0[2]}'.format(glnode.gc), '')[glnode.rand], glfile)
            print(im)
            with open(im, 'r') as hdrfile:
                egrun = Popen(egcmd.split(), stdin = hdrfile, stdout = PIPE, stderr = PIPE)
            
            if imnode != glnode:
                time = datetime.datetime(2019, 1, 1, imnode['coptions']['shour'], 0) + datetime.timedelta(imnode['coptions']['sdoy'] - 1) if imnode['coptions']['anim'] == '0' else \
                    datetime.datetime(2019, 1, 1, int(imnode['coptions']['shour']), int(60*(imnode['coptions']['shour'] - int(imnode['coptions']['shour'])))) + datetime.timedelta(imnode['coptions']['sdoy'] - 1) + datetime.timedelta(hours = int(imnode['coptions']['interval']*i), 
                                    seconds = int(60*(imnode['coptions']['interval']*i - int(imnode['coptions']['interval']*i))))
            else:
                time = datetime.datetime(2019, 1, 1, 1)
            
            with open(os.path.join(svp['viparams']['newdir'], 'images', "temp.glare"), "w") as glaretf:
                for line in egrun.stderr:
                    logentry("Evalglare message: {}".format(line.decode()))
                    if 'perspective' in line.decode():
                        self.report({'ERROR'}, 'Images are not in fisheye format')
                        return {'CANCELLED'}
                
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
                  
            with open('{}.temphdr'.format(os.path.join(svp['viparams']['newdir'], 'images', 'glare')), 'w') as temphdr:
                pcondcmd = "pcond -h+ -u 300 {}.hdr".format(os.path.join(svp['viparams']['newdir'], 'images', '{}-{}'.format(glnode['hdrname'], str(i + svp['liparams']['fs']))))
                Popen(shlex.split(pcondcmd), stdout = temphdr).communicate()

            with open(os.path.join(svp['viparams']['newdir'], 'images', "temp.glare"), "r") as catfile:
                psigncmd = "psign -h {} -cb 0 0 0 -cf 1 1 1".format(int(0.04 * imnode.y))
                psignrun = Popen(shlex.split(psigncmd), stdin = catfile, stdout = PIPE, stderr = PIPE) 

            with open("{}.hdr".format(os.path.join(svp['viparams']['newdir'], 'images', '{}-{}'.format(glnode['hdrname'], str(i + svp['liparams']['fs'])))), 'w') as ghdr:
                pcompcmd = "pcompos {0}.temphdr 0 0 - {1} {2}".format(os.path.join(svp['viparams']['newdir'], 'images', 'glare'), imnode.x, imnode.y*550/800)
                Popen(shlex.split(pcompcmd), stdin = psignrun.stdout, stdout = ghdr).communicate()

            os.remove(os.path.join(svp['viparams']['newdir'], 'images', 'glare.temphdr'))  
            os.remove(os.path.join(svp['viparams']['newdir'], 'images', 'temp.glare')) 
            glnode.postsim()    
                                     
        return {'FINISHED'}

class NODE_OT_Li_Fc(bpy.types.Operator):            
    bl_idname = "node.livifc"
    bl_label = "LiVi False Colour Image"
    bl_description = "False colour an image with falsecolor"
    bl_register = True
    bl_undo = False

    def execute(self, context):
        scene = context.scene
        svp = scene.vi_params
        fcnode = context.node
        fcnode.presim()
        imnode = fcnode.inputs['Image'].links[0].from_node if fcnode.inputs['Image'].links else fcnode
        lmax = '-s {}'.format(fcnode.lmax) if fcnode.lmax else '-s a'
        scaling = '' if fcnode.nscale == '0' else '-log {}'.format(fcnode.decades) 
        mult = '-m {}'.format(eval('{}{}'.format(1, fcnode.multiplier))) if fcnode.multiplier else ''
        legend = '-l {} -lw {} -lh {} {} {} {}'.format(fcnode.unit_name, fcnode.lw, fcnode.lh, lmax, scaling, mult) if fcnode.legend else ''
        bands = '-cb' if fcnode.bands else ''
        contour = '-cl {}'.format(bands) if fcnode.contour else ''
        divisions = '-n {}'.format(fcnode.divisions) if fcnode.divisions != 8 else ''
        
        for i, im in enumerate(imnode['images']): 
            fcim = os.path.join(svp['viparams']['newdir'], 'images', '{}-{}.hdr'.format(fcnode['basename'], i + svp['liparams']['fs']))
            ofile = bpy.path.abspath(fcnode.ofile) if os.path.isfile(bpy.path.abspath(fcnode.ofile)) and fcnode.overlay else bpy.path.abspath(im)
                        
            with open(fcim, 'w') as fcfile:
                if sys.platform == 'win32':
                    temp_file = os.path.join(svp['viparams']['newdir'], 'images', 'temp.hdr')
                    
                    with open(temp_file, 'w') as tfile:
                        pccmd = 'pcond -e {} "{}"'.format(fcnode.disp, os.path.abspath(im))
                        pcrun = Popen(shlex.split(pccmd), stdout = tfile, stderr = PIPE)
                    
                    for line in pcrun.stderr:
                        logentry('Pcond error: {}'.format(line))

                    poverlay = '-p {}'.format(os.path.join(context.scene['viparams']['newdir'], 'images', 'temp.hdr')) if fcnode.contour and fcnode.overlay else ''
                    fccmd = 'falsecolor -i "{}" {} -pal {} {} {} {}'.format(os.path.abspath(im), poverlay, fcnode.coldict[fcnode.colour], legend, contour, divisions) 
                    fcrun = Popen(shlex.split(fccmd), stdout=fcfile, stderr = PIPE) 

                else:
                    poverlay = '-p <(pcond -e {0} "{1}")' .format(fcnode.disp, ofile) if fcnode.contour and fcnode.overlay else ''
                    fccmd = "bash -c 'falsecolor -i \"{}\" {} -pal {} {} {} {}'".format(bpy.path.abspath(im), poverlay, fcnode.coldict[fcnode.colour], legend, contour, divisions) 
                    fcrun = Popen(shlex.split(fccmd), stdout=fcfile, stderr = PIPE)

                logentry('Running falsecolour with the command: {}'.format(fccmd))
               
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
            
class MAT_EnVi_Node(bpy.types.Operator):
    bl_idname = "material.envi_node"
    bl_label = "EnVi Material export"

    def invoke(self, context, event):
        if context.material:            
            cm = context.material
        else:
            cm = bpy.context.active_object.active_material

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

class MAT_EnVi_Node_Remove(bpy.types.Operator):
    bl_idname = "envi_node.remove"
    bl_label = "EnVi Material export"

    def invoke(self, context, event):
        for mat in bpy.data.materials:
            m = 0
            for coll in bpy.data.collections['EnVi Geometry'].children:
                for o in coll.objects:
                    if mat in [ms.material for ms in o.material_slots]:
                        m = 1
            if not m:
                print(mat.name)
#            if mat.envi_nodes and mat not in bpy.data.collections 
#            bpy.data.materials.remove(mat)
        
        return {'FINISHED'}

class NODE_OT_En_Geo(bpy.types.Operator):
    bl_idname = "node.engexport"
    bl_label = "EnVi geometry export"
    bl_context = "scene"

    def invoke(self, context, event):
        objmode()
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

class NODE_OT_En_UV(bpy.types.Operator):
    bl_idname = "node.envi_uv"
    bl_label = "EnVi Material U-Value Calculation"

    def execute(self, context):
        node = context.node

        if node.envi_con_makeup == '1':
            node.ret_uv()

        if node.envi_con_type == 'Window' and node.fclass == '2':
            node.ret_frame_uv()

        return {'FINISHED'}
    
class NODE_OT_En_Con(bpy.types.Operator, ExportHelper):
    bl_idname = "node.encon"
    bl_label = "Export"
    bl_description = "Export the scene to the EnergyPlus file format"
    bl_register = True
    bl_undo = False

    def invoke(self, context, event):
        scene = context.scene
        svp = scene.vi_params
        
        if viparams(self, scene):
            return {'CANCELLED'}
        
        svp['viparams']['vidisp'] = ''
        node = context.node
        (svp['enparams']['fs'], svp['enparams']['fe']) = (node.fs, node.fe) if node.animated else (scene.frame_current, scene.frame_current)
        locnode = node.inputs['Location in'].links[0].from_node
        
        if not os.path.isfile(locnode.weather):
            self.report({'ERROR'}, 'Location node weather file is not valid')
            node.use_custom_color = 1
            return {'CANCELLED'}
        
        node.preexport(scene)
        
        for fi, frame in enumerate(range(node.fs, node.fe + 1)):
            scene.frame_set(frame)
            
            if locnode.outputs['Parameter'].links:
                af = bpy.data.texts[locnode.outputs['Parameter'].links[0].to_node.anim_file].as_string()
                param = locnode.outputs['Parameter'].links[0].to_node.parameter

                for p in locnode.bl_rna.properties:                   
                    if p.is_skip_save:
                        if p.identifier == param:
                            for v in locnode['entries']:
                                if v[1] == af.split('\n')[fi]:
                                    try:
                                        setattr(locnode, param, v[0])
                                    except Exception as e:
                                        self.report({'ERROR'}, 'Error in parametric text file')
                                        return {'CANCELLED'}
                                        
            shutil.copyfile(locnode.weather, os.path.join(svp['viparams']['newdir'], "in{}.epw".format(frame)))
        
        scene.frame_set(node.fs)

        if context.active_object and not context.active_object.visible_get():
            if context.active_object.type == 'MESH':
                bpy.ops.object.mode_set(mode = 'OBJECT')

        enpolymatexport(self, node, locnode, envi_materials(), envi_constructions())
        node.bl_label = node.bl_label[1:] if node.bl_label[0] == '*' else node.bl_label
        node.exported, node.outputs['Context out'].hide = True, False
        node.postexport()
        return {'FINISHED'}

class NODE_OT_En_PVA(bpy.types.Operator):
    bl_idname = "node.pv_area"
    bl_label = "EnVi Material PV area calculation"

    def execute(self, context):
        node = context.node
        node['area'] = bpy.data.materials[node.id_data.name].vi_params['enparams']['area']
        return {'FINISHED'}
    
class NODE_OT_En_PVS(bpy.types.Operator):
    bl_idname = "node.pv_save"
    bl_label = "EnVi Material PV save"

    def execute(self, context):
        node = context.node
        node.save_e1ddict()
#        node['area'] = bpy.data.materials[node.id_data.name].vi_params['enparams']['area']
        return {'FINISHED'}
    
class NODE_OT_En_LayS(bpy.types.Operator):
    bl_idname = "node.lay_save"
    bl_label = "EnVi material save"

    def execute(self, context):
        node = context.node        
        node.save_laydict()
        return {'FINISHED'}
    
class NODE_OT_En_ConS(bpy.types.Operator):
    bl_idname = "node.con_save"
    bl_label = "EnVi construction save"

    def execute(self, context):
        node = context.node        
        node.save_condict()
        return {'FINISHED'}
    
class NODE_OT_En_Sim(bpy.types.Operator):
    bl_idname = "node.ensim"
    bl_label = "Simulate"
    bl_description = "Run EnergyPlus"
    bl_register = True
    bl_undo = False

    def modal(self, context, event):
        if self.pfile.check(self.percent) == 'CANCELLED':                                    
            return {self.terminate('CANCELLED', context)}
        
        while sum([esim.poll() is None for esim in self.esimruns]) < self.processors and self.e < self.lenframes:
            self.esimruns.append(Popen(self.esimcmds[self.e].split(), stderr = PIPE))
            errtext =  self.esimruns[-1].stderr.read().decode()

            if 'Fatal' in errtext:
                logentry('There is something wrong with the Energyplus installation. Check the message below')
                logentry('If using EMS a local installation of EnergyPlus is required')
                logentry(errtext)
                return {self.terminate('CANCELLED', context)}
            
            # if 'EnergyPlus Completed Successfully' not in errtext:
            #     logentry('Energyplus error: {}'.format(errtext))
            #     self.report({'ERROR'}, "Fatal error reported in the EnergyPLus error file. Check the file in Blender's text editor")

            #     for f in range(self.frame, self.frame + len(self.esimruns)):
            #         efilename = "{}{}out.err".format(self.resname, f)

            #         if os.path.isfile(os.path.join(self.nd, efilename)):             
            #             if efilename not in [im.name for im in bpy.data.texts]:
            #                 bpy.data.texts.load(os.path.join(self.nd, efilename))
            #             else:
            #                 bpy.data.texts[efilename].filepath = os.path.join(self.nd, efilename)

            #     return {self.terminate('CANCELLED', context)}
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
                    logentry('There was an error in the EnVi simulation. Check the error log in the text editor')
                
            if all([esim.poll() is not None for esim in self.esimruns]) and self.e == self.lenframes:
                for fname in [fname for fname in os.listdir('.') if fname.split(".")[0] == self.simnode.resname]:
                    os.remove(os.path.join(self.nd, fname))

                for f in range(self.frame, self.frame + self.e):
                    nfns = [fname for fname in os.listdir('.') if fname.split(".")[0] == "{}{}out".format(self.resname, f)]
                    
                    for fname in nfns:
                        os.rename(os.path.join(self.nd, fname), os.path.join(self.nd,fname.replace("eplusout", self.simnode.resname)))
                      
                    efilename = "{}{}out.err".format(self.resname, f)

                    if os.path.isfile(os.path.join(self.nd, efilename)):               
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
                    else:
                        logentry('There was an error in the EnVi simulation. Check the error log in the text editor') 
                        self.report({'ERROR'}, "Fatal error reported in the {} file. Check the file in Blender's text editor".format(efilename))
                        return {self.terminate('CANCELLED', context)}
                
                self.report({'INFO'}, "Calculation is finished.") 
                return {self.terminate('FINISHED', context)}
            else:
                return {'PASS_THROUGH'} 
        else:
            return {'PASS_THROUGH'}   
        
    def invoke(self, context, event):
        scene = context.scene
        svp = scene.vi_params
         
        if viparams(self, scene):
            return {'CANCELLED'}
        
        if shutil.which('energyplus') is None:
            self.report({'ERROR'}, "Energyplus binary is not executable")
            return {'CANCELLED'}
        
        self.frame = svp['enparams']['fs']
        self.frames = range(svp['enparams']['fs'], svp['enparams']['fe'] + 1)
        self.lenframes = len(self.frames)             
        svp['viparams']['visimcontext'] = 'EnVi'
        self.pfile = progressfile(svp['viparams']['newdir'], datetime.datetime.now(), 100)
        self.kivyrun = progressbar(os.path.join(svp['viparams']['newdir'], 'viprogress'), 'EnergyPlus Results')
        wm = context.window_manager
        self._timer = wm.event_timer_add(1, window = context.window)
        wm.modal_handler_add(self)
        self.simnode = context.node
        self.connode = self.simnode.inputs[0].links[0].from_node.name
        self.simnode.presim(context)        
        self.expand = "-x" if svp['viparams'].get('hvactemplate') else ""
        self.resname = (self.simnode.resname, 'eplus')[self.simnode.resname == '']
        os.chdir(svp['viparams']['newdir'])
        self.esimcmds = ["energyplus {0} -w in{1}.epw -p {2} in{1}.idf".format(self.expand, frame, ('{}{}'.format(self.resname, frame))) for frame in self.frames] 
        self.esimruns = []
        self.simnode.run = 1
        self.processors = self.simnode.processors if self.simnode.mp else 1
        self.percent = 0
        self.e = 0
        self.nd = svp['viparams']['newdir']
        return {'RUNNING_MODAL'}
    
    def terminate(self, condition, context):
        scene = context.scene
        svp = scene.vi_params
        print('term')
        self.simnode.postsim(self, condition)
        
        for f in range(self.frame, self.frame + self.e):
            efilename = "{}{}out.err".format(self.resname, f)

            if os.path.isfile(os.path.join(self.nd, efilename)):               
                if efilename not in [im.name for im in bpy.data.texts]:
                    bpy.data.texts.load(os.path.join(self.nd, efilename))
                else:
                    bpy.data.texts[efilename].filepath = os.path.join(self.nd, efilename)

        if condition == 'FINISHED':
            svp['viparams']['resnode'] = '{}@{}'.format(self.simnode.name, self.simnode.id_data.name)
            svp['viparams']['connode'] = '{}@{}'.format(self.connode, self.simnode.id_data.name)
            svp['viparams']['vidisp'] = 'en'

        self.kivyrun.kill() 
             
        for es in self.esimruns:
            if es.poll() is None:
                es.kill()
                
        return condition
    
class OBJECT_OT_VIGridify2(bpy.types.Operator):
    ''''''
    bl_idname = "object.vi_gridify2"
    bl_label = "VI Gridify"
    bl_options = {"REGISTER", 'UNDO'}
    
    rotate: bpy.props.FloatProperty(name = 'Rotation', default = 0, min = 0, max = 360) 
    us: bpy.props.FloatProperty(default = 0.6, min = 0.01) 
    acs: bpy.props.FloatProperty(default = 0.6, min = 0.01) 
    
    @classmethod
    def poll(cls, context):
        obj = context.active_object
        return (obj and obj.type == 'MESH')
    
    def execute(self, context):
        self.o = bpy.context.active_object
        mesh = bmesh.from_edit_mesh(self.o.data)
        mesh.transform(self.o.matrix_world)
        mesh.faces.ensure_lookup_table()
        mesh.verts.ensure_lookup_table()
        fs = [f for f in mesh.faces[:] if f.select]
        self.upv = fs[0].calc_tangent_edge().copy().normalized()
        self.norm = fs[0].normal.copy()
        self.acv = self.upv.copy()
        eul = Euler(radians(-90) * self.norm, 'XYZ')
        self.acv.rotate(eul)
        self.acv = self.upv.cross(self.norm)
        rotation = Euler(radians(self.rotate) * self.norm, 'XYZ')
        self.upv.rotate(rotation)
        self.acv.rotate(rotation)
        vertdots = [Vector.dot(self.upv, vert.co) for vert in fs[0].verts]
        vertdots2 = [Vector.dot(self.acv, vert.co) for vert in fs[0].verts]
        svpos = fs[0].verts[vertdots.index(min(vertdots))].co
        svpos2 = fs[0].verts[vertdots2.index(min(vertdots2))].co
        res1, res2, ngs1, ngs2, gs1, gs2 = 1, 1, self.us, self.acs, self.us, self.acs
        vs = fs[0].verts[:]
        es = fs[0].edges[:]        
        gs = vs + es + fs
          
        while res1:
            res = bmesh.ops.bisect_plane(mesh, geom = gs, dist = 0.001, plane_co = svpos + ngs1 * self.upv, plane_no = self.upv, use_snap_center = 0, clear_outer = 0, clear_inner = 0)
            res1 = res['geom_cut']
            dissvs = [v for v in res1 if isinstance(v, bmesh.types.BMVert) and len(v.link_edges) == 2 and v.calc_edge_angle(1) == 0.0]
            bmesh.ops.dissolve_verts(mesh, verts = dissvs)
                    
            gs = mesh.verts[:] + mesh.edges[:] + [v for v in res['geom'] if isinstance(v, bmesh.types.BMFace)]
            ngs1 += gs1
    
        while res2:
            res = bmesh.ops.bisect_plane(mesh, geom = gs, dist = 0.001, plane_co = svpos2 + ngs2 * self.acv, plane_no = self.acv, use_snap_center = 0, clear_outer = 0, clear_inner = 0)
            res2 = res['geom_cut']
            dissvs = [v for v in res2 if isinstance(v, bmesh.types.BMVert) and len(v.link_edges) == 2 and v.calc_edge_angle(1) == 0.0]
            bmesh.ops.dissolve_verts(mesh, verts = dissvs)
            gs = mesh.verts[:] + mesh.edges[:] + [v for v in res['geom'] if isinstance(v, bmesh.types.BMFace)]
            ngs2 += gs2
        
        mesh.transform(self.o.matrix_world.inverted())
        bmesh.update_edit_mesh(self.o.data)
        mesh.free()
        return {'FINISHED'}
    
class OBJECT_OT_GOct(bpy.types.Operator):
    ''''''
    bl_idname = "object.vi_genoct"
    bl_label = "Octree Generator"
    bl_options = {"REGISTER", 'UNDO'}
    
    @classmethod
    def poll(cls, context):
        obj = context.active_object
        return (obj and obj.type == 'MESH')
    
    def execute(self, context):
        scene = context.scene
        ovp = context.object.vi_params
        gen_octree(scene, context.object, self, ovp.mesh, ovp.triangulate)
        return {'FINISHED'}

class OBJECT_OT_Embod(bpy.types.Operator):
    '''Calculates embodied energy based on object volume'''
    bl_idname = "object.vi_embodied"
    bl_label = "Embodied carbon"
    bl_options = {"REGISTER", 'UNDO'}
    
    @classmethod
    def poll(cls, context):
        obj = context.active_object
        return (obj and obj.type == 'MESH')
    
    def execute(self, context):
        dp = bpy.context.evaluated_depsgraph_get()
        o = context.object
        ovp = o.vi_params
        scene = context.scene
        bm = bmesh.new()
        bm.from_object(o, dp)
        bm.transform(o.matrix_world)
        bm.normal_update()

        if all([e.is_manifold for e in bm.edges]):
            envi_ec = envi_embodied()
            vol = bm.calc_volume()
            ovp['ecdict'] = envi_ec.propdict[ovp.embodiedtype][ovp.embodiedclass][ovp.embodiedmat]
            ovp['ecdict']['ec'] = float(ovp['ecdict']['eckg']) * float(ovp['ecdict']['density']) * vol
            bm.free()
        else:
            self.report({'ERROR'},"You cannot calculate embodied carbon on a non-manifold mesh")
            bm.free()
            return {'CANCELLED'}
        
        return {'FINISHED'}
    
class NODE_OT_Chart(bpy.types.Operator, ExportHelper):
    bl_idname = "node.chart"
    bl_label = "Chart"
    bl_description = "Create a 2D graph from the results file"
    bl_register = True
    bl_undo = True

    def invoke(self, context, event):
        node = context.node
        innodes = list(OrderedDict.fromkeys([inputs.links[0].from_node for inputs in node.inputs if inputs.links]))
        rl = innodes[0]['reslists']
        zrl = list(zip(*rl))
        year = context.scene.vi_params.year
        
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

class NODE_OT_HMChart(bpy.types.Operator, ExportHelper):
    bl_idname = "node.hmchart"
    bl_label = "Heatmap"
    bl_description = "Create a 2D heatmap from the results file"
    bl_register = True
    bl_undo = True



#        maxx, maxy, maxz = max(self.x), max(self.y), max(self.z)
#        minx, miny, minz = min(self.x), min(self.y), min(self.z)
#        print(len(x), len(y), len(self.z))
        # self.fig, self.ax = plt.subplots(figsize=(12, 6))   
        # self.plt.title('Test', size = 20).set_position([.5, 1.025])
        # self.plt.xlabel('Test', size = 18)
        # self.plt.ylabel('Test', size = 18)
        # self.plt.pcolormesh(self.x, self.y, self.z, shading='auto', vmin=minz, vmax=maxz)
        # cbar = self.plt.colorbar(use_gridspec=True, pad = 0.01)
        # cbar.set_label(label='Test',size=18)
        # cbar.ax.tick_params(labelsize=16)
        # self.plt.axis([minx,maxx,miny,maxy])
        # self.plt.xticks(size = 16)
        # self.plt.yticks(size = 16)
        # self.fig.tight_layout()
        # self.plt.show()

    def invoke(self, context, event):
        node = context.node
        node.dupdate(context)
        # innodes = list(OrderedDict.fromkeys([inputs.links[0].from_node for inputs in node.inputs if inputs.links]))
        # rl = innodes[0]['reslists']
        # zrl = list(zip(*rl))
        year = context.scene.vi_params.year
        
        if not mp:
            self.report({'ERROR'},"Matplotlib cannot be found by the Python installation used by Blender")
            return {'CANCELLED'}

        # Sdate = dt.fromordinal(dt(year, 1, 1).toordinal() + node['Start'] - 1)
        # Edate = dt.fromordinal(dt(year, 1, 1).toordinal() + node['End'] - 1)for x in range(20):
        hmchart_disp(self, plt, node, context.scene.vi_params.vi_leg_col)
        return {'FINISHED'}

class NODE_OT_MInfo(bpy.types.Operator):
    bl_idname = "node.metinfo"
    bl_label = "Graphic"
    bl_description = "Creates an Infographic of the choden metric"
    bl_register = True
    bl_undo = False

    def execute(self, context):
        metnode = context.node
        svg_str = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
    <svg width="300" height="300" viewBox="0 0 300 300" id="smile" version="1.1">
        <path
            style="fill:#aaaaff"
            d="M 150,0 A 150,150 0 0 0 0,150 150,150 0 0 0 150,300 150,150 0 0 0 
                300,150 150,150 0 0 0 150,0 Z M 72,65 A 21,29.5 0 0 1 93,94.33 
                21,29.5 0 0 1 72,124 21,29.5 0 0 1 51,94.33 21,29.5 0 0 1 72,65 Z 
                m 156,0 a 21,29.5 0 0 1 21,29.5 21,29.5 0 0 1 -21,29.5 21,29.5 0 0 1 
                -21,-29.5 21,29.5 0 0 1 21,-29.5 z m -158.75,89.5 161.5,0 c 0,44.67 
                -36.125,80.75 -80.75,80.75 -44.67,0 -80.75,-36.125 -80.75,-80.75 z"
        />
    </svg>
    """

        svg_bytes = bytearray(svg_str, encoding='utf-8')
        qimage = QImage.fromData(svg_bytes)
        
        rgba= ndarray(shape = (300, 300, 4), dtype = uint8)

        for x in range(300):
            for y in range(300):
                rgba[299 - y][x] = QColor(qimage.pixel(x, y)).getRgbF()

        imname = "test.png"
        ipheight, ipwidth = 300, 300     
        
        if imname not in [im.name for im in bpy.data.images]:
            bpy.ops.image.new(name=imname, width=ipwidth, height=ipheight, color=(0, 0, 0, 0), alpha=True, generated_type='BLANK', float=False, use_stereo_3d=False)
            im = bpy.data.images[imname]

        else:
            im = bpy.data.images[imname] 
            im.gl_free()
            im.buffers_free()

            if (im.generated_width, im.generated_height) != (ipwidth, ipheight):
                im.generated_width = ipwidth
                im.generated_height = ipheight

            if im.size[:] != (ipwidth, ipheight):
                im.scale(ipwidth, ipheight)
            
        im.pixels.foreach_set(rgba.ravel().astype(float32))

        # Opens new image window
        # area = bpy.context.area
        # t = area.type
        # area.type = 'IMAGE_EDITOR'
        # bpy.ops.screen.area_dupli('INVOKE_DEFAULT')
        #win = bpy.context.window_manager.windows[-1]
        #win.screen.areas[0].spaces[0].show_region_header = 0
        #win.screen.areas[0].spaces[0].show_region_ui = 0
        # area.type = t


class VIEW3D_OT_EnDisplay(bpy.types.Operator):
    bl_idname = "view3d.endisplay"
    bl_label = "EnVi display"
    bl_description = "Display the EnVi results"
    bl_options = {'REGISTER'}
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

# Node utilities from Matalogue    
class TREE_OT_goto_mat(bpy.types.Operator):
    'Show the EnVi nodes for this material'
    bl_idname = 'tree.goto_mat'
    bl_label = 'Go To EnVi Material Group'

    mat: bpy.props.StringProperty(default="")

    def execute(self, context):
        context.space_data.tree_type = 'EnViMatN'
        mat = bpy.data.materials[self.mat]
        context.space_data.node_tree = mat.vi_params.envi_nodes
        context.space_data.node_tree.name = mat.name
        objs_with_mat = 0
        active_set = False
        
        for obj in context.view_layer.objects:
            obj_materials = [slot.material for slot in obj.material_slots]
            if mat in obj_materials:
                objs_with_mat += 1
                obj.select_set(True)
                if not active_set:  # set first object as active
                    active_set = True
                    context.view_layer.objects.active = obj
                    if mat != obj.active_material:
                        for i, x in enumerate(obj.material_slots):
                            if x.material == mat:
                                obj.active_material_index = i
                                break
            else:
                obj.select_set(False)

        if objs_with_mat == 0:
            self.report({'WARNING'}, "No objects in this scene use '" + mat.name + "' material")

        return {'FINISHED'}
    
class TREE_OT_goto_group(bpy.types.Operator):
    'Show the nodes inside this group'
    bl_idname = 'tree.goto_group'
    bl_label = 'Go To Group'

    tree_type: bpy.props.StringProperty(default="")
    tree: bpy.props.StringProperty(default="")

    def execute(self, context):
        try:  # Go up one group as many times as possible - error will occur when the top level is reached
            while True:
                bpy.ops.node.tree_path_parent()
        except:
            pass
        context.space_data.tree_type = self.tree_type
        context.space_data.path.append(bpy.data.node_groups[self.tree])
        context.space_data.node_tree = bpy.data.node_groups[self.tree]
        context.space_data.node_tree.use_fake_user = 1
        return {'FINISHED'}
    
class NODE_OT_CSV(bpy.types.Operator, ExportHelper):
    bl_idname = "node.csvexport"
    bl_label = "Export a CSV file"
    bl_description = "Select the CSV file to export"
    filename = "results"
    filename_ext = ".csv"
    filter_glob: bpy.props.StringProperty(default="*.csv", options={'HIDDEN'})
    bl_register = True
    bl_undo = True

    def draw(self,context):
        layout = self.layout
        row = layout.row()
        row.label(text="Specify the CSV export file with the file browser", icon='WORLD_DATA')

    def execute(self, context):
        node = self.node
        resstring = ''
        resnode = node.inputs['Results in'].links[0].from_node
        rl = resnode['reslists']
        zrl = list(zip(*rl))

        if len(set(zrl[0])) > 1 and node.animated:
            resstring = ''.join(['{} {},'.format(r[2], r[3]) for r in rl if r[0] == 'All']) + '\n'
            metriclist = list(zip(*[r.split() for ri, r in enumerate(zrl[4]) if zrl[0][ri] == 'All']))
        else:
            resstring = ''.join(['{} {} {},'.format(r[0], r[2], r[3]) for r in rl if r[0] != 'All']) + '\n'
            metriclist = list(itertools.zip_longest(*[r.split() for ri, r in enumerate(zrl[4]) if zrl[0][ri] != 'All'], fillvalue = ''))

        for ml in metriclist:
            resstring += ''.join(['{},'.format(m) for m in ml]) + '\n'

        resstring += '\n'

        with open(self.filepath, 'w') as csvfile:
            csvfile.write(resstring)
            
        return {'FINISHED'}

    def invoke(self, context, event):
        self.node = context.node
        svp = context.scene.vi_params

        if self.filepath.split('.')[-1] not in ('csv', 'CSV'):
            self.filepath = os.path.join(svp['viparams']['newdir'], svp['viparams']['filebase'] + '.csv')   

        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}   
    
# Openfoam operators
class NODE_OT_Flo_Case(bpy.types.Operator):
    bl_idname = "node.flovi_case"
    bl_label = "Case export"
    bl_description = "Export an Openfoam case"
    bl_register = True
    bl_undo = False

    def execute(self, context):
        scene = context.scene

        if viparams(self, scene):
            return {'CANCELLED'}

        svp = scene.vi_params
        casenode = context.node# if context.node.bl_label == "FloVi BlockMesh" else context.node.inputs[0].links[0].from_node
#        bmos = [o for o in scene.objects if o.vi_params.vi_type == '2']
        casenode.pre_case(context)
        dobs = [o for o in bpy.data.objects if o.vi_params.vi_type == '2']
        gobs = [o for o in bpy.data.objects if o.vi_params.vi_type == '3']
        obs = dobs + gobs

        for f in os.listdir(svp['flparams']['offilebase']):
            try:
                os.remove(os.path.join(svp['flparams']['offilebase'], f))
            except:
                pass
        
        for f in os.listdir(svp['flparams']['ofcfilebase']):
            try:
                os.remove(os.path.join(svp['flparams']['ofcfilebase'], f))
            except:
                pass
        
        for root, dirs, files in os.walk(os.path.join(svp['flparams']['offilebase'], 'postProcessing')):
            for d in dirs:
                try:
#                    if float(d) != svp['flparams']['st']:
                    shutil.rmtree(os.path.join(root, d))
                except:
                    pass
                # try:
                #     if float(d) != svp['flparams']['st']:
                #         shutil.rmtree(os.path.join(root, d))
                # except:
                #     pass
                
#                    shutil.rmtree()
#        if not os.path.isdir(os.path.join(svp['flparams']['offilebase'], node.stime)):
#            os.makedirs(os.path.join(svp['flparams']['offilebase'], node.stime))

        
        if len(dobs) != 1:
            self.report({'ERROR'},"One, and only one object with the CFD Domain property is allowed")
            return {'CANCELLED'}
        elif [f.material_index for f in dobs[0].data.polygons if f.material_index + 1 > len(dobs[0].data.materials)]:
            self.report({'ERROR'},"Not every domain face has a material attached")
            logentry("Not every face has a material attached")
            return {'CANCELLED'}

        if casenode.turbulence == 'laminar':
            svp['flparams']['solver'] = 'icoFoam' 
            svp['flparams']['solver_type'] = 'if'  
            svp['flparams']['residuals'] = ['p', 'Ux', 'Uy', 'Uz']          
        elif casenode.transience == '0':
            if casenode.buoyancy: 
                svp['flparams']['pref'] = casenode.pabsval               
                svp['flparams']['solver'] = 'buoyantSimpleFoam' 
                if not casenode.buossinesq:
                    svp['flparams']['solver_type'] = 'bsf'
                    if casenode.radiation:
                        svp['flparams']['residuals'] = ['G', 'p_rgh', 'k', 'epsilon', 'h']
                    else:
                        svp['flparams']['residuals'] = ['p_rgh', 'k', 'epsilon', 'h']
                else:
                    svp['flparams']['solver_type'] = 'bbsf'
                    if casenode.radiation:
                        svp['flparams']['residuals'] = ['G', 'Ux', 'Uy', 'Uz', 'k', 'epsilon', 'e', 'p_rgh']
                    else:
                        svp['flparams']['residuals'] = ['Ux', 'Uy', 'Uz', 'k', 'epsilon', 'e', 'p_rgh']
                
            else:
                svp['flparams']['solver'] = 'simpleFoam' 
                svp['flparams']['pref'] = casenode.pnormval
                if casenode.turbulence == 'kEpsilon':
                    svp['flparams']['residuals'] = ['p', 'Ux', 'Uy', 'Uz', 'k', 'epsilon']
                if casenode.turbulence == 'kOmega':  
                    svp['flparams']['residuals'] = ['p', 'Ux', 'Uy', 'Uz', 'k', 'omega']
                if casenode.turbulence == 'SpalartAllmaras':  
                    svp['flparams']['residuals'] = ['p', 'Ux', 'Uy', 'Uz', 'nuTilda']
                svp['flparams']['solver_type'] = 'sf'
               
        elif casenode.transience == '1':
            if casenode.buoyancy:
                svp['flparams']['solver'] = 'buoyantPimpleFoam'
                if not casenode.buossinesq:
                    svp['flparams']['solver_type'] = 'bpf'
                else:
                    svp['flparams']['solver_type'] = 'bbpf'
            else:
                svp['flparams']['solver'] = 'pimpleFoam' 
                svp['flparams']['solver_type'] = 'pf'
        
        svp['flparams']['st'] = casenode.stime
        svp['flparams']['presid'] = casenode.presid
        svp['flparams']['uresid'] = casenode.uresid
        svp['flparams']['keoresid'] = casenode.keoresid
        
        if casenode.parametric:
            frames = range(casenode.frame_start, casenode.frame_end + 1)    
        else:
            frames = [scene.frame_current]

        svp['flparams']['start_frame'] = frames[0]
        svp['flparams']['end_frame'] = frames[-1]
        
        # for frame in frames:
        #     frame_dir = os.path.join(svp['flparams']['offilebase'], str(frame))
        #     frame_system_dir
        #     os.makedirs(frame_dir)
        with open(os.path.join(svp['flparams']['ofsfilebase'], 'controlDict'), 'w') as cdfile:
            cdfile.write(fvcdwrite(casenode.solver, casenode.stime, casenode.dtime, casenode.etime))
        with open(os.path.join(svp['flparams']['ofsfilebase'], 'fvSolution'), 'w') as fvsolfile:
            fvsolfile.write(fvsolwrite(casenode, svp['flparams']['solver_type']))
        with open(os.path.join(svp['flparams']['ofsfilebase'], 'fvSchemes'), 'w') as fvschfile:
            fvschfile.write(fvschwrite(casenode, svp['flparams']['solver_type'])) 

        if casenode.turbulence == 'laminar':
            with open(os.path.join(svp['flparams']['ofcfilebase'], 'transportProperties'), 'w') as tpfile:    
                tpfile.write(fvtpwrite(casenode, svp['flparams']['solver_type']))     
        else:
            with open(os.path.join(svp['flparams']['ofcfilebase'], 'momentumTransport'), 'w') as mtfile:    
                mtfile.write(fvmtwrite(casenode, svp['flparams']['solver_type']))

            if casenode.buoyancy: 
                with open(os.path.join(svp['flparams']['ofcfilebase'], 'pRef'), 'w') as pfile:
                    pfile.write(fvprefwrite(casenode, svp['flparams']['solver_type']))  
                with open(os.path.join(svp['flparams']['ofcfilebase'], 'thermophysicalProperties'), 'w') as tppfile:    
                    tppfile.write(fvtppwrite(casenode, svp['flparams']['solver_type']))
                # with open(os.path.join(svp['flparams']['ofcfilebase'], 'momentumTransport'), 'w') as mtfile:    
                #     mtfile.write(fvmtwrite(casenode, svp['flparams']['solver_type']))
                with open(os.path.join(svp['flparams']['ofcfilebase'], 'g'), 'w') as gfile:    
                    gfile.write(fvgwrite())   
                if casenode.radiation:    
                    with open(os.path.join(svp['flparams']['ofcfilebase'], 'radiationProperties'), 'w') as rpfile:    
                        rpfile.write(fvrpwrite(casenode, svp['flparams']['solver_type']))
            else:
                with open(os.path.join(svp['flparams']['ofcfilebase'], 'transportProperties'), 'w') as tpfile:    
                    tpfile.write(fvtpwrite(casenode, svp['flparams']['solver_type']))

        casenode.post_case()
        return {'FINISHED'}

class NODE_OT_Flo_BM(bpy.types.Operator):
    bl_idname = "node.flovi_bm"
    bl_label = "Blockmesh export"
    bl_description = "Export an Openfoam blockmesh"
    bl_register = True
    bl_undo = False

    def execute(self, context):
        scene = context.scene
        svp = scene.vi_params
        expnode = context.node if context.node.bl_label == "FloVi BlockMesh" else context.node.inputs[0].links[0].from_node
        bmos = [o for o in scene.objects if o.vi_params.vi_type == '2']
        
        if viparams(self, scene):
            return {'CANCELLED'} 
        
        if len(bmos) != 1:
            self.report({'ERROR'},"One and only one object with the CFD Domain property is allowed")
            return {'CANCELLED'}
        elif [f.material_index for f in bmos[0].data.polygons if f.material_index + 1 > len(bmos[0].data.materials)]:
            self.report({'ERROR'},"Not every domain face has a material attached")
            logentry("Not every face has a material attached")
            return {'CANCELLED'}
        with open(os.path.join(svp['flparams']['ofsfilebase'], 'controlDict'), 'w') as cdfile:
            cdfile.write(fvcdwrite("simpleFoam", 0.005, 5))
        with open(os.path.join(svp['flparams']['ofsfilebase'], 'fvSolution'), 'w') as fvsolfile:
            fvsolfile.write(fvsolwrite(expnode))
        with open(os.path.join(svp['flparams']['ofsfilebase'], 'fvSchemes'), 'w') as fvschfile:
            fvschfile.write(fvschwrite(expnode))
        with open(os.path.join(svp['flparams']['ofcpfilebase'], 'blockMeshDict'), 'w') as bmfile:
            bmfile.write(fvbmwrite(bmos[0], expnode))

        call(("blockMesh", "-case", "{}".format(scene['flparams']['offilebase'])))
        fvblbmgen(bmos[0].data.materials, open(os.path.join(scene['flparams']['ofcpfilebase'], 'faces'), 'r'), open(os.path.join(scene['flparams']['ofcpfilebase'], 'points'), 'r'), open(os.path.join(scene['flparams']['ofcpfilebase'], 'boundary'), 'r'), 'blockMesh')
        expnode.export()
        return {'FINISHED'}
    
class NODE_OT_Flo_NG(bpy.types.Operator):
    bl_idname = "node.flovi_ng"
    bl_label = "NetGen export"
    bl_description = "Create a Netgen mesh"
    bl_register = True
    bl_undo = False

    def execute(self, context):
        vi_prefs = bpy.context.preferences.addons[__name__.split('.')[0]].preferences
        scene = context.scene
        svp = scene.vi_params
        vl = context.view_layer
        expnode = context.node
        meshcoll = create_coll(context, 'FloVi Mesh')
        clear_coll(context, meshcoll)
        SetNumThreads(int(svp['viparams']['nproc']))
        maxh = expnode.maxcs
        st = '0'
        totmesh = Mesh()  
        meshes = []
        mats = []
        dp = bpy.context.evaluated_depsgraph_get()
        dobs = [o for o in bpy.data.objects if o.vi_params.vi_type == '2' and o.visible_get()]
        gobs = [o for o in bpy.data.objects if o.vi_params.vi_type == '3' and o.visible_get()]
        obs = dobs + gobs

        for ob in obs:
            bm = bmesh.new()
            bm.from_object(ob, dp)

            if not all([e.is_manifold for e in bm.edges]) or not all([v.is_manifold for v in bm.verts]):
                bm.free()
                logentry('FloVi error: {} is not manifold'.format(ob.name))
                self.report({'ERROR'},'FloVi error: {} is not manifold'.format(ob.name))
                return {'CANCELLED'}

            bm.free()

        if expnode.geo_join and gobs:
            for d in gobs[1:]:
                ubool = gobs[0].modifiers.new(name='union', type='BOOLEAN')
                ubool.object = d
                ubool.operation = 'UNION'
                bpy.ops.object.modifier_apply(modifier=ubool.name)

            gobs = [gobs[0]]
            obs = dobs + gobs

            if expnode.d_diff:
                dbool = dobs[0].modifiers.new(name='diff', type='BOOLEAN')
                dbool.object = gobs[0]
                dbool.operation = 'DIFFERENCE'
                bpy.ops.object.modifier_apply(modifier=dbool.name)
                obs = [dobs[0]]

        omats = []
        mns = [0]
        fds = []
        pdm_error = 0

        for ob in obs:
            mis = empty(len(ob.data.polygons), dtype=uint8)
            ob.data.polygons.foreach_get('material_index', mis)

            try:
                omats.append([ob.material_slots[i].material for i in set(mis)])
            except:
                logentry('FloVi error: {} has missing materials'.format(ob.name))
                self.report({'ERROR'},'FloVi error: {} has missing materials'.format(ob.name))
                return {'CANCELLED'}

            mns.append(len(set(mis)))

        fomats = [item for sublist in omats for item in sublist]
        i= 0
        
        for mis, mats in enumerate(omats):    
            for mi, mat in enumerate(mats):
                if not mis:
                    fd = FaceDescriptor(bc = i, domin = 1, surfnr = i + 1)
                else:
                    fd = FaceDescriptor(bc = i, domin = 0, domout = 1, surfnr = i + 1)

                fd = totmesh.Add(fd)
                fds.append(fd)
                totmesh.SetBCName(fd, mat.name)
                i += 1
                
        for oi, o in enumerate(obs):
            mp = MeshingParameters(maxh=maxh, yangle = expnode.yang, grading = expnode.grading, 
            optsteps2d = expnode.optimisations, optsteps3d = expnode.optimisations, delaunay = True)
            bm = bmesh.new()
            bm.from_object(o, dp)
            bm.transform(o.matrix_world)
            bmesh.ops.recalc_face_normals(bm, faces = bm.faces)
            bm.normal_update()
            tris = bm.calc_loop_triangles()
        
            with open(os.path.join(svp['flparams']['offilebase'], '{}.stl'.format(o.name)), 'w') as stlfile:
                stlfile.write('solid\n')
                
                for tri in tris:
                    stlfile.write('facet normal {0[0]} {0[1]} {0[2]}\nouter loop\n'.format(tri[0].face.normal))
                    
                    for t in tri:
                        stlfile.write('vertex {0[0]} {0[1]} {0[2]}\n'.format(t.vert.co))
                        
                    stlfile.write('end loop\nend facet\n')
                stlfile.write('endsolid\n')
            
            geo = STLGeometry(os.path.join(svp['flparams']['offilebase'], '{}.stl'.format(o.name)))

            for v in bm.verts:
                mp.RestrictH(x=v.co[0],y=v.co[1],z=v.co[2],h=max([o.material_slots[f.material_index].material.vi_params.flovi_ng_max for f in v.link_faces]))
            
            for e in bm.edges:
                if e.calc_length() > 2 * max([o.material_slots[f.material_index].material.vi_params.flovi_ng_max for f in e.link_faces]):
                    segs = int(e.calc_length()/max([o.material_slots[f.material_index].material.vi_params.flovi_ng_max for f in e.link_faces])) + 1
                    
                    for s in range(1, segs):
                        vco = e.verts[0].co + (e.verts[1].co - e.verts[0].co) * s/segs
                        mp.RestrictH(x=vco[0],y=vco[1],z=vco[2],h=min([o.material_slots[f.material_index].material.vi_params.flovi_ng_max for f in v.link_faces]))
                    
            with TaskManager():
                m = geo.GenerateMesh(mp = mp)#'/home/ryan/Store/OneDrive/Blender28/flovi1/Openfoam/meshsize.msz')
    
                for el in m.Elements2D():
                    fpoint = [sum(m[v].p[x]/3 for v in el.vertices) for x in (0, 1, 2)]
                    fnorm = mathutils.geometry.normal([m[v].p for v in el.vertices]) 
                    intersect = 0
                    
                    for face in bm.faces:
                        if bmesh.geometry.intersect_face_point(face, fpoint) and abs(mathutils.geometry.distance_point_to_plane(fpoint, face.calc_center_median(), face.normal)) < expnode.pcorr and abs(fnorm.dot(face.normal)) > expnode.acorr:
                            el.index = omats[oi].index(o.material_slots[face.material_index].material) + 1 + sum(mns[:oi + 1]) 
                            intersect = 1
                            break
                    if not intersect:
                        el.index = 1

                meshes.append(m) 
                m.Save(os.path.join(svp['flparams']['offilebase'], '{}_surface.vol'.format(o.name)))

            bm.free()
            
        for mi, m in enumerate(meshes):   
            pmap1 = { }

            for e in m.Elements2D():
                for v in e.vertices:
                    if (v not in pmap1):
                        pmap1[v] = totmesh.Add(m[v])
                        
                totmesh.Add(Element2D(e.index, [pmap1[v] for v in e.vertices]))

        totmesh.Save(os.path.join(svp['flparams']['offilebase'], 'ng_surf.vol'))

        try:
            with TaskManager():    
                totmesh.GenerateVolumeMesh()

            logentry("Netgen mesh generated")            
            # The below would create a boundary layer but this is nor currently supported in Netgen Python interface
            # totmesh.BoundaryLayer(boundary = 1, thickness = 0.02, material = 't')
            totmesh.Save(os.path.join(svp['flparams']['offilebase'], 'ng.vol'))
            totmesh.Export(os.path.join(svp['flparams']['offilebase'], 'ng.mesh'), format='Neutral Format')

            if sys.platform == 'linux' and os.path.isdir(vi_prefs.ofbin):           
                subprocess.Popen(shlex.split('foamExec netgenNeutralToFoam -case {} {}'.format(svp['flparams']['offilebase'], os.path.join(svp['flparams']['offilebase'], 'ng.mesh')))).wait()
        
                if not os.path.isdir(os.path.join(svp['flparams']['offilebase'], st, 'polyMesh')):
                    os.makedirs(os.path.join(svp['flparams']['offilebase'], st, 'polyMesh'))
            
            elif sys.platform == 'darwin' and os.path.isdir(vi_prefs.ofbin): 
                print("OSX command to open openfoam docker image: {}".format("docker container run -ti --rm -v $PWD:/data -w /data openfoamplus/of_v2012_centos73:release /bin/bash"))
            
            if os.path.isfile(os.path.join(svp['flparams']['offilebase'], 'constant', 'polyMesh', 'boundary')): 
                with open(os.path.join(svp['flparams']['offilebase'], 'constant', 'polyMesh', 'boundary'), 'r') as bfile:
                    nf = []
                    ns = []
        
                    for line in bfile.readlines():
                        if 'nFaces' in line:
                            nf.append(int(line.split()[1].strip(';')))
                        if 'startFace' in line:
                            ns.append(int(line.split()[1].strip(';')))
        
                with open(os.path.join(svp['flparams']['offilebase'], 'constant', 'polyMesh', 'boundary'), 'w') as bfile:
                    bfile.write(ofheader) 
                    cl = 'polyBoundaryMesh' if expnode.bl_label == 'FloVi NetGen' else 'BoundaryMesh'
                    loc = 'constant/polyMesh' if expnode.bl_label == 'FloVi NetGen' else 'Mesh'
                    bfile.write(write_ffile(cl, loc, 'boundary'))
                    bfile.write('// **\n\n{}\n(\n'.format(len(ns)))
                    omi = 0
                    
                    for mi, mats in enumerate(omats):
                        for m in mats:
                            if omi < len(ns):
                                bfile.write(write_bound(obs[mi], m, ns[omi], nf[omi]))
                                omi += 1
                    bfile.write(')\n\n// **\n')
                    
                for file in os.listdir(svp['flparams']['ofcpfilebase']):
                    shutil.copy(os.path.join(svp['flparams']['ofcpfilebase'], file), os.path.join(svp['flparams']['offilebase'], st, 'polyMesh'))
                                    
                if expnode.poly and sys.platform == 'linux' and os.path.isdir(vi_prefs.ofbin):
                    os.chdir(svp['flparams']['offilebase'])
                    pdm = Popen(shlex.split('foamExec polyDualMesh -case {} -noFunctionObjects -noFields -overwrite {}'.format(svp['flparams']['offilebase'], expnode.yang)), stdout = PIPE, stderr = PIPE)

                    for line in pdm.stdout:
                        if 'FOAM aborting' in line.decode():
                            logentry('polyDualMesh error. Check the mesh in Netgen')
                            pdm_error = 1
                    
                    if not pdm_error:
                        Popen(shlex.split('foamExec combinePatchFaces -overwrite -case {} {}'.format(svp['flparams']['offilebase'], expnode.yang))).wait()
                        cm = Popen(shlex.split('foamExec checkMesh -case {}'.format(svp['flparams']['offilebase'])), stdout = PIPE)

                        for line in cm.stdout:
                            if '***Error' in line.decode():
                                logentry('Mesh errors: {}'.format(line))
                    
                    for entry in os.scandir(os.path.join(svp['flparams']['offilebase'], st, 'polyMesh')):
                        if entry.is_file():
                            shutil.copy(os.path.join(svp['flparams']['offilebase'], st, 'polyMesh', entry.name), os.path.join(svp['flparams']['ofcpfilebase']))
                    
                    with open(os.path.join(svp['flparams']['offilebase'], 'constant', 'polyMesh', 'boundary'), 'r') as bfile:
                        nf = []
                        ns = []
            
                        for line in bfile.readlines():
                            if 'nFaces' in line:
                                nf.append(int(line.split()[1].strip(';')))
                            if 'startFace' in line:
                                ns.append(int(line.split()[1].strip(';')))
            
                    with open(os.path.join(svp['flparams']['offilebase'], 'constant', 'polyMesh', 'boundary'), 'w') as bfile:
                        bfile.write(ofheader) 
                        cl = 'polyBoundaryMesh' if expnode.bl_label == 'FloVi NetGen' else 'BoundaryMesh'
                        loc = 'constant/polyMesh' if expnode.bl_label == 'FloVi NetGen' else 'Mesh'
                        bfile.write(write_ffile(cl, loc, 'boundary'))
                        bfile.write('// **\n\n{}\n(\n'.format(len(ns)))
                        omi = 0
                        
                        for mi, mats in enumerate(omats):
                            for m in mats:
                                if omi < len(ns):
                                    bfile.write(write_bound(obs[mi], m, ns[omi], nf[omi]))
                                    omi += 1
                        bfile.write(')\n\n// **\n')
                        
                    for file in os.listdir(svp['flparams']['ofcpfilebase']):
                        shutil.copy(os.path.join(svp['flparams']['ofcpfilebase'], file), os.path.join(svp['flparams']['offilebase'], st, 'polyMesh'))
                    
                    oftomesh(svp['flparams']['offilebase'], vl, fomats, st, ns, nf)
                    expnode.post_export()

                else:
                    oftomesh(svp['flparams']['offilebase'], vl, fomats, st, ns, nf)
                    expnode.post_export()

        except Exception as e:
           logentry("Netgen error: {}".format(e))
           return {'CANCELLED'}
        return {'FINISHED'}
    
class NODE_OT_Flo_Bound(bpy.types.Operator):
    bl_idname = "node.flovi_bound"
    bl_label = "Boundary export"
    bl_description = "Export Openfoam boundaries"
    bl_register = True
    bl_undo = False
    
    def execute(self, context):
        scene = context.scene
        svp = scene.vi_params
        dobs = [o for o in bpy.data.objects if o.vi_params.vi_type == '2']
        gobs = [o for o in bpy.data.objects if o.vi_params.vi_type == '3']
        obs = dobs + gobs
        boundnode = context.node
        meshnode = boundnode.inputs['Mesh in'].links[0].from_node
        casenode = meshnode.inputs['Case in'].links[0].from_node
        fvvarwrite(scene, obs, casenode)

        if boundnode.pv:
            subprocess.Popen(shlex.split('foamExec paraFoam -builtin -case {}'.format(svp['flparams']['offilebase'])))
        else:
            open("{}".format(os.path.join(svp['flparams']['offilebase'], 'project.foam')), "w")
            
        boundnode.post_export()
        return {'FINISHED'}
    
class NODE_OT_Flo_Sim(bpy.types.Operator):
    bl_idname = "node.flovi_sim"
    bl_label = "FloVi simulation"
    bl_description = "Solve an OpenFOAM case"
    bl_register = True
    bl_undo = True
    
    def modal(self, context, event):
        svp = context.scene.vi_params

        if self.run.poll() is None and self.kivyrun.poll() is None:
            with open(self.fpfile, 'r') as fpfile:
                lines = fpfile.readlines()[::-1]
                residict = {}

                for line in lines:
                    if 'Solving for' in line:
                        residict[line.split()[3][:-1]] = abs(float(line.split()[7][:-1]))
                    elif len(line.split()) and line.split()[0] == 'Time' and line.split()[1] == '=':
                        residict[line.split()[0]] = float(line.split()[2])
                    if len(residict) == len(self.residuals) + 1:
                        break

                for var in residict:
                    if var in ('epsilon', 'omega', 'k'):
                        residict[var] = residict[var] - self.econvergence if residict[var] > self.econvergence else 0
                    elif var == 'p':
                        residict[var] = residict[var] - self.pconvergence if residict[var] > self.pconvergence else 0
                    else:
                        residict[var] = residict[var] - self.convergence if residict[var] > self.convergence else 0
 
                if residict:
                    self.pfile.check("\n".join(['{0[0]} {0[1]}'.format(i) for i in residict.items()]))
            return {'PASS_THROUGH'}
        
        elif self.run.poll() is None and self.kivyrun.poll() is not None:
            self.run.kill()
            return {'CANCELLED'}
        else:
            self.kivyrun.kill()
            
            if self.processes > 1:
                Popen(shlex.split("foamExec reconstructPar -case {}".format(svp['flparams']['offilebase']))).wait()
                
            Popen(shlex.split("foamExec postProcess -func writeCellCentres -case {}".format(svp['flparams']['offilebase']))).wait()
            
            if self.pv:
                Popen(shlex.split("foamExec paraFoam -builtin -case {}".format(svp['flparams']['offilebase'])))
            else:
                open("{}".format(os.path.join(svp['flparams']['offilebase'], 'project.foam')), "w")

            reslists = []

            if os.path.isdir(os.path.join(svp['flparams']['offilebase'], 'postProcessing', 'probes', '0')):
                probed = os.path.join(svp['flparams']['offilebase'], 'postProcessing', 'probes', '0')
                resdict = {'p': 'Pressure', 'U': 'Speed', 'T': 'Temperature', 'Ux': 'X velocity', 'Uy': 'Y velocity', 'Uz': 'Z velocity'}

                for f in os.listdir(probed):
                    if f in ('p', 'T'):
                        with open(os.path.join(probed, f), 'r') as resfile:
                            res = []
                            for l in resfile.readlines():
                                if l and l[0] != '#':
                                    res.append(l.split())
                            resarray = array(res)
                            resarray = transpose(resarray)
                                    
                        for ri, r in enumerate(resarray[1:]):
                            reslists.append([str(context.scene.frame_current), 'Zone', svp['flparams']['probes'][ri], resdict[f], ' '.join(['{:5f}'.format(float(res)) for res in r])])

                reslists.append([str(context.scene.frame_current), 'Time', '', 'Steps', ' '.join(['{}'.format(f) for f in resarray[0]])])
                self.simnode['reslists'] = reslists
                self.simnode['frames'] = [context.scene.frame_current]
            
            self.simnode.post_sim()
            return {'FINISHED'}

    def invoke(self, context, event):
        wm = context.window_manager
        scene = context.scene
        svp = scene.vi_params
        self.simnode = context.node
        
        for root, dirs, files in os.walk(svp['flparams']['offilebase']):
            for d in dirs:
                try:
                    if float(d) != svp['flparams']['st'] or 'processor' in d:
                        shutil.rmtree(os.path.join(root, d))
                except:
                    pass
               
        (self.convergence, self.econvergence, self.pconvergence, self.residuals, self.processes)  = (svp['flparams']['uresid'], svp['flparams']['keoresid'], svp['flparams']['presid'], svp['flparams']['residuals'], self.simnode.processes)
        self.fpfile = os.path.join(svp['viparams']['newdir'], 'floviprogress')
        self.pfile = fvprogressfile(svp['viparams']['newdir'])
        self.kivyrun = fvprogressbar(os.path.join(svp['viparams']['newdir'], 'viprogress'), str(self.residuals))
        self.pv = self.simnode.pv
        
        with open(self.fpfile, 'w') as fvprogress:
            if self.processes > 1:
                with open(os.path.join(svp['flparams']['ofsfilebase'], 'decomposeParDict'), 'w') as fvdcpfile:
                    fvdcpfile.write(fvdcpwrite(self.processes))

                Popen(shlex.split("foamExec decomposePar -force -case {}".format(svp['flparams']['offilebase']))).wait()
                self.run = Popen(shlex.split('mpirun --oversubscribe -np {} foamExec {} -parallel -case {}'.format(self.processes, svp['flparams']['solver'], svp['flparams']['offilebase'])), stdout = fvprogress)
            else:
                self.run = Popen(shlex.split('{} {} {} {}'.format('foamExec', svp['flparams']['solver'], "-case", svp['flparams']['offilebase'])), stdout = fvprogress)

        self._timer = wm.event_timer_add(5, window = context.window)
        wm.modal_handler_add(self)        
        return {'RUNNING_MODAL'}
    
    def terminate(self, scene):
        self.run.kill()

#        fvsolvew(casenode, solver)
#        with open(os.path.join(svp['flparams']['ofcpfilebase'], 'blockMeshDict'), 'w') as bmfile:
#            bmfile.write(fvbmwrite(bmos[0], expnode))

#        call(("blockMesh", "-case", "{}".format(scene['flparams']['offilebase'])))
#        fvblbmgen(bmos[0].data.materials, open(os.path.join(scene['flparams']['ofcpfilebase'], 'faces'), 'r'), open(os.path.join(scene['flparams']['ofcpfilebase'], 'points'), 'r'), open(os.path.join(scene['flparams']['ofcpfilebase'], 'boundary'), 'r'), 'blockMesh')
#        expnode.export()

                    # with open(os.path.join(svp['flparams']['offilebase'], 'constant', 'polyMesh', 'boundary'), 'r') as bfile:
                    #     nf = []
                    #     ns = []
                    #     for line in bfile.readlines():
                    #         if 'nFaces' in line:
                    #             nf.append(int(line.split()[1].strip(';')))
                    #         if 'startFace' in line:
                    #             ns.append(int(line.split()[1].strip(';')))
            
                    # with open(os.path.join(svp['flparams']['offilebase'], 'constant', 'polyMesh', 'boundary'), 'w') as bfile:
                    #     bfile.write(ofheader) 
                    #     cl = 'polyBoundaryMesh' if expnode.bl_label == 'FloVi NetGen' else 'BoundaryMesh'
                    #     loc = 'constant/polyMesh' if expnode.bl_label == 'FloVi NetGen' else 'Mesh'
                    #     bfile.write(write_ffile(cl, loc, 'boundary'))
                    #     bfile.write('// **\n\n{}\n(\n'.format(len(ns)))
                    #     omi = 0
                        
                    #     for mi, mats in enumerate(omats):
                    #         for m in mats:
                    #             if omi < len(ns):
                    #                 bfile.write(write_bound(obs[mi], m, ns[omi], nf[omi]))
                    #                 omi += 1
                                
                    #     bfile.write(')\n\n// **\n')
                        
                    # if expnode.poly:
                    #     Popen(shlex.split('foamExec polyDualMesh -case {} -noFields -overwrite {}'.format(svp['flparams']['offilebase'], expnode.yang))).wait()
                    #     Popen(shlex.split('foamExec combinePatchFaces -overwrite -case {} {}'.format(svp['flparams']['offilebase'], expnode.yang))).wait()
                        
                    #     for file in os.listdir(os.path.join(svp['flparams']['offilebase'], st, 'polyMesh')):
                    #         if os.path.isfile(os.path.join(svp['flparams']['offilebase'], st, 'polyMesh', 'file')):
                    #             shutil.copy(os.path.join(svp['flparams']['offilebase'], st, 'polyMesh', file), os.path.join(svp['flparams']['ofcpfilebase']))
                        
                    #     with open(os.path.join(svp['flparams']['offilebase'], '0', 'polyMesh', 'boundary'), 'r') as bfile:
                    #         nf = []
                    #         ns = []
                    #         for line in bfile.readlines():
                    #             if 'nFaces' in line:
                    #                 nf.append(int(line.split()[1].strip(';')))
                    #             if 'startFace' in line:
                    #                 ns.append(int(line.split()[1].strip(';')))



                    # mesh = bpy.data.meshes.new("mesh") 
                    # vcoords = []
                    # findices = []
                    # fi = []
                    # fn = 0
                    # prevline = ''    
                        
                    # with open(os.path.join(svp['flparams']['offilebase'], st, 'polyMesh', 'points'), 'r') as mfile:
                    #     for li, line in enumerate(mfile.readlines()):
                    #         if '(' in line and ')' in line:
                    #             vcoords.append([float(x) for x in line.split('(')[1].split(')')[0].split()])
                                
                    # with open(os.path.join(svp['flparams']['offilebase'], st, 'polyMesh', 'faces'), 'r') as mfile:
                    #     for li, line in enumerate(mfile.readlines()):
                    #         if line:
                    #             if fn:
                    #                 try:
                    #                     fi.append(int(line))
                    #                     fn -= 1
                    #                 except:
                    #                     pass
                                        
                    #             if not fn and fi:
                    #                 findices.append(fi) 
                    #                 fi = []  
                    #                 fn = 0
                    #             elif '(' in line and ')' in line:
                    #                 findices.append([int(x) for x in line.split('(')[1].split(')')[0].split()])
                                    
                    #             else:
                    #                 try:
                    #                     if prevline == '\n' and int(line) < 100:
                    #                         fn = int(line)
                    #                 except:
                    #                     fn = 0
                    #         prevline = line
                    
                    # mesh.from_pydata(vcoords, [], findices)
                    # o = bpy.data.objects.new('Mesh', mesh)
                    # bpy.context.view_layer.active_layer_collection.collection.objects.link(o)
                    # selobj(vl, o)
                    
                    # for mat in fomats:
                    #     bpy.ops.object.material_slot_add()
                    #     o.material_slots[-1].material = mat
                        
                    # bpy.ops.object.material_slot_add()
                    # o.material_slots[-1].material = bpy.data.materials.new("Volume") if 'Volume' not in [m.name for m in bpy.data.materials] else bpy.data.materials["Volume"]
                    
                    # for face in o.data.polygons:
                    #     mi = 0
                    #     for ni, n in enumerate(ns):
                    #         if face.index >= n and face.index <= n + nf[ni]:
                    #             face.material_index = ni
                    #             mi = 1
                    #     if not mi:
                    #         face.material_index = len(ns)

                # with open(os.path.join(svp['flparams']['offilebase'], 'constant', 'polyMesh', 'boundary'), 'r') as bfile:
                #     nf = []
                #     ns = []
                #     for line in bfile.readlines():
                #         if 'nFaces' in line:
                #             nf.append(int(line.split()[1].strip(';')))
                #         if 'startFace' in line:
                #             ns.append(int(line.split()[1].strip(';')))
        
                # with open(os.path.join(svp['flparams']['offilebase'], 'constant', 'polyMesh', 'boundary'), 'w') as bfile:
                #     bfile.write(ofheader) 
                #     cl = 'polyBoundaryMesh' if expnode.bl_label == 'FloVi NetGen' else 'BoundaryMesh'
                #     loc = 'constant/polyMesh' if expnode.bl_label == 'FloVi NetGen' else 'Mesh'
                #     bfile.write(write_ffile(cl, loc, 'boundary'))
                #     bfile.write('// **\n\n{}\n(\n'.format(len(ns)))
                #     omi = 0
                    
                #     for mi, mats in enumerate(omats):
                #         for m in mats:
                #             if omi < len(ns):
                #                 bfile.write(write_bound(obs[mi], m, ns[omi], nf[omi]))
                #                 omi += 1
                            
                #     bfile.write(')\n\n// **\n')