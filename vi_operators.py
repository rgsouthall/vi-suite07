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

import bpy, datetime, mathutils, os, bmesh, shutil, sys, math, shlex
import numpy
from numpy import arange, histogram, array, int8, float16
import bpy_extras.io_utils as io_utils
from subprocess import Popen, PIPE
from collections import OrderedDict
from datetime import datetime as dt
from math import cos, sin, pi, ceil, tan, radians
from time import sleep
from mathutils import Euler, Vector
from xml.dom.minidom import parse, parseString
#from gpu_extras.batch import batch_for_shader
#from bpy_extras import view3d_utils

#from multiprocessing import Pool
from .livi_export import radgexport, createoconv, createradfile
from .livi_calc  import li_calc
from .envi_export import enpolymatexport, pregeo
from .envi_mat import envi_materials, envi_constructions

from .vi_func import selobj, joinobj, solarPosition, viparams, wind_compass, livisimacc

#from .flovi_func import fvcdwrite, fvbmwrite, fvblbmgen, fvvarwrite, fvsolwrite, fvschwrite, fvtppwrite, fvraswrite, fvshmwrite, fvmqwrite, fvsfewrite, fvobjwrite, fvdcpwrite
from .vi_func import ret_plt, logentry, rettree, cmap

from .vi_func import windnum, wind_rose, create_coll, retobjs, progressfile, progressbar
from .vi_func import chunks, clearlayers, clearscene, clearfiles, objmode
from .livi_func import retpmap, radpoints
#from .vi_display import wr_legend, wr_scatter, wr_table, wr_disp
#from .envi_func import processf, retenvires, envizres, envilres, recalculate_text
from .vi_chart import chart_disp

#from .vi_display import wr_legend, results_bar, wr_table, wr_scatter, svf_legend
#from .vi_display import li_display, linumdisplay, ss_legend, ss_scatter, livi_legend#, spnumdisplay, en_air, wr_legend, wr_disp, wr_scatter, wr_table, ss_disp, ss_legend, svf_disp, svf_legend, basic_legend, basic_table, basic_disp, ss_scatter, en_disp, en_pdisp, en_scatter, en_table, en_barchart, comp_table, comp_disp, leed_scatter, cbdm_disp, cbdm_scatter, envals, bsdf, bsdf_disp#, en_barchart, li3D_legend

try:    
    import matplotlib
    matplotlib.use('qt5agg', force = True)
    import matplotlib.cm as mcm
    import matplotlib.colors as mcolors
    mp = 1    
except Exception as e:
#    logentry('Matplotlib error: {}'.format(e))   
    print(e)
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
        sbinvals = arange(0,int(ceil(max(cws))),2)
        dbinvals = arange(-11.25,372.25,22.5)
        dfreq = histogram(awd, bins=dbinvals)[0]
        adfreq = histogram(cwd, bins=dbinvals)[0]
        dfreq[0] = dfreq[0] + dfreq[-1]
        dfreq = dfreq[:-1]
        
        if simnode.wrtype == '0':
            ax.bar(vawd, vaws, bins=sbinvals, normed=True, opening=0.8, edgecolor='white', cmap=mcm.get_cmap(svp.vi_leg_col))
        elif simnode.wrtype == '1':
            ax.box(vawd, vaws, bins=sbinvals, normed=True, cmap=mcm.get_cmap(svp.vi_leg_col))
        elif simnode.wrtype in ('2', '3', '4'):
            ax.contourf(vawd, vaws, bins=sbinvals, normed=True, cmap=mcm.get_cmap(svp.vi_leg_col))
        
        if simnode.max_freq == '1':
            ax.set_rmax(simnode.max_freq_val)
            
        plt.savefig(svp['viparams']['newdir']+'/disp_wind.svg')
        wrme = bpy.data.meshes.new("Wind_rose")   
        wro = bpy.data.objects.new('Wind_rose', wrme) 
        
        if wro.name not in wrcoll.objects:
            wrcoll.objects.link(wro)
            if wro.name in scene.collection.objects:
                scene.collection.objects.unlink(wro)
                
        selobj(context.view_layer, wro)       
        (wro, scale) = wind_rose(wro, simnode['maxres'], svp['viparams']['newdir']+'/disp_wind.svg', simnode.wrtype, mcolors)
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
                 
                    for chunk in chunks(gpoints, int(svp['viparams']['nproc']) * 200):
                        for gp in chunk:
                            pointres = array([(0, 1)[shadtree.ray_cast(posis[g], direc)[3] == None] for direc in valdirecs], dtype = int8)
                            gp[shadres] = (100*(numpy.sum(pointres)/lvaldirecs)).astype(int8)
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

        svp.vi_leg_max, svp.vi_leg_min = 100, 0

        if kivyrun.poll() is None:
            kivyrun.kill()
        
        scene.frame_start, scene.frame_end = svp['liparams']['fs'], svp['liparams']['fe']
        svp['viparams']['vidisp'] = 'svf'
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

class NODE_OT_Li_Geo(bpy.types.Operator):
    bl_idname = "node.ligexport"
    bl_label = "LiVi geometry export"

    def invoke(self, context, event):
        scene = context.scene
        svp = scene.vi_params
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
            if self.kivyrun.poll() is None:
                self.kivyrun.kill() 
            
            if self.o.vi_params.li_bsdf_proxy:
                self.pkgbsdfrun = Popen(shlex.split("pkgBSDF -s {}".format(os.path.join(context.scene.vi_params['viparams']['newdir'], 'bsdfs', '{}.xml'.format(self.mat.name)))), stdin = PIPE, stdout = PIPE)
                self.mat.vi_params['radentry'] = ''.join([line.decode() for line in self.pkgbsdfrun.stdout])
                print(self.mat.vi_params['radentry'])
                
            with open(os.path.join(context.scene.vi_params['viparams']['newdir'], 'bsdfs', '{}.xml'.format(self.mat.name)), 'r') as bsdffile:               
                self.mat.vi_params['bsdf']['xml'] = bsdffile.read()
                bsdf = parseString(self.mat.vi_params['bsdf']['xml'])
                self.mat.vi_params['bsdf']['direcs'] = [path.firstChild.data for path in bsdf.getElementsByTagName('WavelengthDataDirection')]
                self.mat.vi_params['bsdf']['type'] = [path.firstChild.data for path in bsdf.getElementsByTagName('AngleBasisName')][0]
            
            context.scene.vi_params['viparams']['vidisp'] = 'bsdf'
            return {'FINISHED'}
    
    def execute(self, context):
        scene = context.scene
        svp = scene.vi_params
        depsgraph = bpy.context.evaluated_depsgraph_get()
        vl = context.view_layer
        self.o = context.object
        ovp = self.o.vi_params
        ovp.bsdf_running = 1
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
        
        tm = self.o.evaluated_get(depsgraph).to_mesh()
        bm = bmesh.new()    
        bm.from_mesh(tm) 
        self.o.evaluated_get(depsgraph).to_mesh_clear()
        bm.transform(self.o.matrix_world)
        bm.normal_update()
        bsdffaces = [face for face in bm.faces if self.o.material_slots[face.material_index].material.vi_params.radmatmenu == '8']    
        
        if bsdffaces:
            fvec = bsdffaces[0].normal
            mvp['bsdf']['normal'] = '{0[0]:.4f} {0[1]:.4f} {0[2]:.4f}'.format(fvec)
        else:
            self.report({'ERROR'}, '{} does not have a BSDF material associated with any faces'.format(self.o.name))
            return
        
        self.pfile = progressfile(svp['viparams']['newdir'], datetime.datetime.now(), 100)
        self.kivyrun = progressbar(os.path.join(svp['viparams']['newdir'], 'viprogress'), 'BSDF')
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
        mradfile = ''.join([m.vi_params.radmat(scene) for m in self.o.data.materials if m.vi_params.radmatmenu != '8'])                  
        gradfile = radpoints(self.o, [face for face in bm.faces if self.o.material_slots and face.material_index < len(self.o.material_slots) and self.o.material_slots[face.material_index].material.vi_params.radmatmenu != '8'], 0)
        bm.free()  
        bsdfsamp = ovp.li_bsdf_ksamp if ovp.li_bsdf_tensor == ' ' else 2**(int(ovp.li_bsdf_res) * 2) * int(ovp.li_bsdf_tsamp)         
        gbcmd = "genBSDF +geom meter -r '{}' {} {} -c {} {} -n {}".format(ovp.li_bsdf_rcparam,  ovp.li_bsdf_tensor, (ovp.li_bsdf_res, ' ')[ovp.li_bsdf_tensor == ' '], bsdfsamp, ovp.li_bsdf_direc, svp['viparams']['nproc'])

        with open(os.path.join(svp['viparams']['newdir'], 'bsdfs', '{}_mg'.format(self.mat.name)), 'w') as mgfile:
            mgfile.write(mradfile+gradfile)

        with open(os.path.join(svp['viparams']['newdir'], 'bsdfs', '{}_mg'.format(self.mat.name)), 'r') as mgfile: 
            with open(os.path.join(svp['viparams']['newdir'], 'bsdfs', '{}.xml'.format(self.mat.name)), 'w') as bsdffile:
                self.bsdfrun = Popen(shlex.split(gbcmd), stdin = mgfile, stdout = bsdffile)
                
        vl.objects.active = self.o
        wm = context.window_manager
        self._timer = wm.event_timer_add(1, window = context.window)
        wm.modal_handler_add(self)        
        return {'RUNNING_MODAL'}
        
class MATERIAL_OT_Li_LBSDF(bpy.types.Operator, io_utils.ImportHelper):
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
        
class MATERIAL_OT_Li_SBSDF(bpy.types.Operator):
    bl_idname = "material.save_bsdf"
    bl_label = "Save BSDF"
    bl_description = "Save a BSDF for the current selected object"
    bl_register = True

    filename_ext = ".XML;.xml;"
    filter_glob: bpy.props.StringProperty(default="*.XML;*.xml;", options={'HIDDEN'})
    filepath: bpy.props.StringProperty(subtype='FILE_PATH', options={'HIDDEN', 'SKIP_SAVE'})
    
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
        
class NODE_OT_Li_Pre(bpy.types.Operator, io_utils.ExportHelper):
    bl_idname = "node.radpreview"
    bl_label = "LiVi preview"
    bl_description = "Prevew the scene with Radiance"
    bl_register = True
    bl_undo = False
    
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
                self.pfile = progressfile(svp['viparams']['newdir'], datetime.datetime.now(), 100)
                self.kivyrun = progressbar(os.path.join(svp['viparams']['newdir'], 'viprogress'), 'Photon Map')
                amentry, pportentry, cpentry, cpfileentry = retpmap(self.simnode, frame, scene)
                open('{}.pmapmon'.format(svp['viparams']['filebase']), 'w')
                pmcmd = 'mkpmap -t 20 -e {1}.pmapmon -fo+ -bv+ -apD 0.001 {0} -apg {1}-{2}.gpm {3} {4} {5} {1}-{2}.oct'.format(pportentry, svp['viparams']['filebase'], frame, self.simnode.pmapgno, cpentry, amentry)
                logentry('Photon map command: {}'.format(pmcmd))
                pmrun = Popen(pmcmd.split(), stderr = PIPE, stdout = PIPE)

                while pmrun.poll() is None:   
                    sleep(10)
                    with open('{}.pmapmon'.format(svp['viparams']['filebase']), 'r') as vip:
                        for line in vip.readlines()[::-1]:
                            if '%' in line:
                                curres = float(line.split()[6][:-2])
                                break
                                
                    if self.pfile.check(curres) == 'CANCELLED': 
                        pmrun.kill()                                   
                        return {'CANCELLED'}
                
                if self.kivyrun.poll() is None:
                    self.kivyrun.kill()
                        
                with open('{}.pmapmon'.format(svp['viparams']['filebase']), 'r') as pmapfile:
                    for line in pmapfile.readlines():
                        if line in pmerrdict:
                            logentry(line)
                            self.report({'ERROR'}, pmerrdict[line])
                            return {'CANCELLED'}
                                        
                rvucmd = "rvu -w {11} -ap {8} 50 {9} -n {0} -vv {1:.3f} -vh {2} -vd {3[0]:.3f} {3[1]:.3f} {3[2]:.3f} -vp {4[0]:.3f} {4[1]:.3f} {4[2]:.3f} -vu {10[0]:.3f} {10[1]:.3f} {10[2]:.3f} {5} {6}-{7}.oct".format(svp['viparams']['wnproc'], 
                                 vv, cang, vd, cam.location, self.simnode['radparams'], svp['viparams']['filebase'], scene.frame_current, '{}-{}.gpm'.format(svp['viparams']['filebase'], frame), cpfileentry, cam.matrix_world.to_quaternion()@ mathutils.Vector((0, 1, 0)), ('', '-i')[self.simnode.illu])
                
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
            self.processors = simnode['Processors']
            self.processes = simnode.processes
            self.radparams = simnode['radparams']
            self.viewparams = simnode['viewparams']
            self.pmparams = simnode['pmparams']
            self.pmaps = simnode['pmaps']
            self.pmapgnos = simnode['pmapgnos']
            self.pmapcnos = simnode['pmapcnos']
            self.folder = svp['viparams']['newdir']
            self.fb = svp['viparams']['filebase']
            self.mp = 0 if sys.platform == 'win32' else simnode.mp 
            self.basename = simnode['basename']
            
            for frame in range(self.fs, self.fe + 1):
                createradfile(scene, frame, self, simnode)
                createoconv(scene, frame, self, simnode)
                with open('{}-{}'.format(self.pmfile, frame), 'w'):
                    pass
                
            scene.frame_set(svp['liparams']['fs'])
            self.pmcmds = ['mkpmap -t 10 -e {6} -bv+ +fo -apD 0.001 {0} -apg {1}-{2}.gpm {3} {4} {5} {1}-{2}.oct'.format(self.pmparams[str(frame)]['pportentry'], svp['viparams']['filebase'], frame, self.pmapgnos[str(frame)], self.pmparams[str(frame)]['cpentry'], self.pmparams[str(frame)]['amentry'], '{}-{}'.format(self.pmfile, frame)) for frame in range(self.fs, self.fe + 1)]                   
            self.rppmcmds = [('', ' -ap {} {}'.format('{}-{}.gpm'.format(svp['viparams']['filebase'], frame), self.pmparams[str(frame)]['cpfileentry']))[self.pmaps[frame - self.fs]] for frame in range(self.fs, self.fe + 1)]
            self.rpictcmds = ["rpict -t 10 -e {} ".format(self.rpictfile) + ' '.join(['{0[0]} {0[1]}'.format(i) for i in self.viewparams[str(frame)].items()]) + self.rppmcmds[frame - self.fs] + self.radparams + "{0}-{1}.oct".format(svp['viparams']['filebase'], frame, os.path.join(svp['viparams']['newdir'], 'images', self.basename)) for frame in range(self.fs, self.fe + 1)]
            self.rpiececmds = ["rpiece -t 10 -e {} ".format(self.rpictfile) + ' '.join(['{0[0]} {0[1]}'.format(i) for i in self.viewparams[str(frame)].items()]) + self.rppmcmds[frame - self.fs] + self.radparams + "-o {2}-{1}.hdr {0}-{1}.oct".format(svp['viparams']['filebase'], frame, os.path.join(svp['viparams']['newdir'], 'images', self.basename)) for frame in range(self.fs, self.fe + 1)]
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
        scene = context.scene
        svp = scene.vi_params
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
            fcim = os.path.join(svp['viparams']['newdir'], 'images', '{}-{}.hdr'.format(fcnode['basename'], i + svp['liparams']['fs']))
            ofile = bpy.path.abspath(fcnode.ofile) if os.path.isfile(bpy.path.abspath(fcnode.ofile)) and fcnode.overlay else bpy.path.abspath(im)
                        
            with open(fcim, 'w') as fcfile:
                if sys.platform == 'win32':
                    temp_file = os.path.join(svp['viparams']['newdir'], 'images', 'temp.hdr')
                    
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

class NODE_OT_En_UV(bpy.types.Operator):
    bl_idname = "node.envi_uv"
    bl_label = "EnVi Material U-Value Calculation"

    def execute(self, context):
#        resists = []
        node = context.node
        node.ret_uv()
#        lsock = node.inputs['Outer layer']
#                
#        while lsock.links:
#            resists.append(lsock.links[0].from_node.ret_resist())
#            lsock = lsock.links[0].from_node.inputs['Layer']
#
#        node.uv = '{:.3f}'.format(1/(sum(resists) + 0.12 + 0.08))
        return {'FINISHED'}
    
class NODE_OT_En_Con(bpy.types.Operator, io_utils.ExportHelper):
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
        
        for frame in range(node.fs, node.fe + 1):
            scene.frame_set(frame)
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
        self.simnode.postsim(self, condition)
        
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
    
class NODE_OT_Chart(bpy.types.Operator, io_utils.ExportHelper):
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

# Node utilities from Matalogue    
class TREE_OT_goto_mat(bpy.types.Operator):
    'Show the EnVi nodes for this material'
    bl_idname = 'tree.goto_mat'
    bl_label = 'Go To EnVi Material Group'

    mat: bpy.props.StringProperty(default="")

    def execute(self, context):
        context.space_data.tree_type = 'EnViMatN'
#        context.space_data.node_tree = self.mat.vi_params.envi_nodes
#        context.space_data.shader_type = 'OBJECT'
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
#            slot = dummy.material_slots[0]
#            slot.material = mat

        return {'FINISHED'}
    
class TREE_OT_goto_group(bpy.types.Operator):
    'Show the nodes inside this group'
    bl_idname = 'tree.goto_group'
    bl_label = 'Go To Group'

    tree_type: bpy.props.StringProperty(default="")
    tree: bpy.props.StringProperty(default="")

    def execute(self, context):
#        print(self.tree)
        try:  # Go up one group as many times as possible - error will occur when the top level is reached
            while True:
                bpy.ops.node.tree_path_parent()
        except:
            pass
        context.space_data.tree_type = self.tree_type
        context.space_data.path.append(bpy.data.node_groups[self.tree])
        context.space_data.node_tree = bpy.data.node_groups[self.tree]
        context.space_data.node_tree.use_fake_user = 1
#        print(dir(context.space_data))
        return {'FINISHED'}
    
class NODE_OT_CSV(bpy.types.Operator, io_utils.ExportHelper):
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
        node = context.node
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
            self.filepath = os.path.join(context.scene.vi_params['viparams']['newdir'], context.scene.vi_params['viparams']['filebase'] + '.csv')            
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}    