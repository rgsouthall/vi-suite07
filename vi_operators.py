# # ##### BEGIN GPL LICENSE BLOCK #####
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

import bpy, bpy_extras, datetime, mathutils, os, bmesh, shutil, sys, shlex, itertools, inspect, aud, multiprocessing, threading, gc
from pathlib import Path
import subprocess
import numpy
from numpy import arange, histogram, array, int8, int16, int32, float16, empty, uint8, transpose, where, ndarray, place, zeros, average, float32, float64, concatenate, ones, array2string, square
from numpy import sum as nsum
from numpy import max as nmax
from numpy import mean as nmean
from scipy import signal
from scipy.io import wavfile
from scipy import signal
# from bpy.types import Timer
from bpy_extras.io_utils import ExportHelper, ImportHelper
from subprocess import Popen, PIPE, call
from collections import OrderedDict
from datetime import datetime as dt
from math import cos, sin, pi, ceil, tan, radians, log
from time import sleep
from mathutils import Euler, Vector, Matrix
from mathutils.geometry import distance_point_to_plane
from xml.dom.minidom import parseString
from .livi_export import radgexport, createoconv, createradfile, gen_octree, radpoints
from .envi_export import enpolymatexport, pregeo
from .envi_mat import envi_materials, envi_constructions, envi_embodied, envi_eclasstype
from .envi_func import write_ec, write_ob_ec
from .vi_func import selobj, joinobj, solarPosition, viparams, wind_compass
from .flovi_func import ofheader, fvcdwrite, fvvarwrite, fvsolwrite, fvschwrite, fvtpwrite, fvmtwrite
from .flovi_func import fvdcpwrite, write_ffile, write_bound, fvtppwrite, fvgwrite, fvrpwrite, fvprefwrite, oftomesh, fvmodwrite
from .vi_func import ret_plt, logentry, rettree, cmap, fvprogressfile, cancel_window, qtfvprogress
from .vi_func import windnum, wind_rose, create_coll, create_empty_coll, move_to_coll, retobjs, progressfile, progressbar
from .vi_func import chunks, clearlayers, clearscene, clearfiles, objmode, clear_coll, bm_to_stl, qtprogressbar
from .livi_func import retpmap
from .auvi_func import rir2sti
from .vi_chart import chart_disp, hmchart_disp, ec_pie, wlc_line, com_line
from .vi_dicts import rvuerrdict, pmerrdict, flovi_b_dict
from PySide6.QtWidgets import QApplication
import OpenImageIO
OpenImageIO.attribute("missingcolor", "0,0,0")
from OpenImageIO import ImageInput, ImageBuf

if sys.platform != 'win32':
    if multiprocessing.get_start_method() != 'fork':
        multiprocessing.set_start_method('fork', force=True)

try:
    import netgen
    from netgen import occ
    from netgen.meshing import MeshingParameters, FaceDescriptor, Element2D, Mesh, MeshingStep  # , BoundaryLayerParameters
    from pyngcore import SetNumThreads, TaskManager

except Exception as e:
    print(e)

try:
    import matplotlib
    matplotlib.use('qtagg', force=True)
    import matplotlib.cm as mcm
    import matplotlib.colors as mcolors
    from matplotlib import pyplot as plt
    plt.set_loglevel("error")
    from .windrose import WindroseAxes
    mp = 1
except Exception as e:
    print("No matplotlib: {}".format(e))
    mp = 0

try:
    import pyroomacoustics as pra
    pra_rt = True
    pra.constants.set("num_threads", 8)
    ra = 1
except Exception:
    ra = 0

docker_path = 'docker'

if sys.platform == 'darwin':
    if os.path.isfile(f'{Path.home()}/.docker/bin/docker'):
        docker_path = f'{Path.home()}/.docker/bin/docker'

c_freqs = [125, 250, 500, 1000, 2000, 4000, 8000]
pdll_path = os.path.dirname(bpy.app.binary_path)


class NODE_OT_ASCImport(bpy.types.Operator, ImportHelper):
    bl_idname = "node.ascimport"
    bl_label = "Select ESRI Grid file"
    bl_description = "Select the ESRI Grid file to process"
    filename = ""
    filename_ext = ".asc"
    filter_glob: bpy.props.StringProperty(default="*.asc", options={'HIDDEN'})
    bl_register = True
    bl_undo = False

    def draw(self, context):
        layout = self.layout
        row = layout.row()
        row.label(text="Open an asc file with the file browser", icon='WORLD_DATA')

    def execute(self, context):
        scene = context.scene
        asccoll = create_coll(context, 'Terrain')
        startxs, startys, vlen = [], [], 0
        rp = os.path.dirname(os.path.realpath(self.filepath))
        ascfiles = [self.filepath] if self.node.single else [os.path.join(rp, file) for file in os.listdir(rp) if file.endswith('.asc')]
        obs = []
        hd = {'ncols': 0, 'nrows': 0, 'xllcorner': 0, 'yllcorner': 0, 'cellsize': 0, 'NODATA_value': 0}

        for file in ascfiles:
            basename = file.split(os.sep)[-1].split('.')[0]
            me = bpy.data.meshes.new("{} mesh".format(basename))
            bm = bmesh.new()
            li = 0

            with open(file, 'r') as ascfile:
                lines = ascfile.readlines()

                while len(lines[li].split()) == 2:
                    if lines[li].split()[0] in hd:
                        hd[lines[li].split()[0]] = eval(lines[li].split()[1])
                    li += 1

                vlen = hd['nrows'] * hd['ncols']
                startxs.append(hd['xllcorner'])
                startys.append(hd['yllcorner'])
                x, y = 0, hd['nrows']

                for li, line in enumerate(lines[li:]):
                    for zval in line.split():
                        [bm.verts.new((x * hd['cellsize'], y * hd['cellsize'], float(zval)))]
                        x += 1
                    x = 0
                    y -= 1

            bm.verts.ensure_lookup_table()
            faces = [(i + 1, i, i + hd['ncols'], i + hd['ncols'] + 1) for i in range(0, vlen - hd['ncols']) if (i + 1) % hd['ncols']]
            [bm.faces.new([bm.verts[fv] for fv in face]) for face in faces]

            if self.node.clear_nodata == '1':
                bmesh.ops.delete(bm, geom=[v for v in bm.verts if v.co[2] == hd['NODATA_value']], context=1)

            elif self.node.clear_nodata == '0':
                for v in bm.verts:
                    if v.co[2] == hd['NODATA_value']:
                        v.co[2] = 0

            bm.to_mesh(me)
            bm.free()
            ob = bpy.data.objects.new(basename, me)

            if ob.name not in asccoll.objects:
                asccoll.objects.link(ob)
                if ob.name in scene.collection.objects:
                    scene.collection.objects.unlink(ob)

            obs.append(ob)

        minstartx, minstarty = min(startxs), min(startys)

        for o, ob in enumerate(obs):
            ob.location = (startxs[o] - minstartx, startys[o] - minstarty, 0)

        return {'FINISHED'}

    def invoke(self, context, event):
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
        context.view_layer.layer_collection.children[wrcoll.name].exclude = 0
        plt = ret_plt()

        if viparams(self, scene):
            return {'CANCELLED'}

        if not plt:
            self.report({'ERROR'}, "There is something wrong with your matplotlib installation")
            return {'FINISHED'}

        simnode.export()
        locnode = simnode.inputs['Location in'].links[0].from_node
        svp['viparams']['resnode'], svp['viparams']['restree'] = simnode.name, simnode.id_data.name
        svp['viparams']['vidisp'], svp.vi_display = 'wr', 0
        svp['viparams']['visimcontext'] = 'Wind'
        rl = locnode['reslists']
        cdoys = [float(c) for c in [r[4].split() for r in rl if r[0] == '0' and r[1] == 'Time' and r[2] == 'Time' and r[3] == 'DOS'][0]]
        cwd = [float(c) for c in [r[4].split() for r in rl if r[0] == '0' and r[1] == 'Climate' and r[2] == 'Exterior' and r[3] == 'Wind Direction (deg)'][0]]

        if simnode.temp:
            cd = [float(c) for c in [r[4].split() for r in rl if r[0] == '0' and r[1] == 'Climate' and r[2] == 'Exterior' and r[3] == 'Temperature (degC)'][0]]
        else:
            cd = [float(c) for c in [r[4].split() for r in rl if r[0] == '0' and r[1] == 'Climate' and r[2] == 'Exterior' and r[3] == 'Wind Speed (m/s)'][0]]

        doys = list(range(simnode.sdoy, simnode.edoy + 1)) if simnode.edoy > simnode.sdoy else list(range(1, simnode.edoy + 1)) + list(range(simnode.sdoy, 366))
        awd = array([wd for di, wd in enumerate(cwd) if cdoys[di] in doys])
        ad = array([d for di, d in enumerate(cd) if cdoys[di] in doys])
        validdata = where(awd > 0) if max(cwd) == 360 else where(awd > -1)
        vawd = awd[validdata]
        vad = ad[validdata]
        simnode['maxres'], simnode['minres'], simnode['avres'] = max(cd), min(cd), sum(cd) / len(cd)
        sbinvals = arange(0, int(ceil(max(vad))), 2)
        dbinvals = arange(-11.25, 372.25, 22.5)
        dfreq = histogram(awd, bins=dbinvals)[0]
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
            ax.set_rmax(100 * nmax(dfreq) / len(awd) + 0.5)

        plt.savefig(svp['viparams']['newdir'] + '/disp_wind.svg')
        wrme = bpy.data.meshes.new("Wind_rose")
        wro = bpy.data.objects.new('Wind Rose', wrme)

        if wro.name not in wrcoll.objects:
            wrcoll.objects.link(wro)
            if wro.name in scene.collection.objects:
                scene.collection.objects.unlink(wro)

        selobj(context.view_layer, wro)

        (wro, scale) = wind_rose(wro, (simnode['maxres'], simnode.max_freq_val)[simnode.max_freq == '1'], svp['viparams']['newdir'] + '/disp_wind.svg', simnode.wrtype, mcolors)

        wro = joinobj(context.view_layer, wro)
        ovp = wro.vi_params
        ovp['maxres'], ovp['minres'], ovp['avres'], ovp['nbins'], ovp['VIType'] = max(ad), min(ad), sum(ad) / len(ad), len(sbinvals), 'Wind_Plane'
        simnode['maxfreq'] = 100 * nmax(dfreq) / len(awd)
        windnum((100 * nmax(dfreq) / len(awd) + 0.5, simnode.max_freq_val)[simnode.max_freq == '1'], (0, 0, 0), scale, wind_compass((0, 0, 0), scale, wro, wro.data.materials['wr-000000']))
        plt.close()
        ovp['table'] = array([["", 'Minimum', 'Average', 'Maximum'],
                             [('Speed (m/s)', 'Temperature (C)')[simnode.temp], ovp['minres'], '{:.1f}'.format(ovp['avres']), ovp['maxres']],
                             ['Direction (\u00B0)', min(awd), '{:.1f}'.format(sum(awd) / len(awd)), max(awd)]])
        ovp['d'] = ad.reshape(len(doys), 24).T.tolist()
        ovp['wd'] = awd.reshape(len(doys), 24).T.tolist()
        ovp['days'] = array(doys, dtype=float)
        ovp['hours'] = arange(0, 24, dtype=float)
        ovp['maxfreq'] = 100 * nmax(dfreq) / len(awd)
        simnode['nbins'] = len(sbinvals)
        simnode['d'] = array(cd).reshape(365, 24).T.tolist()
        simnode['wd'] = array(cwd).reshape(365, 24).T.tolist()
        simnode['days'] = arange(1, 366, dtype=float)
        simnode['hours'] = arange(0, 24, dtype=float)
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
        context.evaluated_depsgraph_get()

        if viparams(self, scene):
            return {'CANCELLED'}

        shadobs = retobjs('livig')

        if not shadobs:
            self.report({'ERROR'}, "No shading objects with a material attached.")
            return {'CANCELLED'}

        simnode = context.node
        svp['viparams']['restree'] = simnode.id_data.name
        clearscene(context, self)

        for o in scene.objects:
            o.vi_params.vi_type_string = ''

        calcobs = retobjs('ssc')

        if not calcobs:
            self.report({'ERROR'}, "No objects have a light sensor material attached.")
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
            alts = [(rrow + 0.5) * 90 / (2 * 7 + 0.5) for rrow in range(0, 15)]
            azis = (60, 60, 60, 60, 48, 48, 48, 48, 36, 36, 24, 24, 12, 12, 1)

        elif simnode.skypatches == '2':
            alts = [(rrow + 0.5) * 90 / (4 * 7 + 0.5) for rrow in range(0, 29)]
            azis = (120, 120, 120, 120, 120, 120, 120, 120, 96, 96, 96, 96, 96, 96, 96, 96, 72, 72, 72, 72, 48, 48, 48, 48, 24, 24, 24, 24, 1)

        for a, azi in enumerate(azis):
            for az in arange(0, 360, 360 / azi):
                x.append(sin(az * pi / 180) * cos(alts[a] * pi / 180))
                y.append(cos(az * pi / 180) * cos(alts[a] * pi / 180))
                z.append(sin(alts[a] * pi / 180))

        valdirecs = [v for v in zip(x, y, z)]
        lvaldirecs = len(valdirecs)
        calcsteps = len(frange) * sum(len([f for f in o.data.polygons if o.data.materials[f.material_index].vi_params.mattype == '1']) for o in calcobs)
        curres, reslists = 0, []
        pfile = progressfile(svp['viparams']['newdir'], datetime.datetime.now(), calcsteps)
        pb = qtprogressbar(os.path.join(svp['viparams']['newdir'], 'viprogress'), pdll_path, 'Sky View')

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
                shadtree = rettree(scene, [ob for ob in shadobs if ob.visible_get()], ('', '1')[simnode.signore])
                shadres = geom.layers.float['svf{}'.format(frame)]

                if gpoints:
                    posis = [gp.calc_center_median() + gp.normal.normalized() * simnode.offset for gp in gpoints] if simnode.cpoint == '0' else [gp.co + gp.normal.normalized() * simnode.offset for gp in gpoints]

                    for chunk in chunks(gpoints, int(svp['viparams']['nproc']) * 200):
                        for gp in chunk:
                            pointres = array([(0, 1)[not shadtree.ray_cast(posis[g], direc)[3]] for direc in valdirecs], dtype=int8)
                            gp[shadres] = (100 * (nsum(pointres) / lvaldirecs)).astype(int8)
                            g += 1

                        curres += len(chunk)

                        if pfile.check(curres) == 'CANCELLED':
                            return {'CANCELLED'}

                    shadres = [gp[shadres] for gp in gpoints]
                    ovp['omin']['svf{}'.format(frame)], ovp['omax']['svf{}'.format(frame)], ovp['oave']['svf{}'.format(frame)] = min(shadres), max(shadres), sum(shadres) / len(shadres)
                    reslists.append([str(frame), 'Zone spatial', o.name, 'X', ' '.join(['{:.3f}'.format(p[0]) for p in posis])])
                    reslists.append([str(frame), 'Zone spatial', o.name, 'Y', ' '.join(['{:.3f}'.format(p[1]) for p in posis])])
                    reslists.append([str(frame), 'Zone spatial', o.name, 'Z', ' '.join(['{:.3f}'.format(p[2]) for p in posis])])
                    reslists.append([str(frame), 'Zone spatial', o.name, 'SVF', ' '.join(['{:.3f}'.format(sr) for sr in shadres])])
                    avres.append(ovp['oave']['svf{}'.format(frame)])
                    minres.append(ovp['omin']['svf{}'.format(frame)])
                    maxres.append(ovp['omax']['svf{}'.format(frame)])

            reslists.append(['All', 'Frames', 'Frames', 'Frames', ' '.join(['{}'.format(f) for f in frange])])
            reslists.append(['All', 'Zone spatial', o.name, 'Minimum', ' '.join(['{:.3f}'.format(mr) for mr in minres])])
            reslists.append(['All', 'Zone spatial', o.name, 'Average', ' '.join(['{:.3f}'.format(mr) for mr in avres])])
            reslists.append(['All', 'Zone spatial', o.name, 'Maximum', ' '.join(['{:.3f}'.format(mr) for mr in maxres])])
            bm.transform(o.matrix_world.inverted())
            bm.to_mesh(o.data)
            bm.free()
            o.vi_params.vi_type_string = 'LiVi Calc'

        svp.vi_leg_max, svp.vi_leg_min = 100, 0

        if pb.poll() is None:
            pb.kill()

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
        bpy.context.evaluated_depsgraph_get()
        svp.vi_display = 0

        if viparams(self, scene):
            return {'CANCELLED'}

        shadobs = retobjs('livig')

        if not shadobs:
            self.report({'ERROR'}, "No shading objects or none with a material attached.")
            return {'CANCELLED'}

        simnode = context.node
        svp['viparams']['restree'] = simnode.id_data.name
        clearscene(context, self)

        for o in scene.objects:
            o.vi_params.vi_type_string = ''

        calcobs = retobjs('ssc')

        if not calcobs:
            self.report({'ERROR'}, "No objects have a light sensor material attached.")
            return {'CANCELLED'}

        svp['viparams']['visimcontext'] = 'Shadow'

        if not svp.get('liparams'):
            svp['liparams'] = {}

        svp['liparams']['cp'], svp['liparams']['unit'], svp['liparams']['type'] = simnode.cpoint, 'Sunlit (% hrs)', 'VI Shadow'
        simnode.preexport()
        (svp['liparams']['fs'], svp['liparams']['fe']) = (scene.frame_current, scene.frame_current) if simnode.animmenu == 'Static' else (simnode.startframe, simnode.endframe)
        cmap(svp)

        if simnode.starthour > simnode.endhour:
            self.report({'ERROR'}, "End hour is before start hour.")
            return {'CANCELLED'}

        svp['viparams']['resnode'], simnode['Animation'] = simnode.name, simnode.animmenu
        (scmaxres, scminres, scavres) = [[x] * (svp['liparams']['fe'] - svp['liparams']['fs'] + 1) for x in (0, 100, 0)]
        frange = range(svp['liparams']['fs'], svp['liparams']['fe'] + 1)
        time = datetime.datetime(2018, simnode.sdate.month, simnode.sdate.day, simnode.starthour)
        y = 2018 if simnode.edoy >= simnode.sdoy else 2019
        endtime = datetime.datetime(y, simnode.edate.month, simnode.edate.day, simnode.endhour)
        interval = datetime.timedelta(hours=1 / simnode.interval)
        times = [time + interval * t for t in range(int((endtime - time) / interval) + simnode.interval) if simnode.starthour <= (time + interval * t).hour <= simnode.endhour]
        sps = array([solarPosition(t.timetuple().tm_yday, t.hour + t.minute / 60, svp.latitude, svp.longitude)[2:] for t in times])
        valmask = array([sp[0] > 0 for sp in sps], dtype=int8)
        direcs = array([(-sin(sp[1]), -cos(sp[1]), tan(sp[0])) for sp in sps])
        valdirecs = [mathutils.Vector((-sin(sp[1]), -cos(sp[1]), tan(sp[0]))) for sp in sps if sp[0] > 0]

        if not valdirecs:
            self.report({'ERROR'}, "No hours specified with sun above the horizon")
            return {'CANCELLED'}

        lvaldirecs = len(valdirecs)
        ilvaldirecs = 1 / lvaldirecs
        calcsteps = len(frange) * sum(len([f for f in o.data.polygons if o.data.materials[f.material_index].vi_params.mattype == '1']) for o in calcobs)
        curres, reslists = 0, []
        pfile = progressfile(svp['viparams']['newdir'], datetime.datetime.now(), calcsteps)
        # kivyrun = progressbar(os.path.join(scene.vi_params['viparams']['newdir'], 'viprogress'), 'Shadow Map')
        pb = qtprogressbar(os.path.join(scene.vi_params['viparams']['newdir'], 'viprogress'), pdll_path, 'Shadow Map')
        logentry(f'Conducting shadow map calculation with {simnode.interval} samples per hour for {int(len(direcs) / simnode.interval)} total hours and {lvaldirecs} available sun hours')

        for frame in frange:
            reslists.append([str(frame), 'Time', 'Time', 'Month', ' '.join([str(t.month) for t in times])])
            reslists.append([str(frame), 'Time', 'Time', 'Day', ' '.join([str(t.day) for t in times])])
            reslists.append([str(frame), 'Time', 'Time', 'Hour', ' '.join([str(t.hour) for t in times])])
            reslists.append([str(frame), 'Time', 'Time', 'DOS', ' '.join([str(t.timetuple().tm_yday - times[0].timetuple().tm_yday) for t in times])])
            reslists.append([str(frame), 'Time', 'Time', 'DOY', ' '.join([str(t.timetuple().tm_yday) for t in times])])

        for oi, o in enumerate(calcobs):
            ovp = o.vi_params

            for k in [k for k in ovp.keys()]:
                del ovp[k]

            if any([s < 0 for s in o.scale]):
                logentry('Negative scaling on calculation object {}. Results may not be as expected'.format(o.name))
                self.report({'WARNING'}, 'Negative scaling on calculation object {}. Results may not be as expected'.format(o.name))

            ovp['omin'], ovp['omax'], ovp['oave'] = {}, {}, {}

            if simnode.sdoy <= simnode.edoy:
                ovp['days'] = arange(simnode.sdoy, simnode.edoy + 1, dtype=float)
            else:
                ovp['days'] = arange(simnode.sdoy, simnode.edoy + 1, dtype=float)

            ovp['hours'] = arange(simnode.starthour, simnode.endhour + 1, 1 / simnode.interval, dtype=float)
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
                shadtree = rettree(scene, [ob for ob in shadobs if ob.visible_get()], ('', '1')[simnode.signore])
                shadres = geom.layers.float['sm{}'.format(frame)]

                if gpoints:
                    posis = [gp.calc_center_median() + gp.normal.normalized() * simnode.offset for gp in gpoints] if simnode.cpoint == '0' else [gp.co + gp.normal.normalized() * simnode.offset for gp in gpoints]
                    allpoints = zeros((len(gpoints), len(direcs)), dtype=int8)

                    for chunk in chunks(gpoints, int(svp['viparams']['nproc']) * 200):
                        for gp in chunk:
                            pointres = array([(0, 1)[not shadtree.ray_cast(posis[g], direc)[3]] and gp.normal.normalized().dot(direc) > 0 for direc in valdirecs], dtype=int8)
                            place(allpoints[g], valmask == 1, pointres)
                            gp[shadres] = (100 * (nsum(pointres) * ilvaldirecs)).astype(float16)
                            g += 1

                        curres += len(chunk)

                        if pfile.check(curres) == 'CANCELLED':
                            return {'CANCELLED'}

                    ap = average(allpoints, axis=0)
                    shadres = [gp[shadres] for gp in gpoints]
                    hsr = array(100 * ap)
                    ovp['ss{}'.format(frame)] = hsr.reshape(len(ovp['days']), len(ovp['hours'])).T.tolist()
                    reslists.append([str(frame), 'Zone temporal', o.name, 'Sunlit %', ' '.join([str(ss) for ss in hsr])])
                    ovp['omin']['sm{}'.format(frame)] = min(shadres)
                    ovp['omax']['sm{}'.format(frame)] = max(shadres)
                    ovp['oave']['sm{}'.format(frame)] = sum(shadres) / len(shadres)
                    reslists.append([str(frame), 'Zone spatial', o.name, 'X', ' '.join(['{:.3f}'.format(p[0]) for p in posis])])
                    reslists.append([str(frame), 'Zone spatial', o.name, 'Y', ' '.join(['{:.3f}'.format(p[1]) for p in posis])])
                    reslists.append([str(frame), 'Zone spatial', o.name, 'Z', ' '.join(['{:.3f}'.format(p[2]) for p in posis])])
                    reslists.append([str(frame), 'Zone spatial', o.name, 'Sunlit %', ' '.join(['{:.3f}'.format(sr) for sr in shadres])])
                    avres.append(ovp['oave']['sm{}'.format(frame)])
                    minres.append(ovp['omin']['sm{}'.format(frame)])
                    maxres.append(ovp['omax']['sm{}'.format(frame)])

            reslists.append(['All', 'Frames', 'Frames', 'Frames', ' '.join(['{}'.format(f) for f in frange])])
            reslists.append(['All', 'Zone spatial', o.name, 'Min. sunlit %', ' '.join(['{:.3f}'.format(mr) for mr in minres])])
            reslists.append(['All', 'Zone spatial', o.name, 'Ave. sunlit %', ' '.join(['{:.3f}'.format(mr) for mr in avres])])
            reslists.append(['All', 'Zone spatial', o.name, 'Max. sunlit %', ' '.join(['{:.3f}'.format(mr) for mr in maxres])])

            bm.transform(o.matrix_world.inverted())
            bm.to_mesh(o.data)
            bm.free()
            o.vi_params.vi_type_string = 'LiVi Calc'

        svp.vi_leg_max, svp.vi_leg_min = 100, 0

        if pb.poll() is None:
            pb.kill()

        scene.frame_start, scene.frame_end = svp['liparams']['fs'], svp['liparams']['fe']
        simnode['reslists'] = reslists
        simnode['frames'] = [f for f in frange]
        simnode.postexport(scene)
        svp['viparams']['vidisp'] = 'ss'
        return {'FINISHED'}


class OBJECT_OT_EcS(bpy.types.Operator):
    bl_idname = "object.ec_save"
    bl_label = "Embodied material save"

    def execute(self, context):
        '''ID, Quantity, unit, density, weight, ec per unit, ec per kg'''
        ob = context.object if context.object else context.collection
        ovp = ob.vi_params
        envi_ecs = envi_embodied()
        envi_ecs.update()

        if ovp.ec_unit == 'kg':
            weight = f'{ovp.ec_amount:.4f}'
            eckg = ovp.ec_du / ovp.ec_amount
            ecdu = ovp.ec_du
        elif ovp.ec_unit == 'tonnes':
            weight = f'{ovp.ec_amount * 1000:.4f}'
            eckg = ovp.ec_du / ovp.ec_amount * 1000
            ecdu = ovp.ec_du
        elif ovp.ec_unit == 'm2':
            weight = f'{ovp.ec_weight:.4f}'
            eckg = ovp.ec_du / ovp.ec_weight
            ecdu = ovp.ec_du
        elif ovp.ec_unit == 'm3':
            weight = f'{ovp.ec_amount * ovp.ec_density:.4f}'
            eckg = ovp.ec_du / (ovp.ec_amount * ovp.ec_density)
            ecdu = ovp.ec_du
        else:
            weight = f'{ovp.ec_weight:.4f}'
            ecdu = ovp.ec_du
            eckg = ovp.ec_du / ovp.ec_weight

        ec_dict = envi_ecs.get_dat()

        if ovp.ec_class not in ec_dict and ovp.ec_class.upper() in [ec_class.upper() for ec_class in ec_dict]:
            self.report({'ERROR'}, 'A class with this spelling but a different case already exists. Pick a different class name')
            return {'CANCELLED'}

        elif ovp.ec_class not in ec_dict:
            ec_dict[ovp.ec_class] = {}
            ec_dict[ovp.ec_class][ovp.ec_type] = {}
            ec_dict[ovp.ec_class][ovp.ec_type][ovp.ec_name] = {"id": ovp.ec_id,
                                                               "quantity": '{:.4f}'.format(ovp.ec_amount),
                                                               "unit": ovp.ec_unit,
                                                               "density": '{:.4f}'.format(ovp.ec_density),
                                                               "weight": weight,
                                                               "ecdu": '{:.4f}'.format(ecdu),
                                                               "eckg": '{:.4f}'.format(eckg),
                                                               "modules": ovp.ec_mod}

        elif ovp.ec_type not in ec_dict[ovp.ec_class]:
            ec_dict[ovp.ec_class][ovp.ec_type] = {}
            ec_dict[ovp.ec_class][ovp.ec_type][ovp.ec_name] = {"id": ovp.ec_id,
                                                               "quantity": '{:.4f}'.format(ovp.ec_amount),
                                                               "unit": ovp.ec_unit,
                                                               "density": '{:.4f}'.format(ovp.ec_density),
                                                               "weight": weight,
                                                               "ecdu": '{:.4f}'.format(ecdu),
                                                               "eckg": '{:.4f}'.format(eckg),
                                                               "modules": ovp.ec_mod}
        else:
            ec_dict[ovp.ec_class][ovp.ec_type][ovp.ec_name] = {"id": ovp.ec_id,
                                                               "quantity": '{:.4f}'.format(ovp.ec_amount),
                                                               "unit": ovp.ec_unit,
                                                               "density": '{:.4f}'.format(ovp.ec_density),
                                                               "weight": weight,
                                                               "ecdu": '{:.4f}'.format(ecdu),
                                                               "eckg": '{:.4f}'.format(eckg),
                                                               "modules": ovp.ec_mod}

        # if ovp.ec_class not in ec_dict and ovp.ec_class.upper() in [ec_class.upper() for ec_class in ec_dict]:
        #     self.report({'ERROR'}, 'A class with this spelling but a different case already exists. Pick a different class name')
        #     return{'CANCELLED'}

        envi_ecs.set_dat(ec_dict)
        envi_ecs.ec_save()
        ovp.ee.update()
        envi_eclasstype(ovp, context)
        ovp.embodiedclass = ovp.ec_class
        ovp.embodiedtype = ovp.ec_type
        ovp.embodiedname = ovp.ec_name
        return {'FINISHED'}


class OBJECT_OT_EcE(bpy.types.Operator):
    bl_idname = "object.ec_edit"
    bl_label = "Embodied material edit"

    def execute(self, context):
        ob = context.object if context.object else context.collection
        ovp = ob.vi_params
        envi_ecs = envi_embodied()
        envi_ecs.update()
        ec_dict = envi_ecs.get_dat()
        ob_dict = ec_dict[ovp.embodiedclass][ovp.embodiedtype][ovp.embodiedmat]
        ovp.ec_type = ovp.embodiedtype
        ovp.ec_class = ovp.embodiedclass
        ovp.ec_id = ob_dict['id']
        ovp.ec_amount = float(ob_dict['quantity'])
        ovp.ec_density = float(ob_dict['density'])
        ovp.ec_mod = ob_dict['modules']
        ovp.ec_unit = ob_dict['unit']
        ovp.ec_name = ovp.embodiedmat
        ovp.ec_du = float(ob_dict['ecdu'])
        ovp.embodiedclass = 'Custom'
        return {'FINISHED'}


class NODE_OT_EcE(bpy.types.Operator):
    bl_idname = "node.ec_edit"
    bl_label = "Embodied material edit"

    def execute(self, context):
        ob = context.object if context.object else context.collection
        node = context.node
        # ovp = ob.vi_params
        envi_ecs = envi_embodied()
        envi_ecs.update()
        ec_dict = envi_ecs.get_dat()
        node_dict = ec_dict[node.embodiedclass][node.embodiedtype][node.embodiedmat]
        node.ec_class = node.embodiedclass
        node.ec_type = node.embodiedtype
        node.ec_id = node_dict['id']

        try:
            node.ec_amount = float(node_dict['quantity'])
        except Exception:
            node.ec_amount = 0
        try:
            node.ec_density = float(node_dict['density'])
        except Exception:
            node.ec_density = 0
        try:
            node.ec_mod = node_dict['modules']
        except Exception:
            node.ec_mod = 'A1-A3'
        try:
            node.ec_unit = node_dict['unit']
        except Exception:
            node.ec_unit = 'kg'

        node.ec_name = node.embodiedmat
        node.ec_du = float(node_dict['ecdu'])
        node.embodiedclass = 'Custom'
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
        # clearfiles(svp['liparams']['objfilebase'])
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
        self.svp.vi_views = 1
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

    def draw(self, context):
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
            elif self.nodeprop == 'wavname':
                self.node.wavname = self.filepath

        if " " in self.filepath:
            self.report({'ERROR'}, "There is a space either in the filename or its directory location. Remove this space and retry opening the file.")
        return {'FINISHED'}

    def invoke(self, context, event):
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
    filepath: bpy.props.StringProperty(subtype='FILE_PATH', options={'HIDDEN', 'SKIP_SAVE'})
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
                return {'PASS_THROUGH'}
        else:
            self.o.vi_params.bsdf_running = 0
            filepath = os.path.join(context.scene.vi_params['viparams']['newdir'], 'bsdfs', '{}.xml'.format(self.mat.name))

            if self.pb.poll() is None:
                self.pb.kill()

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
        self.pb = qtprogressbar(os.path.join(svp['viparams']['newdir'], 'viprogress'), pdll_path, 'BSDF')
        zvec, yvec = mvp.li_bsdf_up, mathutils.Vector((0, 1, 0))
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
            logentry('Transform unneccesary: {}'.format(e))

#        mvp['bsdf']['up'] = '{0[0]:.4f} {0[1]:.4f} {0[2]:.4f}'.format(mvp.li_bsdf_up)
        vposis = list(zip(*[v.co[:] for v in bm.verts]))
        (maxx, maxy, maxz) = [max(p) for p in vposis]
        (minx, miny, minz) = [min(p) for p in vposis]
        bsdftrans = mathutils.Matrix.Translation(mathutils.Vector((-(maxx + minx) / 2, -(maxy + miny) / 2, -maxz)))
        bm.transform(bsdftrans)
        mradfile = ''.join([m.vi_params.radmat(scene) for m in self.o.data.materials if m.vi_params.radmatmenu != '8'])
        gradfile = radpoints(self.o, [face for face in bm.faces if self.o.material_slots and face.material_index < len(self.o.material_slots) and self.o.material_slots[face.material_index].material.vi_params.radmatmenu != '8'], 0)
        bm.free()
        bsdfsamp = ovp.li_bsdf_ksamp if ovp.li_bsdf_tensor == ' ' else 2**(int(ovp.li_bsdf_res) * 2) * int(ovp.li_bsdf_tsamp)
#        gbcmd = "genBSDF -geom {} -r '{}' {} {} -c {} {} -n {}".format(ovp.li_bsdf_dimen,  ovp.li_bsdf_rcparam,  ovp.li_bsdf_tensor, (ovp.li_bsdf_res, ' ')[ovp.li_bsdf_tensor == ' '], bsdfsamp, ovp.li_bsdf_direc, svp['viparams']['nproc'])
        # Adding MGF geometry does not work (black inner face)
        gbcmd = "genBSDF +geom {} -r '{}' {} {} -c {} {} -n {}".format(ovp.li_bsdf_dimen, ovp.li_bsdf_rcparam, ovp.li_bsdf_tensor, (ovp.li_bsdf_res, ' ')[ovp.li_bsdf_tensor == ' '], bsdfsamp, ovp.li_bsdf_direc, svp['viparams']['nproc'])
        logentry('genBSDF running with command: {}'.format(gbcmd))

        with open(os.path.join(svp['viparams']['newdir'], 'bsdfs', '{}_mg'.format(self.mat.name)), 'w') as mgfile:
            mgfile.write(mradfile + gradfile)

        with open(os.path.join(svp['viparams']['newdir'], 'bsdfs', '{}_mg'.format(self.mat.name)), 'r') as mgfile:
            with open(os.path.join(svp['viparams']['newdir'], 'bsdfs', '{}.xml'.format(self.mat.name)), 'w') as bsdffile:
                self.bsdfrun = Popen(shlex.split(gbcmd), stdin=mgfile, stdout=bsdffile)

        mvp['bsdf']['type'] = 'LBNL/Klems Full' if ovp.li_bsdf_tensor == ' ' else 'Tensor'
        vl.objects.active = self.o
        wm = context.window_manager
        self._timer = wm.event_timer_add(1, window=context.window)
        wm.modal_handler_add(self)
        ovp.bsdf_running = 1
        return {'RUNNING_MODAL'}


class MATERIAL_OT_Li_LBSDF(bpy.types.Operator, ImportHelper):
    bl_idname = "material.load_bsdf"
    bl_label = "Select BSDF file"
    filename_ext = ".XML;.xml;"
    filter_glob: bpy.props.StringProperty(default="*.XML;*.xml;", options={'HIDDEN'})
    filepath: bpy.props.StringProperty(subtype='FILE_PATH', options={'HIDDEN', 'SKIP_SAVE'})

    def draw(self, context):
        layout = self.layout
        row = layout.row()
        row.label(text="Import BSDF XML file with the file browser", icon='WORLD_DATA')
        row = layout.row()

    def execute(self, context):
        mvp = context.material.vi_params
        mvp['bsdf'] = {}
        if not os.path.isfile(self.filepath):
            self.report({'ERROR'}, "No valid XMl file selected.")
            return {'CANCELLED'}
        else:
            with open(self.filepath, 'r') as bsdffile:
                mvp['bsdf']['xml'] = bsdffile.read()
                mvp['bsdf']['filepath'] = self.filepath
                bsdf = parseString(mvp['bsdf']['xml'])
                mvp['bsdf']['direcs'] = [path.firstChild.data for path in bsdf.getElementsByTagName('WavelengthDataDirection')]
                mvp['bsdf']['type'] = [path.firstChild.data for path in bsdf.getElementsByTagName('AngleBasisName')][0]
                mvp['bsdf']['proxied'] = mvp.li_bsdf_proxy_depth
            return {'FINISHED'}

    def invoke(self, context, event):
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}


class MATERIAL_OT_Li_DBSDF(bpy.types.Operator):
    bl_idname = "material.del_bsdf"
    bl_label = "Del BSDF"
    bl_description = "Delete a BSDF for the current selected object"
    bl_register = True
    bl_undo = False

    def execute(self, context):
        del context.material.vi_params['bsdf']
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

    def draw(self, context):
        layout = self.layout
        row = layout.row()
        row.label(text="Save BSDF XML file with the file browser", icon='WORLD_DATA')

    def execute(self, context):
        with open(self.filepath, 'w') as bsdfsave:
            bsdfsave.write(context.material.vi_params['bsdf']['xml'])
        return {'FINISHED'}

    def invoke(self, context, event):
        self.filepath = '{}.xml'.format(context.material.name)
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
            if self.rvurun.poll() is not None:
                self.simnode.run = 0

                for line in self.rvurun.stderr:
                    if b'fatal IO error' not in line and b'events remaining' not in line and b'Broken pipe' not in line and b'explicit kill' not in line:
                        logentry(line)

                    if b'explicit kill' not in line:
                        for rvuerr in rvuerrdict:
                            if rvuerr in line.decode():
                                self.report({'ERROR'}, rvuerrdict[rvuerr])
                                return {'CANCELLED'}
                            else:
                                self.report({'ERROR'}, f'rvu error: {line.decode()}')
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
        elif not all([i.links[0].from_node['Text'] for i in self.simnode.inputs]):
            self.report({'ERROR'}, 'Missing Radiance description. Check Geometry/Context exports')
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
            cang = '180 -vth ' if self.simnode['coptions']['Context'] == 'Basic' and self.simnode['coptions']['Type'] == '1' else cam.data.angle_x * 180 / pi
            vv = 180 if self.simnode['coptions']['Context'] == 'Basic' and self.simnode['coptions']['Type'] == '1' else cam.data.angle_y * 180 / pi
            vd = (0.001, 0, -1 * cam.matrix_world[2][2]) if (round(-1 * cam.matrix_world[0][2], 3), round(-1 * cam.matrix_world[1][2], 3)) == (0.0, 0.0) else [-1 * cam.matrix_world[i][2] for i in range(3)]

            if self.simnode.pmap:
                self.pfile = progressfile(svp['viparams']['newdir'], datetime.datetime.now(), 100)
                self.pb = qtprogressbar(os.path.join(svp['viparams']['newdir'], 'viprogress'), pdll_path, 'Photon Map')
                amentry, pportentry, gpentry, cpentry, gpfileentry, cpfileentry = retpmap(self.simnode, frame, scene)
                open('{}-{}'.format(self.pmfile, frame), 'w')
                pmcmd = 'mkpmap {8} -t 2 -e "{1}" {6} -fo+ -bv{9} -apD 0.1 {0} {3} {4} {5} "{7}-{2}.oct"'.format(pportentry,
                                                                                                                 '{}-{}'.format(self.pmfile, frame),
                                                                                                                 frame, gpentry,
                                                                                                                 cpentry, amentry,
                                                                                                                 ('-n {}'.format(svp['viparams']['wnproc']), '')[sys.platform == 'win32'],
                                                                                                                 svp['viparams']['filebase'], self.simnode.pmapoptions, ('-', '+')[self.simnode.bfv])
                logentry('Photon map command: {}'.format(pmcmd))
                os.chdir(svp['viparams']['newdir'])
                pmrun = Popen(shlex.split(pmcmd), stderr=PIPE, stdout=PIPE)

                for line in pmrun.stderr:
                    self.report({'ERROR'}, f'mkpmap errer: {line.decode()}')

                    if self.pb.poll() is None:
                        self.pb.kill()

                    pmrun.kill()
                    return {'CANCELLED'}

                while pmrun.poll() is None:
                    sleep(2)
                    with open('{}-{}'.format(self.pmfile, frame), 'r') as vip:
                        for line in vip.readlines()[::-1]:
                            if '%' in line:
                                for entry in line.split():
                                    if '%' in entry:
                                        curres = float(entry[:-2]) if float(entry[:-2]) > curres else curres

                    if self.pfile.check(curres) == 'CANCELLED':
                        pmrun.kill()
                        return {'CANCELLED'}

                if self.pb.poll() is None:
                    self.pb.kill()

                with open('{}-{}'.format(self.pmfile, frame), 'r') as pmapfile:
                    for line in pmapfile.readlines():
                        if 'fatal -' in line:
                            for pmerr in pmerrdict:
                                if pmerr in line:
                                    logentry(pmerrdict[pmerr])
                                    self.report({'ERROR'}, pmerrdict[pmerr])
                                    return {'CANCELLED'}

                            logentry(f'Photon map error: {line}')
                            self.report({'ERROR'}, 'Unknown photon map error. Check the VI-Suite log file')
                            return {'CANCELLED'}

                if self.simnode.pmappreview:
                    create_empty_coll(context, 'LiVi Results')
                    gpmbm = bmesh.new()
                    bmesh.ops.create_icosphere(gpmbm, subdivisions=1, radius=0.025, calc_uvs=False)
                    lvo = len(gpmbm.verts)

                    if self.simnode.pmapgno:
                        verts_out, faces_out = [], []

                        for li, line in enumerate(Popen(shlex.split('pmapdump -a -c 0 0 1 {0}-{1}.gpm'.format(svp['viparams']['filebase'], frame)),
                                                        stdout=PIPE, stderr=PIPE).stdout):
                            dl = line.decode().split()
                            matrix = Matrix.Translation(Vector([float(x) for x in dl[:3]]))
                            nbm = gpmbm.copy()
                            bmesh.ops.transform(nbm, matrix=matrix, verts=nbm.verts)
                            verts_out += [v.co.to_tuple() for v in nbm.verts]
                            faces_out += [[j.index + li * len(nbm.verts) for j in i.verts] for i in nbm.faces]
                            nbm.free()

                            if li > self.simnode.pmapgno or li > 100000:
                                break

                        gpm_mesh = bpy.data.meshes.new('gpm_mesh')
                        gpm_mesh.from_pydata(verts_out, [], faces_out)
                        gpmobj = bpy.data.objects.new('GlobalPM', gpm_mesh)
                        gpmobj.vi_params.vi_type_string = 'LiVi Res'
                        scene.collection.objects.link(gpmobj)
                        move_to_coll(bpy.context, 'LiVi Results', gpmobj)

                    if self.simnode.pmapcno:
                        verts_out, faces_out = [], []

                        for li, line in enumerate(Popen(shlex.split('pmapdump -a -c 0 0 1 {0}-{1}.cpm'.format(svp['viparams']['filebase'], frame)),
                                                        stdout=PIPE, stderr=PIPE).stdout):
                            dl = line.decode().split()
                            matrix = Matrix.Translation(Vector([float(x) for x in dl[:3]]))
                            nbm = gpmbm.copy()
                            bmesh.ops.transform(nbm, matrix=matrix, verts=nbm.verts)
                            verts_out += [v.co.to_tuple() for v in nbm.verts]
                            faces_out += [[j.index + li * len(nbm.verts) for j in i.verts] for i in nbm.faces]
                            nbm.free()

                            if li > self.simnode.pmapcno or li > 100000:
                                break

                        cpm_mesh = bpy.data.meshes.new('cpm_mesh')
                        cpm_mesh.from_pydata(verts_out, [], faces_out)
                        cpmobj = bpy.data.objects.new('CausticPM', cpm_mesh)
                        cpmobj.vi_params.vi_type_string = 'LiVi Res'
                        scene.collection.objects.link(cpmobj)
                        move_to_coll(bpy.context, 'LiVi Results', cpmobj)

                    gpmbm.free()

                rvucmd = 'rvu -w {11} {12} {9} -n {0} -vv {1:.3f} -vh {2:.3f} -vd {3[0]:.3f} {3[1]:.3f} {3[2]:.3f} -vp {4[0]:.3f} {4[1]:.3f} {4[2]:.3f} -vu {10[0]:.3f} {10[1]:.3f} {10[2]:.3f} {5} "{6}-{7}.oct"'.format(svp['viparams']['wnproc'], vv, cang, vd, cam.location, self.simnode['rvuparams'], svp['viparams']['filebase'], scene.frame_current, '{}-{}.gpm'.format(svp['viparams']['filebase'], frame),
                cpfileentry, cam.matrix_world.to_quaternion() @ mathutils.Vector((0, 1, 0)), ('', '-i')[self.simnode.illu], gpfileentry)

            else:
                rvucmd = 'rvu -w {9} -n {0} -vv {1:.3f} -vh {2:.3f} -vd {3[0]:.3f} {3[1]:.3f} {3[2]:.3f} -vp {4[0]:.3f} {4[1]:.3f} {4[2]:.3f} -vu {8[0]:.3f} {8[1]:.3f} {8[2]:.3f} {5} "{6}-{7}.oct"'.format(svp['viparams']['wnproc'],
                                 vv, cang, vd, cam.location, self.simnode['rvuparams'], svp['viparams']['filebase'], scene.frame_current, cam.matrix_world.to_quaternion() @ mathutils.Vector((0, 1, 0)), ('', '-i')[self.simnode.illu])

            logentry('Rvu command: {}'.format(rvucmd))
            self.rvurun = Popen(shlex.split(rvucmd), stdout=PIPE, stderr=PIPE)
            context.node.run = 1
            wm = context.window_manager
            self._timer = wm.event_timer_add(1, window=context.window)
            wm.modal_handler_add(self)
            self.simnode.hide = 0
            return {'RUNNING_MODAL'}


class NODE_OT_Li_Sim(bpy.types.Operator):
    bl_idname = "node.livicalc"
    bl_label = "LiVi simulation"
    bl_register = True
    bl_undo = False

    def modal(self, context, event):
        self.simnode.postsim(self.reslists)
        self.report({'INFO'}, "Simulation is finished")
        return {'FINISHED'}

    def invoke(self, context, event):
        scene = context.scene
        vl = context.view_layer
        frame = scene.frame_current
        svp = scene.vi_params
        svp.vi_display = 0
        svp['viparams']['vidisp'] = ''

        if viparams(self, scene):
            return {'CANCELLED'}

        objmode()
        clearscene(context, self)
        self.simnode = context.node

        if not all([i.links[0].from_node['Text'] for i in self.simnode.inputs]):
            self.report({'ERROR'}, 'Missing Radiance description. Check Geometry/Context exports')
            return {'CANCELLED'}

        self.simnode.presim()
        contextdict = {'Basic': 'LiVi Basic', 'CBDM': 'LiVi CBDM'}

        # Set scene parameters
        svp['viparams']['visimcontext'] = contextdict[self.simnode['coptions']['Context']]
        svp['liparams']['fs'] = min((self.simnode['coptions']['fs'], self.simnode['goptions']['fs']))
        svp['liparams']['fe'] = max((self.simnode['coptions']['fe'], self.simnode['goptions']['fe']))
        svp['liparams']['cp'] = self.simnode['goptions']['cp']
        svp['liparams']['unit'] = self.simnode['coptions']['unit']
        svp['liparams']['type'] = self.simnode['coptions']['Type']
        scene.frame_start, scene.frame_end = svp['liparams']['fs'], svp['liparams']['fe']
        self.simnode.sim(scene)
        curres = 0
        rtcmds, rccmds = [], []
        scontext = self.simnode['coptions']['Context']
        subcontext = self.simnode['coptions']['Type']
        patches = self.simnode['coptions']['cbdm_res']
        rh = (146, 578, 0, 2306, 0, 0, 0, 9218).index(patches) + 1
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
                pmfile = os.path.join(svp['viparams']['newdir'], 'pmprogress')
                pfile = progressfile(svp['viparams']['newdir'], datetime.datetime.now(), 100)
                self.pb = qtprogressbar(os.path.join(svp['viparams']['newdir'], 'viprogress'), pdll_path, 'Photon map')
                amentry, pportentry, gpentry, cpentry, gpfileentry, cpfileentry = retpmap(self.simnode, frame, scene)

                if scontext == 'Basic' or (scontext == 'CBDM' and subcontext == '0'):
                    pmcmd = 'mkpmap {7} -t 2 -e "{1}-{3}" -fo+ -bv+ -apD 0.001 {0} {4} {5} {6} "{2}-{3}.oct"'.format(pportentry, pmfile, svp['viparams']['filebase'],
                                                                                                                     frame, gpentry, cpentry, amentry, ('-n {}'.format(svp['viparams']['wnproc']), '')[sys.platform == 'win32'])
                else:
                    pmcmd = 'mkpmap {4} -t 2 -e "{1}-{3}" -fo+ -bv+ -apC "{2}-{3}.copm" {0} "{2}-{3}.oct"'.format(self.simnode.pmapgno, pmfile, svp['viparams']['filebase'],
                                                                                                                  frame, ('-n {}'.format(svp['viparams']['wnproc']), '')[sys.platform == 'win32'])

                logentry('Generating photon map: {}'.format(pmcmd))
                pmrun = Popen(shlex.split(pmcmd), stderr=PIPE, stdout=PIPE)

                for line in pmrun.stderr:
                    self.report({'ERROR'}, f'mkpmap errer: {line.decode()}')

                    if self.pb.poll() is None:
                        self.pb.kill()

                    pmrun.kill()
                    return {'CANCELLED'}

                while pmrun.poll() is None:
                    sleep(10)
                    with open(f'{pmfile}-{frame}', 'r') as vip:
                        for line in vip.readlines()[::-1]:
                            if '%' in line:
                                curres = float(line.split()[6][:-2]) / len(frames)
                                break

                    if pfile.check(curres) == 'CANCELLED':
                        pmrun.kill()
                        return {'CANCELLED'}

                if self.pb.poll() is None:
                    self.pb.kill()

                with open(f'{pmfile}-{frame}', 'r') as pmapfile:
                    for line in pmapfile.readlines():
                        if 'fatal -' in line:
                            for pmerr in pmerrdict:
                                if pmerr in line:
                                    logentry(pmerrdict[pmerr])
                                    self.report({'ERROR'}, pmerrdict[pmerr])
                                    return {'CANCELLED'}

                            logentry(f'Photon map error: {line}')
                            self.report({'ERROR'}, 'Unknown photon map error. Check the VI-Suite log file')
                            return {'CANCELLED'}

            if scontext == 'Basic' or (scontext == 'CBDM' and subcontext == '0'):
                if os.path.isfile("{}-{}.af".format(svp['viparams']['filebase'], frame)):
                    os.remove("{}-{}.af".format(svp['viparams']['filebase'], frame))

                if self.simnode.pmap:
                    rtcmds.append('rtrace -n {0} -w {1} {5} {4} -faa -h -ov -I "{2}-{3}.oct"'.format(svp['viparams']['nproc'], self.simnode['radparams'], svp['viparams']['filebase'], frame, cpfileentry, gpfileentry))
                else:
                    rtcmds.append('rtrace -n {0} -w {1} -faa -h -ov -I "{2}-{3}.oct"'.format(svp['viparams']['nproc'], self.simnode['radparams'], svp['viparams']['filebase'], frame))
            else:
                if self.simnode.pmap:
                    rccmds.append('rcontrib -w -h -I -fo -ap {2}-{3}.copm {0} -n {1} -e MF:{4} -f reinhart.cal -b rbin -bn Nrbins -m sky_glow "{2}-{3}.oct"'.format(self.simnode['radparams'], svp['viparams']['nproc'], svp['viparams']['filebase'], frame, rh))
                else:
                    rccmds.append('rcontrib -w -h -I -fo {} -n {} -e MF:{} -f reinhart.cal -b rbin -bn Nrbins -m sky_glow "{}-{}.oct"'.format(self.simnode['radparams'], svp['viparams']['nproc'], rh, svp['viparams']['filebase'], frame))

        try:
            tpoints = [o.vi_params['rtpnum'] for o in bpy.data.objects if o.vi_params.vi_type_string == 'LiVi Calc']
        except Exception as e:
            self.report({'ERROR'}, 'Re-export the LiVi geometry: {}'.format(e))
            return {'CANCELLED'}

        calcsteps = sum(tpoints) * len(frames)
        pfile = progressfile(svp['viparams']['newdir'], datetime.datetime.now(), calcsteps)
        self.pb = qtprogressbar(os.path.join(svp['viparams']['newdir'], 'viprogress'), pdll_path, 'Lighting')
        self.reslists = []
        obs = [o for o in bpy.data.objects if o.vi_params.vi_type_string == 'LiVi Calc']

        for oi, o in enumerate(obs):
            ovp = o.vi_params
            curres = sum(tpoints[:oi] * len(frames))
            selobj(vl, o)
            ovp['omax'], ovp['omin'], ovp['oave'] = {}, {}, {}

            if scontext == 'Basic':
                bccout = ovp.basiccalcapply(scene, frames, rtcmds, self.simnode, curres, pfile)
                if bccout == 'CANCELLED':
                    if self.pb.poll() is None:
                        self.pb.kill()
                    return {'CANCELLED'}
                else:
                    self.reslists += bccout

            elif scontext == 'CBDM' and subcontext == '0':
                lhout = ovp.lhcalcapply(scene, frames, rtcmds, self.simnode, curres, pfile)

                if lhout == 'CANCELLED':
                    if self.pb.poll() is None:
                        self.pb.kill()
                    return {'CANCELLED'}
                else:
                    self.reslists += lhout

            elif scontext == 'CBDM' and subcontext in ('1', '2'):
                times = [datetime.datetime.strptime(time, "%d/%m/%y %H:%M:%S") for time in self.simnode['coptions']['times']]

                for f, frame in enumerate(frames):
                    self.reslists.append([str(frame), 'Time', 'Time', 'Month', ' '.join([str(t.month) for t in times])])
                    self.reslists.append([str(frame), 'Time', 'Time', 'Day', ' '.join([str(t.day) for t in times])])
                    self.reslists.append([str(frame), 'Time', 'Time', 'Hour', ' '.join([str(t.hour) for t in times])])
                    self.reslists.append([str(frame), 'Time', 'Time', 'DOS', ' '.join([str(t.timetuple().tm_yday - times[0].timetuple().tm_yday) for t in times])])
                    self.reslists.append([str(frame), 'Time', 'Time', 'DOY', ' '.join([str(t.timetuple().tm_yday) for t in times])])
                    self.reslists.append([str(frame), 'Climate', 'Exterior', 'Daylight', self.simnode['coptions']['dl_hours']])

                cbdmout = ovp.udidacalcapply(scene, frames, rccmds, self.simnode, curres, pfile)

                if cbdmout == 'CANCELLED':
                    if self.pb.poll() is None:
                        self.pb.kill()
                    return {'CANCELLED'}
                else:
                    self.reslists += cbdmout

            elif scontext == 'CBDM' and subcontext == '3':
                cbdmout = ovp.adgpcalcapply(scene, frames, rccmds, self.simnode, curres, pfile)

                if cbdmout == 'CANCELLED':
                    if self.pb.poll() is None:
                        self.pb.kill()
                    return {'CANCELLED'}
                else:
                    self.reslists += cbdmout

        if self.pb.poll() is None:
            self.pb.kill()

        svp['viparams']['resnode'] = self.simnode.name
        svp['viparams']['restree'] = self.simnode.id_data.name
        svp['viparams']['vidisp'] = 'li'
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
                self.pmruns.append(Popen(shlex.split(self.pmcmds[self.p]), stderr=PIPE))
            self.p += 1

        if all([pm.poll() is not None for pm in self.pmruns]) and sum(self.pmaps) == self.p and not self.rpruns:
            self.pmfin = 1

            if len(self.pmruns):
                self.pb.kill()

        if self.pmfin:
            if len(self.rpruns) == 0:
                self.pfile = progressfile(self.folder, datetime.datetime.now(), 100)

                if len(self.pmruns):
                    self.pb = qtprogressbar(os.path.join(self.folder, 'viprogress'), pdll_path, 'Radiance Image')

            if self.mp:
                while self.xindex < self.processes and sum([rp.poll() is None for rp in self.rpruns]) < self.processors and self.frame <= self.fe:
                    echo = Popen(['echo', '{}'.format(self.xindex), '0'], stdout=PIPE)
                    self.rpruns.append(Popen(shlex.split(self.rpiececmds[self.frame - self.fs]), stdin=echo.stdout, stderr=PIPE))

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
                        self.rpruns.append(Popen(shlex.split(self.rpictcmds[len(self.rpruns)]), stdout=imfile, stderr=PIPE))

                    logentry('rpict command: {}'.format(self.rpictcmds[self.frame - self.fs]))

                try:
                    if [rp.poll() for rp in self.rpruns][self.frame - self.fs] is not None:
                        if os.path.join(self.folder, 'images', '{}-{}.hdr'.format(self.basename, self.frame)) not in self.images:
                            self.images.append(os.path.join(self.folder, 'images', '{}-{}.hdr'.format(self.basename, self.frame)))

                            if self.frame < self.fe:
                                self.frame += 1

                except Exception as e:
                    print('Frame passing: {}'.format(e))

        if event.type == 'TIMER':
            f = self.frame if self.frame <= self.fe else self.fe

            if self.pmfin:
                self.percent = 0

                with open('{}-{}'.format(self.pmfile, f), 'r') as vip:
                    vip_lines = vip.readlines()

                    for line in vip_lines:
                        for pmerr in pmerrdict:
                            if pmerr in line:
                                self.report({'ERROR'}, pmerrdict[pmerr])
                                self.pb.kill()
                                self.simnode.run = 0
                                return {'CANCELLED'}

            elif os.path.isfile(f'{self.pmfile}-{f}') and not self.rpruns:
                if sum(self.pmaps):
                    self.percent = 0

                    with open('{}-{}'.format(self.pmfile, f), 'r') as vip:
                        for line in vip.readlines()[::-1]:
                            if '% after' in line:
                                perc = [float(ls[:-2]) for ls in line.split() if '%' in ls][0] / sum(self.pmaps)
                                self.percent = perc if perc > self.percent else self.percent
                                break

            if self.pmfin and self.rpruns and all([rp.poll() is not None for rp in self.rpruns]):
                for line in self.rpruns[0].stderr:
                    logentry('Rpict error: {}'.format(line.decode()))

                    for rvuerr in rvuerrdict:
                        if rvuerr in line.decode():
                            self.report({'ERROR'}, rvuerrdict[rvuerr])
                            self.pb.kill()
                            self.simnode.run = 0
                            return {'CANCELLED'}

                self.imupdate(f)

                if os.path.isfile("{}-{}.hdr_temp".format(os.path.join(self.folder, 'images', self.basename), f)):
                    os.remove("{}-{}.hdr".format(os.path.join(self.folder, 'images', self.basename), f))
                    os.rename("{}-{}.hdr_temp".format(os.path.join(self.folder, 'images', self.basename), f), "{}-{}.hdr".format(os.path.join(self.folder, 'images', self.basename), f))

                return {self.terminate()}

            elif self.pmfin and self.mp:
                if self.percent != 100 * sum([r.poll() is not None for r in self.rpruns]) / (self.processes * self.frames):
                    self.percent = 100 * sum([r.poll() is not None for r in self.rpruns]) / (self.processes * self.frames)
                    self.imupdate(self.fs + int(sum([rp.poll() is not None for rp in self.rpruns]) / self.processes))

            elif self.pmfin and os.path.isfile(self.rpictfile):
                lines = [line for line in open(self.rpictfile, 'r') if '% after' in line][::-1]
                if lines:
                    for lineentry in lines[0].split():
                        if '%' in lineentry and self.percent != (float(lineentry.strip('%')) + (f - self.fs) * 100) / self.frames:
                            newpercent = (float(lineentry.strip('%')) * sum([r.poll() is None for r in self.rpruns]) + 100 * sum([r.poll() is not None for r in self.rpruns])) / self.frames

                            if self.percent != newpercent:
                                self.percent = newpercent
                                self.imupdate(f)

            return {'PASS_THROUGH'}
        else:
            return {'PASS_THROUGH'}

    def imupdate(self, f):
        imfp = "{}-{}.hdr".format(os.path.join(self.folder, 'images', self.basename), f)
        inp = ImageInput.open(imfp)

        if not inp:
            gicmd = f'getinfo "{imfp}"'
            girun = Popen(shlex.split(gicmd), stdout=PIPE, stderr=PIPE)

            for line in girun.stdout:
                lid = line.decode()

                if "VIEW=" in lid:
                    if "-vu" in lid:
                        lsplit = lid.split()
                        vui = lsplit.index("-vu")
                        lsplit[vui + 1] = lsplit[vui + 1][:5]
                        lsplit[vui + 2] = lsplit[vui + 2][:5]
                        lsplit[vui + 3] = lsplit[vui + 3][:5]
                        newline = ' '.join(lsplit)

                        with open(imfp, "r") as ifile:
                            gicmd = f'getinfo -r "{newline}"'

                            with open(imfp + "_temp", "w") as nifile:
                                girun = Popen(shlex.split(gicmd), stdin=ifile, stdout=nifile).wait()

                        inp = ImageInput.open(imfp + "_temp")

        if inp:
            spec = inp.spec()
            rgb = zeros((spec.height, spec.width, 3), float32)

            for y in range(spec.height):
                sl = inp.read_scanline(y, 0, "float32")

                if not isinstance(sl, ndarray):
                    pass
                else:
                    if sl.shape[0] == spec.width:
                        rgb[y] = sl

            rgba = concatenate((rgb[::-1, :, :], ones((spec.height, spec.width, 1), float32)), axis=2)

            if 'livi_preview' not in bpy.data.images:
                im = bpy.data.images.new('livi_preview', spec.width, spec.height, alpha=False, float_buffer=True)
                im.file_format = 'HDR'
            else:
                im = bpy.data.images['livi_preview']
                if im.size[:] != (spec.width, spec.height):
                    im.scale(spec.width, spec.height)

            im.pixels.foreach_set(rgba.flatten())
            inp.close()

    def terminate(self):
        self.pb.kill()

        for pm in self.pmruns:
            if pm.poll() is None:
                pm.kill()

                try:
                    os.system("killall -9 mkpmap")
                except Exception:
                    pass

        for rp in self.rpruns:
            if rp.poll() is None:
                rp.kill()

                try:
                    os.system("killall -9 rpict")
                except Exception:
                    pass

        self.simnode.postsim(self.images)

        if os.path.isfile(self.rpictfile):
            try:
                os.remove(self.rpictfile)
            except Exception:
                pass

        return 'FINISHED'

    def execute(self, context):
        scene = context.scene
        svp = scene.vi_params
        self.xindex, self.p = 0, 0
        self.cam = scene.camera
        simnode = context.node
        self.fs, self.fe = simnode.retframes()
        self.simnode = simnode

        if not all([i.links[0].from_node['Text'] for i in self.simnode.inputs]):
            self.report({'ERROR'}, 'Missing Radiance description. Check Geometry/Context exports')
            return {'CANCELLED'}

        if simnode.camera and bpy.data.cameras.get(simnode.camera.lstrip()):
            self.percent = 0
            self.reslists, self.images = [], []
            self.res = []
            self.rpictfile = os.path.join(svp['viparams']['newdir'], 'rpictprogress')

            if os.path.isfile(self.rpictfile):
                os.remove(self.rpictfile)

            self.pmfile = os.path.join(svp['viparams']['newdir'], 'pmprogress')
            simnode.presim()
            svp['liparams']['fs'], svp['liparams']['fe'] = simnode.retframes()
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
            vps = [' '.join(['{0[0]} {0[1]}'.format(i) for i in self.viewparams[str(frame)].items()]) for frame in range(self.fs, self.fe + 1)]
            vps_vwrays = [' '.join(['{0[0]} {0[1]}'.format(i) for i in self.viewparams[str(frame)].items() if i[0] not in ('-X', '-Y', '-i', '-x', '-y')]) for frame in range(self.fs, self.fe + 1)]
            self.pmcmds = ['mkpmap {7} -t 2 -e "{6}" -bv+ +fo -apD 0.001 {0} -apg "{1}-{2}.gpm" {3} {4} {5} "{1}-{2}.oct"'.format(self.pmparams[str(frame)]['pportentry'],
                           self.fb, frame, self.pmapgnos[str(frame)], self.pmparams[str(frame)]['cpentry'], self.pmparams[str(frame)]['amentry'],
                           '{}-{}'.format(self.pmfile, frame), ('-n {}'.format(svp['viparams']['wnproc']), '')[sys.platform == 'win32']) for frame in range(self.fs, self.fe + 1)]

            self.rppmcmds = [('', ' -ap "{}" {}'.format('{}-{}.gpm'.format(self.fb, frame), self.pmparams[str(frame)]['cpfileentry']))[self.pmaps[frame - self.fs]] for frame in range(self.fs, self.fe + 1)]
            self.rpictcmds = ['rpict -u+ -pa 0 -t 10 -e "{}" '.format(self.rpictfile) + vps[frame - self.fs] + self.rppmcmds[frame - self.fs] + self.radparams + '"{0}-{1}.oct"'.format(self.fb, frame) for frame in range(self.fs, self.fe + 1)]
            self.rpiececmds = ['rpiece -u+ -t 10 -af "{}" -e "{}" '.format('{}-{}.amb'.format(self.fb, frame), self.rpictfile) + vps[frame - self.fs] + self.rppmcmds[frame - self.fs] + self.radparams + '-o "{2}-{1}.hdr" "{0}-{1}.oct"'.format(self.fb, frame, os.path.join(self.folder, 'images', self.basename)) for frame in range(self.fs, self.fe + 1)]

            if simnode.normal or simnode.albedo:
                for frame in range(self.fs, self.fe + 1):
                    # res = (int(simnode.x/self.processes) * self.processes, simnode.y) if self.mp else (simnode.x, simnode.y)
                    res = (simnode.x, simnode.y)
                    vwcmd = 'vwrays -pa 0 -x {0[0]} -y {0[1]} -ff {1}'.format(res, vps_vwrays[frame - self.fs])

                    if simnode.normal:
                        rtcmd = f'rtrace -n {self.processors} -on -ffa -ab 0 -ad 0 -aa 0 -ar 8 -as 0 -dr 0 -lr -1 "{self.fb}-{frame}.oct"'
                        vwrun = Popen(shlex.split(vwcmd), stdout=PIPE, stderr=PIPE)
                        rtrun = Popen(shlex.split(rtcmd), stdin=vwrun.stdout, stdout=PIPE, stderr=PIPE)
                        normdata = [line.decode('utf-8').strip('\n').strip('\r').split('\t')[:3] + ['1'] for line in rtrun.stdout]

                        for ni, nline in enumerate(normdata):
                            if not nline[0]:
                                start_data = ni + 1
                                break

                        try:
                            d_list = normdata[start_data:]
                        except Exception:
                            self.report({'ERROR'}, "Missing octree. Re-export the geometry and context")
                            logentry('ERROR: Missing octree. Re-export the geometry and context')
                            return {'CANCELLED'}

                        d_x = [(float(dl[0]) * 0.5 + 0.5) for dl in d_list]
                        d_y = [(float(dl[1]) * 0.5 + 0.5) for dl in d_list]
                        d_z = [(float(dl[2]) * 0.5 + 0.5) for dl in d_list]
                        d_list = [(d_x[li], d_y[li], d_z[li], 1) for li in range(len(d_x))]
                        d_array = numpy.array(d_list, dtype=numpy.float32).reshape(res[1], res[0], 4)[::-1, :, :]

                        if f'{simnode.camera}-norm-{frame}' in bpy.data.images:
                            nmim = bpy.data.images[f'{simnode.camera}-norm-{frame}']
                            nmim.scale(res[0], res[1])
                        else:
                            nmim = bpy.data.images.new(name=f'{simnode.camera}-norm-{frame}', width=res[0], height=res[1], float_buffer=True, alpha=False, is_data=True)

                        nmim.pixels.foreach_set(d_array.flatten())
                        nmim.pack()

                    if simnode.albedo:
                        albcmd = f'rtrace -n {self.processors} -ov -ffa -ab 0 -ad 0 -av 1.0 1.0 1.0 -aa 0 -ar 8 -as 0 -dr 0 -st 1 -lr -1 "{self.fb}-{frame}.oct"'
                        vw2run = Popen(shlex.split(vwcmd), stdout=PIPE, stderr=PIPE)
                        albrun = Popen(shlex.split(albcmd), stdin=vw2run.stdout, stdout=PIPE, stderr=PIPE)
                        albdata = [line.decode('utf-8').strip('\n').split('\t')[:3] + ['1'] for line in albrun.stdout]

                        for ai, aline in enumerate(albdata):
                            if not aline[0]:
                                start_data = ai + 1
                                break

                        alb_array = numpy.array(albdata[start_data:], dtype=numpy.float32).reshape(res[1], res[0], 4)[::-1, :, :]

                        if f'{simnode.camera}-albedo-{frame}' in bpy.data.images:
                            albim = bpy.data.images[f'{simnode.camera}-albedo-{frame}']
                            albim.scale(res[0], res[1])
                        else:
                            albim = bpy.data.images.new(name=f'{simnode.camera}-albedo-{frame}', width=res[0], height=res[1], float_buffer=True, alpha=True, is_data=True)

                        albim.pixels.foreach_set(alb_array.flatten())
                        albim.pack()

            self.starttime = datetime.datetime.now()
            self.pfile = progressfile(self.folder, datetime.datetime.now(), 100)
            (self.pmfin, flag) = (0, 'Photon Maps') if sum(self.pmaps) else (1, 'Radiance Images')
            self.pb = qtprogressbar(os.path.join(self.folder, 'viprogress'), pdll_path, flag)

            if os.path.isfile("{}-{}.hdr".format(os.path.join(self.folder, 'images', self.basename), self.frame)):
                os.remove("{}-{}.hdr".format(os.path.join(self.folder, 'images', self.basename), self.frame))

            wm = context.window_manager
            self._timer = wm.event_timer_add(2, window=context.window)
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
        vffile = bpy.path.abspath(glnode.vffile)
        vfcmd = f'-vf "{vffile}"' if imnode == glnode else ''
        glnode.presim()

        for i, im in enumerate(imnode['images']):
            with open(im, 'rb') as imtext:
                for line in imtext:
                    if line.decode()[0:2] == '-Y':
                        imy = int(line.decode().split()[1])
                        imx = int(line.decode().split()[3])
                        break

            glfile = os.path.join(svp['viparams']['newdir'], 'images', '{}-{}.hdr'.format(glnode['hdrname'], i + svp['liparams']['fs']))
            gicmd = f'getinfo "{im}"'
            girun = Popen(shlex.split(gicmd), stdout=PIPE, stderr=PIPE)

            if any(['VIEW=' in line.decode() for line in girun.stdout]):
                egcmd = 'evalglare {} -c "{}"'.format(('-u {0[0]} {0[1]} {0[2]}'.format(glnode.gc), '')[glnode.rand], glfile)
            elif not glnode.vffile:
                self.report({'ERROR'}, 'A view file is required but no valid file was supplied')
                return {'CANCELLED'}
            else:
                gicmd2 = f'getinfo "{vffile}"'
                girun2 = Popen(shlex.split(gicmd2), stdout=PIPE, stderr=PIPE)

                if any(['VIEW=' in line.decode() for line in girun2.stdout]):
                    egcmd = 'evalglare {} {} -c "{}"'.format(('-u {0[0]} {0[1]} {0[2]}'.format(glnode.gc), '')[glnode.rand], vfcmd, glfile)
                else:
                    self.report({'ERROR'}, 'View file does not contain a valid view specification')
                    return {'CANCELLED'}

            logentry(f"Running evalglare with command: {egcmd}")

            with open(im, 'r') as hdrfile:
                egrun = Popen(shlex.split(egcmd), stdin=hdrfile, stdout=PIPE, stderr=PIPE)

            if imnode != glnode:
                imc = imnode['coptions']
                time = datetime.datetime(2019, 1, 1, imc['shour'], 0) + datetime.timedelta(imc['sdoy'] - 1) if imc['anim'] == '0' else \
                    datetime.datetime(2019, 1, 1, int(imc['shour']),
                                      int(60 * (imc['shour'] - int(imc['shour'])))) + datetime.timedelta(imc['sdoy'] - 1) + datetime.timedelta(hours=int(imc['interval'] * i),
                                                                                                                                             seconds=int(60 * (imc['interval'] * i - int(imc['interval'] * i))))
                gtime = "{0:0>2d}/{1:0>2d} {2:0>2d}:{3:0>2d}\n".format(time.day, time.month, time.hour, time.minute)
            else:
                time = datetime.datetime(2019, 1, 1, 1)
                gtime = ''

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
                        glaretf.write("{0}dgp: {1:.2f}\ndgi: {2:.2f}\nugr: {3:.2f}\nvcp: {4:.2f}\ncgi: {5:.2f}\nLv: {6:.0f}\n".format(gtime, *res))
                        res.append(res)
                        reslists += [[str(i + svp['liparams']['fs']), 'Camera', 'Camera', 'DGP', '{0[0]}'.format(res)],
                                     [str(i + svp['liparams']['fs']), 'Camera', 'Camera', 'DGI', '{0[1]}'.format(res)],
                                     [str(i + svp['liparams']['fs']), 'Camera', 'Camera', 'UGR', '{0[2]}'.format(res)],
                                     [str(i + svp['liparams']['fs']), 'Camera', 'Camera', 'VCP', '{0[3]}'.format(res)],
                                     [str(i + svp['liparams']['fs']), 'Camera', 'Camera', 'CGI', '{[4]}'.format(res)],
                                     [str(i + svp['liparams']['fs']), 'Camera', 'Camera', 'LV', '{[5]}'.format(res)]]

            with open('{}.temphdr'.format(os.path.join(svp['viparams']['newdir'], 'images', 'glare')), 'w') as temphdr:
                pcondcmd = 'pcond -h+ -u 300 "{}"'.format(glfile)
                Popen(shlex.split(pcondcmd), stdout=temphdr).communicate()

            with open(os.path.join(svp['viparams']['newdir'], 'images', "temp.glare"), "r") as catfile:
                psigncmd = "psign -h {} -cb 0 0 0 -cf 1 1 1".format(int(0.04 * imy))
                psignrun = Popen(shlex.split(psigncmd), stdin=catfile, stdout=PIPE, stderr=PIPE)

            with open(glfile, 'w') as ghdr:
                pcompcmd = 'pcompos "{0}.temphdr" 0 0 - {1} {2}'.format(os.path.join(svp['viparams']['newdir'], 'images', 'glare'), imx, imy * 550 / 800)
                Popen(shlex.split(pcompcmd), stdin=psignrun.stdout, stdout=ghdr).communicate()

            try:
                os.remove(os.path.join(svp['viparams']['newdir'], 'images', 'glare.temphdr'))
                os.remove(os.path.join(svp['viparams']['newdir'], 'images', 'temp.glare'))
            except Exception:
                pass

            if glfile not in [i.filepath for i in bpy.data.images]:
                bpy.data.images.load(glfile)
            else:
                [i.reload() for i in bpy.data.images if i.filepath == glfile][0]

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
        mult = '-m {}'.format(eval('{}{}'.format(179, fcnode.multiplier))) if fcnode.multiplier else ''
        legend = '-l "{}" -lw {} -lh {} {} {} {}'.format(fcnode.unit_name, fcnode.lw, fcnode.lh, lmax, scaling, mult) if fcnode.legend else ''
        contour = ('', '-cl', '-cb', '-cp')[int(fcnode.contour)]
        divisions = '-n {}'.format(fcnode.divisions) if fcnode.divisions != 8 else ''

        for i, im in enumerate(imnode['images']):
            fcim = os.path.join(svp['viparams']['newdir'], 'images', '{}-{}.hdr'.format(fcnode['basename'], i + svp['liparams']['fs']))
            ofile = bpy.path.abspath(fcnode.ofile) if os.path.isfile(bpy.path.abspath(fcnode.ofile)) and fcnode.overlay else bpy.path.abspath(im)

            with open(fcim, 'w') as fcfile:
                if sys.platform == 'win32':
                    temp_file = os.path.join(svp['viparams']['newdir'], 'images', 'temp.hdr')

                    with open(temp_file, 'w') as tfile:
                        pccmd = 'pcond -e {} "{}"'.format(fcnode.disp, os.path.abspath(im))
                        pcrun = Popen(shlex.split(pccmd), stdout=tfile, stderr=PIPE)

                    for line in pcrun.stderr:
                        logentry('Pcond error: {}'.format(line))

                    poverlay = '-p "{}"'.format(os.path.join(svp['viparams']['newdir'], 'images', 'temp.hdr')) if fcnode.contour and fcnode.overlay else ''
                    fccmd = 'falsecolor -i "{}" {} -pal {} {} {} {}'.format(os.path.abspath(im), poverlay, fcnode.coldict[fcnode.colour], legend, contour, divisions)
                    fcrun = Popen(shlex.split(fccmd), stdout=fcfile, stderr=PIPE)

                else:
                    poverlay = '-p <(pcond -e {0} "{1}")' .format(fcnode.disp, ofile) if fcnode.contour and fcnode.overlay else ''
                    fccmd = "bash -c 'falsecolor -i \"{}\" {} -pal {} {} {} {}'".format(bpy.path.abspath(im), poverlay, fcnode.coldict[fcnode.colour], legend, contour, divisions)
                    fcrun = Popen(shlex.split(fccmd), stdout=fcfile, stderr=PIPE)

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
            bpy.ops.node.new_node_tree(type='EnViMatN', name=cm.name)
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
            mvp = mat.vi_params
            m = 0

            if mvp.envi_nodes:
                for o in bpy.data.objects:
                    if mat in [ms.material for ms in o.material_slots]:
                        m = 1

                if not m:
                    mat.delete()

        return {'FINISHED'}


class NODE_OT_En_Geo(bpy.types.Operator):
    bl_idname = "node.engexport"
    bl_label = "EnVi geometry export"
    # bl_context = "scene"

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


class NODE_OT_En_EC(bpy.types.Operator):
    bl_idname = "node.envi_ec"
    bl_label = "EnVi Material Embodied Carbon Calculation"

    def execute(self, context):
        node = context.node

        if node.envi_con_makeup == '1':
            node.ret_ec()

        if node.envi_con_type == 'Window' and node.fclass == '2':
            node.ret_frame_ec()

        return {'FINISHED'}


class NODE_OT_En_Con(bpy.types.Operator, ExportHelper):
    bl_idname = "node.encon"
    bl_label = "Export"
    bl_description = "Export the scene to the EnergyPlus file format"
    bl_register = True
    bl_undo = False

    def invoke(self, context, event):
        scene = context.scene

        if viparams(self, scene):
            return {'CANCELLED'}

        if not bpy.data.collections.get('EnVi Geometry'):
            return {'CANCELLED'}
        else:
            eg_coll = bpy.data.collections['EnVi Geometry']

        svp = scene.vi_params
        svp['viparams']['vidisp'] = ''
        reslists = []
        node = context.node
        frames = range(node.fs, node.fe + 1) if node.animated else [scene.frame_current]
        (svp['enparams']['fs'], svp['enparams']['fe']) = (node.fs, node.fe) if node.animated else (scene.frame_current, scene.frame_current)
        locnode = node.inputs['Location in'].links[0].from_node

        if not os.path.isfile(locnode.weather):
            self.report({'ERROR'}, 'Location node weather file is not valid')
            node.use_custom_color = 1
            return {'CANCELLED'}

        node.preexport(scene)

        for fi, frame in enumerate(frames):
            scene.frame_set(frame)

            # if locnode.outputs['Parameter'].links:
            #     af = bpy.data.texts[locnode.outputs['Parameter'].links[0].to_node.anim_file].as_string()
            #     param = locnode.outputs['Parameter'].links[0].to_node.parameter

            #     for p in locnode.bl_rna.properties:
            #         if p.is_skip_save:
            #             if p.identifier == param:
            #                 for v in locnode['entries']:
            #                     if v[1] == af.split('\n')[fi]:
            #                         try:
            #                             setattr(locnode, param, v[0])
            #                         except Exception as e:
            #                             self.report({'ERROR'}, 'Error in parametric text file: {}'.format(e))
            #                             return {'CANCELLED'}

            shutil.copyfile(locnode.weather, os.path.join(svp['viparams']['newdir'], "in{}.epw".format(frame)))

        if context.active_object and not context.active_object.visible_get():
            if context.active_object.type == 'MESH':
                bpy.ops.object.mode_set(mode='OBJECT')

        error = enpolymatexport(self, eg_coll, node, locnode, envi_materials(), envi_constructions())

        if error:
            return {'CANCELLED'}

        svp['ecparams'] = {'ec_text': 'Scenario, Entity, Entity name, ID, Class, Type, Sub-type, Modules, Volume (m3)/Surface (m2), EC (kgCO2e), EC (kgCO2e/y), EC (kgCO2e/m2), EC (kgCO2e/m2/y)\n'}
        reslists = write_ob_ec(scene, eg_coll, frames, reslists)
        reslists = write_ec(scene, eg_coll, frames, reslists)
        node['reslists'] = reslists
        node.bl_label = node.bl_label[1:] if node.bl_label[0] == '*' else node.bl_label
        node.exported, node.outputs['Context out'].hide = True, False
        node.postexport()
        scene.frame_set(node.fs)
        return {'FINISHED'}


class NODE_OT_En_PVA(bpy.types.Operator):
    bl_idname = "node.pv_area"
    bl_label = "EnVi Material PV area calculation"

    def execute(self, context):
        node = context.node
        try:
            node['area'] = bpy.data.materials[node.id_data.name].vi_params['enparams']['pvarea']
        except Exception:
            node['area'] = 0
        return {'FINISHED'}


class NODE_OT_En_PVS(bpy.types.Operator):
    bl_idname = "node.pv_save"
    bl_label = "EnVi Material PV save"

    def execute(self, context):
        node = context.node
        node.save_e1ddict()
        return {'FINISHED'}


class NODE_OT_En_LayS(bpy.types.Operator):
    bl_idname = "node.lay_save"
    bl_label = "EnVi material save"
    bl_description = "Save layer to the material database"

    def execute(self, context):
        node = context.node
        node.save_laydict()
        return {'FINISHED'}


class NODE_OT_En_ConS(bpy.types.Operator):
    bl_idname = "node.con_save"
    bl_label = "EnVi construction save"
    bl_description = "Save material to the construction database"

    def execute(self, context):
        node = context.node
        node.save_condict()
        return {'FINISHED'}


class NODE_OT_En_EcS(bpy.types.Operator):
    bl_idname = "node.ec_save"
    bl_label = "EnVi embodied material save"

    def execute(self, context):
        node = context.node
        node.save_ecdict()
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
            self.esimruns.append(Popen(self.esimcmds[self.e].split(), stderr=PIPE))
            self.e += 1

        if event.type == 'TIMER':
            for esim in self.esimruns:
                if esim.poll() is None:
                    errtext = esim.stderr.read().decode()

                    if 'Fatal' in errtext:
                        logentry('There is something wrong with the Energyplus installation. Check the message below')
                        logentry('If using EMS a local installation of EnergyPlus is required')
                        logentry(errtext)
                        return {self.terminate('CANCELLED', context)}

            if len(self.esimruns) > 1:
                self.percent = 100 * sum([esim.poll() is not None for esim in self.esimruns]) / self.lenframes

            else:
                try:
                    with open(os.path.join(self.nd, f"{self.resname}{self.frame}out.err"), 'r') as errfile:
                        err_lines = errfile.readlines()

                        if not err_lines:
                            self.report({'ERROR'}, f"Illegal entry in the in{self.frame}.idf file. Check the file in Blender's text editor")
                            return {self.terminate('CANCELLED', context)}

                    with open(os.path.join(self.nd, f"{self.resname}{self.frame}out.eso"), 'r') as resfile:
                        res_lines = [line for line in resfile.readlines()[::-1] if line.split(',')[0] == '2' and len(line.split(',')) == 9]

                        if not res_lines:
                            self.report({'ERROR'}, f"Fatal error reported in the in{self.frame}.idf file. Check the file in Blender's text editor")
                            return {self.terminate('CANCELLED', context)}

                        for resline in res_lines:
                            self.percent = 100 * int(resline.split(',')[1]) / (self.simnode.dedoy - self.simnode.dsdoy)
                            break

                except Exception as e:
                    print(e)
                    logentry('There was an error in the EnVi simulation. Check the error log in the text editor')

            if all([esim.poll() is not None for esim in self.esimruns]) and self.e == self.lenframes:
                for fname in [fname for fname in os.listdir('.') if fname.split(".")[0] == self.simnode.resname]:
                    os.remove(os.path.join(self.nd, fname))

                for f in range(self.frame, self.frame + self.e):
                    nfns = [fname for fname in os.listdir('.') if fname.split(".")[0] == f"{self.resname}{f}out"]

                    for fname in nfns:
                        os.rename(os.path.join(self.nd, fname), os.path.join(self.nd, fname.replace("eplusout", self.simnode.resname)))

                    efilename = f"{self.resname}{f}out.err"

                    if os.path.isfile(os.path.join(self.nd, efilename)):
                        if efilename not in [im.name for im in bpy.data.texts]:
                            bpy.data.texts.load(os.path.join(self.nd, efilename))
                        else:
                            bpy.data.texts[efilename].filepath = os.path.join(self.nd, efilename)

                            with open(os.path.join(self.nd, efilename), 'r') as e_file:
                                bpy.data.texts[efilename].from_string(e_file.read())

                        if '**  Fatal  **' in bpy.data.texts[efilename].as_string():
                            self.report({'ERROR'}, "Fatal error reported in the {} file. Check the file in Blender's text editor".format(efilename))
                            return {self.terminate('CANCELLED', context)}

                        if 'EnergyPlus Terminated--Error(s) Detected' in self.esimruns[f - self.frame].stderr.read().decode() or not [f for f in nfns if f.split(".")[1] == "eso"] or self.simnode.run == 0:
                            if not [f for f in nfns if f.split(".")[1] == "eso"]:
                                errtext = "There is no results file. Check you have selected results outputs and that there are no errors in the .err file in the Blender text editor."
                            else:
                                errtext = "There was an error in the input IDF file in{}.idf. Check the *.err file in Blender's text editor.".format(self.frame)

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
        self.pb = qtprogressbar(os.path.join(svp['viparams']['newdir'], 'viprogress'), pdll_path, 'EnergyPlus Results')
        wm = context.window_manager
        self._timer = wm.event_timer_add(1, window=context.window)
        wm.modal_handler_add(self)
        self.simnode = context.node
        self.connode = self.simnode.inputs[0].links[0].from_node.name
        self.simnode.presim(context)
        self.expand = "-x" if svp['viparams'].get('hvactemplate') else ""
        self.resname = (self.simnode.resname, 'eplus')[self.simnode.resname == '']
        os.chdir(svp['viparams']['newdir'])
        self.esimcmds = ["energyplus {0} -w in{1}.epw -p {2} in{1}.idf".format(self.expand, frame, ('{}{}'.format(self.resname, frame))) for frame in self.frames]
        logentry(f"Running EnergyPlus with command: {self.esimcmds[0]}")
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

        self.pb.kill()

        for es in self.esimruns:
            if es.poll() is None:
                es.kill()

        return condition


class OBJECT_OT_VIGridify(bpy.types.Operator):
    ''''''
    bl_idname = "object.vi_gridify"
    bl_label = "VI Gridify"
    bl_options = {"REGISTER", 'UNDO'}

    rotate: bpy.props.FloatProperty(name='Rotation', default=0, min=0, max=360)
    us: bpy.props.FloatProperty(default=0.6, min=0.01)
    acs: bpy.props.FloatProperty(default=0.6, min=0.01)

    @classmethod
    def poll(cls, context):
        obj = context.active_object
        return (obj and obj.type == 'MESH')

    def execute(self, context):
        o = bpy.context.active_object
        mesh = bmesh.from_edit_mesh(o.data)
        mesh.transform(o.matrix_world)
        mesh.normal_update()
        mesh.faces.ensure_lookup_table()
        mesh.verts.ensure_lookup_table()
        fs = [f for f in mesh.faces[:] if f.select]

        if fs:
            self.upv = fs[0].calc_tangent_edge().copy().normalized()
            self.norm = fs[0].normal.copy()
            self.acv = self.upv.copy()
            eul = Euler(radians(-90) * self.norm, 'XYZ')
            self.acv.rotate(eul)
            self.acv = self.upv.cross(self.norm)
            rotation = Euler(radians(self.rotate) * self.norm, 'XYZ')
            self.upv.rotate(rotation)
            self.acv.rotate(rotation)
            allverts = list(set([entry for sublist in [f.verts[:] for f in fs] for entry in sublist]))
            alledges = list(set([entry for sublist in [f.edges[:] for f in fs] for entry in sublist]))
            vertdots = [Vector.dot(self.upv, vert.co) for vert in allverts]
            vertdots2 = [Vector.dot(self.acv, vert.co) for vert in allverts]
            svpos = allverts[vertdots.index(min(vertdots))].co
            svpos2 = allverts[vertdots2.index(min(vertdots2))].co
            res1, res2, ngs1, ngs2, gs1, gs2 = 1, 1, self.us, self.acs, self.us, self.acs
            gs = allverts + alledges + fs

            while res1:
                res = bmesh.ops.bisect_plane(mesh, geom=gs, dist=0.001, plane_co=svpos + ngs1 * self.upv, plane_no=self.upv, use_snap_center=0, clear_outer=0, clear_inner=0)
                res1 = res['geom_cut']
                dissvs = [v for v in res1 if isinstance(v, bmesh.types.BMVert) and len(v.link_edges) == 2 and v.calc_edge_angle(1) == 0.0]
                bmesh.ops.dissolve_verts(mesh, verts=dissvs)
                gs = mesh.verts[:] + mesh.edges[:] + [v for v in res['geom'] if isinstance(v, bmesh.types.BMFace)]
                ngs1 += gs1

            while res2:
                res = bmesh.ops.bisect_plane(mesh, geom=gs, dist=0.001, plane_co=svpos2 + ngs2 * self.acv, plane_no=self.acv, use_snap_center=0, clear_outer=0, clear_inner=0)
                res2 = res['geom_cut']
                dissvs = [v for v in res2 if isinstance(v, bmesh.types.BMVert) and len(v.link_edges) == 2 and v.calc_edge_angle(1) == 0.0]
                bmesh.ops.dissolve_verts(mesh, verts=dissvs)
                gs = mesh.verts[:] + mesh.edges[:] + [v for v in res['geom'] if isinstance(v, bmesh.types.BMFace)]
                ngs2 += gs2

            mesh.transform(o.matrix_world.inverted())
            bmesh.update_edit_mesh(o.data)
            mesh.free()
            bpy.ops.object.editmode_toggle()
            bpy.ops.object.editmode_toggle()
            return {'FINISHED'}
        else:
            self.report({'ERROR'}, "No faces selected")
            return {'CANCELLED'}


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


class NODE_OT_EC(bpy.types.Operator):
    '''Calculates embodied energy based on object volumes
        Writes to reslists volume, ec (kgco2e), ec/m2 floor area'''
    bl_idname = "node.ec_calc"
    bl_label = "Embodied carbon"
    bl_options = {"REGISTER", 'UNDO'}

    def execute(self, context):
        scene = context.scene

        if viparams(self, scene):
            return {'CANCELLED'}

        envi_ec = envi_embodied()
        envi_ec.update()
        dp = bpy.context.evaluated_depsgraph_get()
        node = context.node
        node.presim(context)
        entity = 'Object' if node.entities == '0' else 'Zone'
        reslists = []
        frames = range(node.startframe, node.endframe + 1) if node.parametric else (context.scene.frame_current, )
        # envi_coll = (bpy.data.collections.get('EnVi Geometry') and (o.name not in bpy.data.collections['EnVi Geometry'].all_objects)) or not bpy.data.collections.get('EnVi Geometry')

        if node.entities == '0':
            obs = [o for o in context.scene.objects if o.type == 'MESH' and (o.vi_params.embodied or o.users_collection[0].vi_params.embodied) and o.visible_get() and \
                                                                            ((bpy.data.collections.get('EnVi Geometry') and (o.name not in bpy.data.collections['EnVi Geometry'].all_objects) or not bpy.data.collections.get('EnVi Geometry')))]

            for frame in frames:
                scene.frame_set(frame)
                ecs, vols = [], []

                for o in obs:
                    ovp = o.vi_params if o.vi_params.embodied else o.users_collection[0].vi_params

                    if all((ovp.embodiedclass, ovp.embodiedtype, ovp.embodiedmat)):
                        if ovp.embodiedclass == 'Custom':
                            logentry(f"Object {o.name} has unsaved EC data")
                            self.report({'ERROR'}, f"Object {o.name} has unsaved EC data")
                            return {'CANCELLED'}

                        ecdict = envi_ec.propdict[ovp.embodiedclass][ovp.embodiedtype][ovp.embodiedmat]
                        vol, mass, area = 0, 0, 0

                        if ecdict['unit'] in ('kg', 'm3', 'm2', 'tonnes') or (ecdict['unit'] == 'each' and ovp.ec_rep == '0'):
                            bm = bmesh.new()
                            bm.from_object(o, dp)
                            bm.transform(o.matrix_world)

                            if node.heal:
                                bmesh.ops.remove_doubles(bm, verts=bm.verts, dist=0.001)
                                bmesh.ops.recalc_face_normals(bm, faces=bm.faces)

                            if all([e.is_manifold for e in bm.edges]) and ecdict['unit'] in ('m3', 'kg'):
                                vol = bm.calc_volume()
                                mass = vol * float(ecdict['density'])
                                ec = float(ecdict['eckg']) * float(ecdict['density']) * vol + ovp.ec_amount_mod * vol / (float(ecdict['weight']) / float(ecdict['density']))  # * node.tyears / float(ovp.ec_life)

                            elif ecdict['unit'] == 'each':
                                ec = float(ecdict['ecdu']) / float(ecdict['quantity'])  # * node.tyears / float(ovp.ec_life)

                            elif ecdict['unit'] == 'm2':
                                if ovp.ec_arep == '0':
                                    area = ovp.ec_ma
                                elif ovp.ec_arep == '1':
                                    area = max([f.calc_area() for f in bm.faces])
                                else:
                                    area = sum([f.calc_area() for f in bm.faces])

                                ec = float(ecdict['ecdu']) * area / float(ecdict['quantity'])  # * node.tyears / float(ovp.ec_life)

                            else:
                                logentry(f"{o.name} is not manifold. Embodied energy metrics have not been exported")
                                self.report({'WARNING'}, f"{o.name} is not manifold. Embodied energy metrics have not been exported")
                                continue

                            bm.free()

                        else:
                            ec = (float(ecdict['ecdu']) + ovp.ec_amount_mod) * ovp.ec_items  # if ovp.ec_rep == '1' else float(ecdict['ecdu']) * node.tyears / float(ovp.ec_life) * vol/(float(ecdict['eckg']) * float(ecdict['density']))

                        ecs.append(ec)
                        vols.append(vol)
                        reslists.append([f'{frame}', 'Embodied carbon', o.name, 'ID', ecdict['id']])
                        reslists.append([f'{frame}', 'Embodied carbon', o.name, 'Declared unit', ecdict['unit']])
                        reslists.append([f'{frame}', 'Embodied carbon', o.name, 'DU quantity', ecdict['quantity']])
                        reslists.append([f'{frame}', 'Embodied carbon', o.name, 'Modules', ecdict['modules']])

                        if vol:
                            reslists.append([f'{frame}', 'Embodied carbon', o.name, 'Object volume (m3)', f'{vol:.4f}'])
                        else:
                            reslists.append([f'{frame}', 'Embodied carbon', o.name, 'Object volume (m3)', 'N/A'])

                        if mass:
                            reslists.append([f'{frame}', 'Embodied carbon', o.name, 'Object density (kg/m3)', f'{float(ecdict["density"]):.4f}'])
                            reslists.append([f'{frame}', 'Embodied carbon', o.name, 'Object mass (kg)', f'{mass:.4f}'])
                        else:
                            reslists.append([f'{frame}', 'Embodied carbon', o.name, 'Object density (kg/m3)', 'N/A'])
                            reslists.append([f'{frame}', 'Embodied carbon', o.name, 'Object mass (kg)', 'N/A'])

                        if area:
                            reslists.append([f'{frame}', 'Embodied carbon', o.name, 'Object area (m2)', f'{area:.4f}'])
                        else:
                            reslists.append([f'{frame}', 'Embodied carbon', o.name, 'Object area (m2)', 'N/A'])

                        reslists.append([f'{frame}', 'Embodied carbon', o.name, 'Object EC(kgCO2e/DU)', f'{float(ecdict["ecdu"]):.4f}'])
                        reslists.append([f'{frame}', 'Embodied carbon', o.name, 'Object EC (kgCO2e)', '{:.4f}'.format(ec)])
                        reslists.append([f'{frame}', 'Embodied carbon', o.name, 'Object EC (kgCO2e/y)', '{:.4f}'.format(ec / ovp.ec_life)])
                        reslists.append([f'{frame}', 'Embodied carbon', o.name, 'Object EC (kgCO2e/m2)', '{:.4f}'.format(ec / node.fa)])
                        reslists.append([f'{frame}', 'Embodied carbon', o.name, 'Object EC (kgCO2e/m2/y)', '{:.4f}'.format(ec / (node.fa * ovp.ec_life))])

                    else:
                        self.report({'WARNING'}, f"Object {o.name} has a missing embodied carbon setting")
                        logentry('Object {} has a missing embodied carbon setting'.format(o.name))

        else:
            cobs = {coll.name: [ob for ob in coll.objects if ob.type == 'MESH' and \
            (ob.vi_params.embodied or coll.vi_params.embodied) and ob.visible_get()] for coll in bpy.data.collections if not coll.hide_viewport and \
            ((bpy.data.collections.get('EnVi Geometry') and coll.name not in bpy.data.collections['EnVi Geometry'].children) or not bpy.data.collections.get('EnVi Geometry'))}

            for frame in frames:
                scene.frame_set(frame)
                tvols = []

                for cob in cobs:
                    ecs, ecys, vols, nmobs = [], [], [], []

                    for ob in cobs[cob]:
                        ovp = ob.vi_params if ob.vi_params.embodied else bpy.data.collections[cob].vi_params

                        if all((ovp.embodiedclass, ovp.embodiedtype, ovp.embodiedmat)):
                            if ovp.embodiedclass == 'Custom':
                                logentry(f"Object {o.name} has unsaved EC data")
                                self.report({'ERROR'}, f"Object {o.name} has unsaved EC data")
                                return {'CANCELLED'}

                            ecdict = envi_ec.propdict[ovp.embodiedclass][ovp.embodiedtype][ovp.embodiedmat]

                            if ecdict['unit'] in ('kg', 'm3', 'm2', 'tonnes') or (ecdict['unit'] == 'each' and ovp.ec_rep == '0'):
                                bm = bmesh.new()
                                bm.from_object(ob, dp)
                                bm.transform(ob.matrix_world)

                                if node.heal:
                                    bmesh.ops.remove_doubles(bm, verts=bm.verts, dist=0.001)
                                    bmesh.ops.recalc_face_normals(bm, faces=bm.faces)

                                if all([e.is_manifold for e in bm.edges]):
                                    vol = bm.calc_volume()
                                    ec = float(ecdict['eckg']) * float(ecdict['density']) * vol + ovp.ec_amount_mod * vol / (float(ecdict['weight']) / float(ecdict['density']))
                                    ecy = ec / ovp.ec_life

                                elif ecdict['unit'] == 'm2':
                                    area = sum([f.calc_area() for f in bm.faces])
                                    ec = float(ecdict['ecdu']) * area / float(ecdict['quantity'])
                                    ecy = ec / ovp.ec_life
                                    vol = 0
                                    mass = 0

                                else:
                                    logentry(f"{ob.name} is not manifold. Embodied energy metrics have not been exported")
                                    self.report({'WARNING'}, f"{ob.name} is not manifold. Embodied energy metrics have not been exported")
                                    nmobs.append(ob)

                                bm.free()

                            else:
                                ec = (float(ecdict['ecdu']) + ovp.ec_amount_mod) * ovp.ec_items
                                ecy = ec / ovp.ec_life
                                vol = 0

                            ecs.append(ec)
                            ecys.append(ecy)
                            vols.append(vol)

                        else:
                            self.report({'WARNING'}, f"Object {ob.name} has a missing embodied carbon setting")
                            logentry('Object {} has a missing embodied carbon setting'.format(ob.name))

                    if cobs[cob] and len(nmobs) != len(cobs[cob]):
                        reslists.append([f'{frame}', 'Embodied carbon', cob, 'Zone EC (kgCO2e)', '{:.3f}'.format(sum(ecs))])
                        reslists.append([f'{frame}', 'Embodied carbon', cob, 'Zone EC (kgCO2e/y)', '{:.3f}'.format(sum(ecys))])
                        # reslists.append([f'{frame}', 'Embodied carbon', cob, 'Zone volume (m3)', '{:.3f}'.format(sum(vols))])
                        reslists.append([f'{frame}', 'Embodied carbon', cob, 'Zone EC (kgCO2e/m2)', '{:.3f}'.format(sum(ecs) / node.fa)])
                        reslists.append([f'{frame}', 'Embodied carbon', cob, 'Zone EC (kgCO2e/m2/y)', '{:.3f}'.format(sum(ecys) / node.fa)])

                    elif cobs[cob] and len(nmobs) == len(cobs[cob]):
                        logentry(f"All objects in collection {cob} are non-manifold. Embodied energy metrics have not been exported")
                        self.report({'WARNING'}, f"All objects in collection {c} are non-manifold. Embodied energy metrics have not been exported")

        if len(frames) > 1:
            obs = [ob.name for ob in obs] if node.entities == '0' else [cob for cob in cobs if cobs[cob]]

            if obs:
                reslists.append(['All', 'Frames', 'Frames', 'Frames', ' '.join(['{}'.format(frame) for frame in frames])])

                for ob in obs:
                    reslists.append(['All', 'Embodied carbon', ob, f'{entity} volume (m3)', ' '.join([ec[4] for ec in reslists if ec[2] == ob and ec[3] == f'{entity} volume (m3)'])])
                    reslists.append(['All', 'Embodied carbon', ob, f'{entity}  EC (kgCO2e)', ' '.join([ec[4] for ec in reslists if ec[2] == ob and ec[3] == f'{entity} EC (kgCO2e)'])])
                    reslists.append(['All', 'Embodied carbon', ob, f'{entity}  EC (kgCO2e/y)', ' '.join([ec[4] for ec in reslists if ec[2] == ob and ec[3] == f'{entity} EC (kgCO2e/y)'])])
                    reslists.append(['All', 'Embodied carbon', ob, f'{entity}  EC (kgCO2e/m2)', ' '.join([ec[4] for ec in reslists if ec[2] == ob and ec[3] == f'{entity} EC (kgCO2e/m2)'])])
                    reslists.append(['All', 'Embodied carbon', ob, f'{entity}  EC (kgCO2e/m2/y)', ' '.join([ec[4] for ec in reslists if ec[2] == ob and ec[3] == f'{entity} EC (kgCO2e/m2/y)'])])

        node['reslists'] = reslists
        node.postsim()
        scene.frame_set(frames[0])
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
        bm = bmesh.new()
        bm.from_object(o, dp)
        bm.transform(o.matrix_world)
        bm.normal_update()

        if all([e.is_manifold for e in bm.edges]):
            envi_ec = envi_embodied()
            vol = bm.calc_volume()
            ovp['ecdict'] = envi_ec.propdict[ovp.embodiedclass][ovp.embodiedtype][ovp.embodiedmat]
            ovp['ecdict']['ec'] = float(ovp['ecdict']['eckg']) * float(ovp['ecdict']['density']) * vol
            bm.free()
        else:
            self.report({'ERROR'}, "You cannot calculate embodied carbon on a non-manifold mesh")
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

        try:
            node.inputs['X-axis'].framemenu
        except Exception as e:
            if node.inputs['X-axis'].framemenu not in zrl[0]:
                self.report({'ERROR'}, f"There are no results in the results file. Check the results.err file in Blender's text editor: {e}")
                return {'CANCELLED'}

        if not mp:
            self.report({'ERROR'}, "Matplotlib cannot be found by the Python installation used by Blender")
            return {'CANCELLED'}

        plt.clf()
        Sdate = dt.fromordinal(dt(year, 1, 1).toordinal() + node['Start'] - 1)  # + datetime.timedelta(hours = node.dsh - 1)
        Edate = dt.fromordinal(dt(year, 1, 1).toordinal() + node['End'] - 1) + datetime.timedelta(hours=23)
        chart_disp(self, plt, node, innodes, Sdate, Edate)
        return {'FINISHED'}


class NODE_OT_HMChart(bpy.types.Operator, ExportHelper):
    bl_idname = "node.hmchart"
    bl_label = "Heatmap"
    bl_description = "Create a 2D heatmap from the results file"
    bl_register = True
    bl_undo = True

    def invoke(self, context, event):
        node = context.node
        node.nodeupdate(context)

        if not mp:
            self.report({'ERROR'}, "Matplotlib cannot be found by the Python installation used by Blender")
            return {'CANCELLED'}

        hmchart_disp(self, plt, node, context.scene.vi_params.vi_leg_col)
        return {'FINISHED'}


class NODE_OT_ECPie(bpy.types.Operator, ExportHelper):
    bl_idname = "node.ec_pie"
    bl_label = "Pie Chart"
    bl_description = "Create a pie chart of embodied carbon"
    bl_register = True
    bl_undo = True

    def invoke(self, context, event):
        node = context.node

        if not mp:
            self.report({'ERROR'}, "Matplotlib cannot be found by the Python installation used by Blender")
            return {'CANCELLED'}

        ec_pie(self, plt, node)
        return {'FINISHED'}


class NODE_OT_WLCLine(bpy.types.Operator, ExportHelper):
    bl_idname = "node.wlc_line"
    bl_label = "Line Chart"
    bl_description = "Create a line chart of whole-life carbon"
    bl_register = True
    bl_undo = True

    def invoke(self, context, event):
        node = context.node

        if not mp:
            self.report({'ERROR'}, "Matplotlib cannot be found by the Python installation used by Blender")
            return {'CANCELLED'}

        wlc_line(self, plt, node)
        return {'FINISHED'}


class NODE_OT_COMLine(bpy.types.Operator, ExportHelper):
    bl_idname = "node.com_line"
    bl_label = "Line Chart"
    bl_description = "Create a line chart of over-heating risk"
    bl_register = True
    bl_undo = True

    def invoke(self, context, event):
        node = context.node

        if not mp:
            self.report({'ERROR'}, "Matplotlib cannot be found by the Python installation used by Blender")
            return {'CANCELLED'}

        com_line(self, plt, node)
        return {'FINISHED'}


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

                if not active_set:
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
        try:
            while True:
                bpy.ops.node.tree_path_parent()
        except Exception:
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

    def draw(self, context):
        layout = self.layout
        row = layout.row()
        row.label(text="Specify the CSV export file with the file browser", icon='WORLD_DATA')

    def execute(self, context):
        node = self.node
        resstring = ''
        resnode = node.inputs['Results in'].links[0].from_node
        rl = resnode['reslists']
        zrl = list(zip(*rl))

        if resnode.bl_idname == 'No_Vi_EC':
            if (len(set(zrl[0])) > 1 and node.animated) or set(zrl[0]) == {'All'}:
                htext, rtext = ',', ''
                onames = [r[2] for r in rl if r[0] != 'All']
                obs = sorted(set(onames), key=onames.index)
                resnames = [r[3] for r in rl if r[0] != 'All']
                res = sorted(set(resnames), key=resnames.index)
                resstring = ''.join(['{} {},'.format(r[2], r[3]) for r in rl if r[0] == 'All']) + '\n'
                metriclist = list(zip(*[r.split() for ri, r in enumerate(zrl[4]) if zrl[0][ri] == 'All']))

                for ml in metriclist:
                    resstring += ''.join(['{},'.format(m) for m in ml]) + '\n'

                resstring += '\n'

            else:
                htext, rtext = ',', ''
                fnames = [r[0] for r in rl if r[0] != 'All']
                frames = sorted(set(fnames), key=fnames.index)
                onames = [r[2] for r in rl if r[0] != 'All']
                obs = sorted(set(onames), key=onames.index)
                resnames = [r[3] for r in rl if r[0] != 'All']
                res = sorted(set(resnames), key=resnames.index)

                for f in frames:
                    for o in obs:
                        htext += f'{f} {o},'

                htext += '\n'

                for mi, m in enumerate(res):
                    for r in rl:
                        if r[3] == m and r[0] != 'All':
                            if m not in rtext:
                                rtext += '{}, {},'.format(m, r[4])
                            else:
                                rtext += '{},'.format(r[4])
                    rtext += '\n'

                resstring = htext + rtext

        else:
            if (len(set(zrl[0])) > 1 and node.animated) or set(zrl[0]) == {'All'}:
                resstring = ''.join(['{} {},'.format(r[2], r[3]) for r in rl if r[0] == 'All']) + '\n'
                metriclist = list(zip(*[r.split() for ri, r in enumerate(zrl[4]) if zrl[0][ri] == 'All']))
            else:
                resstring = ''.join(['{} {} {},'.format(r[0], r[2], r[3]) for r in rl if r[0] != 'All']) + '\n'
                metriclist = list(itertools.zip_longest(*[r.split() for ri, r in enumerate(zrl[4]) if zrl[0][ri] != 'All'], fillvalue=''))

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
        dp = bpy.context.evaluated_depsgraph_get()
        scene = context.scene
        svp = scene.vi_params
        dobs = [o for o in bpy.data.objects if o.vi_params.vi_type == '2' and o.visible_get()]
        gobs = [o for o in bpy.data.objects if o.vi_params.vi_type == '2' and o.visible_get()]

        if viparams(self, scene):
            return {'CANCELLED'}

        if len(dobs) != 1:
            self.report({'ERROR'}, "One, and only one object with the CFD Domain property is allowed")
            return {'CANCELLED'}

        if [f.material_index for f in dobs[0].data.polygons if f.material_index + 1 > len(dobs[0].data.materials)]:
            self.report({'ERROR'}, "Not every domain face has a material attached")
            logentry("Not every face has a material attached")
            return {'CANCELLED'}

        if [dobs[0].material_slots[f.material_index].material for f in dobs[0].data.polygons if ' ' in dobs[0].material_slots[f.material_index].material.name]:
            self.report({'ERROR'}, "There is a space in one of the domain boundary material names")
            logentry("There is a space in one of the boundary material names")
            return {'CANCELLED'}

        for gob in gobs:
            if [f.material_index for f in gob.data.polygons if f.material_index + 1 > len(gob.data.materials)]:
                self.report({'ERROR'}, f"Not every face for object {gob.name} has a material attached")
                logentry(f"Not every face of object {gob.name} has a material attached")
                return {'CANCELLED'}

            if [gob.material_slots[f.material_index].material for f in gob.data.polygons if ' ' in gob.material_slots[f.material_index].material.name]:
                self.report({'ERROR'}, f"There is a space in one of the {gob.name} material names")
                logentry(f"There is a space in one of the {gob.name} material names")
                return {'CANCELLED'}

        casenode = context.node
        casenode.pre_case(context)

        if casenode.parametric:
            frames = range(casenode.frame_start, casenode.frame_end + 1)
        else:
            frames = [scene.frame_current]

        svp['flparams']['start_frame'] = frames[0]
        svp['flparams']['end_frame'] = frames[-1]

        for frame in frames:
            scene.frame_set(frame)
            frame_offb = os.path.join(svp['flparams']['offilebase'], str(frame))
            frame_ofcfb = os.path.join(frame_offb, 'constant')
            frame_ofsfb = os.path.join(frame_offb, 'system')

            for ofdir in (frame_offb, frame_ofcfb, frame_ofsfb):
                if not os.path.isdir(ofdir):
                    os.makedirs(ofdir)

            for f in os.listdir(frame_offb):
                try:
                    os.remove(os.path.join(frame_offb, f))
                except Exception:
                    pass

            for f in os.listdir(frame_ofcfb):
                try:
                    os.remove(os.path.join(frame_ofcfb, f))
                except Exception:
                    pass

            for root, dirs, files in os.walk(os.path.join(frame_ofcfb, 'postProcessing')):
                for d in dirs:
                    try:
                        shutil.rmtree(os.path.join(root, d))
                    except Exception:
                        pass

            svp['flparams']['et'] = casenode.etime
            svp['flparams']['features'] = {'turb': {'kEpsilon': 'kE'}}
            svp['flparams']['features']['rad'] = casenode.buoyancy and casenode.radiation
            svp['flparams']['radmodel'] = casenode.radmodel
            svp['flparams']['features']['buoy'] = casenode.buoyancy
            base_residuals = ['Ux', 'Uy', 'Uz']
            turb_residuals = ['k', 'epsilon']
            rad_residuals = ['G'] if svp['flparams']['features']['rad'] else []
            buoy_residuals = ['p_rgh', 'e'] if casenode.buoyancy else ['p']

            if casenode.buoyancy:
                svp['flparams']['pref'] = casenode.pabsval
                svp['flparams']['solver'] = 'buoyantFoam'
                svp['flparams']['solver_type'] = 'bf'
                svp['flparams']['params'] = 'bke'
            else:
                svp['flparams']['solver'] = 'simpleFoam'
                svp['flparams']['pref'] = casenode.pnormval
                svp['flparams']['solver_type'] = 'sf'
                svp['flparams']['params'] = 'ke'

            svp['flparams']['residuals'] = base_residuals + buoy_residuals + turb_residuals + rad_residuals
            svp['flparams']['st'] = 0
            svp['flparams']['presid'] = casenode.presid
            svp['flparams']['uresid'] = casenode.uresid
            svp['flparams']['keoresid'] = casenode.keoresid

            with open(os.path.join(frame_ofsfb, 'controlDict'), 'w') as cdfile:
                cdfile.write(fvcdwrite(svp, casenode, dp))
            with open(os.path.join(frame_ofsfb, 'fvSolution'), 'w') as fvsolfile:
                fvsolfile.write(fvsolwrite(svp, casenode))
            with open(os.path.join(frame_ofsfb, 'fvSchemes'), 'w') as fvschfile:
                fvschfile.write(fvschwrite(svp, casenode))
            with open(os.path.join(frame_ofcfb, 'momentumTransport'), 'w') as mtfile:
                mtfile.write(fvmtwrite(casenode, svp['flparams']['features']))

            if casenode.buoyancy:
                with open(os.path.join(frame_ofcfb, 'pRef'), 'w') as pfile:
                    pfile.write(fvprefwrite(casenode))
                with open(os.path.join(frame_ofcfb, 'physicalProperties'), 'w') as ppfile:
                    ppfile.write(fvtppwrite(svp))
                with open(os.path.join(frame_ofcfb, 'g'), 'w') as gfile:
                    gfile.write(fvgwrite())

                if casenode.radiation:
                    with open(os.path.join(frame_ofcfb, 'radiationProperties'), 'w') as rpfile:
                        rpfile.write(fvrpwrite(casenode))
                    with open(os.path.join(frame_ofcfb, 'fvModels'), 'w') as fvmfile:
                        fvmfile.write(fvmodwrite(casenode))

            else:
                with open(os.path.join(frame_ofcfb, 'physicalProperties'), 'w') as ppfile:
                    ppfile.write(fvtpwrite())

        casenode.post_case()
        return {'FINISHED'}


class NODE_OT_Flo_NG(bpy.types.Operator):
    bl_idname = "node.flovi_ng"
    bl_label = "NetGen export"
    bl_description = "Create a Netgen mesh"
    bl_register = True
    bl_undo = False

    def invoke(self, context, event):
        try:
            bpy.ops.object.mode_set(mode='OBJECT')
        except Exception:
            pass

        self.vi_prefs = bpy.context.preferences.addons[__name__.split('.')[0]].preferences
        addonpath = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
        scene = context.scene
        svp = scene.vi_params
        self.offb = svp['flparams']['offilebase']
        self.vl = context.view_layer
        self.expnode = context.node
        self.expnode.running = 1
        self.surf_complete = 0
        self.vol_complete = 0
        self.vol_running = 0
        case_nodes = [link.from_node for link in self.expnode.inputs['Case in'].links]
        bound_nodes = [link.to_node for link in self.expnode.outputs['Mesh out'].links]

        for node in case_nodes + bound_nodes:
            node.use_custom_color = 1

        self.curcoll = context.collection
        meshcoll = create_coll(context, 'FloVi Mesh')
        clear_coll(context, meshcoll)
        dp = bpy.context.evaluated_depsgraph_get()
        dobs = [o for o in bpy.data.objects if o.vi_params.vi_type == '2' and o.visible_get() and o.name not in meshcoll.objects]

        if not dobs:
            logentry('FloVi requires a domain object but none was found. Check the domain object is not hidden or in the FloVi Mesh collection')
            self.report({'ERROR'}, 'No, or hidden, domain objects')
            self.expnode.running = 0
            return {'CANCELLED'}

        elif len(dobs) > 1:
            self.report({'WARNING'}, 'More then one domain object found. Only the first is exported')

        gobs = [o for o in bpy.data.objects if o.vi_params.vi_type == '3' and o.visible_get() and o.name not in meshcoll.objects]
        self.obs = dobs + gobs

        if any([any(s < 0 for s in o.scale) for o in self.obs]):
            logentry('CFD domain or geometry has a negative scale. Cannot proceed')
            self.report({'ERROR'}, 'CFD domain or geometry has a negative scale. Cannot proceed')
            self.expnode.running = 0
            return {'CANCELLED'}

        for fn in ('ng.mesh', 'ng.vol', 'ng_surf.stl', 'ng_surf.vol'):
            if os.path.isfile(os.path.join(self.offb, fn)):
                os.remove(os.path.join(self.offb, fn))

        if os.environ.get('LD_LIBRARY_PATH'):
            os.environ['LD_LIBRARY_PATH'] += os.pathsep + os.path.join(addonpath, 'Python', sys.platform, 'netgen')
        else:
            os.environ['LD_LIBRARY_PATH'] = os.path.join(addonpath, 'Python', sys.platform, 'netgen')

        mp = MeshingParameters(maxh=self.expnode.maxcs, minh=0.25 * self.expnode.maxcs, grading=self.expnode.grading,
                               optsteps2d=self.expnode.optimisations, optsteps3d=self.expnode.optimisations,
                               delaunay=True, maxoutersteps=self.expnode.maxsteps)

        SetNumThreads(int(svp['viparams']['nproc']))
        mns = [0]
        self.omats = []
        mats, self.matnames, g_geos = [], [], []
        pmap1 = {}
        totmesh = Mesh()
        totmesh.SetMaterial(1, 'air')
        surf_no = 0
        fds = []
        b_mats = []
        e_maxs = {}

        for ob in self.obs:
            fi = 0
            faces = []
            bm = bmesh.new()
            bm.from_object(ob, dp)
            bm.transform(ob.matrix_world)
            min_elen = min([edge.calc_length() for edge in bm.edges])

            if not ob.material_slots:
                logentry(f'{ob.name} has faces with an unspecified material or an empty material slot')
                self.report({'ERROR'}, f'{ob.name} has faces with an unspecified material or an empty material slot')
                self.expnode.running = 0
                return {'CANCELLED'}

            if ' ' in ob.name:
                logentry(f'{ob.name} has a space in the name')
                self.report({'ERROR'}, f'{ob.name} has a space in the name')
                self.expnode.running = 0
                return {'CANCELLED'}

            if len(bm.faces) > 50000:
                bm.free()
                logentry('{} has more than 50000 faces. Simplify the geometry'.format(ob.name))
                self.report({'ERROR'}, '{} has more than 50000 faces. Simplify the geometry'.format(ob.name))
                self.expnode.running = 0
                return {'CANCELLED'}

            if not all([e.is_manifold for e in bm.edges]) or not all([v.is_manifold for v in bm.verts]):
                bmesh.ops.remove_doubles(bm, verts=bm.verts, dist=0.00001)

                if not all([e.is_manifold for e in bm.edges]) or not all([v.is_manifold for v in bm.verts]):
                    bm.free()
                    logentry('FloVi error: {} is not manifold'.format(ob.name))
                    self.report({'ERROR'}, 'FloVi error: {} is not manifold'.format(ob.name))
                    self.expnode.running = 0
                    return {'CANCELLED'}
                else:
                    logentry('{} had double vertices removed to be manifold. Check its mesh geometry'.format(ob.name))
                    self.report({'WARNING'}, '{} had double vertices removed to be manifold. Check its mesh geometry'.format(ob.name))

            for mi, ms in enumerate(ob.material_slots):
                for face in bm.faces:
                    if face.material_index == mi:
                        if ms.material and f'{ob.name}_{ms.material.name}' not in b_mats:
                            b_mats.append(f'{ob.name}_{ms.material.name}')
                        elif not ms.material:
                            logentry(f'{ob.name} has faces with an unspecified material or an empty material slot')
                            self.report({'ERROR'}, f'{ob.name} has faces with an unspecified material or an empty material slot')
                            self.expnode.running = 0
                            bm.free()
                            return {'CANCELLED'}
                        elif ' ' in ms.material.name:
                            logentry(f'The material {ms.material.name} has a space in the name')
                            self.report({'ERROR'}, f'The material {ms.material.name} has a space in the name')
                            self.expnode.running = 0
                            bm.free()
                            return {'CANCELLED'}

                        face.index = fi
                        fi += 1

            bm.faces.sort()
            bm.normal_update()
            bm.faces.ensure_lookup_table()
            ob_mesh = bpy.data.meshes.new("mesh")
            bm.to_mesh(ob_mesh)
            bmesh.ops.split_edges(bm, edges=bm.edges)

            for mi, ms in enumerate(ob.material_slots):
                if ms.material:
                    fd = FaceDescriptor(bc=surf_no, domin=1, surfnr=surf_no + 1)
                    self.matnames.append(ob.name + '_' + ms.material.name)
                    e_maxs[ob.name + '_' + ms.material.name] = ms.material.vi_params.flovi_ng_emax
                    self.omats.append(ms.material)
                    fd.bcname = ms.material.name
                    fd.color = ms.material.diffuse_color[:3]
                    fds.append(fd)
                    totmesh.Add(fd)
                    surf_no += 1

            lbm = len(bm.faces)

            for fi, face in enumerate(bm.faces[:lbm]):
                try:
                    matname = ob.material_slots[face.material_index].material.name

                    if ' ' in matname:
                        logentry(f'FloVi error: Material {matname} has a space in the name')
                        self.report({'ERROR'}, f'{matname} has a space in the name')
                        self.expnode.running = 0
                        bm.free()
                        return {'CANCELLED'}

                except Exception:
                    logentry(f'FloVi error: {ob.name} mesh has faces that reference a non-existant material')
                    self.report({'ERROR'}, f'{ob.name} mesh has faces that reference a non-existant material')
                    self.expnode.running = 0
                    bm.free()
                    return {'CANCELLED'}

                edges = [occ.Segment(occ.gp_Pnt(tuple(loop.vert.co)), occ.gp_Pnt(tuple(loop.link_loop_next.vert.co))) for loop in face.loops]
                wire = occ.Wire(edges)
                f = occ.Face(wire)

                if len(f.edges) > 2:
                    f.name = ob.name + '_' + matname
                    f.mat(matname)
                    f.bc(ob.name + '_' + matname)
                    f.layer = face.index
                    f.maxh = ob.material_slots[face.material_index].material.vi_params.flovi_ng_max
                    faces.append(f)

                else:
                    fc = Vector([fc for fc in face.calc_center_bounds()])
                    fn = Vector([fn for fn in face.normal])
                    evs = [(loop.vert, loop.link_loop_next.vert) for loop in face.loops]
                    vs = [loop.vert for loop in face.loops]

                    for v in vs:
                        dist = distance_point_to_plane(v.co, fc, fn)

                        if abs(dist) < min_elen * 0.4:
                            v.co -= dist * fn

                    edges = [occ.Segment(occ.gp_Pnt(tuple(loop.vert.co)), occ.gp_Pnt(tuple(loop.link_loop_next.vert.co))) for loop in face.loops]
                    wire = occ.Wire(edges)
                    f = occ.Face(wire)

                    if len(f.edges) > 2:
                        f.name = ob.name + '_' + matname
                        f.mat(matname)
                        f.bc(ob.name + '_' + matname)
                        f.layer = face.index
                        f.maxh = ob.material_slots[face.material_index].material.vi_params.flovi_ng_max
                        faces.append(f)
                    else:
                        logentry(f'Object {ob.name} face with index {face.index} had to be triangulated. This could lead to poor mesh quality')
                        t_faces = bmesh.ops.triangulate(bm, faces=[face], quad_method='BEAUTY', ngon_method='BEAUTY')['faces']

                        for ti, tf in enumerate(t_faces):
                            edges = [occ.Edge(occ.Vertex(occ.gp_Pnt(tuple(loop.vert.co))), occ.Vertex(occ.gp_Pnt(tuple(loop.link_loop_next.vert.co)))) for loop in tf.loops]
                            wire = occ.Wire(edges)
                            f = occ.Face(wire)
                            f.name = ob.name + '_' + matname
                            f.mat(matname)
                            f.bc(ob.name + '_' + matname)
                            f.layer = face.index
                            f.maxh = ob.material_slots[face.material_index].material.vi_params.flovi_ng_max
                            faces.append(f)

            bm.free()

            if ob == dobs[0]:
                d_geo = occ.OCCGeometry(occ.Compound(faces))
                fns = [face.name for face in d_geo.shape.faces]
                fms = [face.maxh for face in d_geo.shape.faces]
                d_geo.Heal(tolerance=min_elen * 0.8)

                if None in set([face.name for face in d_geo.shape.faces]):
                    for fi, face in enumerate(d_geo.shape.faces):
                        if face.name is None:
                            face.name = fns[fi]
                            face.maxh = fms[fi]

                if self.expnode.debug_step:
                    d_geo.shape.WriteStep(os.path.join(svp['flparams']['offilebase'], 'empty_domain.step'))

                if len(d_geo.shape.SubShapes(occ.SOLID)) != 1:
                    logentry(f'FloVi error: {ob.name} cannot be converted to a single solid')
                    self.report({'ERROR'}, f'{ob.name} cannot be converted to a single solid')
                    self.expnode.running = 0
                    return {'CANCELLED'}

            else:
                mesh_islands = bpy_extras.mesh_utils.mesh_linked_triangles(ob_mesh)

                if len(mesh_islands) > 1:
                    for mi, mesh_island in enumerate(mesh_islands):
                        g_geo = occ.OCCGeometry(occ.Compound([face for face in faces if face.layer in set([f.polygon_index for f in mesh_island])]))
                        fns = [face.name for face in g_geo.shape.faces]
                        fms = [face.maxh for face in g_geo.shape.faces]
                        fcs = [face.center for face in g_geo.shape.faces]
                        g_geo.Heal(tolerance=min_elen * 0.8)

                        if len(g_geo.shape.SubShapes(occ.SOLID)):
                            for g_geo_solid in g_geo.shape.SubShapes(occ.SOLID):
                                if not all([face.name for face in g_geo_solid.faces]):
                                    for fi, face in enumerate(g_geo_solid.faces):
                                        if face.name is None:
                                            for fci, fc in enumerate(fcs):
                                                if (Vector(face.center) - Vector(fc)).length < 0.001:
                                                    face.name = fns[fci]
                                                    face.maxh = fms[fci]
                                                    break

                                            if face.name is None:
                                                face.name = fns[fi]
                                                face.maxh = fms[fi]

                                g_geos.append(g_geo_solid)

                        else:
                            g_geo.shape.WriteStep(os.path.join(svp['flparams']['offilebase'], f'{ob.name}.step'))
                            g_geo = occ.OCCGeometry(os.path.join(svp['flparams']['offilebase'], f'{ob.name}.step'))
                            g_geo.Heal(tolerance=min_elen * 0.8)

                            for g_geo_solid in g_geo.shape.SubShapes(occ.SOLID):
                                if not all([face.name for face in g_geo_solid.faces]):
                                    for fi, face in enumerate(g_geo_solid.faces):
                                        if face.name is None or face.name == '':
                                            for fci, fc in enumerate(fcs):
                                                if (Vector(face.center) - Vector(fc)).length < 0.001:
                                                    face.name = fns[fci]
                                                    face.maxh = fms[fci]
                                                    break
                                            if face.name is None or face.name == '':
                                                face.name = fns[fi]
                                                face.maxh = fms[fi]

                                g_geos.append(g_geo_solid)

                        if not len(g_geo.shape.SubShapes(occ.SOLID)):
                            logentry(f'FloVi warning: {ob.name} shell {mi} cannot be converted to a solid')
                            self.report({'WARNING'}, f'{ob.name} shell {mi} cannot be converted to a solid')

                        if self.expnode.debug_step:
                            g_geo.shape.WriteStep(os.path.join(svp['flparams']['offilebase'], f'{ob.name}_{mi}.step'))

                else:
                    g_geo = occ.OCCGeometry(occ.Compound(faces))
                    fns = [face.name for face in g_geo.shape.faces]
                    fms = [face.maxh for face in g_geo.shape.faces]
                    fcs = [face.center for face in g_geo.shape.faces]
                    g_geo.Heal(tolerance=min_elen * 0.8)

                    for g_geo_solid in g_geo.shape.SubShapes(occ.SOLID):
                        if not all([face.name for face in g_geo.shape.faces]):
                            for fi, face in enumerate(g_geo_solid.faces):
                                if face.name is None or face.name == '':
                                    for fci, fc in enumerate(fcs):
                                        if (Vector(face.center) - Vector(fc)).length < 0.001:
                                            face.name = fns[fci]
                                            face.maxh = fms[fci]
                                            break
                                    if face.name is None or face.name == '':
                                        face.name = fns[fi]
                                        face.maxh = fms[fi]

                        g_geos.append(g_geo_solid)

                    if not len(g_geo.shape.SubShapes(occ.SOLID)):
                        g_geo.shape.WriteStep(os.path.join(svp['flparams']['offilebase'], f'{ob.name}.step'))
                        g_geo = occ.OCCGeometry(os.path.join(svp['flparams']['offilebase'], f'{ob.name}.step'))
                        fns = [face.name for face in g_geo.shape.faces]
                        fms = [face.maxh for face in g_geo.shape.faces]
                        fcs = [face.center for face in g_geo.shape.faces]
                        g_geo.Heal(tolerance=min_elen * 0.8)

                        for g_geo_solid in g_geo.shape.SubShapes(occ.SOLID):
                            if None in set([face.name for face in g_geo_solid.faces]):
                                for fi, face in enumerate(g_geo_solid.faces):
                                    if face.name is None:
                                        for fci, fc in enumerate(fcs):
                                            if (Vector(face.center) - Vector(fc)).length < 0.001:
                                                face.name = fns[fci]
                                                face.maxh = fms[fci]
                                                break

                                        if face.name is None:
                                            face.name = fns[fi]
                                            face.maxh = fms[fi]

                                g_geos.append(g_geo_solid)

                            if not len(g_geo.shape.SubShapes(occ.SOLID)):
                                logentry(f'FloVi warning: {ob.name} cannot be converted to a solid')
                                self.report({'WARNING'}, f'{ob.name} cannot be converted to a solid')

                    if self.expnode.debug_step:
                        g_geo.shape.WriteStep(os.path.join(svp['flparams']['offilebase'], f'{ob.name}.step'))

        totmesh.Save(os.path.join(svp['flparams']['offilebase'], 'ng.vol'))

        if g_geos:
            for gi, g_solid in enumerate([g for g in g_geos]):
                d_geo = d_geo.shape - g_solid if not gi else d_geo - g_solid

        else:
            d_geo = d_geo.shape

        fns = [face.name for face in d_geo.faces]
        fms = [face.maxh for face in d_geo.faces]
        fcs = [face.center for face in d_geo.faces]
        d_geo = occ.OCCGeometry(d_geo)
        d_geo.Heal(tolerance=0.001)

        if None in set([face.name for face in d_geo.shape.faces]):
            for fi, face in enumerate(d_geo.shape.faces):
                if face.name is None:
                    for fci, fc in enumerate(fcs):
                        if (Vector(face.center) - Vector(fc)).length < 0.001:
                            face.name = fns[fci]
                            face.maxh = fms[fci]
                            break
                    if face.name is None:
                        face.name = fns[fi]
                        face.maxh = fms[fi]

        d_geo.shape.WriteStep(os.path.join(svp['flparams']['offilebase'], 'flovi_geometry.step'))
        self.mis = [self.matnames.index(face.name) for face in d_geo.shape.faces]

        with open(os.path.join(svp['flparams']['offilebase'], 'ngpy.py'), 'w') as ngpyfile:
            ngpyfile.write(inspect.cleandoc('''
            import os, math
            import numpy as np
            from netgen import occ
            from netgen.meshing import MeshingParameters, FaceDescriptor, Element2D, Mesh, MeshingStep, meshsize
            from pyngcore import SetNumThreads, TaskManager
            SetNumThreads({0})
            mp = MeshingParameters(meshsize.fine, maxh={1}, minh={2}, grading={3}, optimize3d='cmdmustm', optsteps3d={4}, delaunay=True, maxoutersteps={5})
            geo = occ.OCCGeometry(os.path.join(r'{6}', 'flovi_geometry.step'))
            e_maxs = {7}

            def unit_vector(vector):
                """ Returns the unit vector of the vector.  """
                return vector / np.linalg.norm(vector)

            def angle_between(v1, v2):
                v1_u = unit_vector(v1)
                v2_u = unit_vector(v2)
                return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

            for edge in geo.shape.edges:
                vecs = []
                e_maxhs = []

                for face in geo.shape.faces:
                    if edge in face.edges:
                        f_dict = eval('{{'+str(face)+'}}')

                        if 'TShape' in f_dict and face.name:
                            f_dir = f_dict['TShape']['Surface']['pos']['Direction']
                            vecs.append(f_dir)
                            e_maxhs.append(e_maxs[face.name] * face.maxh)

                if len(vecs) == 2 and abs(angle_between(vecs[0], vecs[1]) > 0.1 or len(set(e_maxhs)) > 1):
                    edge.maxh = min(e_maxhs)
                    e_len = ((edge.vertices[0].p[0] - edge.vertices[1].p[0])**2 + (edge.vertices[0].p[1] - edge.vertices[1].p[1])**2 + (edge.vertices[0].p[2] - edge.vertices[1].p[2])**2)**0.5

                    if e_len > 2 * edge.maxh:
                        segs = int(e_len/edge.maxh) + 1
                        for s in range(1, segs):
                            vco = [edge.vertices[0].p[i] + (edge.vertices[1].p[i]* s/segs - edge.vertices[0].p[i]* s/segs) for i in range(3)]
                            mp.RestrictH(x=vco[0], y=vco[1], z=vco[2], h=edge.maxh)

                        for v in edge.vertices:
                            mp.RestrictH(x=v.p[0], y=v.p[1], z=v.p[2], h=edge.maxh)

            with TaskManager():
                surf_mesh = geo.GenerateMesh(mp=mp, perfstepsend=MeshingStep.MESHSURFACE)

            #surf_mesh.Refine(adaptive=True)
            surf_mesh.Save(os.path.join(r'{6}', 'ng_surf.vol'))
            surf_mesh.Export(os.path.join(r'{6}', 'ng_surf.stl'), 'STL Format')
            '''.format(int(svp['viparams']['nproc']), self.expnode.maxcs, 0.0, self.expnode.grading,
                       self.expnode.optimisations, self.expnode.maxsteps, svp['flparams']['offilebase'], e_maxs)))

        self.surf_run = Popen(shlex.split('"{}" "{}"'.format(sys.executable, os.path.join(svp['flparams']['offilebase'], 'ngpy.py'))), stdout=PIPE, stderr=PIPE)
        self.surf_cancel = cancel_window(os.path.join(svp['viparams']['newdir'], 'viprogress'), pdll_path, 'Surface Mesh')
        self._timer = context.window_manager.event_timer_add(2, window=context.window)
        context.window_manager.modal_handler_add(self)
        return {'RUNNING_MODAL'}

    def modal(self, context, event):
        scene = context.scene
        svp = scene.vi_params

        if self.surf_run.poll() is None:
            if self.surf_cancel.poll() is not None:
                self.surf_run.kill()
                self.expnode.running = 0
                return {'CANCELLED'}
        else:
            surf_err_lines = []

            for surf_err, line in enumerate(self.surf_run.stderr):
                surf_err_lines.append(line.decode())
                logentry(f"Surface mesh error: {line.decode()}")

            if surf_err_lines:
                self.report({'ERROR'}, "Surface meshing failed. Check the vi-suite-log file in Blender's text editor")

                if self.surf_cancel.poll() is None:
                    self.surf_cancel.kill()

                self.expnode.running = 0
                return {'CANCELLED'}
            else:
                self.surf_cancel.kill()
                self.surf_complete = 1

        if self.surf_complete and not self.vol_complete and not self.vol_running:
            with open(os.path.join(svp['flparams']['offilebase'], 'ngpy.py'), 'w') as ngpyfile:
                ngpyfile.write(inspect.cleandoc('''
                import os
                from netgen.meshing import MeshingParameters, FaceDescriptor, Element2D, Mesh
                from pyngcore import SetNumThreads, TaskManager
                mp = MeshingParameters(maxh={3}, grading={4}, optsteps3d={5}, optimize3d='cmdmustm')
                SetNumThreads({0})
                surf_mesh = Mesh()
                surf_mesh.Load(os.path.join(r'{1}', 'ng_surf.vol'))
                tot_mesh = Mesh()
                tot_mesh.Load(os.path.join(r'{1}', 'ng.vol'))
                mis = {2}
                pmap = {{}}

                for ei, el in enumerate(surf_mesh.Elements2D()):
                    el.index = mis[el.index - 1] + 1

                    for v in el.vertices:
                        if (v not in pmap):
                            pmap[v] = tot_mesh.Add(surf_mesh[v])

                    tot_mesh.Add(Element2D(el.index, [pmap[v] for v in el.vertices]))

                with TaskManager():
                    tot_mesh.GenerateVolumeMesh(mp=mp)

                tot_mesh.Save(os.path.join(r'{1}', 'ng.vol'))
                tot_mesh.Export(os.path.join(r'{1}', 'ng.mesh'), format='Neutral Format')
                '''.format(int(svp['viparams']['nproc']), svp['flparams']['offilebase'], self.mis, self.expnode.maxcs, self.expnode.grading,
                           self.expnode.optimisations)))

            self.vol_running = 1
            self.vol_run = Popen(shlex.split('"{}" "{}"'.format(sys.executable, os.path.join(svp['flparams']['offilebase'], 'ngpy.py'))), stdout=PIPE, stderr=PIPE)
            self.vol_cancel = cancel_window(os.path.join(svp['viparams']['newdir'], 'viprogress'), pdll_path, 'Volume Mesh')

        if self.vol_running == 1:
            if self.vol_run.poll() is None:
                if self.vol_cancel.poll() is not None:
                    self.vol_run.kill()
                    self.expnode.running = 0
                    return {'CANCELLED'}

            else:
                vol_err_lines = []

                for line in self.vol_run.stderr:
                    vol_err_lines.append(line.decode())
                    logentry(f"FloVi volume mesh error: {line.decode()}")

                if vol_err_lines:
                    if self.vol_cancel.poll() is None:
                        self.vol_cancel.kill()

                    self.expnode.running = 0
                    self.report({'ERROR'}, "Volume meshing failed. Check the vi-suite-log file in Blender's text editor")
                    return {'CANCELLED'}
                else:
                    self.vol_cancel.kill()
                    self.vol_complete = 1

        if self.surf_complete and self.vol_complete:
            self.conv_cancel = cancel_window(os.path.join(svp['viparams']['newdir'], 'viprogress'), pdll_path, 'Converting to OpenFOAM')

            for frame in range(svp['flparams']['start_frame'], svp['flparams']['end_frame'] + 1):
                frame_offb = os.path.join(svp['flparams']['offilebase'], str(frame))
                frame_ofcfb = os.path.join(frame_offb, 'constant')
                st = '0'

                if os.path.isfile(os.path.join(frame_offb, st, 'polyMesh', 'points')):
                    os.remove(os.path.join(frame_offb, st, 'polyMesh', 'points'))

                pdm_error = 0
                scene = context.scene
                svp = scene.vi_params

                if os.path.isfile(os.path.join(self.offb, 'ng.mesh')):
                    os.chdir(self.offb)

                    if sys.platform == 'linux' and os.path.isdir(self.vi_prefs.ofbin):
                        nntf_cmd = 'foamExec netgenNeutralToFoam -noFunctionObjects -case {} {}'.format(frame_offb, os.path.join(self.offb, 'ng.mesh'))
                        subprocess.Popen(shlex.split(nntf_cmd)).wait()

                    elif sys.platform in ('darwin', 'win32'):
                        nntf_cmd = '{} run -it --rm -v "{}":/home/openfoam/data dicehub/openfoam:12 "netgenNeutralToFoam -case data/{} {}"'.format(docker_path, self.offb, frame, 'data/ng.mesh')
                        subprocess.Popen(nntf_cmd, shell=True).wait()

                    logentry(f'Running netgenNeutraltoFoam with command: {nntf_cmd}')

                    if not os.path.isdir(os.path.join(frame_offb, st, 'polyMesh')):
                        os.makedirs(os.path.join(frame_offb, st, 'polyMesh'))

                elif not os.path.isfile(os.path.join(svp['flparams']['offilebase'], 'ng.mesh')):
                    logentry('Netgen volume meshing did not complete')
                    self.expnode.running = 0
                    self.conv_cancel.kill()
                    self.report({'ERROR'}, 'Netgen volume meshing did not complete')
                    return {'CANCELLED'}

                if self.conv_cancel.poll() is not None:
                    self.expnode.running = 0
                    return {'CANCELLED'}

                if os.path.isfile(os.path.join(frame_ofcfb, 'polyMesh', 'boundary')):
                    with open(os.path.join(frame_ofcfb, 'polyMesh', 'boundary'), 'r') as bfile:
                        nf = []
                        ns = []

                        for line in bfile.readlines():
                            if 'nFaces' in line:
                                nf.append(int(line.split()[1].strip(';')))
                            if 'startFace' in line:
                                ns.append(int(line.split()[1].strip(';')))

                    with open(os.path.join(frame_ofcfb, 'polyMesh', 'boundary'), 'w') as bfile:
                        bfile.write(ofheader)
                        cl = 'polyBoundaryMesh' if self.expnode.bl_label == 'FloVi NetGen' else 'BoundaryMesh'
                        loc = 'constant/polyMesh' if self.expnode.bl_label == 'FloVi NetGen' else 'Mesh'
                        bfile.write(write_ffile(cl, loc, 'boundary'))
                        bfile.write('// **\n\n{}\n(\n'.format(len(ns)))
                        omi = 0

                        for mi, mat in enumerate(self.omats):
                            if omi < len(ns):
                                bfile.write(write_bound(self.matnames[mi], mat, ns[omi], nf[omi]))
                                omi += 1

                        bfile.write(')\n\n// **\n')

                    for file in os.listdir(os.path.join(frame_ofcfb, 'polyMesh')):
                        shutil.copy(os.path.join(os.path.join(frame_ofcfb, 'polyMesh'), file),
                                    os.path.join(frame_offb, st, 'polyMesh'))

                    if self.expnode.poly:
                        if sys.platform == 'linux' and os.path.isdir(self.vi_prefs.ofbin):
                            pdm = Popen(shlex.split('foamExec polyDualMesh -case ./{} -concaveMultiCells -noFunctionObjects -overwrite {}'.format(frame, 5)), stdout=PIPE, stderr=PIPE)

                        elif sys.platform in ('darwin', 'win32'):
                            pdm_cmd = '{} run -it --rm -v "{}":/home/openfoam/data dicehub/openfoam:12 "polyDualMesh -case data -concaveMultiCells -noFunctionObjects -overwrite {}"'.format(docker_path, frame_offb, 5)
                            pdm = Popen(pdm_cmd, shell=True, stdout=PIPE, stderr=PIPE)

                        for line in pdm.stdout:
                            if 'FOAM aborting' in line.decode():
                                logentry('polyDualMesh error. Check the mesh in Netgen')
                                pdm_error = 1

                    if self.conv_cancel.poll() is not None:
                        self.expnode.running = 0
                        return {'CANCELLED'}

                    if not pdm_error:
                        if sys.platform == 'linux':
                            cpf_cmd = 'foamExec combinePatchFaces -overwrite -noFunctionObjects -case {} {}'.format(frame_offb, 5)
                            Popen(shlex.split(cpf_cmd)).wait()
                        elif sys.platform in ('darwin', 'win32'):
                            cpf_cmd = '{} run -it --rm -v "{}":/home/openfoam/data dicehub/openfoam:12 "combinePatchFaces -overwrite -noFunctionObjects -case data {}"'.format(docker_path, frame_offb, 5)
                            Popen(cpf_cmd, shell=True).wait()

                        # if sys.platform == 'linux':
                        #     cm = Popen(shlex.split('foamExec checkMesh -case {}'.format(frame_offb)), stdout=PIPE)
                        # elif sys.platform in ('darwin', 'win32'):
                        #     cm_cmd = 'docker run -it --rm -v "{}":/home/openfoam/data dicehub/openfoam:12 "checkMesh -case data"'.format(frame_offb)
                        #     cm = Popen(cm_cmd, shell=True, stdout=PIPE)

                        # for line in cm.stdout:
                        #     if '***Error' in line.decode():
                        #         logentry('Mesh errors:{}'.format(line.decode()))
                        #     elif '*Number' in line.decode() and sys.platform == 'linux':
                        #         Popen(shlex.split('foamExec foamToVTK -faceSet nonOrthoFaces -case {}'.format(frame_offb)), stdout=PIPE)
                        #     else:
                        #         print(line.decode())

                        for entry in os.scandir(os.path.join(frame_offb, st, 'polyMesh')):
                            if entry.is_file():
                                shutil.copy(os.path.join(frame_offb, st, 'polyMesh', entry.name), os.path.join(frame_ofcfb, 'polyMesh'))

                        with open(os.path.join(frame_ofcfb, 'polyMesh', 'boundary'), 'r') as bfile:
                            nf = []
                            ns = []

                            for line in bfile.readlines():
                                if 'nFaces' in line:
                                    nf.append(int(line.split()[1].strip(';')))
                                if 'startFace' in line:
                                    ns.append(int(line.split()[1].strip(';')))

                        with open(os.path.join(frame_ofcfb, 'polyMesh', 'boundary'), 'w') as bfile:
                            bfile.write(ofheader)
                            cl = 'polyBoundaryMesh' if self.expnode.bl_label == 'FloVi NetGen' else 'BoundaryMesh'
                            loc = 'constant/polyMesh' if self.expnode.bl_label == 'FloVi NetGen' else 'Mesh'
                            bfile.write(write_ffile(cl, loc, 'boundary'))
                            bfile.write('// **\n\n{}\n(\n'.format(len(ns)))
                            omi = 0

                            for mi, mat in enumerate(self.omats):
                                if omi < len(ns):
                                    bfile.write(write_bound(self.matnames[mi], mat, ns[omi], nf[omi]))
                                    omi += 1

                            bfile.write(')\n\n// **\n')

                        for file in os.listdir(os.path.join(frame_ofcfb, 'polyMesh')):
                            shutil.copy(os.path.join(os.path.join(frame_ofcfb, 'polyMesh'), file), os.path.join(frame_offb, st, 'polyMesh'))

                open("{}".format(os.path.join(frame_offb, '{}.foam'.format(frame))), "w")

            if self.conv_cancel.poll() is not None:
                self.expnode.running = 0
                return {'CANCELLED'}

            if os.path.isfile(os.path.join(frame_offb, st, 'polyMesh', 'points')):
                oftomesh(frame_offb, self.vl, self.omats, st, ns, nf, self.expnode.b_only)
            else:
                logentry('Netgen volume meshing failed:')
                self.report({'ERROR'}, 'Volume meshing failed')
                self.conv_cancel.kill()
                self.expnode.running = 0
                return {'CANCELLED'}

            self.expnode.post_export()
            self.conv_cancel.kill()
            create_coll(context, self.curcoll.name)
            self.expnode.running = 0
            return {'FINISHED'}

        else:
            return {'PASS_THROUGH'}


class NODE_OT_Flo_Bound(bpy.types.Operator):
    bl_idname = "node.flovi_bound"
    bl_label = "Boundary export"
    bl_description = "Export Openfoam boundaries"
    bl_register = True
    bl_undo = False

    def execute(self, context):
        scene = context.scene
        svp = scene.vi_params
        dobs = [o for o in bpy.data.objects if o.visible_get() and o.vi_params.vi_type == '2']
        gobs = [o for o in bpy.data.objects if o.visible_get() and o.vi_params.vi_type == '3']
        obs = dobs + gobs
        boundnode = context.node
        meshnode = boundnode.inputs['Mesh in'].links[0].from_node
        casenode = meshnode.inputs['Case in'].links[0].from_node
        offb = svp['flparams']['offilebase']
        b_dict = fvvarwrite(scene, obs, casenode)

        with open(os.path.join(scene.vi_params['viparams']['newdir'], 'boundary_summary.txt'), 'w') as b_file:
            for mat in b_dict:
                b_file.write(f'{mat}\n')

                for b in b_dict[mat]:
                    b_file.write(f'{b}\n')
                    b_file.write(f'{b_dict[mat][b]}\n')

        for frame in range(svp['flparams']['start_frame'], svp['flparams']['end_frame'] + 1):
            frame_offb = os.path.join(offb, str(frame))
            frame_ofcfb = os.path.join(frame_offb, 'constant')

            if os.path.isfile(os.path.join(frame_ofcfb, 'polyMesh', 'boundary')):
                with open(os.path.join(frame_ofcfb, 'polyMesh', 'boundary'), 'r') as bfile:
                    lines = []

                    for line in bfile.readlines():
                        for ob in obs:
                            for mat in ob.data.materials:
                                if line.strip() in b_dict and line.strip() == f'{ob.name}_{mat.name}':
                                    bound = flovi_b_dict[mat.vi_params.flovi_bmb_type]
                                if line.split() and line.split()[0] == 'type':
                                    line = f"        type            {bound};\n"

                        lines.append(line)
            else:
                logentry('No boundary file found. Re-create the mesh after exporting the case node')
                self.report({'ERROR'}, "No boundary file found. Re-create the mesh after exporting the case node")
                return {'CANCELLED'}

            with open(os.path.join(frame_ofcfb, 'polyMesh', 'boundary'), 'w') as bfile:
                bfile.write(''.join(lines))

        boundnode.post_export()
        return {'FINISHED'}


class NODE_OT_Flo_Sim(bpy.types.Operator):
    bl_idname = "node.flovi_sim"
    bl_label = "FloVi simulation"
    bl_description = "Solve an OpenFOAM case"
    bl_register = True
    bl_undo = True

    def invoke(self, context, event):
        wm = context.window_manager
        scene = context.scene
        svp = scene.vi_params
        self.simnode = context.node
        self.simnode.presim()
        self.convergence = svp['flparams']['uresid']
        self.econvergence = svp['flparams']['keoresid']
        self.pconvergence = svp['flparams']['presid']
        self.residuals = svp['flparams']['residuals']
        self.processes = self.simnode.processes
        self.fpfile = os.path.join(svp['viparams']['newdir'], 'floviprogress')
        self.pfile = fvprogressfile(svp['viparams']['newdir'])
        self.pb = qtfvprogress(os.path.join(svp['viparams']['newdir'], 'viprogress'), pdll_path, svp['flparams']['et'], str(self.residuals), svp['flparams']['start_frame'])
        self.pv = self.simnode.pv
        self.runs = []
        self.reslists = []
        self.o_dict = {}
        self.frames = range(svp['flparams']['start_frame'], svp['flparams']['end_frame'] + 1)
        self.simnode['frames'] = [f for f in self.frames]
        fframe_offb = os.path.join(svp['flparams']['offilebase'], str(svp['flparams']['start_frame']))
        os.chdir(svp['flparams']['offilebase'])

        for frame in self.frames:
            frame_offb = os.path.join(svp['flparams']['offilebase'], str(frame))

            for root, dirs, files in os.walk(frame_offb):
                for d in dirs:
                    if 'processor' in d:
                        shutil.rmtree(os.path.join(root, d))
                    try:
                        if float(d) != svp['flparams']['st']:
                            shutil.rmtree(os.path.join(root, d))
                    except Exception:
                        pass
                    if 'postProcessing' in d:
                        shutil.rmtree(os.path.join(root, d))

            if sys.platform == 'linux':
                pp_cmd = "foamExec foamPostProcess -func writeCellCentres -case {}".format(frame_offb)
                Popen(shlex.split(pp_cmd)).wait()
            elif sys.platform in ('darwin', 'win32'):
                pp_cmd = '{} run -it --rm -v "{}":/home/openfoam/data dicehub/openfoam:12 "foamPostProcess -func writeCellCentres -case data"'.format(docker_path, frame_offb)
                Popen(pp_cmd, shell=True).wait()

            if self.processes > 1:
                with open(os.path.join(frame_offb, 'system', 'decomposeParDict'), 'w') as fvdcpfile:
                    fvdcpfile.write(fvdcpwrite(self.processes))

                if sys.platform == 'linux':
                    dcp_cmd = "foamExec decomposePar -force -case {}".format(frame_offb)
                    Popen(shlex.split(dcp_cmd)).wait()

                elif sys.platform in ('darwin', 'win32'):
                    dcp_cmd = '{} run -it --rm -v "{}":/home/openfoam/data dicehub/openfoam:12 "decomposePar -force -case data"'.format(docker_path, frame_offb)
                    Popen(dcp_cmd, shell=True).wait()

        with open(self.fpfile, 'w') as fvprogress:
            if self.processes > 1:
                if sys.platform == 'linux':
                    self.runs.append(Popen(shlex.split('mpirun --oversubscribe -np {} foamExec {} -parallel -case {}'.format(self.processes,
                                                                                                                             svp['flparams']['solver'],
                                                                                                                             fframe_offb)), stdout=fvprogress))
                elif sys.platform in ('darwin', 'win32'):
                    self.runs.append(Popen('{} run -it --rm -v "{}":/home/openfoam/data dicehub/openfoam:12 "mpirun --oversubscribe -np {} {} -parallel -case data"'.format(docker_path, fframe_offb, self.processes,
                                                                                                                                                                                svp['flparams']['solver']), shell=True, stdout=fvprogress))
            else:
                if sys.platform == 'linux':
                    sol_cmd = '{} {} {} {}'.format('foamExec', svp['flparams']['solver'], "-case", fframe_offb)
                    self.runs.append(Popen(shlex.split(sol_cmd), stderr=PIPE, stdout=fvprogress))

                elif sys.platform in ('darwin', 'win32'):
                    sol_cmd = '{} run -it --rm -v "{}":/home/openfoam/data dicehub/openfoam:12 "{} -case data"'.format(docker_path, fframe_offb, svp['flparams']['solver'])
                    self.runs.append(Popen(sol_cmd, shell=True, stderr=PIPE, stdout=fvprogress))

                logentry('Running solver with command: {}'.format(sol_cmd))

        self._timer = wm.event_timer_add(5, window=context.window)
        wm.modal_handler_add(self)
        return {'RUNNING_MODAL'}

    def terminate(self, scene):
        for run in self.runs:
            run.kill()

        return {'CANCELLED'}

    def modal(self, context, event):
        svp = context.scene.vi_params
        frame_n = svp['flparams']['start_frame'] + len(self.runs)
        frame_c = svp['flparams']['start_frame'] + len(self.runs) - 1
        frame_noffb = os.path.join(svp['flparams']['offilebase'], str(frame_n))
        frame_coffb = os.path.join(svp['flparams']['offilebase'], str(frame_c))

        if self.runs[-1].poll() is None and self.pb.poll() is None:
            with open(self.fpfile, 'r') as fpfile:
                lines = fpfile.readlines()[::-1]
                residict = {}

                for line in lines:
                    line_split = line.split()

                    if 'Initial residual' in line:
                        residict[line_split[3][:-1]] = abs(float(line_split[line_split.index('Initial') + 3][:-1]))
                    elif len(line.split()) and line.split()[0] == 'Time' and line.split()[1] == '=':
                        residict[line.split()[0]] = float(line.split()[2].strip('s'))
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

        elif self.runs[-1].poll() is None and self.pb.poll() is not None:
            self.runs[-1].kill()
            
            if self.processes > 1:
                if sys.platform == 'linux':
                    Popen(shlex.split("foamExec reconstructPar -case {}".format(frame_coffb))).wait()
                elif sys.platform in ('darwin', 'win32'):
                    Popen('{} run -it --rm -v "{}":/home/openfoam/data dicehub/openfoam:12 "foamExec reconstructPar -case data"'.format(docker_path, frame_coffb), shell=True).wait()

            open("{}".format(os.path.join(frame_coffb, '{}.foam'.format(frame_c))), "w")
            self.simnode.running = 0
            logentry('Cancelling FloVi simulation')
            return {'CANCELLED'}

        elif self.pb.poll() is None or self.runs[-1].poll is not None:
            self.pb.kill()
            dline = ['', '']

            if self.runs[-1].stderr:
                for li, line in enumerate(self.runs[-1].stderr):
                    dline[0] = dline[1]
                    dline[1] = line.decode()

                    if 'Unable to set reference cell for field p' in dline[1]:
                        self.runs[-1].kill()
                        self.report({'ERROR'}, "Pressure reference point needs to be supplied or is outside the domain")
                        logentry('ERROR: Pressure reference point needs to be supplied or is outside the domain')
                        self.simnode.running = 0
                        return {'CANCELLED'}
                    elif 'Please supply either pRefCell or pRefPoint' in dline[1]:
                        self.runs[-1].kill()
                        self.report({'ERROR'}, "Pressure reference point needs to be supplied")
                        logentry('ERROR: Pressure reference point needs to be supplied')
                        self.simnode.running = 0
                        return {'CANCELLED'}
                    elif 'You are probably trying to solve for a field with a default boundary condition' in dline[1]:
                        dlist = dline[0].split()
                        pi = dlist.index('patch')
                        fi = dlist.index('field')
                        self.runs[-1].kill()
                        self.report({'ERROR'}, "Change the calculated boundary on {} for field {}".format(dlist[pi + 1], dlist[fi + 1]))
                        logentry('ERROR: Change the calculated boundary on {} for field {}'.format(dlist[pi + 1], dlist[fi + 1]))
                        self.simnode.running = 0
                        return {'CANCELLED'}
                    elif 'Continuity error cannot be removed by adjusting the outflow' in dline[1]:
                        self.report({'ERROR'}, "Mass flow discrepencies cannot be resolved.")
                        logentry('ERROR: Mass flow discrepencies cannot be resolved. This can happen if using fixed velocity on all boundaries and the areas of inflow and outflow boundaries are different.')
                        self.simnode.running = 0
                        return {'CANCELLED'}
                    elif 'not constraint type' in dline[1]:
                        logentry('ERROR: Mesh and material boundaries are out-of-sync. Recreate the mesh')
                        self.report({'ERROR'}, 'Mesh and material boundaries are out-of-sync. Recreate the mesh')
                        self.simnode.running = 0
                        return {'CANCELLED'}
                    elif 'Invalid wall function specification' in dline[1]:
                        logentry('ERROR: Mesh and material boundaries are out-of-sync. Recreate the mesh')
                        self.report({'ERROR'}, 'Mesh and material boundaries are out-of-sync. Recreate the mesh')
                        self.simnode.running = 0
                        return {'CANCELLED'}
                    else:
                        logentry(f'ERROR: {dline[1]}')

            self.runs[-1].kill()
            open("{}".format(os.path.join(frame_coffb, '{}.foam'.format(frame_c))), "w")

            if self.processes > 1:
                if sys.platform == 'linux':
                    Popen(shlex.split("foamExec reconstructPar -case {}".format(frame_coffb))).wait()
                elif sys.platform in ('darwin', 'win32'):
                    Popen('{} run -it --rm -v "{}":/home/openfoam/data dicehub/openfoam:12 "reconstructPar -case data"'.format(docker_path, frame_coffb), shell=True).wait()

            resdict = {'p': 'Pressure', 'U': 'Speed', 'T': 'Temperature', 'Ux': 'X velocity', 'Uy': 'Y velocity', 'Uz': 'Z velocity', 'Q': 'Volumetric flow rate', 'k': 'Turbulent KE', 'epsilon': 'Turbulent dissipation'}

            for oname in svp['flparams']['probes']:
                if os.path.isdir(os.path.join(frame_coffb, 'postProcessing', oname, '0')):
                    probed = os.path.join(frame_coffb, 'postProcessing', oname, '0')

                    if 'p' in os.listdir(probed):
                        if str(frame_c) not in self.o_dict:
                            self.o_dict[str(frame_c)] = {}

                        self.o_dict[str(frame_c)][oname] = {}

                    for f in os.listdir(probed):
                        if f in ('p', 'T', 'k', 'epsilon'):
                            res = []

                            with open(os.path.join(probed, f), 'r') as resfile:
                                for line in resfile.readlines():
                                    if line and line[0] != '#':
                                        res.append(line.split())

                                resarray = array(res)
                                resarray = transpose(resarray)

                            logentry('{} final {} for frame {} at time {} = {}'.format(oname, resdict[f], frame_c, resarray[0][-1], resarray[1:][-1][-1]))
                            self.o_dict[str(frame_c)][oname][f] = float(resarray[1:][-1][-1])

                            for ri, r in enumerate(resarray[1:]):
                                self.reslists.append([str(frame_c), 'Probe', oname, resdict[f], ' '.join(['{:5f}'.format(float(res)) for res in r])])

                        elif f in ('U'):
                            ts = []
                            u_vals = []
                            ux_vals = []
                            uy_vals = []
                            uz_vals = []

                            with open(os.path.join(probed, f), 'r') as resfile:
                                for line in resfile.readlines():
                                    if line and line[0] != '#':
                                        ls = line.split('(')
                                        ts.append(ls[0].strip())
                                        uvec = mathutils.Vector([float(u.strip(')')) for u in ls[1].split()])
                                        u_vals.append(uvec.length)
                                        ux_vals.append(uvec[0])
                                        uy_vals.append(uvec[1])
                                        uz_vals.append(uvec[2])

                            logentry('{} final speed for frame {} at time {} = {}'.format(oname, frame_c, ts[-1], u_vals[-1]))
                            logentry('{} final X velocity for frame {} at time {} = {}'.format(oname, frame_c, ts[-1], ux_vals[-1]))
                            logentry('{} final Y velocity for frame {} at time {} = {}'.format(oname, frame_c, ts[-1], uy_vals[-1]))
                            logentry('{} final Z velocity for frame {} at time {} = {}'.format(oname, frame_c, ts[-1], uz_vals[-1]))
                            self.o_dict[str(frame_c)][oname]['U'] = u_vals[-1]
                            self.o_dict[str(frame_c)][oname]['Ux'] = ux_vals[-1]
                            self.o_dict[str(frame_c)][oname]['Uy'] = uy_vals[-1]
                            self.o_dict[str(frame_c)][oname]['Uz'] = uz_vals[-1]
                            self.reslists.append([str(frame_c), 'Probe', oname, 'Speed', ' '.join(['{:5f}'.format(u) for u in u_vals])])
                            self.reslists.append([str(frame_c), 'Probe', oname, 'X velocity', ' '.join(['{:5f}'.format(u) for u in ux_vals])])
                            self.reslists.append([str(frame_c), 'Probe', oname, 'Y velocity', ' '.join(['{:5f}'.format(u) for u in uy_vals])])
                            self.reslists.append([str(frame_c), 'Probe', oname, 'Z velocity', ' '.join(['{:5f}'.format(u) for u in uz_vals])])

                    self.reslists.append([str(frame_c), 'Timestep', 'Probe', 'Seconds', ' '.join(['{}'.format(f) for f in resarray[0]])])
                    self.simnode['frames'] = [f for f in self.frames]

            for oname in svp['flparams']['s_probes']:
                vfs, times = [], []

                if sys.platform == 'linux':
                    vf_run = Popen(shlex.split('foamExec foamPostProcess -func "triSurfaceVolumetricFlowRate(name={0}, triSurface={0}.stl)" -case {1}'.format(oname, frame_coffb)), stdout=PIPE)
                elif sys.platform in ('darwin', 'win32'):
                    vf_run = Popen(f'{docker_path} run -it --rm -v "{frame_coffb}":/home/openfoam/data dicehub/openfoam:12 "foamPostProcess -func triSurfaceVolumetricFlowRate\\(triSurface="{oname}.stl"\\) -case data"', stdout=PIPE, stderr=PIPE, shell=True)
                    
                if str(frame_c) not in self.o_dict:
                    self.o_dict[str(frame_c)] = {}

                self.o_dict[str(frame_c)][oname] = {}

                for line in vf_run.stdout.readlines()[::-1]:
                    if "U =" in line.decode():
                        vfs.append(line.decode().split()[-1])

                    elif 'Time =' in line.decode():
                        ti = line.decode().split()[-1].strip('s')
                        times.append(ti)

                if vfs and times:
                    logentry('{} final volume flow rate for frame {} at time {} = {}'.format(oname, frame_c, times[0], vfs[0]))

                    if 'Timestep' not in [r[1] for r in self.reslists]:
                        self.reslists.append([str(frame_c), 'Timestep', 'Timestep', 'Seconds', ' '.join(['{}'.format(ti) for ti in times[::-1]])])

                    self.o_dict[str(frame_c)][oname]['Q'] = float(vfs[0])
                    self.reslists.append([str(frame_c), 'Probe', oname, 'Volume flow rate', ' '.join(['{}'.format(vf) for vf in vfs[::-1]])])

            for oname in svp['flparams']['b_probes']:
                if os.path.isdir(os.path.join(frame_coffb, 'postProcessing', oname + '_vf', '0')):
                    probed = os.path.join(frame_coffb, 'postProcessing', oname + '_vf', '0')

                    if 'surfaceFieldValue.dat' in os.listdir(os.path.join(probed)):
                        if str(frame_c) not in self.o_dict:
                            self.o_dict[str(frame_c)] = {}

                        self.o_dict[str(frame_c)][oname] = {}
                        metrics = 'Q'
                        t_res = []
                        q_res = []

                        with open(os.path.join(probed, 'surfaceFieldValue.dat'), 'r') as resfile:
                            for line in resfile.readlines():
                                if line and line[0] != '#':
                                    q_res.append(line.split()[1])
                                    t_res.append(line.split()[0])

                            q_array = array(q_res)
                            logentry('{} final {} for frame {} at time {} = {:.2f}'.format(oname, 'Q', frame_c, t_res[-1], float(q_array[-1])))
                            self.o_dict[str(frame_c)][oname]['Q'] = float(q_array[-1])
                            self.reslists.append([str(frame_c), 'Probe', oname, 'Q', ' '.join(['{:5f}'.format(float(q)) for q in q_array])])

                if os.path.isdir(os.path.join(frame_coffb, 'postProcessing', oname, '0')):
                    probed = os.path.join(frame_coffb, 'postProcessing', oname, '0')

                    if 'surfaceFieldValue.dat' in os.listdir(os.path.join(probed)):
                        if str(frame_c) not in self.o_dict:
                            self.o_dict[str(frame_c)] = {}

                        if oname not in self.o_dict[str(frame_c)]:
                            self.o_dict[str(frame_c)][oname] = {}

                        metrics = ('p', 'Ux', 'Uy', 'Uz', 'U')
                        t_res = []
                        p_res = []
                        ux_res = []
                        uy_res = []
                        uz_res = []
                        u_res = []

                        with open(os.path.join(probed, 'surfaceFieldValue.dat'), 'r') as resfile:
                            for line in resfile.readlines():
                                if '# Time' in line:
                                    p_index = line.split().index('areaAverage(p)') - 1
                                    u_index = line.split().index('areaAverage(U)') - 1
                                if line and line[0] != '#':
                                    t_res.append(line.split()[0])
                                    p_res.append(line.split()[p_index])
                                    ux_res.append(line.split()[u_index].strip('('))
                                    uy_res.append(line.split()[u_index + 1])
                                    uz_res.append(line.split()[u_index + 2].strip(')'))
                                    u_res.append((float(ux_res[-1])**2 + float(uy_res[-1])**2 + float(uz_res[-1])**2)**0.5)

                        for ri, res in enumerate((p_res, ux_res, uy_res, uz_res, u_res)):
                            res_array = array(res)
                            logentry('{} final {} for frame {} at time {} = {:.2f}'.format(oname, resdict[metrics[ri]], frame_c, t_res[-1], float(res_array[-1])))
                            self.o_dict[str(frame_c)][oname][metrics[ri]] = float(res_array[-1])
                            self.reslists.append([str(frame_c), 'Probe', oname, resdict[metrics[ri]], ' '.join(['{:5f}'.format(float(res)) for res in res_array])])

                    if 'Seconds' not in [r[3] for r in self.reslists]:
                        self.reslists.append([str(frame_c), 'Timestep', 'Probe', 'Seconds', ' '.join(['{}'.format(t) for t in t_res])])

            if len(self.runs) < svp['flparams']['end_frame'] - svp['flparams']['start_frame'] + 1:
                self.pb = qtfvprogress(os.path.join(svp['viparams']['newdir'], 'viprogress'), pdll_path, svp['flparams']['et'], str(self.residuals), frame_n)

                with open(self.fpfile, 'w') as fvprogress:
                    if sys.platform == 'linux':
                        if self.processes > 1:
                            self.runs.append(Popen(shlex.split('mpirun --oversubscribe -np {} foamExec {} -parallel -case {}'.format(self.processes,
                                                                                                                                     svp['flparams']['solver'],
                                                                                                                                     frame_noffb)), stderr=PIPE, stdout=fvprogress))
                        else:
                            self.runs.append(Popen(shlex.split('{} {} {} {}'.format('foamExec', svp['flparams']['solver'], "-case", frame_noffb)), stderr=PIPE, stdout=fvprogress))

                    elif sys.platform in ('darwin', 'win32'):
                        if self.processes > 1:
                            self.runs.append(Popen('{} run -it --rm -v {}:/home/openfoam/data dicehub/openfoam:12 "mpirun --oversubscribe -np {} {} -parallel -case data"'.format(docker_path, frame_noffb,
                                                                                                                                                                                      self.processes,
                                                                                                                                                                                      svp['flparams']['solver']), stderr=PIPE, stdout=fvprogress))
                        else:
                            self.runs.append(Popen('{} run -it --rm -v "{}":/home/openfoam/data dicehub/openfoam:12 "{} -case data"'.format(docker_path, frame_noffb, svp['flparams']['solver']), shell=True, stderr=PIPE, stdout=fvprogress))

                return {'PASS_THROUGH'}

            if len(self.frames) > 1:
                if self.o_dict:
                    self.reslists.append(['All', 'Frames', 'Frames', 'Frames', ' '.join(['{}'.format(f) for f in self.frames])])

                    for oname in self.o_dict[str(self.frames[0])]:
                        for param in self.o_dict[str(self.frames[0])][oname]:
                            self.reslists.append(['All', 'Probe', oname, resdict[param], ' '.join(['{:.3f}'.format(self.o_dict[str(f)][oname][param]) for f in self.frames])])

            self.simnode['reslists'] = self.reslists
            self.simnode['frames'] = [f for f in self.frames]
            self.simnode.postsim()

            if self.pv and sys.platform == 'linux':
                Popen(shlex.split("foamExec paraFoam -builtin -case {}".format(frame_coffb)))

            return {'FINISHED'}


class NODE_OT_Au_Rir(bpy.types.Operator):
    bl_idname = "node.rir_sim"
    bl_label = "IR Generator"
    bl_description = "Generates a room impulse response"
    bl_register = True
    bl_undo = False

    def calc_thread(self, room, q_rts, ir_list):
        i = 0

        try:
            room.compute_rir()
            rts = room.measure_rt60(plot=False, decay_db=60)
        except Exception:
            try:
                rts = room.measure_rt60(plot=False, decay_db=30)
            except Exception as e:
                q_rts.put(e)
                return

        try:
            for mi, mrir in enumerate(room.rir):
                for si, srir in enumerate(mrir):
                    ir_list[i] = srir
                    i += 1
        except Exception as e:
            print(e)

        q_rts.put(rts)
        return

    def execute(self, context):
        scene = context.scene
        svp = scene.vi_params
        svp.vi_display = 0

        if viparams(self, scene):
            return {'CANCELLED'}

        if not svp.get('viparams'):
            svp['viparams'] = {}

        if not svp.get('liparams'):
            svp['liparams'] = {}

        clearscene(context, self)

        for o in scene.objects:
            o.vi_params.vi_type_string = ''

        simnode = context.node
        simnode.presim()
        dp = context.evaluated_depsgraph_get()
        empties = [o for o in bpy.data.objects if o.type == 'EMPTY' and o.visible_get()]
        sources = [o for o in empties if o.vi_params.auvi_sl == '0']
        simnode['coptions']['au_sources'] = [s.name for s in sources]
        mics = [o for o in empties if o.vi_params.auvi_sl == '1']
        mic_arrays = [o for o in bpy.data.objects if o.type == 'MESH' and o.material_slots and o.visible_get() and any([o.material_slots[p.material_index].material.vi_params.mattype == '1' for p in o.data.polygons])]

        for o in mic_arrays:
            (o.vi_params['omax'], o.vi_params['omin'], o.vi_params['oave'], o.vi_params['livires']) = ({}, {}, {}, {})

        robs = [o for o in bpy.data.objects if o.type == 'MESH' and o.visible_get() and any([ms.material.vi_params.mattype == '3' for ms in o.material_slots])]
        mats = [mat for mat in bpy.data.materials if mat.vi_params.mattype == '3']
        reslists = []
        frames = [f for f in range(simnode.startframe, simnode.endframe + 1)] if simnode.animated else [scene.frame_current]
        resdict = {}

        for frame in frames:
            amat_abs = {}
            amat_scatts = {}
            amat_dict = {}
            resdict[str(frame)] = {}
            scene.frame_set(frame)

            for mat in mats:
                mvp = mat.vi_params

                if not mvp.am.updated:
                    mvp.am.update()

                if mvp.auvi_abs_class == '0':
                    new_abs = mvp.am.mat_dict['absorption'][mvp.auvi_type_abs][mvp.auvi_mat_abs]['coeffs']

                    if len(new_abs) < 7:
                        new_abs += [new_abs[-1]] * (7 - len(new_abs))

                elif mvp.auvi_abs_flat:
                    new_abs = [mvp.auvi_o1_abs] * 7

                else:
                    new_abs = [mvp.auvi_o1_abs, mvp.auvi_o2_abs, mvp.auvi_o3_abs, mvp.auvi_o4_abs, mvp.auvi_o5_abs, mvp.auvi_o6_abs, mvp.auvi_o7_abs]

                if mvp.auvi_scatt_class == '0':
                    new_scatts = mvp.am.mat_dict['scattering'][mvp.auvi_type_scatt][mvp.auvi_mat_scatt]['coeffs']

                    if len(new_scatts) < 7:
                        new_scatts += [new_scatts[-1]] * (7 - len(new_scatts))

                elif mvp.auvi_scatt_flat:
                    new_scatts = [round(mvp.auvi_o1_scatt, 2)] * 7
                else:
                    new_scatts = [round(mvp.auvi_o1_scatt, 2), round(mvp.auvi_o2_scatt, 2), round(mvp.auvi_o3_scatt, 2), round(mvp.auvi_o4_scatt, 2),
                                  round(mvp.auvi_o5_scatt, 2), round(mvp.auvi_o6_scatt, 2), round(mvp.auvi_o7_scatt, 2)]

                amat_abs = {'coeffs': new_abs, "center_freqs": c_freqs}
                amat_scatts = {'coeffs': new_scatts, "center_freqs": c_freqs}
                amat_dict[mat.name] = pra.Material(energy_absorption=amat_abs, scattering=amat_scatts)

            for rob in robs:
                walls = []
                mic_names = []
                room_bm = bmesh.new()
                room_bm.from_object(rob, dp)
                room_bm.transform(rob.matrix_world)
                bmesh.ops.triangulate(room_bm, faces=room_bm.faces)

                for face in room_bm.faces:
                    mat = rob.material_slots[face.material_index].material

                    if mat.name in amat_dict:
                        poly_vecs = array([v.co for v in face.verts]).T
                        walls.append(
                            pra.wall_factory(
                                poly_vecs,
                                amat_dict[mat.name].energy_absorption["coeffs"],
                                amat_dict[mat.name].scattering["coeffs"],
                            )
                        )

                room_bm.free()

                room = (
                    pra.Room(
                        walls,
                        fs=16000,
                        max_order=simnode.max_order,
                        ray_tracing=pra_rt,
                        air_absorption=False,
                        use_rand_ism=True,
                        max_rand_disp=0.1
                    )
                )

                if pra_rt:
                    room.set_ray_tracing(n_rays=simnode.rt_rays, time_thres=10.0, receiver_radius=simnode.r_radius,
                                         hist_bin_size=0.004, energy_thres=1e-08)

                for source in sources:
                    if room.is_inside(source.location[:]):
                        room.add_source(source.location[:])

                if not room.sources:
                    self.report({'ERROR'}, f'No visible sources inside room {rob.name}')
                    return {'CANCELLED'}

                for mic in mics:
                    if room.is_inside(mic.location[:]):
                        room.add_microphone(mic.location[:])
                        mic_names.append(mic.name)

                for mic_a in mic_arrays:
                    mic_bm = bmesh.new()
                    mic_bm.from_mesh(mic_a.data)
                    mic_bm.transform(mic_a.matrix_world)

                    if not mic_bm.faces.layers.int.get('cindex'):
                        mic_bm.faces.layers.int.new('cindex')

                    bm_ir = mic_bm.faces.layers.int['cindex']
                    mic_bm.faces.ensure_lookup_table()

                    for fi, f in enumerate(mic_bm.faces):
                        if mic_a.material_slots[f.material_index].material.vi_params.mattype == '1':
                            if room.is_inside(f.calc_center_median()[:]):
                                if mic_a.vi_params.vi_type_string != 'LiVi Calc':
                                    mic_a.vi_params.vi_type_string = 'LiVi Calc'

                                f[bm_ir] = 1
                                room = room.add_microphone(f.calc_center_median()[:])
                                mic_names.append(f'{mic_a.name}-{f.index}')
                            else:
                                f[bm_ir] = 0
                        else:
                            f[bm_ir] = 0

                    mic_bm.transform(mic_a.matrix_world.inverted())
                    mic_bm.to_mesh(mic_a.data)
                    mic_bm.free()

                if not room.n_mics:
                    self.report({'ERROR'}, 'No visible listeners inside the room')
                    return {'CANCELLED'}

                Lsf = array([62.9, 62.9, 59.2, 53.2, 47.2, 41.2, 35.2])
                Lsf = Lsf if room.volume < 250 else Lsf + 10
                octave = pra.acoustics.OctaveBandsFactory(base_frequency=125, fs=16000, n_fft=512)

                try:
                    if sys.platform != 'win32':
                        auvi_cancel = cancel_window(os.path.join(scene.vi_params['viparams']['newdir'], pdll_path, 'viprogress'), 'Calculating RTs')
                        p_manager = multiprocessing.Manager()
                        ir_list = p_manager.list()

                        for m in range(room.n_mics * len(sources)):
                            ir_list.append([])

                        q_rts = multiprocessing.Queue()
                        t1 = multiprocessing.Process(target=self.calc_thread, args=(room, q_rts, ir_list))
                        t1.daemon = True
                        t1.start()

                        while t1.is_alive():
                            if auvi_cancel.poll() is not None:
                                t1.kill()
                                return {'CANCELLED'}
                            sleep(0.1)

                        if auvi_cancel.poll() is None:
                            auvi_cancel.kill()

                        t1.join()
                        rts = q_rts.get()

                        if isinstance(rts, Exception):
                            logentry(f'Current settings cannot produce valid RTs. Try increasing reciver radius or RT rays {str(rts)}')
                            self.report({'ERROR'}, 'Current settings cannot produce valid RTs. Try increasing reciver radius or RT rays')
                            return {'CANCELLED'}

                        rirs = ir_list
                        t1.kill()

                    else:
                        room.compute_rir()
                        rirs = []

                        for mrir in room.rir:
                            for srir in mrir:
                                rirs.append(srir)
                        try:
                            rts = room.measure_rt60(plot=False, decay_db=60)
                        except Exception:
                            logentry("Can't get a reliable 60dB reduction. Extrapolating from a 30dB reduction")
                            rts = room.measure_rt60(plot=False, decay_db=30)

                    gc.collect()

                except Exception as e:
                    logentry(str(e))
                    self.report({'ERROR'}, str(e))
                    return {'CANCELLED'}

                i = 0
                for mi, mic_rt in enumerate(rts):
                    if mi < len(mics):
                        for si, source_rt in enumerate(mic_rt):
                            reslists.append([str(frame), 'Probe', f'{mic_names[mi]} - {sources[si].name}', 'Seconds', ' '.join([str(s / 16000) for s in range(len(rirs[i]))])])
                            reslists.append([str(frame), 'Probe', f'{mic_names[mi]} - {sources[si].name}', 'RIR', ' '.join(rirs[i].astype('str'))])
                            reslists.append([str(frame), 'Probe', f'{mic_names[mi]} - {sources[si].name}', 'RT', f'{source_rt:.3f}'])
                            resdict[str(frame)][f'{mic_names[mi]} - {sources[si].name}'] = f'{source_rt:.3f}'
                            i += 1

                for si, source in enumerate(sources):
                    fi = len(mics)

                    for mic_a in mic_arrays:
                        ovp = mic_a.vi_params
                        mic_bm = bmesh.new()
                        mic_bm.from_mesh(mic_a.data)

                        if not mic_bm.faces.layers.float.get(f'{source.name}_rt{frame}'):
                            mic_bm.faces.layers.float.new(f'{source.name}_rt{frame}')
                        if not mic_bm.faces.layers.float.get(f'{source.name}_vol{frame}'):
                            mic_bm.faces.layers.float.new(f'{source.name}_vol{frame}')
                        if not mic_bm.faces.layers.float.get(f'{source.name}_sti{frame}'):
                            mic_bm.faces.layers.float.new(f'{source.name}_sti{frame}')

                        bm_rtres = mic_bm.faces.layers.float[f'{source.name}_rt{frame}']
                        bm_volres = mic_bm.faces.layers.float[f'{source.name}_vol{frame}']
                        bm_stires = mic_bm.faces.layers.float[f'{source.name}_sti{frame}']
                        bm_ir = mic_bm.faces.layers.int['cindex']

                        for face in mic_bm.faces:
                            if face[bm_ir]:
                                try:
                                    face[bm_rtres] = rts[fi][si]
                                    face[bm_volres] = 10 * log(nsum(square(rirs[fi * len(sources) + si]) / 16000) / 6E-07, 10)
                                    face[bm_stires] = rir2sti(rirs[fi * len(sources) + si], room.volume, source.location, mic_a.matrix_world @ face.calc_center_bounds(), octave, 'male', Lsf)

                                except Exception as e:
                                    print(e)

                                fi += 1

                        res_rt = [face[bm_rtres] for face in mic_bm.faces if face[bm_ir]]
                        res_vol = [face[bm_volres] for face in mic_bm.faces if face[bm_ir]]
                        res_sti = [face[bm_stires] for face in mic_bm.faces if face[bm_ir]]

                        if res_rt:
                            ovp['omax'][f'rt{frame}'] = max(res_rt) if not ovp['omax'].get(f'rt{frame}') or max(res_rt) > ovp['omax'][f'rt{frame}'] else ovp['omax'][f'rt{frame}']
                            ovp['oave'][f'rt{frame}'] = sum(res_rt) / len(res_rt)
                            ovp['omin'][f'rt{frame}'] = min(res_rt) if not ovp['omin'].get(f'rt{frame}') or min(res_rt) < ovp['omin'][f'rt{frame}'] else ovp['omin'][f'rt{frame}']
                            ovp['omax'][f'vol{frame}'] = max(res_vol) if not ovp['omax'].get(f'vol{frame}') or max(res_vol) > ovp['omax'][f'vol{frame}'] else ovp['omax'][f'vol{frame}']
                            ovp['oave'][f'vol{frame}'] = sum(res_vol) / len(res_vol)
                            ovp['omin'][f'vol{frame}'] = min(res_vol) if not ovp['omin'].get(f'vol{frame}') or min(res_vol) < ovp['omin'][f'vol{frame}'] else ovp['omin'][f'vol{frame}']
                            ovp['omax'][f'sti{frame}'] = max(res_sti) if not ovp['omax'].get(f'sti{frame}') or max(res_sti) > ovp['omax'][f'sti{frame}'] else ovp['omax'][f'sti{frame}']
                            ovp['oave'][f'sti{frame}'] = sum(res_sti) / len(res_sti)
                            ovp['omin'][f'sti{frame}'] = min(res_sti) if not ovp['omin'].get(f'sti{frame}') or min(res_sti) < ovp['omin'][f'sti{frame}'] else ovp['omin'][f'sti{frame}']
                            ovp['livires'][f'rt{frame}'] = res_rt
                            ovp['livires'][f'vol{frame}'] = res_vol
                            ovp['livires'][f'sti{frame}'] = res_sti
                            mic_bm.to_mesh(mic_a.data)
                        else:
                            self.report({'WARNING'}, f'No results on sensor mesh {mic_a.name}')

                        mic_bm.free()

                print("Room volume:", f'{room.get_volume():.2f}')
                print("RT60 (Simulated):", rts[0, 0])
                print("RT60 (Sabine):", room.rt60_theory(formula='sabine'))
                print("RT60 (Eyring):", room.rt60_theory(formula='eyring'))

        if len(frames) > 1:
            reslists.append(['All', 'Frames', 'Frames', 'Frames', ' '.join(['{}'.format(f) for f in frames])])
            p_dict = {}

            for frame in frames:
                for rl in reslists:
                    if rl[0] == str(frame):
                        if rl[3] == 'RT':
                            if frame == frames[0]:
                                p_dict[rl[2]] = [rl[4]]
                            else:
                                p_dict[rl[2]].append(rl[4])

            for pd in p_dict:
                reslists.append(['All', 'Probe', pd, 'RT', ' '.join(p_dict[pd])])

        if mic_arrays:
            svp['liparams']['unit'] = 'RT60 (s)'
            svp['liparams']['type'] = 'VI Acoustics'
            svp['liparams']['fs'] = frames[0]
            svp['liparams']['fe'] = frames[-1]
            svp['liparams']['cp'] = '0'
            svp['liparams']['offset'] = 0.01
            svp['viparams']['resnode'], svp['viparams']['restree'] = simnode.name, simnode.id_data.name
            svp['viparams']['vidisp'] = 'rt'

        if not simnode.get('goptions'):
            simnode['goptions'] = {}

        simnode['goptions']['offset'] = 0.01
        simnode['reslists'] = reslists
        simnode['resdict'] = resdict
        simnode['frames'] = frames
        svp['liparams']['sources'] = simnode['coptions']['au_sources']
        simnode.postsim()
        scene.frame_set(frames[0])
        return {'FINISHED'}


class NODE_OT_WavSelect(NODE_OT_FileSelect):
    bl_idname = "node.wavselect"
    bl_label = "Select WAV file"
    bl_description = "Select the WAV to be convolved"
    filename_ext = ".WAV;.wav;"
    filter_glob: bpy.props.StringProperty(default="*.WAV;*.wav;", options={'HIDDEN'})
    nodeprop = 'wavname'
    filepath: bpy.props.StringProperty(subtype='FILE_PATH', options={'HIDDEN', 'SKIP_SAVE'})
    fextlist = ("WAV", "wav")


class NODE_OT_Au_Conv(bpy.types.Operator):
    bl_idname = "node.auvi_conv"
    bl_label = "Convolve"
    bl_description = "Generates a convolved sound based on an impulse response"
    bl_register = True
    bl_undo = False

    def execute(self, context):
        convnode = context.node
        ir_node = convnode.inputs[0].links[0].from_node
        fs, audio = wavfile.read(convnode.wavname)
        odt = audio.dtype

        if len(audio.shape) > 1:
            audio = nmean(audio, axis=1).astype(odt)

        if odt == int16:
            audio = audio.astype(float32) / 32768
        elif odt == int32:
            audio = (audio.astype(float32) / 2147483647)
        elif odt == float64:
            audio.astype(float32)

        if fs != 16000:
            samples = int(16000 * len(audio) / fs)
            audio = signal.resample(audio, samples)
            fs = 16000

        for rl in ir_node['reslists']:
            if f'{rl[0]} - {rl[2]}' == convnode.rir and rl[3] == 'RIR':
                ir = array([float(s) for s in rl[4].split()]).astype(float32, order='C')
                break

        convnode['convolved_audio'] = []
        convnode['convolved_audio'] = signal.fftconvolve(audio, ir, mode="full").astype(float32, order='C')
        convnode.postsim()
        return {'FINISHED'}


class NODE_OT_Au_Play(bpy.types.Operator):
    bl_idname = "node.auvi_play"
    bl_label = "Play"
    bl_description = "Plays a wav file"
    bl_register = True
    bl_undo = False

    def execute(self, context):
        scene = context.scene
        self.convnode = context.node
        wm = context.window_manager
        device = aud.Device()

        if not os.path.isfile(self.convnode.wavname):
            self.report({'ERROR'}, 'No file found')
            return {'CANCELLED'}

        sound = aud.Sound(self.convnode.wavname)
        self.handle = device.play(sound)
        self.convnode.play_o = True
        self._timer = wm.event_timer_add(0.1, window=context.window)
        wm.modal_handler_add(self)
        return {'RUNNING_MODAL'}

    def modal(self, context, event):
        scene = context.scene

        if not self.convnode.play_o:
            self.handle.stop()
            return {'FINISHED'}

        if not self.handle.status:
            self.convnode.play_o = False
            return {'FINISHED'}

        return {'PASS_THROUGH'}


class NODE_OT_Au_Stop(bpy.types.Operator):
    bl_idname = "node.auvi_stop"
    bl_label = "Stop"
    bl_description = "Stop playing a wav file"
    bl_register = True
    bl_undo = False

    def execute(self, context):
        scene = context.scene
        convnode = context.node
        convnode.play_o = False
        convnode.play_c = False
        return {'FINISHED'}


class NODE_OT_Au_PlayC(bpy.types.Operator):
    bl_idname = "node.auvi_playc"
    bl_label = "Play convolved file"
    bl_description = "Plays a convolved wav file"
    bl_register = True
    bl_undo = False

    def execute(self, context):
        scene = context.scene
        self.convnode = context.node
        wm = context.window_manager
        device = aud.Device()
        sound = aud.Sound.buffer(array(self.convnode['convolved_audio']).astype(float32), 16000)
        sound = sound.resample(48000, True)
        self.handle = device.play(sound)
        self.convnode.play_c = True
        self._timer = wm.event_timer_add(0.1, window=context.window)
        wm.modal_handler_add(self)
        return {'RUNNING_MODAL'}

    def modal(self, context, event):
        scene = context.scene

        if not self.convnode.play_c:
            self.handle.stop()
            return {'FINISHED'}

        if not self.handle.status:
            self.convnode.play_c = False

            return {'FINISHED'}

        return {'PASS_THROUGH'}


class NODE_OT_Au_Save(bpy.types.Operator, ExportHelper):
    bl_idname = "node.auvi_save"
    bl_label = "Save"
    bl_description = "Save a convolved wav file"
    bl_register = True
    bl_undo = False
    filename_ext = ".wav"

    filter_glob: bpy.props.StringProperty(default="*.wav", options={'HIDDEN'}, maxlen=255)

    def invoke(self, context, event):
        self.node = context.node
        self.filepath = f"{self.node.rir}.wav"
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}

    def execute(self, context):
        sound = aud.Sound.buffer(array(self.node['convolved_audio']).astype(float32), 16000)
        sound.write(self.filepath, 16000, 1, 0, 0, 0, 16, 256)
        return {'FINISHED'}

# class ADDON_OT_PyInstall(bpy.types.Operator):
#     bl_idname = "addon.pyimport"
#     bl_label = "Install Python dependencies"
#     bl_description = "Installs matplotlib, PyQt6, kivy and netgen"
#
#     def execute(self, context):
#         if not os.path.isdir(os.path.join(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))), 'Python', sys.platform, 'pip')):
#             gp_cmd = '{} {} --target {}'.format(sys.executable, os.path.join(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))), 'get-pip.py'),
#                                                 os.path.join(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))), 'Python', sys.platform))
#             Popen(shlex.split(gp_cmd))
#
#         if not os.path.isdir(os.path.join(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))), 'Python', sys.platform, 'kivy')):
#             kivy_cmd = '{} -m pip install kivy --target {}'.format(sys.executable,
#                                                                    os.path.join(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))), 'Python', sys.platform))
#             Popen(shlex.split(kivy_cmd))
#
#         if not os.path.isdir(os.path.join(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))), 'Python', sys.platform, 'PyQt6')):
#             pyqt_cmd = '{} -m pip install PyQt6 --target {}'.format(sys.executable,
#                                                                     os.path.join(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))), 'Python', sys.platform))
#             Popen(shlex.split(pyqt_cmd))
#
#         if not os.path.isdir(os.path.join(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))), 'Python', sys.platform, 'matplotlib')):
#             mp_cmd = '{} -m pip install matplotlib --target {}'.format(sys.executable,
#                                                                        os.path.join(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))), 'Python', sys.platform))
#             Popen(shlex.split(mp_cmd))
#
#         if not os.path.isdir(os.path.join(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))), 'Python', sys.platform, 'netgen')):
#             ng_cmd = '{} -m pip install netgen --target {}'.format(sys.executable,
#                                                                    os.path.join(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))), 'Python', sys.platform))
#             Popen(shlex.split(ng_cmd))
#
#         if not os.path.isdir(os.path.join(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))), 'Python', sys.platform, 'netgen')):
#             ng_cmd = '{} -m pip install pyroomacoustics --target {}'.format(sys.executable,
#                                                                    os.path.join(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))), 'Python', sys.platform))
#             Popen(shlex.split(ng_cmd))
#
#         return {'FINISHED'}
