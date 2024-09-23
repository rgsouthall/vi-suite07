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


import bpy, blf, mathutils, datetime, os, inspect, gpu, bmesh
from gpu_extras.batch import batch_for_shader
from mathutils import Vector
from bpy_extras import view3d_utils
from .vi_func import ret_vp_loc, viewdesc, draw_index, draw_time, blf_props, retcols, retdp
from .vi_func import ret_res_vals, draw_index_distance, selobj, mp2im, move_obs
from .vi_func import logentry, move_to_coll, cmap, retvpvloc, objmode, skframe, clearscene
from .vi_func import solarPosition, solarRiseSet, create_coll, create_empty_coll, compass, joinobj, sunpath, sunpath1
from .livi_func import setscenelivivals, res_interpolate, res_direction
from .auvi_func import setsceneauvivals 
# from .livi_export import spfc
from .vi_dicts import res2unit, unit2res
from . import livi_export
from .vi_svg import vi_info
from math import pi, log10, atan2, sin, cos
from numpy import array, repeat, logspace, multiply, digitize, frombuffer, ubyte, float32, int8
from numpy import min as nmin
from numpy import max as nmax
from numpy import sum as nsum
from numpy import log10 as nlog10
from numpy import append as nappend
from xml.dom.minidom import parseString
# from bpy.app.handlers import persistent
from PySide6.QtGui import QImage, QPdfWriter, QPagedPaintDevice, QPainter, QPageSize
from PySide6.QtPrintSupport import QPrinter
from PySide6.QtSvg import QSvgRenderer
from PySide6.QtCore import QSizeF, QMarginsF
from PySide6.QtCore import QSize
from PySide6.QtWidgets import QApplication


try:
    import matplotlib
    matplotlib.use('qtagg', force=True)
    import matplotlib.pyplot as plt
    plt.ion()
    plt.ioff()
    import matplotlib.cm as mcm
    import matplotlib.colors as mcolors
    from matplotlib.patches import Rectangle
    from matplotlib.collections import PatchCollection
    mp = 1
except Exception as e:
    print("No matplotlib: {}".format(e))
    mp = 0

kfsa = array([0.02391, 0.02377, 0.02341, 0.02738, 0.02933, 0.03496, 0.04787, 0.05180, 0.13552])
kfact = array([0.9981, 0.9811, 0.9361, 0.8627, 0.7631, 0.6403, 0.4981, 0.3407, 0.1294])


def ret_dcoords(context):
    return (context.area.regions[0].height + context.area.regions[1].height, context.area.regions[4].width, context.area.regions[-1].height)


def script_update(self, context):
    print('script_update')
    svp = context.scene.vi_params

    if svp.vi_res_process == '2' and svp.script_file in bpy.data.texts:
        script = bpy.data.texts[svp.script_file]
        exec(script.as_string())

    leg_update(self, context)


def col_update(self, context):
    cmap(context.scene.vi_params)


def leg_update(self, context):
    scene = context.scene
    svp = scene.vi_params
    params = 'liparams'
    type_strings = ('LiVi Res', 'LiVi Calc')
    disp_menu = svp.li_disp_menu
    legmm = leg_min_max(svp)
    frames = range(svp[params]['fs'], svp[params]['fe'] + 1)
    obs = [o for o in scene.objects if o.vi_params.vi_type_string == type_strings[0]]
    cobs = [o for o in scene.objects if o.vi_params.vi_type_string == type_strings[1]]
    increment = 1/svp.vi_leg_levels

    if svp.vi_leg_scale == '0':
        bins = array([increment * i for i in range(1, svp.vi_leg_levels + 1)])

    elif svp.vi_leg_scale == '1':
        slices = logspace(0, 2, svp.vi_leg_levels + 1, True)
        bins = array([(slices[i] - increment * (svp.vi_leg_levels - i))/100 for i in range(svp.vi_leg_levels + 1)])
        bins = array([1 - log10(i)/log10(svp.vi_leg_levels + 1) for i in range(1, svp.vi_leg_levels + 2)][::-1])
        bins = bins[1:-1]

    for o in obs:
        selobj(context.view_layer, o)
        bm = bmesh.new()
        bm.from_mesh(o.data)
        cmap(self)
        
        if svp.vi_disp_process == '2':
            if o.name[:-3] in [cob.name for cob in cobs]:
                cob = cobs[[cob.name for cob in cobs].index(o.name[:-3])]                
                res_interpolate(scene, context.evaluated_depsgraph_get(), cob, o, plt, svp[params]['offset'])
        
        elif svp.vi_disp_process == '3':
            if o.name[:-3] in [cob.name for cob in cobs]:
                cob = cobs[[cob.name for cob in cobs].index(o.name[:-3])] 
                res_direction(scene, cob, o, svp[params]['offset'])

            if len(o.material_slots) != svp.vi_leg_levels:
                for matname in ['{}#{}'.format('vi-suite', i) for i in range(0, svp.vi_leg_levels)]:
                    if bpy.data.materials[matname] not in o.data.materials[:]:
                        bpy.ops.object.material_slot_add()
                        o.material_slots[-1].material = bpy.data.materials[matname]
                while len(o.material_slots) > svp.vi_leg_levels:
                    bpy.ops.object.material_slot_remove()
            
            for f, frame in enumerate(frames):
                if disp_menu == 'aga1v':
                    res_name = 'aga{}v{}'.format(svp.vi_views, frame)
                elif disp_menu == 'ago1v':
                    res_name = 'ago{}v{}'.format(svp.vi_views, frame)
                elif disp_menu == 'rt':
                    res_name = f'{svp.au_sources}_rt{frame}'
                else:
                    res_name = '{}{}'.format(disp_menu, frame)

                if bm.faces.layers.float.get(res_name):
                    vires = bm.faces.layers.float[res_name]
                    ovals = array([f[vires] for f in bm.faces])
                elif bm.verts.layers.float.get(res_name):
                    vires = bm.verts.layers.float[res_name]
                    ovals = array([sum([vert[vires] for vert in f.verts])/len(f.verts) for f in bm.faces])

                ovals = array(ret_res_vals(svp, ovals))

                if legmm[1] > legmm[0]:
                    vals = ovals - legmm[0]
                    vals = vals/(legmm[1] - legmm[0])
                else:
                    vals = array([legmm[1] for f in bm.faces])

                if svp.vi_res_process == '2' and svp.script_file:
                    nmatis = array(ovals).astype(int8)
                else:
                    nmatis = digitize(vals, bins, right=False).clip(0, svp.vi_leg_levels - 1)

                if len(frames) == 1:
                    o.data.polygons.foreach_set('material_index', nmatis)
                    o.data.update()

                elif len(frames) > 1:
                    for fi, fc in enumerate(o.data.animation_data.action.fcurves):
                        fc.keyframe_points[f].co = frame, nmatis[fi]

        else:
            if len(o.material_slots) != svp.vi_leg_levels:
                for matname in ['{}#{}'.format('vi-suite', i) for i in range(0, svp.vi_leg_levels)]:
                    if bpy.data.materials[matname] not in o.data.materials[:]:
                        bpy.ops.object.material_slot_add()
                        o.material_slots[-1].material = bpy.data.materials[matname]
                while len(o.material_slots) > svp.vi_leg_levels:
                    bpy.ops.object.material_slot_remove()

            for f, frame in enumerate(frames):
                if disp_menu == 'aga1v':
                    res_name = 'aga{}v{}'.format(svp.vi_views, frame)
                elif disp_menu == 'ago1v':
                    res_name = 'ago{}v{}'.format(svp.vi_views, frame)
                elif disp_menu == 'rt':
                    res_name = f'{svp.au_sources}_rt{frame}'
                else:
                    res_name = '{}{}'.format(disp_menu, frame)

                if bm.faces.layers.float.get(res_name):
                    vires = bm.faces.layers.float[res_name]
                    ovals = array([f[vires] for f in bm.faces])
                elif bm.verts.layers.float.get(res_name):
                    vires = bm.verts.layers.float[res_name]
                    ovals = array([sum([vert[vires] for vert in f.verts])/len(f.verts) for f in bm.faces])

                ovals = array(ret_res_vals(svp, ovals))

                if legmm[1] > legmm[0]:
                    vals = ovals - legmm[0]
                    vals = vals/(legmm[1] - legmm[0])
                else:
                    vals = array([legmm[1] for f in bm.faces])

                if svp.vi_res_process == '2' and svp.script_file:
                    nmatis = array(ovals).astype(int8)
                else:
                    nmatis = digitize(vals, bins, right=False).clip(0, svp.vi_leg_levels - 1)

                if len(frames) == 1:
                    o.data.polygons.foreach_set('material_index', nmatis)
                    o.data.update()

                elif len(frames) > 1:
                    for fi, fc in enumerate(o.data.animation_data.action.fcurves):
                        fc.keyframe_points[f].co = frame, nmatis[fi]
        bm.free()
    
    try:
        context.space_data.region_3d.view_location[2] += 0.0001
    except Exception as e:
        pass


def leg_min_max(svp):
    try:
        if svp.vi_res_process == '2' and 'resmod' in bpy.app.driver_namespace.keys():
            return bpy.app.driver_namespace['resmod']([svp.vi_leg_min, svp.vi_leg_max])
        elif svp.vi_res_process == '1' and svp.vi_res_mod:
            return (eval('{}{}'.format(svp.vi_leg_min, svp.vi_res_mod)), eval('{}{}'.format(svp.vi_leg_max, svp.vi_res_mod)))
        else:
            return (svp.vi_leg_min, svp.vi_leg_max)
    except Exception as e:
        logentry('Error setting legend values: {}'.format(e))
        return (svp.vi_leg_min, svp.vi_leg_max)


def e_update(self, context):
    scene = context.scene
    svp = scene.vi_params
    maxo, mino = svp.vi_leg_max, svp.vi_leg_min
    odiff = svp.vi_leg_max - svp.vi_leg_min
    
    if context.active_object and context.active_object.mode == 'EDIT':
        return

    if odiff:
        for frame in range(svp['liparams']['fs'], svp['liparams']['fe'] + 1):
            if svp.li_disp_menu == 'aga1v':
                res_name = 'aga{}v{}'.format(svp.vi_views, frame)
            elif svp.li_disp_menu == 'ago1v':
                res_name = 'ago{}v{}'.format(svp.vi_views, frame)
            elif svp.li_disp_menu == 'rt':
                res_name = f'{svp.au_sources}_rt{frame}'
            else:
                res_name = '{}{}'.format(svp.li_disp_menu, frame)

            for o in [obj for obj in bpy.data.objects if obj.vi_params.vi_type_string == 'LiVi Res' and obj.data.shape_keys and str(frame) in [sk.name for sk in obj.data.shape_keys.key_blocks]]:
                ovp = o.vi_params
                bm = bmesh.new()
                bm.from_mesh(o.data)
                bm.transform(o.matrix_world)
                skb = bm.verts.layers.shape['Basis']
                skf = bm.verts.layers.shape[str(frame)]

                if str(frame) in ovp['omax']:
                    if bm.faces.layers.float.get(res_name):
                        extrude = bm.faces.layers.int['extrude']

                        res = bm.faces.layers.float[res_name]  # if context.scene['cp'] == '0' else bm.verts.layers.float['res{}'.format(frame)]
                        faces = [f for f in bm.faces if f[extrude]]
                        fnorms = array([f.normal.normalized() for f in faces]).T
                        fres = array([f[res] for f in faces])
                        extrudes = (0.1 * svp.vi_disp_3dlevel * (nlog10(maxo * (fres + 1 - mino)/odiff)) * fnorms).T if svp.vi_leg_scale == '1' else \
                            multiply(fnorms, svp.vi_disp_3dlevel * ((fres - mino)/odiff)).T

                        for f, face in enumerate(faces):
                            for v in face.verts:
                                v[skf] = v[skb] + mathutils.Vector(extrudes[f])

                    elif bm.verts.layers.float.get(res_name):
                        res = bm.verts.layers.float[res_name]
                        vnorms = array([v.normal.normalized() for v in bm.verts]).T
                        vres = array([v[res] for v in bm.verts])
                        extrudes = multiply(vnorms, svp.vi_disp_3dlevel * ((vres-mino)/odiff)).T if svp.vi_leg_scale == '0' else \
                            [0.1 * svp.vi_disp_3dlevel * (log10(maxo * (v[res] + 1 - mino)/odiff)) * v.normal.normalized() for v in bm.verts]

                        for v, vert in enumerate(bm.verts):
                            vert[skf] = vert[skb] + mathutils.Vector(extrudes[v])

                bm.transform(o.matrix_world.inverted())
                bm.to_mesh(o.data)
                bm.free()


def t_update(self, context):
    for o in [o for o in context.scene.objects if o.type == 'MESH' and 'lightarray' not in o.name and not o.hide_viewport and o.vi_params.vi_type_string == 'LiVi Res']:
        o.show_transparent = 1
    for mat in [bpy.data.materials['{}#{}'.format('vi-suite', index)] for index in range(context.scene.vi_params.vi_leg_levels)]:
        mat.blend_method = 'BLEND'
        mat.diffuse_color[3] = self.id_data.vi_params.vi_disp_trans
    
    cmap(self)


def w_update(self, context):
    o = context.active_object
    if o and o.type == 'MESH':
        (o.show_wire, o.show_all_edges) = (1, 1) if context.scene.vi_params.vi_disp_wire else (0, 0)


def livires_update(self, context):
    setscenelivivals(context.scene)

    for o in [o for o in bpy.data.objects if o.vi_params.vi_type_string == 'LiVi Res']:
        o.vi_params.lividisplay(context.scene)

    e_update(self, context)


def auvires_update(self, context):
    setsceneauvivals(context.scene)

    for o in [o for o in bpy.data.objects if o.vi_params.vi_type_string == 'AuVi Res']:
        o.vi_params.lividisplay(context.scene)

    e_update(self, context)

def rendview(i):
    for scrn in bpy.data.screens:
        for area in scrn.areas:
            if area.type == 'VIEW_3D':
                for space in area.spaces:
                    if space.type == 'VIEW_3D':
                        space.clip_start = 0.1
                        bpy.context.scene['cs'] = space.clip_start


def li_display(context, disp_op, simnode):
    if not [o for o in bpy.data.objects if o.vi_params.vi_type_string == 'LiVi Calc']:
        return 'CANCELLED'

    scene, obreslist, obcalclist = context.scene, [], []
    dp = context.evaluated_depsgraph_get()
    svp = scene.vi_params
    svp.li_disp_menu = unit2res[svp['liparams']['unit']]
    setscenelivivals(scene)

    (rcol, mtype) = ('hot', 'livi') if 'LiVi' in simnode.bl_label else ('grey', 'shad')

    for geo in context.view_layer.objects:
        context.view_layer.objects.active = geo

        if getattr(geo, 'mode') != 'OBJECT':
            bpy.ops.object.mode_set(mode='OBJECT')

    bpy.ops.object.select_all(action='DESELECT')

    if not bpy.app.handlers.frame_change_post:
        bpy.app.handlers.frame_change_post.append(livi_export.cyfc1)

    for o in context.view_layer.objects:
        if o.type == "MESH" and o.vi_params.vi_type_string == 'LiVi Calc' and o.visible_get():
            bpy.ops.object.select_all(action='DESELECT')
            obcalclist.append(o)

    scene.frame_set(svp['liparams']['fs'])
    context.view_layer.objects.active = None
    cmap(svp)

    for i, o in enumerate(obcalclist):
        ovp = o.vi_params
        
        if svp.vi_disp_process == "2":
            mesh = bpy.data.meshes.new(f'{o.name}_res')
            ores = bpy.data.objects.new(f'{o.name}_res', mesh)
            ores.name, ores.show_wire, ores.show_all_edges, ores.display_type, orvp, ores.vi_params.vi_type_string = o.name+"res", 1, 1, 'SOLID', ores.vi_params, 'LiVi Res'
            move_to_coll(context, 'LiVi Results', ores)
            context.view_layer.layer_collection.children['LiVi Results'].exclude = 0
            context.view_layer.objects.active = ores
            ores.visible_diffuse, ores.visible_glossy, ores.visible_transmission, ores.visible_volume_scatter, ores.visible_shadow = 0, 0, 0, 0, 0
            ores.vi_params.vi_type_string == 'LiVi Res'
            orvp['omax'], orvp['omin'], orvp['oave'] = ovp['omax'], ovp['omin'], ovp['oave']
            selobj(context.view_layer, ores)
            res_interpolate(scene, dp, o, ores, plt, simnode['goptions']['offset'])
            ores.vi_params.lividisplay(scene)

        elif svp.vi_disp_process == "3":
            mesh = bpy.data.meshes.new(f'{o.name}_res')
            ores = bpy.data.objects.new(f'{o.name}_res', mesh)
            ores.name, ores.show_wire, ores.show_all_edges, ores.display_type, orvp, ores.vi_params.vi_type_string = o.name+"res", 1, 1, 'SOLID', ores.vi_params, 'LiVi Res'
            move_to_coll(context, 'LiVi Results', ores)
            context.view_layer.layer_collection.children['LiVi Results'].exclude = 0
            selobj(context.view_layer, ores)
            res_direction(scene, o, ores, simnode['goptions']['offset'])
            orvp['omax'], orvp['oave'], orvp['omin'] = ovp['omax'], ovp['omin'], ovp['oave']
            orvp.lividisplay(scene)
            obreslist.append(ores)

        else:
            bm = bmesh.new()
            bm.from_object(o, dp)

            if svp['liparams']['cp'] == '0':
                cindex = bm.faces.layers.int['cindex']

                for f in [f for f in bm.faces if f[cindex] < 1]:
                    bm.faces.remove(f)

                [bm.verts.remove(v) for v in bm.verts if not v.link_faces]

            elif svp['liparams']['cp'] == '1':
                cindex = bm.verts.layers.int['cindex']

                for v in [v for v in bm.verts if v[cindex] < 1]:
                    bm.verts.remove(v)
                for v in bm.verts:
                    v.select = True

            while bm.verts.layers.shape:
                bm.verts.layers.shape.remove(bm.verts.layers.shape[-1])

            for v in bm.verts:
                v.co += mathutils.Vector((nsum([f.normal for f in v.link_faces], axis=0))).normalized() * simnode['goptions']['offset']

            selobj(context.view_layer, o)
            bpy.ops.object.duplicate()

            for face in bm.faces:
                face.select = True

            if not context.active_object:
                disp_op.report({'ERROR'}, "No display object. If in local view switch to global view and/or re-export the geometry")
                return 'CANCELLED'

            ores = context.active_object
            ores.name, ores.show_wire, ores.show_all_edges, ores.display_type, orvp, ores.vi_params.vi_type_string = o.name+"res", 1, 1, 'SOLID', ores.vi_params, 'LiVi Res'
            move_to_coll(context, 'LiVi Results', ores)
            context.view_layer.layer_collection.children['LiVi Results'].exclude = 0
            context.view_layer.objects.active = ores

            while ores.material_slots:
                bpy.ops.object.material_slot_remove()

            while ores.data.shape_keys:
                context.object.active_shape_key_index = 0
                bpy.ops.object.shape_key_remove(all=True)

            ores.visible_diffuse, ores.visible_glossy, ores.visible_transmission, ores.visible_volume_scatter, ores.visible_shadow = 0, 0, 0, 0, 0
            obreslist.append(ores)
            ores.vi_params.vi_type_string == 'LiVi Res'
            orvp['omax'], orvp['omin'], orvp['oave'] = ovp['omax'], ovp['omin'], ovp['oave']
            selobj(context.view_layer, ores)

            for matname in ['{}#{}'.format('vi-suite', i) for i in range(svp.vi_leg_levels)]:
                if bpy.data.materials[matname] not in ores.data.materials[:]:
                    bpy.ops.object.material_slot_add()
                    ores.material_slots[-1].material = bpy.data.materials[matname]

            if svp.vi_disp_process == "1" and svp['liparams']['cp'] == '0':
                bm.faces.layers.int.new('extrude')
                extrude = bm.faces.layers.int['extrude']

                for face in bmesh.ops.extrude_discrete_faces(bm, faces=bm.faces)['faces']:
                    face.select = True
                    face[extrude] = 1

            bm.to_mesh(ores.data)
            bm.free()
            bpy.ops.object.shade_flat()
            ores.vi_params.lividisplay(scene)

            if svp.vi_disp_process == "1" and not ores.data.shape_keys:
                selobj(context.view_layer, ores)
                bpy.ops.object.shape_key_add(from_mix=False)

                for frame in range(svp['liparams']['fs'], svp['liparams']['fe'] + 1):
                    bpy.ops.object.shape_key_add(from_mix=False)
                    ores.active_shape_key.name, ores.active_shape_key.value = str(frame), 1

    skframe('', scene, obreslist, 'liparams')
    bpy.ops.wm.save_mainfile(check_existing=False)
    scene.frame_set(svp['liparams']['fs'])
    rendview(1)

class linumdisplay():
    def __init__(self, disp_op, context):
        scene = context.scene
        svp = scene.vi_params
        self.fn = scene.frame_current - svp['liparams']['fs']
        self.level = svp.vi_disp_3dlevel
        self.disp_op = disp_op
        svp.vi_display_rp = 0
        self.fs = svp.vi_display_rp_fs
        self.fontmult = 1
        self.obreslist = [ob for ob in scene.objects if ob.vi_params.vi_type_string == 'LiVi Res']

        if not svp.vi_display_sel_only:
            self.obd = self.obreslist
        else:
            self.obd = [context.active_object] if context.active_object in self.obreslist else []

        self.omws = [o.matrix_world for o in self.obd]
        mid_x, mid_y, self.width, self.height = viewdesc(context)
        self.view_location = retvpvloc(context)
        objmode()
        self.update(context)

    def draw(self, context):
        self.u = 0
        scene = context.scene
        svp = scene.vi_params
        bcao = bpy.context.active_object
        self.fontmult = 2  # if context.space_data.region_3d.is_perspective else 500

        if not svp.get('viparams') or svp['viparams']['vidisp'] not in ('svf', 'li', 'ss', 'lcpanel', 'rt'):
            svp.vi_display = 0
            return

        if scene.frame_current not in range(svp['liparams']['fs'], svp['liparams']['fe'] + 1):
            self.disp_op.report({'INFO'}, "Outside result frame range")
            return

        if not svp.vi_display_rp or (bcao not in self.obreslist and svp.vi_display_sel_only) or (bcao and bcao.mode == 'EDIT'):
            return

        if (self.width, self.height) != viewdesc(context)[2:]:
            mid_x, mid_y, self.width, self.height = viewdesc(context)
            self.u = 1

        if self.view_location != retvpvloc(context):
            self.view_location = retvpvloc(context)
            self.u = 1

        if not svp.vi_display_sel_only:
            obd = self.obreslist
        else:
            obd = [context.active_object] if context.active_object in self.obreslist else []

        if self.obd != obd:
            self.obd = obd
            self.u = 1

        if self.fn != scene.frame_current - svp['liparams']['fs']:
            self.fn = scene.frame_current - svp['liparams']['fs']
            self.u = 1

        if self.level != svp.vi_disp_3dlevel:
            self.level = svp.vi_disp_3dlevel
            self.u = 1

        blf_props(scene, self.width, self.height)

        if self.u:
            self.update(context)
        else:
            draw_index_distance(self.allpcs, self.allres, self.fontmult * svp.vi_display_rp_fs, svp.vi_display_rp_fc, svp.vi_display_rp_fsh, self.alldepths)

        if svp.vi_display_rp_fs != self.fs:
            self.fs = svp.vi_display_rp_fs

    def update(self, context):
        dp = context.evaluated_depsgraph_get()
        scene = context.scene
        vl = context.view_layer
        svp = scene.vi_params
        self.allpcs, self.alldepths, self.allres = array([]), array([]), array([])

        for ob in self.obd:
            if ob.data.get('shape_keys') and str(self.fn) in [sk.name for sk in ob.data.shape_keys.key_blocks] and ob.active_shape_key.name != str(self.fn):
                ob.active_shape_key_index = [sk.name for sk in ob.data.shape_keys.key_blocks].index(str(self.fn))

        for ob in self.obd:
            res = array([])
            bm = bmesh.new()
            bm.from_object(ob, dp)
            bm.transform(ob.matrix_world)
            bm.normal_update()

            if svp.li_disp_menu == 'aga1v':
                var = 'aga{}v'.format(svp.vi_views)
            elif svp.li_disp_menu == 'ago1v':
                var = 'ago{}v'.format(svp.vi_views)
            elif svp.li_disp_menu == 'rt':
                var = f'{svp.au_sources}_rt'
            else:
                var = svp.li_disp_menu

            if bm.faces.layers.float.get('{}{}'.format(var, scene.frame_current)):
                geom = bm.faces
            elif bm.verts.layers.float.get('{}{}'.format(var, scene.frame_current)):
                geom = bm.verts
            else:
                self.disp_op.report({'ERROR'}, f"No result data on {ob.name}. Re-export LiVi Context and Geometry")
                return 'CANCELLED'

            geom = bm.faces if bm.faces.layers.float.get('{}{}'.format(var, scene.frame_current)) else bm.verts
            geom.ensure_lookup_table()
            livires = geom.layers.float['{}{}'.format(var, scene.frame_current)]

            if bm.faces.layers.float.get('{}{}'.format(var, scene.frame_current)):
                if svp.vi_disp_process == "1":
                    extrude = geom.layers.int['extrude']
                    faces = [f for f in geom if f.select and f[extrude]]
                else:
                    faces = [f for f in geom if f.select]

                distances = [(self.view_location - f.calc_center_median_weighted() + svp.vi_display_rp_off * f.normal.normalized()).length for f in faces]

                if svp.vi_display_vis_only:
                    fcos = [f.calc_center_median_weighted() + svp.vi_display_rp_off * f.normal.normalized() for f in faces]
                    direcs = [self.view_location - f for f in fcos]

                    try:
                        (faces, distances) = map(list, zip(*[[f, distances[i]] for i, f in enumerate(faces) if not scene.ray_cast(vl.depsgraph, fcos[i], direcs[i], distance=distances[i])[0]]))
                    except Exception as e:
                        (faces, distances) = ([], [])

                if faces:
                    face2d = [view3d_utils.location_3d_to_region_2d(context.region, context.region_data, f.calc_center_median_weighted()) for f in faces]

                    try:
                        (faces, pcs, depths) = map(list, zip(*[[f, face2d[fi], distances[fi]] for fi, f in enumerate(faces) if
                                                               face2d[fi] and 0 < face2d[fi][0] < self.width and 0 < face2d[fi][1] < self.height]))
                    except Exception:
                        (faces, pcs, depths) = ([], [], [])

                    res = array([f[livires] for f in faces])
                    res = ret_res_vals(svp, res)

            elif bm.verts.layers.float.get('{}{}'.format(var, scene.frame_current)):
                verts = [v for v in geom if not v.hide and v.select and (context.space_data.region_3d.view_location -
                         self.view_location).dot(v.co + svp.vi_display_rp_off * v.normal.normalized() - self.view_location) /
                         ((context.space_data.region_3d.view_location-self.view_location).length * (v.co + svp.vi_display_rp_off * v.normal.normalized() - self.view_location).length) > 0]
                distances = [(self.view_location - v.co + svp.vi_display_rp_off * v.normal.normalized()).length for v in verts]

                if svp.vi_display_vis_only:
                    vcos = [v.co + svp.vi_display_rp_off * v.normal.normalized() for v in verts]
                    direcs = [self.view_location - v for v in vcos]

                    try:
                        (verts, distances) = map(list, zip(*[[v, distances[i]] for i, v in enumerate(verts) if not scene.ray_cast(vl.depsgraph, vcos[i], direcs[i], distance=distances[i])[0]]))
                    except Exception:
                        (verts, distances) = ([], [])

                if verts:
                    vert2d = [view3d_utils.location_3d_to_region_2d(context.region, context.region_data, v.co) for v in verts]
                    try:
                        (verts, pcs, depths) = map(list, zip(*[[v, vert2d[vi], distances[vi]] for vi, v in enumerate(verts) if
                                                             vert2d[vi] and 0 < vert2d[vi][0] < self.width and 0 < vert2d[vi][1] < self.height]))
                    except Exception:
                        (verts, pcs, depths) = ([], [], [])

                    res = array([v[livires] for v in verts])
                    res = ret_res_vals(svp, res)
            bm.free()

            if len(res):
                self.allpcs = nappend(self.allpcs, array(pcs))
                self.alldepths = nappend(self.alldepths, array(depths))
                self.allres = nappend(self.allres, res)

        if len(self.alldepths):
            self.alldepths = self.alldepths/nmin(self.alldepths)
            draw_index_distance(self.allpcs, self.allres, self.fontmult * svp.vi_display_rp_fs, svp.vi_display_rp_fc, svp.vi_display_rp_fsh, self.alldepths)


class Base_Display():
    def __init__(self, ipos, width, height, xdiff, ydiff):
        self.ispos = ipos
        self.iepos = [ipos[0] + 40, ipos[1] + 40]
        self.xdiff, self.ydiff = xdiff, ydiff
        self.lspos = [ipos[0] - 5, ipos[1] - self.ydiff - 25]
        self.lepos = [ipos[0] - 5 + self.xdiff, self.lspos[1] + self.ydiff]
        self.resize, self.move, self.expand = 0, 0, 0
        self.hl = [1, 1, 1, 1]
        self.cao = None


class results_bar():
    def __init__(self, images):
        self.images = images
        self.rh = 0
        self.xpos = 0
        self.shaders = [gpu.shader.from_builtin('UNIFORM_COLOR'), gpu.shader.from_builtin('UNIFORM_COLOR')]
        self.f_indices = ((0, 1, 2), (2, 3, 0))
        self.tex_coords = ((0, 0), (1, 0), (1, 1), (0, 1))
        self.no = len(images)
        self.yoffset = 10
        self.size = 50
        self.isize = self.size - 10
        self.iyoffset = self.yoffset + (self.size - self.isize)/2
        self.ixoffset = self.isize + 5
        self.iyoffsetb = self.iyoffset + self.isize
        self.ipos = []

        for ii, im in enumerate(images):
            if im not in bpy.data.images:
                bpy.data.images.load(os.path.join(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))), 'Images', im))
            
            self.shaders.append(gpu.shader.from_builtin('IMAGE'))
            pos = self.ret_coords(self.xpos, self.rh, ii)
            self.ipos.append(pos)

    def ret_coords(self, xpos, rh, no):
        return ((xpos + 5 + no * self.size, rh - self.iyoffsetb),
                (xpos + self.ixoffset + no * self.size, rh - self.iyoffsetb),
                (xpos + self.ixoffset + no * self.size, rh - self.iyoffset),
                (xpos + 5 + no * self.size, rh - self.iyoffset))

    def draw(self, xpos, rh):
        if self.rh != rh or xpos != self.xpos - 10:
            self.ipos = []
            self.rh = rh
            self.xpos = xpos

            v_coords = ((self.xpos, rh - self.yoffset - self.size), (self.xpos + self.no * self.size, rh - self.yoffset - self.size),
                        (self.xpos + self.no * self.size, rh - self.yoffset), (self.xpos, rh - self.yoffset), (self.xpos, rh - self.yoffset - self.size))

            self.batches = [batch_for_shader(self.shaders[1], 'TRIS', {"pos": v_coords}, indices=self.f_indices),
                            batch_for_shader(self.shaders[0], 'LINE_STRIP', {"pos": v_coords})]

            for i in range(self.no):
                pos = self.ret_coords(self.xpos, rh, i)
                self.batches.append(batch_for_shader(self.shaders[i + 2], 'TRI_FAN', {"pos": pos, "texCoord": self.tex_coords}))
                self.ipos.append(pos)
        
        for si, s in enumerate(self.shaders):
            s.bind()

            if si == 0:
                s.uniform_float("color", (1, 1, 1, 1))
                self.batches[si].draw(s)
            elif si == 1:
                s.uniform_float("color", (0, 0, 0, 1))
                self.batches[si].draw(s)
            else:
                im = bpy.data.images[self.images[si - 2]]
                texture = gpu.texture.from_image(im)

                if im.gl_load():
                    raise Exception()

                s.uniform_sampler("image", texture)
                #s.uniform_int("image", 0)
                self.batches[si].draw(s)


class draw_bsdf(Base_Display):
    def __init__(self, context, unit, pos, width, height, xdiff, ydiff):
        Base_Display.__init__(self, pos, width, height, xdiff, ydiff)
        self.plt = plt
        self.pw, self.ph = 0.175 * xdiff, 0.35 * ydiff
        self.type_select = 0
        self.patch_hl = 0
        self.scale_select = 'Log'
        self.buttons = {}
        self.num_disp = 0
        self.leg_max, self.leg_min = 100, 0
        self.base_unit = unit
        self.font_id = blf.load(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'Fonts', 'NotoSans-Regular.ttf'))
        self.dpi = 157
        self.v_coords = [(0, 0), (0, 1), (1, 1), (1, 0)]
        self.f_indices = [(0, 1, 2), (2, 3, 0)]
        self.segments = (1, 8, 16, 20, 24, 24, 24, 16, 12)
        self.radii = (6.9, 20.6, 34.4, 48.2, 62, 75.8, 89.5, 103.3, 124)
        self.f_colours = [(1, 1, 1, 1)] * (721 + 8 * 720)
        self.imspos = (self.lspos[0], self.lspos[1])
        self.image = 'bsdfplot.png'
        self.isize = (self.xdiff, self.xdiff - 50)
        self.iimage = 'bsdf_empty.png'
        self.iisize = (250, 270)
        self.type = context.active_object.active_material.vi_params['bsdf']['type']

        if self.iimage not in [im.name for im in bpy.data.images]:
            bpy.data.images.load(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'Images', 'bsdf_empty.png'))

        self.vi_coords = [(0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0)]
        self.tex_coords = ((0, 0), (0, 1), (1, 1), (1, 0))
        self.sr = 0
        self.cr = 0
        self.sseg = 1
        self.cseg = 0
        self.srs = 1
        self.crs = 0
        self.create_batch('all')
        self.get_data(context)

    def get_data(self, context):
        self.mat = context.active_object.active_material
        bsdf = parseString(self.mat.vi_params['bsdf']['xml'])
        self.radtype = [path.firstChild.data for path in bsdf.getElementsByTagName('Wavelength')]
        self.rad_select = self.radtype[0]
        self.dattype = [path.firstChild.data for path in bsdf.getElementsByTagName('WavelengthDataDirection')]
        self.direc = self.dattype[0]
        self.type_select = self.dattype[0].split()[0]
        self.dir_select = self.dattype[0].split()[1]
        self.uthetas = [float(path.firstChild.data) for path in bsdf.getElementsByTagName('UpperTheta')]
        self.phis = [int(path.firstChild.data) for path in bsdf.getElementsByTagName('nPhis')]

        if ',' in bsdf.getElementsByTagName('ScatteringData')[0].firstChild.data:
            self.scatdat = [array([float(nv) for nv in path.firstChild.data.strip('\t').strip('\n').strip().split(',') if nv]) for path in bsdf.getElementsByTagName('ScatteringData')]
        else:
            self.scatdat = [array([float(nv) for nv in path.firstChild.data.strip('\t').strip('\n').strip(',').split(' ') if nv]) for path in bsdf.getElementsByTagName('ScatteringData')]

        self.plot(context)

    def plot(self, context):
        scene = context.scene
        svp = scene.vi_params
        leg_min = svp.vi_bsdfleg_min if svp.vi_bsdfleg_scale == '0' or svp.vi_bsdfleg_min > 0 else svp.vi_bsdfleg_min + 0.01
        self.col = svp.vi_leg_col
        self.centre = (self.lspos[0] + 0.225 * self.xdiff, self.lspos[1] + 0.425 * self.ydiff)
        self.plt.clf()
        self.plt.close()
        self.fig = self.plt.figure(figsize=(4, 3.5), dpi=100)
        ax = self.plt.subplot(111, projection='polar')
        ax.bar(0, 0)
        self.plt.title('{} {}'.format(self.rad_select, svp.vi_bsdf_direc), size=9, y=1.025)
        ax.axis([0, 2 * pi, 0, 1])
        ax.spines['polar'].set_visible(False)
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])

        for dti, dt in enumerate(self.dattype):
            if dt == svp.vi_bsdf_direc:
                self.scat_select = dti
                break

        selectdat = self.scatdat[self.scat_select].reshape(145, 145)  # if self.scale_select == 'Linear' else nlog10((self.scatdat[self.scat_select] + 1).reshape(145, 145))
        widths = [0] + [self.uthetas[w]/90 for w in range(9)]
        patches, p = [], 0
        sa = repeat(kfsa, self.phis)
        act = repeat(kfact, self.phis)
        patchdat = selectdat[self.sseg - 1] * act * sa * 100
        bg = self.plt.Rectangle((0, 0), 2 * pi, 1, color=mcm.get_cmap(svp.vi_leg_col)((0, 0.01)[svp.vi_bsdfleg_scale == '1']), zorder=0)

        for ring in range(1, 10):
            angdiv = pi/self.phis[ring - 1]
            anglerange = range(self.phis[ring - 1], 0, -1)  # if self.type_select == 'Transmission' else range(self.phis[ring - 1])
            ri = widths[ring] - widths[ring-1]

            for wedge in anglerange:
                phi1, phi2 = wedge * 2 * angdiv - angdiv, (wedge + 1) * 2 * angdiv - angdiv
                patches.append(Rectangle((phi1, widths[ring - 1]), phi2 - phi1, ri))

                if self.num_disp:
                    y = 0 if ring == 1 else 0.5 * (widths[ring] + widths[ring-1])
                    self.plt.text(0.5 * (phi1 + phi2), y, ('{:.1f}', '{:.0f}')[patchdat[p] >= 10].format(patchdat[p]), ha="center",
                                  va='center', family='sans-serif', size=self.num_disp)
                p += 1

        pc = PatchCollection(patches, norm=mcolors.LogNorm(vmin=leg_min, vmax=svp.vi_bsdfleg_max), cmap=self.col, linewidths=[0] + 144*[0.5],
                             edgecolors=('black',)) if svp.vi_bsdfleg_scale == '1' else PatchCollection(patches, cmap=self.col, linewidths=[0] + 144*[0.5], edgecolors=('black',))
        pc.set_array(patchdat)
        ax.add_collection(pc)
        ax.add_artist(bg)
        cb = self.plt.colorbar(pc, fraction=0.04, pad=0.02, format='%3g')
        cb.set_label(label='Percentage of incoming flux (%)', size=9)
        cb.ax.tick_params(labelsize=7)
        pc.set_clim(vmin=leg_min + 0.01, vmax=svp.vi_bsdfleg_max)
        self.plt.tight_layout()
        self.save(scene)

    def save(self, scene):
        mp2im(self.fig, self.image)

    def ret_coords(self):
        if self.type == 'LBNL/Klems Full':
            vl_coords, f_indices = [], []
            v_coords = [(0, 0)]

            if self.sseg == 1:
                va_coords = [(0, 0)]
                va_coords += [(self.radii[self.sr] * cos((x*0.5 - 360/(2 * self.segments[self.sr]))*pi/180),
                              self.radii[self.sr] * sin((x*0.5 - 360/(2 * self.segments[self.sr]))*pi/180)) for x in range(720)]
                fa_indices = [(0, x + (1, -719)[x > 0 and not x % 720], x) for x in range(1, 720 + 1)]
            else:
                va_coords = [(-self.radii[self.sr - 1] * cos((x*0.5 - 360/(2 * self.segments[self.sr]))*pi/180),
                              self.radii[self.sr - 1] * sin((x*0.5 - 360/(2 * self.segments[self.sr]))*pi/180)) for x in range(int((self.srs - 1) * 720/self.segments[self.sr]),
                             int((self.srs) * 720/self.segments[self.sr]) + 1)]
                va_coords += [(-self.radii[self.sr] * cos((x*0.5 - 360/(2 * self.segments[self.sr]))*pi/180),
                              self.radii[self.sr] * sin((x*0.5 - 360/(2 * self.segments[self.sr]))*pi/180)) for x in range(int((self.srs - 1) * 720/self.segments[self.sr]),
                              int((self.srs)*720/self.segments[self.sr]) + 1)]

                fa_indices = [(x, x+1, int(x+720/self.segments[self.sr] + 1)) for x in range(int(720/self.segments[self.sr]))] + \
                             [(x-1, x, int(x-720/self.segments[self.sr] - 1)) for x in range(int(720/self.segments[self.sr] + 2), int(1440/self.segments[self.sr] + 2))]

            fa_colours = [(1, 1, 0.0, 1) for _ in range(len(va_coords))]

            for ri, radius in enumerate(self.radii):
                if ri < 8:
                    for si, s in enumerate(range(self.segments[ri + 1])):
                        vl_coords += [(radius * cos((360/self.segments[ri + 1] * si + 360/(2 * self.segments[ri + 1])) * pi/180),
                                      radius * sin((360/self.segments[ri + 1] * si + 360/(2 * self.segments[ri + 1])) * pi/180)),
                                      (self.radii[ri + 1] * cos((360/self.segments[ri + 1] * si + 360/(2 * self.segments[ri + 1])) * pi/180),
                                       self.radii[ri + 1] * sin((360/self.segments[ri + 1] * si + 360/(2 * self.segments[ri + 1])) * pi/180))]

                vl_coords1 = [[radius * cos((x*0.5 - 360/(2 * self.segments[ri]))*pi/180), radius * sin((x*0.5 - 360/(2 * self.segments[ri]))*pi/180)] for x in range(720)]
                vl_coords2 = [[radius * cos(((x*0.5 + 0.5) - 360/(2 * self.segments[ri]))*pi/180), radius * sin(((x*0.5 + 0.5) - 360/(2 * self.segments[ri]))*pi/180)] for x in range(720)]
                vl_coordstot = list(zip(vl_coords1, vl_coords2))
                vl_coords += [item for sublist in vl_coordstot for item in sublist]
                v_coords += vl_coords1

            f_indices = [(0, x + (1, -719)[x > 0 and not x % 720], x) for x in range(1, 720*9 + 1)][::-1]
            return (vl_coords, v_coords, f_indices, va_coords, fa_colours, fa_indices)

    def create_batch(self, sel):
        line_shader = gpu.types.GPUShaderCreateInfo()
        line_shader.push_constant('MAT4', 'ModelViewProjectionMatrix')
        line_shader.push_constant('VEC2', 'spos')
        line_shader.push_constant('FLOAT', 'zpos')
        line_shader.push_constant('VEC2', 'size')
        line_shader.push_constant('VEC4', 'colour')
        line_shader.vertex_in(0, 'VEC2', 'position')
        line_shader.fragment_out(0, 'VEC4', 'FragColour')
        line_shader.vertex_source(
            '''
            void main()
                {
                    float xpos = spos[0] + position[0] * size[0];
                    float ypos = spos[1] + position[1] * size[1];
                    gl_Position = ModelViewProjectionMatrix * vec4(int(xpos), int(ypos), zpos, 1.0f);
                    vec4 pp = gl_Position;
                }
            '''
        )
        line_shader.fragment_source(
            '''
            void main()
                {
                    FragColour = colour;
                }
            '''
        )

        arc_shader = gpu.types.GPUShaderCreateInfo()
        arc_shader.push_constant('MAT4', 'ModelViewProjectionMatrix')
        arc_shader.push_constant('VEC2', 'spos')
        arc_shader.push_constant('VEC2', 'size')
        arc_shader.push_constant('FLOAT', 'zpos')
        arc_shader.vertex_in(0, 'VEC4', 'colour')
        arc_shader.vertex_in(1, 'VEC2', 'position')
        arc_shader_iface = gpu.types.GPUStageInterfaceInfo('arc_interface')
        arc_shader_iface.flat('VEC4', 'a_colour')
        arc_shader.vertex_out(arc_shader_iface)
        arc_shader.fragment_out(0, 'VEC4', 'FragColour')
        arc_shader.vertex_source(
            '''
                void main()
                    {
                    float xpos = spos[0] + position[0] * size[0];
                    float ypos = spos[1] + position[1] * size[1];
                    gl_Position = ModelViewProjectionMatrix * vec4(int(xpos), int(ypos), zpos, 1.0f);
                    a_colour = colour;
                    }
            '''
        )
        arc_shader.fragment_source(
            '''
                void main()
                    {
                        FragColour = a_colour;
                    }
            '''
        )

        (vl_coords, v_coords, f_indices, va_coords, fa_colours, fa_indices) = self.ret_coords()

        if sel == 'all':
            b_coords = [(0, 0), (0, 1), (1, 1), (1, 0)]
            b_indices = [(0, 1, 3), (1, 2, 3)]
            b_colours = [(1, 1, 1, 1), (1, 1, 1, 1), (1, 1, 1, 1), (1, 1, 1, 1)]
            self.back_shader = gpu.shader.create_from_info(arc_shader)
            self.arcline_shader = gpu.shader.create_from_info(line_shader)
            self.image_shader = gpu.shader.from_builtin('IMAGE')
            self.vi2_coords = [(self.lspos[0], self.lspos[1]), (self.lspos[0], self.lspos[1] + self.isize[1]), (self.lspos[0] + self.isize[0], self.lspos[1] + self.isize[1]), (self.lspos[0] + self.isize[0], self.lspos[1])]
            self.image_batch = batch_for_shader(self.image_shader, 'TRI_FAN', {"pos": self.vi2_coords, "texCoord": self.tex_coords})
            self.vi3_coords = [(self.lspos[0] + 50, self.lepos[1] - 280), (self.lspos[0] + 50, self.lepos[1] - 280 + self.iisize[1]), (self.lspos[0] + 50 + self.iisize[0], self.lepos[1] - 280 + self.iisize[1]), (self.lspos[0] + 50 + self.iisize[0], self.lepos[1] - 280)]
            self.iimage_shader = gpu.shader.from_builtin('IMAGE')
            self.iimage_batch = batch_for_shader(self.iimage_shader, 'TRI_FAN', {"pos": self.vi3_coords, "texCoord": self.tex_coords})
            self.back_batch = batch_for_shader(self.back_shader, 'TRIS', {"position": b_coords, "colour": b_colours}, indices=b_indices)

        self.arc_shader = gpu.shader.create_from_info(arc_shader)
        self.arc_batch = batch_for_shader(self.arc_shader, 'TRIS', {"position": va_coords, "colour": fa_colours}, indices=fa_indices)

    def draw(self, context):
        if self.expand:
            (r0h, r2w, r5h) = ret_dcoords(context)
            self.back_shader.bind()
            self.back_shader.uniform_float("size", (400, 650))
            self.back_shader.uniform_float("spos", (self.lspos))
            self.back_shader.uniform_float("zpos", -0.5)
            self.back_batch.draw(self.back_shader)
            gpu.state.depth_test_set('LESS')
            gpu.state.depth_mask_set(False)
            gpu.state.line_width_set(1)
            self.iimage_shader.bind()
            iim = bpy.data.images[self.iimage]
            texture = gpu.texture.from_image(iim)

            if iim.gl_load():
                raise Exception()

            self.iimage_shader.uniform_sampler("image", texture)
            self.vi2_coords = [(self.lspos[0], self.lspos[1]), (self.lspos[0], self.lspos[1] + self.isize[1]), (self.lspos[0] + self.isize[0], self.lspos[1] + self.isize[1]), (self.lspos[0] + self.isize[0], self.lspos[1])]
            self.vi3_coords = [(self.lspos[0] + 50, self.lepos[1] - 280), (self.lspos[0] + 50, self.lepos[1] - 280 + self.iisize[1]), (self.lspos[0] + 50 + self.iisize[0], self.lepos[1] - 280 + self.iisize[1]), (self.lspos[0] + 50 + self.iisize[0], self.lepos[1] - 280)]
            self.iimage_batch = batch_for_shader(self.iimage_shader, 'TRI_FAN', {"pos": self.vi3_coords, "texCoord": self.tex_coords})
            self.iimage_batch.draw(self.iimage_shader)
            self.arc_shader.bind()
            self.arc_shader.uniform_float("size", (1, 1))
            self.arc_shader.uniform_float("spos", (self.lspos[0] + 50 + 125, self.lepos[1] - 155))
            self.arc_shader.uniform_float("zpos", -0.25)
            self.arc_batch.draw(self.arc_shader)
            self.image_shader.bind()
            im = bpy.data.images[self.image]
            texture = gpu.texture.from_image(im)

            if im.gl_load():
                raise Exception()

            self.image_shader.uniform_sampler("image", texture)
            self.image_batch = batch_for_shader(self.image_shader, 'TRI_FAN', {"pos": self.vi2_coords, "texCoord": self.tex_coords})
            self.image_batch.draw(self.image_shader)


class wr_legend(Base_Display):
    def __init__(self, context, unit, pos, width, height, xdiff, ydiff):
        Base_Display.__init__(self, pos, width, height, xdiff, ydiff)
        self.unit = unit
        self.font_id = 0
        self.dpi = 72
        self.update(context)
        self.create_batch()
        self.line_shader.bind()
        self.line_shader.uniform_float("colour", (0, 0, 0, 1))

    def update(self, context):
        scene = context.scene
        svp = scene.vi_params
        simnode = bpy.data.node_groups[svp['viparams']['restree']].nodes[svp['viparams']['resnode']]
        self.cao = context.active_object

        if self.cao and self.cao.vi_params.get('VIType') and self.cao.vi_params['VIType'] == 'Wind_Plane':
            self.levels = self.cao.vi_params['nbins']
            maxres = self.cao.vi_params['maxres']
        else:
            self.levels = simnode['nbins']
            maxres = simnode['maxres']

        self.cols = retcols(mcm.get_cmap(svp.vi_leg_col), self.levels)
        self.colours = [item for item in [self.cols[i] for i in range(self.levels)] for i in range(4)]

        if not svp.get('liparams'):
            svp.vi_display = 0
            return

        self.resvals = ['{0:.0f} - {1:.0f}'.format(2*i, 2*(i+1)) for i in range(simnode['nbins'])]
        self.resvals[-1] = self.resvals[-1][:-int(len('{:.0f}'.format(maxres)))] + "Inf"
        blf.size(self.font_id, 12 * self.dpi/72)
        self.titxdimen = blf.dimensions(self.font_id, self.unit)[0]
        self.resxdimen = blf.dimensions(self.font_id, self.resvals[-1])[0]
        self.mydimen = blf.dimensions(self.font_id, self.unit)[1]

    def draw(self, context):
        (r0h, r2w, r5h) = ret_dcoords(context)

        if self.expand:
            if self.resize:
                self.xdiff = self.lepos[0] - self.lspos[0]
                self.ydiff = self.lepos[1] - self.lspos[1]
            elif self.move:
                self.lspos[1] = self.lepos[1] - self.ydiff
                self.lepos[0] = self.lspos[0] + self.xdiff
            if self.lepos[1] > r5h - r0h:
                self.lspos[1] = r5h - r0h - self.ydiff
                self.lepos[1] = r5h - r0h
            if self.lspos[0] < r2w:
                self.lspos[0] = r2w

            self.base_shader.bind()
            self.base_shader.uniform_float("size", (self.xdiff, self.ydiff))
            self.base_shader.uniform_float("spos", self.lspos)
            self.base_shader.uniform_float("colour", self.hl)
            self.base_batch.draw(self.base_shader)
            self.col_shader.bind()
            self.col_shader.uniform_float("size", (self.xdiff, self.ydiff))
            self.col_shader.uniform_float("spos", self.lspos)
            self.col_batch.draw(self.col_shader)
            self.line_shader.bind()
            self.line_shader.uniform_float("size", (self.xdiff, self.ydiff))
            self.line_shader.uniform_float("spos", self.lspos)
            self.line_shader.uniform_float("colour", (0, 0, 0, 1))
            self.line_batch.draw(self.line_shader)
            fontscale = max(self.titxdimen/(self.xdiff * 0.9), self.resxdimen/(self.xdiff * 0.65), self.mydimen * 1.25/(self.lh * self.ydiff))
            blf.enable(self.font_id, 4)
            blf.enable(self.font_id, 8)
            # blf.shadow(self.font_id, 3, 0.7, 0.7, 0.7, 1)
            blf.size(self.font_id, 12 * int(self.dpi/fontscale*72))
            blf.position(self.font_id, self.lspos[0] + (self.xdiff - blf.dimensions(self.font_id, self.unit)[0]) * 0.45,
                         self.lepos[1] - 0.5 * (self.lh * self.ydiff) - blf.dimensions(self.font_id, self.unit)[1] * 0.3, 0)
            blf.color(self.font_id, 0, 0, 0, 1)
            blf.draw(self.font_id, self.unit)
            # blf.shadow(self.font_id, 3, 0.8, 0.8, 0.8, 1)
            blf.size(self.font_id, 14* int(self.dpi/fontscale*72))

            for i in range(self.levels):
                num = self.resvals[i]
                ndimen = blf.dimensions(self.font_id, "{}".format(num))
                blf.position(self.font_id, int(self.lepos[0] - self.xdiff * 0.05 - ndimen[0]),
                             int(self.lspos[1] + i * self.lh * self.ydiff) + int((self.lh * self.ydiff - ndimen[1])*0.55), 0)
                blf.draw(self.font_id, "{}".format(self.resvals[i]))

            blf.disable(self.font_id, 8)
            blf.disable(self.font_id, 4)

    def create_batch(self):
        base_shader = gpu.types.GPUShaderCreateInfo()
        base_shader.push_constant('MAT4', 'ModelViewProjectionMatrix')
        base_shader.push_constant('VEC2', 'spos')
        base_shader.push_constant('VEC2', 'size')
        base_shader.push_constant('VEC4', 'colour')
        base_shader.vertex_in(0, 'VEC2', 'position')
        base_shader.fragment_out(0, 'VEC4', 'FragColour')
        base_shader.vertex_source(
            '''
            void main()
            {
                float xpos = spos[0] + position[0] * size[0];
                float ypos = spos[1] + position[1] * size[1];
                gl_Position = ModelViewProjectionMatrix * vec4(int(xpos), int(ypos), -0.1f, 1.0f);
            }
            '''
        )
        base_shader.fragment_source(
            '''
            void main()
                {
                    FragColour = colour;
                }
            '''
        )

        col_shader = gpu.types.GPUShaderCreateInfo()
        col_shader.push_constant('MAT4', 'ModelViewProjectionMatrix')
        col_shader.push_constant('VEC2', 'spos')
        col_shader.push_constant('VEC2', 'size')
        # col_shader.push_constant('VEC4', 'colour')
        col_shader.vertex_in(0, 'VEC2', 'position')
        col_shader.vertex_in(1, 'VEC4', 'colour')
        col_shader_iface = gpu.types.GPUStageInterfaceInfo('col_interface')
        col_shader_iface.flat('VEC4', 'f_colour')
        col_shader.vertex_out(col_shader_iface)
        col_shader.fragment_out(0, 'VEC4', 'FragColour')
        col_shader.vertex_source(
            '''
            void main()
            {
                float xpos = spos[0] + position[0] * size[0];
                float ypos = spos[1] + position[1] * size[1];
                gl_Position = ModelViewProjectionMatrix * vec4(int(xpos), int(ypos), 0.0f, 1.0f);
                f_colour = colour;
            }
            '''
        )
        col_shader.fragment_source(
            '''
            void main()
                {
                    FragColour = f_colour;
                }
            '''
        )

        self.base_shader = gpu.shader.create_from_info(base_shader)
        self.line_shader = gpu.shader.create_from_info(base_shader)
        self.col_shader = gpu.shader.create_from_info(col_shader)
        v_coords = [(0, 0), (0, 1), (1, 1), (1, 0)]
        f_indices = [(0, 1, 2), (2, 3, 0)]
        lh = 1/(self.levels + 1.25)
        vl_coords = v_coords
        f_indices = [(0, 1, 2), (2, 3, 0)]
        fl1_indices = [tuple(array((0, 1, 2)) + 4 * i) for i in range(self.levels)]
        fl2_indices = [tuple(array((2, 3, 0)) + 4 * i) for i in range(self.levels)]
        fl_indices = list(fl1_indices) + list(fl2_indices)

        for i in range(0, self.levels):
            vl_coords += [(0, i * lh), (0.4, i * lh), (0.4, (i + 1) * lh), (0, (i + 1) * lh)]

        self.base_batch = batch_for_shader(self.base_shader, 'TRIS', {"position": v_coords}, indices=f_indices)
        self.line_batch = batch_for_shader(self.line_shader, 'LINE_STRIP', {"position": vl_coords})
        self.col_batch = batch_for_shader(self.col_shader, 'TRIS', {"position": vl_coords[4:], "colour": self.colours}, indices=fl_indices)
        self.lh = lh


class wr_table(Base_Display):
    def __init__(self, context, pos, width, height, xdiff, ydiff):
        Base_Display.__init__(self, pos, width, height, xdiff, ydiff)
        self.font_id = 0
        self.dpi = int(0.15 * ydiff)
        self.update(context)
        self.create_batch()
        self.line_shader.bind()
        self.line_shader.uniform_float("colour", (0, 0, 0, 1))

    def update(self, context):
        self.cao = context.active_object

        if self.cao and self.cao.vi_params.get('VIType') == 'Wind_Plane':
            self.rcarray = array(self.cao.vi_params['table'])
        else:
            self.rcarray = array([['Invalid object']])

    def create_batch(self):
        base_shader = gpu.types.GPUShaderCreateInfo()
        base_shader.push_constant('MAT4', 'ModelViewProjectionMatrix')
        base_shader.push_constant('VEC2', 'spos')
        base_shader.push_constant('VEC2', 'size')
        base_shader.push_constant('VEC4', 'colour')
        base_shader.vertex_in(0, 'VEC2', 'position')
        base_shader.fragment_out(0, 'VEC4', 'FragColour')
        base_shader.vertex_source(
            '''
            void main()
            {
                float xpos = spos[0] + position[0] * size[0];
                float ypos = spos[1] + position[1] * size[1];
                gl_Position = ModelViewProjectionMatrix * vec4(int(xpos), int(ypos), 0.0f, 1.0f);
            }
            '''
        )
        base_shader.fragment_source(
            '''
            void main()
                {
                    FragColour = colour;
                }
            '''
        )

        self.base_shader = gpu.shader.create_from_info(base_shader)
        self.line2_shader = gpu.shader.from_builtin('POLYLINE_SMOOTH_COLOR')
        self.line2_shader.uniform_float("viewportSize", gpu.state.viewport_get()[2:])
        self.line2_shader.uniform_float("lineWidth", 5)
        self.line_shader = gpu.shader.create_from_info(base_shader)
        v_coords = [(0, 0), (0, 1), (1, 1), (1, 0)]
        f_indices = [(0, 1, 2), (2, 3, 0)]
        vl_coords = v_coords
        rno = len(self.rcarray)
        cno = len(self.rcarray[0])
        rh = 1/rno
        blf.size(0, 24*300/72)
        ctws = array([int(max([blf.dimensions(0, 'u{}'.format(e))[0] for e in entry])) for entry in self.rcarray.T])
        ctws = ctws/sum(ctws)
        ctws = [sum(ctws[:i]) for i in range(4)] + [1]

        for ci in range(cno):
            for ri in range(rno):
                vl_coords += [(ctws[ci], ri * rh), (ctws[ci + 1], ri * rh), (ctws[ci + 1], (ri + 1) * rh)]  # , (ci * rw, (ri + 1) * rh), (ci * rw, ri * rh)]

        vl_coords += [(0, 1)]
        # vln_coords= [(100, 100), (200, 100), (200, 150), (100, 175)]
        vertex_colors = [(0, 0, 0, 1) for x in range(4)]
        vertex_colors[-2] = (0, 0, 1, 0)
        vertex_colors[-1] = (0, 0, 1, 1)
        self.base_batch = batch_for_shader(self.base_shader, 'TRIS', {"position": v_coords}, indices=f_indices)
        # self.line2_batch = batch_for_shader(self.line2_shader, 'LINE_STRIP', {"pos": vln_coords, "color": vertex_colors})
        self.line_batch = batch_for_shader(self.line_shader, 'LINE_STRIP', {"position": vl_coords})
        # self.line_batch.draw(self.line_shader)

    def draw(self, context):
        (r0h, r2w, r5h) = ret_dcoords(context)

        if self.expand:
            if self.resize:
                self.xdiff = self.lepos[0] - self.lspos[0]
                self.ydiff = self.lepos[1] - self.lspos[1]
            elif self.move:
                self.lspos[1] = self.lepos[1] - self.ydiff
                self.lepos[0] = self.lspos[0] + self.xdiff
            if self.lepos[1] > r5h - r0h:
                self.lspos[1] = r5h - r0h - self.ydiff
                self.lepos[1] = r5h - r0h
            if self.lspos[0] < r2w:
                self.lspos[0] = r2w

            self.base_shader.bind()
            self.base_shader.uniform_float("size", (self.xdiff, self.ydiff))
            self.base_shader.uniform_float("spos", self.lspos)
            self.base_shader.uniform_float("colour", self.hl)
            self.base_batch.draw(self.base_shader)
            self.line_shader.bind()
            self.line_shader.uniform_float("size", (self.xdiff, self.ydiff))
            self.line_shader.uniform_float("spos", self.lspos)
            self.line_shader.uniform_float("colour", (0, 0, 0, 1))
            self.line_batch.draw(self.line_shader)
            # self.line2_batch.draw(self.line2_shader)
            # self.line2_shader.bind()

            fid = self.font_id
            blf.enable(fid, 4)
            blf.enable(fid, 8)
            # blf.shadow(self.font_id, 5, 0.7, 0.7, 0.7, 1)
            blf.size(fid, 24*300/72)
            rcshape = self.rcarray.shape
            [rowno, colno] = self.rcarray.shape
            ctws = array([int(max([blf.dimensions(fid, '{}'.format(e))[0] for e in entry])) for entry in self.rcarray.T])
            ctws = self.xdiff * ctws/sum(ctws)
            ctws = [sum(ctws[:i]) for i in range(4)] + [self.xdiff]
            ctws = [sum(ctws[i:i+2])/2 for i in range(4)]
            coltextwidths = array([int(max([blf.dimensions(fid, '{}'.format(e))[0] for e in entry]) + 0.05 * self.xdiff) for entry in self.rcarray.T])
            colscale = sum(coltextwidths)/(self.xdiff * 0.98)
            maxrowtextheight = max([max([blf.dimensions(fid, '{}'.format(e))[1] for e in entry if e]) for entry in self.rcarray.T])
            rowtextheight = maxrowtextheight + 0.1 * self.ydiff/rowno
            rowscale = (rowno * rowtextheight)/(self.ydiff - self.xdiff * 0.025)
            rowheight = int((self.ydiff - self.xdiff * 0.01)/rowno)
            rowtops = [int(self.lepos[1] - self.xdiff * 0.005 - r * rowheight) for r in range(rowno)]
            rowbots = [int(self.lepos[1] - self.xdiff * 0.005 - (r + 1) * rowheight) for r in range(rowno)]
            rowmids = [0.5 * (rowtops[r] + rowbots[r]) for r in range(rowno)]

            if abs(max(colscale, rowscale) - 1) > 0.05:
                self.fontdpi = int(280/max(colscale, rowscale))

            blf.size(fid, 24 * self.fontdpi/72)
            blf.color(fid, 0, 0, 0, 1)

            for r in range(rcshape[0]):
                for c in range(rcshape[1]):
                    if self.rcarray[r][c]:
                        if c == 0:
                            blf.position(fid, self.lspos[0] + 0.01 * self.xdiff, int(rowmids[r] - 0.4 * blf.dimensions(fid, 'H')[1]), 0)
                        else:
                            blf.position(fid, self.lspos[0] + ctws[c] - int(blf.dimensions(fid, '{}'.format(self.rcarray[r][c]))[0] * 0.5), int(rowmids[r] - 0.5 * blf.dimensions(fid, 'H')[1]), 0)
                        blf.draw(fid, '{}'.format(self.rcarray[r][c]))

            blf.disable(fid, 8)
            blf.disable(fid, 4)


class draw_legend(Base_Display):
    def __init__(self, context, unit, pos, width, height, xdiff, ydiff, levels):
        Base_Display.__init__(self, pos, width, height, xdiff, ydiff)
        self.base_unit = unit
        self.font_id = blf.load(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'Fonts/NotoSans-Light.ttf'))
        self.dpi = 72
        self.levels = levels
        self.v_coords = [(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)]
        self.f_indices = [(0, 1, 2), (2, 3, 0)]
        self.update(context)
        self.create_batch()

    def update(self, context):
        scene = context.scene
        svp = scene.vi_params

        if svp.li_disp_menu != 'None':
            self.levels = svp.vi_leg_levels
            self.lh = 1/(self.levels + 1.5)
            self.cao = context.active_object
            self.cols = retcols(mcm.get_cmap(svp.vi_leg_col), self.levels)
            (self.minres, self.maxres) = leg_min_max(svp)
            self.col, self.scale = svp.vi_leg_col, svp.vi_leg_scale
            self.unit = res2unit[svp.li_disp_menu] if not svp.vi_leg_unit else svp.vi_leg_unit
            self.cols = retcols(mcm.get_cmap(svp.vi_leg_col), self.levels)
            resdiff = self.maxres - self.minres

            if not svp.get('liparams'):
                svp.vi_display = 0
                return

            dplaces = retdp(self.maxres, 1)
            resvals = [format(self.minres + i*(resdiff)/self.levels, '.{}f'.format(dplaces)) for i in range(self.levels + 1)] if self.scale == '0' else \
                      [format(self.minres + (1 - log10(i)/log10(self.levels + 1))*(resdiff), '.{}f'.format(dplaces)) for i in range(1, self.levels + 2)[::-1]]

            if svp.vi_res_process == '2' and 'restext' in bpy.app.driver_namespace.keys():
                self.resvals = [''] + bpy.app.driver_namespace.get('restext')()
                if len(self.resvals) != self.levels:
                    self.resvals = ['{0}'.format(resvals[i]) for i in range(self.levels + 1)]
            else:
                self.resvals = ['{0}'.format(resvals[i]) for i in range(self.levels + 1)]

            self.colours = [item for item in [self.cols[i] for i in range(self.levels)] for i in range(5)][:-1]
            blf.size(self.font_id, 12 * self.dpi/72)
            self.titxdimen = blf.dimensions(self.font_id, self.unit)[0]
            self.resxdimen = blf.dimensions(self.font_id, self.resvals[-1])[0]
            self.mydimen = blf.dimensions(self.font_id, 'M')[1]

    def ret_coords(self):
        lh = 1/(self.levels + 1.5)
        vl_coords = []
        fl1_indices = [tuple(array((0, 1, 2)) + 5 * i) for i in range(self.levels)]
        fl2_indices = [tuple(array((2, 3, 4)) + 5 * i) for i in range(self.levels - 1)]
        fl3_indices = [tuple(array((0, 2, (3, 4)[i < self.levels - 1])) + 5 * i) for i in range(self.levels)]
        fl_indices = list(fl1_indices) + list(fl2_indices) + list(fl3_indices)

        for i in range(0, self.levels):
            if i < self.levels - 1:
                vl_coords += [(0.05, (i+0.2) * lh), (0.4, (i+0.2) * lh), (0.4, ((i+0.2) + 1) * lh), (0.45, ((i+0.2) + 1) * lh), (0.05, ((i+0.2) + 1) * lh)]
            else:
                vl_coords += [(0.05, (i+0.2) * lh), (0.4, (i+0.2) * lh), (0.4, ((i+0.2) + 1) * lh), (0.05, ((i+0.2) + 1) * lh)]

        return (vl_coords, fl_indices)

    def draw(self, context):
        (r0h, r2w, r5h) = ret_dcoords(context)
        svp = context.scene.vi_params

        if self.expand:
            if self.resize:
                self.xdiff = self.lepos[0] - self.lspos[0]
                self.ydiff = self.lepos[1] - self.lspos[1]
            elif self.move:
                self.lspos[1] = self.lepos[1] - self.ydiff
                self.lepos[0] = self.lspos[0] + self.xdiff

            if self.lepos[1] > r5h - r0h:
                self.lspos[1] = r5h - r0h - self.ydiff
                self.lepos[1] = r5h - r0h

            if self.lepos[0] < r2w:
                self.lepos[0] = r2w

            self.base_shader.bind()
            self.base_shader.uniform_float("size", (self.xdiff, self.ydiff))
            self.base_shader.uniform_float("spos", self.lspos)
            self.base_shader.uniform_float("colour", self.hl)
            self.base_batch.draw(self.base_shader)
            self.basel_shader.bind()
            self.basel_shader.uniform_float("size", (self.xdiff, self.ydiff))
            self.basel_shader.uniform_float("spos", self.lspos)
            self.basel_shader.uniform_float("colour", [0, 0, 0, 1])
            self.basel_batch.draw(self.basel_shader)
            self.unit = svp.vi_leg_unit if svp.vi_leg_unit else res2unit[svp.li_disp_menu]

            if self.levels != svp.vi_leg_levels or self.cols != retcols(mcm.get_cmap(svp.vi_leg_col), self.levels) or (self.minres, self.maxres) != leg_min_max(svp):
                self.update(context)
                (vl_coords, fl_indices) = self.ret_coords()
                self.line_batch = batch_for_shader(self.line_shader, 'LINE_LOOP', {"position": vl_coords})
                self.col_batch = batch_for_shader(self.col_shader, 'TRIS', {"position": vl_coords, "colour": self.colours}, indices=fl_indices)

            self.col_shader.bind()
            self.col_shader.uniform_float("size", (self.xdiff, self.ydiff))
            self.col_shader.uniform_float("spos", self.lspos)
            self.col_batch.draw(self.col_shader)
            self.line_shader.bind()
            self.line_shader.uniform_float("size", (self.xdiff, self.ydiff))
            self.line_shader.uniform_float("spos", self.lspos)
            self.line_shader.uniform_float("colour", (0, 0, 0, 1))
            self.line_batch.draw(self.line_shader)
            tfontscale = max(self.titxdimen/(self.xdiff * 0.95), self.mydimen * 1.15/(self.lh * self.ydiff))
            blf.enable(self.font_id, 4)
            blf.enable(self.font_id, 8)
            blf.enable(self.font_id, blf.SHADOW)
            blf.shadow(self.font_id, 3, 0.7, 0.7, 0.7, 1)
            blf.size(self.font_id, int(12/tfontscale) * self.dpi/72)
            blf.position(self.font_id, self.lspos[0] + (self.xdiff - blf.dimensions(self.font_id, self.unit)[0]) * 0.45,
                         self.lepos[1] - 0.55 * (self.lh * self.ydiff) - blf.dimensions(self.font_id, 'M')[1] * 0.5, 0)
            blf.color(self.font_id, 0, 0, 0, 1)
            blf.draw(self.font_id, self.unit)
            # blf.disable(self.font_id, blf.SHADOW)
            lfontscale = max(self.resxdimen/(self.xdiff * 0.45), self.mydimen * 1.15/(self.lh * self.ydiff))
            blf.size(self.font_id, int(11/lfontscale) * self.dpi/72)

            for i in range(1, self.levels):
                num = self.resvals[i]
                ndimen = blf.dimensions(self.font_id, "{}".format(num))
                blf.position(self.font_id, int(self.lepos[0] - self.xdiff * 0.05 - ndimen[0]), int(self.lspos[1] + (i + 0.25) * self.lh * self.ydiff - ndimen[1]*0.5), 0)
                blf.draw(self.font_id, "{}".format(self.resvals[i]))

            blf.disable(self.font_id, blf.SHADOW)
            blf.disable(self.font_id, 8)
            blf.disable(self.font_id, 4)

    def create_batch(self):
        '''New code for cross-platform shaders'''
        base_shader = gpu.types.GPUShaderCreateInfo()
        base_shader.push_constant('MAT4', 'ModelViewProjectionMatrix')
        base_shader.push_constant('VEC2', 'spos')
        base_shader.push_constant('VEC2', 'size')
        base_shader.push_constant('VEC4', 'colour')
        base_shader.vertex_in(0, 'VEC2', 'position')

        base_shader.fragment_out(0, 'VEC4', 'FragColour')
        base_shader.vertex_source(
            '''
            void main()
            {
                float xpos = spos[0] + position[0] * size[0];
                float ypos = spos[1] + position[1] * size[1];
                gl_Position = ModelViewProjectionMatrix * vec4(int(xpos), int(ypos), -0.1f, 1.0f);
            }
            '''
        )
        base_shader.fragment_source(
            '''
            void main()
                {
                    FragColour = colour;
                }
            '''
        )

        col_shader = gpu.types.GPUShaderCreateInfo()
        col_shader.push_constant('MAT4', 'ModelViewProjectionMatrix')
        col_shader.push_constant('VEC2', 'spos')
        col_shader.push_constant('VEC2', 'size')
        col_shader.vertex_in(0, 'VEC2', 'position')
        col_shader.vertex_in(1, 'VEC4', 'colour')
        col_shader_iface = gpu.types.GPUStageInterfaceInfo('col_interface')
        col_shader_iface.flat('VEC4', 'f_colour')
        col_shader.vertex_out(col_shader_iface)
        col_shader.fragment_out(0, 'VEC4', 'FragColour')

        col_shader.vertex_source(
            '''
            void main()
            {
                float xpos = spos[0] + position[0] * size[0];
                float ypos = spos[1] + position[1] * size[1];
                gl_Position = ModelViewProjectionMatrix * vec4(int(xpos), int(ypos), -0.1f, 1.0f);
                f_colour = colour;
            }
            '''
        )
        col_shader.fragment_source(
            '''
            void main()
                {
                    FragColour = f_colour;
                }
            '''
        )

        self.base_shader = gpu.shader.create_from_info(base_shader)
        self.basel_shader = gpu.shader.create_from_info(base_shader)
        self.line_shader = gpu.shader.create_from_info(base_shader)
        self.col_shader = gpu.shader.create_from_info(col_shader)
        (vl_coords, fl_indices) = self.ret_coords()
        self.base_batch = batch_for_shader(self.base_shader, 'TRIS', {"position": self.v_coords}, indices=self.f_indices)
        self.basel_batch = batch_for_shader(self.basel_shader, 'LINE_STRIP', {"position": self.v_coords})
        self.line_batch = batch_for_shader(self.line_shader, 'LINE_LOOP', {"position": vl_coords})
        self.col_batch = batch_for_shader(self.col_shader, 'TRIS', {"position": vl_coords, "colour": self.colours}, indices=fl_indices)


def draw_dhscatter(self, x, y, z, tit, xlab, ylab, zlab, valmin, valmax, col):
    self.plt.close()
    x = [x[0] - 0.5] + [xval + 0.5 for xval in x]
    y = [y[0] - 0.5] + [yval + 0.5 for yval in y]
    self.fig, self.ax = plt.subplots(figsize=(12, 6))
    self.plt.title(tit, size=20).set_position([.5, 1.025])
    self.plt.xlabel(xlab, size=18)
    self.plt.ylabel(ylab, size=18)
    self.plt.pcolormesh(x, y, z, cmap=col, shading='auto', vmin=valmin, vmax=valmax)
    cbar = self.plt.colorbar(use_gridspec=True, pad=0.01)
    cbar.set_label(label=zlab, size=18)
    cbar.ax.tick_params(labelsize=16)
    self.plt.axis([min(x), max(x), min(y), max(y)])
    self.plt.xticks(size=16)
    self.plt.yticks(size=16)
    self.fig.tight_layout()


def save_plot(self, scene, filename):
    fileloc = os.path.join(scene.vi_params['viparams']['newdir'], 'images', filename)
    self.plt.savefig(fileloc, bbox_inches='tight')

    if filename not in [i.name for i in bpy.data.images]:
        self.gimage = filename
        bpy.data.images.load(fileloc)
    else:
        self.gimage = filename
        bpy.data.images[filename].reload()

    bpy.data.images[self.gimage].user_clear()


def show_plot(self, context):
    try:
        self.plt.show()
        self.update(context)

    except Exception as e:
        logentry('Error showing matplotlib graph: {}'.format(e))


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
                if d % 183 == 0:
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
        self.sumnorm = mathutils.Matrix().Rotation(pi/2, 4, 'X')@mathutils.Vector([0] + list((sumcoords[0]-sumcoords[1])[1:])).normalized()
        self.winnorm = mathutils.Matrix().Rotation(pi/2, 4, 'X')@mathutils.Vector([0] + list((wincoords[0]-wincoords[1])[1:])).normalized()

        for a, b in zip(coords[:-1], coords[1:]):
            line_lengths.append(line_lengths[-1] + (a - b).length)

        return (coords, line_lengths, breaks)

    def create_batch(self, scene, node):
        sp_shader = gpu.types.GPUShaderCreateInfo()
        sp_shader.push_constant('MAT4', 'viewProjectionMatrix')
        sp_shader.push_constant('MAT4', 'sp_matrix')
        sp_shader.push_constant('VEC4', 'colour1')
        sp_shader.push_constant('VEC4', 'colour2')
        sp_shader.push_constant('VEC4', 'colour3')
        sp_shader.push_constant('FLOAT', 'dash_density')
        sp_shader.push_constant('FLOAT', 'dash_ratio')
        sp_shader.vertex_in(0, 'VEC3', 'position')
        sp_shader.vertex_in(1, 'FLOAT', 'arcLength')
        sp_shader.vertex_in(2, 'UINT', 'line_break')
        sp_shader_iface = gpu.types.GPUStageInterfaceInfo('sp_interface')
        sp_shader_iface.smooth('VEC4', 'v_colour1')
        sp_shader_iface.smooth('VEC4', 'v_colour2')
        sp_shader_iface.smooth('VEC4', 'v_colour3')
        sp_shader_iface.smooth('FLOAT', 'v_ArcLength')
        sp_shader_iface.smooth('FLOAT', 'zpos')
        sp_shader_iface.flat('UINT', 'lb')
        sp_shader.vertex_out(sp_shader_iface)
        sp_shader.fragment_out(0, 'VEC4', 'FragColour')

        sp_shader.vertex_source(
            '''
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
        )
        sp_shader.fragment_source(
            '''
            void main()
                {
                    if (zpos < 0) {discard;}
                else if (lb == uint(1)) {discard;}
                else if (sin(v_ArcLength * dash_density) > dash_ratio) {FragColour = v_colour1;} else {FragColour = v_colour2;}
                if (lb == uint(2)) {FragColour = v_colour3;}
                }
            '''
        )

        sun_shader = gpu.types.GPUShaderCreateInfo()
        sun_shader.push_constant('MAT4', 'viewProjectionMatrix')
        sun_shader.push_constant('MAT4', 'sp_matrix')
        sun_shader.push_constant('VEC4', 'sun_colour')
        sun_shader.vertex_in(0, 'VEC3', 'position')
        sun_shader.fragment_out(0, 'VEC4', 'FragColour')
        sun_shader.vertex_source(
            '''
            void main()
            {
                gl_Position = viewProjectionMatrix * sp_matrix * vec4(position, 1.0f);
                    gl_Position[2] -= 0.001;
            }
            '''
        )
        sun_shader.fragment_source(
            '''
            void main()
                {
                   vec2 pos = gl_PointCoord - vec2(0.5);
                    if (length(pos) < 0.4) {FragColour = sun_colour;}
                    if (length(pos) <= 0.5) {
                        FragColour = sun_colour;
                        FragColour[3] = (0.5 - length(pos)) * 10;
                        }
                    if (length(pos) > 0.5) {discard;}
                }
            '''
        )

        globe_shader = gpu.types.GPUShaderCreateInfo()
        globe_shader.push_constant('MAT4', 'viewProjectionMatrix')
        globe_shader.push_constant('MAT4', 'sp_matrix')
        globe_shader.push_constant('VEC4', 'colour')
        globe_shader.vertex_in(0, 'VEC3', 'position')
        globe_shader.fragment_out(0, 'VEC4', 'FragColour')
        globe_shader.vertex_source(
            '''
            void main()
            {
                gl_Position = viewProjectionMatrix * sp_matrix * vec4(position, 1.0f);
            }
            '''
        )
        globe_shader.fragment_source(
            '''
            void main()
                {
                   FragColour = colour;
                }
            '''
        )

        range_shader = gpu.types.GPUShaderCreateInfo()
        range_shader.push_constant('MAT4', 'viewProjectionMatrix')
        range_shader.push_constant('MAT4', 'sp_matrix')
        range_shader.vertex_in(0, 'VEC3', 'position')
        range_shader.vertex_in(1, 'VEC3', 'colour')
        range_shader.fragment_out(0, 'VEC4', 'FragColour')
        range_shader_iface = gpu.types.GPUStageInterfaceInfo('range_interface')
        range_shader_iface.smooth('VEC3', 'tri_colour')
        range_shader.vertex_out(range_shader_iface)
        range_shader.vertex_source(
            '''
            void main()
            {
                gl_Position = viewProjectionMatrix * sp_matrix * vec4(position, 1.0f);
                tri_colour = colour;
            }
            '''
        )
        range_shader.fragment_source(
            '''
            void main()
                {
                   FragColour = vec4(tri_colour, 1.0);
                }
            '''
        )

        self.sp_shader = gpu.shader.create_from_info(sp_shader)
        self.sun_shader = gpu.shader.create_from_info(sun_shader)
        self.globe_shader = gpu.shader.create_from_info(globe_shader)
        self.range_shader = gpu.shader.create_from_info(range_shader)
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
        scene = context.scene
        svp = scene.vi_params

        try:
            # Draw lines
            gpu.state.depth_test_set('LESS')
            gpu.state.depth_mask_set(False)
            gpu.state.blend_set('ALPHA')
            gpu.state.line_width_set(svp.sp_line_width)
            gpu.state.point_size_set(svp.sp_sun_size)
            # bgl.glEnable(bgl.GL_DEPTH_TEST)
            # bgl.glDepthFunc(bgl.GL_LESS)
            # bgl.glDepthMask(bgl.GL_FALSE)
            # bgl.glEnable(bgl.GL_BLEND)
            # bgl.glLineWidth(context.scene.vi_params.sp_line_width)
            # bgl.glPointSize(context.scene.vi_params.sp_sun_size)
            # bgl.glEnable(bgl.GL_LINE_SMOOTH)
            # bgl.glEnable(bgl.GL_MULTISAMPLE)

            try:
                self.sp_shader.bind()
            except Exception:
                self.create_batch(context.scene, node)

            matrix = bpy.context.region_data.perspective_matrix
            sp_matrix = scene.objects['SPathMesh'].matrix_world
            sun_pos = [so.location[:] for so in scene.objects if so.type == 'LIGHT' and so.data.type == 'SUN' and not so.hide_viewport]
            self.sp_shader.uniform_float("viewProjectionMatrix", matrix)
            self.sp_shader.uniform_float("sp_matrix", sp_matrix)
            self.sp_shader.uniform_float("colour1", svp.sp_hour_dash)
            self.sp_shader.uniform_float("colour2", svp.sp_hour_main)
            self.sp_shader.uniform_float("colour3", svp.sp_season_main)
            self.sp_shader.uniform_float("dash_ratio", svp.sp_hour_dash_ratio)
            self.sp_shader.uniform_float("dash_density", svp.sp_hour_dash_density)
            self.sun_shader.bind()
            self.sun_shader.uniform_float("viewProjectionMatrix", matrix)
            self.sun_shader.uniform_float("sp_matrix", sp_matrix)
            self.sun_shader.uniform_float("sun_colour", svp.sp_sun_colour)
            self.globe_shader.bind()
            self.globe_shader.uniform_float("viewProjectionMatrix", matrix)
            self.globe_shader.uniform_float("sp_matrix", sp_matrix)
            self.globe_shader.uniform_float("colour", svp.sp_globe_colour)
            self.range_shader.bind()
            self.range_shader.uniform_float("viewProjectionMatrix", matrix)
            self.range_shader.uniform_float("sp_matrix", sp_matrix)

            if self.latitude != svp.latitude or self.longitude != svp.longitude or self.sd != svp.sp_sd or self.sh != svp.sp_sh or self.ss != svp.sp_sun_size:
                (coords, line_lengths, breaks) = self.ret_coords(scene, node)
                self.sp_batch = batch_for_shader(self.sp_shader, 'LINE_STRIP', {"position": coords, "arcLength": line_lengths, "line_break": breaks})
                sun_pos = [so.location[:] for so in scene.objects if so.type == 'LIGHT' and so.data.type == 'SUN' and not so.hide_viewport]
                self.sun_batch = batch_for_shader(self.sun_shader, 'POINTS', {"position": sun_pos})
                globe_v_coords, globe_f_indices = self.ret_globe_geometry(self.latitude, self.longitude)
                self.globe_batch = batch_for_shader(self.globe_shader, 'TRIS', {"position": globe_v_coords}, indices=globe_f_indices)
                range_v_coords, range_f_indices, range_col_indices = self.ret_range_geometry(self.latitude, self.longitude)
                self.range_batch = batch_for_shader(self.range_shader, 'TRIS', {"position": range_v_coords, "colour": range_col_indices})  # , indices=range_f_indices)
                self.latitude = svp.latitude
                self.longitude = svp.longitude
                self.sd = svp.sp_sd
                self.sh = svp.sp_sh
                self.ss = svp.sp_sun_size

            self.range_batch.draw(self.range_shader)
            self.globe_batch.draw(self.globe_shader)
            self.sun_batch.draw(self.sun_shader)
            self.sp_batch.draw(self.sp_shader)
            gpu.state.line_width_set(1)
            gpu.state.point_size_set(1)
            gpu.state.blend_set('NONE')
            gpu.state.depth_mask_set(True)
            gpu.state.depth_test_set('NONE')
            # bgl.glDisable(bgl.GL_MULTISAMPLE)
            # bgl.glDisable(bgl.GL_LINE_SMOOTH)
            # gpu.state.blend_set('NONE')
            # #bgl.glDisable(bgl.GL_BLEND)
            # bgl.glClear(bgl.GL_DEPTH_BUFFER_BIT)
            # bgl.glDisable(bgl.GL_DEPTH_TEST)
            # bgl.glDepthMask(bgl.GL_TRUE)
            # bgl.glPointSize(1)

        except Exception as e:
            logentry('Sun path error: {}'.format(e))
            svp.vi_display = 0

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
                    hs.append(str(int(co.split('-')[1]))+':00')
                    pos.append(view3d_utils.location_3d_to_region_2d(context.region,
                                                                     context.region_data,
                                                                     ob_mat@coords[co]))

                if pos:
                    draw_index(pos, hs, dists, svp.vi_display_rp_fs, svp.vi_display_rp_fc, svp.vi_display_rp_fsh, hour=1)

            if [ob.get('VIType') == 'Sun' for ob in bpy.data.objects if ob.parent == spob] and svp['spparams']['suns'] == '0':
                sobs = [ob for ob in bpy.data.objects if ob.get('VIType') == 'Sun' and ob.parent == spob]

                if sobs and svp.sp_td:
                    sunloc = ob_mat@sobs[0].location
                    solpos = view3d_utils.location_3d_to_region_2d(context.region, context.region_data, sunloc)

                    try:
                        if 0 < solpos[0] < width and 0 < solpos[1] < height and not scene.ray_cast(context.view_layer.depsgraph, sobs[0].location + 0.05 * (vl - sunloc), vl - sunloc)[0]:
                            soltime = datetime.datetime.fromordinal(svp.sp_sd)
                            soltime += datetime.timedelta(hours=svp.sp_sh)
                            sre = sobs[0].rotation_euler
                            blf_props(scene, width, height)
                            sol_text = soltime.strftime('     %d %b %X') + ' alt: {:.1f} azi: {:.1f}'.format(90 - sre[0]*180/pi, (180, -180)[sre[2] < -pi] - sre[2]*180/pi)
                            draw_time(solpos, sol_text, svp.vi_display_rp_fs,
                                      svp.vi_display_rp_fc, svp.vi_display_rp_fsh)

                    except Exception as e:
                        logentry("Sun path number error: {}".format(e))
            blf.disable(0, 4)
        else:
            return

    def ret_sun_geometry(self, dia, suns):
        sun_v_coords, sun_f_indices = [], []

        for sun in suns:
            sunbm = bmesh.new()
            bmesh.ops.create_uvsphere(sunbm, u_segments=12, v_segments=12, radius=dia, matrix=sun.matrix_world, calc_uvs=0)
            bmesh.ops.triangulate(sunbm, faces=sunbm.faces, quad_method='BEAUTY', ngon_method='BEAUTY')
            sun_v_coords += [v.co[:] for v in sunbm.verts]
            sun_f_indices += [[v.index for v in face.verts] for face in sunbm.faces]
            sunbm.free()
        return sun_v_coords, sun_f_indices

    def ret_globe_geometry(self, lat, long):
        midalt = solarPosition(79, 12, lat, long)[2]
        globebm = bmesh.new()
        altrot = mathutils.Matrix().Rotation(midalt, 4, 'X')
        bmesh.ops.create_uvsphere(globebm, u_segments=48, v_segments=48, radius=100,
                                  matrix=altrot, calc_uvs=0)
        bmesh.ops.bisect_plane(globebm, geom=globebm.verts[:] + globebm.edges[:] + globebm.faces[:], dist=0.01,
                               plane_co=(0, 0, 0), plane_no=(0, 0, 1), use_snap_center=False, clear_outer=False, clear_inner=True)
        bmesh.ops.bisect_plane(globebm, geom=globebm.verts[:] + globebm.edges[:] + globebm.faces[:], dist=0.1,
                               plane_co=self.winmid, plane_no=self.winnorm, use_snap_center=False,
                               clear_outer=True, clear_inner=False)
        bmesh.ops.bisect_plane(globebm, geom=globebm.verts[:] + globebm.edges[:] + globebm.faces[:], dist=0.1,
                               plane_co=self.summid, plane_no=self.sumnorm, use_snap_center=False,
                               clear_outer=False, clear_inner=True)
        bmesh.ops.triangulate(globebm, faces=globebm.faces, quad_method='BEAUTY', ngon_method='BEAUTY')
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
                    mornevediff = 360  # if bpy.context.scene.latitude >= 0 else 360

                startset = morn if lat >= 0 else eve
                angrange = [startset + a * 0.0125 * mornevediff for a in range(0, 81)]
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
        objmode()
        node = context.node
        scene.display.shadow_focus = 1
        svp = scene.vi_params
        svp.vi_display = 0
        svp['viparams'] = {}
        svp['spparams'] = {}
        svp['spparams']['suns'] = node.suns
        self.spcoll = create_coll(context, 'SunPath')
        context.view_layer.layer_collection.children[self.spcoll.name].exclude = 0
        sd = 100

        # Set the node colour
        node.export()

        svp['viparams']['resnode'], svp['viparams']['restree'] = node.name, node.id_data.name
        scene.cursor.location = (0.0, 0.0, 0.0)
        suns = [ob for ob in self.spcoll.objects if ob.type == 'LIGHT' and ob.data.type == 'SUN']
        requiredsuns = {'0': 1, '1': 12, '2': 24}[node.suns]
        matdict = {'SPBase': (0, 0, 0, 1), 'SPPlat': (1, 1, 1, 1), 'SPGrey': (0.25, 0.25, 0.25, 1)}

        for mat in [mat for mat in matdict if mat not in bpy.data.materials]:
            bpy.data.materials.new(mat)
            bpy.data.materials[mat].diffuse_color = matdict[mat][:4]
            bpy.data.materials[mat].use_nodes = True
            nodes = bpy.data.materials[mat].node_tree.nodes

            for n in nodes:
                nodes.remove(n)

            node_material = nodes.new(type='ShaderNodeBsdfDiffuse')
            node_material.inputs[0].default_value = matdict[mat]
            node_material.location = 0, 0
            node_output = nodes.new(type='ShaderNodeOutputMaterial')
            node_output.location = 400, 0
            links = bpy.data.materials[mat].node_tree.links
            links.new(node_material.outputs[0], node_output.inputs[0])

        if suns:
            for sun in suns[requiredsuns:]:
                bpy.data.objects.remove(sun, do_unlink=True, do_id_user=True, do_ui_user=True)

            suns = [ob for ob in context.scene.objects if ob.parent and ob.type == 'LIGHT' and ob.data.type == 'SUN' and ob.parent.get('VIType') == "SPathMesh"]
            [sun.animation_data_clear() for sun in suns]

        if not suns or len(suns) < requiredsuns:
            for rs in range(requiredsuns - len(suns)):
                bpy.ops.object.light_add(type='SUN', radius=1, align='WORLD', location=(0, 0, 0))
                suns.append(context.active_object)

        if scene.render.engine == 'CYCLES' and scene.world.get('node_tree') and 'Sky Texture' in [no.bl_label for no in scene.world.node_tree.nodes]:
            scene.world.node_tree.animation_data_clear()

        if bpy.context.active_object and not bpy.context.active_object.hide_viewport:
            if bpy.context.active_object.type == 'MESH':
                bpy.ops.object.mode_set(mode='OBJECT')

        if any(ob.get('VIType') == "SPathMesh" for ob in context.scene.objects):
            spathob = [ob for ob in context.scene.objects if ob.get('VIType') == "SPathMesh"][0]
        else:
            spathob = compass((0, 0, 0.01), sd, bpy.data.materials['SPPlat'], bpy.data.materials['SPBase'], bpy.data.materials['SPGrey'])

            if spathob.name not in self.spcoll.objects:
                self.spcoll.objects.link(spathob)
                if spathob.name in scene.collection.objects:
                    scene.collection.objects.unlink(spathob)

            spathob.location, spathob.name,  spathob['VIType'] = (0, 0, 0), "SPathMesh", "SPathMesh"
            selobj(context.view_layer, spathob)
            joinobj(context.view_layer, [spathob])
            spathob.visible_diffuse, spathob.visible_shadow, spathob.visible_glossy, spathob.visible_transmission, spathob.visible_volume_scatter = [False] * 5
            spathob.show_transparent = True

        for s, sun in enumerate(suns):
            if sun.name not in self.spcoll.objects:
                self.spcoll.objects.link(sun)
                scene.collection.objects.unlink(sun)

            sun.data.shadow_soft_size = 0.01
            sun['VIType'] = 'Sun'
            sun['solhour'], sun['solday'] = svp.sp_sh, svp.sp_sd
            sun.name = sun.data.name = 'Sun{}'.format(s)
            sun.parent = spathob

        if sunpath1 not in bpy.app.handlers.frame_change_post:
            bpy.app.handlers.frame_change_post.append(sunpath1)

        svp['viparams']['vidisp'] = 'sp'
        svp['viparams']['visimcontext'] = 'SunPath'
        sunpath(context)

        if context.screen:
                for a in context.screen.areas:
                    if a.type == 'VIEW_3D':
                        a.spaces[0].shading.shadow_intensity = svp.sp_sun_strength * 0.5

        self.suns = [sun for sun in scene.objects if sun.type == "LIGHT" and sun.data.type == 'SUN']
        self.sp = scene.objects['SPathMesh']
        self.latitude = svp.latitude
        self.longitude = svp.longitude
        self.sd = svp.sp_sd
        self.sh = svp.sp_sh
        self.ss = svp.sp_sun_size
        self.create_batch(scene, node)
        self.draw_handle_sp = bpy.types.SpaceView3D.draw_handler_add(self.draw_sp, (self, context, node), "WINDOW", "POST_VIEW")
        self.draw_handle_spnum = bpy.types.SpaceView3D.draw_handler_add(self.draw_spnum, (self, context), 'WINDOW', 'POST_PIXEL')
        bpy.app.driver_namespace["sp"] = self.draw_handle_sp
        bpy.app.driver_namespace["spnum"] = self.draw_handle_spnum
        svp['viparams']['drivers'] = ['sp', 'spnum']
        context.window_manager.modal_handler_add(self)
        svp.vi_display = 1
        rendview(1)
        return {'RUNNING_MODAL'}

    def modal(self, context, event):
        scene = context.scene
        svp = scene.vi_params

        if context.area:
            context.area.tag_redraw()

        if svp.vi_display == 0 or svp['viparams']['vidisp'] != 'sp' or not context.scene.objects.get('SPathMesh'):
            try:
                bpy.types.SpaceView3D.draw_handler_remove(self.draw_handle_sp, "WINDOW")
                bpy.types.SpaceView3D.draw_handler_remove(self.draw_handle_spnum, 'WINDOW')
            except Exception:
                pass

            svp.vi_display = 0

            for h in bpy.app.handlers.frame_change_post:
                bpy.app.handlers.frame_change_post.remove(h)

            [bpy.data.objects.remove(o, do_unlink=True, do_id_user=True, do_ui_user=True) for o in bpy.data.objects if o.vi_params.get('VIType') and o.vi_params['VIType'] in ('SunMesh', 'SkyMesh')]
            context.view_layer.layer_collection.children[self.spcoll.name].exclude = 1
            return {'CANCELLED'}

        return {'PASS_THROUGH'}


class VIEW3D_OT_WRDisplay(bpy.types.Operator):
    bl_idname = "view3d.wrdisplay"
    bl_label = "Wind rose number display"
    bl_description = "Project the windrose numbers on to the viewport"
    bl_register = True
    bl_undo = False

    def execute(self, context):
        (r0h, r2w, r5h) = ret_dcoords(context)
        scene = context.scene
        svp = scene.vi_params
        svp.vi_display = 1
        svp['viparams']['vidisp'] = 'wr'
        self.wt, self.scatcol, self.scattmax, self.scattmin, self.scattmaxval, self.scattminval, self.scattcol = svp.wind_type, 0, 0, 0, 0, 0, svp.vi_scatt_col
        self.images = ('legend.png', 'table.png')
        self.results_bar = results_bar(self.images)
        simnode = bpy.data.node_groups[svp['viparams']['restree']].nodes[svp['viparams']['resnode']]
        leg_unit = 'Speed (m/s)' if not simnode.temp else 'Temp (C)'
        self.legend = wr_legend(context, leg_unit, [305, r5h - r0h - 80], r2w, r5h - r0h, 100, 400)
        self.table = wr_table(context, [355, r5h - r0h - 80], r2w, r5h - r0h, 400, 60)
        self.cao = [o for o in scene.objects if o.vi_params.get('VIType') == "Wind_Plane"][0] if [o for o in scene.objects if o.vi_params.get('VIType') == "Wind_Plane"] else 0

        if not self.cao:
            return {'CANCELLED'}
        else:
            self.zdata = array(self.cao.vi_params[('d', 'wd')[int(svp.wind_type)]])
            self.zmax = nmax(self.zdata) if svp.vi_scatt_max == '0' else svp.vi_scatt_max_val
            self.zmin = nmin(self.zdata) if svp.vi_scatt_min == '0' else svp.vi_scatt_min_val
            self.xdata = self.cao.vi_params['days']
            self.ydata = self.cao.vi_params['hours']
            self.title = ('Wind Speed', 'Wind Direction')[int(svp.wind_type)]
            self.scatt_legend = ('Speed (m/s)', 'Direction ($^o$ deg from north)')[int(svp.wind_type)]
            self.xtitle = 'Days'
            self.ytitle = 'Hours'

        self.image = 'wr_scatter.png'
        self.zdata = array([])
        self.xtitle = 'Days'
        self.ytitle = 'Hours'
        self.height = r5h - r0h
        self.draw_handle_wrnum = bpy.types.SpaceView3D.draw_handler_add(self.draw_wrnum, (context, ), 'WINDOW', 'POST_PIXEL')
        bpy.app.driver_namespace["wr"] = self.draw_handle_wrnum
        context.area.tag_redraw()
        context.window_manager.modal_handler_add(self)
        return {'RUNNING_MODAL'}

    def modal(self, context, event):
        scene = context.scene
        svp = scene.vi_params
        redraw = 0

        if svp.vi_display == 0 or svp['viparams']['vidisp'] != 'wr' or not context.area:
            svp.vi_display = 0

            if context.area:
                context.area.tag_redraw()

            return {'CANCELLED'}

        if self.cao != context.active_object:
            if context.active_object and context.active_object.vi_params.get('VIType') == 'Wind_Plane':
                self.cao = context.active_object
                self.table.update(context)

        if event.type != 'INBETWEEN_MOUSEMOVE' and context.region and context.area.type == 'VIEW_3D' and context.region.type == 'WINDOW':
            mx, my = event.mouse_region_x, event.mouse_region_y

            for w, window in enumerate((self.legend, self.table)):             
                if self.results_bar.ipos[w][0][0] < mx < self.results_bar.ipos[w][1][0] and self.results_bar.ipos[w][0][1] < my < self.results_bar.ipos[w][2][1]:
                    window.hl = (0.8, 0.8, 0.8, 0.8)
                    redraw = 1

                    if event.type == 'LEFTMOUSE':
                        if event.value == 'RELEASE':
                            window.expand = 0 if window.expand else 1

                elif window.expand and abs(window.lspos[0] - mx) < 10 and abs(window.lepos[1] - my) < 10:
                    window.hl = (0.8, 0.8, 0.8, 0.8)
                    redraw = 1

                    if event.type == 'LEFTMOUSE':
                        if event.value == 'PRESS':
                            window.move = 1
                            window.draw(context)
                            context.area.tag_redraw()
                        elif window.move and event.value == 'RELEASE':
                            window.move = 0

                        return {'RUNNING_MODAL'}

                elif window.expand and abs(window.lepos[0] - mx) < 10 and abs(window.lspos[1] - my) < 10:
                    window.hl = (0.8, 0.8, 0.8, 0.8)
                    context.area.tag_redraw()

                    if event.type == 'LEFTMOUSE':
                        if event.value == 'PRESS':
                            window.resize = 1
                            window.draw(context)
                            context.area.tag_redraw()
                        elif window.resize and event.value == 'RELEASE':
                            window.resize = 0

                        return {'RUNNING_MODAL'}

                elif window.hl == (0.8, 0.8, 0.8, 0.8):
                    window.hl = (1, 1, 1, 1)
                    redraw = 1

                if event.type == 'MOUSEMOVE':
                    if window.move:
                        window.lspos[0], window.lepos[1] = mx, my
                        window.draw(context)
                        redraw = 1

                    elif window.resize:
                        window.lepos[0], window.lspos[1] = mx, my
                        window.draw(context)
                        redraw = 1
        if redraw:
            context.area.tag_redraw()

        return {'PASS_THROUGH'}

    def draw_wrnum(self, context):
        svp = context.scene.vi_params
        (r0h, r2w, r5h) = ret_dcoords(context)

        try:
            self.results_bar.draw(r2w, r5h - r0h)
            self.legend.draw(context)
            self.table.draw(context)

        except Exception as e:
            logentry("Something went wrong with wind rose display: {}".format(e))
            svp.vi_display == 0


class VIEW3D_OT_SVFDisplay(bpy.types.Operator):
    '''Display results legend and stats in the 3D View'''
    bl_idname = "view3d.svfdisplay"
    bl_label = "SVF display"
    bl_description = "Display sky view factor metrics"
    bl_register = True
    bl_undo = False

    def invoke(self, context, event):
        (r0h, r2w, r5h) = ret_dcoords(context)
        svp = context.scene.vi_params
        self.livi_coll = create_empty_coll(context, 'LiVi Results')
        context.view_layer.layer_collection.children[self.livi_coll.name].exclude = 0
        svp.vi_display = 1
        svp['viparams']['vidisp'] = 'svf'
        svp['viparams']['drivers'] = ['svf']
        self.simnode = bpy.data.node_groups[svp['viparams']['restree']].nodes[svp['viparams']['resnode']]
        
        if li_display(context, self, self.simnode) == 'CANCELLED':
            logentry('No result geometry present or visible')
            self.report({'ERROR'}, "No result geometry present or visible")
            svp.vi_display = 0
            return {'CANCELLED'}

        self.results_bar = results_bar(('legend.png',))
        legend_icon_pos = self.results_bar.ret_coords(r2w, r5h - r0h, 0)[0]
        self.legend = draw_legend(context, 'Sky View (%)', legend_icon_pos, r2w, r5h - r0h, 100, 400, 20)
        
        if svp.vi_disp_process != '2':
            self.legend_num = linumdisplay(self, context)
        
        self.height = r5h - r0h
        self.draw_handle_svfnum = bpy.types.SpaceView3D.draw_handler_add(self.draw_svfnum, (context, ), 'WINDOW', 'POST_PIXEL')
        bpy.app.driver_namespace["svf"] = self.draw_handle_svfnum
        self.cao = context.active_object
        context.region.tag_redraw()
        context.window_manager.modal_handler_add(self)
        return {'RUNNING_MODAL'}

    def draw_svfnum(self, context):
        (r0h, r2w, r5h) = ret_dcoords(context)
        svp = context.scene.vi_params

        try:
            self.results_bar.draw(r2w, r5h - r0h)
            self.legend.draw(context)
            
            if svp.vi_disp_process != '2':
                self.legend_num.draw(context)
                
        except Exception:
            svp.vi_display = 0

    def modal(self, context, event):
        scene = context.scene
        svp = scene.vi_params
        redraw = 0

        if svp.vi_display == 0 or svp['viparams']['vidisp'] != 'svf' or not context.area:
            svp.vi_display = 0
            move_obs(context.scene.collection, bpy.data.collections['LiVi Results'], 'LiVi Res')
            
            try:
                context.view_layer.layer_collection.children[self.livi_coll.name].exclude = 1
            except:
                pass

            if context.area:
                context.area.tag_redraw()

            return {'CANCELLED'}

        if event.type != 'INBETWEEN_MOUSEMOVE' and context.region and context.area.type == 'VIEW_3D' and context.region.type == 'WINDOW':
            mx, my = event.mouse_region_x, event.mouse_region_y

            # Legend routine
            if self.results_bar.ipos[0][0][0] < mx < self.results_bar.ipos[0][1][0] and self.results_bar.ipos[0][0][1] < my < self.results_bar.ipos[0][2][1]:
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


class VIEW3D_OT_RTDisplay(bpy.types.Operator):
    '''Display results legend and stats in the 3D View'''
    bl_idname = "view3d.rtdisplay"
    bl_label = "RT display"
    bl_description = "Display reverberation times"
    bl_register = True
    bl_undo = False

    def invoke(self, context, event):
        (r0h, r2w, r5h) = ret_dcoords(context)
        svp = context.scene.vi_params
        self.livi_coll = create_empty_coll(context, 'LiVi Results')
        context.view_layer.layer_collection.children[self.livi_coll.name].exclude = 0
        svp.vi_display = 1
        svp['viparams']['vidisp'] = 'rt'
        svp['viparams']['drivers'] = ['rt']
        self.simnode = bpy.data.node_groups[svp['viparams']['restree']].nodes[svp['viparams']['resnode']]
        
        if li_display(context, self, self.simnode) == 'CANCELLED':
            logentry('No result geometry present or visible')
            self.report({'ERROR'}, "No result geometry present or visible")
            svp.vi_display = 0
            return {'CANCELLED'}

        self.results_bar = results_bar(('legend.png',))
        legend_icon_pos = self.results_bar.ret_coords(r2w, r5h - r0h, 0)[0]
        self.legend = draw_legend(context, 'RT60 (%)', legend_icon_pos, r2w, r5h - r0h, 100, 400, 20)
        
        if svp.vi_disp_process != '2':
            self.legend_num = linumdisplay(self, context)
        
        self.height = r5h - r0h
        self.draw_handle_rtnum = bpy.types.SpaceView3D.draw_handler_add(self.draw_rtnum, (context, ), 'WINDOW', 'POST_PIXEL')
        bpy.app.driver_namespace["rt"] = self.draw_handle_rtnum
        self.cao = context.active_object
        context.region.tag_redraw()
        context.window_manager.modal_handler_add(self)
        return {'RUNNING_MODAL'}

    def draw_rtnum(self, context):
        (r0h, r2w, r5h) = ret_dcoords(context)
        svp = context.scene.vi_params

        try:
            self.results_bar.draw(r2w, r5h - r0h)
            self.legend.draw(context)

            if svp.vi_disp_process != '2':
                self.legend_num.draw(context)
 
        except Exception as e:
            logentry(f'Problem with number display: {e}')
            svp.vi_display = 0

    def modal(self, context, event):
        scene = context.scene
        svp = scene.vi_params
        redraw = 0

        if svp.vi_display == 0 or svp['viparams']['vidisp'] != 'rt' or not context.area:
            svp.vi_display = 0
            move_obs(context.scene.collection, bpy.data.collections['LiVi Results'], 'LiVi Res')
            
            try:
                context.view_layer.layer_collection.children[self.livi_coll.name].exclude = 1
            except:
                pass

            if context.area:
                context.area.tag_redraw()

            return {'CANCELLED'}

        if event.type != 'INBETWEEN_MOUSEMOVE' and context.region and context.area.type == 'VIEW_3D' and context.region.type == 'WINDOW':
            mx, my = event.mouse_region_x, event.mouse_region_y

            # Legend routine
            if self.results_bar.ipos[0][0][0] < mx < self.results_bar.ipos[0][1][0] and self.results_bar.ipos[0][0][1] < my < self.results_bar.ipos[0][2][1]:
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


class VIEW3D_OT_SSDisplay(bpy.types.Operator):
    '''Display results legend and stats in the 3D View'''
    bl_idname = "view3d.ssdisplay"
    bl_label = "Shadow study metric display"
    bl_description = "Display shadow study metrics"
    bl_register = True
    bl_undo = False

    def invoke(self, context, event):
        (r0h, r2w, r5h) = ret_dcoords(context)
        self.scene = context.scene
        region = context.region
        self.height = region.height
        self.width = region.width
        svp = context.scene.vi_params
        svp['viparams']['vidisp'] = 'ss'
        svp['viparams']['drivers'] = ['ss']
        self.livi_coll = create_empty_coll(context, 'LiVi Results')
        self.cao = context.active_object
        self.image = 'ss_scatter.png'
        self.frame = self.scene.frame_current
        self.zmax = 100
        self.zmin = 0
        self.scattmax, self.scattmin, self.scattmaxval, self.scattminval, self.scattcol = 0, 0, 0, 0, svp.vi_scatt_col
        self.title = 'Total Area Sunlit'
        self.scatt_legend = 'Area (%)'
        self.xtitle = 'Days'
        self.ytitle = 'Hours'
        self.simnode = bpy.data.node_groups[svp['viparams']['restree']].nodes[svp['viparams']['resnode']]
        
        if li_display(context, self, self.simnode) == 'CANCELLED':
            logentry('No result geometry present or visible')
            self.report({'ERROR'}, "No result geometry present or visible")
            svp.vi_display = 0
            return {'CANCELLED'}
        
        self.images = ('legend.png', )
        self.results_bar = results_bar(self.images)
        legend_icon_pos = self.results_bar.ret_coords(r2w, r5h - r0h, 0)[0]
        self.legend = draw_legend(context, 'Sunlit (%)', legend_icon_pos, r2w, r5h - r0h, 100, 400, 20)

        if svp.vi_disp_process != '2':
            self.num_display = linumdisplay(self, context)

        svp.vi_disp_wire = 1
        self.draw_handle_ssnum = bpy.types.SpaceView3D.draw_handler_add(self.draw_ssnum, (context, ), 'WINDOW', 'POST_PIXEL')
        bpy.app.driver_namespace["ss"] = self.draw_handle_ssnum
        context.window_manager.modal_handler_add(self)
        svp.vi_display = 1
        return {'RUNNING_MODAL'}

    def modal(self, context, event):
        scene = context.scene
        svp = scene.vi_params
        redraw = 0
        updates = [0 for i in self.images]

        if svp.vi_display == 0 or svp['viparams']['vidisp'] != 'ss' or not [o for o in bpy.data.objects if o.vi_params.vi_type_string == 'LiVi Calc'] or not context.area:
            svp.vi_display = 0
            move_obs(context.scene.collection, bpy.data.collections['LiVi Results'], 'LiVi Res')
            context.view_layer.layer_collection.children[self.livi_coll.name].exclude = 1

            if context.area:
                context.area.tag_redraw()

            return {'CANCELLED'}

        if self.scattcol != svp.vi_leg_col or self.cao != context.active_object or self.frame != svp.vi_frames:
            self.cao = context.active_object
            self.scattcol = svp.vi_leg_col
            self.frame = svp.vi_frames

            if self.cao and self.cao.vi_params.get('ss{}'.format(self.frame)):
                self.zdata = array(context.active_object.vi_params['ss{}'.format(self.frame)])
                self.xdata = array(self.cao.vi_params['days'])
                self.ydata = array(self.cao.vi_params['hours'])

        if event.type != 'INBETWEEN_MOUSEMOVE' and context.region and context.area.type == 'VIEW_3D' and context.region.type == 'WINDOW':
            mx, my = event.mouse_region_x, event.mouse_region_y

            for w, window in enumerate((self.legend, )):
                if self.results_bar.ipos[w][0][0] < mx < self.results_bar.ipos[w][1][0] and self.results_bar.ipos[w][0][1] < my < self.results_bar.ipos[w][2][1]:
                    window.hl = (0.8, 0.8, 0.8, 0.8)
                    redraw = 1
                    if event.type == 'LEFTMOUSE':
                        if event.value == 'RELEASE':
                            window.expand = 0 if window.expand else 1
                            context.area.tag_redraw()
                            return {'RUNNING_MODAL'}

                elif w == 1 and window.expand and window.lspos[0] + 0.1 * window.xdiff < mx < window.lepos[0] - 0.1 * window.xdiff and window.lspos[1] + \
                        0.1 * window.ydiff < my < window.lepos[1] - 0.1 * window.ydiff:
                    window.hl = (0.8, 0.8, 0.8, 0.8)
                    redraw = 1
                    if event.type == 'LEFTMOUSE' and event.value == 'RELEASE':
                        window.show_plot(context)

                elif window.expand and abs(window.lspos[0] - mx) < 10 and abs(window.lepos[1] - my) < 10:
                    window.hl = (0.8, 0.8, 0.8, 0.8)
                    redraw = 1
                    if event.type == 'LEFTMOUSE':
                        if event.value == 'PRESS':
                            window.move = 1
                            window.draw(context)
                            context.area.tag_redraw()
                        elif window.move and event.value == 'RELEASE':
                            window.move = 0
                        return {'RUNNING_MODAL'}

                elif window.expand and abs(window.lepos[0] - mx) < 10 and abs(window.lspos[1] - my) < 10:
                    window.hl = (0.8, 0.8, 0.8, 0.8)
                    context.area.tag_redraw()
                    if event.type == 'LEFTMOUSE':
                        if event.value == 'PRESS':
                            window.resize = 1
                            window.draw(context)
                            context.area.tag_redraw()
                        elif window.resize and event.value == 'RELEASE':
                            window.resize = 0
                        return {'RUNNING_MODAL'}

                elif window.hl == (0.8, 0.8, 0.8, 0.8):
                    window.hl = (1, 1, 1, 1)
                    redraw = 1

                if event.type == 'MOUSEMOVE':
                    if window.move:
                        window.lspos[0], window.lepos[1] = mx, my
                        window.draw(context)
                        context.area.tag_redraw()

                    elif window.resize:
                        window.lepos[0], window.lspos[1] = mx, my
                        window.draw(context)
                        context.area.tag_redraw()

        if redraw:
            context.area.tag_redraw()
            return {'RUNNING_MODAL'}

        return {'PASS_THROUGH'}

    def draw_ssnum(self, context):
        (r0h, r2w, r5h) = ret_dcoords(context)
        svp = context.scene.vi_params

        try:
            self.results_bar.draw(r2w, r5h - r0h)
            self.legend.draw(context)

            if svp.vi_disp_process != '2':
                self.num_display.draw(context)

        except Exception:
            svp.vi_display = 0


class VIEW3D_OT_Li_DBSDF(bpy.types.Operator):
    bl_idname = "view3d.bsdf_display"
    bl_label = "BSDF display"
    bl_description = "Display BSDF"
    bl_register = True
    bl_undo = False

    def modal(self, context, event):
        scene = context.scene
        svp = scene.vi_params
        (r0h, r2w, r5h) = ret_dcoords(context)
        redraw = 0

        if svp.vi_display == 0 or svp['viparams']['vidisp'] != 'bsdf_panel' or not context.area:
            svp['viparams']['vidisp'] = 'bsdf'
            svp.vi_display = 0
            self.bsdf.plt.close()

            if context.area:
                context.area.tag_redraw()

            return {'CANCELLED'}

        if self.bsdf.expand and any((self.bsdf.leg_max != svp.vi_bsdfleg_max,
                                     self.bsdf.leg_min != svp.vi_bsdfleg_min,
                                     self.bsdf.col != svp.vi_leg_col,
                                     self.bsdf.scale_select != svp.vi_bsdfleg_scale,
                                     self.bsdf.direc != svp.vi_bsdf_direc,
                                     self.bsdf.num_disp != svp.vi_bsdf_font)):
            self.bsdf.col = svp.vi_leg_col
            self.bsdf.leg_max = svp.vi_bsdfleg_max
            self.bsdf.leg_min = svp.vi_bsdfleg_min
            self.bsdf.scale_select = svp.vi_bsdfleg_scale
            self.bsdf.direc = svp.vi_bsdf_direc
            self.bsdf.num_disp = svp.vi_bsdf_font
            self.bsdf.plot(context)
            context.region.tag_redraw()
            return {'PASS_THROUGH'}

        if event.type != 'INBETWEEN_MOUSEMOVE' and context.region and context.area.type == 'VIEW_3D' and context.region.type == 'WINDOW':
            mx, my = event.mouse_region_x, event.mouse_region_y
            rbxl = self.results_bar.ret_coords(r2w, r5h - r0h, 0)

            if rbxl[0][0] < mx < rbxl[1][0] and rbxl[0][1] < my < rbxl[2][1]:
                self.bsdf.hl = (0.8, 0.8, 0.8, 0.8)
                redraw = 1

                if event.type == 'LEFTMOUSE':
                    if event.value == 'RELEASE':
                        self.bsdf.expand = 0 if self.bsdf.expand else 1

                context.region.tag_redraw()
                return {'RUNNING_MODAL'}

            elif self.bsdf.expand:
                self.bsdf.lspos = [self.results_bar.ret_coords(r2w, r5h - r0h, 0)[0][0] - 5, self.results_bar.ret_coords(r2w, r5h - r0h, 0)[0][1] - self.bsdf.ydiff - 25]
                self.bsdf.lepos = [self.bsdf.lspos[0] + self.bsdf.xdiff, self.bsdf.lspos[1] + self.bsdf.ydiff]

                if ((mx - (self.bsdf.lspos[0] + 50 + 125))**2 + (my - self.bsdf.lepos[1] + 155)**2)**0.5 <= 124:
                    dist = ((mx - (self.bsdf.lspos[0] + 50 + 125))**2 + (my - self.bsdf.lepos[1] + 155)**2)**0.5

                    for radius in self.bsdf.radii:
                        if dist <= radius:
                            redraw = 1
                            mring = self.bsdf.radii.index(radius)
                            mangle = atan2(-my + (self.bsdf.lepos[1] - 155), mx - (self.bsdf.lspos[0] + 50 + 125)) + pi + pi/self.bsdf.segments[mring]

                            if mring == 0:
                                ms = 1
                            else:
                                ms = int(mangle*self.bsdf.segments[mring]/(2*pi)) + 1 if int(mangle*self.bsdf.segments[mring]/(2*pi)) < self.bsdf.segments[mring] else 1

                            msegment = sum([self.bsdf.segments[ri] for ri in range(mring)]) + ms
                            break

                    self.bsdf.cr, self.bsdf.crs, self.bsdf.cseg = mring, ms, msegment
                    self.bsdf.create_batch('sel')

                    if event.type == 'LEFTMOUSE':
                        redraw = 1

                        if self.bsdf.cseg != self.bsdf.sseg:
                            self.bsdf.sr, self.bsdf.srs, self.bsdf.sseg = mring, ms, msegment
                            self.bsdf.plot(context)

                elif (self.bsdf.imspos[0] < mx < self.bsdf.imspos[0] + 400) and (self.bsdf.imspos[1] < my < self.bsdf.imspos[1] + 400):
                    if event.type == 'LEFTMOUSE':
                        if event.value == 'RELEASE':
                            self.bsdf.plt.show()
                            return {'RUNNING_MODAL'}

                if redraw:
                    context.region.tag_redraw()
                    return {'RUNNING_MODAL'}

        return {'PASS_THROUGH'}

    def invoke(self, context, event):
        cao = context.active_object
        caomvp = cao.active_material.vi_params
        (r0h, r2w, r5h) = ret_dcoords(context)
        scene = context.scene
        svp = scene.vi_params

        if cao and caomvp.get('bsdf') and caomvp['bsdf']['xml'] and caomvp['bsdf']['type'] == 'LBNL/Klems Full':
            bsdf = parseString(caomvp['bsdf']['xml'])
            svp['liparams']['bsdf_direcs'] = [(path.firstChild.data, path.firstChild.data, 'BSDF Direction') for path in bsdf.getElementsByTagName('WavelengthDataDirection')]
            self.images = ['bsdf.png']
            self.results_bar = results_bar(self.images)
            self.bsdf = draw_bsdf(context, '', self.results_bar.ret_coords(r2w, r5h - r0h, 0)[0], r2w, r5h - r0h, 400, 650)
            svp.vi_display = 1
            self.draw_handle_bsdfnum = bpy.types.SpaceView3D.draw_handler_add(self.draw_bsdfnum, (context, ), 'WINDOW', 'POST_PIXEL')
            bpy.app.driver_namespace["bsdf"] = self.draw_handle_bsdfnum
            context.window_manager.modal_handler_add(self)
            svp['viparams']['vidisp'] = 'bsdf_panel'
            svp['viparams']['drivers'] = ['bsdf']
            context.area.tag_redraw()
            return {'RUNNING_MODAL'}
        else:
            self.report({'ERROR'}, "Selected material contains no BSDF information or contains the wrong BSDF type (only Klems is supported)")
            return {'CANCELLED'}

    def draw_bsdfnum(self, context):
        (r0h, r2w, r5h) = ret_dcoords(context)

        try:
            self.results_bar.draw(r2w, r5h - r0h)
            self.bsdf.draw(context)
        except Exception:
            context.scene.vi_params.vi_display = 0


class VIEW3D_OT_Li_BD(bpy.types.Operator):
    '''Display results legend and stats in the 3D View'''
    bl_idname = "view3d.libd"
    bl_label = "LiVi metric display"
    bl_description = "Display lighting metrics"
    bl_register = True
    bl_undo = False

    def modal(self, context, event):
        scene = context.scene
        svp = scene.vi_params
        redraw = 0

        if svp.vi_display == 0 or not context.area or svp['viparams']['vidisp'] != 'li' or not [o for o in context.scene.objects if o.vi_params.vi_type_string == 'LiVi Res']:
            if context.active_object and context.active_object.mode == 'EDIT':
                bpy.ops.object.mode_set(mode='OBJECT')

            svp.vi_display = 0
            move_obs(context.scene.collection, bpy.data.collections['LiVi Results'], 'LiVi Res')

            if context.area:
                context.area.tag_redraw()
            else:
                logentry('You have encountered a Blender bug: "internal error: modal gizmo-map handler has invalid area". \
                          Do not maximise a window while the display operator is running.')

            context.view_layer.layer_collection.children['LiVi Results'].exclude = 1
            return {'CANCELLED'}

        if event.type != 'INBETWEEN_MOUSEMOVE' and context.region and context.area.type == 'VIEW_3D' and context.region.type == 'WINDOW':
            (r0h, r2w, r5h) = ret_dcoords(context)
            mx, my = event.mouse_region_x, event.mouse_region_y

            # Legend routine
            rbxl = self.results_bar.ret_coords(r2w, r5h - r0h, 0)

            if rbxl[0][0] < mx < rbxl[1][0] and rbxl[0][1] < my < rbxl[2][1]:
                self.legend.hl = (0.8, 0.8, 0.8, 0.8)
                redraw = 1

                if event.type == 'LEFTMOUSE':
                    if event.value == 'RELEASE':
                        self.legend.expand = 0 if self.legend.expand else 1
                        if context.area:
                            context.area.tag_redraw()
                        return {'RUNNING_MODAL'}

            elif self.legend.expand and abs(self.legend.lspos[0] - mx) < 10 and abs(self.legend.lepos[1] - my) < 10:
                self.legend.hl = (0.8, 0.8, 0.8, 0.8)
                redraw = 1

                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.legend.move = 1
                        self.legend.draw(context)

                        if context.area:
                            context.area.tag_redraw()

                    elif self.legend.move and event.value == 'RELEASE':
                        self.legend.move = 0

                    return {'RUNNING_MODAL'}

            elif self.legend.expand and abs(self.legend.lepos[0] - mx) < 10 and abs(self.legend.lspos[1] - my) < 10:
                self.legend.hl = (0.8, 0.8, 0.8, 0.8)

                if context.area:
                    context.area.tag_redraw()

                if event.type == 'LEFTMOUSE':
                    if event.value == 'PRESS':
                        self.legend.resize = 1
                        self.legend.draw(context)

                        if context.area:
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
                    redraw = 1
                elif self.legend.resize:
                    self.legend.lepos[0], self.legend.lspos[1] = mx, my
                    self.legend.draw(context)
                    redraw = 1

            if redraw:
                if context.area:
                    context.area.tag_redraw()

        return {'PASS_THROUGH'}

    def invoke(self, context, event):
        (r0h, r2w, r5h) = ret_dcoords(context)
        self.scene = context.scene
        svp = context.scene.vi_params
        svp.vi_display, svp.vi_disp_wire = 1, 1
        clearscene(context, self)
        svp['viparams']['vidisp'] = 'li'
        svp['viparams']['drivers'] = ['li']
        self.simnode = bpy.data.node_groups[svp['viparams']['restree']].nodes[svp['viparams']['resnode']]
        self.images = ['legend.png']
        self.results_bar = results_bar(self.images)
        self.frame = self.scene.frame_current

        if svp.vi_res_process == '2' and svp.script_file in bpy.data.texts and 'resmod' not in bpy.app.driver_namespace.keys():
            script = bpy.data.texts[svp.script_file]
            exec(script.as_string())

        if li_display(context, self, self.simnode) == 'CANCELLED':
            logentry('No result geometry present or visible')
            self.report({'ERROR'}, "No result geometry present or visible")
            svp.vi_display = 0
            return {'CANCELLED'}

        self.legend = draw_legend(context, svp['liparams']['unit'], self.results_bar.ret_coords(r2w, r5h - r0h, 0)[0], r2w, r5h - r0h, 75, 400, 20)
        
        if svp.vi_disp_process != '2':
            self.legend_num = linumdisplay(self, context)
        
        self.draw_handle_linum = bpy.types.SpaceView3D.draw_handler_add(self.draw_linum, (context, ), 'WINDOW', 'POST_PIXEL')
        bpy.app.driver_namespace["li"] = self.draw_handle_linum
        context.window_manager.modal_handler_add(self)
        return {'RUNNING_MODAL'}

    def draw_linum(self, context):
        (r0h, r2w, r5h) = ret_dcoords(context)
        svp = context.scene.vi_params
        
        try:
            self.results_bar.draw(r2w, r5h - r0h)
            self.legend.draw(context)
            
            if svp.vi_disp_process != '2':
                self.legend_num.draw(context)

        except Exception as e:
            logentry('Quitting LiVi display {}'.format(e))
            svp.vi_display = 0


class NODE_OT_Vi_Info(bpy.types.Operator):
    '''Display result infographics'''
    bl_idname = "node.vi_info"
    bl_label = "Inforgraphic display"
    bl_description = "Display infographics"
    bl_register = True
    bl_undo = False

    def execute(self, context):
        scene = context.scene
        svp = scene.vi_params
        dim = 800
        node = context.node

        if node.metric == '1' and node.light_menu == '0':
            dim = (800, 800)
            imname, svg_bytes = vi_info(node, dim, svp, ir=node['res']['ratioDF'], aDF=node['res']['avDF'], rDF=node['res']['rDF'], creds=node['res']['b_creds'], area=node['res']['b_area'], areas=node['res']['b_areas'])
        elif node.metric == '1' and node.light_menu == '2':
            dim = (800, 800)
            imname, svg_bytes = vi_info(node, dim, svp, ir=node['res']['ratioDF'], aDF=node['res']['avDF'])
        elif node.metric == '1' and node.light_menu == '1':
            dim = (600, 800)
            imname, svg_bytes = vi_info(node, dim, svp, sda=node['res']['sda'], sdapass=node['res']['sdapass'],
                                        ase=node['res']['ase'], asepass=node['res']['asepass'], o1=node['res']['o1'],
                                        tc=node['res']['tc'], totarea=node['res']['totarea'], svarea=node['res']['svarea'])

        elif node.metric == '6':
            dim = (1000, 1000)
            imname, svg_bytes = vi_info(node, dim, svp, wlc=node['res']['wl'], ec=node['res']['ec'], noc=node['res']['noc'], oc=node['res']['oc'], of=node['res']['of'])
            # qtsvg = QSvgRenderer()
            # qtsvg.load(svg_bytes)
            # printer = QPdfWriter(os.path.join(svp['viparams']['newdir'], "WLC.pdf"))
            # #printer.setPageSize(QPageSize(QPageSize.PageSizeId.A4))
            # # printer.setPageMargins(QMarginsF(0, 0, 0, 0))
            # printer.setResolution(300)
            # painter = QPainter(printer)
            # qtsvg.render(painter)
            # painter.end()

        with open(os.path.join(svp['viparams']['newdir'], "metric.svg"), 'w') as metric_file:
            metric_file.write(svg_bytes.decode())

        image = QImage.fromData(svg_bytes)
        image = image.convertToFormat(QImage.Format.Format_RGBA8888)
        image = image.mirrored(0, 1)
        buf = image.bits()

        if buf:
            arr = frombuffer(buf, dtype=ubyte).astype(float32)
            ipwidth, ipheight = dim[0], dim[1]

            if imname not in [im.name for im in bpy.data.images]:
                bpy.ops.image.new(name=imname, width=ipwidth, height=ipheight, color=(0, 0, 0, 0), alpha=True,
                                generated_type='BLANK', float=False, use_stereo_3d=False)
                im = bpy.data.images[imname]
                im.colorspace_settings.name = 'sRGB'

            else:
                im = bpy.data.images[imname]
                im.gl_free()
                im.buffers_free()

                if (im.generated_width, im.generated_height) != (ipwidth, ipheight):
                    im.generated_width = ipwidth
                    im.generated_height = ipheight

                if im.size[:] != (ipwidth, ipheight):
                    im.scale(ipwidth, ipheight)

            im.pixels.foreach_set((arr/255))
            im.scale(int(dim[0]), int(dim[1]))
            area = context.area
            t = area.type
            area.type = 'IMAGE_EDITOR'
            area.spaces.active.image = im
            bpy.ops.image.view_all()
            bpy.ops.screen.area_dupli('INVOKE_DEFAULT')
            win = bpy.context.window_manager.windows[-1]
            win.screen.areas[0].spaces[0].show_region_header = 0
            win.screen.areas[0].spaces[0].show_region_toolbar = 0
            win.screen.areas[0].spaces[0].show_region_ui = 0
            win.screen.areas[0].spaces[0].show_gizmo = 0
            win.screen.areas[0].spaces[0].use_realtime_update = 0

            with context.temp_override(window=win, area=win.screen.areas[0]):
                bpy.ops.image.view_zoom_ratio(ratio=1)

            area.type = t
        else:
            self.report({'ERROR'}, "No image data")
            return {'CANCELLED'}
        return {'FINISHED'}


