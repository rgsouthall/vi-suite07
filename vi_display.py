import bpy, blf, mathutils, datetime, os, bgl, inspect, gpu, bmesh, random
from gpu_extras.batch import batch_for_shader
from mathutils import Vector
from bpy_extras import view3d_utils
from .vi_func import ret_vp_loc, viewdesc, draw_index, draw_time, blf_props, retcols, drawloop, drawpoly, leg_min_max, retdp, xy2radial
from .vi_func import ret_res_vals, draw_index_distance, setscenelivivals, selobj
from .vi_func import logentry, move_to_coll, cmap, retvpvloc, objmode, skframe, clearscene
from .vi_func import solarPosition, solarRiseSet, create_coll, compass, joinobj, sunpath
from .livi_export import spfc
from .vi_dicts import unit_dict
from . import livi_export
from math import pi, log10, atan2, sin, cos, atan2
from numpy import array, repeat
from numpy import min as nmin
from numpy import max as nmax
from numpy import sum as nsum
from numpy import append as nappend
from xml.dom.minidom import parse, parseString

try:
    import matplotlib
    matplotlib.use('qt5agg', warn = False, force = True)
    import matplotlib.pyplot as plt
    import matplotlib.cm as mcm  
    import matplotlib.colors as mcolors
    from matplotlib.patches import Rectangle
    from matplotlib.collections import PatchCollection
    mp = 1
except:
    mp = 0

kfsa = array([0.02391, 0.02377, 0.02341, 0.02738, 0.02933, 0.03496, 0.04787, 0.05180, 0.13552])
kfact = array([0.9981, 0.9811, 0.9361, 0.8627, 0.7631, 0.6403, 0.4981, 0.3407, 0.1294])

def rendview(i):
    for scrn in bpy.data.screens:
        if scrn.name == 'Default':
            bpy.context.window.screen = scrn
            for area in scrn.areas:
                if area.type == 'VIEW_3D':
                    for space in area.spaces:
                        if space.type == 'VIEW_3D':
                            space.clip_start = 0.1
                            bpy.context.scene['cs'] = space.clip_start
                            
def li_display(disp_op, simnode):    
    scene, obreslist, obcalclist = bpy.context.scene, [], []
    svp = scene.vi_params
    svp['liparams']['livir'] = []
    setscenelivivals(scene)
    try:
        scene.display_settings.display_device = 'None'
    except:
        pass
    (rcol, mtype) =  ('hot', 'livi') if 'LiVi' in simnode.bl_label else ('grey', 'shad')

    for geo in scene.objects:
        bpy.context.view_layer.objects.active = geo
        
        if getattr(geo, 'mode') != 'OBJECT':
            bpy.ops.object.mode_set(mode = 'OBJECT')

    bpy.ops.object.select_all(action = 'DESELECT')

    if not bpy.app.handlers.frame_change_post:
        bpy.app.handlers.frame_change_post.append(livi_export.cyfc1)
        
    for o in scene.objects:
        if o.type == "MESH" and o.get('licalc') and o.hide == False:
            bpy.ops.object.select_all(action = 'DESELECT')
            obcalclist.append(o)
    
    scene.frame_set(svp['liparams']['fs'])
    bpy.context.view_layer.objects.active = None
    
    for i, o in enumerate([scene.objects[oname] for oname in svp['liparams']['{}c'.format(mtype)]]):        
        bm = bmesh.new()
        tempmesh = o.to_mesh()
        bm.from_mesh(tempmesh)
        o.to_mesh_clear() 
        ovp = o.vi_params
#        bm.normal_update()
#        bm.transform(o.matrix_world)
#        
                 
        if svp['liparams']['cp'] == '0':  
            cindex = bm.faces.layers.int['cindex']
            for f in [f for f in bm.faces if f[cindex] < 1]:
                bm.faces.remove(f)
            [bm.verts.remove(v) for v in bm.verts if not v.link_faces]

        elif svp['liparams']['cp'] == '1':
            cindex =  bm.verts.layers.int['cindex']
            for v in [v for v in bm.verts if v[cindex] < 1]:
                bm.verts.remove(v)
            for v in bm.verts:
                v.select = True
        
        while bm.verts.layers.shape:
            bm.verts.layers.shape.remove(bm.verts.layers.shape[-1])
        
        for v in bm.verts:
            v.co += mathutils.Vector((nsum([f.normal for f in v.link_faces], axis = 0))).normalized()  * simnode['goptions']['offset']
        
        selobj(bpy.context.view_layer, o)
        bpy.ops.object.duplicate() 
        
        for face in bm.faces:
            face.select = True 
        
        if not bpy.context.active_object:
            disp_op.report({'ERROR'},"No display object. If in local view switch to global view and/or re-export the geometry")
            return 'CANCELLED'
            
        ores = bpy.context.active_object
        ores.name, ores.show_wire, ores.display_type, orvp = o.name+"res", 1, 'SOLID', ores.vi_params
        move_to_coll(bpy.context, 'LiVi Results', ores)
        
        while ores.material_slots:
            bpy.ops.object.material_slot_remove()
        
        while ores.data.shape_keys:
            bpy.context.object.active_shape_key_index = 0
            bpy.ops.object.shape_key_remove(all=True)
            
        cv = ores.cycles_visibility
        cv.diffuse, cv.glossy, cv.transmission, cv.scatter, cv.shadow = 0, 0, 0, 0, 0        
        obreslist.append(ores)
        orvp['omax'], orvp['omin'], orvp['oave'] = ovp['omax'], ovp['omin'], ovp['oave'] 
        selobj(bpy.context.view_layer, ores)
        cmap(svp)
        
        for matname in ['{}#{}'.format('vi-suite', i) for i in range(svp.vi_leg_levels + 1)]:
            if bpy.data.materials[matname] not in ores.data.materials[:]:
                bpy.ops.object.material_slot_add()
                ores.material_slots[-1].material = bpy.data.materials[matname]
        
        if svp.vi_disp_3d == 1 and svp['liparams']['cp'] == '0':
            bm.faces.layers.int.new('extrude')
            extrude = bm.faces.layers.int['extrude']
            for face in bmesh.ops.extrude_discrete_faces(bm, faces = bm.faces)['faces']:
                face.select = True
                face[extrude] = 1
                
#        bm.transform(o.matrix_world.inverted())
        bm.to_mesh(ores.data)
        bm.free()
        bpy.ops.object.shade_flat()        
        ores.vi_params.lividisplay(scene)
                
        if svp.vi_disp_3d == 1 and ores.data.shape_keys == None:
            selobj(bpy.context.view_layer, ores)
            bpy.ops.object.shape_key_add(from_mix = False)
            
            for frame in range(svp['liparams']['fs'], svp['liparams']['fe'] + 1):
                bpy.ops.object.shape_key_add(from_mix = False)
                ores.active_shape_key.name, ores.active_shape_key.value = str(frame), 1
                
    svp['liparams']['livir'] = [o.name for o in obreslist]            
    skframe('', scene, obreslist)                                   
    bpy.ops.wm.save_mainfile(check_existing = False)
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
        self.obreslist = [ob for ob in scene.objects if ob.name in svp['liparams']['livir']]

        if svp.vi_display_sel_only == False:
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
        self.fontmult = 2 #if context.space_data.region_3d.is_perspective else 500
        
        if not svp.get('viparams') or svp['viparams']['vidisp'] not in ('svf', 'lipanel', 'ss', 'lcpanel'):
            svp.vi_display = 0
            return
        if scene.frame_current not in range(svp['liparams']['fs'], svp['liparams']['fe'] + 1):
            self.disp_op.report({'INFO'},"Outside result frame range")
            return
        if svp.vi_display_rp != True \
             or (bpy.context.active_object not in self.obreslist and svp.vi_display_sel_only == True)  \
             or (bpy.context.active_object and bpy.context.active_object.mode == 'EDIT'):
             return
        
        if (self.width, self.height) != viewdesc(context)[2:]:
            mid_x, mid_y, self.width, self.height = viewdesc(context)
            self.u = 1
            
        if self.view_location != retvpvloc(context):
            self.view_location = retvpvloc(context)
            self.u = 1
            
        if svp.vi_display_sel_only == False:
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
#            bpy.context.user_preferences.system.window_draw_method = bpy.context.user_preferences.system.window_draw_method
           
    def update(self, context):
        scene = context.scene
        svp = scene.vi_params
        self.allpcs, self.alldepths, self.allres = array([]), array([]), array([])
#        try:
        for ob in self.obd:
            if ob.data.get('shape_keys') and str(self.fn) in [sk.name for sk in ob.data.shape_keys.key_blocks] and ob.active_shape_key.name != str(self.fn):
                ob.active_shape_key_index = [sk.name for sk in ob.data.shape_keys.key_blocks].index(str(self.fn))
#                try:
            omw = ob.matrix_world
            bm = bmesh.new()
            tempmesh = ob.to_mesh()
            bm.from_mesh(tempmesh)
            ob.to_mesh_clear()
#                    bpy.data.meshes.remove(tempmesh)
            bm.transform(omw)
            bm.normal_update()
            geom = bm.faces if bm.faces.layers.float.get('res{}'.format(scene.frame_current)) else bm.verts
            geom.ensure_lookup_table()
            livires = geom.layers.float['res{}'.format(scene.frame_current)]
    
            if bm.faces.layers.float.get('res{}'.format(scene.frame_current)):
                if svp.vi_disp_3d:
                    extrude = geom.layers.int['extrude']                                
                    faces = [f for f in geom if f.select and f[extrude]]
                else:
                    faces = [f for f in geom if f.select]

                distances = [(self.view_location - f.calc_center_median_weighted() + svp.vi_display_rp_off * f.normal.normalized()).length for f in faces]
   
                if svp.vi_display_vis_only:
                    fcos = [f.calc_center_median_weighted() + svp.vi_display_rp_off * f.normal.normalized() for f in faces]
                    direcs = [self.view_location - f for f in fcos]
                    (faces, distances) = map(list, zip(*[[f, distances[i]] for i, f in enumerate(faces) if not scene.ray_cast(context.view_layer, fcos[i], direcs[i], distance=distances[i])[0]]))

                face2d = [view3d_utils.location_3d_to_region_2d(context.region, context.region_data, f.calc_center_median_weighted()) for f in faces]
                (faces, pcs, depths) = map(list, zip(*[[f, face2d[fi], distances[fi]] for fi, f in enumerate(faces) if face2d[fi] and 0 < face2d[fi][0] < self.width and 0 < face2d[fi][1] < self.height]))          
                res = [f[livires] for f in faces] 
                res = ret_res_vals(svp, res)
                
            elif bm.verts.layers.float.get('res{}'.format(scene.frame_current)):                        
                verts = [v for v in geom if not v.hide and v.select and (context.space_data.region_3d.view_location - self.view_location).dot(v.co + svp.vi_display_rp_off * v.normal.normalized() - self.view_location)/((context.space_data.region_3d.view_location-self.view_location).length * (v.co + svp.vi_display_rp_off * v.normal.normalized() - self.view_location).length) > 0]
                distances = [(self.view_location - v.co + svp.vi_display_rp_off * v.normal.normalized()).length for v in verts]

                if svp.vi_display_vis_only:
                    vcos = [v.co + svp.vi_display_rp_off * v.normal.normalized() for v in verts]
                    direcs = [self.view_location - v for v in vcos]
                    (verts, distances) = map(list, zip(*[[v, distances[i]] for i, v in enumerate(verts) if not scene.ray_cast(vcos[i], direcs[i], distance=distances[i])[0]]))
                    
                vert2d = [view3d_utils.location_3d_to_region_2d(context.region, context.region_data, v.co) for v in verts]
                (verts, pcs, depths) = map(list, zip(*[[v, vert2d[vi], distances[vi]] for vi, v in enumerate(verts) if vert2d[vi] and 0 < vert2d[vi][0] < self.width and 0 < vert2d[vi][1] < self.height]))
                res = [v[livires] for v in verts] if not svp.vi_res_mod else [eval('{}{}'.format(v[livires], svp.vi_res_mod)) for v in verts]
                
            bm.free()
            self.allpcs = nappend(self.allpcs, array(pcs))
            self.alldepths = nappend(self.alldepths, array(depths))
            self.allres = nappend(self.allres, array(res))

        self.alldepths = self.alldepths/nmin(self.alldepths)        
        draw_index_distance(self.allpcs, self.allres, self.fontmult * svp.vi_display_rp_fs, svp.vi_display_rp_fc, svp.vi_display_rp_fsh, self.alldepths)    

#class Base_Display():
#    def __init__(self, pos, width, height, iname, xdiff, ydiff):
#        self.pos = pos
#        self.ispos = pos
#        self.iepos = [pos[0] + 40, pos[1] + 40]
#        self.spos = [int(self.pos[0] - 0.025 * width), int(self.pos[1] - 0.0125 * height)]
#        self.epos = [int(self.pos[0] + 0.025 * width), int(self.pos[1] + 0.0125 * height)]  
#        self.lspos = [self.spos[0], self.spos[1] - ydiff]
#        self.lepos = [self.spos[0] + xdiff, self.spos[1]]
#        self.lpos = (self.pos[0] + 0.2 * width, self.pos[1] - 0.2 * height)
#        self.resize = 0
#        self.press = 0
#        self.move = 0
#        self.expand = 0
#        if iname not in bpy.data.images:
#            bpy.data.images.load(os.path.join(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))), 'Images', iname))
#        self.image = iname
#        self.hl = [1, 1, 1, 1]
#        self.cao = None
#        self.xdiff, self.ydiff = xdiff, ydiff
#        
#    def draw(self, context, width, height):
#        self.width, self.height = context.region.width, context.region.height
#        if self.pos[1] > height:
#            self.pos[1] = height
##        self.spos = (int(self.pos[0] - 25), int(self.pos[1] - 15))
##        self.epos = (int(self.pos[0] + 25), int(self.pos[1] + 15))
#        
#        if self.expand == 0:
#            self.drawclosed()
#            
#        if self.expand == 1:
#            self.drawopen(context)
#    
#    def drawclosed(self):
#        draw_icon_new(self) 
        
class Base_Display():
    def __init__(self, ipos, width, height, xdiff, ydiff):
        self.ispos = ipos
        self.iepos = [ipos[0] + 40, ipos[1] + 40]
        self.lspos = [self.ispos[0], self.ispos[1] - ydiff - 20]
        self.lepos = [self.ispos[0] + xdiff, self.ispos[1] - 20]
        self.resize, self.move, self.expand, self.xdiff, self.ydiff = 0, 0, 0, xdiff, ydiff
        self.hl = [1, 1, 1, 1]
        self.cao = None
        self.xdiff, self.ydiff = xdiff, ydiff
        self.ah = height
        self.aw = width
        
class results_bar():
    def __init__(self, images, pos, region):
        self.images = images
        self.pos = pos
        self.rh = 0
        self.rw = region.width
        self.shaders = [gpu.shader.from_builtin('2D_UNIFORM_COLOR'), gpu.shader.from_builtin('2D_UNIFORM_COLOR')]
#        self.height = 0
        self.f_indices = ((0, 1, 2), (2, 3, 0))
        self.tex_coords = ((0, 0), (1, 0), (1, 1), (0, 1))
        self.no = len(images)
        
        for im in images:
            if im not in bpy.data.images:
                bpy.data.images.load(os.path.join(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))), 'Images', im))
            self.shaders.append(gpu.shader.from_builtin('2D_IMAGE'))
        
    def draw(self, rh):
        v_coords = ((self.pos, rh - 85), (self.pos + self.no * 50, rh - 85),
                    (self.pos + self.no * 50, rh - 35), (self.pos, rh - 35))
        
        
        if self.rh != rh:
            self.batches = [batch_for_shader(self.shaders[1], 'TRIS', {"pos": v_coords}, indices = self.f_indices),
                            batch_for_shader(self.shaders[0], 'LINE_LOOP', {"pos": v_coords})]

            for i in range(self.no):
                pos = ((self.pos + 5 + i * 50, rh - 80), (self.pos + 45 + i * 50, rh - 80),(self.pos + 45 + i * 50, rh - 40), (self.pos + 5 + i * 50, rh - 40))
                self.batches.append(batch_for_shader(self.shaders[i + 2], 'TRI_FAN', {"pos": pos, "texCoord": self.tex_coords}))
                self.rh = rh
                
        for si, s in enumerate(self.shaders):  
            s.bind()
            if si == 0:
                s.uniform_float("color", (1, 1, 1, 1))
                self.batches[si].draw(s)
            elif si == 1:
                s.uniform_float("color", (0, 0, 0, 1))
                bgl.glEnable(bgl.GL_BLEND)
                bgl.glEnable(bgl.GL_LINE_SMOOTH)
                self.batches[si].draw(s)
                bgl.glDisable(bgl.GL_LINE_SMOOTH)
                bgl.glDisable(bgl.GL_BLEND)
            else:
                im = bpy.data.images[self.images[si - 2]]
                if im.gl_load():
                    raise Exception()
                bgl.glActiveTexture(bgl.GL_TEXTURE0)
                bgl.glBindTexture(bgl.GL_TEXTURE_2D, im.bindcode)
                s.uniform_int("image", 0)
                self.batches[si].draw(s)
    
def spnumdisplay(disp_op, context):
    scene = context.scene
    svp = scene.vi_params

    if bpy.data.objects.get('SPathMesh'):
        spob = bpy.data.objects['SPathMesh'] 
        ob_mat = spob.matrix_world
        mid_x, mid_y, width, height = viewdesc(context)
        vl = ret_vp_loc(context)
        blf_props(scene, width, height)
        
        if svp.sp_hd:
            pvecs = [ob_mat@mathutils.Vector(p[:]) for p in spob['numpos'].values()]
            pvals = [int(p.split('-')[1]) for p in spob['numpos'].keys()]
            p2ds = [view3d_utils.location_3d_to_region_2d(context.region, context.region_data, p) for p in pvecs]
            vispoints = [pi for pi, p in enumerate(pvals) if p2ds[pi] and 0 < p2ds[pi][0] < width and 0 < p2ds[pi][1] < height and scene.ray_cast(context.view_layer, vl, pvecs[pi] - vl, distance = (pvecs[pi] - vl).length)[4] == spob]
            
            if vispoints:
                hs = [pvals[pi] for pi in vispoints]
                posis = [p2ds[pi] for pi in vispoints]                
                draw_index(posis, hs, svp.vi_display_rp_fs, svp.vi_display_rp_fc, svp.vi_display_rp_fsh)
                
        if [ob.get('VIType') == 'Sun' for ob in bpy.data.objects] and svp['spparams']['suns'] == '0':
            sobs = [ob for ob in bpy.data.objects if ob.get('VIType') == 'Sun']
            
            if sobs and svp.sp_td:
                sunloc = ob_mat@sobs[0].location
                solpos = view3d_utils.location_3d_to_region_2d(context.region, context.region_data, sunloc)
                
                try:
                    if 0 < solpos[0] < width and 0 < solpos[1] < height and not scene.ray_cast(context.view_layer, sobs[0].location + 0.05 * (vl - sunloc), vl - sunloc)[0]:
                        soltime = datetime.datetime.fromordinal(scene.solday)
                        soltime += datetime.timedelta(hours = scene.sp_sh)
                        sre = sobs[0].rotation_euler
                        blf_props(scene, width, height)
                        draw_time(solpos, soltime.strftime('  %d %b %X') + ' alt: {:.1f} azi: {:.1f}'.format(90 - sre[0]*180/pi, (180, -180)[sre[2] < -pi] - sre[2]*180/pi), 
                                   svp.vi_display_rp_fs, svp.vi_display_rp_fc, svp.vi_display_rp_fsh)
                        
                except Exception as e:
                    print(e)
        blf.disable(0, 4)
    else:
        return
    
class bsdf_disp(Base_Display):
    def __init__(self, context, unit, pos, width, height, xdiff, ydiff):
        Base_Display.__init__(self, pos, width, height, xdiff, ydiff)
        self.plt = plt
        self.pw, self.ph = 0.175 * xdiff, 0.35 * ydiff
        self.radii = array([0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1, 1.125])
#        self.patch_select = 0
        self.type_select = 0
        self.patch_hl = 0
        self.scale_select = 'Log'
        self.buttons = {}
        self.num_disp = 0
        self.leg_max, self.leg_min = 100, 0
        self.base_unit = unit
        self.font_id = blf.load(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'Fonts/NotoSans-Regular.ttf'))
        self.dpi = 157      
        self.v_coords = [(0, 0), (0, 1), (1, 1), (1, 0)]
        self.f_indices = [(0, 1, 2), (2, 3, 0)]
        self.segments = (1, 8, 16, 20, 24, 24, 24, 16, 12)            
        self.radii = (10, 30, 50, 70, 90, 110, 130, 150, 170)
        self.f_colours = [(1, 1, 1, 1)] * (721 + 8 * 720)
        self.imspos = (self.lepos[0] - self.ydiff - 75, self.lspos[1])
        self.image = 'bsdfplot.png'
        self.vi_coords = [(0.02, 0.02), (0.02, 0.98), (0.98, 0.98), (0.98, 0.02)] 
        self.tex_coords = ((0, 0), (0, 1), (1, 1), (1, 0))
        self.sr = 0
        self.cr = 0
        self.sseg = 1
        self.cseg = 0
        self.srs = 1
        self.crs = 0
#        self.update(context)
        self.create_batch('all')
        self.get_data(context)
#        self.bsdf = context.active_object.active_material['bsdf']['xml']
    
    def get_data(self, context):
        self.mat = context.active_object.active_material
        bsdf = parseString(self.mat.vi_params['bsdf']['xml'])
#        coltype = [path.firstChild.data for path in bsdf.getElementsByTagName('ColumnAngleBasis')]
#        rowtype = [path.firstChild.data for path in bsdf.getElementsByTagName('RowAngleBasis')]
        self.radtype = [path.firstChild.data for path in bsdf.getElementsByTagName('Wavelength')]
        self.rad_select = self.radtype[0]
        self.dattype = [path.firstChild.data for path in bsdf.getElementsByTagName('WavelengthDataDirection')]
        self.direc = self.dattype[0]
        self.type_select = self.dattype[0].split()[0]
        self.dir_select = self.dattype[0].split()[1]
#        lthetas = [path.firstChild.data for path in bsdf.getElementsByTagName('LowerTheta')]
        self.uthetas = [float(path.firstChild.data) for path in bsdf.getElementsByTagName('UpperTheta')]
        self.phis = [int(path.firstChild.data) for path in bsdf.getElementsByTagName('nPhis')]
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
        self.fig = self.plt.figure(figsize=(7.5, 7.5), dpi = 100)
        ax = self.plt.subplot(111, projection = 'polar')
        ax.bar(0, 0)
        self.plt.title('{} {}'.format(self.rad_select, svp.vi_bsdf_direc), size = 19, y = 1.025)
        ax.axis([0, 2 * pi, 0, 1])
        ax.spines['polar'].set_visible(False)
        ax.xaxis.set_ticks([])
        ax.yaxis.set_ticks([])
        
        for dti, dt in enumerate(self.dattype):
            if dt == svp.vi_bsdf_direc:
                self.scat_select = dti
                break 
            
#        for rdi, raddat in enumerate(['{0[0]} {0[1]}'.format(z) for z in zip(self.radtype, self.dattype)]):
#            if raddat == '{} {} {}'.format(self.rad_select, self.type_select, self.dir_select):
#                self.scat_select = rdi
#                break 

        selectdat = self.scatdat[self.scat_select].reshape(145, 145)# if self.scale_select == 'Linear' else nlog10((self.scatdat[self.scat_select] + 1).reshape(145, 145)) 
        widths = [0] + [self.uthetas[w]/90 for w in range(9)]
        patches, p = [], 0
        sa = repeat(kfsa, self.phis)
        act = repeat(kfact, self.phis)
        patchdat = selectdat[self.sseg - 1] * act * sa * 100
        bg = self.plt.Rectangle((0, 0), 2 * pi, 1, color=mcm.get_cmap(svp.vi_leg_col)((0, 0.01)[svp.vi_bsdfleg_scale == '1']), zorder = 0)      

        for ring in range(1, 10):
            angdiv = pi/self.phis[ring - 1]
            anglerange = range(self.phis[ring - 1], 0, -1)# if self.type_select == 'Transmission' else range(self.phis[ring - 1])
            ri = widths[ring] - widths[ring-1]

            for wedge in anglerange:
                phi1, phi2 = wedge * 2 * angdiv - angdiv, (wedge + 1) * 2 * angdiv - angdiv
                patches.append(Rectangle((phi1, widths[ring - 1]), phi2 - phi1, ri)) 
                if self.num_disp:
                    y = 0 if ring == 1 else 0.5 * (widths[ring] + widths[ring-1])
                    self.plt.text(0.5 * (phi1 + phi2), y, ('{:.1f}', '{:.0f}')[patchdat[p] >= 10].format(patchdat[p]), ha="center", va = 'center', family='sans-serif', size=10)
                p += 1
                
        pc = PatchCollection(patches, norm=mcolors.LogNorm(vmin=leg_min, vmax = svp.vi_bsdfleg_max), cmap=self.col, linewidths = [0] + 144*[0.5], edgecolors = ('black',)) if svp.vi_bsdfleg_scale == '1' else PatchCollection(patches, cmap=self.col, linewidths = [0] + 144*[0.5], edgecolors = ('black',))        
#        pc.set_linewidth(repeat(array([0, 0.5]), array([1, 144])))
        pc.set_array(patchdat)
        
        ax.add_collection(pc)
        ax.add_artist(bg)
        self.plt.colorbar(pc, fraction=0.04, pad=0.02, format = '%3g').set_label(label='Percentage of incoming flux (%)', size=18)
        pc.set_clim(vmin=leg_min + 0.01, vmax= svp.vi_bsdfleg_max)
        self.plt.tight_layout()
        self.save(scene)
                        
    def save(self, scene):
        self.plt.savefig(os.path.join(scene.vi_params['viparams']['newdir'], 'images', 'bsdfplot.png'), bbox_inches='tight')

        if 'bsdfplot.png' not in [i.name for i in bpy.data.images]:
            self.gimage = bpy.data.images.load(os.path.join(scene.vi_params['viparams']['newdir'], 'images', 'bsdfplot.png'))
        else:
            bpy.data.images['bsdfplot.png'].filepath = os.path.join(scene.vi_params['viparams']['newdir'], 'images', 'bsdfplot.png')
            bpy.data.images['bsdfplot.png'].reload()
            self.gimage = bpy.data.images['bsdfplot.png']
            
    def ret_coords(self): 
        if True:#self.type == 'Klems':
            vl_coords, f_indices = [], []
            v_coords = [(0, 0)]
            steps = 720/self.segments[self.sr]
            
            if self.sseg == 1:
                va_coords = [(0, 0)]
                va_coords += [(self.radii[self.sr] * cos((x*0.5 - 360/(2 * self.segments[self.sr]))*pi/180), 
                              self.radii[self.sr] * sin((x*0.5 - 360/(2 * self.segments[self.sr]))*pi/180)) for x in range(720)]
                fa_indices = [(0, x + (1, -719)[x > 0 and not x%720], x) for x in range(1, 720 + 1)]
            else:
#                va_coords = [(self.radii[cr - 1] * cos((x*0.5 - 360/(2 * self.segments[cr]))*pi/180), 
#                              self.radii[cr - 1] * sin((x*0.5 - 360/(2 * self.segments[cr]))*pi/180)) for x in range(int(rs/2 * self.segments[cr]*720), int((rs)/2 * self.segments[cr]*720) + 1)]
                va_coords = [(-self.radii[self.sr - 1] * cos((x*0.5 - 360/(2 * self.segments[self.sr]))*pi/180), 
                              self.radii[self.sr - 1] * sin((x*0.5 - 360/(2 * self.segments[self.sr]))*pi/180)) for x in range(int((self.srs - 1) * 720/self.segments[self.sr]), int((self.srs) * 720/self.segments[self.sr]) + 1)]
                va_coords += [(-self.radii[self.sr] * cos((x*0.5 - 360/(2 * self.segments[self.sr]))*pi/180), 
                              self.radii[self.sr] * sin((x*0.5 - 360/(2 * self.segments[self.sr]))*pi/180)) for x in range(int((self.srs - 1) * 720/self.segments[self.sr]), int((self.srs)*720/self.segments[self.sr]) + 1)]

                fa_indices = [(x, x+1, x+720/self.segments[self.sr] + 1) for x in range(int(720/self.segments[self.sr]))] + [(x-1, x, x-720/self.segments[self.sr] - 1) for x in range(int(720/self.segments[self.sr] + 2), int(1440/self.segments[self.sr]) + 2)]

            fa_colours = [(0.5, 0.5, 0.5, 1) for _ in range(len(va_coords))]
         
            for ri, radius in enumerate(self.radii):                                          
                if ri < 8:
                    for si, s in enumerate(range(self.segments[ri + 1])):
#                        f_colours += [(1, 1, 1, 1)] * int(720/segments[ri + 1])
                        vl_coords += [(radius * cos((360/self.segments[ri + 1] * si + 360/(2 * self.segments[ri + 1])) * pi/180), 
                                      radius * sin((360/self.segments[ri + 1] * si + 360/(2 * self.segments[ri + 1])) * pi/180)), 
                                      (self.radii[ri +1] * cos((360/self.segments[ri + 1] * si + 360/(2 * self.segments[ri + 1])) * pi/180), 
                                       self.radii[ri+1] * sin((360/self.segments[ri + 1] * si + 360/(2 * self.segments[ri + 1])) * pi/180))]

                vl_coords1 = [[radius * cos((x*0.5 - 360/(2 * self.segments[ri]))*pi/180), radius * sin((x*0.5 - 360/(2 * self.segments[ri]))*pi/180)] for x in range(720)]
                vl_coords2 = [[radius * cos(((x*0.5 + 0.5) - 360/(2 * self.segments[ri]))*pi/180), radius * sin(((x*0.5 + 0.5) - 360/(2 * self.segments[ri]))*pi/180)] for x in range(720)]
                vl_coordstot = list(zip(vl_coords1, vl_coords2))
                vl_coords += [item for sublist in vl_coordstot for item in sublist]
                v_coords += vl_coords1
                        
            f_indices = [(0, x + (1, -719)[x > 0 and not x%720], x) for x in range(1, 720*9 + 1)][::-1] 
            return (vl_coords, v_coords, f_indices, va_coords, fa_colours, fa_indices)
        
    def create_batch(self, sel):
        line_vertex_shader = '''
            uniform mat4 ModelViewProjectionMatrix;
            uniform vec2 spos;
            uniform vec2 size;
            in vec2 position;
            
            void main()
                {
                   float xpos = spos[0] + position[0] * size[0];
                   float ypos = spos[1] + position[1] * size[1]; 
                   gl_Position = ModelViewProjectionMatrix * vec4(int(xpos), int(ypos), -0.1f, 1.0f);
                }
            '''
            
        line_fragment_shader = '''
            uniform vec4 colour;
            out vec4 FragColour;
            
            void main()
                {
                    FragColour = colour;
                }
           
            '''  
        arc_vertex_shader = '''
            uniform mat4 ModelViewProjectionMatrix;
            uniform vec2 spos;
            uniform vec2 size;
            in vec2 position;
            in vec4 colour;
            flat out vec4 a_colour;
            
            void main()
                {
                   float xpos = spos[0] + position[0] * size[0];
                   float ypos = spos[1] + position[1] * size[1]; 
                   gl_Position = ModelViewProjectionMatrix * vec4(int(xpos), int(ypos), -0.1f, 1.0f);
                   a_colour = colour;
                }
            '''
            
        arc_fragment_shader = '''
            flat in vec4 a_colour;
            out vec4 FragColour;
            
            void main()
                {
                    FragColour = a_colour;
                }
           
            '''
            
        image_vertex_shader = '''
            uniform mat4 ModelViewProjectionMatrix;
            uniform vec2 spos;
            uniform vec2 size;
            in vec2 texCoord;
            in vec2 position;
            out vec2 texCoord_interp;
            
            void main()
            {
              float xpos = spos[0] + position[0] * size[0];
              float ypos = spos[1] + position[1] * size[1]; 
              gl_Position = ModelViewProjectionMatrix * vec4(int(xpos), int(ypos), 0.0f, 1.0f);
//              gl_Position.z = 1.0;
              texCoord_interp = texCoord;
            }
        '''
        
        image_fragment_shader = '''
            in vec2 texCoord_interp;
            out vec4 fragColor;
            
            uniform sampler2D image;
            
            void main()
            {
              fragColor = texture(image, texCoord_interp);
            }

        '''
        (vl_coords, v_coords, f_indices, va_coords, fa_colours, fa_indices) = self.ret_coords()
        
        if sel == 'all':
            b_coords = [(self.lspos[0], self.lspos[1]), (self.lspos[0], self.lepos[1]), (self.lepos[0], self.lepos[1]), (self.lepos[0], self.lspos[1])]
            b_indices = [(0, 1, 3), (1, 2, 3)]
            b_colours = [(1, 1, 1, 1), (1, 1, 1, 1), (1, 1, 1, 1), (1, 1, 1, 1)]
            self.back_shader =  gpu.types.GPUShader(arc_vertex_shader, arc_fragment_shader) 
            self.arcline_shader = gpu.types.GPUShader(line_vertex_shader, line_fragment_shader) 
            self.image_shader = gpu.types.GPUShader(image_vertex_shader, image_fragment_shader)
            self.image_batch = batch_for_shader(self.image_shader, 'TRI_FAN', {"position": self.vi_coords, "texCoord": self.tex_coords})
#        f_colours = [(random.random(), 0, 0, 1) for x in range(len(v_coords))]
#        f_colours = [item for item in f_colours for i in range(3)]
            
            self.back_batch = batch_for_shader(self.back_shader, 'TRIS', {"position": b_coords, "colour": b_colours}, indices = b_indices)
            self.back_shader.bind()
            self.back_shader.uniform_float("size", (1, 1))
            self.back_shader.uniform_float("spos", (0, 0))
            self.arcline_batch = batch_for_shader(self.arcline_shader, 'LINES', {"position": vl_coords})
            self.arcline_shader.bind()
            self.arcline_shader.uniform_float("size", (1, 1))
            self.arcline_shader.uniform_float("spos", (self.lspos[0] + 175, self.ah - 300))
            self.arcline_shader.uniform_float("colour", (0.6, 0.6, 0.6, 1))
        self.arc_shader = gpu.types.GPUShader(arc_vertex_shader, arc_fragment_shader)
        self.arc_batch = batch_for_shader(self.arc_shader, 'TRIS', {"position": va_coords, "colour": fa_colours}, indices = fa_indices)
#        self.col_batch = batch_for_shader(self.col_shader, 'TRIS', {"position": vl_coords[4:], "colour": self.colours}, indices = fl_indices)
    
    def draw(self, context):
        if self.expand:
            self.ah = context.region.height
            self.aw = context.region.width
            self.back_batch.draw(self.back_shader)
    #        bgl.glEnable(bgl.GL_FLAT)
            self.arc_shader.bind()
            self.arc_shader.uniform_float("size", (1,1))
            self.arc_shader.uniform_float("spos", (self.lspos[0] + 175, self.ah - 300))
    #        self.arc_shader.uniform_float("colour", (1, 1, 1, 1)) 
            self.arc_batch.draw(self.arc_shader)
    #        bgl.glDisable(bgl.GL_FLAT)
             
    ##        bgl.glShadeModel(bgl.GL_SMOOTH);
            bgl.glEnable(bgl.GL_DEPTH_TEST)
            bgl.glDepthFunc(bgl.GL_LESS)
            bgl.glDepthMask(bgl.GL_FALSE)
            bgl.glEnable(bgl.GL_BLEND)
            bgl.glEnable(bgl.GL_LINE_SMOOTH)
    #        bgl.glEnable(bgl.GL_BLEND)
    #        bgl.glDepthMask(False)
    #        bgl.glBlendFunc(bgl.GL_SRC_ALPHA, bgl.GL_ONE_MINUS_SRC_ALPHA);
    #        bgl.glHint(bgl.GL_LINE_SMOOTH_HINT, bgl.GL_DONT_CARE);
            bgl.glLineWidth(1)
    #        bgl.glEnable(bgl.GL_BLEND)
    #        bgl.glEnable(bgl.GL_LINE_SMOOTH)
    #        bgl.glEnable(bgl.GL_MULTISAMPLE)
    #        bgl.glHint(bgl.GL_LINE_SMOOTH_HINT, bgl.GL_NICEST)
            self.arcline_batch.draw(self.arcline_shader)
            bgl.glDepthMask(bgl.GL_TRUE)
    
            bgl.glDisable(bgl.GL_LINE_SMOOTH)
            bgl.glDisable(bgl.GL_BLEND)
            self.image_shader.bind()
            self.image_shader.uniform_float("size", (self.xdiff - self.ydiff, self.ydiff))
            self.image_shader.uniform_float("spos", self.imspos) 
            im = bpy.data.images[self.image]
            if im.gl_load():
                raise Exception()
            bgl.glActiveTexture(bgl.GL_TEXTURE0)
            bgl.glBindTexture(bgl.GL_TEXTURE_2D, im.bindcode)
            self.image_shader.uniform_int("image", 0)
            self.image_batch.draw(self.image_shader)
            blf.enable(0, 4)
            blf.enable(0, 8)
            blf.shadow(self.font_id, 5, 0.7, 0.7, 0.7, 1)
            blf.size(self.font_id, 10, 144)
            blf.position(self.font_id, self.lspos[0] + 175 - blf.dimensions(self.font_id, 'Incoming')[0]*0.5, self.lepos[1] - blf.dimensions(self.font_id, 'Incoming')[1], 0) 
            blf.color(self.font_id, 0, 0, 0, 1)      
            blf.draw(self.font_id, 'Incoming')
            blf.disable(0, 8)  
            blf.disable(0, 4)
        
class svf_legend(Base_Display):
    def __init__(self, context, unit, pos, width, height, xdiff, ydiff):
        Base_Display.__init__(self, pos, width, height, xdiff, ydiff)
        self.unit = unit
        self.font_id = blf.load(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'Fonts/NotoSans-Regular.ttf'))
        self.dpi = 300
        self.levels = 20        
        self.v_coords = [(0, 0), (0, 1), (1, 1), (1, 0)]
        self.f_indices = [(0, 1, 2), (2, 3, 0)]
        self.rh = context.region.height
        self.update(context)
        self.create_batch()        
                        
    def update(self, context):        
        scene = context.scene
        svp = scene.vi_params
        self.levels = svp.vi_leg_levels
        self.lh = 1/(self.levels + 1.25)
        self.cao = context.active_object        
        self.cols = retcols(mcm.get_cmap(svp.vi_leg_col), self.levels)
        (self.minres, self.maxres) = leg_min_max(svp)
        self.col, self.scale = svp.vi_leg_col, svp.vi_leg_scale
        resdiff = self.maxres - self.minres
        
        if not svp.get('liparams'):
            svp.vi_display = 0
            return

        resvals = [format(self.minres + i*(resdiff)/self.levels, '.0f') for i in range(self.levels + 1)] if self.scale == '0' else \
                        [format(self.minres + (1 - log10(i)/log10(self.levels + 1))*(resdiff), '.0f') for i in range(1, self.levels + 2)[::-1]]
        self.resvals = ['{0} - {1}'.format(resvals[i], resvals[i+1]) for i in range(self.levels)]
        self.colours = [item for item in [self.cols[i] for i in range(self.levels)] for i in range(4)]  
        print(self.colours)              
        blf.size(self.font_id, 12, self.dpi)        
        self.titxdimen = blf.dimensions(self.font_id, self.unit)[0]
        self.resxdimen = blf.dimensions(self.font_id, self.resvals[-1])[0]
        self.mydimen = blf.dimensions(self.font_id, self.unit)[1]

    def ret_coords(self):      
        lh = 1/(self.levels + 1.25) 
        vl_coords = self.v_coords[:]
        fl1_indices = [tuple(array((0, 1, 2)) + 4 * i) for i in range(self.levels)]
        fl2_indices = [tuple(array((2, 3, 0)) + 4 * i) for i in range(self.levels)]
        fl_indices = list(fl1_indices) + list(fl2_indices)
        
        for i in range(0, self.levels):
            vl_coords += [(0, i * lh), (0.35, i * lh), (0.35, (i + 1) * lh), (0, (i + 1) * lh)]
        return (vl_coords, fl_indices)
    
    def draw(self, context):
        self.rw = context.region.width
        svp = context.scene.vi_params
        self.cols = retcols(mcm.get_cmap(svp.vi_leg_col), self.levels)
        
        if self.rh != context.region.height:
            self.lepos[1] = context.region.height - (self.rh - self.lepos[1])
            self.lspos[1] = self.lepos[1] - self.ydiff
            self.rh = context.region.height
        
        if self.expand:   
            if self.resize:
                self.xdiff = self.lepos[0] - self.lspos[0]
                self.ydiff = self.lepos[1] - self.lspos[1]
            elif self.move:
                self.lspos[1] = self.lepos[1] - self.ydiff
                self.lepos[0] = self.lspos[0] + self.xdiff
            if self.lepos[1] > self.rh:
                self.lspos[1] = self.rh - self.ydiff 
                self.lepos[1] = self.rh
            if self.lepos[0] > self.rw:
                self.lspos[0] = self.rw - self.xdiff   
                self.lepos[0] = self.rw
                
            self.base_shader.bind()
            self.base_shader.uniform_float("size", (self.xdiff, self.ydiff))
            self.base_shader.uniform_float("spos", self.lspos)
            self.base_shader.uniform_float("colour", self.hl)      
            self.base_batch.draw(self.base_shader)  
            
            if self.levels != svp.vi_leg_levels or self.colours != [item for item in [self.cols[i] for i in range(self.levels)] for i in range(4)] or (self.minres, self.maxres) != leg_min_max(svp):
                self.update(context)
                (vl_coords, fl_indices) = self.ret_coords()
                self.line_batch = batch_for_shader(self.line_shader, 'LINE_LOOP', {"position": vl_coords})
                self.col_batch = batch_for_shader(self.col_shader, 'TRIS', {"position": vl_coords[4:], "colour": self.colours}, indices = fl_indices)
                               
            self.col_shader.bind()
            self.col_shader.uniform_float("size", (self.xdiff, self.ydiff))
            self.col_shader.uniform_float("spos", self.lspos)  
            self.col_batch.draw(self.col_shader)            
            self.line_shader.bind()
            self.line_shader.uniform_float("size", (self.xdiff, self.ydiff))
            self.line_shader.uniform_float("spos", self.lspos)
            self.line_shader.uniform_float("colour", (0, 0, 0, 1))      
            self.line_batch.draw(self.line_shader)
            fontscale = max(self.titxdimen/(self.xdiff * 0.99), self.resxdimen/(self.xdiff * 0.65), self.mydimen * 1.1/(self.lh * self.ydiff))
            blf.enable(0, 4)
            blf.enable(0, 8)
            blf.shadow(self.font_id, 5, 0.7, 0.7, 0.7, 1)
            blf.size(self.font_id, 12, int(self.dpi/fontscale))
            blf.position(self.font_id, self.lspos[0] + (self.xdiff - blf.dimensions(self.font_id, self.unit)[0]) * 0.45, self.lepos[1] - 0.5 * (self.lh * self.ydiff) - blf.dimensions(self.font_id, self.unit)[1] * 0.3, 0) 
            blf.color(self.font_id, 0, 0, 0, 1)      
            blf.draw(self.font_id, self.unit)
            blf.shadow(self.font_id, 5, 0.8, 0.8, 0.8, 1)    
            blf.size(self.font_id, 12, int((self.dpi - 25)/fontscale))
            
            for i in range(self.levels):
                num = self.resvals[i]            
                ndimen = blf.dimensions(self.font_id, "{}".format(num))
                blf.position(self.font_id, int(self.lepos[0] - self.xdiff * 0.05 - ndimen[0]), int(self.lspos[1] + i * self.lh * self.ydiff) + int((self.lh * self.ydiff - ndimen[1])*0.55), 0)
                blf.draw(self.font_id, "{}".format(self.resvals[i]))
                
            blf.disable(0, 8)  
            blf.disable(0, 4)
           
    def create_batch(self):
        base_vertex_shader = '''
            uniform mat4 ModelViewProjectionMatrix;
            uniform vec2 spos;
            uniform vec2 size;
            in vec2 position;
            
            void main()
                {
                   float xpos = spos[0] + position[0] * size[0];
                   float ypos = spos[1] + position[1] * size[1]; 
                   gl_Position = ModelViewProjectionMatrix * vec4(int(xpos), int(ypos), -0.1f, 1.0f);
                }
            '''
            
        base_fragment_shader = '''
            uniform vec4 colour;
            out vec4 FragColour;
            
            void main()
                {
                    FragColour = colour;
                }
           
            '''
            
        col_vertex_shader = '''
            uniform mat4 ModelViewProjectionMatrix;
            uniform vec2 spos;
            uniform vec2 size;
            in vec2 position;
            in vec4 colour;
            flat out vec4 f_colour;
            
            void main()
                {
                   float xpos = spos[0] + position[0] * size[0];
                   float ypos = spos[1] + position[1] * size[1]; 
                   gl_Position = ModelViewProjectionMatrix * vec4(int(xpos), int(ypos), -0.1f, 1.0f);
                   f_colour = colour;
                }
            '''
            
        col_fragment_shader = '''
            uniform vec4 colour;
            out vec4 FragColour;
            flat in vec4 f_colour;
            
            void main()
                {
                    FragColour = f_colour;
                }
           
            '''  
            
        self.base_shader = gpu.types.GPUShader(base_vertex_shader, base_fragment_shader) 
        self.line_shader = gpu.types.GPUShader(base_vertex_shader, base_fragment_shader) 
        self.col_shader = gpu.types.GPUShader(col_vertex_shader, col_fragment_shader)
        (vl_coords, fl_indices) = self.ret_coords()
        self.base_batch = batch_for_shader(self.base_shader, 'TRIS', {"position": self.v_coords}, indices = self.f_indices)
        self.line_batch = batch_for_shader(self.line_shader, 'LINE_LOOP', {"position": vl_coords})
        self.col_batch = batch_for_shader(self.col_shader, 'TRIS', {"position": vl_coords[4:], "colour": self.colours}, indices = fl_indices)

class wr_legend(Base_Display):
    def __init__(self, context, unit, pos, width, height, xdiff, ydiff):
        Base_Display.__init__(self, pos, width, height, xdiff, ydiff)
        self.unit = unit
        self.font_id = 0
        self.dpi = 300
        self.update(context)
        self.create_batch()
        self.line_shader.bind()
        self.line_shader.uniform_float("colour", (0, 0, 0, 1))  
                        
    def update(self, context):
        scene = context.scene
        svp = scene.vi_params
        simnode = bpy.data.node_groups[svp['viparams']['restree']].nodes[svp['viparams']['resnode']]        
        self.cao = context.active_object

        if self.cao and self.cao.get('VIType') and self.cao['VIType'] == 'Wind_Plane':            
            self.levels = self.cao['nbins']
            maxres = self.cao['maxres']
        else:
            self.levels = simnode['nbins']
            maxres = simnode['maxres']
        
        self.cols = retcols(mcm.get_cmap(svp.vi_leg_col), self.levels)
        self.colours = [item for item in [self.cols[i] for i in range(self.levels)] for i in range(4)] 
        
        if not svp.get('liparams'):
            svp.vi_display = 0
            return
        
#        self.cols = retcols(mcm.get_cmap(svp.vi_leg_col), self.levels)
        self.resvals = ['{0:.0f} - {1:.0f}'.format(2*i, 2*(i+1)) for i in range(simnode['nbins'])]
        self.resvals[-1] = self.resvals[-1][:-int(len('{:.0f}'.format(maxres)))] + "Inf"
#        self.colours = [item for item in [self.cols[i] for i in range(self.levels)] for i in range(4)]
        
        blf.size(self.font_id, 12, self.dpi)        
        self.titxdimen = blf.dimensions(self.font_id, self.unit)[0]
        self.resxdimen = blf.dimensions(self.font_id, self.resvals[-1])[0]
        self.mydimen = blf.dimensions(self.font_id, self.unit)[1]
        
    def draw(self, ah, aw):
        self.ah = ah
        self.aw = aw

        if self.expand:
            if self.resize:
                self.xdiff = self.lepos[0] - self.lspos[0]
                self.ydiff = self.lepos[1] - self.lspos[1]
            elif self.move:
                self.lspos[1] = self.lepos[1] - self.ydiff
                self.lepos[0] = self.lspos[0] + self.xdiff
            if self.lepos[1] > ah:
                self.lspos[1] = ah - self.ydiff 
                self.lepos[1] = ah
            if self.lepos[0] > aw:
                self.lspos[0] = aw - self.xdiff   
                self.lepos[0] = aw
                
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
            blf.enable(0, 4)
            blf.enable(0, 8)
            blf.shadow(self.font_id, 5, 0.7, 0.7, 0.7, 1)
            blf.size(self.font_id, 12, int(self.dpi/fontscale))
            blf.position(self.font_id, self.lspos[0] + (self.xdiff - blf.dimensions(self.font_id, self.unit)[0]) * 0.45, self.lepos[1] - 0.5 * (self.lh * self.ydiff) - blf.dimensions(self.font_id, self.unit)[1] * 0.3, 0) 
            blf.color(self.font_id, 0, 0, 0, 1)      
            blf.draw(self.font_id, self.unit)
            blf.shadow(self.font_id, 5, 0.8, 0.8, 0.8, 1)    
            blf.size(self.font_id, 12, int((self.dpi - 50)/fontscale))
            
            for i in range(self.levels):
                num = self.resvals[i]            
                ndimen = blf.dimensions(self.font_id, "{}".format(num))
                blf.position(self.font_id, int(self.lepos[0] - self.xdiff * 0.05 - ndimen[0]), int(self.lspos[1] + i * self.lh * self.ydiff) + int((self.lh * self.ydiff - ndimen[1])*0.55), 0)
                blf.draw(self.font_id, "{}".format(self.resvals[i]))
                
            blf.disable(0, 8)  
            blf.disable(0, 4)
                
    def create_batch(self):
        base_vertex_shader = '''
            uniform mat4 ModelViewProjectionMatrix;
            uniform vec2 spos;
            uniform vec2 size;
            in vec2 position;
            
            void main()
                {
                   float xpos = spos[0] + position[0] * size[0];
                   float ypos = spos[1] + position[1] * size[1]; 
                   gl_Position = ModelViewProjectionMatrix * vec4(int(xpos), int(ypos), 0.0f, 1.0f);
                }
            '''
            
        base_fragment_shader = '''
            uniform vec4 colour;
            out vec4 FragColour;
            
            void main()
                {
                    FragColour = colour;
                }
           
            '''
            
        col_vertex_shader = '''
            uniform mat4 ModelViewProjectionMatrix;
            uniform vec2 spos;
            uniform vec2 size;
            in vec2 position;
            in vec4 colour;
            flat out vec4 f_colour;
            
            void main()
                {
                   float xpos = spos[0] + position[0] * size[0];
                   float ypos = spos[1] + position[1] * size[1]; 
                   gl_Position = ModelViewProjectionMatrix * vec4(int(xpos), int(ypos), 0.0f, 1.0f);
                   f_colour = colour;
                }
            '''
            
        col_fragment_shader = '''
            uniform vec4 colour;
            out vec4 FragColour;
            flat in vec4 f_colour;
            
            void main()
                {
                    FragColour = f_colour;
                }
           
            '''  
            
        self.base_shader = gpu.types.GPUShader(base_vertex_shader, base_fragment_shader) 
        self.line_shader = gpu.types.GPUShader(base_vertex_shader, base_fragment_shader) 
        self.col_shader = gpu.types.GPUShader(col_vertex_shader, col_fragment_shader)
        v_coords = [(0, 0), (0, 1), (1, 1), (1, 0)]
        f_indices = [(0, 1, 2), (2, 3, 0)]        
        lh = 1/(self.levels + 1.25) 
        vl_coords = v_coords
        f_indices = [(0, 1, 2), (2, 3, 0)]
        fl1_indices = [tuple(array((0, 1, 2)) +4 * i) for i in range(self.levels)]
        fl2_indices = [tuple(array((2, 3, 0)) +4 * i) for i in range(self.levels)]
        fl_indices = list(fl1_indices) + list(fl2_indices)
    
        for i in range(0, self.levels):
            vl_coords += [(0, i * lh), (0.4, i * lh), (0.4, (i + 1) * lh), (0, (i + 1) * lh)]

        self.base_batch = batch_for_shader(self.base_shader, 'TRIS', {"position": v_coords}, indices = f_indices)
        self.line_batch = batch_for_shader(self.line_shader, 'LINE_LOOP', {"position": vl_coords})
        self.col_batch = batch_for_shader(self.col_shader, 'TRIS', {"position": vl_coords[4:], "colour": self.colours}, indices = fl_indices)
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
        
        if self.cao and self.cao.vi_params.get('ws'):
            self.rcarray = array(self.cao.vi_params['table']) 
        else:
            self.rcarray = array([['Invalid object']])
            
    def create_batch(self):
        base_vertex_shader = '''
            uniform mat4 ModelViewProjectionMatrix;
            uniform vec2 spos;
            uniform vec2 size;
            in vec2 position;
            
            void main()
                {
                   float xpos = spos[0] + position[0] * size[0];
                   float ypos = spos[1] + position[1] * size[1]; 
                   gl_Position = ModelViewProjectionMatrix * vec4(int(xpos), int(ypos), 0.0f, 1.0f);
                }
        '''
        
        base_fragment_shader = '''
            uniform vec4 colour;
            out vec4 FragColour;
            
            void main()
                {
                    FragColour = colour;
                }
           
            '''
        self.base_shader = gpu.types.GPUShader(base_vertex_shader, base_fragment_shader) 
        self.line_shader = gpu.types.GPUShader(base_vertex_shader, base_fragment_shader)
        v_coords = [(0, 0), (0, 1), (1, 1), (1, 0)]
        f_indices = [(0, 1, 2), (2, 3, 0)]                
        vl_coords = v_coords
        rno = len(self.rcarray)
        cno = len(self.rcarray[0])
        rh = 1/rno 
        blf.size(0, 24, 300)
        ctws = array([int(max([blf.dimensions(0, 'u{}'.format(e))[0] for e in entry])) for entry in self.rcarray.T])
        ctws = ctws/sum(ctws)
        ctws = [sum(ctws[:i]) for i in range(4)] + [1]
                
        for ci in range(cno):
            for ri in range(rno):            
                vl_coords += [(ctws[ci], ri * rh), (ctws[ci + 1], ri * rh), (ctws[ci + 1], (ri + 1) * rh)]#, (ci * rw, (ri + 1) * rh), (ci * rw, ri * rh)]
        vl_coords += [(0, 1)]
        self.base_batch = batch_for_shader(self.base_shader, 'TRIS', {"position": v_coords}, indices = f_indices)
        self.line_batch = batch_for_shader(self.line_shader, 'LINE_LOOP', {"position": vl_coords})
        
    def draw(self, ah, aw):
        self.ah = ah
        self.aw = aw
        
        if self.expand:
            if self.resize:
                self.xdiff = self.lepos[0] - self.lspos[0]
                self.ydiff = self.lepos[1] - self.lspos[1]
            elif self.move:
                self.lspos[1] = self.lepos[1] - self.ydiff
                self.lepos[0] = self.lspos[0] + self.xdiff
            if self.lepos[1] > ah:
                self.lspos[1] = ah - self.ydiff 
                self.lepos[1] = ah
            if self.lepos[0] > aw:
                self.lspos[0] = aw - self.xdiff   
                self.lepos[0] = aw
                
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
            fid = self.font_id
            blf.enable(0, 4)
            blf.enable(0, 8)
            blf.shadow(self.font_id, 5, 0.7, 0.7, 0.7, 1)
            blf.size(fid, 24, 300)
            rcshape = self.rcarray.shape
            [rowno, colno] = self.rcarray.shape            
#            colpos = [int(0.01 * self.xdiff)]
            ctws = array([int(max([blf.dimensions(fid, '{}'.format(e))[0] for e in entry])) for entry in self.rcarray.T])
            ctws = self.xdiff * ctws/sum(ctws) 
            ctws = [sum(ctws[:i]) for i in range(4)] + [self.xdiff]
            ctws = [sum(ctws[i:i+2])/2 for i in range(4)]            
            coltextwidths = array([int(max([blf.dimensions(fid, '{}'.format(e))[0] for e in entry]) + 0.05 * self.xdiff) for entry in self.rcarray.T])
            colscale = sum(coltextwidths)/(self.xdiff * 0.98)
#            colwidths = (coltextwidths/colscale).astype(int)
           
#            for cw in colwidths:
#                colpos.append(cw + colpos[-1])
        
            maxrowtextheight = max([max([blf.dimensions(fid, '{}'.format(e))[1] for e in entry if e])  for entry in self.rcarray.T])
            rowtextheight = maxrowtextheight + 0.1 * self.ydiff/rowno
            rowscale = (rowno * rowtextheight)/(self.ydiff - self.xdiff * 0.025)
            rowheight = int((self.ydiff - self.xdiff * 0.01)/rowno)
        #    rowoffset = 0.5 * maxrowtextheight
            rowtops = [int(self.lepos[1]  - self.xdiff * 0.005 - r * rowheight) for r in range(rowno)]
            rowbots = [int(self.lepos[1]  - self.xdiff * 0.005 - (r + 1) * rowheight) for r in range(rowno)]
            rowmids = [0.5 * (rowtops[r] + rowbots[r]) for r in range(rowno)]
            
            if abs(max(colscale, rowscale) - 1) > 0.05:
                self.fontdpi = int(280/max(colscale, rowscale))
           
            blf.size(fid, 24, self.fontdpi)
            blf.color(fid, 0, 0, 0, 1)       
            
            for r in range(rcshape[0]):
                for c in range(rcshape[1]):
                    if self.rcarray[r][c]:                
                        if c == 0:
                            blf.position(fid, self.lspos[0] + 0.01 * self.xdiff, int(rowmids[r] - 0.4 * blf.dimensions(fid, 'H')[1]), 0)#.format(self.rcarray[r][c]))[1])), 0)#int(self.lepos[1] - rowoffset - rowheight * (r + 0.5)), 0)
                        else:
                            blf.position(fid, self.lspos[0] + ctws[c] - int(blf.dimensions(fid, '{}'.format(self.rcarray[r][c]))[0] * 0.5), int(rowmids[r] - 0.5 * blf.dimensions(fid, 'H')[1]), 0)
                        blf.draw(fid, '{}'.format(self.rcarray[r][c]))

            blf.disable(0, 8)
            blf.disable(0, 4)

class ss_scatter(Base_Display):
    def __init__(self, context, pos, width, height, xdiff, ydiff):
        Base_Display.__init__(self, pos, width, height, xdiff, ydiff)
        self.image = 'ss_scatter.png'
        self.font_id = 0
        self.dpi = int(0.15 * ydiff)
        self.v_coords = [(0, 0), (0, 1), (1, 1), (1, 0)] 
        self.vi_coords = [(0.02, 0.02), (0.02, 0.98), (0.98, 0.98), (0.98, 0.02)] 
        self.f_indices = ((0, 1, 2), (2, 3, 0))
        self.tex_coords = ((0, 0), (0, 1), (1, 1), (1, 0))
        self.update(context)
        self.create_batch()
        self.line_shader.bind()
        self.line_shader.uniform_float("colour", (0, 0, 0, 1)) 
                
    def update(self, context):
        scene = context.scene
        svp = scene.vi_params
        self.cao = context.active_object
        self.col = svp.vi_scatt_col
        
        if self.cao and self.cao.get('ss'):
#            self.unit = svp.wind_type 
            zdata = array(self.cao['ss'])
            zmax = nmax(zdata) if svp.vi_scatt_max == '0' else svp.vi_scatt_max_val
            zmin = nmin(zdata) if svp.vi_scatt_min == '0' else svp.vi_scatt_min_val
            (title, cbtitle) = ('Wind Speed', 'Speed (m/s)') if svp.wind_type == '0' else ('Wind Direction', u'Direction (\u00B0)')
            self.plt = plt
            draw_dhscatter(self, scene, self.cao['days'], self.cao['hours'], zdata, title, 'Days', 'Hours', cbtitle, zmin, zmax)  
            save_plot(self, scene, self.image)
            
    def create_batch(self):
        base_vertex_shader = '''
            uniform mat4 ModelViewProjectionMatrix;
            uniform vec2 spos;
            uniform vec2 size;
            in vec2 position;
            
            void main()
                {
                   float xpos = spos[0] + position[0] * size[0];
                   float ypos = spos[1] + position[1] * size[1]; 
                   gl_Position = ModelViewProjectionMatrix * vec4(int(xpos), int(ypos), 0.0f, 1.0f);
                }
        '''
        
        base_fragment_shader = '''
            uniform vec4 colour;
            out vec4 FragColour;
            
            void main()
                {
                    FragColour = colour;
                }
           
            '''
        image_vertex_shader = '''
            uniform mat4 ModelViewProjectionMatrix;
            uniform vec2 spos;
            uniform vec2 size;
            in vec2 texCoord;
            in vec2 position;
            out vec2 texCoord_interp;
            
            void main()
            {
              float xpos = spos[0] + position[0] * size[0];
              float ypos = spos[1] + position[1] * size[1]; 
              gl_Position = ModelViewProjectionMatrix * vec4(int(xpos), int(ypos), 0.0f, 1.0f);
//              gl_Position.z = 1.0;
              texCoord_interp = texCoord;
            }
        '''
        
        image_fragment_shader = '''
            in vec2 texCoord_interp;
            out vec4 fragColor;
            
            uniform sampler2D image;
            
            void main()
            {
              fragColor = texture(image, texCoord_interp);
            }

        '''
        self.base_shader = gpu.types.GPUShader(base_vertex_shader, base_fragment_shader) 
        self.line_shader = gpu.types.GPUShader(base_vertex_shader, base_fragment_shader)
        self.image_shader = gpu.types.GPUShader(image_vertex_shader, image_fragment_shader)
        self.base_batch = batch_for_shader(self.base_shader, 'TRIS', {"position": self.v_coords}, indices = self.f_indices)
        self.line_batch = batch_for_shader(self.line_shader, 'LINE_LOOP', {"position": self.v_coords})
        self.image_batch = batch_for_shader(self.image_shader, 'TRI_FAN', {"position": self.vi_coords, "texCoord": self.tex_coords})
        
    def draw(self, ah, aw):
        self.ah = ah
        self.aw = aw
        
        if self.expand:
            if self.resize:
                self.xdiff = self.lepos[0] - self.lspos[0]
                self.ydiff = self.lepos[1] - self.lspos[1]
            elif self.move:
                self.lspos[1] = self.lepos[1] - self.ydiff
                self.lepos[0] = self.lspos[0] + self.xdiff
            if self.lepos[1] > ah:
                self.lspos[1] = ah - self.ydiff 
                self.lepos[1] = ah
            if self.lepos[0] > aw:
                self.lspos[0] = aw - self.xdiff   
                self.lepos[0] = aw
                
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
            self.image_shader.bind()
            self.image_shader.uniform_float("size", (self.xdiff, self.ydiff))
            self.image_shader.uniform_float("spos", self.lspos) 
            im = bpy.data.images[self.image]
            if im.gl_load():
                raise Exception()
            bgl.glActiveTexture(bgl.GL_TEXTURE0)
            bgl.glBindTexture(bgl.GL_TEXTURE_2D, im.bindcode)
            self.image_shader.uniform_int("image", 0)
            self.image_batch.draw(self.image_shader)
            
    def show_plot(self, context):
        show_plot(self, context)
            
class ss_legend(Base_Display):
    def __init__(self, context, unit, pos, width, height, xdiff, ydiff):
        Base_Display.__init__(self, pos, width, height, xdiff, ydiff)
        self.base_unit = unit
        self.font_id = 0
        self.dpi = 300
        self.levels = 20        
        self.v_coords = [(0, 0), (0, 1), (1, 1), (1, 0)]
        self.f_indices = [(0, 1, 2), (2, 3, 0)]
        self.update(context)
        self.create_batch()
                                
    def update(self, context):        
        scene = context.scene
        svp = scene.vi_params
        self.levels = svp.vi_leg_levels
        self.lh = 1/(self.levels + 1.25)
        self.cao = context.active_object        
        self.cols = retcols(mcm.get_cmap(svp.vi_leg_col), self.levels)
        (self.minres, self.maxres) = leg_min_max(svp)
        self.col, self.scale = svp.vi_leg_col, svp.vi_leg_scale

        for key, val in unit_dict.items():
            if val == svp.li_disp_basic:
                self.base_unit =  key

        self.unit = self.base_unit if not svp.vi_leg_unit else svp.vi_leg_unit
        self.cols = retcols(mcm.get_cmap(svp.vi_leg_col), self.levels)
        resdiff = self.maxres - self.minres
        
        if not svp.get('liparams'):
            svp.vi_display = 0
            return
        dplaces = retdp(self.maxres, 1)
        resvals = [format(self.minres + i*(resdiff)/self.levels, '.{}f'.format(dplaces)) for i in range(self.levels + 1)] if self.scale == '0' else \
                        [format(self.minres + (1 - log10(i)/log10(self.levels + 1))*(resdiff), '.{}f'.format(dplaces)) for i in range(1, self.levels + 2)[::-1]]

        self.resvals = ['{0} - {1}'.format(resvals[i], resvals[i+1]) for i in range(self.levels)]
        self.colours = [item for item in [self.cols[i] for i in range(self.levels)] for i in range(4)]                
        blf.size(self.font_id, 12, self.dpi)        
        self.titxdimen = blf.dimensions(self.font_id, self.unit)[0]
        self.resxdimen = blf.dimensions(self.font_id, self.resvals[-1])[0]
        self.mydimen = blf.dimensions(self.font_id, self.unit)[1]

    def ret_coords(self):      
        lh = 1/(self.levels + 1.25) 
        vl_coords = self.v_coords[:]
        fl1_indices = [tuple(array((0, 1, 2)) + 4 * i) for i in range(self.levels)]
        fl2_indices = [tuple(array((2, 3, 0)) + 4 * i) for i in range(self.levels)]
        fl_indices = list(fl1_indices) + list(fl2_indices)
        
        for i in range(0, self.levels):
            vl_coords += [(0, i * lh), (0.35, i * lh), (0.35, (i + 1) * lh), (0, (i + 1) * lh)]
        return (vl_coords, fl_indices)
    
    def draw(self, context):
        self.ah = context.area.height
        self.aw = context.area.width
        svp = context.scene.vi_params
        
        if self.expand:
            if self.resize:
                self.xdiff = self.lepos[0] - self.lspos[0]
                self.ydiff = self.lepos[1] - self.lspos[1]
            elif self.move:
                self.lspos[1] = self.lepos[1] - self.ydiff
                self.lepos[0] = self.lspos[0] + self.xdiff
            if self.lepos[1] > self.ah:
                self.lspos[1] = self.ah - self.ydiff 
                self.lepos[1] = self.ah
            if self.lepos[0] > self.aw:
                self.lspos[0] = self.aw - self.xdiff   
                self.lepos[0] = self.aw
                
            self.base_shader.bind()
            self.base_shader.uniform_float("size", (self.xdiff, self.ydiff))
            self.base_shader.uniform_float("spos", self.lspos)
            self.base_shader.uniform_float("colour", self.hl)      
            self.base_batch.draw(self.base_shader)  

            if self.levels != svp.vi_leg_levels or self.cols != retcols(mcm.get_cmap(svp.vi_leg_col), self.levels) or (self.minres, self.maxres) != leg_min_max(svp):
                self.update(context)
                (vl_coords, fl_indices) = self.ret_coords()
                self.line_batch = batch_for_shader(self.line_shader, 'LINE_LOOP', {"position": vl_coords})
                self.col_batch = batch_for_shader(self.col_shader, 'TRIS', {"position": vl_coords[4:], "colour": self.colours}, indices = fl_indices)
                
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
            blf.enable(0, 4)
            blf.enable(0, 8)
            blf.shadow(self.font_id, 5, 0.7, 0.7, 0.7, 1)
            blf.size(self.font_id, 12, int(self.dpi/fontscale))
            blf.position(self.font_id, self.lspos[0] + (self.xdiff - blf.dimensions(self.font_id, self.unit)[0]) * 0.45, self.lepos[1] - 0.6 * (self.lh * self.ydiff) - blf.dimensions(self.font_id, self.unit)[1] * 0.3, 0) 
            blf.color(self.font_id, 0, 0, 0, 1)   
            blf.draw(self.font_id, self.unit)
            blf.shadow(self.font_id, 5, 0.8, 0.8, 0.8, 1)    
            blf.size(self.font_id, 12, int((self.dpi - 50)/fontscale))
            
            for i in range(self.levels):
                num = self.resvals[i]            
                ndimen = blf.dimensions(self.font_id, "{}".format(num))
                blf.position(self.font_id, int(self.lepos[0] - self.xdiff * 0.025 - ndimen[0]), int(self.lspos[1] + i * self.lh * self.ydiff) + int((self.lh * self.ydiff - ndimen[1])*0.55), 0)
                blf.draw(self.font_id, "{}".format(self.resvals[i]))
                
            blf.disable(0, 8)  
            blf.disable(0, 4)
           
    def create_batch(self):
        base_vertex_shader = '''
            uniform mat4 ModelViewProjectionMatrix;
            uniform vec2 spos;
            uniform vec2 size;
            in vec2 position;
            
            void main()
                {
                   float xpos = spos[0] + position[0] * size[0];
                   float ypos = spos[1] + position[1] * size[1]; 
                   gl_Position = ModelViewProjectionMatrix * vec4(int(xpos), int(ypos), -0.1f, 1.0f);
                }
            '''
            
        base_fragment_shader = '''
            uniform vec4 colour;
            out vec4 FragColour;
            
            void main()
                {
                    FragColour = colour;
                }
           
            '''
            
        col_vertex_shader = '''
            uniform mat4 ModelViewProjectionMatrix;
            uniform vec2 spos;
            uniform vec2 size;
            in vec2 position;
            in vec4 colour;
            flat out vec4 f_colour;
            
            void main()
                {
                   float xpos = spos[0] + position[0] * size[0];
                   float ypos = spos[1] + position[1] * size[1]; 
                   gl_Position = ModelViewProjectionMatrix * vec4(int(xpos), int(ypos), -0.1f, 1.0f);
                   f_colour = colour;
                }
            '''
            
        col_fragment_shader = '''
            uniform vec4 colour;
            out vec4 FragColour;
            flat in vec4 f_colour;
            
            void main()
                {
                    FragColour = f_colour;
                }
           
            '''  
            
        self.base_shader = gpu.types.GPUShader(base_vertex_shader, base_fragment_shader) 
        self.line_shader = gpu.types.GPUShader(base_vertex_shader, base_fragment_shader) 
        self.col_shader = gpu.types.GPUShader(col_vertex_shader, col_fragment_shader)
        (vl_coords, fl_indices) = self.ret_coords()
        self.base_batch = batch_for_shader(self.base_shader, 'TRIS', {"position": self.v_coords}, indices = self.f_indices)
        self.line_batch = batch_for_shader(self.line_shader, 'LINE_LOOP', {"position": vl_coords})
        self.col_batch = batch_for_shader(self.col_shader, 'TRIS', {"position": vl_coords[4:], "colour": self.colours}, indices = fl_indices)

class livi_legend(Base_Display):
    def __init__(self, context, unit, pos, width, height, xdiff, ydiff):
        Base_Display.__init__(self, pos, width, height, xdiff, ydiff)
        self.base_unit = unit
        self.font_id = blf.load(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'Fonts/NotoSans-Regular.ttf'))
        self.dpi = 157
        self.levels = 20        
        self.v_coords = [(0, 0), (0, 1), (1, 1), (1, 0)]
        self.f_indices = [(0, 1, 2), (2, 3, 0)]
        self.update(context)
        self.create_batch()
                                
    def update(self, context):        
        scene = context.scene
        svp = scene.vi_params
        self.levels = svp.vi_leg_levels
        self.lh = 1/(self.levels + 1.25)
        self.cao = context.active_object        
        self.cols = retcols(mcm.get_cmap(svp.vi_leg_col), self.levels)
        (self.minres, self.maxres) = leg_min_max(svp)
        self.col, self.scale = svp.vi_leg_col, svp.vi_leg_scale

        for key, val in unit_dict.items():
            if val == svp.li_disp_basic:
                self.base_unit =  key
                
        self.unit = self.base_unit if not svp.vi_leg_unit else svp.vi_leg_unit
        self.cols = retcols(mcm.get_cmap(svp.vi_leg_col), self.levels)
        resdiff = self.maxres - self.minres
        
        if not svp.get('liparams'):
            svp.vi_display = 0
            return
        dplaces = retdp(self.maxres, 1)
        resvals = [format(self.minres + i*(resdiff)/self.levels, '.{}f'.format(dplaces)) for i in range(self.levels + 1)] if self.scale == '0' else \
                        [format(self.minres + (1 - log10(i)/log10(self.levels + 1))*(resdiff), '.{}f'.format(dplaces)) for i in range(1, self.levels + 2)[::-1]]

        self.resvals = ['{0} - {1}'.format(resvals[i], resvals[i+1]) for i in range(self.levels)]
        self.colours = [item for item in [self.cols[i] for i in range(self.levels)] for i in range(4)]                
        blf.size(self.font_id, 12, self.dpi)        
        self.titxdimen = blf.dimensions(self.font_id, self.unit)[0]
        self.resxdimen = blf.dimensions(self.font_id, self.resvals[-1])[0]
        self.mydimen = blf.dimensions(self.font_id, self.unit)[1]

    def ret_coords(self):      
        lh = 1/(self.levels + 1.25) 
        vl_coords = self.v_coords[:]
        fl1_indices = [tuple(array((0, 1, 2)) + 4 * i) for i in range(self.levels)]
        fl2_indices = [tuple(array((2, 3, 0)) + 4 * i) for i in range(self.levels)]
        fl_indices = list(fl1_indices) + list(fl2_indices)
        
        for i in range(0, self.levels):
            vl_coords += [(0, i * lh), (0.35, i * lh), (0.35, (i + 1) * lh), (0, (i + 1) * lh)]
        return (vl_coords, fl_indices)
    
    def draw(self, context):
        self.ah = context.region.height
        self.aw = context.region.width
#        print(self.ah, self.aw)
        svp = context.scene.vi_params
        
        if self.expand:
            if self.resize:
                self.xdiff = self.lepos[0] - self.lspos[0]
                self.ydiff = self.lepos[1] - self.lspos[1]
            elif self.move:
                self.lspos[1] = self.lepos[1] - self.ydiff
                self.lepos[0] = self.lspos[0] + self.xdiff
            if self.lepos[1] > self.ah:
                self.lspos[1] = self.ah - self.ydiff 
                self.lepos[1] = self.ah
            if self.lepos[0] > self.aw:
                self.lspos[0] = self.aw - self.xdiff   
                self.lepos[0] = self.aw
                
            self.base_shader.bind()
            self.base_shader.uniform_float("size", (self.xdiff, self.ydiff))
            self.base_shader.uniform_float("spos", self.lspos)
            self.base_shader.uniform_float("colour", self.hl)      
            self.base_batch.draw(self.base_shader)  
            self.unit = svp.vi_leg_unit if svp.vi_leg_unit else self.unit
            
            if self.levels != svp.vi_leg_levels or self.cols != retcols(mcm.get_cmap(svp.vi_leg_col), self.levels) or (self.minres, self.maxres) != leg_min_max(svp):
                self.update(context)
                (vl_coords, fl_indices) = self.ret_coords()
                self.line_batch = batch_for_shader(self.line_shader, 'LINE_LOOP', {"position": vl_coords})
                self.col_batch = batch_for_shader(self.col_shader, 'TRIS', {"position": vl_coords[4:], "colour": self.colours}, indices = fl_indices)
                               
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
            blf.enable(0, 4)
            blf.enable(0, 8)
            blf.shadow(self.font_id, 5, 0.7, 0.7, 0.7, 1)
            blf.shadow_offset(self.font_id, 1, 1)
            blf.size(self.font_id, int(14/fontscale), self.dpi)
            blf.position(self.font_id, self.lspos[0] + (self.xdiff - blf.dimensions(self.font_id, self.unit)[0]) * 0.45, self.lepos[1] - 0.6 * (self.lh * self.ydiff) - blf.dimensions(self.font_id, self.unit)[1] * 0.3, 0) 
            blf.color(self.font_id, 0, 0, 0, 1)   
            blf.draw(self.font_id, self.unit)
            blf.shadow(self.font_id, 5, 0.8, 0.8, 0.8, 1)    
            blf.size(self.font_id, int(11/fontscale), self.dpi)
            
            for i in range(self.levels):
                num = self.resvals[i]            
                ndimen = blf.dimensions(self.font_id, "{}".format(num))
                blf.position(self.font_id, int(self.lepos[0] - self.xdiff * 0.025 - ndimen[0]), int(self.lspos[1] + i * self.lh * self.ydiff) + int((self.lh * self.ydiff - ndimen[1])*0.55), 0)
                blf.draw(self.font_id, "{}".format(self.resvals[i]))
                
            blf.disable(0, 8)  
            blf.disable(0, 4)
           
    def create_batch(self):
        base_vertex_shader = '''
            uniform mat4 ModelViewProjectionMatrix;
            uniform vec2 spos;
            uniform vec2 size;
            in vec2 position;
            
            void main()
                {
                   float xpos = spos[0] + position[0] * size[0];
                   float ypos = spos[1] + position[1] * size[1]; 
                   gl_Position = ModelViewProjectionMatrix * vec4(int(xpos), int(ypos), -0.1f, 1.0f);
                }
            '''
            
        base_fragment_shader = '''
            uniform vec4 colour;
            out vec4 FragColour;
            
            void main()
                {
                    FragColour = colour;
                }
           
            '''
            
        col_vertex_shader = '''
            uniform mat4 ModelViewProjectionMatrix;
            uniform vec2 spos;
            uniform vec2 size;
            in vec2 position;
            in vec4 colour;
            flat out vec4 f_colour;
            
            void main()
                {
                   float xpos = spos[0] + position[0] * size[0];
                   float ypos = spos[1] + position[1] * size[1]; 
                   gl_Position = ModelViewProjectionMatrix * vec4(int(xpos), int(ypos), -0.1f, 1.0f);
                   f_colour = colour;
                }
            '''
            
        col_fragment_shader = '''
            uniform vec4 colour;
            out vec4 FragColour;
            flat in vec4 f_colour;
            
            void main()
                {
                    FragColour = f_colour;
                }
           
            '''  
            
        self.base_shader = gpu.types.GPUShader(base_vertex_shader, base_fragment_shader) 
        self.line_shader = gpu.types.GPUShader(base_vertex_shader, base_fragment_shader) 
        self.col_shader = gpu.types.GPUShader(col_vertex_shader, col_fragment_shader)
        (vl_coords, fl_indices) = self.ret_coords()
        self.base_batch = batch_for_shader(self.base_shader, 'TRIS', {"position": self.v_coords}, indices = self.f_indices)
        self.line_batch = batch_for_shader(self.line_shader, 'LINE_LOOP', {"position": vl_coords})
        self.col_batch = batch_for_shader(self.col_shader, 'TRIS', {"position": vl_coords[4:], "colour": self.colours}, indices = fl_indices)
        
def draw_legend(self, scene, unit):
    font_id = 0
    blf.enable(0, 4)
    blf.enable(0, 8)
    blf.shadow(font_id, 5, 0.7, 0.7, 0.7, 1)    
    levels = len(self.resvals)
    xdiff = self.lepos[0] - self.lspos[0]
    ydiff = self.lepos[1] - self.lspos[1]
    lh = ydiff/(levels + 1.25)   
    blf.size(font_id, 12, 300)
    titxdimen = blf.dimensions(font_id, unit)[0]
    resxdimen = blf.dimensions(font_id, self.resvals[-1])[0]
    mydimen = blf.dimensions(font_id, unit)[1]
    fontscale = max(titxdimen/(xdiff * 0.9), resxdimen/(xdiff * 0.6), mydimen * 1.25/lh)
    blf.size(font_id, 12, int(300/fontscale))

    if not self.resize:
        self.lspos = [self.spos[0], self.spos[1] - ydiff]
        self.lepos = [self.lspos[0] + xdiff, self.spos[1]]            
    else:
        self.lspos = [self.spos[0], self.lspos[1]]
        self.lepos = [self.lepos[0], self.spos[1]]
    
    bgl.glLineWidth(1)
    self.legl_shader = gpu.shader.from_builtin('2D_UNIFORM_COLOR')
    self.legf_shader = gpu.shader.from_builtin('2D_UNIFORM_COLOR')
    self.legfc_shader = gpu.shader.from_builtin('2D_FLAT_COLOR')
    colours = [item for item in [self.cols[i] for i in range(levels)] for i in range(4)]
    v_coords = [(self.lspos[0], self.lspos[1]), (self.lspos[0], self.lepos[1]), (self.lepos[0], self.lepos[1]), (self.lepos[0], self.lspos[1])]
    vl_coords = v_coords
    f_indices = [(0, 1, 2), (2, 3, 0)]
    fl1_indices = [tuple(array((0, 1, 2)) +4 * i) for i in range(levels)]
    fl2_indices = [tuple(array((2, 3, 0)) +4 * i) for i in range(levels)]
    fl_indices = list(fl1_indices) + list(fl2_indices)

    for i in range(0, levels):
        vl_coords += [(self.lspos[0], int(self.lspos[1] + i * lh)), (int(self.lspos[0] + xdiff * 0.4), int(self.lspos[1] + i * lh)), (int(self.lspos[0] + xdiff * 0.4), int(self.lspos[1] + (i + 1) * lh)), (self.lspos[0], int(self.lspos[1] + (i + 1) * lh))]

    self.legl_batch = batch_for_shader(self.legl_shader, 'LINE_LOOP', {"pos": vl_coords})
    self.legf_batch = batch_for_shader(self.legf_shader, 'TRIS', {"pos": v_coords}, indices = f_indices)
    self.legfc_batch = batch_for_shader(self.legfc_shader, 'TRIS', {"pos": vl_coords[4:], "color": colours}, indices = fl_indices)
    bgl.glEnable(bgl.GL_BLEND)
    self.legf_shader.bind()
    self.legf_shader.uniform_float("color", (self.hl))
    self.legf_batch.draw(self.legf_shader)
    bgl.glDisable(bgl.GL_BLEND)
    
    self.legfc_shader.bind()
    self.legfc_batch.draw(self.legfc_shader)
    
    self.legl_shader.bind()
    self.legl_shader.uniform_float("color", (0, 0, 0, 1))
    self.legl_batch.draw(self.legl_shader)

    blf.position(font_id, self.lspos[0] + (xdiff - blf.dimensions(font_id, unit)[0]) * 0.45, self.spos[1] - 0.5 * lh - blf.dimensions(font_id, unit)[1] * 0.3, 0) 
    blf.color(font_id, 0, 0, 0, 1)      
    blf.draw(font_id, unit)
#    blf.enable(0, blf.SHADOW)
#    blf.enable(0, blf.KERNING_DEFAULT)
#    blf.shadow(0, 5, 0, 0, 0, 0.7)
    
#    bgl.glColor4f(*scene.vi_display_rp_fc)

    blf.shadow(font_id, 5, 0.8, 0.8, 0.8, 1)
    
    blf.size(font_id, 12, int(250/fontscale))
    bgl.glDisable(bgl.GL_BLEND)
    
#    self.legl_shader = gpu.shader.from_builtin('2D_UNIFORM_COLOR')
   
    for i in range(levels):
        num = self.resvals[i]
#        rgba = self.cols[i]
        bgl.glHint(bgl.GL_LINE_SMOOTH_HINT, bgl.GL_NICEST)
        ndimen = blf.dimensions(font_id, "{}".format(num))
        blf.position(font_id, int(self.lepos[0] - xdiff * 0.05 - ndimen[0]), int(self.lspos[1] + i * lh) + int((lh - ndimen[1])*0.5), 0)
#        bgl.glColor4f(0, 0, 0, 1)
        blf.draw(font_id, "{}".format(self.resvals[i]))
    
    bgl.glLineWidth(1)
#    bgl.glColor4f(0, 0, 0, 1)
    blf.disable(0, 8)  
    blf.disable(0, 4)
    
def draw_icon_new(self):    
#    IMAGE_NAME = "Untitled"
    image = bpy.data.images[self.image]

    
    shader = gpu.shader.from_builtin('2D_IMAGE')
    batch = batch_for_shader(
        shader, 'TRI_FAN',
        {
            "pos": ((305, self.height - 80), (345, self.height - 80), (345, self.height - 40), (305, self.height - 40)),
            "texCoord": ((0, 0), (1, 0), (1, 1), (0, 1)),
        },
    )
    
    if image.gl_load():
        raise Exception()

    
    def draw():
        bgl.glActiveTexture(bgl.GL_TEXTURE0)
        bgl.glBindTexture(bgl.GL_TEXTURE_2D, image.bindcode)
    
        shader.bind()
        shader.uniform_int("image", 0)
        batch.draw(shader)

    draw()
#    bpy.types.SpaceView3D.draw_handler_add(draw, (), 'WINDOW', 'POST_PIXEL')
    
def draw_icon(self):
    drawpoly(self.spos[0], self.spos[1], self.epos[0], self.epos[1], *self.hl)        
    drawloop(self.spos[0], self.spos[1], self.epos[0], self.epos[1])
    bgl.glEnable(bgl.GL_BLEND)
    bpy.data.images[self.image].gl_load(bgl.GL_NEAREST, bgl.GL_NEAREST)
    bgl.glBindTexture(bgl.GL_TEXTURE_2D, bpy.data.images[self.image].bindcode[0])
    bgl.glTexParameteri(bgl.GL_TEXTURE_2D,
                            bgl.GL_TEXTURE_MAG_FILTER, bgl.GL_LINEAR)
    bgl.glTexParameteri(bgl.GL_TEXTURE_2D,
                            bgl.GL_TEXTURE_MIN_FILTER, bgl.GL_LINEAR)
    bgl.glEnable(bgl.GL_TEXTURE_2D)
    bgl.glColor4f(1, 1, 1, 1)
    bgl.glBegin(bgl.GL_QUADS)
    bgl.glTexCoord2i(0, 0)
    bgl.glVertex2f(self.spos[0] + 1, self.spos[1] + 1)
    bgl.glTexCoord2i(1, 0)
    bgl.glVertex2f(self.epos[0] - 1, self.spos[1] + 1)
    bgl.glTexCoord2i(1, 1)
    bgl.glVertex2f(self.epos[0] - 1, self.epos[1] - 1)
    bgl.glTexCoord2i(0, 1)
    bgl.glVertex2f(self.spos[0] + 1, self.epos[1] - 1)
    bgl.glEnd()
    bgl.glDisable(bgl.GL_TEXTURE_2D)
    bgl.glDisable(bgl.GL_BLEND)
    bgl.glFlush()
    
def draw_image(self, topgap):
    draw_icon(self)
    self.xdiff = self.lepos[0] - self.lspos[0]
    self.ydiff = self.lepos[1] - self.lspos[1]
    if not self.resize:
        self.lspos = [self.spos[0], self.spos[1] - self.ydiff]
        self.lepos = [self.lspos[0] + self.xdiff, self.spos[1]]            
    else:
        self.lspos = [self.spos[0], self.lspos[1]]
        self.lepos = [self.lepos[0], self.spos[1]]

    bpy.data.images[self.gimage].reload()
    drawpoly(self.lspos[0], self.lspos[1], self.lepos[0], self.lepos[1], 1, 1, 1, 1)        
    drawloop(self.lspos[0], self.lspos[1], self.lepos[0], self.lepos[1])
    bgl.glEnable(bgl.GL_BLEND)
    bpy.data.images[self.gimage].gl_load(bgl.GL_NEAREST, bgl.GL_NEAREST)
    bgl.glBindTexture(bgl.GL_TEXTURE_2D, bpy.data.images[self.gimage].bindcode[0])
    bgl.glTexParameteri(bgl.GL_TEXTURE_2D,
                            bgl.GL_TEXTURE_MAG_FILTER, bgl.GL_LINEAR)
    bgl.glTexParameteri(bgl.GL_TEXTURE_2D,
                            bgl.GL_TEXTURE_MIN_FILTER, bgl.GL_LINEAR)
    bgl.glEnable(bgl.GL_TEXTURE_2D)
    bgl.glColor4f(1, 1, 1, 1)
    bgl.glBegin(bgl.GL_QUADS)
    bgl.glTexCoord2i(0, 0)
    bgl.glVertex2f(self.lspos[0] + 5, self.lspos[1] + 5)
    bgl.glTexCoord2i(1, 0)
    bgl.glVertex2f(self.lepos[0] - 5, self.lspos[1] + 5)
    bgl.glTexCoord2i(1, 1)
    bgl.glVertex2f(self.lepos[0] - 5, self.lepos[1] - topgap)
    bgl.glTexCoord2i(0, 1)
    bgl.glVertex2f(self.lspos[0] + 5, self.lepos[1] - topgap)
    bgl.glEnd()
    bgl.glDisable(bgl.GL_TEXTURE_2D)
    bgl.glFlush()
    
#def draw_dhscatter(self, scene, x, y, z, tit, xlab, ylab, zlab, valmin, valmax):
#    self.plt.close()
#    self.col = scene.vi_leg_col
#    x = [x[0] - 0.5] + [xval + 0.5 for xval in x] 
#    y = [y[0] - 0.5] + [yval + 0.5 for yval in y]
#    self.plt.figure(figsize=(6 + len(x)/len(y), 6))
#    
#    self.plt.title(tit, size = 20).set_position([.5, 1.025])
#    self.plt.xlabel(xlab, size = 18)
#    self.plt.ylabel(ylab, size = 18)
#    self.plt.pcolor(x, y, z, cmap=self.col, vmin=valmin, vmax=valmax)#, norm=plt.matplotlib.colors.LogNorm())#, edgecolors='b', linewidths=1, vmin = 0, vmax = 4000)
#    cbar = self.plt.colorbar(use_gridspec=True)
#    cbar.set_label(label=zlab,size=18)
#    cbar.ax.tick_params(labelsize=16)
#    self.plt.axis([min(x),max(x),min(y),max(y)])
#    self.plt.xticks(size = 16)
#    self.plt.yticks(size = 16)
#    self.plt.tight_layout(rect=[0, 0, 1 + ((len(x)/len(y)) - 1) * 0.005, 1])

def draw_dhscatter(self, scene, x, y, z, tit, xlab, ylab, zlab, valmin, valmax):
    self.plt.close()
    x = [x[0] - 0.5] + [xval + 0.5 for xval in x] 
    y = [y[0] - 0.5] + [yval + 0.5 for yval in y]
    self.plt.figure(figsize=(4 + len(x)/len(y), 6))    
    self.plt.title(tit, size = 20).set_position([.5, 1.025])
    self.plt.xlabel(xlab, size = 18)
    self.plt.ylabel(ylab, size = 18)
    self.plt.pcolor(x, y, z, cmap=self.col, vmin=valmin, vmax=valmax)#, norm=plt.matplotlib.colors.LogNorm())#, edgecolors='b', linewidths=1, vmin = 0, vmax = 4000)
    cbar = self.plt.colorbar(use_gridspec=True)
    cbar.set_label(label=zlab,size=18)
    cbar.ax.tick_params(labelsize=16)
    self.plt.axis([min(x),max(x),min(y),max(y)])
    self.plt.xticks(size = 16)
    self.plt.yticks(size = 16)
    self.plt.tight_layout(rect=[0, 0, 1 + ((len(x)/len(y)) - 1) * 0.005, 1])  
    
def draw_table(self):
    draw_icon(self) 
    font_id = 0
    blf.enable(0, 4)
    blf.enable(0, 8)
    blf.shadow(font_id, 5, 0.9, 0.9, 0.9, 1)
    blf.size(font_id, 44, self.fontdpi)
    rcshape = self.rcarray.shape
    [rowno, colno] = self.rcarray.shape
    
    self.xdiff = self.lepos[0] - self.lspos[0]
    self.ydiff = self.lepos[1] - self.lspos[1]
    colpos = [int(0.01 * self.xdiff)]
    
    if not self.resize:
        self.lspos = [self.spos[0], self.spos[1] - self.ydiff]
        self.lepos = [self.lspos[0] + self.xdiff, self.spos[1]]            
    else:
        self.lspos = [self.spos[0], self.lspos[1]]
        self.lepos = [self.lepos[0], self.spos[1]]
        
    coltextwidths = array([int(max([blf.dimensions(font_id, '{}'.format(e))[0] for e in entry]) + 0.05 * self.xdiff) for entry in self.rcarray.T])
    colscale = sum(coltextwidths)/(self.xdiff * 0.98)
    colwidths = (coltextwidths/colscale).astype(int)
   
    for cw in colwidths:
        colpos.append(cw + colpos[-1])

    maxrowtextheight = max([max([blf.dimensions(font_id, '{}'.format(e))[1] for e in entry if e])  for entry in self.rcarray.T])
    rowtextheight = maxrowtextheight + 0.1 * self.ydiff/rowno
    rowscale = (rowno * rowtextheight)/(self.ydiff - self.xdiff * 0.025)
    rowheight = int((self.ydiff - self.xdiff * 0.01)/rowno)
#    rowoffset = 0.5 * maxrowtextheight
    rowtops = [int(self.lepos[1]  - self.xdiff * 0.005 - r * rowheight) for r in range(rowno)]
    rowbots = [int(self.lepos[1]  - self.xdiff * 0.005 - (r + 1) * rowheight) for r in range(rowno)]
    rowmids = [0.5 * (rowtops[r] + rowbots[r]) for r in range(rowno)]
    
    if abs(max(colscale, rowscale) - 1) > 0.05:
        self.fontdpi = int(self.fontdpi/max(colscale, rowscale))
   
    blf.size(font_id, 48, self.fontdpi)
    drawpoly(self.lspos[0], self.lspos[1], self.lepos[0], self.lepos[1], 1, 1, 1, 1)        
    drawloop(self.lspos[0], self.lspos[1], self.lepos[0], self.lepos[1])       
    bgl.glEnable(bgl.GL_BLEND)
    bgl.glBlendFunc(bgl.GL_SRC_ALPHA, bgl.GL_ONE_MINUS_SRC_ALPHA)
    
    for r in range(rcshape[0]):
        for c in range(rcshape[1]):
            if self.rcarray[r][c]:                
                if c == 0:
                    blf.position(font_id, self.lspos[0] + colpos[c] + 0.005 * self.xdiff, int(rowmids[r] - 0.5 * blf.dimensions(font_id, 'H')[1]), 0)#.format(self.rcarray[r][c]))[1])), 0)#int(self.lepos[1] - rowoffset - rowheight * (r + 0.5)), 0)
                else:
                    blf.position(font_id, self.lspos[0] + colpos[c] + colwidths[c] * 0.5 - int(blf.dimensions(font_id, '{}'.format(self.rcarray[r][c]))[0] * 0.5), int(rowmids[r] - 0.5 * blf.dimensions(font_id, 'H')[1]), 0)
                drawloop(int(self.lspos[0] + colpos[c]), rowtops[r], self.lspos[0] + colpos[c + 1], rowbots[r])                
                if self.rcarray[r][c] == 'Pass':
                    bgl.glColor3f(0.0, 0.6, 0.0)
                elif self.rcarray[r][c] == 'Fail':
                    bgl.glColor3f(0.6, 0.0, 0.0)
                else:
                    bgl.glColor3f(0.0, 0.0, 0.0)
                blf.draw(font_id, '{}'.format(self.rcarray[r][c]))
#    else:
#        for r in range(rcshape[0]):
#            for c in range(rcshape[1]):
#                if self.rcarray[r][c]:
#                    if c == 0:
#                        blf.position(font_id, self.lspos[0] + colpos[c] + 0.01 * self.xdiff, self.lepos[1] -  0.01 * self.xdiff - int(rowheight * (r + 0.25)) - int(blf.dimensions(font_id, '{}'.format(self.rcarray[1][1]))[1]), 0)
#                    else:
#                        blf.position(font_id, self.lspos[0] + colpos[c] + colwidths[c] * 0.5 - int(blf.dimensions(font_id, '{}'.format(self.rcarray[r][c]))[0] * 0.5), self.lepos[1] -  0.01 * self.xdiff - int(rowheight * (r + 0.25)) - int(blf.dimensions(font_id, '{}'.format(self.rcarray[1][1]))[1]), 0)
#                    drawloop(int(self.lspos[0] + colpos[c]), int(self.lepos[1] - 0.01 * self.xdiff - r * rowheight), self.lspos[0] + colpos[c + 1], int(self.lepos[1] - 0.01 * self.xdiff - (r + 1) * rowheight))                
#                    blf.draw(font_id, '{}'.format(self.rcarray[r][c]))
    bgl.glDisable(bgl.GL_BLEND) 
    blf.disable(0, 8)
    blf.disable(0, 4)
    bgl.glEnd()
    bgl.glFlush()

def save_plot(self, scene, filename):
    fileloc = os.path.join(scene.vi_params['viparams']['newdir'], 'images', filename)
    self.plt.savefig(fileloc, pad_inches = 0.1)
    
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
        
#class wr_legend(Base_Display):
#    def __init__(self, pos, width, height, iname, xdiff, ydiff):
#        Base_Display.__init__(self, pos, width, height, iname, xdiff, ydiff)
#        
#    def update(self, context):
#        scene = context.scene
#        simnode = bpy.data.node_groups[scene.vi_params['viparams']['restree']].nodes[scene.vi_params['viparams']['resnode']]        
#        self.cao = context.active_object
#        covp = self.cao.vi_params
#
#        if self.cao and covp.get('VIType') and covp['VIType'] == 'Wind_Plane':            
#            levels = covp['nbins']
#            maxres = covp['maxres']
#        else:
#            levels = simnode['nbins']
#            maxres = simnode['maxres']
#        self.cols = retcols(mcm.get_cmap(scene.vi_leg_col), levels)
#        
#        if not scene.get('liparams'):
#            scene.vi_display = 0
#            return
#
#        self.resvals = ['{0:.0f} - {1:.0f}'.format(2*i, 2*(i+1)) for i in range(simnode['nbins'])]
#        self.resvals[-1] = self.resvals[-1][:-int(len('{:.0f}'.format(maxres)))] + "Inf"  
#        
#    def drawopen(self, context):
#        draw_legend(self, context.scene, 'Speed (m/s)')
        
#class wr_scatter(Base_Display):
#    def __init__(self, pos, width, height, iname, xdiff, ydiff):
#        Base_Display.__init__(self, pos, width, height, iname, xdiff, ydiff)
#        self.unit = '0'
#        
#    def update(self, context):
#        self.cao = context.active_object
#        covp = self.cao.vi_params
#        if self.cao and covp.get('ws'):
#            self.unit = context.scene.wind_type 
#            zdata = array(covp['ws']) if context.scene.wind_type == '0' else array(covp['wd'])
#            (title, cbtitle) = ('Wind Speed', 'Speed (m/s)') if context.scene.wind_type == '0' else ('Wind Direction', u'Direction (\u00B0)')
#            self.plt = plt
#            draw_dhscatter(self, context.scene, covp['days'], covp['hours'], zdata, title, 'Days', 'Hours', cbtitle, nmin(zdata), nmax(zdata))  
#            save_plot(self, context.scene, 'scatter.png')
#        
#    def drawopen(self, context):
#        draw_image(self, 0)
#        
#    def show_plot(self):
#        show_plot(self)

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
        self.sumnorm = mathutils.Matrix().Rotation(pi/2, 4, 'X')@mathutils.Vector([0] + list((sumcoords[0]-sumcoords[1])[1:])).normalized()   
        self.winnorm = mathutils.Matrix().Rotation(pi/2, 4, 'X')@mathutils.Vector([0] + list((wincoords[0]-wincoords[1])[1:])).normalized()    
                
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
#        sun2_vertex_shader = '''
#            uniform mat4 viewProjectionMatrix;
#            uniform mat4 sp_matrix;
#            in vec3 position;
#            in vec3 sun_position;
#//            uniform float sun_radius;
#//            out vec4 sun_position;
#//            out mat4 sp_matrix;
#//            out mat4 viewProjectionMatrix;
#            out float sun_alpha;
#            out vec4 sun_Position;
#            out vec3 sp;
#            
#            void main()
#                {
#                    gl_Position = viewProjectionMatrix * sp_matrix * vec4(position, 1.0f);
#                    sun_Position = viewProjectionMatrix * sp_matrix * vec4(sun_position, 1.0f); 
#//                    sun_Position = ftransform(vec4(sun_position, 1.0f));
#//                    sun_alpha = length(gl_Position.xy - sun_Position.xy) * 0.1;
#//                    sp = sun_position;
#                }
#            '''
#        sun2_fragment_shader = '''
#            uniform vec4 sun_colour;
#            uniform vec4 viewport;
#            out vec4 FragColour;
#            in vec4 sun_Position;
#            
#            void main()
#                {
#                    vec3 ndc = sun_Position.xyz/sun_Position.w;
#                    vec2 vp_coord = ndc.xy * 0.5 + 0.5;
#                    vec2 vpp_coord = vp_coord * viewport.zw;
#                    float pos = length(gl_FragCoord.xy - vpp_coord);
#                    float radius = sun_Position.z * 100;
#                    if (pos > radius) 
#                        {discard;}
#                    if (pos >= radius - 10) 
#                        {
#                            FragColour = sun_colour;
#                            FragColour[3] = 0.1*(radius - 10 - pos);
#                        }
#                    if (pos < radius -10) 
#                        {FragColour = sun_colour;}
#                }
#           
#            '''
            
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
        svp = scene.vi_params
       
        if context.area:
            context.area.tag_redraw()
            
        if svp.vi_display == 0 or svp['viparams']['vidisp'] != 'sp':
            bpy.types.SpaceView3D.draw_handler_remove(self.draw_handle_sp, "WINDOW")
            bpy.types.SpaceView3D.draw_handler_remove(self.draw_handle_spnum, 'WINDOW')
            svp['viparams']['vidisp'] = ''
            
            for h in bpy.app.handlers.frame_change_post:
                bpy.app.handlers.frame_change_post.remove(h)
                
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
        node = context.node
        scene.display.shadow_focus = 1
#        context.space_data.shading.light = 'FLAT'
        svp = scene.vi_params
        svp['viparams'] = {}
        svp['spparams'] = {}
        svp['spparams']['suns'] = node.suns
        spcoll = create_coll(context, 'SunPath')            
        sd = 100
        
        # Set the node colour
        node.export()
        
        svp['viparams']['resnode'], svp['viparams']['restree'] = node.name, node.id_data.name
        scene.cursor.location = (0.0, 0.0, 0.0)
        suns = [ob for ob in scene.objects if ob.parent and ob.type == 'LIGHT' and ob.data.type == 'SUN' and ob.parent.get('VIType') == "SPathMesh" ]
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
            sun['solhour'], sun['solday'] = svp.sp_sh, svp.sp_sd
            sun.name = sun.data.name ='Sun{}'.format(s)
            sun.parent = spathob

        if spfc not in bpy.app.handlers.frame_change_post:
            bpy.app.handlers.frame_change_post.append(spfc)

        svp['viparams']['vidisp'] = 'sp'
        svp['viparams']['visimcontext'] = 'SunPath'
        sunpath(scene)
#        node = context.node
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
        context.window_manager.modal_handler_add(self)
        svp.vi_display = 1
        return {'RUNNING_MODAL'}
        
class wr_scatter(Base_Display):
    def __init__(self, context, pos, width, height, xdiff, ydiff):
        Base_Display.__init__(self, pos, width, height, xdiff, ydiff)
        self.image = 'wind_scatter.png'
        self.font_id = 0
        self.dpi = int(0.15 * ydiff)
        self.v_coords = [(0, 0), (0, 1), (1, 1), (1, 0)] 
        self.vi_coords = [(0.02, 0.02), (0.02, 0.98), (0.98, 0.98), (0.98, 0.02)] 
        self.f_indices = ((0, 1, 2), (2, 3, 0))
        self.tex_coords = ((0, 0), (0, 1), (1, 1), (1, 0))
        self.update(context)
        self.create_batch()
        self.line_shader.bind()
        self.line_shader.uniform_float("colour", (0, 0, 0, 1)) 
        
        
    def update(self, context):
        scene = context.scene
        svp = scene.vi_params
        self.cao = context.active_object
        self.col = svp.vi_scatt_col
        
        if self.cao and self.cao.vi_params.get('ws'):
#            self.unit = svp.wind_type 
            zdata = array(self.cao.vi_params['ws']) if svp.wind_type == '0' else array(self.cao.vi_params['wd'])
            zmax = nmax(zdata) if svp.vi_scatt_max == '0' else svp.vi_scatt_max_val
            zmin = nmin(zdata) if svp.vi_scatt_min == '0' else svp.vi_scatt_min_val
            (title, cbtitle) = ('Wind Speed', 'Speed (m/s)') if svp.wind_type == '0' else ('Wind Direction', u'Direction (\u00B0)')
            self.plt = plt
            draw_dhscatter(self, scene, self.cao.vi_params['days'], self.cao.vi_params['hours'], zdata, title, 'Days', 'Hours', cbtitle, zmin, zmax)  
            save_plot(self, scene, self.image)
            
    def create_batch(self):
        base_vertex_shader = '''
            uniform mat4 ModelViewProjectionMatrix;
            uniform vec2 spos;
            uniform vec2 size;
            in vec2 position;
            
            void main()
                {
                   float xpos = spos[0] + position[0] * size[0];
                   float ypos = spos[1] + position[1] * size[1]; 
                   gl_Position = ModelViewProjectionMatrix * vec4(int(xpos), int(ypos), 0.0f, 1.0f);
                }
        '''
        
        base_fragment_shader = '''
            uniform vec4 colour;
            out vec4 FragColour;
            
            void main()
                {
                    FragColour = colour;
                }
           
            '''
        image_vertex_shader = '''
            uniform mat4 ModelViewProjectionMatrix;
            uniform vec2 spos;
            uniform vec2 size;
            in vec2 texCoord;
            in vec2 position;
            out vec2 texCoord_interp;
            
            void main()
            {
              float xpos = spos[0] + position[0] * size[0];
              float ypos = spos[1] + position[1] * size[1]; 
              gl_Position = ModelViewProjectionMatrix * vec4(int(xpos), int(ypos), 0.0f, 1.0f);
//              gl_Position.z = 1.0;
              texCoord_interp = texCoord;
            }
        '''
        
        image_fragment_shader = '''
            in vec2 texCoord_interp;
            out vec4 fragColor;
            
            uniform sampler2D image;
            
            void main()
            {
              fragColor = texture(image, texCoord_interp);
            }

        '''
        self.base_shader = gpu.types.GPUShader(base_vertex_shader, base_fragment_shader) 
        self.line_shader = gpu.types.GPUShader(base_vertex_shader, base_fragment_shader)
        self.image_shader = gpu.types.GPUShader(image_vertex_shader, image_fragment_shader)
        self.base_batch = batch_for_shader(self.base_shader, 'TRIS', {"position": self.v_coords}, indices = self.f_indices)
        self.line_batch = batch_for_shader(self.line_shader, 'LINE_LOOP', {"position": self.v_coords})
        self.image_batch = batch_for_shader(self.image_shader, 'TRI_FAN', {"position": self.vi_coords, "texCoord": self.tex_coords})
        
    def draw(self, ah, aw):
        self.ah = ah
        self.aw = aw
        
        if self.expand:
            if self.resize:
                self.xdiff = self.lepos[0] - self.lspos[0]
                self.ydiff = self.lepos[1] - self.lspos[1]
            elif self.move:
                self.lspos[1] = self.lepos[1] - self.ydiff
                self.lepos[0] = self.lspos[0] + self.xdiff
            if self.lepos[1] > ah:
                self.lspos[1] = ah - self.ydiff 
                self.lepos[1] = ah
            if self.lepos[0] > aw:
                self.lspos[0] = aw - self.xdiff   
                self.lepos[0] = aw
                
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
            self.image_shader.bind()
            self.image_shader.uniform_float("size", (self.xdiff, self.ydiff))
            self.image_shader.uniform_float("spos", self.lspos) 
            im = bpy.data.images[self.image]
            if im.gl_load():
                raise Exception()
            bgl.glActiveTexture(bgl.GL_TEXTURE0)
            bgl.glBindTexture(bgl.GL_TEXTURE_2D, im.bindcode)
            self.image_shader.uniform_int("image", 0)
            self.image_batch.draw(self.image_shader)
            
    def show_plot(self, context):
        show_plot(self, context)
        

        
#class wr_table(Base_Display):
#    def __init__(self, pos, width, height, iname, xdiff, ydiff):
#        Base_Display.__init__(self, pos, width, height, iname, xdiff, ydiff)
#        self.fontdpi = int(0.15 * ydiff)
#        
#    def update(self, context):
#        self.cao = context.active_object
#        covp = self.cao.vi_params
#        if self.cao and covp.get('ws'):
#            self.rcarray = array(covp['table'])  
#        
#    def drawopen(self, context):
#        draw_table(self)
        
#def wr_disp(self, context, simnode):
#    try:
#        if self._handle_wr_disp:
#            width, height = context.region.width, context.region.height
#            self.legend.draw(context, width, height)
#            self.dhscatter.draw(context, width, height)
#            self.table.draw(context, width, height)
#    except:
#        pass

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
        self.results_bar = results_bar(('legend.png', 'table.png', 'scatter.png'), 300, area)
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
                
            if svp.vi_disp_refresh:
                self.dhscatter.update(context)
                self.dhscatter.draw(ah, aw)
                svp.vi_disp_refresh = 0
                context.area.tag_redraw()
            
            if redraw:
                context.area.tag_redraw()        
        return {'PASS_THROUGH'}
    
    def draw_wrnum(self, context):
        region = context.region
        rh = region.height  
        rw = region.width
        self.results_bar.draw(rh)
        self.legend.draw(rh, rw)
        self.table.draw(rh, rw)
        self.dhscatter.draw(rh, rw)
        
class VIEW3D_OT_SVFDisplay(bpy.types.Operator):
    '''Display results legend and stats in the 3D View'''
    bl_idname = "view3d.svfdisplay"
    bl_label = "Shadow study metric display"
    bl_description = "Display shadow study metrics"
    bl_register = True
    bl_undo = False
    
    def invoke(self, context, event):        
        region = context.region
        svp = context.scene.vi_params
        try:
            bpy.types.SpaceView3D.draw_handler_remove(self.draw_handle_svfnum, 'WINDOW')
#            bpy.types.SpaceView3D.draw_handler_remove(self._handle_ss_disp, 'WINDOW')
        except:
            pass
        res_coll = bpy.data.collections.get('LiVi Results')
        if res_coll:
            for o in res_coll.objects:
                bpy.data.objects.remove(o)
        svp.vi_display = 1
        svp['viparams']['vidisp'] = 'svf'
        self.simnode = bpy.data.node_groups[svp['viparams']['restree']].nodes[svp['viparams']['resnode']]
        li_display(self, self.simnode)
        self.results_bar = results_bar(('legend.png',), 300, region)
        self.legend = svf_legend(context, 'Sky View (%)', [305, region.height - 80], region.width, region.height, 100, 400)
        self.legend_num = linumdisplay(self, context)
        self.height = region.height
        self.draw_handle_svfnum = bpy.types.SpaceView3D.draw_handler_add(self.draw_svfnum, (context, ), 'WINDOW', 'POST_PIXEL')
        self.cao = context.active_object
        region.tag_redraw()
        context.window_manager.modal_handler_add(self)
        return {'RUNNING_MODAL'}
    
    def draw_svfnum(self, context):  
        region = context.region    
        self.results_bar.draw(region.height)
        self.legend.draw(context)
        self.legend_num.draw(context)
        
    def modal(self, context, event):    
        scene = context.scene
        svp = scene.vi_params
        redraw = 0
           
        if svp.vi_display == 0 or svp['viparams']['vidisp'] != 'svf' or event.type == 'ESC':
            svp.vi_display = 0
            bpy.types.SpaceView3D.draw_handler_remove(self.draw_handle_svfnum, 'WINDOW')
            context.area.tag_redraw()
            return {'CANCELLED'}

        if event.type != 'INBETWEEN_MOUSEMOVE' and context.region and context.area.type == 'VIEW_3D' and context.region.type == 'WINDOW':            
            mx, my = event.mouse_region_x, event.mouse_region_y 

#           for ui_element in (self.ui_elements):
                
            # Legend routine 
            if self.legend.ispos[0] < mx < self.legend.iepos[0] and self.legend.ispos[1] < my < self.legend.iepos[1]:
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
        region = context.region
        self.height = region.height
        self.width = region.width
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
        self.results_bar = results_bar(('legend.png', 'scatter.png'), 300, region)
        self.legend = ss_legend(context, 'Sunlit (%)', [305, self.height - 80], self.width, self.height, 125, 400)
        self.num_display = linumdisplay(self, context)
        self.dhscatter = ss_scatter(context, [355, self.height - 80], self.width, self.height, 600, 200)
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
        
class VIEW3D_OT_Li_DBSDF(bpy.types.Operator):
    bl_idname = "view3d.bsdf_display"
    bl_label = "BSDF display"
    bl_description = "Display BSDF"
    bl_register = True
    bl_undo = False
        
    def modal(self, context, event):
        scene = context.scene
        svp = scene.vi_params
        ah, aw = context.region.height, context.region.width
        redraw = 0
#        self.bsdf.f_colours = [(1, 1, 1, 1)] * (721 + 8 * 720)   
        if svp.vi_display == 0 or svp['viparams']['vidisp'] != 'bsdf_panel' or event.type == 'ESC':
            svp['viparams']['vidisp'] = 'bsdf'
            svp.vi_display = 0
            bpy.types.SpaceView3D.draw_handler_remove(self._handle_bsdfnum, 'WINDOW')
            context.area.tag_redraw()
            return {'CANCELLED'}

        if self.bsdf.expand and any((self.bsdf.leg_max != svp.vi_bsdfleg_max, 
                                     self.bsdf.leg_min != svp.vi_bsdfleg_min, 
                                     self.bsdf.col != svp.vi_leg_col, 
                                     self.bsdf.scale_select != svp.vi_bsdfleg_scale,
                                     self.bsdf.direc != svp.vi_bsdf_direc)):
            self.bsdf.col = svp.vi_leg_col
            self.bsdf.leg_max = svp.vi_bsdfleg_max
            self.bsdf.leg_min = svp.vi_bsdfleg_min
            self.bsdf.scale_select = svp.vi_bsdfleg_scale
            self.bsdf.direc = svp.vi_bsdf_direc
            self.bsdf.plot(context)
            context.region.tag_redraw()
            print('passing', (self.bsdf.leg_max != svp.vi_bsdfleg_max, self.bsdf.leg_min != svp.vi_bsdfleg_min, self.bsdf.col != svp.vi_leg_col, self.bsdf.scale_select != svp.vi_bsdfleg_scale))
            return {'PASS_THROUGH'}
        
        if event.type != 'INBETWEEN_MOUSEMOVE' and context.region and context.area.type == 'VIEW_3D' and context.region.type == 'WINDOW':            
            mx, my = event.mouse_region_x, event.mouse_region_y            
            # Legend routine 
            if self.bsdf.ispos[0] < mx < self.bsdf.iepos[0] and self.bsdf.ispos[1] < my < self.bsdf.iepos[1]:
                self.bsdf.hl = (0.8, 0.8, 0.8, 0.8) 
                redraw = 1
                if event.type == 'LEFTMOUSE':
                    if event.value == 'RELEASE':
                        self.bsdf.expand = 0 if self.bsdf.expand else 1

                context.region.tag_redraw()
                
            elif self.bsdf.expand and ((mx - self.bsdf.lspos[0] - 175)**2 + (my - ah + 300)**2)**0.5 <= 170:
                dist = ((mx - self.bsdf.lspos[0] - 175)**2 + (my - ah + 300)**2)**0.5
                for radius in self.bsdf.radii:
                    if dist <= radius:
                        redraw = 1
                        mring = self.bsdf.radii.index(radius)
#                        angles = [0] + [si*2*pi/self.bsdf.segments[ring] for si in range(1, self.bsdf.segments[mring])]
                        mangle = atan2((-my + ah - 300), (mx - self.bsdf.lspos[0] - 175)) + pi + pi/self.bsdf.segments[mring]
#                        phi = atan2(-my + self.bsdf.centre[1], mx - self.bsdf.centre[0]) + pi
#                        for a in angles:
#                            if abs(a - angle) < 2*pi/self.bsdf.segments[ring]:
                        if mring == 0:
                            ms = 1
#                            self.bsdf.f_colours[:721] = [(0.5, 0.5, 0.5, 1)] * 721
                        else:
                            ms = int(mangle*self.bsdf.segments[mring]/(2*pi)) + 1 if int(mangle*self.bsdf.segments[mring]/(2*pi)) < self.bsdf.segments[mring] else 1 
#                            self.bsdf.f_colours[int(720*mring + ms*720/self.bsdf.segments[mring] + 1) + 1: int(720*mring + (ms+1)*2*pi/self.bsdf.segments[mring]) + 2] = [(0.5, 0.5, 0.5, 1)] * ((int(720*mring + (ms+1)*2*pi/self.bsdf.segments[mring]) + 2) - (int(720*mring + ms*2*pi/self.bsdf.segments[mring] + 1) + 1))
                        msegment = sum([self.bsdf.segments[ri] for ri in range(mring)]) + ms
                        break
                self.bsdf.cr, self.bsdf.crs, self.bsdf.cseg = mring, ms, msegment
                self.bsdf.create_batch('sel')
                
                if event.type == 'LEFTMOUSE':
                    
#                    self.bsdf.f_colours[int(720*mring + s*720/self.bsdf.segments[mring] + 1) + 1: int(720*mring + (s+1)*720/self.bsdf.segments[mring]) + 2] = [(0.5, 0.5, 0.5, 1)] * ((int(720*mring + (s+1)*2*pi/self.bsdf.segments[mring]) + 2) - (int(720*mring + s*2*pi/self.bsdf.segments[mring] + 1) + 1))

#                    if self.bsdf.expand:
                    redraw = 1
                    self.bsdf.create_batch('arc')
                    if self.bsdf.cseg != self.bsdf.sseg:
                        self.bsdf.sr, self.bsdf.srs, self.bsdf.sseg = mring, ms, msegment
                        self.bsdf.plot(context)

                
                if redraw:
                    context.region.tag_redraw()
                return {'RUNNING_MODAL'} 
        return {'PASS_THROUGH'}                    
#                self.bsdf.draw(context)

#        if event.type != 'INBETWEEN_MOUSEMOVE' and context.region and context.area.type == 'VIEW_3D' and context.region.type == 'WINDOW':  
#            if context.scene['viparams']['vidisp'] != 'bsdf_panel':
#                self.remove(context)
#                context.scene['viparams']['vidisp'] = self.olddisp
#                return {'CANCELLED'}
#
#            mx, my = event.mouse_region_x, event.mouse_region_y
#            if self.bsdf.spos[0] < mx < self.bsdf.epos[0] and self.bsdf.spos[1] < my < self.bsdf.epos[1]:
#                self.bsdf.hl = (0, 1, 1, 1)  
#                
#                if event.type == 'LEFTMOUSE':
#                    if event.value == 'PRESS':
#                        self.bsdfpress = 1
#                        self.bsdfmove = 0
#                        return {'RUNNING_MODAL'}
#                    elif event.value == 'RELEASE':
#                        if not self.bsdfmove:
#                            self.bsdf.expand = 0 if self.bsdf.expand else 1
#                        self.bsdfpress = 0
#                        self.bsdfmove = 0
#                        context.area.tag_redraw()
#                        return {'RUNNING_MODAL'}
#                        
#                elif event.type == 'ESC':
#                    self.remove(context)
#                    context.scene['viparams']['vidisp'] = self.olddisp
#                    return {'CANCELLED'}                   
#                elif self.bsdfpress and event.type == 'MOUSEMOVE':
#                     self.bsdfmove = 1
#                     self.bsdfpress = 0
#                            
#            elif abs(self.bsdf.lepos[0] - mx) < 10 and abs(self.bsdf.lspos[1] - my) < 10:
#                self.bsdf.hl = (0, 1, 1, 1) 
#                if event.type == 'LEFTMOUSE':
#                    if event.value == 'PRESS':
#                        self.bsdf.resize = 1
#                    if self.bsdf.resize and event.value == 'RELEASE':
#                        self.bsdf.resize = 0
#                    return {'RUNNING_MODAL'}  
#            
#            elif all((self.bsdf.expand, self.bsdf.lspos[0] + 0.45 * self.bsdf.xdiff < mx < self.bsdf.lspos[0] + 0.8 * self.bsdf.xdiff, self.bsdf.lspos[1] + 0.06 * self.bsdf.ydiff < my < self.bsdf.lepos[1] - 5)):
#                if event.type == 'LEFTMOUSE' and event.value == 'PRESS':
#                    self.bsdf.plt.show()
#            
#            else:
#                for butrange in self.bsdf.buttons:
#                    if self.bsdf.buttons[butrange][0] - 0.0075 * self.bsdf.xdiff < mx < self.bsdf.buttons[butrange][0] + 0.0075 * self.bsdf.xdiff and self.bsdf.buttons[butrange][1] - 0.01 * self.bsdf.ydiff < my < self.bsdf.buttons[butrange][1] + 0.01 * self.bsdf.ydiff:
#                        if event.type == 'LEFTMOUSE' and event.value == 'PRESS' and self.bsdf.expand:
#                            if butrange in ('Front', 'Back'):
#                                self.bsdf.dir_select = butrange
#                            elif butrange in ('Visible', 'Solar', 'Discrete'):
#                                self.bsdf.rad_select = butrange
#                            elif butrange in ('Transmission', 'Reflection'):
#                                self.bsdf.type_select = butrange
#                            self.bsdf.plot(context)
#
#                self.bsdf.hl = (1, 1, 1, 1)
#                                
#            if event.type == 'MOUSEMOVE':                
#                if self.bsdfmove:
#                    self.bsdf.pos = [mx, my]
#                    context.area.tag_redraw()
#                    return {'RUNNING_MODAL'}
#                if self.bsdf.resize:
#                    self.bsdf.lepos[0], self.bsdf.lspos[1] = mx, my
#            
#            if self.bsdf.expand and self.bsdf.lspos[0] < mx < self.bsdf.lepos[0] and self.bsdf.lspos[1] < my < self.bsdf.lepos[1]:
#                theta, phi = xy2radial(self.bsdf.centre, (mx, my), self.bsdf.pw, self.bsdf.ph)
#                phi = atan2(-my + self.bsdf.centre[1], mx - self.bsdf.centre[0]) + pi
#
#                if theta < self.bsdf.radii[-1]:
#                    for ri, r in enumerate(self.bsdf.radii):
#                        if theta < r:
#                            break
#
#                    upperangles = [p * 2 * pi/self.bsdf.phis[ri] + pi/self.bsdf.phis[ri]  for p in range(int(self.bsdf.phis[ri]))]
#                    uai = 0
#
#                    if ri > 0:
#                        for uai, ua in enumerate(upperangles): 
#                            if phi > upperangles[-1]:
#                                uai = 0
#                                break
#                            if phi < ua:
#                                break
#
#                    self.bsdf.patch_hl = sum(self.bsdf.phis[0:ri]) + uai
#                    if event.type in ('LEFTMOUSE', 'RIGHTMOUSE')  and event.value == 'PRESS':                        
#                        self.bsdf.num_disp = 1 if event.type == 'RIGHTMOUSE' else 0    
#                        self.bsdf.patch_select = sum(self.bsdf.phis[0:ri]) + uai
#                        self.bsdf.plot(context)
#                        context.area.tag_redraw()
#                        return {'RUNNING_MODAL'}
#                        
#                else:
#                    self.bsdf.patch_hl = None
#                    
            

                
    def invoke(self, context, event):
        cao = context.active_object
        region = context.region
        scene = context.scene
        svp = scene.vi_params
#        width, height = area.width, area.height
        if cao and cao.active_material.get('bsdf') and cao.active_material.vi_params['bsdf']['xml'] and cao.active_material.vi_params['bsdf']['type'] == 'LBNL/Klems Full':
            bsdf = parseString(cao.active_material.vi_params['bsdf']['xml'])
#        coltype = [path.firstChild.data for path in bsdf.getElementsByTagName('ColumnAngleBasis')]
#        rowtype = [path.firstChild.data for path in bsdf.getElementsByTagName('RowAngleBasis')]

            svp['liparams']['bsdf_direcs'] = [(path.firstChild.data, path.firstChild.data, 'BSDF Direction') for path in bsdf.getElementsByTagName('WavelengthDataDirection')]
#            svp.vi_bsdf_direc = bpy.props.EnumProperty(items = direcs, name = "", description = "BSDf display direction")
#            width, height = context.region.width, context.region.height
            self.images = ['bsdf.png']
            self.results_bar = results_bar(self.images, 300, region)
            self.bsdf = bsdf_disp(context, '', [305, region.height - 80], region.width, region.height, 800, 375)
#            context.scene.bsdf = bsdf_disp(context, '', [305, region.height - 80], region.width, region.height, 800, 375)
#            self.bsdf.update(context)
#            self.bsdfpress, self.bsdfmove, self.bsdfresize = 0, 0, 0
            svp.vi_display = 1
            self._handle_bsdfnum = bpy.types.SpaceView3D.draw_handler_add(self.draw_bsdfnum, (context, ), 'WINDOW', 'POST_PIXEL')
#            self.olddisp = context.scene['viparams']['vidisp']
            context.window_manager.modal_handler_add(self)
#            if not context.scene.get('viparams'):
#                context.scene['viparams'] = {}
            svp['viparams']['vidisp'] = 'bsdf_panel'
            context.area.tag_redraw()  
#            self.results_bar = results_bar(('legend.png',), 300, area)
#            self.legend = svf_legend(context, 'Sky View (%)', [305, area.height - 80], area.width, area.height, 125, 400)
#            self.legend_num = linumdisplay(self, context)
#            self.height = area.height
#            self.draw_handle_svfnum = bpy.types.SpaceView3D.draw_handler_add(self.draw_svfnum, (context, ), 'WINDOW', 'POST_PIXEL')
#            self.cao = context.active_object
#            area.tag_redraw()
#            context.window_manager.modal_handler_add(self)
            return {'RUNNING_MODAL'}
            #return {'RUNNING_MODAL'}
        else:
            self.report({'ERROR'},"Selected material contains no BSDF information or contains the wrong BSDF type (only Klems is supported)")
            return {'CANCELLED'}
            
    def remove(self, context):
        self.bsdf.plt.close()
        bpy.types.SpaceView3D.draw_handler_remove(self._handle_bsdfnum, 'WINDOW')
        context.scene.vi_params['viparams']['vidisp'] = 'bsdf'
        bpy.data.images.remove(self.bsdf.gimage)
        context.area.tag_redraw()
    
    def draw_bsdfnum(self, context):
        self.results_bar.draw(context.region.height)
        self.bsdf.draw(context)

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
        region = context.region
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

#        lnd = linumdisplay(self, context)
   #     self._handle_pointres = bpy.types.SpaceView3D.draw_handler_add(lnd.draw, (context, ), 'WINDOW', 'POST_PIXEL')
        self.legend = livi_legend(context, svp['liparams']['unit'], [305, region.height - 80], region.width, region.height, 125, 400)

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

