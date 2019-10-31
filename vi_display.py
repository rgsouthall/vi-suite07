import bpy, blf, mathutils, datetime, os, bgl, inspect, gpu
from gpu_extras.batch import batch_for_shader
from bpy_extras import view3d_utils
from .vi_func import ret_vp_loc, viewdesc, draw_index, draw_time, blf_props, retcols, drawloop, drawpoly
from math import pi
from numpy import array
from numpy import min as nmin
from numpy import max as nmax

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
    
def spnumdisplay(disp_op, context):
    scene = context.scene

    if bpy.data.objects.get('SPathMesh'):
        spob = bpy.data.objects['SPathMesh'] 
        ob_mat = spob.matrix_world
        mid_x, mid_y, width, height = viewdesc(context)
        vl = ret_vp_loc(context)
        blf_props(scene, width, height)
        
        if scene.vi_params.sp_hd:
            pvecs = [ob_mat@mathutils.Vector(p[:]) for p in spob['numpos'].values()]
            pvals = [int(p.split('-')[1]) for p in spob['numpos'].keys()]
            p2ds = [view3d_utils.location_3d_to_region_2d(context.region, context.region_data, p) for p in pvecs]
            vispoints = [pi for pi, p in enumerate(pvals) if p2ds[pi] and 0 < p2ds[pi][0] < width and 0 < p2ds[pi][1] < height and scene.ray_cast(context.view_layer, vl, pvecs[pi] - vl, distance = (pvecs[pi] - vl).length)[4] == spob]
            
            if vispoints:
                hs = [pvals[pi] for pi in vispoints]
                posis = [p2ds[pi] for pi in vispoints]                
                draw_index(posis, hs, scene.vi_params.display_rp_fs, scene.vi_params.display_rp_fc, scene.vi_params.display_rp_fsh)
                
        if [ob.get('VIType') == 'Sun' for ob in bpy.data.objects] and scene.vi_params['spparams']['suns'] == '0':
            sobs = [ob for ob in bpy.data.objects if ob.get('VIType') == 'Sun']
            
            if sobs and scene.vi_params.sp_td:
                sunloc = ob_mat@sobs[0].location
                solpos = view3d_utils.location_3d_to_region_2d(context.region, context.region_data, sunloc)
                
                try:
                    if 0 < solpos[0] < width and 0 < solpos[1] < height and not scene.ray_cast(context.view_layer, sobs[0].location + 0.05 * (vl - sunloc), vl - sunloc)[0]:
                        soltime = datetime.datetime.fromordinal(scene.solday)
                        soltime += datetime.timedelta(hours = scene.solhour)
                        sre = sobs[0].rotation_euler
                        blf_props(scene, width, height)
                        draw_time(solpos, soltime.strftime('  %d %b %X') + ' alt: {:.1f} azi: {:.1f}'.format(90 - sre[0]*180/pi, (180, -180)[sre[2] < -pi] - sre[2]*180/pi), 
                                   scene.vi_params.display_rp_fs, scene.vi_params.display_rp_fc, scene.vi_params.display_rp_fsh)
                        
                except Exception as e:
                    print(e)
        blf.disable(0, 4)
    else:
        return

class Base_Display():
    def __init__(self, pos, width, height, iname, xdiff, ydiff):
        self.pos = pos
        self.ispos = pos
        self.iepos = [pos[0] + 40, pos[1] + 40]
        self.spos = [int(self.pos[0] - 0.025 * width), int(self.pos[1] - 0.0125 * height)]
        self.epos = [int(self.pos[0] + 0.025 * width), int(self.pos[1] + 0.0125 * height)]  
        self.lspos = [self.spos[0], self.spos[1] - ydiff]
        self.lepos = [self.spos[0] + xdiff, self.spos[1]]
        self.lpos = (self.pos[0] + 0.2 * width, self.pos[1] - 0.2 * height)
        self.resize = 0
        self.press = 0
        self.move = 0
        self.expand = 0
        if iname not in bpy.data.images:
            bpy.data.images.load(os.path.join(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))), 'Images', iname))
        self.image = iname
        self.hl = [1, 1, 1, 1]
        self.cao = None
        self.xdiff, self.ydiff = xdiff, ydiff
        
    def draw(self, context, width, height):
        self.width, self.height = context.region.width, context.region.height
        if self.pos[1] > height:
            self.pos[1] = height
#        self.spos = (int(self.pos[0] - 25), int(self.pos[1] - 15))
#        self.epos = (int(self.pos[0] + 25), int(self.pos[1] + 15))
        
        if self.expand == 0:
            self.drawclosed()
            
        if self.expand == 1:
            self.drawopen(context)
    
    def drawclosed(self):
        draw_icon_new(self) 

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
        rgba = self.cols[i]
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
    IMAGE_NAME = "Untitled"
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
    
def draw_dhscatter(self, scene, x, y, z, tit, xlab, ylab, zlab, valmin, valmax):
    self.plt.close()
    self.col = scene.vi_leg_col
    x = [x[0] - 0.5] + [xval + 0.5 for xval in x] 
    y = [y[0] - 0.5] + [yval + 0.5 for yval in y]
    self.plt.figure(figsize=(6 + len(x)/len(y), 6))
    
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
    fileloc = os.path.join(scene['viparams']['newdir'], 'images', filename)
    self.plt.savefig(fileloc, pad_inches = 0.1)
    
    if filename not in [i.name for i in bpy.data.images]:
        self.gimage = filename
        bpy.data.images.load(fileloc)
    else:
        self.gimage = filename
        bpy.data.images[filename].reload()
        
    bpy.data.images[self.gimage].user_clear()

def show_plot(self):
    try:
        self.plt.show()
    except:
        pass
        
class wr_legend(Base_Display):
    def __init__(self, pos, width, height, iname, xdiff, ydiff):
        Base_Display.__init__(self, pos, width, height, iname, xdiff, ydiff)
        
    def update(self, context):
        scene = context.scene
        simnode = bpy.data.node_groups[scene.vi_params['viparams']['restree']].nodes[scene.vi_params['viparams']['resnode']]        
        self.cao = context.active_object
        covp = self.cao.vi_params

        if self.cao and covp.get('VIType') and covp['VIType'] == 'Wind_Plane':            
            levels = covp['nbins']
            maxres = covp['maxres']
        else:
            levels = simnode['nbins']
            maxres = simnode['maxres']
        self.cols = retcols(mcm.get_cmap(scene.vi_leg_col), levels)
        
        if not scene.get('liparams'):
            scene.vi_display = 0
            return

        self.resvals = ['{0:.0f} - {1:.0f}'.format(2*i, 2*(i+1)) for i in range(simnode['nbins'])]
        self.resvals[-1] = self.resvals[-1][:-int(len('{:.0f}'.format(maxres)))] + "Inf"  
        
    def drawopen(self, context):
        draw_legend(self, context.scene, 'Speed (m/s)')
        
class wr_scatter(Base_Display):
    def __init__(self, pos, width, height, iname, xdiff, ydiff):
        Base_Display.__init__(self, pos, width, height, iname, xdiff, ydiff)
        self.unit = '0'
        
    def update(self, context):
        self.cao = context.active_object
        covp = self.cao.vi_params
        if self.cao and covp.get('ws'):
            self.unit = context.scene.wind_type 
            zdata = array(covp['ws']) if context.scene.wind_type == '0' else array(covp['wd'])
            (title, cbtitle) = ('Wind Speed', 'Speed (m/s)') if context.scene.wind_type == '0' else ('Wind Direction', u'Direction (\u00B0)')
            self.plt = plt
            draw_dhscatter(self, context.scene, covp['days'], covp['hours'], zdata, title, 'Days', 'Hours', cbtitle, nmin(zdata), nmax(zdata))  
            save_plot(self, context.scene, 'scatter.png')
        
    def drawopen(self, context):
        draw_image(self, 0)
        
    def show_plot(self):
        show_plot(self)
        
class wr_table(Base_Display):
    def __init__(self, pos, width, height, iname, xdiff, ydiff):
        Base_Display.__init__(self, pos, width, height, iname, xdiff, ydiff)
        self.fontdpi = int(0.15 * ydiff)
        
    def update(self, context):
        self.cao = context.active_object
        covp = self.cao.vi_params
        if self.cao and covp.get('ws'):
            self.rcarray = array(covp['table'])  
        
    def drawopen(self, context):
        draw_table(self)
        
def wr_disp(self, context, simnode):
    try:
        if self._handle_wr_disp:
            width, height = context.region.width, context.region.height
            self.legend.draw(context, width, height)
            self.dhscatter.draw(context, width, height)
            self.table.draw(context, width, height)
    except:
        pass
    

