import bpy, os, inspect, gpu, bgl, blf, bmesh, mathutils
from gpu_extras.batch import batch_for_shader
from bpy_extras import view3d_utils
from .vi_func import retcols, logentry, setscenelivivals, selobj, cmap, skframe, retvpvloc, viewdesc, objmode, ret_res_vals, draw_index_distance, blf_props, leg_min_max, retdp, move_to_coll
from .vi_dicts import unit_dict
from . import livi_export
from numpy import array
from numpy import sum as nsum
from numpy import min as nmin
from numpy import max as nmax
from numpy import append as nappend
from math import log10


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
    def __init__(self, images, pos, area):
        self.images = images
        self.pos = pos
        self.ah = area.height
        self.aw = area.width
        self.shaders = [gpu.shader.from_builtin('2D_UNIFORM_COLOR'), gpu.shader.from_builtin('2D_UNIFORM_COLOR')]
        self.height = 0
        self.f_indices = ((0, 1, 2), (2, 3, 0))
        self.tex_coords = ((0, 0), (1, 0), (1, 1), (0, 1))
        self.no = len(images)
        
        for im in images:
            if im not in bpy.data.images:
                bpy.data.images.load(os.path.join(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))), 'Images', im))
            self.shaders.append(gpu.shader.from_builtin('2D_IMAGE'))
        
    def draw(self, ah):
        v_coords = ((self.pos, ah - 85), (self.pos + self.no * 50, ah - 85),
                    (self.pos + self.no * 50, ah - 35), (self.pos, ah - 35))
        
        
        if self.height != ah:
            self.batches = [batch_for_shader(self.shaders[1], 'TRIS', {"pos": v_coords}, indices = self.f_indices),
                            batch_for_shader(self.shaders[0], 'LINE_LOOP', {"pos": v_coords})]

            for i in range(self.no):
                pos = ((self.pos + 5 + i * 50, ah - 80), (self.pos + 45 + i * 50, ah - 80),(self.pos + 45 + i * 50, ah - 40), (self.pos + 5 + i * 50, ah - 40))
                self.batches.append(batch_for_shader(self.shaders[i + 2], 'TRI_FAN', {"pos": pos, "texCoord": self.tex_coords}))
                self.height = ah
                
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
        self.cols = retcols(mcm.get_cmap(svp.vi_scatt_col), self.levels)
        
        if not scene.get('liparams'):
            svp.vi_display = 0
            return

        self.resvals = ['{0:.0f} - {1:.0f}'.format(2*i, 2*(i+1)) for i in range(simnode['nbins'])]
        self.resvals[-1] = self.resvals[-1][:-int(len('{:.0f}'.format(maxres)))] + "Inf"
        self.colours = [item for item in [self.cols[i] for i in range(self.levels)] for i in range(4)]
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
        
class svf_legend(Base_Display):
    def __init__(self, context, unit, pos, width, height, xdiff, ydiff):
        Base_Display.__init__(self, pos, width, height, xdiff, ydiff)
        self.unit = unit
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
        self.cols = retcols(mcm.get_cmap(svp.vi_leg_col), self.levels)
        resdiff = self.maxres - self.minres
        
        if not svp.get('liparams'):
            svp.vi_display = 0
            return

        resvals = [format(self.minres + i*(resdiff)/self.levels, '.2f') for i in range(self.levels + 1)] if self.scale == '0' else \
                        [format(self.minres + (1 - log10(i)/log10(self.levels + 1))*(resdiff), '.2f') for i in range(1, self.levels + 2)[::-1]]
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
            vl_coords += [(0, i * lh), (0.45, i * lh), (0.45, (i + 1) * lh), (0, (i + 1) * lh)]
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
        
        if self.cao and self.cao.get('ws'):
            self.rcarray = array(self.cao['table']) 
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
                            blf.position(fid, self.lspos[0] + 0.005 * self.xdiff, int(rowmids[r] - 0.4 * blf.dimensions(fid, 'H')[1]), 0)#.format(self.rcarray[r][c]))[1])), 0)#int(self.lepos[1] - rowoffset - rowheight * (r + 0.5)), 0)
                        else:
                            blf.position(fid, self.lspos[0] + ctws[c] - int(blf.dimensions(fid, '{}'.format(self.rcarray[r][c]))[0] * 0.5), int(rowmids[r] - 0.5 * blf.dimensions(fid, 'H')[1]), 0)
                        blf.draw(fid, '{}'.format(self.rcarray[r][c]))

            blf.disable(0, 8)
            blf.disable(0, 4)

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
        
        if self.cao and self.cao.get('ws'):
#            self.unit = svp.wind_type 
            zdata = array(self.cao['ws']) if svp.wind_type == '0' else array(self.cao['wd'])
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

def show_plot(self, context):
    try:
        self.plt.show()
        self.update(context)
       
    except Exception as e:
        logentry('Error showing matplotlib graph: {}'.format(e))
        
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
        move_to_coll(bpy.context, 'Livi Results', ores)
        
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
            print('stop', svp['viparams']['vidisp'])
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
            bpy.context.user_preferences.system.window_draw_method = bpy.context.user_preferences.system.window_draw_method
           
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
            
#                except Exception as e:
#                    print(e)
#                    pcs, depths, res = [], [], []
                
            self.allpcs = nappend(self.allpcs, array(pcs))
            self.alldepths = nappend(self.alldepths, array(depths))
            self.allres = nappend(self.allres, array(res))

        self.alldepths = self.alldepths/nmin(self.alldepths)        
        draw_index_distance(self.allpcs, self.allres, self.fontmult * svp.vi_display_rp_fs, svp.vi_display_rp_fc, svp.vi_display_rp_fsh, self.alldepths)    

#        except Exception as e:
#            logentry('Error in LiVi number display: {}'.format(e))
            
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