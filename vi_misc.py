import bpy, os, inspect, gpu, bgl, blf
from gpu_extras.batch import batch_for_shader
from .vi_func import retcols
from numpy import array

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

#def draw_icon(self, sxpos): 
#    size = 40
#    yoffset = 40
#    image = bpy.data.images[self.image] 
#    shader = gpu.shader.from_builtin('2D_IMAGE')
#    batch = batch_for_shader(shader, 'TRI_FAN', 
#                             {"pos": ((sxpos, self.height - (yoffset + size)), 
#                                      (sxpos + size, self.height - (yoffset + size)), 
#                                      (sxpos + size, self.height - yoffset), 
#                                      (sxpos, self.height - yoffset)),
#                            "texCoord": ((0, 0), (1, 0), (1, 1), (0, 1))})
#    
#    if image.gl_load():
#        raise Exception()
#
#    bgl.glActiveTexture(bgl.GL_TEXTURE0)
#    bgl.glBindTexture(bgl.GL_TEXTURE_2D, image.bindcode)    
#    shader.bind()
#    shader.uniform_int("image", 0)
#    batch.draw(shader)

class Base_Display():
    def __init__(self, ipos, width, height, iname, xdiff, ydiff):
        self.ispos = ipos
        self.iepos = [ipos[0] + 40, ipos[1] + 40]
        self.lspos = [self.ispos[0], self.ispos[1] - ydiff - 20]
        self.lepos = [self.ispos[0] + xdiff, self.ispos[1] - 20]
        self.resize, self.move, self.expand, self.xdiff, self.ydiff = 0, 0, 0, xdiff, ydiff
        
        if iname not in bpy.data.images:
            bpy.data.images.load(os.path.join(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))), 'Images', iname))
            
        self.image = iname
        self.hl = [1, 1, 1, 1]
        self.cao = None
        self.xdiff, self.ydiff = xdiff, ydiff
        self.ah = height
        self.aw = width
        
#    def draw(self, context, width, height):
#        self.width, self.height = width, height
#        if self.lepos[1] > height:
#            self.lepos[1] = height
##        self.spos = (int(self.pos[0] - 25), int(self.pos[1] - 15))
##        self.epos = (int(self.pos[0] + 25), int(self.pos[1] + 15))
#        
#        if self.expand == 0:
#            self.drawclosed()
#            
#        if self.expand == 1:
#            self.drawopen(context)
#    
#    def draw_closed(self):
#        draw_icon(self, 305) 
        
class results_bar():
    def __init__(self, images, pos, area):
        self.images = images
        self.pos = pos
        self.ah = area.height
        self.aw = area.width
        self.shaders = [gpu.shader.from_builtin('2D_UNIFORM_COLOR'), gpu.shader.from_builtin('2D_UNIFORM_COLOR')]
        self.height = 0
        for im in images:
            if im not in bpy.data.images:
                bpy.data.images.load(os.path.join(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))), 'Images', im))
            self.shaders.append(gpu.shader.from_builtin('2D_IMAGE'))
#        self.wrbl_shader = gpu.shader.from_builtin('2D_UNIFORM_COLOR')
#        self.wrbf_shader = gpu.shader.from_builtin('2D_UNIFORM_COLOR')
#        self.leg_shader = gpu.shader.from_builtin('2D_IMAGE')
#        self.tab_shader = gpu.shader.from_builtin('2D_IMAGE')
#        self.hm_shader = gpu.shader.from_builtin('2D_IMAGE')
#        print(self.shaders)
        
    def draw(self, ah, no):
        v_coords = ((self.pos, ah - 85), (self.pos + no * 50, ah - 85),
                    (self.pos + no * 50, ah - 35), (self.pos, ah - 35))
        f_indices = ((0, 1, 2), (2, 3, 0))
        tex_coords = ((0, 0), (1, 0), (1, 1), (0, 1))
        
        if self.height != ah:
            self.batches = [batch_for_shader(self.shaders[1], 'TRIS', {"pos": v_coords}, indices = f_indices),
                            batch_for_shader(self.shaders[0], 'LINE_LOOP', {"pos": v_coords})]

            for i in range(no):
#                print(self.shaders[no+2])
                pos = ((self.pos + 5 + i * 50, ah - 80), (self.pos + 45 + i * 50, ah - 80),(self.pos + 45 + i * 50, ah - 40), (self.pos + 5 + i * 50, ah - 40))
                self.batches.append(batch_for_shader(self.shaders[i + 2], 'TRI_FAN', {"pos": pos, "texCoord": tex_coords}))
#                self.leg_batch = batch_for_shader(self.leg_shader, 'TRI_FAN', {"pos": leg_pos, "texCoord": tex_coords})
#                tab_pos = ((355, ah - 80), (395, ah - 80),(395, ah - 40), (355, ah - 40))
#                self.tab_batch = batch_for_shader(self.tab_shader, 'TRI_FAN', {"pos": tab_pos, "texCoord": tex_coords})
#                hm_pos = ((405, ah - 80), (445, ah - 80),(445, ah - 40), (405, ah - 40))
#                self.hm_batch = batch_for_shader(self.hm_shader, 'TRI_FAN', {"pos": hm_pos, "texCoord": tex_coords})
                self.height = ah
        for si, s in enumerate(self.shaders):            
            if si == 0:
                s.bind()
                s.uniform_float("color", (1, 1, 1, 1))
                self.batches[si].draw(s)
            elif si == 1:
                s.bind()
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
                s.bind()
                s.uniform_int("image", 0)
                self.batches[si].draw(s)
                
#        self.wrbf_shader.bind()
#        self.wrbf_shader.uniform_float("color", (1, 1, 1, 1))
#        self.wrbf_batch.draw(self.wrbf_shader)
#        self.wrbl_shader.bind()
#        self.wrbl_shader.uniform_float("color", (0, 0, 0, 1))
#        bgl.glEnable(bgl.GL_BLEND)
#        bgl.glEnable(bgl.GL_LINE_SMOOTH)
#        bgl.glLineWidth(1)        
#        self.wrbl_batch.draw(self.wrbl_shader)
#        bgl.glLineWidth(1)        
#        bgl.glDisable(bgl.GL_LINE_SMOOTH)
#        bgl.glDisable(bgl.GL_BLEND)
#        leg_image = bpy.data.images[self.images[0]]
#         
#        if leg_image.gl_load():
#            raise Exception()
#            
#        bgl.glActiveTexture(bgl.GL_TEXTURE0)
#        bgl.glBindTexture(bgl.GL_TEXTURE_2D, leg_image.bindcode)
#        self.leg_shader.bind()
#        self.leg_shader.uniform_int("image", 0)
#        self.leg_batch.draw(self.leg_shader)
#        
#        tab_image = bpy.data.images[self.images[1]]
#        if tab_image.gl_load():
#            raise Exception()
#            
#        bgl.glActiveTexture(bgl.GL_TEXTURE0)
#        bgl.glBindTexture(bgl.GL_TEXTURE_2D, tab_image.bindcode)
#        self.tab_shader.bind()
#        self.tab_shader.uniform_int("image", 0)
#        self.tab_batch.draw(self.tab_shader)
#        
#        hm_image = bpy.data.images[self.images[2]]
#        
#        if hm_image.gl_load():
#            raise Exception()
#            
#        bgl.glActiveTexture(bgl.GL_TEXTURE0)
#        bgl.glBindTexture(bgl.GL_TEXTURE_2D, hm_image.bindcode)
#        self.hm_shader.bind()
#        self.hm_shader.uniform_int("image", 0)
#        self.hm_batch.draw(self.hm_shader)


class wr_legend2(Base_Display):
    def __init__(self, context, unit, pos, width, height, iname, xdiff, ydiff):
        Base_Display.__init__(self, pos, width, height, iname, xdiff, ydiff)
        self.unit = unit
        self.font_id = 0
        self.update(context)
        self.create_batch()
        self.line_shader.bind()
        self.line_shader.uniform_float("colour", (0, 0, 0, 1))  
                
    def update(self, context):
        font_id = 0
        scene = context.scene
        simnode = bpy.data.node_groups[scene.vi_params['viparams']['restree']].nodes[scene.vi_params['viparams']['resnode']]        
        self.cao = context.active_object

        if self.cao and self.cao.get('VIType') and self.cao['VIType'] == 'Wind_Plane':            
            self.levels = self.cao['nbins']
            maxres = self.cao['maxres']
        else:
            self.levels = simnode['nbins']
            maxres = simnode['maxres']
        self.cols = retcols(mcm.get_cmap(scene.vi_leg_col), self.levels)
        
        if not scene.get('liparams'):
            scene.vi_display = 0
            return

        self.resvals = ['{0:.0f} - {1:.0f}'.format(2*i, 2*(i+1)) for i in range(simnode['nbins'])]
        self.resvals[-1] = self.resvals[-1][:-int(len('{:.0f}'.format(maxres)))] + "Inf"
        self.colours = [item for item in [self.cols[i] for i in range(self.levels)] for i in range(4)]
        blf.size(font_id, 12, 300)        
        self.titxdimen = blf.dimensions(font_id, self.unit)[0]
        self.resxdimen = blf.dimensions(font_id, self.resvals[-1])[0]
        self.mydimen = blf.dimensions(font_id, self.unit)[1]
        
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
            blf.size(self.font_id, 12, int(300/fontscale))
            blf.position(self.font_id, self.lspos[0] + (self.xdiff - blf.dimensions(self.font_id, self.unit)[0]) * 0.45, self.lepos[1] - 0.5 * (self.lh * self.ydiff) - blf.dimensions(self.font_id, self.unit)[1] * 0.3, 0) 
            blf.color(self.font_id, 0, 0, 0, 1)      
            blf.draw(self.font_id, self.unit)
            blf.shadow(self.font_id, 5, 0.8, 0.8, 0.8, 1)    
            blf.size(self.font_id, 12, int(250/fontscale))
            bgl.glDisable(bgl.GL_BLEND)
            bgl.glHint(bgl.GL_LINE_SMOOTH_HINT, bgl.GL_NICEST)
            
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
            vl_coords += [(0, i * lh), (0.45, i * lh), (0.45, (i + 1) * lh), (0, (i + 1) * lh)]

        self.base_batch = batch_for_shader(self.base_shader, 'TRIS', {"position": v_coords}, indices = f_indices)
        self.line_batch = batch_for_shader(self.line_shader, 'LINE_LOOP', {"position": vl_coords})
        self.col_batch = batch_for_shader(self.col_shader, 'TRIS', {"position": vl_coords[4:], "colour": self.colours}, indices = fl_indices)
        self.lh = lh
        
#def draw_legend(self, scene, unit):
#    font_id = 0
#    blf.enable(0, 4)
#    blf.enable(0, 8)
#    blf.shadow(font_id, 5, 0.7, 0.7, 0.7, 1)    
#    levels = len(self.resvals)
#    xdiff = self.lepos[0] - self.lspos[0]
#    ydiff = self.lepos[1] - self.lspos[1]
#    lh = ydiff/(levels + 1.25)   
#    blf.size(font_id, 12, 300)
#    titxdimen = blf.dimensions(font_id, unit)[0]
#    resxdimen = blf.dimensions(font_id, self.resvals[-1])[0]
#    mydimen = blf.dimensions(font_id, unit)[1]
#    fontscale = max(titxdimen/(xdiff * 0.9), resxdimen/(xdiff * 0.6), mydimen * 1.25/lh)
#    blf.size(font_id, 12, int(300/fontscale))
#
#    if not self.resize:
#        self.lspos = [self.spos[0], self.spos[1] - ydiff]
#        self.lepos = [self.lspos[0] + xdiff, self.spos[1]]            
#    else:
#        self.lspos = [self.spos[0], self.lspos[1]]
#        self.lepos = [self.lepos[0], self.spos[1]]
#    
#    bgl.glLineWidth(1)
#    self.legl_shader = gpu.shader.from_builtin('2D_UNIFORM_COLOR')
#    self.legf_shader = gpu.shader.from_builtin('2D_UNIFORM_COLOR')
#    self.legfc_shader = gpu.shader.from_builtin('2D_FLAT_COLOR')
#    colours = [item for item in [self.cols[i] for i in range(levels)] for i in range(4)]
#    v_coords = [(self.lspos[0], self.lspos[1]), (self.lspos[0], self.lepos[1]), (self.lepos[0], self.lepos[1]), (self.lepos[0], self.lspos[1])]
#    vl_coords = v_coords
#    f_indices = [(0, 1, 2), (2, 3, 0)]
#    fl1_indices = [tuple(array((0, 1, 2)) +4 * i) for i in range(levels)]
#    fl2_indices = [tuple(array((2, 3, 0)) +4 * i) for i in range(levels)]
#    fl_indices = list(fl1_indices) + list(fl2_indices)
#
#    for i in range(0, levels):
#        vl_coords += [(self.lspos[0], int(self.lspos[1] + i * lh)), (int(self.lspos[0] + xdiff * 0.4), int(self.lspos[1] + i * lh)), (int(self.lspos[0] + xdiff * 0.4), int(self.lspos[1] + (i + 1) * lh)), (self.lspos[0], int(self.lspos[1] + (i + 1) * lh))]
#
#    self.legl_batch = batch_for_shader(self.legl_shader, 'LINE_LOOP', {"pos": vl_coords})
#    self.legf_batch = batch_for_shader(self.legf_shader, 'TRIS', {"pos": v_coords}, indices = f_indices)
#    self.legfc_batch = batch_for_shader(self.legfc_shader, 'TRIS', {"pos": vl_coords[4:], "color": colours}, indices = fl_indices)
#    bgl.glEnable(bgl.GL_BLEND)
#    self.legf_shader.bind()
#    self.legf_shader.uniform_float("color", (self.hl))
#    self.legf_batch.draw(self.legf_shader)
#    bgl.glDisable(bgl.GL_BLEND)
#    
#    self.legfc_shader.bind()
#    self.legfc_batch.draw(self.legfc_shader)
#    
#    self.legl_shader.bind()
#    self.legl_shader.uniform_float("color", (0, 0, 0, 1))
#    self.legl_batch.draw(self.legl_shader)
#
#    blf.position(font_id, self.lspos[0] + (xdiff - blf.dimensions(font_id, unit)[0]) * 0.45, self.spos[1] - 0.5 * lh - blf.dimensions(font_id, unit)[1] * 0.3, 0) 
#    blf.color(font_id, 0, 0, 0, 1)      
#    blf.draw(font_id, unit)
##    blf.enable(0, blf.SHADOW)
##    blf.enable(0, blf.KERNING_DEFAULT)
##    blf.shadow(0, 5, 0, 0, 0, 0.7)
#    
##    bgl.glColor4f(*scene.vi_display_rp_fc)
#
#    blf.shadow(font_id, 5, 0.8, 0.8, 0.8, 1)
#    
#    blf.size(font_id, 12, int(250/fontscale))
#    bgl.glDisable(bgl.GL_BLEND)
#    
##    self.legl_shader = gpu.shader.from_builtin('2D_UNIFORM_COLOR')
#   
#    for i in range(levels):
#        num = self.resvals[i]
#        rgba = self.cols[i]
#        bgl.glHint(bgl.GL_LINE_SMOOTH_HINT, bgl.GL_NICEST)
#        ndimen = blf.dimensions(font_id, "{}".format(num))
#        blf.position(font_id, int(self.lepos[0] - xdiff * 0.05 - ndimen[0]), int(self.lspos[1] + i * lh) + int((lh - ndimen[1])*0.5), 0)
##        bgl.glColor4f(0, 0, 0, 1)
#        blf.draw(font_id, "{}".format(self.resvals[i]))
#    
#    bgl.glLineWidth(1)
##    bgl.glColor4f(0, 0, 0, 1)
#    blf.disable(0, 8)  
#    blf.disable(0, 4)