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


import bpy, glob, os, inspect, datetime, shutil, time, math, mathutils, sys
from bpy.props import EnumProperty, FloatProperty, IntProperty, BoolProperty, StringProperty, FloatVectorProperty
from bpy.types import NodeTree, Node, NodeSocket
from nodeitems_utils import NodeCategory, NodeItem
from subprocess import Popen
from .vi_func import socklink, socklink2, uvsocklink, newrow, epwlatilongi, nodeid, nodeinputs, remlink, rettimes, sockhide, selobj, cbdmhdr, cbdmmtx
from .vi_func import hdrsky, nodecolour, facearea, retelaarea, iprop, bprop, eprop, fprop, sunposlivi, retdates, validradparams, retpmap
#from .envi_func import retrmenus, resnameunits, enresprops, epentry, epschedwrite, processf, get_mat, get_con_node
#from .livi_export import livi_sun, livi_sky, livi_ground, hdrexport
#from .envi_mat import envi_materials, envi_constructions, envi_layer, envi_layertype, envi_con_list


#envi_mats = envi_materials()
#envi_cons = envi_constructions()

class ViNetwork(NodeTree):
    '''A node tree for VI-Suite analysis.'''
    bl_idname = 'ViN'
    bl_label = 'VI-Suite Nodes'
    bl_icon = 'NODETREE'
    viparams = {}
        
class ViNodes:
    @classmethod
    def poll(cls, ntree):
        return ntree.bl_idname == 'ViN'

class ViLoc(Node, ViNodes):
    '''Node describing a geographical location manually or with an EPW file'''
    bl_idname = 'ViLoc'
    bl_label = 'VI Location'
    bl_icon = 'FORCE_WIND'

    def updatelatlong(self, context):
        context.space_data.edit_tree == ''
#        print(NodeTree.get_from_context(context))
#        print(self.id_data)
        scene = context.scene
        nodecolour(self, self.ready())
        reslists = []

        if self.loc == '1':
            entries = []
            addonfolder = os.path.basename(os.path.dirname(os.path.abspath(__file__)))
            vi_prefs = bpy.context.preferences.addons['{}'.format(addonfolder)].preferences
        
            if vi_prefs and os.path.isdir(bpy.path.abspath(vi_prefs.epweath)):
                epwpath = bpy.path.abspath(vi_prefs.epweath)
            else:
                epwpath = os.path.dirname(os.path.abspath(__file__)) + '/EPFiles/Weather/'
    
            for wfile in glob.glob(epwpath+"/*.epw"):
                with open(wfile, 'r') as wf:
                    for wfl in wf.readlines():
                        if wfl.split(',')[0].upper() == 'LOCATION':
                            entries.append((wfile, '{} - {}'.format(wfl.split(',')[3], wfl.split(',')[1]), 'Weather Location'))
                            break
            self['entries'] = entries if entries else [('None', 'None', 'None')]
            
            if os.path.isfile(self.weather):            
                with open(self.weather, 'r') as epwfile:                
                    self['frames'] = ['0']
                    epwlines = epwfile.readlines()[8:]
                    epwcolumns = list(zip(*[epwline.split(',') for epwline in epwlines]))
                    self['year'] = 2015 if len(epwlines) == 8760 else 2016
                    times = ('Month', 'Day', 'Hour', 'DOS')
                    
                    for t, ti in enumerate([' '.join(epwcolumns[c]) for c in range(1,4)] + [' '.join(['{}'.format(int(d/24) + 1) for d in range(len(epwlines))])]):
                        reslists.append(['0', 'Time', '', times[t], ti])
                        
                    for c in {"Temperature ("+ u'\u00b0'+"C)": 6, 'Humidity (%)': 8, "Direct Solar (W/m"+u'\u00b2'+")": 14, "Diffuse Solar (W/m"+u'\u00b2'+")": 15,
                              'Wind Direction (deg)': 20, 'Wind Speed (m/s)': 21}.items():
                        reslists.append(['0', 'Climate', '', c[0], ' '.join([cdata for cdata in list(epwcolumns[c[1]])])])
    
                    self.outputs['Location out']['epwtext'] = epwfile.read()
                    self.outputs['Location out']['valid'] = ['Location', 'Vi Results']
            else:
                self.outputs['Location out']['epwtext'] = ''
                self.outputs['Location out']['valid'] = ['Location']

        socklink2(self.outputs['Location out'], self.id_data)
        self['reslists'] = reslists
        (scene.latitude, scene.longitude) = epwlatilongi(context.scene, self) if self.loc == '1' and self.weather != 'None' else (scene.latitude, scene.longitude)

        for node in [l.to_node for l in self.outputs['Location out'].links]:
            node.update()
                
    def retentries(self, context):
        try:
            return [tuple(e) for e in self['entries']]
        except:
            return [('None', 'None','None' )]
                  
    weather: EnumProperty(name='', items=retentries, update=updatelatlong)
    loc: EnumProperty(items=[("0", "Manual", "Manual location"), ("1", "EPW ", "Get location from EPW file")], name = "", description = "Location", default = "0", update = updatelatlong)
    maxws: FloatProperty(name="", description="Max wind speed", min=0, max=90, default=0)
    minws: FloatProperty(name="", description="Min wind speed", min=0, max=90, default=0)
    avws: FloatProperty(name="", description="Average wind speed", min=0, max=0, default=0)
    dsdoy: IntProperty(name="", description="", min=1, max=365, default=1)
    dedoy: IntProperty(name="", description="", min=1, max=365, default=365)

    def init(self, context):
        self['nodeid'] = nodeid(self)    
#        self.id_data.use_fake_user = True
#        bpy.data.node_groups[nodeid(self).split('@')[1]].use_fake_user = True
        
        self.outputs.new('ViLoc', 'Location out')
        self['year'] = 2015
        self['entries'] = [('None', 'None', 'None')] 
        NodeTree.get_from_context(context).use_fake_user = True

    def update(self):
        socklink2(self.outputs['Location out'], self.id_data)
        nodecolour(self, self.ready())
        
    def draw_buttons(self, context, layout):
        newrow(layout, "Source:", self, 'loc')
        
        if self.loc == "1":
            newrow(layout, "Weather file:", self, 'weather')
        else:
            newrow(layout, 'Latitude', context.scene.vi_params, "latitude")
            newrow(layout, 'Longitude', context.scene.vi_params, "longitude")
            
    def ready(self):
        if self.loc == '1' and not self.weather:
            return 1
        if any([link.to_node.bl_label in ('LiVi CBDM', 'EnVi Export') and self.loc != "1" for link in self.outputs['Location out'].links]):
            return 1
        return 0
        
class ViSPNode(Node, ViNodes):
    '''Node describing a VI-Suite sun path'''
    bl_idname = 'ViSPNode'
    bl_label = 'VI Sun Path'
#    bl_icon = 'LAMP'
    
    def nodeupdate(self, context):
        nodecolour(self, self['exportstate'] != [str(x) for x in (self.suns)])
    
    suns: EnumProperty(items = [('0', 'Single', 'Single sun'), ('1', 'Monthly', 'Monthly sun for chosen time'), ('2', 'Hourly', 'Hourly sun for chosen date')], name = '', description = 'Sunpath sun type', default = '0', update=nodeupdate)
#    res: FloatProperty(name="", description="Calc point offset", min=1, max=10, default=6, update = nodeupdate)

    def init(self, context):
#        self['nodeid'] = nodeid(self)
        self.inputs.new('ViLoc', 'Location in')
        self['exportstate'] = '0'
        nodecolour(self, 1)

    def draw_buttons(self, context, layout):
        if self.inputs['Location in'].links:
            newrow(layout, 'Suns:', self, 'suns')
            row = layout.row()
            row.operator("node.sunpath", text="Create Sun Path")#.nodeid = self['nodeid']
        else:
            row = layout.row()
            row.label(text="Connect location node")

    def export(self):
        nodecolour(self, 0)
        self['exportstate'] = [str(x) for x in (self.suns)]
        
    def update(self):
        pass

class ViWRNode(bpy.types.Node, ViNodes):
    '''Node describing a VI-Suite wind rose generator'''
    bl_idname = 'ViWRNode'
    bl_label = 'VI Wind Rose'
    bl_icon = 'FORCE_WIND'

    def nodeupdate(self, context):
        nodecolour(self, self['exportstate'] != [str(x) for x in (self.wrtype, self.sdoy, self.edoy)])

    wrtype: EnumProperty(items = [("0", "Hist 1", "Stacked histogram"), ("1", "Hist 2", "Stacked Histogram 2"), ("2", "Cont 1", "Filled contour"), ("3", "Cont 2", "Edged contour"), ("4", "Cont 3", "Lined contour")], name = "", default = '0', update = nodeupdate)
    sdoy: IntProperty(name = "", description = "Day of simulation", min = 1, max = 365, default = 1, update = nodeupdate)
    edoy: IntProperty(name = "", description = "Day of simulation", min = 1, max = 365, default = 365, update = nodeupdate)
    max_freq: EnumProperty(items = [("0", "Data", "Max frequency taken from data"), ("1", "Specified", "User entered value")], name = "", default = '0', update = nodeupdate)
    max_freq_val: FloatProperty(name = "", description = "Max frequency", min = 1, max = 100, default = 20, update = nodeupdate)

    def init(self, context):
        self.inputs.new('ViLoc', 'Location in')
        self['exportstate'] = ''
        nodecolour(self, 1)

    def draw_buttons(self, context, layout):
        if nodeinputs(self) and self.inputs[0].links[0].from_node.loc == '1':
            (sdate, edate) = retdates(self.sdoy, self.edoy, self.inputs[0].links[0].from_node['year'])
            newrow(layout, 'Type:', self, "wrtype")
            newrow(layout, 'Start day {}/{}:'.format(sdate.day, sdate.month), self, "sdoy")
            newrow(layout, 'End day {}/{}:'.format(edate.day, edate.month), self, "edoy")
            newrow(layout, 'Colour:', context.scene.vi_params, 'vi_leg_col')
            newrow(layout, 'Max frequency:', self, 'max_freq')
            if self.max_freq == '1':
               newrow(layout, 'Frequency:', self, 'max_freq_val') 
            row = layout.row()
            row.operator("node.windrose", text="Create Wind Rose")
        else:
            row = layout.row()
            row.label(text = 'Location node error')

    def export(self):
        nodecolour(self, 0)
        self['exportstate'] = [str(x) for x in (self.wrtype, self.sdoy, self.edoy, self.max_freq, self.max_freq_val)]
        
    def update(self):
        pass

class ViSVFNode(Node, ViNodes):
    '''Node for sky view factor analysis'''
    bl_idname = 'ViSVFNode'
    bl_label = 'VI SVF'
    bl_icon = 'COLOR'
    
    def nodeupdate(self, context):
        nodecolour(self, self['exportstate'] != [str(x) for x in (self.startframe, self.endframe, self.cpoint, self.offset, self.animmenu)])
    
    animtype = [('Static', "Static", "Simple static analysis"), ('Geometry', "Geometry", "Animated geometry analysis")]
    animmenu: EnumProperty(name="", description="Animation type", items=animtype, default = 'Static', update = nodeupdate)
    startframe: IntProperty(name = '', default = 0, min = 0, max = 1024, description = 'Start frame')
    endframe: IntProperty(name = '', default = 0, min = 0, max = 1024, description = 'End frame')
    cpoint: EnumProperty(items=[("0", "Faces", "Export faces for calculation points"),("1", "Vertices", "Export vertices for calculation points"), ],
            name="", description="Specify the calculation point geometry", default="0", update = nodeupdate)
    offset: FloatProperty(name="", description="Calc point offset", min=0.001, max=10, default=0.01, update = nodeupdate)
    signore: BoolProperty(name = '', default = 0, description = 'Ignore sensor surfaces', update = nodeupdate)
    skytype = [('0', "Tregenza", "145 Tregenza sky patches"), ('1', "Reinhart 577", "577 Reinhart sky patches"), ('2', 'Reinhart 2305', '2305 Reinhart sky patches')]
    skypatches: EnumProperty(name="", description="Animation type", items=skytype, default = '0', update = nodeupdate)
    
    def init(self, context):
#        self['nodeid'] = nodeid(self)
        self['goptions'] = {}
        self.outputs.new('ViR', 'Results out')
        nodecolour(self, 1)

    def draw_buttons(self, context, layout):
        newrow(layout, 'Ignore sensor:', self, "signore")
        newrow(layout, 'Animation:', self, "animmenu")
        if self.animmenu != 'Static': 
            row = layout.row(align=True)
            row.alignment = 'EXPAND'
            row.label('Frames:')
            row.prop(self, 'startframe')
            row.prop(self, 'endframe')
        newrow(layout, 'Sky patches:', self, "skypatches")
        newrow(layout, 'Result point:', self, "cpoint")
        newrow(layout, 'Offset:', self, 'offset')
        row = layout.row()
        row.operator("node.svf", text="Sky View Factor")#.nodeid = self['nodeid']
     
    def preexport(self):
        self['goptions']['offset'] = self.offset
        
    def postexport(self, scene):
        nodecolour(self, 0)
        self.outputs['Results out'].hide = False if self.get('reslists') else True            
        self['exportstate'] = [str(x) for x in (self.startframe, self.endframe, self.cpoint, self.offset, self.animmenu)]

    def update(self):
        socklink2(self.outputs['Results out'], self.id_data)
        
class ViNodeCategory(NodeCategory):
    @classmethod
    def poll(cls, context):
        return context.space_data.tree_type == 'ViN'

class ViLocSock(NodeSocket):
    '''Vi Location socket'''
    bl_idname = 'ViLoc'
    bl_label = 'Location socket'
    valid = ['Location']

    def draw(self, context, layout, node, text):
        layout.label(text = text)

    def draw_color(self, context, node):
        return (0.45, 1.0, 0.45, 1.0)
    
    def ret_valid(self, node):
        return ['Location']
    
class ViR(NodeSocket):
    '''Vi results socket'''
    bl_idname = 'ViR'
    bl_label = 'VI results'

    valid = ['Vi Results']
    link_limit = 1

    def draw(self, context, layout, node, text):
        layout.label(text = text)

    def draw_color(self, context, node):
        return (0.0, 1.0, 0.0, 0.75)
    
    def ret_valid(self, node):
        return ['Vi Results']
        
####################### Vi Nodes Categories ##############################

viexnodecat = []
                
vifilenodecat = []
vinodecat = [NodeItem("ViSPNode", label="Sun Path"), NodeItem("ViWRNode", label="Wind Rose"), NodeItem("ViSVFNode", label="Sky View")]

vigennodecat = []

vidisnodecat = []
vioutnodecat = []
viimnodecat = []
viinnodecat = [NodeItem("ViLoc", label="VI Location")]

vinode_categories = [ViNodeCategory("Output", "Output Nodes", items=vioutnodecat), 
                     ViNodeCategory("Edit", "Edit Nodes", items=vifilenodecat), 
                     ViNodeCategory("Image", "Image Nodes", items=viimnodecat), 
                     ViNodeCategory("Display", "Display Nodes", items=vidisnodecat), 
                     ViNodeCategory("Generative", "Generative Nodes", items=vigennodecat), 
                     ViNodeCategory("Analysis", "Analysis Nodes", items=vinodecat), 
                     ViNodeCategory("Process", "Process Nodes", items=viexnodecat), 
                     ViNodeCategory("Input", "Input Nodes", items=viinnodecat)]



