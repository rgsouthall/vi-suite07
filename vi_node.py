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


import bpy, glob, os, inspect, datetime, shutil, time, math, mathutils, sys, json
from bpy.props import EnumProperty, FloatProperty, IntProperty, BoolProperty, StringProperty, FloatVectorProperty
from bpy.types import NodeTree, Node, NodeSocket
from nodeitems_utils import NodeCategory, NodeItem
from subprocess import Popen
from .vi_func import socklink, socklink2, uvsocklink, uvsocklink2, newrow, epwlatilongi, nodeinputs, remlink, rettimes, sockhide, selobj
from .vi_func import nodecolour, facearea, retelaarea, iprop, bprop, eprop, fprop, sunposlivi, retdates
from .vi_func import delobj, logentry
from .livi_func import hdrsky, cbdmhdr, cbdmmtx, retpmap, validradparams
from .envi_func import retrmenus, resnameunits, enresprops, epentry, epschedwrite, processf, get_mat, get_con_node, get_con_node2
from .livi_export import livi_sun, livi_sky, livi_ground, hdrexport
from .envi_mat import envi_materials, envi_constructions, envi_layer, envi_layertype, envi_con_list
#from .vi_dicts import e1ddict

envi_mats = envi_materials()
envi_cons = envi_constructions()

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

# Input nodes
        
class No_Loc(Node, ViNodes):
    '''Node describing a geographical location manually or with an EPW file'''
    bl_idname = 'No_Loc'
    bl_label = 'VI Location'
    bl_icon = 'FORCE_WIND'

    def updatelatlong(self, context):
        context.space_data.edit_tree == ''
        scene = context.scene
        svp = scene.vi_params
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
        (svp.latitude, svp.longitude) = epwlatilongi(context.scene, self) if self.loc == '1' and self.weather != 'None' else (svp.latitude, svp.longitude)

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
#        self['nodeid'] = nodeid(self)    
#        self.id_data.use_fake_user = True
#        bpy.data.node_groups[nodeid(self).split('@')[1]].use_fake_user = True
        
        self.outputs.new('So_Vi_Loc', 'Location out')
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

# Export Nodes
class No_Li_Geo(Node, ViNodes):
    '''Node describing a LiVi geometry export node'''
    bl_idname = 'No_Li_Geo'
    bl_label = 'LiVi Geometry'

    def nodeupdate(self, context):
        nodecolour(self, self['exportstate'] != [str(x) for x in (self.animated, self.startframe, self.endframe, self.cpoint, self.offset, self.fallback)])

    cpoint: EnumProperty(items=[("0", "Faces", "Export faces for calculation points"),("1", "Vertices", "Export vertices for calculation points"), ],
            name="", description="Specify the calculation point geometry", default="0", update = nodeupdate)
    offset: FloatProperty(name="", description="Calc point offset", min = 0.001, max = 1, default = 0.01, update = nodeupdate)
    animated: BoolProperty(name="", description="Animated analysis", default = 0, update = nodeupdate)
    startframe: IntProperty(name="", description="Start frame for animation", min = 0, default = 0, update = nodeupdate)
    endframe: IntProperty(name="", description="End frame for animation", min = 0, default = 0, update = nodeupdate)
    fallback: BoolProperty(name="", description="Enforce simple geometry export", default = 0, update = nodeupdate)
    
    def init(self, context):
        self['exportstate'] = ''
#        self['nodeid'] = nodeid(self)
        self.outputs.new('So_Li_Geo', 'Geometry out')
        nodecolour(self, 1)

    def draw_buttons(self, context, layout):
        newrow(layout, 'Fallback:', self, 'fallback')
        newrow(layout, 'Animated:', self, 'animated')
        if self.animated:
            row = layout.row()
            row.label(text = 'Frames:')
            col = row.column()
            subrow = col.row(align=True)
            subrow.prop(self, 'startframe')
            subrow.prop(self, 'endframe')

        newrow(layout, 'Result point:', self, 'cpoint')
        newrow(layout, 'Offset:', self, 'offset')
        row = layout.row()
        row.operator("node.ligexport", text = "Export")

    def update(self):
        socklink(self.outputs['Geometry out'], self.id_data.name)

    def preexport(self, scene):
        self['Text'] = {}
        self['Options'] = {'offset': self.offset, 'fs': (scene.frame_current, self.startframe)[self.animated], 'fe': (scene.frame_current, self.endframe)[self.animated], 'cp': self.cpoint, 'anim': self.animated}
        
    def postexport(self, scene):
        self.id_data.use_fake_user = 1
        self['exportstate'] = [str(x) for x in (self.animated, self.startframe, self.endframe, self.cpoint, self.offset, self.fallback)]
        nodecolour(self, 0)

class No_Li_Con(Node, ViNodes):
    '''Node for creating a LiVi context'''
    bl_idname = 'No_Li_Con'
    bl_label = 'LiVi Context'

    def nodeupdate(self, context):
        scene = context.scene
        nodecolour(self, self['exportstate'] != [str(x) for x in (self.contextmenu, self.spectrummenu, self.canalysismenu, self.cbanalysismenu, 
                   self.animated, self.skymenu, self.shour, self.sdoy, self.startmonth, self.endmonth, self.damin, self.dasupp, self.dalux, self.daauto,
                   self.ehour, self.edoy, self.interval, self.hdr, self.hdrname, self.skyname, self.resname, self.turb, self.mtxname, self.cbdm_start_hour,
                   self.cbdm_end_hour, self.bambuildmenu)])
        if self.edoy < self.sdoy:
            self.edoy = self.sdoy
        if self.edoy == self.sdoy:
            if self.ehour < self.shour:
                self.ehour = self.shour
        
        self['skynum'] = int(self.skymenu)         
        suns = [ob for ob in scene.objects if ob.type == 'LIGHT' and ob.data.type == 'SUN'] 
                
        if self.contextmenu == 'Basic' and ((self.skyprog == '0' and self['skynum'] < 2) or (self.skyprog == '1' and self.epsilon > 1)):
            starttime = datetime.datetime(2015, 1, 1, int(self.shour), int((self.shour - int(self.shour))*60)) + datetime.timedelta(self.sdoy - 1) if self['skynum'] < 3 else datetime.datetime(2013, 1, 1, 12)                                       
            self['endframe'] = self.startframe + int(((24 * (self.edoy - self.sdoy) + self.ehour - self.shour)/self.interval)) if self.animated else [scene.frame_current]
            frames = range(self.startframe, self['endframe'] + 1) if self.animated else [scene.frame_current]
            scene.frame_start, scene.frame_end = self.startframe, frames[-1]
            
            if suns:
                sun = suns[0]
                sun['VIType'] = 'Sun'
                [delobj(bpy.context.view_layer, sun) for sun in suns[1:]]
#                bpy.ops.object.delete(use_global=False)
#                bpy.data.objects.remove(sun, do_unlink=True, do_id_user=True, do_ui_user=True)
#                scene.collection.objects.unlink(sun)
                
#                [scene.objects.unlink(o) for o in suns[1:]]
            else:
                bpy.ops.object.light_add(type='SUN')
                sun = bpy.context.object
                sun['VIType'] = 'Sun'

            if self.inputs['Location in'].links and suns:
                sunposlivi(scene, self, frames, sun, starttime)
        else:
            for so in suns:
                selobj(context.view_layer, so)
                bpy.ops.object.delete()
                
    spectrumtype =  [('0', "Visible", "Visible radiation spectrum calculation"), ('1', "Full", "Full radiation spectrum calculation")]                           
#    banalysistype = [('0', "Illu/Irrad/DF", "Illumninance/Irradiance/Daylight Factor Calculation"), ('1', "Glare", "Glare Calculation")]
    skylist = [("0", "Sunny", "CIE Sunny Sky description"), ("1", "Partly Coudy", "CIE Sunny Sky description"),
               ("2", "Coudy", "CIE Partly Cloudy Sky description"), ("3", "DF Sky", "Daylight Factor Sky description")]

    contexttype = [('Basic', "Basic", "Basic analysis"), ('Compliance', "Compliance", "Compliance analysis"), ('CBDM', "CBDM", "Climate based daylight modelling")]
    contextmenu: EnumProperty(name="", description="Contexttype type", items=contexttype, default = 'Basic', update = nodeupdate)
    animated: BoolProperty(name="", description="Animated sky", default=False, update = nodeupdate)
    offset: FloatProperty(name="", description="Calc point offset", min=0.001, max=1, default=0.01, update = nodeupdate)
    spectrummenu: EnumProperty(name="", description = "Visible/full radiation spectrum selection", items = spectrumtype, default = '0', update = nodeupdate)
    skyprog: EnumProperty(name="", items=[('0', "Gensky", "Basic sky creation"), ('1', "Gendaylit", "Perez sky creation"),
                                                     ("2", "HDR Sky", "HDR file sky"), ("3", "Radiance Sky", "Radiance file sky"), ("4", "None", "No Sky")], description="Specify sky creation", default="0", update = nodeupdate)
    epsilon: FloatProperty(name="", description="Hour of simulation", min=1, max=8, default=6.3, update = nodeupdate)
    delta: FloatProperty(name="", description="Hour of simulation", min=0.05, max=0.5, default=0.15, update = nodeupdate)
    skymenu: EnumProperty(name="", items=skylist, description="Specify the type of sky for the simulation", default="0", update = nodeupdate)
    gref: FloatProperty(name="", description="Ground reflectance", min=0.0, max=1.0, default=0.0, update = nodeupdate)
    gcol:FloatVectorProperty(size = 3, name = '', description="Ground colour", attr = 'Color', default = [0, 1, 0], subtype = 'COLOR', update = nodeupdate)
    shour: FloatProperty(name="", description="Hour of simulation", min=0, max=23.99, default=12, subtype='TIME', unit='TIME', update = nodeupdate)
    sdoy: IntProperty(name="", description="Day of simulation", min=1, max=365, default=1, update = nodeupdate)
    ehour: FloatProperty(name="", description="Hour of simulation", min=0, max=23.99, default=12, subtype='TIME', unit='TIME', update = nodeupdate)
    edoy: IntProperty(name="", description="Day of simulation", min=1, max=365, default=1, update = nodeupdate)
    interval: FloatProperty(name="", description="Site Latitude", min=1/60, max=24, default=1, update = nodeupdate)
    hdr: BoolProperty(name="", description="Export HDR panoramas", default=False, update = nodeupdate)
    skyname: StringProperty(name="", description="Name of the radiance sky file", default="", subtype="FILE_PATH", update = nodeupdate)
    mtxname: StringProperty(name="", description="Name of the radiance sky file", default="", subtype="FILE_PATH", update = nodeupdate)
    resname: StringProperty()
    turb: FloatProperty(name="", description="Sky Turbidity", min=1.0, max=5.0, default=2.75, update = nodeupdate)
    canalysistype = [('0', "BREEAM", "BREEAM HEA1 calculation"), ('1', "CfSH", "Code for Sustainable Homes calculation"), ('2', "Green Star", "Green Star Calculation"), ('3', "LEED", "LEED v4 Daylight calculation")]
    bambuildtype = [('0', "School", "School lighting standard"), ('1', "Higher Education", "Higher education lighting standard"), ('2', "Healthcare", "Healthcare lighting standard"), ('3', "Residential", "Residential lighting standard"), ('4', "Retail", "Retail lighting standard"), ('5', "Office & other", "Office and other space lighting standard")]
    lebuildtype = [('0', "Office/Education/Commercial", "Office/Education/Commercial lighting standard"), ('1', "Healthcare", "Healthcare lighting standard")]
    canalysismenu: EnumProperty(name="", description="Type of analysis", items = canalysistype, default = '0', update = nodeupdate)
    bambuildmenu: EnumProperty(name="", description="Type of building", items=bambuildtype, default = '0', update = nodeupdate)
    lebuildmenu: EnumProperty(name="", description="Type of building", items=lebuildtype, default = '0', update = nodeupdate)
    cusacc: StringProperty(name="", description="Custom Radiance simulation parameters", default="", update = nodeupdate)
    buildstorey: EnumProperty(items=[("0", "Single", "Single storey building"),("1", "Multi", "Multi-storey building")], name="", description="Building storeys", default="0", update = nodeupdate)
    cbanalysistype = [('0', "Exposure", "LuxHours/Irradiance Exposure Calculation"), ('1', "Hourly irradiance", "Irradiance for each simulation time step"), ('2', "DA/UDI/SDA/ASE", "Climate based daylighting metrics")]
    cbanalysismenu: EnumProperty(name="", description="Type of lighting analysis", items = cbanalysistype, default = '0', update = nodeupdate)
#    leanalysistype = [('0', "Light Exposure", "LuxHours Calculation"), ('1', "Radiation Exposure", "kWh/m"+ u'\u00b2' + " Calculation"), ('2', "Daylight Autonomy", "DA (%) Calculation")]
    sourcetype = [('0', "EPW", "EnergyPlus weather file"), ('1', "HDR", "HDR sky file")]
    sourcetype2 = [('0', "EPW", "EnergyPlus weather file"), ('1', "VEC", "Generated vector file")]
    sourcemenu: EnumProperty(name="", description="Source type", items=sourcetype, default = '0', update = nodeupdate)
    sourcemenu2: EnumProperty(name="", description="Source type", items=sourcetype2, default = '0', update = nodeupdate)
    hdrname: StringProperty(name="", description="Name of the composite HDR sky file", default="vi-suite.hdr", update = nodeupdate)
    hdrmap: EnumProperty(items=[("0", "Polar", "Polar to LatLong HDR mapping"),("1", "Angular", "Light probe or angular mapping")], name="", description="Type of HDR panorama mapping", default="0", update = nodeupdate)
    hdrangle: FloatProperty(name="", description="HDR rotation (deg)", min=0, max=360, default=0, update = nodeupdate)
    hdrradius: FloatProperty(name="", description="HDR radius (m)", min=0, max=5000, default=1000, update = nodeupdate)
    mtxname: StringProperty(name="", description="Name of the calculated vector sky file", default="", subtype="FILE_PATH", update = nodeupdate)
    weekdays: BoolProperty(name = '', default = False, update = nodeupdate)
    cbdm_start_hour:  IntProperty(name = '', default = 8, min = 1, max = 24, update = nodeupdate)
    cbdm_end_hour:  IntProperty(name = '', default = 20, min = 1, max = 24, update = nodeupdate)
    dalux:  IntProperty(name = '', default = 300, min = 1, max = 2000, update = nodeupdate)
    damin: IntProperty(name = '', default = 100, min = 1, max = 2000, update = nodeupdate)
    dasupp: IntProperty(name = '', default = 300, min = 1, max = 2000, update = nodeupdate)
    daauto: IntProperty(name = '', default = 3000, min = 1, max = 5000, update = nodeupdate)
    sdamin: IntProperty(name = '', default = 300, min = 1, max = 2000, update = nodeupdate) 
    asemax: IntProperty(name = '', default = 1000, min = 1, max = 2000, update = nodeupdate)
    startmonth: IntProperty(name = '', default = 1, min = 1, max = 12, description = 'Start Month', update = nodeupdate)
    endmonth: IntProperty(name = '', default = 12, min = 1, max = 12, description = 'End Month', update = nodeupdate)
    startframe: IntProperty(name = '', default = 0, min = 0, description = 'Start Frame', update = nodeupdate)

    def init(self, context):
        self['exportstate'], self['skynum'] = '', 0
#        self['nodeid'] = nodeid(self)
        self['whitesky'] = "void glow sky_glow \n0 \n0 \n4 1 1 1 0 \nsky_glow source sky \n0 \n0 \n4 0 0 1 180 \nvoid glow ground_glow \n0 \n0 \n4 1 1 1 0 \nground_glow source ground \n0 \n0 \n4 0 0 -1 180\n\n"
        self.outputs.new('So_Li_Con', 'Context out')
        self.inputs.new('So_Vi_Loc', 'Location in')
        nodecolour(self, 1)
        self.hdrname = ''
        self.skyname = ''
        self.mtxname = ''

    def draw_buttons(self, context, layout):
        newrow(layout, 'Context:', self, 'contextmenu')
        (sdate, edate) = retdates(self.sdoy, self.edoy, 2015)
        if self.contextmenu == 'Basic':            
            newrow(layout, "Program:", self, 'skyprog')
            if self.skyprog == '0':
                newrow(layout, "Sky type:", self, 'skymenu')
                newrow(layout, "Ground ref:", self, 'gref')
                newrow(layout, "Ground col:", self, 'gcol')
                
                if self.skymenu in ('0', '1', '2'):
                    newrow(layout, "Start hour:", self, 'shour')
                    newrow(layout, 'Start day {}/{}:'.format(sdate.day, sdate.month), self, "sdoy")
                    newrow(layout, "Animation;", self, 'animated')
                    
                    if self.animated:
                        newrow(layout, "Start frame:", self, 'startframe')
                        row = layout.row()
                        row.label(text = 'End frame:')
                        row.label(text = '{}'.format(self['endframe']))
                        newrow(layout, "End hour:", self, 'ehour')
                        newrow(layout, 'End day {}/{}:'.format(edate.day, edate.month), self, "edoy")
                        newrow(layout, "Interval (hours):", self, 'interval')
                    newrow(layout, "Turbidity", self, 'turb')
                
            elif self.skyprog == '1':
                newrow(layout, "Spectrum:", self, 'spectrummenu')
                newrow(layout, "Epsilon:", self, 'epsilon')
                newrow(layout, "Delta:", self, 'delta')
                newrow(layout, "Ground ref:", self, 'gref')
                newrow(layout, "Ground col:", self, 'gcol')
                newrow(layout, "Start hour:", self, 'shour')
                newrow(layout, 'Start day {}/{}:'.format(sdate.day, sdate.month), self, "sdoy")
                newrow(layout, "Animation;", self, 'animated')
                
                if self.animated:
                    newrow(layout, "Start frame:", self, 'startframe')
                    row = layout.row()
                    row.label(text = 'End frame:')
                    row.label(text = '{}'.format(self['endframe']))
                    newrow(layout, "End hour:", self, 'ehour')
                    newrow(layout, 'End day {}/{}:'.format(edate.day, edate.month), self, "edoy")
                    newrow(layout, "Interval (hours):", self, 'interval')
                
            elif self.skyprog == '2':
                row = layout.row()
                row.operator('node.hdrselect', text = 'HDR select')
                row.prop(self, 'hdrname')
                newrow(layout, "HDR format:", self, 'hdrmap')
                newrow(layout, "HDR rotation:", self, 'hdrangle')
                newrow(layout, "HDR radius:", self, 'hdrradius')

            elif self.skyprog == '3':
                row = layout.row()
                
                row.operator('node.skyselect', text = 'Sky select')
                row.prop(self, 'skyname')
            row = layout.row()

            if self.skyprog in ("0", "1"):
                newrow(layout, 'HDR:', self, 'hdr')

        elif self.contextmenu == 'Compliance':
            newrow(layout, "Standard:", self, 'canalysismenu')
            if self.canalysismenu == '0':
                newrow(layout, "Building type:", self, 'bambuildmenu')
                newrow(layout, "Storeys:", self, 'buildstorey')
            if self.canalysismenu == '2':
                newrow(layout, "Building type:", self, 'bambuildmenu')
            if self.canalysismenu == '3':
                newrow(layout, "Building type:", self, 'lebuildmenu')
                newrow(layout, 'Weekdays only:', self, 'weekdays')
                newrow(layout, 'Start hour:', self, 'cbdm_start_hour')
                newrow(layout, 'End hour:', self, 'cbdm_end_hour')
                newrow(layout, 'Source file:', self, 'sourcemenu2') 
                if self.sourcemenu2 == '1':
                    row = layout.row()
                    row.operator('node.mtxselect', text = 'Select MTX')
                    row = layout.row()
                    row.prop(self, 'mtxname')
            newrow(layout, 'HDR:', self, 'hdr')
                
        elif self.contextmenu == 'CBDM':
            newrow(layout, 'Type:', self, 'cbanalysismenu')
            if self.cbanalysismenu == '0':
                newrow(layout, "Spectrum:", self, 'spectrummenu')
            newrow(layout, 'Start day {}/{}:'.format(sdate.day, sdate.month), self, "sdoy")
            newrow(layout, 'End day {}/{}:'.format(edate.day, edate.month), self, "edoy")
            newrow(layout, 'Weekdays only:', self, 'weekdays')
            newrow(layout, 'Start hour:', self, 'cbdm_start_hour')
            newrow(layout, 'End hour:', self, 'cbdm_end_hour')
            row = layout.row()
            row.label(text = "--")
            
            if self.cbanalysismenu == '2':
                newrow(layout, '(s)DA Min lux:', self, 'dalux')
                row = layout.row()
                row.label(text = "--")
                newrow(layout, 'UDI Fell short (Max):', self, 'damin')
                newrow(layout, 'UDI Supplementry (Max):', self, 'dasupp')
                newrow(layout, 'UDI Autonomous (Max):', self, 'daauto')
                row = layout.row()
                row.label(text = "--")
                newrow(layout, 'ASE Lux level:', self, 'asemax')
                   
            if self.cbanalysismenu == '0':
                newrow(layout, 'Source file:', self, 'sourcemenu')
            else:
                newrow(layout, 'Source file:', self, 'sourcemenu2')
            row = layout.row()

            if self.sourcemenu2 == '1' and self.cbanalysismenu in ('1', '2'):
                newrow(layout, "MTX file:", self, 'mtxname')

            if self.sourcemenu == '1' and self.cbanalysismenu == '0':
                newrow(layout, "HDR file:", self, 'hdrname')
            else:
                newrow(layout, 'HDR:', self, 'hdr')
        
        if self.contextmenu == 'Basic':
            if int(self.skymenu) > 2 or (int(self.skymenu) < 3 and self.inputs['Location in'].links):
                row = layout.row()
                row.operator("node.liexport", text = "Export")
        elif self.contextmenu == 'Compliance' and self.canalysismenu != '3':
            row = layout.row()
            row.operator("node.liexport", text = "Export")
        elif (self.contextmenu == 'CBDM' and self.cbanalysismenu == '0' and self.sourcemenu2 == '1') or \
            (self.contextmenu == 'CBDM' and self.cbanalysismenu != '0' and self.sourcemenu == '1'):         
            row = layout.row()
            row.operator("node.liexport", text = "Export")  
        elif self.inputs['Location in'].links and self.inputs['Location in'].links[0].from_node.loc == '1' and self.inputs['Location in'].links[0].from_node.weather != 'None':
            row = layout.row()
            row.operator("node.liexport", text = "Export")

    def update(self):
        socklink(self.outputs['Context out'], self.id_data.name)
        if self.inputs.get('Location in'):
            self.nodeupdate(bpy.context) 
    
    def preexport(self):
        (interval, shour, ehour) = (1, self.cbdm_start_hour, self.cbdm_end_hour - 1) if self.contextmenu == 'CBDM' or (self.contextmenu == 'Compliance' and self.canalysismenu =='3') else (round(self.interval, 3), self.shour, self.ehour)        
        starttime = datetime.datetime(2015, 1, 1, 0) + datetime.timedelta(days = self.sdoy - 1) + datetime.timedelta(hours = shour)

        if self.contextmenu == 'CBDM' or (self.contextmenu == 'Basic' and self.animated):            
            endtime = datetime.datetime(2015, 1, 1, 0) + datetime.timedelta(days = self.edoy - 1)  + datetime.timedelta(hours = ehour)
        elif self.contextmenu == 'Compliance' and self.canalysismenu == '3':
            starttime = datetime.datetime(2015, 1, 1, 0) + datetime.timedelta(hours = shour)
            endtime = datetime.datetime(2015, 1, 1, 0) + datetime.timedelta(days = 364)  + datetime.timedelta(hours = ehour)
        else:
            endtime = starttime

        times = [starttime]
        ctime = starttime
        
        while ctime < endtime:
            ctime += datetime.timedelta(hours = interval)
            if (self.contextmenu == 'Compliance' and self.canalysismenu == '3') or self.contextmenu == 'CBDM':
                if shour <= ctime.hour <= ehour:
                    times.append(ctime)
            else:
                times.append(ctime)
               
        self.times = times 
        self.starttime = times[0]
        self.endtime = times[-1]
        self['skynum'] = int(self.skymenu)
        self['hours'] = 0 if not self.animated or int(self.skymenu) > 2  else (self.endtime-self.starttime).seconds/3600
        self['epwbase'] = os.path.splitext(os.path.basename(self.inputs['Location in'].links[0].from_node.weather)) if self.inputs['Location in'].links else ''
        self['Text'], self['Options'] = {}, {}
        self['watts'] = 0#1 if self.contextmenu == "CBDM" and self.cbanalysismenu in ('1', '2') else 0
        
    def export(self, scene, export_op):   
        svp = scene.vi_params
        self.startframe = self.startframe if self.animated and self.contextmenu == 'Basic' else scene.frame_current 
        self['endframe'] = self.startframe + int(((24 * (self.edoy - self.sdoy) + self.ehour - self.shour)/self.interval)) if self.contextmenu == 'Basic' and self.animated else scene.frame_current
        self['mtxfile'] = ''
        self['preview'] = 0
        
        if self.contextmenu == "Basic":  
            self['preview'] = 1
            
            if self.skyprog in ('0', '1'):
                self['skytypeparams'] = ("+s", "+i", "-c", "-b 22.86 -c")[self['skynum']] if self.skyprog == '0' else "-P {} {} -O {}".format(self.epsilon, self.delta, int(self.spectrummenu))

                for f, frame in enumerate(range(self.startframe, self['endframe'] + 1)):                  
                    skytext = livi_sun(scene, self, f) + livi_sky(self['skynum']) + livi_ground(*self.gcol, self.gref)
                    
                    if self['skynum'] < 2 or (self.skyprog == '1' and self.epsilon > 1):
                        if frame == self.startframe:
                            if 'SUN' in [ob.data.type for ob in scene.objects if ob.type == 'LIGHT' and ob.get('VIType')]:
                                sun = [ob for ob in scene.objects if ob.get('VIType') == 'Sun'][0]
                            else:
                                bpy.ops.object.light_add(type='SUN')
                                sun = bpy.context.object
                                sun['VIType'] = 'Sun'
 
                    if self.hdr:
                        hdrexport(scene, f, frame, self, skytext)                        
                    
                    self['Text'][str(frame)] = skytext

            elif self.skyprog == '2':
                if self.hdrname and os.path.isfile(self.hdrname):
                    if self.hdrname not in bpy.data.images:
                        bpy.data.images.load(self.hdrname)
                    self['Text'][str(scene.frame_current)] = hdrsky(self.hdrname, self.hdrmap, self.hdrangle, self.hdrradius)
                else:
                    export_op.report({'ERROR'}, "Not a valid HDR file")
                    return 'Error'
            
            elif self.skyprog == '3':
                if self.skyname and os.path.isfile(self.skyname):
                    shutil.copyfile(self.skyname, "{}-0.sky".format(svp['viparams']['filebase']))
                    
                    with open(self.skyname, 'r') as radfiler:
                        self['Text'][str(scene.frame_current)] = radfiler.read()
                        if self.hdr:
                            hdrexport(scene, 0, scene.frame_current, self, radfiler.read())
                else:
                    export_op.report({'ERROR'}, "Not a valid Radiance sky file")
                    return 'Error'

            elif self.skyprog == '4':
                self['Text'][str(scene.frame_current)] = ''
        
        elif self.contextmenu == "CBDM":
            if (self.cbanalysismenu =='0' and self.sourcemenu == '0') or (self.cbanalysismenu != '0' and self.sourcemenu2 == '0'):
                self['mtxfile'] = cbdmmtx(self, scene, self.inputs['Location in'].links[0].from_node, export_op)
            elif self.cbanalysismenu != '0' and self.sourcemenu2 == '1':
                self['mtxfile'] = self.mtxname

            if self.cbanalysismenu == '0' :
                self['preview'] = 1
                self['Text'][str(scene.frame_current)] = cbdmhdr(self, scene)
            else:
                self['Text'][str(scene.frame_current)] = "void glow sky_glow \n0 \n0 \n4 1 1 1 0 \nsky_glow source sky \n0 \n0 \n4 0 0 1 180 \nvoid glow ground_glow \n0 \n0 \n4 1 1 1 0 \nground_glow source ground \n0 \n0 \n4 0 0 -1 180\n\n"

                if self.sourcemenu2 == '0':
                    with open("{}.mtx".format(os.path.join(svp['viparams']['newdir'], self['epwbase'][0])), 'r') as mtxfile:
                        self['Options']['MTX'] = mtxfile.read()
                else:
                    with open(self.mtxname, 'r') as mtxfile:
                        self['Options']['MTX'] = mtxfile.read()
                if self.hdr:
                    self['Text'][str(scene.frame_current)] = cbdmhdr(self, scene)

        elif self.contextmenu == "Compliance":
            if self.canalysismenu in ('0', '1', '2'):            
                self['skytypeparams'] = ("-b 22.86 -c", "-b 22.86 -c", "-b 18 -u")[int(self.canalysismenu)]
                skyentry = livi_sun(scene, self, 0) + livi_sky(3)
                
                if self.canalysismenu in ('0', '1'):
                    self.starttime = datetime.datetime(2015, 1, 1, 12)
                    self['preview'] = 1
                    if self.hdr:
                        hdrexport(scene, 0, scene.frame_current, self, skyentry)
                else:
                    self.starttime = datetime.datetime(2015, 9, 11, 9)
                self['Text'][str(scene.frame_current)] = skyentry
            else:
                if self.sourcemenu2 == '0':
                    self['mtxfile'] = cbdmmtx(self, scene, self.inputs['Location in'].links[0].from_node, export_op)
                elif self.sourcemenu2 == '1':
                    self['mtxfile'] = self.mtxname
                
                self['Text'][str(scene.frame_current)] = "void glow sky_glow \n0 \n0 \n4 1 1 1 0 \nsky_glow source sky \n0 \n0 \n4 0 0 1 180 \nvoid glow ground_glow \n0 \n0 \n4 1 1 1 0 \nground_glow source ground \n0 \n0 \n4 0 0 -1 180\n\n"

                if self.sourcemenu2 == '0':
                    with open("{}.mtx".format(os.path.join(svp['viparams']['newdir'], self['epwbase'][0])), 'r') as mtxfile:
                        self['Options']['MTX'] = mtxfile.read()
                else:
                    with open(self.mtxname, 'r') as mtxfile:
                        self['Options']['MTX'] = mtxfile.read()
                if self.hdr:
                    self['Text'][str(scene.frame_current)] = cbdmhdr(self, scene)
                
    def postexport(self):  
        typedict = {'Basic': '0', 'Compliance': self.canalysismenu, 'CBDM': self.cbanalysismenu}
        unitdict = {'Basic': ("Lux", 'W/m2 (f)')[self.skyprog == '1' and self.spectrummenu =='1'], 'Compliance': ('DF (%)', 'DF (%)', 'DF (%)', 'sDA (%)')[int(self.canalysismenu)], 'CBDM': (('Mlxh', 'kWh')[int(self.spectrummenu)], 'kWh', 'DA (%)')[int(self.cbanalysismenu)]}
        btypedict = {'0': self.bambuildmenu, '1': '', '2': self.bambuildmenu, '3': self.lebuildmenu}
#        ('Luuminance, Irradiance, )
        self['Options'] = {'Context': self.contextmenu, 'Preview': self['preview'], 'Type': typedict[self.contextmenu], 'fs': self.startframe, 'fe': self['endframe'],
                    'anim': self.animated, 'shour': self.shour, 'sdoy': self.sdoy, 'ehour': self.ehour, 'edoy': self.edoy, 'interval': self.interval, 'buildtype': btypedict[self.canalysismenu], 'canalysis': self.canalysismenu, 'storey': self.buildstorey,
                    'bambuild': self.bambuildmenu, 'cbanalysis': self.cbanalysismenu, 'unit': unitdict[self.contextmenu], 'damin': self.damin, 'dalux': self.dalux, 'dasupp': self.dasupp, 'daauto': self.daauto, 'asemax': self.asemax, 'cbdm_sh': self.cbdm_start_hour, 
                    'cbdm_eh': self.cbdm_end_hour, 'weekdays': (7, 5)[self.weekdays], 'sourcemenu': (self.sourcemenu, self.sourcemenu2)[self.cbanalysismenu not in ('2', '3', '4', '5')],
                    'mtxfile': self['mtxfile'], 'times': [t.strftime("%d/%m/%y %H:%M:%S") for t in self.times]}
        nodecolour(self, 0)
        self['exportstate'] = [str(x) for x in (self.contextmenu, self.spectrummenu, self.canalysismenu, self.cbanalysismenu, 
                   self.animated, self.skymenu, self.shour, self.sdoy, self.startmonth, self.endmonth, self.damin, self.dasupp, self.dalux, self.daauto,
                   self.ehour, self.edoy, self.interval, self.hdr, self.hdrname, self.skyname, self.resname, self.turb, self.mtxname, self.cbdm_start_hour,
                   self.cbdm_end_hour, self.bambuildmenu)]

class No_Li_Im(Node, ViNodes):
    '''Node describing a LiVi image generation'''
    bl_idname = 'No_Li_Im'
    bl_label = 'LiVi Im'
    
    def nodeupdate(self, context):
        self["_RNA_UI"] = {"Processors": {"min": 1, "max": int(context.scene.vi_params['viparams']['nproc']), "name": ""}}
        nodecolour(self, self['exportstate'] != [str(x) for x in (self.camera, self.basename, self.illu, self.fisheye, self.fov,
                   self.mp, self['Processors'], self.processes, self.cusacc, self.simacc, self.pmap, self.pmapgno, self.pmapcno,
                   self.x, self.y)])
        if bpy.data.objects.get(self.camera):
            context.scene.camera = bpy.data.objects[self.camera]
        
        if self.simacc == '3':
            self.validparams = validradparams(self.cusacc)
            
    startframe: IntProperty(name = '', default = 0)
    endframe: IntProperty(name = '', default = 0)
    cusacc: StringProperty(
            name="", description="Custom Radiance simulation parameters", default="", update = nodeupdate)
    simacc: EnumProperty(items=[("0", "Low", "Low accuracy and high speed (preview)"),("1", "Medium", "Medium speed and accuracy"), ("2", "High", "High but slow accuracy"), 
                                           ("3", "Custom", "Edit Radiance parameters")], name="", description="Simulation accuracy", default="0", update = nodeupdate)
    rpictparams = (("-ab", 2, 3, 4), ("-ad", 256, 1024, 4096), ("-as", 128, 512, 2048), ("-aa", 0, 0, 0), ("-dj", 0, 0.7, 1), 
                   ("-ds", 0.5, 0.15, 0.15), ("-dr", 1, 3, 5), ("-ss", 0, 2, 5), ("-st", 1, 0.75, 0.1), ("-lw", 0.0001, 0.00001, 0.0000002), ("-lr", 3, 3, 4))
    pmap: BoolProperty(name = '', default = False, update = nodeupdate)
    pmapgno: IntProperty(name = '', default = 50000)
    pmapcno: IntProperty(name = '', default = 0)
    x: IntProperty(name = '', min = 1, max = 10000, default = 2000, update = nodeupdate)
    y: IntProperty(name = '', min = 1, max = 10000, default = 1000, update = nodeupdate)
    basename: StringProperty(name="", description="Base name for image files", default="", update = nodeupdate)
    run: BoolProperty(name = '', default = False) 
    illu: BoolProperty(name = '', default = True, update = nodeupdate)
    validparams: BoolProperty(name = '', default = True)
    mp: BoolProperty(name = '', default = False, update = nodeupdate)
    camera: StringProperty(description="Textfile to show", update = nodeupdate)
    fisheye: BoolProperty(name = '', default = 0, update = nodeupdate)
    fov: FloatProperty(name = '', default = 180, min = 1, max = 180, update = nodeupdate)
    processes: IntProperty(name = '', default = 1, min = 1, max = 1000, update = nodeupdate)
    
    def retframes(self):
        try:
            return min([c['fs'] for c in (self.inputs['Context in'].links[0].from_node['Options'], self.inputs['Geometry in'].links[0].from_node['Options'])]),\
                    max([c['fe'] for c in (self.inputs['Context in'].links[0].from_node['Options'], self.inputs['Geometry in'].links[0].from_node['Options'])])            
        except:
            return 0, 0
            
    def init(self, context):
        self['exportstate'] = ''
        self.inputs.new('So_Li_Geo', 'Geometry in')
        self.inputs.new('So_Li_Con', 'Context in')
        self.outputs.new('So_Li_Im', 'Image')
        self['Processors'] = 1
        
    def draw_buttons(self, context, layout):       
        sf, ef = self.retframes()
        row = layout.row()
        row.label(text = 'Frames: {} - {}'.format(sf, ef))
        layout.prop_search(self, 'camera', bpy.data, 'cameras', text='Camera*', icon='NONE')
        
        if all([sock.links for sock in self.inputs]) and self.camera:
            newrow(layout, 'Base name:', self, 'basename')        
            newrow(layout, 'Illuminance*:', self, 'illu')
            newrow(layout, 'Fisheye*:', self, 'fisheye')
    
            if self.fisheye:
                newrow(layout, 'FoV*:', self, 'fov')
            
            newrow(layout, 'Accuracy:', self, 'simacc')
    
            if self.simacc == '3':
                newrow(layout, "Radiance parameters:", self, 'cusacc')
            newrow(layout, 'Photon map*:', self, 'pmap')
    
            if self.pmap:
               newrow(layout, 'Global photons*:', self, 'pmapgno')
               newrow(layout, 'Caustic photons*:', self, 'pmapcno')

               
            if self.simacc != '3' or (self.simacc == '3' and self.validparams) and not self.run:
                row = layout.row()
                row.operator("node.radpreview", text = 'Preview') 
            newrow(layout, 'X resolution*:', self, 'x')
            newrow(layout, 'Y resolution*:', self, 'y')            
            
            if sys.platform != 'win32':
                newrow(layout, 'Multi-thread:', self, 'mp')
                if self.mp:
                    row = layout.row()
                    row.prop(self, '["Processors"]')
                    newrow(layout, 'Processes:', self, 'processes')
            if (self.simacc != '3' or (self.simacc == '3' and self.validparams)) and not self.run:
                row = layout.row()
                row.operator("node.radimage", text = 'Image')
        
    def update(self):        
        self.run = 0
        
    def presim(self):
        scene = bpy.context.scene
        if sys.platform == 'win32':
            self.mp = 0
        pmaps = []
        sf, ef, = self.retframes()
        self['frames'] = range(sf, ef + 1)
        self['viewparams'] = {str(f): {} for f in self['frames']}
        self['pmparams'] = {str(f): {} for f in self['frames']}
        self['pmaps'], self['pmapgnos'], self['pmapcnos'] = {}, {}, {}
        self['coptions'] = self.inputs['Context in'].links[0].from_node['Options']
        self['goptions'] = self.inputs['Geometry in'].links[0].from_node['Options']
        self['radfiles'], self['reslists'] = {}, [[]]
        self['radparams'] = self.cusacc if self.simacc == '3' else (" {0[0]} {1[0]} {0[1]} {1[1]} {0[2]} {1[2]} {0[3]} {1[3]} {0[4]} {1[4]} {0[5]} {1[5]} {0[6]} {1[6]} {0[7]} {1[7]} {0[8]} {1[8]} {0[9]} {1[9]} {0[10]} {1[10]} ".format([n[0] for n in self.rpictparams], [n[int(self.simacc)+1] for n in self.rpictparams]))
        self['basename'] = self.basename if self.basename else 'image'
        
        for frame in self['frames']:
            scene.frame_set(frame)
            scene.camera = bpy.data.objects[self.camera]
            cam = bpy.data.objects[self.camera]
            cang = cam.data.angle*180/math.pi if not self.fisheye else self.fov
            vh = cang if self.x >= self.y else cang * self.x/self.y 
            vv = cang if self.x < self.y else cang * self.y/self.x 
            vd = (0.001, 0, -1*cam.matrix_world[2][2]) if (round(-1*cam.matrix_world[0][2], 3), round(-1*cam.matrix_world[1][2], 3)) == (0.0, 0.0) else [-1*cam.matrix_world[i][2] for i in range(3)]
            pmaps.append(self.pmap)
            self['pmapgnos'][str(frame)] = self.pmapgno
            self['pmapcnos'][str(frame)] = self.pmapcno
            self['pmparams'][str(frame)]['amentry'], self['pmparams'][str(frame)]['pportentry'], self['pmparams'][str(frame)]['cpentry'], self['pmparams'][str(frame)]['cpfileentry'] = retpmap(self, frame, scene)
                
            if self.fisheye and self.fov == 180:
                self['viewparams'][str(frame)]['-vth'] = ''
                
            (self['viewparams'][str(frame)]['-vh'], self['viewparams'][str(frame)]['-vv']) = (self.fov, self.fov) if self.fisheye else (vh, vv)
            self['viewparams'][str(frame)]['-vd'] = ' '.join(['{:.3f}'.format(v) for v in vd])
            self['viewparams'][str(frame)]['-x'], self['viewparams'][str(frame)]['-y'] = self.x, self.y
            if self.mp:
                self['viewparams'][str(frame)]['-X'], self['viewparams'][str(frame)]['-Y'] = self.processes, 1
            self['viewparams'][str(frame)]['-vp'] = '{0[0]:.3f} {0[1]:.3f} {0[2]:.3f}'.format(cam.location)
            self['viewparams'][str(frame)]['-vu'] = '{0[0]:.3f} {0[1]:.3f} {0[2]:.3f}'.format(cam.matrix_world.to_quaternion()@mathutils.Vector((0, 1, 0)))
            
            if self.illu:
                self['viewparams'][str(frame)]['-i'] = ''
        self['pmaps'] = pmaps
        self.run = 1
        nodecolour(self, 1)
                
    def postsim(self, images):
        self['images'] = images
        self.run = 0
        self['exportstate'] = [str(x) for x in (self.camera, self.basename, self.illu, self.fisheye, self.fov,
            self.mp, self['Processors'], self.processes, self.cusacc, self.simacc, self.pmap, self.pmapgno, self.pmapcno,
            self.x, self.y)]
        nodecolour(self, 0)  

class No_Li_Gl(Node, ViNodes):
    '''Node describing a LiVi glare analysis'''
    bl_idname = 'No_Li_Gl'
    bl_label = 'LiVi Glare analysis' 

    def nodeupdate(self, context):
        nodecolour(self, self['exportstate'] != [str(x) for x in (self.hdrname)])

    hdrname: StringProperty(name="", description="Base name of the Glare image", default="", update = nodeupdate)    
    gc: FloatVectorProperty(size = 3, name = '', attr = 'Color', default = [1, 0, 0], subtype = 'COLOR', update = nodeupdate)
    rand: BoolProperty(name = '', default = True, update = nodeupdate)

    def init(self, context):
        self['exportstate'] = ''
        self.inputs.new('So_Li_Im', 'Image')
                
    def draw_buttons(self, context, layout):
        if self.inputs['Image'].links and os.path.isfile(bpy.path.abspath(self.inputs['Image'].links[0].from_node['images'][0])):
            newrow(layout, 'Base name:', self, 'hdrname')
            newrow(layout, 'Random:', self, 'rand')
            if not self.rand:
                newrow(layout, 'Colour:', self, 'gc')
            row = layout.row()
            row.operator("node.liviglare", text = 'Glare')
    
    def presim(self):
        self['hdrname'] = self.hdrname if self.hdrname else 'glare'

    def sim(self):
        for im in self['images']:
            with open(im+'glare', 'w') as glfile:
                Popen('evalglare {}'.format(im), stdout = glfile)

    def postsim(self):
        self['exportstate'] = [str(x) for x in (self.hdrname, self.colour, self.lmax, self.unit, self.nscale, self.decades, 
                   self.legend, self.lw, self.lh, self.contour, self.overlay, self.bands)]
        nodecolour(self, 0)

class No_Li_Fc(Node, ViNodes):
    '''Node describing a LiVi false colour image generation'''
    bl_idname = 'No_Li_Fc'
    bl_label = 'LiVi False Colour Image' 

    def nodeupdate(self, context):
        nodecolour(self, self['exportstate'] != [str(x) for x in (self.basename, self.colour, self.lmax, self.unit, self.nscale, self.decades, 
                   self.legend, self.lw, self.lh, self.contour, self.overlay, self.bands, self.ofile, self.hdrfile)])

    basename: StringProperty(name="", description="Base name of the falsecolour image(s)", default="", update = nodeupdate)    
    colour: EnumProperty(items=[("0", "Default", "Default color mapping"), ("1", "Spectral", "Spectral color mapping"), ("2", "Thermal", "Thermal colour mapping"), ("3", "PM3D", "PM3D colour mapping"), ("4", "Eco", "Eco color mapping")],
            name="", description="Simulation accuracy", default="0", update = nodeupdate)             
    lmax: IntProperty(name = '', min = 0, max = 100000, default = 1000, update = nodeupdate)
    unit: EnumProperty(items=[("0", "Lux", "Spectral color mapping"),("1", "Candelas", "Thermal colour mapping"), ("2", "DF", "PM3D colour mapping"), ("3", "Irradiance(v)", "PM3D colour mapping")],
            name="", description="Unit", default="0", update = nodeupdate)
    nscale: EnumProperty(items=[("0", "Linear", "Linear mapping"),("1", "Log", "Logarithmic mapping")],
            name="", description="Scale", default="0", update = nodeupdate)
    decades: IntProperty(name = '', min = 1, max = 5, default = 2, update = nodeupdate)
    unitdict = {'0': 'Lux', '1': 'cd/m2', '2': 'DF', '3': 'W/m2'}
    unitmult = {'0': 179, '1': 179, '2': 1.79, '3': 1}
    legend: BoolProperty(name = '', default = True, update = nodeupdate)
    lw: IntProperty(name = '', min = 1, max = 1000, default = 100, update = nodeupdate)
    lh: IntProperty(name = '', min = 1, max = 1000, default = 200, update = nodeupdate)
    contour: BoolProperty(name = '', default = False, update = nodeupdate)
    overlay: BoolProperty(name = '', default = False, update = nodeupdate)
    bands: BoolProperty(name = '', default = False, update = nodeupdate)
    coldict = {'0': 'def', '1': 'spec', '2': 'hot', '3': 'pm3d', '4': 'eco'}
    divisions: IntProperty(name = '', min = 1, max = 50, default = 8, update = nodeupdate)
    ofile: StringProperty(name="", description="Location of the file to overlay", default="", subtype="FILE_PATH", update = nodeupdate)
    hdrfile: StringProperty(name="", description="Location of the file to overlay", default="", subtype="FILE_PATH", update = nodeupdate)
    disp: FloatProperty(name = '', min = 0.0001, max = 10, default = 1, precision = 4, update = nodeupdate)
    
    def init(self, context):
        self['exportstate'] = ''
        self.inputs.new('ViLiI', 'Image')
                
    def draw_buttons(self, context, layout):
        if not self.inputs['Image'].links or not self.inputs['Image'].links[0].from_node['images'] or not os.path.isfile(bpy.path.abspath(self.inputs['Image'].links[0].from_node['images'][0])):
            row = layout.row()
            row.prop(self, 'hdrfile')
        if (self.inputs['Image'].links and self.inputs['Image'].links[0].from_node['images'] and os.path.isfile(bpy.path.abspath(self.inputs['Image'].links[0].from_node['images'][0]))) or os.path.isfile(self.hdrfile): 
            newrow(layout, 'Base name:', self, 'basename')
            newrow(layout, 'Unit:', self, 'unit')
            newrow(layout, 'Colour:', self, 'colour')
            newrow(layout, 'Divisions:', self, 'divisions')
            newrow(layout, 'Legend:', self, 'legend')
            
            if self.legend:
                newrow(layout, 'Scale:', self, 'nscale')
                if self.nscale == '1':
                    newrow(layout, 'Decades:', self, 'decades')
                newrow(layout, 'Legend max:', self, 'lmax')
                newrow(layout, 'Legend width:', self, 'lw')
                newrow(layout, 'Legend height:', self, 'lh')
            
            newrow(layout, 'Contour:', self, 'contour')
            
            if self.contour:
               newrow(layout, 'Overlay:', self, 'overlay') 
               if self.overlay:
                   newrow(layout, 'Overlay file:', self, 'ofile') 
                   newrow(layout, 'Overlay exposure:', self, 'disp')
               newrow(layout, 'Bands:', self, 'bands') 
   
            if self.inputs['Image'].links and os.path.isfile(self.inputs['Image'].links[0].from_node['images'][0]):
                row = layout.row()
                row.operator("node.livifc", text = 'Process')
            
    def presim(self):
        self['basename'] = self.basename if self.basename else 'fc'
        
    def postsim(self):
        self['exportstate'] = [str(x) for x in (self.basename, self.colour, self.lmax, self.unit, self.nscale, self.decades, 
                   self.legend, self.lw, self.lh, self.contour, self.overlay, self.bands)]
        nodecolour(self, 0)
        

        
class No_Li_Sim(Node, ViNodes):
    '''Node describing a LiVi simulation'''
    bl_idname = 'No_Li_Sim'
    bl_label = 'LiVi Simulation'

    def nodeupdate(self, context):
        nodecolour(self, self['exportstate'] != [str(x) for x in (self.cusacc, self.simacc, self.csimacc, self.pmap, self.pmapcno, self.pmapgno)])
        if self.simacc == '3':
            self.validparams = validradparams(self.cusacc)
        
    simacc: EnumProperty(items=[("0", "Low", "Low accuracy and high speed (preview)"),("1", "Medium", "Medium speed and accuracy"), ("2", "High", "High but slow accuracy"),("3", "Custom", "Edit Radiance parameters"), ],
            name="", description="Simulation accuracy", default="0", update = nodeupdate)
    csimacc: EnumProperty(items=[("0", "Custom", "Edit Radiance parameters"), ("1", "Initial", "Initial accuracy for this metric"), ("2", "Final", "Final accuracy for this metric")],
            name="", description="Simulation accuracy", default="1", update = nodeupdate)
    cusacc: StringProperty(
            name="", description="Custom Radiance simulation parameters", default="", update = nodeupdate)
    rtracebasic = (("-ab", 2, 3, 4), ("-ad", 256, 1024, 4096), ("-as", 128, 512, 2048), ("-aa", 0, 0, 0), ("-dj", 0, 0.7, 1), ("-ds", 0, 0.5, 0.15), ("-dr", 1, 3, 5), ("-ss", 0, 2, 5), ("-st", 1, 0.75, 0.1), ("-lw", 0.0001, 0.00001, 0.000002), ("-lr", 2, 3, 4))
    rtraceadvance = (("-ab", 3, 5), ("-ad", 4096, 8192), ("-as", 512, 1024), ("-aa", 0.0, 0.0), ("-dj", 0.7, 1), ("-ds", 0.5, 0.15), ("-dr", 2, 3), ("-ss", 2, 5), ("-st", 0.75, 0.1), ("-lw", 1e-4, 1e-5), ("-lr", 3, 5))
    rvubasic = (("-ab", 2, 3, 4), ("-ad", 256, 1024, 4096), ("-as", 128, 512, 2048), ("-aa", 0, 0, 0), ("-dj", 0, 0.7, 1), ("-ds", 0.5, 0.15, 0.15), ("-dr", 1, 3, 5), ("-ss", 0, 2, 5), ("-st", 1, 0.75, 0.1), ("-lw", 0.0001, 0.00001, 0.0000002), ("-lr", 3, 3, 4))
    rvuadvance = (("-ab", 3, 5), ("-ad", 4096, 8192), ("-as", 1024, 2048), ("-aa", 0.0, 0.0), ("-dj", 0.7, 1), ("-ds", 0.5, 0.15), ("-dr", 2, 3), ("-ss", 2, 5), ("-st", 0.75, 0.1), ("-lw", 1e-4, 1e-5), ("-lr", 3, 5))
    pmap: BoolProperty(name = '', default = False)
    pmapgno: IntProperty(name = '', default = 50000)
    pmapcno: IntProperty(name = '', default = 0)
    run: IntProperty(default = 0)
    validparams: BoolProperty(name = '', default = True)
    illu: BoolProperty(name = '', default = False)

    def init(self, context):
        self['simdict'] = {'Basic': 'simacc', 'Compliance':'csimacc', 'CBDM':'csimacc'}
        self.inputs.new('So_Li_Geo', 'Geometry in')
        self.inputs.new('So_Li_Con', 'Context in')
        self.outputs.new('So_Vi_Res', 'Results out')
        nodecolour(self, 1)
        self['maxres'], self['minres'], self['avres'], self['exportstate'], self['year'] = {}, {}, {}, '', 2015
        
    def draw_buttons(self, context, layout): 
        scene = context.scene
        svp = scene.vi_params
        try:
            row = layout.row()
            row.label(text = 'Frames: {} - {}'.format(min([c['fs'] for c in (self.inputs['Context in'].links[0].from_node['Options'], self.inputs['Geometry in'].links[0].from_node['Options'])]), max([c['fe'] for c in (self.inputs['Context in'].links[0].from_node['Options'], self.inputs['Geometry in'].links[0].from_node['Options'])])))
            cinnode = self.inputs['Context in'].links[0].from_node
            newrow(layout, 'Photon map:', self, 'pmap')
            if self.pmap:
               newrow(layout, 'Global photons:', self, 'pmapgno')
               newrow(layout, 'Caustic photons:', self, 'pmapcno')
            row = layout.row()
            row.label(text = "Accuracy:")            
            row.prop(self, self['simdict'][cinnode['Options']['Context']])
            
            if (self.simacc == '3' and cinnode['Options']['Context'] == 'Basic') or (self.csimacc == '0' and cinnode['Options']['Context'] in ('Compliance', 'CBDM')):
               newrow(layout, "Radiance parameters:", self, 'cusacc')
            if not self.run and self.validparams:
                if cinnode['Options']['Preview']:
                    row = layout.row()
                    row.operator("node.radpreview", text = 'Preview')
#                if cinnode['Options']['Context'] == 'Basic' and cinnode['Options']['Type'] == '1' and not self.run:
#                    row.operator("node.liviglare", text = 'Calculate').nodeid = self['nodeid']
                if [o for o in scene.objects if o.name in svp['liparams']['livic']]:
                    row.operator("node.livicalc", text = 'Calculate')
        except Exception as e:
            logentry('Problem with LiVi simulation: {}'.format(e))

    def update(self):
        if self.outputs.get('Results out'):
            socklink(self.outputs['Results out'], self.id_data.name)
        self.run = 0
    
    def presim(self):
        self['coptions'] = self.inputs['Context in'].links[0].from_node['Options']
        self['goptions'] = self.inputs['Geometry in'].links[0].from_node['Options']
        self['radfiles'], self['reslists'] = {}, [[]]
        if self['coptions']['Context'] == 'Basic':
            self['radparams'] = self.cusacc if self.simacc == '3' else (" {0[0]} {1[0]} {0[1]} {1[1]} {0[2]} {1[2]} {0[3]} {1[3]} {0[4]} {1[4]} {0[5]} {1[5]} {0[6]} {1[6]} {0[7]} {1[7]} {0[8]} {1[8]} {0[9]} {1[9]} {0[10]} {1[10]} ".format([n[0] for n in self.rtracebasic], [n[int(self.simacc)+1] for n in self.rtracebasic]))
        else:
            self['radparams'] = self.cusacc if self.csimacc == '0' else (" {0[0]} {1[0]} {0[1]} {1[1]} {0[2]} {1[2]} {0[3]} {1[3]} {0[4]} {1[4]} {0[5]} {1[5]} {0[6]} {1[6]} {0[7]} {1[7]} {0[8]} {1[8]} {0[9]} {1[9]} {0[10]} {1[10]} ".format([n[0] for n in self.rtraceadvance], [n[int(self.csimacc)] for n in self.rtraceadvance]))
    
    def sim(self, scene):
        svp = scene.vi_params
        self['frames'] = range(svp['liparams']['fs'], svp['liparams']['fe'] + 1)
        
    def postsim(self):
        self['exportstate'] = [str(x) for x in (self.cusacc, self.simacc, self.csimacc, self.pmap, self.pmapcno, self.pmapgno)]
        nodecolour(self, 0)
        
class ViSPNode(Node, ViNodes):
    '''Node describing a VI-Suite sun path'''
    bl_idname = 'ViSPNode'
    bl_label = 'VI Sun Path'
    
    def nodeupdate(self, context):
        nodecolour(self, self['exportstate'] != [str(x) for x in (self.suns)])
    
    suns: EnumProperty(items = [('0', 'Single', 'Single sun'), ('1', 'Monthly', 'Monthly sun for chosen time'), ('2', 'Hourly', 'Hourly sun for chosen date')], name = '', description = 'Sunpath sun type', default = '0', update=nodeupdate)
#    res: FloatProperty(name="", description="Calc point offset", min=1, max=10, default=6, update = nodeupdate)

    def init(self, context):
#        self['nodeid'] = nodeid(self)
        self.inputs.new('So_Vi_Loc', 'Location in')
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

class ViWRNode(Node, ViNodes):
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
        self.inputs.new('So_Vi_Loc', 'Location in')
        self['exportstate'] = ''
        nodecolour(self, 1)

    def draw_buttons(self, context, layout):
        if nodeinputs(self) and self.inputs[0].links[0].from_node.loc == '1':
            (sdate, edate) = retdates(self.sdoy, self.edoy, self.inputs[0].links[0].from_node['year'])
            newrow(layout, 'Type:', self, "wrtype")
            newrow(layout, 'Start day {}/{}:'.format(sdate.day, sdate.month), self, "sdoy")
            newrow(layout, 'End day {}/{}:'.format(edate.day, edate.month), self, "edoy")
            newrow(layout, 'Colour:', context.scene.vi_params, 'vi_scatt_col')
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
        self['goptions'] = {}
        self.outputs.new('So_Vi_Res', 'Results out')
        nodecolour(self, 1)

    def draw_buttons(self, context, layout):
        newrow(layout, 'Ignore sensor:', self, "signore")
        newrow(layout, 'Animation:', self, "animmenu")
        if self.animmenu != 'Static': 
            row = layout.row(align=True)
            row.alignment = 'EXPAND'
            row.label(text = 'Frames:')
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
        
class ViSSNode(Node, ViNodes):
    '''Node to create a VI-Suite shadow map'''
    bl_idname = 'ViSSNode'
    bl_label = 'VI Shadow Map'

    def nodeupdate(self, context):
        nodecolour(self, self['exportstate'] != [str(x) for x in (self.animmenu, self.sdoy, self.edoy, self.starthour, self.endhour, self.interval, self.cpoint, self.offset)])

    animtype = [('Static', "Static", "Simple static analysis"), ('Geometry', "Geometry", "Animated geometry analysis")]
    animmenu: EnumProperty(name="", description="Animation type", items=animtype, default = 'Static', update = nodeupdate)
    startframe: IntProperty(name = '', default = 0, min = 0, max = 1024, description = 'Start frame')
    endframe: IntProperty(name = '', default = 0, min = 0, max = 1024, description = 'End frame')
    starthour: IntProperty(name = '', default = 1, min = 1, max = 24, description = 'Start hour')
    endhour: IntProperty(name = '', default = 24, min = 1, max = 24, description = 'End hour')
    interval: IntProperty(name = '', default = 1, min = 1, max = 60, description = 'Number of simulation steps per hour')
    sdoy: IntProperty(name = '', default = 1, min = 1, max = 365, description = 'Start Day', update = nodeupdate)
    edoy: IntProperty(name = '', default = 365, min = 1, max = 365, description = 'End Day', update = nodeupdate)
    cpoint: EnumProperty(items=[("0", "Faces", "Export faces for calculation points"),("1", "Vertices", "Export vertices for calculation points"), ],
            name="", description="Specify the calculation point geometry", default="0", update = nodeupdate)
    offset: FloatProperty(name="", description="Calc point offset", min=0.001, max=10, default=0.01, update = nodeupdate)
    signore: BoolProperty(name = '', default = 0, description = 'Ignore sensor surfaces', update = nodeupdate)
    
    def init(self, context):
#        self['nodeid'] = nodeid(self)
        self.inputs.new('So_Vi_Loc', 'Location in')
        self.outputs.new('ViR', 'Results out')
        self.outputs['Results out'].hide = True
        self['exportstate'] = ''
        self['goptions'] = {}
        nodecolour(self, 1)

    def draw_buttons(self, context, layout):
        if nodeinputs(self):
            (sdate, edate) = retdates(self.sdoy, self.edoy, self.inputs[0].links[0].from_node['year'])
            newrow(layout, 'Ignore sensor:', self, "signore")
            newrow(layout, 'Animation:', self, "animmenu")
            if self.animmenu != 'Static':            
                row = layout.row(align=True)
                row.alignment = 'EXPAND'
                row.label(text = 'Frames:')
                row.prop(self, 'startframe')
                row.prop(self, 'endframe')
            newrow(layout, 'Start day {}/{}:'.format(sdate.day, sdate.month), self, "sdoy")
            newrow(layout, 'End day {}/{}:'.format(edate.day, edate.month), self, "edoy")
            newrow(layout, 'Start hour:', self, "starthour")
            newrow(layout, 'End hour:', self, "endhour")
            newrow(layout, 'Hour steps:', self, "interval")
            newrow(layout, 'Result point:', self, "cpoint")
            newrow(layout, 'Offset:', self, 'offset')
            row = layout.row()
            row.operator("node.shad", text = 'Calculate')

    def preexport(self):
        (self.sdate, self.edate) = retdates(self.sdoy, self.edoy, self.inputs[0].links[0].from_node['year'])
        self['goptions']['offset'] = self.offset

    def postexport(self, scene):
        nodecolour(self, 0)
        self.outputs['Results out'].hide = False if self.get('reslists') else True            
        self['exportstate'] = [str(x) for x in (self.animmenu, self.sdoy, self.edoy, self.starthour, self.endhour, self.interval, self.cpoint, self.offset)]
    
    def update(self):
        if self.outputs.get('Results out'):
            socklink(self.outputs['Results out'], self.id_data.name)

# Edit nodes
class No_Text(Node, ViNodes):
    '''Text Export Node'''
    bl_idname = 'No_Text'
    bl_label = 'VI Text Edit'
    
    contextmenu: StringProperty(name = '')

    def init(self, context):
        self['bt'] = ''
        self.outputs.new('So_Text', 'Text out')
        self.inputs.new('So_Text', 'Text in')
        self.outputs['Text out']['Text'] = {}
        self.outputs['Text out']['Options'] = {}
        
    def draw_buttons(self, context, layout):
        if self.inputs['Text in'].links:
            inodename = self.inputs['Text in'].links[0].from_node.name
            row = layout.row()
            row.label(text = 'Text name: {}'.format(inodename))            
            if inodename in [im.name for im in bpy.data.texts] and self['bt'] != bpy.data.texts[inodename].as_string():
                row = layout.row()
                row.operator('node.textupdate', text = 'Update')

    def update(self):
        socklink(self.outputs['Text out'], self.id_data.name)
        if self.inputs and self.inputs['Text in'].links:
            self['Options'] = self.inputs['Text in'].links[0].from_node['Options']
            self['Text'] = self.inputs['Text in'].links[0].from_node['Text']
            inodename = self.inputs['Text in'].links[0].from_node.name
            sframes = sorted([int(frame) for frame in self.inputs['Text in'].links[0].from_node['Text'].keys()])
            t = ''.join(['# Frame {}\n{}\n\n'.format(f, self.inputs['Text in'].links[0].from_node['Text'][str(f)]) for f in sframes])
            bt = bpy.data.texts.new(inodename) if inodename not in [im.name for im in bpy.data.texts] else bpy.data.texts[inodename]
            bt.from_string(t)
            self['bt'] = bt.as_string()
        else:
            self['Text'] = {}

    def textupdate(self, bt):
        inodename = self.inputs['Text in'].links[0].from_node.name
        bt = bpy.data.texts.new(inodename) if inodename not in [im.name for im in bpy.data.texts] else bpy.data.texts[inodename]
        btlines = [line.body for line in bt.lines]
        self['bt'] = bt.as_string()
        btheads = [line for line in btlines if '# Frame' in line]
        btstring = ''.join([self['bt'].replace(bth, '***') for bth in btheads])
        btbodies = btstring.split('***\n')[1:]
        btframes = [head.split()[2] for head in btheads]
        self['Text'] = {bthb[0]:bthb[1] for bthb in zip(btframes, btbodies)}
        
class No_En_Geo(Node, ViNodes):
    '''Node describing an EnVi Geometry Export'''
    bl_idname = 'No_En_Geo'
    bl_label = 'EnVi Geometry'
    
    geo_offset = FloatVectorProperty(name="", description="", default=(0.0, 0.0, 0.0), min=sys.float_info.min, max=sys.float_info.max, soft_min=sys.float_info.min, soft_max=sys.float_info.max, step=3, precision=1, subtype='TRANSLATION', unit='NONE', size=3, update=None, get=None, set=None)
#    def nodeupdate(self, context):
#        nodecolour(self, self['exportstate'] != [str(x) for x in (self.animmenu)])
    
    def init(self, context):
        self.outputs.new('So_En_Geo', 'Geometry out')
#        self['exportstate'] = ''
        nodecolour(self, 1)

    def draw_buttons(self, context, layout):
        newrow(layout, 'Offset', self, 'geo_offset')
        row = layout.row()
        row.operator("node.engexport", text = "Export")

    def update(self):
        socklink(self.outputs['Geometry out'], self.id_data.name)
        
    def preexport(self, scene):
         pass
               
    def postexport(self):
        nodecolour(self, 0)
        
class No_En_Con(Node, ViNodes):
    '''Node describing an EnergyPlus export'''
    bl_idname = 'No_En_Con'
    bl_label = 'EnVi Export'

    def nodeupdate(self, context):
        nodecolour(self, self['exportstate'] != [str(x) for x in (self.loc, self.terrain, self.timesteps, self.animated, self.fs, self.fe, self.sdoy, self.edoy)])

    animated: BoolProperty(name="", description="Animated analysis", update = nodeupdate)
    fs: IntProperty(name="", description="Start frame", default = 0, min = 0, update = nodeupdate)
    fe: IntProperty(name="", description="End frame", default = 0, min = 0, update = nodeupdate)
    loc: StringProperty(name="", description="Identifier for this project", default="", update = nodeupdate)
    terrain: EnumProperty(items=[("0", "City", "Towns, city outskirts, centre of large cities"),
                   ("1", "Urban", "Urban, Industrial, Forest"),("2", "Suburbs", "Rough, Wooded Country, Suburbs"),
                    ("3", "Country", "Flat, Open Country"),("4", "Ocean", "Ocean, very flat country")],
                    name="", description="Specify the surrounding terrain", default="0", update = nodeupdate)

    addonpath = os.path.dirname(inspect.getfile(inspect.currentframe()))
    matpath = addonpath+'/EPFiles/Materials/Materials.data'
    sdoy: IntProperty(name = "", description = "Day of simulation", min = 1, max = 365, default = 1, update = nodeupdate)
    edoy: IntProperty(name = "", description = "Day of simulation", min = 1, max = 365, default = 365, update = nodeupdate)
    timesteps: IntProperty(name = "", description = "Time steps per hour", min = 1, max = 60, default = 1, update = nodeupdate)
    restype: EnumProperty(items = [("0", "Zone Thermal", "Thermal Results"), ("1", "Comfort", "Comfort Results"), 
                                   ("2", "Zone Ventilation", "Zone Ventilation Results"), ("3", "Ventilation Link", "Ventilation Link Results"), 
                                   ("4", "Thermal Chimney", "Thermal Chimney Results"), ("5", "Power", "Power Production Results")],
                                   name="", description="Specify the EnVi results category", default="0", update = nodeupdate)

    (resaam, resaws, resawd, resah, resasm, restt, resh, restwh, restwc, reswsg, rescpp, rescpm, resvls, resvmh, resim, resiach, resco2, resihl, resl12ms,
     reslof, resmrt, resocc, resh, resfhb, ressah, ressac, reshrhw, restcvf, restcmf, restcot, restchl, restchg, restcv, restcm, resldp, resoeg,
     respve, respvw, respvt, respveff) = resnameunits()
     
    def init(self, context):
        self.inputs.new('So_En_Geo', 'Geometry in')
        self.inputs.new('So_Vi_Loc', 'Location in')
        self.outputs.new('So_En_Con', 'Context out')
        self['exportstate'] = ''
        self['year'] = 2015
        nodecolour(self, 1)

    def draw_buttons(self, context, layout):
        (sdate, edate) = retdates(self.sdoy, self.edoy, self['year'])
        row = layout.row()
        row.label(text = 'Animation:')
        row.prop(self, 'animated')
        if self.animated:
            newrow(layout, 'Start frame:', self, 'fs')
            newrow(layout, 'End frame:', self, 'fe')
        newrow(layout, "Name/location", self, "loc")
        row = layout.row()
        row.label(text = 'Terrain:')
        col = row.column()
        col.prop(self, "terrain")
        newrow(layout, 'Start day {}/{}:'.format(sdate.day, sdate.month), self, "sdoy")
        newrow(layout, 'End day {}/{}:'.format(edate.day, edate.month), self, "edoy")
        newrow(layout, 'Time-steps/hour', self, "timesteps")
        row = layout.row()
        row.label(text = 'Results Category:')
        col = row.column()
        col.prop(self, "restype")
        resdict = enresprops('')
        
        for rprop in resdict[self.restype]:
            if not rprop:
                row = layout.row()
            else:
                row.prop(self, rprop)
                
        if all([s.links for s in self.inputs]) and not any([s.links[0].from_node.use_custom_color for s in self.inputs]):
            row = layout.row()
            row.operator("node.encon", text = 'Export')

    def update(self):
        if self.inputs.get('Location in') and self.outputs.get('Context out'):
            socklink(self.outputs['Context out'], self.id_data.name)
            self['year'] = self.inputs['Location in'].links[0].from_node['year'] if self.inputs['Location in'].links else 2015
    
    def preexport(self, scene):
        (self.fs, self.fe) = (self.fs, self.fe) if self.animated else (scene.frame_current, scene.frame_current)
        scene.vi_params['enparams']['fs'], scene.vi_params['enparams']['fe'] = self.fs, self.fe
        (self.sdate, self.edate) = retdates(self.sdoy, self.edoy, self['year'])
        
    def postexport(self):
        nodecolour(self, 0)
        self['exportstate'] = [str(x) for x in (self.loc, self.terrain, self.timesteps, self.animated, self.fs, self.fe, self.sdoy, self.edoy)]

class No_En_Sim(Node, ViNodes):
    '''Node describing an EnergyPlus simulation'''
    bl_idname = 'No_En_Sim'
    bl_label = 'EnVi Simulation'

    def init(self, context):
        self.inputs.new('So_En_Con', 'Context in')
        self.outputs.new('So_Vi_Res', 'Results out')
        self['exportstate'] = ''
        self['Start'], self['End'] = 1, 365
        self['AStart'], self['AEnd'] = 0, 0
        nodecolour(self, 1)

    def nodeupdate(self, context):
        nodecolour(self, self['exportstate'] != [self.resname])

    resname: StringProperty(name="", description="Base name for the results files", default="results", update = nodeupdate)
    resfilename: StringProperty(name = "", default = 'results')
    dsdoy: IntProperty()
    dedoy: IntProperty() 
    run: IntProperty(min = -1, default = -1)
    processors: IntProperty(name = '', min = 1, default = 4)#max = bpy.context.scene['viparams']['nproc'], default = bpy.context.scene['viparams']['nproc'])
    mp: BoolProperty(name = "", default = False)

    def draw_buttons(self, context, layout):
        scene = context.scene
        svp = scene.vi_params
        if self.run > -1:
            row = layout.row()
            row.label(text = 'Calculating {}%'.format(self.run))
        elif self.inputs['Context in'].links and not self.inputs['Context in'].links[0].from_node.use_custom_color:
            if svp['enparams']['fe'] > svp['enparams']['fs']:
                newrow(layout, 'Multi-core:', self, 'mp')
                if self.mp:
                    newrow(layout, 'Processors:', self, 'processors')
            newrow(layout, 'Results name:', self, 'resname')
            row = layout.row()
            row.operator("node.ensim", text = 'Calculate')

    def update(self):
        if self.outputs.get('Results out'):
            socklink(self.outputs['Results out'], self.id_data.name)

    def presim(self, context):
        innode = self.inputs['Context in'].links[0].from_node
        self['frames'] = range(context.scene.vi_params['enparams']['fs'], context.scene.vi_params['enparams']['fe'] + 1)
        self.resfilename = os.path.join(context.scene.vi_params['viparams']['newdir'], self.resname+'.eso')
        self['year'] = innode['year']
        self.dsdoy = innode.sdoy # (locnode.startmonthnode.sdoy
        self.dedoy = innode.edoy
        self["_RNA_UI"] = {"Start": {"min":innode.sdoy, "max":innode.edoy}, "End": {"min":innode.sdoy, "max":innode.edoy}, 
            "AStart": {"name": '', "min": context.scene.vi_params['enparams']['fs'], "max": context.scene.vi_params['enparams']['fe']}, 
            "AEnd": {"min":context.scene.vi_params['enparams']['fs'], "max":context.scene.vi_params['enparams']['fe']}}
        self['Start'], self['End'] = innode.sdoy, innode.edoy
#        self["_RNA_UI"] = {"AStart": {"min":context.scene['enparams']['fs'], "max":context.scene['enparams']['fe']}, "AEnd": {"min":context.scene['enparams']['fs'], "max":context.scene['enparams']['fe']}}
        self['AStart'], self['AEnd'] = context.scene.vi_params['enparams']['fs'], context.scene.vi_params['enparams']['fe']
     
    def postsim(self, sim_op, condition):
#        scene = bpy.context.scene
        nodecolour(self, 0)
        self.run = -1
        if condition == 'FINISHED':
            processf(sim_op, self)

class No_Vi_Chart(Node, ViNodes):
    '''Node for 2D results plotting'''
    bl_idname = 'No_Vi_Chart'
    bl_label = 'VI Chart'
    
    def aupdate(self, context):
        self.update()

    def pmitems(self, context):
        return [tuple(p) for p in self['pmitems']]

    ctypes = [("0", "Line/Scatter", "Line/Scatter Plot")]
    charttype: EnumProperty(items = ctypes, name = "Chart Type", default = "0")
    timemenu: EnumProperty(items=[("0", "Hourly", "Hourly results"),("1", "Daily", "Daily results"), ("2", "Monthly", "Monthly results")],
                name="Period", description="Results frequency", default="0")
    parametricmenu: EnumProperty(items=pmitems, name="", description="Parametric result display", update=aupdate)    
    bl_width_max = 800
    dpi: IntProperty(name = 'DPI', description = "DPI of the shown figure", default = 92, min = 92)
           
    def init(self, context):
        self.inputs.new("ViEnRXIn", "X-axis")
        self.inputs.new("ViEnRY1In", "Y-axis 1")
        self.inputs["Y-axis 1"].hide = True
        self.inputs.new("ViEnRY2In", "Y-axis 2")
        self.inputs["Y-axis 2"].hide = True
        self.inputs.new("ViEnRY3In", "Y-axis 3")
        self.inputs["Y-axis 3"].hide = True
        self['Start'], self['End'] = 1, 365
        self['pmitems'] = [("0", "Static", "Static results")]
        self.update()

    def draw_buttons(self, context, layout):
        if self.inputs['X-axis'].links:
            innode = self.inputs['X-axis'].links[0].from_node
            if innode.get('reslists'):
                newrow(layout, 'Animated:', self, 'parametricmenu')
                if self.parametricmenu == '0':                
                    (sdate, edate) = retdates(self['Start'], self['End'], innode['year']) 
                    label = "Start/End Day: {}/{} {}/{}".format(sdate.day, sdate.month, edate.day, edate.month)
                else:
                    row = layout.row()
                    label = "Frame"
     
                row = layout.row()    
                row.label(text = label)
                row.prop(self, '["Start"]')
                row.prop(self, '["End"]')
                    
                if self.parametricmenu == '0':
                    row = layout.row()
                    row.prop(self, "charttype")
                    row.prop(self, "timemenu")
                    row.prop(self, "dpi")
    
                if self.inputs['Y-axis 1'].links and 'NodeSocketUndefined' not in [sock.bl_idname for sock in self.inputs if sock.links]:
                    row = layout.row()
                    row.operator("node.chart", text = 'Create plot')
                    row = layout.row()
                    row.label(text = "------------------")

    def update(self):
        try:
            if not self.inputs['X-axis'].links or not self.inputs['X-axis'].links[0].from_node['reslists']:
                class ViEnRXIn(So_En_ResU):
                    '''Energy geometry out socket'''
                    bl_idname = 'ViEnRXIn'
                    bl_label = 'X-axis'    
                    valid = ['Vi Results']
                    
            else:
                innode = self.inputs['X-axis'].links[0].from_node
                rl = innode['reslists']
                zrl = list(zip(*rl))
                try:
                    if len(set(zrl[0])) > 1:
                        self['pmitems'] = [("0", "Static", "Static results"), ("1", "Parametric", "Parametric results")]
                    else:
                        self['pmitems'] = [("0", "Static", "Static results")]
                except:
                    self['pmitems'] = [("0", "Static", "Static results")]
                
                time.sleep(0.1)
    
                if self.parametricmenu == '1' and len(set(zrl[0])) > 1:
                    frames = [int(k) for k in set(zrl[0]) if k != 'All']
                    startframe, endframe = min(frames), max(frames)
                    self["_RNA_UI"] = {"Start": {"min":startframe, "max":endframe}, "End": {"min":startframe, "max":endframe}}
                    self['Start'], self['End'] = startframe, endframe
                else:
                    if 'Month' in zrl[3]:
                        startday = datetime.datetime(int(innode['year']), int(zrl[4][zrl[3].index('Month')].split()[0]), int(zrl[4][zrl[3].index('Day')].split()[0])).timetuple().tm_yday
                        endday = datetime.datetime(int(innode['year']), int(zrl[4][zrl[3].index('Month')].split()[-1]), int(zrl[4][zrl[3].index('Day')].split()[-1])).timetuple().tm_yday
                        self["_RNA_UI"] = {"Start": {"min":startday, "max":endday}, "End": {"min":startday, "max":endday}}
                        self['Start'], self['End'] = startday, endday
    
                if self.inputs.get('Y-axis 1'):
                    self.inputs['Y-axis 1'].hide = False
    
                class ViEnRXIn(So_En_Res):
                    '''Energy geometry out socket'''
                    bl_idname = 'ViEnRXIn'
                    bl_label = 'X-axis'
                                    
#                    if innode['reslists']:
                    (valid, framemenu, statmenu, rtypemenu, climmenu, zonemenu, 
                     zonermenu, linkmenu, linkrmenu, enmenu, enrmenu, chimmenu, 
                     chimrmenu, posmenu, posrmenu, cammenu, camrmenu, powmenu, powrmenu, multfactor) = retrmenus(innode, self, 'X-axis')
                        
            bpy.utils.register_class(ViEnRXIn)
    
            if self.inputs.get('Y-axis 1'):
                if not self.inputs['Y-axis 1'].links or not self.inputs['Y-axis 1'].links[0].from_node['reslists']:
                    class ViEnRY1In(So_En_ResU):
                        '''Energy geometry out socket'''
                        bl_idname = 'ViEnRY1In'
                        bl_label = 'Y-axis 1'
    
                    if self.inputs.get('Y-axis 2'):
                        self.inputs['Y-axis 2'].hide = True
                else:
                    innode = self.inputs['Y-axis 1'].links[0].from_node
    
                    class ViEnRY1In(So_En_Res):
                        '''Energy geometry out socket'''
                        bl_idname = 'ViEnRY1In'
                        bl_label = 'Y-axis 1'
                        (valid, framemenu, statmenu, rtypemenu, climmenu, zonemenu, 
                         zonermenu, linkmenu, linkrmenu, enmenu, enrmenu, chimmenu, 
                         chimrmenu, posmenu, posrmenu, cammenu, camrmenu, powmenu, powrmenu, multfactor) = retrmenus(innode, self, 'Y-axis 1')
    
                    self.inputs['Y-axis 2'].hide = False
                bpy.utils.register_class(ViEnRY1In)
    
            if self.inputs.get('Y-axis 2'):
                if not self.inputs['Y-axis 2'].links or not self.inputs['Y-axis 2'].links[0].from_node['reslists']:
                    class ViEnRY2In(So_En_ResU):
                        '''Energy geometry out socket'''
                        bl_idname = 'ViEnRY2In'
                        bl_label = 'Y-axis 2'
    
                    if self.inputs.get('Y-axis 3'):
                        self.inputs['Y-axis 3'].hide = True
                else:
                    innode = self.inputs[2].links[0].from_node
    
                    class ViEnRY2In(So_En_Res):
                        '''Energy geometry out socket'''
                        bl_idname = 'ViEnRY2In'
                        bl_label = 'Y-axis 2'
    
                        (valid, framemenu, statmenu, rtypemenu, climmenu, zonemenu, 
                         zonermenu, linkmenu, linkrmenu, enmenu, enrmenu, chimmenu, 
                         chimrmenu, posmenu, posrmenu, cammenu, camrmenu, powmenu, powrmenu, multfactor) = retrmenus(innode, self, 'Y-axis 2')
    
                    self.inputs['Y-axis 3'].hide = False
    
                bpy.utils.register_class(ViEnRY2In)
    
            if self.inputs.get('Y-axis 3'):
                if not self.inputs['Y-axis 3'].links or not self.inputs['Y-axis 3'].links[0].from_node['reslists']:
                    class ViEnRY3In(So_En_ResU):
                        '''Energy geometry out socket'''
                        bl_idname = 'ViEnRY3In'
                        bl_label = 'Y-axis 3'
                else:
                    innode = self.inputs[3].links[0].from_node
    
                    class ViEnRY3In(So_En_Res):
                        '''Energy geometry out socket'''
                        bl_idname = 'ViEnRY3In'
                        bl_label = 'Y-axis 3'
    
                        (valid, framemenu, statmenu, rtypemenu, climmenu, zonemenu, 
                         zonermenu, linkmenu, linkrmenu, enmenu, enrmenu, chimmenu, 
                         chimrmenu, posmenu, posrmenu, cammenu, camrmenu, powmenu, powrmenu, multfactor) = retrmenus(innode, self, 'Y-axis 3')
    
                bpy.utils.register_class(ViEnRY3In)
        except Exception as e:
            print('Chart node update failure 2 ', e)

class No_Vi_Metrics(Node, ViNodes):
    '''Node for result metrics'''
    bl_idname = 'No_Vi_Metrics'
    bl_label = 'VI Metrics' 
    
    def zupdate(self, context):
        self.update()
        
    def zitems(self, context):
        if self.inputs[0].links:
            rl = self.inputs[0].links[0].from_node['reslists']
#            zrl = list(zip(*rl))
            znames = set([z[2] for z in rl if z[1] == 'Zone'])
            return [(zn, zn, 'Zone name') for zn in znames] + [('All', 'All', 'All zones')]
        else:
            return [('None', 'None', 'None')]
    
    metric: EnumProperty(items=[("0", "Energy", "Energy results"), ("1", "Lighting", "Lighting results")],
                name="", description="Results type", default="0", update=zupdate)   
    energy_menu: EnumProperty(items=[("0", "SAP", "SAP results")],
                name="", description="Results metric", default="0", update=zupdate)
    light_menu: EnumProperty(items=[("0", "BREEAM", "BREEAM HEA1 results"),
                                    ("1", "LEED", "LEED v4 results")],
                name="", description="Results metric", default="0", update=zupdate)
    zone_menu: EnumProperty(items=zitems,
                name="", description="Zone results", update=zupdate)
    mod: FloatProperty(name="kWh", description="Energy modifier (kWh)", update=zupdate)
    
    def init(self, context):
        self['res'] = {}
        self.inputs.new('So_Vi_Res', 'Results in') 
        
    def draw_buttons(self, context, layout):
        newrow(layout, 'Type:', self, "metric")
        if self.metric == '0':
            newrow(layout, 'Metric:', self, "energy_menu")
        else:
            newrow(layout, 'Metric:', self, "light_menu")
        newrow(layout, 'Zone', self, "zone_menu")
        
        if self.metric == '0':
            if self.energy_menu == '0':
                if self['res'] and self['res'].get('hkwh'):
                    newrow(layout, 'Modifier', self, 'mod')
                    row = layout.row()
                    pvkwh = self['res']['pvkwh'] if self['res']['pvkwh'] == 'N/A' else "{:.2f}".format(self['res']['pvkwh'])
                    row.label(text = "PV (kWh): {}".format(pvkwh))
                    pva = "{:.2f}".format(self['res']['pvkwh']/self['res']['fa']) if self['res']['fa'] != 'N/A' and self['res']['fa'] > 0 else 'N/A' 
                    row = layout.row()
                    row.label(text = "PV (kWh/m2): {}".format(pva))
                    row = layout.row()
                    hkwh = self['res']['hkwh'] if self['res']['hkwh'] == 'N/A' else "{:.2f}".format(self['res']['hkwh'] + self['res']['ahkwh'])
                    row.label(text = "Heating (kWh): {}".format(hkwh))
                    
                    row = layout.row()
                    ha = "{:.2f}".format((self['res']['hkwh'] + self['res']['ahkwh'])/self['res']['fa']) if self['res']['fa'] != 'N/A' and self['res']['fa'] > 0 else 'N/A' 
                    row.label(text = "Heating (kWh/m2): {}".format(ha))
                    row = layout.row()
                    ckwh = self['res']['pvkwh'] if self['res']['pvkwh'] == 'N/A' else "{:.2f}".format(self['res']['ckwh'])
                    row.label(text = "Cooling (kWh): {}".format(ckwh))
                    row = layout.row()
                    ca = "{:.2f}".format(self['res']['ckwh']/self['res']['fa']) if self['res']['fa'] != 'N/A' and self['res']['fa'] > 0 else 'N/A' 
                    row.label(text = "Cooling (kWh/m2): {}".format(ca))
                    
                    if self.zone_menu == 'All':
                        row = layout.row()
                        wkwh = self['res']['wkwh'] if self['res']['wkwh'] == 'N/A' else "{:.2f}".format(self['res']['wkwh'])
                        row.label(text = "Hot water (kWh): {}".format(wkwh))
                        row = layout.row()
                        ecf = "{:.2f}".format(self['res']['ECF']) if self['res']['ECF'] != 'N/A' else 'N/A' 
                        row.label(text = "ECF: {}".format(ecf))
                        row = layout.row()
                        epc = "{:.0f}".format(self['res']['EPC']) if self['res']['EPC'] != 'N/A' else 'N/A' 
                        row.label(text = "EPC: {} ({})".format(epc, self['res']['EPCL']))
                    
    def update(self):
        if self.inputs[0].links:
            if self.metric == '0' and bpy.data.collections.get('EnVi Geometry'):
                rl = self.inputs[0].links[0].from_node['reslists']
                self['res']['pvkwh'] = 0
                self['res']['hkwh'] = 0
                self['res']['ahkwh'] = 0
                self['res']['ckwh'] = 0            
                self['res']['fa'] = sum([c.vi_params['enparams']['floorarea'] for c in bpy.data.collections['EnVi Geometry'].children]) if self.zone_menu == 'All' else bpy.data.collections['EnVi Geometry'].children[self.zone_menu].vi_params['enparams']['floorarea']

                if self.energy_menu == '0':
                    if self['res']['fa'] > 13.9:
                        occ = 1 + 1.76*(1 - math.exp(-0.000349 * (self['res']['fa']-13.9)**2)) + 0.0013 * (self['res']['fa'] - 13.9)
                    else:
                        occ = 1
        
                    Vda = 25 * occ + 36
                    md = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
                    ff = (1.10, 1.06, 1.02, 0.98, 0.94, 0.90, 0.90, 0.94, 0.98, 1.02, 1.06, 1.10, 1.00)
                    dtm = (41.2, 41.4, 40.1, 37.6, 36.4, 33.9, 30.4, 33.4, 33.5, 36.3, 39.4, 39.9, 37.0)
                    self['res']['wkwh'] = 1.15 * sum([4.18/3600 * Vda * z[0] * z[1] * z[2] for z in zip(md, ff, dtm)])
            
                    for r in rl:
                        if self.zone_menu == 'All':
        #                    self['res']['fa'] = sum([c.vi_params['enparams']['floorarea'] for c in bpy.data.collections['EnVi Geometry'].children]) 
                            if r[3] == 'PV Power (W)':
                                self['res']['pvkwh'] += sum(float(p) for p in r[4].split()) * 0.001
                            elif r[3] == 'Heating (W)':
                                self['res']['hkwh'] += sum(float(p) for p in r[4].split()) * 0.001
                            elif r[3] == 'Air heating (W)':
                                self['res']['ahkwh'] += sum(float(p) for p in r[4].split()) * 0.001
                            elif r[3] == 'Cooling (W)':
                                self['res']['ckwh'] += sum(float(p) for p in r[4].split()) * 0.001
                             
                            self['res']['totkwh'] = self['res']['hkwh'] + self['res']['ahkwh'] + self['res']['ckwh'] + self.mod + self['res']['wkwh'] - self['res']['pvkwh']
                            self['res']['ECF'] = 0.42*(54 + self['res']['totkwh'] * 0.1319)/(self['res']['fa'] + 45) 
                            self['res']['EPC'] = 100 - 13.95 * self['res']['ECF'] if self['res']['ECF'] < 3.5 else 117 - 121 * math.log10(self['res']['ECF'])
                            epcletts = ('A', 'B', 'C', 'D', 'E', 'F','G')
                            epcnum = (92, 81, 69, 55, 39, 21, 1)
                            
                            for ei, en in enumerate(epcnum):
                                if self['res']['EPC'] > en:
                                    self['res']['EPCL'] = epcletts[ei]
                                    break
                        
                        else:
                            self['res']['fa'] = bpy.data.collections[self.zone_menu].vi_params['enparams']['floorarea']
                            if r[2] == self.zone_menu:
                                if r[3] == 'Heating (W)':
                                    self['res']['hkwh'] = sum(float(p) for p in r[4].split()) * 0.001
                                elif r[3] == 'Cooling (W)':        
                                    self['res']['ckwh'] = sum(float(p) for p in r[4].split()) * 0.001
                            elif r[1] == 'Power' and 'EN_' + r[2].split('_')[1] == self.zone_menu and r[3] == 'PV Power (W)':
                                    self['res']['pvkwh'] += sum(float(p) for p in r[4].split()) * 0.001
            elif self.metric == '1':
                self['res']['pvkwh'] = 'N/A'
                self['res']['hkwh'] = 'N/A'
                self['res']['ckwh'] = 'N/A'
                self['res']['fa'] = 'N/A'
                self['res']['ECF'] = 'N/A'
#    def ret_metrics(self):
#        if self.inputs['Results in'].links:
#            if self.metric_menu == '0':
#                resdict = self.inputs['Results in'].links[0].from_node['resdict']
#                print(resdict.keys())

class No_CSV(Node, ViNodes):
    '''CSV Export Node'''
    bl_idname = 'No_CSV'
    bl_label = 'VI CSV Export'
    
    animated: BoolProperty(name = '', description = 'Animated results', default = 0)

    def init(self, context):
        self.inputs.new('So_Vi_Res', 'Results in')

    def draw_buttons(self, context, layout):
        try:
            rl = self.inputs['Results in'].links[0].from_node['reslists']
            zrl = list(zip(*rl))
            if len(set(zrl[0])) > 1:
                newrow(layout, 'Animated:', self, 'animated')
            row = layout.row()
            row.operator('node.csvexport', text = 'Export CSV file')
        except:
            pass
        
    def update(self):
        pass    
      
class ViNodeCategory(NodeCategory):
    @classmethod
    def poll(cls, context):
        return context.space_data.tree_type == 'ViN'

class So_Vi_Loc(NodeSocket):
    '''Vi Location socket'''
    bl_idname = 'So_Vi_Loc'
    bl_label = 'Location socket'
    valid = ['Location']

    def draw(self, context, layout, node, text):
        layout.label(text = text)

    def draw_color(self, context, node):
        return (0.45, 1.0, 0.45, 1.0)
    
    def ret_valid(self, node):
        return ['Location']
    
class So_Li_Geo(NodeSocket):
    '''Lighting geometry socket'''
    bl_idname = 'So_Li_Geo'
    bl_label = 'Geometry'

    valid = ['LiVi Geometry', 'text']
    link_limit = 1

    def draw(self, context, layout, node, text):
        layout.label(text = text)

    def draw_color(self, context, node):
        return (0.3, 0.17, 0.07, 0.75)

class So_Li_Con(NodeSocket):
    '''Lighting context in socket'''
    bl_idname = 'So_Li_Con'
    bl_label = 'Context'

    valid = ['LiVi Context', 'text']
    link_limit = 1

    def draw(self, context, layout, node, text):
        layout.label(text = text)

    def draw_color(self, context, node):
        return (1.0, 1.0, 0.0, 0.75)
    
class So_Text(NodeSocket):
    '''VI text socket'''
    bl_idname = 'So_Text'
    bl_label = 'VI text export'

    valid = ['text']
    link_limit = 1

    def draw(self, context, layout, node, text):
        layout.label(text = text)

    def draw_color(self, context, node):
        return (0.2, 1.0, 0.0, 0.75)
    
class So_Vi_Res(NodeSocket):
    '''Vi results socket'''
    bl_idname = 'So_Vi_Res'
    bl_label = 'VI results'

    valid = ['Vi Results']
    link_limit = 1

    def draw(self, context, layout, node, text):
        layout.label(text = text)

    def draw_color(self, context, node):
        return (0.0, 1.0, 0.0, 0.75)
    
    def ret_valid(self, node):
        return ['Vi Results']
    
class So_Li_Im(NodeSocket):
    '''LiVi image socket'''
    bl_idname = 'So_Li_Im'
    bl_label = 'Image'

    valid = ['image']
    link_limit = 1

    def draw(self, context, layout, node, text):
        layout.label(text = text)

    def draw_color(self, context, node):
        return (0.5, 1.0, 0.0, 0.75)
    
class So_En_Geo(NodeSocket):
    '''EnVi geometry out socket'''
    bl_idname = 'So_En_Geo'
    bl_label = 'EnVi Geometry Socket'

    valid = ['EnVi Geometry']
    link_limit = 1

    def draw(self, context, layout, node, text):
        layout.label(text = text)

    def draw_color(self, context, node):
        return (0.0, 0.0, 1.0, 0.75)
    
class So_En_Con(NodeSocket):
    '''EnVi context socket'''
    bl_idname = 'So_En_Con'
    bl_label = 'EnVi context'

    valid = ['EnVi Context']
    link_limit = 1

    def draw(self, context, layout, node, text):
        layout.label(text = text)

    def draw_color(self, context, node):
        return (0.0, 1.0, 1.0, 0.75)
    
class So_En_Res(NodeSocket):
    '''Results socket'''
    bl_idname = 'ViEnRIn'
    bl_label = 'Results axis'
    valid = ['Vi Results']

    def draw(self, context, layout, node, text):
        typedict = {"Time": [], "Frames": [], "Climate": ['climmenu'], "Zone": ("zonemenu", "zonermenu"), 
                    "Linkage":("linkmenu", "linkrmenu"), "External":("enmenu", "enrmenu"), 
                    "Chimney":("chimmenu", "chimrmenu"), "Position":("posmenu", "posrmenu"), 
                    "Camera":("cammenu", "camrmenu"), "Power":("powmenu", "powrmenu")}
        row = layout.row()

        if self.links and self.links[0].from_node.get('frames'):
            if len(self.links[0].from_node['frames']) > 1 and node.parametricmenu == '0': 
                row.prop(self, "framemenu", text = text)
                row.prop(self, "rtypemenu")
            else:
                row.prop(self, "rtypemenu", text = text)

            for rtype in typedict[self.rtypemenu]:
                row.prop(self, rtype)
            if self.node.timemenu in ('1', '2') and self.rtypemenu !='Time' and node.parametricmenu == '0':
                row.prop(self, "statmenu")
            if self.rtypemenu != 'Time':
                row.prop(self, 'multfactor')
        else:
            row.label('No results')

    def draw_color(self, context, node):
        return (0.0, 1.0, 0.0, 0.75)
        
class So_En_ResU(NodeSocket):
    '''Vi unlinked results socket'''
    bl_idname = 'ViEnRInU'
    bl_label = 'Axis'
    valid = ['Vi Results']

    def draw_color(self, context, node):
        return (0.0, 1.0, 0.0, 0.75)

    def draw(self, context, layout, node, text):
        layout.label(text = self.bl_label)
    
        
####################### Vi Nodes Categories ##############################

vi_process = [NodeItem("No_Li_Geo", label="LiVi Geometry"), NodeItem("No_Li_Con", label="LiVi Context"), 
              NodeItem("No_En_Geo", label="EnVi Geometry"), NodeItem("No_En_Con", label="EnVi Context")]
                
vi_edit = [NodeItem("No_Text", label="Text Edit")]
vi_analysis = [NodeItem("ViSPNode", label="Sun Path"), NodeItem("ViWRNode", label="Wind Rose"), 
             NodeItem("ViSVFNode", label="Sky View"), NodeItem("ViSSNode", label="Shadow map"),
             NodeItem("No_Li_Sim", label="LiVi Simulation"), NodeItem("No_En_Sim", label="EnVi Simulation")]

vi_gen = []

vi_display = [NodeItem("No_Vi_Chart", label="Chart"), NodeItem("No_Vi_Metrics", label="Metrics")]
vi_out = [NodeItem("No_CSV", label="CSV")]
vi_image = [NodeItem("No_Li_Im", label="LiVi Image"), NodeItem("No_Li_Gl", label="LiVi Glare"), NodeItem("No_Li_Fc", label="LiVi False-colour")]
vi_input = [NodeItem("No_Loc", label="VI Location")]

vinode_categories = [ViNodeCategory("Output", "Output Nodes", items=vi_out), 
                     ViNodeCategory("Edit", "Edit Nodes", items=vi_edit), 
                     ViNodeCategory("Image", "Image Nodes", items=vi_image), 
                     ViNodeCategory("Display", "Display Nodes", items=vi_display), 
                     ViNodeCategory("Generative", "Generative Nodes", items=vi_gen), 
                     ViNodeCategory("Analysis", "Analysis Nodes", items=vi_analysis), 
                     ViNodeCategory("Process", "Process Nodes", items=vi_process), 
                     ViNodeCategory("Input", "Input Nodes", items=vi_input)]

class EnViNetwork(NodeTree):
    '''A node tree for the creation of EnVi advanced networks.'''
    bl_idname = 'EnViN'
    bl_label = 'EnVi Network'
    bl_icon = 'FORCE_WIND'
    nodetypes = {}

class EnViNodes:
    @classmethod
    def poll(cls, ntree):
        return ntree.bl_idname == 'EnViN'
    
class So_En_Net_Bound(NodeSocket):
    '''A plain zone boundary socket'''
    bl_idname = 'So_En_Net_Bound'
    bl_label = 'Plain zone boundary socket'
    bl_color = (1.0, 1.0, 0.2, 0.5)

    valid = ['Boundary']
    sn: StringProperty()
    uvalue: StringProperty()

    def draw(self, context, layout, node, text):
        layout.label(text = text)

    def draw_color(self, context, node):
        return (0.5, 0.2, 0.0, 0.75)

    def ret_valid(self, node):
        return ['Boundary']
    
class So_En_Sched(NodeSocket):
    '''Fraction schedule socket'''
    bl_idname = 'So_En_Sched'
    bl_label = 'Schedule socket'
    bl_color = (1.0, 1.0, 0.0, 0.75)

    valid = ['Schedule']
    schedule = ['Fraction']

    def draw(self, context, layout, node, text):
        layout.label(text = text)

    def draw_color(self, context, node):
        return (1.0, 1.0, 0.0, 0.75)

class So_En_Net_TSched(NodeSocket):
    '''Temperature schedule socket'''
    bl_idname = 'So_En_Net_TSched'
    bl_label = 'Schedule socket'
    bl_color = (1.0, 1.0, 0.0, 0.75)

    valid = ['Schedule']
    schedule = ['Temperature']

    def draw(self, context, layout, node, text):
        layout.label(text = text)

    def draw_color(self, context, node):
        return (1.0, 1.0, 0.0, 0.75)

class So_En_Net_SSFlow(NodeSocket):
    '''A sub-surface flow socket'''
    bl_idname = 'So_En_Net_SSFlow'
    bl_label = 'Sub-surface flow socket'

    sn: StringProperty()
    valid = ['Sub-surface']

    def draw(self, context, layout, node, text):
        layout.label(text = text)

    def draw_color(self, context, node):
        return (0.1, 1.0, 0.2, 0.75)

    def ret_valid(self, node):
        return ['Sub-surface']
    
class So_En_Net_SFlow(NodeSocket):
    '''A surface flow socket'''
    bl_idname = 'So_En_Net_SFlow'
    bl_label = 'Surface flow socket'

    sn: StringProperty()
    valid = ['Surface']

    def draw(self, context, layout, node, text):
        layout.label(text = text)

    def draw_color(self, context, node):
        return (1.0, 0.2, 0.2, 0.75)

    def ret_valid(self, node):
        return ['Surface']
    
class So_En_Net_SSSFlow(NodeSocket):
    '''A surface or sub-surface flow socket'''
    bl_idname = 'So_En_Net_SSSFlow'
    bl_label = '(Sub-)Surface flow socket'

    sn: StringProperty()
    valid = ['(Sub)Surface']

    def draw(self, context, layout, node, text):
        layout.label(text = text)

    def draw_color(self, context, node):
        return (1.0, 1.0, 0.2, 0.75)
    
    def ret_valid(self, node):
        return ['(Sub)Surface']

class So_En_Net_CRef(NodeSocket):
    '''A plain zone airflow component socket'''
    bl_idname = 'So_En_Net_CRef'
    bl_label = 'Plain zone airflow component socket'

    sn = StringProperty()

    def draw(self, context, layout, node, text):
        layout.label(text = text)

    def draw_color(self, context, node):
        return (1.0, 0.4, 0.0, 0.75)

class So_En_Net_Occ(NodeSocket):
    '''An EnVi zone occupancy socket'''
    bl_idname = 'So_En_Net_Occ'
    bl_label = 'Zone occupancy socket'

    sn: StringProperty()
    valid = ['Occupancy']

    def draw(self, context, layout, node, text):
        layout.label(text = text)

    def draw_color(self, context, node):
        return (1.0, 0.2, 0.2, 0.75)

class So_En_Net_Eq(NodeSocket):
    '''An EnVi zone equipment socket'''
    bl_idname = 'So_En_Net_Eq'
    bl_label = 'Zone equipment socket'

    sn: StringProperty()
    valid = ['Equipment']

    def draw(self, context, layout, node, text):
        layout.label(text = text)

    def draw_color(self, context, node):
        return (1.0, 0.2, 0.2, 0.75)

class So_En_Net_Inf(NodeSocket):
    '''An EnVi zone infiltration socket'''
    bl_idname = 'So_En_Net_Inf'
    bl_label = 'Zone infiltration socket'

    valid = ['Infiltration']

    def draw(self, context, layout, node, text):
        layout.label(text = text)

    def draw_color(self, context, node):
        return (1.0, 0.2, 0.2, 0.75)

class So_En_Net_Hvac(NodeSocket):
    '''An EnVi zone HVAC socket'''
    bl_idname = 'So_En_Net_Hvac'
    bl_label = 'Zone HVAC socket'

    sn: StringProperty()
    valid = ['HVAC']

    def draw(self, context, layout, node, text):
        layout.label(text = text)

    def draw_color(self, context, node):
        return (1.0, 0.2, 0.2, 0.75)

class So_En_Net_WPC(NodeSocket):
    '''An EnVi external node WPC socket'''
    bl_idname = 'So_En_Net_WPC'
    bl_label = 'External node WPC'

    sn: StringProperty()
    valid = ['WPC']

    def draw(self, context, layout, node, text):
        layout.label(text = text)

    def draw_color(self, context, node):
        return (0.2, 0.2, 0.2, 0.75)

class So_En_Net_Act(NodeSocket):
    '''An EnVi actuator socket'''
    bl_idname = 'So_En_Net_Act'
    bl_label = 'EnVi actuator socket'

    sn: StringProperty()
    valid = ['Actuator']

    def draw(self, context, layout, node, text):
        layout.label(text = text)

    def draw_color(self, context, node):
        return (0.2, 0.9, 0.9, 0.75)
    
    def ret_valid(self, node):
        return ['Actuator']
    
class So_En_Net_Sense(NodeSocket):
    '''An EnVi sensor socket'''
    bl_idname = 'So_En_Net_Sense'
    bl_label = 'EnVi sensor socket'
    
    sn: StringProperty()
    valid = ['Sensor']

    def draw(self, context, layout, node, text):
        layout.label(text = text)

    def draw_color(self, context, node):
        return (0.9, 0.9, 0.2, 0.75)
    
    def ret_valid(self, node):
        return ['Sensor']
    
#class So_En_Net_SFlow(NodeSocket):
#    '''A surface flow socket'''
#    bl_idname = 'So_En_Net_SFlow'
#    bl_label = 'Surface flow socket'
#
#    sn = StringProperty()
#    valid = ['Surface']
#
#    def draw(self, context, layout, node, text):
#        layout.label(text)
#
#    def draw_color(self, context, node):
#        return (1.0, 0.2, 0.2, 0.75)
#
#class So_En_Net_SSFlow(NodeSocket):
#    '''A surface or sub-surface flow socket'''
#    bl_idname = 'So_En_Net_SSFlow'
#    bl_label = '(Sub-)Surface flow socket'
#
#    sn = StringProperty()
#    valid = ['(Sub)Surface']
#
#    def draw(self, context, layout, node, text):
#        layout.label(text)
#
#    def draw_color(self, context, node):
#        return (1.0, 1.0, 0.2, 0.75)
    
class No_En_Net_Zone(Node, EnViNodes):
    '''Node describing a simulation zone'''
    bl_idname = 'No_En_Net_Zone'
    bl_label = 'Zone'
    bl_icon = 'SOUND'

    def zupdate(self, context):
        self.afs = 0
        col = bpy.data.collections[self.zone]

        for obj in col.objects:
            odm = [m.material for m in obj.material_slots]
            olinks = [(o.name, o.links[0].to_node.name, o.links[0].to_socket.name) for o in self.outputs if o.links and o.bl_idname in ('So_En_Net_Bound', 'So_En_Net_SFlow', 'So_En_Net_SSFlow')]
            ilinks = [(i.name, i.links[0].from_node.name, i.links[0].from_socket.name) for i in self.inputs if i.links and i.bl_idname in ('So_En_Net_Bound', 'So_En_Net_SFlow', 'So_En_Net_SSFlow')]
            bfacelist = sorted([face for face in obj.data.polygons if get_con_node(odm[face.material_index].vi_params).envi_con_con == 'Zone'], key=lambda face: -face.center[2])
    #        buvals = [retuval(odm[face.material_index]) for face in bfacelist]
            
            bsocklist = ['{}_{}_b'.format(odm[face.material_index].name, face.index) for face in bfacelist]
            sfacelist = sorted([face for face in obj.data.polygons if get_con_node(odm[face.material_index].vi_params).envi_afsurface == 1 and get_con_node(odm[face.material_index].vi_params).envi_con_type not in ('Window', 'Door')], key=lambda face: -face.center[2])
            ssocklist = ['{}_{}_s'.format(odm[face.material_index].name, face.index) for face in sfacelist]
            ssfacelist = sorted([face for face in obj.data.polygons if get_con_node(odm[face.material_index].vi_params).envi_afsurface == 1 and get_con_node(odm[face.material_index].vi_params).envi_con_type in ('Window', 'Door')], key=lambda face: -face.center[2])
            sssocklist = ['{}_{}_ss'.format(odm[face.material_index].name, face.index) for face in ssfacelist]
            [self.outputs.remove(oname) for oname in self.outputs if oname.bl_idname in ('So_En_Net_Bound', 'So_En_Net_SFlow', 'So_En_Net_SSFlow')]# and oname not in bsocklist + ssocklist + sssocklist]
            [self.inputs.remove(iname) for iname in self.inputs if iname.bl_idname in ('So_En_Net_Bound', 'So_En_Net_SFlow', 'So_En_Net_SSFlow')]# and iname not in bsocklist + ssocklist + sssocklist]
    
            for sock in bsocklist:
                self.outputs.new('So_En_Net_Bound', sock)
                self.outputs[0].sn = sock.split('_')[-2]
                self.inputs.new('So_En_Net_Bound', sock).sn = sock.split('_')[-2]
            for sock in ssocklist:
                self.afs += 1
                self.outputs.new('So_En_Net_SFlow', sock).sn = sock.split('_')[-2]
                self.inputs.new('So_En_Net_SFlow', sock).sn = sock.split('_')[-2]
            for sock in sssocklist:
                self.afs += 1
                self.outputs.new('So_En_Net_SSFlow', sock).sn = sock.split('_')[-2]
                self.inputs.new('So_En_Net_SSFlow', sock).sn = sock.split('_')[-2]                
    #        for s, sock in enumerate(bsocklist):
    #            self.outputs[sock].uvalue = '{:.4f}'.format(buvals[s])    
    #            self.inputs[sock].uvalue = '{:.4f}'.format(buvals[s]) 
            for olink in olinks:
                try:
                    self.id_data.links.new(self.outputs[olink[0]], self.id_data.nodes[olink[1]].inputs[olink[2]])
                except:
                    pass
            for ilink in ilinks:
                try:
                    self.id_data.links.new(self.id_data.nodes[ilink[1]].outputs[ilink[2]], self.inputs[ilink[0]])
                except:
                    pass
        self.vol_update(context)
        
    def vol_update(self, context):
        coll = bpy.data.collections[self.zone] 
        
        for obj in coll.objects:
            obj['volume'] = obj['auto_volume'] if obj.get('auto_volume') and self.volcalc == '0' else self.zonevolume
            self['volume'] = obj['auto_volume'] if obj.get('auto_volume') else 0
            
    def tspsupdate(self, context):
        if self.control != 'Temperature' and self.inputs['TSPSchedule'].links:
            remlink(self, self.inputs['TSPSchedule'].links)
        self.inputs['TSPSchedule'].hide = False if self.control == 'Temperature' else True
        self.update()
                
    zone: StringProperty(name = '', update = zupdate)
    controltype = [("NoVent", "None", "No ventilation control"), ("Constant", "Constant", "From vent availability schedule"), ("Temperature", "Temperature", "Temperature control")]
    control: EnumProperty(name="", description="Ventilation control type", items=controltype, default='NoVent', update=tspsupdate)
    volcalc: EnumProperty(name="", description="Volume calculation type", items=[('0', 'Auto', 'Automatic calculation (check EnVi error file)'), ('1', 'Manual', 'Manual volume')], default='0', update=vol_update)
    zonevolume: FloatProperty(name = '', min = 0, default = 100, update=vol_update)
    mvof: FloatProperty(default = 0, name = "", min = 0, max = 1)
    lowerlim: FloatProperty(default = 0, name = "", min = 0, max = 100)
    upperlim: FloatProperty(default = 50, name = "", min = 0, max = 100)
    afs: IntProperty(default = 0, name = "")
    alllinked: BoolProperty(default = 0, name = "")
    envi_oca: eprop([("0", "Default", "Use the system wide convection algorithm"), ("1", "Simple", "Use the simple convection algorithm"), ("2", "TARP", "Use the detailed convection algorithm"), ("3", "DOE-2", "Use the Trombe wall convection algorithm"), ("4", "MoWitt", "Use the adaptive convection algorithm"), ("5", "Adaptive", "Use the adaptive convection algorithm")], "", "Specify the EnVi zone outside convection algorithm", "0")
    envi_ica: eprop([("0", "Default", "Use the system wide convection algorithm"), ("1", "Simple", "Use the simple convection algorithm"), ("2", "Detailed", "Use the detailed convection algorithm"), ("3", "Trombe", "Use the Trombe wall convection algorithm"), ("4", "Adaptive", "Use the adaptive convection algorithm")], "", "Specify the EnVi zone inside convection algorithm", "0")
    envi_hab: BoolProperty(default = 0, name = "")
#    zone = StringProperty(name = '', default = "en_Chimney")
#    tcsched = EnumProperty(name="", description="Ventilation control type", items=[('On', 'On', 'Always on'), ('Off', 'Off', 'Always off'), ('Sched', 'Schedule', 'Scheduled operation')], default='On', update = supdate)
#    waw = FloatProperty(name = '', min = 0.001, default = 1)
#    ocs = FloatProperty(name = '', min = 0.001, default = 1)
#    odc = FloatProperty(name = '', min = 0.001, default = 0.6)
    
    def init(self, context):
        self['tsps'] = 1
        self['volume'] = 10
        self.inputs.new('So_En_Net_Hvac', 'HVAC')
        self.inputs.new('So_En_Net_Occ', 'Occupancy')
        self.inputs.new('So_En_Net_Eq', 'Equipment')
        self.inputs.new('So_En_Net_Inf', 'Infiltration')
        self.inputs.new('So_En_Sched', 'TSPSchedule')
        self.inputs.new('So_En_Sched', 'VASchedule')

    def update(self):
        sflowdict = {'So_En_Net_SFlow': 'Envi surface flow', 'So_En_Net_SSFlow': 'Envi sub-surface flow'}
        [bi, si, ssi, bo, so , sso] = [1, 1, 1, 1, 1, 1]
                
        try:            
            for inp in [inp for inp in self.inputs if inp.bl_idname in ('So_En_Net_Bound', 'So_En_Net_SFlow', 'So_En_Net_SSFlow')]:
                self.outputs[inp.name].hide = True if inp.links and self.outputs[inp.name].bl_idname == inp.bl_idname else False
    
            for outp in [outp for outp in self.outputs if outp.bl_idname in ('So_En_Net_Bound', 'So_En_Net_SFlow', 'So_En_Net_SSFlow')]:
                if self.inputs.get(outp.name):
                    self.inputs[outp.name].hide = True if outp.links and self.inputs[outp.name].bl_idname == outp.bl_idname else False
    
            for inp in [inp for inp in self.inputs if inp.bl_idname in ('So_En_Net_Bound', 'So_En_Net_SFlow', 'So_En_Net_SSFlow')]:
                if inp.bl_idname == 'So_En_Bound' and not inp.hide and not inp.links:
                    bi = 0
                elif inp.bl_idname in sflowdict:
                    if (not inp.hide and not inp.links) or (inp.links and inp.links[0].from_node.bl_label != sflowdict[inp.bl_idname]):
                        si = 0
                        if inp.links:
                            remlink(self, [inp.links[0]])    
            
            for outp in [outp for outp in self.outputs if outp.bl_idname in ('So_En_Net_Bound', 'So_En_Net_SFlow', 'So_En_Net_SSFlow')]:
                if outp.bl_idname == 'So_En_Bound' and not outp.hide and not outp.links:
                    bo = 0
                elif outp.bl_idname  in sflowdict:
                    if (not outp.hide and not outp.links) or (outp.links and outp.links[0].to_node.bl_label != sflowdict[outp.bl_idname]):
                        so = 0
                        if outp.links:
                            remlink(self, [outp.links[0]])

        except Exception as e:
            logentry("There was a problem an EnVi Zone {} node socket change: {}".format(self.zone, e))
        
        for sock in self.outputs:
            socklink2(sock, self.id_data)
            
        self.alllinked = 1 if all((bi, si, ssi, bo, so, sso)) else 0
        nodecolour(self, self.errorcode())
        
    def uvsockupdate(self):
        for sock in self.outputs:
            socklink2(sock, self.id_data)
            if sock.bl_idname == 'So_En_Net_Bound':
                uvsocklink2(sock, self.id_data)
    
    def errorcode(self):
        if self.afs == 1:
            return 'Too few air-flow surfaces'
        if self.control == 'Temperature' and not self.inputs['TSPSchedule'].is_linked:
            return 'TSPSchedule not linked'
        elif self.alllinked == 0:
            return 'Unlinked air-flow/boundary socket'
        else:
            return ''
                    
    def draw_buttons(self, context, layout):
        if self.errorcode():
            row = layout.row()
            row.label(text = 'Error - {}'.format(self.errorcode()))
        newrow(layout, 'Zone:', self, 'zone')
        if bpy.data.collections.get(self.zone):
#            cvp = bpy.data.collections[self.zone].vi_params
            newrow(layout, "Habitable:", self, 'envi_hab')
            newrow(layout, "Inside convection:", self, 'envi_ica')
            newrow(layout, "Outside convection:", self, 'envi_oca')
        yesno = (1, self.control == 'Temperature', self.control == 'Temperature', self.control == 'Temperature')
        vals = (("Control type:", "control"), ("Minimum OF:", "mvof"), ("Lower:", "lowerlim"), ("Upper:", "upperlim"))
        newrow(layout, 'Volume calc:', self, 'volcalc')
        
        if self.volcalc == '0':
            row = layout.row()
            row.label(text = 'Auto volume: {:.1f}'.format(self['volume']))
        else:
            newrow(layout, 'Volume:', self, 'zonevolume')
            
        [newrow(layout, val[0], self, val[1]) for v, val in enumerate(vals) if yesno[v]]
        
    def epwrite(self):
        (tempschedname, mvof, lowerlim, upperlim) = (self.zone + '_tspsched', self.mvof, self.lowerlim, self.upperlim) if self.inputs['TSPSchedule'].is_linked else ('', '', '', '')
        vaschedname = self.zone + '_vasched' if self.inputs['VASchedule'].is_linked else ''
        params = ('Zone Name',
        'Ventilation Control Mode', 'Ventilation Control Zone Temperature Setpoint Schedule Name',
        'Minimum Venting Open Factor (dimensionless)',
        'Indoor and Outdoor Temperature Diffeence Lower Limit for Maximum Venting Opening Factor (deltaC)',
        'Indoor and Outdoor Temperature Diffeence Upper Limit for Minimum Venting Opening Factor (deltaC)',
        'Indoor and Outdoor Enthalpy Difference Lower Limit For Maximum Venting Open Factor (deltaJ/kg)',
        'Indoor and Outdoor Enthalpy Difference Upper Limit for Minimum Venting Open Factor (deltaJ/kg)',
        'Venting Availability Schedule Name')

        paramvs = (self.zone, self.control, tempschedname, mvof, lowerlim, upperlim, '0.0', '300000.0', vaschedname)
        return epentry('AirflowNetwork:MultiZone:Zone', params, paramvs)
    
#    def tc_epwrite(self):
#        params =  ('Name of Thermal Chimney System', 'Name of Thermal Chimney Zone',
#                   'Availability Schedule Name', 'Width of the Absorber Wall',
#                   'Cross Sectional Area of Air Channel Outlet', 'Discharge Coefficient',
#                   'Zone Name 1', 'Distance from the Top of the Thermal Chimney to Inlet 1',
#                   'Relative Ratios of Air Flow Rates Passing through Zone 1', 'Cross Sectional Areas of Air Channel Inlet 1',
#                   'Zone Name 2', 'Distance from the Top of the Thermal Chimney to Inlet 2',
#                   'Relative Ratios of Air Flow Rates Passing through Zone 2', ' Cross Sectional Areas of Air Channel Inlet 2')

class No_En_Net_TC(Node, EnViNodes):
    '''Zone Thermal Chimney node'''
    bl_idname = 'No_En_Net_TC'
    bl_label = 'Chimney'
    bl_icon = 'SOUND'

    def zupdate(self, context):
        zonenames= []
        obj = bpy.data.objects[self.zone]
        odm = obj.data.materials
        bsocklist = ['{}_{}_b'.format(odm[face.material_index].name, face.index) for face in obj.data.polygons if get_con_node(odm[face.material_index].vi_params).envi_con_con  == 'Zone' and odm[face.material_index].name not in [outp.name for outp in self.outputs if outp.bl_idname == 'So_En_Net_Bound']]

        for oname in [outputs for outputs in self.outputs if outputs.name not in bsocklist and outputs.bl_idname == 'EnViBoundSocket']:
            self.outputs.remove(oname)
            
        for iname in [inputs for inputs in self.inputs if inputs.name not in bsocklist and inputs.bl_idname == 'EnViBoundSocket']:
            self.inputs.remove(iname)
            
        for sock in sorted(set(bsocklist)):
            if not self.outputs.get(sock):
                self.outputs.new('EnViBoundSocket', sock).sn = sock.split('_')[-2]
            if not self.inputs.get(sock):
                self.inputs.new('EnViBoundSocket', sock).sn = sock.split('_')[-2]
                
        for sock in (self.inputs[:] + self.outputs[:]):
            if sock.bl_idname == 'EnViBoundSocket' and sock.links:
                zonenames += [(link.from_node.zone, link.to_node.zone)[sock.is_output] for link in sock.links]

        nodecolour(self, all([get_con_node(mat.vi_params).envi_con_type != 'Window' for mat in bpy.data.objects[self.zone].data.materials if mat]))
        self['zonenames'] = zonenames

    def supdate(self, context):
        self.inputs.new['Schedule'].hide = False if self.sched == 'Sched' else True

    zone: StringProperty(name = '', default = "en_Chimney")
    sched: EnumProperty(name="", description="Ventilation control type", items=[('On', 'On', 'Always on'), ('Off', 'Off', 'Always off'), ('Sched', 'Schedule', 'Scheduled operation')], default='On', update = supdate)
    waw: FloatProperty(name = '', min = 0.001, default = 1)
    ocs: FloatProperty(name = '', min = 0.001, default = 1)
    odc: FloatProperty(name = '', min = 0.001, default = 0.6)

    def init(self, context):
        self.inputs.new('So_En_Sched', 'Schedule')
        self['zonenames'] = []

    def draw_buttons(self, context, layout):
        newrow(layout, 'Zone:', self, 'zone')
        newrow(layout, 'Schedule:', self, 'sched')
        newrow(layout, 'Width Absorber:', self, 'waw')
        newrow(layout, 'Outlet area:', self, 'ocs')
        newrow(layout, 'Outlet DC:', self, 'odc')

        for z, zn in enumerate(self['zonenames']):
            row=layout.row()
            row.label(zn)
            row=layout.row()
            row.prop(self, '["Distance {}"]'.format(z))
            row=layout.row()
            row.prop(self, '["Relative Ratio {}"]'.format(z))
            row=layout.row()
            row.prop(self, '["Cross Section {}"]'.format(z))

    def update(self):
        bi, bo = 1, 1
        zonenames, fheights, fareas = [], [], []
        for inp in [inp for inp in self.inputs if inp.bl_idname == 'So_En_Net_Bound']:
            self.outputs[inp.name].hide = True if inp.is_linked and self.outputs[inp.name].bl_idname == inp.bl_idname else False

        for outp in [outp for outp in self.outputs if outp.bl_idname in 'So_En_Net_Bound']:
            self.inputs[outp.name].hide = True if outp.is_linked and self.inputs[outp.name].bl_idname == outp.bl_idname else False

        if [inp for inp in self.inputs if inp.bl_idname == 'So_En_Net_Bound' and not inp.hide and not inp.links]:
            bi = 0
                
        if [outp for outp in self.outputs if outp.bl_idname == 'So_En_Net_Bound' and not outp.hide and not outp.links]:
            bo = 0
        
        nodecolour(self, not all((bi, bo)))
        
        for sock in [sock for sock in self.inputs[:] + self.outputs[:] if sock.bl_idname == 'So_En_Net_Bound']:
            if sock.links and self.zone in [o.name for o in bpy.data.objects]:
                zonenames += [link.to_node.zone for link in sock.links]
                fheights += [max([(bpy.data.objects[self.zone].matrix_world * vert.co)[2] for vert in bpy.data.objects[self.zone].data.vertices]) - (bpy.data.objects[link.to_node.zone].matrix_world * bpy.data.objects[link.to_node.zone].data.polygons[int(link.to_socket.sn)].center)[2] for link in sock.links]
                fareas += [facearea(bpy.data.objects[link.to_node.zone], bpy.data.objects[link.to_node.zone].data.polygons[int(link.to_socket.sn)]) for link in sock.links]
    
            self['zonenames'] = zonenames
            for z, zn in enumerate(self['zonenames']):
                self['Distance {}'.format(z)] = fheights[z]
                self['Relative Ratio {}'.format(z)] = 1.0
                self['Cross Section {}'.format(z)] = fareas[z]
                
        for sock in self.outputs:
            socklink(sock, self['nodeid'].split('@')[1])

    def uvsockupdate(self):
        for sock in self.outputs:
            socklink(sock, self['nodeid'].split('@')[1])
            
            if sock.bl_idname == 'EnViBoundSocket':
                uvsocklink(sock, self['nodeid'].split('@')[1])
                
    def epwrite(self):
        scheduled = 1 if self.inputs['Schedule'].links and not self.inputs['Schedule'].links[0].to_node.use_custom_color else 0
        paramvs = ('{}_TC'.format(self.zone), self.zone, ('', '{}_TCSched'.format(self.zone))[scheduled], self.waw, self.ocs, self.odc)
        params = ('Name of Thermal Chimney System', 'Name of Thermal Chimney Zone', 'Availability Schedule Name', 'Width of the Absorber Wall',
                  'Cross Sectional Area of Air Channel Outlet', 'Discharge Coefficient')

        for z, zn in enumerate(self['zonenames']):
            params += (' Zone Name {}'.format(z + 1), 'Distance from the Top of the Thermal Chimney to Inlet {}'.format(z + 1), 'Relative Ratios of Air Flow Rates Passing through Zone {}'.format(z + 1),
                       'Cross Sectional Areas of Air Channel Inlet {}'.format(z + 1))
            paramvs += (zn, self['Distance {}'.format(z)], self['Relative Ratio {}'.format(z)], self['Cross Section {}'.format(z)])

        return epentry('ZoneThermalChimney', params, paramvs)


class No_En_Net_Hvac(Node, EnViNodes):
    '''Zone HVAC node'''
    bl_idname = 'No_En_Net_Hvac'
    bl_label = 'HVAC'

    def hupdate(self, context):
        self.h = 1 if self.envi_hvachlt != '4' else 0
        self.c = 1 if self.envi_hvacclt != '4' else 0
        self['hc'] = ('', 'SingleHeating', 'SingleCooling', 'DualSetpoint')[(not self.h and not self.c, self.h and not self.c, not self.h and self.c, self.h and self.c).index(1)]
        
    envi_hvact: bprop("", "", False)
    envi_hvacht: fprop(u'\u00b0C', "Heating temperature:", 1, 99, 50)
    envi_hvacct: fprop(u'\u00b0C', "Cooling temperature:", -10, 20, 13)
    envi_hvachlt: EnumProperty(items = [('0', 'LimitFlowRate', 'LimitFlowRate'), ('1', 'LimitCapacity', 'LimitCapacity'), ('2', 'LimitFlowRateAndCapacity', 'LimitFlowRateAndCapacity'), ('3', 'NoLimit', 'NoLimit'), ('4', 'None', 'No heating')], name = '', description = "Heating limit type", default = '4', update = hupdate)
    envi_hvachaf: FloatProperty(name = u'm\u00b3/s', description = "Heating air flow rate", min = 0, max = 60, default = 1, precision = 4)
    envi_hvacshc: fprop("W", "Sensible heating capacity", 0, 10000, 1000)
    envi_hvacclt: EnumProperty(items = [('0', 'LimitFlowRate', 'LimitFlowRate'), ('1', 'LimitCapacity', 'LimitCapacity'), ('2', 'LimitFlowRateAndCapacity', 'LimitFlowRateAndCapacity'), ('3', 'NoLimit', 'NoLimit'), ('4', 'None', 'No cooling')], name = '', description = "Cooling limit type", default = '4', update = hupdate)
    envi_hvaccaf: FloatProperty(name = u'm\u00b3/s', description = "Cooling air flow rate", min = 0, max = 60, default = 1, precision = 4)
    envi_hvacscc: fprop("W", "Sensible cooling capacity", 0, 10000, 1000)
    envi_hvacoam: eprop([('0', 'None', 'None'), ('1', 'Flow/Zone', 'Flow/Zone'), ('2', 'Flow/Person', 'Flow/Person'), ('3', 'Flow/Area', 'Flow/Area'), ('4', 'Sum', 'Sum'), ('5', 'Maximum ', 'Maximum'), ('6', 'ACH/Detailed', 'ACH/Detailed')], '', "Cooling limit type", '2')
    envi_hvacfrp: fprop(u'm\u00b3/s/p', "Flow rate per person", 0, 1, 0.008)
    envi_hvacfrzfa: fprop("", "Flow rate per zone area", 0, 1, 0.008)
    envi_hvacfrz: FloatProperty(name = u'm\u00b3/s', description = "Flow rate per zone", min = 0, max = 100, default = 0.1, precision = 4)
    envi_hvacfach: fprop("", "ACH", 0, 10, 1)
    envi_hvachr: eprop([('0', 'None', 'None'), ('1', 'Sensible', 'Flow/Zone')], '', "Heat recovery type", '0')
    envi_hvachre: fprop("", "Heat recovery efficiency", 0, 1, 0.7)
    h: iprop('', '', 0, 1, 0)
    c: iprop('', '', 0, 1, 0)
    actlist = [("0", "Air supply temp", "Actuate an ideal air load system supply temperature"), ("1", "Air supply flow", "Actuate an ideal air load system flow rate"),
               ("2", "Outdoor Air supply flow", "Actuate an ideal air load system outdoor air flow rate")]
    acttype: EnumProperty(name="", description="Actuator type", items=actlist, default='0')
    compdict = {'0': 'AirFlow Network Window/Door Opening'}
    actdict =  {'0': ('Venting Opening Factor', 'of')}
    envi_heat: BoolProperty(name = "Heating", description = 'Turn on zone heating', default = 0)
    envi_htsp: FloatProperty(name = u'\u00b0C', description = "Temperature", min = 0, max = 50, default = 20)
    envi_cool: BoolProperty(name = "Cooling", description = "Turn on zone cooling", default = 0)
    envi_ctsp: FloatProperty(name = u'\u00b0'+"C", description = "Temperature", min = 0, max = 50, default = 20)

    def init(self, context):
        self['hc'] = ''
        self['ctdict'] = {'DualSetpoint': 4, 'SingleHeating': 1, 'SingleCooling': 2}
        self['limittype'] = {'0': 'LimitFlowRate', '1': 'LimitCapacity', '2': 'LimitFlowRateAndCapacity', '3': 'NoLimit', '4': ''}
        self.outputs.new('So_En_Net_Hvac', 'HVAC')
        self.inputs.new('So_En_Sched', 'Schedule')
        self.inputs.new('So_En_Net_TSched', 'HSchedule')
        self.inputs.new('So_En_Net_TSched', 'CSchedule')

    def draw_buttons(self, context, layout):
        row = layout.row()
        row.label(text = 'HVAC Template:')
        row.prop(self, 'envi_hvact')
        row = layout.row()
        row.label(text = 'Heating -----------')
        newrow(layout, 'Heating limit:', self, 'envi_hvachlt')
        if self.envi_hvachlt != '4':
            newrow(layout, 'Heating temp:', self, 'envi_hvacht')
            if self.envi_hvachlt in ('0', '2',):
                newrow(layout, 'Heating airflow:', self, 'envi_hvachaf')
            if self.envi_hvachlt in ('1', '2'):
                newrow(layout, 'Heating capacity:', self, 'envi_hvacshc')
            if not self.inputs['HSchedule'].links:
                newrow(layout, 'Thermostat level:', self, 'envi_htsp')
            newrow(layout, 'Heat recovery:', self, 'envi_hvachr')
            if self.envi_hvachr != '0':
                newrow(layout, 'HR eff.:', self, 'envi_hvachre')

        row = layout.row()
        row.label(text = 'Cooling ------------')
        newrow(layout, 'Cooling limit:', self, 'envi_hvacclt')
        if self.envi_hvacclt != '4':
            newrow(layout, 'Cooling temp:', self, 'envi_hvacct')
            if self.envi_hvacclt in ('0', '2'):
                newrow(layout, 'Cooling airflow:', self, 'envi_hvaccaf')
            if self.envi_hvacclt in ('1', '2'):
                newrow(layout, 'Cooling capacity:', self, 'envi_hvacscc')
            if not self.inputs['CSchedule'].links:
                newrow(layout, 'Thermostat level:', self, 'envi_ctsp')

        if (self.envi_hvachlt, self.envi_hvacclt) != ('4', '4'):
            row = layout.row()
            row.label(text = 'Outdoor air --------------')
            newrow(layout, 'Outdoor air:', self, 'envi_hvacoam')
            if self.envi_hvacoam in ('2', '4', '5'):
                newrow(layout, 'Flow/person:', self, 'envi_hvacfrp')
            if self.envi_hvacoam in ('1', '4', '5'):
                newrow(layout, 'Zone flow:', self, 'envi_hvacfrz')
            if self.envi_hvacoam in ('3', '4', '5'):
                newrow(layout, 'Flow/area:', self, 'envi_hvacfrzfa')
            if self.envi_hvacoam in ('4', '5', '6') and not self.envi_hvact:
                newrow(layout, 'ACH', self, 'envi_hvacfach')

    def update(self):
        for sock in self.outputs:
            socklink(sock, self.id_data.name)
    
    def eptcwrite(self, zn):
        if self['hc'] in ('SingleHeating', 'SingleCooling', 'DualSetpoint'):
            return epschedwrite(zn + '_thermocontrol', 'Control Type', ['Through: 12/31'], [['For: Alldays']], [[[['Until: 24:00,{}'.format(self['ctdict'][self['hc']])]]]])
        else:
            return ''
            
    def eptspwrite(self, zn):
        params = ['Name', 'Setpoint Temperature Schedule Name']
        if self['hc'] ==  'DualSetpoint':
            params += ['Setpoint Temperature Schedule Name 2']
            paramvs = [zn +'_tsp', zn + '_htspsched', zn + '_ctspsched']
        elif self['hc'] == 'SingleHeating':
            paramvs = [zn +'_tsp', zn + '_htspsched']
        elif self['hc'] == 'SingleCooling':
            paramvs = [zn +'_tsp', zn + '_ctspsched']

        if self['hc'] in ('SingleHeating', 'SingleCooling', 'DualSetpoint'):
            params2 = ('Name', 'Zone or Zonelist Name', 'Control Type Schedule Name', 'Control 1 Object Type', 'Control 1 Name')
            paramvs2 = (zn+'_thermostat', zn, zn +'_thermocontrol', 'ThermostatSetpoint:{}'.format(self['hc']), zn + '_tsp')
            return epentry('ThermostatSetpoint:{}'.format(self['hc']), params, paramvs) + epentry('ZoneControl:Thermostat', params2, paramvs2)
        else:
            return ''

    def ephwrite(self, zn):
        params = ('Name', 'Availability Schedule Name', 'Zone Supply Air Node Name', 'Zone Exhaust Air Node Name', 'System Inlet Air Node Name',
              "Maximum Heating Supply Air Temperature (degC)", "Minimum Cooling Supply Air Temperature (degC)",
              'Maximum Heating Supply Air Humidity Ratio (kgWater/kgDryAir)', 'Minimum Cooling Supply Air Humidity Ratio (kgWater/kgDryAir)',
              'Heating Limit', 'Maximum Heating Air Flow Rate (m3/s)', 'Maximum Sensible Heating Capacity (W)',
              'Cooling limit', 'Maximum Cooling Air Flow Rate (m3/s)', 'Maximum Total Cooling Capacity (W)', 'Heating Availability Schedule Name',
              'Cooling Availability Schedule Name', 'Dehumidification Control Type', 'Cooling Sensible Heat Ratio (dimensionless)', 'Humidification Control Type',
              'Design Specification Outdoor Air Object Name', 'Outdoor Air Inlet Node Name', 'Demand Controlled Ventilation Type', 'Outdoor Air Economizer Type',
              'Heat Recovery Type', 'Sensible Heat Recovery Effectiveness (dimensionless)', 'Latent Heat Recovery Effectiveness (dimensionless)')
        paramvs = ('{}_Air'.format(zn), zn + '_hvacsched', '{}_supairnode'.format(zn), '', '', self.envi_hvacht, self.envi_hvacct, 0.015, 0.009, self['limittype'][self.envi_hvachlt],
                   '{:.4f}'.format(self.envi_hvachaf) if self.envi_hvachlt in ('0', '2') else '', self.envi_hvacshc if self.envi_hvachlt in ('1', '2') else '', self['limittype'][self.envi_hvacclt],
                   '{:.4f}'.format(self.envi_hvaccaf) if self.envi_hvacclt in ('0', '2') else '', self.envi_hvacscc if self.envi_hvacclt in ('1', '2') else '',
                   '', '', 'ConstantSupplyHumidityRatio', '', 'ConstantSupplyHumidityRatio', (zn + ' Outdoor Air', '')[self.envi_hvacoam == '0'], '', '', '', ('None', 'Sensible')[int(self.envi_hvachr)], self.envi_hvachre, '')
        entry = epentry('ZoneHVAC:IdealLoadsAirSystem', params, paramvs)

        if self.envi_hvacoam != '0':
            oam = {'0':'None', '1':'Flow/Zone', '2':'Flow/Person', '3':'Flow/Area', '4':'Sum', '5':'Maximum', '6':'AirChanges/Hour'}
            params2 = ('Name', 'Outdoor Air  Method', 'Outdoor Air Flow per Person (m3/s)', 'Outdoor Air Flow per Zone Floor Area (m3/s-m2)', 'Outdoor Air  Flow per Zone',
            'Outdoor Air Flow Air Changes per Hour', 'Outdoor Air Flow Rate Fraction Schedule Name')
            paramvs2 =(zn + ' Outdoor Air', oam[self.envi_hvacoam], '{:.4f}'.format(self.envi_hvacfrp) if self.envi_hvacoam in ('2', '4', '5') else '',
                        '{:.4f}'.format(self.envi_hvacfrzfa) if self.envi_hvacoam in ('3', '4', '5') else '', '{:.4f}'.format(self.envi_hvacfrz) if self.envi_hvacoam in ('1', '4', '5') else '',
                        '{:.4f}'.format(self.envi_hvacfach) if self.envi_hvacoam in ('4', '5', '6') else '', '')
            entry += epentry('DesignSpecification:OutdoorAir', params2, paramvs2)
        return entry

    def hvactwrite(self, zn):
        self.hupdate()
        oam = {'0':'None', '1':'Flow/Zone', '2':'Flow/Person', '3':'Flow/Area', '4':'Sum', '5':'Maximum', '6':'DetailedSpecification'}
        params = ('Zone Name' , 'Thermostat Name', 'System Availability Schedule Name', 'Maximum Heating Supply Air Temperature', 'Minimum Cooling Supply Air Temperature',
                'Maximum Heating Supply Air Humidity Ratio (kgWater/kgDryAir)', 'Minimum Cooling Supply Air Humidity Ratio (kgWater/kgDryAir)', 'Heating Limit', 'Maximum Heating Air Flow Rate (m3/s)',
                'Maximum Sensible Heating Capacity (W)', 'Cooling Limit', 'Maximum Cooling Air Flow Rate (m3/s)', 'Maximum Total Cooling Capacity (W)', 'Heating Availability Schedule Name',
                'Cooling Availability Schedule Name', 'Dehumidification Control Type', 'Cooling Sensible Heat Ratio', 'Dehumidification Setpoint (percent)', 'Humidification Control Type',
                'Humidification Setpoint (percent)', 'Outdoor Air Method', 'Outdoor Air Flow Rate per Person (m3/s)', 'Outdoor Air Flow Rate per Zone Floor (m3/s-m2)', 'Outdoor Air Flow Rate per Zone (m3/s)',
                'Design Specification Outdoor Air Object', 'Demand Controlled Ventilation Type', 'Outdoor Air Economizer Type', 'Heat Recovery Type', 'Sensible Heat Recovery Effectiveness',
                'Latent Heat Recovery Effectiveness')
        paramvs = (zn, '', zn + '_hvacsched', self.envi_hvacht, self.envi_hvacct, 0.015, 0.009, self['limittype'][self.envi_hvachlt], self.envi_hvachaf if self.envi_hvachlt in ('0', '2') else '',
                   self.envi_hvacshc if self.envi_hvachlt in ('1', '2') else '', self['limittype'][self.envi_hvacclt], self.envi_hvaccaf if self.envi_hvacclt in ('0', '2') else '',
                    self.envi_hvacscc if self.envi_hvacclt in ('1', '2') else '', '', '', 'None', '', '', 'None', '', oam[self.envi_hvacoam], '{:.4f}'.format(self.envi_hvacfrp) if self.envi_hvacoam in ('2', '4', '5') else '',
                    '{:.4f}'.format(self.envi_hvacfrzfa) if self.envi_hvacoam in ('3', '4', '5') else '', '{:.4f}'.format(self.envi_hvacfrz) if self.envi_hvacoam in ('1', '4', '5') else '', '', 'None', 'NoEconomizer', ('None', 'Sensible')[int(self.envi_hvachr)], self.envi_hvachre, 0.65)
        bpy.context.scene['enparams']['hvactemplate'] = 1
        return epentry('HVACTemplate:Zone:IdealLoadsAirSystem', params, paramvs)

    def epewrite(self, zn):
        params = ('Zone Name', 'Zone Conditioning Equipment List Name', 'Zone Air Inlet Node or NodeList Name', 'Zone Air Exhaust Node or NodeList Name',
                  'Zone Air Node Name', 'Zone Return Air Node Name')
        paramvs = (zn, zn + '_Equipment', zn + '_supairnode', '', zn + '_airnode', zn + '_retairnode')
        params2 = ('Name', 'Load Distribution Scheme', 'Zone Equipment 1 Object Type', 'Zone Equipment 1 Name', 'Zone Equipment 1 Cooling Sequence', 'Zone Equipment 1 Heating or No-Load Sequence')
        paramvs2 = (zn + '_Equipment', 'SequentialLoad', 'ZoneHVAC:IdealLoadsAirSystem', zn + '_Air', 1, 1)
        return epentry('ZoneHVAC:EquipmentConnections', params, paramvs) + epentry('ZoneHVAC:EquipmentList', params2, paramvs2)

    def schedwrite(self, zn):
        pass
    
class No_En_Net_Occ(Node, EnViNodes):
    '''Zone occupancy node'''
    bl_idname = 'No_En_Net_Occ'
    bl_label = 'Occupancy'
    bl_icon = 'SOUND'

    envi_occwatts: IntProperty(name = "W/p", description = "Watts per person", min = 1, max = 800, default = 90)
    envi_weff: FloatProperty(name = "", description = "Work efficiency", min = 0, max = 1, default = 0.0)
    envi_airv: FloatProperty(name = "", description = "Average air velocity", min = 0, max = 1, default = 0.1)
    envi_cloth: FloatProperty(name = "", description = "Clothing level", min = 0, max = 10, default = 0.5)
    envi_occtype: EnumProperty(items = [("0", "None", "No occupancy"),("1", "Occupants", "Actual number of people"), ("2", "Person/m"+ u'\u00b2', "Number of people per squared metre floor area"),
                                              ("3", "m"+ u'\u00b2'+"/Person", "Floor area per person")], name = "", description = "The type of zone occupancy specification", default = "0")
    envi_occsmax: FloatProperty(name = "", description = "Maximum level of occupancy that will occur in this schedule", min = 1, max = 500, default = 1)
    envi_comfort: BoolProperty(name = "", description = "Enable comfort calculations for this space", default = False)
    envi_co2: BoolProperty(name = "", description = "Enable CO2 concentration calculations", default = False)

    def init(self, context):
        self.outputs.new('So_En_Net_Occ', 'Occupancy')
        self.inputs.new('So_En_Sched', 'OSchedule')
        self.inputs.new('So_En_Sched', 'ASchedule')
        self.inputs.new('So_En_Sched', 'WSchedule')
        self.inputs.new('So_En_Sched', 'VSchedule')
        self.inputs.new('So_En_Sched', 'CSchedule')

    def draw_buttons(self, context, layout):
        newrow(layout, 'Type:', self, "envi_occtype")
        if self.envi_occtype != '0':
            newrow(layout, 'Max level:', self, "envi_occsmax")
            if not self.inputs['ASchedule'].links:
                newrow(layout, 'Activity level:', self, 'envi_occwatts')
            newrow(layout, 'Comfort calc:', self, 'envi_comfort')
            if self.envi_comfort:
                if not self.inputs['WSchedule'].links:
                    newrow(layout, 'Work efficiency:', self, 'envi_weff')
                if not self.inputs['VSchedule'].links:
                    newrow(layout, 'Air velocity:', self, 'envi_airv')
                if not self.inputs['CSchedule'].links:
                    newrow(layout, 'Clothing:', self, 'envi_cloth')
                newrow(layout, 'CO2:', self, 'envi_co2')

    def update(self):
        if self.inputs.get('CSchedule'):
            for sock in  self.outputs:
                socklink(sock, self.id_data.name)

    def epwrite(self, zn):
        pdict = {'0': '', '1':'People', '2': 'People/Area', '3': 'Area/Person'}
        plist = ['', '', '']
        plist[int(self.envi_occtype) - 1] = self.envi_occsmax
        params =  ['Name', 'Zone or ZoneList Name', 'Number of People Schedule Name', 'Number of People Calculation Method', 'Number of People', 'People per Zone Floor Area (person/m2)',
        'Zone Floor Area per Person (m2/person)', 'Fraction Radiant', 'Sensible Heat Fraction', 'Activity Level Schedule Name']
        paramvs = [zn + "_occupancy", zn, zn + '_occsched', pdict[self.envi_occtype]] + plist + [0.3, '', zn + '_actsched']
        if self.envi_comfort:
            params += ['Carbon Dioxide Generation Rate (m3/s-W)', 'Enable ASHRAE 55 Comfort Warnings',
                       'Mean Radiant Temperature Calculation Type', 'Surface Name/Angle Factor List Name', 'Work Efficiency Schedule Name', 'Clothing Insulation Calculation Method', 'Clothing Insulation Calculation Method Schedule Name',
                       'Clothing Insulation Schedule Name', 'Air Velocity Schedule Name', 'Thermal Comfort Model 1 Type']
            paramvs += [3.82E-8, 'No', 'zoneaveraged', '', zn + '_wesched', 'ClothingInsulationSchedule', '', zn + '_closched', zn + '_avsched', 'FANGER']
        return epentry('People', params, paramvs)

class No_En_Net_Eq(Node, EnViNodes):
    '''EnVi equipment node'''
    bl_idname = 'No_En_Net_Eq'
    bl_label = 'Equipment'
    bl_icon = 'SOUND'

    envi_equiptype: EnumProperty(items = [("0", "None", "No equipment"),("1", "EquipmentLevel", "Overall equpiment gains"), ("2", "Watts/Area", "Equipment gains per square metre floor area"),
                                              ("3", "Watts/Person", "Equipment gains per occupant")], name = "", description = "The type of zone equipment gain specification", default = "0")
    envi_equipmax: FloatProperty(name = "", description = "Maximum level of equipment gain", min = 1, max = 50000, default = 1)

    def init(self, context):
        self.outputs.new('So_En_Net_Eq', 'Equipment')
        self.inputs.new('So_En_Sched', 'Schedule')

    def draw_buttons(self, context, layout):
        newrow(layout, 'Type:', self, "envi_equiptype")
        if self.envi_equiptype != '0':
            newrow(layout, 'Max level:', self, "envi_equipmax")

    def update(self):
        for sock in self.outputs:
            socklink(sock, self.id_data.name)

    def oewrite(self, zn):
        edict = {'0': '', '1':'EquipmentLevel', '2': 'Watts/Area', '3': 'Watts/Person'}
        elist = ['', '', '']
        elist[int(self.envi_equiptype) - 1] = self.envi_equipmax
        params = ('Name', 'Fuel type', 'Zone Name', 'SCHEDULE Name', 'Design Level calculation method', 'Design Level (W)', 'Power per Zone Floor Area (Watts/m2)', 'Power per Person (Watts/person)', \
        'Fraction Latent', 'Fraction Radiant', 'Fraction Lost')
        paramvs = [zn + "_equip", 'Electricity', zn, zn + "_eqsched", edict[self.envi_equiptype]] + elist + ['0', '0', '0']
        return epentry('OtherEquipment', params, paramvs)

class No_En_Net_Inf(Node, EnViNodes):
    '''EnVi infiltration node'''
    bl_idname = 'No_En_Net_Inf'
    bl_label = 'Infiltration'
    bl_icon = 'SOUND'

    envi_inftype: EnumProperty(items = [("0", "None", "No infiltration"), ("1", 'Flow/Zone', "Absolute flow rate in m{}/s".format(u'\u00b3')), ("2", "Flow/Area", 'Flow in m{}/s per m{} floor area'.format(u'\u00b3', u'\u00b2')),
                                 ("3", "Flow/ExteriorArea", 'Flow in m{}/s per m{} external surface area'.format(u'\u00b3', u'\u00b2')), ("4", "Flow/ExteriorWallArea", 'Flow in m{}/s per m{} external wall surface area'.format(u'\u00b3', u'\u00b2')),
                                 ("5", "ACH", "ACH flow rate")], name = "", description = "The type of zone infiltration specification", default = "0")
    unit = {'0':'', '1': '(m{}/s)'.format(u'\u00b3'), '2': '(m{}/s.m{})'.format(u'\u00b3', u'\u00b2'), '3': '(m{}/s per m{})'.format(u'\u00b3', u'\u00b2'), '4': '(m{}/s per m{})'.format(u'\u00b3', u'\u00b2'), "5": "(ACH)"}
    envi_inflevel: FloatProperty(name = "", description = "Level of Infiltration", min = 0, max = 500, default = 0.001)

    def init(self, context):
        self.outputs.new('So_En_Net_Inf', 'Infiltration')
        self.inputs.new('So_En_Sched', 'Schedule')

    def draw_buttons(self, context, layout):
        newrow(layout, 'Type:', self, "envi_inftype")
        if self.envi_inftype != '0':
            newrow(layout, 'Level {}:'.format(self.unit[self.envi_inftype]), self, "envi_inflevel")

    def update(self):
        for sock in self.outputs:
            socklink(sock, self.id_data.name)

    def epwrite(self, zn):
        infildict = {'0': '', '1': 'Flow/Zone', '2': 'Flow/Area', '3': 'Flow/ExteriorArea', '4': 'Flow/ExteriorWallArea',
                          '5': 'AirChanges/Hour', '6': 'Flow/Zone'}
        inflist = ['', '', '', '']
        infdict = {'1': '0', '2': '1', '3':'2', '4':'2', '5': '3', '6': '0'}
        inflist[int(infdict[self.envi_inftype])] = '{:.4f}'.format(self.envi_inflevel)
        params = ('Name', 'Zone or ZoneList Name', 'Schedule Name', 'Design Flow Rate Calculation Method', 'Design Flow Rate {m3/s}', 'Flow per Zone Floor Area {m3/s-m2}',
               'Flow per Exterior Surface Area {m3/s-m2}', 'Air Changes per Hour {1/hr}', 'Constant Term Coefficient', 'Temperature Term Coefficient',
                'Velocity Term Coefficient', 'Velocity Squared Term Coefficient')
        paramvs = [zn + '_infiltration', zn, zn + '_infsched', infildict[self.envi_inftype]] + inflist + [1, 0, 0, 0]
        return epentry('ZoneInfiltration:DesignFlowRate', params, paramvs)

#class No_En_Net_TC(Node, EnViNodes):
#    '''EnVi Thermal Chimney node'''
#    bl_idname = 'No_En_Net_TC'
#    bl_label = 'Chimney'
#    bl_icon = 'SOUND'
#
#    def zupdate(self, context):
#        zonenames= []
#        obj = bpy.data.objects[self.zone]
#        odm = obj.data.materials
#        bsocklist = ['{}_{}_b'.format(odm[face.material_index].name, face.index) for face in obj.data.polygons if get_con_node(odm[face.material_index]).envi_boundary == 1 and odm[face.material_index].name not in [outp.name for outp in self.outputs if outp.bl_idname == 'EnViBoundSocket']]
#
#        for oname in [outputs for outputs in self.outputs if outputs.name not in bsocklist and outputs.bl_idname == 'EnViBoundSocket']:
#            self.outputs.remove(oname)
#            
#        for iname in [inputs for inputs in self.inputs if inputs.name not in bsocklist and inputs.bl_idname == 'EnViBoundSocket']:
#            self.inputs.remove(iname)
#            
#        for sock in sorted(set(bsocklist)):
#            if not self.outputs.get(sock):
#                self.outputs.new('So_En_Net_Bound', sock).sn = sock.split('_')[-2]
#            if not self.inputs.get(sock):
#                self.inputs.new('So_En_Net_Bound', sock).sn = sock.split('_')[-2]
#                
#        for sock in (self.inputs[:] + self.outputs[:]):
#            if sock.bl_idname == 'So_En_Net_Bound' and sock.links:
#                zonenames += [(link.from_node.zone, link.to_node.zone)[sock.is_output] for link in sock.links]
#
#        nodecolour(self, all([get_con_node(mat).envi_con_type != 'Window' for mat in bpy.data.objects[self.zone].data.materials if mat]))
#        self['zonenames'] = zonenames
#
#    def supdate(self, context):
#        self.inputs.new['Schedule'].hide = False if self.sched == 'Sched' else True
#
#    zone: StringProperty(name = '', default = "en_Chimney")
#    sched: EnumProperty(name="", description="Ventilation control type", items=[('On', 'On', 'Always on'), ('Off', 'Off', 'Always off'), ('Sched', 'Schedule', 'Scheduled operation')], default='On', update = supdate)
#    waw: FloatProperty(name = '', min = 0.001, default = 1)
#    ocs: FloatProperty(name = '', min = 0.001, default = 1)
#    odc: FloatProperty(name = '', min = 0.001, default = 0.6)
#
#    def init(self, context):
#        self.inputs.new('So_En_Sched', 'Schedule')
#        self['zonenames'] = []
#
#    def draw_buttons(self, context, layout):
#        newrow(layout, 'Zone:', self, 'zone')
#        newrow(layout, 'Schedule:', self, 'sched')
#        newrow(layout, 'Width Absorber:', self, 'waw')
#        newrow(layout, 'Outlet area:', self, 'ocs')
#        newrow(layout, 'Outlet DC:', self, 'odc')
#
#        for z, zn in enumerate(self['zonenames']):
#            row=layout.row()
#            row.label(zn)
#            row=layout.row()
#            row.prop(self, '["Distance {}"]'.format(z))
#            row=layout.row()
#            row.prop(self, '["Relative Ratio {}"]'.format(z))
#            row=layout.row()
#            row.prop(self, '["Cross Section {}"]'.format(z))
#
#    def update(self):
#        bi, bo = 1, 1
#        zonenames, fheights, fareas = [], [], []
#        for inp in [inp for inp in self.inputs if inp.bl_idname == 'So_En_Net_Bound']:
#            self.outputs[inp.name].hide = True if inp.is_linked and self.outputs[inp.name].bl_idname == inp.bl_idname else False
#
#        for outp in [outp for outp in self.outputs if outp.bl_idname in 'So_En_Net_Bound']:
#            self.inputs[outp.name].hide = True if outp.is_linked and self.inputs[outp.name].bl_idname == outp.bl_idname else False
#
#        if [inp for inp in self.inputs if inp.bl_idname == 'So_En_Net_Bound' and not inp.hide and not inp.links]:
#            bi = 0
#                
#        if [outp for outp in self.outputs if outp.bl_idname == 'So_En_Net_Bound' and not outp.hide and not outp.links]:
#            bo = 0
#        
#        nodecolour(self, not all((bi, bo)))
#        
#        for sock in [sock for sock in self.inputs[:] + self.outputs[:] if sock.bl_idname == 'So_En_Net_Bound']:
#            if sock.links and self.zone in [o.name for o in bpy.data.objects]:
#                zonenames += [link.to_node.zone for link in sock.links]
#                fheights += [max([(bpy.data.objects[self.zone].matrix_world * vert.co)[2] for vert in bpy.data.objects[self.zone].data.vertices]) - (bpy.data.objects[link.to_node.zone].matrix_world * bpy.data.objects[link.to_node.zone].data.polygons[int(link.to_socket.sn)].center)[2] for link in sock.links]
#                fareas += [facearea(bpy.data.objects[link.to_node.zone], bpy.data.objects[link.to_node.zone].data.polygons[int(link.to_socket.sn)]) for link in sock.links]
#    
#            self['zonenames'] = zonenames
#            for z, zn in enumerate(self['zonenames']):
#                self['Distance {}'.format(z)] = fheights[z]
#                self['Relative Ratio {}'.format(z)] = 1.0
#                self['Cross Section {}'.format(z)] = fareas[z]
#                
#        for sock in self.outputs:
#            socklink(sock, self.id_data.name)
#
#    def uvsockupdate(self):
#        for sock in self.outputs:
#            socklink(sock, self.id_data.name)
#            
#            if sock.bl_idname == 'EnViBoundSocket':
#                uvsocklink(sock, self.id_data.name)
#                
#    def epwrite(self):
#        scheduled = 1 if self.inputs['Schedule'].links and not self.inputs['Schedule'].links[0].to_node.use_custom_color else 0
#        paramvs = ('{}_TC'.format(self.zone), self.zone, ('', '{}_TCSched'.format(self.zone))[scheduled], self.waw, self.ocs, self.odc)
#        params = ('Name of Thermal Chimney System', 'Name of Thermal Chimney Zone', 'Availability Schedule Name', 'Width of the Absorber Wall',
#                  'Cross Sectional Area of Air Channel Outlet', 'Discharge Coefficient')
#
#        for z, zn in enumerate(self['zonenames']):
#            params += (' Zone Name {}'.format(z + 1), 'Distance from the Top of the Thermal Chimney to Inlet {}'.format(z + 1), 'Relative Ratios of Air Flow Rates Passing through Zone {}'.format(z + 1),
#                       'Cross Sectional Areas of Air Channel Inlet {}'.format(z + 1))
#            paramvs += (zn, self['Distance {}'.format(z)], self['Relative Ratio {}'.format(z)], self['Cross Section {}'.format(z)])
#
#        return epentry('ZoneThermalChimney', params, paramvs)

class No_En_Net_SSFlow(Node, EnViNodes):
    '''Sub-surface airflow node'''
    bl_idname = 'No_En_Net_SSFlow'
    bl_label = 'Envi sub-surface flow'
    bl_icon = 'SOUND'

    def supdate(self, context):
        if self.linkmenu in ('Crack', 'EF', 'ELA') or self.controls != 'Temperature':
            if self.inputs['TSPSchedule'].is_linked:
                remlink(self, self.inputs['TSPSchedule'].links)
        if self.linkmenu in ('Crack', 'EF', 'ELA') or self.controls in ('ZoneLevel', 'NoVent'):
            if self.inputs['VASchedule'].is_linked:
                remlink(self, self.inputs['VASchedule'].links)

        self.inputs['TSPSchedule'].hide = False if self.linkmenu in ('SO', 'DO', 'HO') and self.controls == 'Temperature' else True
        self.inputs['VASchedule'].hide = False if self.linkmenu in ('SO', 'DO', 'HO') else True
        self.legal()

    linktype = [("SO", "Simple Opening", "Simple opening element"),("DO", "Detailed Opening", "Detailed opening element"),
        ("HO", "Horizontal Opening", "Horizontal opening element"),("Crack", "Crack", "Crack aperture used for leakage calculation"),
        ("ELA", "ELA", "Effective leakage area")]

    linkmenu: EnumProperty(name="Type", description="Linkage type", items=linktype, default='SO', update = supdate)

    wdof1: FloatProperty(default = 0.1, min = 0.001, max = 1, name = "", description = 'Opening Factor 1 (dimensionless)')
    controltype = [("ZoneLevel", "ZoneLevel", "Zone level ventilation control"), ("NoVent", "None", "No ventilation control"),
                   ("Constant", "Constant", "From vent availability schedule"), ("Temperature", "Temperature", "Temperature control")]
    controls: EnumProperty(name="", description="Ventilation control type", items=controltype, default='ZoneLevel', update = supdate)
    mvof: FloatProperty(default = 0, min = 0, max = 1, name = "", description = 'Minimium venting open factor')
    lvof: FloatProperty(default = 0, min = 0, max = 100, name = "", description = 'Indoor and Outdoor Temperature Difference Lower Limit For Maximum Venting Open Factor (deltaC)')
    uvof: FloatProperty(default = 1, min = 1, max = 100, name = "", description = 'Indoor and Outdoor Temperature Difference Upper Limit For Minimum Venting Open Factor (deltaC)')
    amfcc: FloatProperty(default = 0.001, min = 0.00001, max = 1, precision = 5, name = "", description = 'Air Mass Flow Coefficient When Opening is Closed (kg/s-m)')
    amfec: FloatProperty(default = 0.65, min = 0.5, max = 1, name = '', description =  'Air Mass Flow Exponent When Opening is Closed (dimensionless)')
    lvo: EnumProperty(items = [('NonPivoted', 'NonPivoted', 'Non pivoting opening'), ('HorizontallyPivoted', 'HPivoted', 'Horizontally pivoting opening')], name = '', default = 'NonPivoted', description = 'Type of Rectanguler Large Vertical Opening (LVO)')
    ecl: FloatProperty(default = 0.0, min = 0, name = '', description = 'Extra Crack Length or Height of Pivoting Axis (m)')
    noof: IntProperty(default = 2, min = 2, max = 4, name = '', description = 'Number of Sets of Opening Factor Data')
    spa: IntProperty(default = 90, min = 0, max = 90, name = '', description = 'Sloping Plane Angle')
    dcof: FloatProperty(default = 0.7, min = 0.01, max = 1, name = '', description = 'Discharge Coefficient')
    ddtw: FloatProperty(default = 0.001, min = 0, max = 10, name = '', description = 'Minimum Density Difference for Two-way Flow')
    amfc: FloatProperty(min = 0.001, max = 1, default = 0.01, precision = 5, name = "")
    amfe: FloatProperty(min = 0.5, max = 1, default = 0.65, precision = 3, name = "")
    dlen: FloatProperty(default = 2, name = "")
    dhyd: FloatProperty(default = 0.1, name = "")
    dcs: FloatProperty(default = 0.1, name = "")
    dsr: FloatProperty(default = 0.0009, name = "")
    dlc: FloatProperty(default = 1.0, name = "")
    dhtc: FloatProperty(default = 0.772, name = "")
    dmtc: FloatProperty(default = 0.0001, name = "")
    fe: FloatProperty(default = 0.6, min = 0, max = 1, name = "")
    rpd: FloatProperty(default = 4, min = 0.1, max = 50, name = "")
    of1:FloatProperty(default = 0.0, min = 0.0, max = 0, name = '', description = 'Opening Factor {} (dimensionless)')
    (of2, of3, of4) =  [FloatProperty(default = 1.0, min = 0.01, max = 1, name = '', description = 'Opening Factor {} (dimensionless)'.format(i)) for i in range(3)]
    (dcof1, dcof2, dcof3, dcof4) = [FloatProperty(default = 0.65, min = 0.01, max = 1, name = '', description = 'Discharge Coefficient for Opening Factor {} (dimensionless)'.format(i)) for i in range(4)]
    (wfof1, wfof2, wfof3, wfof4) = [FloatProperty(default = 0.0, min = 0, max = 1, name = '', description = 'Width Factor for Opening Factor {} (dimensionless)'.format(i)) for i in range(4)]
    (hfof1, hfof2, hfof3, hfof4) = [FloatProperty(default = 0.0, min = 0, max = 1, name = '', description = 'Height Factor for Opening Factor {} (dimensionless)'.format(i)) for i in range(4)]
    (sfof1, sfof2, sfof3, sfof4) = [FloatProperty(default = 0.0, min = 0, max = 1, name = '', description = 'Start Height Factor for Opening Factor {} (dimensionless)'.format(i)) for i in range(4)]
    dcof: FloatProperty(default = 0.65, min = 0.01, max = 1, name = '', description = 'Discharge Coefficient')
    extnode:  BoolProperty(default = 0)
    actlist = [("0", "Opening factor", "Actuate the opening factor")]
    acttype: EnumProperty(name="", description="Actuator type", items=actlist, default='0')
    compdict = {'0': 'AirFlow Network Window/Door Opening'}
    actdict =  {'0': ('Venting Opening Factor', 'of')}
    adict = {'Window': 'win', 'Door': 'door'}

    def init(self, context):
        self['init'] = 1
        self['ela'] = 1.0
        self.inputs.new('So_En_Sched', 'VASchedule')
        self.inputs.new('So_En_Sched', 'TSPSchedule')
        self.inputs['TSPSchedule'].hide = True
        self.inputs.new('So_En_Net_SSFlow', 'Node 1', identifier = 'Node1_s')
        self.inputs.new('So_En_Net_SSFlow', 'Node 2', identifier = 'Node2_s')
        self.outputs.new('So_En_Net_SSFlow', 'Node 1', identifier = 'Node1_s')
        self.outputs.new('So_En_Net_SSFlow', 'Node 2', identifier = 'Node2_s')
        self.color = (1.0, 0.3, 0.3)
        self['layoutdict'] = {'SO':(('Closed FC', 'amfcc'), ('Closed FE', 'amfec'), ('Density diff', 'ddtw'), ('DC', 'dcof')), 'DO':(('Closed FC', 'amfcc'), ('Closed FE', 'amfec'),
                           ('Opening type', 'lvo'), ('Crack length', 'ecl'), ('OF Number', 'noof'), ('OF1', 'of1'), ('DC1', 'dcof1'), ('Width OF1', 'wfof1'), ('Height OF1', 'hfof1'),
                            ('Start height OF1', 'sfof1'), ('OF2', 'of2'), ('DC2', 'dcof2'), ('Width OF2', 'wfof2'), ('Height OF2', 'hfof2'), ('Start height OF2', 'sfof2')),
                            'OF3': (('OF3', 'of3'), ('DC3', 'dcof3'), ('Width OF3', 'wfof3'), ('Height OF3', 'hfof3'), ('Start height OF3', 'sfof3')),
                            'OF4': (('OF4', 'of4'), ('DC4', 'dcof4'), ('Width OF4', 'wfof4'), ('Height OF4', 'hfof4'), ('Start height OF4', 'sfof4')),
                            'HO': (('Closed FC', 'amfcc'), ('Closed FE', 'amfec'), ('Slope', 'spa'), ('DC', 'dcof')), 'Crack': (('Coefficient', 'amfc'), ('Exponent', 'amfe'), ('Factor', 'of1')),
                            'ELA': (('ELA', '["ela"]'), ('DC', 'dcof'), ('PA diff', 'rpd'), ('FE', 'fe'))}

    def update(self):
        if self.get('layoutdict'):
            for sock in self.outputs:
                socklink(sock, self.id_data.name)
            if self.linkmenu == 'ELA':
                retelaarea(self)
            self.extnode = 0
            for sock in self.inputs[:] + self.outputs[:]:
                for l in sock.links:
                    if (l.from_node, l.to_node)[sock.is_output].bl_idname == 'No_En_Net_Ext':
                        self.extnode = 1
            if self.outputs.get('Node 2'):
                sockhide(self, ('Node 1', 'Node 2'))
            self.legal()

    def draw_buttons(self, context, layout):
        layout.prop(self, 'linkmenu')
        if self.linkmenu in ('SO', 'DO', 'HO'):
            newrow(layout, 'Win/Door OF:', self, 'wdof1')
            newrow(layout, "Control type:", self, 'controls')
            if self.linkmenu in ('SO', 'DO') and self.controls == 'Temperature':
                newrow(layout, "Limit OF:", self, 'mvof')
                newrow(layout, "Lower OF:", self, 'lvof')
                newrow(layout, "Upper OF:", self, 'uvof')

        row = layout.row()
        row.label(text = 'Component options:')

        for vals in self['layoutdict'][self.linkmenu]:
            newrow(layout, vals[0], self, vals[1])
        if self.noof > 2:
            for of3vals in self['layoutdict']['OF3']:
                newrow(layout, of3vals[0], self, of3vals[1])
            if self.noof > 3:
                for of4vals in self['layoutdict']['OF4']:
                    newrow(layout, of4vals[0], self, of4vals[1])

    def epwrite(self, exp_op, enng):
        surfentry, en, snames = '', '', []
        tspsname = '{}_tspsched'.format(self.name) if self.inputs['TSPSchedule'].is_linked and self.linkmenu in ('SO', 'DO', 'HO') and self.controls == 'Temperature' else ''
        vasname = '{}_vasched'.format(self.name) if self.inputs['VASchedule'].is_linked and self.linkmenu in ('SO', 'DO', 'HO') else ''
        
        for sock in (self.inputs[:] + self.outputs[:]):
            for link in sock.links:
                othernode = (link.from_node, link.to_node)[sock.is_output]
                if othernode.bl_idname == 'No_En_Net_Ext' and enng['enviparams']['wpca'] == 1:
                    en = othernode.name

        if self.linkmenu == 'DO':
            cfparams = ('Name', 'Air Mass Flow Coefficient When Opening is Closed (kg/s-m)', 'Air Mass Flow Exponent When Opening is Closed (dimensionless)',
                       'Type of Rectanguler Large Vertical Opening (LVO)', 'Extra Crack Length or Height of Pivoting Axis (m)', 'Number of Sets of Opening Factor Data',
                        'Opening Factor 1 (dimensionless)', 'Discharge Coefficient for Opening Factor 1 (dimensionless)', 'Width Factor for Opening Factor 1 (dimensionless)',
                        'Height Factor for Opening Factor 1 (dimensionless)', 'Start Height Factor for Opening Factor 1 (dimensionless)', 'Opening Factor 2 (dimensionless)',
                        'Discharge Coefficient for Opening Factor 2 (dimensionless)', 'Width Factor for Opening Factor 2 (dimensionless)', 'Height Factor for Opening Factor 2 (dimensionless)',
                        'Start Height Factor for Opening Factor 2 (dimensionless)', 'Opening Factor 3 (dimensionless)', 'Discharge Coefficient for Opening Factor 3 (dimensionless)',
                        'Width Factor for Opening Factor 3 (dimensionless)', 'Height Factor for Opening Factor 3 (dimensionless)', 'Start Height Factor for Opening Factor 3 (dimensionless)',
                        'Opening Factor 4 (dimensionless)', 'Discharge Coefficient for Opening Factor 4 (dimensionless)', 'Width Factor for Opening Factor 4 (dimensionless)',
                        'Height Factor for Opening Factor 4 (dimensionless)', 'Start Height Factor for Opening Factor 4 (dimensionless)')
            cfparamsv = ('{}_{}'.format(self.name, self.linkmenu), self.amfcc, self.amfec, self.lvo, self.ecl, self.noof, '{:.3f}'.format(self.of1), self.dcof1,self.wfof1, self.hfof1, self.sfof1,
                         self.of2, self.dcof2,self.wfof2, self.hfof2, self.sfof2, self.of3, self.dcof3,self.wfof3, self.hfof3, self.sfof3, self.of4, self.dcof4,self.wfof4, self.hfof4, self.sfof4)

        elif self.linkmenu == 'SO':
            cfparams = ('Name', 'Air Mass Flow Coefficient When Opening is Closed (kg/s-m)', 'Air Mass Flow Exponent When Opening is Closed (dimensionless)', 'Minimum Density Difference for Two-Way Flow (kg/m3)', 'Discharge Coefficient (dimensionless)')
            cfparamsv = ('{}_{}'.format(self.name, self.linkmenu), '{:.5f}'.format(self.amfcc), '{:.3f}'.format(self.amfec), '{:.3f}'.format(self.ddtw), '{:.3f}'.format(self.dcof))

        elif self.linkmenu == 'HO':
            if not (self.inputs['Node 1'].is_linked or self.inputs['Node 2'].is_linked and self.outputs['Node 1'].is_linked or self.outputs['Node 2'].is_linked):
                exp_op.report({'ERROR'}, 'All horizontal opening surfaces must sit on the boundary between two thermal zones')

            cfparams = ('Name', 'Air Mass Flow Coefficient When Opening is Closed (kg/s-m)', 'Air Mass Flow Exponent When Opening is Closed (dimensionless)', 'Sloping Plane Angle (deg)', 'Discharge Coefficient (dimensionless)')
            cfparamsv = ('{}_{}'.format(self.name, self.linkmenu), '{:.5f}'.format(self.amfcc), '{:.2f}'.format(self.amfec), '{:.2f}'.format(self.spa), '{:.2f}'.format(self.dcof))

        elif self.linkmenu == 'ELA':
            cfparams = ('Name', 'Effective Leakage Area (m2)', 'Discharge Coefficient (dimensionless)', 'Reference Pressure Difference (Pa)', 'Air Mass Flow Exponent (dimensionless)')
            cfparamsv = ('{}_{}'.format(self.name, self.linkmenu), '{:.5f}'.format(self['ela']), '{:.2f}'.format(self.dcof), '{:.1f}'.format(self.rpd), '{:.3f}'.format(self.amfe))

        elif self.linkmenu == 'Crack':
            crname = 'ReferenceCrackConditions' if enng['enviparams']['crref'] == 1 else ''
            cfparams = ('Name', 'Air Mass Flow Coefficient at Reference Conditions (kg/s)', 'Air Mass Flow Exponent (dimensionless)', 'Reference Crack Conditions')
            cfparamsv = ('{}_{}'.format(self.name, self.linkmenu), self.amfc, self.amfe, crname)

        cftypedict = {'DO':'Component:DetailedOpening', 'SO':'Component:SimpleOpening', 'HO':'Component:HorizontalOpening', 'Crack':'Surface:Crack', 'ELA':'Surface:EffectiveLeakageArea'}
        cfentry = epentry('AirflowNetwork:MultiZone:{}'.format(cftypedict[self.linkmenu]), cfparams, cfparamsv)

        for sock in (self.inputs[:] + self.outputs[:]):
            for link in sock.links:
                othersock = (link.from_socket, link.to_socket)[sock.is_output]
                othernode = (link.from_node, link.to_node)[sock.is_output]
                
                if sock.bl_idname == 'So_En_Net_SSFlow' and othernode.bl_idname == 'No_En_Net_Zone':
                    # The conditional below checks if the airflow surface is also on a boundary. If so only the surface belonging to the outputting zone node is written.
                    if (othersock.name[0:-2]+'b' in [s.name for s in othernode.outputs] and othernode.outputs[othersock.name[0:-2]+'b'].links) or othersock.name[0:-2]+'b' not in [s.name for s in othernode.outputs]:
                        zn = othernode.zone
                        sn = othersock.sn
                        snames.append(('win-', 'door-')[get_con_node(bpy.data.materials[othersock.name[:-len(sn)-4]].vi_params).envi_con_type == 'Door']+zn+'_'+sn)
                        params = ('Surface Name', 'Leakage Component Name', 'External Node Name', 'Window/Door Opening Factor')
                        paramvs = (snames[-1], '{}_{}'.format(self.name, self.linkmenu), en, self.wdof1)
                        if self.linkmenu in ('SO', 'DO'):
                            params += ('Ventilation Control Mode', 'Vent Temperature Schedule Name', 'Limit  Value on Multiplier for Modulating Venting Open Factor (dimensionless)', \
                            'Lower Value on Inside/Outside Temperature Difference for Modulating the Venting Open Factor (deltaC)', 'Upper Value on Inside/Outside Temperature Difference for Modulating the Venting Open Factor (deltaC)',\
                            'Lower Value on Inside/Outside Enthalpy Difference for Modulating the Venting Open Factor (J/kg)', 'Upper Value on Inside/Outside Enthalpy Difference for Modulating the Venting Open Factor (J/kg)', 'Venting Availability Schedule Name')
                            paramvs += (self.controls, tspsname, '{:.2f}'.format(self.mvof), self.lvof, self.uvof, '', '', vasname)
                        elif self.linkmenu =='HO':
                            params += ('Ventilation Control Mode', 'Ventilation Control Zone Temperature Setpoint Schedule Name')
                            paramvs += (self.controls, '')
                        surfentry += epentry('AirflowNetwork:MultiZone:Surface', params, paramvs)
        self['sname'] = snames
        self.legal()
        return surfentry + cfentry

    def legal(self):
        nodecolour(self, 1) if (self.controls == 'Temperature' and not self.inputs['TSPSchedule'].is_linked) or (self.id_data['enviparams']['wpca'] and not self.extnode) else nodecolour(self, 0)
        for sock in self.inputs[:] + self.outputs[:]:
            sock.hide = sock.hide
        self.id_data.interface_update(bpy.context)

class No_En_Net_Ext(Node, EnViNodes):
    '''Node describing an EnVi external node'''
    bl_idname = 'No_En_Net_Ext'
    bl_label = 'Envi External Node'
    bl_icon = 'SOUND'

    height: FloatProperty(default = 1.0)
    (wpc1, wpc2, wpc3, wpc4, wpc5, wpc6, wpc7, wpc8, wpc9, wpc10, wpc11, wpc12) = [FloatProperty(name = '', default = 0, min = -1, max = 1) for x in range(12)]
    enname: StringProperty()

    def init(self, context):
        self.inputs.new('So_En_Net_SSFlow', 'Sub surface')
        self.inputs.new('So_En_Net_SFlow', 'Surface')
        self.outputs.new('So_En_Net_SSFlow', 'Sub surface')
        self.outputs.new('So_En_Net_SFlow', 'Surface')

    def draw_buttons(self, context, layout):
        layout.prop(self, 'height')
        row= layout.row()
        row.label(text = 'WPC Values')
        
        for w in range(1, 13):
            row = layout.row()
            row.prop(self, 'wpc{}'.format(w))

    def update(self):
        for sock in self.outputs:
            socklink(sock, self.id_data.name)
        sockhide(self, ('Sub surface', 'Surface'))

    def epwrite(self, enng):
        enentry, wpcname, wpcentry = '', '', ''
        for sock in self.inputs[:] + self.outputs[:]:
            for link in sock.links:
                wpcname = self.name+'_wpcvals'
                wpcs = (self.wpc1, self.wpc2, self.wpc3, self.wpc4, self.wpc5, self.wpc6, self.wpc7, self.wpc8, self.wpc9, self.wpc10, self.wpc11, self.wpc12)
                wparams = ['Name', 'AirflowNetwork:MultiZone:WindPressureCoefficientArray Name'] + ['Wind Pressure Coefficient Value {} (dimensionless)'.format(w + 1) for w in range(enng['enviparams']['wpcn'])]
                wparamvs =  ['{}_wpcvals'.format(self.name), 'WPC Array'] + [wpcs[wp] for wp in range(len(wparams))]
                wpcentry = epentry('AirflowNetwork:MultiZone:WindPressureCoefficientValues', wparams, wparamvs)
                params = ['Name', 'External Node Height (m)', 'Wind Pressure Coefficient Values Object Name']
                paramvs = [self.name, self.height, wpcname]
                enentry = epentry('AirflowNetwork:MultiZone:ExternalNode', params, paramvs)
        return enentry + wpcentry

class No_En_Net_SFlow(Node, EnViNodes):
    '''Surface airflow node'''
    bl_idname = 'No_En_Net_SFlow'
    bl_label = 'Envi surface flow'
    bl_icon = 'SOUND'

    linktype = [("Crack", "Crack", "Crack aperture used for leakage calculation"),
        ("ELA", "ELA", "Effective leakage area")]

    linkmenu: EnumProperty(name="Type", description="Linkage type", items=linktype, default='ELA')
    of: FloatProperty(default = 0.1, min = 0.001, max = 1, name = "", description = 'Opening Factor 1 (dimensionless)')
    ecl: FloatProperty(default = 0.0, min = 0, name = '', description = 'Extra Crack Length or Height of Pivoting Axis (m)')
    dcof: FloatProperty(default = 0.7, min = 0, max = 1, name = '', description = 'Discharge Coefficient')
    amfc: FloatProperty(min = 0.001, max = 1, default = 0.01, name = "", precision = 5, description = 'Flow coefficient')
    amfe: FloatProperty(min = 0.5, max = 1, default = 0.65, name = "", precision = 3, description = 'Flow exponent')
    dlen: FloatProperty(default = 2, name = "")
    dhyd: FloatProperty(default = 0.1, name = "")
    dcs: FloatProperty(default = 0.1, name = "")
    dsr: FloatProperty(default = 0.0009, name = "")
    dlc: FloatProperty(default = 1.0, name = "")
    dhtc: FloatProperty(default = 0.772, name = "")
    dmtc: FloatProperty(default = 0.0001, name = "")
    cf: FloatProperty(default = 1, min = 0, max = 1, name = "")
    rpd: FloatProperty(default = 4, min = 0.1, max = 50, name = "")
    fe: FloatProperty(default = 4, min = 0.1, max = 1, name = "", description = 'Fan Efficiency')
    pr: IntProperty(default = 500, min = 1, max = 10000, name = "", description = 'Fan Pressure Rise')
    mf: FloatProperty(default = 0.1, min = 0.001, max = 5, name = "", description = 'Maximum Fan Flow Rate (m3/s)')
    extnode: BoolProperty(default = 0)

    def init(self, context):
        self['ela'] = 1.0
        self.inputs.new('So_En_Net_SFlow', 'Node 1')
        self.inputs.new('So_En_Net_SFlow', 'Node 2')
        self.outputs.new('So_En_Net_SFlow', 'Node 1')
        self.outputs.new('So_En_Net_SFlow', 'Node 2')

    def update(self):
        for sock in self.outputs:
            socklink(sock, self.id_data.name)
        if self.linkmenu == 'ELA':
            retelaarea(self)
        self.extnode = 0
        for sock in self.inputs[:] + self.outputs[:]:
            for l in sock.links:
                if (l.from_node, l.to_node)[sock.is_output].bl_idname == 'No_En_Net_Ext':
                    self.extnode = 1
        if self.outputs.get('Node 2'):
            sockhide(self, ('Node 1', 'Node 2'))
        self.legal()

    def draw_buttons(self, context, layout):
        layout.prop(self, 'linkmenu')
        layoutdict = {'Crack':(('Coefficient', 'amfc'), ('Exponent', 'amfe'), ('Factor', 'of')), 'ELA':(('ELA (m^2)', '["ela"]'), ('DC', 'dcof'), ('PA diff (Pa)', 'rpd'), ('FE', 'amfe')),
        'EF':(('Off FC', 'amfc'), ('Off FE', 'amfe'), ('Efficiency', 'fe'), ('PA rise (Pa)', 'pr'), ('Max flow', 'mf'))}

        for vals in layoutdict[self.linkmenu]:
            newrow(layout, '{}:'.format(vals[0]), self, vals[1])

    def epwrite(self, exp_op, enng):
        fentry, crentry, zn, en, surfentry, crname, snames = '', '', '', '', '', '', []
        for sock in (self.inputs[:] + self.outputs[:]):
            for link in sock.links:
                othernode = (link.from_node, link.to_node)[sock.is_output]
                if othernode.bl_idname == 'No_En_Net_Ext' and enng['enviparams']['wpca'] == 1:
                    en = othernode.name

        if self.linkmenu == 'ELA':
            cfparams = ('Name', 'Effective Leakage Area (m2)', 'Discharge Coefficient (dimensionless)', 'Reference Pressure Difference (Pa)', 'Air Mass Flow Exponent (dimensionless)')
            cfparamvs = ('{}_{}'.format(self.name, self.linkmenu), '{:.5f}'.format(self['ela']), self.dcof, self.rpd, '{:.5f}'.format(self.amfe))

        elif self.linkmenu == 'Crack':
            crname = 'ReferenceCrackConditions' if enng['enviparams']['crref'] == 1 else ''
            cfparams = ('Name', 'Air Mass Flow Coefficient at Reference Conditions (kg/s)', 'Air Mass Flow Exponent (dimensionless)', 'Reference Crack Conditions')
            cfparamvs = ('{}_{}'.format(self.name, self.linkmenu), '{:.5f}'.format(self.amfc), '{:.5f}'.format(self.amfe), crname)

        elif self.linkmenu == 'EF':
            cfparams = ('Name', 'Air Mass Flow Coefficient When the Zone Exhaust Fan is Off at Reference Conditions (kg/s)', 'Air Mass Flow Exponent When the Zone Exhaust Fan is Off (dimensionless)')
            cfparamvs = ('{}_{}'.format(self.name, self.linkmenu), self.amfc, self.amfe)
            schedname = self.inputs['Fan Schedule'].links[0].from_node.name if self.inputs['Fan Schedule'].is_linked else ''

            for sock in [inp for inp in self.inputs if 'Node' in inp.name and inp.is_linked] + [outp for outp in self.outputs if 'Node' in outp.name and outp.is_linked]:
                zname = (sock.links[0].from_node, sock.links[0].to_node)[sock.is_output].zone

            fparams = ('Name', 'Availability Schedule Name', 'Fan Efficiency', 'Pressure Rise (Pa)', 'Maximum Flow Rate (m3/s)', 'Air Inlet Node Name', 'Air Outlet Node Name', 'End-Use Subcategory')
            fparamvs = ('{}_{}'.format(self.name,  self.linkmenu), schedname, self.fe, self.pr, self.mf, '{} Exhaust Node'.format(zname), '{} Exhaust Fan Outlet Node'.format(zname), '{} Exhaust'.format(zname))
            fentry = epentry('Fan:ZoneExhaust', fparams, fparamvs)

        cftypedict = {'Crack':'Surface:Crack', 'ELA':'Surface:EffectiveLeakageArea', 'EF': 'Component:ZoneExhaustFan'}
        cfentry = epentry('AirflowNetwork:MultiZone:{}'.format(cftypedict[self.linkmenu]), cfparams, cfparamvs)

        for sock in self.inputs[:] + self.outputs[:]:
            for link in sock.links:
                othersock = (link.from_socket, link.to_socket)[sock.is_output]
                othernode = (link.from_node, link.to_node)[sock.is_output]
                if sock.bl_idname == 'So_En_Net_SFlow' and othernode.bl_idname == 'No_En_Net_Zone':
                    # The conditional below checks if the airflow surface is also on a boundary. If so only the surface belonging to the outputting zone node is written.
                    if (othersock.name[0:-1]+'b' in [s.name for s in othernode.outputs[:]] and othernode.outputs[othersock.name[0:-1]+'b'].links) or othersock.name[0:-1]+'b' not in [s.name for s in othernode.outputs]:                        
                        sn = othersock.sn
                        zn = othernode.zone
                        snames.append(zn+'_'+sn)
                        params = ('Surface Name', 'Leakage Component Name', 'External Node Name', 'Window/Door Opening Factor, or Crack Factor (dimensionless)')
                        paramvs = (snames[-1], '{}_{}'.format(self.name, self.linkmenu), en, '{:.5f}'.format(self.of))
                        surfentry += epentry('AirflowNetwork:MultiZone:Surface', params, paramvs)

        self['sname'] = snames
        self.legal()
        return surfentry + cfentry + crentry + fentry

    def legal(self):
        try:
            nodecolour(self, 1) if not self.extnode and self.id_data['enviparams']['wpca'] else nodecolour(self, 0)
            self.id_data.interface_update(bpy.context)
        except:
            nodecolour(self, 1)

class No_En_Net_ACon(Node, EnViNodes):
    '''Node defining the overall airflow network simulation'''
    bl_idname = 'No_En_Net_ACon'
    bl_label = 'AFN Control'
    bl_icon = 'SOUND'

    def wpcupdate(self, context):
        if self.wpctype == 'SurfaceAverageCalculation':
            if self.inputs['WPC Array'].is_linked:
                remlink(self, self.inputs['WPC Array'].links)
            self.inputs['WPC Array'].hide = True
        elif self.wpctype == 'Input':
            self.inputs['WPC Array'].hide = False
        self.legal()

    afnname: StringProperty(name = '')
    afntype: EnumProperty(items = [('MultizoneWithDistribution', 'MultizoneWithDistribution', 'Include a forced airflow system in the model'),
                                              ('MultizoneWithoutDistribution', 'MultizoneWithoutDistribution', 'Exclude a forced airflow system in the model'),
                                              ('MultizoneWithDistributionOnlyDuringFanOperation', 'MultizoneWithDistributionOnlyDuringFanOperation', 'Apply forced air system only when in operation'),
                                              ('NoMultizoneOrDistribution', 'NoMultizoneOrDistribution', 'Only zone infiltration controls are modelled')], name = "", default = 'MultizoneWithoutDistribution')

    wpctype: EnumProperty(items = [('SurfaceAverageCalculation', 'SurfaceAverageCalculation', 'Calculate wind pressure coefficients based on oblong building assumption'),
                                              ('Input', 'Input', 'Input wind pressure coefficients from an external source')], name = "", default = 'SurfaceAverageCalculation', update = wpcupdate)
    wpcaname: StringProperty()
    wpchs: EnumProperty(items = [('OpeningHeight', 'OpeningHeight', 'Calculate wind pressure coefficients based on opening height'),
                                              ('ExternalNode', 'ExternalNode', 'Calculate wind pressure coefficients based on external node height')], name = "", default = 'OpeningHeight')
    buildtype: EnumProperty(items = [('LowRise', 'Low Rise', 'Height is less than 3x the longest wall'),
                                              ('HighRise', 'High Rise', 'Height is more than 3x the longest wall')], name = "", default = 'LowRise')

    maxiter: IntProperty(default = 500, description = 'Maximum Number of Iterations', name = "")

    initmet: EnumProperty(items = [('ZeroNodePressures', 'ZeroNodePressures', 'Initilisation type'),
                                              ('LinearInitializationMethod', 'LinearInitializationMethod', 'Initilisation type')], name = "", default = 'ZeroNodePressures')
    rcontol: FloatProperty(default = 0.0001, description = 'Relative Airflow Convergence Tolerance', name = "")
    acontol: FloatProperty(min = 0.000001, max = 0.1, default = 0.000001, description = 'Absolute Airflow Convergence Tolerance', name = "")
    conal: FloatProperty(default = -0.1, max = 1, min = -1, description = 'Convergence Acceleration Limit', name = "")
    aalax: IntProperty(default = 0, max = 180, min = 0, description = 'Azimuth Angle of Long Axis of Building', name = "")
    rsala: FloatProperty(default = 1, max = 1, min = 0, description = 'Ratio of Building Width Along Short Axis to Width Along Long Axis', name = "")

    def init(self, context):
        self.inputs.new('So_En_Net_WPC', 'WPC Array')

    def draw_buttons(self, context, layout):
        yesno = (1, 1, 1, self.wpctype == 'Input', self.wpctype != 'Input' and self.wpctype == 'SurfaceAverageCalculation', 1, 1, 1, 1, 1, self.wpctype == 'SurfaceAverageCalculation', self.wpctype == 'SurfaceAverageCalculation')
        vals = (('Name:', 'afnname'), ('Type:', 'afntype'), ('WPC type:', 'wpctype'), ('WPC height', 'wpchs'), ('Build type:', 'buildtype'), ('Max iter:','maxiter'), ('Init method:', 'initmet'),
         ('Rel Converge:', 'rcontol'), ('Abs Converge:', 'acontol'), ('Converge Lim:', 'conal'), ('Azimuth:', 'aalax'), ('Axis ratio:', 'rsala'))
        [newrow(layout, val[0], self, val[1]) for v, val in enumerate(vals) if yesno[v]]

    def epwrite(self, exp_op, enng):
        wpcaentry = ''
        if self.wpctype == 'Input' and not self.inputs['WPC Array'].is_linked:
            exp_op.report({'ERROR'},"WPC array input has been selected in the control node, but no WPC array node is attached")
            return 'ERROR'

#        wpcaname = 'WPC Array' if not self.wpcaname else self.wpcaname
        self.afnname = 'default' if not self.afnname else self.afnname
        wpctype = 1 if self.wpctype == 'Input' else 0
        paramvs = (self.afnname, self.afntype,
                     self.wpctype, ("", self.wpchs)[wpctype], (self.buildtype, "")[wpctype], self.maxiter, self.initmet,
                    '{:.3E}'.format(self.rcontol), '{:.3E}'.format(self.acontol), '{:.3E}'.format(self.conal), (self.aalax, "")[wpctype], (self.rsala, "")[wpctype])

        params = ('Name', 'AirflowNetwork Control', 'Wind Pressure Coefficient Type', \
        'Height Selection for Local Wind Pressure Calculation', 'Building Type', 'Maximum Number of Iterations (dimensionless)', 'Initialization Type', \
        'Relative Airflow Convergence Tolerance (dimensionless)', 'Absolute Airflow Convergence Tolerance (kg/s)', 'Convergence Acceleration Limit (dimensionless)', \
        'Azimuth Angle of Long Axis of Building (deg)', 'Ratio of Building Width Along Short Axis to Width Along Long Axis')

        simentry = epentry('AirflowNetwork:SimulationControl', params, paramvs)

        if self.inputs['WPC Array'].is_linked:
            (wpcaentry, enng['enviparams']['wpcn']) = self.inputs['WPC Array'].links[0].from_node.epwrite() if wpctype == 1 else ('', 0)
            enng['enviparams']['wpca'] = 1
        self.legal()
        return simentry + wpcaentry

    def update(self):
        self.legal()

    def legal(self):
        try:
            bpy.data.node_groups[self['nodeid'].split('@')[1]]['enviparams']['wpca'] = 1 if self.wpctype == 'Input' and self.inputs['WPC Array'].is_linked else 0
            nodecolour(self, self.wpctype == 'Input' and not self.inputs['WPC Array'].is_linked)
            for node in [node for node in bpy.data.node_groups[self['nodeid'].split('@')[1]].nodes if node.bl_idname in ('EnViSFlow', 'EnViSSFlow')]:
                node.legal()
        except:
            pass
        
class No_En_Net_Sched(Node, EnViNodes):
    '''Node describing a schedule'''
    bl_idname = 'No_En_Net_Sched'
    bl_label = 'Schedule'
    bl_icon = 'SOUND'

    def tupdate(self, context):
        try:
            err = 0
            if self.t2 <= self.t1 and self.t1 < 365:
                self.t2 = self.t1 + 1
                if self.t3 <= self.t2 and self.t2 < 365:
                    self.t3 = self.t2 + 1
                    if self.t4 != 365:
                        self.t4 = 365

            tn = (self.t1, self.t2, self.t3, self.t4).index(365) + 1
            if max((self.t1, self.t2, self.t3, self.t4)[:tn]) != 365:
                err = 1
            if any([not f for f in (self.f1, self.f2, self.f3, self.f4)[:tn]]):
                err = 1
            if any([not u or '. ' in u or len(u.split(';')) != len((self.f1, self.f2, self.f3, self.f4)[i].split(' ')) for i, u in enumerate((self.u1, self.u2, self.u3, self.u4)[:tn])]):
                err = 1

            for f in (self.f1, self.f2, self.f3, self.f4)[:tn]:
                for fd in f.split(' '):
                    if not fd or (fd and fd.upper() not in ("ALLDAYS", "WEEKDAYS", "WEEKENDS", "MONDAY", "TUESDAY", "WEDNESDAY", "THURSDAY", "FRIDAY", "SATURDAY", "SUNDAY", "ALLOTHERDAYS")):
                        err = 1

            for u in (self.u1, self.u2, self.u3, self.u4)[:tn]:
                for uf in u.split(';'):
                    for ud in uf.split(','):
                        if len(ud.split()[0].split(':')) != 2 or int(ud.split()[0].split(':')[0]) not in range(1, 25) or len(ud.split()[0].split(':')) != 2 or not ud.split()[0].split(':')[1].isdigit() or int(ud.split()[0].split(':')[1]) not in range(0, 60):
                            err = 1
            nodecolour(self, err)

        except:
            nodecolour(self, 1)

    file: EnumProperty(name = '', items = [("0", "None", "No file"), ("1", "Select", "Select file"), ("2", "Generate", "Generate file")], default = '0')
    select_file: StringProperty(name="", description="Name of the variable file", default="", subtype="FILE_PATH")
    cn: IntProperty(name = "", default = 1, min = 1)
    rtsat: IntProperty(name = "", default = 0, min = 0)
    hours: IntProperty(name = "", default = 8760, min = 1, max = 8760)
    delim: EnumProperty(name = '', items = [("Comma", "Comma", "Comma delimiter"), ("Space", "Space", "space delimiter")], default = 'Comma')
    generate_file: StringProperty(default = "", name = "")
    (u1, u2, u3, u4) =  [StringProperty(name = "", description = "Valid entries (; separated for each 'For', comma separated for each day, space separated for each time value pair)", update = tupdate)] * 4
    (f1, f2, f3, f4) =  [StringProperty(name = "", description = "Valid entries (space separated): AllDays, Weekdays, Weekends, Monday, Tuesday, Wednesday, Thursday, Friday, Saturday, Sunday, AllOtherDays", update = tupdate)] * 4
    (t1, t2, t3, t4) = [IntProperty(name = "", default = 365, min = 1, max = 365, update = tupdate)] * 4

    def init(self, context):
        self.outputs.new('So_En_Sched', 'Schedule')
        self['scheddict'] = {'TSPSchedule': 'Any Number', 'VASchedule': 'Fraction', 'Fan Schedule': 'Fraction', 'HSchedule': 'Temperature', 'CSchedule': 'Temperature'}
        self.tupdate(context)
        nodecolour(self, 1)

    def draw_buttons(self, context, layout):
        uvals, u = (1, self.u1, self.u2, self.u3, self.u4), 0
        tvals = (0, self.t1, self.t2, self.t3, self.t4)
        newrow(layout, 'From file', self, 'file')
        
        if self.file == "1":
            newrow(layout, 'Select', self, 'select_file')
            newrow(layout, 'Columns', self, 'cn')
            newrow(layout, 'Skip rows', self, 'rtsat')
            newrow(layout, 'Delimiter', self, 'delim')
        elif self.file == "2":
            newrow(layout, 'Generate', self, 'generate_file')

        if self.file != "1":        
            while uvals[u] and tvals[u] < 365:
                [newrow(layout, v[0], self, v[1]) for v in (('End day {}:'.format(u+1), 't'+str(u+1)), ('Fors:', 'f'+str(u+1)), ('Untils:', 'u'+str(u+1)))]
                u += 1

    def update(self):
        for sock in self.outputs:
            socklink(sock, self.id_data.name)
        self.id_data.interface_update(bpy.context)

    def epwrite(self, name, stype):
        schedtext, ths = '', []
        for tosock in [link.to_socket for link in self.outputs['Schedule'].links]:
            if not schedtext:
                for t in (self.t1, self.t2, self.t3, self.t4):
                    ths.append(t)
                    if t == 365:
                        break
#                ths = [self.t1, self.t2, self.t3, self.t4]
                fos = [fs for fs in (self.f1, self.f2, self.f3, self.f4) if fs]
                uns = [us for us in (self.u1, self.u2, self.u3, self.u4) if us]
                ts, fs, us = rettimes(ths, fos, uns)
                
#                if self.file == '0':
                schedtext = epschedwrite(name, stype, ts, fs, us)
        return schedtext
            
    def epwrite_sel_file(self, name):               
        params = ('Name', 'ScheduleType', 'Name of File', 'Column Number', 'Rows to Skip at Top', 'Number of Hours of Data', 'Column Separator')
        paramvs = (name, 'Any number', os.path.abspath(self.select_file), self.cn, self.rtsat, 8760, self.delim) 
        schedtext = epentry('Schedule:File', params, paramvs)
        '''    Schedule:File,
        elecTDVfromCZ01res, !- Name
        Any Number, !- ScheduleType
        TDV_kBtu_CTZ01.csv, !- Name of File
        2, !- Column Number
        4, !- Rows to Skip at Top
        8760, !- Number of Hours of Data
        Comma; !- Column Separator'''
        return schedtext
    
    def epwrite_gen_file(self, name, data, newdir):
        schedtext, ths = '', []
        for tosock in [link.to_socket for link in self.outputs['Schedule'].links]:
            if not schedtext:
                for t in (self.t1, self.t2, self.t3, self.t4):
                    ths.append(t)
                    if t == 365:
                        break
#                ths = [self.t1, self.t2, self.t3, self.t4]
                fos = [fs for fs in (self.f1, self.f2, self.f3, self.f4) if fs]
                uns = [us for us in (self.u1, self.u2, self.u3, self.u4) if us]
                ts, fs, us = rettimes(ths, fos, uns)
        for t in ts:
            for f in fs:
                for u in us:
                    for hi, h in enumerate((datetime.datetime(2015, 1, 1, 0, 00) - datetime.datetime(2014, 1, 1, 0, 00)).hours):
                        if h.day <= self.ts:
#                            if f == 'Weekday'
                            data[hi] = 1
        with open(os.path.join(newdir, name), 'w') as sched_file:
            sched_file.write(',\n'.join([d for d in data]))
        params = ('Name', 'ScheduleType', 'Name of File', 'Column Number', 'Rows to Skip at Top', 'Number of Hours of Data', 'Column Separator')
        paramvs = (name, 'Any number', os.path.abspath(self.select_file), self.cn, self.rtsat, 8760, self.delim) 
        schedtext = epentry('Schedule:File', params, paramvs)    
        return schedtext

class No_En_Net_Prog(Node, EnViNodes):
    '''Node describing an EMS Program'''
    bl_idname = 'No_En_Net_Prog'
    bl_label = 'Envi Program'
    bl_icon = 'SOUND'

    text_file: StringProperty(description="Textfile to show")

    def init(self, context):
        self.outputs.new('So_En_Net_Sense', 'Sensor')
        self.outputs.new('So_En_Net_Act', 'Actuator')
        nodecolour(self, 1)

    def draw_buttons(self, context, layout):
        layout.prop_search(self, 'text_file', bpy.data, 'texts', text='File', icon='TEXT')

    def update(self):
        for sock in self.outputs:
            socklink2(sock, self.id_data.name)
        nodecolour(self, not all([sock.links for sock in self.outputs]) and any([sock.links for sock in self.outputs]))

    def epwrite(self):
        if not (self.outputs['Sensor'].links and self.outputs['Actuator'].links):
            return ''
        else:
            sentries, aentries = '', ''
            for slink in self.outputs['Sensor'].links:
                snode = slink.to_node
                sparams = ('Name', 'Output:Variable or Output:Meter Index Key Name', 'EnergyManagementSystem:Sensor')
                
                if snode.bl_idname == 'No_En_Net_EMSZone':
                    sparamvs = ('{}_{}'.format(snode.emszone, snode.sensordict[snode.sensortype][0]), '{}'.format(snode.emszone), snode.sensordict[snode.sensortype][1])
    
                elif snode.bl_label == 'No_En_Net_Occ':
                    for zlink in snode.outputs['Occupancy'].links:
                        znode = zlink.to_node
                        sparamvs = ('{}_{}'.format(znode.zone, snode.sensordict[snode.sensortype][0]), '{}'.format(znode.zone), snode.sensordict[snode.sensortype][1])
                sentries += epentry('EnergyManagementSystem:Sensor', sparams, sparamvs)
    
            for alink in self.outputs['Actuator'].links:
                anode, asocket = alink.to_node, alink.to_socket
                aparams = ('Name', 'Actuated Component Unique Name', 'Actuated Component Type', 'Actuated Component Control Type')
                aparamvs = (asocket.name, asocket.sn, anode.compdict[anode.acttype], anode.actdict[anode.acttype][0])
                aentries += epentry('EnergyManagementSystem:Actuator', aparams, aparamvs)
    
            cmparams = ('Name', 'EnergyPlus Model Calling Point', 'Program Name 1')
            cmparamvs = (self.name.replace(' ', '_'), 'BeginTimestepBeforePredictor', '{}_controller'.format(self.name.replace(' ', '_')))
            cmentry = epentry('EnergyManagementSystem:ProgramCallingManager', cmparams, cmparamvs)
            pparams = ['Name'] + ['line{}'.format(l) for l, line in enumerate(bpy.data.texts[self.text_file].lines) if line.body and line.body.strip()[0] != '!']
            pparamvs = ['{}_controller'.format(self.name.replace(' ', '_'))] + [line.body.strip() for line in bpy.data.texts[self.text_file].lines if line.body and line.body.strip()[0] != '!']
            pentry = epentry('EnergyManagementSystem:Program', pparams, pparamvs)
            return sentries + aentries + cmentry + pentry

class No_En_Net_EMSZone(Node, EnViNodes):
    '''Node describing a simulation zone'''
    bl_idname = 'No_En_Net_EMSZone'
    bl_label = 'EMS Zone'
    bl_icon = 'SOUND'

    def supdate(self, context):
        self.inputs[0].name = '{}_{}'.format(self.emszone, self.sensordict[self.sensortype][0])

    def zupdate(self, context):
        adict = {'Window': 'win', 'Door': 'door'}
        self.supdate(context)
        sssocklist = []
        
        try:            
            obj = bpy.data.collections[self.emszone].objects[0]
            odm = [ms.material for ms in obj.material_slots]
            
            for face in obj.data.polygons:
                mat = odm[face.material_index]
                
                for emnode in mat.vi_params.envi_nodes.nodes:
                    if emnode.bl_idname == 'No_En_Mat_Con' and emnode.active and emnode.envi_afsurface and emnode.envi_con_type in ('Window', 'Door'):
                        sssocklist.append('{}_{}_{}_{}'.format(adict[emnode.envi_con_type], self.emszone, face.index, self.actdict[self.acttype][1]))         
            
            self.inputs[0].hide = False
            nodecolour(self, 0)
        except:
            self.inputs[0].hide = True
            nodecolour(self, 1)

        for iname in [inputs for inputs in self.inputs if inputs.name not in sssocklist and inputs.bl_idname == 'So_En_Net_Act']:
            try: self.inputs.remove(iname)
            except: pass

        for sock in sorted(set(sssocklist)):
            if not self.inputs.get(sock):
                try: 
                    self.inputs.new('So_En_Net_Act', sock).sn = sock.split('_')[0] + '-' + '_'.join(sock.split('_')[1:-1])
                except Exception as e: print('3190', e)

    emszone: StringProperty(name = '', update = zupdate)
    sensorlist = [("0", "Zone Temperature", "Sense the zone temperature"), ("1", "Zone Humidity", "Sense the zone humidity"), ("2", "Zone CO2", "Sense the zone CO2"),
                  ("3", "Zone Occupancy", "Sense the zone occupancy"), ("4", "Zone Equipment", "Sense the equipment level")]
    sensortype: EnumProperty(name="", description="Linkage type", items=sensorlist, default='0', update = supdate)
    sensordict = {'0':  ('Temp', 'Zone Mean Air Temperature'), '1': ('RH', 'Zone Air Relative Humidity'), 
                  '2': ('CO2', 'AFN Node CO2 Concentration'), '3': ('Occ', 'Zone Occupancy'), '4': ('Equip', 'Zone Equipment')}
    actlist = [("0", "Opening factor", "Actuate the opening factor"), ("1", "Air supply temp", "Actuate an ideal air load system supply temperature"),
               ("2", "Air supply flow", "Actuate an ideal air load system flow rate"), ("3", "Outdoor Air supply flow", "Actuate an ideal air load system outdoor air flow rate")]
    acttype: EnumProperty(name="", description="Actuator type", items=actlist, default='0')
    compdict = {'0': 'AirFlow Network Window/Door Opening'}
    actdict =  {'0': ('Venting Opening Factor', 'of')}

    def init(self, context):
        self.inputs.new('So_En_Net_Sense', 'Sensor')
        self.inputs[0].hide = True
        self.inputs[0].display_shape = 'SQUARE'
        nodecolour(self, 1)

    def draw_buttons(self, context, layout):
        newrow(layout, 'Zone:', self, "emszone")
        if self.emszone in [o.name for o in bpy.data.objects]:
            newrow(layout, 'Sensor', self, 'sensortype')
        if len(self.inputs) > 1:
            newrow(layout, 'Actuator', self, 'acttype')
    
class EnViNodeCategory(NodeCategory):
    @classmethod
    def poll(cls, context):
        return context.space_data.tree_type == 'EnViN'

envi_zone = [NodeItem("No_En_Net_Zone", label="Zone"), NodeItem("No_En_Net_Occ", label="Occupancy"),
             NodeItem("No_En_Net_Hvac", label="HVAC"), NodeItem("No_En_Net_Eq", label="Equipment"),
             NodeItem("No_En_Net_Inf", label="Infiltration"), NodeItem("No_En_Net_TC", label="Thermal Chimney")]
envi_sched = [NodeItem("No_En_Net_Sched", label="Schedule")]
envi_airflow = [NodeItem("No_En_Net_SFlow", label="Surface Flow"), NodeItem("No_En_Net_SSFlow", label="Sub-surface Flow"),
                NodeItem("No_En_Net_Ext", label="External Air")]
envi_ems = [NodeItem("No_En_Net_EMSZone", label="EMS Zone"), NodeItem("No_En_Net_Prog", label="EMS Program")]

envinode_categories = [EnViNodeCategory("Zone", "Zone Nodes", items=envi_zone), 
                       EnViNodeCategory("Schedule", "Schedule Nodes", items=envi_sched),
                       EnViNodeCategory("Airflow", "Airflow Nodes", items=envi_airflow),
                       EnViNodeCategory("EMS", "EMS Nodes", items=envi_ems)]

        
#        EnViNodeCategory("Control", "Control Node", items=[NodeItem("AFNCon", label="Control Node"), NodeItem("EnViWPCA", label="WPCA Node"), NodeItem("EnViCrRef", label="Crack Reference")]),
#        EnViNodeCategory("Nodes", "Zone Nodes", items=[NodeItem("EnViZone", label="Zone Node"), NodeItem("EnViExt", label="External Node"), NodeItem("EnViOcc", label="Ocupancy Node")
#        , NodeItem("EnViEq", label="Equipment Node"), NodeItem("EnViHvac", label="HVAC Node"), NodeItem("EnViInf", label="Infiltration Node"), NodeItem("EnViTC", label="Thermal Chimney Node")]),
#        EnViNodeCategory("LinkNodes", "Airflow Link Nodes", items=[
#            NodeItem("EnViSSFlow", label="Sub-surface Flow Node"), NodeItem("EnViSFlow", label="Surface Flow Node")]),
#        EnViNodeCategory("SchedNodes", "Schedule Nodes", items=[NodeItem("EnViSched", label="Schedule")]),
#        EnViNodeCategory("EMSNodes", "EMS Nodes", items=[NodeItem("EnViProg", label="Program"), NodeItem("EnViEMSZone", label="Zone")])]


class EnViMatNodes:
    @classmethod
    def poll(cls, ntree):
        return ntree.bl_idname == 'EnViMatN'

class EnViMatNetwork(NodeTree):
    '''A node tree for the creation of EnVi materials.'''
    bl_idname = 'EnViMatN'
    bl_label = 'EnVi Material'
    bl_icon = 'IMGDISPLAY'
    nodetypes = {}
    
class EnViMatNodeCategory(NodeCategory):
    @classmethod
    def poll(cls, context):
        return context.space_data.tree_type == 'EnViMatN'

class So_En_Mat_Op(NodeSocket):
    '''EnVi opaque layer socket'''
    bl_idname = 'So_En_Mat_Op'
    bl_label = 'Opaque layer socket'
    valid = ['OLayer']

    def draw(self, context, layout, node, text):
        layout.label(text = text)

    def draw_color(self, context, node):
        return (0.65, 0.16, 0.16, 1)

    def ret_valid(self, node):
        return ['OLayer']

class So_En_Mat_Ou(NodeSocket):
    '''EnVi outer layer socket'''
    bl_idname = 'So_En_Mat_Ou'
    bl_label = 'Outer layer socket'
    valid = ['OLayer', 'Tlayer', 'ScreenLayer']

    def draw(self, context, layout, node, text):
        layout.label(text = text)

    def draw_color(self, context, node):
        if node.envi_con_type == 'Window':
            return (0.65, 0.65, 1, 1)  
        else:
            return (0.65, 0.16, 0.16, 1)
    
    def ret_valid(self, node):
        if node.envi_con_type == 'Window':
            return ['TLayer', 'ScreenLayer', 'BlindLayer', 'ShadeLayer']
        else:
            return ['OLayer']
    
class So_En_Mat_Tr(NodeSocket):
    '''EnVi transparent layer socket'''
    bl_idname = 'So_En_Mat_Tr'
    bl_label = 'Transparent layer socket'
    valid = ['TLayer']

    def draw(self, context, layout, node, text):
        layout.label(text = text)

    def draw_color(self, context, node):
        return (0.65, 0.65, 1, 1.0)

    def ret_valid(self, node):
        return ['TLayer']
    
class So_En_Mat_Fr(NodeSocket):
    '''EnVi frame socket'''
    bl_idname = 'So_En_Mat_Fr'
    bl_label = 'Window frame socket'
    valid = ['Olayer']

    def draw(self, context, layout, node, text):
        layout.label(text = text)

    def draw_color(self, context, node):
        return (1, 0.1, 1, 1.0)
    
    def ret_valid(self, node):
        return ['OLayer']
    
class So_En_Mat_Gas(NodeSocket):
    '''EnVi gas layer socket'''
    bl_idname = 'So_En_Mat_Gas'
    bl_label = 'Gas layer socket'
    valid = ['GLayer']

    def draw(self, context, layout, node, text):
        layout.label(text = text)

    def draw_color(self, context, node):
        return (1, 1, 1, 1.0)
    
    def ret_valid(self, node):
        return ['GLayer']
    
class So_En_Mat_Sh(NodeSocket):
    '''EnVi shade layer socket'''
    bl_idname = 'So_En_Mat_Sh'
    bl_label = 'Shade layer socket'
    valid = ['GLayer', 'Tlayer']

    def draw(self, context, layout, node, text):
        layout.label(text)

    def draw_color(self, context, node):
        return (0, 0, 0, 1.0)

    def ret_valid(self, node):
        return ['GLayer', 'Tlayer']
    
class So_En_Mat_Sc(NodeSocket):
    '''EnVi screen layer socket'''
    bl_idname = 'So_En_Mat_Sc'
    bl_label = 'External screen layer socket'
    valid = ['ScreenLayer']

    def draw(self, context, layout, node, text):
        layout.label(text)

    def draw_color(self, context, node):
        return (0.65, 0.65, 1, 1.0)
    
    def ret_valid(self, node):
        return ['ScreenLayer']

class So_En_Mat_Sw(NodeSocket):
    '''EnVi switchable glazing layer socket'''
    bl_idname = 'So_En_Mat_Sw'
    bl_label = 'Switchable glazing layer socket'
#    valid = ['SGLayer']

    def draw(self, context, layout, node, text):
        layout.label(text)

    def draw_color(self, context, node):
        return (0.65, 0.65, 1, 1.0)
    
    def ret_valid(self, node):
        return ['TLayer']
    
class So_En_Mat_ShC(NodeSocket):
    '''EnVi shade control socket'''
    bl_idname = 'So_En_Mat_ShC'
    bl_label = 'Shade control socket'
    valid = ['SControl']

    def draw(self, context, layout, node, text):
        layout.label(text)

    def draw_color(self, context, node):
        return (0, 0, 0, 1.0)
    
class So_En_Mat_PV(NodeSocket):
    '''EnVi Photovoltaic socket'''
    bl_idname = 'So_En_Mat_PV'
    bl_label = 'PV socket'
    valid = ['PV']

    def draw(self, context, layout, node, text):
        layout.label(text = text)

    def draw_color(self, context, node):
        return (0, 0, 0, 1.0)
    
class So_En_Mat_PVG(NodeSocket):
    '''EnVi Photovoltaic generator socket'''
    bl_idname = 'So_En_Mat_PVG'
    bl_label = 'PV Generator socket'
    valid = ['PVG']

    def draw(self, context, layout, node, text):
        layout.label(text = text)

    def draw_color(self, context, node):
        return (0.9, 0.9, 0, 1.0)
    
class No_En_Mat_Con(Node, EnViMatNodes):
    '''Node defining the EnVi material construction'''
    bl_idname = 'No_En_Mat_Con'
    bl_label = 'EnVi Construction'
    bl_icon = 'FORCE_WIND'
    
    def con_update(self, context):
        if len(self.inputs) == 3:
#            if not self.pv:
#                remlink(self, self.inputs['PV'].links)
#                self.inputs['PV'].hide = True
            if self.envi_con_type == 'Shading':
                self.inputs['Schedule'].hide = False
            else:
                for link in self.inputs['Schedule'].links:
                    self.id_data.links.remove(link)
                self.inputs['Schedule'].hide = True
                
            if self.envi_con_makeup != "1" or self.envi_con_type in ('Shading', 'None'):
                for link in self.inputs['Outer layer'].links:
                    self.id_data.links.remove(link)
                self.inputs['Outer layer'].hide = True                
            else:
                self.inputs['Outer layer'].hide = False
#            else:
#                self.inputs['Outer layer'].hide = False
#                if self.pp != '0':
#                    remlink(self.id_data, self.inputs['PV'].links)
#                    self.inputs['PV'].hide = True
#                else:
#                    self.inputs[''].hide = False
                
            [link.from_node.update() for link in self.inputs['Outer layer'].links]
            get_mat(self, 0).vi_params.envi_type = self.envi_con_type   
            self.pv_update()
            self.update()
    
#    def frame_update(self, context):
#        if self.fclass in ("0", "1"):
#            for link in self.inputs['Outer frame layer'].links:
#                self.id_data.links.remove(link)
#            self.inputs['Outer frame layer'].hide = True
#        else:
#            self.inputs['Outer frame layer'].hide = False
#
#        self.update()
    def pv_update(self):
        if (self.envi_con_type in ('Wall', 'Roof') and self.envi_con_con != 'Thermal mass') or self.envi_con_type == 'Shading':
            self.inputs['PV'].hide = False
        else:
            remlink(self.id_data, self.inputs['PV'].links)
            self.inputs['PV'].hide = True
            
    def active_update(self, context):
        if self.active:
            for node in [n for n in self.id_data.nodes if n.bl_idname == 'No_En_Mat_Con' and n != self]:
                node.active = False
                
    def bc_update(self, context):
        if self.envi_con_type in ("Wall", "Roof"):
            return [("External", "External", "External boundary"),
             ("Zone", "Zone", "Zone boundary"),
             ("Thermal mass", "Thermal mass", "Adiabatic")]
        elif self.envi_con_type in ("Door", "Window"):
            return [("External", "External", "External boundary"),
             ("Zone", "Zone", "Zone boundary")]
        elif self.envi_con_type == "Floor":
            return [("Ground", "Ground", "Ground boundary"), ("External", "External", "External boundary"),
             ("Zone", "Zone", "Zone boundary"), ("Thermal mass", "Thermal mass", "Adiabatic")]
        else:
            return [("None", "None", "None")]
        
    def uv_update(self, context):
        pstcs, resists = [], []
        
        if self.envi_con_type in ('Wall', 'Floor', 'Roof'):
            if self.envi_con_makeup == '0':
                ecs = envi_constructions()
                ems = envi_materials()
                con_layers = ecs.propdict[self.envi_con_type][self.envi_con_list]
                thicks = [0.001 * tc for tc in [self.lt0, self.lt1, self.lt2, self.lt3, 
                                                self.lt4, self.lt5, self.lt6, self.lt7, self.lt8, self.lt9][:len(con_layers)]]

                for p, psmat in enumerate(con_layers):
                    pi = 2 if psmat in ems.gas_dat else 1
                    pstcs.append(float(ems.matdat[psmat][pi]))
                    resists.append((thicks[p]/float(ems.matdat[psmat][pi]), float(ems.matdat[psmat][pi]))[ems.matdat[psmat][0] == 'Gas'])
                uv = 1/(sum(resists) + 0.12 + 0.08)
                self.uv = '{:.3f}'.format(uv)                
        else:
            self.uv = "N/A"
                                      
    con_name: StringProperty(name = "", description = "", default = '')
    envi_con_type: EnumProperty(items = [("Wall", "Wall", "Wall construction"),
                                            ("Floor", "Floor", "Ground floor construction"),
                                            ("Roof", "Roof", "Roof construction"),
                                            ("Window", "Window", "Window construction"), 
                                            ("Door", "Door", "Door construction"),
                                            ("Shading", "Shading", "Shading material"),
                                            ("None", "None", "Surface to be ignored")], 
                                            name = "", 
                                            description = "Specify the construction type", 
                                            default = "None", update = con_update)
    envi_con_makeup: EnumProperty(items = [("0", "Pre-set", "Construction pre-set"),
                                            ("1", "Layers", "Custom layers"),
                                            ("2", "Dummy", "Adiabatic")], 
                                            name = "", 
                                            description = "Pre-set construction of custom layers", 
                                            default = "0", update = con_update)
    envi_con_con: EnumProperty(items = bc_update, 
                                            name = "", 
                                            description = "Construction context", update = con_update)
    envi_simple_glazing: BoolProperty(name = "", description = "Flag to signify whether to use a EP simple glazing representation", default = False)
    envi_sg_uv: FloatProperty(name = "W/m^2.K", description = "Window U-Value", min = 0.01, max = 10, default = 2.4)
    envi_sg_shgc: FloatProperty(name = "", description = "Window Solar Heat Gain Coefficient", min = 0, max = 1, default = 0.7)
    envi_sg_vt: FloatProperty(name = "", description = "Window Visible Transmittance", min = 0, max = 1, default = 0.8)
    envi_afsurface: BoolProperty(name = "", description = "Flag to signify whether the material represents an airflow surface", default = False)
    lt0: FloatProperty(name = "mm", description = "Layer thickness (mm)", min = 0.1, default = 100, update = uv_update)
    lt1: FloatProperty(name = "mm", description = "Layer thickness (mm)", min = 0.1, default = 100, update = uv_update)
    lt2: FloatProperty(name = "mm", description = "Layer thickness (mm)", min = 0.1, default = 100, update = uv_update)
    lt3: FloatProperty(name = "mm", description = "Layer thickness (mm)", min = 0.1, default = 100, update = uv_update)
    lt4: FloatProperty(name = "mm", description = "Layer thickness (mm)", min = 0.1, default = 100, update = uv_update)
    lt5: FloatProperty(name = "mm", description = "Layer thickness (mm)", min = 0.1, default = 100, update = uv_update)
    lt6: FloatProperty(name = "mm", description = "Layer thickness (mm)", min = 0.1, default = 100, update = uv_update)
    lt7: FloatProperty(name = "mm", description = "Layer thickness (mm)", min = 0.1, default = 100, update = uv_update)
    lt8: FloatProperty(name = "mm", description = "Layer thickness (mm)", min = 0.1, default = 100, update = uv_update)
    lt9: FloatProperty(name = "mm", description = "Layer thickness (mm)", min = 0.1, default = 100, update = uv_update)
    uv: StringProperty(name = "", description = "Construction U-Value", default = "N/A")
    envi_con_list: EnumProperty(items = envi_con_list, name = "", description = "Database construction")
    active: BoolProperty(name = "", description = "Active construction", default = False, update = active_update)
    
    # Frame parameters
    fclass: EnumProperty(items = [("0", "Simple spec.", "Simple frame designation"),
                                   ("1", "Detailed spec.", "Advanced frame designation"),
                                   ("2", "Layers", "Layered frame designation")], 
                                    name = "", 
                                    description = "Window frame specification", 
                                    default = "0")#, update = frame_update)
    
    fmat: EnumProperty(items = [("0", "Wood", "Wooden frame"),
                                   ("1", "Aluminium", "Aluminium frame"),
                                   ("2", "Plastic", "uPVC frame"),
                                   ("3", "Layers", "Layered frame")], 
                                    name = "", 
                                    description = "Frame material", 
                                    default = "0", update = con_update)
    
    fthi: FloatProperty(name = "m", description = "Frame thickness", min = 0.001, max = 10, default = 0.05)
    farea: FloatProperty(name = "%", description = "Frame area percentage", min = 0.001, max = 100, default = 10)
    fw: FloatProperty(name = "m", description = "Frame Width", min = 0.0, max = 10, default = 0.2)
    fop: FloatProperty(name = "m", description = "Frame Outside Projection", min = 0.01, max = 10, default = 0.1)
    fip: FloatProperty(name = "m", description = "Frame Inside Projection", min = 0.01, max = 10, default = 0.1)
    ftc: FloatProperty(name = "W/m.K", description = "Frame Conductance", min = 0.01, max = 10, default = 0.1)
    fratio: FloatProperty(name = "", description = "Ratio of Frame-Edge Glass Conductance to Center-Of-Glass Conductance", min = 0.1, max = 10, default = 1)
    fsa: FloatProperty(name = "", description = "Frame Solar Absorptance", min = 0.01, max = 1, default = 0.7)
    fva: FloatProperty(name = "", description = "Frame Visible Absorptance", min = 0.01, max = 1, default = 0.7)
    fte: FloatProperty(name = "", description = "Frame Thermal Emissivity", min = 0.01, max = 1, default = 0.7)
    dt: EnumProperty(items = [("0", "None", "None"), ("1", "DividedLite", "Divided lites"), ("2", "Suspended", "Suspended divider")], 
                                        name = "", description = "Type of divider", default = "0")
    dw: FloatProperty(name = "m", description = "Divider Width", min = 0.001, max = 10, default = 0.01)
    dhd: IntProperty(name = "", description = "Number of Horizontal Dividers", min = 0, max = 10, default = 0)
    dvd: IntProperty(name = "", description = "Number of Vertical Dividers", min = 0, max = 10, default = 0)
    dop: FloatProperty(name = "m", description = "Divider Outside Projection", min = 0.0, max = 10, default = 0.01)
    dip: FloatProperty(name = "m", description = "Divider Inside Projection", min = 0.0, max = 10, default = 0.01)
    dtc: FloatProperty(name = "W/m.K", description = "Divider Conductance", min = 0.001, max = 10, default = 0.1)
    dratio: FloatProperty(name = "", description = "Ratio of Divider-Edge Glass Conductance to Center-Of-Glass Conductance", min = 0.1, max = 10, default = 1)
    dsa: FloatProperty(name = "", description = "Divider Solar Absorptance", min = 0.01, max = 1, default = 0.7)
    dva: FloatProperty(name = "", description = "Divider Visible Absorptance", min = 0.01, max = 1, default = 0.7)
    dte: FloatProperty(name = "", description = "Divider Thermal Emissivity", min = 0.01, max = 1, default = 0.7)
    orsa: FloatProperty(name = "", description = "Outside Reveal Solar Absorptance", min = 0.01, max = 1, default = 0.7)
    isd: FloatProperty(name = "m", description = "Inside Sill Depth (m)", min = 0.0, max = 10, default = 0.1)
    issa: FloatProperty(name = "", description = "Inside Sill Solar Absorptance", min = 0.01, max = 1, default = 0.7)
    ird: FloatProperty(name = "m", description = "Inside Reveal Depth (m)", min = 0.0, max = 10, default = 0.1)
    irsa: FloatProperty(name = "", description = "Inside Reveal Solar Absorptance", min = 0.01, max = 1, default = 0.7)
    resist: FloatProperty(name = "", description = "U-value", min = 0.01, max = 10, default = 0.7)
    
    def init(self, context):
        self.inputs.new('So_En_Mat_PV', 'PV')
        self.inputs['PV'].hide = True
        self.inputs.new('So_En_Mat_Ou', 'Outer layer')
        self.inputs['Outer layer'].hide = True
        self.inputs.new('So_En_Sched', 'Schedule')
        self.inputs['Schedule'].hide = True
#        self.inputs.new('SO_EN_Mat_Fr', 'Outer frame layer')
#        self.inputs['Outer frame layer'].hide = True
        
    def draw_buttons(self, context, layout):
        newrow(layout, 'Active:', self, 'active')
        newrow(layout, 'Type:', self, "envi_con_type")
       
        if self.envi_con_type not in ('None', 'Shading'):
            newrow(layout, 'Boundary:', self, "envi_con_con")
            
            if self.envi_con_type != "Shading":
                if self.envi_con_con in ('External', 'Zone') and not self.inputs['PV'].links:
                    newrow(layout, 'Air-flow:', self, "envi_afsurface")
                newrow(layout, 'Specification:', self, "envi_con_makeup")
                
                if self.envi_con_makeup == '0':                    
                    if self.envi_con_type == 'Window':
                        newrow(layout, 'Simple glazing:', self, "envi_simple_glazing")
    
                        if self.envi_simple_glazing:
                            newrow(layout, 'U-Value:', self, "envi_sg_uv")
                            newrow(layout, 'SHGC:', self, "envi_sg_shgc")
                            newrow(layout, 'Vis trans.:', self, "envi_sg_vt")
                    
                    if self.envi_con_type != 'Window' or not self.envi_simple_glazing:
                        row = layout.row()                
                        row.prop(self, 'envi_con_list')
                        
                        con_type = {'Roof': 'Ceiling', 'Floor': 'Internal floor', 'Wall': 'Internal wall'}[self.envi_con_type] if self.envi_con_con in ('Thermal mass', 'Zone') and self.envi_con_type in ('Roof', 'Wall', 'Floor') else self.envi_con_type
        
                        for l, layername in enumerate(envi_cons.propdict[con_type][self.envi_con_list]):    
                            row = layout.row()
                            
                            if layername in envi_mats.wgas_dat:
                                row.label(text = '{} ({})'.format(layername, "14mm"))
                                row.prop(self, "lt{}".format(l))
    
                            elif layername in envi_mats.gas_dat:
                                row.label(text = '{} ({})'.format(layername, "20-50mm"))
                                row.prop(self, "lt{}".format(l))
    
                            elif layername in envi_mats.glass_dat:
                                row.label(text = '{} ({})'.format(layername, "{}mm".format(float(envi_mats.matdat[layername][3])*1000)))
                                row.prop(self, "lt{}".format(l))
    
                            else:
                                row.label(text = '{} ({})'.format(layername, "{}mm".format(envi_mats.matdat[layername][7])))
                                row.prop(self, "lt{}".format(l))
                   
            if self.envi_con_type in ('Window', 'Door'):
                newrow(layout, 'Frame:', self, "fclass")
                if self.fclass == '0' or self.envi_con_type == 'Door':
                    newrow(layout, 'Material:', self, "fmat")
                    newrow(layout, '% frame area:', self, "farea")
                    
                elif self.fclass == '1':
                    newrow(layout, "Width:", self, "fw")
                    newrow(layout, "Outer p:", self, "fop")
                    newrow(layout, "Inner p:", self, "fip")
                    newrow(layout, "Conductivity:", self, "ftc") 
                    newrow(layout, "Cond. ratio:", self, "fratio")
                    newrow(layout, "Solar absorp.:", self, "fsa")
                    newrow(layout, "Visible trans.:", self, "fva")
                    newrow(layout, "Thermal emmis.:", self, "fte")
                    newrow(layout, "Divider type:", self, "dt")
        
                    if self.dt != '0':
                        row = layout.row()
                        row.label('--Divider--')
                        newrow(layout, "    Width:", self, "dw")
                        newrow(layout, "    No. (h):", self, "dhd")
                        newrow(layout, "    No. (v):", self, "dvd")
                        newrow(layout, "    Outer proj.:", self, "dop")
                        newrow(layout, "    Inner proj.:", self, "dip")
                        newrow(layout, "    Conductivity:", self, "dtc")
                        newrow(layout, "    Cond. ratio:", self, "dratio")
                        newrow(layout, "    Sol. abs.:", self, "dsa")
                        newrow(layout, "    Vis. abs.:", self, "dva")
                        newrow(layout, "    Emissivity:", self, "dte")
                    newrow(layout, "Reveal sol. abs.:", self, "orsa")
                    newrow(layout, "Sill depth:", self, "isd")
                    newrow(layout, "Sill sol. abs.:", self, "issa")
                    newrow(layout, "Reveal depth:", self, "ird")
                    newrow(layout, "Inner reveal sol. abs:", self, "irsa")
                else:
                    newrow(layout, '% area:', self, "farea")
            
            elif self.envi_con_type in ('Wall', 'Floor', 'Roof'):
                row = layout.row()
                if self.envi_con_makeup == '0':
                    try:                        
                        row.label(text = 'U-value  = {} W/m2.K'.format(self.uv)) 
                    except: 
                        row.label(text = 'U-value  = N/A') 
                elif self.envi_con_makeup == '1':
                    row.operator('node.envi_uv', text = "UV Calc")
                    try:                        
                        row.label(text = 'U-value  = {} W/m2.K'.format(self.uv)) 
                    except: 
                        row.label(text = 'U-value  = N/A')
                        
        if self.envi_con_makeup == '1' and self.envi_con_type != 'Shading':
            newrow(layout, 'Name:', self, "con_name")
            if self.con_name:
                row = layout.row()
                row.operator('node.con_save', text = "Save")
#            elif self.envi_con_type == 'PV' and self.envi_con_makeup == '0':
#                newrow(layout, "Series in parallel:", self, "ssp")
#                newrow(layout, "Modules in series:", self, "mis")
#                newrow(layout, "Area:", self, "fsa")
#                newrow(layout, "Efficiency:", self, "eff")
        
    def update(self):
        if len(self.inputs) == 3:
            self.valid()
    
    def valid(self):
        if (self.envi_con_makeup == '1' and not self.inputs['Outer layer'].links and self.envi_con_type != 'Shading'):
            nodecolour(self, 1)
        else:
            nodecolour(self, 0)
    
    def ret_uv(self):
        if self.envi_con_makeup == '1':
            resists = []
            lsock = self.inputs['Outer layer']
            
            while lsock.links:
                resists.append(lsock.links[0].from_node.ret_resist())
                lsock = lsock.links[0].from_node.inputs['Layer']   

            self.uv = '{:.3f}'.format(1/(sum(resists) + 0.12 + 0.08))
        return self.uv
    
    def ret_nodes(self):
        nodes = [self]
        lsock = self.inputs['Outer layer']
            
        while lsock.links:
            nodes.append(lsock.links[0].from_node)
            lsock = lsock.links[0].from_node.inputs['Layer']
        
        return nodes
    
    def save_condict(self):
        lks = self.inputs['Outer layer'].links
        lay_names = [lks[0].from_node.lay_name] if lks[0].from_node.layer == '1' else [lks[0].from_node.material]
        
        while lks:
            lks = lks[0].from_node.inputs[0].links
            if lks:
                lay_name = lks[0].from_node.lay_name if lks[0].from_node.layer == '1' else lks[0].from_node.material
                lay_names.append(lay_name)
            
        envi_cons.get_dat('{} - {}'.format(self.envi_con_type, self.envi_con_con))[self.con_name] = lay_names   
        envi_cons.get_dat('{} - {}'.format(self.envi_con_type, self.envi_con_con))['{} (reversed)'.format(self.con_name)] = lay_names[::-1] 
        envi_cons.con_save()
               
    def pv_ep_write(self, sn):
        self['matname'] = get_mat(self, 1).name
        
        params = ('Name', 'Surface Name', 'Photovoltaic Performance Object Type', 
                  'Module Performance Name', 'Heat Transfer Integration Mode', 
                  'Number of Series Strings in Parallel', 'Number of Modules in Series')
                
        paramvs = ['{}-pv'.format(sn), sn, 
                   ('PhotovoltaicPerformance:Simple', 'PhotovoltaicPerformance:EquivalentOne-Diode', 'PhotovoltaicPerformance:Sandia')[int(self.pp)], '{}-pv-performance'.format(sn),
                   self.hti, self.ssp, self.mis]

        ep_text = epentry('Generator:Photovoltaic', params, paramvs)
        
        if self.pp == '0':
            params = ('Name', 'Fraction of Surface Area with Active Solar Cell', 'Conversion Efficiency Input Mode', 'Value for Cell Efficiency if Fixed', 'Efficiency Schedule Name')
            paramvs = ('{}-pv-performance'.format(sn), self.pvsa * 0.01, ('Fixed', 'Scheduled')[len(self.inputs['PV Schedule'].links)], self.eff * 0.01, ('', '{}-pv-performance-schedule'.format(sn))[len(self.inputs['PV Schedule'].links)])
            ep_text += epentry('PhotovoltaicPerformance:Simple', params, paramvs)
            
            if self.inputs['PV Schedule'].links:
                ep_text += self.inputs['PV Schedule'].links[0].from_node.epwrite('{}-pv-performance-schedule'.format(sn), 'Fraction')
            
        elif self.pp == '1':
            params = ('Name', 'Cell type', 'Number of Cells in Series', 'Active Area (m2)', 'Transmittance Absorptance Product',
                      'Semiconductor Bandgap (eV)', 'Shunt Resistance (ohms)', 'Short Circuit Current (A)', 'Open Circuit Voltage (V)',
                      'Reference Temperature (C)', 'Reference Insolation (W/m2)', 'Module Current at Maximum Power (A)', 
                      'Module Voltage at Maximum Power (V)', 'Temperature Coefficient of Short Circuit Current (A/K)',
                      'Temperature Coefficient of Open Circuit Voltage (V/K)', 'Nominal Operating Cell Temperature Test Ambient Temperature (C)',
                      'Nominal Operating Cell Temperature Test Cell Temperature (C)', 'Nominal Operating Cell Temperature Test Insolation (W/m2)',
                      'Module Heat Loss Coefficient (W/m2-K)', 'Total Heat Capacity (J/m2-K)')
            paramvs = ('{}-pv-performance'.format(sn), ('CrystallineSilicon', 'AmorphousSilicon')[int(self.ct)], (self.cis, self.e1ddict[self.e1menu][6])[self.e1menu != 'Custom'], (self.aa, self.e1ddict[self.e1menu][8])[self.e1menu != 'Custom'],
                       self.tap, self.sbg, self.sr, (self.scc, self.e1ddict[self.e1menu][0])[self.e1menu != 'Custom'], (self.ocv, self.e1ddict[self.e1menu][1])[self.e1menu != 'Custom'],
                       self.rt, self.ri, (self.mc, self.e1ddict[self.e1menu][3])[self.e1menu != 'Custom'], (self.mv, self.e1ddict[self.e1menu][2])[self.e1menu != 'Custom'],
                       (self.tcscc, self.e1ddict[self.e1menu][4])[self.e1menu != 'Custom'], (self.tcocv, self.e1ddict[self.e1menu][5])[self.e1menu != 'Custom'],
                       self.atnoct, (self.ctnoct, self.e1ddict[self.e1menu][7] - 273.14)[self.e1menu != 'Custom'], self.inoct, self.hlc, self.thc)
            ep_text += epentry('PhotovoltaicPerformance:EquivalentOne-Diode', params, paramvs)        
        return ep_text
     
    def ep_write(self, mn):
        self['matname'] = get_mat(self, 1).name
        print('1', self['matname'], mn)
        con_type = {'Roof': 'Ceiling', 'Floor': 'Internal floor', 'Wall': 'Internal wall'}[self.envi_con_type] if self.envi_con_con in ('Thermal mass', 'Zone') and self.envi_con_type in ('Roof', 'Wall', 'Floor') else self.envi_con_type
        envi_mats = envi_materials()
        
        if self.envi_con_makeup == '0':
            if self.envi_con_type == 'Window' and self.envi_simple_glazing:
                params = ['Name', 'Outside layer']
                paramvs = [mn, mn + '_sg']
                ep_text = epentry('Construction', params, paramvs)
                params = ('Name', 'U-Factor', 'Solar Heat Gain Coefficient', 'Visible Transmittance')
                paramvs = [self['matname'] + '_sg'] + ['{:.3f}'.format(p) for p in (self.envi_sg_uv, self.envi_sg_shgc, self.envi_sg_vt)]                
                ep_text += epentry("WindowMaterial:SimpleGlazingSystem", params, paramvs)                          
            else:
                self.thicklist = [self.lt0, self.lt1, self.lt2, self.lt3, self.lt4, self.lt5, self.lt6, self.lt7, self.lt8, self.lt9]
                mats = envi_cons.propdict[con_type][self.envi_con_list]
                params = ['Name', 'Outside layer'] + ['Layer {}'.format(i + 1) for i in range(len(mats) - 1)]        
                paramvs = [mn] + ['{}-layer-{}'.format(mn, mi) for mi, m in enumerate(mats)]
                ep_text = epentry('Construction', params, paramvs)
                
                for pm, presetmat in enumerate(mats):  
                    matlist = list(envi_mats.matdat[presetmat])
                    layer_name = '{}-layer-{}'.format(mn, pm)
                    
                    if envi_mats.namedict.get(presetmat) == None:
                        envi_mats.namedict[presetmat] = 0
                        envi_mats.thickdict[presetmat] = [self.thicklist[pm]/1000]
                    else:
                        envi_mats.namedict[presetmat] = envi_mats.namedict[presetmat] + 1
                        envi_mats.thickdict[presetmat].append(self.thicklist[pm]/1000)
                    
                    if self.envi_con_type in ('Wall', 'Floor', 'Roof', 'Ceiling', 'Door') and presetmat not in envi_mats.gas_dat:
                        self.resist += self.thicklist[pm]/1000/float(matlist[1])
                        params = ('Name', 'Roughness', 'Thickness (m)', 'Conductivity (W/m-K)', 'Density (kg/m3)', 'Specific Heat Capacity (J/kg-K)', 'Thermal Absorptance', 'Solar Absorptance', 'Visible Absorptance')                    
                        paramvs = ['{}-layer-{}'.format(mn, pm), matlist[0], str(self.thicklist[pm]/1000)] + matlist[1:8]                    
                        ep_text += epentry("Material", params, paramvs)
    
                        if presetmat in envi_mats.pcmd_datd:
                            stringmat = envi_mats.pcmd_datd[presetmat]
                            params = ('Name', 'Temperature Coefficient for Thermal Conductivity (W/m-K2)')
                            paramvs = ('{}-layer-{}'.format(self['matname'], pm), stringmat[0])
                            
                            for i, te in enumerate(stringmat[1].split()):
                                params += ('Temperature {} (C)'.format(i), 'Enthalpy {} (J/kg)'.format(i))
                                paramvs +=(te.split(':')[0], te.split(':')[1])
                                
                            ep_text += epentry("MaterialProperty:PhaseChange", params, paramvs)
                            pcmparams = ('Name', 'Algorithm', 'Construction Name')
                            pcmparamsv = ('{} CondFD override'.format(mn), 'ConductionFiniteDifference', mn)
                            ep_text += epentry('SurfaceProperty:HeatTransferAlgorithm:Construction', pcmparams, pcmparamsv)
    
                    elif presetmat in envi_mats.gas_dat:
                        params = ('Name', 'Resistance')
                        paramvs = ('{}-layer-{}'.format(mn, pm), matlist[2])
                        ep_text += epentry("Material:AirGap", params, paramvs)
                    
                    elif self.envi_con_type =='Window':
                        if envi_mats.matdat[presetmat][0] == 'Glazing':
                            params = ('Name', 'Optical Data Type', 'Window Glass Spectral Data Set Name', 'Thickness (m)', 'Solar Transmittance at Normal Incidence', 'Front Side Solar Reflectance at Normal Incidence',
                          'Back Side Solar Reflectance at Normal Incidence', 'Visible Transmittance at Normal Incidence', 'Front Side Visible Reflectance at Normal Incidence', 'Back Side Visible Reflectance at Normal Incidence',
                          'Infrared Transmittance at Normal Incidence', 'Front Side Infrared Hemispherical Emissivity', 'Back Side Infrared Hemispherical Emissivity', 'Conductivity (W/m-K)',
                          'Dirt Correction Factor for Solar and Visible Transmittance', 'Solar Diffusing')
                            paramvs = ['{}-layer-{}'.format(mn, pm)] + matlist[1:3] + [self.thicklist[pm]] + ['{:.3f}'.format(float(sm)) for sm in matlist[4:-1]] + [1, ('No', 'Yes')[matlist[-1]]]
                            ep_text += epentry("WindowMaterial:{}".format(matlist[0]), params, paramvs)
                    
                        elif envi_mats.matdat[presetmat][0] == 'Gas':
                            params = ('Name', 'Gas Type', 'Thickness')
                            paramvs = [layer_name] + [matlist[1]] + [self.thicklist[pm]]
                            ep_text += epentry("WindowMaterial:Gas", params, paramvs)
                    
        elif self.envi_con_makeup == '1':            
            in_sock = self.inputs['Outer layer']# if self.envi_con_type == "Window" else self.inputs[0]
            n = 0
            params = ['Name']
            paramvs = [mn]
            ep_text = ''
            self.resist = 0
            get_mat(self, 1).vi_params.envi_shading = 0

            while in_sock.links:
                node = in_sock.links[0].from_node

                if node.bl_idname not in ('envi_sl_node', 'envi_bl_node', 'envi_screen_node', 'envi_sgl_node'):                    
                    paramvs.append('{}-layer-{}'.format(mn, n)) 
                    params.append(('Outside layer', 'Layer {}'.format(n))[n > 0])
                    ep_text += node.ep_write(n, mn)  
                    self.resist += node.resist
                else:
                    get_mat(self, 1).vi_params.envi_shading = 1
                    
                in_sock = node.inputs['Layer']
                n += 1
                
            ep_text += epentry('Construction', params, paramvs)
            
            if get_mat(self, 1).vi_params.envi_shading:
                in_sock = self.inputs['Outer layer']
                n = 0
                params = ['Name'] 
                paramvs = ['{}-shading'.format(mn)]
                
                while in_sock.links:
                    node = in_sock.links[0].from_node
                    
                    if node.outputs['Layer'].links[0].to_node.bl_idname != 'envi_sgl_node':
                        paramvs.append('{}-layer-{}'.format(mn, n)) 
                        params.append(('Outside layer', 'Layer {}'.format(n))[n > 0])
                    
                    in_sock = node.inputs['Layer']

                    if node.bl_idname in ('envi_sl_node', 'envi_bl_node', 'envi_screen_node', 'envi_sgl_node'):
                        ep_text += node.ep_write(n, mn)
                    
                    n += 1
                ep_text += epentry('Construction', params, paramvs)
                
        if self.envi_con_type in ('Window', 'Door'): 
            if self.fclass == '0' or self.envi_con_type == 'Door':
                params = ('Name', 'Roughness', 'Thickness (m)', 'Conductivity (W/m-K)', 'Density (kg/m3)', 'Specific Heat (J/kg-K)', 'Thermal Absorptance', 'Solar Absorptance', 'Visible Absorptance', 'Name', 'Outside Layer')
                paramvs = ('{}-frame-layer{}'.format(mn, 0), 'Smooth', '0.12', '0.1', '1400.00', '1000', '0.9', '0.6', '0.6', '{}-frame'.format(mn), '{}-frame-layer{}'.format(mn, 0))
                ep_text += epentry('Material', params[:-2], paramvs[:-2])
                ep_text += epentry('Construction', params[-2:], paramvs[-2:])
            
            elif self.fclass == '1':
                fparams = ('Frame/Divider Name', 'Frame Width', 'Frame Outside Projection', 'Frame Inside Projection', 'Frame Conductance', 
                           'Ratio of Frame-Edge Glass Conductance to Center-Of-Glass Conductance', 'Frame Solar Absorptance', 'Frame Visible Absorptance', 
                           'Frame Thermal Emissivity', 'Divider Type', 'Divider Width', 'Number of Horizontal Dividers', 'Number of Vertical Dividers',
                           'Divider Outside Projection', 'Divider Inside Projection', 'Divider Conductance', 'Ratio of Divider-Edge Glass Conductance to Center-Of-Glass Conductance',
                           'Divider Solar Absorptance', 'Divider Visible Absorptance', 'Divider Thermal Emissivity', 'Outside Reveal Solar Absorptance',
                           'Inside Sill Depth (m)', 'Inside Sill Solar Absorptance', 'Inside Reveal Depth (m)', 'Inside Reveal Solar Absorptance')
                fparamvs = ['{}-fad'.format(self['matname'])] +  ['{:.3f}'.format(p) for p in (self.fw, self.fop, self.fip, self.ftc, self.fratio, self.fsa, self.fva, self.fte)] +\
                            [('', 'DividedLite', 'Suspended')[int(self.dt)]] + ['{:.3f}'.format(p) for p in (self.dw, self.dhd, self.dvd, self.dop, self.dip, self.dtc, self.dratio, self.dsa, self.dva, self.dte, self.orsa, self.isd, 
                            self.issa, self.ird, self.irsa)]
                ep_text += epentry('WindowProperty:FrameAndDivider', fparams, fparamvs)
                
            elif self.fclass == '2':
                ep_text += self.layer_write(self.inputs['Outer frame layer'], mn)
        
        return ep_text
    
    def layer_write(self, in_sock, matname):
        ep_text = ''
        n = 0
        params = ['Name']
        paramvs = ['{}-frame'.format(matname)]
        
        while in_sock.links:
            node = in_sock.links[0].from_node
            paramvs.append('{}-frame-layer-{}'.format(matname, n)) 
            params.append(('Outside layer', 'Layer {}'.format(n))[n > 0])
            ep_text += node.ep_write(n)                    
            in_sock = node.inputs['Layer']
            n += 1
            
        ep_text += epentry('Construction', params, paramvs)
        return ep_text
    
class No_En_Mat_Op(Node, EnViMatNodes):
    '''Node defining the EnVi opaque material layer'''
    bl_idname = 'No_En_Mat_Op'
    bl_label = 'EnVi opaque layer'
    
    def lay_update(self, context):
        if self.layer == '1' and self.lay_name == '':
            nodecolour(self, 1) 
        elif self.layer == '0' and not self.material:
            nodecolour(self, 1) 
        else:
            nodecolour(self, 0)
             
    
    lay_name: StringProperty(name = '', description = 'Custom layer name', update = lay_update)
    layer: EnumProperty(items = [("0", "Database", "Select from database"), 
                                        ("1", "Custom", "Define custom material properties")], 
                                        name = "", description = "Class of layer", default = "0", update = lay_update)

    materialtype: EnumProperty(items = envi_layertype, name = "", description = "Layer material type", update = lay_update)
    material: EnumProperty(items = envi_layer, name = "", description = "Layer material", update = lay_update)
    thi: FloatProperty(name = "mm", description = "Thickness (mm)", min = 0.1, max = 10000, default = 100, options={'ANIMATABLE'})
    tc: FloatProperty(name = "W/m.K", description = "Thickness (mm)", min = 0.001, max = 10, default = 0.5)
    rough: EnumProperty(items = [("VeryRough", "VeryRough", "Roughness"), 
                                  ("Rough", "Rough", "Roughness"), 
                                  ("MediumRough", "MediumRough", "Roughness"),
                                  ("MediumSmooth", "MediumSmooth", "Roughness"), 
                                  ("Smooth", "Smooth", "Roughness"), 
                                  ("VerySmooth", "VerySmooth", "Roughness")],
                                  name = "", 
                                  description = "Specify the material roughness for convection calculations", 
                                  default = "Rough")
    
    rho: FloatProperty(name = "kg/m^3", description = "Density", min = 0.001, max = 10000, default = 800)
    shc: FloatProperty(name = "J/kg", description = "Thickness (mm)", min = 0.01, max = 10000, default = 800)
    tab: FloatProperty(name = "", description = "Thermal absorptance", min = 0, max = 1, precision = 2, default = 0.7)
    sab: FloatProperty(name = "", description = "Solar absorptance", min = 0, max = 1, precision = 2, default = 0.7)
    vab: FloatProperty(name = "", description = "Visual absorptance", min = 0, max = 1, precision = 2, default = 0.7)
    pcm: BoolProperty(name = "", description = "Phase Change Material", default = 0)
    tctc: FloatProperty(name = "", description = "Temp. coeff. for thermal conductivity (W/m-K2)", min = 0, max = 50, default = 0)
    tempemps: StringProperty(name = "", description = "Temperature/empalthy pairs (e.g. T1:E1 T2:E2)", default = "")
    resist: FloatProperty(name = "", description = "", min = 0, default = 0) 
    envi_con_type: StringProperty(name = "", description = "Name")
    
    def init(self, context):
        self.outputs.new('So_En_Mat_Op', 'Layer')
        self.inputs.new('So_En_Mat_Op', 'Layer')
        
    def draw_buttons(self, context, layout):
        newrow(layout, "Type:", self, "materialtype")
        newrow(layout, "Class:", self, "layer")
        if self.layer == '0':
            newrow(layout, "Material:", self, "material")
            newrow(layout, "Thickness:", self, "thi") 
        else:
            newrow(layout, "Name:", self, "lay_name")
            newrow(layout, "Conductivity:", self, "tc")
            newrow(layout, "Thickness:", self, "thi")
            newrow(layout, "Roughness:", self, "rough")
            newrow(layout, "Density:", self, "rho")
            newrow(layout, "SHC:", self, "shc")            
            newrow(layout, "Therm absorb:", self, "tab")
            newrow(layout, "Solar absorb:", self, "sab")
            newrow(layout, "Vis absorb:", self, "vab")
            
            if self.materialtype == '8':
                newrow(layout, "Temps:Emps", self, "tempemps")
                newrow(layout, "TCTC:", self, "tctc")
            
            row = layout.row()
            row.operator('node.lay_save', text = "Layer Save")
    
    def ret_resist(self):
        if self.layer == '0':
            matlist = list(envi_mats.matdat[self.material])
            
            if self.materialtype != '6': 
                self.resist = self.thi * 0.001/float(matlist[1])
            else:
                self.resist = float(matlist[2])            
        else:
            self.resist = self.thi * 0.001/self.tc
        
        return self.resist
    
    def save_laydict(self):
        '''Roughness, Conductivity {W/m-K}, Density {kg/m3}, Specific Heat {J/kg-K}, Thermal Absorbtance, 
            Solar Absorbtance, Visible Absorbtance, Default thickness'''
        envi_mats.get_dat(self.materialtype)[self.lay_name] = [self.rough, '{:.4f}'.format(self.tc), 
                                                             '{:.2f}'.format(self.rho), '{:.2f}'.format(self.shc), 
                                                             '{:.2f}'.format(self.tab), '{:.2f}'.format(self.sab), 
                                                             '{:.2f}'.format(self.vab), self.thi]       
        if self.materialtype == '8':
            envi_mats.get_dat('9')[self.lay_name] = [self.tctc, self.tempemps]
            
        envi_mats.lay_save()
        
    def update(self):
        socklink2(self.outputs['Layer'], self.id_data)
        if self.outputs['Layer'].links:
            self.envi_con_type = self.outputs['Layer'].links[0].to_node.envi_con_type if self.outputs['Layer'].links[0].to_socket.bl_idname != 'envi_f_sock' else 'Frame'
        self.valid()
    
    def valid(self):
        if not self.outputs["Layer"].links:
            nodecolour(self, 1)
        else:
            nodecolour(self, 0)
    
    def pcm_write(self):
        params = ['Name', 'Temperature Coefficient for Thermal Conductivity (W/m-K2)']
        
        if self.layer == '0':
            paramvs = [self['layer_name'], envi_mats.pcmd_dat[self.material][0]]
            mtempemps = envi_mats.pcmd_dat[self.material][1]
        else:
            paramvs = [self['layer_name'], self.tctc]
            mtempemps = self.tempemps
            
        for i, te in enumerate(mtempemps.split()):
            params += ('Temperature {} (C)'.format(i), 'Enthalpy {} (J/kg)'.format(i))
            paramvs +=(te.split(':')[0], te.split(':')[1])

        pcmparams = ('Name', 'Algorithm', 'Construction Name')
        pcmparamsv = ('{} CondFD override'.format(self['matname']), 'ConductionFiniteDifference', self['matname'])
    
        return epentry("MaterialProperty:PhaseChange", params, paramvs) + epentry('SurfaceProperty:HeatTransferAlgorithm:Construction', pcmparams, pcmparamsv)
        
    def ep_write(self, ln, mn):
        for material in bpy.data.materials:
            if self.id_data == material.vi_params.envi_nodes:
                break
        self['matname'] = get_mat(self, 1).name
        self['layer_name'] = '{}-layer-{}'.format(mn, ln) if self.envi_con_type != 'Frame' else '{}-frame-layer-{}'.format(mn, ln)
        
        if self.materialtype != '6':
            params = ('Name', 'Roughness', 'Thickness (m)', 'Conductivity (W/m-K)', 'Density (kg/m3)', 'Specific Heat Capacity (J/kg-K)', 'Thermal Absorptance', 'Solar Absorptance', 'Visible Absorptance')
            header = 'Material'
        else:
            params = ('Name', 'Resistance')
            header = 'Material:AirGap'

        if self.layer == '0':
            matlist = list(envi_mats.matdat[self.material])
            
            if self.materialtype != '6':
                paramvs = [self['layer_name'], matlist[0], '{:.3f}'.format(self.thi * 0.001)] + matlist[1:8]  
            else:
                paramvs = [self['layer_name'], matlist[2]]
            
        else:
            paramvs = ['{}-layer-{}'.format(mn, ln), self.rough, '{:.3f}'.format(self.thi * 0.001), '{:.3f}'.format(self.tc), '{:.3f}'.format(self.rho), '{:.3f}'.format(self.shc), '{:.3f}'.format(self.tab), 
                       '{:.3f}'.format(self.sab), '{:.3f}'.format(self.vab)]
            
        ep_text = epentry(header, params, paramvs)
        
        if self.materialtype == '8':
            ep_text += self.pcm_write()
        
        return ep_text
            
class No_En_Mat_Tr(Node, EnViMatNodes):
    '''Node defining the EnVi transparent material layer'''
    bl_idname = 'No_En_Mat_Tr'
    bl_label = 'EnVi transparent layer'
    
    lay_name: StringProperty(name = '', description = 'Custom layer name')
    layer: EnumProperty(items = [("0", "Database", "Select from database"), 
                                        ("1", "Custom", "Define custom material properties")], 
                                        name = "", description = "Composition of the layer", default = "0")
    materialtype: EnumProperty(items = envi_layertype, name = "", description = "Layer material type")
    material: EnumProperty(items = envi_layer, name = "", description = "Layer material")
    thi: FloatProperty(name = "mm", description = "Thickness (mm)", min = 0.1, max = 1000, default = 6)
    tc: FloatProperty(name = "W/m.K", description = "Thermal Conductivity (W/m.K)", precision = 3, min = 0.1, max = 10, default = 0.8)
    stn: FloatProperty(name = "", description = "Solar normal transmittance", precision = 3, min = 0, max = 1, default = 0.7)
    fsn: FloatProperty(name = "", description = "Solar front normal reflectance", precision = 3, min = 0, max = 1, default = 0.07)
    bsn: FloatProperty(name = "", description = "Solar back normal reflectance", precision = 3, min = 0, max = 1, default = 0.07)
    vtn: FloatProperty(name = "", description = "Visible Transmittance at Normal Incidence", precision = 3, min = 0, max = 1, default = 0.89)
    fvrn: FloatProperty(name = "", description = "Front Side Visible Reflectance at Normal Incidence", precision = 3, min = 0, max = 1, default = 0.07)
    bvrn: FloatProperty(name = "", description = "Back Side Visible Reflectance at Normal Incidence", precision = 3, min = 0, max = 1, default = 0.07)
    itn: FloatProperty(name = "", description = "Infrared Transmittance at Normal Incidence", precision = 3, min = 0, max = 1, default = 0.0)
    fie: FloatProperty(name = "", description = "Front Side Infrared Hemispherical Emissivity", precision = 3, min = 0, max = 1, default = 0.84)
    bie: FloatProperty(name = "", description = "Back Side Infrared Hemispherical Emissivity", precision = 3, min = 0, max = 1, default = 0.84)
    diff: BoolProperty(name = "", description = "Diffusing", default = 0)
    envi_con_type: StringProperty(name = "", description = "Name")
    resist: FloatProperty(name = "", description = "", min = 0, default = 0)
    
    def init(self, context):
        self.outputs.new('So_En_Mat_Tr', 'Layer')
        self.inputs.new('So_En_Mat_Gas', 'Layer')
        
    def draw_buttons(self, context, layout):        
        newrow(layout, "Class:", self, "layer")
        if self.layer == '0':
            newrow(layout, "Material:", self, "material")
            newrow(layout, "Thickness:", self, "thi") 
        else:
            newrow(layout, "Name:", self, "lay_name")
            newrow(layout, "Conductivity:", self, "tc")
            newrow(layout, "Thickness:", self, "thi")
            newrow(layout, "STN:", self, "stn")
            newrow(layout, "FSN:", self, "fsn")
            newrow(layout, "BSN:", self, "bsn")
            newrow(layout, "VTN:", self, "vtn")
            newrow(layout, "FVRN:", self, "fvrn")
            newrow(layout, "BVRN:", self, "bvrn")
            newrow(layout, "ITN:", self, "itn")
            newrow(layout, "FIE:", self, "fie")
            newrow(layout, "BIE:", self, "bie")
            newrow(layout, "Diffuse:", self, "diff")
            
            if self.lay_name:
                row = layout.row()
                row.operator('node.lay_save', text = "Layer Save")
            
    def update(self):
        socklink2(self.outputs['Layer'], self.id_data)

        if self.outputs['Layer'].links:
            self.envi_con_type = self.outputs['Layer'].links[0].to_node.envi_con_type

        self.valid()
    
    def valid(self):
        if not self.outputs["Layer"].links:
            nodecolour(self, 1)
        else:
            nodecolour(self, 0)
    
    def save_laydict(self):
        '''Name', 'Optical Data Type', 'Window Glass Spectral Data Set Name', 'Thickness (m)', 'Solar Transmittance at Normal Incidence', 'Front Side Solar Reflectance at Normal Incidence',
                  'Back Side Solar Reflectance at Normal Incidence', 'Visible Transmittance at Normal Incidence', 'Front Side Visible Reflectance at Normal Incidence', 'Back Side Visible Reflectance at Normal Incidence',
                  'Infrared Transmittance at Normal Incidence', 'Front Side Infrared Hemispherical Emissivity', 'Back Side Infrared Hemispherical Emissivity', 'Conductivity (W/m-K)',
                  'Dirt Correction Factor for Solar and Visible Transmittance', 'Solar Diffusing'''

        envi_mats.get_dat('Glass')[self.lay_name] = ['Glazing', 'SpectralAverage', '', '{:.4f}'.format(self.thi * 0.001), '{:.4f}'.format(self.stn), '{:.4f}'.format(self.fsn), '{:.4f}'.format(self.bsn), '{:.4f}'.format(self.vtn), 
                         '{:.4f}'.format(self.fvrn), '{:.4f}'.format(self.bvrn), '{:.4f}'.format(self.itn), '{:.4f}'.format(self.fie), '{:.4f}'.format(self.bie), '{:.4f}'.format(self.diff)]       
        envi_mats.lay_save()

    def ret_resist(self):
        if self.layer == '0':
            matlist = list(envi_mats.matdat[self.material])            
            self.resist = float(matlist[13])            
        else:
            self.resist = self.thi * 0.001/self.tc
            
        return self.resist
        
    def ep_write(self, ln, mn):
        for material in bpy.data.materials:
            if self.id_data == material.vi_params.envi_nodes:
                break

        layer_name = '{}-layer-{}'.format(mn, ln)
        params = ('Name', 'Optical Data Type', 'Window Glass Spectral Data Set Name', 'Thickness (m)', 'Solar Transmittance at Normal Incidence', 'Front Side Solar Reflectance at Normal Incidence',
                  'Back Side Solar Reflectance at Normal Incidence', 'Visible Transmittance at Normal Incidence', 'Front Side Visible Reflectance at Normal Incidence', 'Back Side Visible Reflectance at Normal Incidence',
                  'Infrared Transmittance at Normal Incidence', 'Front Side Infrared Hemispherical Emissivity', 'Back Side Infrared Hemispherical Emissivity', 'Conductivity (W/m-K)',
                  'Dirt Correction Factor for Solar and Visible Transmittance', 'Solar Diffusing')

        if self.layer == '0':
            matlist = list(envi_mats.matdat[self.material])
            paramvs = [layer_name] + matlist[1:3] + [self.thi] + ['{:.3f}'.format(float(sm)) for sm in matlist[4:-1]] + [1, ('No', 'Yes')[matlist[-1]]]
            
        else:
            paramvs = ['{}-layer-{}'.format(mn, ln), 'SpectralAverage', '', self.thi/1000, '{:.3f}'.format(self.stn), '{:.3f}'.format(self.fsn), '{:.3f}'.format(self.bsn), 
                       '{:.3f}'.format(self.vtn), '{:.3f}'.format(self.fvrn), '{:.3f}'.format(self.bvrn), '{:.3f}'.format(self.itn),
                       '{:.3f}'.format(self.fie), '{:.3f}'.format(self.bie), '{:.3f}'.format(self.tc), 1, ('No', 'Yes')[self.diff]]

        return epentry("WindowMaterial:Glazing", params, paramvs)
    
    
class No_En_Mat_Gas(Node, EnViMatNodes):
    '''Node defining the EnVi transparent gas layer'''
    bl_idname = 'No_En_Mat_Gas'
    bl_label = 'EnVi gas layer'
    
    layer: EnumProperty(items = [("0", "Database", "Select from database"), 
                                        ("1", "Custom", "Define custom material properties")], 
                                        name = "", description = "Composition of the layer", default = "0")
    materialtype: EnumProperty(items = envi_layertype, name = "", description = "Layer material type")
    material: EnumProperty(items = envi_layer, name = "", description = "Layer material")
    thi: FloatProperty(name = "mm", description = "Thickness (mm)", min = 0.1, max = 1000, default = 14)
    ccA: FloatProperty(name = "W/m.K", description = "Conductivity coefficient A", min = 0.1, max = 10, default = 0.003, precision = 5)
    ccB: FloatProperty(name = "W/m.K^2", description = "Conductivity coefficient B", min = 0.0, max = 10, default = 0.00008, precision = 5)
    ccC: FloatProperty(name = "W/m.K^3", description = "Conductivity coefficient C", min = 0.0, max = 10, default = 0, precision = 5)
    vcA: FloatProperty(name = "kg/m.s", description = "Viscosity coefficient A", min = 0.1, max = 10000, default = 800)
    vcB: FloatProperty(name = "kg/m.s.K", description = "Viscosity coefficient B", min = 0, max = 1, default = 0.7)
    vcC: FloatProperty(name = "kg/m.s.K^2", description = "Viscosity coefficient C", min = 0, max = 1, default = 0.7)
    shcA: FloatProperty(name = "J/kg.K", description = "Specific heat coefficient A", min = 0, max = 1, default = 0.7)
    shcB: FloatProperty(name = "J/kg.K^2", description = "Specific heat coefficient A", min = 0, max = 1, default = 0.7)
    shcC: FloatProperty(name = "J/kg.K^3", description = "Specific heat coefficient A", min = 0, max = 1, default = 0.7)
    mw: FloatProperty(name = "kg/kmol", description = "Molecular weight", min = 20, max = 100, default = 20)
    shr: FloatProperty(name = "", description = "Specific heat ratio", min = 1, max = 10, default = 2)
    resist: FloatProperty(name = "", description = "", min = 0, default = 0)
    envi_con_type: StringProperty(name = "", description = "Name")
    
    def init(self, context):
        self.outputs.new('So_En_Mat_Gas', 'Layer')
        self.inputs.new('So_En_Mat_Tr', 'Layer')
        
    def draw_buttons(self, context, layout):
        if self.outputs['Layer'].links:
            newrow(layout, "Class:", self, "layer")
            if self.layer == '0':
                newrow(layout, "Material:", self, "material")
                newrow(layout, "Thickness:", self, "thi") 
            else:
                newrow(layout, "Thickness:", self, "thi")
                newrow(layout, "Coeff A:", self, "ccA")
                newrow(layout, "Coeff B:", self, "ccB")
                newrow(layout, "Coeff C:", self, "ccC")
                newrow(layout, "Viscosity A:", self, "vcA")
                newrow(layout, "Viscosity A:", self, "vcB")
                newrow(layout, "Viscosity A:", self, "vcC")
                newrow(layout, "SHC A:", self, "shcA")
                newrow(layout, "SHC A:", self, "shcB")
                newrow(layout, "SHC A:", self, "shcC")
                newrow(layout, "Mol Weight:", self, "mw")
                newrow(layout, "SHR:", self, "shr")

    def update(self):
        socklink2(self.outputs['Layer'], self.id_data)
        if self.outputs['Layer'].links:
            self.envi_con_type = self.outputs['Layer'].links[0].to_node.envi_con_type
        self.valid()
    
    def valid(self):
        if not self.outputs["Layer"].links or not self.inputs["Layer"].links:
            nodecolour(self, 1)
        else:
            nodecolour(self, 0)

    def ret_resist(self):
        if self.layer == '0':
            matlist = list(envi_mats.matdat[self.material])
            self.resist = self.thi * 0.001/float(matlist[4])          
        else:
            self.tc = self.ccA + self.ccB * 293.14
            self.resist = self.thi * 0.001/self.tc
            
        return self.resist
    
    def ep_write(self, ln, mn):
        for material in bpy.data.materials:
            if self.id_data == material.vi_params.envi_nodes:
                break
        if self.layer == '0':
            params = ('Name', 'Gas Type', 'Thickness')
            paramvs = ['{}-layer-{}'.format(mn, ln), self.material, self.thi]
            
        else:
            params = ('gap name', 'type', 'thickness', 'Conductivity Coefficient A', 'Conductivity Coefficient B', 'Conductivity Coefficient C', 
                      'Conductivity Viscosity A', 'Conductivity Viscosity B', 'Conductivity Viscosity C', 'Specific Heat Coefficient A',
                      'Specific Heat Coefficient B', 'Specific Heat Coefficient C', 'Molecular Weight', 'Specific Heat Ratio')
            paramvs = ['{}-layer-{}'.format(mn, ln), 'Custom', '{:.3f}'.format(self.thi), '{:.3f}'.format(self.ccA), '{:.3f}'.format(self.ccB), 
                       '{:.3f}'.format(self.ccC), '{:.3f}'.format(self.vcA), '{:.3f}'.format(self.vcB), '{:.3f}'.format(self.vcC), '{:.3f}'.format(self.shcA),
                       '{:.3f}'.format(self.shcB), '{:.3f}'.format(self.shcC), '{:.3f}'.format(self.mw), '{:.3f}'.format(self.shr)]
   
        return epentry("WindowMaterial:Gas", params, paramvs)

class No_En_Mat_Sh(Node, EnViMatNodes):
    '''Node defining an EnVi window shader'''
    bl_idname = 'No_En_Mat_Sh'
    bl_label = 'EnVi shade'
        
    st: FloatProperty(name = "", description = "Solar transmittance", min = 0.0, max = 1, default = 0.05)
    sr: FloatProperty(name = "", description = "Solar reflectance", min = 0.0, max = 1, default = 0.3)
    vt: FloatProperty(name = "", description = "Visible transmittance", min = 0.0, max = 1, default = 0.05)
    vr: FloatProperty(name = "", description = "Visible reflectance", min = 0.0, max = 1, default = 0.3)
    ihe: FloatProperty(name = "", description = "Infrared Hemispherical Emissivity", min = 0.0, max = 1, default = 0.9)
    it: FloatProperty(name = "", description = "Infrared Transmittance", min = 0.0, max = 1, default = 0.0)
    thi: FloatProperty(name = "mm", description = "Thickness", min = 0.1, max = 1000, default = 5)
    tc: FloatProperty(name = "W/m.K", description = "Conductivity", min = 0.0001, max = 10, default = 0.1)
    sgd: FloatProperty(name = "mm", description = "Shade to glass distance", min = 0.1, max = 1000, default = 50)
    tom: FloatProperty(name = "", description = "Top opening multiplier", min = 0.0, max = 1, default = 0.5)
    bom: FloatProperty(name = "", description = "Bottom opening multiplier", min = 0.0, max = 1, default = 0.5)
    lom: FloatProperty(name = "", description = "Left-side opening multiplier", min = 0.0, max = 1, default = 0.5)
    rom: FloatProperty(name = "", description = "Right-side opening multiplier", min = 0.0, max = 1, default = 0.5)
    afp: FloatProperty(name = "", description = "Air flow permeability", min = 0.0, max = 1, default = 0.)
    envi_con_type: StringProperty(name = "", description = "Name")
    
    def init(self, context):
        self.outputs.new('envi_sl_sock', 'Layer')
        self.inputs.new('envi_sc_sock', 'Control')        
        self.inputs.new('envi_sl_sock', 'Layer')
                
    def draw_buttons(self, context, layout):
        if self.outputs['Layer'].links:
            newrow(layout, "Solar trans.:", self, "st")
            newrow(layout, "Solar reflec.:", self, "sr")
            newrow(layout, "Vis. trans.:", self, "vt")
            newrow(layout, "Vis. reflec.:", self, "vr") 
            newrow(layout, "IHE:", self, "ihe")
            newrow(layout, "Infra. trans.:", self, "it")
            newrow(layout, "Thickness:", self, "thi")
            newrow(layout, "Conductivity:", self, "tc")
            newrow(layout, "Glass distance:", self, "sgd")
            newrow(layout, "Top mult.:", self, "tom")
            newrow(layout, "Bottom mult.:", self, "bom")
            newrow(layout, "Left milt.:", self, "lom")
            newrow(layout, "Right mult.:", self, "rom")
            newrow(layout, "Air perm.:", self, "afp")
    
    def valid(self):
        if not self.outputs["Layer"].links or not self.inputs["Layer"].links:
            nodecolour(self, 1)
        else:
            nodecolour(self, 0)
            
    def update(self):
        socklink2(self.outputs['Layer'], self.id_data)
        self.valid()
        
    def ep_write(self, ln, mn):
#        for material in bpy.data.materials:
#            if self.id_data == material.envi_nodes:
#                break
        params = ('Name', 'Solar transmittance', 'Solar Reflectance', 'Visible reflectance', 'Infrared Hemispherical Emissivity', 'Infrared Transmittance', 'Thickness {m}',
                  'Conductivity {W/m-K}', 'Shade to glass distance {m}', 'Top opening multiplier', 'Top opening multiplier', 'Bottom opening multiplier', 'Left-side opening multiplier',
                  'Right-side opening multiplier', 'Air flow permeability')
        paramvs = ['{}-layer-{}'.format(mn, ln)] + ['{:.3f}'.format(p) for p in (self.st, self.sr, self.vt, self.vr, self.ihe, self.it, 0.001 * self.thi, self.tc, 0.001 * self.sgd,
                   self.tom, self.bom, self.lom, self.rom, self.afp)]
  
        return epentry('WindowMaterial:Shade', params, paramvs) + self.inputs['Control'].links[0].from_node.ep_write(ln)
        
class No_En_Mat_Sc(Node, EnViMatNodes):
    '''Node defining an EnVi external screen'''
    bl_idname = 'No_En_Mat_Sc'
    bl_label = 'EnVi screen'
    
    rb: EnumProperty(items = [("DoNotModel", "DoNotModel", "Do not model reflected beam component"), 
                               ("ModelAsDirectBeam", "ModelAsDirectBeam", "Model reflectred beam as beam"),
                               ("ModelAsDiffuse", "ModelAsDiffuse", "Model reflected beam as diffuse")], 
                                name = "", description = "Composition of the layer", default = "ModelAsDiffuse")
    ta: EnumProperty(items = [("0", "0", "Angle of Resolution for Screen Transmittance Output Map"), 
                               ("1", "1", "Angle of Resolution for Screen Transmittance Output Map"),
                               ("2", "2", "Angle of Resolution for Screen Transmittance Output Map"),
                               ("3", "3", "Angle of Resolution for Screen Transmittance Output Map"),
                               ("5", "5", "Angle of Resolution for Screen Transmittance Output Map")], 
                                name = "", description = "Angle of Resolution for Screen Transmittance Output Map", default = "0")

    dsr: FloatProperty(name = "", description = "Diffuse solar reflectance", min = 0.0, max = 0.99, default = 0.5)
    vr: FloatProperty(name = "", description = "Visible reflectance", min = 0.0, max = 1, default = 0.6)
    the: FloatProperty(name = "", description = "Thermal Hemispherical Emissivity", min = 0.0, max = 1, default = 0.9)
    tc: FloatProperty(name = "W/m.K", description = "Conductivity", min = 0.0001, max = 10, default = 0.1)
    sme: FloatProperty(name = "mm", description = "Screen Material Spacing", min = 1, max = 1000, default = 50)
    smd: FloatProperty(name = "mm", description = "Screen Material Diameter", min = 1, max = 1000, default = 25)
    sgd: FloatProperty(name = "mm", description = "Screen to glass distance", min = 1, max = 1000, default = 25)
    tom: FloatProperty(name = "", description = "Top opening multiplier", min = 0.0, max = 1, default = 0.5)
    bom: FloatProperty(name = "", description = "Bottom opening multiplier", min = 0.0, max = 1, default = 0.5)
    lom: FloatProperty(name = "", description = "Left-side opening multiplier", min = 0.0, max = 1, default = 0.5)
    rom: FloatProperty(name = "", description = "Right-side opening multiplier", min = 0.0, max = 1, default = 0.5)
    
    def init(self, context):
        self.outputs.new('envi_screen_sock', 'Outer Layer')
        self.inputs.new('envi_sc_sock', 'Control')
        self.inputs.new('envi_tl_sock', 'Layer')
        
    def draw_buttons(self, context, layout):
        if self.outputs['Outer Layer'].links:
            newrow(layout, "Reflected beam:", self, "rb")
            newrow(layout, "Diffuse reflectance:", self, "dsr")
            newrow(layout, "Visible reflectance:", self, "vr") 
            newrow(layout, "Thermal emmisivity:", self, "the")
            newrow(layout, "Conductivity:", self, "tc")
            newrow(layout, "Material spacing:", self, "sme")
            newrow(layout, "Material diameter:", self, "smd")
            newrow(layout, "Distance:", self, "sgd")
            newrow(layout, "Top mult.:", self, "tom")
            newrow(layout, "Bottom mult.:", self, "bom")
            newrow(layout, "Left milt.:", self, "lom")
            newrow(layout, "Right mult.:", self, "rom")
            newrow(layout, "Resolution angle:", self, "ta")
    
    def valid(self):
        if not self.outputs["Outer Layer"].links or not self.inputs["Layer"].links or not self.inputs["Control"].links:
            nodecolour(self, 1)
        else:
            nodecolour(self, 0)
            
    def update(self):
        socklink2(self.outputs['Outer Layer'], self.id_data)
        self.valid()
        
    def ep_write(self, ln):
        for material in bpy.data.materials:
            if self.id_data == material.envi_nodes:
                break
            
        params = ('Name', 'Reflected Beam Transmittance Accounting Method', 'Diffuse Solar Reflectance', 'Diffuse Visible Reflectance',
                  'Thermal Hemispherical Emissivity', 'Conductivity (W/m-K)', 'Screen Material Spacing (m)', 'Screen Material Diameter (m)',
                  'Screen-to-Glass Distance (m)', 'Top Opening Multiplier', 'Bottom Opening Multiplier', 'Left-Side Opening Multiplier', 
                  'Right-Side Opening Multiplier', 'Angle of Resolution for Output Map (deg)')
        
        paramvs = ['{}-layer-{}'.format(material.name, ln), self.rb] + ['{:.3f}'.format(p) for p in (self.dsr, self.vr, self.the, self.tc, 0.001 * self.sme, 0.001 * self.smd,
                   0.001 * self.sgd, self.tom, self.bom, self.lom, self.rom)] + [self.ta]
  
        return epentry('WindowMaterial:Screen', params, paramvs) + self.inputs['Control'].links[0].from_node.ep_write(ln)
    
class No_En_Mat_Bl(Node, EnViMatNodes):
    '''Node defining an EnVi window blind'''
    bl_idname = 'No_En_Mat_Bl'
    bl_label = 'EnVi blind'
    
    so: EnumProperty(items = [("0", "Horizontal", "Select from database"), 
                                ("1", "Vertical", "Define custom material properties")],
                                name = "", description = "Slat orientation", default = '0')
    sw: FloatProperty(name = "mm", description = "Slat width", min = 0.1, max = 1000, default = 25)
    ss: FloatProperty(name = "mm", description = "Slat separation", min = 0.1, max = 1000, default = 20)
    st: FloatProperty(name = "mm", description = "Slat thickness", min = 0.1, max = 1000, default = 50)
    sa: FloatProperty(name = "deg", description = "Slat angle", min = 0.0, max = 90, default = 45)
    stc: FloatProperty(name = "W/m.K", description = "Slat conductivity", min = 0.01, max = 100, default = 10)
    sbst: FloatProperty(name = "", description = "Slat beam solar transmittance", min = 0.0, max = 1, default = 0.0)
    fbst: FloatProperty(name = "", description = "Front Side Slat beam solar reflectance", min = 0.0, max = 1, default = 0.8)
    bbst: FloatProperty(name = "", description = "Back Side Slat beam solar reflectance", min = 0.0001, max = 10, default = 0.8)
    sdst: FloatProperty(name = "m", description = "Slat diffuse solar transmittance", min = 0.0, max = 1, default = 0.0)
    fdsr: FloatProperty(name = "", description = "Front Side Slat diffuse solar reflectance", min = 0.0, max = 1, default = 0.8)
    bdsr: FloatProperty(name = "", description = "Back Side Slat diffuse solar reflectance", min = 0.0, max = 1, default = 0.8)
    sbvt: FloatProperty(name = "", description = "Slat beam visible transmittance", min = 0.0, max = 1, default = 0.0)
    fbvr: FloatProperty(name = "", description = "Front Side Slat beam visible reflectance", min = 0.0, max = 1, default = 0.7)
    bbvr: FloatProperty(name = "", description = "Back Side Slat beam visible reflectance", min = 0.0, max = 1, default = 0.7)
    sdvt: FloatProperty(name = "", description = "Slat diffuse visible transmittance", min = 0.0, max = 1, default = 0.0)
    fdvr: FloatProperty(name = "", description = "Front Side Slat diffuse visible reflectance", min = 0.0, max = 1, default = 0.7)
    bdvr: FloatProperty(name = "", description = "Back Side Slat diffuse visible reflectance", min = 0.0, max = 1, default = 0.7)
    sit: FloatProperty(name = "", description = "Slat Infrared hemispherical transmittance", min = 0.0, max = 1, default = 0.0)
    sfie: FloatProperty(name = "", description = "Front Side Slat Infrared hemispherical emissivity", min = 0.0, max = 1, default = 0.9)
    sbie: FloatProperty(name = "", description = "Back Side Slat Infrared hemispherical emissivity", min = 0.0, max = 1, default = 0.9)
    bgd: FloatProperty(name = "mm", description = "Blind-to-glass distance", min = 0.1, max = 1000, default = 50)
    tom: FloatProperty(name = "", description = "Blind top opening multiplier", min = 0.0, max = 1, default = 0.0)
    bom: FloatProperty(name = "", description = "Blind bottom opening multiplier", min = 0.0, max = 1, default = 0.0)
    lom: FloatProperty(name = "", description = "Blind left-side opening multiplier", min = 0.0, max = 1, default = 0.5)
    rom: FloatProperty(name = "", description = "Blind right-side opening multiplier", min = 0.0, max = 1, default = 0.5)
    minsa: FloatProperty(name = "deg", description = "Minimum slat angle", min = 0.0, max = 90, default = 0.0)
    maxsa: FloatProperty(name = "deg", description = "Maximum slat angle", min = 0.0, max = 90, default = 0.0)
    
    def init(self, context):
        self.outputs.new('envi_sl_sock', 'Layer')
        self.inputs.new('envi_sc_sock', 'Control')
        self.inputs.new('envi_sl_sock', 'Layer')
        
    def draw_buttons(self, context, layout):
        if self.outputs['Layer'].links:
            newrow(layout, "Slat orient.:", self, "so")
            newrow(layout, "Slat width:", self, "sw")
            newrow(layout, "Slat sep.:", self, "ss")
            newrow(layout, "Slat thick.:", self, "st")
            newrow(layout, "Slat angle:", self, "sa") 
            newrow(layout, "Slat cond.:", self, "stc")
            newrow(layout, "Slat beam trans.:", self, "sbst")
            newrow(layout, "Front beam trans.:", self, "fbst")
            newrow(layout, "Back beam trans.:", self, "bbst")
            newrow(layout, "Slat diff. trans.:", self, "sdst")
            newrow(layout, "Front diff. reflec.:", self, "fdsr")
            newrow(layout, "Back diff. reflec.:", self, "bdsr")
            newrow(layout, "Slat beam trans.:", self, "sbvt")
            newrow(layout, "Front beam vis. reflec.:", self, "fbvr")
            newrow(layout, "Back beam vis. reflec.:", self, "bbvr")
            newrow(layout, "Slat diff. vis. trans.:", self, "sdvt")
            newrow(layout, "Front diff. vis. reflec.:", self, "fdvr")
            newrow(layout, "Back diff. vis. reflec.:", self, "bdvr")
            newrow(layout, "Infra. trans.:", self, "sit") 
            newrow(layout, "Front infra. emiss.:", self, "sfie")
            newrow(layout, "Back infra. emiss.:", self, "sbie")
            newrow(layout, "Glass dist.:", self, "bgd")
            newrow(layout, "Top mult.:", self, "tom")
            newrow(layout, "Bottom mult.:", self, "bom")
            newrow(layout, "Left mult.:", self, "lom")
            newrow(layout, "Right mult.:", self, "rom")
            newrow(layout, "Min ang.:", self, "minsa")
            newrow(layout, "Max ang.:", self, "maxsa")

    def valid(self):
        if not self.outputs["Layer"].links or not self.inputs["Layer"].links:
            nodecolour(self, 1)
        else:
            nodecolour(self, 0)
            
    def update(self):
        socklink2(self.outputs['Layer'], self.id_data)
        self.valid()
        
    def ep_write(self, ln):
        for material in bpy.data.materials:
            if self.id_data == material.envi_nodes:
                break
        params = ('Name', 'Slat orientation', 'Slat width (m)', 'Slat separation (m)', 'Slat thickness (m)', 'Slat angle (deg)', 'Slat conductivity (W/m.K)',
                  'Slat beam solar transmittance', 'Front Side Slat beam solar reflectance', 'Back Side Slat beam solar reflectance', 'Slat diffuse solar transmittance', 
                  'Front Side Slat diffuse solar reflectance', 'Back Side Slat diffuse solar reflectance', 'Slat beam visible transmittance', 'Front Side Slat beam visible reflectance',
                  'Back Side Slat beam visible reflectance', 'Slat diffuse visible transmittance', "Front Side Slat diffuse visible reflectance", "Back Side Slat diffuse visible reflectance",
                  "Slat Infrared hemispherical transmittance", "Front Side Slat Infrared hemispherical emissivity", "Back Side Slat Infrared hemispherical emissivity", "Blind-to-glass distance",
                  "Blind top opening multiplier", "Blind bottom opening multiplier", "Blind left-side opening multiplier", "Blind right-side opening multiplier", "Minimum slat angle", "Maximum slat angle")
        paramvs = ['{}-layer-{}'.format(material.name, ln), ('Horizontal', 'Vertical')[int(self.so)]] + ['{:.3f}'.format(p) for p in (0.001 * self.sw, 0.001 * self.ss, 0.001 * self.st, self.sa, self.stc, self.sbst, self.fbst, self.bbst, self.sdst, self.fdsr, self.bdsr, self.sbvt,
                   self.fbvr, self.bbvr, self.sdvt, self.fdvr, self.bdvr, self.sit, self.sfie, self.sbie, 0.001 * self.bgd, self.tom, self.bom, self.lom, self.rom, self.minsa, self.maxsa)]
  
        return epentry('WindowMaterial:Blind', params, paramvs) + self.inputs['Control'].links[0].from_node.ep_write(ln)
    
class No_En_Mat_SG(Node, EnViMatNodes):
    '''Node defining the EnVi switchable glazing layer'''
    bl_idname = 'No_En_Mat_SG'
    bl_label = 'EnVi switchable glazing layer'
    
    layer: EnumProperty(items = [("0", "Database", "Select from database"), 
                                        ("1", "Custom", "Define custom material properties")], 
                                        name = "", description = "Composition of the layer", default = "0")
    materialtype: EnumProperty(items = envi_layertype, name = "", description = "Layer material type")
    mats = [((mat, mat, 'Layer material')) for mat in envi_materials().glass_dat.keys()]
    material: EnumProperty(items = mats, name = "", description = "Glass material")
    thi: FloatProperty(name = "mm", description = "Thickness (mm)", min = 0.1, max = 10000, default = 100)
    tc: FloatProperty(name = "W/m.K", description = "Thermal Conductivity (W/m.K)", min = 0.1, max = 10000, default = 100)
    stn: FloatProperty(name = "", description = "Solar normal transmittance", min = 0, max = 1, default = 0.7)
    fsn: FloatProperty(name = "", description = "Solar front normal reflectance", min = 0, max = 1, default = 0.7)
    bsn: FloatProperty(name = "", description = "Solar back normal reflectance", min = 0, max = 1, default = 0.7)
    vtn: FloatProperty(name = "", description = "Visible Transmittance at Normal Incidence", min = 0, max = 1, default = 0.7)
    fvrn: FloatProperty(name = "", description = "Front Side Visible Reflectance at Normal Incidence", min = 0, max = 1, default = 0.7)
    bvrn: FloatProperty(name = "", description = "Back Side Visible Reflectance at Normal Incidence", min = 0, max = 1, default = 0.7)
    itn: FloatProperty(name = "", description = "Infrared Transmittance at Normal Incidence", min = 0, max = 1, default = 0.7)
    fie: FloatProperty(name = "", description = "Front Side Infrared Hemispherical Emissivity'", min = 0, max = 1, default = 0.7)
    bie: FloatProperty(name = "", description = "Back Side Infrared Hemispherical Emissivity", min = 0, max = 1, default = 0.7)
    diff: BoolProperty(name = "", description = "Diffusing", default = 0)
    
    def init(self, context):
        self.inputs.new('envi_sc_sock', 'Control')
        self.inputs.new('envi_sgl_sock', 'Layer')
        self.outputs.new('envi_sgl_sock', 'Layer')
        
    def draw_buttons(self, context, layout):
        if self.inputs['Layer'].links:
            newrow(layout, "Class:", self, "layer")
            if self.layer == '0':
                newrow(layout, "Material:", self, "material")
                newrow(layout, "Thickness:", self, "thi") 
            else:
                newrow(layout, "Conductivity:", self, "tc")
                newrow(layout, "Thickness:", self, "thi")
                newrow(layout, "STN:", self, "stn")
                newrow(layout, "FSN:", self, "fsn")
                newrow(layout, "BSN:", self, "bsn")
                newrow(layout, "VTN:", self, "vtn")
                newrow(layout, "FVRN:", self, "fvrn")
                newrow(layout, "BVRN:", self, "bvrn")
                newrow(layout, "ITN:", self, "itn")
                newrow(layout, "FIE:", self, "fie")
                newrow(layout, "BIE:", self, "bie")
                newrow(layout, "Diffuse:", self, "diff")

    def update(self):
        self.valid()
    
    def valid(self):
        if not self.outputs["Layer"].links:
            nodecolour(self, 1)
        else:
            nodecolour(self, 0)
     
    def ep_write(self, ln):
        for material in bpy.data.materials:
            if self.id_data == material.envi_nodes:
                break

        layer_name = '{}-layer-{}'.format(material.name, ln)
        params = ('Name', 'Optical Data Type', 'Window Glass Spectral Data Set Name', 'Thickness (m)', 'Solar Transmittance at Normal Incidence', 'Front Side Solar Reflectance at Normal Incidence',
                  'Back Side Solar Reflectance at Normal Incidence', 'Visible Transmittance at Normal Incidence', 'Front Side Visible Reflectance at Normal Incidence', 'Back Side Visible Reflectance at Normal Incidence',
                  'Infrared Transmittance at Normal Incidence', 'Front Side Infrared Hemispherical Emissivity', 'Back Side Infrared Hemispherical Emissivity', 'Conductivity (W/m-K)',
                  'Dirt Correction Factor for Solar and Visible Transmittance', 'Solar Diffusing')

        if self.layer == '0':
            matlist = list(envi_mats.matdat[self.material])
            paramvs = [layer_name] + matlist[1:3] + [self.thi] + ['{:.3f}'.format(float(sm)) for sm in matlist[4:-1]] + [1, ('No', 'Yes')[matlist[-1]]]
            
        else:
            paramvs = ['{}-layer-{}'.format(material.name, ln), 'SpectralAverage', '', self.thi/1000, '{:.3f}'.format(self.stn), '{:.3f}'.format(self.fsn), '{:.3f}'.format(self.bsn), 
                       '{:.3f}'.format(self.vtn), '{:.3f}'.format(self.fvrn), '{:.3f}'.format(self.bvrn), '{:.3f}'.format(self.itn),
                       '{:.3f}'.format(self.fie), '{:.3f}'.format(self.bie), '{:.3f}'.format(self.tc), 1, ('No', 'Yes')[self.diff]]

        return epentry("WindowMaterial:Glazing", params, paramvs)   + self.inputs['Control'].links[0].from_node.ep_write(ln) 
          
class No_En_Mat_ShC(Node, EnViMatNodes):
    '''Node defining an EnVi window shade control'''
    bl_idname = 'No_En_Mat_ShC'
    bl_label = 'EnVi shade control'
        
    ttuple = ("Alwayson", "Alwaysoff", "OnIfScheduleAllows", "OnIfHighSolarOnWindow", "OnIfHighHorizontalSolar", 
              "OnIfHighOutdoorAirTemperature", 
              "OnIfHighZoneAirTemperature", "OnIfHighZoneCooling", "OnIfHighGlare", "MeetDaylightIlluminanceSetpoint",
              "OnNightIfLowOutdoorTempAndOffDay", "OnNightIfLowInsideTempAndOffDay", "OnNightIfHeatingAndOffDay",
              "OnNightIfLowOutdoorTempAndOnDayIfCooling", "OnNightIfHeatingAndOnDayIfCooling", 
              "OffNightAndOnDayIfCoolingAndHighSolarOnWindow", "OnNightAndOnDayIfCoolingAndHighSolarOnWindow", 
              "OnIfHighOutdoorAirTempAndHighSolarOnWindow", "OnIfHighOutdoorAirTempAndHighHorizontalSolar")
    
    def type_menu(self, context):
        try:
            if self.outputs['Control'].links[0].to_node.bl_idname == 'envi_screen_node':
                return [(self.ttuple[t], self.ttuple[t], self.ttuple[t]) for t in (0, 1, 2)]
            elif self.outputs['Control'].links[0].to_node.bl_idname in ('envi_bl_node', 'envi_sl_node'):
                return [(self.ttuple[t], self.ttuple[t], self.ttuple[t]) for t in (0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18)]
            else:
                return [(t, t, t) for t in self.ttuple]
        except Exception as e:
            logentry('Shade control error {}'.format(e))
            return [('None', 'None', 'None')]

    ctype: EnumProperty(items = type_menu, name = "", description = "Shading device")
    sp: FloatProperty(name = "", description = "Setpoint (W/m2, W or deg C)", min = 0.0, max = 1000, default = 20)
    sac: EnumProperty(items = [("FixedSlatAngle", "Always on", "Shading component"), 
                                ("ScheduledSlatAngle", "OnIfHighOutdoorAirTempAndHighSolarOnWindow", "Switchable glazing component"),
                                ("BlockBeamSolar", "OnIfHighOutdoorAirTempAndHighHorizontalSolar", "Switchable glazing component")
                                ],
                                name = "", description = "Shading device", default = 'FixedSlatAngle')
    sp2: FloatProperty(name = "", description = "Setpoint 2 (W/m2, W or deg C)", min = 0.0, max = 1000, default = 20)
      
    def init(self, context):
        self.outputs.new('envi_sgl_sock', 'Layer')
        self.outputs['Layer'].hide = True
        self.outputs.new('So_En_Mat_ShC', 'Control')
        self.inputs.new('So_En_Sched', 'Schedule')

    def draw_buttons(self, context, layout):
        newrow(layout, "Shading device:", self, 'ctype')
        if self.ctype not in ('Always on', 'Always off', 'OnIfScheduleAllows', 'OnIfHighGlare', 'DaylightIlluminance'):
            newrow(layout, "Set-point", self, 'sp')
        if self.outputs['Control'].links and self.outputs['Control'].links[0].to_node.bl_idname == 'envi_blind_node':
            newrow(layout, 'Slat angle:', self, 'sac')
           
    def valid(self):
        if not self.outputs["Control"].links:
            nodecolour(self, 1)
        else:
            nodecolour(self, 0)
            
    def update(self):
        socklink2(self.outputs['Control'], self.id_data)
#        socklink2(self.outputs['Layer'], self.id_data)
#        if self.outputs['Control'].links and self.outputs['Control'].links[0].to_node.bl_idname == 'envi_tl_node':
#            self.outputs['Layer'].hide = False
#        else:
#            for link in self.outputs['Shade'].links:
#                self.id_data.links.remove(link)
#            self.outputs['Layer'].hide = True
        self.valid()
        
    def ep_write(self, ln):
        for material in bpy.data.materials:
            if self.id_data == material.envi_nodes:
                break
            
        if self.outputs['Control'].links[0].to_node.bl_idname == 'envi_screen_node':
            st = 'ExteriorScreen' 
        elif self.outputs['Control'].links[0].to_node.bl_idname == 'envi_bl_node':
            if ln == 0:
                st = 'ExteriorBlind'
            elif self.outputs['Control'].links[0].to_node.inputs['Layer'].links:
                st = 'BetweenGlassBlind'
            else:
                st = 'InteriorBlind'
        elif self.outputs['Control'].links[0].to_node.bl_idname == 'envi_sl_node':
            if ln == 0:
                st = 'ExteriorShade'
            elif self.outputs['Control'].links[0].to_node.inputs['Layer'].links:
                st = 'BetweenGlassShade'
            else:
                st = 'InteriorShade'
        else:
            st = 'SwitchableGlazing'
        
        (scs, sn) = ('Yes', '{}-shading-schedule-{}'.format(material.name, ln)) if self.inputs['Schedule'].links else ('No', '')
                
        params = ('Name', 'Shading Type', 'Construction with Shading Name', 'Shading Control Type', 'Schedule Name', 'Setpoint (W/m2, W or deg C)', 'Shading Control Is Scheduled',
                  'Glare Control Is Active', 'Shading Device Material Name', 'Type of Slat Angle Control for Blinds', 'Slat Angle Schedule Name')
        paramvs = ('{}-shading-control'.format(material.name), st, '{}-shading'.format(material.name), self.ctype, sn, self.sp, scs, 'No', '', self.sac, '')
        return epentry('WindowProperty:ShadingControl', params, paramvs)

class No_En_Mat_PV(Node, EnViMatNodes):
    '''Node defining an EnVi photovoltaic module'''
    bl_idname = 'No_En_Mat_PV'
    bl_label = 'EnVi PV'
    
    def pv_update(self, context):
        pass
        
    sandia_dict = {}
    l = -40
    
    def ret_e1dmenu(self, context):  
        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'EPFiles', '{}'.format('PV_database.json')), 'r') as pv_json:
            e1ddict = json.loads(pv_json.read())
        return [(p, p, '{} module'.format(p)) for p in e1ddict]
            
    def ret_sandiamenu(self, context):  
        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'EPFiles', '{}'.format('SandiaPVdata.json')), 'r') as sandia_json:
            sandiadict = json.loads(sandia_json.read())
        return [(p, p, '{} module'.format(p)) for p in sandiadict]
    
    def ret_e1ddict(self):
        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'EPFiles', '{}'.format('PV_database.json')), 'r') as pv_json:
            e1ddict = json.loads(pv_json.read())
            return e1ddict
        
    def save_e1ddict(self):
        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'EPFiles', '{}'.format('PV_database.json')), 'r') as pv_json:
            e1ddict = json.loads(pv_json.read())
            
        e1ddict[self.pv_name] = [self.scc, self.ocv, self.mv, self.mc, self.tcscc, self.tcocv, self.mis, self.ctnoct, self.mod_area]
        self.e1ddict = e1ddict
#        self.e1dmenu = self.pv_name
        
        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'EPFiles', '{}'.format('PV_database.json')), 'w') as e1d_jfile:
            e1d_jfile.write(json.dumps(e1ddict))  
    
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'EPFiles', '{}'.format('PV_database.json')), 'r') as pv_json:
        e1ddict = json.loads(pv_json.read())
        
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'EPFiles', '{}'.format('SandiaPVdata.json')), 'r') as sandia_json:
        sandiadict = json.loads(sandia_json.read())
                
    ct: EnumProperty(items = [("0", "Crystalline", "Do not model reflected beam component"), 
                               ("1", "Amorphous", "Model reflectred beam as beam")], 
                                name = "", description = "Photovoltaic Type", default = "0")

    hti: EnumProperty(items = [("Decoupled", "Decoupled", "Decoupled"), 
                               ("DecoupledUllebergDynamic", "Ulleberg", "DecoupledUllebergDynamic"),
                               ("IntegratedSurfaceOutsideFace", "SurfaceOutside", "IntegratedSurfaceOutsideFace"),
                               ("IntegratedTranspiredCollector", "Transpired", "IntegratedTranspiredCollector"),
                               ("IntegratedExteriorVentedCavity", "ExteriorVented", "IntegratedExteriorVentedCavity"),
                               ("PhotovoltaicThermalSolarCollector", "PVThermal", "PhotovoltaicThermalSolarCollector")], 
                                name = "", description = "Conversion Efficiency Input Mode'", default = "Decoupled")
        
    pp: EnumProperty(items = [("0", "Simple", "Do not model reflected beam component"), 
                               ("1", "One-Diode", "Model reflectred beam as beam"),
                               ("2", "Sandia", "Model reflected beam as diffuse")], 
                                name = "", description = "Photovoltaic Performance Object Type", default = "0", update = pv_update)
    pv_name: StringProperty(name = '', description = 'Name of the custom PV model', default = '')
    e1dmenu: EnumProperty(items = ret_e1dmenu, name = "", description = "Module type")
    smenu: EnumProperty(items = ret_sandiamenu, name = "", description = "Module type")
    mod_area: FloatProperty(name = "m2", description = "PV module area", min = 0.1, default = 5)
    pvsa: FloatProperty(name = "%", description = "Fraction of Surface Area with Active Solar Cells", min = 50, max = 100, default = 90)
    aa: FloatProperty(name = "m2", description = "Active area", min = 0.1, max = 10000, default = 5)
    eff: FloatProperty(name = "%", description = "Visible reflectance", min = 0.0, max = 100, default = 20)
    ssp: IntProperty(name = "", description = "Number of series strings in parallel", min = 1, max = 100, default = 1)
    mis: IntProperty(name = "", description = "Number of modules in series", min = 1, max = 100, default = 1)
    cis: IntProperty(name = "", description = "Number of cells in series", min = 1, max = 100, default = 36) 
    tap: FloatProperty(name = "", description = "Transmittance absorptance product", min = -1, max = 1, default = 0.9)
    sbg: FloatProperty(name = "eV", description = "Semiconductor band-gap", min = 0.1, max = 5, default = 1.12)
    sr: FloatProperty(name = "W", description = "Shunt resistance", min = 1, default = 1000000)
    scc: FloatProperty(name = "Amps", description = "Short circuit current", min = 1, max = 1000, default = 25)
    sgd: FloatProperty(name = "mm", description = "Screen to glass distance", min = 1, max = 1000, default = 25)
    ocv: FloatProperty(name = "V", description = "Open circuit voltage", min = 0.0, max = 100, default = 60)
    rt: FloatProperty(name = "C", description = "Reference temperature", min = 0, max = 40, default = 25)
    ri: FloatProperty(name = "W/m2", description = "Reference insolation", min = 100, max = 2000, default = 1000)
    mc: FloatProperty(name = "Amps", description = "Module current at maximum power", min = 1, max = 10, default = 5.6)
    mv: FloatProperty(name = "V", description = "Module voltage at maximum power", min = 0.0, max = 75, default = 17)
    tcscc: FloatProperty(name = "A/K", description = "Temperature Coefficient of Short Circuit Current", precision = 5, min = 0.00001, max = 0.01, default = 0.002)
    tcocv: FloatProperty(name = "V/K", description = "Temperature Coefficient of Open Circuit Voltage", precision = 5, min = -0.5, max = 0, default = -0.1)
    atnoct: FloatProperty(name = "C", description = "Reference ambient temperature", min = 0, max = 40, default = 20)
    ctnoct: FloatProperty(name = "C", description = "Nominal Operating Cell Temperature Test Cell Temperature", min = 0, max = 60, default = 45)
    inoct: FloatProperty(name = "W/m2", description = "Nominal Operating Cell Temperature Test Insolation", min = 100, max = 2000, default = 800)
    hlc: FloatProperty(name = "W/m2.K", description = "Module heat loss coefficient", min = 0.0, max = 50, default = 30)
    thc: FloatProperty(name = " J/m2.K", description = " Total Heat Capacity", min = 10000, max = 100000, default = 50000)
    
    def init(self, context):
        mat = bpy.data.materials[self.id_data.name]
        if not mat.vi_params.get('enparams'):
            mat.vi_params['enparams'] = {'area' : -1} 
        self['area'] = 0
        self.inputs.new('So_En_Mat_PVG', 'PV Generator')
        self.inputs.new('So_En_Sched', 'PV Schedule')
#        self.inputs['PV Schedule'].hide = True
        self.outputs.new('So_En_Mat_PV', 'PV')
        
    def draw_buttons(self, context, layout):
        mat = bpy.data.materials[self.id_data.name]
        row = layout.row()
#        print(self['area'], mat.vi_params['enparams']['area'])
#        if self['area'] != mat.vi_params['enparams']['area']:            
        row.operator('node.pv_area', text = "Area Calc")
        try:                        
            row.label(text = '{:.2f} m2'.format(self['area'])) 
        except: 
            row.label(text = 'Area  = N/A')
        #else:
#            row.label(text = 'Area: {:.2f}m2'.format(self['area']))
        newrow(layout, "Heat transfer:", self, "hti")
        newrow(layout, 'Type:', self, 'pp')
        
        if self.pp == '0':
            newrow(layout, "PV area ratio:", self, "pvsa")
            newrow(layout, "Efficiency:", self, "eff")
            
        elif self.pp == '1':
            newrow(layout, "Model:", self, "e1dmenu")

            if self.e1dmenu == 'Custom':
                newrow(layout, "Name:", self, "pv_name")
                newrow(layout, "Module area:", self, "mod_area")
                newrow(layout, "Cell type:", self, "ct")
                newrow(layout, "Cells:", self, "cis")                   
                newrow(layout, "SCC:", self, "scc") 
                newrow(layout, "OCV:", self, "ocv")                
                newrow(layout, "Max power I:", self, "mc")
                newrow(layout, "Max power V:", self, "mv")
                newrow(layout, "TCSCC:", self, "tcscc")
                newrow(layout, "TCOCV:", self, "tcocv")                
                newrow(layout, "Nominal temp.:", self, "ctnoct")
                if self.pv_name:
                    row=layout.row()
                    row.operator('node.pv_save', text = "PV Save")
            newrow(layout, "Trans*absorp:", self, "tap")
            newrow(layout, "Band gap:", self, "sbg")
            newrow(layout, "Shunt:", self, "sr")
            newrow(layout, "Ref. temp.:", self, "rt")
            newrow(layout, "Ref. insol.:", self, "ri")
            newrow(layout, "Ambient temp.:", self, "atnoct")
            newrow(layout, "Test Insolation:", self, "inoct")
            newrow(layout, "Heat loss coeff.:", self, "hlc")
            newrow(layout, "Heat capacity:", self, "thc")
                             
        if self.pp == '2':
            newrow(layout, 'Model:', self, 'smenu')
    
    def update(self):
        if len(self.outputs) + len(self.inputs) == 3:
            if self.outputs['PV'].links and not self.inputs['PV Generator'].links:
                nodecolour(self, 1)         
            else:
                nodecolour(self, 0)                     
    
    def ep_write(self, sn, area):
        self['matname'] = get_mat(self, 1).name
        marea = (area, (self.mod_area, self.e1ddict[self.e1dmenu][8])[self.e1dmenu != 'Custom'], self.sandiadict[self.smenu][0])[int(self.pp)]
        
        if marea /area > 1.01:
            logentry("No PV data exported as the face {} area is smaller than the module area".format(sn))
            return '! No PV data exported as the face {} area is smaller than the module area\n\n'.format(sn)
            
        params = ('Name', 'Surface Name', 'Photovoltaic Performance Object Type', 
                  'Module Performance Name', 'Heat Transfer Integration Mode', 
                  'Number of Series Strings in Parallel', 'Number of Modules in Series')
                
        paramvs = ['{}-pv'.format(sn), sn, 
                   ('PhotovoltaicPerformance:Simple', 'PhotovoltaicPerformance:EquivalentOne-Diode', 'PhotovoltaicPerformance:Sandia')[int(self.pp)], '{}-pv-performance'.format(sn),
                   self.hti, self.ssp, area/marea]

        ep_text = epentry('Generator:Photovoltaic', params, paramvs)
        
        if self.pp == '0':
            params = ('Name', 'Fraction of Surface Area with Active Solar Cell', 'Conversion Efficiency Input Mode', 'Value for Cell Efficiency if Fixed', 'Efficiency Schedule Name')
            paramvs = ('{}-pv-performance'.format(sn), self.pvsa * 0.01, ('Fixed', 'Scheduled')[len(self.inputs['PV Schedule'].links)], self.eff * 0.01, ('', '{}-pv-performance-schedule'.format(sn))[len(self.inputs['PV Schedule'].links)])
            ep_text += epentry('PhotovoltaicPerformance:Simple', params, paramvs)
            
            if self.inputs['PV Schedule'].links:
                ep_text += self.inputs['PV Schedule'].links[0].from_node.epwrite('{}-pv-performance-schedule'.format(sn), 'Fraction')
            
        elif self.pp == '1':
            params = ('Name', 'Cell type', 'Number of Cells in Series', 'Active Area (m2)', 'Transmittance Absorptance Product',
                      'Semiconductor Bandgap (eV)', 'Shunt Resistance (ohms)', 'Short Circuit Current (A)', 'Open Circuit Voltage (V)',
                      'Reference Temperature (C)', 'Reference Insolation (W/m2)', 'Module Current at Maximum Power (A)', 
                      'Module Voltage at Maximum Power (V)', 'Temperature Coefficient of Short Circuit Current (A/K)',
                      'Temperature Coefficient of Open Circuit Voltage (V/K)', 'Nominal Operating Cell Temperature Test Ambient Temperature (C)',
                      'Nominal Operating Cell Temperature Test Cell Temperature (C)', 'Nominal Operating Cell Temperature Test Insolation (W/m2)',
                      'Module Heat Loss Coefficient (W/m2-K)', 'Total Heat Capacity (J/m2-K)')
            paramvs = ('{}-pv-performance'.format(sn), ('CrystallineSilicon', 'AmorphousSilicon')[int(self.ct)], (self.cis, self.e1ddict[self.e1dmenu][6])[self.e1dmenu != 'Custom'], (self.mod_area, self.e1ddict[self.e1dmenu][8])[self.e1dmenu != 'Custom'],
                       self.tap, self.sbg, self.sr, (self.scc, self.e1ddict[self.e1dmenu][0])[self.e1dmenu != 'Custom'], (self.ocv, self.e1ddict[self.e1dmenu][1])[self.e1dmenu != 'Custom'],
                       self.rt, self.ri, (self.mc, self.e1ddict[self.e1dmenu][3])[self.e1dmenu != 'Custom'], (self.mv, self.e1ddict[self.e1dmenu][2])[self.e1dmenu != 'Custom'],
                       (self.tcscc, self.e1ddict[self.e1dmenu][4])[self.e1dmenu != 'Custom'], (self.tcocv, self.e1ddict[self.e1dmenu][5])[self.e1dmenu != 'Custom'],
                       self.atnoct, (self.ctnoct, self.e1ddict[self.e1dmenu][7])[self.e1dmenu != 'Custom'], self.inoct, self.hlc, self.thc)
            ep_text += epentry('PhotovoltaicPerformance:EquivalentOne-Diode', params, paramvs)   
            
        elif self.pp == '2':
            params = ('Name', 'Field Active Area Acoll (m2), single module', 'NcellSer (unitless)', 'NparSerCells',
                      'Isc0 (Amps)', 'Voc0 (Volts)', 'Imp0 (Amps)', 'Vmp0 (Volts)', 'aIsc (1/ degC)', 'aImp (1/ degC)',
                      'C0 (unitless)', 'C1 (unitless)', 'BVoc0 (Volts/degC)', 'mBVoc (Volts/degC)', 'BVmp0 (Volts/degC)',
                      'mBVmp (Volts/degC)', 'Diode Factor (n) (Unitless)', 'C2 (Unitless)', 'C3 (Unitless)',
                      'A0 (Unitless)', 'A1 (Unitless)', 'A2 (Unitless)', 'A3 (Unitless)', 'A4 (Unitless)', 'B0 (Unitless)',
                      'B1 (Unitless)', 'B2 (Unitless)', 'B3 (Unitless)', 'B4 (Unitless)', 'B5 (Unitless)', 'dT0 (degC)',
                      'fd (Unitless)','a (Unitless)','b (Unitless)','C4 (Unitless)','C5 (Unitless)','Ix0 (Amps)','Ixx0 (Amps)',
                      'C6 (Unitless)','C7 (Unitless)')
            paramvs = ['{}-pv-performance'.format(sn)] + self.sandiadict[self.smenu]
            ep_text += epentry('PhotovoltaicPerformance:Sandia', params, paramvs)
        return ep_text                    
    
class No_En_Mat_PVG(Node, EnViMatNodes):
    '''Node defining an EnVi photovoltaic generator'''
    bl_idname = 'No_En_Mat_PVG'
    bl_label = 'EnVi PV Generator'    
    
    it: EnumProperty(items = [('0', 'Simple', 'Simple Inverter')], name = "", description = "Inverter type")
    ie: FloatProperty(name = "%", description = "Inverter efficiency (%)", min = 0.0, max = 100, default = 90)
    rf: FloatProperty(name = "", description = "Fraction", min = 0, max = 1, default = 0.1)
    maxload: IntProperty(name = "", description = "Max Power (W)", max = 10, default = 1000)
    
#    def elc_write()

    def init(self, context):
        self.inputs.new('So_En_Sched', 'Schedule')
#        self.inputs['PV Schedule'].hide = True
        self.outputs.new('So_En_Mat_PVG', 'PV Gen')
        
    def draw_buttons(self, context, layout):
        newrow(layout, "Inverter type:", self, "it")
        newrow(layout, "Inverter efficiency:", self, "ie")
        newrow(layout, 'Radiative fraction:', self, 'rf')
        
class No_En_Mat_Sched(Node, EnViNodes):
    '''Node describing a schedule'''
    bl_idname = 'No_En_Mat_Sched'
    bl_label = 'Schedule'
    bl_icon = 'SOUND'

    def tupdate(self, context):
        try:
            err = 0
            if self.t2 <= self.t1 and self.t1 < 365:
                self.t2 = self.t1 + 1
                if self.t3 <= self.t2 and self.t2 < 365:
                    self.t3 = self.t2 + 1
                    if self.t4 != 365:
                        self.t4 = 365

            tn = (self.t1, self.t2, self.t3, self.t4).index(365) + 1
            if max((self.t1, self.t2, self.t3, self.t4)[:tn]) != 365:
                err = 1
            if any([not f for f in (self.f1, self.f2, self.f3, self.f4)[:tn]]):
                err = 1
            if any([not u or '. ' in u or len(u.split(';')) != len((self.f1, self.f2, self.f3, self.f4)[i].split(' ')) for i, u in enumerate((self.u1, self.u2, self.u3, self.u4)[:tn])]):
                err = 1

            for f in (self.f1, self.f2, self.f3, self.f4)[:tn]:
                for fd in f.split(' '):
                    if not fd or (fd and fd.upper() not in ("ALLDAYS", "WEEKDAYS", "WEEKENDS", "MONDAY", "TUESDAY", "WEDNESDAY", "THURSDAY", "FRIDAY", "SATURDAY", "SUNDAY", "ALLOTHERDAYS")):
                        err = 1

            for u in (self.u1, self.u2, self.u3, self.u4)[:tn]:
                for uf in u.split(';'):
                    for ud in uf.split(','):
                        if len(ud.split()[0].split(':')) != 2 or int(ud.split()[0].split(':')[0]) not in range(1, 25) or len(ud.split()[0].split(':')) != 2 or not ud.split()[0].split(':')[1].isdigit() or int(ud.split()[0].split(':')[1]) not in range(0, 60):
                            err = 1
            nodecolour(self, err)

        except:
            nodecolour(self, 1)

    file: EnumProperty(name = '', items = [("0", "None", "No file"), ("1", "Select", "Select file"), ("2", "Generate", "Generate file")], default = '0')
    select_file: StringProperty(name="", description="Name of the variable file", default="", subtype="FILE_PATH")
    cn: IntProperty(name = "", default = 1, min = 1)
    rtsat: IntProperty(name = "", default = 0, min = 0)
    hours: IntProperty(name = "", default = 8760, min = 1, max = 8760)
    delim: EnumProperty(name = '', items = [("Comma", "Comma", "Comma delimiter"), ("Space", "Space", "space delimiter")], default = 'Comma')
    generate_file: StringProperty(default = "", name = "")
    (u1, u2, u3, u4) =  [StringProperty(name = "", description = "Valid entries (; separated for each 'For', comma separated for each day, space separated for each time value pair)", update = tupdate)] * 4
    (f1, f2, f3, f4) =  [StringProperty(name = "", description = "Valid entries (space separated): AllDays, Weekdays, Weekends, Monday, Tuesday, Wednesday, Thursday, Friday, Saturday, Sunday, AllOtherDays", update = tupdate)] * 4
    (t1, t2, t3, t4) = [IntProperty(name = "", default = 365, min = 1, max = 365, update = tupdate)] * 4

    def init(self, context):
        self.outputs.new('So_En_Sched', 'Schedule')
        self['scheddict'] = {'TSPSchedule': 'Any Number', 'VASchedule': 'Fraction', 'Fan Schedule': 'Fraction', 'HSchedule': 'Temperature', 'CSchedule': 'Temperature'}
        self.tupdate(context)
        nodecolour(self, 1)

    def draw_buttons(self, context, layout):
        uvals, u = (1, self.u1, self.u2, self.u3, self.u4), 0
        tvals = (0, self.t1, self.t2, self.t3, self.t4)
        newrow(layout, 'From file', self, 'file')
        
        if self.file == "1":
            newrow(layout, 'Select', self, 'select_file')
            newrow(layout, 'Columns', self, 'cn')
            newrow(layout, 'Skip rows', self, 'rtsat')
            newrow(layout, 'Delimiter', self, 'delim')
        elif self.file == "2":
            newrow(layout, 'Generate', self, 'generate_file')

        if self.file != "1":        
            while uvals[u] and tvals[u] < 365:
                [newrow(layout, v[0], self, v[1]) for v in (('End day {}:'.format(u+1), 't'+str(u+1)), ('Fors:', 'f'+str(u+1)), ('Untils:', 'u'+str(u+1)))]
                u += 1

    def update(self):
        for sock in self.outputs:
            socklink(sock, self.id_data.name)
        self.id_data.interface_update(bpy.context)

    def epwrite(self, name, stype):
        schedtext, ths = '', []
        for tosock in [link.to_socket for link in self.outputs['Schedule'].links]:
            if not schedtext:
                for t in (self.t1, self.t2, self.t3, self.t4):
                    ths.append(t)
                    if t == 365:
                        break
#                ths = [self.t1, self.t2, self.t3, self.t4]
                fos = [fs for fs in (self.f1, self.f2, self.f3, self.f4) if fs]
                uns = [us for us in (self.u1, self.u2, self.u3, self.u4) if us]
                ts, fs, us = rettimes(ths, fos, uns)
                
#                if self.file == '0':
                schedtext = epschedwrite(name, stype, ts, fs, us)
        return schedtext
            
    def epwrite_sel_file(self, name):               
        params = ('Name', 'ScheduleType', 'Name of File', 'Column Number', 'Rows to Skip at Top', 'Number of Hours of Data', 'Column Separator')
        paramvs = (name, 'Any number', os.path.abspath(self.select_file), self.cn, self.rtsat, 8760, self.delim) 
        schedtext = epentry('Schedule:File', params, paramvs)
        '''    Schedule:File,
        elecTDVfromCZ01res, !- Name
        Any Number, !- ScheduleType
        TDV_kBtu_CTZ01.csv, !- Name of File
        2, !- Column Number
        4, !- Rows to Skip at Top
        8760, !- Number of Hours of Data
        Comma; !- Column Separator'''
        return schedtext
    
    def epwrite_gen_file(self, name, data, newdir):
        schedtext, ths = '', []
        for tosock in [link.to_socket for link in self.outputs['Schedule'].links]:
            if not schedtext:
                for t in (self.t1, self.t2, self.t3, self.t4):
                    ths.append(t)
                    if t == 365:
                        break
#                ths = [self.t1, self.t2, self.t3, self.t4]
                fos = [fs for fs in (self.f1, self.f2, self.f3, self.f4) if fs]
                uns = [us for us in (self.u1, self.u2, self.u3, self.u4) if us]
                ts, fs, us = rettimes(ths, fos, uns)
        for t in ts:
            for f in fs:
                for u in us:
                    for hi, h in enumerate((datetime.datetime(2015, 1, 1, 0, 00) - datetime.datetime(2014, 1, 1, 0, 00)).hours):
                        if h.day <= self.ts:
#                            if f == 'Weekday'
                            data[hi] = 1
        with open(os.path.join(newdir, name), 'w') as sched_file:
            sched_file.write(',\n'.join([d for d in data]))
        params = ('Name', 'ScheduleType', 'Name of File', 'Column Number', 'Rows to Skip at Top', 'Number of Hours of Data', 'Column Separator')
        paramvs = (name, 'Any number', os.path.abspath(self.select_file), self.cn, self.rtsat, 8760, self.delim) 
        schedtext = epentry('Schedule:File', params, paramvs)    
        return schedtext
   
envi_mat_con = [NodeItem("No_En_Mat_Con", label="Construction Node")]
envi_mat_lay = [NodeItem("No_En_Mat_Op", label="Opaque layer"),
                NodeItem("No_En_Mat_Tr", label="Transparency layer"),
                NodeItem("No_En_Mat_Gas", label="Gas layer")]
envi_mat_sha = [NodeItem("No_En_Mat_Sh", label="Shading layer"),
                NodeItem("No_En_Mat_Bl", label="Blind layer"),
                NodeItem("No_En_Mat_Sc", label="Screen layer"),
                NodeItem("No_En_Mat_SG", label="Switchable layer"),
                NodeItem("No_En_Mat_ShC", label="Shading Control Node")]
envi_mat_sch = [NodeItem("No_En_Mat_Sched", label="Schedule")]
envi_mat_pv = [NodeItem("No_En_Mat_PV", label="PV"),
               NodeItem("No_En_Mat_PVG", label="PV Generator")] 

envimatnode_categories = [
        EnViMatNodeCategory("Type", "Type Node", items=envi_mat_con),
        EnViMatNodeCategory("Layer", "Layer Node", items=envi_mat_lay), 
        EnViMatNodeCategory("Shading", "Shading Node", items=envi_mat_sha),
        EnViMatNodeCategory("Schedule", "Schedule Node", items=envi_mat_sch),
        EnViMatNodeCategory("Power", "PV Node", items=envi_mat_pv)]
            
            

#    eff: FloatProperty(name = "%", description = "Efficiency (%)", min = 0.0, max = 100, default = 20)
#    ssp: IntProperty(name = "", description = "Number of series strings in parallel", min = 1, max = 100, default = 5)
#    mis: IntProperty(name = "", description = "Number of modules in series", min = 1, max = 100, default = 5)
#    pv: BoolProperty(name = "", description = "Photovoltaic shader", default = False, update = con_update)
#    pp: EnumProperty(items = [("0", "Simple", "Simple PV calculation"), 
#                               ("1", "One-Diode", "One-diode PV calculation"),
#                               ("2", "Sandia", "MSandia PV database")], 
#                                name = "", description = "Photovoltaic Performance Object Type", default = "0", update = con_update)
#    ct: EnumProperty(items = [("0", "Crystalline", "Do not model reflected beam component"), 
#                               ("1", "Amorphous", "Model reflectred beam as beam")], 
#                                name = "", description = "Photovoltaic Type", default = "0")
#    hti: EnumProperty(items = [("Decoupled", "Decoupled", "Decoupled"), 
#                               ("DecoupledUllebergDynamic", "Ulleberg", "DecoupledUllebergDynamic"),
#                               ("IntegratedSurfaceOutsideFace", "SurfaceOutside", "IntegratedSurfaceOutsideFace"),
#                               ("IntegratedTranspiredCollector", "Transpired", "IntegratedTranspiredCollector"),
#                               ("IntegratedExteriorVentedCavity", "ExteriorVented", "IntegratedExteriorVentedCavity"),
#                               ("PhotovoltaicThermalSolarCollector", "PVThermal", "PhotovoltaicThermalSolarCollector")], 
#                                name = "", description = "Conversion Efficiency Input Mode'", default = "Decoupled")

#    e1ddict = {'ASE 300-DFG/50': (6.2, 60, 50.5, 5.6, 0.001, -0.0038, 216, 318, 2.43), 
#               'BPsolar 275': (4.75,  21.4, 17, 4.45, 0.00065, -0.08, 36, 320, 0.63),
#               'BPsolar 3160': (4.8, 44.2, 35.1, 4.55, 0.00065, -0.16, 72, 320, 1.26),
#               'BPsolar 380': (4.8, 22.1, 17.6, 4.55, 0.00065, -0.08, 36, 320, 0.65),
#               'BPsolar 4160': (4.9, 44.2, 35.4, 4.52, 0.00065, -0.16, 72, 320, 1.26),
#               'BPsolar 5170': (5, 44.2, 36, 4.72, 0.00065, -0.16, 72, 320, 1.26),
#               'BPsolar 585': (5, 22.1, 18, 4.72, 0.00065, -0.08, 36, 320, 0.65),
#               'Shell SM110-12': (6.9, 21.7, 17.5, 6.3, 0.0028, -0.076, 36, 318, 0.86856),
#               'Shell SM110-24': (3.45, 43.5, 35, 3.15, 0.0014, -0.152, 72, 318, 0.86856),
#               'Shell SP70': (4.7, 21.4, 16.5, 4.25, 0.002, -0.076, 36, 318, 0.6324),
#               'Shell SP75': (4.8, 21.7, 17, 4.4, 0.002, -0.076, 36, 318, 0.6324),
#               'Shell SP140': (4.7, 42.8, 33, 4.25, 0.002, -0.152, 72, 318, 1.320308),
#               'Shell SP150': (4.8, 43.4, 34, 4.4, 0.002, -0.152, 72, 318, 1.320308),
#               'Shell S70': (4.5, 21.2, 17, 4, 0.002, -0.076, 36, 317, 0.7076),
#               'Shell S75': (4.7, 21.6, 17.6, 4.2, 0.002, -0.076, 36, 317, 0.7076),
#               'Shell S105': (4.5, 31.8, 25.5, 3.9, 0.002, -0.115, 54, 317, 1.037),
#               'Shell S115': (4.7, 32.8, 26.8, 4.2, 0.002, -0.115, 54, 317, 1.037),
#               'Shell ST40': (2.68, 23.3, 16.6, 2.41, 0.00035, -0.1, 16, 320, 0.424104),
#               'UniSolar PVL-64': (4.8, 23.8, 16.5, 3.88, 0.00065, -0.1, 40, 323, 0.65),
#               'UniSolar PVL-128': (4.8, 47.6, 33, 3.88, 0.00065, -0.2, 80, 323, 1.25),
#               'Custom': (None, None, None, None, None, None, None, None, None)}
#    
#    e1items = [(p, p, '{} module'.format(p)) for p in e1ddict]
#    e1menu: EnumProperty(items = e1items, name = "", description = "Module type", default = 'ASE 300-DFG/50')
#    pvsa: FloatProperty(name = "%", description = "Fraction of Surface Area with Active Solar Cells", min = 10, max = 100, default = 90)
#    aa: FloatProperty(name = "m2", description = "Active area", min = 0.1, max = 10000, default = 5)
#    eff: FloatProperty(name = "%", description = "Visible reflectance", min = 0.0, max = 100, default = 20)
#    ssp: IntProperty(name = "", description = "Number of series strings in parallel", min = 1, max = 100, default = 5)
#    mis: IntProperty(name = "", description = "Number of modules in series", min = 1, max = 100, default = 5)
#    cis: IntProperty(name = "", description = "Number of cells in series", min = 1, max = 100, default = 36) 
#    tap: FloatProperty(name = "", description = "Transmittance absorptance product", min = -1, max = 1, default = 0.9)
#    sbg: FloatProperty(name = "eV", description = "Semiconductor band-gap", min = 0.1, max = 5, default = 1.12)
#    sr: FloatProperty(name = "W", description = "Shunt resistance", min = 1, default = 1000000)
#    scc: FloatProperty(name = "Amps", description = "Short circuit current", min = 1, max = 1000, default = 25)
#    sgd: FloatProperty(name = "mm", description = "Screen to glass distance", min = 1, max = 1000, default = 25)
#    ocv: FloatProperty(name = "V", description = "Open circuit voltage", min = 0.0, max = 100, default = 60)
#    rt: FloatProperty(name = "C", description = "Reference temperature", min = 0, max = 40, default = 25)
#    ri: FloatProperty(name = "W/m2", description = "Reference insolation", min = 100, max = 2000, default = 1000)
#    mc: FloatProperty(name = "", description = "Module current at maximum power", min = 1, max = 10, default = 5.6)
#    mv: FloatProperty(name = "", description = "Module voltage at maximum power", min = 0.0, max = 75, default = 17)
#    tcscc: FloatProperty(name = "", description = "Temperature Coefficient of Short Circuit Current", min = 0.00001, max = 0.01, default = 0.002)
#    tcocv: FloatProperty(name = "", description = "Temperature Coefficient of Open Circuit Voltage", min = -0.5, max = 0, default = -0.1)
#    atnoct: FloatProperty(name = "C", description = "Reference ambient temperature", min = 0, max = 40, default = 20)
#    ctnoct: FloatProperty(name = "C", description = "Nominal Operating Cell Temperature Test Cell Temperature", min = 0, max = 60, default = 45)
#    inoct: FloatProperty(name = "W/m2", description = "Nominal Operating Cell Temperature Test Insolation", min = 100, max = 2000, default = 800)
#    hlc: FloatProperty(name = "W/m2.K", description = "Module heat loss coefficient", min = 0.0, max = 50, default = 30)
#    thc: FloatProperty(name = " J/m2.K", description = " Total Heat Capacity", min = 10000, max = 100000, default = 50000)
            
#            if self.envi_con_type in ("Wall", "Floor", "Roof", "Window", "Door"):
                
                
#                    newrow(layout, 'PV:', self, "pv")
#                    self.inputs['PV'].hide = False
                    
                    
#                    if self.pv:
#                        newrow(layout, "Heat transfer:", self, "hti")
#                        newrow(layout, "Photovoltaic:", self, "pp")
#                                                    
#                        if self.pp == '0':
#                            newrow(layout, "PV area ratio:", self, "pvsa")
#                            newrow(layout, "Efficiency:", self, "eff")
#                            
#                        elif self.pp == '1':
#                            newrow(layout, "Model:", self, "e1menu")
#                            newrow(layout, "Series in parallel:", self, "ssp")
#                            newrow(layout, "Modules in series:", self, "mis")
#                            
#                            if self.e1menu == 'Custom':
#                                newrow(layout, "Cell type:", self, "ct")
#                                newrow(layout, "Silicon:", self, "mis")
#                                newrow(layout, "Area:", self, "pvsa")
#                                newrow(layout, "Trans*absorp:", self, "tap")
#                                newrow(layout, "Band gap:", self, "sbg")
#                                newrow(layout, "Shunt:", self, "sr")    
#                                newrow(layout, "Short:", self, "scc") 
#                                newrow(layout, "Open:", self, "ocv")
#                                newrow(layout, "Ref. temp.:", self, "rt")
#                                newrow(layout, "Ref. insol.:", self, "ri")
#                                newrow(layout, "Max current:", self, "mc")
#                                newrow(layout, "Max voltage:", self, "mv")
#                                newrow(layout, "Max current:", self, "tcscc")
#                                newrow(layout, "Max voltage:", self, "tcocv")
#                                newrow(layout, "Ambient temp.:", self, "atnoct")
#                                newrow(layout, "Cell temp.:", self, "ctnoct")
#                                newrow(layout, "Insolation:", self, "inoct")
#                            else:
#                                newrow(layout, "Trans*absorp:", self, "tap")
#                                newrow(layout, "Band gap:", self, "sbg")
#                                newrow(layout, "Shunt:", self, "sr")    
#                                newrow(layout, "Ref. temp.:", self, "rt")
#                                newrow(layout, "Ref. insol.:", self, "ri")
#                                newrow(layout, "Test ambient:", self, "atnoct")
#                                newrow(layout, "Test Insolation:", self, "inoct")
#                                newrow(layout, "Heat loss coeff.:", self, "hlc")
#                                newrow(layout, "Heat capacity:", self, "thc")