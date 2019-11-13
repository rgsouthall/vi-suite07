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


import bpy, glob, os, inspect, datetime, shutil, time, math, mathutils, sys, shlex
from bpy.props import EnumProperty, FloatProperty, IntProperty, BoolProperty, StringProperty, FloatVectorProperty
from bpy.types import NodeTree, Node, NodeSocket
from nodeitems_utils import NodeCategory, NodeItem
from subprocess import Popen, PIPE
from .vi_func import socklink, socklink2, uvsocklink, newrow, epwlatilongi, nodeinputs, remlink, rettimes, sockhide, selobj, cbdmhdr, cbdmmtx
from .vi_func import hdrsky, nodecolour, facearea, retelaarea, iprop, bprop, eprop, fprop, sunposlivi, retdates, validradparams, retpmap
from .vi_func import delobj, logentry
from .envi_func import retrmenus, resnameunits, enresprops, epentry, epschedwrite, processf, get_mat, get_con_node
from .livi_export import livi_sun, livi_sky, livi_ground, hdrexport
from .envi_mat import envi_materials, envi_constructions, envi_layer, envi_layertype, envi_con_list


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
        self.startframe = self.startframe if self.animated and self.contextmenu == 'Basic' else scene.frame_current 
        self['endframe'] = self.startframe + int(((24 * (self.edoy - self.sdoy) + self.ehour - self.shour)/self.interval)) if self.contextmenu == 'Basic' and self.animated else scene.frame_current
        self['mtxfile'] = ''
        self['preview'] = 0
        
        if self.contextmenu == "Basic":  
            self['preview'] = 1
            
            if self.skyprog in ('0', '1'):
                self['skytypeparams'] = ("+s", "+i", "-c", "-b 22.86 -c")[self['skynum']] if self.skyprog == '0' else "-P {} {} -O {}".format(self.epsilon, self.delta, int(self.spectrummenu))
                print(self['skytypeparams'])
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
                    shutil.copyfile(self.skyname, "{}-0.sky".format(scene['viparams']['filebase']))
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
                    with open("{}.mtx".format(os.path.join(scene['viparams']['newdir'], self['epwbase'][0])), 'r') as mtxfile:
                        self['Options']['MTX'] = mtxfile.read()
                else:
                    with open(self.mtxname, 'r') as mtxfile:
                        self['Options']['MTX'] = mtxfile.read()
                if self.hdr:
                    self['Text'][str(scene.frame_current)] = cbdmhdr(self, scene)

        elif self.contextmenu == "Compliance":
            if self.canalysismenu in ('0', '1', '2'):            
                self['skytypeparams'] = ("-b 22.86 -c", "-b 22.86 -c", "-b 18 -u")[int(self.canalysismenu)]
                skyentry = livi_sun(scene, self, 0, 0) + livi_sky(3)
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
                    with open("{}.mtx".format(os.path.join(scene['viparams']['newdir'], self['epwbase'][0])), 'r') as mtxfile:
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
        self["_RNA_UI"] = {"Processors": {"min": 1, "max": int(context.scene['viparams']['nproc']), "name": ""}}
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
            newrow(layout, 'Multi-thread:', self, 'mp')
            
            if self.mp and sys.platform != 'win32':
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
        self.outputs.new('ViR', 'Results out')
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
        self.inputs.new('ViLoc', 'Location in')
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
        self.outputs.new('ViR', 'Results out')
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
    
#    def nodeupdate(self, context):
#        nodecolour(self, self['exportstate'] != [str(x) for x in (self.animmenu)])
    
    def init(self, context):
        self.outputs.new('So_En_Geo', 'Geometry out')
#        self['exportstate'] = ''
        nodecolour(self, 1)

    def draw_buttons(self, context, layout):
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
        scene['enparams']['fs'], scene['enparams']['fe'] = self.fs, self.fe
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
        self['frames'] = range(context.scene['enparams']['fs'], context.scene['enparams']['fe'] + 1)
        self.resfilename = os.path.join(context.scene['viparams']['newdir'], self.resname+'.eso')
        self['year'] = innode['year']
        self.dsdoy = innode.sdoy # (locnode.startmonthnode.sdoy
        self.dedoy = innode.edoy
        self["_RNA_UI"] = {"Start": {"min":innode.sdoy, "max":innode.edoy}, "End": {"min":innode.sdoy, "max":innode.edoy}, "AStart": {"name": '', "min":context.scene['enparams']['fs'], "max":context.scene['enparams']['fe']}, "AEnd": {"min":context.scene['enparams']['fs'], "max":context.scene['enparams']['fe']}}
        self['Start'], self['End'] = innode.sdoy, innode.edoy
#        self["_RNA_UI"] = {"AStart": {"min":context.scene['enparams']['fs'], "max":context.scene['enparams']['fe']}, "AEnd": {"min":context.scene['enparams']['fs'], "max":context.scene['enparams']['fe']}}
        self['AStart'], self['AEnd'] = context.scene['enparams']['fs'], context.scene['enparams']['fe']
     
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

vigennodecat = []

vi_display = [NodeItem("No_Vi_Chart", label="Chart")]
vioutnodecat = []
vi_image = [NodeItem("No_Li_Im", label="LiVi Image"), NodeItem("No_Li_Gl", label="LiVi Glare"), NodeItem("No_Li_Fc", label="LiVi False-colour")]
vi_input = [NodeItem("No_Loc", label="VI Location")]

vinode_categories = [ViNodeCategory("Output", "Output Nodes", items=vioutnodecat), 
                     ViNodeCategory("Edit", "Edit Nodes", items=vi_edit), 
                     ViNodeCategory("Image", "Image Nodes", items=vi_image), 
                     ViNodeCategory("Display", "Display Nodes", items=vi_display), 
                     ViNodeCategory("Generative", "Generative Nodes", items=vigennodecat), 
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
        layout.label(text)

    def draw_color(self, context, node):
        return (0.5, 0.2, 0.0, 0.75)

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
        layout.label(text)

    def draw_color(self, context, node):
        return (0.2, 0.2, 0.2, 0.75)

class So_En_Net_Act(NodeSocket):
    '''An EnVi actuator socket'''
    bl_idname = 'So_En_Net_Act'
    bl_label = 'EnVi actuator socket'

    sn: StringProperty()
    valid = ['Actuator']

    def draw(self, context, layout, node, text):
        layout.label(text)

    def draw_color(self, context, node):
        return (0.2, 0.9, 0.9, 0.75)

class So_En_Net_Sense(NodeSocket):
    '''An EnVi sensor socket'''
    bl_idname = 'So_En_Net_Sense'
    bl_label = 'EnVi sensor socket'

    sn: StringProperty()
    valid = ['Sensor']

    def draw(self, context, layout, node, text):
        layout.label(text)

    def draw_color(self, context, node):
        return (0.9, 0.9, 0.2, 0.75)
    
class No_En_Net_Zone(Node, EnViNodes):
    '''Node describing a simulation zone'''
    bl_idname = 'No_En_Net_Zone'
    bl_label = 'Zone'
    bl_icon = 'SOUND'

    def zupdate(self, context):
        self.afs = 0
        col = bpy.data.collections[self.zone]
        
        for obj in col.objects:
            odm = obj.data.materials
            bfacelist = sorted([face for face in obj.data.polygons if get_con_node(odm[face.material_index].vi_params).envi_con_con == 'Zone'], key=lambda face: -face.center[2])
    #        buvals = [retuval(odm[face.material_index]) for face in bfacelist]
            bsocklist = ['{}_{}_b'.format(odm[face.material_index].name, face.index) for face in bfacelist]
            sfacelist = sorted([face for face in obj.data.polygons if get_con_node(odm[face.material_index].vi_params).envi_afsurface == 1 and get_con_node(odm[face.material_index]).envi_con_type not in ('Window', 'Door')], key=lambda face: -face.center[2])
            ssocklist = ['{}_{}_s'.format(odm[face.material_index].name, face.index) for face in sfacelist]
            ssfacelist = sorted([face for face in obj.data.polygons if get_con_node(odm[face.material_index].vi_params).envi_afsurface == 1 and get_con_node(odm[face.material_index]).envi_con_type in ('Window', 'Door')], key=lambda face: -face.center[2])
            sssocklist = ['{}_{}_ss'.format(odm[face.material_index].name, face.index) for face in ssfacelist]
    
            [self.outputs.remove(oname) for oname in self.outputs if oname.bl_idname in ('So_En_Net_Bound', 'So_En_Net_SFlow', 'So_En_Net_SSFlow')]
            [self.inputs.remove(iname) for iname in self.inputs if iname.bl_idname in ('So_En_Net_Bound', 'So_En_Net_SFlow', 'So_En_Net_SSFlow')]
    
            for sock in bsocklist:
                self.outputs.new('So_En_Net_Bound', sock).sn = sock.split('_')[-2]
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
        self.vol_update(context)
        
    def vol_update(self, context):
        coll = bpy.data.collections[self.zone] 
        for obj in coll.objects:
            obj['volume'] = obj['auto_volume'] if self.volcalc == '0' else self.zonevolume
            self['volume'] = obj['auto_volume']
            
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
            print('hello')
        except Exception as e:
            print("Don't panic")
            
        self.alllinked = 1 if all((bi, si, ssi, bo, so, sso)) else 0
        nodecolour(self, self.errorcode())
        
    def uvsockupdate(self):
        for sock in self.outputs:
            socklink(sock, self['nodeid'].split('@')[1])
            if sock.bl_idname == 'So_En_Net_Bound':
                uvsocklink(sock, self['nodeid'].split('@')[1])
    
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
    
class EnViNodeCategory(NodeCategory):
    @classmethod
    def poll(cls, context):
        return context.space_data.tree_type == 'EnViN'

envi_zone = [NodeItem("No_En_Net_Zone", label="EnVi Zone"), NodeItem("No_En_Net_Occ", label="EnVi Occupancy")]
envinode_categories = [EnViNodeCategory("Nodes", "Zone Nodes", items=envi_zone)]

        
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
        layout.label(text)

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
        layout.label(text)

    def draw_color(self, context, node):
        return (1, 0.1, 1, 1.0)
    
    def ret_valid(self, node):
        return ['OLayer']
    
class So_En_Mat_G(NodeSocket):
    '''EnVi gas layer socket'''
    bl_idname = 'So_En_Mat_G'
    bl_label = 'Gas layer socket'
    valid = ['GLayer']

    def draw(self, context, layout, node, text):
        layout.label(text)

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
    
class No_En_Mat_Con(Node, EnViMatNodes):
    '''Node defining the EnVi material construction'''
    bl_idname = 'No_En_Mat_Con'
    bl_label = 'EnVi Construction'
    bl_icon = 'FORCE_WIND'
    
    def con_update(self, context):
        if len(self.inputs) == 2:
            if not self.pv:
                remlink(self, self.inputs['PV Schedule'].links)
                self.inputs['PV Schedule'].hide = True
                if self.envi_con_makeup != "1" or self.envi_con_type in ('Shading', 'None'):
                    for link in self.inputs['Outer layer'].links:
                        self.id_data.links.remove(link)
                    self.inputs['Outer layer'].hide = True
                else:
                    self.inputs['Outer layer'].hide = False
            else:
                self.inputs['Outer layer'].hide = False
                if self.pp != '0':
                    remlink(self, self.inputs['PV Schedule'].links)
                    self.inputs['PV Schedule'].hide = True
                else:
                    self.inputs['PV Schedule'].hide = False
                
            [link.from_node.update() for link in self.inputs['Outer layer'].links]
            get_mat(self, 0).vi_params.envi_type = self.envi_con_type        
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
            return [("", "", "")]
        
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
                                      
    matname: StringProperty(name = "", description = "", default = '')
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

    eff: FloatProperty(name = "%", description = "Efficiency (%)", min = 0.0, max = 100, default = 20)
    ssp: IntProperty(name = "", description = "Number of series strings in parallel", min = 1, max = 100, default = 5)
    mis: IntProperty(name = "", description = "Number of modules in series", min = 1, max = 100, default = 5)
    pv: BoolProperty(name = "", description = "Photovoltaic shader", default = False, update = con_update)
    pp: EnumProperty(items = [("0", "Simple", "Do not model reflected beam component"), 
                               ("1", "One-Diode", "Model reflectred beam as beam"),
                               ("2", "Sandia", "Model reflected beam as diffuse")], 
                                name = "", description = "Photovoltaic Performance Object Type", default = "0", update = con_update)
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

    e1ddict = {'ASE 300-DFG/50': (6.2, 60, 50.5, 5.6, 0.001, -0.0038, 216, 318, 2.43), 
               'BPsolar 275': (4.75,  21.4, 17, 4.45, 0.00065, -0.08, 36, 320, 0.63),
               'BPsolar 3160': (4.8, 44.2, 35.1, 4.55, 0.00065, -0.16, 72, 320, 1.26),
               'BPsolar 380': (4.8, 22.1, 17.6, 4.55, 0.00065, -0.08, 36, 320, 0.65),
               'BPsolar 4160': (4.9, 44.2, 35.4, 4.52, 0.00065, -0.16, 72, 320, 1.26),
               'BPsolar 5170': (5, 44.2, 36, 4.72, 0.00065, -0.16, 72, 320, 1.26),
               'BPsolar 585': (5, 22.1, 18, 4.72, 0.00065, -0.08, 36, 320, 0.65),
               'Shell SM110-12': (6.9, 21.7, 17.5, 6.3, 0.0028, -0.076, 36, 318, 0.86856),
               'Shell SM110-24': (3.45, 43.5, 35, 3.15, 0.0014, -0.152, 72, 318, 0.86856),
               'Shell SP70': (4.7, 21.4, 16.5, 4.25, 0.002, -0.076, 36, 318, 0.6324),
               'Shell SP75': (4.8, 21.7, 17, 4.4, 0.002, -0.076, 36, 318, 0.6324),
               'Shell SP140': (4.7, 42.8, 33, 4.25, 0.002, -0.152, 72, 318, 1.320308),
               'Shell SP150': (4.8, 43.4, 34, 4.4, 0.002, -0.152, 72, 318, 1.320308),
               'Shell S70': (4.5, 21.2, 17, 4, 0.002, -0.076, 36, 317, 0.7076),
               'Shell S75': (4.7, 21.6, 17.6, 4.2, 0.002, -0.076, 36, 317, 0.7076),
               'Shell S105': (4.5, 31.8, 25.5, 3.9, 0.002, -0.115, 54, 317, 1.037),
               'Shell S115': (4.7, 32.8, 26.8, 4.2, 0.002, -0.115, 54, 317, 1.037),
               'Shell ST40': (2.68, 23.3, 16.6, 2.41, 0.00035, -0.1, 16, 320, 0.424104),
               'UniSolar PVL-64': (4.8, 23.8, 16.5, 3.88, 0.00065, -0.1, 40, 323, 0.65),
               'UniSolar PVL-128': (4.8, 47.6, 33, 3.88, 0.00065, -0.2, 80, 323, 1.25),
               'Custom': (None, None, None, None, None, None, None, None, None)}
    
    e1items = [(p, p, '{} module'.format(p)) for p in e1ddict]
    e1menu: EnumProperty(items = e1items, name = "", description = "Module type", default = 'ASE 300-DFG/50')
    pvsa: FloatProperty(name = "%", description = "Fraction of Surface Area with Active Solar Cells", min = 10, max = 100, default = 90)
    aa: FloatProperty(name = "m2", description = "Active area", min = 0.1, max = 10000, default = 5)
    eff: FloatProperty(name = "%", description = "Visible reflectance", min = 0.0, max = 100, default = 20)
    ssp: IntProperty(name = "", description = "Number of series strings in parallel", min = 1, max = 100, default = 5)
    mis: IntProperty(name = "", description = "Number of modules in series", min = 1, max = 100, default = 5)
    cis: IntProperty(name = "", description = "Number of cells in series", min = 1, max = 100, default = 36) 
    tap: FloatProperty(name = "", description = "Transmittance absorptance product", min = -1, max = 1, default = 0.9)
    sbg: FloatProperty(name = "eV", description = "Semiconductor band-gap", min = 0.1, max = 5, default = 1.12)
    sr: FloatProperty(name = "W", description = "Shunt resistance", min = 1, default = 1000000)
    scc: FloatProperty(name = "Amps", description = "Short circuit current", min = 1, max = 1000, default = 25)
    sgd: FloatProperty(name = "mm", description = "Screen to glass distance", min = 1, max = 1000, default = 25)
    ocv: FloatProperty(name = "V", description = "Open circuit voltage", min = 0.0, max = 100, default = 60)
    rt: FloatProperty(name = "C", description = "Reference temperature", min = 0, max = 40, default = 25)
    ri: FloatProperty(name = "W/m2", description = "Reference insolation", min = 100, max = 2000, default = 1000)
    mc: FloatProperty(name = "", description = "Module current at maximum power", min = 1, max = 10, default = 5.6)
    mv: FloatProperty(name = "", description = "Module voltage at maximum power", min = 0.0, max = 75, default = 17)
    tcscc: FloatProperty(name = "", description = "Temperature Coefficient of Short Circuit Current", min = 0.00001, max = 0.01, default = 0.002)
    tcocv: FloatProperty(name = "", description = "Temperature Coefficient of Open Circuit Voltage", min = -0.5, max = 0, default = -0.1)
    atnoct: FloatProperty(name = "C", description = "Reference ambient temperature", min = 0, max = 40, default = 20)
    ctnoct: FloatProperty(name = "C", description = "Nominal Operating Cell Temperature Test Cell Temperature", min = 0, max = 60, default = 45)
    inoct: FloatProperty(name = "W/m2", description = "Nominal Operating Cell Temperature Test Insolation", min = 100, max = 2000, default = 800)
    hlc: FloatProperty(name = "W/m2.K", description = "Module heat loss coefficient", min = 0.0, max = 50, default = 30)
    thc: FloatProperty(name = " J/m2.K", description = " Total Heat Capacity", min = 10000, max = 100000, default = 50000)
    
    def init(self, context):
        self.inputs.new('So_En_Sched', 'PV Schedule')
        self.inputs['PV Schedule'].hide = True
        self.inputs.new('So_En_Mat_Ou', 'Outer layer')
        self.inputs['Outer layer'].hide = True
#        self.inputs.new('SO_EN_Mat_Fr', 'Outer frame layer')
#        self.inputs['Outer frame layer'].hide = True
        
    def draw_buttons(self, context, layout):
        newrow(layout, 'Active:', self, 'active')
        newrow(layout, 'Type:', self, "envi_con_type")
       
        if self.envi_con_type != "None":
            newrow(layout, 'Boundary:', self, "envi_con_con")

            if self.envi_con_type in ("Wall", "Floor", "Roof", "Window", "Door"):
                
                if self.envi_con_type in ("Wall", "Roof") and self.envi_con_con == 'External':
                    newrow(layout, 'PV:', self, "pv")
                    
                    if self.pv:
                        newrow(layout, "Heat transfer:", self, "hti")
                        newrow(layout, "Photovoltaic:", self, "pp")
                                                    
                        if self.pp == '0':
                            newrow(layout, "PV area ratio:", self, "pvsa")
                            newrow(layout, "Efficiency:", self, "eff")
                            
                    elif self.pp == '1':
                        newrow(layout, "Model:", self, "e1menu")
                        newrow(layout, "Series in parallel:", self, "ssp")
                        newrow(layout, "Modules in series:", self, "mis")
                        
                        if self.e1menu == 'Custom':
                            newrow(layout, "Cell type:", self, "ct")
                            newrow(layout, "Silicon:", self, "mis")
                            newrow(layout, "Area:", self, "pvsa")
                            newrow(layout, "Trans*absorp:", self, "tap")
                            newrow(layout, "Band gap:", self, "sbg")
                            newrow(layout, "Shunt:", self, "sr")    
                            newrow(layout, "Short:", self, "scc") 
                            newrow(layout, "Open:", self, "ocv")
                            newrow(layout, "Ref. temp.:", self, "rt")
                            newrow(layout, "Ref. insol.:", self, "ri")
                            newrow(layout, "Max current:", self, "mc")
                            newrow(layout, "Max voltage:", self, "mv")
                            newrow(layout, "Max current:", self, "tcscc")
                            newrow(layout, "Max voltage:", self, "tcocv")
                            newrow(layout, "Ambient temp.:", self, "atnoct")
                            newrow(layout, "Cell temp.:", self, "ctnoct")
                            newrow(layout, "Insolation:", self, "inoct")
                        else:
                            newrow(layout, "Trans*absorp:", self, "tap")
                            newrow(layout, "Band gap:", self, "sbg")
                            newrow(layout, "Shunt:", self, "sr")    
                            newrow(layout, "Ref. temp.:", self, "rt")
                            newrow(layout, "Ref. insol.:", self, "ri")
                            newrow(layout, "Test ambient:", self, "atnoct")
                            newrow(layout, "Test Insolation:", self, "inoct")
                            newrow(layout, "Heat loss coeff.:", self, "hlc")
                            newrow(layout, "Heat capacity:", self, "thc")
                            
            if self.envi_con_type != "Shading" and not self.pv:
                if self.envi_con_con in ('External', 'Zone'):
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
                   
            if self.envi_con_type == 'Window':
                newrow(layout, 'Frame:', self, "fclass")
                if self.fclass == '0':
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
            
            elif self.envi_con_type in ('Wall', 'Floor', 'Roof') and not self.pv:
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
            elif self.envi_con_type == 'PV' and self.envi_con_makeup == '0':
                newrow(layout, "Series in parallel:", self, "ssp")
                newrow(layout, "Modules in series:", self, "mis")
                newrow(layout, "Area:", self, "fsa")
                newrow(layout, "Efficiency:", self, "eff")
        
    def update(self):
        if len(self.inputs) == 3:
            self.valid()
    
    def valid(self):
        if ((self.envi_con_makeup == '1' or self.pv) and not self.inputs['Outer layer'].links and self.envi_con_type != 'Shading') or \
            (not self.inputs['Outer frame layer'].links and not self.inputs['Outer frame layer'].hide):
            nodecolour(self, 1)
        else:
            nodecolour(self, 0)
            
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
     
    def ep_write(self):
        self['matname'] = get_mat(self, 1).name
        con_type = {'Roof': 'Ceiling', 'Floor': 'Internal floor', 'Wall': 'Internal wall'}[self.envi_con_type] if self.envi_con_con in ('Thermal mass', 'Zone') and self.envi_con_type in ('Roof', 'Wall', 'Floor') else self.envi_con_type
   
        if self.envi_con_makeup == '0':
            if self.envi_con_type == 'Window' and self.envi_simple_glazing:
                params = ['Name', 'Outside layer']
                paramvs = [self['matname'], self['matname'] + '_sg']
                ep_text = epentry('Construction', params, paramvs)
                params = ('Name', 'U-Factor', 'Solar Heat Gain Coefficient', 'Visible Transmittance')
                paramvs = [self['matname'] + '_sg'] + ['{:.3f}'.format(p) for p in (self.envi_sg_uv, self.envi_sg_shgc, self.envi_sg_vt)]                
                ep_text += epentry("WindowMaterial:SimpleGlazingSystem", params, paramvs)                          
            else:
                self.thicklist = [self.lt0, self.lt1, self.lt2, self.lt3, self.lt4, self.lt5, self.lt6, self.lt7, self.lt8, self.lt9]
                mats = envi_cons.propdict[con_type][self.envi_con_list]
                params = ['Name', 'Outside layer'] + ['Layer {}'.format(i + 1) for i in range(len(mats) - 1)]        
                paramvs = [self['matname']] + ['{}-layer-{}'.format(self['matname'], mi) for mi, m in enumerate(mats)]
                ep_text = epentry('Construction', params, paramvs)
                
                for pm, presetmat in enumerate(mats):  
                    matlist = list(envi_mats.matdat[presetmat])
                    layer_name = '{}-layer-{}'.format(self['matname'], pm)
                    
                    if envi_mats.namedict.get(presetmat) == None:
                        envi_mats.namedict[presetmat] = 0
                        envi_mats.thickdict[presetmat] = [self.thicklist[pm]/1000]
                    else:
                        envi_mats.namedict[presetmat] = envi_mats.namedict[presetmat] + 1
                        envi_mats.thickdict[presetmat].append(self.thicklist[pm]/1000)
                    
                    if self.envi_con_type in ('Wall', 'Floor', 'Roof', 'Ceiling', 'Door') and presetmat not in envi_mats.gas_dat:
                        self.resist += self.thicklist[pm]/1000/float(matlist[1])
                        params = ('Name', 'Roughness', 'Thickness (m)', 'Conductivity (W/m-K)', 'Density (kg/m3)', 'Specific Heat Capacity (J/kg-K)', 'Thermal Absorptance', 'Solar Absorptance', 'Visible Absorptance')                    
                        paramvs = ['{}-layer-{}'.format(self['matname'], pm), matlist[0], str(self.thicklist[pm]/1000)] + matlist[1:8]                    
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
                            pcmparamsv = ('{} CondFD override'.format(self['matname']), 'ConductionFiniteDifference', self['matname'])
                            ep_text += epentry('SurfaceProperty:HeatTransferAlgorithm:Construction', pcmparams, pcmparamsv)
    
                    elif presetmat in envi_mats.gas_dat:
                        params = ('Name', 'Resistance')
                        paramvs = ('{}-layer-{}'.format(self['matname'], pm), matlist[2])
                        ep_text += epentry("Material:AirGap", params, paramvs)
                    
                    elif self.envi_con_type =='Window':
                        if envi_mats.matdat[presetmat][0] == 'Glazing':
                            params = ('Name', 'Optical Data Type', 'Window Glass Spectral Data Set Name', 'Thickness (m)', 'Solar Transmittance at Normal Incidence', 'Front Side Solar Reflectance at Normal Incidence',
                          'Back Side Solar Reflectance at Normal Incidence', 'Visible Transmittance at Normal Incidence', 'Front Side Visible Reflectance at Normal Incidence', 'Back Side Visible Reflectance at Normal Incidence',
                          'Infrared Transmittance at Normal Incidence', 'Front Side Infrared Hemispherical Emissivity', 'Back Side Infrared Hemispherical Emissivity', 'Conductivity (W/m-K)',
                          'Dirt Correction Factor for Solar and Visible Transmittance', 'Solar Diffusing')
                            paramvs = ['{}-layer-{}'.format(self['matname'], pm)] + matlist[1:3] + [self.thicklist[pm]] + ['{:.3f}'.format(float(sm)) for sm in matlist[4:-1]] + [1, ('No', 'Yes')[matlist[-1]]]
                            ep_text += epentry("WindowMaterial:{}".format(matlist[0]), params, paramvs)
                    
                        elif envi_mats.matdat[presetmat][0] == 'Gas':
                            params = ('Name', 'Gas Type', 'Thickness')
                            paramvs = [layer_name] + [matlist[1]] + [self.thicklist[pm]]
                            ep_text += epentry("WindowMaterial:Gas", params, paramvs)
                    
        elif self.envi_con_makeup == '1':            
            in_sock = self.inputs['Outer layer']# if self.envi_con_type == "Window" else self.inputs[0]
            n = 0
            params = ['Name']
            paramvs = [self['matname']]
            ep_text = ''
            self.resist = 0
            get_mat(self, 1).vi_params.envi_shading = 0

            while in_sock.links:
                node = in_sock.links[0].from_node

                if node.bl_idname not in ('envi_sl_node', 'envi_bl_node', 'envi_screen_node', 'envi_sgl_node'):                    
                    paramvs.append('{}-layer-{}'.format(self['matname'], n)) 
                    params.append(('Outside layer', 'Layer {}'.format(n))[n > 0])
                    ep_text += node.ep_write(n)  
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
                paramvs = ['{}-shading'.format(self['matname'])]
                
                while in_sock.links:
                    node = in_sock.links[0].from_node
                    
                    if node.outputs['Layer'].links[0].to_node.bl_idname != 'envi_sgl_node':
                        paramvs.append('{}-layer-{}'.format(self['matname'], n)) 
                        params.append(('Outside layer', 'Layer {}'.format(n))[n > 0])
                    
                    in_sock = node.inputs['Layer']

                    if node.bl_idname in ('envi_sl_node', 'envi_bl_node', 'envi_screen_node', 'envi_sgl_node'):
                        ep_text += node.ep_write(n)
                    
                    n += 1
                ep_text += epentry('Construction', params, paramvs)
                
        if self.envi_con_type =='Window': 
            if self.fclass == '0':
                params = ('Name', 'Roughness', 'Thickness (m)', 'Conductivity (W/m-K)', 'Density (kg/m3)', 'Specific Heat (J/kg-K)', 'Thermal Absorptance', 'Solar Absorptance', 'Visible Absorptance', 'Name', 'Outside Layer')
                paramvs = ('{}-frame-layer{}'.format(self['matname'], 0), 'Rough', '0.12', '0.1', '1400.00', '1000', '0.9', '0.6', '0.6', '{}-frame'.format(self['matname']), '{}-frame-layer{}'.format(self['matname'], 0))
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
                ep_text += self.layer_write(self.inputs['Outer frame layer'], self['matname'])
        
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
    

    layer: EnumProperty(items = [("0", "Database", "Select from database"), 
                                        ("1", "Custom", "Define custom material properties")], 
                                        name = "", description = "Class of layer", default = "0")

    materialtype: EnumProperty(items = envi_layertype, name = "", description = "Layer material type")
    material: EnumProperty(items = envi_layer, name = "", description = "Layer material")
    thi: FloatProperty(name = "mm", description = "Thickness (mm)", min = 0.1, max = 10000, default = 100)
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
    shc: FloatProperty(name = "J/kg", description = "Thickness (mm)", min = 0.001, max = 10000, default = 800)
    tab: FloatProperty(name = "", description = "Thickness (mm)", min = 0, max = 1, default = 0.7)
    sab: FloatProperty(name = "", description = "Thickness (mm)", min = 0, max = 1, default = 0.7)
    vab: FloatProperty(name = "", description = "Thickness (mm)", min = 0, max = 1, default = 0.7)
    pcm: BoolProperty(name = "", description = "Phase Change Material", default = 0)
    tctc: FloatProperty(name = "", description = "Temp. coeff. for thermal conductivity (W/m-K2)", min = 0, max = 50, default = 0)
    tempemps: StringProperty(name = "", description = "Temperature/empalthy pairs", default = "")
    resist: FloatProperty(name = "", description = "", min = 0, default = 0) 
    envi_con_type: StringProperty(name = "", description = "Name")
    
    def init(self, context):
        self.outputs.new('So_En_Mat_Op', 'Layer')
        self.inputs.new('So_En_Mat_Op', 'Layer')
        
    def draw_buttons(self, context, layout):
        newrow(layout, "Class:", self, "layer")
        if self.layer == '0':
            newrow(layout, "Type:", self, "materialtype")
            newrow(layout, "Material:", self, "material")
            newrow(layout, "Thickness:", self, "thi") 
        else:
            newrow(layout, "Conductivity:", self, "tc")
            newrow(layout, "Thickness:", self, "thi")
            newrow(layout, "Roughness:", self, "rough")
            newrow(layout, "Density:", self, "rho")
            newrow(layout, "SHC:", self, "shc")
            newrow(layout, "Therm absorb:", self, "tab")
            newrow(layout, "Solar absorb:", self, "sab")
            newrow(layout, "Vis absorb:", self, "vab")
            newrow(layout, "PCM:", self, "pcm")

            if self.pcm:
                newrow(layout, "TCTC:", self, "tctc")
                newrow(layout, "Temps:Emps", self, "tempemps")
    
    def ret_resist(self):
        if self.layer == '0':
            matlist = list(envi_mats.matdat[self.material])
            
            if self.materialtype != '6': 
                self.resist = self.thi * 0.001/float(matlist[1])
            else:
                self.resist = float(matlist[2])
            
        else:
            self.resist = self.thi/self.tc
        
        return self.resist
            
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
        pcmparamsv = ('{} CondFD override'.format(self['layer_name']), 'ConductionFiniteDifference', self['layer_name'])
    
        return epentry("MaterialProperty:PhaseChange", params, paramvs) + epentry('SurfaceProperty:HeatTransferAlgorithm:Construction', pcmparams, pcmparamsv)
        
    def ep_write(self, ln):
        for material in bpy.data.materials:
            if self.id_data == material.vi_params.envi_nodes:
                break
        self['matname'] = get_mat(self, 1).name
        self['layer_name'] = '{}-layer-{}'.format(self['matname'], ln) if self.envi_con_type != 'Frame' else '{}-frame-layer-{}'.format(self['matname'], ln)
        
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
            paramvs = ['{}-layer-{}'.format(material.name, ln), self.rough, '{:.3f}'.format(self.thi * 0.001), '{:.3f}'.format(self.tc), '{:.3f}'.format(self.rho), '{:.3f}'.format(self.shc), '{:.3f}'.format(self.tab), 
                       '{:.3f}'.format(self.sab), '{:.3f}'.format(self.vab)]
            
        ep_text = epentry(header, params, paramvs)
        
        if self.pcm and self.material in envi_mats.pcmd_datd:
            ep_text += self.pcm_write()
        
        return ep_text
            
class No_En_Mat_Tr(Node, EnViMatNodes):
    '''Node defining the EnVi transparent material layer'''
    bl_idname = 'No_En_Mat_Tr'
    bl_label = 'EnVi transparent layer'
    
    layer: EnumProperty(items = [("0", "Database", "Select from database"), 
                                        ("1", "Custom", "Define custom material properties")], 
                                        name = "", description = "Composition of the layer", default = "0")
    materialtype: EnumProperty(items = envi_layertype, name = "", description = "Layer material type")
    material: EnumProperty(items = envi_layer, name = "", description = "Layer material")
    thi: FloatProperty(name = "mm", description = "Thickness (mm)", min = 0.1, max = 1000, default = 6)
    tc: FloatProperty(name = "W/m.K", description = "Thermal Conductivity (W/m.K)", min = 0.1, max = 10, default = 0.8)
    stn: FloatProperty(name = "", description = "Solar normal transmittance", min = 0, max = 1, default = 0.7)
    fsn: FloatProperty(name = "", description = "Solar front normal reflectance", min = 0, max = 1, default = 0.7)
    bsn: FloatProperty(name = "", description = "Solar back normal reflectance", min = 0, max = 1, default = 0.7)
    vtn: FloatProperty(name = "", description = "Visible Transmittance at Normal Incidence", min = 0, max = 1, default = 0.7)
    fvrn: FloatProperty(name = "", description = "Front Side Visible Reflectance at Normal Incidence", min = 0, max = 1, default = 0.1)
    bvrn: FloatProperty(name = "", description = "Back Side Visible Reflectance at Normal Incidence", min = 0, max = 1, default = 0.1)
    itn: FloatProperty(name = "", description = "Infrared Transmittance at Normal Incidence", min = 0, max = 1, default = 0.1)
    fie: FloatProperty(name = "", description = "Front Side Infrared Hemispherical Emissivity'", min = 0, max = 1, default = 0.7)
    bie: FloatProperty(name = "", description = "Back Side Infrared Hemispherical Emissivity", min = 0, max = 1, default = 0.7)
    diff: BoolProperty(name = "", description = "Diffusing", default = 0)
    envi_con_type: StringProperty(name = "", description = "Name")
    
    def init(self, context):
        self.outputs.new('envi_tl_sock', 'Layer')
        self.inputs.new('envi_gl_sock', 'Layer')
        
    def draw_buttons(self, context, layout):
        if self.outputs['Layer'].links:
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
        socklink2(self.outputs['Layer'], self.id_data)

        if self.outputs['Layer'].links:
            self.envi_con_type = self.outputs['Layer'].links[0].to_node.envi_con_type

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
    mw: FloatProperty(name = "kg/kmol", description = "Molecular weight", min = 0, max = 1, default = 0.7)
    shr: FloatProperty(name = "", description = "Specific heat ratio", min = 0, max = 1, default = 0.7)
    envi_con_type: StringProperty(name = "", description = "Name")
    
    def init(self, context):
        self.outputs.new('envi_gl_sock', 'Layer')
        self.inputs.new('envi_tl_sock', 'Layer')
        
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

    def ep_write(self, ln):
        for material in bpy.data.materials:
            if self.id_data == material.envi_nodes:
                break
        if self.layer == '0':
            params = ('Name', 'Gas Type', 'Thickness')
            paramvs = ['{}-layer-{}'.format(material.name, ln), self.material, self.thi]
            
        else:
            params = ('gap name', 'type', 'thickness', 'Conductivity Coefficient A', 'Conductivity Coefficient B', 'Conductivity Coefficient C', 
                      'Conductivity Viscosity A', 'Conductivity Viscosity B', 'Conductivity Viscosity C', 'Specific Heat Coefficient A',
                      'Specific Heat Coefficient B', 'Specific Heat Coefficient C', 'Molecular Weight', 'Specific Heat Ratio')
            paramvs = ['{}-layer-{}'.format(material.name, ln), 'Custom', '{:.3f}'.format(self.thi), '{:.3f}'.format(self.ccA), '{:.3f}'.format(self.ccB), 
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
        
    def ep_write(self, ln):
        for material in bpy.data.materials:
            if self.id_data == material.envi_nodes:
                break
        params = ('Name', 'Solar transmittance', 'Solar Reflectance', 'Visible reflectance', 'Infrared Hemispherical Emissivity', 'Infrared Transmittance', 'Thickness {m}',
                  'Conductivity {W/m-K}', 'Shade to glass distance {m}', 'Top opening multiplier', 'Top opening multiplier', 'Bottom opening multiplier', 'Left-side opening multiplier',
                  'Right-side opening multiplier', 'Air flow permeability')
        paramvs = ['{}-layer-{}'.format(material.name, ln)] + ['{:.3f}'.format(p) for p in (self.st, self.sr, self.vt, self.vr, self.ihe, self.it, 0.001 * self.thi, self.tc, 0.001 * self.sgd,
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
#            print('Shade control error {}'.format(e))
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

    
envimatnode_categories = [
        EnViMatNodeCategory("Type", "Type Node", items=[NodeItem("No_En_Mat_Con", label="Construction Node"),
                                                     NodeItem("No_En_Mat_PV", label="PV Node")]),
        EnViMatNodeCategory("Layer", "Layer Node", items=[NodeItem("No_En_Mat_Op", label="Opaque layer"),
                                                       NodeItem("No_En_Mat_Tr", label="Transparency layer"),
                                                       NodeItem("No_En_Mat_G", label="Gas layer")]), 
        EnViMatNodeCategory("Shading", "Shading Node", items=[NodeItem("No_En_Mat_Sh", label="Shading layer"),
                                                       NodeItem("No_En_Mat_Bl", label="Blind layer"),
                                                       NodeItem("No_En_Mat_Sc", label="Screen layer"),
                                                       NodeItem("No_En_Mat_SG", label="Switchable layer"),
                                                       NodeItem("No_En_Mat_ShC", label="Shading Control Node")]),
        EnViMatNodeCategory("Schedule", "Schedule Node", items=[NodeItem("No_En_Sched", label="Schedule")])]
