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
        self.inputs.new('ViLoc', 'Location in')
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
        self.inputs.new('ViLoc', 'Location in')
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
    
    def nodeupdate(self, context):
        nodecolour(self, self['exportstate'] != [str(x) for x in (self.animmenu)])
    
    def init(self, context):
        self.outputs.new('So_En_Geo', 'Geometry out')
        self['exportstate'] = ''
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
    '''Energy geometry out socket'''
    bl_idname = 'So_En_Geo'
    bl_label = 'EnVi Geometry'

    valid = ['EnVi Geometry']
    link_limit = 1

    def draw(self, context, layout, node, text):
        layout.label(text)

    def draw_color(self, context, node):
        return (0.0, 0.0, 1.0, 0.75)
        
####################### Vi Nodes Categories ##############################

vi_process = [NodeItem("No_Li_Geo", label="LiVi Geometry"), NodeItem("No_Li_Con", label="LiVi Context"), 
              NodeItem("No_En_Geo", label="EnVi Geometry")]
                
vi_edit = [NodeItem("No_Text", label="Text Edit")]
vi_analysis = [NodeItem("ViSPNode", label="Sun Path"), NodeItem("ViWRNode", label="Wind Rose"), 
             NodeItem("ViSVFNode", label="Sky View"), NodeItem("ViSSNode", label="Shadow map"),
             NodeItem("No_Li_Sim", label="LiVi Simulation")]

vigennodecat = []

vidisnodecat = []
vioutnodecat = []
vi_image = [NodeItem("No_Li_Im", label="LiVi Image"), NodeItem("No_Li_Gl", label="LiVi Glare"), NodeItem("No_Li_Fc", label="LiVi False-colour")]
vi_input = [NodeItem("No_Loc", label="VI Location")]

vinode_categories = [ViNodeCategory("Output", "Output Nodes", items=vioutnodecat), 
                     ViNodeCategory("Edit", "Edit Nodes", items=vi_edit), 
                     ViNodeCategory("Image", "Image Nodes", items=vi_image), 
                     ViNodeCategory("Display", "Display Nodes", items=vidisnodecat), 
                     ViNodeCategory("Generative", "Generative Nodes", items=vigennodecat), 
                     ViNodeCategory("Analysis", "Analysis Nodes", items=vi_analysis), 
                     ViNodeCategory("Process", "Process Nodes", items=vi_process), 
                     ViNodeCategory("Input", "Input Nodes", items=vi_input)]

class EnViNetwork(NodeTree):
    '''A node tree for the creation of EnVi advanced ventilation networks.'''
    bl_idname = 'EnViN'
    bl_label = 'EnVi Network'
    bl_icon = 'FORCE_WIND'
    nodetypes = {}

class EnViNodes:
    @classmethod
    def poll(cls, ntree):
        return ntree.bl_idname == 'EnViN'
    
class No_En_Net_Zone(Node, EnViNodes):
    '''Node describing a simulation zone'''
    bl_idname = 'No_En_Net_Zone'
    bl_label = 'Zone'
    bl_icon = 'SOUND'

    def zupdate(self, context):
        self.afs = 0
        obj = bpy.data.objects[self.zone]
        odm = obj.data.materials
        bfacelist = sorted([face for face in obj.data.polygons if get_con_node(odm[face.material_index]).envi_con_con == 'Zone'], key=lambda face: -face.center[2])
#        buvals = [retuval(odm[face.material_index]) for face in bfacelist]
        bsocklist = ['{}_{}_b'.format(odm[face.material_index].name, face.index) for face in bfacelist]
        sfacelist = sorted([face for face in obj.data.polygons if get_con_node(odm[face.material_index]).envi_afsurface == 1 and get_con_node(odm[face.material_index]).envi_con_type not in ('Window', 'Door')], key=lambda face: -face.center[2])
        ssocklist = ['{}_{}_s'.format(odm[face.material_index].name, face.index) for face in sfacelist]
        ssfacelist = sorted([face for face in obj.data.polygons if get_con_node(odm[face.material_index]).envi_afsurface == 1 and get_con_node(odm[face.material_index]).envi_con_type in ('Window', 'Door')], key=lambda face: -face.center[2])
        sssocklist = ['{}_{}_ss'.format(odm[face.material_index].name, face.index) for face in ssfacelist]

        [self.outputs.remove(oname) for oname in self.outputs if oname.bl_idname in ('EnViBoundSocket', 'EnViSFlowSocket', 'EnViSSFlowSocket')]
        [self.inputs.remove(iname) for iname in self.inputs if iname.bl_idname in ('EnViBoundSocket', 'EnViSFlowSocket', 'EnViSSFlowSocket')]

        for sock in bsocklist:
            self.outputs.new('EnViBoundSocket', sock).sn = sock.split('_')[-2]
            self.inputs.new('EnViBoundSocket', sock).sn = sock.split('_')[-2]
        for sock in ssocklist:
            self.afs += 1
            self.outputs.new('EnViSFlowSocket', sock).sn = sock.split('_')[-2]
            self.inputs.new('EnViSFlowSocket', sock).sn = sock.split('_')[-2]
        for sock in sssocklist:
            self.afs += 1
            self.outputs.new('EnViSSFlowSocket', sock).sn = sock.split('_')[-2]
            self.inputs.new('EnViSSFlowSocket', sock).sn = sock.split('_')[-2]                
#        for s, sock in enumerate(bsocklist):
#            self.outputs[sock].uvalue = '{:.4f}'.format(buvals[s])    
#            self.inputs[sock].uvalue = '{:.4f}'.format(buvals[s]) 
        self.vol_update(context)
        
    def vol_update(self, context):
        obj = bpy.data.objects[self.zone]      
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
        self.inputs.new('EnViHvacSocket', 'HVAC')
        self.inputs.new('EnViOccSocket', 'Occupancy')
        self.inputs.new('EnViEqSocket', 'Equipment')
        self.inputs.new('EnViInfSocket', 'Infiltration')
        self.inputs.new('EnViSchedSocket', 'TSPSchedule')
        self.inputs.new('EnViSchedSocket', 'VASchedule')

    def update(self):
        sflowdict = {'EnViSFlowSocket': 'Envi surface flow', 'EnViSSFlowSocket': 'Envi sub-surface flow'}
        [bi, si, ssi, bo, so , sso] = [1, 1, 1, 1, 1, 1]
                
        try:
            for inp in [inp for inp in self.inputs if inp.bl_idname in ('EnViBoundSocket', 'EnViSFlowSocket', 'EnViSSFlowSocket')]:
                self.outputs[inp.name].hide = True if inp.links and self.outputs[inp.name].bl_idname == inp.bl_idname else False
    
            for outp in [outp for outp in self.outputs if outp.bl_idname in ('EnViBoundSocket', 'EnViSFlowSocket', 'EnViSSFlowSocket')]:
                self.inputs[outp.name].hide = True if outp.links and self.inputs[outp.name].bl_idname == outp.bl_idname else False
    
            for inp in [inp for inp in self.inputs if inp.bl_idname in ('EnViBoundSocket', 'EnViSFlowSocket', 'EnViSSFlowSocket')]:
                if inp.bl_idname == 'EnViBoundSocket' and not inp.hide and not inp.links:
                    bi = 0
                elif inp.bl_idname in sflowdict:
                    if (not inp.hide and not inp.links) or (inp.links and inp.links[0].from_node.bl_label != sflowdict[inp.bl_idname]):
                        si = 0
                        if inp.links:
                            remlink(self, [inp.links[0]])    
            
            for outp in [outp for outp in self.outputs if outp.bl_idname in ('EnViBoundSocket', 'EnViSFlowSocket', 'EnViSSFlowSocket')]:
                if outp.bl_idname == 'EnViBoundSocket' and not outp.hide and not outp.links:
                    bo = 0
                elif outp.bl_idname  in sflowdict:
                    if (not outp.hide and not outp.links) or (outp.links and outp.links[0].to_node.bl_label != sflowdict[outp.bl_idname]):
                        so = 0
                        if outp.links:
                            remlink(self, [outp.links[0]])

        except Exception as e:
            print("Don't panic")
            
        self.alllinked = 1 if all((bi, si, ssi, bo, so, sso)) else 0
        nodecolour(self, self.errorcode())
        
    def uvsockupdate(self):
        for sock in self.outputs:
            socklink(sock, self['nodeid'].split('@')[1])
            if sock.bl_idname == 'EnViBoundSocket':
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
            row.label('Error - {}'.format(self.errorcode()))
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
    
class No_En_Mat_Con(Node, EnViMatNodes):
    '''Node defining the EnVi material construction'''
    bl_idname = 'No_En_Mat_Con'
    bl_label = 'EnVi Construction'
    bl_icon = 'FORCE_WIND'
    
    def con_update(self, context):
        if len(self.inputs) == 3:
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
            get_mat(self, 0).envi_type = self.envi_con_type        
            self.update()
    
    def frame_update(self, context):
        if self.fclass in ("0", "1"):
            for link in self.inputs['Outer frame layer'].links:
                self.id_data.links.remove(link)
            self.inputs['Outer frame layer'].hide = True
        else:
            self.inputs['Outer frame layer'].hide = False

        self.update()
        
    def active_update(self, context):
        if self.active:
            for node in [n for n in self.id_data.nodes if n.bl_idname == 'EnViCon' and n != self]:
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
                                    default = "0", update = frame_update)
    
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
        self.inputs.new('EnViSchedSocket', 'PV Schedule')
        self.inputs['PV Schedule'].hide = True
        self.inputs.new('envi_ol_sock', 'Outer layer')
        self.inputs['Outer layer'].hide = True
        self.inputs.new('envi_f_sock', 'Outer frame layer')
        self.inputs['Outer frame layer'].hide = True
        
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
            get_mat(self, 1).envi_shading = 0

            while in_sock.links:
                node = in_sock.links[0].from_node

                if node.bl_idname not in ('envi_sl_node', 'envi_bl_node', 'envi_screen_node', 'envi_sgl_node'):                    
                    paramvs.append('{}-layer-{}'.format(self['matname'], n)) 
                    params.append(('Outside layer', 'Layer {}'.format(n))[n > 0])
                    ep_text += node.ep_write(n)  
                    self.resist += node.resist
                else:
                    get_mat(self, 1).envi_shading = 1
                    
                in_sock = node.inputs['Layer']
                n += 1
                
            ep_text += epentry('Construction', params, paramvs)
            
            if get_mat(self, 1).envi_shading:
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
#            if self.envi_simple_glazing:
                
                
                
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
    
envimatnode_categories = [
        EnViMatNodeCategory("Type", "Type Node", items=[NodeItem("EnViCon", label="Construction Node"),
                                                     NodeItem("envi_frame_node", label="Frame Node"),
                                                     NodeItem("envi_pv_node", label="PV Node")]),
        EnViMatNodeCategory("Layer", "Layer Node", items=[NodeItem("envi_ol_node", label="Opaque layer"),
                                                       NodeItem("envi_tl_node", label="Transparency layer"),
                                                       NodeItem("envi_gl_node", label="Gas layer")]), 
        EnViMatNodeCategory("Shading", "Shading Node", items=[NodeItem("envi_sl_node", label="Shading layer"),
                                                       NodeItem("envi_bl_node", label="Blind layer"),
                                                       NodeItem("envi_screen_node", label="Screen layer"),
                                                       NodeItem("envi_sgl_node", label="Switchable layer"),
                                                       NodeItem("envi_sc_node", label="Shading Control Node")]),
        EnViMatNodeCategory("Schedule", "Schedule Node", items=[NodeItem("EnViSched", label="Schedule")])]
