# -*- coding: utf-8 -*-socklink2
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
from collections import OrderedDict
from bpy.props import EnumProperty, FloatProperty, IntProperty, BoolProperty, StringProperty, FloatVectorProperty
from bpy.types import NodeTree, Node, NodeSocket
from nodeitems_utils import NodeCategory, NodeItem
from subprocess import Popen, PIPE
from .vi_func import socklink, socklink2, uvsocklink, uvsocklink2, newrow, epwlatilongi, nodeinputs, remlink, rettimes, sockhide, selobj
from .vi_func import nodecolour, facearea, retelaarea, iprop, bprop, eprop, fprop, retdates
from .vi_func import delobj, logentry, ret_camera_menu, ret_param, ret_empty_menu, ret_datab, epentry
from .livi_func import hdrsky, cbdmhdr, cbdmmtx, retpmap, validradparams, sunposlivi
from .envi_func import retrmenus, enresprops, epschedwrite, processf, get_mat, get_con_node
from .livi_export import livi_sun, livi_sky, livi_ground, hdrexport
from .envi_mat import envi_materials, envi_constructions, envi_embodied, envi_layer, envi_layertype, envi_elayertype, envi_eclasstype, envi_emattype
from numpy import array, stack, where, unique
from numpy import sum as nsum
from .vi_dicts import rpictparams, rvuparams, rtraceparams, rtracecbdmparams
import matplotlib
matplotlib.use('qt5agg', force=True)
cur_dir = os.getcwd()

try:
    # addonpath = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    # os.chdir(os.path.join(addonpath, 'Python', sys.platform, 'netgen'))
    import netgen
    from netgen.meshing import MeshingParameters, FaceDescriptor, Element2D, Mesh
    from netgen.stl import STLGeometry
    from pyngcore import SetNumThreads, TaskManager
    ng = 1
except Exception as e:
    print('Problem with Netgen installation: {}'.format(e))
    ng = 0

ofoam = 0

if sys.platform in ('darwin', 'win32'):
    dck_run = Popen('docker images --quiet', shell=True, stdout=PIPE)

    for line in dck_run.stdout.readlines():
        if 'dce34e9a03e1' in line.decode():
            ofoam = 1

flo_libs = [ng, ofoam]
os.chdir(cur_dir)
envi_mats = envi_materials()
envi_cons = envi_constructions()
# envi_ecs = envi_embodied()


class ViNetwork(NodeTree):
    '''A node tree for VI-Suite analysis.'''
    bl_idname = 'ViN'
    bl_label = 'VI-Suite Nodes'
    bl_icon = 'NODETREE'
#    viparams = {}


class ViNodes:
    @classmethod
    def poll(cls, ntree):
        return ntree.bl_idname == 'ViN'

# Parametric nodes


class No_Anim(Node, ViNodes):
    '''Node to automate changes in parameters'''
    bl_idname = 'No_Anim'
    bl_label = 'VI Parametric'
    bl_icon = 'ANIM'

    def retparams(self, context):
        if self.inputs[0].links:
            return [(p.identifier, p.description, p.identifier) for p in self.inputs[0].links[0].from_node.bl_rna.properties if p.is_skip_save]
        else:
            return [('None', 'None', 'None')]

    parameter: EnumProperty(name='', description='Parameter to be animated', items=retparams)
    anim_file: StringProperty(name='')

    def init(self, context):
        self.inputs.new('So_Anim', 'Parameter')

    def draw_buttons(self, context, layout):
        newrow(layout, "Parameter:", self, 'parameter')
        layout.prop_search(self, 'anim_file', bpy.data, 'texts', text='File', icon='TEXT')


# Input nodes

class No_Loc(Node, ViNodes):
    '''Node describing a geographical location manually or with an EPW file'''
    bl_idname = 'No_Loc'
    bl_label = 'VI Location'
    bl_icon = 'WORLD'

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
                    try:
                        for wfl in wf.readlines():
                            if wfl.split(',')[0].upper() == 'LOCATION':
                                entries.append((wfile, '{} - {}'.format(wfl.split(',')[3], wfl.split(',')[1]), 'Weather Location'))
                                break

                            elif wfl.split(',')[0].upper() == "B'LOCATION":
                                with open(wfile, 'rb') as wfb:
                                    for wfbl in wfb.readlines():
                                        wfl = wfbl.decode()

                                        if wfl.split(',')[0].upper()[2:] == 'LOCATION':
                                            entries.append((wfile, '{} - {}'.format(wfl.split(',')[3], wfl.split(',')[1]), 'Weather Location'))
                                            break

                                logentry("Byte formatting found in file {}. Attempting to read byte format. If it fails remove leading b', end ' and all /r line endings".format(wfile))

                    except (TypeError, UnicodeDecodeError):
                        logentry(f'Non-unicode character found in {wfile}')

            self['entries'] = entries if entries else [('None', 'None', 'None')]

            if os.path.isfile(self.weather):
                with open(self.weather, 'r') as epwfile:
                    self['frames'] = ['0']
                    llist = epwfile.readlines()
                    epwlines = llist[8:] if llist[0][:2] != "b'" else llist[8:-1]
                    epwcolumns = list(zip(*[epwline.split(',') for epwline in epwlines]))
                    svp.year = 2019 if len(epwlines) == 8760 else 2020
                    times = ('Month', 'Day', 'Hour', 'DOS')

                    for t, ti in enumerate([' '.join(epwcolumns[c]) for c in range(1, 4)] + [' '.join(['{}'.format(int(d/24) + 1) for d in range(len(epwlines))])]):
                        reslists.append(['0', 'Time', 'Time', times[t], ti])

                    for c in {"Temperature (degC)": 6, 'Humidity (%)': 8, "Direct Solar (W/m^2)": 14, "Diffuse Solar (W/m^2)": 15,
                              'Wind Direction (deg)': 20, 'Wind Speed (m/s)': 21}.items():
                        reslists.append(['0', 'Climate', 'Exterior', c[0], ' '.join([cdata for cdata in list(epwcolumns[c[1]])])])

                    self.outputs['Location out']['epwtext'] = epwfile.read()
                    self.outputs['Location out']['valid'] = ['Location', 'Vi Results']
            else:
                self.outputs['Location out']['epwtext'] = ''
                self.outputs['Location out']['valid'] = ['Location']

        socklink(self.outputs['Location out'], self.id_data.name)
        self['reslists'] = reslists
        (svp.latitude, svp.longitude) = epwlatilongi(context.scene, self) if self.loc == '1' and self.weather != 'None' else (svp.latitude, svp.longitude)

        for node in [link.to_node for link in self.outputs['Location out'].links]:
            node.nodeupdate(context)

    def retentries(self, context):
        try:
            return [tuple(e) for e in self['entries']]
        except Exception:
            return [('None', 'None', 'None')]

    weather: EnumProperty(name='', description="Weather file", items=retentries, options={'SKIP_SAVE'}, update=updatelatlong)
    weather_anim: BoolProperty(name='', description='Animate weather file', default=False)
    weather_anim_file: StringProperty(name='')
    loc: EnumProperty(items=[("0", "Manual", "Manual location"), ("1", "EPW ", "Get location from EPW file")], name="", description="Location", default="0", update=updatelatlong)
    maxws: FloatProperty(name="", description="Max wind speed", min=0, max=90, default=0)
    minws: FloatProperty(name="", description="Min wind speed", min=0, max=90, default=0)
    avws: FloatProperty(name="", description="Average wind speed", min=0, max=0, default=0)
    dsdoy: IntProperty(name="", description="", min=1, max=365, default=1)
    dedoy: IntProperty(name="", description="", min=1, max=365, default=365)

    def init(self, context):
        self.outputs.new('So_Vi_Loc', 'Location out')
        self.outputs.new('So_Anim', 'Parameter')
        self['entries'] = [('None', 'None', 'None')]

        try:
            NodeTree.get_from_context(context).use_fake_user = True
        except Exception:
            pass

    def update(self):
        for sock in self.outputs:
            socklink(sock, self.id_data.name)

        nodecolour(self, self.ready())

    def draw_buttons(self, context, layout):
        newrow(layout, "Source:", self, 'loc')

        if self.loc == "1":
            newrow(layout, "Weather file:", self, 'weather')

            if self.weather_anim:
                layout.prop_search(self, 'weather_anim_file', bpy.data, 'texts', text='File', icon='TEXT')

        else:
            newrow(layout, 'Latitude', context.scene.vi_params, "latitude")
            newrow(layout, 'Longitude', context.scene.vi_params, "longitude")

    def ready(self):
        if self.loc == '1' and not self.weather:
            return 1
        if any([link.to_node.bl_label in ('LiVi CBDM', 'EnVi Export') and self.loc != "1" for link in self.outputs['Location out'].links]):
            return 1
        return 0


class No_ASC_Import(Node, ViNodes):
    '''Node describing a LiVi geometry export node'''
    bl_idname = 'No_ASC_Import'
    bl_label = 'VI ASC Import'
    bl_icon = 'GRID'

    single: BoolProperty(name='', default=False)
    ascfile: StringProperty()
    clear_nodata: EnumProperty(name="", description="Deal with no data", items=[('0', 'Zero', 'Make no data zero'), ('1', 'Delete', 'Delete no data')], default='0')

    def draw_buttons(self, context, layout):
        newrow(layout, 'Single file:', self, 'single')
        newrow(layout, 'No data:', self, 'clear_nodata')
        row = layout.row()
        row.operator('node.ascimport', text='Import ASC')

# Export Nodes


class No_Li_Geo(Node, ViNodes):
    '''Node describing a LiVi geometry export node'''
    bl_idname = 'No_Li_Geo'
    bl_label = 'LiVi Geometry'
    bl_icon = 'OBJECT_DATA'

    def ret_params(self):
        return [str(x) for x in (self.animated, self.startframe, self.endframe, self.cpoint, self.offset, self.mesh)]

    def nodeupdate(self, context):
        context.scene.vi_params.vi_nodes = self.id_data
        nodecolour(self, self['exportstate'] != self.ret_params())

    cpoint: EnumProperty(items=[("0", "Faces", "Export faces for calculation points"), ("1", "Vertices", "Export vertices for calculation points")],
                         name="", description="Specify the calculation point geometry", default="0", update=nodeupdate)
    offset: FloatProperty(name="", description="Calc point offset", min=0.001, max=1, default=0.02, precision=3, update=nodeupdate)
    animated: BoolProperty(name="", description="Animated analysis", default=0, update=nodeupdate)
    startframe: IntProperty(name="", description="Start frame for animation", min=0, default=0, update=nodeupdate)
    endframe: IntProperty(name="", description="End frame for animation", min=0, default=0, update=nodeupdate)
    mesh: BoolProperty(name="", description="Radiance mesh geometry export", default=0, update=nodeupdate)
    triangulate: BoolProperty(name="", description="Triangulate mesh geometry for export", default=0, update=nodeupdate)

    def init(self, context):
        self['exportstate'] = ''
        self.outputs.new('So_Li_Geo', 'Geometry out')
        nodecolour(self, 1)

    def draw_buttons(self, context, layout):
        newrow(layout, 'Triangulate:', self, 'triangulate')
        newrow(layout, 'Mesh:', self, 'mesh')
        newrow(layout, 'Animated:', self, 'animated')

        if self.animated:
            row = layout.row()
            row.label(text='Frames:')
            col = row.column()
            subrow = col.row(align=True)
            subrow.prop(self, 'startframe')
            subrow.prop(self, 'endframe')

        newrow(layout, 'Result point:', self, 'cpoint')
        newrow(layout, 'Offset:', self, 'offset')
        row = layout.row()
        row.operator("node.ligexport", text="Export")

    def update(self):
        for sock in self.outputs:
            socklink(sock, self.id_data.name)

    def preexport(self, scene):
        self['Text'] = {}
        self['Options'] = {'offset': self.offset, 'fs': (scene.frame_current, self.startframe)[self.animated],
                           'fe': (scene.frame_current, self.endframe)[self.animated], 'cp': self.cpoint, 'anim': self.animated}

    def postexport(self, scene):
        self.id_data.use_fake_user = 1
        self['exportstate'] = self.ret_params()
        nodecolour(self, 0)


class No_Li_Sen(Node, ViNodes):
    '''Node for creating LiVi sensing geometry'''
    bl_idname = 'No_Li_Sen'
    bl_label = 'LiVi Sense'
    bl_icon = 'TEXTURE'

    def draw_buttons(self, context, layout):
        pass


class No_Li_Con(Node, ViNodes):
    '''Node for creating a LiVi context'''
    bl_idname = 'No_Li_Con'
    bl_label = 'LiVi Context'
    bl_icon = 'LIGHT_SUN'

    def ret_params(self):
        return ['{}'.format(x) for x in (self.contextmenu, self.spectrummenu, self.cbanalysismenu,
                self.animated, self.skymenu, self.shour, self.sdoy, self.startmonth, self.endmonth, self.damin, self.dasupp, self.dalux, self.daauto,
                self.ehour, self.edoy, self.interval, self.hdr, self.hdrname, self.skyname, self.resname, self.turb, self.mtxname, self.cbdm_start_hour,
                self.cbdm_end_hour, self.leed4, self.colour, self.cbdm_res, self.ay, self.sp)]

    def nodeupdate(self, context):
        scene = context.scene
        scene.vi_params.vi_nodes = self.id_data
        nodecolour(self, self['exportstate'] != self.ret_params())

        if self.edoy < self.sdoy:
            self.edoy = self.sdoy
        if self.edoy == self.sdoy:
            if self.ehour < self.shour:
                self.ehour = self.shour

        self['skynum'] = int(self.skymenu)
        suns = [ob for ob in scene.objects if ob.type == 'LIGHT' and ob.data.type == 'SUN' and ob.visible_get()]

        if self.contextmenu == 'Basic' and ((self.skyprog == '0' and self['skynum'] < 3) or (self.skyprog == '1' and self.epsilon > 1)):
            starttime = datetime.datetime(2015, 1, 1, int(self.shour), int((self.shour - int(self.shour))*60)) + \
                datetime.timedelta(self.sdoy - 1) if self['skynum'] < 3 else datetime.datetime(2013, 1, 1, 12)
            self['endframe'] = self.startframe + int(((24 * (self.edoy - self.sdoy) + self.ehour - self.shour)/self.interval)) if self.animated else [scene.frame_current]
            frames = range(self.startframe, self['endframe'] + 1) if self.animated else [scene.frame_current]
            scene.frame_start, scene.frame_end = self.startframe, frames[-1]

            if suns:
                sun = suns[0]
                sun['VIType'] = 'Sun'
                [delobj(bpy.context.view_layer, sun) for sun in suns[1:]]
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

        if sys.platform == 'win32' and (self.hdr or self.cbanalysismenu == '0') and self.cbdm_res == 3:
            self.cbdm_res = 2

    spectrumtype = [('0', "Visible", "Visible radiation spectrum calculation"), ('1', "Full", "Full radiation spectrum calculation")]
    skylist = [("0", "Sunny", "CIE Sunny Sky description"), ("1", "Partly Coudy", "CIE Sunny Sky description"),
               ("2", "Coudy", "CIE Partly Cloudy Sky description"), ("3", "DF Sky", "Daylight Factor Sky description")]

    contexttype = [('Basic', "Basic", "Basic analysis"), ('CBDM', "CBDM", "Climate based daylight modelling")]
    contextmenu: EnumProperty(name="", description="Context type", items=contexttype, default='Basic', update=nodeupdate)
    animated: BoolProperty(name="", description="Animated sky", default=False, update=nodeupdate)
    offset: FloatProperty(name="", description="Calc point offset", min=0.001, max=1, default=0.01, update=nodeupdate)
    spectrummenu: EnumProperty(name="", description="Visible/full radiation spectrum selection", items=spectrumtype, default='0', update=nodeupdate)
    skyprog: EnumProperty(name="", items=[('0', "Gensky", "Basic sky creation"), ('1', "Gendaylit", "Perez sky creation"),
                                          ("2", "HDR Sky", "HDR file sky"), ("3", "Radiance Sky", "Radiance file sky"), ("4", "None", "No Sky")],
                          description="Specify sky creation", default="0", update=nodeupdate)
    epsilon: FloatProperty(name="", description="Sky epsilon", min=1, max=8, default=6.3, options={'SKIP_SAVE'}, update=nodeupdate)
    delta: FloatProperty(name="", description="Sky delta", min=0.05, max=0.5, default=0.15, options={'SKIP_SAVE'}, update=nodeupdate)
    skymenu: EnumProperty(name="", items=skylist, description="Specify the type of sky for the simulation", default="0", update=nodeupdate)
    gref: FloatProperty(name="", description="Ground reflectance", min=0.0, max=1.0, default=0.0, options={'SKIP_SAVE'}, update=nodeupdate)
    gcol: FloatVectorProperty(size=3, name='', description="Ground colour", attr='Color', default=[0, 1, 0], subtype='COLOR', options={'SKIP_SAVE'}, update=nodeupdate)
    sdist: FloatProperty(name="", description="Blender sun distance", min=0.0, default=50.0, update=nodeupdate)
    shour: FloatProperty(name="", description="Start hour of simulation", min=0, max=23.99, default=12, subtype='TIME', unit='TIME', options={'ANIMATABLE'}, update=nodeupdate)
    sdoy: IntProperty(name="", description="Start day of simulation", min=1, max=365, default=1, update=nodeupdate)
    ehour: FloatProperty(name="", description="End hour of simulation", min=0, max=23.99, default=12, subtype='TIME', unit='TIME', update=nodeupdate)
    edoy: IntProperty(name="", description="End day of simulation", min=1, max=365, default=1, update=nodeupdate)
    interval: FloatProperty(name="", description="Site Latitude", min=1/60, max=24, default=1, update=nodeupdate)
    hdr: BoolProperty(name="", description="Export HDR panoramas", default=False, update=nodeupdate)
    skyname: StringProperty(name="", description="Name of the radiance sky file", default="", subtype="FILE_PATH", update=nodeupdate)
    resname: StringProperty()
    turb: FloatProperty(name="", description="Sky Turbidity", min=1.0, max=5.0, default=2.75, update=nodeupdate)
    cusacc: StringProperty(name="Custom parameters", description="Custom Radiance simulation parameters", default="", update=nodeupdate)
    buildstorey: EnumProperty(items=[("0", "Single", "Single storey building"), ("1", "Multi", "Multi-storey building")],
                              name="", description="Building storeys", default="0", update=nodeupdate)
    cbanalysistype = [('0', "Exposure", "LuxHours/Irradiance Exposure Calculation"), ('1', "Hourly irradiance", "Irradiance for each simulation time step"),
                      ('2', "DA/UDI/SDA/ASE", "Climate based daylighting metrics")]
    cbanalysismenu: EnumProperty(name="", description="Type of lighting analysis", items=cbanalysistype, default='0', update=nodeupdate)
    sourcetype = [('0', "EPW", "EnergyPlus weather file"), ('1', "HDR", "HDR sky file")]
    sourcetype2 = [('0', "EPW", "EnergyPlus weather file"), ('1', "VEC", "Generated vector file")]
    sourcemenu: EnumProperty(name="", description="Source type", items=sourcetype, default='0', update=nodeupdate)
    sourcemenu2: EnumProperty(name="", description="Source type", items=sourcetype2, default='0', update=nodeupdate)
    hdrname: StringProperty(name="", description="Name of the composite HDR sky file", default="", subtype="FILE_PATH", update=nodeupdate)
    hdrmap: EnumProperty(items=[("0", "Polar", "Latitude/longitude (equirectangular) format"), ("1", "Angular", "Light probe or angular map format")],
                         name="", description="Type of HDR panorama mapping", default="0", update=nodeupdate)
    hdrangle: FloatProperty(name="", description="HDR rotation (deg)", min=0, max=360, default=0, update=nodeupdate)
    hdrradius: FloatProperty(name="", description="HDR radius (m)", min=0, max=5000, default=1000, update=nodeupdate)
    mtxname: StringProperty(name="", description="Name of the calculated vector sky file", default="", subtype="FILE_PATH", update=nodeupdate)
    weekdays: BoolProperty(name='', default=False, update=nodeupdate)
    cbdm_start_hour:  IntProperty(name='', description="Hour of the day (1 being the first hour: midnight to 1am)", default=8, min=1, max=24, update=nodeupdate)
    cbdm_end_hour:  IntProperty(name='', description="Hour of the day (24 being the last hour: 11pm to midnight)", default=20, min=1, max=24, update=nodeupdate)
    cbdm_res: IntProperty(name='', default=1, min=1, max=3, update=nodeupdate)
    dalux:  IntProperty(name='lux', default=300, min=1, max=2000, update=nodeupdate)
    damin: IntProperty(name='lux', default=100, min=1, max=2000, update=nodeupdate)
    dasupp: IntProperty(name='lux', default=300, min=1, max=2000, update=nodeupdate)
    daauto: IntProperty(name='lux', default=3000, min=1, max=5000, update=nodeupdate)
    sdamin: IntProperty(name='lux', default=300, min=1, max=2000, update=nodeupdate)
    asemax: IntProperty(name='lux', default=1000, min=1, max=2000, update=nodeupdate)
    startmonth: IntProperty(name='', default=1, min=1, max=12, description='Start Month', update=nodeupdate)
    endmonth: IntProperty(name='', default=12, min=1, max=12, description='End Month', update=nodeupdate)
    startframe: IntProperty(name='', default=0, min=0, description='Start Frame', update=nodeupdate)
    leed4: BoolProperty(name='', description='LEED v4 Compliance',  default=False, update=nodeupdate)
    ay: BoolProperty(name='', description='All year simulation',  default=False, update=nodeupdate)
    colour: BoolProperty(name='', description='Coloured Gendaylit sky',  default=False, update=nodeupdate)
    sp: BoolProperty(name='', description='Split channels',  default=False, update=nodeupdate)

    def init(self, context):
        self['exportstate'], self['skynum'] = '', 0
        self['whitesky'] = ("void glow sky_glow \n0 \n0 \n4 1 1 1 0 \nsky_glow source sky \n0 \n0 \n4 0 0 1 180 \n"
                            "void glow ground_glow \n0 \n0 \n4 1 1 1 0 \nground_glow source ground \n0 \n0 \n4 0 0 -1 180\n\n")
        self.outputs.new('So_Li_Con', 'Context out')
        self.inputs.new('So_Vi_Loc', 'Location in')
        self.outputs['Context out'].hide = True
        self.outputs.new('So_Anim', 'Parameter')
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
                    newrow(layout, "Sun distance:", self, 'sdist')
                    newrow(layout, "Start hour {}:{}:".format(int(self.shour), int((self.shour*60) % 60)), self, 'shour')
                    newrow(layout, 'Start day {}/{}:'.format(sdate.day, sdate.month), self, "sdoy")
                    newrow(layout, "Animation;", self, 'animated')

                    if self.animated:
                        newrow(layout, "Start frame:", self, 'startframe')
                        row = layout.row()
                        row.label(text='End frame:')
                        row.label(text='{}'.format(self['endframe']))
                        newrow(layout, "End hour {}:{}:".format(int(self.ehour), int((self.ehour*60) % 60)), self, 'ehour')
                        newrow(layout, 'End day {}/{}:'.format(edate.day, edate.month), self, "edoy")
                        newrow(layout, "Interval (hours):", self, 'interval')

                    newrow(layout, "Turbidity", self, 'turb')

            elif self.skyprog == '1':
                newrow(layout, "Spectrum:", self, 'spectrummenu')
                newrow(layout, 'Colour sky', self, "colour")
                newrow(layout, "Epsilon:", self, 'epsilon')
                newrow(layout, "Delta:", self, 'delta')
                newrow(layout, "Ground ref:", self, 'gref')
                newrow(layout, "Ground col:", self, 'gcol')
                newrow(layout, "Start hour {}:{}:".format(int(self.shour), int((self.shour*60) % 60)), self, 'shour')
                newrow(layout, 'Start day {}/{}:'.format(sdate.day, sdate.month), self, "sdoy")
                newrow(layout, "Animation;", self, 'animated')

                if self.animated:
                    newrow(layout, "Start frame:", self, 'startframe')
                    row = layout.row()
                    row.label(text='End frame:')
                    row.label(text='{}'.format(self['endframe']))
                    newrow(layout, "End hour {}:{}:".format(int(self.ehour), int((self.ehour*60) % 60)), self, 'ehour')
                    newrow(layout, 'End day {}/{}:'.format(edate.day, edate.month), self, "edoy")
                    newrow(layout, "Interval (hours):", self, 'interval')

            elif self.skyprog == '2':
                row = layout.row()
                row.operator('node.hdrselect', text='HDR select')
                row.prop(self, 'hdrname')
                newrow(layout, "HDR format:", self, 'hdrmap')
                newrow(layout, "HDR rotation:", self, 'hdrangle')
                newrow(layout, "HDR radius:", self, 'hdrradius')

            elif self.skyprog == '3':
                row = layout.row()
                row.prop_search(self, 'skyname', bpy.data, 'texts', text='', icon='NONE')

            row = layout.row()

            if self.skyprog in ("0", "1", "3"):
                newrow(layout, 'HDR:', self, 'hdr')

            newrow(layout, 'Split channels:', self, 'sp')

        elif self.contextmenu == 'CBDM':
            newrow(layout, 'Type:', self, 'cbanalysismenu')

            if self.cbanalysismenu == '0' and self.sourcemenu != '1':
                newrow(layout, "Spectrum:", self, 'spectrummenu')

            newrow(layout, 'All year:', self, 'ay')
            newrow(layout, 'Weekdays only:', self, 'weekdays')

            if self.cbanalysismenu == '2':
                newrow(layout, 'LEED v4:', self, 'leed4')

            if self.cbanalysismenu in ('0', '1') or (self.cbanalysismenu == '2' and not self.leed4):
                if not self.ay:
                    newrow(layout, 'Start day {}/{}:'.format(sdate.day, sdate.month), self, "sdoy")
                    newrow(layout, 'End day {}/{}:'.format(edate.day, edate.month), self, "edoy")
                    newrow(layout, 'Start hour:', self, 'cbdm_start_hour')
                    newrow(layout, 'End hour:', self, 'cbdm_end_hour')

                if self.cbanalysismenu == '2':
                    row = layout.row()
                    row.label(text="--")
                    newrow(layout, '(s)DA (Min):', self, 'dalux')
                    newrow(layout, 'UDI Low (Max):', self, 'damin')
                    newrow(layout, 'UDI Supplementry (Max):', self, 'dasupp')
                    newrow(layout, 'UDI Autonomous (Max):', self, 'daauto')
                    newrow(layout, 'ASE level:', self, 'asemax')
                    row = layout.row()
                    row.label(text="--")

            elif self.cbanalysismenu == '2' and self.leed4:
                newrow(layout, 'Start hour:', self, 'cbdm_start_hour')
                newrow(layout, 'End hour:', self, 'cbdm_end_hour')

            if self.cbanalysismenu == '0':
                newrow(layout, 'Source file:', self, 'sourcemenu')
            else:
                newrow(layout, 'Source file:', self, 'sourcemenu2')

            row = layout.row()

            if self.sourcemenu2 == '1' and self.cbanalysismenu in ('1', '2'):
                newrow(layout, "MTX file:", self, 'mtxname')

            elif self.sourcemenu == '1' and self.cbanalysismenu == '0':
                newrow(layout, "HDR file:", self, 'hdrname')
            else:
                newrow(layout, 'Resolution:', self, 'cbdm_res')

            if self.cbanalysismenu != '0':
                newrow(layout, 'HDR:', self, 'hdr')

        if self.contextmenu == 'Basic':
            if int(self.skymenu) > 2 or int(self.skyprog) > 1 or (int(self.skymenu) < 3 and self.inputs['Location in'].links):
                row = layout.row()
                row.operator("node.liexport", text="Export")

        elif (self.contextmenu == 'CBDM' and self.cbanalysismenu == '0' and self.sourcemenu == '1'):
            if os.path.isfile(bpy.path.abspath(self.hdrname)):
                row = layout.row()
                row.operator("node.liexport", text="Export")
            else:
                row = layout.row()
                row.label(text="ERROR: No valid HDR file")

        elif (self.contextmenu == 'CBDM' and self.cbanalysismenu != '0' and self.sourcemenu2 == '1'):
            if os.path.isfile(bpy.path.abspath(self.mtxname)):
                row = layout.row()
                row.operator("node.liexport", text="Export")
            else:
                row = layout.row()
                row.label(text="ERROR: No valid MTX file")

        elif self.inputs['Location in'].links and self.inputs['Location in'].links[0].from_node.loc == '1' and self.inputs['Location in'].links[0].from_node.weather != 'None':
            row = layout.row()
            row.operator("node.liexport", text="Export")
        else:
            row = layout.row()
            row.label(text="ERROR: No valid EPW file")

    def update(self):
        for sock in self.outputs:
            socklink(sock, self.id_data.name)

        if self.inputs.get('Location in'):
            self.nodeupdate(bpy.context)

    def preexport(self):
        if self.contextmenu == 'CBDM' and self.leed4:
            self.asemax = 1000
            self.sdamin = 300
            self.daauto = 3000
            self.sdoy = 1
            self.edoy = 365

            if self.cbdm_start_hour > 9:
                self.cbdm_start_hour = 9
            if self.cbdm_end_hour < 16:
                self.cbdm_end_hour = 16

        if self.contextmenu == 'Basic':
            (shour, ehour) = (self.shour, self.ehour)
            (sdoy, edoy) = (self.sdoy, self.edoy)

        elif self.contextmenu == 'CBDM':
            if self.cbanalysismenu == '2' and self.leed4:
                (shour, ehour) = (self.cbdm_start_hour - 1, self.cbdm_end_hour - 1)
            elif self.ay:
                (shour, ehour) = (0, 23)
            else:
                (shour, ehour) = (self.cbdm_start_hour - 1, self.cbdm_end_hour - 1)
            (sdoy, edoy) = (self.sdoy, self.edoy) if not self.ay else (1, 365)

        interval = 1
        starttime = datetime.datetime(2015, 1, 1, 0) + datetime.timedelta(days=sdoy - 1) + datetime.timedelta(hours=shour)

        if (self.contextmenu == 'CBDM' and not self.leed4) or (self.contextmenu == 'Basic' and self.animated):
            endtime = datetime.datetime(2015, 1, 1, 0) + datetime.timedelta(days=edoy - 1) + datetime.timedelta(hours=ehour)
        elif self.contextmenu == 'CBDM':
            starttime = datetime.datetime(2015, 1, 1, 0) + datetime.timedelta(hours=shour)
            endtime = datetime.datetime(2015, 1, 1, 0) + datetime.timedelta(days=edoy - 1) + datetime.timedelta(hours=ehour)
        else:
            endtime = starttime

        times = [starttime]
        ctime = starttime

        while ctime < endtime:
            ctime += datetime.timedelta(hours=interval)

            if self.contextmenu == 'CBDM':

                if shour <= ctime.hour <= ehour:
                    times.append(ctime)
            else:
                times.append(ctime)

        self.times = [t for t in times if t.weekday() <= (6, 4)[self.weekdays]]
        self.starttime = times[0]
        self.endtime = times[-1]
        self['skynum'] = int(self.skymenu)
        self['hours'] = 0 if not self.animated or int(self.skymenu) > 2 else (self.endtime-self.starttime).seconds/3600
        self['epwbase'] = os.path.splitext(os.path.basename(self.inputs['Location in'].links[0].from_node.weather)) if self.inputs['Location in'].links else ''
        self['Text'], self['Options'] = {}, {}
        self['watts'] = 1 if self.contextmenu == "CBDM" and ((self.cbanalysismenu == '0' and self.spectrummenu == '1') or self.cbanalysismenu == '1') else 0

    def export(self, scene, export_op):
        self._valid = 1
        svp = scene.vi_params
        self.startframe = self.startframe if self.animated and self.contextmenu == 'Basic' else scene.frame_current
        self['endframe'] = self.startframe + int(((24 * (self.edoy - self.sdoy) + self.ehour - self.shour)/self.interval)) if self.contextmenu == 'Basic' and self.animated else scene.frame_current
        self['mtxfile'] = ''
        self['mtxfilens'] = ''
        self['preview'] = 0

        if self.contextmenu == "Basic":
            self['preview'] = 1

            if self.skyprog in ('0', '1'):
                self['skytypeparams'] = ("+s", "+i", "-c", "-b 22.86 -c")[self['skynum']] if self.skyprog == '0' else "-P {} {} -O {} {}".format(self.epsilon,
                                                                                                                                                 self.delta,
                                                                                                                                                 int(self.spectrummenu),
                                                                                                                                                 ('', '-C')[self.colour])

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
                hdr_loc = bpy.path.abspath(self.hdrname)
                if self.hdrname and os.path.isfile(hdr_loc):
                    if self.hdrname not in bpy.data.images:
                        bpy.data.images.load(self.hdrname)

                    self['Text'][str(scene.frame_current)] = hdrsky(hdr_loc, self.hdrmap, self.hdrangle, self.hdrradius)
                else:
                    export_op.report({'ERROR'}, "Not a valid HDR file")
                    return 'Error'

            elif self.skyprog == '3':
                if self.skyname and bpy.data.texts.get(self.skyname):
                    self['Text'][str(scene.frame_current)] = bpy.data.texts[self.skyname].as_string()

                    if self.hdr:
                        hdrexport(scene, 0, scene.frame_current, self, bpy.data.texts[self.skyname].as_string())
                else:
                    export_op.report({'ERROR'}, "Not a valid Radiance sky file")
                    return 'Error'

            elif self.skyprog == '4':
                self['Text'][str(scene.frame_current)] = ''

        elif self.contextmenu == "CBDM":
            if (self.cbanalysismenu == '0' and self.sourcemenu == '0') or (self.cbanalysismenu != '0' and self.sourcemenu2 == '0'):
                (self['mtxfile'], self['mtxfilens']) = cbdmmtx(self, scene, self.inputs['Location in'].links[0].from_node, export_op)

            elif self.cbanalysismenu != '0' and self.sourcemenu2 == '1':
                self['mtxfile'] = bpy.path.abspath(self.mtxname)
                matrix_ns = '.'.join(self['mtxfile'].split('.')[:-2] + [self['mtxfile'].split('.')[-2] + 'ns'] + [self['mtxfile'].split('.')[-1]])
                self['mtxfilens'] = matrix_ns

            if self.cbanalysismenu == '0':
                self['preview'] = 1
                self['Text'][str(scene.frame_current)] = cbdmhdr(self, scene, export_op)
            else:
                self['Text'][str(scene.frame_current)] = ("void glow sky_glow \n0 \n0 \n4 1 1 1 0 \nsky_glow source sky \n0 \n0 \n4 0 0 1 180 \n"
                                                          "void glow ground_glow \n0 \n0 \n4 1 1 1 0 \nground_glow source ground \n0 \n0 \n4 0 0 -1 180\n\n")

                if self.sourcemenu2 == '0':
                    with open("{}.mtx".format(os.path.join(svp['viparams']['newdir'], self['epwbase'][0])), 'r') as mtxfile:
                        self['Options']['MTX'] = mtxfile.read()
                else:
                    with open(bpy.path.abspath(self.mtxname), 'r') as mtxfile:
                        for line in mtxfile.readlines():
                            if line.split('=')[0] == 'NROWS':
                                self.cbdm_res = (146, 578, 2306).index(int(line.split('=')[1])) + 1
                                break
                            # elif line.split('=')[0] == 'NCOLS':
                            #     if len(self.times) != int(line.split('=')[1]):
                            #         export_op.report({'ERROR'}, "Outdated MTX file")
                            #         self._valid = 0
                            #         return

                        self['Options']['MTX'] = mtxfile.read()

                if self.hdr:
                    cbdmhdr(self, scene, export_op)

    # def check_mtx(self):
    #     if self.cbanalysismenu != '0' and self.sourcemenu2 == '1':

    def postexport(self):
        (csh, ceh) = (self.cbdm_start_hour, self.cbdm_end_hour) if not self.ay or (self.cbanalysismenu == '2' and self.leed4) else (1, 24)
        (sdoy, edoy) = (self.sdoy, self.edoy) if self.contextmenu == '0' or not self.ay else (1, 365)
        typedict = {'Basic': '0', 'CBDM': self.cbanalysismenu}
        basic_unit = 'W/m2' if self.sp else (("Lux", "DF (%)")[self.skyprog == '0' and self.skymenu == '3'], 'W/m2 (f)')[self.skyprog == '1' and self.spectrummenu == '1']
        unitdict = {'Basic': basic_unit,
                    'CBDM': (('klxh', 'kWh (f)')[int(self.spectrummenu)], 'kWh (f)', 'DA (%)')[int(self.cbanalysismenu)]}
        self['Options'] = {'Context': self.contextmenu, 'Preview': self['preview'], 'Type': typedict[self.contextmenu],
                           'fs': self.startframe, 'fe': self['endframe'], 'anim': self.animated, 'shour': self.shour,
                           'sdoy': self.sdoy, 'ehour': self.ehour, 'edoy': self.edoy, 'interval': self.interval,
                           'cbanalysis': self.cbanalysismenu, 'unit': unitdict[self.contextmenu],
                           'damin': self.damin, 'dalux': self.dalux, 'dasupp': self.dasupp,
                           'daauto': self.daauto, 'asemax': self.asemax, 'cbdm_sh': csh,
                           'cbdm_eh': ceh, 'cbdm_sd': sdoy, 'cbdm_ed': edoy, 'weekdays': (7, 5)[self.weekdays],
                           'sourcemenu': (self.sourcemenu, self.sourcemenu2)[self.cbanalysismenu not in ('2', '3', '4', '5')],
                           'mtxfile': self['mtxfile'], 'mtxfilens': self['mtxfilens'], 'times': [t.strftime("%d/%m/%y %H:%M:%S") for t in self.times],
                           'leed4': self.leed4, 'colour': self.colour, 'cbdm_res': (146, 578, 2306)[self.cbdm_res - 1],
                           'sm': self.skymenu, 'sp': self.skyprog, 'ay': self.ay}

        # if self._valid:
        nodecolour(self, 0)
        self.outputs['Context out'].hide = False
        self['exportstate'] = self.ret_params()


class No_Vi_Im(Node, ViNodes):
    '''Image node'''
    bl_idname = 'No_Vi_Im'
    bl_label = 'Image'
    bl_icon = 'IMAGE'

    def nodeupdate(self, context):
        self['images'] = [bpy.data.images[self.image].filepath]
        nodecolour(self, self['exportstate'] != [str(x) for x in (self.image)])

    image: StringProperty(description="Select image", update=nodeupdate)

    def init(self, context):
        self['exportstate'] = ''
        self.outputs.new('So_Li_Im', 'Image')

    def draw_buttons(self, context, layout):
        layout.prop_search(self, 'image', bpy.data, 'images', text='Image', icon='NONE')

    def update(self):
        for sock in self.outputs:
            socklink(sock, self.id_data.name)

# Image nodes


class No_Li_Im(Node, ViNodes):
    '''Node describing a LiVi image generation'''
    bl_idname = 'No_Li_Im'
    bl_label = 'LiVi Image'
    bl_icon = 'IMAGE'

    def ret_params(self):
        return [str(x) for x in (self.camera, self.basename, self.illu, self.fisheye, self.fov, self.mp, self.processors,
                self.processes, self.cusacc, self.simacc, self.pmap, self.pmapgno, self.pmapcno, self.x, self.y)]

    def nodeupdate(self, context):
        if self.processors > int(context.scene.vi_params['viparams']['nproc']):
            self.processors = int(context.scene.vi_params['viparams']['nproc'])

        nodecolour(self, self['exportstate'] != self.ret_params())

        if bpy.data.objects.get(self.camera):
            context.scene.camera = bpy.data.objects[self.camera]

        if self.simacc == '3':
            self.validparams = validradparams(self.cusacc)

    startframe: IntProperty(name='', default=0)
    endframe: IntProperty(name='', default=0)
    cusacc: StringProperty(name="Custom parameters", description="Custom Radiance simulation parameters", default="", update=nodeupdate)
    simacc: EnumProperty(items=[("0", "Low", "Low accuracy and high speed (preview)"), ("1", "Medium", "Medium speed and accuracy"),
                                ("2", "High", "High but slow accuracy"), ("3", "Custom", "Edit Radiance parameters")],
                         name="", description="Simulation accuracy", default="0", update=nodeupdate)
    pmap: BoolProperty(name='', default=False, update=nodeupdate)
    pmapgno: IntProperty(name='', description="Number of global photons", default=50000)
    pmapcno: IntProperty(name='', description="Number of caustic photons", default=0)
    pmapvno: IntProperty(name='', description="Number of visualised photons", min=100, max=5000, default=500)
    pmapoptions: StringProperty(name="", description="Additional pmap parameters", default="", update=nodeupdate)
    pmappreview: BoolProperty(name='', default=0, update=nodeupdate)
    bfv: BoolProperty(name='', description="Turn on back face visibility (may cause light bleed but deals with planar geometry)", default=True, update=nodeupdate)
    x: IntProperty(name='', min=1, max=10000, default=2000, update=nodeupdate)
    y: IntProperty(name='', min=1, max=10000, default=1000, update=nodeupdate)
    basename: StringProperty(name="", description="Base name for image files", default="", update=nodeupdate)
    run: BoolProperty(name='', default=False)
    illu: BoolProperty(name='', default=True, update=nodeupdate)
    validparams: BoolProperty(name='', default=True)
    mp: BoolProperty(name='', default=False, update=nodeupdate)
    camera: StringProperty(description="Select camera", update=nodeupdate)
    fisheye: BoolProperty(name='', default=0, update=nodeupdate)
    fov: FloatProperty(name='', default=180, min=1, max=180, update=nodeupdate)
    processors: IntProperty(name='', default=1, min=1, max=128, update=nodeupdate)
    processes: IntProperty(name='', default=1, min=1, max=1000, update=nodeupdate)
    normal: BoolProperty(name='', description="Generate denoising normal map", default=False)
    albedo: BoolProperty(name='', description="Generate denoising albedo map", default=False)

    def retframes(self):
        try:
            return min([c['fs'] for c in (self.inputs['Context in'].links[0].from_node['Options'], self.inputs['Geometry in'].links[0].from_node['Options'])]),\
                    max([c['fe'] for c in (self.inputs['Context in'].links[0].from_node['Options'], self.inputs['Geometry in'].links[0].from_node['Options'])])
        except Exception:
            return 0, 0

    def init(self, context):
        self['exportstate'] = ''
        self.inputs.new('So_Li_Geo', 'Geometry in')
        self.inputs.new('So_Li_Con', 'Context in')
        self.outputs.new('So_Li_Im', 'Image')

    def draw_buttons(self, context, layout):
        sf, ef = self.retframes()
        row = layout.row()
        row.label(text='Frames: {} - {}'.format(sf, ef))
        layout.prop_search(self, 'camera', bpy.data, 'cameras', text='Camera', icon='NONE')

        if all([sock.links for sock in self.inputs]) and self.camera:
            newrow(layout, 'Base name:', self, 'basename')
            newrow(layout, 'Illuminance:', self, 'illu')
            newrow(layout, 'Fisheye:', self, 'fisheye')

            if self.fisheye:
                newrow(layout, 'FoV:', self, 'fov')

            newrow(layout, 'Accuracy:', self, 'simacc')
            newrow(layout, 'Normal:', self, 'normal')
            newrow(layout, 'Albedo:', self, 'albedo')

            if self.simacc == '3':
                row = layout.row()
                row.prop(self, 'cusacc')

            newrow(layout, 'Photon map:', self, 'pmap')

            if self.pmap:
                newrow(layout, 'Global photons:', self, 'pmapgno')
                newrow(layout, 'Caustic photons:', self, 'pmapcno')
                newrow(layout, 'Back face visability:', self, 'bfv')
                newrow(layout, 'Photons options:', self, 'pmapoptions')
                newrow(layout, 'Preview photons:', self, 'pmappreview')

                if self.pmappreview:
                    newrow(layout, 'Visualised photons:', self, 'pmapvno')

            if self.simacc != '3' or (self.simacc == '3' and self.validparams) and not self.run:
                row = layout.row()
                row.operator("node.radpreview", text='Preview')

            newrow(layout, 'X resolution:', self, 'x')
            newrow(layout, 'Y resolution:', self, 'y')

            if sys.platform != 'win32':
                newrow(layout, 'Multi-thread:', self, 'mp')

                if self.mp:
                    newrow(layout, 'Processors:', self, 'processors')
                    newrow(layout, 'Processes:', self, 'processes')

            if (self.simacc != '3' or (self.simacc == '3' and self.validparams)) and not self.run:
                row = layout.row()
                row.operator("node.radimage", text='Image')

    def update(self):
        for sock in self.outputs:
            socklink(sock, self.id_data.name)

        self.run = 0

    def presim(self):
        self.time = datetime.datetime.now()
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
        self['rpictparams'] = ' {} '.format(self.cusacc) if self.simacc == '3' else ''.join([' {} {} '.format(k, rpictparams[k][int(self.simacc)]) for k in rpictparams])
        self['rvuparams'] = ' {} '.format(self.cusacc) if self.simacc == '3' else ''.join([' {} {} '.format(k, rvuparams[k][int(self.simacc)]) for k in rvuparams])
        self['basename'] = self.basename if self.basename else 'image'

        for frame in self['frames']:
            sf = str(frame)
            scene.frame_set(frame)
            scene.camera = bpy.data.objects[self.camera.lstrip()]
            cam = bpy.data.objects[self.camera.lstrip()]
            cang = cam.data.angle if not self.fisheye else self.fov
            2 * math.atan((0.5 * self.y) / (0.5 * self.x / math.tan(cang / 2)))
            vh = 180/math.pi * cang if self.x >= self.y else 180/math.pi * 2 * math.atan((0.5 * self.x) / (0.5 * self.y / math.tan(cang / 2)))
            vv = 180/math.pi * cang if self.x < self.y else 180/math.pi * 2 * math.atan((0.5 * self.y) / (0.5 * self.x / math.tan(cang / 2)))
            vd = (0.001, 0, -1*cam.matrix_world[2][2]) if (round(-1*cam.matrix_world[0][2], 3), round(-1*cam.matrix_world[1][2], 3)) == (0.0, 0.0) else [-1*cam.matrix_world[i][2] for i in range(3)]
            pmaps.append(self.pmap)
            self['pmapgnos'][sf] = self.pmapgno
            self['pmapcnos'][sf] = self.pmapcno
            self['pmparams'][sf]['amentry'], self['pmparams'][sf]['pportentry'], self['pmparams'][sf]['gpentry'], self['pmparams'][sf]['cpentry'], self['pmparams'][sf]['gpfileentry'], self['pmparams'][sf]['cpfileentry'] = retpmap(self, frame, scene)
            # amentry, pportentry, gpentry, cpentry, gpfileentry, cpfileentry
            if self.fisheye and self.fov == 180:
                self['viewparams'][sf]['-vth'] = ''

            (self['viewparams'][sf]['-vh'], self['viewparams'][sf]['-vv']) = (self.fov, self.fov) if self.fisheye else ('{:.3f}'.format(vh), '{:.3f}'.format(vv))
            self['viewparams'][sf]['-vd'] = ' '.join(['{:.3f}'.format(v) for v in vd])
            self['viewparams'][sf]['-x'], self['viewparams'][sf]['-y'] = self.x, self.y

            if self.mp:
                self['viewparams'][sf]['-X'], self['viewparams'][sf]['-Y'] = self.processes, 1

            self['viewparams'][sf]['-vp'] = '{0[0]:.3f} {0[1]:.3f} {0[2]:.3f}'.format(cam.location)
            self['viewparams'][sf]['-vu'] = '{0[0]:.3f} {0[1]:.3f} {0[2]:.3f}'.format(cam.matrix_world.to_quaternion()@mathutils.Vector((0, 1, 0)))

            if self.illu:
                self['viewparams'][sf]['-i'] = ''

        self['pmaps'] = pmaps
        self.run = 1
        nodecolour(self, 1)

    def postsim(self, images):
        self['images'] = images
        self.run = 0
        self['exportstate'] = self.ret_params()
        logentry('Time to render: {}'.format(datetime.datetime.now() - self.time))
        nodecolour(self, 0)


class No_Li_Gl(Node, ViNodes):
    '''Node describing a LiVi glare analysis'''
    bl_idname = 'No_Li_Gl'
    bl_label = 'LiVi Glare'
    bl_icon = 'IMAGE'

    def ret_params(self):
        return [str(x) for x in (self.hdrname, self.rand, self.gc)]

    def nodeupdate(self, context):
        self['images'] = [bpy.path.abspath(self.hdrfile)]
        nodecolour(self, self['exportstate'] != self.ret_params())

    hdrname: StringProperty(name="", description="Base name of the Glare image", default="", update=nodeupdate)
    hdrfile: StringProperty(name="", description="Location of the HDR image", default="", subtype="FILE_PATH", update=nodeupdate)
    vffile: StringProperty(name="", description="Location of the view file image", default="", subtype="FILE_PATH", update=nodeupdate)
    gc: FloatVectorProperty(size=3, name='', attr='Color', default=[1, 0, 0], subtype='COLOR', update=nodeupdate)
    rand: BoolProperty(name='', default=True, update=nodeupdate)
    x: IntProperty(name='', min=1, max=10000, default=2000, update=nodeupdate)
    y: IntProperty(name='', min=1, max=10000, default=1000, update=nodeupdate)
    month: IntProperty(name='', min=1, max=12, default=1, update=nodeupdate)
    day: IntProperty(name='', min=1, max=31, default=1, update=nodeupdate)
    hour: IntProperty(name='', min=1, max=24, default=1, update=nodeupdate)
    minutes: IntProperty(name='', min=0, max=59, default=0, update=nodeupdate)

    def init(self, context):
        self['exportstate'] = ''
        self.inputs.new('So_Li_Im', 'Image')

    def draw_buttons(self, context, layout):
        # These are here to provide ancilliary info for a regular HDR but Blender does not save HDRs with the required header information
        # if not self.inputs['Image'].links or not self.inputs['Image'].links[0].from_node['images'] or not os.path.isfile(bpy.path.abspath(self.inputs['Image'].links[0].from_node['images'][0])):
        #     row = layout.row()
        #     row.prop(self, 'hdrfile')
        #     newrow(layout, 'Month:', self, 'month')
        #     newrow(layout, 'Day:', self, 'day')
        #     newrow(layout, 'Hour:', self, 'hour')
        #     newrow(layout, 'Minutes:', self, 'minutes')
        #     newrow(layout, 'X dimen.:', self, 'x')
        #     newrow(layout, 'Y dimen.:', self, 'y')
        # if not sim.links or not sim.links[0].from_node['images'] or not os.path.isfile(bpy.path.abspath(sim.links[0].from_node['images'][0])):
        #     row = layout.row()
        #     row.prop(self, 'hdrfile')

        if not self.inputs['Image'].links or not self.inputs['Image'].links[0].from_node['images'] or not os.path.isfile(bpy.path.abspath(self.inputs['Image'].links[0].from_node['images'][0])):
            newrow(layout, 'Base image', self, 'hdrfile')

            if os.path.isfile(bpy.path.abspath(self.hdrfile)):
                newrow(layout, 'View image', self, 'vffile')

        if (self.inputs['Image'].links and self.inputs['Image'].links[0].from_node['images']
                and os.path.isfile(bpy.path.abspath(self.inputs['Image'].links[0].from_node['images'][0]))) or os.path.isfile(bpy.path.abspath(self.hdrfile)):
            newrow(layout, 'Base name:', self, 'hdrname')
            newrow(layout, 'Random:', self, 'rand')

            if not self.rand:
                newrow(layout, 'Colour:', self, 'gc')

            row = layout.row()
            row.operator("node.liviglare", text='Glare')

    def presim(self):
        self['hdrname'] = self.hdrname if self.hdrname else 'glare'

    def sim(self):
        for im in self['images']:
            with open(im+'glare', 'w') as glfile:
                Popen('evalglare {}'.format(im), stdout=glfile)

    def postsim(self):
        self['exportstate'] = self.ret_params()
        nodecolour(self, 0)


class No_Li_Fc(Node, ViNodes):
    '''Node describing a LiVi false colour image generation'''
    bl_idname = 'No_Li_Fc'
    bl_label = 'LiVi False Colour'
    bl_icon = 'IMAGE'

    def nodeupdate(self, context):
        nodecolour(self, self['exportstate'] != [str(x) for x in (self.basename, self.colour, self.lmax, self.unit, self.nscale, self.decades,
                   self.legend, self.lw, self.lh, self.contour, self.overlay, self.ofile, self.hdrfile)])
        self['images'] = [bpy.path.abspath(self.hdrfile)]

    basename: StringProperty(name="", description="Base name of the falsecolour image(s)", default="", update=nodeupdate)
    colour: EnumProperty(items=[("0", "Default", "Default color mapping"), ("1", "Spectral", "Spectral color mapping"), ("2", "Thermal", "Thermal colour mapping"),
                                ("3", "PM3D", "PM3D colour mapping"), ("4", "Eco", "Eco color mapping")], name="", description="Simulation accuracy", default="0", update=nodeupdate)
    lmax: IntProperty(name='', description="Legend max: 0 for auto", min=0, max=100000, default=1000, update=nodeupdate)
    unit: EnumProperty(items=[("0", "Lux", "Spectral color mapping"), ("1", "Candelas", "Thermal colour mapping"), ("2", "DF", "PM3D colour mapping"), ("3", "Irradiance(v)", "PM3D colour mapping")],
                       name="", description="Unit", default="0", update=nodeupdate)
    nscale: EnumProperty(items=[("0", "Linear", "Linear mapping"), ("1", "Log", "Logarithmic mapping")],
                         name="", description="Scale", default="0", update=nodeupdate)
    decades: IntProperty(name='', min=1, max=5, default=2, update=nodeupdate)
    unitdict = {'0': 'Lux', '1': 'cd/m2', '2': 'DF', '3': 'W/m2', '4': 'W/m2.sr'}
    unitmult = {'0': 179, '1': 179, '2': 1.79, '3': 1, '4': 1}
    legend: BoolProperty(name='', default=True, update=nodeupdate)
    lw: IntProperty(name='', min=1, max=1000, default=100, update=nodeupdate)
    lh: IntProperty(name='', min=1, max=1000, default=200, update=nodeupdate)
    contour: EnumProperty(items=[("0", "None", "No contours"), ("1", "Contour", "Line contours"), ("2", "Bands", "Banded contours"), ("3", "Poster", "Posterised contours")],
                          name="", description="Countour type", default="0", update=nodeupdate)
    overlay: BoolProperty(name='', default=False, update=nodeupdate)
    coldict = {'0': 'def', '1': 'spec', '2': 'hot', '3': 'pm3d', '4': 'eco'}
    divisions: IntProperty(name='', min=1, max=50, default=8, update=nodeupdate)
    ofile: StringProperty(name="", description="Location of the file to overlay", default="", subtype="FILE_PATH", update=nodeupdate)
    hdrfile: StringProperty(name="", description="Location of the file to overlay", default="", subtype="FILE_PATH", update=nodeupdate)
    unit_name: StringProperty(name="", description="Legend unit", default="", update=nodeupdate)
    multiplier: StringProperty(name="", description="Unit multiplier", default="", update=nodeupdate)
    disp: FloatProperty(name='', min=0.0001, max=10, default=1, precision=4, update=nodeupdate)

    def init(self, context):
        self['exportstate'] = ''
        self.inputs.new('So_Li_Im', 'Image')

    def draw_buttons(self, context, layout):
        sim = self.inputs['Image']

        if not sim.links or not sim.links[0].from_node['images'] or not os.path.isfile(bpy.path.abspath(sim.links[0].from_node['images'][0])):
            row = layout.row()
            row.prop(self, 'hdrfile')

        if (sim.links and sim.links[0].from_node['images'] and os.path.isfile(bpy.path.abspath(sim.links[0].from_node['images'][0]))) or os.path.isfile(bpy.path.abspath(self.hdrfile)):
            newrow(layout, 'Base name:', self, 'basename')
            newrow(layout, 'Unit:', self, 'unit_name')

            if self.unit_name:
                newrow(layout, 'Multiplier:', self, 'multiplier')
                newrow(layout, 'Colour:', self, 'colour')
                newrow(layout, 'Legend:', self, 'legend')

                if self.legend:
                    newrow(layout, 'Scale:', self, 'nscale')

                    if self.nscale == '1':
                        newrow(layout, 'Decades:', self, 'decades')

                    newrow(layout, 'Divisions:', self, 'divisions')
                    newrow(layout, 'Legend max:', self, 'lmax')
                    newrow(layout, 'Legend width:', self, 'lw')
                    newrow(layout, 'Legend height:', self, 'lh')

                newrow(layout, 'Contour:', self, 'contour')

                if self.contour != '0':
                    newrow(layout, 'Overlay:', self, 'overlay')

                    if self.overlay:
                        newrow(layout, 'Overlay file:', self, 'ofile')
                        newrow(layout, 'Overlay exposure:', self, 'disp')

                row = layout.row()
                row.operator("node.livifc", text='Process')

    def presim(self):
        self['basename'] = self.basename if self.basename else 'fc'

    def postsim(self):
        self['exportstate'] = [str(x) for x in (self.basename, self.colour, self.lmax, self.unit, self.nscale, self.decades,
                               self.legend, self.lw, self.lh, self.contour, self.overlay, self.ofile, self.hdrfile)]
        nodecolour(self, 0)

# Analysis nodes


class No_Li_Sim(Node, ViNodes):
    '''Node describing a LiVi simulation'''
    bl_idname = 'No_Li_Sim'
    bl_label = 'LiVi Simulation'

    def ret_params(self):
        return [str(x) for x in (self.cusacc, self.simacc, self.csimacc, self.pmap, self.pmapcno, self.pmapgno, self.pmapoptions, self.pmappreview)]

    def nodeupdate(self, context):
        nodecolour(self, self['exportstate'] != self.ret_params())

        if self.simacc == '3' or self.csimacc == "0":
            self.validparams = validradparams(self.cusacc)

    simacc: EnumProperty(items=[("0", "Low", "Low accuracy and high speed (preview)"), ("1", "Medium", "Medium speed and accuracy"), ("2", "High", "High but slow accuracy"),
                                ("3", "Custom", "Edit Radiance parameters")], name="", description="Simulation accuracy", default="0", update=nodeupdate)
    csimacc: EnumProperty(items=[("0", "Custom", "Edit Radiance parameters"), ("1", "Initial", "Initial accuracy for this metric"),
                                 ("2", "Final", "Final accuracy for this metric")], name="", description="Simulation accuracy", default="1", update=nodeupdate)
    cusacc: StringProperty(name="Custom parameters", description="Custom Radiance simulation parameters", default="", update=nodeupdate)

    pmap: BoolProperty(name='', description="Turn on photon mapping", default=False, update=nodeupdate)
    bfv: BoolProperty(name='', description="Turn on back face visibility (may cause light bleed but deals with planar geometry", default=True, update=nodeupdate)
    pmapgno: IntProperty(name='', default=50000, update=nodeupdate)
    pmapcno: IntProperty(name='', default=0, update=nodeupdate)
    pmapoptions: StringProperty(name="", description="Additional pmap parameters", default="", update=nodeupdate)
    pmappreview: BoolProperty(name='', default=0, update=nodeupdate)
    run: IntProperty(default=0)
    validparams: BoolProperty(name='', default=True)
    illu: BoolProperty(name='', default=False)
    camera: EnumProperty(items=ret_camera_menu, name='', description='Camera')
    new_res: BoolProperty(name='', default=False)

    def init(self, context):
        self['simdict'] = {'Basic': 'simacc', 'CBDM': 'csimacc'}
        self.inputs.new('So_Li_Geo', 'Geometry in')
        self.inputs.new('So_Li_Con', 'Context in')
        self.outputs.new('So_Vi_Res', 'Results out')
        nodecolour(self, 1)
        self['maxres'], self['minres'], self['avres'], self['exportstate'] = {}, {}, {}, ''

    def draw_buttons(self, context, layout):
        scene = context.scene
        svp = scene.vi_params

        if self.inputs['Context in'].links and self.inputs['Geometry in'].links:
            cinnode = self.inputs['Context in'].links[0].from_node
            ginnode = self.inputs['Geometry in'].links[0].from_node

            if cinnode.get('Options'):
                row = layout.row()
                row.label(text='Frames: {} - {}'.format(min([c['fs'] for c in (cinnode['Options'], ginnode['Options'])]), max([c['fe'] for c in (cinnode['Options'], ginnode['Options'])])))
                newrow(layout, 'Photon map:', self, 'pmap')

                if self.pmap:
                    newrow(layout, 'Global photons:', self, 'pmapgno')
                    newrow(layout, 'Caustic photons:', self, 'pmapcno')
                    newrow(layout, 'Back face visability:', self, 'bfv')
                    newrow(layout, 'Photon options:', self, 'pmapoptions')
                    newrow(layout, 'Preview photons:', self, 'pmappreview')

                row = layout.row()
                row.label(text="Accuracy:")
                row.prop(self, self['simdict'][cinnode['Options']['Context']])

                if (self.simacc == '3' and cinnode['Options']['Context'] == 'Basic') or (self.csimacc == '0' and cinnode['Options']['Context'] == 'CBDM'):
                    row = layout.row()
                    row.prop(self, 'cusacc')

                if not self.run and (self.simacc != '3' or self.validparams):
                    if cinnode['Options']['Preview']:
                        row = layout.row()
                        row.prop(self, "camera")

                        if self.camera != 'None':
                            row.operator("node.radpreview", text='Preview')

                    if [o for o in scene.objects if o.vi_params.vi_type_string == 'LiVi Calc']:
                        row = layout.row()
                        row.operator("node.livicalc", text='Calculate')

    def update(self):
        for sock in self.outputs:
            socklink(sock, self.id_data.name)
        self.run = 0

    def presim(self):
        self['coptions'] = self.inputs['Context in'].links[0].from_node['Options']
        self['goptions'] = self.inputs['Geometry in'].links[0].from_node['Options']
        self['radfiles'], self['reslists'] = {}, [[]]

        if self['coptions']['Context'] == 'Basic':
            self['radparams'] = ' {} '.format(self.cusacc) if self.simacc == '3' else ''.join([' {} {} '.format(k, rtraceparams[k][int(self.simacc)]) for k in rtraceparams])
            self['rvuparams'] = ' {} '.format(self.cusacc) if self.simacc == '3' else ''.join([' {} {} '.format(k, rvuparams[k][int(self.simacc)]) for k in rvuparams])
        else:
            self['radparams'] = ' {} '.format(self.cusacc) if self.csimacc == '0' else ''.join([' {} {} '.format(k, rtracecbdmparams[k][int(self.simacc) - 1]) for k in rtracecbdmparams])
            self['rvuparams'] = ' {} '.format(self.cusacc) if self.csimacc == '0' else ''.join([' {} {} '.format(k, rvuparams[k][int(self.simacc) - 1]) for k in rvuparams])

    def sim(self, scene):
        svp = scene.vi_params
        self['frames'] = range(svp['liparams']['fs'], svp['liparams']['fe'] + 1)

    def postsim(self, calcout):
        self['exportstate'] = self.ret_params()
        self['reslists'] = calcout

        if self.outputs[0].links:
            for dnode in set([li.to_node for li in self.outputs[0].links]):
                if dnode.bl_idname == 'No_Vi_Metrics':
                    dnode.update()
                elif dnode.bl_idname == 'No_Vi_HMChart':
                    dnode.update()

        nodecolour(self, 0)


class No_Vi_SP(Node, ViNodes):
    '''Node describing a VI-Suite sun path'''
    bl_idname = 'No_Vi_SP'
    bl_label = 'VI Sun Path'

    def nodeupdate(self, context):
        nodecolour(self, self['exportstate'] != [str(x) for x in (self.suns)])

    suns: EnumProperty(items=[('0', 'Single', 'Single sun'), ('1', 'Monthly', 'Monthly sun for chosen time'), ('2', 'Hourly', 'Hourly sun for chosen date')],
                       name='', description='Sunpath sun type', default='0', update=nodeupdate)

    def init(self, context):
        self.inputs.new('So_Vi_Loc', 'Location in')
        self['exportstate'] = '0'
        nodecolour(self, 1)

    def draw_buttons(self, context, layout):
        if self.inputs['Location in'].links:
            newrow(layout, 'Suns:', self, 'suns')
            row = layout.row()
            row.operator("node.sunpath", text="Create Sun Path")
        else:
            row = layout.row()
            row.label(text="Connect location node")

    def export(self):
        nodecolour(self, 0)
        self['exportstate'] = [str(x) for x in (self.suns)]

    def update(self):
        pass


class No_Vi_WR(Node, ViNodes):
    '''Node describing a VI-Suite wind rose generator'''
    bl_idname = 'No_Vi_WR'
    bl_label = 'VI Wind Rose'
    bl_icon = 'FORCE_WIND'

    def nodeupdate(self, context):
        nodecolour(self, self['exportstate'] != [str(x) for x in (self.wrtype, self.sdoy, self.edoy, self.max_freq, self.max_freq_val, self.temp)])

    wrtype: EnumProperty(items=[("0", "Hist 1", "Stacked histogram"), ("1", "Hist 2", "Stacked Histogram 2"), ("2", "Cont 1", "Filled contour"),
                                ("3", "Cont 2", "Edged contour"), ("4", "Cont 3", "Lined contour")], name="", default='0', update=nodeupdate)
    sdoy: IntProperty(name="", description="Day of simulation", min=1, max=365, default=1, update=nodeupdate)
    edoy: IntProperty(name="", description="Day of simulation", min=1, max=365, default=365, update=nodeupdate)
    max_freq: EnumProperty(items=[("0", "Data", "Max frequency taken from data"), ("1", "Specified", "User entered value")], name="", default='0', update=nodeupdate)
    max_freq_val: FloatProperty(name="", description="Max frequency", min=1, max=100, default=20, update=nodeupdate)
    temp: BoolProperty(name="", description="Plot temperature", default=0, update=nodeupdate)

    def init(self, context):
        self.inputs.new('So_Vi_Loc', 'Location in')
        self['exportstate'] = ''
        nodecolour(self, 1)

    def draw_buttons(self, context, layout):
        if nodeinputs(self) and self.inputs[0].links[0].from_node.loc == '1':
            (sdate, edate) = retdates(self.sdoy, self.edoy, context.scene.vi_params.year)
            newrow(layout, 'Type:', self, "wrtype")
            newrow(layout, 'Temperature:', self, "temp")
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
            row.label(text='Location node error')

    def export(self):
        nodecolour(self, 0)
        self['exportstate'] = [str(x) for x in (self.wrtype, self.sdoy, self.edoy, self.max_freq, self.max_freq_val, self.temp)]

    def update(self):
        pass


class No_Vi_SVF(Node, ViNodes):
    '''Node for sky view factor analysis'''
    bl_idname = 'No_Vi_SVF'
    bl_label = 'VI SVF'
    bl_icon = 'MOD_SOFT'

    def nodeupdate(self, context):
        context.scene.vi_params.vi_nodes = self.id_data
        nodecolour(self, self['exportstate'] != [str(x) for x in (self.startframe, self.endframe, self.cpoint, self.offset, self.animmenu)])

    animtype = [('Static', "Static", "Simple static analysis"), ('Geometry', "Geometry", "Animated geometry analysis")]
    animmenu: EnumProperty(name="", description="Animation type", items=animtype, default='Static', update=nodeupdate)
    startframe: IntProperty(name='', default=0, min=0, max=1024, description='Start frame')
    endframe: IntProperty(name='', default=0, min=0, max=1024, description='End frame')
    cpoint: EnumProperty(items=[("0", "Faces", "Export faces for calculation points"), ("1", "Vertices", "Export vertices for calculation points")],
                         name="", description="Specify the calculation point geometry", default="0", update=nodeupdate)
    offset: FloatProperty(name="", description="Calc point offset", min=0.001, max=10, default=0.01, precision=3, update=nodeupdate)
    signore: BoolProperty(name='', default=0, description='Ignore sensor surfaces', update=nodeupdate)
    skytype = [('0', "Tregenza", "145 Tregenza sky patches"), ('1', "Reinhart 577", "577 Reinhart sky patches"), ('2', 'Reinhart 2305', '2305 Reinhart sky patches')]
    skypatches: EnumProperty(name="", description="Animation type", items=skytype, default='0', update=nodeupdate)

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
            row.label(text='Frames:')
            row.prop(self, 'startframe')
            row.prop(self, 'endframe')

        newrow(layout, 'Sky patches:', self, "skypatches")
        newrow(layout, 'Result point:', self, "cpoint")
        newrow(layout, 'Offset:', self, 'offset')
        row = layout.row()
        row.operator("node.svf", text="Sky View Factor")

    def preexport(self):
        self['goptions']['offset'] = self.offset

    def postexport(self, scene):
        nodecolour(self, 0)
        self.outputs['Results out'].hide = False if self.get('reslists') else True
        self['exportstate'] = [str(x) for x in (self.startframe, self.endframe, self.cpoint, self.offset, self.animmenu)]

    def update(self):
        for sock in self.outputs:
            socklink(sock, self.id_data.name)


class No_Vi_SS(Node, ViNodes):
    '''Node to create a VI-Suite shadow map'''
    bl_idname = 'No_Vi_SS'
    bl_label = 'VI Shadow Map'

    def nodeupdate(self, context):
        context.scene.vi_params.vi_nodes = self.id_data
        nodecolour(self, self['exportstate'] != [str(x) for x in (self.animmenu, self.sdoy, self.edoy, self.starthour, self.endhour, self.interval, self.cpoint, self.offset)])

    animtype = [('Static', "Static", "Simple static analysis"), ('Geometry', "Geometry", "Animated geometry analysis")]
    animmenu: EnumProperty(name="", description="Animation type", items=animtype, default='Static', update=nodeupdate)
    startframe: IntProperty(name='', default=0, min=0, max=1024, description='Start frame')
    endframe: IntProperty(name='', default=0, min=0, max=1024, description='End frame')
    starthour: IntProperty(name='', default=1, min=1, max=24, description='Start hour')
    endhour: IntProperty(name='', default=24, min=1, max=24, description='End hour')
    interval: IntProperty(name='', default=1, min=1, max=60, description='Number of simulation steps per hour')
    sdoy: IntProperty(name='', default=1, min=1, max=365, description='Start Day', update=nodeupdate)
    edoy: IntProperty(name='', default=365, min=1, max=365, description='End Day', update=nodeupdate)
    cpoint: EnumProperty(items=[("0", "Faces", "Export faces for calculation points"), ("1", "Vertices", "Export vertices for calculation points")],
                         name="", description="Specify the calculation point geometry", default="0", update=nodeupdate)
    offset: FloatProperty(name="", description="Calc point offset", min=0.001, max=10, default=0.01, update=nodeupdate)
    signore: BoolProperty(name='', default=0, description='Ignore sensor surfaces', update=nodeupdate)

    def init(self, context):
        self.inputs.new('So_Vi_Loc', 'Location in')
        self.outputs.new('So_Vi_Res', 'Results out')
        self.outputs['Results out'].hide = True
        self['exportstate'] = ''
        self['goptions'] = {}
        nodecolour(self, 1)

    def draw_buttons(self, context, layout):
        if nodeinputs(self):
            (sdate, edate) = retdates(self.sdoy, self.edoy, context.scene.vi_params.year)
            newrow(layout, 'Ignore sensor:', self, "signore")
            newrow(layout, 'Animation:', self, "animmenu")

            if self.animmenu != 'Static':
                row = layout.row(align=True)
                row.alignment = 'EXPAND'
                row.label(text='Frames:')
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
            row.operator("node.shad", text='Calculate')

    def preexport(self):
        (self.sdate, self.edate) = retdates(self.sdoy, self.edoy, bpy.context.scene.vi_params.year)
        self['goptions']['offset'] = self.offset

    def postexport(self, scene):
        nodecolour(self, 0)
        self.outputs['Results out'].hide = False if self.get('reslists') else True
        self['exportstate'] = [str(x) for x in (self.animmenu, self.sdoy, self.edoy, self.starthour, self.endhour, self.interval, self.cpoint, self.offset)]

    def update(self):
        for sock in self.outputs:
            socklink(sock, self.id_data.name)


class No_Vi_EC(Node, ViNodes):
    '''Node to calculate embodied carbon'''
    bl_idname = 'No_Vi_EC'
    bl_label = 'VI Embodied Carbon'

    def ret_params(self):
        return [str(x) for x in (self.entities, self.parametric, self.startframe, self.endframe, self.tyears, self.heal)]

    def nodeupdate(self, context):
        nodecolour(self, self['exportstate'] != self.ret_params())

    entities: EnumProperty(name="", description="Entity type", items=[('0', "Objects", "Calculate EC based on objects"), ('1', "Collections", "Calculate EC based on collections")], default='0', update=nodeupdate)
    parametric: BoolProperty(name='', default=0, description='Ignore sensor surfaces', update=nodeupdate)
    heal: BoolProperty(name='', default=0, description='Attempt to heal meshes', update=nodeupdate)
    fa: FloatProperty(name='m2', default=50, min=1, max=1024, description='Floor area', update=nodeupdate)
    startframe: IntProperty(name='', default=0, min=0, max=1024, description='Start frame', update=nodeupdate)
    endframe: IntProperty(name='', default=0, min=0, max=1024, description='End frame', update=nodeupdate)
    tyears: IntProperty(name='Years', default=60, min=1, max=120, description='Total timeframe', update=nodeupdate)

    def init(self, context):
        self.outputs.new('So_Vi_Res', 'Results out')
        self.outputs['Results out'].hide = True
        self['exportstate'] = ''
        nodecolour(self, 1)

    def draw_buttons(self, context, layout):
        newrow(layout, 'Entities:', self, 'entities')
        newrow(layout, 'Parametric:', self, 'parametric')

        if self.parametric:
            row = layout.row(align=True)
            row.alignment = 'EXPAND'
            row.label(text='Frames:')
            row.prop(self, 'startframe')
            row.prop(self, 'endframe')

        newrow(layout, 'Heal:', self, 'heal')
        newrow(layout, 'Floor area:', self, 'fa')
        newrow(layout, 'Timeframe:', self, 'tyears')
        row = layout.row()
        row.operator("node.ec_calc", text='Calculate')

    def presim(self, context):
        self['frames'] = [context.scene.frame_current] if not self.parametric else range(self.startframe, self.endframe + 1)

    def postsim(self):
        self.outputs['Results out'].hide = False
        nodecolour(self, 0)
        self['exportstate'] = self.ret_params()

        if self.outputs[0].links:
            for dnode in set([link.to_node for link in self.outputs[0].links]):
                if dnode.bl_idname == 'No_Vi_Metrics':
                    dnode.update()


# Edit nodes


class No_Text(Node, ViNodes):
    '''Text Export Node'''
    bl_idname = 'No_Text'
    bl_label = 'VI Text Edit'

    contextmenu: StringProperty(name='')

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
            row.label(text='Text name: {}'.format(inodename))
            if inodename in [im.name for im in bpy.data.texts] and self['bt'] != bpy.data.texts[inodename].as_string():
                row = layout.row()
                row.operator('node.textupdate', text='Update')

    def update(self):
        for sock in self.outputs:
            socklink(sock, self.id_data.name)

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
        self['Text'] = {bthb[0]: bthb[1] for bthb in zip(btframes, btbodies)}


class No_En_Geo(Node, ViNodes):
    '''Node describing an EnVi Geometry Export'''
    bl_idname = 'No_En_Geo'
    bl_label = 'EnVi Geometry'

    geo_offset: FloatVectorProperty(name="", description="", default=(0.0, 0.0, 0.0), min=sys.float_info.min,
                                    max=sys.float_info.max, soft_min=sys.float_info.min, soft_max=sys.float_info.max, step=3, precision=1,
                                    subtype='TRANSLATION', unit='NONE', size=3, update=None, get=None, set=None)

    def init(self, context):
        self.outputs.new('So_En_Geo', 'Geometry out')
        nodecolour(self, 1)

    def draw_buttons(self, context, layout):
        newrow(layout, 'Offset', self, 'geo_offset')
        row = layout.row()
        row.operator("node.engexport", text="Export")

    def update(self):
        for sock in self.outputs:
            socklink(sock, self.id_data.name)

    def preexport(self, scene):
        pass

    def postexport(self):
        nodecolour(self, 0)


class No_En_Con(Node, ViNodes):
    '''Node describing an EnergyPlus export'''
    bl_idname = 'No_En_Con'
    bl_label = 'EnVi Export'

    def nodeupdate(self, context):
        context.scene.vi_params.vi_nodes = self.id_data
        nodecolour(self, self['exportstate'] != [str(x) for x in (self.loc, self.terrain, self.timesteps, self.animated, self.fs, self.fe, self.sdoy, self.edoy)])

    animated: BoolProperty(name="", description="Animated analysis", update=nodeupdate)
    fs: IntProperty(name="", description="Start frame", default=0, min=0, update=nodeupdate)
    fe: IntProperty(name="", description="End frame", default=0, min=0, update=nodeupdate)
    loc: StringProperty(name="", description="Identifier for this project", default="", update=nodeupdate)
    terrain: EnumProperty(items=[("0", "City", "Towns, city outskirts, centre of large cities"),
                                 ("1", "Urban", "Urban, Industrial, Forest"), ("2", "Suburbs", "Rough, Wooded Country, Suburbs"),
                                 ("3", "Country", "Flat, Open Country"), ("4", "Ocean", "Ocean, very flat country")],
                          name="", description="Exposure context", default="0", options={'SKIP_SAVE'}, update=nodeupdate)
    solar: EnumProperty(items=[("0", "MinimalShadowing", ""), ("1", "FullExterior", ""), ("2", "FullInteriorAndExterior", ""),
                               ("3", "FullExteriorWithReflections", ""), ("4", "FullInteriorAndExteriorWithReflections", "")],
                        name="", description="Solar calcs", default="2", update=nodeupdate)
    addonpath = os.path.dirname(inspect.getfile(inspect.currentframe()))
    matpath = addonpath+'/EPFiles/Materials/Materials.data'
    sdoy: IntProperty(name="", description="Start day of simulation", min=1, max=365, default=1, update=nodeupdate)
    edoy: IntProperty(name="", description="End day of simulation", min=1, max=365, default=365, update=nodeupdate)
    timesteps: IntProperty(name="", description="Time steps per hour", min=1, max=60, default=4, update=nodeupdate)
    shadow_calc: EnumProperty(items=[("0", "CPU", "CPU based shadow calculations"), ("1", "GPU", "GPU based shadow calculations")],
                              name="", description="Specify the EnVi results category", default="0", update=nodeupdate)
    restype: EnumProperty(items=[("0", "Zone Thermal", "Thermal Results"), ("1", "Comfort", "Comfort Results"),
                                 ("2", "Zone Ventilation", "Zone Ventilation Results"), ("3", "Ventilation Link", "Ventilation Link Results"),
                                 ("4", "Power", "Power Production Results")],
                          name="", description="Specify the EnVi results category", default="0", update=nodeupdate)

    resaam: bpy.props.BoolProperty(name='Air', description='Ambient air metrics', default=False)
    resaws: bpy.props.BoolProperty(name='Wind speed', description='Ambient wind speed', default=False)
    resawd: bpy.props.BoolProperty(name='Wind direction', description='Ambient wind direction', default=False)
    resah: bpy.props.BoolProperty(name='Humidity', description='Ambient humidity', default=False)
    resasm: bpy.props.BoolProperty(name='Solar', description='Ambient solar metrics', default=False)
    restt: bpy.props.BoolProperty(name='Temperature', description='Zone temperature (degC)', default=False)
    resh: bpy.props.BoolProperty(name='Humidity', description='Zone humidity (%)', default=False)
    restwh: bpy.props.BoolProperty(name='Heating', description='Zone heating (W)', default=False)
    restwc: bpy.props.BoolProperty(name='Cooling', description='Zone cooling (W)', default=False)
    reswsg: bpy.props.BoolProperty(name='Solar gain', description='Window solar gain (W)', default=False)
    rescpp: bpy.props.BoolProperty(name='PPD', description='Percentage People Dissatisfied', default=False)
    rescpm: bpy.props.BoolProperty(name='PMV', description='Predicted Mean Vote', default=False)
    resvls: bpy.props.BoolProperty(name='Ventilation (l/s)', description='Zone ventilation rate (l/s)', default=False)
    resvmh: bpy.props.BoolProperty(name='Ventilation (m\u00b3/h)', description='Zone ventilation rate (m\u00b3/h)', default=False)
    resim: bpy.props.BoolProperty(name='Infiltration (m\u00b3/h)', description='Zone infiltration rate (m\u00b3/h)', default=False)
    resiach: bpy.props.BoolProperty(name='Infiltration (ACH)', description='Zone infiltration rate (ACH)', default=False)
    resco2: bpy.props.BoolProperty(name='CO2 (ppm)', description='Zone CO2 concentration (ppm)', default=False)
    resihl: bpy.props.BoolProperty(name='Exfiltration (W)', description='Exfiltration heat loss (W)', default=False)
    resl12ms: bpy.props.BoolProperty(name=u'Flow (m\u00b3/s)', description=u'Linkage flow (m\u00b3/s)', default=False)
    reslof: bpy.props.BoolProperty(name='Opening factor', description='Linkage opening factor', default=False)
    resmrt: bpy.props.BoolProperty(name='MRT', description='Mean radiant temperature (degC)', default=False)
    resocc: bpy.props.BoolProperty(name='Occupancy', description='Zone occupancy', default=False)
    resh: bpy.props.BoolProperty(name='Humidity', description='Zone humidity (%)', default=False)
    resfhb: bpy.props.BoolProperty(name='Heat balance', description='Fabric heat balance (W)', default=False)
    ressah: bpy.props.BoolProperty(name='Air heating', description='Air heating (W)', default=False)
    ressac: bpy.props.BoolProperty(name='Air cooling', description='Air cooling (W)', default=False)
    reshrhw: bpy.props.BoolProperty(name='Heat recovery', description='Heat recovery heating (W)', default=False)
    resldp: bpy.props.BoolProperty(name='Delta P', description='Pressure difference (Pa)', default=False)
    resoeg: bpy.props.BoolProperty(name='Equipment', description='Other equipment heat gains (W)', default=False)
    respve: bpy.props.BoolProperty(name='PV Energy', description='PV energy (J)', default=False)
    respvw: bpy.props.BoolProperty(name='PV Power', description='PV power (W)', default=False)
    respvt: bpy.props.BoolProperty(name='PV Temperature', description='PV Temperature (degC)', default=False)
    respveff: bpy.props.BoolProperty(name='PV Efficiency', description='PV efficiency (%)', default=False)

    def init(self, context):
        self.inputs.new('So_En_Geo', 'Geometry in')
        self.inputs.new('So_Vi_Loc', 'Location in')
        self.outputs.new('So_En_Con', 'Context out')
        self.outputs.new('So_Anim', 'Parameter')
        self['exportstate'] = ''
        nodecolour(self, 1)

    def draw_buttons(self, context, layout):
        (sdate, edate) = retdates(self.sdoy, self.edoy, context.scene.vi_params.year)
        newrow(layout, 'Shadow:', self, 'shadow_calc')
        row = layout.row()
        row.label(text='Animation:')
        row.prop(self, 'animated')

        if self.animated:
            newrow(layout, 'Start frame:', self, 'fs')
            newrow(layout, 'End frame:', self, 'fe')

        newrow(layout, "Name/location:", self, "loc")
        newrow(layout, "Solar:", self, "solar")
        row = layout.row()
        row.label(text = 'Exposure:')
        col = row.column()
        col.prop(self, "terrain")
        newrow(layout, 'Start day {}/{}:'.format(sdate.day, sdate.month), self, "sdoy")
        newrow(layout, 'End day {}/{}:'.format(edate.day, edate.month), self, "edoy")
        newrow(layout, 'Time-steps/hour', self, "timesteps")
        row = layout.row()
        row.label(text='Results Category:')
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
            row.operator("node.encon", text='Export')

    def update(self):
        for sock in self.outputs:
            socklink(sock, self.id_data.name)

    def preexport(self, scene):
        (self.fs, self.fe) = (self.fs, self.fe) if self.animated else (scene.frame_current, scene.frame_current)
        scene.vi_params['enparams']['fs'], scene.vi_params['enparams']['fe'] = self.fs, self.fe
        (self.sdate, self.edate) = retdates(self.sdoy, self.edoy, scene.vi_params.year)

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
    resfilename: StringProperty(name="", default='results')
    dsdoy: IntProperty()
    dedoy: IntProperty()
    run: IntProperty(min=-1, default=-1)
    processors: IntProperty(name='', min=1, default=4)  #max = bpy.context.scene['viparams']['nproc'], default = bpy.context.scene['viparams']['nproc'])
    mp: BoolProperty(name="", default=False)
    # finished: BoolProperty(name="", default=False)

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
        for sock in self.outputs:
            socklink(sock, self.id_data.name)

        self.run = -1

    def presim(self, context):
        innode = self.inputs['Context in'].links[0].from_node
        self['frames'] = range(context.scene.vi_params['enparams']['fs'], context.scene.vi_params['enparams']['fe'] + 1)
        self.resfilename = os.path.join(context.scene.vi_params['viparams']['newdir'], self.resname+'.eso')
        self.dsdoy = innode.sdoy
        self.dedoy = innode.edoy
        self["_RNA_UI"] = {"Start": {"min":innode.sdoy, "max":innode.edoy}, "End": {"min":innode.sdoy, "max":innode.edoy},
            "AStart": {"name": '', "min": context.scene.vi_params['enparams']['fs'], "max": context.scene.vi_params['enparams']['fe']},
            "AEnd": {"min":context.scene.vi_params['enparams']['fs'], "max":context.scene.vi_params['enparams']['fe']}}
        self['Start'], self['End'] = innode.sdoy, innode.edoy
        self['AStart'], self['AEnd'] = context.scene.vi_params['enparams']['fs'], context.scene.vi_params['enparams']['fe']

    def postsim(self, sim_op, condition):
        innode = self.inputs['Context in'].links[0].from_node
        nodecolour(self, 0)
        self.run = -1

        if condition == 'FINISHED':
            processf(sim_op, self, innode)

            if self.outputs[0].links:
                for dnode in set([l.to_node for l in self.outputs[0].links]):
                    if dnode.bl_idname in ('No_Vi_Metrics', 'No_Vi_HMChart', 'No_Vi_Chart'):
                        dnode.update()


class No_En_IF(Node, ViNodes):
    '''Node for EnergyPlus input file selection'''
    bl_idname = 'ViEnInNode'
    bl_label = 'EnVi Input File'

    def nodeupdate(self, context):
        context.scene['enparams']['fs'] = context.scene['enparams']['fe'] = context.scene.frame_current
        shutil.copyfile(self.idffilename, os.path.join(context.scene['viparams']['newdir'], 'in{}.idf'.format(context.scene.frame_current)))
        locnode = self.inputs['Location in'].links[0].from_node
        shutil.copyfile(locnode.weather, os.path.join(context.scene['viparams']['newdir'], 'in{}.epw'.format(context.scene.frame_current)))

        with open(self.idffilename, 'r', errors='ignore') as idff:
            idfflines = idff.readlines()
            for l, line in enumerate(idfflines):
                if line.split(',')[0].lstrip(' ').upper() == 'RUNPERIOD':
                    self.sdoy = datetime.datetime(context.scene.vi_params.year, int(idfflines[l+2].split(',')[0].lstrip(' ')), int(idfflines[l+3].split(',')[0].lstrip(' '))).timetuple().tm_yday
                    self.edoy = datetime.datetime(context.scene.vi_params.year, int(idfflines[l+4].split(',')[0].lstrip(' ')), int(idfflines[l+5].split(',')[0].lstrip(' '))).timetuple().tm_yday
                    self.outputs['Context out'].hide = False
                    break
            nodecolour(self, 0)

    idfname = StringProperty(name="", description="Name of the EnVi results file", default="", update=nodeupdate)
    sdoy = IntProperty(name='', default=1, min=1, max=365)
    edoy = IntProperty(name='', default=365, min=1, max=365)
    newdir = StringProperty()

    def init(self, context):
        self.inputs.new('ViLoc', 'Location in')
        self.outputs.new('ViEnC', 'Context out')
        self.outputs['Context out'].hide = True
        self['nodeid'] = nodeid(self)
        nodecolour(self, 1)

    def draw_buttons(self, context, layout):
        row = layout.row()

        if self.inputs['Location in'].links:
            row = layout.row()
            row.operator('node.idfselect', text='IDF select').nodeid = self['nodeid']
            row.prop(self, 'idfname')
        else:
            row.label('Connect Location node')

    def update(self):
        if self.outputs.get('Context out'):
            (self.outputs['Context out'], self['nodeid'].split('@')[1])
        if not self.inputs['Location in'].links:
            nodecolour(self, 1)


class No_En_RF(Node, ViNodes):
    '''Node for EnergyPlus results file selection'''
    bl_idname = 'ViEnRFNode'
    bl_label = 'EnVi Results File'

    def nodeupdate(self, context):
        nodecolour(self, self['exportstate'] != [self.resfilename])
        self['frames'] = [context.scene.frame_current]

    esoname = StringProperty(name="", description="Name of the EnVi results file", default="", update=nodeupdate)
    filebase = StringProperty(name="", description="Name of the EnVi results file", default="")
    dsdoy, dedoy = IntProperty(), IntProperty()

    def init(self, context):
        self.outputs.new('ViR', 'Results out')
        self['exportstate'] = ''
        nodecolour(self, 1)

    def draw_buttons(self, context, layout):
        row = layout.row()
        row.operator('node.esoselect', text='ESO select').nodeid = self['nodeid']
        row.prop(self, 'esoname')
        row = layout.row()
        row.operator("node.fileprocess", text='Process file').nodeid = self['nodeid']

    def update(self):
        for sock in self.outputs:
            socklink(sock, self.id_data.name)

    def export(self):
        self['exportstate'] = [self.resfilename]
        nodecolour(self, 0)


class So_En_ResU(NodeSocket):
    '''Vi unlinked results socket'''
    bl_idname = 'So_En_ResU'
    bl_label = 'Axis'
    valid = ['Vi Results']

    r_len: IntProperty()

    def draw_color(self, context, node):
        return (0.0, 1.0, 0.0, 0.75)

    def draw(self, context, layout, node, text):
        layout.label(text=self.bl_label)


class ViEnRIn(So_En_ResU):
    '''Chart input socket'''
    bl_idname = 'ViEnRIn'
    bl_label = 'Axis'

    def f_menu(self, context):
        try:
            frs = [(f'{frame}', f'{frame}', f'Frame {frame}') for frame in self['resdict'].keys() if frame]

            if frs:
                return frs
            else:
                return [('None', 'None', 'None')]
        except:
            return [('None', 'None', 'None')]

    def r_menu(self, context):
        try:
            ress = [(f'{res}', f'{res}', f'Frame {res}') for res in self['resdict'][self.framemenu].keys() if res]
            if ress:
                return ress
            else:
                return [('None', 'None', 'None')]
        except Exception:
            return [('None', 'None', 'None')]

    def z_menu(self, context):
        try:
            zres = [(f'{res}', f'{res}', f'Zone {res}') for res in self['resdict'][self.framemenu][self.resultmenu].keys() if res]

            if zres:
                return zres
            else:
                return [('None', 'None', 'None')]
        except Exception:
            try:
                r = list(self['resdict'][self.framemenu].keys())[0]
                z = list(self['resdict'][self.framemenu][r].keys())[0]
                self.zonemenu = list(self['resdict'][self.framemenu][r].keys())[0]
                return [(f'{res}', f'{res}', f'Zone {res}') for res in self['resdict'][self.framemenu][r].keys() if res]
            except Exception:
                return [('None', 'None', 'None')]

    def m_menu(self, context):
        try:
            m_entries =  [(f'{res}', f'{res}', f'Frame {res}') for res in self['resdict'][self.framemenu][self.resultmenu][self.zonemenu] if res]

            if m_entries:
                return m_entries
            else:
                return [('None', 'None', 'None')]
        except Exception:
            # try:
            #     r = list(self['resdict'][self.framemenu].keys())[0]
            #     z = list(self['resdict'][self.framemenu][r].keys())[0]
            #     return [(f'{res}', f'{res}', f'Frame {res}') for res in self['resdict'][self.framemenu][r][z] if res]
            # except Exception:
            return [('None', 'None', 'None')]

    def f_update(self, context):
        if not self.framemenu:
            self.framemenu = self.f_menu(context)[0][0]

        self.r_update(context)

    def r_update(self, context):
        try:
            innode = self.links[0].from_node
            rl = innode['reslists']
            zrl = list(zip(*rl))
            self.r_len = [len(res[4].split()) for res in rl if res[0] == self.framemenu and res[1] == self.resultmenu and res[2] == self.zonemenu and res[3] == self.metricmenu][0]
        except Exception:
            self.r_len = 0

        if not self.resultmenu:
            self.resultmenu = self.r_menu(context)[0][0]

        self.z_update(context)

        if self.name == 'X-axis':
            if self.resultmenu == 'Time':
                startday = datetime.datetime(bpy.context.scene.vi_params.year, int(zrl[4][zrl[3].index('Month')].split()[0]), int(zrl[4][zrl[3].index('Day')].split()[0])).timetuple().tm_yday
                endday = datetime.datetime(bpy.context.scene.vi_params.year, int(zrl[4][zrl[3].index('Month')].split()[-1]), int(zrl[4][zrl[3].index('Day')].split()[-1])).timetuple().tm_yday
                self.node["_RNA_UI"] = {"Start": {"min":startday, "max":endday}, "End": {"min":startday, "max":endday}}
                self.node['Start'], self.node['End'] = startday, endday

            elif self.resultmenu == 'Frames':
                frames = [int(k) for k in set(zrl[0]) if k != 'All']
                startframe, endframe = min(frames), max(frames)
                frame = 'All'
                self.node["_RNA_UI"] = {"Start": {"min":startframe, "max":endframe}, "End": {"min":startframe, "max":endframe}}
                self.node['Start'], self.node['End'] = startframe, endframe

            else:
                xs = range(1, 1 + self.r_len)
                startx, endx = min(xs), max(xs)
                self.node["_RNA_UI"] = {"Start": {"min":startx, "max":endx}, "End": {"min":startx, "max":endx}}
                self.node['Start'], self.node['End'] = startx, endx

    def z_update(self, context):
        try:
            innode = self.links[0].from_node
            rl = innode['reslists']
            self.r_len = [len(res[4].split()) for res in rl if res[0] == self.framemenu and res[1] == self.resultmenu and res[2] == self.zonemenu and res[3] == self.metricmenu][0]
        except Exception:
            self.r_len = 0

        if not self.zonemenu:
            self.zonemenu = self.z_menu(context)[0][0]

        self.m_update(context)

    def m_update(self, context):
        try:
            innode = self.links[0].from_node
            rl = innode['reslists']
            self.r_len = [len(res[4].split()) for res in rl if res[0] == self.framemenu and res[1] == self.resultmenu and res[2] == self.zonemenu and res[3] == self.metricmenu][0]
        except Exception:
            self.r_len = 0

        if self.name == 'X-axis':
            if not self.metricmenu:
                self.metricmenu = self.m_menu(context)[0][0]

            self.update_menus(('Y-axis 1', 'Y-axis 2', 'Y-axis 3'))
        else:
            if not self.metricmenu:
                self.metricmenu = self.m_menu(context)[0][0]

    def update_menus(self, axes):
        rl = self.node.inputs['X-axis'].links[0].from_node['reslists']

        if rl:
            try:
                self.r_len = [len(res[4].split()) for res in rl if res[0] == self.node.inputs['X-axis'].framemenu and res[1] == self.node.inputs['X-axis'].resultmenu and res[2] == self.node.inputs['X-axis'].zonemenu and res[3] == self.node.inputs['X-axis'].metricmenu][0]
            except:
                self.r_len = 0

            for ax in axes:
                if self.node.inputs[ax].links:
                    rl = self.node.inputs[ax].links[0].from_node['reslists']
                    resdict = {}

                    for r in rl:
                        if r and len(r[4].split()) == self.node.inputs['X-axis'].r_len and r[1] not in ('Time', 'Frames'):
                            resdict[r[0]] = {} if r[0] not in resdict else resdict[r[0]]

                            try:
                                resdict[r[0]][r[1]] = {} if r[1] not in resdict[r[0]] else resdict[r[0]][r[1]]
                                resdict[r[0]][r[1]][r[2]] = [] if r[2] not in resdict[r[0]][r[1]] else resdict[r[0]][r[1]][r[2]]
                                resdict[r[0]][r[1]][r[2]].append(r[3])
                            except Exception:
                                pass

                    self.node.inputs[ax]['resdict'] = resdict
                    self.node.inputs[ax].r_len = self.node.inputs['X-axis'].r_len if resdict else 0

            # self.update_entries(axes)

    def update_entries(self, context):
        axes = ('X-axis', 'Y-axis 1', 'Y-axis 2', 'Y-axis 3')

        for ax in axes:
            if self.node.inputs[ax].framemenu not in [m[0] for m in self.node.inputs[ax].f_menu(context)]:
                self.node.inputs[ax].framemenu = [m[0] for m in self.node.inputs[ax].f_menu(context)][0]
            if self.node.inputs[ax].resultmenu not in [m[0] for m in self.node.inputs[ax].r_menu(context)]:
                self.node.inputs[ax].resultmenu = [m[0] for m in self.node.inputs[ax].r_menu(context)][0]
            if self.node.inputs[ax].zonemenu not in [m[0] for m in self.node.inputs[ax].z_menu(context)]:
                self.node.inputs[ax].zonemenu = [m[0] for m in self.node.inputs[ax].z_menu(context)][0]
            if self.node.inputs[ax].metricmenu not in [m[0] for m in self.node.inputs[ax].m_menu(context)]:
                self.node.inputs[ax].metricmenu = [m[0] for m in self.node.inputs[ax].m_menu(context)][0]


    framemenu: EnumProperty(items=f_menu, name="", description="Frame number", update=f_update)
    resultmenu: EnumProperty(items=r_menu, name="", description="Result type", update=r_update)
    zonemenu: EnumProperty(items=z_menu, name="", description="Location display", update=z_update)
    metricmenu: EnumProperty(items=m_menu, name="", description="Metric display", update=m_update)
    multfactor: FloatProperty(name="", description="Result multiplication factor", min=-10000, max=10000, default=1)
    statmenu: EnumProperty(items=[('Average', 'Average', 'Average Value'), ('Maximum', 'Maximum', 'Maximum Value'), ('Minimum', 'Minimum', 'Minimum Value'), ('Sum', 'Sum', 'Sum Value')],
                                      name="", description="Zone result", default='Average')

    def draw(self, context, layout, node, text):
        typedict = {"Time": [], "Frames": [], "Climate": ['climmenu'],
                    "Zone spatial": ("zonemenu", "zonermenu"), "Zone temporal": ("zonemenu", "zonermenu"),
                    "Embodied carbon": ("ecmenu", "ecrmenu"), "Linkage": ("linkmenu", "linkrmenu"),
                    "External": ("enmenu", "enrmenu"), "Position": ("posmenu", "posrmenu"),
                    "Camera": ("cammenu", "camrmenu"), "Power": ("powmenu", "powrmenu"),
                    "Probe": ("probemenu", "probermenu")}
        row = layout.row()

        if self.links and self.links[0].from_node.get('frames'):
            row.prop(self, "framemenu", text=text)
            row.prop(self, "resultmenu")

            if self.resultmenu not in ('None', 'Time', 'Frames', 'Climate'):
                row.prop(self, "zonemenu")

            if self.zonemenu not in ('None', 'Time'):
                row.prop(self, "metricmenu")

            if self.resultmenu not in ('Time', 'Frames') and self.node.timemenu in ('1', '2') and self.node.inputs['X-axis'].resultmenu == 'Time':
                row.prop(self, "statmenu")

            if self.resultmenu not in ('Time', 'Frames'):
                row.prop(self, 'multfactor')

    def draw_color(self, context, node):
        return (0.0, 1.0, 0.0, 0.75)

class ViEnRY1In(So_En_ResU):
    '''Chart y1 input socket'''
    bl_idname = 'ViEnRY1In'
    bl_label = 'Y-axis 1'


class ViEnRY2In(So_En_ResU):
    '''Chart y2 input socket'''
    bl_idname = 'ViEnRY2In'
    bl_label = 'Y-axis 2'


class ViEnRY3In(So_En_ResU):
    '''Chart y3 input socket'''
    bl_idname = 'ViEnRY3In'
    bl_label = 'Y-axis 3'


class No_Vi_Chart(Node, ViNodes):
    '''Node for 2D results plotting'''
    bl_idname = 'No_Vi_Chart'
    bl_label = 'VI Chart'

    def nodeupdate(self, context):
        self.update()

    def pmitems(self, context):
        return [tuple(p) for p in self['pmitems']]

    ctypes = [("0", "Line/Scatter", "Line/Scatter Plot")]
    charttype: EnumProperty(items=ctypes, name="", default="0")
    timemenu: EnumProperty(items=[("0", "Hourly", "Hourly results"),("1", "Daily", "Daily results"), ("2", "Monthly", "Monthly results")],
                           name="", description="Results frequency", default="0")
    bl_width_max = 800
    dpi: IntProperty(name='DPI', description="DPI of the shown figure", default=92, min=92)

    def init(self, context):
        self.inputs.new("ViEnRIn", "X-axis")
        self.inputs.new("ViEnRIn", "Y-axis 1")
        self.inputs["Y-axis 1"].hide = True
        self.inputs.new("ViEnRIn", "Y-axis 2")
        self.inputs["Y-axis 2"].hide = True
        self.inputs.new("ViEnRIn", "Y-axis 3")
        self.inputs["Y-axis 3"].hide = True
        self['Start'], self['End'] = 1, 365
        self['pmitems'] = [("0", "Static", "Static results")]
        self.update()

    def draw_buttons(self, context, layout):
        if self.inputs['X-axis'].links:
            innode = self.inputs['X-axis'].links[0].from_node
            rsx = self.inputs['X-axis']

            if innode.get('reslists') and len(innode['reslists']) > 1:
                if rsx.resultmenu == 'Time':
                    (sdate, edate) = retdates(self['Start'], self['End'], context.scene.vi_params.year)
                    label = "Start/End Day: {}/{} {}/{}".format(sdate.day, sdate.month, edate.day, edate.month)
                elif rsx.resultmenu == 'Frames':
                    label = "Frame"
                else:
                    label = 'X range'

                row = layout.row()
                row.label(text=label)
                row.prop(self, '["Start"]')
                row.prop(self, '["End"]')
                row = layout.row()
                row.label(text="Chart type:")
                row.prop(self, "charttype")

                if rsx.resultmenu == 'Time':
                    row.prop(self, "timemenu")

                row.prop(self, "dpi")

                r_lens = [sock.r_len for sock in self.inputs if sock.links]

                if len(r_lens) > 1 and all([r == r_lens[0] for r in r_lens]):
                    row = layout.row()
                    row.operator("node.chart", text = 'Create plot')
                    row = layout.row()
                    row.label(text = "------------------")
            else:
                row = layout.row()
                row.label(text = "No results")

    def update(self):
        if self.inputs.get('X-axis'):
            if self.inputs['X-axis'].links and self.inputs['X-axis'].links[0].from_node.get('reslists') and self.inputs['X-axis'].links[0].from_node['reslists']:
                rsx = self.inputs['X-axis']
                innode = rsx.links[0].from_node
                rl = innode['reslists']
                zrl = list(zip(*rl))
                rsx['resdict'] = {}
                resdict = {}

                for r in rl:
                    if r:
                        resdict[r[0]] = {} if r[0] not in resdict else resdict[r[0]]

                        try:
                            resdict[r[0]][r[1]] = {} if r[1] not in resdict[r[0]] else resdict[r[0]][r[1]]
                            resdict[r[0]][r[1]][r[2]] = [] if r[2] not in resdict[r[0]][r[1]] else resdict[r[0]][r[1]][r[2]]
                            resdict[r[0]][r[1]][r[2]].append(r[3])
                        except Exception:
                            pass

                rsx['resdict'] = resdict

                if resdict:
                    i = 0

                    while rsx.framemenu not in [r[0] for r in rsx.f_menu(0)] and i < 5:
                        try:
                            rsx.framemenu = rsx.f_menu(self)[0][0]
                        except Exception:
                            i += 1

                    i = 0

                    while rsx.resultmenu not in [r[0] for r in rsx.r_menu(0)] and i < 5:
                        try:
                            rsx.resultmenu = rsx.r_menu(0)[0][0]
                        except Exception:
                            i += 1

                    i = 0

                    while rsx.zonemenu not in [r[0] for r in rsx.z_menu(0)] and i < 5:
                        try:
                            rsx.zonemenu = rsx.z_menu(0)[0][0]
                        except Exception:
                            i += 1

                    i = 0

                    while rsx.metricmenu not in [r[0] for r in rsx.m_menu(0)] and i < 5:
                        try:
                            rsx.metricmenu = rsx.m_menu(0)[0][0]
                        except Exception:
                            i += 1

                    rsx.metricmenu = rsx.metricmenu

                    if rsx.resultmenu == 'Time':
                        startday = datetime.datetime(bpy.context.scene.vi_params.year, int(zrl[4][zrl[3].index('Month')].split()[0]), int(zrl[4][zrl[3].index('Day')].split()[0])).timetuple().tm_yday
                        endday = datetime.datetime(bpy.context.scene.vi_params.year, int(zrl[4][zrl[3].index('Month')].split()[-1]), int(zrl[4][zrl[3].index('Day')].split()[-1])).timetuple().tm_yday
                        self["_RNA_UI"] = {"Start": {"min":startday, "max":endday}, "End": {"min":startday, "max":endday}}
                        self['Start'], self['End'] = startday, endday

                    elif rsx.resultmenu == 'Frames':
                        frames = [int(k) for k in set(zrl[0]) if k != 'All']
                        startframe, endframe = min(frames), max(frames)
                        frame = 'All'
                        self["_RNA_UI"] = {"Start": {"min":startframe, "max":endframe}, "End": {"min":startframe, "max":endframe}}
                        self['Start'], self['End'] = startframe, endframe

                    else:
                        xs = range(1, 1 + [len(res[4].split()) for res in rl if res[0] == rsx.framemenu and res[1] == rsx.resultmenu and res[2] == rsx.zonemenu and res[3] == rsx.metricmenu][0])
                        (startx, endx) = (min(xs), max(xs)) if len(xs) > 1 else (xs[0], xs[0])
                        self["_RNA_UI"] = {"Start": {"min":startx, "max":endx}, "End": {"min":startx, "max":endx}}
                        self['Start'], self['End'] = startx, endx

                    if self.inputs.get('Y-axis 1'):
                        self.inputs['Y-axis 1'].hide = False

                    if self.inputs.get('Y-axis 1'):
                        if self.inputs['Y-axis 1'].links and self.inputs['Y-axis 1'].links[0].from_node['reslists']:
                            rsy1 = self.inputs['Y-axis 1']
                            rsy1.update_menus(('Y-axis 1', ))

                            i = 0

                            while rsy1.framemenu not in [r[0] for r in rsy1.f_menu(0)] and i < 5:
                                try:
                                    rsy1.framemenu = rsy1.f_menu(0)[0][0]
                                except Exception:
                                    i += 1
                            i = 0

                            while rsy1.resultmenu not in [r[0] for r in rsy1.r_menu(0)]  and i < 5:
                                try:
                                    rsy1.resultmenu = rsy1.r_menu(0)[0][0]
                                except Exception:
                                    i += 1
                            i = 0
                            while rsy1.zonemenu not in [r[0] for r in rsy1.z_menu(0)]  and i < 5:
                                try:
                                    rsy1.zonemenu = rsy1.z_menu(0)[0][0]
                                except Exception:
                                    i += 1
                            i = 0

                            while rsy1.metricmenu not in [r[0] for r in rsy1.m_menu(0)]  and i < 5:
                                try:
                                    rsy1.metricmenu = rsy1.m_menu(0)[0][0]
                                except Exception:
                                    i += 1

                            rsy1.framemenu = rsy1.framemenu
                            rsy1.resultmenu = rsy1.resultmenu
                            rsy1.zonemenu = rsy1.zonemenu
                            rsy1.metricmenu = rsy1.metricmenu

                            if self.inputs.get('Y-axis 2'):
                                self.inputs['Y-axis 2'].hide = False

                        else:
                            if self.inputs.get('Y-axis 2'):
                                self.inputs['Y-axis 2'].hide = True

                    if self.inputs.get('Y-axis 2'):
                        if self.inputs['Y-axis 2'].links and self.inputs['Y-axis 2'].links[0].from_node['reslists']:
                            rsy2 = self.inputs['Y-axis 2']
                            rsy2.update_menus(('Y-axis 2', ))

                            i = 0

                            while rsy2.framemenu not in [r[0] for r in rsy2.f_menu(0)] and i < 5:
                                try:
                                    rsy2.framemenu = rsy2.f_menu(0)[0][0]
                                except Exception:
                                    i += 1
                            i = 0

                            while rsy2.resultmenu not in [r[0] for r in rsy2.r_menu(0)]  and i < 5:
                                try:
                                    rsy2.resultmenu = rsy2.r_menu(0)[0][0]
                                except Exception:
                                    i += 1
                            i = 0
                            while rsy2.zonemenu not in [r[0] for r in rsy2.z_menu(0)]  and i < 5:
                                try:
                                    rsy2.zonemenu = rsy2.z_menu(0)[0][0]
                                except Exception:
                                    i += 1
                            i = 0

                            while rsy2.metricmenu not in [r[0] for r in rsy2.m_menu(0)]  and i < 5:
                                try:
                                    rsy2.metricmenu = rsy2.m_menu(0)[0][0]
                                except Exception:
                                    i += 1

                            rsy2.framemenu = rsy2.framemenu
                            rsy2.resultmenu = rsy2.resultmenu
                            rsy2.zonemenu = rsy2.zonemenu
                            rsy2.metricmenu = rsy2.metricmenu

                            if self.inputs.get('Y-axis 3'):
                                self.inputs['Y-axis 3'].hide = False
                        else:
                            self.inputs['Y-axis 3'].hide = True

                    if self.inputs.get('Y-axis 3'):
                        if self.inputs['Y-axis 3'].links and self.inputs['Y-axis 3'].links[0].from_node['reslists']:
                            rsy3 = self.inputs['Y-axis 3']
                            rsy3.update_menus(('Y-axis 3', ))

                            i = 0

                            while rsy3.framemenu not in [r[0] for r in rsy3.f_menu(0)] and i < 5:
                                try:
                                    rsy3.framemenu = rsy3.f_menu(0)[0][0]
                                except Exception:
                                    i += 1
                            i = 0

                            while rsy3.resultmenu not in [r[0] for r in rsy3.r_menu(0)]  and i < 5:
                                try:
                                    rsy3.resultmenu = rsy3.r_menu(0)[0][0]
                                except Exception:
                                    i += 1
                            i = 0
                            while rsy3.zonemenu not in [r[0] for r in rsy3.z_menu(0)]  and i < 5:
                                try:
                                    rsy3.zonemenu = rsy3.z_menu(0)[0][0]
                                except Exception:
                                    i += 1
                            i = 0

                            while rsy3.metricmenu not in [r[0] for r in rsy3.m_menu(0)]  and i < 5:
                                try:
                                    rsy3.metricmenu = rsy3.m_menu(0)[0][0]
                                except Exception:
                                    i += 1

                            rsy3.framemenu = rsy3.framemenu
                            rsy3.resultmenu = rsy3.resultmenu
                            rsy3.zonemenu = rsy3.zonemenu
                            rsy3.metricmenu = rsy3.metricmenu

                else:
                    self.inputs['Y-axis 1']['resdict'], self.inputs['Y-axis 2']['resdict'], self.inputs['Y-axis 3']['resdict'] = {}, {}, {}
                    if self.inputs.get('Y-axis 1'):
                        self.inputs['Y-axis 1'].hide = True
                    if self.inputs.get('Y-axis 2'):
                        self.inputs['Y-axis 2'].hide = True
                    if self.inputs.get('Y-axis 3'):
                        self.inputs['Y-axis 3'].hide = True
            else:
                if self.inputs.get('Y-axis 1'):
                    self.inputs['Y-axis 1'].hide = True


class No_Vi_HMChart(Node, ViNodes):
    '''Node for 2D results plotting'''
    bl_idname = 'No_Vi_HMChart'
    bl_label = 'VI Heatmap'

    def update(self):
        resdict = {}
        if self.inputs['Results in'].links:
            innode = self.inputs['Results in'].links[0].from_node

            for rl in innode['reslists']:
                if rl and rl[0] != 'All':
                    resdict[rl[0]] = {} if rl[0] not in resdict else resdict[rl[0]]

                try:
                    if rl and rl[1] in ('Zone temporal', 'Power', 'Climate'):
                        resdict[rl[0]][rl[1]] = {} if rl[1] not in resdict[rl[0]] else resdict[rl[0]][rl[1]]
                        resdict[rl[0]][rl[1]][rl[2]] = [] if rl[2] not in resdict[rl[0]][rl[1]] else resdict[rl[0]][rl[1]][rl[2]]
                        resdict[rl[0]][rl[1]][rl[2]].append(rl[3])
                except Exception:
                    pass

            self['resdict'] = resdict
            self['times'] = list(resdict.keys())

            i = 0

            while self.framemenu not in [r[0] for r in self.f_menu(0)] and i < 5:
                try:
                    self.framemenu = self.f_menu(self)[0][0] # list(resdict.keys())[0]
                except Exception:
                    i += 1

            i = 0

            while self.resmenu not in [r[0] for r in self.r_menu(0)] and i < 5:
                try:
                    self.resmenu = self.r_menu(0)[0][0]
                except Exception:
                    i += 1

            i = 0

            while self.locmenu not in [r[0] for r in self.z_menu(0)] and i < 5:
                try:
                    self.locmenu = self.z_menu(0)[0][0]
                except Exception:
                    i += 1

            i = 0

            while self.metricmenu not in [r[0] for r in self.m_menu(0)] and i < 5:
                try:
                    self.metricmenu = self.m_menu(0)[0][0]
                except Exception:
                    i += 1

    def f_menu(self, context):
        try:
            frames = [(f'{frame}', f'{frame}', f'Frame {frame}') for frame in self['resdict'].keys()]

            if frames:
                return frames
            else:
                return [('None', 'None', 'None')]

        except Exception:
            return [('None', 'None', 'None')]

    def r_menu(self, context):
        try:
            ress = [(f'{res}', f'{res}', f'Frame {res}') for res in self['resdict'][self.framemenu].keys()]

            if ress:
                return ress
            else:
                return [('None', 'None', 'None')]

        except Exception:
            return [('None', 'None', 'None')]

    def z_menu(self, context):
        try:
            return [(f'{res}', f'{res}', f'Zone {res}') for res in self['resdict'][self.framemenu][self.resmenu].keys()]

        except Exception:
            try:
                r = list(self['resdict'][self.framemenu].keys())[0]
                z = list(self['resdict'][self.framemenu][r].keys())[0]
                self.locmenu = list(self['resdict'][self.framemenu][r].keys())[0]
                return [(f'{res}', f'{res}', f'Zone {res}') for res in self['resdict'][self.framemenu][r].keys()]
            except Exception:
                return [('None', 'None', 'None')]


    def m_menu(self, context):
        try:
            return [(f'{res}', f'{res}', f'Frame {res}') for res in self['resdict'][self.framemenu][self.resmenu][self.locmenu]]
        except Exception:
            try:
                r = list(self['resdict'][self.framemenu].keys())[0]
                z = list(self['resdict'][self.framemenu][r].keys())[0]
                return [(f'{res}', f'{res}', f'Frame {res}') for res in self['resdict'][self.framemenu][r][z]]
            except Exception:
                return [('None', 'None', 'None')]

    def mupdate(self, context):
        if self.locmenu not in [l[0] for l in self.z_menu(context)]:
            self.locmenu = self.z_menu(context)[0][0]

        if self.metricmenu not in [l[0] for l in self.m_menu(context)]:
            self.metricmenu = self.m_menu(context)[0][0]

    def nodeupdate(self, context):
        if self.inputs['Results in'].links:
            innode = self.inputs['Results in'].links[0].from_node
            rl = innode['reslists']
            self.mupdate(context)

            for r in rl:
                if r[0] == self.framemenu:
                    if r[1] == 'Time':
                        if r[3] == 'DOS':
                            self.x = (array([float(r) for r in r[4].split()]) + 0.5)
                            dno = len(unique(self.x))

                        elif r[3] == 'Hour':
                            self.y = (array([float(r) for r in r[4].split()]) + 0.5)
                            hno = len(unique(self.y))

                    elif r[1] == self.resmenu:
                        if self.resmenu == 'Climate':
                            if r[3] == self.metricmenu:
                                self.z = array([float(r) for r in r[4].split()])

                        elif r[3] == self.metricmenu and r[2] == self.locmenu:
                            self.z = array([float(r) for r in r[4].split()])
            try:
                self.x = self.x.reshape(dno, hno)
                self.y = self.y.reshape(dno, hno)
                self.z = self.z.reshape(dno, hno)

            except Exception:
                logentry('Mis-match in result length. Try reconnecting the Heatmap chart node')

    dpi: IntProperty(name='DPI', description="DPI of the shown figure", default=92, min=92)
    framemenu: EnumProperty(items=f_menu, name="", description="Frame number")
    resmenu: EnumProperty(items=r_menu, name="", description="Result type", update=mupdate)
    locmenu: EnumProperty(items=z_menu, name="", description="Result location", update=mupdate)
    metricmenu: EnumProperty(items=m_menu, name="", description="Result metric", update=mupdate)
    metricrange: EnumProperty(items=[('0', 'Auto', 'Automatic range based on max/min values'), ('1', 'Custom', 'Custom range based on input values')], name="", description="Result metric")
    cf: BoolProperty(name="", description="Contour fill", default=0)
    cl: BoolProperty(name="", description="Contour fill", default=0)
    lvals: StringProperty(name="", description="Space separated contour values", default="")
    lw: FloatProperty(name="", description="Line width", min=0.0, max=10, default=0.1)
    clevels: IntProperty(name='', description="Number of contour levels", default=10, min=1)
    daystart: IntProperty(name='', description="Start day", default=1, min=1, max=365)
    dayend: IntProperty(name='', description="End day", default=365, min=1, max=365)
    hourstart: IntProperty(name='', description="Start hour", default=1, min=1, max=365)
    hourend: IntProperty(name='', description="End hour", default=24, min=1, max=365)
    varmin: IntProperty(name='', description="Variable minimum", default=0)
    varmax: IntProperty(name='', description="Varaible maximum", default=20)
    grid: BoolProperty(name="", description="Grid", default=0)
    x, y, z = array([]), array([]), array([])

    def init(self, context):
        self['times'] = [('', '', "")]
        self['rtypes'] = [('', '', "")]
        self['locs'] = [('', '', "")]
        self['metrics'] = [('', '', "")]
        self.inputs.new("So_Vi_Res", "Results in")

    def draw_buttons(self, context, layout):
        svp = context.scene.vi_params

        if self.inputs['Results in'].links:
            innode = self.inputs['Results in'].links[0].from_node

            if innode.get('reslists'):
                row = layout.row()
                row.prop(self, "framemenu")
                row.prop(self, "dpi")
                row = layout.row()
                row.prop(self, "resmenu")
                row.prop(self, "locmenu")
                row.prop(self, "metricmenu")
                row = layout.row()
                row.label(text = 'Days:')
                row.prop(self, 'daystart')
                row.prop(self, 'dayend')
                row = layout.row()
                row.label(text = 'Hours:')
                row.prop(self, 'hourstart')
                row.prop(self, 'hourend')
                row = layout.row()
                newrow(layout, 'Metric range:', self, "metricrange")

                if self.metricrange == '1':
                    row = layout.row()
                    row.label(text='Range:')
                    row.prop(self, 'varmin')
                    row.prop(self, 'varmax')

                newrow(layout, 'Colour map:', svp, "vi_leg_col")
                newrow(layout, 'Contour lines:', self, "cl")

                if self.cl:
                    newrow(layout, 'Line width:', self, "lw")
                    newrow(layout, 'Contour values:', self, "lvals")

                newrow(layout, 'Contour fill:', self, "cf")

                if self.cf:
                    newrow(layout, 'Grid display:', self, "grid")

                if self.cl or self.cf:
                    newrow(layout, 'Contour levels:', self, "clevels")

                if self.framemenu and self.metricmenu != 'None':
                    row = layout.row()
                    row.operator("node.hmchart", text='Create heatmap')

class No_Vi_Metrics(Node, ViNodes):
    '''Node for result metrics'''
    bl_idname = 'No_Vi_Metrics'
    bl_label = 'VI Metrics'

    def zupdate(self, context):
        self.update()

    def zitems(self, context):
        if self.inputs[0].links:
            try:
                return [tuple(z) for z in self['znames']]
            except Exception:
                return [('None', 'None', 'None')]

        # elif self.metric == '3' and self.em_menu == '0':
        #     if [o for o in bpy.context.visible_objects if o.vi_params.get('ecdict')]:
        #         return [(o.name, o.name, o.name) for o in bpy.context.visible_objects if o.vi_params.get('ecdict')]  + [('All', 'All', 'All objects')]
        #     else:
        #         return [('None', 'None', 'None')]
        else:
            return [('None', 'None', 'None')]

    def frames(self, context):
        if self.inputs[0].links:
            try:
                # if self.frame_menu not in [f[0] for f in self['frames']]:
                #     print('frame', self.frame_menu, [f[0] for f in self['frames']])
                    #self.frame_menu = self['frames'][0][0]
                return [tuple(f) for f in self['frames']]

            except Exception:
                #if self.frame_menu != 'None':
                    #self.frame_menu = 'None'
                return [('None', 'None', 'None')]
        else:
            #if self.frame_menu != 'None':
                #self.frame_menu = 'None'
            return [('None', 'None', 'None')]

    def probes(self, context):
        if self.inputs[0].links and self['rl']:
            probes = set([z[3] for z in self['rl'] if z[3] in ('Temperature', 'Pressure', 'Velocity')])

            if probes:
                return [(m.lower(), m, 'Probe metric') for m in probes if m != 'Steps'] + [('wpc', 'WPC', 'Probe menu')]
            else:
                return [('None', 'None', 'None')]
        else:
            return [('None', 'None', 'None')]

    def ec_types(self, context):
        if self.inputs[0].links and self['rl']:
            ec_typemenu = []
            try:
                ec_types = set([z[3] for z in self['rl']])
                if 'Object EC (kgCO2e)' in ec_types:
                    ec_typemenu.append('Object')
                if 'Surface EC (kgCO2e/y)' in ec_types:
                    ec_typemenu.append('Surface')
                if 'Zone EC (kgCO2e/y)' in ec_types:
                    ec_typemenu.append('Zone')

                return [(ect, ect, 'EC type') for ect in ec_typemenu] if ec_typemenu else [('None', 'None', 'None')]
            except Exception:
                return [('None', 'None', 'None')]
        else:
            return [('None', 'None', 'None')]

    metric: EnumProperty(items=[("0", "Energy", "Energy results"),
                                ("1", "Lighting", "Lighting results"),
                                ("2", "Flow", "Flow results"),
                                ("3", "Embodied carbon", "Embodied carbon results"),
                                ("4", "Comfort", "Comfort results"),
                                ("5", "IAQ", "Internal air quality results"),
                                ("6", "WLC", "WHole life carbon results")],
                                name="", description="Results type", default="0", update=zupdate)
    energy_menu: EnumProperty(items=[("0", "SAP", "SAP results"),
                                    ("1", "RIBA 2030", "RIBA 2030 results"),
                                    ("2", "PassivHaus", "PassivHaus reults")],
                                    name="", description="Results metric", default="0", update=zupdate)
    light_menu: EnumProperty(items=[("0", "BREEAM", "BREEAM HEA1 results"),
                                    ("1", "LEED", "LEED v4 results"),
                                    ("2", "RIBA 2030", "RIBA 2030 results")],
                                    name="", description="Results metric", default="0", update=zupdate)
    em_menu: EnumProperty(items=ec_types, name="", description="Results metric", update=zupdate)
    leed_menu: BoolProperty(name="", description="LEED space type", default=0)
    riba_menu: EnumProperty(items=[("0", "Domestic", "Domestic scenario"),
                                    ("1", "Office", "Office scenario"),
                                    ("2", "School", "School scenario")],
                                    name="", description="RIBA building class", default="0", update=zupdate)
    breeam_menu: EnumProperty(items=[("0", "Education", "Education scenario"),
                                    ("1", "Healthcare", "Healthcare scenario"),
                                    ("2", "Multi-residential", "Multi-residential scenario"),
                                    ("3", "Retail", "Retail scenario"),
                                    ("4", "Prison", "Prison scenario"),
                                    ("5", "Office", "Office scenario"),
                                    ("6", "Creche", "Creche scenario"),
                                    ("7", "Other", "Other scenario")],
                                    name="", description="BREEAM space type", default="0", update=zupdate)
    breeam_edumenu: EnumProperty(items=[("0", "School", "School context"),
                                        ("1", "Higher education", "Higher education scenario")],
                                        name="", description="BREEAM education space type", default="0", update=zupdate)
    breeam_healthmenu: EnumProperty(items=[("0", "Staff/public", "Staff/public context"),
                                        ("1", "Patient", "Patient scenario")],
                                        name="", description="BREEAM healthcare space type", default="0", update=zupdate)
    breeam_multimenu: EnumProperty(items=[("0", "Kitchen", "Staff/public context"),
                                        ("1", "Living", "Patient scenario"),
                                        ("2", "Communal", "Patient scenario")],
                                        name="", description="BREEAM multi-residential space type", default="0", update=zupdate)
    breeam_retailmenu: EnumProperty(items=[("0", "Sales", "Staff/public context"),
                                        ("1", "Other", "Patient scenario")],
                                    name="", description="BREEAM retail space type", default="0", update=zupdate)
    breeam_prisonmenu: EnumProperty(items=[("0", "Cells", "Custody cells context"),
                                          ("2", "Atrium", "Atrium scenario"),
                                          ("3", "Patient", "Patient context"),
                                          ("4", "Lecture", "Lecture context")],
                                   name="", description="BREEAM other space type", default="0", update=zupdate)
    g_roof: BoolProperty(name="", description="Glazed roof", default=0, update=zupdate)
    com_menu: EnumProperty(items=[("0", "Overheating", "Overheating analysis")],
                                name="", description="Comfort type", default="0", update=zupdate)
    iaq_menu: EnumProperty(items=[("0", "RIBA 2030", "RIBA 2030 CO2 Criteria")],
                                name="", description="IAQ type", default="0", update=zupdate)
    ec_menu: EnumProperty(items=[("0", "RIBA 2030", "RIBA 2030 Embodied Carbon Criteria")],
                                name="", description="Embodied carbon standard", default="0", update=zupdate)
    zone_menu: EnumProperty(items=zitems,
                name="", description="Zone results", update=zupdate)
    frame_menu: EnumProperty(items=frames,
                name="", description="Frame results", update=zupdate)
    mod: FloatProperty(name="kWh/m2/y", description="Energy modifier (kWh/m2/y)", update=zupdate)
    hwmod: FloatProperty(name="kWh/m2/y", description="Hot water consumption (kWh/m2/y)", update=zupdate)
    heat_type: EnumProperty(items=[('0', 'Gas', 'Gas based heating'), ('1', 'Electric', 'Electric based hheating')], name="", description="Heating energy source", default='0', update=zupdate)
    elec_cop: FloatProperty(name="", description="Coefficient-of-performance of the electrical heating system", default=1.0, update=zupdate)
    hw_cop: FloatProperty(name="", description="Coefficient-of-performance of the electrical hot water system", default=1.0, update=zupdate)
    ac_cop: FloatProperty(name="", description="Coefficient-of-performance of the electrical cooling system", default=3.0, update=zupdate)
    gas_eff: FloatProperty(name="%", description="Efficiency of the gas heating system", default=90, min=1.0, max=100.0, update=zupdate)
    carb_fac: FloatProperty(name="kgCO2/kWh", description="Electrical grid carbon factor", default=0.21, min=0.01, max=1, update=zupdate)
    carb_annc: FloatProperty(name="%", description="Annual change in carbon factor", default=-0.5, min=-100, max=100, update=zupdate)
    ec_years: IntProperty(name="Years", description="Timespan for embodied carbon calculations", default=60, min=1, max=100, update=zupdate)
    probe_menu: EnumProperty(items=probes, name="", description="Probe results", update=zupdate)
    ws: FloatProperty(name="m/s", description="Freesteam wind speed", update=zupdate)
    occ: BoolProperty(name="", description="Only occupied hours", update=zupdate)

    def init(self, context):
        self['res'] = {}
        self.inputs.new('So_Vi_Res', 'Results in')
        self.outputs.new('So_Vi_Res', 'Results out')
        # self.outputs['Results out'].hide = True
        self['riba_en'] = {'0': 35, '1': 55, '2': 60}

    def draw_buttons(self, context, layout):
        if self.inputs[0].links:
            if self.inputs[0].links[0].from_node.get('new_res'):
                row = layout.row()
                row.label(text = "**Reconnect Node**")

            else:
                newrow(layout, 'Type:', self, "metric")

                if self.metric == '0':
                    newrow(layout, 'Metric:', self, "energy_menu")
                elif self.metric == '1':
                    newrow(layout, 'Metric:', self, "light_menu")
                # elif self.metric == '3':
                #     newrow(layout, 'Embodied standard:', self, "ec_menu")

                newrow(layout, 'Frame', self, "frame_menu")
                if self.frame_menu and self.frame_menu != 'None':

                    if self.metric == '3':
                        newrow(layout, 'Embodied type:', self, "em_menu")

                    newrow(layout, ('Zone:', f'{self.em_menu}:')[self.metric == '3'], self, "zone_menu")

                    if self.metric == '0':
                        if self['res'] and self['res'].get('hkwh'):
                            newrow(layout, 'Heating type:', self, 'heat_type')

                            if self.heat_type == '0':
                                newrow(layout, 'Heating efficiency:', self, 'gas_eff')
                            else:
                                newrow(layout, 'Heating COP:', self, 'elec_cop')
                                newrow(layout, 'Hot water COP:', self, 'hw_cop')

                            if self['res'].get('ckwh'):
                                newrow(layout, 'Air-con COP:', self, 'ac_cop')

                            # newrow(layout, 'Hot water:', self, 'hwmod')
                            newrow(layout, 'Unregulated:', self, 'mod')
                            newrow(layout, 'Carbon factor:', self, 'carb_fac')

                        if self.energy_menu == '0':
                            if self['res'].get('fa'):
                                row = layout.row()
                                pvkwh = self['res']['pvkwh'] if self['res']['pvkwh'] == 'N/A' else "{:.2f}".format(self['res']['pvkwh'])
                                row.label(text="PV (kWh): {}".format(pvkwh))
                                pva = "{:.2f}".format(self['res']['pvkwh']/self['res']['fa']) if self['res']['fa'] != 'N/A' and self['res']['fa'] > 0 else 'N/A'
                                row = layout.row()
                                row.label(text="PV (kWh/m2): {}".format(pva))
                                row = layout.row()
                                hkwh = self['res']['hkwh'] if self['res']['hkwh'] == 'N/A' else "{:.2f}".format(self['res']['hkwh'])
                                row.label(text="Heating (kWh): {}".format(hkwh))
                                row = layout.row()
                                ha = "{:.2f}".format(self['res']['hkwh']/self['res']['fa']) if self['res']['fa'] != 'N/A' and self['res']['fa'] > 0 else 'N/A'
                                row.label(text="Heating (kWh/m2): {}".format(ha))
                                row = layout.row()
                                ckwh = self['res']['pvkwh'] if self['res']['pvkwh'] == 'N/A' else "{:.2f}".format(self['res']['ckwh'])
                                row.label(text="Cooling (kWh): {}".format(ckwh))
                                row = layout.row()
                                ca = "{:.2f}".format(self['res']['ckwh']/self['res']['fa']) if self['res']['fa'] != 'N/A' and self['res']['fa'] > 0 else 'N/A'
                                row.label(text="Cooling (kWh/m2): {}".format(ca))

                                if self.zone_menu == 'All':
                                    row = layout.row()
                                    wkwh = self['res']['wkwh'] if self['res']['wkwh'] == 'N/A' else "{:.2f}".format(self['res']['wkwh'])
                                    row.label(text="Hot water (kWh): {}".format(wkwh))
                                    row = layout.row()
                                    wkwhm2 = "{:.2f}".format(self['res']['wkwh']/self['res']['fa']) if self['res']['fa'] != 'N/A' and self['res']['fa'] > 0 else 'N/A'
                                    row.label(text="Hot water (kWh/m2): {}".format(wkwhm2))
                                    row = layout.row()
                                    ecf = "{:.2f}".format(self['res']['ECF']) if self['res']['ECF'] != 'N/A' else 'N/A'
                                    row.label(text="ECF: {}".format(ecf))
                                    row = layout.row()
                                    epc = "{:.0f}".format(self['res']['EPC']) if self['res']['EPC'] != 'N/A' else 'N/A'
                                    row.label(text="EPC: {} ({})".format(epc, self['res']['EPCL']))

                        elif self.energy_menu == '1':
                            if self['res'].get('fa'):
                                if self['res'].get('totkwh'):
                                    newrow(layout, 'Type', self, 'riba_menu')
                                    tar = self['riba_en'][self.riba_menu]
                                    epass = '(FAIL kWh/m2 > {})'.format(tar) if self['res']['totkwh']/self['res']['fa'] > tar else '(PASS kWh/m2 <= {})'.format(tar)
    #                                shpass = '(FAIL kWh/m2 > {})'.format(20) if self['res']['totkwh']/self['res']['fa'] > 20 else '(PASS kWh/m2 <= {})'.format(20)
                                    row = layout.row()
                                    row.label(text = "Space heating (kWh/m2): {:.1f}".format(self['res']['hkwh']/self['res']['fa']))
                                    row = layout.row()
                                    row.label(text = "Space cooling (kWh/m2): {:.1f}".format(self['res']['ckwh']/self['res']['fa']))
                                    row = layout.row()
                                    row.label(text = "Power production (kWh/m2): {:.1f}".format(self['res']['pvkwh']/self['res']['fa']))
                                    row = layout.row()
                                    row.label(text="Gross operational (kWh/m2): {:.1f} {}".format(self['res']['totkwh']/self['res']['fa'], epass))
                                    row = layout.row()
                                    row.label(text="Net operational (kWh/m2): {:.1f}".format(self['res']['netkwh']/self['res']['fa']))

                                elif self['res'].get('hkwh'):
                                    row = layout.row()
                                    row.label(text="Space heating (kWh): {:.1f}".format(self['res']['hkwh']))
                                    row = layout.row()
                                    row.label(text="Space heating (kWh/m2): {:.1f}".format(self['res']['hkwh']/self['res']['fa']))

                                    if self['res'].get('ckwh'):
                                        row = layout.row()
                                        row.label(text="Space cooling (kWh): {:.1f}".format(self['res']['ckwh']))
                                        row = layout.row()
                                        row.label(text="Space cooling (kWh/m2): {:.1f}".format(self['res']['ckwh']/self['res']['fa']))
                            else:
                                row = layout.row()
                                row.label(text="No floor area")

                        elif self.energy_menu == '2':
                            if self['res'].get('fa'):
                                if self['res'].get('hkwh'):
                                    epass = '(FAIL kWh/m2 > {})'.format(15) if self['res']['hkwh']/self['res']['fa'] > 15 else '(PASS kWh/m2 <= {})'.format(15)
                                    row = layout.row()
                                    row.label(text="Space heating (kWh/m2): {:.1f} {}".format(self['res']['hkwh']/self['res']['fa'], epass))

                    elif self.metric == '1':
                        if self.light_menu == '0':
                            areaDF = 'N/A' if self['res']['areaDF'] < 0 else self['res']['areaDF']
                            avDF = 'N/A' if self['res']['avDF'] < 0 else self['res']['avDF']
                            minDF = 'N/A' if self['res']['minDF'] < 0 else self['res']['minDF']
                            ratioDF = 'N/A' if self['res']['ratioDF'] < 0 else self['res']['ratioDF']
                            sda = 'N/A' if self['res']['sda'] < 0 else self['res']['sda']
                            credits = 0
                            newrow(layout, 'Space:', self, "breeam_menu")

                            if self.breeam_menu == '0':
                                # if avDF >= 2 and ratioDF >= 0.3:
                                newrow(layout, 'Education space:', self, "breeam_edumenu")

                            elif self.breeam_menu == '1':
                                    # if ratioDF >= 0.3:
                                    newrow(layout, 'Health space:', self, "breeam_healthmenu")

                            elif self.breeam_menu == '2':
                                # if ratioDF >= 0.3:
                                newrow(layout, 'Multi-res space:', self, "breeam_multimenu")

                            elif self.breeam_menu == '3':
                                newrow(layout, 'Retail space:', self, "breeam_retailmenu")

                            elif self.breeam_menu == '4':
                                newrow(layout, 'Prison space:', self, "breeam_prisonmenu")

                            newrow(layout, 'Glazed roof:', self, "g_roof")

                            if 'N/A' not in (areaDF, avDF):
                                if self.breeam_menu == '0':
                                    # if avDF >= 2 and ratioDF >= 0.3:
                                    # newrow(layout, 'Education space:', self, "breeam_edumenu")
                                    if self.breeam_edumenu == '0':
                                        if areaDF >= 80:
                                            credits = 2
                                    elif self.breeam_edumenu == '1':
                                        if areaDF >= 60:
                                            credits = 1
                                        if areaDF >= 80:
                                            credits = 2

                                elif self.breeam_menu == '1':
                                    # if ratioDF >= 0.3:
                                    # newrow(layout, 'Health space:', self, "breeam_healthmenu")
                                    if self.breeam_healthmenu == '0':
                                        if areaDF >= 80:
                                            credits = 2

                                    elif self.breeam_healthmenu == '1':
                                        if areaDF >= 80:
                                            credits = 1 if avDF < 3 else 2

                                elif self.breeam_menu == '2':
                                    # if ratioDF >= 0.3:
                                    # newrow(layout, 'Multi-res space:', self, "breeam_multimenu")
                                    if areaDF >= 80:
                                        credits = 2

                                elif self.breeam_menu == '3':
                                    # newrow(layout, 'Retail space:', self, "breeam_retailmenu")
                                    if self.breeam_retailmenu == '0':
                                        if areaDF >= 35:
                                            credits = 1

                                    elif self.breeam_retailmenu == '1':
                                        if areaDF >= 80:
                                            credits = 1

                                elif self.breeam_menu == '4':
                                    # newrow(layout, 'Other space:', self, "breeam_othermenu")
                                    if self.breeam_prisonmenu == '0':
                                        if areaDF >= 80:
                                            credits = 2

                                    elif self.breeam_prisonmenu == '1':
                                        if areaDF >= 80:
                                            credits = 2

                                    elif self.breeam_prisonmenu == '2':
                                        if areaDF >= 80:
                                            credits = 2

                                    elif self.breeam_prisonmenu == '3':
                                        if areaDF >= 80:
                                            credits = 2

                                    elif self.breeam_prisonmenu == '4':
                                        if areaDF >= 80:
                                            credits = 2

                                elif self.breeam_menu == '5':
                                    if areaDF >= 80:
                                        credits = 2

                                elif self.breeam_menu == '6':
                                    if areaDF >= 80:
                                        credits = 2

                                elif self.breeam_menu == '7':
                                    if areaDF >= 80:
                                        credits = 1

                                row = layout.row()
                                row.label(text="Average DF: {}%".format(avDF))
                                row = layout.row()
                                row.label(text="Minimum DF: {}%".format(minDF))
                                row = layout.row()
                                row.label(text="Uniformity: {}".format(ratioDF))
                                row = layout.row()
                                row.label(text="Compliant area: {}%".format(areaDF))
                                row = layout.row()
                                row.label(text="Credits: {}".format(credits))

                            elif sda != 'N/A':
                                row = layout.row()
                                row.label(text=f"Hours above target: {sda}")
                                # row = layout.row()
                                # row.label(text=f"Area above target: {sdaarea}")

                        elif self.light_menu == '1':
                            newrow(layout, 'Healthcare', self, 'leed_menu')
                            (self['res']['sdapass'], self['res']['tc']) = ((75, 90), 2) if self.leed_menu else ((55, 75), 3)

                            if self['res'] and self['res'].get('ase') > -1:
                                if self['res']['ase'] < 0:
                                    (sda, ase, o1) = ('sDA300 (%): N/A', 'ASE1000 (hours): N/A', 'Total credits: N/A')
                                else:
                                    sdares = self['res']['sdapa'] if self.leed_menu else self['res']['sda']

                                    if self['res']['ase'] <= 10:
                                        if round(sdares, 3) >= (55, 75)[self.leed_menu]:
                                            self['res']['o1'] = (2, 1)[self.leed_menu]
                                        else:
                                            self['res']['o1'] = 0

                                        if round(sdares, 3) >= (75, 90)[self.leed_menu]:
                                            self['res']['o1'] = (3, 2)[self.leed_menu]
                                    else:
                                        self['res']['o1'] = 0

                                    (sda, ase, o1) = ('sDA300/50% (% area): {0:.1f} | > ({1[0]}, {1[1]}) | {2}'.format(sdares, self['res']['sdapass'], ('Pass', 'Fail')[round(sdares, 3) < self['res']['sdapass'][0]]),
                                                    'ASE1000/250 (% area): {:.1f} | < 10 | {}'.format(self['res']['ase'], ('Pass', 'Fail')[self['res']['ase'] > 10]),
                                                    'Total credits: {}'.format(self['res']['o1']))

                                if self.leed_menu:
                                    row = layout.row()
                                    row.label(text='Perimeter area: {:.1f} | {}'.format(self['res']['sv'], ('Pass', 'Fail')[self['res']['sv'] < 90]))

                                row = layout.row()
                                row.label(text=sda)
                                row = layout.row()
                                row.label(text=ase)
                                row = layout.row()
                                row.label(text=o1)

                        elif self.light_menu == '2':
                            if self['res']['avDF'] < 0:
                                (dfpass, udfpass, avDF, uDF) = ('', '', 'N/A', 'N/A')
                            else:
                                dfpass = '(FAIL DF < 2)' if self['res']['avDF'] < 2 else '(PASS DF >= 2)'
                                udfpass = '(FAIL UDF < 0.4)' if self['res']['ratioDF'] < 0.4 else '(PASS UDF >= 0.4)'
                                avDF = self['res']['avDF']
                                uDF = self['res']['ratioDF']

                            row = layout.row()
                            row.label(text="Average DF: {} {}".format(avDF, dfpass))
                            row = layout.row()
                            row.label(text="Uniformity: {} {}".format(uDF, udfpass))

                    elif self.metric == '2' and self.probe_menu != 'None':
                        newrow(layout, 'Wind speed', self, "ws")

                        if self['res']['pressure']:
                            if self.zone_menu == 'All':
                                newrow(layout, 'Metric', self, "probe_menu")

                        if self.zone_menu == 'All':
                            for z in self['res'][self.probe_menu]:
                                row = layout.row()
                                row.label(text="{}: {}".format(z, self['res'][self.probe_menu][z]))
                        else:
                            if self['res']:
                                for m in self['res']:
                                    if self['res'][m].get(self.zone_menu):
                                        row = layout.row()
                                        row.label(text="{} {}: {}".format(self.zone_menu, m, self['res'][m][self.zone_menu]))

                    elif self.metric == '3':
                        if self['res']['ec'] and self.em_menu in ('Object', 'Surface', 'Zone') and self.frame_menu != 'All':
                            # newrow(layout, 'Type', self, 'riba_menu')
                            newrow(layout, 'Timespan:', self, "ec_years")
                            row = layout.row()

                            if self.zone_menu == 'All':
                                row.label(text='Total EC: {:.2f} kgCO2e'.format(float(self['res']['ec']['All'])))

                                if self['res']['ecm2']:
                                    row = layout.row()
                                    row.label(text='Total EC/m2: {:.2f} kgCO2e/m2'.format(float(self['res']['ecm2']['All'])))
                                if self['res']['ecm2y']:
                                    row = layout.row()
                                    row.label(text='Total EC/m2/y: {:.2f} kgCO2e/m2/y'.format(float(self['res']['ecm2y']['All'])))
                                if self['res']['vol']:
                                    row = layout.row()
                                    row.label(text='Total volume: {:.2f} m3'.format(float(self['res']['vol']['All'])))
                                if self['res']['area']:
                                    row = layout.row()
                                    row.label(text='Total area: {:.2f} m2'.format(self['res']['area'][self.zone_menu]))
                            else:
                                row.label(text='{} EC: {:.2f} kgCO2e'.format(self.em_menu, self['res']['ec'][self.zone_menu]))

                                if self['res']['ecm2'] and self.zone_menu != 'None':
                                    row = layout.row()
                                    row.label(text='{} EC/m2: {:.2f} kgCO2e/m2 (floor area)'.format(self.em_menu, self['res']['ecm2'][self.zone_menu]))
                                if self['res']['ecm2y'] and self.zone_menu != 'None':
                                    row = layout.row()
                                    row.label(text='{} EC/m2/y: {:.2f} kgCO2e/m2/y (floor area)'.format(self.em_menu, self['res']['ecm2y'][self.zone_menu]))
                                if self['res']['vol'] and self.zone_menu != 'None':
                                    row = layout.row()
                                    row.label(text='{} volume: {:.2f} m3'.format(self.em_menu, self['res']['vol'][self.zone_menu]))
                                if self.em_menu != 'Object':
                                    if self['res']['area'] and self.zone_menu != 'None':
                                        row = layout.row()
                                        row.label(text='{} area: {:.2f} m2'.format(self.em_menu, self['res']['area'][self.zone_menu]))

                    elif self.metric == '4':
                        newrow(layout, 'Comfort type:', self, "com_menu")

                        if self.com_menu == '0':
                            if self['res']['ooh'] >= 0:
                                newrow(layout, 'Occupied', self, "occ")

                            (r1, r2) = ('ooh', 'ooh2') if self.occ else ('oh', 'oh2')

                            if self['res'][r2] >= 0:
                                row = layout.row()
                                row.label(text="Time above 28degC: {:.1f}%".format(self['res'][r2]))

                            if self['res'][r1] >= 0:
                                row = layout.row()
                                row.label(text="Time between 25-28degC: {:.1f}%".format(self['res'][r1]))

                            elif self.frame_menu != 'All':
                                row = layout.row()
                                row.label(text="Overheating data not available")

                    elif self.metric == '5':
                        newrow(layout, 'IAQ type:', self, "iaq_menu")

                        if self.iaq_menu == '0':
                            row = layout.row()
                            if self['res']['co2'] >= 0:
                                row.label(text="Time CO2 below 900ppm: {:.1f}%".format(self['res']['co2']))
                            else:
                                row.label(text="CO2 data not available")

                    elif self.metric == '6':
                        if self['res'].get('wlc'):
                            newrow(layout, 'Timespan:', self, "ec_years")
                            newrow(layout, 'Heating type:', self, 'heat_type')

                            if self.heat_type == '0':
                                newrow(layout, 'Heating efficiency:', self, 'gas_eff')

                                if self['res']['ckwh']:
                                    newrow(layout, 'Cooling COP:', self, 'ac_cop')
                                    newrow(layout, 'Carbon factor:', self, 'carb_fac')
                                    newrow(layout, 'Annual change:', self, 'carb_annc')
                            else:
                                newrow(layout, 'Heating COP:', self, 'elec_cop')
                                newrow(layout, 'Hot water COP:', self, 'hw_cop')

                                if self['res']['ckwh']:
                                    newrow(layout, 'Cooling EER:', self, 'ac_cop')

                                newrow(layout, 'Carbon factor:', self, 'carb_fac')
                                newrow(layout, 'Annual change:', self, 'carb_annc')

                            newrow(layout, 'Unregulated', self, 'mod')
                            newrow(layout, 'Hot water:', self, 'hwmod')

                            if self['res']['owlc']:
                                row = layout.row()
                                row.label(text="Gross operational carbon: {:.1f} kgCO2e".format(self['res']['owlc']))
                            if self['res']['ofwlc']:
                                row = layout.row()
                                row.label(text="Carbon offset: {:.1f} kgCO2e".format(self['res']['ofwlc']))
                            if self['res']['owlc']:
                                row = layout.row()
                                row.label(text="Net operational carbon: {:.1f} kgCO2e".format(self['res']['owlc'] - self['res']['ofwlc']))
                            if self['res']['ecwlc']:
                                row = layout.row()
                                row.label(text="Total embodied carbon: {:.1f} kgCO2e".format(self['res']['ecwlc']))
                            if self['res']['wlc']:
                                row = layout.row()
                                row.label(text="Whole-life carbon: {:.1f} kgCO2e".format(self['res']['wlc']))

                            if self.frame_menu != 'All':
                                row = layout.row()
                                row.operator('node.vi_info', text='Infographic')

                    if self.metric == '1':
                        if self.light_menu == '2':
                            if self['res']['ratioDF'] >= 0:
                                row = layout.row()
                                row.operator('node.vi_info', text='Infographic')
                        elif self.light_menu == '1':
                            if self['res']['sda'] >= 0:
                                row = layout.row()
                                row.operator('node.vi_info', text='Infographic')
                    elif self.metric == '3':
                        if self['res']['ec'] and self.frame_menu != 'All':
                            row = layout.row()
                            row.operator('node.ec_pie', text='Pie chart')
                        # elif self.frame_menu == 'All':
                        #     row = layout.row()
                        #     row.operator('node.ec_line', text='Line chart')
                    elif self.metric == '4' and self.frame_menu == 'All':
                        if self['res'].get('alloh1s'):
                            row = layout.row()
                            row.operator('node.com_line', text='Line chart')

                    elif self.metric == '6' and self.frame_menu == 'All':
                        if self['res']['wl']:
                            row = layout.row()
                            row.operator('node.wlc_line', text='Line chart')


    def update(self):
        if self.inputs[0].links:
            self['rl'] = self.inputs[0].links[0].from_node['reslists']

            if self['rl'] and len(self['rl'][0]):
                frames = list(dict.fromkeys([z[0] for z in self['rl']]))
                self['frames'] =  [(f, f, 'Frame') for f in frames]

                if self.metric == '3':
                    # if  == 'Surface':
                    znames = sorted(list(dict.fromkeys([z[2] for z in self['rl'] if z[1] == 'Embodied carbon' and z[3] == f"{self.em_menu} EC (kgCO2e/y)"])))
                    # elif self.em_menu == 'Zone':
                    #     znames = sorted(list(dict.fromkeys([z[2] for z in self['rl'] if z[1] == 'Embodied carbon' and z[3] == "Zone EC (kgCO2e/y)"])))
                    # elif self.em_menu == 'Object':
                    #     znames = sorted(list(dict.fromkeys([z[2] for z in self['rl'] if z[1] == 'Embodied carbon' and z[3] == "Object EC (kgCO2e/y)"])))

                    if any ([z[3] == "Total EC (kgCO2e/y)" for z in self['rl']]):
                        self['znames'] = [(zn, zn, 'Zone name') for zn in znames] + [('All', 'All', 'All entities')]
                    else:
                        self['znames'] = [(zn, zn, 'Zone name') for zn in znames]

                elif self.metric == '0':
                    znames = sorted(list(dict.fromkeys([z[2] for z in self['rl'] if z[1] == 'Zone temporal'])))
                    self['znames'] = [(zn, zn, 'Zone name') for zn in znames] + [('All', 'All', 'All entities')]
                else:
                    znames = sorted(list(dict.fromkeys([z[2] for z in self['rl'] if z[1] in ('Zone spatial', 'Zone temporal')])))
                    self['znames'] = [(zn, zn, 'Zone name') for zn in znames]

                self.inputs[0].links[0].from_node.new_res = 0

            else:
                self['rl'] = []
                self['frames'] = [('None', 'None', 'None')]
                self['znames'] = [('None', 'None', 'None')]
        else:
            self['rl'] = []
            self['frames'] = [('None', 'None', 'None')]
            self['znames'] = [('None', 'None', 'None')]

        if self.frame_menu == '' or self.frame_menu not in [sf[0] for sf in self['frames']]:
            self.frame_menu = self['frames'][0][0]

        if self.zone_menu == '' or self.zone_menu not in [szn[0] for szn in self['znames']]:
            self.zone_menu = self['znames'][0][0]

        # if self.em_menu == '' or self.em_menu not in ('Object', 'Surface', 'Zone'):
            # self.em_menu = 'None'
        # if self.probe_menu == '' or self.probe_menu not in [spr[0] for spr in self['probes']]:
        #     self.probe_menu = self['probes'][0][0]

        self.res_update()

    def res_update(self):
        self['res'] = {}
        if self.metric == '0' and bpy.data.collections.get('EnVi Geometry') and 'Heating (W)' in [metric[0][3] for metric in zip(self['rl'])]:
            self['res']['pvkwh'] = 0
            self['res']['hkwh'] = 0
            self['res']['ahkwh'] = 0
            self['res']['ckwh'] = 0
            self['res']['ackwh'] = 0
            geo_coll = bpy.data.collections['EnVi Geometry']
            heat_mod = self.gas_eff/100 if self.heat_type == '0' else self.elec_cop
            hw_mod = self.gas_eff/100 if self.heat_type == '0' else self.hw_cop

            if self.zone_menu == 'All':
                if geo_coll.vi_params['enparams'].get('floorarea'):
                    self['res']['fa'] = geo_coll.vi_params['enparams']['floorarea'][str(self.frame_menu)]

            elif self.zone_menu != 'None':
                if geo_coll.children[self.zone_menu].vi_params['enparams'].get('floorarea'):
                    self['res']['fa'] = geo_coll.children[self.zone_menu].vi_params['enparams']['floorarea'][str(self.frame_menu)]

            if self['res']['fa'] > 13.9:
                occ = 1 + 1.76*(1 - math.exp(-0.000349 * (self['res']['fa'] -13.9)**2)) + 0.0013 * (self['res']['fa'] - 13.9)
            else:
                occ = 1

            Vda = 25 * occ + 36
            md = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
            ff = (1.10, 1.06, 1.02, 0.98, 0.94, 0.90, 0.90, 0.94, 0.98, 1.02, 1.06, 1.10, 1.00)
            dtm = (41.2, 41.4, 40.1, 37.6, 36.4, 33.9, 30.4, 33.4, 33.5, 36.3, 39.4, 39.9, 37.0)
            self['res']['wkwh'] = 1.15 * sum([4.18/3600 * Vda * z[0] * z[1] * z[2] / hw_mod for z in zip(md, ff, dtm)])

            if self.energy_menu == '0':
                for r in self['rl']:
                    if r[0] == self.frame_menu and self.zone_menu == 'All':
                        if r[3] == 'PV power (W)' and r[2] == 'All':
                            self['res']['pvkwh'] += sum(float(p) for p in r[4].split()) * 0.001
                        elif r[3] == 'Heating (W)':
                            self['res']['hkwh'] += sum(float(p) for p in r[4].split()) * 0.001
                        elif r[3] == 'Air heating (W)':
                            self['res']['ahkwh'] += sum(float(p) for p in r[4].split()) * 0.001
                        elif r[3] == 'Cooling (W)':
                            self['res']['ckwh'] += sum(float(p) for p in r[4].split()) * 0.001
                        elif r[3] == 'Air cooling (W)':
                            self['res']['ackwh'] += sum(float(p) for p in r[4].split()) * 0.001

                    elif r[0] == self.frame_menu:
                        if r[2] == self.zone_menu:
                            if r[3] == 'Heating (W)':
                                self['res']['hkwh'] = sum(float(p) for p in r[4].split()) * 0.001
                            if r[3] == 'Air heating (W)':
                                self['res']['ahkwh'] = sum(float(p) for p in r[4].split()) * 0.001
                            elif r[3] == 'Cooling (W)':
                                self['res']['ckwh'] = sum(float(p) for p in r[4].split()) * 0.001
                            elif r[3] == 'Air cooling (W)':
                                self['res']['ackwh'] = sum(float(p) for p in r[4].split()) * 0.001
                        elif r[1] == 'Power' and '_'.join(r[2].split('_')[:-1]) == self.zone_menu and r[3] == 'PV power (W)':
                                self['res']['pvkwh'] += sum(float(p) for p in r[4].split()) * 0.001

                self['res']['hkwh'] = self['res']['hkwh'] / heat_mod + self['res']['ahkwh'] - self['res']['hkwh'] if self['res']['ahkwh'] > self['res']['hkwh'] else self['res']['hkwh'] / heat_mod
                self['res']['ckwh'] = self['res']['ckwh'] / self.ac_cop + self['res']['ackwh'] - self['res']['ckwh'] if self['res']['ackwh'] > self['res']['ckwh'] else self['res']['ckwh'] / self.ac_cop
                self['res']['totkwh'] = self['res']['hkwh'] + self['res']['ckwh'] + self.mod * self['res']['fa'] + self['res']['wkwh'] - self['res']['pvkwh']
                self['res']['ECF'] = 0.42*(54 + self['res']['totkwh'] * 0.1319)/(self['res']['fa'] + 45)
                self['res']['EPC'] = 100 - 13.95 * self['res']['ECF'] if self['res']['ECF'] < 3.5 else 117 - 121 * math.log10(self['res']['ECF'])
                epcletts = ('A', 'B', 'C', 'D', 'E', 'F','G')
                epcnum = (92, 81, 69, 55, 39, 21, 1)

                for ei, en in enumerate(epcnum):
                    if self['res']['EPC'] < epcnum[-1]:
                        self['res']['EPCL'] = 'U'
                    elif self['res']['EPC'] > en:
                        self['res']['EPCL'] = epcletts[ei]
                        break

            elif self.energy_menu == '1':
                for r in self['rl']:
                    if r[0] == self.frame_menu and self.zone_menu == 'All':
                        if r[3] == 'PV power (W)' and r[2] == 'All':
                            self['res']['pvkwh'] += sum(float(p) for p in r[4].split()) * 0.001
                        elif r[3] == 'Heating (W)':
                            self['res']['hkwh'] += sum(float(p) for p in r[4].split()) * 0.001
                        elif r[3] == 'Air heating (W)':
                            self['res']['ahkwh'] += sum(float(p) for p in r[4].split()) * 0.001
                        elif r[3] == 'Cooling (W)':
                            self['res']['ckwh'] += sum(float(p) for p in r[4].split()) * 0.001
                        elif r[3] == 'Air cooling (W)':
                            self['res']['ackwh'] += sum(float(p) for p in r[4].split()) * 0.001

                    elif r[0] == self.frame_menu:
                        if r[2] == self.zone_menu :
                            if r[3] == 'Heating (W)':
                                self['res']['hkwh'] = sum(float(p) for p in r[4].split()) * 0.001
                            elif r[3] == 'Air heating (W)':
                                self['res']['ahkwh'] = sum(float(p) for p in r[4].split()) * 0.001
                            elif r[3] == 'Cooling (W)':
                                self['res']['ckwh'] = sum(float(p) for p in r[4].split()) * 0.001
                            elif r[3] == 'Air cooling (W)':
                                self['res']['ackwh'] = sum(float(p) for p in r[4].split()) * 0.001
                        elif r[2] != 'All' and r[1] == 'Power' and 'EN_' + r[2].split('_')[1] == self.zone_menu and r[3] == 'PV power (W)':
                                self['res']['pvkwh'] += sum(float(p) for p in r[4].split()) * 0.001

                self['res']['hkwh'] = self['res']['hkwh'] / heat_mod + self['res']['ahkwh'] - self['res']['hkwh'] if self['res']['ahkwh'] > self['res']['hkwh'] else self['res']['hkwh'] / heat_mod
                self['res']['ckwh'] = self['res']['ckwh'] / self.ac_cop + self['res']['ackwh'] - self['res']['ckwh'] if self['res']['ackwh'] > self['res']['ckwh'] else self['res']['ckwh'] / self.ac_cop
                self['res']['totkwh'] = (self['res']['hkwh'] + self['res']['ckwh'] + self.mod * self['res']['fa'])
                self['res']['netkwh'] = (self['res']['hkwh'] + self['res']['ckwh'] - self['res']['pvkwh'] + self.mod * self['res']['fa'])

            elif self.energy_menu == '2':
                for r in self['rl']:
                    if r[0] == self.frame_menu and self.zone_menu == 'All':
                        if r[3] == 'PV power (W)' and r[2] == 'All':
                            self['res']['pvkwh'] += sum(float(p) for p in r[4].split()) * 0.001
                        elif r[3] == 'Heating (W)':
                            self['res']['hkwh'] += sum(float(p) for p in r[4].split()) * 0.001
                        elif r[3] == 'Air heating (W)':
                            self['res']['ahkwh'] += sum(float(p) for p in r[4].split()) * 0.001
                        elif r[3] == 'Cooling (W)':
                            self['res']['ckwh'] += sum(float(p) for p in r[4].split()) * 0.001
                        elif r[3] == 'Air cooling (W)':
                            self['res']['ackwh'] += sum(float(p) for p in r[4].split()) * 0.001

                    elif r[0] == self.frame_menu:
                        if r[2] == self.zone_menu :
                            if r[3] == 'Heating (W)':
                                self['res']['hkwh'] = sum(float(p) for p in r[4].split()) * 0.001
                            elif r[3] == 'Air heating (W)':
                                self['res']['ahkwh'] = sum(float(p) for p in r[4].split()) * 0.001
                            elif r[3] == 'Cooling (W)':
                                self['res']['ckwh'] = sum(float(p) for p in r[4].split()) * 0.001
                            elif r[3] == 'Air cooling (W)':
                                self['res']['ackwh'] = sum(float(p) for p in r[4].split()) * 0.001
                        elif r[2] != 'All' and r[1] == 'Power' and 'EN_' + r[2].split('_')[1] == self.zone_menu and r[3] == 'PV power (W)':
                                self['res']['pvkwh'] += sum(float(p) for p in r[4].split()) * 0.001

                self['res']['hkwh'] = self['res']['hkwh'] / heat_mod + self['res']['ahkwh'] - self['res']['hkwh'] if self['res']['ahkwh'] > self['res']['hkwh'] else self['res']['hkwh'] / heat_mod

        elif self.metric == '1':
            self['res']['avDF'] = -1
            self['res']['ratioDF'] = -1
            self['res']['ase'] = -1
            self['res']['sda'] = -1
            self['res']['sdaarea'] = -1
            self['res']['asepass'] = -1
            self['res']['sdapass'] = -1
            self['res']['sv'] = -1
            self['res']['auto'] = -1
            self['res']['o1'] = -1
            self['res']['areaDF'] = -1
            self['res']['minDF'] = -1
            mir = 0.7 if self.g_roof else 0.3

            if self.light_menu == '0':
                if self.breeam_menu == '0':
                    mDFs = (2, )
                    mAs = (0.8, )
                    sdah = 2000
                    sdal = 300
                    sdaml = 90
                    cred = 2

                    if self.breeam_edumenu == '1':
                        mAs = (0.6, 0.8)
                    #     cred = 1

                elif self.breeam_menu == '1':
                    cred = 2
                    mAs = (0.8, )
                    mDFs = (2, )
                    sdah = 2650
                    sdal = 300
                    sdaml = 90

                    if self.breeam_healthmenu == '1':
                        mDF = (3, 2)

                elif self.breeam_menu == '2':
                    cred = 1
                    mDFs = (2, )
                    mAs = (0.8, )
                    msdaA = 1
                    sdah = 3450
                    sdal = 100
                    sdaml = 30

                    if self.breeam_multimenu == '2':
                        sdah = 2650
                        sdal = 200
                        sdaml = 60
                        msdaA = 0.8

                elif self.breeam_menu == '3':
                    cred = 1
                    mDFs = (2, )
                    mAs = (0.8, )

                    if self.breeam_retailmenu == '0':
                        mAs = (0.35, )

                elif self.breeam_menu == '4':
                    cred = 2

                    if self.breeam_prisonmenu == '0':
                        mAs = (0.8, )
                        mDFs = (1.5, )
                    elif self.breeam_prisonmenu == '1':
                        mAs = (0.8, )
                        mDFs = (3, )
                    elif self.breeam_prisonmenu == '2':
                        mAs = (0.8, )
                        mDFs = (3, )
                    elif self.breeam_prisonmenu == '3':
                        mDFs = (2, )
                        mAs = (0.8, )

                elif self.breeam_menu == '5':
                    cred = 2
                    mDFs = (2, )
                    mAs = (0.8, )

                elif self.breeam_menu == '6':
                    cred = 2
                    mDFs = (2, )
                    mAs = (0.8, )

                elif self.breeam_menu == '7':
                    cred = 1
                    mDFs = (2, )
                    mAs = (0.8, )

                for r in self['rl']:
                    if r[0] == self.frame_menu:
                        if r[2] == self.zone_menu:
                            if r[3] == 'Areas (m2)':
                                dfareas = array([float(p) for p in r[4].split()])
                            elif r[3] == 'DF (%)':
                                df = array([float(p) for p in r[4].split()])
                            elif r[3] == 'Spatial Daylight Autonomy (% area)':
                                sareas = array([float(p) for p in r[4].split()])

                try:
                    tdf = stack((df, dfareas), axis=1)
                    stdf = tdf[tdf[:, 0].argsort()][::-1]
                    aDF = stdf[0][0]
                    rarea = stdf[0][1]
                    ir = stdf[0][0]/(nsum(df)/len(df))
                    i = 1

                    for mDF in mDFs:
                        for mA in mAs:
                            while (aDF >= mDF and i < len(df) and ir >= mir):
                                aDF = nsum(stdf[0: i + 1, 0] * dfareas[0: i + 1])/nsum(dfareas[0: i + 1])
                                minDF = stdf[i][0]
                                ir = minDF/aDF
                                rarea += stdf[i][1]
                                if rarea/nsum(dfareas) > mA:
                                    break
                                i += 1

                    self['res']['avDF'] = round(aDF, 2)
                    self['res']['areaDF'] = round(100 * rarea/nsum(dfareas), 2)
                    self['res']['ratioDF'] = round(ir, 2)
                    self['res']['minDF'] = round(minDF, 2)

                except Exception:
                    pass

                try:
                    vsareas = sareas[where(sareas > mA * 100)]
                    self['res']['sda'] = len(vsareas)

                    if self['res']['sda'] > sdah:
                        credits = cred

                except Exception:
                    pass

            elif self.light_menu == '2':
                for r in self['rl']:
                    if r[0] == self.frame_menu:
                        if r[2] == self.zone_menu:
                            if r[3] == 'Areas (m2)':
                                dfareas = array([float(p) for p in r[4].split()])
                            elif r[3] == 'DF (%)':
                                df = array([float(p) for p in r[4].split()])

                try:
                    self['res']['avDF'] = round(sum(df * dfareas)/sum(dfareas), 2)
                    self['res']['ratioDF'] = round(min(df)/self['res']['avDF'], 2)
                except Exception:
                    pass

            elif self.light_menu == '1':
                if 'Annual Sunlight Exposure (% area)' in [r[3] for r in self['rl']]:
                    for r in self['rl']:
                        if r[0] == self.frame_menu:
                            if r[2] == self.zone_menu:
                                if self.zone_menu in bpy.context.scene.objects:
                                    res_ob = bpy.context.scene.objects[self.zone_menu]

                                    if res_ob.vi_params['livires'].get('totarea{}'.format(self.frame_menu)):
                                        self['res']['totarea'] = res_ob.vi_params['livires']['totarea{}'.format(self.frame_menu)]
                                        self['res']['svarea'] = res_ob.vi_params['livires']['svarea{}'.format(self.frame_menu)]

                                        if r[3] == 'Annual Sunlight Exposure (% area)':
                                            self['res']['ase'] = 100 * res_ob.vi_params['livires']['ase{}'.format(self.frame_menu)]
                                            self['res']['asepass'] = 10
                                        elif r[3] == 'Spatial Daylight Autonomy (% area)':
                                            self['res']['sda'] = 100 * res_ob.vi_params['livires']['sda{}'.format(self.frame_menu)]
                                        elif r[3] == 'Spatial Daylight Autonomy (% perimeter area)':
                                            self['res']['sdapa'] = 100 * res_ob.vi_params['livires']['sdapa{}'.format(self.frame_menu)]
                                        elif r[3] == 'UDI-a Area (%)':
                                            udiaareas = array([float(p) for p in r[4].split()])
                                            im = self.inputs[0].links[0].from_node['coptions']['times'].index('20/03/15 09:00:00')
                                            ie = self.inputs[0].links[0].from_node['coptions']['times'].index('20/03/15 15:00:00')
                                            self['res']['udiam'] = udiaareas[im]
                                            self['res']['udiae'] = udiaareas[ie]

                                        self['res']['sv'] = 100 * res_ob.vi_params['livires']['svarea{}'.format(self.frame_menu)]/res_ob.vi_params['livires']['totarea{}'.format(self.frame_menu)]

        elif self.metric == '2':
            self['res']['pressure'] = {}
            self['res']['speed'] = {}
            self['res']['temperature'] = {}
            self['res']['xvelocity'] = {}
            self['res']['zvelocity'] = {}
            self['res']['yvelocity'] = {}
            self['res']['wpc'] = {}
            znames = set([z[2] for z in self['rl'] if z[1] == 'Zone spatial'])

            for zn in znames:
                for r in self['rl']:
                    if r[2] == zn:
                        if r[3] == 'Pressure':
                            self['res']['pressure'][zn] = float(r[4].split()[-1])
                            self['res']['wpc'][zn] = round((float(r[4].split()[-1]) - bpy.context.scene.vi_params['flparams']['pref'])/(0.5*1.225*(self.ws**2)), 3)
                        elif r[3] == 'Speed':
                            self['res']['speed'][zn] = float(r[4].split()[-1])
                        elif r[3] == 'Temperature':
                            self['res']['temperature'][zn] = float(r[4].split()[-1])

        elif self.metric == '3':
            self['res']['ec'] = {}
            self['res']['ecm2'] = {}
            self['res']['ecm2y'] = {}
            self['res']['vol'] = {}
            self['res']['area'] = {}

            for r in self['rl']:
                if r[0] == self.frame_menu and self.frame_menu != 'All':
                    if self.em_menu == 'Object':
                        if r[3] == 'Object EC (kgCO2e/y)':
                            self['res']['ec'][r[2]] = float(r[4]) * self.ec_years
                        elif r[3] == 'Object EC (kgCO2e/m2/y)':
                            self['res']['ecm2'][r[2]] = float(r[4]) * self.ec_years
                            self['res']['ecm2y'][r[2]] = float(r[4])
                        elif r[3] == 'Object volume (m3)':
                            self['res']['vol'][r[2]] = float(r[4])
                    elif self.em_menu == 'Surface':
                        if r[3] == 'Surface EC (kgCO2e/y)':
                            self['res']['ec'][r[2]] = float(r[4]) * self.ec_years
                        elif r[3] == 'Surface EC (kgCO2e/m2/y)':
                            self['res']['ecm2'][r[2]] = float(r[4]) * self.ec_years
                            self['res']['ecm2y'][r[2]] = float(r[4])
                        elif r[3] == 'Surface volume (m3)':
                            self['res']['vol'][r[2]] = float(r[4])
                        elif r[3] == 'Surface area (m2)':
                            self['res']['area'][r[2]] = float(r[4])
                    elif self.em_menu == 'Zone':
                        if r[3] == 'Zone EC (kgCO2e/y)':
                            self['res']['ec'][r[2]] = float(r[4]) * self.ec_years
                        elif r[3] == 'Zone EC (kgCO2e/m2/y)':
                            self['res']['ecm2'][r[2]] = float(r[4]) * self.ec_years
                            self['res']['ecm2y'][r[2]] = float(r[4])
                        elif r[3] == 'Surface area (m2)':
                            self['res']['area'][r[2]] = float(r[4])

            self['res']['ec']['All'] = sum([self['res']['ec'][ec] for ec in self['res']['ec']])
            self['res']['ecm2']['All'] = sum([self['res']['ecm2'][ec] for ec in self['res']['ecm2']])
            self['res']['ecm2y']['All'] = sum([self['res']['ecm2y'][ec] for ec in self['res']['ecm2y']])
                    # if r[2] == 'All' and r[3] == f'{self.em_menu} EC (kgCO2e/y)':
                    #     self['res']['ec']['All'] = float(r[4]) * self.ec_years
                    # elif r[2] == 'All' and r[3] == f'{self.em_menu} EC (kgCO2e/m2/y)':
                    #     self['res']['ecm2']['All'] = float(r[4]) * self.ec_years
                    #     self['res']['ecm2y']['All'] = float(r[4])
                    # elif r[2] == 'All' and r[3] == f'{self.em_menu} volume (m3)':
                    #     self['res']['vol']['All'] = float(r[4])
                    # elif r[2] == 'All' and r[3] == 'Object surface area (m2)':
                    #     if self.em_menu != 'Object':
                    #         self['res']['area']['All'] = float(r[4])

                # elif self.frame_menu == 'All':
                #     if r[0] == 'All' and
                #     self['res'][ec]

        elif self.metric == '4':
            self['res']['oh'] = -1
            self['res']['oh2'] = -1
            self['res']['ooh'] = -1
            self['res']['ooh2'] = -1
            temps = array([])
            occs = array([])
            alltemps = []
            alloccs = []

            if self.frame_menu != 'All':
                for r in self['rl']:
                    if r[0] == self.frame_menu:
                        if r[2] == self.zone_menu:
                            if r[3] == 'Temperature (degC)':
                                temps = array([float(p) for p in r[4].split()])
                            elif r[3] == 'Occupancy':
                                occs = array([float(p) for p in r[4].split()])

                if len(temps):
                    self['res']['oh'] = 100 * len(where(temps > 25)[0])/len(temps)
                    self['res']['oh2'] = 100 * len(where(temps > 28)[0])/len(temps)

                    if len(occs):
                        self['res']['ooh'] = 100 * len(where(where(occs > 0, temps, 0) > 25)[0])/len(where(occs > 0)[0])
                        self['res']['ooh2'] = 100 * len(where(where(occs > 0, temps, 0) > 28)[0])/len(where(occs > 0)[0])
            else:
                for r in self['rl']:
                    if r[2] == self.zone_menu:
                        if r[3] == 'Temperature (degC)':
                            alltemps.append(array([float(p) for p in r[4].split()]))
                        elif r[3] == 'Occupancy':
                            alloccs.append(array([float(p) for p in r[4].split()]))

                if alltemps:
                    self['res']['alloh1s'] = [100 * len(where(temps > 25)[0])/len(temps) for temps in alltemps]
                    self['res']['alloh2s'] = [100 * len(where(temps > 28)[0])/len(temps) for temps in alltemps]

                    if alloccs:
                        self['res']['allooh1s'] = [100 * len(where(where(alloccs[ti] > 0, temps, 0) > 25)[0])/len(where(alloccs[ti] > 0)[0]) for ti, temps in enumerate(alltemps)]
                        self['res']['allooh2s'] = [100 * len(where(where(alloccs[ti] > 0, temps, 0) > 28)[0])/len(where(alloccs[ti] > 0)[0]) for ti, temps in enumerate(alltemps)]
                # print(alloh1s, alloh2s)

        elif self.metric == '5':
            self['res']['co2'] = -1

            for r in self['rl']:
                if r[0] == self.frame_menu:
                    if r[2] == self.zone_menu:
                        if r[3] == 'CO2 (ppm)':
                            if self.iaq_menu == '0':
                                co2s = array([float(p) for p in r[4].split()])
                                self['res']['co2'] = 100 * len(where(co2s < 900)[0])/len(co2s)

        elif self.metric == '6':
            self['res']['wlc'] = 0
            self['res']['owlc'] = 0
            self['res']['ecwlc'] = 0
            self['res']['wlc'] = 0
            self['res']['ofwlc'] = 0
            self['res']['ckwh'] = 0

            if bpy.data.collections.get('EnVi Geometry') and self.zone_menu != 'None':
                hours = 8760
                ofcs = {}
                owlcs, ecwlcs, ofwlcs, wlcs, reslists = [], [], [], [], []
                cop = 1/self.elec_cop if self.heat_type == '1' else 1/(self.gas_eff * 0.01)
                hw_cop = 1/self.hw_cop if self.heat_type == '1' else 1/(self.gas_eff * 0.01)
                geo_coll = bpy.data.collections['EnVi Geometry']
                frames = [f[0] for f in self['frames'] if f[0] != 'All']

                if self.zone_menu == 'All':
                    if geo_coll.vi_params['enparams'].get('floorarea'):
                        self['res']['fa'] = list(geo_coll.vi_params['enparams']['floorarea'].values())

                elif self.zone_menu != '':
                    if geo_coll.children[self.zone_menu].vi_params['enparams'].get('floorarea'):
                        self['res']['fa'] = list(geo_coll.children[self.zone_menu].vi_params['enparams']['floorarea'].values())

                for frame in frames:
                    pv_kwh = 0
                    heat_kwh = 0
                    aheat_kwh = 0
                    cool_kwh = 0
                    acool_kwh = 0
                    airheat_kwh = 0
                    aairheat_kwh = 0
                    aircool_kwh = 0
                    aaircool_kwh = 0
                    ec_kgco2e = 0

                    if self.frame_menu and self.zone_menu:
                        for r in self['rl']:
                            if r[0] == frame:
                                if r[2] == self.zone_menu:
                                    if r[3] == 'Heating (W)':
                                        heat_kwh = sum([float(h) for h in r[4].split()]) * 0.001
                                        aheat_kwh += heat_kwh
                                        hours = len(r[4].split())
                                    if r[3] == 'Air heating (W)':
                                        airheat_kwh = sum([float(h) for h in r[4].split()]) * 0.001
                                        aairheat_kwh += airheat_kwh
                                    if r[3] == 'Cooling (W)':
                                        cool_kwh = sum([float(h) for h in r[4].split()]) * 0.001
                                        acool_kwh += cool_kwh
                                    if r[3] == 'Air cooling (W)':
                                        aircool_kwh = sum([float(h) for h in r[4].split()]) * 0.001
                                        aaircool_kwh += aaircool_kwh
                                    elif r[3] == '{} EC (kgCO2e/y)'.format(('Zone', 'Total')[self.zone_menu == 'All']):
                                        ec_kgco2e = float(r[4])
                                    elif r[3] == 'PV power (W)':
                                        pv_kwh += sum(float(p) for p in r[4].split()) * 0.001

                                elif r[1] == 'Power' and '_'.join(r[2].split('_')[:-1]) == self.zone_menu and r[3] == 'PV power (W)':
                                    pv_kwh += sum(float(p) for p in r[4].split()) * 0.001
                                elif r[1] == 'Zone temporal' and self.zone_menu == 'All' and r[3] == 'Heating (W)':
                                    aheat_kwh += sum(float(p) for p in r[4].split()) * 0.001
                                elif r[1] == 'Zone temporal' and self.zone_menu == 'All' and r[3] == 'Cooling (W)':
                                    acool_kwh += sum(float(p) for p in r[4].split()) * 0.001

                        heat_kwh = airheat_kwh if airheat_kwh > heat_kwh else heat_kwh
                        cool_kwh = aircool_kwh if aircool_kwh > cool_kwh else cool_kwh
                        aheat_kwh = aairheat_kwh if aairheat_kwh > aheat_kwh else aheat_kwh
                        acool_kwh = aaircool_kwh if aaircool_kwh > acool_kwh else acool_kwh
                        (heat_kwh, cool_kwh) = (aheat_kwh, acool_kwh) if self.zone_menu == 'All' else (heat_kwh, cool_kwh)
                        o_kwh = heat_kwh * 8760/hours * cop + cool_kwh * 8760/hours * self.ac_cop + (self.mod + self.hwmod * hw_cop) * self['res']['fa'][int(frame) - int(frames[0])]
                        owlc = o_kwh * sum([self.carb_fac * (1 + (self.carb_annc * 0.01))**i for i in range(self.ec_years)])
                        ecwlc = ec_kgco2e * self.ec_years
                        ofwlc = pv_kwh * sum([self.carb_fac * (1 + (self.carb_annc * 0.01))**i for i in range(self.ec_years)])
                        wlc = owlc + ecwlc - ofwlc

                        if self.frame_menu == frame:
                            self['res']['owlc'] = owlc
                            self['res']['ecwlc'] = ecwlc
                            self['res']['ofwlc'] = ofwlc
                            self['res']['wlc'] = wlc

                        owlcs.append(owlc)
                        ecwlcs.append(ecwlc)
                        ofwlcs.append(ofwlc)
                        wlcs.append(wlc)

                    reslists.append([str(frame), 'Carbon', self.zone_menu, 'Carbon offset (kgCO2e)', '{:.3f}'.format(ofwlc)])
                    reslists.append([str(frame), 'Carbon', self.zone_menu, 'Gross operational carbon (kgCO2e)', '{:.3f}'.format(owlc)])
                    reslists.append([str(frame), 'Carbon', self.zone_menu, 'Embodied carbon (kgCO2e)', '{:.3f}'.format(ecwlc)])
                    reslists.append([str(frame), 'Carbon', self.zone_menu, 'Net operational carbon (kgCO2e)', '{:.3f}'.format(owlc - ofwlc)])
                    reslists.append([str(frame), 'Carbon', self.zone_menu, 'Whole-life carbon (kgCO2e)', '{:.3f}'.format(wlc)])

                self['res']['ec'] = ecwlcs
                self['res']['of'] = ofwlcs
                self['res']['wl'] = wlcs
                self['res']['oc'] = owlcs
                self['res']['noc'] = [owlcs[i] - ofwlc for i, ofwlc in enumerate(ofwlcs)]
                reslists.append(['All', 'Frames', 'Frames', 'Frames', ' '.join(['{}'.format(f[0]) for f in self['frames']])])
                reslists.append(['All', 'Carbon', self.zone_menu, 'Carbon offset (kgCO2e)', ' '.join(['{:.3f}'.format(ofwlc) for ofwlc in ofwlcs])])
                reslists.append(['All', 'Carbon', self.zone_menu, 'Gross operational carbon (kgCO2e)', ' '.join(['{:.3f}'.format(owlc) for owlc in owlcs])])
                reslists.append(['All', 'Carbon', self.zone_menu, 'Embodied carbon (kgCO2e)', ' '.join(['{:.3f}'.format(ecwlc) for ecwlc in ecwlcs])])
                reslists.append(['All', 'Carbon', self.zone_menu, 'Net operational carbon (kgCO2e)', ' '.join(['{:.3f}'.format(owlcs[i] - ofwlc) for i, ofwlc in enumerate(ofwlcs)])])
                reslists.append(['All', 'Carbon', self.zone_menu, 'Whole-life carbon (kgCO2e)', ' '.join(['{:.3f}'.format(wlc) for wlc in wlcs])])
                self['reslists'] = reslists

    # def ret_reslists(self, zones):
    #     reslists = []

    #     for frame in [f[0] for f in self['frames']]:
    #         for zone in zones:
    #             for r in self['rl']:
    #                 if r[2] == zone:
    #                     if r[3] == 'Heating (W)':
    #                         heat_kwh = sum([float(h) for h in r[4].split()]) * 0.001
    #                         aheat_kwh += heat_kwh
    #                         hours = len(r[4].split())
    #                     if r[3] == 'Cooling (W)':
    #                         cool_kwh = sum([float(h) for h in r[4].split()]) * 0.001
    #                         acool_kwh += cool_kwh
    #                     elif r[3] == '{} EC (kgCO2e/y)'.format(('Zone', 'Total')[self.zone_menu == 'All']):
    #                         ec_kgco2e = float(r[4])
    #                     elif r[3] == 'PV power (W)':
    #                         pv_kwh += sum(float(p) for p in r[4].split()) * 0.001

    #                 elif r[1] == 'Power' and '_'.join(r[2].split('_')[:-1]) == zone and r[3] == 'PV power (W)':
    #                     pv_kwh += sum(float(p) for p in r[4].split()) * 0.001
    #                 elif r[1] == 'Power' and zone == 'All' and r[3] == 'PV power (W)':
    #                     pv_kwh += sum(float(p) for p in r[4].split()) * 0.001
    #                 elif r[1] == 'Zone temporal' and zone == 'All' and r[3] == 'Heating (W)':
    #                     aheat_kwh += sum(float(p) for p in r[4].split()) * 0.001
    #                 elif r[1] == 'Zone temporal' and zone == 'All' and r[3] == 'Cooling (W)':
    #                     acool_kwh += sum(float(p) for p in r[4].split()) * 0.001

    #             (heat_kwh, cool_kwh) = (aheat_kwh, acool_kwh) if self.zone_menu == 'All' else (heat_kwh, cool_kwh)
    #             o_kwh = heat_kwh * 8760/hours * cop + cool_kwh * 8760/hours * self.ac_cop + (self.mod + self.hwmod * hw_cop) * self['res']['fa']
    #             owlc = o_kwh * self.ec_years * self.carb_fac * (1 + (self.carb_annc * 0.01))**self.ec_years
    #             ecwlc = ec_kgco2e * self.ec_years
    #             ofwlc = pv_kwh * self.ec_years * self.carb_fac * (1 + (self.carb_annc * 0.01))**self.ec_years
    #             wlc = owlc + ecwlc - ofwlc

    #             owlcs.append(owlc)
    #             ecwlcs.append(ecwlc)
    #             ofwlcs.append(ofwlc)
    #             wlcs.append(wlc)
    #             # reslists.append(['All', 'Frames', 'Frames', 'Frames', ' '.join(['{}'.format(f[0]) for f in self['frames']])])
    #             # reslists.append(['All', 'Carbon', zone, 'Embodied carbon (kgCO2e)', ' '.join(['{:.3f}'.format()])
    #             # reslists.append(['All', 'Carbon', zone, 'Embodied carbon (kgCO2e)', ' '.join(['{:.3f}'.format()])])
    #             # reslists.append(['All', 'Carbon', zone, 'Embodied carbon (kgCO2e)', ])

    def ret_metrics(self):
        if self.inputs['Results in'].links:
            reslist = self.inputs['Results in'].links[0].from_node['reslists']


class No_CSV(Node, ViNodes):
    '''CSV Export Node'''
    bl_idname = 'No_CSV'
    bl_label = 'VI CSV Export'
    bl_icon = 'SPREADSHEET'

    animated: BoolProperty(name='', description='Animated results', default=0)

    def init(self, context):
        self.inputs.new('So_Vi_Res', 'Results in')

    def draw_buttons(self, context, layout):
        try:
            rl = self.inputs['Results in'].links[0].from_node['reslists']
            zrl = list(zip(*rl))

            if len(set(zrl[0])) > 1:
                newrow(layout, 'Animated:', self, 'animated')

            row = layout.row()
            row.operator('node.csvexport', text='Export CSV file')
        except Exception:
            pass

    def update(self):
        pass


class ViNodeCategory(NodeCategory):
    @classmethod
    def poll(cls, context):
        return context.space_data.tree_type == 'ViN'


class So_Anim(NodeSocket):
    '''Vi Animation socket'''
    bl_idname = 'So_Anim'
    bl_label = 'Animation socket'

    valid = ['Animation']

    def draw(self, context, layout, node, text):
        layout.label(text=text)

    def draw_color(self, context, node):
        return (1, 1.0, 0.45, 1.0)

    def ret_valid(self, node):
        return ['Animation']


class So_Vi_Loc(NodeSocket):
    '''Vi Location socket'''
    bl_idname = 'So_Vi_Loc'
    bl_label = 'Location socket'

    valid = ['Location']

    def draw(self, context, layout, node, text):
        layout.label(text=text)

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
        layout.label(text=text)

    def draw_color(self, context, node):
        return (0.3, 0.17, 0.07, 0.75)


class So_Li_Con(NodeSocket):
    '''Lighting context socket'''
    bl_idname = 'So_Li_Con'
    bl_label = 'Context'

    valid = ['LiVi Context', 'text']
    link_limit = 1

    def draw(self, context, layout, node, text):
        layout.label(text=text)

    def draw_color(self, context, node):
        return (1.0, 1.0, 0.0, 0.75)


class So_Text(NodeSocket):
    '''VI text socket'''
    bl_idname = 'So_Text'
    bl_label = 'VI text export'

    valid = ['text']
    link_limit = 1

    def draw(self, context, layout, node, text):
        layout.label(text=text)

    def draw_color(self, context, node):
        return (0.2, 1.0, 0.0, 0.75)


class So_Vi_Res(NodeSocket):
    '''Vi results socket'''
    bl_idname = 'So_Vi_Res'
    bl_label = 'VI results'

    valid = ['Vi Results']
    link_limit = 1

    def draw(self, context, layout, node, text):
        layout.label(text=self.bl_label)

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
        layout.label(text=text)

    def draw_color(self, context, node):
        return (0.5, 1.0, 0.0, 0.75)


class So_En_Geo(NodeSocket):
    '''EnVi geometry out socket'''
    bl_idname = 'So_En_Geo'
    bl_label = 'EnVi Geometry Socket'

    valid = ['EnVi Geometry']
    link_limit = 1

    def draw(self, context, layout, node, text):
        layout.label(text=text)

    def draw_color(self, context, node):
        return (0.0, 0.0, 1.0, 0.75)


class So_En_Con(NodeSocket):
    '''EnVi context socket'''
    bl_idname = 'So_En_Con'
    bl_label = 'EnVi context'

    valid = ['EnVi Context']
    link_limit = 1

    def draw(self, context, layout, node, text):
        layout.label(text=text)

    def draw_color(self, context, node):
        return (0.0, 1.0, 1.0, 0.75)


class So_En_Res(NodeSocket):
    '''Results socket'''
    bl_idname = 'So_En_Res'
    bl_label = 'Results axis'
    valid = ['Vi Results']

    def draw(self, context, layout, node, text):
        typedict = {"Time": [], "Frames": [], "Climate": ['climmenu'],
                    "Zone spatial": ("zonemenu", "zonermenu"), "Zone temporal": ("zonemenu", "zonermenu"),
                    "Embodied carbon": ("ecmenu", "ecrmenu"), "Linkage": ("linkmenu", "linkrmenu"),
                    "External": ("enmenu", "enrmenu"), "Position": ("posmenu", "posrmenu"),
                    "Camera": ("cammenu", "camrmenu"), "Power": ("powmenu", "powrmenu"),
                    "Probe": ("probemenu", "probermenu")}
        row = layout.row()

        if self.links and self.links[0].from_node.get('frames'):
            if len(self.links[0].from_node['frames']) > 1 and node.parametricmenu == '0':
                row.prop(self, "framemenu", text=text)
                row.prop(self, "rtypemenu")
            else:
                row.prop(self, "rtypemenu", text=text)

            if self.rtypemenu in typedict:
                for rtype in typedict[self.rtypemenu]:
                    row.prop(self, rtype)

            if self.node.timemenu in ('1', '2') and self.rtypemenu != 'Time' and node.parametricmenu == '0':
                row.prop(self, "statmenu")

            if self.rtypemenu != 'Time':
                row.prop(self, 'multfactor')
        else:
            layout.label(text=self.bl_label)

    def draw_color(self, context, node):
        return (0.0, 1.0, 0.0, 0.75)

# Openfoam nodes


class No_Flo_Case(Node, ViNodes):
    '''Openfoam case export node'''
    bl_idname = 'No_Flo_Case'
    bl_label = 'FloVi Case'
    bl_icon = 'FILE_FOLDER'

    def ret_params(self):
        return [str(x) for x in (self.scenario, self.buoyancy, self.parametric, self.frame_start, self.frame_end,
                                 self.dtime, self.etime, self.pnormval, self.pabsval, self.uval, self.tval, self.nutval, self.nutildaval,
                                 self.kval, self.epval, self.oval, self.presid, self.uresid, self.keoresid, self.aval, self.p_rghval,
                                 self.Gval, self.radmodel, self.solar, self.sun, self.comfort, self.clo, self.met, self.rh, self.age)]

    def nodeupdate(self, context):
        nodecolour(self, self['exportstate'] != self.ret_params())

        if self.scenario in ('2', '3') and self.buoyancy == 0:
            self.buoyancy = 1
            if self.scenario == '3' and self.age == 1:
                self.age = 0
        elif self.scenario in ('0', '1') and self.buoyancy == 1:
            self.buoyancy = 0

        context.scene.vi_params['flparams']['scenario'] = self.scenario

    scenario: EnumProperty(items=[('0', 'External flow', 'Wind induced flow'), ('1', 'Internal flow', 'Internal forced flow'), ('2', 'Forced convection', 'Forced convection'),
                                  ('3', 'Free convection', 'Free convection'), ('4', 'Custom', 'Custom scenario')], name='', description='Scenario type', default=0, update=nodeupdate)
    parametric: BoolProperty(name='', description='Parametric simulation', default=0, update=nodeupdate)
    frame_start: IntProperty(name="", description="Start frame", min=0, default=0, update=nodeupdate)
    frame_end: IntProperty(name="", description="End frame", min=0, default=0, update=nodeupdate)
    # solver: EnumProperty(name='', items=[('simpleFoam', 'SimpleFoam', 'SimpleFoam solver')], description='Solver selection', default='simpleFoam')
    # transience: EnumProperty(name='', items=[('0', 'Steady', 'Steady state simulation'),
    #                                              ('1', 'Transient', 'Transient simulation')], description='Transience selection', default='0', update=nodeupdate)

    # turbulence: EnumProperty(items=[('laminar', 'Laminar', 'Steady state turbulence solver'),
    #                                   ('kEpsilon', 'k-Epsilon', 'Transient laminar solver'),
    #                                   ('kOmega', 'k-Omega', 'Transient turbulence solver'),
    #                                   ('SpalartAllmaras', 'Spalart-Allmaras', 'Spalart-Allmaras turbulence solver')], name="",
    #                                     default='kEpsilon', update=nodeupdate)
    buoyancy: BoolProperty(name='', description='Thermal', default=0, update=nodeupdate)
    radiation: BoolProperty(name='', description='Radiation', default=0, update=nodeupdate)
    solar: BoolProperty(name='', description='Radiation', default=0, update=nodeupdate)
    sun: StringProperty(name="", description="Sun for solar radiation analysis", default="", update=nodeupdate)
    # buossinesq: BoolProperty(name='', description='Buossinesq approximation', default=0, update=nodeupdate)
    # stime: FloatProperty(name='', description='Simulation start time', min=0, max=10, default=0)
    dtime: FloatProperty(name='', description='False time step', min=0.001, max=10, default=0.5, precision=4, update=nodeupdate)
    etime: FloatProperty(name='', description='Simulation end time', min=1, default=5, update=nodeupdate)
    w_step: EnumProperty(items=[('0', 'None', 'No reference pressure'), ('1', 'Relative', 'Relative reference pressure'), ('2', 'Absolute', 'Absolute reference pressure')],
                        name='', description='Write step', default=0, update=nodeupdate)
    w_int: IntProperty(name="", description="Write interval in time steps", min=1, default=10, update=nodeupdate)
    pval: FloatProperty(name="", description="Field pressure (relative)", min=-500, max=500, default=0.0, update=nodeupdate)
    pnormval: FloatProperty(name="", description="Field pressure (normalised)", min=-500, max=500, default=0.0, update=nodeupdate)
    pabsval: IntProperty(name="", description="Field pressure (absolute)", min=0, max=10000000, default=100000, update=nodeupdate)
    p_ref: EnumProperty(items=[('0', 'None', 'No reference pressure'), ('1', 'Relative', 'Relative reference pressure'), ('2', 'Absolute', 'Absolute reference pressure')],
                        name='', description='Reference pressure', default=0, update=nodeupdate)
    p_ref_point: EnumProperty(items=ret_empty_menu, name='', description='Camera', update=nodeupdate)
    p_ref_val: FloatProperty(name="", description="Reference pressure value", min=-5000000, max=5000000, default=0.0, update=nodeupdate)
    uval: FloatVectorProperty(size=3, name='', attr='Velocity', default=[0, 0, 0], unit='VELOCITY', subtype='VELOCITY', min=-100, max=100, update=nodeupdate)
    uval_type: EnumProperty(name='', items=[('0', 'Vector', 'Air dirction and speed by vector'),
                                                 ('1', 'Azimuth', 'Transient simulation')], description='Velocity type', default='0', update=nodeupdate)
    uval_azi: FloatProperty(name="", description="Air direction azimuth (degrees from north)", min=0, max=360, default=0.0, update=nodeupdate)
    umag: FloatProperty(name="m/s", description="Air speed (m/s)", min=0, max=100, default=5.0, update=nodeupdate)
    tval: FloatProperty(name="K", description="Field Temperature (K)", min=0.0, max=500, default=293.14, update=nodeupdate)
    nutval: FloatProperty(name="", description="Nut domain value", min=0.0, max=500, default=0.0, update=nodeupdate)
    nutildaval: FloatProperty(name="", description="NuTilda domain value", min=0.0, max=500, default=0.0, update=nodeupdate)
    kval: FloatProperty(name="", description="k domain value", min=0.001, max=500, default=0.8, update=nodeupdate)
    epval: FloatProperty(name="", description="Epsilon domain value", min=0.001, max=500, default=0.1, update=nodeupdate)
    oval: FloatProperty(name="", description="Omega domain value", min=0.1, max=500, default=0.1, update=nodeupdate)
#    hval: FloatProperty(name="", description="Enthalpy domain value", min=0.1, max=500, default=0.1, update=nodeupdate)
   #  enval: FloatProperty(name="", description="Enthalpy domain value", min=0.1, max=500, default=0.1, update=nodeupdate)
    presid: FloatProperty(name="", description="p convergence criteria", precision=6, min=0.000001, max=0.01, default=0.0001, update=nodeupdate)
    uresid: FloatProperty(name="", description="U convergence criteria", precision=6, min=0.000001, max=0.5, default=0.0001, update=nodeupdate)
    keoresid: FloatProperty(name="", description="k/e/o convergence criteria", precision=6, min=0.000001, max=0.5, default=0.0001, update=nodeupdate)
    enresid: FloatProperty(name="", description="Enthalpy convergence criteria", precision=6, min=0.000001, max=0.5, default=0.0001, update=nodeupdate)
    aval: FloatProperty(name="", description="alphat value", min=0.0, max=500, default=0.1, update=nodeupdate)
    p_rghval: FloatProperty(name="", description="p_rgh value", min=-1000000, max=1000000, default=0.0, update=nodeupdate)
    Gval: FloatProperty(name="", description="Field radiation value", min=0.0, max=500, default=0.0, update=nodeupdate)
    radmodel: EnumProperty(name='', items=[('0', 'P1', 'P1 radiation model'), ('1', 'fvDOM', 'fvDOM radiation model')], description='Radiation model selection', default='0', update=nodeupdate)
    comfort: BoolProperty(name='', description='Comfort parameters', default=0, update=nodeupdate)
    clo: FloatProperty(name="clo", description="Clothing level", min=0.25, max=1.5, default=0.75, update=nodeupdate)
    met: FloatProperty(name="met", description="Metabolic rate", min=0.5, max=5, default=1, update=nodeupdate)
    rh: FloatProperty(name="%", description="Relative humidity", min=0, max=100, default=60, update=nodeupdate)
    age: BoolProperty(name='', description='Air age', default=0, update=nodeupdate)

    def init(self, context):
        self['exportstate'] = ''
        self.outputs.new('So_Flo_Case', 'Case out')
        nodecolour(self, 1)

        if context:
            context.scene.vi_params['flparams']['params'] = 'l'

    def draw_buttons(self, context, layout):
        newrow(layout, 'Scenrio:', self, 'scenario')
        newrow(layout, 'Parametric:', self, 'parametric')

        if self.parametric:
            row = layout.row()
            row.label(text = 'Frames:')
            col = row.column()
            subrow = col.row(align=True)
            subrow.prop(self, 'frame_start')
            subrow.prop(self, 'frame_end')

        newrow(layout, 'Time step:', self, 'dtime')
        newrow(layout, 'End time:', self, 'etime')
        newrow(layout, 'Write interval:', self, 'w_int')

        if self.scenario in ('2', '3', '4'):
            newrow(layout, 'Buoyancy:', self, 'buoyancy')

        if self.buoyancy:
            newrow(layout, 'Comfort:', self, 'comfort')

            if self.comfort:
                newrow(layout, 'Clothing', self, 'clo')
                newrow(layout, 'Metabolic', self, 'met')
                newrow(layout, 'Humidity', self, 'rh')

            newrow(layout, 'Age:', self, 'age')
            newrow(layout, 'Radiation:', self, 'radiation')
            newrow(layout, 'Field pressure:', self, 'pabsval')
            newrow(layout, 'Field p_rgh:', self, 'p_rghval')
        else:
            newrow(layout, 'Field pressure:', self, 'pnormval')

        newrow(layout, 'Reference pressure:', self, 'p_ref')

        if self.p_ref != '0':
            newrow(layout, 'Reference pressure:', self, 'p_ref_val')
            newrow(layout, 'Reference point:', self, 'p_ref_point')

        newrow(layout, 'Field k:', self, 'kval')
        newrow(layout, 'Field epsilon:', self, 'epval')

        if self.buoyancy:
            newrow(layout, 'Field T:', self, 'tval')
            newrow(layout, 'Field alphat:', self, 'aval')

            if self.radiation:
                newrow(layout, 'Solar:', self, 'solar')

                if self.solar:
                    layout.prop_search(self, 'sun', bpy.data, 'lights', text='Sun', icon='NONE')

                newrow(layout, 'Rad model:', self, 'radmodel')
                newrow(layout, 'G value:', self, 'Gval')

        newrow(layout, 'Field nut:', self, 'nutval')
        newrow(layout, 'Field velocity type:', self, 'uval_type')

        if self.uval_type == '0':
            newrow(layout, 'Field velocity:', self, 'uval')
        else:
            newrow(layout, 'Field flow azi:', self, 'uval_azi')
            newrow(layout, 'Field flow speed:', self, 'umag')

        newrow(layout, 'p Residual:', self, 'presid')
        newrow(layout, 'U Residual:', self, 'uresid')
        newrow(layout, 'k/epsilon Residual:', self, 'keoresid')

        if self.buoyancy:
            newrow(layout, 'e residual:', self, 'enresid')

        row = layout.row()
        row.operator("node.flovi_case", text = "Export")

    def update(self):
        for sock in self.outputs:
            socklink(sock, self.id_data.name)

    def pre_case(self, context):
        self.nodeupdate(context)

    def post_case(self):
        self['exportstate'] = self.ret_params()
        nodecolour(self, 0)

class No_Flo_NG(Node, ViNodes):
    '''Netgen mesh export'''
    bl_idname = 'No_Flo_NG'
    bl_label = 'FloVi NetGen'
    bl_icon = 'MESH_ICOSPHERE'

    def ret_params(self):
        return [str(x) for x in (self.poly, self.pcorr, self.acorr, self.maxcs, self.yang, self.grading, self.optimisations, self.fang)]

    def nodeupdate(self, context):
        nodecolour(self, self['exportstate'] != self.ret_params())

    poly: BoolProperty(name = '', description = 'Create polygonal mesh', default = 0, update = nodeupdate)
    pcorr: FloatProperty(name = "m", description = "Maximum distance for position correspondance", min = 0, max = 1, default = 0.1, update = nodeupdate)
    acorr: FloatProperty(name = "", description = "Minimum cosine for angular correspondance", min = 0, max = 1, default = 0.9, update = nodeupdate)
    maxcs: FloatProperty(name = "m", description = "Max global cell size", min = 0, max = 100, default = 1, update = nodeupdate)
    yang: FloatProperty(name = "deg", description = "Minimum angle for separate faces", min = 0, max = 90, default = 1, update = nodeupdate)
    grading: FloatProperty(name = "", description = "Small to large cell inflation", min = 0, max = 5, default = 0.3, update = nodeupdate)
    processors: IntProperty(name = "", description = "Number of processers", min = 0, max = 32, default = 1, update = nodeupdate)
    optimisations: IntProperty(name = "", description = "Number of optimisation steps", min = 0, max = 32, default = 3, update = nodeupdate)
    maxsteps: IntProperty(name = "", description = "Number of attempts", min = 0, max = 10, default = 3, update = nodeupdate)
    fang: FloatProperty(name = "deg", description = "Minimum angle for separate faces", min = 0, max = 90, default = 30, update = nodeupdate)
    geo_join: BoolProperty(name = '', description = 'Join Geometries', default = 0, update = nodeupdate)
    d_diff: BoolProperty(name = '', description = 'Extract geometries from domain', default = 0, update = nodeupdate)
    running: BoolProperty(name = '', description = '', default = 0)

    def init(self, context):
        self['exportstate'] = ''
        self.inputs.new('So_Flo_Case', 'Case in')
        self.outputs.new('So_Flo_Mesh', 'Mesh out')
        nodecolour(self, 1)

    def draw_buttons(self, context, layout):
        if sys.platform == 'linux':
            addonfolder = os.path.basename(os.path.dirname(os.path.abspath(__file__)))
            vi_prefs = bpy.context.preferences.addons['{}'.format(addonfolder)].preferences

            if os.path.isdir(vi_prefs.ofbin) and os.path.isfile(os.path.join(vi_prefs.ofbin, 'foamExec')):
                flo_libs[1] = 1

        if self.inputs and self.inputs['Case in'].links:
            if all(flo_libs):
                # newrow(layout, 'Join geometries:', self, 'geo_join')
                # if self.geo_join:
                #     newrow(layout, 'Domain extraction:', self, 'd_diff')
                newrow(layout, 'Cell size:', self, 'maxcs')
                newrow(layout, 'Position corr:', self, 'pcorr')
                newrow(layout, 'Angular corr:', self, 'acorr')
                newrow(layout, 'Distinction angle:', self, 'yang')
                newrow(layout, 'Inflation:', self, 'grading')
                newrow(layout, 'Optimisations:', self, 'optimisations')
                newrow(layout, 'Attempts:', self, 'maxsteps')
                newrow(layout, 'Polygonal:', self, 'poly')

                if not self.running:
                    row = layout.row()
                    row.operator("node.flovi_ng", text="Generate")

            elif not flo_libs[0]:
                row = layout.row()
                row.label(text = 'Netgen not found')

            elif not flo_libs[1]:
                if sys.platform == 'linux':
                    row = layout.row()
                    row.label(text='No OpenFOAM directory set')
                else:
                    row = layout.row()
                    row.label(text='No docker container running')

    def update(self):
        for sock in self.outputs:
            socklink(sock, self.id_data.name)

        self.running = 0

    def post_export(self):
        self['exportstate'] = self.ret_params()
        nodecolour(self, 0)

class No_Flo_Bound(Node, ViNodes):
    '''Openfoam boundary export'''
    bl_idname = 'No_Flo_Bound'
    bl_label = 'FloVi Boundary'
    bl_icon = 'MESH_ICOSPHERE'

    def nodeupdate(self, context):
        nodecolour(self, self['exportstate'] != [str(self.pv)])

    pv: BoolProperty(name='', description='Open Paraview', default=0, update=nodeupdate)

    def init(self, context):
        self['exportstate'] = ''
        self.inputs.new('So_Flo_Mesh', 'Mesh in')
        self.outputs.new('So_Flo_Con', 'Context out')
        nodecolour(self, 1)

    def draw_buttons(self, context, layout):
        if self.inputs and self.inputs['Mesh in'].links:
            newrow(layout, 'Paraview:', self, 'pv')
            row = layout.row()
            row.operator("node.flovi_bound", text="Generate")

    def update(self):
        for sock in self.outputs:
            socklink(sock, self.id_data.name)

    def post_export(self):
        self['exportstate'] = [str(self.pv)]
        nodecolour(self, 0)

class So_Flo_Con(NodeSocket):
    '''FloVi context socket'''
    bl_idname = 'So_Flo_Con'
    bl_label = 'FloVi context socket'

    valid = ['FloVi context']
    link_limit = 1

    def draw(self, context, layout, node, text):
        layout.label(text = text)

    def draw_color(self, context, node):
        return (0.0, 1.0, 1.0, 1.0)

class So_Flo_Case(NodeSocket):
    '''FloVi case socket'''
    bl_idname = 'So_Flo_Case'
    bl_label = 'FloVi Case socket'

    valid = ['FloVi case']
    link_limit = 1

    def draw(self, context, layout, node, text):
        layout.label(text = text)

    def draw_color(self, context, node):
        return (1.0, 1.0, 1.0, 1.0)


class So_Flo_Mesh(NodeSocket):
    '''FloVi mesh socket'''
    bl_idname = 'So_Flo_Mesh'
    bl_label = 'FloVi Mesh socket'

    valid = ['FloVi mesh']
    link_limit = 1

    def draw(self, context, layout, node, text):
        layout.label(text = text)

    def draw_color(self, context, node):
        return (0.5, 1.0, 0.0, 0.75)

    def draw_buttons(self, context, layout):
        row = layout.row()
        row.operator("node.blockmesh", text = "Export")

class No_Flo_BMesh(Node, ViNodes):
    '''Openfoam blockmesh export node'''
    bl_idname = 'No_Flo_BMesh'
    bl_label = 'FloVi BlockMesh'
    bl_icon = 'GRID'

    def nodeupdate(self, context):
        nodecolour(self, self['exportstate'] != [str(x) for x in (self.bm_xres, self.bm_yres, self.bm_zres, self.bm_xgrad, self.bm_ygrad, self.bm_zgrad)])

    bm_xres: IntProperty(name = "X", description = "Blockmesh X resolution", min = 0, max = 1000, default = 10, update = nodeupdate)
    bm_yres: IntProperty(name = "Y", description = "Blockmesh Y resolution", min = 0, max = 1000, default = 10, update = nodeupdate)
    bm_zres: IntProperty(name = "Z", description = "Blockmesh Z resolution", min = 0, max = 1000, default = 10, update = nodeupdate)
    bm_xgrad: FloatProperty(name = "X", description = "Blockmesh X simple grading", min = 0, max = 10, default = 1, update = nodeupdate)
    bm_ygrad: FloatProperty(name = "Y", description = "Blockmesh Y simple grading", min = 0, max = 10, default = 1, update = nodeupdate)
    bm_zgrad: FloatProperty(name = "Z", description = "Blockmesh Z simple grading", min = 0, max = 10, default = 1, update = nodeupdate)

    def init(self, context):
        self['exportstate'] = ''
        self.outputs.new('So_Flo_Mesh', 'Mesh out')
        nodecolour(self, 1)

    def draw_buttons(self, context, layout):
        split = layout.split()
        col = split.column(align=True)
        col.label(text="Cell resolution:")
        col.prop(self, "bm_xres")
        col.prop(self, "bm_yres")
        col.prop(self, "bm_zres")
        col = split.column(align=True)
        col.label(text="Cell grading:")
        col.prop(self, "bm_xgrad")
        col.prop(self, "bm_ygrad")
        col.prop(self, "bm_zgrad")
        row = layout.row()
        row.operator("node.flovi_bm", text = "Export")

    def update(self):
        for sock in self.outputs:
            socklink(sock, self.id_data.name)

    def export(self):
        self.exportstate = [str(x) for x in (self.bm_xres, self.bm_yres, self.bm_zres, self.bm_xgrad, self.bm_ygrad, self.bm_zgrad)]
        nodecolour(self, 0)

class No_Flo_Sim(Node, ViNodes):
    '''Openfoam simulation node'''
    bl_idname = 'No_Flo_Sim'
    bl_label = 'FloVi Simulation'
    bl_icon = 'FORCE_HARMONIC'

    def nodeupdate(self, context):
        nodecolour(self, self['exportstate'] != [str(x) for x in (self.processes, self.pv)])

    processes: IntProperty(name="", description="Number of processors", min=1, max=1000, default=1, update=nodeupdate)
    pv: BoolProperty(name="", description="Open paraview", default=0, update=nodeupdate)

    def init(self, context):
        self['exportstate'] = ''
        self.inputs.new('So_Flo_Con', 'Context in')
        self.outputs.new('So_Vi_Res', 'Results out')
        nodecolour(self, 1)

    def draw_buttons(self, context, layout):
        newrow(layout, 'Processes:', self, 'processes')
        newrow(layout, 'Paraview:', self, 'pv')
        row = layout.row()
        row.operator("node.flovi_sim", text = "Calculate")

    def update(self):
        for sock in self.outputs:
            socklink(sock, self.id_data.name)

    def presim(self):
        expnode = self.inputs['Context in'].links[0].from_node
        return (expnode.convergence, expnode.econvergence, expnode['residuals'], expnode.processes, expnode.solver)

    def post_sim(self):
        self['exportstate'] = [str(x) for x in (self.processes, self.pv)]
        nodecolour(self, 0)

####################### Vi Nodes Categories ##############################

vi_process = [NodeItem("No_Li_Geo", label="LiVi Geometry"), NodeItem("No_Li_Con", label="LiVi Context"),
              NodeItem("No_En_Geo", label="EnVi Geometry"), NodeItem("No_En_Con", label="EnVi Context"), NodeItem("No_Flo_Case", label="FloVi Case"),
              NodeItem("No_Flo_NG", label="FloVi NetGen"), NodeItem("No_Flo_Bound", label="FloVi Boundary")]

vi_edit = [NodeItem("No_Text", label="Text Edit")]

vi_analysis = [NodeItem("No_Vi_SP", label="Sun Path"), NodeItem("No_Vi_WR", label="Wind Rose"),
             NodeItem("No_Vi_SVF", label="Sky View"), NodeItem("No_Vi_SS", label="Shadow map"),
             NodeItem("No_Li_Sim", label="LiVi Simulation"), NodeItem("No_En_Sim", label="EnVi Simulation"),
             NodeItem("No_Flo_Sim", label="FloVi Simulation"), NodeItem("No_Vi_EC", label="Embodied Carbon")]

vi_anim = [NodeItem("No_Anim", label="Parametric")]
vi_display = [NodeItem("No_Vi_Chart", label="Chart"), NodeItem("No_Vi_HMChart", label="Heatmap"), NodeItem("No_Vi_Metrics", label="Metrics")]
vi_out = [NodeItem("No_CSV", label="CSV")]
vi_image = [NodeItem("No_Li_Im", label="LiVi Image"),
            NodeItem("No_Li_Gl", label="LiVi Glare"), NodeItem("No_Li_Fc", label="LiVi False-colour")]
vi_input = [NodeItem("No_Loc", label="VI Location"), NodeItem("No_ASC_Import", label="ASC Import")]

# Names must be unique
vinode_categories = [ViNodeCategory("Output", "Output Nodes", items=vi_out),
                     ViNodeCategory("ViPara", "Parametric Nodes", items=vi_anim),
                     ViNodeCategory("Edit", "Edit Nodes", items=vi_edit),
                     ViNodeCategory("Image", "Image Nodes", items=vi_image),
                     ViNodeCategory("Display", "Display Nodes", items=vi_display),
                     ViNodeCategory("Analysis", "Analysis Nodes", items=vi_analysis),
                     ViNodeCategory("Process", "Process Nodes", items=vi_process),
                     ViNodeCategory("Input", "Input Nodes", items=vi_input)]


# EnVi Network Nodes

class EnViNetwork(NodeTree):
    '''A node tree for the creation of EnVi advanced networks.'''
    bl_idname = 'EnViN'
    bl_label = 'EnVi Network'
    bl_icon = 'FORCE_WIND'


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
    sn: IntProperty()
    uvalue: StringProperty()
    viuid: StringProperty()
    link_limit = 1

    def draw(self, context, layout, node, text):
        layout.label(text = text)

    def draw_color(self, context, node):
        return (0.5, 0.2, 0.0, 0.75)

    def ret_valid(self, node):
        return ['Boundary']

class So_En_Net_Sched(NodeSocket):
    '''Fraction schedule socket'''
    bl_idname = 'So_En_Net_Sched'
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

    sn: IntProperty()
    valid = ['Sub-surface']
    viuid: StringProperty()
    link_limit = 1

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

    sn: IntProperty()
    valid = ['Surface']
    viuid: StringProperty()
    link_limit = 1

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

class No_En_Net_Zone(Node, EnViNodes):
    '''Node describing a simulation zone'''
    bl_idname = 'No_En_Net_Zone'
    bl_label = 'EnVi Zone'
    bl_icon = 'CUBE'

    def zupdate(self, context):
        self.afs = 0
        col = bpy.data.collections[self.zone]

        for obj in [o for o in col.objects if not o.vi_params.embodied]:
            odm = [m.material for m in obj.material_slots]
            bfacelist = sorted([face for face in obj.data.polygons if get_con_node(odm[face.material_index].vi_params).envi_con_con == 'Zone'], key=lambda face: -face.center[2])
            sfacelist = sorted([face for face in obj.data.polygons if get_con_node(odm[face.material_index].vi_params).envi_afsurface == 1 and get_con_node(odm[face.material_index].vi_params).envi_con_type not in ('Window', 'Door')], key=lambda face: -face.center[2])
            ssfacelist = sorted([face for face in obj.data.polygons if get_con_node(odm[face.material_index].vi_params).envi_afsurface == 1 and get_con_node(odm[face.material_index].vi_params).envi_con_type in ('Window', 'Door')], key=lambda face: -face.center[2])
            [self.outputs.remove(oname) for oname in self.outputs if oname.bl_idname in ('So_En_Net_Bound', 'So_En_Net_SFlow', 'So_En_Net_SSFlow')]# and oname not in bsocklist + ssocklist + sssocklist]
            [self.inputs.remove(iname) for iname in self.inputs if iname.bl_idname in ('So_En_Net_Bound', 'So_En_Net_SFlow', 'So_En_Net_SSFlow')]# and iname not in bsocklist + ssocklist + sssocklist]

            for bface in bfacelist:
                self.outputs.new('So_En_Net_Bound', '{}_{}_b'.format(odm[bface.material_index].name, bface.index)).sn = bface.index
                self.outputs[-1].viuid = '{}#{}'.format(obj.name, obj.data.polygon_layers_int["viuid"].data[bface.index].value)
                self.outputs[-1].link_limit = 1
                self.inputs.new('So_En_Net_Bound', '{}_{}_b'.format(odm[bface.material_index].name, bface.index)).sn = bface.index
                self.inputs[-1].viuid = '{}#{}'.format(obj.name, obj.data.polygon_layers_int["viuid"].data[bface.index].value)
                self.inputs[-1].link_limit = 1
            for sface in sfacelist:
                self.afs += 1
                self.outputs.new('So_En_Net_SFlow', '{}_{}_s'.format(odm[sface.material_index].name, sface.index)).sn = sface.index
                self.outputs[-1].viuid = '{}#{}'.format(obj.name, obj.data.polygon_layers_int["viuid"].data[sface.index].value)
                self.outputs[-1].link_limit = 1
                self.inputs.new('So_En_Net_SFlow', '{}_{}_s'.format(odm[sface.material_index].name, sface.index)).sn = sface.index
                self.inputs[-1].viuid = '{}#{}'.format(obj.name, obj.data.polygon_layers_int["viuid"].data[sface.index].value)
                self.inputs[-1].link_limit = 1
            for ssface in ssfacelist:
                self.afs += 1
                self.outputs.new('So_En_Net_SSFlow', '{}_{}_s'.format(odm[ssface.material_index].name, ssface.index)).sn = ssface.index
                self.outputs[-1].viuid = '{}#{}'.format(obj.name, obj.data.polygon_layers_int["viuid"].data[ssface.index].value)
                self.outputs[-1].link_limit = 1
                self.inputs.new('So_En_Net_SSFlow', '{}_{}_s'.format(odm[ssface.material_index].name, ssface.index)).sn = ssface.index
                self.inputs[-1].viuid = '{}#{}'.format(obj.name, obj.data.polygon_layers_int["viuid"].data[ssface.index].value)
                self.inputs[-1].link_limit = 1

        self.vol_update(context)
        self['nbound'] = len(bfacelist)
        self['nsflow'] = len(sfacelist)
        self['nssflow'] = len(ssfacelist)

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

    def init(self, context):
        self['nbound'] = 0
        self['nsflow'] = 0
        self['nssflow'] = 0
        self['tsps'] = 1
        self['volume'] = 10
        self.inputs.new('So_En_Net_Hvac', 'HVAC')
        self.inputs.new('So_En_Net_Occ', 'Occupancy')
        self.inputs.new('So_En_Net_Eq', 'Equipment')
        self.inputs.new('So_En_Net_Inf', 'Infiltration')
        self.inputs.new('So_En_Net_Sched', 'TSPSchedule')
        self.inputs.new('So_En_Net_Sched', 'VASchedule')

    def update(self):
        sflowdict = {'So_En_Net_SFlow': 'Envi surface flow', 'So_En_Net_SSFlow': 'Envi sub-surface flow'}
        [bi, si, ssi, bo, so , sso] = [1, 1, 1, 1, 1, 1]

        if len(self.inputs) > 5 and len(self.outputs) == self['nbound'] + self['nsflow'] + self['nssflow']:
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
                            remlink(self.id_data, [inp.links[0]])

            for outp in [outp for outp in self.outputs if outp.bl_idname in ('So_En_Net_Bound', 'So_En_Net_SFlow', 'So_En_Net_SSFlow')]:
                if outp.bl_idname == 'So_En_Bound' and not outp.hide and not outp.links:
                    bo = 0
                elif outp.bl_idname in sflowdict:
                    if (not outp.hide and not outp.links) or (outp.links and outp.links[0].to_node.bl_label != sflowdict[outp.bl_idname]):
                        so = 0
                        if outp.links:
                            remlink(self.id_data, [outp.links[0]])

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




class No_En_Net_Hvac(Node, EnViNodes):
    '''Zone HVAC node'''
    bl_idname = 'No_En_Net_Hvac'
    bl_label = 'HVAC'

    def hupdate(self, context):
        self.h = 1 if self.envi_hvachlt != '4' else 0
        self.c = 1 if self.envi_hvacclt != '4' else 0
        self['hc'] = ('', 'SingleHeating', 'SingleCooling', 'DualSetpoint')[(not self.h and not self.c, self.h and not self.c, not self.h and self.c, self.h and self.c).index(1)]

    def hupdate_nc(self):
        self.h = 1 if self.envi_hvachlt != '4' else 0
        self.c = 1 if self.envi_hvacclt != '4' else 0
        self['hc'] = ('', 'SingleHeating', 'SingleCooling', 'DualSetpoint')[(not self.h and not self.c, self.h and not self.c, not self.h and self.c, self.h and self.c).index(1)]

    envi_hvact: bprop("", "", False)
    envi_hvacht: fprop(u'\u00b0C', "Heating temperature:", 1, 99, 50)
    envi_hvacct: fprop(u'\u00b0C', "Cooling temperature:", -10, 20, 13)
    envi_hvachlt: EnumProperty(items = [('0', 'LimitFlowRate', 'LimitFlowRate'), ('1', 'LimitCapacity', 'LimitCapacity'), ('2', 'LimitFlowRateAndCapacity', 'LimitFlowRateAndCapacity'), ('3', 'NoLimit', 'NoLimit'), ('4', 'None', 'No heating')], name = '', description = "Heating limit type", default = '4', update = hupdate)
    envi_hvachaf: FloatProperty(name = u'm\u00b3/s', description = "Heating air flow rate", min = 0, max = 60, default = 1, precision = 4, options={'ANIMATABLE'})
    envi_hvacshc: fprop("W", "Sensible heating capacity", 0, 10000, 1000)
    envi_hvacclt: EnumProperty(items = [('0', 'LimitFlowRate', 'LimitFlowRate'), ('1', 'LimitCapacity', 'LimitCapacity'), ('2', 'LimitFlowRateAndCapacity', 'LimitFlowRateAndCapacity'), ('3', 'NoLimit', 'NoLimit'), ('4', 'None', 'No cooling')], name = '', description = "Cooling limit type", default = '4', update = hupdate)
    envi_hvaccaf: FloatProperty(name = u'm\u00b3/s', description = "Cooling air flow rate", min = 0, max = 60, default = 1, precision = 4)
    envi_hvacscc: fprop("W", "Sensible cooling capacity", 0, 10000, 1000)
    envi_hvacoam: eprop([('0', 'None', 'None'), ('1', 'Flow/Zone', 'Flow/Zone'), ('2', 'Flow/Person', 'Flow/Person'), ('3', 'Flow/Area', 'Flow/Area'), ('4', 'Sum', 'Sum'), ('5', 'Maximum ', 'Maximum'), ('6', 'ACH/Detailed', 'ACH/Detailed')], '', "Cooling limit type", '0')
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
    envi_ctsp: FloatProperty(name = u'\u00b0'+"C", description = "Temperature", min = 0, max = 50, default = 24)

    def init(self, context):
        self['hc'] = ''
        self['ctdict'] = {'DualSetpoint': 4, 'SingleHeating': 1, 'SingleCooling': 2}
        self['limittype'] = {'0': 'LimitFlowRate', '1': 'LimitCapacity', '2': 'LimitFlowRateAndCapacity', '3': 'NoLimit', '4': ''}
        self.outputs.new('So_En_Net_Hvac', 'HVAC')
        self.inputs.new('So_En_Net_Sched', 'Schedule')
        self.inputs.new('So_En_Net_TSched', 'HSchedule')
        self.inputs.new('So_En_Net_TSched', 'CSchedule')

    def draw_buttons(self, context, layout):
        # row = layout.row()
        # row.label(text = 'HVAC Template:')
        # row.prop(self, 'envi_hvact')
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
                   '', '', 'ConstantSupplyHumidityRatio', '', 'ConstantSupplyHumidityRatio', (zn + ' Outdoor Air', '')[self.envi_hvacoam == '0'], '', '', '', ('None', 'Sensible')[int(self.envi_hvachr)], '{:.2f}'.format(self.envi_hvachre), '')
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
        self.hupdate_nc()
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
        bpy.context.scene.vi_params['enparams']['hvactemplate'] = 1
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
    bl_icon = 'ARMATURE_DATA'

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
        self.inputs.new('So_En_Net_Sched', 'OSchedule')
        self.inputs.new('So_En_Net_Sched', 'ASchedule')
        self.inputs.new('So_En_Net_Sched', 'WSchedule')
        self.inputs.new('So_En_Net_Sched', 'VSchedule')
        self.inputs.new('So_En_Net_Sched', 'CSchedule')

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
        self.inputs.new('So_En_Net_Sched', 'Schedule')

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
        self.inputs.new('So_En_Net_Sched', 'Schedule')

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

        if self.envi_inftype != '0':
            inflist[int(infdict[self.envi_inftype])] = '{:.4f}'.format(self.envi_inflevel)
            params = ('Name', 'Zone or ZoneList Name', 'Schedule Name', 'Design Flow Rate Calculation Method', 'Design Flow Rate {m3/s}', 'Flow per Zone Floor Area {m3/s-m2}',
                'Flow per Exterior Surface Area {m3/s-m2}', 'Air Changes per Hour {1/hr}', 'Constant Term Coefficient', 'Temperature Term Coefficient',
                    'Velocity Term Coefficient', 'Velocity Squared Term Coefficient')
            paramvs = [zn + '_infiltration', zn, zn + '_infsched', infildict[self.envi_inftype]] + inflist + [1, 0, 0, 0]
            return epentry('ZoneInfiltration:DesignFlowRate', params, paramvs)
        else:
            return ''

class No_En_Net_SSFlow(Node, EnViNodes):
    '''Sub-surface airflow node'''
    bl_idname = 'No_En_Net_SSFlow'
    bl_label = 'Envi sub-surface flow'
    bl_icon = 'FORCE_WIND'

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
    uvof: FloatProperty(default = 100, min = 1, max = 100, name = "", description = 'Indoor and Outdoor Temperature Difference Upper Limit For Minimum Venting Open Factor (deltaC)')
    amfcc: FloatProperty(default = 0.001, min = 0.00001, max = 1, precision = 5, name = "", description = 'Air Mass Flow Coefficient When Opening is Closed (kg/s-m)')
    amfec: FloatProperty(default = 0.65, min = 0.5, max = 1, name = '', description =  'Air Mass Flow Exponent When Opening is Closed (dimensionless)')
    lvo: EnumProperty(items = [('NonPivoted', 'NonPivoted', 'Non pivoting opening'), ('HorizontallyPivoted', 'HPivoted', 'Horizontally pivoting opening')], name = '', default = 'NonPivoted', description = 'Type of Rectanguler Large Vertical Opening (LVO)')
    ecl: FloatProperty(default = 0.0, min = 0, name = '', description = 'Extra Crack Length or Height of Pivoting Axis (m)')
    noof: IntProperty(default = 2, min = 2, max = 4, name = '', description = 'Number of Sets of Opening Factor Data')
    spa: IntProperty(default = 90, min = 0, max = 90, name = '', description = 'Sloping Plane Angle')
    dcof: FloatProperty(default = 0.7, min = 0.01, max = 1, name = '', description = 'Discharge Coefficient')
    ddtw: FloatProperty(default = 0.001, min = 0.001, max = 10, name = '', description = 'Minimum Density Difference for Two-way Flow')
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
    of1: FloatProperty(default = 0.0, min = 0.0, max = 0, name = '', description = 'Opening Factor 1 (dimensionless)')
    of2: FloatProperty(default = 0.0, min = 0.0, max = 0, name = '', description = 'Opening Factor 2 (dimensionless)')
    of3: FloatProperty(default = 0.0, min = 0.0, max = 0, name = '', description = 'Opening Factor 3 (dimensionless)')
    of4: FloatProperty(default = 1.0, min = 0.01, max = 1, name = '', description = 'Opening Factor 4 (dimensionless)')
    dcof1: FloatProperty(default = 0.65, min = 0.01, max = 1, name = '', description = 'Discharge Coefficient for Opening Factor 1 (dimensionless)')
    dcof2: FloatProperty(default = 0.65, min = 0.01, max = 1, name = '', description = 'Discharge Coefficient for Opening Factor 2 (dimensionless)')
    dcof3: FloatProperty(default = 0.65, min = 0.01, max = 1, name = '', description = 'Discharge Coefficient for Opening Factor 3 (dimensionless)')
    dcof4: FloatProperty(default = 0.65, min = 0.01, max = 1, name = '', description = 'Discharge Coefficient for Opening Factor 4 (dimensionless)')
    wfof1: FloatProperty(default = 0.0, min = 0, max = 1, name = '', description = 'Width Factor for Opening Factor 1 (dimensionless)')
    wfof2: FloatProperty(default = 0.0, min = 0, max = 1, name = '', description = 'Width Factor for Opening Factor 2 (dimensionless)')
    wfof3: FloatProperty(default = 0.0, min = 0, max = 1, name = '', description = 'Width Factor for Opening Factor 3 (dimensionless)')
    wfof4: FloatProperty(default = 0.0, min = 0, max = 1, name = '', description = 'Width Factor for Opening Factor 4 (dimensionless)')
    hfof1: FloatProperty(default = 0.0, min = 0, max = 1, name = '', description = 'Height Factor for Opening Factor 1 (dimensionless)')
    hfof2: FloatProperty(default = 0.0, min = 0, max = 1, name = '', description = 'Height Factor for Opening Factor 2 (dimensionless)')
    hfof3: FloatProperty(default = 0.0, min = 0, max = 1, name = '', description = 'Height Factor for Opening Factor 3 (dimensionless)')
    hfof4: FloatProperty(default = 0.0, min = 0, max = 1, name = '', description = 'Height Factor for Opening Factor 4 (dimensionless)')
    sfof1: FloatProperty(default = 0.0, min = 0, max = 1, name = '', description = 'Start Height Factor for Opening Factor 1 (dimensionless)')
    sfof2: FloatProperty(default = 0.0, min = 0, max = 1, name = '', description = 'Start Height Factor for Opening Factor 2 (dimensionless)')
    sfof3: FloatProperty(default = 0.0, min = 0, max = 1, name = '', description = 'Start Height Factor for Opening Factor 3 (dimensionless)')
    sfof4: FloatProperty(default = 0.0, min = 0, max = 1, name = '', description = 'Start Height Factor for Opening Factor 4 (dimensionless)')
    dcof: FloatProperty(default = 0.65, min = 0.01, max = 1, name = '', description = 'Discharge Coefficient')
    extnode: BoolProperty(default = 0)
    actlist = [("0", "Opening factor", "Actuate the opening factor")]
    acttype: EnumProperty(name="", description="Actuator type", items=actlist, default='0')
    compdict = {'0': 'AirFlow Network Window/Door Opening'}
    actdict =  {'0': ('Venting Opening Factor', 'of')}
    adict = {'Window': 'win', 'Door': 'door'}

    def init(self, context):
        self['init'] = 1
        self['ela'] = 1.0
        self.inputs.new('So_En_Net_Sched', 'VASchedule')
        self.inputs.new('So_En_Net_Sched', 'TSPSchedule')
        self.inputs['TSPSchedule'].hide = True
        self.inputs.new('So_En_Net_SSFlow', 'Node 1', identifier='Node1_s').link_limit = 1
        self.inputs.new('So_En_Net_SSFlow', 'Node 2', identifier='Node2_s').link_limit = 1
        self.outputs.new('So_En_Net_SSFlow', 'Node 1', identifier='Node1_s').link_limit = 1
        self.outputs.new('So_En_Net_SSFlow', 'Node 2', identifier='Node2_s').link_limit = 1
        self.outputs.new('So_Anim', 'Parameter')
        self.color = (1.0, 0.3, 0.3)
        self['layoutdict'] = {'SO':(('Closed FC', 'amfcc'), ('Closed FE', 'amfec'), ('Density diff', 'ddtw'), ('DC', 'dcof')), 'DO':(('Closed FC', 'amfcc'), ('Closed FE', 'amfec'),
                           ('Opening type', 'lvo'), ('Crack length', 'ecl'), ('OF Number', 'noof'), ('OF1', 'of1'), ('DC1', 'dcof1'), ('Width OF1', 'wfof1'), ('Height OF1', 'hfof1'),
                            ('Start height OF1', 'sfof1'), ('OF2', 'of2'), ('DC2', 'dcof2'), ('Width OF2', 'wfof2'), ('Height OF2', 'hfof2'), ('Start height OF2', 'sfof2')),
                            'OF3': (('OF3', 'of3'), ('DC3', 'dcof3'), ('Width OF3', 'wfof3'), ('Height OF3', 'hfof3'), ('Start height OF3', 'sfof3')),
                            'OF4': (('OF4', 'of4'), ('DC4', 'dcof4'), ('Width OF4', 'wfof4'), ('Height OF4', 'hfof4'), ('Start height OF4', 'sfof4')),
                            'HO': (('Closed FC', 'amfcc'), ('Closed FE', 'amfec'), ('Slope', 'spa'), ('DC', 'dcof')), 'Crack': (('Coefficient', 'amfc'), ('Exponent', 'amfe'), ('Factor', 'of1')),
                            'ELA': (('ELA', '["ela"]'), ('DC', 'dcof'), ('PA diff', 'rpd'), ('FE', 'fe'))}

    def update(self):
        self.inputs[2].viuid = '{}#{}'.format(self.name, 1)
        self.inputs[3].viuid = '{}#{}'.format(self.name, 2)
        self.outputs[0].viuid = '{}#{}'.format(self.name, 3)
        self.outputs[1].viuid = '{}#{}'.format(self.name, 4)

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

            if self.outputs.get('Parameter'):
                sockhide(self, ('Node 1', 'Node 2'))

            self.legal()

    def draw_buttons(self, context, layout):
        layout.prop(self, 'linkmenu')
        if self.linkmenu in ('SO', 'DO', 'HO'):
            newrow(layout, 'Win/Door OF:', self, 'wdof1')
            newrow(layout, "Control type:", self, 'controls')
            if self.linkmenu in ('SO', 'DO') and self.controls == 'Temperature':
                newrow(layout, "Limit OF:", self, 'mvof')
                newrow(layout, "Lower deltaT:", self, 'lvof')
                newrow(layout, "Upper deltaT:", self, 'uvof')

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
                        snames.append('{}{}_{}'.format(('win-', 'door-')[get_con_node(bpy.data.materials[othersock.name[:-len(str(sn))-3]].vi_params).envi_con_type == 'Door'], zn, sn))
                        params = ('Surface Name', 'Leakage Component Name', 'External Node Name', 'Window/Door Opening Factor')
                        paramvs = (snames[-1], '{}_{}'.format(self.name, self.linkmenu), en, '{:.5f}'.format(self.wdof1))

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
    bl_label = 'Envi External'
    bl_icon = 'FORCE_WIND'

    height: FloatProperty(default = 1.0)
    wpc1: FloatProperty(name = '', default = 0, min = -1, max = 1)
    wpc2: FloatProperty(name = '', default = 0, min = -1, max = 1)
    wpc3: FloatProperty(name = '', default = 0, min = -1, max = 1)
    wpc4: FloatProperty(name = '', default = 0, min = -1, max = 1)
    wpc5: FloatProperty(name = '', default = 0, min = -1, max = 1)
    wpc6: FloatProperty(name = '', default = 0, min = -1, max = 1)
    wpc7: FloatProperty(name = '', default = 0, min = -1, max = 1)
    wpc8: FloatProperty(name = '', default = 0, min = -1, max = 1)
    wpc9: FloatProperty(name = '', default = 0, min = -1, max = 1)
    wpc10: FloatProperty(name = '', default = 0, min = -1, max = 1)
    wpc11: FloatProperty(name = '', default = 0, min = -1, max = 1)
    wpc12: FloatProperty(name = '', default = 0, min = -1, max = 1)
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
        direcs = ('N', 'NNE', 'NE', 'ENE', 'E', 'ESE', 'SE', 'SSE', 'S', 'SSW', 'SW', 'WSW', 'W', 'WNW', 'NW', 'NNW')

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
    bl_icon = 'FORCE_WIND'

    linktype = [("Crack", "Crack", "Crack aperture used for leakage calculation"),
        ("ELA", "ELA", "Effective leakage area")]

    linkmenu: EnumProperty(name="Type", description="Linkage type", items=linktype, default='ELA')
    ela: FloatProperty(default = 0.1, min = 0.00001, max = 1, name = "", description = 'Effective leakage area', options = {'SKIP_SAVE'})
    of: FloatProperty(default = 0.1, min = 0.001, max = 1, name = "", description = 'Opening Factor 1 (dimensionless)')
    ecl: FloatProperty(default = 0.0, min = 0, name = '', description = 'Extra Crack Length or Height of Pivoting Axis (m)')
    dcof: FloatProperty(default = 0.7, min = 0, max = 1, name = '', description = 'Discharge Coefficient')
    amfc: FloatProperty(min = 0.001, max = 1, default = 0.01, name = "", precision = 5, description = 'Flow coefficient', options = {'SKIP_SAVE'})
    amfe: FloatProperty(min = 0.5, max = 1, default = 0.65, name = "", precision = 3, description = 'Flow exponent', options = {'SKIP_SAVE'})
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
        self.outputs.new('So_Anim', 'Parameter')

    def update(self):
        self.inputs[0].viuid = '{}#{}'.format(self.name, 1)
        self.inputs[1].viuid = '{}#{}'.format(self.name, 2)
        self.outputs[0].viuid = '{}#{}'.format(self.name, 3)
        self.outputs[1].viuid = '{}#{}'.format(self.name, 4)

        for sock in self.outputs:
            socklink(sock, self.id_data.name)

        if self.linkmenu == 'ELA':
            retelaarea(self)

        self.extnode = 0

        for sock in self.inputs[:] + self.outputs[:]:
            for l in sock.links:
                if (l.from_node, l.to_node)[sock.is_output].bl_idname == 'No_En_Net_Ext':
                    self.extnode = 1

        if self.outputs.get('Parameter'):
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
        paradict = {}

        for p in self.bl_rna.properties:
            if p.is_skip_save and p.identifier in [l.to_node.parameter for l in self.outputs['Parameter'].links]:
                for l in self.outputs['Parameter'].links:
                    if p.identifier == l.to_node.parameter and l.to_node.anim_file:
                        tf = bpy.data.texts[l.to_node.anim_file]
                        setattr(self, p.identifier, ret_param(getattr(self, p.identifier), tf.as_string().split('\n')[bpy.context.scene.frame_current - bpy.context.scene.vi_params['enparams']['fs']]))

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
                        snames.append('{}_{}'.format(zn, sn))
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
        except Exception:
            nodecolour(self, 1)

class No_En_Net_ACon(Node, EnViNodes):
    '''Node defining the overall airflow network simulation'''
    bl_idname = 'No_En_Net_ACon'
    bl_label = 'EnVi AFN'
    bl_icon = 'FORCE_WIND'

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
        except Exception:
            pass

class No_En_Net_WPC(Node, EnViNodes):
    '''Node describing Wind Pressure Coefficient array'''
    bl_idname = 'No_En_Net_WPC'
    bl_label = 'EnVi WPC'
    bl_icon = 'FORCE_WIND'

    (ang1, ang2, ang3, ang4, ang5, ang6, ang7, ang8, ang9, ang10, ang11, ang12) = [IntProperty(name = '', default = 0, min = 0, max = 360) for x in range(12)]

    def init(self, context):
        self.outputs.new('So_En_Net_WPC', 'WPC values')

    def draw_buttons(self, context, layout):
        row = layout.row()
        row.label(text = 'WPC Angles')

        for w in range(1, 13):
            row = layout.row()
            row.prop(self, 'ang{}'.format(w))

    def update(self):
        for sock in self.outputs:
            socklink(sock, self.id_data.name)

        self.id_data.interface_update(bpy.context)

    def epwrite(self):
        angs = (self.ang1,self.ang2, self.ang3, self.ang4, self.ang5, self.ang6, self.ang7, self.ang8, self.ang9, self.ang10, self.ang11, self.ang12)
        aparamvs = ['WPC Array'] + [wd for w, wd in enumerate(angs) if wd not in angs[:w]]
        aparams = ['Name'] + ['Wind Direction {} (deg)'.format(w + 1) for w in range(len(aparamvs) - 1)]
        return (epentry('AirflowNetwork:MultiZone:WindPressureCoefficientArray', aparams, aparamvs), len(aparamvs) - 1)

class No_En_Net_Sched(Node, EnViNodes):
    '''Node describing a schedule'''
    bl_idname = 'No_En_Net_Sched'
    bl_label = 'Net Schedule'
    bl_icon = 'TIME'

    def tupdate(self, context):
        try:
            if self.source == '1':
                if os.path.isfile(bpy.path.abspath(self.select_file)):
                    nodecolour(self, 0)
                else:
                    nodecolour(self, 1)
            else:
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

        except Exception:
            nodecolour(self, 1)

    source: EnumProperty(name = '', items = [("0", "Node", "Generate schedule within the node"), ("1", "File", "Select schedule file")], default = '0', update = tupdate)
    select_file: StringProperty(name="", description="Name of the variable file", default="", subtype="FILE_PATH", update = tupdate)
    cn: IntProperty(name = "", default = 1, min = 1)
    rtsat: IntProperty(name = "", default = 0, min = 0)
    hours: IntProperty(name = "", default = 8760, min = 1, max = 8760)
    delim: EnumProperty(name = '', items = [("Comma", "Comma", "Comma delimiter"), ("Space", "Space", "space delimiter")], default = 'Comma')
#    generate_file: StringProperty(default = "", name = "")
    u1: StringProperty(name = "", description = "Valid entries (; separated for each 'For', comma separated for each day, space separated for each time value pair)", update = tupdate)
    u2: StringProperty(name = "", description = "Valid entries (; separated for each 'For', comma separated for each day, space separated for each time value pair)", update = tupdate)
    u3: StringProperty(name = "", description = "Valid entries (; separated for each 'For', comma separated for each day, space separated for each time value pair)", update = tupdate)
    u4: StringProperty(name = "", description = "Valid entries (; separated for each 'For', comma separated for each day, space separated for each time value pair)", update = tupdate)
    f1: StringProperty(name = "", description = "Valid entries (space separated): AllDays, Weekdays, Weekends, Monday, Tuesday, Wednesday, Thursday, Friday, Saturday, Sunday, AllOtherDays", update = tupdate)
    f2: StringProperty(name = "", description = "Valid entries (space separated): AllDays, Weekdays, Weekends, Monday, Tuesday, Wednesday, Thursday, Friday, Saturday, Sunday, AllOtherDays", update = tupdate)
    f3: StringProperty(name = "", description = "Valid entries (space separated): AllDays, Weekdays, Weekends, Monday, Tuesday, Wednesday, Thursday, Friday, Saturday, Sunday, AllOtherDays", update = tupdate)
    f4: StringProperty(name = "", description = "Valid entries (space separated): AllDays, Weekdays, Weekends, Monday, Tuesday, Wednesday, Thursday, Friday, Saturday, Sunday, AllOtherDays", update = tupdate)
    t1: IntProperty(name = "", default = 365, min = 1, max = 365, update = tupdate)
    t2: IntProperty(name = "", default = 365, min = 1, max = 365, update = tupdate)
    t3: IntProperty(name = "", default = 365, min = 1, max = 365, update = tupdate)
    t4: IntProperty(name = "", default = 365, min = 1, max = 365, update = tupdate)

    def init(self, context):
        self.outputs.new('So_En_Net_Sched', 'Schedule')
        self['scheddict'] = {'TSPSchedule': 'Any Number', 'VASchedule': 'Fraction', 'Fan Schedule': 'Fraction', 'HSchedule': 'Temperature', 'CSchedule': 'Temperature'}
        self.tupdate(context)
        nodecolour(self, 1)

    def draw_buttons(self, context, layout):
        uvals, u = (1, self.u1, self.u2, self.u3, self.u4), 0
        tvals = (0, self.t1, self.t2, self.t3, self.t4)
        newrow(layout, 'Source', self, 'source')

        if self.source == "1":
            newrow(layout, 'Select', self, 'select_file')
            newrow(layout, 'Columns', self, 'cn')
            newrow(layout, 'Skip rows', self, 'rtsat')
            newrow(layout, 'Delimiter', self, 'delim')

        if self.source != "1":
            while uvals[u] and tvals[u] < 365:
                [newrow(layout, v[0], self, v[1]) for v in (('End day {}:'.format(u+1), 't'+str(u+1)), ('Fors:', 'f'+str(u+1)), ('Untils:', 'u'+str(u+1)))]
                u += 1

    def update(self):
        for sock in self.outputs:
            socklink(sock, self.id_data.name)

        self.id_data.interface_update(bpy.context)

    def epwrite(self, name, stype):
        if self.source == '0':
            schedtext, ths = '', []

            for tosock in [link.to_socket for link in self.outputs['Schedule'].links]:
                if not schedtext:
                    for t in (self.t1, self.t2, self.t3, self.t4):
                        ths.append(t)
                        if t == 365:
                            break

                    fos = [fs for fs in (self.f1, self.f2, self.f3, self.f4) if fs]
                    uns = [us for us in (self.u1, self.u2, self.u3, self.u4) if us]
                    ts, fs, us = rettimes(ths, fos, uns)
                    schedtext = epschedwrite(name, stype, ts, fs, us)
            return schedtext
        else:
            params = ('Name', 'ScheduleType', 'Name of File', 'Column Number', 'Rows to Skip at Top', 'Number of Hours of Data', 'Column Separator')
            paramvs = (name, 'Any number', bpy.path.abspath(self.select_file), self.cn, self.rtsat, 8760, self.delim)
            schedtext = epentry('Schedule:File', params, paramvs)
            return schedtext

    def epwrite_sel_file(self, name):
        params = ('Name', 'ScheduleType', 'Name of File', 'Column Number', 'Rows to Skip at Top', 'Number of Hours of Data', 'Column Separator')
        paramvs = (name, 'Any number', os.path.abspath(self.select_file), self.cn, self.rtsat, 8760, self.delim)
        schedtext = epentry('Schedule:File', params, paramvs)
        return schedtext

    def epwrite_gen_file(self, name, data, newdir):
        schedtext, ths = '', []

        for tosock in [link.to_socket for link in self.outputs['Schedule'].links]:
            if not schedtext:
                for t in (self.t1, self.t2, self.t3, self.t4):
                    ths.append(t)
                    if t == 365:
                        break

                fos = [fs for fs in (self.f1, self.f2, self.f3, self.f4) if fs]
                uns = [us for us in (self.u1, self.u2, self.u3, self.u4) if us]
                ts, fs, us = rettimes(ths, fos, uns)

        for t in ts:
            for f in fs:
                for u in us:
                    for hi, h in enumerate((datetime.datetime(2015, 1, 1, 0, 00) - datetime.datetime(2014, 1, 1, 0, 00)).hours):
                        if h.day <= self.ts:
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
    bl_icon = 'TEXT'

    text_file: StringProperty(description="Textfile to show")

    def init(self, context):
        self.outputs.new('So_En_Net_Sense', 'Sensor')
        self.outputs.new('So_En_Net_Act', 'Actuator')
        nodecolour(self, 1)

    def draw_buttons(self, context, layout):
        layout.prop_search(self, 'text_file', bpy.data, 'texts', text='File', icon='TEXT')

    def update(self):
        for sock in self.outputs:
            socklink(sock, self.id_data.name)

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
            pparams = ['Name'] + ['line{}'.format(l) for l, line in enumerate(bpy.data.texts[self.text_file.lstrip()].lines) if line.body and line.body.strip()[0] != '!']
            pparamvs = ['{}_controller'.format(self.name.replace(' ', '_'))] + [line.body.strip() for line in bpy.data.texts[self.text_file.lstrip()].lines if line.body and line.body.strip()[0] != '!']
            pentry = epentry('EnergyManagementSystem:Program', pparams, pparamvs)
            return sentries + aentries + cmentry + pentry


class No_En_Net_EMSZone(Node, EnViNodes):
    '''Node describing a Energy management System routine'''
    bl_idname = 'No_En_Net_EMSZone'
    bl_label = 'EMS Zone'
    bl_icon = 'LIGHTPROBE_CUBEMAP'

    def zonelist(self, context):
        return [(c.name, c.name, c.name) for c in bpy.data.collections['EnVi Geometry'].children]

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
        except Exception:
            self.inputs[0].hide = True
            nodecolour(self, 1)

        for iname in [inputs for inputs in self.inputs if inputs.name not in sssocklist and inputs.bl_idname == 'So_En_Net_Act']:
            try: self.inputs.remove(iname)
            except Exception: pass

        for sock in sorted(set(sssocklist)):
            if not self.inputs.get(sock):
                try:
                    self.inputs.new('So_En_Net_Act', sock).sn = sock.split('_')[0] + '-' + '_'.join(sock.split('_')[1:-1])
                except Exception as e: print('3190', e)

#    emszone: StringProperty(name = '', update = zupdate)
    emszone: EnumProperty(name="", description="Zone name", items=zonelist, update = zupdate)
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

class No_En_Net_EMSPy(Node, EnViNodes):
    '''Node identifying a Python interface to the EMS'''
    bl_idname = 'No_En_Net_EMSPy'
    bl_label = 'EMS'
    bl_icon = 'FILE_TEXT'

    py_mod: StringProperty(description="Python file", options = {'SKIP_SAVE'})
    py_class: StringProperty(name = '', description="Textfile to show")

    def init(self, context):
        self.outputs.new('So_Anim', 'Parameter')

    def draw_buttons(self, context, layout):
        layout.prop_search(self, 'py_mod', bpy.data, 'texts', text='Module', icon='TEXT')
        newrow(layout, "Class:", self, 'py_class')

    def epwrite(self):
        for p in self.bl_rna.properties:
            if p.is_skip_save and p.identifier in [l.to_node.parameter for l in self.outputs['Parameter'].links]:
                for l in self.outputs['Parameter'].links:
                    if p.identifier == l.to_node.parameter and l.to_node.anim_file:
                        tf = bpy.data.texts[l.to_node.anim_file]
    #                    paradict[p.identifier] = tf.as_string().split('\n')[bpy.context.scene.frame_current]
                        setattr(self, p.identifier, ret_param(getattr(self, p.identifier), tf.as_string().split('\n')[bpy.context.scene.frame_current - bpy.context.scene.vi_params['enparams']['fs']]))

        if self.py_mod and self.py_class:
            pysparams = ('Name', 'Add Current Working Directory to Search Path', 'Add Input File Directory to Search Path', 'Search Path N')
            pysparamvs = ('Module_search_path', 'No', 'Yes', '')
            pyparams = ('Name', 'Run During Warmup Days', 'Python Module Name', 'Plugin Class Name')
            pyparamvs = ('Python_module', 'Yes', self.py_mod, self.py_class)
            return epentry('PythonPlugin:SearchPaths', pysparams, pysparamvs) + epentry('PythonPlugin:Instance', pyparams, pyparamvs)
        else:
            return ''

class No_En_Net_Anim(Node, EnViNodes):
    '''Node to automate changes in parameters'''
    bl_idname = 'No_En_Net_Anim'
    bl_label = 'VI Animation'
    bl_icon = 'ANIM'

    def retparams(self, context):
        if self.inputs[0].links:
            return [(p.identifier, p.description, p.identifier) for p in self.inputs[0].links[0].from_node.bl_rna.properties if p.is_skip_save]
        else:
            return [('None', 'None', 'None')]

    parameter: EnumProperty(name='', description = 'Parameter to be animated', items=retparams)
    anim_file: StringProperty(name = '')

    def init(self, context):
        self.inputs.new('So_Anim', 'Parameter')

    def draw_buttons(self, context, layout):
        newrow(layout, "Parameter:", self, 'parameter')
        layout.prop_search(self, 'anim_file', bpy.data, 'texts', text='File', icon='TEXT')

class EnViNodeCategory(NodeCategory):
    @classmethod
    def poll(cls, context):
        return context.space_data.tree_type == 'EnViN'

envi_zone = [NodeItem("No_En_Net_Zone", label="Zone"), NodeItem("No_En_Net_Occ", label="Occupancy"),
             NodeItem("No_En_Net_Hvac", label="HVAC"), NodeItem("No_En_Net_Eq", label="Equipment"),
             NodeItem("No_En_Net_Inf", label="Infiltration"), NodeItem("No_En_Net_TC", label="Thermal Chimney")]

envi_sched = [NodeItem("No_En_Net_Sched", label="Schedule Net")]

envi_airflow = [NodeItem("No_En_Net_SFlow", label="Surface Flow"), NodeItem("No_En_Net_SSFlow", label="Sub-surface Flow"),
                NodeItem("No_En_Net_Ext", label="External Air"), NodeItem("No_En_Net_WPC", label="WPC")]

envi_ems = [NodeItem("No_En_Net_EMSZone", label="EMS Zone"), NodeItem("No_En_Net_Prog", label="EMS Program"),
            NodeItem("No_En_Net_EMSPy", label="EMS Python")]

envi_para = [NodeItem("No_En_Net_Anim", label="Parametric")]

envinode_categories = [EnViNodeCategory("Zone", "Zone Nodes", items=envi_zone),
                       EnViNodeCategory("Schedule_Net", "Schedule Nodes", items=envi_sched),
                       EnViNodeCategory("Airflow", "Airflow Nodes", items=envi_airflow),
                       EnViNodeCategory("EMS", "EMS Nodes", items=envi_ems),
                       EnViNodeCategory("EnNPara", "Parametric Nodes", items=envi_para)]

class EnViMatNodes:
    @classmethod
    def poll(cls, ntree):
        return ntree.bl_idname == 'EnViMatN'

class EnViMatNetwork(NodeTree):
    '''A node tree for the creation of EnVi materials.'''
    bl_idname = 'EnViMatN'
    bl_label = 'EnVi Material'
    bl_icon = 'IMGDISPLAY'
#    nodetypes = {}

class EnViMatNodeCategory(NodeCategory):
    @classmethod
    def poll(cls, context):
        return context.space_data.tree_type == 'EnViMatN'

class So_En_Mat_Sched(NodeSocket):
    '''Fraction schedule socket'''
    bl_idname = 'So_En_Mat_Sched'
    bl_label = 'Schedule socket'
    bl_color = (1.0, 1.0, 0.0, 0.75)

    valid = ['Schedule']
    schedule = ['Fraction']

    def draw(self, context, layout, node, text):
        layout.label(text = text)

    def draw_color(self, context, node):
        return (1.0, 1.0, 0.0, 0.75)

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

    valid = ['OLayer', 'TLayer', 'ScreenLayer', 'BlindLayer', 'ShadeLayer']

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
        layout.label(text=text)

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

    valid = ['Shade']

    def draw(self, context, layout, node, text):
        layout.label(text = text)

    def draw_color(self, context, node):
        return (0, 0, 0, 1.0)

    def ret_valid(self, node):
        return ['Shade']

class So_En_Mat_Sc(NodeSocket):
    '''EnVi screen layer socket'''
    bl_idname = 'So_En_Mat_Sc'
    bl_label = 'External screen layer socket'

    valid = ['ScreenLayer']

    def draw(self, context, layout, node, text):
        layout.label(text = text)

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
        layout.label(text = text)

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
        layout.label(text = text)

    def draw_color(self, context, node):
        return (0, 0, 0, 1.0)

    def ret_valid(self, node):
        return ['SControl']

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
    bl_icon = 'NODE_COMPOSITING'

    def envi_con_list(self, context):
        if not self.ec.updated:
            self.ec.update()

        return [(mat, mat, 'Construction') for mat in ((self.ec.wall_con, self.ec.iwall_con)[self.envi_con_con in ("Zone", "Thermal mass")],
                                                        (self.ec.roof_con, self.ec.ceil_con)[self.envi_con_con in ("Zone", "Thermal mass")],
                                                        (self.ec.floor_con, self.ec.ifloor_con)[self.envi_con_con in ("Zone", "Thermal mass")],
                                                        self.ec.door_con, self.ec.glaze_con, self.ec.pv_con)[("Wall", "Roof", "Floor", "Door", "Window", "PV").index(self.envi_con_type)]]

    def con_update(self, context):
        if not self.em.updated:
            self.em.update()

        if not self.ec.updated:
            self.ec.update()

        if len(self.inputs) == 4:
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

            [link.from_node.update() for link in self.inputs['Outer layer'].links]
            #get_mat(self, 0).vi_params.envi_type = self.envi_con_type
            self.pv_update()
            self.update()

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

    def con_type(self, ect):
        if self.envi_con_type == 'Wall':
            envi_con_type = 'Internal wall' if self.envi_con_con == 'Zone' else 'Wall'
        elif self.envi_con_type == 'Floor':
            envi_con_type = 'Internal floor' if self.envi_con_con == 'Zone' else 'Floor'
        elif self.envi_con_type == 'Roof':
            envi_con_type = 'Ceiling' if self.envi_con_con == 'Zone' else 'Roof'
        else:
            envi_con_type = ect

        return envi_con_type

    def uv_update(self, context):
        pstcs, resists = [], []

        if self.envi_con_type in ('Wall', 'Floor', 'Roof'):
            if self.envi_con_makeup == '0':
                con_layers = self.ec.propdict[self.con_type(self.envi_con_type)][self.envi_con_list]
                thicks = [0.001 * tc for tc in [self.lt0, self.lt1, self.lt2, self.lt3,
                                                self.lt4, self.lt5, self.lt6, self.lt7, self.lt8, self.lt9][:len(con_layers)]]

                for p, psmat in enumerate(con_layers):
                    pi = 2 if psmat in self.em.gas_dat else 1
                    pstcs.append(float(self.em.matdat[psmat][pi]))
                    resists.append((thicks[p]/float(self.em.matdat[psmat][pi]), float(self.em.matdat[psmat][pi]))[self.em.matdat[psmat][0] == 'Gas'])

                uv = 1/(sum(resists) + 0.12 + 0.08)
                self.cuv = '{:.3f}'.format(uv)
        else:
            self.cuv = "N/A"

    def frame_update(self, context):
        if self.fclass == '2':
            self.inputs['Outer frame layer'].hide = False
        else:
            for link in self.inputs['Outer frame layer'].links:
                self.id_data.links.unlink(link)

            self.inputs['Outer frame layer'].hide = True


    con_name: StringProperty(name="", description="", default='')
    envi_con_type: EnumProperty(items=[("Wall", "Wall", "Wall construction"),
                                       ("Floor", "Floor", "Ground floor construction"),
                                       ("Roof", "Roof", "Roof construction"),
                                       ("Window", "Window", "Window construction"),
                                       ("Door", "Door", "Door construction"),
                                       ("Shading", "Shading", "Shading material"),
                                       ("None", "None", "Surface to be ignored")],
                                name="",
                                description="Specify the construction type",
                                default="None", update=con_update)
    envi_con_makeup: EnumProperty(items=[("0", "Pre-set", "Construction pre-set"),
                                            ("1", "Layers", "Custom layers")],
                                            name="",
                                            description="Pre-set construction of custom layers",
                                            default="0", update=con_update)
    envi_con_con: EnumProperty(items=bc_update,
                                            name="",
                                            description="Construction context", update=con_update)
    envi_simple_glazing: BoolProperty(name="", description="Flag to signify whether to use a EP simple glazing representation", default=False)
    envi_sg_uv: FloatProperty(name="W/m^2.K", description="Window U-Value", min=0.01, max=10, default=2.4)
    envi_sg_shgc: FloatProperty(name="", description="Window Solar Heat Gain Coefficient", min=0, max=1, default=0.7)
    envi_sg_vt: FloatProperty(name="", description="Window Visible Transmittance", min=0, max=1, default=0.8)
    envi_afsurface: BoolProperty(name="", description="Flag to signify whether the material represents an airflow surface", default=False)
    lt0: FloatProperty(name="mm", description="Layer thickness (mm)", min=0.1, default=100, update=uv_update)
    lt1: FloatProperty(name="mm", description="Layer thickness (mm)", min=0.1, default=100, update=uv_update)
    lt2: FloatProperty(name="mm", description="Layer thickness (mm)", min=0.1, default=100, update=uv_update)
    lt3: FloatProperty(name="mm", description="Layer thickness (mm)", min=0.1, default=100, update=uv_update)
    lt4: FloatProperty(name="mm", description="Layer thickness (mm)", min=0.1, default=100, update=uv_update)
    lt5: FloatProperty(name="mm", description="Layer thickness (mm)", min=0.1, default=100, update=uv_update)
    lt6: FloatProperty(name="mm", description="Layer thickness (mm)", min=0.1, default=100, update=uv_update)
    lt7: FloatProperty(name="mm", description="Layer thickness (mm)", min=0.1, default=100, update=uv_update)
    lt8: FloatProperty(name="mm", description="Layer thickness (mm)", min=0.1, default=100, update=uv_update)
    lt9: FloatProperty(name="mm", description="Layer thickness (mm)", min=0.1, default=100, update=uv_update)
    cuv: StringProperty(name="", description="Construction U-Value", default="N/A")
    cec: StringProperty(name="", description="Construction embodied carbon", default="N/A")
    cecy: StringProperty(name="", description="Construction embodied carbon per year", default="N/A")
    frame_cuv: StringProperty(name="", description="Frame U-Value", default="N/A")
    envi_con_list: EnumProperty(items=envi_con_list, name="", description="Database construction")
    active: BoolProperty(name="", description="Active construction", default=False, update=active_update)
    em = envi_materials()
    ec = envi_constructions()

    # Frame parameters
    fclass: EnumProperty(items=[("0", "Simple spec.", "Simple frame designation"),
                                   ("1", "Detailed spec.", "Advanced frame designation"),
                                   ("2", "Layers", "Layered frame designation")],
                                    name="",
                                    description="Window frame specification",
                                    default="0", update=frame_update)

    fmat: EnumProperty(items=[("0", "Wood", "Wooden frame"),
                                   ("1", "Aluminium", "Aluminium frame"),
                                   ("2", "Plastic", "uPVC frame"),
                                   ("3", "Layers", "Layered frame")],
                                    name="",
                                    description="Frame material",
                                    default="0", update=con_update)

    fthi: FloatProperty(name="m", description="Frame thickness", min=0.001, max=10, default=0.05)
    farea: FloatProperty(name="%", description="Frame area percentage", min=0.01, max=100, default=10)
    fw: FloatProperty(name="m", description="Frame Width", min=0.0, max=10, default=0.2)
    fop: FloatProperty(name="m", description="Frame Outside Projection", min=0.01, max=10, default=0.1)
    fip: FloatProperty(name="m", description="Frame Inside Projection", min=0.01, max=10, default=0.1)
    ftc: FloatProperty(name="W/m.K", description="Frame Conductance", min=0.01, max=10, default=0.1)
    fratio: FloatProperty(name="", description="Ratio of Frame-Edge Glass Conductance to Center-Of-Glass Conductance", min=0.1, max=10, default=1)
    fsa: FloatProperty(name="", description="Frame Solar Absorptance", min=0.01, max=1, default=0.7)
    fva: FloatProperty(name="", description="Frame Visible Absorptance", min=0.01, max=1, default=0.7)
    fte: FloatProperty(name="", description="Frame Thermal Emissivity", min=0.01, max=1, default=0.7)
    dt: EnumProperty(items=[("0", "None", "None"), ("1", "DividedLite", "Divided lites"), ("2", "Suspended", "Suspended divider")],
                                        name="", description="Type of divider", default="0")
    dw: FloatProperty(name="m", description="Divider Width", min=0.001, max=10, default=0.01)
    dhd: IntProperty(name="", description="Number of Horizontal Dividers", min=0, max=10, default=0)
    dvd: IntProperty(name="", description="Number of Vertical Dividers", min=0, max=10, default=0)
    dop: FloatProperty(name="m", description="Divider Outside Projection", min=0.0, max=10, default=0.01)
    dip: FloatProperty(name="m", description="Divider Inside Projection", min=0.0, max=10, default=0.01)
    dtc: FloatProperty(name="W/m.K", description="Divider Conductance", min=0.001, max=10, default=0.1)
    dratio: FloatProperty(name="", description="Ratio of Divider-Edge Glass Conductance to Center-Of-Glass Conductance", min=0.1, max=10, default=1)
    dsa: FloatProperty(name="", description="Divider Solar Absorptance", min=0.01, max=1, default=0.7)
    dva: FloatProperty(name="", description="Divider Visible Absorptance", min=0.01, max=1, default=0.7)
    dte: FloatProperty(name="", description="Divider Thermal Emissivity", min=0.01, max=1, default=0.7)
    orsa: FloatProperty(name="", description="Outside Reveal Solar Absorptance", min=0.01, max=1, default=0.7)
    isd: FloatProperty(name="m", description="Inside Sill Depth (m)", min=0.0, max=10, default=0.1)
    issa: FloatProperty(name="", description="Inside Sill Solar Absorptance", min=0.01, max=1, default=0.7)
    ird: FloatProperty(name="m", description="Inside Reveal Depth (m)", min=0.0, max=10, default=0.1)
    irsa: FloatProperty(name="", description="Inside Reveal Solar Absorptance", min=0.01, max=1, default=0.7)
    resist: FloatProperty(name="", description="U-value", min=0.01, max=10, default=0.7)

    def init(self, context):
        self.inputs.new('So_En_Mat_PV', 'PV')
        self.inputs['PV'].hide = True
        self.inputs.new('So_En_Mat_Ou', 'Outer layer')
        self.inputs['Outer layer'].hide = True
        self.inputs.new('So_En_Mat_Sched', 'Schedule')
        self.inputs['Schedule'].hide = True
        self.inputs.new('So_En_Mat_Fr', 'Outer frame layer')
        self.inputs['Outer frame layer'].hide = True

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

                        for l, layername in enumerate(self.ec.propdict[con_type][self.envi_con_list]):
                            row = layout.row()

                            if layername in self.em.wgas_dat:
                                row.label(text = '{} ({})'.format(layername, "14mm"))
                                row.prop(self, "lt{}".format(l))

                            elif layername in self.em.gas_dat:
                                row.label(text = '{} ({})'.format(layername, "20-50mm"))
                                row.prop(self, "lt{}".format(l))

                            elif layername in self.em.glass_dat:
                                row.label(text = '{} ({})'.format(layername, "{}mm".format(float(self.em.matdat[layername][3])*1000)))
                                row.prop(self, "lt{}".format(l))

                            else:
                                row.label(text = '{} ({})'.format(layername, "{}mm".format(self.em.matdat[layername][7])))
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

                elif self.fclass == '2':
                    newrow(layout, '% area:', self, "farea")
                    row = layout.row()
                    row.operator('node.envi_uv', text = "UV Calc")

                    try:
                        row.label(text = 'U-value  = {} W/m2.K'.format(self.frame_cuv))
                    except Exception:
                        row.label(text = 'U-value  = N/A')

                if self.envi_con_makeup == '1':
                    row = layout.row()
                    row.operator('node.envi_ec', text = "EC Calc")

                    try:
                        row.label(text = 'EC  = {} kgCO2e/m2'.format(self.cec))
                    except Exception:
                        row.label(text = 'EC  = N/A')


            elif self.envi_con_type in ('Wall', 'Floor', 'Roof'):
                row = layout.row()

                if self.envi_con_makeup == '0':
                    try:
                        row.label(text = 'U-value  = {} W/m2.K'.format(self.cuv))
                    except Exception:
                        row.label(text = 'U-value  = N/A')

                elif self.envi_con_makeup == '1':
                    row.operator('node.envi_uv', text = "UV Calc")

                    try:
                        row.label(text = 'U-value  = {} W/m2.K'.format(self.cuv))
                    except Exception:
                        row.label(text = 'U-value  = N/A')

                    row = layout.row()
                    row.operator('node.envi_ec', text = "EC Calc")

                    try:
                        row.label(text = 'EC  = {} kgCO2e/m2'.format(self.cec))
                    except Exception:
                        row.label(text = 'EC  = N/A')

        if self.envi_con_makeup == '1' and self.envi_con_type != 'Shading':
            newrow(layout, 'Name:', self, "con_name")

            if self.con_name and self.inputs['Outer layer'].links:
                row = layout.row()
                row.operator('node.con_save', text = "Save")

    def update(self):
        if not self.em.updated:
            self.em.update()

        if len(self.inputs) == 4:
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
                if lsock.links[0].from_node.bl_idname not in ("No_En_Mat_Sh",):
                    resists.append(lsock.links[0].from_node.ret_resist())
                lsock = lsock.links[0].from_node.inputs['Layer']

            self.cuv = '{:.3f}'.format(1/(sum(resists) + 0.12 + 0.08))
        return self.cuv

    def ret_frame_uv(self):
        if self.envi_con_type == 'Window' and self.fclass == '2':
            resists = []
            lsock = self.inputs['Outer frame layer']

            while lsock.links:
                resists.append(lsock.links[0].from_node.ret_resist())
                lsock = lsock.links[0].from_node.inputs['Layer']

            self.frame_cuv = '{:.3f}'.format(1/(sum(resists) + 0.12 + 0.08))
        return self.frame_cuv

    def ret_ec(self):
        if self.envi_con_makeup == '1':
            ecss, ecys = [], []
            lsock = self.inputs['Outer layer']

            while lsock.links:
                l_node = lsock.links[0].from_node

                if l_node.embodied:
                    node_ec = l_node.ret_ec()
                    ecss.append(node_ec[0])
                    ecys.append(node_ec[1])

                lsock = lsock.links[0].from_node.inputs['Layer']

            self.cec = '{:.3f}'.format(sum(ecss)) if ecss else 'N/A'
            self.cecy = '{:.3f}'.format(sum(ecys)) if ecys else 'N/A'

        return (self.cec, self.cecy)

    def ret_frame_ec(self):
        if self.envi_con_type == 'Window' and self.fclass == '2':
            ecs = []
            lsock = self.inputs['Outer frame layer']

            while lsock.links:
                l_node = lsock.links[0].from_node

                if l_node.embodied:
                    ecs.append(lsock.links[0].from_node.ret_ec())

                lsock = lsock.links[0].from_node.inputs['Layer']

            self.frame_cec = '{:.3f}'.format(sum(ecs))

        return self.frame_cec

    def ret_nodes(self):
        nodes = [self]
        lsock = self.inputs['Outer layer']

        while lsock.links:
            nodes.append(lsock.links[0].from_node)
            lsock = lsock.links[0].from_node.inputs['Layer']

        return nodes

    def save_condict(self):
        if not self.ec.updated:
            self.ec.update()

        self.em.update()

        lks = self.inputs['Outer layer'].links
        lay_names = [lks[0].from_node.lay_name] if lks[0].from_node.layer == '1' else [lks[0].from_node.material]

        while lks:
            lks = lks[0].from_node.inputs[0].links

            if lks:
                lay_name = lks[0].from_node.lay_name if lks[0].from_node.layer == '1' else lks[0].from_node.material
                lay_names.append(lay_name)


        self.ec.get_dat('{} - {}'.format(self.envi_con_type, self.envi_con_con))[self.con_name] = lay_names
        self.ec.get_dat('{} - {}'.format(self.envi_con_type, self.envi_con_con))['{} (reversed)'.format(self.con_name)] = lay_names[::-1]
        self.ec.con_save()

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

    def ep_write(self, mn, ln):
        self['matname'] = get_mat(self, 1).name

        con_type = {'Roof': 'Ceiling', 'Floor': 'Internal floor', 'Wall': 'Internal wall'}[self.envi_con_type] if self.envi_con_con in ('Thermal mass', 'Zone') and self.envi_con_type in ('Roof', 'Wall', 'Floor') else self.envi_con_type

        if not self.ec.updated:
            self.ec.update()

        if not self.em.updated:
            self.em.update()

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
                mats = self.ec.propdict[con_type][self.envi_con_list]
                params = ['Name', 'Outside layer'] + ['Layer {}'.format(i + 1) for i in range(len(mats) - 1)]
                paramvs = [mn] + ['{}-layer-{}'.format(ln, mi) for mi, m in enumerate(mats)]
                ep_text = epentry('Construction', params, paramvs)

                for pm, presetmat in enumerate(mats):
                    matlist = list(self.em.matdat[presetmat])
                    layer_name = '{}-layer-{}'.format(mn, pm)

                    if self.em.namedict.get(presetmat) == None:
                        self.em.namedict[presetmat] = 0
                        self.em.thickdict[presetmat] = [self.thicklist[pm] * 0.001]
                    else:
                        self.em.namedict[presetmat] = self.em.namedict[presetmat] + 1
                        self.em.thickdict[presetmat].append(self.thicklist[pm] * 0.001)

                    if self.envi_con_type in ('Wall', 'Floor', 'Roof', 'Ceiling', 'Door') and presetmat not in self.em.gas_dat:
                        self.resist += self.thicklist[pm]* 0.001/float(matlist[1])
                        params = ('Name', 'Roughness', 'Thickness (m)', 'Conductivity (W/m-K)', 'Density (kg/m3)', 'Specific Heat Capacity (J/kg-K)', 'Thermal Absorptance', 'Solar Absorptance', 'Visible Absorptance')
                        paramvs = ['{}-layer-{}'.format(mn, pm), matlist[0], str(self.thicklist[pm] * 0.001)] + matlist[1:8]
                        ep_text += epentry("Material", params, paramvs)

                        if presetmat in self.em.pcmd_datd:
                            stringmat = self.em.pcmd_datd[presetmat]
                            params = ('Name', 'Temperature Coefficient for Thermal Conductivity (W/m-K2)')
                            paramvs = ('{}-layer-{}'.format(self['matname'], pm), stringmat[0])

                            for i, te in enumerate(stringmat[1].split()):
                                params += ('Temperature {} (C)'.format(i), 'Enthalpy {} (J/kg)'.format(i))
                                paramvs +=(te.split(':')[0], te.split(':')[1])

                            ep_text += epentry("MaterialProperty:PhaseChange", params, paramvs)
                            pcmparams = ('Name', 'Algorithm', 'Construction Name')
                            pcmparamsv = ('{} CondFD override'.format(mn), 'ConductionFiniteDifference', mn)
                            ep_text += epentry('SurfaceProperty:HeatTransferAlgorithm:Construction', pcmparams, pcmparamsv)

                    elif presetmat in self.em.gas_dat:
                        params = ('Name', 'Resistance')
                        paramvs = ('{}-layer-{}'.format(mn, pm), matlist[2])
                        ep_text += epentry("Material:AirGap", params, paramvs)

                    elif self.envi_con_type =='Window':
                        if self.em.matdat[presetmat][0] == 'Glazing':
                            params = ('Name', 'Optical Data Type', 'Window Glass Spectral Data Set Name', 'Thickness (m)', 'Solar Transmittance at Normal Incidence', 'Front Side Solar Reflectance at Normal Incidence',
                          'Back Side Solar Reflectance at Normal Incidence', 'Visible Transmittance at Normal Incidence', 'Front Side Visible Reflectance at Normal Incidence', 'Back Side Visible Reflectance at Normal Incidence',
                          'Infrared Transmittance at Normal Incidence', 'Front Side Infrared Hemispherical Emissivity', 'Back Side Infrared Hemispherical Emissivity', 'Conductivity (W/m-K)',
                          'Dirt Correction Factor for Solar and Visible Transmittance', 'Solar Diffusing')
                            paramvs = ['{}-layer-{}'.format(mn, pm)] + matlist[1:3] + [self.thicklist[pm] * 0.001] + ['{:.3f}'.format(float(sm)) for sm in matlist[4:-1]] + [1, ('No', 'Yes')[matlist[-1]]]
                            ep_text += epentry("WindowMaterial:{}".format(matlist[0]), params, paramvs)

                        elif self.em.matdat[presetmat][0] == 'Gas':
                            params = ('Name', 'Gas Type', 'Thickness')
                            paramvs = [layer_name] + [matlist[1]] + [self.thicklist[pm] * 0.001]
                            ep_text += epentry("WindowMaterial:Gas", params, paramvs)

        elif self.envi_con_makeup == '1':
            in_sock = self.inputs['Outer layer']# if self.envi_con_type == "Window" else self.inputs[0]
            n = 0
            params = ['Name']
            paramvs = [mn]
            ep_text = ''
            self.resist, ecm2 = 0, 0
            get_mat(self, 1).vi_params.envi_shading = 0

            while in_sock.links:
                node = in_sock.links[0].from_node
                paramvs.append('{}-layer-{}'.format(ln, n))
                params.append(('Outside layer', 'Layer {}'.format(n))[n > 0])
                ep_text += node.ep_write(n, mn)
                self.resist += node.resist

                if node.get('embodied') and node['ecm2']:
                    ecm2 += float(node['ecm2'])

                if (node.inputs.get('Shade') and (node.inputs['Shade'].links or node.outputs['Shade'].links)) or (node.outputs.get('Shade') and (node.outputs['Shade'].links or node.inputs['Shade'].links)):
                    get_mat(self, 1).vi_params.envi_shading = 1

                    if node.inputs['Layer'].links and node.inputs['Layer'].links[0].from_node.bl_idname == 'No_En_Mat_Gas' and node.inputs['Shade'].links and node.inputs['Shade'].links[0].from_node.bl_idname != 'No_En_Mat_SG':
                        g_t = node.inputs['Layer'].links[0].from_node.thi
                        s_t = node.inputs['Shade'].links[0].from_node.thi
                        ep_text += node.inputs['Layer'].links[0].from_node.ep_write(n + 1, mn + '_split', tmod = (g_t - s_t)/(2 * g_t))

                    elif node.outputs['Layer'].links and node.outputs['Layer'].links[0].to_node.bl_idname == 'No_En_Mat_Gas' and node.outputs['Shade'].links:
                        g_t = node.outputs['Layer'].links[0].to_node.thi
                        s_t = node.outputs['Shade'].links[0].to_node.thi
                        ep_text += node.outputs['Layer'].links[0].to_node.ep_write(n - 1, mn + '_split', tmod = (g_t - s_t)/(2 * g_t))

                get_mat(self, 1).vi_params['enparams']['ecm2'] = ecm2
                in_sock = node.inputs['Layer']
                n += 1

            if mn == self.id_data.name:
                ep_text += epentry('Construction', params, paramvs)
            else:
                ep_text += epentry('Construction', params, [paramvs[0]] + paramvs[1:][::-1])

            if get_mat(self, 1).vi_params.envi_shading:
                in_sock = self.inputs['Outer layer']
                n = 0
                params = ['Name']
                paramvs = ['{}-shading'.format(mn)]

                while in_sock.links:
                    node = in_sock.links[0].from_node

                    if node.bl_idname == 'No_En_Mat_Tr' and node.inputs['Shade'].links and node.inputs['Shade'].links[0].from_node.bl_idname == 'No_En_Mat_SG':
                        paramvs.append('{}-shading-{}'.format(mn, n))
                        params.append(('Outer shader', 'Shading layer {}'.format(n))[n > 0])
                        ep_text += node.inputs['Shade'].links[0].from_node.ep_write(n, mn)

                    elif node.bl_idname == 'No_En_Mat_Tr' and node.outputs['Layer'].links[0].to_node.bl_idname != 'No_En_Mat_Gas':
                        if node.outputs['Shade'].links:
                            paramvs.append('{}-shading-{}'.format(mn, n))
                            params.append(('Outer shader', 'Shading layer {}'.format(n))[n > 0])
                            ep_text += node.outputs['Shade'].links[0].to_node.ep_write(n, mn)

                        paramvs.append('{}-layer-{}'.format(mn, n))
                        params.append(('Outside layer', 'Layer {}'.format(n))[n > 0])

                        if node.inputs['Shade'].links:
                            paramvs.append('{}-shading-{}'.format(mn, n))
                            params.append(('Outer shader', 'Shading layer {}'.format(n))[n > 0])
                            ep_text += node.inputs['Shade'].links[0].from_node.ep_write(n, mn)

                    elif node.bl_idname == 'No_En_Mat_Gas':
                        shade_bool = len(node.outputs['Layer'].links[0].to_node.inputs['Shade'].links) or (node.inputs['Layer'].links and len(node.inputs['Layer'].links[0].from_node.outputs['Shade'].links))
                        if shade_bool:
                            if node.outputs['Layer'].links[0].to_node.inputs['Shade'].links and node.outputs['Layer'].links[0].to_node.inputs['Shade'].links[0].from_node.bl_idname == 'No_En_Mat_SG':
                                shade_bool = 0

                        if not shade_bool:
                            paramvs.append('{}-layer-{}'.format(mn, n))
                            params.append(('Outside layer', 'Layer {}'.format(n))[n > 0])
                        else:
                            paramvs.append('{}-layer-{}'.format(mn + '_split', n))
                            params.append(('Outside layer', 'Layer {}'.format(n))[n > 0])
                            paramvs.append('{}-shading-{}'.format(mn, n))
                            params.append(('Outside layer', 'Layer {}'.format(n))[n > 0])

                            if node.outputs['Layer'].links and node.outputs['Layer'].links[0].to_node.inputs['Shade'].links:
                                shade_node = node.outputs['Layer'].links[0].to_node.inputs['Shade'].links[0].from_node
                            else:
                                shade_node = node.inputs['Layer'].links[0].from_node.outputs['Shade'].links[0].to_node

                            ep_text += shade_node.ep_write(n, mn)
                            paramvs.append('{}-layer-{}'.format(mn + '_split', n))
                            params.append(('Outside layer', 'Layer {}'.format(n))[n > 0])

                    elif not node.inputs['Layer'].links:
                        if not node.inputs['Shade'].links:
                            paramvs.append('{}-layer-{}'.format(mn, n))
                            params.append(('Outside layer', 'Layer {}'.format(n))[n > 0])

                        if node.inputs['Shade'].links:
                            paramvs.append('{}-shading-{}'.format(mn, n))
                            params.append(('Outer shader', 'Shading layer {}'.format(n))[n > 0])
                            ep_text += node.inputs['Shade'].links[0].from_node.ep_write(n, mn)

                    else:
                        paramvs.append('{}-layer-{}'.format(mn, n))
                        params.append(('Outside layer', 'Layer {}'.format(n))[n > 0])

                    n += 1
                    in_sock = node.inputs['Layer']

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
                in_sock = self.inputs['Outer frame layer']
                n = 0
                params = ['Name']
                paramvs = ['{}-frame'.format(mn)]

                while in_sock.links:
                    paramvs.append('{}-frame-layer-{}'.format(mn, n))
                    params.append(('Outside layer', 'Layer {}'.format(n))[n > 0])
                    node = in_sock.links[0].from_node
                    ep_text += node.ep_write(n, mn)
                    in_sock = node.inputs['Layer']
                    n += 1

                ep_text += epentry('Construction', params, paramvs)

        return ep_text

    # def layer_write(self, in_sock, matname):
    #     ep_text = ''
    #     n = 0

    #     while in_sock.links:
    #         node = in_sock.links[0].from_node
    #         paramvs.append('{}-frame-layer-{}'.format(matname, n))
    #         params.append(('Outside layer', 'Layer {}'.format(n))[n > 0])
    #         ep_text += node.ep_write(n, mn)
    #         in_sock = node.inputs['Layer']
    #         n += 1

    #     ep_text += epentry('Construction', params, paramvs)
    #     return ep_text

class No_En_Mat_Op(Node, EnViMatNodes):
    '''Node defining the EnVi opaque material layer'''
    bl_idname = 'No_En_Mat_Op'
    bl_label = 'EnVi opaque layer'

    def lay_update(self, context):
        if not self.em.updated:
            self.em.update()

        c_vals = self.em.get_dat(self.materialtype)[self.material]
        self.rough = c_vals[0]
        self.tc = float(c_vals[1])
        self.rho = float(c_vals[2])
        self.shc = float(c_vals[3])
        self.tab = float(c_vals[4])
        self.sab = float(c_vals[5])
        self.vab = float(c_vals[6])

        if self.layer == '1' and self.lay_name == '':
            nodecolour(self, 1)
        elif self.layer == '0' and not self.material:
            nodecolour(self, 1)
        else:
            nodecolour(self, 0)

    def ec_update(self, context):
        self['ecentries'] = []

        if not self.ee.updated:
            self.ee.update()

        if self.embodiedtype in self.ee.propdict.keys():
            if self.embodiedclass not in self.ee.propdict[self.embodiedtype]:
                self.embodiedclass = list(self.ee.propdict[self.embodiedtype].keys())[0]

            if self.embodiedmat not in self.ee.propdict[self.embodiedtype][self.embodiedclass]:
                self.embodiedmat = list(self.ee.propdict[self.embodiedtype][self.embodiedclass])[0]

            try:
                self['ecdict'] = self.ee.propdict[self.embodiedtype][self.embodiedclass][self.embodiedmat]
                self['ecm2'] = '{:.3f}'.format(float(self['ecdict']['eckg']) * float(self['ecdict']['density']) * self.thi * 0.001)
                self['ecm2y'] = '{:.3f}'.format(float(self['ecdict']['eckg']) * float(self['ecdict']['density']) * self.thi * 0.001/self.ec_life)
                self['ecentries'] = [(k, self.ee.propdict[self.embodiedtype][self.embodiedclass][self.embodiedmat][k]) for k in self['ecdict'].keys()]

            except Exception as e:
                self['ecm2'] = 'N/A'
                self['ecm2y'] = 'N/A'

    lay_name: StringProperty(name='', description='Custom layer name', update=lay_update)
    layer: EnumProperty(items=[("0", "Database", "Select from database"),
                                        ("1", "Custom", "Define custom material properties")],
                                        name="", description="Class of layer", default="0", update=lay_update)

    materialtype: EnumProperty(items=envi_layertype, name="", description="Layer material type", update=lay_update)
    embodied: BoolProperty(name="", description="Embodied carbon", default=0, update=ec_update)
    embodiedtype: EnumProperty(items=envi_elayertype, name="", description="Layer embodied material class", update=ec_update)
    embodiedclass: EnumProperty(items=envi_eclasstype, name="", description="Layer embodied class", update=ec_update)
    embodiedmat: EnumProperty(items=envi_emattype, name="", description="Layer embodied material", update=ec_update)
    ecm2: FloatProperty(name="kgCo2e/m2", description="Embodied carbon per metre squared", min=0.0, default=0.0)
    material: EnumProperty(items=envi_layer, name="", description="Layer material", update=lay_update)
    thi: FloatProperty(name="mm", description="Thickness (mm)", min=0.1, max=10000, default=100, update=ec_update)
    tc: FloatProperty(name="W/m.K", description="Thermal conductivity (W/m.K)", min=0.001, max=10, precision=3, default=0.5)
    rough: EnumProperty(items=[("VeryRough", "VeryRough", "Roughness"),
                                  ("Rough", "Rough", "Roughness"),
                                  ("MediumRough", "MediumRough", "Roughness"),
                                  ("MediumSmooth", "MediumSmooth", "Roughness"),
                                  ("Smooth", "Smooth", "Roughness"),
                                  ("VerySmooth", "VerySmooth", "Roughness")],
                                  name="",
                                  description="Specify the material roughness for convection calculations",
                                  default="Rough")

    rho: FloatProperty(name="kg/m^3", description="Density", min=0.001, max=10000, default=800)
    shc: FloatProperty(name="J/kg", description="Thickness (mm)", min=0.01, max=10000, default=800)
    tab: FloatProperty(name="", description="Thermal absorptance", min=0, max=1, precision=2, default=0.7)
    sab: FloatProperty(name="", description="Solar absorptance", min=0, max=1, precision=2, default=0.7)
    vab: FloatProperty(name="", description="Visual absorptance", min=0, max=1, precision=2, default=0.7)
    pcm: BoolProperty(name="", description="Phase Change Material", default=0)
    tctc: FloatProperty(name="", description="Temp. coeff. for thermal conductivity (W/m-K2)", min=0, max=50, default=0)
    tempemps: StringProperty(name="", description="Temperature/empalthy pairs (e.g. T1:E1 T2:E2)", default="")
    resist: FloatProperty(name="", description="", min=0, default=0)
    envi_con_type: StringProperty(name="", description="Name")
    ec_id: StringProperty(name="", description="Embodied id")
    ec_type: StringProperty(name="", description="Embodied type")
    ec_class: StringProperty(name="", description="Embodied class")
    ec_name: StringProperty(name="", description="Embodied name")
    ec_unit:EnumProperty(items=[("kg", "kg", "per kilogram"),
                                  ("m2", "m2", "per square metre"),
                                  ("m3", "m3", "per cubic metre"),
                                  ("unit", "Unit", "per unit")],
                                  name="",
                                  description="Embodied carbon unit",
                                  default="kg")
    ec_amount: FloatProperty(name="", description="", min=0.001, precision=3, default=1)
    ec_kgco2e: FloatProperty(name="", description="Embodied carbon per kg amount", precision=3, default=100)
    # ec_m2: FloatProperty(name="", description="Embodied carbon per area amount", default=100)
    ec_density: FloatProperty(name="kg/m^3", description="Material density", default=1000)
    ec_mod: StringProperty(name="", description="Embodied modules")
    ec_life: IntProperty(name="Years", description="Service life in years", min=1, max=100, default=60, update=ec_update)
    em = envi_materials()
    ee = envi_embodied()

    def init(self, context):
        self.outputs.new('So_En_Mat_Op', 'Layer')
        self.inputs.new('So_En_Mat_Op', 'Layer')
        self.outputs.new('So_Anim', 'Parameter')

    def draw_buttons(self, context, layout):
        newrow(layout, "Type:", self, "materialtype")
        newrow(layout, "Specification:", self, "layer")

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

            if self.lay_name:
                row = layout.row()
                row.operator('node.lay_save', text="Layer Save")

        newrow(layout, "Embodied:", self, "embodied")

        if self.embodied:
            newrow(layout, "Embodied type:", self, "embodiedtype")

            if self.embodiedtype != 'Custom':
                newrow(layout, "Embodied class:", self, "embodiedclass")
                newrow(layout, "Embodied material:", self, "embodiedmat")
                newrow(layout, "Service life:", self, "ec_life")

                try:
                    for ec in self['ecentries']:
                        row = layout.row()
                        row.label(text='{}: {}'.format(ec[0], ec[1]))

                    row = layout.row()
                    row.label(text='ec/m2: {}'.format(self['ecm2']))
                    row = layout.row()
                    row.label(text='ec/m2/y: {}'.format(self['ecm2y']))

                except Exception:
                    pass

            else:
                newrow(layout, "Embodied id:", self, "ec_id")
                newrow(layout, "Embodied type:", self, "ec_type")
                newrow(layout, "Embodied class:", self, "ec_class")
                newrow(layout, "Embodied name:", self, "ec_name")
                newrow(layout, "Embodied modules:", self, "ec_mod")
                newrow(layout, "Embodied unit:", self, "ec_unit")

                if self.ec_unit in ('kg', 'm2', 'm3'):
                    newrow(layout, "Embodied amount:", self, "ec_amount")

                newrow(layout, "GWP per amount:", self, "ec_kgco2e")
                newrow(layout, "Embodied density:", self, "ec_density")
                newrow(layout, "Service life:", self, "ec_life")

                if all((self.ec_id, self.ec_name, self.ec_type, self.ec_class, self.ec_mod)):
                    row = layout.row()
                    row.operator('node.ec_save', text="Embodied Save")

    def ret_resist(self):
        if self.layer == '0':
            if not self.material:
                self.update()
            matlist = list(self.em.matdat[self.material])

            if self.materialtype != '6':
                self.resist = self.thi * 0.001/float(matlist[1])
            else:
                self.resist = float(matlist[2])
        else:
            self.resist = self.thi * 0.001/self.tc

        return self.resist

    def ret_ec(self):
        self.ec_update(0)
        try:
            return (float(self['ecm2']), float(self['ecm2y']))
        except Exception:
            return (0, 0)

    def save_laydict(self):
        '''Roughness, Conductivity {W/m-K}, Density {kg/m3}, Specific Heat {J/kg-K}, Thermal Absorbtance,
            Solar Absorbtance, Visible Absorbtance, Default thickness'''
        self.em.get_dat(self.materialtype)[self.lay_name] = [self.rough, '{:.4f}'.format(self.tc),
                                                             '{:.2f}'.format(self.rho), '{:.2f}'.format(self.shc),
                                                             '{:.2f}'.format(self.tab), '{:.2f}'.format(self.sab),
                                                             '{:.2f}'.format(self.vab), '{:.2f}'.format(self.thi)]
        if self.materialtype == '8':
            self.em.get_dat('9')[self.lay_name] = [self.tctc, self.tempemps]

        self.em.lay_save()

    def save_ecdict(self):
        '''ID, Quantity, unit, density, weight, ec per unit, ec per kg, modules'''
        if self.ec_unit == 'kg':
            weight = self.ec_amount
            eckg = self.ec_kgco2e/weight
            ecdu = self.ec_kgco2e
        elif self.ec_unit == 'm2':
            weight = self.ec_amount * self.ec_density * self.thi * 0.001
            eckg = self.ec_kgco2e/weight
            ecdu = self.ec_kgco2e
        elif self.ec_unit == 'm3':
            weight = self.ec_amount * self.ec_density
            eckg = self.ec_kgco2e/weight
            ecdu = self.ec_kgco2e
        else:
            weight = ''
            ecdu = self.ec_kgco2e
            eckg = self.ec_kgco2e

        ec_dict = self.ee.get_dat()

        if self.ec_type not in ec_dict:
            ec_dict[self.ec_type] = {}
            ec_dict[self.ec_type][self.ec_class] = {}
            ec_dict[self.ec_type][self.ec_class][self.ec_name] = {"id": self.ec_id,
                                                                  "quantity": '{:.4f}'.format(self.ec_amount),
                                                                  "unit": self.ec_unit,
                                                                  "density": '{:.4f}'.format(self.ec_density),
                                                                  "weight": '{:.4f}'.format(weight),
                                                                  "ecdu": '{:.4f}'.format(ecdu),
                                                                  "eckg": '{:.4f}'.format(eckg),
                                                                  "modules": self.ec_mod}
        elif self.ec_class not in ec_dict[self.ec_type]:
            ec_dict[self.ec_type][self.ec_class] = {}
            ec_dict[self.ec_type][self.ec_class][self.ec_name] = {"id": self.ec_id,
                                                                  "quantity": '{:.4f}'.format(self.ec_amount),
                                                                  "unit": self.ec_unit,
                                                                  "density": '{:.4f}'.format(self.ec_density),
                                                                  "weight": '{:.4f}'.format(weight),
                                                                  "ecdu": '{:.4f}'.format(ecdu),
                                                                  "eckg": '{:.4f}'.format(eckg),
                                                                  "modules": self.ec_mod}
        else:
            ec_dict[self.ec_type][self.ec_class][self.ec_name] = {"id": self.ec_id,
                                                                  "quantity": '{:.4f}'.format(self.ec_amount),
                                                                  "unit": self.ec_unit,
                                                                  "density": '{:.4f}'.format(self.ec_density),
                                                                  "weight": '{:.4f}'.format(weight),
                                                                  "ecdu": '{:.4f}'.format(ecdu),
                                                                  "eckg": '{:.4f}'.format(eckg),
                                                                  "modules": self.ec_mod}
        self.ee.set_dat(ec_dict)
        self.ee.ec_save()

    def update(self):
        if not self.material:
            print('no material')

        for sock in self.outputs:
            socklink2(sock, self.id_data)

        if self.outputs['Layer'].links:
            ect_node = self.outputs['Layer'].links[0].to_node if self.outputs['Layer'].links[0].to_node.envi_con_type else self.outputs['Layer'].links[0].to_node.outputs['Layer'].links[0].to_node
            self.envi_con_type = ect_node.envi_con_type if self.outputs['Layer'].links[0].to_socket.bl_idname != 'So_En_Mat_Fr' else 'Frame'

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
            paramvs += (te.split(':')[0], te.split(':')[1])

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
            matlist = list(self.em.matdat[self.material])

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

    def lay_update(self, context):
        if not self.em.updated:
            self.em.update()

        if self.layer == '1' and self.lay_name == '':
            nodecolour(self, 1)
        elif self.layer == '0' and not self.material:
            nodecolour(self, 1)
        else:
            nodecolour(self, 0)

    def ec_update(self, context):
        self['ecentries'] = []

        if not self.ee.updated:
            self.ee.update()

        if self.embodiedtype in self.ee.propdict.keys():
            if self.embodiedclass not in self.ee.propdict[self.embodiedtype]:
                self.embodiedclass = list(self.ee.propdict[self.embodiedtype].keys())[0]

            if self.embodiedmat not in self.ee.propdict[self.embodiedtype][self.embodiedclass]:
                self.embodiedmat = list(self.ee.propdict[self.embodiedtype][self.embodiedclass])[0]

            try:
                self['ecdict'] = self.ee.propdict[self.embodiedtype][self.embodiedclass][self.embodiedmat]
                self['ecm2'] = '{:.3f}'.format(float(self['ecdict']['eckg']) * float(self['ecdict']['density']) * self.thi * 0.001)
                self['ecm2y'] = '{:.3f}'.format(float(self['ecdict']['eckg']) * float(self['ecdict']['density']) * self.thi * 0.001/self.ec_life)
                self['ecentries'] = [(k, self.ee.propdict[self.embodiedtype][self.embodiedclass][self.embodiedmat][k]) for k in self['ecdict'].keys()]

            except Exception as e:
                self['ecm2'] = 'N/A'
                self['ecm2y'] = 'N/A'

    lay_name: StringProperty(name='', description='Custom layer name')
    layer: EnumProperty(items=[("0", "Database", "Select from database"),
                                        ("1", "Custom", "Define custom material properties")],
                                        name="", description="Composition of the layer", default="0", update=lay_update)
    materialtype: EnumProperty(items=envi_layertype, name="", description="Layer material type")
    material: EnumProperty(items=envi_layer, name="", description="Layer material", update=lay_update)
    thi: FloatProperty(name="mm", description="Thickness (mm)", min=0.1, max=1000, default=6)
    tc: FloatProperty(name="W/m.K", description="Thermal Conductivity (W/m.K)", precision=3, min=0.1, max=10, default=0.8)
    stn: FloatProperty(name="", description="Solar normal transmittance", precision=3, min=0, max=1, default=0.7)
    fsn: FloatProperty(name="", description="Solar front normal reflectance", precision=3, min=0, max=1, default=0.07)
    bsn: FloatProperty(name="", description="Solar back normal reflectance", precision=3, min=0, max=1, default=0.07)
    vtn: FloatProperty(name="", description="Visible Transmittance at Normal Incidence", precision=3, min=0, max=1, default=0.89)
    fvrn: FloatProperty(name="", description="Front Side Visible Reflectance at Normal Incidence", precision=3, min=0, max=1, default=0.07)
    bvrn: FloatProperty(name="", description="Back Side Visible Reflectance at Normal Incidence", precision=3, min=0, max=1, default=0.07)
    itn: FloatProperty(name="", description="Infrared Transmittance at Normal Incidence", precision=3, min=0, max=1, default=0.0)
    fie: FloatProperty(name="", description="Front Side Infrared Hemispherical Emissivity", precision=3, min=0, max=1, default=0.84)
    bie: FloatProperty(name="", description="Back Side Infrared Hemispherical Emissivity", precision=3, min=0, max=1, default=0.84)
    diff: BoolProperty(name="", description="Diffusing", default=0)
    envi_con_type: StringProperty(name="", description="Name")
    resist: FloatProperty(name="", description="", min=0, default=0)
    embodied: BoolProperty(name="", description="Embodied carbon", default=0, update=ec_update)
    embodiedtype: EnumProperty(items=envi_elayertype, name="", description="Layer embodied material class", update=ec_update)
    embodiedclass: EnumProperty(items=envi_eclasstype, name="", description="Layer embodied class", update=ec_update)
    embodiedmat: EnumProperty(items=envi_emattype, name="", description="Layer embodied material", update=ec_update)
    ecm2: FloatProperty(name="kgCo2e/m2", description="Embodied carbon per metre squared", min=0.0, default=0.0)
    ec_id: StringProperty(name="", description="Embodied id")
    ec_type: StringProperty(name="", description="Embodied type")
    ec_class: StringProperty(name="", description="Embodied class")
    ec_name: StringProperty(name="", description="Embodied name")
    ec_unit:EnumProperty(items=[("kg", "kg", "per kilogram"),
                                  ("m2", "m2", "per square metre"),
                                  ("m3", "m3", "per cubic metre"),
                                  ("unit", "Unit", "per unit")],
                                  name="",
                                  description="Embodied carbon unit",
                                  default="kg")
    ec_amount: FloatProperty(name="", description="", min=0.001, precision=3, default=1)
    ec_kgco2e: FloatProperty(name="", description="Embodied carbon per kg amount", precision=3, default=100)
    ec_density: FloatProperty(name="kg/m^3", description="Material density", default=1000)
    ec_life: IntProperty(name="y", description="Service life in years", min=1, max=100, default=60, update=ec_update)
    ec_mod: StringProperty(name="", description="Embodied modules")
    em = envi_materials()
    ee = envi_embodied()

    def init(self, context):
        self.outputs.new('So_En_Mat_Tr', 'Layer')
        self.inputs.new('So_En_Mat_Gas', 'Layer')
        self.outputs.new('So_En_Mat_Sh', 'Shade')
        self.outputs['Shade'].link_limit = 1
        self.inputs.new('So_En_Mat_Sh', 'Shade')

    def draw_buttons(self, context, layout):
        newrow(layout, "Specification:", self, "layer")

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
                row.operator('node.lay_save', text="Layer Save")

        newrow(layout, "Embodied:", self, "embodied")

        if self.embodied:
            newrow(layout, "Embodied type:", self, "embodiedtype")

            if self.embodiedtype != 'Custom':
                newrow(layout, "Embodied class:", self, "embodiedclass")
                newrow(layout, "Embodied material:", self, "embodiedmat")
                newrow(layout, "Service life:", self, "ec_life")

                try:
                    for ec in self['ecentries']:
                        row = layout.row()
                        row.label(text='{}: {}'.format(ec[0], ec[1]))

                    row = layout.row()
                    row.label(text='ec/m2: {}'.format(self['ecm2']))
                    row = layout.row()
                    row.label(text='ec/m2/y: {}'.format(self['ecm2y']))

                except Exception as e:
                    print(e)

            else:
                newrow(layout, "Embodied id:", self, "ec_id")
                newrow(layout, "Embodied type:", self, "ec_type")
                newrow(layout, "Embodied class:", self, "ec_class")
                newrow(layout, "Embodied name:", self, "ec_name")
                newrow(layout, "Embodied modules:", self, "ec_mod")
                newrow(layout, "Embodied unit:", self, "ec_unit")

                if self.ec_unit in ('kg', 'm2', 'm3'):
                    newrow(layout, "Embodied amount:", self, "ec_amount")

                newrow(layout, "GWP per amount:", self, "ec_kgco2e")
                newrow(layout, "Embodied density:", self, "ec_density")
                newrow(layout, "Service life:", self, "ec_life")

                if all((self.ec_id, self.ec_name, self.ec_type, self.ec_class, self.ec_mod)):
                    row = layout.row()
                    row.operator('node.ec_save', text="Embodied Save")

    def update(self):
        for sock in self.outputs:
            socklink2(sock, self.id_data)

        if self.outputs['Layer'].links:
            self.envi_con_type = self.outputs['Layer'].links[0].to_node.envi_con_type

        if self.outputs['Shade'].links and self.outputs['Shade'].links[0].to_node.bl_idname == 'No_En_Mat_Sc' and self.outputs['Layer'].links and self.outputs['Layer'].links[0].to_node.bl_idname != 'No_En_Mat_Con':
            self.id_data.links.remove(self.outputs['Shade'].links[0])

        if self.outputs['Shade'].links and self.outputs['Shade'].links[0].to_node.bl_idname in ('No_En_Mat_Bl', 'No_En_Mat_Sh') and self.inputs['Layer'].links and self.outputs['Layer'].links[0].to_node.bl_idname != 'No_En_Mat_Con':
            self.id_data.links.remove(self.outputs['Shade'].links[0])

        if self.inputs['Shade'].links and self.inputs['Shade'].links[0].from_node.bl_idname in ('No_En_Mat_Bl', 'No_En_Mat_Sh') and self.inputs['Layer'].links and self.inputs['Layer'].links[0].from_node.inputs['Layer'].links:
            if self.inputs['Layer'].links[0].from_node.inputs['Layer'].links[0].from_node.inputs['Layer'].links:
                self.id_data.links.remove(self.inputs['Shade'].links[0])

        # if self.inputs['Shade'].links and self.inputs['Shade'].links[0].from_node.bl_idname == 'No_En_Mat_SG' and self.inputs['Layer'].links:
        #     self.id_data.links.remove(self.inputs['Shade'].links[0])

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
                         '{:.4f}'.format(self.fvrn), '{:.4f}'.format(self.bvrn), '{:.4f}'.format(self.itn), '{:.4f}'.format(self.fie), '{:.4f}'.format(self.bie), int(self.diff)]
        envi_mats.lay_save()

    def ret_resist(self):
        if not self.em.updated:
            self.em.update()

        if self.layer == '0':
            matlist = list(self.em.matdat[self.material])
            self.resist = float(matlist[13])
        else:
            self.resist = self.thi * 0.001/self.tc

        return self.resist

    def ret_ec(self):
        self.ec_update(0)
        try:
            return (float(self['ecm2']), float(self['ecm2y']))
        except Exception:
            return (0, 0)

    def ep_write(self, ln, mn):
        self.em.update()

        for material in bpy.data.materials:
            if self.id_data == material.vi_params.envi_nodes:
                break

        layer_name = '{}-layer-{}'.format(mn, ln)
        params = ('Name', 'Optical Data Type', 'Window Glass Spectral Data Set Name', 'Thickness (m)', 'Solar Transmittance at Normal Incidence', 'Front Side Solar Reflectance at Normal Incidence',
                  'Back Side Solar Reflectance at Normal Incidence', 'Visible Transmittance at Normal Incidence', 'Front Side Visible Reflectance at Normal Incidence', 'Back Side Visible Reflectance at Normal Incidence',
                  'Infrared Transmittance at Normal Incidence', 'Front Side Infrared Hemispherical Emissivity', 'Back Side Infrared Hemispherical Emissivity', 'Conductivity (W/m-K)',
                  'Dirt Correction Factor for Solar and Visible Transmittance', 'Solar Diffusing')

        if self.layer == '0':
            matlist = list(self.em.matdat[self.material])
            paramvs = [layer_name] + matlist[1:3] + [self.thi * 0.001] + ['{:.3f}'.format(float(sm)) for sm in matlist[4:-1]] + [1, ('No', 'Yes')[matlist[-1]]]

        else:
            paramvs = ['{}-layer-{}'.format(mn, ln), 'SpectralAverage', '', self.thi * 0.001, '{:.3f}'.format(self.stn), '{:.3f}'.format(self.fsn), '{:.3f}'.format(self.bsn),
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
    thi: FloatProperty(name = "mm", description = "Thickness (mm)", min = 0.1, max = 5000, default = 14)
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
    embodied: BoolProperty(name = "", description = "Embodied carbon", default = 0)
    em = envi_materials()
    ee = envi_embodied()

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
                newrow(layout, "Viscosity B:", self, "vcB")
                newrow(layout, "Viscosity C:", self, "vcC")
                newrow(layout, "SHC A:", self, "shcA")
                newrow(layout, "SHC B:", self, "shcB")
                newrow(layout, "SHC C:", self, "shcC")
                newrow(layout, "Mol Weight:", self, "mw")
                newrow(layout, "SHR:", self, "shr")

    def update(self):
        for sock in self.outputs:
            socklink(sock, self.id_data.name)

        if self.outputs['Layer'].links:
            # ect_node = self.outputs['Layer'].links[0].to_node if self.outputs['Layer'].links[0].to_node.envi_con_type else self.outputs['Layer'].links[0].to_node.outputs['Layer'].links[0].to_node
            self.envi_con_type = self.outputs['Layer'].links[0].to_node.envi_con_type

        self.valid()

    def valid(self):
        if not self.outputs["Layer"].links or not self.inputs["Layer"].links:
            nodecolour(self, 1)
        else:
            nodecolour(self, 0)

    def ret_resist(self):
        if not self.em.updated:
            self.em.update()

        if self.layer == '0':
            matlist = list(self.em.matdat[self.material])
            self.resist = self.thi * 0.001/float(matlist[4])
        else:
            self.tc = self.ccA + self.ccB * 293.14
            self.resist = self.thi * 0.001/self.tc

        return self.resist

    def ret_ec(self):
        return (0, 0)

    def ep_write(self, ln, mn, **kwargs):
        tmod = kwargs['tmod'] if kwargs.get('tmod') else 1

        if self.layer == '0':
            params = ('Name', 'Gas Type', 'Thickness')
            paramvs = ['{}-layer-{}'.format(mn, ln), self.material, self.thi * tmod * 0.001]

        else:
            params = ('gap name', 'type', 'thickness', 'Conductivity Coefficient A', 'Conductivity Coefficient B', 'Conductivity Coefficient C',
                      'Conductivity Viscosity A', 'Conductivity Viscosity B', 'Conductivity Viscosity C', 'Specific Heat Coefficient A',
                      'Specific Heat Coefficient B', 'Specific Heat Coefficient C', 'Molecular Weight', 'Specific Heat Ratio')
            paramvs = ['{}-layer-{}'.format(mn, ln), 'Custom', '{:.3f}'.format(self.thi * 0.001), '{:.3f}'.format(self.ccA), '{:.3f}'.format(self.ccB),
                       '{:.3f}'.format(self.ccC), '{:.3f}'.format(self.vcA), '{:.3f}'.format(self.vcB), '{:.3f}'.format(self.vcC), '{:.3f}'.format(self.shcA),
                       '{:.3f}'.format(self.shcB), '{:.3f}'.format(self.shcC), '{:.3f}'.format(self.mw), '{:.3f}'.format(self.shr)]

        return epentry("WindowMaterial:Gas", params, paramvs)

class No_En_Mat_Sh(Node, EnViMatNodes):
    '''Node defining an EnVi window shader'''
    bl_idname = 'No_En_Mat_Sh'
    bl_label = 'EnVi shade'

    st: FloatProperty(name="", description="Solar transmittance", min=0.0, max=1, default=0.05)
    sr: FloatProperty(name="", description="Solar reflectance", min=0.0, max=1, default=0.3)
    vt: FloatProperty(name="", description="Visible transmittance", min=0.0, max=1, default=0.05)
    vr: FloatProperty(name="", description="Visible reflectance", min=0.0, max=1, default=0.3)
    ihe: FloatProperty(name="", description="Infrared Hemispherical Emissivity", min=0.0, max=1, default=0.9)
    it: FloatProperty(name="", description="Infrared Transmittance", min=0.0, max=1, default=0.0)
    thi: FloatProperty(name="mm", description="Thickness", min=0.1, max=1000, default=5)
    tc: FloatProperty(name="W/m.K", description="Conductivity", min=0.0001, max=10, precision=3, default=0.1)
    sgd: FloatProperty(name="mm", description="Shade to glass distance", min=0.1, max=1000, default=50)
    tom: FloatProperty(name="", description="Top opening multiplier", min=0.0, max=1, default=0.5)
    bom: FloatProperty(name="", description="Bottom opening multiplier", min=0.0, max=1, default=0.5)
    lom: FloatProperty(name="", description="Left-side opening multiplier", min=0.0, max=1, default=0.5)
    rom: FloatProperty(name="", description="Right-side opening multiplier", min=0.0, max=1, default=0.5)
    afp: FloatProperty(name="", description="Air flow permeability", min=0.0, max=1, default=0.)
    envi_con_type: StringProperty(name="", description="Name")
    resist: FloatProperty(name = "", description = "", min = 0, default = 0)

    def init(self, context):
        self.outputs.new('So_En_Mat_Sh', 'Shade')
        self.outputs['Shade'].link_limit = 1
        self.inputs.new('So_En_Mat_Sh', 'Shade')
        self.inputs.new('So_En_Mat_ShC', 'Control')

    def draw_buttons(self, context, layout):
        if self.outputs["Shade"].links or self.inputs["Shade"].links:
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
        if self.outputs["Shade"].links or self.inputs["Shade"].links:
            nodecolour(self, 0)
        else:
            nodecolour(self, 1)

    def ret_resist(self):
        self.resist = 0
        return 0

    def ret_ec(self):
        return (0, 0)

    def update(self):
        for sock in self.outputs:
            socklink2(sock, self.id_data)

        if self.outputs['Shade'].links:
            if self.outputs['Shade'].links[0].to_node.bl_idname != 'No_En_Mat_Tr':
                self.id_data.links.remove(self.outputs['Shade'].links[0])

        if self.outputs["Shade"].links:
            self.inputs["Shade"].hide = True
        elif self.inputs["Shade"].links:
            self.outputs["Shade"].hide = True
        else:
            (self.inputs["Shade"].hide, self.outputs["Shade"].hide) = (False, False)

        self.valid()

    def ep_write(self, ln, mn):
        params = ('Name', 'Solar transmittance', 'Solar Reflectance', 'Visible transmittance', 'Visible reflectance', 'Infrared Hemispherical Emissivity', 'Infrared Transmittance', 'Thickness {m}',
                  'Conductivity {W/m-K}', 'Shade to glass distance {m}', 'Top opening multiplier', 'Bottom opening multiplier', 'Left-side opening multiplier',
                  'Right-side opening multiplier', 'Air flow permeability')
        paramvs = ['{}-shading-{}'.format(mn, ln)] + ['{:.3f}'.format(p) for p in (self.st, self.sr, self.vt, self.vr, self.ihe, self.it, 0.001 * self.thi, self.tc, 0.001 * self.sgd,
                   self.tom, self.bom, self.lom, self.rom, self.afp)]

        return epentry('WindowMaterial:Shade', params, paramvs) # + self.inputs['Control'].links[0].from_node.ep_write(ln, mn, zn, sn)

class No_En_Mat_Sc(Node, EnViMatNodes):
    '''Node defining an EnVi external screen'''
    bl_idname = 'No_En_Mat_Sc'
    bl_label = 'EnVi screen'

    rb: EnumProperty(items=[("DoNotModel", "DoNotModel", "Do not model reflected beam component"),
                               ("ModelAsDirectBeam", "ModelAsDirectBeam", "Model reflectred beam as beam"),
                               ("ModelAsDiffuse", "ModelAsDiffuse", "Model reflected beam as diffuse")],
                                name="", description="Composition of the layer", default="ModelAsDiffuse")
    ta: EnumProperty(items=[("0", "0", "Angle of Resolution for Screen Transmittance Output Map"),
                               ("1", "1", "Angle of Resolution for Screen Transmittance Output Map"),
                               ("2", "2", "Angle of Resolution for Screen Transmittance Output Map"),
                               ("3", "3", "Angle of Resolution for Screen Transmittance Output Map"),
                               ("5", "5", "Angle of Resolution for Screen Transmittance Output Map")],
                                name="", description="Angle of Resolution for Screen Transmittance Output Map", default="0")

    dsr: FloatProperty(name="", description="Diffuse solar reflectance", min=0.0, max=0.99, default=0.5)
    vr: FloatProperty(name="", description="Visible reflectance", min=0.0, max=1, default=0.6)
    the: FloatProperty(name="", description="Thermal Hemispherical Emissivity", min=0.0, max=1, default=0.9)
    tc: FloatProperty(name="W/m.K", description="Conductivity", min=0.0001, max=10, precision=3, default=0.1)
    sme: FloatProperty(name="mm", description="Screen Material Spacing", min=1, max=1000, default=50)
    smd: FloatProperty(name="mm", description="Screen Material Diameter", min=1, max=1000, default=25)
    sgd: FloatProperty(name="mm", description="Screen to glass distance", min=1, max=1000, default=25)
    tom: FloatProperty(name="", description="Top opening multiplier", min=0.0, max=1, default=0.0)
    bom: FloatProperty(name="", description="Bottom opening multiplier", min=0.0, max=1, default=0.0)
    lom: FloatProperty(name="", description="Left-side opening multiplier", min=0.0, max=1, default=0.0)
    rom: FloatProperty(name="", description="Right-side opening multiplier", min=0.0, max=1, default=0.0)
    resist: FloatProperty(name = "", description = "", min = 0, default = 0)

    def init(self, context):
        self.inputs.new('So_En_Mat_Sh', 'Shade')
        self.inputs.new('So_En_Mat_ShC', 'Control')

    def draw_buttons(self, context, layout):
        if self.inputs['Shade'].links:
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
        if not self.inputs["Shade"].links or not self.inputs["Control"].links:
            nodecolour(self, 1)
        else:
            nodecolour(self, 0)

    def update(self):
        for sock in self.outputs:
            socklink(sock, self.id_data.name)

        if self.outputs["Shade"].links:
            self.inputs["Shade"].hide = True
        elif self.inputs["Shade"].links:
            self.outputs["Shade"].hide = True
        else:
            (self.inputs["Shade"].hide, self.outputs["Shade"].hide) = (False, False)

        self.valid()

    def ret_resist(self):
        self.resist = 0
        return 0

    def ret_ec(self):
        return (0, 0)

    def ep_write(self, ln, mn):
        params = ('Name', 'Reflected Beam Transmittance Accounting Method', 'Diffuse Solar Reflectance', 'Diffuse Visible Reflectance',
                  'Thermal Hemispherical Emissivity', 'Conductivity (W/m-K)', 'Screen Material Spacing (m)', 'Screen Material Diameter (m)',
                  'Screen-to-Glass Distance (m)', 'Top Opening Multiplier', 'Bottom Opening Multiplier', 'Left-Side Opening Multiplier',
                  'Right-Side Opening Multiplier', 'Angle of Resolution for Output Map (deg)')

        paramvs = ['{}-shading-{}'.format(mn, ln), self.rb] + ['{:.3f}'.format(p) for p in (self.dsr, self.vr, self.the, self.tc, 0.001 * self.sme, 0.001 * self.smd,
                   0.001 * self.sgd, self.tom, self.bom, self.lom, self.rom)] + [self.ta]

        return epentry('WindowMaterial:Screen', params, paramvs)

class No_En_Mat_Bl(Node, EnViMatNodes):
    '''Node defining an EnVi window blind'''
    bl_idname = 'No_En_Mat_Bl'
    bl_label = 'EnVi blind'

    so: EnumProperty(items=[("0", "Horizontal", "Select from database"),
                                ("1", "Vertical", "Define custom material properties")],
                                name="", description="Slat orientation", default='0')
    sw: FloatProperty(name="mm", description="Slat width", min=0.1, max=1000, default=12)
    ss: FloatProperty(name="mm", description="Slat separation", min=0.1, max=1000, default=20)
    st: FloatProperty(name="mm", description="Slat thickness", min=0.1, max=1000, default=2)
    sa: FloatProperty(name="deg", description="Slat angle", min=0.0, max=90, default=45)
    stc: FloatProperty(name="W/m.K", description="Slat conductivity", min=0.01, max=100, default=10)
    sbst: FloatProperty(name="", description="Slat beam solar transmittance", min=0.0, max=1, default=0.0)
    fbsr: FloatProperty(name="", description="Front Side Slat beam solar reflectance", min=0.0, max=1, default=0.8)
    bbsr: FloatProperty(name="", description="Back Side Slat beam solar reflectance", min=0.0001, max=10, default=0.8)
    sdst: FloatProperty(name="", description="Slat diffuse solar transmittance", min=0.0, max=1, default=0.0)
    fdsr: FloatProperty(name="", description="Front Side Slat diffuse solar reflectance", min=0.0, max=1, default=0.8)
    bdsr: FloatProperty(name="", description="Back Side Slat diffuse solar reflectance", min=0.0, max=1, default=0.8)
    sbvt: FloatProperty(name="", description="Slat beam visible transmittance", min=0.0, max=1, default=0.0)
    fbvr: FloatProperty(name="", description="Front Side Slat beam visible reflectance", min=0.0, max=1, default=0.7)
    bbvr: FloatProperty(name="", description="Back Side Slat beam visible reflectance", min=0.0, max=1, default=0.7)
    sdvt: FloatProperty(name="", description="Slat diffuse visible transmittance", min=0.0, max=1, default=0.0)
    fdvr: FloatProperty(name="", description="Front Side Slat diffuse visible reflectance", min=0.0, max=1, default=0.7)
    bdvr: FloatProperty(name="", description="Back Side Slat diffuse visible reflectance", min=0.0, max=1, default=0.7)
    sit: FloatProperty(name="", description="Slat Infrared hemispherical transmittance", min=0.0, max=1, default=0.0)
    sfie: FloatProperty(name="", description="Front Side Slat Infrared hemispherical emissivity", min=0.0, max=1, default=0.9)
    sbie: FloatProperty(name="", description="Back Side Slat Infrared hemispherical emissivity", min=0.0, max=1, default=0.9)
    bgd: FloatProperty(name="mm", description="Blind-to-glass distance", min=0.1, max=1000, default=50)
    tom: FloatProperty(name="", description="Blind top opening multiplier", min=0.0, max=1, default=0.0)
    bom: FloatProperty(name="", description="Blind bottom opening multiplier", min=0.0, max=1, default=0.0)
    lom: FloatProperty(name="", description="Blind left-side opening multiplier", min=0.0, max=1, default=0.5)
    rom: FloatProperty(name="", description="Blind right-side opening multiplier", min=0.0, max=1, default=0.5)
    minsa: FloatProperty(name="deg", description="Minimum slat angle", min=0.0, max=90, default=0.0)
    maxsa: FloatProperty(name="deg", description="Maximum slat angle", min=0.0, max=90, default=90.0)
    resist: FloatProperty(name = "", description = "", min = 0, default = 0)
    thi: FloatProperty(name = "", description = "", min = 0, default = 0)

    def init(self, context):
        self.outputs.new('So_En_Mat_Sh', 'Shade')
        self.outputs['Shade'].link_limit = 1
        self.inputs.new('So_En_Mat_Sh', 'Shade')
        self.inputs.new('So_En_Mat_ShC', 'Control')

    def draw_buttons(self, context, layout):
        if self.outputs['Shade'].links or self.inputs['Shade'].links:
            newrow(layout, "Slat orient.:", self, "so")
            newrow(layout, "Slat width:", self, "sw")
            newrow(layout, "Slat sep.:", self, "ss")
            newrow(layout, "Slat thick.:", self, "st")
            newrow(layout, "Slat angle:", self, "sa")
            newrow(layout, "Slat cond.:", self, "stc")
            newrow(layout, "Slat beam trans.:", self, "sbst")
            newrow(layout, "Front beam reflec.:", self, "fbsr")
            newrow(layout, "Back beam reflec.:", self, "bbsr")
            newrow(layout, "Slat diff. reflec.:", self, "sdst")
            newrow(layout, "Front diff. reflec.:", self, "fdsr")
            newrow(layout, "Back diff. reflec.:", self, "bdsr")
            newrow(layout, "Slat beam vis. trans.:", self, "sbvt")
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
        if (self.outputs["Shade"].links or self.inputs["Shade"].links) and self.inputs["Control"].links:
            nodecolour(self, 0)
        else:
            nodecolour(self, 1)

    def update(self):
        for sock in self.outputs:
            socklink(sock, self.id_data.name)

        if self.outputs["Shade"].links:
            self.inputs["Shade"].hide = True
        elif self.inputs["Shade"].links:
            self.outputs["Shade"].hide = True
        else:
            (self.inputs["Shade"].hide, self.outputs["Shade"].hide) = (False, False)

        self.valid()

    def ret_resist(self):
        self.resist = 0
        return 0

    def ret_ec(self):
        return (0, 0)

    def ep_write(self, ln, mn):
        params = ('Name', 'Slat orientation', 'Slat width (m)', 'Slat separation (m)', 'Slat thickness (m)', 'Slat angle (deg)', 'Slat conductivity (W/m.K)',
                  'Slat beam solar transmittance', 'Front Side Slat beam solar reflectance', 'Back Side Slat beam solar reflectance', 'Slat diffuse solar transmittance',
                  'Front Side Slat diffuse solar reflectance', 'Back Side Slat diffuse solar reflectance', 'Slat beam visible transmittance', 'Front Side Slat beam visible reflectance',
                  'Back Side Slat beam visible reflectance', 'Slat diffuse visible transmittance', "Front Side Slat diffuse visible reflectance", "Back Side Slat diffuse visible reflectance",
                  "Slat Infrared hemispherical transmittance", "Front Side Slat Infrared hemispherical emissivity", "Back Side Slat Infrared hemispherical emissivity", "Blind-to-glass distance",
                  "Blind top opening multiplier", "Blind bottom opening multiplier", "Blind left-side opening multiplier", "Blind right-side opening multiplier", "Minimum slat angle", "Maximum slat angle")
        paramvs = ['{}-shading-{}'.format(mn, ln), ('Horizontal', 'Vertical')[int(self.so)]] + ['{:.3f}'.format(p) for p in (0.001 * self.sw, 0.001 * self.ss, 0.001 * self.st, self.sa, self.stc, self.sbst,
                                                                                                                             self.fbsr, self.bbsr, self.sdst, self.fdsr, self.bdsr, self.sbvt,
                                                                                                                             self.fbvr, self.bbvr, self.sdvt, self.fdvr, self.bdvr, self.sit, self.sfie,
                                                                                                                             self.sbie, 0.001 * self.bgd, self.tom, self.bom, self.lom, self.rom, self.minsa, self.maxsa)]

        return epentry('WindowMaterial:Blind', params, paramvs)

class No_En_Mat_SG(Node, EnViMatNodes):
    '''Node defining the EnVi switchable glazing layer'''
    bl_idname = 'No_En_Mat_SG'
    bl_label = 'EnVi switchable glazing layer'

    def lay_update(self, context):
        if not self.em.updated:
            self.em.update()

        if self.layer == '1' and self.lay_name == '':
            nodecolour(self, 1)
        elif self.layer == '0' and not self.material:
            nodecolour(self, 1)
        else:
            nodecolour(self, 0)

    def ec_update(self, context):
        self['ecentries'] = []

        if not self.ee.updated:
            self.ee.update()

        if self.embodiedtype in self.ee.propdict.keys():
            if self.embodiedclass not in self.ee.propdict[self.embodiedtype]:
                self.embodiedclass = list(self.ee.propdict[self.embodiedtype].keys())[0]

            if self.embodiedmat not in self.ee.propdict[self.embodiedtype][self.embodiedclass]:
                self.embodiedmat = list(self.ee.propdict[self.embodiedtype][self.embodiedclass])[0]

            try:
                self['ecdict'] = self.ee.propdict[self.embodiedtype][self.embodiedclass][self.embodiedmat]
                self['ecm2'] = '{:.3f}'.format(float(self['ecdict']['eckg']) * float(self['ecdict']['density']) * self.thi * 0.001)
                self['ecm2y'] = '{:.3f}'.format(float(self['ecdict']['eckg']) * float(self['ecdict']['density']) * self.thi * 0.001/self.ec_life)
                self['ecentries'] = [(k, self.ee.propdict[self.embodiedtype][self.embodiedclass][self.embodiedmat][k]) for k in self['ecdict'].keys()]

            except Exception as e:
                self['ecm2'] = 'N/A'
                self['ecm2y'] = 'N/A'

    layer: EnumProperty(items=[("0", "Database", "Select from database"),
                                        ("1", "Custom", "Define custom material properties")],
                                        name="", description="Composition of the layer", default="0")
    materialtype: EnumProperty(items=envi_layertype, name="", description="Layer material type")
    material: EnumProperty(items=envi_layer, name="", description="Glass material")
    thi: FloatProperty(name="mm", description="Thickness (mm)", min=0.1, max=10000, default=100)
    tc: FloatProperty(name="W/m.K", description="Thermal Conductivity (W/m.K)", min=0.1, max=10000, precision=3, default=0.9)
    stn: FloatProperty(name="", description="Solar normal transmittance", min=0, max=1, default=0.837)
    fsn: FloatProperty(name="", description="Solar front normal reflectance", min=0, max=1, default=0.075)
    bsn: FloatProperty(name="", description="Solar back normal reflectance", min=0, max=1, default=0.075)
    vtn: FloatProperty(name="", description="Visible Transmittance at Normal Incidence", min=0, max=1, default=0.898)
    fvrn: FloatProperty(name="", description="Front Side Visible Reflectance at Normal Incidence", min=0, max=1, default=0.081)
    bvrn: FloatProperty(name="", description="Back Side Visible Reflectance at Normal Incidence", min=0, max=1, default=0.081)
    itn: FloatProperty(name="", description="Infrared Transmittance at Normal Incidence", min=0, max=1, default=0.0)
    fie: FloatProperty(name="", description="Front Side Infrared Hemispherical Emissivity'", min=0, max=1, default=0.84)
    bie: FloatProperty(name="", description="Back Side Infrared Hemispherical Emissivity", min=0, max=1, default=0.84)
    diff: BoolProperty(name="", description="Diffusing", default=0)
    envi_con_type: StringProperty(name="", description="Name")
    resist: FloatProperty(name="", description="", min=0, default=0)
    embodied: BoolProperty(name="", description="Embodied carbon", default=0, update=ec_update)
    embodiedtype: EnumProperty(items=envi_elayertype, name="", description="Layer embodied material class", update=ec_update)
    embodiedclass: EnumProperty(items=envi_eclasstype, name="", description="Layer embodied class", update=ec_update)
    embodiedmat: EnumProperty(items=envi_emattype, name="", description="Layer embodied material", update=ec_update)
    ecm2: FloatProperty(name="kgCo2e/m2", description="Embodied carbon per metre squared", min=0.0, default=0.0)
    ec_id: StringProperty(name="", description="Embodied id")
    ec_type: StringProperty(name="", description="Embodied type")
    ec_class: StringProperty(name="", description="Embodied class")
    ec_name: StringProperty(name="", description="Embodied name")
    ec_unit:EnumProperty(items=[("kg", "kg", "per kilogram"),
                                  ("m2", "m2", "per square metre"),
                                  ("m3", "m3", "per cubic metre"),
                                  ("unit", "Unit", "per unit")],
                                  name="",
                                  description="Embodied carbon unit",
                                  default="kg")
    ec_amount: FloatProperty(name="", description="", min=0.001, precision=3, default=1)
    ec_kgco2e: FloatProperty(name="", description="Embodied carbon per kg amount", precision=3, default=100)

    ec_density: FloatProperty(name="kg/m^3", description="Material density", default=1000)
    ec_life: IntProperty(name="y", description="Service life in years", min=1, max=100, default=60, update=ec_update)
    ec_mod: StringProperty(name="", description="Embodied modules")
    em = envi_materials()
    ee = envi_embodied()

    def init(self, context):
        self.outputs.new('So_En_Mat_Sh', 'Shade')
        self.outputs['Shade'].link_limit = 1
        self.inputs.new('So_En_Mat_Sh', 'Shade')
        self.inputs['Shade'].hide = True
        self.inputs.new('So_En_Mat_ShC', 'Control')

    def draw_buttons(self, context, layout):
        if self.outputs['Shade'].links:
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
        for sock in self.outputs:
            socklink(sock, self.id_data.name)

        if self.outputs["Shade"].links:
            self.envi_con_type = self.outputs['Shade'].links[0].to_node.envi_con_type
            # self.inputs["Shade"].hide = True
        # elif self.inputs["Shade"].links:
        #     self.envi_con_type = self.inputs['Shade'].links[0].to_node.envi_con_type
        #     self.outputs["Shade"].hide = True
        # else:
        #     (self.inputs["Shade"].hide, self.outputs["Shade"].hide) = (False, False)

        self.valid()

    def valid(self):
        if not self.outputs["Shade"].links:
            nodecolour(self, 1)
        else:
            nodecolour(self, 0)

    def ret_resist(self):
        self.resist = 0
        return 0

    def ret_ec(self):
        return (0, 0)

    def ep_write(self, ln, mn):
        if not self.em.updated:
            self.em.update()

        layer_name = '{}-shading-{}'.format(mn, ln)
        params = ('Name', 'Optical Data Type', 'Window Glass Spectral Data Set Name', 'Thickness (m)', 'Solar Transmittance at Normal Incidence', 'Front Side Solar Reflectance at Normal Incidence',
                  'Back Side Solar Reflectance at Normal Incidence', 'Visible Transmittance at Normal Incidence', 'Front Side Visible Reflectance at Normal Incidence', 'Back Side Visible Reflectance at Normal Incidence',
                  'Infrared Transmittance at Normal Incidence', 'Front Side Infrared Hemispherical Emissivity', 'Back Side Infrared Hemispherical Emissivity', 'Conductivity (W/m-K)',
                  'Dirt Correction Factor for Solar and Visible Transmittance', 'Solar Diffusing')

        if self.layer == '0':
            matlist = list(self.em.matdat[self.material])
            paramvs = [layer_name] + matlist[1:3] + [self.thi] + ['{:.3f}'.format(float(sm)) for sm in matlist[4:-1]] + [1, ('No', 'Yes')[matlist[-1]]]

        else:
            paramvs = ['{}-shading-{}'.format(mn, ln), 'SpectralAverage', '', self.thi * 0.001, '{:.3f}'.format(self.stn), '{:.3f}'.format(self.fsn), '{:.3f}'.format(self.bsn),
                       '{:.3f}'.format(self.vtn), '{:.3f}'.format(self.fvrn), '{:.3f}'.format(self.bvrn), '{:.3f}'.format(self.itn),
                       '{:.3f}'.format(self.fie), '{:.3f}'.format(self.bie), '{:.3f}'.format(self.tc), 1, ('No', 'Yes')[self.diff]]

        return epentry("WindowMaterial:Glazing", params, paramvs)

class No_En_Mat_ShC(Node, EnViMatNodes):
    '''Node defining an EnVi window shade control'''
    bl_idname = 'No_En_Mat_ShC'
    bl_label = 'EnVi shade control'

    ttuple = ("AlwaysOn", "AlwaysOff", "OnIfScheduleAllows", "OnIfHighSolarOnWindow", "OnIfHighHorizontalSolar",
              "OnIfHighOutdoorAirTemperature",
              "OnIfHighZoneAirTemperature", "OnIfHighZoneCooling", "OnIfHighGlare", "MeetDaylightIlluminanceSetpoint",
              "OnNightIfLowOutdoorTempAndOffDay", "OnNightIfLowInsideTempAndOffDay", "OnNightIfHeatingAndOffDay",
              "OnNightIfLowOutdoorTempAndOnDayIfCooling", "OnNightIfHeatingAndOnDayIfCooling",
              "OffNightAndOnDayIfCoolingAndHighSolarOnWindow", "OnNightAndOnDayIfCoolingAndHighSolarOnWindow",
              "OnIfHighOutdoorAirTempAndHighSolarOnWindow", "OnIfHighOutdoorAirTempAndHighHorizontalSolar",
              "OnIfHighZoneAirTempAndHighSolarOnWindow", "OnIfHighZoneAirTempAndHighHorizontalSolar")

    def type_menu(self, context):
        try:
            if self.outputs['Control'].links[0].to_node.bl_idname == 'No_En_Mat_Sc':
                return [(self.ttuple[t], self.ttuple[t], self.ttuple[t]) for t in (0, 1, 2)]
            elif self.outputs['Control'].links[0].to_node.bl_idname in ('No_En_Mat_Bl', 'No_En_Mat_Sh', 'No_En_Mat_SG'):
                return [(self.ttuple[t], self.ttuple[t], self.ttuple[t]) for t in (0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20)]
            else:
                return [(t, t, t) for t in self.ttuple]
        except Exception:
            return [('None', 'None', 'None')]

    def schupdate(self, context):
        if self.ctype not in ("AlwaysOn", "AlwaysOff", "OnIfHighGlare"):
            self.inputs['Schedule'].hide = False
        else:
            if self.inputs['Schedule'].links:
                self.id_data.links.remove(self.inputs['Schedule'].links[0])
            self.inputs['Schedule'].hide = True

        if self.ctype == 'OnIfScheduleAllows':
            if not self.inputs['Schedule'].links:
                nodecolour(self, 1)

    def slatupdate(self, context):
        if self.sac != 'ScheduledSlatAngle':
            if self.inputs['Slat schedule'].links:
                self.id_data.links.remove(self.inputs['Slat schedule'].links[0])
            self.inputs['Slat schedule'].hide = True
        else:
            self.inputs['Slat schedule'].hide = False

            if not self.inputs['Slat schedule'].links:
                nodecolour(self, 1)


    ctype: EnumProperty(items=type_menu, name="", description="Shading device", update=schupdate)
    sp: FloatProperty(name="", description="Setpoint (W/m2, W or deg C)", min=0.0, max=1000, default=20)
    sac: EnumProperty(items=[("FixedSlatAngle", "FixedSlatAngle", "Constant slat angle"),
                                ("ScheduledSlatAngle", "ScheduledSlatAngle", "Scheduled slat angle"),
                                ("BlockBeamSolar", "BlockBeamSolar", "Block beam solar")
                                ],
                                name="", description="Blind slat angle control", default='FixedSlatAngle', update=slatupdate)
    sp2: FloatProperty(name="", description="Setpoint 2 (W/m2, W or deg C)", min=0.0, max=1000, default=20)
    resist: FloatProperty(name="", description="", min=0, default=0)

    def init(self, context):
        self.outputs.new('So_En_Mat_ShC', 'Control')
        self.outputs['Control'].link_limit = 1
        self.inputs.new('So_En_Mat_Sched', 'Schedule')
        self.inputs.new('So_En_Mat_Sched', 'Slat schedule')
        self.inputs['Schedule'].hide = True
        self.inputs['Slat schedule'].hide = True

    def draw_buttons(self, context, layout):
        if self.outputs.get('Control'):
            newrow(layout, "Shading device:", self, 'ctype')

            if self.ctype not in ('AlwaysOn', 'AlwaysOff', 'OnIfScheduleAllows', 'OnIfHighGlare', 'DaylightIlluminance'):
                newrow(layout, "Set-point", self, 'sp')
            if self.outputs['Control'].links and self.outputs['Control'].links[0].to_node.bl_idname == 'No_En_Mat_Bl':
                newrow(layout, 'Slat angle:', self, 'sac')
            if self.ctype in ("OnIfHighOutdoorAirTempAndHighSolarOnWindow", "OnIfHighOutdoorAirTempAndHighHorizontalSolar", "OnIfHighZoneAirTempAndHighSolarOnWindow", "OnIfHighZoneAirTempAndHighHorizontalSolar"):
                newrow(layout, "Set-point 2", self, 'sp2')

    def valid(self):
        if not self.outputs["Control"].links or (self.ctype == "OnIfScheduleAllows" and not self.inputs['Schedule'].links):
            nodecolour(self, 1)
        else:
            nodecolour(self, 0)

    def update(self):
        for sock in self.outputs:
            socklink(sock, self.id_data.name)

        if not self.outputs['Control'].links:
            self.ctype = 'None'
        else:
            if self.outputs['Control'].links[0].to_node.bl_idname != 'No_En_Mat_Bl':
                if self.inputs['Slat schedule'].links:
                    self.id_data.links.remove(self.inputs['Slat schedule'].links[0])
                self.inputs['Slat schedule'].hide = True
            else:
                self.inputs['Slat schedule'].hide = False

        self.valid()
        self.ctype = self.ctype

    def ret_resist(self):
        self.resist = 0
        return 0

    def ret_ec(self):
        return (0, 0)

    def ep_write(self, ln, mn, zn, sn):
        shade_node = self.outputs['Control'].links[0].to_node

        if shade_node.bl_idname == 'No_En_Mat_Sc':
            st = 'ExteriorScreen'
        elif shade_node.bl_idname == 'No_En_Mat_Bl':
            if shade_node.inputs['Shade'].links and shade_node.inputs['Shade'].links[0].from_node.outputs['Layer'].links[0].to_node.bl_idname == 'No_En_Mat_Con':
                st = 'ExteriorBlind'
            elif (shade_node.outputs['Shade'].links and shade_node.outputs['Shade'].links[0].to_node.inputs['Layer'].links) or (shade_node.inputs['Shade'].links and shade_node.inputs['Shade'].links[0].from_node.outputs['Layer'].links):
                st = 'BetweenGlassBlind'
            else:
                st = 'InteriorBlind'
        elif shade_node.bl_idname == 'No_En_Mat_Sh':
            if shade_node.inputs['Shade'].links and shade_node.inputs['Shade'].links[0].from_node.outputs['Layer'].links[0].to_node.bl_idname == 'No_En_Mat_Con':
                st = 'ExteriorShade'
            elif (shade_node.outputs['Shade'].links and shade_node.outputs['Shade'].links[0].to_node.inputs['Layer'].links) or (shade_node.inputs['Shade'].links and shade_node.inputs['Shade'].links[0].from_node.outputs['Layer'].links):
                st = 'BetweenGlassShade'
            else:
                st = 'InteriorShade'
        else:
            st = 'SwitchableGlazing'

        (scs, scn) = ('Yes', '{}-shading-schedule'.format(sn)) if self.inputs['Schedule'].links else ('No', '')
        slcn = '{}-slat-schedule'.format(sn) if self.inputs['Slat schedule'].links else ''

        params = ('Name', 'Zone Name', 'Shading Control Sequence Number', 'Shading Type', 'Construction with Shading Name', 'Shading Control Type', 'Schedule Name', 'Setpoint (W/m2, W or deg C)', 'Shading Control Is Scheduled',
                  'Glare Control Is Active', 'Shading Device Material Name', 'Type of Slat Angle Control for Blinds', 'Slat Angle Schedule Name', 'Setpoint 2 (W/m2, deg C or cd/m2)', 'Daylighting Control Object Name',
                  'Multiple Surface Control Type', 'Fenestration Surface 1 Name')
        paramvs = ('{}-shading-control'.format(sn), zn, 1, st, '{}-shading'.format(mn), self.ctype, scn, self.sp, scs, 'No', '', self.sac, slcn, '', '', 'Sequential', sn)

        ss_text = self.inputs['Slat schedule'].links[0].from_node.ep_write(slcn, 'Any number') if self.inputs['Slat schedule'].links else ''
        return epentry('WindowShadingControl', params, paramvs) + ss_text

class No_En_Mat_PV(Node, EnViMatNodes):
    '''Node defining an EnVi photovoltaic module'''
    bl_idname = 'No_En_Mat_PV'
    bl_label = 'EnVi PV'

    def pv_update(self, context):
        pass

    sandia_dict = {}
    l = -40

    def ret_e1dmenu(self, context):
        with open(ret_datab('PV_database.json', 'r'), 'r') as e1d_jfile:
            e1d_dict = json.loads(e1d_jfile.read())
        return [(p, p, '{} module'.format(p)) for p in e1d_dict]

    def ret_sandiamenu(self, context):
        # with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'EPFiles', '{}'.format('SandiaPVdata.json')), 'r') as sandia_json:
        #     sandiadict = json.loads(sandia_json.read())
        return [(p, p, '{} module'.format(p)) for p in self.sandiadict]

    # def ret_e1ddict(self):
    #     with open(ret_datab('PV_database.json', 'r'), 'r') as e1d_jfile:
    #         e1d_dict = json.loads(e1d_jfile.read())
    #         return e1d_dict

    def save_e1ddict(self):
        with open(ret_datab('PV_database.json', 'r'), 'r') as e1d_jfile:
            e1d_dict = json.loads(e1d_jfile.read())

        e1d_dict[self.pv_name] = [f'{self.scc:.2f}', f'{self.ocv:.2f}', f'{self.mv:.2f}', f'{self.mc:.2f}', f'{self.tcscc:.3f}', f'{self.tcocv:.3f}', f'{self.cis:.2f}', f'{self.ctnoct:.2f}', f'{self.mod_area:.2f}']
        #self.e1d_dict = e1d_dict

        with open(ret_datab('PV_database.json', 'w'), 'w') as e1d_jfile:
            e1d_jfile.write(json.dumps(e1d_dict, indent=2))

    # with open(ret_datab('PV_database.json', 'r'), 'r') as e1d_jfile:
    #     e1d_dict = json.loads(e1d_jfile.read())

    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'EPFiles', '{}'.format('SandiaPVdata.json')), 'r') as sandia_json:
        sandiadict = json.loads(sandia_json.read())

    ct: EnumProperty(items=[("0", "Crystalline", "Do not model reflected beam component"),
                               ("1", "Amorphous", "Model reflectred beam as beam")],
                                name="", description="Photovoltaic Type", default="0")

    hti: EnumProperty(items=[("Decoupled", "Decoupled", "Decoupled"),
                               ("DecoupledUllebergDynamic", "Ulleberg", "DecoupledUllebergDynamic")],
                                name="", description="Conversion Efficiency Input Mode'", default="Decoupled")

    pp: EnumProperty(items=[("0", "Simple", "Do not model reflected beam component"),
                               ("1", "One-Diode", "Model reflectred beam as beam"),
                               ("2", "Sandia", "Model reflected beam as diffuse")],
                                name="", description="Photovoltaic Performance Object Type", default="0", update=pv_update)
    pv_name: StringProperty(name='', description='Name of the custom PV model', default='')
    e1dmenu: EnumProperty(items=ret_e1dmenu, name="", description="Module type")
    smenu: EnumProperty(items=ret_sandiamenu, name="", description="Module type")
    mod_area: FloatProperty(name="m2", description="PV module area", min=0.1, default=5)
    pvsa: FloatProperty(name="%", description="Fraction of Surface Area with Active Solar Cells", min=50, max=100, default=90)
    aa: FloatProperty(name="m2", description="Active area", min=0.1, max=10000, default=5)
    eff: FloatProperty(name="%", description="Visible reflectance", min=0.0, max=100, default=20)
    ssp: IntProperty(name="", description="Number of series strings in parallel", min=1, max=100, default=1)
    mis: IntProperty(name="", description="Number of modules in series", min=1, max=100, default=1)
    cis: IntProperty(name="", description="Number of cells in series", min=1, max=100, default=36)
    tap: FloatProperty(name="", description="Transmittance absorptance product", min=-1, max=1, precision=3, default=0.9)
    sbg: FloatProperty(name="eV", description="Semiconductor band-gap", min=0.1, max=5, default=1.12)
    sr: FloatProperty(name="W", description="Shunt resistance", min=1, default=1000000)
    scc: FloatProperty(name="Amps", description="Short circuit current", min=1, max=1000, default=25)
    sgd: FloatProperty(name="mm", description="Screen to glass distance", min=1, max=1000, default=25)
    ocv: FloatProperty(name="V", description="Open circuit voltage", min=0.0, max=100, default=60)
    rt: FloatProperty(name="C", description="Reference temperature", min=0, max=40, default=25)
    ri: FloatProperty(name="W/m2", description="Reference insolation", min=100, max=2000, default=1000)
    mc: FloatProperty(name="Amps", description="Module current at maximum power", min=1, max=10, precision=3, default=5.6)
    mv: FloatProperty(name="V", description="Module voltage at maximum power", min=0.0, max=75, precision=3, default=17)
    tcscc: FloatProperty(name="A/K", description="Temperature Coefficient of Short Circuit Current", precision=5, min=0.00001, max=0.01, default=0.002)
    tcocv: FloatProperty(name="V/K", description="Temperature Coefficient of Open Circuit Voltage", precision=5, min=-0.5, max=0, default=-0.1)
    atnoct: FloatProperty(name="C", description="Reference ambient temperature", min=0, max=40, default=20)
    ctnoct: FloatProperty(name="C", description="Nominal Operating Cell Temperature Test Cell Temperature", min=0, max=60, default=45)
    inoct: FloatProperty(name="W/m2", description="Nominal Operating Cell Temperature Test Insolation", min=100, max=2000, default=800)
    hlc: FloatProperty(name="W/m2.K", description="Module heat loss coefficient", min=0.0, max=50, default=30)
    thc: FloatProperty(name=" J/m2.K", description=" Total Heat Capacity", min=10000, max=100000, default=50000)

    def init(self, context):
        mat = bpy.data.materials[self.id_data.name]

        if not mat.vi_params.get('enparams'):
            mat.vi_params['enparams'] = {'area' : -1}

        self['area'] = 0
        self.inputs.new('So_En_Mat_PVG', 'PV Generator')
        self.inputs.new('So_En_Mat_Sched', 'PV Schedule')
        self.outputs.new('So_En_Mat_PV', 'PV')

    def draw_buttons(self, context, layout):
        mat = bpy.data.materials[self.id_data.name]

        if mat.vi_params['enparams'].get('pvarea'):
            row = layout.row()
            # row.operator('node.pv_area', text = "Area Calc")

            # try:
            row = layout.row()
            row.label(text = 'Area = {:.2f}m2'.format(mat.vi_params['enparams'].get('pvarea')))

        else:
            row = layout.row()
            row.label(text = 'Area = N/A')

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
        for sock in self.outputs:
            socklink(sock, self.id_data.name)

        if len(self.outputs) + len(self.inputs) == 3:
            if self.outputs['PV'].links and not self.inputs['PV Generator'].links:
                nodecolour(self, 1)
            else:
                nodecolour(self, 0)

    def ep_write(self, sn, area):
        with open(ret_datab('PV_database.json', 'r'), 'r') as e1d_jfile:
            e1d_dict = json.loads(e1d_jfile.read())

        self['matname'] = get_mat(self, 1).name
        marea = float((area, (self.mod_area, e1d_dict[self.e1dmenu][8])[self.e1dmenu != 'Custom'], self.sandiadict[self.smenu][0])[int(self.pp)])

        if marea/area > 1.01:
            logentry("No PV data exported as the face {} area is smaller than the module area".format(sn))
            return '! No PV data exported as the face {} area is smaller than the module area\n\n'.format(sn)

        params = ('Name', 'Surface Name', 'Photovoltaic Performance Object Type',
                  'Module Performance Name', 'Heat Transfer Integration Mode',
                  'Number of Series Strings in Parallel', 'Number of Modules in Series')

        paramvs = ['{}-pv'.format(sn), sn,
                   ('PhotovoltaicPerformance:Simple', 'PhotovoltaicPerformance:EquivalentOne-Diode', 'PhotovoltaicPerformance:Sandia')[int(self.pp)], '{}-pv-performance'.format(sn),
                   self.hti, self.ssp, int(area/marea)]

        ep_text = epentry('Generator:Photovoltaic', params, paramvs)

        if self.pp == '0':
            params = ('Name', 'Fraction of Surface Area with Active Solar Cell', 'Conversion Efficiency Input Mode', 'Value for Cell Efficiency if Fixed', 'Efficiency Schedule Name')
            paramvs = ('{}-pv-performance'.format(sn), self.pvsa * 0.01, ('Fixed', 'Scheduled')[len(self.inputs['PV Schedule'].links)], self.eff * 0.01, ('', '{}-pv-performance-schedule'.format(sn))[len(self.inputs['PV Schedule'].links)])
            ep_text += epentry('PhotovoltaicPerformance:Simple', params, paramvs)

            if self.inputs['PV Schedule'].links:
                ep_text += self.inputs['PV Schedule'].links[0].from_node.ep_write('{}-pv-performance-schedule'.format(sn), 'Fraction')

        elif self.pp == '1':
            params = ('Name', 'Cell type', 'Number of Cells in Series', 'Active Area (m2)', 'Transmittance Absorptance Product',
                      'Semiconductor Bandgap (eV)', 'Shunt Resistance (ohms)', 'Short Circuit Current (A)', 'Open Circuit Voltage (V)',
                      'Reference Temperature (C)', 'Reference Insolation (W/m2)', 'Module Current at Maximum Power (A)',
                      'Module Voltage at Maximum Power (V)', 'Temperature Coefficient of Short Circuit Current (A/K)',
                      'Temperature Coefficient of Open Circuit Voltage (V/K)', 'Nominal Operating Cell Temperature Test Ambient Temperature (C)',
                      'Nominal Operating Cell Temperature Test Cell Temperature (C)', 'Nominal Operating Cell Temperature Test Insolation (W/m2)',
                      'Module Heat Loss Coefficient (W/m2-K)', 'Total Heat Capacity (J/m2-K)')
            paramvs = ('{}-pv-performance'.format(sn), ('CrystallineSilicon', 'AmorphousSilicon')[int(self.ct)], (self.cis, e1d_dict[self.e1dmenu][6])[self.e1dmenu != 'Custom'], (self.mod_area, e1d_dict[self.e1dmenu][8])[self.e1dmenu != 'Custom'],
                       f'{self.tap:.3f}', f'{self.sbg:.3f}', self.sr, (f'{self.scc:.3f}', e1d_dict[self.e1dmenu][0])[self.e1dmenu != 'Custom'], (f'{self.ocv:.3f}', e1d_dict[self.e1dmenu][1])[self.e1dmenu != 'Custom'],
                       self.rt, self.ri, (f'{self.mc:.3f}', e1d_dict[self.e1dmenu][3])[self.e1dmenu != 'Custom'], (f'{self.mv:.3f}', e1d_dict[self.e1dmenu][2])[self.e1dmenu != 'Custom'],
                       (f'{self.tcscc:.3f}', e1d_dict[self.e1dmenu][4])[self.e1dmenu != 'Custom'], (f'{self.tcocv:.3f}', e1d_dict[self.e1dmenu][5])[self.e1dmenu != 'Custom'],
                       self.atnoct, (self.ctnoct, e1d_dict[self.e1dmenu][7])[self.e1dmenu != 'Custom'], self.inoct, self.hlc, self.thc)
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

    def init(self, context):
        self.inputs.new('So_En_Mat_Sched', 'Schedule')
        self.outputs.new('So_En_Mat_PVG', 'PV Gen')

    def draw_buttons(self, context, layout):
        newrow(layout, "Inverter type:", self, "it")
        newrow(layout, "Inverter efficiency:", self, "ie")
        newrow(layout, 'Radiative fraction:', self, 'rf')

class No_En_Mat_Sched(Node, EnViMatNodes):
    '''Node describing a schedule'''
    bl_idname = 'No_En_Mat_Sched'
    bl_label = 'Schedule'
    bl_icon = 'TIME'

    def tupdate(self, context):
        try:
            if self.source == '1':
                if os.path.isfile(bpy.path.abspath(self.select_file)):
                    nodecolour(self, 0)
                else:
                    nodecolour(self, 1)
            else:
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

        except Exception:
            nodecolour(self, 1)

    source: EnumProperty(name = '', items = [("0", "None", "No file"), ("1", "Select", "Select file")], default = '0')
    select_file: StringProperty(name="", description="Name of the variable file", default="", subtype="FILE_PATH")
    cn: IntProperty(name = "", default = 1, min = 1)
    rtsat: IntProperty(name = "", default = 0, min = 0)
    hours: IntProperty(name = "", default = 8760, min = 1, max = 8760)
    delim: EnumProperty(name = '', items = [("Comma", "Comma", "Comma delimiter"), ("Space", "Space", "space delimiter")], default = 'Comma')
    generate_file: StringProperty(default = "", name = "")
    u1: StringProperty(name = "", description = "Valid entries (; separated for each 'For', comma separated for each day, space separated for each time value pair)", update = tupdate)
    u2: StringProperty(name = "", description = "Valid entries (; separated for each 'For', comma separated for each day, space separated for each time value pair)", update = tupdate)
    u3: StringProperty(name = "", description = "Valid entries (; separated for each 'For', comma separated for each day, space separated for each time value pair)", update = tupdate)
    u4: StringProperty(name = "", description = "Valid entries (; separated for each 'For', comma separated for each day, space separated for each time value pair)", update = tupdate)
    f1: StringProperty(name = "", description = "Valid entries (space separated): AllDays, Weekdays, Weekends, Monday, Tuesday, Wednesday, Thursday, Friday, Saturday, Sunday, AllOtherDays", update = tupdate)
    f2: StringProperty(name = "", description = "Valid entries (space separated): AllDays, Weekdays, Weekends, Monday, Tuesday, Wednesday, Thursday, Friday, Saturday, Sunday, AllOtherDays", update = tupdate)
    f3: StringProperty(name = "", description = "Valid entries (space separated): AllDays, Weekdays, Weekends, Monday, Tuesday, Wednesday, Thursday, Friday, Saturday, Sunday, AllOtherDays", update = tupdate)
    f4: StringProperty(name = "", description = "Valid entries (space separated): AllDays, Weekdays, Weekends, Monday, Tuesday, Wednesday, Thursday, Friday, Saturday, Sunday, AllOtherDays", update = tupdate)
    t1: IntProperty(name = "", default = 365, min = 1, max = 365, update = tupdate)
    t2: IntProperty(name = "", default = 365, min = 1, max = 365, update = tupdate)
    t3: IntProperty(name = "", default = 365, min = 1, max = 365, update = tupdate)
    t4: IntProperty(name = "", default = 365, min = 1, max = 365, update = tupdate)

    def init(self, context):
        self.outputs.new('So_En_Net_Sched', 'Schedule')
        self['scheddict'] = {'TSPSchedule': 'Any Number', 'VASchedule': 'Fraction', 'Fan Schedule': 'Fraction', 'HSchedule': 'Temperature', 'CSchedule': 'Temperature'}
        self.tupdate(context)
        nodecolour(self, 1)

    def draw_buttons(self, context, layout):
        uvals, u = (1, self.u1, self.u2, self.u3, self.u4), 0
        tvals = (0, self.t1, self.t2, self.t3, self.t4)
        newrow(layout, 'Source:', self, 'source')

        if self.source == "1":
            newrow(layout, 'Select:', self, 'select_file')
            newrow(layout, 'Columns:', self, 'cn')
            newrow(layout, 'Skip rows:', self, 'rtsat')
            newrow(layout, 'Delimiter:', self, 'delim')
#        elif self.source == "2":
#            newrow(layout, 'Generate', self, 'generate_file')

        if self.source != "1":
            while uvals[u] and tvals[u] < 365:
                [newrow(layout, v[0], self, v[1]) for v in (('End day {}:'.format(u+1), 't'+str(u+1)), ('Fors:', 'f'+str(u+1)), ('Untils:', 'u'+str(u+1)))]
                u += 1

    def update(self):
        for sock in self.outputs:
            socklink(sock, self.id_data.name)

        self.id_data.interface_update(bpy.context)

    def ep_write(self, name, stype):
        if self.source == '0':
            schedtext, ths = '', []

            for tosock in [link.to_socket for link in self.outputs['Schedule'].links]:
                if not schedtext:
                    for t in (self.t1, self.t2, self.t3, self.t4):
                        ths.append(t)
                        if t == 365:
                            break

                    fos = [fs for fs in (self.f1, self.f2, self.f3, self.f4) if fs]
                    uns = [us for us in (self.u1, self.u2, self.u3, self.u4) if us]
                    ts, fs, us = rettimes(ths, fos, uns)
                    schedtext = epschedwrite(name, stype, ts, fs, us)
            # return schedtext
        else:
            params = ('Name', 'ScheduleType', 'Name of File', 'Column Number', 'Rows to Skip at Top', 'Number of Hours of Data', 'Column Separator')
            paramvs = (name, 'Any number', bpy.path.abspath(self.select_file), self.cn, self.rtsat, 8760, self.delim)
            schedtext = epentry('Schedule:File', params, paramvs)
            # return schedtext

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

                fos = [fs for fs in (self.f1, self.f2, self.f3, self.f4) if fs]
                uns = [us for us in (self.u1, self.u2, self.u3, self.u4) if us]
                ts, fs, us = rettimes(ths, fos, uns)
        for t in ts:
            for f in fs:
                for u in us:
                    for hi, h in enumerate((datetime.datetime(2015, 1, 1, 0, 00) - datetime.datetime(2014, 1, 1, 0, 00)).hours):
                        if h.day <= self.ts:
                            data[hi] = 1

        with open(os.path.join(newdir, name), 'w') as sched_file:
            sched_file.write(',\n'.join([d for d in data]))

        params = ('Name', 'ScheduleType', 'Name of File', 'Column Number', 'Rows to Skip at Top', 'Number of Hours of Data', 'Column Separator')
        paramvs = (name, 'Any number', os.path.abspath(self.select_file), self.cn, self.rtsat, 8760, self.delim)
        schedtext = epentry('Schedule:File', params, paramvs)
        return schedtext

class No_En_Mat_Anim(Node, EnViNodes):
    '''Node to automate changes in parameters'''
    bl_idname = 'No_En_Mat_Anim'
    bl_label = 'VI Animation'
    bl_icon = 'ANIM'

    def retparams(self, context):
        if self.inputs[0].links:
            return [(p.identifier, p.description, p.identifier) for p in self.inputs[0].links[0].from_node.bl_rna.properties if p.is_skip_save]
        else:
            return [('None', 'None', 'None')]

    parameter: EnumProperty(name='', description = 'Parameter to be animated', items=retparams)
    anim_file: StringProperty(name = '')

    def init(self, context):
        self.inputs.new('So_Anim', 'Parameter')

    def draw_buttons(self, context, layout):
        newrow(layout, "Parameter:", self, 'parameter')
        layout.prop_search(self, 'anim_file', bpy.data, 'texts', text='File', icon='TEXT')

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
envi_mat_para = [NodeItem("No_En_Mat_Anim", label="Animation")]

envimatnode_categories = [
        EnViMatNodeCategory("Type", "Type Node", items=envi_mat_con),
        EnViMatNodeCategory("Layer", "Layer Node", items=envi_mat_lay),
        EnViMatNodeCategory("Shading", "Shading Node", items=envi_mat_sha),
        EnViMatNodeCategory("Schedule", "Schedule Node", items=envi_mat_sch),
        EnViMatNodeCategory("Power", "PV Node", items=envi_mat_pv),
        EnViMatNodeCategory("Parametric", "Parametric Node", items=envi_mat_para)]



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

#        elif self.contextmenu == "Compliance":
#            if self.canalysismenu in ('0', '1', '2'):
#                self['skytypeparams'] = ("-b 22.86 -c", "-b 22.86 -c", "-b 18 -u")[int(self.canalysismenu)]
#                skyentry = livi_sun(scene, self, 0) + livi_sky(3)
#
#                if self.canalysismenu in ('0', '1'):
#                    self.starttime = datetime.datetime(2015, 1, 1, 12)
#                    self['preview'] = 1
#                    if self.hdr:
#                        hdrexport(scene, 0, scene.frame_current, self, skyentry)
#                else:
#                    self.starttime = datetime.datetime(2015, 9, 11, 9)
#                self['Text'][str(scene.frame_current)] = skyentry
#            else:
#                if self.sourcemenu2 == '0':
#                    self['mtxfile'] = cbdmmtx(self, scene, self.inputs['Location in'].links[0].from_node, export_op)
#                elif self.sourcemenu2 == '1':
#                    self['mtxfile'] = self.mtxname
#
#                self['Text'][str(scene.frame_current)] = "void glow sky_glow \n0 \n0 \n4 1 1 1 0 \nsky_glow source sky \n0 \n0 \n4 0 0 1 180 \nvoid glow ground_glow \n0 \n0 \n4 1 1 1 0 \nground_glow source ground \n0 \n0 \n4 0 0 -1 180\n\n"
#
#                if self.sourcemenu2 == '0':
#                    with open("{}.mtx".format(os.path.join(svp['viparams']['newdir'], self['epwbase'][0])), 'r') as mtxfile:
#                        self['Options']['MTX'] = mtxfile.read()
#                else:
#                    with open(self.mtxname, 'r') as mtxfile:
#                        self['Options']['MTX'] = mtxfile.read()
#                if self.hdr:
#                    self['Text'][str(scene.frame_current)] = cbdmhdr(self, scene)

       # params = ''

        # if self.buoyancy:
        #     params += 't'

        #     if self.buossinesq:
        #         params += 'b'
        #     if self.radiation:
        #         if self.radmodel == '0':
        #             params += 'p'
        #         else:
        #             params += 'f'

        # if self.turbulence == 'laminar':
        #     params += 'l'
        # elif self.turbulence == 'kEpsilon':
        #     params += 'k'
        # elif self.turbulence == 'kOmega':
        #     params += 'o'
        # elif self.turbulence == 'SpalartAllmaras':
        #     params += 's'

        # if context.scene.vi_params.get('flparams') and context.scene.vi_params['flparams'].get('solver_type'):
        #     context.scene.vi_params['flparams']['params'] = params

        # class No_En_Net_TC(Node, EnViNodes):
#     '''Zone Thermal Chimney node'''
#     bl_idname = 'No_En_Net_TC'
#     bl_label = 'Chimney'
#     bl_icon = 'SOUND'

#     def zupdate(self, context):
#         zonenames= []
#         obj = bpy.data.objects[self.zone]
#         odm = obj.data.materials
#         bsocklist = ['{}_{}_b'.format(odm[face.material_index].name, face.index) for face in obj.data.polygons if get_con_node(odm[face.material_index].vi_params).envi_con_con  == 'Zone' and odm[face.material_index].name not in [outp.name for outp in self.outputs if outp.bl_idname == 'So_En_Net_Bound']]

#         for oname in [outputs for outputs in self.outputs if outputs.name not in bsocklist and outputs.bl_idname == 'So_En_Net_Bound']:
#             self.outputs.remove(oname)

#         for iname in [inputs for inputs in self.inputs if inputs.name not in bsocklist and inputs.bl_idname == 'So_En_Net_Bound']:
#             self.inputs.remove(iname)

#         for sock in sorted(set(bsocklist)):
#             if not self.outputs.get(sock):
#                 self.outputs.new('So_En_Net_Bound', sock).sn = sock.split('_')[-2]
#             if not self.inputs.get(sock):
#                 self.inputs.new('So_En_Net_Bound', sock).sn = sock.split('_')[-2]

#         for sock in (self.inputs[:] + self.outputs[:]):
#             if sock.bl_idname == 'So_En_Net_Bound' and sock.links:
#                 zonenames += [(link.from_node.zone, link.to_node.zone)[sock.is_output] for link in sock.links]

#         nodecolour(self, all([get_con_node(mat.vi_params).envi_con_type != 'Window' for mat in bpy.data.objects[self.zone].data.materials if mat and mat.vi_params.envi_nodes]))
#         self['zonenames'] = zonenames

#     def supdate(self, context):
#         self.inputs.new['Schedule'].hide = False if self.sched == 'Sched' else True

#     zone: StringProperty(name = '', default = "en_Chimney")
#     sched: EnumProperty(name="", description="Ventilation control type", items=[('On', 'On', 'Always on'), ('Off', 'Off', 'Always off'), ('Sched', 'Schedule', 'Scheduled operation')], default='On', update = supdate)
#     waw: FloatProperty(name = '', min = 0.001, default = 1)
#     ocs: FloatProperty(name = '', min = 0.001, default = 1)
#     odc: FloatProperty(name = '', min = 0.001, default = 0.6)

#     def init(self, context):
#         self.inputs.new('So_En_Net_Sched', 'Schedule')
#         self['zonenames'] = []

#     def draw_buttons(self, context, layout):
#         newrow(layout, 'Zone:', self, 'zone')
#         newrow(layout, 'Schedule:', self, 'sched')
#         newrow(layout, 'Width Absorber:', self, 'waw')
#         newrow(layout, 'Outlet area:', self, 'ocs')
#         newrow(layout, 'Outlet DC:', self, 'odc')

#         for z, zn in enumerate(self['zonenames']):
#             row=layout.row()
#             row.label(zn)
#             row=layout.row()
#             row.prop(self, '["Distance {}"]'.format(z))
#             row=layout.row()
#             row.prop(self, '["Relative Ratio {}"]'.format(z))
#             row=layout.row()
#             row.prop(self, '["Cross Section {}"]'.format(z))

#     def update(self):
#         bi, bo = 1, 1
#         zonenames, fheights, fareas = [], [], []

#         for inp in [inp for inp in self.inputs if inp.bl_idname == 'So_En_Net_Bound']:
#             self.outputs[inp.name].hide = True if inp.is_linked and self.outputs[inp.name].bl_idname == inp.bl_idname else False

#         for outp in [outp for outp in self.outputs if outp.bl_idname in 'So_En_Net_Bound']:
#             self.inputs[outp.name].hide = True if outp.is_linked and self.inputs[outp.name].bl_idname == outp.bl_idname else False

#         if [inp for inp in self.inputs if inp.bl_idname == 'So_En_Net_Bound' and not inp.hide and not inp.links]:
#             bi = 0

#         if [outp for outp in self.outputs if outp.bl_idname == 'So_En_Net_Bound' and not outp.hide and not outp.links]:
#             bo = 0

#         nodecolour(self, not all((bi, bo)))

#         for sock in [sock for sock in self.inputs[:] + self.outputs[:] if sock.bl_idname == 'So_En_Net_Bound']:
#             if sock.links and self.zone in [o.name for o in bpy.data.objects]:
#                 zonenames += [link.to_node.zone for link in sock.links]
#                 fheights += [max([(bpy.data.objects[self.zone].matrix_world * vert.co)[2] for vert in bpy.data.objects[self.zone].data.vertices]) - (bpy.data.objects[link.to_node.zone].matrix_world * bpy.data.objects[link.to_node.zone].data.polygons[int(link.to_socket.sn)].center)[2] for link in sock.links]
#                 fareas += [facearea(bpy.data.objects[link.to_node.zone], bpy.data.objects[link.to_node.zone].data.polygons[int(link.to_socket.sn)]) for link in sock.links]

#             self['zonenames'] = zonenames

#             for z, zn in enumerate(self['zonenames']):
#                 self['Distance {}'.format(z)] = fheights[z]
#                 self['Relative Ratio {}'.format(z)] = 1.0
#                 self['Cross Section {}'.format(z)] = fareas[z]

#         for sock in self.outputs:
#             socklink(sock, self.id_data.name)

#     def uvsockupdate(self):
#         for sock in self.outputs:
#             socklink(sock, self.id_data.name)

#             if sock.bl_idname == 'EnViBoundSocket':
#                 uvsocklink(sock, self.id_data.name)

#     def epwrite(self):
#         scheduled = 1 if self.inputs['Schedule'].links and not self.inputs['Schedule'].links[0].to_node.use_custom_color else 0
#         paramvs = ('{}_TC'.format(self.zone), self.zone, ('', '{}_TCSched'.format(self.zone))[scheduled], self.waw, self.ocs, self.odc)
#         params = ('Name of Thermal Chimney System', 'Name of Thermal Chimney Zone', 'Availability Schedule Name', 'Width of the Absorber Wall',
#                   'Cross Sectional Area of Air Channel Outlet', 'Discharge Coefficient')

#         for z, zn in enumerate(self['zonenames']):
#             params += (' Zone Name {}'.format(z + 1), 'Distance from the Top of the Thermal Chimney to Inlet {}'.format(z + 1), 'Relative Ratios of Air Flow Rates Passing through Zone {}'.format(z + 1),
#                        'Cross Sectional Areas of Air Channel Inlet {}'.format(z + 1))
#             paramvs += (zn, self['Distance {}'.format(z)], self['Relative Ratio {}'.format(z)], self['Cross Section {}'.format(z)])

#         return epentry('ZoneThermalChimney', params, paramvs)