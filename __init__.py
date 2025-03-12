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

bl_info = {
    "name": "VI-Suite",
    "author": "Ryan Southall",
    "version": (0, 7, 0),
    "blender": (4, 3),
    "location": "Node Editor & 3D View > Properties Panel",
    "description": "Radiance/EnergyPlus/OpenFOAM exporter and results visualiser",
    "warning": "This is a beta script. Some functionality is buggy",
    "tracker_url": "https://github.com/rgsouthall/vi-suite07/issues",
    "doc_url": "http://blogs.brighton.ac.uk/visuite/",
    "category": "Import-Export"}

if "bpy" in locals():
    import imp
    imp.reload(vi_operators)
    imp.reload(vi_ui)
    imp.reload(vi_func)
    imp.reload(vi_node)
    imp.reload(envi_mat)
else:
    import sys, os, inspect, shlex, bpy, requests, shutil, glob
    from subprocess import Popen, call
    import nodeitems_utils
    from bpy.app.handlers import persistent
    from bpy.props import StringProperty, EnumProperty, IntProperty, FloatProperty, BoolProperty
    from bpy.types import AddonPreferences, Image, Material
    addonpath = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

    if sys.platform == 'linux':
        os.environ['SDL_VIDEO_X11_NET_WM_BYPASS_COMPOSITOR'] = '0'
        os.environ['PYTHON_INCLUDE_DIR'] = os.path.join(addonpath, 'Python', sys.platform, 'include')

    #os.environ["KIVY_NO_CONSOLELOG"] = "1"
    install_fails = []
    sys_install = 0
    
    def install_libs():
        install_fails = []

        try:
            from PySide6.QtGui import QImage
            install_fails.append(0)
        except:
            install_fails.append(1)

        try:
            import matplotlib
            matplotlib.use('qtagg', force=True)
            import matplotlib.pyplot as plt
            plt.text(0, 0, 'dummy')
            plt.clf()
            install_fails.append(0)
        except Exception as e:
            install_fails.append(1)

        # try:
        #     from kivy.config import Config
        #     Config.set('kivy', 'log_level', 'error')
        #     Config.write()
        #     from kivy.app import App
        #     install_fails.append(0)
        # except:
        #     install_fails.append(0)

        try:
            import scipy
            install_fails.append(0)
        except:
            install_fails.append(1)

        return install_fails
    
    if not all(install_libs()):
        sys_install = 1

    if any(install_libs()):
        print('Some system libraries were not found. Checking built-in libraries')
        sys.path.insert(0, os.path.join(addonpath, 'Python', sys.platform, '{}ib'.format(('l', 'L')[sys.platform == 'win32']),
                                        ('python{}.{}'.format(sys.version_info.major, sys.version_info.minor), '')[sys.platform == 'win32'],
                                        'site-packages'))
        sys.path.insert(0, os.path.join(addonpath, 'Python', sys.platform))

        if os.environ.get('PYTHONPATH'):
            if os.path.join(addonpath, 'Python', sys.platform) not in os.environ['PYTHONPATH']:
                os.environ['PYTHONPATH'] = os.path.join(addonpath, 'Python', sys.platform) + os.pathsep + os.environ['PYTHONPATH']
                
                if sys.platform != 'linux':
                    os.environ['PYTHONPATH'] = os.path.join(addonpath, 'Python', sys.platform) + os.pathsep + os.path.join(addonpath, 'Python', sys.platform, '{}ib'.format(('l', 'L')[sys.platform == 'win32']),
                                                                        ('python{}.{}'.format(sys.version_info.major, sys.version_info.minor), '')[sys.platform == 'win32'],
                                                                        'site-packages') + os.pathsep + os.environ['PYTHONPATH']
        else:
            os.environ['PYTHONPATH'] = os.path.join(addonpath, 'Python', sys.platform)
            if sys.platform != 'linux':
                os.environ['PYTHONPATH'] += os.pathsep + os.path.join(addonpath, 'Python', sys.platform, '{}ib'.format(('l', 'L')[sys.platform == 'win32']),
                                                                    ('python{}.{}'.format(sys.version_info.major, sys.version_info.minor), '')[sys.platform == 'win32'],
                                                                    'site-packages')
        
        if sys.platform == 'linux':
            if not os.environ.get('LD_LIBRARY_PATH'):
                os.environ['LD_LIBRARY_PATH'] = os.path.join(addonpath, 'Python', sys.platform)

            elif os.path.join(addonpath, 'Python', sys.platform) not in os.environ['LD_LIBRARY_PATH']:
                os.environ['LD_LIBRARY_PATH'] += os.pathsep + os.path.join(addonpath, 'Python', sys.platform)
            
        elif sys.platform == 'darwin':
            if not os.environ.get('DYLD_LIBRARY_PATH'):
                os.environ['DYLD_LIBRARY_PATH'] = os.path.join(addonpath, 'Python', sys.platform)
        
        if os.environ.get('PATH'):
            if os.path.join(addonpath, 'Python', sys.platform, 'bin') not in os.environ['PATH']:
                os.environ['PATH'] += os.pathsep + os.path.join(addonpath, 'Python', sys.platform, 'bin')
        else:
            os.environ['PATH'] = os.path.join(addonpath, 'Python', sys.platform, 'bin')

        install_fails = install_libs()

        if any(install_fails):
            print('Some built-in libraries were not found. Installing.')
            pyqt_cmd = '"{}" -m pip install -r "{}" --target "{}"'.format(sys.executable, os.path.join(addonpath, 'Python', 'requirements.txt'), os.path.join(addonpath, 'Python', sys.platform))
            Popen(shlex.split(pyqt_cmd)).wait()
            print('Installed local libraries')

        else:
            print('Found built-in libraries')
    
    else:
        print('VI-Suite: Using system libraries')

    try:
        import netgen
    except:
        if sys.platform == 'darwin':
            ngocc_cmd = '"{0}" -m pip install --upgrade --force --prefix "{1}" netgen-occt==7.8.1'.format(sys.executable, os.path.join(addonpath, 'Python', sys.platform))
            Popen(shlex.split(ngocc_cmd)).wait()
            ng_cmd = '"{0}" -m pip install --upgrade --force --target "{1}" netgen-mesher==6.2.2406.post119.dev1'.format(sys.executable, os.path.join(addonpath, 'Python', sys.platform))
            Popen(shlex.split(ng_cmd)).wait()
            src_path = os.path.join(addonpath, 'Python', sys.platform, 'lib')
            dest_path = os.path.join(addonpath, 'Python')
            files = os.listdir(src_path)

            for f in files:
                if sys.platform == 'darwin' and f.endswith('.dylib') or sys.platform == 'linux' and f.startswith('libT'):
                    src_file = os.path.join(src_path, f)
                    dest_file = os.path.join(dest_path, f)
                    shutil.move(src_file, dest_file)

        elif sys.platform == 'win32':
            ngocc_cmd = '"{0}" -m pip install --target "{1}" netgen-occt==7.8.1'.format(sys.executable, addonpath)
            Popen(shlex.split(ngocc_cmd)).wait()
            ng_cmd = '"{0}" -m pip install --target "{1}" netgen-mesher==6.2.2406.post119.dev1'.format(sys.executable, os.path.join(addonpath, 'Python', sys.platform))
            Popen(shlex.split(ng_cmd)).wait()
        
        elif sys.platform == 'linux':
            if not sys_install:
                ng_cmd = '"{0}" -m pip install --upgrade netgen-mesher==6.2.2406.post119.dev1'.format(sys.executable)
                Popen(shlex.split(ng_cmd)).wait()
            else:
                print("Blender is using the system's Python installation, but no system installation of netgen was found")
        try:
            import netgen
        except Exception as e:
            print('Netgen installation failed and is disabled')

    try:
        import pyroomacoustics as pra
    except:
        if sys.platform == 'linux':
            print('For pyroomacoustics functionality on linux, a system install of Blender, PySide6, Matplotlib, Netgen and pyroomacoustics is required')
        else:
            try:
                import pyroomacoustics as pra
            except Exception as e:
                print('pyroomacoustics installation failed and is disabled')
    
    if sys.platform in ('linux', 'darwin'):
        ep_path = os.path.join(addonpath, 'EPFiles', sys.platform)
        rad_path = os.path.join(addonpath, 'RadFiles', sys.platform, 'bin')

        for fn in ('cnt', 'epw2wea', 'evalglare', 'falsecolor', 'genBSDF', 'gendaylit', 'gendaymtx', 'gensky',
                   'getbbox', 'getinfo', 'ies2rad', 'mkpmap', 'obj2mesh', 'oconv', 'pcomb', 'pcompos', 'pcond',
                   'pfilt', 'pkgBSDF', 'pmapdump', 'psign', 'rad2mgf', 'rcalc', 'rcontrib', 'rfluxmtx', 'rmtxop',
                   'rpict', 'rpiece', 'rtrace', 'rttree_reduce', 'rvu', 'vwrays', 'wrapBSDF', 'xform', 'dcglare'):
            try:
                if not os.access(os.path.join(rad_path, fn), os.X_OK):
                    os.chmod(os.path.join(rad_path, fn), 0o775)
            except Exception:
                print('{} not found'.format(fn))

        for fn in ('energyplus-22.1.0', 'ExpandObjects'):
            try:
                if not os.access(os.path.join(ep_path, fn), os.X_OK):
                    os.chmod(os.path.join(ep_path, fn), 0o775)
            except Exception:
                print('{} not found'.format(fn))
        

        if not os.path.islink(os.path.join(ep_path, 'energyplus')):
            os.symlink(os.path.join(ep_path, 'energyplus-22.1.0'), os.path.join(ep_path, 'energyplus'))
        elif not os.path.isfile(os.path.join(ep_path, 'energyplus')):
            os.remove(os.path.join(ep_path, 'energyplus'))
            os.symlink(os.path.join(ep_path, 'energyplus-22.1.0'), os.path.join(ep_path, 'energyplus'))

        if not os.path.islink(os.path.join(ep_path, 'libenergyplusapi.{}'.format(('so', 'dylib')[sys.platform == 'darwin']))):
            os.symlink(os.path.join(ep_path, 'libenergyplusapi{}.22.1.0{}'.format(('.so', '')[sys.platform == 'darwin'], ('', '.dylib')[sys.platform == 'darwin'])),
                                    os.path.join(ep_path, 'libenergyplusapi.{}'.format(('so', 'dylib')[sys.platform == 'darwin'])))
        elif not os.path.isfile(os.path.join(ep_path, 'libenergyplusapi.{}'.format(('so', 'dylib')[sys.platform == 'darwin']))):
            os.remove(os.path.join(ep_path, 'libenergyplusapi.{}'.format(('so', 'dylib')[sys.platform == 'darwin'])))
            os.symlink(os.path.join(ep_path, 'libenergyplusapi{}.22.1.0{}'.format(('.so', '')[sys.platform == 'darwin'], ('', '.dylib')[sys.platform == 'darwin'])),
                                    os.path.join(ep_path, 'libenergyplusapi.{}'.format(('so', 'dylib')[sys.platform == 'darwin'])))

    
    from .vi_node import vinode_categories, envinode_categories, envimatnode_categories, ViNetwork, No_Loc, So_Vi_Loc
    from .vi_node import No_Vi_SP, No_Vi_WR, No_Vi_SVF, So_Vi_Res, No_Vi_SS
    from .vi_node import No_Li_Geo, No_Li_Con, No_Li_Sen, So_Li_Geo, So_Li_Con, No_Text, So_Text, No_CSV
    from .vi_node import No_Vi_Im, No_Li_Im, So_Li_Im, No_Li_Gl, No_Li_Fc
    from .vi_node import No_Li_Sim, No_ASC_Import, So_Flo_Mesh
    from .vi_node import No_En_Net_Zone, No_En_Net_Occ, So_En_Net_Eq, So_En_Net_Inf, So_En_Net_Hvac, No_En_Net_Hvac
    from .vi_node import No_En_Geo, So_En_Geo, EnViNetwork, EnViMatNetwork, No_En_Con, So_En_Con
    from .vi_node import No_En_Mat_Con, No_En_Mat_Sc, No_En_Mat_Sh, No_En_Mat_ShC, No_En_Mat_Bl, No_En_Mat_SG, No_En_Mat_Op, No_En_Mat_Tr, No_En_Mat_Gas, So_En_Mat_Ou, So_En_Mat_Op, No_En_Mat_Sched
    from .vi_node import So_En_Net_Occ, So_En_Net_Sched, So_En_Mat_Sched, No_En_Net_Sched, No_En_Sim, No_Vi_Chart, No_Vi_HMChart, So_En_Res, So_En_ResU, So_En_Net_TSched, No_En_Net_Eq, No_En_Net_Inf
    from .vi_node import No_En_Net_SFlow, No_En_Net_SSFlow, So_En_Net_SFlow, So_En_Net_SSFlow, So_En_Mat_PV, No_En_Mat_PV
    from .vi_node import So_En_Mat_PVG, No_En_Mat_PVG, No_Vi_Metrics, So_En_Mat_Tr, So_En_Mat_Gas, So_En_Mat_Fr, So_En_Net_Bound, No_En_Net_ACon, No_En_Net_Ext
    from .vi_node import No_En_Net_EMSZone, No_En_Net_Prog, No_En_Net_EMSPy, So_En_Net_Act, So_En_Net_Sense, No_Flo_Case, So_Flo_Case, No_Flo_NG, So_Flo_Con, No_Flo_Bound, No_Flo_Sim
    from .vi_node import No_En_IF, No_En_RF, So_En_Net_WPC, No_En_Net_Azi, No_Anim, So_Anim, No_En_Net_Anim, No_En_Mat_Anim, ViEnRIn
    from .vi_node import So_En_Mat_Sh, So_En_Mat_ShC, So_En_Mat_Sc, No_Vi_EC
    from .vi_node import No_Au_Sim, No_Au_Conv, So_Au_IR, So_Au_Scene
    from .vi_func import iprop, bprop, eprop, fprop, sprop, fvprop, sunpath1
    from .vi_func import lividisplay, logentry, ob_to_stl, ec_update
    from .livi_func import rtpoints, lhcalcapply, udidacalcapply, basiccalcapply, adgpcalcapply, radmat, retsv
    from .envi_func import enunits, enpunits, enparametric, resnameunits, aresnameunits
    from .envi_mat import envi_elayertype, envi_eclasstype, envi_emattype, envi_embodied
    from .flovi_func import fvmat, ret_fvbp_menu, ret_fvbu_menu, ret_fvbnut_menu, ret_fvbk_menu, ret_fvbepsilon_menu
    from .flovi_func import ret_fvb_menu, ret_fvbomega_menu, ret_fvbt_menu, ret_fvba_menu, ret_fvbprgh_menu, ret_fvrad_menu, ret_fvi_menu
    from .vi_operators import NODE_OT_WindRose, NODE_OT_SVF, NODE_OT_En_Con, NODE_OT_En_Sim, NODE_OT_TextUpdate
    from .vi_operators import MAT_EnVi_Node, NODE_OT_Shadow, NODE_OT_CSV, NODE_OT_ASCImport, NODE_OT_FileSelect, NODE_OT_HdrSelect
    from .vi_operators import NODE_OT_Li_Geo, NODE_OT_Li_Con, NODE_OT_Li_Pre, NODE_OT_Li_Sim, NODE_OT_EC, OBJECT_OT_EcS, OBJECT_OT_EcE, NODE_OT_ECPie, NODE_OT_WLCLine, NODE_OT_COMLine
    from .vi_operators import NODE_OT_Li_Im, NODE_OT_Li_Gl, NODE_OT_Li_Fc, NODE_OT_En_Geo, OBJECT_OT_VIGridify, OBJECT_OT_Embod, NODE_OT_En_UV, NODE_OT_En_EC, MAT_EnVi_Node_Remove
    from .vi_operators import NODE_OT_Chart, NODE_OT_HMChart, NODE_OT_En_PVA, NODE_OT_En_PVS, NODE_OT_En_LayS, NODE_OT_En_EcS, NODE_OT_En_ConS, TREE_OT_goto_mat, TREE_OT_goto_group
    from .vi_operators import OBJECT_OT_Li_GBSDF, OBJECT_OT_GOct, MATERIAL_OT_Li_LBSDF, MATERIAL_OT_Li_SBSDF, MATERIAL_OT_Li_DBSDF, NODE_OT_EcE
    from .vi_operators import NODE_OT_Flo_Case, NODE_OT_Flo_NG, NODE_OT_Flo_Bound, NODE_OT_Flo_Sim, NODE_OT_Au_Rir, NODE_OT_WavSelect, NODE_OT_Au_Conv, NODE_OT_Au_Play, NODE_OT_Au_Stop
    from .vi_operators import NODE_OT_Au_PlayC, NODE_OT_Au_Save
    from .vi_display import VIEW3D_OT_WRDisplay, VIEW3D_OT_SVFDisplay, VIEW3D_OT_Li_BD, VIEW3D_OT_Li_DBSDF, VIEW3D_OT_SSDisplay, VIEW3D_OT_RTDisplay, NODE_OT_SunPath, NODE_OT_Vi_Info
    from .vi_display import script_update, col_update, leg_update, w_update, t_update, livires_update, e_update
    from .vi_ui import VI_PT_3D, VI_PT_Mat, VI_PT_Ob, VI_PT_Col, VI_PT_Gridify, TREE_PT_envim, TREE_PT_envin, TREE_PT_vi
    from .vi_dicts import colours
    from .auvi_mat import auvi_type_abs, auvi_type_scatt, auvi_abs, auvi_scatt, auvi_materials


def abspath(self, context):
    if self.radbin != bpy.path.abspath(self.radbin):
        self.radbin = bpy.path.abspath(self.radbin)
    if self.radlib != bpy.path.abspath(self.radlib):
        self.radlib = bpy.path.abspath(self.radlib)
    if self.epbin != bpy.path.abspath(self.epbin):
        self.epbin = bpy.path.abspath(self.epbin)
    if self.epweath != bpy.path.abspath(self.epweath):
        self.epweath = bpy.path.abspath(self.epweath)
    if self.ofbin != bpy.path.abspath(self.ofbin):
        self.ofbin = bpy.path.abspath(self.ofbin)
    if self.oflib != bpy.path.abspath(self.oflib):
        self.oflib = bpy.path.abspath(self.oflib)
    if self.ofetc != bpy.path.abspath(self.ofetc):
        self.ofetc = bpy.path.abspath(self.ofetc)
    if self.datab != bpy.path.abspath(self.datab):
        self.datab = bpy.path.abspath(self.datab)
    
    if context.scene:
        svp = context.scene.vi_params

        if not svp.get('viparams'):
            svp['viparams'] = {}

        if self.datab:
            svp['viparams']['datab'] = bpy.path.abspath(self.datab)

    path_update()


def flovi_levels(self, context):
    if self.flovi_slmin > self.flovi_slmax:
        self.flovi_slmin -= 1


def unititems(self, context):
    try:
        scene = context.scene
        svp = scene.vi_params

        if svp['liparams']['unit'] == 'W/m2 (f)':
            return [('firradm2', 'W/m2 (f)', 'Full spectrum irradiance per metre square'),
                    ('firrad', 'W (f)', 'Full spectrum irradiance')]
        elif svp['liparams']['unit'] == 'Lux':
            return [('illu', 'Lux', 'Illuminance'),
                    ('virradm2', 'W/m2 (v)', 'Visible spectrum irradiance')]
        elif svp['liparams']['unit'] == 'DF (%)':
            return [('df', 'DF (%)', 'Daylight factor'),
                    ('virradm2', 'W/m2 (v)', 'Visible spectrum irradiance')]
        elif svp['liparams']['unit'] == 'klxh':
            return [('illuh', 'klux-hours', 'kilolux-hours'), ('virradh', 'kWh (v)', 'kilo-Watt hours (visible spectrum)'), ('virradhm2', 'kWh/m2 (v)', 'kilo-Watt hours per square metre (visible spectrum)')]
        elif svp['liparams']['unit'] == 'kWh (f)':
            return [('firradh', 'kWh (f)', 'kilo-Watt hours (solar spectrum)'),
                    ('firradhm2', 'kWh/m2 (f)', 'kilo-Watt hours per square metre (solar spectrum)')]
        elif svp['liparams']['unit'] == 'DA (%)':
            return [("da", "DA", "Daylight Autonomy"),
                    ("udilow", "UDI (low)", "Useful daylight illuminance (low)"),
                    ("udisup", "UDI (supp)", "Useful daylight illuminance (supplemented)"),
                    ("udiauto", "UDI (auto)", "Useful daylight illuminance (autonomous)"),
                    ("udihi", "UDI (high)", "Useful daylight illuminance (high)"),
                    ("maxlux", "Lux level (max)", "Maximum lux level"),
                    ("avelux", "Lux level (ave)", "Average lux level"),
                    ("minlux", "Lux level (min)", "Minimum lux level")]
        elif svp['liparams']['unit'] == 'sDA (%)':
            return [("sda", "sDA", "Spatial Daylight Autonomy"),
                    ("ase", "ASE", "Annual Sunlight Exposure"),
                    ("sv", "Perimeter area", "Perimeter")]
        elif svp['liparams']['unit'] == 'EN17037_100 (% hrs)':
            return [("en100", "EN17037_100", "EN17037 100 lux"),
                    ("en300", "EN17037_300", "EN17037 300 lux")]
        elif svp['liparams']['unit'] == 'Sunlit (% hrs)':
            return [('sm', '% Sunlit', '% of time sunlit')]
        elif svp['liparams']['unit'] == 'SVF (%)':
            return [('svf', 'SVF (%)', '% of sky visible')]
        elif svp['liparams']['unit'] == 'W/m2':
            return [('firradm2', 'W/m2', 'Combined irradiance per metre square'),
                    ('firrad', 'W', 'Combined irradiance'),
                    ('firradrm2', 'W/m2 (red)', 'Red irradiance per metre square'),
                    ('firradr', 'W (red)', 'Red irradiance'),
                    ('firradgm2', 'W/m2 (green)', 'Green irradiance per metre square'),
                    ('firradg', 'W (green)', 'Green irradiance'),
                    ('firradbm2', 'W/m2 (blue)', 'Blue irradiance per metre square'),
                    ('firradb', 'W (blue)', 'Blue irradiance')]
        elif svp['liparams']['unit'] == 'GO (%)':
            return [('ago1v', 'GO (%)', 'Glare occurance'),
                    ('aga1v', 'GA (%)', 'Glare autonomy')]
        elif svp['liparams']['unit'] == 'RT60 (s)':
            return [('rt', 'RT60 (s)', 'Reverberation time (60db)'),
                    ('sti', 'STI', 'Speech Transmission Index'),
                    ('vol', 'TSL (dB)', 'Total sound level (db)')]
        else:
            return [('None', 'None', 'None')]

    except Exception as e:
        logentry(f'Error: {e}. Check the LiVi Result Type menu')
        return [('None', 'None', 'None')]


def ret_im(self, context):
    return [('None', 'None', 'No image')] + [(im.name, im.name, 'Blender image') for im in bpy.data.images]


def ret_envi_mats(self, context):
    return [(mat.name, mat.name, mat.name) for mat in bpy.data.materials if mat.vi_params.envi_nodes]


def bsdf_direcs(self, context):
    try:
        return [tuple(i) for i in context.scene.vi_params['liparams']['bsdf_direcs']]
    except Exception:
        return [('None', 'None', 'None')]

def get_sources(self, context):
    try:
        return [(i, i, i) for i in context.scene.vi_params['liparams']['sources']]
    except Exception:
        return [('None', 'None', 'None')]

def ga_views(self, context):
    try:
        return [(str(i + 1), f'View {i + 1}', f'View {i + 1}') for i in range(context.scene.vi_params['liparams']['views'])]
    except Exception:
        return [('None', 'None', 'None')]


class VIPreferences(AddonPreferences):
    bl_idname = __name__
    radbin: StringProperty(name='', description='Radiance binary directory location', default='', subtype='DIR_PATH', update=abspath)
    radlib: StringProperty(name='', description='Radiance library directory location', default='', subtype='DIR_PATH', update=abspath)
    epbin: StringProperty(name='', description='EnergyPlus binary directory location', default='', subtype='DIR_PATH', update=abspath)
    epweath: StringProperty(name='', description='EnergyPlus weather directory location', default='', subtype='DIR_PATH', update=abspath)
    ofbin: StringProperty(name='', description='OpenFOAM binary directory location', default='', subtype='DIR_PATH', update=abspath)
    oflib: StringProperty(name='', description='OpenFOAM library directory location', default='', subtype='DIR_PATH', update=abspath)
    ofetc: StringProperty(name='', description='OpenFOAM etc directory location', default='', subtype='DIR_PATH', update=abspath)
    datab: StringProperty(name='', description='Database directory', default='', subtype='DIR_PATH', update=abspath)
    ui_dict = {"Radiance bin directory:": 'radbin', "Radiance lib directory:": 'radlib', "EnergyPlus bin directory:": 'epbin',
               "EnergyPlus weather directory:": 'epweath', 'Database directory': 'datab', 'OpenFOAM bin directory': 'ofbin'}

    def draw(self, context):
        layout = self.layout

        # if not os.path.isdir(os.path.join(addonpath, 'Python', sys.platform, 'numpy')):
        #     row = layout.row()
        #     row.label(text='You do not seem to have an internet connection. Cannot download Python libraries')

        for entry in self.ui_dict:
            if entry == 'OpenFOAM bin directory' and sys.platform != 'linux':
                pass
            else:
                row = layout.row()
                row.label(text=entry)
                row.prop(self, self.ui_dict[entry])


def d_update(self, context):
    csvp = context.scene.vi_params
    
    if not self.vi_display and csvp.get('viparams'):
        dns = bpy.app.driver_namespace
        try:
            for d in csvp['viparams'].get('drivers'):
                if d in dns:
                    bpy.types.SpaceView3D.draw_handler_remove(dns[d], 'WINDOW')
                    logentry('Stopping {} display'.format(d))
        except:
            pass

        csvp['viparams']['drivers'] = []


class VI_Params_Scene(bpy.types.PropertyGroup):
    def get_frame(self):
        return self.id_data.frame_current

    def set_frame(self, value):
        if value > self['liparams']['fe']:
            self['vi_frames'] = self['liparams']['fe']
        elif value < self['liparams']['fs']:
            self['vi_frames'] = self['liparams']['fs']
        else:
            self.id_data.frame_set(value)
            self['vi_frames'] = value
        
        
    def update_frame(self, context):
        self.id_data.vi_params.vi_leg_max = self.id_data.vi_params.vi_leg_max

    def get_view(self):
        return self.get('vi_views', self['liparams']['views'])

    def set_view(self, value):
        if self['liparams'].get('views'):
            if value > self['liparams']['views']:
                self['vi_views'] = self['liparams']['views']
            elif value < 1:
                self['vi_views'] = 1
            else:
                self['vi_views'] = value
        else:
            if value < 1:
                self['vi_views'] = 1
            else:
                self['vi_views'] = value


    def disp_options(self, context):
        svp = context.scene.vi_params
        try:
            if svp['liparams']['unit'] == 'GO (%)':
                return [("0", "None", "No processing"), ("1", "3D", "3D results display"), 
                        ("2", "Interpolate", "Interpolate results display"), ("3", "Direction", "Directional results display")]
            else:
                return [("0", "None", "No processing"), ("1", "3D", "3D results display"), 
                        ("2", "Interpolate", "Interpolate results display")]
        except:
            return [("0", "None", "No processing"), ("1", "3D", "3D results display"), 
                    ("2", "Interpolate", "Interpolate results display")]

    vi_name = sprop("", "VI-Suite addon directory name", 1024, "")
    year: iprop("", 'Year', 2019, 2020, 2019)
    vipath: sprop("VI Path", "Path to files included with the VI-Suite ", 1024, addonpath)
    vi_frames: IntProperty(name="", description="Frame of parametric run", get=get_frame, set=set_frame, update=update_frame)
    solday: IntProperty(name="", description="Day of year", min=1, max=365, default=1, update=sunpath1)
    solhour: bpy.props.FloatProperty(name="", description="Time of day", subtype='TIME', unit='TIME', min=0, max=24, default=12, update=sunpath1)
    sp_hour_dash: fvprop(4, "", 'Main colour of the hour lines', [1.0, 0.0, 0.0, 0.0], 'COLOR', 0, 1)
    sp_hour_main: fvprop(4, "", 'Dash colour of the hour lines', [1.0, 1.0, 0.0, 1.0], 'COLOR', 0, 1)
    sp_season_main: fvprop(4, "", 'Main colour of the season lines', [1.0, 0.0, 0.0, 1.0], 'COLOR', 0, 1)
    sp_season_dash: fvprop(4, "", 'Dash colour of the season lines', [1.0, 1.0, 1.0, 0.0], 'COLOR', 0, 1)
    sp_sun_colour: fvprop(4, "", 'Sun colour', [1.0, 1.0, 0.0, 1.0], 'COLOR', 0, 1)
    sp_globe_colour: fvprop(4, "", 'Sun colour', [0.0, 0.0, 1.0, 0.25], 'COLOR', 0, 1)
    sp_sun_angle: FloatProperty(name="", description="Sun size", min=0, max=1, default=0.01, update=sunpath1)
    sp_sun_size: iprop("", 'Sun size', 1, 50, 10)
    sp_sun_strength: FloatProperty(name="", description="Sun strength", min=0, max=100, default=1.0, update=sunpath1)
    sp_season_dash_ratio: fprop("", "Ratio of line to dash of season lines", 0, 5, 0)
    sp_hour_dash_ratio: fprop("", "Ratio of line to dash of hour lines", -1, 1, 0.5)
    sp_hour_dash_density: fprop("", "Ratio of line to dash of hour lines", 0, 5, 1)
    sp_line_width: iprop("", "Sun path line width", 0, 50, 2)
    latitude: FloatProperty(name="", description="Site decimal latitude (N is positive)",
                            min=-89.99, max=89.99, default=52.0, update=sunpath1)
    longitude: FloatProperty(name="", description="Site decimal longitude (E is positive)",
                             min=-180, max=180, default=0.0, update=sunpath1)
    sp_suns: EnumProperty(items=[('0', 'Single', 'Single sun'),
                                 ('1', 'Monthly', 'Monthly sun for chosen time'),
                                 ('2', 'Hourly', 'Hourly sun for chosen date')],
                          name='', description='Sunpath sun type', default='0', update=sunpath1)
    sp_sst: FloatProperty(name="", description="Sun strength", min=0, max=100, default=0.1, update=sunpath1)
    sp_ssi: FloatProperty(name="", description="Sun size", min=0, max=1, default=0.01, update=sunpath1)
    sp_sd: IntProperty(name="", description="Day of year", min=1, max=365, default=1, options={'ANIMATABLE'}, update=sunpath1)
    sp_sh: FloatProperty(name="", description="Time of day", subtype='TIME', unit='TIME',
                         min=0, max=24, default=12, options={'ANIMATABLE'}, update=sunpath1)
    sp_hd: bprop("", "", 0)
    sp_up: bprop("", "", 0)
    sp_td: bprop("", "", 0)
    li_disp_panel: iprop("Display Panel", "Shows the Display Panel", -1, 2, 0)
    li_disp_menu: EnumProperty(items=unititems, name="", description="LiVi metric selection", update=livires_update)
    vi_display_rp_fsh: fvprop(4, "", "Font shadow", [0.0, 0.0, 0.0, 1.0], 'COLOR', 0, 1)
    vi_display_rp_fs: iprop("", "Point result font size", 4, 64, 24)
    vi_display_rp_fc: fvprop(4, "", "Font colour", [0.0, 0.0, 0.0, 1.0], 'COLOR', 0, 1)
    vi_display_rp_sh: bprop("", "Toggle for font shadow display",  False)
    vi_display: BoolProperty(name="", description="Toggle results display", default=0, update=d_update)
    vi_disp_process: EnumProperty(items=disp_options,  name="", description="Result processing")
    vi_disp_pos: EnumProperty(items=[("0", "Positive", "Positive face placement"), ("1", "Negative", "Negative face placement")],  name="", 
                              description="Face placement to avoid z-fighting", default="0", update=leg_update)
    vi_leg_unit: sprop("", "Legend unit", 1024, "")
    vi_leg_max: FloatProperty(name="", description="Legend maximum", min=0, max=1000000, default=1000, update=leg_update)
    vi_leg_min: FloatProperty(name="", description="Legend minimum", min=0, max=1000000, default=0, update=leg_update)
    vi_leg_col: EnumProperty(items=colours, name="", description="Legend colours", default='rainbow', update=col_update)
    vi_leg_levels: IntProperty(name="", description="Day of year", min=2, max=100, default=20, update=leg_update)
    vi_arrow_size: FloatProperty(name="", description="Arrow size", min=0.01, max=10, default=0.5, update=leg_update)
    vi_leg_scale: EnumProperty(items=[('0', 'Linear', 'Linear scale'), ('1', 'Log', 'Logarithmic scale')], name="", description="Legend scale", default='0', update=leg_update)
    wind_type: eprop([("0", "Speed", "Wind Speed (m/s)"), ("1", "Direction", "Wind Direction (deg. from North)")], "", "Wind metric", "0")
    vi_disp_trans: FloatProperty(name="", description="Sensing material transparency", min=0, max=1, default=0, update=t_update)
    vi_disp_wire: BoolProperty(name="", description="Draw wire frame", default=0, update=w_update)
    vi_disp_mat: BoolProperty(name="", description="Turn on/off result material emission", default=0, update=col_update)
    vi_disp_ems: FloatProperty(name="", description="Emissive strength", default=1, min=0, update=col_update)
    vi_scatt_col: EnumProperty(items=colours, name="", description="Scatter colour", default='rainbow')
    vi_disp_refresh: bprop("", "Refresh display",  False)
    vi_res_mod: sprop("", "Result modifier", 1024, "")
    vi_res_process: EnumProperty(items=[("0", "None", ""), ("1", "Modifier", ""), ("2", "Script", "")], name="", description="Specify the type of data processing", default="0", update=script_update)
    script_file: StringProperty(description="Text file to show", update=script_update)
    ss_disp_panel: iprop("Display Panel", "Shows the Display Panel", -1, 2, 0)
    vi_display_rp: bprop("", "", False)
    vi_display_rp_off: fprop("", "Surface offset for number display", 0, 5, 0.001)
    vi_display_sel_only: bprop("", "", False)
    vi_display_vis_only: bprop("", "", False)
    vi_disp_3dlevel: FloatProperty(name="", description="Level of 3D result plane extrusion", min=0, max=500, default=0, update=e_update)
    vi_bsdf_direc: EnumProperty(items=bsdf_direcs, name="", description="BSDf display direction")
    vi_bsdfleg_max: FloatProperty(name="", description="Legend maximum", min=0, max=1000000, default=100)
    vi_bsdfleg_min: FloatProperty(name="", description="Legend minimum", min=0, max=1000000, default=0)
    vi_bsdfleg_scale: EnumProperty(items=[('0', 'Linear', 'Linear scale'), ('1', 'Log', 'Logarithmic scale')], name="", description="Legend scale", default='0')
    vi_bsdf_font: iprop("", "Font size for BSDF numerical display", 0, 48, 0)
    vi_views: IntProperty(name="", description="GA views", get=get_view, set=set_view, update=livires_update)
    en_disp_type: EnumProperty(items=enparametric, name="", description="Type of EnVi display")
    vi_nodes: bpy.props.PointerProperty(type=bpy.types.NodeTree)
    envi_nodes: bpy.props.PointerProperty(type=bpy.types.NodeTree)
    resas_disp: bprop("", "", False)
    reszt_disp: bprop("", "", False)
    reszh_disp: bprop("", "", False)
    resaa_disp: bprop("", "", False)

    (resaa_disp, resaws_disp, resawd_disp, resah_disp, resas_disp, reszt_disp, reszh_disp, reszhw_disp, reszcw_disp, reszsg_disp, reszppd_disp,
     reszpmv_disp, resvls_disp, resvmh_disp, resim_disp, resiach_disp, reszco_disp, resihl_disp, reszlf_disp, reszof_disp, resmrt_disp,
     resocc_disp, resh_disp, resfhb_disp, reszahw_disp, reszacw_disp, reshrhw_disp, restcvf_disp, restcmf_disp, restcot_disp, restchl_disp,
     restchg_disp, restcv_disp, restcm_disp, resldp_disp, resoeg_disp, respve_disp, respvw_disp, respveff_disp, respvt_disp) = resnameunits()

    (resazmaxt_disp, resazmint_disp, resazavet_disp,
     resazmaxhw_disp, resazminhw_disp, resazavehw_disp,
     resazth_disp, resazthm_disp,
     resazmaxcw_disp, resazmincw_disp, resazavecw_disp,
     resaztc_disp, resaztcm_disp,
     resazmaxco_disp, resazaveco_disp, resazminco_disp,
     resazlmaxf_disp, resazlminf_disp, resazlavef_disp,
     resazmaxshg_disp, resazminshg_disp, resazaveshg_disp,
     resaztshg_disp, resaztshgm_disp) = aresnameunits()

    au_sources: EnumProperty(items=get_sources, name="", description="Au Sources", update=livires_update)


class VI_Params_Object(bpy.types.PropertyGroup):
    # VI-Suite object definitions
    vi_type_string: sprop("", "VI Suite object type", 1024, "")
    vi_type: eprop([("0", "None", "Not a VI-Suite specific object"),
                    ("1", "EnVi Surface", "Designates an EnVi surface"),
                    ("2", "CFD Domain", "Specifies an OpenFoam BlockMesh"),
                    ("3", "CFD Geometry", "Specifies an OpenFoam geometry"),
                    ("4", "Light Array", "Specifies a LiVi lighting array"),
                    ("5", "Complex Fenestration", "Specifies complex fenestration for BSDF generation"),
                    ("6", "CFD Probe", "Specifies a mesh as a CFD probe")],
                   "", "Specify the type of VI-Suite zone", "0")

    # LiVi object properties
    livi_merr: bprop("LiVi simple mesh export", "Boolean for simple mesh export", False)
    ies_name: bpy.props.StringProperty(name="", description="Name of the IES file", default="", subtype="FILE_PATH")
    ies_strength: fprop("", "IES brightness multiplier", 0, 1, 1)
    ies_unit: eprop([("m", "Meters", ""), ("c", "Centimeters", ""), ("f", "Feet", ""), ("i", "Inches", "")], "", "Specify the IES file measurement unit", "m")
    ies_colmenu: eprop([("0", "RGB", ""), ("1", "Temperature", "")], "", "Specify the IES colour type", "0")
    ies_rgb: fvprop(3, "", 'IES Colour', [1.0, 1.0, 1.0], 'COLOR', 0, 1)
    ies_ct: iprop("", "Colour temperature in Kelven", 0, 12000, 4700)
    limerr: bprop("", "", False)
    manip: bprop("", "", False)
    bsdf_proxy: bprop("", "", False)
    basiccalcapply = basiccalcapply
    rtpoints = rtpoints
    udidacalcapply = udidacalcapply
    lividisplay = lividisplay
    lhcalcapply = lhcalcapply
    adgpcalcapply = adgpcalcapply
    li_bsdf_direc: EnumProperty(items=[('+b -f', 'Backwards', 'Backwards BSDF'),
                                       ('+f -b', 'Forwards', 'Forwards BSDF'),
                                       ('+b +f', 'Bi-directional', 'Bi-directional BSDF')], name='', description='BSDF direction', default='+b -f')
    li_bsdf_proxy: bprop("", "Include proxy geometry in the BSDF", False)
    li_bsdf_dimen: EnumProperty(items=[('meter', 'Meters', 'BSDF length unit'),
                                       ('millimeter', 'Millimeters', 'BSDF length unit'),
                                       ('centimeter', 'Centimeters', 'BSDF length unit'),
                                       ('foot', 'Feet', 'BSDF length unit'),
                                       ('inch', 'Inches', 'BSDF length unit')], name='', description='BSDF length unit', default='meter')
    li_bsdf_tensor: EnumProperty(items=[(' ', 'Klems', 'Uniform Klems sample'),
                                        ('-t3', 'Symmentric', 'Symmetric Tensor BSDF'),
                                        ('-t4', 'Assymmetric', 'Asymmetric Tensor BSDF')], name='', description='BSDF tensor', default=' ')
    li_bsdf_res: EnumProperty(items=[('1', '2x2', '2x2 sampling resolution'),
                                     ('2', '4x4', '4x4 sampling resolution'),
                                     ('3', '8x8', '8x8 sampling resolution'),
                                     ('4', '16x16', '16x16 sampling resolution'),
                                     ('5', '32x32', '32x32 sampling resolution'),
                                     ('6', '64x64', '64x64 sampling resolution'),
                                     ('7', '128x128', '128x128 sampling resolution')], name='', description='BSDF resolution', default='4')
    li_bsdf_tsamp: IntProperty(name='', description='Tensor samples per region', min=1, max=10000, default=4)
    li_bsdf_ksamp: IntProperty(name='', description='Klem samples', min=1, default=2000)
    li_bsdf_rcparam: sprop("", "rcontrib parameters", 1024, "")
    bsdf_running: bprop("", "Running BSDF calculation", False)
    retsv = retsv
    # envi_type: eprop([("0", "Construction", "Thermal Construction"), ("1", "Shading", "Shading Object")], "", "Specify the EnVi surface type", "0")
    flovi_solver: EnumProperty(items=[('icoFoam', 'IcoFoam', 'Transient laminar solver'), ('simpleFoam', 'SimpleFoam', 'Transient turbulent solver'),
                                      ('bBSimpleFoam', 'buoyantBoussinesqSimpleFoam', 'Bouyant Boussinesq Turbulent solver'), ('bSimpleFoam', 'buoyantSimpleFoam', 'Bouyant Turbulent solver')],
                               name="", default='icoFoam')
    flovi_turb: EnumProperty(items=[('kEpsilon', 'K-Epsilon', 'K-Epsion turbulence model'), ('kOmega', 'K-Omega', 'K-Omega turbulence model'),
                                    ('SpalartAllmaras', 'SpalartAllmaras', 'SpalartAllmaras turbulence model')],
                             name="", default='kEpsilon')
    flovi_fl: IntProperty(name='', description='SnappyHexMesh object features levels', min=1, max=20, default=4)
    flovi_slmax: IntProperty(name='', description='SnappyHexMesh surface maximum levels', min=1, max=20, default=4, update=flovi_levels)
    flovi_slmin: IntProperty(name='', description='SnappyHexMesh surface minimum levels', min=1, max=20, default=3, update=flovi_levels)
    flovi_sl: iprop('', 'SnappyHexMesh surface minimum levels', 0, 20, 3)
    mesh: bprop("", "Radiance mesh geometry export", 1)
    triangulate: bprop("", "Triangulate mesh geometry for export", 0)
    flovi_ufield: fvprop(3, '', 'Velocity field value', [0, 0, 0], 'VELOCITY', -100, 100)
    flovi_pfield: fprop("", "p field value", 0, 500, 0)
    flovi_nutfield: fprop("", "nut field value", 0, 500, 0)
    flovi_kfield: fprop("", "k field value", 0, 500, 1.5)
    flovi_efield: fprop("", "e field value", 0, 500, 0.03)
    flovi_ofield: fprop("", "o field value", 0, 500, 0.03)
    flovi_probe: bprop("", "OpenFoam probe", False)
    embodied: BoolProperty(name="", description="Embodied carbon", default=0, update=ec_update)
    embodiedtype: EnumProperty(items=envi_elayertype, name="", description="Layer embodied material type", update=ec_update)
    embodiedclass: EnumProperty(items=envi_eclasstype, name="", description="Layer embodied class", update=ec_update)
    embodiedmat: EnumProperty(items=envi_emattype, name="", description="Layer embodied material", update=ec_update)
    ec_id: StringProperty(name="", description="Embodied id (unique indentifier")
    ec_type: StringProperty(name="", description="Embodied type e.g. Insulation")
    ec_class: StringProperty(name="", description="Embodied class (class of type e.g. phenolic foam)")
    ec_name: StringProperty(name="", description="Embodied name")
    ec_unit: EnumProperty(items=[("kg", "kg", "per kilogram"),
                                  ("m2", "m2", "per square metre"),
                                  ("m3", "m3", "per cubic metre"),
                                  ("each", "Each", "per item")],
                                  name="",
                                  description="Embodied carbon unit",
                                  default="kg")
    ec_amount: FloatProperty(name="", description="EC amount of the declared unit", min=0.001, default=1, precision=3)
    ec_amount_mod: FloatProperty(name="", description="EC amount modifier per declared unit", min=0.00, default=0, precision=3)
    ec_du: FloatProperty(name="", description="Embodied carbon per declared unit", default=100, precision=3)
    ec_ma: FloatProperty(name="", description="Area of object in embodied carbon terms", default=1, precision=3)
    ec_weight: FloatProperty(name="kg", description="Weight", default=1)
    ec_density: FloatProperty(name="kg/m^3", description="Material density", default=1000)
    ec_life: iprop("y", "Lifespan in years", 1, 100, 60)
    ec_rep: eprop([("0", "Volume", "Object represents a material volume"), ("1", "Item", "Object represents an item")], "", "Specify what the object represents in EC terms", "0")
    ec_arep: eprop([("0", "Manual area", "Area defined manually"), ("1", "Object area", "Area defined by all object faces")], "", "Specify what the object represents in EC terms", "0")
    ec_items: FloatProperty(name="", description="Number of items the object represents", min=0.001, default=1, precision=3)
    ec_mod: StringProperty(name="", description="Embodied modules reported")
    ee = envi_embodied()
    write_stl = ob_to_stl
    auvi_sl: eprop([("0", "Source", "Object represents a sound source"), ("1", "Listener", "Object represents an item")], "", "Specify what the object represents in AuVi terms", "0")
    

class VI_Params_Material(bpy.types.PropertyGroup):
    radtex: bprop("", "Flag to signify whether the material has a texture associated with it", False)
    amask: StringProperty(name="", description="Name of the alpha mask file", default="", subtype="FILE_PATH")
    radnorm: bprop("", "Flag to signify whether the material has a normal map associated with it", False)
    nu: fvprop(3, '', 'Image up vector', [0, 0, 1], 'XYZ', -1, 1)
    nside: fvprop(3, '', 'Image side vector', [-1, 0, 0], 'XYZ', -1, 1)
    radcolour: fvprop(3, "Material Reflectance", 'Material Reflectance', [0.8, 0.8, 0.8], 'COLOR', 0, 1)
    radcolmenu: eprop([("0", "RGB", "Specify colour temperature"), ("1", "Temperature", "Specify colour temperature")], "Colour type:", "Specify the colour input", "0")
    radrough: fprop("Roughness", "Material roughness", 0, 1, 0.1)
    radspec: fprop("Specularity", "Material specular reflection", 0, 1, 0.0)
    radtrans: fvprop(3, "Transmission", 'Material transmission', [0.8, 0.8, 0.8], 'COLOR', 0, 1)
    radtransmit: fprop("Transmittance", "Material transmittance", 0, 1, 0.9)
    radtransmenu: eprop([("0", "Transmission", "RGB transmission"), ("1", "Transmittance", "Transmittance")], "Trans type", "Specify the material transmission", "0")
    radtransdiff: fprop("Transmission", "Material diffuse transmission", 0, 1, 0.1)
    radtranspec: fprop("Trans spec", "Material specular transmission", 0, 1, 0.1)
    radior: fprop("IOR", "Material index of refractionn", 0, 5, 1.5)
    radct: iprop("Temperature (K)", "Colour temperature in Kelvin", 0, 12000, 4700)
    radintensity: fprop("Intensity", u"Material radiance (W/sr/m\u00b2)", 0, 10000, 1)
    radfile: sprop("", "Radiance file material description", 1024, "")
    vi_shadow: bprop("VI Shadow", "Flag to signify whether the material represents a VI Shadow sensing surface", False)
    livi_sense: bprop("LiVi Sensor", "Flag to signify whether the material represents a LiVi sensing surface", False)
    livi_compliance: bprop("LiVi Compliance Surface", "Flag to signify whether the material represents a LiVi compliance surface", False)
    gl_roof: bprop("Glazed Roof", "Flag to signify whether the communal area has a glazed roof", False)
    hspacetype = [('0', 'Public/Staff', 'Public/Staff area'), ('1', 'Patient', 'Patient area')]
    rspacetype = [('0', "Kitchen", "Kitchen space"), ('1', "Living/Dining/Study", "Living/Dining/Study area"), ('2', "Communal", "Non-residential or communal area")]
    respacetype = [('0', "Sales", "Sales space"), ('1', "Occupied", "Occupied space")]
    lespacetype = [('0', "Healthcare", "Healthcare space"), ('1', "Other", "Other space")]
    hspacemenu: eprop(hspacetype, "", "Type of healthcare space", '0')
    brspacemenu: eprop(rspacetype, "", "Type of residential space", '0')
    crspacemenu: eprop(rspacetype[:2], "", "Type of residential space", '0')
    respacemenu: eprop(respacetype, "", "Type of retail space", '0')
    lespacemenu: eprop(lespacetype, "", "Type of space", '0')
    BSDF: bprop("", "Flag to signify a BSDF material", False)
    mattype: eprop([("0", "Geometry", "Geometry"), ("1", 'Light sensor', "LiVi sensing material"), 
                    ("2", "FloVi boundary", 'FloVi blockmesh boundary'), ("3", "AuVi material", 'Acoustic material')], "", "VI-Suite material type", "0")
    envi_nodes: bpy.props.PointerProperty(type=bpy.types.NodeTree)
    envi_rev_enum: EnumProperty(items=ret_envi_mats, name='', description='EnVi material')
    envi_type: sprop("", "EnVi Material type", 64, "None")
    envi_shading: bprop("", "Flag to signify whether the material contains shading elements", False)
    envi_boundary: bprop("", "Flag to signify whether the material represents a zone boundary", False)
    envi_export: bprop("Material Export", "Flag to tell EnVi to export this material", False)
    pport: bprop("", "Flag to signify whether the material represents a Photon Port", False)
    radtypes = [('0', 'Plastic', 'Plastic Radiance material'), ('1', 'Glass', 'Glass Radiance material'), ('2', 'Dielectric', 'Dialectric Radiance material'),
                ('3', 'Translucent', 'Translucent Radiance material'), ('4', 'Mirror', 'Mirror Radiance material'), ('5', 'Light', 'Emission Radiance material'),
                ('6', 'Metal', 'Metal Radiance material'), ('7', 'Anti-matter', 'Antimatter Radiance material'), ('8', 'BSDF', 'BSDF Radiance material'), ('9', 'Custom', 'Custom Radiance material')]
    radmatmenu: eprop(radtypes, "", "Type of Radiance material", '0')
    radmatdict = {'0': ['radcolour', 0, 'radrough', 'radspec'], '1': ['radtransmenu', 0, 'radtrans', 0, 'radtransmit'], '2': ['radtrans', 0, 'radior'], '3': ['radcolour', 0, 'radspec', 'radrough', 0, 'radtransdiff',  'radtranspec'], '4': ['radcolour'],
    '5': ['radcolmenu', 0, 'radcolour', 0, 'radct',  0, 'radintensity'], '6': ['radcolour', 0, 'radrough', 'radspec'], '7': [], '8': [], '9': []}
    radmat = radmat
    li_bsdf_proxy_depth: fprop("", "Depth of proxy geometry", -10, 10, 0)
    li_bsdf_up: fvprop(3, '', 'BSDF up vector', [0, 0, 1], 'XYZ', -1, 1)
    li_tex: bpy.props.PointerProperty(name="", type=Image, description='Texture map')
    li_tex_black: bpy.props.PointerProperty(name="", type=Material, description='Texture map (black)')
    li_am: bpy.props.PointerProperty(name="", type=Image, description='Alpha mask')
    li_norm: bpy.props.PointerProperty(name="", type=Image, description='Normal map')
    li_norm_strength: fprop("", "Normal strength", 0, 1000, 1.0)
    li_dirt: bprop("", "Applies a noise pattern to reduce reflectivity", False)
    li_dirt_spacing: fprop("m", "Noise dimension", 0, 100, 0.1)
    li_dirt_level: fprop("", "Level of reflectivity reduction", 0, 1, 0.5)

    # FloVi Materials
    flovi_bmb_type: EnumProperty(items=ret_fvb_menu, name="", description="FloVi boundary")
    flovi_mat = fvmat
    flovi_bmbp_subtype: EnumProperty(items=ret_fvbp_menu, name="", description="FloVi boundary sub-type")
    flovi_bmbp_val: fprop("", "Pressure value", -1000, 1000, 0.0)
    flovi_p_field: bprop("", "Take boundary velocity from the field velocity", False)
    flovi_bmbp_p0val: fprop("", "Pressure value", -1000, 1000, 0)
    flovi_bmbp_gamma: fprop("", "Pressure value", -1000, 1000, 1.4)

    flovi_bmbu_subtype: EnumProperty(items=ret_fvbu_menu, name="", description="FloVi sub-type boundary")
    flovi_bmbu_val: fvprop(3, '', 'Vector value', [0, 0, 0], 'VELOCITY', -100, 100)
    flovi_bmbu_inval: fvprop(3, '', 'Vector value', [0, 0, 0], 'VELOCITY', -100, 100)
    flovi_u_type: eprop([("0", "Vector", "Specify velocity as a vector"), ("1", 'Azimuth', "Specify velocity as an azimuth and speed")], "", "Velocity type", "0")
    flovi_u_field: bprop("", "Take boundary velocity from the field velocity", False)
    flovi_u_azi: fprop("deg.", "Flow azimuth direction", 0, 360, 0)
    flovi_u_speed: fprop("m/s", "Flow speed", 0, 1000, 1)
    flovi_u_uref: fprop("m/s", "Reference speed", 0, 500, 0)
    flovi_u_zref: fprop("m", "Reference z height", 10, 500, 0)
    flovi_u_zdir: fvprop(3, '', 'Up direction', [0, 0, 1], 'XYZ', -100, 100)
    flovi_u_fdir: fvprop(3, '', 'Flow direction', [1, 0, 0], 'XYZ', -100, 100)
    flovi_u_z0: fprop("m", "	Surface roughness length", 0.1, 500, 0)
    flovi_u_zground: fprop("m", "Ground height in global Z", 0, 500, 0)
    flovi_u_d: fprop("", "Displacement value", 0, 500, 0)

    flovi_bmbnut_subtype: EnumProperty(items=ret_fvbnut_menu, name="", description="FloVi sub-type boundary")
    flovi_bmbnut_val: fprop("", "Nut value", -1000, 1000, 0.0)
    flovi_nut_field: bprop("", "Take boundary nut from the field nut", False)

    flovi_k_subtype: EnumProperty(items=ret_fvbk_menu, name="", description="FloVi k sub-type boundary")
    flovi_k_val: FloatProperty(name="", description="k value", min=0.0001, max=1, default=0.001, precision=4)
    flovi_k_intensity: fprop("", "k value", 0.001, 1000, 0.14)
    flovi_k_field: bprop("", "Take boundary k from the field k", False)

    flovi_bmbe_subtype: EnumProperty(items=ret_fvbepsilon_menu, name="", description="FloVi epsilon sub-type boundary")
    flovi_bmbe_val: FloatProperty(name="", description="Epsilon value", min=0.0001, max=1, default=0.001, precision=4)
    flovi_eml_val: FloatProperty(name="", description="Mixing length", min=0.001, max=1, default=0.005, precision=4)
    flovi_e_field: bprop("", "Take boundary epsilon from the field epsilon", False)

    flovi_bmbt_subtype: EnumProperty(items = ret_fvbt_menu, name = "", description = "FloVi sub-type boundary")
    flovi_bmbt_val: fprop("", "T value", 0, 1000, 300)
    flovi_bmbti_val: fprop("", "T inlet/outlet value", 0, 1000, 300)
    flovi_t_field: bprop("", "Take boundary t from the field t", False)

    flovi_a_subtype: EnumProperty(items = ret_fvba_menu, name = "", description = "FloVi sub-type boundary")
    flovi_a_val: fprop("", "T value", -1000, 1000, 0.0)
    flovi_a_field: bprop("", "Take boundary alphat from the field alphat", False)

    flovi_prgh_subtype: EnumProperty(items = ret_fvbprgh_menu, name = "", description = "FloVi sub-type boundary")
    flovi_prgh_val: fprop("", "p_rgh value", -100000, 100000, 0.0)
    flovi_prgh_p: fprop("", "p_rgh p", -100000, 100000, 0.0)
    flovi_prgh_field: bprop("", "Take boundary p_rgh from the field p_rgh", True)
    flovi_prgh_p0: fprop("", "p_rgh p0", -100000, 100000, 0.0)
    flovi_prgh_gamma: fprop("", "p_rgh gamma", 0, 10, 1.4)

    flovi_ng_max: fprop("", "Netgen max cell size (face defined)", 0.001, 100, 0.5)
    flovi_ng_emax: fprop("", "Netgen max cell size (edge defined)", 0.001, 1, 1)
    flovi_rad_subtype: EnumProperty(items=ret_fvrad_menu, name="", description="FloVi sub-type boundary")
    flovi_rad_em: eprop([('lookup', 'Lookup', 'Lookup emissivity')], "", "Emissivity mode", 'lookup')
    flovi_rad_e: fprop("", "Emissivity value", 0, 1, 0.5)
    flovi_rad_val: fprop("", "Radiation value", 0, 10000, 0)

    flovi_i_subtype: EnumProperty(items=ret_fvi_menu, name="", description="FloVi sub-type boundary")
    flovi_i_em: eprop([('lookup', 'Lookup', 'Lookup emissivity')], "", "Emissivity mode", 'lookup')
    flovi_i_e: fprop("", "Emissivity value", 0, 1, 0.5)
    flovi_i_val: fprop("", "Radiation value", 0, 10000, 0)
    flovi_probe: bprop("", "Turn on pressure monitoring", False)
    flovi_htc: bprop("", "Turn on heat transfer coefficient calculation", False)

    # AuVi materials
    auvi_abs_class: eprop([('0', 'Database', 'Database material'), ('1', 'Custom', 'Custom material')], "", "AuVi material class", '0')
    auvi_scatt_class: eprop([('0', 'Database', 'Database material'), ('1', 'Custom', 'Custom material')], "", "AuVi material class", '1')
    # auvi_scatter: bprop("", "Turn on scattering for ray tracing", False)
    auvi_type_abs: EnumProperty(items=auvi_type_abs, name="", description="Layer embodied material type")
    auvi_mat_abs: EnumProperty(items=auvi_abs, name="", description="Layer embodied material type")
    auvi_type_scatt: EnumProperty(items=auvi_type_scatt, name="", description="Layer embodied material type")
    auvi_mat_scatt: EnumProperty(items=auvi_scatt, name="", description="Layer embodied material type")
    auvi_o1_abs: fprop("", "125 Hz absorption value", 0, 1, 0.1)
    auvi_o1_scatt: fprop("", "125 Hz scatter value", 0, 1, 0.05)
    auvi_o2_abs: fprop("", "250 Hz absorption value", 0, 1, 0.1)
    auvi_o2_scatt: fprop("", "250 Hz scatter value", 0, 1, 0.05)
    auvi_o3_abs: fprop("", "500 Hz absorption value", 0, 1, 0.1)
    auvi_o3_scatt: fprop("", "500 Hz scatter value", 0, 1, 0.05)
    auvi_o4_abs: fprop("", "1000 Hz absorption value", 0, 1, 0.1)
    auvi_o4_scatt: fprop("", "1000 Hz scatter value", 0, 1, 0.05)
    auvi_o5_abs: fprop("", "2000 Hz absorption value", 0, 1, 0.1)
    auvi_o5_scatt: fprop("", "2000 Hz scatter value", 0, 1, 0.05)
    auvi_o6_abs: fprop("", "4000 Hz absorption value", 0, 1, 0.1)
    auvi_o6_scatt: fprop("", "4000 Hz scatter value", 0, 1, 0.05)
    auvi_o7_abs: fprop("", "8000 Hz absorption value", 0, 1, 0.1)
    auvi_o7_scatt: fprop("", "8000 Hz scatter value", 0, 1, 0.05)
    auvi_abs_flat: bprop("", "Turns on equal absorbtion for all frequency bands", True)
    auvi_scatt_flat: bprop("", "Turns on equal scattering for all frequency bands", True)
    am = auvi_materials()

class VI_Params_Collection(bpy.types.PropertyGroup):
    envi_collection: bprop("", "Flag to tell EnVi to export this collection", False)
    envi_zone: bprop("", "Flag to tell EnVi this is a geometry collection", False)
    envi_shade: bprop("", "Flag to tell EnVi this is a shader collection", False)
    embodied: BoolProperty(name="", description="Embodied carbon", default=0)
    embodiedtype: EnumProperty(items=envi_elayertype, name="", description="Layer embodied material class", update=ec_update)
    embodiedclass: EnumProperty(items=envi_eclasstype, name="", description="Layer embodied class", update=ec_update)
    embodiedmat: EnumProperty(items=envi_emattype, name="", description="Layer embodied material", update=ec_update)
    ec_id: StringProperty(name="", description="Embodied id (unique indentifier")
    ec_type: StringProperty(name="", description="Embodied type e.g. Insulation")
    ec_class: StringProperty(name="", description="Embodied class (class of type e.g. phenolic foam)")
    ec_name: StringProperty(name="", description="Embodied name")
    ec_unit: EnumProperty(items=[("kg", "kg", "per kilogram"),
                                  ("m2", "m2", "per square metre"),
                                  ("m3", "m3", "per cubic metre"),
                                  ("each", "Each", "per item")],
                                  name="",
                                  description="Embodied carbon unit",
                                  default="kg")
    ec_amount: FloatProperty(name="", description="Amount of the declared unit", min=0.001, default=1, precision=3)
    ec_amount_mod: FloatProperty(name="", description="EC amount modifier per amount of DU", min=0.00, default=0, precision=3)
    ec_du: FloatProperty(name="", description="Embodied carbon per declared unit", default=100)
    ec_weight: FloatProperty(name="kg", description="Weight", default=1)
    ec_density: FloatProperty(name="kg/m^3", description="Material density", default=1000)
    ec_life: iprop("y", "Lifespan in years", 1, 100, 60)
    ec_mod: StringProperty(name="", description="Embodied modules reported")
    ee = envi_embodied()

class VI_Params_Link(bpy.types.PropertyGroup):
    vi_uid: iprop("ID", "Unique ID", 0, 10000, 0)

@persistent
def update_chart_node(dummy):
    try:
        for ng in [ng for ng in bpy.data.node_groups if ng.bl_idname == 'ViN']:
            [node.update() for node in ng.nodes if node.bl_label == 'VI Chart']
    except Exception as e:
        print('Chart node cannot update:', e)

@persistent
def update_dir(dummy):
    if bpy.context.scene.vi_params.get('viparams'):
        fp = bpy.data.filepath
        bpy.context.scene.vi_params['viparams']['newdir'] = os.path.join(os.path.dirname(fp), os.path.splitext(os.path.basename(fp))[0])

@persistent
def display_off(dummy):
    if bpy.context.scene.vi_params.get('viparams') and bpy.context.scene.vi_params['viparams'].get('vidisp'):
        ifdict = {'sspanel': 'ss', 'lipanel': 'li', 'enpanel': 'en', 'bsdf_panel': 'bsdf'}

        if bpy.context.scene.vi_params['viparams']['vidisp'] in ifdict:
            bpy.context.scene.vi_params['viparams']['vidisp'] = ifdict[bpy.context.scene.vi_params['viparams']['vidisp']]

        bpy.context.scene.vi_params.vi_display = 0

@persistent
def display_off_load(dummy):
    if bpy.context.scene.vi_params.get('vi_display'):
        bpy.context.scene.vi_params.vi_display = 0

bpy.app.handlers.load_post.append(display_off_load)

@persistent
def select_nodetree(dummy):
    for space in getViEditorSpaces():
        vings = [ng for ng in bpy.data.node_groups if ng.bl_idname == 'ViN']
        if vings:
            space.node_tree = vings[0]
            space.node_tree.use_fake_user = 1

    for space in getEnViEditorSpaces():
        envings = [ng for ng in bpy.data.node_groups if ng.bl_idname == 'EnViN']
        if envings:
            space.node_tree = envings[0]

    for space in getEnViMaterialSpaces():
        try:
            if space.node_tree != bpy.context.active_object.active_material.vi_params.envi_nodes:
                envings = [ng for ng in bpy.data.node_groups if ng.bl_idname == 'EnViMatN' and ng == bpy.context.active_object.active_material.vi_params.envi_nodes]
                if envings:
                    space.node_tree = envings[0]
        except Exception as nt_error:
            pass

bpy.app.handlers.depsgraph_update_post.append(select_nodetree)

def getViEditorSpaces():
    if bpy.context.screen:
        return [area.spaces.active for area in bpy.context.screen.areas if area and area.type == "NODE_EDITOR" and area.spaces.active.tree_type == "ViN" and not area.spaces.active.edit_tree]
    else:
        return []

def getEnViEditorSpaces():
    if bpy.context.screen:
        return [area.spaces.active for area in bpy.context.screen.areas if area and area.type == "NODE_EDITOR" and area.spaces.active.tree_type == "EnViN" and not area.spaces.active.edit_tree]
    else:
        return []

def getEnViMaterialSpaces():
    if bpy.context.screen:
        return [area.spaces.active for area in bpy.context.screen.areas if area and area.type == "NODE_EDITOR" and area.spaces.active.tree_type == "EnViMatN"]
    else:
        return []

def path_update():
    vi_prefs = bpy.context.preferences.addons[__name__].preferences
    epdir = vi_prefs.epbin if vi_prefs and vi_prefs.epbin and os.path.isdir(vi_prefs.epbin) else os.path.join('{}'.format(addonpath), 'EPFiles', str(sys.platform))
    radldir = vi_prefs.radlib if vi_prefs and os.path.isdir(vi_prefs.radlib) else os.path.join('{}'.format(addonpath), 'RadFiles', 'lib')
    radbdir = vi_prefs.radbin if vi_prefs and os.path.isdir(vi_prefs.radbin) else os.path.join('{}'.format(addonpath), 'RadFiles', str(sys.platform), 'bin')
    ofbdir = os.path.abspath(vi_prefs.ofbin) if vi_prefs and os.path.isdir(vi_prefs.ofbin) else os.path.join('{}'.format(addonpath), 'OFFiles', str(sys.platform), 'bin')
    # datadir = os.path.abspath(vi_prefs.datab) if vi_prefs and os.path.isdir(vi_prefs.datab) else os.path.join('{}'.format(addonpath), 'EPFiles')

    if not os.environ.get('RAYPATH') or radldir not in os.environ['RAYPATH'] or radbdir not in os.environ['PATH'] or epdir not in os.environ['PATH'] or ofbdir not in os.environ['PATH']:
        if vi_prefs and os.path.isdir(vi_prefs.radlib):
            os.environ["RAYPATH"] = '{0}{1}{2}'.format(radldir, os.pathsep, os.path.join(addonpath, 'RadFiles', 'lib'))
        else:
            os.environ["RAYPATH"] = radldir

        if radbdir not in os.environ["PATH"]:
            native_path = os.path.join('{}'.format(addonpath), 'RadFiles', str(sys.platform), 'bin')
            if native_path in os.environ["PATH"]:
                os.environ["PATH"].replace(native_path, radbdir)
            else:
                os.environ["PATH"] += '{0}{1}'.format(os.pathsep, radbdir)

        os.environ["PATH"] += "{0}{1}{0}{2}{0}{3}".format(os.pathsep, epdir, ofbdir, os.path.join('{}'.format(addonpath), 'Python', str(sys.platform), 'bin'))
        sys.path.append(ofbdir)

classes = (VIPreferences, ViNetwork, No_Loc, So_Vi_Loc, No_Vi_SP, NODE_OT_SunPath, NODE_OT_TextUpdate, NODE_OT_FileSelect, NODE_OT_HdrSelect,
           VI_PT_3D, VI_Params_Scene, VI_Params_Object, VI_Params_Material, VI_Params_Collection, No_Vi_WR, No_Vi_SVF, NODE_OT_WindRose, VIEW3D_OT_WRDisplay,
           NODE_OT_SVF, So_Vi_Res, VI_PT_Mat, VIEW3D_OT_SVFDisplay, MAT_EnVi_Node, No_Vi_SS, NODE_OT_Shadow, VIEW3D_OT_SSDisplay, VIEW3D_OT_RTDisplay,
           No_Li_Geo, No_Li_Con, No_Li_Sen, So_Li_Geo, NODE_OT_Li_Geo, So_Li_Con, NODE_OT_Li_Con, No_Text, So_Text,
           No_Vi_Im, No_Li_Im, So_Li_Im, NODE_OT_Li_Im, NODE_OT_Li_Pre, No_Li_Sim, NODE_OT_Li_Sim, VIEW3D_OT_Li_BD,
           No_Li_Gl, No_Li_Fc, NODE_OT_Li_Gl, NODE_OT_Li_Fc, No_En_Geo, VI_PT_Ob, NODE_OT_En_Geo, EnViNetwork, No_En_Net_Zone,
           EnViMatNetwork, No_En_Mat_Con, VI_PT_Gridify, OBJECT_OT_VIGridify, No_En_Mat_Sc, No_En_Mat_Sh, No_En_Mat_ShC, No_En_Mat_Bl, No_En_Mat_SG,
           NODE_OT_En_UV, NODE_OT_En_EC, No_En_Net_Occ, So_En_Net_Occ, So_En_Net_Sched, So_En_Mat_Sched, So_En_Net_Inf, So_En_Net_Hvac, So_En_Net_Eq,
           No_En_Mat_Op, No_En_Mat_Tr, So_En_Mat_Ou, So_En_Mat_Fr, So_En_Mat_Op, So_En_Mat_Tr, So_En_Mat_Gas, No_En_Con,
           So_En_Mat_Sh, So_En_Mat_ShC, So_En_Mat_Sc, No_Vi_EC, NODE_OT_EC, OBJECT_OT_EcS, OBJECT_OT_EcE, NODE_OT_ECPie, NODE_OT_WLCLine, NODE_OT_COMLine,
           So_En_Con, So_En_Geo, NODE_OT_En_Con, No_En_Sim, NODE_OT_En_Sim, No_En_Mat_Gas,
           No_Vi_Chart, No_Vi_HMChart, So_En_Res, So_En_ResU, NODE_OT_Chart, NODE_OT_HMChart, No_En_Net_Hvac, So_En_Net_TSched, No_En_Net_Eq, No_En_Net_Sched, No_En_Net_Inf,
           No_En_Net_SFlow, No_En_Net_SSFlow, So_En_Net_SFlow, So_En_Net_SSFlow, So_En_Mat_PV, No_En_Mat_PV, No_En_Mat_Sched,
           So_En_Mat_PVG, No_En_Mat_PVG, NODE_OT_En_PVA, No_Vi_Metrics, NODE_OT_En_PVS, NODE_OT_En_LayS, NODE_OT_En_EcS, NODE_OT_En_ConS, So_En_Net_Bound,
           No_En_Net_ACon, No_En_Net_Ext, No_En_Net_EMSZone, No_En_Net_Prog, No_En_Net_EMSPy, So_En_Net_Act, So_En_Net_Sense,
           TREE_PT_vi, TREE_PT_envin, TREE_PT_envim,  TREE_OT_goto_mat, TREE_OT_goto_group,
           OBJECT_OT_Li_GBSDF, MATERIAL_OT_Li_LBSDF, MATERIAL_OT_Li_SBSDF, OBJECT_OT_GOct, OBJECT_OT_Embod, MATERIAL_OT_Li_DBSDF, VIEW3D_OT_Li_DBSDF, NODE_OT_CSV, No_CSV,
           NODE_OT_ASCImport, No_ASC_Import, So_Flo_Mesh, No_Flo_Case, So_Flo_Case, NODE_OT_Flo_Case, No_Flo_NG, NODE_OT_Flo_NG,
           So_Flo_Con, No_Flo_Bound, NODE_OT_Flo_Bound, No_Flo_Sim, NODE_OT_Flo_Sim, No_En_IF, No_En_RF, So_En_Net_WPC, No_En_Net_Azi, MAT_EnVi_Node_Remove, No_Anim, So_Anim,
           No_En_Net_Anim, No_En_Mat_Anim, VI_PT_Col, NODE_OT_Vi_Info, ViEnRIn, NODE_OT_EcE, So_Au_Scene, So_Au_IR, No_Au_Sim, No_Au_Conv, NODE_OT_Au_Conv,
           NODE_OT_Au_Rir, NODE_OT_WavSelect, NODE_OT_Au_Play, NODE_OT_Au_Stop, NODE_OT_Au_PlayC, NODE_OT_Au_Save)


def register():
    for cl in classes:
        bpy.utils.register_class(cl)

    Object, Scene, Material, Collection, Link = bpy.types.Object, bpy.types.Scene, bpy.types.Material, bpy.types.Collection, bpy.types.NodeLink
    Scene.vi_params = bpy.props.PointerProperty(type=VI_Params_Scene)
    Object.vi_params = bpy.props.PointerProperty(type=VI_Params_Object)
    Material.vi_params = bpy.props.PointerProperty(type=VI_Params_Material)
    Collection.vi_params = bpy.props.PointerProperty(type=VI_Params_Collection)
    Link.vi_params = bpy.props.PointerProperty(type=VI_Params_Link)
    nodeitems_utils.register_node_categories("Vi Nodes", vinode_categories)
    nodeitems_utils.register_node_categories("EnVi Nodes", envinode_categories)
    nodeitems_utils.register_node_categories("EnVi Material Nodes", envimatnode_categories)

    if update_chart_node not in bpy.app.handlers.load_post:
        bpy.app.handlers.load_post.append(update_chart_node)

    if display_off not in bpy.app.handlers.load_pre:
        bpy.app.handlers.load_pre.append(display_off)

    if update_dir not in bpy.app.handlers.load_post:
        bpy.app.handlers.load_post.append(update_dir)

    path_update()

def unregister():
    nodeitems_utils.unregister_node_categories("EnVi Material Nodes")
    nodeitems_utils.unregister_node_categories("EnVi Nodes")
    nodeitems_utils.unregister_node_categories("Vi Nodes")

    for cl in reversed(classes):
        try:
            bpy.utils.unregister_class(cl)
        except Exception:
            pass

    if update_chart_node in bpy.app.handlers.load_post:
        bpy.app.handlers.load_post.remove(update_chart_node)

    if display_off in bpy.app.handlers.load_post:
        bpy.app.handlers.load_post.remove(display_off)

    if update_dir in bpy.app.handlers.load_post:
        bpy.app.handlers.load_post.remove(update_dir)
