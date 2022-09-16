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
#

# EnVi materials databases
import os, json
from collections import OrderedDict
from .envi_func import epentry


class envi_materials(object):
    '''Defines materials with a comma separated dictionary, with material name as key, giving
    (Roughness, Conductivity {W/m-K}, Density {kg/m3}, Specific Heat {J/kg-K}, Thermal Absorbtance,
    Solar Absorbtance, Visible Absorbtance, Default thickness)'''

    def __init__(self):
        self.update()

    def update(self):
        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'EPFiles', '{}'.format('Material_database.json')), 'r') as mat_jfile:
            mat_dict = json.loads(mat_jfile.read())

        self.metal_datd = mat_dict['Metal']
        self.brick_datd = mat_dict['Brick']
        self.cladding_datd = mat_dict['Cladding']
        self.concrete_datd = mat_dict['Concrete']
        self.wood_datd = mat_dict['Wood']
        self.stone_datd = mat_dict['Stone']
        self.gas_datd = mat_dict['Gas']
        self.wgas_datd = mat_dict['WGas']
        self.glass_datd = mat_dict['Glass']
        self.insulation_datd = mat_dict['Insulation']
        self.pcm_datd = mat_dict['PCM']
        self.pcmd_datd = mat_dict['PCMD']
        self.plastic_datd = mat_dict['Plastic']
        self.metal_dat = OrderedDict(sorted(self.metal_datd.items()))
        self.brick_dat = OrderedDict(sorted(self.brick_datd.items()))
        self.cladding_dat = OrderedDict(sorted(self.cladding_datd.items()))
        self.concrete_dat = OrderedDict(sorted(self.concrete_datd.items()))
        self.wood_dat = OrderedDict(sorted(self.wood_datd.items()))
        self.stone_dat = OrderedDict(sorted(self.stone_datd.items()))
        self.gas_dat = OrderedDict(sorted(self.gas_datd.items()))
        self.wgas_dat = OrderedDict(sorted(self.wgas_datd.items()))
        self.glass_dat = OrderedDict(sorted(self.glass_datd.items()))
        self.insulation_dat = OrderedDict(sorted(self.insulation_datd.items()))
        self.pcm_dat = OrderedDict(sorted(self.pcm_datd.items()))
        self.pcmd_dat = OrderedDict(sorted(self.pcmd_datd.items()))
        self.plastic_dat = OrderedDict(sorted(self.plastic_datd.items()))
        self.pv_datd = {'Default PV': ('Rough', '0.035', '29', '1213', '0.90', '0.5', '0.5', '50')}
        self.pv_dat = OrderedDict(sorted(self.pv_datd.items()))

        self.namedict = OrderedDict()
        self.thickdict = OrderedDict()
        self.i = 0
        self.matdat = OrderedDict()

        for dat in (self.brick_dat, self.cladding_dat, self.concrete_dat, self.gas_dat, self.insulation_dat, self.metal_dat,
                    self.stone_dat, self.wood_dat, self.glass_dat, self.wgas_dat, self.pcm_dat, self.pv_dat, self.plastic_dat):
            self.matdat.update(dat)

    def get_dat(self, mat_type):
        mat_dict = {'Glass': self.glass_datd, '3': self.metal_datd, '0': self.brick_datd, '1': self.cladding_datd,
                   '2': self.concrete_datd, '5': self.wood_datd, '4': self.stone_datd, '6': self.gas_datd,
                   'WGas': self.wgas_datd, '7': self.insulation_datd, '8': self.pcm_datd, '9': self.pcmd_datd, '10': self.plastic_datd}
        return mat_dict[mat_type]

    def lay_save(self):
        mat_dict = {'Glass': self.glass_datd, 'Metal': self.metal_datd, 'Brick': self.brick_datd, 'Cladding': self.cladding_datd,
                   'Concrete': self.concrete_datd, 'Wood': self.wood_datd, 'Stone': self.stone_datd, 'Gas': self.gas_datd,
                   'WGas': self.wgas_datd, 'Insulation': self.insulation_datd, 'PCM': self.pcm_datd, 'PCMD': self.pcmd_datd, 'Plastic': self.plastic_datd}
        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'EPFiles', '{}'.format('Material_database.json')), 'w') as mat_jfile:
            mat_jfile.write(json.dumps(mat_dict))
        self.update()

    def omat_write(self, idf_file, name, stringmat, thickness):
        params = ('Name', 'Roughness', 'Thickness (m)', 'Conductivity (W/m-K)', 'Density (kg/m3)', 'Specific Heat Capacity (J/kg-K)', 'Thermal Absorptance', 'Solar Absorptance', 'Visible Absorptance')
        paramvs = [name, stringmat[0], thickness] + stringmat[1:8]
        idf_file.write(epentry("Material", params, paramvs))

    def amat_write(self, idf_file, name, stringmat):
        params = ('Name', 'Resistance')
        paramvs = (name, stringmat[0])
        idf_file.write(epentry("Material:AirGap", params, paramvs))

    def tmat_write(self, idf_file, name, stringmat, thickness):
        params = ('Name', 'Optical Data Type', 'Window Glass Spectral Data Set Name', 'Thickness (m)', 'Solar Transmittance at Normal Incidence', 'Front Side Solar Reflectance at Normal Incidence',
                  'Back Side Solar Reflectance at Normal Incidence', 'Visible Transmittance at Normal Incidence', 'Front Side Visible Reflectance at Normal Incidence', 'Back Side Visible Reflectance at Normal Incidence',
                  'Infrared Transmittance at Normal Incidence', 'Front Side Infrared Hemispherical Emissivity', 'Back Side Infrared Hemispherical Emissivity', 'Conductivity (W/m-K)',
                  'Dirt Correction Factor for Solar and Visible Transmittance', 'Solar Diffusing')
        paramvs = [name] + stringmat[1:3] + [thickness] + ['{:.3f}'.format(float(sm)) for sm in stringmat[4:-1]] + [1, ('No', 'Yes')[stringmat[-1]]]
        idf_file.write(epentry("WindowMaterial:{}".format(stringmat[0]), params, paramvs))

    def gmat_write(self, idf_file, name, stringmat, thickness):
        params = ('Name', 'Gas Type', 'Thickness')
        paramvs = [name] + [stringmat[1]] + [thickness]
        idf_file.write(epentry("WindowMaterial:Gas", params, paramvs))

    def pcmmat_write(self, idf_file, name, stringmat):
        params = ('Name', 'Temperature Coefficient for Thermal Conductivity (W/m-K2)')
        paramvs = (name, stringmat[0])
        for i, te in enumerate(stringmat[1].split()):
            params += ('Temperature {} (C)'.format(i), 'Enthalpy {} (J/kg)'.format(i))
            paramvs +=(te.split(':')[0], te.split(':')[1])
        idf_file.write(epentry("MaterialProperty:PhaseChange", params, paramvs))

    def sg_write(self, idf_file, name, uv, shgc, vt):
        params = ('Name', 'U-Factor', 'Solar Heat Gain Coefficient', 'Visible Transmittance')
        paramvs = [name] + ['{:.3f}'.format(p) for p in (uv, shgc, vt)]
        idf_file.write(epentry("WindowMaterial:SimpleGlazingSystem", params, paramvs))

class envi_constructions(object):
    def __init__(self):
        self.update()

    def update(self):
        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'EPFiles', '{}'.format('Construction_database.json')), 'r') as con_jfile:
            con_dict = json.loads(con_jfile.read())

        self.wall_cond = con_dict['Wall']
        self.iwall_cond = con_dict['Party wall']
        self.ceil_cond = con_dict['Ceiling']
        self.ifloor_cond = con_dict['Internal floor']
        self.floor_cond = con_dict['Floor']
        self.roof_cond = con_dict['Roof']
        self.door_cond = con_dict['Door']
        self.glaze_cond = con_dict['Glazing']
        self.wall_con = OrderedDict(sorted(self.wall_cond.items()))
        self.iwall_con = OrderedDict(sorted(self.iwall_cond.items()))
        self.ceil_con = OrderedDict(sorted(self.ceil_cond.items()))
        self.ifloor_con = OrderedDict(sorted(self.ifloor_cond.items()))
        self.floor_con = OrderedDict(sorted(self.floor_cond.items()))
        self.roof_con = OrderedDict(sorted(self.roof_cond.items()))
        self.door_con = OrderedDict(sorted(self.door_cond.items()))
        self.pv_cond = {'Simple PV': ['Default PV']}
        self.pv_con = OrderedDict(sorted(self.pv_cond.items()))
        self.glaze_con = OrderedDict(sorted(self.glaze_cond.items()))
        self.p = 0
        self.propdict = {'Wall': self.wall_con, 'Floor': self.floor_con, 'Roof': self.roof_con, 'Ceiling': self.ceil_con, 'Door': self.door_con,
                                'Window': self.glaze_con, 'PV': self.pv_con, 'Internal floor': self.ifloor_con, 'Internal wall': self.iwall_con}

    def con_write(self, idf_file, contype, name, nl, mn, cln):
        params = ['Name', 'Outside layer'] + ['Layer {}'.format(i + 1) for i in range(len(cln) - 1)]
        paramvs = [mn] + cln#'{}-{}'.format(con[name][0], nl)]
        idf_file.write(epentry('Construction', params, paramvs))

    def get_dat(self, mat_type):
        mat_dict = {'Wall - External': self.wall_cond, 'Wall - Zone': self.iwall_cond, 'Wall - Adiabatic': self.iwall_cond,
                    'Floor - Ground': self.floor_cond, 'Floor - Zone': self.ifloor_cond, 'Roof - External': self.roof_cond,
                    'Roof - Zone': self.ceil_cond,
                    'Door - External': self.door_cond, 'Door - Zone': self.door_cond, 'Window - External': self.glaze_cond,
                    'Window - Zone': self.glaze_cond}
        return mat_dict[mat_type]

    def con_save(self):
        con_dict = {'Wall': self.wall_cond, 'Party wall': self.iwall_cond,
                    'Internal floor': self.ifloor_cond, 'Floor': self.floor_cond, 'Ceiling': self.ceil_cond,
                    'Door': self.door_cond, 'Glazing': self.glaze_cond, 'Roof': self.roof_cond}
        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'EPFiles', '{}'.format('Construction_database.json')), 'w') as con_jfile:
            con_jfile.write(json.dumps(con_dict))
        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'EPFiles', '{}'.format('Construction_database_backup.json')), 'w') as con_bujfile:
            con_bujfile.write(json.dumps(con_dict))
        self.update()


def retmatdict(ect, t, l):
    print(ect, t, l)
    if ect in ('Wall', 'Roof', 'Floor', 'Door', 'Ceiling', 'Frame'):
        typelist = [("0", "Brick", "Choose a material from the brick database"),("1", "Cladding", "Choose a material from the cladding database"),
                    ("2", "Concrete", "Choose a material from the concrete database"),("3", "Metal", "Choose a material from the metal database"),
                    ("4", "Stone", "Choose a material from the stone database"),("5", "Wood", "Choose a material from the wood database"),
                    ("6", "Gas", "Choose a material from the gas database"),("7", "Insulation", "Choose a material from the insulation database"),
                    ("8", "PCM", "Choose a material from the phase change database"), ("9", "PV", "Choose a material from the photovoltaic database"),
                    ("10", "Plastic", "Choose a material from the plastic database")]
        matdict = {'0': envi_materials().brick_dat.keys(), '1': envi_materials().cladding_dat.keys(), '2': envi_materials().concrete_dat.keys(),
                   '3': envi_materials().metal_dat.keys(), '4': envi_materials().stone_dat.keys(), '5': envi_materials().wood_dat.keys(),
                   '6': envi_materials().gas_dat.keys(), '7': envi_materials().insulation_dat.keys(), '8': envi_materials().pcm_dat.keys(),
                   '9': envi_materials().pv_dat.keys(), '10': envi_materials().plastic_dat.keys()}

    elif ect == 'Window':
        if not l % 2:
            typelist = [("0", "Glass", "Choose a material from the glass database")]
            matdict = {'0': envi_materials().glass_dat.keys()}
        else:
            typelist = [("0", "Gas", "Choose a material from the gas database")]
            matdict = {'0': envi_materials().wgas_dat.keys()}
    else:
        typelist = [('0', 'None', 'None')]
        matdict = {'0': ['None']}
    if t:
        return typelist
    else:
        return matdict


def envi_con_list(self, context):
    ec = envi_constructions()
    return [(mat, mat, 'Construction') for mat in ((ec.wall_con, ec.iwall_con)[self.envi_con_con in ("Zone", "Thermal mass")], (ec.roof_con, ec.ceil_con)[self.envi_con_con in ("Zone", "Thermal mass")], (ec.floor_con, ec.ifloor_con)[self.envi_con_con in ("Zone", "Thermal mass")], ec.door_con, ec.glaze_con, ec.pv_con)[("Wall", "Roof", "Floor", "Door", "Window", "PV").index(self.envi_con_type)]]


def retuval(mat):
    if mat.envi_con_type not in ('None', 'Shading', 'Aperture', 'Window'):
        resists, em, ec = [], envi_materials(), envi_constructions()
        thicks = [0.001 * tc for tc in (mat.envi_export_lo_thi, mat.envi_export_l1_thi, mat.envi_export_l2_thi, mat.envi_export_l3_thi, mat.envi_export_l4_thi)]
        laymats = (mat.envi_material_lo, mat.envi_material_l1, mat.envi_material_l2, mat.envi_material_l3, mat.envi_material_l4)
        pstcs = []

        if mat.envi_con_makeup == '1':
            lays = (mat.envi_layero, mat.envi_layer1, mat.envi_layer2, mat.envi_layer3, mat.envi_layer4)
            ctcs = (mat.envi_export_lo_tc, mat.envi_export_l1_tc, mat.envi_export_l2_tc, mat.envi_export_l3_tc, mat.envi_export_l4_tc)

            for l, lay in enumerate(lays):
                if lay == '1':
                    dtc = em.matdat[laymats[l]][2] if em.matdat[laymats[l]][0] == 'Gas' else em.matdat[laymats[l]][1]
                    resists.append((thicks[l]/float(dtc), float(dtc))[em.matdat[laymats[l]][0] == 'Gas'])
                if lay == '2':
                    resists.append(thicks[l]/ctcs[l])

        elif mat.envi_con_makeup == '0':
            for p, psmat in enumerate(ec.propdict[mat.envi_con_type][mat.envi_con_list]):
                pi = 2 if psmat in em.gas_dat else 1
                pstcs.append(float(em.matdat[psmat][pi]))
                resists.append((thicks[p]/float(em.matdat[psmat][pi]), float(em.matdat[psmat][pi]))[em.matdat[psmat][0] == 'Gas'])
        uv = 1/(sum(resists) + 0.12 + 0.08)
        mat.envi_material_uv = '{:.3f}'.format(uv)
        return uv
    else:
        return 1.0


def envi_layertype(self, context):
    return retmatdict(self.envi_con_type, 1, self.bl_idname == 'No_En_Mat_Gas')


def envi_elayertype(self, context):
    return [(k, k, '{} type'.format(k)) for k in envi_embodied().propdict.keys()]


def envi_layer(self, context):
    if self.materialtype:
        return [((mat, mat, 'Layer material')) for mat in list(retmatdict(self.envi_con_type, 0, self.bl_idname == 'No_En_Mat_Gas')[self.materialtype])]
    else:
        return [('', '', '')]


def envi_eclasstype(self, context):
    if self.embodiedtype:
        return [((mat, mat, 'Embodied class material')) for mat in envi_embodied().propdict[self.embodiedtype]]
    else:
        return [('', '', '')]


def envi_emattype(self, context):
    if self.embodiedtype:
        return [((mat, mat, 'Embodied material')) for mat in envi_embodied().propdict[self.embodiedtype][self.embodiedclass]]
    else:
        return [('', '', '')]


class envi_embodied(object):
    '''Defines materials with a comma separated dictionary, with material name as key, giving
    embodied energy metrics from ICE databsae v3.0'''

    def __init__(self):
        self.update()

    def update(self):
        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'EPFiles', '{}'.format('EC_database.json')), 'r') as ec_jfile:
            ec_dict = json.loads(ec_jfile.read())
            self.agg_ecd = ec_dict['aggregatesand']
            self.alu_ecd = ec_dict['aluminium']
            self.asp_ecd = ec_dict['asphalt']
            self.bit_ecd = ec_dict['bitumen']
            self.cement_ecd = ec_dict['cement']
            self.clay_ecd = ec_dict['clay']
            self.con_ecd = ec_dict['concrete']
            self.glass_ecd = ec_dict['glass']
            self.steel_ecd = ec_dict['steel']
            self.timber_ecd = ec_dict['timber']
            self.agg_ec = OrderedDict(sorted(self.agg_ecd.items()))
            self.alu_ec = OrderedDict(sorted(self.alu_ecd.items()))
            self.asp_ec = OrderedDict(sorted(self.asp_ecd.items()))
            self.bit_ec = OrderedDict(sorted(self.bit_ecd.items()))
            self.cement_ec = OrderedDict(sorted(self.cement_ecd.items()))
            self.clay_ec = OrderedDict(sorted(self.clay_ecd.items()))
            self.con_ec = OrderedDict(sorted(self.con_ecd.items()))
            self.glass_ec = OrderedDict(sorted(self.glass_ecd.items()))
            self.steel_ec = OrderedDict(sorted(self.steel_ecd.items()))
            self.timber_ec = OrderedDict(sorted(self.timber_ecd.items()))
            self.propdict = {'Aggregate': self.agg_ec, 'Aluminium': self.alu_ec, 'Asphalt': self.asp_ec, 'Bitumen': self.bit_ec, 'Cement': self.cement_ec,
                             'Clay': self.clay_ec, 'Concrete': self.con_ec, 'Glass': self.glass_ec, 'Steel': self.steel_ec, 'Timber': self.timber_ec}


#            print(ec_dict.keys())


#def retmatdict(self, t, l):
#    if self.envi_con_type in ('Wall', 'Roof', 'Floor', 'Door', 'Ceiling'):
#        typelist = [("0", "Brick", "Choose a material from the brick database"),("1", "Cladding", "Choose a material from the cladding database"), ("2", "Concrete", "Choose a material from the concrete database"),("3", "Metal", "Choose a material from the metal database"),
#                   ("4", "Stone", "Choose a material from the stone database"),("5", "Wood", "Choose a material from the wood database"),
#                   ("6", "Gas", "Choose a material from the gas database"),("7", "Insulation", "Choose a material from the insulation database"),
#                    ("8", "PCM", "Choose a material from the phase change database")]
#        matdict = {'0': envi_materials().brick_dat.keys(), '1': envi_materials().cladding_dat.keys(), '2': envi_materials().concrete_dat.keys(), '3': envi_materials().metal_dat.keys(), '4': envi_materials().stone_dat.keys(),
#                   '5': envi_materials().wood_dat.keys(), '6': envi_materials().gas_dat.keys(), '7': envi_materials().insulation_dat.keys(), '8': envi_materials().pcm_dat.keys()}
#
#    elif self.envi_con_type == 'Window':
#        if not l % 2:
#            typelist = [("0", "Glass", "Choose a material from the glass database")]
#            matdict = {'0': envi_materials().glass_dat.keys()}
#        else:
#            typelist = [("0", "Gas", "Choose a material from the gas database")]
#            matdict = {'0': envi_materials().wgas_dat.keys()}
#    else:
#        typelist = [('', '', '')]
#    if t:
#        return typelist
#    else:
#        return matdict
#
# #def envi_layerotype(self, context):
#    return retmatdict(self, 1, 0)
#
#def envi_layer1type(self, context):
#    return retmatdict(self, 1, 1)
#
#def envi_layer2type(self, context):
#    return retmatdict(self, 1, 2)
#
#def envi_layer3type(self, context):
#    return retmatdict(self, 1, 3)
#
#def envi_layer4type(self, context):
#    return retmatdict(self, 1, 4)
#
#def envi_layero(self, context):
#    return [((mat, mat, 'Layer material')) for mat in list(retmatdict(self, 0, 0)[self.envi_type_lo])]
#
#def envi_layer1(self, context):
#    return [((mat, mat, 'Layer material')) for mat in list(retmatdict(self, 0, 1)[self.envi_type_l1])]
#
#def envi_layer2(self, context):
#    return [((mat, mat, 'Layer material')) for mat in list(retmatdict(self, 0, 2)[self.envi_type_l2])]
#
#def envi_layer3(self, context):
#    return [((mat, mat, 'Layer material')) for mat in list(retmatdict(self, 0, 3)[self.envi_type_l3])]
#
#def envi_layer4(self, context):
#    return [((mat, mat, 'Layer material')) for mat in list(retmatdict(self, 0, 4)[self.envi_type_l4])]

