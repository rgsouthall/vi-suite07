# EnVi materials database
import os, json
from collections import OrderedDict
from .envi_func import epentry

class envi_materials(object):
    '''Defines materials with a comma separated dictionary, with material name as key, giving 
    (Roughness, Conductivity {W/m-K}, Density {kg/m3}, Specific Heat {J/kg-K}, Thermal Absorbtance, 
    Solar Absorbtance, Visible Absorbtance, Default thickness)'''

    def __init__(self):
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
        
#        self.metal_datd = {'Copper': ('Smooth', '200', '8900.0', '418.00', '0.72', '0.65', '0.65', '5'),
#                        'Steel': ('Smooth', '50', '7800.0', '502.0', '0.12', '0.2', '0.2', '5'),
#                        'Aluminium': ('Smooth', '210', '2700', '880.00', '0.22', '0.2', '0.2', '5'),
#                        'Lead': ('Smooth', '35.3', '11340', '128.00', '0.05', '0.05', '0.05', '5')}
        self.metal_dat = OrderedDict(sorted(self.metal_datd.items()))

#        self.brick_datd = {'Standard Brick': ('Rough', '0.8', '1800', '900.00', '0.900000', '0.600000', '0.600000', '100'),
#                        'Inner brick': ('Rough', '0.62', '1800', '840.00', '0.93', '0.700000', '0.700000', '100'),
#                        'Outer brick': ('Rough', '0.96', '2000', '650.00', '0.90', '0.930000', '0.930000', '100'),
#                        'Vermiculite insulating brick': ('Rough', '0.27', '700', '840.00', '0.90', '0.650000', '0.650000', '100'),
#                        'Honeycomb brick': ('Rough', '0.27', '1700', '1000.00', '0.90', '0.7', '0.7', '102'),
#                        'Hollow terracota': ('Rough', '0.6', '845', '800', '0.90', '0.7', '0.7', '102')}
        self.brick_dat = OrderedDict(sorted(self.brick_datd.items()))

#        self.cladding_datd = {'Stucco': ('Smooth', '0.692', '1858', '836.00', '0.900000', '0.9200000', '0.920000', '25'),
#                              'Plaster board': ('Smooth', '0.7264', '1602', '836.00', '0.400000', '0.400000', '0.40000', '20'),
#                              'Plaster': ('Smooth', '1.5', '1900', '840.00', '0.300000', '0.300000', '0.3000', '5'),
#                              'Roof tiles': ('Smooth', '0.84', '1900', '840.00', '0.800000', '0.800000', '0.80000', '20')}
        self.cladding_dat = OrderedDict(sorted(self.cladding_datd.items()))

#        self.concrete_datd = {'Light mix concrete': ('MediumRough', '0.38', '1200.0', '653', '0.9', '0.65', '0.65', '100'),
#                        'Aerated concrete block': ('Rough', '0.24', '750.0', '1000', '0.9', '0.65', '0.65', '100'),
#                        'Inner concrete block': ('Rough', '0.51', '1400.0', '1000', '0.9', '0.65', '0.65', '100'),
#                        'Heavy mix concrete': ('Rough', '1.4', '2100.0', '840.0', '0.90', '0.65', '0.65', '100'),
#                        'Concrete Floor slab': ('MediumRough', '1.73', '2242.6', '836.0', '0.90', '0.65', '0.65', '100'),
#                        'Hemcrete': ('Rough', '0.09', '330.0', '2100', '0.900000', '0.600000', '0.600000', '50'),
#                        'Screed': ('MediumRough', '0.41', '1200.0', '2100', '0.900000', '0.600000', '0.600000', '50')}
        self.concrete_dat = OrderedDict(sorted(self.concrete_datd.items()))

#        self.wood_datd = {'Wood flooring': ('MediumSmooth', '0.14', '600.0', '1210.0', '0.91', '0.65', '0.65', '25'),
#                        'Parquet flooring': ('MediumSmooth', '0.17', '740.0', '2000.0', '0.90', '0.65', '0.65', '12'),
#                        'Medium hardboard': ('MediumSmooth', '0.17', '740.0', '2000.0', '0.90', '0.65', '0.65', '12'),
#                        'Plywood': ('MediumSmooth', '0.15', '700.0', '1420.0', '0.90', '0.65', '0.65', '25'),
#                        'Chipboard': ('MediumSmooth', '0.15', '800.0', '2093.0', '0.91', '0.65', '0.65', '25'),
#                        'OSB': ('MediumSmooth', '0.13', '640.0', '800.0', '0.91', '0.65', '0.65', '15'),
#                        'Hardwood': ('MediumSmooth', '0.16', '720.8', '1255.2', '0.90', '0.78', '0.78', '50')}
        self.wood_dat = OrderedDict(sorted(self.wood_datd.items()))

#        self.stone_datd = {'Sandstone': ('MediumSmooth', '1.83', '2200.0', '712.0', '0.90', '0.6', '0.6', '200'),
#                          'Limestone': ('MediumSmooth', '1.3', '2180.0', '720.0', '0.90', '0.6', '0.6', '200'),
#                          'Clay tile': ('MediumSmooth', '0.85', '1900.0', '837.0', '0.90', '0.6', '0.6', '6'),
#                          'Common earth': ('Rough', '1.28', '1460.0', '879.0', '0.90', '0.85', '0.85', '200'),
#                          'Gravel': ('Rough', '1.28', '1460.0', '879.0', '0.90', '0.85', '0.85', '200'),
#                          'Tuff': ('MediumRough', '0.4', '1400.0', '800.0', '0.90', '0.65', '0.65', '200'),
#                          'Rammed earth': ('Rough', '1.25', '1540.0', '1260.0', '0.90', '0.65', '0.65', '250'),
#                          'Sand': ('Rough', '0.2', '1500.0', '700.0', '0.90', '0.65', '0.65', '250')}
                          
        self.stone_dat = OrderedDict(sorted(self.stone_datd.items()))

#        self.gas_datd = {'Air 20-50mm': ('Gas', 'Air', '0.17'),
#                        'Horizontal Air 20-50mm Heat Down': ('Gas', 'Air', '0.21'),
#                        'Horizontal Air 20-50mm Heat Up': ('Gas', 'Air', '0.17')}
        self.gas_dat = OrderedDict(sorted(self.gas_datd.items()))

#        self.wgas_datd = {'Argon': ('Gas', 'Argon', '', '0.2', '0.016'),
#                        'Krypton':('Gas', 'Krypton', '', '0.22', '0.00943'),
#                        'Xenon':('Gas', 'Xenon', '', '0.25', '0.00565'),
#                        'Air': ('Gas', 'Air', '', '0.17', '0.024')}
        self.wgas_dat = OrderedDict(sorted(self.wgas_datd.items()))

#        self.glass_datd = {'Clear 6mm': ('Glazing', 'SpectralAverage', '', '0.006', '0.775', '0.071', '0.071', '0.881', '0.080', '0.080', '0.0', '0.84', '0.84', '0.9', 0),
#                          'Clear 4mm': ('Glazing', 'SpectralAverage', '', '0.004', '0.837', '0.075', '0.075', '0.898', '0.081', '0.081', '0.0', '0.84', '0.84', '0.9', 0),
#                          'Clear 3mm': ('Glazing', 'SpectralAverage', '', '0.003', '0.837', '0.075', '0.075', '0.898', '0.081', '0.081', '0.0', '0.84', '0.84', '0.9', 0),
#                          'Clear 6mm Soft LoE': ('Glazing', 'SpectralAverage', '', '0.006', '0.600', '0.0170', '0.220', '0.840', '0.055', '0.078', '0.0', '0.05', '0.84', '0.9', 0),
#                          'Clear 4mm Soft LoE': ('Glazing', 'SpectralAverage', '', '0.004', '0.600', '0.0170', '0.220', '0.840', '0.055', '0.078', '0.0', '0.05', '0.84', '0.9', 0),
#                          'Clear 3mm Soft LoE': ('Glazing', 'SpectralAverage', '', '0.003', '0.630', '0.190', '0.220', '0.850', '0.056', '0.079', '0.0', '0.05', '0.84', '0.9', 0),
#                          'Clear 6mm Hard LoE': ('Glazing', 'SpectralAverage', '', '0.006', '0.600', '0.0170', '0.220', '0.840', '0.055', '0.078', '0.0', '0.15', '0.84', '0.9', 0),
#                          'Clear 4mm Hard LoE': ('Glazing', 'SpectralAverage', '', '0.004', '0.600', '0.0170', '0.220', '0.840', '0.055', '0.078', '0.0', '0.15', '0.84', '0.9', 0),
#                          'Clear 3mm Hard LoE': ('Glazing', 'SpectralAverage', '', '0.003', '0.630', '0.190', '0.220', '0.850', '0.056', '0.079', '0.0', '0.15', '0.84', '0.9', 0)}
        self.glass_dat = OrderedDict(sorted(self.glass_datd.items()))

#        self.insulation_datd = {'Glass fibre quilt': ('Rough', '0.04', '12.0', '840.0', '0.9', '0.65', '0.65', '100'),
#                        'EPS': ('MediumSmooth', '0.035', '15', '1000.0', '0.90', '0.7', '0.7', '100'),
#                        'Cavity wall insul': ('Rough', '0.037', '300.0', '1000.0', '0.90', '0.6', '0.6', '100'),
#                        'Roofing felt': ('Rough', '0.19', '960.0', '837.0', '0.90', '0.9', '0.9', '6'),
#                        'Wilton wool carpet': ('Rough', '0.06', '186.0', '1360.0', '0.90', '0.60', '0.60', '5'),
#                        'Thermawall TW50': ('MediumSmooth', '0.022', '32.000', '1500', '0.900000', '0.600000', '0.600000', '200'),
#                        'Stramit': ('Rough', '0.1', '380.0', '2100', '0.900000', '0.600000', '0.600000', '50'),
#                        'Straw bale': ('Rough', '0.07', '110.0', '2000', '0.900000', '0.600000', '0.600000', '50'),
#                        'Foamglass': ('MediumSmooth', '0.04', '120.0', '840', '0.900000', '0.600000', '0.600000', '50'),
#                        'Calsitherm': ('Rough', '0.059', '220.0', '1500', '0.900000', '0.600000', '0.600000', '50'),
#                        'Cellulose (attic)': ('Rough', '0.04', '25.0', '1600', '0.900000', '0.600000', '0.600000', '200'),
#                        'Thermafloor TF70': ('Smooth', '0.022', '32.0', '1500', '0.100000', '0.100000', '0.100000', '250'),
#                        'Aerogel insulation': ('Smooth', '0.015', '2.0', '840', '0.100000', '0.100000', '0.100000', '60')}
   
        self.insulation_dat = OrderedDict(sorted(self.insulation_datd.items()))
        
#        self.pcm_datd = {'PCM plaster board': ('Smooth', '0.7264', '1602', '836.00', '0.90', '0.92', '0.92', '20'),
#                         'DuPont Energain': ('Smooth', '0.16', '850', '2500', '0.90', '0.92', '0.92', '5')}
        self.pcm_dat = OrderedDict(sorted(self.pcm_datd.items()))
        
#        self.pcmd_datd = {'PCM plaster board': ('0.0', '-20.0:0.1 22:18260 22.1:32000 60:71000'),
#                        'DuPont Energain': ('0.0', '-9.0:0.001 15.0:93760 26.0:191185 80.0:332460')}
        self.pcmd_dat = OrderedDict(sorted(self.pcmd_datd.items()))
        
        self.pv_datd = {'Default PV': ('Rough', '0.035', '29', '1213', '0.90', '0.5', '0.5', '50')}
        self.pv_dat = OrderedDict(sorted(self.pv_datd.items()))
        
        self.namedict = OrderedDict()
        self.thickdict = OrderedDict()
        self.i = 0
        self.matdat = OrderedDict()
        
        for dat in (self.brick_dat, self.cladding_dat, self.concrete_dat, self.gas_dat, self.insulation_dat, self.metal_dat, 
                    self.stone_dat, self.wood_dat, self.glass_dat, self.wgas_dat, self.pcm_dat, self.pv_dat):
            self.matdat.update(dat)
#        matdict = {'Glass': self.glass_datd, 'Metal': self.metal_datd, 'Brick': self.brick_datd, 'Cladding': self.cladding_datd,
#                   'Concrete': self.concrete_datd, 'Wood': self.wood_datd, 'Stone': self.stone_datd, 'Gas': self.gas_datd,
#                   'WGas': self.wgas_datd, 'Insulation': self.insulation_datd, 'PCM': self.pcm_datd, 'PCMD': self.pcmd_datd}    
#        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'EPFiles', '{}'.format('Material_database.json')), 'w') as mat_jfile:
#            mat_jfile.write(json.dumps(matdict))
#        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'EPFiles', '{}'.format('Material_database.json')), 'r') as mat_jfile:
#            mat_jfile.write(json.dumps(matdict))
            
    def get_dat(self, mat_type):
        mat_dict = {'Glass': self.glass_datd, '3': self.metal_datd, '0': self.brick_datd, '1': self.cladding_datd,
                   '2': self.concrete_datd, '5': self.wood_datd, '4': self.stone_datd, '6': self.gas_datd,
                   'WGas': self.wgas_datd, '7': self.insulation_datd, '8': self.pcm_datd, '9': self.pcmd_datd} 
        return mat_dict[mat_type]
    
    def lay_save(self):
        mat_dict = {'Glass': self.glass_datd, 'Metal': self.metal_datd, 'Brick': self.brick_datd, 'Cladding': self.cladding_datd,
                   'Concrete': self.concrete_datd, 'Wood': self.wood_datd, 'Stone': self.stone_datd, 'Gas': self.gas_datd,
                   'WGas': self.wgas_datd, 'Insulation': self.insulation_datd, 'PCM': self.pcm_datd, 'PCMD': self.pcmd_datd}    
        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'EPFiles', '{}'.format('Material_database.json')), 'w') as mat_jfile:
            mat_jfile.write(json.dumps(mat_dict))
        
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
#        self.wall_cond = {'External Wall 1': ('Standard Brick', 'Thermawall TW50', 'Inner concrete block'), 'Kingston PH 1': ('Plywood', 'EPS', 'Plywood'),
#        'SIP': ('OSB', 'EPS', 'OSB')}
        self.wall_con = OrderedDict(sorted(self.wall_cond.items()))
        
#        self.iwall_cond = {'Party Wall 1': ('Plaster board', 'Standard Brick', 'Plaster board')}
        self.iwall_con = OrderedDict(sorted(self.iwall_cond.items()))
        
#        self.ceil_cond = {'Ceiling 1': ('Chipboard', 'EPS', 'Plaster board')}
        self.ceil_con = OrderedDict(sorted(self.ceil_cond.items()))
        
#        self.ifloor_cond = {'Internal floor 1': ('Plaster board', 'EPS', 'Chipboard')}
        self.ifloor_con = OrderedDict(sorted(self.ifloor_cond.items()))
        
#        self.floor_cond = {'Ground Floor 1': ('Common earth', 'Gravel', 'Heavy mix concrete', 'Horizontal Air 20-50mm Heat Down', 'Chipboard'),
#                           'Kingston PH 1': ('Common earth', 'Gravel', 'EPS', 'Heavy mix concrete')}
        self.floor_con = OrderedDict(sorted(self.floor_cond.items()))
        
#        self.roof_cond = {'Roof 1': ('Clay tile', 'Roofing felt', 'Plywood')}
        self.roof_con = OrderedDict(sorted(self.roof_cond.items()))

#        self.door_cond = {'Internal Door 1': ('Chipboard', 'Hardwood', 'Chipboard')}
        self.door_con = OrderedDict(sorted(self.door_cond.items()))
        
        self.pv_cond = {'Simple PV': ['Default PV']}
        self.pv_con = OrderedDict(sorted(self.pv_cond.items()))
        
#        self.glaze_cond = {'Standard Double Glazing': ('Clear 3mm', 'Air', 'Clear 3mm'), 'Low-E Double Glazing': ('Clear 3mm', 'Air', 'Clear 3mm Hard LoE'), 
#                           'PassivHaus': ('Clear 3mm', 'Argon', 'Clear 3mm Soft LoE', 'Argon', 'Clear 3mm Soft LoE'), 'Velfac 200 Double': ('Clear 4mm', 'Argon', 'Clear 4mm Soft LoE'),
#                           'Velfac 200 Triple': ('Clear 4mm', 'Argon', 'Clear 4mm Soft LoE', 'Argon', 'Clear 4mm Soft LoE')}
        self.glaze_con = OrderedDict(sorted(self.glaze_cond.items()))
        
        self.p = 0
        self.propdict = {'Wall': self.wall_con, 'Floor': self.floor_con, 'Roof': self.roof_con, 'Ceiling': self.ceil_con, 'Door': self.door_con, 
                                'Window': self.glaze_con, 'PV': self.pv_con, 'Internal floor': self.ifloor_con, 'Internal wall': self.iwall_con} 
#        condict = {'Wall': self.wall_cond, 'Party wall': self.iwall_cond, 'Ceiling': self.ceil_cond, 'Internal floor': self.ifloor_cond,
#                   'Floor': self.floor_cond, 'Roof': self.roof_cond, 'Door': self.door_cond, 'Glazing': self.glaze_cond}
#
#        with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'EPFiles', '{}'.format('Construction_database.json')), 'w') as con_jfile:
#            con_jfile.write(json.dumps(condict))
            
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
        
def retmatdict(ect, t, l): 
    if ect in ('Wall', 'Roof', 'Floor', 'Door', 'Ceiling', 'Frame'):
        typelist = [("0", "Brick", "Choose a material from the brick database"),("1", "Cladding", "Choose a material from the cladding database"), 
                    ("2", "Concrete", "Choose a material from the concrete database"),("3", "Metal", "Choose a material from the metal database"),
                   ("4", "Stone", "Choose a material from the stone database"),("5", "Wood", "Choose a material from the wood database"),
                   ("6", "Gas", "Choose a material from the gas database"),("7", "Insulation", "Choose a material from the insulation database"),
                    ("8", "PCM", "Choose a material from the phase change database"), ("9", "PV", "Choose a material from the photovoltaic database")]
        matdict = {'0': envi_materials().brick_dat.keys(), '1': envi_materials().cladding_dat.keys(), '2': envi_materials().concrete_dat.keys(), '3': envi_materials().metal_dat.keys(), '4': envi_materials().stone_dat.keys(),
                   '5': envi_materials().wood_dat.keys(), '6': envi_materials().gas_dat.keys(), '7': envi_materials().insulation_dat.keys(), '8': envi_materials().pcm_dat.keys(), '9': envi_materials().pv_dat.keys()}

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
            
def envi_layerotype(self, context):   
    return retmatdict(self, 1, 0) 
    
def envi_layer1type(self, context):   
    return retmatdict(self, 1, 1) 
    
def envi_layer2type(self, context):   
    return retmatdict(self, 1, 2) 
    
def envi_layer3type(self, context):   
    return retmatdict(self, 1, 3) 
    
def envi_layer4type(self, context):   
    return retmatdict(self, 1, 4) 
                    
def envi_layero(self, context):   
    return [((mat, mat, 'Layer material')) for mat in list(retmatdict(self, 0, 0)[self.envi_type_lo])]

def envi_layer1(self, context):
    return [((mat, mat, 'Layer material')) for mat in list(retmatdict(self, 0, 1)[self.envi_type_l1])]

def envi_layer2(self, context):
    return [((mat, mat, 'Layer material')) for mat in list(retmatdict(self, 0, 2)[self.envi_type_l2])]

def envi_layer3(self, context):
    return [((mat, mat, 'Layer material')) for mat in list(retmatdict(self, 0, 3)[self.envi_type_l3])]

def envi_layer4(self, context):
    return [((mat, mat, 'Layer material')) for mat in list(retmatdict(self, 0, 4)[self.envi_type_l4])]
    
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
    return retmatdict(self.envi_con_type, 1, self.bl_idname == 'envi_gl_node')   
          
def envi_layer(self, context):  
    if self.materialtype:
        return [((mat, mat, 'Layer material')) for mat in list(retmatdict(self.envi_con_type, 0, self.bl_idname == 'envi_gl_node')[self.materialtype])]  
    else:
        return [('', '', '')]    

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
     