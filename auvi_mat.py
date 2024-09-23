import os, json
from collections import OrderedDict
from .vi_func import ret_datab

class auvi_materials(object):
    '''Defines auvi materials with a comma separated dictionary, with material name as key, giving
    (Roughness, Conductivity {W/m-K}, Density {kg/m3}, Specific Heat {J/kg-K}, Thermal Absorbtance,
    Solar Absorbtance, Visible Absorbtance, Default thickness)'''
    updated = 0

    def __init__(self):
        self.updated = 0

    def update(self):
        with open(ret_datab('AuVi_database.json', 'r'), 'r') as mat_jfile:
            self.mat_dict = json.loads(mat_jfile.read())

        self.updated = 1


def auvi_type_abs(self, context):
    if not self.am.updated:
        self.am.update()

    return [((mat, mat, 'Acoustic material')) for mat in self.am.mat_dict['absorption'].keys()]

def auvi_abs(self, context):
    if not self.am.updated:
        self.am.update()

    if self.auvi_type_abs:
        return [((mat, mat, 'Acoustic material')) for mat in self.am.mat_dict['absorption'][self.auvi_type_abs].keys()]
    else:
        return [('', '', '')]

def auvi_type_scatt(self, context):
    if not self.am.updated:
        self.am.update()

    return [((mat, mat, 'Acoustic material')) for mat in self.am.mat_dict['scattering'].keys()]


def auvi_scatt(self, context):
    if not self.am.updated:
        self.am.update()

    if self.auvi_type_scatt:
        return [((mat, mat, 'Acoustic material')) for mat in self.am.mat_dict['scattering'][self.auvi_type_scatt].keys()]
    else:
        return [('', '', '')]