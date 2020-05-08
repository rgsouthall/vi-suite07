unit2res = {'Lux': 'illu', 'DF (%)': 'df', u'W/m\u00b2 (v)': 'virrad', u'W/m\u00b2 (f)': 'firrad', 'SVF': 'svf', 'sDA (%)': 'sda', 'ASE (hrs)': 'ase',
             'lxh': 'illuh','kWh (f)': 'firradh', 'kWh/m2 (f)': 'firradhm2', 'kWh (v)': 'virradh', u'kWh/m\u00b2 (v)': 'virradhm2', 'DA (%)': 'da', 
             'UDI-f (%)': 'udilow', 'UDI-s (%)': 'udisup', '% Sunlit': 'sm',
             'UDI-a (%)': 'udiauto', 'UDI-e (%)': 'udihi', 'kWh': 'kwh', 'kWh/m2': 'kwhm2', 'Lux (max)': 'maxlux', 'Lux (min)': 'minlux', 'Lux (ave)': 'avelux'}

res2unit = {unit2res[u]: u for u in unit2res}

colours = [('rainbow', 'Rainbow', 'Rainbow colour scale'), ('gray', 'Grey', 'Grey colour scale'), ('hot', 'Hot', 'Hot colour scale'),
           ('CMRmap', 'CMR', 'CMR colour scale'), ('jet', 'Jet', 'Jet colour scale'), ('plasma', 'Plasma', 'Plasma colour scale'), 
           ('hsv', 'HSV', 'HSV colour scale'), ('viridis', 'Viridis', 'Viridis colour scale')]

coldict = {'0': 'rainbow', '1': 'gray', '2': 'hot', '3': 'CMRmap', '4': 'jet', '5': 'plasma'}

e1ddict = {'BPsolar 275': (4.75,  21.4, 17, 4.45, 0.00065, -0.08, 36, 320, 0.63),
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