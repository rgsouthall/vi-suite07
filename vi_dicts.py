rpictparams = {"-ab": (2, 4, 8), "-ad": (256, 1024, 4096), "-as": (128, 512, 2048), "-aa": (0, 0, 0), "-ar": (8, 32, 128),
               "-ss": (0, 2, 5), "-st": (1, 0.75, 0.1), "-lw": (0.0005, 0.0001, 0.00002), "-lr": (2, 4, 8),
               "-pj": (0, 0, 0), "-pt": (0.15, 0.05, 0.0), "-ps": (1, 1, 1),
               "-dj": (0.0, 0.7, 1), "-ds": (0, 0.5, 0.15), "-dr": (0, 1, 3), "-dt": (0.5, 0.05, 0.0), "-dc": (0.25, 0.5, 1)}

rvuparams = {"-ab": (2, 3, 4), "-ad": (256, 1024, 4096), "-as": (128, 512, 2048), "-aa": (0, 0, 0), "-ar": (8, 32, 128),
             "-dj": (0, 0.7, 1), "-ds": (0.5, 0.15, 0.15), "-dr": (1, 3, 5),
             "-ss": (0, 2, 5), "-st": (1, 0.75, 0.1), "-lw": (0.0005, 0.0001, 0.00002),
             "-lr": (2, 3, 4), "-ps": (4, 1, 1)}

rtraceparams = {"-ab": (2, 4, 8), "-ad": (256, 1024, 4096), "-as": (128, 512, 2048), "-aa": (0, 0, 0), "-ar": (8, 32, 128),
                "-dj": (0, 0.7, 1), "-ds": (0, 0.5, 0.15), "-dr": (1, 3, 5), "-ss": (0, 2, 5), "-st": (1, 0.75, 0.1),
                "-lw": (0.0001, 0.00001, 0.000002), "-lr": (2, 3, 4)}

rtracecbdmparams = {"-ab": (4, 8), "-ad": (4096, 8192), "-as": (512, 1024), "-aa": (0, 0), "-ar": (8, 32, 128),
                    "-dj": (0.7, 1), "-ds": (0.5, 0.15), "-dr": (2, 3), "-ss": (2, 5), "-st": (0.75, 0.1),
                    "-ss": (2, 5), "-st": (0.75, 0.1),
                    "-lw": (1e-4, 1e-7), "-lr": (3, 5)}

unit2res = {'Lux': 'illu', 'DF (%)': 'df', 'W/m2 (v)': 'virradm2', 'W/m2 (f)': 'firradm2', 'W (f)': 'firrad',
            'W/m2': 'firradm2', 'W': 'firrad', 'W/m2 (red)': 'firradrm2', 'W (red)': 'firradr',
            'W/m2 (green)': 'firradgm2', 'W (green)': 'firradg', 'W/m2 (blue)': 'firradbm2', 'W (blue)': 'firradb',
            'SVF (%)': 'svf', 'sDA (%)': 'sda', 'ASE (hrs)': 'ase', 'Perimeter area': 'sv',
            'klxh': 'illuh', 'kWh (f)': 'firradh', 'kWh/m2 (f)': 'firradhm2', 'kWh (v)': 'virradh',
            'kWh/m2 (v)': 'virradhm2', 'DA (%)': 'da', 'UDI-f (%)': 'udilow', 'UDI-s (%)': 'udisup',
            'Sunlit time (%)': 'sm', 'UDI-a (%)': 'udiauto', 'UDI-e (%)': 'udihi', 'kWh': 'kwh',
            'kWh/m2': 'kwhm2', 'Lux (max)': 'maxlux', 'Lux (min)': 'minlux', 'Lux (ave)': 'avelux'}

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

rvuerrdict = {'view up parallel to view direction': "Camera cannot point directly upwards",
              ' x11': "No X11 display server found. You may need to install XQuartz",
              'source center': "A light source has concave faces. Use mesh - cleanup - split concave faces"}

pmerrdict = {'fatal - too many prepasses, no global photons stored\n': "Too many prepasses have occurred. Make sure light sources can see your geometry",
             'fatal - too many prepasses, no global photons stored, no caustic photons stored\n': "Too many prepasses have occurred. Turn off caustic photons and encompass the scene",
             'fatal - zero flux from light sources\n': "No light flux, make sure there is a light source and that photon port normals point inwards",
             'fatal - no light sources in distribPhotons\n': "No light sources. Photon mapping does not work with HDR skies",
             'fatal - no valid photon ports found\n': 'Make sure photon ports are valid',
             'fatal - failed photon distribution\n': 'Do the lights see enough geometry?'}

rnu = {'0': ("Air", "Ambient air metrics"), '1': ("Wind Speed", "Ambient Wind Speed (m/s)"), '2': ("Wind Direction", "Ambient Wind Direction (degrees from North)"),
       '3': ("Humidity", "Ambient Humidity"), '4': ("Solar", 'Ambient solar metrics'), '5': ("Temperature", "Zone Temperature"), '6': ("Humidity", "Zone Humidity"),
       '7': ("Heating Watts", "Zone Heating Requirement (Watts)"), '8': ("Cooling Watts", "Zone Cooling Requirement (Watts)"),
       '9': ("Solar Gain", "Window Solar Gain (Watts)"), '10': ("PPD", "Percentage Proportion Dissatisfied"), '11': ("PMV", "Predicted Mean Vote"),
       '12': ("Ventilation (l/s)", "Zone Ventilation rate (l/s)"), '13': (u'Ventilation (m\u00b3/h)', u'Zone Ventilation rate (m\u00b3/h)'),
       '14': (u'Infiltration (m\u00b3/hr)',  u'Zone Infiltration (m\u00b3/hr)'), '15': ('Infiltration (ACH)', 'Zone Infiltration rate (ACH)'), '16': ('CO2 (ppm)', 'Zone CO2 concentration (ppm)'),
       '17': ("Heat loss (W)", "Ventilation Heat Loss (W)"), '18': (u'Flow (m\u00b3/s)', u'Linkage flow (m\u00b3/s)'), '19': ('Opening factor', 'Linkage Opening Factor'),
       '20': ("MRT (K)", "Mean Radiant Temperature (K)"), '21': ('Occupancy', 'Occupancy count'), '22': ("Humidity", "Zone Humidity"),
       '23': ("Fabric HB (W)", "Fabric convective heat balance"), '24': ("Air Heating", "Zone air heating"), '25': ("Air Cooling", "Zone air cooling"),
       '26': ("HR Heating", "Heat recovery heating (W)"), '27': ("Volume flow", "Thermal chimney volume flow rate (m3/2)"), '28': ("Mass flow", "Thermal chmimney mass flow rate (kg/s"),
       '29': ("Out temp.", "Thermal chimney outlet temperature (C)"), '30': ("Heat loss", "Thermal chimney heat loss (W)"), '31': ("Heat gain", "Thermal chimney heat gain (W)"),
       '32': ("Volume", "Thermal chimnwey volume (m3)"), '33': ("Mass", "Thermal chimney mass (kg)"), '34': ('delta P', 'Linkage Pressure Differential (Pa)'),
       '35': ('Equipment', 'Other equipment heat gains (W)'), '36': ('PV Energy', 'PV energy (J)'), '37': ('PV Power', 'PV power (W)'),
       '38': ('PV Temp.', 'PV Temperature (C)'), '39': ('PV Eff.', 'PV efficiency (%)')}

arnu = {'0': (u"Max temp (\u2103)", "Maximum zone temperature"), '1': (u"Min temp (\u2103)", "Minimum zone temperature"), '2': (u"Avg temp (\u2103)", "Average zone temperature"),
        '3': ("Max heating (W)", "Max Zone heating"), '4': ("Min heating (W)", "Min Zone heating"), '5': ("Avg heating (W)", "Avg Zone heating"),
        '6': ("Total heating (kWh)", "Total zone heating"), '7': (u"Total heating (kWh/m\u00b2)", "Total zone heating per floor area"),
        '8': ("Max cooling (W)", "Max Zone cooling"), '9': ("Min cooling (W)", "Min Zone cooling"), '10': ("Avg cooling (W)", "Avg Zone colling"),
        '11': ("Total cooling (kWh)", "Total zone cooling"), '12': (u"Total cooling (kWh/m\u00b2)", "Total zone cooling per floor area"),
        '13': ("Max CO2 (ppm)", u"Maximum zone CO\u2082 level"), '14': ("Avg CO2 (ppm)", u"Average zone CO\u2082 level"), '15': ("Min CO2 (ppm)", u"Minimum zone CO\u2082 level"),
        '16': (u"Max flow in (m\u00b3/s)", u"Maximum linkage flow level"), '17': (u"Min flow in (m\u00b3/s)", u"Minimum linkage flow level"),
        '18': (u"Avg flow in (m\u00b3/s)", u"Average linkage flow level"),
        '19': ('Max SHG (W)', 'Maximum Solar Heat Gain'), '20': ('Min SHG (W)', 'Minimum Solar Heat Gain'), '21': ('Avg SHG (W)', 'Average Solar Heat Gain'),
        '22': ('Total SHG (kWh)', 'Total solar heat gain'), '23': ('Total SHG (kWh/m2)', 'Total Solar heat gain per floor area')}

envdict = {'Site Outdoor Air Drybulb Temperature [C] !Hourly': "Temperature (degC)",
           'Site Outdoor Air Relative Humidity [%] !Hourly': 'Humidity (%)',
           'Site Wind Direction [deg] !Hourly': 'Wind Direction (deg)',
           'Site Wind Speed [m/s] !Hourly': 'Wind Speed (m/s)',
           'Site Diffuse Solar Radiation Rate per Area [W/m2] !Hourly': "Diffuse Solar (W/m^2)",
           'Site Direct Solar Radiation Rate per Area [W/m2] !Hourly': "Direct Solar (W/m^2)"}

zresdict = {'Zone Air Temperature [C] !Hourly': "Temperature (degC)",
            'Zone Air Relative Humidity [%] !Hourly': 'Humidity (%)',
            'Zone Air System Sensible Heating Rate [W] !Hourly': 'Heating (W)',
            'Zone Air System Sensible Cooling Rate [W] !Hourly': 'Cooling (W)',
            'Zone Ideal Loads Supply Air Sensible Heating Rate [W] !Hourly': 'Air heating (W)',
            'Zone Ideal Loads Heat Recovery Sensible Heating Rate [W] !Hourly': 'HR heating (W)',
            'Zone Ideal Loads Supply Air Sensible Cooling Rate [W] !Hourly': 'Air cooling (W)',
            'Zone Windows Total Transmitted Solar Radiation Rate [W] !Hourly': 'Solar gain (W)',
            'AFN Zone Infiltration Volume [m3] !Hourly': 'Infiltration (m3/hr)',
            'AFN Zone Infiltration Air Change Rate [ach] !Hourly': 'Infiltration (ACH)',
            'Zone Infiltration Current Density Volume [m3] !Hourly': 'Infiltration (m3/hr)',
            'Zone Infiltration Air Change Rate [ach] !Hourly': 'Infiltration (ACH)',
            'Zone Exfiltration Sensible Heat Transfer Rate [W] !Hourly': 'Exfiltration heat loss (W)',
            'Zone Mean Air Temperature [C] ! Hourly': 'Mean Temperature (degC)',
            'Zone Thermal Comfort Fanger Model PPD [%] !Hourly': 'PPD (%)',
            'Zone Thermal Comfort Fanger Model PMV [] !Hourly': 'PMV',
            'AFN Node CO2 Concentration [ppm] !Hourly': 'CO2 (ppm)',
            'Zone Air CO2 Concentration [ppm] !Hourly': 'CO2 (ppm)',
            'Zone Mean Radiant Temperature [C] !Hourly': 'MRT',
            'Zone People Occupant Count [] !Hourly': 'Occupancy',
            'Zone Air Heat Balance Surface Convection Rate [W] !Hourly': 'Heat balance (W)',
            'Zone Air Heat Balance Internal Convective Heat Gain Rate [W] !Hourly': 'Internal heat balance (W)',
            'Zone Thermal Chimney Current Density Air Volume Flow Rate [m3/s] !Hourly': 'Volume flow (m3/s)',
            'Zone Thermal Chimney Mass Flow Rate [kg/s] !Hourly': 'Mass flow (kg/s)',
            'Zone Thermal Chimney Outlet Temperature [C] !Hourly': 'Outlet temperature (C)',
            'Zone Thermal Chimney Heat Loss Energy [J] !Hourly': 'TC heat loss (J)',
            'Zone Thermal Chimney Heat Gain Energy [J] !Hourly': 'TC heat gain (J)',
            'Zone Thermal Chimney Volume [m3] !Hourly': 'TC VOLUME (m3)',
            'Zone Thermal Chimney Mass [kg] !Hourly': 'TC mass(kg)',
            'Zone Other Equipment Total Heating Rate [W] !Hourly': 'Equipment (W)'}

enresdict = {'AFN Node CO2 Concentration [ppm] !Hourly': 'CO2 (ppm)',
             'AFN Node Wind Pressure [Pa] !Hourly': 'WP (PA)'}

lresdict = {'AFN Linkage Node 1 to Node 2 Volume Flow Rate [m3/s] !Hourly': 'Linkage Flow out',
            'AFN Linkage Node 2 to Node 1 Volume Flow Rate [m3/s] !Hourly': 'Linkage Flow in',
            'AFN Surface Venting Window or Door Opening Factor [] !Hourly': 'Opening Factor',
            'AFN Linkage Node 1 to Node 2 Pressure Difference [Pa] !Hourly': 'delta P (Pa)'}

presdict = {'Generator Produced DC Electricity Energy [J] !Hourly': 'PV energy (J)',
            'Generator Produced DC Electricity Rate [W] !Hourly': 'PV power (W)',
            'Generator PV Array Efficiency [] !Hourly': 'PV efficiency',
            'Generator PV Cell Temperature [C] !Hourly': 'PV temperature (C)'}

hdict = {}
