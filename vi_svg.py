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

import os
from math import sin, cos, pi


def vi_info(node, dim, svp, **kwargs):
    if node.metric == '0' and node.energy_menu == '0':
        pass
#        x = 50
#        y = 100
#        svg_str = '''<path style="fill:#f8fb00;fill-opacity:1;stroke:#000000;stroke-width:0.999989;stroke-linecap:butt;stroke-linejoin:miter;stroke-dasharray:none;stroke-opacity:1"
#                    d="M 90.457354,0.35500172 1.214155,88.866548 H 74.819729 L 36.551591,164.64998 126.80966,72.577488 H 56.922653 Z"
#                    id="path3397" />'''
#     elif node.metric == '1' and node.light_menu == '3':
#         ir = kwargs['ir']
#         aDF = kwargs['aDF']
#         irscaled = ir if ir < 0.7 else 0.7
#         adfscaled = aDF if aDF < 3.2 else 3.2
#         adfxpos = 400 - (300 * sin(adfscaled*pi/2.0))
#         adfypos = 400 + (300 * cos(adfscaled*pi/2.0))
#         (adffill, adfsweep) = ("255, 0, 0", 0) if aDF < 2.0 else ("0, 255, 0", 1)
#         irxpos = 400 - (250 * sin(irscaled*pi/0.4))
#         irypos = 400 + (250 * cos(irscaled*pi/0.4))
#         (irfill, irsweep) = ("255, 0, 0", 0) if ir < 0.4 else ("0, 255, 0", 1)
#         imname = "RIBA_lighting_{}".format(node.zone_menu)
#         svg_str = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
#         <svg
#         id="svg5"
#         version="1.1"
#         viewBox="0 0 {0} {0}"
#         height="{0}"
#         width="{0}"
#         xmlns="http://www.w3.org/2000/svg"
#         xmlns:svg="http://www.w3.org/2000/svg">

#         <rect style="fill:rgb(255, 255, 255)" width="{0}" height="{0}"/>
#         <path style="fill:rgb({1})" d="M 400 700
#             A 300 300, 1, {2}, 1, {3:.0f} {4:.0f}
#             L 400 400 Z"/>
#         <circle  style="fill:rgb(255, 255, 255)" cx="400" cy="400" r="260"/>
#         <path style="fill:rgb({5})" d="M 400 650
#             A 250 250, 1, {6}, 1, {7:.0f} {8:.0f}
#             L 400 400 Z"/>
#         <circle  style="fill:rgb(255, 255, 255)" cx="400" cy="400" r="210"/>
#         <text text-anchor="middle" x="400" y="400" style="font-size: 48px">RIBA 2030</text>
#         <text text-anchor="middle" x="400" y="450" style="font-size: 48px">Lighting</text>
#         <text x="20" y="770" style="font-size: 48px">Sensor: {9}</text>
#         </svg>
#         """.format(dim, adffill, adfsweep, adfxpos, adfypos, irfill, irsweep, irxpos, irypos, node.zone_menu)

#         return imname, bytearray(svg_str, encoding='utf-8')

#     elif node.metric == '1' and node.light_menu == 'old':
#         ir = kwargs['ir']
#         aDF = kwargs['aDF']
#         adfpos = 700 - aDF * 150 if aDF < 4 else 100
#         adfheight = 700 - adfpos
#         irpos = 700 - ir * 750 if ir < 600/750 else 100
#         irheight = 700 - irpos
#         adffill = "255, 0, 0" if aDF < 2.0 else "0, 255, 0"
#         irfill = "255, 0, 0" if ir < 0.4 else "0, 255, 0"
#         imname = "RIBA_lighting_{}".format(node.zone_menu)
#         svg_str = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
#         <svg
#         id="svg5"
#         version="1.1"
#         viewBox="0 0 {0[0]} {0[1]}"
#         width="{0[0]}"
#         height="{0[1]}"
#         xmlns="http://www.w3.org/2000/svg"
#         xmlns:svg="http://www.w3.org/2000/svg">
#         <defs
#         id="defs2" />
#         <rect style="fill:rgb(255, 255, 255)" width="{0[0]}" height="{0[1]}"/>
#         <rect style="fill:rgb({1})" ry="5" x="50" y="{2}" width="100" height="{3}"/>
#         <rect style="fill:rgb({4})" ry="5" x="175" y="{5}" width="100" height="{6}"/>

#         <text text-anchor="middle" x="400" y="50" style="font-size: 48px">RIBA 2030 Lighting</text>
#         <text text-anchor="middle" x="100" y="{7}" style="font-size: 36px">{8}</text>
#         <text text-anchor="middle" x="225" y="{9}" style="font-size: 36px">{10}</text>
#         <text text-anchor="middle" x="100" y="750" style="font-size: 48px">aDF</text>
#         <text text-anchor="middle" x="225" y="750" style="font-size: 48px">IR</text>
#         <text x="275" y="750" style="font-size: 48px">Sensor: {11}</text>
#         </svg>
#         """.format(dim, adffill, adfpos, adfheight, irfill, irpos, irheight, adfpos - 25, aDF, irpos - 25, ir, node.zone_menu)


#         return imname, bytearray(svg_str, encoding='utf-8')

    elif node.metric == '1' and node.light_menu == '0':
        if kwargs.get('aDF'):
            ir = kwargs['ir']
            aDF = kwargs['aDF']
            adfpos = 335 - aDF * 67.5 if aDF < 2 else 100 + 200/aDF
            adfheight = 335 - adfpos
            irpos = 335 - ir * 337.5 if ir < 0.3 else 100 + 25/ir
            irheight = 335 - irpos
            adffill = "255, 0, 0" if aDF < 2.0 else "0, 255, 0"
            irfill = "255, 0, 0" if ir < 0.3 else "0, 255, 0"
            imname = "BREEAM_lighting_{}".format(node.zone_menu)
            svg_str = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
            <svg
            id="svg5"
            version="1.1"
            viewBox="0 0 400 400"
            width="{0[0]}"
            height="{0[1]}"
            xmlns="http://www.w3.org/2000/svg"
            xmlns:svg="http://www.w3.org/2000/svg">
            <defs>
                <linearGradient id="lGadf" x1="60" x2="190" y1="335" y2="{7:.3f}" gradientUnits="userSpaceOnUse">
                    <stop style="stop-color:rgb({1})" offset=".26316"/>
                    <stop style="stop-color:rgb({1});stop-opacity:.6392" offset="1"/>
                </linearGradient>
                <linearGradient id="lGir" x1="210" x2="340" y1="335" y2="{9:.3f}" gradientUnits="userSpaceOnUse">
                    <stop style="stop-color:rgb({4})" offset=".26316"/>
                    <stop style="stop-color:rgb({4});stop-opacity:.6392" offset="1"/>
                </linearGradient>
            </defs>

            <rect style="fill:rgb(255, 255, 255)" width="{0[0]}" height="{0[1]}"/>
            <text text-anchor="middle" x="200" y="32" style="font-size:18px;font-family:arial">BREEAM HEA 01 - Daylighting</text>
            <text text-anchor="middle" x="125" y="60" style="font-size:26px;font-family:arial">DF</text>
            <text text-anchor="middle" x="275" y="60" style="font-size:26px;font-family:arial">UR</text>
            <text text-anchor="middle" x="125" y="75" style="font-size:13px;font-family:arial">Average Daylight Factor</text>
            <text text-anchor="middle" x="275" y="75" style="font-size:13px;font-family:arial">Uniformity Ratio</text>
            <rect ry="4" x="60" y="{2:.3f}" width="130" height="{3:.3f}" style="fill:rgb({1});fill-rule:evenodd;fill:url(#lGadf);stroke-width:0.5;stroke:#000000"/>
            <rect ry="4" x="210" y="{5:.3f}" width="130" height="{6:.3f}" style="fill:rgb({4});fill-rule:evenodd;fill:url(#lGir);stroke-width:0.5;stroke:#000000"/>
            <text text-anchor="middle" x="30" y="209" style="font-size:24px;font-family:arial">2</text>
            <text text-anchor="middle" x="370" y="209" style="font-size:24px;font-family:arial">0.3</text>
            <path d="m50 200h300" style="fill:none;stroke-dasharray:5, 5;stroke-width:1;stroke:#4d4d4d"/>
            <path d="m50 340h300" style="fill:none;stroke-width:1;stroke:#000000"/>
            <path d="m200 80v270" style="fill:none;stroke-width:1;stroke:#000000"/>
            <text text-anchor="middle" x="125" y="365" style="font-size:24px;font-family:arial">{8}</text>
            <text text-anchor="middle" x="275" y="365" style="font-size:24px;font-family:arial">{10}</text>
            <text x="10" y="380" style="font-size:14px;font-family:arial">Sensor location: {11}</text>
            <text x="10" y="395" style="font-size:14px;font-family:arial">Compliant area: {12}%</text>
            </svg>
            """.format(dim, adffill, adfpos, adfheight, irfill, irpos, irheight, 335 - adfheight, aDF, 335 - irheight, ir, node.zone_menu, int(node['res']['b_area'][0] * 100))

            with open(os.path.join(svp['viparams']['newdir'], 'images', 'BREEAM_{}_light.svg'.format(node.zone_menu)), 'w') as svg_file:
                svg_file.write(svg_str)

        return imname, bytearray(svg_str, encoding='utf-8')

    elif node.metric == '1' and node.light_menu == '2':
        ir = kwargs['ir']
        aDF = kwargs['aDF']
        adfpos = 335 - aDF * 67.5 if aDF < 2 else 100 + 200/aDF
        adfheight = 335 - adfpos
        irpos = 335 - ir * 337.5 if ir < 0.4 else 100 + 25/ir
        irheight = 335 - irpos
        adffill = "255, 0, 0" if aDF < 2.0 else "0, 255, 0"
        irfill = "255, 0, 0" if ir < 0.4 else "0, 255, 0"
        imname = "RIBA_lighting_{}".format(node.zone_menu)
        svg_str = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
        <svg
        id="svg5"
        version="1.1"
        viewBox="0 0 400 400"
        width="{0[0]}"
        height="{0[1]}"
        xmlns="http://www.w3.org/2000/svg"
        xmlns:svg="http://www.w3.org/2000/svg">
        <defs>
            <linearGradient id="lGadf" x1="60" x2="190" y1="335" y2="{7:.3f}" gradientUnits="userSpaceOnUse">
                <stop style="stop-color:rgb({1})" offset=".26316"/>
                <stop style="stop-color:rgb({1});stop-opacity:.6392" offset="1"/>
            </linearGradient>
            <linearGradient id="lGir" x1="210" x2="340" y1="335" y2="{9:.3f}" gradientUnits="userSpaceOnUse">
                <stop style="stop-color:rgb({4})" offset=".26316"/>
                <stop style="stop-color:rgb({4});stop-opacity:.6392" offset="1"/>
            </linearGradient>
        </defs>

        <rect style="fill:rgb(255, 255, 255)" width="{0[0]}" height="{0[1]}"/>
        <text text-anchor="middle" x="200" y="32" style="font-size:18px;font-family:arial">RIBA 2030 Climate Challenge - Daylighting</text>
        <text text-anchor="middle" x="125" y="60" style="font-size:26px;font-family:arial">DF</text>
        <text text-anchor="middle" x="275" y="60" style="font-size:26px;font-family:arial">UR</text>
        <text text-anchor="middle" x="125" y="75" style="font-size:13px;font-family:arial">Average Daylight Factor</text>
        <text text-anchor="middle" x="275" y="75" style="font-size:13px;font-family:arial">Uniformity Ratio</text>
        <rect ry="4" x="60" y="{2:.3f}" width="130" height="{3:.3f}" style="fill:rgb({1});fill-rule:evenodd;fill:url(#lGadf);stroke-width:0.5;stroke:#000000"/>
        <rect ry="4" x="210" y="{5:.3f}" width="130" height="{6:.3f}" style="fill:rgb({4});fill-rule:evenodd;fill:url(#lGir);stroke-width:0.5;stroke:#000000"/>
        <text text-anchor="middle" x="30" y="209" style="font-size:24px;font-family:arial">2.0</text>
        <text text-anchor="middle" x="370" y="209" style="font-size:24px;font-family:arial">0.4</text>
        <path d="m50 200h300" style="fill:none;stroke-dasharray:5, 5;stroke-width:1;stroke:#4d4d4d"/>
        <path d="m50 340h300" style="fill:none;stroke-width:1;stroke:#000000"/>
        <path d="m200 80v270" style="fill:none;stroke-width:1;stroke:#000000"/>
        <text text-anchor="middle" x="125" y="365" style="font-size:24px;font-family:arial">{8}</text>
        <text text-anchor="middle" x="275" y="365" style="font-size:24px;font-family:arial">{10}</text>
        <text x="10" y="390" style="font-size:16px;font-family:arial">Sensor location: {11}</text>
        </svg>
        """.format(dim, adffill, adfpos, adfheight, irfill, irpos, irheight, 335 - adfheight, aDF, 335 - irheight, ir, node.zone_menu)

        with open(os.path.join(svp['viparams']['newdir'], 'images', 'RIBA_{}_light.svg'.format(node.zone_menu)), 'w') as svg_file:
            svg_file.write(svg_str)

        return imname, bytearray(svg_str, encoding='utf-8')

    elif node.metric == '1' and node.light_menu == '1':
        comps = ['&lt;', '>=', '>=', '>=']
        sda = kwargs['sda']
        sdapass = kwargs['sdapass']
        ase = kwargs['ase']
        asepass = kwargs['asepass']
        comps[0] = '&lt;' if ase < asepass else '>='
        tcredits = kwargs['tc']
        credits = kwargs['o1']
        asepoints = tcredits if ase < asepass else '0'
        totarea = kwargs['totarea']
        svarea = kwargs['svarea']
        comps[3] = '>=' if svarea/totarea >= 0.9 else '&lt;'
        imname = "LEED_lighting_{}".format(node.zone_menu)

        if node.leed_menu:
            sda1points = 1 if sda >= sdapass[0] else 0
        else:
            sda1points = 2 if sda >= sdapass[0] else 0

        comps[1] = '>=' if sda >= sdapass[0] else '&lt;'
        sda2points = 1 if sda >= sdapass[1] else 0
        comps[2] = '>=' if sda >= sdapass[1] else '&lt;'
        '''<text x="350" y="115" style="font-size: 20px;font-family:arial">Perimeter: {:.2f}m\u00b2</text>'''
        hc_svg = '''
        <text x="350" y="130" style="font-size: 22px;font-family:arial">Perimeter: {:.2f}%</text>
        <text x="350" y="165" style="font-size: 22px;font-family:arial">Perimeter {} 90%: {}</text>
        '''.format(100 * svarea/totarea, comps[3], ('Fail', 'Pass')[svarea/totarea >= 0.9]) if node.leed_menu else ""

        svg_str = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
        <svg
        id="svg5"
        version="1.1"
        viewBox="0 0 {0[0]} {0[1]}"
        width="{0[0]}"
        height="{0[1]}"
        xmlns="http://www.w3.org/2000/svg"
        xmlns:svg="http://www.w3.org/2000/svg">

        <defs
            id="defs111">
            <linearGradient
                id="linearGradient1">
                <stop
                    style="stop-color:#e2e2e2;stop-opacity:1;"
                    offset="0"
                    id="stop4186" />
                <stop
                    style="stop-color:#ffffff;stop-opacity:1;"
                    offset="1"
                    id="stop4188" />
            </linearGradient>
            <linearGradient
                xlink:href="#linearGradient1"
                id="linearGradient2"
                x1="0"
                y1="400"
                x2="800"
                y2="400"
                gradientUnits="userSpaceOnUse" />
        </defs>

        <rect style="fill:url(#linearGradient2);stroke:rgb(0,0,0);stroke-width:0" width="{0[0]}" height="{0[1]}"/>
        <text text-anchor="middle" x="300" y="50" style="font-size: 32px;font-family:arial">LEED v4 (Option 1) Daylighting Analysis</text>
        <text x="25" y="100" style="font-size: 26px;font-family:arial">Sensor location: {1}</text>
        <text x="25" y="130" style="font-size: 26px;font-family:arial">Sensor area: {2:.2f}m\u00b2</text>{3}
        <text x="350" y="270" style="font-size: 28px;font-family:arial">ASE {12[0]} {4}%: {5}</text>
        <text x="350" y="300" style="font-size: 24px;font-family:arial">({13} Points achievable)</text>
        <text x="350" y="555" style="font-size: 28px;font-family:arial">sDA {12[1]} {6}%: {7}</text>
        <text x="410" y="585" style="font-size: 24px;font-family:arial">({8} Points)</text>
        <text x="350" y="645" style="font-size: 28px;font-family:arial">sDA {12[2]} {9}%: {10}</text>
        <text x="410" y="675" style="font-size: 24px;font-family:arial">({11} Points)</text>
        <text text-anchor="middle" x="175" y="750" style="font-size: 20px;font-family:arial">Spatial Daylight Autonomy (300/50%)</text>
        <text text-anchor="middle" x="175" y="425" style="font-size: 20px;font-family:arial">Annual Sunlight Exposure (1000/250)</text>
        """.format(dim, node.zone_menu, totarea, hc_svg, asepass, ('Pass', 'Fail')[ase >= asepass], sdapass[0], ('Fail', 'Pass')[sda >= sdapass[0]],
                   sda1points, sdapass[1], ('Fail', 'Pass')[sda >= sdapass[1]], sda2points, comps, asepoints)

        for b in range(20):
            bfill = "0, 0, 255" if (b + 1) * 5 <= sdapass[1] else "0, 255, 0"
            bfill = "255, 0, 0" if (b + 1) * 5 <= sdapass[0] else bfill
            alpha = 0.95 if -5 <= sda - ((b + 1) * 5) <= 0 else 0.45
            svg_str += '        <rect style="fill:rgb({})" fill-opacity="{}" stroke="rgb(0, 0, 0)" stroke-width="1" x="{}" y="{}" width="{}" height="{}"/>\n'.format(bfill,
                                                                                                                                                                     alpha,
                                                                                                                                                                     25 + int(b%4) * 75,
                                                                                                                                                                     675 - int(b/4) * 50,
                                                                                                                                                                     75,
                                                                                                                                                                     50)

            if alpha == 0.95:
                svg_str += '        <text text-anchor="middle" x="{}" y="{}" style="font-size:24px;font-family:arial">{:.1f}</text>'.format(65 + int(b % 4) * 75, 708 - int(b/4) * 50, sda)

        for b in range(20):
            bfill = "255, 0, 0" if (b + 1) * 5 > asepass else "0, 255, 0"
            alpha = 1.0 if -5 <= ase - ((b + 1) * 5) <= 0 else 0.45
            svg_str += '        <rect style="fill:rgb({})" fill-opacity="{}" stroke="rgb(0, 0, 0)" stroke-width="1" x="{}" y="{}" width="{}" height="{}"/>\n'.format(bfill,
                                                                                                                                                                     alpha, 25 + int(b%4) * 75,
                                                                                                                                                                     350 - int(b/4) * 50, 75, 50)

            if alpha == 1.0:
                svg_str += '        <text text-anchor="middle" x="{}" y="{}" style="font-size:24px;font-family:arial">{:.1f}</text>'.format(65 + int(b % 4) * 75, 383 - int(b/4) * 50, ase)

        svg_str += '        <text x="360" y="780" style="font-size:26px;font-family:arial">Total points: {} of {}</text>'.format(credits, tcredits)

        svg_str += "</svg>"

        with open(os.path.join(svp['viparams']['newdir'], 'images', 'LEED_{}_light.svg'.format(node.zone_menu)), 'w') as svg_file:
            svg_file.write(svg_str)

        return imname, bytearray(svg_str, encoding='utf-8')

    elif node.metric == '0' and node.energy_menu == '0':
        pass

    elif node.metric == '0' and node.energy_menu == '1':
        for met in ('Heating (W)', 'Cooling (W)', 'Temperatature (C)', 'CO2 (ppm)', 'Ventilation (heat (W)', 'Fabric heat (W)', 'Power (W)'):
            pass

    elif node.metric == '6':
        if node.frame_menu != 'All':
            imname = "Whole_life_carbon"
            wlc = kwargs['wlc']
            ec = kwargs['ec']
            noc = kwargs['noc']
            min_res = min((noc[0], ec[0], wlc[0]))
            max_res = max((noc[0], ec[0], wlc[0]))
            min_res = round(min_res, -2)
            max_res = round(max_res, -2)

            if min_res < 0:
                min_res -= 100

            l_interval = int(round((max_res - min_res)/7, -2))
            l_range = range(int(min_res), int(max_res) + l_interval, l_interval)
            l_diff = 400/len(l_range)
            svg_str = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
            <svg
            id="svg5"
            version="1.1"
            viewBox="0 0 {0[0]} {0[0]}"
            width="{0[0]}"
            height="{0[1]}"
            xmlns="http://www.w3.org/2000/svg"
            xmlns:svg="http://www.w3.org/2000/svg">\n""".format(dim)
            svg_str += '''<defs>
            <radialGradient id="WLC_grad" gradientUnits="userSpaceOnUse"
            cx="500" cy="500" r="300" fx="500" fy="500">
            <stop offset="0%" stop-color="rgb(64, 255, 64)" stop-opacity="0.75" />
            <stop offset="50%" stop-color="rgb(64, 64, 255)"  stop-opacity="0.75"/>
            <stop offset="100%" stop-color="rgb(255, 64, 64)"  stop-opacity="0.75"/>
            </radialGradient>
            </defs>'''
            svg_str += '<rect x="0" y="0" width="{0[0]}" height="1000" style="fill:white;stroke:black;stroke-width:5"/>\n'.format(dim)

            for ii, i in enumerate(l_range):
                svg_str += '<rect x="{0}" y="{1}" width="{2}" height="{3}" style="fill:none;stroke:rgb({4}, {4}, {4});stroke-width:2"/>\n'.format(450 - l_diff*(ii), 450 - l_diff*(ii),
                                                                                                                                              100 + 2 * l_diff*(ii), 100 + 2 * l_diff*(ii),
                                                                                                                                              (200 - 20 * ii))
                svg_str += '<rect x="460" y="{}" width="80" height="20" style="fill:white;stroke:none"/>\n'.format(540 + l_diff*(ii))
                svg_str += '<text x="500" y="{}" text-anchor="middle" style="font-size:28px;font-family:arial">{}</text>\n'.format(560 + l_diff*(ii), i)

            svg_str += '<text x="500" y="75" text-anchor="middle" style="font-size:48px;font-family:arial">Whole-life Carbon</text>\n'
            svg_str += '<text x="75" y="500" transform="rotate(-90,75,500)" text-anchor="middle" style="font-size:48px;font-family:arial">Net operational carbon</text>\n'
            svg_str += '<text x="925" y="500" transform="rotate(90,925,500)" text-anchor="middle" style="font-size:48px;font-family:arial">Embodied carbon</text>\n'
            svg_str += '<text x="500" y="975" text-anchor="middle" style="font-size:48px;font-family:arial">kgCO2e</text>\n'
            wlc_pos = 450 - ((wlc[0] - min_res)/(max_res - min_res)) * (900-550)
            noc_pos = 450 - ((noc[0] - min_res)/(max_res - min_res)) * (900-550)
            ec_pos = 550 + ((ec[0] - min_res)/(max_res - min_res)) * (900-550)
            svg_str += '<polygon points="{},{} {},{} {},{}" style="fill:url(#WLC_grad);stroke:black;stroke-width:2"/>\n'.format(noc_pos, 500, 500, wlc_pos, ec_pos, 500)
            svg_str += '<circle cx="500" cy="{0}" r="15" style="fill:yellow;stroke:black;stroke-width:2"/>\n'.format(wlc_pos)
            svg_str += '<circle cx="{}" cy="500" r="15" style="fill:yellow;stroke:black;stroke-width:2"/>\n'.format(noc_pos)
            svg_str += '<circle cx="{}" cy="500" r="15" style="fill:yellow;stroke:black;stroke-width:2"/>\n'.format(ec_pos)
            svg_str += "</svg>"

            return imname, bytearray(svg_str, encoding='utf-8')
        else:
            imname = f"Whole_life_carbon_{node.zone_menu}"
            svg_str = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
            <svg
            id="svg5"
            version="1.1"
            viewBox="0 0 {0[0]} {0[1]}"
            width="{0[0]}"
            height="{0[1]}"
            xmlns="http://www.w3.org/2000/svg"
            xmlns:svg="http://www.w3.org/2000/svg">\n""".format(dim)

            max_val = max([max(res) for res in (kwargs['wlc'], kwargs['noc'], kwargs['ec'])])
            min_val = min([min(res) for res in (kwargs['wlc'], kwargs['noc'], kwargs['ec'])])
            res_new = (dim[1] - 250)/(max_val-min_val)
            res_yscale = 275/max((abs(max_val), abs(min_val)))
            valsay = [min_val + y * round(max_val - min_val, -2) * 0.2 for y in range(6)]
            res_xscale = (dim[0] - 220)/(len(kwargs['wlc']) - 1)
            cols = ('#440154', '#20A387', '#FDE725')
            pointsx = [120 + r * res_xscale for r in range(len(kwargs['wlc']))]
            svg_str += '<rect width="{0[0]}" height="{0[1]}" style="fill:white;stroke:none;stroke-width:5"/>\n'.format(dim)

            if min_val <= 0:
                svg_str += '<line x1="120" y1="{0:.3f}" x2="900" y2="{0:.3f}" style="stroke:green;stroke-width:4"/>\n'.format(650 + min_val * res_new)

            svg_str += '<rect x="120" y="100" width="780" height="550" style="fill:none;stroke:grey;stroke-width:4"/>\n'.format(650 + min_val * res_new)
            svg_str += '<text x="30" y="375" text-anchor="middle" transform="rotate(-90,30,375)" style="font-size:36px;font-family:arial">kgCO\u2082e</text>\n'.format(650 + min_val * res_new)
            svg_str += '<text x="{}" y="725" text-anchor="middle" style="font-size:36px;font-family:arial">Scenario</text>\n'.format(dim[0] * 0.5)
            svg_str += '<text x="{}" y="40" text-anchor="middle" style="font-size:38px;font-family:arial">Whole-life Carbon</text>\n'.format(dim[0] * 0.5)
            svg_str += '<text x="{}" y="80" text-anchor="middle" style="font-size:36px;font-family:arial">Zone: {}</text>\n'.format(dim[0] * 0.5, node.zone_menu)

            for ipx, pointx in enumerate(pointsx):
                svg_str += '<line x1="{0}" y1="100" x2="{0}" y2="660" style="stroke:grey;stroke-width:2"/>\n'.format(pointx)
                svg_str += '<text x="{}" y="690" text-anchor="middle" style="font-size:30px;font-family:arial">{}</text>\n'.format(pointx, ipx + 1)

            for ivy, valy in enumerate(valsay):
                svg_str += '<line x1="120" y1="{0}" x2="900" y2="{0}" style="stroke:grey;stroke-width:2"/>\n'.format((dim[1] - 150) - (valy - min_val) * res_new)
                svg_str += '<text x="75" y="{}" text-anchor="middle" style="font-size:30px;font-family:arial">{:.0f}</text>\n'.format((dim[1] - 150) - (valy - min_val) * res_new + 10, valy)

            for ri, res in enumerate((kwargs['wlc'], kwargs['noc'], kwargs['ec'])):
                pointsy = [(dim[1] - 150) - (r - min_val) * res_new for r in res]
                pointsz = list(zip(pointsx, pointsy))
                points = ' '.join(['{:.2f},{:.2f}'.format(p[0], p[1]) for p in pointsz])
                svg_str += '<polyline points="{}" style="fill:none;stroke:{};stroke-width:5"/>\n'.format(points, cols[ri])
                svg_str += '<rect x="{}" y="740" width="50" height="40" style="fill:{};stroke:black;stroke-width:2"/>\n'.format(50 + ri * 225, cols[ri])
                svg_str += '<text x="{}" y="770" style="font-size:30px;font-family:arial">{}</text>'.format(110 + ri * 225, ('Whole-life', 'Operational', 'Embodied')[ri])
                # svg_str += '<text x="700" y="770" style="font-size:30px;font-family:arial">Zone: {}</text>'.format(node.zone_menu)

                for point in pointsz:
                    svg_str += '<circle cx="{0[0]}" cy="{0[1]}" r="10" style="fill:{1};stroke:black;stroke-width:1"/>\n'.format(point, cols[ri])

            svg_str += "</svg>"
            return imname, bytearray(svg_str, encoding='utf-8')
