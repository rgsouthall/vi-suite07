from math import sin, cos, pi

def vi_info(node, dim, **kwargs):
    if node.metric == '1' and node.light_menu == '3':
        ir = kwargs['ir']
        aDF = kwargs['aDF']
        irscaled = ir if ir < 0.7 else 0.7
        adfscaled = aDF if aDF < 3.2 else 3.2
        adfxpos = 400 - (300 * sin(adfscaled*pi/2.0))
        adfypos = 400 + (300 * cos(adfscaled*pi/2.0))
        (adffill, adfsweep) = ("255, 0, 0", 0) if aDF < 2.0 else ("0, 255, 0", 1)
        irxpos = 400 - (250 * sin(irscaled*pi/0.4))
        irypos = 400 + (250 * cos(irscaled*pi/0.4))
        (irfill, irsweep) = ("255, 0, 0", 0) if ir < 0.4 else ("0, 255, 0", 1)


        imname = "RIBA_lighting_{}".format(node.zone_menu)
        svg_str = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
        <svg
        id="svg5"
        version="1.1"
        viewBox="0 0 {0} {0}"
        height="{0}"
        width="{0}"
        xmlns="http://www.w3.org/2000/svg"
        xmlns:svg="http://www.w3.org/2000/svg">
        <defs
        id="defs2" />
        <rect style="fill:rgb(255, 255, 255)" width="{0}" height="{0}"/>
        <path style="fill:rgb({1})" d="M 400 700
            A 300 300, 1, {2}, 1, {3:.0f} {4:.0f}
            L 400 400 Z"/>
        <circle  style="fill:rgb(255, 255, 255)" cx="400" cy="400" r="260"/>
        <path style="fill:rgb({5})" d="M 400 650
            A 250 250, 1, {6}, 1, {7:.0f} {8:.0f}
            L 400 400 Z"/>
        <circle  style="fill:rgb(255, 255, 255)" cx="400" cy="400" r="210"/>
        <text text-anchor="middle" x="400" y="400" style="font-size: 48px">RIBA 2030</text>
        <text text-anchor="middle" x="400" y="450" style="font-size: 48px">Lighting</text>
        <text x="20" y="770" style="font-size: 48px">Sensor: {9}</text>
        </svg>
        """.format(dim, adffill, adfsweep, adfxpos, adfypos, irfill, irsweep, irxpos, irypos, node.zone_menu)
        print(svg_str)
        return imname, bytearray(svg_str, encoding='utf-8')

    elif node.metric == '1' and node.light_menu == '2':
        ir = kwargs['ir']
        aDF = kwargs['aDF']
        adfpos = 700 - aDF * 150 if aDF < 4 else 100
        adfheight = 700 - adfpos
        irpos = 700 - ir * 750 if ir < 600/750 else 100
        irheight = 700 - irpos
        adffill = "255, 0, 0" if aDF < 2.0 else "0, 255, 0"
        irfill = "255, 0, 0" if ir < 0.4 else "0, 255, 0"
        imname = "RIBA_lighting_{}".format(node.zone_menu)
        svg_str = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
        <svg
        id="svg5"
        version="1.1"
        viewBox="0 0 {0} {0}"
        height="{0}"
        width="{0}"
        xmlns="http://www.w3.org/2000/svg"
        xmlns:svg="http://www.w3.org/2000/svg">
        <defs
        id="defs2" />
        <rect style="fill:rgb(255, 255, 255)" width="{0}" height="{0}"/>
        <rect style="fill:rgb({1})" ry="5" x="50" y="{2}" width="100" height="{3}"/>
        <rect style="fill:rgb({4})" ry="5" x="175" y="{5}" width="100" height="{6}"/>
        
        <text text-anchor="middle" x="400" y="50" style="font-size: 48px">RIBA 2030 Lighting</text>
        <text text-anchor="middle" x="100" y="{7}" style="font-size: 36px">{8}</text>
        <text text-anchor="middle" x="225" y="{9}" style="font-size: 36px">{10}</text>
        <text text-anchor="middle" x="100" y="750" style="font-size: 48px">aDF</text>
        <text text-anchor="middle" x="225" y="750" style="font-size: 48px">IR</text>
        <text x="275" y="750" style="font-size: 48px">Sensor: {11}</text>
        </svg>
        """.format(dim, adffill, adfpos, adfheight, irfill, irpos, irheight, adfpos - 25, aDF, irpos - 25, ir, node.zone_menu)

        print(svg_str)

        return imname, bytearray(svg_str, encoding='utf-8')

    elif node.metric == '1' and node.light_menu == '1':
        sda = kwargs['sda']
        sdapass = kwargs['sdapass']
        ase = kwargs['ase']
        asepass = kwargs['asepass']
        credits = kwargs['o1']
        imname = "LEED_lighting_{}".format(node.zone_menu)

        svg_str = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
        <svg
        id="svg5"
        version="1.1"
        viewBox="0 0 {0} {0}"
        height="{0}"
        width="{0}"
        xmlns="http://www.w3.org/2000/svg"
        xmlns:svg="http://www.w3.org/2000/svg">
        <rect style="fill:rgb(255, 255, 255)" width="{0}" height="{0}"/>
        <text text-anchor="middle" x="400" y="50" style="font-size: 32px">LEED v4 (Option 1) Daylighting Analysis</text>
        <text x="450" y="100" style="font-size: 20px">Sensor name: {1}</text>
        <text x="450" y="125" style="font-size: 20px">Floor area: {1}</text>
        <text text-anchor="middle" x="250" y="375" style="font-size: 24px">Spatial Daylight Autonomy</text>
        <text text-anchor="middle" x="250" y="725" style="font-size: 24px">Annual Sunlight Exposure</text>
        """.format(dim, node.zone_menu)

        for b in range(20):
            bfill = "255, 128, 128" if (b + 1) * 5 <= sdapass[0] else "128, 255, 128" 
            alpha = 1.0 if -5 <= sda - ((b + 1) * 5) < 0 else 0.25

            svg_str += '        <rect style="fill:rgb({})" fill-opacity="{}" stroke="rgb(0, 0, 0)" stroke_width="1" x="{}" y="{}" width="{}" height="{}"/>\n'.format(bfill, alpha, 100 + int(b%4) * 75, 300 - int(b/4) * 50, 75, 50)
            
            if alpha == 1.0:
                svg_str += '        <text text-anchor="middle" x="{}" y="{}" style="font-size: 24px">{:.1f}</text>'.format(137.5 + int(b%4) * 75, 333 - int(b/4) * 50, sda)

            
        
        for b in range(20):
            bfill = "255, 128, 128" if (b + 1) * 5 > asepass else "128, 255, 128" 
            alpha = 1.0 if -5 <= ase - ((b + 1) * 5) < 0 else 0.25
            svg_str += '        <rect style="fill:rgb({})" fill-opacity="{}" stroke="rgb(0, 0, 0)" stroke_width="1" x="{}" y="{}" width="{}" height="{}"/>\n'.format(bfill, alpha, 100 + int(b%4) * 75, 650 - int(b/4) * 50, 75, 50)
            
            if alpha == 1.0:
                svg_str += '        <text text-anchor="middle" x="{}" y="{}" style="font-size: 24px">{:.1f}</text>'.format(137.5 + int(b%4) * 75, 683 - int(b/4) * 50, ase)
        
        svg_str += '        <text text-anchor="middle" x="600" y="775" style="font-size: 32px">Credits: {}</text>'.format(credits)
        svg_str += "</svg>"

        print(svg_str)
        return imname, bytearray(svg_str, encoding='utf-8')