def vi_info(node, dim):
    imname = "test.png"
    svg_str = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
    <svg width="{0}" height="{0}" viewBox="0 0 {0} {0}" id="smile" version="1.1">
    <rect width="{0}" height="{0}" style="fill:rgb(255, 255, 255)"/>
        <path
            style="fill:rgb(0, 0, 255);shape-rendering:geometricPrecision"
            d="M 150,0 A 150,150 0 0 0 0,150 150,150 0 0 0 150,300 150,150 0 0 0 
                300,150 150,150 0 0 0 150,0 Z M 72,65 A 21,29.5 0 0 1 93,94.33 
                21,29.5 0 0 1 72,124 21,29.5 0 0 1 51,94.33 21,29.5 0 0 1 72,65 Z 
                m 156,0 a 21,29.5 0 0 1 21,29.5 21,29.5 0 0 1 -21,29.5 21,29.5 0 0 1 
                -21,-29.5 21,29.5 0 0 1 21,-29.5 z m -158.75,89.5 161.5,0 c 0,44.67 
                -36.125,80.75 -80.75,80.75 -44.67,0 -80.75,-36.125 -80.75,-80.75 z"
        />
    </svg>
    """.format(dim)

    return imname, bytearray(svg_str, encoding='utf-8')