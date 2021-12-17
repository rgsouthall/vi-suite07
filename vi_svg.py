def vi_info(node, dim, **kwargs):
    for key, value in kwargs.items():
        print("{} = {}".format(key, value))
    if node.metric == '1' and node.light_menu == '2':
        imname = "RIBA_lighting"
        svg_str = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
        <svg
        id="svg5"
        version="1.1"
        viewBox="0 0 800 800"
        height="800"
        width="800"
        xmlns="http://www.w3.org/2000/svg"
        xmlns:svg="http://www.w3.org/2000/svg">
        <defs
        id="defs2" />
        <rect style="fill:rgb(255, 255, 255)" width="800" height="800"/>
        <path style="fill:rgb(255, 0, 0)" d="M 400 700
            A 300 300, 1, 0, 1, 100 400
            L 400 400 Z"/>
        <circle  style="fill:rgb(255, 255, 255)" cx="400" cy="400" r="260"/>
        <path style="fill:rgb(0, 255, 0)" d="M 400 650
            A 250 250, 1, 1, 1, 650 400
            L 400 400 Z"/>
        <circle  style="fill:rgb(255, 255, 255)" cx="400" cy="400" r="210"/>
        </svg>
        """.format(dim)

        return imname, bytearray(svg_str, encoding='utf-8')