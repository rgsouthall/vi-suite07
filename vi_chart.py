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

import bpy
import sys
from .envi_func import retmenu
from numpy import amax, amin, linspace, uint8, sin, cos, sign, deg2rad


def label(dnode, metric, axis, variant):
    catdict = {'clim': 'Ambient', 'zone spatial': 'Zone spatial', 'zone temporal': 'Zone temporal',
               'Embodied carbon': 'Embodied carbon', 'Linkage': 'Linkage', 'External': 'External', 'Frames': 'Frame',
               'metric': dnode.inputs[axis].resultmenu + ' metric', 'type': metric}
    animdict = {'metric': dnode.inputs[axis].resultmenu, 'type': metric}

    if dnode.inputs['X-axis'].resultmenu == 'Frames':
        return animdict[variant]
    else:
        return catdict[variant]


def llabel(dnode, metric, axis, variant):
    rdict = {'Climate': 'Ambient',
             'Zone spatial': dnode.inputs[axis].zonemenu,
             'Zone temporal': dnode.inputs[axis].zonemenu,
             'Embodied carbon': dnode.inputs[axis].zonemenu,
             'Carbon': dnode.inputs[axis].zonemenu,
             'Linkage': dnode.inputs[axis].zonemenu,
             'Frames': 'Frames', 'Camera': dnode.inputs[axis].zonemenu,
             'Position': dnode.inputs[axis].zonemenu,
             'External': dnode.inputs[axis].zonemenu,
             'Power': dnode.inputs[axis].zonemenu,
             'Probe': dnode.inputs[axis].zonemenu,
             'Surface': dnode.inputs[axis].zonemenu,
             'Time': dnode.inputs[axis].zonemenu}
    ldict = {'type': rdict[dnode.inputs[axis].resultmenu], 'metric': metric, }
    return ldict[variant]


def statdata(res, stat):
    if stat == 'Average':
        return ([sum(r)/len(r) for r in res])
    elif stat == 'Maximum':
        return ([max(r) for r in res])
    elif stat == 'Minimum':
        return ([min(r) for r in res])
    elif stat == 'Sum':
        return ([sum(r) for r in res])


def rvariant(dnode):
    axes = ('Y-axis 1', 'Y-axis 2', 'Y-axis 3')
    zones = [dnode.inputs[axis].zonemenu for axis in axes if dnode.inputs[axis].links and dnode.inputs[axis].resultmenu == 'Zone spatial']
    clims = [dnode.inputs[axis].metricmenu for axis in axes if dnode.inputs[axis].links and dnode.inputs[axis].resultmenu == 'Climate']

    if zones and len(set(zones)) + len(set(clims)) == len(zones + clims):
        return 'type'
    else:
        return 'metric'


def timedata(datastring, timetype, stattype, months, days, dos, dnode, Sdate, Edate):
    if dnode.inputs['X-axis'].resultmenu != 'Time' or timetype == '0':
        return datastring
    else:
        if timetype == '1' and dnode.inputs['X-axis'].resultmenu == 'Time':
            res = [[] for d in range(len(set(dos)))]

            for h, val in enumerate(datastring):
                res[dos[h] - dos[0]].append(val)

        elif timetype == '2' and dnode.inputs['X-axis'].resultmenu == 'Time':
            res = [[] for m in range(len(set(months)))]

            for h, val in enumerate(datastring):
                res[months[h] - months[0]].append(val)

        return (statdata(res, stattype))


def retframe(axis, dnode, frames):
    if len(set(frames)) > 1 and dnode.inputs['X-axis'].resultmenu == 'Frames':
        return 'All'
    elif len(set(frames)) > 1:
        return dnode.inputs[axis].framemenu
    else:
        return frames[0]


def chart_disp(chart_op, plt, dnode, rnodes, Sdate, Edate):
    plt.close('all')
    plt.style.use('bmh')
    fig, ax = plt.subplots(dpi=dnode.dpi)
    variant = rvariant(dnode)
    rsx = dnode.inputs['X-axis']
    rnx = rsx.links[0].from_node
    rlx = rnx['reslists']
    rzlx = list(zip(*rlx))
    mdata = [rx[4].split() for rx in rlx if rx[0] == rsx.framemenu and rx[1] == 'Time' and rx[2] == 'Time' and rx[3] == 'Month']
    ddata = [rx[4].split() for rx in rlx if rx[0] == rsx.framemenu and rx[1] == 'Time' and rx[2] == 'Time' and rx[3] == 'Day']
    sdata = [rx[4].split() for rx in rlx if rx[0] == rsx.framemenu and rx[1] == 'Time' and rx[2] == 'Time' and rx[3] == 'DOS']
    hdata = [rx[4].split() for rx in rlx if rx[0] == rsx.framemenu and rx[1] == 'Time' and rx[2] == 'Time' and rx[3] == 'Hour']
    tdata = [rx[4].split() for rx in rlx if rx[0] == rsx.framemenu and rx[1] == 'Time' and rx[2] == 'Time' and rx[3] == 'Steps']

    if rsx.resultmenu != 'Time':
        si, ei = dnode["Start"] - dnode.id_properties_ui("Start").as_dict()['min'], dnode["End"] - dnode.id_properties_ui("Start").as_dict()['min']
    else:
        sm, sd, em, ed = Sdate.month, Sdate.day, Edate.month, Edate.day

        if mdata:
            (dm, dd) = ([int(x) for x in mdata[0]], [int(x) for x in ddata[0]])

            for i in range(len(hdata[0])):
                if sm == dm[i] and sd == dd[i]:  # and sh == dh[i] - 1:
                    si = i
                    break

            for i in range(len(hdata[0])):
                ei = i

                if em == dm[i] and ed == dd[i]:
                    ei += 23
                    break

            mdata = [int(m) for m in mdata[0]][si:ei + 1]
            ddata = [int(d) for d in ddata[0]][si:ei + 1]
            sdata = [int(s) for s in sdata[0]][si:ei + 1]
        else:
            si, ei = 0, -2

    linestyles = ('solid', '--', ':')
    colors = ('k', 'k', 'k')

    if rsx.resultmenu == 'Time':
        if dnode.timemenu == '0':
            xdata = range(1, ei-si + 2)
            xlabel = 'Time (hours)'
        elif dnode.timemenu == '1':
            xdata = range(dnode['Start'], dnode['End'] + 1)
            xlabel = 'Time (day of year)'
        elif dnode.timemenu == '2':
            xdata = range(Sdate.month, Edate.month + 1)
            xlabel = 'Time (months)'
        if rnx.bl_label == 'FloVi Simulation':
            xdata = [float(s) for s in tdata[0]][si:ei]
            xlabel = 'False time'
    else:
        # menus = retmenu(dnode, 'X-axis', rsx.resultmenu)
        # print(rsx.framemenu, rsx.resultmenu, rsx.zonemenu, rsx.metricmenu)
        data = [rx[4].split()[si:ei + 1] for rx in rlx if rx[0] == rsx.framemenu and rx[1] == rsx.resultmenu and rx[2] == rsx.zonemenu and rx[3] == rsx.metricmenu][0]
        xdata = timedata([rsx.multfactor * float(xd) for xd in data], dnode.timemenu, rsx.statmenu, mdata, ddata, sdata, dnode, Sdate, Edate)
        xlabel = rsx.metricmenu

    rny1 = dnode.inputs['Y-axis 1'].links[0].from_node
    rly1 = rny1['reslists']
    rzly1 = list(zip(*rly1))
    framey1 = retframe('Y-axis 1', dnode, rzly1[0])
    #menusy1 = retmenu(dnode, 'Y-axis 1', dnode.inputs['Y-axis 1'].resultmenu)

    try:
        #print(framey1, dnode.inputs['Y-axis 1'].resultmenu, menusy1[0], menusy1[1])
        y1d = [ry1[4].split()[si:ei + 1] for ry1 in rly1 if ry1[0] == framey1 and ry1[1] == dnode.inputs['Y-axis 1'].resultmenu and ry1[2] == dnode.inputs['Y-axis 1'].zonemenu and ry1[3] == dnode.inputs['Y-axis 1'].metricmenu][0]
    except Exception as e:
        chart_op.report({'ERROR'}, 'Invalid data on the y1 axis: {}'.format(e))
        return

    y1data = timedata([dnode.inputs['Y-axis 1'].multfactor * float(y) for y in y1d], dnode.timemenu, dnode.inputs['Y-axis 1'].statmenu, mdata, ddata, sdata, dnode, Sdate, Edate)
    ylabel = label(dnode, dnode.inputs['Y-axis 1'].metricmenu, 'Y-axis 1', variant)
    drange = checkdata(chart_op, xdata, y1data)
    line, = ax.plot(xdata[:drange], y1data[:drange], color=colors[0], ls=linestyles[0], linewidth=1, label=llabel(dnode, dnode.inputs['Y-axis 1'].metricmenu, 'Y-axis 1', variant))

    if dnode.inputs['Y-axis 2'].links:
        rny2 = dnode.inputs['Y-axis 2'].links[0].from_node
        rly2 = rny2['reslists']
        rzly2 = list(zip(*rly2))
        framey2 = retframe('Y-axis 2', dnode, rzly2[0])
        # menusy2 = retmenu(dnode, 'Y-axis 2', dnode.inputs['Y-axis 2'].resultmenu)

        try:
            y2d = [ry2[4].split()[si:ei + 1] for ry2 in rly2 if ry2[0] == framey2 and ry2[1] == dnode.inputs['Y-axis 2'].resultmenu and ry2[2] == dnode.inputs['Y-axis 2'].zonemenu and ry2[3] == dnode.inputs['Y-axis 2'].metricmenu][0]
        except Exception as e:
            chart_op.report({'ERROR'}, 'Invalid data on the y2 axis: {}'.format(e))
            return

        y2data = timedata([dnode.inputs['Y-axis 2'].multfactor * float(y) for y in y2d], dnode.timemenu, dnode.inputs['Y-axis 2'].statmenu, mdata, ddata, sdata, dnode, Sdate, Edate)
        drange = checkdata(chart_op, xdata, y2data)
        line, = ax.plot(xdata[:drange], y2data[:drange], color=colors[1], ls=linestyles[1], linewidth=1, label=llabel(dnode, dnode.inputs['Y-axis 2'].metricmenu, 'Y-axis 2', variant))

    if dnode.inputs['Y-axis 3'].links:
        rny3 = dnode.inputs['Y-axis 3'].links[0].from_node
        rly3 = rny3['reslists']
        rzly3 = list(zip(*rly3))
        framey3 = retframe('Y-axis 3', dnode, rzly3[0])
        # menusy3 = retmenu(dnode, 'Y-axis 3', dnode.inputs['Y-axis 3'].resultmenu)

        try:
            y3d = [ry3[4].split()[si:ei + 1] for ry3 in rly3 if ry3[0] == framey3 and ry3[1] == dnode.inputs['Y-axis 3'].resultmenu and ry3[2] == dnode.inputs['Y-axis 3'].zonemenu and ry3[3] == dnode.inputs['Y-axis 3'].metricmenu][0]
        except Exception as e:
            chart_op.report({'ERROR'}, 'Invalid data on the y3 axis: {}'.format(e))
            return

        y3data = timedata([dnode.inputs['Y-axis 3'].multfactor * float(y) for y in y3d], dnode.timemenu, dnode.inputs['Y-axis 3'].statmenu, mdata, ddata, sdata, dnode, Sdate, Edate)
        drange = checkdata(chart_op, xdata, y3data)
        line, = ax.plot(xdata[:drange], y3data[:drange], color=colors[2], ls=linestyles[2], linewidth=1, label=llabel(dnode, dnode.inputs['Y-axis 3'].metricmenu, 'Y-axis 3', variant))

    try:
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.legend()
        plt.grid(True)
        # plt.show(block=str(sys.platform) != 'win32')

        if sys.platform == 'darwin':
            plt.ion()

        plt.show()
    except Exception as e:
        chart_op.report({'ERROR'}, '{} Invalid data for this component'.format(e))

    def plot_graph(*args):
        args[0][0].plot()
        args[0][0].show()


def checkdata(chart_op, x, y):
    if len(x) != len(y):
        chart_op.report({'WARNING'}, 'X ({} points) and Y ({} points) data are not the same length'.format(len(x), len(y)))
        drange = min(len(x), len(y))
    else:
        drange = len(x)
    return drange


def hmchart_disp(chart_op, plt, dnode, col):
    x, y, z, var = dnode.x, dnode.y, dnode.z, dnode.metricmenu
    xmin = dnode.daystart if dnode.daystart > amin(x) else amin(x)
    xmax = dnode.dayend if dnode.dayend < amax(x) else amax(x)
    ymin = dnode.hourstart if dnode.hourstart > amin(y) else amin(y)
    ymax = dnode.hourend if dnode.hourend < amax(y) else amax(y)
    zmin = dnode.varmin if dnode.metricrange == '1' else amin(z)
    zmax = dnode.varmax if dnode.metricrange == '1' else amax(z)

    if zmax <= zmin and dnode.cf:
        chart_op.report({'ERROR'}, 'Maximum metric value is not greater than minimum')
        return

    plt.clf()
    plt.close()
    fig, ax = plt.subplots(figsize=(12, 6), dpi=dnode.dpi)
    plt.xlabel('Days', size=16)
    plt.ylabel('Hours', size=16)
    low_extend = 'lower' if zmin > amin(z) else 0
    up_extend = 'upper' if zmax < amax(z) else 0

    if all((low_extend, up_extend)):
        bar_extend = 'both'
    elif not any((low_extend, up_extend)):
        bar_extend = 'neither'
    else:
        bar_extend = ('', 'min')[low_extend == 'lower'] + ('', 'max')[up_extend == 'upper']

    # if dnode.cf:
    #     plt.contourf(x - 0.5, y, z, linspace(zmin, zmax, num=dnode.clevels + 1), levels=[zmin + (i + 1) * (zmax - zmin)/(dnode.clevels) for i in range(dnode.clevels - 1)], cmap=col, extend='both')
    # else:
    #     plt.pcolormesh(x, y, z, cmap=col, shading='auto', vmin=zmin, vmax=zmax, edgecolors='k', linewidths=0.075, snap=True, antialiased=True)

    # cbar = plt.colorbar(use_gridspec=True, pad=0.01, extend='neither')

    # if dnode.cl:
    #     try:
    #         ls = dnode.clevels + 1 if not dnode.lvals else [float(lev) for lev in dnode.lvals.split(" ")]
    #         cp = plt.contour(x - 0.5, y, z, linspace(zmin, zmax, num=dnode.clevels + 1), levels=ls, colors='Black', linewidths=dnode.lw)
    #         plt.clabel(cp, inline=True, fontsize=10)
    #     except Exception as e:
    #         print('except', linspace(zmin, zmax, num=dnode.clevels + 1))
    #         cp = plt.contour(x - 0.5, y, z, linspace(zmin, zmax, num=dnode.clevels + 1), levels=[zmin + i * (zmax - zmin)/(dnode.clevels) for i in range(dnode.clevels)], colors='Black', linewidths=dnode.lw)

    # if dnode.grid and dnode.cf:
    #     ax.grid(True, which='both', zorder=10)

    # cbar.set_label(label=var, size=16)
    # cbar.ax.tick_params(labelsize=14)

    # if dnode.inputs[0].links[0].from_node.bl_idname == 'No_Loc':
    #     plt.axis([xmin - 0.5, xmax + 0.5, ymin - 0.5, ymax + 0.5])
    # else:
    #     plt.axis([xmin - 0.5, xmax + 0.5, ymin - 0.5, ymax + 0.5])

    if dnode.cf:
        plt.contourf(x + 0.5, y + 0.5, z, linspace(zmin, zmax, num=dnode.clevels + 1),
                     levels=[zmin + i * (zmax - zmin)/(dnode.clevels) for i in range(dnode.clevels + 1)], cmap=col, extend=bar_extend)
        plt.axis([xmin + 0.5, xmax + 0.5, ymin + 0.5, ymax + 0.5])
    else:
        plt.axis([xmin - 0.5, xmax + 0.5, ymin - 0.5, ymax + 0.5])
        plt.pcolormesh(x, y, z, cmap=col, shading='auto', vmin=zmin, vmax=zmax, edgecolors='k', linewidths=0.075, snap=True, antialiased=True)

    cbar = plt.colorbar(use_gridspec=True, pad=0.01, extend=bar_extend)
    cbar.set_label(label=var, size=16)
    cbar.ax.tick_params(labelsize=14)

    if dnode.cl:
        try:
            ls = dnode.clevels + 1 if not dnode.lvals else [float(lev) for lev in dnode.lvals.split(" ")]
            cp = plt.contour(x, y, z, linspace(zmin, zmax, num=dnode.clevels), levels=ls, colors='Black', linewidths=dnode.lw)
            plt.clabel(cp, inline=True, fontsize=10)
        except Exception as e:
            cp = plt.contour(x + 0.5, y + 0.5, z, linspace(zmin, zmax, num=dnode.clevels), levels=[zmin + i * (zmax - zmin)/(dnode.clevels) for i in range(dnode.clevels + 1)][1:], colors='Black', linewidths=dnode.lw)

    if dnode.grid and dnode.cf:
        ax.grid(True, which='both', zorder=10)

    plt.xticks(size=14)
    plt.yticks(size=14)
    fig.tight_layout()
    
    if sys.platform == 'darwin':
        plt.ion()

    plt.show()


def ec_pie(chart_op, plt, node):
    plt.clf()
    plt.close()
    fig, ax = plt.subplots(figsize=(8, 6), subplot_kw=dict(aspect="equal"))
    labels = ['{}\n{:.1f} kgCO$_2$e'.format(k, node['res']['ec'][k]) for k in node['res']['ec'].keys() if k != 'All' and node['res']['ec'][k] > 0]
    values = [node['res']['ec'][k] for k in node['res']['ec'].keys() if k != 'All' and node['res']['ec'][k] > 0]
    wedge_properties = {"width": 0.3, "edgecolor": "w", 'linewidth': 2}
    cmap = plt.get_cmap('viridis')
    colors = [list(cmap(i)[:3]) + [0.7] for i in linspace(0, 1, len(values))]
    wedges, texts = ax.pie(values, wedgeprops=wedge_properties, startangle=0, shadow=False, colors=colors)
    bbox_props = dict(boxstyle="round,pad=0.3,rounding_size=0.1", fc="w", ec="grey", lw=0.72)
    kw = dict(arrowprops=dict(arrowstyle="-"), bbox=bbox_props, zorder=0, va="baseline")

    for i, p in enumerate(wedges):
        ang = (p.theta2 - p.theta1)/2. + p.theta1
        y = sin(deg2rad(ang))
        x = cos(deg2rad(ang))
        horizontalalignment = {-1: "right", 1: "left"}[int(sign(x))]
        connectionstyle = "arc, angleA=0, angleB={}, armA=20, armB=40, rad=5".format(ang)
        kw["arrowprops"].update({"connectionstyle": connectionstyle})

        ax.annotate(labels[i], xy=(x, y), xytext=(1.15*sign(x), 1.35*y),
                    horizontalalignment=horizontalalignment, verticalalignment="baseline", **kw)

    ax.annotate('Total kgCO$_2$e\n{:.1f}\nTotal kgCO$_2$e/m$^2$\n{:.1f}\nTotal kgCO$_2$e/m$^2$/y\n{:.1f}'.format(float(node['res']['ec']['All']), float(node['res']['ecm2']['All']), float(node['res']['ecm2y']['All'])),
                xy=(0, 0), xytext=(0, 0), horizontalalignment='center', va="center", size=14)
    
    if sys.platform == 'darwin':
        plt.ion()
    
    plt.show()


def wlc_line(chart_op, plt, node):
    plt.clf()
    plt.close()
    plt.style.use('bmh')
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    plt.subplots_adjust(bottom=0.2)
    data = [node['res']['wl'], node['res']['ec'], node['res']['noc']]
    xdata = [int(f[0]) for f in node['frames'] if f[0] != 'All']
    # cols = ('#440154', '#20A387', '#FDE725')
    cols = ('gold', 'turquoise', 'lime')
    lines = ('o-', '^-', 'd-')
    symbols = ('o', '^', 'd')
    labels = ('Whole-life', 'Embodied', 'Operational')

    for di, d in enumerate(data):
        ax.plot(xdata, list(d), lines[di], color=cols[di], alpha=0.5)
        ax.plot(xdata, list(d), symbols[di], color='k', markerfacecolor=(0, 0, 0, 0), markeredgecolor='grey')
        props = dict(boxstyle='round', facecolor=cols[di], alpha=0.5, edgecolor='k')
        plt.gcf().text(0.0125 * (di + 1) * 17, 0.05, labels[di], fontsize=14, bbox=props)

    ax.set_xticks(xdata)
    ax.set_title('Whole-life Carbon Analysis')
    plt.xlabel('Scenario')
    plt.ylabel('Carbon kgCO$_2$e')
    plt.grid(True)

    if sys.platform == 'darwin':
        plt.ion()

    plt.show()


def com_line(chart_op, plt, node):
    plt.clf()
    plt.close()
    plt.style.use('bmh')
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    plt.subplots_adjust(bottom=0.2)
    data = [node['res']['alloh1s'], node['res']['alloh2s']] if not node.occ or not node['res'].get('allooh1s') else [node['res']['allooh1s'], node['res']['allooh2s']]
    xdata = [int(f[0]) for f in node['frames'] if f[0] != 'All']
    # cols = ('#440154', '#20A387', '#FDE725')
    cols = ('orangered', 'gold')
    lines = ('o-', '^-', 'd-')
    symbols = ('o', '^', 'd')
    labels = ('25 - 28$^o$C', '> 28$^o$C')

    for di, d in enumerate(data):
        ax.plot(xdata, list(d), lines[di], color=cols[di], alpha=0.5)
        ax.plot(xdata, list(d), symbols[di], color='k', markerfacecolor=(0, 0, 0, 0), markeredgecolor='grey')
        props = dict(boxstyle='round', facecolor=cols[di], alpha=0.5, edgecolor='k')
        plt.gcf().text(0.0125 * (di + 1) * 17, 0.05, labels[di], fontsize=14, bbox=props)

    ax.set_xticks(xdata)
    ax.set_title('Over-heating Analysis')
    plt.xlabel('Scenario')
    occ = 'Occupied ' if node.occ and node['res'].get('allooh1s') else ''
    plt.ylabel(f'% {occ}hours')
    plt.grid(True)

    if sys.platform == 'darwin':
        plt.ion()

    plt.show()

# def ec_line(chart_op, plt, node):
#     plt.clf()
#     plt.close()
#     fig, ax = plt.subplots(figsize=(8, 6), subplot_kw=dict(aspect="equal"))

# def ec_wlc(chart_op, plt, node):