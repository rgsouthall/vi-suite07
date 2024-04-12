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

import bpy, bmesh, os, datetime, shlex, sys, math, pickle, shutil, time
from math import sin, cos, pi, log10
from mathutils import Vector
from subprocess import Popen, PIPE, STDOUT
from numpy import array, where, in1d, transpose, savetxt, int8, float16, float32, float64, digitize, zeros, choose, inner, average, amax, amin, concatenate, logical_and, genfromtxt, logspace
from numpy import sum as nsum
from numpy import max as nmax
from numpy import min as nmin
from numpy import mean as nmean
from numpy import append as nappend
from .vi_func import vertarea, logentry, ct2RGB, clearlayers, chunks, selobj, solarPosition, sunapply
from .vi_dicts import res2unit, unit2res


def sunposlivi(scene, skynode, frames, sun, stime):
    svp = scene.vi_params

    if skynode['skynum'] < 3 or (skynode.skyprog == '1' and skynode.epsilon > 1):
        times = [stime + frame*datetime.timedelta(seconds=3600*skynode.interval) for frame in range(len(frames))]
        solposs = [solarPosition(t.timetuple()[7], t.hour + (t.minute)*0.016666, svp.latitude, svp.longitude) for t in times]
        beamvals = [(0, 3)[solposs[t][0] > 0] for t in range(len(times))] if skynode['skynum'] < 2  or (skynode.skyprog == '1' and skynode.epsilon > 1) else [0 for t in range(len(times))]
        skyvals = [1 for t in range(len(times))]

    elif skynode['skynum'] == 3 and skynode.skyprog == '0':
        times = [datetime.datetime(2015, 3, 20, 12, 0)]
        solposs = [solarPosition(t.timetuple()[7], t.hour + (t.minute)*0.016666, 0, 0) for t in times]
        beamvals = [0 for t in range(len(times))]
        skyvals = [1 for t in range(len(times))]

    shaddict = {'0': 0.01, '1': 1, '2': 1, '3': 1}
    values = list(zip([shaddict[str(skynode['skynum'])] for t in range(len(times))], beamvals, skyvals))
    sunapply(scene, sun, values, solposs, frames, skynode.sdist)


def face_bsdf(o, m, mname, f):
    if m.vi_params.get('bsdf'):
        uv = '{0[0]:.4f} {0[1]:.4f} {0[2]:.4f}'.format(m.vi_params.li_bsdf_up)
        #MGF geometry does not work (black inner face)
        # if m.vi_params.li_bsdf_proxy_depth and m.vi_params['bsdf'].get('proxied'):
        #     trans = f.calc_center_median() - Vector([float(p) for p in m.vi_params['bsdf']['pos'].split()])
        #     rot = f.normal.rotation_difference(Vector((0, 0, 1))).to_euler()
        #     rot_deg = [180*r/math.pi for r in (rot.x, rot.y, rot.z)]
        #     print(rot_deg)
        #     upv = Vector((0, 1, 0))
        #     upv.rotate(rot)
        #     rot2 = upv.rotation_difference(Vector((0, 0, 1))).to_euler()
        #     rot2_deg = [180*r/math.pi for r in (rot2.x, rot2.y, rot2.z)]
        #     print('rot2', upv, rot2_deg)
        #     rot3 = f.normal.rotation_difference(Vector([-float(p) for p in m.vi_params['bsdf']['normal'].split()])).to_euler()
        #     rot3_deg = [180*r/math.pi for r in (rot3.x, rot3.y, rot3.z)]
        #     print(f.normal, m.vi_params['bsdf']['normal'])
        #     print(rot3_deg)
        #     pkgbsdfrun = Popen(shlex.split("pkgBSDF -s {}".format(os.path.join(bpy.context.scene.vi_params['viparams']['newdir'], 'bsdfs', '{}.xml'.format(m.name)))), stdout = PIPE)
        #     xformrun = Popen(shlex.split("xform -rx {0[0]} -ry {0[1]} -rz {0[2]} -t {1[0]} {1[1]} {1[2]}".format(rot_deg, trans)), stdin = pkgbsdfrun.stdout, stdout = PIPE)
        #     radentry = ''.join([line.decode() for line in xformrun.stdout])
        #     radentry = radentry.replace('m_{}_f'.format(m.name), mname)
        #     radentry = radentry.replace(' {}.xml '.format(m.name), ' {} '.format(os.path.join(bpy.context.scene.vi_params['viparams']['newdir'], 'bsdfs', '{}.xml'.format(m.name))))
        # else:
        radentry = 'void BSDF {0}\n6 {3} "{1}" {2} .\n0\n0\n\n'.format(mname, m.vi_params['bsdf']['filepath'], uv, m.vi_params.li_bsdf_proxy_depth)
        return radentry
    else:
        return ''

def res_interpolate(scene, dp, o, ores, plt, offset):
    bm = bmesh.new()
    bm.from_object(o, dp)
    bmesh.ops.triangulate(bm, faces=bm.faces[:])
    ovp = o.vi_params
    svp = scene.vi_params

    if svp.vi_leg_levels < 3:
        svp.vi_leg_levels = 3

    if svp.li_disp_menu == 'aga1v':
        var = 'aga{}v'.format(svp.vi_views)
    elif svp.li_disp_menu == 'ago1v':
        var = 'ago{}v'.format(svp.vi_views)
    else:
        var = svp.li_disp_menu

    if bm.faces.layers.float.get('{}{}'.format(var, scene.frame_current)):
        geom = bm.faces
    elif bm.verts.layers.float.get('{}{}'.format(var, scene.frame_current)):
        geom = bm.verts
    else:
        logentry(f"No result data on {o.name}. Re-export LiVi Context and Geometry")
        return

    if svp.vi_leg_scale == '0':
        levels = [svp.vi_leg_min + i * (svp.vi_leg_max - svp.vi_leg_min)/svp.vi_leg_levels for i in range(1, svp.vi_leg_levels)]

    elif svp.vi_leg_scale == '1':
        slices = logspace(0, 2, svp.vi_leg_levels + 1, True)
        bins = array([1 - log10(i)/log10(svp.vi_leg_levels + 1) for i in range(1, svp.vi_leg_levels + 2)][::-1])
        levels = svp.vi_leg_min + (svp.vi_leg_max - svp.vi_leg_min) * bins[1:-1]

    poss = [v.co for v in bm.verts]
    xs = [p[0] for p in poss]
    ys = [p[1] for p in poss]
    tris = [[v.index for v in face.verts] for face in bm.faces]
    res_lay = geom.layers.float['{}{}'.format(var, scene.frame_current)]

    if isinstance(geom[:][0], bmesh.types.BMVert):
        ress = array([v[res_lay] for v in bm.verts])
    else:
        ress = array([sum([f[res_lay] for f in v.link_faces])/len(v.link_faces) for v in bm.verts])

    bm.free()
    CS = plt.tricontourf(xs, ys, tris, ress, levels=levels, extend="both")
    vi, vcos, eis, fis, mis, v_start = 0, [], [], [], [], 0 
    plt.gca().set_aspect('equal')

    for csi, ca in enumerate(CS.allsegs):
        for ci, c_cos in enumerate(ca):
            vcos += [o.matrix_world@Vector(c + [offset + (1, -1)[svp.vi_disp_pos == "1"] * 0.0001 * csi]) for c in c_cos.tolist()]
            fis.append([vi for vi in range(v_start, v_start + len(c_cos))])      
            mis.append(csi)
            v_start += len(c_cos)

    while ores.material_slots:
        bpy.ops.object.material_slot_remove()

    ores.data.clear_geometry()
    ores.data.from_pydata(vcos, [], fis)   
    ores.data.validate()
    ores.data.update(calc_edges=True)

    for fi, face in enumerate(ores.data.polygons):
        face.material_index = mis[fi]
    
    for matname in ['{}#{}'.format('vi-suite', i) for i in range(svp.vi_leg_levels)]:
        if bpy.data.materials[matname] not in ores.data.materials[:]:
            bpy.ops.object.material_slot_add()
            ores.material_slots[-1].material = bpy.data.materials[matname]


def res_direction(scene, o, ores, offset):
    svp = scene.vi_params
    views = len([attrib for attrib in o.data.attributes if 'aga' == attrib.name[:3]])
    vf = o.data.attributes[f'aga{views}v{scene.frame_current}'].domain
    v_cos = []
    f_vs = []
    osizel = 0.2
    osizeh = 0.5
    azi = o.vi_params['azi'] + (svp.vi_views - 1) * 360/views
    v_angles = [azi, azi + 120, azi + 240]
    angxs = [sin(ang*pi/180) for ang in v_angles]
    angys = [cos(ang*pi/180) for ang in v_angles]
    rps = o.data.polygons if vf == 'FACE' else o.data.vertices
    maxval = max([v.value for v in o.data.attributes[f'{svp.li_disp_menu}{svp.vi_frames}'].data])
    minval = min([v.value for v in o.data.attributes[f'{svp.li_disp_menu}{svp.vi_frames}'].data])

    for rpi, rp in enumerate(rps):
        val = o.data.attributes[f'{svp.li_disp_menu}{svp.vi_frames}'].data[rpi].value
        size = osizel + (osizeh - osizel) * (val - minval)/(maxval - minval)
        size = svp.vi_arrow_size
        pc = o.matrix_world@rp.center if vf == 'FACE' else o.matrix_world@rp.co
        v_cos.append([pc[0] + angxs[0] * size, pc[1] + angys[0] * size, pc[2] + offset])
        v_cos.append([pc[0] + angxs[1] * size * 0.5, pc[1] + angys[1] * size * 0.5, pc[2] + offset])
        v_cos.append([pc[0] - angxs[0] * size * 0.5, pc[1] - angys[0] * size * 0.5, pc[2] + offset])
        v_cos.append([pc[0] + angxs[2] * size * 0.5, pc[1] + angys[2] * size * 0.5, pc[2] + offset])
        f_vs.append([rpi*4, rpi*4+1, rpi*4+2, rpi*4+3])
    
    while ores.material_slots:
        bpy.ops.object.material_slot_remove()

    ores.data.clear_geometry()   
    ores.data.from_pydata(vertices=v_cos, edges=[], faces=f_vs)
    bm = bmesh.new()
    bm.from_mesh(ores.data)

    for attrib in o.data.attributes:
        if attrib.data_type == 'FLOAT':
            bm.faces.layers.float.new(attrib.name)
            res = bm.faces.layers.float[attrib.name]
            
            for fi, face in enumerate(bm.faces):
                face[res] = o.data.attributes[attrib.name].data[fi].value
                    
        elif attrib.data_type == 'INT':
            bm.faces.layers.int.new(attrib.name)
            res = bm.faces.layers.int[attrib.name]
            
            for fi, face in enumerate(bm.faces):
                face[res] = o.data.attributes[attrib.name].data[fi].value

    if svp.vi_leg_scale == '0':
        levels = [svp.vi_leg_min + i * (svp.vi_leg_max - svp.vi_leg_min)/svp.vi_leg_levels for i in range(1, svp.vi_leg_levels)]

    elif svp.vi_leg_scale == '1':
        slices = logspace(0, 2, svp.vi_leg_levels + 1, True)
        bins = array([1 - log10(i)/log10(svp.vi_leg_levels + 1) for i in range(1, svp.vi_leg_levels + 2)][::-1])
        levels = svp.vi_leg_min + (svp.vi_leg_max - svp.vi_leg_min) * bins[1:-1]

    bm.to_mesh(ores.data)
    bm.free()

    for matname in ['{}#{}'.format('vi-suite', i) for i in range(svp.vi_leg_levels)]:
        if bpy.data.materials[matname] not in ores.data.materials[:]:
            bpy.ops.object.material_slot_add()
            ores.material_slots[-1].material = bpy.data.materials[matname]


def rtpoints(self, bm, offset, cp, frame):
    geom = bm.verts if cp == '1' else bm.faces
    cindex = geom.layers.int['cindex']
    rt = geom.layers.string['rt{}'.format(frame)]

    for gp in geom:
        gp[cindex] = 0

    geom.ensure_lookup_table()
    resfaces = [face for face in bm.faces if face.material_index <= len(self.id_data.data.materials) and self.id_data.data.materials[face.material_index] and self.id_data.data.materials[face.material_index].vi_params.mattype == '1']
    self['cfaces'] = [face.index for face in resfaces]

    if cp == '0':
        gpoints = resfaces
        gpcos = [gp.calc_center_median() for gp in gpoints]
        self['cverts'], self['lisenseareas'][frame] = [], [f.calc_area() for f in gpoints]

    elif cp == '1':
        gis = sorted(set([item.index for sublist in [face.verts[:] for face in resfaces] for item in sublist]))
        gpoints = [geom[gi] for gi in gis]
        gpcos = [gp.co for gp in gpoints]
        self['cverts'], self['lisenseareas'][frame] = gp.index, [vertarea(bm, gp) for gp in gpoints]

    for g, gp in enumerate(gpoints):
        gp[rt] = '{0[0]:.4f} {0[1]:.4f} {0[2]:.4f} {1[0]:.4f} {1[1]:.4f} {1[2]:.4f}'.format([gpcos[g][i] + offset * gp.normal.normalized()[i] for i in range(3)], gp.normal[:]).encode('utf-8')
        gp[cindex] = g + 1

    self['rtpnum'] = g + 1


def setscenelivivals(scene):
    svp = scene.vi_params
    svp['liparams']['maxres'], svp['liparams']['minres'], svp['liparams']['avres'] = {}, {}, {}
    res = svp.li_disp_menu
    olist = [o for o in bpy.data.objects if o.vi_params.vi_type_string == 'LiVi Calc']

    for frame in range(svp['liparams']['fs'], svp['liparams']['fe'] + 1):
        svp['liparams']['maxres'][str(frame)] = max([o.vi_params['omax']['{}{}'.format(res, frame)] for o in olist])
        svp['liparams']['minres'][str(frame)] = min([o.vi_params['omin']['{}{}'.format(res, frame)] for o in olist])
        svp['liparams']['avres'][str(frame)] = sum([o.vi_params['oave']['{}{}'.format(res, frame)] for o in olist])/len([o.vi_params['oave']['{}{}'.format(res, frame)] for o in olist])

    svp.vi_leg_max = max(svp['liparams']['maxres'].values())
    svp.vi_leg_min = min(svp['liparams']['minres'].values())


def validradparams(params):
    valids = ('-ps', '-pt', '-pj', '-dj', '-ds', '-dt', '-dc', '-dr', '-dp', '-ss', '-st', '-sj', '-ab',
              '-av', '-aa',	'-ar', '-ad', '-as', '-lr', '-lw', '-u+')

    for p, param in enumerate(params.split()):
        if not p % 2 and (param not in valids):
            return 0
        elif p % 2:
            try:
                float(param)
            except Exception:
                return 0
    return 1


def ret_radentry(self, radname, mod):
    if self.radmatmenu == '8':
        radentry = ''
        # if self.get('bsdf'):
        #     radentry = 'void BSDF {0}\n6 {1} "{2}" {3[0]} {3[1]} {3[2]} .\n0\n0\n'.format(radname, self.li_bsdf_proxy_depth, self['bsdf']['filepath'], self.li_bsdf_up)
        # else:
        #     logentry(f'{self.id_data.name} has no BSDF data. A simple plastic material has been exported instead')
        #     radentry = '# dummy material\nvoid plastic {}\n0\n0\n5 0.8 0.8 0.8 0.1 0.1\n\n'.format(radname)
    elif self.radmatmenu == '9':
        radentry = bpy.data.texts[self.radfile].as_string()+'\n\n' if self.radfile in [t.name for t in bpy.data.texts] else '# dummy material\nvoid plastic {}\n0\n0\n5 0.8 0.8 0.8 0.1 0.1\n\n'.format(radname)
    else:
        if self.radtransmenu == '0':
            tn = self.radtrans
        else:
            tn = (((0.8402528435 + 0.0072522239 * self.radtransmit * self.radtransmit) ** 0.5) - 0.9166530661)/(0.0036261119 * self.radtransmit)
            tn = (tn, tn, tn)

        radentry = ('# ' + ('plastic', 'glass', 'dielectric', 'translucent', 'mirror', 'light', 'metal', 'antimatter', 'BSDF', 'custom')[int(self.radmatmenu)] + ' material\n' +
                    '{} {} {}\n'.format(mod, ('plastic', 'glass', 'dielectric', 'trans', 'mirror', 'light', 'metal', 'antimatter', 'BSDF', 'custom')[int(self.radmatmenu)], radname) +
                    {'0': '0\n0\n5 {0[0]:.3f} {0[1]:.3f} {0[2]:.3f} {1:.3f} {2:.3f}\n'.format(self.radcolour, self.radspec, self.radrough),
                     '1': '0\n0\n3 {0[0]:.3f} {0[1]:.3f} {0[2]:.3f}\n'.format(tn),
                     '2': '0\n0\n5 {0[0]:.3f} {0[1]:.3f} {0[2]:.3f} {1:.3f} 0\n'.format(self.radtrans, self.radior),
                     '3': '0\n0\n7 {0[0]:.3f} {0[1]:.3f} {0[2]:.3f} {1:.3f} {2:.3f} {3:.3f} {4:.3f}\n'.format(self.radcolour, self.radspec, self.radrough, self.radtransdiff, self.radtranspec),
                     '4': '0\n0\n3 {0[0]:.3f} {0[1]:.3f} {0[2]:.3f}\n'.format(self.radcolour),
                     '5': '0\n0\n3 {0[0]:.3f} {0[1]:.3f} {0[2]:.3f}\n'.format([c * self.radintensity for c in (self.radcolour, ct2RGB(self.radct))[self.radcolmenu == '1']]),
                     '6': '0\n0\n5 {0[0]:.3f} {0[1]:.3f} {0[2]:.3f} {1:.3f} {2:.3f}\n'.format(self.radcolour, self.radspec, self.radrough),
                     '7': '1 void\n0\n0\n', '8': '1 void\n0\n0\n', '9': '1 void\n0\n0\n'}[self.radmatmenu] + '\n')
    return radentry


def radmat(self, scene):
    svp = scene.vi_params
    radname = self.id_data.name.replace(" ", "_")
    radname = radname.replace(",", "")
    self['radname'] = radname
    dirt_entry = f'void brightfunc {radname}_dirt\n4 dirt dirt.cal -s {3.3 * self.li_dirt_spacing:.2f}\n0\n1 {self.li_dirt_level:.2f}\n\n' if self.li_dirt else ''
    mod = f'{radname}_dirt' if self.li_dirt else 'void'

    if self.mattype == '0' and self.radmatmenu in ('0', '1', '2', '3', '6') and any((self.li_tex, self.li_am, self.li_norm)):

        try:
            if self.li_tex and self.radtex:
                fd, fn = os.path.dirname(bpy.data.filepath), os.path.splitext(os.path.basename(bpy.data.filepath))[0]
                nd = os.path.join(fd, fn)
                svp['liparams']['texfilebase'] = os.path.join(nd, 'textures')
                teximage = self.li_tex
                teximageloc = os.path.join(svp['liparams']['texfilebase'], '{}.hdr'.format(radname))
                off = scene.render.image_settings.file_format
                scene.render.image_settings.file_format = 'HDR'
                teximage.save_render(teximageloc)
                scene.render.image_settings.file_format = off
                (w, h) = teximage.size
                ar = ('*{}'.format(w/h), '') if w >= h else ('', '*{}'.format(h/w))
                dirt_entry = f'void brightfunc {radname}_dirt\n4 dirt dirt.cal -s {3.3 * self.li_dirt_spacing:.2f}\n0\n1 {self.li_dirt_level:.2f}\n\n' if self.li_dirt else ''
                dirt_mod = f'{radname}_dirt' if self.li_dirt else 'void'
                radentries = ["{}{} colorpict {}_tex\n7 red green blue '{}' . frac(Lu){} frac(Lv){}\n0\n0\n\n".format(dirt_entry, dirt_mod, radname, teximageloc, ar[0], ar[1])]
                mod = '{}_tex'.format(radname)
                radentries.append(ret_radentry(self, radname, mod))

            else:
                mod = 'void'
                radentries = [dirt_entry, ret_radentry(self, radname, mod)]

            if self.li_am and self.radtex:
                t_mod = mod
                mod = f'{radname}_im'
                b_mod = 'void' if not self.li_tex_black else self.li_tex_black.name
                amim = self.li_am
                amloc = os.path.join(svp['liparams']['texfilebase'], '{}_am.hdr'.format(radname))
                off = scene.render.image_settings.file_format
                scene.render.image_settings.file_format = 'HDR'
                amim.save_render(amloc)
                scene.render.image_settings.file_format = off
                (w, h) = amim.size
                ar = ('*{}'.format(w/h), '') if w >= h else ('', '*{}'.format(h/w))
                radentries[1] = ret_radentry(self, mod, t_mod)
                radentries.append("# alpha mapped material\nvoid mixpict {0}\n7 {5} '{4}' grey '{1}' . frac(Lu){2} frac(Lv){3}\n0\n0\n\n".format(f'{radname}', amloc, ar[0], ar[1], b_mod, mod))
                mod = '{}'.format(radname)

            elif self.li_norm and self.radtex:
                t_mod = mod
                mod = '{}_norm'.format(radname)
                norm = self.li_norm
                (w, h) = norm.size
                ar = ('*{}'.format(w/h), '') if w >= h else ('', '*{}'.format(h/w))
                normpixels = zeros(norm.size[0] * norm.size[1] * 4, dtype='float32')
                norm.pixels.foreach_get(normpixels)
                header = '2\n0 1 {}\n0 1 {}\n'.format(norm.size[1], norm.size[0])
                xdat = -1 + 2 * normpixels[:][0::4].reshape(norm.size[0], norm.size[1])
                ydat = -1 + 2 * normpixels[:][1::4].reshape(norm.size[0], norm.size[1])
                savetxt(os.path.join(svp['liparams']['texfilebase'], '{}.ddx'.format(radname)), xdat, fmt='%.2f', header=header, comments='')
                savetxt(os.path.join(svp['liparams']['texfilebase'], '{}.ddy'.format(radname)), ydat, fmt='%.2f', header=header, comments='')
                radentries.append("{7} texdata {0}\n9 ddx ddy ddz '{1}.ddx' '{1}.ddy' '{1}.ddy' nm.cal frac(Lv){2} frac(Lu){3}\n0\n7 {4} {5[0]} {5[1]} {5[2]} {6[0]} {6[1]} {6[2]}\n\n".format(mod,
                                    os.path.join(svp['viparams']['newdir'], 'textures', radname), ar[1], ar[1], self.li_norm_strength, self.nu, self.nside, t_mod))
                radentries.append(ret_radentry(self, radname, mod))


        except Exception as e:
            logentry('Problem with texture export {}'.format(e))
            return ''
    else:
        radentries = [dirt_entry, ret_radentry(self, radname, mod)]

    self['radentry'] = ''.join(radentries)
    return self['radentry']


def cbdmmtx(self, scene, locnode, export_op):
    svp = scene.vi_params
    res = (1, 2, 4, 8)[self.cbdm_res - 1]
    os.chdir(svp['viparams']['newdir'])
    (csh, ceh) = (self.cbdm_start_hour, self.cbdm_end_hour + 1) if not self.ay or (self.cbanalysismenu == '2' and self.leed4) else (0, 24)
    (sdoy, edoy) = (self.sdoy, self.cbdm_edoy) if not self.ay else (1, 365)

    if self['epwbase'][1] in (".epw", ".EPW"):
        with open(locnode.weather, "r") as epwfile:
            epwlines = epwfile.readlines()
            self['epwyear'] = epwlines[8].split(",")[0]

        Popen(("epw2wea", locnode.weather, "{}.wea".format(os.path.join(svp['viparams']['newdir'], self['epwbase'][0])))).wait()

        with open("{}.wea".format(os.path.join(svp['viparams']['newdir'], self['epwbase'][0])), 'r') as weafile:
            weadata = weafile.readlines()

        with open("{}.wea".format(os.path.join(svp['viparams']['newdir'], self['epwbase'][0])), 'w') as weafile:
            for line in weadata:
                ls = line.split()

                if len(ls) != 5:
                    weafile.write(line)
                elif csh <= float(ls[2]) <= ceh and sdoy <= datetime.datetime(svp['year'], int(ls[0]), int(ls[1])).timetuple().tm_yday <= edoy and datetime.datetime(svp['year'], int(ls[0]), int(ls[1])).weekday() <= (6, 4)[self.weekdays]:
                    weafile.write(line)

        gdmcmd = ('gendaymtx -m {} {} "{}"'.format(res, ('-O0', '-O1')[self['watts']],
                  "{0}.wea".format(os.path.join(svp['viparams']['newdir'], self['epwbase'][0]))))
        gdmcmdns = ('gendaymtx -d -m {} {} "{}"'.format(res, ('-O0', '-O1')[self['watts']],
                    "{0}.wea".format(os.path.join(svp['viparams']['newdir'], self['epwbase'][0]))))

        with open("{}.mtx".format(os.path.join(svp['viparams']['newdir'], self['epwbase'][0])), 'w') as mtxfile:
            Popen(shlex.split(gdmcmd), stdout=mtxfile, stderr=STDOUT).communicate()

        with open("{}ns.mtx".format(os.path.join(svp['viparams']['newdir'], self['epwbase'][0])), 'w') as mtxfile:
            Popen(shlex.split(gdmcmdns), stdout=mtxfile, stderr=STDOUT).communicate()
        
        with open("{}-whitesky.oct".format(svp['viparams']['filebase']), 'w') as wsfile:
            oconvcmd = "oconv -w -"
            Popen(shlex.split(oconvcmd), stdin=PIPE, stdout=wsfile).communicate(input=self['whitesky'].encode(sys.getfilesystemencoding()))
        
        return ("{}.mtx".format(os.path.join(svp['viparams']['newdir'], self['epwbase'][0])), "{}ns.mtx".format(os.path.join(svp['viparams']['newdir'], self['epwbase'][0])))

    else:
        export_op.report({'ERROR'}, "Not a valid EPW file")
        return ('', '')


def cbdmhdr(node, scene, exp_op):
    patches = (146, 578, 2306, 9218)[node.cbdm_res - 1]
    cbdm_res = (146, 578, 0, 2306, 9218).index(patches) + 1
    svp = scene.vi_params
    svpnd = svp['viparams']['newdir']
    targethdr = os.path.join(svpnd, node['epwbase'][0]+"{}.hdr".format(('l', 'w')[node['watts']]))
    temphdr = os.path.join(svpnd, "temp.hdr")
    latlonghdr = os.path.join(svpnd, node['epwbase'][0]+"{}p.hdr".format(('l', 'w')[node['watts']]))
    skyentry = hdrsky(node.hdrname, '1', 0, 1000) if node.sourcemenu == '1' and node.cbanalysismenu == '0' else hdrsky(targethdr, '1', 0, 1000)

    if node.sourcemenu != '1' or node.cbanalysismenu == '2':
        mtxlines = open(node['mtxfile'], 'r').readlines()
        
        # for line in mtxlines:
        #     if line.split('=')[0] == 'NROWS':
        #         patches = int(line.split('=')[1])
        #         cbdm_res = (146, 578, 0, 2306).index(patches) + 1
        #     elif line.split('=')[0] == 'NCOLS':
        #         mtxhours = int(line.split('=')[1])

        #         if mtxhours != len(node.times):
        #             exp_op.report({'ERROR'}, "Outdated MTX file")
        #             node._valid = 0
        #             return

        vecvals, vals = mtx2vals(mtxlines, datetime.datetime(svp['year'], 1, 1).weekday(), node, node.times)
        
        pcombfiles = ''.join(["{} ".format(os.path.join(svpnd, 'ps{}.hdr'.format(i))) for i in range(patches)])
        vwcmd = 'vwrays -ff -x 600 -y 600 -vta -vp 0 0 0 -vd 0 1 0 -vu 0 0 1 -vh 360 -vv 360 -vo 0 -va 0 -vs 0 -vl 0'
        rcontribcmd = 'rcontrib -bn {} -fo -ab 0 -ad 1 -n {} -ffc -x 600 -y 600 -ld- -V+ -e MF:{} -f reinhart.cal -b rbin -o "{}" -m sky_glow "{}-whitesky.oct"'.format(patches, svp['viparams']['nproc'],
                                                                                                                                                                 cbdm_res,
                                                                                                                                                                 os.path.join(svpnd, 'p%d.hdr'),
                                                                                                                                                                 os.path.join(svpnd,
                                                                                                                                                                 svp['viparams']['filename']))

        vwrun = Popen(shlex.split(vwcmd), stdout=PIPE)
        rcrun = Popen(shlex.split(rcontribcmd), stderr=PIPE, stdin=vwrun.stdout)
        rcrun.wait()


        for line in rcrun.stderr:
            logentry('HDR generation error: {}'.format(line))

        for j in range(patches):
            with open(os.path.join(svpnd, "ps{}.hdr".format(j)), 'w') as psfile:
                Popen(shlex.split('pcomb -h -s {} "{}"'.format(vals[j], os.path.join(svpnd, 'p{}.hdr'.format(j)))), stdout=psfile).wait()

            if not j:
                shutil.copyfile(os.path.join(svpnd, 'ps0.hdr'), os.path.join(svpnd, 'temp.hdr'))
            else:
                with open(os.path.join(svpnd, "running.hdr"), 'w') as runhdr:
                    Popen(shlex.split('pcomb -h "{}" "{}"'.format(os.path.join(svpnd, 'temp.hdr'), os.path.join(svpnd, 'ps{}.hdr'.format(j)))), stdout=runhdr).wait()

                shutil.copyfile(os.path.join(svpnd, 'running.hdr'), os.path.join(svpnd, 'temp.hdr'))

        shutil.copyfile(os.path.join(svpnd, 'temp.hdr'), targethdr)
        os.remove(os.path.join(svpnd, 'temp.hdr'))
        os.remove(os.path.join(svpnd, 'running.hdr'))

        # with open(targethdr, 'w') as epwhdr:
        #     if sys.platform == 'win32':
        #         Popen("pcomb -h {}".format(pcombfiles), stdout=epwhdr).wait()
        #     else:
        #         Popen(shlex.split('pcomb -h {}'.format(pcombfiles)), stdout=epwhdr).wait()

        [os.remove(os.path.join(svpnd, 'p{}.hdr'.format(i))) for i in range(patches)]
        [os.remove(os.path.join(svpnd, 'ps{}.hdr'.format(i))) for i in range(patches)]
        node.hdrname = targethdr

        if node.hdr:
            with open('{}.oct'.format(os.path.join(svpnd, node['epwbase'][0])), 'w') as hdroct:
                Popen(shlex.split('oconv -w - '), stdin=PIPE, stdout=hdroct, stderr=STDOUT).communicate(input=skyentry.encode(sys.getfilesystemencoding()))

            cntrun = Popen('cnt 750 1500'.split(), stdout=PIPE)
            rccmd = 'rcalc -f "{}" -e XD=1500;YD=750;inXD=0.000666;inYD=0.001333'.format(os.path.join(svp.vipath, 'RadFiles', 'lib', 'latlong.cal'))
            logentry('Running rcalc: {}'.format(rccmd))
            rcalcrun = Popen(shlex.split(rccmd), stdin=cntrun.stdout, stdout=PIPE)

            with open(latlonghdr, 'w') as panohdr:
                rtcmd = 'rtrace -n {} -x 1500 -y 750 -fac "{}.oct"'.format(svp['viparams']['nproc'], os.path.join(svpnd, node['epwbase'][0]))
                logentry('Running rtrace: {}'.format(rtcmd))
                Popen(shlex.split(rtcmd), stdin=rcalcrun.stdout, stdout=panohdr)

    return skyentry


def mtx2vals(mtxlines, fwd, node, times):
    for m, mtxline in enumerate(mtxlines):
        if 'NROWS' in mtxline:
            patches = int(mtxline.split('=')[1])

        elif mtxline == '\n':
            startline = m + 1
            break

    tothours = len(times)
    hours = [t.hour for t in times]
    mtxlarray = array([0.333 * sum([float(lv) for lv in fval.split(" ")]) for fval in mtxlines[startline:] if fval != '\n'], dtype=float)
    mtxshapearray = mtxlarray.reshape(patches, int(len(mtxlarray)/patches))
    vals = nsum(mtxshapearray, axis=1)
    vvarray = transpose(mtxshapearray)
    vvlist = vvarray.tolist()
    vecvals = [[hours[x], (fwd+int(hours[x]/24)) % 7, *vvlist[x]] for x in range(tothours)]
    return (vecvals, vals)


def hdrsky(hdrfile, hdrmap, hdrangle, hdrradius):
    hdrangle = '1 {:.3f}'.format(hdrangle * math.pi/180) if hdrangle else '1 0'
    hdrfn = {'0': 'sphere2latlong', '1': 'sphere2angmap'}[hdrmap]
    return ("# Sky material\nvoid colorpict hdr_env\n7 red green blue '{}' {}.cal sb_u sb_v\n0\n{}\n\nhdr_env glow env_glow\n0\n0\n4 1 1 1 0\n\nenv_glow bubble sky\n0\n0\n4 0 0 0 {}\n\n".format(hdrfile, hdrfn, hdrangle, hdrradius))


def retpmap(node, frame, scene):
    svp = scene.vi_params
    pportmats = ' '.join([mat.name.replace(" ", "_") for mat in bpy.data.materials if mat.vi_params.pport and mat.vi_params.get('radentry')])
    ammats = ' '.join([mat.name.replace(" ", "_") for mat in bpy.data.materials if mat.vi_params.mattype == '1' and mat.vi_params.radmatmenu == '7' and mat.vi_params.get('radentry')])
    pportentry = ' '.join(['-apo {}'.format(ppm) for ppm in pportmats.split()]) if pportmats else ''
    amentry = '-aps {}'.format(ammats) if ammats else ''
    gpentry = '-apg "{}-{}.gpm" {}'.format(svp['viparams']['filebase'], frame, node.pmapgno) if node.pmapgno else ''
    cpentry = '-apc "{}-{}.cpm" {}'.format(svp['viparams']['filebase'], frame, node.pmapcno) if node.pmapcno else ''
    gpfileentry = '-ap "{}-{}.gpm" 50'.format(svp['viparams']['filebase'], frame) if node.pmapgno else ''
    cpfileentry = '-ap "{}-{}.cpm" 50'.format(svp['viparams']['filebase'], frame) if node.pmapcno else ''
    return amentry, pportentry, gpentry, cpentry, gpfileentry, cpfileentry


def retsv(self, scene, frame, rtframe, chunk, rt):
    svcmd = "rcontrib -w -I -n {} {} -m sky_glow {}-{}.oct ".format(scene.vi_params['viparams']['nproc'], '-ab 1 -ad 8192 -aa 0 -ar 512 -as 1024 -lw 0.0002 ',
                                                                    scene.vi_params['viparams']['filebase'], frame)
    rtrun = Popen(svcmd.split(), stdin=PIPE, stdout=PIPE, stderr=STDOUT, universal_newlines=True).communicate(input='\n'.join([c[rt].decode('utf-8') for c in chunk]))
    reslines = nsum(array([[float(rv) for rv in r.split('\t')[:3]] for r in rtrun[0].splitlines()[10:]]), axis=1)
    reslines[reslines > 0] = 1
    return reslines.astype(int8)


def basiccalcapply(self, scene, frames, rtcmds, simnode, curres, pfile):
    svp = scene.vi_params
    reslists = []
    ll = svp.vi_leg_levels
    increment = 1/ll
    bm = bmesh.new()
    bm.from_mesh(self.id_data.data)
    bm.transform(self.id_data.matrix_world)
    self['omax'], self['omin'], self['oave'], self['livires'] = {}, {}, {}, {}
    clearlayers(bm, 'f')
    geom = bm.verts if svp['liparams']['cp'] == '1' else bm.faces
    cindex = geom.layers.int['cindex']

    for f, frame in enumerate(frames):
        self['res{}'.format(frame)] = {}
        if svp['liparams']['unit'] == 'Lux':
            geom.layers.float.new('illu{}'.format(frame))
            geom.layers.float.new('virradm2{}'.format(frame))
            illures = geom.layers.float['illu{}'.format(frame)]
            virradm2res = geom.layers.float['virradm2{}'.format(frame)]

        elif svp['liparams']['unit'] == 'DF (%)':
            geom.layers.float.new('df{}'.format(frame))
            geom.layers.float.new('virradm2{}'.format(frame))
            dfres = geom.layers.float['df{}'.format(frame)]
            virradm2res = geom.layers.float['virradm2{}'.format(frame)]

        elif svp['liparams']['unit'] == 'W/m2 (f)':
            geom.layers.float.new('firrad{}'.format(frame))
            geom.layers.float.new('firradm2{}'.format(frame))
            firradres = geom.layers.float['firrad{}'.format(frame)]
            firradm2res = geom.layers.float['firradm2{}'.format(frame)]

        elif svp['liparams']['unit'] == 'W/m2':
            geom.layers.float.new('firrad{}'.format(frame))
            geom.layers.float.new('firradm2{}'.format(frame))
            firradres = geom.layers.float['firrad{}'.format(frame)]
            firradm2res = geom.layers.float['firradm2{}'.format(frame)]
            geom.layers.float.new('firradr{}'.format(frame))
            geom.layers.float.new('firradrm2{}'.format(frame))
            firradrres = geom.layers.float['firradr{}'.format(frame)]
            firradrm2res = geom.layers.float['firradrm2{}'.format(frame)]
            geom.layers.float.new('firradg{}'.format(frame))
            geom.layers.float.new('firradgm2{}'.format(frame))
            firradgres = geom.layers.float['firradg{}'.format(frame)]
            firradgm2res = geom.layers.float['firradgm2{}'.format(frame)]
            geom.layers.float.new('firradb{}'.format(frame))
            geom.layers.float.new('firradbm2{}'.format(frame))
            firradbres = geom.layers.float['firradb{}'.format(frame)]
            firradbm2res = geom.layers.float['firradbm2{}'.format(frame)]

        geom.layers.float.new('res{}'.format(frame))

        if geom.layers.string.get('rt{}'.format(frame)):
            rtframe = frame
        else:
            kints = [int(k[2:]) for k in geom.layers.string.keys()]
            rtframe = max(kints) if frame > max(kints) else min(kints)

        rt = geom.layers.string['rt{}'.format(rtframe)]
        logentry('Running rtrace: {}'.format(rtcmds[f]))

        for chunk in chunks([g for g in geom if g[rt]], int(svp['viparams']['nproc']) * 500):
            rtrun = Popen(shlex.split(rtcmds[f]), stdin=PIPE, stdout=PIPE, stderr=PIPE, universal_newlines=True).communicate(input='\n'.join([c[rt].decode('utf-8') for c in chunk]))

            if rtrun[1]:
                logentry('rtrun error: {}'.format(rtrun[1]))
                pfile.check('CANCELLED')
                bm.free()
                return 'CANCELLED'
            else:
                xyzirrad = array([[float(v) for v in sl.split('\t')[:3]] for sl in rtrun[0].splitlines()])

                if svp['liparams']['unit'] == 'W/m2 (f)':
                    firradm2 = nsum(xyzirrad * array([0.333, 0.333, 0.333]), axis=1)
                elif svp['liparams']['unit'] in ('Lux', 'DF (%)'):
                    virradm2 = nsum(xyzirrad * array([0.333, 0.333, 0.333]), axis=1)

                    if svp['liparams']['unit'] == 'Lux':
                        illu = nsum(xyzirrad * array([0.265, 0.67, 0.065]), axis=1) * 179
                    elif svp['liparams']['unit'] == 'DF (%)':
                        df = nsum(xyzirrad * array([0.265, 0.67, 0.065]), axis=1) * 1.79

                elif svp['liparams']['unit'] == 'W/m2':
                    firradm2 = nsum(xyzirrad * array([0.333, 0.333, 0.333]), axis=1)
                    firradrm2 = xyzirrad[:, 0] * 0.333
                    firradgm2 = xyzirrad[:, 1] * 0.333
                    firradbm2 = xyzirrad[:, 2] * 0.333

                for gi, gp in enumerate(chunk):
                    gparea = gp.calc_area() if svp['liparams']['cp'] == '0' else vertarea(bm, gp)

                    if svp['liparams']['unit'] == 'W/m2 (f)':
                        gp[firradm2res] = firradm2[gi].astype(float32)
                        gp[firradres] = (firradm2[gi] * gparea).astype(float32)

                    elif svp['liparams']['unit'] in ('Lux', 'DF (%)'):
                        gp[virradm2res] = virradm2[gi].astype(float32)

                        if svp['liparams']['unit'] == 'Lux':
                            gp[illures] = illu[gi].astype(float32)
                        elif svp['liparams']['unit'] == 'DF (%)':
                            gp[dfres] = df[gi].astype(float32)

                    elif svp['liparams']['unit'] == 'W/m2':
                        gp[firradm2res] = firradm2[gi].astype(float32)
                        gp[firradres] = (firradm2[gi] * gparea).astype(float32)
                        gp[firradrm2res] = firradrm2[gi].astype(float32)
                        gp[firradrres] = (firradrm2[gi] * gparea).astype(float32)
                        gp[firradgm2res] = firradgm2[gi].astype(float32)
                        gp[firradgres] = (firradgm2[gi] * gparea).astype(float32)
                        gp[firradbm2res] = firradbm2[gi].astype(float32)
                        gp[firradbres] = (firradbm2[gi] * gparea).astype(float32)

            curres += len(chunk)

            if pfile.check(curres) == 'CANCELLED':
                bm.free()
                return 'CANCELLED'

        if svp['liparams']['unit'] == 'W/m2 (f)':
            oirradm2 = array([g[firradm2res] for g in geom]).astype(float64)
            maxoirradm2, minoirradm2, aveoirradm2 = nmax(oirradm2), nmin(oirradm2), nmean(oirradm2)
            oirrad = array([g[firradres] for g in geom]).astype(float64)
            maxoirrad, minoirrad, aveoirrad = nmax(oirrad), nmin(oirrad), nmean(oirrad)

        elif svp['liparams']['unit'] in ('Lux', 'DF (%)'):
            ovirrad = array([g[virradm2res] for g in geom]).astype(float64)
            maxovirrad, minovirrad, aveovirrad = nmax(ovirrad), nmin(ovirrad), nmean(ovirrad)

            if svp['liparams']['unit'] == 'Lux':
                oillu = array([g[illures] for g in geom]).astype(float64)
                maxoillu, minoillu, aveoillu = nmax(oillu), nmin(oillu), nmean(oillu)
            elif svp['liparams']['unit'] == 'DF (%)':
                odf = array([g[dfres] for g in geom]).astype(float64)
                maxodf, minodf, aveodf = nmax(odf), nmin(odf), nmean(odf)

        elif svp['liparams']['unit'] == 'W/m2':
            oirradm2 = array([g[firradm2res] for g in geom]).astype(float64)
            maxoirradm2, minoirradm2, aveoirradm2 = nmax(oirradm2), nmin(oirradm2), nmean(oirradm2)
            oirrad = array([g[firradres] for g in geom]).astype(float64)
            maxoirrad, minoirrad, aveoirrad = nmax(oirrad), nmin(oirrad), nmean(oirrad)
            oirradrm2 = array([g[firradrm2res] for g in geom]).astype(float64)
            maxoirradrm2, minoirradrm2, aveoirradrm2 = nmax(oirradrm2), nmin(oirradrm2), nmean(oirradrm2)
            oirradr = array([g[firradrres] for g in geom]).astype(float64)
            maxoirradr, minoirradr, aveoirradr = nmax(oirradr), nmin(oirradr), nmean(oirradr)
            oirradgm2 = array([g[firradgm2res] for g in geom]).astype(float64)
            maxoirradgm2, minoirradgm2, aveoirradgm2 = nmax(oirradgm2), nmin(oirradgm2), nmean(oirradgm2)
            oirradg = array([g[firradgres] for g in geom]).astype(float64)
            maxoirradg, minoirradg, aveoirradg = nmax(oirradg), nmin(oirradg), nmean(oirradg)
            oirradbm2 = array([g[firradbm2res] for g in geom]).astype(float64)
            maxoirradbm2, minoirradbm2, aveoirradbm2 = nmax(oirradbm2), nmin(oirradbm2), nmean(oirradbm2)
            oirradb = array([g[firradbres] for g in geom]).astype(float64)
            maxoirradb, minoirradb, aveoirradb = nmax(oirradb), nmin(oirradb), nmean(oirradb)

        if svp['liparams']['unit'] == 'W/m2 (f)':
            self['omax']['firrad{}'.format(frame)] = maxoirrad
            self['oave']['firrad{}'.format(frame)] = aveoirrad
            self['omin']['firrad{}'.format(frame)] = minoirrad
            self['omax']['firradm2{}'.format(frame)] = maxoirradm2
            self['oave']['firradm2{}'.format(frame)] = aveoirradm2
            self['omin']['firradm2{}'.format(frame)] = minoirradm2

        elif svp['liparams']['unit'] in ('Lux', 'DF (%)'):
            self['omax']['virradm2{}'.format(frame)] = maxovirrad
            self['oave']['virradm2{}'.format(frame)] = aveovirrad
            self['omin']['virradm2{}'.format(frame)] = minovirrad

            if svp['liparams']['unit'] == 'Lux':
                self['omax']['illu{}'.format(frame)] = maxoillu
                self['oave']['illu{}'.format(frame)] = aveoillu
                self['omin']['illu{}'.format(frame)] = minoillu

            elif svp['liparams']['unit'] == 'DF (%)':
                self['omax']['df{}'.format(frame)] = maxodf
                self['omin']['df{}'.format(frame)] = minodf
                self['oave']['df{}'.format(frame)] = aveodf

        if svp['liparams']['unit'] == 'W/m2':
            self['omax']['firrad{}'.format(frame)] = maxoirrad
            self['oave']['firrad{}'.format(frame)] = aveoirrad
            self['omin']['firrad{}'.format(frame)] = minoirrad
            self['omax']['firradm2{}'.format(frame)] = maxoirradm2
            self['oave']['firradm2{}'.format(frame)] = aveoirradm2
            self['omin']['firradm2{}'.format(frame)] = minoirradm2
            self['omax']['firradr{}'.format(frame)] = maxoirradr
            self['oave']['firradr{}'.format(frame)] = aveoirradr
            self['omin']['firradr{}'.format(frame)] = minoirradr
            self['omax']['firradrm2{}'.format(frame)] = maxoirradrm2
            self['oave']['firradrm2{}'.format(frame)] = aveoirradrm2
            self['omin']['firradrm2{}'.format(frame)] = minoirradrm2
            self['omax']['firradg{}'.format(frame)] = maxoirradg
            self['oave']['firradg{}'.format(frame)] = aveoirradg
            self['omin']['firradg{}'.format(frame)] = minoirradg
            self['omax']['firradgm2{}'.format(frame)] = maxoirradgm2
            self['oave']['firradgm2{}'.format(frame)] = aveoirradgm2
            self['omin']['firradgm2{}'.format(frame)] = minoirradgm2
            self['omax']['firradb{}'.format(frame)] = maxoirradb
            self['oave']['firradb{}'.format(frame)] = aveoirradb
            self['omin']['firradb{}'.format(frame)] = minoirradb
            self['omax']['firradbm2{}'.format(frame)] = maxoirradbm2
            self['oave']['firradbm2{}'.format(frame)] = aveoirradbm2
            self['omin']['firradbm2{}'.format(frame)] = minoirradbm2

        posis = [v.co for v in bm.verts if v[cindex] > 0] if svp['liparams']['cp'] == '1' else [f.calc_center_median() for f in bm.faces if f[cindex] > 0]
        rgeom = [g for g in geom if g[cindex] > 0]
        rareas = [gp.calc_area() for gp in geom] if svp['liparams']['cp'] == '0' else [vertarea(bm, gp) for gp in geom]
        reslists.append([str(frame), 'Zone spatial', self.id_data.name, 'X', ' '.join(['{:.3f}'.format(p[0]) for p in posis])])
        reslists.append([str(frame), 'Zone spatial', self.id_data.name, 'Y', ' '.join(['{:.3f}'.format(p[1]) for p in posis])])
        reslists.append([str(frame), 'Zone spatial', self.id_data.name, 'Z', ' '.join(['{:.3f}'.format(p[2]) for p in posis])])
        reslists.append([str(frame), 'Zone spatial', self.id_data.name, 'Areas (m2)', ' '.join(['{:.3f}'.format(ra) for ra in rareas])])

        if svp['liparams']['unit'] == 'W/m2 (f)':
            firradbinvals = [self['omin']['firrad{}'.format(frame)] + (self['omax']['firrad{}'.format(frame)] - self['omin']['firrad{}'.format(frame)])/ll * (i + increment) for i in range(ll)]
            self['livires']['valbins'] = firradbinvals
            reslists.append([str(frame), 'Zone spatial', self.id_data.name, 'Full irradiance (W)', ' '.join(['{:.3f}'.format(g[firradres]) for g in rgeom])])
            reslists.append([str(frame), 'Zone spatial', self.id_data.name, 'Full irradiance (W/m2)', ' '.join(['{:.3f}'.format(g[firradm2res]) for g in rgeom])])

        elif svp['liparams']['unit'] in ('Lux', 'DF (%)'):
            reslists.append([str(frame), 'Zone spatial', self.id_data.name, 'Visible irradiance (W/m2)', ' '.join(['{:.3f}'.format(g[virradm2res]) for g in rgeom])])

            if svp['liparams']['unit'] == 'Lux':
                reslists.append([str(frame), 'Zone spatial', self.id_data.name, 'Illuminance (lux)', ' '.join(['{:.3f}'.format(g[illures]) for g in rgeom])])
            elif svp['liparams']['unit'] == 'DF (%)':
                reslists.append([str(frame), 'Zone spatial', self.id_data.name, 'DF (%)', ' '.join(['{:.3f}'.format(g[dfres]) for g in rgeom])])

        elif svp['liparams']['unit'] == 'W/m2':
            firradbinvals = [self['omin']['firrad{}'.format(frame)] + (self['omax']['firrad{}'.format(frame)] - self['omin']['firrad{}'.format(frame)])/ll * (i + increment) for i in range(ll)]
            self['livires']['valbins'] = firradbinvals
            reslists.append([str(frame), 'Zone spatial', self.id_data.name, 'Full irradiance (W)', ' '.join(['{:.3f}'.format(g[firradres]) for g in rgeom])])
            reslists.append([str(frame), 'Zone spatial', self.id_data.name, 'Full irradiance (W/m2)', ' '.join(['{:.3f}'.format(g[firradm2res]) for g in rgeom])])
            firradrbinvals = [self['omin']['firradr{}'.format(frame)] + (self['omax']['firradr{}'.format(frame)] - self['omin']['firradr{}'.format(frame)])/ll * (i + increment) for i in range(ll)]
            self['livires']['rvalbins'] = firradrbinvals
            reslists.append([str(frame), 'Zone spatial', self.id_data.name, 'Red irradiance (W)', ' '.join(['{:.3f}'.format(g[firradrres]) for g in rgeom])])
            reslists.append([str(frame), 'Zone spatial', self.id_data.name, 'Red irradiance (W/m2)', ' '.join(['{:.3f}'.format(g[firradrm2res]) for g in rgeom])])
            firradbbinvals = [self['omin']['firradb{}'.format(frame)] + (self['omax']['firradb{}'.format(frame)] - self['omin']['firradb{}'.format(frame)])/ll * (i + increment) for i in range(ll)]
            self['livires']['bvalbins'] = firradbbinvals
            reslists.append([str(frame), 'Zone spatial', self.id_data.name, 'Blue irradiance (W)', ' '.join(['{:.3f}'.format(g[firradbres]) for g in rgeom])])
            reslists.append([str(frame), 'Zone spatial', self.id_data.name, 'Blue irradiance (W/m2)', ' '.join(['{:.3f}'.format(g[firradbm2res]) for g in rgeom])])
            firradgbinvals = [self['omin']['firradg{}'.format(frame)] + (self['omax']['firradg{}'.format(frame)] - self['omin']['firradg{}'.format(frame)])/ll * (i + increment) for i in range(ll)]
            self['livires']['bvalbins'] = firradgbinvals
            reslists.append([str(frame), 'Zone spatial', self.id_data.name, 'Green irradiance (W)', ' '.join(['{:.3f}'.format(g[firradgres]) for g in rgeom])])
            reslists.append([str(frame), 'Zone spatial', self.id_data.name, 'Green irradiance (W/m2)', ' '.join(['{:.3f}'.format(g[firradgm2res]) for g in rgeom])])

    if len(frames) > 1:
        reslists.append(['All', 'Frames', 'Frames', 'Frames', ' '.join([str(f) for f in frames])])

        if svp['liparams']['unit'] == 'W/m2 (f)':
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Average full irradiance (W/m2)', ' '.join(['{:.3f}'.format(self['oave']['firradm2{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Maximum full irradiance (W/m2)', ' '.join(['{:.3f}'.format(self['omax']['firradm2{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Minimum full irradiance (W/m2)', ' '.join(['{:.3f}'.format(self['omin']['firradm2{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Average full irradiance (W)', ' '.join(['{:.3f}'.format(self['oave']['firrad{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Maximum full irradiance (W)', ' '.join(['{:.3f}'.format(self['omax']['firrad{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Minimum full irradiance (W)', ' '.join(['{:.3f}'.format(self['omin']['firrad{}'.format(frame)]) for frame in frames])])

        elif svp['liparams']['unit'] in ('Lux', 'DF (%)'):
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Average visible irradiance (W/m2)', ' '.join(['{:.3f}'.format(self['oave']['virradm2{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Maximum visible irradiance (W/m2)', ' '.join(['{:.3f}'.format(self['omax']['virradm2{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Minimum visible irradiance (W/m2)', ' '.join(['{:.3f}'.format(self['omin']['virradm2{}'.format(frame)]) for frame in frames])])

            if svp['liparams']['unit'] == 'Lux':
                reslists.append(['All', 'Zone spatial', self.id_data.name, 'Average illuminance (lux)', ' '.join(['{:.3f}'.format(self['oave']['illu{}'.format(frame)]) for frame in frames])])
                reslists.append(['All', 'Zone spatial', self.id_data.name, 'Maximum illuminance (lux)', ' '.join(['{:.3f}'.format(self['omax']['illu{}'.format(frame)]) for frame in frames])])
                reslists.append(['All', 'Zone spatial', self.id_data.name, 'Minimum illuminance (lux)', ' '.join(['{:.3f}'.format(self['omin']['illu{}'.format(frame)]) for frame in frames])])

            elif svp['liparams']['unit'] == 'DF (%)':
                reslists.append(['All', 'Zone spatial', self.id_data.name, 'Average DF (%)', ' '.join(['{:.3f}'.format(self['oave']['df{}'.format(frame)]) for frame in frames])])
                reslists.append(['All', 'Zone spatial', self.id_data.name, 'Maximum DF (%)', ' '.join(['{:.3f}'.format(self['omax']['df{}'.format(frame)]) for frame in frames])])
                reslists.append(['All', 'Zone spatial', self.id_data.name, 'Minimum DF (%)', ' '.join(['{:.3f}'.format(self['omin']['df{}'.format(frame)]) for frame in frames])])

            ir = []

            if self['oave'].get('illu{}'.format(frames[0])) or self['oave'].get('df{}'.format(frames[0])):
                for frame in frames:
                    dfillu = 'illu' if self['oave'].get('illu{}'.format(frame)) else 'df'
                    if self['oave'][f'{dfillu}{frame}'] > 0:
                        ir.append('{:.3f}'.format(self['omin'][f'{dfillu}{frame}']/self['oave'][f'{dfillu}{frame}']))
                    else:
                        ir.append('0')

                reslists.append(['All', 'Zone spatial', self.id_data.name, 'Uniformity ratio', ' '.join(ir)])

        elif svp['liparams']['unit'] == 'W/m2':
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Average full irradiance (W/m2)', ' '.join(['{:.3f}'.format(self['oave']['firradm2{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Maximum full irradiance (W/m2)', ' '.join(['{:.3f}'.format(self['omax']['firradm2{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Minimum full irradiance (W/m2)', ' '.join(['{:.3f}'.format(self['omin']['firradm2{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Average full irradiance (W)', ' '.join(['{:.3f}'.format(self['oave']['firrad{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Maximum full irradiance (W)', ' '.join(['{:.3f}'.format(self['omax']['firrad{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Minimum full irradiance (W)', ' '.join(['{:.3f}'.format(self['omin']['firrad{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Average red irradiance (W/m2)', ' '.join(['{:.3f}'.format(self['oave']['firradrm2{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Maximum red irradiance (W/m2)', ' '.join(['{:.3f}'.format(self['omax']['firradrm2{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Minimum red irradiance (W/m2)', ' '.join(['{:.3f}'.format(self['omin']['firradrm2{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Average red irradiance (W)', ' '.join(['{:.3f}'.format(self['oave']['firradr{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Maximum red irradiance (W)', ' '.join(['{:.3f}'.format(self['omax']['firradr{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Minimum red irradiance (W)', ' '.join(['{:.3f}'.format(self['omin']['firradr{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Average green irradiance (W/m2)', ' '.join(['{:.3f}'.format(self['oave']['firradgm2{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Maximum green irradiance (W/m2)', ' '.join(['{:.3f}'.format(self['omax']['firradgm2{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Minimum green irradiance (W/m2)', ' '.join(['{:.3f}'.format(self['omin']['firradgm2{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Average green irradiance (W)', ' '.join(['{:.3f}'.format(self['oave']['firradg{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Maximum green irradiance (W)', ' '.join(['{:.3f}'.format(self['omax']['firradg{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Minimum green irradiance (W)', ' '.join(['{:.3f}'.format(self['omin']['firradg{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Average blue irradiance (W/m2)', ' '.join(['{:.3f}'.format(self['oave']['firradbm2{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Maximum blue irradiance (W/m2)', ' '.join(['{:.3f}'.format(self['omax']['firradbm2{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Minimum blue irradiance (W/m2)', ' '.join(['{:.3f}'.format(self['omin']['firradbm2{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Average blue irradiance (W)', ' '.join(['{:.3f}'.format(self['oave']['firradb{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Maximum blue irradiance (W)', ' '.join(['{:.3f}'.format(self['omax']['firradb{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Minimum blue irradiance (W)', ' '.join(['{:.3f}'.format(self['omin']['firradb{}'.format(frame)]) for frame in frames])])

    bm.transform(self.id_data.matrix_world.inverted())
    bm.to_mesh(self.id_data.data)
    bm.free()
    return reslists


def lhcalcapply(self, scene, frames, rtcmds, simnode, curres, pfile):
    reslists = []
    svp = scene.vi_params
    bm = bmesh.new()
    bm.from_mesh(self.id_data.data)
    bm.transform(self.id_data.matrix_world)
    self['omax'], self['omin'], self['oave'] = {}, {}, {}
    clearlayers(bm, 'f')
    geom = bm.faces if svp['liparams']['cp'] == '0' else bm.verts
    cindex = geom.layers.int['cindex']

    for f, frame in enumerate(frames):
        geom.layers.float.new('res{}'.format(frame))

        if simnode['coptions']['unit'] == 'klxh':
            geom.layers.float.new('virradhm2{}'.format(frame))
            geom.layers.float.new('virradh{}'.format(frame))
            geom.layers.float.new('illuh{}'.format(frame))
            virradm2res = geom.layers.float['virradhm2{}'.format(frame)]
            virradres = geom.layers.float['virradh{}'.format(frame)]
            illures = geom.layers.float['illuh{}'.format(frame)]
        elif simnode['coptions']['unit'] == 'kWh (f)':
            geom.layers.float.new('firradhm2{}'.format(frame))
            geom.layers.float.new('firradh{}'.format(frame))
            firradm2res = geom.layers.float['firradhm2{}'.format(frame)]
            firradres = geom.layers.float['firradh{}'.format(frame)]

        if geom.layers.string.get('rt{}'.format(frame)):
            rtframe = frame
        else:
            kints = [int(k[2:]) for k in geom.layers.string.keys()]
            rtframe = max(kints) if frame > max(kints) else min(kints)

        rt = geom.layers.string['rt{}'.format(rtframe)]
        gps = [g for g in geom if g[rt]]
        areas = array([g.calc_area() for g in gps] if svp['liparams']['cp'] == '0' else [vertarea(bm, g) for g in gps])

        for chunk in chunks(gps, int(svp['viparams']['nproc']) * 20):
            careas = array([c.calc_area() if svp['liparams']['cp'] == '0' else vertarea(bm, c) for c in chunk])
            rtrun = Popen(shlex.split(rtcmds[f]), stdin=PIPE, stdout=PIPE, stderr=PIPE, universal_newlines=True).communicate(input='\n'.join([c[rt].decode('utf-8') for c in chunk]))
            logentry('Running rtrace with command: {}'.format(rtcmds[f]))

            if rtrun[1]:
                logentry('Rtrace error: {}'.format(rtrun[1]))
                return 'CANCELLED'

            xyzirrad = array([[float(v) for v in sl.split('\t')[:3]] for sl in rtrun[0].splitlines()])

            if simnode['coptions']['unit'] == 'klxh':
                virradm2 = nsum(xyzirrad * array([0.333, 0.333, 0.333]), axis=1) * 1e-3
                illu = nsum(xyzirrad * array([0.265, 0.67, 0.065]), axis=1) * 0.179
                virrad = virradm2 * careas

            elif simnode['coptions']['unit'] == 'kWh (f)':
                firradm2 = nsum(xyzirrad * array([0.333, 0.333, 0.333]), axis=1) * 1e-3
                firrad = firradm2 * careas

            for gi, gp in enumerate(chunk):
                if simnode['coptions']['unit'] == 'klxh':
                    gp[virradres] = virrad[gi].astype(float32)
                    gp[virradm2res] = virradm2[gi].astype(float32)
                    gp[illures] = illu[gi].astype(float32)
                elif simnode['coptions']['unit'] == 'kWh (f)':
                    gp[firradm2res] = firradm2[gi].astype(float32)
                    gp[firradres] = firrad[gi].astype(float32)

            curres += len(chunk)

            if pfile.check(curres) == 'CANCELLED':
                bm.free()
                return {'CANCELLED'}

        if simnode['coptions']['unit'] == 'klxh':
            oillu = array([g[illures] for g in gps])
            ovirrad = array([g[virradres] for g in gps])
            ovirradm2 = array([g[virradm2res] for g in gps])
            maxoillu = nmax(oillu)
            maxovirrad = nmax(ovirrad)
            maxovirradm2 = nmax(ovirradm2)
            minovirradm2 = nmin(ovirradm2)
            minovirrad = nmin(ovirrad)
            minoillu = nmin(oillu)
            aveovirradm2 = nmean(ovirradm2)
            aveovirrad = nmean(ovirrad)
            aveoillu = nmean(oillu)
            self['omax']['virradh{}'.format(frame)] = maxovirrad
            self['omin']['virradh{}'.format(frame)] = minovirrad
            self['oave']['virradh{}'.format(frame)] = aveovirrad
            self['omax']['virradhm2{}'.format(frame)] = maxovirradm2
            self['omin']['virradhm2{}'.format(frame)] = minovirradm2
            self['oave']['virradhm2{}'.format(frame)] = aveovirradm2
            self['omax']['illuh{}'.format(frame)] = maxoillu
            self['omin']['illuh{}'.format(frame)] = minoillu
            self['oave']['illuh{}'.format(frame)] = aveoillu

        elif simnode['coptions']['unit'] == 'kWh (f)':
            ofirradm2 = array([g[firradm2res] for g in gps])
            ofirrad = array([g[firradres] for g in gps])
            maxofirradm2 = nmax(ofirradm2)
            maxofirrad = nmax(ofirrad)
            minofirradm2 = nmin(ofirradm2)
            minofirrad = nmin(ofirrad)
            aveofirradm2 = nmean(ofirradm2)
            aveofirrad = nmean(ofirrad)
            self['omax']['firradh{}'.format(frame)] = maxofirrad
            self['omin']['firradh{}'.format(frame)] = minofirrad
            self['oave']['firradh{}'.format(frame)] = aveofirrad
            self['omax']['firradhm2{}'.format(frame)] = maxofirradm2
            self['omin']['firradhm2{}'.format(frame)] = minofirradm2
            self['oave']['firradhm2{}'.format(frame)] = aveofirradm2

        posis = [v.co for v in bm.verts if v[cindex] > 0] if svp['liparams']['cp'] == '1' else [f.calc_center_median() for f in bm.faces if f[cindex] > 1]
        reslists.append([str(frame), 'Zone spatial', self.id_data.name, 'X', ' '.join([str(p[0]) for p in posis])])
        reslists.append([str(frame), 'Zone spatial', self.id_data.name, 'Y', ' '.join([str(p[1]) for p in posis])])
        reslists.append([str(frame), 'Zone spatial', self.id_data.name, 'Z', ' '.join([str(p[2]) for p in posis])])
        reslists.append([str(frame), 'Zone spatial', self.id_data.name, 'Area', ' '.join([str(a) for a in areas])])

        if simnode['coptions']['unit'] == 'klxh':
            reslists.append([str(frame), 'Zone spatial', self.id_data.name, 'Visible irradiance (kwh)', ' '.join([str(g[virradres]) for g in geom if g[cindex] > 0])])
            reslists.append([str(frame), 'Zone spatial', self.id_data.name, 'Illuminance (klxh)', ' '.join([str(g[illures]) for g in geom if g[cindex] > 0])])
            reslists.append([str(frame), 'Zone spatial', self.id_data.name, 'Visible irradiance (kwh/m2)', ' '.join([str(g[virradm2res]) for g in geom if g[cindex] > 0])])
        elif simnode['coptions']['unit'] == 'kWh (f)':
            reslists.append([str(frame), 'Zone spatial', self.id_data.name, 'Full irradiance (kwh)', ' '.join([str(g[firradres]) for g in geom if g[cindex] > 0])])
            reslists.append([str(frame), 'Zone spatial', self.id_data.name, 'Full irradiance (kwh/m2)', ' '.join([str(g[firradm2res]) for g in geom if g[cindex] > 0])])

    if len(frames) > 1:
        reslists.append(['All', 'Frames', 'Frames', 'Frames', ' '.join([str(f) for f in frames])])

        if simnode['coptions']['unit'] == 'klxh':
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Average klxh', ' '.join(['{:.3f}'.format(self['oave']['illuh{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Maximum klxh', ' '.join(['{:.3f}'.format(self['omax']['illuh{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Minimum klxh', ' '.join(['{:.3f}'.format(self['omin']['illuh{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Average visual radiation exposure (kWh)', ' '.join(['{:.3f}'.format(self['oave']['virradh{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Maximum visual radiation exposure (kWh)', ' '.join(['{:.3f}'.format(self['omax']['virradh{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Minimum visual radiation exposure (kWh)', ' '.join(['{:.3f}'.format(self['omin']['virradh{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Average visual radiation exposure (kWh/m2)', ' '.join(['{:.3f}'.format(self['oave']['virradhm2{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Maximum visual radiation exposure (kWh/m2)', ' '.join(['{:.3f}'.format(self['omax']['virradhm2{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Minimum visual radiation exposure (kWh/m2)', ' '.join(['{:.3f}'.format(self['omin']['virradhm2{}'.format(frame)]) for frame in frames])])

        elif simnode['coptions']['unit'] == 'kWh (f)':
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Average solar exposure (kWh)', ' '.join(['{:.3f}'.format(self['oave']['firradh{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Maximum solar exposure (kWh)', ' '.join(['{:.3f}'.format(self['omax']['firradh{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Minimum solar exposure (kWh)', ' '.join(['{:.3f}'.format(self['omin']['firradh{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Average solar exposure (kWh/m2)', ' '.join(['{:.3f}'.format(self['oave']['firradhm2{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Maximum solar exposure (kWh/m2)', ' '.join(['{:.3f}'.format(self['omax']['firradhm2{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Minimum solar exposure (kWh/m2)', ' '.join(['{:.3f}'.format(self['omin']['firradhm2{}'.format(frame)]) for frame in frames])])

    bm.transform(self.id_data.matrix_world.inverted())
    bm.to_mesh(self.id_data.data)
    bm.free()
    return reslists


def udidacalcapply(self, scene, frames, rccmds, simnode, curres, pfile):
    svp = scene.vi_params
    self['livires'] = {}
    reslists = []
    selobj(bpy.context.view_layer, self.id_data)
    bm = bmesh.new()
    bm.from_mesh(self.id_data.data)
    bm.transform(self.id_data.matrix_world)
    bm.normal_update()
    clearlayers(bm, 'f')
    geom = bm.faces if svp['liparams']['cp'] == '0' else bm.verts
    self['omax'], self['omin'], self['oave'] = {}, {}, {}
    mtxlines = open(simnode.inputs['Context in'].links[0].from_node['Options']['mtxfile'], 'r').readlines()
    mtxlinesns = open(simnode.inputs['Context in'].links[0].from_node['Options']['mtxfilens'], 'r').readlines()

    for line in mtxlines:
        if line.split("=")[0] == 'NROWS':
            patches = int(line.split("=")[1])
            break

    if self.get('wattres'):
        del self['wattres']

    illumod = array((47.4, 120, 11.6)).astype(float32)
    wattmod = array((0.265, 0.67, 0.065)).astype(float32)
    times = [datetime.datetime.strptime(time, "%d/%m/%y %H:%M:%S") for time in simnode['coptions']['times']]
    vecvals, vals = mtx2vals(mtxlines, datetime.datetime(2015, 1, 1).weekday(), simnode, times)
    vecvalsns, valsns = mtx2vals(mtxlinesns, datetime.datetime(2015, 1, 1).weekday(), simnode, times)
    cbdm_days = list(set([t.timetuple().tm_yday for t in times]))
    cbdm_hours = [h for h in range(simnode['coptions']['cbdm_sh'], simnode['coptions']['cbdm_eh'] + 1)]
    dno, hno = len(cbdm_days), len(cbdm_hours)
    sdah = 10 if simnode['coptions']['ay'] else hno
    (luxmin, luxmax) = (simnode['coptions']['dalux'], simnode['coptions']['asemax'])
    vecvals = array([vv[2:] for vv in vecvals if vv[1] < simnode['coptions']['weekdays']]).astype(float32)
    vecvalsns = array([vv[2:] for vv in vecvalsns if vv[1] < simnode['coptions']['weekdays']]).astype(float32)
    hours = vecvals.shape[0]
    logentry(f'Running CBDM calculation over {hours} hours')
    hour_array = array([t.hour for t in times])
    restypes = ('da', 'sda', 'sv', 'ase', 'res', 'udilow', 'udisup', 'udiauto', 'udihi', 'firradh', 'firradhm2', 'maxlux', 'minlux', 'avelux')
    self['livires']['cbdm_days'] = cbdm_days
    self['livires']['cbdm_hours'] = cbdm_hours

    for f, frame in enumerate(frames):
        reslists.append([str(frame), 'Time', 'Time', 'Month', ' '.join([str(t.month) for t in times])])
        reslists.append([str(frame), 'Time', 'Time', 'Day', ' '.join([str(t.day) for t in times])])
        reslists.append([str(frame), 'Time', 'Time', 'Hour', ' '.join([str(t.hour) for t in times])])
        reslists.append([str(frame), 'Time', 'Time', 'DOS', ' '.join([str(t.timetuple().tm_yday - times[0].timetuple().tm_yday) for t in times])])

        for restype in restypes:
            geom.layers.float.new('{}{}'.format(restype, frame))

        (resda, ressda, ressv, resase, res, resudilow, resudisup, resudiauto, resudihi, firrad, firradm2, maxillu, minillu, aveillu) = [geom.layers.float['{}{}'.format(r, frame)] for r in restypes]

        if geom.layers.string.get('rt{}'.format(frame)):
            rtframe = frame
        else:
            kints = [int(k[2:]) for k in geom.layers.string.keys()]
            rtframe = max(kints) if frame > max(kints) else min(kints)

        rt = geom.layers.string['rt{}'.format(rtframe)]
        rgeom = [f for f in bm.faces if f[rt]] if svp['liparams']['cp'] == '0' else [v for v in bm.verts if v[rt]]
        reslen = len(rgeom)
        areas = array([g.calc_area() for g in rgeom] if svp['liparams']['cp'] == '0' else [vertarea(bm, g) for g in rgeom])
        totarea = sum(areas)

        for ch, chunk in enumerate(chunks([g for g in rgeom], int(svp['viparams']['nproc']) * 20)):
            if not ch:
                logentry(f"Running rcontrib with the command: {rccmds[f]}")

            sensrun = Popen(shlex.split(rccmds[f]), stdin=PIPE, stdout=PIPE, stderr=PIPE, universal_newlines=True).communicate(input='\n'.join([c[rt].decode('utf-8') for c in chunk]))
            resarray = array([[float(v) for v in sl.strip('\n').strip('\r\n').split('\t') if v] for sl in sensrun[0].splitlines()]).reshape(len(chunk), patches, 3).astype(float32)
            chareas = array([c.calc_area() for c in chunk]) if svp['liparams']['cp'] == '0' else array([vertarea(bm, c) for c in chunk])
            sensarray = nsum(resarray*illumod, axis=2).astype(float32)
            finalillu = inner(sensarray, vecvals).astype(float64)

            if svp['viparams']['visimcontext'] == 'LiVi CBDM' and simnode['coptions']['cbanalysis'] == '1':
                wattarray = nsum(resarray*wattmod, axis=2).astype(float32)
                firradm2array = inner(wattarray, vecvals).astype(float32)
                firradarray = (firradm2array.T * chareas).T.astype(float32)
                kwh = 0.001 * nsum(firradarray, axis=1)
                kwhm2 = 0.001 * nsum(firradm2array, axis=1)

                if not ch:
                    totfinalwatt = nsum(firradarray, axis=0)
                    totfinalwattm2 = average(firradm2array, axis=0)
                    finalkwh = kwh
                    finalkwhm2 = kwhm2
                else:
                    totfinalwatt += nsum(firradarray, axis=0)
                    totfinalwattm2 += average(firradm2array, axis=0)
                    finalkwh = concatenate((finalkwh, kwh))
                    finalkwhm2 = concatenate((finalkwhm2, kwhm2))

                for gi, gp in enumerate(chunk):
                    gp[firrad] = kwh[gi]
                    gp[firradm2] = kwhm2[gi]

            elif svp['viparams']['visimcontext'] == 'LiVi CBDM' and simnode['coptions']['cbanalysis'] == '2':
                rclist = shlex.split(rccmds[f], posix=False)
                rccmd = ' '.join(rclist[:4] + ['-ab 1 -ad 8192 -lw 0.0001 -lr 0'] + rclist[-13:])

                if not ch:
                    logentry(f"Running rcontrib (no sky) with the command: {rccmd}")

                sensrunns = Popen(shlex.split(rccmd), stdin=PIPE, stdout=PIPE, stderr=PIPE, universal_newlines=True).communicate(input='\n'.join([c[rt].decode('utf-8') for c in chunk]))
                resarrayd = array([[float(v) for v in sl.strip('\n').strip('\r\n').split('\t') if v] for sl in sensrunns[0].splitlines()]).reshape(len(chunk), patches, 3).astype(float32)
                sensarrayns = nsum(resarrayd*illumod, axis=2).astype(float32)
                finalilluns = inner(sensarrayns, vecvalsns).astype(float32)
                dabool = choose(finalillu >= luxmin, [0, 1]).astype(int8)
                asebool = choose(finalilluns >= luxmax, [0, 1]).astype(int8)
                udilbool = choose(finalillu < simnode['coptions']['damin'], [0, 1]).astype(int8)
                udisbool = choose(finalillu < simnode['coptions']['dasupp'], [0, 1]).astype(int8) - udilbool
                udiabool = choose(finalillu < simnode['coptions']['daauto'], [0, 1]).astype(int8) - udilbool - udisbool
                udihbool = choose(finalillu >= simnode['coptions']['daauto'], [0, 1]).astype(int8)
                svbool = choose(nsum(nsum(resarrayd, axis=1), axis=1) > 0.0, [0, 1]).astype(int8)
                sdafinalillu = finalillu[:, logical_and(8 <= hour_array, hour_array < 18)] if simnode['coptions']['ay'] else finalillu
                sdabool = choose(sdafinalillu >= 300, [0, 1]).astype(int8)
                daareares = (dabool.T*chareas).T
                udilareares = (udilbool.T*chareas).T
                udisareares = (udisbool.T*chareas).T
                udiaareares = (udiabool.T*chareas).T
                udihareares = (udihbool.T*chareas).T
                aseareares = (asebool.T*chareas).T
                sdaareares = (sdabool.T*chareas).T
                sdaarearespa = (sdabool.T*chareas*svbool).T
                dares = dabool.sum(axis=1)*100/hours
                udilow = udilbool.sum(axis=1)*100/hours
                udisup = udisbool.sum(axis=1)*100/hours
                udiauto = udiabool.sum(axis=1)*100/hours
                udihi = udihbool.sum(axis=1)*100/hours
                sdares = sdabool.sum(axis=1)*100/(dno * sdah)
                aseres = asebool.sum(axis=1)*1.0

                if not ch:
                    totfinalillu = finalillu
                    totdaarea = nsum(100 * daareares/totarea, axis=0)
                    totudiaarea = nsum(100 * udiaareares/totarea, axis=0)
                    totudisarea = nsum(100 * udisareares/totarea, axis=0)
                    totudilarea = nsum(100 * udilareares/totarea, axis=0)
                    totudiharea = nsum(100 * udihareares/totarea, axis=0)
                    totsdaarea = nsum(sdaareares, axis=0)
                    totsdaareapa = nsum(sdaarearespa, axis=0)
                    totasearea = nsum(aseareares, axis=0)
                    svarea = nsum(chareas * svbool)
                else:
                    nappend(totfinalillu, finalillu)
                    totdaarea += nsum(100 * daareares/totarea, axis=0)
                    totudiaarea += nsum(100 * udiaareares/totarea, axis=0)
                    totudilarea += nsum(100 * udilareares/totarea, axis=0)
                    totudisarea += nsum(100 * udisareares/totarea, axis=0)
                    totudiharea += nsum(100 * udihareares/totarea, axis=0)
                    totsdaarea += nsum(sdaareares, axis=0)
                    totsdaareapa += nsum(sdaarearespa, axis=0)
                    totasearea += nsum(aseareares, axis=0)
                    svarea += nsum(chareas * svbool)

                for gi, gp in enumerate(chunk):
                    gp[resda] = dares[gi]
                    gp[res] = dares[gi]
                    gp[resudilow] = udilow[gi]
                    gp[resudisup] = udisup[gi]
                    gp[resudiauto] = udiauto[gi]
                    gp[resudihi] = udihi[gi]
                    gp[maxillu] = max(finalillu[gi])
                    gp[minillu] = min(finalillu[gi])
                    gp[aveillu] = nmean(finalillu[gi])
                    gp[ressda] = sdares[gi]
                    gp[resase] = aseres[gi]
                    gp[ressv] = svbool[gi]

            curres += len(chunk)

            if pfile.check(curres) == 'CANCELLED':
                bm.free()
                return {'CANCELLED'}

        if svp['viparams']['visimcontext'] == 'LiVi CBDM' and simnode['coptions']['cbanalysis'] == '1':
            self['omax']['firradh{}'.format(frame)] = nmax(finalkwh).astype(float64)
            self['omin']['firradh{}'.format(frame)] = nmin(finalkwh).astype(float64)
            self['oave']['firradh{}'.format(frame)] = nmean(finalkwh).astype(float64)
            self['omax']['firradhm2{}'.format(frame)] = nmax(finalkwhm2).astype(float64)
            self['omin']['firradhm2{}'.format(frame)] = nmin(finalkwhm2).astype(float64)
            self['oave']['firradhm2{}'.format(frame)] = nmean(finalkwhm2).astype(float64)
            self['livires']['firradh{}'.format(frame)] = (0.001*totfinalwatt).reshape(dno, hno).transpose().tolist()
            self['livires']['firradhm2{}'.format(frame)] = (0.001*totfinalwattm2).reshape(dno, hno).transpose().tolist()
            reslists.append([str(frame), 'Zone temporal', self.id_data.name, 'sum kW', ' '.join([str(p) for p in 0.001 * totfinalwatt])])
            reslists.append([str(frame), 'Zone temporal', self.id_data.name, 'ave kW/m2', ' '.join([str(p) for p in 0.001 * totfinalwattm2])])

        elif svp['viparams']['visimcontext'] == 'LiVi CBDM' and simnode['coptions']['cbanalysis'] == '2':
            dares = [gp[resda] for gp in rgeom]
            udilow = [gp[resudilow] for gp in rgeom]
            udisup = [gp[resudisup] for gp in rgeom]
            udiauto = [gp[resudiauto] for gp in rgeom]
            udihi = [gp[resudihi] for gp in rgeom]
            sdas = [gp[ressda] for gp in rgeom]
            ases = [gp[resase] for gp in rgeom]
            svres = [gp[ressv] for gp in rgeom]
            self['omax']['udilow{}'.format(frame)] = max(udilow)
            self['omin']['udilow{}'.format(frame)] = min(udilow)
            self['oave']['udilow{}'.format(frame)] = nmean(udilow)
            self['omax']['udisup{}'.format(frame)] = max(udisup)
            self['omin']['udisup{}'.format(frame)] = min(udisup)
            self['oave']['udisup{}'.format(frame)] = nmean(udisup)
            self['omax']['udiauto{}'.format(frame)] = max(udiauto)
            self['omin']['udiauto{}'.format(frame)] = min(udiauto)
            self['oave']['udiauto{}'.format(frame)] = sum(udiauto)
            self['omax']['udihi{}'.format(frame)] = max(udihi)
            self['omin']['udihi{}'.format(frame)] = min(udihi)
            self['oave']['udihi{}'.format(frame)] = nmean(udihi)
            self['omax']['da{}'.format(frame)] = max(dares)
            self['omin']['da{}'.format(frame)] = min(dares)
            self['oave']['da{}'.format(frame)] = nmean(dares)
            self['omax']['maxlux{}'.format(frame)] = max(nmax(totfinalillu, axis=1).astype(float64))
            self['omin']['maxlux{}'.format(frame)] = min(nmax(totfinalillu, axis=1).astype(float64))
            self['oave']['maxlux{}'.format(frame)] = nmean(nmax(totfinalillu, axis=1).astype(float64))
            self['omax']['minlux{}'.format(frame)] = max(nmin(totfinalillu, axis=1).astype(float64))
            self['omin']['minlux{}'.format(frame)] = min(nmin(totfinalillu, axis=1).astype(float64))
            self['oave']['minlux{}'.format(frame)] = nmean(nmin(totfinalillu, axis=1).astype(float64))
            self['omax']['avelux{}'.format(frame)] = max(nmean(totfinalillu, axis=1).astype(float64))
            self['omin']['avelux{}'.format(frame)] = min(nmean(totfinalillu, axis=1).astype(float64))
            self['oave']['avelux{}'.format(frame)] = nmean(nmean(totfinalillu, axis=1).astype(float64))
            self['livires']['dhilluave{}'.format(frame)] = average(totfinalillu, axis=0).flatten().reshape(dno, hno).transpose().tolist()
            self['livires']['dhillumin{}'.format(frame)] = amin(totfinalillu, axis=0).reshape(dno, hno).transpose().tolist()
            self['livires']['dhillumax{}'.format(frame)] = amax(totfinalillu, axis=0).reshape(dno, hno).transpose().tolist()
            self['livires']['daarea{}'.format(frame)] = totdaarea.reshape(dno, hno).transpose().tolist()
            self['livires']['udiaarea{}'.format(frame)] = totudiaarea.reshape(dno, hno).transpose().tolist()
            self['livires']['udisarea{}'.format(frame)] = totudisarea.reshape(dno, hno).transpose().tolist()
            self['livires']['udilarea{}'.format(frame)] = totudilarea.reshape(dno, hno).transpose().tolist()
            self['livires']['udiharea{}'.format(frame)] = totudiharea.reshape(dno, hno).transpose().tolist()
            self['omax']['sv{}'.format(frame)] = max(svres)
            self['omin']['sv{}'.format(frame)] = min(svres)
            self['oave']['sv{}'.format(frame)] = nmean(svres)
            reslists.append([str(frame), 'Zone temporal', self.id_data.name, 'Daylight Autonomy Area (%)', ' '.join([f'{p:.2f}' for p in totdaarea])])
            reslists.append([str(frame), 'Zone temporal', self.id_data.name, 'UDI-a Area (%)', ' '.join([f'{p:.2f}' for p in totudiaarea])])
            reslists.append([str(frame), 'Zone temporal', self.id_data.name, 'UDI-s Area (%)', ' '.join([f'{p:.2f}' for p in totudisarea])])
            reslists.append([str(frame), 'Zone temporal', self.id_data.name, 'UDI-l Area (%)', ' '.join([f'{p:.2f}' for p in totudilarea])])
            reslists.append([str(frame), 'Zone temporal', self.id_data.name, 'UDI-h Area (%)', ' '.join([f'{p:.2f}' for p in totudiharea])])
            overallsdaareapa = sum([g.calc_area() for g in rgeom if g[ressv] == 1.0]) if svp['liparams']['cp'] == '0' else sum([vertarea(bm, g) for g in rgeom if g[ressv] == 1.0])
            overallsdaarea = totarea
            self['omax']['sda{}'.format(frame)] = max(sdas)
            self['omin']['sda{}'.format(frame)] = min(sdas)
            self['oave']['sda{}'.format(frame)] = sum(sdas)/reslen
            self['omax']['ase{}'.format(frame)] = max(ases)
            self['omin']['ase{}'.format(frame)] = min(ases)
            self['oave']['ase{}'.format(frame)] = sum(ases)/reslen
            self['livires']['asearea{}'.format(frame)] = (100 * totasearea/totarea).reshape(dno, hno).transpose().tolist()
            self['livires']['sdaarea{}'.format(frame)] = (100 * totsdaarea/overallsdaarea).reshape(dno, sdah).transpose().tolist()
            self['livires']['sdaareapa{}'.format(frame)] = (100 * totsdaareapa/svarea).reshape(dno, sdah).transpose().tolist()
            self['livires']['totarea{}'.format(frame)] = totarea
            self['livires']['svarea{}'.format(frame)] = svarea
            self['livires']['ase{}'.format(frame)] = sum(areas[array(ases) > 250])/totarea
            sumsdaareas = sum(areas[array(sdas) >= 50])
            self['livires']['sda{}'.format(frame)] = sumsdaareas/totarea
            self['livires']['sdapa{}'.format(frame)] = sumsdaareas/svarea
            reslists.append([str(frame), 'Zone temporal', self.id_data.name, 'Spatial Daylight Autonomy (% area)', ' '.join([f'{p:.2f}' for p in 100 * totsdaarea/totarea])])

            if svarea:
                reslists.append([str(frame), 'Zone temporal', self.id_data.name, 'Spatial Daylight Autonomy (% perimeter area)', ' '.join([f'{p:.2f}' for p in 100 * totsdaareapa/svarea])])
            else:
                logentry(f'Sensor object {self.id_data.name} cannot see any sky.')

            reslists.append([str(frame), 'Zone temporal', self.id_data.name, 'Annual Sunlight Exposure (% area)', ' '.join([f'{p:.2f}' for p in 100 * totasearea/totarea])])

    if len(frames) > 1:
        reslists.append(['All', 'Frames', 'Frames', 'Frames', ' '.join([str(f) for f in frames])])
        if svp['viparams']['visimcontext'] == 'LiVi CBDM' and simnode['coptions']['cbanalysis'] == '1':
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Average solar exposure (kWh)', ' '.join(['{:.3f}'.format(self['oave']['firradh{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Maximum solar exposure (kWh)', ' '.join(['{:.3f}'.format(self['omax']['firradh{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Minimum solar exposure (kWh)', ' '.join(['{:.3f}'.format(self['omin']['firradh{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Average solar exposure (kWh/m2)', ' '.join(['{:.3f}'.format(self['oave']['firradhm2{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Maximum solar exposure (kWh/m2)', ' '.join(['{:.3f}'.format(self['omax']['firradhm2{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Minimum solar exposure (kWh/m2)', ' '.join(['{:.3f}'.format(self['omin']['firradhm2{}'.format(frame)]) for frame in frames])])

        elif svp['viparams']['visimcontext'] == 'LiVi CBDM' and simnode['coptions']['cbanalysis'] == '2':
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Average sDA (% area)', ' '.join(['{:.3f}'.format(self['oave']['sda{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Maximum sDA (% area)', ' '.join(['{:.3f}'.format(self['omax']['sda{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Minimum sDA (% area)', ' '.join(['{:.3f}'.format(self['omin']['sda{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Average DA (% area)', ' '.join(['{:.3f}'.format(self['oave']['da{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Maximum DA (% area)', ' '.join(['{:.3f}'.format(self['omax']['da{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Minimum DA (% area)', ' '.join(['{:.3f}'.format(self['omin']['da{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Average ASE (hours)', ' '.join(['{:.3f}'.format(self['oave']['ase{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Maximum ASE (hours)', ' '.join(['{:.3f}'.format(self['omax']['ase{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Minimum ASE (hours)', ' '.join(['{:.3f}'.format(self['omin']['ase{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Average UDI-l (% area)', ' '.join(['{:.3f}'.format(self['oave']['udilow{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Maximum UDI-l (% area)', ' '.join(['{:.3f}'.format(self['omax']['udilow{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Minimum UDI-l (% area)', ' '.join(['{:.3f}'.format(self['omin']['udilow{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Average UDI-s (% area)', ' '.join(['{:.3f}'.format(self['oave']['udisup{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Maximum UDI-s (% area)', ' '.join(['{:.3f}'.format(self['omax']['udisup{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Minimum UDI-s (% area)', ' '.join(['{:.3f}'.format(self['omin']['udisup{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Average UDI-a (% area)', ' '.join(['{:.3f}'.format(self['oave']['udiauto{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Maximum UDI-a (% area)', ' '.join(['{:.3f}'.format(self['omax']['udiauto{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Minimum UDI-a (% area)', ' '.join(['{:.3f}'.format(self['omin']['udiauto{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Average UDI-h (% area)', ' '.join(['{:.3f}'.format(self['oave']['udihi{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Maximum UDI-h (% area)', ' '.join(['{:.3f}'.format(self['omax']['udihi{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone spatial', self.id_data.name, 'Minimum UDI-h (% area)', ' '.join(['{:.3f}'.format(self['omin']['udihi{}'.format(frame)]) for frame in frames])])

    bm.transform(self.id_data.matrix_world.inverted())
    bm.to_mesh(self.id_data.data)
    bm.free()
    return reslists

def adgpcalcapply(self, scene, frames, rccmds, simnode, curres, pfile):
    svp = scene.vi_params
    self['livires'] = {}
    self['azi'] = simnode['coptions']['dgp_azi']
    reslists = []
    selobj(bpy.context.view_layer, self.id_data)
    bm = bmesh.new()
    bm.from_mesh(self.id_data.data)
    bm.transform(self.id_data.matrix_world)
    bm.normal_update()
    clearlayers(bm, 'f')
    geom = bm.faces if svp['liparams']['cp'] == '0' else bm.verts
    self['omax'], self['omin'], self['oave'] = {}, {}, {}
    mtxlines = open(simnode.inputs['Context in'].links[0].from_node['Options']['mtxfile'], 'r').readlines()
    mtxlinesns = open(simnode.inputs['Context in'].links[0].from_node['Options']['mtxfilens'], 'r').readlines()
    cbdm_res = (1, 2, 4, 8)[(146, 578, 2306, 9218).index(simnode['coptions']['cbdm_res'])]

    for line in mtxlines:
        if line.split("=")[0] == 'NROWS':
            patches = int(line.split("=")[1])
            break

    if self.get('wattres'):
        del self['wattres']

    times = [datetime.datetime.strptime(time, "%d/%m/%y %H:%M:%S") for time in simnode['coptions']['times']]
    vecvals, vals = mtx2vals(mtxlines, datetime.datetime(2015, 1, 1).weekday(), simnode, times)
    vecvalsns, valsns = mtx2vals(mtxlinesns, datetime.datetime(2015, 1, 1).weekday(), simnode, times)
    cbdm_days = list(set([t.timetuple().tm_yday for t in times]))
    cbdm_hours = [h for h in range(simnode['coptions']['cbdm_sh'], simnode['coptions']['cbdm_eh'] + 1)]
    dno, hno = len(cbdm_days), len(cbdm_hours)
    vecvals = array([vv[2:] for vv in vecvals if vv[1] < simnode['coptions']['weekdays']]).astype(float32)
    vecvalsns = array([vv[2:] for vv in vecvalsns if vv[1] < simnode['coptions']['weekdays']]).astype(float32)
    hours = vecvals.shape[0]
    logentry(f'Running annual glare calculation over {hours} hours')
    hour_array = array([t.hour for t in times])
    self['livires']['cbdm_days'] = cbdm_days
    self['livires']['cbdm_hours'] = cbdm_hours

    for f, frame in enumerate(frames):
        if simnode['coptions']['dgp_hourly']:
            reslists.append([str(frame), 'Time', 'Time', 'Month', ' '.join([str(t.month) for t in times])])
            reslists.append([str(frame), 'Time', 'Time', 'Day', ' '.join([str(t.day) for t in times])])
            reslists.append([str(frame), 'Time', 'Time', 'Hour', ' '.join([str(t.hour) for t in times])])
            reslists.append([str(frame), 'Time', 'Time', 'DOS', ' '.join([str(t.timetuple().tm_yday - times[0].timetuple().tm_yday) for t in times])])

        agas, agos = [], []
        
        for agv in range(simnode['coptions']['ga_views']):
            geom.layers.float.new('{}{}'.format(f'aga{agv + 1}v', frame))
            geom.layers.float.new('{}{}'.format(f'ago{agv + 1}v', frame))
            agas.append(geom.layers.float['{}{}'.format(f'aga{agv + 1}v', frame)])
            agos.append(geom.layers.float['{}{}'.format(f'ago{agv + 1}v', frame)])

            if not agv:
                logentry('Calculating view {} with azimuth {}'.format(agv + 1, simnode['coptions']['dgp_azi'])) 
            else:
                logentry('Calculating view {} with azimuth {}'.format(agv + 1, simnode['coptions']['dgp_azi'] + agv * 360/simnode['coptions']['ga_views'])) 
        
        if geom.layers.string.get('rt{}'.format(frame)):
            rtframe = frame
        else:
            kints = [int(k[2:]) for k in geom.layers.string.keys()]
            rtframe = max(kints) if frame > max(kints) else min(kints)

        rt = geom.layers.string['rt{}'.format(rtframe)]
        
        if not geom.layers.string.get('rv{}'.format(frame)):
            geom.layers.string.new('rv{}'.format(rtframe))
        
        rv = geom.layers.string['rv{}'.format(rtframe)]

        rgeom = [f for f in bm.faces if f[rt]] if svp['liparams']['cp'] == '0' else [v for v in bm.verts if v[rt]]
        reslen = len(rgeom)
        views_text = ''

        for r in rgeom:
            for view in range(simnode['coptions']['ga_views']):
                if not view:
                    r[rv] = '{} {:.3f} {:.3f} 0'.format(' '.join(r[rt].decode('utf-8').split()[:3]), sin(simnode['coptions']['dgp_azi']*pi/180), cos(simnode['coptions']['dgp_azi'] * pi/180)).encode('utf-8')
                    views_text += '{} {:.3f} {:.3f} 0\n'.format(' '.join(r[rt].decode('utf-8').split()[:3]), 
                                                                        sin(simnode['coptions']['dgp_azi']*pi/180), 
                                                                        cos(simnode['coptions']['dgp_azi'] * pi/180))    
                else:
                    r[rv] = '{} {:.3f} {:.3f} 0'.format(' '.join(r[rt].decode('utf-8').split()[:3]), 
                                                                sin((simnode['coptions']['dgp_azi'] + view * 360/simnode['coptions']['ga_views'])*pi/180),
                                                                cos((simnode['coptions']['dgp_azi'] + view * 360/simnode['coptions']['ga_views'])*pi/180)).encode('utf-8')
                    views_text += '{} {:.3f} {:.3f} 0\n'.format(' '.join(r[rt].decode('utf-8').split()[:3]), 
                                                                        sin((simnode['coptions']['dgp_azi'] + view * 360/simnode['coptions']['ga_views'])*pi/180),
                                                                        cos((simnode['coptions']['dgp_azi'] + view * 360/simnode['coptions']['ga_views'])*pi/180))
                    
        with open(os.path.join(svp['viparams']['newdir'], 'views_{}.vf'.format(frame)), 'w') as vf_file:
            vf_file.write(views_text)
        
        rareas = array([g.calc_area() for g in rgeom] if svp['liparams']['cp'] == '0' else [vertarea(bm, g) for g in rgeom])
        rpos = array([g.calc_center_median() for g in rgeom] if svp['liparams']['cp'] == '0' else [g.co for g in rgeom])
        totarea = sum(rareas)
        dc1_cmd = 'rcontrib -e MF:{} -f reinhart.cal -b rbin -bn Nrbins -n {} -m sky_glow -I+ -ab 1 -ad 50000 -lw .00002 -lr -10 -faf "{}-{}.oct"'.format(cbdm_res, svp['viparams']['nproc'], svp['viparams']['filebase'], frame)
        dc8_cmd = 'rcontrib -e MF:{} -f reinhart.cal -b rbin -bn Nrbins -n {} -m sky_glow -I+ -ab 8 -ad 50000 -lw .00002 -lr -10 -faf "{}-{}.oct"'.format(cbdm_res, svp['viparams']['nproc'], svp['viparams']['filebase'], frame)
        
        with open("{}-{}.dc1".format(svp['viparams']['filebase'], frame), 'w') as dc1_file:
            with open(os.path.join(svp['viparams']['newdir'], 'views_{}.vf'.format(frame)), 'r') as vf_file:
                logentry(f'Running rcontrib with command: {dc1_cmd}')
                dc1_run = Popen(shlex.split(dc1_cmd), stdin=vf_file, stdout=dc1_file, stderr=PIPE, universal_newlines=True).communicate(input='\n'.join([r[rv].decode('utf-8') for r in rgeom]))
        
        curres += 0.33 * reslen

        if pfile.check(curres) == 'CANCELLED':
            bm.free()
            return {'CANCELLED'}
        
        with open("{}-{}.dc8".format(svp['viparams']['filebase'], frame), 'w') as dc8_file:
            with open(os.path.join(svp['viparams']['newdir'], 'views_{}.vf'.format(frame)), 'r') as vf_file:
                logentry(f'Running rcontrib with command: {dc8_cmd}')
                dc8_run = Popen(shlex.split(dc8_cmd), stdin=vf_file, stdout=dc8_file, stderr=PIPE, universal_newlines=True).communicate(input='\n'.join([r[rv].decode('utf-8') for r in rgeom]))
        
        curres += 0.33 * reslen

        if pfile.check(curres) == 'CANCELLED':
            bm.free()
            return {'CANCELLED'}

        dgp_level = '-l {:.2f}'.format(0.01 * simnode['coptions']['dgp_thresh']) if not simnode['coptions']['dgp_hourly'] else ''
        # dgc_cmds = []
        # with open(os.path.join(svp['viparams']['newdir'], 'views_{}.vf'.format(frame)), 'r') as vf_file:
        #         vf_lines = vf_file.readlines()
        #         n_lines = len(vf_lines)

        #         for ci, chunk in enumerate(chunks(vf_lines, int(svp['viparams']['nproc']))):
        #             with open(os.path.join(svp['viparams']['newdir'], 'views_{}.vf{}'.format(frame, ci)), 'w') as pvf_file:
        #                 pvf_file.write(''.join(chunk))
        #             with open(os.path.join(svp['viparams']['newdir'], 'views_{}.vf{}'.format(frame, ci)), 'r') as pvf_pfile:
        #                 # dc1_cmd = 'rcontrib -e MF:{} -f reinhart.cal -b rbin -bn Nrbins -n {} -m sky_glow -I+ -ab 1 -ad 50000 -lw .00002 -lr -10 -faf "{}-{}.oct"'.format(cbdm_res, svp['viparams']['nproc'], svp['viparams']['filebase'], frame)
        #                 # dc1_run = Popen(shlex.split(dc1_cmd), stdin=vf_pfile, stdout=dc1_file, stderr=PIPE, universal_newlines=True).communicate(input='\n'.join([r[rv].decode('utf-8') for r in rgeom]))
        #                 # dc8_cmd = 'rcontrib -e MF:{} -f reinhart.cal -b rbin -bn Nrbins -n {} -m sky_glow -I+ -ab 8 -ad 50000 -lw .00002 -lr -10 -faf "{}-{}.oct"'.format(cbdm_res, svp['viparams']['nproc'], svp['viparams']['filebase'], frame)
        #                 # dc8_run = Popen(shlex.split(dc8_cmd), stdin=vf_pfile, stdout=dc8_file, stderr=PIPE, universal_newlines=True).communicate(input='\n'.join([r[rv].decode('utf-8') for r in rgeom])) 
        #                 dgc_cmds.append('dcglare -h {4} -vf "{0}" "{1}-{3}.dc1" "{1}-{3}.dc8" "{2}"'.format(os.path.join(svp['viparams']['newdir'], 'views_{}.vf{}'.format(frame, ci)), 
        #                                                                                 svp['viparams']['filebase'],
        #                                                                                 simnode.inputs['Context in'].links[0].from_node['Options']['mtxfile'], 
        #                                                                                 frame, dgp_level))
                
        #         dgc_runs = [Popen(shlex.split(dgc_cmd), stdout=PIPE) for dgc_cmd in dgc_cmds]                                                             

        dcg_cmd = 'dcglare -h {4} -vf "{0}" "{1}-{3}.dc1" "{1}-{3}.dc8" "{2}"'.format(os.path.join(svp['viparams']['newdir'], 'views_{}.vf'.format(frame)), 
                                                                    svp['viparams']['filebase'],
                                                                    simnode.inputs['Context in'].links[0].from_node['Options']['mtxfile'], 
                                                                    frame, dgp_level)
        logentry(f'Running dcglare with command {dcg_cmd}')
        dgc_run = Popen(shlex.split(dcg_cmd), stdout=PIPE)
        # while not all([dgc_run.poll() is not None for dgc_run in dgc_runs]):
        #     time.sleep(0.2)
        
        #for ci, chunk in enumerate(chunks(rgeom, int(svp['viparams']['nproc']))):

        if not simnode['coptions']['dgp_hourly']:
            resarray = array([float(line.decode()) for line in dgc_run.stdout])
            ri = 0

            # for rg in chunk:
            for rg in rgeom:
                for vi in range(simnode['coptions']['ga_views']):
                    rg[agas[vi]] = 100 * resarray[ri]
                    rg[agos[vi]] = (1 - resarray[ri]) * 100
                    ri += 1

                # reslists.append([str(frame), 'Zone spatial', self.id_data.name, 'X', ' '.join(['{:.3f}'.format(p[0]) for p in rpos])])
                # reslists.append([str(frame), 'Zone spatial', self.id_data.name, 'Y', ' '.join(['{:.3f}'.format(p[1]) for p in rpos])])
                # reslists.append([str(frame), 'Zone spatial', self.id_data.name, 'Z', ' '.join(['{:.3f}'.format(p[2]) for p in rpos])])
                # reslists.append([str(frame), 'Zone spatial', self.id_data.name, 'Areas (m2)', ' '.join(['{:.3f}'.format(ra) for ra in rareas])])
                
                # for vi in range(simnode['coptions']['ga_views']):
                #     reslists.append([str(frame), 'Zone spatial', self.id_data.name, f'GA view {vi + 1}', ' '.join([str(rg[agas[vi]]) for rg in rgeom])])
                #     reslists.append([str(frame), 'Zone spatial', self.id_data.name, f'GO view {vi + 1}', ' '.join([str(rg[agos[vi]]) for rg in rgeom])])

        else:
            resarray = genfromtxt(dgc_run.stdout, dtype=float)

            for vi in range(simnode['coptions']['ga_views']):
                dgp_res = nsum(((resarray[vi::simnode['coptions']['ga_views']] <= 0.01 * simnode['coptions']['dgp_thresh']).T * rareas).T, axis=0)/nsum(rareas)
                reslists.append([str(frame), 'Zone temporal', self.id_data.name, f'GA view {vi + 1} (% area)', ' '.join([f'{p:.2f}' for p in 100 * dgp_res])])
                reslists.append([str(frame), 'Zone temporal', self.id_data.name, f'GO view {vi + 1} (% area)', ' '.join([f'{p:.2f}' for p in 100 * (1 - dgp_res)])])
                
                for ri, rg in enumerate(rgeom):
                    rg[agas[vi]] = 100 * sum(resarray[vi::simnode['coptions']['ga_views']][ri] <= 0.01 * simnode['coptions']['dgp_thresh'])/len(resarray[ri])
                    rg[agos[vi]] = 100 * (1 - sum(resarray[vi::simnode['coptions']['ga_views']][ri] <= 0.01 * simnode['coptions']['dgp_thresh'])/len(resarray[ri]))
                    ri += 1

        for vi in range(simnode['coptions']['ga_views']):
            self['omax']['aga{}v{}'.format(vi + 1, frame)] = max(rg[agas[vi]] for rg in rgeom)
            self['oave']['aga{}v{}'.format(vi + 1, frame)] = sum(rg[agas[vi]] for rg in rgeom)/len(rgeom)
            self['omin']['aga{}v{}'.format(vi + 1, frame)] = min(rg[agas[vi]] for rg in rgeom)
            self['omax']['ago{}v{}'.format(vi + 1, frame)] = max(rg[agos[vi]] for rg in rgeom)
            self['oave']['ago{}v{}'.format(vi + 1, frame)] = sum(rg[agos[vi]] for rg in rgeom)/len(rgeom)
            self['omin']['ago{}v{}'.format(vi + 1, frame)] = min(rg[agos[vi]] for rg in rgeom)

        self['omax']['aga1v{}'.format(frame)] = max([self['omax']['aga{}v{}'.format(vi + 1, frame)] for vi in range(simnode['coptions']['ga_views'])])
        self['oave']['aga1v{}'.format(frame)] = sum([self['oave']['aga{}v{}'.format(vi + 1, frame)] for vi in range(simnode['coptions']['ga_views'])])/simnode['coptions']['ga_views']
        self['omin']['aga1v{}'.format(frame)] = min([self['omin']['aga{}v{}'.format(vi + 1, frame)] for vi in range(simnode['coptions']['ga_views'])])
        self['omax']['ago1v{}'.format(frame)] = max([self['omax']['ago{}v{}'.format(vi + 1, frame)] for vi in range(simnode['coptions']['ga_views'])])
        self['oave']['ago1v{}'.format(frame)] = sum([self['oave']['ago{}v{}'.format(vi + 1, frame)] for vi in range(simnode['coptions']['ga_views'])])/simnode['coptions']['ga_views']
        self['omin']['ago1v{}'.format(frame)] = min([self['omin']['ago{}v{}'.format(vi + 1, frame)] for vi in range(simnode['coptions']['ga_views'])])
        svp['liparams']['views'] = simnode['coptions']['ga_views']
    
    bm.transform(self.id_data.matrix_world.inverted())
    bm.to_mesh(self.id_data.data)
    bm.free()
    return reslists



