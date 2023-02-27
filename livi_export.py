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

import bpy, os, math, subprocess, datetime, bmesh, shlex
from mathutils import Vector
from math import pi
from subprocess import PIPE, Popen, TimeoutExpired
from .vi_func import clearscene, solarPosition, retobjs, clearlayers, ct2RGB, logentry, sunpath2
from .livi_func import face_bsdf
from numpy import array, where, in1d


def radpoints(o, faces, sks):
    fentries = ['']*len(faces)
    valid_fmis = [msi for msi, ms in enumerate(o.material_slots) if o.material_slots[msi].material]
    o_mats = [o.material_slots[fmi].material for fmi in valid_fmis]
    mns = [m.name.replace(" ", "_").replace(",", "") for m in o.data.materials if m]
    on = o.name.replace(" ", "_")

    if sks:
        (skv0, skv1, skl0, skl1) = sks

    for f, face in enumerate(faces):
        if face.material_index in valid_fmis:
            fmi = valid_fmis.index(face.material_index)
            m = o_mats[fmi]
            mvp = m.vi_params
            mname = mns[fmi] if not mvp.get('bsdf') else '{}_{}_{}'.format(mns[fmi], o.name, face.index)
            mentry = face_bsdf(o, m, mname, face)
            fentry = "# Polygon \n{} polygon poly_{}_{}\n0\n0\n{}\n".format(mname, on, face.index, 3*len(face.verts))

            if sks:
                ventries = ''.join([" {0[0]:.6f} {0[1]:.6f} {0[2]:.6f}\n".format((o.matrix_world@Vector((v[skl0][0]+(v[skl1][0]-v[skl0][0])*skv1,
                                                                                                         v[skl0][1]+(v[skl1][1]-v[skl0][1])*skv1,
                                                                                                         v[skl0][2]+(v[skl1][2]-v[skl0][2])*skv1)))) for v in face.verts])
            else:
                ventries = ''.join([" {0[0]:.6f} {0[1]:.6f} {0[2]:.6f}\n".format(v.co) for v in face.verts])

            fentries[f] = ''.join((mentry, fentry, ventries+'\n'))

        else:
            logentry(f'{o.name} face {face.index} has no material defined. Check that mesh modifiers have been applied')

    return ''.join(fentries)


def bmesh2mesh(scene, obmesh, o, frame, tmf, m_export, tri):
    svp = scene.vi_params
    ftext, gradfile, vtext = '', '', ''
    bm = obmesh.copy()

    if tri:
        bmesh.ops.triangulate(bm, faces=[f for f in bm.faces if not o.material_slots[f.material_index].material.vi_params.pport])
    if not m_export:
        gradfile += radpoints(o, bm.faces, 0)
    else:
        valid_fmis = [msi for msi, ms in enumerate(o.material_slots) if o.material_slots[msi].material]
        o_mats = [o.material_slots[fmi].material for fmi in valid_fmis]
        mesh_faces = [f for f in bm.faces if f.material_index in valid_fmis]
        mrms = array([m.vi_params.radmatmenu for m in o_mats if m])
        mpps = array([not m.vi_params.pport for m in o_mats if m])
        mnpps = where(mpps, 0, 1)
        mmrms = in1d(mrms, array(('0', '1', '2', '3', '6', '9')))
        fmrms = in1d(mrms, array(('0', '1', '2', '3', '6', '7', '9')), invert=True)
        mfaces = [f for f in mesh_faces if o.material_slots[f.material_index].material and (mmrms * mpps)[valid_fmis.index(f.material_index)]]
        ffaces = [f for f in mesh_faces if o.material_slots[f.material_index].material and (fmrms + mnpps)[valid_fmis.index(f.material_index)]]
        mmats = [mat for mat in o_mats if mat and mat.vi_params.radmatmenu in ('0', '1', '2', '3', '6', '9')]
        otext = 'o {}\n'.format(o.name)
        vtext = ''.join(['v {0[0]:.6f} {0[1]:.6f} {0[2]:.6f}\n'.format(v.co) for v in bm.verts])

        if o.data.polygons and o.data.polygons[0].use_smooth:
            vtext += ''.join(['vn {0[0]:.6f} {0[1]:.6f} {0[2]:.6f}\n'.format(v.normal.normalized()) for v in bm.verts])

        if not o.data.uv_layers:
            if mfaces:
                for mat in mmats:
                    matname = mat.vi_params['radname']
                    ftext += "usemtl {}\n".format(matname) + ''.join(['f {}\n'.format(' '.join(('{0}', '{0}//{0}')[f.smooth].format(v.index + 1) for v in f.verts)) for f in mfaces if o.data.materials[f.material_index] == mat])
        else:
            uv_layer = bm.loops.layers.uv.values()[0]
            bm.faces.ensure_lookup_table()
            vtext += ''.join([''.join(['vt {0[0]} {0[1]}\n'.format(loop[uv_layer].uv) for loop in face.loops]) for face in bm.faces])
            li = 1

            for face in bm.faces:
                for loop in face.loops:
                    loop.index = li
                    li += 1

            if mfaces:
                for mat in mmats:
                    matname = mat.vi_params['radname']
                    ftext += "usemtl {}\n".format(matname) + ''.join(['f {}\n'.format(' '.join(('{0}/{1}'.format(loop.vert.index + 1, loop.index), '{0}/{1}/{0}'.format(loop.vert.index + 1, loop.index))[f.smooth] for loop in f.loops)) for f in mfaces if o.data.materials[f.material_index] == mat])

        if ffaces:
            gradfile += radpoints(o, ffaces, 0)

        if ftext:
            mfile = os.path.join(svp['viparams']['newdir'], 'obj', '{}-{}.mesh'.format(o.name.replace(' ', '_'), frame))
            ofile = os.path.join(svp['viparams']['newdir'], 'obj', '{}-{}.obj'.format(o.name.replace(' ', '_'), frame))

            with open(ofile, 'w') as objfile:
                objfile.write(otext + vtext + ftext)

            with open(mfile, 'w') as mesh:
                logentry('Running obj2mesh: {}'.format('obj2mesh -w -a {} '.format(tmf)))
                o2mrun = Popen(shlex.split("obj2mesh -w -a '{}' ".format(tmf)), stdout=mesh, stdin=PIPE, stderr=PIPE, universal_newlines=True).communicate(input=(otext + vtext + ftext))

            if os.path.getsize(mfile) and not o2mrun[1]:
                gradfile += "void mesh id \n1 '{}'\n0\n0\n\n".format(mfile)

            else:
                if o2mrun[1]:
                    logentry('Obj2mesh error: {}. Using mesh geometry export on {}. Try triangulating the Radiance mesh export'.format(o2mrun[1], o.name))

                gradfile += radpoints(o, mfaces, 0)

    bm.free()
    return gradfile


def radgexport(export_op, node):
    dp = bpy.context.evaluated_depsgraph_get()
    mats = bpy.data.materials
    scene = bpy.context.scene
    svp = scene.vi_params
    clearscene(bpy.context, export_op)
    frames = range(node['Options']['fs'], node['Options']['fe'] + 1)
    svp['liparams']['cp'] = node.cpoint
    geooblist, caloblist, lightlist = retobjs('livig'), retobjs('livic'), retobjs('livil')

    for o in caloblist:
        if any([s < 0 for s in o.scale]):
            logentry('Negative scaling on calculation object {}. Results may not be as expected'.format(o.name))
            export_op.report({'WARNING'}, 'Negative scaling on calculation object {}. Results may not be as expected'.format(o.name))

    if not node.mesh:
        for m in mats:
            if m.vi_params.radtex:
                logentry('Mesh option not selected so the texture on material {} will not be applied'.format(m.name))
                export_op.report({'WARNING'}, 'Mesh export has not been selected so the texture on material {} will not be applied'.format(m.name))

    svp['liparams']['livig'], svp['liparams']['livic'], svp['liparams']['livil'] = [o.name for o in geooblist], [o.name for o in caloblist], [o.name for o in lightlist]
    eolist = set(geooblist + caloblist)

    for o in eolist:
        ovp = o.vi_params
        ovt = ovp.vi_type

        if not node.animated:
            o.animation_data_clear()
            o.data.animation_data_clear()

        for k in [k for k in ovp.keys()]:
            del ovp[k]

        if o in caloblist:
            ovp['rtpoints'] = {}
            ovp['lisenseareas'] = {}
            ovp.vi_type_string = 'LiVi Calc'

        o.vi_params.vi_type = ovt

    for frame in frames:
        scene.frame_set(frame)
        mradfile = "".join([m.vi_params.radmat(scene) for m in mats if m])
        gradfile = "# Geometry \n\n"
        lradfile = "# Lights \n\n"
        bpy.ops.object.select_all(action='DESELECT')
        tempmatfilename = svp['viparams']['filebase']+".tempmat"

        with open(tempmatfilename, "w") as tempmatfile:
            tempmatfile.write(mradfile)

        for o in [ob for ob in eolist]:
            ovp = o.vi_params
            bm = bmesh.new()
            bm.from_object(o, dp)
            bm.transform(o.matrix_world)
            bm.normal_update()

            gradfile += bmesh2mesh(scene, bm, o, frame, tempmatfilename, node.mesh, node.triangulate)

            if o in caloblist:
                geom = (bm.faces, bm.verts)[int(node.cpoint)]

                if frame == frames[0]:
                    clearlayers(bm, 'a')
                    geom.layers.int.new('cindex')
                    o.vi_params['cpoint'] = node.cpoint

                geom.layers.string.new('rt{}'.format(frame))
                ovp.rtpoints(bm, node.offset, str(frame))
                bm.transform(o.matrix_world.inverted())
                bm.to_mesh(o.data)
                o.vi_params.vi_type_string = 'LiVi Calc'

            bm.free()

            if o.particle_systems:
                ps = o.particle_systems[0]
                hl = ps.settings.hair_length
                ps.settings.hair_length = 0
#                dp = bpy.context.evaluated_depsgraph_get()
                ps = o.evaluated_get(dp).particle_systems[0]
                particles = ps.particles
                dob = ps.settings.instance_object
                (t, r, s) = dob.matrix_world.decompose()
                dob_axis = Vector((0.0, 1.0, 0.0))
                dob_axis_glo = r@dob_axis
                dobs = [dob] if dob else []
                dobs = ps.settings.dupli_group.objects if not dobs else dobs

                for dob in dobs:
                    if not os.path.isfile(os.path.join(svp['viparams']['newdir'], 'octrees', '{}.oct'.format(dob.name.replace(' ', '_')))):
                        logentry('Octree for object {} not found in {}'.format(dob.name, os.path.join(svp['viparams']['newdir'], 'octrees')))
#                        gen_octree(scene, dob, export_op, node.fallback)
#                        if node.fallback or dob.hide_get(): # Should check visibility with dob.hide_get() but it's not working
                    else:
                        # dovp = dob.vi_params
                        # bm = bmesh.new()
                        # tempmesh = dob.evaluated_get(depsgraph).to_mesh()
                        # bm.from_mesh(tempmesh)
                        # bm.transform(dob.matrix_world)
                        # bm.normal_update()
                        # dob.to_mesh_clear()
                        # bmesh2mesh(scene, bm, dob, frame, tempmatfilename, 0, node.triangulate)
                        # bm.free()

                        for p, part in enumerate(particles):
                            nv = Vector((1, 0, 0))
                            nv.rotate(part.rotation)
                            rdiff = dob_axis_glo.rotation_difference(nv).to_euler()
                            gradfile += 'void instance {7}\n17 {6} -t {2[0]:.4f} {2[1]:.4f} {2[2]:.4f} -s {4:.3f} -rx {5[0]:.4f} -ry {5[1]:.4f} -rz {5[2]:.4f} -t {3[0]:.4f} {3[1]:.4f} {3[2]:.4f} \n0\n0\n\n'.format(dob.name,
                                        p, [-p for p in dob.location], part.location, part.size, [180.0 * r/math.pi for r in (rdiff.x, rdiff.y, rdiff.z)],
                                        os.path.join(svp['viparams']['newdir'], 'octrees', '{}.oct'.format(dob.name.replace(' ', '_'), frame)), '{}_copy_{}'.format(o.name, p))

                o.particle_systems[0].settings.hair_length = hl

    # Lights export routine
        for o in [ob for ob in lightlist if ob.visible_get()]:
            ovp = o.vi_params

            # if ' ' in bpy.data.filepath:
            #     logentry('There is a space in the Blender file name or directory path - re-save with no spaces in the filename/directory path')
            #     export_op.report({'ERROR'}, 'There is a space in the Blender file name or directory path - re-save with no spaces in the filename/directory path')
            # elif ' ' in ovp.ies_name:
            #     logentry('There is a space in the {} IES file name or directory path - move/rename it'.format(o.name))
            #     export_op.report({'ERROR'}, 'There is a space in the {} IES file name or directory path - move/rename it'.format(o.name))
            # else:
            ab_ies_path = bpy.path.abspath(ovp.ies_name)
            iesname = os.path.splitext(os.path.basename(ab_ies_path))[0]

            if os.path.isfile(ab_ies_path):
                iescmd = "ies2rad -t default -m {0} -c {1[0]:.4f} {1[1]:.4f} {1[2]:.4f} -p '{2}' -d{3} -o '{4}-{5}' '{6}'".format(ovp.ies_strength, (ovp.ies_rgb, ct2RGB(ovp.ies_ct))[ovp.ies_colmenu == '1'], svp['liparams']['lightfilebase'], ovp.ies_unit, iesname, frame, ab_ies_path)
                logentry('Running ies2rad with command: {}'.format(iescmd))
                subprocess.call(shlex.split(iescmd))

                with open(os.path.join(svp['liparams']['lightfilebase'], '{}-{}.rad'.format(iesname, frame)), 'r') as dat_file:
                    dat_str = dat_file.read()
                    dat_str = dat_str.replace(os.path.join(svp['liparams']['lightfilebase'], '{}-{}.dat'.format(iesname, frame)), '"{}"'.format(os.path.join(svp['liparams']['lightfilebase'], '{}-{}.dat'.format(iesname, frame))))

                    for suf in (f'-{frame}_dist', f'-{frame}_light', f'-{frame}.u', f'-{frame}.s'):
                        dat_str = dat_str.replace(iesname+suf, f'"{iesname}{suf}"')

                    dat_str = dat_str.replace(f' {iesname}-{frame}.d', f' "{iesname}-{frame}.d"')

                with open(os.path.join(svp['liparams']['lightfilebase'], '{}-{}.rad'.format(iesname, frame)), 'w') as dat_file:
                    dat_file.write(dat_str)

                if o.type == 'LIGHT':
                    if o.parent:
                        o = o.parent
                    lradfile += u'!xform -rx {0[0]:.4f} -ry {0[1]:.4f} -rz {0[2]:.4f} -t {1[0]:.4f} {1[1]:.4f} {1[2]:.4f} "{2}.rad"\n\n'.format([(180/pi)*o.rotation_euler[i] for i in range(3)], o.location, os.path.join(svp['liparams']['lightfilebase'], iesname+"-{}".format(frame)))

                elif o.type == 'MESH':
                    tm = o.to_mesh()
                    bm = bmesh.new()
                    bm.from_mesh(tm)
                    o.to_mesh_clear()
                    bm.transform(o.matrix_world)
                    bm.normal_update()
                    bm.faces.ensure_lookup_table()

                    for f in bm.faces:
                        lrot = Vector.rotation_difference(Vector((0, 0, -1)), f.normal).to_euler('XYZ')
                        lradfile += u'!xform -rx {0[0]:.4f} -ry {0[1]:.4f} -rz {0[2]:.4f} -t {1[0]:.4f} {1[1]:.4f} {1[2]:.4f} "{2}"{3}'.format([math.degrees(lr) for lr in lrot], f.calc_center_median(), os.path.join(svp['liparams']['lightfilebase'], iesname+"-{}.rad".format(frame)), ('\n', '\n\n')[f == bm.faces[-1]])
                    bm.free()

            elif iesname:
                export_op.report({'ERROR'}, 'The IES file associated with {} cannot be found'.format(o.name))

        sradfile = "# Sky \n\n"
        node['Text'][str(frame)] = mradfile+gradfile+lradfile+sradfile


def gen_octree(scene, o, op, mesh, tri):
    dg = bpy.context.evaluated_depsgraph_get()
    bm = bmesh.new()
    tempmesh = o.evaluated_get(dg).to_mesh()
    bm.from_mesh(tempmesh)
    bm.transform(o.matrix_world)
    bm.normal_update()
    o.to_mesh_clear()
    nd = scene.vi_params['viparams']['newdir']

    if not os.path.isdir(os.path.join(nd, 'octrees')):
        os.makedirs(os.path.join(nd, 'octrees'))

    mats = o.data.materials
    mradfile = "".join([m.vi_params.radmat(scene) for m in mats if m])
    bpy.ops.object.select_all(action='DESELECT')
    mf = os.path.join(nd, 'octrees', '{}.mat'.format(o.name))

    with open(mf, "w") as tempmatfile:
        tempmatfile.write(mradfile)

    gradfile = bmesh2mesh(scene, bm, o, scene.frame_current, mf, mesh, tri)

    with open(os.path.join(nd, 'octrees', '{}.oct'.format(o.name)), "wb") as octfile:
        try:
            ocrun =  Popen("oconv -w -".split(), stdin=PIPE, stderr=PIPE, stdout=octfile, universal_newlines=True)
            err = ocrun.communicate(input=mradfile + gradfile, timeout=600)[1]

            if err:
                logentry('Oconv conversion error: {}'.format(err))
                return 'CANCELLED'

        except TimeoutExpired:
            ocrun.kill()
            errmsg = 'Oconv conversion taking too long. Try joining/simplfying geometry or using non-mesh geometry export'
            op.report({'ERROR'}, errmsg)
            logentry('Oconv error: {}'.format(errmsg))
            return 'CANCELLED'

    for line in err:
        logentry('Oconv error: {}'.format(line))
    if err and 'fatal -' in err:
        op.report({'ERROR'}, 'Oconv conversion failure: {}'.format(err))
        return 'CANCELLED'
    elif err and 'set overflow' in err:
        op.report({'ERROR'}, 'Ratio of largest to smallest geometry is too large. Clean up mesh geometry or decrease the radius of any HDR panorama')
        return 'CANCELLED'

def livi_sun(scene, node, frame):
    svp = scene.vi_params

    if node.skyprog in ('0', '1') and node.contextmenu == 'Basic':
        simtime = node.starttime + frame*datetime.timedelta(seconds=3600*node.interval)
        solalt, solazi, beta, phi = solarPosition(simtime.timetuple()[7], simtime.hour + (simtime.minute)*0.016666, svp.latitude, svp.longitude)

        if node.skyprog == '0':
            gscmd = "gensky -ang {:.3f} {:.3f} {} -t {} -g {}".format(solalt, solazi, node['skytypeparams'], node.turb, node.gref)
        else:
            gscmd = "gendaylit -ang {:.3f} {:.3f} {} -g {}".format(solalt, solazi, node['skytypeparams'], node.gref)
    else:
        gscmd = "gensky -ang {:.3f} {:.3f} {} -g {}".format(45, 0, node['skytypeparams'], node.gref)

    logentry('Generating sky with the command: {}'.format(gscmd))
    gsrun = Popen(gscmd.split(), stdout=PIPE)
    return gsrun.stdout.read().decode()

def hdrexport(scene, f, frame, node, skytext):
    svp = scene.vi_params

    with open('{}-{}sky.oct'.format(svp['viparams']['filebase'], frame), 'w') as skyoct:
        Popen('oconv -w -'.split(), stdin=PIPE, stdout=skyoct).communicate(input=skytext.encode('utf-8'))

    with open(os.path.join(svp['viparams']['newdir'], str(frame)+".hdr"), 'w') as hdrfile:
        rpictcmd = 'rpict -vta -vp 0 0 0 -vd 0 1 0 -vu 0 0 1 -vh 360 -vv 360 -x 1500 -y 1500 "{}-{}sky.oct"'.format(svp['viparams']['filebase'], frame)
        Popen(shlex.split(rpictcmd), stdout=hdrfile).communicate()

    cntrun = Popen('cnt 750 1500'.split(), stdout=PIPE)
    rcalccmd = 'rcalc -f "{}" -e XD=1500;YD=750;inXD=0.000666;inYD=0.001333'.format(os.path.join(svp.vipath, 'RadFiles', 'lib', 'latlong.cal'))
    rcalcrun = Popen(shlex.split(rcalccmd), stdin=cntrun.stdout, stdout=PIPE)
    rtracecmd = 'rtrace -n {} -x 1500 -y 750 -fac "{}-{}sky.oct"'.format(svp['viparams']['nproc'], svp['viparams']['filebase'], frame)

    with open('{}p.hdr'.format(os.path.join(svp['viparams']['newdir'], str(frame))), 'w') as hdrim:
        Popen(shlex.split(rtracecmd), stdin=rcalcrun.stdout, stdout=hdrim).communicate()

    if '{}p.hdr'.format(frame) not in bpy.data.images:
        bpy.data.images.load(os.path.join(svp['viparams']['newdir'], "{}p.hdr".format(frame)))
    else:
        bpy.data.images['{}p.hdr'.format(frame)].reload()

def livi_sky(sn):
    skytext = "4 .8 .8 1 0\n\n" if sn < 3 else "4 1 1 1 0\n\n"
    return "\nskyfunc glow sky_glow\n0\n0\n" + skytext + "sky_glow source sky\n0\n0\n4 0 0 1  180\n\n"

def livi_ground(r, g, b, ref):
    fac = ref/(r * 0.265 + g * 0.670 + b * 0.065)
    if ref:
        return "skyfunc glow ground_glow\n0\n0\n4 {0[0]:.3f} {0[1]:.3f} {0[2]:.3f} 0\n\nground_glow source ground\n0\n0\n4 0 0 -1 180\n\n".format([c*fac for c in (r, g, b)])
    else:
        return ''

def createradfile(scene, frame, export_op, simnode):
    radtext = ''
    links = (list(simnode.inputs['Geometry in'].links[:]) + list(simnode.inputs['Context in'].links[:]))

    for link in links:
        if str(frame) in link.from_node['Text']:
            radtext += link.from_node['Text'][str(frame)]
        elif frame < min([int(k) for k in link.from_node['Text'].keys()]):
            radtext += link.from_node['Text'][str(min([int(k) for k in link.from_node['Text'].keys()]))]
        elif frame > max([int(k) for k in link.from_node['Text'].keys()]):
            radtext += link.from_node['Text'][str(max([int(k) for k in link.from_node['Text'].keys()]))]

    simnode['radfiles'][str(frame)] = radtext

def createoconv(scene, frame, sim_op, simnode, **kwargs):
    svp = scene.vi_params
    fbase = "{0}-{1}".format(svp['viparams']['filebase'], frame)

    with open("{}.oct".format(fbase), "wb") as octfile:
        try:
            ocrun = Popen("oconv -w -".split(), stdin=PIPE, stderr=PIPE, stdout=octfile, universal_newlines=True)
            err = ocrun.communicate(input=simnode['radfiles'][str(frame)], timeout=600)[1]

            if err:
                logentry('Oconv conversion error: {}'.format(err))
                return 'CANCELLED'

        except TimeoutExpired:
            ocrun.kill()
            errmsg = 'Oconv conversion taking too long. Try joining/simplfying geometry or using non-mesh geometry export'
            sim_op.report({'ERROR'}, errmsg)
            logentry('Oconv error: {}'.format(errmsg))
            return 'CANCELLED'

    for line in err:
        logentry('Oconv error: {}'.format(line))
    if err and 'fatal -' in err:
        sim_op.report({'ERROR'}, 'Oconv conversion failure: {}'.format(err))
        return 'CANCELLED'
    elif err and 'set overflow' in err:
        sim_op.report({'ERROR'}, 'Ratio of largest to smallest geometry is too large. Clean up mesh geometry or decrease the radius of any HDR panorama')
        return 'CANCELLED'


def spfc(self):
    sunpath2()


def cyfc1(self):
    scene = bpy.context.scene
    svp = scene.vi_params

    if 'LiVi' in svp['viparams']['resnode'] or 'Shadow' in svp['viparams']['resnode']:
        for material in [m for m in bpy.data.materials if m.use_nodes and m.vi_params.mattype == '1']:
            try:
                if any([node.bl_label == 'Attribute' for node in material.node_tree.nodes]):
                    material.node_tree.nodes["Attribute"].attribute_name = str(scene.frame_current)
            except Exception as e:
                print(e, 'Something wrong with changing the material attribute name')
