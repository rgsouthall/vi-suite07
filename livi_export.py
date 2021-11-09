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

import bpy, os, math, subprocess, datetime, bmesh, mathutils, shlex
from math import sin, cos, pi
from subprocess import PIPE, Popen, TimeoutExpired
from .vi_func import clearscene, solarPosition, retobjs, clearlayers, viparams, ct2RGB, logentry, sunpath2
from .livi_func import face_bsdf
from numpy import array, where, in1d

def radpoints(o, faces, sks):
    fentries = ['']*len(faces) 
    mns = [m.name.replace(" ", "_").replace(",", "") for m in o.data.materials]
    on = o.name.replace(" ", "_")
    
    if sks:
        (skv0, skv1, skl0, skl1) = sks

    for f, face in enumerate(faces):
        fmi = face.material_index
        m = o.data.materials[fmi]
        mvp = m.vi_params
        mname = mns[fmi] if not mvp.get('bsdf') else '{}_{}_{}'.format(mns[fmi], o.name, face.index)
        mentry = face_bsdf(o, m, mname, face)
        fentry = "# Polygon \n{} polygon poly_{}_{}\n0\n0\n{}\n".format(mname, on, face.index, 3*len(face.verts))
        
        if sks:
            ventries = ''.join([" {0[0]:.6f} {0[1]:.6f} {0[2]:.6f}\n".format((o.matrix_world@mathutils.Vector((v[skl0][0]+(v[skl1][0]-v[skl0][0])*skv1, v[skl0][1]+(v[skl1][1]-v[skl0][1])*skv1, v[skl0][2]+(v[skl1][2]-v[skl0][2])*skv1)))) for v in face.verts])
        else:
            ventries = ''.join([" {0[0]:.6f} {0[1]:.6f} {0[2]:.6f}\n".format(v.co) for v in face.verts])
        
        fentries[f] = ''.join((mentry, fentry, ventries+'\n'))        
    return ''.join(fentries)

def bmesh2mesh(scene, obmesh, o, frame, tmf, m, tri):
    svp = scene.vi_params
    ftext, gradfile, vtext = '', '', ''
    bm = obmesh.copy()

    if tri:
        bmesh.ops.triangulate(bm, faces = [f for f in bm.faces if not o.material_slots[f.material_index].material.vi_params.pport])
    if not m:
#        gradfile += radpoints(o, [f for f in bm.faces if f.calc_area() >= 0.0], 0)
        gradfile += radpoints(o, bm.faces, 0)
#        for f in [f for f in bm.faces if f.calc_area() == 0.0]:
#            logentry('Object {} face {} has not been exported as it is too small'.format(o.name, f.index))
    else:
        mrms = array([m.vi_params.radmatmenu for m in o.data.materials])
        mpps = array([not m.vi_params.pport for m in o.data.materials])        
        mnpps = where(mpps, 0, 1)        
        mmrms = in1d(mrms, array(('0', '1', '2', '3', '6', '9')))        
        fmrms = in1d(mrms, array(('0', '1', '2', '3', '6', '7', '9')), invert = True)
        mfaces = [f for f in bm.faces if (mmrms * mpps)[f.material_index]]
        ffaces = [f for f in bm.faces if (fmrms + mnpps)[f.material_index]]        
        mmats = [mat for mat in o.data.materials if mat.vi_params.radmatmenu in ('0', '1', '2', '3', '6', '9')]
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
                    li +=1
                    
            if mfaces:
                for mat in mmats:
                    matname = mat.vi_params['radname']
                    ftext += "usemtl {}\n".format(matname) + ''.join(['f {}\n'.format(' '.join(('{0}/{1}'.format(loop.vert.index + 1, loop.index), '{0}/{1}/{0}'.format(loop.vert.index + 1, loop.index))[f.smooth]  for loop in f.loops)) for f in mfaces if o.data.materials[f.material_index] == mat])
              
        if ffaces:
            gradfile += radpoints(o, ffaces, 0)
    
        if ftext:   
            mfile = os.path.join(svp['viparams']['newdir'], 'obj', '{}-{}.mesh'.format(o.name.replace(' ', '_'), frame))
            ofile = os.path.join(svp['viparams']['newdir'], 'obj', '{}-{}.obj'.format(o.name.replace(' ', '_'), frame))
            
            with open(ofile, 'w') as objfile:
                objfile.write(otext + vtext + ftext)
                        
            with open(mfile, 'w') as mesh:
                logentry('Running obj2mesh: {}'.format('obj2mesh -w -a {} '.format(tmf)))
                o2mrun = Popen('obj2mesh -w -a {} '.format(tmf).split(), stdout = mesh, stdin = PIPE, stderr = PIPE, universal_newlines=True).communicate(input = (otext + vtext + ftext))
                
                               
            if os.path.getsize(mfile) and not o2mrun[1]:
                gradfile += "void mesh id \n1 {}\n0\n0\n\n".format(mfile)
    
            else:
                if o2mrun[1]:
                    logentry('Obj2mesh error: {}. Using mesh geometry export on {}. Try triangulating the Radiance mesh export'.format(o2mrun[1], o.name))
    
                gradfile += radpoints(o, mfaces, 0)
#    print(gradfile)
           
    bm.free()       
    return gradfile

def radgexport(export_op, node, **kwargs):
    dp = bpy.context.evaluated_depsgraph_get()
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
            
    svp['liparams']['livig'], svp['liparams']['livic'], svp['liparams']['livil'] = [o.name for o in geooblist], [o.name for o in caloblist], [o.name for o in lightlist]
    eolist = set(geooblist + caloblist)
    mats = bpy.data.materials

    for o in eolist:  
        ovp = o.vi_params
        ovt = ovp.vi_type

        if not node.animated:
            o.animation_data_clear()
            o.data.animation_data_clear()

        for k in [k for k in ovp.keys()]:
            del ovp[k]
        
        if o in caloblist:
            o.vi_params['rtpoints'] = {}
            o.vi_params['lisenseareas'] = {}
        
        o.vi_params.vi_type = ovt

    for frame in frames:
        scene.frame_set(frame)
        mradfile =  "".join([m.vi_params.radmat(scene) for m in mats if m])
        gradfile = "# Geometry \n\n"
        lradfile = "# Lights \n\n"
        bpy.ops.object.select_all(action='DESELECT')
        tempmatfilename = svp['viparams']['filebase']+".tempmat"
        
        with open(tempmatfilename, "w") as tempmatfile:
            tempmatfile.write(mradfile)

        for o in eolist:
            ovp = o.vi_params
            bm = bmesh.new()
            bm.from_object(o, dp)
#            tempmesh = o.evaluated_get(dp).to_mesh()
#            bm.from_mesh(tempmesh)
            bm.transform(o.matrix_world)
            bm.normal_update() 
            o.to_mesh_clear()
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
                dob_axis = mathutils.Vector((0.0, 1.0, 0.0))
                dob_axis_glo = r@dob_axis
                dobs = [dob] if dob else []
                dobs = ps.settings.dupli_group.objects if not dobs else dobs
                
                for dob in dobs: 
                    if not os.path.isfile(os.path.join(svp['viparams']['newdir'], 'octrees', '{}.oct'.format(dob.name.replace(' ', '_')))):
                        logentry('Octree for object {} not found in {}'.format(dob.name, os.path.join(svp['viparams']['newdir'], 'octrees')))
#                        gen_octree(scene, dob, export_op, node.fallback)
                    
#                        if node.fallback or dob.hide_get(): # Should check visibility with dob.hide_get() but it's not working
                    else:
#                        dovp = dob.vi_params
                        # bm = bmesh.new()
                        # tempmesh = dob.evaluated_get(depsgraph).to_mesh()
                        # bm.from_mesh(tempmesh)
                        # bm.transform(dob.matrix_world)
                        # bm.normal_update() 
                        # dob.to_mesh_clear()
                        # bmesh2mesh(scene, bm, dob, frame, tempmatfilename, 0, node.triangulate)
                        # bm.free()
        
                        for p, part in enumerate(particles):
                           nv = mathutils.Vector((1, 0, 0))
                           nv.rotate(part.rotation)
                           rdiff = dob_axis_glo.rotation_difference(nv).to_euler()
                           gradfile += 'void instance {7}\n17 {6} -t {2[0]:.4f} {2[1]:.4f} {2[2]:.4f} -s {4:.3f} -rx {5[0]:.4f} -ry {5[1]:.4f} -rz {5[2]:.4f} -t {3[0]:.4f} {3[1]:.4f} {3[2]:.4f} \n0\n0\n\n'.format(dob.name, 
                                        p, [-p for p in dob.location], part.location, part.size, [180.0 * r/math.pi for r in (rdiff.x, rdiff.y, rdiff.z)], 
                                        os.path.join(svp['viparams']['newdir'], 'octrees', '{}.oct'.format(dob.name.replace(' ', '_'), frame)), '{}_copy_{}'.format(o.name, p))

                o.particle_systems[0].settings.hair_length = hl

    # Lights export routine        
        for o in lightlist:
            ovp = o.vi_params
            
            if ' ' in ovp.ies_name:
                logentry('There is a space in the {} IES file name - rename it'.format(o.name))
                export_op.report({'ERROR'}, 'There is a space in the {} IES file name - rename it'.format(o.name))
            else:
                ab_ies_path = bpy.path.abspath(ovp.ies_name)
                iesname = os.path.splitext(os.path.basename(ab_ies_path))[0]

                if os.path.isfile(ab_ies_path):
                    iescmd = "ies2rad -t default -m {0} -c {1[0]:.4f} {1[1]:.4f} {1[2]:.4f} -p '{2}' -d{3} -o {4}-{5} '{6}'".format(ovp.ies_strength, (ovp.ies_rgb, ct2RGB(ovp.ies_ct))[ovp.ies_colmenu == '1'], svp['liparams']['lightfilebase'], ovp.ies_unit, iesname, frame, ab_ies_path)
                    logentry('Running ies2rad with command: {}'.format(iescmd))
                    subprocess.call(shlex.split(iescmd))
    
                    if o.type == 'LIGHT':
                        if o.parent:
                            o = o.parent
                        lradfile += u'!xform -rx {0[0]:.4f} -ry {0[1]:.4f} -rz {0[2]:.4f} -t {1[0]:.4f} {1[1]:.4f} {1[2]:.4f} "{2}.rad"\n\n'.format([(180/pi)*o.rotation_euler[i] for i in range(3)], o.location, os.path.join(svp['liparams']['lightfilebase'], iesname+"-{}".format(frame)))
    
                    elif o.type == 'MESH':
                        tm = o.to_mesh()
                        bm = bmesh.new()
                        bm.from_mesh(tm)
                        o.to_mesh.clear()
                        bm.transform(o.matrix_world)
                        bm.normal_update()
                        bm.faces.ensure_lookup_table()
                        
                        for f in bm.faces: 
                            lrot = mathutils.Vector.rotation_difference(mathutils.Vector((0, 0, -1)), f.normal).to_euler('XYZ')
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
    mradfile =  "".join([m.vi_params.radmat(scene) for m in mats if m])        
    bpy.ops.object.select_all(action='DESELECT')
    mf = os.path.join(nd, 'octrees', '{}.mat'.format(o.name))
        
    with open(mf, "w") as tempmatfile:
        tempmatfile.write(mradfile)  
        
    gradfile = bmesh2mesh(scene, bm, o, scene.frame_current, mf, mesh, tri)
    print(gradfile + mradfile)
    with open(os.path.join(nd, 'octrees', '{}.oct'.format(o.name)), "wb") as octfile:
        try:
            ocrun =  Popen("oconv -w -".split(), stdin = PIPE, stderr = PIPE, stdout = octfile, universal_newlines=True)
            err = ocrun.communicate(input = mradfile + gradfile, timeout = 600)[1]

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
        simtime = node.starttime + frame*datetime.timedelta(seconds = 3600*node.interval)
        solalt, solazi, beta, phi = solarPosition(simtime.timetuple()[7], simtime.hour + (simtime.minute)*0.016666, svp.latitude, svp.longitude)
        
        if node.skyprog == '0':
            gscmd = "gensky -ang {} {} {} -t {} -g {}".format(solalt, solazi, node['skytypeparams'], node.turb, node.gref)            
        else:
            gscmd = "gendaylit -ang {} {} {} -g {}".format(solalt, solazi, node['skytypeparams'], node.gref)
    else:
        gscmd = "gensky -ang {} {} {} -g {}".format(45, 0, node['skytypeparams'], node.gref)
        
    logentry('Generating sky with the command: {}'.format(gscmd))
    gsrun = Popen(gscmd.split(), stdout = PIPE) 
    return gsrun.stdout.read().decode()

def hdrexport(scene, f, frame, node, skytext):
    svp = scene.vi_params
    
    with open('{}-{}sky.oct'.format(svp['viparams']['filebase'], frame), 'w') as skyoct:
        Popen('oconv -w -'.split(), stdin = PIPE, stdout = skyoct).communicate(input = skytext.encode('utf-8'))

    with open(os.path.join(svp['viparams']['newdir'], str(frame)+".hdr"), 'w') as hdrfile:
        rpictcmd = 'rpict -vta -vp 0 0 0 -vd 0 1 0 -vu 0 0 1 -vh 360 -vv 360 -x 1500 -y 1500 "{}-{}sky.oct"'.format(svp['viparams']['filebase'], frame)
        Popen(shlex.split(rpictcmd), stdout = hdrfile).communicate()

    cntrun = Popen('cnt 750 1500'.split(), stdout = PIPE)
    rcalccmd = 'rcalc -f "{}" -e XD=1500;YD=750;inXD=0.000666;inYD=0.001333'.format(os.path.join(svp.vipath, 'RadFiles', 'lib', 'latlong.cal'))
    rcalcrun = Popen(shlex.split(rcalccmd), stdin = cntrun.stdout, stdout = PIPE)
    rtracecmd = 'rtrace -n {} -x 1500 -y 750 -fac "{}-{}sky.oct"'.format(svp['viparams']['nproc'], svp['viparams']['filebase'], frame)

    with open('{}p.hdr'.format(os.path.join(svp['viparams']['newdir'], str(frame))), 'w') as hdrim:
        Popen(shlex.split(rtracecmd), stdin = rcalcrun.stdout, stdout = hdrim).communicate()

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
            ocrun =  Popen("oconv -w -".split(), stdin = PIPE, stderr = PIPE, stdout = octfile, universal_newlines=True)
            err = ocrun.communicate(input = simnode['radfiles'][str(frame)], timeout = 600)[1]

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
    # scene = bpy.context.scene
    # svp = scene.vi_params

    # if not svp['viparams'].get('newframe'):
    #     svp['viparams']['newframe'] = 1
    # else:
    #     svp['viparams']['newframe'] = 0
    #     scene.frame_set(scene.frame_current)
        
    # if svp['viparams']['resnode'] == 'VI Sun Path':
    #     spoblist = {ob.get('VIType'):ob for ob in scene.objects if ob.get('VIType') in ('Sun', 'SPathMesh')}
    #     beta, phi = solarPosition(svp.sp_sd, svp.sp_sh, svp.latitude, svp.longitude)[2:]
    #     print(beta, phi)
    #     if scene.world.use_nodes == False:
    #         scene.world.use_nodes = True
    #     nt = bpy.data.worlds[0].node_tree

    #     if nt and nt.nodes.get('Sky Texture'):
    #         scene.world.node_tree.nodes['Sky Texture'].sun_direction = -sin(phi), -cos(phi), sin(beta)

    #     for ob in scene.objects:
    #         if ob.get('VIType') == 'Sun':
    #             ob.rotation_euler = pi * 0.5 - beta, 0, -phi 
    #             ob.location.z = spoblist['SPathMesh'].location.z + 100 * sin(beta)                
    #             ob.location.x = spoblist['SPathMesh'].location.x -(100**2 - (spoblist['Sun'].location.z-spoblist['SPathMesh'].location.z)**2)**0.5 * sin(phi)
    #             ob.location.y = spoblist['SPathMesh'].location.y -(100**2 - (spoblist['Sun'].location.z-spoblist['SPathMesh'].location.z)**2)**0.5 * cos(phi)
                
    #             if ob.data.node_tree:
    #                 for blnode in [blnode for blnode in ob.data.node_tree.nodes if blnode.bl_label == 'Blackbody']:
    #                     blnode.inputs[0].default_value = 2500 + 3000*sin(beta)**0.5
    #                 for emnode in [emnode for emnode in ob.data.node_tree.nodes if emnode.bl_label == 'Emission']:
    #                     emnode.inputs[1].default_value = 10 * sin(beta)

    #         elif ob.get('VIType') == 'SunMesh':
    #             ob.location = (0, 0, 0)
    #             if ob.data.materials[0].node_tree:
    #                 for smblnode in [smblnode for smblnode in ob.data.materials[0].node_tree.nodes if ob.data.materials and smblnode.bl_label == 'Blackbody']:
    #                     smblnode.inputs[0].default_value = 2500 + 3000*sin(beta)**0.5
    # else:
    #     return
    
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
        
# def genbsdf(scene, export_op, o):     
#     if viparams(export_op, scene):
#         return
    
#     svp = scene.vi_params
#     bsdfmats = [mat for mat in o.data.materials if mat.radmatmenu == '8']

#     if bsdfmats:
#         mat = bsdfmats[0]
#         mat['bsdf'] = {} 
#     else:
#         export_op.report({'ERROR'}, '{} does not have a BSDF material attached'.format(o.name))
    
#     tm = o.to_mesh(scene = scene, apply_modifiers = True, settings = 'PREVIEW')
#     bm = bmesh.new()    
#     bm.from_mesh(tm) 
#     bpy.data.meshes.remove(tm)
#     bm.transform(o.matrix_world)
#     bm.normal_update()
#     bsdffaces = [face for face in bm.faces if o.data.materials[face.material_index].radmatmenu == '8']    
    
#     if bsdffaces:
#         fvec = bsdffaces[0].normal
#         mat['bsdf']['normal'] = '{0[0]:.5f} {0[1]:.5f} {0[2]:.5f}'.format(fvec)
#         mat['bsdf']['pos'] == '{0[0]:.5f} {0[1]:.5f} {0[2]:.5f}'.format(bsdffaces[0].calc_centre_median())
#     else:
#         export_op.report({'ERROR'}, '{} does not have a BSDF material associated with any faces'.format(o.name))
#         return
    
#     zvec, xvec = mathutils.Vector((0, 0, 1)), mathutils.Vector((1, 0, 0))
#     svec = mathutils.Vector.cross(fvec, zvec)
#     bm.faces.ensure_lookup_table()
#     bsdfrotz = mathutils.Matrix.Rotation(mathutils.Vector.angle(fvec, zvec), 4, svec)
#     bm.transform(bsdfrotz)
#     bsdfrotx = mathutils.Matrix.Rotation(math.pi + mathutils.Vector.angle_signed(mathutils.Vector(xvec[:2]), mathutils.Vector(svec[:2])), 4, zvec)#mathutils.Vector.cross(svec, xvec))
#     bm.transform(bsdfrotx)
#     vposis = list(zip(*[v.co[:] for v in bm.verts]))
#     (maxx, maxy, maxz) = [max(p) for p in vposis]
#     (minx, miny, minz) = [min(p) for p in vposis]
#     bsdftrans = mathutils.Matrix.Translation(mathutils.Vector((-(maxx + minx)/2, -(maxy + miny)/2, -maxz)))
#     bm.transform(bsdftrans)
#     mradfile = ''.join([m.radmat(scene) for m in o.data.materials if m.radmatmenu != '8'])                  
#     gradfile = radpoints(o, [face for face in bm.faces if o.data.materials and face.material_index < len(o.data.materials) and o.data.materials[face.material_index].radmatmenu != '8'], 0)
#     bm.free()  
#     bsdfsamp = o.li_bsdf_ksamp if o.li_bsdf_tensor == ' ' else 2**(int(o.li_bsdf_res) * 2) * int(o.li_bsdf_tsamp) 
#     gbcmd = "genBSDF +geom {} -r '{}' {} {} -c {} {} -n {}".format(o.li_bsdf_dimen,  o.li_bsdf_rcparam,  o.li_bsdf_tensor, (o.li_bsdf_res, ' ')[o.li_bsdf_tensor == ' '], bsdfsamp, o.li_bsdf_direc, svp['viparams']['nproc'])
#     mat['bsdf']['xml'] = Popen(shlex.split(gbcmd), stdin = PIPE, stdout = PIPE).communicate(input = (mradfile+gradfile).encode('utf-8'))[0].decode()
#     svp['viparams']['vidisp'] = 'bsdf'
#     mat['bsdf']['type'] = o.li_bsdf_tensor
    
        