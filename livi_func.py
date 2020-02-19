import mathutils, bpy, bmesh, os, datetime, shlex, sys, math
from subprocess import Popen, PIPE, STDOUT
from numpy import array, where, in1d, transpose, savetxt, int8, float16, float32, float64, digitize, zeros, choose, inner, average, amax, amin
from numpy import sum as nsum
from numpy import max as nmax
from numpy import min as nmin
from numpy import mean as nmean
from numpy import append as nappend
#from numpy import delete as ndelete
from .vi_func import vertarea, logentry, ct2RGB, clearlayers, chunks, selobj

def radpoints(o, faces, sks):
    fentries = ['']*len(faces) 
    mns = [m.name.replace(" ", "_").replace(",", "") for m in o.data.materials]
    on = o.name.replace(" ", "_")
    
    if sks:
        (skv0, skv1, skl0, skl1) = sks

    for f, face in enumerate(faces):
        fmi = face.material_index
        mname = mns[fmi]
        fentry = "# Polygon \n{} polygon poly_{}_{}\n0\n0\n{}\n".format(mname, on, face.index, 3*len(face.verts))
        if sks:
            ventries = ''.join([" {0[0]:.4f} {0[1]:.4f} {0[2]:.4f}\n".format((o.matrix_world*mathutils.Vector((v[skl0][0]+(v[skl1][0]-v[skl0][0])*skv1, v[skl0][1]+(v[skl1][1]-v[skl0][1])*skv1, v[skl0][2]+(v[skl1][2]-v[skl0][2])*skv1)))) for v in face.verts])
        else:
            ventries = ''.join([" {0[0]:.4f} {0[1]:.4f} {0[2]:.4f}\n".format(v.co) for v in face.verts])
        fentries[f] = ''.join((fentry, ventries+'\n'))        
    return ''.join(fentries)

def radbsdf(self, radname, fi, rot, trans):
    fmat = self.data.materials[self.data.polygons[fi].material_index]
    pdepth = fmat['bsdf']['proxy_depth'] if self.bsdf_proxy else 0 
    bsdfxml = self.data.materials[self.data.polygons[fi].material_index]['bsdf']['xml']
    radname = '{}_{}_{}'.format(fmat.name, self.name, fi)
    radentry = 'void BSDF {0}\n16 {4:.4f} {1} 0 0 1 . -rx {2[0]:.4f} -ry {2[1]:.4f} -rz {2[2]:.4f} -t {3[0]:.4f} {3[1]:.4f} {3[2]:.4f}\n0\n0\n\n'.format(radname, bsdfxml, rot, trans, pdepth)
    return radentry
           
def rtpoints(self, bm, offset, frame):    
    geom = bm.verts if self['cpoint'] == '1' else bm.faces 
    cindex = geom.layers.int['cindex']
    rt = geom.layers.string['rt{}'.format(frame)]
    for gp in geom:
        gp[cindex] = 0 
    geom.ensure_lookup_table()
    resfaces = [face for face in bm.faces if self.id_data.data.materials[face.material_index].vi_params.mattype == '1']
    self['cfaces'] = [face.index for face in resfaces]
       
    if self['cpoint'] == '0': 
        gpoints = resfaces
        gpcos =  [gp.calc_center_median_weighted() for gp in gpoints]
        self['cverts'], self['lisenseareas'][frame] = [], [f.calc_area() for f in gpoints]       

    elif self['cpoint'] == '1': 
        gis = sorted(set([item.index for sublist in [face.verts[:] for face in resfaces] for item in sublist]))
        gpoints = [geom[gi] for gi in gis]
        gpcos = [gp.co for gp in gpoints]
        self['cverts'], self['lisenseareas'][frame] = gp.index, [vertarea(bm, gp) for gp in gpoints]    
    
    for g, gp in enumerate(gpoints):
        gp[rt] = '{0[0]:.4f} {0[1]:.4f} {0[2]:.4f} {1[0]:.4f} {1[1]:.4f} {1[2]:.4f}'.format([gpcos[g][i] + offset * gp.normal.normalized()[i] for i in range(3)], gp.normal[:]).encode('utf-8')
        gp[cindex] = g + 1
        
    self['rtpnum'] = g + 1
    
def validradparams(params):
    valids = ('-ps', '-pt', '-pj', '-pj', '-dj', '-ds', '-dt', '-dc', '-dr', '-dp',	'-ss',	'-st',	'-st', '-ab',	'-av',	'-aa',	'-ar',	'-ad',	'-as',	'-lr',	'-lw')
    for p, param in enumerate(params.split()):
        if not p%2 and param not in valids:
            return 0
        elif  p%2:
            try: float(param)
            except: return 0   
    return 1    

def bmesh2mesh(scene, obmesh, o, frame, tmf, fb):
    svp = scene.vi_params
    ftext, gradfile, vtext = '', '', ''

#    try:
    bm = obmesh.copy()
    bmesh.ops.remove_doubles(bm, verts = bm.verts, dist = 0.0001)
    bmesh.ops.dissolve_limit(bm, angle_limit = 0.0001, use_dissolve_boundaries = False, verts = bm.verts, edges = bm.edges, delimit = {'NORMAL'})
    bmesh.ops.connect_verts_nonplanar(bm, angle_limit = 0.0001, faces = bm.faces)
    mrms = array([m.vi_params.radmatmenu for m in o.data.materials])
    mpps = array([not m.vi_params.pport for m in o.data.materials])        
    mnpps = where(mpps, 0, 1)        
    mmrms = in1d(mrms, array(('0', '1', '2', '3', '6', '9')))        
    fmrms = in1d(mrms, array(('0', '1', '2', '3', '6', '7', '9')), invert = True)
    mfaces = [f for f in bm.faces if (mmrms * mpps)[f.material_index]]
    ffaces = [f for f in bm.faces if (fmrms + mnpps)[f.material_index]]        
    mmats = [mat for mat in o.data.materials if mat.vi_params.radmatmenu in ('0', '1', '2', '3', '6', '9')]
#    for mm in mmats:
#        if not mm.get()
    otext = 'o {}\n'.format(o.name)
    vtext = ''.join(['v {0[0]:.6f} {0[1]:.6f} {0[2]:.6f}\n'.format(v.co) for v in bm.verts])
    
    if o.data.polygons[0].use_smooth:
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
        
        with open(mfile, 'w') as mesh:
            o2mrun = Popen('obj2mesh -w -a {} '.format(tmf).split(), stdout = mesh, stdin = PIPE, stderr = PIPE, universal_newlines=True).communicate(input = (otext + vtext + ftext))
                           
        if os.path.getsize(mfile) and not o2mrun[1] and not fb:
            gradfile += "void mesh id \n1 {}\n0\n0\n\n".format(mfile)

        else:
            if o2mrun[1]:
                logentry('Obj2mesh error: {}. Using geometry export fallback on {}'.format(o2mrun[1], o.name))

            gradfile += radpoints(o, mfaces, 0)

        with open(ofile, 'w') as objfile:
            objfile.write(otext + vtext + ftext)

    bm.free()
        
    return gradfile
    
#    except Exception as e:
#        logentry('LiVi mesh export error for {}: {}'.format(o.name, e))
#        return gradfile
    
def radmat(self, scene):
#    mvp = self.vi_params
    svp = scene.vi_params
    radname = self.id_data.name.replace(" ", "_")
    radname = radname.replace(",", "")
    print(self.id_data.name)
    self['radname'] = radname
    radtex = ''
    mod = 'void' 
    
    if self.radmatmenu in ('0', '1', '2', '3', '6') and self.radtex:
        try:
            teximage = self.node_tree.nodes['Material Output'].inputs['Surface'].links[0].from_node.inputs['Color'].links[0].from_node.image
            teximageloc = os.path.join(svp['liparams']['texfilebase'],'{}.hdr'.format(radname))
            off = scene.render.image_settings.file_format 
            scene.render.image_settings.file_format = 'HDR'
            teximage.save_render(teximageloc, scene)
            scene.render.image_settings.file_format = off
            (w, h) = teximage.size
            ar = ('*{}'.format(w/h), '') if w >= h else ('', '*{}'.format(h/w))
            radtex = 'void colorpict {}_tex\n7 red green blue {} . frac(Lu){} frac(Lv){}\n0\n0\n\n'.format(radname, '{}'.format(teximageloc), ar[0], ar[1])
            mod = '{}_tex'.format(radname)
            
            try:
                if self.radnorm:             
                    normimage = self.id_data.node_tree.nodes['Material Output'].inputs['Surface'].links[0].from_node.inputs['Normal'].links[0].from_node.inputs['Color'].links[0].from_node.image
                    header = '2\n0 1 {}\n0 1 {}\n'.format(normimage.size[1], normimage.size[0])
                    xdat = -1 + 2 * array(normimage.pixels[:][0::4]).reshape(normimage.size[0], normimage.size[1])
                    ydat = -1 + 2 * array(normimage.pixels[:][1::4]).reshape(normimage.size[0], normimage.size[1])# if self.gup == '0' else 1 - 2 * array(normimage.pixels[:][1::4]).reshape(normimage.size[0], normimage.size[1])
                    savetxt(os.path.join(svp['liparams']['texfilebase'],'{}.ddx'.format(radname)), xdat, fmt='%.2f', header = header, comments='')
                    savetxt(os.path.join(svp['liparams']['texfilebase'],'{}.ddy'.format(radname)), ydat, fmt='%.2f', header = header, comments='')
                    radtex += "{0}_tex texdata {0}_norm\n9 ddx ddy ddz {1}.ddx {1}.ddy {1}.ddy nm.cal frac(Lv){2} frac(Lu){3}\n0\n7 {4} {5[0]} {5[1]} {5[2]} {6[0]} {6[1]} {6[2]}\n\n".format(radname, os.path.join(svp['viparams']['newdir'], 'textures', radname), ar[1], ar[1], self.ns, self.nu, self.nside)
                    mod = '{}_norm'.format(radname)
                    
            except Exception as e:
                print('Problem with normal export {}'.format(e))
                
        except Exception as e:
            print('Problem with texture export {}'.format(e))
         
    radentry = '# ' + ('plastic', 'glass', 'dielectric', 'translucent', 'mirror', 'light', 'metal', 'antimatter', 'bsdf', 'custom')[int(self.radmatmenu)] + ' material\n' + \
            '{} {} {}\n'.format(mod, ('plastic', 'glass', 'dielectric', 'trans', 'mirror', 'light', 'metal', 'antimatter', 'bsdf', 'custom')[int(self.radmatmenu)], radname) + \
           {'0': '0\n0\n5 {0[0]:.3f} {0[1]:.3f} {0[2]:.3f} {1:.3f} {2:.3f}\n'.format(self.radcolour, self.radspec, self.radrough), 
            '1': '0\n0\n3 {0[0]:.3f} {0[1]:.3f} {0[2]:.3f}\n'.format(self.radcolour), 
            '2': '0\n0\n5 {0[0]:.3f} {0[1]:.3f} {0[2]:.3f} {1:.3f} 0\n'.format(self.radcolour, self.radior),
            '3': '0\n0\n7 {0[0]:.3f} {0[1]:.3f} {0[2]:.3f} {1:.3f} {2:.3f} {3:.3f} {4:.3f}\n'.format(self.radcolour, self.radspec, self.radrough, self.radtrans, self.radtranspec), 
            '4': '0\n0\n3 {0[0]:.3f} {0[1]:.3f} {0[2]:.3f}\n'.format(self.radcolour),
            '5': '0\n0\n3 {0[0]:.3f} {0[1]:.3f} {0[2]:.3f}\n'.format([c * self.radintensity for c in (self.radcolour, ct2RGB(self.radct))[self.radcolmenu == '1']]), 
            '6': '0\n0\n5 {0[0]:.3f} {0[1]:.3f} {0[2]:.3f} {1:.3f} {2:.3f}\n'.format(self.radcolour, self.radspec, self.radrough), 
            '7': '1 void\n0\n0\n', '8': '1 void\n0\n0\n', '9': '1 void\n0\n0\n'}[self.radmatmenu] + '\n'

    if self.radmatmenu == '8' and self.get('bsdf') and self['bsdf'].get('xml'):
        bsdfxml = os.path.join(svp['viparams']['newdir'], 'bsdfs', '{}.xml'.format(radname))
        
        
        with open(bsdfxml, 'w') as bsdffile:
            bsdffile.write(self['bsdf']['xml'])
        radentry = 'void BSDF {0}\n6 {1:.4f} {2} 0 0 1 .\n0\n0\n\n'.format(radname, self.li_bsdf_proxy_depth, bsdfxml)
        
    elif self.radmatmenu == '9':
        radentry = bpy.data.texts[self.radfile].as_string()+'\n\n' if self.radfile in [t.name for t in bpy.data.texts] else '# dummy material\nvoid plastic {}\n0\n0\n5 0.8 0.8 0.8 0.1 0.1\n\n'.format(radname)
                        
    self['radentry'] = radtex + radentry
    return(radtex + radentry)    
    
def cbdmmtx(self, scene, locnode, export_op):
    svp = scene.vi_params
    os.chdir(svp['viparams']['newdir'])  
     
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
                elif self.cbdm_start_hour <= float(ls[2]) <= self.cbdm_end_hour and self.sdoy <= datetime.datetime(2015, int(ls[0]), int(ls[1])).timetuple().tm_yday <= self.edoy:
                    weafile.write(line)
                
        gdmcmd = ("gendaymtx -m 1 {} {}".format(('', '-O1')[self['watts']], 
                  "{0}.wea".format(os.path.join(svp['viparams']['newdir'], self['epwbase'][0]))))
        with open("{}.mtx".format(os.path.join(svp['viparams']['newdir'], self['epwbase'][0])), 'w') as mtxfile:
            Popen(gdmcmd.split(), stdout = mtxfile, stderr=STDOUT).communicate()
        with open("{}-whitesky.oct".format(svp['viparams']['filebase']), 'w') as wsfile:
            oconvcmd = "oconv -w -"
            Popen(shlex.split(oconvcmd), stdin = PIPE, stdout = wsfile).communicate(input = self['whitesky'].encode(sys.getfilesystemencoding()))
        return "{}.mtx".format(os.path.join(svp['viparams']['newdir'], self['epwbase'][0]))
    else:
        export_op.report({'ERROR'}, "Not a valid EPW file")
        return ''
    
def cbdmhdr(node, scene):
    svp = scene.vi_params
    targethdr = os.path.join(svp['viparams']['newdir'], node['epwbase'][0]+"{}.hdr".format(('l', 'w')[node['watts']]))
    latlonghdr = os.path.join(svp['viparams']['newdir'], node['epwbase'][0]+"{}p.hdr".format(('l', 'w')[node['watts']]))
    skyentry = hdrsky(node.hdrname, '1', 0, 1000) if node.sourcemenu == '1' and  node.cbanalysismenu == '0' else hdrsky(targethdr, '1', 0, 1000)

    if node.sourcemenu != '1' or node.cbanalysismenu == '2':
        vecvals, vals = mtx2vals(open(node['mtxfile'], 'r').readlines(), datetime.datetime(2015, 1, 1).weekday(), node, node.times)
        pcombfiles = ''.join(["{} ".format(os.path.join(svp['viparams']['newdir'], 'ps{}.hdr'.format(i))) for i in range(146)])
        vwcmd = "vwrays -ff -x 600 -y 600 -vta -vp 0 0 0 -vd 0 1 0 -vu 0 0 1 -vh 360 -vv 360 -vo 0 -va 0 -vs 0 -vl 0"
        rcontribcmd = "rcontrib -bn 146 -fo -ab 0 -ad 1 -n {} -ffc -x 600 -y 600 -ld- -V+ -f tregenza.cal -b tbin -o {} -m sky_glow {}-whitesky.oct".format(svp['viparams']['nproc'], 
                                                           os.path.join(svp['viparams']['newdir'], 'p%d.hdr'), 
                                                           os.path.join(svp['viparams']['newdir'], 
                                                                        svp['viparams']['filename']))
        vwrun = Popen(vwcmd.split(), stdout = PIPE)
        rcrun = Popen(rcontribcmd.split(), stderr = PIPE, stdin = vwrun.stdout)
        for line in rcrun.stderr:
            logentry('HDR generation error: {}'.format(line))
    
        for j in range(146):
            with open(os.path.join(svp['viparams']['newdir'], "ps{}.hdr".format(j)), 'w') as psfile:
                Popen("pcomb -s {} {}".format(vals[j], os.path.join(svp['viparams']['newdir'], 'p{}.hdr'.format(j))).split(), stdout = psfile).wait()
        with open(targethdr, 'w') as epwhdr:
            Popen("pcomb -h {}".format(pcombfiles).split(), stdout = epwhdr).wait()
        
        [os.remove(os.path.join(svp['viparams']['newdir'], 'p{}.hdr'.format(i))) for i in range (146)]
        [os.remove(os.path.join(svp['viparams']['newdir'], 'ps{}.hdr'.format(i))) for i in range (146)]
        node.hdrname = targethdr
    
        if node.hdr:
            with open('{}.oct'.format(os.path.join(svp['viparams']['newdir'], node['epwbase'][0])), 'w') as hdroct:
                Popen(shlex.split("oconv -w - "), stdin = PIPE, stdout=hdroct, stderr=STDOUT).communicate(input = skyentry.encode(sys.getfilesystemencoding()))
            cntrun = Popen('cnt 750 1500'.split(), stdout = PIPE)
            rcalcrun = Popen('rcalc -f {} -e XD=1500;YD=750;inXD=0.000666;inYD=0.001333'.format(os.path.join(svp.vipath, 'Radfiles', 'lib', 'latlong.cal')).split(), stdin = cntrun.stdout, stdout = PIPE)
            with open(latlonghdr, 'w') as panohdr:
                rtcmd = 'rtrace -n {} -x 1500 -y 750 -fac {}.oct'.format(svp['viparams']['nproc'], os.path.join(svp['viparams']['newdir'], node['epwbase'][0]))
                Popen(rtcmd.split(), stdin = rcalcrun.stdout, stdout = panohdr)
    return skyentry

def mtx2vals(mtxlines, fwd, node, times):    
    for m, mtxline in enumerate(mtxlines):
        if 'NROWS' in mtxline:
            patches = int(mtxline.split('=')[1])
            
        elif mtxline == '\n':
            startline = m + 1
            break

#    sdoy = (times[0] - datetime.datetime(2015, 1, 1)).days
#    shour = times[0].hour
#    edoy = (times[-1] - datetime.datetime(2015, 1, 1)).days + 1
#    ehour = times[-1].hour
    tothours = len(times)
    hours = [t.hour for t in times]
    
#    invalidhours = [h for h in range(8760) if h < sdoy * 24 or h > edoy  * 24 or h%24 < shour or h%24 > ehour] 
    mtxlarray = array([0.333 * sum([float(lv) for lv in fval.split(" ")]) for fval in mtxlines[startline:] if fval != '\n'], dtype=float)
    mtxshapearray = mtxlarray.reshape(patches, int(len(mtxlarray)/patches))
#    mtxshapearray = ndelete(mtxshapearray, invalidhours, 1)
    vals = nsum(mtxshapearray, axis = 1)
    vvarray = transpose(mtxshapearray)
    vvlist = vvarray.tolist()
    vecvals = [[hours[x], (fwd+int(hours[x]/24))%7, *vvlist[x]] for x in range(tothours)]
    return(vecvals, vals)
    
def hdrsky(hdrfile, hdrmap, hdrangle, hdrradius):
    hdrangle = '1 {:.3f}'.format(hdrangle * math.pi/180) if hdrangle else '1 0'
    hdrfn = {'0': 'sphere2latlong', '1': 'sphere2angmap'}[hdrmap]
    return("# Sky material\nvoid colorpict hdr_env\n7 red green blue '{}' {}.cal sb_u sb_v\n0\n{}\n\nhdr_env glow env_glow\n0\n0\n4 1 1 1 0\n\nenv_glow bubble sky\n0\n0\n4 0 0 0 {}\n\n".format(hdrfile, hdrfn, hdrangle, hdrradius))

def retpmap(node, frame, scene):
    svp = scene.vi_params
    pportmats = ' '.join([mat.name.replace(" ", "_") for mat in bpy.data.materials if mat.vi_params.pport and mat.vi_params.get('radentry')])
    ammats = ' '.join([mat.name.replace(" ", "_") for mat in bpy.data.materials if mat.vi_params.mattype == '1' and mat.vi_params.radmatmenu == '7' and mat.vi_params.get('radentry')])
    pportentry = ' '.join(['-apo {}'.format(ppm) for ppm in pportmats.split()]) if pportmats else ''
    amentry = '-aps {}'.format(ammats) if ammats else ''
    cpentry = '-apc {}-{}.cpm {}'.format(svp['viparams']['filebase'], frame, node.pmapcno) if node.pmapcno else ''
    cpfileentry = '-ap {}-{}.cpm 50'.format(svp['viparams']['filebase'], frame) if node.pmapcno else ''  
    return amentry, pportentry, cpentry, cpfileentry   

def retsv(self, scene, frame, rtframe, chunk, rt):
    svcmd = "rcontrib -w -I -n {} {} -m sky_glow {}-{}.oct ".format(scene.vi_params['viparams']['nproc'], '-ab 1 -ad 8192 -aa 0 -ar 512 -as 1024 -lw 0.0002 ', scene.vi_params['viparams']['filebase'], frame)    
    rtrun = Popen(svcmd.split(), stdin = PIPE, stdout=PIPE, stderr=STDOUT, universal_newlines=True).communicate(input = '\n'.join([c[rt].decode('utf-8') for c in chunk]))                
    reslines = nsum(array([[float(rv) for rv in r.split('\t')[:3]] for r in rtrun[0].splitlines()[10:]]), axis = 1)
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
    geom = bm.verts if self['cpoint'] == '1' else bm.faces
    cindex = geom.layers.int['cindex']
    totarea = sum([gp.calc_area() for gp in geom if gp[cindex] > 0]) if self['cpoint'] == '0' else sum([vertarea(bm, gp) for gp in geom])
    
    for f, frame in enumerate(frames):
        self['res{}'.format(frame)] = {}
        if svp['liparams']['unit'] in ('DF (%)', 'Lux'):
            geom.layers.float.new('virrad{}'.format(frame))
            geom.layers.float.new('illu{}'.format(frame))
            virradres = geom.layers.float['virrad{}'.format(frame)]
            illures = geom.layers.float['illu{}'.format(frame)]
        if svp['liparams']['unit'] == 'DF (%)':
            geom.layers.float.new('df{}'.format(frame))
            dfres = geom.layers.float['df{}'.format(frame)]
        elif svp['liparams']['unit'] == 'W/m2 (f)':
            geom.layers.float.new('firrad{}'.format(frame))
            firradres = geom.layers.float['firrad{}'.format(frame)]

        geom.layers.float.new('res{}'.format(frame))
        res =  geom.layers.float['res{}'.format(frame)]
        
        if geom.layers.string.get('rt{}'.format(frame)):
            rtframe = frame
        else:
            kints = [int(k[2:]) for k in geom.layers.string.keys()]
            rtframe = max(kints) if frame > max(kints) else min(kints)
        
        rt =  geom.layers.string['rt{}'.format(rtframe)]
            
        for chunk in chunks([g for g in geom if g[rt]], int(svp['viparams']['nproc']) * 500):
            rtrun = Popen(rtcmds[f].split(), stdin = PIPE, stdout=PIPE, stderr=PIPE, universal_newlines=True).communicate(input = '\n'.join([c[rt].decode('utf-8') for c in chunk]))   
            xyzirrad = array([[float(v) for v in sl.split('\t')[:3]] for sl in rtrun[0].splitlines()])
            if svp['liparams']['unit'] == 'W/m2 (f)':
                firrad = nsum(xyzirrad * array([0.333, 0.333, 0.333]), axis = 1)
            elif svp['liparams']['unit'] in ('DF (%)', 'Lux'):
                virrad = nsum(xyzirrad * array([0.26, 0.67, 0.065]), axis = 1)
                illu = virrad * 179
#            firrad = virrad * 1.64
                if svp['liparams']['unit'] == 'DF (%)':
                    df = illu * 0.01
            
            for gi, gp in enumerate(chunk):                
                if svp['liparams']['unit'] == 'W/m2 (f)':
                    gp[firradres] = firrad[gi].astype(float32)
                    gp[res] = firrad[gi].astype(float32)
                elif svp['liparams']['unit'] in ('DF (%)', 'Lux'):   
                    gp[virradres] = virrad[gi].astype(float32)
                    gp[illures] = illu[gi].astype(float32)                    
                    if svp['liparams']['unit'] == 'DF':
                        gp[dfres] = df[gi].astype(float16)
                    gp[res] = illu[gi].astype(float32)
                
            curres += len(chunk)
            if pfile.check(curres) == 'CANCELLED':
                bm.free()
                return {'CANCELLED'}

        oirrad = array([g[virradres] for g in geom]).astype(float64) if svp['liparams']['unit'] in ('DF (%)', 'Lux') else array([g[firradres] for g in geom]).astype(float64)
        maxoirrad, minoirrad, aveoirrad = nmax(oirrad), nmin(oirrad), nmean(oirrad)
        self['livires'][str(frame)] = (maxoirrad, minoirrad, aveoirrad)
        if svp['liparams']['unit'] == 'W/m2 (f)':
            self['omax']['firrad{}'.format(frame)] =  maxoirrad
            self['oave']['firrad{}'.format(frame)] = aveoirrad            
            self['omin']['firrad{}'.format(frame)] = minoirrad
            if self['omax']['firrad{}'.format(frame)] > self['omin']['firrad{}'.format(frame)]:
                vals = [(gp[res] - self['omin']['firrad{}'.format(frame)])/(self['omax']['firrad{}'.format(frame)] - self['omin']['firrad{}'.format(frame)]) for gp in geom]
            else:
                vals = [1 for gp in geom]
        elif svp['liparams']['unit'] in ('DF (%)', 'Lux'):
            self['omax']['virrad{}'.format(frame)] = maxoirrad
            self['omax']['illu{}'.format(frame)] =  maxoirrad * 179
            self['omin']['virrad{}'.format(frame)] = minoirrad
            self['oave']['illu{}'.format(frame)] = aveoirrad * 179
            self['oave']['virrad{}'.format(frame)] = aveoirrad
            self['omin']['illu{}'.format(frame)] = minoirrad * 179
            
            if self['omax']['illu{}'.format(frame)] > self['omin']['illu{}'.format(frame)]:
                vals = [(gp[res] - self['omin']['illu{}'.format(frame)])/(self['omax']['illu{}'.format(frame)] - self['omin']['illu{}'.format(frame)]) for gp in geom]
            else:
                vals = [1 for gp in geom]
                
            if svp['liparams']['unit'] == 'DF':        
                self['omax']['df{}'.format(frame)] =  maxoirrad * 1.79
                self['oave']['df{}'.format(frame)] = aveoirrad * 1.79       
                self['omin']['df{}'.format(frame)] = minoirrad * 1.79
        
        tableheaders = [["", 'Minimum', 'Average', 'Maximum']]
        posis = [v.co for v in bm.verts if v[cindex] > 0] if self['cpoint'] == '1' else [f.calc_center_bounds() for f in bm.faces if f[cindex] > 1]
#        illubinvals = [self['omin']['illu{}'.format(frame)] + (self['omax']['illu{}'.format(frame)] - self['omin']['illu{}'.format(frame)])/ll * (i + increment) for i in range(ll)]
        bins = array([increment * i for i in range(1, ll)])
        
#        if self['omax']['illu{}'.format(frame)] > self['omin']['illu{}'.format(frame)]:
#            vals = [(gp[res] - self['omin']['illu{}'.format(frame)])/(self['omax']['illu{}'.format(frame)] - self['omin']['illu{}'.format(frame)]) for gp in geom]         
        
            
        ais = digitize(vals, bins)
        rgeom = [g for g in geom if g[cindex] > 0]
        rareas = [gp.calc_area() for gp in geom] if self['cpoint'] == '0' else [vertarea(bm, gp) for gp in geom]
        sareas = zeros(ll)
        
        for ai in range(ll):
            sareas[ai] = sum([rareas[gi]/totarea for gi in range(len(rgeom)) if ais[gi] == ai])
                
        self['livires']['areabins'] = sareas
        
        reslists.append([str(frame), 'Zone', self.id_data.name, 'X', ' '.join(['{:.3f}'.format(p[0]) for p in posis])])
        reslists.append([str(frame), 'Zone', self.id_data.name, 'Y', ' '.join(['{:.3f}'.format(p[1]) for p in posis])])
        reslists.append([str(frame), 'Zone', self.id_data.name, 'Z', ' '.join(['{:.3f}'.format(p[2]) for p in posis])])
        reslists.append([str(frame), 'Zone', self.id_data.name, 'Areas (m2)', ' '.join(['{:.3f}'.format(ra) for ra in rareas])])
        
        if svp['liparams']['unit'] == 'W/m2 (f)':
            firradbinvals = [self['omin']['firrad{}'.format(frame)] + (self['omax']['firrad{}'.format(frame)] - self['omin']['firrad{}'.format(frame)])/ll * (i + increment) for i in range(ll)]
            self['livires']['valbins'] = firradbinvals
            self['tablefi{}'.format(frame)] = array(tableheaders + [['Full Irradiance (W/m2)', '{:.1f}'.format(self['omin']['firrad{}'.format(frame)]), '{:.1f}'.format(self['oave']['firrad{}'.format(frame)]), '{:.1f}'.format(self['omax']['firrad{}'.format(frame)])]])
            reslists.append([str(frame), 'Zone', self.id_data.name, 'Full Irradiance (W/m2)', ' '.join(['{:.3f}'.format(g[firradres]) for g in rgeom])])

        elif svp['liparams']['unit'] in ('DF (%)', 'Lux'):
            illubinvals = [self['omin']['illu{}'.format(frame)] + (self['omax']['illu{}'.format(frame)] - self['omin']['illu{}'.format(frame)])/ll * (i + increment) for i in range(ll)]
            self['livires']['valbins'] = illubinvals
            self['tableillu{}'.format(frame)] = array(tableheaders + [['Illuminance (lux)', 
            '{:.1f}'.format(self['omin']['illu{}'.format(frame)]), '{:.1f}'.format(self['oave']['illu{}'.format(frame)]), 
            '{:.1f}'.format(self['omax']['illu{}'.format(frame)])]])
            self['tablevi{}'.format(frame)] = array(tableheaders + [['Visual Irradiance (W/m2)', '{:.1f}'.format(self['omin']['virrad{}'.format(frame)]), '{:.1f}'.format(self['oave']['virrad{}'.format(frame)]), '{:.1f}'.format(self['omax']['virrad{}'.format(frame)])]])
            reslists.append([str(frame), 'Zone', self.id_data.name, 'Illuminance (lux)', ' '.join(['{:.3f}'.format(g[illures]) for g in rgeom])])
            reslists.append([str(frame), 'Zone', self.id_data.name, 'Visible Irradiance (W/m2)', ' '.join(['{:.3f}'.format(g[virradres]) for g in rgeom])])

            if svp['liparams']['unit'] == 'DF': 
                dfbinvals = [self['omin']['df{}'.format(frame)] + (self['omax']['df{}'.format(frame)] - self['omin']['df{}'.format(frame)])/ll * (i + increment) for i in range(ll)]
                self['livires']['valbins'] = dfbinvals
                self['tabledf{}'.format(frame)] = array(tableheaders + [['DF (%)', '{:.1f}'.format(self['omin']['df{}'.format(frame)]), '{:.1f}'.format(self['oave']['df{}'.format(frame)]), '{:.1f}'.format(self['omax']['df{}'.format(frame)])]])
                reslists.append([str(frame), 'Zone', self.id_data.name, 'DF (%)', ' '.join(['{:.3f}'.format(g[dfres]) for g in rgeom])])

    if len(frames) > 1:
        reslists.append(['All', 'Frames', '', 'Frames', ' '.join([str(f) for f in frames])])
        if svp['liparams']['unit'] == 'W/m2 (f)':
            reslists.append(['All', 'Zone', self.id_data.name, 'Average irradiance (W/m2)', ' '.join(['{:.3f}'.format(self['oave']['firrad{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone', self.id_data.name, 'Maximum irradiance (W/m2)', ' '.join(['{:.3f}'.format(self['omax']['firrad{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone', self.id_data.name, 'Minimum irradiance (W/m2)', ' '.join(['{:.3f}'.format(self['omin']['firrad{}'.format(frame)]) for frame in frames])])
        elif svp['liparams']['unit'] in ('DF (%)', 'Lux'):            
            reslists.append(['All', 'Zone', self.id_data.name, 'Average illuminance (lux)', ' '.join(['{:.3f}'.format(self['oave']['illu{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone', self.id_data.name, 'Maximum illuminance (lux)', ' '.join(['{:.3f}'.format(self['omax']['illu{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone', self.id_data.name, 'Minimum illuminance (lux)', ' '.join(['{:.3f}'.format(self['omin']['illu{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone', self.id_data.name, 'Average irradiance (W/m2)', ' '.join(['{:.3f}'.format(self['oave']['virrad{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone', self.id_data.name, 'Maximum irradiance (W/m2)', ' '.join(['{:.3f}'.format(self['omax']['virrad{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone', self.id_data.name, 'Minimum irradiance (W/m2)', ' '.join(['{:.3f}'.format(self['omin']['virrad{}'.format(frame)]) for frame in frames])])
            reslists.append(['All', 'Zone', self.id_data.name, 'Illuminance ratio', ' '.join(['{:.3f}'.format(self['omin']['illu{}'.format(frame)]/self['oave']['illu{}'.format(frame)]) for frame in frames])])
            if svp['liparams']['unit'] == 'DF': 
                reslists.append(['All', 'Zone', self.id_data.name, 'Average DF (lux)', ' '.join(['{:.3f}'.format(self['oave']['df{}'.format(frame)]) for frame in frames])])
                reslists.append(['All', 'Zone', self.id_data.name, 'Maximum DF (lux)', ' '.join(['{:.3f}'.format(self['omax']['df{}'.format(frame)]) for frame in frames])])
                reslists.append(['All', 'Zone', self.id_data.name, 'Minimum DF (lux)', ' '.join(['{:.3f}'.format(self['omin']['df{}'.format(frame)]) for frame in frames])])
   
    bm.transform(self.id_data.matrix_world.inverted())
    bm.to_mesh(self.id_data.data)
    bm.free()
    return reslists
    
def lhcalcapply(self, scene, frames, rtcmds, simnode, curres, pfile):
    reslists = []
    svp = scene.vi_params
    bm = bmesh.new()
    bm.from_mesh(self.id_data.data)
    self['omax'], self['omin'], self['oave'] = {}, {}, {}
    clearlayers(bm, 'f')
    geom = bm.verts if self['cpoint'] == '1' else bm.faces
    cindex = geom.layers.int['cindex']
    
    for f, frame in enumerate(frames): 
        geom.layers.float.new('firradm2{}'.format(frame))
        geom.layers.float.new('virradm2{}'.format(frame))
        geom.layers.float.new('firrad{}'.format(frame))
        geom.layers.float.new('virrad{}'.format(frame))
        geom.layers.float.new('illu{}'.format(frame))
        geom.layers.float.new('res{}'.format(frame))
        firradm2res = geom.layers.float['firradm2{}'.format(frame)]
        virradm2res = geom.layers.float['virradm2{}'.format(frame)]
        firradres = geom.layers.float['firrad{}'.format(frame)]
        virradres = geom.layers.float['virrad{}'.format(frame)]
        illures = geom.layers.float['illu{}'.format(frame)]
         
        if geom.layers.string.get('rt{}'.format(frame)):
            rtframe = frame
        else:
            kints = [int(k[2:]) for k in geom.layers.string.keys()]
            rtframe  = max(kints) if frame > max(kints) else  min(kints)
        
        rt = geom.layers.string['rt{}'.format(rtframe)]
        gps = [g for g in geom if g[rt]]
        areas = array([g.calc_area() for g in gps] if self['cpoint'] == '0' else [vertarea(bm, g) for g in gps])

        for chunk in chunks(gps, int(svp['viparams']['nproc']) * 200):
            careas = array([c.calc_area() if self['cpoint'] == '0' else vertarea(bm, c) for c in chunk])
            rtrun = Popen(rtcmds[f].split(), stdin = PIPE, stdout=PIPE, stderr=PIPE, universal_newlines=True).communicate(input = '\n'.join([c[rt].decode('utf-8') for c in chunk]))   
            xyzirrad = array([[float(v) for v in sl.split('\t')[:3]] for sl in rtrun[0].splitlines()])
            virradm2 = nsum(xyzirrad * array([0.26, 0.67, 0.065]), axis = 1) * 1e-3
            virrad = virradm2 * careas
            firradm2 = virradm2 * 1.64
            firrad = firradm2 * careas
            illu = virradm2 * 179e-3

            for gi, gp in enumerate(chunk):
                gp[firradm2res] = firradm2[gi].astype(float32)
                gp[virradm2res] = virradm2[gi].astype(float32)
                gp[firradres] = firrad[gi].astype(float32)
                gp[virradres] = virrad[gi].astype(float32)
                gp[illures] = illu[gi].astype(float32)
            
            curres += len(chunk)
            if pfile.check(curres) == 'CANCELLED':
                bm.free()
                return {'CANCELLED'}
                
        ovirradm2 = array([g[virradm2res] for g in gps])
        ovirrad = array([g[virradres] for g in gps])
        maxovirradm2 = nmax(ovirradm2)
        maxovirrad = nmax(ovirrad)
        minovirradm2 = nmin(ovirradm2)
        minovirrad = nmin(ovirrad)
        aveovirradm2 = nmean(ovirradm2)
        aveovirrad = nmean(ovirrad)
        self['omax']['firrad{}'.format(frame)] = maxovirrad * 1.64
        self['omin']['firrad{}'.format(frame)] = minovirrad * 1.64
        self['oave']['firrad{}'.format(frame)] = aveovirrad * 1.64
        self['omax']['firradm2{}'.format(frame)] = maxovirradm2  * 1.64
        self['omin']['firradm2{}'.format(frame)] = minovirradm2  * 1.64
        self['oave']['firradm2{}'.format(frame)] = aveovirradm2  * 1.64
        self['omax']['virrad{}'.format(frame)] = maxovirrad
        self['omin']['virrad{}'.format(frame)] = minovirrad
        self['oave']['virrad{}'.format(frame)] = aveovirrad
        self['omax']['virradm2{}'.format(frame)] = maxovirradm2
        self['omin']['virradm2{}'.format(frame)] = minovirradm2
        self['oave']['virradm2{}'.format(frame)] = aveovirradm2
        self['omax']['illu{}'.format(frame)] = maxovirradm2 * 178e-3
        self['omin']['illu{}'.format(frame)] = minovirradm2 * 178e-3
        self['oave']['illu{}'.format(frame)] = aveovirradm2 * 178e-3
        self['tablemlxh{}'.format(frame)] = array([["", 'Minimum', 'Average', 'Maximum'], 
            ['Luxhours (Mlxh)', '{:.1f}'.format(self['omin']['illu{}'.format(frame)]), '{:.1f}'.format(self['oave']['illu{}'.format(frame)]), '{:.1f}'.format(self['omax']['illu{}'.format(frame)])]])
        self['tablefim2{}'.format(frame)] = array([["", 'Minimum', 'Average', 'Maximum'], 
            ['Full Irradiance (kWh/m2)', '{:.1f}'.format(self['omin']['firradm2{}'.format(frame)]), '{:.1f}'.format(self['oave']['firradm2{}'.format(frame)]), '{:.1f}'.format(self['omax']['firradm2{}'.format(frame)])]])
        self['tablevim2{}'.format(frame)] = array([["", 'Minimum', 'Average', 'Maximum'], 
            ['Visual Irradiance (kWh/m2)', '{:.1f}'.format(self['omin']['virradm2{}'.format(frame)]), '{:.1f}'.format(self['oave']['virradm2{}'.format(frame)]), '{:.1f}'.format(self['omax']['virradm2{}'.format(frame)])]])
        self['tablefi{}'.format(frame)] = array([["", 'Minimum', 'Average', 'Maximum'], 
            ['Full Irradiance (kWh)', '{:.1f}'.format(self['omin']['firrad{}'.format(frame)]), '{:.1f}'.format(self['oave']['firrad{}'.format(frame)]), '{:.1f}'.format(self['omax']['firrad{}'.format(frame)])]])
        self['tablevi{}'.format(frame)] = array([["", 'Minimum', 'Average', 'Maximum'], 
            ['Visual Irradiance (kWh)', '{:.1f}'.format(self['omin']['virrad{}'.format(frame)]), '{:.1f}'.format(self['oave']['virrad{}'.format(frame)]), '{:.1f}'.format(self['omax']['virrad{}'.format(frame)])]])

        posis = [v.co for v in bm.verts if v[cindex] > 0] if self['cpoint'] == '1' else [f.calc_center_bounds() for f in bm.faces if f[cindex] > 1]
        reslists.append([str(frame), 'Zone', self.id_data.name, 'X', ' '.join([str(p[0]) for p in posis])])
        reslists.append([str(frame), 'Zone', self.id_data.name, 'Y', ' '.join([str(p[0]) for p in posis])])
        reslists.append([str(frame), 'Zone', self.id_data.name, 'Z', ' '.join([str(p[0]) for p in posis])])
        reslists.append([str(frame), 'Zone', self.id_data.name, 'Area', ' '.join([str(a) for a in areas])])
        reslists.append([str(frame), 'Zone', self.id_data.name, 'Full irradiance', ' '.join([str(g[firradres]) for g in geom if g[cindex] > 0])])
        reslists.append([str(frame), 'Zone', self.id_data.name, 'Visible irradiance', ' '.join([str(g[virradres]) for g in geom if g[cindex] > 0])])
        reslists.append([str(frame), 'Zone', self.id_data.name, 'Illuminance (Mlxh)', ' '.join([str(g[illures]) for g in geom if g[cindex] > 0])])
    bm.to_mesh(self.id_data.data)
    bm.free()
    return reslists
                    
def compcalcapply(self, scene, frames, rtcmds, simnode, curres, pfile):  
    svp = scene.vi_params
    pfs, epfs = [[] for f in frames], [[] for f in frames]
    self['compmat'] = [slot.material.name for slot in self.id_data.material_slots if slot.material.vi_params.mattype == '1'][0]
    self['omax'], self['omin'], self['oave'] = {}, {}, {}
    self['crit'], self['ecrit'], spacetype = retcrits(simnode, self['compmat'])    
    comps, ecomps =  {str(f): [] for f in frames}, {str(f): [] for f in frames}
    crits, dfpass, edfpass = [], {str(f): 0 for f in frames}, {str(f): 0 for f in frames} 
    selobj(bpy.context.view_layer, self.id_data)
    bm = bmesh.new()
    bm.from_mesh(self.id_data.data)
    clearlayers(bm, 'f')
    geom = bm.verts if simnode['goptions']['cp'] == '1' else bm.faces
    reslen = len(geom)
    cindex = geom.layers.int['cindex']
    pf = ('Fail', 'Pass')

    for f, frame in enumerate(frames):
        reslists, scores, escores, metric, emetric = [], [], [], [], []
        geom.layers.float.new('sv{}'.format(frame))
        geom.layers.float.new('df{}'.format(frame))
        geom.layers.float.new('res{}'.format(frame))
        dfres = geom.layers.float['df{}'.format(frame)]
        svres = geom.layers.float['sv{}'.format(frame)]
        res = geom.layers.float['res{}'.format(frame)]
        
        if geom.layers.string.get('rt{}'.format(frame)):
            rtframe = frame
        else:
            kints = [int(k[2:]) for k in geom.layers.string.keys()]
            rtframe  = max(kints) if frame > max(kints) else  min(kints)
        
        rt = geom.layers.string['rt{}'.format(rtframe)]
        
        for chunk in chunks([g for g in geom if g[rt]], int(svp['viparams']['nproc']) * 50):
            rtrun = Popen(rtcmds[f].split(), stdin = PIPE, stdout=PIPE, stderr=PIPE, universal_newlines=True).communicate(input = '\n'.join([c[rt].decode('utf-8') for c in chunk]))   
            xyzirrad = array([[float(v) for v in sl.split('\t')[:3]] for sl in rtrun[0].splitlines()])
            virrad = nsum(xyzirrad * array([0.26, 0.67, 0.065]), axis = 1)
            illu = virrad * 179
            df = illu * 0.01
            sv = self.retsv(scene, frame, rtframe, chunk, rt)

            for gi, gp in enumerate(chunk):
                gp[dfres] = df[gi].astype(float16)
                gp[svres] = sv[gi].astype(int8)
                gp[res] = illu[gi].astype(float32)
            
            curres += len(chunk)
            if pfile.check(curres) == 'CANCELLED':
                bm.free()
                return {'CANCELLED'}

        resillu = array([gp[res] for gp in geom if gp[cindex] > 0], dtype = float32)
        resdf = array([gp[dfres] for gp in geom if gp[cindex] > 0], dtype = float32)
        ressv = array([gp[svres] for gp in geom if gp[cindex] > 0], dtype = int8)
        
        self['omax']['df{}'.format(frame)] = nmax(resdf).astype(float64)
        self['omin']['df{}'.format(frame)] = nmin(resdf).astype(float64)
        self['oave']['df{}'.format(frame)] = nmean(resdf).astype(float64)
        self['omax']['sv{}'.format(frame)] =  1.0
        self['omin']['sv{}'.format(frame)] = 0.0
        self['oave']['sv{}'.format(frame)] = nmean(ressv)
        self['tabledf{}'.format(frame)] = array([["", 'Minimum', 'Average', 'Maximum'], 
            ['DF (%)', '{:.1f}'.format(self['omin']['df{}'.format(frame)]), '{:.1f}'.format(self['oave']['df{}'.format(frame)]), '{:.1f}'.format(self['omax']['df{}'.format(frame)])]])
        self['tablesv{}'.format(frame)] = array([['', '% area with', '% area without'], ['Sky View', '{:.1f}'.format(100 * nsum(ressv)/len(ressv)), '{:.1f}'.format(100 - 100 * nsum(ressv)/(len(ressv)))]])

        posis = [v.co for v in bm.verts if v[cindex] > 0] if self['cpoint'] == '1' else [f.calc_center_bounds() for f in bm.faces if f[cindex] > 1]
        reslists.append([str(frame), 'Zone', self.name, 'X', ' '.join([str(p[0]) for p in posis])])
        reslists.append([str(frame), 'Zone', self.name, 'Y', ' '.join([str(p[0]) for p in posis])])
        reslists.append([str(frame), 'Zone', self.name, 'Z', ' '.join([str(p[0]) for p in posis])])
        resdict = {'DF': resdf, 'Illuminance': resillu, 'Sky View': ressv}
        
        for unit in resdict:
            reslists.append([str(frame), 'Zone', self.name, unit, ' '.join([str(r) for r in resdict[unit]])])
        
        dftotarea, dfpassarea, edfpassarea, edftotarea = 0, 0, 0, 0
        oareas = self['lisenseareas'][str(frame)]
        oarea = sum(oareas)
        passarea = 0

        for c in self['crit']:
            if c[0] == 'Average':
                if c[2] == 'DF':
                    dfpass[str(frame)] = 1
                    dfpassarea = dfpassarea + oarea if sum(resdf)/reslen > float(c[3]) else dfpassarea
                    comps[str(frame)].append((0, 1)[sum(resdf)/reslen > float(c[3])])
                    comps[str(frame)].append(sum(resdf)/reslen)
                    dftotarea += oarea
                    metric.append(['Average DF', c[3], '{:.1f}'.format(comps[str(frame)][-1]), pf[comps[str(frame)][-2]]])
                    
            elif c[0] == 'Min':
                comps[str(frame)].append((0, 1)[min(resdf) > float(c[3])])
                comps[str(frame)].append(min(resdf))
                metric.append(['Minimum DF', c[3], '{:.1f}'.format(comps[str(frame)][-1]), pf[comps[str(frame)][-2]]])
    
            elif c[0] == 'Ratio':
                comps[str(frame)].append((0, 1)[min(resdf)/(sum(resdf)/reslen) >= float(c[3])])
                comps[str(frame)].append(min(resdf)/(sum(resdf)/reslen))
                metric.append(['Ratio of minimum to average DF', c[3], '{:.1f}'.format(comps[str(frame)][-1]), pf[comps[str(frame)][-2]]])
            
            elif c[0] == 'Percent':
                if c[2] == 'PDF':
                    dfpass[str(frame)] = 1
                    dfpassarea = sum([area for p, area in enumerate(oareas) if resdf[p] > int(c[3])])
                    comps[str(frame)].append((0, 1)[dfpassarea > float(c[1])*oarea/100])
                    comps[str(frame)].append(100*dfpassarea/oarea)
                    dftotarea += oarea
                    metric.append(['% area with Point DF > {}'.format(c[3]), c[1], '{:.1f}'.format(comps[str(frame)][-1]), pf[comps[str(frame)][-2]]])                    
                elif c[2] == 'Skyview': 
                    passarea = sum([area for p, area in enumerate(oareas) if ressv[p] > 0])
                    comps[str(frame)].append((0, 1)[passarea >= float(c[1])*oarea/100])
                    comps[str(frame)].append(100*passarea/oarea)
                    passarea = 0
                    metric.append(['% area with sky view', c[1], '{:.1f}'.format(comps[str(frame)][-1]), pf[comps[str(frame)][-2]]])
                elif c[2] == 'DF':  
                    passareapc = 100 * sum([area for p, area in enumerate(oareas) if resdf[p] > float(c[3])])/oarea
                    comps[str(frame)].append((0, 1)[sum([area * resdf[p] for p, area in enumerate(oareas)])/oarea > float(c[3])])
                    comps[str(frame)].append(sum([area * resdf[p] for p, area in enumerate(oareas)])/oarea)
                    metric.append(['% area with DF > {}'.format(c[3]), c[1], '{:.1f}'.format(passareapc), pf[passareapc >= float(c[1])]])
            scores.append(c[4])  

        passfails = [m[-1] for m in metric]

        if simnode['coptions']['canalysis'] == '0':
            if 'Pass' not in passfails:
                opf = 'FAIL'
            elif 'Fail' not in passfails:
                opf = 'PASS'
            elif 'Fail' in [c for i, c in enumerate(passfails) if scores[i] == '1']:
                opf = 'FAIL'
            elif 'Pass' not in [c for i, c in enumerate(passfails) if scores[i] == '0.75'] and len([c for i, c in enumerate(list(zip(metric))[-1]) if scores[i] == '0.75']) > 0:
                if 'Pass' not in [c for i, c in enumerate(passfails) if scores[i] == '0.5'] and len([c for i, c in enumerate(list(zip(metric))[-1]) if scores[i] == '0.5']) > 0:
                    opf = 'FAIL'
                else:
                    opf = 'PASS'
            else:
                opf = 'PASS'
            pfs[f].append(opf)
        
        elif simnode['coptions']['canalysis'] == '1':
            for met in metric:
                if self.data.materials[self['compmat']].crspacemenu == '0': # Kitchen Space
                
                    if met[0] == 'Average DF':
                        pfs[f].append((1, met[-1]))
                    else:
                        pfs[f].append((2, met[-1]))
                else: # Living Space
                    pfs[f].append((0, met[-1]))
                        
        elif simnode['coptions']['canalysis'] == '2':
            pfs[f] = [m[-1] for m in metric]

        if self['ecrit']:
            emetric = [['', '', '', ''], ['Exemplary requirements: ', '', '', '']]
        
            for e in self['ecrit']:
                if e[0] == 'Percent':
                    if e[2] == 'DF':
                        epassareapc = 100 * sum([area for p, area in enumerate(oareas) if resdf[p] > float(e[3])])/oarea
                        ecomps[str(frame)].append((0, 1)[sum([area * resdf[p] for p, area in enumerate(oareas)])/oarea > float(e[3])])
                        ecomps[str(frame)].append(sum([area * resdf[p] for p, area in enumerate(oareas)])/oarea)
                        emetric.append(['% area with DF > {}'.format(e[3]), e[1], '{:.1f}'.format(epassareapc), pf[epassareapc >= float(e[1])]])
                        
                    if e[2] == 'PDF':
                        edfpass[str(frame)] = 1
                        edfpassarea = sum([area for p, area in enumerate(oareas) if resdf[p] > float(e[3])])      
                        ecomps[str(frame)].append((0, 1)[dfpassarea > float(e[1])*oarea/100])
                        ecomps[str(frame)].append(100*edfpassarea/oarea)
                        edftotarea += oarea
                        emetric.append(['% area with Point DF > {}'.format(e[3]), e[1], '{:.1f}'.format(ecomps[str(frame)][-1]), pf[ecomps[str(frame)][-2]]])
        
                    elif e[2] == 'Skyview':
                        passarea = sum([area for p, area in enumerate(oareas) if ressv[p] > 0])
                        ecomps[str(frame)].append((0, 1)[passarea >= int(e[1]) * oarea/100])
                        ecomps[str(frame)].append(100*passarea/oarea)
                        passarea = 0
                        emetric.append(['% area with sky view', e[1], '{:.1f}'.format(ecomps[str(frame)][-1]), pf[ecomps[str(frame)][-2]]])
        
                elif e[0] == 'Min':
                    ecomps[str(frame)].append((0, 1)[min(resdf) > float(e[3])])
                    ecomps[str(frame)].append(min(resdf))
                    emetric.append(['Minimum DF', e[3], '{:.1f}'.format(ecomps[str(frame)][-1]), pf[ecomps[str(frame)][-2]]])
        
                elif e[0] == 'Ratio':
                    ecomps[str(frame)].append((0, 1)[min(resdf)/(sum(resdf)/reslen) >= float(e[3])])
                    ecomps[str(frame)].append(min(resdf)/(sum(resdf)/reslen))
                    emetric.append(['Ratio of minimum to average DF', e[3], '{:.1f}'.format(ecomps[str(frame)][-1]), pf[ecomps[str(frame)][-2]]])
        
                elif e[0] == 'Average':
                    ecomps[str(frame)].append((0, 1)[sum(resdf)/reslen > float(e[3])])
                    ecomps[str(frame)].append(sum(resdf)/reslen)
                    emetric.append(['% area with Average DF > {}'.format(e[3]), e[1], ecomps[str(frame)][-1], pf[ecomps[str(frame)][-2]]])
                crits.append(self['crit'])
                escores.append(e[4])
                
            epassfails = [em[-1] for em in emetric[2:]]

            if 'Pass' not in epassfails:
                epf = 'FAIL'     
            if 'Fail' not in epassfails:
                epf = 'PASS' 
            elif 'Fail' in [c for i, c in enumerate(epassfails) if escores[i] == '1']:
                epf = 'FAIL'
            elif 'Pass' not in [c for i, c in enumerate(epassfails) if escores[i] == '0.75'] and len([c for i, c in enumerate(list(zip(emetric))[-1]) if escores[i] == '0.75']) > 0:
                if 'Pass' not in [c for i, c in enumerate(epassfails) if escores[i] == '0.5'] and len([c for i, c in enumerate(list(zip(emetric))[-1]) if escores[i] == '0.5']) > 0:
                    epf = 'FAIL'
                else:
                    epf = 'EXEMPLARY'
            else:
                epf = 'EXEMPLARY'

            epfs[f].append(epf)
    
        if dfpass[str(frame)] == 1:
            dfpass[str(frame)] = 2 if dfpassarea/dftotarea >= (0.8, 0.35)[simnode['coptions']['canalysis'] == '0' and simnode['coptions']['buildtype'] == '4'] else dfpass[str(frame)]
        if edfpass[str(frame)] == 1:
            edfpass[str(frame)] = 2 if edfpassarea/edftotarea >= (0.8, 0.5)[simnode['coptions']['canalysis'] == '0' and simnode['coptions']['buildtype'] == '4'] else edfpass[str(frame)]
        
        smetric = [['Standard: {}'.format(('BREEAM HEA1', 'CfSH', 'Green Star', 'LEED EQ8.1')[int(simnode['coptions']['Type'])]), '', '', ''], 
                    ['Space type: {}'.format(spacetype), '', '', ''], ['', '', '', ''], ['Standard requirements:', 'Target', 'Result', 'Pass/Fail']] + metric
        self['tablecomp{}'.format(frame)] = smetric if not self['ecrit'] else smetric + emetric

    bm.to_mesh(self.id_data.data)
    bm.free()
    return (pfs, epfs, reslists)
    
def udidacalcapply(self, scene, frames, rccmds, simnode, curres, pfile):
    self['livires'] = {}
    self['compmat'] = [slot.material.name for slot in self.id_data.material_slots if slot.material.vi_params.mattype == '1'][0]
    selobj(scene, self)
    bm = bmesh.new()
    bm.from_mesh(self.id_data.data)
    bm.transform(self.matrix_world)
    clearlayers(bm, 'f')
    geom = bm.verts if self['cpoint'] == '1' else bm.faces
    reslen = len(geom)

    if self.get('wattres'):
        del self['wattres']
        
    illuarray = array((47.4, 120, 11.6)).astype(float32)
    vwattarray = array((0.265, 0.67, 0.065)).astype(float32)
    fwattarray = vwattarray * 1.64
    times = [datetime.datetime.strptime(time, "%d/%m/%y %H:%M:%S") for time in simnode['coptions']['times']]                          
    vecvals, vals = mtx2vals(open(simnode.inputs['Context in'].links[0].from_node['Options']['mtxfile'], 'r').readlines(), datetime.datetime(2010, 1, 1).weekday(), simnode, times)
    cbdm_days = [d for d in range(simnode['coptions']['sdoy'], simnode['coptions']['edoy'] + 1)] if scene['viparams']['visimcontext'] == 'LiVi CBDM' else [d for d in range(1, 366)]
    cbdm_hours = [h for h in range(simnode['coptions']['cbdm_sh'], simnode['coptions']['cbdm_eh'] + 1)]
    dno, hno = len(cbdm_days), len(cbdm_hours)    
    (luxmin, luxmax) = (simnode['coptions']['dalux'], simnode['coptions']['asemax']) if scene['viparams']['visimcontext'] != 'LiVi Compliance' else (300, 1000)
    vecvals = array([vv[2:] for vv in vecvals if vv[1] < simnode['coptions']['weekdays']]).astype(float32)
    hours = vecvals.shape[0]
    restypes = ('da', 'sda', 'ase', 'res', 'udilow', 'udisup', 'udiauto', 'udihi', 'kW', 'kW/m2', 'illu')
    self['livires']['cbdm_days'] = cbdm_days
    self['livires']['cbdm_hours'] = cbdm_hours

    for f, frame in enumerate(frames):        
        reslists = [[str(frame), 'Time', '', 'Month', ' '.join([str(t.month) for t in times])]]
        reslists.append([str(frame), 'Time', '', 'Day', ' '.join([str(t.day) for t in times])])
        reslists.append([str(frame), 'Time', '', 'Hour', ' '.join([str(t.hour) for t in times])])
        reslists.append([str(frame), 'Time', '', 'DOS', ' '.join([str(t.timetuple().tm_yday - times[0].timetuple().tm_yday) for t in times])])

        for restype in restypes:
            geom.layers.float.new('{}{}'.format(restype, frame))
        (resda, ressda, resase, res, resudilow, resudisup, resudiauto, resudihi, reskw, reskwm2, resillu) = [geom.layers.float['{}{}'.format(r, frame)] for r in restypes]

        if simnode['coptions']['buildtype'] == '1':
            geom.layers.float.new('sv{}'.format(frame))
            ressv = geom.layers.float['sv{}'.format(frame)]
        
        if geom.layers.string.get('rt{}'.format(frame)):
            rtframe = frame
        else:
            kints = [int(k[2:]) for k in geom.layers.string.keys()]
            rtframe  = max(kints) if frame > max(kints) else  min(kints)
        
        rt = geom.layers.string['rt{}'.format(rtframe)]
        totarea = sum([g.calc_area() for g in geom if g[rt]]) if self['cpoint'] == '0' else sum([vertarea(bm, g) for g in geom if g[rt]])
                
        for ch, chunk in enumerate(chunks([g for g in geom if g[rt]], int(scene['viparams']['nproc']) * 40)):
            sensrun = Popen(rccmds[f].split(), stdin=PIPE, stdout=PIPE, universal_newlines=True).communicate(input = '\n'.join([c[rt].decode('utf-8') for c in chunk]))
#            resarray = array([[float(v) for v in sl.split('\t') if v] for sl in sensrun[0].splitlines() if sl not in ('\n', '\r\n')]).reshape(len(chunk), 146, 3).astype(float32)
            resarray = array([[float(v) for v in sl.strip('\n').strip('\r\n').split('\t') if v] for sl in sensrun[0].splitlines()]).reshape(len(chunk), 146, 3).astype(float32)
            chareas = array([c.calc_area() for c in chunk]) if self['cpoint'] == '0' else array([vertarea(bm, c) for c in chunk]).astype(float32)
            sensarray = nsum(resarray*illuarray, axis = 2).astype(float32)
            wsensearray  = nsum(resarray*fwattarray, axis = 2).astype(float32)
            finalillu = inner(sensarray, vecvals).astype(float64)
            
            if scene['viparams']['visimcontext'] != 'LiVi Compliance':            
                finalwattm2 = inner(wsensearray, vecvals).astype(float32)
                wsensearraym2 = (wsensearray.T * chareas).T.astype(float32)
                finalwatt = inner(wsensearraym2, vecvals).astype(float32)  
                dabool = choose(finalillu >= simnode['coptions']['dalux'], [0, 1]).astype(int8)
                udilbool = choose(finalillu < simnode['coptions']['damin'], [0, 1]).astype(int8)
                udisbool = choose(finalillu < simnode['coptions']['dasupp'], [0, 1]).astype(int8) - udilbool
                udiabool = choose(finalillu < simnode['coptions']['daauto'], [0, 1]).astype(int8) - udilbool - udisbool
                udihbool = choose(finalillu >= simnode['coptions']['daauto'], [0, 1]).astype(int8)                       
                daareares = (dabool.T*chareas).T             
                udilareares = (udilbool.T*chareas).T
                udisareares = (udisbool.T*chareas).T
                udiaareares = (udiabool.T*chareas).T
                udihareares = (udihbool.T*chareas).T
                dares = dabool.sum(axis = 1)*100/hours
                udilow = udilbool.sum(axis = 1)*100/hours
                udisup = udisbool.sum(axis = 1)*100/hours
                udiauto = udiabool.sum(axis = 1)*100/hours
                udihi = udihbool.sum(axis = 1)*100/hours
                kwh = 0.001 * nsum(finalwatt, axis = 1)
                kwhm2 = 0.001 * nsum(finalwattm2, axis = 1)
            
            if scene['viparams']['visimcontext'] == 'LiVi Compliance' and simnode['coptions']['buildtype'] == '1':
                svres = self.retsv(scene, frame, rtframe, chunk, rt)
                sdaareas = where([sv > 0 for sv in svres], chareas, 0)
            else:
                sdaareas = chareas
            sdabool = choose(finalillu >= luxmin, [0, 1]).astype(int8)
            asebool = choose(finalillu >= luxmax, [0, 1]).astype(int8)
            aseareares = (asebool.T*chareas).T
            sdaareares = (sdabool.T*sdaareas).T            
            sdares = sdabool.sum(axis = 1)*100/hours
            aseres = asebool.sum(axis = 1)*1.0
                                    
            for gi, gp in enumerate(chunk):
                if scene['viparams']['visimcontext'] != 'LiVi Compliance':
                    gp[resda] = dares[gi]                
                    gp[res] = dares[gi]
                    gp[resudilow] = udilow[gi]
                    gp[resudisup] = udisup[gi]
                    gp[resudiauto] = udiauto[gi]
                    gp[resudihi] = udihi[gi]
                    gp[reskw] = kwh[gi]
                    gp[reskwm2] = kwhm2[gi]
                    gp[resillu] = max(finalillu[gi])
                
                elif simnode['coptions']['buildtype'] == '1':
                    gp[ressv] = svres[gi]
                gp[ressda] = sdares[gi]
                gp[resase] = aseres[gi]

            if not ch:
                if scene['viparams']['visimcontext'] != 'LiVi Compliance':
                    totfinalillu = finalillu
                    totdaarea = nsum(100 * daareares/totarea, axis = 0)
                    totudiaarea = nsum(100 * udiaareares/totarea, axis = 0)
                    totudisarea = nsum(100 * udisareares/totarea, axis = 0)
                    totudilarea = nsum(100 * udilareares/totarea, axis = 0)
                    totudiharea = nsum(100 * udihareares/totarea, axis = 0)                
                    
                if scene['viparams']['visimcontext'] == 'LiVi CBDM'  and simnode['coptions']['cbanalysis'] == '1':
                    totfinalwatt = nsum(finalwatt, axis = 0)#nsum(inner(sensarray, vecvals), axis = 0)
                    totfinalwattm2 = average(finalwattm2, axis = 0)
                else:
                    totsdaarea = nsum(sdaareares, axis = 0)
                    totasearea = nsum(aseareares, axis = 0)
            else:
                if scene['viparams']['visimcontext'] != 'LiVi Compliance':
                    nappend(totfinalillu, finalillu)
                    totdaarea += nsum(100 * daareares/totarea, axis = 0)
                    totudiaarea += nsum(100 * udiaareares/totarea, axis = 0)
                    totudilarea += nsum(100 * udilareares/totarea, axis = 0)
                    totudisarea += nsum(100 * udisareares/totarea, axis = 0)
                    totudiharea += nsum(100 * udihareares/totarea, axis = 0)
                if scene['viparams']['visimcontext'] == 'LiVi CBDM'  and simnode['coptions']['cbanalysis'] == '1':
                    totfinalwatt += nsum(finalwatt, axis = 0)#nsum(inner(sensarray, vecvals), axis = 0)
                    totfinalwattm2 += average(finalwattm2, axis = 0)
                else:
                    totsdaarea += nsum(sdaareares, axis = 0)
                    totasearea += nsum(aseareares, axis = 0)
              
            curres += len(chunk)
            if pfile.check(curres) == 'CANCELLED':
                bm.free()
                return {'CANCELLED'}

        if scene['viparams']['visimcontext'] != 'LiVi Compliance':
            dares = [gp[resda] for gp in geom] 
            udilow = [gp[resudilow] for gp in geom] 
            udisup = [gp[resudisup] for gp in geom]
            udiauto = [gp[resudiauto] for gp in geom]
            udihi = [gp[resudihi] for gp in geom]
            kwh = [gp[reskw] for gp in geom]
            kwhm2 = [gp[reskwm2] for gp in geom]
            self['omax']['udilow{}'.format(frame)] = max(udilow)
            self['omin']['udilow{}'.format(frame)] = min(udilow)
            self['oave']['udilow{}'.format(frame)] = sum(udilow)/reslen
            self['omax']['udisup{}'.format(frame)] = max(udisup)
            self['omin']['udisup{}'.format(frame)] = min(udisup)
            self['oave']['udisup{}'.format(frame)] = sum(udisup)/reslen
            self['omax']['udiauto{}'.format(frame)] = max(udiauto)
            self['omin']['udiauto{}'.format(frame)] = min(udiauto)
            self['oave']['udiauto{}'.format(frame)] = sum(udiauto)/reslen
            self['omax']['udihi{}'.format(frame)] = max(udihi)
            self['omin']['udihi{}'.format(frame)] = min(udihi)
            self['oave']['udihi{}'.format(frame)] = sum(udihi)/reslen
            self['omax']['da{}'.format(frame)] = max(dares)
            self['omin']['da{}'.format(frame)] = min(dares)
            self['oave']['da{}'.format(frame)] = sum(dares)/reslen
            self['omax']['illu{}'.format(frame)] = amax(totfinalillu)
            self['omin']['illu{}'.format(frame)] = amin(totfinalillu)
            self['oave']['illu{}'.format(frame)] = nmean(totfinalillu)/reslen
            self['livires']['dhilluave{}'.format(frame)] = average(totfinalillu, axis = 0).flatten().reshape(dno, hno).transpose().tolist()
            self['livires']['dhillumin{}'.format(frame)] = amin(totfinalillu, axis = 0).reshape(dno, hno).transpose().tolist()
            self['livires']['dhillumax{}'.format(frame)] = amax(totfinalillu, axis = 0).reshape(dno, hno).transpose().tolist()
            self['livires']['daarea{}'.format(frame)] = totdaarea.reshape(dno, hno).transpose().tolist()
            self['livires']['udiaarea{}'.format(frame)] = totudiaarea.reshape(dno, hno).transpose().tolist()
            self['livires']['udisarea{}'.format(frame)] = totudisarea.reshape(dno, hno).transpose().tolist()
            self['livires']['udilarea{}'.format(frame)] = totudilarea.reshape(dno, hno).transpose().tolist()
            self['livires']['udiharea{}'.format(frame)] = totudiharea.reshape(dno, hno).transpose().tolist()
            
            self['tableudil{}'.format(frame)] = array([["", 'Minimum', 'Average', 'Maximum'], 
                ['UDI-l (% area)', '{:.1f}'.format(self['omin']['udilow{}'.format(frame)]), '{:.1f}'.format(self['oave']['udilow{}'.format(frame)]), '{:.1f}'.format(self['omax']['udilow{}'.format(frame)])]])
            self['tableudis{}'.format(frame)] = array([["", 'Minimum', 'Average', 'Maximum'], 
                ['UDI-s (% area)', '{:.1f}'.format(self['omin']['udisup{}'.format(frame)]), '{:.1f}'.format(self['oave']['udisup{}'.format(frame)]), '{:.1f}'.format(self['omax']['udisup{}'.format(frame)])]])
            self['tableudia{}'.format(frame)] = array([["", 'Minimum', 'Average', 'Maximum'], 
                ['UDI-a (% area)', '{:.1f}'.format(self['omin']['udiauto{}'.format(frame)]), '{:.1f}'.format(self['oave']['udiauto{}'.format(frame)]), '{:.1f}'.format(self['omax']['udiauto{}'.format(frame)])]])
            self['tableudie{}'.format(frame)] = array([["", 'Minimum', 'Average', 'Maximum'], 
                ['UDI-e (% area)', '{:.1f}'.format(self['omin']['udihi{}'.format(frame)]), '{:.1f}'.format(self['oave']['udihi{}'.format(frame)]), '{:.1f}'.format(self['omax']['udihi{}'.format(frame)])]])
            self['tableillu{}'.format(frame)] = array([["", 'Minimum', 'Average', 'Maximum'], 
                ['Illuminance (lux)', '{:.1f}'.format(self['omin']['illu{}'.format(frame)]), '{:.1f}'.format(self['oave']['illu{}'.format(frame)]), '{:.1f}'.format(self['omax']['illu{}'.format(frame)])]])
            self['tableda{}'.format(frame)] = array([["", 'Minimum', 'Average', 'Maximum'], 
                ['Daylight availability (% time)', '{:.1f}'.format(self['omin']['da{}'.format(frame)]), '{:.1f}'.format(self['oave']['da{}'.format(frame)]), '{:.1f}'.format(self['omax']['da{}'.format(frame)])]])
            
            reslists.append([str(frame), 'Zone', self.name, 'Daylight Autonomy Area (%)', ' '.join([str(p) for p in totdaarea])])
            reslists.append([str(frame), 'Zone', self.name, 'UDI-a Area (%)', ' '.join([str(p) for p in totudiaarea])])
            reslists.append([str(frame), 'Zone', self.name, 'UDI-s Area (%)', ' '.join([str(p) for p in totudisarea])])
            reslists.append([str(frame), 'Zone', self.name, 'UDI-l Area (%)', ' '.join([str(p) for p in totudilarea])])
            reslists.append([str(frame), 'Zone', self.name, 'UDI-h Area (%)', ' '.join([str(p) for p in totudiharea])])
        
        if scene['viparams']['visimcontext'] == 'LiVi CBDM' and simnode['coptions']['cbanalysis'] == '1': 
            self['omax']['kW{}'.format(frame)] = max(kwh)
            self['omin']['kW{}'.format(frame)] = min(kwh)
            self['oave']['kW{}'.format(frame)] = sum(kwh)/reslen
            self['omax']['kW/m2{}'.format(frame)] = max(kwhm2)
            self['omin']['kW/m2{}'.format(frame)] = min(kwhm2)
            self['oave']['kW/m2{}'.format(frame)] = sum(kwhm2)/reslen
            self['livires']['kW{}'.format(frame)] =  (0.001*totfinalwatt).reshape(dno, hno).transpose().tolist()
            self['livires']['kW/m2{}'.format(frame)] =  (0.001*totfinalwattm2).reshape(dno, hno).transpose().tolist()
            self['tablekwh{}'.format(frame)] = array([["", 'Minimum', 'Average', 'Maximum'], 
                ['Irradiance (kW)', '{:.1f}'.format(self['omin']['kW{}'.format(frame)]), '{:.1f}'.format(self['oave']['kW{}'.format(frame)]), '{:.1f}'.format(self['omax']['kW{}'.format(frame)])]])
            self['tablekwhm2{}'.format(frame)] = array([["", 'Minimum', 'Average', 'Maximum'], 
                ['Irradiance (kW/m2)', '{:.1f}'.format(self['omin']['kW/m2{}'.format(frame)]), '{:.1f}'.format(self['oave']['kW/m2{}'.format(frame)]), '{:.1f}'.format(self['omax']['kW/m2{}'.format(frame)])]])
            reslists.append([str(frame), 'Zone', self.name, 'kW', ' '.join([str(p) for p in 0.001 * totfinalwatt])])
            reslists.append([str(frame), 'Zone', self.name, 'kW/m2', ' '.join([str(p) for p in 0.001 * totfinalwattm2])])
        else:
            sdares = [gp[ressda] for gp in geom]
            aseres = [gp[resase] for gp in geom]
            if scene['viparams']['visimcontext'] == 'LiVi Compliance' and simnode['coptions']['buildtype'] == '1':
                overallsdaarea = sum([g.calc_area() for g in geom if g[rt] and g[ressv]]) if self['cpoint'] == '0' else sum([vertarea(bm, g) for g in geom if g[rt] and g[ressv]]) 
            else:
                overallsdaarea = totarea
            self['omax']['sda{}'.format(frame)] = max(sdares)
            self['omin']['sda{}'.format(frame)] = min(sdares)
            self['oave']['sda{}'.format(frame)] = sum(sdares)/reslen
            self['omax']['ase{}'.format(frame)] = max(aseres)
            self['omin']['ase{}'.format(frame)] = min(aseres)
            self['oave']['ase{}'.format(frame)] = sum(aseres)/reslen
            self['livires']['asearea{}'.format(frame)] = (100 * totasearea/totarea).reshape(dno, hno).transpose().tolist()
            self['livires']['sdaarea{}'.format(frame)] = (100 * totsdaarea/overallsdaarea).reshape(dno, hno).transpose().tolist()
            self['tablesda{}'.format(frame)] = array([["", 'Minimum', 'Average', 'Maximum'], 
                ['sDA (% hours)', '{:.1f}'.format(self['omin']['sda{}'.format(frame)]), '{:.1f}'.format(self['oave']['sda{}'.format(frame)]), '{:.1f}'.format(self['omax']['sda{}'.format(frame)])]])
            self['tablease{}'.format(frame)] = array([["", 'Minimum', 'Average', 'Maximum'], 
                ['ASE (hrs)', '{:.1f}'.format(self['omin']['ase{}'.format(frame)]), '{:.1f}'.format(self['oave']['ase{}'.format(frame)]), '{:.1f}'.format(self['omax']['ase{}'.format(frame)])]])
            reslists.append([str(frame), 'Zone', self.name, 'Annual Sunlight Exposure (% area)', ' '.join([str(p) for p in 100 * totasearea/totarea])])
            reslists.append([str(frame), 'Zone', self.name, 'Spatial Daylight Autonomy (% area)', ' '.join([str(p) for p in 100 * totsdaarea/overallsdaarea])])
            
        metric, scores, pf = [], [], ('Fail', 'Pass')

        if scene['viparams']['visimcontext'] == 'LiVi Compliance':
            self['crit'], self['ecrit'], spacetype = retcrits(simnode, self['compmat'])
            sdapassarea, asepassarea, comps = 0, 0, {str(f): [] for f in frames}
            oareas = self['lisenseareas'][str(frame)]
            oarea = sum(oareas)
            geom.ensure_lookup_table()
            hoarea = sum([oa for o, oa in enumerate(oareas) if geom[o][ressv] > 0]) if simnode['coptions']['buildtype'] == '3' else oarea
            aoarea = hoarea if simnode['coptions']['buildtype'] == '1' else oarea     
            self['oarea'] = aoarea

            for c in self['crit']:
                if c[0] == 'Percent':        
                    if c[2] == 'SDA':
                        sdapassarea = sum([area for p, area in enumerate(oareas) if sdares[p] >= 50 and svres[p] > 0]) if simnode['coptions']['buildtype'] == '1' else sum([area for p, area in enumerate(oareas) if sdares[p] >= 50])
                        comps[str(frame)].append((0, 1)[sdapassarea >= float(c[1])*oarea/100])
                        comps[str(frame)].append(100*sdapassarea/aoarea)
                        self['sdapassarea'] = sdapassarea
                        metric.append(['% area with SDA', c[1], '{:.1f}'.format(comps[str(frame)][-1]), pf[comps[str(frame)][-2]]])
                    
                    elif c[2] == 'ASE':
                        asepassarea = sum([area for p, area in enumerate(oareas) if aseres[p] > 250 and svres[p] > 0]) if simnode['coptions']['buildtype'] == '1' else sum([area for p, area in enumerate(oareas) if aseres[p] > 250])
                        comps[str(frame)].append((0, 1)[asepassarea <= float(c[1])*aoarea/100])
                        comps[str(frame)].append(100*asepassarea/aoarea)
                        self['asepassarea'] = asepassarea
                        metric.append(['% area with ASE', c[1], '{:.1f}'.format(comps[str(frame)][-1]), pf[comps[str(frame)][-2]]])
                    scores.append(c[4])

            self['comps'] = comps
            self['tablecomp{}'.format(frame)] = [['Standard: {}'.format('LEEDv4 EQ8.1'), '', '', ''], 
                ['Space type: {}'.format(spacetype), '', '', ''], ['', '', '', ''], ['Standard requirements:', 'Target', 'Result', 'Pass/Fail']] + metric
        
    bm.transform(self.matrix_world.inverted())        
    bm.to_mesh(self.data)
    bm.free()
    return [m[-1] for m in metric], scores, reslists

def retcrits(simnode, matname):
    ecrit = []
    mat = bpy.data.materials[matname]
    if simnode['coptions']['canalysis'] == '0':
        if simnode['coptions']['buildtype'] in ('0', '5'):
            if not mat.vi_params.gl_roof:
                crit = [['Percent', 80, 'DF', 2, '1'], ['Ratio', 100, 'Uni', 0.4, '0.5'], ['Min', 100, 'PDF', 0.8, '0.5'], ['Percent', 80, 'Skyview', 1, '0.75']]
                ecrit = [['Percent', 80, 'DF', 4, '1'], ['Min', 100, 'PDF', 1.6, '0.75']] if simnode['coptions']['storey'] == '0' else [['Percent', 80, 'DF', 3, '1'], ['Min', 100, 'PDF', 1.2, '0.75']] 
            else:
                crit = [['Percent', 80, 'DF', 2, '1'], ['Ratio', 100, 'Uni', 0.7, '0.5'], ['Min', 100, 'PDF', 1.4, '0.5'], ['Percent', 100, 'Skyview', 1, '0.75']]
                ecrit = [['Percent', 80, 'DF', 4, '1'], ['Min', 100, 'PDF', 2.8, '0.75']] if simnode['coptions']['storey'] == '0' else [['Percent', 80, 'DF', 3, '1'], ['Min', 100, 'PDF', 2.1, '0.75']]
            spacetype = 'School' if simnode['coptions']['buildtype'] == '0' else 'Office & Other'
        elif simnode['coptions']['buildtype'] == '1':
            if not mat.vi_params.gl_roof:
                crit = [['Percent', 80, 'DF', 2, '1'], ['Ratio', 100, 'Uni', 0.4, '0.5'], ['Min', 100, 'PDF', 0.8, '0.5'], ['Percent', 80, 'Skyview', 1, '0.75']]
                ecrit = [['Percent', 80, 'DF', 4, '1'], ['Min', 100, 'PDF', 1.6, '0.75']] if simnode['coptions']['storey'] == '0' else [['Percent', 80, 'DF', 3, '1'], ['Min', 100, 'PDF', 1.2, '0.75']]
            else:
                crit = [['Percent', 80, 'DF', 2, '1'], ['Ratio', 100, 'Uni', 0.7, '0.5'], ['Min', 100, 'PDF', 1.4, '0.5'], ['Percent', 100, 'Skyview', 1, '0.75']]
                ecrit= [['Percent', 80, 'DF', 4, '1'], ['Min', 100, 'PDF', 2.8, '0.75']] if simnode['coptions']['storey'] == '0' else [['Percent', 80, 'DF', 3, '1'], ['Min', 100, 'PDF', 2.1, '0.75']]
            spacetype = 'Higher Education'
        elif simnode['coptions']['buildtype'] == '2':
            crit = [['Percent', 80, 'DF', 2, '1']] if mat.hspacemenu == '0' else [['Percent', 80, 'DF', 3, '2']]
            ecrit = [['Percent', 80, 'DF', 4, '1'], ['Min', 100, 'PDF', 1.6, '0.75']] if simnode['coptions']['storey'] == '0' else [['Min', 100, 'PDF', 1.6, '0.75'], ['Min', 100, 'PDF', 1.2, '0.75']]
            spacetype = 'Healthcare - Patient' if mat.hspacemenu == '0' else 'Healthcare - Public'
        elif simnode['coptions']['buildtype'] == '3':
            if mat.brspacemenu == '0':
                crit = [['Percent', 80, 'DF', 2, '1'], ['Percent', 100, 'Skyview', 1, '0.75']]
                ecrit = [['Percent', 80, 'DF', 4, '1'], ['Min', 100, 'PDF', 1.6, '0.75']] if simnode['coptions']['storey'] == '0' else [['Percent', 80, 'DF', 3, '1'], ['Min', 100, 'PDF', 1.2, '0.75']]
                spacetype = 'Residential - Kitchen'
            elif mat.brspacemenu == '1':
                crit = [['Percent', 80, 'DF', 1.5, '1'], ['Percent', 100, 'Skyview', 1, '0.75']]
                ecrit = [['Percent', 80, 'DF', 4, '1'], ['Min', 100, 'PDF', 1.6, '0.75']] if simnode['coptions']['storey'] == '0' else [['Percent', 80, 'DF', 3, '1'], ['Min', 100, 'PDF', 1.2, '0.75']]
                spacetype = 'Residential - Living/Dining/Study'
            elif mat.brspacemenu == '2':
                if not mat.gl_roof:
                    crit = [['Percent', 80, 'DF', 2, '1'], ['Ratio', 100, 'Uni', 0.4, '0.5'], ['Min', 100, 'PDF', 0.8, '0.5'], ['Percent', 80, 'Skyview', 1, '0.75']]
                    ecrit = [['Percent', 80, 'DF', 4, '1'], ['Min', 100, 'PDF', 1.6, '0.75']] if simnode['coptions']['storey'] == '0' else [['Percent', 80, 'DF', 3, '1'], ['Min', 100, 'PDF', 1.2, '0.75']]
                else:
                    crit = [['Percent', 80, 'DF', 2, '1'], ['Ratio', 100, 'Uni', 0.7, '0.5'],['Min', 100, 'PDF', 1.4, '0.5'], ['Percent', 100, 'Skyview', 1, '0.75']] 
                    ecrit = [['Percent', 80, 'DF', 4, '1'], ['Min', 100, 'PDF', 2.8, '0.75']] if simnode['coptions']['storey'] == '0' else [['Percent', 80, 'DF', 3, '1'], ['Min', 100, 'PDF', 2.1, '0.75']]
                spacetype = 'Residential - Communal'

        elif simnode['coptions']['buildtype'] == '4':
            if mat.respacemenu == '0':
                crit = [['Percent', 35, 'PDF', 2, '1']]
                ecrit = [['Percent', 50, 'PDF', 2, '1']]
                spacetype = 'Retail - Sales'
            elif mat.respacemenu == '1':
                if not mat.gl_roof:
                    crit = [['Percent', 80, 'DF', 2, '1'], ['Ratio', 100, 'Uni', 0.4, '0.5'], ['Min', 100, 'PDF', 0.8, '0.5'], ['Percent', 80, 'Skyview', 1, '0.75']] 
                    ecrit = [['Percent', 80, 'DF', 4, '1'], ['Min', 100, 'PDF', 1.6, '0.75']] if simnode['coptions']['storey'] == '0' else [['Percent', 80, 'DF', 3, '1'], ['Min', 100, 'PDF', 1.2, '0.75']]   
                else:
                    crit = [['Percent', 80, 'DF', 2, '1'], ['Ratio', 100, 'Uni', 0.7, '0.5'], ['Min', 100, 'PDF', 1.4, '0.5'], ['Percent', 100, 'Skyview', 1, '0.75']]
                    ecrit = [['Percent', 80, 'DF', 4, '1'], ['Min', 100, 'PDF', 2.8, '0.75']] if simnode['coptions']['storey'] == '0' else [['Percent', 80, 'DF', 3, '1'],['Min', 100, 'PDF', 2.1, '0.75']] 
                spacetype = 'Retail - Occupied'
    
    elif simnode['coptions']['canalysis'] == '1':
        crit = [['Average', 100, 'DF', 2, '1'], ['Percent', 80, 'Skyview', 1, '0.75']] if mat.crspacemenu == '0' else [['Average', 100, 'DF', 1.5, '1'], ['Percent', 80, 'Skyview', 1, '0.75']]
        spacetype = 'Residential - Kitchen' if mat.crspacemenu == '0' else 'Residential - Living/Dining/Study'

    elif simnode['coptions']['canalysis'] == '2':
        if simnode['coptions']['buildtype'] == '0':
            crit = [['Percent', 30, 'DF', 2, '1'], ['Percent', 60, 'DF', 2, '1'], ['Percent', 90, 'DF', 2, '1'], ['Percent', 50, 'DF', 4, '1']]
            spacetype = 'School'
        if simnode['coptions']['buildtype'] == '1':
            crit = [['Percent', 30, 'DF', 2, '1'], ['Percent', 60, 'DF', 2, '1'], ['Percent', 90, 'DF', 2, '1']]
            spacetype = 'Higher Education'
        if simnode['coptions']['buildtype'] == '2':
            crit = [['Percent', 30, 'DF', 2.5, '1'], ['Percent', 60, 'DF', 2.5, '1'], ['Percent', 90, 'DF', 2.5, '1']] if mat.hspacemenu == '0' else [['Percent', 30, 'DF', 3, '1'], ['Percent', 60, 'DF', 3, '1'], ['Percent', 90, 'DF', 3, '1']]
            spacetype = 'Healthcare'
        if simnode['coptions']['buildtype'] == '3':
            crit = [['Percent', 60, 'DF', 2.0, '1'], ['Percent', 90, 'DF', 2.0, '1']] if mat.brspacemenu == '0' else [['Percent', 60, 'DF', 1.5, '1'], ['Percent', 90, 'DF', 1.5, '1']]   
            spacetype = 'Residential'             
        if simnode['coptions']['buildtype'] in ('4', '5'):
            crit = [['Percent', 30, 'DF', 2, '1'], ['Percent', 60, 'DF', 2, '1'], ['Percent', 90, 'DF', 2, '1']]
            spacetype = 'Retail/Office/Public' 
            
    elif simnode['coptions']['canalysis'] == '3':
        spacetype = ('Office/Education/Commercial', 'Healthcare')[int(simnode['coptions']['buildtype'])]
        if simnode['coptions']['buildtype'] == '0':
            crit = [['Percent', 55, 'SDA', 300, '2'], ['Percent', 75, 'SDA', 300, '1']]
        else:
            crit = [['Percent', 75, 'SDA', 300, '1'], ['Percent', 90, 'SDA', 300, '1']]
#            spacetype = 'Healthcare'
        crit.append(['Percent', 10, 'ASE', 1000, '1', 250])
                   
    return [[c[0], str(c[1]), c[2], str(c[3]), c[4]] for c in crit[:]], [[c[0], str(c[1]), c[2], str(c[3]), c[4]] for c in ecrit[:]], spacetype