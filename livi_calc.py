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

import bpy, os, datetime, shlex
from subprocess import Popen, PIPE
from time import sleep
from . import livi_export
from .vi_func import selobj, progressbar, progressfile, logentry
from .livi_func import retpmap

def radfexport(scene, export_op, connode, geonode, frames):
    for frame in frames:
        livi_export.fexport(scene, frame, export_op, connode, geonode, pause = 1)

def li_calc(calc_op, simnode, simacc, **kwargs): 
    scene = bpy.context.scene
    svp = scene.vi_params
    vl = bpy.context.view_layer
    pfs, epfs, curres = [], [], 0
    context = simnode['coptions']['Context']
    subcontext = simnode['coptions']['Type']
    patches = simnode['coptions']['cbdm_res']
    svp['liparams']['maxres'], svp['liparams']['minres'], svp['liparams']['avres'] = {}, {}, {}
    frames = range(svp['liparams']['fs'], svp['liparams']['fe'] + 1) if not kwargs.get('genframe') else [kwargs['genframe']]
    os.chdir(svp['viparams']['newdir'])
    rtcmds, rccmds = [], []
    builddict = {'0': ('School', 'Higher Education', 'Healthcare', 'Residential', 'Retail', 'Office & Other'), '2': ('School', 'Higher Education', 'Healthcare', 'Residential', 'Retail', 'Office & Other'), '3': ('Office/Education/Commercial', 'Healthcare')}
    
    for f, frame in enumerate(frames):
        if simnode.pmap:
            pmappfile = open(os.path.join(svp['viparams']['newdir'], 'viprogress'), 'w')
            pmappfile.close()
            pfile = progressfile(svp['viparams']['newdir'], datetime.datetime.now(), 100)
            kivyrun = progressbar(os.path.join(svp['viparams']['newdir'], 'viprogress'), 'Photon map')
            errdict = {'fatal - too many prepasses, no global photons stored\n': "Too many prepasses have occurred. Make sure light sources can see your geometry",
            'fatal - too many prepasses, no global photons stored, no caustic photons stored\n': "Too many prepasses have occurred. Turn off caustic photons and encompass the scene",
           'fatal - zero flux from light sources\n': "No light flux, make sure there is a light source and that photon port normals point inwards",
           'fatal - no light sources\n': "No light sources. Photon mapping does not work with HDR skies",
           'fatal - no valid photon ports found\n': 'Re-export the geometry'}
            amentry, pportentry, cpentry, cpfileentry = retpmap(simnode, frame, scene)
            open('{}.pmapmon'.format(svp['viparams']['filebase']), 'w')
            
            if context == 'Basic' or (context == 'CBDM' and subcontext == '0'):
                pmcmd = 'mkpmap -n {6} -t 10 -e "{1}.pmapmon" -fo+ -bv+ -apD 0.001 {0} -apg "{1}-{2}.gpm" {3} {4} {5} "{1}-{2}.oct"'.format(pportentry, svp['viparams']['filebase'], frame, simnode.pmapgno, cpentry, amentry, svp['viparams']['wnproc'])
            else:
                pmcmd = 'mkpmap -n {3} -t 10 -e "{1}.pmapmon" -fo+ -bv+ -apC "{1}.cpm" {0} "{1}-{2}.oct"'.format(simnode.pmapgno, svp['viparams']['filebase'], frame, svp['viparams']['wnproc'])
            
            logentry('Generating photon map: {}'.format(pmcmd))
            pmrun = Popen(shlex.split(pmcmd), stderr = PIPE, stdout = PIPE)
                
            while pmrun.poll() is None:   
                sleep(10)
                with open('{}.pmapmon'.format(svp['viparams']['filebase']), 'r') as vip:
                    for line in vip.readlines()[::-1]:
                        if '%' in line:
                            curres = float(line.split()[6][:-2])/len(frames)
                            break
                                
                if pfile.check(curres) == 'CANCELLED': 
                    pmrun.kill()                                   
                    return 'CANCELLED'
            
            if kivyrun.poll() is None:
                kivyrun.kill()
                    
            with open('{}.pmapmon'.format(svp['viparams']['filebase']), 'r') as pmapfile:
                pmlines = pmapfile.readlines()
                if pmlines:
                    for line in pmlines:
                        if line in errdict:
                            calc_op.report({'ERROR'}, errdict[line])
                            return 'CANCELLED'
                        if 'fatal - ' in line:
                            calc_op.report({'ERROR'}, line)
                            return 'CANCELLED'
                else:
                    calc_op.report({'ERROR'}, 'There is a problem with pmap generation. Check there are no non-ascii characters in the project directory file path')
                    return 'CANCELLED'
                
        if context == 'Basic' or (context == 'CBDM' and subcontext == '0'):# or (context == 'Compliance' and int(subcontext) < 3):
            if os.path.isfile("{}-{}.af".format(svp['viparams']['filebase'], frame)):
                os.remove("{}-{}.af".format(svp['viparams']['filebase'], frame))
            if simnode.pmap:
                rtcmds.append('rtrace -n {0} -w {1} -ap "{2}-{3}.gpm" 50 {4} -faa -h -ov -I "{2}-{3}.oct"'.format(svp['viparams']['nproc'], simnode['radparams'], svp['viparams']['filebase'], frame, cpfileentry)) #+" | tee "+lexport.newdir+lexport.fold+self.simlistn[int(lexport.metric)]+"-"+str(frame)+".res"
            else:
                rtcmds.append('rtrace -n {0} -w {1} -faa -h -ov -I "{2}-{3}.oct"'.format(svp['viparams']['nproc'], simnode['radparams'], svp['viparams']['filebase'], frame)) #+" | tee "+lexport.newdir+lexport.fold+self.simlistn[int(lexport.metric)]+"-"+str(frame)+".res"
        else:
            if simnode.pmap:
                rccmds.append('rcontrib -w  -h -I -fo -ap {2}.cpm -bn {4} {0} -n {1} -f tregenza.cal -b tbin -m sky_glow "{2}-{3}.oct"'.format(simnode['radparams'], svp['viparams']['nproc'], svp['viparams']['filebase'], frame, patches))
            else:   
                rccmds.append('rcontrib -w  -h -I -fo -bn {} {} -n {} -f tregenza.cal -b tbin -m sky_glow "{}-{}.oct"'.format(patches, simnode['radparams'], svp['viparams']['nproc'], svp['viparams']['filebase'], frame))

    try:
        tpoints = [o.vi_params['rtpnum'] for o in bpy.data.objects if o.name in svp['liparams']['livic']]
    except:
        calc_op.report({'ERROR'}, 'Re-export the LiVi geometry')
        return 'CANCELLED'

    calcsteps = sum(tpoints) * len(frames)
    pfile = progressfile(svp['viparams']['newdir'], datetime.datetime.now(), calcsteps)
    kivyrun = progressbar(os.path.join(svp['viparams']['newdir'], 'viprogress'), 'Lighting')
    reslists = []
    obs = [o for o in bpy.data.objects if o.name in svp['liparams']['livic']]

    for oi, o in enumerate(obs):
        ovp = o.vi_params
        curres = sum(tpoints[:oi] * len(frames))
        selobj(vl, o)
        ovp['omax'], ovp['omin'], ovp['oave']  = {}, {}, {}
        
        if context == 'Basic':
            bccout = ovp.basiccalcapply(scene, frames, rtcmds, simnode, curres, pfile)
            if bccout == 'CANCELLED':
                if kivyrun.poll() is None:
                    kivyrun.kill()
                return 'CANCELLED'
            else:
                reslists += bccout
                
        elif context == 'CBDM' and subcontext == '0':
            lhout = ovp.lhcalcapply(scene, frames, rtcmds, simnode, curres, pfile)
            if lhout  == 'CANCELLED':
                if kivyrun.poll() is None:
                    kivyrun.kill()
                return 'CANCELLED'
            else:
                reslists += lhout
        
        elif (context == 'CBDM' and subcontext in ('1', '2')):# or (context == 'Compliance' and subcontext == '3'):
            cbdmout = ovp.udidacalcapply(scene, frames, rccmds, simnode, curres, pfile)
            if cbdmout == 'CANCELLED':
                if kivyrun.poll() is None:
                    kivyrun.kill()
                return 'CANCELLED'
            else:
                reslists += cbdmout
           
    if kivyrun.poll() is None:
        kivyrun.kill()
        
    return reslists
            

