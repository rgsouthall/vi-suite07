import bpy, blf, mathutils, datetime
from bpy_extras import view3d_utils
from .vi_func import ret_vp_loc, viewdesc, draw_index, draw_time, blf_props
from math import pi

def spnumdisplay(disp_op, context):
    scene = context.scene

    if bpy.data.objects.get('SPathMesh'):
        spob = bpy.data.objects['SPathMesh'] 
        ob_mat = spob.matrix_world
        mid_x, mid_y, width, height = viewdesc(context)
        vl = ret_vp_loc(context)
        blf_props(scene, width, height)
        
        if scene.vi_params.sp_hd:
            pvecs = [ob_mat@mathutils.Vector(p[:]) for p in spob['numpos'].values()]
            pvals = [int(p.split('-')[1]) for p in spob['numpos'].keys()]
            p2ds = [view3d_utils.location_3d_to_region_2d(context.region, context.region_data, p) for p in pvecs]
            vispoints = [pi for pi, p in enumerate(pvals) if p2ds[pi] and 0 < p2ds[pi][0] < width and 0 < p2ds[pi][1] < height and scene.ray_cast(context.view_layer, vl, pvecs[pi] - vl, distance = (pvecs[pi] - vl).length)[4] == spob]
            
            if vispoints:
                hs = [pvals[pi] for pi in vispoints]
                posis = [p2ds[pi] for pi in vispoints]                
                draw_index(posis, hs, scene.vi_params.display_rp_fs, scene.vi_params.display_rp_fc, scene.vi_params.display_rp_fsh)
                
        if [ob.get('VIType') == 'Sun' for ob in bpy.data.objects] and scene['spparams']['suns'] == '0':
            sobs = [ob for ob in bpy.data.objects if ob.get('VIType') == 'Sun']
            
            if sobs and scene.vi_params.sp_td:
                sunloc = ob_mat@sobs[0].location
                solpos = view3d_utils.location_3d_to_region_2d(context.region, context.region_data, sunloc)
                
                try:
                    if 0 < solpos[0] < width and 0 < solpos[1] < height and not scene.ray_cast(context.view_layer, sobs[0].location + 0.05 * (vl - sunloc), vl - sunloc)[0]:
                        soltime = datetime.datetime.fromordinal(scene.solday)
                        soltime += datetime.timedelta(hours = scene.solhour)
                        sre = sobs[0].rotation_euler
                        blf_props(scene, width, height)
                        draw_time(solpos, soltime.strftime('  %d %b %X') + ' alt: {:.1f} azi: {:.1f}'.format(90 - sre[0]*180/pi, (180, -180)[sre[2] < -pi] - sre[2]*180/pi), 
                                   scene.vi_params.display_rp_fs, scene.vi_params.display_rp_fc, scene.vi_params.display_rp_fsh)
                        
                except Exception as e:
                    print(e)
        blf.disable(0, 4)
    else:
        return