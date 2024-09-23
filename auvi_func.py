import bpy

def setsceneauvivals(scene):
    svp = scene.vi_params
    svp['auparams']['maxres'], svp['auparams']['minres'], svp['auparams']['avres'] = {}, {}, {}
    res = svp.au_disp_menu
    olist = [o for o in bpy.data.objects if o.vi_params.vi_type_string == 'AuVi Calc']

    for frame in range(svp['auparams']['fs'], svp['auparams']['fe'] + 1):
        svp['auparams']['maxres'][str(frame)] = max([o.vi_params['omax']['{}{}'.format(res, frame)] for o in olist])
        svp['auparams']['minres'][str(frame)] = min([o.vi_params['omin']['{}{}'.format(res, frame)] for o in olist])
        svp['auparams']['avres'][str(frame)] = sum([o.vi_params['oave']['{}{}'.format(res, frame)] for o in olist])/len([o.vi_params['oave']['{}{}'.format(res, frame)] for o in olist])

    svp.vi_leg_max = max(svp['auparams']['maxres'].values())
    svp.vi_leg_min = min(svp['auparams']['minres'].values())