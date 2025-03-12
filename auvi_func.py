import bpy, math
import numpy as np
try:
    import pyroomacoustics as pra
except:
    pass

Qf = 1  # source directivity factor
c = 343.0  # sound speed (m/s)
Lsf = [62.9, 62.9, 59.2, 53.2, 47.2, 41.2, 35.2]  # UNI EN ISO 9921:2014 Table H.2 (speech spectrum at 1 m in front of the mouth of a male speaker (dB))
Lnf = [41.0, 43.0, 50.0, 47.0, 42.0, 42.0, 39.0]  # UNI EN ISO 9921:2014 Table H.1 (ambient noise spectrum (dB))
Lsf = np.array(Lsf)
Lnf = np.array(Lnf)
F = [0.63, 0.8, 1, 1.25, 1.6, 2, 2.5, 3.15, 4, 5, 6.3, 8, 10, 12.5]  # BS EN 60268-16:2011 par.A.2.2 (STI modulation frequencies (Hz))
alpha_male = [0.085, 0.127, 0.230, 0.233, 0.309, 0.224, 0.173]  # BS EN 60268-16:2011 Table A.3 (MTI octave band weighting factors)
beta_male = [0.085, 0.078, 0.065, 0.011, 0.047, 0.095, 0]  # BS EN 60268-16:2011 Table A.3 (MTI octave band weighting factors)
alpha_female = [0, 0.117, 0.223, 0.216, 0.328, 0.250, 0.194]  # BS EN 60268-16:2011 Table A.3 (MTI octave band weighting factors)
beta_female = [0, 0.099, 0.066, 0.062, 0.025, 0.076, 0]  # BS EN 60268-16:2011 Table A.3 (MTI octave band weighting factors)

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

def modulation_depth_reduction_factor(r, Qf, rcf, F, Tf, Lsf, Lnf):
    A = (Qf / r**2) + (1 / rcf**2) * (1 + ((2 * math.pi * F * Tf) / 13.8)**2)**-1
    # print("A", A)
    B = ((2 * math.pi * F * Tf) / (13.8 * rcf**2)) * (1 + (2 * math.pi * F * Tf / 13.8)**2)**-1
    # print("B", B)
    C = (Qf / r**2) + (1 / rcf**2) + Qf * 10**((-Lsf + Lnf) / 10)
    # print("C", C)
    mfF = (math.sqrt(A**2 + B**2)) / C
    return mfF


def Speech_Transmssion_Index(MTI, alpha, beta):
    MTI_alpha_pond = np.multiply(MTI, alpha)
    MTI_beta_pond = np.multiply(MTI, beta)
    sum_MTI_alpha_pond = np.sum(MTI_alpha_pond)
    sum_MTI_beta_pond = np.sum(MTI_beta_pond)
    STI = sum_MTI_alpha_pond - sum_MTI_beta_pond
    return STI


def rir2sti(rir, vol, s_pos, m_pos, octave, speaker_gender, Lsf):
    # global Lsf
    # #print(Lsf)
    # #octave = pra.acoustics.OctaveBandsFactory(base_frequency=125, fs=fs, n_fft=512)
    # if vol < 250:
    #     Lsf = Lsf
    # elif vol >= 250:
    #     Lsf = Lsf + 10
    n_octave_bands = len(octave.centers)
    r = np.linalg.norm(s_pos - m_pos)

    rt = []
    for j in range(n_octave_bands):
        rir_octave_bands = octave.analysis(rir, band=j)
        rt_octave_bands = pra.experimental.rt60.measure_rt60(rir_octave_bands, fs=16000, decay_db=60, plot=False, rt60_tgt=None)
        rt.append(rt_octave_bands)
    #print("\nSchroeder RT60 in octave bands at position {} (s):\n".format(1), np.round(rt, 2))

    # Critical Distances (7)
    cost = (2 * vol / c)
    rt = np.array(rt)
    rcf = np.sqrt(np.multiply(cost, rt))
    #print("\nCritical distances at position {} (m):\n".format(i+1), np.round(rcf, 2))

    # Modulation Transfer Function matrix (14 x 7 = 98)
    MTF = []
    for k in range(len(F)):
        for l in range(len(rt)):
            mfF = modulation_depth_reduction_factor(r, Qf, rcf[l], F[k], rt[l], Lsf[l], Lnf[l])
            MTF.append(mfF)
    MTF = np.reshape(MTF, (14, 7))
    #print("\nMTF calculated by Schroeder equation IEC 60268-16 at position {}:\n".format(1), np.round(MTF, 3))

    # Modulation Transfer Indexes (7)
    MTI = np.mean(MTF, axis=0)
    #print("\nModulation Transfer Indexes at position {}:\n".format(1), np.round(MTI, 3))

    # Speech Transmission Index
    if speaker_gender == 'male':
        STI_simul = Speech_Transmssion_Index(MTI, alpha_male, beta_male)
    elif speaker_gender == 'female':
        STI_simul = Speech_Transmssion_Index(MTI, alpha_female, beta_female)
    #print(f'STI: {STI_simul}')
    return STI_simul
    # STI_simul_list.append(STI_simul)

    # STI_P1_simul = STI_simul_list[0]
    # print("\nSTI P1 =", np.round(STI_P1_simul, 2))
    # STI_P2_simul = STI_simul_list[1]
    # print("STI P2 =", np.round(STI_P2_simul, 2))
    # STI_P3_simul = STI_simul_list[2]
    # print("STI P3 =", np.round(STI_P3_simul, 2))
    # STI_P4_simul = STI_simul_list[3]
    # print("STI P4 =", np.round(STI_P4_simul, 2))

    # STI_avg_simul = (sum(STI_simul_list))/(len(STI_simul_list))