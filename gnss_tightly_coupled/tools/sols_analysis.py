import os
import math as m
import matplotlib.pyplot as plt
from parse_sol import parse_sol
import cmn_func as cmn


def analysis_pos(dir, site, data, ana_handle):
    dir_site = os.path.join(dir, site)
    if not os.path.exists(dir_site):
        os.makedirs(dir_site)

    t = []
    e = []
    n = []
    u = []
    enu = []
    for row in data:
        i = 0
        for l in row['epoch']:
            if row['epoch'][4] == 'NONE':
                break
            i = i + 1
            if i == 1:
                t.append(l)
            else:
                continue

    for row in data:
        i = 0
        for l in row['pos']:
            i = i + 1
            if i == 1:
                e.append(l)
            elif i == 2:
                n.append(l)
            elif i == 3:
                u.append(l)
            else:
                continue

    if e.__len__() <= 0 or n.__len__() <= 0 or u.__len__() <= 0 or t.__len__() <= 0:
        return

    for i in range(1, e.__len__()):
        enu.append(m.sqrt(e[i]*e[i]+n[i]*n[i]+u[i]*u[i]))

    [time_e, rms_e, std_e] = cmn.seri_cover(e, t, 0.10)
    [time_n, rms_n, std_n] = cmn.seri_cover(n, t, 0.10)
    [time_u, rms_u, std_u] = cmn.seri_cover(u, t, 0.10)
    [time_3d, rms_3d, std_3d] = cmn.seri_cover(enu, t, 0.10)
    buff = '{:4s}, {:.2f}, {:.2f}, {:.2f}, {:.2f}, {:.1f}, {:.1f}, {:.1f}, {:.1f}\n'.format(site, rms_e*100, rms_n*100, rms_u*100, rms_3d*100,
                                                                                            time_e, time_n, time_u, time_3d)
    if rms_3d > 0.15 or time_3d > 100:
        print('', buff)
        return

    if ana_handle is not '':

        ana_handle.write(buff)
    else:
        label_x = 'E(cm): RMS={:.2f} STD={:.2f}, COV_T(min): {:.1f}'.format(rms_e*100, std_e*100, time_e)
        label_y = 'N(cm): RMS={:.2f} STD={:.2f}, COV_T(min): {:.1f}'.format(rms_n*100, std_n*100, time_n)
        label_z = 'U(cm): RMS={:.2f} STD={:.2f}, COV_T(min): {:.1f}'.format(rms_u*100, std_u*100, time_u)
        print(label_x)
        print(label_y)
        print(label_z)


def walk_files(dir, ana_dir):
    if not os.path.exists(ana_dir):
        os.makedirs(ana_dir)

    for root, dirs, files in os.walk(dir):
        for d in dirs:
            if d == 'analysis':
                continue
            if len(d) == 4:
                continue
            pos_dir = os.path.join(root, d)
            ana_path = os.path.join(ana_dir, d+'_analysis.txt')
            ana_handle = open(ana_path, mode='w')
            print(pos_dir)
            for rs, ds, fs in os.walk(pos_dir):
                for f in fs:
                    pos_path = os.path.join(pos_dir, f)
                    # print(pos_path)
                    if pos_path[-4:] == '.pos':
                        site = pos_path[-8:-4]
                        data = parse_sol(pos_dir, site)
                        analysis_pos(pos_dir, site, data, ana_handle)
                    else:
                        continue
    ana_handle.close()


if __name__ == '__main__':
    BATCH_ANA = 1
    main_dir = '/home/cc/dataset/data_batch/2020/2092/041/result_wum/PPP_KINE'
    sub_dir = ''
    if sub_dir[-4:] == '.pos':
        BATCH_ANA = 0
    else:
        BATCH_ANA = 1
    ana_dir = os.path.join(main_dir, 'analysis')
    sol_dir = os.path.join(main_dir, sub_dir)
    if BATCH_ANA:
        walk_files(sol_dir, ana_dir)
    else:
        pos_dir = ''
        site = ''
        if sol_dir[-4:] == '.pos':
            pos_dir = sol_dir[:-8]
            site = sol_dir[-8:-4]
            data = parse_sol(pos_dir, site)
            analysis_pos(pos_dir, site, data, '')
