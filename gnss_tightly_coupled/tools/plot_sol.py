# coding=utf-8
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.dates as mdate
from parse_sol import parse_sol
import cmn_func as cmn

COLOR1 = (77/256, 133/256, 189/256)
COLOR2 = (247/256, 144/256, 61/256)
COLOR3 = (89/256, 169/256, 90/256)
MARK_SIZE = 1.5

START_POS = -1
ID = START_POS+1
STAT = START_POS+2
OBS = START_POS+3
AZ = START_POS+4
EL = START_POS+5
RES_P1 = START_POS+6
RES_P2 = START_POS+7
RES_P3 = START_POS+8
RES_L1 = START_POS+9
RES_L2 = START_POS+10
RES_L3 = START_POS+11
AMB1 = START_POS+12
AMB2 = START_POS+13
AMB3 = START_POS+14
MW_AMB = START_POS+15
TRPH = START_POS+16
TROW = START_POS+17
MAPH = START_POS+18
MAPW = START_POS+19
ION = START_POS+20
LOCK = START_POS+21
OUTC = START_POS+22


def extract_sat_info(data, index, sat_id):
    temp_t = []
    temp_el = []
    temp_info = []
    for row in data:
        i = 0
        for l in row['epoch']:
            # if row['epoch'][4] == 'NONE':
            #     break
            i = i + 1
            if i == 1:
                temp_t.append(l)
            else:
                continue

    for row in data:
        for l in row['sat']:
            if sat_id == l[0]:
                j = 0
                for s in l:
                    j = j + 1
                    if j == index:
                        temp_info.append(s)
                    elif j == EL:
                        temp_el.append(s)
                    else:
                        continue
            else:
                continue

    t = []
    info = []
    el = []
    for i in range(0, temp_info.__len__()):
        if temp_info[i] is not None:
            info.append(float(temp_info[i]))
            t.append(temp_t[i])
            el.append(float(temp_el[i]))

    return [t, el, info]


def dop_plot(dir, site, data):
    dir_site = os.path.join(dir, site)
    if not os.path.exists(dir_site):
        os.makedirs(dir_site)

    t = []
    pdop = []
    total_sat_num = []
    use_sat_num = []

    for row in data:
        i = 0
        for l in row['epoch']:
            i = i + 1
            if i == 1:
                t.append(l)
            elif i == 6:
                total_sat_num.append(l)
            elif i == 7:
                use_sat_num.append(l)
            elif i == 8:
                pdop.append(l)
            else:
                continue

    fig = plt.figure(figsize=(4, 3))
    ax = fig.add_subplot(111)
    ax.xaxis.set_major_formatter(mdate.DateFormatter('%H'))

    plt.plot_date(t, total_sat_num, '-', color=COLOR1, label='total')
    plt.plot_date(t, use_sat_num, '-', color=COLOR2, label='used')
    plt.legend(loc='best')

    plt.xlabel('Time [h]')
    plt.ylabel('Satellite Number')

    ax1 = ax.twinx()
    ax1.xaxis.set_major_formatter(mdate.DateFormatter('%H'))
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.plot_date(t, pdop, '-', color=COLOR3, label='PDOP')
    ax1.tick_params(axis='y', colors=COLOR3)
    plt.ylabel('PDOP')

    # plt.show()
    save_path = os.path.join(dir_site, 'dop.png')
    fig.savefig(save_path, bbox_inches='tight', dpi=300)
    plt.close(fig)
    plt.clf()


def enu_plot(dir, site, data, pre):
    dir_site = os.path.join(dir, site)
    if not os.path.exists(dir_site):
        os.makedirs(dir_site)

    t = []
    e = []
    n = []
    u = []
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

    [time_e, rms_e, std_e] = cmn.seri_cover(e, t, 0.1)
    [time_n, rms_n, std_n] = cmn.seri_cover(n, t, 0.1)
    [time_u, rms_u, std_u] = cmn.seri_cover(u, t, 0.1)
    # label_x = 'E(cm): RMS={:.2f} STD={:.2f}, COV_T(min): {:.1f}'.format(rms_e*100, std_e*100, time_e)
    # label_y = 'N(cm): RMS={:.2f} STD={:.2f}, COV_T(min): {:.1f}'.format(rms_n*100, std_n*100, time_n)
    # label_z = 'U(cm): RMS={:.2f} STD={:.2f}, COV_T(min): {:.1f}'.format(rms_u*100, std_u*100, time_u)
    label_x = 'E(cm): RMS={:.2f}, COV_T(min): {:.1f}'.format(rms_e*100, time_e)
    label_y = 'N(cm): RMS={:.2f}, COV_T(min): {:.1f}'.format(rms_n*100, time_n)
    label_z = 'U(cm): RMS={:.2f}, COV_T(min): {:.1f}'.format(rms_u*100, time_u)

    fig = plt.figure(figsize=(2.6, 2.4))
    ax = fig.add_subplot(111)
    ax.grid(axis='x', linestyle='-.', linewidth=0.2)
    ax.grid(axis='y', linestyle='-.', linewidth=0.2)
    ax.xaxis.set_major_formatter(mdate.DateFormatter('%H'))
    plt.plot_date(t, e, color=COLOR1, ms=0.5, label=label_x)
    plt.plot_date(t, n, color=COLOR2, ms=0.5, label=label_y)
    plt.plot_date(t, u, color=COLOR3, ms=0.5, label=label_z)

    if 'SPP' in dir:
        plt.ylim([-3.0, 3.0])
    elif "PPP" in dir and "SF" in dir:
        plt.ylim([-1.5, 1.5])
        ax.axhline(y=0.3, color='y', linestyle='-', linewidth=0.4)
        ax.axhline(y=-0.3, color='y', linestyle='-', linewidth=0.4)
    else:
        plt.ylim([-0.05, 0.05])
        ax.axhline(y=0.1, color='y', linestyle='-', linewidth=0.4)
        ax.axhline(y=-0.1, color='y', linestyle='-', linewidth=0.4)
    legend_font = {'size': 6}
    plt.legend(prop=legend_font)
    plt.xlabel('Time [h]')
    plt.ylabel('Position Error [m]')
    if pre == '':
        tlt = site
        save_name = 'pos.png'
    else:
        tlt = site + ' ' + pre
        save_name = pre + '_pos.png'
    plt.title(tlt)
    # plt.show()

    save_path = os.path.join(dir_site, save_name)
    fig.savefig(save_path, bbox_inches='tight', dpi=300)
    plt.close(fig)
    plt.clf()


def plot_isb():
    a = 1


def trp_plot(dir, site, data):
    dir_site = os.path.join(dir, site)
    if not os.path.exists(dir_site):
        os.makedirs(dir_site)

    t = []
    dry = []
    wet = []
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
        for l in row['trp']:
            i = i + 1
            if i == 1:
                dry.append(l)
            elif i == 2:
                wet.append(l)
            else:
                continue

    fig = plt.figure(figsize=(4, 3))
    ax = fig.add_subplot(111)
    ax.xaxis.set_major_formatter(mdate.DateFormatter('%H'))
    ax.grid(axis='x', linestyle='-.', linewidth=0.2)
    ax.grid(axis='y', linestyle='-.', linewidth=0.2)
    plt.plot_date(t, dry, '-', color=COLOR1)

    plt.xlabel('Time [h]')
    plt.ylabel('Dry [m]')

    ax1 = ax.twinx()
    ax1.xaxis.set_major_formatter(mdate.DateFormatter('%H'))
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.plot_date(t, wet, '-', color=COLOR3)
    ax1.tick_params(axis='y', colors=COLOR3)
    plt.ylabel('Wet [m]')

    # plt.show()
    save_path = os.path.join(dir_site, 'trp.png')
    fig.savefig(save_path, bbox_inches='tight', dpi=300)
    plt.close(fig)
    plt.clf()


def sat_plot(dir, site, data, index, sat_id):
    dir_site = os.path.join(dir, site)
    if not os.path.exists(dir_site):
        os.makedirs(dir_site)

    pre = ''
    if index == RES_P1:
        pre = '_res_P1'
        ylabel = 'P1 residual [m]'
    elif index == RES_P2:
        pre = '_res_P2'
        ylabel = 'P2 residual [m]'
    elif index == RES_P3:
        pre = '_res_P3'
        ylabel = 'P3 residual [m]'
    elif index == RES_L1:
        pre = '_res_L1'
        ylabel = 'L1 residual [m]'
    elif index == RES_L2:
        pre = '_res_L2'
        ylabel = 'L2 residual [m]'
    elif index == RES_L3:
        pre = '_res_L3'
        ylabel = 'L3 residual [m]'
    elif index == AMB1:
        pre = '_amb1'
        ylabel = 'L1 ambiguity [m]'
    elif index == AMB2:
        pre = '_amb2'
        ylabel = 'L2 ambiguity [m]'
    elif index == AMB3:
        pre = '_amb3'
        ylabel = 'L3 ambiguity [m]'
    elif index == MW_AMB:
        pre = '_mw'
        ylabel = 'MW ambiguity [m]'

    [t, el, info] = extract_sat_info(data, index, sat_id)
    fig = plt.figure(figsize=(4, 3))
    ax = fig.add_subplot(111)
    ax.xaxis.set_major_formatter(mdate.DateFormatter('%H'))
    ax.grid(axis='x', linestyle='-.', linewidth=0.2)
    ax.grid(axis='y', linestyle='-.', linewidth=0.2)
    plt.plot_date(t, info, ms=MARK_SIZE, color=COLOR1)
    plt.ylabel(ylabel)

    ax1 = ax.twinx()
    ax1.xaxis.set_major_formatter(mdate.DateFormatter('%H'))
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.plot_date(t, el, ms=MARK_SIZE, color=COLOR3)
    ax1.tick_params(axis='y', colors=COLOR3)
    plt.ylabel('Elevation [deg]')
    plt.xlabel('Time [h]')

    # plt.show()

    save_path = os.path.join(dir_site, sat_id + pre + '.png')
    fig.savefig(save_path, bbox_inches='tight', dpi=300)
    plt.close(fig)
    plt.clf()


def walk_files(dir):
    for root, dirs, files in os.walk(dir):
        for d in dirs:
            pos_dir = os.path.join(root, d)
            print(pos_dir)
            for rs, ds, fs in os.walk(pos_dir):
                for f in fs:
                    pos_path = os.path.join(pos_dir, f)
                    # print(pos_path)
                    if pos_path[-4:] == '.pos':
                        site = pos_path[-8:-4]
                        data = parse_sol(pos_dir, site)
                        if ENU_PLOT:
                            enu_plot(pos_dir, site, data, d)
                        if DOP_PLOT:
                            dop_plot(pos_dir, site, data)
                        if TRP_PLOT:
                            trp_plot(pos_dir, site, data)
                        if SAT_PLOT:
                            sat_plot(pos_dir, site, data, RES_P1, SAT_ID)
                            sat_plot(pos_dir, site, data, RES_L1, SAT_ID)
                            sat_plot(pos_dir, site, data, AMB1, SAT_ID)
                            sat_plot(pos_dir, site, data, MW_AMB, SAT_ID)
                    else:
                        continue


if __name__ == '__main__':
    BATCH_PLOT = 1
    # 指定路径,若路径后４位为'.pos',则处理单文件，否则处理该路径下的所有.pos文件
    sol_dir = '/home/wlzhang/software/PPPLib/dataset/data_mgex/result_wum/PPP_STATIC/G_DF_IF/jfng_F.pos'
    # sol_dir = '/home/cc/dataset/data_batch/2020/2092/041/result_wum/PPP_KINE'
    # sol_dir = '/home/cc/dataset/data_batch/2020/2092/041/result_wum/PPP_KINE/G_DF_IF/abmf.pos'

    # sol_dir = 'E:\\PhdWorks\\DataProcess\\PPPLib-Dataset\\data_batch\\2020\\2091\\038\\result_wum\\PPP_STATIC'
    if sol_dir[-4:] == '.pos':
        BATCH_PLOT = 0
    else:
        BATCH_PLOT = 1

    # 指定绘制类型
    DOP_PLOT = 0  # 绘制卫星数量和pdop
    ENU_PLOT = 1  # 绘制ENU位置误差
    TRP_PLOT = 0  # 绘制对流层延迟
    # 绘制卫星有关，指定卫星id后绘制单颗卫星，不指定则不绘制 　
    SAT_PLOT = 0
    SAT_ID = 'G16'

    if BATCH_PLOT:
        walk_files(sol_dir)
    else:
        pos_dir = ''
        site = ''
        if sol_dir[-4:] == '.pos':
            pos_dir = sol_dir[:-10]
            site = sol_dir[-10:-4]
            data = parse_sol(pos_dir, site)
        if ENU_PLOT:
            enu_plot(pos_dir, site, data, '')
        if DOP_PLOT:
            dop_plot(pos_dir, site, data)
        if TRP_PLOT:
            trp_plot(pos_dir, site, data)
        if SAT_PLOT:
            sat_plot(pos_dir, site, data, RES_P1, SAT_ID)
            sat_plot(pos_dir, site, data, RES_L1, SAT_ID)
            sat_plot(pos_dir, site, data, AMB1,   SAT_ID)
            sat_plot(pos_dir, site, data, MW_AMB, SAT_ID)




