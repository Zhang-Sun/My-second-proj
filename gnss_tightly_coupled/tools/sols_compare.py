import os
import math as m
import matplotlib.pyplot as plt
import matplotlib.dates as mdate
from parse_sol import parse_sol
import cmn_func as cmn


def get_pos(data):
    t1 = []
    e1 = []
    n1 = []
    u1 = []
    for row in data:
        i = 0
        for l in row['epoch']:
            if row['epoch'][4] == 'NONE':
                break
            i = i + 1
            if i == 1:
                t1.append(l)
            else:
                continue

    for row in data:
        i = 0
        for l in row['pos']:
            i = i + 1
            if i == 1:
                e1.append(l)
            elif i == 2:
                n1.append(l)
            elif i == 3:
                u1.append(l)
            else:
                continue

    if e1.__len__() <= 0 or n1.__len__() <= 0 or u1.__len__() <= 0 or t1.__len__() <= 0:
        return

    return [t1, e1, n1, u1]


def compare_pos(site, data, n, save_path, sub1, sub2):
    labs = [sub1, sub2]
    pos = []
    for i in range(0, n):
        pos.append(get_pos(data[i]))

    fig = plt.figure(figsize=(4, 3))
    ax1 = fig.add_subplot(311)
    plt.title(site)
    ax1.xaxis.set_major_formatter(mdate.DateFormatter('%H'))
    plt.ylim([-0.2, 0.2])
    plt.ylabel('North [m]')
    plt.xticks([])
    for i in range(0, n):
        ax1.plot_date(pos[i][0], pos[i][1], ms=0.5, label=labs[i])
    ax2 = fig.add_subplot(312)
    ax2.xaxis.set_major_formatter(mdate.DateFormatter('%H'))
    plt.ylim([-0.2, 0.2])
    plt.ylabel('East [m]')
    plt.xticks([])
    for i in range(0, n):
        ax2.plot_date(pos[i][0], pos[i][2], ms=0.5, label=labs[i])
    ax3 = fig.add_subplot(313)
    ax3.xaxis.set_major_formatter(mdate.DateFormatter('%H'))
    plt.ylim([-0.2, 0.2])
    plt.xlabel('Time [h]')
    plt.ylabel('Up [m]')
    for i in range(0, n):
        ax3.plot_date(pos[i][0], pos[i][3], ms=0.5, label=labs[i])

    legend_font = {'size': 6}
    plt.legend(prop=legend_font)

    save_name = save_path + '.png'
    fig.savefig(save_name, bbox_inches='tight', dpi=300)
    plt.close(fig)
    plt.clf()


if __name__ == '__main__':
    # main_dir = '/home/cc/dataset/data_batch/2020/2092/041/result_wum/PPP_KINE'
    main_dir = '/home/cc/dataset/data_mgex/result_wum/PPP_KINE'
    sub_dir1 = 'GR_DF_IF'
    sub_dir2 = 'GR_DF_UC'
    site = 'wuh2'
    data = []
    n_data = 0

    comp_dir1 = os.path.join(main_dir, sub_dir1)
    comp_dir2 = os.path.join(main_dir, sub_dir2)
    comp_save_dir = os.path.join(main_dir, 'compare')

    if not os.path.exists(comp_save_dir):
        os.makedirs(comp_save_dir)
    save_path = os.path.join(comp_save_dir, sub_dir1 + '-' + sub_dir2 + '-' + site)

    for root, dirs, files in os.walk(comp_dir1):
        pos_dir1 = comp_dir1
        for f in files:
            if site in f:
                pos_path1 = os.path.join(comp_dir1, f)
                data1 = parse_sol(pos_dir1, site)
                n_data += 1
                data.append(data1)
                break
            else:
                continue

    for root, dirs, files in os.walk(comp_dir2):
        pos_dir2 = comp_dir2
        for f in files:
            if site in f:
                pos_path2 = os.path.join(comp_dir2, f)
                data2 = parse_sol(pos_dir2, site)
                n_data += 1
                data.append(data2)
                break
            else:
                continue

    compare_pos(site, data, n_data, save_path, sub_dir1, sub_dir2)
