#! python3
# -*- coding: UTF-8 -*-

# @Project ：corr.py
# @File ：corr_new.py
# @Author ：Jiangcg
# Email: successjiang@tju.edu.cn
# @Date ：2021-06-26 19:05 
# @Desc :


import os
import sys
from itertools import product
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colorbar import ColorbarBase
import pandas as pd
from plotnine import *


def corr_plot(x_cof, x_name=None, left_down=None, right_top=None, threshold_left=0, threshold_right=0.9,
              title=None, label_axis="off", front_raito=1):
    """
    plot corr
    Parameters
    ----------
    x_cof:np.ndarray
        correlation coefficient matrix
    x_name:list,None
        feature names
    left_down:None,"pie","fill","text","circle"
        type for left_down
    right_top:None,"pie","fill","text","circle"
        type for right_top
    threshold_left:int
        threshold for show.
    threshold_right:int
        threshold for show.
    title:str
        picture title
    label_axis:"left","right","off"
        label_axis kwargs
    front_raito:int
        front scare for show
    """
    x_cof = np.round(x_cof, 2)
    name = x_name
    size = x_cof
    or_size = np.nan_to_num((abs(size) / size) * (1 - abs(size)))
    n = size.shape[0]
    explode = (0, 0)
    gs = gridspec.GridSpec(n, n, width_ratios=[1]*16, height_ratios=[1.5]*16)
    gs.update(wspace=0, hspace=0)
    cmap = plt.get_cmap("bwr")  # args
    fill_colors = cmap(size / 2 + 0.5)  # args
    fig = plt.figure(figsize=(6, 6), frameon=True)  # args
    title_fontsize = round(15 * front_raito)  # c_args
    ax_fontsize = round(15 * front_raito)
    score_fontsize = round(15 * front_raito)
    circle_size = round(400 * front_raito)
    fig.text(0.5, 0.05, title, fontsize=title_fontsize, horizontalalignment='center',
             verticalalignment='center')  # zou, xia
    for i, j in product(range(n), range(n)):
        if j < i and abs(size[i, j]) >= threshold_left:
            types = left_down
        elif j > i and abs(size[i, j]) >= threshold_right:
            types = right_top
        else:
            types = None
        if types is "pie":
            ax = plt.subplot(gs[i, j])
            ax.pie((size[i, j], or_size[i, j]), explode=explode, labels=None, autopct=None, shadow=True,
                   startangle=90,
                   colors=[fill_colors[i, j], 'w'], wedgeprops=dict(width=1, edgecolor='black', linewidth=0.5),
                   counterclock=False,
                   frame=False, center=(0, 0))
            ax.set_xlim(-1, 1)
            ax.axis('equal')
        elif types is "fill":
            ax = plt.subplot(gs[i, j])
            ax.set_facecolor(fill_colors[i, j])
            [ax.spines[_].set_color('w') for _ in ['right', 'top', 'left', 'bottom']]
            ax.text(0.5, 0.5, size[i, j],
                    fontdict={"color": "black"},  # args
                    fontsize=score_fontsize,  # c_arg
                    horizontalalignment='center', verticalalignment='center')
            ax.set_xticks([])
            ax.set_yticks([])
        elif types is "text":
            ax = plt.subplot(gs[i, j])
            if size[i, j] > 0:
                ax.text(0.5, 0.5, size[i, j],
                        fontdict={"color": "r"},  # args
                        fontsize=score_fontsize,  # c_arg
                        horizontalalignment='center', verticalalignment='center')
            else:
                ax.text(0.5, 0.5, size[i, j],
                        fontdict={"color": "b"},  # args
                        fontsize=score_fontsize,  # c_arg
                        horizontalalignment='center', verticalalignment='center')
            # ax.axis('equal')
            # ax.set_xlim(-1, 1)
            ax.set_xticks([])
            ax.set_yticks([])
            # plt.axis('off')
        elif types is "circle":
            ax = plt.subplot(gs[i, j])
            ax.axis('equal')
            ax.set_xlim(-1, 1)
            ax.scatter(0, 0, color=fill_colors[i, j], s=circle_size * abs(size[i, j]) ** 2)
            ax.set_xticks([])
            ax.set_yticks([])
            # plt.axis('off')
        else:
            pass
    for k in range(n):
        ax = plt.subplot(gs[k, k])
        ax.text(0.5, 0.5, name[k], fontsize=ax_fontsize, horizontalalignment='center', verticalalignment='center')
        ax.set_xticks([])
        ax.set_yticks([])
        if label_axis is "left":
            color = ["w", "w", "b", "b"]
            [ax.spines[i].set_color(j) for i, j in zip(['right', 'top', 'left', 'bottom'], color)]
        elif label_axis is "right":
            color = ["b", "b", "w", "w"]
            [ax.spines[i].set_color(j) for i, j in zip(['right', 'top', 'left', 'bottom'], color)]
        else:
            plt.axis('off')
    fig.subplots_adjust(right=0.75)
    cbar_ax = fig.add_axes([0.75, 0.70, 0.025, 0.15])
    norm = mpl.colors.Normalize(vmin=-1, vmax=1)
    ColorbarBase(cbar_ax, cmap=cmap, norm=norm)
    fig.set_size_inches(3, 3, forward=True)
    plt.show()


df = pd.read_csv('features.csv')
features_corr = df.corr().values
name0 = df.T.index[1:]
corr_plot(features_corr, name0, left_down="circle", right_top="text", threshold_right=0.9, label_axis="off", front_raito=0.5)
