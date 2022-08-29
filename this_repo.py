import sys, os
import numpy as np
from matplotlib import pyplot as plt

# Edit the font, font size, and axes width
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.size'] = 8
plt.rcParams['figure.figsize'] = [4.0, 3.0]
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams["image.interpolation"] = None

print("Executed: `asnlib_mplset.py`")


def mark_line(ax, pos, orientation_dim, pdict={'color':'red'}):
    """[summary]

    Args:
        ax ([type]): [description]
        pos ([tuple] of len 2, containing arrays): [description]
        pdict (dict, optional): [description]. Defaults to {'color':'red'}.

    Returns:
        [type]: [description]
    """

    xlim0 = ax.get_xlim()
    ylim0 = ax.get_ylim()
    ax.set_xlim(xlim0)
    ax.set_ylim(ylim0)

    positions = pos[orientation_dim]
    line_ends = pos[int(not bool(orientation_dim))] #has to be arry of len2 or None

    if orientation_dim == 0: #horizontal line
        vpos = positions[0]
        if (line_ends is not None):
            xlim0 = tuple(line_ends)
        ax.plot(np.array(xlim0), np.ones(2)*vpos, **pdict)
    else: #vertical line
        hpos = positions[0]
        if (line_ends is not None):
            ylim0 = tuple(line_ends)
        ax.plot(np.ones(2)*hpos, np.array(ylim0), **pdict)

    return ax

