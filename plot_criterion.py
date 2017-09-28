# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 12:47:36 2017

@author: wangronin
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 14:15:50 2017

@author: wangronin
"""

import pdb

import os
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

from rpy2.robjects import r
from rpy2.robjects import pandas2ri

from matplotlib import rcParams
    
pandas2ri.activate()

# plot settings
plt.ioff()
fig_width = 22
fig_height = 22 * 9 / 16
_color = ['r', 'm', 'b', 'g', 'c']

rcParams['font.size'] = 15
rcParams['xtick.direction'] = 'out'
rcParams['ytick.direction'] = 'out'
rcParams['font.size'] = 15
rcParams['legend.numpoints'] = 1 
rcParams['xtick.labelsize'] = 13
rcParams['ytick.labelsize'] = 13
rcParams['xtick.major.size'] = 7
rcParams['xtick.major.width'] = 1
rcParams['ytick.major.size'] = 7
rcParams['ytick.major.width'] = 1

# setup working directory
os.chdir(os.path.expanduser('~')  + '/Desktop/Pandora')

# ----------------------- Loading the data set from R ------------------------
r['load']('./grad.RData')

X = r['X']
Y = r['Y']
X1 = r['X1']
Y1 = r['Y1']
values = r['values']
#gradient = r['gradient']

# diagostic plots of gradient field 
def plot_contour_gradient(ax, x_lb, x_ub, title='f', is_log=False, n_level=30):
    
    fig = ax.figure
    foo = 1e-30
    
    fitness = values
    try:
        fitness = (fitness - np.min(fitness)) / (np.max(fitness) - np.min(fitness)) + foo
        if is_log:
            fitness = np.log(fitness)
        CS = ax.contour(X, Y, fitness, n_level, cmap=plt.cm.winter, linewidths=1)
#        plt.clabel(CS, inline=1, fontsize=5)
#        fig.colorbar(CS, ax=ax, fraction=0.046, pad=0.04)
        fig.colorbar(CS, ax=ax)
    except:
        pdb.set_trace()
    
#    dx = gradient
#    n_X1 = X1.shape[1]
#    dx_norm = np.sqrt(np.sum(dx ** 2.0, axis=1)) # in case of zero gradients
#    dx /= dx_norm.reshape(-1, 1)
#    dx1 = dx[:, 0].reshape(-1, n_X1)
#    dx2 = dx[:, 1].reshape(-1, n_X1)
    
#    dx = gradient
    n_X1 = X1.shape[1]
    dx1, dx2 = r['dx'], r['dy']
    
    dx_norm = np.sqrt(dx1 ** 2. + dx2 ** 2.) # in case of zero gradients
    dx1 /= dx_norm
    dx2 /= dx_norm
  
    CS = ax.quiver(X1, Y1, dx1, dx2, dx_norm, cmap=plt.cm.jet, 
                   #norm=colors.LogNorm(vmin=1e-100, vmax=dx_norm.max()),
                   headlength=5)
   
#        fig.colorbar(CS, ax=ax)
    
    ax.set_xlabel('$x_1$')
    ax.set_ylabel('$x_2$')
    ax.grid(True)
    ax.set_title(title)
    ax.set_xlim(x_lb[0], x_ub[0])
    ax.set_ylim(x_lb[1], x_ub[1])
    
    plt.show()
    
fig0, ax = plt.subplots(1, 1,
                     figsize=(fig_width, fig_height), 
                     subplot_kw={'aspect':'equal'}, dpi=100)

plot_contour_gradient(ax, (-5, 0), (10, 15), 'MGF(0.8)')