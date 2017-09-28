# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 16:06:58 2017

@author: wangronin
"""

import pdb
from GaussianProcess import GaussianProcess
from GaussianProcess.trend import constant_trend, linear_trend

from criteria import EI
from scipy.stats import norm

import numpy as np
from numpy import sin

import matplotlib.pyplot as plt
from matplotlib import rcParams

# ------------------------------ matplotlib settings -----------------------------
fig_width = 22
fig_height = fig_width * 9 / 16

plt.ioff()
rcParams['legend.numpoints'] = 1 
rcParams['xtick.labelsize'] = 20
rcParams['ytick.labelsize'] = 20
rcParams['xtick.major.size'] = 7
rcParams['xtick.major.width'] = 1
rcParams['ytick.major.size'] = 7
rcParams['ytick.major.width'] = 1
rcParams['axes.labelsize'] = 50
rcParams['font.size'] = 20
plt.style.use('ggplot')

f = lambda x: x* sin(x) + 1
ub = 2.2
lb = -3

np.random.seed(1)

# ------------------------------ The GP model -----------------------------
thetaL = 1e-5 * (ub - lb) * np.ones(1)
thetaU = 10 * (ub - lb) * np.ones(1)
theta0 = np.random.rand(1) * (thetaU - thetaL) + thetaL

mean = constant_trend(1, beta=0)
model = GaussianProcess(mean=mean,
                        corr='matern',
                        theta0=theta0,
                        thetaL=thetaL,
                        thetaU=thetaU,
                        nugget=None,
                        noise_estim=False,
                        optimizer='BFGS',
                        wait_iter=5,
                        random_start=15,
                        likelihood='concentrated',
                        eval_budget=100)
                        
x_train = np.array([-2.5, -1.2])
y_train = f(x_train)

x = np.linspace(lb, ub, 500)
y = f(x)

add_on = np.array([1.8, -0.4,  0.12, -0.07, 0.001])
for i, p in enumerate(add_on):
    
    fig, ax0 = plt.subplots(1, 1, subplot_kw={'aspect':'auto'}, figsize=(fig_width, fig_height))

    ax0.grid(True)              
    ax0.hold(True)
    ax0.set_xlabel('$x$', fontsize=28)
    ax0.set_ylabel('$f$', fontsize=28)
    ax0.set_xlim([lb, ub])
    ax0.set_ylim([-0.1, 3.5])
    
    lines = []
    yp = f([p])

    lines += ax0.plot(x_train, y_train, ls='none', ms=10, marker='o', mfc='k', mec='k', alpha=0.8)

    x_train = np.r_[x_train, p]
    y_train = np.r_[y_train, yp]
    model.fit(x_train.reshape(-1, 1), y_train)
    
    if i > 0:
        ax0.plot(p, yp, ls='none', ms=10, marker='o', mfc='r', mec='r', alpha=0.8)
    else:
        ax0.plot(p, yp, ls='none', ms=10, marker='o', mfc='k', mec='k', alpha=0.8)

        scale = [1e-1 / 2., 1e-1 * 0.7]
        _min = np.min(y_train)

        xx1 = np.linspace(lb, ub, 500)
        xx2 = np.linspace(-0.6, -0.1, 100)
        ax0.plot(xx1, _min*np.ones(np.size(xx1)), 'k--', lw=1.5)
        
        for k, loc in enumerate([-0.8, -0.4]):
            y_hat, mse = model.predict([[loc]], eval_MSE=True)

            loc_mean = y_hat.flatten()
            loc_sigma = np.sum(np.sqrt(mse))
            x1 = np.linspace(1, loc_mean, 100)
            x2 = np.linspace(loc_mean, 3, 100)
            y1 = scale[k] * norm.pdf(x1, loc=loc_mean, scale=loc_sigma)
            y2 = scale[k] * norm.pdf(x2, loc=loc_mean, scale=loc_sigma)

            max_y = max(y1)
            # ax0.plot((k, max_y), (loc_mean, loc_mean), 'k-', lw=1.5)
            
            truncate_point = scale[k] * norm.pdf(_min, loc=loc_mean, scale=loc_sigma)
            x_truncate = y1[y1 <= truncate_point] + loc
            y_truncate = x1[y1 <= truncate_point]

            ax0.plot(y1 + loc, x1, 'k-', lw=2)
            ax0.plot(y2 + loc, x2, 'k-', lw=2)
            ax0.plot((loc, loc), (1, 3), 'k-', lw=2)

            ax0.fill_between(x_truncate, y_truncate, 
                             _min* np.ones(np.size(x_truncate)), color='r', alpha=.5)
    
    y_hat, mse = model.predict(x.reshape(-1, 1), eval_MSE=True)
    y_hat = y_hat.flatten()
    sd = np.sqrt(mse)
    
    y_up = y_hat + sd * 1.96
    y_down = y_hat - sd * 1.96
    
    acquisition_func = EI(model, None, minimize=True)

    func = acquisition_func(x.reshape(-1, 1))
    func = func / max(func)

    lines += ax0.plot(x, y, ls='-', lw=1.5, color='k')
    lines += ax0.plot(x, y_hat, ls='--', lw=1.5, color='k')
    lines.append(ax0.fill_between(x, y_down, y_up, facecolor='#2673A2', 
                                  interpolate=True, alpha=0.4))
    
    lines += ax0.plot(x, func, ls='-', lw=2, color='r')
    
    ax0.set_title('iteration N={}'.format(i + 1), fontsize=25)
    
    plt.legend(lines, ['data points', 'objective', 'model prediction', '95% CI',
                       'Acquisition function'], fontsize=20)
    
    plt.tight_layout()
    plt.savefig('N{}.png'.format(i + 1), dpi=300)
    # plt.show()

