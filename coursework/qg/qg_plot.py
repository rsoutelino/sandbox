#!/usr/local/epd-7.2-1-rh5-x86_64/bin/python
# -*- coding:utf-8 -*-

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.mlab import griddata
from mpl_toolkits.basemap import Basemap
import seawater.csiro as sw
import scipy.io as sp
import romslab

################################################################################



def plot_instability_conditions(v, dvdz, dqds, depth, figname):
    # making shaded regions for BC and NBUC
    f = romslab.near(v, 0)[0][0]
    a = depth[f] * -1
    iy = depth[:f+1] * -1
    ix = v[:f+1]
    verts = [(0, 0)] + zip(ix, iy) + [(0, a)]
    poly1 = Polygon(verts, facecolor='r', edgecolor='k', alpha=0.2)

    iy = depth[f:] * -1
    ix = v[f:]
    verts = [(0, a)] + zip(ix, iy) + [(0, depth.max()*-1)]
    poly2 = Polygon(verts, facecolor='b', edgecolor='k', alpha=0.2)

    fig, axarr = plt.subplots(1, 3, sharey=True)
    fig.set_size_inches(12,8)
    fig.set_facecolor('w')
    p1, p2, p3 = axarr[0], axarr[1], axarr[2]
    # fig.tight_layout()

    p1.add_patch(poly1)
    p1.add_patch(poly2)
    p1.grid()
    p1.set_axis_bgcolor('0.95')
    p1.plot([0, 0], [-1600, 0], 'k')
    p1.plot(v, -depth, 'k', linewidth=2)
    p1.set_title(r"$\bar{v}$", fontsize=22)
    p1.set_ylabel("z [m]", fontsize=12)
    p1.set_xlabel(r"cm s$^{-1}$", fontsize=12)
    p1.text(-0.08*100, -100, "BC")
    p1.text(0.03*100, -600, "NBUC")

    p2.set_axis_bgcolor('0.95')
    p2.plot([0, 0], [-1600, 0], 'k')
    p2.grid()
    p2.plot(dvdz, -depth, linewidth=2)
    p2.axis([-0.0013, 0.00091, -1600, 0])
    tit = p2.set_title(r"$\frac{\partial \bar{v}}{\partial z}$", fontsize=22)
    tit.set_position((0.5, 1.02))
    p2.set_xlabel(r"s$^{-1}$", fontsize=12)
    p2.set_xticklabels([0, -4, -2, 0, 2, 4])
    p2.text(0.0004, -1685, "1e-4")

    p3.set_axis_bgcolor('0.95')
    p3.grid()
    p3.plot([0, 0], [-1600, 0], 'k')
    p3.plot(dqds, -depth, linewidth=2)
    p3.axis([-6.9597e-10, 6.49504e-10, -1600, 0])
    tit = p3.set_title(r"$\frac{\partial \bar{q}}{\partial n}$", fontsize=22)
    tit.set_position((0.5, 1.02))
    p3.set_xlabel(r"m$^{-1}$ s$^{-1}$", fontsize=12)
    plt.savefig(figname)


def plot_average_fields(lon, lat, Ubar, Vbar, VbarS, h, zlevs, depths, lims, M, figname):
    matplotlib.rcParams.update({'font.size': 8})
    mlon, mlat = M(lon, lat)
    fig, axarr = plt.subplots(1, 3)
    fig.set_size_inches(12,8)
    fig.set_facecolor('w')
    p1, p2, p3 = axarr[0], axarr[1], axarr[2]
    sc = 3
    d = 4
    # horiz plots
    for k in range(len(depths)):
        a = k+1
        # con = M.contourf(mlon, mlat, h, 30, cmap=plt.cm.Blues, alpha=0.2, ax=axarr[a])
        vel = np.sqrt(Ubar[k]**2 + Vbar[k]**2)
        cmap = plt.cm.Reds if k == 0 else plt.cm.Blues
        M.pcolormesh(mlon, mlat, vel, cmap=cmap, ax=axarr[a])
        # M.contour(mlon, mlat, h, [200, 1000], colors='k', ax=axarr[a], alpha=0.5)
        M.quiver(mlon[::d, ::d], mlat[::d, ::d], Ubar[k][::d, ::d]*sc, Vbar[k][::d, ::d]*sc, scale=10, ax=axarr[a])
        vmax = np.sqrt(Ubar[k]**2 + Vbar[k]**2).max()
        text = '%.0fcm/s' %(vmax*100)
        ax, ay = M( lims[0] + 0.2, lims[-1] - 1.11 ) 
        M.quiver(ax, ay, vmax*sc,0,scale=10,color='k',zorder=10, headwidth=1, headlength=2, ax=axarr[a])
        ax, ay = M( lims[0] + 0.2, lims[-1] - 1.0 ) 
        axarr[a].text(ax, ay, text, color='k', fontsize=8, fontweight='bold')
        text = "Z: %sm" %depths[k]
        ax, ay = M( lims[0] + 0.2, lims[-1] - 1.4 )
        axarr[a].text(ax, ay, text, color='k', fontsize=8, fontweight='bold')
        M.drawcoastlines(zorder=5, ax=axarr[a])
        M.fillcontinents(zorder=3, ax=axarr[a])
        ax, ay = M([-38.7, -36.1], [-19, -19])
        M.plot(ax, ay, 'k', ax=axarr[a], linewidth=3, alpha=0.5)
        M.drawparallels(np.arange(lims[2], lims[3], 1), labels=[1, 0, 0, 0],
                                  dashes=[1,1000], zorder=6, ax=axarr[a] )
        M.drawmeridians(np.arange(lims[0], lims[1], 1), labels=[0, 0, 0, 1],
                                  dashes=[1,1000], zorder=7, ax=axarr[a] )
    
    # vertical section
    a, b, c = VbarS.shape
    xsec = lon[0,...]
    xsec = xsec.repeat(a, axis=0)
    xsec.shape = (c, a)
    xsec = xsec.transpose()

    sec = axarr[0].contourf(xsec, zlevs[:,10,:], VbarS[:, 10, :], 
          np.arange(-0.23, 0.23+0.02, 0.02), cmap=plt.cm.RdBu)
    axarr[0].set_xlim(-38.7, -36.1)
    axarr[0].set_xlabel("Longitude")
    axarr[0].set_ylabel("Z[m]")
    axarr[0].set_aspect(0.0024)
    axarr[0].set_axis_bgcolor('0.8')
    
    ax = fig.add_axes([0.14, 0.35, 0.015, 0.3])
    cbar = plt.colorbar(sec,cax=ax,orientation='vertical')
    cbar.set_label('m/s')
    plt.savefig(figname)


def plot_rossby(lon, lat, M, lims, h, Ro, figname):
    matplotlib.rcParams.update({'font.size': 8})
    mlon, mlat = M(lon, lat)
    fig, axarr = plt.subplots(1, 2)
    fig.set_size_inches(8,6)
    fig.set_facecolor('w')
    for k in range(2):
        cmap = plt.cm.Greens
        pc = M.pcolormesh(mlon, mlat, Ro[k], vmin=0, vmax=0.06, cmap=cmap, ax=axarr[k])
        # M.contour(mlon, mlat, h, [200, 1000], colors='k', ax=axarr[k], alpha=0.5)
        M.drawcoastlines(zorder=5, ax=axarr[k])
        M.fillcontinents(zorder=3, ax=axarr[k])
        ax, ay = M( lims[0] + 0.2, lims[-1] - 0.5 ) 
        text = "BC" if k == 0 else "NBUC"
        axarr[k].text(ax, ay, text, color='k', fontsize=8, fontweight='bold')
        M.drawparallels(np.arange(lims[2], lims[3], 1), labels=[1, 0, 0, 0],
                                  dashes=[1,1000], zorder=6, ax=axarr[k] )
        M.drawmeridians(np.arange(lims[0], lims[1], 1), labels=[0, 0, 1, 0],
                                  dashes=[1,1000], zorder=7, ax=axarr[k] )

    ax = fig.add_axes([0.23, 0.13, 0.54, 0.015])
    cbar = plt.colorbar(pc, cax=ax, orientation='horizontal')
    cbar.set_label('Ro')
    plt.savefig(figname)


def plot_energy_conversion(BC, figname):
    matplotlib.rcParams.update({'font.size': 10})
    fig = plt.figure(figsize=(8,4), facecolor='w')
    plt.plot([0, 40],[0, 0], 'k--')
    plt.plot([10, 10], [-1, 1], 'k', [20, 20], [-1, 1], 'k', [30, 30], [-1, 1],
             'k', linewidth=20, alpha=0.2 )
    plt.gca().set_axis_bgcolor('0.95')
    plt.grid()
    plt.plot(BC[0], 'r', linewidth=2, label='BC')
    plt.plot(BC[1], 'b', linewidth=2, label='NBUC')
    plt.gca().set_xlim([8, 32])
    plt.gca().set_ylim([-3e-12, 3e-12])
    plt.xlabel('Time [days]')
    plt.ylabel('[m$^2$s$^{-3}$]')
    plt.title('Horizontally-averaged Baroclinic Conversion Rate')
    plt.legend(loc=2)
    plt.savefig(figname)