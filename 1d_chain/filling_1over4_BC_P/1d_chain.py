#!/usr/bin/env python

# coding:utf-8
from __future__ import print_function
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def ene(kx):
    return -2.0*np.cos(kx)

def calc_k_ene(Lx,BCx):
    if BCx == 'AP' or BCx == 'antiperiodic':
        xshift = 0.5
    elif BCx == 'P' or BCx == 'periodic':
        xshift = 0.0
    else:
        xshift = 0.0
    list_kx = np.array([2.0*np.pi*((x+xshift)/Lx-float(Lx//2)/Lx) for x in range(Lx)])
    list_enekx = np.array([ene(kx) for kx in list_kx])
    list_intkx = np.array([x for x in range(Lx)])
    return list_enekx, list_intkx, xshift

def calc_shell_cond(Lx,BCx,filling_numer,filling_denom):
    filling = float(filling_numer)/float(filling_denom)
    numel = Lx*filling_numer//filling_denom
    list_enekx, list_intkx, xshift = calc_k_ene(Lx,BCx)
    list_ind = np.argsort(list_enekx)
    list_sorted_enekx = list_enekx[list_ind]
    list_sorted_intkx = list_intkx[list_ind]
    chemipo = 0.5*(list_sorted_enekx[numel] + list_sorted_enekx[numel-1])
    totene = np.sum(list_sorted_enekx[:numel])
    gap = list_sorted_enekx[numel] - list_sorted_enekx[numel-1]
    if np.abs(gap)>1e-10:
        shellcond = 'closed'
    else:
        shellcond = 'open'
    return filling, numel, chemipo, totene, gap, shellcond, \
        list_sorted_enekx, list_sorted_intkx, xshift

def main():
    BCx = 'P'
    filling_numer = 1
    filling_denom = 4
    list_Lx = []
    list_enedens = []
    file = open("dat_1d_chain",'w')
    file.write("# Lx filling(=n/2) BCx num_electrons(=nup=ndown) chemi_potential ene ene_dens gap shell_cond\n")
    for Lx in range(4,60,4):
        filling, numel, chemipo, totene, gap, shellcond, \
            list_enekx, list_intkx, xshift = \
            calc_shell_cond(Lx,BCx,filling_numer,filling_denom)
        list_Lx.append(Lx)
        list_enedens.append(totene/Lx)
        file.write("{} {} {} {} {} {} {} {} {}\n".format(\
            Lx,filling,BCx,numel,chemipo,totene,totene/Lx,gap,shellcond))
    file.close()

    list_Lx = np.array(list_Lx)
    list_enedens = np.array(list_enedens)
    plt.xlabel("1/Lx")
    plt.ylabel("E/Lx")
    plt.plot(1.0/list_Lx,list_enedens,color='blue',marker='o',markerfacecolor='white')
#    plt.axhline(y=-2.0/np.pi,color='red')
    plt.savefig("fig_1d_chain_enedens.png")
    plt.cla()
    plt.clf()

    Lx = 36
    filling, numel, chemipo, totene, gap, shellcond, \
        list_enekx, list_intkx, xshift = \
        calc_shell_cond(Lx,BCx,filling_numer,filling_denom)
    list_kx = (list_intkx.astype(np.float64)+xshift)/Lx-float(Lx//2)/Lx
    list_ky = np.array([0 for x in range(Lx)])
    plt.xlabel("kx/pi")
    plt.ylabel("")
    plt.xticks([-0.5,-0.25,0,0.25,0.5])
    plt.yticks([])
    plt.xlim(-0.55,0.55)
## https://stackoverflow.com/questions/17990845/how-to-equalize-the-scales-of-x-axis-and-y-axis-in-python-matplotlib
#    plt.axis('equal')
#    plt.gca().set_aspect('equal',adjustable='box')
    plt.gca().set_aspect(0.5,adjustable='box')
    plt.plot(list_kx,list_ky,color='blue',marker='o',\
        markerfacecolor='white',linestyle='None')
    plt.plot(list_kx[:numel],list_ky[:numel],color='blue',marker='o',\
        markerfacecolor='blue',linestyle='None')
    plt.savefig("fig_1d_chain_fermisurface.png")
    plt.cla()
    plt.clf() 

    Lx = 2**18
    nbins = int(np.sqrt(Lx))//2
    filling, numel, chemipo, totene, gap, shellcond, \
        list_enekx, list_intkx, xshift = \
        calc_shell_cond(Lx,BCx,filling_numer,filling_denom)
    plt.xlabel("E")
    plt.ylabel("DOS")
    plt.hist(list_enekx-chemipo,bins=nbins,density=True)
    plt.savefig("fig_1d_chain_dos.png")

if __name__ == "__main__":
    main()
