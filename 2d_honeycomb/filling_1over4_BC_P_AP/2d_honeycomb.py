#!/usr/bin/env python

# coding:utf-8
from __future__ import print_function
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def ene(orb,kx,ky):
    enecmplx = -(1.0 + np.exp(-1j*kx) + np.exp(-1j*ky))
    ene = np.abs(enecmplx)
    sign = 1.0-2.0*orb
    return sign*ene

def calc_k_ene(Lorb,Lx,Ly,BCx,BCy):
    if BCx == 'AP' or BCx == 'antiperiodic':
        xshift = 0.5
    elif BCx == 'P' or BCx == 'periodic':
        xshift = 0.0
    else:
        xshift = 0.0
    if BCy == 'AP' or BCy == 'antiperiodic':
        yshift = 0.5
    elif BCy == 'P' or BCy == 'periodic':
        yshift = 0.0
    else:
        yshift = 0.0
    list_kx = np.array([2.0*np.pi*((x+xshift)/Lx-float(Lx//2)/Lx) for x in range(Lx)])
    list_ky = np.array([2.0*np.pi*((y+yshift)/Ly-float(Ly//2)/Ly) for y in range(Ly)])
    list_enekxky = np.array([ene(orb,kx,ky) for ky in list_ky for kx in list_kx for orb in range(Lorb)])
    list_intkxky = np.array([Lorb*(Lx*y+x)+orb for y in range(Ly) for x in range(Lx) for orb in range(Lorb)])
    return list_enekxky, list_intkxky, xshift, yshift

def calc_shell_cond(Lorb,Lx,Ly,BCx,BCy,filling_numer,filling_denom):
    filling = float(filling_numer)/float(filling_denom)
    numel = Lx*Ly*Lorb*filling_numer//filling_denom
    list_enekxky, list_intkxky, xshift, yshift = calc_k_ene(Lorb,Lx,Ly,BCx,BCy)
    list_ind = np.argsort(list_enekxky)
    list_sorted_enekxky = list_enekxky[list_ind]
    list_sorted_intkxky = list_intkxky[list_ind]
    chemipo = 0.5*(list_sorted_enekxky[numel] + list_sorted_enekxky[numel-1])
    totene = np.sum(list_sorted_enekxky[:numel])
    gap = list_sorted_enekxky[numel] - list_sorted_enekxky[numel-1]
    if np.abs(gap)>1e-10:
        shellcond = 'closed'
    else:
        shellcond = 'open'
    return filling, numel, chemipo, totene, gap, shellcond, \
        list_sorted_enekxky, list_sorted_intkxky, xshift, yshift

def main():
    Lorb = 2
    BCx = 'P'
    BCy = 'AP'
#    BCy = 'P'
    filling_numer = 1
#    filling_denom = 2 # half fill
    filling_denom = 4 # quarter fill
    list_L = []
    list_enedens = []
    file = open("dat_2d_honeycomb",'w')
    file.write("# L filling(=n/2) BCx BCy num_electrons(=nup=ndown) chemi_potential ene ene_dens gap shell_cond\n")
    for L in range(2,50,2):
        Lx = L
        Ly = L
        filling, numel, chemipo, totene, gap, shellcond, \
            list_enekxky, list_intkxky, xshift, yshift = \
            calc_shell_cond(Lorb,Lx,Ly,BCx,BCy,filling_numer,filling_denom)
        list_L.append(L)
        list_enedens.append(totene/Lx/Ly/Lorb)
        file.write("{} {} {} {} {} {} {} {} {} {}\n".format(\
            L,filling,BCx,BCy,numel,chemipo,totene,totene/Lx/Ly/Lorb,gap,shellcond))
    file.close()

    list_L = np.array(list_L)
    list_enedens = np.array(list_enedens)
    plt.xlabel("1/L^2")
    plt.ylabel("E/L^2/Lorb")
    plt.plot(1.0/list_L**2,list_enedens,color='blue',marker='o',markerfacecolor='white')
    plt.savefig("fig_2d_honeycomb_enedens.png")
    plt.cla()
    plt.clf()

    L = 30
#    L = 32
    Lx = L
    Ly = L
    filling, numel, chemipo, totene, gap, shellcond, \
        list_enekxky, list_intkxky, xshift, yshift = \
        calc_shell_cond(Lorb,Lx,Ly,BCx,BCy,filling_numer,filling_denom)
    list_intorb = list_intkxky%Lorb
    list_intkx = (list_intkxky//Lorb)%Lx
    list_intky = (list_intkxky//Lorb)//Lx
    list_kx = (list_intkx.astype(np.float64)+xshift)/Lx-float(Lx//2)/Lx
    list_ky = (list_intky.astype(np.float64)+yshift)/Ly-float(Ly//2)/Ly
    plt.xlabel("kx/pi")
    plt.ylabel("ky/pi")
    plt.xticks([-0.5,-0.25,0,0.25,0.5])
    plt.yticks([-0.5,-0.25,0,0.25,0.5])
    plt.xlim(-0.55,0.55)
    plt.ylim(-0.55,0.55)
## https://stackoverflow.com/questions/17990845/how-to-equalize-the-scales-of-x-axis-and-y-axis-in-python-matplotlib
#    plt.axis('equal')
    plt.gca().set_aspect('equal',adjustable='box')
    plt.plot(list_kx,list_ky,color='blue',marker='o',\
        markerfacecolor='white',linestyle='None')
    if numel > Lx*Ly:
        plt.plot(list_kx[Lx*Ly:numel],list_ky[Lx*Ly:numel],color='blue',marker='o',\
            markerfacecolor='blue',linestyle='None')
    else:
        plt.plot(list_kx[:numel],list_ky[:numel],color='blue',marker='o',\
            markerfacecolor='blue',linestyle='None')
    plt.savefig("fig_2d_honeycomb_fermisurface.png")
    plt.cla()
    plt.clf()

    L = 2**9
    Lx = L
    Ly = L
    nbins = L//2
    filling, numel, chemipo, totene, gap, shellcond, \
        list_enekxky, list_intkxky, xshift, yshift = \
        calc_shell_cond(Lorb,Lx,Ly,BCx,BCy,filling_numer,filling_denom)
    plt.xlabel("E")
    plt.ylabel("DOS")
    plt.hist(list_enekxky-chemipo,bins=nbins,density=True)
    plt.savefig("fig_2d_honeycomb_dos.png")

if __name__ == "__main__":
    main()
