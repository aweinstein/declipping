#!/usr/bin/env python
import sys
import os.path
import shutil
import subprocess
import argparse

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from scipy.io import loadmat
from matplotlib.backends.backend_pdf import PdfPages

params = {'backend': 'Agg',
          'axes.labelsize': 22,
          'text.fontsize': 22,
          'legend.fontsize': 14,
          'xtick.labelsize': 18,
          'ytick.labelsize': 18,
          'savefig.dpi' : 600,
          'ps.usedistiller' : 'xpdf',
          'text.usetex' : True,
          'font.family': 'serif',
          'font.serif' : ['Times'],
          #'axes.linewidth' : 0.5,
          #'xtick.major.size' : 2,
          #'ytick.major.size' : 2,
          #'font.size' : 9}
          }
mpl.rcParams.update(params)

bottom_margin=0.15
markers = 'ovs^p*Dx<>1234hH+d|_'
s_d = 'arc,angleA=0,armA=20,angleB=-90,armB=15,rad=7'
s_u = 'angle3,angleA=0,angleB=-90'
# Map between the name of the methods and the name used in the representations
method_name = {'bp':'BP',
               'bpcc':'BPCC', 
               'rw_ell1':r'R$\ell_1$CC',
               'omp': 'OMPn',
               'omp_n':'OMP',
               'ls_omp':'LS-OMP',
               'astart-omp':'A*-OMP',
               'c_omp': 'C-OMP',
               'isd': 'ISD',
               'rand': 'CS',
               'tpcc': 'TPCC',
               'cvp_omp': 'cvpOMP',
               'consOMP': 'consOMP'}


# List of methods to Plot
plot_methods = [
    'bp',
    'omp_n',
    'isd',
    'rw_ell1',
    'c_omp',
    'rand',
    'tpcc',
    'cvp_omp',
    'bpcc',
    'consOMP',]

def open_file(file_name):
    """Open the Matlab file file_name. Return a dictionary with the data."""
    
    dir_name, file_name = os.path.split(file_name)
    print 'Opening data file', file_name
    return  loadmat(dir_name +'/'+ file_name, struct_as_record=True)

def save_fig(fig, file_name):
    fn = 'pdfs/' + os.path.split(file_name)[1].split('.')[0] + '.pdf'
    if isinstance(fig, list):
        pdf = PdfPages(fn)
        for f in fig:
            pdf.savefig(f)
        pdf.close()
    else:
        fig.savefig(fn, format='pdf', dpi=300)
    print 'File %s created' % fn
    return fn

def footnote(file_name):
    fn = os.path.split(file_name)[1].split('.')[0]
    plt.figtext(1, 0, fn, fontsize=6, horizontalalignment='right')
    
def experiment_1(file_name, args):
    data = open_file(file_name)
    K = data['K'].flatten()
    methods = data['methods'].flatten()
    Ms = data['Ms']
    Cls = data['Cls']
    maxs = data['maxs'].flatten()

    grey_colors = [str(g) for g in np.linspace(0, 0.6, len(methods))]
    style = [s_u, s_u, s_u, s_d, s_u, s_u, s_d, s_d, s_d]
    d = [1, 1, 1, -1, 1, 1, -1, -1, -1]

    figs = []
    # Plot Mmin versus K
    figs.append(plt.figure())
    ax = plt.subplot(111)
    k = 1
    K = 2*K
    for i in np.arange(len(methods)):
        if methods[i][0] not in plot_methods:
            continue
        p = plt.plot(K, Ms[:,i], markers[i] + '-',
                     markersize=9,
                     linewidth=2,
                     color = grey_colors[i],
                     label=method_name[methods[i][0]])
        m = method_name[methods[i][0]]
        #ax.annotate(m, xy=(k+2, Ms[k+1,i]),  xycoords='data',
        #            xytext=(-40, d[i]*30), textcoords='offset points',
        #            arrowprops=dict(arrowstyle="->",
        #                            connectionstyle=style[i]),) 
        k += 1

    plt.legend()
    plt.ylim([0,100])
    plt.xlabel(r'Sparsity level $K$')
    plt.ylabel(r'$M_{min}$')
    if args.footnote:
        footnote(file_name)
    plt.gcf().subplots_adjust(bottom=bottom_margin)

    fn = save_fig(figs, file_name)

    #plt.show()
    return fn

def experiment_2(file_name, args):
    data = open_file(file_name)
   
    K = data['K'].flatten()
    methods = data['methods'].flatten()
    p_recovery = data['p_recovery']
    M = data['M'].flatten()
    grey_colors = [str(g) for g in np.linspace(0, 0.6, len(methods))]
    fig = plt.figure()
    ax = plt.subplot(111)
    
    k = 0
    K = 2*K
    for i in np.arange(len(methods)):
        if methods[i][0] not in plot_methods:
            continue
        p = plt.plot(K, p_recovery[i,:], markers[i] + '-',
                     markersize=9,
                     linewidth=2,
                     color = grey_colors[i])
        m = method_name[methods[i][0]]
        xy_ann = (K[k+6], p_recovery[i,k+6])
        if k == 0:
            ax.annotate(m, xy=xy_ann,  xycoords='data',
                        xytext=(-40, 20), textcoords='offset points',
                        arrowprops=dict(arrowstyle="->",
                                        connectionstyle=s_u,))
        else:
            ax.annotate(m, xy=xy_ann,  xycoords='data',
                        xytext=(-40, -30), textcoords='offset points',
                        arrowprops=dict(arrowstyle="->",
                                        connectionstyle=s_u,))
        k += 1
        
    plt.xlabel(r'Sparsity level $K$')
    plt.ylabel('Probability of recovery')
    #plt.title('Number of measurements, M=%d, n_trials=%d' %
    #          (data['M'], data['n_trials'])) 
    plt.ylim([0, 1.15])
    if args.footnote:
        footnote(file_name)
    plt.gcf().subplots_adjust(bottom=bottom_margin)
    fn = save_fig(fig, file_name)
    #plt.show()
    return fn

def experiment_3(file_name, args):
    data = open_file(file_name)
    
    K = data['K'].flatten()
    M = data['M'].flatten()
    p_recovery = data['p_recovery']
    grey_colors = [str(g) for g in np.linspace(0, 0.6, len(M))]
    fig = plt.figure()
    ax = plt.subplot(111)

    name = method_name[data['method'][0]]
    k = 0
    offsets = [7,9,12,14,17]
    K = 2*K
    for i in np.arange(len(M)):
        p = plt.plot(K, p_recovery[i,:], markers[i] + '-',
                     markersize=9,
                     linewidth=2,
                     color = grey_colors[i],
                     label='M=%d' % M[i])

        offset = offsets[k]
        k += 1
    plt.legend()
    plt.xlabel(r'Sparsity level $K$')
    plt.ylabel('Probability of recovery')
    #plt.title('Number of measurements, M=%d, n_trials=%d' %
    #          (data['M'], data['n_trials'])) 
    plt.ylim([0, 1.1])
    plt.xlim(xmin=min(K))
    if args.footnote:
        footnote(file_name)
    plt.gcf().subplots_adjust(bottom=bottom_margin)
    fn = save_fig(fig, file_name)
    return fn
    
def mangasarian(file_name, args):
    d = open_file(file_name)
    uniqueness = d['uniqueness']
    recovery = d['recovery']
    K = d['K'].flatten()
    M = d['M'].flatten()
    fig = plt.figure()
    ax = plt.subplot(111)
    # BP m=0, CBP m=1, CS-BP m = 2 
    # Lets plot only BP (m = 0)
    m = 0 # Index corresponding to BP
    cm = ('0.1', '0.5', '0.9')
    x_offset = (0,0.38,0)
    y_offset = (0,0,4)
    ms = 50 # marker size
    for m in (0,1,2):
        r0u0_x, r0u0_y, r0u1_x, r0u1_y = [], [], [], []
        r1u0_x, r1u0_y, r1u1_x, r1u1_y = [], [], [], []
        for i in range(uniqueness.shape[0]):
            for j in range(uniqueness.shape[1]):
                u = uniqueness[i,j,m]
                r = recovery[i,j,m]
                if r==0 and u ==0:
                    r0u0_x.append(K[i] + x_offset[m])
                    r0u0_y.append(M[j] + y_offset[m])
                elif r==0 and u==1:
                    r0u1_x.append(K[i] + x_offset[m])
                    r0u1_y.append(M[j] + y_offset[m])
                elif r==1 and u==0:
                    r1u0_x.append(K[i] + x_offset[m])
                    r1u0_y.append(M[j] + y_offset[m])
                elif r==1 and u==1:
                    r1u1_x.append(K[i] + x_offset[m])
                    r1u1_y.append(M[j] + y_offset[m])

        if len(r0u0_x) > 0:
            plt.scatter(r0u0_x, r0u0_y, s=ms, marker='s', c=cm[m],
                        edgecolor='0', label='R=0, U=0')
        if len(r0u1_x) > 0:
            plt.scatter(r0u1_x, r0u1_y, s=ms, marker='^', c=cm[m],
                        edgecolor='0', label='R=0, U=1')
        if len(r1u0_x) > 0:
            plt.scatter(r1u0_x, r1u0_y, s=ms, marker='o', c=cm[m],
                        edgecolor='0', label='R=1, U=0')
        if len(r1u1_x) > 0:
            plt.scatter(r1u1_x, r1u1_y, s=ms, marker='d', c=cm[m],
                        edgecolor='0', label='R=1, U=1')
            

    # Annotate the lower left corner with the names BP, CBP and CS-BP
    x = max(K)
    y = min(M)
    ax.annotate('$(P1)$', (x,y), textcoords='offset points', xytext=(30,-20),
                arrowprops=dict(arrowstyle='->',
                                connectionstyle=('arc3,rad=-0.3')))
    ax.annotate('$(P2)$', (x+x_offset[1],y+y_offset[1]),
                textcoords='offset points', xytext=(30,0),
                arrowprops=dict(arrowstyle='->',
                                connectionstyle=('arc3,rad=-0.3')))
    ax.annotate('CS', (x+x_offset[2],y+y_offset[2]),
                textcoords='offset points', xytext=(30,20),
                arrowprops=dict(arrowstyle='->',
                                connectionstyle=('arc3,rad=-0.3')))
    plt.xlabel('Sparsity level $K$')
    plt.ylabel('Number of measurements $M$')
    plt.xlim([0.8, x+2])
    #plt.legend(loc='upper right',scatterpoints=1, markerscale=0.5)
    if args.footnote:
        footnote(file_name)
    fn = save_fig(fig, file_name)
    plt.show()
    return fn

def test_1_real(file_name, args):
    d = open_file(file_name)
    x = d['x']
    x_bp = d['x_hat_l1_wo']
    x_cbp = d['x_hat_l1']
    x_rw = d['x_hat_rew']
    cl = d['cl']
    
    fig = plt.figure()
    ax = plt.subplot(111)

    alpha = 1

    plt.plot(x_bp, c='0.2', lw=5, alpha=alpha, label='BP')
    plt.plot(x_cbp, c='0.4', lw=5, alpha=alpha, label='BPCC')
    plt.plot(x_rw, c='0.6', lw=5, alpha=alpha, label=r'R$\ell_1$CC')
    plt.plot(x, '--',  c='0.0', lw=1.5, label=r'$x$')
    # We set the fontsize explicitly here because we use the figure in a two
    # figures per colum figure.
    plt.legend(loc='upper right',
               prop=mpl.font_manager.FontProperties(size=20))
    plt.ylim((-1.1, 1.1))
    plt.axhline(cl, c='k', ls=':')
    plt.axhline(-cl, c='k', ls=':')
    ax.set_yticks((-cl,cl))
    ax.set_yticklabels(('$C_l$','$C_u$'), fontsize=30)
    ax.set_xticks([])
    if args.footnote:
        footnote(file_name)
    #plt.gcf().subplots_adjust(bottom=bottom_margin)
    plt.gcf().subplots_adjust(bottom=0.01, top=0.99)

    #plt.show()
    fn = save_fig(fig, file_name)
    return fn

def test_1_real_color(file_name, args):
    d = open_file(file_name)
    x = d['x']
    x_bp = d['x_hat_l1_wo']
    x_cbp = d['x_hat_l1']
    x_rw = d['x_hat_rew']
    cl = d['cl']
    
    fig = plt.figure()
    ax = plt.subplot(111)

    plt.plot(x, c='0.0', lw=1, label='x')
    plt.plot(x_bp, c='r', lw=4, alpha=0.4, label='BP')
    #plt.plot(x_cbp, c='0.4', lw=4, alpha=0.4, label='C-BP')
    plt.plot(x_rw, c='b', lw=4, alpha=0.4, label='Rew')
    plt.legend(loc='upper right')
    plt.ylim((-1.1, 1.1))
    plt.axhline(cl, c='k', ls=':')
    plt.axhline(-cl, c='k', ls=':')
    ax.set_yticks((-cl,cl))
    ax.set_yticklabels(('$C_l$','$C_u$'))
    ax.set_xticks([])
    if args.footnote:
        footnote(file_name)
    #plt.show()
    plt.gcf().subplots_adjust(bottom=bottom_margin)
    fn = save_fig(fig, 'fig_1_b')
    return fn
    
def realization_sparse_signals(file_name, args):
    d = open_file(file_name)
    xs_det = d['xs_det']
    xs_rand = d['xs_rand']

    n = 1 # Number of realizations to plot
    
    fig = plt.figure(figsize=(12,3))
    axs_1 = [fig.add_subplot(n,2,2*i+1) for i in range(n)]
    axs_2 = [fig.add_subplot(n,2,2*i+2) for i in range(n)]
    for  i in range(n):
        axs_1[i].plot(xs_det[:,i], lw=2, c='k')
        axs_2[i].plot(xs_rand[:,i], lw=2, c='k')
    for ax in axs_1 + axs_2:
        ax.set_xticks([])
        ax.set_yticks([-5, 0, 5])
        ax.set_xlim([0, xs_det.shape[0]])
    axs_1[0].set_title('Deterministic support')
    axs_2[0].set_title('Random support')
    if args.footnote:
        footnote(file_name)
    #plt.show()
    plt.gcf().subplots_adjust(bottom=bottom_margin)
    fn = save_fig(fig, file_name)
    return fn

def latex_table(file_names, args):
    d1 = open_file(file_names[0])
    d2 = open_file(file_names[1])
    K1 = d1['K'].flatten()
    K2 = d2['K'].flatten()
    m1 = [m[0] for m in d1['methods'].flatten()]
    m2 = [m[0] for m in d2['methods'].flatten()]
    Ms = d1['Ms']
    Cls = d1['Cls']
    maxs = d1['maxs']
    p_recovery = d2['p_recovery']
    
    if set(K1) != set(K2):
        print "K1 and K2 are not equal. Don't know what to do"
        return locals()
    if m1 != m2:
        print "m1 and m2 are not equal. Processing only the common methods."
    K = K1
    methods = set(m1).intersection(m2)
    n_K = len(K)
    
    table = [r'\begin{tabular}{%s}' % ('c' * (n_K + 2),) ]
    table.append(r'\toprule')
    table.append(r'\multirow{2}{*}{Method}  & & \multicolumn{%d}{c}{Sparsity K} \\' %
                 (n_K,) ) 
    table.append('& ' + ' '.join([' & ' + str(k)  for k in K]) + r' \\ \hline')

    for m in methods:
        i_m1 = m1.index(m)
        i_m2 = m2.index(m)
        s1 = r'\multirow{3}{*}{%s} & $M_{min}$' % method_name[m]
        s2 = ' '.join([' & ' + str(k)  for k in Ms[:,i_m1]]) + r' \\ '
        s3 = r'\cline{3-%d}' % (n_K + 2, ) 
        table.append(s1 + s2 + s3)
        x = Cls[:,i_m1] / maxs.T
        x = x.flatten()
        s1 = '& $Cu_{min}$ '  + \
             ' '.join([' & ' + ('%.2f' % k)  for k in x]) + r' \\ '
        table.append(s1 + s3)
        s1 = '& p  '
        s2 = ' '.join([' & ' + '%.2f' % k  for k in p_recovery[i_m2,:]]) + r' \\ '
        table.append(s1 + s2 + r'\hline')

    # Remove the last \hline
    table[-1] = table[-1][:table[-1].find('\\hline')]
    table.append(r'\bottomrule')
    table.append(r'\end{tabular}')

    f = open('../reports/tech_report/table.tex','w')
    f.write('\n'.join(table))
    f.close()
    return locals()

def interpol_declipp(file_name, args):
    d = open_file(file_name)
    ss = d['ss'].flatten()
    y = d['y'].flatten()
    y_hats = d['y_hats']

    figs = []
    for i in range(len(ss)):
        figs.append(plt.figure())
        ax = plt.subplot(111)
        plt.plot(y, c='0.0', label='$x$')
        plt.plot(y_hats[i,:], c='0.2', lw=3, alpha=0.7, label='$\hat{x}$')
        plt.xlim((500, 750))
        plt.ylim((0,1.25))
        cl = ss[i]
        plt.axhline(cl, c='k', ls=':')
        ax.set_yticks((cl,))
        ax.set_yticklabels(('$C_u$',))
        ax.set_xticks([])
        plt.legend()
        if args.footnote:
            footnote(file_name)
    
    #plt.show()
    fn = save_fig(figs, file_name)
    return fn

def supp_estimation(file_name, args):
    d = open_file(file_name)
    alpha = np.abs(d['alpha'])
    alpha_c = np.abs(d['alpha_c'])
    x = d['x']
    x_c = d['x_c']
    n = np.double(d['n'].flatten())
    cl = d['cl']

    figs = []
    figs.append(plt.figure(figsize=(6,3)))
    ax = plt.subplot(111)
    plt.plot(x, c='0.4', lw=2, label='$x$')
    plt.plot(x_c, c='0.0', lw=1, alpha=1, label= '$x_c$')
    plt.axhline(cl, c='k', ls=':')
    plt.axhline(-cl, c='k', ls=':')
    ax.set_yticks((-cl,cl))
    # We set the fontsize explicitly here because we use the figure in a two
    # figures per colum figure.
    ax.set_yticklabels(('$-C_u$','$C_u$'), fontsize=22)
    ax.set_xticks([])
    plt.xlabel('Sample (n)')
    plt.gcf().subplots_adjust(bottom=bottom_margin, top=0.99)
    plt.legend(prop=mpl.font_manager.FontProperties(size=16))
    
    figs.append(plt.figure(figsize=(6,3)))
    ax = plt.subplot(111)
    n2 = len(n) / 2
    n = n[:n2]
    markerline, stemlines, baseline =  plt.stem(n, alpha[:n2], '0.5')
    plt.setp(markerline, 'markerfacecolor', '0.5')
    markerline, stemlines, baseline = plt.stem(n, alpha_c[:n2], '0.0')
    plt.setp(markerline, 'markerfacecolor', '0.0')
    plt.ylabel(r'$|\alpha(k)|$')
    #plt.legend(('a','', 'c'))
    ax.set_xticks([])
    ax.set_yticks([])
    plt.xlabel('Frequency (k)')    
    plt.gcf().subplots_adjust(bottom=bottom_margin, top=0.99)
    if args.footnote:
        footnote(file_name)
    fn = save_fig(figs, file_name)
    return fn

def test_8(file_name, args):
    d = open_file(file_name)
    x = d['x']
    x_c = d['x_c']
    x_hat = d['x_hat']
    cl = d['clip_level'].flatten()[0]
    
    figs = []
    figs.append(plt.figure(figsize=(6,3)))
    #figs.append(plt.figure(figsize=(12,6)))
    ax = plt.subplot(111)
    plt.plot(x, c='0.4', lw=4, label='$x$')
    plt.plot(x_hat, c='0.0', lw=0.5, alpha=1, label='$\hat{x}$')
    plt.legend()
    cl2 = 0.2
    plt.axhline(cl, c='k', ls=':')
    plt.axhline(-cl, c='k', ls=':')
    plt.axhline(cl2, c='k', ls=':')
    plt.axhline(-cl2, c='k', ls=':')
    ax.set_yticks((-cl, -cl2, cl2, cl))
    ax.set_yticklabels((r'$%.1f$' % -cl,
                        r'$%.1f$' % -cl2,
                        r'$%.1f$' % cl2,
                        r'$%.1f$' % cl))
    ax.set_xticks([])
    ax.set_xlim((0,128))
    plt.xlabel('$n$')
    if args.footnote:
        footnote(file_name)
    plt.gcf().subplots_adjust(bottom=bottom_margin)
    fn = save_fig(figs, file_name)
    return fn

def reconstruct_new(file_name, args):
    d = open_file(file_name)
    y_c = d['y_c']
    y_hat = d['y_hat']
    cu = d['clipu'][0,0]
    cl = d['clipl'][0,0]
    
    fig = plt.figure(figsize=(6,3))
    ax = plt.subplot(111)
    plt.plot(y_c, c='0.0', lw=0.5, label='$x_c$')
    plt.plot(y_hat, c='0.0', lw=1.5, alpha=0.5, label='$\hat{x}$')
    plt.legend()
    plt.axhline(cu, c='k', ls=':')
    plt.axhline(cl, c='k', ls=':')
    ax.set_yticks((cl,cu))
    ax.set_yticklabels(('$C_l$','$C_u$'))
    ax.set_xticks([])
    #ax.set_xlim((0,128))
    plt.xlabel('$n$')
    if args.footnote:
        footnote(file_name)
    plt.gcf().subplots_adjust(bottom=bottom_margin)
    fn = save_fig(fig, file_name)
    
    return fn

def tpcc(file_name, args):
    d = open_file(file_name)
    x = d['x']
    x_hat = d['x_hat']
    x_c = d['x_c']
    cl = d['clip_level']

    print d['Phi'].shape
    
    figs = []
    figs.append(plt.figure(figsize=(6,3)))
    ax = plt.subplot(111)
    plt.plot(x, c='0.0', lw=0.5, label='$x+z$')
    plt.plot(x_hat, c='0.0', lw=2, alpha=0.5, label='$\widehat{x}$')
    #plt.legend(loc=2)
    plt.axhline(cl, c='k', ls=':')
    plt.axhline(-cl, c='k', ls=':')
    ax.set_yticks((-cl,cl))
    ax.set_yticklabels(('$-C_u$','$C_u$'))
    ax.set_xticks([])
    ax.set_xlim((0,128))
    #ax.set_ylim((-4,4))
    plt.xlabel('$n$')

    n = 104
    xy_ann = (n, x[n]+0.12)
    ax.annotate(r'$x+z$', xy=xy_ann,  xycoords='data',
                xytext=(-40, -30), textcoords='offset points',
                arrowprops=dict(arrowstyle="->",
                                connectionstyle=s_d,))
    n = 55
    xy_ann = (n, x_hat[n]+0.06)
    ax.annotate(r'$\widehat{x}$', xy=xy_ann,  xycoords='data',
                xytext=(-29, -30), textcoords='offset points',
                arrowprops=dict(arrowstyle="->",
                                connectionstyle=s_d,))
    if args.footnote:
        footnote(file_name)
    plt.gcf().subplots_adjust(bottom=bottom_margin)
    fn = save_fig(figs, file_name)
    
    #plt.show()
    return fn

def dispatcher(file_name, args):
    if file_name.count('experiment_1'):
        d = experiment_1(file_name, args)
    elif file_name.count('experiment_2'):
        d = experiment_2(file_name, args)
    elif file_name.count('experiment_3'):
        d = experiment_3(file_name, args)
    elif file_name.count('mangasarian'):
        d = mangasarian(file_name, args)
    elif file_name.count('test_1_real'):
        d = test_1_real(file_name, args)
        #d = test_1_real_color(file_name)
    elif file_name.count('realization_sparse_signals'):
        d = realization_sparse_signals(file_name, args)
    elif file_name.count('interpol_declipp'):
        d = interpol_declipp(file_name, args)
    elif file_name.count('supp_estimation'):
        d = supp_estimation(file_name, args)
    elif file_name.count('test_8'):
        d = test_8(file_name, args)
    elif file_name.count('reconstruct_new'):
        d = reconstruct_new(file_name, args)
    elif file_name.count('pOMP'):
        d = tpcc(file_name, args)
    else:
        print "Don't know what to do"
    return d

    

if __name__ == '__main__':
    # Dictionary with 'key = file name of the PDF used in the latex document'
    # and 'value = file name of the corresponding mat file' (both without the 
    # file extension)
    plot_list = {'exp1_TPCC': 'experiment_1_20111014T223757',
                 'exp2_TPCC': 'experiment_2_20111015T220756',
                 'exp3_TPCC': 'experiment_3_20111016T102024',
                 'supp_estimation': 'supp_estimation_20101217T132559',
                 'test_1_a': 'test_1_real_20101122T093959',
                 'test_1_b': 'test_1_real_20101122T100157',
                 'test_8': 'test_8_20101221T090224'}

    matfiles_dir = '../matlab'
    latex_dir = '../reports/figures/'

    parser = argparse.ArgumentParser(description='Generate PDF files with the '
                                     'plots for the declipping paper.')
    parser.add_argument('-footnote', action='store_true', default=False,
                        help='Add footprint with the file name to the figures.')
    parser.add_argument('-latex', action='store_true', default=False,
                        help='Compile the latex document')
    parser.add_argument('-all', action='store_true', default=False,
                        help='Create all the figures.')
    parser.add_argument('-print_list', action='store_true', default=False,
                        help='Print the list of plots and their filenames')
    parser.add_argument('files', nargs='*')
    args = parser.parse_args()

    if args.print_list:
        for latex_name, file_name in plot_list.iteritems():
            print '%s -> %s' % (latex_name, file_name)
        sys.exit(0)
        
    if len(args.files) == 0 and args.all is False:
        print 'Error: Specify at least one file or use -all option.'
        sys.exit(0)

    for name in args.files:
        if name not in plot_list:
            print 'Error: %s not in the list of plots.' % name
            sys.exit(0)

    for latex_name, file_name in plot_list.iteritems():
        if latex_name in args.files or args.all:
            fn = os.path.join(matfiles_dir, file_name + '.mat')
            print '----- Creating plot for file', fn, '-----'
            pdf_fn = dispatcher(fn, args)
            latex_pdf_name = os.path.join(latex_dir, latex_name + '.pdf')
            print 'Copying %s to %s' % (pdf_fn, latex_pdf_name)
            shutil.copy(pdf_fn, latex_pdf_name)

    if args.latex:
        doc_dir = '../reports/paper'
        latex_file = 'declipping_paper.tex'
        print 'Compiling', latex_file, '...'
        devnull = open(os.devnull, 'w')
        subprocess.call(('pdflatex', latex_file), cwd=doc_dir, stdout=devnull)
    
    print 'Done!'
