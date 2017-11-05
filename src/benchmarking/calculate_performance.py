import pandas as pd
import numpy as np
from sklearn import metrics
import matplotlib.pyplot as plt



def calculate_perfromance(folder_path):
    sumauc = []
    conlist = ['no_constraints', 'GC_constraints', 'aa_constraints', 'length_constraints',\
               'GCAA_constraints', 'AAlength_constraints', 'GClength_constraints', \
               'ecoli_constraints', 'bsubtilus_constraints', 'scoelicolor_constraints', \
               'combined']

    for constraint in conlist:
        if constraint == 'combined':
            df_ecoli = pd.read_csv('%secoli_constraints.csv' %folder_path, index_col=0)
            df_bsub = pd.read_csv('%sbsubtilus_constraints.csv' %folder_path, index_col=0)
            df_scoeli = pd.read_csv('%sscoelicolor_constraints.csv' %folder_path, index_col=0)

        else:
            df = pd.read_csv('%s%s.csv'%(folder_path, constraint), index_col=0)

        TPR = []
        FPR = []
        hm = []

        for i in range(20, 62):
            hm.append([])
            for j in range(20, 62):
                if j > i:
                    if constraint == 'combined':
                        highlist = list(df_ecoli[df_ecoli.CUB == j]['calc_CUB']) + \
                                    list(df_bsub[df_bsub.CUB == j]['calc_CUB']) + \
                                    list(df_scoeli[df_scoeli.CUB == j]['calc_CUB'])
                        lowlist = list(df_ecoli[df_ecoli.CUB == i]['calc_CUB']) + \
                                    list(df_bsub[df_bsub.CUB == i]['calc_CUB']) + \
                                    list(df_scoeli[df_scoeli.CUB == i]['calc_CUB'])
                    else:
                        highlist = list(df[df.CUB == j]['calc_CUB'])
                        lowlist = list(df[df.CUB == i]['calc_CUB'])

                    y = np.array([0]*len(highlist) + [1]*len(lowlist))
                    scores = np.array(highlist + lowlist)
                    fpr, tpr, thresholds = metrics.roc_curve(y, scores, pos_label=0)
                    auc = metrics.auc(fpr, tpr, reorder=True)
                    if auc < 0.5:
                        auc = 1-auc

                    hm[-1].append(auc)
                else:
                    hm[-1].append(1)
        hm = np.array(hm).transpose()
        avgAUC = (sum(sum(hm))-903)/861
        sumauc.append(avgAUC)
    return np.array(sumauc)*2-1

def compare_performance(newlist):

    barray_short = pd.read_csv('./data/metric_performance.csv', index_col=0)

    MetricList = ['Info', 'Nov', 'RCBS', 'SCUO', 'Wright']
    Metriclist = ['iCUB', r'$N_C^{\prime}$', 'RCBS', 'SCUO', r'$N_C$']
    critlist = ['EqualCodon', \
                'GCVariation_ecoli','AAVariaton_ecoli','LengthVariation_ecoli',\
                'GCAAVariation_ecoli', 'AALengthVariaton_ecoli', \
                'GCLengthVariation_ecoli', \
                'Xseq_ecoli', 'Xseq_bacsub', \
                'Xseq_scoeli','Xseq_comb']
    critlist1 = ['Simple', 'GC', 'AA', 'Length', \
                 'AA+GC', 'AA+Length', 'GC+Length', \
                'E. coli', 'B. sub.', \
                'S. coeli.', 'E. coli + \nB. sub. + \nS. coeli']
    cmap = plt.cm.seismic
    colorlist = ['c'] + [cmap(i) for i in np.linspace(0,1,6)][0:4]
    linewidthlist = [1,1,1,1,1]

    criteraname = dict(zip(critlist, critlist1))
    mnameconverstion = dict(zip(MetricList, Metriclist))
    colordict = dict(zip(MetricList, colorlist))
    linewidthdict = dict(zip(MetricList, linewidthlist))

    fig = plt.figure(figsize=(8, 3))
    temp = fig.add_axes([0,0,1,1])
    temp.set_yticks([])
    temp.set_xticks([])
    temp.spines['right'].set_color('none')
    temp.spines['top'].set_color('none')
    temp.spines['bottom'].set_color('none')
    temp.spines['left'].set_color('none')

    xind = 0.09
    width = 1-xind
    ayind = 0.06
    aheight = 1-ayind
    wid = 0.4
    # ind = np.arange(len(barray_short))
    x = np.arange(1,12)
    sub1 = fig.add_axes([xind, ayind, width, aheight])
    sub1.yaxis.grid(color='1', linestyle='-', zorder=3)
    sub1.xaxis.grid(color='1', linestyle='-', zorder=3)
    b=0.5
    for M in MetricList:
        sub1.plot(x, barray_short[M], label=mnameconverstion[M], color=colordict[M],
                  marker='.',linewidth=linewidthdict[M], zorder=13)
    sub1.plot(x, newlist, label='New', color='red', marker='.', linewidth=2)
    sub1.fill_between([0.5, 1.5], [b, b], [1, 1], color='0.9')
    sub1.fill_between([1.5, 4.5], [b, b], [1, 1], color='0.8')
    sub1.fill_between([4.5, 7.5], [b, b], [1, 1], color='0.7')
    sub1.fill_between([7.5, 10.5], [b, b], [1, 1], color='0.6')
    sub1.fill_between([10.5, 11.5], [b, b], [1, 1], color='0.5')
    sub1.plot([1.5, 1.5], [b, 1.3], linestyle='--', color='black', linewidth=0.5)
    sub1.plot([4.5, 4.5], [b, 1.3], linestyle='--', color='black', linewidth=0.5)
    sub1.plot([7.5, 7.5], [b, 1.3], linestyle='--', color='black', linewidth=0.5)
    sub1.plot([10.5, 10.5], [b, 1], linestyle='--', color='black', linewidth=0.5)
    sub1.set_xticks(np.arange(1,12))
    sub1.set_xticklabels(critlist1, rotation=20, ha='right', rotation_mode='anchor',
                    va='top')
    sub1.set_ylim(b, 1)
    sub1.set_xlim(0.5, 11.5)
    sub1.set_ylabel('Average Score')
    sub1.tick_params(length=0)
    y = b
    sub1.text(3, y, '1 constraint', va='bottom', ha='center')
    sub1.text(6, y, '2 constraints', va='bottom', ha='center')
    sub1.text(9, y, 'Organismal Genomes \n (3 constraints)', va='bottom', ha='center')
    sides = ['right', 'top', 'left', 'bottom']
    for line in sides:
        sub1.spines[line].set_color('none')

    legend = sub1.legend(frameon=False, loc=9, bbox_to_anchor=(0.5, 1.17)
                        ,fontsize=12, ncol=6, columnspacing=3)

    plt.show()
