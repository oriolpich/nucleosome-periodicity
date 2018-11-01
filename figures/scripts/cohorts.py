import collections
import math

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from statsmodels.sandbox.stats.multicomp import multipletests

from nucperiod import plot as nuceriod_plt, increase


TTYPES = {'ALL-US':'Acute Lymphoblastic Leukemia',
    'AML-US':'Acute Myeloid Leukemia',
    'BLCA-CN':'Bladder Urothelial Cancer',
    'BLCA-US':'Bladder Urothelial Cancer',
    'BOCA-FR':'Bone Cancer',
    'BOCA-UK':'Bone Cancer',
    'BRCA-CN':'Breast Cancer',
    'BRCA-EU':'Breast Cancer',
    'BRCA-FR':'Breast Cancer',
    'BRCA-KR':'Breast Cancer',
    'BRCA-MX':'Breast Cancer',
    'BRCA-UK':'Breast Cancer',
    'BRCA-US':'Breast Cancer',
    'BTCA-JP':'Biliary Tract Cancer',
    'BTCA-SG':'Biliary Tract Cancer',
    'CCSK-US':'Clear Cell Sarcomas of the Kidney',
    'CESC-US':'Cervical Squamous Cell Carcinoma',
    'CLLE-ES':'Chronic Lymphocytic Leukemia',
    'CMDI-UK':'Chronic Myeloid Disorders',
    'COAD-US':'Colorectal Cancer',
    'CRC':'Colorectal Cancer',
    'COCA-CN':'Colorectal Cancer',
    'DLBC-US':'Lymphoid Neoplasm Diffuse Large B-cell Lymphoma',
    'EOPC-DE':'Prostate Adenocarcinoma',
    'ESAD-UK':'Esophageal Cancer',
    'ESCA-CN':'Esophageal Cancer',
    'GACA-CN':'Gastric Cancer',
    'GBM-CN':'Glioblastoma Multiforme',
    'GBM-US':'Glioblastoma Multiforme',
    'HNCA-MX':'Head and Neck Cancer',
    'HNSC-US':'Head and Neck Cancer',
    'KICH-US':'Kidney Chromophobe',
    'KIRC-US':'Kidney Clear Cell Carcinoma',
    'KIRP-US':'Kidney RPCC',
    'LAML-CN':'Acute Myeloid Leukemia',
    'LAML-KR':'Acute Myeloid Leukemia',
    'LAML-US':'Acute Myeloid Leukemia',
    'LGG-US':'Brain Lower Grade Glioma',
    'LIAD-FR':'Benign Liver Tumour',
    'LICA-CN':'Liver Cancer',
    'LICA-FR':'Liver Cancer',
    'LIHC-US':'Liver Cance',
    'LIHM-FR':'Liver Cancer',
    'LINC-JP':'Liver Cancer',
    'LIRI-JP':'Liver Cancer',
    'LMS-FR':'Leiomyosarcoma',
    'LUAD-US':'Lung Adenocarcinoma',
    'LUSC-CN':'Lung Squamous Cell Carcinoma',
    'LUSC-KR':'Lung Squamous Cell Carcinoma',
    'LUSC-US':'Lung Squamous Cell Carcinoma',
    'MALY-DE':'Malignant Lymphoma',
    'MELA-AU':'Skin Cancer',
    'NACA-CN':'Nasopharyngeal cancer',
    'NBL-US':'Neuroblastoma',
    'NHLY-MX':'Non Hodgkin Lymphoma',
    'NKTL-SG':'NK T Cell Lymphoma',
    'ORCA-IN':'Oral Cancer',
    'OS-US':'Osteosarcoma',
    'OV-AU':'Ovarian Cancer',
    'OV-CN':'Ovarian Cancer',
    'OV-US':'Ovarian Cancer',
    'PAAD-US':'Pancreatic Cancer',
    'PACA-AU':'Pancreatic Endocrine',
    'PACA-CA':'Pancreatic Cancer',
    'PACA-CN':'Pancreatic Ductal',
    'PAEN-AU':'Pancreatic Endocrine',
    'PAEN-IT':'Pancreatic Endocrine',
    'PBCA-DE':'Pediatric Brain Cancer',
    'PEME-CA':'Pediatric Medulloblastoma',
    'PRAD-CA':'Prostate Adenocarcinoma',
    'PRAD-FR':'Prostate Adenocarcinoma',
    'PRAD-CN':'Prostate Adenocarcinoma',
    'PRAD-UK':'Prostate Adenocarcinoma',
    'PRAD-US':'Prostate Adenocarcinoma',
    'PRCA-FR':'Prostate Adenocarcinoma',
    'READ-US':'Rectum Adenocarcinoma',
    'RECA-CN':'Kidney Clear Cell Carcinoma',
    'RECA-EU':'Kidney Clear Cell Carcinoma',
    'RT-US':'Rhabdoid Tumors',
    'RTBL-FR':'Retinoblastoma',
    'SARC-US':'Sarcoma',
    'SKCA-BR':'Skin Cancer',
    'SKCM-US':'Skin Cancer',
    'STAD-US':'Gastric Adenocarcinoma',
    'THCA-CN':'Thyroid Cancer',
    'THCA-SA':'Thyroid Cancer',
    'THCA-US':'Thyroid Cancer',
    'UTCA-FR':'Uterine Cancer',
    'UCEC-US':'Uterine Cancer',
    'WT-US':'Wilms Tumor',
    'Eyelid': 'Normal Skin'}
TTYPES.update({k.split('-')[0]: v for (k, v) in TTYPES.items()})


COLORS = {'Colorectal Cancer':'#191970',
    'Acute Myeloid Leukemia':'#CD6600',
    'Lung Adenocarcinoma':'#efd014ff',
    'Benign Liver Tumour':'grey',
    'Lung Squamous Cell Carcinoma':'#FDF5E6',
    'Kidney Clear Cell Carcinoma':'#FF4500',
    'Brain Lower Grade Glioma':'#D8BFD8',
    'Kidney Chromophobe':'#FF4500',
    'Uterine Cancer':'#FF8C69',
    'Pediatric Brain Cancer':'blue',
    'Glioblastoma Multiforme':'#3D3D3D',
    'Oral Cancer':'yellow',
    'Head and Neck Cancer':'#8B2323',
    'Bladder Urothelial Cancer':'#EEAD0E',
    'Liver Cancer':'#006400',
    'Kidney RCCC':'#FF4500',
    'Pancreatic Cancer':'#7A378B',
    'Ovarian Cancer':'#008B8B',
    'Pancreatic Endocrine':'#E066FF',
    'Prostate Adenocarcinoma':'#87CEFA',
    'Malignant Lymphoma':'#DDCDCD',
    'Chronic Myeloid Disorders':'#DDCDCD',
    'Chronic Lymphocytic Leukemia':'#F4A35D',
    'Breast Cancer':'#CD6090',
    'Leiomyosarcoma':'#FFEC8B',
    'Head and Neck Thyroid Carcinoma':'pink',
    'Bone Cancer':'#FFD700',
    'Skin Cancer':'#000000',
    'Biliary Tract Cancer':'#00CD66',
    'NK T Cell Lymphoma':'#698B22',
    'Gastric Cancer':'#BFEFFF',
    'Esophageal Cancer':'#1E90FF',
    'Thyroid Cancer':'brown',
    'Normal Skin': '#000000'}


def find(files, min_muts=None):
    df = increase.load_d(files)
    df['ttype'] = df['name'].map(TTYPES)

    if min_muts is not None:
        df = df[df['nmuts_whole_nucleosome'] > min_muts]

    df.sort_values(by='nmuts_whole_nucleosome', ascending=False, inplace=True)

    ttypes = set()
    names = []
    for i, row in df.iterrows():
        if row['ttype'] not in ttypes:
            names.append(row['name'])
            ttypes.add(row['ttype'])
    return names


def _load(cohorts, tumors):

    df = increase.load_d(cohorts)

    df['increase_in'] = df['observed_in'] - df['expected_in']
    df['prop_increase_in'] = df['increase_in'] / df['expected_in']
    df['ttype'] = df['name'].map(TTYPES)

    df = df[df['name'].isin(tumors)].sort_values(by='prop_increase_in', ascending=True)

    # add Q-value
    df['empirical_pvalue_snr'] = df['empirical_pvalue_snr'].fillna(1)
    df['empirical_pvalue_snr'].replace(0, 0.001, inplace=True)

    qvals = multipletests(df['empirical_pvalue_snr'].tolist(), method='fdr_bh')
    df['qvals_snr'] = qvals[1]

    return df


def zoomout(cohorts, tumors):
    nuceriod_plt.config_params(14)

    toplot = _load(cohorts, tumors)

    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(8, 7), sharex=True)

    for i, row in toplot.iterrows():
        edgecolor = None
        snr = row['snr']
        if row['qvals_snr'] < 0.05 and snr > 8:
            edgecolor = 'darkred'
        if row['cross_validation_max'] > 0:

            ax[0].scatter(row['peak'], np.log2(snr), s=200, linewidth=2, edgecolor=edgecolor, color='white',
                          label=row['name'], alpha=1)

            ax[0].scatter(row['peak'], np.log2(snr), s=80, edgecolor='grey', linewidth=0.5, color=COLORS[row['ttype']],
                          label=row['name'], alpha=1)

            if edgecolor == 'darkred':
                ax[0].text(row['peak'] + 10, np.log2(snr) - 0.15, row['ttype'])

    for i, row in toplot.iterrows():
        edgecolor = None
        snr = row['snr']
        if row['qvals_snr'] < 0.05 and snr > 8:
            edgecolor = 'darkred'
        if row['cross_validation_max'] < 0:

            ax[1].scatter(row['peak'], -np.log2(snr), s=200, linewidth=2, edgecolor=edgecolor, color='white',
                          label=row['name'], alpha=1)

            ax[1].scatter(row['peak'], -np.log2(snr), s=80, edgecolor='grey', linewidth=0.5, color=COLORS[row['ttype']],
                          label=row['name'], alpha=1)

            if edgecolor == 'darkred':
                ax[1].text(row['peak'] + 10, -np.log2(snr) - 0.15, row['ttype'])

    plt.ylabel('log2(SNR)')
    plt.xlabel('Period')

    xlim = [i for i in range(50, 260, 20)]
    ax[0].set_xticks(xlim)
    ax[1].set_xticks(xlim)
    yvals = [i for i in range(2, 10, 2)]
    ax[0].set_yticks(yvals)

    yvals = [i for i in range(-8, 0, 2)]
    ax[1].set_yticks(yvals)
    ylabels = [str(2 ** abs(i)) for i in range(2, 10, 2)]
    ax[0].set_yticklabels(ylabels)

    ylabels = [str(2 ** abs(i)) for i in range(-8, 0, 2)]
    ax[1].set_yticklabels(ylabels)
    ax[0].spines['right'].set_visible(False)
    ax[0].spines['top'].set_visible(False)

    ax[1].xaxis.set_ticks_position('top')
    ax[1].spines['bottom'].set_visible(False)
    ax[1].spines['right'].set_visible(False)

    ax[0].set_ylim(1.5, 10)
    ax[1].set_ylim(-10, -1.5)

    plt.tight_layout()


def zoomin(cohorts, tumors):
    nuceriod_plt.config_params(14)

    toplot = _load(cohorts, tumors)

    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(8, 7), sharex=True)

    toplot.sort_values(by='qvals_snr', ascending=False, inplace=True)

    for i, row in toplot.iterrows():
        edgecolor = None
        snr = row['snr']
        if (row['qvals_snr'] < 0.05)&(snr>8):
            edgecolor = 'darkred'

        if row['cross_validation_max'] < 0:

            ax[1].scatter(row['peak'], -np.log2(snr), s=200, linewidth=2, edgecolor=edgecolor, color='white',
                          label=row['name'], alpha=1)

            ax[1].scatter(row['peak'], -np.log2(snr), s=80, edgecolor='grey', linewidth=0.5, color=COLORS[row['ttype']],
                          label=row['name'], alpha=1)

            if edgecolor == 'darkred':
                ax[1].text(row['peak'] + 1, -np.log2(snr) - 0.15, row['ttype'])

    for i, row in toplot.iterrows():
        edgecolor = None
        snr = row['snr']
        if row['qvals_snr'] < 0.05 and snr > 8:
            edgecolor = 'darkred'
        if row['cross_validation_max'] > 0:
            ax[0].scatter(row['peak'], np.log2(snr), s=200, linewidth=2, edgecolor=edgecolor, color='white',
                          label=row['name'], alpha=1)

            ax[0].scatter(row['peak'], np.log2(snr), s=80, edgecolor='grey', linewidth=0.5, color=COLORS[row['ttype']],
                          label=row['name'], alpha=1)

            if edgecolor == 'darkred':
                ax[0].text(row['peak'] + 1, np.log2(snr) - 0.15, row['ttype'])

    plt.xlabel('Period')

    xlim = [i for i in range(8, 22, 2)]
    ax[0].set_xticks(xlim)
    ax[1].set_xticks(xlim)
    yvals = [i for i in range(2, 10, 2)]
    ax[0].set_yticks(yvals)

    yvals = [i for i in range(-8, 0, 2)]
    ax[1].set_yticks(yvals)
    ylabels = [str(2 ** abs(i)) for i in range(2, 10, 2)]
    ax[0].set_yticklabels(ylabels)

    ylabels = ['{}'.format(str(2 ** abs(i))) for i in range(-8, 0, 2)]
    ax[1].set_yticklabels(ylabels)

    ax[0].spines['right'].set_visible(False)
    ax[0].spines['top'].set_visible(False)

    ax[1].xaxis.set_ticks_position('top')
    ax[1].spines['bottom'].set_visible(False)
    ax[1].spines['right'].set_visible(False)

    ax[0].set_ylim(1.5, 10)
    ax[1].set_ylim(-10, -1.5)

    plt.tight_layout()


def compare(cohorts_5mer, cohorts_3mer, cohorts_linker, tumors=None):
    nuceriod_plt.config_params(14)

    df_5mer = increase.load_d(cohorts_5mer)
    df_5mer['control'] = 'mer5'
    df_3mer = increase.load_d(cohorts_3mer)
    df_3mer['control'] = 'mer3'
    df_linker = increase.load_d(cohorts_linker)
    df_linker['control'] = 'linker'

    df = pd.concat([df_5mer, df_3mer, df_linker])

    df['increase_in'] = df['observed_in'] - df['expected_in']
    df['prop_increase_in'] = df['increase_in'] / df['expected_in']

    df['ttype'] = df['name'].map(TTYPES)
    if tumors is not None:
        df = df[df['name'].isin(tumors)]

    toplot = df.sort_values(by='prop_increase_in', ascending=True)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 5))

    labels = []
    colors = []
    red_patch = mpatches.Patch(color='red', label='3-mers')
    green_patch = mpatches.Patch(color='green', label='no nucleosomes in context')
    orange_patch = mpatches.Patch(color='orange', label='5-mers')

    for ix, (ttype, data) in enumerate(toplot.sort_values(by='snr', ascending=True).groupby(by='ttype', sort=False)):
        snr1 = data[data['control'] == 'mer3']['snr'].tolist()[0]
        snr2 = data[data['control'] == 'linker']['snr'].tolist()[0]
        snr3 = data[data['control'] == 'mer5']['snr'].tolist()[0]
        colors.append('red')
        colors.append('green')
        colors.append('orange')
        ax.scatter(ix, math.log2(snr1), color='red', s=15, alpha=0.8)
        ax.scatter(ix, math.log2(snr2), color='green', s=15, alpha=0.8)
        ax.scatter(ix, math.log2(snr3), color='orange', s=15, alpha=0.8)

        labels.append(ttype)

    plt.xticks([i for i in range(ix + 1)], labels, rotation=90)
    tick = [2, 4, 6, 8]
    plt.yticks(tick, [str(2 ** t) for t in tick])
    plt.ylabel('log2(SNR)')
    plt.legend(handles=[red_patch, green_patch, orange_patch])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()


def rotational(cohorts_high, cohorts_low, tumors=None):
    nuceriod_plt.config_params(14)

    df_high = increase.load_d(cohorts_high)
    df_high['control'] = 'high'
    df_low = increase.load_d(cohorts_low)
    df_low['control'] = 'low'

    df = pd.concat([df_high, df_low])

    df['increase_in'] = df['observed_in'] - df['expected_in']
    df['prop_increase_in'] = df['increase_in'] / df['expected_in']

    df['ttype'] = df['name'].map(TTYPES)
    if tumors is not None:
        df = df[df['name'].isin(tumors)]

    toplot = df.sort_values(by='prop_increase_in', ascending=True)

    order = ['low', 'high']
    count = 0
    xvals = []
    yvals = []
    colors = []
    dic_t = collections.defaultdict(dict)
    labels = []
    for sig, data in toplot.sort_values(by='snr').groupby(by='ttype', sort=False):
        if sig in COLORS:
            labels.append('{}'.format(sig))
            for i in order:
                val = data[data['control'] == i]['snr'].tolist()[0]
                dic_t['Sign {}'.format(sig)][i] = val
                yvals.append(math.log2(val))
                xvals.append(count)
                colors.append(COLORS[sig])
                count += 1

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 1.5))
    ax.bar(xvals, yvals, color=colors, label=labels)
    ax.set_ylabel('SNR')
    plt.xticks(np.arange(0.5, 54, 2), labels, rotation=90, fontsize=13)
    tick = [2, 4, 6, 8]
    plt.yticks(tick, [str(2 ** t) for t in tick])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def _get_project(name):
    if '-' in name:
        return name.split('-')[1]
    else:
        return '505'


def generate_table(cohorts, tumors):
    df = _load(cohorts, tumors)
    df['project'] = df['name'].apply(_get_project)
    table_out = df[['ttype', 'name', 'project', 'peak', 'snr', 'empirical_pvalue_snr', 'qvals_snr', 'cross_validation_max']].copy()
    table_out.columns = ['Tumor Name', 'Cohort', 'Project', 'Peak', 'SNR', 'P-value', 'Q-value', 'Phase']
    table_out.sort_values(by='Tumor Name', inplace=True)
    return table_out
