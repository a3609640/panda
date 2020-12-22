#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 14:48:10 2020

@author: suwu
"""
from collections import Iterable
import igraph
import markov_clustering as mc
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles
import mygene
import networkx as nx
import numpy as np
import os
import pandas as pd
from pandas import DataFrame
import pygraphviz as pgv
from networkx.drawing.nx_agraph import graphviz_layout
import pylab
import seaborn as sbn
from scipy.cluster.hierarchy import dendrogram, linkage
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as ssd
from scipy.stats import pearsonr
from sklearn.cluster import KMeans 
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import random

# set current direcotry
os.chdir('/home/suwu/Documents/Protein_compl_Output')
os.getcwd()

###### retrieve proteomics data ######

def ccle_pro():
    ## download the proteomics data 
    ## https://gygi.hms.harvard.edu/data/ccle/protein_quant_current_normalized.csv.gz
    CCLE_PRO = pd.read_csv("~/Downloads/protein_quant_current_normalized.csv") 
    # concatenate the Gene_Symbol (with duplicate names) and the Uniprot_Acc
    CCLE_PRO["Gene_Symbol"] = CCLE_PRO["Gene_Symbol"].fillna('') + ' ' + CCLE_PRO["Uniprot_Acc"]
    CCLE_PRO = CCLE_PRO.set_index('Gene_Symbol')
    CCLE_PRO_subset = CCLE_PRO[CCLE_PRO.columns.drop(list(CCLE_PRO.filter(regex='Peptides')))]
    CCLE_PRO_subset = CCLE_PRO_subset.drop(columns=['Protein_Id', 
                                                    'Description',
                                                    'Group_ID',
                                                    'Uniprot',
                                                    'Uniprot_Acc']).T
    #CCLE_PRO_subset.index = CCLE_PRO_subset.index.str.split('_').str[0]
    CCLE_PRO_subset.index.names = ['ccle_name']
    return CCLE_PRO_subset
CCLE_PRO_subset = ccle_pro()
#CCLE_PRO_subset[CCLE_PRO_subset.duplicated()]


EIF =  pd.Index(['EIF4G1 Q04637-9',
                 'EIF4A1 P60842',
                 'EIF4E P06730-2',
                 'EIF4EBP1 Q13541'])

def reg_coef(x,y,label=None,color=None,**kwargs):
    ax = plt.gca()
    r,p = pearsonr(x,y)
    ax.annotate('r = {:.2f}'.format(r), xy=(0.5,0.5), xycoords='axes fraction', ha='center')
    ax.set_axis_off()
    
def eif_ccle_scatter(data, proteins):
    CCLE_EIF_PRO = CCLE_PRO_subset[proteins]
    CCLE_EIF_PRO['ccle'] = CCLE_EIF_PRO.index
    CCLE_EIF_PRO['ccle'] = CCLE_EIF_PRO['ccle'].str.split("_").str[1] 
    #print(CCLE_EIF_PRO.dtypes)
    #scatter plot of two protein expression across cell lines
    CCLE_EIF_PRO.dropna(inplace = True)
    g = sbn.PairGrid(CCLE_EIF_PRO)
    #g.map(sbn.scatterplot)
    g.map_diag(sbn.histplot)
    #g.map_lower(sbn.scatterplot)
    g.map_lower(sbn.regplot)
    g.map_upper(reg_coef)
    g.add_legend()
    #scatter plot of two protein expression across cell lines with color
    g = sbn.PairGrid(CCLE_EIF_PRO, hue="ccle",corner=True)
    g.map_lower(sbn.scatterplot)
    g.map_diag(sbn.histplot)
    g.add_legend()
eif_ccle_scatter(CCLE_PRO_subset, EIF)    


###### retrieve proteomics data of EIF4F complex ######
def eif_corr_pro(eif):
    y = pd.DataFrame() 
    for gene in eif:
        print(gene)
        x =  CCLE_PRO_subset.corrwith(CCLE_PRO_subset[gene])
        x.name = gene
       # x.index.is_unique
       # x.index = x.index.where(~x.index.duplicated(), x.index + '_dp')
        y = pd.concat([y, x], axis=1)
    y.dropna(inplace=True)    
    return y

EIF_COR_PRO = eif_corr_pro(EIF)

EIF_COR_PRO_sig = EIF_COR_PRO.loc[(EIF_COR_PRO['EIF4G1 Q04637-9'] >= 0.5) | (EIF_COR_PRO['EIF4G1 Q04637-9'] <= -0.5)
                                  #|(EIF_COR_PRO['EIF4G1 Q04637-8'] >= 0.3) |  (EIF_COR_PRO['EIF4G1 Q04637-8'] <= -0.3) 
                                  | (EIF_COR_PRO['EIF4A1 P60842'] >= 0.5) |  (EIF_COR_PRO['EIF4A1 P60842'] <= -0.5) 
                                  | (EIF_COR_PRO['EIF4E P06730-2'] >= 0.5) |  (EIF_COR_PRO['EIF4E P06730-2'] <= -0.5) 
                                  | (EIF_COR_PRO['EIF4EBP1 Q13541'] >= 0.5) |  (EIF_COR_PRO['EIF4EBP1 Q13541'] <= -0.5) 
                                  #| (EIF_COR_PRO['MKNK2 Q9HBH9'] >= 0.5) |  (EIF_COR_PRO['MKNK2 Q9HBH9'] <= -0.5) 
                                  #| (EIF_COR_PRO['MKNK1 Q9BUB5'] >= 0.5) |  (EIF_COR_PRO['MKNK1 Q9BUB5'] <= -0.5) 
                                  ] 

def eif_corr_scatter(data):
    data.dropna(inplace = True)
    g = sbn.PairGrid(data)
    #g.map(sbn.scatterplot)
    g.map_diag(sbn.histplot)
    #g.map_lower(sbn.scatterplot)
    g.map_lower(sbn.regplot)
    g.map_upper(reg_coef)
#CCLE_PRO[CCLE_PRO['Gene_Symbol'].str.contains("MYC")]['Gene_Symbol']
eif_corr_scatter(EIF_COR_PRO)
eif_corr_scatter(EIF_COR_PRO_sig)

def plot_EIF_Venn(data): 
    pos_eIF4G1 = pd.Series(data.index[data['EIF4G1 Q04637-9'] >= 0.5])
    pos_eIF4G1 = pos_eIF4G1.apply(lambda x: x.split(' ')[0])
    pos_eIF4A1 = pd.Series(data.index[data['EIF4A1 P60842'] >= 0.5])
    pos_eIF4A1 = pos_eIF4A1.apply(lambda x: x.split(' ')[0])
    pos_eIF4E = pd.Series(data.index[data['EIF4E P06730-2'] >= 0.5])
    pos_eIF4E = pos_eIF4E.apply(lambda x: x.split(' ')[0])
    fig = plt.figure(figsize = (12, 10))
    pos_A = set(pos_eIF4G1)
    pos_B = set(pos_eIF4A1)
    pos_C = set(pos_eIF4E)
    v = venn3([pos_A, pos_B, pos_C], ('pos_eIF4G1', 'pos_eIF4A1', 'pos_eIF4E'))
    
    neg_eIF4G1 = pd.Series(data.index[data['EIF4G1 Q04637-9'] <= -0.5])
    neg_eIF4G1 = neg_eIF4G1.apply(lambda x: x.split(' ')[0])
    neg_eIF4A1 = pd.Series(data.index[data['EIF4A1 P60842'] <= -0.5])
    neg_eIF4A1 = neg_eIF4A1.apply(lambda x: x.split(' ')[0])
    neg_eIF4E = pd.Series(data.index[data['EIF4E P06730-2'] <= -0.5])
    neg_eIF4E = neg_eIF4E.apply(lambda x: x.split(' ')[0])
    fig = plt.figure(figsize = (12, 10))
    neg_eIF4G1 = set(neg_eIF4G1)
    neg_eIF4A1 = set(neg_eIF4A1)
    neg_eIF4E = set(neg_eIF4E)
    v = venn3([neg_eIF4G1, neg_eIF4A1, neg_eIF4E], ('neg_eIF4G1', 'neg_eIF4A1', 'neg_eIF4E'))
plot_EIF_Venn(EIF_COR_PRO)


#EIF =  pd.Index(['EIF4G1 Q04637-9', 
#                 'EIF4A1 P60842',
#                 'EIF4E P06730-2',
#                 'EIF4EBP1 Q13541',
#                 'MKNK1 Q9BUB5',
#                 'MKNK2 Q9HBH9'])


def plot_cor_PCA(data): 
    def corr_pca(data):
        pca = PCA(n_components = 4)
        penguins_pca = pca.fit_transform(data)
        coeff = np.transpose(pca.components_[0:2, :])
        n = coeff.shape[0]
        pc_df = pd.DataFrame(data = penguins_pca, 
                             columns = ['PC1', 'PC2','PC3', 'PC4'])
        var_explained = pca.explained_variance_ratio_*100
        def plot_2d_pca(data):
            plt.figure(figsize=(12,10))
            with sbn.plotting_context("notebook",font_scale=1.25):
                sbn.scatterplot(x    = "PC1", 
                                y    = "PC2",
                                data = pc_df, 
                                #hue="Species",
                                #style="Sex",
                                s=100)
                for i in range(n):
            #plot as arrows the variable scores (each variable has a score for PC1 and one for PC2)
                    plt.arrow(0, 
                              0, 
                              coeff[i,0], 
                              coeff[i,1], 
                              color       = 'k', 
                              alpha       = 0.9, 
                              head_width  = 0.02, 
                              head_length = 0.05,
                              linestyle   = '-',
                              linewidth   = 1.5, 
                              overhang    = 0.2)
                    plt.text(coeff[i,0]* 1.15, 
                             coeff[i,1] * 1.15, 
                             list(data.columns.values)[i], 
                             color = 'k', 
                             ha = 'center', 
                             va = 'center',
                             fontsize = 10)  
            plt.xlabel("PC1: "+f'{var_explained[0]:.0f}'+"%")
            plt.ylabel("PC2: "+f'{var_explained[1]:.0f}'+"%")
            plt.show()
        plot_2d_pca(data)
        
        def plot_class_pca(data):
            PC = range(1, pca.n_components_+1)
            plt.bar(PC, pca.explained_variance_ratio_)
            plt.xlabel('Principal Components')
            plt.ylabel('Variance %')
            plt.xticks(PC)
            plt.show()
            # Putting components in a dataframe for later
            PCA_components = pd.DataFrame(penguins_pca)
            inertias = []
            # Creating 10 K-Mean models while varying the number of clusters (k)
            for k in range(1,10):
                model = KMeans(n_clusters=k)
                # Fit model to samples
                model.fit(PCA_components.iloc[:,:3])                
                # Append the inertia to the list of inertias
                inertias.append(model.inertia_)                
            plt.plot(range(1,10), inertias, '-p')
            plt.xlabel('number of clusters, k')
            plt.ylabel('inertia')
            plt.show()
            model = KMeans(n_clusters = 4)
            model.fit(PCA_components.iloc[:,:2])
            labels = model.predict(PCA_components.iloc[:,:2])
            
            plt.figure(figsize = (12, 10))
            with sbn.plotting_context("notebook",font_scale = 1.25):
                sbn.scatterplot(x       = "PC1", 
                                y       = "PC2",
                                data    = pc_df, 
                                #c    = labels,
                                hue     = labels,
                                palette = "deep",
                                s       = 100)
                for i in range(n):
            #plot as arrows the variable scores (each variable has a score for PC1 and one for PC2)
                    plt.arrow(0, 
                              0, 
                              coeff[i,0], 
                              coeff[i,1], 
                              color = 'k', 
                              alpha = 0.9,
                              linestyle = '-',
                              linewidth = 1.5, head_width = 0.02, head_length=0.05,
                              overhang = 0.2)
                    plt.text(coeff[i,0] * 1.15, 
                             coeff[i,1] * 1.15, 
                             list(data.columns.values)[i], 
                             color = 'k', 
                             ha = 'center', 
                             va = 'center',
                             fontsize = 10)  
            plt.xlabel("PC1: "+f'{var_explained[0]:.0f}'+"%")
            plt.ylabel("PC2: "+f'{var_explained[1]:.0f}'+"%")
        plot_class_pca(data)
    corr_pca(data)
    
    def plot_3d_pca(data): 
        # Get the iris dataset
        sbn.set_style("white")
        # Run The PCA
        pca = PCA(n_components = 3)
        pca.fit(data)
        # Store results of PCA in a data frame
        result = pd.DataFrame(pca.transform(data), 
                            columns = ['PCA%i' % i for i in range(3)], 
                            index = data.index)
        # Plot initialisation
        fig = plt.figure(dpi = 600)
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(result['PCA0'], 
                   result['PCA1'], 
                   result['PCA2'], 
                   #c=my_color, 
                   #cmap="Set2_r", 
                   s=2)
         
        # make simple, bare axis lines through space:
        xAxisLine = ((min(result['PCA0']), max(result['PCA0'])), (0, 0), (0,0))
        ax.plot(xAxisLine[0], xAxisLine[1], xAxisLine[2], 'r')
        yAxisLine = ((0, 0), (min(result['PCA1']), max(result['PCA1'])), (0,0))
        ax.plot(yAxisLine[0], yAxisLine[1], yAxisLine[2], 'r')
        zAxisLine = ((0, 0), (0,0), (min(result['PCA2']), max(result['PCA2'])))
        ax.plot(zAxisLine[0], zAxisLine[1], zAxisLine[2], 'r')
         
        # label the axes
        var_explained = pca.explained_variance_ratio_*100
        ax.set_xlabel("PC1: "+f'{var_explained[0]:.0f}'+"%")
        ax.set_ylabel("PC2: "+f'{var_explained[1]:.0f}'+"%")
        ax.set_zlabel("PC3: "+f'{var_explained[2]:.0f}'+"%")
        ax.set_title("PCA on the EIF-corr genes")
        #plt.show()
    plot_3d_pca(data)
plot_cor_PCA(EIF_COR_PRO_sig)


# clustering of PCA plot and retrieve the co-expression genes
def label_cluster(data):
    pca = PCA(n_components = 4)
    penguins_pca = pca.fit_transform(data)
    coeff = np.transpose(pca.components_[0:2, :])
    n = coeff.shape[0]
    pc_df = pd.DataFrame(data = penguins_pca, 
                         columns = ['PC1', 'PC2','PC3', 'PC4'])
    var_explained = pca.explained_variance_ratio_*100
    PC = range(1, pca.n_components_+1)
    pc_df = pd.DataFrame(data = penguins_pca, 
                         columns = ['PC1', 'PC2','PC3', 'PC4'])
    plt.bar(PC, pca.explained_variance_ratio_)
    plt.xlabel('Principal Components')
    plt.ylabel('Variance %')
    plt.xticks(PC)
    plt.show()
    # Putting components in a dataframe for later
    PCA_components = pd.DataFrame(penguins_pca)
    inertias = []
    # Creating 10 K-Mean models while varying the number of clusters (k)
    for k in range(1,10):
        model = KMeans(n_clusters=k)
        # Fit model to samples
        model.fit(PCA_components.iloc[:,:3])                
        # Append the inertia to the list of inertias
        inertias.append(model.inertia_)                
    plt.plot(range(1,10), inertias, '-p')
    plt.xlabel('number of clusters, k')
    plt.ylabel('inertia')
    plt.show()
    model = KMeans(n_clusters = 4)
    model.fit(PCA_components.iloc[:,:2])
    labels = model.predict(PCA_components.iloc[:,:2])
    
    plt.figure(figsize = (12, 10))
    with sbn.plotting_context("notebook",font_scale = 1.25):
        sbn.scatterplot(x       = "PC1", 
                        y       = "PC2",
                        data    = pc_df, 
                        #c    = labels,
                        hue     = labels,
                        palette = "deep",
                        s       = 100)
        for i in range(n):
    #plot as arrows the variable scores (each variable has a score for PC1 and one for PC2)
            plt.arrow(0, 
                      0, 
                      coeff[i,0], 
                      coeff[i,1], 
                      color       = 'k', 
                      alpha       = 0.9,
                      linestyle   = '-',
                      linewidth   = 1.5, 
                      head_width  = 0.02, 
                      head_length = 0.05,
                      overhang = 0.2)
            plt.text(coeff[i,0] * 1.15, 
                     coeff[i,1] * 1.15, 
                     list(data.columns.values)[i], 
                     color = 'k', 
                     ha = 'center', 
                     va = 'center',
                     fontsize = 10)  
    plt.xlabel("PC1: "+f'{var_explained[0]:.0f}'+"%")
    plt.ylabel("PC2: "+f'{var_explained[1]:.0f}'+"%")
        
    df2 = pd.DataFrame(data.index)
    df2["label"] = labels
    #df2['idx'] = df2.groupby('label').cumcount()
    cluster0 = pd.Series(df2['Gene_Symbol'][df2['label'] == 0])
    cluster1 = pd.Series(df2['Gene_Symbol'][df2['label'] == 1])
    cluster2 = pd.Series(df2['Gene_Symbol'][df2['label'] == 2])
    cluster3 = pd.Series(df2['Gene_Symbol'][df2['label'] == 3])
    return (cluster0, cluster1, cluster2, cluster3)
cluster = label_cluster(EIF_COR_PRO_sig)
cluster0 = cluster[0].apply(lambda x: x.split(' ')[0])
cluster0.to_csv('cluster0.csv', index = False) 
cluster1 = cluster[1].apply(lambda x: x.split(' ')[0])
cluster1.to_csv('cluster1.csv', index = False) 
cluster2 = cluster[2].apply(lambda x: x.split(' ')[0])
cluster2.to_csv('cluster2.csv', index = False) 
cluster3 = cluster[3].apply(lambda x: x.split(' ')[0])
cluster3.to_csv('cluster3.csv', index = False) 

EIF4G_COEXP = EIF_COR_PRO_sig[EIF_COR_PRO_sig.index.isin(cluster[1])]
# keep gene symbol only
EIF4G_COEXP.index =  EIF4G_COEXP.index.str.split(' ').str[0]


##############################
###### network analysis ###### 
##############################
# prepare the reference network data
def protein_interaction_reference ():
    #HuRI = pd.read_csv("~/Downloads/HI-union.tsv", sep='\t', header = None)
    ## download the reference file 
    ## https://stringdb-static.org/download/protein.physical.links.detailed.v11.0/9606.protein.physical.links.detailed.v11.0.txt.gz
    #HuRI = pd.read_csv("~/Downloads/9606.protein.physical.links.detailed.v11.0.txt", sep=' ')
    HuRI = pd.read_csv("~/Downloads/9606.protein.links.detailed.v11.0.txt", sep=' ')
    HuRI.head
    HuRI['protein1'] = HuRI['protein1'].str.split('\.').str[-1].str.strip()
    HuRI['protein2'] = HuRI['protein2'].str.split('\.').str[-1].str.strip()
    codes1, uniques1 = pd.factorize(HuRI['protein1'])
    codes2, uniques2 = pd.factorize(HuRI['protein2'])
    # Mapping ensembl gene ids to gene symbols¶
    mg = mygene.MyGeneInfo()
    node1 = mg.querymany(uniques1, 
                        scopes  = 'ensembl.protein', 
                        fields  = 'symbol', 
                        species = 'human',
                        as_dataframe = True)
    node2 = mg.querymany(uniques2, 
                        scopes  = 'ensembl.protein', 
                        fields  = 'symbol', 
                        species = 'human',
                        as_dataframe = True)
    dict1 = pd.Series(node1.symbol.values,index = node1.index).to_dict()
    HuRI['protein1'] = HuRI['protein1'].map(dict1)
    HuRI['protein2'] = HuRI['protein2'].map(dict1)
    return (HuRI)
HuRI = protein_interaction_reference()
# HuRI.to_csv('HuRI.csv', index = False) 

def EIF4F_interaction_reference (x):
    # HuRI = protein_interaction_reference()
    ## construct network for all proteins interacting eIF4G1
    ## find all eIF4G1 interacting proteins from reference data
    EIF4F_RI = HuRI.loc[HuRI['protein1'].isin(x)]
    EIF4F_RI = EIF4F_RI[(EIF4F_RI['experimental'] > 0) & (EIF4F_RI['database'] > 0)] 
    #EIF4F_RI = EIF4F_RI[EIF4F_RI['database'] > 400] 
    EIF4F_RI = EIF4F_RI.dropna()  
    ## construct a list of proteins with interaction to eIF4G1
    EIF4F_RI_List = EIF4F_RI['protein2'].append(pd.Series(x))
    ## construct a network file recording all protein interactions within EIF4G_RI_list
    EIF4F_InterPro_Net = HuRI[HuRI['protein1'].isin(EIF4F_RI_List) & HuRI['protein2'].isin(EIF4F_RI_List)]
    EIF4F_InterPro_Net.reset_index(inplace = True)
    EIF4F_InterPro_Net = EIF4F_InterPro_Net.drop('index', 1)
    EIF4F_InterPro_Net = EIF4F_InterPro_Net[['protein1','protein2','experimental', 'database']]
    # Remove reverse duplicates from dataframe
    cols = ['protein1','protein2']
    EIF4F_InterPro_Net[cols] = np.sort(EIF4F_InterPro_Net[cols].values, axis=1)
    EIF4F_InterPro_Net = EIF4F_InterPro_Net.drop_duplicates()
    #EIF4F_InterPro_Net = EIF4F_InterPro_Net[(EIF4F_InterPro_Net['experimental'] != 0) & 
    #                                        (EIF4F_InterPro_Net['database'] != 0)] 
    melted_df = pd.melt(EIF4F_InterPro_Net, 
                        id_vars=['protein1','protein2'],
                        value_vars=['experimental','database'],
                        var_name='sources',value_name='scores')
    melted_df.loc[melted_df.sources == 'experimental', 'color'] = "r"
    melted_df.loc[melted_df.sources == 'database', 'color'] = "b"
    melted_df_r = melted_df[melted_df['color'] == "r" ] 
    melted_df_r = melted_df_r[(melted_df_r['scores'] != 0)] 
    melted_df_b = melted_df[melted_df['color'] == "b" ] 
    melted_df_b = melted_df_b[(melted_df_b['scores'] != 0)] 
    return (EIF4F_InterPro_Net, melted_df_r, melted_df_b)    
uniques1, uniques2, EIF4F_InterPro_Net = EIF4F_interaction_reference (["EIF4G1"])   

def plot_EIF4F_ref_network (x):
    # uniques1, uniques2, EIF4F_InterPro_Net = EIF4F_interaction_reference (x)
    ## plot nodes interacting eIF4G1 by experiments and database
    ## plot edges for experimental interactions
    EIF4F_InterPro_Net, melted_df_r, melted_df_b = EIF4F_interaction_reference (x)   
    G = nx.from_pandas_edgelist(EIF4F_InterPro_Net,
                                'protein1',
                                'protein2', 
                                edge_attr = ['experimental','database'])
    G.add_nodes_from(nodes_for_adding = EIF4F_InterPro_Net.protein1.tolist())

    labels = [i for i in dict(G.nodes).keys()]
    labels = {i:i for i in dict(G.nodes).keys()}
    
    fig, ax = plt.subplots(figsize = (20,20))
    pos = nx.kamada_kawai_layout(G)
    #pos = nx.spring_layout(G, seed = 50)
    nx.draw_networkx_nodes(G, pos, ax = ax, label =True, node_size = 2000)
    nx.draw_networkx_edges(G, pos, ax=ax)
    _ = nx.draw_networkx_labels(G, pos, labels, ax=ax)
   
    ## plot nodes interacting eIF4G1 by experiment and database
    ## plot edges of experimental interaction.

    G = nx.from_pandas_edgelist(melted_df_r, 
                                'protein1',
                                'protein2', 
                                edge_attr = ['sources','color'],
                                create_using = nx.MultiGraph())
    G.add_nodes_from(nodes_for_adding = EIF4F_InterPro_Net.protein1.tolist())

    #weights = nx.get_edge_attributes(G,'weight').values()
    labels = [i for i in dict(G.nodes).keys()]
    labels = {i:i for i in dict(G.nodes).keys()}
    #pos = nx.kamada_kawai_layout(G)
    #pos = nx.spring_layout(G, seed = 100)
    pos= graphviz_layout(G, prog = "neato")
    fig, ax = plt.subplots(figsize = (20,20))
    nx.draw_networkx_nodes(G, pos, ax = ax, label =True, node_size = 2000)
    nx.draw_networkx_edges(G, pos, edge_color = "r", ax=ax)
    _ = nx.draw_networkx_labels(G, pos, labels, ax=ax)
    
plot_EIF4F_ref_network (["EIF4A1","EIF4G1"]) 


def EIF4F_inter_coexp (x):
    ## construct a network file recording all protein-protein interactions within cluster0
    EIF4F_InterPro_Net = HuRI[HuRI['protein1'].isin(x) & HuRI['protein2'].isin(x)]
    EIF4F_InterPro_Net.reset_index(inplace = True)
    EIF4F_InterPro_Net = EIF4F_InterPro_Net.drop('index', 1)
    EIF4F_InterPro_Net =  EIF4F_InterPro_Net[['protein1',
                                              'protein2',
                                              'experimental',
                                              'database']]
    # Remove reverse duplicates from dataframe
    cols = ['protein1','protein2']
    EIF4F_InterPro_Net[cols] = np.sort(EIF4F_InterPro_Net[cols].values, axis=1)
    EIF4F_InterPro_Net = EIF4F_InterPro_Net.drop_duplicates()
    EIF4F_InterPro_Net = EIF4F_InterPro_Net[(EIF4F_InterPro_Net['experimental'] != 0) & 
                                            (EIF4F_InterPro_Net['database'] != 0)] 
    melted_df = pd.melt(EIF4F_InterPro_Net, 
                        id_vars=['protein1','protein2'],
                        value_vars=['experimental','database'],
                        var_name='sources',value_name='scores')
    melted_df.loc[melted_df.sources == 'experimental', 'color'] = "r"
    melted_df.loc[melted_df.sources == 'database', 'color'] = "b"
    melted_df_r = melted_df[melted_df['color'] == "r" ] 
    melted_df_r = melted_df_r[(melted_df_r['scores'] != 0)] 
    melted_df_b = melted_df[melted_df['color'] == "b" ] 
    melted_df_b = melted_df_b[(melted_df_b['scores'] != 0)]
    return (EIF4F_InterPro_Net, melted_df_r, melted_df_b)    
#uniques1, uniques2, EIF4F_InterPro_Net = EIF4F_inter_coexp (cluster2)   

def plot_EIF4F_coexp_network (x):
    EIF4F_InterPro_Net, melted_df_r, melted_df_b = EIF4F_inter_coexp (x)
    ## plot nodes interacting eIF4G1 by experiments and database
    ## plot edges for experimental interactions
    G = nx.from_pandas_edgelist(EIF4F_InterPro_Net,
                                'protein1',
                                'protein2', 
                                edge_attr = ['experimental','database'])
    G.add_nodes_from(nodes_for_adding = EIF4F_InterPro_Net.protein1.tolist())
    protein = list(EIF4F_InterPro_Net.protein1.unique())

    labels = [i for i in dict(G.nodes).keys()]
    labels = {i:i for i in dict(G.nodes).keys()}
    
    fig, ax = plt.subplots(figsize = (20,20))
    pos = nx.kamada_kawai_layout(G)
    #pos = nx.spring_layout(G, seed = 50)
    # Draw every protein
    nx.draw_networkx_nodes(G, 
                           pos, 
                           ax = ax, 
                           label =True, 
                           node_color='#cccccc', 
                           node_size=100)
    #nx.draw_networkx_nodes(G, pos, 
    #                       nodelist = ["EIF4A1","EIF4G1"], 
    #                       node_color='orange', 
    #                       node_size=100)
    # Draw POPULAR protein
    popular_protein = [item for item in protein if G.degree(item) > 20]
    nx.draw_networkx_nodes(G, pos, 
                           nodelist = popular_protein, 
                           node_color = 'orange', 
                           node_size = 100)
    nx.draw_networkx_edges(G, pos, ax=ax, 
                           width=1, 
                           edge_color="#cccccc")
    _ = nx.draw_networkx_labels(G, pos, labels, ax=ax)

    #Compute largest connected component of the network  (LC)
    #lsit the components in network (g)
    components = nx.connected_components(G)
    #compare among components and find the one having maximun length(LC)
    largest_component = max(components, key=len)
    #largest_component
    # Q1.draw LC
    subgraph = G.subgraph(largest_component)
    #pos = nx.spring_layout(subgraph) # force nodes to separte 
    #pos= graphviz_layout(G, prog = "neato")
    pos = nx.kamada_kawai_layout(G)
    betCent = nx.betweenness_centrality(subgraph, normalized=True, endpoints=True)
    node_color = [20000.0 * G.degree(v) for v in subgraph]
    node_size =  [v * 10000 for v in betCent.values()]
    plt.figure(figsize=(20,15))
    nx.draw_networkx(subgraph, pos = pos, with_labels = False,
                     node_color = node_color,
                     node_size = node_size)
    plt.axis('off')
    
    ## plot nodes interacting eIF4G1 by experiment and database
    ## plot edges of experimental interaction.
    G = nx.from_pandas_edgelist(melted_df_r, 
                                'protein1',
                                'protein2', 
                                edge_attr = ['sources','color'],
                                create_using = nx.MultiGraph())
    G.add_nodes_from(nodes_for_adding = melted_df_r.protein1.tolist())

    #weights = nx.get_edge_attributes(G,'weight').values()
    labels = [i for i in dict(G.nodes).keys()]
    labels = {i:i for i in dict(G.nodes).keys()}
    pos = nx.kamada_kawai_layout(G)
    #pos = nx.spring_layout(G, seed = 50)
    #pos= graphviz_layout(G, prog = "neato")
    fig, ax = plt.subplots(figsize = (20,20))
    nx.draw_networkx_nodes(G, pos, ax = ax, 
                           label =True, 
                           node_size = 100, 
                           cmap=plt.cm.Blues)
    nx.draw_networkx_edges(G, pos, edge_color = "r", ax=ax)
    _ = nx.draw_networkx_labels(G, pos, labels, ax=ax)
           
    fig, ax = plt.subplots(figsize = (20,20))
    # number of nodes to use
    numnodes = 200
    # generate random positions as a dictionary where the key is the node id and the value
    # is a tuple containing 2D coordinates
    positions = {i:(random.random() * 2 - 1, random.random() * 2 - 1) for i in range(numnodes)}
    # use networkx to generate the graph
    network = nx.random_geometric_graph(numnodes, 0.3, pos=positions)
    # then get the adjacency matrix (in sparse form)
    matrix = nx.to_scipy_sparse_matrix(network)
    # run the MCL algorithm on the adjacency matrix and retrieve the clusters
    result = mc.run_mcl(matrix)           # run MCL with default parameters
    clusters = mc.get_clusters(result)    # get clusters
    mc.draw_graph(matrix, clusters, pos=positions, node_size=50, with_labels=False, edge_color="silver")
                  
plot_EIF4F_coexp_network(cluster2)
list(map(plot_EIF4F_coexp_network, [cluster0,cluster1,cluster2,cluster3]))


def EIF4F_inter_coexp_sub (x, y):
    ## construct a network file recording all protein-protein interactions within cluster0
    EIF4F_InterPro = HuRI[HuRI['protein1'].isin(x) & HuRI['protein2'].isin(y)]
    EIF4F_InterPro = EIF4F_InterPro.dropna()  
    ## construct a list of proteins with interaction to eIF4G1
    EIF4F_InterPro_List = EIF4F_InterPro['protein2'].append(pd.Series(x))
    ## construct a network file recording all protein interactions within EIF4G_RI_list
    EIF4F_InterPro_Net = HuRI[HuRI['protein1'].isin(EIF4F_InterPro_List) & 
                              HuRI['protein2'].isin(EIF4F_InterPro_List)]
    
    EIF4F_InterPro_Net.reset_index(inplace = True)
    EIF4F_InterPro_Net = EIF4F_InterPro_Net.drop('index', 1)
    EIF4F_InterPro_Net =  EIF4F_InterPro_Net[['protein1',
                                              'protein2',
                                              'experimental',
                                              'database']]
    # Remove reverse duplicates from dataframe
    cols = ['protein1','protein2']
    EIF4F_InterPro_Net[cols] = np.sort(EIF4F_InterPro_Net[cols].values, axis=1)
    EIF4F_InterPro_Net = EIF4F_InterPro_Net.drop_duplicates()
    EIF4F_InterPro_Net = EIF4F_InterPro_Net[(EIF4F_InterPro_Net['experimental'] != 0) & 
                                            (EIF4F_InterPro_Net['database'] != 0)] 

    melted_df = pd.melt(EIF4F_InterPro_Net, 
                        id_vars=['protein1','protein2'],
                        value_vars=['experimental','database'],
                        var_name='sources',value_name='scores')
    melted_df.loc[melted_df.sources == 'experimental', 'color'] = "r"
    melted_df.loc[melted_df.sources == 'database', 'color'] = "b"
    melted_df_r = melted_df[melted_df['color'] == "r" ] 
    melted_df_r = melted_df_r[(melted_df_r['scores'] != 0)] 
    melted_df_b = melted_df[melted_df['color'] == "b" ] 
    melted_df_b = melted_df_b[(melted_df_b['scores'] != 0)]
    return (EIF4F_InterPro_Net, melted_df_r, melted_df_b)    


def plot_EIF4F_coexp_sub_network (x, y):
    #EIF4F_InterPro_Net, melted_df_r, melted_df_b = EIF4F_inter_coexp (cluster2)
    EIF4F_InterPro_Net, melted_df_r, melted_df_b = EIF4F_inter_coexp_sub (x, y)
    ## plot nodes interacting eIF4G1 by experiments and database
    ## plot edges for experimental interactions
    G = nx.from_pandas_edgelist(EIF4F_InterPro_Net,
                                'protein1',
                                'protein2', 
                                edge_attr = ['experimental','database'])
    G.add_nodes_from(nodes_for_adding = EIF4F_InterPro_Net.protein1.tolist())
    protein = list(EIF4F_InterPro_Net.protein1.unique())

    labels = [i for i in dict(G.nodes).keys()]
    labels = {i:i for i in dict(G.nodes).keys()}
    
    fig, ax = plt.subplots(figsize = (20,20))
    #pos = nx.spring_layout(G,iterations=50)
    pos = nx.kamada_kawai_layout(G)
    #pos = nx.spring_layout(G, seed = 50)
    # Draw every protein
    
    nx.draw_networkx_nodes(G, 
                           pos, 
                           ax = ax, 
                           label =True, 
                           node_color='#cccccc', 
                           node_size=100)
    nx.draw_networkx_nodes(G, pos, 
                           nodelist = x, 
                           node_color='orange', 
                           node_size=100)
    # Draw POPULAR protein
    popular_protein = [item for item in protein if G.degree(item) > 10]
    nx.draw_networkx_nodes(G, pos, 
                           nodelist = popular_protein, 
                           node_color = 'orange', 
                           node_size = 100)
    nx.draw_networkx_edges(G, pos, ax=ax, 
                           width=1, 
                           edge_color="#cccccc")
    _ = nx.draw_networkx_labels(G, pos, labels, ax=ax)

    #Compute largest connected component of the network  (LC)
    #lsit the components in network (g)
    components = nx.connected_components(G)
    #compare among components and find the one having maximun length(LC)
    largest_component = max(components, key=len)
    #largest_component
    # Q1.draw LC
    subgraph = G.subgraph(largest_component)
    #pos = nx.spring_layout(subgraph) # force nodes to separte 
    #pos= graphviz_layout(G, prog = "neato")
    pos = nx.kamada_kawai_layout(G)
    betCent = nx.betweenness_centrality(subgraph, normalized=True, endpoints=True)
    node_color = [20000.0 * G.degree(v) for v in subgraph]
    node_size =  [v * 10000 for v in betCent.values()]
    plt.figure(figsize=(20,15))
    nx.draw_networkx(subgraph, pos = pos, with_labels = False,
                     node_color = node_color,
                     node_size = node_size)
    plt.axis('off')
    
    ## plot nodes interacting eIF4G1 by experiment and database
    ## plot edges of experimental interaction.
    G = nx.from_pandas_edgelist(melted_df_r, 
                                'protein1',
                                'protein2', 
                                edge_attr = ['sources','color'],
                                create_using = nx.MultiGraph())
    G.add_nodes_from(nodes_for_adding = melted_df_r.protein1.tolist())
    labels = [i for i in dict(G.nodes).keys()]
    labels = {i:i for i in dict(G.nodes).keys()}
    pos = nx.kamada_kawai_layout(G)
    #pos = nx.spring_layout(G, seed = 50)
    #pos= graphviz_layout(G, prog = "neato")
    fig, ax = plt.subplots(figsize = (20,20))
    nx.draw_networkx_nodes(G, pos, ax = ax, 
                           label =True, 
                           node_size = 100, 
                           cmap=plt.cm.Blues)
    nx.draw_networkx_edges(G, pos, edge_color = "r", ax=ax)
    _ = nx.draw_networkx_labels(G, pos, labels, ax=ax)              
plot_EIF4F_coexp_sub_network(["EIF4G1"],cluster2)
plot_EIF4F_coexp_sub_network(["EIF4E"],cluster1)


### IRES ###
def plot_ires_heatmap():
    IRES_list = pd.read_csv("~/Downloads/human_IRES_info.csv")
    EIF_COR_PRO_sig['Gene_Symbol'] = EIF_COR_PRO_sig.index
    EIF_COR_PRO_sig['Gene_Symbol'] = EIF_COR_PRO_sig['Gene_Symbol'].apply(lambda x: x.split(' ')[0])
    EIF_COR_PRO_sig.set_index('Gene_Symbol', inplace=True)

    ires_merge = pd.merge(EIF_COR_PRO_sig, IRES_list, how='inner', on=['Gene_Symbol'])
    ires_merge.set_index('Gene_Symbol', inplace=True)

    sbn.clustermap(ires_merge.iloc[:, 0:4],
                       method='centroid',               
                       metric='euclidean',
                       tree_kws=dict(linewidths=.5, colors=(0.2, 0.2, 0.4)),
                       cmap="coolwarm")
plot_ires_heatmap()


def plot_tsne (data): 
    tsne_obj = TSNE(n_components = 2, 
                    random_state = 0).fit_transform(data)
    #tsne_obj
    #tsne_em = TSNE(n_components = 2, perplexity   = 30.0, n_iter       = 1000, verbose      =1).fit_transform(data)
    tsne_df = pd.DataFrame({'X':tsne_obj[:,0],
                            'Y':tsne_obj[:,1]})
    tsne_df.head()
    plt.figure(figsize = (16,10))
    sbn.scatterplot(x = "X", 
                    y = "Y",
                    data = tsne_df)
plot_tsne(EIF_COR_PRO_sig)


# Displaying dataframe as an heatmap  
# with diverging colourmap as coolwarm 
D = sch.distance.pdist(EIF_COR_PRO_sig, metric='euclidean')
L = sch.linkage(D, method='centroid')

sch.dendrogram(sch.linkage(D, method='centroid'), 
               orientation='top',  p=5, truncate_mode='level', 
               color_threshold=.665)

h = sbn.clustermap(EIF_COR_PRO_sig,
                   method='centroid',               
                   metric='euclidean',
                   tree_kws=dict(linewidths=.5, colors=(0.2, 0.2, 0.4)),
                   cmap="coolwarm")


def plot_heatmap (data):
    # clustering colums
    data_1D_X = ssd.pdist(data.T, 'euclidean')
    X = sch.linkage(data_1D_X, method='centroid')
    # clustering rows
    data_1D_Y = ssd.pdist(data, 'euclidean')
    Y = sch.linkage(data_1D_Y, method='centroid')
    #plot first dendrogram
    fig = plt.figure(figsize=(8, 8))
    
    #sch.set_link_color_palette(['m', 'c', 'y', 'k','g','r'])
    ax1 = fig.add_axes([0.09, 0.1, 0.2, 0.6])
    Z1 = sch.dendrogram(Y, orientation='left',color_threshold= 0.65)
    ax1.set_xticks([])
    ax1.set_yticks([])
    
    # second dendrogram.
    ax2 = fig.add_axes([0.3, 0.71, 0.6, 0.1])
    Z2 = sch.dendrogram(X)
    ax2.set_xticks([])
    ax2.set_yticks([])
    
    # plot matrix
    axmatrix = fig.add_axes([0.3, 0.1, 0.6, 0.6])
    # sorts based of clustering
    idx1 = Z1['leaves']
    idx2 = Z2['leaves']
    D = data.values[idx1, :]
    D = D[:, idx2]
    im = axmatrix.matshow(D, aspect='auto', origin='lower', cmap="coolwarm")
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])
plot_heatmap(EIF_COR_PRO_sig)





###########################################
## load RNA-seq data for CCLE cell lines ##
###########################################
## https://portals.broadinstitute.org/ccle/data
def ccle_rna():
    ## https://depmap.org/portal/download/
    CCLE_RNA = pd.read_csv("~/Downloads/CCLE_expression_full.csv")
    ccle = CCLE_RNA.set_index('Unnamed: 0')
    ccle.index.names = ['DepMap_ID']
    return ccle
CCLE_RNA = ccle_rna()

def eif_corr_rna(eif):
    y = pd.DataFrame() 
    for gene in eif:
        print(gene)
        x =  CCLE_RNA.corrwith(CCLE_RNA[gene])
        x.name = gene
        y = pd.concat([y, x], axis=1)
    return y
#gee = ['EIF4G1 (ENSG00000114867)']
#EIF_COR_RNA = eif_corr_rna(gee) 
def plot_eif_corr_rna(gene):
    EIF_COR_RNA = eif_corr_rna(gene)
    EIF_COR_RNA_sig = EIF_COR_RNA.loc[(EIF_COR_RNA['EIF4G1 (ENSG00000114867)'] >= 0.3) 
                                      | (EIF_COR_RNA['EIF4G1 (ENSG00000114867)'] <= -0.3)
                                      | (EIF_COR_RNA['EIF4A1 (ENSG00000161960)'] >= 0.3) 
                                      | (EIF_COR_RNA['EIF4A1 (ENSG00000161960)'] <= -0.3) 
                                      | (EIF_COR_RNA['EIF4E (ENSG00000151247)'] >= 0.3) 
                                      | (EIF_COR_RNA['EIF4E (ENSG00000151247)'] <= -0.3) 
                                      | (EIF_COR_RNA['EIF4EBP1 (ENSG00000187840)'] >= 0.3) 
                                      | (EIF_COR_RNA['EIF4EBP1 (ENSG00000187840)'] <= -0.3) 
                                      ] 
    EIF_COR_RNA_sig.dropna(inplace=True)
    EIF_COR_RNA_sig.dtypes
    #EIF_COR_PRO['index1'] = EIF_COR_PRO.index
    # Displaying dataframe as an heatmap  
    # with diverging colourmap as coolwarm 
    sbn.clustermap(EIF_COR_RNA_sig,               
                   #metric='correlation',
                   #standard_scale=1, 
                   cmap="coolwarm")
EIF =  pd.Index(['EIF4G1 (ENSG00000114867)', 
                 'EIF4A1 (ENSG00000161960)',
                 'EIF4E (ENSG00000151247)',
                 'EIF4EBP1 (ENSG00000187840)'])
plot_eif_corr_rna(EIF)

#pd.Series(CCLE_PRO["Gene_Symbol"]).is_unique
#pd.Series(CCLE_PRO_subset.index).is_unique

#ids = CCLE_PRO["Gene_Symbol"]
#du = CCLE_PRO[ids.isin(ids[ids.duplicated()])].sort_values("Gene_Symbol")

def gdsc_drugrespose():
## http://www.cancerrxgene.org/downloads
    GDSC_drugrespose = pd.read_csv("~/Downloads/PANCANCER_IC_Mon Nov  9 16_46_39 2020.csv")
    GDSC_drugrespose.dtypes # The data type of each column.
    # convert "Drug name" from object to category 
    #GDSC_drugrespose['Drug name'] = GDSC_drugrespose['Drug name','Cell line name'].astype('category')
    for col in ['Drug name', 'Drug Id', 'Cell line name', 'Tissue']:
        GDSC_drugrespose[col] = GDSC_drugrespose[col].astype('category')
    return GDSC_drugrespose
GDSC_drugrespose = gdsc_drugrespose()
drug_list = GDSC_drugrespose['Drug name'].cat.categories


def cpd_corr_pro(a_drug_list, protein):
    CCLE_PRO_EIF = CCLE_PRO_subset[EIF]
    CCLE_PRO_EIF.index = CCLE_PRO_EIF.index.str.split('_').str[0]
    CCLE_PRO_EIF.index.names = ['Cell line name']
    y = pd.DataFrame() 
    for drug in a_drug_list:
        print(drug)
## select one drug type
        GDSC_drugrespose_subset = GDSC_drugrespose.loc[GDSC_drugrespose['Drug name'] == drug]
        GDSC_drugrespose_subset['Cell line name'] = GDSC_drugrespose_subset['Cell line name'].str.replace('-', 
                                                                                                          '', 
                                                                                                        regex=True)
        GDSC_drugrespose_subset.set_index('Cell line name', inplace=True)
        GDSC_drugrespose_subset = GDSC_drugrespose_subset.drop(GDSC_drugrespose_subset.columns[[0,1,2,3,4,5]], 
                                                               axis=1)
        GDSC_drugrespose_subset = GDSC_drugrespose_subset.drop(['Max conc', 
                                                                'RMSE', 
                                                                'Dataset version', 
                                                                'Z score', 
                                                                'AUC'], 1)
        # select one gene expression
        # CCLE_expression_subset = CCLE_expression.loc[CCLE_expression['Description'] == 'EIF4G1']
        result = pd.merge(GDSC_drugrespose_subset, CCLE_PRO_EIF, on='Cell line name')
        x =  result.corrwith(result['IC50'])
        x.name = drug
        y = pd.concat([y, x], axis=1)
    return y
#cpd_pro_correlation = cpd_corr_pro(drug, CCLE_PRO_EIF)  
CPD_COR_EIF = cpd_corr_pro(drug_list, CCLE_PRO_subset).T
# CPD_COR_EIF['indx'] = CPD_COR_EIF.index
CPD_COR_EIF = CPD_COR_EIF.drop(['IC50'], 1)
CPD_COR_EIF_sig = CPD_COR_EIF.loc[(CPD_COR_EIF['EIF4G1 Q04637-9'] >= 0.3) | (CPD_COR_EIF['EIF4G1 Q04637-9'] <= -0.3)
                                  #|(EIF_COR_PRO['EIF4G1 Q04637-8'] >= 0.3) |  (EIF_COR_PRO['EIF4G1 Q04637-8'] <= -0.3) 
                                  | (CPD_COR_EIF['EIF4A1 P60842'] >= 0.3) |  (CPD_COR_EIF['EIF4A1 P60842'] <= -0.3) 
                                  | (CPD_COR_EIF['EIF4E P06730-2'] >= 0.3) |  (CPD_COR_EIF['EIF4E P06730-2'] <= -0.3) 
                                  | (CPD_COR_EIF['EIF4EBP1 Q13541'] >= 0.3) |  (CPD_COR_EIF['EIF4EBP1 Q13541'] <= -0.3) 
                                  #| (CPD_COR_EIF['MKNK2 Q9HBH9'] >= 0.3) |  (CPD_COR_EIF['MKNK2 Q9HBH9'] <= -0.3) 
                                  #| (CPD_COR_EIF['MKNK1 Q9BUB5'] >= 0.3) |  (CPD_COR_EIF['MKNK1 Q9BUB5'] <= -0.3) 
                                  ] 
CPD_COR_EIF_sig.dropna(inplace=True)
sbn.clustermap(CPD_COR_EIF_sig,
                   method='centroid',               
                   metric='euclidean',
                   tree_kws=dict(linewidths=.5, colors=(0.2, 0.2, 0.4)),
                   cmap="coolwarm")







def ccle_rna_annotation():
    ## load annontation data for CCLE cell lines
    ## https://depmap.org/portal/download/
    CCLE_annotation = pd.read_csv("~/Downloads/sample_info.csv")
    CCLE_annotation_subset = CCLE_annotation[['DepMap_ID', 'stripped_cell_line_name']]
    CCLE_annotation_subset.set_index('DepMap_ID', inplace=True)
    CCLE_RNA_annotation = pd.merge(CCLE_RNA, CCLE_annotation_subset,on = 'DepMap_ID')
    CCLE_RNA_annotation.set_index('stripped_cell_line_name', inplace=True)
    CCLE_RNA_annotation.index.names = ['Cell line name']
    return CCLE_RNA_annotation
CCLE_RNA_annotation = ccle_rna_annotation()
CCLE_RNA_EIF = CCLE_RNA_annotation[['EIF4G1 (ENSG00000114867)', 
                                    'EIF4A1 (ENSG00000161960)',
                                    'EIF4E (ENSG00000151247)',
                                    'EIF4EBP1 (ENSG00000187840)']]

def cpd_corr_rna(a_drug_list,rna):
    y = pd.DataFrame() 
    for drug in a_drug_list:
        print(drug)
        GDSC_drugrespose_subset = GDSC_drugrespose.loc[GDSC_drugrespose['Drug name'] == drug]
        GDSC_drugrespose_subset['Cell line name'] = GDSC_drugrespose_subset['Cell line name'].str.replace('-', '', regex=True)
        GDSC_drugrespose_subset.set_index('Cell line name', inplace=True)
        GDSC_drugrespose_subset = GDSC_drugrespose_subset.drop(GDSC_drugrespose_subset.columns[[0,1,2,3,4,5]], axis=1)
        GDSC_drugrespose_subset = GDSC_drugrespose_subset.drop(['Max conc', 'RMSE', 'Dataset version'], 1)
        # select one gene expression
        # CCLE_expression_subset = CCLE_expression.loc[CCLE_expression['Description'] == 'EIF4G1']
        result = pd.merge(GDSC_drugrespose_subset, rna, on='Cell line name')
        x =  result.corrwith(result['IC50'])
        x.name = drug
        y = pd.concat([y, x], axis=1)
    return y
cpd_rna_correlation = cpd_corr_rna(drug_list, CCLE_RNA_EIF).T    
cpd_rna_correlation = cpd_corr_rna(['Camptothecin'], CCLE_RNA_EIF).T    


EIF4F_corr = cpd_rna_correlation.loc[cpd_rna_correlation.index.isin(EIF4F["Name"])]
EIF4F_corr_T = EIF4F_corr.T
cpd_rna_correlation['index1'] = cpd_rna_correlation.index




# select one drug type
GDSC_drugrespose_subset = GDSC_drugrespose.loc[GDSC_drugrespose['Drug name'] == 'Camptothecin']
GDSC_drugrespose_subset['Cell line name'] = GDSC_drugrespose_subset['Cell line name'].str.replace('-', '', regex=True)
GDSC_drugrespose_subset.set_index('Cell line name', inplace=True)
GDSC_drugrespose_subset = GDSC_drugrespose_subset.drop(GDSC_drugrespose_subset.columns[[0,1,2,3,4,5]], axis=1)
GDSC_drugrespose_subset = GDSC_drugrespose_subset.drop(['Max conc', 'RMSE', 'Dataset version'], 1)

# select one gene expression
# CCLE_expression_subset = CCLE_expression.loc[CCLE_expression['Description'] == 'EIF4G1']
CCLE_expression = CCLE_expression.set_index('Name')
CCLE_expression = CCLE_expression.drop(['Description'], 1)
CCLE_expression_t = CCLE_expression.T 
CCLE_expression_t.dtypes
CCLE_expression_t.index = CCLE_expression_t.index.str.split('_').str[0]
CCLE_expression_t.index.names = ['Cell line name']
result = pd.merge(GDSC_drugrespose_subset, CCLE_expression_t, on='Cell line name')
y =  result.corrwith(result['IC50'])
y.name = "Camptothecin"
y =  y.to_frame()












new_header = CCLE_expression_transposed.iloc[0] #grab the first row for the header
CCLE_expression_transposed = CCLE_expression_transposed[1:] #take the data less the header row
CCLE_expression_transposed.columns = new_header #set the header row as the df header
CCLE_expression_transposed_drop = CCLE_expression_transposed.drop(index='Description')
CCLE_expression_transposed_drop_reset = CCLE_expression_transposed_drop.reset_index()
CCLE_expression_transposed_drop_reset['Cell line name'] = CCLE_expression_transposed_drop_reset['index'].str.split('_').str[0]
CCLE_expression_transposed_drop_reset_rename = CCLE_expression_transposed_drop_reset.drop(['index'], axis=1)
x = CCLE_expression_transposed_drop_reset_rename['Cell line name']
y = list(CCLE_expression_transposed_drop.index.values)

result = pd.merge(GDSC_drugrespose_subset, CCLE_expression_transposed_drop_reset_rename, on='Cell line name')
result_1 = result.drop(result.columns[[0,1,2,3,4,5,6]], axis=1)
result_1 = result_1.drop(['Max conc', 'RMSE','Dataset version'], 1)

print (result_1[['ENSG00000210195.2', 'ENSG00000210196.2']].corr(result_1['IC50']))

result_1.dtypes
