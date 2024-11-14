#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to generate the Fig4 of the manuscript. 
Code tested using python 3.10.6
"""

#%% import modules

import pandas as pd # 1.5.1
import numpy as np # 1.26.4
import math
from skbio.diversity import alpha_diversity # 0.5.8
import plotly.express as px # 5.10.0
import plotly.io as pio
pio.renderers.default='svg'

#%% Load tables

# Sheet 1: Contains HLA genotypes from present study + Immel et al 2021 (Additional file 1: Data S1) plus calls 
# from the 1KG populations (https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HLA_types/20181129_HLA_types_full_1000_Genomes_Project_panel.txt).
# Sheet 2: contains the modern German allele frequencies (https://www.allelefrequencies.net/pop6001c.asp?pop_id=3767)
hla_file =  '../data/hla_calls.xlsx'

# directory to output plots and tables
outdir = './'

# load HLA genotypes
dfo = pd.read_excel(hla_file,  header = 0, sheet_name="genotypes")
dfo = dfo.replace('..', np.nan) 


# load modern German reference frequencies
ref = pd.read_excel(hla_file,  header = 0, sheet_name="allelefrequencies.net")
ref["Counts"] = ref["Allele Frequency"] * ref['Sample Size']
ref["Counts"] = ref["Counts"].apply(lambda x: round(x))
ref["Population"] = "GER" # relabel the reference for plot
ref["HLA"] = "HLA-" + ref["HLA"] # adjust the column to subset

## define frequency filter to use (do not plot below this freq)
freq_filter = 0.01

#%% adjust tables

df = dfo.copy()
# add the 1000 genomes with the pop instead of superpop
invert_1kg = df[df["Group"].isin([ 'EUR', 'EAS', 'SAS', 'AFR', 'AMR'])]
invert_1kg = invert_1kg.rename({"Group":"Site", "Site":"Group"}, axis=1)
df = pd.concat([df,invert_1kg])
print("Size of populations with HLA genotypes:\n", df.Group.value_counts())

#%% determine the dict for the columns

dict_cols = {
    "HLA-A":['HLA-A (1)', 'HLA-A (2)'], 
    "HLA-B":['HLA-B (1)', 'HLA-B (2)'],
    "HLA-C":['HLA-C (1)', 'HLA-C (2)'], 
    "HLA-DPB1":['HLA-DPB1 (1)', 'HLA-DPB1 (2)'],
    "HLA-DQB1":['HLA-DQB1 (1)', 'HLA-DQB1 (2)'], 
    "HLA-DRB1":['HLA-DRB1 (1)', 'HLA-DRB1 (2)'],
    }

res = 2 # digits resolution to use for frequencies

#%% 

# plot order of using the groups and the color for plots
plot_order = {

     'EF' : px.colors.qualitative.Plotly[1],
     "LF": px.colors.qualitative.Plotly[0],
     "GER": px.colors.qualitative.Plotly[4],
       
      "CEU": px.colors.qualitative.Safe[4],
      "GBR": px.colors.qualitative.Safe[5],
      "FIN": px.colors.qualitative.Safe[3],
      "IBS": px.colors.qualitative.Safe[7],
      "TSI": px.colors.qualitative.Safe[9],

    }

# specific the groups where only allele frequency is available and not the 
# actual genotypes -> leads to different way to get shannon diversity
ref_groups = ["GER"]

## filter df to only include groups to plot
group_col = "Group"
df = df[df[group_col].isin(plot_order.keys())]


#%% adjust resolution

def isNaN(num): # function to check if is nan
    return num != num

def res_filter(x): # function to remove genotypes with lower res
    if isNaN(x) == False:
        nfields = len(x.split(":")[0:res])
        test = "PASS" if nfields == res else "FAIL"
    else:
        test = "FAIL"
    return x if test == "PASS" else np.nan


cols_list = [item for sublist in dict_cols.values() for item in sublist]
for col in cols_list:
    # adjust digit resolution to res (if more, then decrease)
    df[col] = df[col].apply(lambda x: ":".join(x.split(":")[0:res]) if isNaN(x) == False else x) 
    # remove lower res
    df[col] = df[col].apply(lambda x: res_filter(x) )

#%% get dict with dataframes with all freqs per site

dict_dfs = {}

for key, l in dict_cols.items():
    dfA = pd.concat([
        df[[group_col, 'Sample ID', l[0] ]].rename({l[0]:key},axis=1) ,  
        df[[group_col, 'Sample ID', l[1] ]].rename({l[1]:key},axis=1)
        ]).sort_values([group_col,"Sample ID"]).groupby([group_col,key]).size().reset_index(name='Counts')
    
    dfA["Frequency"] = dfA.apply(lambda x: round(x["Counts"]/((dfA[dfA[group_col]==x[group_col]]["Counts"].sum())),3), axis=1) # AC
    dfA["Sample Size"] = dfA.apply(lambda x: round(((dfA[dfA[group_col]==x[group_col]]["Counts"].sum())),3), axis=1) # "Sample  size" here refers to the total number of alleles- NCHROBS
    
    # merge with data from allelefrequencies.net
    dfA = pd.concat([
        dfA,
        ref[ref["HLA"]==key].rename({'Allele':key,"Population":group_col,'Allele Frequency':'Frequency'}, axis=1)[[group_col,key,"Counts","Frequency", "Sample Size"]] # all alleles in reference
        ])
    
    dfA["Allele"] = key
    dfA.rename({key:"HLA"},axis=1, inplace=True) # adjust column name before concat 
    
    dict_dfs[key] = dfA

dfu = pd.concat( dict_dfs.values() ) # merge dfs

#%% Diversity: calculate Shannon Index (H')

# check pop sizes
# df[group_col].value_counts()
                    

# generate Shannon index for diversity (H')
shannon_l = []

for k, i in dict_cols.items():
    # identify the pop with smallest sample size for each allele
    dfs = df[[group_col,"Sample ID"]+i]
    dfs = dfs[ (~dfs[i[0]].isna() ) & (~dfs[i[1]].isna()) ] # remove nans
    gsize = dfs.groupby("Group").agg(counts=("Sample ID", "count")).sort_values(["counts"]).reset_index()
    resample_size = gsize["counts"][0]
    smalles_n_g = gsize[group_col][0] 
    print("="*50)
    print(f"Group with smallest sample size for {k}: {smalles_n_g} (n={resample_size})")
    
    # get shannon index for each group
    for g in plot_order.keys():
        dfs0 = dfu[(dfu[group_col]==g) & (dfu["Allele"]==k)]
        if len(dfs0) > 0:
            for b in range(100): # bootstrapping
                dfs0b = dfs0.sample(n=resample_size*2, replace=True, weights="Frequency")
                dfs2b = pd.DataFrame(dfs0b.HLA.value_counts())
                counts_b = list(dfs2b.HLA)
                H_b = alpha_diversity('shannon', counts_b, base=math.e).values[0]
                shannon_l.append([k,g,H_b,f"Bootstrap{b}",resample_size])

        

shannon = pd.DataFrame(shannon_l, columns=["Locus", group_col, "Shannon Index (H')", "Resample #", "Resample size"])


#%% plot shannon index (Fig 4)

# adjust the name of each facet
shannon["facet_col"] = "<b>" + shannon["Locus"] + "</b>" + " (n=" + shannon["Resample size"].astype(str)  + ")"

fig = px.box(shannon, 
             x=group_col, 
             y="Shannon Index (H')", 
             facet_col="facet_col",
             color=group_col,
             color_discrete_sequence=list(plot_order.values()),
             category_orders={"Group": plot_order.keys()} , 
             facet_col_wrap=3,
             labels={"Group":"Population"}
             )



fig.for_each_annotation(lambda a: a.update(text=a.text.replace("facet_col=","")))


fig.update_layout(
    # autosize=False,
    width=700,  # 800
    height=500, # 1300
    
    # font_size=12,
    
    # xaxis_title="",
    # yaxis_title="",
    
    
    plot_bgcolor= "white", 
    paper_bgcolor= "white" ,
    
    
    legend=dict(
        orientation="v", 
        yanchor="top",  
        y=1, 
        xanchor="right", 
        x=1.2
    )
    )
    

fig.update_xaxes(title_text="", showticklabels=False)
fig.update_yaxes(title_text="")


fig.update_xaxes(showline = True, linecolor = 'black', linewidth = 1,  mirror = True,
                 # showgrid=True, gridwidth=0.3, gridcolor='lightgrey',
                 )
fig.update_yaxes(showline = True, linecolor = 'black', linewidth = 1,  mirror = True,
                  showgrid=True, gridwidth=0.3, gridcolor='lightgrey', 
                   range=[-0.05, 3.4],  dtick=1, 
                 )

fig.update_traces(
    # quartilemethod="exclusive", # or "inclusive", or "linear"
    marker=dict(
                size=3, 
                opacity=0.5,
                # line=dict(width=0.19, color='black', )
                )
    )

# fig.show()
# fig.show(renderer="browser")

fig.write_image(f"{outdir}/Fig4.pdf", engine='kaleido')
fig.write_image(f"{outdir}/Fig4.svg", engine='kaleido')
