#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to generate the Fig2 of the manuscript. 
Code tested using python 3.10.6
"""

#%% import modules

import pandas as pd # 1.5.1
import numpy as np # 1.26.4
from scipy import stats # 1.13.1
from statsmodels.sandbox.stats.multicomp import multipletests # 0.14.2
import itertools
# import random
import plotly.express as px # 5.10.0
import plotly.io as pio
pio.renderers.default='svg'
# random.seed(42) 

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

# plot order if using the groups and the color for plots
plot_order = {

     'EF' : px.colors.qualitative.Plotly[1],
     "LF": px.colors.qualitative.Plotly[0],
     "GER": px.colors.qualitative.Plotly[4],
       
      #  "CEU": "rgb(156, 156, 94)",       ### 1000 genomes pops of interest (EUR)
      #  "GBR": "rgb(102,102,102)",
      #  "TSI": "rgb(179,179,179)",
      #  "IBS": "#8C564B",
      #  "FIN": "rgb(102, 17, 0)",
    
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

dfu = pd.concat( dict_dfs.values() ) # Data S10



#%% only keep loci for plot that have more than 1% in at least one of the groups
alleles_tokeep = dfu[dfu["Frequency"] > freq_filter].HLA.tolist()
dff = dfu[dfu["HLA"].isin(alleles_tokeep)].reset_index(drop=True)

# recreate dict with fitlered alleles
dict_dfs = {}
for al in dff.Allele.unique():
    subset = dff[dff["Allele"]==al]
    dict_dfs[al] = subset
      

#%% individual plots used to create multi-panel Figure S14

# for key, dfA in dict_dfs.items():
#     fig = px.bar(dfA,  x="HLA", y="Frequency", title=key,
#                   color=group_col, barmode='group',
#                   # text="Counts",
#                   category_orders = {group_col:plot_order.keys(), "HLA":sorted(list(dfA["HLA"].unique()))},
#                   height=350, width=15*len(dfA),
#                   color_discrete_sequence=list(plot_order.values()),
#                   labels = {"HLA":""} )
#     fig.update_xaxes(tickangle=45)
#     fig.update_traces(textfont_size=12, textangle=0, textposition="outside", cliponaxis=False)
#     fig.update_yaxes(tick0 = 0, dtick = 0.10, showline=True, linewidth=2, linecolor='black', gridcolor='lightgrey',
#                       range=[0,round(dff.Frequency.max()+0.1,1)]
#                       ) 

    
#     fig.show()
#     # To save plot
#     # fig.write_image(f"{outdir}/FigS14.panel-{key}.svg", engine='kaleido') 
#     # fig.write_image(f"{outdir}/FigS14.panel-{key}.pdf", engine='kaleido')



#%% pairwise Fisher tests to check if there are significant difference between pops

fisher_results = []
for allele in dff["HLA"].unique():
    df_subset = dff[dff["HLA"]==allele]
    
    # test if more than two pops for the calculations
    if len(df_subset) >= 2:
        list_groups = df_subset["Group"].tolist()
        group_combinations = list(itertools.combinations(list_groups, r=2))
        for combination in group_combinations:
            try:
                data = df_subset[df_subset["Group"].isin(combination)][["Counts","Sample Size"]]
                data["OtherAlleleAC"] = data["Sample Size"] - data["Counts"] 
                data = data[["Counts", "OtherAlleleAC"]] 
                _listMAFs_ = list(df_subset[df_subset["Group"].isin(combination)]["Frequency"])
                difference = abs(float(_listMAFs_[0]) - float(_listMAFs_[1]))
                oddsratio, p_value = stats.fisher_exact(table=data, alternative="two-sided")
                fisher_results.append([allele, combination[0], combination[1],p_value, oddsratio, difference])
            except:
                print("WARNING: Unable to perform Fisher test for SNP %s: %s, %s" % (allele, combination[0], combination[1]))
                fisher_results.append([allele, combination[0], combination[1],float("nan"), float("nan"), float("nan")])



df_results_fisher = pd.DataFrame(fisher_results, columns=["HLA allele","Pop1", "Pop2", "p-value (Fisher)","OR (Fisher)","Freq. difference"])

# correct for multiple testing - using FDR TSBH
pvals = df_results_fisher["p-value (Fisher)"].tolist()

reject, pvals_corrected, alphacSidak, alphacBonf = multipletests(
    pvals, 
    alpha=0.05, 
    method='fdr_tsbh', 
    is_sorted=False, 
    returnsorted=False
    )

df_results_fisher["p-value (fdr_tsbh)"] = pvals_corrected

#%% filter tests for p<0.05 and diff of freqs > 0.1 

dft = df_results_fisher[
    (df_results_fisher["p-value (fdr_tsbh)"]< 0.05)
    & (df_results_fisher["Freq. difference"] >= 0.1)
    ] 

# annotated dft
def annotate_p_value(p_value):
    if p_value <= 0.001:
        return '***'
    elif p_value <= 0.01:
        return '**'
    elif p_value <= 0.05:
        return '*'
    else:
        return 'ns'
    
dft["anno"] = dft["p-value (fdr_tsbh)"].apply(lambda x: annotate_p_value(x)) # Table S4


#%% # Save all tables in xlsx
dfs_dict = { "Frequencies":dfu,  "Sig. different":dft  }
writer = pd.ExcelWriter(f"{outdir}/results_HLA_freqs.xlsx", )
for name, dfi in dfs_dict.items():
    try:
        dfi.to_excel(writer, sheet_name=name, index=False)
    except:
        pass
writer.close()


#%% Plot facet plot for alleles showing significant changes - Fig2

test_df = dff[dff["HLA"].isin(dft["HLA allele"].unique())].sort_values("HLA").reset_index(drop=True)


fig = px.bar(test_df, x="Group", y="Frequency", facet_col="HLA", facet_col_wrap=6,
             color="Group", color_discrete_sequence=list(plot_order.values()),
             category_orders={"Group": plot_order.keys()} , 
             )

fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[1]))

dict_facet = {hla:i for i,hla in  enumerate([x.text for x in fig.layout.annotations], start=1)} # make dict for annotating correctly


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
        orientation="h", 
        yanchor="top",  
        y=-0.02, 
        xanchor="center", 
        x=0.5 
    ),
    
    bargap=0.5,  
    # margin=dict(l=20, r=20, t=40, b=20),  


)

fig.update_xaxes(showline = True, linecolor = 'black', linewidth = 1,  mirror = True,
                 # showgrid=True, gridwidth=0.3, gridcolor='lightgrey',
                 )
fig.update_yaxes(showline = True, linecolor = 'black', linewidth = 1,  mirror = True,
                  showgrid=True, gridwidth=0.3, gridcolor='lightgrey', 
                   range=[-0.05, 1.05],  dtick=0.25
                 )

fig.update_xaxes(title_text="", showticklabels=False)
# fig.update_yaxes(title_text="")

# add statistical annotations
l_y = []
for _, row in dft.sort_values("HLA allele").reset_index(drop=True).iterrows():

    # try:
        allele = row['HLA allele']
        group1 = row['Pop1']
        group2 = row['Pop2']
        p_value = row['p-value (fdr_tsbh)']
        annotation_text = row['anno']
        
        # get x positions of the groups
        x0 = list(plot_order.keys()).index(group1)
        x1 = list(plot_order.keys()).index(group2)
        
        # get facet number corresponding to the allele
        facet_number = dict_facet[allele] 
        
        # adjust l_y
        adjust_y = l_y.count(facet_number) * 0.15
        l_y.append(facet_number) 
        
        
        # determine the y position 
        y_position = test_df[test_df['HLA'] == allele]['Frequency'].max() + 0.05 + adjust_y
        
        # line for the comparison
        fig.add_shape(
            type="line",
            x0=x0, x1=x1,
            y0=y_position, y1=y_position,
            line=dict(color="black", width=1),
            xref=f"x{facet_number}", yref=f"y{facet_number}"
        )
        
        # asterisks or p-value
        fig.add_annotation(
            text=annotation_text,
            x=(x0 + x1) / 2, 
            y=y_position + 0.05,  # slightly above the line
            xref=f"x{facet_number}", yref=f"y{facet_number}",
            showarrow=False,
            font=dict(color="black", size=12)
        )
        
        
        # deebug
        # print(f"{allele}\n{facet_number}\n{x0}, {x1}\n{y_position}")
    
    # except:
    #     # print(_,row)
    #     pass



# fig.show()
fig.write_image(f"{outdir}/Fig2_significantHLA.pdf", engine='kaleido')
fig.write_image(f"{outdir}/Fig2_significantHLA.svg", engine='kaleido')



