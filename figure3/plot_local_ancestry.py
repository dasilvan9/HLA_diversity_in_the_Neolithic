#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to generate the Fig3 of the manuscript. 
Code tested using python 3.10.6
"""

import pandas as pd # 1.5.1
import os
from scipy.stats import mannwhitneyu # 1.13.1
import statsmodels.stats.multitest as smm # 0.14.2
from glob import glob
import plotly.express as px  # 5.10.0
import plotly.io as pio
pio.renderers.default = 'svg'


mydir = "/home/nicolas/hamburg" # mount folder, leave as "" if data in local pc

## If running RFmix v2 again, please uncomment lines below to load the original *fb.tsv files
# fb_files = glob(f"{mydir}/data110/dasilva/Project_Neolithic_Pops/1240K_run/RFmix2/run2/LF*G22*fb.tsv")
# fb_files_dict = {int(x.split("/")[-1].split(".")[1].replace("chr","")): x for x in fb_files }


outdir = './'  # output directory
# outname = "Results_LF.WHG.G22" # name used in output files

ancestry = "WHG" # ancestry name to subset from RFmix output

supdataf = "../data/AdditionalFile1.xlsx"

#%% Concat data and calculate mean ancestry per SNP

# if os.path.exists(f"{outdir}/{outname}.tsv"):
#     df = pd.read_csv(f"{outdir}/{outname}.tsv.zip", sep="\t")
# else:
#     dfl = []
#     for chrom in range(1,23):
#         data = pd.read_csv(fb_files_dict[chrom], sep="\t", header=0, skiprows=1)    
#         ancestry_columns = [col for col in data.columns if f":::{ancestry}" in col]
#         data['mean_painting'] = data[ancestry_columns].mean(axis=1)
#         dfl.append(data)
    
#     df = pd.concat(dfl)
    
#     # add genome wide z-scores
#     df['zscore'] = (df['mean_painting'] - df['mean_painting'].mean()) / df['mean_painting'].std()
    
#     df.to_csv(f"{outdir}/{outname}.tsv.zip",sep = "\t", index=False)

df = pd.read_csv("../data/Results_LF.WHG.G22.tsv.zip", sep="\t")

#%% Plot mean ancestry for HLA region - Fig 3A

# subset whole MHC region "chr6:28477797-33448354" + BUFFER
buffer = 5000000 # 5 Mb
hla_plus_buffer = df[(df["chromosome"]==6) 
         & (df["physical_position"]>=(28477797 - buffer)) 
         & (df["physical_position"]<=(33448354 + buffer)) 
         ]

fig = px.scatter(hla_plus_buffer, 
                 x="physical_position", 
                  y="mean_painting",
                 color="zscore",
                 labels={
                     "physical_position": "Position",
                     "mean_painting": f"Mean {ancestry} ancestry (per variant)",
                     "zscore": "Z-score",
                     } ,
                 color_continuous_scale=px.colors.sequential.Inferno[::-1]  
                 )



## add mean WHG ancestry
fig.add_hline(y=df["mean_painting"].mean(), line_width=3, line_dash="dash",line_color="gray", 
              name="Genomewide mean", 
              annotation_text="Genomewide mean", 
               annotation_position="top left"
              )



# ## highlight MHC region
fig.add_vrect(
    x0=28477797, x1=33448354,
    fillcolor=px.colors.qualitative.Plotly[7],
    opacity=0.5,
    layer="below", line_width=0,
)


fig.update_xaxes(showline = True, linecolor = 'black', linewidth = 1,  mirror = True,
                 showgrid=True, gridwidth=0.3, gridcolor='lightgrey', 
                 )
fig.update_yaxes(showline = True, linecolor = 'black', linewidth = 1,  mirror = True,
                 showgrid=True, gridwidth=0.3, gridcolor='lightgrey', 
                 )


fig.update_traces(marker=dict(size=10, opacity=0.5,line=dict(width=0.19, color='black', )))


fig.update_layout(
                   height=500, width=500,
                  plot_bgcolor="white",
                  legend=dict(
                      x=0,
                      y=1,
                      bordercolor="Black",
                      borderwidth=0.8
                      )
                  )




# fig.show(renderer="browser")
# fig.show()

# for fmt in ["svg", "png", "pdf"]:
#     fig.write_image(f"{outdir}/{outname}.HLA_region.{fmt}", engine='kaleido')
# fig.write_html(f"{outdir}/{outname}.HLA_region.html")

for fmt in ["svg", "png", "pdf"]:
    fig.write_image(f"{outdir}/Fig3A.{fmt}", engine='kaleido')
fig.write_html(f"{outdir}/Fig3A.html")

#%% Calculate mean ancestry for each individual (across MHC region)
mhc_subset = df[(df["chromosome"]==6) 
         & (df["physical_position"]>=(28477797)) 
         & (df["physical_position"]<=(33448354)) 
         ]


ancestry_columns = [col for col in mhc_subset.columns if f":::{ancestry}" in col]
mean_ancestry_perindividual = []
for individual in set(col.split(':::')[0] for col in ancestry_columns):
    hap1_col = f"{individual}:::hap1:::{ancestry}"
    hap2_col = f"{individual}:::hap2:::{ancestry}"
    mean_ancestry_ind = pd.concat([mhc_subset[hap1_col], mhc_subset[hap2_col]], axis=0).mean()
    mean_ancestry_perindividual.append([individual, mean_ancestry_ind]) 


mhc_perind = pd.DataFrame(mean_ancestry_perindividual, columns=["Laboratory ID", f"Mean {ancestry} anc [MHC]"])

#%% Load table with genotpyes


genotypes = pd.read_excel(supdataf, sheet_name="DataS1", skiprows=1, header=0)
genotypes.columns = [x.replace("\n", "") for x in genotypes.columns]
genotypes["Laboratory ID"] = genotypes["Laboratory ID"].str.replace(";","_") # adjust ids to match as it is reported in fb files
lf_samples = [col.split(":::")[0] for col in df.columns if f":::hap1:::{ancestry}" in col]   # high covered LF samples used for RFmix            
genotypes = genotypes[genotypes["Laboratory ID"].isin(lf_samples)] # filter for the ones we have data
genotypes = genotypes.merge(mhc_perind, how="left", on="Laboratory ID") # merge with mhc_perind

# adjust resolution (limit to two field)
field_res = 2
def res_filter(x, res): 
    if x != "..":
        nfields = len(x.split(":")[0:res])
        test = "PASS" if nfields == res else "FAIL"
    else:
        test = "FAIL"
    return x if test == "PASS" else None


for col in [x for x in genotypes.columns if "HLA-" in x]:
    genotypes[col] = genotypes[col].apply(lambda x: ":".join(x.split(":")[0:field_res]) if x != ".." else x) 
    genotypes[col] = genotypes[col].apply(lambda x: res_filter(x, field_res) )

#%% do statistical test between carriers and non carriers of selected alleles
sig_alleles = ["A*02:01", "DRB1*08:01", "DQB1*04:02"] # significant alleles more frequent in LF (Fig. 2)

tests = []
boxplot_dfs = []
for sig_allele in sig_alleles:
    locus = sig_allele.split("*")[0]
    hla_cols = (f"HLA-{locus} (1)", f"HLA-{locus} (2)")
    genotypes_flt = genotypes[(genotypes[hla_cols[0]].notna()) & (genotypes[hla_cols[1]].notna())]

    df_present = genotypes_flt[((genotypes_flt[hla_cols[0]] == sig_allele) | (genotypes_flt[hla_cols[1]] == sig_allele)) ]
    df_present["ALLELE"] = sig_allele
    df_present["CARRIER_STATUS"] = "+"
    
    df_absent = genotypes_flt[~((genotypes_flt[f"HLA-{locus} (1)"] == sig_allele) | (genotypes_flt[f"HLA-{locus} (2)"] == sig_allele))]
    df_absent["ALLELE"] = sig_allele
    df_absent["CARRIER_STATUS"] = "-"
    
    anc_present = df_present[f"Mean {ancestry} anc [MHC]"].values
    anc_absent = df_absent[f"Mean {ancestry} anc [MHC]"].values
    if len(anc_present)>0 and len(anc_absent)>0:
        statistic, p_value = mannwhitneyu(anc_present, anc_absent, alternative='two-sided')
        tests.append([sig_allele, locus, statistic, p_value])
        boxplot_dfs.append(df_present)
        boxplot_dfs.append(df_absent)
        

tests_df = pd.DataFrame(tests, columns=["allele", "locus", "statistic", "p-value"])
reject, corrected_p_values, _, _ = smm.multipletests(tests_df["p-value"].values, alpha=0.05, method='holm-sidak') # multiple test correction
tests_df["p-value (corrected)"] = corrected_p_values

def annotate_p_value(p_value):
    if p_value <= 0.001:
        return '***'
    elif p_value <= 0.01:
        return '**'
    elif p_value <= 0.05:
        return '*'
    else:
        return 'ns'

tests_df["p-value annotation"] = tests_df["p-value (corrected)"].apply(lambda x: annotate_p_value(x) )

boxplot_df = pd.concat(boxplot_dfs)



#%% plot boxplot - Fig3 B

fig = px.box(boxplot_df, 
             x="ALLELE", 
             y=f"Mean {ancestry} anc [MHC]", 
             color="CARRIER_STATUS",
             points='all',
             hover_data = ["Laboratory ID", "Site"],
             color_discrete_map={
                 "+":px.colors.qualitative.Dark2[0],
                 "-":px.colors.qualitative.Dark2[1]
                     },
             labels={
                 f"Mean {ancestry} anc [MHC]": f"Mean {ancestry} ancestry (per sample)",
                 "ALLELE": "",
                 "CARRIER_STATUS": "Carrier<br>status"} 
             )



fig.update_xaxes(showline = True, linecolor = 'black', linewidth = 1,  mirror = True,
                 showgrid=True, gridwidth=0.3, gridcolor='lightgrey', 
                 )
fig.update_yaxes(showline = True, linecolor = 'black', linewidth = 1,  mirror = True,
                 showgrid=True, gridwidth=0.3, gridcolor='lightgrey', 
                 )


fig.update_traces(
    quartilemethod="exclusive", # or "inclusive", or "linear"
    marker=dict(
                # size=7, 
                opacity=0.5,line=dict(width=0.19, color='black', )))


for idx, row in tests_df.iterrows():
    allele = row['allele']
    p_val = row['p-value (corrected)']
    p_anno = row['p-value annotation']
    
    # x positions for the annotation line 
    x_position = boxplot_df['ALLELE'].unique().tolist().index(allele)
    x_start = x_position - 0.2
    x_end = x_position + 0.2

    # line above the box plots
    fig.add_shape(
        type="line",
        x0=x_start, y0=fig.data[0].y.max() + 0.05, 
        x1=x_end, y1=fig.data[0].y.max() + 0.05,
        line=dict(color="black", width=1)
    )

    # Add annotation with the corrected p-value
    fig.add_annotation(
        x=x_position, 
        y=fig.data[0].y.max() + 0.04,
        # text=f"p = {p_val:.3f}",
        text=f"{p_anno}",
        showarrow=False,
        yshift=10
    )
    

fig.update_layout(
                   height=500, width=500,
                  plot_bgcolor="white",
                  # legend=dict(
                  #     # x=0,
                  #     # y=1,
                  #     bordercolor="Black",
                  #     borderwidth=0.8
                  #     )
                  )




# fig.show(renderer="browser") 
# fig.show()

# for fmt in ["svg", "png", "pdf"]:
#     fig.write_image(f"{outdir}/{outname}.HLA_region.boxplot.{fmt}", engine='kaleido')
# fig.write_html(f"{outdir}/{outname}.HLA_region.boxplot.html")

for fmt in ["svg", "png", "pdf"]:
    fig.write_image(f"{outdir}/Fig3B.{fmt}", engine='kaleido')
fig.write_html(f"{outdir}/Fig3B.html")