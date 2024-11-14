#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to generate the Fig1 of the manuscript. 
Code tested using python 3.10.6
"""


import pandas as pd # 1.5.1
import plotly.express as px # 5.10.0
import plotly.graph_objects as go
import plotly.io as pio
import numpy as np # 1.26.4
from scipy.spatial import ConvexHull  # 1.13.1
from scipy import stats
pio.renderers.default='svg'

# Input file is a merge of tables containing sample metadata table (.anno), 
# PCA results (evec from smartpca) and results from supervised ADMXITURE (.Q).
dfo = pd.read_csv("../data/PCA_data.tsv", sep="\t",)

# Additional file 1: Data S4 with qpadm results
datas4 = pd.read_excel("../data/AdditionalFile1.xlsx", sheet_name='DataS4', skiprows=1)

# directory to output plots 
outdir = './'

# Simplify some labels
dfo["Group ID"] = dfo["Group ID"].str.replace("_published2", "")
dfo["Group ID"] = dfo["Group ID"].str.replace("_published", "")
dfo["Group ID"] = dfo["Group ID"].str.replace("_Afanasievo.SG", "_Afanasievo")

## Remove a published outlier from the plot
dfo = dfo[~dfo["Version ID"].isin([
    "I15940_published" # 'Italy_Sardinia_C_oAfrica'
    ])].reset_index(drop=True)

# check dtypes are correct
dfo = dfo.replace("..", np.nan) 
for col in ["PC1","PC2", 'Steppe', 'WHG', 'Farmer', 'Steppe(SE)', 'WHG(SE)', 'Farmer(SE)', "Lat.", "Long."]:
    dfo[col] = dfo[col].astype(float)

for col in ["Date BP"]:
    dfo[col] = dfo[col].astype(int)
    
# invert pc1 and pc2 to look like the "standard" PCA from other publications (eg: Papac2021)
dfo["PC1"] = -dfo["PC1"]
dfo["PC2"] = -dfo["PC2"]


#%% Manage data

# Define group labels for plotting specific groups 
WHG_l = ["Italy_North_Villabruna_HG", 'Luxembourg_Loschbour_published.DG', 'Spain_ElMiron', 'Spain_HG', 'Hungary_EN_HG_Koros','Italy_Sicily_HG_OrienteC', 'Switzerland_Bichon.SG', 'Croatia_Mesolithic_HG', 'Serbia_IronGates_Mesolithic', 'Romania_IronGates_Mesolithic.SG', 'Latvia_HG',]
CHG_l = ['Georgia_Kotias.SG',   'Georgia_Satsurblia.SG', ] 
EHG_l = ["Russia_HG_Karelia",   "Russia_HG_Samara",    "Russia_Popovo_HG",]
Steppe_l = ['Russia_Samara_EBA_Yamnaya',   'Russia_Samara_EBA_Yamnaya_published',    "Russia_Samara_EBA_Yamnaya_published2",  'Russia_Afanasievo.SG',     "Russia_Afanasievo.DG",    "Russia_Afanasievo",]
Anatolia_l = ['Turkey_N', ]
pca_ref_anc_l = WHG_l + CHG_l + EHG_l + Steppe_l + Anatolia_l
pca_ref_anc = dfo[dfo["Group ID"].isin(pca_ref_anc_l)]


# focus on central european countries to plot
countries = [
       'Czech Republic', 'Hungary', 'France', 'Italy', 'Denmark', 
       'Poland','Luxembourg', 'Croatia', 'Germany', 'Belgium', 
       'Austria', 'Switzerland', 'Netherlands','Slovakia', "Slovenia",
       ]

pca = dfo[
    (dfo["Country"].isin(countries)) &
    (dfo["Date BP"]<7500) &
    (dfo["Date BP"]>4000) 
    & (~dfo["Group Label"].isin(pca_ref_anc_l))
 ]

EF_l = ["Hinkelstein", "Grossgartach", "Fellbachoeffingen", "Niederpoering"]
LF_l = ["Altendorf_PhaseI",  "Niedertiefenbach", "Rimbeck", "Warburg"]


#%% possible markers and colors

symbolVec = [
# main symbols
 'square', 'circle','cross',   'pentagon', 'star', "x", 'diamond', 'triangle-up', 'hexagram', 'hexagon', # 'octagon',

# - cross
 'square-cross', 'circle-cross', "diamond-cross",
 
# - open
 'square-open',  'circle-open', 'cross-open',   'pentagon-open', 'star-open', "x-open", 'diamond-open', 'triangle-up-open', 'hexagram-open', 'hexagon-open',
 
# -x
  'square-x', 'circle-x', 'diamond-x', 
  
# -dot
  'square-dot', 'circle-dot', 'cross-dot',   'pentagon-dot', 'star-dot', "x-dot", 'diamond-dot', 'triangle-up-dot', 'hexagram-dot', 'hexagon-dot', 
  
# -open-dot
  'square-open-dot', 'circle-open-dot', 'cross-open-dot',  'pentagon-open-dot', 'star-open-dot', "x-open-dot", 'diamond-open-dot', 'triangle-up-open-dot', 'hexagram-open-dot',  'hexagon-open-dot', 
  
# other
'diamond-tall', 'diamond-wide', 'star-triangle-up',
'diamond-tall-open', 'diamond-wide-open', 'star-triangle-down',
'diamond-tall-dot', 'diamond-wide-dot', "star-square",
'diamond-tall-open-dot', 'diamond-wide-open-dot', "star-diamond",
'triangle-left',  'triangle-right', 'triangle-down', "star-triangle-up-open",
'triangle-left-open',  'triangle-right-open', 'triangle-down-open', 
"star-triangle-down-open", "star-square-open", "star-diamond-open",
'triangle-down-dot', "star-triangle-down-dot", "star-square-dot", 
"star-diamond-dot",'triangle-down-open-dot', "star-triangle-down-open-dot", "star-square-open-dot", "star-diamond-open-dot", "arrow", "arrow-wide",
"arrow-bar-left", "arrow-bar-right", "arrow-bar-up", "arrow-bar-down", 
# "arrow-open",
  "arrow-bar-left-open", 
"arrow-bar-right-open", "arrow-bar-up-open", "arrow-bar-down-open", 
  

 ]

# only open
symbolVec_open = [x for x in symbolVec if "open" in x]

colorVec = px.colors.qualitative.Dark24

def get_markers(data_l):
    unique_strings = list(set(data_l))
    string_to_index = {string: index  for index, string in enumerate(unique_strings)}
    symbols = [ symbolVec[string_to_index[string]] for string in data_l]
    return symbols


def color_symbol_comb(l, colorVec=colorVec, symbolVec=symbolVec):
    # l = pca_BA["Group Label"].to_list()
    unique_l = sorted(list(set(l)))
    
    comb_l = []
    for c in colorVec:
        for s in symbolVec:
            comb_l.append([c,s])
    
    c_dict = {}
    s_dict = {}
    for i, pop in enumerate(unique_l):
        c_dict[pop] = comb_l[i][0]
        s_dict[pop] = comb_l[i][1]
    
    return [c_dict[pop] for pop in l], [s_dict[pop] for pop in l]


#%% prepare dfs to plot

##%% published EN  7500-6400 BP
pca_EN_range = [7500, 6400]
pca_EN = dfo[
    (dfo["Country"].isin(countries)) &
    ((dfo["Date BP"]>pca_EN_range[1]) 
    & (dfo["Date BP"]<=pca_EN_range[0]))
    & (~dfo["Group ID"].isin(pca_ref_anc_l + EF_l + LF_l))
    ]

pca_EN["group"] = f"{pca_EN_range[0]}-{pca_EN_range[1]} BP"
pca_EN["legendrank"] = 101
pca_EN["color"] = "#FFA15A"
# pca_EN["marker"] = get_markers(pca_EN["Group ID"].tolist())
pca_EN["marker2"] = "square-open"
pca_EN["color2"] = "#FFA15A"
pca_EN["size"] = 10
pca_EN["opacity"] = 0.7
pca_EN = pca_EN.sort_values("Group ID").reset_index(drop=True)


##%% published LN  6400-5000 BP
pca_LN_range = [6400, 5000]
pca_LN = dfo[
    (dfo["Country"].isin(countries)) &
    ((dfo["Date BP"]>pca_LN_range[1]) 
    & (dfo["Date BP"]<=pca_LN_range[0]))
    & (~dfo["Group ID"].isin(pca_ref_anc_l + EF_l + LF_l))
    ]

pca_LN["group"] = f"{pca_LN_range[0]}-{pca_LN_range[1]} BP"
pca_LN["legendrank"] = 102
pca_LN["color"] = "#17BECF" 
# pca_LN["marker"] = get_markers(pca_LN["Group ID"].tolist())
pca_LN["marker2"] = "circle-open"
pca_LN["color2"] = "#17BECF" 
pca_LN["size"] = 10
pca_LN["opacity"] = 0.7
pca_LN = pca_LN.sort_values("Group ID").reset_index(drop=True)


##%% published BA  4600-4000 BP
pca_BA_range = [5000, 4000]
pca_BA = dfo[
    (dfo["Country"].isin(countries)) &
    ((dfo["Date BP"]>pca_BA_range[1]) 
    & (dfo["Date BP"]<=pca_BA_range[0]))
    & (~dfo["Group ID"].isin(pca_ref_anc_l + EF_l + LF_l))
    ]

pca_BA["group"] = f"{pca_BA_range[0]}-{pca_BA_range[1]} BP"
pca_BA["legendrank"] = 103
pca_BA["color"] = "rgb(15, 133, 84)"
# pca_BA["marker"] = get_markers(pca_BA["Group ID"].tolist())
pca_BA["marker2"] = "diamond-open-dot"
pca_BA["color2"] = px.colors.qualitative.Plotly[8]
pca_BA["size"] = 10
pca_BA["opacity"] = 0.7
pca_BA = pca_BA.sort_values("Group ID").reset_index(drop=True)


##%% prepare modern ref (HO populations)

pca_ref_mod = dfo[dfo["Date BP"]==0]
pca_ref_mod["group"] = "Modern references"
pca_ref_mod["legendrank"] = 1001
pca_ref_mod["color"] = "rgb(102,102,102)"
pca_ref_mod["marker"] = get_markers(pca_ref_mod["Group ID"].tolist())
pca_ref_mod["marker2"] = "diamond-cross-open"
pca_ref_mod["color2"] = "rgb(102,102,102)"
pca_ref_mod["size"] = 9
pca_ref_mod["opacity"] = 0.3
pca_ref_mod = pca_ref_mod.sort_values("Group ID").reset_index(drop=True)

##%%  WHG 

pca_whg = dfo[dfo["Group ID"].isin(WHG_l)]
pca_whg["group"] = "Western hunter-gatherers"
pca_whg["legendrank"] = 201
pca_whg["color"] = "#990099"
# pca_whg["marker"] = get_markers(pca_whg["Group ID"].tolist())
pca_whg["marker2"] = "cross-open"
pca_whg["color2"] = "#990099"
pca_whg["size"] = 10
pca_whg["opacity"] = 0.7
pca_whg = pca_whg.sort_values("Group ID").reset_index(drop=True)

##%% prepare EHG 

pca_ehg = dfo[dfo["Group ID"].isin(EHG_l)]
pca_ehg["group"] = "Eastern hunter-gatherers"
pca_ehg["legendrank"] = 202
pca_ehg["color"] = "#8C564B"
# pca_ehg["marker"] = get_markers(pca_ehg["Group ID"].tolist())
pca_ehg["marker2"] = "x-open-dot"
pca_ehg["color2"] = "#8C564B"
pca_ehg["size"] = 10
pca_ehg["opacity"] = 0.7
pca_ehg = pca_ehg.sort_values("Group ID").reset_index(drop=True)

##%% prepare CHG 

pca_chg = dfo[dfo["Group ID"].isin(CHG_l)]
pca_chg["group"] = "Caucasus hunter-gatherers"
pca_chg["legendrank"] = 203
pca_chg["color"] = "#72B7B2"
# pca_chg["marker"] = get_markers(pca_chg["Group ID"].tolist())
pca_chg["marker2"] = "pentagon-open"
pca_chg["color2"] = "#72B7B2"
pca_chg["size"] = 10
pca_chg["opacity"] = 0.7
pca_chg = pca_chg.sort_values("Group ID").reset_index(drop=True)


##%% prepare steppe

pca_steppe = dfo[dfo["Group ID"].isin(Steppe_l)]
pca_steppe["group"] = "Steppe herders"
pca_steppe["legendrank"] = 204
pca_steppe["color"] = px.colors.qualitative.Pastel[4] 
# pca_steppe["marker"] = get_markers(pca_steppe["Group ID"].tolist())
pca_steppe["marker2"] = "star-open-dot"
pca_steppe["color2"] = px.colors.qualitative.Dark24[5] 
pca_steppe["size"] = 10
pca_steppe["opacity"] = 0.7
pca_steppe = pca_steppe.sort_values("Group ID").reset_index(drop=True)

##%% prepare anatolia

pca_anatolia = dfo[dfo["Group ID"].isin(Anatolia_l)]
pca_anatolia["group"] = "Anatolian farmers"
pca_anatolia["legendrank"] = 200
pca_anatolia["color"] = "#FECB52"
# pca_anatolia["marker"] = get_markers(pca_anatolia["Group ID"].tolist())
pca_anatolia["marker2"] = "square-x-open"
pca_anatolia["color2"] = px.colors.qualitative.Set1[6]
pca_anatolia["size"] = 10
pca_anatolia["opacity"] = 0.7
pca_anatolia = pca_anatolia.sort_values("Group ID").reset_index(drop=True)

##%% prepare EF

pca_EF = dfo[dfo["Group ID"].isin(EF_l)]
pca_EF["group"] = "Early Farmers (EF)"
pca_EF["legendrank"] = 1
pca_EF["color"] = "#EF553B"
pca_EF["color2"] = "#EF553B"
pca_EF["marker"] = get_markers(pca_EF["Group ID"].tolist())
pca_EF["marker2"] = get_markers(pca_EF["Group ID"].tolist())
pca_EF["size"] = 12
pca_EF["opacity"] = 0.9
pca_EF = pca_EF.sort_values("Group ID").reset_index(drop=True)
pca_EF["Group ID"] = pca_EF["Group ID"].map({'Hinkelstein':'Trebur/Hinkelstein', 
'Grossgartach':'Trebur/Großgartach', 
'Fellbachoeffingen':'Fellbach-Oeffingen', 
'Niederpoering':'Niederpöring',})


# add jitter to dates to improve vizualization
def add_jitter(x):
    jitter = np.random.randint(-50, 50)
    return x + jitter


pca_EF["Date BP"] = pca_EF["Date BP"].map(add_jitter)

##%% prepare LF

pca_LF = dfo[dfo["Group ID"].isin(LF_l)]
pca_LF["group"] = "Late Farmers (LF)"
pca_LF["legendrank"] = 2
pca_LF["color"] = "#636EFA"
pca_LF["color2"] = "#636EFA"
pca_LF["marker"] = get_markers(pca_LF["Group ID"].tolist())
pca_LF["marker2"] = get_markers(pca_LF["Group ID"].tolist())
pca_LF["size"] = 12
pca_LF["opacity"] = 0.9
pca_LF = pca_LF.sort_values("Group ID").reset_index(drop=True)
pca_LF["Group ID"] = pca_LF["Group ID"].map({'Altendorf_PhaseI':'Altendorf',}).fillna(pca_LF["Group ID"])


pca_LF["Date BP"] = pca_LF["Date BP"].map(add_jitter)

#%% plot pca (Fig 1A)

fig = go.Figure()

# add convex hull 
for pca_hull, hull_c in [
                [pca_EF, '#EF553B'], 
                [pca_LF, '#636EFA']
                 ]:
    
    highlighted_data = pca_hull[["PC1","PC2"]]
    
    
    # plot convex hull
    hull = ConvexHull(highlighted_data)
    hull_vertices = np.append(hull.vertices, hull.vertices[0])  # close the hull
    
    
    
    fig.add_trace(
        go.Scatter(
            x=pca_hull.loc[hull_vertices, "PC1"],
            y=pca_hull.loc[hull_vertices, "PC2"],
            mode='lines',
            line=dict(color=hull_c, width=0),
            fill='toself',
            fillcolor=hull_c,
            opacity=0.2,
            showlegend=False
        )
    )




################### add the published (EN, LN, BA), refs ( HGs, Steppe, Anatolia) + EF and LF



for df0 in [
        pca_BA,  pca_EN, pca_LN, 
        pca_whg, pca_ehg, # pca_chg, 
        pca_anatolia, pca_steppe, 
        ]:
    d = df0.copy()
    fig.add_trace(
        go.Scatter(
            x=        d["PC1"], 
            y=        d["PC2"],
            
            text = d.apply(
                lambda row: f"<i>Group Label:</i> {row['Group Label']}<br>" +
                            f"<i>Version ID:</i> {row['Version ID']}<br>" +
                            f"<i>Publication:</i> {row['Publication']}<br>" +
                            f"<i>SNPs:</i> {row['SNPs']}<br>" +
                            f"<i>Date:</i> {row['Date BP']} BP<br>" +  
                            f"<i>Steppe anc.:</i> {row['Steppe']}<br>" +
                            f"<i>WHG anc.:</i> {row['WHG']}<br>" +
                            f"<i>Farmer anc.:</i> {row['Farmer']}" 
                            ,
                axis=1
            ),
            marker=dict(
                    size=d["size"],
                    color=d["color2"],
                    symbol=d["marker2"],
                    opacity=d["opacity"],
                    ),        
            mode='markers',
            name=d.group.unique()[0],
            legendgroup="Published individuals",
            legendgrouptitle_text="Published individuals",
            legendrank=d.legendrank.unique()[0]
            )
        )
        
        
for df0 in [ pca_EF, pca_LF]:
    for g in df0["Group ID"].unique():
        d = df0[df0["Group ID"]==g]
        fig.add_trace(
            go.Scatter(
                x=        d["PC1"], 
                y=        d["PC2"],
                
                text = d.apply(
                    lambda row: f"<i>Group Label:</i> {row['Group Label']}<br>" +
                                f"<i>Version ID:</i> {row['Version ID']}<br>" +
                                f"<i>Publication:</i> {row['Publication']}<br>" +
                                f"<i>SNPs:</i> {row['SNPs']}<br>" +
                                f"<i>Date:</i> {row['Date BP']} BP<br>" +  
                                f"<i>Steppe anc.:</i> {row['Steppe']}<br>" +
                                f"<i>WHG anc.:</i> {row['WHG']}<br>" +
                                f"<i>Farmer anc.:</i> {row['Farmer']}" 
                                ,
                    axis=1
                ),
                marker=dict(
                        size=d["size"],
                        color=d["color"],
                        symbol=d["marker"],
                        opacity=d["opacity"],
                        ),        
                mode='markers',
                name=g,
                legendgroup=d.group.unique()[0],
                legendgrouptitle_text=d.group.unique()[0],
                legendrank=d.legendrank.unique()[0]
                )
            )
        

# make toggle in items of legend group    
fig.update_layout(legend=dict(groupclick="toggleitem"))

fig.update_layout(
    autosize=False,
    width=600,  # 800
    height=800, # 1300
    
    font_size=12,
    
    xaxis_title="PC1",
    yaxis_title="PC2",
    
    legend=dict(orientation="h",
                yanchor="top",
                y=-0.2,
                xanchor="left",
                x=0,
                # font=dict(size=10),
                # itemwidth=30
                ),
    
    plot_bgcolor= "white",
    paper_bgcolor= "white" ,

)

fig.update_xaxes(showline = True, linecolor = 'black', linewidth = 1,  mirror = True)
fig.update_yaxes(showline = True, linecolor = 'black', linewidth = 1,  mirror = True)

# fig.show()
# fig.show(renderer="browser")

# # save
for fmt in ["svg", "pdf"]: # 'png'
    fig.write_image(f"{outdir}/Fig1B.{fmt}", engine='kaleido')
fig.write_html(f"{outdir}/Fig1B.html")

#%% plot with modern references - Fig S1

# adjust ancient ref
pca_ancient_ref = pd.concat([
    pca_EN, pca_LN, pca_BA, 
    pca_whg, pca_ehg, pca_chg, 
    pca_anatolia, pca_steppe, 
    ])

pca_ancient_ref = pca_ancient_ref.sort_values(["Group Label"]).reset_index(drop=True)

# adjust label (to simplify a bit)
pca_ancient_ref["Group Label"] = pca_ancient_ref["Group Label"].str.replace("_published2","").str.replace("_published","")

pca_ancient_ref["color3"], pca_ancient_ref["marker3"]  = color_symbol_comb(
    l = pca_ancient_ref["Group Label"].to_list(),
    colorVec=px.colors.qualitative.Plotly[2:], 
    symbolVec=symbolVec_open
    )



pca_ancient_ref["GROUP"] = pca_ancient_ref["Group Label"]

# adjust modern ref
pca_ref_mod["GROUP"] = pca_ref_mod["Group Label"]
pca_ref_mod["color3"] = pca_ref_mod["color"]
pca_ref_mod["marker3"] = pca_ref_mod["marker"]

# adjust our ef and lf
pca_EF_LF = pd.concat([pca_EF, pca_LF])
pca_EF_LF["color3"] = pca_EF_LF["color"]
pca_EF_LF["marker3"] = pca_EF_LF["marker"]


# 1=reference with same marker: 2=references with different markers

fig = go.Figure()

###################  add convex hull 
for pca_hull, hull_c, hull_name in [
                [pca_EF, '#EF553B', "EF space"], 
                [pca_LF, '#636EFA', "LF space"],                
                 ]:
    # Filter data for the chosen cluster
    highlighted_data = pca_hull[["PC1","PC2"]]
    
    
    # Plot Convex Hull
    hull = ConvexHull(highlighted_data)
    hull_vertices = np.append(hull.vertices, hull.vertices[0])  # Close the hull    
    
    fig.add_trace(
        go.Scatter(
            x=pca_hull.loc[hull_vertices, "PC1"],
            y=pca_hull.loc[hull_vertices, "PC2"],
            mode='lines',
            line=dict(color=hull_c, width=0),
            fill='toself',
            fillcolor=hull_c,
            opacity=0.2,
            name=hull_name,
            legendgroup="Spaces",
            legendgrouptitle_text='Spaces',
            legendrank=1001,
            showlegend=True,
        )
    )


###################  add the modern reference
# plot each group with a differnt marker
for g in pca_ref_mod["Group ID"].unique():
    d = pca_ref_mod[pca_ref_mod["Group ID"]==g]
    fig.add_trace(
        go.Scatter(
            x=        d["PC1"], 
            y=        d["PC2"],
            text = d.apply(
                lambda row: f"<i>Group Label:</i> {row['Group Label']}<br>" +
                            f"<i>Version ID:</i> {row['Version ID']}<br>" +
                            f"<i>Publication:</i> {row['Publication']}<br>" +
                            f"<i>Date:</i> {row['Date BP']} BP" 
                            ,
                axis=1
            ),
            marker=dict(
                    size=d["size"],
                    color=d["color3"],
                    symbol=d["marker3"],
                    opacity=d["opacity"],
                    ),        
            mode='markers',
            name=g,
            legendgroup="Modern references",
            legendgrouptitle_text='Modern references',
            legendrank=1001,
            showlegend=True,
            )
        )


################### add the published (EN, LN, BA) refs ( HGs, Steppe, Anatolia) + EF and LF


for g in pca_ancient_ref["Group Label"].unique():
    d = pca_ancient_ref[pca_ancient_ref["Group Label"]==g]
# for g in df0["Group ID"].unique():
#     d = df0[df0["Group ID"]==g]
    fig.add_trace(
        go.Scatter(
            x=        d["PC1"], 
            y=        d["PC2"],
            
            text = d.apply(
                lambda row: f"<i>Group Label:</i> {row['Group Label']}<br>" +
                            f"<i>Version ID:</i> {row['Version ID']}<br>" +
                            f"<i>Master ID:</i> {row['Master ID']}<br>" +
                            f"<i>Publication:</i> {row['Publication']}<br>" +
                            f"<i>SNPs:</i> {row['SNPs']}<br>" +
                            f"<i>Date:</i> {row['Date BP']} BP<br>" +  
                            f"<i>Steppe anc.:</i> {row['Steppe']}<br>" +
                            f"<i>WHG anc.:</i> {row['WHG']}<br>" +
                            f"<i>Farmer anc.:</i> {row['Farmer']}" 
                            ,
                axis=1
            ),
            marker=dict(
                    size=d["size"],
                    color=d["color3"],
                    symbol=d["marker3"],
                    opacity=d["opacity"],
                    ),        
            mode='markers',
            name=g,
            legendgroup="Ancient references",
            legendgrouptitle_text="Ancient references",
            legendrank=d.legendrank.unique()[0],
            showlegend=True,
            )
        )
        
### plot EF LF

for g in pca_EF_LF["Group Label"].unique():
    d = pca_EF_LF[pca_EF_LF["Group Label"]==g]
    fig.add_trace(
        go.Scatter(
            x=        d["PC1"], 
            y=        d["PC2"],
            
            text = d.apply(
                lambda row: f"<i>Group Label:</i> {row['Group Label']}<br>" +
                            f"<i>Version ID:</i> {row['Version ID']}<br>" +
                            f"<i>Master ID:</i> {row['Master ID']}<br>" +
                            f"<i>Publication:</i> {row['Publication']}<br>" +
                            f"<i>SNPs:</i> {row['SNPs']}<br>" +
                            f"<i>Date:</i> {row['Date BP']} BP<br>" +  
                            f"<i>Steppe anc.:</i> {row['Steppe']}<br>" +
                            f"<i>WHG anc.:</i> {row['WHG']}<br>" +
                            f"<i>Farmer anc.:</i> {row['Farmer']}" 
                            ,
                axis=1
            ),
            marker=dict(
                    size=d["size"],
                    color=d["color"],
                    symbol=d["marker"],
                    opacity=d["opacity"],
                    ),        
            mode='markers',
            name=g,
            legendgroup=d.group.unique()[0],
            legendgrouptitle_text=d.group.unique()[0],
            legendrank=d.legendrank.unique()[0],
            showlegend=True,
            )
        )
    
    
# make toggle in items of legend group    
fig.update_layout(legend=dict(groupclick="toggleitem"))

fig.update_layout(
    autosize=False,
    width=1100,
    height=1600,
    
    font_size=12,
    
    xaxis_title="PC1",
    yaxis_title="PC2",
    
    legend=dict(orientation="h",
                yanchor="top",
                y=-0.10,
                xanchor="left",
                x=0,
                # font=dict(size=10),
                # itemwidth=30
                ),
    
    plot_bgcolor= "white", #'rgba(0, 0, 0, 0)',
    paper_bgcolor= "white" ,#'rgba(0, 0, 0, 0)', # background color
    
)

fig.update_xaxes(showline = True, linecolor = 'black', linewidth = 1,  mirror = True)
fig.update_yaxes(showline = True, linecolor = 'black', linewidth = 1,  mirror = True)

# fig.show()
# fig.show(renderer="browser")

# ## save
for fmt in ["svg", "pdf"]:
    fig.write_image(f"{outdir}/FigS1.{fmt}", engine='kaleido')
fig.write_html(f"{outdir}/FigS1.html")

#%% plot WHG ancestry - Fig 1C

fig = go.Figure()

################### add the published (EN, LN, BA) refs ( HGs, Steppe, Anatolia) + EF and LF

for df0 in [
        pca_EN, pca_LN, pca_BA, 

        ]:
    d = df0.copy()
    fig.add_trace(
        go.Scatter(
            x=        -d["Date BP"], 
            y=        d["WHG"],
            
            text = d.apply(
                lambda row: f"<i>Group Label:</i> {row['Group Label']}<br>" +
                            f"<i>Version ID:</i> {row['Version ID']}<br>" +
                            f"<i>Publication:</i> {row['Publication']}<br>" +
                            f"<i>SNPs:</i> {row['SNPs']}<br>" +
                            f"<i>Date:</i> {row['Date BP']} BP<br>" +  
                            f"<i>Steppe anc.:</i> {row['Steppe']}<br>" +
                            f"<i>WHG anc.:</i> {row['WHG']}<br>" +
                            f"<i>Farmer anc.:</i> {row['Farmer']}" 
                            ,
                axis=1
            ),
            marker=dict(
                    size=d["size"],
                    color=d["color2"],
                    symbol=d["marker2"],
                    opacity=d["opacity"],
                    ),        
            mode='markers',
            name=d.group.unique()[0],
            legendgroup="Published individuals",
            legendgrouptitle_text="Published individuals",
            legendrank=d.legendrank.unique()[0]
            )
        )
        
        
for df0 in [ pca_EF, pca_LF]:
    for g in df0["Group ID"].unique():
        d = df0[df0["Group ID"]==g]
        fig.add_trace(
            go.Scatter(
                x=        -d["Date BP"], 
                y=        d["WHG"],
                
                text = d.apply(
                    lambda row: f"<i>Group Label:</i> {row['Group Label']}<br>" +
                                f"<i>Version ID:</i> {row['Version ID']}<br>" +
                                f"<i>Publication:</i> {row['Publication']}<br>" +
                                f"<i>SNPs:</i> {row['SNPs']}<br>" +
                                f"<i>Date:</i> {row['Date BP']} BP<br>" +  
                                f"<i>Steppe anc.:</i> {row['Steppe']}<br>" +
                                f"<i>WHG anc.:</i> {row['WHG']}<br>" +
                                f"<i>Farmer anc.:</i> {row['Farmer']}" 
                                ,
                    axis=1
                ),
                marker=dict(
                        size=d["size"],
                        color=d["color"],
                        symbol=d["marker"],
                        opacity=d["opacity"],
                        ),        
                mode='markers',
                name=g,
                legendgroup=d.group.unique()[0],
                legendgrouptitle_text=d.group.unique()[0],
                legendrank=d.legendrank.unique()[0]
                )
            )
    

# make toggle in items of legend group    
fig.update_layout(
    autosize=False,
    width=500,  # 800
    height=550, # 1300
    
    font_size=12,
    
    xaxis_title="Date (BP)",
    yaxis_title="Western Hunter-Gatherer ancestry",
    
    legend=dict(groupclick="toggleitem",
                orientation="h",
                yanchor="top",
                y=-0.20,
                xanchor="left",
                x=0,
                # font=dict(size=10),
                # itemwidth=30
                ),
    
    plot_bgcolor= "white", #'rgba(0, 0, 0, 0)',
    paper_bgcolor= "white" ,#'rgba(0, 0, 0, 0)', # background color

)

fig.update_xaxes(showline = True, linecolor = 'black', linewidth = 1,  mirror = True)
fig.update_yaxes(showline = True, linecolor = 'black', linewidth = 1,  mirror = True)

# fig.show()
# fig.show(renderer="browser")




# ## save
for fmt in ["svg",  "pdf"]:
    fig.write_image(f"{outdir}/Fig1C.{fmt}", engine='kaleido')
fig.write_html(f"{outdir}/Fig1C.html")


    
#%% Map of samples - Fig 1A

# cropping the map
lon_min=0
lon_max=20

lat_min=45
lat_max=60

        
fig = go.Figure()

################### add the published (EN, LN, BA) refs ( HGs, Steppe, Anatolia) + EF and LF


for df0 in [
        pca_BA, pca_EN, pca_LN,  
        pca_whg, 
        ]:
    
    #filter by the lat/lon to be displayed
    df0 = df0[
        (df0["Lat."]>=lat_min) & (df0["Lat."]<=lat_max) &
        (df0["Long."]>=lon_min) & (df0["Long."]<=lon_max)
        ]
    
    
    d = df0.copy()
    
    fig.add_trace(
        go.Scattergeo(
            lat=        d["Lat."], 
            lon=        d["Long."],
            
            text = d.apply(
                lambda row: f"<i>Group Label:</i> {row['Group Label']}<br>" +
                            f"<i>Version ID:</i> {row['Version ID']}<br>" +
                            f"<i>Publication:</i> {row['Publication']}<br>" +
                            f"<i>SNPs:</i> {row['SNPs']}<br>" +
                            f"<i>Date:</i> {row['Date BP']} BP<br>" +  
                            f"<i>Steppe anc.:</i> {row['Steppe']}<br>" +
                            f"<i>WHG anc.:</i> {row['WHG']}<br>" +
                            f"<i>Farmer anc.:</i> {row['Farmer']}" 
                            ,
                axis=1
            ),
            marker=dict(
                    size=d["size"]*.9,
                    color=d["color2"],
                    symbol=d["marker2"],
                    opacity=d["opacity"],
                    ),        
            mode='markers',
            name=d.group.unique()[0],
            legendgroup="Published ancient",
            legendgrouptitle_text="Published ancient",
            legendrank=d.legendrank.unique()[0]
            )
        )
        



for df0 in [
        pca_EF, pca_LF]:
    
    #filter by the lat/lon to be displayed
    df0 = df0[
        (df0["Lat."]>=lat_min) & (df0["Lat."]<=lat_max) &
        (df0["Long."]>=lon_min) & (df0["Long."]<=lon_max)
        ]
    
    
    for g in df0["Group ID"].unique():
        d = df0[df0["Group ID"]==g]
        
        fig.add_trace(
            go.Scattergeo(
                lat=        d["Lat."], 
                lon=        d["Long."],
                
                text = d.apply(
                    lambda row: f"<i>Group Label:</i> {row['Group Label']}<br>" +
                                f"<i>Version ID:</i> {row['Version ID']}<br>" +
                                f"<i>Publication:</i> {row['Publication']}<br>" +
                                f"<i>SNPs:</i> {row['SNPs']}<br>" +
                                f"<i>Date:</i> {row['Date BP']} BP<br>" +  
                                f"<i>Steppe anc.:</i> {row['Steppe']}<br>" +
                                f"<i>WHG anc.:</i> {row['WHG']}<br>" +
                                f"<i>Farmer anc.:</i> {row['Farmer']}" 
                                ,
                    axis=1
                ),
                marker=dict(
                        size=d["size"]*.9,
                        color=d["color"],
                        symbol=d["marker"],
                        opacity=d["opacity"],
                        ),        
                mode='markers',
                name=g,
                legendgroup=d.group.unique()[0],
                legendgrouptitle_text=d.group.unique()[0],
                legendrank=d.legendrank.unique()[0]
                )
            )
        

# make toggle in items of legend group    
fig.update_layout(legend=dict(groupclick="toggleitem"))

fig.update_layout(
    autosize=False,
    width=400,  # 800
    height=900, # 1300
    
    font_size=12,
    
    
    legend=dict(orientation="h",
                yanchor="top",
                y=-0.04,
                xanchor="left",
                x=0,
                # font=dict(size=10),
                # itemwidth=30
                ),
    
    plot_bgcolor= "white", #'rgba(0, 0, 0, 0)',
    paper_bgcolor= "white" ,#'rgba(0, 0, 0, 0)', # background color
    
    geo = dict(
        showland = True, landcolor="wheat",
        showocean=True, oceancolor="lightsteelblue",
        showlakes=True, lakecolor="lightsteelblue",
        showcountries=True, countrycolor="Black", countrywidth=1,
        showsubunits=True, subunitcolor="Blue",
        showcoastlines=True, coastlinecolor="black", coastlinewidth=1,
        showrivers=True, rivercolor="lightsteelblue", riverwidth=0.5,
        showframe=True, framecolor="black", framewidth=1,
        
        
        lonaxis_range=[lon_min, lon_max],  # cropping
        lataxis_range=[lat_min, lat_max],  # cropping
        
        resolution = 50,
        
        projection = dict(
            type = 'mercator',
        ),
        )
    

)

# fig.show()
# fig.show(renderer="browser")


# # ## save
for fmt in ["svg", "pdf"]:
    fig.write_image(f"{outdir}/Fig1A.{fmt}", engine='kaleido')
fig.write_html(f"{outdir}/Fig1A.html")

#%% adjust df for plot qpadm results - Fig 1D


# adjust colum for plots
datas4 = datas4.rename({x:x.replace("\n", "_").replace("_estimated", "").replace("SE","stderr").replace(" ","") for x in datas4.columns}, axis=1)

# annotate with ef or lf
datas4.loc[datas4["target"].isin(['Fellbachoeffingen',     
                            'Niederpoering',     
                            'Trebur_Grossgartach',    
                            'Trebur_Hinkelstein',]), "group"] = "EF"

    
datas4.loc[datas4["target"].isin([    'Altendorf',
                            'Rimbeck',
                            'Warburg',
                            'Niedertiefenbach',]), "group"] = "LF"
    

datas4["target"] = datas4["target"].map({'Fellbachoeffingen':'Fellbach-Oeffingen',     
                            'Niederpoering':'Niederpöring',     
                            'Trebur_Grossgartach': 'Trebur/Großgartach',    
                            'Trebur_Hinkelstein': 'Trebur/Hinkelstein'}).fillna(datas4["target"])

# addd *2 stderr + proportion
datas4["source1_2stderr"] = abs(datas4["source1_stderr"]) * 2 + abs(datas4["source1_proportion"] )
datas4["source2_2stderr"] = abs(datas4["source2_stderr"]) * 2 + abs(datas4["source2_proportion"] )

def sum_sqr(stderrs):
    return np.sum(stderrs**2)

def weighted_stderr(sum_sqr, n_models):
    return np.sqrt(sum_sqr / n_models)


# filter two-way models
datas4_2 = datas4[
    (datas4['tail_probability'] > 0.05) &
    (datas4['feasibility'] != "infeasible") &
    (datas4['source3'] == "..") &
    ((datas4['source1_proportion'] >= 0) & (datas4['source1_proportion'] <= 1)) &
    ((datas4['source2_proportion'] >= 0) & (datas4['source2_proportion'] <= 1)) 
    ]


# get proportions averaged across feasible models
datas4_2_stats_target = datas4_2.groupby("target").agg(
    
    # mean median
    source1_proportion_mean=("source1_proportion", "mean"),
    source1_proportion_median=("source1_proportion", "median"),
    source2_proportion_mean=("source2_proportion", "mean"),
    source2_proportion_median=("source2_proportion", "median"),
    # min max
    source1_proportion_min=("source1_proportion", "min"),
    source1_proportion_max=("source1_proportion", "max"),
    source2_proportion_min=("source2_proportion", "min"),
    source2_proportion_max=("source2_proportion", "max"),
    
    # info to calculate weighted stderr
    n_models=("target", "count"),
    sum_sqr1=("source1_stderr", sum_sqr),
    sum_sqr2=("source2_stderr", sum_sqr),
    
    source1_stderr=("source1_stderr", stats.sem),
    source2_stderr=("source2_stderr", stats.sem),
    ).reset_index()



datas4_2_stats_target["source1_weighted_stderr"] = datas4_2_stats_target.apply(lambda row: 
                        weighted_stderr(row["sum_sqr1"], row["n_models"])
                        , axis=1)
    
datas4_2_stats_target["source2_weighted_stderr"] = datas4_2_stats_target.apply(lambda row: 
                        weighted_stderr(row["sum_sqr2"], row["n_models"])
                        , axis=1)

    
    
    
    
    
# get the mean across the EF and LF group
datas4_2_stats_group = datas4_2.groupby("group").agg(
    # mean median
    source1_proportion_mean=("source1_proportion", "mean"),
    source1_proportion_median=("source1_proportion", "median"),
    source2_proportion_mean=("source2_proportion", "mean"),
    source2_proportion_median=("source2_proportion", "median"),
    # min max
    source1_proportion_min=("source1_proportion", "min"),
    source1_proportion_max=("source1_proportion", "max"),
    source2_proportion_min=("source2_proportion", "min"),
    source2_proportion_max=("source2_proportion", "max")).reset_index()



#%% plot order

plot_order = [
    # EF
    'Fellbach-Oeffingen',
    'Niederpöring',
    'Trebur/Großgartach',
    'Trebur/Hinkelstein',
    
    #LF
    'Altendorf',
    'Rimbeck',
    'Warburg',
    'Niedertiefenbach',
    ]




#%% two way plots
   
fig = go.Figure(data=[

    
    go.Bar(
        name='WHG', 
           y=datas4_2_stats_target['target'].tolist(), 
           x=datas4_2_stats_target['source2_proportion_mean'].tolist(),
           base=datas4_2_stats_target['source1_proportion_mean'].tolist(),
           orientation='h',
           textposition="outside",
           textangle=0,
           marker=dict(color=px.colors.qualitative.Bold[0] ),
           ),
    
    
    go.Bar(
        name="Farmer", 
           y=datas4_2_stats_target['target'].tolist(), 
           x=datas4_2_stats_target['source1_proportion_mean'].tolist(),
            error_x=dict(type="data", array=datas4_2_stats_target['source1_stderr'].tolist()),
           orientation='h',
           marker=dict(color=px.colors.qualitative.Bold[1] ),
           ),
    
])



fig.update_xaxes(
                 showgrid=False, 
                 range=[0, 1.1],
                 dtick = 0.25,
                 )


fig.update_yaxes(showline = False,
                 showgrid=False, 
                 categoryorder='array', categoryarray= list(reversed(plot_order))
                 )

# Change the bar mode
fig.update_layout(barmode='stack', 
                  width=500,
                  height=350,
                  # margin=dict(t=100, b=100, l=100, r=50),
                  # autosize=False ,
                  plot_bgcolor="white",
                  bargap=0.03)


# fig.show(renderer="browser")
# fig.show()


# # ## save
for fmt in ["svg", "pdf"]:
    fig.write_image(f"{outdir}/Fig1D.{fmt}", engine='kaleido')
# fig.write_html(f"{outdir}/Fig1D.html")