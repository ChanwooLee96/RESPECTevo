import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import openpyxl
import seaborn.objects as so
from datetime import datetime
from scipy import stats
############### Global variance setting ###############
Date=str(datetime.today().strftime("%Y%m%d"))[2:]
main_title='Suppression rate test'
subtitle1=''
subtitle2=''
ylabel='Suppression rate'
xlabel='Spacer'
"""sorting_var  ### x value in excel"""
type_name='Cycle'  ###hue name in excel
value_name_insrt='MSC (ug/mL)' ### y value in excel
font_label= {'size': 12,
             'weight' : 'bold',
             'color'  : 'Gray',
             'verticalalignment': 'baseline',
             'horizontalalignment': 'center'}
colors6=["#727171","#16803a","#bf930a","#689429"]

def suppression_freq(i_file_path,o_file_name,
                     sheet="suppression summary",width=5.5,height=3,legend_num=3,sorting_var="Spacer",type_name="Cycle",datatype="png",ymax=2,ymin=10**-7,log=True,point=False,swarm_line=0.2,alpha_=0.7):
    df=pd.read_excel(i_file_path,sheet_name=sheet)
    File_name=f"{Date}_{o_file_name}.{datatype}"
    plt.title(main_title)
    fig, ax = plt.subplots(figsize=(width,height))
    
    
    color=colors6[0:legend_num]
    sns.set_palette(sns.color_palette(color))

    
    sns.set_palette(sns.color_palette(color))
    
    if point:
        point_graph=sns.pointplot(data=df,
                        x=sorting_var,
                        y=value_name_insrt,
                        hue=type_name,
                        linestyle="None",
                        dodge=0.6,
                        marker="_",
                        markersize=15,
                        markeredgewidth=2,
                        errorbar=("ci",95),
                        err_kws={'linewidth':0.5},
                        # native_scale=False,
                        capsize=0.05,
                        palette='dark:black'
                        # palette=color
                        )
        swarm=sns.swarmplot(data=df,
                        x=sorting_var,
                        y=value_name_insrt,
                        hue=type_name,
                        palette=color,
                        dodge=True,
                        size=3.5,
                        linewidth=swarm_line,
                        edgecolor="auto",
                        legend=False,
                        alpha=alpha_)

    else:
        swarm=sns.swarmplot(data=df,
            x=sorting_var,
            y=value_name_insrt,
            hue=type_name,
            palette=color,
            # color="black",
            dodge=True,
            size=2.5,
            linewidth=swarm_line,
            edgecolor="gray",
            alpha=alpha_)    
        if log:
            bar=so.Plot(df,x=sorting_var,y=value_name_insrt,color=type_name).add(so.Bar(baseline=ymin,edgewidth=0.7),so.Agg(),so.Dodge()).scale(y="log",color=color).add(so.Range(color="gray",linewidth=0.5), so.Est(),so.Dodge())
        else:
            bar=so.Plot(df,x=sorting_var,y=value_name_insrt,color=type_name).add(so.Bar(baseline=ymin,edgewidth=0.7),so.Agg(),so.Dodge()).scale(color=color).add(so.Range(color="gray",linewidth=0.5), so.Est(),so.Dodge())
        
    # bar=sns.barplot(data=df,
    #             x=sorting_var,
    #             y=value_name_insrt,
    #             errcolor="gray",
    #             errwidth=0.5,
    #             capsize=0,
    #             hue=type_name
    #             )

    handles,labels=ax.get_legend_handles_labels()
    lgd = ax.legend(handles[0:legend_num], labels[0:legend_num],
               loc='center',
               bbox_to_anchor=(1.2,0.5),
               handletextpad=0.5,
               )
    for i in range(0,legend_num):
        lgd.legendHandles[i]._sizes = [40]


    ax.set_xlabel(xlabel,fontdict=font_label, labelpad=8)
    ax.set_ylabel(ylabel,fontdict=font_label, labelpad=10)
    if log:
        ax.set(yscale="log")
    ax.set_ylim(ymin,ymax)

    if point != True:
        bar.on(ax).save(File_name,dpi=300)
    fig.tight_layout()
    plt.savefig(File_name,dpi=300)
    print("graph done")
    plt.close()


def two_sample_stats_independent(i_file_path,Sorting_var="target",hue_var="",value_name="MIC (ug/mL)",sheet="Sheet1",sample1="FolA1+2",sample2="FolA1",mod="ttest",hue1="",hue2="",hypothesis_mod="greater",var=False):
    df=pd.read_excel(i_file_path,sheet_name=sheet)
    population1=[]
    df1=df[df[Sorting_var]==sample1]
    if hue1!="":
        df1=df1[df1[hue_var]==hue1]
    population1=df1[value_name].to_list()  ##Get only sample1's value to list
    population2=[]
    df2=df[df[Sorting_var]==sample2]
    if hue2!="":
        df2=df2[df2[hue_var]==hue2]
    # print(df1)
    # print(df2)
    population2=df2[value_name].to_list() ##Get only sample2's value to list
    if mod=="ttest":
        result=stats.ttest_ind(population1,population2,alternative=hypothesis_mod,equal_var=var)
    if mod=="mann-whitney":
        result=stats.mannwhitneyu(population1, population2, alternative=hypothesis_mod)
    print(f"Sheet name: {sheet}, Sample1: {sample1}_{hue1}, Sample2: {sample2}_{hue2} Test: {mod}")
    print(f"{sample1}_{hue1} is {hypothesis_mod} than {sample2}_{hue2}, in p-value of {result.pvalue:.5f}")
    



basicpath="RESPECTevo-Code_Rawdata/2_Plotting_Statistics/4_TMP_selection_conc/"


file_name="_Fig5_TMP_conc._by_selection.xlsx"


suppression_freq(basicpath+file_name,"Set39_Maximal_survival_passage_summary_huediff_box",
                 sheet="Summary_maximal_survival",width=6,height=3,datatype="pdf", legend_num=4,
                 sorting_var="Passage",type_name="target",ymin=0.05,ymax=3*10**3,point=True,alpha_=0.8,swarm_line=0)

sheet_type=["Sel1","Sel2","Sel3"]
sample_type=["FolA1","FolA2","W3110"]
for i in sheet_type:
    for j in sample_type:
        two_sample_stats_independent(basicpath+file_name,Sorting_var="target",value_name="MSC (ug/mL)",
                                    sheet=i,sample1="FolA1+2",sample2=j,mod="mann-whitney",hypothesis_mod="two-sided")
sample_type=["FolA1","FolA2","FolA1+2"]
for i in sheet_type:
    for j in sample_type:
        two_sample_stats_independent(basicpath+file_name,Sorting_var="target",value_name="MSC (ug/mL)",
                                    sheet=i,sample1="W3110",sample2=j,mod="mann-whitney",hypothesis_mod="two-sided")