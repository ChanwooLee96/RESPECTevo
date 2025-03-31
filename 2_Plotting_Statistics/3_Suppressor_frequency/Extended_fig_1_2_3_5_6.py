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
"""sorting_var  """### x value in excel
type_name='Cycle'  ### hue value in excel
value_name_insrt='Suppression rate' ### y value in excel 
font_label= {'size': 12,
             'weight' : 'bold',
             'color'  : 'Gray',
             'verticalalignment': 'baseline',
             'horizontalalignment': 'center'}
colors6=["#EB697F","#FF8E7A","#986D81","#FFB37A","#FFB3DA","#FF5EB1"]
colors=["#EB697F","#FF8E7A","#986D81","#FFB37A","#FFB3DA"]

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
        


    handles,labels=ax.get_legend_handles_labels()
    lgd = ax.legend(handles[0:legend_num], labels[0:legend_num],
               loc='center',
               bbox_to_anchor=(1.2,0.5),
               handletextpad=0.5,
               )
    for i in range(0,legend_num):
        lgd.legendHandles[i]._sizes = [80]


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
    population2=df2[value_name].to_list() ##Get only sample2's value to list
    if mod=="ttest":
        result=stats.ttest_ind(population1,population2,alternative=hypothesis_mod,equal_var=var)
    if mod=="mann-whitney":
        result=stats.mannwhitneyu(population1, population2, alternative=hypothesis_mod)
    print(f"Sheet name: {sheet}, Sample1: {sample1}_{hue1}, Sample2: {sample2}_{hue2} Test: {mod}")
    print(f"{sample1}_{hue1} is {hypothesis_mod} than {sample2}_{hue2}, in p-value of {result.pvalue:.5f}")
    


########### Data processing ###########
basicpath="RESPECTevo-Code_Rawdata/2_Plotting_Statistics/3_Suppressor_frequency/rawdata_zip/"


"""Extended fig1 Mock cycle"""
file_name="_1_Ext_fig1_Mockcycle.xlsx"
## 1st panel pm_16_CasB##
suppression_freq(basicpath+file_name,"PmCasB_cyc1_bf_aft_mock",
                 sheet="1stpanel",width=4,height=3.6,legend_num=2,datatype="pdf",ymin=10**-6,
                 sorting_var="Types",type_name="mock")
## 2nd panel pm_16_CasB ungKO##
suppression_freq(basicpath+file_name,"ung_cyc1_bf_aft_mock",
                 sheet="2ndpanel",width=4,height=3.6,legend_num=2,datatype="pdf",ymin=10**-6,
                 sorting_var="Types",type_name="mock")
## 2nd panel TadDE##
suppression_freq(basicpath+file_name,"TadDE_cyc1_bf_aft_mock",
                 sheet="3rdpanel",width=5,height=3.6,legend_num=2,datatype="pdf",ymin=10**-6,
                 sorting_var="Types",type_name="mock")
    
"""Mock cycle suppression stats"""
sample_type=["Pm_L16_CasB_33nt","Pm_L16_CasB_57nt"]
hue_type=["bf.mock","aft.mock"]
for i in range(0,len(sample_type)):
    two_sample_stats_independent(basicpath+file_name,Sorting_var="Types",hue_var="mock",value_name="Suppression rate",
                                    sheet="1stpanel",sample1=sample_type[i],hue1=hue_type[0],sample2=sample_type[i],hue2=hue_type[1],mod="ttest",hypothesis_mod="two-sided")
    
sample_type=["Pm_L16_CasB","Pm_L16_CasB UNG KO"]
hue_type=["bf.mock","aft.mock"]
for i in range(0,len(sample_type)):
    two_sample_stats_independent(basicpath+file_name,Sorting_var="Types",hue_var="mock",value_name="Suppression rate",
                                    sheet="2ndpanel",sample1=sample_type[i],hue1=hue_type[0],sample2=sample_type[i],hue2=hue_type[1],mod="ttest",hypothesis_mod="two-sided")  
      

sample_type=["TadDE_L16_CasB","TadDE_L32_CasB","TadDE_L64_CasB","TadDE_L16_CasC","TadDE_L32_CasC","TadDE_L64_CasC"]
hue_type=["bf.mock","aft.mock"]
for i in range(0,len(sample_type)):
    two_sample_stats_independent(basicpath+file_name,Sorting_var="Types",hue_var="mock",value_name="Suppression rate",
                                    sheet="3rdpanel",sample1=sample_type[i],hue1=hue_type[0],sample2=sample_type[i],hue2=hue_type[1],mod="ttest",hypothesis_mod="two-sided")  








"""Extended fig2 Pm_L16_CasB +-mutator """
file_name="_2_Ext_fig2_+-mutator_pmcda.xlsx"
suppression_freq(basicpath+file_name,"pmCasB_cascade_105bp_Cycle_supprfreq",sorting_var="Type",type_name="Cycle",
                 sheet="suppression",width=5,height=3.6,legend_num=5,datatype="pdf",ymin=5*10**-8)







"""Extended fig2 UNG KO"""
file_name="_3_Ext_fig2_PmCasB_C-ung-KO.xlsx"
suppression_freq(basicpath+file_name,"PmCasB_105bp_UNGKO_supprfreq",sorting_var="Type",type_name="Cycle",
                 sheet="suppression",width=5,height=3.6,legend_num=4,datatype="pdf",ymin=5*10**-8)


""" UNG KO stats"""
sample_type=["Pm-CasB 105nt","Pm-CasB 105nt ung","pm-CasB (-)crRNA","pm-CasB (-)crRNA ung"]
hue_type=["Cycle4"]
for i in [0,1]:
    two_sample_stats_independent(basicpath+file_name,Sorting_var="Type",hue_var="Cycle",value_name="Suppression rate",
                                    sheet="suppression",sample1=sample_type[2*i],hue1=hue_type[0],sample2=sample_type[2*i+1],hue2=hue_type[0],mod="ttest",hypothesis_mod="two-sided")





"""Extended fig3 TadA-8e"""
file_name="_4_Ext_fig3_TadA_optimization.xlsx"
suppression_freq(basicpath+file_name,"TadA_optimization_cyc4",
                 sheet="suppression_Cyc4",width=5,height=3.6,legend_num=3,datatype="pdf",ymin=5*10**-8,
                 sorting_var="Types",type_name="linker")

suppression_freq(basicpath+file_name,"TadA_optimization_cyc0",
                 sheet="suppression_Cyc0",width=5,height=3.6,legend_num=3,datatype="pdf",ymin=5*10**-8,
                 sorting_var="Types",type_name="linker")





"""Extended fig3 PmCDA-APOBEC"""
file_name="_5_Ext_fig3_pmCDA_rAPOBEC_optimization.xlsx"
suppression_freq(basicpath+file_name,"Pm_APO_optimization_cyc4_CasB",
                 sheet="Cyc4_CasB",width=4,height=3.6,legend_num=3,datatype="pdf",ymin=5*10**-8,
                 sorting_var="Cascade type",type_name="linker")
suppression_freq(basicpath+file_name,"Pm_APO_optimization_cyc4_CasC",
                 sheet="Cyc4_CasC",width=4,height=3.6,legend_num=3,datatype="pdf",ymin=5*10**-8,
                 sorting_var="Cascade type",type_name="linker")

suppression_freq(basicpath+file_name,"Pm_APO_optimization_cyc0_CasB",
                 sheet="Cyc0_CasB",width=4,height=3.6,legend_num=3,datatype="pdf",ymin=5*10**-8,
                 sorting_var="Cascade type",type_name="linker")
suppression_freq(basicpath+file_name,"Pm_APO_optimization_cyc0_CasC",
                 sheet="Cyc0_CasC",width=4,height=3.6,legend_num=3,datatype="pdf",ymin=5*10**-8,
                 sorting_var="Cascade type",type_name="linker")






"""Extended fig3 TadDE"""
file_name="_6_Ext_fig3_TadDE_optimization.xlsx"
suppression_freq(basicpath+file_name,"Set24_TadDE_optimization_cyc4",
                 sheet="suppression_Cyc4",width=5,height=3.6,legend_num=3,datatype="pdf",ymin=5*10**-8,
                 sorting_var="Types",type_name="Linker")

suppression_freq(basicpath+file_name,"Set24_TadDE_optimization_cyc0",
                 sheet="suppression_Cyc0",width=5,height=3.6,legend_num=3,datatype="pdf",ymin=5*10**-8,
                 sorting_var="Types",type_name="Linker")





"""Extended fig5 TadDE proximity"""
file_name="_7_Ext_fig5_TadDE proximity.xlsx"
suppression_freq(basicpath+file_name,"Set29_Proximity_effect_Cycle_supprfreq",sorting_var="Type",type_name="Cycle",
                 sheet="suppression summary",width=4.5,height=3.6,legend_num=4,datatype="pdf",ymin=5*10**-8)


"""Extended fig5 TadDE proximity """
colors6=["#EB697F","#FFB37A","#FF8E7A","#986D81","#FFB3DA","#FF5EB1"]
file_name="_8_Ext_fig5_RESPECTevo_rifampicin.xlsx"
suppression_freq(basicpath+file_name,"_RESPECTevo_rifampicin_Cycle_supprfreq",sorting_var="Type",type_name="Cycle",
                 sheet="suppression summary",width=3,height=3.6,legend_num=2,datatype="pdf",ymin=5*10**-8)


sample_type=["TadDE_L16_casB 105nt","TadDE_L16_casB (-)crRNA"]
hue_type=["Cycle4","Cycle0"]
for i in [0,1]:
    two_sample_stats_independent(basicpath+file_name,Sorting_var="Type",hue_var="Cycle",value_name="Suppression rate",
                                    sheet="suppression summary",sample1=sample_type[i],hue1=hue_type[0],sample2=sample_type[i],hue2=hue_type[1],mod="ttest",hypothesis_mod="two-sided")


"""Extended fig6 spacer length"""
file_name="_9_Ext_fig6_spacer_length.xlsx"
suppression_freq(basicpath+file_name,"spacer length",
                 sheet="suppression summary",width=7,height=3.6,legend_num=4,datatype="pdf",ymin=10**-6,
                 sorting_var="Spacer",type_name="Cycle")





