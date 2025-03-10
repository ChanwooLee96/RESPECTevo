import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import openpyxl
from datetime import datetime
import seaborn.objects as so
from scipy import stats
############### Global variance setting ###############
Date=str(datetime.today().strftime("%Y%m%d"))[2:]
main_title='Tranformation efficiency test'
subtitle1=''
subtitle2=''
ylabel='Log(Interference efficiency)'
xlabel='Type'
sorting_var='Spacer'   ### x value in excel
value_name_insrt='Interference rate' ### y value in excel
font_label= {'size': 12,
             'weight' : 'bold',
             'color'  : 'Gray',
             'verticalalignment': 'baseline',
             'horizontalalignment': 'center'}

gray1=["#ADB7C0"]
skyblue1=["#2B9FD5"]

####### function definition: draw both swarm and box plot ######
####### x= excel data, y=sorting name of the output file ###########
def Cas3_interference(i_file_path,o_file_name,sheet="Sheet1",
                      quantative=False,width=5.5,height=3,legend=False,legend_num=6,sorting_var="Spacer",datatype="png",log=True,ymax=10**3,ymin=3*10**-1):
    df=pd.read_excel(i_file_path,sheet_name=sheet)
    File_name=f"{Date}_{o_file_name}.{datatype}"
    plt.title(main_title)
    fig, ax = plt.subplots(figsize=(width,height))
    # df[value_name_insrt]=df[value_name_insrt].apply(lambda x: 0 if x==0 else np.log10(x)) 
    color=skyblue1[0:legend_num]
    sns.set_palette(sns.color_palette(color))
    

    if quantative:
        scatter=sns.scatterplot(data=df,
                    x=sorting_var,
                    y=value_name_insrt,
                    # palette=color,
                    s=10,   
                    linewidth=0.5,
                    edgecolor="black",
                    alpha=0.7,
                    legend=False)
    else:
        swarm=sns.swarmplot(data=df,
                x=sorting_var,
                y=value_name_insrt,
                # palette=color,
                # color="black",
                dodge=False,
                size=2.5,
                linewidth=0.2,
                edgecolor="black",
                alpha=0.7)
    sns.set_palette(sns.color_palette(color))
    # plt.bar(x=df.index,height=value_name_insrt,color=color,bottom=1,ecolor="gray",data=df)
    # bar=so.Plot(df,x=sorting_var,y=value_name_insrt).add(so.Bar(baseline=ymin),so.Agg()).add(so.Range(), so.Est()).scale(y="log")
    if log==True:
        bar=so.Plot(df,x=sorting_var,y=value_name_insrt).add(so.Bar(baseline=ymin,color=skyblue1[0]),so.Agg()).add(so.Range(color="gray"), so.Est()).scale(y="log")
    if log==False:
        bar=so.Plot(df,x=sorting_var,y=value_name_insrt).add(so.Bar(baseline=ymin,color=skyblue1[0]),so.Agg()).add(so.Range(color="gray"), so.Est())
    # bar=sns.catplot(data=df,
    #             x=sorting_var,
    #             y=value_name_insrt,
    #             errcolor="gray",
    #             errwidth=0.5,
    #             palette=color,
    #             capsize=0,
    #             kind='bar'

    #             )
    # ax.axhline(y=0,color="black",alpha=0.9,linewidth=0.75)
    if legend:
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
    if log ==True:
        ax.set(yscale="log")
    ax.set_ylim(ymin,ymax)
    fig.tight_layout()
    bar.on(ax).save(File_name,dpi=300)
    # plt.savefig(File_name,dpi=300)
    print("graph done")


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
 
 
 

basicpath="RESPECTevo-Code_Rawdata/2_Plotting_Statistics/2_Cas3_interference/rawdata_zip/"

"""Endogenous Cas3 activity barplot (FigS1)"""
file_name="_Cas3 interference_endogenous_cas3.xlsx"
Cas3_interference(basicpath+file_name,"Endogenous_Cas3_expression_check",
                 sheet="Sheet1",width=3,height=2.4,datatype="pdf",
                 sorting_var="Cas3 Type")
file_name="_Cas3 interference_endogenous_cas3.xlsx"
Cas3_interference(basicpath+file_name,"Endogenous_Cas3_expression_check",
                 sheet="Sheet1",width=3,height=2.4,datatype="pdf",
                 sorting_var="Cas3 Type")

"""Endogenous cas3 statistics (Extended fig1)"""
file_name="_Cas3 interference_endogenous_cas3.xlsx"
sample_type=["WT_Cascade(+)","WT_Cascade(-)","D75A_Cascade(+)","D75A_Cascade(-)","Empty_Cascade(+)"]
for i in range(0,len(sample_type)):
    two_sample_stats_independent(basicpath+file_name,Sorting_var="Cas3 Type",value_name="Interference rate",
                                    sheet="Sheet1",sample1=sample_type[i],sample2="Empty_Cascade(-)",mod="ttest",hypothesis_mod="two-sided")


"""spacer length barplot (Extended fig6)"""
file_name="_Cas3 interference_spacer_length.xlsx"
Cas3_interference(basicpath+file_name,"Set26_length_Interference",
                 sheet="33-417nt",width=5,height=3,datatype="pdf",
                 sorting_var="Spacer",ymax=2*10**5, quantative=True)

Cas3_interference(basicpath+file_name,"Set26_length_Interference_nocrRNA",
                 sheet="nocrRNA",width=1.03,height=3,datatype="pdf",
                 sorting_var="Spacer",ymax=2*10**5, quantative=False)



"""Mismatch Effect barplot (Extended fig7)"""
file_name="_Cas3_interference_Mismatch_effect.xlsx"
Cas3_interference(basicpath+file_name,"Set27_Mismatch_Interference",
                 sheet="Sheet1",width=3.5,height=3,datatype="pdf",
                 sorting_var="Spacer",ymax=5*10**4)



"""Mismatch effect statistics (Extended fig7)""" 
file_name="_Cas3_interference_Mismatch_effect.xlsx"
sample_type=["No mismatch","mismatch+4","mismatch+34","mismatch+70","mismatch+34&70","mismatch+34to39",]
for i in range(0,len(sample_type)):
    two_sample_stats_independent(basicpath+file_name,Sorting_var="Spacer",value_name="Interference rate",
                                    sheet="Sheet1",sample1=sample_type[i],sample2="no crRNA",mod="ttest",hypothesis_mod="two-sided")
for x in ["mismatch+4","mismatch+34to39"]:
    two_sample_stats_independent(basicpath+file_name,Sorting_var="Spacer",value_name="Interference rate",
                                    sheet="Sheet1",sample1=x,sample2="No mismatch",mod="ttest",hypothesis_mod="two-sided") 