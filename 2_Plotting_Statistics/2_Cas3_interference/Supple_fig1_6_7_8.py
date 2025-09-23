import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import openpyxl
from datetime import datetime
import seaborn.objects as so
from scipy import stats

import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
import statsmodels.api as sm
############### Global variance setting ###############
Date=str(datetime.today().strftime("%Y%m%d"))[2:]
main_title='Tranformation efficiency test'
subtitle1=''
subtitle2=''
ylabel='Log(Interference efficiency)'
xlabel='Type'
sorting_var='Spacer'   ### x value in excel
value_name_insrt='Interference_rate' ### y value in excel
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


def two_sample_stats_independent_v2(i_file_path,namelist=["FolA1","FolA2","FolA1+2","W3110"],o_filename:str="",
                                 Sorting_var:str="target",value_name:str="MSC_ugpmL",sheet:str="Sel3",
                                 omnibus:str="Permuted WTS",pairwise_num=[(),()],statsmod="permuted_brunner-munzel",post_hoc_method:str="holm-sidak",
                                 data_transformation:str="log",date=Date,equal_var=False):
    ### R package setting ###
    utils = importr("utils")
    utils.chooseCRANmirror(ind=49) ##korea cran mirror
    df=pd.read_excel(i_file_path,sheet_name=sheet)
    
    test_result_list=[]
    test_column=["group 1","group 2","g1 shapiro-wilk","g2 shapiro-wilk","Omnibus_test_type","Omnibus_pvalue","Omnibus_stats","paired_test_type","paired_p_value" ,"paired_stats"] ##set column title

    if data_transformation=="log":
        df = df.copy()  
        df.loc[:, value_name] = np.log10(df[value_name].to_numpy())
    GFD = importr("GFD")
    # Convert pandas -> R data.frame and bind into R's global env ---
    with localconverter(ro.default_converter + pandas2ri.converter):
        r_df = ro.conversion.py2rpy(df)
    ro.globalenv["dat"] = r_df
    ro.globalenv["form"] = ro.r(f'{value_name} ~ {Sorting_var}')
    permutation_number=20000

    ro.r(f'omnibus <- GFD::GFD(form, data = dat, nperm = {permutation_number})')  ##conduct omnibus test 
    if omnibus=="ANOVA":
        omnibus_result_pval = float(ro.r('as.numeric(omnibus$ATS[4])')[0])
        omnibus_result_stat = float(ro.r('as.numeric(omnibus$ATS[1])')[0])
    elif omnibus=="WTS":
        omnibus_result_pval = float(ro.r('as.numeric(omnibus$WTS[3])')[0])
        omnibus_result_stat = float(ro.r('as.numeric(omnibus$WTS[1])')[0])

    elif omnibus=="Permuted WTS":
        omnibus_result_pval = float(ro.r('as.numeric(omnibus$WTS[4])')[0])
        omnibus_result_stat = float(ro.r('as.numeric(omnibus$WTS[1])')[0])
    else:
        omnibus="N/A"

    if omnibus_result_pval==0:
        
        omnibus_result_pval=1/(permutation_number+1)
        print(f"omnibus p value is found to be 0, 1/permutation+1 would be more proper:{omnibus_result_pval} [Belinda Phipson & Gordon K Smyth, 2010]")
        
    
    for pair in pairwise_num:
        group1_name=namelist[pair[0]]
        group2_name=namelist[pair[1]]
        group1_data=df.loc[df[Sorting_var]==group1_name,value_name].squeeze().to_numpy()
        group2_data=df.loc[df[Sorting_var]==group2_name,value_name].squeeze().to_numpy()       
        shapiro_g1=stats.shapiro(group1_data).pvalue
        shapiro_g2=stats.shapiro(group2_data).pvalue

        if statsmod=="ttest":
            result=stats.ttest_ind(group1_data,group2_data,alternative="two-sided",equal_var=equal_var)
        if statsmod=="mann-whitney":
            result=stats.mannwhitneyu(group1_data,group2_data, alternative="two-sided")

        if statsmod=="permuted_brunner-munzel": # if N<10, asympotic brunner munzel is preferred
            importr("brunnermunzel")
            ro.globalenv["g1"] = FloatVector(list(group1_data))
            ro.globalenv["g2"] = FloatVector(list(group2_data))
            # Run permuted Brunner–Munzel (studentized permutation; robust to variance/shape)
            ro.r(f'bm<-brunnermunzel.test(g1,g2,perm=TRUE)')

            # Extract scalar components back into Python
            bm_stat = "N/A"
            bm_pval = float(ro.r('as.numeric(bm$p.value)')[0])
            bm_est  = float(ro.r('as.numeric(bm$estimate)')[0])  # default = P(X<Y)+.5*P(X=Y)
            class result:
                statistic= bm_stat,
                pvalue= bm_pval,
                bm_effect_estimate= bm_est
        if statsmod=="asympotic_brunner-munzel":
            importr("brunnermunzel")
            ro.globalenv["g1"] = FloatVector(list(group1_data))
            ro.globalenv["g2"] = FloatVector(list(group2_data))

            # if N>10, asympotic brunner munzel is enough
            ro.r(f'bm<-brunnermunzel.test(g1,g2)')

            # Extract scalar components back into Python

            bm_stat = float(ro.r('as.numeric(bm$statistic)')[0])
            bm_pval = float(ro.r('as.numeric(bm$p.value)')[0])
            bm_est  = float(ro.r('as.numeric(bm$estimate)')[0])  # default = P(X<Y)+.5*P(X=Y)
            class result:
                statistic= bm_stat,
                pvalue= bm_pval,
                bm_effect_estimate= bm_est
                    
        test_result_list.append([group1_name, group2_name, shapiro_g1,shapiro_g2, omnibus, omnibus_result_pval, omnibus_result_stat,statsmod,result.pvalue,result.statistic]) 
        # make column which contain ["group 1","group 2","g1 shapiro-wilk","g2 shapiro-wilk","Omnibus_test_type","Omnibus_pvalue","Omnibus_stats","paired_test_type","paired_p_value" ,"paired_stats"]
    
    

    test_result_df=pd.DataFrame(test_result_list,columns=test_column)

    test_result_df.reset_index(inplace=True)

    ### adjusted p-value ###

    pvals = np.array(list(test_result_df.loc[:, "paired_p_value"]), dtype=float).squeeze()
    _, p_corr, _, _=sm.stats.multipletests(pvals, method=post_hoc_method)
    test_result_df.loc[:, "adjusted_pval"]=p_corr

    ####  generating excel ####
    writer=pd.ExcelWriter(f"{date}_{statsmod}_{o_filename}.xlsx",engine='xlsxwriter')
    test_result_df.to_excel(writer, sheet_name="Sheet1") # test result
    writer.close()
    print("statistics done")
 
 
 

basicpath="RESPECTevo-Code_Rawdata/2_Plotting_Statistics/2_Cas3_interference/rawdata_zip/"

"""Endogenous Cas3 activity barplot (FigS1)"""
file_name="_Cas3 interference_endogenous_cas3.xlsx"
Cas3_interference(basicpath+file_name,"Endogenous_Cas3_expression_check",
                 sheet="Sheet1",width=3,height=2.4,datatype="pdf",
                 sorting_var="Cas3_Type")

"""Endogenous cas3 statistics (Supple fig1)"""
file_name="_Cas3 interference_endogenous_cas3.xlsx"
sample_type=["WT_Cascade(+)","WT_Cascade(-)","D75A_Cascade(+)","D75A_Cascade(-)","Empty_Cascade(+)","Empty_Cascade(-)"]

two_sample_stats_independent_v2(basicpath+file_name,Sorting_var="Cas3_Type",value_name="Interference_rate",pairwise_num=[(5,0),(5,1),(5,2),(5,3),(5,4)],
                                sheet="Sheet1",namelist=sample_type,statsmod="ttest",omnibus="ANOVA",o_filename="Cas3_endogenous_log",data_transformation="log")




"""spacer length barplot (Supple fig6)"""
file_name="_Cas3 interference_spacer_length.xlsx"
Cas3_interference(basicpath+file_name,"Set26_length_Interference",
                 sheet="33-417nt",width=5,height=3,datatype="pdf",
                 sorting_var="Spacer",ymax=2*10**5, quantative=True)

Cas3_interference(basicpath+file_name,"Set26_length_Interference_nocrRNA",
                 sheet="nocrRNA",width=1.03,height=3,datatype="pdf",
                 sorting_var="Spacer",ymax=2*10**5, quantative=False)

"""pheS site2 spacer length barplot (supple fig7)"""
file_name="_Cas3_interference_PheSsite2_spacer_length.xlsx"
Cas3_interference(basicpath+file_name,"Set45_PheSsite2_length_interference",
                 sheet="Sheet1",width=4.5,height=3,datatype="pdf",
                 sorting_var="Spacer",ymax=2*10**4, quantative=True)

Cas3_interference(basicpath+file_name,"Set45_PheSsite2_length_interference_nocrRNA",
                 sheet="Sheet2",width=1.1,height=3,datatype="pdf",
                 sorting_var="Spacer",ymax=2*10**4, quantative=False)


"""Mismatch Effect barplot (Supple fig8)"""
file_name="_Cas3_interference_Mismatch_effect.xlsx"
Cas3_interference(basicpath+file_name,"Set27_Mismatch_Interference",
                 sheet="Sheet1",width=3.5,height=3,datatype="pdf",
                 sorting_var="Spacer",ymax=5*10**4)



"""Mismatch effect statistics (Extended fig7)""" 
file_name="_Cas3_interference_Mismatch_effect.xlsx"
sample_type=["No mismatch","mismatch+4","mismatch+34","mismatch+70","mismatch+34&70","mismatch+34to39","no crRNA"]
two_sample_stats_independent_v2(basicpath+file_name,Sorting_var="Spacer",value_name="Interference_rate",pairwise_num=[(6,0),(6,1),(6,2),(6,3),(6,4),(6,5),(0,1),(0,5)],
                                sheet="Sheet1",namelist=sample_type,statsmod="ttest",omnibus="ANOVA",o_filename="Cas3_Mismatch_log",data_transformation="log")
