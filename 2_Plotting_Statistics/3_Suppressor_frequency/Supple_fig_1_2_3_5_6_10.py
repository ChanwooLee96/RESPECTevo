import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import openpyxl
import seaborn.objects as so
from datetime import datetime
from scipy import stats

import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
import statsmodels.api as sm
############### Global variance setting ###############
Date=str(datetime.today().strftime("%Y%m%d"))[2:]
main_title='Suppression rate test'
subtitle1=''
subtitle2=''
ylabel='Suppression_rate'
xlabel='Spacer'
"""sorting_var  """### x value in excel
type_name='Cycle'  ### hue value in excel
value_name_insrt='Suppression_rate' ### y value in excel 
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
            edgecolor="auto",
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
        lgd.legend_handles[i]._sizes = [80]


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


def two_sample_stats_independent_v2(i_file_path,namelist=[],o_filename:str="",
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
    


########### Data processing ###########
basicpath="RESPECTevo-Code_Rawdata/2_Plotting_Statistics/3_Suppressor_frequency/rawdata_zip/"


"""Supplementary fig1 Mock cycle"""
file_name="_1_Supple_fig1_Mockcycle.xlsx"
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
sample_type=["Pm_L16_CasB_33nt_bf_mock","Pm_L16_CasB_33nt_aft_mock","Pm_L16_CasB_57nt_bf_mock","Pm_L16_CasB_57nt_aft_mock"]
two_sample_stats_independent_v2(basicpath+file_name,Sorting_var="Types_mock",value_name="Suppression_rate",pairwise_num=[(0,1),(2,3)],
                                sheet="1stpanel",namelist=sample_type,statsmod="ttest",omnibus="ANOVA",o_filename="Mock_cycle_1st_log",data_transformation="log")
    
sample_type=["Pm_L16_CasB_bf_mock","Pm_L16_CasB_aft_mock","Pm_L16_CasB UNG KO_bf_mock","Pm_L16_CasB UNG KO_aft_mock"]
two_sample_stats_independent_v2(basicpath+file_name,Sorting_var="Types_mock",value_name="Suppression_rate",pairwise_num=[(0,1),(2,3)],
                                sheet="2ndpanel",namelist=sample_type,statsmod="ttest",omnibus="ANOVA",o_filename="Mock_cycle_2nd_log",data_transformation="log")

sample_type=["TadDE_L16_CasB","TadDE_L32_CasB","TadDE_L64_CasB","TadDE_L16_CasC","TadDE_L32_CasC","TadDE_L64_CasC"]
hue_type=["bf_mock","aft_mock"]
sample_hue_type=[]
for i in sample_type:
    for j in hue_type:
        sample_hue_type.append(i+"_"+j)
two_sample_stats_independent_v2(basicpath+file_name,Sorting_var="Types_mock",value_name="Suppression_rate",pairwise_num=[(0,1),(2,3),(4,5),(6,7),(8,9),(10,11)],
                                sheet="3rdpanel",namelist=sample_hue_type,statsmod="ttest",omnibus="ANOVA",o_filename="Mock_cycle_3rd_log",data_transformation="log")






"""Supplementary fig2 Pm_L16_CasB +-mutator """
file_name="_2_Supple_fig2_+-mutator_pmcda.xlsx"
suppression_freq(basicpath+file_name,"pmCasB_cascade_105bp_Cycle_supprfreq",sorting_var="Type",type_name="Cycle",
                 sheet="suppression",width=5,height=3.6,legend_num=5,datatype="pdf",ymin=5*10**-8)







"""Supplementary fig2 UNG KO"""
file_name="_3_Supple_fig2_PmCasB_C-ung-KO.xlsx"
suppression_freq(basicpath+file_name,"PmCasB_105bp_UNGKO_supprfreq",sorting_var="Type",type_name="Cycle",
                 sheet="suppression",width=5,height=3.6,legend_num=4,datatype="pdf",ymin=5*10**-8)


""" UNG KO stats"""
sample_type=["Pm-CasB 105nt_Cycle4","Pm-CasB 105nt ung_Cycle4","pm-CasB (-)crRNA_Cycle4","pm-CasB (-)crRNA ung_Cycle4"]
two_sample_stats_independent_v2(basicpath+file_name,Sorting_var="Type_Cycle",value_name="Suppression_rate",pairwise_num=[(0,1),(2,3)],
                                sheet="suppression",namelist=sample_type,statsmod="ttest",omnibus="ANOVA",o_filename="UNGKO_log",data_transformation="log")




"""Supplementary fig3 TadA-8e"""
file_name="_4_Supple_fig3_TadA_optimization.xlsx"
suppression_freq(basicpath+file_name,"TadA_optimization_cyc4",
                 sheet="suppression_Cyc4",width=5,height=3.6,legend_num=3,datatype="pdf",ymin=5*10**-8,
                 sorting_var="Types",type_name="linker")

suppression_freq(basicpath+file_name,"TadA_optimization_cyc0",
                 sheet="suppression_Cyc0",width=5,height=3.6,legend_num=3,datatype="pdf",ymin=5*10**-8,
                 sorting_var="Types",type_name="linker")





"""Supplementary fig3 PmCDA-APOBEC"""
file_name="_5_Supple_fig3_pmCDA_rAPOBEC_optimization.xlsx"
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






"""Supplementary fig3 TadDE"""
file_name="_6_Supple_fig3_TadDE_optimization.xlsx"
suppression_freq(basicpath+file_name,"Set24_TadDE_optimization_cyc4",
                 sheet="suppression_Cyc4",width=5,height=3.6,legend_num=3,datatype="pdf",ymin=5*10**-8,
                 sorting_var="Types",type_name="Linker")

suppression_freq(basicpath+file_name,"Set24_TadDE_optimization_cyc0",
                 sheet="suppression_Cyc0",width=5,height=3.6,legend_num=3,datatype="pdf",ymin=5*10**-8,
                 sorting_var="Types",type_name="Linker")





"""Supplementary fig5 TadDE proximity"""
file_name="_7_Supple_fig5_TadDE proximity.xlsx"
suppression_freq(basicpath+file_name,"Set29_Proximity_effect_Cycle_supprfreq",sorting_var="Type",type_name="Cycle",
                 sheet="suppression summary",width=4.5,height=3.6,legend_num=4,datatype="pdf",ymin=5*10**-8)




"""Supplementary fig6 spacer length"""
file_name="_8_Supple_fig6_spacer_length.xlsx"
suppression_freq(basicpath+file_name,"spacer length",
                 sheet="suppression summary",width=7,height=3.6,legend_num=4,datatype="pdf",ymin=10**-6,
                 sorting_var="Spacer",type_name="Cycle")


"""Supplementary fig10 off target effect"""

colors6=["#EB697F","#FFB37A","#FF8E7A","#986D81","#FFB3DA","#FF5EB1"]
file_name="_9_supple_fig10_off_target_PCP_rif.xlsx"
suppression_freq(basicpath+file_name,"Set42_off_target_PCPA",
                 sheet="Suppression summary_PCPA",width=5,height=3,datatype="pdf", legend_num=2,
                 sorting_var="Spacer",type_name="Cycle",ymin=10**-7,ymax=0.75)
suppression_freq(basicpath+file_name,"Set42_off_target_Rifampicin",
                 sheet="Suppression summary_Rif",width=5,height=3,datatype="pdf", legend_num=2,
                 sorting_var="Spacer",type_name="Cycle",ymin=10**-7,ymax=0.75)


sample_type=["pheS 105","LacZ B1","LacZ B2","LacZ B3","LacZ T1","RpoB 1", "RpoB 2", "Empty"]
two_sample_stats_independent_v2(basicpath+file_name,Sorting_var="Spacer",value_name="Suppression_rate",pairwise_num=[(7,0),(7,1),(7,2),(7,3),(7,4),(7,5),(7,6)],
                                sheet="PCPA_cycle4",namelist=sample_type,statsmod="ttest",omnibus="ANOVA",o_filename="offtarget_PCPA",data_transformation="log")

two_sample_stats_independent_v2(basicpath+file_name,Sorting_var="Spacer",value_name="Suppression_rate",pairwise_num=[(7,0),(7,1),(7,2),(7,3),(7,4),(7,5),(7,6)],
                                sheet="Rifampicin_cycle4",namelist=sample_type,statsmod="ttest",omnibus="ANOVA",o_filename="offtarget_Rif",data_transformation="log")

