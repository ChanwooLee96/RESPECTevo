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
colors_ten=["#727171","#16803a","#16803a","#bf930a","#bf930a","#bf930a","#689429","#689429","#689429","#689429"]


def ODgraph(i_file_path,o_file_name, symlog=True, marker=True, error_style="bars", error="ci", alpha_=0.85,legend=True, linelist=[], sheet_name="Summary",
                 huename="type",xname="TMP",yname="O.D.600", date=Date, width=5, height=2.8, color=colors_ten, marker_size=7, outtype="pdf"):
    
    File_name=f"{date}_{o_file_name}-Line_OD_growth.{outtype}"
    # plt.title(main_title)
    fig, ax = plt.subplots(figsize=(width,height))
    colorlist=color
    sns.set_palette(sns.color_palette(colorlist))
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    # sns.set_style("darkgrid")

    df=pd.read_excel(i_file_path,sheet_name)        
    line=sns.lineplot(
                data=df,
                x=xname,
                y=yname,
                hue=huename,
                style=huename,
                markers=True,
                err_style=error_style,
                ax=ax,
                dashes=False,
                errorbar=error,
                # hue_order=,
                alpha=alpha_,
                markersize=marker_size,
                linewidth=1,
                palette=color
                )
    if linelist != []:
        for i in linelist:
            ax.axvline(x=i,linestyle='dashed',color="#C2A239",alpha=1,linewidth=1)
                
    
    for lines in line.get_lines():
        lines.set_markeredgecolor("gray")
        lines.set_markeredgewidth('0')
        # lines.set_markerfacecolor("none")

    # handles,labels=ax.get_legend_handles_labels()
    if legend:
        lgd = ax.legend(
                        # handles=handles[0:len_mutlist], 
                        # labels=labels[0:len_mutlist],
                        loc='center',
                        bbox_to_anchor=(1.4,0.5),###location of legend box (x,y,width,height)
                        # borderaxespad=1.1, ##pad between legend border and axes
                        handletextpad=0.15,
                        frameon=True)
    else:
        ax.get_legend().remove()
    
    # a=graphrange[0]
    # b=graphrange[1]
    # ax.set_xticks(range(a,b,gap))    
    # ax.set_xticklabels(range(a,b,gap))
    
    ax.tick_params(axis='both', which='major')
    ax.set_xlabel("position",fontdict=font_label, labelpad=10)
    ax.set_ylabel(yname,fontdict=font_label, labelpad=10)           
    if symlog:
        ax.set(xscale="log")
        ax.set_xlabel("LOG("+xname+")",fontdict=font_label, labelpad=10)

    fig.tight_layout()
    plt.savefig(File_name,dpi=600)



"""Set39 TMP OD graph"""
basicpath="RESPECTevo-Code_Rawdata/2_Plotting_Statistics/5_Growth_curve/" 

ODgraph(basicpath+"_Fig5D_plasmid_expressed_folA_variants.xlsx","plasmid-FolA",sheet_name="Sheet1",color=colors_ten,alpha_=0.75,marker_size=7.5)