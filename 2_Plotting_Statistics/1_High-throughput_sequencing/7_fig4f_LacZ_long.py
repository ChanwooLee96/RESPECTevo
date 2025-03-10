import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from datetime import datetime
from scipy import stats
Date=str(datetime.today().strftime("%Y%m%d"))[2:]

####Global variables ####
main_title='NGS data result'
ylabel='Frequency'
ylabel2='Log(Frquency ratio)'
ylabel_ratio='Log(Frquency'
xlabel='Mutation type'
type_name2="Base"
type_name3="Spacer"


font_label= {'size': 12,
             'weight' : 'bold',
             'color'  : 'Gray',
             'verticalalignment': 'baseline',
             'horizontalalignment': 'center'}


colors6=["#EB697F","#FF8E7A","#986D81","#FFB37A","#FFB3DA","#FF5EB1"]
bluered4=["#357FC2","#D23064","#76CAF2","#E37298"]
bluered2_1=["#357FC2","#D23064"] ##deep blue, red
bluered2_2=["#76CAF2","#E37298"] ##skyblue, pale red
blue1=["#357FC2"]
blue2=["#76CAF2"]
gray1=["#ADB7C0"]
red1=["#D23064"]
red2=["#E37298"]

A_mut_list=["AtoG","AtoT","AtoC"]
G_mut_list=["GtoA","GtoC","GtoT"]
T_mut_list=["TtoA","TtoG","TtoC"]
C_mut_list=["CtoA","CtoG","CtoT"]
All_mut_list=["AtoT","AtoG","AtoC","TtoA","TtoG","TtoC","GtoA","GtoT","GtoC","CtoA","CtoT","CtoG"]
freq_limit=(10**-6)
sorting_var='Type'   ### x value column setting in xlsx file
value_name_insrt='frequency' ### y value column setting in xlsx file
type_name='Cycle'  ### hue column setting in xlsx file
ylim_min=10**-5
ylim_max=10**0
outlier=[159,162,165,168,171,201,204,207,210,213,216,219,234,543,546,549,552,553]  ### codon changed position in Phes_codon and will be eliminated 

def NGS_lineplot(mutlist,i_merging_name,merging_filelist,o_filename,
                 minmax_range=["321","562"], log=True, marker=True, error_style="band", error=lambda x: (x.min(), x.max()), linelist=[], graphrange=[], gap=25, alpha_=1,legend=True, ylim_min=2*10**-5, ylim_max=6*10**-1,
                 outlierlist=outlier,typename=type_name, date=Date, width=5.4, height=2.1, color=bluered4,marker_size=5, outtype="png"):
    """
    Lineplot generating function
        merging_filelist : should be given as list of file paths list. e.g. [[a,b,c],[d,e,f]]
                    -> a,b,c and d,e,f would be merged by group respectively.
        i_merging_name : should be given as list of names. e.g. ["1","2"]
                    -> a,b,c file and d,e,f file would be named as 1 and 2 respectivley
    """
    min_range=int(minmax_range[0]) ###setting range of analysis
    max_range=int(minmax_range[1])
    if mutlist==All_mut_list: ###setting naming for mutation list I would use
        muttype="All"
    else:
        muttype=""
        count=1
        for x in mutlist: 
            muttype=muttype+x
            count=count+1
            if count==4:
                break
    len_namelist=len(i_merging_name)     
    if len(merging_filelist)!=len_namelist:
        print("the number of given files is not matched")
        exit()
        
    #### Xposition frequency in files, would be merged in one file,with given mutation type#####
    merged_df=pd.DataFrame()
    for count in range(0,len_namelist):
        ### Melting, Merging files in merginglist######

        for filepath in merging_filelist[count]:
            df=pd.read_excel(filepath,"Sheet1")
            ## Select proper range ####
            for i in outlierlist:
                df=df[df["position"]!=i]
            df=df[df["position"]>min_range-1] 
            df=df[df["position"]<max_range+1]
            
            Mut_merged_df=pd.DataFrame()
            for mut in mutlist:
                columnlist=["position",mut]
                temp_df=df[columnlist] ### Extract given columns from dataframe ###
                temp_df=pd.melt(temp_df,id_vars="position",var_name=typename, value_name=value_name_insrt)
                ### Merging1: Merging by Mutation type  ####
                Mut_merged_df=pd.concat([temp_df,Mut_merged_df])
            Mut_merged_df=Mut_merged_df.dropna()
            Mut_merged_df.reset_index(inplace=True)

            #### Naming #####
            temp_list=[] 
            for x in Mut_merged_df.index:
                temp_list.append(i_merging_name[count])    
            Mut_merged_df["Cycle"] = temp_list
            
            
        #### Merging2: Merging by Experiment type#####
            merged_df=pd.concat([Mut_merged_df,merged_df])
    df_ref=merged_df[merged_df["Cycle"]==i_merging_name[0]]
    average_ref=df_ref["frequency"].mean()
    merged_df_noref=merged_df[merged_df["Cycle"]!=i_merging_name[0]]

    ### making graph ###
    File_name=f"{date}_{o_filename}-merged_line_{min_range}-{max_range}_{muttype}.{outtype}"
    fig, ax = plt.subplots(figsize=(width,height))
    colorlist=color[0:len(mutlist)-1]
    sns.set_palette(sns.color_palette(colorlist))
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)

            
    line=sns.lineplot(
                data=merged_df_noref,
                x="position",
                y=value_name_insrt,
                hue=typename,
                style="Cycle",
                markers=["o"],
                err_style=error_style,
                ax=ax,
                errorbar=error,
                hue_order=mutlist,
                alpha=alpha_,
                markersize=marker_size,
                linewidth=0,
                palette=color
                )
    if linelist != []:
        for i in linelist:
            ax.axvline(x=i,linestyle='dashed',color="#C2A239",alpha=1,linewidth=1)
            
    ax.axhline(y=average_ref,linestyle=':',color="black",alpha=0.9,linewidth=0.75)
    
    
    for lines in line.get_lines():

        lines.set_markeredgecolor("gray")
        lines.set_markeredgewidth('0')
        # lines.set_markerfacecolor("none")

    len_mutlist=len(mutlist)

    if legend:
        lgd = ax.legend(
                        loc='center',
                        bbox_to_anchor=(1.1,0.5),###location of legend box (x,y,width,height)
                        handletextpad=0.15,
                        frameon=True)
    else:
        ax.get_legend().remove()
    a=graphrange[0]
    b=graphrange[1]

    ax.set_xticks(range(a,b,gap))    
    ax.set_xticklabels(range(a,b,gap))
    ax.tick_params(axis='both', which='major')
    ax.set_xlabel("position",fontdict=font_label, labelpad=10)
    ax.set_ylabel(ylabel,fontdict=font_label, labelpad=10)           
    ax.set_xlim(list(map(int,minmax_range)))
    if log:
        ax.set(yscale="log")
        ax.set_ylabel("LOG("+ylabel+")",fontdict=font_label, labelpad=10)
        ax.set_ylim(ylim_min,ylim_max)
    fig.tight_layout()
    plt.savefig(File_name,dpi=600)


    print("graph done")



""" LacZ B1B0 B1B4 B1B5  """
basicpath="RESPECTevo-Code_Rawdata/2_Plotting_Statistics/1_High-throughput_sequencing/7_fig4f_LacZ_long_rawdata/"               

Wholelist_LacZ=[[["1"],["4"]], #B1B0
                [["2"],["5"]], #B1B4
                [["3"],["6"]]  #B1B5
                ]


namingdict={0:"B1B0_Cyc4",1:"B1B4_Cyc4",2:"B1B5_Cyc4"}
minmax_range_dict={0:["247","2947"],1:["247","2947"],2:["247","2947"]}
lacZrange=[250,3000]
graphrange_dict={0:lacZrange,1:lacZrange,2:lacZrange}
linelist_dict={0:[1094,1198,1466,1570],1:[1466,1570,1862,1966],2:[1466,1570,2266,2370]}
namelist=["Cycle0","Cycle4"]
for count0 in range(0,len(Wholelist_LacZ)):
    filelist_temp=Wholelist_LacZ[count0]
    for count in range(0,len(filelist_temp)):
        for count2 in range(0,len(filelist_temp[count])):
            filelist_temp[count][count2]=basicpath+f"_Exp18-{filelist_temp[count][count2]}_0.xlsx"
            
    NGS_lineplot(["CtoT","AtoG"],namelist,filelist_temp,f"Set40_{namingdict[count0]}",
                log=True,minmax_range=minmax_range_dict[count0],color=bluered2_1,graphrange=graphrange_dict[count0],linelist=linelist_dict[count0],
                typename="Mutation",marker=True,alpha_=0.9,error_style="bars",legend=False,gap=250,outlierlist=[],marker_size=3.2,outtype="pdf",width=7.7,height=2.1)
    plt.close()
