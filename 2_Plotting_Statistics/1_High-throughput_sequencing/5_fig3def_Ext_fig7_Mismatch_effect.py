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
    colorlist=color[0:len(mutlist)]
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
    ax.set_xlim(list(map(int,graphrange)))
    if log:
        ax.set(yscale="log")
        ax.set_ylabel("LOG("+ylabel+")",fontdict=font_label, labelpad=10)
        ax.set_ylim(ylim_min,ylim_max)
    fig.tight_layout()
    plt.savefig(File_name,dpi=600)


    print("graph done")


def NGSmutation_frequency_v2(mutlist, namelist,i_filelist,o_filename,rangelist,rangemod="same",mergemod="average",outlierlist=outlier,typename=type_name,width=9,height=3,date=Date,outtype="png",
                             graphmod="Pointbar",startpoint="547",color_list=colors6,legend=True,vlinelist=[],ymin=ylim_min,ymax=ylim_max):
    """ Get mutation frequency from excel file lists

    Args:
        mutlist (list or str): Desired mutation list. e.g.[GtoA, CtoT] would give you both GtoA and CtoT mutations in data. If ="All", all possible mutation would be analyzed. 
        namelist (list): Descriptive names of each file, orders need to be matched with i_filelist. Data would be named by given list.
        i_filelist (list of list): List of File pathways. To merge multiple data at once, File path should be given in list of list. e.g. [[A],[B,C],[C,D,E]]-> B,C & C,D,E data would be merged in one data respectively
        o_filename (str): _description_
        rangelist (list): If rangemod="same", List of start point and end point of data. the number of list should be even number. e.g. [a1,a2,b1,b2] would extract a1 to a2 and b1 to b2
        rangemod (str): "same" and "different" is possible. "same" -> rangelist: list, "different -> list(list). Defaults to "same"
        mergemod (str): "raw" and "average" is possible. "raw" ->merge frequency data w/o change. "average" ->show average frequency data. Defaults to "average"
        typename (_type_, optional): _description_. Defaults to type_name.
        width (int, optional): _description_. Defaults to 9.
        height (int, optional): _description_. Defaults to 3.
        date (_type_, optional): _description_. Defaults to Date.
        outtype (str, optional): _description_. Defaults to "png".
        graphmod (str, optional): _description_. Defaults to "Pointbar".
        startpoint (str, optional): _description_. Defaults to "547".
        color_list (_type_, optional): _description_. Defaults to colors6.
        legend (bool, optional): _description_. Defaults to True.
    """

    if mutlist==All_mut_list: ###setting naming for mutation list I would use
        muttype="All"
    else:
        muttype=""
        count=1
        for x in mutlist: 
            muttype=muttype+x
            count=count+1
            if count==5:
                break
    count=0
    if rangemod=="same":
        if len(rangelist)%2!=0:
            print("the number of rangelist is not matched")
            exit()
    if rangemod=="different":
        for z in rangelist:
            if len(z)%2!=0:
                print("the number of rangelist is not matched")
                exit()
            count=count+1
    len_namelist=len(namelist)
    ### 1.frequency###
    if len(namelist)!=len(i_filelist):
        print("the number of given files is not matched")
        exit()
    
    ## e.g.) i_file list=[[a,g],[b,c],[d,e,f]] 
    count=0 ### count is address of elements of i_file_list
    merged_df=pd.DataFrame() ###DataFrame that would merge the data of a and g  in [a,g]
    merged_df_2=pd.DataFrame() ### DataFrame that would merge the data of merged_df
    average_df=pd.DataFrame() 
    for i_name in namelist: #get file name and file, the number of elements in i_file_list is same with the number of elements namelist  
        count0=0 
        df_freq_temp=pd.DataFrame()
        for y in i_filelist[count]: #e.g. until the reading of [a,g] ends 
            i_file=i_filelist[count][count0] ## reading each a and g
            df=pd.read_excel(i_file,"Sheet1")
            df_temp=pd.DataFrame()
            df0=pd.DataFrame()
            ### sorting data to have proper range
            if rangemod=="different":
                rangenum=len(rangelist[count])/2 #check the total number of min,max pair
                count1=0
                while count1<rangenum:
                    min_range=rangelist[count][2*count1]
                    max_range=rangelist[count][2*count1+1]
                    df0=df
                    for i in outlierlist:
                        df0=df0[df0["position"]!=i]
                    df0=df0[df0["position"]>min_range-1]
                    df0=df0[df0["position"]<max_range+1] # slice by its range
                    df_temp=pd.concat([df_temp,df0])
                    count1=count1+1
            if rangemod=="same":
                rangenum=len(rangelist)/2 #check the total number of min,max pair
                count1=0
                while count1<rangenum:
                    min_range=rangelist[2*count1]
                    max_range=rangelist[2*count1+1]
                    df0=df
                    for i in outlierlist:
                        df0=df0[df0["position"]!=i]
                    df0=df0[df0["position"]>min_range-1]
                    df0=df0[df0["position"]<max_range+1] # slice by its range
                    df_temp=pd.concat([df_temp,df0])
                    count1=count1+1                    
            df1=df_temp[mutlist] #get mutation type such as GtoT
            df1=pd.melt(df1, var_name=sorting_var, value_name=value_name_insrt) ## metling and setting its column name (e.g. frequency, type)
            ###adding i_name to column in melted dataframe1#####
            temp_list=[]
            for x in df1.index:
                temp_list.append(i_name)
            df1[typename]=temp_list ###by this a and g in [a,g] will have same name
            ### merging ####
            merged_df=pd.concat([df1,merged_df]) ### merge downside of the data
            df_freq_temp=pd.concat([df_freq_temp,df1["frequency"]],axis=1) ### merge only frequency column to the right side 
            count0=count0+1   ## get back to read next element in list, if first was a, then g would be read
        
        ### elements from a and g in [a,g] are already merged. Before reading next list [b,c], save the merged data of a and g 
        merged_df_2=pd.concat([merged_df_2,merged_df]) ### merged the data, After all 'for' command ended, all data like a,g,b,c,d,e,f would be merged with different naming in typename column.
        df1["frequency"]=df_freq_temp.mean(axis="columns") ## making average data of [a,g]. 
        average_df=pd.concat([average_df,df1])
        count=count+1

    if mergemod=="raw":
        merged_df=merged_df.dropna()
        merged_df=merged_df[merged_df["frequency"]>freq_limit]
    if mergemod=="average":
        average_df=average_df.dropna()
        average_df=average_df[average_df["frequency"]>freq_limit]
        merged_df=average_df
        
    ### now this dataframe contains frequency those have mutation type(e.g.GtoA) and name(e.g.cycle number) as information
    File_name=f"{date}_{o_filename}_{muttype}.{outtype}"
    
    ###### Making graph #######
    plt.title(main_title)
    fig, ax = plt.subplots(figsize=(width,height))
    color_repeat=int(len(namelist)/len(color_list))
    color=color_list*color_repeat
    sns.set_palette(sns.color_palette(color))
    if graphmod=="Box":
        box=sns.boxplot(data=merged_df,
                    x=sorting_var,
                    y=value_name_insrt,
                    hue=typename,
                    hue_order=namelist,
                    boxprops=dict(alpha=0.75),
                    showmeans=True,
                    fliersize=0,
                    meanprops={'marker': 'd', 'markeredgecolor': 'black', 'lw': 0.5,"markersize":"5","alpha":1},
                    palette=color)
    elif graphmod=="Pointbar":
        point=sns.boxplot(data=merged_df,
                    x=sorting_var,
                    y=value_name_insrt,
                    hue=typename,
                    showmeans=True,
                    meanline=True,
                    meanprops={'color': 'k', 'ls': '-', 'lw': 1},
                    medianprops={'visible': False},
                    whiskerprops={'visible': False},
                    showfliers=False,
                    showbox=False,
                    showcaps=False,
                    color="white",
                    orient="h",
                    hue_order=namelist
                    )
    strip=sns.stripplot(data=merged_df,
            x=sorting_var,
            y=value_name_insrt,
            hue=typename,
            dodge=True,
            jitter=0.25,
            alpha=0.6,
            palette=color,
            edgecolor="black",
            linewidth=0.05,
            size=3,
            hue_order=namelist
            )
    for i in vlinelist:
        ax.axhline(y=i,linestyle=':',color="black",alpha=0.9,linewidth=0.75)
    if legend:
        handles,labels=ax.get_legend_handles_labels()
        lgd = ax.legend(handles[0:len(namelist)], labels[0:len(namelist)],
                loc='center',
                bbox_to_anchor=(1.2,0.5),
                handletextpad=0.3)
        count=0
        for x in namelist:
            lgd.legendHandles[count]._sizes = [40]
            count=count+1
    else:
        ax.get_legend().remove()    
    ax.set_xlabel(xlabel,fontdict=font_label, labelpad=8)
    ax.set_ylabel(ylabel,fontdict=font_label, labelpad=10)
    ax.set(yscale="log")
    ax.set_ylim(ymin,ymax)
    fig.tight_layout()
    plt.savefig(File_name,dpi=300)
    print("graph done")


def NGSstats_v2(mutlist, namelist,i_filelist,o_filename,rangelist,rangemod="same",mergemod="average",outlierlist=outlier,typename=type_name,date=Date, ref_num=0,
                            statsmod="mann-whitney",hypothesis_mod="greater",log="False",var=False):    
        
        """ Get mutation frequency from excel file lists. 1st list of files as reference, statistics of given file can be calculated.

    Args:
        mutlist (list or str): Desired mutation list. e.g.[GtoA, CtoT] would give you both GtoA and CtoT mutations in data. If ="All", all possible mutation would be analyzed. 
        namelist (list): Descriptive names of each file, orders need to be matched with i_filelist. Data would be named by given list.
        i_filelist (list of list): List of File pathways. To merge multiple data at once, File path should be given in list of list. e.g. [[A],[B,C],[C,D,E]]-> B,C & C,D,E data would be merged in one data respectively
        o_filename (str): _description_
        rangelist (list): If rangemod="same", List of start point and end point of data. the number of list should be even number. e.g. [a1,a2,b1,b2] would extract a1 to a2 and b1 to b2
        rangemod (str): "same" and "different" is possible. "same" -> rangelist: list, "different -> list(list). Defaults to "same"
        mergemod (str): "raw" and "average" is possible. "raw" ->merge frequency data w/o change. "average" ->show average frequency data. Defaults to "average"
        typename (_type_, optional): _description_. Defaults to type_name.
    
    """
        if mutlist==All_mut_list: ###setting naming for mutation list which would be used
            muttype="All"
        else:
            muttype=""
            count=1
            for x in mutlist: 
                muttype=muttype+x
                count=count+1
                if count==5:
                    break
        count=0
        if rangemod=="same":
            if len(rangelist)%2!=0:
                print("the number of rangelist is not matched")
                exit()
        if rangemod=="different":
            for z in rangelist:
                if len(z)%2!=0:
                    print("the number of rangelist is not matched")
                    exit()
                count=count+1
        len_namelist=len(namelist)
        ### 1.frequency###
        if len(namelist)!=len(i_filelist):
            print("the number of given files is not matched")
            exit()
        ##e.g. i_file list=[[a,g],[b,c],[d,e,f]]
        count=0 ### count is address of element in filelist
        merged_df=pd.DataFrame() ### this dataframe would be used to merge [a,g]
        merged_df_2=pd.DataFrame() ### this dataframe would be used to merge result of merged_df
        average_df=pd.DataFrame() 
        for i_name in namelist: #get file name and file, len(i_file_list)=len(namelist).  
            count0=0 
            df_freq_temp=pd.DataFrame()
            for y in i_filelist[count]: #e.g. i_file list=[[a,g],[b,c],[d,e,f]] --> until [a,g] ends
                i_file=i_filelist[count][count0] ## read a in 1st "for"cycle and g in 2nd "for"cycle
                df=pd.read_excel(i_file,"Sheet1")
                df_temp=pd.DataFrame()
                df0=pd.DataFrame()
                ### Classify data then slice by its range 
                if rangemod=="different":
                    rangenum=len(rangelist[count])/2 #get the number of min max pair
                    count1=0
                    while count1<rangenum:
                        min_range=rangelist[count][2*count1]
                        max_range=rangelist[count][2*count1+1]
                        df0=df
                        for i in outlierlist:
                            df0=df0[df0["position"]!=i]
                        df0=df0[df0["position"]>min_range-1]
                        df0=df0[df0["position"]<max_range+1] # slice by its range
                        df_temp=pd.concat([df_temp,df0])
                        count1=count1+1
                if rangemod=="same":
                    rangenum=len(rangelist)/2  #get the number of min max pair
                    count1=0
                    while count1<rangenum:
                        min_range=rangelist[2*count1]
                        max_range=rangelist[2*count1+1]
                        df0=df
                        for i in outlierlist:
                            df0=df0[df0["position"]!=i]
                        df0=df0[df0["position"]>min_range-1]
                        df0=df0[df0["position"]<max_range+1] # slice by its range
                        df_temp=pd.concat([df_temp,df0])
                        count1=count1+1                    
                df1=df_temp[mutlist] #get mutation type such as GtoT
                df1=pd.melt(df1, var_name=sorting_var, value_name=value_name_insrt) ## metling and setting its column name (e.g. frequency, type)
                ###adding i_name to column in melted dataframe1#####
                temp_list=[]
                for x in df1.index:
                    temp_list.append(i_name)
                df1[typename]=temp_list ###[a,g] get the same typename column
                ### merging ####
                merged_df=pd.concat([df1,merged_df]) 
                df_freq_temp=pd.concat([df_freq_temp,df1["frequency"]],axis=1) ###gather only frequency to merge tothe right side
                count0=count0+1   ## #get back to the command then read next. if 1st is a, then next would be g
            
            ### data from a, g is already merged. save the data before starting next cycle which read [b,c]
            merged_df_2=pd.concat([merged_df_2,merged_df]) ###After all "for" commands end,  all a,g,b,c,d,e,f are gathered but with different typename column
            df1["frequency"]=df_freq_temp.mean(axis="columns") 
            average_df=pd.concat([average_df,df1])
            count=count+1

        if mergemod=="raw":
            merged_df=merged_df.dropna()
            merged_df=merged_df[merged_df["frequency"]>freq_limit]
        if mergemod=="average":
            average_df=average_df.dropna()
            average_df=average_df[average_df["frequency"]>freq_limit]
            merged_df=average_df               
        
        ref_name=namelist[ref_num]  # set given number of string in namelist as standard(ref), ref_df will be extracted from merged_df
        ref_df=merged_df[merged_df[type_name]==ref_name]
        extracted_df=merged_df[merged_df[type_name]!=ref_name]
        
        test_result_list=[]
        test_column=[f"Ref: {ref_name}","mutation type","Result_p_value" ,"Reuslt_stats"] ##set column title
        
        comparison_list=list(range(0,len_namelist))
        comparison_list.remove(ref_num) ## except given reference, other elements in list would be analyzed
        for count in comparison_list:
            i_name=namelist[count]
            temp_df=extracted_df[extracted_df[type_name]==i_name]  #set temp_df as dataframe which will be used for comparison

            for mut in mutlist:
                ref=ref_df[ref_df[sorting_var]==mut][value_name_insrt].squeeze().to_numpy()
                comparison=temp_df[temp_df[sorting_var]==mut][value_name_insrt].squeeze().to_numpy()
                if log:
                    ref=np.log(ref)
                    comparison=np.log(comparison)
                if statsmod=="ttest":
                    result=stats.ttest_ind(ref,comparison,alternative=hypothesis_mod,equal_var=var)
                if statsmod=="mann-whitney":
                    result=stats.mannwhitneyu(ref,comparison, alternative=hypothesis_mod)
                test_result_list.append([i_name,mut,result.pvalue,result.statistic]) # make column which contain [sample type, mutation type, pvalue, statistics(like t or u value)]
                
        test_result_df=pd.DataFrame(test_result_list,columns=test_column)
        test_result_df.sort_values(by=["mutation type"],inplace=True)  ## 
        ####  generating excel ####
        writer=pd.ExcelWriter(f"{date}_{statsmod}_{hypothesis_mod}_{o_filename}.xlsx",engine='xlsxwriter')
        merged_df.to_excel(writer, sheet_name="Raw data") # Raw data
        test_result_df.to_excel(writer, sheet_name="Sheet1") # test result
        writer.close()          


def NGSstats_average_mutfreq(mutlist, namelist,i_filelist,o_filename,rangelist,rangemod="same",mergemod="average",outlierlist=outlier,typename=type_name,date=Date):  ##measuring average frequency
    if mutlist==All_mut_list: ###setting naming for mutation list I would use
        muttype="All"
    else:
        muttype=""
        count=1
        for x in mutlist: 
            muttype=muttype+x
            count=count+1
            if count==5:
                break
    count=0
    if rangemod=="same":
        if len(rangelist)%2!=0:
            print("the number of rangelist is not matched")
            exit()
    if rangemod=="different":
        for z in rangelist:
            if len(z)%2!=0:
                print("the number of rangelist is not matched")
                exit()
            count=count+1
    len_namelist=len(namelist)
    ### 1.frequency###
    if len(namelist)!=len(i_filelist):
        print("the number of given files is not matched")
        exit()
    ##e.g. i_file list=[[a,g],[b,c],[d,e,f]]
    count=0 ### count is address of element in filelist
    merged_df=pd.DataFrame() ### this dataframe would be used to merge [a,g]
    merged_df_2=pd.DataFrame() ### this dataframe would be used to merge result of merged_df
    average_df=pd.DataFrame() 
    for i_name in namelist: #get file name and file, len(i_file_list)=len(namelist).  
        count0=0 
        df_freq_temp=pd.DataFrame()
        for y in i_filelist[count]: #e.g. i_file list=[[a,g],[b,c],[d,e,f]] --> until [a,g] ends
            i_file=i_filelist[count][count0] ## read a in 1st "for"cycle and g in 2nd "for"cycle
            df=pd.read_excel(i_file,"Sheet1")
            df_temp=pd.DataFrame()
            df0=pd.DataFrame()
            ### Classify data then slice by its range 
            if rangemod=="different":
                rangenum=len(rangelist[count])/2 #get the number of min max pair
                count1=0
                while count1<rangenum:
                    min_range=rangelist[count][2*count1]
                    max_range=rangelist[count][2*count1+1]
                    df0=df
                    for i in outlierlist:
                        df0=df0[df0["position"]!=i]
                    df0=df0[df0["position"]>min_range-1]
                    df0=df0[df0["position"]<max_range+1] # slice by its range
                    df_temp=pd.concat([df_temp,df0])
                    count1=count1+1
            if rangemod=="same":
                rangenum=len(rangelist)/2 #get the number of min max pair
                count1=0
                while count1<rangenum:
                    min_range=rangelist[2*count1]
                    max_range=rangelist[2*count1+1]
                    df0=df
                    for i in outlierlist:
                        df0=df0[df0["position"]!=i]
                    df0=df0[df0["position"]>min_range-1]
                    df0=df0[df0["position"]<max_range+1] # slice by its range
                    df_temp=pd.concat([df_temp,df0])
                    count1=count1+1                    
            df1=df_temp[mutlist] #get mutation type such as GtoT
            df1=pd.melt(df1, var_name=sorting_var, value_name=value_name_insrt) ## metling and setting its column name (e.g. frequency, type)
            ###adding i_name to column in melted dataframe1#####
            temp_list=[]
            for x in df1.index:
                temp_list.append(i_name)
            df1[typename]=temp_list ###[a,g] get the same typename column
            ### merging ####
            merged_df=pd.concat([df1,merged_df]) 
            df_freq_temp=pd.concat([df_freq_temp,df1["frequency"]],axis=1) ###gather only frequency to merge tothe right side
            count0=count0+1   ## #get back to the command then read next. if 1st is a, then next would be g
        
        ### data from a, g is already merged. save the data before starting next cycle which read [b,c]
        merged_df_2=pd.concat([merged_df_2,merged_df]) ###After all "for" commands end,  all a,g,b,c,d,e,f are gathered but with different typename column
        df1["frequency"]=df_freq_temp.mean(axis="columns") 
        average_df=pd.concat([average_df,df1])
        count=count+1

    if mergemod=="raw":
        merged_df=merged_df.dropna()
        merged_df=merged_df[merged_df["frequency"]>freq_limit]
    if mergemod=="average":
        average_df=average_df.dropna()
        average_df=average_df[average_df["frequency"]>freq_limit]
        merged_df=average_df
    
    #### count n number and  average mutation frequency from each sample###
    temp_column=[f"{sorting_var}",{typename},"n number" ,"Average mutation frequency", "standard deviation"] ##set column title
    temp_result_list=[]
    for mut_name in mutlist:
        for ii_name in namelist:
             ### extract frequency data having given mutation type and cycle number
            freq_array=np.array(merged_df[(merged_df[sorting_var]==mut_name)&(merged_df[typename]==ii_name)][value_name_insrt])
            n_number=len(freq_array)
            avr_freq=np.mean(freq_array)
            std_freq=np.std(freq_array, ddof=1)
            temp_result_list.append([mut_name,ii_name,n_number,avr_freq,std_freq]) # make column which contain [mutation type, sample type, n number, average, std]. final product would be list of list
    
    temp_result_df=pd.DataFrame(temp_result_list,columns=temp_column)
    # temp_result_df.sort_values(by=["mutation type"],inplace=True)
            
    ####  generating excel ####
    writer=pd.ExcelWriter(f"{date}_Averge mutation frequency_{o_filename}.xlsx",engine='xlsxwriter')
    merged_df.to_excel(writer, sheet_name="Raw data") # Raw data
    temp_result_df.to_excel(writer, sheet_name="average statistics") # average result
    writer.close()
    print("statistics done")  
    
       
basicpath="RESPECTevo-Code_Rawdata/2_Plotting_Statistics/1_High-throughput_sequencing/5_fig3def_Ext_fig7_Mismatch_effect_rawdata/"       

"""### Exp11, Set27 Mismatch ###"""
      
Wholelist_MutgRNA=[[["4_4"],["1_1","1_2","1_3"]],  ##Mis4
                   [["4_5"],["1_4","1_5","1_6"]],  ##Mis34
                   [["4_6"],["2_1","2_2","2_3"]],  ##Mis70
                   [["5_1"],["2_4","2_5","2_6"]],  ##Mis34to39
                   [["5_2"],["3_1","3_2","3_3"]],  ##Mut34_70
                   [["5_3"],["3_4","3_5","3_6"]],  ##NoMis
                   [["5_4"],["4_1","4_2","4_3"]]   ##NocrRNA
                   ]  

### point plot fig3e extended fig7###
namingdict={0:"Mis4",1:"Mis34",2:"Mis70",3:"Mis34to39",4:"Mut34_70",5:"NoMis",6:"NocrRNA"}
minmax_range_dict={0:["321","562"],1:["321","562"],2:["321","562"],3:["321","562"],4:["321","562"],5:["321","562"],6:["321","562"]}
pheSrange=[300,575]

graphrange_dict={0:pheSrange,1:pheSrange,2:pheSrange,3:pheSrange,4:pheSrange,5:pheSrange,6:pheSrange}
linelist_dict={0:[443,547,544],1:[443,547,514],2:[443,547,478],3:[443,547,514,509],4:[443,547,514,478],5:[443,547],6:[]}
namelist=["Cycle0","Cycle4"]
for count0 in range(0,len(Wholelist_MutgRNA)):
    filelist_temp=Wholelist_MutgRNA[count0]
    for count in range(0,len(filelist_temp)):
        for count2 in range(0,len(filelist_temp[count])):
            filelist_temp[count][count2]=basicpath+f"_Exp11-{filelist_temp[count][count2]}.xlsx"
            
    NGS_lineplot(["CtoT","AtoG"],namelist,filelist_temp,f"Set27_gRNA_mutation_{namingdict[count0]}",
                log=True,minmax_range=minmax_range_dict[count0],color=bluered2_1,graphrange=graphrange_dict[count0],linelist=linelist_dict[count0],
                typename="Mutation",marker=True,alpha_=0.9,error_style="bars",legend=False,outtype="pdf")
    plt.close()
    NGS_lineplot(["GtoA","TtoC"],namelist,filelist_temp,f"Set27_gRNA_mutation_{namingdict[count0]}",
                log=True,minmax_range=minmax_range_dict[count0],color=bluered2_2,graphrange=graphrange_dict[count0],linelist=linelist_dict[count0],
                typename="Mutation",marker=True,alpha_=0.9,error_style="bars",legend=False,outtype="pdf")
    plt.close()



### Boxplot targeted region fig3f extended fig7###
filelist_temp=[["3_4","3_5","3_6"],["1_1","1_2","1_3"],["1_4","1_5","1_6"],["2_1","2_2","2_3"],["3_1","3_2","3_3"],["2_4","2_5","2_6"],["4_1","4_2","4_3"]]
range__list=[[443,547],[443,547],[443,547],[443,547],[443,547],[443,547],[321,562]]
num=0
for i in filelist_temp:
    num1=0
    for j in filelist_temp[num]:
        filelist_temp[num][num1]=basicpath+f"_Exp11-{filelist_temp[num][num1]}.xlsx"
        num1=num1+1
    num=num+1   
namelist=["NoMis","Mis4","Mis34","Mis70","Mis34_70","Mis34to39","nocrRNA"]        

NGSmutation_frequency_v2(["CtoT"],namelist,filelist_temp,"Set27_mismtach_summary_average_summary_average_mutrange",range__list,width=3,graphmod="Box",color_list=blue1,legend=False,rangemod="different",outtype="pdf")
NGSmutation_frequency_v2(["GtoA"],namelist,filelist_temp,"Set27_mismtach_summary_average_summary_average_mutrange",range__list,width=3,graphmod="Box",color_list=blue2,legend=False,rangemod="different",outtype="pdf")
NGSmutation_frequency_v2(["AtoG"],namelist,filelist_temp,"Set27_mismtach_summary_average_summary_average_mutrange",range__list,width=3,graphmod="Box",color_list=red1,legend=False,rangemod="different",outtype="pdf")
NGSmutation_frequency_v2(["TtoC"],namelist,filelist_temp,"Set27_mismtach_summary_average_summary_average_mutrange",range__list,width=3,graphmod="Box",color_list=red2,legend=False,rangemod="different",outtype="pdf")
NGSstats_v2(["CtoT","GtoA","TtoC","AtoG"],namelist,filelist_temp,"Set27_mismtach_cyc4_mutrange_ref0",rangelist=range__list,rangemod="different",hypothesis_mod="two-sided",ref_num=0)
NGSstats_average_mutfreq(["CtoT","GtoA","TtoC","AtoG"],namelist,filelist_temp,"Set27_mismtach_cyc4_mutrange",rangelist=range__list,rangemod="different")


### Boxplot untargeted region extended fig7###
filelist_temp=[["3_4","3_5","3_6"],["1_1","1_2","1_3"],["1_4","1_5","1_6"],["2_1","2_2","2_3"],["3_1","3_2","3_3"],["2_4","2_5","2_6"]]

num=0
for i in filelist_temp:
    num1=0
    for j in filelist_temp[num]:
        filelist_temp[num][num1]=basicpath+f"_Exp11-{filelist_temp[num][num1]}.xlsx"
        num1=num1+1
    num=num+1   
namelist=["NoMis","Mis4","Mis34","Mis70","Mis34_70","Mis34to39"]  

NGSmutation_frequency_v2(["CtoT"],namelist,filelist_temp,"Set27_mismtach_summary_average_summary_average_mutrange_X",[321,442,548,562],width=3,graphmod="Box",color_list=blue1,legend=False,rangemod="same",outtype="pdf")
NGSmutation_frequency_v2(["GtoA"],namelist,filelist_temp,"Set27_mismtach_summary_average_summary_average_mutrange_X",[321,442,548,562],width=3,graphmod="Box",color_list=blue2,legend=False,rangemod="same",outtype="pdf")
NGSmutation_frequency_v2(["AtoG"],namelist,filelist_temp,"Set27_mismtach_summary_average_summary_average_mutrange_X",[321,442,548,562],width=3,graphmod="Box",color_list=red1,legend=False,rangemod="same",outtype="pdf")
NGSmutation_frequency_v2(["TtoC"],namelist,filelist_temp,"Set27_mismtach_summary_average_summary_average_mutrange_X",[321,442,548,562],width=3,graphmod="Box",color_list=red2,legend=False,rangemod="same",outtype="pdf")
NGSstats_average_mutfreq(["CtoT","GtoA","TtoC","AtoG"],namelist,filelist_temp,"Set27_mismtach_cyc4_mutrange_X",[321,442,548,562],rangemod="same")

