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
bluered4=["#357FC2","#D23064","#07abe3","#e35b72"]
bluered2_1=["#357FC2","#D23064"] ##deep blue, red
bluered2_2=["#07abe3","#e35b72"] ##skyblue_2, pale red_2

blue1=["#357FC2"]
blue2=["#07abe3"]
gray1=["#ADB7C0"]
red1=["#D23064"]
red2=["#e35b72"]

A_mut_list=["AtoG","AtoT","AtoC"]
G_mut_list=["GtoA","GtoC","GtoT"]
T_mut_list=["TtoA","TtoG","TtoC"]
C_mut_list=["CtoA","CtoG","CtoT"]
All_mut_list=["AtoT","AtoG","AtoC","TtoA","TtoG","TtoC","GtoA","GtoT","GtoC","CtoA","CtoT","CtoG"]
freq_limit=(10**-6)
sorting_var='Type'   ### x value column setting in xlsx file
value_name_insrt='frequency' ### y value column setting in xlsx file
value_name_insrt_ratio='frequency ratio'
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

def NGSstats_average_mutfreq(mutlist, namelist,i_filelist,o_filename,rangelist,rangemod="same",mergemod="average",outlierlist=outlier,typename=type_name,date=Date
                             , intype="xlsx"):
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
            if intype=="xlsx":
                df=pd.read_excel(i_file,"Sheet1")
            elif intype=="csv":
                df=pd.read_csv(i_file)
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
    temp_column=[f"{sorting_var}",{typename},"n number" ,"Average mutation frequency", "standard deviation","normality_pvalue","log normality_pvalue"] ##set column title
    temp_result_list=[]
    for mut_name in mutlist:
        for ii_name in namelist:
             ### extract frequency data having given mutation type and cycle number
            freq_array=np.array(merged_df[(merged_df[sorting_var]==mut_name)&(merged_df[typename]==ii_name)][value_name_insrt])
            n_number=len(freq_array)
            avr_freq=np.mean(freq_array)
            std_freq=np.std(freq_array, ddof=1)
            shapiro=stats.shapiro(freq_array) #shapiro-wilk test for normality validation in original data
            logshapiro=stats.shapiro(np.log(freq_array))  #shapiro-wilk test for normality validation in log-transformed data
            normal_pvalue=shapiro.pvalue
            lognormal_pvalue=logshapiro.pvalue
            temp_result_list.append([mut_name,ii_name,n_number,avr_freq,std_freq,normal_pvalue,lognormal_pvalue]) # make column which contain [mutation type, sample type, n number, average, std, shapiro test pvalue, log form shapiro test pvalue]. final product would be list of list
    
    temp_result_df=pd.DataFrame(temp_result_list,columns=temp_column)
    # temp_result_df.sort_values(by=["mutation type"],inplace=True)
            
    ####  generating excel ####
    writer=pd.ExcelWriter(f"{date}_Averge mutation frequency_{o_filename}.xlsx",engine='xlsxwriter')
    merged_df.to_excel(writer, sheet_name="Raw data") # Raw data
    temp_result_df.to_excel(writer, sheet_name="average statistics") # average result
    writer.close()
    print("statistics done")    
    
def NGS_lineplot_csv(mutlist,i_merging_name,merging_filelist,o_filename,
                 minmax_range=["321","562"], log=True, marker=True, error_style="band", error=lambda x: (x.min(), x.max()), linelist=[], graphrange=[], gap=25, alpha_=1,legend=True, ylim_min=2*10**-5, ylim_max=6*10**-1,
                 outlierlist=outlier,typename="Mutation", date=Date, width=5.4, height=2.1, color=bluered4,marker_size=5,include_ref=False,ratio=False,intype="xlsx", outtype="png"):
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
    ratio_df=pd.DataFrame()
    for count in range(0,len_namelist):
        ### Melting, Merging files in merginglist######

        for filepath in merging_filelist[count]:
            if intype=="xlsx":
                df=pd.read_excel(filepath,"Sheet1")
            if intype=="csv":
                df=pd.read_csv(filepath)
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
    
    if ratio:

        for mut in mutlist:
            df_ref=merged_df[(merged_df["Cycle"]==i_merging_name[0]) & (merged_df[typename]==mut)]
            df_noref=merged_df[(merged_df["Cycle"]==i_merging_name[1]) & (merged_df[typename]==mut)]
            df_noref=df_noref.rename(columns={value_name_insrt:"temp"+value_name_insrt}) 
            df_ref=df_ref.rename(columns={value_name_insrt:"ref"+value_name_insrt})
            temp_df=pd.DataFrame()        
            temp_df=pd.merge(df_ref,df_noref,on="position",how="outer")
            temp_df=temp_df.dropna(axis=0)
            temp_df=temp_df[temp_df["ref"+value_name_insrt]>0]    ####frequency가 0으로 나누는 일이 없도록 limit 이하의 값을 가진건 버림####
            temp_df=temp_df[temp_df["temp"+value_name_insrt]>0]   ####frequency의 일관성을 위해 limit 이하의 값을 가진건 버림####
            temp_df[value_name_insrt_ratio]=temp_df["temp"+value_name_insrt]/temp_df["ref"+value_name_insrt] ##divide by reference frequency        
            temp_df=temp_df.rename(columns={typename+"_x":typename}) ## 둘을 합치는 도중에 만들어진 이름 변경을 원래대로 돌려놓기
            temp_df1=temp_df[["position",typename,value_name_insrt_ratio,"Cycle_x"]]
            ratio_df=pd.concat([temp_df1,ratio_df])
    
    if not ratio:
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
    if ratio:
        print(ratio_df.head())
        print(ratio_df.tail())
        print(ratio_df.dtypes)
        line=sns.lineplot(
                data=ratio_df,
                x="position",
                y=value_name_insrt_ratio,
                hue=typename,
                style="Cycle_x",
                markers=True,
                err_style=error_style,
                ax=ax,
                hue_order=mutlist,
                alpha=alpha_,
                markersize=marker_size,
                linewidth=0,
                palette=color
                )
        for lines in line.get_lines():

            lines.set_markeredgecolor("gray")
            lines.set_markeredgewidth('0')
    if not ratio:
        if include_ref:
            ref_line=sns.lineplot(
                    data=df_ref,
                    x="position",
                    y=value_name_insrt,
                    hue=typename,
                    style="Cycle",
                    markers=["X"],
                    err_style=error_style,
                    ax=ax,
                    errorbar=error,
                    hue_order=mutlist,
                    alpha=alpha_,
                    markersize=marker_size,
                    linewidth=0,
                    palette=["k"]
                    )
            for lines in ref_line.get_lines():

                lines.set_markeredgecolor("gray")
                lines.set_markeredgewidth('0')
        
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
        for lines in line.get_lines():

            lines.set_markeredgecolor("gray")
            lines.set_markeredgewidth('0')


    if linelist != []:
        for i in linelist:
            ax.axvline(x=i,linestyle='dashed',color="#C2A239",alpha=1,linewidth=1)
    if not ratio:        
        ax.axhline(y=average_ref,linestyle=':',color="black",alpha=0.9,linewidth=0.75)
    
    




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

def NGSmutation_frequency_csv_v2(mutlist, namelist,i_filelist,o_filename,rangelist,rangemod="same",mergemod="average",outlierlist=outlier,typename=type_name,width=9,height=3,date=Date,outtype="png",
                             graphmod="Pointbar",startpoint="547",color_list=colors6,legend=True,vlinelist=[],intype="csv",stripplot=True,log=True,ymin=ylim_min,ymax=ylim_max):
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
            if intype=="xlsx":
                df=pd.read_excel(i_file,"Sheet1")
            if intype=="csv":
                df=pd.read_csv(i_file)
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
    elif graphmod=="violin":
        violin=sns.violinplot(data=merged_df,
                    x=sorting_var,
                    y=value_name_insrt,
                    hue=typename,
                    palette=color,
                    inner="box",
                    hue_order=namelist)
    if stripplot:
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
            lgd.legend_handles[count]._sizes = [40]
            count=count+1
    else:
        ax.get_legend().remove()    
    ax.set_xlabel(xlabel,fontdict=font_label, labelpad=8)
    ax.set_ylabel(ylabel,fontdict=font_label, labelpad=10)
    if log:
        ax.set(yscale="log")
    ax.set_ylim(ymin,ymax)
    fig.tight_layout()
    plt.savefig(File_name,dpi=300)
    print("graph done")
 
       
basicpath="RESPECTevo-Code_Rawdata/2_Plotting_Statistics/1_High-throughput_sequencing/4_fig3abc_Supple_fig6_7_spacer_length_rawdata/"       


"""### spacer length###"""
        
Wholelist_gRNAlength=[[["13-1_3"],["10-1_3","13-1_1","13-1_2"]],  ##33bp
                      [["13-1_6"],["10-1_4","13-1_4","13-1_5"]],  ##57bp
                      [["13-2_3"],["10-1_5","13-2_1","13-2_2"]],  ##105bp
                      [["13-2_6"],["10-1_6","13-2_4","13-2_5"]],  ##129bp
                      [["13-3_3"],["10-2_1","13-3_1","13-3_2"]],  ##153bp
                      [["13-3_6"],["10-2_2","13-3_4","13-3_5"]],  ##177bp
                      [["13-4_3"],["10-2_3","13-4_1","13-4_2"]],  ##189bp
                      [["13-4_6"],["10-2_4","13-4_4","13-4_5"]],   ##201bp
                      [["-225bp-Cyc0"],["-225bp-1","-225bp-2","-225bp-3"]], ##225bp
                      [["-309bp-Cyc0"],["-309bp-1","-309bp-2","-309bp-3"]], ##309bp
                      [["-417bp-Cyc0"],["-417bp-1","-417bp-2","-417bp-3"]], ##417bp
                      [["-pETBAD-Cyc0"],["-pETBAD-1","-pETBAD-2","-pETBAD-3"]] ### no crRNA
                      ]  


### point plot (fig3b supple fig6) ####
namingdict={0:"33bp",1:"57bp",2:"105bp",3:"129bp",4:"153bp",5:"177bp",6:"189bp",7:"201bp",8:"225bp",9:"309bp",10:"417bp",11:"Empty"}
minmax_range_dict={0:["321","562"],1:["321","562"],2:["321","562"],3:["321","562"],4:["321","562"],5:["321","562"],6:["321","562"],7:["321","562"],8:["83","562"],9:["83","562"],10:["83","562"],11:["83","562"]}
pheSrange=[300,575]
pheSrange2=[75,575]
graphrange_dict={0:pheSrange,1:pheSrange,2:pheSrange,3:pheSrange,4:pheSrange,5:pheSrange,6:pheSrange,7:pheSrange,8:pheSrange2,9:pheSrange2,10:pheSrange2,11:pheSrange2}
gap1=25
gap2=50
gap_dict={0:gap1,1:gap1,2:gap1,3:gap1,4:gap1,5:gap1,6:gap1,7:gap1,8:gap2,9:gap2,10:gap2,11:gap2}
linelist_dict={0:[515,547],1:[491,547],2:[443,547],3:[419,547],4:[395,547],5:[371,547],6:[359,547],7:[347,547],8:[323,547],9:[239,547],10:[131,547],11:[]}
namelist=["Cycle 0","Cycle 4"]
for count0 in range(0,len(Wholelist_gRNAlength)):
    filelist_temp=Wholelist_gRNAlength[count0]
    for count in range(0,len(filelist_temp)):
        for count2 in range(0,len(filelist_temp[count])):
            filelist_temp[count][count2]=basicpath+f"_Exp{filelist_temp[count][count2]}.xlsx"
            
    NGS_lineplot(["CtoT","AtoG"],namelist,filelist_temp,f"Set26_gRNA_length_{namingdict[count0]}",
                log=True,minmax_range=minmax_range_dict[count0],color=bluered2_1,graphrange=graphrange_dict[count0],linelist=linelist_dict[count0],gap=gap_dict[count0],
                typename="Mutation",marker=True,alpha_=0.9,legend=False,error_style="bars",outtype="pdf")
    plt.close()
    NGS_lineplot(["GtoA","TtoC"],namelist,filelist_temp,f"Set26_gRNA_length_{namingdict[count0]}",
                log=True,minmax_range=minmax_range_dict[count0],color=bluered2_2,graphrange=graphrange_dict[count0],linelist=linelist_dict[count0],gap=gap_dict[count0],
                typename="Mutation",marker=True,alpha_=0.9,legend=False,error_style="bars",outtype="pdf")
    plt.close()



#### targeted Boxplot fig3c supple fig6 ####
filelist_temp=[["10-1_3","13-1_1","13-1_2"],["10-1_4","13-1_4","13-1_5"],["10-1_5","13-2_1","13-2_2"],["10-1_6","13-2_4","13-2_5"],
               ["10-2_1","13-3_1","13-3_2"],["10-2_2","13-3_4","13-3_5"],["10-2_3","13-4_1","13-4_2"],["10-2_4","13-4_4","13-4_5"],
               ["-225bp-1","-225bp-2","-225bp-3"],["-309bp-1","-309bp-2","-309bp-3"],["-417bp-1","-417bp-2","-417bp-3"],["-pETBAD-1","-pETBAD-2","-pETBAD-3"]]
range__list=[[515,547],[491,547],[443,547],[419,547],[395,547],[371,547],[359,547],[347,547],[323,547],[239,547],[131,547],[83,562]]
num=0
for i in filelist_temp:
    num1=0
    for j in filelist_temp[num]:
        filelist_temp[num][num1]=basicpath+f"_Exp{filelist_temp[num][num1]}.xlsx"
        num1=num1+1
    num=num+1
namelist=["33bp","57bp","105bp","129bp","153bp","177bp","189bp","201bp","225bp","309bp","417bp","Empty"]
    
NGSmutation_frequency_v2(["CtoT"],namelist,filelist_temp,"Set26_spacerlength_cyc4_summary_average_mutrange",range__list,width=3.6,graphmod="Box",color_list=blue1,legend=False,rangemod="different",outtype="pdf")
NGSmutation_frequency_v2(["GtoA"],namelist,filelist_temp,"Set26_spacerlength_cyc4_summary_average_mutrange",range__list,width=3.6,graphmod="Box",color_list=blue2,legend=False,rangemod="different",outtype="pdf")
NGSmutation_frequency_v2(["AtoG"],namelist,filelist_temp,"Set26_spacerlength_cyc4_summary_average_mutrange",range__list,width=3.6,graphmod="Box",color_list=red1,legend=False,rangemod="different",outtype="pdf")
NGSmutation_frequency_v2(["TtoC"],namelist,filelist_temp,"Set26_spacerlength_cyc4_summary_average_mutrange",range__list,width=3.6,graphmod="Box",color_list=red2,legend=False,rangemod="different",outtype="pdf")
NGSstats_average_mutfreq(["CtoT","GtoA","TtoC","AtoG"],namelist,filelist_temp,"Set26_spacerlength_cyc4_summary_average_mutrange",rangelist=range__list,rangemod="different")



#### untargeted Boxplot fig3c supple fig6 ####
filelist_temp=[["10-1_3","13-1_1","13-1_2"],["10-1_4","13-1_4","13-1_5"],["10-1_5","13-2_1","13-2_2"],["10-1_6","13-2_4","13-2_5"],
               ["10-2_1","13-3_1","13-3_2"],["10-2_2","13-3_4","13-3_5"],["10-2_3","13-4_1","13-4_2"],["10-2_4","13-4_4","13-4_5"],
               ["-225bp-1","-225bp-2","-225bp-3"],["-309bp-1","-309bp-2","-309bp-3"],["-417bp-1","-417bp-2","-417bp-3"]]
range__list=[[321,514,548,562],[321,490,548,562],[321,442,548,562],[321,418,548,562],[321,394,548,562],[321,370,548,562],[321,358,548,562],[321,346,548,562],
             [83,322,548,562],[83,238,548,562],[83,130,548,562]]
num=0
for i in filelist_temp:
    num1=0
    for j in filelist_temp[num]:
        filelist_temp[num][num1]=basicpath+f"_Exp{filelist_temp[num][num1]}.xlsx"
        num1=num1+1
    num=num+1
namelist=["33bp","57bp","105bp","129bp","153bp","177bp","189bp","201bp","225bp","309bp","417bp"]
    
NGSmutation_frequency_v2(["CtoT"],namelist,filelist_temp,"Set26_spacerlength_cyc4_summary_average_mutrange_X",range__list,width=3.6,graphmod="Box",color_list=blue1,legend=False,rangemod="different",outtype="pdf")
NGSmutation_frequency_v2(["GtoA"],namelist,filelist_temp,"Set26_spacerlength_cyc4_summary_average_mutrange_X",range__list,width=3.6,graphmod="Box",color_list=blue2,legend=False,rangemod="different",outtype="pdf")
NGSmutation_frequency_v2(["AtoG"],namelist,filelist_temp,"Set26_spacerlength_cyc4_summary_average_mutrange_X",range__list,width=3.6,graphmod="Box",color_list=red1,legend=False,rangemod="different",outtype="pdf")
NGSmutation_frequency_v2(["TtoC"],namelist,filelist_temp,"Set26_spacerlength_cyc4_summary_average_mutrange_X",range__list,width=3.6,graphmod="Box",color_list=red2,legend=False,rangemod="different",outtype="pdf")
NGSstats_average_mutfreq(["CtoT","GtoA","TtoC","AtoG"],namelist,filelist_temp,"Set26_spacerlength_cyc4_summary_average_mutrangeX",rangelist=range__list,rangemod="different")



#### PheS site2 targeting with different length of spacer, supplementary data figure
Wholelist_pheSsite2=[[["5_4"],["1_1","1_2","1_3"]],  ##57
                    [["5_5"],["1_4","1_5","1_6"]],  ##105
                    [["5_6"],["2_1","2_2","2_3"]],  ##153 
                    [["6_1"],["2_4","2_5","2_6"]],  ##177
                    [["6_2"],["3_1","3_2","3_3"]],  ##201
                    [["6_3"],["3_4","3_5","3_6"]],  ##261
                    [["6_4"],["4_1","4_2","4_3"]],  ##309
                    [["6_5"],["4_4","4_5","4_6"]],  ##417
                    [["6_6"],["5_1","5_2","5_3"]]   ##nocrRNA
                    ]  

namingdict={0:"57nt",1:"105nt",2:"153nt",3:"177nt",4:"201nt",5:"261nt",6:"309nt",7:"417nt",8:"nocrRNA"}
minmax_range_dict={0:["-32","207"],1:["-32","207"],2:["-32","207"],3:["-32","207"],4:["-32","207"],5:["-32","207"],6:["-32","207"],7:["-32","207"],8:["-32","207"]}
pheS3range=[-50,225]

graphrange_dict={0:pheS3range,1:pheS3range,2:pheS3range,3:pheS3range,4:pheS3range,5:pheS3range,6:pheS3range,7:pheS3range,8:pheS3range}
linelist_dict={0:[151,207],1:[103,207],2:[55,207],3:[31,207],4:[7,207],5:[207],6:[207],7:[207],8:[]}
namelist=["Cycle0","Cycle4"]
for count0 in range(0,len(Wholelist_pheSsite2)):
    filelist_temp=Wholelist_pheSsite2[count0]
    for count in range(0,len(filelist_temp)):
        for count2 in range(0,len(filelist_temp[count])):
            filelist_temp[count][count2]=basicpath+f"_Exp20-{filelist_temp[count][count2]}_mutation_pos_changed.csv"
            
    NGS_lineplot_csv(["CtoT","AtoG"],namelist,filelist_temp,f"Set45_PheS_site2_{namingdict[count0]}",
                log=True,minmax_range=minmax_range_dict[count0],color=bluered2_1,graphrange=graphrange_dict[count0],linelist=linelist_dict[count0],
                typename="Mutation",marker=True,alpha_=0.9,error_style="bars",legend=False,outtype="pdf",outlierlist=[],intype="csv")
    plt.close()
    NGS_lineplot_csv(["GtoA","TtoC"],namelist,filelist_temp,f"Set45_PheS_site2_{namingdict[count0]}",
                log=True,minmax_range=minmax_range_dict[count0],color=bluered2_2,graphrange=graphrange_dict[count0],linelist=linelist_dict[count0],
                typename="Mutation",marker=True,alpha_=0.9,error_style="bars",legend=False,outtype="pdf",outlierlist=[],intype="csv")
    plt.close()    



filelist_temp=[["1_1","1_2","1_3"],["1_4","1_5","1_6"],["2_1","2_2","2_3"],["2_4","2_5","2_6"],["3_1","3_2","3_3"],["3_4","3_5","3_6"],["4_1","4_2","4_3"],["4_4","4_5","4_6"],["5_1","5_2","5_3"]]
range__list=[[151,207],[103,207],[55,207],[31,207],[7,207],[-32,207],[-32,207],[-32,207],[-32,207]]

num2=0
for j in filelist_temp:
    num1=0
    for i in filelist_temp[num2]:
        filelist_temp[num2][num1]=basicpath+f"_Exp20-{filelist_temp[num2][num1]}_mutation_pos_changed.csv"
        num1=num1+1
    num2=num2+1
    
namelist=["57nt","105nt","153nt","177nt","201nt","261nt","309nt","417nt","nocrRNA"]
NGSmutation_frequency_csv_v2(["CtoT"],namelist,filelist_temp,"Set45_pheSsite2_log_BOX",range__list,width=2.7,graphmod="Box",stripplot=True,color_list=blue1,legend=False,rangemod="different",outtype="pdf",outlierlist=[],log=True)
NGSmutation_frequency_csv_v2(["GtoA"],namelist,filelist_temp,"Set45_pheSsite2_log_BOX",range__list,width=2.7,graphmod="Box",stripplot=True,color_list=blue2,legend=False,rangemod="different",outtype="pdf",outlierlist=[],log=True)
NGSmutation_frequency_csv_v2(["AtoG"],namelist,filelist_temp,"Set45_pheSsite2_log_BOX",range__list,width=2.7,graphmod="Box",stripplot=True,color_list=red1,legend=False,rangemod="different",outtype="pdf",outlierlist=[],log=True)
NGSmutation_frequency_csv_v2(["TtoC"],namelist,filelist_temp,"Set45_pheSsite2_log_BOX",range__list,width=2.7,graphmod="Box",stripplot=True,color_list=red2,legend=False,rangemod="different",outtype="pdf",outlierlist=[],log=True)
NGSstats_average_mutfreq(["CtoT","GtoA","TtoC","AtoG"],namelist,filelist_temp,"Set45_pheSsite2",rangelist=range__list,rangemod="different",intype="csv")
   