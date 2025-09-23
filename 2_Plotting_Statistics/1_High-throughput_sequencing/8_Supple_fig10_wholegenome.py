import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from datetime import datetime
from scipy import stats
import matplotlib as mpl
mpl.rcParams['agg.path.chunksize'] = 1000000
import cProfile

Date=str(datetime.today().strftime("%Y%m%d"))[2:]

A_mut_list=["AtoG","AtoT","AtoC"]
G_mut_list=["GtoA","GtoC","GtoT"]
T_mut_list=["TtoA","TtoG","TtoC"]
C_mut_list=["CtoA","CtoG","CtoT"]
All_mut_list=["AtoT","AtoG","AtoC","TtoA","TtoG","TtoC","GtoA","GtoT","GtoC","CtoA","CtoT","CtoG"]
freq_limit=(10**-6)
freq_limit=(-10**6)

main_title='NGS data result'
subtitle1=''
subtitle2=''
ylabel='Frequency'
ylabel2='Log(Frquency ratio)'
ylabel_ratio='Log(Frquency'
xlabel='Mutation type'
sorting_var='Type'   
type_name='Cycle' 
type_name2="Base"
type_name3="Spacer"
value_name_insrt='frequency' 
value_name_insrt2="freq_ratio"
value_name_insrt_ratio="freq_ratio"
font_label= {'size': 12,
             'weight' : 'bold',
             'color'  : 'Gray',
             'verticalalignment': 'baseline',
             'horizontalalignment': 'center'}
ylim_min=10**-5
ylim_max=10**0
bluered4=["#357FC2","#D23064","#76CAF2","#E37298"]
bluered2_1=["#357FC2","#D23064"] ##deep blue, red
bluered2_2=["#07abe3","#e35b72"] ##skyblue_2, pale red_2
redblue2_1=["#D23064","#357FC2"]
redblue2_2=["#e35b72","#07abe3"]
red12=["#D23064","#e35b72"]
black=["#000000"]
outlier=[]



def WGS_quartile_freq_lineplot(mutlist,i_merging_name,merging_filelist,o_filename,minmax_range=[], analysis_range=[50,50], outlierlist=[],typename="Mutation", date=Date, intype="csv", percentile_ref:float=50,
                    log=True, marker=True, error_style="band", error=lambda x: (x.min(), x.max()), linelist=[], graphrange=[], gap=25, alpha_=1,legend=True, ylim_min=2*10**-5, ylim_max=6*10**-1,
                    width=5.4, height=2.1, color=bluered4,marker_size=5,include_ref=False, rasterize=False,outtype="png"):
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
            if intype=="xlsx":
                df=pd.read_excel(filepath,"Sheet1")
            if intype=="csv":
                df=pd.read_csv(filepath)
            ## Select proper range ####
            df = df[df["position"].between(min_range, max_range) & ~df["position"].isin(outlierlist)]
            
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
            Mut_merged_df["Cycle"] = i_merging_name[count]
            
            
        #### Merging2: Merging by Experiment type#####
            merged_df=pd.concat([Mut_merged_df,merged_df])

    def window_mean_count_irregular(df, W, val_col="value", pos_col="position"):
        pos = df[pos_col].to_numpy(np.int64)
        val = df[val_col].to_numpy(np.float64)

        # vectorized boundaries
        left  = np.searchsorted(pos, pos - W, side="left")  # in each position[i], find the index of position[i]-W -->array of resulting index
        right = np.searchsorted(pos, pos + W, side="right") # in each position[i], find the index of position[i]+W -->array of resulting index

        # prefix sums for O(1) range sums
        csum = np.concatenate(([0.0], np.cumsum(val)))  # cumulative sume from start to denoted index
        counts = right - left
        sums = csum[right] - csum[left]
        means = sums / np.maximum(counts, 1)

        out = pd.DataFrame({
            "nnumber": counts.astype(np.int32),
            "mean": means
        }, index=df.index)
        return out

    W=analysis_range[0]

    parts = []
    for _, sub in merged_df.groupby(["Cycle", typename], sort=False, observed=True):   ### category name in _, each divided dataframe in sub
        stats = window_mean_count_irregular(sub, W, val_col=value_name_insrt, pos_col="position")  ##in each divided dataframe which have same cylce and typename, calculate mean value in each position -> array of means
        parts.append(pd.concat([sub, stats], axis=1)) 

    # merged_df2 = pd.concat(parts, axis=0).sort_index()
    merged_df2 = pd.concat(parts, axis=0)
    # merged_df2=merged_df2.sort_values(by=["Cycle", typename])
    percentile_y = "mean" if percentile_ref == 0 else f"{percentile_ref} percentile"
    if percentile_ref == 0:
        merged_df2.rename(columns={"mean": percentile_y}, inplace=True)
    else:
        # exact percentile : slice each window and compute quantile
        merged_df2=merged_df2.sort_values(by=["Cycle", typename])
        q = percentile_ref / 100.0
        def rolling_percentile_exact(df):
            pos = df["position"].to_numpy(np.int64)
            val = df[value_name_insrt].to_numpy(np.float64)
            left  = np.searchsorted(pos, pos - W, side="left")
            right = np.searchsorted(pos, pos + W, side="right")
            out = np.empty(len(df), dtype=np.float64)
            for i in range(len(df)):
                out[i] = np.quantile(val[left[i]:right[i]], q, method="linear")
            return out
        
        groupby= merged_df2.groupby(["Cycle", typename], sort=False, observed=True) ##as sort is false, resulting groupby would have exact same key with original merged_df2
        percent_df=groupby.apply(lambda s: pd.Series(rolling_percentile_exact(s)))

        merged_df2[percentile_y] = percent_df.reset_index(level=[0,1], drop=True).to_numpy()   ##exclude index of percent_df, resulting in series of percentile through index


    
    df_ref=merged_df2[merged_df2["Cycle"]==i_merging_name[0]]
    average_ref=df_ref[percentile_y].mean()
    merged_df_noref=merged_df2[merged_df2["Cycle"]!=i_merging_name[0]]   
    
    ###save to csv ###
    merged_df2.to_csv(f"{date}_{o_filename}-{percentile_ref}_percentile_line_{min_range}-{max_range}_{muttype}.csv") 


    ### making graph ###

    # 1) Fast Matplotlib settings

    mpl.style.use('fast')
    mpl.rcParams['path.simplify'] = True
    mpl.rcParams['path.simplify_threshold'] = 1.0

    File_name=f"{date}_{o_filename}-{percentile_ref}_percentile_line_{min_range}-{max_range}_{muttype}.{outtype}"
    fig, ax = plt.subplots(figsize=(width,height))
    colorlist=color[0:len(mutlist)]
    sns.set_palette(sns.color_palette(colorlist))
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    
    
    if include_ref:
        ref_line=sns.lineplot(
                data=df_ref,
                x="position",
                y=percentile_y,
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
                y=percentile_y,
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
    if rasterize:
        for x in ax.collections + ax.lines:
            x.set_rasterized(True)
    ax.set_xticks(range(a,b,gap))    
    ax.set_xticklabels(range(a,b,gap))
    ax.tick_params(axis='both', which='major')
    ax.set_xlabel("position",fontdict=font_label, labelpad=10)
    ax.set_ylabel(f"{percentile_y}",fontdict=font_label, labelpad=10)           
    ax.set_xlim(list(map(int,graphrange)))
    if log:
        ax.set(yscale="log")
        ax.set_ylabel("LOG("+ylabel+")",fontdict=font_label, labelpad=10)
    ax.set_ylim(ylim_min,ylim_max)
    fig.tight_layout()
    plt.savefig(File_name,dpi=600)
    print("graph done")    


def WGS_mut_finding(mutlist,i_merging_name,merging_filelist,o_filename,minmax_range=[], analysis_range=[50,50], outlierlist=[],typename="Mutation", date=Date, intype="csv",
                    q2limit: float=3*10**-4, q1limit: float=10**-4):
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
            if intype=="xlsx":
                df=pd.read_excel(filepath,"Sheet1")
            if intype=="csv":
                df=pd.read_csv(filepath)
            ## Select proper range ####
            df = df[df["position"].between(min_range, max_range) & ~df["position"].isin(outlierlist)]
            
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
            Mut_merged_df["Cycle"] = i_merging_name[count]
            
            
        #### Merging2: Merging by Experiment type#####
            merged_df=pd.concat([Mut_merged_df,merged_df])
    def window_mean_count_irregular(df, W, val_col="value", pos_col="position"):
        pos = df[pos_col].to_numpy(np.int64)
        val = df[val_col].to_numpy(np.float64)

        left  = np.searchsorted(pos, pos - W, side="left")     # indices for [pos-W, pos+W]
        right = np.searchsorted(pos, pos + W, side="right")

        csum = np.concatenate(([0.0], np.cumsum(val)))         # prefix sums
        counts = right - left
        sums = csum[right] - csum[left]
        means = sums / np.maximum(counts, 1)

        return pd.DataFrame({"nnumber": counts.astype(np.int32),
                            "mean":    means}, index=df.index)

    # ---- 1) Build base with mean/count per (Cycle, typename) group ----
    W = analysis_range[0]  # window half-width, assumes symmetric ±W
    parts = []
    for _, sub in merged_df.sort_values(["Cycle", typename, "position"], kind="mergesort") \
                        .groupby(["Cycle", typename], sort=False, observed=True):
        stats = window_mean_count_irregular(sub, W, val_col=value_name_insrt, pos_col="position")
        parts.append(pd.concat([sub, stats], axis=1))

    merged_df2 = pd.concat(parts, axis=0)  # all stats preserved

    # ---- 2) Exact rolling percentiles (p25/p50/p75) per group ----
    def rolling_percentiles_exact(df, W, val_col="value", pos_col="position"):
        pos = df[pos_col].to_numpy(np.int64)
        val = df[val_col].to_numpy(np.float64)
        left  = np.searchsorted(pos, pos - W, side="left")
        right = np.searchsorted(pos, pos + W, side="right")
        n = len(df)
        p25 = np.empty(n, dtype=np.float64)
        p50 = np.empty(n, dtype=np.float64)
        p75 = np.empty(n, dtype=np.float64)
        for i in range(n):
            w = val[left[i]:right[i]]
            if w.size:
                p25[i] = np.quantile(w, 0.25, method="linear")
                p50[i] = np.quantile(w, 0.50, method="linear")
                p75[i] = np.quantile(w, 0.75, method="linear")
            else:
                p25[i] = p50[i] = p75[i] = np.nan
        # IMPORTANT: return with the SAME INDEX so apply() aligns rows correctly
        return pd.DataFrame({"25 percentile": p25,
                            "50 percentile": p50,
                            "75 percentile": p75}, index=df.index)

    grouby = merged_df2.groupby(["Cycle", typename], sort=False, observed=True)
    per_df = grouby.apply(lambda s: rolling_percentiles_exact(s, W, value_name_insrt, "position"))
    # apply stacks results with a MultiIndex on rows (group keys + original index)
    per_df = per_df.reset_index(level=[0,1], drop=True)   # flatten → align 1:1 with merged_df2 rows

    # Attach percentiles back
    merged_df2 = pd.concat([merged_df2, per_df], axis=1)

    # ---- 3) Filter & sort by percentile standards ----
    p50_min = q2limit  # <-- set your threshold
    p25_min = q1limit  # <-- set your threshold

    filtered = merged_df2[
        (merged_df2["50 percentile"] >= p50_min) &
        (merged_df2["25 percentile"] >= p25_min)
    ]

    # Stable block order + deterministic inside-block order
    filtered = filtered.sort_values(
        ["50 percentile", "25 percentile", "Cycle", typename, "position"],
        ascending=[False, False, True, True, True],
        kind="mergesort"  # stable; keeps within-group relative order
    )

    # ---- 4) Save CSVs  ----
    # filtered/sorted subset
    filtered.to_csv(f"{date}_{o_filename}_Mutated_region_{muttype}.csv", index=False, float_format="%.6g")
    print(f"file created. total {len(filtered)} number of range detected")     


def WGS_depth_finding_csv(mutlist,i_merging_name,merging_filelist,o_filename, minmax_range=["321","562"], depth_standard:float=3*10**4, outlierlist=[],typename="Mutation", date=Date,ratio=False,intype="csv"):
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
            columnlist=["position","q30_depth"]
            temp_df=df[columnlist] ### Extract given columns from dataframe ###
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
            df_ref=merged_df[(merged_df["Cycle"]==i_merging_name[0])]
            df_noref=merged_df[(merged_df["Cycle"]==i_merging_name[1])]
            df_noref=df_noref.rename(columns={"q30_depth":"temp"+"q30_depth"}) 
            df_ref=df_ref.rename(columns={"q30_depth":"ref"+"q30_depth"})
            temp_df=pd.DataFrame()        
            temp_df=pd.merge(df_ref,df_noref,on="position",how="outer")
            temp_df=temp_df.dropna(axis=0)
            temp_df=temp_df[temp_df["ref"+"q30_depth"]>0]    ####frequency가 0으로 나누는 일이 없도록 limit 이하의 값을 가진건 버림####
            temp_df=temp_df[temp_df["temp"+"q30_depth"]>0]   ####frequency의 일관성을 위해 limit 이하의 값을 가진건 버림####
            temp_df["q30_depth_ratio"]=temp_df["temp"+"q30_depth"]/temp_df["ref"+"q30_depth"] ##divide by reference frequency        
            temp_df1=temp_df[["position","q30_depth_ratio","Cycle_x"]]
            ratio_df=pd.concat([temp_df1,ratio_df])
        ratio_df=ratio_df[(ratio_df["q30_depth_ratio"]>depth_standard)]
        ratio_df=ratio_df[["position","q30_depth_ratio","Cycle_x"]]
        ratio_df.to_csv(f"{date}_{o_filename}_depthratio.csv") 
        print(f"file created. total {len(ratio_df)} number of position detected")   
    if not ratio:
        merged_df=merged_df[(merged_df["q30_depth"]>depth_standard)]
        merged_df=merged_df[["position","q30_depth","Cycle"]]
        merged_df.to_csv(f"{date}_{o_filename}_depth.csv") 
        print(f"file created. total {len(merged_df)} number of position detected")  
  
def NGS_depth_lineplot_csv(mutlist,i_merging_name,merging_filelist,o_filename,
                 minmax_range=["321","562"], log=True, marker=True, error_style="band", error=lambda x: (x.min(), x.max()), linelist=[], graphrange=[], gap=25, alpha_=1,legend=True, ylim_min=2*10**-5, ylim_max=6*10**-1,
                 outlierlist=outlier,typename="Mutation", date=Date, width=5.4, height=2.1, color=bluered4,marker_size=5,include_ref=False,ratio=False,intype="xlsx", rasterize=True,outtype="pdf"):
    """
    Lineplot generating function
        merging_filelist : should be given as list of file paths list. e.g. [[a,b,c],[d,e,f]]
                    -> a,b,c and d,e,f would be merged by group respectively.
        i_merging_name : should be given as list of names. e.g. ["1","2"]
                    -> a,b,c file and d,e,f file would be named as 1 and 2 respectivley
    """
    min_range=int(minmax_range[0]) ###setting range of analysis
    max_range=int(minmax_range[1])

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
            columnlist=["position","q30_depth"]
            temp_df=df[columnlist] ### Extract given columns from dataframe ###
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

        df_ref=merged_df[(merged_df["Cycle"]==i_merging_name[0])]
        df_noref=merged_df[(merged_df["Cycle"]==i_merging_name[1])]
        df_noref=df_noref.rename(columns={"q30_depth":"temp"+"q30_depth"}) 
        df_ref=df_ref.rename(columns={"q30_depth":"ref"+"q30_depth"})
        temp_df=pd.DataFrame()        
        temp_df=pd.merge(df_ref,df_noref,on="position",how="outer")
        temp_df=temp_df.dropna(axis=0)
        temp_df=temp_df[temp_df["ref"+"q30_depth"]>0]    
        temp_df=temp_df[temp_df["temp"+"q30_depth"]>0]   
        temp_df["q30_depth ratio"]=temp_df["temp"+"q30_depth"]/temp_df["ref"+"q30_depth"] ##divide by reference frequency        
        temp_df1=temp_df[["position","q30_depth ratio","Cycle_x"]]
        ratio_df=pd.concat([temp_df1,ratio_df])
    
    if not ratio:
        df_ref=merged_df[merged_df["Cycle"]==i_merging_name[0]]
        average_ref=df_ref["q30_depth"].mean()
        merged_df_noref=merged_df[merged_df["Cycle"]!=i_merging_name[0]]
    
    
    ### making graph ###
    File_name=f"{date}_{o_filename}-merged_line_{min_range}-{max_range}.{outtype}"
    fig, ax = plt.subplots(figsize=(width,height))
    colorlist=color[0:len(mutlist)]

    mpl.style.use('fast')
    mpl.rcParams['path.simplify'] = True
    mpl.rcParams['path.simplify_threshold'] = 1.0

    sns.set_palette(sns.color_palette(colorlist))
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    if ratio:
        line=sns.lineplot(
                data=ratio_df,
                x="position",
                y="q30_depth ratio",
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
                    y="q30_depth",
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
                    y="q30_depth",
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
    if rasterize:
        for x in ax.collections + ax.lines:
            x.set_rasterized(True)
    ax.set_xticks(range(a,b,gap))    
    ax.set_xticklabels(range(a,b,gap))
    ax.tick_params(axis='both', which='major')
    ax.set_xlabel("position",fontdict=font_label, labelpad=10)
    ax.set_ylabel("depth",fontdict=font_label, labelpad=10)           
    ax.set_xlim(list(map(int,graphrange)))
    if log:
        ax.set(yscale="log")
        ax.set_ylabel("LOG("+ylabel+")",fontdict=font_label, labelpad=10)
    ax.set_ylim(ylim_min,ylim_max)
    fig.tight_layout()
    plt.savefig(File_name,dpi=600)


    print("graph done")    


"""Supplementary fig 10 LacZ w/o riboJ in B1-B5 region ecoli WGS"""
basicpath="RESPECTevo-Code_Rawdata/2_Plotting_Statistics/1_High-throughput_sequencing/8_Supple_fig10_wholegenome_rawdata/"

###lacZ gene###
Wholelist_Ecoli=[[["1"],["2"]],  
                 [["1"],["1"]]
                ]
namingdict={0:"B1B5_Cyc4",1:"B1B5_Cyc0"}
minmax_range_dict={0:["362455","365529"],1:["362455","365529"]} #lacZ
Ecolirange=[362300,365600] #lacZ
graphrange_dict={0:Ecolirange,1:Ecolirange}
linelist_dict={0:[],1:[]}
namelist=["Cycle0","Cycle4"]
filelist_temp=[]
for count0 in range(0,len(Wholelist_Ecoli)):
    filelist_temp=Wholelist_Ecoli[count0]
    for count in range(0,len(filelist_temp)):
        for count2 in range(0,len(filelist_temp[count])):
            filelist_temp[count][count2]=basicpath+f"_Exp19-{filelist_temp[count][count2]}_0_mutation.csv"
            
    WGS_quartile_freq_lineplot(["GtoA","TtoC"],namelist,filelist_temp,f"Set40_{namingdict[count0]}_linear_lacZonly",
                log=False,minmax_range=minmax_range_dict[count0] ,color=bluered2_1,graphrange=graphrange_dict[count0],linelist=linelist_dict[count0], percentile_ref=50,
                typename="Mutation",marker=True,alpha_=0.9,error_style="bars",legend=False,gap=500,outlierlist=[],marker_size=2.5,intype="csv",outtype="pdf",width=5.1,height=2.1
                ,ylim_max=2.5*10**-3,ylim_min=-10**-4)
    plt.close()

    WGS_quartile_freq_lineplot(["CtoT","AtoG"],namelist,filelist_temp,f"Set40_{namingdict[count0]}_linear_lacZonly",
                log=False,minmax_range=minmax_range_dict[count0] ,color=bluered2_2,graphrange=graphrange_dict[count0],linelist=linelist_dict[count0], percentile_ref=50,
                typename="Mutation",marker=True,alpha_=0.9,error_style="bars",legend=False,gap=500,outlierlist=[],marker_size=2.5,intype="csv",outtype="pdf",width=5.1,height=2.1
                ,ylim_max=2.5*10**-3,ylim_min=-10**-4)
    plt.close()    

###candidate1 (ybdA)
Wholelist_Ecoli=[[["1"],["2"]],  
                 [["1"],["1"]]
                ]
namingdict={0:"B1B5_Cyc4",1:"B1B5_Cyc0"}
minmax_range_dict={0:["1467000","1470000"],1:["1467000","1470000"]} #candidate1
Ecolirange=[1467000,1470000] #candidate1
graphrange_dict={0:Ecolirange,1:Ecolirange}
linelist_dict={0:[],1:[]}
namelist=["Cycle0","Cycle4"]
filelist_temp=[]
for count0 in range(0,len(Wholelist_Ecoli)):
    filelist_temp=Wholelist_Ecoli[count0]
    for count in range(0,len(filelist_temp)):
        for count2 in range(0,len(filelist_temp[count])):
            filelist_temp[count][count2]=basicpath+f"_Exp19-{filelist_temp[count][count2]}_0_mutation.csv"
            
    WGS_quartile_freq_lineplot(["GtoA","TtoC"],namelist,filelist_temp,f"Set40_{namingdict[count0]}_linear_candidate_only",
                log=False,minmax_range=minmax_range_dict[count0] ,color=bluered2_1,graphrange=graphrange_dict[count0],linelist=linelist_dict[count0],
                typename="Mutation",marker=True,alpha_=0.9,error_style="bars",legend=False,gap=500,outlierlist=[],marker_size=2.5,intype="csv",outtype="pdf",width=5.1,height=2.1
                ,ylim_max=2.5*10**-3,ylim_min=-10**-4)
    plt.close()

###candidate2 (hcp)
Wholelist_Ecoli=[[["1"],["2"]],  
                 [["1"],["1"]]
                ]
namingdict={0:"B1B5_Cyc4",1:"B1B5_Cyc0"}
minmax_range_dict={0:["912700","915700"],1:["912700","915700"]} #candidate2
Ecolirange=[912700,915700] #candidate2
graphrange_dict={0:Ecolirange,1:Ecolirange}
linelist_dict={0:[],1:[]}
namelist=["Cycle0","Cycle4"]
filelist_temp=[]
for count0 in range(0,len(Wholelist_Ecoli)):
    filelist_temp=Wholelist_Ecoli[count0]
    for count in range(0,len(filelist_temp)):
        for count2 in range(0,len(filelist_temp[count])):
            filelist_temp[count][count2]=basicpath+f"_Exp19-{filelist_temp[count][count2]}_0_mutation.csv"
            
    WGS_quartile_freq_lineplot(["CtoT","AtoG"],namelist,filelist_temp,f"Set40_{namingdict[count0]}_linear_candidate2_only",
                log=False,minmax_range=minmax_range_dict[count0] ,color=bluered2_1,graphrange=graphrange_dict[count0],linelist=linelist_dict[count0],
                typename="Mutation",marker=True,alpha_=0.9,error_style="bars",legend=False,gap=500,outlierlist=[],marker_size=2.5,intype="csv",outtype="pdf",width=5.1,height=2.1
                ,ylim_max=2.5*10**-3,ylim_min=-10**-4)
    plt.close()

#### whole genome ###
transversion_list=["AtoT","AtoC","TtoA","TtoG","GtoT","GtoC","CtoA","CtoG"]
black_list=["#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000"]
Wholelist_Ecoli=[[["1"],["2"]],  
                 [["1"],["1"]]
                ]
namingdict={0:"B1B5_Cyc4",1:"B1B5_Cyc0"}
minmax_range_dict={0:["1","4646332"],1:["1","4646332"]} ##whole genome 
Ecolirange=[0,4800000]
graphrange_dict={0:Ecolirange,1:Ecolirange}
linelist_dict={0:[],1:[]}
namelist=["Cycle0","Cycle4"]
filelist_temp=[]
for count0 in range(0,len(Wholelist_Ecoli)):
    filelist_temp=Wholelist_Ecoli[count0]
    for count in range(0,len(filelist_temp)):
        for count2 in range(0,len(filelist_temp[count])):
            filelist_temp[count][count2]=basicpath+f"_Exp19-{filelist_temp[count][count2]}_0_mutation.csv"

    WGS_quartile_freq_lineplot(["TtoC","GtoA"],namelist,filelist_temp,f"Set40_wholegenome_{namingdict[count0]}_linear_raterized",
            log=False,minmax_range=minmax_range_dict[count0] ,color=redblue2_1,graphrange=graphrange_dict[count0],linelist=linelist_dict[count0],
            typename="Mutation",marker=True,alpha_=0.9,error_style="bars",legend=False,gap=500000,outlierlist=[],marker_size=0.8,intype="csv",outtype="pdf",width=7,height=2.1, rasterize=True
            ,ylim_max=2.5*10**-3,ylim_min=-10**-4)    
    plt.close()

    WGS_quartile_freq_lineplot(["AtoG","CtoT"],namelist,filelist_temp,f"Set40_wholegenome_{namingdict[count0]}_linear_raterized",
            log=False,minmax_range=minmax_range_dict[count0] ,color=redblue2_2,graphrange=graphrange_dict[count0],linelist=linelist_dict[count0],
            typename="Mutation",marker=True,alpha_=0.9,error_style="bars",legend=False,gap=500000,outlierlist=[],marker_size=0.8,intype="csv",outtype="pdf",width=7,height=2.1, rasterize=True
            ,ylim_max=2.5*10**-3,ylim_min=-10**-4)    
    plt.close()


    WGS_quartile_freq_lineplot(transversion_list,namelist,filelist_temp,f"Set40_wholegenome_black_{namingdict[count0]}_linear_raterized",
            log=False,minmax_range=minmax_range_dict[count0] ,color=black_list,graphrange=graphrange_dict[count0],linelist=linelist_dict[count0],
            typename="Mutation",marker=True,alpha_=0.9,error_style="bars",legend=False,gap=500000,outlierlist=[],marker_size=0.8,intype="csv",outtype="pdf",width=7,height=2.1, rasterize=True
            ,ylim_max=2.5*10**-3,ylim_min=-10**-4)    
    plt.close()


#### Depth plotting ####
Wholelist_Ecoli=[[["1"],["2"]]]
namingdict={0:"B1B5_Cyc4"}
minmax_range_dict={0:["1","4646332"]} ##whole genome 
Ecolirange=[0,4800000]
graphrange_dict={0:Ecolirange}
namelist=["Cycle0","Cycle4"]
filelist_temp=[]
for count0 in range(0,len(Wholelist_Ecoli)):
    filelist_temp=Wholelist_Ecoli[count0]
    for count in range(0,len(filelist_temp)):
        for count2 in range(0,len(filelist_temp[count])):
            filelist_temp[count][count2]=basicpath+f"_Exp19-{filelist_temp[count][count2]}_0_mutation.csv"    

    NGS_depth_lineplot_csv(["CtoT"],namelist,filelist_temp,f"Set40_depth_{namingdict[count0]}",
            log=True,minmax_range=minmax_range_dict[count0],color=black,graphrange=graphrange_dict[count0],linelist=[], ylim_min=10**2, ylim_max=10**6,
            typename="Mutation",marker=True,alpha_=0.9,error_style="bars",legend=False,gap=500000,outlierlist=[],marker_size=0.8,intype="csv",outtype="pdf",width=7,height=2.1)
    plt.close()
    NGS_depth_lineplot_csv(["CtoT"],namelist,filelist_temp,f"Set40_depth_{namingdict[count0]}",
            log=True,minmax_range=minmax_range_dict[count0],color=black,graphrange=graphrange_dict[count0],linelist=[], ylim_min=10**2, ylim_max=10**6,
            typename="Mutation",marker=True,alpha_=0.9,error_style="bars",legend=False,gap=500000,outlierlist=[],marker_size=0.8,intype="csv",outtype="pdf",width=7,height=2.1)
    plt.close()
    NGS_depth_lineplot_csv(["CtoT"],namelist,filelist_temp,f"Set40_depthratio_{namingdict[count0]}",
            log=True,minmax_range=minmax_range_dict[count0],color=black,graphrange=graphrange_dict[count0],linelist=[], ratio=True, ylim_min=10**-2, ylim_max=10**2,
            typename="Mutation",marker=True,alpha_=0.9,error_style="bars",legend=False,gap=500000,outlierlist=[],marker_size=0.8,intype="csv",outtype="pdf",width=7,height=2.1)
    plt.close()
    
#### finding point exceeding criteria ###    
    WGS_depth_finding_csv(["CtoT"],namelist,filelist_temp,f"Set40_{namingdict[count0]}",minmax_range=minmax_range_dict[count0])
    WGS_depth_finding_csv(["CtoT"],namelist,filelist_temp,f"Set40_{namingdict[count0]}",minmax_range=minmax_range_dict[count0],ratio=True, depth_standard=10)

    WGS_mut_finding(All_mut_list,namelist,filelist_temp,f"Set40_criteria_{namingdict[count0]}_q2_4e-4",minmax_range=minmax_range_dict[count0],q1limit=0,q2limit=4*10**-4)

