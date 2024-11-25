# MNanalyser 231122
# input: python MNanalyser.py {GNPS cytoscape result}
#       {samplelist} -db {accurate mass list}
import argparse
import pandas as pd
import os
import glob
from datetime import datetime
from function import calc_ratio, get_tops, get_compound_name, classify_mn, classify_single, get_isotope_intensity, get_first_decimal
import sys

version = 1.0

parser = argparse.ArgumentParser()

parser.add_argument("gnps_result_folder", type=str, help="Directory containing input GNPS cytoscape result files is required!!")
parser.add_argument("sample_list", type=str, help="Input Analyzed Sample list (tab separated format, eg. tsv file, txt file) is required!! You must create this file (tab separated format). ex. GNPS_sample_list.tsv")
parser.add_argument(
    "-mse", type=float, default="0.05", help="ms search error value (default;+- 0.05)")
    
args = parser.parse_args()


samplelist = args.sample_list
if os.path.isfile(samplelist) is False:
    print("ERROR!: The SECOND ARGUMENT must be the sample list")
    print("Make it by yourself as comma or tab separated text, as follow")
    print("Sample\tStrain\tMedium")
    print("or")
    print("Sample\tStrain\tMedium\tCondition")
    print("Without header")
else:
    df_sample = pd.read_table(samplelist, sep="\t")
    for i in [".wiff", ".abf", "wiffscan"]:
        df_sample.iloc[:, 0] = df_sample.apply(
            lambda x: x.str.replace(i, "", regex=False)).iloc[:, 0]

columns = ["Alignment ID", "Average Rt(min)", "Average Mz",
           "Metabolite name", "Adduct type",
           "MS/MS assigned", "S/N average",
           "MS1 isotopic spectrum", "MS/MS spectrum"]

gnpsdir = args.gnps_result_folder
if os.path.isdir(gnpsdir) is False:
    print(
        "ERROR!: The FIRST ARGUMENT must be /PATH/FOR/YOUR/MN/ANALYSIS/RESULT/FILE!!"
    )
    sys.exit
else:
    df_quant = pd.read_table(
        glob.glob(f"{gnpsdir}/quantification_table/*.txt")[0],
        sep="\t",
        header=3,
        low_memory=False,
        usecols=lambda x: x in (columns +
                                df_sample.iloc[:, 0].values.tolist()))
    """""""""
    Concatenate the result files obtained by GNPS analysis.
    Analyzed MS data, DB search result and MN information are mearged.
    """""""""
    df_dbresult = pd.read_table(
        glob.glob(f"{gnpsdir}/DB_result/*.tsv")[0], sep="\t",
        usecols=["#Scan#", "Compound_Name", "Smiles", "INCHI"])
    df_dbresult= df_dbresult.rename(
        {"#Scan#": "Alignment ID"}, axis=1)
    df_clusterinfo = pd.read_table(
        glob.glob(f"{gnpsdir}/clusterinfo_summary/*.tsv")[0], sep="\t",
        usecols=["cluster index", "componentindex"]).rename(
        {"cluster index": "Alignment ID"}, axis=1)

df_quant = pd.merge(df_quant, df_clusterinfo,
                    on="Alignment ID", how="left")
df_quant = pd.merge(df_quant, df_dbresult,
                    on="Alignment ID", how="left")


"""""""""
Comparison of compound production based on ion intensity across all samples.

In this step, information about the ion intensity values and the sample name 
of the sample that produced the most compounds and the sample that produced 
the second most compounds is obtained.

Top1: A certain compound with the strongest ion intensity among all samples
Top1_sample: Sample name of Top1
Top2: A certain compound with the second strongest ion intensity among all samples
Top2_sample: Sample name of Top2
"""""""""

grouped = df_sample.groupby(df_sample.columns[1])
groups, cols = [], []
for group, data in grouped:
    samples = data.iloc[:, 0].to_list()
    df_quant[f"{group}_max"] = df_quant[samples].max(axis=1)
    df_quant[f"{group}_max_sample"] = df_quant[samples].idxmax(axis=1)
    cols.append(f"{group}_max")
    groups.append(group)

df_quant[["First_Ratio", "Max_Ratio", "Ratios", "Max_Ratio_Num"]
            ] = df_quant.apply(calc_ratio, columns=cols, axis=1)
df_quant[["Top1", "Top1_sample", "Top2", "Top2_sample"]
            ] = df_quant.apply(get_tops, columns=cols, axis=1)


df_quant["Top1_group"] = df_quant["Top1_sample"].apply(
    lambda x: df_sample[df_sample.iloc[:, 0] == x].iloc[0, 1]
)

controlname = "control"
control_samples = df_sample[df_sample.iloc[:, 1]
                            == controlname].iloc[:, 0].to_list()
df_quant["Medium"] = False
df_quant.loc[
    df_quant["Metabolite name"].apply(
        lambda x: any(substring in x for substring in control_samples)
    )
    | df_quant["Top1_sample"].apply(
        lambda x: x in control_samples
    ), "Medium"
] = True

df_quant["Medium(MSMS)"] = False
df_quant.loc[(
    df_quant["Metabolite name"].apply(
        lambda x: any(substring in x for substring in control_samples)
    )
    & ~ df_quant["Metabolite name"].str.contains("w/o MS2"))
    | df_quant["Top1_sample"].apply(
        lambda x: x in control_samples
), "Medium(MSMS)"
] = True


df_quant["GNPS_Hit"] = df_quant["Compound_Name"].notna()


df_quant["Local_Hit"] = False
df_quant.loc[
    (~ df_quant["Medium"] & (df_quant["Metabolite name"] != "Unknown")), "Local_Hit"] = True
df_quant["Local_Hit(MSMS)"] = False
df_quant.loc[
    (~ df_quant["Medium"] &
     (df_quant["Metabolite name"] != "Unknown") &
     (df_quant["MS/MS spectrum"].notna()) &
     ~ df_quant["Metabolite name"].str.contains("w/o MS2")), "Local_Hit"
] = True


df_quant["Hit_Compound"] = df_quant.apply(get_compound_name, axis=1)


df_quant["Unknown_Metabolite"] = True
df_quant.loc[
    (df_quant["GNPS_Hit"]
     | df_quant["Local_Hit"]
     | df_quant["Medium"]
     ), "Unknown_Metabolite"
] = False

df_quant["Unknown_Metabolite(MSMS)"] = True
df_quant.loc[
    (df_quant["GNPS_Hit"]
     | df_quant["Local_Hit(MSMS)"]
     | df_quant["Medium(MSMS)"]
     ), "Unknown_Metabolite(MSMS)"
] = False



for i in ["", "(MSMS)"]:
    series = df_quant.groupby(
        "componentindex").apply(classify_mn, accuracy=i)
    series.name = f"MN_type{i}"
    df_quant = pd.merge(df_quant, series, on="componentindex", how="inner")
    df_quant.loc[df_quant["componentindex"] == -1, f"MN_type{i}"] = (
        df_quant.loc[df_quant["componentindex"] == -1, :].apply(
            classify_single, accuracy=i, axis=1
        )
    )



for target in ["Max_Ratio", "Max_Ratio_Num", "First_Ratio"]:
    for func in ["mean", "max", "min"]:
        df_quant[f"Net_{func}_({target})"] = df_quant.groupby(
            "componentindex")[target].transform(func)
        df_quant.loc[
            df_quant["componentindex"] == -1, f"Net_{func}_({target})"
        ] = (
            df_quant.loc[df_quant["componentindex"] == -1, target]
        )


df_quant["MN_Group_Variety"] = (
    df_quant.groupby("componentindex")["Top1_group"].transform("nunique")
)
df_quant.loc[
    df_quant["componentindex"] == -1, "MN_Group_Variety"
] = 1


now = datetime.now()

formatted_datetime = now.strftime("%y%m%d%H%M%S")
print("========================================================================\n")
print("Output analysis file is" + f"{formatted_datetime}_mndataMNA{version}.txt\n")
print("========================================================================")
df_quant.to_csv(f"{formatted_datetime}_mndataMNA{version}.txt",
                sep="\t", index=False)