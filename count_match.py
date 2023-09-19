import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
import glob
import re

def residue_finder(res):

    res_list = res.strip("[]").split(", ")
    res_dic = {}
    der_dic = {}
    for one_res in res_list:
        resn = one_res.split(" ")[0][1:]
        at = resn.split("-")[-1][:2]
        perc = float(re.search(r"\d+.\d+", one_res.split(" ")[1]).group())
        if (resn[:3] in ["GLU", "ASP", "ARG"]) and (at in ["OE", "OD", "NH"]):
            new_resn = resn[:-1]
            if new_resn not in der_dic.keys():
                der_dic[new_resn] = perc
            elif new_resn in der_dic.keys():
                der_dic[new_resn]+=perc
        else:
            if perc >= 50.0:
                res_dic[resn] = perc

    if len(der_dic) != 0:
        for der in der_dic.keys():
            if der_dic[der] >= 50.0:
                res_dic[der] = der_dic[der]

    return res_dic


def match_calc(df):
    prot_lig = df['Prot - Lig']
    HS_idx = df['HS_idx']
    Distance = df['Distance \n(Lig - HS)']
    Lig_feat = df['Ligand \nFeature']
    Wat_feat = df['Water \nFeature']
    Don_Prot = df['Protein atoms \n(Prot to HS %)']
    Acc_Prot = df['Protein atoms \n(HS to Prot %)']
    total_num_int = len(prot_lig)

    gmatch=0
    for hs in HS_idx:
        if math.isnan(hs) == False:
            gmatch+=1

    fcnl_grp_match_perc = round((gmatch/total_num_int)*100, 2)

    fmatch=0
    feat_match_idx = []
    for i in range(total_num_int):
        lig_feat = Lig_feat[i]
        wat_feat = Wat_feat[i]
        wat_prot_don = Don_Prot[i]
        wat_prot_acc = Acc_Prot[i]

        if (type(wat_prot_acc) == str) and ("NONE" not in wat_prot_acc):
            wat_prot_acc_dic = residue_finder(wat_prot_acc)
        if (type(wat_prot_don) == str) and ("NONE" not in wat_prot_don):
            wat_prot_don_dic = residue_finder(wat_prot_don)

        if type(lig_feat) == str and type(wat_feat) == str:
            if (lig_feat == wat_feat) or (lig_feat in wat_feat) or (wat_feat in lig_feat):
                lig_prot = prot_lig[i].split(" - ")[-1]
                lig_prot_at = prot_lig[i].split(" - ")[-1].split("-")[-1][:2]
                if (lig_prot[:3] in ["GLU", "ASP", "ARG"]) and (lig_prot_at in ["OE", "OD", "NH"]):
                    lig_prot = lig_prot[:-3]
                if lig_feat == "Acceptor":
                    if lig_prot in wat_prot_don_dic.keys():
                        fmatch += 1
                        feat_match_idx.append(i)
                elif lig_feat == "Donor":
                    if lig_prot in wat_prot_acc_dic.keys():
                        fmatch += 1
                        feat_match_idx.append(i)
                elif lig_feat == "Acceptor/Donor":
                    if (lig_prot in wat_prot_don_dic.keys()) or (lig_prot in wat_prot_acc_dic.keys()):
                        fmatch += 1
                        feat_match_idx.append(i)

    feat_fcnl_grp_match_perc = round((fmatch / total_num_int) * 100, 2)
    # print(feat_match_idx)
    lmatch=0
    for i in feat_match_idx:
        dist = Distance[i]
        # print(i, dist)
        lig_feat = Lig_feat[i]
        wat_feat = Wat_feat[i]
        if math.isnan(dist) == False:
            if dist <= 1.4:
                if type(lig_feat) == str and type(wat_feat) == str:
                    if (lig_feat == wat_feat) or (lig_feat in wat_feat) or (wat_feat in lig_feat):
                        lmatch+=1

    same_loc_feat_match_perc = round((lmatch/total_num_int)*100, 2)

    lmatch2=0
    for i in range(total_num_int):
        dist = Distance[i]
        # print(i, dist)
        lig_feat = Lig_feat[i]
        wat_feat = Wat_feat[i]
        if math.isnan(dist) == False:
            if dist <= 1.4:
                lmatch2 += 1
    same_loc_match_perc = round((lmatch2/total_num_int)*100, 2)


    return total_num_int, fcnl_grp_match_perc, feat_fcnl_grp_match_perc, same_loc_feat_match_perc
    # return total_num_int, same_loc_match_perc

file_path = "/yeonji-data/yeonji/Dropbox/myfolder_data/wbp_results/"
sim = "ahr"
data_table_files = sorted(glob.glob(file_path+sim+"_match2/fcnl_grp_close/*csv"))
# data_table_files = sorted(glob.glob(file_path+"MATCH/"+sim+"/fcnl_grp_close/*csv"))
out_df_inverse = pd.DataFrame(columns=range(len(data_table_files)), index=['System', 'Total_num_int', 'Functional_Grp_Match', 'Feature_Functional_Grp_Match', 'Position_Feature_Match'])
# out_df_inverse2 = pd.DataFrame(columns=range(len(data_table_files)), index=['System', 'Total_num_int', 'Position_Match'])

for n, file in enumerate(data_table_files):
    sys_name = file.split("/")[-1].split("_")[0]
    print(sys_name)
    data_table = pd.read_csv(file)
    match_perc = match_calc(data_table)
    if match_perc[0] >= 3:
        # out_df_inverse2[n] = [sys_name, match_perc[0], match_perc[1]]
        out_df_inverse[n] = [sys_name, match_perc[0], match_perc[1], match_perc[2], match_perc[3]]

    # print(match_perc[0], match_perc[1], match_perc[2], match_perc[3])
out_df = out_df_inverse.T
out_df.dropna(subset=["System"], inplace=True)
out_df = out_df.reset_index(drop=True)

out_df.to_csv(file_path+sim+"_match2/"+sim+"_position_match.csv")
# out_df.to_csv(file_path+"MATCH/"+sim+"/match.csv")


range_perc = np.arange(0.0, 110.0, 10)
fcn_grp_match = (out_df["Functional_Grp_Match"].value_counts(bins=range_perc, sort=False)).tolist()
feat_match = (out_df["Feature_Functional_Grp_Match"].value_counts(bins=10, sort=False)).tolist()
out_df["Feature"] = pd.cut(out_df["Feature_Functional_Grp_Match"], bins=range_perc, right=True)
counts = (out_df["Feature"].value_counts().sort_index()).tolist()

loc_match = (out_df["Position_Feature_Match"].value_counts(bins=range_perc, sort=False)).tolist()


max_fcn_score = max(fcn_grp_match)
max_feat_score = max(feat_match)
max_loc_score = max(loc_match)


fig, axes = plt.subplots(1, 3, figsize=(12, 4))

axes[0].bar([x - 4.6 for x in range_perc[1:11]], fcn_grp_match, align="center", width=7, color='#18c983')
axes[0].set_xticks(range_perc[0:11])
axes[0].set_yticks(np.arange(0, max_fcn_score + 1, 1))
axes[0].set_xlabel("%")
axes[0].set_ylabel(r"Frequency")
axes[0].text(0, 11, "  Atom matches \n                 / \nProt-Lig Interaction", style='italic', fontsize=6.4, bbox={'facecolor': '#e1f5f7', 'alpha': 0.5})
axes[0].grid(color='gray', alpha=0.5)
# axes[0].set_title("Atom match")

axes[1].bar([x - 5 for x in range_perc[1:11]], feat_match, align="center", width=7, color='#0c88e8')
axes[1].set_xticks(range_perc[0:11])
axes[1].set_yticks(np.arange(0, max_feat_score + 1, 1))
axes[1].set_xlabel("%")
axes[1].set_ylabel(r"Frequency")
axes[1].text(0, 6.0, "Feature matches \n                / \n Prot-Lig Fetures", style='italic', fontsize=6.4, bbox={'facecolor': '#e1f5f7', 'alpha': 0.5})
axes[1].grid(color='gray', alpha=0.5)
axes[1].set_title("Feature match")


axes[2].bar([x - 5 for x in range_perc[1:11]], loc_match, align="center", width=7, color='#fa7537')
axes[2].set_xticks(range_perc[0:11])
axes[2].set_yticks(np.arange(0, max_loc_score + 1, 1))
axes[2].set_xlabel("%")
axes[2].set_ylabel(r"Frequency")
axes[1].text(0, 6.0, "Feature matches \n                / \n Prot-Lig Fetures", style='italic', fontsize=6.4, bbox={'facecolor': '#e1f5f7', 'alpha': 0.5})
axes[2].grid(color='gray', alpha=0.5)

plt.tight_layout()
fig.savefig(file_path+sim+"_match/"+sim+"_match_analysis.png")
# plt.show()