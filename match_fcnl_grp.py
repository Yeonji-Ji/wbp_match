import pandas as pd
import math
import scipy
import mdtraj as md
import numpy as np
import re
import glob
import pymol

########################################################################################################################
# Generate prot-lig interaction and prot-water interaction matches based on the functional groups.
# change the path of the files from line 525.
########################################################################################################################
# class and functions are defined below #

class Find_WBP:
    def __init__(self, topo, traj, ligand, cent_crds, hs_summary):

        self.topo = topo
        self.traj = traj
        self.ligand = ligand  # hydration site water
        self.cent_crds = cent_crds
        self.hs_summary = hs_summary
        # self.mpo = mpo

        self.lig_don_h = {}      # {H-index: don_N_index}
        self.lig_acc = []        # {index}
        self.lig_hv_at = {}      # {index: element}
        self.prot_don_h = {}     # {H-index: don_N_index}
        self.prot_acc = []       # {index}
        self.prot_hv_at = {}     # {index: element}

        self.PL_hbond = []       # Pairs of [P, L] within the criteria of distance and angle.
        self.PL_no_angle = []    # Pairs of [P, L] within the distance criteria.
        self.prot_res_list = []

        self.res_match_dic = {}

    # Function to return pairs within distance of cutoff.
    def find_pair(self, pair, criteria=1.2):

        pair_list = np.array([[pair[0], pair[1]]])
        dist = md.compute_distances(self.traj, pair_list)
        roundup = round(dist[0][0]*10, 1)

        # if dist * 10 <= criteria:
        if roundup <= criteria:
            return pair


    def ligand_fcn(self):

        lig_hydrogen = {}

        # From complex file, extract ligand N, O, H index.
        for at in self.ligand:
            element = self.topo.atom(at).element.symbol
            if element in ["N", "O", "S"]:
                self.lig_hv_at[at] = element
            if element == "H":
                lig_hydrogen[at] = element

        # Possible complex of ligand N, O and H.
        lig_at_h = [[x, y] for x in lig_hydrogen.keys() for y in self.lig_hv_at.keys()]

        # Use find_pair function to find bonded N,O,S-H pairs. Fill the dictionary, lig_don_h.
        for pair in lig_at_h:
            bond_pair = self.find_pair(pair)
            if bond_pair:
                self.lig_don_h[pair[0]] = pair[1]

        # Among heavy atoms in lig_hv_at dictionary, if the atom does not have a bonded hydrogen, add to the list, lig_acc.
        for at in self.lig_hv_at:
            if at not in self.lig_don_h.values():
                self.lig_acc.append(at)


        return self.lig_don_h, self.lig_acc


    def protein_fcn(self):

        prot_hydrogen = {}

        solutes = md.compute_neighbors(self.traj, 0.4, query_indices=self.lig_hv_at)  # List of proteins within 4 A of ligand heavy atoms.
        solutes = [idx for idx in solutes[0] if idx not in self.ligand]  # Exclude indices of ligand.

        # From complex file, extract protein N, O, H index.
        for at in solutes:
            element = self.topo.atom(at).element.symbol

            if element in ["N", "O", "S"]:
                self.prot_hv_at[at] = element
            if element == "H":
                prot_hydrogen[at] = element

        # Possible complex of ligand N, O and H.
        prot_at_h = [[x, y] for x in prot_hydrogen.keys() for y in self.prot_hv_at.keys()]

        # Use find_pair function to find bonded N-H or O-H pairs. Fill the dictionary, lig_don_h.
        for pair in prot_at_h:
            bond_pair = self.find_pair(pair, 1.2)
            if bond_pair:
                self.prot_don_h[pair[0]] = pair[1]

        # Among heavy atoms in lig_hv_at dictionary, if the atom does not have a bonded hydrogen, add to the list, lig_acc.
        for at in self.prot_hv_at:
            if at not in self.prot_don_h.values():
                self.prot_acc.append(at)

        # for don in self.prot_don_h.keys():
        #     if self.prot_don_h[don] == 3560:
        #         print(self.prot_don_h[don], self.topo.atom(3560).name, self.topo.atom(3560).residue, don)

        return self.prot_don_h, self.prot_acc


    def prot_lig_interaction(self, distance=3.6, angle=30):

        don_acc_element = ["oxygen", "sulfur"]
        prot_resn_list = []
        PL_within_d = []        # Including pairs of {distance : [protein heavy atom, ligand heavy atom]} within cutoff. Always [P, L].

        # Possible complex of protein heavy atoms and ligand heavy atoms. Always [P, L]
        prot_lig_hv_at = [[x, y] for x in self.prot_hv_at for y in self.lig_hv_at]
        # print(prot_lig_hv_at)

        # for at in [4821, 4822, 4823, 4825, 4827, 4829, 4831]:
        #     print(self.topo.atom(at))

        # pair_list = np.array([[3560, 4831]])
        # dist = md.compute_distances(self.traj, pair_list)
        # print(round(dist[0][0]*10, 1))

        # Use find_pair function to find protein heavy atom and ligand heavy atom pairs within input distance (default=3.6 A). Add to the list, PL_within_d.
        for pair in prot_lig_hv_at:
            hb_dist_pair = self.find_pair(pair, distance)
            if hb_dist_pair:
                PL_within_d.append(hb_dist_pair)

        # From possible P-distance-L pairs, if either of a pair is in *_don_h dictionary, find the angle A-D-H.
        # If donor atom is found, the other is always counted as an acceptor atom.
        for pair in PL_within_d:
            if pair[0] in self.prot_don_h.values():  # Check protein (donor) atom.
                if pair[1] in self.lig_acc or (str(self.topo.atom(pair[1]).element) in don_acc_element):
                    hydrogens = [k for k, v in self.prot_don_h.items() if v == pair[0]]
                    if len(hydrogens) == 1:
                        h1 = hydrogens[0]
                        triplet_pd = np.array([[pair[1], pair[0], h1]])
                        tri_angle = md.compute_angles(self.traj, triplet_pd)
                        deg = round(math.degrees(tri_angle), 3)
                        self.PL_no_angle.append([deg, pair[0], pair[1], "pd"])
                        if deg <= angle:
                            self.PL_hbond.append([deg, pair[0], pair[1], h1, "pd"])
                            prot_res = int(self.topo.atom(pair[0]).residue.index)
                            prot_resn_list.append(str(self.topo.atom(pair[0]).residue.name) + str(prot_res))
                            self.prot_res_list.append(prot_res+1)

                    elif len(hydrogens) > 1:
                        for n in range(len(hydrogens)):
                            h1 = hydrogens[n]
                            triplet_pd = np.array([[pair[1], pair[0], h1]])
                            tri_angle = md.compute_angles(self.traj, triplet_pd)
                            deg = round(math.degrees(tri_angle), 3)
                            self.PL_no_angle.append([deg, pair[0], pair[1], "pd"])
                            if deg <= angle:
                                self.PL_hbond.append([deg, pair[0], pair[1], h1, "pd"])
                                prot_res = int(self.topo.atom(pair[0]).residue.index)
                                prot_resn_list.append(str(self.topo.atom(pair[0]).residue.name) + str(prot_res))
                                self.prot_res_list.append(prot_res+1)


            if pair[1] in self.lig_don_h.values():  # Check Ligand (donor) atom.
                if pair[0] in self.prot_acc or (str(self.topo.atom(pair[0]).element) in don_acc_element):
                    hydrogens = [k for k, v in self.lig_don_h.items() if v == pair[1]]
                    if len(hydrogens) == 1:
                        h1 = hydrogens[0]
                        triplet_ld = np.array([[pair[0], pair[1], h1]])
                        tri_angle = md.compute_angles(self.traj, triplet_ld)
                        deg = round(math.degrees(tri_angle), 3)
                        self.PL_no_angle.append([deg, pair[0], pair[1], "ld"])
                        if deg <= angle:
                            self.PL_hbond.append([deg, pair[0], pair[1], h1, "ld"])
                            prot_res = int(self.topo.atom(pair[0]).residue.index)
                            prot_resn_list.append(str(self.topo.atom(pair[0]).residue.name) + str(prot_res))
                            self.prot_res_list.append(prot_res+1)

                    elif len(hydrogens) > 1:
                        for n in range(len(hydrogens)):
                            h1 = hydrogens[n]
                            triplet_ld = np.array([[pair[0], pair[1], h1]])
                            deg = md.compute_angles(self.traj, triplet_ld)
                            deg = round(math.degrees(deg), 3)
                            self.PL_no_angle.append([deg, pair[0], pair[1], "ld"])
                            if deg <= angle:
                                self.PL_hbond.append([deg, pair[0], pair[1], h1, "ld"])
                                prot_res = int(self.topo.atom(pair[0]).residue.index)
                                prot_resn_list.append(str(self.topo.atom(pair[0]).residue.name) + str(prot_res))
                                self.prot_res_list.append(prot_res+1)

        self.prot_res_list = list(set(self.prot_res_list))

        return self.PL_hbond, self.PL_no_angle, self.prot_res_list


    def assign_wat_feature(self, acc_rate, don_rate, hbRate=0.6, greaseRate=0.2):

        if acc_rate >= hbRate and don_rate < hbRate:
            wat_feature = "Acceptor"
        elif acc_rate < hbRate and don_rate >= hbRate:
            wat_feature = "Donor"
        elif acc_rate >= hbRate and don_rate >= hbRate:
            wat_feature = "Acceptor/Donor"
        elif acc_rate < greaseRate and don_rate < greaseRate:
            wat_feature = "Grease"
        elif (greaseRate <= acc_rate < hbRate and don_rate < greaseRate) or (greaseRate <= don_rate < hbRate and acc_rate < greaseRate) \
                or (greaseRate <= don_rate < hbRate and greaseRate <= acc_rate < hbRate):
            wat_feature = "N/A"

        return wat_feature


    def wat_feature_uniorm(self, feat):
        if feat == "acceptor":
            new_feat = "Acceptor"
        elif feat == "donor":
            new_feat = "Donor"
        elif feat == "both":
            new_feat = "Acceptor/Donor"
        elif feat == "grease":
            new_feat = "Grease"

        return new_feat


    def col_rename_and_extract(self):
        fixed_col = ['HS_idx', 'Acc_sw', 'Don_sw', 'wat_to_prot_res', 'prot_to_wat_res']
        column_list = list(self.hs_summary.columns)
        column_renamed = {}
        for col in column_list:
            if col in ["index", "Hydration_Site"]:
                column_renamed[col] = "HS_idx"
            if col in ["Acc_sw", "Don_sw"]:
                column_renamed[col] = col
            if col in ["solute_acceptors", "protein_acceptor_percentage"]:
                column_renamed[col] = "wat_to_prot_res"
            if col in ["solute_donors", "protein_donor_percentage"]:
                column_renamed[col] = "prot_to_wat_res"
            if col in ["nwat", "occupancy", "Feature"]:
                column_renamed[col] = col

        to_add = list(col for col in column_renamed.values() if col not in fixed_col)
        df_out = self.hs_summary.rename(columns=column_renamed)[fixed_col + to_add]

        return df_out


    def format_res_for_match(self, res):
        if (type(res) == str) and (res != "NONE"):
            if ("[" or "]") in res:
                res_list = re.findall("\[.*?\]", res)
                replaced_res = [one.replace("[", "").replace("]", "") for one in res_list]
                new_format = []
                for r in replaced_res:
                    res = r.split(", ")
                    res_name_idx = res[0].split("-")[0]
                    res_name = re.search(r'[a-zA-Z]+', res_name_idx).group()
                    res_idx = int(re.search(r'\d+', res_name_idx).group())+1
                    at_name = res[0].split("-")[1]
                    perc = res[1]
                    new_name = res_name + str(res_idx) + "-" + at_name + " (" + perc + ")"
                    new_format.append(new_name)
            else:
                res_list = res.split(",")
                new_format = []
                for r in res_list:
                    res = r.split('-')
                    res_name = re.search(r'[a-zA-Z]+', res[0]).group()
                    res_idx = int(re.search(r'\d+', res[0]).group())+1
                    at_name = res[1]
                    perc = str("{:.2f}".format(float(res[2])*100))
                    new_name = res_name  + str(res_idx) + "-" + at_name + " (" + perc + ")"
                    new_format.append(new_name)
        else:
            new_format = ["NONE"]
        # print(new_format)
        return new_format


    def glu_asp_process_for_hs(self, res):
        res_dic = {}
        # print(res)
        for one in res:
            if one != "NONE":
                res_at = one.split(" ")[0]      # GLH92-OE1
                resn = res_at.split("-")[0][:3] # GLH
                atn = res_at.split("-")[1][:2]  # OE1
                if (resn in ["GLU"] and "OE" in atn) or (resn in ["ASP"] and "OD" in atn) or (resn in ["ARG"] and "NH" in atn):
                    res_at = res_at[:-1]
                elif resn == "GLH":
                    res_at = "GLU"+res_at[3:-1]
                elif resn == "ASH":
                    res_at == "ASP"+res_at[3:-1]
                perc = float(one.split(" ")[1].strip("()"))
                if res_at not in res_dic.keys():
                    res_dic[res_at] = perc
                else:
                    res_dic[res_at] += perc

            new_res = [res + " (" + str(round(res_dic[res], 3)) + ")" for res in res_dic.keys()]
        return new_res


    def find_match(self):

        feature_list = []
        prot_at_list = [str(self.topo.atom(li[1])) for li in self.PL_hbond]
        self.res_match_dic = {}
        for key in prot_at_list:
            atn = key.split("-")[1][:2]
            if (key[:3] in ["ASP", "GLU", "ARG"]) and (atn in ["OE", "OD", "NH"]):
                key = key[:-1]
            self.res_match_dic[key] = []

        hs_summary_new = self.col_rename_and_extract()

        hs_summary_new['wat_to_prot_res'] = hs_summary_new['wat_to_prot_res'].apply(self.format_res_for_match)
        hs_summary_new['prot_to_wat_res'] = hs_summary_new['prot_to_wat_res'].apply(self.format_res_for_match)
        res_list = hs_summary_new['wat_to_prot_res'].apply(self.glu_asp_process_for_hs) + hs_summary_new['prot_to_wat_res'].apply(self.glu_asp_process_for_hs)
        num_hs = len(hs_summary_new)
        # print(hs_summary_new)
        # hs_summary_new.to_csv("/home/yeonji/akt1_bbr_summary_new.csv")

        for idx in range(num_hs):
            res_dict = {}
            for res in res_list[idx]:
                if res != 'NONE':
                    resi = res.split(' ')[0]
                    perc = float(re.findall(r'\((\d+.\d+)\)', res)[0])
                    if resi not in res_dict.keys():
                        res_dict[resi] = perc
                    elif resi in res_dict.keys():
                        res_dict[resi] += perc
            # print(res_dict)
            # print(self.res_match_dic)
            # print(hs_summary_new.loc[idx, 'wat_to_prot_res'])
            for res in res_dict.keys():
                if res in self.res_match_dic.keys():
                    rate = res_dict[res]
                    if rate >= 50.0:
                        # print(res)
                        if idx not in self.res_match_dic[res]:
                            self.res_match_dic[res].append(idx)

            acc_rate = hs_summary_new.loc[idx, 'Acc_sw']
            don_rate = hs_summary_new.loc[idx, 'Don_sw']
            feature = self.assign_wat_feature(acc_rate=acc_rate, don_rate=don_rate)
            feature_list.append(feature)

        if 'Feature' in hs_summary_new.columns:
            hs_summary_new['Feature'] = hs_summary_new["Feature"].apply(self.wat_feature_uniorm)
        elif 'Feature' not in hs_summary_new.columns:
            hs_summary_new['Feature'] = feature_list



        # print(self.res_match_dic)
        return hs_summary_new, self.res_match_dic


    def summarize_data(self, hs_summary_new, distance=1.4, hbRate=0.6, greaseRate=0.2):
        global column_list
        lig_at_list = []
        two_int = {}
        for n, pair in enumerate(self.PL_hbond):
            if pair[2] not in lig_at_list:
                lig_at_list.append(pair[2])
            elif pair[2] in lig_at_list:
                two_int[pair[2]] = n
        # print(self.PL_hbond)
        # print(lig_at_list)
        PL_list = []
        # for at in lig_at_list:
        #     print(self.topo.atom(at))
        for lig_at in lig_at_list:

            prot_at_list = [pair[1] for pair in self.PL_hbond if pair[2] == lig_at]
            prot_res_list = [str(self.topo.atom(prot_at)) for prot_at in prot_at_list]
            # print(prot_res_list)
            # prot_res_set = set(prot_res_list)
            # print(prot_res_set)
            prot_res = ", ".join(prot_res_list)

            # print(prot_res)

            lig_res = str(self.topo.atom(lig_at).residue.name)+"-"+str(self.topo.atom(lig_at)).split("-")[-1]
            # print(lig_res)

            angle_PL = []
            for pair in self.PL_hbond:
                if pair[2] == lig_at:
                    angle = pair[0]
                    angle_feat = pair[4]
                    angle_at = self.topo.atom(pair[1])
                    resi = self.topo.atom(pair[1]).residue.index+1
                    if angle_feat == "pd":
                        angle_str = str(angle) + " (" + str(angle_at) + " to LIG)"
                    elif angle_feat == "ld":
                        angle_str = str(angle) + " (LIG to " + str(angle_at) + ")"
                    angle_PL.append(angle_str)
            # print(angle_PL, prot_res)
            PL_list.append([angle_PL, lig_res, prot_res, lig_at])

        summary_all_HS_dic = {}
        a=0
        summary_closest_HS_dic = {}
        c=0
        # print(PL_list)
        for PL in sorted(PL_list):
            lig_at = PL[-1]
            prot_at = PL[2].split(", ")
            res_DER_shortened = []
            res_DER_renamed = []

            for p_at in prot_at:
                atn = p_at.split("-")[1][:2]
                if (p_at[:3] in ["GLU", "ASP", "ARG"]) and (atn in ("OE", "OD", "NH")):
                    p_at = p_at[:-1]
                    p_at_for_data = p_at+"1,2"
                else:
                    p_at = p_at
                    p_at_for_data = p_at
                if p_at not in res_DER_shortened:
                    res_DER_shortened.append(p_at)
                if p_at not in res_DER_renamed:
                    res_DER_renamed.append(p_at_for_data)

            res_DER_shortened = set(res_DER_shortened)
            res_DER_renamed = list(set(res_DER_renamed))

            for i, residue in enumerate(res_DER_shortened):
                # print(res_DER_renamed[i])

                if ("nwat" and "occupancy") in hs_summary_new.columns:
                    column_list = ['Prot - Lig', 'HS_idx', 'Distance \n(Lig - HS)', 'n_WAT', 'Occupancy',
                                         'Ligand \nFeature', 'Water \nFeature', 'HS Acc rate \n(Prot to HS)',
                                         'HS Don rate \n(HS to Prot)', 'Protein atoms \n(Prot to HS %)',
                                         'Protein atoms \n(HS to Prot %)', 'Prot_Lig_angle']
                elif ("nwat" and "occupancy") not in hs_summary_new.columns:
                    column_list = ['Prot - Lig', 'HS_idx', 'Distance \n(Lig - HS)', 'Ligand \nFeature',
                                         'Water \nFeature', 'HS Acc rate \n(Prot to HS)',
                                         'HS Don rate \n(HS to Prot)', 'Protein atoms \n(Prot to HS %)',
                                         'Protein atoms \n(HS to Prot %)', 'Prot_Lig_angle']
                column_length = len(column_list)

                hs_list = self.res_match_dic[residue]

                prot_lig_int = "LIG" + PL[1][3:] + " - " + res_DER_renamed[i]
                # print(PL[1][3:])
                # print(prot_lig_int)

                prot_lig_angle_list = [PL[0][j] for j, angle in enumerate(PL[0]) if residue in angle]
                if len(prot_lig_angle_list) == 1:
                    prot_lig_angle = prot_lig_angle_list[0]
                elif len(prot_lig_angle_list) >= 2:
                    prot_lig_angle = ", ".join(prot_lig_angle_list)
                # print(prot_lig_angle)
                if ("LIG to" in prot_lig_angle) and ("to LIG" not in prot_lig_angle):
                    lig_feature = "Donor"
                elif ("LIG to" not in prot_lig_angle) and ("to LIG" in prot_lig_angle):
                    lig_feature = "Acceptor"
                elif ("LIG to" in prot_lig_angle) and ("to LIG" in prot_lig_angle):
                    lig_feature = "Acceptor/Donor"



                if len(hs_list) >= 1:
                    lig_hb_cords = self.traj.xyz[0, lig_at, :] * 10
                    dist_list = scipy.spatial.distance.cdist([lig_hb_cords], self.cent_crds[hs_list], "euclidean")[0]
                    closest_HS = hs_list[np.argmin(dist_list)]

                    for m, hsite in enumerate(hs_list):

                        # print(prot_lig_int)
                        # print(prot_lig_angle)

                        HS_feature = hs_summary_new.loc[hsite, 'Feature']
                        distance = round(scipy.spatial.distance.cdist([lig_hb_cords], [self.cent_crds[hsite]], "euclidean")[0][0], 3)
                        acc_rate = hs_summary_new.loc[hsite, 'Acc_sw']
                        Prot2Wat = hs_summary_new.loc[hsite, 'prot_to_wat_res']
                        don_rate = hs_summary_new.loc[hsite, 'Don_sw']
                        Wat2Prot = hs_summary_new.loc[hsite, 'wat_to_prot_res']

                        if ("nwat" and "occupancy") in hs_summary_new.columns:
                            n_wat = hs_summary_new.loc[hsite, 'nwat']
                            occupancy = hs_summary_new.loc[hsite, 'occupancy']
                            match_HS_list = [prot_lig_int, str(hsite), str(distance), n_wat, occupancy, lig_feature, HS_feature, acc_rate, don_rate, Prot2Wat, Wat2Prot, prot_lig_angle]

                        if ("nwat" and "occupancy") not in hs_summary_new.columns:
                            match_HS_list = [prot_lig_int, str(hsite), str(distance), lig_feature, HS_feature, acc_rate, don_rate, Prot2Wat, Wat2Prot, prot_lig_angle]
                            #print(match_HS_list)

                        summary_all_HS_dic[a] = match_HS_list
                        a+=1
                        if hsite == closest_HS:
                            summary_closest_HS_dic[c] = match_HS_list
                            c+=1

                elif len(hs_list) == 0:
                    if column_length == 12:
                        one_int_list = [prot_lig_int, None, None, None, None, lig_feature, None, None, None, None, None, prot_lig_angle]
                    if column_length == 10:
                        one_int_list = [prot_lig_int, None, None, lig_feature, None, None, None, None, None, prot_lig_angle]
                        # print(one_int_list)
                    summary_all_HS_dic[a] = one_int_list
                    a+=1
                    summary_closest_HS_dic[c] = one_int_list
                    c+=1
        #         if match_HS_list:
        #             print(match_HS_list)
        #         elif one_int_list:
        #             print(one_int_list)

        # print(summary_all_HS_dic)
        summary_closest_data = pd.DataFrame.from_dict(summary_closest_HS_dic, orient='index', columns=column_list)
        # print(summary_closest_data.T)
        summary_closest_data = summary_closest_data.reset_index()
        # print(summary_closest_data)
        summary_all_data = pd.DataFrame.from_dict(summary_all_HS_dic, orient='index', columns=column_list)
        summary_all_data = summary_all_data.reset_index()


        ### PYMOL ###

        PL_Dist_list_PyMOL = []                     # THIS LIST IS FOR VISUALIZATION
        for n, pair in enumerate(self.PL_hbond):
            if pair[-1] == "pd":
                dist_pair = [pair[2], pair[3]]      #(acc-Lig, donH)
            elif pair[-1] == "ld":
                dist_pair = [pair[1], pair[3]]      #(acc-Prot, donH)
            PL_Dist_list_PyMOL.append(dist_pair)

        ### hs_summary ###
        PyMOL_HS_dic = {'acc': [], 'don': [], 'both':[], 'grease': [], 'na': []}
        for i in range(len(hs_summary_new)):
            HS_idx = str(hs_summary_new.loc[i, 'HS_idx'])
            if hs_summary_new.loc[i, 'Feature'] == "Donor":
                PyMOL_HS_dic['don'].append(HS_idx)
            if hs_summary_new.loc[i, 'Feature'] == "Acceptor":
                PyMOL_HS_dic['acc'].append(HS_idx)
            if hs_summary_new.loc[i, 'Feature'] == "Acceptor/Donor":
                PyMOL_HS_dic['both'].append(HS_idx)
            if hs_summary_new.loc[i, 'Feature'] == "Grease":
                PyMOL_HS_dic['grease'].append(HS_idx)
            if hs_summary_new.loc[i, 'Feature'] == "N/A":
                PyMOL_HS_dic['na'].append(HS_idx)


        return summary_closest_data, summary_all_data, PL_Dist_list_PyMOL, PyMOL_HS_dic



########################################################################################################################
### Run the class and functions here ###

######################################### P A T H ############################################
complex_path = "/Users/yeonji/Dropbox (Lehman College)/myfolder_data/wbp_results/ahr_complex/"
file_path = "/Users/yeonji/Dropbox (Lehman College)/myfolder_data/wbp_results/"
##############################################################################################

complex_files = sorted(glob.glob(complex_path + "sahh*pdb"))


for file in complex_files:
    target = file.split("/")[-1].split("_")[0]
    print(target)

    if target:

        ### Input files path ###
        sim = "ahr"
        hs_summary_file = pd.read_csv(file_path+sim+"_summary/csv/"+target+"_min_"+sim+"_hsa_summary.csv")
        clustercenter_file = md.load(file_path+sim+"_cc/"+target+"_"+sim+"_cc.pdb")
        mpo_file = file_path+sim+"_mpo/"+target+"_"+sim+"_unshifted_mpo.pdb"

        ### Load the input files ###
        prot_traj = md.load(file)
        prot_top = prot_traj.topology
        prot = prot_top.select("is_protein")
        lig = prot_top.select("resname UNK")
        crds_cent = clustercenter_file.xyz[0, :, :] * 10


        ### run class and functions ###
        WBP = Find_WBP(prot_top, prot_traj, lig, crds_cent, hs_summary_file)
        WBP.ligand_fcn()
        WBP.protein_fcn()
        pl_int = WBP.prot_lig_interaction()                      # Default: distance=3.6A and angle=30
        df_new = WBP.find_match()[0]
        match_dic = WBP.find_match()[1]
        data = WBP.summarize_data(df_new)
        data_table_closest = data[0]
        data_table_all = data[1]
        print(data_table_closest)

        ### save the results ###
        data_table_closest.to_csv(file_path+sim+"_match/fcnl_grp_close/"+target+"_"+sim+"_match_data_close.csv")
        data_table_all.to_csv(file_path+sim+"_match/fcnl_grp_all/"+target+"_"+sim+"_match_data_all.csv")



        ### Generate PyMOL pse file ###
        PyMOL_pl_int_resi_list = [str(at) for at in pl_int[2]]
        PyMOL_pl_pair = data[2]
        PyMOL_hs_dict = data[3]

        match_hs = [hs for hs in set(data[0]['HS_idx']) if type(hs)==str]

        pl_int_resi_sele_str = '+'.join(PyMOL_pl_int_resi_list)
        hs_acc_sele_str = '+'.join(PyMOL_hs_dict['acc'])
        hs_don_sele_str = '+'.join(PyMOL_hs_dict['don'])
        hs_acc_don_sele_str = "+".join(PyMOL_hs_dict['both'])
        hs_grease_sele_str = "+".join(PyMOL_hs_dict['grease'])
        hs_na_sele_str = "+".join(PyMOL_hs_dict['na'])
        if len(match_hs) >=1:
            match_hs_str = "+".join(match_hs)

        pymol.cmd.load(file, "comp")
        pymol.cmd.load(mpo_file, "mpo")

        pymol.cmd.create("prot", "comp and not organic")
        pymol.cmd.create("lig", "comp and organic and resn UNK")

        pymol.cmd.color("grey70", "prot and element c")
        pymol.cmd.color("cyan", "lig and element c")

        pymol.cmd.select("side_chains", "prot and resi " + pl_int_resi_sele_str)
        pymol.cmd.show("sticks", "side_chains")
        pymol.cmd.hide("everything", "mpo")

        pymol.cmd.set("stick_radius", 0.2, "prot")
        pymol.cmd.set("stick_radius", 0.1, "lig")
        pymol.cmd.set("stick_h_scale", 0.8, "lig")



        if hs_acc_sele_str:
            pymol.cmd.select("acc", "mpo and resi " + hs_acc_sele_str)
            pymol.preset.ball_and_stick("acc", mode=1)
            pymol.cmd.color("red", "acc and element o")
        if hs_don_sele_str:
            pymol.cmd.select("don", "mpo and resi " + hs_don_sele_str)
            pymol.preset.ball_and_stick("don", mode=1)
            pymol.cmd.color("blue", "don and element o")
        if hs_acc_don_sele_str:
            pymol.cmd.select("acc_don", "mpo and resi " + hs_acc_don_sele_str)
            pymol.preset.ball_and_stick("acc_don", mode=1)
            pymol.cmd.color("purple", "acc_don and element o")
        if hs_grease_sele_str:
            pymol.cmd.select("grease", "mpo and resi " + hs_grease_sele_str)
            pymol.preset.ball_and_stick("grease", mode=1)
            pymol.cmd.color("green", "grease and element o")

        pymol.cmd.disable("mpo")
        if match_hs_str:
            pymol.cmd.create("match_HS", "mpo and resi " + match_hs_str)
            pymol.preset.ball_and_stick("match_HS", mode=1)


        for n, pair in enumerate(PyMOL_pl_pair):
            at1 = "comp and index " + str(pair[0] + 1)
            at2 = "comp and index " + str(pair[1] + 1)
            pymol.cmd.distance("pl" + str(n), at1, at2)

        pymol.cmd.hide("label")
        pymol.cmd.group("prot_lig_int", "pl*")
        pymol.cmd.disable("comp")

        pymol.cmd.set("bg_rgb", "[1, 1, 1]")
        pymol.cmd.set("depth_cue", "0")
        pymol.cmd.set("spec_reflect", 0)
        pymol.cmd.set("spec_power", 800)
        pymol.cmd.set("ambient", 0.5)
        pymol.cmd.set("dash_width", 2)
        pymol.cmd.set("dash_gap", 0.3)
        pymol.cmd.set("antialias", 2)
        pymol.cmd.set("light_count", 1)
        pymol.cmd.set("cartoon_transparency", 0.9)
        pymol.cmd.set("ray_trace_mode", 1)
        pymol.cmd.set("ray_trace_gain", 1)
        pymol.cmd.set("ray_texture", 1)

        ### Save pymol sessions in the path ###
        pymol.cmd.save(file_path+sim+"_pymol/"+target + "_"+sim+"_vis.pse")
        pymol.cmd.delete("*")
