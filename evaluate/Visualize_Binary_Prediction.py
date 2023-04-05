# Publication:  "Emap2sec+:  Detecting Protein and DNA/RNA Structures in Cryo-EM Maps of Intermediate Resolution Using Deep Learning", Xiao Wang, Eman Alnabati, Tunde W. Aderinwale, Sai Raghavendra Maddhuri Venkata Subramaniya, Genki Terashi, and Daisuke Kihara, BioRxiv (2020)

# Emap2sec+ is a computational tool using deep learning that can accurately identify structures, alpha helices, beta sheets, other(coils/turns) and DNA/RNA, in cryo-Electron Microscopy (EM) maps of medium to low resolution.
# Copyright (C) 2020 Xiao Wang, Eman Alnabati, Tunde W Aderinwale, Sai Raghavendra Maddhuri, Genki Terashi, Daisuke Kihara, and Purdue University.
# License: GPL v3 for academic use. (For commercial use, please contact us for different licensing.)
# Contact: Daisuke Kihara (dkihara@purdue.edu)


#

# This program is free software: you can redistribute it and/or modify

# it under the terms of the GNU General Public License as published by

# the Free Software Foundation, version 3.

#

# This program is distributed in the hope that it will be useful,

# but WITHOUT ANY WARRANTY; without even the implied warranty of

# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the

# GNU General Public License V3 for more details.

#

# You should have received a copy of the GNU v3.0 General Public License

# along with this program.  If not, see https://www.gnu.org/licenses/gpl-3.0.en.html.

import os
import numpy as np
from ops.os_operation import mkdir

def Visualize_Binary_Prediction(save_path,map_name,Final_Predict_file,factor,output_str):
    Final_Predict_Dict={}
    with open(Final_Predict_file,'r') as file:
        line=file.readline()
        while line:
            line=line.strip()
            split_result=line.split()
            key=split_result[0]
            tmp_label = int(split_result[1])
            tmp_label = 0 if tmp_label<=2 else 1
            Final_Predict_Dict[key]=tmp_label
            line=file.readline()
    save_path=os.path.join(save_path,output_str)
    mkdir(save_path)
    tmp_visual_pred_path = os.path.join(save_path, map_name +output_str+ "_pred.pdb")
    natm = 1
    chain_dict = Build_Chain_ID()
    Nstep = 11
    back = int((float(Nstep) - 1) / 2)
    with open(tmp_visual_pred_path, 'w') as predfile:
        for key in Final_Predict_Dict.keys():
            tmp_label = Final_Predict_Dict[key]
            tmp_chain = chain_dict[tmp_label]
            coordinate=key.split(",")
            tmp_coordinate=[]
            for tmp_loc_idx in range(3):
                tmp_loc=coordinate[tmp_loc_idx]
                tmp_coordinate.append(float(tmp_loc)*factor+back)
            line = "ATOM%7d  %3s %3s%2s%4d    " % (natm, "CA ", "ALA", " " + tmp_chain, natm)
            line += "%8.3f%8.3f%8.3f%6.2f\n" % (tmp_coordinate[0], tmp_coordinate[1], tmp_coordinate[2], 1)
            predfile.write(line)
    #color it based on our definition
    #using pymol -u *.pml to open it
    tmp_visual_script_path = os.path.join(save_path, map_name + output_str + "_pred.pml")
    with open(tmp_visual_script_path,'w') as file:
        file.write("load "+map_name +output_str+ "_pred.pdb\n")
        current_obj_name=map_name +output_str+ "_pred"
        file.write("show spheres, "+current_obj_name+"\n")
        file.write("set sphere_scale, 0.5\n")
        file.write("color red, chain A and "+current_obj_name+"\n")
        file.write("color cyan, chain B and " + current_obj_name + "\n")
        file.write("select protein, chain A and "+current_obj_name+"\n")
        file.write("select DNA_RNA, chain B and " + current_obj_name + "\n")
        file.write("bg_color 0\n")

def Visualize_Binary_Confident_Prediction(save_path,map_name,Final_Predict_file,factor,output_str):
    Final_Predict_Dict={}
    Final_Prob_Dict = {}
    with open(Final_Predict_file,'r') as file:
        line=file.readline()
        while line:
            line=line.strip()
            split_result=line.split()
            key=split_result[0]
            tmp_label = int(split_result[1])
            tmp_label = 0 if tmp_label<=2 else 1
            if tmp_label==0:
                tmp_prob=0
                for k in range(3):
                    tmp_prob += float(split_result[2 + tmp_label])
            else:
                tmp_prob=float(split_result[5])

            Final_Prob_Dict[key]=tmp_prob
            Final_Predict_Dict[key]=tmp_label
            line=file.readline()
    save_path=os.path.join(save_path,output_str)
    mkdir(save_path)
    tmp_visual_pred_path = os.path.join(save_path, map_name +output_str+ "_predC.pdb")
    natm = 1
    chain_dict = Build_Chain_ID()
    Nstep = 11
    back = int((float(Nstep) - 1) / 2)
    with open(tmp_visual_pred_path, 'w') as predfile:
        for key in Final_Predict_Dict.keys():
            tmp_label = Final_Predict_Dict[key]
            tmp_prob = Final_Prob_Dict[key]
            if tmp_prob < 0.8:
                continue
            tmp_chain = chain_dict[tmp_label]
            coordinate=key.split(",")
            tmp_coordinate=[]
            for tmp_loc_idx in range(3):
                tmp_loc=coordinate[tmp_loc_idx]
                tmp_coordinate.append(float(tmp_loc)*factor+back)
            line = "ATOM%7d  %3s %3s%2s%4d    " % (natm, "CA ", "ALA", " " + tmp_chain, natm)
            line += "%8.3f%8.3f%8.3f%6.2f\n" % (tmp_coordinate[0], tmp_coordinate[1], tmp_coordinate[2], 1)
            predfile.write(line)
    #color it based on our definition
    #using pymol -u *.pml to open it
    tmp_visual_script_path = os.path.join(save_path, map_name + output_str + "_pred.pml")
    with open(tmp_visual_script_path,'w') as file:
        file.write("load "+map_name +output_str+ "_pred.pdb\n")
        current_obj_name=map_name +output_str+ "_pred"
        file.write("show spheres, "+current_obj_name+"\n")
        file.write("set sphere_scale, 0.5\n")
        file.write("color red, chain A and "+current_obj_name+"\n")
        file.write("color cyan, chain B and " + current_obj_name + "\n")
        file.write("select protein, chain A and "+current_obj_name+"\n")
        file.write("select DNA_RNA, chain B and " + current_obj_name + "\n")
        file.write("bg_color 0\n")


def Build_Chain_ID():
    chain_dict={}
    chain_dict[0]='A'
    chain_dict[1]='B'
    chain_dict[2] = 'C'
    chain_dict[3] = 'D'
    chain_dict[4] = 'E'
    chain_dict[-1] = 'E'
    return chain_dict