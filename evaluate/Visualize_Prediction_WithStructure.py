

import os
import numpy as np
from ops.os_operation import mkdir
def Visualize_Prediction_WithStructure(save_path,map_name,Final_Predict_file,factor,real_loc_ref,output_str):
    Final_Predict_Dict={}
    with open(Final_Predict_file,'r') as file:
        line=file.readline()
        while line:
            line=line.strip()
            split_result=line.split()
            key=split_result[0]
            Final_Predict_Dict[key]=int(split_result[1])
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
            if key not in real_loc_ref:
                continue#because it's background, ignore predictions
            tmp_coordinate=real_loc_ref[key]
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
        file.write("color green, chain A and "+current_obj_name+"\n")
        file.write("color yellow, chain B and " + current_obj_name + "\n")
        file.write("color red, chain C and " + current_obj_name + "\n")
        file.write("color cyan, chain D and " + current_obj_name + "\n")
        file.write("select coil, chain A and "+current_obj_name+"\n")
        file.write("select beta, chain B and " + current_obj_name + "\n")
        file.write("select alpha, chain C and " + current_obj_name + "\n")
        file.write("select DNA_RNA, chain D and " + current_obj_name + "\n")
        file.write("bg_color 0\n")


def Visualize_Confident_Prediction_WithStructure(save_path,map_name,Final_Predict_file,factor,real_loc_ref,output_str):
    Final_Predict_Dict={}
    Final_Prob_Dict = {}
    n_class=4
    with open(Final_Predict_file,'r') as file:
        line=file.readline()
        while line:
            line=line.strip()
            split_result=line.split()
            key=split_result[0]
            tmp_label=int(split_result[1])
            Final_Predict_Dict[key]=tmp_label
            Final_Prob_Dict[key]=float(split_result[2+tmp_label])#the prob of the predicted class by our Model
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
            tmp_prob=Final_Prob_Dict[key]
            if tmp_prob<0.9:
                continue
            tmp_chain = chain_dict[tmp_label]
            if key not in real_loc_ref:
                continue#because it's background, ignore predictions
            tmp_coordinate = real_loc_ref[key]
            line = "ATOM%7d  %3s %3s%2s%4d    " % (natm, "CA ", "ALA", " " + tmp_chain, natm)
            line += "%8.3f%8.3f%8.3f%6.2f\n" % (tmp_coordinate[0], tmp_coordinate[1], tmp_coordinate[2], 1)
            predfile.write(line)
    #color it based on our definition
    #using pymol -u *.pml to open it
    tmp_visual_script_path = os.path.join(save_path, map_name + output_str + "_predC.pml")
    with open(tmp_visual_script_path,'w') as file:
        file.write("load "+map_name +output_str+ "_predC.pdb\n")
        current_obj_name=map_name +output_str+ "_predC"
        file.write("show spheres, "+current_obj_name+"\n")
        file.write("set sphere_scale, 0.5\n")
        file.write("color green, chain A and "+current_obj_name+"\n")
        file.write("color yellow, chain B and " + current_obj_name + "\n")
        file.write("color red, chain C and " + current_obj_name + "\n")
        file.write("color cyan, chain D and " + current_obj_name + "\n")
        file.write("select coil, chain A and " + current_obj_name + "\n")
        file.write("select beta, chain B and " + current_obj_name + "\n")
        file.write("select alpha, chain C and " + current_obj_name + "\n")
        file.write("select DNA_RNA, chain D and " + current_obj_name + "\n")
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