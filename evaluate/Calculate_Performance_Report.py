#Publication:  "Emap2sec+:  Detecting Protein and DNA/RNA Structures in Cryo-EM Maps of Intermediate Resolution Using Deep Learning", Xiao Wang, Eman Alnabati, Tunde W. Aderinwale, Sai Raghavendra Maddhuri Venkata Subramaniya, Genki Terashi, and Daisuke Kihara, BioRxiv (2020)

# Emap2sec+ is a computational tool using deep learning that can accurately identify structures, alpha helices, beta sheets, other(coils/turns) and DNA/RNA, in cryo-Electron Microscopy (EM) maps of medium to low resolution.
# Copyright (C) 2020 Xiao Wang, Eman Alnabati, Tunde W Aderinwale, Sai Raghavendra Maddhuri, Genki Terashi, Daisuke Kihara, and Purdue University.
# License: GPL v3 for academic use. (For commercial use, please contact us for different licensing.)
# Contact: Daisuke Kihara (dkihara@purdue.edu)


#

# This program is free software: you can redistribute it and/or modify

# it under the terms of the GNU General Public License as published by

# the Free Software Foundation, either version 3 of the License, or

# (at your option) any later version.

#

# This program is distributed in the hope that it will be useful,

# but WITHOUT ANY WARRANTY; without even the implied warranty of

# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the

# GNU General Public License V3 for more details.

#

# You should have received a copy of the GNU v3.0 General Public License

# along with this program.  If not, see https://www.gnu.org/licenses/gpl-3.0.en.html.

import os
from sklearn.metrics import classification_report

def Calculate_Performance_Report(save_path, map_name,input_path,phase_pred_file,save_str):
    #first get the prediciton dict
    input_size=11
    Final_Predict_Dict = {}
    with open(phase_pred_file, 'r') as file:
        line = file.readline()
        while line:
            line = line.strip()
            split_result = line.split()
            key = split_result[0]
            Final_Predict_Dict[key] = int(split_result[1])
            line = file.readline()
    Label_Dict={}
    with open(input_path,'r') as file:
        line = file.readline()
        while line:
            line = line.strip('\n')
            result = line.split(',')
            if len(result) < 100 or len(result) != input_size ** 3 + 4:  # no map info
                line = file.readline()
                continue
            assert len(result) == input_size ** 3 + 4
            tmp_position = ''
            for k in range(3):
                tmp_position += str(result[k]) + ','
            tmp_label=str(result[3])
            Label_Dict[tmp_position]=tmp_label
            line=file.readline()
    Include_Type_Set=set()
    Predict_Array=[]
    True_Array=[]
    for key in Final_Predict_Dict:
        if key in Label_Dict:
            tmp_label=Label_Dict[key]
            if tmp_label=="-1":#background
                continue
            cur_pred = Final_Predict_Dict[key]
            Predict_Array.append(cur_pred)
            possible_list = tmp_label.split(":")
            first_correct = int(possible_list[0])
            possible_set = set()
            for item in possible_list:
                possible_set.add(int(item))
            if int(cur_pred) not in possible_set:
                True_Array.append(first_correct)
                Include_Type_Set.add(first_correct)
            else:
                True_Array.append(cur_pred)
            Include_Type_Set.add(cur_pred)


    #write report
    final_path=os.path.join(save_path,map_name+save_str+"_report.txt")
    file=open(final_path,'w')
    target_names=['Coil','Beta','Alpha','DNA/RNA']
    final_target_names=[x for k,x in enumerate(target_names) if k in Include_Type_Set]#sometimes no beta/alpha for input map
    print(classification_report(True_Array,Predict_Array,target_names=final_target_names),file=file)
    file.close()
