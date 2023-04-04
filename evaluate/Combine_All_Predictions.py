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
from collections import defaultdict
import numpy as np
from scipy.special import softmax
def Combine_All_Predictions(save_path, map_name, All_Output_File):
    Final_Predict_Dict = defaultdict(list)
    num_classes=4

    for Predict_file in All_Output_File:
        with open(Predict_file, 'r') as file:
            line = file.readline()
            while line:
                line = line.strip()
                split_result = line.split()
                key = split_result[0]
                pred_list=[]
                for k in range(num_classes):
                    pred_list.append(float(split_result[k+2]))
                Final_Predict_Dict[key].append(pred_list)
                line = file.readline()
    final_pred_path = os.path.join(save_path, map_name + "_final_pred.txt")
    with open(final_pred_path,'w') as file:
        for key in Final_Predict_Dict:
            pred_list=Final_Predict_Dict[key]
            pred_record=np.zeros(4)
            for tmp_list in pred_list:
                for k in range(num_classes):
                    pred_record[k]+=tmp_list[k]
            pred_record=softmax(pred_record)
            pred_label=int(np.argmax(pred_record))
            file.write(key+" "+str(pred_label)+" ")
            for prob in pred_record:
                file.write(str(prob)+" ")
            file.write("\n")
            file.flush()
    return final_pred_path
