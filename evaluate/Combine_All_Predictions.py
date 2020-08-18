
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