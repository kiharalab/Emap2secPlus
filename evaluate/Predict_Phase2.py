import os
import torch
from ops.os_operation import mkdir
import numpy as np
from torch.autograd import Variable
from torch import nn
import torch.nn.functional as F

def Predict_Phase2(save_path,map_name,Phase1_Final_Prediction_Dict,indicate,fold,batch_size):
    input_size = 11
    dimension = input_size
    n_classes = 4
    final_classes = 4
    model_path = os.path.join(os.getcwd(), 'best_model')
    model_path = os.path.join(model_path, indicate)
    if indicate == "REAL":
        model_path = os.path.join(model_path, 'Fold' + str(fold))
    p2_path = os.path.join(model_path, 'Phase2_Model.pkl')
    # actually phase1 Model
    p2_model = torch.load(p2_path)
    p2_model = p2_model.cuda()
    p2_model = nn.DataParallel(p2_model, device_ids=None)
    p2_model.eval()
    Neighbor_count = 7
    tmp_result_path = os.path.join(save_path, map_name + "_phase2pred.txt")
    if os.path.exists(tmp_result_path):

        return tmp_result_path
    Phase2_Input = []
    Key_List = []
    Phase2_Final_Prediction_Dict = {}
    for key in Phase1_Final_Prediction_Dict.keys():
        tmp_values = key.split(',')
        x = int(tmp_values[0])
        y = int(tmp_values[1])
        z = int(tmp_values[2])
        input_example = np.zeros([n_classes, Neighbor_count, Neighbor_count, Neighbor_count])
        for i in range(Neighbor_count):
            for j in range(Neighbor_count):
                for k in range(Neighbor_count):
                    cur_x = x + i - int(np.floor(Neighbor_count / 2))
                    cur_y = y + j - int(np.floor(Neighbor_count / 2))
                    cur_z = z + k - int(np.floor(Neighbor_count / 2))
                    tmp_refer_key = str(cur_x) + ',' + str(cur_y) + ',' + str(cur_z) + ','
                    if tmp_refer_key in Phase1_Final_Prediction_Dict.keys():
                        tmp_phase2_input = Phase1_Final_Prediction_Dict[
                            tmp_refer_key]  # a num_Class dimension vector
                        for tmp_cls in range(n_classes):
                            input_example[tmp_cls, i, j, k] = tmp_phase2_input[tmp_cls]
        Phase2_Input.append(input_example)
        Key_List.append(key)
        if len(Phase2_Input) == batch_size:
            input_use = np.array(Phase2_Input)
            input_data = torch.from_numpy(input_use).float()
            # print(input_data.size())
            input_data = Variable(input_data)
            with torch.no_grad():
                output, p1, p2 = p2_model(input_data)
                output = F.softmax(output, dim=1)
            cb_array = output.cpu()
            cb_array = cb_array.detach().numpy()
            for i in range(len(Key_List)):
                tmp_key = Key_List[i]
                tmp_phase2_input = cb_array[i]
                Phase2_Final_Prediction_Dict[tmp_key] = tmp_phase2_input
            Phase2_Input = []
            Key_List = []
    if len(Phase2_Input) > 0:
        input_use = np.array(Phase2_Input)
        input_data = torch.from_numpy(input_use).float()
        # print(input_data.size())
        input_data = Variable(input_data)
        with torch.no_grad():
            output, p1, p2 = p2_model(input_data)
            output = F.softmax(output, dim=1)
        cb_array = output.cpu()
        cb_array = cb_array.detach().numpy()
        for i in range(len(Key_List)):
            tmp_key = Key_List[i]
            tmp_phase2_input = cb_array[i]
            Phase2_Final_Prediction_Dict[tmp_key] = tmp_phase2_input
        print("phase 2 final predictions:", tmp_phase2_input)
        Phase2_Input = []
        Key_List = []
    with open(tmp_result_path, 'w') as file:
        for tmp_key in Phase2_Final_Prediction_Dict:
            cur_prob = Phase2_Final_Prediction_Dict[tmp_key]
            cur_pred = int(np.argmax(cur_prob))
            file.write(str(tmp_key) + " " + str(int(cur_pred))+" ")
            #write all the probability values in the
            for prob in cur_prob:
                file.write(str(prob)+" ")
            file.write("\n")
            file.flush()
    Final_P2T_Dict={}
    #then forbid some transitions as we stated in the paper.
    tmp_result_path = os.path.join(save_path, map_name + "_phase2predT.txt")
    with open(tmp_result_path, 'w') as file:
        for key in Phase1_Final_Prediction_Dict.keys():
            pred1_prob = Phase1_Final_Prediction_Dict[key]
            pred2_prob = Phase2_Final_Prediction_Dict[key]
            pred1_label=int(np.argmax(pred1_prob))
            pred2_label = int(np.argmax(pred2_prob))
            # label1=true_dict2[key]
            # label2=true_dict1[key]
            if pred1_label == 1:  # beta to any
                pred2_label = pred2_label
            elif pred2_label == 3:  # any to DNA/RNA
                pred2_label = pred2_label
            else:
                pred2_label = pred1_label
            # get label through stride and trimmap considering here we do not have all label info
            file.write(str(key) + " " +str(int(pred2_label))+" ")
            for prob in pred2_prob:
                file.write(str(prob)+" ")
            file.write("\n")
            file.flush()
    return tmp_result_path
