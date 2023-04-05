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
import sys
import torch
import numpy as np
from torch.autograd import Variable
from torch import nn
import torch.nn.functional as F
from scipy.special import softmax

# Importing error codes from main
current_dir = os.path.dirname(os.path.abspath(__file__))
main_dir = os.path.join(current_dir, '../')
sys.path.append(main_dir)
from main import ErrorCodes

def controlExecution(function, *args):
    """
    This function attepts to run the given function with arguments, and
    cotrols the possible exception ocurred during such execution.
    Color red is added to the exception message.
    """
    separator = ErrorCodes.CODE_MESSAGE_SEPARATOR.value
    try:
        return function(*args)
    except RuntimeError as rte:
        # Capturing RuntimeError (most likely CUDA out of memory)
        raise RuntimeError(str(ErrorCodes.CUDA_OUT_OF_MEMORY_CODE.value) + separator + "\033[31m" + str(rte) + "\033[0m")
    except Exception as e:
        # Capturing rest of exceptions
        raise Exception(str(ErrorCodes.DEFULT_ERROR_CODE.value) + separator + "\033[31m" + str(e) + "\033[0m")

def Predict_Phase1(save_path,map_name,input_path,indicate,fold,batch_size,params):
    input_size = 11
    dimension = input_size
    n_classes = 4
    #model_path=os.path.join(os.getcwd(),'best_model')
    model_path = os.path.abspath(params['M'])
    model_path=os.path.join(model_path,indicate)
    if indicate=="REAL":
        model_path=os.path.join(model_path,'Fold'+str(fold))
    cmodel_path = os.path.join(model_path, "Coil_Model.pkl")
    coil_model = torch.load(cmodel_path)
    coil_model = coil_model.cuda()
    coil_model = nn.DataParallel(coil_model, device_ids=None)
    coil_model.eval()
    # beta Model
    bmodel_path = os.path.join(model_path, "Beta_Model.pkl")
    beta_model = torch.load(bmodel_path)
    beta_model = beta_model.cuda()
    beta_model = nn.DataParallel(beta_model, device_ids=None)
    beta_model.eval()
    # alpha Model
    amodel_path = os.path.join(model_path, "Alpha_Model.pkl")
    alpha_model = torch.load(amodel_path)
    alpha_model = alpha_model.cuda()
    alpha_model = nn.DataParallel(alpha_model, device_ids=None)
    alpha_model.eval()
    # drna Model
    dmodel_path = os.path.join(model_path, "DRNA_Model.pkl")
    drna_model = torch.load(dmodel_path)
    drna_model = drna_model.cuda()
    drna_model = nn.DataParallel(drna_model, device_ids=None)
    drna_model.eval()
    #multi class Model
    combine_path = os.path.join(model_path, "All_Model.pkl")
    all_model = torch.load(combine_path)
    all_model = all_model.cuda()
    all_model = nn.DataParallel(all_model, device_ids=None)
    all_model.eval()

    combine_path = os.path.join(model_path, 'Combine_Modelall.pkl')
    # actually phase1 Model
    cb_model = torch.load(combine_path)
    cb_model = cb_model.cuda()
    cb_model = nn.DataParallel(cb_model, device_ids=None)
    cb_model.eval()

    Neighbor_count = 7
    tmp_result_path0 = os.path.join(save_path, map_name + "_phase1pred.txt")
    tmp_result_path = os.path.join(save_path, map_name + "_phase2pred.txt")
    if os.path.exists(tmp_result_path) and os.path.exists(tmp_result_path):
        #phase 1 prediction dict for phase 2 using
        Prediction_Dict={}
        with open(tmp_result_path,'r') as file:
            line=file.readline()
            while line:
                line=line.strip()
                split_result=line.split()
                tmp_key=split_result[0]
                tmp_pred=[]
                for k in range(n_classes):
                    tmp_pred.append(float(split_result[k+2]))
                Prediction_Dict[tmp_key]=np.array(tmp_pred)
                line=file.readline()

        return Prediction_Dict,tmp_result_path,tmp_result_path0
    Input_file = input_path
    Prediction_Dict = {}  # record the prediction with a dict
    #Output_Dict = {}
    Position_Order_List = []  # used to calculate the
    with open(Input_file, 'r') as file:
        line = file.readline()
        count = 0
        input_use = np.zeros([batch_size, 1, input_size, input_size, input_size])
        position_list = []
        outout_list = []
        batch_count=0
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
            position_list.append(tmp_position)
            outout_list.append(str(result[3]))  # real label,it may include two labels at the same time
            # print()

            for k in range(4, len(result)):
                index = k - 4
                x = index % (dimension)
                y = int(((index - x) % (dimension ** 2)) / dimension)
                z = int((index - x - y * dimension) / (dimension ** 2))
                input_use[count, 0, x, y, z] = float(result[k])

            count += 1
            if count == batch_size:
                batch_count+=1
                if batch_count%10==0:
                    print("Step 1 Phase1 Predicted %d batches, in total %d voxels"%(batch_count,batch_count*batch_size))
                # predict when a batch formed
                input_data = torch.from_numpy(input_use).float()
                # print(input_data.size())
                input_data = Variable(input_data)
                with torch.no_grad():
                    output = controlExecution(coil_model, input_data)[0]
                    # print(output.size())
                    output = F.softmax(output, dim=1)
                coil_array = output.cpu()
                coil_array = coil_array.detach().numpy()
                with torch.no_grad():
                    output = beta_model(input_data)[0]
                    # print(output.size())
                    output = F.softmax(output, dim=1)
                beta_array = output.cpu()
                beta_array = beta_array.detach().numpy()

                with torch.no_grad():
                    output = alpha_model(input_data)[0]
                    # print(output.size())
                    output = F.softmax(output, dim=1)
                alpha_array = output.cpu()
                alpha_array = alpha_array.detach().numpy()

                with torch.no_grad():
                    output = drna_model(input_data)[0]
                    # print(output.size())
                    output = F.softmax(output, dim=1)
                drna_array = output.cpu()
                drna_array = drna_array.detach().numpy()
                with torch.no_grad():
                    output, p1, p2 = all_model(input_data)
                    # print(output.size())
                    output = F.softmax(output, dim=1)
                output2array = output.cpu()
                output2array = output2array.detach().numpy()
                for i in range(len(position_list)):
                    tmp_key = position_list[i]
                    tmp_phase2_input = []
                    tmp_phase2_input.append(float(coil_array[i, 1]))
                    tmp_phase2_input.append(float(beta_array[i, 1]))
                    tmp_phase2_input.append(float(alpha_array[i, 1]))
                    tmp_phase2_input.append(float(drna_array[i, 1]))
                    tmp_all_output = output2array[i]
                    tmp_phase2_input += list(tmp_all_output)
                    Position_Order_List.append(tmp_key)
                    Prediction_Dict[tmp_key] = tmp_phase2_input
                    #Output_Dict[tmp_key] = tmp_phase2_output
                # print(output2array[0])
                position_list = []
                outout_list = []
                input_use = np.zeros([batch_size, 1, input_size, input_size, input_size])
                count = 0
            line = file.readline()
        # predict when a batch formed
        if len(position_list) > 0:
            input_use = input_use[0:len(position_list)]
            input_data = torch.from_numpy(input_use).float()
            input_data = Variable(input_data)
            with torch.no_grad():
                output, p1, p2 = coil_model(input_data)
                # print(output.size())
                output = F.softmax(output, dim=1)
            coil_array = output.cpu()
            coil_array = coil_array.detach().numpy()
            with torch.no_grad():
                output, p1, p2 = beta_model(input_data)
                # print(output.size())
                output = F.softmax(output, dim=1)
            beta_array = output.cpu()
            beta_array = beta_array.detach().numpy()

            with torch.no_grad():
                output, p1, p2 = alpha_model(input_data)
                # print(output.size())
                output = F.softmax(output, dim=1)
            alpha_array = output.cpu()
            alpha_array = alpha_array.detach().numpy()

            with torch.no_grad():
                output, p1, p2 = drna_model(input_data)
                # print(output.size())
                output = F.softmax(output, dim=1)
            drna_array = output.cpu()
            drna_array = drna_array.detach().numpy()
            with torch.no_grad():
                output, p1, p2 = all_model(input_data)
                # print(output.size())
                output = F.softmax(output, dim=1)
            output2array = output.cpu()
            output2array = output2array.detach().numpy()

            for i in range(len(position_list)):
                tmp_key = position_list[i]

                tmp_phase2_output = outout_list[i]
                tmp_phase2_input = []
                tmp_phase2_input.append(float(coil_array[i, 1]))
                tmp_phase2_input.append(float(beta_array[i, 1]))
                tmp_phase2_input.append(float(alpha_array[i, 1]))
                tmp_phase2_input.append(float(drna_array[i, 1]))
                tmp_all_output = output2array[i]
                tmp_phase2_input += list(tmp_all_output)
                Prediction_Dict[tmp_key] = tmp_phase2_input
                #Output_Dict[tmp_key] = tmp_phase2_output
                Position_Order_List.append(tmp_key)
        # print(output2array[0])
    print("after Phase1, we have %d predictions" % len(Prediction_Dict))
    Phase1_Input = []
    Key_List = []
    Phase1_Final_Prediction_Dict = {}
    batch_count=0
    # construct the input from Prediction_Dict for this purpose
    for key in Prediction_Dict.keys():
        tmp_values = key.split(',')
        x = int(tmp_values[0])
        y = int(tmp_values[1])
        z = int(tmp_values[2])
        input_example = np.zeros([n_classes * 2, Neighbor_count, Neighbor_count, Neighbor_count])
        for i in range(Neighbor_count):
            for j in range(Neighbor_count):
                for k in range(Neighbor_count):
                    cur_x = x + i - int(np.floor(Neighbor_count / 2))
                    cur_y = y + j - int(np.floor(Neighbor_count / 2))
                    cur_z = z + k - int(np.floor(Neighbor_count / 2))
                    tmp_refer_key = str(cur_x) + ',' + str(cur_y) + ',' + str(cur_z) + ','
                    if tmp_refer_key in Prediction_Dict.keys():
                        tmp_phase2_input = Prediction_Dict[tmp_refer_key]  # a num_Class dimension vector
                        for tmp_cls in range(n_classes * 2):
                            input_example[tmp_cls, i, j, k] = tmp_phase2_input[tmp_cls]
        # process the output

        Phase1_Input.append(input_example)
        Key_List.append(key)
        if len(Phase1_Input) == batch_size:
            batch_count+=1
            if batch_count % 10 == 0:
                print("Step 2 Phase1 Predicted %d batches, in total %d/%d voxels" % (
                batch_count, batch_count * batch_size,len(Prediction_Dict)))
            input_use = np.array(Phase1_Input)
            input_data = torch.from_numpy(input_use).float()
            # print(input_data.size())
            input_data = Variable(input_data)
            with torch.no_grad():
                output, p1, p2 = cb_model(input_data)
                output = F.softmax(output, dim=1)
            cb_array = output.cpu()
            cb_array = cb_array.detach().numpy()
            for i in range(batch_size):
                tmp_key = Key_List[i]
                tmp_phase2_input = cb_array[i]
                Phase1_Final_Prediction_Dict[tmp_key] = tmp_phase2_input

            Phase1_Input = []
            Key_List = []
    if len(Phase1_Input) > 0:
        input_use = np.array(Phase1_Input)
        input_data = torch.from_numpy(input_use).float()
        # print(input_data.size())
        input_data = Variable(input_data)
        with torch.no_grad():
            output, p1, p2 = cb_model(input_data)
            output = F.softmax(output, dim=1)
        cb_array = output.cpu()
        cb_array = cb_array.detach().numpy()
        for i in range(len(Key_List)):
            tmp_key = Key_List[i]
            tmp_phase2_input = cb_array[i]
            Phase1_Final_Prediction_Dict[tmp_key] = tmp_phase2_input

        Phase1_Input = []
        Key_List = []
    print("After Phase 2, we have %d predictions" % (
     len(Phase1_Final_Prediction_Dict)))
    with open(tmp_result_path0, 'w') as file:
        for tmp_key in Prediction_Dict.keys():
            cur_prob = Prediction_Dict[tmp_key]
            use_pred = []
            for k in range(n_classes):
                tmp_prob1 = cur_prob[k]
                tmp_prob2 = cur_prob[k + n_classes]
                use_pred.append(tmp_prob1 + tmp_prob2)
            use_pred = np.array(use_pred)
            cur_prob =softmax(use_pred)
            cur_pred = int(np.argmax(cur_prob))
            file.write(str(tmp_key) + " " + str(int(cur_pred)) + " ")
            # write all the probability values in the
            for prob in cur_prob:
                file.write(str(prob) + " ")
            file.write("\n")
            file.flush()
    # write a txt file
    with open(tmp_result_path, 'w') as file:
        for tmp_key in Phase1_Final_Prediction_Dict:
            cur_prob = Phase1_Final_Prediction_Dict[tmp_key]
            cur_pred = int(np.argmax(cur_prob))
            file.write(str(tmp_key) +  " " + str(int(cur_pred)) +" ")
            #write all the probability values in the
            for prob in cur_prob:
                file.write(str(prob)+" ")
            file.write("\n")
            file.flush()

    return Phase1_Final_Prediction_Dict,tmp_result_path,tmp_result_path0
