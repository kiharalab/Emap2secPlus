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
from ops.argparser import argparser
from ops.os_operation import mkdir
import shutil
from enum import Enum

# Defining constants
RESULT_DEFAULT_FOLDER = 'Predict_Result'
RESULT_WITH_PDB_DEFAULT_FOLDER = 'Predict_Result_WithPDB'

# Error codes
class ErrorCodes(Enum):
    CODE_MESSAGE_SEPARATOR = "<-->" # Something that shouldn't be found in the exception string by chance
    DEFULT_ERROR_CODE = 1
    CUDA_OUT_OF_MEMORY_CODE = 2
    WRONG_PARAMETER_CODE = 3

def getOutputPath(custom_path, mode):
    """ This function returns the output path for the results. """
    if custom_path:
        return custom_path
    relative_default_path = RESULT_WITH_PDB_DEFAULT_FOLDER if (mode == 1 or mode == 3) else RESULT_DEFAULT_FOLDER
    return os.path.join(os.getcwd(), relative_default_path)

def splitException(exception):
    """
    This function extracts the error message and the error code from a custom error coded exception.
    Color red is added to the exception message.
    """
    exception = str(exception)
    splitIndex = exception.find(ErrorCodes.CODE_MESSAGE_SEPARATOR.value)
    if splitIndex == -1:
        return 1, ErrorCodes.DEFULT_ERROR_CODE.value, "\033[31m" + exception + "\033[0m"
    return int(exception[0:splitIndex]), "\033[31m" + exception[splitIndex + len(ErrorCodes.CODE_MESSAGE_SEPARATOR.value):] + "\033[0m"

def execute_with_exceptions(function, *args):
    """
    This function executes a function and captures the custom exception, prints it through terminal,
    and exits with the generated code.
    """
    try:
        return function(*args)
    except Exception as e:
        # Formatting and printing the exception
        code, message = splitException(e)
        print(message)
        exit(code)

if __name__ == "__main__":
    params = argparser()
    if params['mode']==0:
        input_map=params['F']
        input_map=os.path.abspath(input_map)
        type=params['type']
        choose = params['gpu']
        os.environ["CUDA_VISIBLE_DEVICES"] = choose
        if type==0:
            indicate='SIMU6'
        elif type==1:
            indicate='SIMU10'
        elif type==2:
            indicate= 'SIMU_MIX'
        elif type==3:
            indicate='REAL'
        else:
            print("we only have 4 type predictions: simulated(0,1,2) and experimental map(3)")
            exit(ErrorCodes.WRONG_PARAMETER_CODE)
        factor = 2  # reduce 4 to 2 to get more data
        save_path=getOutputPath(params["output_folder"], 0)
        mkdir(save_path)
        if not params['output_folder']:
            save_path = os.path.join(save_path, indicate)
            mkdir(save_path)
            fold = params['fold']  # specify use which fold Model based on real map
            if type==3:
                save_path = os.path.join(save_path, "Fold%d_Model_Result"%fold)
                mkdir(save_path)
            name_split=os.path.split(input_map)
            map_name=name_split[1]
            map_name=map_name.split(".")[0]
            save_path=os.path.join(save_path,map_name)
            mkdir(save_path)
        else:
            map_name="input"
        from process_map.Unify_Map import Unify_Map
        input_map = Unify_Map(input_map,os.path.join(save_path,map_name+"_unified.mrc"))
        # reform the map voxel size to 1A instead of experimental voxel size
        from process_map.Reform_Map_Voxel import Reform_Map_Voxel,Reform_Map_Voxel_Final
        output_map=os.path.join(save_path,map_name+".mrc")
        if type==3:
            input_map = Reform_Map_Voxel(input_map, output_map)
        else:
            shutil.copy(input_map,output_map)
        from process_map.Build_Map import Build_Map
        contour_level=params['contour']
        trimmap_path = Build_Map(save_path,map_name,input_map, type, factor, contour_level, compile=params["no_compilation"])
        #prepare input for using Model to predict
        from prepare_data.Prepare_Input import Prepare_Input
        input_path = Prepare_Input(save_path,map_name,trimmap_path,factor)
        #use the input to predict output

        batch_size=params['batch_size']
        from evaluate.Predict_Phase1 import Predict_Phase1
        phase1_pred_dict,phase1_pred_file,step1_pred_file=\
            execute_with_exceptions(Predict_Phase1, save_path,map_name,input_path,indicate,fold,batch_size,params)
        #visualize phase 1
        from evaluate.Visualize_Prediction import Visualize_Prediction,Visualize_Confident_Prediction
        Visualize_Prediction(save_path, map_name, step1_pred_file,  factor, 'Phase1')
        Visualize_Confident_Prediction(save_path, map_name, step1_pred_file, factor, 'Phase1')
        Visualize_Prediction(save_path,map_name,phase1_pred_file,factor,'Phase2')
        Visualize_Confident_Prediction(save_path,map_name,phase1_pred_file,factor,'Phase2')


    elif params['mode']==1:
        #predict maps with pdb structure and output predictions as well as voxel based accu and F1
        input_map = params['F']
        input_map = os.path.abspath(input_map)
        input_pdb=params['P']
        pdb_path=os.path.abspath(input_pdb)
        type = params['type']
        choose = params['gpu']
        os.environ["CUDA_VISIBLE_DEVICES"] = choose
        if type == 0:
            indicate = 'SIMU6'
        elif type == 1:
            indicate = 'SIMU10'
        elif type == 2:
            indicate = 'SIMU_MIX'
        elif type == 3:
            indicate = 'REAL'
        else:
            print("we only have 4 type predictions: simulated(0,1,2) and experimental map(3)")
            exit(ErrorCodes.WRONG_PARAMETER_CODE)
        factor = 2  # reduce 4 to 2 to get more data
        save_path = getOutputPath(params["output_folder"], 1)
        mkdir(save_path)
        if not params['output_folder']:
            save_path = os.path.join(save_path, indicate)
            mkdir(save_path)
            fold = params['fold']  # specify use which fold Model based on real map
            if type == 3:
                save_path = os.path.join(save_path, "Fold%d_Model_Result" % fold)
                mkdir(save_path)
            name_split = os.path.split(input_map)
            map_name = name_split[1]
            map_name = map_name.split(".")[0]
            save_path = os.path.join(save_path, map_name)
            mkdir(save_path)
        else:
            map_name="input"
        # reform the map voxel size to 1A instead of experimental voxel size
        from process_map.Unify_Map import Unify_Map
        input_map = Unify_Map(input_map,os.path.join(save_path,map_name+"_unified.mrc"))
        from process_map.Reform_Map_Voxel import Reform_Map_Voxel,Reform_Map_Voxel_Final

        output_map = os.path.join(save_path, map_name + ".mrc")
        if type==3:
            input_map = Reform_Map_Voxel(input_map, output_map)
        else:
            shutil.copy(input_map,output_map)
        from process_map.Build_Map import Build_Map_WithStructure
        contour_level = params['contour']
        trimmap_path = Build_Map_WithStructure(save_path, map_name, input_map, type, factor, contour_level,pdb_path, compile=params["no_compilation"])
        from prepare_data.Gen_Stride import Gen_Stride
        stride_path = Gen_Stride(save_path,map_name,pdb_path)
        from prepare_data.Prepare_Input import Prepare_Input_WithStructure
        input_path = Prepare_Input_WithStructure(save_path, map_name, trimmap_path, factor,pdb_path,stride_path)

        from evaluate.Visualize_Prediction_WithStructure import Visualize_Prediction_WithStructure, \
            Visualize_Confident_Prediction_WithStructure
        from prepare_data.Build_Refer_Loc_Dict import Build_Refer_Loc_Dict

        real_loc_refer = Build_Refer_Loc_Dict(trimmap_path, factor)
        from prepare_data.Build_label_file import Build_label_file
        real_label_path=Build_label_file(save_path, map_name, input_path)
        Visualize_Prediction_WithStructure(save_path, map_name, real_label_path, factor, real_loc_refer, "REAL")
        batch_size = params['batch_size']
        from evaluate.Predict_Phase1 import Predict_Phase1

        phase1_pred_dict, phase1_pred_file,step1_pred_file = \
            execute_with_exceptions(Predict_Phase1, save_path, map_name, input_path, indicate, fold, batch_size,params)

        Visualize_Prediction_WithStructure(save_path, map_name, step1_pred_file, factor, real_loc_refer, 'Phase1')
        Visualize_Confident_Prediction_WithStructure(save_path, map_name, step1_pred_file, factor, real_loc_refer,
                                                     'Phase1')
        Visualize_Prediction_WithStructure(save_path, map_name, phase1_pred_file, factor, real_loc_refer,'Phase2')
        Visualize_Confident_Prediction_WithStructure(save_path, map_name, phase1_pred_file, factor,real_loc_refer, 'Phase2')

        from evaluate.Calculate_Performance_Report import Calculate_Performance_Report
        Calculate_Performance_Report(save_path, map_name,input_path,step1_pred_file,'Phase1')
        Calculate_Performance_Report(save_path, map_name, input_path, phase1_pred_file, 'Phase2')
    elif params['mode']==2:
        #augmented predictions by all 4 models for experimental maps
        input_map = params['F']
        input_map = os.path.abspath(input_map)
        type = params['type']
        choose = params['gpu']
        os.environ["CUDA_VISIBLE_DEVICES"] = choose
        assert type==3
        indicate="REAL"
        factor = 2  # reduce 4 to 2 to get more data
        save_path0 = getOutputPath(params["output_folder"], 2)
        mkdir(save_path0)
        if not params['output_folder']:
            save_path0 = os.path.join(save_path0, indicate)
            mkdir(save_path0)
            name_split = os.path.split(input_map)
            map_name = name_split[1]
            map_name = map_name.split(".")[0]
            save_path0 = os.path.join(save_path0, map_name)
            mkdir(save_path0)
        else:
            map_name="input"
        from process_map.Unify_Map import Unify_Map
        input_map = Unify_Map(input_map,os.path.join(save_path0,map_name+"_unified.mrc"))
        from process_map.Reform_Map_Voxel import Reform_Map_Voxel, Reform_Map_Voxel_Final

        output_map = os.path.join(save_path0, map_name + ".mrc")
        if type==3:
            input_map = Reform_Map_Voxel(input_map, output_map)
        else:
            shutil.copy(input_map,output_map)

        from process_map.Build_Map import Build_Map

        contour_level = params['contour']
        trimmap_path = Build_Map(save_path0, map_name, input_map,type, factor, contour_level, compile=params["no_compilation"])
        # prepare input for using Model to predict
        from prepare_data.Prepare_Input import Prepare_Input

        input_path = Prepare_Input(save_path0, map_name, trimmap_path, factor)
        All_Output_File=[]
        for fold in range(1,5):
            save_path = os.path.join(save_path0, "Fold%d_Model_Result" % fold)
            mkdir(save_path)
            # reform the map voxel size to 1A instead of experimental voxel size
            # use the input to predict output

            batch_size = params['batch_size']
            from evaluate.Predict_Phase1 import Predict_Phase1

            phase1_pred_dict, phase1_pred_file,step1_pred_file = \
                execute_with_exceptions(Predict_Phase1, save_path, map_name, input_path, indicate, fold, batch_size,params)
            # visualize phase 1
            from evaluate.Visualize_Prediction import Visualize_Prediction, Visualize_Confident_Prediction

            Visualize_Prediction(save_path, map_name, step1_pred_file, factor, 'Phase1')
            Visualize_Confident_Prediction(save_path, map_name, phase1_pred_file, factor, 'Phase1')
            Visualize_Prediction(save_path, map_name, step1_pred_file, factor, 'Phase2')
            Visualize_Confident_Prediction(save_path, map_name, phase1_pred_file, factor, 'Phase2')
            All_Output_File.append(phase1_pred_file)
        from evaluate.Combine_All_Predictions import Combine_All_Predictions
        final_pred_path=Combine_All_Predictions(save_path0, map_name, All_Output_File)
        Visualize_Prediction(save_path0, map_name, final_pred_path, factor, 'Final')
        Visualize_Confident_Prediction(save_path0, map_name,final_pred_path, factor, 'Final')

    elif params['mode']==3:
        # predict maps with pdb structure and output predictions as well as voxel based accu and F1
        input_map = params['F']
        input_map = os.path.abspath(input_map)
        input_pdb = params['P']
        pdb_path = os.path.abspath(input_pdb)
        type = params['type']
        choose = params['gpu']
        os.environ["CUDA_VISIBLE_DEVICES"] = choose
        assert type == 3
        indicate = "REAL"
        factor = 2  # reduce 4 to 2 to get more data
        save_path0 = getOutputPath(params["output_folder"], 3)
        mkdir(save_path0)
        if not params['output_folder']:
            save_path0 = os.path.join(save_path0, indicate)
            mkdir(save_path0)
            name_split = os.path.split(input_map)
            map_name = name_split[1]
            map_name = map_name.split(".")[0]
            save_path0 = os.path.join(save_path0, map_name)
            mkdir(save_path0)
        else:
            map_name="input"
        from process_map.Unify_Map import Unify_Map
        input_map = Unify_Map(input_map,os.path.join(save_path0,map_name+"_unified.mrc"))
        from process_map.Reform_Map_Voxel import Reform_Map_Voxel, Reform_Map_Voxel_Final

        output_map = os.path.join(save_path0, map_name + ".mrc")
        if type==3:
            input_map = Reform_Map_Voxel(input_map, output_map)
        else:
            shutil.copy(input_map,output_map)
        from process_map.Build_Map import Build_Map_WithStructure

        contour_level = params['contour']
        trimmap_path = Build_Map_WithStructure(save_path0, map_name, input_map, type, factor, contour_level,
                                               pdb_path, compile=params["no_compilation"])
        from prepare_data.Gen_Stride import Gen_Stride

        stride_path = Gen_Stride(save_path0, map_name, pdb_path)
        from prepare_data.Prepare_Input import Prepare_Input_WithStructure

        input_path = Prepare_Input_WithStructure(save_path0, map_name, trimmap_path, factor, pdb_path, stride_path)
        from prepare_data.Build_Refer_Loc_Dict import Build_Refer_Loc_Dict

        real_loc_refer = Build_Refer_Loc_Dict(trimmap_path, factor)
        All_Output_File = []
        for fold in range(1,5):
            save_path = os.path.join(save_path0, "Fold%d_Model_Result" % fold)
            mkdir(save_path)
            # reform the map voxel size to 1A instead of experimental voxel size
            # use the input to predict output

            batch_size = params['batch_size']
            from evaluate.Predict_Phase1 import Predict_Phase1

            phase1_pred_dict, phase1_pred_file,step1_pred_file = \
                execute_with_exceptions(Predict_Phase1, save_path, map_name, input_path, indicate, fold,batch_size,params)
            # visualize phase 1
            from evaluate.Visualize_Prediction_WithStructure import Visualize_Prediction_WithStructure, \
                Visualize_Confident_Prediction_WithStructure

            Visualize_Prediction_WithStructure(save_path, map_name, step1_pred_file, factor, real_loc_refer, 'Phase1')
            Visualize_Confident_Prediction_WithStructure(save_path, map_name, step1_pred_file, factor, real_loc_refer,
                                                         'Phase1')
            Visualize_Prediction_WithStructure(save_path, map_name, phase1_pred_file, factor, real_loc_refer, 'Phase2')
            Visualize_Confident_Prediction_WithStructure(save_path, map_name, phase1_pred_file, factor, real_loc_refer,
                                                         'Phase2')
            from evaluate.Calculate_Performance_Report import Calculate_Performance_Report

            Calculate_Performance_Report(save_path, map_name, input_path, step1_pred_file, 'Phase1')
            Calculate_Performance_Report(save_path, map_name, input_path, phase1_pred_file, 'Phase2')
            All_Output_File.append(phase1_pred_file)

        from evaluate.Combine_All_Predictions import Combine_All_Predictions

        final_pred_path = Combine_All_Predictions(save_path0, map_name, All_Output_File)
        Visualize_Prediction_WithStructure(save_path0, map_name, final_pred_path, factor, real_loc_refer, 'Final')
        Visualize_Confident_Prediction_WithStructure(save_path0, map_name, final_pred_path, factor, real_loc_refer,
                                                     'Final')
        from prepare_data.Build_label_file import Build_label_file

        real_label_path = Build_label_file(save_path0, map_name, input_path)
        Visualize_Prediction_WithStructure(save_path0, map_name, real_label_path, factor, real_loc_refer, 'REAL')
        Calculate_Performance_Report(save_path0, map_name, input_path,final_pred_path, 'Final')

    elif params['mode']==4:
        #make protein+DNA/RNA structure binary prediction
        input_map = params['F']
        input_map = os.path.abspath(input_map)
        type = params['type']
        choose = params['gpu']
        os.environ["CUDA_VISIBLE_DEVICES"] = choose
        assert type == 3
        indicate = "REAL"
        factor = 2  # reduce 4 to 2 to get more data
        if params['output_folder']:
            name_split = os.path.split(input_map)
            map_name = name_split[1]
            map_name = map_name.split(".")[0]
            save_path0 = os.path.join(getOutputPath(params["output_folder"], 4), "Binary", indicate, map_name)
            mkdir(save_path0)
        else:
            save_path0 = getOutputPath(params["output_folder"], 4)
            mkdir(save_path0)
            map_name="input"
        from process_map.Unify_Map import Unify_Map
        input_map = Unify_Map(input_map,os.path.join(save_path0,map_name+"_unified.mrc"))
        from process_map.Reform_Map_Voxel import Reform_Map_Voxel, Reform_Map_Voxel_Final

        output_map = os.path.join(save_path0, map_name + ".mrc")
        input_map = Reform_Map_Voxel(input_map, output_map)

        from process_map.Build_Map import Build_Map

        contour_level = params['contour']
        trimmap_path = Build_Map(save_path0, map_name, input_map, type, factor, contour_level, compile=params["no_compilation"])

        from prepare_data.Prepare_Input import Prepare_Input

        input_path = Prepare_Input(save_path0, map_name, trimmap_path, factor)
        All_Output_File = []
        from evaluate.Visualize_Binary_Prediction import Visualize_Binary_Prediction, \
            Visualize_Binary_Confident_Prediction
        for fold in range(1, 5):
            save_path = os.path.join(save_path0, "Fold%d_Model_Result" % fold)
            mkdir(save_path)
            # reform the map voxel size to 1A instead of experimental voxel size
            # use the input to predict output

            batch_size = params['batch_size']
            from evaluate.Predict_Phase1 import Predict_Phase1

            phase1_pred_dict, phase1_pred_file, step1_pred_file = \
                execute_with_exceptions(Predict_Phase1, save_path, map_name, input_path,indicate, fold, batch_size,params)
            # visualize phase 1

            Visualize_Binary_Prediction(save_path, map_name, step1_pred_file, factor, 'Phase1')
            Visualize_Binary_Confident_Prediction(save_path, map_name, phase1_pred_file, factor, 'Phase1')
            Visualize_Binary_Prediction(save_path, map_name, step1_pred_file, factor, 'Phase2')
            Visualize_Binary_Confident_Prediction(save_path, map_name, phase1_pred_file, factor, 'Phase2')
            All_Output_File.append(phase1_pred_file)
        from evaluate.Combine_All_Predictions import Combine_All_Predictions

        final_pred_path = Combine_All_Predictions(save_path0, map_name, All_Output_File)
        Visualize_Binary_Prediction(save_path0, map_name, final_pred_path, factor, 'Final')
        Visualize_Binary_Confident_Prediction(save_path0, map_name, final_pred_path, factor, 'Final')
