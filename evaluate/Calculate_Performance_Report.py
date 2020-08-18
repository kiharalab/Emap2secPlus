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