

import os


def Build_label_file(save_path, map_name, input_path):
    tmp_result_path = os.path.join(save_path, map_name + "_label.txt")
    input_size=11
    with open(input_path, 'r') as file:
        with open(tmp_result_path,'w') as wfile:
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
                wfile.write(tmp_position+" ")
                label_info=str(result[3])
                label_info=label_info.split(":")
                first_correct=label_info[0]
                wfile.write(str(first_correct)+"\n")
                line = file.readline()
    return tmp_result_path