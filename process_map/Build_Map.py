import os
from ops.os_operation import mkdir
def Build_Map(save_path,split_name,input_path,type,factor,contour_level):
    """
    :param input_path: path for map
    :param type: simulated map or real map
    :return:
    """

    output_path=os.path.join(save_path,split_name+'.trimmap')
    #code path
    code_path0=os.path.join(os.getcwd(),'process_map')
    code_path=os.path.join(code_path0,'bin')
    ##make the file automatically
    root_path=os.getcwd()
    os.chdir(code_path0)
    os.system("make clean")
    os.system("make")
    os.chdir(root_path)
    code_path=os.path.join(code_path,'map2train')
    if type==3:
        commandline = code_path + ' ' + input_path + '  -r 3.0 -c ' + str(
            contour_level) + ' -sstep ' + str(factor) + ' >' + output_path
    else:
        commandline = code_path + ' ' + input_path + ' -r 3.0 -c ' + str(
            contour_level) + ' -sstep ' + str(factor) + ' -ignorestart >' + output_path
    print("Extracting Trimmap:",commandline)#add ignorestart or not based on different maps
    os.system(commandline)
    return output_path


def Build_Map_WithStructure(save_path,split_name,input_path,type,factor,contour_level,pdb_path):
    """
    :param input_path: path for map
    :param type: simulated map or real map
    :return:
    """

    output_path=os.path.join(save_path,split_name+'.trimmap')
    if os.path.exists(output_path) and os.path.getsize(output_path)>100000:
        return output_path
    #code path
    code_path0=os.path.join(os.getcwd(),'process_map')
    code_path=os.path.join(code_path0,'bin')
    ##make the file automatically
    root_path=os.getcwd()
    os.chdir(code_path0)
    os.system("make clean")
    os.system("make")
    os.chdir(root_path)
    code_path=os.path.join(code_path,'map2train')
    if type==3:
        commandline = code_path + ' ' + input_path +' -P ' + pdb_path +  '  -r 3.0 -c ' + str(
            contour_level) + ' -sstep ' + str(factor) + ' >' + output_path
    else:
        commandline = code_path + ' ' + input_path + ' -P ' + pdb_path + ' -r 3.0 -c ' + str(
            contour_level) + ' -sstep ' + str(factor) + ' -ignorestart >' + output_path
    print("Extracting Trimmap:",commandline)#add ignorestart or not based on different maps
    os.system(commandline)
    return output_path

