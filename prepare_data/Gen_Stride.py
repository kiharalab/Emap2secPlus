
import os
from ops.os_operation import mkdir
def Gen_Stride(save_path,map_name,pdb_path):
    code_path=os.path.join(os.getcwd(),'ops')
    code_path=os.path.join(code_path,'stride')
    os.system("chmod 777 "+code_path)#give full permission to this software
    root_path=os.getcwd()
    os.chdir(save_path)
    output_path = os.path.join(save_path, map_name+'.stride')
    os.system("chmod 777 " + save_path)
    os.system("chmod 777 " + output_path)  # give full permission to this path
    output_name=map_name+'.stride'
    os.system(code_path+' -f'+output_name+' '+pdb_path)
    os.chdir(root_path)
    return output_path
