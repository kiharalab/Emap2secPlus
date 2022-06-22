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

def Build_Map(save_path,split_name,input_path,type,factor,contour_level,compile=True):
    """
    :param input_path: path for map
    :param type: simulated map or real map
    :return:
    """

    output_path=os.path.join(save_path,split_name+'.trimmap')
    #code path
    code_path0=os.path.join(os.getcwd(),'process_map')
    code_path=os.path.join(code_path0,'bin')
    ##make the file automatically unless skipping compilation is specified
    if compile:
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


def Build_Map_WithStructure(save_path,split_name,input_path,type,factor,contour_level,pdb_path,compile=True):
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
    ##make the file automatically unless skipping compilation is specified
    if compile:
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

