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
