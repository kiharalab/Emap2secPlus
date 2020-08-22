# Publication:  "Emap2sec+:  Detecting Protein and DNA/RNA Structures in Cryo-EM Maps of Intermediate Resolution Using Deep Learning", Xiao Wang, Eman Alnabati, Tunde W. Aderinwale, Sai Raghavendra Maddhuri Venkata Subramaniya, Genki Terashi, and Daisuke Kihara, BioRxiv (2020)

# Emap2sec+ is a computational tool using deep learning that can accurately identify structures, alpha helices, beta sheets, other(coils/turns) and DNA/RNA, in cryo-Electron Microscopy (EM) maps of medium to low resolution.
# Copyright (C) 2020 Xiao Wang, Eman Alnabati, Tunde W Aderinwale, Sai Raghavendra Maddhuri, Genki Terashi, Daisuke Kihara, and Purdue University.
# License: GPL v3 for academic use. (For commercial use, please contact us for different licensing.)
# Contact: Daisuke Kihara (dkihara@purdue.edu)


#

# This program is free software: you can redistribute it and/or modify

# it under the terms of the GNU General Public License as published by

# the Free Software Foundation, either version 3 of the License, or

# (at your option) any later version.

#

# This program is distributed in the hope that it will be useful,

# but WITHOUT ANY WARRANTY; without even the implied warranty of

# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the

# GNU General Public License V3 for more details.

#

# You should have received a copy of the GNU v3.0 General Public License

# along with this program.  If not, see https://www.gnu.org/licenses/gpl-3.0.en.html.

def Extract_PDB_Refer(tmp_pdb_path):
    count_ca = 0
    Refer_Dict = {}
    Refer_Dict2={}
    with open(tmp_pdb_path, 'r') as file:
        line = file.readline()
        # residue_id=line[22:26]

        while line:
            if line[:4] == "ATOM":
                split_result=line.strip().split()
                atom_id = line[13:16]
                if atom_id == "CA ":
                    count_ca += 1
                    residue_id =str(line[22:27].strip()) #str(split_result[5])
                    chain_id=line[21]
                    residue_name=line[17:20]
                    tmp_key=chain_id+","+residue_name+","+str(residue_id)
                    Refer_Dict[count_ca]=tmp_key
                    Refer_Dict2[count_ca]=residue_name+","+str(residue_id)
            elif line[:6]=="ENDMDL":#for pdb with multiple models
                break
            line = file.readline()
    print("CA count %d, REFER DICT size %d" % (count_ca,len(Refer_Dict)))
    return Refer_Dict,Refer_Dict2
