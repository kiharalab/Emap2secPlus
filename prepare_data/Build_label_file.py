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
