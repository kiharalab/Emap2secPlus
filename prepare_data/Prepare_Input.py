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
def Prepare_Input(save_path,map_name,trimmap_path,factor):

    output_path = os.path.join(save_path, map_name + '.input')
    if os.path.exists(output_path) and os.path.getsize(output_path)>1000:
        return output_path
    Gen_Secondary_Eval(trimmap_path,output_path,factor)
    return output_path

def Prepare_Input_WithStructure(save_path,map_name,trimmap_path,factor,pdb_path,stride_path):

    output_path = os.path.join(save_path, map_name + '.input')
    if os.path.exists(output_path) and os.path.getsize(output_path)>1000:
        return output_path
    from prepare_data.Extract_PDB_Refer import Extract_PDB_Refer
    PDB_refer_dict, PDB_refer_dict2 = Extract_PDB_Refer(pdb_path)#building dict for some cases stride output is strange and does not have the order
    Gen_Secondary_Mark(trimmap_path, stride_path, output_path, PDB_refer_dict, PDB_refer_dict2,
                       factor)  #generate input with labels then it's easier to evaluate
    return output_path


def Gen_Secondary_Mark(trimmap_path, stride_path, output_path, PDB_refer_dict, PDB_refer_dict2,
                       factor):
    fil1 = open(trimmap_path)  # Trimap generated by the Process_map code
    fil3 = open(output_path, 'w')
    fil4 = open(stride_path)
    lines1 = fil1.readlines()
    lines4 = fil4.readlines();

    # Get the secondary structure label from stride files
    Label_Dict = {}  #
    Label_Dict2 = {}  # do not include chain information, just locate based on very naive residue_id+residue_type
    j = 0;
    acount = 0
    bcount = 0
    ccount = 0
    use_chain_flag = True
    actual_acount=0
    actual_bcount=0
    actual_ccount=0
    for line in lines4:
        if (line.split()[0] == 'ASG'):
            residue_id = str(line.split()[3].strip())
            chain_id = line.split()[2]
            residue_name = line.split()[1]
            if chain_id == '-':
                use_chain_flag = False
            tmp_key = chain_id + "," + residue_name + "," + str(residue_id)
            tmp_key1 = residue_name + "," + str(residue_id)
            if tmp_key1 in Label_Dict2.keys() and use_chain_flag == False:
                print("can not be used found two same key in stride file %s/%s", tmp_key1, tmp_key)
                exit()
            if (line.split()[5] == 'H' or line.split()[5] == 'G' or line.split()[5] == 'I'):
                if tmp_key not in Label_Dict:
                    actual_acount+=1
                Label_Dict[tmp_key] = 2
                Label_Dict2[tmp_key1] = 2
                acount += 1
            elif (line.split()[5] == 'B' or line.split()[5] == 'E'):
                if tmp_key not in Label_Dict:
                    actual_bcount+=1
                Label_Dict[tmp_key] = 1
                Label_Dict2[tmp_key1] = 1
                bcount += 1
            else:
                if tmp_key not in Label_Dict:
                    actual_ccount+=1
                Label_Dict[tmp_key] = 0
                Label_Dict2[tmp_key1] = 0
                ccount += 1
    print("alhpa : %d, beta : %d, none : %d\n" % (acount, bcount, ccount))
    print("actual alhpa : %d, beta : %d, none : %d\n" % (actual_acount, actual_bcount, actual_ccount))
    if acount!=actual_acount or bcount!=actual_bcount or ccount!=actual_ccount:
        print("!!!Stride file exists problems!!!")
        print("Predictions can be trusted, while the evaluation results can't be trusted!!!")
    print("Residue %d " % len(Label_Dict))
    print("we used "+str(use_chain_flag)+" chain search dict")
    max_resi_id = len(Label_Dict)
    count = 0
    coords = []
    for line in lines1:
        flag = 0;
        if (line.rstrip() == ''):
            continue
        elif (line.startswith("-2") or line.startswith('#C: Res= -2')):
            continue
        elif (line.startswith('#C:')):
            equ = line.split('=')
            coords = equ[len(equ) - 1].split(' ')[1:4]
            prev_line = line
            # print(coords)
            continue

        if (line.startswith("-1")):
            fil3.write(str(int(int(coords[0]) / factor)) + "," + str(int(int(coords[1]) / factor)) + "," + str(
                int(int(coords[2].rstrip()) / factor)) + ',' + line)
            del coords
            # fil3.write(line)
            # fil3.write('\n')
            continue

        elif (line.startswith('#Base') or line.startswith('#Steps')):
            continue

        elif (line.startswith("#Voxel")):
            print('!!!!!!!!!1here!!!!!!!!!!!!!!!!!')
            fil3.write(line.split()[2] + "," + line.split()[3] + "," + line.split()[4])
            fil3.write('\n')
            continue
        elif (line.startswith("#dmax")):
            continue
        # debug info
        # if int(int(coords[0]) / factor)==28 and int(int(coords[1]) / factor)==17 and int(int(coords[2].rstrip()) / factor)==3:
        #     print(prev_line)
        #     print(line)
        #     continue
        # if (count == 1):
        #    print(line)
        # print(labelArray1[(int(line.split(',')[0])-1)])

        li = line.split(',');

        # print(len(labelArray1))
        for i in li:
            if (flag == 0):
                # First write label
                # print(int(i)-1);
                lbs = i.split(';')  # For the residue and RNA mark part, we use ; to divide them
                # ones = [0, 0, 0, 0]
                ones = []
                for j in lbs:
                    # print(lbs)
                    # print((int(j)-1))
                    j = int(j)
                    if j == -999:
                        ones.append(3)
                        continue
                    if j > max_resi_id:  # remove those target because some simulated maps have 2 structures which will cause problems.
                        #print("exceeding the maximum residue number!!!")
                        continue
                    if j not in PDB_refer_dict:
                        #print(lbs)
                        #print("%d residue no refer in pdb, in total %d residue" % (j, len(PDB_refer_dict)))
                        continue
                    if use_chain_flag == False:
                        tmp_key = PDB_refer_dict2[int(j)]
                        if tmp_key in Label_Dict2:
                            tmp_label = Label_Dict2[tmp_key]
                        else:
                            print("no predictions in stride for %s" % tmp_key)
                    else:
                        tmp_key = PDB_refer_dict[int(j)]
                        if tmp_key in Label_Dict:
                            tmp_label = Label_Dict[tmp_key]
                        else:
                            print("no predictions in stride for %s" % tmp_key)

                    ones.append(tmp_label)

                if (len(ones) == 0):
                    label = "-1"
                else:
                    # labelA = [str(k) for k, l in enumerate(ones) if l != 0]
                    labelA = [str(k) for k in ones]
                    label = ':'.join(labelA)

                fil3.write(str(int(int(coords[0]) / factor)) + "," + str(int(int(coords[1]) / factor)) + "," + str(
                    int(int(coords[2].rstrip()) / factor)) + ',')
                fil3.write(label)
                # print(label);
                flag = 1;

            else:
                # Then write the input
                # print(i);
                fil3.write(',');
                fil3.write(i);
        #    fil3.write('\n');
        del coords
        count += 1



def Gen_Secondary_Eval(dataFile,outputFile,factor):
    """
    :param dataFile: trimap file generated by our code in process_map
    :param strideFile: stride file generated by ops/stride, mark secondary structure
    :param outputFile:
    :param pdb_id:
    :param factor:
    :return:
    """
    fil1 = open(dataFile)#Trimap generated by the Process_map code
    fil3 = open(outputFile, 'w')

    lines1 = fil1.readlines()
    count = 0
    coords = []
    for line in lines1:
        if count%1000==0:
            print("Finishing %d/%d lines of trimmap to input"%(count,len(lines1)))
        if (line.rstrip() == ''):
            continue
        elif (line.startswith("-2") or line.startswith('#C: Res= -2')):
            continue
        elif (line.startswith('#C:')):
            equ = line.split('=')
            coords = equ[-1].split(' ')[1:4]#use real coord
            # print(coords)
            continue

        if (line.startswith("-1")):
            fil3.write(str(int(float(coords[0])/ factor)) + "," + str(int(float(coords[1])/ factor)) + "," + str(
                int(float(coords[2].rstrip())/ factor )) + ',' + line)
            # fil3.write(line)
            # fil3.write('\n')
            continue

        elif (line.startswith('#Base') or line.startswith('#Steps')):
            continue

        elif (line.startswith("#Voxel")):
            print('!!!!!!!!!1here!!!!!!!!!!!!!!!!!')
            fil3.write(line.split()[2] + "," + line.split()[3] + "," + line.split()[4])
            fil3.write('\n')
            continue
        elif (line.startswith("#dmax")):
            continue
        # print(labelArray1[(int(line.split(',')[0])-1)])
        li = line.split(',');
        flag = 0;
    # print(len(labelArray1))
        for i in li:
            if (flag == 0):
                #First write label
                # print(int(i)-1);
                label="0"
                fil3.write(str(int(float(coords[0]) / factor)) + "," + str(int(float(coords[1])/ factor )) + "," + str(
                int(float(coords[2].rstrip())/ factor )) + ',')
                fil3.write(label)
                # print(label);
                flag = 1;

            else:
                #Then write the input
                # print(i);
                fil3.write(',');
                fil3.write(i);
        #    fil3.write('\n');
        count += 1
