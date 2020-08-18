
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