

def Build_Refer_Loc_Dict(tmp_trimmap_path,factor):
    """
    mapping *.input location to real locations
    :param tmp_trimmap_path:
    :param factor:
    :return:
    """
    trimmap_input_location_refer_dict = {}  # map the input location to voxel id for DNA/RNA
    with open(tmp_trimmap_path, 'r') as file:
        line = file.readline()
        while line:
            if (line.rstrip() == ''):
                line = file.readline()
                continue
            if (line.startswith("-2") or line.startswith('#C: Res= -2')):
                line = file.readline()
                continue

            if (line.startswith("-1")):
                line = file.readline()
                continue
            if (line.startswith("#C: Res= -1")):
                line = file.readline()
                continue

            if (line.startswith('#C:')):
                current_label = int(line[8:].split()[0])
                split_list = line.split()
                bx = float(split_list[10])
                by = float(split_list[11])
                bz = float(split_list[12])
                Nstep = split_list[6]
                back = int((float(Nstep) - 1) / 2)  # reasonable considering griding edge fact
                cx = bx + back
                cy = by + back
                cz = bz + back
                tmp_coordinate = []
                tmp_coordinate.append(cx)
                tmp_coordinate.append(cy)
                tmp_coordinate.append(cz)

                # map input position to the count_drna
                real_pos = ""
                for k in range(14, 17):
                    real_pos += str(int(int(split_list[k]) / factor)) + ","
                trimmap_input_location_refer_dict[real_pos] = tmp_coordinate

                line = file.readline()

                continue
            if (line.startswith('#Base') or line.startswith('#Steps')):
                line = file.readline()
                continue

            if (line.startswith("#Voxel")):
                line = file.readline()
                continue
            if (line.startswith("#dmax")):
                line = file.readline()
                continue
            line = file.readline()
    return trimmap_input_location_refer_dict