
import os
import mrcfile
import numpy as np
from multiprocessing import Pool, Lock
from multiprocessing.sharedctypes import Value, Array
#global size,data,iterator,s
from numba import jit
#https://github.com/jglaser/interp3d
import shutil


@jit(nopython=True,nogil=True)
def interpolate_fast(data,data_new,size,iterator1,iterator2,iterator3,prev_voxel_size):
    for i in range(1, iterator1, 1):
        if(i%10==0):
            print("Finished %d",i)
        for j in range(1, iterator2, 1):
            for k in range(1, iterator3, 1):
                count = [int(i / prev_voxel_size), int(j / prev_voxel_size), int(k / prev_voxel_size)]
                e1 = count[0] + 1
                e2 = count[1] + 1
                e3 = count[2] + 1

                if (count[0] >= size[0] - 1):  # or count[1]>=size[1]-1 or count[2]>=size[2]-1 ):
                    # print(count)
                    e1 = count[0]
                    continue
                if (count[1] >= size[1] - 1):
                    e2 = count[1]
                    continue
                if (count[2] >= size[2] - 1):
                    e3 = count[2]
                    continue
                diff1 = [i - count[0] * prev_voxel_size, j - count[1] *prev_voxel_size, k - count[2] * prev_voxel_size]
                diff2=[e1*prev_voxel_size-i,e2*prev_voxel_size-j,e3*prev_voxel_size-k]
                # print(diff)
                val1 = data[count[0], count[1], count[2]]
                val2 = data[e1, count[1], count[2]]
                val3 = data[e1, e2, count[2]]
                val4 = data[count[0], e2, count[2]]
                val5 = data[count[0], count[1], e3]
                val6 = data[e1, count[1], e3]
                val7 = data[e1, e2, e3]
                val8 = data[count[0], e2, e3]
                #val = (val1 + diff[0] * (val2 - val1) + diff[1] * (val4 - val1) + diff[2] * (val5 - val1) + diff[0] *
                #       diff[1] * (val1 - val2 + val3 - val4) + diff[0] * diff[2] * (val1 - val2 - val5 + val6) + diff[
                #           1] * diff[2] * (
                #               val1 - val4 - val5 + val8) + diff[0] * diff[1] * diff[2] * (
                #               val8 - val7 + val6 - val5 + val4 - val3 + val2 - val1))
                u1=diff1[0]
                u2=diff2[0]
                v1=diff1[1]
                v2=diff2[1]
                w1=diff1[2]
                w2=diff2[2]
                val=(w2*(v1*(u1*val3+u2*val4)+v2*(u1*val2+u2*val1))+w1*(v1*(u1*val7+u2*val8)+v2*(u1*val6+u2*val5)))/(w1+w2)/(v1+v2)/(u1+u2)
                data_new[i,j,k]=val
    return data_new#np.float32(data_new)

@jit(nopython=True,nogil=True)
def interpolate_fast_general(data,data_new,size,iterator1,iterator2,iterator3,
                             prev_voxel_size1,prev_voxel_size2,prev_voxel_size3):
    for i in range(1, iterator1, 1):
        if(i%10==0):
            print("Finished %d",i)
        for j in range(1, iterator2, 1):
            for k in range(1, iterator3, 1):
                count = [int(i / prev_voxel_size1), int(j / prev_voxel_size2), int(k / prev_voxel_size3)]
                e1 = count[0] + 1
                e2 = count[1] + 1
                e3 = count[2] + 1

                if (count[0] >= size[0] - 1):  # or count[1]>=size[1]-1 or count[2]>=size[2]-1 ):
                    # print(count)
                    e1 = count[0]
                    continue
                if (count[1] >= size[1] - 1):
                    e2 = count[1]
                    continue
                if (count[2] >= size[2] - 1):
                    e3 = count[2]
                    continue
                diff1 = [i - count[0] * prev_voxel_size1, j - count[1] *prev_voxel_size2, k - count[2] * prev_voxel_size3]
                diff2 = [e1*prev_voxel_size1-i,e2*prev_voxel_size2-j,e3*prev_voxel_size3-k]
                # print(diff)
                val1 = data[count[0], count[1], count[2]]
                val2 = data[e1, count[1], count[2]]
                val3 = data[e1, e2, count[2]]
                val4 = data[count[0], e2, count[2]]
                val5 = data[count[0], count[1], e3]
                val6 = data[e1, count[1], e3]
                val7 = data[e1, e2, e3]
                val8 = data[count[0], e2, e3]
                #val = (val1 + diff[0] * (val2 - val1) + diff[1] * (val4 - val1) + diff[2] * (val5 - val1) + diff[0] *
                #       diff[1] * (val1 - val2 + val3 - val4) + diff[0] * diff[2] * (val1 - val2 - val5 + val6) + diff[
                #           1] * diff[2] * (
                #               val1 - val4 - val5 + val8) + diff[0] * diff[1] * diff[2] * (
                #               val8 - val7 + val6 - val5 + val4 - val3 + val2 - val1))
                u1=diff1[0]
                u2=diff2[0]
                v1=diff1[1]
                v2=diff2[1]
                w1=diff1[2]
                w2=diff2[2]
                val=(w2*(v1*(u1*val3+u2*val4)+v2*(u1*val2+u2*val1))+w1*(v1*(u1*val7+u2*val8)+v2*(u1*val6+u2*val5)))/(w1+w2)/(v1+v2)/(u1+u2)
                data_new[i,j,k]=val
    return data_new#np.float32(data_new)

def interpolate_slow(data,data_new,size,iterator1,iterator2,iterator3,prev_voxel_size):
    for i in range(1, iterator1, 1):
        if(i%10==0):
            print("Finished %d",i)
        for j in range(1, iterator2, 1):
            for k in range(1, iterator3, 1):
                count = [int(i / prev_voxel_size), int(j / prev_voxel_size), int(k / prev_voxel_size)]
                e1 = count[0] + 1
                e2 = count[1] + 1
                e3 = count[2] + 1

                if (count[0] >= size[0] - 1):  # or count[1]>=size[1]-1 or count[2]>=size[2]-1 ):
                    # print(count)
                    e1 = count[0]
                    continue
                if (count[1] >= size[1] - 1):
                    e2 = count[1]
                    continue
                if (count[2] >= size[2] - 1):
                    e3 = count[2]
                    continue
                diff1 = [i - count[0] * prev_voxel_size, j - count[1] *prev_voxel_size, k - count[2] * prev_voxel_size]
                diff2=[e1*prev_voxel_size-i,e2*prev_voxel_size-j,e3*prev_voxel_size-k]
                # print(diff)
                val1 = data[count[0], count[1], count[2]]
                val2 = data[e1, count[1], count[2]]
                val3 = data[e1, e2, count[2]]
                val4 = data[count[0], e2, count[2]]
                val5 = data[count[0], count[1], e3]
                val6 = data[e1, count[1], e3]
                val7 = data[e1, e2, e3]
                val8 = data[count[0], e2, e3]
                #val = (val1 + diff[0] * (val2 - val1) + diff[1] * (val4 - val1) + diff[2] * (val5 - val1) + diff[0] *
                #       diff[1] * (val1 - val2 + val3 - val4) + diff[0] * diff[2] * (val1 - val2 - val5 + val6) + diff[
                #           1] * diff[2] * (
                #               val1 - val4 - val5 + val8) + diff[0] * diff[1] * diff[2] * (
                #               val8 - val7 + val6 - val5 + val4 - val3 + val2 - val1))
                u1=diff1[0]
                u2=diff2[0]
                v1=diff1[1]
                v2=diff2[1]
                w1=diff1[2]
                w2=diff2[2]
                val=(w2*(v1*(u1*val3+u2*val4)+v2*(u1*val2+u2*val1))+w1*(v1*(u1*val7+u2*val8)+v2*(u1*val6+u2*val5)))/(w1+w2)/(v1+v2)/(u1+u2)
                data_new[i,j,k]=val
    return data_new



def Reform_Map_Voxel_Final(map_path,new_map_path):
    from scipy.interpolate import RegularGridInterpolator
    if not os.path.exists(new_map_path):
        with mrcfile.open(map_path,permissive=True) as mrc:
            prev_voxel_size = mrc.voxel_size
            prev_voxel_size_x = float(prev_voxel_size['x'])
            prev_voxel_size_y = float(prev_voxel_size['y'])
            prev_voxel_size_z = float(prev_voxel_size['z'])
            nx, ny, nz, nxs, nys, nzs, mx, my, mz = \
                mrc.header.nx, mrc.header.ny, mrc.header.nz, \
                mrc.header.nxstart, mrc.header.nystart, mrc.header.nzstart, \
                mrc.header.mx, mrc.header.my, mrc.header.mz
            orig = mrc.header.origin
            print("Origin:", orig)
            print("Previous voxel size:", prev_voxel_size)
            print("nx, ny, nz",nx,ny,nz)
            print("nxs,nys,nzs",nxs,nys,nzs)
            print("mx,my,mz",mx,my,mz)
            data = mrc.data
            data = np.swapaxes(data, 0, 2)
            size = np.shape(data)
            x = np.arange(size[0])
            y = np.arange(size[1])
            z = np.arange(size[2])
            my_interpolating_function = RegularGridInterpolator((x, y, z), data)
            it_val_x = int(np.floor(size[0] * prev_voxel_size_x))
            it_val_y = int(np.floor(size[1] * prev_voxel_size_y))
            it_val_z = int(np.floor(size[2] * prev_voxel_size_z))
            print("Previouse size:", size, " Current map size:", [it_val_x,it_val_y,it_val_z])
            data_new = np.zeros([it_val_x, it_val_y, it_val_z])
            for i in range(it_val_x):
                if i%10==0:
                    print("Resizing finished %d/%d"%(i,it_val_x))
                for j in range(it_val_y):
                    for k in range(it_val_z):
                        if ( i/prev_voxel_size_x>= size[0] - 1):
                            x_query=size[0] - 1
                        else:
                            x_query = i/prev_voxel_size_x

                        if ( j/prev_voxel_size_y>= size[1] - 1):
                            y_query=size[1] - 1
                        else:
                            y_query = j/prev_voxel_size_y
                        if ( k/prev_voxel_size_z>= size[2] - 1):
                            z_query=size[2] - 1
                        else:
                            z_query = k/prev_voxel_size_z
                        current_query=np.array([x_query,y_query,z_query])
                        current_value=float(my_interpolating_function(current_query))
                        data_new[i,j,k]=current_value

            data_new = np.swapaxes(data_new, 0, 2)
            data_new = np.float32(data_new)
            mrc_new = mrcfile.new(new_map_path, data=data_new, overwrite=True)
            vsize = mrc_new.voxel_size
            vsize.flags.writeable = True
            vsize.x = 1.0
            vsize.y = 1.0
            vsize.z = 1.0
            mrc_new.voxel_size = vsize
            mrc_new.update_header_from_data()
            mrc_new.header.nxstart = nxs * prev_voxel_size_x
            mrc_new.header.nystart = nys * prev_voxel_size_y
            mrc_new.header.nzstart = nzs * prev_voxel_size_z
            mrc_new.header.mapc = mrc.header.mapc
            mrc_new.header.mapr = mrc.header.mapr
            mrc_new.header.maps = mrc.header.maps
            mrc_new.header.origin = orig
            mrc_new.update_header_stats()
            mrc.print_header()
            mrc_new.print_header()
            mrc_new.close()
            mrc.close()
            del data
            del data_new
    return new_map_path


def Reform_Map_Voxel(map_path,new_map_path):
    if not os.path.exists(new_map_path):
        with mrcfile.open(map_path,permissive=True) as mrc:
            prev_voxel_size=mrc.voxel_size
            #assert len(prev_voxel_size)==3

            if not(prev_voxel_size['x']==prev_voxel_size['y'] and prev_voxel_size['x']==prev_voxel_size['z']):
                mrc.close()
                print("Grid size of different axis is different, resizing automatically.")
                return Reform_Map_Voxel_Final(map_path, new_map_path)
            prev_voxel_size=float(prev_voxel_size['x'])
            nx, ny, nz, nxs, nys, nzs, mx, my, mz =\
                mrc.header.nx, mrc.header.ny, mrc.header.nz, \
                mrc.header.nxstart, mrc.header.nystart, mrc.header.nzstart,\
                mrc.header.mx, mrc.header.my, mrc.header.mz
            orig = mrc.header.origin
            print("Origin:",orig)
            print("Previous voxel size:",prev_voxel_size)
            data = mrc.data
            data = np.swapaxes(data, 0, 2)
            size = np.shape(data)
            if (prev_voxel_size==1):
                shutil.copy(map_path,new_map_path)
                return new_map_path
            if (prev_voxel_size < 1):
                mrc.close()
                print("Grid size is smaller than 1, resizing automatically.")
                return Reform_Map_Voxel_Final(map_path, new_map_path)
            it_val1 = int(np.floor(size[0] * prev_voxel_size))
            it_val2 = int(np.floor(size[1] * prev_voxel_size))
            it_val3 = int(np.floor(size[2] * prev_voxel_size))
            print("Previouse size:",size," Current map size:",it_val1,it_val2,it_val3)
            data_new = np.zeros([it_val1,it_val2,it_val3])
            data_new[0, 0, 0] = data[0, 0, 0]
            data_new[it_val1 - 1, it_val2 - 1, it_val3 - 1] = data[
                size[0] - 1, size[1] - 1, size[2] - 1]
            #iterator = Value('i', it_val)
            #s = Value('d', float(prev_voxel_size))
            #pool = Pool(3)
            #out_1d = pool.map(interpolate,enumerate(np.reshape(data_new, (iterator.value * iterator.value * iterator.value,))))
            #data_new = np.array(out_1d).reshape(iterator.value, iterator.value, iterator.value)
            try:
                data_new=interpolate_fast(data,data_new,size,it_val1,it_val2,it_val3,prev_voxel_size)
            except:
                data_new = np.zeros([it_val1, it_val2, it_val3])
                data_new[0, 0, 0] = data[0, 0, 0]
                data_new[it_val1 - 1, it_val2 - 1, it_val3 - 1] = data[
                    size[0] - 1, size[1] - 1, size[2] - 1]
                data_new = interpolate_slow(data, data_new, size, it_val1,it_val2,it_val3, prev_voxel_size)
            data_new = np.swapaxes(data_new, 0, 2)
            data_new=np.float32(data_new)
            mrc_new = mrcfile.new(new_map_path, data=data_new, overwrite=True)
            vsize = mrc_new.voxel_size
            vsize.flags.writeable = True
            vsize.x = 1.0
            vsize.y = 1.0
            vsize.z = 1.0
            mrc_new.voxel_size = vsize
            mrc_new.update_header_from_data()
            mrc_new.header.nxstart = nxs * prev_voxel_size
            mrc_new.header.nystart = nys * prev_voxel_size
            mrc_new.header.nzstart = nzs * prev_voxel_size
            mrc_new.header.mapc = mrc.header.mapc
            mrc_new.header.mapr = mrc.header.mapr
            mrc_new.header.maps = mrc.header.maps
            mrc_new.header.origin = orig
            mrc_new.update_header_stats()
            mrc.print_header()
            mrc_new.print_header()
            mrc_new.close()
            mrc.close()
            del data
            del data_new
           # del out_1d
    return new_map_path
