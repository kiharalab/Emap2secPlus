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


import argparse

def argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-F',type=str, required=True,help='map path')#File path for our MAINMAST code
    parser.add_argument('--mode',type=int,required=True,help='0: Predict structures for EM MAP\n'
                                                             '1: Predict structures for EM maps with pdb structure\n'
                                                             '2: Predict structure for experimental maps with 4 models\n'
                                                             '3: Predict and evaluate structure for experimental maps with 4 models\n' )
    parser.add_argument('--resize',type=int,default=0,help="0: resizing maps with numba optimized (some maps size are not supported);\n"
                                                           "1: resizing maps with scipy (relatively slow but support almost all maps).")
    parser.add_argument('-P',type=str,default="",help="PDB path for evaluating Model's performance")
    parser.add_argument('-M',type=str,default="best_model",help="Trained model path which saved all the trained models")
    parser.add_argument('--type',type=int,help='0:simulated map at 6 A 1: simulated map at 10 A 2: simulated map at 6-10 A 3:experimental map')
    parser.add_argument('--gpu',type=str,default='0',help='gpu id choose for training')
    parser.add_argument('--class', type=int, default='4', help='number of classes')
    parser.add_argument('--batch_size', type=int, default=256, help='batch size for training')
    parser.add_argument('--contour', type=float, default=0.0, help='Contour level for real map')
    parser.add_argument('--cardinality', default=32, type=int, help='ResNeXt cardinality')
    parser.add_argument('--drop_rate',type=float,default=0.3,help="Drop out rate for the phase2 Model")
    parser.add_argument('--fold',type=int,default=1,help='specify the fold Model used for predicting the real map')
    parser.add_argument('--output_folder',type=str,default="",help='specify a custom folder where results will be stored')
    parser.add_argument('--no_compilation', action='store_false', help='skip automatic compilation before execution')
    args = parser.parse_args()
    # try:
    #     import ray,socket
    #     rayinit()
    # except:
    #     print('ray need to be installed')#We do not need this since GAN can't be paralleled.
    params = vars(args)
    return params
