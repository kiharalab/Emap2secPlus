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

import torch
from torch import nn
from torch import optim
import math
from functools import partial
import torch.nn.functional as F
class Phase2_CNNMLP(nn.Module):

    def __init__(self,input_channel=4,num_classes=5,dropout_rate=0,neighbor=3):
        """
        :param num_classes:num of classes that we need to predict
        """
        super(Phase2_CNNMLP, self).__init__()
        self.conv1 = nn.Conv3d(
            input_channel,
            32,
            kernel_size=2,
            stride=(1, 1, 1),
            padding=0,
            bias=False)  # Out input channel size is 1, our map
        self.bn = nn.BatchNorm3d(32)
        self.relu = nn.ReLU(inplace=True)
        #then fully connected layer
        self.fc1 = nn.Linear(32*((neighbor-1)**3), 128)
        self.bn1 = nn.BatchNorm1d(128)
        self.relu1 = nn.ReLU(inplace=True)

        self.fc2=nn.Linear(128,64)
        self.bn2 = nn.BatchNorm1d(64)
        self.relu2 = nn.ReLU(inplace=True)

        self.fc3=nn.Linear(64,32)
        self.bn3 = nn.BatchNorm1d(32)
        self.relu3 = nn.ReLU(inplace=True)

        self.fc4=nn.Linear(32,num_classes)
        self.dropout_rate=dropout_rate
    def forward(self, x):
        x = self.conv1(x)
        x = self.bn(x)
        x = self.relu(x)
        x=x.view(x.size(0), -1)
        x = self.fc1(x)
        x = self.bn1(x)
        x = self.relu1(x)
        if self.dropout_rate>0:
            x = F.dropout(x, p=self.dropout_rate, training=self.training)
        x = self.fc2(x)
        x = self.bn2(x)
        x = self.relu2(x)
        if self.dropout_rate > 0:
            x = F.dropout(x, p=self.dropout_rate, training=self.training)
        x = self.fc3(x)
        x = self.bn3(x)
        x = self.relu3(x)
        if self.dropout_rate > 0:
            x = F.dropout(x, p=self.dropout_rate, training=self.training)
        x=self.fc4(x)
        return x,None,None#to easily make use of train a val a



class Phase2_CNNMLP2(nn.Module):

    def __init__(self,input_channel=8,num_classes=4,dropout_rate=0,neighbor=7):
        """
        :param num_classes:num of classes that we need to predict
        """
        super(Phase2_CNNMLP2, self).__init__()
        self.conv1 = nn.Conv3d(
            input_channel,
            32,
            kernel_size=3,
            stride=(2, 2, 2),
            padding=0,
            bias=False)  # Out input channel size is 1, our map
        self.bn = nn.BatchNorm3d(32)
        self.relu = nn.ReLU(inplace=True)
        #then fully connected layer
        self.fc1 = nn.Linear(32*(int((neighbor-3)/2+1)**3), 128)
        self.bn1 = nn.BatchNorm1d(128)
        self.relu1 = nn.ReLU(inplace=True)

        self.fc2=nn.Linear(128,64)
        self.bn2 = nn.BatchNorm1d(64)
        self.relu2 = nn.ReLU(inplace=True)

        self.fc3=nn.Linear(64,32)
        self.bn3 = nn.BatchNorm1d(32)
        self.relu3 = nn.ReLU(inplace=True)

        self.fc4=nn.Linear(32,num_classes)
        self.dropout_rate=dropout_rate
    def forward(self, x):
        x = self.conv1(x)
        x = self.bn(x)
        x = self.relu(x)
        x=x.view(x.size(0), -1)
        x = self.fc1(x)
        x = self.bn1(x)
        x = self.relu1(x)
        if self.dropout_rate>0:
            x = F.dropout(x, p=self.dropout_rate, training=self.training)
        x = self.fc2(x)
        x = self.bn2(x)
        x = self.relu2(x)
        if self.dropout_rate > 0:
            x = F.dropout(x, p=self.dropout_rate, training=self.training)
        x = self.fc3(x)
        x = self.bn3(x)
        x = self.relu3(x)
        if self.dropout_rate > 0:
            x = F.dropout(x, p=self.dropout_rate, training=self.training)
        x=self.fc4(x)
        return x,None,None#to easily make use of train a val a
