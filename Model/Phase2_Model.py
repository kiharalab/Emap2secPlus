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