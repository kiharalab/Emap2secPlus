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

import torch.nn as nn
import math
from functools import partial
from Model.resnetxt import ResNeXtBottleneck
from Model.resnetxt import downsample_basic_block
from Model.Attention import Attention3D
class Attention_ResNeXt(nn.Module):

    def __init__(self,
                 block,
                 attention_block,
                 layers,
                 sample_size,
                 sample_duration,
                 shortcut_type='B',
                 cardinality=32,
                 num_classes=400):
        """
        :param block: choose bloce
        :param layers: 4-length list: combine 4 differen layer with different specifications
        :param sample_size:
        :param sample_duration:
        :param shortcut_type:
        :param cardinality:
        :param num_classes:num of classes that we need to predict
        """
        self.inplanes = 64
        super(Attention_ResNeXt, self).__init__()
        self.attention_block1=attention_block(1,'relu')
        self.attention_block2 = attention_block(2048, 'relu')
        self.conv1 = nn.Conv3d(
            1,
            64,
            kernel_size=3,
            stride=(2, 2, 2),
            padding=1,
            bias=False)#Out input channel size is 1, our map
        self.bn1 = nn.BatchNorm3d(64)
        self.relu = nn.ReLU(inplace=True)
        self.maxpool = nn.MaxPool3d(kernel_size=(3, 3, 3), stride=2, padding=1)
        self.layer1 = self._make_layer(block, 128, layers[0], shortcut_type,
                                       cardinality)
        self.layer2 = self._make_layer(
            block, 256, layers[1], shortcut_type, cardinality, stride=2)
        self.layer3 = self._make_layer(
            block, 512, layers[2], shortcut_type, cardinality, stride=2)
        self.layer4 = self._make_layer(
            block, 1024, layers[3], shortcut_type, cardinality, stride=2)
        last_duration = int(math.ceil(sample_duration / 16))
        last_size = int(math.ceil(sample_size / 32))
        self.avgpool = nn.AvgPool3d(
            (last_duration, last_size, last_size), stride=1)
        self.fc = nn.Linear(cardinality * 32 * block.expansion, num_classes)

        for m in self.modules():
            if isinstance(m, nn.Conv3d):
                m.weight = nn.init.kaiming_normal(m.weight, mode='fan_out')#Kaiming init for weight distribution
            elif isinstance(m, nn.BatchNorm3d):
                m.weight.data.fill_(1)
                m.bias.data.zero_()

    def _make_layer(self,
                    block,
                    planes,
                    blocks,
                    shortcut_type,
                    cardinality,
                    stride=1):
        """
        :param block: residue block define before
        :param planes: number of filters
        :param blocks: number of times use this block
        :param shortcut_type: specify use downsample or not
        :param cardinality:
        :param stride: stride in CNN
        :return:
        """
        downsample = None
        if stride != 1 or self.inplanes != planes * block.expansion:
            if shortcut_type == 'A':
                downsample = partial(
                    downsample_basic_block,
                    planes=planes * block.expansion,
                    stride=stride)
            else:
                downsample = nn.Sequential(
                    nn.Conv3d(
                        self.inplanes,
                        planes * block.expansion,
                        kernel_size=1,
                        stride=stride,
                        bias=False), nn.BatchNorm3d(planes * block.expansion))

        layers = []
        layers.append(
            block(self.inplanes, planes, cardinality, stride, downsample))
        self.inplanes = planes * block.expansion
        for i in range(1, blocks):
            layers.append(block(self.inplanes, planes, cardinality))

        return nn.Sequential(*layers)

    def forward(self, x):
        x,p1=self.attention_block1(x)

        x = self.conv1(x)
        x = self.bn1(x)
        x = self.relu(x)
        x = self.maxpool(x)

        x = self.layer1(x)
        x = self.layer2(x)
        x = self.layer3(x)
        x = self.layer4(x)

        x = self.avgpool(x)
        x,p2=self.attention_block2(x)
        x = x.view(x.size(0), -1)
        x = self.fc(x)
        #x =F.softmax(x, dim=1)

        return x,p1,p2

class NoAttention_ResNeXt(nn.Module):

    def __init__(self,
                 block,
                 attention_block,
                 layers,
                 sample_size,
                 sample_duration,
                 shortcut_type='B',
                 cardinality=32,
                 num_classes=400):
        """
        :param block: choose bloce
        :param layers: 4-length list: combine 4 differen layer with different specifications
        :param sample_size:
        :param sample_duration:
        :param shortcut_type:
        :param cardinality:
        :param num_classes:num of classes that we need to predict
        """
        self.inplanes = 64
        super(NoAttention_ResNeXt, self).__init__()
        self.attention_block1=attention_block(1,'relu')
        self.attention_block2 = attention_block(2048, 'relu')
        self.conv1 = nn.Conv3d(
            1,
            64,
            kernel_size=3,
            stride=(2, 2, 2),
            padding=1,
            bias=False)#Out input channel size is 1, our map
        self.bn1 = nn.BatchNorm3d(64)
        self.relu = nn.ReLU(inplace=True)
        self.maxpool = nn.MaxPool3d(kernel_size=(3, 3, 3), stride=2, padding=1)
        self.layer1 = self._make_layer(block, 128, layers[0], shortcut_type,
                                       cardinality)
        self.layer2 = self._make_layer(
            block, 256, layers[1], shortcut_type, cardinality, stride=2)
        self.layer3 = self._make_layer(
            block, 512, layers[2], shortcut_type, cardinality, stride=2)
        self.layer4 = self._make_layer(
            block, 1024, layers[3], shortcut_type, cardinality, stride=2)
        last_duration = int(math.ceil(sample_duration / 16))
        last_size = int(math.ceil(sample_size / 32))
        self.avgpool = nn.AvgPool3d(
            (last_duration, last_size, last_size), stride=1)
        self.fc = nn.Linear(cardinality * 32 * block.expansion, num_classes)

        for m in self.modules():
            if isinstance(m, nn.Conv3d):
                m.weight = nn.init.kaiming_normal(m.weight, mode='fan_out')#Kaiming init for weight distribution
            elif isinstance(m, nn.BatchNorm3d):
                m.weight.data.fill_(1)
                m.bias.data.zero_()

    def _make_layer(self,
                    block,
                    planes,
                    blocks,
                    shortcut_type,
                    cardinality,
                    stride=1):
        """
        :param block: residue block define before
        :param planes: number of filters
        :param blocks: number of times use this block
        :param shortcut_type: specify use downsample or not
        :param cardinality:
        :param stride: stride in CNN
        :return:
        """
        downsample = None
        if stride != 1 or self.inplanes != planes * block.expansion:
            if shortcut_type == 'A':
                downsample = partial(
                    downsample_basic_block,
                    planes=planes * block.expansion,
                    stride=stride)
            else:
                downsample = nn.Sequential(
                    nn.Conv3d(
                        self.inplanes,
                        planes * block.expansion,
                        kernel_size=1,
                        stride=stride,
                        bias=False), nn.BatchNorm3d(planes * block.expansion))

        layers = []
        layers.append(
            block(self.inplanes, planes, cardinality, stride, downsample))
        self.inplanes = planes * block.expansion
        for i in range(1, blocks):
            layers.append(block(self.inplanes, planes, cardinality))

        return nn.Sequential(*layers)

    def forward(self, x):
        #x,p1=self.attention_block1(x)

        x = self.conv1(x)
        x = self.bn1(x)
        x = self.relu(x)
        x = self.maxpool(x)

        x = self.layer1(x)
        x = self.layer2(x)
        x = self.layer3(x)
        x = self.layer4(x)

        x = self.avgpool(x)
        #x,p2=self.attention_block2(x)
        x = x.view(x.size(0), -1)
        x = self.fc(x)
        #x =F.softmax(x, dim=1)

        return x,0,0


class ResNeXt_Combine(nn.Module):

    def __init__(self,
                 block,
                 layers,
                 sample_size,
                 sample_duration,
                 shortcut_type='B',
                 cardinality=32,
                 num_classes=400,input_dim=1):
        """
        :param block: choose bloce
        :param layers: 4-length list: combine 4 differen layer with different specifications
        :param sample_size:
        :param sample_duration:
        :param shortcut_type:
        :param cardinality:
        :param num_classes:num of classes that we need to predict
        """
        self.inplanes = 64
        super(ResNeXt_Combine, self).__init__()
        self.conv1 = nn.Conv3d(
            input_dim,
            64,
            kernel_size=3,
            stride=(2, 2, 2),
            padding=1,
            bias=False)#Out input channel size is 1, our map
        self.bn1 = nn.BatchNorm3d(64)
        self.relu = nn.ReLU(inplace=True)
        self.maxpool = nn.MaxPool3d(kernel_size=(3, 3, 3), stride=2, padding=1)
        self.layer1 = self._make_layer(block, 128, layers[0], shortcut_type,
                                       cardinality)
        self.layer2 = self._make_layer(
            block, 256, layers[1], shortcut_type, cardinality, stride=2)
        self.layer3 = self._make_layer(
            block, 512, layers[2], shortcut_type, cardinality, stride=2)
        self.layer4 = self._make_layer(
            block, 1024, layers[3], shortcut_type, cardinality, stride=2)
        last_duration = int(math.ceil(sample_duration / 16))
        last_size = int(math.ceil(sample_size / 32))
        self.avgpool = nn.AvgPool3d(
            (last_duration, last_size, last_size), stride=1)
        self.fc = nn.Linear(cardinality * 32 * block.expansion, num_classes)

        for m in self.modules():
            if isinstance(m, nn.Conv3d):
                m.weight = nn.init.kaiming_normal(m.weight, mode='fan_out')#Kaiming init for weight distribution
            elif isinstance(m, nn.BatchNorm3d):
                m.weight.data.fill_(1)
                m.bias.data.zero_()

    def _make_layer(self,
                    block,
                    planes,
                    blocks,
                    shortcut_type,
                    cardinality,
                    stride=1):
        """
        :param block: residue block define before
        :param planes: number of filters
        :param blocks: number of times use this block
        :param shortcut_type: specify use downsample or not
        :param cardinality:
        :param stride: stride in CNN
        :return:
        """
        downsample = None
        if stride != 1 or self.inplanes != planes * block.expansion:
            if shortcut_type == 'A':
                downsample = partial(
                    downsample_basic_block,
                    planes=planes * block.expansion,
                    stride=stride)
            else:
                downsample = nn.Sequential(
                    nn.Conv3d(
                        self.inplanes,
                        planes * block.expansion,
                        kernel_size=1,
                        stride=stride,
                        bias=False), nn.BatchNorm3d(planes * block.expansion))

        layers = []
        layers.append(
            block(self.inplanes, planes, cardinality, stride, downsample))
        self.inplanes = planes * block.expansion
        for i in range(1, blocks):
            layers.append(block(self.inplanes, planes, cardinality))

        return nn.Sequential(*layers)

    def forward(self, x):
        x = self.conv1(x)
        x = self.bn1(x)
        x = self.relu(x)
        x = self.maxpool(x)
        x = self.layer1(x)
        x = self.layer2(x)
        x = self.layer3(x)
        x = self.layer4(x)
        x = self.avgpool(x)
        x = x.view(x.size(0), -1)
        x = self.fc(x)

        return x,0,0


def Aresnet50(**kwargs):
    """Constructs a ResNet-50 Model. We will use this
    """
    model = Attention_ResNeXt(ResNeXtBottleneck,Attention3D, [3, 4, 6, 3], **kwargs)
    return model

def Resnet50(**kwargs):
    """Constructs a ResNet-50 Model. We will use this
    """
    model = NoAttention_ResNeXt(ResNeXtBottleneck,Attention3D, [3, 4, 6, 3], **kwargs)
    return model

def Aresnet20(**kwargs):
    """Constructs a ResNet-20 Model. We will use this，
    actually 6*3+2=20
    """
    model = Attention_ResNeXt(ResNeXtBottleneck,Attention3D,  [1, 2, 2, 1], **kwargs)
    return model


def Aresnet101(**kwargs):
    """Constructs a ResNet-101 Model.
    """
    model = Attention_ResNeXt(ResNeXtBottleneck,Attention3D,  [3, 4, 23, 3], **kwargs)
    return model


def Aresnet152(**kwargs):
    """Constructs a ResNet-101 Model.
    """
    model = Attention_ResNeXt(ResNeXtBottleneck,Attention3D,  [3, 8, 36, 3], **kwargs)
    return model

def Resnet_Phase1_Combine(**kwargs):
    model=ResNeXt_Combine(ResNeXtBottleneck,[1, 2, 2, 1], **kwargs)
    return model
