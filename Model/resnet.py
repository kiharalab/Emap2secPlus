import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.autograd import Variable
import math
from functools import partial
def conv3x3x3(in_planes, out_planes, stride=1):
    # 3x3x3 convolution with padding
    return nn.Conv3d(
        in_planes,
        out_planes,
        kernel_size=3,
        stride=stride,
        padding=1,
        bias=False)


class ResidualBlock(nn.Module):
    def __init__(self, in_channels=None, out_channels=32):
        super(ResidualBlock, self).__init__()
        if in_channels is None:
            in_channels = out_channels
        self.conv1 = nn.Conv3d(in_channels, out_channels, kernel_size=3, padding=1)
        self.bn1 = nn.BatchNorm3d(out_channels)
        self.prelu = nn.PReLU()
        self.conv2 = nn.Conv3d(out_channels, out_channels, kernel_size=3, padding=1)
        self.bn2 = nn.BatchNorm3d(out_channels)

    def forward(self, x):
        #         print("Res input:"+str(x.size()))
        residual = self.conv1(x)
        residual = self.bn1(residual)
        residual = self.prelu(residual)
        residual = self.conv2(residual)
        residual = self.bn2(residual)
        #         print("Res out:"+str(residual.size()))
        return x + residual

def normal_init(m, mean, std):
    if isinstance(m, nn.ConvTranspose2d) or isinstance(m, nn.Conv2d):
        m.weight.data.normal_(mean, std)
        m.bias.data.zero_()
class Discriminator(nn.Module):
    # initializers
    def __init__(self, d, no_of_blocks):
        super(Discriminator, self).__init__()

        model_sequence = []
        model_sequence += [ResidualBlock(1, d)]
        self.no_of_blocks = no_of_blocks
        for i in range(self.no_of_blocks - 1):
            model_sequence += [ResidualBlock(d)]
        # model_sequence += [nn.Conv3d(d, 1, kernel_size=3, stride=1,padding=1)]
        self.model = nn.Sequential(*model_sequence)

    def weight_init(self, mean=0, std=0.01):
        for m in self._modules:
            normal_init(self._modules[m], mean, std)

            # forward method

    def forward(self, input):
        #         print("D imnput:"+str(input.size()))
        x = self.model(input)
        #         print("D output:"+str(x.size()))
        x = torch.sigmoid(x)
        return x