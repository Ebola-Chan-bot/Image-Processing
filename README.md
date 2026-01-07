埃博拉酱的图像处理工具包，提供一系列 MATLAB Image Processing Toolbox 内置函数所欠缺，但却常用的图像处理功能。

本包中所有函数均在ImageProcessing命名空间下，使用前需import。
```MATLAB
import ImageProcessing.*
```
函数文件内有详细文档，可用doc命令查看。
```MATLAB
%内置dct2的升级版，支持多维数组和GPU数组
function arg1=Dct2(varargin)

%内置edge的升级版，支持任意维度图像数组（但只在前两维计算边缘），支持GPU数组
function BW=Edge(I)

%内置fitgeotform2d的升级版，根据输入点数自动选择合适的变换方案
function TForm2D = FitGeoTForm2D(Points)

%根据变换前后的点XY坐标计算出变换矩阵
function LeftMatrix = FixedPointTform2D(FromPoints,ToPoints)

%对隐函数在指定范围内作图
function [Line,Scatter]=FPlot(Func,MinX,MaxX,MinY,MaxY,Step,options)

%将GIF多帧图像读入为RGB视频
function [Video,FrameRate] = Gif2Video(GifPath)

%内置idct2的升级版，支持多维数组和GPU数组
function a = IDct2(varargin)

%内置imwarp的升级版，支持3D图像批量变换
function varargout=ImWarp(Image,varargin)

%内置ind2rgb的升级版，支持任意维度图像（但第3维必须单一，留给RGB通道）
function RGB = Ind2Rgb(Index,Map)

%将两张图叠成一张PNG幻影坦克
function [Color,Alpha] = MirageTank(BlackImage,WhiteImage,varargin)

%内置normxcorr2的升级版，支持多图批量操作
function C=NormXCorr2(template,A,Partial)

%对NormXCorr2的模板进行预处理
function FDN = Nxc2TPreprocess(Template,SizeA)

%使用Nxc2TPreprocess的预处理结果计算NormXCorr2
function C = Nxc2APostprocess(A,SizeT,FDN,Partial)

%将视频转换为GIF动图
function Video2Gif(VideoPath,GifPath,options)
```