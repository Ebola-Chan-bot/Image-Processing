埃博拉酱的图像处理工具包，提供一系列MATLAB Image Processing Toolbox内置函数所欠缺，但却常用的图像处理功能。依赖[埃博拉酱的MATLAB扩展](https://ww2.mathworks.cn/matlabcentral/fileexchange/96344-matlab-extension)

本项目的发布版本号遵循[语义化版本](https://semver.org/lang/zh-CN/)规范。开发者认为这是一个优秀的规范，并向每一位开发者推荐遵守此规范。
# 目录
本包中所有函数均在ImageProcessing命名空间下，使用前需import。使用命名空间是一个好习惯，可以有效防止命名冲突，避免编码时不必要的代码提示干扰。
- [VideoBatchMeasure](#VideoBatchMeasure) 根据ROI批量测量视频像素值
- [VideoDrawROI](#VideoDraoROI) 给视频画ROI用于测量
# VideoBatchMeasure
根据ROI批量测量视频像素值
## 输入参数
Flags(1,1)string，重复，可选设置以下功能旗帜：
- Gpu，使用GPU加速计算
- Verbose，输出进度信息
- Parallel，使用并行计算。如果还指定了Gpu，将只在一个进程中使用GPU，其他进程仍使用CPU。

VideoPaths(:,1)string，名称值，视频路径，默认打开文件选择对话框要求用户手动选择

Masks(:,1)cell，名称值，测量遮罩。默认调用ImageProcessing.VideoDraoROI要求用户手动绘制。
## 返回值
Measurements(:,1)cell，每个元胞对应一个视频，元胞内是(:,:)single，第1维ROI，第2维时间帧，每个数值代表该ROI内所有像素在该时间帧内的平均值
Masks(:,1)cell，同输入参数中的Masks。如果你的遮罩是临时手绘的，可以用这个返回值取得遮罩。
# VideoDrawROI
给视频画ROI用于测量
## 输入参数
NoROIs(1,1)uint8，必需，每个视频要画多少个ROI

VideoPaths(:,1)string，可选，所有视频文件路径。默认打开文件选择对话框要求用户手动选择
## 返回值
Masks(:,1)cell，每个视频的测量遮罩，每个元胞对应一个视频。元胞内是(:,:,:)logical，前两维是视频的高宽，第3维是不同的ROI。