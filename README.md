埃博拉酱的图像处理工具包，提供一系列MATLAB Image Processing Toolbox内置函数所欠缺，但却常用的图像处理功能

本项目的发布版本号遵循[语义化版本](https://semver.org/lang/zh-CN/)规范。开发者认为这是一个优秀的规范，并向每一位开发者推荐遵守此规范。
# 目录
本包中所有函数均在ImageProcessing命名空间下，使用前需import。使用命名空间是一个好习惯，可以有效防止命名冲突，避免编码时不必要的代码提示干扰。
- [VideoPixelMeasure](#VideoPixelMeasure)
# VideoPixelMeasure
视频像素测量

本函数用ROI遮罩测量视频平面上像素值随视频帧的时间变化，彩色计算灰度。
## 名称值参数
VideoPath(1,1)string，视频文件路径。默认打开文件选择对话框供用户手动选择。

Masks(:,:,:)logical，测量遮罩。第1维高度，第2维宽度，第3维不同ROI。默认要求用户手动画圈。
## 返回值
Measurements(:,:)double，测量值。第1维ROI，第2维时间。每个数值是某时刻某ROI的全像素平均值。