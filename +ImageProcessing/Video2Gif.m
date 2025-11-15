%[text] 将视频转换为GIF动图
%[text] ## 语法
%[text] ```matlabCodeExample
%[text] ImageProcessing.Video2Gif(VideoPath,GifPath,Name=Value);
%[text] ```
%[text] ## 输入参数
%[text] VideoPath(1,1)string，输入视频路径
%[text] GifPath(1,1)string，输出GIF路径
%[text] ### 名称值参数
%[text] Scale(1,1)缩放倍率，小于1则缩小，大于1则放大。视频的宽高和帧数都会同步缩放，因此生成文件尺寸会变成缩放倍率的立方倍。
%[text] CropRange(1,3)cell，裁剪范围。前两个元胞里是(1,:)，分别表示要截取的高、宽方向像素索引；第3个元胞里是(1,2)，表示要截取的首帧和尾帧。此外，三个元胞里也都可以设为':'，表示不作截取。截取范围以Scale缩放之前的图像尺寸为准。
%[text] **See also** [imresize3](<matlab:doc imresize3>) [rgb2ind](<matlab:doc rgb2ind>) [VideoReader](<matlab:doc VideoReader>) [imwrite](<matlab:doc imwrite>) [ImageProcessing.Ind2Rgb](<matlab:doc ImageProcessing.Ind2Rgb>) [ImageProcessing.Gif2Video](<matlab:doc ImageProcessing.Gif2Video>)
function Video2Gif(VideoPath,GifPath,options)
arguments
	VideoPath
	GifPath
	options.Scale=1
	options.CropRange={':',':',':'}
end
if options.CropRange{3}==':'
	Range={};
else
	Range=options.CropRange(3);
end
VR=VideoReader(VideoPath);
Frames=VR.read(Range{:},'native');
Frames=imresize3(permute(Frames(options.CropRange{1:2},:,:),[1,2,4,3]),options.Scale);
[ImageHeight,ImageWidth,NumFrames]=size(Frames,1,2,3);
[Frames,ColorMap]=rgb2ind(reshape(Frames,ImageHeight,[],3),256);
imwrite(reshape(Frames,ImageHeight,ImageWidth,1,NumFrames),ColorMap,GifPath,'gif',DelayTime=1/(options.Scale*VR.FrameRate),LoopCount=Inf);
%经测试，跟一张一张写并手动做重复像素优化大小无差异
end

%[appendix]{"version":"1.0"}
%---
