%[text] 内置fitgeotform2d的升级版，根据输入点数自动选择合适的变换方案
%[text] 此函数调用fitgeotform2d，但用户无需关心更多复杂的变换方案，只需提供输入输出点对，自动选择变换方案。如果结果不满意，可以再尝试内置fitgeotform2d。
%[text] 具体来说，变换方案选择策略是：
%[text] 1. 1对点，平移变换
%[text] 2. 2对点，相似变换，额外考虑旋转和各向同性缩放，
%[text] 3. 3对点，仿射变换，额外考虑各向异性缩放
%[text] 4. 4~5对点，影射变换，额外考虑投影
%[text] 5. 6对及以上，多项式变换，允许非线性局部畸变 \
%[text] ## 语法
%[text] ```matlabCodeExample
%[text] TForm2D = ImageProcessing.FitGeoTForm2D(Points);
%[text] ```
%[text] ## 输入参数
%[text] Points(:,2,2)，输入点对张量。第1维是不同的点对，第2维是XY，第3维是变换前后。输入点对越多，变换越精确。
%[text] ## 返回值
%[text] TForm2D(1,1)images.geotrans.internal.GeometricTransformation，可用于imwarp
%[text] **See also** [fitgeotform2d](<matlab:doc fitgeotform2d>) [imwarp](<matlab:doc imwarp>)
function TForm2D = FitGeoTForm2D(Points)
NumPoints=height(Points);
if NumPoints<2
	TForm2D=simtform2d([1,0,0;0,1,0;Points(:,:,2)-Points(:,:,1),1]);
elseif NumPoints<3
	TForm2D=fitgeotform2d(Points(:,:,1),Points(:,:,2),'similarity');
elseif NumPoints<4
	TForm2D=fitgeotform2d(Points(:,:,1),Points(:,:,2),'affine');
elseif NumPoints<6
	TForm2D=fitgeotform2d(Points(:,:,1),Points(:,:,2),'projective');
elseif NumPoints<10
	TForm2D=fitgeotform2d(Points(:,:,1),Points(:,:,2),'polynomial',2);
elseif NumPoints<15
	TForm2D=fitgeotform2d(Points(:,:,1),Points(:,:,2),'polynomial',3);
else
	TForm2D=fitgeotform2d(Points(:,:,1),Points(:,:,2),'polynomial',4);
end
end

%[appendix]{"version":"1.0"}
%---
