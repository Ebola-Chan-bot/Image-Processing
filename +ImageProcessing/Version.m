function V = Version
V.Me='v3.7.0';
V.MatlabExtension='v18.2.0';
V.MATLAB='R2022b';
persistent NewVersion
try
	if isempty(NewVersion)
		NewVersion=TextAnalytics.CheckUpdateFromGitHub('https://github.com/Silver-Fang/Image-Processing/releases','埃博拉酱的图像处理工具箱',V.Me);
	end
catch ME
	if ME.identifier~="MATLAB:undefinedVarOrClass"
		ME.rethrow;
	end
end