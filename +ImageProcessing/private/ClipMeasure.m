function Measurement = ClipMeasure(VR,GpuAvailable,VideoMasks,varargin)
Data=VR.read(varargin{:});
if GpuAvailable
	Data=gpuArray(Data);
end
Data=single(Data);
if size(Data,3)>1
	Data=0.2989*Data(:,:,1,:)+0.5870*Data(:,:,2,:)+0.1140*Data(:,:,3,:);
end
Data=reshape(Data,[],size(Data,4));
Measurement=gather(VideoMasks*Data);