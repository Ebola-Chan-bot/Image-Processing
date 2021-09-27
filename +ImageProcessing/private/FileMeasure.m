function Measurements = FileMeasure(Path,VideoMasks,FreeMemory,GpuAvailable)
VR=VideoReader(Path);
Data=VR.readFrame;
VideoMasks=reshape(VideoMasks,[],size(VideoMasks,3))';
VideoMasks=VideoMasks./sum(VideoMasks,2);
FrameSize=numel(Data)*4;
NoFrames=VR.NumFrames;
if FrameSize*NoFrames<FreeMemory
	Measurements=ClipMeasure(VR,GpuAvailable,VideoMasks);
else
	ActualBlocks=ceil(NoFrames/floor(FreeMemory/FrameSize));
	Frames1=floor(NoFrames/ActualBlocks);
	Frames2=Frames1+1;
	Blocks1=ActualBlocks*Frames2-NoFrames;
	VideoMeasurement=cell(1,ActualBlocks);
	ReadRange=[1 Frames1];
	for B=1:Blocks1
		VideoMeasurement{B}=ClipMeasure(VR,GpuAvailable,VideoMasks,ReadRange);
		ReadRange=ReadRange+Frames1;
	end
	ReadRange(2)=ReadRange(2)+1;
	for B=Blocks1+1:ActualBlocks
		VideoMeasurement{B}=ClipMeasure(VR,GpuAvailable,VideoMasks,ReadRange);
		ReadRange=ReadRange+Frames2;
	end
	Measurements=[VideoMeasurement{:}];
end