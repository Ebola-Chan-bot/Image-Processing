function [outputImage,outputRef] = CImWarp(varargin)

narginchk(2,inf);

isDisplacementFieldSyntax = isnumeric(varargin{2});
if isDisplacementFieldSyntax
    % Handle displacement field syntaxes as a completely separate case.
    [parsedInputs,catConverter,isInputCategorical] = parseInputsDisplacementFieldSyntax(varargin{:});
    method = parsedInputs.InterpolationMethod;
    fillValues = parsedInputs.FillValues;
    D = parsedInputs.DisplacementField;
    
    outputImage = images.geotrans.internal.applyDisplacementField(parsedInputs.InputImage,...
        D,method,fillValues, parsedInputs.SmoothEdges);
    
    sizeD = size(D);
    if ndims(D) == 3
        outputRef = imref2d(sizeD(1:2));
    else
        outputRef = imref3d(sizeD(1:3));
    end
    
    if isInputCategorical
        outputImage = catConverter.numeric2Categorical(outputImage);
    end

    return;
end

[R_A, varargin] = images.geotrans.internal.preparseSpatialReferencingObjects(varargin{:});
[parsedInputs,catConverter,isInputCategorical] = parseInputs(varargin{:});

method = parsedInputs.InterpolationMethod;
fillValues = parsedInputs.FillValues;
SmoothEdges = parsedInputs.SmoothEdges;

% Some older geometric transformation types are automatically converted by
% parseInputs to the corresponding types introduced in R2022b:
% 
% rigid2d -> rigidtform2d
% rigid3d -> rigidtform3d
% affine2d -> affinetform2d
% affine3d -> affinetform3d
% projective2d -> projtform2d
%
% Subsequent code can assume the tform variable below does not contain one
% of the older types.
tform = parsedInputs.GeometricTransform;

% Convert translation, similarity, and rigid transformations to affine so
% that code paths optimized for affine transformations can be used later.
if isa(tform,'transltform2d') || isa(tform,'rigidtform2d') || isa(tform,'simtform2d')
    tform = affinetform2d(tform);
elseif isa(tform,'transltform3d') || isa(tform,'rigidtform3d') || isa(tform,'simtform3d')
    tform = affinetform3d(tform);
end

fillValues = cast(fillValues, 'like', parsedInputs.InputImage);

% Check agreement of input image with dimensionality of tform
% images.geotrans.internal.checkImageAgreementWithTform(parsedInputs.InputImage,tform);

inputSpatialReferencingNotSpecified = isempty(R_A);
if inputSpatialReferencingNotSpecified
    if isa(R_A,'imref3d')
        R_A = imref3d(size(parsedInputs.InputImage,1:3));
    else
        R_A = imref2d(size(parsedInputs.InputImage));
    end
else
    % Check agreement of input spatial referencing object with input image.
    images.spatialref.internal.checkSpatialRefAgreementWithInputImage(parsedInputs.InputImage,R_A);
end

% check agreement of fillValues with dimensionality of problem
images.internal.checkFillValues(fillValues,parsedInputs.InputImage,tform.Dimensionality);

% If the 'OutputView' was not specified, we have to determine the world
% limits and the image size of the output from the input spatial
% referencing information and the geometric transform.
if isempty(parsedInputs.OutputView)
    outputRef = calculateOutputSpatialReferencing(R_A,tform);
else
    outputRef = parsedInputs.OutputView;
    checkOutputViewAgreementWithTform(outputRef,tform);
end

outputImage = remapPointsAndResample(parsedInputs.InputImage,R_A,tform,outputRef,method,fillValues, SmoothEdges);

if isInputCategorical
    outputImage = catConverter.numeric2Categorical(outputImage);
else
    outputImage = cast(outputImage,'like',parsedInputs.InputImage);
end

function outputImage = remapPointsAndResample(inputImage,R_A,tform,outputRef,method,fillValues, SmoothEdges)

if tform.Dimensionality ==2
    
    useIPPLibrary = images.internal.useIPPLibrary();
    useIPPAffine = useIPPLibrary && isa(tform,'affinetform2d') && ~isProblemSizeTooBig(inputImage) && ~isvector(inputImage);

    if useIPPAffine
        outputImage = ippWarpAffine(inputImage,R_A,tform,outputRef,method,fillValues, SmoothEdges);
    elseif isa(tform,'affinetform2d') || isa(tform,'projtform2d')
        outputImage = remapAndResampleInvertible2d(inputImage,R_A,tform,outputRef,method,fillValues, SmoothEdges);
    else
        outputImage = remapAndResampleGeneric2d(inputImage,R_A,tform,outputRef,method,fillValues, SmoothEdges);
    end

else
    %3d transformation
    if isa(tform,'affinetform3d')
        outputImage = remapAndResampleInvertible3d(inputImage,R_A,tform,outputRef,method,fillValues, SmoothEdges);
    else
        outputImage = remapAndResampleGeneric3d(inputImage,R_A,tform,outputRef,method,fillValues, SmoothEdges);
    end
end

function checkOutputViewAgreementWithTform(Rout,tform)

if (tform.Dimensionality == 3) && ~isa(Rout,'imref3d') || ((tform.Dimensionality==2) && isa(Rout,'imref3d'))
    error(message('images:imwarp:outputViewTformDimsMismatch','''OutputView'''));
end


function R_out = calculateOutputSpatialReferencing(R_A,tform)
% Applies geometric transform to input spatially referenced grid to figure
% out the resolution and world limits after application of the forward
% transformation.
R_out = images.spatialref.internal.applyGeometricTransformToSpatialRef(R_A,tform);

function [parsedOutputs,catConverter,isInputCategorical] = parseInputsDisplacementFieldSyntax(varargin)

isInputCategorical = iscategorical(varargin{1});
catConverter = [];

parser = inputParser();
parser.addRequired('InputImage',@validateInputImage);
parser.addRequired('DisplacementField',@validateDField);
parser.addOptional('InterpolationMethod','',@validateInterpMethod);
parser.addParameter('SmoothEdges',false,@validateSmoothEdges);

if isInputCategorical
    catConverter = images.internal.utils.CategoricalConverter(categories(varargin{1}));
    parser.addParameter('FillValues',missing,@(fillVal)validateCategoricalFillValues(fillVal,catConverter.Categories));
else
    parser.addParameter('FillValues',0,@validateFillValues);
end

varargin = remapPartialParamNamesImwarp(varargin{:});

parser.parse(varargin{:});

parsedOutputs = parser.Results;

displacementFieldDataMismatch = ndims(parsedOutputs.DisplacementField) == 4 &&...
    ismatrix(parsedOutputs.InputImage);
if displacementFieldDataMismatch
    error(message('images:imwarp:displacementField3dData2d'));
end

if isInputCategorical
    % For categorical inputs,
    % 1. 'nearest' is the only interpolation method.
    % 2.  Default Fill values will be set to missing, which corresponds to
    % '<undefined>' label in the output categorical result. Any other
    % FillValue results in an error.
   if ~(isempty(parsedOutputs.InterpolationMethod)||...
                       isequal(parsedOutputs.InterpolationMethod,'nearest'))
       error(message('MATLAB:images:validate:badMethodForCategorical'));
   end
   
   % Get corresponding numeric value for the classname
   parsedOutputs.FillValues = catConverter.getNumericValue(parsedOutputs.FillValues);
   
   % Set default interpolation to 'nearest' for categorical inputs
   parsedOutputs.InterpolationMethod = 'nearest';
   
   parsedOutputs.InputImage = catConverter.categorical2Numeric(parsedOutputs.InputImage);
   
else
    % Set default interpolation to 'linear'
    if(isempty(parsedOutputs.InterpolationMethod))
        parsedOutputs.InterpolationMethod = 'linear';
    end
end

method = postProcessMethodString(parsedOutputs.InterpolationMethod);
parsedOutputs.InterpolationMethod = method;

function [parsedOutputs,catConverter,isInputCategorical] = parseInputs(varargin)

isInputCategorical = iscategorical(varargin{1});
catConverter = [];

parser = inputParser();
parser.addRequired('InputImage',@validateInputImage);
parser.addRequired('GeometricTransform',@validateTform);
parser.addOptional('InterpolationMethod','',@validateInterpMethod);
parser.addParameter('SmoothEdges',false,@validateSmoothEdges);
parser.addParameter('OutputView',[],@(ref) isa(ref,'imref2d') || isa(ref,'imref3d'));

if isInputCategorical
    parser.addParameter('FillValues',missing,@(fillVal)validateCategoricalFillValues(fillVal,categories(varargin{1})));
else
    parser.addParameter('FillValues',0,@validateFillValues);
end

varargin = remapPartialParamNamesImwarp(varargin{:});

parser.parse(varargin{:});

parsedOutputs = parser.Results;

% Convert certain older geometric transformation types to the types
% introduced in R2022b.
tform_in = parsedOutputs.GeometricTransform;
if isa(tform_in,'rigid2d')
    tform_converted = rigidtform2d(tform_in.T');
elseif isa(tform_in,'affine2d')
    tform_converted = affinetform2d(tform_in.T');
elseif isa(tform_in,'projective2d')
    tform_converted = projtform2d(tform_in.T');
elseif isa(tform_in,'rigid3d')
    tform_converted = rigidtform3d(tform_in.T');
elseif isa(tform_in,'affine3d')
    tform_converted = affinetform3d(tform_in.T');
else
    tform_converted = tform_in;
end
parsedOutputs.GeometricTransform = tform_converted;

if isInputCategorical
    % For categorical inputs,
    % 1. 'nearest' is the only interpolation method.
    % 2.  Default Fill values will be set to missing, which corresponds to
    % '<undefined>' label in the output categorical result. Any other
    % FillValue results in an error.
   if ~(isempty(parsedOutputs.InterpolationMethod)||...
                       isequal(parsedOutputs.InterpolationMethod,'nearest'))
       error(message('MATLAB:images:validate:badMethodForCategorical'));
   end
   
   catConverter = images.internal.utils.CategoricalConverter(categories(parsedOutputs.InputImage));
   
   % Get corresponding numeric value for the classname
   parsedOutputs.FillValues = catConverter.getNumericValue(parsedOutputs.FillValues);
   
   % Set default interpolation to 'nearest' for categorical inputs
   parsedOutputs.InterpolationMethod = 'nearest'; 
   
   parsedOutputs.InputImage = catConverter.categorical2Numeric(parsedOutputs.InputImage);
   
else
    % Set default interpolation to 'linear'
    if(isempty(parsedOutputs.InterpolationMethod))
        parsedOutputs.InterpolationMethod = 'linear';
    end
end

method = postProcessMethodString(parsedOutputs.InterpolationMethod);
parsedOutputs.InterpolationMethod = method;

function TF = validateInterpMethod(method)

validatestring(method,...
    {'nearest','linear','cubic','bilinear','bicubic'}, ...
    'imwarp', 'InterpolationMethod');

TF = true;

function TF = validateInputImage(img)

allowedTypes = {'logical','uint8', 'uint16', 'uint32', 'int8','int16','int32','single','double','categorical'};
validateattributes(img,allowedTypes,...
    {'nonempty','nonsparse'},'imwarp','A',1);

TF = true;

function TF = validateFillValues(fillVal)

validateattributes(fillVal,{'numeric'},...
        {'nonempty','nonsparse'},'imwarp','FillValues');

TF = true;

function TF = validateCategoricalFillValues(fillVal,cats)

% FillVal for categorical input can be of the valid category in form of
% char vector, or missing
if ischar(fillVal) && any(strcmp(fillVal,cats))
    TF = true;
elseif ~ischar(fillVal) && ~isnumeric(fillVal) && isscalar(fillVal) && ismissing(fillVal)
    TF = true;
else
    error(message('MATLAB:images:validate:badFillValueForCategorical'));
end


function TF = validateSmoothEdges(SmoothEdges)
validateattributes(SmoothEdges,{'logical'},...
    {'nonempty','scalar'},'imwarp','SmoothEdges');

TF = true;


function TF = validateTform(t)

validateattributes(t,{'images.geotrans.internal.GeometricTransformation'},{'scalar','nonempty'},'imwarp','tform');

TF = true;

function methodOut = postProcessMethodString(methodIn)

methodIn = validatestring(methodIn,...
    {'nearest','linear','cubic','bilinear','bicubic'});
% We commonly use bilinear and bicubic in IPT, so both names should work
% for 2-D and 3-D input. This is consistent with interp2 and interp3 in
% MATLAB.

keys   = {'nearest','linear','cubic','bilinear','bicubic'};
values = {'nearest', 'linear','cubic','linear','cubic'};
methodMap = containers.Map(keys,values);
methodOut = methodMap(methodIn);

function TF = validateDField(D)

allowedTypes = {'logical','uint8', 'uint16', 'uint32', 'int8','int16','int32','single','double'};
validateattributes(D, allowedTypes, {'nonempty','nonsparse','finite'},'imwarp','D');

sizeD = size(D);
valid2DCase = (ndims(D) == 3) && (sizeD(3) == 2);
valid3DCase = (ndims(D) == 4) && (sizeD(4) == 3);
if ~(valid2DCase || valid3DCase)
    error(message('images:imwarp:invalidDSize'));
end

TF = true;

function varargin_out = remapPartialParamNamesImwarp(varargin)

varargin_out = varargin;
if (nargin > 2)
    % Parse input, replacing partial name matches with the canonical form.
    varargin_out(3:end) = images.internal.remapPartialParamNames({'OutputView','FillValues','SmoothEdges'}, ...
        varargin{3:end});
end


function [paddedImage,Rpadded] = padImage(A,RA,fillVal)

pad = 2;

if isscalar(fillVal)
    paddedImage = padarray(A,[pad pad],fillVal);
else
    sizeInputImage = size(A);
    sizeOutputImage = sizeInputImage;
    sizeOutputImage(1:2) = sizeOutputImage(1:2) + [2*pad 2*pad];
    if islogical(A)
        paddedImage = false(sizeOutputImage);
    else
        paddedImage = zeros(sizeOutputImage,'like',A);
    end
    [~,~,numPlanes] = size(A);
    for i = 1:numPlanes
        paddedImage(:,:,i) = padarray(A(:,:,i),[pad pad],fillVal(i));
    end

end

Rpadded = imref2d(size(paddedImage), RA.PixelExtentInWorldX*[-pad pad]+RA.XWorldLimits,...
    RA.PixelExtentInWorldY*[-pad pad]+RA.YWorldLimits);


function outputImage = remapAndResampleInvertible3d(inputImage,Rin,tform,Rout,method,fillValues, SmoothEdges)

% Define affine transformation that maps from intrinsic system of
% output image to world system of output image.
tIntrinsictoWorldOutput = Rout.TransformIntrinsicToWorld;

% Define affine transformation that maps from world system of
% input image to intrinsic system of input image.
tWorldToIntrinsicInput = Rin.TransformWorldToIntrinsic;

% Form transformation to go from output intrinsic to input intrinsic.
% NOTE: tComp is in post-multiply form.
tComp = tIntrinsictoWorldOutput / tform.A' * tWorldToIntrinsicInput;
% Find the transform that takes from input intrinsic to output intrinsic
tComp(:,4)=[0 0 0 1]; % avoid round off issues due to inversion above

if ~SmoothEdges && (string(method)~="cubic") &&...
        (isa(inputImage,'uint8') || isa(inputImage,'int16') || isa(inputImage,'uint16')|| isfloat(inputImage))&&ndims(inputImage)<=3
             
    % Fast common case
    outputImage = images.internal.builtins.warp3d(inputImage, double(tComp),  ...
                        Rout.ImageSize, fillValues, method);
else
    tformComposite = affinetform3d(tComp');
    % Form plaid grid of intrinsic points in output image.
    [dstXIntrinsic,dstYIntrinsic,dstZIntrinsic] = meshgrid(1:Rout.ImageSize(2),...
        1:Rout.ImageSize(1),...
        1:Rout.ImageSize(3));
    [srcXIntrinsic,srcYIntrinsic, srcZIntrinsic] = ...
        tformComposite.transformPointsForward(dstXIntrinsic,dstYIntrinsic, dstZIntrinsic);
    clear dstXIntrinsic dstYIntrinsic dstZIntrinsic;
    outputImage = CInterp3D(inputImage,srcXIntrinsic,srcYIntrinsic,srcZIntrinsic,method,fillValues, SmoothEdges);
end

function outputImage = remapAndResampleInvertible2d(inputImage,Rin,tform,Rout,method,fillValues, SmoothEdges)
[srcXIntrinsic,srcYIntrinsic] = images.geotrans.internal.getSourceMappingInvertible2d(Rin,tform,Rout);
% Mimics syntax of interp2. Has different edge behavior that uses 'fill'
outputImage = images.internal.interp2d(inputImage,srcXIntrinsic,srcYIntrinsic,method,fillValues, SmoothEdges);

function outputImage = remapAndResampleGeneric3d(inputImage,R_A,tform,outputRef,method,fillValues, SmoothEdges)

% Form plaid grid of intrinsic points in output image.
[X,Y,Z] = meshgrid(1:outputRef.ImageSize(2),...
    1:outputRef.ImageSize(1),...
    1:outputRef.ImageSize(3));

% Find location of pixel centers of destination image in world coordinates
% as the starting point for reverse mapping.
[X, Y, Z] = outputRef.intrinsicToWorldAlgo(X,Y,Z);

% Reverse map pixel centers from destination image to source image via
% inverse transformation.
[X,Y,Z] = tform.transformPointsInverse(X,Y,Z);

% Find srcX, srcY, srcZ in intrinsic coordinates to use when
% interpolating.
[X,Y,Z] = R_A.worldToIntrinsicAlgo(X,Y,Z);

% Mimics syntax of interp3. Has different edge behavior that uses
% 'fill'
outputImage = CInterp3D(inputImage,X,Y,Z,method,fillValues, SmoothEdges);

function outputImage = remapAndResampleGeneric2d(inputImage,R_A,tform,outputRef,method,fillValues, SmoothEdges)

% Form plaid grid of intrinsic points in output image.
[X,Y] = meshgrid(1:outputRef.ImageSize(2),1:outputRef.ImageSize(1));

% Find location of pixel centers of destination image in world coordinates
% as the starting point for reverse mapping.
[X, Y] = outputRef.intrinsicToWorldAlgo(X,Y);

% Reverse map pixel centers from destination image to source image via
% inverse transformation.
[X,Y] = tform.transformPointsInverse(X,Y);

% Find srcX srcY in intrinsic coordinates to use when interpolating.
% remapmex only knows how to work in intrinsic coordinates, interp2
% supports intrinsic or world.
[X,Y] = R_A.worldToIntrinsicAlgo(X,Y);

% Mimics syntax of interp2. Has different edge behavior that uses 'fill'
outputImage = images.internal.interp2d(inputImage,X,Y,method,fillValues, SmoothEdges);

function outputImage = ippWarpAffine(inputImage,Rin,tform,Rout,interp,fillVal, SmoothEdges)

fillVal = manageFillValue(inputImage,fillVal);

% To achieve desired edge behavior, pad input image with fill values so
% that fill values will be interpolated with source image values at the
% edges. Account for this effect by also including the added extents in the
% spatial referencing object associated with inputImage, since we've added to the
% world extent of inputImage.
if(SmoothEdges)
    [inputImage,Rin] = padImage(inputImage,Rin,fillVal);
end

% NOTE: tComp is in post-multiply form here.
tComp = composeTransformation(tform, Rin, Rout);
T = inv(tComp);
T = T(1:3,1:2);

% Convert types to match IPP support
[inputImage,origClass] = castImageForIPPUse(inputImage);

% Handle complex inputs by simply calling into IPP twice with the real and
% imaginary parts.

fillVal = double(fillVal);
if isreal(inputImage)
    outputImage = images.internal.builtins.ippgeotrans(inputImage,double(T),Rout.ImageSize,interp,fillVal);
else
    outputImage = complex( images.internal.builtins.ippgeotrans(real(inputImage),double(T),Rout.ImageSize,interp,real(fillVal)),...
        images.internal.builtins.ippgeotrans(imag(inputImage),double(T),Rout.ImageSize,interp,imag(fillVal)));
end

outputImage = cast(outputImage,origClass);

function TF = isProblemSizeTooBig(inputImage)
% IPP cannot handle double-precision inputs that are too big. Switch to
% using MATLAB's interp2 when the image is double-precision and is too big.

imageIsDoublePrecision = isa(inputImage,'double');

padSize = 3;
numel2DInputImage = (size(inputImage,1) + 2*padSize) * (size(inputImage,2) + 2*padSize);

% The size threshold is double(intmax('int32'))/8. The double-precision
% IPP routine can only handle images that have fewer than this many pixels.
% This is hypothesized to be because they use an int to hold a pointer
% offset for the input image. This overflows when the offset becomes large
% enough that ptrOffset*sizeof(double) exceeds intmax.
sizeThreshold = 2.6844e+08;
TF = imageIsDoublePrecision && (numel2DInputImage>=sizeThreshold);

function [intermediateImage, origClass] = castImageForIPPUse(inputImage)

origClass = class(inputImage);
intermediateImage = inputImage;
if(islogical(inputImage))
    intermediateImage = uint8(inputImage);
elseif(isa(inputImage,'int8') || isa(inputImage,'int16'))
    intermediateImage = single(inputImage);
elseif(isa(inputImage,'uint32') || isa(inputImage,'int32'))
    intermediateImage = double(inputImage);
end

function fillVal = manageFillValue(inputImage,fillVal)
% Always represent fillValue as a Mx1 vector if the inputImage is
% multi-channel.

if (~ismatrix(inputImage) && isscalar(fillVal))
    % If we are doing plane at time behavior, make sure fillValues
    % always propagates through code as a matrix of size determine by
    % dimensions 3:end of inputImage.
    sizeInputImage = size(inputImage);
    if (ndims(inputImage)==3)
        % This must be handled as a special case because repmat(X,N)
        % replicates a scalar X as a NxN matrix. We want a Nx1 vector.
        sizeVec = [sizeInputImage(3) 1];
    else
        sizeVec = sizeInputImage(3:end);
    end
    fillVal = repmat(fillVal,sizeVec);
end

function tComp = composeTransformation(tform, Rin, Rout)
% Form composite transformation that defines the forward transformation
% from intrinsic points in the input image in the Intel intrinsic system to
% intrinsic points in the output image in the Intel intrinsic system. This
% composite transform accounts for the spatial referencing of the input and
% output images, and differences between the MATLAB and Intel intrinsic
% systems.
%
% NOTE: tComp is returned in post-multiply form.

% The intrinsic coordinate system of IPP is 0 based. 0,0 is the location of
% the center of the first pixel in IPP. We must translate by 1 in each
% dimension to account for this.
tIntelIntrinsicToMATLABIntrinsic = [1 0 0; 0 1 0; 1 1 1];

% Define affine transformation that maps from intrinsic system of
% output image to world system of output image.
tIntrinsictoWorldOutput = Rout.TransformIntrinsicToWorld;

% Define affine transformation that maps from world system of
% input image to intrinsic system of input image.
tWorldToIntrinsicInput = Rin.TransformWorldToIntrinsic;

% Transform from intrinsic system of MATLAB to intrinsic system of Intel.
tMATLABIntrinsicToIntelIntrinsic = [1 0 0; 0 1 0; -1 -1 1];

% Form composite transformation that defines the forward transformation
% from intrinsic points in the input image in the Intel intrinsic system to
% intrinsic points in the output image in the Intel intrinsic system. This
% composite transform accounts for the spatial referencing of the input and
% output images, and differences between the MATLAB and Intel intrinsic
% systems.
tComp = tIntelIntrinsicToMATLABIntrinsic*tIntrinsictoWorldOutput / tform.A' * tWorldToIntrinsicInput*tMATLABIntrinsicToIntelIntrinsic;

% Copyright 2012-2022 The MathWorks, Inc.