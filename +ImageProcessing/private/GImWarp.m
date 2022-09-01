function [outputImage,outputRef] = GImWarp(varargin)
%IMWARP Apply geometric transformation to image.
%   B = IMWARP(A,TFORM) transforms the image A according to the geometric
%   transformation defined by TFORM, which is a geometric transformation
%   object. B is the output image. TFORM can be a 2-D or 3-D geometric
%   transformation. If TFORM is 2-D and ndims(A) > 2, such as for an RGB
%   image, then the same 2-D transformation is automatically applied to all
%   2-D planes along the higher dimensions. If TFORM is 3-D, then A must be
%   a 3-D image volume.
%
%   [B, RB] = IMWARP(A,RA,TFORM) transforms a spatially referenced
%   image specified by the image data A and the associated spatial
%   referencing object RA. When TFORM is a 2-D geometric transformation,
%   RA must be a 2-D spatial referencing object. When TFORM is a 3-D
%   geometric transformation, RA must be a 3-D spatial referencing
%   object. The output is a spatially referenced image specified by the
%   image data B and the associated spatial referencing object RB.
%
%   B = IMWARP(...,INTERP) specifies the form of interpolation to
%   use.  INTERP can be one of the strings 'nearest', or 'linear'.
%
%   [B,RB] = IMWARP(...,PARAM1,VAL1,PARAM2,VAL2,...)
%   specifies parameters that control various aspects of the geometric
%   transformation. Parameter names can be abbreviated, and case does not
%   matter.
%
%   B = IMWARP(A,D) transforms the input image A according to the
%   displacement field defined by D. B is the output image. D is an MxNx2
%   numeric matrix when the input image A is 2-D. D is an MxNxPx3 matrix
%   when the input image A is 3-D. Plane at a time behavior is also
%   supported when A is MxNxP and the input displacement field is MxNx2, in
%   which case D is applied to A one plane at a time. The first plane of
%   the displacement field, D(:,:,1) describes the X component of additive
%   displacement that is added to column,row locations in D to produce
%   remapped locations in A. Similarly, D(:,:,2) describes the Y component
%   of additive displacement values and in the 3-D case, D(:,:,:,3) describes
%   the Z component of additive displacement. The unit of displacement
%   values in D is pixels. D defines the grid size and location of the
%   output image. It is assumed that D is referenced to the default
%   intrinsic coordinate system. Both A and D are gpuArray type.
%
%   Parameters include:
%
%   'OutputView'        An imref2d or imref3d object. The ImageSize,
%                       XWorldLimits, and YWorldLimits properties of the
%                       specified spatial referencing object define the
%                       size of the output image and the location of the
%                       output image in the world coordinate system.
%
%   'FillValues'        An array containing one or several fill values.
%                       Fill values are used for output pixels when the
%                       corresponding inverse transformed location in the
%                       input image is completely outside the input image
%                       boundaries.
%
%                       If A is 2-D then 'FillValues' must be a
%                       scalar. If A is 3-D and the geometric
%                       transformation is 3-D, then 'FillValues' must be a
%                       scalar. If A is N-D and the geometric
%                       transformation is 2-D, then 'FillValues' may be
%                       either scalar or an array whose size matches
%                       dimensions 3 to N of A. For example, if A is a
%                       uint8 RGB image that is 200-by-200-by-3, then
%                       'FillValues' can be a scalar or a 3-by-1 array. In
%                       this RGB image example, possibilities for
%                       'FillValues' include:
%
%                           0                 - fill with black
%                           [0;0;0]           - also fill with black
%                           255               - fill with white
%                           [255;255;255]     - also fill with white
%                           [0;0;255]         - fill with blue
%                           [255;255;0]       - fill with yellow
%
%                       If A is 4-D with size 200-by-200-by-3-by-10, then
%                       'FillValues' can be a scalar or a 3-by-10 array.
%
%   Class Support
%   -------------
%   A can be of any nonsparse numeric class. A can also be logical.  The
%   class of B is the same as the class of A. TFORM is a geometric
%   transformation object. RA and RB are spatial referencing objects of
%   class imref2d or imref3d. D is displacement field of nonsparse 
%   numeric class.
%
%   Notes
%   -----
%   - The gpuArray version of imwarp only supports 'nearest' and 'linear'
%   interpolation.
%
%   - The gpuArray version of imwarp always uses 'SmoothEdges' false.
%
%   Example 1
%   ---------
%   % Apply a horizontal shear to an intensity image.
%
%       I = gpuArray(imread('cameraman.tif'));
%       tform = affine2d([1 0 0; .5 1 0; 0 0 1]);
%       J = imwarp(I,tform);
%       figure
%       imshow(J)
%
%   Example 2
%   ---------
%   % Apply a rotation transformation to a 3-D MRI dataset.
%
%       s = load('mri');
%       mriVolume = squeeze(s.D);
%
%       % Form rotation transformation about Y axis
%       theta = pi/8;
%       t = [cos(theta)  0      -sin(theta)   0
%           0             1              0     0
%           sin(theta)    0       cos(theta)   0
%           0             0              0     1]
%
%       tform = affine3d(t);
%       mriVolumeRotated = imwarp(gpuArray(mriVolume),tform);
%
%       figure, volshow(mriVolume);
%       figure, volshow(gather(mriVolumeRotated));
%
%   See also AFFINE2D, AFFINE3D, PROJECTIVE2D, IMREF2D, IMREF3D,
%   IMREGTFORM, GEOMETRICTRANSFORM2D, GEOMETRICTRANSFORM3D

%   Copyright 2019-2020 The MathWorks, Inc.

narginchk(2,inf);

isDisplacementFieldSyntax = isnumeric(varargin{2});

if isDisplacementFieldSyntax
%Handle displacement field syntaxes as a separate case.

    if ~isa(varargin{2},"gpuArray")
        varargin{2} = gpuArray(varargin{2});
    end
    
    parsedInputs = images.geotrans.internal.imwarpDisplacementParseInputs(varargin{:});
    manageUnsupportedGPUOptions(parsedInputs);

    method = parsedInputs.InterpolationMethod;
    fillValues = parsedInputs.FillValues;
    displacementField = cast(parsedInputs.DisplacementField,'double');
    fillValues = cast(fillValues, 'like', parsedInputs.InputImage);
    
    % check 2d/3d senario
    if ndims(displacementField)==4
        R_A = imref3d(size(parsedInputs.InputImage));
    else
        R_A = imref2d(size(parsedInputs.InputImage));
    end
    
    %check agreement of fillValues with dimensionality of problem
    dimensionality = ndims(displacementField) - 1;
    images.internal.checkFillValues(fillValues,parsedInputs.InputImage,dimensionality);    
    
    sizeD = size(displacementField);
    if ndims(displacementField) == 3
        outputRef = imref2d(sizeD(1:2));
    else
        outputRef = imref3d(sizeD(1:3));
    end
    
    %set the output image size same as input image.
    outputImage = remapPointsAndResample(parsedInputs.InputImage,R_A,displacementField,outputRef,method,fillValues,true);

    outputImage = cast(outputImage,'like',parsedInputs.InputImage);
    
    return;
end

[R_A, varargin] = images.geotrans.internal.preparseSpatialReferencingObjects(varargin{:});

parsedInputs = images.geotrans.internal.imwarpParseInputs(varargin{:});

manageUnsupportedGPUOptions(parsedInputs);

method = parsedInputs.InterpolationMethod;
fillValues = parsedInputs.FillValues;
tform = parsedInputs.GeometricTransform;

if isa(tform,'rigid2d')
    tform = affine2d(tform.T);
elseif isa(tform,'rigid3d')
    tform = affine3d(tform.T);
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

outputImage = remapPointsAndResample(parsedInputs.InputImage,R_A,tform,outputRef,method,fillValues,false);

outputImage = cast(outputImage,'like',parsedInputs.InputImage);

end

function manageUnsupportedGPUOptions(parsedInputs)

if parsedInputs.SmoothEdges
    error(message('images:imwarp:gpuSmoothEdgesNotSupported'));
end

if string(parsedInputs.InterpolationMethod) == "cubic"
    error(message('images:imwarp:gpuCubicInterpolationNotSupported'));
end

end


function outputImage = remapPointsAndResample(inputImage,R_A,tform,Rout,method,fillValues,isDisplacementFieldSyntax)

if isa(tform,'affine3d') || ndims(tform)==4||tform.Dimensionality==3
    outputImage = remapAndResampleInvertible3d(inputImage,R_A,tform,Rout,method,fillValues,isDisplacementFieldSyntax);
else
    outputImage = remapAndResampleInvertible2d(inputImage,R_A,tform,Rout,method,fillValues,isDisplacementFieldSyntax);
end

end

function checkOutputViewAgreementWithTform(Rout,tform)

if (tform.Dimensionality == 3) && ~isa(Rout,'imref3d') || ((tform.Dimensionality==2) && isa(Rout,'imref3d'))
    error(message('images:imwarp:outputViewTformDimsMismatch','''OutputView'''));
end

end

function R_out = calculateOutputSpatialReferencing(R_A,tform)
% Applies geometric transform to input spatially referenced grid to figure
% out the resolution and world limits after application of the forward
% transformation.
R_out = images.spatialref.internal.applyGeometricTransformToSpatialRef(R_A,tform);

end


function B = remapAndResampleInvertible2d(in,Rin,tform,Rout,method,fillValues,isDisplacementFieldSyntax)

% Form an input grid on the GPU, down-sampling the first two dimensions but
% keeping all of the third.

x = gpuArray.colon(1,1,Rout.ImageSize(2));
y = gpuArray.colon(1,1,Rout.ImageSize(1))';

if ndims(in) > 2
    lengthZ = prod(size(in,3:ndims(in)));
else
    lengthZ = 1;
end
z = reshape(gpuArray.colon(1,1,lengthZ),1,1,[]);

% Manage plane at a time behavior for dims 4-N.
origSize = size(in);
if ndims(in) > 3
    in = reshape(in,[size(in,1),size(in,2),lengthZ]);
end

% Collapse spatial referencing and specified transformation into a single
% composite transformation matrix. Not necessary for displacement field.

origClass = string(underlyingType(in));
numRows = size(in,1);
numCols = size(in,2);

if isscalar(fillValues) && (size(in,3) > 1)
   fillValues = repmat(fillValues,1,size(in,3)); 
end

% Take displacement field case totally seperate
if isDisplacementFieldSyntax
    if origClass == "double"
        fillValue = double(fillValues);
    else
        in = single(in);
        fillValue = single(fillValues);
    end
    
    switch method
        case 'linear'
            B = arrayfun(@warpElemDisplacementLinear, x, y, z);
        case 'nearest'
            B = arrayfun(@warpElemDisplacementNearest, x, y, z);
        otherwise
            assert(0, sprintf('Unexpected error: Incorrect method %s',method));
    end
    
    B = cast(B,origClass);

    % Manage plane at a time behavior for dims 4-N.
    if length(origSize) > 3
        B = reshape(B,[size(B,1),size(B,2),origSize(3:end)]);
    end
    
    return;
end

tComp = Rout.TransformIntrinsicToWorld / tform.T * Rin.TransformWorldToIntrinsic; 
tformMatrix = tComp;

if (origClass ~= "double")
    in = single(in);
    fillValue = single(fillValues);
    tformMatrix = single(tformMatrix);
else
    fillValue = double(fillValues);
    tformMatrix = double(tformMatrix);
end

switch(method)
        
    case 'linear'
        if isa(tform,'affine2d')
            B = arrayfun(@warpElemAffineLinear, x, y, z);
        else
            B = arrayfun(@warpElemProjectiveLinear, x, y, z);
        end
    case 'nearest'
        if isa(tform,'affine2d')
            B = arrayfun(@warpElemAffineNearest, x, y, z);
        else
            B = arrayfun(@warpElemProjectiveNearest, x, y, z);
        end
        
    otherwise
        assert(0, sprintf('Unexpected error: Incorrect method %s',method));
end
        
B = cast(B,origClass);

% Manage plane at a time behavior for dims 4-N.
if length(origSize) > 3
    B = reshape(B,[size(B,1),size(B,2),origSize(3:end)]);
end
    
    function [u,v] = mapPointDisplacement(xLoc,yLoc)
        u = xLoc + tform(yLoc,xLoc,1);
        v = yLoc + tform(yLoc,xLoc,2);
    end
    
    function [u,v] = mapPointAffine(xLoc,yLoc)
        u = xLoc*tformMatrix(1,1)+yLoc*tformMatrix(2,1)+tformMatrix(3,1);
        v = xLoc*tformMatrix(1,2)+yLoc*tformMatrix(2,2)+tformMatrix(3,2);
    end

    function [u,v] = mapPointProjective(xLoc,yLoc)
        u = xLoc*tformMatrix(1,1)+yLoc*tformMatrix(2,1)+tformMatrix(3,1);
        v = xLoc*tformMatrix(1,2)+yLoc*tformMatrix(2,2)+tformMatrix(3,2);
        w = xLoc*tformMatrix(1,3)+yLoc*tformMatrix(2,3)+tformMatrix(3,3);
        u = u/w;
        v = v/w;
    end

    function TF = outOfBounds(u,v)
        TF = (u < 1) || (u > (numCols)) ||...
             (v < 1) || (v > (numRows));
    end

    function vOut = interp2dLinear(u,v,zLoc)
        
        if outOfBounds(u,v)
            vOut = fillValue(zLoc);
        else
            vNorth = floor(v);
            vSouth = ceil(v);
            uWest = floor(u);
            uEast = ceil(u);
            
            weightNorth = triangle(vNorth-v);
            weightSouth = 1-weightNorth;
            valueWest = weightNorth*in(vNorth,uWest,zLoc) + weightSouth*in(vSouth,uWest,zLoc);
            valueEast = weightNorth*in(vNorth,uEast,zLoc) + weightSouth*in(vSouth,uEast,zLoc);
            
            weightWest = triangle(uWest-u);
            weightEast = 1-weightWest;
            
            vOut = weightWest*valueWest + weightEast*valueEast;
            
        end
        
    end

    function vOut = interp2dNearest(u,v,zLoc)
        
        if outOfBounds(u,v)
            vOut = fillValue(zLoc);
        else
            vOut = in(round(v),round(u),zLoc);
        end
        
    end

    function vOut = warpElemDisplacementLinear(xLoc,yLoc,zLoc)
        [u,v] = mapPointDisplacement(xLoc,yLoc);
        vOut = interp2dLinear(u,v,zLoc);
    end

    function vOut = warpElemDisplacementNearest(xLoc,yLoc,zLoc)
        [u,v] = mapPointDisplacement(xLoc,yLoc);
        vOut = interp2dNearest(u,v,zLoc);
    end

    function vOut = warpElemAffineNearest(xLoc,yLoc,zLoc)
        [u,v] = mapPointAffine(xLoc,yLoc);
        vOut = interp2dNearest(u,v,zLoc);
    end

    function vOut = warpElemAffineLinear(xLoc,yLoc,zLoc)
        [u,v] = mapPointAffine(xLoc,yLoc);
        vOut = interp2dLinear(u,v,zLoc);
    end

    function vOut = warpElemProjectiveNearest(xLoc,yLoc,zLoc)
        [u,v] = mapPointProjective(xLoc,yLoc);
        vOut = interp2dNearest(u,v,zLoc);
    end

    function vOut = warpElemProjectiveLinear(xLoc,yLoc,zLoc)
        [u,v] = mapPointProjective(xLoc,yLoc);
        vOut = interp2dLinear(u,v,zLoc);
    end

end


function B = remapAndResampleInvertible3d(in,Rin,tform,Rout,method,fillValues,isDisplacementFieldSyntax)

% Form an input grid on the GPU, down-sampling the first two dimensions but
% keeping all of the third.
x = gpuArray.colon(1,1,Rout.ImageSize(2));
y = gpuArray.colon(1,1,Rout.ImageSize(1))';
z = reshape(gpuArray.colon(1,1,Rout.ImageSize(3)),1,1,[]);

HigherSize=size(in,4:max(ndims(in),4));
t=reshape(gpuArray.colon(1,1,prod(HigherSize)),1,1,1,[]);
if isscalar(fillValues)
	fillValues=repmat(fillValues,1,numel(t));
end

% Collapse spatial referencing and specified transformation into a single
% composite transformation matrix.

origClass = string(underlyingType(in));
numRows = size(in,1);
numCols = size(in,2);
numPlanes = size(in,3);

% take displacement field input totally seperate
if isDisplacementFieldSyntax
    if origClass == "double"
        fillValue = double(fillValues);
    else
        in = single(in);
        fillValue = single(fillValues);
    end
    
    switch method
        case 'linear'
            B = arrayfun(@warpElemDisplacementLinear, x, y, z,t);
        case 'nearest'
            B = arrayfun(@warpElemDisplacementNearest, x, y, z,t);
        otherwise
            assert(0, sprintf('Unexpected error: Incorrect method %s',method));
    end
    
    B = cast(B,origClass);
    return;
end
    
tComp = Rout.TransformIntrinsicToWorld / tform.T * Rin.TransformWorldToIntrinsic;
tformMatrix = tComp;

if (origClass ~= "double")
    in = single(in);
    fillValue = single(fillValues);
    tformMatrix = single(tformMatrix);
else
    fillValue = double(fillValues);
    tformMatrix = double(tformMatrix);
end

switch(method)
    
    case 'linear'
        B = arrayfun(@warpElemLinear, x, y, z,t);
    case 'nearest'
        B = arrayfun(@warpElemNearest, x, y, z,t);
    otherwise
        error('');
end

B = reshape(cast(B,origClass),[size(B,1:3),HigherSize]);

    function [u,v,w] = mapPoint(xLoc,yLoc,zLoc)
        u = xLoc*tformMatrix(1,1)+yLoc*tformMatrix(2,1)+zLoc*tformMatrix(3,1)+tformMatrix(4,1);
        v = xLoc*tformMatrix(1,2)+yLoc*tformMatrix(2,2)+zLoc*tformMatrix(3,2)+tformMatrix(4,2);
        w = xLoc*tformMatrix(1,3)+yLoc*tformMatrix(2,3)+zLoc*tformMatrix(3,3)+tformMatrix(4,3);
    end
    
    function [u,v,w] = mapPointDisplacement(xLoc,yLoc,zLoc)
        u = xLoc + tform(yLoc,xLoc,zLoc,1);
        v = yLoc + tform(yLoc,xLoc,zLoc,2);
        w = zLoc + tform(yLoc,xLoc,zLoc,3);
    end

    function TF = outOfBounds(u,v,w)
        TF = (u < 1) || (u > (numCols)) ||...
             (v < 1) || (v > (numRows)) ||...
             (w < 1) || (w > (numPlanes));
    end

    function vOut = interp3dLinear(u,v,w,tLoc)
        
        if outOfBounds(u,v,w)
            vOut = fillValue(tLoc);
        else
            vNorth = floor(v);
            vSouth = ceil(v);
            uWest = floor(u);
            uEast = ceil(u);
            wUp = floor(w);
            wDown = ceil(w);
            
            weightNorth = triangle(vNorth-v);
            weightSouth = 1-weightNorth;
            weightWest = triangle(uWest-u);
            weightEast = 1-weightWest;
            weightUp = triangle(wUp-w);
            weightDown = 1-weightUp;
            
            valueWest = weightNorth*in(vNorth,uWest,wUp,tLoc) + weightSouth*in(vSouth,uWest,wUp,tLoc);
            valueEast = weightNorth*in(vNorth,uEast,wUp,tLoc) + weightSouth*in(vSouth,uEast,wUp,tLoc);
            
            v1 = weightWest*valueWest + weightEast*valueEast;
            
            valueWest = weightNorth*in(vNorth,uWest,wDown,tLoc) + weightSouth*in(vSouth,uWest,wDown,tLoc);
            valueEast = weightNorth*in(vNorth,uEast,wDown,tLoc) + weightSouth*in(vSouth,uEast,wDown,tLoc);
            
            v2 = weightWest*valueWest + weightEast*valueEast;
            
            vOut = weightUp*v1 + weightDown * v2;

        end
    
    end

    function vOut = interp3dNearest(u,v,w,tLoc)
        
        if outOfBounds(u,v,w)
            vOut = fillValue(tLoc);
        else
            vOut = in(round(v),round(u),round(w),tLoc);            
        end
        
    end

    function vOut = warpElemDisplacementLinear(xLoc,yLoc,zLoc,tLoc)
        
        [u,v,w] = mapPointDisplacement(xLoc,yLoc,zLoc);    
        vOut = interp3dLinear(u,v,w,tLoc);
        
    end

    function vOut = warpElemDisplacementNearest(xLoc,yLoc,zLoc,tLoc)
        
        [u,v,w] = mapPointDisplacement(xLoc,yLoc,zLoc);    
        vOut = interp3dNearest(u,v,w,tLoc);
        
    end
    
    function vOut = warpElemLinear(xLoc,yLoc,zLoc,tLoc)
        
        [u,v,w] = mapPoint(xLoc,yLoc,zLoc);
        vOut = interp3dLinear(u,v,w,tLoc);
        
    end

    function vOut = warpElemNearest(xLoc,yLoc,zLoc,tLoc)
        
        [u,v,w] = mapPoint(xLoc,yLoc,zLoc);
        vOut = interp3dNearest(u,v,w,tLoc);
        
    end
    
end

function f = triangle(x)

f = (x+1) .* ((-1 <= x) & (x < 0)) + (1-x) .* ((0 <= x) & (x <= 1));

end