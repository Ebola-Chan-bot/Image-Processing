function C = GNormXCorr2(T,A)
A = shiftData(gpuArray(A));
T = shiftData(gpuArray(T));
[m,n] = size(T,1,2);
outsize = [m,n] +  size(A,1:2) - 1;
if (numel(T)<2500)
	xcorr_TA = conv2(rot90(T,2),A);
else
	xcorr_TA = ifft2(fft2(rot90(T,2),outsize(1),outsize(2)) .* fft2(A,         outsize(1),outsize(2)), 'symmetric');
end
mn    = m*n;
[numerator,denom] = arrayfun(@computeNumDen,local_sum(A,m,n),local_sum(A.*A,m,n),xcorr_TA,sum(T,[1,2]),sqrt(mn-1)*std(T,0,[1,2]));
	function [num,den] = computeNumDen(local_sum_a,local_sum_a2,xcorr_ta,sumT,denom_T)
		den            = denom_T*sqrt(max(local_sum_a2 - (local_sum_a^2)/mn,0));
		num            = xcorr_ta - local_sum_a*sumT/mn;
	end
tol = sqrt( eps( max(abs(reshape(denom,numel(denom),1)))) );
C   = arrayfun(@computeNCC,numerator,denom);
	function c = computeNCC(num,den)
		if den>tol
			c = num/den;
			c = (abs(c)-1<=sqrt(eps(1)))*c;
		else
			c=0;
		end
	end
end
function local_sum_A = local_sum(A,m,n)
local_sum_A = cumsum(padarray(A,[m n]),1);
local_sum_A = cumsum(local_sum_A(1+m:end-1,:,:)-local_sum_A(1:end-m-1,:,:),2);
local_sum_A = local_sum_A(:,1+n:end-1,:)-local_sum_A(:,1:end-n-1,:);
local_sum_A=reshape(local_sum_A,[size(local_sum_A,1:2),size(A,3:ndims(A))]);
end
function A = shiftData(A)
if ~ ismember(underlyingType(A), {'uint8','uint16','uint32'})
	min_A = min(A(:));
	if min_A < 0
		A = A - min_A;
	end
end
A = double(A);
end