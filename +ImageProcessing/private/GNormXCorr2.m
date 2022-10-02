function C = GNormXCorr2(T,A)
A = ShiftData(gpuArray(A));
T = ShiftData(gpuArray(T));
[m,n] = size(T,1,2);
outsize = [m,n] +  size(A,1:2) - 1;
if (numel(T)<2500)
	xcorr_TA = conv2(rot90(T,2),A);
else
	xcorr_TA = ifft2(fft2(rot90(T,2),outsize(1),outsize(2)) .* fft2(A,         outsize(1),outsize(2)), 'symmetric');
end
mn    = m*n;
[numerator,denom] = arrayfun(@computeNumDen,LocalSum(A,m,n),LocalSum(A.*A,m,n),xcorr_TA,sum(T,[1,2]),sqrt(mn-1)*std(T,0,[1,2]));
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