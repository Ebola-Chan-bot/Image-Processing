function C = CNormXCorr2(T,A)
A = shiftData(A);
T = shiftData(T);
[m,n] = size(T,1,2);
[p,q] = size(A,1,2);
R = m+p-1;
S = n+q-1;
K_fft = 3.3e-7;
Tr = K_fft*R*log(R);
if S==R
	Ts = Tr;
else
	Ts = K_fft*S*log(S);
end
CTM=coder.target('MATLAB');
if (prod([2.7e-8,m,n,p,q]) < 3*(S*Tr + R*Ts))
	xcorr_TA = conv2(rot90(T,2),A);
else
	F=fft2(rot90(T,2),R,S).*fft2(A,R,S);
	if CTM
		xcorr_TA = ifft2(F,'symmetric');
	else
		xcorr_TA = ifft2(F,'nonsymmetric');
	end
end
mn = m*n;
local_sum_A = local_sum(A,m,n);
denom = sqrt(mn-1)*std(T,0,1:2).*sqrt( max( local_sum(A.*A,m,n) - (local_sum_A.^2)/mn,0) );
numerator = (xcorr_TA - local_sum_A.*sum(T,1:2)/mn );
if CTM
	C = zeros(size(numerator));
else
	C = zeros(size(numerator),'like',xcorr_TA);
end
i_nonzero = denom > sqrt( eps( max(abs(denom(:)))) );
C(i_nonzero) = numerator(i_nonzero) ./ denom(i_nonzero);
C( ( abs(C) - 1 ) > sqrt(eps(1)) ) = 0;
if ~CTM
	C = real(C);
end
function local_sum_A = local_sum(A,m,n)
local_sum_A = cumsum(padarray(A,[m n]),1);
local_sum_A = cumsum(local_sum_A(1+m:end-1,:,:)-local_sum_A(1:end-m-1,:,:),2);
local_sum_A = local_sum_A(:,1+n:end-1,:)-local_sum_A(:,1:end-n-1,:);
local_sum_A=reshape(local_sum_A,[size(local_sum_A,1:2),size(A,3:ndims(A))]);
function A = shiftData(A)
if ~(isa(A,'uint8') || isa(A,'uint16') || isa(A,'uint32'))
	min_A = min(A(:));
	if min_A < 0
		A = A - min_A;
	end
end
A=double(A);