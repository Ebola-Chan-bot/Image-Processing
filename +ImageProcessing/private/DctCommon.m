function arg1=DctCommon(Fun,arg1,mrows,ncols)
[m, n] = size(arg1,1,2);
% Basic algorithm.
if (nargin == 2)
	if m > 1 && n > 1
		arg1 = pagetranspose(Fun(pagetranspose(Fun(arg1))));
		return;
	else
		mrows = m;
		ncols = n;
	end
end

% Padding for vector input.
if nargin==3, ncols = mrows(2); mrows = mrows(1); end

if m == 1 && mrows > m, arg1(2, 1,:) = 0; m = 2; end
if n == 1 && ncols > n, arg1(1, 2,:) = 0; n = 2; end
if m == 1, mrows = ncols; ncols = 1; end   % For row vector.

% Transform.

arg1 = Fun(arg1, mrows);
if m > 1 && n > 1, arg1 = pagetranspose(Fun(pagetranspose(arg1), ncols)); end