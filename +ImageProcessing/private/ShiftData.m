function A = ShiftData(A)
if ~ ismember(underlyingType(A), {'uint8','uint16','uint32'})
	min_A = min(A(:));
	if min_A < 0
		A = A - min_A;
	end
end
A = double(A);