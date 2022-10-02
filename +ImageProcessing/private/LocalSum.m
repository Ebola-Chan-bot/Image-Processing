function S = LocalSum(A,m,n)
S = cumsum(padarray(A,[m n]),1);
S = cumsum(S(1+m:end-1,:,:)-S(1:end-m-1,:,:),2);
S = S(:,1+n:end-1,:)-S(:,1:end-n-1,:);
S=reshape(S,[size(S,1:2),size(A,3:ndims(A))]);