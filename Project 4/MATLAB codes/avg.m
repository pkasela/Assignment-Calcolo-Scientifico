function B=avg(A)


if nargin<2, k = 1; end
if size(A,1)==1, A = A'; end
if k<2, B = (A(2:end,:)+A(1:end-1,:))/2; else, B = avg(A,k-1); end
if size(A,2)==1, B = B'; end



%if(size(A,1)==1), A=A'; end %se A riga allora la traspongo in colonna


%B=(A(1:end-1,:)+A(2:end,:))/2; %media tra due elementi consecutivi per colonna


%if(size(A,2)==1), B=B'; end %se A era riga, la metto come prima 

return
%It work with martix, column and row

%% EXAMPLE
% x=0:0.1:1;
% xa=avg(x);

%% OSS

%U_a=(avg(U'))' for an average per row, where U is matrix
%else it will to an average on column