% [AN,NF]=normalize(A)
% AN = normalized version of A
% NF = normalization factor
% It does not normalize columns
% with norms less that 1e-11
% Chris Siefert 2/15/04
function [AN,NF]=normalize(A)
SZ=size(A,2);
NF=zeros(SZ,1);
AN=0*A;
for I=1:SZ,
  NF(I)=norm(A(:,I));
  if(NF(I) < 1e-11)
    NF(I)=1;AN(:,I)=A(:,I);
  else
   AN(:,I)=A(:,I)/NF(I);
  end
end
