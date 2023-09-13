function [Unflattened] = BlowUp(Vec,nfine,ncoarse,PInds);

nn = length(Vec);
BigVec = sparse(PInds,ones(nn,1),Vec,nfine*ncoarse,1,nn);
Unflattened = reshape(BigVec,ncoarse,nfine)';

% test
%Tvec = rand(nnz(Ppattern),1);
%Unflattened = BlowUp(Tvec,nfine,ncoarse,PInds);
%Pvec = Flatten(Unflattened,PInds);
%nnz(Pvec-Tvec),

