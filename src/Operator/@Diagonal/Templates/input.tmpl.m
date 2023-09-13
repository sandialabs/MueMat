#ifndef INV
     AMAT = A.GetMatrixData();
     if isempty(AMAT)
       ST=dbstack; error('%s:: Missing matrix data.\n', ST(1).name);
     end
#else
     FactorBlkDiag(A); % This do nothing if the factorization is already computed
     AMAT = A.GetDiagonalInvData();
#endif

