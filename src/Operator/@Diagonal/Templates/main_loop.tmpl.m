         OutVecPart = FUNCNAME(FUNC,SubFunc)(AMAT, Vec, [fGlobalRow fLocalRow], [lGlobalRow lLocalRow], BlkSize, VecType);
         %TEST: OutVecTEST(fLocalRow:lLocalRow,:) = OutVecPart;

         [I,J,V] = find(OutVecPart);
         I = I + fLocalRow - 1; % shift == OutVec(fLocalRow:lLocalRow,:) = DSubFuncBlkApply(...)
         OutVecIJV.Add(i,I,J,V);
