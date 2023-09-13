function [Afiltered] = TestSchurDrop(A, Level, tol)

  Afiltered = spones(Level.Get('Arr').GetMatrixData) .* SAdrop(A, Level, tol);

