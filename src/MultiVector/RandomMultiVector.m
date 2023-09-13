function [mv] = RandomMultiVector(Length,numVectors,isSparse)
  % Utility for creating random multivectors.
  % See also MultiVector.
  if ~varexist('isSparse'), isSparse = false; end
  mv = MultiVector([],Length,numVectors,isSparse);
  mv.Random();
end
