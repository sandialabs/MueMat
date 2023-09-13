function [Pvec] = Flatten(Pmat, PatternInds)
  Pvec           = reshape(Pmat',[],1);
  Pvec           = Pvec(PatternInds);
end %Flatten()
