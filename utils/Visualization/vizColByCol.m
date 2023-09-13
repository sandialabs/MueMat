function vizColByCol(MgHierarchy)
  % Visualize the difference between a fine nullspace vector
  % and the interpolated coarse grid nullspace vector.
  fine = MgHierarchy.GetLevel(1);
  coarse = MgHierarchy.GetLevel(2);

  eminP = coarse.Get('P').GetMatrixData();
  CNS = coarse.Get('NullSpace');
  FNS = fine.Get('NullSpace');

  dif = eminP*CNS - FNS;
  fine = fine.Set('NullSpace', dif);
  MgHierarchy = MgHierarchy.SetLevel(fine,1);
  vizNS(MgHierarchy,'emin',1);
  fprintf('column by column norm(FNull - Pfinal*CNull) = ');
  for jj=1:size(dif,2), fprintf('%10.2e ',norm(dif(:,jj))); end
  fprintf('\n');
end
