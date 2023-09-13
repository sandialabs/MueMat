function vizNS(MgHierarchy,MgType,level)
  % Produces plots of the 2D nullspace vectors associated with level "level" on a mesh
  %
  % input:
  %   MgHierarchy -- an AMG class object
  %   MgType      -- string indicating type of AMG, e.g., 'sa'.  Used in plot titles.  Defaults to "unknown".
  %   level       -- integer level number, defaults to 2.
  %
  % output:
  %   one 3D plot for each nullspace vector
  %
  % Example usage:
  %   vizCoarseNS(MgHierarchy,'emin');
  %
  % Comments:
  %   For each vector, different degrees of freedom will appear in unique colors in the same plot.
  %   Degrees of freedom that are identically zero are not plotted.
  if isa(MgHierarchy,'AMG') ~= true
    error('The first argument to ''vizCoarseNS'' must be an ''AMG'' class.');
  end
  if ~varexist('MgType') MgType = 'unknown'; end
  if isa(MgType,'char') ~= true
    error('The second argument to ''vizCoarseNS'' must be a character array.');
  end
  if ~varexist('level') level = 2; end
  levelBucket = MgHierarchy.GetLevel(level);
  numCDOF = levelBucket.Get('A').GetRowMap().ConstBlkSize();
  CNS = levelBucket.Get('NullSpace');
  numCNS = size(CNS,2);
  X = levelBucket.Get('xcoords');
  Y = levelBucket.Get('ycoords');
  c = 'brk';
  for ii=1:numCNS
    figure
    labels = '';
    for jj=1:numCDOF
      component = CNS(jj:numCDOF:size(CNS,1),ii);
      if nnz(component) > 0
        plot3(X, Y, component, [c(jj) '.']);
        labels = [labels '''' 'dof ' num2str(jj) ''',' ];
      end
      hold on;
    end
    h = eval( ['legend(' labels(1:length(labels)-1) ');'] ); %get rid of last comma
    set(h,'FontSize',15);
    title([MgType ' RBM #' num2str(ii) ' on level ' num2str(level)],'FontSize',17)
    xlabel('X','FontSize',15); ylabel('Y','FontSize',15); zlabel('Z','FontSize',15);
    grid on; box on;
  end
end %function vizCoarseNS()
