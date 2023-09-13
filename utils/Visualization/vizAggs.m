function vizAggs(MgHierarchy,level,markerSize)
  % Produces plots of the aggregates (2D)
  %
  % input:
  %   MgHierarchy -- an AMG class object
  %   level       -- integer level number, defaults to 2.
  %   markerSize  -- specifies the size of the marker in points.
  % output:
  %   plot of the aggregates
  %
  % Example usage:
  %   vizAggs(MgHierarchy,1,8);
  if isa(MgHierarchy,'Hierarchy') ~= true
    error('The first argument to ''vizAggs'' must be an ''Hierarchy'' class.');
  end
  if ~varexist('level'), level = 2; end
  if ~varexist('markerSize'), markerSize = 8; end

  CurrentLevel = MgHierarchy.GetLevel(level);

  if ~CurrentLevel.IsAvailable('Aggregates')
    error('The Aggregate structure is not available. Did you save the aggregates? Does the hierarchy have enough levels?');
  end
  
  % Data
  Aggregates = CurrentLevel.Get('Aggregates');
  X = CurrentLevel.Get('xcoords');
  Y = CurrentLevel.Get('ycoords');

  % Associate each aggregate with a color (randomly)
  c = 'byrkmgc';
  numColors = length(c);
  
  numAggs = size(Aggregates.NodesInAgg,1);
  for ii=1:numAggs,
    aggColors(ii) = randi(numColors);
  end

  numNodes = length(Aggregates.AggId);
  for ii=1:numNodes
    aggColor(ii) = aggColors(Aggregates.AggId(ii));
  end

  % Plot aggregates
  Z = zeros(length(X),1); % 2D
  plot3k([X Y Z],aggColor,'s',[],markerSize);
  hold on
  
  % Distinguish root points
  plot(X(Aggregates.Roots),Y(Aggregates.Roots),'wo','MarkerSize',6,'MarkerFaceColor','w');;
  % plot(X(Aggregates.Roots),Y(Aggregates.Roots),'w.','MarkerSize',markerSize-2);

  view(2);

end %function vizAggs()
