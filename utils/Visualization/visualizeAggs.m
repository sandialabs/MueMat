function visualizeAggs(amg,coords)

    % make sure coords are "row vectors"
    % first row = x coords
    % second row = y coords
    if size(coords,1) > size(coords,2)
        coords = coords';
    end
    subplot(1,length(amg.Levels_),1);
    fcol = [0.5,0.5+0.5*rand(1,1),0.5*rand(1,1)];
    plot(coords(1,:),coords(2,:),'Marker','.','MarkerEdgeColor',fcol,'LineStyle','none');

    for i = 1:length(amg.Levels_)-1
        subplot(1,length(amg.Levels_),i+1);
        aggInfo = amg.GetLevel(i).Get('Aggregates');
        visualizeLevelAggs(aggInfo,coords);
        coords = [coords(1,aggInfo.Roots);coords(2,aggInfo.Roots)];
    end




end

%% visualizes aggregates using aggInfo and coords
% for current level
function visualizeLevelAggs(aggInfo,coords)
    hold on

    % loop over all aggregates
    for aggId = 1 : size(aggInfo.NodesInAgg,1)
        fcol = [0.5,0.5+0.5*rand(1,1),0.5*rand(1,1)];
        nodeIds = find(aggInfo.NodesInAgg(aggId,:));
        aggCoords = coords(:,nodeIds);

        if(size(aggCoords,2) > 1)
           if length(nodeIds) > 2,
              K = convhull(aggCoords(1,:),aggCoords(2,:));
              fill(aggCoords(1,K),aggCoords(2,K),fcol);
           end
        else
           plot(aggCoords(1,:),aggCoords(2,:),'rx');
           alpha(0.5);
        end
    end

    plot(coords(1,aggInfo.Roots),coords(2,aggInfo.Roots),'x');

    hold off
end
