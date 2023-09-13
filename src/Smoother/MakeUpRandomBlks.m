%
% Utility function for testing things
% It makes an arbitrary number of
% large blocks and fills this with arbitrary blocks from Amat.
% FIXME Because this function creates a new random stream each time it's called,
% FIXME I suspect that it will always return the same Collection.
%
%
function [Collection] = MakeUpRandomBlks(Amat,type)

    Collection  = [];

    PrevStream = RandStream.setGlobalStream(RandStream.create('mrg32k3a','NumStreams',1));
    RowMap = Amat.GetRowMap();
    NBlks = RowMap.NNodes();
    if strcmp(type,'NonOverlapping'),

      %
      % generate a random integer for each block row
      % numbers with the same random number are grouped
      % into the same block.
      %
      temp = ceil(rand(1,NBlks)*NBlks/4);

      [aa,bb] = unique(temp);
      Collection.NSubsets = length(aa);
      for i=1:Collection.NSubsets
         Collection.Subsets(i) = CreateDOFSubset(RowMap, 'Scattered', -1, -1, ...
                                          find(temp==aa(i)));
      end
    elseif strcmp(type,'Overlapping'),
        Collection.NSubsets = ceil(NBlks/3);
      InABlock = zeros(NBlks,1);
      for i=1:Collection.NSubsets-1
         BlkRows = find( rand(NBlks,1) < .1);
         Collection.Subsets(i) = CreateDOFSubset(RowMap, 'Scattered', -1, -1,BlkRows);
         InABlock(BlkRows) = 1;
      end
      BlkRows = find( rand(NBlks,1) < .1);
      InABlock(BlkRows) = 1;
      BlkRows = sort([ BlkRows ; find(InABlock == 0)]);
      Collection.Subsets(Collection.NSubsets) = ...
                     CreateDOFSubset(RowMap, 'Scattered', -1, -1,BlkRows);
    else
       fprintf('MakeUpRandomBlks:: Unrecognized type %s\n',type);
       keyboard;
    end
    RandStream.setGlobalStream(PrevStream);

