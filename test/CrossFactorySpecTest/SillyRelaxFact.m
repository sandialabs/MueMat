% MUST BE UPDATED with the new design of smootherfactory

classdef SillyRelaxFact < SmootherFactoryBase
% This class is only used to demonstrate an aggregate-need by
% a smoother using the CrossFactorySpecs class.
%
   methods
    function [this] = SillyRelaxFact()
    % constructor makes request for aggregates.
        Need.SaveAggregates = 'all';
        this.AddNeeds(Need);
        fprintf('SillyRelaxFact: Need aggregates\n');
    end

    function [PreSmoo,PostSmoo,Amat] = Build(this, Level, Specs)
    % Build attempts to retrieve aggregates and complains if they are not found.

         this.TempOutputLevel(Specs);
         Amat = Level.Get('A');
         PreSmoo = 7; PostSmoo = 5;
         AggInfo = Level.Get('Aggregates');
         if isempty(AggInfo), 
            fprintf('SillyRelaxFact: Aggregates not found on level %d\n',...
                     Level.GetLevelId());
         elseif this.GetOutputLevel() > 5,
            fprintf('SillyRelaxFact: Setting smoother on level %d\n',...
                     Level.GetLevelId());
         end
         this.RestoreOutputLevel();
    end
   end
end % classdef
