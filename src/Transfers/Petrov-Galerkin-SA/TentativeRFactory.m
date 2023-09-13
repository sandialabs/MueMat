%% RFactory
% small example for a concrete simple implementation of the RFactory
% interface. class creates a tentative restriction operator using the
% tentative prolongator.
%%

classdef TentativeRFactory < RFactory
    % class demonstrates the use of SavePtent need
    % usage:
    % PRfact = GenericPRFactory(SaPFactory(), TentativeRFactory());
    % alternatively you can use:
    % Rfact = TransPFactory();
    % Rfact.usePtent(true);
    % PRfact = GenericPRFactory(SaPFactory(), Rfact);
    methods
        function this = TentativeRFactory()
            % since we don't want to recompute the tentative prolongator
            % for the restrictor, we need the tentative prolongator stored
            % on all levels
% DEPRECATED
%            Need.SavePtent ='all';
            if this.GetOutputLevel() > 5,
                fprintf('TransPFactory: Need Tentative Prolongator\n');
            end
%            this.AddNeeds(Need);
        end

        function flag = BuildR(this,FineLevel,CoarseLevel,Specs)
            flag = true;

            this.TempOutputLevel(Specs);

            Pmat = CoarseLevel.Get2('Ptent');
            if isempty(Pmat)
                error('we need Ptent from the CoarseLevel!');
            end

            CoarseLevel.Set('R', Pmat');

            this.RestoreOutputLevel();
        end
    end

    methods (Access = protected)
        function Copy_(this,src,mc)
            %COPY_
            %
            %   SYNTAX   obj.Copy_(src, mc);
            %
            %     src - Object to copy
            %     mc  - MATLAB Metaclass
            [cmd, data, mc] = this.CopyCmd_(src,mc);
            eval(cmd);
        end
    end
end
