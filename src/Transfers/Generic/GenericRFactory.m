% This restriction factory wraps a prolongation factory.
%
% The restriction operator is build by applying
% this.PFactory_.Build() to A' and by transposing the resulting P.
%
% This factory works with both Level and Level.
%
% Example (Unsymmetric Smoothed Aggregation):
%
% Pfact = SaPFactory(CoalesceDropFactory(),AggregationFactory());
% Rfact = GenericRFactory(Pfact);
%
% In this example, the RFact will do:
%  - fakeA = A'
%  - Get a prolongator 'fakeP' from the SAPFactory applied on fakeA
%  - R = fakeP'
%
% So, R = fakeP' = (S . fakePtent)' = fakePtent' (I - fakeA' D^{-1} w)
%                                   = Rtent      (I -  w  A  D^{-1})
%
% See also (formula 2.5):
% M. Sala and R. Tuminaro in the paper "A new Petrov-Galerkin
% smoothed aggregation preconditioner for nonsymmetric linear
% systems", SIAM Journal on Scientific Computing, 2008, 31, 143-166

classdef GenericRFactory < RFactory

    properties % TODO (Access = private)
        PFactory_
    end

    methods
        function [this] = GenericRFactory(PFactory)
            % Copy constructor
            if nargin == 1 && isa(PFactory, class(this)), this.Copy_(PFactory,[]); return; end
            %

            if ~PFactory.SupportsRestrictionMode()
                error('GenericRFactory: PFactory has no support for a restriciton mode!\n');
            end

            this.PFactory_ = PFactory;

            % Make sense to reuse aggregates by default
            %this.ReUseAggregates(true);
            % no, wen don't want to aggregate the same thing twice!

        end

        function flag = Build(this,FineLevel,CoarseLevel)

            %% switch given PFactory to restrictor mode!
            pmode = this.PFactory_.ProlongationMode(false);

            %% adapt given PFactory
            %reUseR = this.PFactory_.ReUseP(this.reUseR_);
            %JG: TODO
            %reUseRtent = this.PFactory_.ReUsePtent(this.reUseRtent_);
            %reUseAggregates = this.PFactory_.ReUseAggregates(this.reUseAggregates_);
            %reUseGraph = this.PFactory_.ReUseGraph(this.reUseGraph_);

            flag = this.PFactory_.Build(FineLevel, CoarseLevel);

            %% reset given PFactory
            %this.PFactory_.ReUseP(reUseR);
            %this.PFactory_.ReUsePtent(reUseRtent);
            %this.PFactory_.ReUseAggregates(reUseAggregates);
            %this.PFactory_.ReUseGraph(reUseGraph);

            %% switch given PFactory (back) to initial mode
            this.PFactory_.ProlongationMode(pmode);

            %% check if 'R' was successfully generated with this.PFactory_
            if ~CoarseLevel.IsAvailable('R',this.PFactory_),
                error('GenericRFactory: R could not be set by PFactory_');
            end

            % redeclare generating factory to this GenericRFactory
            CoarseLevel.Set('R', CoarseLevel.Get('R', this.PFactory_), this);
            CoarseLevel.Release('R', this.PFactory_);
        end
        function SetNeeds(this, FineLevel, CoarseLevel)
        % Obtain any cross factory specifications
        % TODO add requests for AggFactory and...
          %FineLevel.Request('NullSpace');
          %CoarseLevel.Request('NullSpace'); % the result from initialP!

          pmode = this.PFactory_.ProlongationMode(false);
          this.PFactory_.SetNeeds(FineLevel, CoarseLevel);
          this.PFactory_.ProlongationMode(pmode);

          % request for result
          CoarseLevel.Request('R',this.PFactory_);
        end
    end %public methods

    methods (Access = protected)

        function Copy_(this, src, mc)
            %COPY_
            %
            %   SYNTAX   obj.Copy_(src, mc);
            %
            %     src - Object to copy
            %     mc  - MATLAB Metaclass
            [cmd, data, mc] = this.CopyCmd_(src,mc);
            eval(cmd);
        end

    end % methods

end

% TODO: - Left and Right NullSpace ?
%       - Smart Operator Transpose: don't copy the matrix A but only change the multiply method.
%                                   + fix transpose() method of Operator for multiple views.
%       - How to know if we can destroy aggregates after RFact.Build() ?
%         possible option -> if this.SaveAggregate is false, we destroy it ?
%       - ReUseP(): error if not SA/Emin (=GenericDataBucket)
