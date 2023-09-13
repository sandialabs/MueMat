%
% CrossFactorySpecTest.m demonstrates how parameters affecting several 
% factories are set and passed. Potential situations include
%
%    1) Needing a factory to save data (which is not normally saved) 
%       so that a class, function, or factory can use it. Some examples are ...
%
%       a) aggregates saved for plotting.
%       b) tentative prolongators saved by prolongator factory for use within
%          restrictor factory.
%       c) prolongator sparsity pattern saved by prolongator factory so that a
%          related linear solve can reuse the pattern during its multigrid setup. 
%
%    2) Specify output or debugging level for all factories..
%
% One important point is that 'Needs' can come from users or be generated 
% internally. For example, saving aggregates for plotting is a user-need while
% keeping aggregates for a restriction factory is an internal-need.  In the 
% later case, the user may have no idea that a particular algorithm requires 
% aggregates to be saved.
% 
% Key Concept 1
% =============
% Needs to a factory are made via 
%
%          MyFactory.AddNeeds(MyNeeds)
%
% MyNeeds is a structure with fields used like keys.  In principle, these can 
% be anything which some function expects.  While anything is possible, we 
% highly recommend adhering to the following convections when saving/reusing:
%      
%   MyNeeds.Savexxx           = 'startlevel:endlevel' or 'all'
%   MyNeeds.PersistentSavexxx = 'startlevel:endlevel' or 'all'
%   MyNeeds.ReUsexxx          = 'startlevel:endlevel' or 'all'
%
% where 'xxx' refers to a specific object (e.g. 'Aggregates'), 'all' saves or
% retrieves on all levels while startlevel/endlevel provide level specific
% control. 'Persistent' means that the object is not removed by Hierarchy.Clean()
% after solving. This allows for data reuse between multiple linear solves.
% 'Save' without 'Persistent' means that an object is only saved during a single
% linear solve and is then removed by Hierarchy.Clean(). This allows one factory to
% use something created by another factory.
%
% Key Concept 2
% =============
% The need mechanism is only for parameters that affect multiple factories.
% Thus, a factory specific set command (e.g., BlkSmootherFactory.SetOmega(...))
% should be used to set something particular to just one factory.
%
% The only exception is reuse. Reuse is both cross-factory (restrictor algorithm
% wants a tentative prolongator) and single factory (prolongator scheme reuses
% the same aggregates for different linear solves). Theoretically, a factory
% specific command (something like SaPFactory.SaveAggregates()) could be used.
% However, since there is already a cross-factory mechanism for 'needing' 
% aggregates, it is employed by a factory for reuse by itself.
%
% Key Concept 3
% =============
% Individual needs are not directly accessed by a factory's build method.
% Instead, MyFactory.Build( ....) is passed a CrossFactory object that 
% defines the state. This is because a factory's needs do not include other
% factory needs.  CrossFactory is created by Hierarchy's hierarchy populating
% functions (FillHierarchy, FullPopulate, SetSmoothers). These functions combine
% all factory needs as any specific needs given directly to the Hierarchy populator. 
%
% THUS, INDIVIDUAL FACTORY NEEDS ARE GLOBAL AND PASSED TO ALL FACTORIES INVOLVED
% HIERARCHY POPULATING. This is the main reason for not using the cross-factory
% mechansim for factory specific parameters.  For example, 
%
%          MyNeed.OutputLevel = 10
%          MyFactory.AddNeeds(MyNeeds)
%         
% actually set the output level to 10 in all factories during the build process.
%
% Key Concept 4
% =============
% CrossFactory is not stored permanently in the Hierarchy class. It must be 
% recreated from individual needs (which are stored) each time we populate.
% The main reason is that several populates might be invoked with different 
% factories to build hybrid methods (2 geometric MG levels followed by several
% multigrid levels).  Needs for these different factory sets may vary and so
% they should not reside permanently within the Hierarchy object. They are instead
% connected to specific build/populate.
%
% Note: the merged need list is returned and can be stored by a user.
%
%
% Two CrossFactory Functions
% ===============================
%
%  1) MergeNeeds(A,B) merges needs as best as it can figure out. For example,
%          A.SavePtent = '1:2';
%          B.SavePtent = '4:5'; 
%          C = MergeNeeds(A,B);
%     sets the value of C.SavePtent to '1:5'  (min(1,4) and max(2,5)).
%
%  2) TrueOrFalse(key, LevelId) checks if a key exists. If it is of the form
%     'x:y' it returns true only if x <= LevelId <= y. If it is not of the form
%     'x:y' (including 'all'), it always returns true. If the key does not
%     exist it returns false.

% Most cases are simple, though we highlight tricky stuff.
clear all;

SetHomeDir;
P=path; path([MUEMAT_ROOT_DIR,'/test/CrossFactorySpecTest'],P);

Amat        = BuildLaplace1DBlk(1,-1,81);
Finest      = Level();
Finest.Set('A', Amat);
Finest.Set('NullSpace', ones(Amat.GetRowMap().NDOFs(),1));
MgHierarchy = Hierarchy();
MgHierarchy.SetOutputLevel(1);
MgHierarchy.SetLevel(Finest,1);
fprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
fprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
fprintf('Test 1: Saving Tentative Prolongator for use in Restrictor Factory\n');
fprintf('        Rfact internally stores need for someone to save a tentative\n');
fprintf('        prolongator. FullPopulate() merges needs and passes a\n'); 
fprintf('        CrossFactory to build methods including the one in Pfact.\n');
fprintf('        Note: Verbose output not requested, so nothing appears\n');
AmalgamateDropFact= CoalesceDropFactory();
AggFact  = AggregationFactory();
Pfact       = SaPFactory(AmalgamateDropFact,AggFact);
Rfact       = TransPFactory();
Rfact.UsePtent();
PRfact = GenericPRFactory(Pfact,Rfact);
Acfact      = RAPFactory();
Sfact       = SmootherFactory(Smoother('GaussSeidel',1,1.));
NewMg = MgHierarchy.Copy();
NewMg.FullPopulate(PRfact, Acfact, Sfact, 1, 4);
fprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
fprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n\n');
%
%
%
fprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
fprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
fprintf('Test 2: User repeats Test 1 but requests verbose output\n');
fprintf('        It should be noted that the same effect is achieved (BUT HIGHLY\n');
fprintf('        DISCOURAGED) by removing Need from the arguments for FullPopulate\n');
fprintf('        and invoking Pfact.AddNeed(Need); just before FullPopulate().\n'); 
fprintf('        If one wants only Pfact to be verbose, use instead\n');
fprintf('             Pfact = Pfact.SetOutputLevel(10);\n');
Pfact.SetOutputLevel(10);
PRfact = GenericPRFactory(Pfact);
NewMg = MgHierarchy.Copy();
NewMg.FullPopulate(PRfact, Acfact, Sfact, 1, 4);
fprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
fprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n\n');
%
%
%
fprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
fprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
fprintf('Test 3: Aggregation factory requests graph be saved.\n');
fprintf('        SillyAggFact indicates that it needs a matrix graph. A\n'); 
fprintf('        A CrossFactory is built using Pfact.GetNeeds()\n');
fprintf('        which invokes AggFact.GetNeeds() so that these get\n');
fprintf('        passed up to the Hierarchy populator.\n');
AggFact  = SillyAggFact();
Pfact       = SaPFactory(AmalgamateDropFact,AggFact);
Rfact       = TransPFactory();
PRfact = GenericPRFactory(Pfact,Rfact);
NewMg = MgHierarchy.Copy();
NewMg.FullPopulate(PRfact, Acfact, Sfact, 1, 4);
fprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
fprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n\n');
%
%
%
fprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
fprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
fprintf('Test 4: Smoother requests aggregates be saved.\n');
AggFact  = AggregationFactory();
Pfact       = SaPFactory(AmalgamateDropFact,AggFact);
PRfact = GenericPRFactory(Pfact);
Sfact       = SillyRelaxFact();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n==> 1st attempt fails because smoother needs not given to FillHierarchy\n\n');
NewMg = MgHierarchy.Copy();
NewMg.FillHierarchy(PRfact, Acfact, 1, 4);
NewMg = NewMg.SetSmoothers(Sfact, 1, 4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n==> 2nd attempts almost succeeds. Smoother needs given to FullPopulate,\n');
fprintf('==> but fails on coarsest level as there are no aggregates there.\n\n');
NewMg = MgHierarchy.Copy();
NewMg.FullPopulate(PRfact, Acfact, Sfact, 1, 4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n==> 3rd attempt succeeds by turning off coarsest level smoother.\n\n');
NewMg = MgHierarchy.Copy();
NewMg.FullPopulate(PRfact, Acfact, Sfact, 1, 4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n==> 4th attempt succeeds by manually passing smoother needs into FillHierarchy.\n\n');
NewMg = MgHierarchy.Copy();
NewMg.FillHierarchy(PRfact, Acfact, 1, 4,...
                                    CrossFactory.MergeNeeds(Sfact.GetNeeds()));
NewMg.SetSmoothers(Sfact, 1, 4, ...
                                    CrossFactory.MergeNeeds(Sfact.GetNeeds()));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
fprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n\n');
%
%
%
fprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
fprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n');
fprintf('Test 5: User-requested reuse\n');
fprintf('        1st FillHierarchy() saves tentative prolongators\n');
fprintf('        while the 2nd FillHierarchy() uses them\n');
Sfact            = SmootherFactory(Smoother('GaussSeidel',1,1.));
Need.PersistentSavePtent = 'all';
Need.PersistentPtent = 'all';
Need.Ptent = 'all';
NewMg = MgHierarchy.Copy();
NewMg.FillHierarchy(PRfact, Acfact, 1, 4, Need);
Need.ReUsePtent = 'all';
NewMg.FillHierarchy(PRfact, Acfact, 1, 4, Need);

path(P);
