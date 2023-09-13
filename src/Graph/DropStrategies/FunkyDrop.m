%
% New funky dropping strategy. The basic idea is that we have two
% criteria:
%
%        1) percentage of largest off-diag : roughly, drop stuff below a
%                                            certain percentage of the 
%                                            largest nonzero offdiagonal
%                                            in row.
% and 
%
%        2) percentage of off-diagonal sum : roughly, don't drop stuff
%                                            that would lower sum(abs(row))row))
%                                            lower than a certain percentage.
%
% The details are kind of nasty, but I put more comments in the code below
%
function [matrix] = FunkyDrop(matrix, MyLevel, widget)

n = size(matrix,1);
%
% Grab parameters that define the dropping specifications
%
SizePercent     = widget.SizePercent;
SizeDelta       = widget.SizeDelta;
LowerSumPercent = widget.LowerSumPercent;
SizeWeight      = widget.SizeWeight;
SumWeight       = widget.SumWeight;
%
% Symmetrically scale matrix, then strip diag, and take absolute values
%
original = matrix;
invDiag  = sparse(1:n,1:n,abs(1./sqrt(abs(diag(matrix)))));
matrix   = invDiag'*matrix*invDiag;     
matrix   = abs(matrix - spdiags(diag(matrix),0,n,n));
%
% Grab nonzeros and sort so that rows are consecutive. Additionally,
% sort within rows so that large entries appear first.
%
[Rows,Cols,Values] = find(matrix);
nnn = length(Rows);
[dummy,iii] = sort(Rows);
Rows    = Rows(iii);
Cols    = Cols(iii);
Values  = Values(iii);
RowStart= zeros(n+1,1);   % RowStart(i) gives the start location of row i
                           % within the arrays Rows,Cols,Values.
Start      = 1;
RowStart(1)= Start;
Count      = 1;
Rows = [Rows ; -7]; % Append crud so code below works for last row.
for Row= 1:n+1
   flag = 1;
   if Count > nnn, flag = 0; end;
   while flag == 1,
      if Rows(Count) ~= Row,
         Finish = Count-1;
         [newValues, iii] = sort(-abs(Values(Start:Finish)));
         temp = Cols(Start:Finish);
         Cols(Start:Finish) = temp(iii);
         temp = Values(Start:Finish);
         Values(Start:Finish) = temp(iii);
         Start = Count;
         RowStart(Row+1) = Start;
         flag  = 0;
      end
      Count = Count + 1;
   end
end
Rows  = Rows(1:end-1); % Remove crud from end of Rows
NumNzs= length(Rows);

% Now compute a bunch of row-wise quantities. Some of these quantities
% are stored in an n-vector while others are stored in a vector of
% length NumNzs. That is, even though it is a row-quantity, we store
% repeated entries of the row corresponding to each nonzero within that row.
% This makes some Matlab stuff a bit more efficient. To do this we 
% make use of RowToNnzs

RowToNnzs      = sparse( (1:NumNzs), Rows, ones(NumNzs,1), NumNzs, n);
MaxValuePerRow = max(matrix')';
RowSums        = matrix*ones(n,1);
MyMaxValue     = RowToNnzs*MaxValuePerRow;
TotalSum       = RowToNnzs*RowSums;
%
%  Within each row we look for a place to cut-off (between values kept and 
%  values tossed).  The actual cut-off spot will be based on three scores: 
%  SizeScore, Gradient, and SumScore. When determining SizeScore, 
%  focus on those whose values are between MinSizes and MaxSizes where
%  these are based on the lower percentage of SizePercent-SizeDelta
%  and the upper percentage of SizePercent+SizeDelta. The actual value of
%  SizeScore favors those values that are right at MyMaxValue*SizePercent
%  and trails off as we approach MinSizes and MaxSizes. SizeWeights determines
%  how steep this trail-off is. When determining SumScore, we use MinSums
%  which is the percentage of TotalSum = sum(abs(row)) in conjunction 
%  with the following formula
%
%    SumScore  = 1 - SumWeight*(TotalSum - MyNnzSum)./(TotalSum - MinSums);
%
%  where MyNnzSum = sum(abs(myvalue + row values bigger than me)). Thus,
%  when MyNnzSum = MinSums, the fraction is just equal to SumWeight and
%  when MyNnzSum = 0 then it is just SumWeight*TotalSum/(TotalSum-MinSums)

SizeScore  = 1 - SizeWeight*abs((Values./MyMaxValue) - SizePercent)/SizeDelta;
MinSizes   = (SizePercent - SizeDelta)*MaxValuePerRow;
MaxSizes   = (SizePercent + SizeDelta)*MaxValuePerRow;
SizeScore(abs(Values) < RowToNnzs*MinSizes) = 0.;
SizeScore(abs(Values) > RowToNnzs*MaxSizes) = 0.;
MinSums   = LowerSumPercent*TotalSum;
MyNnzSum = zeros(NumNzs,1);
for i=1:n
     MySum = 0;
     for j=RowStart(i):RowStart(i+1)-1
        MySum = MySum + Values(j);
        MyNnzSum(j) = MySum;
     end
end
SumScore  = 1 - SumWeight*(TotalSum - MyNnzSum)./(TotalSum - MinSums);
%
% Any sum below MinSums is not considered as a possible cut-off
%
ZeroOut = find((MyNnzSum < MinSums) .* (SumWeight ~= 0));
SumScore(ZeroOut) = 0;

% compute some kind of difference to measure how steep values within a
% row are decreasing. Use this as well in the weighting.
%
Gradient      = [ 0 ; Values(1:end-1)-Values(2:end)];
CombinedScore = Gradient .* SumScore .* SizeScore;

% Convert the CombinedScore into a CutOff ... essentially taking the
% largest CombinedScore within a row.
%
CutOff = RowStart(2:end);
for i=1:n
   [MaxScore,ii] = max(CombinedScore(RowStart(i)+1:RowStart(i+1)-1));
   if MaxScore ~= 0, 
        CutOff(i) = RowStart(i)+ii;
   else  
        % if all the CombinedScores are 0, use the MinSizes as the 
        % determining factor for the cut-off
        for j=RowStart(i+1)-1:-1:RowStart(i)+1
           if Values(j) > MinSizes(i), CutOff(i) = j+1; break; end
        end
   end
end
inds  = find( (1:NumNzs)' < (RowToNnzs*CutOff) );
matrix = sparse(Rows(inds),Cols(inds),Values(inds),n,n);
boo = matrix*ones(n,1);
aaa = find(boo == 0);
if ~isempty(aaa),
   fprintf('There should not be any empties\n');
   keyboard;
end
matrix = original .* (spones(abs(matrix) + abs(matrix'))+speye(n,n));
