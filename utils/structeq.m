%%% Test if two struct are equal.
%
% This function is not recursive and only test the first level of
% struct fields. Values must be string or valid for ~= operator.

% Warning: if struct1.blabla=default and blabla field don't exist,
% struct1 and 2 are not equal.

function [bool] = structeq(struct1, struct2)

% struct1.test='a';
% struct1.test2=2;

% struct2.test2=3;
% struct2.test='a';
% %struct2.test3='a';

s1 = orderfields(struct1);
s2 = orderfields(struct2);

f1 = fieldnames(s1);
f2 = fieldnames(s2);

bool = 0;

% Compare length
if length(f1) ~= length(f2), return; end

% Compare fields
for i=1:length(f1)
  if ~strcmp(f1{i}, f2{i}), return; end
end

% Compare fields values
for i=1:length(f1)
  a = getfield(s1,f1{i});
  b = getfield(s2,f2{i});
  
  if ischar(a) && ischar(b) && ~strcmp(a,b)
    return;
  elseif a~=b
    return;
  end
end  
  
bool = 1;

end