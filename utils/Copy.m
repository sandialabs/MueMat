function newvar = Copy(var)

  if iscell(var) 
    if ndims(var) ~= 2, error('Copy.m: TODO 1'); end
    [n,m] = size(var);
    if n ~= 1, error('Copy.m: TODO 2'); end
   
    newvar = cell(size(var));
    for i=1:m
        newvar{i} = Copy(var{i});
    end
    
  elseif strcmp(class(var),'containers.Map')
   
    newvar = containers.Map();  % TODO: match KeyType/ValueType of var
    for k=keys(var)
      newvar(k{1}) = Copy(var(k{1})); % why k{1} ?
    end
       
  elseif isstruct(var)
      
    % also allows for vectors of struct objects...  
    newvar = [];
    for i=1:length(var)  % copy vector of length(var) structure objects
        newstruc = struct;
        names = fieldnames(var(i));
        for k=1:length(names)
            data = getfield(var(i),names{k});
            data2 = Copy(data);
           newstruc = setfield(newstruc,names{k},data2);
        end
        newvar = [newvar, newstruc];
    end
    
  else
    newvar = var;
    
  end
  
  %TODO if isvector() % attention! avoid recursion!
  
end