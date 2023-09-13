% This script converts operators objects contained in .mat files to matlab matrices.

%cd old
filelist=dir('*.mat');

for i=1:size(filelist,1)
  filename = filelist(i).name;
  fprintf('%s: ',filename);
  
  load(filename);
  if isa(Amat, 'Operator')
    Amat=Amat.GetMatrixData();
    save(filename, 'Amat', 'nullspace', 'coords');
    fprintf('updated\n');
  else
    fprintf('nothing to do\n');
  end

end

%cd .. % old