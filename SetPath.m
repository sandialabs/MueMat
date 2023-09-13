function SetPath()
  %SETPATH Setup search path for MueMat
  %
  %   SYNTAX   SetPath();
  
  addFolders({   'src',
                 'gallery',
                 'example',
                 'test',
                 'experiment', 
                 'utils', 
                 %'thirdparty',
             });

  addpath(fullfile(pwd,'doc'));
  addpath(pwd);
           
end

function addFolders(dirList)
  %SETPATH Add folders and subfolders to MATLAB search path
  %
  %   SYNTAX   addFolders(dirList);
  %
  %   EXAMPLE
  %
  %     addFolders({'src', 'example'}); path;
  %
  for i=size(dirList):-1:1 % in reverse order
    addpath(genpath(fullfile(pwd, dirList{i})));
  end
end
