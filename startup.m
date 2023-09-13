function startup(logo)
  
  if ~exist('logo','var'), logo = true; end;
    
  opengl neverselect
  
  % Setup search path
  SetPath
  
  % Create 'SetHomeDir.m'
  fid = fopen('SetHomeDir.m','w');
  if fid > 0
    fprintf(fid,'%% This file is automatically generated by ''startup.m''.\n');
    fprintf(fid,'MUEMAT_ROOT_DIR=''%s'';\n',pwd);
    fclose(fid);
    clear fid
  end
  
  % logo
  if logo, mue; end
  
end
