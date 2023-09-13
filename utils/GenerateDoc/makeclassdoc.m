SetHomeDir
cd(MUEMAT_ROOT_DIR);
%!utils/makeclassdoc.sh
!utils/checkContents-m.sh

%cd([MUEMAT_ROOT_DIR '/doc/source']);
opts.outputDir = [MUEMAT_ROOT_DIR '/doc/html'];

filelist = dir([MUEMAT_ROOT_DIR '/src/']);

for i=3:size(filelist,1) % loop skips . and ..
  if (filelist(i).isdir)
    %makedoc(['ref_' filelist(i).name '.m']);

    ContentsM = [MUEMAT_ROOT_DIR '/src/' filelist(i).name '/Contents.m'];
    if exist(ContentsM, 'file')
      makedoc2(ContentsM, opts, ['ref_' filelist(i).name '.html']);
    end
  end
end

% more...
makedoc2([MUEMAT_ROOT_DIR '/example/Contents.m'], opts, ['Tutorial.html']);