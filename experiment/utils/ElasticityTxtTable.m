% Charts for Elasticity Test

%fd = file descriptor

%TODO: - add an option to select printed columns (as in ElasticityTexTable)
%      - interface of ElasticityTexTable and ElasticityTxtTable should be
%        identical (file descriptor as output etc.)

function ElasticityTxtTable(dim, stretch, DATE, fd)

if ~varexist('fd'), fd=1; end % default = stdout

% Load File: (data.ITERS, LABEL, NLEVELS, NLIST, data.OC, RESID)
if(stretch==1), data = load(sprintf('etest%dd-filter3_304.%s.mat',dim,DATE));
else data = load(sprintf('etest%dd-filter3_304.s%d.%s.mat',dim,stretch,DATE));end

fprintf(fd, '\n');
fprintf(fd, '\n');
fprintf(fd, 'dim=%d, stretch=%d', dim, stretch);
fprintf(fd, '\n');

% fprintf(fd, '        |  SA-NR  |ML-SA-NR ||    SA   |  ML-SA  ||   EMIN  |\n');
% fprintf(fd, ' SZ LVL | ITS  OC | ITS  OC || ITS  OC | ITS  OC || ITS  OC |\n');
% fprintf(fd, '-------------------------------------------------------------\n');
% fprintf(fd, '%3d %2d  | %2d %4.2f | %2d %4.2f || %2d %4.2f | %2d %4.2f || %2d %4.2f |\n',[data.NLIST', data.NLEVELS(:,1)...
%     data.ITERS(:,1),data.OC(:,1),... % SA-NR
%     data.ITERS(:,5),data.OC(:,5),... % ML-SA-NR
%     data.ITERS(:,2),data.OC(:,2),... % SA
%     data.ITERS(:,4),data.OC(:,4),... % ML-SA
%     data.ITERS(:,3),data.OC(:,3)]'); % EMIN

% fprintf(fd, '        |  SA-NR  |ML-SA-NR ||    SA   |  ML-SA  ||   EMIN  || EMIN-NR | EMIN-SA || EMIN-NR2| EMIN-SA2|| SA-NF   | SA-NR-NF|| NOSMOO  |NOSMOO-NR|\n');
% fprintf(fd, ' SZ LVL | ITS  OC | ITS  OC || ITS  OC | ITS  OC || ITS  OC || ITS  OC | ITS  OC || ITS  OC | ITS  OC || ITS  OC | ITS  OC || ITS  OC | ITS  OC |\n');
% fprintf(fd, '-------------------------------------------------------------------------------------------------------------------------------------------------\n');
% fprintf(fd, '%3d %2d  | %2d %4.2f | %2d %4.2f || %2d %4.2f | %2d %4.2f || %2d %4.2f || %2d %4.2f | %2d %4.2f || %2d %4.2f | %2d %4.2f || %2d %4.2f | %2d %4.2f || %2d %4.2f | %2d %4.2f |\n',...
%         [data.NLIST',data.NLEVELS(:,1),...
%          data.ITERS(:,1),data.OC(:,1),... % SA-NR
%          data.ITERS(:,5),data.OC(:,5),... % ML-SA-NR
%          data.ITERS(:,2),data.OC(:,2),... % SA
%          data.ITERS(:,4),data.OC(:,4),... % ML-SA
%          data.ITERS(:,3),data.OC(:,3),... % EMIN
%          data.ITERS(:,7),data.OC(:,7),... % EMIN-SANR
%          data.ITERS(:,6),data.OC(:,6),...; % EMIN-SA
%          data.ITERS(:,9),data.OC(:,9),... % EMIN-SANR2
%          data.ITERS(:,8),data.OC(:,8),... % EMIN-SA2
%          data.ITERS(:,10),data.OC(:,10),...% SA-NF
%          data.ITERS(:,11),data.OC(:,11),... % SA-NR-NF
%          data.ITERS(:,13),data.OC(:,13),...% NoSmoo
%          data.ITERS(:,12),data.OC(:,12)]'); % NoSmoo-NR
     
fprintf(fd, '        |               No Rotation                 ||  ******  ||             With Rotation                 |\n');
fprintf(fd, '        |  NOSMOO  |    SA    |   SA-NF  |   EMIN   ||   EMIN   ||  NOSMOO  |    SA    |  SA-NF   |   EMIN   |\n');
fprintf(fd, ' SZ LVL | ITS  OC  | ITS  OC  | ITS  OC  | ITS  OC  || ITS  OC  || ITS  OC  | ITS  OC  | ITS  OC  | ITS  OC  |\n');
fprintf(fd, '--------------------------------------------------------------------------------------------------------------\n');
fprintf(fd, '%3d %2d  | %3d %4.2f | %3d %4.2f | %3d %4.2f | %3d %4.2f || %3d %4.2f || %3d %4.2f | %3d %4.2f | %3d %4.2f | %3d %4.2f |\n',...
        [data.NLIST', data.NLEVELS(:,1),...
         data.ITERS(:,12), data.OC(:,12),... % NoSmoo-NR 
         data.ITERS(:,1),  data.OC(:,1) ,... % SA-NR
         data.ITERS(:,11), data.OC(:,11),... % SA-NR-NF
         data.ITERS(:,7),  data.OC(:,7) ,... % EMIN-SANR
         data.ITERS(:,3),  data.OC(:,3) ,... % EMIN
         data.ITERS(:,13), data.OC(:,13),... % NoSmoo
         data.ITERS(:,2),  data.OC(:,2) ,... % SA
         data.ITERS(:,10), data.OC(:,10),... % SA-NF
         data.ITERS(:,6),  data.OC(:,6) ,... % EMIN-SA
         ]'); 



