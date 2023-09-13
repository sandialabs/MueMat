function Abar = TestNoCrossCrackDrop(A, Level, X)
% Function to drop entries in A depending on Crack-Crossing
% X contains xCr and COORD information
if ~isfield(X,'droptol'), X.droptol = 1e-3; end
    fprintf('Predropping crack-cross...');
    Pattern = DontCrossCrackMat(X,Level.Get('Arr').GetMatrixData,Level.Get('dof2node'));
    Abar = spones(Pattern).*SAdrop(Level.Get('Arr').GetMatrixData(), Level,X.droptol);

    fprintf('Done!\n')
    % keyboard
end



function Pattern = DontCrossCrackMat(X,A,dof2node)
% X needs the following properties/methods
% xCr = {ncr x 1}(2x2)
% COORD = [ndof x 2] where ndof is size of A(ndof x ndof)
numrows = size(A,1);
Pattern = sparse(A~=0);
for i = 1:numrows
    dofcouple = setdiff(find(A(:,i)~=0), i);
    nodei = X.COORD(dof2node(i),:);
    nodecouple = X.COORD(dof2node(dofcouple),:);
    ncr = length(X.xCr);
    for k = 1:ncr
        lsetki = MyLevelSets(nodei,X.xCr,k);
        nlseti = lsetki(1);
        tlseti = lsetki(2);
        
        lsetkcouple = MyLevelSets(nodecouple,X.xCr,k);
        nlsetcouple = lsetkcouple(:,1);
        tlsetcouple = lsetkcouple(:,2);
        
        if (tlseti <= 0)  % tangent crossing
           ncrossers = dofcouple(find(nlsetcouple * nlseti < 0)); %normal cross
           Pattern(i,ncrossers) = 0;
           Pattern(ncrossers,i) = 0;
        end;
    end
    
end
end


function plset = MyLevelSets(pt,xCr,k)
% Function to generate levelsets at given point for given crack
numx = size(pt,1); 
plset = zeros(numx,2);

% The point coordinates
x = pt(:,1); 
y = pt(:,2);       
% Crack tips
x0 = xCr{k}(1,1); y0 = xCr{k}(1,2); % tip1
x1 = xCr{k}(2,1); y1 = xCr{k}(2,2); % tip2
% Crack segment, normalized and length
seg = xCr{k}(2,:)-xCr{k}(1,:);
t = 1/norm(seg)*seg;
l = sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)) ;

% Normal Levelset
phi = (y0-y1)*x + (x1-x0)*y + (x0*y1-x1*y0);
plset(:,1) = phi/l;                                    % normal LS

% Tangent Levelset
tls1 = ([x-xCr{k}(1,1) y-xCr{k}(1,2)])*t';
tls2 = ([x-xCr{k}(2,1) y-xCr{k}(2,2)])*t';
C1 = find(abs(tls1)<abs(tls2));
C2 = find(abs(tls1)>=abs(tls2));
plset(C1,2) = -tls1(C1);
plset(C2,2) =  tls2(C2);

end
