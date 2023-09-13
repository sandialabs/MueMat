function [AggId, roots] = SquareAggs(Amat, options)

A     = Amat.GetMatrixData();
dim   = options.Dim;

if dim == 1,
fprintf('SquareAggs currently only works for 2D problems\n');
keyboard;
   nx    = options.Nx;
   Diamx = options.Diamx;
   if mod(nx,Diamx) ~= 0,
      error('MakeSquareAggs: nx must be divisible by options.Diamx');
   end
   AggId = ceil( (1:nx)/Diamx - .001);

elseif dim == 2,
   nx  = options.Nx;
   ny  = options.Ny;
   Diamx = options.Diamx;
   Diamy = options.Diamy;
   if mod(nx,Diamx) ~= 0,
      error('MakeSquareAggs: nx must be divisible by options.Diamx');
   end
   if mod(ny,Diamy) ~= 0,
      error('MakeSquareAggs: ny must be divisible by options.Diamy');
   end
   tx = ceil( (1:nx)/Diamx - .001);
   ty = ceil( (1:ny)/Diamy - .001);

   k = 1; roots = [];
   for j=1:length(ty)
      for i=1:length(tx)
         rotx = Diamx*ceil(i/Diamx) - floor(Diamx/2);
         roty = Diamy*ceil(j/Diamy) - floor(Diamy/2);
         ttt = (roty-1)*nx + rotx;
         if (k == 1) roots(k) = ttt; k = k+1;
         elseif (ttt > roots(k-1)), roots(k) = ttt; k=k+1; end
      end
   end
   temp = (ty-1)*max(tx);
   AggId = reshape(ones(length(tx),1)*temp+tx'*ones(1,length(ty)),nx*ny,1);
else
fprintf('SquareAggs currently only works for 2D problems\n');
keyboard;
   nx  = options.Nx;
   ny  = options.Ny;
   nz  = options.Nz;
   Diamx = options.Diamx;
   Diamy = options.Diamy;
   Diamz = options.Diamz;
   if mod(nx,Diamx) ~= 0,
      error('MakeSquareAggs: nx must be divisible by options.Diamx');
   end
   if mod(ny,Diamy) ~= 0,
      error('MakeSquareAggs: ny must be divisible by options.Diamy');
   end
   if mod(nz,Diamz) ~= 0,
      error('MakeSquareAggs: nz must be divisible by options.Diamz');
   end
   tx = ceil( (1:nx)/Diamx - .001);
   ty = ceil( (1:ny)/Diamy - .001);
   tz = ceil( (1:nz)/Diamz - .001);
   AggId = zeros(nx*ny*nz,1);
   MaxTx = max(tx);  MaxTy = max(ty);
   for i=1:nx, for j=1:ny, for k=1:nz
      AggId(i+(j-1)*nx+(k-1)*nx*ny) =  ...
             (tz(k)-1)*MaxTx*MaxTy + (ty(j)-1)*MaxTx + tx(i);
           end; end; end
end

%ndofs = Amat.GetRowMap().NDOFs();
%NodesInAgg = sparse(max(AggId), ndofs);
%for i=1:ndofs, NodesInAgg(AggId(i),i) = 1; end

