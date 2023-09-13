% The IJVCellMatrix class allows to store IJV matrices by parts.
%
% IJVCellMatrix(nBlk)  : nBlk = number of parts
% M.Add(iBlk, I, J, V) : replace the part iBlkm by [I,J,V]
% M.GetIJV()           : get the matrix in IJV format.
%
% Example :
%
% M = [I,J,V] with I=[1 2 3 4 5 6], J=[1 1 3 3 4 5]; V=[6 5 4 3 2 1];
%
% I1 = [1 2 3 4]; I2 = [5 6]; 
% J1 = [1 1 3 3]; J2 = [4 5];
% V1 = [6 5 4 3]; V2 = [2 1];
%
% myMatrix = IJVCellMatrix(2);
% myMatrix = myMatrix.Add(1, I1,J1,V1);
% myMatrix = myMatrix.Add(2, I2,J2,V2);
% [I,J,V]  = myMatrix.GetIJV();
%
% TODO: in which directory I can put this file ?
%
classdef IJVCellMatrix < CopiableHandle
   properties (SetAccess = private)
      Cell_
   end
   methods
      function this = IJVCellMatrix(nBlk)
         % Copy constructor
         if nargin == 1 && isa(nBlk, class(this)), this.Copy_(nBlk,[]); return; end
         % 
      
         this.Cell_ = cell(3,nBlk); % 3 = I,J,V
      end
      function Add(this, iBlk, I, J, V)
         this.Cell_{1,iBlk} = I;
        this.Cell_{2,iBlk} = J;
        this.Cell_{3,iBlk} = V;
      end
      function [II, JJ, VV] = GetIJV(this) % concat I,J,V subvectors
         dim = 1; if (size(this.Cell_{1,1},1) == 1), dim=2; end % why ?

         II = cat(dim,this.Cell_{1,:}); 
         JJ = cat(dim,this.Cell_{2,:});
         VV = cat(dim,this.Cell_{3,:}); 
      end
   end % methods
   
   methods (Access = protected)  

    function Copy_(this, src, mc)
      %COPY_ 
      % 
      %   SYNTAX   obj.Copy_(src, mc);
      % 
      %     src - Object to copy
      %     mc  - MATLAB Metaclass
      [cmd, data, mc] = this.CopyCmd_(src,mc);
      eval(cmd);    
    end
    
  end % methods
end % classdef