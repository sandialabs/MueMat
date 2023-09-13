function vizAggsFromP(P, Pt, Ncpn1, Ncpn2, coords)
  % Graph aggregates from the prolongator (2D)
  %
  % input:
  %   P       -- Prolongator (matlab matrix)
  %   Pt      -- Tentative Prolongator (matlab matrix)
  %   Ncpn1   -- DOF in the fine grid
  %   Ncpn2   -- DOF in the coarse grid
  %   coords  -- Location of the nodes
  % output:
  %   plot of the aggregates, one by one
  %
  % Example usage:
  %  P =MgHierarchy.Levels_{2}.Get('P');
  %  Pt=MgHierarchy.Levels_{2}.Get('Ptent');
  %  % Elasticitity2D, SA: Ncpn1=2; Ncpn2=3;
  %  vizAggsFromP(P,GetMatrixData(), Pt.GetMatrixData(), 2, 3, coords);
  
  % Create nodal P/Ptent
  Pnode=sparse(size(P,1)/Ncpn1,size(P,2)/Ncpn2);
  Ptnode=sparse(size(Pt,1)/Ncpn1,size(Pt,2)/Ncpn2);
  
  for ii=1:Ncpn1, 
    for jj=1:Ncpn2, 
      Pnode  = Pnode  + (abs(P (ii:Ncpn1:end,jj:Ncpn2:end))> 0);
      Ptnode = Ptnode + (abs(Pt(ii:Ncpn1:end,jj:Ncpn2:end))> 0);
    end
  end
  
  % Create agg list from P
  for I=1:size(Pnode,2),
    agg{I}=find(Pnode(:,I));
    tagg{I}=find(Ptnode(:,I));
  end
  
  %plot(coords(:,1),coords(:,2),'rx');
  for IDX=1:size(Pnode,2),
    plot(coords(:,1),coords(:,2),'r.',coords(tagg{IDX},1),coords(tagg{IDX},2),'bx',coords(agg{IDX},1),coords(agg{IDX},2),'gs');
    title(sprintf('aggregate %d',IDX));
    pause;
  end

end