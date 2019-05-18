% Returns the element information
% def.nNod - number of nodes in the element
% def.nDof - number of degrees of freedom in the element
% def.Dof - Degrees of freedom correspond to the nodes 
% def.eltEdge - nodes of each element edge


function def=eltdef(typ,eltID,dofs)


typName = elttypname(typ,eltID);

switch lower(typName)
    case 'volu4'    
        def.nNod=4;
        def.nDof=4*dofs;  % 4 for 1 degree of freedom, 8 for 2 degrees of freedom, 12 for 3 degrees of freedom
        switch dofs
            case 1
                def.Dof=[1.03;2.03;3.03;4.03];% 1,2,3,4 represents the nodes and 01-x direction, 03-z direction
            case 2
                def.Dof=[1.01;1.03;2.01;2.03;3.01;3.03;4.01;4.03]; 
            case 3
                def.Dof = [1.01;1.02;1.03;2.01;2.02;2.03;3.01;3.02;3.03;4.01;4.02;4.03];
        end
        def.eltEdge=[1 2;2 3;3 4;4 1];
        def.eltSide=[1 2 3 4];
  
    case 'volu7'
        def.nNod=7;
        def.nDof=7*dofs;
        switch dofs
            case 1 
                def.Dof=[1.03;2.03;3.03;4.03;5.03;6.03;7.03];
            case 2
                def.Dof=[1.01;1.03;2.01;2.03;3.01;3.03;4.01;4.03;5.01;5.03;6.01;6.03;7.01;7.03];
            case 3
                def.Dof=[1.01;1.02;1.03;2.01;2.02;2.03;3.01;3.02;3.03;4.01;4.02;4.03;5.01;5.02;5.03;6.01;6.02;6.03;7.01;7.02;7.03];
        end

        
    case 'volu8'
        def.nNod=8;
        def.nDof=8*dofs;
        switch dofs
            case 1 
                def.Dof=[1.03;2.03;3.03;4.03;5.03;6.03;7.03;8.03];
            case 2
                def.Dof=[1.01;1.03;2.01;2.03;3.01;3.03;4.01;4.03;5.01;5.03;6.01;6.03;7.01;7.03;8.01;8.03];
            case 3
                def.Dof=[1.01;1.02;1.03;2.01;2.02;2.03;3.01;3.02;3.03;4.01;4.02;4.03;5.01;5.02;5.03;6.01;6.02;6.03;7.01;7.02;7.03;8.01;8.02;8.03];
        end 
        
    case 'beam1'     
         def.nNod=1;
         def.nDof=1*dofs;
         switch dofs
             case 1
                 def.Dof=[1.03];
             case 2
                def.Dof=[1.02;1.03];
             case 3
               def.Dof=[1.01;1.02;1.03];
         end
    case 'mass1'     
        def.nNod=1;
        def.nDof=1;
        def.Dof=[1.03];  
    case 'link2'   
        def.nNod=2;
        def.nDof=2;
        def.Dof=[1.03;2.03];
        def.eltEdge=[1 2];  
otherwise
  error('Unknown element type');
end

function typName = elttypname(typ,eltID)
nTyp=size(typ,1);
typName='';
for iTyp=1:nTyp
  if typ{iTyp,1}==eltID, typName=typ{iTyp,2}; end
end
