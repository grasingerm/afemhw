function [K]=tangent_stiffness(S,t,Ue,x,y,C)
%
%Function for obtaining the tangent stiffness matrix for a Q4 element in a
%finite strain finite element model
%Assume: Given x1, x2, x3, x4, y1, y2, y3, y4 (mapped nodal coords) for the
%element.  Also given {S} at each integration point, and {Ue}.  Also given
%t=element thickness.  Also given [C] such that {Sdot}=[C]{Edot}.
%
%Created 10/20/05 John Brigham
%Last Updated 10/20/05 JCB
%Comments:  Assumed all numbering (nodes and int points) starts at bottom
%left and counts counter-clockwise
%
%
%Define the integration points for a Q4 parent element
int_point=1/sqrt(3)*[-1, -1; 1, -1; 1, 1; -1, 1];
%initialize Kg and Km to zero
Kg=zeros(8,8);
Km=zeros(8,8);
%For all 4 integration points of a Q4
for i=1:4
    %define the coordinates of xi and eta from int point
    xi=int_point(i,1);
    eta=int_point(i,2);
    %Form the 4 shape functions N1, N2, N3, and N4
    N=1/4*[(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)];
    %Form derivatives of N1, N2, N3, and N4 w.r.t. parent coords
        %as in [dN/dxi], [dN/deta]
    dNdxi=1/4*[(eta-1), (1-eta), (1+eta), (-1-eta)];
    dNdeta=1/4*[(xi-1), (-1-xi), (1+xi), (1-xi)];
    %Form the jacobian mapping matrix [J]
        %[J]=[dx/dxi, dy/dxi; dx/deta, dy/deta]=[dN/dxi; dN/deta]*[x,y]
    J=[dNdxi;dNdeta]*[x,y];
    %Calculate [dNj/dx; dNj/dy] from [J]^-1*[dNj/dxi; dNj/deta] for
        %j=1,2,3,4
    Dummy=inv(J)*[dNdxi;dNdeta];
    dNdx=Dummy(1,:);
    dNdy=Dummy(2,:);
    %Form [Bl] from [dNi/dx] and [dNi/dy]
    Bl=[dNdx,0,0,0,0; 0,0,0,0,dNdy; dNdy,dNdx];
    %Form [G1] and [G2] from [dNi/dx] and [dNi/dy]
    G1=[dNdx,0,0,0,0;0,0,0,0,dNdx];
    G2=[dNdy,0,0,0,0;0,0,0,0,dNdy];
    G=[G1;G2];
    %Form {theta1} and {theta2} from {thetai}=[Gi]*{Ue}
    theta1=G1*Ue;
    theta2=G2*Ue;
    %Form [A] from {theta1} and {theta2}
    A=[theta1',0,0;0,0,theta2';theta2',theta1'];
    %Calculate [Bnl] from [A]*[G]
    Bnl=A*G;
    %Calculate [Bhat] from [Bl]+[Bnl]
    Bhat=Bl+Bnl;
    %Form [Ms] from {S} and [I2]
    Ms=[S(1,i),0,S(3,i),0;0,S(1,i),0,S(3,i);S(3,i),0,S(2,i),0;0,S(3,i),0,S(2,i)];
    %Calculate [Kgpoint]=[G]'*[Ms]*[G]
    Kgpoint=G'*Ms*G;
    %Calculate [Kmpoint]=[Bhat]'*[C]*[Bhat]
    Kmpoint=Bhat'*C*Bhat;
    %Calculate |[J]|=det[J]
    det_J=det(J);
    %Calculate [Kg] by the weighted sum of [Kgpoint]*t*|[J]| (for a Q4
        %element the weights are 1 for the 4 integration points)
    Kg=Kg+Kgpoint*t*det_J;
    %Calculate [Km] by the weighted sum of [Kmpoint]*t*|[J]| (for a Q4
        %element the weights are 1 for the 4 integration points)
    Km=Km+Kmpoint*t*det_J;
%end loop for each integration point
end
%
%Cacluate [K] from [Kg]+[Km]
K=Kg+Km;