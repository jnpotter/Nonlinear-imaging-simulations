function [Fx Fy]=non_linear_forcing_computation1(linear_ux, linear_uy)

global mu rho K dx dy n m A B C gpu_on

global duxdx duxdy duydx duydy dduxdxdx dduxdxdy dduxdydy dduydydy dduydxdy dduydxdx
% 
 ux=linear_ux;uy=linear_uy;
 
 
%  %1st order derivatives (6th order accuracy)
% duxdx(:,4:end-3)=(-1/60*ux(:,(4-3):(end-3-3)) + 3/20*ux(:,(4-2):(end-3-2)) - 3/4*ux(:,(4-1):(end-3-1)) +  3/4*ux(:,(4+1):(end-3+1)) - 3/20*ux(:,(4+2):(end-3+2)) +1/60*ux(:,(4+3):(end-3+3)))/dx; 
% duydy(4:end-3,:)=(-1/60*uy((4-3):(end-3-3),:) + 3/20*uy((4-2):(end-3-2),:) - 3/4*uy((4-1):(end-3-1),:) +  3/4*uy((4+1):(end-3+1),:) - 3/20*uy((4+2):(end-3+2),:) +1/60*uy((4+3):(end-3+3),:))/dy; 
% duxdy(4:end-3,:)=(-1/60*ux((4-3):(end-3-3),:) + 3/20*ux((4-2):(end-3-2),:) - 3/4*ux((4-1):(end-3-1),:) +  3/4*ux((4+1):(end-3+1),:) - 3/20*ux((4+2):(end-3+2),:) +1/60*ux((4+3):(end-3+3),:))/dy; 
% duydx(:,4:end-3)=(-1/60*uy(:,(4-3):(end-3-3)) + 3/20*uy(:,(4-2):(end-3-2)) - 3/4*uy(:,(4-1):(end-3-1)) +  3/4*uy(:,(4+1):(end-3+1)) - 3/20*uy(:,(4+2):(end-3+2)) +1/60*uy(:,(4+3):(end-3+3)))/dx; 
% 
% %2nd order derivatives (6th order accuracy)
% dduxdxdx(:,4:end-3)=(1/90*ux(:,(4-3):(end-3-3)) - 3/20*ux(:,(4-2):(end-3-2)) + 3/2*ux(:,(4-1):(end-3-1)) -49/18*ux(:,4:end-3) +  3/2*ux(:,(4+1):(end-3+1)) - 3/20*ux(:,(4+2):(end-3+2)) +1/90*ux(:,(4+3):(end-3+3)))/dx; 
% dduxdydy(4:end-3,:)=(1/90*ux((4-3):(end-3-3),:) - 3/20*ux((4-2):(end-3-2),:) + 3/2*ux((4-1):(end-3-1),:) -49/18*ux(4:end-3,:) +  3/2*ux((4+1):(end-3+1),:) - 3/20*ux((4+2):(end-3+2),:) +1/90*ux((4+3):(end-3+3),:))/dy; 
% dduydydy(4:end-3,:)=(1/90*uy((4-3):(end-3-3),:) - 3/20*uy((4-2):(end-3-2),:) + 3/2*uy((4-1):(end-3-1),:) -49/18*uy(4:end-3,:) +  3/2*uy((4+1):(end-3+1),:) - 3/20*uy((4+2):(end-3+2),:) +1/90*uy((4+3):(end-3+3),:))/dy; 
% dduydxdx(:,4:end-3)=(1/90*uy(:,(4-3):(end-3-3)) - 3/20*uy(:,(4-2):(end-3-2)) + 3/2*uy(:,(4-1):(end-3-1)) -49/18*uy(:,4:end-3) +  3/2*uy(:,(4+1):(end-3+1)) - 3/20*uy(:,(4+2):(end-3+2)) +1/90*uy(:,(4+3):(end-3+3)))/dx; 
% 
% 
% dduxdxdy(4:end-3,4:end-3)=...
% +(-1/60)*(-1/60*ux((4-3):(end-3-3),(4-3):(end-3-3)) + 3/20*ux((4-3):(end-3-3),(4-2):(end-3-2)) - 3/4*ux((4-3):(end-3-3),(4-1):(end-3-1)) +  3/4*ux((4-3):(end-3-3),(4+1):(end-3+1)) - 3/20*ux((4-3):(end-3-3),(4+2):(end-3+2)) +1/60*ux((4-3):(end-3-3),(4+3):(end-3+3)))/dx/dy...
% +(3/20)*(-1/60*ux((4-2):(end-3-2),(4-3):(end-3-3)) + 3/20*ux((4-2):(end-3-2),(4-2):(end-3-2)) - 3/4*ux((4-2):(end-3-2),(4-1):(end-3-1)) +  3/4*ux((4-2):(end-3-2),(4+1):(end-3+1)) - 3/20*ux((4-2):(end-3-2),(4+2):(end-3+2)) +1/60*ux((4-2):(end-3-2),(4+3):(end-3+3)))/dx/dy... 
% +(-3/4)*(-1/60*ux((4-1):(end-3-1),(4-3):(end-3-3)) + 3/20*ux((4-1):(end-3-1),(4-2):(end-3-2)) - 3/4*ux((4-1):(end-3-1),(4-1):(end-3-1)) +  3/4*ux((4-1):(end-3-1),(4+1):(end-3+1)) - 3/20*ux((4-1):(end-3-1),(4+2):(end-3+2)) +1/60*ux((4-1):(end-3-1),(4+3):(end-3+3)))/dx/dy...
% +(3/4)*(-1/60*ux((4+1):(end-3+1),(4-3):(end-3-3)) + 3/20*ux((4+1):(end-3+1),(4-2):(end-3-2)) - 3/4*ux((4+1):(end-3+1),(4-1):(end-3-1)) +  3/4*ux((4+1):(end-3+1),(4+1):(end-3+1)) - 3/20*ux((4+1):(end-3+1),(4+2):(end-3+2)) +1/60*ux((4+1):(end-3+1),(4+3):(end-3+3)))/dx/dy...
% +(-3/20)*(-1/60*ux((4+2):(end-3+2),(4-3):(end-3-3)) + 3/20*ux((4+2):(end-3+2),(4-2):(end-3-2)) - 3/4*ux((4+2):(end-3+2),(4-1):(end-3-1)) +  3/4*ux((4+2):(end-3+2),(4+1):(end-3+1)) - 3/20*ux((4+2):(end-3+2),(4+2):(end-3+2)) +1/60*ux((4+2):(end-3+2),(4+3):(end-3+3)))/dx/dy...
% +(1/60)*(-1/60*ux((4+3):(end-3+3),(4-3):(end-3-3)) + 3/20*ux((4+3):(end-3+3),(4-2):(end-3-2)) - 3/4*ux((4+3):(end-3+3),(4-1):(end-3-1)) +  3/4*ux((4+3):(end-3+3),(4+1):(end-3+1)) - 3/20*ux((4+3):(end-3+3),(4+2):(end-3+2)) +1/60*ux((4+3):(end-3+3),(4+3):(end-3+3)))/dx/dy; 
% 
% 
% dduydxdy(4:end-3,4:end-3)=...
% +(-1/60)*(-1/60*uy((4-3):(end-3-3),(4-3):(end-3-3)) + 3/20*uy((4-3):(end-3-3),(4-2):(end-3-2)) - 3/4*uy((4-3):(end-3-3),(4-1):(end-3-1)) +  3/4*uy((4-3):(end-3-3),(4+1):(end-3+1)) - 3/20*uy((4-3):(end-3-3),(4+2):(end-3+2)) +1/60*uy((4-3):(end-3-3),(4+3):(end-3+3)))/dx/dy...
% +(3/20)*(-1/60*uy((4-2):(end-3-2),(4-3):(end-3-3)) + 3/20*uy((4-2):(end-3-2),(4-2):(end-3-2)) - 3/4*uy((4-2):(end-3-2),(4-1):(end-3-1)) +  3/4*uy((4-2):(end-3-2),(4+1):(end-3+1)) - 3/20*uy((4-2):(end-3-2),(4+2):(end-3+2)) +1/60*uy((4-2):(end-3-2),(4+3):(end-3+3)))/dx/dy... 
% +(-3/4)*(-1/60*uy((4-1):(end-3-1),(4-3):(end-3-3)) + 3/20*uy((4-1):(end-3-1),(4-2):(end-3-2)) - 3/4*uy((4-1):(end-3-1),(4-1):(end-3-1)) +  3/4*uy((4-1):(end-3-1),(4+1):(end-3+1)) - 3/20*uy((4-1):(end-3-1),(4+2):(end-3+2)) +1/60*uy((4-1):(end-3-1),(4+3):(end-3+3)))/dx/dy...
% +(3/4)*(-1/60*uy((4+1):(end-3+1),(4-3):(end-3-3)) + 3/20*uy((4+1):(end-3+1),(4-2):(end-3-2)) - 3/4*uy((4+1):(end-3+1),(4-1):(end-3-1)) +  3/4*uy((4+1):(end-3+1),(4+1):(end-3+1)) - 3/20*uy((4+1):(end-3+1),(4+2):(end-3+2)) +1/60*uy((4+1):(end-3+1),(4+3):(end-3+3)))/dx/dy...
% +(-3/20)*(-1/60*uy((4+2):(end-3+2),(4-3):(end-3-3)) + 3/20*uy((4+2):(end-3+2),(4-2):(end-3-2)) - 3/4*uy((4+2):(end-3+2),(4-1):(end-3-1)) +  3/4*uy((4+2):(end-3+2),(4+1):(end-3+1)) - 3/20*uy((4+2):(end-3+2),(4+2):(end-3+2)) +1/60*uy((4+2):(end-3+2),(4+3):(end-3+3)))/dx/dy...
% +(1/60)*(-1/60*uy((4+3):(end-3+3),(4-3):(end-3-3)) + 3/20*uy((4+3):(end-3+3),(4-2):(end-3-2)) - 3/4*uy((4+3):(end-3+3),(4-1):(end-3-1)) +  3/4*uy((4+3):(end-3+3),(4+1):(end-3+1)) - 3/20*uy((4+3):(end-3+3),(4+2):(end-3+2)) +1/60*uy((4+3):(end-3+3),(4+3):(end-3+3)))/dx/dy; 

% % 
%1st order derivatives (4th order accuracy)
duxdx(:,3:end-2)=(1/12*ux(:,(3-2):(end-2-2)) - 2/3*ux(:,(3-1):(end-2-1)) +  2/3*ux(:,(3+1):(end-2+1)) - 1/12*ux(:,(3+2):(end-2+2)))/dx; 
duydy(3:end-2,:)=(1/12*uy((3-2):(end-2-2),:) - 2/3*uy((3-1):(end-2-1),:) +  2/3*uy((3+1):(end-2+1),:) - 1/12*uy((3+2):(end-2+2),:))/dy;
duxdy(3:end-2,:)=(1/12*ux((3-2):(end-2-2),:) - 2/3*ux((3-1):(end-2-1),:) +  2/3*ux((3+1):(end-2+1),:) - 1/12*ux((3+2):(end-2+2),:))/dy;
duydx(:,3:end-2)=(1/12*uy(:,(3-2):(end-2-2)) - 2/3*uy(:,(3-1):(end-2-1)) +  2/3*uy(:,(3+1):(end-2+1)) - 1/12*uy(:,(3+2):(end-2+2)))/dx;

%2nd order derivatives (4th order accuracy)

dduxdxdx(:,3:end-2)=(-1/12*ux(:,(3-2):(end-2-2)) + 4/3*ux(:,(3-1):(end-2-1)) - 5/2*ux(:,3:end-2) + 4/3*ux(:,(3+1):(end-2+1)) - 1/12*ux(:,(3+2):(end-2+2)))/dx^2;
dduxdydy(3:end-2,:)=(-1/12*ux((3-2):(end-2-2),:) + 4/3*ux((3-1):(end-2-1),:) - 5/2*ux(3:end-2,:) + 4/3*ux((3+1):(end-2+1),:) - 1/12*ux((3+2):(end-2+2),:))/dy^2;
dduydydy(3:end-2,:)=(-1/12*uy((3-2):(end-2-2),:) + 4/3*uy((3-1):(end-2-1),:) - 5/2*uy(3:end-2,:) + 4/3*uy((3+1):(end-2+1),:) - 1/12*uy((3+2):(end-2+2),:))/dy^2;
dduydxdx(:,3:end-2)=(-1/12*uy(:,(3-2):(end-2-2)) + 4/3*uy(:,(3-1):(end-2-1)) - 5/2*uy(:,3:end-2) + 4/3*uy(:,(3+1):(end-2+1)) - 1/12*uy(:,(3+2):(end-2+2)))/dx^2;

dduxdxdy(3:end-2,3:end-2)=...
+(1/12)*(1/12*ux((3-2):(end-2-2),(3-2):(end-2-2)) - 2/3*ux((3-2):(end-2-2),(3-1):(end-2-1)) + 2/3*ux((3-2):(end-2-2),(3+1):(end-2+1)) - 1/12*ux((3-2):(end-2-2),(3+2):(end-2+2)))/dx/dy...
+(-2/3)*(1/12*ux((3-1):(end-2-1),(3-2):(end-2-2)) - 2/3*ux((3-1):(end-2-1),(3-1):(end-2-1)) + 2/3*ux((3-1):(end-2-1),(3+1):(end-2+1)) - 1/12*ux((3-1):(end-2-1),(3+2):(end-2+2)))/dx/dy...
+(2/3)*(1/12*ux((3+1):(end-2+1),(3-2):(end-2-2)) - 2/3*ux((3+1):(end-2+1),(3-1):(end-2-1)) + 2/3*ux((3+1):(end-2+1),(3+1):(end-2+1)) - 1/12*ux((3+1):(end-2+1),(3+2):(end-2+2)))/dx/dy...
+(-1/12)*(1/12*ux((3+2):(end-2+2),(3-2):(end-2-2)) - 2/3*ux((3+2):(end-2+2),(3-1):(end-2-1)) + 2/3*ux((3+2):(end-2+2),(3+1):(end-2+1)) - 1/12*ux((3+2):(end-2+2),(3+2):(end-2+2)))/dx/dy;

dduydxdy(3:end-2,3:end-2)=...
+(1/12)*(1/12*uy((3-2):(end-2-2),(3-2):(end-2-2)) - 2/3*uy((3-2):(end-2-2),(3-1):(end-2-1)) + 2/3*uy((3-2):(end-2-2),(3+1):(end-2+1)) - 1/12*uy((3-2):(end-2-2),(3+2):(end-2+2)))/dx/dy...
+(-2/3)*(1/12*uy((3-1):(end-2-1),(3-2):(end-2-2)) - 2/3*uy((3-1):(end-2-1),(3-1):(end-2-1)) + 2/3*uy((3-1):(end-2-1),(3+1):(end-2+1)) - 1/12*uy((3-1):(end-2-1),(3+2):(end-2+2)))/dx/dy...
+(2/3)*(1/12*uy((3+1):(end-2+1),(3-2):(end-2-2)) - 2/3*uy((3+1):(end-2+1),(3-1):(end-2-1)) + 2/3*uy((3+1):(end-2+1),(3+1):(end-2+1)) - 1/12*uy((3+1):(end-2+1),(3+2):(end-2+2)))/dx/dy...
+(-1/12)*(1/12*uy((3+2):(end-2+2),(3-2):(end-2-2)) - 2/3*uy((3+2):(end-2+2),(3-1):(end-2-1)) + 2/3*uy((3+2):(end-2+2),(3+1):(end-2+1)) - 1/12*uy((3+2):(end-2+2),(3+2):(end-2+2)))/dx/dy;

% %first spatial derivatives
% 
% duxdx(:,2:end-1)=(linear_ux(:,3:end)-linear_ux(:,1:end-2))/dx/2;
% duxdx(:,1)=duxdx(:,2);
% duxdx(:,end)=duxdx(:,end-1);
% duydy(2:end-1,:)=(linear_uy(3:end,:)-linear_uy(1:end-2,:))/dy/2;
% duydy(1,:)=duydy(2,:);
% duydy(end,:)=duydy(end-1,:);
% duxdy(2:end-1,:)=(linear_ux(3:end,:)-linear_ux(1:end-2,:))/dy/2;
% duxdy(1,:)=duxdy(2,:);
% duxdy(end,:)=duxdy(end-1,:);
% duydx(:,2:end-1)=(linear_uy(:,3:end)-linear_uy(:,1:end-2))/dx;
% duydx(1,:)=duydx(2,:);
% duydx(end,:)=duydx(end-1,:);
% 
% %2nd spatial derivatives
% 
% dduxdxdx(:,2:end-1)=(duxdx(:,3:end)-duxdx(:,1:end-2))/dx/2;
% dduxdxdx(:,1)=dduxdxdx(:,2);
% dduxdxdx(:,end)=dduxdxdx(:,end-1);
% 
% dduxdxdy(2:end-1,:)=(duxdx(3:end,:)-duxdx(1:end-2,:))/dy/2;
% dduxdxdy(1,:)=dduxdxdy(2,:);
% dduxdxdy(end,:)=dduxdxdy(end-1,:);
% 
% dduxdydy(2:end-1,:)=(duxdy(3:end,:)-duxdy(1:end-2,:))/dy/2;
% dduxdydy(1,:)=dduxdydy(2,:);
% dduxdydy(end,:)=dduxdydy(end-1,:);
% 
% dduydydy(2:end-1,:)=(duydy(3:end,:)-duydy(1:end-2,:))/dy/2;
% dduydydy(1,:)=dduydydy(2,:);
% dduydydy(end,:)=dduydydy(end-1,:);
% 
% dduydxdy(2:end-1,:)=(duydx(3:end,:)-duydx(1:end-2,:))/dy/2;
% dduydxdy(1,:)=dduydxdy(2,:);
% dduydxdy(end,:)=dduydxdy(end-1,:);
% 
% dduydxdx(:,2:end-1)=(duydx(:,3:end)-duydx(:,1:end-2))/dx/2;
% dduydxdx(:,1)=dduydxdx(:,2);
% dduydxdx(:,end)=dduydxdx(:,end-1);
% 


% %%
% duxdx(2:end-1,:)=(ux(3:end,:)-ux(1:end-2,:))/dx/2;
% duydy(:,2:end-1)=(uy(:,3:end)-uy(:,1:end-2))/dy/2;
% duxdy(:,2:end-1)=(ux(:,3:end)-ux(:,1:end-2))/dy/2;
% duydx(2:end-1,:)=(uy(3:end,:)-uy(1:end-2,:))/dx/2;
% 
% %2nd spatial derivatives
% dduxdxdx(2:end-1,:)=(duxdx(3:end,:)-duxdx(1:end-2,:))/dx/2;
% dduxdxdy(:,2:end-1)=(duxdx(:,3:end)-duxdx(:,1:end-2))/dy/2;
% dduxdydy(:,2:end-1)=(duxdy(:,3:end)-duxdy(:,1:end-2))/dy/2;
% dduydydy(:,2:end-1)=(duydy(:,3:end)-duydy(:,1:end-2))/dy/2;
% dduydxdy(:,2:end-1)=(duydx(:,3:end)-duydx(:,1:end-2))/dy/2;
% dduydxdx(2:end-1,:)=(duydx(3:end,:)-duydx(1:end-2,:))/dx/2;



Fx=(mu+A/4)*(dduxdxdx.*duxdx +dduxdydy.*duxdx +dduydxdx.*duydx +dduydydy.*duydx + dduxdxdx.*duxdx +dduxdydy.*duxdx +dduydxdx.*duxdy +dduydydy.*duxdy +2*(dduxdxdx.*duxdx +dduxdxdy.*duxdy +dduxdxdy.*duydx +dduxdydy.*duydy))...
    +(K+mu/3 +A/4 +B)*(dduxdxdx.*duxdx +dduxdxdy.*duxdy +dduydxdx.*duydx+ dduydxdy.*duydy +dduxdxdx.*duxdx +dduydxdy.*duxdx +dduxdxdy.*duxdy +dduydydy.*duxdy)...
    +(K-2*mu/3 +B)*(dduxdxdx.*duxdx +dduxdydy.*duxdx +dduxdxdx.*duydy +dduxdydy.*duydy)...
    +(A/4+B)*(dduxdxdx.*duxdx +dduydxdy.*duxdx +dduxdxdy.*duydx +dduydydy.*duydx + dduxdxdx.*duxdx +dduxdxdy.*duydx +dduydxdx.*duxdy +dduydxdy.*duydy)...
    +(B+2*C)*(dduxdxdx.*duxdx +dduydxdy.*duxdx +dduxdxdx.*duydy +dduydxdy.*duydy);

Fy=(mu+A/4)*(dduxdxdx.*duxdy +dduxdydy.*duxdy +dduydxdx.*duydy +dduydydy.*duydy +dduxdxdx.*duydx +dduxdydy.*duydx +dduydxdx.*duydy +dduydydy.*duydy +2*(dduydxdx.*duxdx +dduydxdy.*duxdy +dduydxdy.*duydx +dduydydy.*duydy))...
    +(K+mu/3 +A/4 +B)*(dduxdxdy.*duxdx +dduxdydy.*duxdy +dduydxdy.*duydx +dduydydy.*duydy +dduxdxdx.*duydx +dduydxdy.*duydx +dduxdxdy.*duydy +dduydydy.*duydy)...
    +(K-2*mu/3 +B)*(dduydxdx.*duxdx +dduydydy.*duxdx +dduydxdx.*duydy +dduydydy.*duydy)...
    +(A/4+B)*(dduxdxdx.*duxdy +dduydxdy.*duxdy +dduxdxdy.*duydy +dduydydy.*duydy +dduxdxdy.*duxdx +dduxdydy.*duydx +dduydxdy.*duxdy +dduydydy.*duydy)...
    +(B+2*C)*(dduxdxdy.*duxdx +dduydydy.*duxdx +dduxdxdy.*duydy +dduydydy.*duydy);

end