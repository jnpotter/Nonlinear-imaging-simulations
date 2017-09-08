function [ux_dotdot uy_dotdot] = nonlin_wave_eqs_high_order(ux, uy, ux_dot, uy_dot)

global mu rho K dx dy n m gpu_on A B C 

global duxdxh duxdyh duydxh duydyh dduxdxdxh dduxdxdyh dduxdydyh dduydydyh dduydxdyh dduydxdxh

%1st order derivatives (4th order accuracy)
duxdxh(:,3:end-2)=(1/12*ux(:,(3-2):(end-2-2)) - 2/3*ux(:,(3-1):(end-2-1)) +  2/3*ux(:,(3+1):(end-2+1)) - 1/12*ux(:,(3+2):(end-2+2)))/dx; 
duydyh(3:end-2,:)=(1/12*uy((3-2):(end-2-2),:) - 2/3*uy((3-1):(end-2-1),:) +  2/3*uy((3+1):(end-2+1),:) - 1/12*uy((3+2):(end-2+2),:))/dy;
duxdyh(3:end-2,:)=(1/12*ux((3-2):(end-2-2),:) - 2/3*ux((3-1):(end-2-1),:) +  2/3*ux((3+1):(end-2+1),:) - 1/12*ux((3+2):(end-2+2),:))/dy;
duydxh(:,3:end-2)=(1/12*uy(:,(3-2):(end-2-2)) - 2/3*uy(:,(3-1):(end-2-1)) +  2/3*uy(:,(3+1):(end-2+1)) - 1/12*uy(:,(3+2):(end-2+2)))/dx;

%2nd order derivatives (4th order accuracy)

dduxdxdxh(:,3:end-2)=(-1/12*ux(:,(3-2):(end-2-2)) + 4/3*ux(:,(3-1):(end-2-1)) - 5/2*ux(:,3:end-2) + 4/3*ux(:,(3+1):(end-2+1)) - 1/12*ux(:,(3+2):(end-2+2)))/dx^2;
dduxdydyh(3:end-2,:)=(-1/12*ux((3-2):(end-2-2),:) + 4/3*ux((3-1):(end-2-1),:) - 5/2*ux(3:end-2,:) + 4/3*ux((3+1):(end-2+1),:) - 1/12*ux((3+2):(end-2+2),:))/dy^2;
dduydydyh(3:end-2,:)=(-1/12*uy((3-2):(end-2-2),:) + 4/3*uy((3-1):(end-2-1),:) - 5/2*uy(3:end-2,:) + 4/3*uy((3+1):(end-2+1),:) - 1/12*uy((3+2):(end-2+2),:))/dy^2;
dduydxdxh(:,3:end-2)=(-1/12*uy(:,(3-2):(end-2-2)) + 4/3*uy(:,(3-1):(end-2-1)) - 5/2*uy(:,3:end-2) + 4/3*uy(:,(3+1):(end-2+1)) - 1/12*uy(:,(3+2):(end-2+2)))/dx^2;

dduxdxdyh(3:end-2,3:end-2)=...
+(1/12)*(1/12*ux((3-2):(end-2-2),(3-2):(end-2-2)) - 2/3*ux((3-2):(end-2-2),(3-1):(end-2-1)) + 2/3*ux((3-2):(end-2-2),(3+1):(end-2+1)) - 1/12*ux((3-2):(end-2-2),(3+2):(end-2+2)))/dx/dy...
+(-2/3)*(1/12*ux((3-1):(end-2-1),(3-2):(end-2-2)) - 2/3*ux((3-1):(end-2-1),(3-1):(end-2-1)) + 2/3*ux((3-1):(end-2-1),(3+1):(end-2+1)) - 1/12*ux((3-1):(end-2-1),(3+2):(end-2+2)))/dx/dy...
+(2/3)*(1/12*ux((3+1):(end-2+1),(3-2):(end-2-2)) - 2/3*ux((3+1):(end-2+1),(3-1):(end-2-1)) + 2/3*ux((3+1):(end-2+1),(3+1):(end-2+1)) - 1/12*ux((3+1):(end-2+1),(3+2):(end-2+2)))/dx/dy...
+(-1/12)*(1/12*ux((3+2):(end-2+2),(3-2):(end-2-2)) - 2/3*ux((3+2):(end-2+2),(3-1):(end-2-1)) + 2/3*ux((3+2):(end-2+2),(3+1):(end-2+1)) - 1/12*ux((3+2):(end-2+2),(3+2):(end-2+2)))/dx/dy;

dduydxdyh(3:end-2,3:end-2)=...
+(1/12)*(1/12*uy((3-2):(end-2-2),(3-2):(end-2-2)) - 2/3*uy((3-2):(end-2-2),(3-1):(end-2-1)) + 2/3*uy((3-2):(end-2-2),(3+1):(end-2+1)) - 1/12*uy((3-2):(end-2-2),(3+2):(end-2+2)))/dx/dy...
+(-2/3)*(1/12*uy((3-1):(end-2-1),(3-2):(end-2-2)) - 2/3*uy((3-1):(end-2-1),(3-1):(end-2-1)) + 2/3*uy((3-1):(end-2-1),(3+1):(end-2+1)) - 1/12*uy((3-1):(end-2-1),(3+2):(end-2+2)))/dx/dy...
+(2/3)*(1/12*uy((3+1):(end-2+1),(3-2):(end-2-2)) - 2/3*uy((3+1):(end-2+1),(3-1):(end-2-1)) + 2/3*uy((3+1):(end-2+1),(3+1):(end-2+1)) - 1/12*uy((3+1):(end-2+1),(3+2):(end-2+2)))/dx/dy...
+(-1/12)*(1/12*uy((3+2):(end-2+2),(3-2):(end-2-2)) - 2/3*uy((3+2):(end-2+2),(3-1):(end-2-1)) + 2/3*uy((3+2):(end-2+2),(3+1):(end-2+1)) - 1/12*uy((3+2):(end-2+2),(3+2):(end-2+2)))/dx/dy;

%really crude finite difference gradient estimate (improve later)
%1st spatial derivatives

% duxdxh(:,2:end-1)=(ux(:,3:end)-ux(:,1:end-2))/dx/2; %central difference 1st order
% duxdxh(:,1)=(ux(:,2)-ux(:,1))/dx;
% duxdxh(:,end)=(ux(:,end)-ux(:,end-1))/dx;
% 
% duydyh(2:end-1,:)=(uy(3:end,:)-uy(1:end-2,:))/dy/2;
% duydyh(1,:)=(uy(2,:)-uy(1,:))/dy;
% duydyh(end,:)=(uy(end,:)-uy(end-1,:))/dy;
% 
% duxdyh(2:end-1,:)=(ux(3:end,:)-ux(1:end-2,:))/dy/2;
% duxdyh(1,:)=(ux(2,:)-ux(1,:))/dy;
% duxdyh(end,:)=(ux(end,:)-ux(end-1,:))/dy;
% 
% duydxh(:,2:end-1)=(uy(:,3:end)-uy(:,1:end-2))/dx/2;
% duydxh(:,1)=(ux(:,2)-ux(:,1))/dx;
% duydxh(:,end)=(ux(:,end)-ux(:,end-1))/dx;

% %2nd spatial derivatives
% dduxdxdxh(:,2:end-1)=(duxdxh(:,3:end)-duxdxh(:,1:end-2))/dx/2;
% dduxdxdxh(:,1)=(duxdxh(:,2)-duxdxh(:,1))/dx;
% dduxdxdxh(:,end)=(duxdxh(:,end)-duxdxh(:,end-1))/dx;
% 
% dduxdxdyh(2:end-1,:)=(duxdxh(3:end,:)-duxdxh(1:end-2,:))/dy/2;
% dduxdxdyh(1,:)=(duxdxh(2,:)-duxdxh(1,:))/dy;
% dduxdxdyh(end,:)=(duxdxh(end,:)-duxdxh(end-1,:))/dy;
% 
% dduxdydyh(2:end-1,:)=(duxdyh(3:end,:)-duxdyh(1:end-2,:))/dy/2;
% dduxdydyh(1,:)=(duxdyh(2,:)-duxdyh(1,:))/dy;
% dduxdydyh(end,:)=(duxdyh(end,:)-duxdyh(end-1,:))/dy;
% 
% dduydydyh(2:end-1,:)=(duydyh(3:end,:)-duydyh(1:end-2,:))/dy/2;
% dduydydyh(1,:)=(duydyh(2,:)-duydyh(1,:))/dy;
% dduydydyh(end,:)=(duydyh(end,:)-duydyh(end-1,:))/dy;
% 
% dduydxdyh(2:end-1,:)=(duydxh(3:end,:)-duydxh(1:end-2,:))/dy/2;
% dduydxdyh(1,:)=(duydxh(2,:)-duydxh(1,:))/dy;
% dduydxdyh(end,:)=(duydxh(end,:)-duydxh(end-1,:))/dy;
% 
% dduydxdxh(:,2:end-1)=(duydxh(:,3:end)-duydxh(:,1:end-2))/dx/2;
% dduydxdxh(:,1)=(duydxh(:,2)-duydxh(:,1))/dx;
% dduydxdxh(:,end)=(duydxh(:,end)-duydxh(:,end-1))/dx;


%Nonlinear terms
Fx=(mu+A/4)*(dduxdxdxh.*duxdxh +dduxdydyh.*duxdxh +dduydxdxh.*duydxh +dduydydyh.*duydxh + dduxdxdxh.*duxdxh +dduxdydyh.*duxdxh +dduydxdxh.*duxdyh +dduydydyh.*duxdyh +2*(dduxdxdxh.*duxdxh +dduxdxdyh.*duxdyh +dduxdxdyh.*duydxh +dduxdydyh.*duydyh))...
    +(K+mu/3 +A/4 +B)*(dduxdxdxh.*duxdxh +dduxdxdyh.*duxdyh +dduydxdxh.*duydxh+ dduydxdyh.*duydyh +dduxdxdxh.*duxdxh +dduydxdyh.*duxdxh +dduxdxdyh.*duxdyh +dduydydyh.*duxdyh)...
    +(K-2*mu/3 +B)*(dduxdxdxh.*duxdxh +dduxdydyh.*duxdxh +dduxdxdxh.*duydyh +dduxdydyh.*duydyh)...
    +(A/4+B)*(dduxdxdxh.*duxdxh +dduydxdyh.*duxdxh +dduxdxdyh.*duydxh +dduydydyh.*duydxh + dduxdxdxh.*duxdxh +dduxdxdyh.*duydxh +dduydxdxh.*duxdyh +dduydxdyh.*duydyh)...
    +(B+2*C)*(dduxdxdxh.*duxdxh +dduydxdyh.*duxdxh +dduxdxdxh.*duydyh +dduydxdyh.*duydyh);

Fy=(mu+A/4)*(dduxdxdxh.*duxdyh +dduxdydyh.*duxdyh +dduydxdxh.*duydyh +dduydydyh.*duydyh +dduxdxdxh.*duydxh +dduxdydyh.*duydxh +dduydxdxh.*duydyh +dduydydyh.*duydyh +2*(dduydxdxh.*duxdxh +dduydxdyh.*duxdyh +dduydxdyh.*duydxh +dduydydyh.*duydyh))...
    +(K+mu/3 +A/4 +B)*(dduxdxdyh.*duxdxh +dduxdydyh.*duxdyh +dduydxdyh.*duydxh +dduydydyh.*duydyh +dduxdxdxh.*duydxh +dduydxdyh.*duydxh +dduxdxdyh.*duydyh +dduydydyh.*duydyh)...
    +(K-2*mu/3 +B)*(dduydxdxh.*duxdxh +dduydydyh.*duxdxh +dduydxdxh.*duydyh +dduydydyh.*duydyh)...
    +(A/4+B)*(dduxdxdxh.*duxdyh +dduydxdyh.*duxdyh +dduxdxdyh.*duydyh +dduydydyh.*duydyh +dduxdxdyh.*duxdxh +dduxdydyh.*duydxh +dduydxdyh.*duxdyh +dduydydyh.*duydyh)...
    +(B+2*C)*(dduxdxdyh.*duxdxh +dduydydyh.*duxdxh +dduxdxdyh.*duydyh +dduydydyh.*duydyh);


%PDES
ux_dotdot=(mu/rho)*(dduxdxdxh +dduxdydyh)+((K+mu/3)/rho)*(dduxdxdxh +dduydxdyh) +Fx/rho;
uy_dotdot=(mu/rho)*(dduydxdxh +dduydydyh)+((K+mu/3)/rho)*(dduxdxdyh +dduydydyh) +Fy/rho;

% %edge boundary correction
% ux_dotdot(:,end)=ux_dotdot(:,end-1); ux_dotdot(1,:)=ux_dotdot(2,:); ux_dotdot(end,:)=ux_dotdot(end-1,:); ux_dotdot(:,1)=ux_dotdot(:,2);
% uy_dotdot(:,end)=uy_dotdot(:,end-1); uy_dotdot(1,:)=uy_dotdot(2,:); uy_dotdot(end,:)=uy_dotdot(end-1,:); uy_dotdot(:,1)=uy_dotdot(:,2);

%ux_dotdot(1,1)=(ux_dotdot(2,1)+ux_dotdot(1,2)+ux_dotdot(2,2))/3; ux_dotdot(end,1)=(ux_dotdot(end,2)+ux_dotdot(end-1,1)+ux_dotdot(end-1,2))/3;
%ux_dotdot(1,end)=(ux_dotdot(2,end)+ux_dotdot(1,end-1)+ux_dotdot(2,end-1))/3; ux_dotdot(end,end)=(ux_dotdot(end-1,end)+ux_dotdot(end,end-1)+ux_dotdot(end-1,end-1))/3;

end

