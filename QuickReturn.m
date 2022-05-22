% This program solves kinematics of a two link robot
clear all
close all
clc
%remember to calculate the derivatives outside the for loop
nb=5; %number of bodies
% q=sym('q',[1 3*nb]);
Rx=sym('Rx',[1 nb]);
Ry=sym('Ry',[1 nb]);
Th=sym('Th',[1 nb]);
Rxd=sym('Rxd',[1 nb]);
Ryd=sym('Ryd',[1 nb]);
Thd=sym('Thd',[1 nb]);

qNum=[0 0 0 0 0.1 pi/2 0 0-0.05 pi/2 0.2 0.15 0 0.3 0.2 0];
t0=0;dt=0.05;tf=4;

for i=1:nb
    q(3*(i-1)+1)=Rx(i);
    q(3*(i-1)+2)=Ry(i);
    q(3*(i-1)+3)=Th(i);
    qd(3*(i-1)+1)=Rxd(i);
    qd(3*(i-1)+2)=Ryd(i);
    qd(3*(i-1)+3)=Thd(i);
end

syms L1 L2 L3 L4 L5 omg2 Th20 t

CSym1=[Rx(1);
   Ry(1);      
   Th(1);
   Rx(2)-0.5*L2*cos(Th(2));
   Ry(2)-0.5*L2*sin(Th(2));
   Rx(3)-0.5*L3*cos(Th(3));
   Ry(3)-0.5*L3*sin(Th(3))+0.3;
   Rx(3)+0.5*L3*cos(Th(3))-Rx(4)+0.5*L4*cos(Th(4));
   Ry(3)+0.5*L3*sin(Th(3))-Ry(4)+0.5*L4*sin(Th(4));
   Rx(4)+0.5*L4*cos(Th(4))-Rx(5);
   Ry(4)+0.5*L4*sin(Th(4))-Ry(5);
   Th(5);
   Ry(5)-0.2;
   (Rx(3)-Rx(2)-0.5*L2*cos(Th(2)))*(-sin(Th(3)))+(Ry(3)-Ry(2)-0.5*L2*sin(Th(2)))*(cos(Th(3)));
   Th(2)-omg2*t-Th20;];

% C=subs(C1,qb,q)
% Cq=jacobian(C,q)

CqSym=jacobian(CSym1,q);

C1=subs(CSym1,{L1 L2 L3 L4 omg2  Th20},{0 0.2 0.5 0.3 2 pi/2});
Cq1=subs(CqSym,{L1 L2 L3 L4 omg2  Th20},{0 0.2 0.5 0.3 2 pi/2});
ErrorTol=1e-5;ErrorTol2=1e-3;
Idx=1;
Ct1=diff(C1,t);
Cqt1=diff(Cq1,t);

Cq_qd=(Cq1*qd.');
term1_1=jacobian(Cq_qd,q);

%Solutin for C, Cq, Qd
CRes=subs(CSym1,{L1 L2 L3 L4 L5},{0 0.2 0.5 0.3 0});
CqRes=subs(CqSym, {L1 L2 L3 L4 L5},{0 0.2 0.5 0.3 0});
Ctt1=diff(Ct1,t);

QdResSem=-(term1_1*qd.'+2*Cqt1*qd.'+Ctt1);
QdRes=subs(QdResSem,{L1 L2 L3 L4 L5},{0 0.2 0.5 0.3 0});
%QdRes=simple(QdRes);

for tn=t0:dt:tf
% Position analysis
if tn>t0
    qNum=qNum+dt*qdNum';
end
C=subs(C1,t,tn);
Cq=subs(Cq1,t,tn);
qNum=eval(NewtonRaphsonSol(qNum,q,C,Cq,ErrorTol,ErrorTol2,1));

% Velocity analysis
Ct=subs(Ct1,t,tn);
Cq=subs(Cq1,t,tn);
Cqnum=subs(Cq,q,qNum);
qdNum=eval(-Cqnum\Ct);

%Acceleration analysis
Ctt1=diff(Ct1,t);
Ctt=subs(Ctt1,t,tn);

Cqt=subs(Cqt1,t,tn);
Cqtnum=subs(Cqt,q,qNum);
term2=2*Cqtnum*qdNum;

term1_2=subs(term1_1,t,tn);
term1_3=subs(term1_2,q,qNum);
term1_4=subs(term1_3,qd,qdNum');

term1=term1_4*qdNum; 
Qd=-(term1+term2+Ctt);
qddNum=eval(Cqnum\Qd);

Res(Idx,:)=[tn,qNum,qdNum',qddNum'];
Idx=Idx+1;
end
