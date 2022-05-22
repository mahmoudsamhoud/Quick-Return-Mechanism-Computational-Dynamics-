function [qnum2]=NewtonRaphsonSol(qnum,q,C,Cq,ErrorTol1,ErrorTol2,Iflag)

error1=1; error2=1;
    Cnum=subs(C,q,qnum);
    Cqnum=subs(Cq,q,qnum);
NRcounter=0;
    
while  ((error1>ErrorTol1 ) && (error2>ErrorTol2)) && NRcounter<9
    dq=Cqnum\Cnum;
    qnum=qnum-dq';
    Cnum=subs(C,q,qnum);
    if Iflag==1
        Cqnum=subs(Cq,q,qnum);
    end
    error1=norm(dq);
    error2=norm(Cnum);
    
    NRcounter=NRcounter+1;
    
end
qnum2=qnum;