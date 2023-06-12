clear all
close all

V = 5;
mu = 5;
F = 1;
qp1 = 4;
qp2 = 4;
yp1s = 0.51;
yp2s = 0.51;
sF = 5;



%dVdt = @(x,s,p1,p2,V) F;

%finding steady state values for different values of sF

for j = 1:1:90
    clear J;
    n=100;
    dxdt = @(x,s,p1,p2) (mu - F/V)*x;
    dsdt = @(x,s,p1,p2) -(qp1/yp1s + qp2/yp2s)*x + (F/V)*(sF-s);
    dp1dt = @(x,s,p1,p2) qp1*x - (F/V)*p1;
    dp2dt = @(x,s,p1,p2) qp2*x - (F/V)*p2;
    x_in = [10;10;10;10];
    for i=1:n
       
        dxdi = Jacobian_m(dxdt,x_in);
        dsdi = Jacobian_m(dsdt,x_in);
        dp1di = Jacobian_m(dp1dt,x_in);
        dp2di = Jacobian_m(dp2dt,x_in);
        %dVdi = Jacobian_m(dVdt,x_in);
    
        J = [dxdi;dsdi;dp1di;dp2di];
    
        dxdt_val = dxdt(x_in(1),x_in(2),x_in(3),x_in(4));
        dsdt_val = dsdt(x_in(1),x_in(2),x_in(3),x_in(4));
        dp1dt_val = dp1dt(x_in(1),x_in(2),x_in(3),x_in(4));
        dp2dt_val = dp2dt(x_in(1),x_in(2),x_in(3),x_in(4));
       %dVdt_val = dVdt(x_in(1),x_in(2),x_in(3),x_in(4));
      
    
        Fvec = [dxdt_val;dsdt_val;dp1dt_val;dp2dt_val];
        xnew = x_in - inv(J)*Fvec;
        

        if (abs((xnew(1) - x_in(1))) < 10^-6) && (abs((xnew(2) - x_in(2))) < 10^-6) && (abs((xnew(3) - x_in(3))) < 10^-6) && (abs((xnew(4) - x_in(4))) < 10^-6) 
            xroot = xnew;
            xroot_vec(j,1) = xroot(1);
            xroot_vec(j,2) = xroot(2);
            xroot_vec(j,3) = xroot(3);
            xroot_vec(j,4) = xroot(4);
            break
        end
    
        x_in = xnew;
    end
    sF = 5;
    sF = sF+j;
end

figure(1) %Steady state values for different values of sF
sF_vec = linspace(5,95,90);
subplot(221),plot(sF_vec,xroot_vec(:,1))
ylabel('xss')
xlabel('sF')

subplot(222),plot(sF_vec,xroot_vec(:,2))
ylabel('Sss')
xlabel('sF')

subplot(223),plot(sF_vec,xroot_vec(:,3))
ylabel('p1ss')
xlabel('sF')

subplot(224),plot(sF_vec,xroot_vec(:,4))
ylabel('p2ss')
xlabel('sF')

%finding stability of roots

sF = 5;
for k=1:length(xroot_vec)
    ssr1 = [xroot_vec(k,1);xroot_vec(k,2);xroot_vec(k,3);xroot_vec(k,4)];
    xt = @(x,s,p1,p2) (mu - F/V)*x;
    st = @(x,s,p1,p2) -(qp1/yp1s + qp2/yp2s)*x + (F/V)*(sF - s);
    p1t = @(x,s,p1,p2) qp1*x - (F/V)*p1;
    p2t = @(x,s,p1,p2) qp2*x - (F/V)*p2;
    
    jss1 = Jacobian_m(xt,ssr1);
    jss2 = Jacobian_m(st,ssr1);
    jss3 = Jacobian_m(p1t,ssr1);
    jss4 = Jacobian_m(p2t,ssr1);
    
    Jss = [jss1;jss2;jss3;jss4];
    %Jss
    [Vec,lambda] = eig(Jss);
    %lambda
    poseig = 0;
    negeig = 0;
    zeroeig = 0;
    for i=1:4
        if lambda(i,i) > 0
            poseig = poseig + 1 ;
        end
        if  lambda(i,i) < 0
            negeig = negeig + 1 ;
        end
        if lambda(i,i) == 0 
            zeroeig = zeroeig + 1 ;
        end
    end
    if poseig > 0 && negeig == 0
        natureofss(k) = "Unstable node";
    end
    if poseig == 0 && negeig >0
        natureofss(k) = "Stable node";
    end
    if poseig > 0 && negeig > 0
        natureofss(k) = "Saddle point";
    end
    if poseig > 0 && zeroeig > 0
        natureofss(k) = "Stable node";
    end
    if negeig > 0 && zeroeig > 0
        natureofss(k) = "Unstable node";
    end
    if negeig == 0 && zeroeig > 0 && poseig == 0
        natureofss(k) = "Unstable node";
    end
    sF = 5;
    sF = sF + k;
end

%phase plane plot

for m=1:length(xroot_vec)
    if abs(xroot_vec(m,1)) < 1e-10
        xroot_vec(m,1) = 0;
    end
    if abs(xroot_vec(m,2)) < 1e-10
        xroot_vec(m,2) = 0;
    end
    if abs(xroot_vec(m,3)) < 1e-10
        xroot_vec(m,3) = 0;
    end
    if abs(xroot_vec(m,4)) < 1e-10
        xroot_vec(m,4) = 0;
    end
end

figure(2) %saddle point

plot(0,0,'*')
ylabel('deviation variable 1')
xlabel('deviation variable 2')
hold on

xin = [-1000;-1000;-1000;-1000];
for i = 1:20
    c = inv(Vec)*xin;

    xhat = @(t) c(1)*exp(lambda(1,1)*t)*Vec(1,1) + c(2)*exp(lambda(2,2)*t)*Vec(1,2) + c(3)*exp(lambda(3,3)*t)*Vec(1,3) + c(4)*exp(lambda(4,4)*t)*Vec(1,4) ;
    shat = @(t) c(1)*exp(lambda(1,1)*t)*Vec(2,1) + c(2)*exp(lambda(2,2)*t)*Vec(2,2) + c(3)*exp(lambda(3,3)*t)*Vec(2,3) + c(4)*exp(lambda(4,4)*t)*Vec(2,4) ;
    p1hat = @(t) c(1)*exp(lambda(1,1)*t)*Vec(3,1) + c(2)*exp(lambda(2,2)*t)*Vec(3,2) + c(3)*exp(lambda(3,3)*t)*Vec(3,3) + c(4)*exp(lambda(4,4)*t)*Vec(3,4) ;
    p2hat = @(t) c(1)*exp(lambda(1,1)*t)*Vec(4,1) + c(2)*exp(lambda(2,2)*t)*Vec(4,2) + c(3)*exp(lambda(3,3)*t)*Vec(4,3) + c(4)*exp(lambda(4,4)*t)*Vec(4,4) ;
    
    for h = 1:15
        xhatvec(h) = xhat(h);
        shatvec(h) = shat(h);
        p1hatvec(h) = p1hat(h);
        p2hatvec(h) = p2hat(h);
    end
    plot(shatvec,xhatvec)
    plot(p1hatvec,xhatvec)
    plot(p2hatvec,xhatvec)
    plot(shatvec,p1hatvec)
    plot(shatvec,p2hatvec)
    plot(p1hatvec,p2hatvec)
    xin = [-1000;-1000;-1000;-1000];
    xin(1) = xin(1) + 100*i;
    xin(2) = xin(2) + 100*i;
    xin(3) = xin(3) + 100*i;
    xin(4) = xin(4) + 100*i;

end



