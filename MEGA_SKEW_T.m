%Ashton Pohlmann lab10

clear variables;
close all; 
%[P,H,T,Td,u,v]=textread('OAX_12_May_2014_00z.dat','%f %f %f %f %f %f');
[P,H,T,Td,u,v]=textread('TBW_sounding.dat','%f %f %f %f %f %f');
%Used to find the values for the table
%T(1,1)=T(1,1)+5
%T(1,1)=T(1,1)-5
%Td(1,1)=Td(1,1)+5
%Td(1,1)=Td(1,1)-5
%Get the Cape/Cin Values that we need to calculate
[cape,cin] = getcape(83,P,T,Td);

%Totals Totals Function Below
TotalsTotals(T,Td,P)


es2=satvape(Td(1,1));

%SKEW T 
T = c2k(T);
Td= c2k(Td);

TLCL=(2840/(3.5*log(T(1,1))-log(es2)-4.805))+55;
PLCL = P(1,1)*(TLCL/T(1,1))^(0.286);




%pressure levels
plevel50 = [1050:-50:100]';
%temp levels
tlevel10 = -80:10:150;
tlevel03 = -40:10:110;

%temp/temp dew to x(skew'd)
At=tp2x(T,P);
Atd=tp2x(Td,P);

%this stuff creates the chart
skewt = figure(3);
hold on

%pressure lines
yline(plevel50)
set(gca,'yscale','log','ydir','reverse')
set(gca,'ytick',100:50:1050)

%skewed temp lines
T_skew=(tlevel10-(40.*log(0.001.*plevel50))+1.9516);
plot(T_skew,plevel50,'black');

%dry adiabats and mixing ratio lines

%constants things
z=0:.5:25;
b=1050*exp(-z/10);
mxcon=[0.4, 1, 2, 3, 5, 8, 12, 16, 20];
thecontours=-40:10:110;

%calculate dry adiabats and sat stuff for the second part of the lab
nj=95;
nz=length(z);
for i=1:nj
    for j=1:nz
        %another new temperature matrix
        NTM(i,j)=i-2*j-40;
        %dry adiabats/constant potential temp
        pt(i,j)=PotTemp(NTM(i,j),b(j));
        %sat vapor pressure
        es(i,j) = satvape(NTM(i,j));
        %sat mixing ratio
        ws(i,j)=mixingsat(es(i,j),b(j));
        %moist adiabat
        th(i,j)=equ_pot_temp_sat(NTM(i,j),b(j),ws(i,j));
    end 
end

for i=1:nj
    for j=1:nz
        ntm2(i,j)=TLCL;
the1=equ_pot_temp_sat(ntm2(i,j),b(j),ws(i,j));
    end
end






%plotting the lines for lab 5
ptline=contour(NTM(:,1),b,pt'-273.15,tlevel03,'r-');
mline=contour(NTM(:,1),b,ws',mxcon,'color',[0 0.5 1],'linestyle','--');
contour(NTM(:,1),b,th'-273.15,thecontours,'g-')
clabel(mline, 'color',[0 0.5 1])

%plot data
plot(At,P,'r','linewidth',1.5);
plot(Atd,P,'b','linewidth',1.5)
xlim([-40 40])
ylim([100 1050])
title("SkewT LogP")
xlabel("Temperature")
ylabel("Pressure")
hold off

%funtions used above
function x=tp2x(t,p)
slope=40;
x=t-slope.*log(p);
end
function CK = c2k(x)
 CK = (x+273.15);
end
function pt = PotTemp(x,y)
    pt = (x+273.15).*(1000./y).^0.286;
end
%method one equations to make mixing ratio lines.
function x=satvape(y)
    x = 6.112*exp((17.67*y)/(y+243.5));
end
function x = mixingsat(a,b)
    x = ((0.622.*a)./(b-a))*1000;
    %a is es and b is pressure
end
function eqts = equ_pot_temp_sat(t,p,ws)
    t = t+273.15;
    eqts = (t.*(1000./p)^(0.2854.*(1-(0.28.*10^-3).*ws))).*exp(ws.*((3.376./t)-0.00254).*(1+ws.*(0.81.*10^-3)));
end

function TT = TotalsTotals(T,Td,P);
ind = find(P==850);
ind2 = find(P==500);
VT = T(ind)-T(ind2);
CT = Td(ind)-T(ind2);
TT = VT+CT;
end

function e = w2e(p,w)
    e = (p.*w)/(0.622+w);
end




function [cape, cin]=getcape(nk, p_in, t_in, td_in);
%-----------------------------------------------------------------------

% Original code in fortran by G. Bryan translated into MATLAB
% Input variables: nk - number of levels in the sounding (integer)
%        p_in - one-dimensional array of pressure (mb) (real)
%        t_in - one-dimensional array of temperature (C) (real)
%        td_in - one-dimensional array of dewpoint temperature (C) (real)
%
%  Output variables:  cape - Convective Available Potential Energy (J/kg) (real)
%                     cin - Convective Inhibition (J/kg) (real)

%  getcape - a fortran90 subroutine to calculate Convective Available
%            Potential Energy (CAPE) from a sounding.
%
%  Version 1.02                           Last modified:  10 October 2008
%
%  Author:  George H. Bryan
%           Mesoscale and Microscale Meteorology Division
%           National Center for Atmospheric Research
%           Boulder, Colorado, USA
%           gbryan@ucar.edu
%
%  Disclaimer:  This code is made available WITHOUT WARRANTY.
%
%  References:  Bolton (1980, MWR, p. 1046) (constants and definitions)
%               Bryan and Fritsch (2004, MWR, p. 2421) (ice processes)
%
%-----------------------------------------------------------------------
%
%  Input:     nk - number of levels in the sounding (integer)
%
%           p_in - one-dimensional array of pressure (mb) (real)
%
%           t_in - one-dimensional array of temperature (C) (real)
%
%          td_in - one-dimensional array of dewpoint temperature (C) (real)
%
%  Output:  cape - Convective Available Potential Energy (J/kg) (real)
%
%            cin - Convective Inhibition (J/kg) (real)
%
%-----------------------------------------------------------------------

clear global; 
clear functions;

persistent adiabat avgqv avgth b1 b2 cloud converge cp cpdg cpdrd cpi cpl cpm cpv debug_level doit dp dz eps fice fliq frac g i ice k kmax lfc lhf lhs lhv ls1 ls2 lv1 lv2 maxthe ml_depth n narea nloop not_converged orec p p00 p1 p2 parea pb pc pi pi1 pi2 pinc pn ps pt ptv q qi1 qi2 qibar ql1 ql2 qlbar qt qv1 qv2 qvbar rd rddcp reps rm rp00 rv source t t0 t1 t2 tbar td th th1 th2 the thlast thv thv1 thv2 xls xlv z ; 

%% Choose Accuracy

pinc = 10.0;
% (smaller number yields more accurate
%  results,larger number makes code
%  go faster)

%% Constants

g     = 9.81;
p00   = 100000.0;
cp    = 1005.7;
rd    = 287.04;
rv    = 461.5;
xlv   = 2501000.0;
xls   = 2836017.0;
t0    = 273.15;
cpv   = 1875.0;
cpl   = 4190.0;
cpi   = 2118.636;
lv1   = xlv+(cpl-cpv).*t0;
lv2   = cpl-cpv;
ls1   = xls+(cpi-cpv).*t0;
ls2   = cpi-cpv;
rp00  = 1.0./p00;
eps   = rd./rv;
reps  = rv./rd;
rddcp = rd./cp;
cpdrd = cp./rd;
cpdg  = cp./g;
converge = 0.0002;
debug_level =   0;

%% Convert p,t,td to mks units; get pi,q,th,thv

for k=1:nk;
p(k) = 100.0.*p_in(k);
t(k) = 273.15+t_in(k);
td(k) = 273.15+td_in(k);
pi(k) =(p(k).*rp00).^rddcp;
[q(k) ,p(k),td(k)]=getqvs(p(k),td(k));
th(k) = t(k)./pi(k);
thv(k) = th(k).*(1.0+reps.*q(k))./(1.0+q(k));
end; k=fix(nk+1);

%% Get height using the hydrostatic equation

z(1) = 0.0;
for k=2:nk;
dz = -cpdg.*0.5.*(thv(k)+thv(k-1)).*(pi(k)-pi(k-1));
z(k) = z(k-1) + dz;
end; k=fix(nk+1);

%% Define parcel properties at initial location

kmax = 1;
narea = 0.0;
k    = fix(kmax);
th2  = th(kmax);
pi2  = pi(kmax);
p2   = p(kmax);
t2   = t(kmax);
thv2 = thv(kmax);
qv2  = q(kmax);
b2   = 0.0;
ql2 = 0.0;
qi2 = 0.0;
qt  = qv2;
cape = 0.0;
cin  = 0.0;
lfc  = 0.0;
doit = true;
cloud = false;
ice = false;
t2_orig=t2;    
[the ,p2,t2,dumvar4,qv2]=getthe(p2,t2,t2,qv2);    
t2(dumvar4~=t2_orig)=dumvar4(dumvar4~=t2_orig);

%% The Parcel Ascends!

while( doit & (k<nk) );
k = fix(k+1);
b1 =  b2;
dp = p(k-1)-p(k);

if( dp<pinc )
nloop = 1;
else;
nloop = fix(1 + fix( dp./pinc ));
dp = dp./(nloop);
end;

for n=1:nloop;
p1 =  p2;
t1 =  t2;
pi1 = pi2;
th1 = th2;
qv1 = qv2;
ql1 = ql2;
qi1 = qi2;
thv1 = thv2;
p2 = p2 - dp;
pi2 =(p2.*rp00).^rddcp;
thlast = th1;
i = 0;
not_converged = true;

while( not_converged );
i = fix(i + 1);
t2 = thlast.*pi2;
fliq = 1.0;
fice = 0.0;
qv2 = min( qt , fliq.*getqvs(p2,t2) + fice.*getqvi(p2,t2) );
qi2 = max( fice.*(qt-qv2) , 0.0 );
ql2 = max( qt-qv2-qi2 , 0.0 );
tbar  = 0.5.*(t1+t2);
qvbar = 0.5.*(qv1+qv2);
qlbar = 0.5.*(ql1+ql2);
qibar = 0.5.*(qi1+qi2);
lhv = lv1-lv2.*tbar;
lhs = ls1-ls2.*tbar;
lhf = lhs-lhv;
rm=rd+rv.*qvbar;
cpm=cp+cpv.*qvbar+cpl.*qlbar+cpi.*qibar;
th2=th1.*exp(  lhv.*(ql2-ql1)./(cpm.*tbar)+lhs.*(qi2-qi1)./(cpm.*tbar)+(rm./cpm-rd./cp).*log(p2./p1) );

if(i>100)
fprintf('Error: lack of convergence')
return
end;

if( abs(th2-thlast)>converge )
thlast=thlast+0.3.*(th2-thlast);
else;
not_converged = false;
end;

end;

% Latest pressure increment is complete.  Calculate some important stuff:

if( ql2>=1.0e-10 )
cloud = true;
end;

% pseudoadiabat
qt  = qv2;
ql2 = 0.0;
qi2 = 0.0;
end; n=fix(nloop+1);

thv2 = th2.*(1.0+reps.*qv2)./(1.0+qv2+ql2+qi2);
b2 = g.*( thv2-thv(k) )./thv(k);
dz = -cpdg.*0.5.*(thv(k)+thv(k-1)).*(pi(k)-pi(k-1));
t2_orig=t2;    
[the ,p2,t2,dumvar4,qv2]=getthe(p2,t2,t2,qv2);    
t2(dumvar4~=t2_orig)=dumvar4(dumvar4~=t2_orig);

%% Get contributions to CAPE and CIN:

if((b2>=0.0) && (b1<0.0))
% first trip into positive area
ps = p(k-1)+(p(k)-p(k-1)).*(0.0-b1)./(b2-b1);
frac = b2./(b2-b1);
parea =  0.5.*b2.*dz.*frac;
narea = narea-0.5.*b1.*dz.*(1.0-frac);
cin  = cin  + narea;
narea = 0.0;

elseif((b2<0.0) && (b1>0.0));
% first trip into neg area
ps = p(k-1)+(p(k)-p(k-1)).*(0.0-b1)./(b2-b1);
frac = b1./(b1-b2);
parea =  0.5.*b1.*dz.*frac;
narea = -0.5.*b2.*dz.*(1.0-frac);

elseif(b2<0.0);
% still collecting negative buoyancy
parea =  0.0;
narea = narea-0.5.*dz.*(b1+b2);

else;
% still collecting positive buoyancy
parea =  0.5.*dz.*(b1+b2);
narea =  0.0;
end;

cape = cape + max(0.0,parea);

if((p(k)<=10000.0)&&(b2<0.0))
% stop if b < 0 and p < 100 mb
doit = false;
end;
end;


%% Done.
return;
end %subroutine getcape

%% Functions

function [getqvsresult,p,t]=getqvs(p,t);
getqvsresult=[];
persistent eps es ; 
if isempty(es), es=0; end;
if isempty(eps), eps = 287.04./461.5; end;
es = 611.2.*exp(17.67.*(t-273.15)./(t-29.65));
getqvsresult = eps.*es./(p-es);
return;
end %function getqvs

function [getqviresult,p,t]=getqvi(p,t);
getqviresult=[];
persistent eps es ; 
if isempty(es), es=0; end;
if isempty(eps), eps = 287.04./461.5; end;
es = 611.2.*exp(21.8745584.*(t-273.15)./(t-7.66));
getqviresult = eps.*es./(p-es);
return;
end %function getqvi

function [gettheresult,p,t,td,q]=getthe(p,t,td,q);
gettheresult=[];
persistent tlcl ; 
if isempty(tlcl), tlcl=0; end;
if((td-t)>=-0.1 )
tlcl = t;
else;
tlcl = 56.0 +((td-56.0).^(-1) + 0.00125.*log(t./td) ).^(-1);
end;
gettheresult=t.*((100000.0./p).^(0.2854.*(1.0-0.28.*q)) ).*exp(((3376.0./tlcl)-2.54).*q.*(1.0+0.81.*q) );
return;
end %function getthe

    