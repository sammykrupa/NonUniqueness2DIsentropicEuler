% This code searches for solutions to the algebraic equations making up the
% condition of being a subsolution fan, with contact discontinuity initial data.
% The code also symbolically verifies the solution.
% This code is associated with the paper "Contact discontinuities for 2-D isentropic Euler are unique in 1-D but wildly non-unique otherwise"
% by Sam G. Krupa and László Székelyhidi Jr.
% To use the constants which are given in the appendix of the paper, uncomment
% lines 266-303 (inclusive) of this code. 


%number of subsolution fans (pizza slices)
N=3; %N>1

%the epsilon determines the strictness in strict inequalities
epsilon_fixed=1;


%Setup the numerical solver to find a collection of subsolution fans, and run the solver

exitflag=0;


%keep looping (with different random initial conditions for the solver) until the solver algorithm finds a viable solution
while exitflag~=1
    




%setup the optimization variables   
alpha=optimvar('alpha',N);
beta=optimvar('beta',N);
gamma=optimvar('gamma',N);
delta=optimvar('delta',N);

%make v_{-} and v_{+}
vminusone=optimvar('vminusone',1,1);
vplusone=optimvar('vplusone',1,1);
vminustwo=optimvar('vminustwo',1,1);
vplustwo=optimvar('vplustwo',1,1);
%a vector with (in this order) \rho_{-}, \rho_1,\rho_2,...,\rho_{N}, \rho_{+}
rho=optimvar('rho',N+2,'LowerBound',epsilon_fixed);
%corresponding energies
energy=optimvar('energy',1,N+2);
%corresponding energy DERIVATIVES
energyprime=optimvar('energyprime',1,N+2,'LowerBound',epsilon_fixed);

C=optimvar('C',N);

%make \nu_{-} and \nu_{+}
numinus=optimvar('numinus',1,1);
nuplus=optimvar('nuplus',1,1);
nu=optimvar('nu',N-1);
%define the optimization problem

prob=optimproblem;

%define the constraints 

%Rankine-Hugoniot conditions on the left interface
prob.Constraints.RH_left_1=numinus*(rho(1)-rho(2))==rho(1)*vminustwo-rho(2)*beta(1);
prob.Constraints.RH_left_2=numinus*(rho(1)*vminusone-rho(2)*alpha(1))==rho(1)*vminusone*vminustwo-rho(2)*delta(1);
prob.Constraints.RH_left_3=numinus*(rho(1)*vminustwo-rho(2)*beta(1))==rho(1)*vminustwo*vminustwo+rho(2)*gamma(1)+energyprime(1)*rho(1)*rho(1)-energyprime(2)*rho(2)*rho(2)-rho(2)*C(1)*(1/2);

%Rankine-Hugoniot conditions on the interface i

for i=1:(N-1)
constraint = strcat('RH_1','_',num2str(i));
prob.Constraints.(constraint)=nu(i)*(rho(i+1)-rho(i+2))==rho(i+1)*beta(i)-rho(i+2)*beta(i+1);
 
constraint = strcat('RH_2','_',num2str(i));
prob.Constraints.(constraint)=nu(i)*(rho(i+1)*alpha(i)-rho(i+2)*alpha(i+1))==rho(i+1)*delta(i)-rho(i+2)*delta(i+1);


constraint = strcat('RH_3','_',num2str(i));
prob.Constraints.(constraint)=nu(i)*(rho(i+1)*beta(i)-rho(i+2)*beta(i+1))==-rho(i+1)*gamma(i)+rho(i+2)*gamma(i+1)+energyprime(i+1)*rho(i+1)*rho(i+1)-energyprime(i+2)*rho(i+2)*rho(i+2)+rho(i+1)*C(i)*(1/2)-rho(i+2)*C(i+1)*(1/2);

end


%Rankine-Hugoniot conditions on the right interface
prob.Constraints.RH_right_1=nuplus*(rho(N+1)-rho(N+2))==rho(N+1)*beta(N)-rho(N+2)*vplustwo;
prob.Constraints.RH_right_2=nuplus*(rho(N+1)*alpha(N)-rho(N+2)*vplusone)==rho(N+1)*delta(N)-rho(N+2)*vplusone*vplustwo;
prob.Constraints.RH_right_3=nuplus*(rho(N+1)*beta(N)-rho(N+2)*vplustwo)==-rho(N+1)*gamma(N)-rho(N+2)*vplustwo*vplustwo+energyprime(N+1)*rho(N+1)*rho(N+1)-energyprime(N+2)*rho(N+2)*rho(N+2)+rho(N+1)*C(N)*(1/2);

%subsolution conditions
for i=1:N
constraint = strcat('sub1','_',num2str(i));
prob.Constraints.(constraint)=alpha(i)*alpha(i)+beta(i)*beta(i)<=C(i)-epsilon_fixed;

constraint = strcat('sub2','_',num2str(i));
prob.Constraints.(constraint)=((1/2)*C(i)-alpha(i)*alpha(i)+gamma(i))*((1/2)*C(i)-beta(i)*beta(i)-gamma(i))-(delta(i)-alpha(i)*beta(i))^2>=epsilon_fixed;
end

%admissibility condition on the left interface
prob.Constraints.adm_left=numinus*(rho(1)*energy(1)-rho(2)*energy(2))+numinus*(rho(1)*(vminusone*vminusone+vminustwo*vminustwo)*(1/2)-rho(2)*C(1)*(1/2)) <= (rho(1)*energy(1)+energyprime(1)*rho(1)*rho(1))*vminustwo-(rho(2)*energy(2)+energyprime(2)*rho(2)*rho(2))*beta(1)+rho(1)*vminustwo*(vminusone*vminusone+vminustwo*vminustwo)*(1/2)-rho(2)*beta(1)*C(1)*(1/2);


%admissibility condition on inteface i

for i=1:(N-1)
 constraint = strcat('adm','_',num2str(i));
 prob.Constraints.(constraint)=nu(i)*(rho(i+1)*energy(i+1)-rho(i+2)*energy(i+2))+nu(i)*(rho(i+1)*C(i)*(1/2)-rho(i+2)*(C(i+1))*(1/2))<=(rho(i+1)*energy(i+1)+energyprime(i+1)*rho(i+1)*rho(i+1))*beta(i)-(rho(i+2)*energy(i+2)+energyprime(i+2)*rho(i+2)*rho(i+2))*beta(i+1)+rho(i+1)*beta(i)*C(i)*(1/2)-rho(i+2)*beta(i+1)*C(i+1)*(1/2);
end

%admissibility condition on the right interface
prob.Constraints.adm_right=nuplus*(rho(N+1)*energy(N+1)-rho(N+2)*energy(N+2))+nuplus*(rho(N+1)*C(N)*(1/2)-rho(N+2)*(vplusone*vplusone+vplustwo*vplustwo)*(1/2))<=(rho(N+1)*energy(N+1)+energyprime(N+1)*rho(N+1)*rho(N+1))*beta(N)-(rho(N+2)*energy(N+2)+energyprime(N+2)*rho(N+2)*rho(N+2))*vplustwo+rho(N+1)*beta(N)*C(N)*(1/2)-rho(N+2)*vplustwo*(vplusone*vplusone+vplustwo*vplustwo)*(1/2);

%we want the fan parition to actually make sense
prob.Constraints.fan_condition_First=numinus<=nu(1)-epsilon_fixed;
prob.Constraints.fan_condition_Last=nu(end)<=nuplus-epsilon_fixed;

for i=1:(N-2)
constraint = strcat('fan_condition','_',num2str(i));
prob.Constraints.(constraint)=nu(i)<=nu(i+1)-epsilon_fixed;
end

%let's find a contact discontinuity
prob.Constraints.contact1=vplustwo==vminustwo;
prob.Constraints.contact2=rho(1)==rho(N+2);
prob.Constraints.contact3=energy(1)==energy(N+2);
prob.Constraints.contact4=energyprime(1)==energyprime(N+2);







%we want energy function to be strictly increasing and strictly convex
%(note the lower bound on the optimization variable energyprime which is defined
%above)
for i=1:(N+1)
for j=[1:(i-1) (i+1):(N+1)]
constraint = strcat('ConstrConvexityCondition','_',num2str(i),'_',num2str(j));
prob.Constraints.(constraint)=energy(j)-energy(i)-energyprime(i)*(rho(j)-rho(i))>=epsilon_fixed;
end
end


%initialize the solver with random initial points
data.numinus=randi([-1000,0],1,1);
data.nu=zeros(N-1,1);
for i=1:(N-1)
    if i==1
        data.nu(1,1)=data.numinus+randi([1,100],1,1);
    else
        data.nu(i,1)=data.nu(i-1,1)+randi([1,100],1,1);
    end
end
data.nuplus=data.nu(end,1)+randi([0,1000],1,1);
data.alpha = randi([-1000,1000],N,1);
data.beta = randi([-1000,1000],N,1);
data.vminusone=randi([-1000,1000],1,1);
data.vminustwo=randi([-1000,1000],1,1);
data.vplusone=randi([-1000,1000],1,1);
data.vplustwo=randi([-1000,1000],1,1);
data.rho=randi([1,1000],N+2,1);
data.delta=randi([-1000,1000],N,1);
data.gamma=randi([-1000,1000],N,1);
data.C=randi([1,1000],N,1);
data.energy=randi([-1000,1000],N+2,1);
data.energyprime=randi([1,1000],N+2,1);

%let's find a contact discontinuity
data.rho(1)=data.rho(N+2);
data.vminustwo=data.vplustwo;
data.energy(1)=data.energy(N+2);
data.energyprime(1)=data.energyprime(N+2);

%or, this initial data for the solver is already known to find a solution
data.numinus=-245;
data.nu=[-220;-175];
data.nuplus=513;
data.alpha=[-282;473;-211];
data.beta=[367;408;-115];
data.vminusone=-961;
data.vminustwo=-460;
data.vplusone=-151;
data.vplustwo=-460;
data.rho=[392;822;430;888;392];
data.delta=[538;-207;617];
data.gamma=[510;-245;-568];
data.C=[791;950;328];
data.energy=[-666;-123;667;538;-666];
data.energyprime=[589;990;515;885;589];

%run this code after finding an approximate solution to make the solution
%we find even closer to exact solution (this sets the starting point of the
%iteration to be the solution found from the previous iteration)
%data=sol;

%solve!
disp('Start the solver')
%run this code after finding an approximate solution to make the solution
%we find even closer to exact solutoin (the solver might converge to an
%infeasible point, but the solution found may still be very good!)
%[sol,fval,exitflag,output] = solve(prob,data,'Options',optimoptions(prob,'EnableFeasibilityMode',true,'SubproblemAlgorithm','cg','MaxFunctionEvaluations',20000000,'MaxIterations',10000,'ConstraintTolerance',.0000000000000000000000001,'Display','iter'));
[sol,fval,exitflag,output] = solve(prob,data,'Options',optimoptions(prob,'MaxFunctionEvaluations',20000000,'MaxIterations',30000,'ConstraintTolerance',.00000000000000001,'Display','iter'));

end %suitable solution found, end the while loop



%save the original solution (before we start modifying it). Just in case
%we want to examine it later.
original_sol=sol;


%convert numerical solution from solver to symbolic

alpha=sym(sol.alpha,'f');
beta=sym(sol.beta,'f');


vminusone=sym(sol.vminusone,'f');
vplusone=sym(sol.vplusone,'f');
vminustwo=sym(sol.vminustwo,'f');
vplustwo=sym(sol.vplustwo,'f');
rho=sym(sol.rho,'f');
%corresponding energies
energy=sym(sol.energy,'f');
%corresponding energy DERIVATIVES
energyprime=sym(sol.energyprime,'f');
gamma=sym(sol.gamma,'f');
delta=sym(sol.delta,'f');

C=sym(sol.C,'f');

numinus=sym(sol.numinus,'f');
nuplus=sym(sol.nuplus,'f');
nu=sym(sol.nu,'f');



%enforce contact discontinuity (correcting small errors from the numerical
%solver)
rho(1)=rho(N+2);
vminustwo=vplustwo;
energy(1)=energy(N+2);
energyprime(1)=energyprime(N+2);


%correct the various quantities to make equalities exact

%correct the Rankine-Hugoniot conditions on the left interface
numinus=(rho(1)*vminustwo-rho(2)*beta(1))/(rho(1)-rho(2));
delta(1)=(numinus*(rho(1)*vminusone-rho(2)*alpha(1))-rho(1)*vminusone*vminustwo)/(-rho(2));
energyprime(2)=(-(numinus*(rho(1)*vminustwo-rho(2)*beta(1)))+(rho(1)*vminustwo*vminustwo+rho(2)*gamma(1)+energyprime(1)*rho(1)*rho(1)-rho(2)*C(1)*(1/2)))/(rho(2)*rho(2));

%correct the Rankine-Hugoniot conditions on the interface i (for
%i=1,...,N-2)

for i=1:(N-2)
nu(i)=(rho(i+1)*beta(i)-rho(i+2)*beta(i+1))/(rho(i+1)-rho(i+2));
 
delta(i+1)=(nu(i)*(rho(i+1)*alpha(i)-rho(i+2)*alpha(i+1))-rho(i+1)*delta(i))/(-rho(i+2));

energyprime(i+2)=(-(nu(i)*(rho(i+1)*beta(i)-rho(i+2)*beta(i+1)))+(-rho(i+1)*gamma(i)+rho(i+2)*gamma(i+1)+energyprime(i+1)*rho(i+1)*rho(i+1)+rho(i+1)*C(i)*(1/2)-rho(i+2)*C(i+1)*(1/2)))/(rho(i+2)*rho(i+2));
end



%Here, we include the exact values of the constants given in the appendix of the paper
%Uncomment these lines to use exactly the constants as given in the
%paper.
% alpha(1)=sym('-8177336068870495/140737488355328');
% alpha(2)=sym('-4833381446756075/562949953421312');
% alpha(3)=sym('3121572020159473/562949953421312');
% beta(1)=sym('-2536561643647751/140737488355328');
% beta(2)=sym('-1114286601116939/70368744177664');
% beta(3)=sym('-6219197795695073/562949953421312');
% gamma(1)=sym('841617150350781/549755813888');
% gamma(2)=sym('-2850833975067331/17592186044416');
% gamma(3)=sym('-8954832877447991/140737488355328');
% delta(1)=sym('28872176135415855785280523654056524019908546591/27627619078169805047324756605549438692229120');
% delta(2)=sym('-867454945412067709200232997995952542374584537720982074207241594308599438377/32748846874784971211058574285222379723549486466626634273206947431338475520');
% delta(3)=sym('-2871256077954219/35184372088832');
% rho(1)=sym('2708112612978501/281474976710656');
% rho(2)=sym('6811063536043807/562949953421312');
% rho(3)=sym('2057060350258899/562949953421312');
% rho(4)=sym('3062207031116133/281474976710656');
% rho(5)=sym('2708112612978501/281474976710656');
% energy(1)=sym('-5041529442624971/2199023255552');
% energy(2)=sym('-5015532875605977/2199023255552');
% energy(3)=sym('-5073206593829053/2199023255552');
% energy(4)=sym('-2515400677054201/1099511627776');
% energy(5)=sym('-5041529442624971/2199023255552');
% energyprime(1)=sym('1676289422169645/562949953421312');
% energyprime(2)=sym('60006068216738166351756926651195471316209797920723479281182349/9106752347169134708406278810788569756749285046122185710632960');
% energyprime(3)=sym('1831278218949756891087541424381121781924211556364089293813566516592411003/1359848686641341079894728790850171965963388817622134940140753492461486080');
% energyprime(4)=sym('5400383921383283/1125899906842624');
% energyprime(5)=sym('1676289422169645/562949953421312');
% numinus=sym('-6486283176597739958874052307549/196306040423407104692364247040');
% nu(1)=sym('-1153852086001065889673487658885/60824224363690518566334889984');
% nu(2)=sym('-4856156003780791/562949953421312');
% nuplus=sym('7162856387903725/562949953421312');
% vminusone=sym('-4098844157247653/70368744177664');
% vplusone=sym('3603433899522037/562949953421312');
% vminustwo=sym('-996118042660627/70368744177664');
% vplustwo=sym('-996118042660627/70368744177664');
% C(1)=sym('510415269881361/137438953472');
% C(2)=sym('1515879700153707/2199023255552');
% C(3)=sym('1855257252703141/8796093022208');


%calculate D\Gamma at the point (\hat\alpha_N,\hat\beta_N,\hat\delta_N,\hat\rho_N,\hat\nu_+,\hat\nu_{N-1})

pprime = 2*rho(N+1)*energyprime(N+1); %p'(\rho_{N})
DGamma=sym(zeros(6,6));
DGamma=[0 -rho(N+1) 0 -beta(N)+nu(N-1) 0 -(rho(N)-rho(N+1));rho(N+1)*nu(N-1) 0 -rho(N+1) -delta(N)+alpha(N)*nu(N-1) 0 -(rho(N)*alpha(N-1)-rho(N+1)*alpha(N));0 nu(N-1)*rho(N+1) 0 gamma(N)-pprime-.5*C(N)+beta(N)*nu(N-1) 0 -(rho(N)*beta(N-1)-rho(N+1)*beta(N));0 rho(N+1) 0 beta(N)-nuplus -(rho(N+1)-rho(N+2)) 0;-nuplus*rho(N+1) 0 rho(N+1) delta(N)-nuplus*alpha(N) -(rho(N+1)*alpha(N)-rho(N+2)*vplusone) 0;0 -nuplus*rho(N+1) 0 -gamma(N)+pprime+.5*C(N)-nuplus*beta(N) -(rho(N+1)*beta(N)-rho(N+2)*vplustwo) 0];

disp('Singular values of DGamma')
vpa(svd(DGamma))

%calculate some maximum of various quantities  to figure out how
%svd(DGamma) will vary locally (using 1-Lipschitz-ness of singular values)
%we use variable-precision arithmetic:
%see https://www.mathworks.com/help/symbolic/vpa.html
disp('maximum of alpha values')
vpa(max(abs(alpha)))

disp('maximum of beta values')
vpa(max(abs(beta)))

disp('maximum of gamma values')
vpa(max(abs(gamma)))

disp('maximum of delta values')
vpa(max(abs(delta)))

disp('maximum of nu values')
vpa(max(abs(nu)))

disp('maximum of C values')
vpa(max(abs(C)))

disp('maximum of rho values')
vpa(max(abs(rho)))

disp('maximum of energy values')
vpa(max(abs(energy)))

disp('maximum of energyprime values')
vpa(max(abs(energyprime)))

disp('numinus')
vpa(abs(numinus))

disp('nuplus')
vpa(abs(nuplus))

disp('vminusone')
vpa(abs(vminusone))

disp('vminustwo')
vpa(abs(vminustwo))

disp('vplusone')
vpa(abs(vplusone))

disp('vplustwo')
vpa(abs(vplustwo))


%reduce epsilon_fixed (to make some room, when we symbolically check the
%constraints)
epsilon_fixed=sym((1/3)*epsilon_fixed);

%now: symbolically, check the constraints
%see https://www.mathworks.com/help/symbolic/check-symbolic-statements.html
disp('Symbolically, check the constraints')

%Rankine-Hugoniot conditions on the left interface
logical(numinus*(rho(1)-rho(2))==rho(1)*vminustwo-rho(2)*beta(1))
logical(numinus*(rho(1)*vminusone-rho(2)*alpha(1))==rho(1)*vminusone*vminustwo-rho(2)*delta(1))
logical(numinus*(rho(1)*vminustwo-rho(2)*beta(1))==rho(1)*vminustwo*vminustwo+rho(2)*gamma(1)+energyprime(1)*rho(1)*rho(1)-energyprime(2)*rho(2)*rho(2)-rho(2)*C(1)*(1/2))

%Rankine-Hugoniot conditions on the interface i

for i=1:(N-1)
logical(nu(i)*(rho(i+1)-rho(i+2))==rho(i+1)*beta(i)-rho(i+2)*beta(i+1))
logical(nu(i)*(rho(i+1)*alpha(i)-rho(i+2)*alpha(i+1))==rho(i+1)*delta(i)-rho(i+2)*delta(i+1))
logical(nu(i)*(rho(i+1)*beta(i)-rho(i+2)*beta(i+1))==-rho(i+1)*gamma(i)+rho(i+2)*gamma(i+1)+energyprime(i+1)*rho(i+1)*rho(i+1)-energyprime(i+2)*rho(i+2)*rho(i+2)+rho(i+1)*C(i)*(1/2)-rho(i+2)*C(i+1)*(1/2))
end


%Rankine-Hugoniot conditions on the right interface
logical(nuplus*(rho(N+1)-rho(N+2))==rho(N+1)*beta(N)-rho(N+2)*vplustwo)
logical(nuplus*(rho(N+1)*alpha(N)-rho(N+2)*vplusone)==rho(N+1)*delta(N)-rho(N+2)*vplusone*vplustwo)
logical(nuplus*(rho(N+1)*beta(N)-rho(N+2)*vplustwo)==-rho(N+1)*gamma(N)-rho(N+2)*vplustwo*vplustwo+energyprime(N+1)*rho(N+1)*rho(N+1)-energyprime(N+2)*rho(N+2)*rho(N+2)+rho(N+1)*C(N)*(1/2))

%subsolution conditions
for i=1:N
logical(alpha(i)*alpha(i)+beta(i)*beta(i)<=C(i)-epsilon_fixed)
logical(((1/2)*C(i)-alpha(i)*alpha(i)+gamma(i))*((1/2)*C(i)-beta(i)*beta(i)-gamma(i))-(delta(i)-alpha(i)*beta(i))^2>=epsilon_fixed)
end

%admissibility condition on the left interface
logical(numinus*(rho(1)*energy(1)-rho(2)*energy(2))+numinus*(rho(1)*(vminusone*vminusone+vminustwo*vminustwo)*(1/2)-rho(2)*C(1)*(1/2)) <= (rho(1)*energy(1)+energyprime(1)*rho(1)*rho(1))*vminustwo-(rho(2)*energy(2)+energyprime(2)*rho(2)*rho(2))*beta(1)+rho(1)*vminustwo*(vminusone*vminusone+vminustwo*vminustwo)*(1/2)-rho(2)*beta(1)*C(1)*(1/2)-epsilon_fixed)


%admissibility condition on inteface i

for i=1:(N-1)
logical(nu(i)*(rho(i+1)*energy(i+1)-rho(i+2)*energy(i+2))+nu(i)*(rho(i+1)*C(i)*(1/2)-rho(i+2)*(C(i+1))*(1/2))<=(rho(i+1)*energy(i+1)+energyprime(i+1)*rho(i+1)*rho(i+1))*beta(i)-(rho(i+2)*energy(i+2)+energyprime(i+2)*rho(i+2)*rho(i+2))*beta(i+1)+rho(i+1)*beta(i)*C(i)*(1/2)-rho(i+2)*beta(i+1)*C(i+1)*(1/2)-epsilon_fixed)
end

%admissibility condition on the right interface
logical(nuplus*(rho(N+1)*energy(N+1)-rho(N+2)*energy(N+2))+nuplus*(rho(N+1)*C(N)*(1/2)-rho(N+2)*(vplusone*vplusone+vplustwo*vplustwo)*(1/2))<=(rho(N+1)*energy(N+1)+energyprime(N+1)*rho(N+1)*rho(N+1))*beta(N)-(rho(N+2)*energy(N+2)+energyprime(N+2)*rho(N+2)*rho(N+2))*vplustwo+rho(N+1)*beta(N)*C(N)*(1/2)-rho(N+2)*vplustwo*(vplusone*vplusone+vplustwo*vplustwo)*(1/2)-epsilon_fixed)

%we want the fan parition to actually make sense
logical(numinus<=nu(1)-epsilon_fixed)
logical(nu(end)<=nuplus-epsilon_fixed)

for i=1:(N-2)
logical(nu(i)<=nu(i+1)-epsilon_fixed)
end

%we want a contact discontinuity
logical(vplustwo==vminustwo)
logical(rho(1)==rho(N+2))
logical(energy(1)==energy(N+2))
logical(energyprime(1)==energyprime(N+2))





%we want energy function to be strictly increasing and strictly convex

for i=1:(N+1)
for j=[1:(i-1) (i+1):(N+1)]
logical(energy(j)-energy(i)-energyprime(i)*(rho(j)-rho(i))>=epsilon_fixed)
end
end

%check positivity of density, energy derivative
for i=1:(N+2)
    logical(rho(i)>=epsilon_fixed)
    logical(energyprime(i)>=epsilon_fixed)
end

disp('Check the violations on the last two interfaces')

%Rankine-Hugoniot conditions on the interface i=N-1

i=(N-1)
vpa(nu(i)*(rho(i+1)-rho(i+2))-(rho(i+1)*beta(i)-rho(i+2)*beta(i+1)))
 
vpa(nu(i)*(rho(i+1)*alpha(i)-rho(i+2)*alpha(i+1))-(rho(i+1)*delta(i)-rho(i+2)*delta(i+1)))

vpa(nu(i)*(rho(i+1)*beta(i)-rho(i+2)*beta(i+1))-(-rho(i+1)*gamma(i)+rho(i+2)*gamma(i+1)+energyprime(i+1)*rho(i+1)*rho(i+1)-energyprime(i+2)*rho(i+2)*rho(i+2)+rho(i+1)*C(i)*(1/2)-rho(i+2)*C(i+1)*(1/2)))


%Rankine-Hugoniot conditions on the right interface
vpa(nuplus*(rho(N+1)-rho(N+2))-(rho(N+1)*beta(N)-rho(N+2)*vplustwo))
vpa(nuplus*(rho(N+1)*alpha(N)-rho(N+2)*vplusone)-(rho(N+1)*delta(N)-rho(N+2)*vplusone*vplustwo))
vpa(nuplus*(rho(N+1)*beta(N)-rho(N+2)*vplustwo)-(-rho(N+1)*gamma(N)-rho(N+2)*vplustwo*vplustwo+energyprime(N+1)*rho(N+1)*rho(N+1)-energyprime(N+2)*rho(N+2)*rho(N+2)+rho(N+1)*C(N)*(1/2)))

%show that these violations are always less than 10^(-11) in magnitude
UpperBound=str2sym('10^(-11)');

disp('The next 6 variable-precision numbers should all be negative')

%Rankine-Hugoniot conditions on the interface i=N-1

i=(N-1)
vpa(abs(nu(i)*(rho(i+1)-rho(i+2))-(rho(i+1)*beta(i)-rho(i+2)*beta(i+1)))-UpperBound)
 
vpa(abs(nu(i)*(rho(i+1)*alpha(i)-rho(i+2)*alpha(i+1))-(rho(i+1)*delta(i)-rho(i+2)*delta(i+1)))-UpperBound)

vpa(abs(nu(i)*(rho(i+1)*beta(i)-rho(i+2)*beta(i+1))-(-rho(i+1)*gamma(i)+rho(i+2)*gamma(i+1)+energyprime(i+1)*rho(i+1)*rho(i+1)-energyprime(i+2)*rho(i+2)*rho(i+2)+rho(i+1)*C(i)*(1/2)-rho(i+2)*C(i+1)*(1/2)))-UpperBound)


%Rankine-Hugoniot conditions on the right interface
vpa(abs(nuplus*(rho(N+1)-rho(N+2))-(rho(N+1)*beta(N)-rho(N+2)*vplustwo))-UpperBound)
vpa(abs(nuplus*(rho(N+1)*alpha(N)-rho(N+2)*vplusone)-(rho(N+1)*delta(N)-rho(N+2)*vplusone*vplustwo))-UpperBound)
vpa(abs(nuplus*(rho(N+1)*beta(N)-rho(N+2)*vplustwo)-(-rho(N+1)*gamma(N)-rho(N+2)*vplustwo*vplustwo+energyprime(N+1)*rho(N+1)*rho(N+1)-energyprime(N+2)*rho(N+2)*rho(N+2)+rho(N+1)*C(N)*(1/2)))-UpperBound)


