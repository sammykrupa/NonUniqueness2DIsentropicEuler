% This code searches for solutions to the algebraic equations making up the
% condition of being a subsolution fan for
% the 2-D isentropic Euler equations.
% The Riemann initial data are random values, which are chosen each time
% the program is run.
% This code is associated with the paper "Contact discontinuities for 2-D isentropic Euler are unique in 1-D but wildly non-unique otherwise"
% by Sam G. Krupa and László Székelyhidi Jr.


TotalSuccess=0; %keep track of how many times the numerical solver converges successfully


%loop 100 times (with different random Riemann problems)
for TrialLoop=1:100

%print our loop number
TrialLoop

clearvars -except TotalSuccess TrialLoop %start from a blank slate on each loop

%number of pizza slices in the fan partition
N=3; %This code only works with N=3

%the epsilon determines the strictness in strict inequalities
epsilon_fixed=.0001;

%Setup the numerical solver to find a subsolution fan, and run the solver

%setup the optimization variables   
alpha=optimvar('alpha',N);
beta=optimvar('beta',N);
gamma=optimvar('gamma',N);
delta=optimvar('delta',N);


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

%hold v_{-} and v_{+} values fixed: random perturbations of the values used
%for the contact discontinuity non-uniqueness example in the paper.
%see https://www.mathworks.com/help/matlab/ref/rand.html

%choose how large to make perturbations
perturbationSize=1;

vminusone=-4098844157247653/70368744177664+perturbationSize*(rand-.5);
vplusone=3603433899522037/562949953421312+perturbationSize*(rand-.5);
vminustwo=-996118042660627/70368744177664+perturbationSize*(rand-.5);
vplustwo=-996118042660627/70368744177664+perturbationSize*(rand-.5);

%hold rho_{+} and rho_{-} values fixed: random perturbations of the values used
%for the contact discontinuity non-uniqueness example in the paper.

rho1=2708112612978501/281474976710656+perturbationSize*(rand-.5);
rhoN2=2708112612978501/281474976710656+perturbationSize*(rand-.5);



%define the optimization problem

prob=optimproblem;

%define the constraints 

%Rankine-Hugoniot conditions on the left interface
prob.Constraints.RH_left_1=numinus*(rho1-rho(2))==rho1*vminustwo-rho(2)*beta(1);
prob.Constraints.RH_left_2=numinus*(rho1*vminusone-rho(2)*alpha(1))==rho1*vminusone*vminustwo-rho(2)*delta(1);
prob.Constraints.RH_left_3=numinus*(rho1*vminustwo-rho(2)*beta(1))==rho1*vminustwo*vminustwo+rho(2)*gamma(1)+energyprime(1)*rho1*rho1-energyprime(2)*rho(2)*rho(2)-rho(2)*C(1)*(1/2);

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
prob.Constraints.RH_right_1=nuplus*(rho(N+1)-rhoN2)==rho(N+1)*beta(N)-rhoN2*vplustwo;
prob.Constraints.RH_right_2=nuplus*(rho(N+1)*alpha(N)-rhoN2*vplusone)==rho(N+1)*delta(N)-rhoN2*vplusone*vplustwo;
prob.Constraints.RH_right_3=nuplus*(rho(N+1)*beta(N)-rhoN2*vplustwo)==-rho(N+1)*gamma(N)-rhoN2*vplustwo*vplustwo+energyprime(N+1)*rho(N+1)*rho(N+1)-energyprime(N+2)*rhoN2*rhoN2+rho(N+1)*C(N)*(1/2);

%subsolution conditions
for i=1:N
constraint = strcat('sub1','_',num2str(i));
prob.Constraints.(constraint)=alpha(i)*alpha(i)+beta(i)*beta(i)<=C(i)-epsilon_fixed;

constraint = strcat('sub2','_',num2str(i));
prob.Constraints.(constraint)=((1/2)*C(i)-alpha(i)*alpha(i)+gamma(i))*((1/2)*C(i)-beta(i)*beta(i)-gamma(i))-(delta(i)-alpha(i)*beta(i))^2>=epsilon_fixed;
end

%admissibility condition on the left interface
prob.Constraints.adm_left=numinus*(rho1*energy(1)-rho(2)*energy(2))+numinus*(rho1*(vminusone*vminusone+vminustwo*vminustwo)*(1/2)-rho(2)*C(1)*(1/2)) <= (rho1*energy(1)+energyprime(1)*rho1*rho1)*vminustwo-(rho(2)*energy(2)+energyprime(2)*rho(2)*rho(2))*beta(1)+rho1*vminustwo*(vminusone*vminusone+vminustwo*vminustwo)*(1/2)-rho(2)*beta(1)*C(1)*(1/2);


%admissibility condition on inteface i

for i=1:(N-1)
 constraint = strcat('adm','_',num2str(i));
 prob.Constraints.(constraint)=nu(i)*(rho(i+1)*energy(i+1)-rho(i+2)*energy(i+2))+nu(i)*(rho(i+1)*C(i)*(1/2)-rho(i+2)*(C(i+1))*(1/2))<=(rho(i+1)*energy(i+1)+energyprime(i+1)*rho(i+1)*rho(i+1))*beta(i)-(rho(i+2)*energy(i+2)+energyprime(i+2)*rho(i+2)*rho(i+2))*beta(i+1)+rho(i+1)*beta(i)*C(i)*(1/2)-rho(i+2)*beta(i+1)*C(i+1)*(1/2);
end

%admissibility condition on the right interface
prob.Constraints.adm_right=nuplus*(rho(N+1)*energy(N+1)-rhoN2*energy(N+2))+nuplus*(rho(N+1)*C(N)*(1/2)-rhoN2*(vplusone*vplusone+vplustwo*vplustwo)*(1/2))<=(rho(N+1)*energy(N+1)+energyprime(N+1)*rho(N+1)*rho(N+1))*beta(N)-(rhoN2*energy(N+2)+energyprime(N+2)*rhoN2*rhoN2)*vplustwo+rho(N+1)*beta(N)*C(N)*(1/2)-rhoN2*vplustwo*(vplusone*vplusone+vplustwo*vplustwo)*(1/2);

%we want the fan parition to actually make sense
prob.Constraints.fan_condition_First=numinus<=nu(1)-epsilon_fixed;
prob.Constraints.fan_condition_Last=nu(end)<=nuplus-epsilon_fixed;

for i=1:(N-2)
constraint = strcat('fan_condition','_',num2str(i));
prob.Constraints.(constraint)=nu(i)<=nu(i+1)-epsilon_fixed;
end










%we want energy function to be strictly increasing and strictly convex
%(note the lower bound on the optimization variable energyprime which is defined
%above)
for i=1:(N+2)
for j=[1:(i-1) (i+1):(N+2)]
    if i==1
        if j~=(N+2) & j~=1
            constraint = strcat('ConstrConvexityCondition','_',num2str(i),'_',num2str(j));
            prob.Constraints.(constraint)=energy(j)-energy(i)-energyprime(i)*(rho(j)-rho1)>=epsilon_fixed*(rho(j)-rho1)^2;
        end
    end

     if i==(N+2)
        if j~=(N+2) & j~=1
            constraint = strcat('ConstrConvexityCondition','_',num2str(i),'_',num2str(j));
            prob.Constraints.(constraint)=energy(j)-energy(i)-energyprime(i)*(rho(j)-rhoN2)>=epsilon_fixed*(rho(j)-rhoN2)^2;
        end
     end

     if i==1
        if j==(N+2)
            constraint = strcat('ConstrConvexityCondition','_',num2str(i),'_',num2str(j));
            prob.Constraints.(constraint)=energy(j)-energy(i)-energyprime(i)*(rhoN2-rho1)>=epsilon_fixed*(rhoN2-rho1)^2;
        end
     end
     if i==(N+2)
        if j==1
            constraint = strcat('ConstrConvexityCondition','_',num2str(i),'_',num2str(j));
            prob.Constraints.(constraint)=energy(j)-energy(i)-energyprime(i)*(rho1-rhoN2)>=epsilon_fixed*(rho1-rhoN2)^2;
        end
     end

      if j==1
        if i~=(N+2) & i~=1
            constraint = strcat('ConstrConvexityCondition','_',num2str(i),'_',num2str(j));
            prob.Constraints.(constraint)=energy(j)-energy(i)-energyprime(i)*(rho1-rho(i))>=epsilon_fixed*(rho1-rho(i))^2;
        end
      end

      if j==N+2
        if i~=(N+2) & i~=1
            constraint = strcat('ConstrConvexityCondition','_',num2str(i),'_',num2str(j));
            prob.Constraints.(constraint)=energy(j)-energy(i)-energyprime(i)*(rhoN2-rho(i))>=epsilon_fixed*(rhoN2-rho(i))^2;
        end
     end

     if i~=1 & i~=(N+2) & j~=1 & j~=(N+2)
constraint = strcat('ConstrConvexityCondition','_',num2str(i),'_',num2str(j));
prob.Constraints.(constraint)=energy(j)-energy(i)-energyprime(i)*(rho(j)-rho(i))>=epsilon_fixed*(rho(j)-rho(i))^2;
     end
end
end


%set the initial point for the numerical solver: Here, we include the exact values of the constants given in the appendix
%of the paper, which correspond to contact discontinuity initial data.

data.alpha(1)=-8177336068870495/140737488355328;
data.alpha(2)=-4833381446756075/562949953421312;
data.alpha(3)=3121572020159473/562949953421312;
data.beta(1)=-2536561643647751/140737488355328;
data.beta(2)=-1114286601116939/70368744177664;
data.beta(3)=-6219197795695073/562949953421312;
data.gamma(1)=841617150350781/549755813888;
data.gamma(2)=-2850833975067331/17592186044416;
data.gamma(3)=-8954832877447991/140737488355328;
data.delta(1)=28872176135415855785280523654056524019908546591/27627619078169805047324756605549438692229120;
data.delta(2)=-867454945412067709200232997995952542374584537720982074207241594308599438377/32748846874784971211058574285222379723549486466626634273206947431338475520;
data.delta(3)=-2871256077954219/35184372088832;
data.rho(1)=rho1;
data.rho(2)=6811063536043807/562949953421312;
data.rho(3)=2057060350258899/562949953421312;
data.rho(4)=3062207031116133/281474976710656;
data.rho(5)=rhoN2;
data.energy(1)=-5041529442624971/2199023255552;
data.energy(2)=-5015532875605977/2199023255552;
data.energy(3)=-5073206593829053/2199023255552;
data.energy(4)=-2515400677054201/1099511627776;
data.energy(5)=-5041529442624971/2199023255552;
data.energyprime(1)=1676289422169645/562949953421312;
data.energyprime(2)=60006068216738166351756926651195471316209797920723479281182349/9106752347169134708406278810788569756749285046122185710632960;
data.energyprime(3)=1831278218949756891087541424381121781924211556364089293813566516592411003/1359848686641341079894728790850171965963388817622134940140753492461486080;
data.energyprime(4)=5400383921383283/1125899906842624;
data.energyprime(5)=1676289422169645/562949953421312;
data.numinus=-6486283176597739958874052307549/196306040423407104692364247040;
data.nu(1)=-1153852086001065889673487658885/60824224363690518566334889984;
data.nu(2)=-4856156003780791/562949953421312;
data.nuplus=7162856387903725/562949953421312;
data.C(1)=510415269881361/137438953472;
data.C(2)=1515879700153707/2199023255552;
data.C(3)=1855257252703141/8796093022208;




%solve!
disp('Start the solver')
[sol,fval,exitflag,output] = solve(prob,data,'Options',optimoptions(prob,'MaxFunctionEvaluations',20000000,'MaxIterations',2000,'ConstraintTolerance',.00000000000001,'Display','iter'));
%remark that even if the solver converges to an infeasible point, the
%solution may still be valid !

%save the original solution (before we start modifying it). Just in case
%we want to examine it later.
original_sol=sol;


%convert numerical solution from solver to symbolic

alpha=sym(sol.alpha,'f');
beta=sym(sol.beta,'f');


vminusone=sym(vminusone,'f');
vplusone=sym(vplusone,'f');
vminustwo=sym(vminustwo,'f');
vplustwo=sym(vplustwo,'f');
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

%we have fixed values for rho_{-} and rho_{+}
rho(1)=sym(rho1,'f');
rho(N+2)=sym(rhoN2,'f');



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





%calculate D\Gamma at the point (\hat\alpha_N,\hat\beta_N,\hat\delta_N,\hat\rho_N,\hat\nu_+,\hat\nu_{N-1})

pprime = 2*rho(N+1)*energyprime(N+1); %p'(\rho_{N})
DGamma=sym(zeros(6,6));
DGamma=[0 -rho(N+1) 0 -beta(N)+nu(N-1) 0 -(rho(N)-rho(N+1));rho(N+1)*nu(N-1) 0 -rho(N+1) -delta(N)+alpha(N)*nu(N-1) 0 -(rho(N)*alpha(N-1)-rho(N+1)*alpha(N));0 nu(N-1)*rho(N+1) 0 gamma(N)-pprime-.5*C(N)+beta(N)*nu(N-1) 0 -(rho(N)*beta(N-1)-rho(N+1)*beta(N));0 rho(N+1) 0 beta(N)-nuplus -(rho(N+1)-rho(N+2)) 0;-nuplus*rho(N+1) 0 rho(N+1) delta(N)-nuplus*alpha(N) -(rho(N+1)*alpha(N)-rho(N+2)*vplusone) 0;0 -nuplus*rho(N+1) 0 -gamma(N)+pprime+.5*C(N)-nuplus*beta(N) -(rho(N+1)*beta(N)-rho(N+2)*vplustwo) 0];

disp('Singular values of DGamma')
vpa(svd(DGamma))

%calculate some maximum of various quantities to figure out how
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

%keep track of the responses when checking the constraints
tally=[];

%Rankine-Hugoniot conditions on the left interface
tally(end+1)=logical(numinus*(rho(1)-rho(2))==rho(1)*vminustwo-rho(2)*beta(1))
tally(end+1)=logical(numinus*(rho(1)*vminusone-rho(2)*alpha(1))==rho(1)*vminusone*vminustwo-rho(2)*delta(1))
tally(end+1)=logical(numinus*(rho(1)*vminustwo-rho(2)*beta(1))==rho(1)*vminustwo*vminustwo+rho(2)*gamma(1)+energyprime(1)*rho(1)*rho(1)-energyprime(2)*rho(2)*rho(2)-rho(2)*C(1)*(1/2))

%Rankine-Hugoniot conditions on the interface i

for i=1:(N-1)
tally(end+1)=logical(nu(i)*(rho(i+1)-rho(i+2))==rho(i+1)*beta(i)-rho(i+2)*beta(i+1))
tally(end+1)=logical(nu(i)*(rho(i+1)*alpha(i)-rho(i+2)*alpha(i+1))==rho(i+1)*delta(i)-rho(i+2)*delta(i+1))
tally(end+1)=logical(nu(i)*(rho(i+1)*beta(i)-rho(i+2)*beta(i+1))==-rho(i+1)*gamma(i)+rho(i+2)*gamma(i+1)+energyprime(i+1)*rho(i+1)*rho(i+1)-energyprime(i+2)*rho(i+2)*rho(i+2)+rho(i+1)*C(i)*(1/2)-rho(i+2)*C(i+1)*(1/2))
end


%Rankine-Hugoniot conditions on the right interface
tally(end+1)=logical(nuplus*(rho(N+1)-rho(N+2))==rho(N+1)*beta(N)-rho(N+2)*vplustwo)
tally(end+1)=logical(nuplus*(rho(N+1)*alpha(N)-rho(N+2)*vplusone)==rho(N+1)*delta(N)-rho(N+2)*vplusone*vplustwo)
tally(end+1)=logical(nuplus*(rho(N+1)*beta(N)-rho(N+2)*vplustwo)==-rho(N+1)*gamma(N)-rho(N+2)*vplustwo*vplustwo+energyprime(N+1)*rho(N+1)*rho(N+1)-energyprime(N+2)*rho(N+2)*rho(N+2)+rho(N+1)*C(N)*(1/2))

%subsolution conditions
for i=1:N
tally(end+1)=logical(alpha(i)*alpha(i)+beta(i)*beta(i)<=C(i)-epsilon_fixed)
tally(end+1)=logical(((1/2)*C(i)-alpha(i)*alpha(i)+gamma(i))*((1/2)*C(i)-beta(i)*beta(i)-gamma(i))-(delta(i)-alpha(i)*beta(i))^2>=epsilon_fixed)
end

%admissibility condition on the left interface
tally(end+1)=logical(numinus*(rho(1)*energy(1)-rho(2)*energy(2))+numinus*(rho(1)*(vminusone*vminusone+vminustwo*vminustwo)*(1/2)-rho(2)*C(1)*(1/2)) <= (rho(1)*energy(1)+energyprime(1)*rho(1)*rho(1))*vminustwo-(rho(2)*energy(2)+energyprime(2)*rho(2)*rho(2))*beta(1)+rho(1)*vminustwo*(vminusone*vminusone+vminustwo*vminustwo)*(1/2)-rho(2)*beta(1)*C(1)*(1/2)-epsilon_fixed)


%admissibility condition on inteface i

for i=1:(N-1)
tally(end+1)=logical(nu(i)*(rho(i+1)*energy(i+1)-rho(i+2)*energy(i+2))+nu(i)*(rho(i+1)*C(i)*(1/2)-rho(i+2)*(C(i+1))*(1/2))<=(rho(i+1)*energy(i+1)+energyprime(i+1)*rho(i+1)*rho(i+1))*beta(i)-(rho(i+2)*energy(i+2)+energyprime(i+2)*rho(i+2)*rho(i+2))*beta(i+1)+rho(i+1)*beta(i)*C(i)*(1/2)-rho(i+2)*beta(i+1)*C(i+1)*(1/2)-epsilon_fixed)
end

%admissibility condition on the right interface
tally(end+1)=logical(nuplus*(rho(N+1)*energy(N+1)-rho(N+2)*energy(N+2))+nuplus*(rho(N+1)*C(N)*(1/2)-rho(N+2)*(vplusone*vplusone+vplustwo*vplustwo)*(1/2))<=(rho(N+1)*energy(N+1)+energyprime(N+1)*rho(N+1)*rho(N+1))*beta(N)-(rho(N+2)*energy(N+2)+energyprime(N+2)*rho(N+2)*rho(N+2))*vplustwo+rho(N+1)*beta(N)*C(N)*(1/2)-rho(N+2)*vplustwo*(vplusone*vplusone+vplustwo*vplustwo)*(1/2)-epsilon_fixed)

%we want the fan parition to actually make sense
tally(end+1)=logical(numinus<=nu(1)-epsilon_fixed)
tally(end+1)=logical(nu(end)<=nuplus-epsilon_fixed)

for i=1:(N-2)
tally(end+1)=logical(nu(i)<=nu(i+1)-epsilon_fixed)
end




%we want energy function to be strictly increasing and strictly convex

for i=1:(N+2)
for j=[1:(i-1) (i+1):(N+2)]
tally(end+1)=logical(energy(j)-energy(i)-energyprime(i)*(rho(j)-rho(i))>=epsilon_fixed*(rho(j)-rho(i))^2)
end
end

%check positivity of density, energy derivative
for i=1:(N+2)
    tally(end+1)=logical(rho(i)>=epsilon_fixed)
    tally(end+1)=logical(energyprime(i)>=epsilon_fixed)
end


%check symbolically that the violation is not too large on the last two interfaces

UpperBound=str2sym('10^(-9)');

%Rankine-Hugoniot conditions on the interface i=N-1

i=(N-1)

tally(end+1)=logical(abs(nu(i)*(rho(i+1)-rho(i+2))-(rho(i+1)*beta(i)-rho(i+2)*beta(i+1)))-UpperBound<=0)
 
tally(end+1)=logical(abs(nu(i)*(rho(i+1)*alpha(i)-rho(i+2)*alpha(i+1))-(rho(i+1)*delta(i)-rho(i+2)*delta(i+1)))-UpperBound<=0)

tally(end+1)=logical(abs(nu(i)*(rho(i+1)*beta(i)-rho(i+2)*beta(i+1))-(-rho(i+1)*gamma(i)+rho(i+2)*gamma(i+1)+energyprime(i+1)*rho(i+1)*rho(i+1)-energyprime(i+2)*rho(i+2)*rho(i+2)+rho(i+1)*C(i)*(1/2)-rho(i+2)*C(i+1)*(1/2)))-UpperBound<=0)


%Rankine-Hugoniot conditions on the right interface
tally(end+1)=logical(abs(nuplus*(rho(N+1)-rho(N+2))-(rho(N+1)*beta(N)-rho(N+2)*vplustwo))-UpperBound<=0)
tally(end+1)=logical(abs(nuplus*(rho(N+1)*alpha(N)-rho(N+2)*vplusone)-(rho(N+1)*delta(N)-rho(N+2)*vplusone*vplustwo))-UpperBound<=0)
tally(end+1)=logical(abs(nuplus*(rho(N+1)*beta(N)-rho(N+2)*vplustwo)-(-rho(N+1)*gamma(N)-rho(N+2)*vplustwo*vplustwo+energyprime(N+1)*rho(N+1)*rho(N+1)-energyprime(N+2)*rho(N+2)*rho(N+2)+rho(N+1)*C(N)*(1/2)))-UpperBound<=0)

ViolatedConditions=length(tally)-sum(tally);

if ViolatedConditions==6 %solution found
TotalSuccess=TotalSuccess+1;
end

end

disp('Total successes:')
TotalSuccess