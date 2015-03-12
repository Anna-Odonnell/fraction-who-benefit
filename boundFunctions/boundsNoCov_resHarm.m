function [l,u,lflag, uflag] = boundsNoCov_resHarm(YT, YC, d)
%This function calculates the bounds on the fraction who benefit, with a restriction on harm but no
%baseline variables.

%This code assumes the higher the outcome, the better.

%Inputs:
%YT is the vector of outcomes for those assigned to treatment
%YC is the vector of outcomes for those assigned to control
%d applies to the restriction on harm, which is Harm=YC - YT <= d

%Outputs:
%l and u are the lower and upper bounds
%lflag and uflag are flags if there were errors (1 if no error; see exitflag of Matlab linprog for other values)

YT_sort = unique(YT); %returns the values of the vector, sorted from smallest to largest
YC_sort = unique(YC);

mT = length(YT_sort);
mC = length(YC_sort);

%number of variables
varCount = mT * mC;
%number of constraints from marginal pmf's and restriction on harm
constCount = mT + mC + 1; 

%OBJECTIVE FUNCTION
f = zeros(mC, mT); %creates a mC x mT array of zeros
for r = 1:mC
    for c = 1:mT
        if YC_sort(r) < YT_sort(c)
            f(r,c) = 1; %1 if the element corresponds to benefit, 0 otherwise
        end
    end
end

f = reshape(f, varCount, 1); %reshapes f to a vector (varCount x 1)

%LOWER BOUND ON UNKNOWNS
lb = zeros(varCount,1);

%beq 
scaling = max([length(YT);length(YC)]); 
margPMF_C = zeros(mC, 1);
for i = 1:mC
    margPMF_C(i) = sum(YC == YC_sort(i));
end
margPMF_C = margPMF_C/length(YC)*scaling; %scaling the problem up might make the optimization more accurate

margPMF_T = zeros(mT, 1);
for i = 1:mT
    margPMF_T(i) = sum(YT == YT_sort(i));
end
margPMF_T = margPMF_T/length(YT)*scaling;

beq = [margPMF_T; margPMF_C; 0]; %0 is appended for restriction on harm

%Aeq

%i1,i2,j1,j2 are for marginal pmf's
%i3 and j3 are for restriction on harm
j3 = zeros(mC, mT);
for r = 1:mC
    for c = 1:mT
        if (YC_sort(r) - YT_sort(c)) > d
            j3(r,c) = 1;
        end
    end
end

j3 = reshape(j3, varCount, 1);
j3 = find(j3);
i3 = ones(length(j3),1)*constCount;


s = 1;

i1 = sort(repmat((1:mT)', mC, 1)); %1, 2,.., m_T each mC times
i2 = sort(repmat(((mT + 1):(mT+mC))',mT, 1)); %mT+1, mT+2, mT+mC each mT times 
i = [i1; i2; i3];


j1 = (1:1:varCount)';
j2 = transpose(reshape(j1, [mC, mT]));
j2 = reshape(j2, [varCount, 1]);
j = [j1;j2;j3];


Aeq = sparse(i,j,s, constCount, varCount);

%DOING THE OPTIMIZATION

% lower bound on fraction
[x,fval,exitflag] = linprog(f,[],[],Aeq,beq,lb);
% optimSoln_lower = reshape(x, m, m)
if exitflag == 1
    l = fval/scaling;
    lflag = exitflag;
else
    l = 0;
    lflag = exitflag;
end

% upper bound on fraction
f = f * -1;
[x,fval,exitflag] = linprog(f,[],[],Aeq,beq,lb);
% optimSoln_upper = reshape(x, m, m)
if exitflag == 1
    u = -fval/scaling;
    uflag = exitflag;
else
    u = 1;
    uflag = exitflag;
end



