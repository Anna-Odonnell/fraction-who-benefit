function [xl, xu, l, u, eps] = boundsNoCov_res(YT, YC, maxBen, maxHarm)


options = optimoptions(@linprog,'Algorithm','interior-point');

nT = length(YT);
nC = length(YC);

YT_sort = unique(YT); %returns the unique values of the vector, sorted from smallest to largest
YC_sort = unique(YC);

mT = length(YT_sort);
mC = length(YC_sort);

%number of pi_{i,j}'s
varCount = mT * mC;

scale = lcm(nT,nC);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate marginal cdf's% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cdf_C = zeros((mC-1), 1);
for i = 1:(mC-1)
    cdf_C(i) = sum(YC <= YC_sort(i));
end
cdf_C = cdf_C/nC; 

cdf_T = zeros((mT-1), 1);
for i = 1:(mT-1)
    cdf_T(i) = sum(YT <= YT_sort(i));
end
cdf_T = cdf_T/nT;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Which pi_{i,j}'s are affected by the restrictions?% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
res = zeros(mC, mT);

for r = 1:mC
    for c = 1:mT
        if (YT_sort(c) - YC_sort(r) > maxBen) || (YC_sort(r) - YT_sort(c) > maxHarm)
            res(r,c) = 1;
        end
    end
end

res = reshape(res, varCount, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define temp to be used when defining A% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp = zeros(varCount,mC+mT-2);

temp2 = zeros(mC, mT);
for i = 1:(mC-1)
    temp2(i,:) = 1;
    temp(:,i) = reshape(temp2, varCount, 1);
end

temp2 = zeros(mC, mT);
for i = 1:(mT-1)
    temp2(:,i) = 1;
    temp(:,mC-1+i) = reshape(temp2, varCount, 1);
end

temp = transpose(temp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate epsilon% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%OBJECTIVE FUNCTION
f = [zeros(varCount,1); 1];

%LOWER BOUND ON UNKNOWNS
lb = zeros((varCount+1),1);

%b
b = scale * [1; 0; 1; -1; cdf_C; cdf_T; -cdf_C; -cdf_T];

%A
A = zeros((2*mC+2*mT),(varCount+1)); %A has dimensions # of Constraints x # of Unknowns
A(1,:) = transpose(f);
A(2,:) = [transpose(res) 0];
A(3,:) = [ones(1,varCount) 0];
A(4,:) = -A(3,:);
A(5:end,end) = -1;
A(5:(2+mC+mT),1:(end-1)) = temp;
A((3+mC+mT):end,1:(end-1)) = -temp;

%Calculate epsilon
[~,fval,exitflag] = linprog(f,A,b,[],[],lb,[],[],options);
fval = fval/scale;
if exitflag == 1
    if fval < 10 ^ (-10)
        eps = 0;
    else 
        eps = fval;
    end
else
    eps = 10000; %there should be no problem with feasibility; this is to catch coding errors
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate l and u, using epsilon% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars temp temp2 res f b lb nC nT

%OBJECTIVE FUNCTION
f = zeros(mC, mT); %creates a mC x mT array of zeros
for r = 1:mC
    for c = 1:mT
        if YC_sort(r) < YT_sort(c)
            f(r,c) = 1; %1 if the element corresponds to benefit, 0 otherwise
        end
    end
end

f = reshape(f, varCount, 1); 

%LOWER BOUND ON UNKNOWNS
lb = zeros(varCount,1);

%b
b = scale * [0; 1; -1; (cdf_C + eps); (cdf_T + eps); (-cdf_C + eps); (-cdf_T + eps)];

%A 
A = A((2:end),1:(end-1));

% lower bound on fraction
[x,fval,exitflag] = linprog(f,A,b,[],[],lb,[],[],options);
xl = reshape(x,mC,mT)/scale;
if exitflag == 1
    l = fval/scale;
else
    l = 10000; %this is for debugging
end

% upper bound on fraction
f = f * -1;
[x,fval,exitflag] = linprog(f,A,b,[],[],lb,[],[],options);
xu = reshape(x,mC,mT)/scale;
if exitflag == 1
    u = -fval/scale;
else
    u = 10000; %this is for debugging
end
