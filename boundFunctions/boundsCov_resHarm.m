function [l, u, l1, u1, l2, u2, flag] = boundsCov_resHarm(YT1, YT2, YC1, YC2, d1, d2)
%This function calculates the bounds on the fraction who benefit, with a
%baseline variable and restriction on harm.
%This code handles the case that the baseline variable divides the total 
%population into Subpopulation 1 and Subpopulation 2, as in the
%data application of the paper.
%It outputs the bounds on the fraction who benefit for the
%population (i.e., l and u). It also outputs the bounds on the fraction who
%benefit for each subpopulation (i.e., l1 and u1; l2 and u2).

%This code assumes the higher the outcome, the better.

%Inputs:
%YT1 is the vector of outcomes for those assigned to treatment and in 
%Subpopulation 1 
%YC1 is the vector of outcomes for those assigned to control and in
%Subpopulation 1
%YT2 is the vector of outcomes for those assigned to treatment and in 
%Subpopulation 2 
%YC2 is the vector of outcomes for those assigned to control and in
%Subpopulation 2
%d1 applies to the restriction on harm for Subpopulation 1 
%(i.e., Harm =  YC -YT <= d1 for Subpopulation 1)
%d2 applies to the restriction on harm for Subpopulation 2 
%(i.e., Harm = YC- YT <= d2 for Subpopulation 2)

%Outputs:
%l and u are the lower and upper bounds for the total population
%l1 and u1 are the lower and upper bounds for Subpopulation 1
%l2 and u2 are the lower and upper bounds for Subpopulation 2
%lflag and uflag are 0 if there was no error, and 1 otherwise
    [l1, u1, lflag1, uflag1] = boundsNoCov_mTmC_SPARSE_resHarm(YT1, YC1, d1);
    [l2, u2, lflag2, uflag2] = boundsNoCov_mTmC_SPARSE_resHarm(YT2, YC2, d2);
    weight1 = (length(YT1)+length(YC1))/(length(YT1)+length(YC1)+length(YT2)+length(YC2));
    l = l1 * weight1 + l2 * (1-weight1);
    u = u1 * weight1 + u2 * (1-weight1);
    if (lflag1==1) && (lflag2==1) && (uflag1==1) && (uflag2==1)
        flag = 0;
    else
        flag = 1;
    end
end