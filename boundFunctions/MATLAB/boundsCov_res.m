function [l, u, l1, u1, eps1, l2, u2, eps2] = boundsCov_res(YT1, YT2, YC1, YC2, maxBen, maxHarm)
    [~, ~, l1, u1, eps1] = boundsNoCov_res(YT1, YC1, maxBen, maxHarm);
    [~, ~, l2, u2, eps2] = boundsNoCov_res(YT2, YC2, maxBen, maxHarm);
    weight1 = (length(YT1)+length(YC1))/(length(YT1)+length(YC1)+length(YT2)+length(YC2));
    l = l1 * weight1 + l2 * (1-weight1);
    u = u1 * weight1 + u2 * (1-weight1);
end