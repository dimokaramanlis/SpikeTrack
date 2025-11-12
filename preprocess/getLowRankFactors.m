function [U, V] = getLowRankFactors(template, N)
% Returns factors U [Nchan, N] and V [Nt, N]
% such that template = U * V'
[aa, bb, cc] = svds(double(template), N);
U = single(aa * bb);
V = single(cc); % svds already gives V' as cc', so V is cc
end