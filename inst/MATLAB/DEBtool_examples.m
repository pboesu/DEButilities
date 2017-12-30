%DEBtool_M calculations that are used in the DEButilities R package to
%check for numerical equivalence of function outputs

%get_ue0 example
g = 6; k = 6; kap = .8; uHb = .001; vHb = uHb/ (1 - kap);
pars = [g, k, vHb];
[uE0, lb, info] = get_ue0(pars);

%beta0 example
betazero = beta0(0.1,0.2);