function En = RMS_err_metric(ini, new)
% ini : true object wavefront.
% new : calculated object wavefront.
% En : normalized RMS error metric
% Andrew M. Maiden and J. M. Rodenburg, Ultramicroscopy, 109, 1256 (2009).

% The parameter gamma allows 1) for the multiplication of object by a
% constant, and 2) for a constant phase offset.
gamma = sqrt(sum(ini.*conj(new), 'all')/sum(abs(new).^2, 'all'));
En = sum(abs(ini-gamma*new).^2, 'all')/sum(abs(ini).^2, 'all');
