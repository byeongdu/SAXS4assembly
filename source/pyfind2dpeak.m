function ndt = pyfind2dpeak(dt)
% convert matlab array into numpy array
dtnp = py.numpy.array(dt);
% compute image analysis
dt = py.find2dpeaks.detect_peaks(dtnp);
% convert python array into matlab logical array.
ndt = logical(dt);