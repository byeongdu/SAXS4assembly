function Is = smearwl(Iexp, wlpar)
% calculation of the effect of wavelength distribtuion on SAXS pattern.
% distribution function is Gaussian.
% input parameter
%       Iexp : [q, I] - experimental data composed of q and intensity.
%       wlpar : [center, width] - experimental data composed of center, width.

wl = (wlpar(1)-3*wlpar(2)):wlpar(2)/10:(wlpar(1)+3*wlpar(2));
wld = normpdf(wl, wlpar(1), wlpar(2));wld = wld/sum(wld);

%figure
%hold on

%plot(Iexp(:,1), Iexp(:,2), 'b');

Is = [];
temp = zeros(size(Iexp(:,2)));
Is(:,1) = Iexp(:,1);

for i = 1:length(wl);
    %plot(Iexp(:,1)*wl(i)/wlpar(1), Iexp(:,2)*wld(i), 'r')
    temp = temp + spline(Iexp(:,1)*wl(i)/wlpar(1), Iexp(:,2)*wld(i), Is(:,1));
end    

Is(:,2) = temp;