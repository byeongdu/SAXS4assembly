function [gk, chi2, err, gk0] = solve_densitymap(varargin)
% [gk, chi2] = solve_densitymap(hkl, Iq)
% [gk, chi2] = solve_densitymap(hkl, Iq, 'groupN', gn, 'histogram', ht,...
%               'display', true, 'trial_phase', ph)
% hkl : unique hkl (only one for the symmetry equivalent hkls)
%   ex: [2, 0, 0] only for {2, 0, 0} family.
hkl = varargin{1};

Iq = varargin{2};
willDisplay = false;
sg= [];
for i=1:numel(varargin)
    if ~ischar(varargin{i})
        continue
    end
    switch varargin{i}
        case 'groupN'
            groupN = varargin{i+1};
        case 'histogram'
            nh = varargin{i+1};
        case 'display'
            willDisplay = varargin{i+1};
        case 'trial_phase'
            trial_phase = varargin{i+1};
        case 'spacegroup'
            sg = varargin{i+1};
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta = 0.9;         % parameter in HIO algorithm
Iterations = 100;  % number of iterations, typically 200 iterations are enough
N = 51;             % matrix size.
Restart_allowance = 5;
Restart_threshold_error_value = 2E-4;
x = linspace(0, 1, N);%x(end) = [];
[frac_x, frac_y, frac_z] = ndgrid(x);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hkl_op = [];
Iq_op = [];
sym_groupN = [];
groupN_op = [];
if ~isempty(sg)
    for m = 1:size(hkl,1)
        hkls = symmetryoperate(sg, hkl(m,:), 'hkl', true);
        hkl_op = [hkl_op;hkls];
        Iq_op = [Iq_op;Iq(m)*ones(size(hkls, 1), 1)];
        sym_groupN = [sym_groupN;m*ones(size(hkls, 1), 1)];
        if exist('groupN', 'var')
            groupN_op = [groupN_op;groupN(m)*ones(size(hkls, 1), 1)];
        end
    end
    hkl = hkl_op;
    Iq = Iq_op;
    if exist('groupN', 'var')
        groupN = groupN_op;
    end
end


% 
% x = linspace(0, 1, N-1);x(end) = [];
% [frac_x, frac_y, frac_z] = ndgrid(x);





for i = 1:size(hkl, 1)
    if hkl(i,4)>0
        continue
    end
    current_HKL = hkl(i, :);
    % Check if the negated HKL exists in the array
    negated_HKL = -current_HKL;
    % Find rows that match the negated HKL
    match_idx = ismember(hkl, negated_HKL, 'rows');
    if any(match_idx)
        hkl(match_idx, 4) = i;
    end
end
hkl_Fridel_Pair = hkl(:,4);
hkl(:,4) = [];
h = hkl(:,1);
k = hkl(:,2);
l = hkl(:,3);

q = sqrt(sum(hkl.^2, 2));
N2 = size(Iq, 1);
nonobserved = Iq == -1;
if exist('trial_phase', 'var') % if initla phase does not exist...
    phase = trial_phase;
    %Iq(nonobserved) = trial_amp(nonobserved);
else
    % generate initial random phase 
    phase = zeros(size(h));
    phase(hkl_Fridel_Pair==0) = (2*rand(sum(hkl_Fridel_Pair==0),1) - 1)*pi;
    for i=1:numel(hkl_Fridel_Pair)
        if hkl_Fridel_Pair(i)>0
            phase(i) = -phase(hkl_Fridel_Pair(i));
        end
    end
    Iq(nonobserved) = 0;
end
Fq_amplitude_exp = sqrt(Iq);
[maxv, maxind] = max(Fq_amplitude_exp);
Fq_amplitude_exp = Fq_amplitude_exp/maxv;
Fq_trial = Fq_amplitude_exp.*exp(sqrt(-1)*phase);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Getting initial object distribution
ind_out = [];
%gk = fourier_synthesis(hkl, Fq_trial, N);
[tmpmatrix,~,~] = ndgrid(1:N+1,1:N+1,1:N+1);
if ~isempty(sg)
%    gk = symmetryoperate2voxel(gk, sg, 0);
    ind_out = symmetry_asymunit_position(tmpmatrix, sg);
    ind_out_index = unique(ind_out);
    tmp = ind_out==0;
    gk = fourier_synthesis(hkl, Fq_trial, ...
        frac_x(tmp), frac_y(tmp), frac_z(tmp));
    gk = fillgk(gk, ind_out, ind_out_index);
    gk = real(gk);
    gk(end,:,:)=[];
    gk(:,end,:)=[];
    gk(:,:,end)=[];
else
    gk = fourier_synthesis(hkl, Fq_trial, N);
    gk = real(gk);
end


er = ones(Iterations, 1)*nan;

if willDisplay
    figure(5);clf
end

histc_counter = 0;
ermin = 9999999;
ii = 1;
restartC = 0;
err_turnedup = 0;
gk_1 = gk;
while ii<Iterations
    %if error < preverr
    if err_turnedup
        if ermin > Restart_threshold_error_value
            if restartC < Restart_allowance
                % if first 30 is not good, try it again.
                % make random phases (but, phi_hkl = -phi_-h-k-l)
                phase = zeros(size(h));
                phase(hkl_Fridel_Pair==0) = (2*rand(sum(hkl_Fridel_Pair==0),1) - 1)*pi;
                for i=1:numel(hkl_Fridel_Pair)
                    if hkl_Fridel_Pair(i)>0
                        phase(i) = -phase(hkl_Fridel_Pair(i));
                    end
                end


                Fq_trial = Fq_amplitude_exp.*exp(sqrt(-1)*phase);
                %gk = real(fourier_synthesis(hkl, Fq_trial, N));
                if ~isempty(sg)
                    %gk = symmetryoperate2voxel(gk, sg, 0);
                    %gk = fillgk(gk, ind_out_index);
                    gk = real(fourier_synthesis(hkl, Fq_trial, ...
                        frac_x(tmp), frac_y(tmp), frac_z(tmp)));
                    gk = fillgk(gk, ind_out, ind_out_index);
                else
                    gk = real(fourier_synthesis(hkl, Fq_trial, N));
                end

                ii = 1;
                restartC = restartC + 1;
                err_turnedup = 0;
                %ermin = 9999999;
            end
        end
    end

%        numnh = numel(nh);
    if exist('nh', 'var')
        if (histc_counter > 50) & (Iterations - ii > 200) % Once in a while redistribute overlapped intensities.
            if isempty(sg)
                Fq_cnstr = calc_Fq(gk, h(:), k(:), l(:));
            else
                Fq_cnstr = calc_Fq({frac_x(tmp), frac_y(tmp), frac_z(tmp), gk}, h(:), k(:), l(:));
            end
            aFq = abs(Fq_cnstr);           % repartitioning experimental Fq.....
            if exist('groupN', 'var')
                gi = 1;
                ngi = groupN == gi;
                while sum(ngi) > 0
                    Fcnstr = aFq(ngi);
                    sumIq = sum(Fq_amplitude_exp(ngi).^2);
                    sumIqcnstr = sum(Fcnstr.^2);
                    Fq_amplitude_exp(ngi) = sqrt(Fcnstr.^2/sumIqcnstr*sumIq);
                    gi = gi + 1;
                    ngi = groupN == gi;
                end
                fprintf('histogram distribution done.\n');
            end
            if sum(nonobserved) > 0
                Fq_amplitude_exp(nonobserved) = aFq(nonobserved)/aFq(maxind)*max(Fq_amplitude_exp);
                fprintf('Assigning intensities to nonobserved ones.\n');
            end
            %
            histc_counter = 1;
        end
    end
    if isempty(sg)
        Fq_cnstr = calc_Fq(gk, h(:), k(:), l(:));
    else
        Fq_cnstr = calc_Fq({frac_x, frac_y, frac_z, gk}, h(:), k(:), l(:));
    end
    %error = RMS_err_metric(Fq_amplitude_exp, abs(Fq_cnstr));
    aFq = abs(Fq_cnstr);%aFq = aFq/aFq(maxind);
    % error = sum((Fq_amplitude_exp-aFq).^2);
    

    
    
    phase = angle(Fq_cnstr);
    % inversion...
    for i=1:numel(hkl_Fridel_Pair)
        if hkl_Fridel_Pair(i)>0
            phase(i) = -phase(hkl_Fridel_Pair(i));
        end
    end

    %gk_prime = real(fourier_synthesis(hkl, Fq_amplitude_exp.*exp(sqrt(-1)*phase), N));
    if ~isempty(sg)
        gk_prime = real(fourier_synthesis(hkl, Fq_amplitude_exp.*exp(sqrt(-1)*phase), ...
            frac_x(tmp), frac_y(tmp), frac_z(tmp)));
        gk_prime = fillgk(gk_prime, ind_out, ind_out_index);
        gk_prime(end,:,:)=[];
        gk_prime(:,end,:)=[];
        gk_prime(:,:,end)=[];
    else
        %gk = real(fourier_synthesis(hkl, Fq_trial, N));
        gk_prime = real(fourier_synthesis(hkl, Fq_amplitude_exp.*exp(sqrt(-1)*phase), N));
    end
    gk = gk_prime;
    % if ~isempty(sg)
    %     %if rem(ii,20)==0
    %         fprintf('Applying symmetry operation for %i\n', ii)
    %         %gk_prime = symmetryoperate2voxel(gk_prime, sg, 0);
    %         gk_prime = fillgk(gk_prime, ind_out_index);
    %     %end
    % end   

    %error = error_function_Fienup_RMS(Fq_amplitude_exp,aFq);
    error = chi_squared(Fq_amplitude_exp./q, aFq./q,0)/numel(Fq_amplitude_exp);
    %error = chi_squared(Fq_amplitude_exp./q, aFq./q/aFq(maxind),0)/numel(Fq_amplitude_exp);
    % if (error < 1E-10) & (ii ==1)
    %     error = NaN;
    % end
    
    er(ii) = error;

    if error < ermin
        gk0 = gk_prime;
        ermin = er(ii);
        %gk_1 = gk0;
    end
    
    if ii>15
        if er(ii) > er(ii-1)
            err_turnedup = 1;
        end
%        if abs(mean(diff(er(ii-10:ii)))) < 1E-10
        if abs(er(ii-1)-er(ii))/er(ii-1) < 1E-8 % when there is only 0.1% change.
            fprintf('Converged.\n');
            break
        end
        if error>ermin
            gk_1 = gk0;
        end
    end
    
    % Make a change to the density...
            % adjust histogram every time.
    if exist('nh', 'var')
        if ii>1
            gk_prime = histmatch2(gk_prime, nh(:));
        end
    end
    % 
    % if mean(gk_prime(:))<0
    %     gk_prime = -gk_prime;
    % end
    % % Error reduction.
    % gk = gk_prime;
    % t = gk<0;
    % gk(t) = 0;
    % 
  
    
    % HIO algorithm.
    maxgk = max(gk, [], 'all');
    mingk = min(gk, [], 'all');
    diffgk = maxgk - mingk;
    %t = abs(gk) > mingk + diffgk*0.3;
    t = abs(gk) < maxgk - diffgk*0.2;
    gk(t) = gk_1(t) - beta*gk_prime(t); %HIO



    if willDisplay
        figure(5); subplot(2, 2, 1)
        semilogy(er, 'ro-');
        title(sprintf('%i: itr, %d; Err = %0.2e / Err-min %0.2e.\n', restartC, ii, er(ii), ermin))
        figure(5);
        subplot(2,2,2);cla;
        semilogy(q, Fq_amplitude_exp./q, 'ro', q, aFq./q, 'bv')
        %subplot(1, 3, 3);cla;
        %histogram(gk_prime, numnh);

        %figure(6);
        subplot(2,2,3)
%        if isempty(sg)
            gk_tmp = gk_prime;
%        else
%            gk_tmp = fillgk(gk_prime, ind_out, ind_out_index);
%        end
        if ii==1
            im = imagesc(gk_tmp(:,:,1));
            axis image;
        else
            set(im, 'CData', gk_tmp(:,:,1));
        end
        subplot(2,2,4)
        histogram(gk_prime, 100)
        title(sprintf('Itr: %d, Mean density = %0.3e.\n', ii, mean(gk_prime(:))))
    end
    
    gk_1 = gk_prime;

    histc_counter = histc_counter +1;
    drawnow
    ii = ii + 1;
end
err = er;
gk = gk_prime;
gk = check_mirror(gk);

%figure(6);
if willDisplay
    set(im, 'CData', gk(:,:,1));
end
aFq = abs(Fq_cnstr);


Iqc = abs(aFq).^2;
% Iqc = Iqc/max(Iqc);
% Iq = Iq/max(Iq);

t = Iq<1E-10;
Iq(t) = [];
Iqc(t) = [];

chi2 = chi_squared(Iq,Iqc,0);      

end

function map = fillgk(gk, ind_out, ind_out_index)
    N = round(numel(ind_out)^(1/3));
    map = zeros(size(ind_out));
    map(ind_out==0) = gk;
    for index_sym=2:numel(ind_out_index)
        idxfrom = ind_out_index(index_sym);
        map(ind_out==idxfrom)=map(ind_out_index(index_sym));
    end
    map = reshape(map, [N,N,N]);
end
