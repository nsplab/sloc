vols=[];

% 10-20 electrodes
for casenum=3:9,
    casenum
    xs = [];
    ys = [];
    zs = [];
    for i=1:num,
        i
        filename = ['/home/user/Documents/pubs/NER2013/datafiles/head.cost_at_grid_pts.10_20.noise_case' int2str(casenum) '.' int2str(i)];
        filename
        costs = dlmread(filename,' ');
        size(costs(:,7))
        costs_mtx = reshape(costs(:,7), 11,11,11);
        [val indx] = min(costs_mtx(:));
        [n m t] = ind2sub(size(costs_mtx),indx);
        xs(end+1) = n-6;
        ys(end+1) = m-6;
        zs(end+1) = t-6;
    end
    %cov_mtx = cov([xs; ys; zs]');
    %vol = 4/3*pi*sqrt(det(cov_mtx));
    %vols(casenum-2,1) = vol*1.96;
    dist = sqrt(xs.^2+ys.^2+zs.^2);
    vols(casenum-2,1) = mean(dist);
end



% + EV1
for casenum=3:9,
    casenum
    xs = [];
    ys = [];
    zs = [];
    for i=1:num,
        i
        filename = ['/home/user/Documents/pubs/NER2013/datafiles/head.cost_at_grid_pts.10_20.ev1.noise_case' int2str(casenum) '.' int2str(i)];
        filename
        costs = dlmread(filename,' ');
        size(costs(:,7))
        costs_mtx = reshape(costs(:,7), 11,11,11);
        [val indx] = min(costs_mtx(:));
        [n m t] = ind2sub(size(costs_mtx),indx);
        xs(end+1) = n-6;
        ys(end+1) = m-6;
        zs(end+1) = t-6;
    end
        %cov_mtx = cov([xs; ys; zs]');
        %vol = 4/3*pi*sqrt(det(cov_mtx));
        %vols(casenum-2,2) = vol*1.96;
        dist = sqrt(xs.^2+ys.^2+zs.^2);
        vols(casenum-2,2) = mean(dist);
end

plot(vols);
ylabel('Mean distance from the true dipole location (cm)');
xlabel('SNR');
set(gca,'XTickLabel',{'5000','1000','500','100','50','10','5'})
