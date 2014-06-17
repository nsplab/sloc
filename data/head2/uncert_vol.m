
diff_vol=[];
for casenum=3:9,
    casenum
    xs = [];
    ys = [];
    zs = [];
    vols=[];
    for num=99:100,
        num
        for i=1:num,
            i
            filename = ['/home/user/Documents/pubs/NER2013/datafiles/head.cost_at_grid_pts.10_20.noise_case' int2str(casenum) '.' int2str(i)];
            filename
            costs = dlmread(filename,' ');
            size(costs(:,7))
            costs_mtx = reshape(costs(:,7), 11,11,11);
            [val indx] = min(costs_mtx(:));
            [n m t] = ind2sub(size(costs_mtx),indx);
            xs(end+1) = n;
            ys(end+1) = m;
            zs(end+1) = t;
        end
        cov_mtx = cov([xs; ys; zs]');
        vol = 4/3*pi*sqrt(det(cov_mtx));
        vols(num-98) = vol;
    end
    diff_vol(end+1) = (vols(1) - vols(2))/vols(2);
end