function angs = find_angles_from_rotation(fn, R, x0)
%     ff = @(angs) reshape(fn(angs(1), angs(2), angs(3)) - R, 9,1);
    ff = @(angs) reshape(fn(angs) - R, 9,1);
    options = optimoptions('fsolve','Display','none');
    [angs, fval, exitFlag] = fsolve(ff, x0, options);
    if exitFlag < 0
        error('No solution found.')
    end
end
