function best_sol = gridsearch_LPSVD(t, x, Ls, rats)
    best_sol = struct()
    best_ssd = Inf;
    for L = Ls
        for rat = rats
            sol = LPSVD(t, x, L, rat);
            if sol.ssd < best_ssd
                best_sol = sol;
                best_ssd = sol.ssd;
            end
        end
    end
end