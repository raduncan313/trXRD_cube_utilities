function d = clear_masks(d)
    if isfield(d, 'masks')
        d = rmfield(d, 'masks');
    end
end