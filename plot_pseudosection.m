function plot_pseudosection(ele, ABMN, d)
    % Function that takes BERT data struct and plots pseudosection.
    %
    % SYNTAX
    %   plot_pseudosection(ele, ABMN, d)
    %
    % INPUT PARAMETER
    %   ele  ... Matrix [n_ele x 2] of electrode coordinates.
    %   ABMN ... Matrix [n_obs x 4] indexing ele for survey configuration
    %            definitions 
    %   d    ... Vector [n_obs x 1] of respective data
    %
    % Sascha Weit 2022

    % Generate grid coordinates.
    n_ele = max(ABMN(:));
    ref_pos = mean([ele(ABMN(:,1),1), ele(ABMN(:,4),1)].').';
    separation = abs(ABMN(:, 3) - ABMN(:, 2));
    diff = get_diffs();

    % Put together patch cells.
    X_patch = [ref_pos-diff(:,1)/2, ref_pos-diff(:,1)/2, ref_pos+diff(:,2)/2, ref_pos+diff(:,2)/2];
    Y_patch = [separation - 0.5, separation + 0.5, separation + 0.5, separation - 0.5];

    % Plot.
    patch(X_patch.', Y_patch.', d(:).');
    axis equal;
    xlim([min(X_patch(:)), max(X_patch(:))]);
    ylim([min(Y_patch(:)), max(Y_patch(:))]);
    set(gca, 'Ydir', 'reverse');
    xlabel('profile meter');
    ylabel('separation');
    colorbar();

    % Helper.
    function diffs = get_diffs()
        % Calculates horizontal cell sizes.
        diffs = zeros(length(ABMN), 2);
        ele_starts = ABMN(:,1) == 1;
        ele_ends =  ABMN(:,4) == n_ele;
        [~, i_next] = ismember(ABMN+1, ABMN,'rows');
        [~, i_prev] = ismember(ABMN-1, ABMN,'rows');
        i_next(i_next == 0) = [];
        i_prev(i_prev == 0) = [];
        diffs(~ele_starts, 1) = abs(ref_pos(~ele_starts) - ref_pos(i_prev));
        diffs(~ele_ends, 2) = abs(ref_pos(~ele_ends) - ref_pos(i_next));
        diffs(ele_ends & ele_starts, :) = 1;
        diffs(ele_starts, 1) = abs(ele(ABMN(ele_starts,1),1) - ele(ABMN(ele_starts,4),1))./(separation(ele_starts)+2);
        diffs(ele_ends, 2) = abs(ele(ABMN(ele_ends,1),1) - ele(ABMN(ele_ends,4),1))./(separation(ele_ends)+2);
    end
end
