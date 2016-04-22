% get_interval_overlap(foreground, background)
%
% Imagine an interval in the foreground moving over a stationary background
% interval. This function calculates the degree to which the foreground
% interval overlaps the background window, as a fraction of the width of
% the foreground window.
%
% For example:
%
% for i=1:10
%   disp(get_interval_overlap([i, i+4], [5, 9]))
% end
%
%     0
%     0
%     0.3333
%     0.6667
%     1
%     1
%     0.6667
%     0.3333
%     0
%     0
%
% Corresponds to:
%
% 0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16
%                |     bg    |
%                .           .
%  > |  fg01  | (0)          .
%       |  fg02  | (0)       .
%          |  fg03  | (0.3333)
%             |  fg04  | (0.6667)
%                |  fg05  | (1)
%                .  |  fg06  | (1)
%                .     |  fg07  | 0.6667)
%                .        |  fg08  | (0.333)
%                .           |  fg09  | (0)
%                .           .  |  fg10  | (0)
%                .           .
function fraction = get_interval_overlap(foreground, background)

    % Input validations
    assert(is_interval(foreground));
    assert(is_interval(background));

    % To track errors.
    % `overlap` should never end up as nan.
    overlap = NaN;

    % Find out which case we're in.
    if foreground(1) <= background(1)
        % `foreground` starts to the left of `background`
        if foreground(2) <= background(1)
            % `foreground` does not intersect with `background` by more
            % than a point.
            overlap = 0;
        else
            % There is some overlap. `foreground peaks over the left side
            % of background.
            if foreground(2) < background(2)
                % `foreground` hasn't covered all of `background`
                overlap = foreground(2) - background(1);
            else
                % `foreground` has covered all of `background`
                overlap = background(2) - background(1);
            end
        end
    else
        % `foreground` starts to the right of the left-edge of `background`
        if foreground(2) <= background(2)
            % `foreground` is entirely inside `background`
            overlap = foreground(2) - foreground(1);
        else
            % `foreground` ends to the right of `background`
            if foreground(1) >= background(2)
                % `foreground` does not intersect with `background` by more
                % than a point.
                overlap = 0;
            else
                % There is some overlap. `foreground` peaks over the right
                % side of `background`.
                overlap = background(2) - foreground(1);
            end
        end
    end
    
    % Checking that all cases were covered and that a value was assigned.
    assert(~isnan(overlap));
    
    % Avoid dividing by zero
    width = foreground(2) - foreground(1);
    if width == 0
        fraction = 0;
    else
        fraction = overlap / width;
    end

end

%% Subfunctions

function result = is_interval(X)
    result = (numel(X) == 2 && issorted(X));
end
