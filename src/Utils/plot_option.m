function plot_option()
    colors = [
    1 0 0;    % red
    1 1 0;    % yellow
    0 1 0;    % green 
    0 1 1;    % cyan
    0 0 1     % blue
    ];
    n_colors = 256; % resolution
    custom_map = interp1(linspace(0,1,5), colors, linspace(0,1,n_colors));

    custom_map = flipud(custom_map); % reverse

    colormap(custom_map);

end