function plot_stress_Q8(model)
    % PLOT_STRESS:  PLOT THE STRESS CONTOURS
    
    %% Plot Stress contour
    plot_stress_contour_Q8(model, 1); % s_xx
    plot_stress_contour_Q8(model, 2); % s_yy
    plot_stress_contour_Q8(model, 3); % s_xy
    plot_stress_contour_Q8(model, 4); % von mises

end