function model = plot_stress_T3(model)
    % PLOT_STRESS:  PLOT THE STRESS CONTOURS
    
     %% Plot Stress contour
     plot_stress_contour_T3(model, 1); % s_xx
     plot_stress_contour_T3(model, 2); % s_yy
     plot_stress_contour_T3(model, 3); % s_xy
     plot_stress_contour_T3(model, 4); % von mises


end