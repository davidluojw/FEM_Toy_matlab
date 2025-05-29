function [L2_error, H1_error, log_L2_error, log_H1_error, L2_rate_error, H1_rate_error] = convergence_analysis_Q8()
    figure;
    h = [1/4,1/8, 1/12, 1/16, 1/20];
    log_h = log(h);
    
    
    
    % =========================================================================
    % deg=1
    L2_error = [0.004828426207519, 5.422460501647681e-04, 1.574966811557829e-04, 6.596494341468588e-05, 3.365886498163686e-05];
    log_L2_error = log(L2_error);
    H1_error = [0.031889839015168, 0.005989080912743, 0.002504472615763, 0.001379592470969, 8.746777132887221e-04];
    log_H1_error = log(H1_error);
    
    % =========================================================================
    
    
    
    h1 = plot(log_h, log_L2_error,'b-', 'LineWidth', 2);
    hold on;
    h2 = plot(log_h, log_H1_error,'r-', 'LineWidth', 2);
    xlabel("log-h");
    ylabel("log-error");
    legend('L2-error', 'H1-error','Location', 'Best', 'FontSize', 14, 'Box', 'on');
    
    % rate of error
    L2_rate_error = zeros(length(log_h)-1,1);
    H1_rate_error = zeros(length(log_h)-1,1);
    for ii = 1:length(log_h)-1
        L2_rate_error(ii) = (log_L2_error(ii+1) - log_L2_error(ii))/(log_h(ii+1)-log_h(ii));
        H1_rate_error(ii) = (log_H1_error(ii+1) - log_H1_error(ii))/(log_h(ii+1)-log_h(ii));
    end

    % add slope text
    for seg = 1:length(log_h)-1
        x_mid = mean(log_h(seg:seg+1));
        y_mid_L2 = mean(log_L2_error(seg:seg+1)) + 0.2;
        y_mid_H1 = mean(log_H1_error(seg:seg+1)) + 0.2;
        text(x_mid, y_mid_L2, sprintf('Slope: %.4f',L2_rate_error(seg)),...
            'Color', 'b',...
            'HorizontalAlignment', 'center',...
            'VerticalAlignment', 'bottom',...
            'FontSize', 12);
        text(x_mid, y_mid_H1, sprintf('Slope: %.4f',H1_rate_error(seg)),...
            'Color', 'r',...
            'HorizontalAlignment', 'center',...
            'VerticalAlignment', 'bottom',...
            'FontSize', 12);
    end
    legend('show'); 
    hold off;

    fprintf('\n                           Convergence\n');
    fprintf('-------------------------------------------------------------------------------\n');
    fprintf('\tgrid\teL2\t\teH1\t\tlog(eL1)\tlog(eH1)\n');
    
    for i = 1:length(log_h)
        eL2 = L2_error(i);
        eH1 = H1_error(i);
        logeL2 = log_L2_error(i);  
        logeH1 = log_H1_error(i);
        
        % print
        line = sprintf('\t%d\t%.8f\t%.8f\t%.8f\t%.8f\n', i, eL2, eH1, logeL2, logeH1);
        fprintf('%s\n', line);
    end


end