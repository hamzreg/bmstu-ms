function lab_02
    X = [14.90,14.40,13.56,15.55,13.97,16.33,14.37,13.46,15.51,14.69,...
         13.41,14.24,15.65,14.54,13.55,13.15,14.32,15.04,13.27,14.60,...
         13.83,13.93,14.11,14.15,15.48,15.96,14.46,13.87,13.67,15.30,...
         13.95,16.08,18.25,14.93,15.37,14.38,15.56,13.92,14.23,12.80,...
         13.16,13.89,14.24,13.90,12.82,13.20,13.89,13.50,13.44,16.13,...
         14.68,15.27,13.35,13.62,16.16,16.46,13.83,14.13,15.68,15.22,...
         12.59,12.94,13.09,16.54,14.61,14.63,14.17,13.34,16.74,16.30,...
         13.74,15.02,14.96,15.87,16.03,12.87,14.32,14.48,14.57,14.43,...
         12.61,14.52,15.29,12.07,14.58,11.74,14.97,14.31,12.94,12.82,...
         14.13,14.48,12.25,14.39,15.08,12.87,14.25,15.12,15.35,12.27,...
         14.43,13.85,13.16,16.77,14.47,14.89,14.95,14.55,12.80,15.26,...
         13.32,14.92,13.44,13.48,12.81,15.01,13.19,14.68,14.44,14.89];
    n = length(X);

    mu = sum(X) / n;
    S2 = sum((X-mu).^2) / (n - 1);
    fprintf("Выборочное среднее = %.4f\n", mu);
    fprintf("Исправленная выборочная дисперсия = %.4f\n", S2);
    
    gamma = 0.9;
    mu_bottom = get_mu_bottom(gamma, n, mu, S2);
    mu_top = get_mu_top(gamma, n, mu, S2);
    fprintf("\n%.1f-доверительный интервал для математического ожидания: (%.4f, %.4f)\n", gamma, mu_bottom, mu_top);
    
    S2_bottom = get_S2_bottom(gamma, n, S2);
    S2_top = get_S2_top(gamma, n, S2);
    fprintf("\n%.1f-доверительный интервал для дисперсии: (%.4f, %.4f)\n", gamma, S2_bottom, S2_top);

    mu_arr = zeros(1, n);
    S2_arr = zeros(1, n);
    mu_bottom_arr = zeros(1, n);
    mu_top_arr = zeros(1, n);
    S2_bottom_arr = zeros(1, n);
    S2_top_arr = zeros(1, n);
    Rm_arr = zeros(1, n);
    Rs_arr = zeros(1, n);

    for i = 1 : n
      mu_arr(i) = sum(X(1:i)) / i;
      
      if i == 1
        S2_arr(i) = 0;
      else
        S2_arr(i) = sum((X(1:i)-mu_arr(i)).^2) / (i - 1);
      end
      
      Rm_arr(i) = get_Rm(gamma, i, S2_arr(i));
      Rs_arr(i) = get_Rs(gamma, i, S2_arr(i));
      
      mu_bottom_arr(i) = get_mu_bottom(gamma, i, mu_arr(i), S2_arr(i));
      mu_top_arr(i) = get_mu_top(gamma, i, mu_arr(i), S2_arr(i));
      
      S2_bottom_arr(i) = get_S2_bottom(gamma, i, S2_arr(i));
      S2_top_arr(i) = get_S2_top(gamma, i, S2_arr(i));
    end
    
    clf;
    start = 10;
    plot(start:n, zeros(n - start + 1) + mu, "k");
    hold on;
    plot(start:n, mu_arr(start:n), "r");
    plot(start:n, mu_bottom_arr(start:n), "g");
    plot(start:n, mu_top_arr(start:n), "b");
    xlabel("n");
    ylabel("y");
    %xlim([10, 120]);
    set(gca, 'FontSize', 16);
    l = legend('y = \mu^{\^}(x_N)', 'y = \mu^{\^}(x_n)', 'y = \mu_{bottom}(x_n)', 'y = \mu_{top}(x_n)');
    set(l, "interpreter", "tex");
    grid on;
    figure;
    plot(start:n, zeros(n - start + 1) + S2, "k");
    hold on;
    plot(start:n, S2_arr(start:n), "r");
    plot(start:n, S2_bottom_arr(start:n), "g");
    plot(start:n, S2_top_arr(start:n), "b");
    xlabel("n");
    ylabel("z");
    %xlim([10, 120]);
    set(gca, 'FontSize', 16);
    l = legend('z = S^2(x_N)', 'z = S^2(x_n)', 'z = \sigma^2_{bottom}(x_n)', 'z = \sigma^2_{top}(x_n)');
    set(l, "interpreter", "tex");
    grid on;
    figure;
    plot(start:n, Rm_arr(start:n), "k");
    hold on;
    plot(start:n, Rs_arr(start:n), "r");
    xlabel("n");
    ylabel("R");
    %xlim([10, 120]);
    set(gca, 'FontSize', 16);
    l = legend('R = R_m(n)', 'R = R_\sigma^2(n)');
    set(l, "interpreter", "tex");
    grid on;
endfunction

function mu_bottom = get_mu_bottom(gamma, n, mu, S2)
    t = tinv((1 + gamma)/2, n - 1);
    mu_bottom = mu - sqrt(S2) * t / sqrt(n);
end

function mu_top = get_mu_top(gamma, n, mu, S2)
    t = tinv((1 + gamma)/2, n - 1);
    mu_top = mu + sqrt(S2) * t / sqrt(n);
end

function S2_bottom = get_S2_bottom(gamma, n, S2)
    h_bottom = chi2inv((1 + gamma)/2, n - 1);
    S2_bottom = (n - 1) * S2 / h_bottom;
end

function S2_top = get_S2_top(gamma, n, S2)
    h_top = chi2inv((1 - gamma)/2, n - 1);
    S2_top = (n - 1) * S2 / h_top;
end

function Rm = get_Rm(gamma, n, S2)
    t = tinv((1 + gamma)/2, n - 1);
    Rm = 2 * sqrt(S2) * t / sqrt(n);
endfunction

function Rs = get_Rs(gamma, n, S2)
    h_bottom = chi2inv((1 + gamma)/2, n - 1);
    h_top = chi2inv((1 - gamma)/2, n - 1);
    Rs = (n - 1) * S2 * (h_bottom - h_top) / (h_bottom * h_top);
endfunction
