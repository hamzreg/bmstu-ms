function lab_01
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

    X = sort(X);
    n = length(X);

    M_max = max(X);
    M_min = min(X);
    fprintf("\nМаксимальное значение = %.4f\n", M_max);
    fprintf("Минимальное значение = %.4f\n", M_min);

    R = M_max - M_min;
    fprintf("\nРазмах выборки = %.4f\n", R);

    mu = sum(X) / n;
    S2 = sum((X-mu).^2) / (n - 1);
    fprintf("\nОценка математического ожидания = %.4f\n", mu);
    fprintf("Оценка дисперсии = %.4f\n", S2);

    m = floor(log2(n)) + 2;
    fprintf("\nЧисло интервалов = %d\n\n", m);

    d = (X(n) - X(1)) / m;

    J = [];

    for i = 1 : m
      J(i, 1) = X(1) + (i - 1) * d;
      J(i, 2) = X(1) + i * d;
    end
 
    N = zeros(m);
    
    for i = 1 : n
      for j = 1 : m
        if X(i) >= J(j, 1) && X(i) < J(j, 2)
          N(j)++;
        end
      end
    end
    
    N(m)++;
    
    for i = 1 : m - 1
      fprintf("[%.4f,%.4f) - %d\n", J(i, 1), J(i, 2), N(i));
    end
    
    fprintf("[%.4f,%.4f] - %d\n", J(m, 1), J(m, 2), N(m));

    f = [];
    
    for i = 1 : m
      f(i, 1) = (J(i, 1) + J(i, 2)) / 2;
      f(i, 2) = N(i) / (n * d);
    end
    
    args = X(1):1e-4:X(n);
    sigma = sqrt(S2);
    f_normal = normpdf(args, mu, sigma);

    bar(f(:,1), f(:,2), 1, "b");
    hold on;
    plot(args, f_normal, 'g', 'LineWidth', 2);
    grid;

    F = [];
    
    F(1, 2) = 0;
    F(1, 1) = X(1) - 1;
    
    for i = 2 : (n + 1)
      cnt = 0;
      for j = 1 : n
        if X(j) <= X(i - 1)
          cnt++;
        end
      end
      
      F(i, 1) = X(i - 1);
      F(i, 2) = cnt / n;
    end

    F(n + 2, 2) = 1;
    F(n + 2, 1) = X(n) + 1;

    F_normal = normcdf(args, mu, sigma);

    figure;
    stairs(F(:,1), F(:,2), "b");
    hold on;
    plot(args, F_normal, "r", 'LineWidth', 2);
    grid;
endfunction
