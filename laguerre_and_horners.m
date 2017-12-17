% the main function here is laguerre, below will be several helper
% functions I needed to complete evaluations or find polished roots. I'm
% merging everything into one .m file but in order to run, please seperate
% each helper function into their own .m file and run test.m
function[root] = laguerre(func)
    errorTol = 10^-6; % get error tolerance
    n = feval(symengine,'degree',func); % find leading degree
    roots = []; % create a vector of all the roots
    Q = func; % Q(x)
    space = linspace(-10,10,161); % evenly spaced points
    for i = 1:n % run everything for n times
        handle = matlabFunction(Q);
        arr = arrayfun(handle,space); % create an array of all function evaluations
        root = min(arr); % find min
        count = 1;
        while (abs(horners(Q, root)) > errorTol && count < 10) % using count to have it not run forever
            G = horners(diff(Q,1),root) / horners(Q,root); % evaluate G and H
            H = G^2 - (horners(diff(Q,2),root))/horners(Q,root); 
            test_0 = G+sqrt((n-1)*(n*H-G^2)); % figure out which makes the abs value of denom greatest
            test_1 = G-sqrt((n-1)*(n*H-G^2));
            if (abs(test_0) >= abs(test_1))
               a = n / test_0;
            else
                a = n / test_1;
            end
            root = root - a; % adjust root
            count = count + 1;
        end
        if (count ~= 10) % if it didn't flow over then we have our desired root!
            roots = [roots root]; % add this found root into our vector of roots!      
            num = sym2poly(Q); % transform into poly
            den = [];
            den(1) = 1;
            den(2) = -1*root;
            [quot, rem] = deconv(num, den); % a way of finding Q using division
            Q = poly2sym(quot); % discard remainder
        else
            break;
        end
    end
    disp("Original roots: ");
    disp(roots);
    disp("Polished roots: ");
    disp(polished(func, roots));
end

%%% POLISHED FUNCTION %%%
function[pol_roots] = polished(P, roots)
    pol_roots = []; % create vector of polished roots
    for i = 1:length(roots)
        root = roots(i);
        % As mentioned in the email, here is where things can go wrong
        % since it will be very difficult to achieve evaluation <= 10^-24
        while (abs(horners(P,root)) > 10^-24) % repeat same process
            G = horners(diff(P,1),root) / horners(P,root); 
            H = G^2 - (horners(diff(P,2),root))/horners(P,root); 
            test_0 = G+sqrt((n-1)*(n*H-G^2));
            test_1 = G-sqrt((n-1)*(n*H-G^2));
            if (abs(test_0) >= abs(test_1))
               a = n / test_0;
            else
                a = n / test_1;
            end
            root = root - a;
        end
        pol_roots = [pol_roots root];
    end
end

%%% HORNERS FUNCTION %%%
function[sum] = horners(func, p)
    coeff_arr = sym2poly(func); % create array of coefficients
    sum = coeff_arr(1);
    i = 1;
    n = feval(symengine, 'degree', func); % find leading degree
    while (n > 0)
            sum = sum * p + coeff_arr(i+1); % b_k = b_k+1 * p + a_k
            n = n - 1;
            i = i + 1;
    end
end

%%% TEST FILE %%%
syms x;
fun = x^2 - 2*x - 3;
laguerre(fun);

    
    
    
    
    
