# Matlab-Code--ODEs
## This is a code written in Matlab, It gives the solution of a Second Order Homogenous Differential Equation with Constant Coefficients and Initial Conditions.

%%
% Prompting user to input values of a, b and c
a = input ('Input value of coefficient "a":  ');
b = input ('Input value of coefficient "b":  ');
c = input ('Input value of coefficient "c":  ');

% Promting user to input initial conditions
m = input ('Input value of "y(0)":');
n = input('Input value of "y_prime(0)":');


% Quadratic Equation
x1 = (-b + sqrt(b.^2-4.*a.*c))/(2.*a)
x2 = (-b - sqrt(b.^2-4.*a.*c))/(2.*a)

% Checking the nature of the roots

    % Case 1: Real and Distinct Roots y(t)= A*e^(x1*t)+ B*e^(x2*t)
    if  b^2 > 4*a*c
       
        % Solving simultaneous equations
        % y(0) = Ae^(x1*0)+ Be^(x2*0) = m and
        % y'(0) = (x1*Ae^(x1*0))+ (x2*Be^(x2*0)) = n
        
        syms A B
        eqn1 = A + B == m;
        eqn2 = (x1*A) + (x2*B) == n;
        
        sol = solve([eqn1,eqn2], [A, B]);
        Asol = sol.A
        Bsol = sol.B
        
        yt = sprintf ('%fe^%ft + %fe^%ft', Asol, x1, Bsol, x2)
           
    % Case 2: Repeated Distinct Roots y(t)= Ae^(x1t)+ Bte^(x2t)
    elseif b^2 == 4*a*c 
         % Solving simultaneous equations
         % y(0) = Ae^(x1*0)+ Bte^(x2*0)= m
         % y'(0) = (x1*A*e^x1*0)+ (B*e^x1*0) + (x1*B*0*e^x1*0) = n
        
        syms A B
        eqn1 = A == m;
        eqn2 = (x1*A) + (B) == n;
        
        sol = solve([eqn1,eqn2], [A, B]);
        Asol = sol.A
        Bsol = sol.B
        
        yt = sprintf ('e^%ft (%f + %ft)', x1, Asol, Bsol)
        
    % Case 3: Complex Conjugate Roots y(t)= e^(α*t)[x1cos(β*t)+ x2sin(β*t)]
    else % (b^2) < 4*a*c
        
         % Solving simultaneous equations
         % y(0) = (A*e^(P*0)*cos(Q*0) + (B*e^(P*0)*sin(Q*0) = m and
         % y'(0) = e^(P*0)(-A*Q*sin(Q*0) + B*Q*cos(Q*0)) + P*e^(P*0)(A*cos(Q*0)+ B*sin(Q*0))= n
         
      % Separating real and imaginary parts of a complex number
         P = real(x1);
         Q = imag (x1);
         
        syms A B
        eqn1 = A == m;
        eqn2 = (B*Q) + (A*P) == n;
        
        sol = solve([eqn1,eqn2], [A, B]);
        Asol = sol.A
        Bsol = sol.B
    
        yt = sprintf ('e^%ft (%f cos(%ft) + %f sin (%ft))' , P, Asol, Q, Bsol, Q)
        
    end
    
   %%
   % Solving using dsolve
   
    syms y(t)
    dy = diff (y);
    eqtn = a*diff (y , t, 2) + b*diff (y,t) + c*y == 0;
 
    % Initial conditions
    conditions = [y(0)== m, dy(0)== n];
    ySol(t) = dsolve (eqtn, conditions)
