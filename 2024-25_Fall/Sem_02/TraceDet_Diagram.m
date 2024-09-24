%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         3. practice         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Nonlinear dynamics systems  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         2019.10.01.         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Trace-Determinant diagram  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

% plot trajectories from given points,
% change the matrix, the phase portrait also changes (eigenvalue, det, trace also shown)
x = -4:1:4;
y = -4:1:4; 
[X Y] = meshgrid(x,y);

K = 500;
fix = zeros(K,length(x)*length(y)); % initial points for ode45
fiy = zeros(K,length(x)*length(y));

figure;
% plot trace-determinant curve
subplot(1,2,2);
tr=(-8:0.1:8);
plot(tr, tr.^2./4, 'b', tr, zeros(size(tr)), 'b', zeros(size(tr)), tr, 'b');
xlabel('trace');
ylabel('det');
hold on;
%
% change parameter 'a'
border=6;
for a = -border:0.2:border
    % the system:
    % ddx-2a*dx+(4-a^2)*x=0
    % ddx=(a^2-4)*x+2a*dx
    % dx=y;
    % dy=(a^2-4)x+2a*y;
    % Jacobi matrix:
    A = [0 1; -4+a^2/4 a] % matrix form of the equations
   
    se=eig(A); % eigenvalues
    % solve the diffequ for each system:
    n=1;
    for i = 1:length(x)
        for j = 1:length(y)
            %
            if trace(A)>0 % in unstable cases solutions run far from the initial points, hard to see.., because of that:          
               A1=-A; % reverse time (-1*matrix): the initial points become the end points of the trajectories
            else
                A1=A;
            end
            %}
            [t,fi] = ode45(@(t,y) A1*y, linspace(0,5,K), [x(i) y(j)]); %calculate ode45 in given points
            fix(:,n) = fi(:,1); % store the results in the given rpws of 'fix' and 'fiy'
            fiy(:,n) = fi(:,2);
            n = n+1;
        end
    end
    
    %plot trajectories from given points:
    subplot(1,2,1);
    plot(fix,fiy);
    hold on;
    plot(X,Y,'x');
    xlim([2*min(x) 2*max(x)]);
    ylim([2*min(y) 2*max(y)]);
    title(['eigenvalues: ', num2str(se(1)), ', ', num2str(se(2))]);
    xlabel('x');
    ylabel('y');
    hold off
    
    %trace det plot:
    subplot(1,2,2);
    plot(trace(A), det(A), 'r*');
    title(['trace(A): ', num2str(trace(A)), ', det(A):', num2str(det(A))]);
    
    waitforbuttonpress;
    %pause(0.1)

end