function p = psychofun(x,theta)
%PSYCHOFUN Evaluate psychometric function with parameter vector THETA at X.

% Psychometric function parameters
psycho_mu = theta(1);
psycho_sigma = theta(2);
psycho_gammalo = theta(3);
psycho_gammahi = theta(4);

% Psychometric curve for choosing RIGHT
p = 0.5*(1 + erf((x-psycho_mu)/(sqrt(2)*psycho_sigma)));

% Add left (low) and right (high) lapses
p = psycho_gammalo + max(0,(1 - psycho_gammalo - psycho_gammahi))*p;

end
