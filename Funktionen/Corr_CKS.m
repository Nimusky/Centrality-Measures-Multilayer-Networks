% This function computes the correlation and the kendall coefficients
% for 4x2 given centralities.

function [T] = Corr_CKS(ESC, RSC, TC, KATZ, ESC_2, RSC_2, TC_2, KATZ_2)
    
    correlations = zeros(4, 1);
    kendallCoefficients = zeros(4, 1);

    correlations(1) = corr(ESC, ESC_2);
    correlations(2) = corr(RSC, RSC_2);
    correlations(3) = corr(TC, TC_2);
    correlations(4) = corr(KATZ, KATZ_2);

    kendallCoefficients(1) = corr(ESC, ESC_2, 'type', 'Kendall');
    kendallCoefficients(2) = corr(RSC, RSC_2, 'type', 'Kendall');
    kendallCoefficients(3) = corr(TC, TC_2, 'type', 'Kendall');
    kendallCoefficients(4) = corr(KATZ, KATZ_2, 'type', 'Kendall');

   
    T = table(correlations,kendallCoefficients);
    T.Properties.VariableNames = ["Correlation","Kendall Coefficient"];
end

