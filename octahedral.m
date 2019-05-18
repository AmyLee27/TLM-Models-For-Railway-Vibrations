function [oct] = octahedral(epsilon,alpha)

oct = alpha*sqrt((epsilon(1)-epsilon(2)).^2+(epsilon(1)-epsilon(3)).^2+(epsilon(2)-epsilon(3)).^2+6*(epsilon(4).^2+epsilon(5).^2+epsilon(6).^2))./3;
            