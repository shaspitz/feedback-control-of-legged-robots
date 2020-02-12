function [value,isterminal,direction] = two_link_event(t,x)

value = x(2) - 2*x(1); % detect when phi - 2*theta == 0 (approx)
isterminal = 1 ; % stop integration when value == 0
direction = 1 ; % detect zero when function is increasing

end