function Lfy = Lfy_gen(in1)
%LFY_GEN
%    LFY = LFY_GEN(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    20-Oct-2019 13:56:53

dq1 = in1(8,:);
dq2 = in1(9,:);
dq3 = in1(10,:);
Lfy = [dq3;dq1+dq2+dq3.*2.0];