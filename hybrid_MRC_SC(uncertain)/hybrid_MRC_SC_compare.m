clear all; clc;

hw1_4c(2);
 hold on 
hw1_4c(4);
hw1_4c(6);

% picture -----------------------------------------------------------------

title('Compare different number of branch of hybrid MRC');

L=legend('2RX hybrid MRC', '2RX MRC', '4RX hybrid MRC' ...
    , '4RX MRC', '6RX hybrid MRC', '6RX MRC');

set(L,'Fontsize', 12);