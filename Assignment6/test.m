close all
clear all
% figure;
importfile('uComp000.dat');
importfile('uExact000.dat');
% plot(uComp000(:,2));
% hold on
% plot(uExact000(:,2), 'red');
% hold off
% figure;
importfile('uComp001.dat');
importfile('uExact001.dat');
plot(uComp001(:,2));
% figure;
hold on
plot(uExact001(:,2), 'red');
