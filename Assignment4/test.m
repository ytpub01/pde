close all
clear all
importfile('uComp001.dat');
importfile('uExact001.dat');
plot(uComp001(:,2));
hold on
plot(uExact001(:,2), 'red')