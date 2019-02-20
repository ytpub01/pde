close all
clear all
importFile('uComp001.dat');
importFile('uExact001.dat');
plot(uComp001(:,2));
hold on
plot(uExact001(:,2), 'red')