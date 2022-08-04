function [f0,BW,n,nz]=Parasinput()
f0=input('f0(center frequency/MHz)=')*1e6;
BW=input('BW(bandwidth/MHz)=')*1e6;
n=input('n(order)=');
nz=input('nz(number of transmission zeros)=');

