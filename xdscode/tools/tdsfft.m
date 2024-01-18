function [w, y] = tdsfft(t,X,varargin)
% TDSFFT: Does the FFT of matrix X along the time axis (cols) and returns
% the frequencies and complex FT.
% the frequencies are in THz if t is in ps.
% t is assumed to be a 1D column vector.
% X is a matrix where time runs along columns
% the length of the FFT can be specified, in which case it is padded with
% zeros.
if ~isempty(varargin)
    N = varargin{1};
else
    N = size(X,1);
end
dt = (t(2) - t(1));         % ps; then freq is THz
maxFreq = 1/dt;
%N = 2^(nextpow2(size(X,2))+2);	% finds the next power of 2 for FFT
w = linspace(0, maxFreq, N)';
y = (fft(X,N,1)/N);

w = w(1:round(N/2));
y = y(1:round(N/2),:);
