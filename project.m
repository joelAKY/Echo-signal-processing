
clc; clear all ; close all

%function norm_corrfunctin

%sys = tf([1],[1 1 10]);   % transfer function
%sys = tf([1], [0.5  0.5 2]);   % transfer function

FilterOrder = 11;
h1 =  fir1(FilterOrder-1, 0.9)';  % imulse filter response return FilterOrder  coefficients

plot(h1) ; grid ;

[Y filter_delay] = max(abs(h1)) ;

D = 100;

h = [zeros(D,1); h1] ;    % adding delay to the filter
plot(h) ; grid ;

N = 5000; % signal length

noise = randn(N,1);
xref = wgn(N,1,0);  %reference signal,  white gaussion noise  
xref = xref./max(abs(xref)) ; %normalize amplitude range -1 to 1 
xinp = filter (h, 1, xref) ; % echo signal
noise = noise * std(xinp)/(10*std(noise));

xinp = xinp + noise; % adding some noise to input

figure(1) ;
L = min(length(xref), length(xinp)) ;
subplot(211); plot(1:L, xref(1:L), 'b') ; grid ;
xlabel('Time (samples)'); ylabel('Amplitude')
legend ('xref') ;
subplot(212); plot(1:L, xinp(1:L), 'b') ; grid ;

xlabel('Time (samples)'); ylabel('Amplitude')
legend ( 'xinp') ;

 %finding delay on xinp signal
 %{
[acor,lag] = xcorr(xref,xinp);
[cor,I] = max(abs(acor));
lagDiff = -lag(I);
%}
 
alpha = 0.1; 
one_minus_alpha = 1-alpha ;

% finding normalized cross correlation
N=length(xinp); %// length of signal
M=512;  %// length of window
max_lag_supported = 200 ;
%M_2 = M/2 ; %overalp window size 50%
%{%
N = N - M ;
%xf=framing(xinp,M,M_2); %this function frames the signal i will get xf(128,14)
%max_lag_supported= floor(N/M);
M0 = 1;  
%ccor=NaN(1,win_num);
Px= 0 ;
Py = 0 ; Pxy = 0 ;
tCxy_max = 0 ;
idx = 1 ;
for m=1:max_lag_supported
    x1 = xref(1:M) ;
    Px = Px *alpha + one_minus_alpha * (x1'*x1) ;
     y1 = xinp( (m-1) + (1:M)) ;
     Py =  Py * alpha + one_minus_alpha * (y1' * y1) ;
     Pxy = Pxy * alpha + one_minus_alpha * (x1'*y1) ;
     tPxy = abs(Pxy/sqrt(Px + Py)) ;
     Cxy (idx) = tPxy ;
     if tCxy_max < tPxy
         tCxy_max = tPxy ;
         index= m ;
     end
     idx = idx + 1 ;
 end
figure(2) ;
plot(Cxy) ; grid ;
title('Correlation coefficients') ;
%}%
lagDiff = index  ; 
lagDiff = lagDiff - filter_delay ; % compansate the fir filter delay
L = length(xref) ;
xinp_delayed = xinp (lagDiff+1: L) ; % delay compansated input signal
L = min(length(xref), length(xinp_delayed)) ;

[xout w1] =  lms_filter (xref, xinp_delayed, FilterOrder) ;  % lms filter

L = min(min(length(xref), length(xinp)), length(xout)) ;

figure (3);
subplot (211); plot(1:FilterOrder, h1(1:FilterOrder), 'b') ; grid
legend('Actual filter weights') ;
title('Comparison of the actual/estimated filter weights') ;
subplot(212); plot( 1:FilterOrder, w1(1:FilterOrder), 'r') ; grid
legend('Estimated filter weights') ;


figure(4) ;
subplot(311); plot(1:L, xref(1:L), 'b') ; grid ;
xlabel('Time (samples)'); ylabel('Amplitude')
title('ref/inp/out signals') ;
subplot(312); plot(1:L, xinp(1:L), 'b') ; grid ;
xlabel('Time (samples)'); ylabel('Amplitude')
subplot(313); plot(1:L, xout(1:L),'b') ; grid ;
axis([0  L -1 1]) ;
xlabel('Time (samples)');  ylabel('Amplitude')



return ;

% lms filter algorithma with filter order =N=11
%xref = reference signal
%xinp = echo signal
%N = filter order
%e = output (error) signal ;
%w= estimated (adaptive) filter weights 
%e = y - w' * u;
%w = w + mu * u./u^2 * e


function [e w] =  lms_filter (xref, xinp, N) ;

%begin of algorithm
L = min(length(xref), length(xinp)) ;
 mu=0.25; % step size parameter
w = zeros ( N  , 1 ) ;
for n = N : L 
    u = xref(n:-1:n-N+1) ;
    y(n)= w' * u;
    e(n) = xinp(n) - y(n) ;
	w = w + mu * u./(u'*u) * e(n) ;
end 
%plot(w) ; grid ;
end

