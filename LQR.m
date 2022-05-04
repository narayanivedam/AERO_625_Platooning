function [k,Qd,Rd,Nd,s,e] = lqrdjv(a,b,q,r,nn,Ts)
%%%%%%function [k,s,e,Qd,Rd,Nd] = lqrdjv(a,b,q,r,nn,Ts)
%LQRDJV  Discrete linear quadratic regulator design from continuous 
%       cost function.
%
%   ********************************************************************
%
%           =====>  modified by J. Valasek, 5 Dec 94   <=====
% 
%           This routine now sends back the Q^, R^, and M^ matrices.
%
%           NOTE: the order of the passed-back arguments is NOT the
%                 same as the original MATLAB version.
%
%   ********************************************************************
%
%	[K,S,E] = LQRD(A,B,Q,R,Ts) calculates the optimal feedback gain 
%	matrix K such that the discrete feedback law  u[n] = -K x[n] 
%	minimizes a discrete cost function equivalent to the continuous 
%	cost function
%		J = Integral {x'Qx + u'Ru} dt
%                                                       .
%	subject to the continuous constraint equation:  x = Ax + Bu
%
%	Also returned is S, the discrete Riccati equation solution, and 
%	the closed loop eigenvalues E = EIG(Ad-Bd*K).
%
%	The gain matrix is determined by discretizing the continuous plant
%	(A,B,C,D) and continuous weighting matrices (Q,R) using the sample
%	time Ts and the zero order hold approximation. The gain matrix is
%	then calculated using DLQR.
%
%	[K,S,E] = LQRD(A,B,Q,R,N,Ts) includes the cross-term N that 
%	relates u to x in the cost function.
%		J = Integral {x'Qx + u'Ru + 2*x'Nu}
%
%	See also: C2D, LQED, DLQR, and LQR.

%	Clay M. Thompson 7-16-90
%	Copyright (c) 1986-93 by the MathWorks, Inc.

% Reference: This routine is based on the routine JDEQUIV.M by Franklin, 
% Powell and Workman and is described on pp. 439-441 of "Digital Control
% of Dynamic Systems".

error(nargchk(5,6,nargin));
error(abcdchk(a,b));
[nx,na] = size(a); 
[nb,nu] = size(b);

[nq,mq] = size(q);
if (nx ~= nq) | (nx ~= mq), error('A and Q must be the same size.'); end
[nr,mr] = size(r);
if (mr ~= nr) | (nu ~= mr), error('B and R must be consistent.'); end

if nargin==5,
  Ts = nn;
  nn = zeros(nb,nu);
else
  [nnn,mn] = size(nn);
  if (nnn ~= nx) | (mn ~= nu), error('N must be consistent with Q and R.'); end
end

% Check if q is positive semi-definite and symmetric
if any(eig(q) < -eps) | (norm(q'-q,1)/norm(q,1) > eps)
  disp('Warning: Q is not symmetric and positive semi-definite');
end
% Check if r is positive definite and symmetric
if any(eig(r) <= -eps) | (norm(r'-r,1)/norm(r,1) > eps)
  disp('Warning: R is not symmetric and positive definite');
end

% Discretize the state-space system.
[ad,bd] = c2d(a,b,Ts);

% --- Determine discrete equivalent of continuous cost function ---
n = nx+nu;
Za = zeros(nx); Zb = zeros(nx,nu); Zu = zeros(nu);
M = [ -a' Zb   q  nn
      -b' Zu  nn'  r
      Za  Zb   a   b
      Zb' Zu  Zb' Zu];
phi = expm(M*Ts);
phi12 = phi(1:n,n+1:2*n);
phi22 = phi(n+1:2*n,n+1:2*n);
QQ = phi22'*phi12;
QQ = (QQ+QQ')/2;		% Make sure QQ is symmetric
Qd = QQ(1:nx,1:nx) ;
Rd = QQ(nx+1:n,nx+1:n) ;
Nd = QQ(1:nx,nx+1:n) ;

% Design the gain matrix using the discrete plant and discrete cost function
[k,s,e] = dlqr(ad,bd,Qd,Rd,Nd);

