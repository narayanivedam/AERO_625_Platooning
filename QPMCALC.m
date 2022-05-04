%  qpmcalc.m
%  This M-file assembles the Quad Partition Matrix (QPM) and returns
%  the sub-matrices X12 and X22 from passed A, B, C, and D.
%
%  =================================================================
%
%  written by:   J. Valasek
%                WMU Aircraft Design and Control Laboratory
%                19 January 1996
%                
%  =================================================================
%       DIGITAL FLIGHT CONTROL SYSTEMS: Analysis and Design
%                               by
%                          David R. Downing
%                          John Valasek

%   The University of Kansas - Division of Continuing Education

%                          1 September 2003
%  =================================================================
function [X12, X22] = QPMCALC(A, B, C, D) ;

%
% .. determine dimenstions of passed matrices and vectors
%

[rowa, cola] = size(A) ;

[rowb, colb] = size(B) ;

[rowc, colc] = size(C) ;

[rowd, cold] = size(D) ;

%
% .. quad partition matrix and its inverse
%

qpm  = [A, B ; C, D]   ;

qpmi = inv(qpm)       ;

%
% .. break out the X12 and X22 quadrants
%

X12 = qpmi(1:rowa, cola+1:cola+colb)           ;

X22 = qpmi(rowa+1:rowa+rowc, cola+1:cola+colb) ;
