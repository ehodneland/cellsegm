% MAKESTR Making a string
% STR = MAKESTR(STR,TH)  of STR that is TH long.
% Useful when writing to file in columns.
% NB: Must first use NUMSTR to create the string
%
function [str] = makestr(str,th)

str = [str blanks(th-max(size(str)))];