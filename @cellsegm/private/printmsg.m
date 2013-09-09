% PRINTMSG Print message to screen
% 
% PRINTMSG(MSG,FLAG) prints the message in MSG with uppercase if FLAG == 1,
% else lowercase
%
function [] = printmsg(varargin)

msg = varargin{1};

flag = 0;
if nargin == 2
    flag = varargin{2};
end;
if flag == 1
    msg = upper(msg);
end;

level = 1;
if nargin == 3
    level = varargin{3};
end;
    
for i = 1 : level-1
    msg = ['    ' msg];
end;
disp(sprintf('%s',msg))
