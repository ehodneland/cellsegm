function prm = prmfilenucleus
%   =======================================================================================
%   Copyright (C) 2013  Erlend Hodneland
%   Email: erlend.hodneland@biomed.uib.no 
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%   =======================================================================================

% ridge filtering
prm.segmsurf.filterridges = 0;

% smoothing method
prm.segmsurf.smoothim.method = 'dirced';

% minima method
prm.segmsurf.getminima.method = 'nucleus';

% threshold in minima
prm.segmsurf.getminima.nucleus.segmct.thrs.th = 0.8;

% split cells in minima
prm.segmsurf.getminima.nucleus.segmct.split = 1;

% threshold for splitting cells in minima
prm.segmsurf.getminima.nucleus.segmct.splitth = 1;

% classification method
prm.segmsurf.classifycells.method = 'minimacell';

% segmentation plane start
prm.segmplane = 3;

% segmentation channel
prm.segmch = 2;

% nucleus channel
prm.nucleusch = 1;

% segmentation method
prm.method = 'segmsurf';

% subtract nucleus channel from surface stain
prm.subtract = 1;