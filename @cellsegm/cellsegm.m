%    CELLSEGM Class definition file
%
%     =======================================================================================
%     Copyright (C) 2013  Erlend Hodneland
%     Email: erlend.hodneland@biomed.uib.no 
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%     =======================================================================================
%
classdef cellsegm
       
   methods (Static)
      % main functions
      prop = cellprop(im,wat,cellv,propname,h,meanintbck);
      cellsegmentation(varargin);
      [minvol,minvolvox,maxvol,maxvolvox] = cellsize(minvol,maxvol,h,just,O);
      [cellbw,infocells] = classifycells(varargin);
      [minima,minimacell,prm] = getminima(varargin);
      readbioformat(varargin);
      readbioformatmatlab(varargin);
      [cellbw,wat,imsegm,prmout] = segmct(varargin);
      [cellbw,wat,imsegm,prmout] = segmneuron(varargin);      
      [cellbw,wat,imsegm,minima,minimacell,info] = segmsurf(varargin);
      [im] = smoothim(varargin);
      [cellbw] = splitcells(cellbw,splitth,splitvolvox,h);
      viewsegm(varargin);
      [faser] = mergefragments(faser,im,thint,thconv,optlog,optint);
      [faser] = resegmentwat(faser,im);
      [v] = ridgeenhhessian(varargin);
      show(im,fignum);
      showall(varargin);
      show4D(im);
      printcell(fid,A,thmin);
      A = mat2celldirect(A);
      panelstruct(varargin);
      [u] = imresize3d(im,g,method);
      [u] = cohenhdiff(varargin);
      [u] = edgeenhdiff(varargin);
      [u] = dircohenh(varargin);
      
   end 
      
end