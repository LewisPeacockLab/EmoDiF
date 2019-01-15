function savesamesize(h, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: This a third party file written by Richard Quist
% Url: http://www.mathworks.com/matlabcentral/fileexchange/17868-savesamesize
% Copyright (c) 2009, The MathWorks, Inc.
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%     * Neither the name of the The MathWorks, Inc. nor the names 
%       of its contributors may be used to endorse or promote products derived 
%       from this software without specific prior written permission.
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SAVESAMESIZE print/export figure h at same size as on screen.
%
% from http://www.mathworks.com/matlabcentral/fileexchange/17868
%
% SAVESAMESIZE will adjust printing/exporting related properties so the
% figure (h) will print/export the same size as the onscreen display. It 
% will print/export the figure and then restore the changed properties to
% the original values
%
% saveSameSize accepts some optional parameter-value pairs. If used you must
% specify BOTH the parameter name (e.g. 'format') AND the parameter value 
% (e.g. 'png'):
%     parameter      value/meaning
%     ----------     --------------------------------------
%      'format'      output format to generate. 
%                    Must be one of the formats supported by the print command
% 
%                    If not specified, jpeg is used as the default
%
%      'file'        path/name of output file to create
%
%                    If not specified, 'figure<N>.<ext>' will be used 
%                    where  <N> is the figure handle and ext 
%                           <ext> is the extension corresponding to the 
%                                 output format chosen
%
%      'renderer'    Renderer to use when creating the output. 
%                    Must be one renderers supported by MATLAB
%
%                    If not specified, the figure's current renderer will be used 
%
% Example usage:
%     save current figure to jpeg file named figure<n>.jpg using the 
%     figure's current renderer
%         saveSameSize(gcf); 
%
%     save current figure to png file named test.png using the 
%     OpenGL renderer
%         saveSameSize(gcf, 'format', 'png', 'renderer', 'opengl',... 
%                        'file', 'test.png'); 
%
% See also PRINT 

%   Copyright 2007 The MathWorks, Inc.

if ishandle(h) 
    setprops = {'PaperType', 'PaperSize', 'PaperUnits', 'PaperPosition', 'PaperPositionMode'};
    try 
        settings = get(h, setprops);
    catch 
        error('saveSameSize unable to get figure properties');
    end
    
    oldSettings = cell2struct(settings, setprops, 2); % to use when resetting values
    set(h, 'PaperPositionMode', 'auto');
    
    paramSettings = processArgs(h, varargin{:});
    if isstruct(paramSettings) 
        try 
            print(h, '-r0', ...
                 paramSettings.renderer, ...
                 paramSettings.device, ...
                 paramSettings.outfile);
        catch
            % print had a problem - restore settings and rethrow
            set(h, oldSettings);
            rethrow(lasterror);
        end
        set(h, oldSettings);
    end
else
    error('saveSameSize requires a handle to a figure');
end

end

% helper function to deal with parsing command line arguments
%   expecting   ..., 'format', ['jpeg', 'eps', 'png', ...]
%                    'file', 'pathToOutputFileToCreate'
%                    'renderer', 'opengl|painters|zbuffer'
function param_settings = processArgs(h, varargin) 

    if mod(length(varargin),2) ~= 0
       error('saveSameSize:InvalidOptionalParameters', ...
             'saveSameSize optional parameters must be specified as parameter-value pairs');
    end
    %setup defaults
    param_settings.outfile = ''; 
    ext = 'jpg';
    param_settings.device = '-djpeg'; 
    param_settings.renderer = lower(['-' get(h, 'Renderer')]);

    % GJD 090428: for some unfathomable reason, the renderer keeps
    % getting set to '-none', even though it's 'painters' if you run
    % things interactively...???
    if strcmp(param_settings.renderer,'-none')
      param_settings.renderer = '-painters';
    end

    for i = 1 : 2 : length(varargin)-1 
        param_arg = varargin{i};
        % caller specified output device to use
        if ischar(param_arg) && strcmpi(param_arg, 'format')
           val_arg = varargin{i+1};
           % remove '-d' if caller specified it
           if length(val_arg) >= 3 && strcmpi(val_arg(1:2),'-d')
               val_arg = val_arg(3:end); 
           end
           % check device against list of supported devices
           [device, ext] = localIsGoodDevice(val_arg);
           param_settings.device = ['-d' device];
        end 

        % caller specified output file to use 
        if ischar(param_arg) && strcmpi(param_arg, 'file')
           val_arg = varargin{i+1};
           if ischar(val_arg) 
    		  param_settings.outfile = val_arg;
           else 
              error('saveSameSize:InvalidOutputFile', ...
                     'saveSameSize invalid output file specified');
           end
        end

        % caller specified renderer to use 
        if ischar(param_arg) && strcmpi(param_arg, 'renderer')
           val_arg = varargin{i+1};
           if ischar(val_arg) 
              if strcmpi(val_arg, 'OpenGL')
                 renderer = 'opengl';
              elseif strcmpi(val_arg, 'painters')
                 renderer = 'painters';
              elseif strcmpi(val_arg, 'zbuffer')
                 renderer = 'zbuffer';
              else 
                 error('saveSameSize:InvalidRenderer', ...
                     'saveSameSize invalid renderer <%s> specified', val_arg);
              end   
              param_settings.renderer = ['-' renderer];
           elseif isnumeric(val_arg)
               error('saveSameSize:InvalidRenderer', ...
                     'saveSameSize invalid renderer <%s> specified', num2str(val_arg));
           else 
               error('saveSameSize:InvalidRenderer', ...
                     'saveSameSize invalid renderer specified');
           end
        end
    end
    if isempty(param_settings.outfile) 
       param_settings.outfile = ['figure' num2str(h) '.' ext];
    end
end

% helper function to check specified device (dev) against the list of 
% supported devices.
% Returns device/format string for the supported device and the extension
% that corresponds to output of that format.
function [device, ext] = localIsGoodDevice(dev) 
  pj = printjob; %#ok
  device = '';
  [options, devices, extensions ] = printtables( pj ); %#ok
  devIndex = strmatch( dev, devices, 'exact' );
  if length(devIndex) == 1 % one exact match
     device = dev;
  elseif length(devIndex) > 1 
     devIndex = strmatch( dev, devices);
     if length(devIndex) == 1 % one partial match
        device = devices{devIndex};
     end
  end
  if isempty(device) 
    error('saveSameSize:InvalidFormat', ...
          'saveSameSize invalid format specified ');
  end
  ext = extensions{devIndex};
  if isempty(ext) 
    error('saveSameSize:NotFileFormat', ...
          'saveSameSize format specified does not support exporting to a file ');
  end
end 

