function [sparsityvec,changevec,decomplevel] = decomplevelselect(data,wname)
% DECOMPLEVELSELECT Optimal sparsity-based decomposition level selection
% for noise thresholding in wavelet denoising.
% 
%   [SPARSITYVEC,CHANGEVEC,DECOMPLEVEL] = DECOMPLEVELSELECT(DATA,'wname')
%   returns the vector of sparsity component values for each decomposition
%   level, the vector of change in sparsity component values for each
%   decomposition level, and the optimal sparsity-based decomposition level,
%   respectively, of 1-D signal data using the wavelet named in the string
%   or character array 'wname' (see WFILTERS for more information).
%   DATA must be a 1-D array with finite, numeric elements.
%
%   W. Bekerman and M. Srivastava, 01-Dec-2020.
%   Last Revision: 16-Jul-2021.
%   See also WAVEDEC, DETCOEF.

%   Copyright 2021, William J. Bekerman
%   
%   Licensed under the Apache License, Version 2.0 (the "License");
%   you may not use this file except in compliance with the License.
%   You may obtain a copy of the License at
%   
%       http://www.apache.org/licenses/LICENSE-2.0
%   
%   Unless required by applicable law or agreed to in writing, software
%   distributed under the License is distributed on an "AS IS" BASIS,
%   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%   See the License for the specific language governing permissions and
%   limitations under the License.

    narginchk(2,2);
    validateattributes(data,{'numeric'},{'vector','finite','real'},'decomplevelselect','DATA');
    validateattributes(wname,{'string','char'},{},'decomplevelselect','WNAME');
    
    f = data;
    N = floor(log2(length(f))); %% MAXIMUM POSSIBLE DECOMPOSITION LEVEL
    
    [C,L] = wavedec(f,N,wname);
    c = detcoef(C,L,1:N);
    filter_len = length(wfilters(wname));
    
    %{
    MAXIMUM DECOMPOSITION LEVEL BASED ON RATIO BETWEEN LENGTH OF DETAIL
    COMPONENET AND FILTER_LEN
    %}
    N_ratio = 0;
    ratio = Inf;
    while ratio > 1.5 && N_ratio <= N
        % ONLY CONSIDER DECOMPOSITION LEVELS WITH RATIO > 1.5
        N_ratio = N_ratio + 1;
        ratio = length(c{N_ratio})/filter_len;
    end
    N = N_ratio-1;
    
    sparsityvec = zeros(1,N);
    for level = 1:N    
        cd = c{level};
        sparsityvec(level) = max(abs(cd))/sum(abs(cd));     
    end

    changevec = zeros(1,N);
    for level = 2:N    
        changevec(level) = sparsityvec(level)-sparsityvec(level-1);
    end
    
    decomplevel = -1;
    for level = 2:N
        if changevec(level) > 0.05 %% FIVE PERCENT CHANGE IN SPARSITY THRESHOLD
            decomplevel = level-1;
            break;
        end
    end
    
    if decomplevel == -1
        error('Error: No change in sparsity of at least five percent detected.')
    end
    
end


