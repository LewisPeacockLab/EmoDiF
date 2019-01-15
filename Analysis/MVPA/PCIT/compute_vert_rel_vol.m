function [] = compute_vert_rel_vol()
  
  total = 0;
  for k = 1:3
    total = total + compute_vert_rel_vol_per_branch(k); % 6.1676e+16
  end
  
  resolution = 0.0001;
  total_possible = length(-1:resolution:1)^4;
  
  fprintf('\nTOTAL: total_entries = %d (%.2f%% of total)\n\n', ...
    total, 100*total/total_possible);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[super_total_entries] = compute_vert_rel_vol_per_branch(branch)
  
  % resolution = 0.0001;
  resolution = 0.0001;
  num_chunks = 8;
  save_dir = '~/fmri/repref/p-cit/vol_comp';
  
  switch branch
    case 1
      
      %   Branch I: y1 defines the starting point, y2 defines the dip and y3 defines the rise
      %   -1 <= y2 < y1, y2 is the dip so it must fall below the starting point (y1)
      %   y2 < y3 <= 1, y3 is the rise so it must fall above the dip (y2)
      %   -1 <= y4 <= y3, y4 can hold any value that is below the rise (y3)
      %   y2 < y1 <= 1, y1 can hold any value that is above the dip (y2)
      
      y = -1:resolution:1;
      num_y2s = length(y);
      total_entries = nan(num_y2s,1);
      
      for r = 1:num_y2s
        
        if ~mod(r,floor(0.01*num_y2s))
          fprintf('%d%% ', ceil(100*r/num_y2s));
        end
        
        my_y2 = y(r);
        
        my_y1  = y(y > my_y2);
        my_y3 = y(y > my_y2);
        
        num_y1 = length(my_y1);
        num_y3 = length(my_y3);
        if ~num_y1; continue; end
        
        num_y3_and_y4 = sum( (r+1):num_y2s );
        total_entries(r) = num_y1 * num_y3_and_y4;
      end
      
      super_total_entries = nansum(nansum((total_entries)));
      fprintf('\nbranch 1: total_entries = %d (%.2f%% of total)\n\n', ...
        super_total_entries, 100*super_total_entries/num_y2s^4);
      
      
    case 2
      
      %   Branch II: y1 defines the starting point, y2 defines the dip and y4 defines the rise
      %   -1 <= y2 < y1, y2 is the dip so it must fall below the starting point (y1)
      %   y2 < y4 <= 1, y4 is the rise so it must fall above the dip (y2)
      %   y2 < y3 <= y4, y3 can hold any value between the dip (y2) and the rise (y4)
      %   y2 < y1 <= 1, y1 can hold any value that is above the dip (y2)
      
      y = -1:resolution:1;
      num_y2s = length(y);
      total_entries = nan(num_y2s,1);
      
      for r = 1:num_y2s
        
        if ~mod(r,floor(0.01*num_y2s))
          fprintf('%d%% ', ceil(100*r/num_y2s));
        end
        
        my_y2 = y(r);
        
        my_y1  = y(y > my_y2);
        my_y4 = y(y > my_y2);
        
        num_y1 = length(my_y1);
        num_y4 = length(my_y4);
        if ~num_y1; continue; end
        
        num_y3_and_y4 = sum( 0:(num_y2s-r-1) );
        
        %         num_y3 = nan(num_y4,1);
        %         for rr = 1:num_y4
        %           num_y3(rr) = length( y(y < my_y4(rr) & y > my_y2) );
        %           % NOTE: doing " < my_y4 " instead of " <= my_y4 " to avoid duplicate
        %           % curves from "case 1" where y3==y4
        %         end
        %         total_entries(r) = num_y1 * sum(num_y3);
        
        total_entries(r) = num_y1 * num_y3_and_y4;
        
      end
      
      super_total_entries = nansum(nansum((total_entries)));
      fprintf('\nbranch 2: total_entries = %d (%.2f%% of total)\n\n', ...
        super_total_entries, 100*super_total_entries/num_y2s^4);
      
      
    case 3
      
      %  Branch III: y1 defines the starting point, y3 defines the dip and y4 defines the rise
      %       -1 <= y3 < y4, y3 is the dip so it must fall below the rise (y4)
      %       y3 < y4 <= 1, y4 is the rise so it must fall above the dip (y3)
      %       y3 < y1 <= 1, y1 can hold any value that is above the dip (y3)
      %       y3 <= y2 <= 1, y2 can hold any value that is above the dip (y3)
      
      y = -1:resolution:1;
      num_y3s = length(y);
      total_entries = nan(num_y3s,1);
      
      for r = 1:num_y3s
        
        if ~mod(r,floor(0.01*num_y3s))
          fprintf('%d%% ', ceil(100*r/num_y3s));
        end
        
        my_y3 = y(r);
        
        my_y1  = y(y > my_y3);
        my_y2 = y(y > my_y3);
        my_y4 = y(y > my_y3);
        
        num_y1 = length(my_y1);
        num_y2 = length(my_y2);
        num_y4 = length(my_y4);
        if ~num_y1; continue; end
        
        total_entries(r) = num_y1 * num_y2 * num_y4;
        
      end
      
      super_total_entries = nansum(nansum((total_entries)));
      fprintf('\nbranch 3: total_entries = %d (%.2f%% of total)\n\n', ...
        super_total_entries, 100*super_total_entries/num_y3s^4);
      
      
    otherwise, error('Not a valid branch');
  end
  
end


