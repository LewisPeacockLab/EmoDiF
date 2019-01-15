function [] = compute_horz_indpt_vol_tw() %%% analog to JAL's compute_vert_rel_vol.m script.  replacing 'count particles' from family_of_curves%
%script output JLP 6 parameter calc to find chance criteron
%branch 1: total_entries = 22505647854382832 (14.06% of total)
%branch 2: total_entries = 30008265714771040 (18.75% of total)
%branch 3: total_entries = 23339111399098412 (14.58% of total)

%WRONG - USE JLF
  
  total = 0;
  
  total = total + family_of_curves('horz_indpnt','get_nParams')
  for k = 1:3
    total = total + compute_horz_indpt_vol_per_branch(k); % 6.1676e+16
  end
  
  resolution = 0.0001;
  total_possible = length(-1:resolution:1)^4;
  
  fprintf('\nTOTAL: total_entries = %d (%.2f%% of total)\n\n', ...
    total, 100*total/total_possible);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[super_total_entries] = compute_horz_indpt_vol_per_branch(branch)
  
  % resolution = 0.0001;
  resolution = 0.0001;
  num_chunks = 8;
  save_dir = '~/Dropbox (LewPeaLab)/STUDY/IMDiF/Results/P-CIT/vol_comp';
  
  switch branch
    case 1
      
% 		Branch I: y2 defines the dip and y3 defines the rise
% 		-1 <= y2 < 0, y2 is the dip so it must fall below zero
% 		0 < y3 <= 1, y3 is the rise so it must fall above zero
% 		-1 <= y4 <= y3, y4 can hold any value that is below the rise (y3)
% 		y2 < y1 <= 1, y1 can hold any value that is above the dip (y2)

      
      y = -1:resolution:1; %resolution of y axis
      total_curve_num = length(y);
      abovey = 0:resolution:1;
      belowy = -1:resolution:0;
      
      base_num = length(belowy);
      
      %num_y2s = length(y); %number of Y values
      total_entries = nan(base_num,1); %creating matrix
      
      for r = 1:length(belowy) %y2 is constrained as below 0
        
        if ~mod(r,floor(0.01*base_num)) %this gives out percentage calculated
          fprintf('%d%% ', ceil(100*r/base_num));
        end
    
        my_y2 = belowy(r); %gets values where y1 > y2
        my_y1 = y(y> my_y2);
        
        %my_y3 = y(y>0)
        
        nums_y3_y4 = [];
        
        for t = 1:length(abovey)
            
        my_y3 = abovey(t); % gets values where y3>y2 (with the prior line, this creates a dip)
        my_y4 = y(y < my_y3); %gets values where -1 <= y4 <= y3
        

        num_y3 = length(my_y3);
        num_y4 = length(my_y4);
        
        num_y3_y4 = num_y3 * num_y4;
        
        nums_y3_y4 = horzcat(num_y3_y4, nums_y3_y4);
        
        end
                num_y1 = length(my_y1);
        if ~num_y1; continue; end 
%         
%         num_y3_and_y4 = sum( (r+1):base_num );
        total_entries(r) = num_y1 * sum(nums_y3_y4);
      end
      
      super_total_entries = nansum(nansum((total_entries)));
      fprintf('\nbranch 1: total_entries = %d (%.2f%% of total)\n\n', ...
        super_total_entries, 100*super_total_entries/total_curve_num^4);
  
      
    case 2
      
% 		Branch II: y2 defines the dip and y4 defines the rise
% 		-1 <= y2 < 0, y2 is the dip so it must fall below zero
% 		0 < y4 <= 1, y4 is the rise so it must fall above zero
% 		y2 <= y3 <= y4, y3 can hold any value between the dip and the rise
% 		y2 < y1 <= 1, y1 can hold any value that is above the dip (y2)
      
    y = -1:resolution:1; %resolution of y axis
      total_curve_num = length(y);
      abovey = 0:resolution:1;
      belowy = -1:resolution:0;
      
      base_num = length(belowy);
      
      %num_y2s = length(y); %number of Y values
      total_entries = nan(length(belowy),1); %creating matrix
      
      for r = 1:length(belowy) %y2 is constrained as below 0
        
        if ~mod(r,floor(0.01*base_num)) %this gives out percentage calculated
          fprintf('%d%% ', ceil(100*r/base_num));
        end
    
        my_y2 = belowy(r); %gets values where y1 > y2
        my_y1 = y(y> my_y2); % 		y2 < y1 <= 1, y1 can hold any value that is above the dip (y2)
        
        %my_y3 = y(y>0)
        
        my_y3 = y; % 		y2 <= y3 <= y4, y3 can hold any value between the dip and the rise
        
        my_y4 = abovey;% rise gets value above 0
        
       

        num_y3 = length(my_y3);
        num_y4 = length(my_y4);
        
        num_y3_y4 = num_y3 * num_y4;
        
                num_y1 = length(my_y1);
                
        if ~num_y1; continue; end 
%         
%         num_y3_and_y4 = sum( (r+1):base_num );

        total_entries(r) = num_y1 * (num_y3 * num_y4);
  
      end
      
      super_total_entries = nansum(nansum((total_entries)));
      fprintf('\nbranch 2: total_entries = %d (%.2f%% of total)\n\n', ...
        super_total_entries, 100*super_total_entries/total_curve_num^4);
      
      
    case 3
      
% 		Branch III: y3 defines the dip and y4 defines the rise
% 		-1 <= y3 < 0, y3 is the dip so it must fall below zero
% 		0 < y4 <= 1, y4 is the rise so it must fall above zero
% 		y3 < y1 <= 1, y1 can hold any value that is above the dip (y3)
% 		y3 <= y2 <= 1, y2 can hold any value that is above the dip (y3)
% 	
%       
 y = -1:resolution:1; %resolution of y axis
      total_curve_num = length(y);
      abovey = 0:resolution:1;
      belowy = -1:resolution:0;
      
      
      %num_y2s = length(y); %number of Y values
      total_entries = nan(length(belowy),1); %creating matrix
      
      
      for r = 1:length(belowy)
        
        if ~mod(r,floor(0.01*length(belowy)))
          fprintf('%d%% ', ceil(100*r/length(belowy)));
        end
        
        my_y3 = belowy(r);
        
        my_y1  = y(y > my_y3);
        my_y2 = y(y > my_y3);
        my_y4 = abovey;
        
        num_y1 = length(my_y1);
        num_y2 = length(my_y2);
        num_y4 = length(my_y4);
        if ~num_y1; continue; end
        
        total_entries(r) = num_y1 * num_y2 * num_y4;
        
      end
      
      super_total_entries = nansum(nansum((total_entries)));
      fprintf('\nbranch 3: total_entries = %d (%.2f%% of total)\n\n', ...
        super_total_entries, 100*super_total_entries/total_curve_num^4);
      
      
    
  end
  
end


