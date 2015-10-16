function output = day_stamp( byr, bmon, bday, eyr, emon, eday, ops )
%   Create time stamp at daily step
%   input: byr  - first year  [1x1]
%          bmon - first month [1x1]
%          bday - first day   [1x1]
%          eyr  - last year   [1x1]
%          emon - last month  [1x1]
%          eday - last day    [1x1]
%          ops  - 1 => output = date number 
%                 2 => output = date vector [yy mm dd]

% Create daily array
tnum1= datenum(byr,bmon,bday);
tnum2= datenum(eyr,emon,eday);
dnum_day = (tnum1:1:tnum2)';
[yrs, mon ,day] = datevec(dnum_day);

if ops == 1
    output = dnum_day;
elseif ops ==2
    output = [yrs, mon ,day];
end

end
