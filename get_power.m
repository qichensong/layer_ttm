function a = get_power(fname)
s = importdata(fname,'\t',1);
s = split(s);
a = str2double(s{length(s)-1})*0.001;
end