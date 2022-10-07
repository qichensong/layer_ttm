function a = get_data(fname)
s = importdata(fname,'\t',2);
a = s.data;
end