function[t_out]= time_to_necking(directory,tc)

files = dir(fullfile(directory, 'geo*.mat'));
tf = int32(length(files));
t = int32(tf/10);

will = [];
for t = tf-100:1:tf
    geoj = load(directory + sprintf('geo%i.mat',t));
    geo = Geometry(geoj.M,geoj.P);
    will = [will,geo.willmore_energy(1) - 8*pi];
end



dwill = will(2:end)-will(1:end-1)/length(will);
t_out = (tf-tc) + will(1)/mean(dwill);
end
%t_doubled = willmore(t_cup)/((willmore(t_cup) - willmore(t_cup + N))/N);


