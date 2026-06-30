clear;clc;close;
directory = "./data/fs_batch_data/";


Sd = 1000; gamy = 0;


Da1 = 1e0; %less permeable
Da2 = 1e4; %more permeable

for setup =1
    run_tag1 = sprintf('Sd_%.2e_Da_%.2e_gamy_%+.2e', Sd, Da1, gamy);
    run_tag1 = strrep(run_tag1, '+', 'p');
    run_tag1 = strrep(run_tag1, '-', 'm');
    folder1 = directory + run_tag1;
    
    run_tag2 = sprintf('Sd_%.2e_Da_%.2e_gamy_%+.2e', Sd, Da2, gamy);
    run_tag2 = strrep(run_tag2, '+', 'p');
    run_tag2 = strrep(run_tag2, '-', 'm');
    folder2 = directory + run_tag2;
end

for n = 1:150
    geoj1 = load(fullfile(folder1, sprintf("geo%d.mat", n)));
    geo1 = Geometry(geoj1.M, geoj1.P);

    geoj2 = load(fullfile(folder2, sprintf("geo%d.mat", n)));
    geo2 = Geometry(geoj2.M, geoj2.P);

    fprintf("Volume difference: %.9f\n", geo2.volume-geo1.volume);
end



