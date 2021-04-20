function process()
    % read a file 
    lower_limit = 0; 
    upper_limit = 23;
    disp('[INFO]: Scanning dir for *.txt data files');
    files = dir('./src/*.txt'); % count the number of file in the direcotry
    disp('[INFO]: Reading the first file in the directory');
    file_path = fullfile(files(1).folder, files(1).name); % first file 
    data_table = readFile(file_path, lower_limit, upper_limit);% read the first file
    disp('[INFO]; reading the rest of the files a for loop...');
    for i = 2:length(files)% start the loop at 2 since you have already read the first file.
        lower_limit = upper_limit + 1;
        upper_limit = upper_limit + 24; % increment by 24
        file_path = fullfile(files(i).folder, files(i).name);
        new_data_table = readFile(file_path, lower_limit, upper_limit);
        data_table = [data_table; new_data_table]; % append the new table to the existing table
    end
    disp('... Finished reading, starting  computation');
    % compute zhd_saas using the pressure values(Var3)
    zhd_saas = zhd_Saastamoinen(data_table.Var3);
    % compute zhd_hopj using the hopfield model
    zhd_hopj = zhd_Hopfield(data_table.Var4, data_table.Var3);
    % compute zhd_black using the black model
    zhd_black = zhd_Black(data_table.Var4, data_table.Var3);
    % compute for the rms values of the above models
    [zhd_saas_ms,zhd_hopj_ms, zhd_black_ms] = zhd_MS(zhd_saas, zhd_hopj, zhd_black);
    % compute the difference from zhd values
    [zhd_saas_e, zhd_hopj_e, zhd_black_e] = zhdDifference(data_table.Var2, zhd_saas, zhd_hopj, zhd_black);
    % compute for the value of e
    e_val = computeEval(data_table.Var4);
    % compute for the value of e
    zwd = computeZWD(e_val, data_table.Var3, data_table.Var4);
    
    disp('[INFO]: Done computing, plotting...')
    % plot the values of zhd from the models
    zwdPlots(data_table.data_3, zhd_saas, zhd_hopj, zhd_black);
    % plot rms values of zhd
    rmsPlots(data_table.data_3, zhd_saas_ms, zhd_hopj_ms, zhd_black_ms);
    disp('[END]: Saving the data in a file "output.txt" inside "./psd/"');
    % writing the data to a file
    new_table = table(zhd_saas,zhd_hopj, zhd_black, zhd_saas_ms,zhd_hopj_ms, zhd_black_ms,zhd_saas_e, zhd_hopj_e, zhd_black_e, zwd);
    data_table = [data_table new_table];
    writetable(data_table, './psd/output.txt','Delimiter','tab');
end

function data_table = readFile(file_path, lower_limit, upper_limit)
    f_id = fopen(file_path);
    data = textscan(f_id,"%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",24,'HeaderLines',4);
    data_3_temp = lower_limit:1:upper_limit;
    data_3 = data_3_temp.'; % transpose
    data_table = table(data_3,data{6}, data{11}, data{12});
    fclose(f_id);
end

function zhd_saas = zhd_Saastamoinen(pressure)
    % param: pressure is an double array  containing pressure values
    disp('[INFO]: Saastamoinen model - computing ZHD(saas)');
    zhd_saas = zeros(length(pressure), 1);
    % define constants 
    H = -0.0209202; % constant height
    lat =  -0.052288317664; % constant latitude
    % compute the value of zhd
    for i = 1:length(pressure)
        zhd_saas(i) = (0.0022767 * pressure(i)) / (1 - 0.00266* cosd(2 * lat) - 0.00028 * H);
    end
end

function zhd_hopj =  zhd_Hopfield(T_s, P_s)
    % param: T_s is an array of temperature values
    %        P_s  is an array of pressure values
    disp('[INFO]: Hopfield model - computing ZHD(hopj)');
    zhd_hopj = zeros(length(T_s), 1);
    % compute the value zhd_hopj using the Hopfield Model
    for i = 1:length(zhd_hopj)
        zhd_hopj(i) = (((10^(-6))/ 5)*(40136 + 148.72 * (T_s(i) - 273.16)) * 77.64 * (P_s(i)/ T_s(i)));
    end
end
function zhd_black = zhd_Black(T_s, P_s)
    % param: T_s is an array of temperature values
    % param: P_s is an array of pressure values
    disp('[INFO]: Black model - computing ZHD(black)');
    zhd_black = zeros(length(T_s),1);
    % compute the value of the zhd_black using the Black Model
    for i = 1:length(T_s)
        zhd_black(i) = ((0.0022767 * (T_s(i) - 4.12) * P_s(i))/ T_s(i));
    end
end

function [zhd_saas_ms,zhd_hopj_ms, zhd_black_ms] = zhd_MS(zhd_saas, zhd_hopj, zhd_black)
    % compure the MS value of the zhd_saas column
    disp('[INFO]: RMS values- computing ZHD RMS values for each model');
    % param: zhd_saas is a whole column
    syms n
    N_1 = length(zhd_saas);
    zhd_saas_ms  = sqrt(1/N_1) * symsum((zhd_saas.^2),n,1,N_1);
    N_2 = length(zhd_hopj);
    zhd_hopj_ms  = sqrt(1/N_2) * symsum((zhd_hopj.^2),n,1,N_2);
    N_3 = length(zhd_black);
    zhd_black_ms  = sqrt(1/N_3) * symsum((zhd_hopj.^2),n,1,N_3);
end

function [zhd_saas_e, zhd_hopj_e, zhd_black_e] = zhdDifference(zhd, zhd_saas, zhd_hopj, zhd_black)
    % using column  find the difference with the other values of zhd 
    disp('[INFO]: Computing difference(ZHD - ZHD(model))');
    zhd_saas_e = zhd - zhd_saas;
    zhd_hopj_e = zhd - zhd_hopj;
    zhd_black_e  = zhd - zhd_black;
end

function e_val = computeEval(T_s)
    % param: Ts: array of temperature values
    disp('[INFO]: Computing value e');
    e_val = zeros(length(T_s), 1);
    % compute the value of e using the value of Ts
    for i = 1:length(e_val)
        e_val(i) = 0.611 * exp((17.27 * T_s(i))/(T_s(i) + 273.3));
    end
end

function zwd = computeZWD(e_val, P_s,T_s)
    % param:  e_val: array of values of e
    % param : P_s; array of pressure values
    zwd = zeros(length(P_s), 1);
    disp('[INFO]; Computing zwd');
    % using Ifadis  model the value of ZWD can be computed using the value e 
    for i = 1:length(P_s)
        zwd(i) = ((0.554 * 10^(-2)) - (0.880*10^(-4))*(P_s(i) - 1000) + (0.272 * 10^(-4)) * e_val(i) + 2.771 * (e_val(i)/T_s(i)));
    end
end

function zwdPlots(time, zhd_saas, zhd_hopj, zhd_black)
    % plot the value of zhd from the three models
    % 1. Saastomoinen model
    % 2. Hopfield model
    % 3. Black model
    zhd_fig = figure('Name','ZHD');
    plot(time,zhd_saas,'r-', time, zhd_hopj, 'g-',time, zhd_black, 'b-');
    title("ZHD computed from 3 models");
    xlabel("time");
    ylabel("zhd from the models")
    legend("Saastomoinen","Hopfield","Black");
    saveas(zhd_fig,'./psd/zhd_fig.png');
end

function rmsPlots(time, zhd_saas_ms, zhd_hopj_ms, zhd_black_ms)
    % Plotting the RMS value of the ZHD computed from the three models
    rms_fig = figure('Name','RMS');
    plot(time, zhd_saas_ms, 'r-', time, zhd_hopj_ms, 'g-', time, zhd_black_ms,'b-');
    title("RMS values of ZHD computed from the 3 models");
    xlabel("Time hrs");
    ylabel("RMS values of ZHD");
    legend("Saastomoinen","Hopfield","Black");
    saveas(rms_fig, './psd/rms_fig.png');
end