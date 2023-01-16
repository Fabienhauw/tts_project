res_dir = '/network/lustre/iss02/cohen/data/Fabien_official/SYNESTHEX/second_level/Aud/loc/FIR_analyses_rh';
cd(res_dir)
results = dir('fir_analysis*.mat');
for tmp_res = 1 : length(results)
    load(results(tmp_res).name)
    counter = 0;
    clear res_mat
    for tmp_con = 1 : size(save_sfiry,2)
        for tmp_group = 1 : size(save_sfiry,1)
            for tmp_sub = 1 : size(save_sfiry{tmp_group,tmp_con},1)
                for tmp_time = 1 : size(save_sfiry{tmp_group,tmp_con},2)
                    counter = counter + 1;
                    res_mat(counter).firy = save_sfiry{tmp_group,tmp_con}(tmp_sub, tmp_time);
                    res_mat(counter).time = timeaxis(tmp_time);
                    res_mat(counter).subj = tmp_sub;
                    res_mat(counter).cond = tmp_con;
                    res_mat(counter).group = tmp_group;
                end
            end
        end
    end
    
    Firy = [res_mat.firy]'; Time = [res_mat.time]'; Subj_Id = [res_mat.subj]';
    Cond = [res_mat.cond]'; Group_Id = [res_mat.group]';
    table_result = table(Firy, Time, Subj_Id, Cond, Group_Id);
    % Write data to text file
    writetable(table_result, fullfile(res_dir, sprintf('mat_for_stats_%s.txt', results(tmp_res).name(1:end-4))))
end