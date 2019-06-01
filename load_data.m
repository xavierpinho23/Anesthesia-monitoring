% Data loading
%zeus = load('G8Data/zeus_data_pat5_v2.mat');
trc = trc_file('G8Data/EEG_2.TRC');

% Data information
% SQI - quality of the signal

self = trc;
self.get_electrode_info;
self_electrodes =  self.a_file_elec_cell;
self_electrodes_frontal = self_electrodes([1:8 20]);
EEG = self.def_data_access(5,5, self_electrodes_frontal);

while self.a_is_last_segment ~=1
    EEG = [EEG self.get_next_window()];
end

save('EEG_frontal.mat','EEG')


