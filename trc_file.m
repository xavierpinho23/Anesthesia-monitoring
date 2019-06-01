classdef trc_file < handle
    %======================================================================
    %Matlab class that defines a "trc_file" object, enabling TRC parameter
    %and data access. For data access, in first place the method
    %"def_data_access" must be called. This method defines the size of
    %the data window, the time step between data windows, and the
    %considered channels. Then the method "get_next_window" can be called
    %in order to retrive the next data window in the file. A method that
    %retireves only the next non-overlaping portion of data was also
    %develloped. This method is called "get_next_segment".The parameters
    %defined previously in "def_data_access" are considered in the 
    %"get_next_window" and "get_next_segment" methods.
    %
    %EU FP7 Grant 211713 (EPILEPSIAE)
    %
    %César A. D. Teixeira
    %CISUC, FCTUC, University of Coimbra
    %September 2009
    %======================================================================
    
    
     properties
        
        a_n_chan; % number of channels
        a_samp_freq; % sampling rate
        a_n_bytes; % number of bytes per sample in the file
        a_n_samples; % total number of samples per channel
        a_n_data_secs; %  total number of seconds per channel
        a_duration_ts; % dd HH:MM:SS.fff
        a_start_ts; % yyyy-mm-dd HH:MM:SS.fff
        a_stop_ts; % yyyy-mm-dd HH:MM:SS.fff
        a_file_elec_cell; % cell array containing the names of all electrodes presented in the file
        a_conv_factor; % conversion factor
        a_note_struct; % structure array containing the notes and the related sample number

        a_step; % step between data windows (in seconds)
        a_wsize; % window size in seconds
        a_channs_cell; % cell array containing the names of electrodes to be processed. The name can contain channels in MONOPOLAR, BIPOLAR or AVERAGE montage
        a_channs_ave_cell; % cell array containing the names of electrodes to be included in the average. Default value: empty
        a_first_data_bytes_offset; % offset in bytes of the first sample to be read
        a_first_data_samp_offset; % offset in samples of the first sample to be read
        a_curr_bytes_offset; % offset in bytes from the beginning of the file to the beginning of the current window
        a_last_bytes_offset; % offset in bytes from the beginning of the file to the end of the current window
        a_curr_samp_offset; % offset in samples from the beginning of the file to the beginning of the current window
        a_last_samp_offset; % offset in samples from the beginning of the file to the end of the current window
        a_is_first_segment; % flag to indicate first complete segment
        a_is_last_segment; % flag to indicate last complete segment
        a_last_data_mat; % matrix containing the last data window
        a_last_flagg=0;
        
        
    end%End propreties
    
    properties (SetAccess = protected)
        a_file_name;
        a_step_samp; % step between data windows (in samples)
        a_channs_ave_ind; % cell array containing the names of electrodes to be included in the average
        
        a_acq_unit;%Acquisition unit ID
        a_avg_ar_len;
        a_avg_ar_start;
        a_b_imped_ar_len;
        a_b_imped_ar_start;
        a_code_ar_len;
        a_code_ar_start;
        a_comp;
        a_comp_ar_len;
        a_comp_ar_start;
        a_day_rec;
        a_desc;
        a_dvid_ar_len
        a_dvid_ar_start
        a_dvid_start
        a_e_imped_ar_len
        a_e_imped_ar_start
        a_elect_ar_len
        a_elect_ar_start
        a_elect_list
        a_ev_a_ar_len
        a_ev_a_ar_start
        a_ev_b_ar_len
        a_ev_b_ar_start
        
        a_fd
        a_file_type
        
        a_flag_ar_len
        a_flag_ar_start
        a_head_type
        a_hist_ar_len
        a_hist_ar_start
        a_hour_rec
        a_lab
        a_min_rec
        a_montage
        a_montage_ar_len
        a_mont_ar_start
        a_mouth_rec
        a_mpeg_dsync
        a_mux
        
        
        a_n_rec_sec
        a_note_ar_len
        a_note_ar_start
        a_order
        a_red_ar_len
        a_red_ar_start
        
       
        
        a_sec_rec
        a_trig_ar_len
        a_trig_ar_start
        a_year_rec
        a_mes_unit_cell
        a_unit_conv_cell
        a_data_map
        a_dtype
       
        a_abs_file_pt
        a_n_bytes_ahead
        
        a_n_bytes_win
        a_is_subset
        
        
    end%End propreties
    
    
    
    
    
    
    
    
   
    methods
        function self=trc_file(f_name)
            %==============================================================
            %"trc_file" class constructor
            %In order to create an instance of trc_file class, only the TRC
            %file name is required.
            %
            %Input:
            %   f_name-->TRC file name
            %
            %Output:
            %   self-->trc_file class object
            %
            %
            %EU FP7 Grant 211713 (EPILEPSIAE)
            %
            %César A. D. Teixeira
            %CISUC, FCTUC, University of Coimbra
            %September 2009
            %==============================================================

            self.a_file_name=f_name;
            self.a_fd=fopen(self.a_file_name,'r');
            self.a_desc=reshape(fread(self.a_fd,32,'*char'),1,32);%Get TRC description
            self.a_lab=reshape(fread(self.a_fd,32,'*char'),1,32);%Get lab description
            fseek(self.a_fd,64,'cof');%Go 64 bytes ahead
            %===================Get Record date==========================
            self.a_day_rec=fread(self.a_fd,1,'uchar');
            self.a_mouth_rec=fread(self.a_fd,1,'uchar');
            self.a_year_rec=fread(self.a_fd,1,'uchar')+1900;
            %===================Get Record time==========================
            self.a_hour_rec=fread(self.a_fd,1,'uchar');
            self.a_min_rec=fread(self.a_fd,1,'uchar');
            self.a_sec_rec=fread(self.a_fd,1,'uchar');
            %============Compute the number of recorded seconds==========
            self.a_n_rec_sec=self.a_sec_rec+(self.a_min_rec*60)+(self.a_hour_rec*3600);
            %================Get Acquisition Unit description============
            self.a_acq_unit=fread(self.a_fd,1,'uint16');
            %======================Get File Type=========================
            self.a_file_type=fread(self.a_fd,1,'uint16');
            %====================Read Data start offset==================
            self.a_first_data_bytes_offset=fread(self.a_fd,1,'uint32');
            %=============Get Number of Memorised Channels===============
            self.a_n_chan=fread(self.a_fd,1,'uint16');
            %=============Get Distance Between Channel Data==============
            self.a_mux=fread(self.a_fd,1,'uint16');
            %==============Get Minimum Sampling Frequency================
            self.a_samp_freq=fread(self.a_fd,1,'uint16');
            %==========Get Number of Bytes per sample channel============
            self.a_n_bytes=fread(self.a_fd,1,'uint16');
            %=============Get Type of Data Compression===================
            self.a_comp=fread(self.a_fd,1,'uint16');
            %============Get Number of Specific Montages=================
            self.a_montage=fread(self.a_fd,1,'uint16');
            %========Get the Starting sample of digital video============
            self.a_dvid_start=fread(self.a_fd,1,'uint32');
            %Get Number of frames per hour of de-synchronization in MPEG acq.
            self.a_mpeg_dsync=fread(self.a_fd,1,'uint16');
            %=================Move 15 bytes ahead=======================
            fseek(self.a_fd,15,'cof');
            %==================Get Header type==========================
            self.a_head_type=fread(self.a_fd,1,'uint8');

            %==================Get area blocks and sizes=================

            area_blk_names={'self.a_code_ar_start','self.a_elect_ar_start','self.a_note_ar_start','self.a_flag_ar_start','self.a_red_ar_start',...
                'self.a_b_imped_ar_start','self.a_e_imped_ar_start','self.a_mont_ar_start','self.a_comp_ar_start','self.a_avg_ar_start',...
                'self.a_hist_ar_start','self.a_dvid_ar_start','self.a_ev_a_ar_start','self.a_ev_b_ar_start','self.a_trig_ar_start'};

            area_blk_len={'self.a_code_ar_len','self.a_elect_ar_len','self.a_note_ar_len','self.a_flag_ar_len','self.a_red_ar_len',...
                'self.a_b_imped_ar_len','self.a_e_imped_ar_len','self.a_montage_ar_len','self.a_comp_ar_len','self.a_avg_ar_len',...
                'self.a_hist_ar_len','self.a_dvid_ar_len','self.a_ev_a_ar_len','self.a_ev_b_ar_len','self.a_trig_ar_len'};

            for i=1:size(area_blk_names,2)
                fseek(self.a_fd,8,'cof');
                eval([area_blk_names{i},'=fread(self.a_fd,1,''uint32'');']);
                eval([area_blk_len{i},'=fread(self.a_fd,1,''uint32'');']);
            end

            %===================Go 224 bytes ahead======================
            fseek(self.a_fd,224,'cof');

            %==================Get electrode order======================
            %self.order=fread(self.a_fd,2*self.a_n_chan,'uint8');
            %self.order=self.order(1:end-1,1);

            %=============Memmapfile object for data access=============
            self.a_data_map=memmapfile(self.a_file_name,'Offset',self.a_first_data_bytes_offset,'Repeat',1);
            self.a_is_last_segment=0;
            self.a_conv_factor=1;
            
            %==================Compute data length========================
            f_inf=dir(self.a_file_name);
            data_size=(f_inf.bytes)-self.a_first_data_bytes_offset;
            self.a_n_samples=data_size/self.a_n_bytes/self.a_n_chan;
            
            self.a_abs_file_pt=self.a_first_data_bytes_offset+(self.a_n_samples*self.a_n_bytes*self.a_n_chan);
            
            self.a_n_data_secs=self.a_n_samples/self.a_samp_freq;
            
            self.a_start_ts=datestr([self.a_year_rec,self.a_mouth_rec,self.a_day_rec,self.a_hour_rec,self.a_min_rec,self.a_sec_rec],'yyyy-mm-dd HH:MM:SS');
            seq_date_i=datenum(self.a_year_rec,self.a_mouth_rec,self.a_day_rec,self.a_hour_rec,self.a_min_rec,self.a_sec_rec);
            seq_inc=datenum(0,0,0,0,0,self.a_n_data_secs);
            seq_date_fin=seq_date_i+seq_inc;
            self.a_stop_ts=datestr(seq_date_fin,'yyyy-mm-dd HH:MM:SS');
            self.a_duration_ts=datestr(seq_inc,'DD HH:MM:SS');
            
            %================Determine if the TRC file is a subset from
            %other big file
            
            fseek(self.a_fd,self.a_red_ar_start,'bof');
            if fread(self.a_fd,1,'ulong')~=0
                self.a_is_subset=1;
            else
                self.a_is_subset=0;
            end
            %================Set data access flags==================
            self.a_is_last_segment=0;
            self.a_is_first_segment=1;
            self.a_curr_bytes_offset=self.a_first_data_bytes_offset;%Define the current data offset equal to the offset of the first data sample
            self.a_curr_samp_offset=0;
            
        end%End Constructor

        
        function get_notes(self)
            %==============================================================
            %"get_notes" method
            %This method acquires all the information related with the anotations
            %This information is:
            %   -sample where the event occurs
            %   -text describing the associated event
            %
            %This information is stored in the "a_note_struct" attribute that
            %is a struct array with fields 'sample' and 'note'
            %
            %
            %EU FP7 Grant 211713 (EPILEPSIAE)
            %
            %César A. D. Teixeira
            %CISUC, FCTUC, University of Coimbra
            %September 2009
            %==============================================================
            fseek(self.a_fd,self.a_note_ar_start,'bof');%Go to the note area in the file
            self.a_note_struct=[];%Initialize the 'a_note_struct' attribute
            while 1
                sample=fread(self.a_fd,1,'uint32');%Get the current sample number
                if sample~=0%A null sample means that there is no more notes
                    note=reshape(fread(self.a_fd,40,'*char'),1,40);%Get the current event note
                    ixc = isstrprop(note, 'graphic');%Look for non-graphic characters at the end of the note string and eliminate them
                    ixc=find(ixc);
                    if size(ixc,2)~=0
                        ixc=ixc(end);
                        note=note(1:ixc);
                    end
                    self.a_note_struct=[self.a_note_struct,struct('sample',sample,'note',note)];%Pack sample number and note in
                    %a struct and append it to the 'a_note_struct' array
                else
                    break;
                end
            end
        end%End get_notes method
        
        function get_elect_order(self)
            %==============================================================
            %"get_elect_order" method
            %This method acquires the information related with the
            %electrode order
            %This method returned and array, which contains the coded list
            %of the electrodes stored in the file
            %This information is stored in the "order" attribute
            %
            %
            %EU FP7 Grant 211713 (EPILEPSIAE)
            %
            %César A. D. Teixeira
            %CISUC, FCTUC, University of Coimbra
            %September 2009
            %==============================================================
            fseek(self.a_fd,self.a_code_ar_start,'bof');%Go to the code area in the file
            self.a_order=zeros(1,self.a_n_chan);
            for i=1:self.a_n_chan
                self.a_order(i)=fread(self.a_fd,1,'uint16');%Get the current order
            end
        end%End get_elect_order method
        function get_electrode_info(self)
            %==============================================================
            %"get_electrode_info" method
            %This method acquires all the information related with
            %the electrodes.
            %This information is:
            %   -status of the electrode, i.e. On/Off
            %   -type
            %   -positive impedance label
            %   -negative impedance label
            %   -logical minimum
            %   -logical maximum
            %   -logical ground
            %   -physical minimum
            %   -physical maximum
            %   -measurement unit (nV, uV, etc)
            %   -High-pass pre-filtering cutt-off frequency
            %   -High-pass prefiletering type 
            %   -Low-pass pre-filtering cutt-off frequency
            %   -Low-pass prefiletering type
            %   -Sampling rate coeficient 1=minimum, 2=2*minimum, and 4=4*minimum
            %   -Maps latitude (radius=1)
            %   -Maps longitude (radius=1)
            %   -Is the electrode present in map
            %   -Is the electrode used for AVG calculation
            %   -Electrode description
            %   -x coordinate for 3D map analysis
            %   -y coordinate for 3D map analysis
            %   -z coordinate for 3D map analysis
            %
            %This information is retrieved for each channel, and will be
            %used for data conversion in the data access methods.
            %
            %
            %EU FP7 Grant 211713 (EPILEPSIAE)
            %
            %César A. D. Teixeira
            %CISUC, FCTUC, University of Coimbra
            %September 2009
            %==============================================================
            self.get_elect_order();
            

            %self.a_elect_list=[]#List that will contain a sequence of electrode class instances
            self.a_mes_unit_cell=cell(7,2);
            self.a_unit_conv_cell=cell(7,2);
            self.a_mes_unit_cell={'i-1','nVolt';'i0','uVolt';'i1','mVolt';'i2','Volt';'i100','%';'i101','bpm';'i102','Adim'};%cell of measurement units
            self.a_unit_conv_cell={'i-1','1e-3';'i0','1';'i1','1000';'i2','1e6';'i100','1';'i101','1';'i102','1'};%cell of measurement units
            self.a_elect_list=[];
            
         

            for i=1:self.a_n_chan
                off=self.a_elect_ar_start+(self.a_order(i)*128);%Upgrade the file offset for the related channel order
                fseek(self.a_fd,off,'bof');%Point to the region where information related with electrode i is stored
                status=fread(self.a_fd,1,'uint8');%Get electrode i status
                type=fread(self.a_fd,1,'uint8');%Get electrode i Type
                pos_imp_label=reshape(fread(self.a_fd,6,'*char'),1,6);%Get positive impedance label
                neg_imp_label=reshape(fread(self.a_fd,6,'*char'),1,6);%Get negative impedance label
                
                ixc = isstrprop(pos_imp_label, 'graphic');%Look for non-graphic characters at the end of the positive impedance label string and eliminate them
                ixc=find(ixc);
                if size(ixc,2)~=0
                     ixc=ixc(end);
                     pos_imp_label=pos_imp_label(1:ixc);
                end
                
                ixc = isstrprop(neg_imp_label, 'graphic');%Look for non-graphic characters at the end of the negative impedance label string and eliminate them
                ixc=find(ixc);
                if size(ixc,2)~=0
                     ixc=ixc(end);
                     neg_imp_label=neg_imp_label(1:ixc);
                end
                
                l_min=fread(self.a_fd,1,'int32');%Get logic minimum
                l_max=fread(self.a_fd,1,'int32');%Get logic maximum
                l_grd=fread(self.a_fd,1,'int32');%Get logic ground
                p_min=fread(self.a_fd,1,'int32');%Get physical minimum
                p_max=fread(self.a_fd,1,'int32');%Get physical maximum
                bin_unit=fread(self.a_fd,1,'uint16');%Get measure self.eeg_conv_factorment Unit and conversion factor
                mes_unit=self.a_mes_unit_cell{find(ismember(self.a_mes_unit_cell,['i',num2str(bin_unit)])==1),2};%Define the measurement unit
                conv_unit=str2num(self.a_unit_conv_cell{find(ismember(self.a_unit_conv_cell, ['i',num2str(bin_unit)])==1),2});%Define the decimal factor, acourding to the related unit
                pre_high_pass_lim=fread(self.a_fd,1,'uint16')/1000;%High-pass pre-filtering cutt-off frequency
                pre_high_pass_type=fread(self.a_fd,1,'uint16');%High-pass prefiletering type 
                pre_low_pass_lim=fread(self.a_fd,1,'uint16');%Low-pass pre-filtering cutt-off frequency
                pre_low_pass_type=fread(self.a_fd,1,'uint16');%Low-pass prefiletering type
                samp_rate_coef=fread(self.a_fd,1,'uint16');%Sampling rate coeficient 1=minimum, 2=2*minimum, and 4=4*minimum
                fseek(self.a_fd,2,'cof');%Go 2 bytes ahead (reserved area)
                latitudine=fread(self.a_fd,1,'float32');%Maps latitude (radius=1)
                longitudine=fread(self.a_fd,1,'float32');%Maps longitude(radius=1)
                is_present_in_map=fread(self.a_fd,1,'uchar');%Presence of electrode in map
                is_used_in_avg=fread(self.a_fd,1,'uchar');%Used for AVG calculation
                description=reshape(fread(self.a_fd,32,'*char'),1,32);%Electrode description
                x_coord=fread(self.a_fd,1,'float32');%x coordinate for 3D map analysis
                y_coord=fread(self.a_fd,1,'float32');%y coordinate for 3D map analysis
                z_coord=fread(self.a_fd,1,'float32');%z coordinate for 3D map analysis
                coord_type=fread(self.a_fd,1,'uint16');%Get coordinate type
                if coord_type==0
                    coord_type='Polar map (Lat & long)';
                else
                    coord_type='X Y Z (3D System)';
                end
                elect=struct('status',status,'type',type,'pos_imp_label',pos_imp_label,'neg_imp_label',neg_imp_label,...
                    'l_min',l_min,'l_max',l_max,'l_grd',l_grd,'p_min',p_min,'p_max',p_max,'mes_unit',mes_unit,'conv_unit',conv_unit,...
                    'pre_high_pass_lim',pre_high_pass_lim,'pre_high_pass_type',pre_high_pass_type,'pre_low_pass_lim',pre_low_pass_lim,'pre_low_pass_type',pre_low_pass_type,...
                    'samp_rate_coef',samp_rate_coef,'latitudine',latitudine,'longitudine',longitudine,'is_present_in_map',is_present_in_map,'is_used_in_avg',is_used_in_avg,...
                     'description',description,'x_coord',x_coord,'y_coord',y_coord,'z_coord',z_coord,'coord_type',coord_type);%Integrate all the electrode parameters in a struct
                self.a_elect_list=[self.a_elect_list,elect];%Pack the information of all electrodes in a struct array
            end
            self.a_file_elec_cell={self.a_elect_list(:).pos_imp_label};%Put electrode labels in a cell array

        end%End get_electrode_info method
        
        function set_elect_ave(self, ave_elec_cell)
            %%
            %==============================================================
            % This functions sets the indices of electrodes to be
            % considered for an optional common average montage
            %
            % Input:
            %   ave_elec_cell-->cell array containing the names of
            %   electrodes to be included in the average
            %
            % Output:
            %
            %%EU FP7 Grant 211713 (EPILEPSIAE)
            %
            %César A. D. Teixeira
            %CISUC, FCTUC, University of Coimbra
            %November 2009
            %==============================================================
                  
            
            self.a_channs_ave_cell = ave_elec_cell;
            self.a_channs_ave_ind = zeros(1, length(self.a_channs_ave_cell));
            
            for i=1:length(self.a_channs_ave_cell)
                
                    idx=find(ismember(self.a_file_elec_cell,self.a_channs_ave_cell{i})==1);
                    
                    if ~isempty(idx)
                        self.a_channs_ave_ind(i)=idx;
                    end
                    
            end
            
            
        end
        
        function [cdata]=combine_data(self,data)
            %"combine_data" method
            %This method combine data channels acourding to the cell array stored
            %by the a_channs_cell attribute
            %   The a_channs_cell attribute ia a cell array containing the name of the channels
            %       to be processed. The name can contain channels in
            %       MONOPOLAR, BIPOLAR or AVERAGE montage. A combination from
            %       some of them is also possible.
            %       Example:
            %           {'FP1';'FP2';'FPZ-FZ'}
            %            |__________|_______|
            %              Mono      Bipolar --> 2 channels Monopolar, 1 channel Bipolar
            %       For AVERAGE montage, a cell array containing the name of
            %       the channels to be included in the average must be
            %       declared before through the method set_elect_ave
            %
            %       
            %Outputs:
            %   data-->Matrix containing the related data
            %
            %
            %EU FP7 Grant 211713 (EPILEPSIAE)
            %
            %César A. D. Teixeira
            %CISUC, FCTUC, University of Coimbra
            %November 2009
            %==============================================================
            
            nsamples=size(data,1);
            
            cdata=zeros(nsamples,length(self.a_channs_cell));
            
            
            for i=1:length(self.a_channs_cell)
                [sig1 sig2] = strtok(self.a_channs_cell{i}, '-');

                sig1_idx=find(ismember(self.a_file_elec_cell,sig1)==1);
                if isempty(sig1_idx)
                    cdata=[];
                else

                    if ~isempty(sig2)
                        if ~strcmpi(sig2(2:end),'avg')
                            
                            sig2_idx=find(ismember(lower(self.a_file_elec_cell),lower(sig2(2:end)))==1);

                            if ~isempty(sig2_idx)
                                cdata(:,i)=data(:,sig1_idx)-data(:,sig2_idx);
                            else
                                cdata(:,i)=data(:,sig1_idx);
                            end

                        else
                            length(self.a_channs_ave_ind)
                            
                            data_m=sum(data(:,self.a_channs_ave_ind),2)./length(self.a_channs_ave_ind);
                            cdata(:,i)=data(:,sig1_idx)-data_m;
                        end
                    else
                        cdata(:,i)=data(:,sig1_idx);
                    end
                end
            end
                
            
            
            
            
            
        end
        
        

        function [data,time]=def_data_access(self, wsize, step, channs_cell, offset)
            %==============================================================
            %"def_data_access" method
            %This method defines the way as EEG/ECG data will be accessed.
            %Parameters like data window size, time steb between data
            %windows and the considered channels are defined by this
            %method.
            %This method returns a matrix with the first data window. The
            %next data windows are obtained by the "get_next_window"
            %method.
            %
            %Inputs:
            %   step-->Time in seconds between windows. A step equal to the
            %   window size resulted in non-overlaped windows.
            %   wsize-->Window size in seconds
            %   channs_cell-->Cell array containing the name of the channels
            %       to be processed. The name can contain channels in
            %       MONOPOLAR, BIPOLAR or AVERAGE montage. A combination from
            %       some of them is also possible.
            %       Example:
            %           {'FP1';'FP2';'FPZ-FZ'}
            %            |__________|_______|
            %              Mono  Bipolar --> 2 channels Monopolar, 1 channel Bipolar
            %       For AVERAGE montage, a cell array containing the name of
            %       the channels to be included in the average must be
            %       declared before through the method set_elect_ave
            %
            %       
            %Outputs:
            %   data-->Matrix containing the related data
            %
            %
            %EU FP7 Grant 211713 (EPILEPSIAE)
            %
            %César A. D. Teixeira
            %CISUC, FCTUC, University of Coimbra
            %September 2009
            %==============================================================
            
            
            if ~exist('wsize', 'var') || isempty(wsize)
                return;
            end
            
            self.a_wsize = wsize;
            if ~exist('step', 'var') || isempty(step)
                self.a_step = wsize;
            else
                self.a_step = step;
            end       
            if ~exist('offset', 'var') || isempty(offset)
                self.a_first_data_samp_offset = 0;
            else
                self.a_first_data_samp_offset = offset;
            end
            if ~exist('channs_cell', 'var') || isempty(channs_cell)
                self.a_channs_cell = self.a_file_elec_cell;
            else
                self.a_channs_cell = channs_cell;
            end        
            
            
            
            
            
            
            
            
            self.a_n_bytes_win=floor(self.a_wsize*self.a_samp_freq)*self.a_n_chan*self.a_n_bytes;
            
            
            self.a_step_samp=floor(self.a_step*self.a_samp_freq);%Compute the step in samples
            
            self.a_n_bytes_ahead=self.a_step_samp*self.a_n_chan*self.a_n_bytes;%Compute the step in bytes
            
            
            data_types={'1','uint8';'2','uint16';'4','uint32'};%Data type.
            
            
            self.get_electrode_info();%Call the "Get electrode" method for logical and physical gains retrieval

            self.a_conv_factor=(self.a_elect_list(1).p_max-self.a_elect_list(1).p_min)/(self.a_elect_list(1).l_max-self.a_elect_list(1).l_min+1)...
                *self.a_elect_list(1).conv_unit;% Compute EEG conversion factor
            
             nsamples=floor(self.a_wsize*self.a_samp_freq);%Compute the number of samples to be acquired
            
            
            
            %self.a_curr_bytes_offset=self.a_curr_bytes_offset+self.a_n_bytes_win;%Define the offset of the last data sample in the first window
            
            self.a_curr_bytes_offset=self.a_first_data_bytes_offset+(self.a_first_data_samp_offset *self.a_n_bytes * self.a_n_chan);
            self.a_curr_samp_offset=self.a_first_data_samp_offset;
            
            
            
            
            
            self.a_last_bytes_offset=self.a_curr_bytes_offset+self.a_n_bytes_win;%Define the offset of the last data sample in the first window
            self.a_last_samp_offset=self.a_curr_samp_offset+nsamples;
            
            
            self.a_data_map.Offset=self.a_curr_bytes_offset;%Point the memmapfile object to the current data offset
            
            
           
            

            dtype_idx=find(ismember(data_types,num2str(self.a_n_bytes))==1);%Get the appropriate data type
            self.a_dtype=data_types{dtype_idx,2};
            

            self.a_data_map.Format={self.a_dtype [self.a_n_chan nsamples] 'dt'};%Define the data format that the memmapfile object should consider

            try%Try if it is possible to acquire the requested amount of data
                data=(self.a_data_map.Data(1).dt)';
                %data=data(:,self.a_channs_cell);
                self.a_is_last_segment=0;
                
                if self.a_curr_samp_offset==0
                    self.a_is_first_segment=1;
                else
                    self.a_is_first_segment=0;
                end
                
%                 self.curr_off=self.curr_off+nsamples*(self.a_n_chan*self.a_n_bytes);%Bruno

            
            catch%If an error is raised, then redifine the memmapfile object parameters to adjust to the remaining data in file.
                self.a_data_map.Format=self.a_dtype;%Consider independent data samples, not arranged acourding to the number of channels
                self.a_data_map.Repeat=Inf;%Repeat the data format to the end of file
                last_data=self.a_data_map.Data;%Get data
                n_last_samp=size(last_data,1);%Compute how many individual samples were readed
                data=reshape(last_data,self.a_n_chan,n_last_samp/self.a_n_chan)';%Now reshape the unidimensional data matrix acourding to the number of channels
                %data=data(:,self.a_channs_cell);%Consider only the defined channels
                
%                 self.last_samp_off=self.curr_off+(n_last_samp*2);%Bruno
                
                self.a_is_last_segment=1;%Set the last segment flag
            end
            
            
            %data=self.combine_data(data);
            
            
            
            
            if self.a_step<=self.a_wsize
                self.a_data_map.Format={self.a_dtype [self.a_n_chan self.a_step_samp] 'dt'};
            else
                self.a_data_map.Format={self.a_dtype [self.a_n_chan nsamples] 'dt'};
            end
            
            data=((double(data)-self.a_elect_list(1).l_grd).*self.a_conv_factor);%Convert data to the nV, uV, etc range
            data=self.combine_data(data)';
            
            curr_sample=(self.a_curr_bytes_offset-self.a_first_data_bytes_offset)/self.a_n_bytes/self.a_n_chan;
            time=[curr_sample:curr_sample+size(data,2)-1].*(1/self.a_samp_freq);
            
            
            
            self.a_last_data_mat=data;
            
        end%End def_data_access method

        
        
        function data=get_next_window(self)
            %==============================================================
            %"get_next_segment" method
            %This method returns the next non-overlaping data
            %portion.
            % If the end of file was reachedk, in the next call a empty matrix
            %is returned.
            %
            %
            %EU FP7 Grant 211713 (EPILEPSIAE)
            %
            %César A. D. Teixeira
            %CISUC, FCTUC, University of Coimbra
            %September 2009
            %==============================================================
            %if ~self.a_is_last_segment%If the last segment flag was previously set then return an empty matrix
            
            
            if self.a_last_bytes_offset<self.a_abs_file_pt
                self.a_last_flagg=0;
                self.a_data_map.Offset=self.a_last_bytes_offset;%Point to the last sample of the last window
                %self.a_data_map.Format={self.a_dtype [self.a_n_chan self.a_step_samp] 'dt'};

                
                %data=(self.a_data_map.Data(1).dt)';
                %data=data(:,self.a_channs_cell);

                %disp(self.a_n_bytes_ahead)
                
                
                try%Try if it is possible to acquire the requested amount of data
                    data=(self.a_data_map.Data(1).dt)';
                    %data=data(:,self.a_channs_cell);
                    self.a_is_last_segment=0;
                    self.a_is_first_segment=0;
                
                    %((double(data(2,1:10))-self.a_elect_list(1).l_grd).*self.a_conv_factor)
                catch%If an error is raised, then redifine the memmapfile object parameters to adjust to the remaining data in file.
                    self.a_data_map.Format=self.a_dtype;%Consider independent data samples, not arranged acourding to the number of channels
                    self.a_data_map.Repeat=Inf;%Repeat the data format to the end of file
                    last_data=self.a_data_map.Data;%Get data
                    n_last_samp=size(last_data,1);%Compute how many individual samples were readed
                    
                                        
                    
                    data=reshape(last_data,self.a_n_chan,n_last_samp/self.a_n_chan)';%Now reshape the unidimensional data matrix acourding to the number of channels
                     %((double(data(2,1:10))-self.a_elect_list(1).l_grd).*self.a_conv_factor)              
                    %self.a_is_last_segment=1;%Set the last segment flag
                end
                
                
                
                
                self.a_last_bytes_offset=self.a_last_bytes_offset+self.a_n_bytes_ahead;
                self.a_curr_bytes_offset=self.a_curr_bytes_offset+self.a_n_bytes_ahead;
                
                self.a_last_samp_offset=self.a_last_samp_offset+self.a_step_samp;
                self.a_curr_samp_offset=self.a_curr_samp_offset+self.a_step_samp;
                
                 if (self.a_curr_bytes_offset+(3*self.a_n_bytes_ahead))>self.a_abs_file_pt
                    self.a_is_last_segment=1;
                    %disp(self.a_is_last_segment)
                 end
%                 
%                 
%                 self.a_is_first_segment=0;
                

                data=((double(data)-self.a_elect_list(1).l_grd).*self.a_conv_factor);%Convert data to the nV, uV, etc range
                
                data=self.combine_data(data)';
                
                if self.a_step<self.a_wsize%If the step is smaller than the window
                        
                        %Update the original data by appending the new
                        %retrieved data portion with the existing one,
                        %obviously eliminating the data out of the time
                        %window
                        
                                        
                        data=[self.a_last_data_mat(:,self.a_step_samp+1:end),data];
                        
                        
                        
                        
                 end
                
                self.a_last_data_mat=data;
                                          
                
            else
                if ~self.a_last_flagg
                    data=self.a_last_data_mat;
                    self.a_last_flagg=1;
                else
                    data=[];
                end
                
            end


        end%End get_next_segment method
        
        
        
        function data=get_prev_window(self)
            %==============================================================
            %"get_prev_segment" method
            %This method returns the previous non-overlaping data
            %portion.
            % If the beggin of data was reached, in the next call a empty matrix
            %is returned.
            %
            %
            %EU FP7 Grant 211713 (EPILEPSIAE)
            %
            %César A. D. Teixeira
            %CISUC, FCTUC, University of Coimbra
            %September 2009
            %==============================================================
            if ~self.a_is_first_segment%If the last segment flag was previously set then return an empty matrix
                
                
                
                
                
                self.a_curr_bytes_offset=self.a_curr_bytes_offset-self.a_n_bytes_ahead;
                self.a_last_bytes_offset=self.a_last_bytes_offset-self.a_n_bytes_ahead;
                
                self.a_curr_samp_offset=self.a_curr_samp_offset-self.a_step_samp;
                self.a_last_samp_offset=self.a_last_samp_offset-self.a_step_samp;
                
                

                self.a_data_map.Offset=self.a_curr_bytes_offset;%Point to the first sample to acquire
                %self.a_data_map.Format={self.a_dtype [self.a_n_chan self.a_step_samp] 'dt'};

                if (self.a_curr_bytes_offset-(self.a_n_bytes_ahead))<self.a_first_data_bytes_offset
                    self.a_is_first_segment=1;

                end
                data=(self.a_data_map.Data(1).dt)';
                %data=data(:,self.a_channs_cell);

                self.a_is_last_segment=0;

                data=((double(data)-self.a_elect_list(1).l_grd).*self.a_conv_factor);%Convert data to the nV, uV, etc range
                data=self.combine_data(data)';
                
                
                if self.a_step<=self.a_wsize%If the step is smaller than the window
                        
                        %Update the original data by appending the new
                        %retrieved data portion with the existing one,
                        %obviously eliminating the data out of the time
                        %window
                        
                      
                                               
                        data=[data,self.a_last_data_mat(:,1:end-self.a_step_samp)];
                        
                        self.a_last_data_mat=data;
                       
                    end
                
                
                
                
            else
                data=[];%No data to return, i.e. end-of-file
            end


        end%End get_next_segment method
        
        function self = redefine_data_access(self,wsize,step, offset)
            %==============================================================
            %"redefine_trc_par" method
            %This method redefines the parameters of the reading function
            % 
            %
            %
            %EU FP7 Grant 211713 (EPILEPSIAE)
            %
            %César A. D. Teixeira & Bruno Leitao
            %CISUC, FCTUC, University of Coimbra
            %October 2009
            %==============================================================
            
            self.a_step=step;
            self.a_wsize=wsize;
            nsamples=floor(self.a_wsize*self.a_samp_freq);%Compute the number of samples to be acquired
            
            
            
            self.a_n_bytes_win=floor(self.a_wsize*self.a_samp_freq)*self.a_n_chan*self.a_n_bytes;
            
            
            self.a_step_samp=floor(self.a_step*self.a_samp_freq);%Compute the step in samples
            
            self.a_n_bytes_ahead=self.a_step_samp*self.a_n_chan*self.a_n_bytes;%Compute the step in bytes
            
            
            self.a_data_map.Format={self.a_dtype [self.a_n_chan nsamples] 'dt'};%Bruno
            
            
               
                tmp=self.a_curr_bytes_offset;
                self.a_curr_bytes_offset = self.a_first_data_bytes_offset+(offset*self.a_n_bytes * self.a_n_chan);
                                          
                self.a_last_bytes_offset=self.a_curr_bytes_offset+self.a_n_bytes_win;%Define the offset of the last data sample in the first window;
            
                if self.a_last_bytes_offset<=self.a_abs_file_pt 
                    self.a_curr_bytes_offset =tmp;
                    self.a_last_bytes_offset=tmp+self.a_n_bytes_win;
               
                
                
                self.a_curr_samp_offset = offset;
            
                self.a_last_samp_offset=self.a_curr_samp_offset+nsamples;%Define the offset of the last data sample in the first window;
                
                self.a_data_map.Offset=self.a_curr_bytes_offset;
                
                self.a_is_last_segment=0;
                
                
               
                
            
            end

        end%End "redefine_data_access" method
        
        
        
%         function get_montages(self)
%             %==============================================================
%             %"get_montages" method
%             %This method acquires the information related with the
%             %electrode order
%             %This method returned and array, which contains the coded list
%             %of the electrodes stored in the file
%             %This information is stored in the "order" attribute
%             %
%             %
%             %EU FP7 Grant 211713 (EPILEPSIAE)
%             %
%             %César A. D. Teixeira
%             %CISUC, FCTUC, University of Coimbra
%             %September 2009
%             %==============================================================
%             fseek(self.a_fd,self.a_mont_ar_start,'bof');%Go to the code area in the file
%             lines=fread(self.a_fd,1,'uint16')
%             sectors=fread(self.a_fd,1,'uint16')
%             base_time=fread(self.a_fd,1,'uint16')
%             notch=fread(self.a_fd,1,'uint16')
%             color=fread(self.a_fd,32,'uchar')
%         end%End get_montages method
        
        
    end%End Methods
end%End Class
