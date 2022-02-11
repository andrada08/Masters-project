%% Plot behavior (for checking ongoing)

% animals = {'AP110','AP111','AP112'};
% animals = {'AP098','AP099','AP112'};
animals = {'AP108'};

protocol = 'AP_stimWheel';
flexible_name = true;
bhv = struct;

for curr_animal = 1:length(animals)
    
    preload_vars = who;
    
    animal = animals{curr_animal};
    experiments = AP_find_experiments(animal,protocol,flexible_name);
    
    if isempty(experiments)
        disp(['No behavior data: ' animal]);
        continue
    end
    
    disp([animal ', day:'])
    for curr_day = 1:length(experiments)
        
        day = experiments(curr_day).day;
        experiment_num = experiments(curr_day).experiment;
        
        % If multiple experiments, only use the last one
        % (usually multiple happens if mess ups and final one is good)
        for curr_experiment = length(experiment_num)
            
            experiment = experiment_num(curr_experiment);
            
            % Load (behavior/timeline only)
            load_parts.imaging = false;
            AP_load_experiment;
            
            % Get protocol name
            [~,curr_protocol] = fileparts(block.expDef);
            
            % Time of session (in minutes)
            session_duration = block.duration/60;
            
            % Trial counts (use only full trials with a response)
            total_water = sum(block.outputs.rewardValues);
            
            % Hit fraction
            frac_hit = nanmean(trial_outcome == 1);
            
            % (probability of any movement around stim)
            stim_surround_t_centers = -10:0.1:10;
            stim_surround_times = stimOn_times + stim_surround_t_centers;
            stim_surround_move = interp1(Timeline.rawDAQTimestamps, ...
                +wheel_move,stim_surround_times,'previous');
            
            % Resetting quiescence period for each trial
            if isfield(block.events,'trialQuiescenceValues')
                quiescence_t = block.events.trialQuiescenceValues(1:n_trials);
            else
                quiescence_t = 0.5*ones(size(block.events.stimOnTimes));
            end
            
            % ITI time (including quiescence resets)
            iti_t = block.events.stimOnTimes(2:n_trials) - ...
                block.events.stimOffTimes(1:n_trials-1);
            
            % Wheel movements/biases
            left_wheel_velocity = abs(wheel_velocity.*(wheel_velocity < 0));
            right_wheel_velocity = abs(wheel_velocity.*(wheel_velocity > 0));
            wheel_bias = (nansum(right_wheel_velocity)-nansum(left_wheel_velocity))/ ...
                (nansum(right_wheel_velocity)+nansum(left_wheel_velocity));
            
            % Store in behavior structure
            % (general timings)
            bhv(curr_animal).animal = animal;
            bhv(curr_animal).day{curr_day,1} = day;
            bhv(curr_animal).protocol{curr_day,1} = curr_protocol;
            bhv(curr_animal).session_duration(curr_day,1) = session_duration;
            bhv(curr_animal).n_trials(curr_day,1) = n_trials;
            bhv(curr_animal).frac_hit(curr_day,1) = frac_hit;
            bhv(curr_animal).total_water(curr_day,1) = total_water;
            bhv(curr_animal).wheel_velocity(curr_day,1) = nansum(abs(wheel_velocity));
            bhv(curr_animal).stim_move_t{curr_day,1} = stim_to_move;
            bhv(curr_animal).stim_feedback_t{curr_day,1} = stim_to_feedback;
            bhv(curr_animal).quiescence_t{curr_day,1} = quiescence_t';
            bhv(curr_animal).iti_t{curr_day,1} = iti_t';
            bhv(curr_animal).wheel_bias(curr_day,1) = wheel_bias;
            
            % (stim-aligned movement)
            bhv(curr_animal).stim_surround_t = stim_surround_t_centers;
            bhv(curr_animal).stim_surround_wheel{curr_day,1} = stim_surround_move;
            
            AP_print_progress_fraction(curr_day,length(experiments));
        end
        
    end
    
    % Plot summary
    day_num = cellfun(@(x) datenum(x),{experiments.day}');
    day_labels = cellfun(@(day,protocol) [day(6:end)], ...
        {experiments.day}',bhv(curr_animal).protocol,'uni',false);
    
    [unique_protocols,~,protocol_idx] = unique(bhv(curr_animal).protocol);
    protocol_col = hsv(length(unique_protocols));
    
    figure('Name',animal)
    
    % Trials and water
    subplot(2,3,1); hold on;
    yyaxis left
    % plot(day_num,bhv(curr_animal).n_trials./bhv(curr_animal).session_duration,'linewidth',2);
    % ylabel('Trials/min');
    plot(day_num,bhv(curr_animal).n_trials./bhv(curr_animal).session_duration,'linewidth',2);
    ylabel('Trials/min');
    yyaxis right
    plot(day_num,bhv(curr_animal).total_water,'linewidth',2);
    protocol_plot = gscatter(day_num,zeros(size(day_num)),[bhv(curr_animal).protocol]);
    legend off;
    ylabel('Total water');
    xlabel('Day');
    set(gca,'XTick',day_num);
    set(gca,'XTickLabel',day_labels);
    set(gca,'XTickLabelRotation',90);
    
    % Wheel movement and bias
    subplot(2,3,2);
    yyaxis left
    % plot(day_num,bhv(curr_animal).wheel_velocity./bhv(curr_animal).session_duration,'linewidth',2);
    % ylabel('Wheel movement / min');
    plot(day_num,bhv(curr_animal).wheel_velocity,'linewidth',2);
    ylabel('Wheel movement');
    yyaxis right
    plot(day_num,bhv(curr_animal).wheel_bias,'linewidth',2);
    ylim([-1,1]);
    line(xlim,[0,0]);
    ylabel('Wheel bias');
    xlabel('Day');
    set(gca,'XTick',day_num);
    set(gca,'XTickLabel',day_labels);
    set(gca,'XTickLabelRotation',90);
    
    % Stim to move/reward time
    subplot(2,3,3); hold on;
    plot(day_num,cellfun(@nanmedian,bhv(curr_animal).stim_move_t),'k','linewidth',2);
    plot(day_num,cellfun(@nanmedian,bhv(curr_animal).stim_feedback_t),'b','linewidth',2);
    ylabel('Stim to (move/reward) time (s)');
    xlabel('Session')
    set(gca,'YScale','log')
    set(gca,'XTick',day_num);
    set(gca,'XTickLabel',day_labels);
    set(gca,'XTickLabelRotation',90);
    
    % Stim-aligned move fraction
    stim_surround_wheel_cat = ...
        cell2mat(cellfun(@(x) nanmean(x,1),bhv(curr_animal).stim_surround_wheel,'uni',false));
    subplot(2,4,5);
    imagesc(bhv(curr_animal).stim_surround_t,[],stim_surround_wheel_cat);
    colormap(gca,brewermap([],'Greys'));
    caxis([0,1]);
    xlabel('Time from stim (s)');
    set(gca,'YTick',1:length(day_num));
    set(gca,'YTickLabel',day_labels);
    
    subplot(2,4,6); hold on;
    set(gca,'ColorOrder',copper(size(stim_surround_wheel_cat,1)));
    plot(bhv(curr_animal).stim_surround_t,stim_surround_wheel_cat');
    ylim([0,1]);
    ylabel('Fraction move');
    xlabel('Time from stim (s)');
    
    % Reaction time histogram
    rxn_bins = [0:0.05:1];
    rxn_bin_centers = rxn_bins(1:end-1) + diff(rxn_bins)./2;
    rxn_hist = cell2mat(cellfun(@(x) ...
        histcounts(x,rxn_bins,'normalization','probability'), ...
        bhv(curr_animal).stim_move_t,'uni',false));
    
    subplot(2,4,7);
    imagesc(rxn_bin_centers,[],rxn_hist);
    colormap(gca,brewermap([],'Greys'));
    caxis([0,0.5]);
    set(gca,'YTick',1:length(day_num));
    set(gca,'YTickLabel',day_labels);
    xlabel('Reaction time');
    
    subplot(2,4,8); hold on;
    set(gca,'ColorOrder',copper(size(stim_surround_wheel_cat,1)));
    plot(rxn_bin_centers,rxn_hist');
    xlabel('Reaction time');
    ylabel('Fraction');
    
    drawnow;
    
    clearvars('-except',preload_vars{:});
    
end