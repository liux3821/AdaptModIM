function poppop(delta_ort,sndInfo)

if delta_ort > 0
    
    if delta_ort > sndInfo.numSnds
        
        s = PsychPortAudio('GetStatus', sndInfo.pahandle);
        if s.State == 0 % s.Active == 0
            % Schedule finished, engine stopped. Before adding new
            % slots we first must delete the old ones, ie., reset the
            % schedule:
            PsychPortAudio('UseSchedule', sndInfo.pahandle, 2);
        end
        
        PsychPortAudio('AddToSchedule', sndInfo.pahandle, sndInfo.blipBuffer(end));
        
        if s.State == 0 % s.Active == 0
            PsychPortAudio('Start', sndInfo.pahandle, 1, 0);
        end
        
    else
        
        s = PsychPortAudio('GetStatus', sndInfo.pahandle);
        if s.State == 0 % s.Active == 0
            % Schedule finished, engine stopped. Before adding new
            % slots we first must delete the old ones, ie., reset the
            % schedule:
            PsychPortAudio('UseSchedule', sndInfo.pahandle, 2);
        end
        
        PsychPortAudio('AddToSchedule', sndInfo.pahandle, sndInfo.blipBuffer(delta_ort));
        
        if s.State == 0 % s.Active == 0
            PsychPortAudio('Start', sndInfo.pahandle, 1, 0);
        end
        
    end
end