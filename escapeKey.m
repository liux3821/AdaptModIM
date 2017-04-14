function escapeKey(keyIsDown,keyCode)
if keyIsDown == 1
    if keyCode(KbName('ESCAPE')) %
        Priority(0);
        Screen('Close');
        RestoreCluts;
        Screen('CloseAll');
        return
    end
end
end