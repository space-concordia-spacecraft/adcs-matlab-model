try
start_time = datestr([year month day hour minute second],'dd mmm yyyy HH:MM:SS.FFF');

fileID = fopen('AttitudeDataSCODIN.a','w');

fprintf(fileID,['stk.v.11.0\n'...
'BEGIN Attitude\n' ...
'NumberOfAttitudePoints		' num2str(size(q_EB,1)) '\n'...
    'BlockingFactor          20\n'...
    'InterpolationOrder      1\n'...
'CentralBody			Earth\n'...
'ScenarioEpoch			' start_time '\n'...
'CoordinateAxes		Fixed\n'...
'AttitudeTimeQuaternions\n']);

fprintf(fileID,'%4.2f %1.7f %1.7f %1.7f %1.7f\n',transpose(q_EB));

fprintf(fileID, 'END Attitude\n');
fclose(fileID);

Attitudefile = dir('AttitudeDataSCODIN.a');

durationoffile = seconds(q_EB(end,1));
msgbox({'STK Attitude data exported to file!';'AttitudeDataSCODIN.a';...
    ['Start time STK : ' start_time] ;...
    ['Stop time STK : ' datestr(addtodate(datenum(start_time,'dd mmm yyyy HH:MM:SS.FFF'),q_EB(end,1),'second'),'dd mmm yyyy HH:MM:SS.FFF')] ;...
    ['Duration : ' datestr(durationoffile,'dd') ' day  ' datestr(durationoffile,'hh') 'h ' datestr(durationoffile,'MM') 'm ' datestr(durationoffile,'ss') 's'] ;...
    ['Data : ' num2str(size(q_EB,1)) ' Attitude Points'];...
    ['Size : ' num2str(Attitudefile.bytes/1024, '%.0f') ' KB']},'SC-ODIN STK','help');

catch
    msgbox('ERROR : Attitude data not exported to file','SC-ODIN STK','help');
end
    