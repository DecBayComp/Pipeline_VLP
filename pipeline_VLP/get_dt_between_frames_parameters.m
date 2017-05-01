function parameters = get_dt_between_frames_parameters(movie_per_frame,parameters)

%% time between frames
ddtt = [];
parameters.dt_theo                       = movie_per_frame(2).t - movie_per_frame(1).t;
for i = 1 : length ( movie_per_frame)-1
    ddtt = [ddtt;movie_per_frame(i+1).t - movie_per_frame(i).t ];
end
parameters.dt_theo = min(ddtt);

end