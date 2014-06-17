vol = [];
vol.c = [3.3 10.0 0.042 3.3];
vol.r = [0.79 0.81 0.85 0.88];
elec = dlmread('/home/user/code/pelops/karakum/build/elec.txt',' ');
eeg_leadfield4([0 0 0], elec, vol);
