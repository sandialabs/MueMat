function mue()
  d = date;
  if (all(str2num(int2str(d(1:6))) == [49 57 45 65 112 114]))
    fprintf('              __ \n             (..)  "MueMat"\n      /------/||\\\n     / |     ||\n    *  ||----||\n       ~~    ~~\n');
  else
    fprintf('             (__)\n             (oo)  "MueMat"\n      /-------\\/ \n     / |     ||\n    *  ||----||\n       ~~    ~~\n');
  end
end
