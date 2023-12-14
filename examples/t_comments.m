%
%
interpretation( 3, [foo], [ % more comments
% foo
  function(*(_,_), [
    2, 1, 0,
    0, 2, 1, % baz
    1, 0, 2])]). 
% bar
interpretation( 4, [foo baz], [ % more comments
% foo
  function(*(_,_), [
    2, 1, 0, 3, 0, 0, 1, 3, % baz
    1, 0, 2, 3, 3, 3, 3, 3  % ******* *****
    ]) %%%%%%%%%%%%%%%%%%%
  ]). 
