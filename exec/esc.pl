/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Enumerate 3+3 trial outcomes for exact matrix calculations

   This file contains the main logic. It runs with every ISO Prolog system
   that provides a few conforming extensions such as definite clause
   grammars (DCGs) and finite domain constraints.

   Tested with Scryer Prolog and SICStus Prolog.
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

% Set the system up so that strings in double quotes are lists of characters.
% In Scryer Prolog and other recent systems, this is already the default.
:- set_prolog_flag(double_quotes, chars).

% Prefix op * for 'generalizing away' goals (https://www.metalevel.at/prolog/debugging)
:- op(920, fy, *). *_.  % *Goal always succeeds

% DCG esc(D, Lo..Hi) describes a list of 3+3 cohorts FOLLOWING a dose D in Lo..Hi.
% One can read esc(D, Lo..Hi) as the DECISION to escalate *from* D to min(D+1,Hi).
% For example, esc(0, 0..5) *initiates* a trial that has 5 prespecified doses,
% and enrolls the first cohort at D=1.

% - - - - - begin inset for paper - - - - -

tox(T) :- T in 0..3,
	  indomain(T).

% Mnemonic: * is ^ that 'splatted' on dose ceiling.
esc(Hi,Lo..Hi) --> [Hi * T], { tox(T) },
		   (  {T #=< 1}, [mtd_notfound(Hi)]
		   ;  {T #>= 2}, des(Hi, Lo)
		   ).
esc(D, Lo..Hi) --> { D #< Hi, D1 #= D + 1 },
		   [D1 ^ T], { tox(T) },
		   (  {T #= 0}, esc(D1, Lo..Hi)
		   ;  {T #= 1}, sta(D1, Lo..Hi)
		   ;  {T #> 1}, des(D1, Lo)
		   ).

sta(D,  _..D ) --> [D - 0], [mtd_notfound(D)].
sta(D, Lo..Hi) --> { D #< Hi, D in Lo..Hi },
		   [D - 0],
		   esc(D, D..Hi).
sta(D, Lo.._ ) --> [D - T], { tox(T), T #> 0 },
		   des(D, Lo).

% As a mirror image of esc//2, des(D, Lo) moves
% downward FROM D, to max(D-1,Lo).
% NB: De-escalation to D-1 clamps Hi #= D - 1.
des(D, Lo) --> { D_1 #= D - 1 },
	       (  {D_1 #= Lo}, [declare_mtd(Lo)]
	       ;  {D_1 #> Lo}, [D_1 : T], {tox(T)},
		  (  {T #=< 1}, [declare_mtd(D_1)]
		  ;  {T #>= 2}, des(D_1, Lo)
		  )
	       ).

% - - - - - end inset for paper - - - - -

% n_trials(+Drange, DN)
n_trials(Drange, DN) :-
    Dmax in Drange, indomain(Dmax),
    findall(Tr, phrase(esc(0, 0..Dmax), Tr), Trials),
    length(Trials, N),
    DN = (Dmax, N).

%% Here are D & J columns of Table 1 in the 'WWTT' paper:

%?- n_trials(1..10, DN).
%@    DN = (1,10)
%@ ;  DN = (2,46)
%@ ;  DN = (3,154)
%@ ;  DN = (4,442)
%@ ;  DN = (5,1162)
%@ ;  DN = (6,2890)
%@ ;  DN = (7,6922)
%@ ;  DN = (8,16138)
%@ ;  DN = (9,36874)
%@ ;  DN = (10,82954)
%@ ;  false.

%% See all 10 paths in a 1-dose 3+3 trial:

%?- phrase(esc(0, 0..1), Tr).
%@    Tr = [1^0,1*0,mtd_notfound(1)]
%@ ;  Tr = [1^0,1*1,mtd_notfound(1)]
%@ ;  Tr = [1^0,1*2,declare_mtd(0)]
%@ ;  Tr = [1^0,1*3,declare_mtd(0)]
%@ ;  Tr = [1^1,1-0,mtd_notfound(1)]
%@ ;  Tr = [1^1,1-1,declare_mtd(0)]
%@ ;  Tr = [1^1,1-2,declare_mtd(0)]
%@ ;  Tr = [1^1,1-3,declare_mtd(0)]
%@ ;  Tr = [1^2,declare_mtd(0)]
%@ ;  Tr = [1^3,declare_mtd(0)]
%@ ;  false.

%% See all 46 paths in a 2-dose 3+3 trial:

%?- phrase(esc(0, 0..2), Tr).
%@    Tr = [1^0,2^0,2*0,mtd_notfound(2)]
%@ ;  Tr = [1^0,2^0,2*1,mtd_notfound(2)]
%@ ;  Tr = [1^0,2^0,2*2,1:0,declare_mtd(1)]
%@ ;  Tr = [1^0,2^0,2*2,1:1,declare_mtd(1)]
%@ ;  Tr = [1^0,2^0,2*2,1:2,declare_mtd(0)]
%@ ;  Tr = [1^0,2^0,2*2,1:3,declare_mtd(0)]
%@ ;  Tr = [1^0,2^0,2*3,1:0,declare_mtd(1)]
%@ ;  Tr = [1^0,2^0,2*3,1:1,declare_mtd(1)]
%@ ;  Tr = [1^0,2^0,2*3,1:2,declare_mtd(0)]
%@ ;  Tr = [1^0,2^0,2*3,1:3,declare_mtd(0)]
%@ ;  Tr = [1^0,2^1,2-0,mtd_notfound(2)]
%@ ;  Tr = [1^0,2^1,2-1,1:0,declare_mtd(1)]
%@ ;  Tr = [1^0,2^1,2-1,1:1,declare_mtd(1)]
%@ ;  Tr = [1^0,2^1,2-1,1:2,declare_mtd(0)]
%@ ;  Tr = [1^0,2^1,2-1,1:3,declare_mtd(0)]
%@ ;  Tr = [1^0,2^1,2-2,1:0,declare_mtd(1)]
%@ ;  Tr = [1^0,2^1,2-2,1:1,declare_mtd(1)]
%@ ;  Tr = [1^0,2^1,2-2,1:2,declare_mtd(0)]
%@ ;  Tr = [1^0,2^1,2-2,1:3,declare_mtd(0)]
%@ ;  Tr = [1^0,2^1,2-3,1:0,declare_mtd(1)]
%@ ;  Tr = [1^0,2^1,2-3,1:1,declare_mtd(1)]
%@ ;  Tr = [1^0,2^1,2-3,1:2,declare_mtd(0)]
%@ ;  Tr = [1^0,2^1,2-3,1:3,declare_mtd(0)]
%@ ;  Tr = [1^0,2^2,1:0,declare_mtd(1)]
%@ ;  Tr = [1^0,2^2,1:1,declare_mtd(1)]
%@ ;  Tr = [1^0,2^2,1:2,declare_mtd(0)]
%@ ;  Tr = [1^0,2^2,1:3,declare_mtd(0)]
%@ ;  Tr = [1^0,2^3,1:0,declare_mtd(1)]
%@ ;  Tr = [1^0,2^3,1:1,declare_mtd(1)]
%@ ;  Tr = [1^0,2^3,1:2,declare_mtd(0)]
%@ ;  Tr = [1^0,2^3,1:3,declare_mtd(0)]
%@ ;  Tr = [1^1,1-0,2^0,2*0,mtd_notfound(2)]
%@ ;  Tr = [1^1,1-0,2^0,2*1,mtd_notfound(2)]
%@ ;  Tr = [1^1,1-0,2^0,2*2,declare_mtd(1)]
%@ ;  Tr = [1^1,1-0,2^0,2*3,declare_mtd(1)]
%@ ;  Tr = [1^1,1-0,2^1,2-0,mtd_notfound(2)]
%@ ;  Tr = [1^1,1-0,2^1,2-1,declare_mtd(1)]
%@ ;  Tr = [1^1,1-0,2^1,2-2,declare_mtd(1)]
%@ ;  Tr = [1^1,1-0,2^1,2-3,declare_mtd(1)]
%@ ;  Tr = [1^1,1-0,2^2,declare_mtd(1)]
%@ ;  Tr = [1^1,1-0,2^3,declare_mtd(1)]
%@ ;  Tr = [1^1,1-1,declare_mtd(0)]
%@ ;  Tr = [1^1,1-2,declare_mtd(0)]
%@ ;  Tr = [1^1,1-3,declare_mtd(0)]
%@ ;  Tr = [1^2,declare_mtd(0)]
%@ ;  Tr = [1^3,declare_mtd(0)]
%@ ;  false.

%% Transform dose-escalation path lists to the arrays T(c,d,j).

path_matrix(P, D, M) :-
    phrase(pm_(M), P),
    M = C1-C2,
    length(C1,D),
    length(C2,D),
    maplist(ground_or_na, C1),
    maplist(ground_or_na, C2).

pm_(C1-C2) -->
    [D^T], % ascribe to D's 1st cohort
    { nth1(D, C1, T) },
    pm_(C1-C2).

pm_(C1-C2) -->
    (	[D-T]
    ;	[D*T]
    ;	[D:T]
    ), % ascribe to D's 2nd cohort
    { nth1(D, C2, T) },
    pm_(C1-C2).

pm_(_-_) -->
    (	[declare_mtd(_)]
    ;	[mtd_notfound(_)]
    ).

ground_or_na(Term) :- ground(Term).
ground_or_na('NA').

columns_format(1) --> "~w~n".
columns_format(N) --> "~w\t",
		      { N #> 1,
			N1 #= N - 1 },
		      columns_format(N1).

%% Write out tab-delimited input files T<D>.tab
write_T(D) :-
    findall(P, phrase(esc(0, 0..D), P), Paths),
    length(Paths, Len),
    format("D = ~d ~t~8+J = ~d~n", [D, Len]), % feedback to console
    length(Ds, Len), maplist(=(D), Ds),
    maplist(path_matrix, Paths, Ds, Ms) ->
	write_esc_array(D, Ms).

write_esc_array(D, Ms) :-
    phrase(format_("T~d.tab", [D]), Filename), atom_chars(File, Filename),
    open(File, write, OS),
    phrase(columns_format(D), Format),
    phrase(write_esc_array_(OS, Format), Ms) -> close(OS).

write_esc_array_(OS, Format) -->
    [C1-C2],
    { format(OS, Format, C1),
      format(OS, Format, C2) },
    write_esc_array_(OS, Format).
write_esc_array_(_, _) --> [].

:- initialization((maplist(write_T, [2,3,4,5,6,7,8]) -> halt)).
