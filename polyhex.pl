
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                    %%
%% http://www.research.att.com/~njas/sequences/polyhexp.txt           %%
%%                                                                    %%
%% A set of Prolog-definitions that illustrate how the first terms    %%
%% of A0xxxxx are produced.                                           %%
%%                                                                    %%
%% Written by Antti Karttunen, 2004, http://www.iki.fi/kartturi/      %%
%%                                                                    %%
%% Edited:                                                            %%
%%   2004 10 09.                                                      %%
%%   2006 02 27 (Added the Gilbert, Dorfman, Gaspard reference)       %%
%%   2012 04 28-29. (Added maxrotation, etc. Improved selfavoiding.)  %%
%%                                                                    %%
%%                                                                    %%
%% This works with GNU prolog:                                        %%
%% http://www.gnu.org/software/gprolog/gprolog.html                   %%
%%                                                                    %%
%% Load as:                                                           %%
%% consult('/karttu/prolog/polyhexp.txt').                            %%
%%                                                                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Encoding like this is used for example in:
%% T. Gilbert, J. R. Dorfman, and P. Gaspard,
%% "Fractal dimensions of the hydrodynamic modes of diffusion"

%% http://arxiv.org/pdf/nlin.CD/0007008
%% or: http://arxiv.org/PS_cache/nlin/pdf/0007/0007008.pdf

%% On page 4:
%% 
%% Alternatively, a convenient representation of phase space is to consider binary sequences
%% of zeros and ones labeling respectively hops to the left and right. A trajectory of the random
%% walker is then given by an integer coordinate representing, say, the position at time zero
%% and an arbitrary binary sequence, coding all its (past and future) displacements.
%% 
%% and later:
%% 
%% Similar functions can be observed for the 2-adic and 5-adic totally symmetric random walks
%% for k_c = Pi/3. In the 2-adic case, the set is a regular triangular lattice; in the 5-adic case,
%% a regular hexagonal lattice.
%% 

%% This encoding was also independently discovered by Y. Atias.


%% For rotating list left (or right).

rot([],[]).

rot([X|Xs],Y) :-
  append(Xs,[X],Y).

add(X,Y,Z) :-
  length(LX,X),
  length(LY,Y),
  append(LX,LY,LZ),
  length(LZ,Z),
  !.

equlengths([],[]).

equlengths([_|X],[_|Y]) :-
  equlengths(X,Y).

nzeros([],0) :-
  !.

nzeros([0|X],N) :-
  M is N-1,
  nzeros(X,M),
  !.


%% P = Phase: 0-5, *60 degrees. (e.g. 5 = 300 degrees).
%% D = Direction: 0 = turn 60 degrees right (CW). 1 = turn 60 deg. left (CCW).
%% NP = New Phase = P + 2*D - 1 mod 6.

newphase(P,D,NP) :-
  XNP is P + 2*D - 1,
  NP is XNP mod 6.

%% First arg = Phase, 0-5, *60 degrees. (e.g. 5 = 300 degrees).
%% X & Y = Old X & Y.
%% NX & NY = New X & Y.

newcoords(0,X,Y,NX,NY) :-
  NX is X+1,
  NY is Y.

newcoords(1,X,Y,NX,NY) :-
  NX is X,
  NY is Y+1.

newcoords(2,X,Y,NX,NY) :-
  NX is X-1,
  NY is Y+1.

newcoords(3,X,Y,NX,NY) :-
  NX is X-1,
  NY is Y.

newcoords(4,X,Y,NX,NY) :-
  NX is X,
  NY is Y-1.

newcoords(5,X,Y,NX,NY) :-
  NX is X+1,
  NY is Y-1.


uniqelems([]) :-
  !.

uniqelems([CP|CPS]) :-
  memberchk(CP,CPS) -> fail ; uniqelems(CPS).


discardduplicates([],[]) :-
  !.

% If X is found in the tail, then just discard it and check the rest.
discardduplicates([X|Xs],Ys) :-
  memberchk(X,Xs),
  !, % Beware of red cut!
  discardduplicates(Xs,Ys).

% If X was not found in the tail, then include it in Ys obtained by tailrecursing the rest:
discardduplicates([X|Xs],[X|Ys]) :-
  discardduplicates(Xs,Ys).


intlistlessthan([],[]) :-
  fail.

intlistlessthan([X|_],[Y|_]) :-
  X < Y,
  !.

intlistlessthan([X|Xs],[Y|Ys]) :-
  X > Y -> fail;
  intlistlessthan(Xs,Ys).


maxrotationAux(_,_,0) :-
  !.

maxrotationAux(X,Y,N) :-
  intlistlessthan(X,Y) -> fail;
  rot(Y,Z),
  M is N-1,
  maxrotationAux(X,Z,M).


maxrotation(X) :-
  length(X,L),
  L1 is L-1,
  rot(X,X1),
  maxrotationAux(X,X1,L1).



%% When we have taken all the turns,
%% we should have returned back to origo [0|0]:
%% I.e. self-avoiding and self-returning!
selfXavoiding([],_,_,_,[[0|0]|_]) :-
  !.


selfXavoiding([B|Bs],P,X,Y,[C|CPS]) :-
  memberchk(C,CPS) -> fail ; %% It is enough to check that the first coord is not among the rest.
  newphase(P,B,NP), %% Compute the new phase from old phase and given direction
  newcoords(NP,X,Y,NX,NY), %% Compute the new coordinates in the  lattice.
  selfXavoiding(Bs,NP,NX,NY,[[NX|NY]|[C|CPS]]). %% And loop.


%% Foolish:
%% selfXavoiding([B|Bs],P,X,Y,CPS) :-
%%   uniqelems(CPS),   %% All coordinate pairs unique so far?
%%   newphase(P,B,NP), %% Compute the new phase from old phase and given direction
%%   newcoords(NP,X,Y,NX,NY), %% Compute the new coordinates in the  lattice.
%%   selfXavoiding(Bs,NP,NX,NY,[[NX|NY]|CPS]). %% And loop.


%% We start from origo [0|0] with phase 0.
%% If Bs starts with 1 (turn left), then the next phase
%% is 1, and the coords will be [0|1].
%% If Bs starts with 0 (turn right),
%% then the next phase is 5, and the coords will be [1|-1].

selfavoiding(Bs) :-
  selfXavoiding(Bs,0,0,0,[[0|0]]).

selfreaching(Bs) :-
  selfavoiding(Bs) -> fail ;
  !.


%% If we can repeatedly apply these two Post-productions:
%%   A011110B -> A11B
%% and
%%   A01110B  -> A101B
%% until the string of six 1's (111111) results,
%% then the original binary string was also a valid
%% holeless polyhex. (bullshit!)

%% Note that we can ignore the productions
%% A00B -> A100001B, A010B -> A10001B, A0110B -> A1001B
%% because every valid non-monic polyhex must contain at least
%% two instances of one of the convex paths 0110, 01110 or 011110
%% somewhere on its edge (e.g. 1111011110 or 111011101110)
%% and although not necessarily in the beginning of the sequence,
%% at least one instance is somewhere as whole, not broken by wrap-over.
%% This is fortunate, because including those two
%% additional productions would introduce other problems.


ispolyhex(X,N) :-
  ithexrewrite(X,[1,1,1,1,1,1],L),
  length([0|L],N),
  !.


%% ithexrewrite(X,Y,L) is true if X can be rewritten to Y
%% with (length L) polyhex-rewriting rules:
ithexrewrite([1,1,1,1,1,1],[1,1,1,1,1,1],[]) :-
  !.

ithexrewrite(X,Y,[0|L]) :-
  hexrewrite(X,Z),
  ithexrewrite(Z,Y,L).

%% Note the execution order of the clauses.
%% We should not call hexrewrite with two
%% uninstantiated arguments!

%% polyhexnhexes(X,N) is true if X can be rewritten to [1,1,1,1,1,1]
%% with N-1 polyhex-rewriting rules.

polyhexnhexes(X,N) :-
  M is N-1,
  it2hexrewrite([1,1,1,1,1,1],X,M).


polyhexnhexes1selfavoiding(X,N) :-
  M is N-1,
  it2hexrewrite([1,1,1,1,1,1],X,M),
  selfavoiding(X).


polyhexnhexes1maxrot(X,N) :-
  M is N-1,
  it2hexrewrite([1,1,1,1,1,1],X,M),
  maxrotation(X).


polyhexnhexes1selfavoiding1maxrot(X,N) :-
  M is N-1,
  it2hexrewrite([1,1,1,1,1,1],X,M),
  maxrotation(X),
  selfavoiding(X).


polyhexnhexes1selfreaching1maxrot(X,N) :-
  M is N-1,
  it2hexrewrite([1,1,1,1,1,1],X,M),
  selfreaching(X),
  maxrotation(X).


%% By the size (number of hexes).
%% it2hexrewrite(X,Y,N) matches when X can be rewritten
%% (expanded) to Y with N-1 rewrites:
it2hexrewrite(X,X,0) :-
  !.

it2hexrewrite(X,Y,N) :-
  M is N-1,
  hexrewrite(Z,X),
  it2hexrewrite(Z,Y,M).

%% polyhexperim(X,N) is true if X can be rewritten to [1,1,1,1,1,1]
%% (i.e. [1,1,1,1,1,1] can be rewritten to X)
%% with its length contracting by N-6 edges.
%% (I.e., its original length being N).
polyhexperim(X,N) :-
  M is N-6,
  it3hexrewrite([1,1,1,1,1,1],X,M).


polyhexperim1selfavoiding(X,N) :-
  M is N-6,
  it3hexrewrite([1,1,1,1,1,1],X,M),
  selfavoiding(X).

polyhexperim1maxrot(X,N) :-
  M is N-6,
  it3hexrewrite([1,1,1,1,1,1],X,M),
  maxrotation(X).

polyhexperim1selfavoiding1maxrot(X,N) :-
  M is N-6,
  it3hexrewrite([1,1,1,1,1,1],X,M),
  maxrotation(X),
  selfavoiding(X).



%% By the perimeter length:
%% it3hexrewrite(X,Y,N) matches when X can be rewritten
%% to Y and X is N edges shorter than Y.
%% I.e. the third argument is the limit how many
%% edges we can still add to the perimeter.
it3hexrewrite(X,X,0) :-
  !.

it3hexrewrite(X,Y,N) :-
  length(X,A),
  hexrewrite(Z,X),
  length(Z,B),
  D is B-A,
  D =< N,
  M is N-D,
  it3hexrewrite(Z,Y,M).


apolhex(N,Y) :-
  findall(L,polyhexnhexes(L,N),RESULTS),
  length(RESULTS,Y).

apolhexper(N,Y) :-
  M is N*2,
  findall(L,polyhexperim(L,M),RESULTS),
  length(RESULTS,Y).


apolhex1sa(N,Y) :-
  findall(L,polyhexnhexes1selfavoiding(L,N),RESULTS),
  length(RESULTS,Y).

apolhexper1sa(N,Y) :-
  M is N*2,
  findall(L,polyhexperim1selfavoiding(L,M),RESULTS),
  length(RESULTS,Y).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Begins as 1,1,3,10,33,146,618,2803,12824,
%% Compare to A006535 Number of one-sided hexagonal polyominoes with n cells.
%% and A038144 Number of planar n-hexes, or polyhexes (in the sense of A000228,
%% so rotations and reflections count as the same shape) with at least one hole.
%S A006535 1,1,3,10,33,147,620,2821,12942,60639,286190,1364621,6545430,31586358,
%S A038144 0,0,0, 0, 0,  1,  2,  13,   67,  404,  2323,  13517,  76570,  429320,2373965,13004323,


apolhex1sauptorot(N,Y) :-
  findall(L,polyhexnhexes1selfavoiding1maxrot(L,N),RESULTS),
  discardduplicates(RESULTS,RESULTS2),
  length(RESULTS2,Y).


apolhex1sruptorot(N,Y) :- %% Self-Reaching ones.
  findall(L,polyhexnhexes1selfreaching1maxrot(L,N),RESULTS),
  discardduplicates(RESULTS,RESULTS2),
  length(RESULTS2,Y).


%% Without self-avoiding condition, thus allowing one-sided polyhexes whose ends
%% touch with each other. (First such case occurs with size N=6).
%% Begins as 1,1,3,10,33,147,632,2935,13866,
%% case n=7: 632 (= 618 + 14) at least seems fine to me.

apolhex1uptorot(N,Y) :-
  findall(L,polyhexnhexes1maxrot(L,N),RESULTS),
  discardduplicates(RESULTS,RESULTS2),
  length(RESULTS2,Y).


%% Begins as 0,0,1,0,1,1,3,3,15,22,76,174,536,
%% Compare to:
%S A057779   0,0,1,0,1,1,3,2,12,14,50,98,
%N A057779 Hexagonal polyominoes (or polyhexes, A000228) with perimeter 2n.
%% and:
%N A130623 Number of polyhexes with perimeter at most 2n.


apolhexper1sauptorot(N,Y) :-
  M is N*2,
  findall(L,polyhexperim1selfavoiding1maxrot(L,M),RESULTS),
  discardduplicates(RESULTS,RESULTS2),
  length(RESULTS2,Y).


%% Begins as 0,0,1,0,1,1,3,3,15,22,76,174,537,

% Without self-avoiding condition, thus allowing one-sided polyhexes whose ends
% touch with each other. (First such case occurs with N=13, perimeter=26).
apolhexper1uptorot(N,Y) :-
  M is N*2,
  findall(L,polyhexperim1maxrot(L,M),RESULTS),
  discardduplicates(RESULTS,RESULTS2),
  length(RESULTS2,Y).



%% When checking, i.e. rewriting in X -> '111111' direction
%% iterate as long as the Y is not '111111' and the length
%% keeps decreasing.
%% When generating, i.e. rewriting in X <- '111111' direction
%% iterate as long as the result's length is less than given N.

%% Test:
%% Not only ispolyhex([0,1,1,1,1,0,1,0,1,1,1,1,0,1],N).
%% but also ispolyhex([1,1,1,1,0,1,0,1,1,1,1,0,1,0],N). (15738; 32122)
%% should return true, with N=3.
%% ispolyhex([1,1,1,1,0,1,0,1,1,1,1,0,1,0],N). --> N = 3
%% ispolyhex([0,1,1,0,1,1,1,1,0,0,1,1,1,1],N). --> N = 3
%% ispolyhex([1,1,1,1,0,1,0,1,0,1,1,1,1,0,1,0,1,0],N). --> N = 4
%% ispolyhex([1,1,1,1,0,1,0,1,0,1,0,1,1,1,1,0,1,0,1,0,1,0],N). --> N = 5
%% ispolyhex([1,1,1,1,0,1,0,1,0,1,0,1,0,1,1,1,1,0,1,0,1,0,1,0,1,0],N). --> N = 6
%% ispolyhex([1,1,0,1,1,0,1,1,1,0,1,1,0,1],N). --> N = 4
%% This works ONLY if the third rewrite-rule 0110 -> 1001 is present:
%% ispolyhex([1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0],N). --> N = 7 (224694; 486838).


%% One hex is already lopped off:
%% ispolyhex([1,1,0,1,1,0,1,1,1,0,0,1,1,1,0,1,1,0],N). --> N = 6

%% This should work:
%% ispolyhex([1,1,1,1,0,1,1,0,1,1,0,1,1,0,1,1,1,1,0,0,0,0],N). --> N = 5

%% This should not, in no circumstances:
%% ispolyhex([1,1,1,1,0,1,1,1,0,0,1,1,1,0,1,1,1,1,0,0,0,0],N).

%% Neither this should match. However, cases like these prove that
%% the polyhex-language is NOT context-free: (CHECK YOUR TERMINOLOGY!)
%% ispolyhex([1,1,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,1,1,0,0,0,0,0],N).
%% (Currently this returns N=6.)
%% Note that it would not be enough to just filter off patterns with ,0,0,0,0,0, sub-pattern:
%% ispolyhex([1,1,1,1,0,1,1,0,1,0,1,1,0,1,1,0,1,1,0,1,0,1,1,1,1,0,1,0,0,0,0,1,0,0],N).
%% (Currently this returns N=8.)


%% The rewrite-rules follow.
%% hexrewrite(X,Y). will match
%% if X can be rewritten to Y.
%% Fails if no rewriting is possible
%% at any position.

%% 011110B -> 11B (perimeter contracted by four edges)
hexrewrite([0,1,1,1,1,0|X],[1,1|X]).

%% 01110B  -> 101B (perimeter contracted by two edges)
hexrewrite([0,1,1,1,0|X],[1,0,1|X]).

%% 0110B  -> 1001B (perimeter stays the same, one hex lopped off)
hexrewrite([0,1,1,0|X],[1,0,0,1|X]).

%% Search a position where to apply one of the three rules given above.
%% This fails if no such position is found.
hexrewrite([A|X],[A|Y]) :-
  hexrewrite(X,Y).


