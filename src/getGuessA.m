function A_guess=getGuessA(degree, interval)

%Note that these guesses are in [0,1]
guess=load(strcat('guesses/solutionDeg',num2str(degree),'.mat'));
A_guess=guess.A;

if(interval=="m11")
    A_guess=convertAFrom00toM11(A_guess);
end

end