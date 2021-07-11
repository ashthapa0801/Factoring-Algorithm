%% Only works for 4 orders properly!!
% funcT = mod(T,M) // remainder = mod(a^n,m) // b = a^n - floor(a^n./m)*m
% 7^t mod (100)
clc;
clear;

A = 7;                      % Seed T
n = 1024;                   % input("Please Enter the Bit size: ");  P                 % Bit 128, 256, 512, 1024
B = 6268263152346622;       % input('Please Enter the Number to be Factored: ');  M    % Number to be factorized

t = sym(1):1:sym(n);                    % Signal duration/length
z = A.^t;
fncT = (z) - floor(z/B)*B;              % Alternative version to mod(a,b) //powermod(z,B,n);

% Fast fourier Transforms 1-D
% Rearranges 'fncT' by shifting the zero-frequency component to the centre of the array
FTsignal = fft(double(fncT) - mean(double(fncT)))/numel(double(fncT));   
Z = fftshift(FTsignal);                 % Computes descrete fourier transform of Z using fft algorithm
Z_mag = abs(Z.*Z);                      % Returns absolute value of each element 

% Plotting Graphs
title('Register in superposition state')
tiledlayout(4,1)                        % 4 rows meaning 4 different graphs
Z_ax2 = nexttile;                       % Creates an axes object and inputs the graph into the next empty tile of the layout
plot(Z_ax2,fncT)
xlim([0 n+1])                           % Sets the x axis limit
ylabel('Amplitute')
Z_ax3 = nexttile;
stem(Z_ax3,fncT)                        % Displays only the discrete values of the points 
xlim([0 n+1])
ylabel('Amplitute')
xlabel('Signal Duration')
Z_ax = nexttile;
plot(Z_ax,Z_mag);
xlim([0 n+1])
ylabel('Large Zero frequency term')
Z_ax1 = nexttile;
stem(Z_ax1,Z_mag);
xlim([0 n+1])
ylabel('Large Zero frequency term')
xlabel('Signal Duration')
figure(1);

% Check if the number is prime and if composite, calculate the potential divisors 
prime=1;
divisors=[];                        % Empty vector to store potential divisors
for i=2:sqrt(B)
    if rem(B,i)==0
        prime=0;
        divisors = [divisors,i];    % Store the divisors of B
    end
end

if prime==1
    disp([num2str(B),' is prime and it only has factors of 1 and ',num2str(B),'.']);
    return;
else
    disp('The number is Composite.');
end 

% The main part of the algorithm that finds the more factors/divisors 
[pks1,pks2,pks3] = findpeaks(Z_mag);    % Finds the peaks of graphs 3 or 4

if (pks2(1) == 33)                      % For order of 4
    x1 = sqrt((pks2(1)-1)/2);           % (-1) because the graph displays 33 instead of 32. 32 is the correct one it should show.   
    eqn = (A^(x1/2)) + 1;
    eqn1 = (A^(x1/2)) - 1;
    G = divisors;
    emp = zeros(n);                     % n-by-n matrix of zeros
    
    if rem(x1,1)==0                                 % If x1 is an integer proceed with this part
        % Eculideans Algorithm & Greatest Common Divisor
        Ba = B;
        G1 = gcd(eqn,Ba);
        G2 = gcd(eqn1,Ba);
        G3 = sqrt(x1);
        G3a = sqrt(Ba);
        G3s = [G3 G3a];
    end
        
        if rem(G3a,1)==0                                % If G3a is a integer then do not change the number.
            G3a = G3a;
        else
            G3a = [1];                              % If G3a is not a int then assign a value of 1
        end
                                            
        G3s = G3s(randperm(length(G3s)));           % Randomize the order
        G3s(fix(G3s)~=G3s) = [];                    % Remove non integers
        
        G4 = G3s(randperm(length(G3s)));            % Randomize the order
        G4(fix(G4)~=G4) = [];
        G4 = [G4,Ba./G4];                           % Division
        G4 = G4(randperm(length(G4)));
        G4(fix(G4)~=G4) = [];
        
        maxlen5 = max(length(B), length(G4));       % Make cells G3s and G4 the same length
        Ba = [zeros(1,maxlen5 - length(Ba)), Ba];
        Ba = Ba(randperm(length(Ba)));
        G4 = [zeros(1,maxlen5 - length(G4)), G4];
        G5 = gcd(Ba,G4);
        G5(fix(G5)~=G5) = [];
        G5 = G5(randperm(length(G5)));
        
        maxlen = max(length(G3s), length(G4));      % Make cells G3s and G4 the same length
        G3s = [zeros(1,maxlen - length(G3s)), G3s];
        G4 = [zeros(1,maxlen - length(G4)), G4];
        G6 = gcd(G3s,G4);
        G6(fix(G6)~=G6) = [];
        G6 = G6(randperm(length(G6)));
        
        maxlen1 = max(length(G5), length(G6));      % Make cells G5 and G6 the same length
        G5 = [zeros(1,maxlen1 - length(G5)), G5];
        G6 = [zeros(1,maxlen1 - length(G6)), G6];
        G7 = gcd(G5,G6);
        G7(fix(G7)~=G7) = [];   
        G7 = G7(randperm(length(G7)));
        
        maxlen2 = max(length(G6), length(G7));      % Make cells G6 and G7 the same length
        G7 = [zeros(1,maxlen2 - length(G7)), G7];
        G8 = gcd(G6,G7);
        G8(fix(G8)~=G8) = [];
        G8 = G8(randperm(length(G8)));
        
        maxlen3 = max(length(G7), length(G8));      % Make cells G7 and G8 the same length
        G8 = [zeros(1,maxlen3 - length(G8)), G8];
        G9 = gcd(G7,G8);
        G9(fix(G9)~=G9) = [];
        G9 = G9(randperm(length(G9)));
        
        maxlen4 = max(length(G7), length(G9));      % Make cells G5 and G6 the same length
        G7 = [zeros(1,maxlen4 - length(G7)), G7];
        G9 = [zeros(1,maxlen4 - length(G9)), G9];
        G10 = gcd(G7,G9);
        
        G11 = gcd(100,G5); 
        G11a = G11(randperm(length(G11)));
        G11a(fix(G11a)~=G11a) = [];
        Ca = G11a;
        
        maxlen6 = max(length(Ba), length(Ca));      % Make cells Ba and G6 the same length
        Ba = [zeros(1,maxlen6 - length(B)), Ba];
        Ca = [zeros(1,maxlen6 - length(Ca)), Ca];
        G12 = [Ca, B./Ca];
        
        maxlen7 = max(length(G), length(G3s));      % Make cells G and G3s the same length
        G = [zeros(1,maxlen7 - length(G)), G];
        G3s = [zeros(1,maxlen7 - length(G3s)), G3s];
        G13 = gcd(G,G3s);
        G13 = G13(randperm(length(G13)));           % Randomize the order
        G13(fix(G13)~=G13) = [];                    % Remove Non integers
        
        maxlen8 = max(length(G), length(G4));       % Make cells G and G4 the same length
        G = [zeros(1,maxlen8 - length(G)), G];
        G4 = [zeros(1,maxlen8 - length(G4)), G4];
        G14 = gcd(G,G4);
        G14 = G14(randperm(length(G14)));           % Randomize the order
        G14(fix(G14)~=G14) = [];                    % Remove Non integers
        
        maxlen9 = max(length(G), length(G5));       % Make cells G and G5 the same length
        G = [zeros(1,maxlen9 - length(G)), G];
        G5 = [zeros(1,maxlen9 - length(G5)), G5];
        G15 = gcd(G,G5);
        G15 = G15(randperm(length(G15)));           % Randomize the order
        G15(fix(G15)~=G15) = [];                    % Remove Non integers
        
        maxlen9 = max(length(G), length(G6));       % Make cells G and G6 the same length
        G = [zeros(1,maxlen9 - length(G)), G];
        G6 = [zeros(1,maxlen9 - length(G6)), G6];
        G16 = gcd(G,G6);
        G16 = G16(randperm(length(G16)));           % Randomize the order
        G16(fix(G16)~=G16) = [];                    % Remove Non integers
        
        maxlen10 = max(length(G), length(G7));      % Make cells G and G7 the same length
        G = [zeros(1,maxlen10 - length(G)), G];
        G7 = [zeros(1,maxlen10 - length(G7)), G7];
        G17 = gcd(G,G7);
        G17 = G17(randperm(length(G17)));           % Randomize the order
        G17(fix(G17)~=G17) = [];                    % Remove Non integers
        
        maxlen11 = max(length(G), length(G9));      % Make cells G and G7 the same length
        G = [zeros(1,maxlen11 - length(G)), G];
        G9 = [zeros(1,maxlen11 - length(G9)), G9];
        G18 = gcd(G,G9);
        G18 = G18(randperm(length(G18)));           % Randomize the order
        G18(fix(G18)~=G18) = [];                    % Remove Non integers
        
        maxlen12 = max(length(G), length(G10));      % Make cells G and G7 the same length
        G = [zeros(1,maxlen12 - length(G)), G];
        G10 = [zeros(1,maxlen12 - length(G10)), G10];
        G19 = gcd(G,G9);
        G19 = G19(randperm(length(G19)));           % Randomize the order
        G19(fix(G19)~=G19) = [];                    % Remove Non integers
        
        maxlen12 = max(length(G), length(G11));      % Make cells G and G7 the same length
        G = [zeros(1,maxlen12 - length(G)), G];
        G11 = [zeros(1,maxlen12 - length(G11)), G11];
        G20 = gcd(G,G11);
        G20 = G20(randperm(length(G20)));           % Randomize the order
        G20(fix(G20)~=G20) = [];    
        
        G21 = [Ca, B./Ca];
        G22 = [G, B./G];
        
elseif pks2(1) == 17                            % Note 17 is 17-16 because of the error in FFT shift
    x1 = sqrt((pks2(1)-1));
    eqn = (A^(x1/2)) + 1;
    eqn1 = (A^(x1/2)) - 1;
    G = divisors;
    
    if rem(x1,1)==0                             % If x1 is an integer proceed with this part
        % Eculideans Algorithm & Greatest Common Divisor
        Ba = B;
        G1 = gcd(eqn,Ba);
        G2 = gcd(eqn1,Ba);
        G3 = sqrt(x1);
        G3a = sqrt(Ba);
        G3s = [G3 G3a];
    end
        
        if rem(G3a,1)==0                            % If G3a is a integer then do not change the number.
            G3a = G3a;
        else
            G3a = [1];                              % If G3a is not a int then assign a value of 1
        end
                                            
        G3s = G3s(randperm(length(G3s)));           % Randomize the order
        G3s(fix(G3s)~=G3s) = [];                    % Remove non integers
        
        G4 = G3s(randperm(length(G3s)));            % Randomize the order
        G4(fix(G4)~=G4) = [];
        G4 = [G4,Ba./G4];                           % Division
        G4 = G4(randperm(length(G4)));
        G4(fix(G4)~=G4) = [];
        
        maxlen5 = max(length(B), length(G4));       % Make cells G3s and G4 the same length
        Ba = [zeros(1,maxlen5 - length(Ba)), Ba];
        Ba = Ba(randperm(length(Ba)));
        G4 = [zeros(1,maxlen5 - length(G4)), G4];
        G5 = gcd(Ba,G4);
        G5(fix(G5)~=G5) = [];
        G5 = G5(randperm(length(G5)));
        
        maxlen = max(length(G3s), length(G4));      % Make cells G3s and G4 the same length
        G3s = [zeros(1,maxlen - length(G3s)), G3s];
        G4 = [zeros(1,maxlen - length(G4)), G4];
        G6 = gcd(G3s,G4);
        G6(fix(G6)~=G6) = [];
        G6 = G6(randperm(length(G6)));
        
        maxlen1 = max(length(G5), length(G6));      % Make cells G5 and G6 the same length
        G5 = [zeros(1,maxlen1 - length(G5)), G5];
        G6 = [zeros(1,maxlen1 - length(G6)), G6];
        G7 = gcd(G5,G6);
        G7(fix(G7)~=G7) = [];   
        G7 = G7(randperm(length(G7)));
        
        maxlen2 = max(length(G6), length(G7));      % Make cells G6 and G7 the same length
        G7 = [zeros(1,maxlen2 - length(G7)), G7];
        G8 = gcd(G6,G7);
        G8(fix(G8)~=G8) = [];
        G8 = G8(randperm(length(G8)));
        
        maxlen3 = max(length(G7), length(G8));      % Make cells G7 and G8 the same length
        G8 = [zeros(1,maxlen3 - length(G8)), G8];
        G9 = gcd(G7,G8);
        G9(fix(G9)~=G9) = [];
        G9 = G9(randperm(length(G9)));
        
        maxlen4 = max(length(G7), length(G9));      % Make cells G5 and G6 the same length
        G7 = [zeros(1,maxlen4 - length(G7)), G7];
        G9 = [zeros(1,maxlen4 - length(G9)), G9];
        G10 = gcd(G7,G9);
        
        G11 = gcd(100,G5); 
        G11a = G11(randperm(length(G11)));
        G11a(fix(G11a)~=G11a) = [];
        Ca = G11a;
        
        maxlen6 = max(length(Ba), length(Ca));      % Make cells Ba and G6 the same length
        Ba = [zeros(1,maxlen6 - length(B)), Ba];
        Ca = [zeros(1,maxlen6 - length(Ca)), Ca];
        G12 = [Ca, B./Ca];
        
        maxlen7 = max(length(G), length(G3s));      % Make cells G and G3s the same length
        G = [zeros(1,maxlen7 - length(G)), G];
        G3s = [zeros(1,maxlen7 - length(G3s)), G3s];
        G13 = gcd(G,G3s);
        G13 = G13(randperm(length(G13)));           % Randomize the order
        G13(fix(G13)~=G13) = [];                    % Remove Non integers
        
        maxlen8 = max(length(G), length(G4));       % Make cells G and G4 the same length
        G = [zeros(1,maxlen8 - length(G)), G];
        G4 = [zeros(1,maxlen8 - length(G4)), G4];
        G14 = gcd(G,G4);
        G14 = G14(randperm(length(G14)));           % Randomize the order
        G14(fix(G14)~=G14) = [];                    % Remove Non integers
        
        maxlen9 = max(length(G), length(G5));       % Make cells G and G5 the same length
        G = [zeros(1,maxlen9 - length(G)), G];
        G5 = [zeros(1,maxlen9 - length(G5)), G5];
        G15 = gcd(G,G5);
        G15 = G15(randperm(length(G15)));           % Randomize the order
        G15(fix(G15)~=G15) = [];                    % Remove Non integers
        
        maxlen9 = max(length(G), length(G6));       % Make cells G and G6 the same length
        G = [zeros(1,maxlen9 - length(G)), G];
        G6 = [zeros(1,maxlen9 - length(G6)), G6];
        G16 = gcd(G,G6);
        G16 = G16(randperm(length(G16)));           % Randomize the order
        G16(fix(G16)~=G16) = [];                    % Remove Non integers
        
        maxlen10 = max(length(G), length(G7));      % Make cells G and G7 the same length
        G = [zeros(1,maxlen10 - length(G)), G];
        G7 = [zeros(1,maxlen10 - length(G7)), G7];
        G17 = gcd(G,G7);
        G17 = G17(randperm(length(G17)));           % Randomize the order
        G17(fix(G17)~=G17) = [];                    % Remove Non integers
        
        maxlen11 = max(length(G), length(G9));      % Make cells G and G7 the same length
        G = [zeros(1,maxlen11 - length(G)), G];
        G9 = [zeros(1,maxlen11 - length(G9)), G9];
        G18 = gcd(G,G9);
        G18 = G18(randperm(length(G18)));           % Randomize the order
        G18(fix(G18)~=G18) = [];                    % Remove Non integers
        
        maxlen12 = max(length(G), length(G10));      % Make cells G and G7 the same length
        G = [zeros(1,maxlen12 - length(G)), G];
        G10 = [zeros(1,maxlen12 - length(G10)), G10];
        G19 = gcd(G,G9);
        G19 = G19(randperm(length(G19)));           % Randomize the order
        G19(fix(G19)~=G19) = [];                    % Remove Non integers
        
        maxlen12 = max(length(G), length(G11));      % Make cells G and G7 the same length
        G = [zeros(1,maxlen12 - length(G)), G];
        G11 = [zeros(1,maxlen12 - length(G11)), G11];
        G20 = gcd(G,G11);
        G20 = G20(randperm(length(G20)));           % Randomize the order
        G20(fix(G20)~=G20) = [];    
        
        G21 = [Ca, B./Ca];
        G22 = [G, B./G];
        
elseif pks2(3) == 25                            % For orders other than 7 
        x1 = sqrt((pks2(1)));
        eqn = (A^(x1/3)) + 1;
        eqn1 = (A^(x1/3)) - 1;
        G = divisors;
        
        if rem(x1,1)==0                             % If x1 is an integer proceed with this part
        % Eculideans Algorithm & Greatest Common Divisor
        Ba = B;
        G1 = gcd(eqn,Ba);
        G2 = gcd(eqn1,Ba);
        G3 = sqrt(x1);
        G3a = sqrt(Ba);
        G3s = [G3 G3a];
        end
        
        if rem(G3a,1)==0                            % If G3a is a integer then do not change the number.
            G3a = G3a;
        else
            G3a = [1];                              % If G3a is not a int then assign a value of 1
        end
                                            
        G3s = G3s(randperm(length(G3s)));           % Randomize the order
        G3s(fix(G3s)~=G3s) = [];                    % Remove non integers
        
        G4 = G3s(randperm(length(G3s)));            % Randomize the order
        G4(fix(G4)~=G4) = [];
        G4 = [G4,Ba./G4];                           % Division
        G4 = G4(randperm(length(G4)));
        G4(fix(G4)~=G4) = [];
        
        maxlen5 = max(length(B), length(G4));       % Make cells G3s and G4 the same length
        Ba = [zeros(1,maxlen5 - length(Ba)), Ba];
        Ba = Ba(randperm(length(Ba)));
        G4 = [zeros(1,maxlen5 - length(G4)), G4];
        G5 = gcd(Ba,G4);
        G5(fix(G5)~=G5) = [];
        G5 = G5(randperm(length(G5)));
        
        maxlen = max(length(G3s), length(G4));      % Make cells G3s and G4 the same length
        G3s = [zeros(1,maxlen - length(G3s)), G3s];
        G4 = [zeros(1,maxlen - length(G4)), G4];
        G6 = gcd(G3s,G4);
        G6(fix(G6)~=G6) = [];
        G6 = G6(randperm(length(G6)));
        
        maxlen1 = max(length(G5), length(G6));      % Make cells G5 and G6 the same length
        G5 = [zeros(1,maxlen1 - length(G5)), G5];
        G6 = [zeros(1,maxlen1 - length(G6)), G6];
        G7 = gcd(G5,G6);
        G7(fix(G7)~=G7) = [];   
        G7 = G7(randperm(length(G7)));
        
        maxlen2 = max(length(G6), length(G7));      % Make cells G6 and G7 the same length
        G7 = [zeros(1,maxlen2 - length(G7)), G7];
        G8 = gcd(G6,G7);
        G8(fix(G8)~=G8) = [];
        G8 = G8(randperm(length(G8)));
        
        maxlen3 = max(length(G7), length(G8));      % Make cells G7 and G8 the same length
        G8 = [zeros(1,maxlen3 - length(G8)), G8];
        G9 = gcd(G7,G8);
        G9(fix(G9)~=G9) = [];
        G9 = G9(randperm(length(G9)));
        
        maxlen4 = max(length(G7), length(G9));      % Make cells G5 and G6 the same length
        G7 = [zeros(1,maxlen4 - length(G7)), G7];
        G9 = [zeros(1,maxlen4 - length(G9)), G9];
        G10 = gcd(G7,G9);
        
        G11 = gcd(100,G5); 
        G11a = G11(randperm(length(G11)));
        G11a(fix(G11a)~=G11a) = [];
        Ca = G11a;
        
        maxlen6 = max(length(Ba), length(Ca));      % Make cells Ba and G6 the same length
        Ba = [zeros(1,maxlen6 - length(B)), Ba];
        Ca = [zeros(1,maxlen6 - length(Ca)), Ca];
        G12 = [Ca, B./Ca];
        
        maxlen7 = max(length(G), length(G3s));      % Make cells G and G3s the same length
        G = [zeros(1,maxlen7 - length(G)), G];
        G3s = [zeros(1,maxlen7 - length(G3s)), G3s];
        G13 = gcd(G,G3s);
        G13 = G13(randperm(length(G13)));           % Randomize the order
        G13(fix(G13)~=G13) = [];                    % Remove Non integers
        
        maxlen8 = max(length(G), length(G4));       % Make cells G and G4 the same length
        G = [zeros(1,maxlen8 - length(G)), G];
        G4 = [zeros(1,maxlen8 - length(G4)), G4];
        G14 = gcd(G,G4);
        G14 = G14(randperm(length(G14)));           % Randomize the order
        G14(fix(G14)~=G14) = [];                    % Remove Non integers
        
        maxlen9 = max(length(G), length(G5));       % Make cells G and G5 the same length
        G = [zeros(1,maxlen9 - length(G)), G];
        G5 = [zeros(1,maxlen9 - length(G5)), G5];
        G15 = gcd(G,G5);
        G15 = G15(randperm(length(G15)));           % Randomize the order
        G15(fix(G15)~=G15) = [];                    % Remove Non integers
        
        maxlen9 = max(length(G), length(G6));       % Make cells G and G6 the same length
        G = [zeros(1,maxlen9 - length(G)), G];
        G6 = [zeros(1,maxlen9 - length(G6)), G6];
        G16 = gcd(G,G6);
        G16 = G16(randperm(length(G16)));           % Randomize the order
        G16(fix(G16)~=G16) = [];                    % Remove Non integers
        
        maxlen10 = max(length(G), length(G7));      % Make cells G and G7 the same length
        G = [zeros(1,maxlen10 - length(G)), G];
        G7 = [zeros(1,maxlen10 - length(G7)), G7];
        G17 = gcd(G,G7);
        G17 = G17(randperm(length(G17)));           % Randomize the order
        G17(fix(G17)~=G17) = [];                    % Remove Non integers
        
        maxlen11 = max(length(G), length(G9));      % Make cells G and G7 the same length
        G = [zeros(1,maxlen11 - length(G)), G];
        G9 = [zeros(1,maxlen11 - length(G9)), G9];
        G18 = gcd(G,G9);
        G18 = G18(randperm(length(G18)));           % Randomize the order
        G18(fix(G18)~=G18) = [];                    % Remove Non integers
        
        maxlen12 = max(length(G), length(G10));      % Make cells G and G7 the same length
        G = [zeros(1,maxlen12 - length(G)), G];
        G10 = [zeros(1,maxlen12 - length(G10)), G10];
        G19 = gcd(G,G9);
        G19 = G19(randperm(length(G19)));           % Randomize the order
        G19(fix(G19)~=G19) = [];                    % Remove Non integers
        
        maxlen12 = max(length(G), length(G11));      % Make cells G and G7 the same length
        G = [zeros(1,maxlen12 - length(G)), G];
        G11 = [zeros(1,maxlen12 - length(G11)), G11];
        G20 = gcd(G,G11);
        G20 = G20(randperm(length(G20)));           % Randomize the order
        G20(fix(G20)~=G20) = [];    
        
        G21 = [Ca, B./Ca];
        G22 = [G, B./G];
        
else
        G = divisors;
        % Eculideans Algorithm & Greatest Common Divisor
        G1 = 0;
        G2 = 0;
        Ba = B;
        G3 = sqrt(pks2(1));
        G3a = sqrt(Ba);
        G3s = [G3 G3a];
        
        if rem(G3a,1)==0                            % If G3a is a integer then do not change the number.
            G3a = G3a;
        else
            G3a = [1];                              % If G3a is not a int then assign a value of 1
        end
                                            
        G3s = G3s(randperm(length(G3s)));           % Randomize the order
        G3s(fix(G3s)~=G3s) = [];                    % Remove non integers
        
        G4 = G3s(randperm(length(G3s)));            % Randomize the order
        G4(fix(G4)~=G4) = [];
        G4 = [G4,Ba./G4];                           % Division
        G4 = G4(randperm(length(G4)));
        G4(fix(G4)~=G4) = [];
        
        maxlen5 = max(length(B), length(G4));       % Make cells G3s and G4 the same length
        Ba = [zeros(1,maxlen5 - length(Ba)), Ba];
        Ba = Ba(randperm(length(Ba)));
        G4 = [zeros(1,maxlen5 - length(G4)), G4];
        G5 = gcd(Ba,G4);
        G5(fix(G5)~=G5) = [];
        G5 = G5(randperm(length(G5)));
        
        maxlen = max(length(G3s), length(G4));      % Make cells G3s and G4 the same length
        G3s = [zeros(1,maxlen - length(G3s)), G3s];
        G4 = [zeros(1,maxlen - length(G4)), G4];
        G6 = gcd(G3s,G4);
        G6(fix(G6)~=G6) = [];
        G6 = G6(randperm(length(G6)));
        
        maxlen1 = max(length(G5), length(G6));      % Make cells G5 and G6 the same length
        G5 = [zeros(1,maxlen1 - length(G5)), G5];
        G6 = [zeros(1,maxlen1 - length(G6)), G6];
        G7 = gcd(G5,G6);
        G7(fix(G7)~=G7) = [];   
        G7 = G7(randperm(length(G7)));
        
        maxlen2 = max(length(G6), length(G7));      % Make cells G6 and G7 the same length
        G7 = [zeros(1,maxlen2 - length(G7)), G7];
        G8 = gcd(G6,G7);
        G8(fix(G8)~=G8) = [];
        G8 = G8(randperm(length(G8)));
        
        maxlen3 = max(length(G7), length(G8));      % Make cells G7 and G8 the same length
        G8 = [zeros(1,maxlen3 - length(G8)), G8];
        G9 = gcd(G7,G8);
        G9(fix(G9)~=G9) = [];
        G9 = G9(randperm(length(G9)));
        
        maxlen4 = max(length(G7), length(G9));      % Make cells G5 and G6 the same length
        G7 = [zeros(1,maxlen4 - length(G7)), G7];
        G9 = [zeros(1,maxlen4 - length(G9)), G9];
        G10 = gcd(G7,G9);
        
        G11 = gcd(100,G5); 
        G11a = G11(randperm(length(G11)));
        G11a(fix(G11a)~=G11a) = [];
        Ca = G11a;
        
        maxlen6 = max(length(Ba), length(Ca));      % Make cells Ba and G6 the same length
        Ba = [zeros(1,maxlen6 - length(B)), Ba];
        Ca = [zeros(1,maxlen6 - length(Ca)), Ca];
        G12 = [Ca, B./Ca];
        
        maxlen7 = max(length(G), length(G3s));      % Make cells G and G3s the same length
        G = [zeros(1,maxlen7 - length(G)), G];
        G3s = [zeros(1,maxlen7 - length(G3s)), G3s];
        G13 = gcd(G,G3s);
        G13 = G13(randperm(length(G13)));           % Randomize the order
        G13(fix(G13)~=G13) = [];                    % Remove Non integers
        
        maxlen8 = max(length(G), length(G4));       % Make cells G and G4 the same length
        G = [zeros(1,maxlen8 - length(G)), G];
        G4 = [zeros(1,maxlen8 - length(G4)), G4];
        G14 = gcd(G,G4);
        G14 = G14(randperm(length(G14)));           % Randomize the order
        G14(fix(G14)~=G14) = [];                    % Remove Non integers
        
        maxlen9 = max(length(G), length(G5));       % Make cells G and G5 the same length
        G = [zeros(1,maxlen9 - length(G)), G];
        G5 = [zeros(1,maxlen9 - length(G5)), G5];
        G15 = gcd(G,G5);
        G15 = G15(randperm(length(G15)));           % Randomize the order
        G15(fix(G15)~=G15) = [];                    % Remove Non integers
        
        maxlen9 = max(length(G), length(G6));       % Make cells G and G6 the same length
        G = [zeros(1,maxlen9 - length(G)), G];
        G6 = [zeros(1,maxlen9 - length(G6)), G6];
        G16 = gcd(G,G6);
        G16 = G16(randperm(length(G16)));           % Randomize the order
        G16(fix(G16)~=G16) = [];                    % Remove Non integers
        
        maxlen10 = max(length(G), length(G7));      % Make cells G and G7 the same length
        G = [zeros(1,maxlen10 - length(G)), G];
        G7 = [zeros(1,maxlen10 - length(G7)), G7];
        G17 = gcd(G,G7);
        G17 = G17(randperm(length(G17)));           % Randomize the order
        G17(fix(G17)~=G17) = [];                    % Remove Non integers
        
        maxlen11 = max(length(G), length(G9));      % Make cells G and G7 the same length
        G = [zeros(1,maxlen11 - length(G)), G];
        G9 = [zeros(1,maxlen11 - length(G9)), G9];
        G18 = gcd(G,G9);
        G18 = G18(randperm(length(G18)));           % Randomize the order
        G18(fix(G18)~=G18) = [];                    % Remove Non integers
        
        maxlen12 = max(length(G), length(G10));      % Make cells G and G7 the same length
        G = [zeros(1,maxlen12 - length(G)), G];
        G10 = [zeros(1,maxlen12 - length(G10)), G10];
        G19 = gcd(G,G9);
        G19 = G19(randperm(length(G19)));           % Randomize the order
        G19(fix(G19)~=G19) = [];                    % Remove Non integers
        
        maxlen12 = max(length(G), length(G11));      % Make cells G and G7 the same length
        G = [zeros(1,maxlen12 - length(G)), G];
        G11 = [zeros(1,maxlen12 - length(G11)), G11];
        G20 = gcd(G,G11);
        G20 = G20(randperm(length(G20)));           % Randomize the order
        G20(fix(G20)~=G20) = [];    
        
        G21 = [Ca, B./Ca];
        G22 = [G, B./G];
        
end

        
AllCell = [G1 G2 G3 G3a G4 G5 G6 G7 G8 G9 G10 G11 G11a G12 G13 G14 G15 G16 G17 G18 G19 G20 G21 G22 G Ca Ba]; % Put into array
arr = AllCell;                             
sor_arr = [sort(arr)];                      % Order numbers in that array               
sor_arr(isinf(sor_arr)|isnan(sor_arr)) = 0; % Replaces inf or nan with 0 to remove it on the next line of code
sor_arr(sor_arr<=1)=1;                      % Remove any thing below integer 1
sor_arr(fix(sor_arr)~=sor_arr) = [];        % Removes any non integer
uni_no = unique(sor_arr);                   % Remove duplicates
% Display factors
disp(['The factors of ',num2str(B),' are ', num2str(uni_no)])
    

