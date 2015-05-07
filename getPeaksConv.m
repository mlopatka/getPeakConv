function p = getPeaksConv(varargin)
% getPeaksCov - Probabilistic peak detection of first order chromatographic
% or time-series data. For a provided pair of vectors (time and intensity)
% the function uses probabilistic reasoning to evaluate an exhaustive set 
% of (realistic) peak models in a windowed region. Prior probabilities are
% evaluated by the method published by Davis and Giddings 1983.
% Probabilistic reasoning is used to estimate the posterior probability
% that a given point is affected by a chromatographic peak.
%
% Input:
% p = getPeaksConv(x, y, bandWidth, ResSigma, ALPHA, MAX_PEAKS_PER_ROW)
% 
% x : Time index of each measurement
% y : Intensity value of each measurement
% bandWidth: The standard deviation of a prototypical Gaussian peak for the
% chromatographic system expressed in terms of retention time
% ResSigma: Standard deviation of a blank measurement (baseline noise)
% MAX_PEAKS_PER_ROW (optional): The boound on number of peaks allowed to overlap in a
% chromatographic space define by 2 plates in the chromatographic system.
% Default = 2
% ALPHA (optional): The saturation of the chromatogram [1 > ALPHA > 0]
% Default = 0.25
%
% Output:
% p : Estiamted posterior probibility that the point at p is affected by a
% chromatographic peak
%
% Example usage:
% p = getPeaksConv(retention_time_vector, intenstiy_vector, 2.5, 1.2e-5, 0.3, 3);
%
% If this software is useful to your academic work, please cite our
% publication in lieu of thanks:
%
%%% Lopatka M., Vivó-Truyols G., Sjerps M.J., "Probabilistic peak detection
%%% for first-order chromatographic data." Analitica Chimica Acta. 2014.
%%% DOI: 10.1016/j.aca.2014.02.015
%
% Author: Martin Lopatka <martin.lopatka@gmail.com> Created: 29th August, 2013
% Gabriel Vivó-Truyols <g.vivotruyols@uva.nl> Revised: 23rd April, 2015

%% parse inputs
if numel(varargin) > 6, error('too many input arguments'); end
if numel(varargin) > 5, MAX_PEAKS_PER_ROW = varargin{6}; else MAX_PEAKS_PER_ROW = 2; end
if numel(varargin) > 4, ALPHA = varargin{5}; else ALPHA = 0.25; end
if numel(varargin) < 4
    error('too few input arguments'); 
else
    x = varargin{1}(:); % time
    y = varargin{2}(:); % intensity 
    y = y./max(y);
    bandWidth = varargin{3};
    ResSigma = varargin{4};
end

if ~isequal(numel(x), numel(y))
    error('time and intensity vectors must have an equal number of points');
end

if ~and(isscalar(varargin{3}), isscalar(varargin{4}))
    error('bandwidth and residual noise parameters must be scalar');
end
p = zeros(2,numel(y)); % preallocate 
rt_Delta=mean(diff(x)); % this is the time step constant
n = 4*ceil(bandWidth/rt_Delta); % round to integer values, this is the width of one peak in data points
snMax = exp((-.5*((0/ResSigma).^2))./(ResSigma*sqrt(2*pi))); %optional different threshold than 0 
bandWidth = (bandWidth/rt_Delta); %recast bandwidth in terms of the points rather than in terms of Rt
if mod(n,2) > 0, n = n+1; end % we want n to be even so that the points around i are symetrical
% peaks affect  n+1 points, n points on either side of the centre point i. 
% H_p_i: point i is "in a peak"
% H_d_i: point i is not "in a peak"
% if point i is "in a peak" then the peak centre must be in the interval [i-n/2:i+n/2]

%% precalcualte the possible peak configurations insode
h = waitbar(0.0,'Generating Vm Matrix');
Vm = genVm(((2*n)+1), MAX_PEAKS_PER_ROW); % size of the window is define here as ((2*n)+1) points

%% precalculate the models and store them in a cell
temp = size(Vm,1);
modelHolder = cell(temp,5);
midPoint = round((2*n+1)/2);

for i = 1:temp
    [modelHolder{i,4},modelHolder{i,1}] = Vm2Fun(Vm(1,:), bandWidth); % put the model into a cell array as an anonymous function 
%     modelHolder{i,4} lists the peak centers asumed by the model
%     modelHolder{i,1} holds the model as an annonymous function
    if sum(Vm(1,n/2:(round(1.5*n)))) > 0
       modelHolder{i,2} = true; % this Vm row supports H_p
    else
       modelHolder{i,2} = false;% this Vm row supports H_d
    end
%     modelHolder{i,2} contains true if a peak is present in the inner window and flase otherwise
    modelHolder{i,3} = sum(Vm(1,:));
%     modelHolder{i,3} holds the number of peaks specific to that model.
    if numel(modelHolder{i,4})>1 % checks overlapping peaks
        numPeaks = sum(or(pdist(modelHolder{i,4}', 'cityblock')<=round(n/2), (abs(modelHolder{i,4}-midPoint)<round(midPoint/2))));
        modelHolder{i,5} = ALPHA*dngp(numPeaks,ALPHA);
        if numPeaks == 0
            modelHolder{i,5} = ALPHA*dngp(1,ALPHA); % multiple peaks not overlapping so singlet prior is applied
        end
    else
        modelHolder{i,5} = ALPHA*dngp(1,ALPHA); % overlap not possible
    end
    % application of Davis-Giddings prior  
%     modelHolder{i,5} holds the the davis and giddings prior value 
    waitbar((0.05+((i/temp)*0.95)),h);
    Vm(1,:)=[]; %delete the first index in Vm so that we are clearing space as we process the matric row by row
end

modelHolder{1,5} = (ALPHA*dngp(0,ALPHA)); % zero peaks model prior

% housekeeping
clear temp;
if isempty(Vm), clear Vm; else error('indexing problem, not all VmRows considered'); end
close(h); %clean up the waitbar

indM = size(modelHolder,1);
x_t = ([midPoint-n:midPoint+n]-midPoint)';%linspace(-1,1,2*n+1)';%;
% x_t = 1:(2*n+1);
yy=repmat(y',numel(x_t),1);
i=1:numel(x_t);
%temp=fliplr((numel(x_t)-1)./2-(i-1));
temp=(numel(x_t)-1)./2-(i-1);
for i=1:numel(x_t)
   yy(i,:)=circshift(yy(i,:),[0 temp(i)]);
end

%% normalize similar models proportionately!
% This is new strategy implemented after the publication!
prior_vals = unique([modelHolder{:,5}]);            
for i = 1:numel(prior_vals)
%     temp = sum([modelHolder{[modelHolder{:,3}]==i, 5}]); % the sum of priors
    temp = sum([modelHolder{:,5}]==prior_vals(i)); % number of similar models
    modelHolder([modelHolder{:,5}]==prior_vals(i), 5) = cellfun(@(v) v./temp,modelHolder([modelHolder{:,5}]==prior_vals(i), 5) ,'UniformOutput', false); %force return as a cellend
end 

%normalize to sum to unity
temp = sum([modelHolder{:,5}]);
modelHolder(:,5) = cellfun(@(v) v./temp,modelHolder(:,5),'UniformOutput', false); %force return as a cellend

PRIOR_PEAK = sum([modelHolder{2:end,5}]);
PRIOR_NOPEAK = 1-PRIOR_PEAK;
%preallocate
p_of_d_given_Hp = zeros([1, numel(y)]);
p_of_d_given_Hd = zeros([1, numel(y)]);

h = waitbar(0.0,'performing pointwise model evaluation');
mask_holder = [];

for m = 1:size(modelHolder,1)%iterate over all hypotheses
    mH = modelHolder(m,:);

    X = exportXMatrix(x_t, mH{3}, mH{4},bandWidth);
    p = pinv(X'*X)*X';
    
    b=ones(size(X,2),numel(y));
    for i=1:size(p,1)
        b(i,:)=conv(y',fliplr(p(i,:)),'same');
    end
    y_hat=X*b;
    res = yy-y_hat;
    
    %logResProb = log(normpdf(res,0,ResSigma));
    logResProb = log(1/(sqrt(2*pi)*ResSigma))-0.5.*((res./ResSigma).^2);
       
    if m==1
        mask = true(size(y'));
    else
        mask = sum((b(3:end,:)*snMax<0),1)==0;
    end
    
    mask_holder(m,1:numel(b(:,422))) = b(:,422)';
                % PRIOR          % LIKELIHOOD          % MASK
    model_prob = modelHolder{m,5}*exp(sum(logResProb,1).*mask); 
    % remove models fitting a negative (or super tiny if we want) 
    % gaussian as it is impossible for chromatography

    if mH{2}==true %hp supporting hypothesis
        p_of_d_given_Hp = p_of_d_given_Hp(~isnan(p_of_d_given_Hp)) + model_prob;
    else
        p_of_d_given_Hd = p_of_d_given_Hd(~isnan(p_of_d_given_Hd)) + model_prob;
    end

    waitbar(m/indM,h);
end     

delete(get(h, 'children'))
close(h)

p = (p_of_d_given_Hp*PRIOR_PEAK)./((p_of_d_given_Hp*PRIOR_PEAK)+(p_of_d_given_Hd*(PRIOR_NOPEAK)));%posterior probability

p(isnan(p)) = 0; % clean up NaN vlaues
p(1:n) = 0; % unreliable results in the edges are set to zeros
p(end-n:end) = 0;

end

function Vm = genVm(N, maxP)
    % Vn = [true false]; 
    % N = window size in indeces;
    % first we need to preallocate based on our maxP, this is a factorial
    % N^C_maxP = N!/(N-maxP)!maxP!
    % otherwise this goes totally insane where n>22 leads to out of memory
    b = 0;
    for i = 1:maxP
        b = b + factorial(N)/((factorial(N-i)*factorial(i))); % calculate the number of permutations for exactly each number of peaks
    end

    Vm = false([uint64(b)+1,N]); %preallocate - this speeds things up very nicely.
    counter = 2; % ensure that first row in Vm is zeros (no peaks)

    for i = 1:maxP
        indX = nchoosek([1:N],i); %generate the INDECES that will be set to '1' on this iteration of the inner loop
        for j = 1:size(indX,1)
            Vm(counter, [indX(j,:)]) = true; %operate on the rows by logical indexing regardless of the size of indX(j,:) this is fast
            counter = counter+1; %incremetn the counter -- this is guarenteed never to exceed b
        end
        indX=[];
    end
end

function p = dngp(n,alpha)
%davis and giddings prior propbability for peak configurations given a
%particular saturation parameter.
% n is the overlapping numer of peaks i.e. 2 for a double, 3 for a triplet
% alpha is an estimate of the saturation of the chromatographic space
    e = exp(1);
    p = (e^(-2*alpha))*(1-e^(-alpha))^(n-1);
end

function [mu,yHatFun] = Vm2Fun(vMRow, bandWidth)
clear yHatFun;
    baseFun = 'yHatFun = @(a,b';
    mu = find(vMRow); % mu is a vector of means in the VMRow, 
                      % each corresponds the mean of a gaussian 
                      % component or Peak centre in the hypothesis 
                      % defioned by VmRow.
    for i = 1:sum(vMRow) % adding Gaussian components 
        baseFun = [baseFun ,',c',num2str(i)];
    end
    baseFun = [baseFun,',x) (a+(b*x)']; % add linear baseline to the model.
    
    if sum(vMRow) > 0
        for i = 1:sum(vMRow) % adding a variable number of gaussians components 
                             % and their respective c terms depending on the 
                             % specific hypothesis depicted in VMRow.
             baseFun = [baseFun, '+(c',num2str(i),'*(exp(-0.5*((x-',num2str(mu(i)),')/',num2str(bandWidth),').^2)./(',num2str(bandWidth),'*sqrt(2*pi))))'];   
        end
    else
        baseFun =  'yHatFun = @(a,b,x) (a + (b*x)'; %in case there are no peaks
    end
        baseFun = [baseFun,');'];
        eval(baseFun); % yHatFun is now an anonymous function that will be fit to our data         
end

function A = exportXMatrix(j, model, mus, band)
    % use the defined bandwidth of a peak
    snoggy = @(x, mu, sig)(exp((-.5*((x - mu)/sig).^2))./(sig*sqrt(2*pi))); % snoggy is gaussian
    model = model+2;
    A = repmat([1:numel(j)]', [1, model]); %Initialize the system matrix
    A(:,1) = 1; % populate the parameter column of the system matrix for the intercept term
    for i = 1:numel(mus)
        A(:,i+2) = snoggy([1:numel(j)]',mus(i),band); %populate the parameter column for each Gaussian componnent as a black box
    end % (s)noggy in honour of the Nuge
    %numerical solution for the least squares fit problem
end

