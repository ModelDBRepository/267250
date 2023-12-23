%%% This code calculates the phase differences of a user specified 
%%% light pulse at a chosen lux and plots the resultant PRC. The simulation
%%% runs in LD for 25 days, followed by free running in DD for 50 days, and
%%% then the signal is given. The first few days after the signal are
%%% ignored. 

%%% Specify the lux and the light signal for the simulation. E.g. [0 2] 
%%% means a single light pulse signal of 2 hours in length. Note this 
%%% simulation correlates the circadian time with the start of the light
%%% pulse, and not the midpoint of the light pulse. 

lux = 250;
tvec = [0 2]';

param = getParam();

%%% LD
RT = 24;
numdays = 25;
tf = RT*numdays;
numpts=100*numdays;
numconts=64;
lightmax = 150;
dayperiod = 24;
daylength = 12;
light = @(t) lightmax*(mod(t,dayperiod)< daylength);
v0 = [0.0055934725745285;-0.1940044085122700;1.0437461312032239];
options = odeset('Events',@(t,v) peakevent(t,v,param, light));
[t,v,te,ve,ie] = ode23s( @(t,v) rhsKronauer(t,v,param,light), linspace(0,tf,numpts), v0, options);

%%% DD
light = @(t) 0;
v0 = [v(end,1), v(end,2), v(end,3)];
numdays = 50;
tf = param.taux*numdays;
numpts=100*numdays;
options = odeset('Events',@(t,v) peakevent(t,v,param, light), 'Reltol',1e-10,'AbsTol',1e-11);
[t,v_DD,te,ve,ie] = ode45( @(t,v) rhsKronauer(t,v,param,light), linspace(0,tf,numpts), v0, options);
tignore = param.taux*(numdays - 14);

%%% Initial phases 
CBT_initvec = (mod(te(te>tignore), param.taux));
CBT_init = mean(CBT_initvec);

%%% Construct vector of CTs
CTvec = linspace(0,24,250);
CTvec = CTvec*param.taux/24;

%%% Calculate the phase differences
phdiff = zeros(1,length(CTvec));
parfor i = 1:length(CTvec)
    phdiff(i) = pulse(CTvec(i)+tvec,v_DD, param, CBT_init, lux );
end
    
%%% Plotting 
figure; hold on
box on
xaxis = CTvec*24/param.taux;
plot(xaxis, phdiff, 'Color', 'k', 'LineWidth', 1);
set(gca, 'XTick', 0:6:24);
axis([0 24 -4 2])
xlabel('CT (h)', 'FontName', 'SansSerif');
ylabel('Phase difference (h)', 'FontName', 'SansSerif');

function phdiff = pulse(tvec, v_DD, param, initCBT, lux)

    numdays = 11;
    tf = param.taux*numdays;
    tignore = param.taux*(numdays - 7);
    v0 = [v_DD(end,1), v_DD(end,2), v_DD(end,3)];

    light = @(t) reshape(lux*( mod(sum(t(:) > tvec(:)',2),2)==1 ),size(t)); 

    options = odeset('Events',@(t,v) peakevent(t,v,param, light),'Reltol',1e-8,'AbsTol',1e-9);
    tcur = 0;
    tt = [];
    tevec = [];
    vvec = [];
    while tcur<tf
        tfcur = min([tvec(tvec>tcur);tf]);
        [t,v,te,ve,ie] = ode45( @(t,v) rhsKronauer(t,v,param,light), linspace(tcur,tfcur,250), v0, options);
        v0 = v(end,:);
        tevec = [tevec;te];
        tt = [tt; t];
        tcur = tfcur;
        vvec = [vvec; v];
    end

    CBT_finvec = (mod(tevec(tevec>tignore), param.taux));
    CBT_fin = mean(CBT_finvec);
    phdiff = -mean(CBT_fin - initCBT) * 24/param.taux;

end

% the right-hand side of the Kronauer model
function dvdt = rhsKronauer(t,v,param,light)
    I = light(t);
    n = v(1); C = v(2); Cc = v(3); 
    dvdt = NaN(size(v));
    
    alpha = param.alpha0*(I/param.I0)^param.p; 
    B = param.G*(1-n)*alpha*(1-param.b*(C))*(1-param.b*(Cc));
    dvdt(1) = param.lambda*( alpha*(1-n) - param.beta*n );
    dvdt(2) = (pi/12)*(Cc + param.gamma*(C/3 + 4/3*C^3 - 256/105*C^7) + B);
    dvdt(3) = (pi/12)*(B*Cc/3- C*((24/(param.taux*param.f))^2 + param.k*B));

end

function param = getParam
    param.taux = 23.6;
    param.I0 = 9500;
    param.lambda = 60;
    param.f = 0.9974145;
    param.gamma = 0.13;
    param.beta = 0.005;
    param.alpha0 = 1.8;
    param.G = 52;
    param.k =  0.20;
    param.b = 0.59;
    param.p = 0.64;
end

function [value,isterminal,direction] = peakevent(t,v,param,light)
    dvdt = rhsKronauer(t,v,param,light);
    value = dvdt(2);    % the value that we want to be zero
    isterminal = 0;     % halt integration 
    direction = 1;      % -1, dvdt decreasing means max
end