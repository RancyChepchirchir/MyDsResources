figure(1)
model = model1;

Model_name  = model.name;
Model_param = model.param;

model  = Model(Model_name,Model_param);
Noptions = size(Data,1);

Pmodel = NaN(Noptions,1);
for i = 1:Noptions
    K  = Data(i,2);
    T  = Data(i,3);
    r  = Data(i,4);
    q  = Data(i,5);
    Opt_type = Data(i,6);
    
    method = Method('Conv',[1, 1000, -Opt_type]);
    gain_fun = @(S,K) max(Opt_type*(S-K),0);
    
    option = Option(S0,r,q,K,T,1,gain_fun,@(S,K) 0.*S,0,1000,'F',1); 

    Pmodel(i) = Pricer(option,model,method);
end
%%
P = Data(:,1);
K = Data(:,2);
T = Data(:,3);
scatter3(K/S0,T,P,'r','marker','o')
hold on
scatter3(K/S0,T,Pmodel,'b','marker','*')
view(5,20)
xlabel('\fontsize{12} Moneyness K/S')
ylabel('\fontsize{12} Maturities T')
zlabel('\fontsize{12} Prices')
title('\fontsize{14} Balck-Scholes prices')

%% Jump models
figure(2)
model = model2;
Noptions = size(Data,1);

Pmodel = NaN(Noptions,1);
for i = 1:Noptions
    K  = Data(i,2);
    T  = Data(i,3);
    r  = Data(i,4);
    q  = Data(i,5);
    Opt_type = Data(i,6);
    
    method = Method('Conv',[1, 1000, -Opt_type]);
    gain_fun = @(S,K) max(Opt_type*(S-K),0);
    
    option = Option(S0,r,q,K,T,1,gain_fun,@(S,K) 0.*S,0,1000,'F',1); 

    Pmodel(i) = Pricer(option,model,method);
end
%%
P = Data(:,1);
K = Data(:,2);
T = Data(:,3);
%subplot(2,2,1)
scatter3(K/S0,T,P,'r','marker','o')
hold on
scatter3(K/S0,T,Pmodel,'b','marker','*')
view(5,20)
xlabel('\fontsize{12} Moneyness K/S')
ylabel('\fontsize{12} Maturities T')
zlabel('\fontsize{12} Prices')
title('\fontsize{14} Merton Calibration')

%% ------------
model = model3;

Model_name  = model.name;
Model_param = model.param;

model  = Model(Model_name,Model_param);
Noptions = size(Data,1);

Pmodel = NaN(Noptions,1);
for i = 1:Noptions
    K  = Data(i,2);
    T  = Data(i,3);
    r  = Data(i,4);
    q  = Data(i,5);
    Opt_type = Data(i,6);
    
    method = Method('Conv',[1, 1000, -Opt_type]);
    gain_fun = @(S,K) max(Opt_type*(S-K),0);
    
    option = Option(S0,r,q,K,T,1,gain_fun,@(S,K) 0.*S,0,1000,'F',1); 

    Pmodel(i) = Pricer(option,model,method);
end
%%
P = Data(:,1);
K = Data(:,2);
T = Data(:,3);
subplot(2,2,2)
scatter3(K/S0,T,P,'r','marker','o')
hold on
scatter3(K/S0,T,Pmodel,'b','marker','*')
view(5,20)
xlabel('\fontsize{12} Moneyness K/S')
ylabel('\fontsize{12} Maturities T')
zlabel('\fontsize{12} Prices')
title('\fontsize{14} Kou Calibration')
%% -------------
model = model4;

Model_name  = model.name;
Model_param = model.param;

model  = Model(Model_name,Model_param);
Noptions = size(Data,1);

Pmodel = NaN(Noptions,1);
for i = 1:Noptions
    K  = Data(i,2);
    T  = Data(i,3);
    r  = Data(i,4);
    q  = Data(i,5);
    Opt_type = Data(i,6);
    
    method = Method('Conv',[1, 1000, -Opt_type]);
    gain_fun = @(S,K) max(Opt_type*(S-K),0);
    
    option = Option(S0,r,q,K,T,1,gain_fun,@(S,K) 0.*S,0,1000,'F',1); 

    Pmodel(i) = Pricer(option,model,method);
end
%%
P = Data(:,1);
K = Data(:,2);
T = Data(:,3);
subplot(2,2,3)
scatter3(K/S0,T,P,'r','marker','o')
hold on
scatter3(K/S0,T,Pmodel,'b','marker','*')
view(5,20)
xlabel('\fontsize{12} Moneyness K/S')
ylabel('\fontsize{12} Maturities T')
zlabel('\fontsize{12} Prices')
title('\fontsize{14} NIG Calibration')
%% ----------
model = model5;

Model_name  = model.name;
Model_param = model.param;

model  = Model(Model_name,Model_param);
Noptions = size(Data,1);

Pmodel = NaN(Noptions,1);
for i = 1:Noptions
    K  = Data(i,2);
    T  = Data(i,3);
    r  = Data(i,4);
    q  = Data(i,5);
    Opt_type = Data(i,6);
    
    method = Method('Conv',[1, 1000, -Opt_type]);
    gain_fun = @(S,K) max(Opt_type*(S-K),0);
    
    option = Option(S0,r,q,K,T,1,gain_fun,@(S,K) 0.*S,0,1000,'F',1); 
    Pmodel(i) = Pricer(option,model,method);
end
%%
P = Data(:,1);
K = Data(:,2);
T = Data(:,3);
T_unique = unique(T);
subplot(2,2,4)
scatter3(K/S0,T,P,'r','marker','o')
hold on
scatter3(K/S0,T,Pmodel,'b','marker','*')
view(5,20)
xlabel('\fontsize{12} Moneyness K/S')
ylabel('\fontsize{12} Maturities T')
zlabel('\fontsize{12} Prices')
title('\fontsize{14} VG Calibration')

        