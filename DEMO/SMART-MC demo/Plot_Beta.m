% BetaIndex : 2 = age_at_diag, 3 = disease_duration, 4 = Sex",
%             5 = white, 6 = black.
function Plot_Beta(BetaIndex,BetaCellFinal,NonNullPositions)
N = size(BetaCellFinal,2);
HeatData = nan((N+1),N);
for i = 1:(N+1)
    for j = 1:N
        if (NonNullPositions(i,j) == 1)
            HeatData(i,j) = BetaCellFinal{i,j}(BetaIndex);
        end
    end
end

HeatDataISV = HeatData(1,:);
HeatDataTM = HeatData(2:end,:);


figure;
subplot(2, 1, 1);
Data = round(HeatDataISV, 2);
h = heatmap(Data);
h.Position=[0.130 0.800 0.7179 0.12]; % [X,Y, width, length]
if (BetaIndex == 2)
    h.Title ='Coeff. of Age (transformed)';
else
    if (BetaIndex == 3)
        h.Title = 'Coeff. of Disease duration (transformed)';
    else
        if (BetaIndex == 4)
            h.Title = 'Coeff. of Sex (female = 1)';
        else
            if (BetaIndex == 5)
                h.Title = 'Coeff. of Race White (white = 1, black/others =0)';
            else
                if (BetaIndex == 6)
                    h.Title = 'Coeff. of Race Black (black = 1, white/others =0)';
                end
            end
        end
    end
end
mymap = [
    1 0 0
    0 0 1];
colormap(h, mymap);
h.XDisplayLabels = {'BcD', 'GA', 'IB','DF','Nat' 'Fin','Ter','Cyc',...
    'Mit','Ale'};
% h.YDisplayLabels = repmat({''}, size(Data, 2), 1);
h.YDisplayLabels = {''};
h.ColorLimits = [-1, 1];
h.XLabel = 'TO-Treatment';          % Set the X-axis label
h.ColorbarVisible = 'off';
h.FontSize = 12;                % Controls all fonts





subplot(2, 1, 2);
set(gca, 'XAxisLocation', 'top');
Data = round(HeatDataTM, 2);
h = heatmap(Data);
h.Position=[0.130 0.100 0.7179 0.55]; % [X,Y, width, length]
mymap = [
    1 0 0
    0 0 1];
colormap(h, mymap);
h.XDisplayLabels = {'BcD', 'GA', 'IB','DF','Nat' 'Fin','Ter','Cyc',...
    'Mit','Ale'};
h.YDisplayLabels = {'BcD', 'GA', 'IB','DF','Nat' 'Fin','Ter','Cyc',...
    'Mit','Ale'};
h.ColorLimits = [-1, 1];
h.YLabel = 'FROM-Treatment';          % Set the Y-axis label
h.ColorbarVisible = 'on';           % Show colorbar
h.FontSize = 12;
print(['Plots/Beta_values_',num2str(BetaIndex)], '-dpng', '-r300')
end