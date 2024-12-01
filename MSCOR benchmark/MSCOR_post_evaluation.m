values = [7,2,3,6]; % 2 = Griewank, 3 = Negative Sum Squares, 6 = Rastrigin, 7 = Ackley.
B = 5;
M = 5;
NumExp = 100;

which_function = values(1);
for funIdx = 1:4
    which_function = values(funIdx);
    filename = sprintf('Summary_funVals_%d_B_%d_M_%d_NumExp_%d.csv', which_function,B,M,NumExp);
    FunVals = readmatrix(filename);
    filename2 = sprintf('Summary_times_%d_B_%d_M_%d_NumExp_%d.csv', which_function,B,M,NumExp);
    Times = readmatrix(filename2);
    if(funIdx == 1)
        Table = [min(abs(FunVals))', (std(abs(FunVals))/sqrt(length(FunVals)))',...
            round(mean(Times)',2), round((std(abs(Times))/sqrt(length(Times)))',3)];
    else
        Table = [Table; [min(abs(FunVals))', (std(abs(FunVals))/sqrt(length(FunVals)))',...
            round(mean(Times)',2), round((std(abs(Times))/sqrt(length(Times)))',3)]];
    end
end

filename = sprintf('Summary_FINAL_B_%d_M_%d.csv',B,M);
writematrix(Table, filename);


%%%% Figure (only for B = 10, M= 20 scenario) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(B == 10 && M == 20)
    
    figure;
    t = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    for funIdx = 1:4
        which_function = values(funIdx);
        filename = sprintf('Summary_funVals_%d_B_%d_M_%d_NumExp_%d.csv', which_function,B,M,NumExp);
        FunVals = readmatrix(filename);
        
        nexttile;
        boxplot(log10(abs(FunVals)), 'Labels',{'MSCOR','GA','SA','IP','SQP','AS'})
        set(gca, 'FontSize', 13,'FontWeight', 'bold');
        ylabel('log_{10}(f_{soln})', 'Interpreter','tex', 'FontWeight', 'bold', 'FontSize', 14)
        if(which_function == 2)
            title('Griewank (modified)','FontWeight', 'bold', 'FontSize', 20)
        end
        if(which_function == 3)
            title('Neg. sum. square (modified)','FontWeight', 'bold', 'FontSize', 20)
        end
        if(which_function == 6)
            title('Rastrigin (modified)','FontWeight', 'bold', 'FontSize', 20)
        end
        if(which_function == 7)
            title('Ackley (modified)','FontWeight', 'bold', 'FontSize', 20)
        end
        
    end
    
    % Set figure size (in inches) and DPI
    factor = 1.2;
    width = 10*factor;  % Width in inches
    height = 5*factor; % Height in inches
    dpi = 300;  % DPI (dots per inch)
    
    % Adjust the figure size
    set(gcf, 'Units', 'inches', 'Position', [0, 0, width, height]);
    
    % Save the figure with the specified DPI
    print(gcf, 'Benchmark_comparison', '-dpng', ['-r' num2str(dpi)]);
    hold off
end





