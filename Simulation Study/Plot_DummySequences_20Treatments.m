function Plot_DummySequences_20Treatments(DummySequences,XDummyPatient_Raw)

NumDummyClasses = size(DummySequences,3);

figure; % imp treatments 1,3,4,5,6
t = tiledlayout(3, 4, 'Padding', 'compact', 'TileSpacing', 'compact');

for DummyPatientNum = 1:NumDummyClasses
    
    nexttile;
    %%%% Extracting patient characteristics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    age = XDummyPatient_Raw(DummyPatientNum, 2);
    
    if(XDummyPatient_Raw(DummyPatientNum, 4) == 1)
        gender = 'F';
    else
        gender = 'M';
    end
    
    if(XDummyPatient_Raw(DummyPatientNum, 5) == 1)
        race = 'White';
    else
        if(XDummyPatient_Raw(DummyPatientNum, 6) == 1)
            race = 'Black';
        else
            race = 'Others';
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    Mat2bePlotted = DummySequences(:,:,DummyPatientNum);
    imagesc(Mat2bePlotted);
    customColormap = [
        1 0 0;    % Red
        0 1 1;    % Cyan
        0 0 1;    % Blue
        1 1 0;    % Yellow
        0 0 0;    % Black
        1 0.5 0;  % Orange
        1 0 1;    % Magenta
        0.5 0.5 0.5;  % Gray
        0.5 0 0.5;  % Purple
        0.5 1 0;  % Light Green
        ];
    colormap(customColormap);
    c = colorbar;
    caxis([1 10]);
    set(c, 'YDir', 'reverse');
    c.Ticks = 1.45:.9:9.55;  % Positions the ticks at each color level
    c.TickLabels = {'BcD', 'GA', 'IB','DF','Nat','Fin','Ter','Cyc','Mit','Ale'};
    xlabel('Doctor visits', 'FontWeight', 'bold', 'FontSize', 14);
    ylabel('Patients', 'FontWeight', 'bold', 'FontSize', 14);
    title([num2str(age),' yrs, ',gender,', ',race], 'FontWeight', 'bold', 'FontSize', 18);
    
    hold on; % Keep the current plot
    for i = 1:size(Mat2bePlotted, 1)
        for j = 1:size(Mat2bePlotted, 2)
            % Draw a rectangle for each tile
            rectangle('Position', [j-1+0.5, i-1+0.5, 2, 2], 'EdgeColor', 'w', 'LineWidth', 0.5);
        end
    end
   
    
    %print(['Plots/Generated_patient_seq_',num2str(DummyPatientNum)], '-dpng', '-r300')
end

% Set figure size (in inches) and DPI
factor = 2;
width = 8*factor;  % Width in inches
height = 5*factor; % Height in inches
dpi = 300;  % DPI (dots per inch)

% Adjust the figure size
set(gcf, 'Units', 'inches', 'Position', [0, 0, width, height]);

% Save the figure with the specified DPI
print(gcf, 'Plots/Dummy_sequences', '-dpng', ['-r' num2str(dpi)]);

 hold off;
end

