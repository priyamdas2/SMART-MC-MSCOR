function Plot_DummySequences_20Treatments_7_trts(DummySequences,XDummyPatient_Raw)

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
    treatment_order = {'BcD', 'GA', 'IB','DF','Nat','S1P','AL'};
    treatment_colors = [
        1 0 0;    % BcD (Red)
        0 1 1;    % GA (Cyan)
        0 0 1;    % IB (Blue)
        1 1 0;    % DF (Yellow)
        0 0 0;    % Nat (Black)
        1 0.5 0;  % S1P (Orange)
        1 0 1;    % AL (Magenta)
        ];
    colormap(treatment_colors);
    c = colorbar;
    caxis([1 7]);
    set(c, 'YDir', 'reverse');
    c.Ticks = 1.45:.9:6.85;  % Positions the ticks at each color level
    c.TickLabels = treatment_order;
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

