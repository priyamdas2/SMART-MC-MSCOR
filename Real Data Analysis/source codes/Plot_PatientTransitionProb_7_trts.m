
function Plot_PatientTransitionProb_7_trts(MFinalCell, XRaw, PatientNumbers)
TreatmentNumbers = [1 2 3 4 5 6 7];
% -- TreatmentNamesFull = {'B-cell depletion', 'glatiramer acetate', 'Interferon-beta',...
%     'dimethyl fumarate','natalizumab', 'S1P modulators (fingolimod + teriflunomide)',...
%     'Aggressive / Legacy therapies (cyclophosphamide, mitoxantrone, alemtuzumab)'};
% -- Covariates HERE: c("age_at_diag", "disease_duration", "gender","white", "black")
TreatmentNamesShort = {'BcD', 'GA', 'IB','DF','Nat','S1P','AL'};

NumPatients = length(PatientNumbers);

for k = 1:NumPatients
    
    % Gender conversion
    Covariates = XRaw(PatientNumbers(k),2:end);
    if(Covariates(3) == 1)
        gender = "F";
    else
        gender = "M";
    end
    
    
    figure('Position', [100, 100, 800, 800]);
    patNum = PatientNumbers(k);
    
    Mat2disp = MFinalCell(:,:,patNum);
    Mat2dispISV = Mat2disp(1,:);
    Mat2dispTM = Mat2disp(2:end,:);
    % Subplot 1: Initial State Vector
    subplot(2, 1, 1);
    set(gca, 'XAxisLocation', 'top');
    imagesc(Mat2dispISV);
    colormap(hot);
    caxis([0 1]);
    xlabel('TO-Treatments','FontWeight', 'bold', 'FontSize', 16);
    ylabel('');
    % ylabel('Initial SV','FontWeight', 'bold','FontSize', 14);
    xticks(TreatmentNumbers);
    xticklabels(TreatmentNamesShort);
    title(sprintf('Patient number %d', patNum),'FontWeight', 'bold', 'FontSize', 28);
    s = subtitle(sprintf('(Age = %d y, Dis. dur. = %d m, Sex = %s, White = %d, Black = %d)', ...
        Covariates(1), Covariates(2), gender, Covariates(4),Covariates(5)),...
        'FontSize', 16);
    s.Color = 'blue'; 
    set(gca, 'YTick', []);
    
    ax1 = gca; % Get current axes handle
    ax1.Position(4) = ax1.Position(4) * 0.15;
    hold on; % Keep the current plot
    for i = 1:size(Mat2dispISV, 1)
        for j = 1:size(Mat2dispISV, 2)
            % Draw a rectangle for each tile
            rectangle('Position', [j-1+0.5, i-1+0.5, 2, 2], 'EdgeColor', 'w', 'LineWidth', 2);
        end
    end
    hold off;
    
    % Subplot 2: Transition matrix
    subplot(2, 1, 2);
    imagesc(Mat2dispTM);
    colormap(hot);
    caxis([0 1]);
    colorbar;
    ylabel('FROM-Treatments','FontWeight', 'bold', 'FontSize', 16);
    xticks(TreatmentNumbers);
    xticklabels(TreatmentNamesShort);
    yticks(TreatmentNumbers);
    yticklabels(TreatmentNamesShort);
    set(gca, 'XAxisLocation', 'top');
    
    hold on; % Keep the current plot
    for i = 1:size(Mat2dispTM, 1)
        for j = 1:size(Mat2dispTM, 2)
            % Draw a rectangle for each tile
            rectangle('Position', [j-1+0.5, i-1+0.5, 2, 2], 'EdgeColor', 'w', 'LineWidth', 2);
        end
    end
    hold off;
    print(['Plots/Patient_Transition_Mat_SerialNo_',num2str(patNum)], '-dpng', '-r300')
    
end

end
