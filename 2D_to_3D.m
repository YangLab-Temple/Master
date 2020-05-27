function [] = Sasha()

    % parameters
    global file working_dir filename smooth_window d;
    % indicate the name of the file with your data here
    file = 'nup 58 1100 points n term.csv';
    % indicate the name of the folder where you want your results to be
    % saved (create a folder and specify the path in MatLab)
    working_dir = 'nup58 N term';
    filename = [working_dir filesep file];
    % here specify the level of smoothing
    smooth_window = 2;
    % 
    d = 360;
    ranges = { [50 90] };
    for j = 1:length(ranges)
        Chi = [];
       for i = 7 %dr loop, here indicate bin size you want to use. If bin size e.g 3 set it to i = 5
            Chi(end+1) = 1 - MainCalculations(i,ranges{j});
       end
        [best,index] = min(Chi);
        fprintf('Range %i %i best Chi %f with dr %i\n', ranges{j}(1),ranges{j}(2),1-Chi(index),index); 
    end
end

function [Chi] = MainCalculations(dr,range)

    global file;
    global working_dir;
    global filename;
    global smooth_window;
    global d;

    R = (floor(100/dr))*dr;
    bin_size = dr;
    % open csv file
    M = csvread(filename);
    % sort by col 1
    M = sortrows(M,1);
    % choose interval in col 1, take corresponding values from col 2
    range = sort(range);
    a = find(M(:,1)>range(1));
    b = find(M(:,1)<range(2));
    X = intersect(a,b);
    Y = M(X,2);
    % histogram count
    edges = -R : bin_size : R; % -100 : 5 : +100
    bin_center = -R+bin_size/2 : bin_size : R-bin_size/2; % e.g -97.5 : 5 : 97.5
    n_elements = histc(Y,edges);
    n_elements(end-1) = n_elements(end-1) + n_elements(end); % the last element in n_elements
    n_elements(end) = []; % contains number of elements equal to the maximum bin_center
    % smooth
    smoothed_n_elements = smooth(n_elements,smooth_window);

    flip = flipud(smoothed_n_elements);
    for i = 1:length(flip(:,1))
        cum(i,1) = flip(i,1) + smoothed_n_elements(i,1);
    end

    Z = [reshape(bin_center,[],1) cum];
    Z(Z(:,1)<0,:) = [];

    sm = zeros(R/dr + 1,1);
    for i = 1:length(Z)
        sm(i,1) = Z(i,2);
    end
    % here you matrix will be generated depending on your bin size and dr
    Matrix = generateMatrix(R,dr,d);
    mr = inv(Matrix) * sm;
    msm = smooth(mr,smooth_window);
    file_id = sprintf('R_%i_dr_%i_range_%i-%i',R,dr,range(1),range(2));
    bar(Z(:,1),msm(1:length(Z)));
    % the following commands will display the title on each graph generated
    title (['R=',num2str(R),blanks(10),'dr=',num2str(dr),blanks(10),num2str(range(1)),' to ', num2str(range(2))],'FontSize',32);
    saveas(gcf,[working_dir filesep file_id '_msm.png']);
    bar(Z(:,1),sm(1:length(Z)));
    title (['R=',num2str(R),blanks(10),'dr=',num2str(dr),blanks(10),num2str(range(1)),' to ', num2str(range(2))],'FontSize',32)
    saveas(gcf,[working_dir filesep file_id '_sm.png']);
    bar(Z(:,1),mr(1:length(Z)));
    title (['R=',num2str(R),blanks(10),'dr=',num2str(dr),blanks(10),num2str(range(1)),' to ', num2str(range(2))],'FontSize',32)
    saveas(gcf,[working_dir filesep file_id '_m.png']);
    close(gcf);
    csvwrite([working_dir filesep file_id '.csv'],[sm mr msm]);

   % if ~isempty(find(msm<0,1))
     %   disp('Negative values detected');
        msm(msm<0) = 0;
        reverse = Matrix * msm;
        nominator = 0;
        for i = 1:length(reverse)
            nominator = nominator + (sm(i) - reverse(i))^2;
        end
       Chi = 1 - nominator/(sum(sm(:)))^2;
       %Chi = 1 - nominator/(sum(sm(:).^2));
       %Chi = 1 - nominator/sum(power(sm,2));
        fid = fopen([working_dir filesep file_id '_Chi_' num2str(Chi) '.txt'],'w');
        fprintf(fid,'%f',Chi);
        fclose(fid);
        fprintf('Chi square: %f\n',Chi);
   % end

%    keyboard;

end

function [S] = smooth(Data, window)

    if mod(window,2) == 0
        range = window/2;
        denominator = window + 1;
    else
        range = (window-1)/2;
        denominator = window;
    end
    S = Data;

    for i = 1+range : length(Data(:,1))-range
        S(i,1) = sum(sum(Data(i-range:i+range,1))) / denominator;
    end

end

function [M] = generateMatrix(R,dr,d)

    if mod(R,dr) ~= 0
        error('R/dr must be integer');
    end

    N = R/dr + 1;
    M = zeros(N,N);

    for i = 1:N-1
        for j = 1:N-1
            k = i - 0.5;
            if i < j
                M(i,j) = sqrt(j^2 - k^2) - sqrt((j-1)^2 - k^2);
            elseif i == j
                M(i,j) = sqrt(j^2 - k^2);
            end

        end
    end
    for i = 1:N
        k = i - 0.5;
        M(i,N) = real((d/2)/dr - sqrt((N-1)^2 - k^2));
    end
%     disp(M);

end
