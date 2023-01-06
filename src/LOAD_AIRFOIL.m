function [dataX,dataY] = LOAD_AIRFOIL(fileName)

fidAirfoil = fopen(['./Airfoil_DAT/' fileName '.dat']);
dataCell = textscan(fidAirfoil,'%s','Delimiter','\n');
fclose(fidAirfoil);
hdrlns = nan;
for k = 1:1:length(dataCell{1})
    charArr = cell2mat(dataCell{1,1}(k));
    try
        C = strsplit(charArr,' ');
        if (str2double(C(1)) <= 1 && str2double(C(1)) >= 0)
            fprintf('\tHeader Lines: %i\n',k-1);
            hdrlns = k-1;
            break;
        end
    catch
    end
    if (isnan(hdrlns))
        C = strsplit(charArr);
        if (str2double(C(1)) <= 1 && str2double(C(1)) >= 0)
            fprintf('\tHeader Lines: %i\n',k-1);
            hdrlns = k-1;
            break;
        end
    end
end

firstLine = dataCell{1}(hdrlns+1);
split = strsplit(cell2mat(firstLine),' ');

if (length(split) == 1)
    dataX = zeros(length(dataCell{1})-hdrlns,1);
    dataY = zeros(length(dataCell{1})-hdrlns,1);
    
    for i = 1:1:length(dataCell{1})-hdrlns
        if (length(dataCell{1}{i+hdrlns}) == 3)
            vals = cell2mat(dataCell{1}(i+hdrlns));
            dataX(i,1) = str2double(vals(1));
            dataY(i,1) = str2double(vals(3));
        elseif (length(dataCell{1}{i+hdrlns}) == 0)
            
        else
            vals = strsplit(dataCell{1}{i+hdrlns});
            dataX(i,1) = str2double(vals(1));
            dataY(i,1) = str2double(vals(2));
        end
    end
else
    fidAirfoil = fopen([',/Airfoil_DAT/' fileName '.dat']);
    dataBuffer = textscan(fidAirfoil,'%f %f','CollectOutput',1,...
                                                            'Delimiter','',...
                                                            'Whitespace',' ',...
                                                            'HeaderLines',hdrlns);
                                                    
    dataX = dataBuffer{1}(:,1);
    dataY = dataBuffer{1}(:,2);
    fclose(fidAirfoil);
end

if (strcmpi('bacj',fileName))
    dataX = dataX(1:end-1);
    dataY = dataY(1:end-1);
end

indMax = find(dataX == max(dataX));
if (length(indMax) < 2)
    dataX2 = dataX;
    dataX2(indMax(1)) = 0;
    secInd = find(dataX2 == max(dataX2));
    if (min(indMax > secInd) == 1)
        indMax(2) = secInd(1);
    else
        indMax(2) = secInd(secInd > indMax(1));
    end
end

if (strcmp(fileName,'ste87151') || strcmp(fileName,'ste87391') || strcmp(fileName,'stf86361'))
    indMax = [indMax(1); length(dataX)];
end

indMax = sort(indMax);

if (indMax(1) == 1 && indMax(2) == length(dataX))
    dataXX = dataX;
    dataYY = dataY;
elseif (indMax(1) ==  1 && indMax(2) ~= length(dataX))
    dataXX = [dataX(1:indMax(2)-1); flipud(dataX(indMax(2):end))];
    dataYY = [dataY(1:indMax(2)-1); flipud(dataY(indMax(2):end))];
elseif (indMax(1) ~= 1 && indMax(2) == length(dataX))
    dataXX = [flipud(dataX(1:indMax(1))); dataX(indMax(1)+1:end)];
    dataYY = [flipud(dataY(1:indMax(1))); dataY(indMax(1)+1:end)];
elseif (indMax(1) ~= 1 && indMax(2) ~= length(dataX))
    dataXX = [dataX(indMax(2):end); dataX(1:indMax(1))];
    dataYY = [dataY(indMax(2):end); dataY(1:indMax(1))];
end
dataX = dataXX;
dataY = dataYY;
dataArr = [dataX dataY];

[~,ia,~] = unique(dataArr,'rows','stable');
i = true(size(dataArr,1),1);
i(ia) = false;
dataArr(i,:) = [];
dataX = dataArr(:,1);
dataY = dataArr(:,2);

if (dataY(1) ~= dataY(end))
    dataX(end+1) = dataX(1);
    dataY(end+1) = dataY(1);
end
    
    