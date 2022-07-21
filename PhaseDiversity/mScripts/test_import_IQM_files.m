fileFolderIn = '..\..\..\computAO\Data\20200820_Beads_PR\Progressive_batch_50\';
fileNameBase1 = 'beads_IQM_InitialAber_';
fileNameBase2 = 'beads_IQM_Corrected_';
raw1 = zeros(1,2,50);
raw2 = zeros(4,2,50);
iqmGT = zeros(50,1);
iqmAbe = zeros(50,1);
iqmCorr = zeros(50,4);
iqmCorrBetter = zeros(50,4);
iqmCorr2step = zeros(50,1);
iqmCorr2stepBetter = zeros(50,1);
for i = 1:50
    a1 = importdata([fileFolderIn, fileNameBase1, num2str(i-1), '.txt']);
    a2 = importdata([fileFolderIn, fileNameBase2, num2str(i-1), '.txt']);
    raw1(:,:,i) = a1;
    raw2(:,:,i) = a2;
    iqmGT(i) = a1(1);
    iqmAbe(i) = a1(2);
    iqmCorr(i,:) = a2(:,2);
    % count the better cases
    if(a2(1,2)>a1(2))
        iqmCorrBetter(i,1) = 1;
        iqmCorr2step(i) = a2(1,2);
    else
        iqmCorr2step(i) = a2(2,2);
    end
    if(a2(2,2)>a1(2))
        iqmCorrBetter(i,2) = 1;
    end
    if(a2(3,2)>a1(2))
        iqmCorrBetter(i,3) = 1;
    end
    if(a2(4,2)>a1(2))
        iqmCorrBetter(i,4) = 1;
    end
    if(iqmCorr2step(i)>a1(2))
        iqmCorr2stepBetter(i) = 1;
    end
end

iqmGT_Ave = mean(iqmGT(:));
iqmGT_SD = sqrt(var(iqmGT(:)));
iqmAbe_Ave = mean(iqmAbe(:));
iqmAbe_SD = sqrt(var(iqmAbe(:)));
iqmCorr_Ave = mean(iqmCorr, 1);
iqmCorr_SD = sqrt(var(iqmCorr,1));
iqmCorr2step_Ave = mean(iqmCorr2step(:));
iqmCorr2step_SD = sqrt(var(iqmCorr2step(:)));
iqmCorrBetter_Total = sum(iqmCorrBetter,1);
iqmCorr2stepBetter_Total = sum(iqmCorr2stepBetter(:));