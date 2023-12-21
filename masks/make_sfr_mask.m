%% make Schaefer mask



mask = 'Schaefer2018_400Parcels_7Networks_order_FSLMNI152_2mm.nii'
nii = load_nii(mask);


tbl = readtable('Schaefer2018_400Parcels_7Networks_order_FSLMNI152_2mm.Centroid_RAS.csv');


% whole networks
wholeNames = {...
    'Sfr_Vis',...
    'Sfr_DorsAttn',...
    'Sfr_SalVentAttn',...
    'Sfr_Cont',...
    'Sfr_Default',...
    }


% sub networks
subNames = {
    ...
    'Sfr_Cont',...
    'Sfr_SalVentAttn',...
    }

subParcels = {
    ...
    [148, 360, 361, 143, 340],...
    [107, 108, 110, 311, 312, 314, 99, 101, 303, 306],...
    }




%% overall networks
for ff = 1:length(wholeNames)


    parcel = strsplit(wholeNames{ff}, {'_', '.'});
    rois   = tbl.ROILabel(cellfun(@any, regexpi(tbl.ROIName, parcel{2})));

    parcelNii = nii;
    parcelNii.img = ismember(parcelNii.img,  rois);

    save_nii(parcelNii, [wholeNames{ff}, '.nii']);

end


%% sub-networks
for ff = 1:length(subNames)

    for ss = 1:length(subParcels{ff})

        parcelNii = nii;
        parcelNii.img = ismember(parcelNii.img,  subParcels{ff}(ss));

        save_nii(parcelNii, [subNames{ff}, '_', num2str(subParcels{ff}(ss)), '.nii']);
    end

end



