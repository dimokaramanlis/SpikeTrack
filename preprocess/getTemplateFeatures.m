function resfeats = getTemplateFeatures(edata, coords)
%GETTEMPLATEFEATURES Summary of this function goes here
%   Detailed explanation goes here

%--------------------------------------------------------------------------
% copy properties from edata
resfeats = struct();
resfeats.spikecv = mean(edata.spikecv, 2, "omitmissing", "Weights", edata.unitspikes);
spikelv  = edata.spikelv;
spikelv(isinf(spikelv)) = 0;
resfeats.spikelv = mean(spikelv, 2, "omitmissing", "Weights", edata.unitspikes);
resfeats.t2ptime = fixWaveformTimes(edata.t2ptime, edata.unitspikes, 1e-3);
resfeats.trepol  = fixWaveformTimes(edata.trepol, edata.unitspikes,  1e-3);
resfeats.t2pval  = fixWaveformTimes(edata.t2pval, edata.unitspikes,  1);
%--------------------------------------------------------------------------
% extract baseline firing rate from ACG
[Nunits, Nsamps, Ndays] = size(edata.autoCorrs);
dacg              = median(diff(edata.autoCorrsTime));
allacg            = edata.autoCorrs./reshape(edata.unitspikes,[Nunits 1 Ndays])/dacg;
bslrates          = squeeze(mean(allacg(:,end-10:end, :), 2));
resfeats.bslfr = mean(bslrates, 2, "omitmissing", "Weights", edata.unitspikes);
%--------------------------------------------------------------------------
tfs = median(diff(edata.templateTimes));
templatemat = edata.stimTemplatesMean;

maxRexp       = 200;
Nsitesprop    = 13;
Ncompsdenoise = 3;
dttrough      = -20:40;
[Nunits, Nchan, Nt, Ndays] = size(templatemat);
medrem              = median(templatemat,[2 3]);
templatemat         = templatemat - medrem;

twts      = squeeze(range(templatemat, 3));
%--------------------------------------------------------------------------
% initialize collections
resfeats.unitcents   = nan(Nunits, 2,  'single');
resfeats.avgwforms   = nan(Nunits, Nt, 'single');
resfeats.wfmdecays   = nan(Nunits, 1,  'single');
resfeats.gausscents  = nan(Nunits, 2,  'single');
resfeats.gausssigma  = nan(Nunits, 1,  'single');
resfeats.tfitprofile = nan(Nunits, 27, 'single');
resfeats.cfitprofile = nan(Nunits, 27, 'single');
resfeats.unitamps    = nan(Nunits, 1,  'single');
resfeats.speedtop    = nan(Nunits, 1,  'single');
resfeats.speedbottom = nan(Nunits, 1,  'single');
resfeats.axonval     = nan(Nunits, 1,  'single');
%--------------------------------------------------------------------------
for iunit = 1:Nunits 
    unitdayspikes     = edata.unitspikes(iunit, :);
    daycents          = nan(Ndays, 2);
    daywvfm           = nan(Ndays, Nt);
    daydecays         = nan(Ndays, 1);
    dayscentsgauss    = nan(Ndays, 2);
    daysigma          = nan(Ndays, 1);
    daytfitprofile    = nan(Ndays, 27);
    daycfitprofile    = nan(Ndays, 27);
    dayamps           = nan(Ndays, 1);
    dayspeedtop       = nan(Ndays, 1);
    dayspeedbot       = nan(Ndays, 1);
    dayaxonvalues     = nan(Ndays, 1);

    for iday = 1:Ndays
        utemp      = squeeze(templatemat(iunit,:,:, iday));
        if all(isnan(utemp), 'all') || all(utemp == 0, "all") || (unitdayspikes(iday) < 10)
            continue;
        end
        [~, imsite]  = max(twts(iunit,:, iday));
        [~, itrough] = min(utemp(imsite,:));
    %     cwts  = max(abs(utemp(:, dtcell)), [], 2);
    %     cwts  = cwts - quantile(cwts,0.05);
        allds = sqrt(sum((coords - coords(imsite,:)).^2,2));
        cwts  = twts(iunit, :, iday)';
        ichankeep = allds<maxRexp & allds>0 & cwts > 0;
        pfit  = polyfit(allds(ichankeep), log(cwts(ichankeep)), 1);
        %---------------------------------------------------------------------
        % extract decay-based center and avg waveform
    
        xmax     = log(10)/abs(pfit(1));
        iptsuse  = allds<xmax;
        currcent = cwts(iptsuse)'*coords(iptsuse,:)/sum(cwts(iptsuse));
    
        alldsuse = sqrt(sum((coords(iptsuse,:) - currcent).^2,2));
        wvfwts   = 1 - alldsuse/xmax;
        wvfwts   = wvfwts.*(wvfwts>0);
        wvfwts   = wvfwts/sum(wvfwts);
        avgwfm   = wvfwts'*utemp(iptsuse, :);
    
        dayamps(  iday)    = twts(iunit, iptsuse, iday) * wvfwts;
        daycents( iday, :) = currcent;
        daywvfm(  iday, :) = avgwfm;
        daydecays(iday)    = mean((cwts(imsite)-cwts(iptsuse))./allds(iptsuse), 'omitnan');
        %---------------------------------------------------------------------
        % extract Gaussian parameters
    
        [~, isortd] = sort(alldsuse,'descend');
        imgfit      = cwts(iptsuse);
        imgfit      = imgfit - mean(imgfit(isortd(1:6)));
        imgfit      = imgfit/max(imgfit);
        fitparams   = fitCircGaussForMEA(coords(iptsuse,:), double(imgfit));
    
        dayscentsgauss(iday, :) = fitparams(1:2);
        daysigma(iday) = fitparams(3);
        fullparams  = [fitparams(1:2) fitparams([3 3]) 0 1];
        cel = getEllipseFromParams(fullparams, 2, 50);
        sitesel = inpolygon(coords(:,1), coords(:,2), cel(1,:), cel(2,:));
        elsiteidx = find(sitesel);
        axvals    = zeros(size(elsiteidx));
        for isite = 1:nnz(sitesel)
            sctemp = utemp(elsiteidx(isite), :);
            [~, iminsc] = min(sctemp);
            normfac = abs(min(sctemp(1:iminsc)));
            if normfac > 0
                axvals(isite) = range(sctemp(1:iminsc))/normfac - 1;
            end
        end
        dayaxonvalues(iday) = cwts(elsiteidx)'*axvals/sum(cwts(elsiteidx));
        %---------------------------------------------------------------------
        % extract spike propagation
        % sitesuse = imsite + (-Nsitesprop:Nsitesprop);
        sitecoords = coords(:, 2) - coords(imsite, 2);
        [sortedsitecoords, isitesort] = sort(sitecoords, 'ascend');
        sitesuse = isitesort(abs(sortedsitecoords) < 100);

        indsave  = 1:numel(sitesuse);
    
    
        timesfit  = NaN(numel(sitesuse), 1);
        valsfit   = zeros(numel(sitesuse), 1);
        wvfsearch = utemp(sitesuse, :);
        wvfsearch = wvfsearch - median(wvfsearch(:, 1:5), 2);
        dtsearch  = itrough + dttrough;
        dtsearch(dtsearch < 1 | dtsearch > Nt) = [];
    
        wvfsmall  = wvfsearch(:, dtsearch);
        
        % svd model of waveform
        [aa, bb, cc] = svd(wvfsmall);
        svdmodel =  aa(:,1:Ncompsdenoise)*bb(1:Ncompsdenoise,1:Ncompsdenoise)*cc(:,1:Ncompsdenoise)';
    %     clf;      
    %     subplot(1,2,1)
    %     plot(wvfsmall' + (1:numel(sitesuse))*10)
    %     subplot(1,2,2)
    %     plot(svdmodel' + (1:numel(sitesuse))*10)
    %     
    
        for ipt = 1:numel(sitesuse)
            wvfmcurr = svdmodel(ipt,:);

            [~,   imincurr] = min(wvfmcurr);
            [~, inegpeaks]  = findpeaks(-wvfmcurr,'NPeaks',3,'SortStr','descend',...
                'MinPeakDistance',6,'MinPeakProminence',0.1);
            
            inegpeaks(inegpeaks < 6 | inegpeaks > numel(wvfmcurr)-6 ) = [];
            if ~isempty(inegpeaks)
                imincurr = inegpeaks(1);
            end
    %         % replace minimum point 
    %         if abs(minval)<posval
    %             iafter      = (inegpeaks - ipospeak)>0;
    %             if sum(iafter)==1
    %                 imincurr =  inegpeaks(iafter);
    %             else
    % %                 [~, ifar] = max(abs((inegpeaks - ipospeak)));
    % %                 iafter = true(size(iafter));
    % %                 iafter(ifar) = 0;
    %             end
    %         end
    
            timesfit(ipt) = imincurr;
            % valsfit(ipt)  = svdmodel(ipt,imincurr);
        end
    %   
        %---------------------------------------------------------------------
        % triage propagation
    %     line(timesfit, valsfit +(1:numel(sitesuse))'*10,...
    %         'LineStyle','none','Color','k','Marker','o')
    
        % cfit    = coords(sitesuse,2) - coords(imsite,2);
        cfit    = sitecoords(sitesuse);
        daycfitprofile(iday, indsave) = cfit;
        fitpoly2=fit(cfit,timesfit,'poly2', 'Weights',range(wvfsmall,2));
    
        tpred = feval(fitpoly2, cfit);
        
        properrs = timesfit - tpred;
        irem     = abs(properrs) > mad(properrs, 1) * 1.4 * 6;
        timesfit(irem) = nan;
    %     line(timesfit, valsfit +(1:numel(sitesuse))'*10,...
    %         'LineStyle','none','Color','r','Marker','o')
        daytfitprofile(iday, indsave) = timesfit;
        %---------------------------------------------------------------------
        [cfit, ~, icun] = unique(cfit);
        tfit = accumarray(icun, timesfit,[],@nanmean);
        
        ispeedtop = 0;
        iuse = cfit>=0 & ~isnan(tfit);
        if nnz(iuse) > 4
            coefftop = polyfit(cfit(iuse), tfit(iuse), 1);
            ispeedtop = abs(coefftop(1) * tfs)*1e6;
        end
        ispeedbot = 0;
        iuse = cfit<=0 & ~isnan(tfit);
        if  nnz(iuse)  > 4
            coeffbot = polyfit(cfit(iuse), tfit(iuse), 1);
            ispeedbot = abs(coeffbot(1) * tfs)*1e6;
        end
        
        dayspeedtop(iday) = ispeedtop;
        dayspeedbot(iday) = ispeedbot;
        %---------------------------------------------------------------------
    end


    resfeats.unitcents(iunit, :)   = mean(daycents, 1, "omitmissing", "Weights", unitdayspikes);
    resfeats.gausscents(iunit, :)  = mean(dayscentsgauss, 1, "omitmissing", "Weights", unitdayspikes);
    resfeats.unitamps(iunit, :)    = mean(dayamps, 1, "omitmissing", "Weights", unitdayspikes);
    resfeats.speedtop(iunit, :)    = mean(dayspeedtop, 1, "omitmissing", "Weights", unitdayspikes);
    resfeats.speedbottom(iunit, :) = mean(dayspeedbot, 1, "omitmissing", "Weights", unitdayspikes);
    resfeats.axonval(iunit, :)     = mean(dayaxonvalues, 1, "omitmissing", "Weights", unitdayspikes);
    resfeats.wfmdecays(iunit, :)   = mean(daydecays, 1, "omitmissing", "Weights", unitdayspikes);
    resfeats.gausssigma(iunit, :)  = mean(daysigma, 1, "omitmissing", "Weights", unitdayspikes);
    resfeats.avgwforms(iunit, :)   = mean(daywvfm, 1, "omitmissing", "Weights", unitdayspikes);

    % this averaging is not that proper, check later...
    resfeats.tfitprofile(iunit, :) = mean(daytfitprofile, 1, "omitmissing", "Weights", unitdayspikes);
    resfeats.cfitprofile(iunit, :) = mean(daycfitprofile, 1, "omitmissing", "Weights", unitdayspikes);
    

end

% resfeatures.wfmamps   = max(abs(avgwfmall), [], 2);

end



function newtimes = fixWaveformTimes(wvtimes, spikes, thresmax)
        
wvtimes(wvtimes < 0 | wvtimes > thresmax) = nan;

medval = median(wvtimes, 2, 'omitmissing');
t2ptimeall = wvtimes(:, 1, :);
ireplace   = isnan(t2ptimeall);
t2ptimeall(ireplace) = medval(ireplace);

newtimes = mean(squeeze(t2ptimeall), 2, "omitmissing", "Weights",spikes);

end