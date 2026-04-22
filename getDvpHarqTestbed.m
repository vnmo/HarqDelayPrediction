function [dvp] = getDvpHarqTestbed(d, delta, pHarq, f, T, M, Qmax, mcs)
% getDvpHarqTestbed
% Load queue-length PMF from pure_queue_len_outputs/comparison_queue_pmf.csv
% and use the MCS-specific column as waitpdf.
%
% Example CSV columns:
%   queue_length, resMCS5, resMCS9, resMCS10, resMCS12, resMCS14, resMCS20

    %#ok<*INUSD>
    % f, Qmax are kept only to minimize changes in the caller.

    pmfFile = fullfile("pure_queue_len_outputs", "comparison_queue_pmf.csv");

    if ~isfile(pmfFile)
        error("PMF file not found: " + pmfFile);
    end

    tbl = readtable(pmfFile);

    colName = "resMCS" + string(mcs);

    if ~any(strcmp(tbl.Properties.VariableNames, colName))
        error("Column " + colName + " not found in " + pmfFile);
    end

    % Read PMF column for this MCS
    waitpdf = tbl.(colName);
    waitpdf = waitpdf(:);

    % Remove NaN/Inf just in case
    waitpdf(~isfinite(waitpdf)) = 0;

    % Clamp tiny numerical negatives
    waitpdf(waitpdf < 0) = 0;

    % Normalize
    s = sum(waitpdf);
    if s <= 0
        error("Invalid PMF for " + colName + ": sum(waitpdf) <= 0");
    end
    waitpdf = waitpdf / s;

    % DVP calculation
    dvp = 0;

    for k = 0:length(waitpdf)-1
        k_d = min(M, floor(((d / T) - k + delta) / (delta + 2)));
        k_d_CSO = min(M, floor(((d / T) - k + delta - 1) / (delta + 2)));

        k_d = max(k_d, 0);
        k_d_CSO = max(k_d_CSO, 0);

        if k_d == 0
            harq_term_1 = 1;
        else
            harq_term_1 = prod(pHarq(1:k_d));
        end

        if k_d_CSO == 0
            harq_term_2 = 1;
        else
            harq_term_2 = prod(pHarq(1:k_d_CSO));
        end

        dvp = dvp ...
            + (1 - k_d * 0) * waitpdf(k+1) * harq_term_1 ...
            + (k_d * 0)     * waitpdf(k+1) * harq_term_2;
    end
end