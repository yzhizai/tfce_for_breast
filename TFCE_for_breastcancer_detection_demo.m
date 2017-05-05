function TFCE_for_breastcancer_detection_demo

filenames = cellstr(spm_select(Inf, 'image'));

for aa = 1:numel(filenames)
    filename = filenames{aa};
    [pat, tit, ext, ~] = spm_fileparts(filename);
    V = spm_vol(filename);
    Img = spm_read_vols(V);
    Img_tfce = TFCE_for_breastcancer_detection(Img, 0.5, 2, 6);
    V.fname = fullfile(pat, [tit, '_tfce', ext]);
    V = spm_create_vol(V);
    spm_write_vol(V, Img_tfce);  
end
