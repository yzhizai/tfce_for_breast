function Img_tfce = TFCE_for_breastcancer_detection(Img, E, H, adj)
if nargin < 2
    E = 0.5;
    H = 2;
    adj = 6;
elseif nargin < 4
    error('you should provide 1 or 4 args');
end

[dim_x, dim_y, dim_z] = size(Img);
Img_temp = Img;
Img_temp(Img == 0) = [];
Img_tfce = zeros(size(Img));


% init_threshold = quantile(Img_temp(:), 0.5);% discarding the small values.
init_threshold = 0.1;


for cc = 1:dim_z
    for bb = 1:dim_y
        for aa = 1:dim_x
            adj_vox = [];
            vox_hp = Img(aa, bb, cc);
            if vox_hp > init_threshold
                try
                    adj_vox = Img(aa - 1:aa + 1, bb - 1:bb + 1, cc - 1:cc + 1); % adjacent voxels
                end
                if ~isempty(adj_vox)
                    tfce_for_vox = get_tfce_val(adj_vox, vox_hp, E, H, adj);
                else
                    tfce_for_vox = 0;
                end
                Img_tfce(aa, bb, cc) = tfce_for_vox;
            end
        end
    end
end


function val = get_tfce_val(adj_vox, vox_hp, E, H, adj)
if adj == 6
    subindx = [5, 11, 13, 14, 15, 17, 23];
elseif adj == 26
    subindx = 1:27;
else
    error('you should provide 6 or 26 for "adj" argument');
end

adj_vox_val = adj_vox(subindx);

% vox_hp_array = 0:dh:vox_hp;
vox_hp_array = linspace(0, vox_hp, 10);
temp = zeros(numel(vox_hp_array), 2);
for aa = 1:numel(vox_hp_array)
    temp(aa, :) = [sum(adj_vox_val >= vox_hp_array(aa)), vox_hp_array(aa)];
end
val = sum(temp(:, 1).^E.*temp(:, 2).^H);